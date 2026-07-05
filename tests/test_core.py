from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from genome_assembly import AssemblyConfig, assemble_short_reads, n50, nx, simulate_reads
from genome_assembly.benchmark import parse_backend_list, run_benchmark
from genome_assembly.io import FastaRecord, FastqRecord, read_fasta, read_fastq, write_fasta, write_fastq
from genome_assembly.kmers import canonical_kmer, iter_kmers, reverse_complement
from genome_assembly.metrics import assembly_stats
from genome_assembly.cython_backend import cython_available, require_cython_backend
from genome_assembly.native import native_available, require_native


class KmerTests(unittest.TestCase):
    def test_iter_kmers_skips_ambiguous(self):
        self.assertEqual(list(iter_kmers("ACNTG", 3)), [])
        self.assertEqual(list(iter_kmers("ACGT", 3)), ["ACG", "CGT"])

    def test_reverse_complement_and_canonical(self):
        self.assertEqual(reverse_complement("ACGTTA"), "TAACGT")
        self.assertEqual(canonical_kmer("TAACGT"), "ACGTTA")


class AssemblyTests(unittest.TestCase):
    def test_linear_reads_generate_single_contig(self):
        result = assemble_short_reads(["ACGTTGA"], AssemblyConfig(k=3))
        self.assertEqual([contig.sequence for contig in result.contigs], ["ACGTTGA"])
        self.assertEqual(result.summary.nodes, 5)
        self.assertEqual(result.summary.edges, 4)

    def test_overlapping_reads_merge_into_contig(self):
        result = assemble_short_reads(["ACGTTA", "GTTACC"], AssemblyConfig(k=3))
        self.assertEqual(result.contigs[0].sequence, "ACGTTACC")

    def test_min_abundance_filters_edges(self):
        with self.assertRaises(ValueError):
            assemble_short_reads(["ACGTTA"], AssemblyConfig(k=3, min_abundance=2))


class MetricTests(unittest.TestCase):
    def test_nx_metrics(self):
        contigs = ["A" * 100, "C" * 50, "G" * 25]
        self.assertEqual(n50(contigs), 100)
        self.assertEqual(nx(contigs, 90), 25)

    def test_assembly_stats(self):
        report = assembly_stats(["ACGT", "GGGG"], reference_length=10)
        self.assertEqual(report["contigs"], 2)
        self.assertEqual(report["total_bp"], 8)
        self.assertEqual(report["n50"], 4)
        self.assertEqual(report["ng50"], 4)
        self.assertEqual(report["ng90"], 0)
        self.assertEqual(report["coverage_percent"], 80.0)


class SimulationAndIoTests(unittest.TestCase):
    def test_simulation_is_deterministic(self):
        first = simulate_reads("ACGTACGTACGT", read_length=4, coverage=2, seed=1)
        second = simulate_reads("ACGTACGTACGT", read_length=4, coverage=2, seed=1)
        self.assertEqual(first.reads, second.reads)
        self.assertEqual(len(first.reads), 6)

    def test_fasta_and_fastq_roundtrip(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fasta = tmp / "records.fasta"
            fastq = tmp / "records.fastq"

            write_fasta([FastaRecord("seq1", "ACGT")], fasta)
            write_fastq([FastqRecord("read1", "ACGT", "IIII")], fastq)

            self.assertEqual(read_fasta(fasta)[0].sequence, "ACGT")
            self.assertEqual(read_fastq(fastq)[0].quality, "IIII")


class GraphCleaningTests(unittest.TestCase):
    def test_tip_clipping_removes_error_branch(self):
        # Strong backbone plus a short erroneous dead-end branch off node "TAC".
        reads = ["GATTACAGG"] * 3 + ["TTACG"]
        dirty = assemble_short_reads(reads, AssemblyConfig(k=3))
        self.assertGreater(len(dirty.contigs), 1)  # branch splits the backbone

        cleaned = assemble_short_reads(reads, AssemblyConfig(k=3, tip_length=5))
        self.assertEqual([c.sequence for c in cleaned.contigs], ["GATTACAGG"])
        self.assertEqual(cleaned.summary.tips_removed, 1)

    def test_tip_clipping_preserves_clean_linear_assembly(self):
        reads = ["GATTACAGG"] * 3
        cleaned = assemble_short_reads(reads, AssemblyConfig(k=3, tip_length=10, bubble_length=10))
        self.assertEqual([c.sequence for c in cleaned.contigs], ["GATTACAGG"])
        self.assertEqual(cleaned.summary.tips_removed, 0)
        self.assertEqual(cleaned.summary.bubble_edges_removed, 0)

    def test_bubble_popping_keeps_highest_coverage_path(self):
        # Two parallel paths CAT->GAT differing by a middle substitution.
        reads = ["CATCGAT"] * 3 + ["CATTGAT"]
        dirty = assemble_short_reads(reads, AssemblyConfig(k=3))
        self.assertEqual(len(dirty.contigs), 2)

        cleaned = assemble_short_reads(reads, AssemblyConfig(k=3, bubble_length=8))
        self.assertEqual([c.sequence for c in cleaned.contigs], ["CATCGAT"])
        self.assertEqual(cleaned.summary.bubble_edges_removed, 4)

    def test_cleaning_is_backend_agnostic(self):
        reads = ["GATTACAGG"] * 3 + ["TTACG"]
        python_result = assemble_short_reads(reads, AssemblyConfig(k=3, tip_length=5, backend="python"))
        for backend in ("cython", "native"):
            available = cython_available() if backend == "cython" else native_available()
            if not available:
                continue
            result = assemble_short_reads(reads, AssemblyConfig(k=3, tip_length=5, backend=backend))
            self.assertEqual(
                [c.sequence for c in result.contigs],
                [c.sequence for c in python_result.contigs],
                msg=f"{backend} cleaning diverged from python",
            )
            self.assertEqual(result.summary.tips_removed, python_result.summary.tips_removed)


class BenchmarkTests(unittest.TestCase):
    def test_parse_backend_list(self):
        self.assertEqual(parse_backend_list("python, cython,native"), ["python", "cython", "native"])
        with self.assertRaises(ValueError):
            parse_backend_list(" , ")

    def test_run_python_benchmark(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta = Path(tmpdir) / "reference.fasta"
            write_fasta([FastaRecord("ref", "ACGTACGTACGT")], fasta)

            report = run_benchmark(
                fasta,
                backends=["python"],
                k=3,
                read_length=4,
                coverage=2,
                seed=1,
            )

        self.assertEqual(report["benchmark"]["reference_name"], "ref")
        self.assertEqual(report["benchmark"]["read_count"], 6)
        self.assertEqual(report["runs"][0]["backend"], "python")
        self.assertEqual(report["runs"][0]["status"], "ok")
        self.assertIn("wall_time_seconds", report["runs"][0])
        self.assertIn("stats", report["runs"][0])


class NativeBridgeTests(unittest.TestCase):
    def test_native_available_returns_bool(self):
        self.assertIsInstance(native_available(), bool)

    def test_missing_native_extension_has_actionable_error(self):
        if native_available():
            self.assertIsNotNone(require_native())
            return

        with self.assertRaisesRegex(RuntimeError, "maturin develop --release"):
            require_native()

    def test_native_backend_matches_python_when_available(self):
        if not native_available():
            self.skipTest("native extension is not installed")

        reads = ["ACGTTA", "GTTACC"]
        python_result = assemble_short_reads(reads, AssemblyConfig(k=3, backend="python"))
        native_result = assemble_short_reads(reads, AssemblyConfig(k=3, backend="native"))

        self.assertEqual(python_result.summary, native_result.summary)
        self.assertEqual(
            [contig.sequence for contig in python_result.contigs],
            [contig.sequence for contig in native_result.contigs],
        )

    def test_native_multithread_output_matches_single_thread(self):
        if not native_available():
            self.skipTest("native extension is not installed")

        reads = simulate_reads("ACGTACGTTAGGCCATATCGATCGTAGCTAG" * 4, read_length=12, coverage=8, seed=3).reads
        self.assertGreater(len(reads), 4)

        single = assemble_short_reads(reads, AssemblyConfig(k=5, backend="native", threads=1))
        multi = assemble_short_reads(reads, AssemblyConfig(k=5, backend="native", threads=4))

        self.assertEqual(single.summary, multi.summary)
        self.assertEqual(
            [contig.sequence for contig in single.contigs],
            [contig.sequence for contig in multi.contigs],
        )


class CythonBridgeTests(unittest.TestCase):
    def test_cython_available_returns_bool(self):
        self.assertIsInstance(cython_available(), bool)

    def test_missing_cython_extension_has_actionable_error(self):
        if cython_available():
            self.assertIsNotNone(require_cython_backend())
            return

        with self.assertRaisesRegex(RuntimeError, "build_ext --inplace"):
            require_cython_backend()

    def test_cython_import_failure_has_actionable_error(self):
        with patch("genome_assembly.cython_backend.import_module", side_effect=ImportError):
            with self.assertRaisesRegex(RuntimeError, "build_ext --inplace"):
                require_cython_backend()

    def test_cython_backend_matches_python_when_available(self):
        if not cython_available():
            self.skipTest("Cython extension is not built")

        reads = ["ACGTTA", "GTTACC"]
        python_result = assemble_short_reads(reads, AssemblyConfig(k=3, backend="python"))
        cython_result = assemble_short_reads(reads, AssemblyConfig(k=3, backend="cython"))

        self.assertEqual(python_result.summary, cython_result.summary)
        self.assertEqual(
            [contig.sequence for contig in python_result.contigs],
            [contig.sequence for contig in cython_result.contigs],
        )


if __name__ == "__main__":
    unittest.main()
