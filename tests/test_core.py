from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from typer.testing import CliRunner

from genome_assembly import (
    AssemblyConfig,
    assemble_short_reads,
    generate_random_genome,
    minimizer_bucket,
    minimizers,
    n50,
    nx,
    simulate_reads,
    syncmers,
)
from genome_assembly.cli import app
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


class SketchTests(unittest.TestCase):
    def setUp(self):
        self.seq = generate_random_genome(400, seed=11)

    def test_minimizers_are_deterministic_and_ordered(self):
        first = minimizers(self.seq, w=10, m=8)
        second = minimizers(self.seq, w=10, m=8)
        self.assertEqual(first, second)
        self.assertTrue(first)
        positions = [pos for pos, _ in first]
        self.assertEqual(positions, sorted(set(positions)))  # strictly increasing
        for pos, mmer in first:
            self.assertEqual(self.seq[pos : pos + 8], mmer)

    def test_short_sequence_returns_global_minimizer(self):
        result = minimizers("ACGTA", w=50, m=3)
        self.assertEqual(len(result), 1)

    def test_syncmers_have_bounded_spacing(self):
        syncs = syncmers(self.seq, k=15, s=6, t=0)
        self.assertTrue(syncs)
        for pos, kmer in syncs:
            self.assertEqual(len(kmer), 15)
            self.assertEqual(self.seq[pos : pos + 15], kmer)
        # Open syncmers guarantee at least one seed every (k - s + 1) bases.
        positions = [pos for pos, _ in syncs]
        for earlier, later in zip(positions, positions[1:]):
            self.assertLessEqual(later - earlier, 15 - 6 + 1)

    def test_syncmer_offset_out_of_range_raises(self):
        with self.assertRaises(ValueError):
            syncmers(self.seq, k=10, s=4, t=99)

    def test_minimizer_bucket_stable_and_in_range(self):
        buckets = {minimizer_bucket(self.seq[i : i + 21], m=7, num_buckets=16) for i in range(50)}
        self.assertTrue(all(0 <= b < 16 for b in buckets))
        self.assertEqual(
            minimizer_bucket("ACGTACGTACGT", m=5, num_buckets=8),
            minimizer_bucket("ACGTACGTACGT", m=5, num_buckets=8),
        )

    def test_native_sketch_matches_python(self):
        if not native_available():
            self.skipTest("native extension is not installed")
        from genome_assembly import native

        self.assertEqual(minimizers(self.seq, 10, 8), native.minimizers(self.seq, 10, 8))
        self.assertEqual(syncmers(self.seq, 15, 6, 0), native.syncmers(self.seq, 15, 6, 0))
        for i in range(30):
            kmer = self.seq[i : i + 21]
            self.assertEqual(
                minimizer_bucket(kmer, 7, 16), native.minimizer_bucket(kmer, 7, 16)
            )


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


class CLITests(unittest.TestCase):
    def setUp(self):
        self.runner = CliRunner()

    def test_quickstart_runs(self):
        result = self.runner.invoke(app, ["quickstart"])
        self.assertEqual(result.exit_code, 0)
        self.assertIn("three steps", result.output.lower())

    def test_help_shows_examples(self):
        result = self.runner.invoke(app, ["--help"])
        self.assertEqual(result.exit_code, 0)
        self.assertIn("ga quickstart", result.output)

    def test_assemble_bad_params_exit_cleanly(self):
        with tempfile.TemporaryDirectory() as tmp:
            reads = Path(tmp) / "reads.fastq"
            write_fastq([FastqRecord("r1", "ACGT", "IIII")], reads)
            result = self.runner.invoke(
                app, ["assemble", str(reads), "-k", "31", "-o", str(Path(tmp) / "out")]
            )
        # k > read length yields no edges: should exit via a clean typer.Exit
        # (SystemExit), not a leaked ValueError traceback.
        self.assertEqual(result.exit_code, 1)
        self.assertIsInstance(result.exception, SystemExit)


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
        run = report["runs"][0]
        self.assertEqual(run["backend"], "python")
        self.assertEqual(run["status"], "ok")
        self.assertIn("wall_time_seconds", run)
        self.assertIn("stats", run)
        self.assertIn("peak_rss_bytes", run)  # native memory accounting present

    def test_synthetic_genome_benchmark(self):
        report = run_benchmark(
            synthetic_genome_size=5000,
            backends=["python"],
            k=15,
            read_length=100,
            coverage=5,
            seed=3,
        )
        self.assertEqual(report["benchmark"]["reference_name"], "synthetic_5000bp")
        self.assertEqual(report["benchmark"]["reference_length"], 5000)
        self.assertEqual(report["runs"][0]["status"], "ok")

    def test_benchmark_requires_exactly_one_source(self):
        with self.assertRaises(ValueError):
            run_benchmark(backends=["python"])  # neither reference nor size


class GenomeGenerationTests(unittest.TestCase):
    def test_random_genome_is_deterministic(self):
        first = generate_random_genome(500, seed=42)
        second = generate_random_genome(500, seed=42)
        self.assertEqual(first, second)
        self.assertEqual(len(first), 500)
        self.assertTrue(set(first).issubset(set("ACGT")))

    def test_random_genome_seed_changes_output(self):
        self.assertNotEqual(generate_random_genome(500, seed=1), generate_random_genome(500, seed=2))


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
