from pathlib import Path
import tempfile
import unittest

from genome_assembly import AssemblyConfig, assemble_short_reads, n50, nx, simulate_reads
from genome_assembly.io import FastaRecord, FastqRecord, read_fasta, read_fastq, write_fasta, write_fastq
from genome_assembly.kmers import canonical_kmer, iter_kmers, reverse_complement
from genome_assembly.metrics import assembly_stats


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


if __name__ == "__main__":
    unittest.main()
