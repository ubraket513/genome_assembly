"""Small compatibility demo for the package implementation.

Prefer the production CLI:

    ga simulate genomic.fna reads.fastq
    ga assemble reads.fastq --k 31 --outdir assembly_out
"""

from genome_assembly import AssemblyConfig, assemble_short_reads, n50, nx, read_fasta, simulate_reads


def read_fasta_sequence(fasta_file):
    return read_fasta(fasta_file)[0].sequence


def generate_reads(sequence, read_length=200, coverage=2, seed=7):
    simulated = simulate_reads(sequence, read_length=read_length, coverage=coverage, seed=seed)
    return simulated.reads, simulated.coverage


def generate_Nx_stat(contig_array, total_genome_length: int, x):
    return nx(contig_array, x, reference_length=total_genome_length)


def generate_contigs(graph):
    return [contig.sequence for contig in graph.compact_contigs()]


if __name__ == "__main__":
    sequence = read_fasta_sequence("genomic.fna")
    reads, coverage_array = generate_reads(sequence=sequence, coverage=3, seed=7)
    result = assemble_short_reads(reads, AssemblyConfig(k=20))
    contigs = [contig.sequence for contig in result.contigs]
    print(f"Simulated {len(reads)} reads")
    print(f"Mean coverage: {sum(coverage_array) / len(coverage_array):.2f}x")
    print(f"Assembled {len(contigs)} contigs")
    print(f"N50: {n50(contigs, reference_length=len(sequence))}")
