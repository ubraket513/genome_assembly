import matplotlib.pyplot as plt
import numpy as np
import de_bruijn as db
import random as rd
import multiprocessing as mp
import networkx as nx
import warnings, sys
from loguru import logger

warnings.filterwarnings("ignore")

# configure Loguru logger
logger.remove()
logger.add(
    sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} - {level} - {message}", level="INFO"
)

def read_fasta(fasta_file):
    with open(fasta_file, "r") as f:
        lines = f.readlines()
    sequence = np.array(list("".join(line.strip() for line in lines[1:])))  # retrieve the whole genomic sequence from FASTA
    return sequence

def generate_reads(sequence, read_length=200, coverage=3):
    size = len(sequence)
    reads = []
    coverage_array = np.zeros_like(list(sequence), dtype=int)
    # repeat until every dna base in the sequence is covered at least 3 times by generated reads
    itr = 0
    while not np.all(np.greater_equal(coverage_array, coverage*np.ones_like(coverage_array))):
        # if the sequence is not entirely 3x covered by reads
        read_index = rd.randint(0, size-read_length)
        # random read of length 200 starts here
        candidate_coverage = coverage_array[read_index:read_index + read_length]
        candidate_read = sequence[read_index:read_index + read_length]
        # if this random read covers any dna base that has yet 3x covered
        if np.any(np.greater_equal(coverage*np.ones_like(candidate_coverage), candidate_coverage)):
            reads.append(candidate_read)
            coverage_array[read_index:read_index + read_length] += 1 # update the coverage
        if itr % 10000 == 0:
            logger.info(f"{itr} iterations")
        itr += 1
    logger.success(f"simulated {len(reads)} total reads from the sequence")
    return reads, coverage_array

def visualize_coverage(coverage_array):
    plt.plot(coverage_array, 'go', markersize=3)  # Plot as green dots
    plt.xlabel("Position in the Genome")
    plt.ylabel("Coverage")
    plt.title("Coverage Plot of NCBI COVID-19")
    plt.axhline(y=3, color='r', linestyle='-')  # Add red horizontal line at y=3
    plt.show()

def generate_Nx_stat(contig_array, total_genome_length: int, x):
    """
    Construct a Nx-statistics of the contig array with the value x specified
    """
    contig_len_array = sorted([len(contig) for contig in contig_array], reverse=True)
    nX = 0
    index = 0
    while (nX < total_genome_length * x / 100) and (index < len(contig_len_array)):
        nX += contig_len_array[index]
        index += 1
    return nX

def visualize_Nx_stat(contig_array, total_genome_length):
    """
    Visualize the Nx statistics for x=[1,...,100]
    """
    x_lst = np.linspace(1, 100, endpoint=True)
    Nx_stat = [generate_Nx_stat(contig_array, total_genome_length, x) for x in x_lst]
    plt.plot(x_lst, Nx_stat, 'ro-', markersize=3)
    plt.xlabel("Value of x")
    plt.ylabel("Nx Statistics")
    plt.title("Nx Plot") 
    plt.show()

def generate_contigs(graph):
    """
    Generate contigs from de bruijn graph
    """
    pass

if __name__ == '__main__':
    fasta_file = "genomic.fna"
    sequence = read_fasta(fasta_file)
    total_genome_length = len(sequence)
    # simulate reads of coverage 3x
    reads, coverage_array = generate_reads(sequence=sequence)
    visualize_coverage(coverage_array)
    # initialize the de bruijn graph
    de_bruijn = db.create(reads=reads, k=20)
    # compressed_de_bruijn = db.compress(de_bruijn, k=20)
    logger.info("Generating contigs from the de Bruijn graph")
    contig_array = generate_contigs(de_bruijn)
    logger.info("Computing Nx statistics")
    visualize_Nx_stat(contig_array, total_genome_length)