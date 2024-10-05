import matplotlib.pyplot as plt
import numpy as np
import de_bruijn as db
import random as rd
import multiprocessing as mp
import warnings, sys
from loguru import logger

warnings.filterwarnings("ignore")

# configure Loguru logger
logger.remove()
logger.add(
    sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} - {level} - {message}", level="INFO"
)

def generate_reads(fasta_file, read_length=200, coverage=3):
    with open(fasta_file, "r") as f:
        lines = f.readlines()
    sequence = "".join(line.strip() for line in lines[1:])  # retrieve the whole genomic sequence from FASTA
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
    plt.plot(coverage_array)
    plt.xlabel("Position in the Genome")
    plt.ylabel("Coverage")
    plt.title("Coverage Plot of NCBI COVID-19")
    plt.show()


if __name__ == '__main__':
    fasta_file = "genomic.fna"
    reads, coverage_array = generate_reads(fasta_file)
    genome_length = len(reads[0]) * len(reads)  # Adjust if needed
    dbgraphs = []
    visualize_coverage(coverage_array)
    
    def func(read):
        return db.create(read=read, k=20)
    
    with mp.Pool(processes=150) as pool:
        results = pool.map_async(func, reads)  # Use map_async
        # Wait for the processes to finish
        results.wait()
        dbgraphs = results.get()  # Get the actual results
    de_bruijn = db.merge(dbgraphs)
    logger.info(f"Merging de bruijn graphs. {de_bruijn}")
    # print(db.compress_debruijn_graph(de_bruijn))
    # visualize_coverage(reads, genome_length)