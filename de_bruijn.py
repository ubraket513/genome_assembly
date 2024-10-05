import networkx as nx
import numpy as np
from tqdm import trange
import warnings, sys
from loguru import logger

warnings.filterwarnings("ignore")

# configure Loguru logger
logger.remove()
logger.add(
    sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} - {level} - {message}", level="INFO"
)

def get_kmers(read, k:int):
  """
  Generates all possible k-mers from a given read.

  Args:
    read: The input read (string).
    k: The length of k-mers.

  Returns:
    A NumPy array of k-mers.
  """
  indices = np.arange(len(read) - k + 1)[:, np.newaxis] + np.arange(k) # Generate a sliding window of indices
  kmers_np = np.take(read, indices) # Extract k-mers using advanced indexing
  kmers = [''.join(kmer) for kmer in kmers_np] # Convert k-mer array back to a list of strings
  return kmers

def create(reads, k: int):
  """
  Creates a De Bruijn graph from a list of k-mers, removing duplicate edges.

  Args:
    kmers: A list of k-mers.

  Returns:
    A dictionary representing the De Bruijn graph where keys are nodes
    (k-1-mers) and values are lists of their unique neighbors.
  """
  # Create k-mers
  kmers = []
  for index in trange(0, len(reads)): # iterate through all reads
    kmers.extend(get_kmers(reads[index], k+1)) # generate k+1 mers from the read for edges

  # Create a directed graph
  graph = nx.DiGraph()
  for kmer in kmers:
    prefix, suffix = kmer[:-1], kmer[1:]
    graph.add_edge(prefix, suffix)  # Add edges to the graph
  logger.success("Generated de bruijn graph")
  return graph

def should_remove_node(graph, node):
    return graph.in_degree(node) == 1 and graph.out_degree(node) == 1


def compress(G, k:int):  
  for node in list(G.nodes()):
      if should_remove_node(G, node):
          in_node = list(G.in_edges(node))[0][0]
          out_node = list(G.out_edges(node))[0][1]
          G.add_edge(in_node, out_node)
          G.remove_node(node)