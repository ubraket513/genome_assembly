from collections import defaultdict
import numpy as np

def get_kmers(read:str, k:int):
  """
  Generates all possible k-mers from a given read.

  Args:
    read: The input read (string).
    k: The length of k-mers.

  Returns:
    A NumPy array of k-mers.
  """
  read_np = np.array(list(read)) # Convert the read to a NumPy array for efficient slicing

  indices = np.arange(len(read) - k + 1)[:, np.newaxis] + np.arange(k) # Generate a sliding window of indices

  kmers_np = np.take(read_np, indices) # Extract k-mers using advanced indexing

  kmers = [''.join(kmer) for kmer in kmers_np] # Convert k-mer array back to a list of strings

  return kmers


def create(read:str, k: int):
  """
  Creates a De Bruijn graph from a list of k-mers, removing duplicate edges.

  Args:
    kmers: A list of k-mers.

  Returns:
    A dictionary representing the De Bruijn graph where keys are nodes
    (k-1-mers) and values are lists of their unique neighbors.
  """
  kmers = get_kmers(read, k+1) # generate k+1 mers from the read

  graph = defaultdict(set)  # base class for a de bruijn graph

  for kmer in kmers:
    # generate two candidate de bruijn nodes of a particular k+1 mer
    prefix, suffix = kmer[:-1], kmer[1:]
    graph[prefix].add(suffix)  # add an edge between two candidate nodes

  graph = {k: list(v) for k, v in graph.items()}
  print("created a de bruijn graph.")
  return graph


def merge(graphs: list):
  """
  Merges multiple De Bruijn graphs into a single graph.

  Args:
    graphs: A list of De Bruijn graphs (dictionaries).

  Returns:
    A dictionary representing the merged De Bruijn graph.
  """
  merged_graph = defaultdict(set)
  for graph in graphs:
    for node, neighbors in graph.items():
      merged_graph[node].update(neighbors)  # Use update to add multiple neighbors
  # Convert sets back to lists for consistency
  merged_graph = {k: list(v) for k, v in merged_graph.items()}
  return merged_graph


def compress(graph):
    """
    Compresses a De Bruijn graph by merging non-branching paths.

    Args:
      graph: A dictionary representing the De Bruijn graph.

    Returns:
      A dictionary representing the compressed De Bruijn graph.
    """
    compressed_graph = {}
    for node in graph:
        if len(graph[node]) == 1 and len(graph.get(graph[node][0], [])) == 1:  # Non-branching path
            start_node = node
            path = [node]
            while len(graph[path[-1]]) == 1 and len(graph.get(graph[path[-1]][0], [])) == 1:
                path.append(graph[path[-1]][0])
            compressed_graph[start_node] = [path[-1]]
        elif node not in compressed_graph:  # Branching node
            compressed_graph[node] = graph[node]
    print("compressed a de bruijn graph.")
    return compressed_graph