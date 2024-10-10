import networkx as nx
import numpy as np
import warnings, sys
import matplotlib.pyplot as plt
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
  Creates a De Bruijn graph from a list of reads with (k+1)-mer edge labels.

  Args:
      reads: A list of reads (strings).
      k: The length of k-mers to use for constructing the graph 
          (edges will be k+1-mers).

  Returns:
      A NetworkX DiGraph representing the De Bruijn graph.
  """
  edges = get_kmers(reads, k+1)
  graph = nx.DiGraph()
  
  for kmer in edges:
      prefix, suffix = kmer[:-1], kmer[1:]
      # Add edge with the (k+1)-mer as the label
      graph.add_edge(prefix, suffix, label=kmer)  

  edge_labels = {(u, v): kmer for kmer in get_kmers(reads, k + 1) for u, v in [(kmer[:-1], kmer[1:])]}
  nx.set_edge_attributes(graph, edge_labels, "contig")
  logger.success("Generated De Bruijn graph")
  return graph

def is_node_1_to_1(graph, node):
    """
    Check if the node in the graph has exactly one outgoing edge and one incoming edge
    """
    return graph.in_degree(node) == 1 and graph.out_degree(node) == 1

def generate_contigs(graph):
  """
  Contig is a maximal non-branching path in a de Bruijn graph. This method returns a list of contigs and a list of node tuples that corresponds to each contig.
  """
  contigs = [] # list of tuples representing compressed edge list
  edge_labels = {} # label of each compressed edge
  for v in graph.nodes():
    if not is_node_1_to_1(v): # compression begins if v is not a 1-in-1-out node
      if graph.in_degree(v) > 0: # starting node must not be compressed out
        for edge in graph.out_edges(v): # each outgoing edge (v,w) from v
          w = edge[1] # extract w from edge (v, w)
          non_branching_path = [v, w] # nodes of the path starting from v
          compressed_edge_label = graph.get_edge_data(*edge)["contig"] # edge label of the compressed edge
          while is_node_1_to_1(w): # check if there is a longer non-branching path
            u = graph.out_edges(w)[1]
            non_branching_path.append(u) # extend NonBranchingPath by the outgoing edge (w, u) from w 
            compressed_edge_label += u[-1] # append last string of the path to the existing string
            w = u # w ← u
          compressed_path = tuple(non_branching_path[0], non_branching_path[-1])
          edge_labels[compressed_path] = compressed_edge_label
          contigs.append()
  return contigs, edge_labels

def compress(graph):
  """
  Compresses an existing De Bruijn graph by collapsing maximal non-branching paths.
  """
  compressed_edge_lst, contigs_label = generate_contigs(graph)
  compressed_graph = nx.MultiDiGraph(graph.nodes()) # generate compressed graph from existing k-mer nodes
  compressed_graph.add_edges_from(compressed_edge_lst) # add edges from the compressed list
  nx.set_edge_attributes(graph, contigs_label, "contig") # add contigs as edges
  compressed_graph.remove_nodes_from(list(nx.isolates(compressed_graph))) # remove nodes that are compressed out

"""
def plot(graph):
   # Draw the graph
   pos = nx.spring_layout(graph)  # Choose a layout for the graph
   nx.draw(graph, pos, with_labels=True, node_size=1500, node_color="skyblue")
   # Draw edge labels
   nx.draw_networkx_edge_labels(graph, pos, edge_labels=nx.get_edge_attributes(graph, "label"))
   plt.show()
"""