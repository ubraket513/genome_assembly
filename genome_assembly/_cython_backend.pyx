# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True

cdef str _normalize_sequence(object sequence):
    return "".join(str(sequence).split()).upper()


cdef bint _is_unambiguous_dna(str sequence):
    cdef str base
    for base in sequence:
        if base != "A" and base != "C" and base != "G" and base != "T":
            return False
    return True


cpdef list count_kmers(object reads, Py_ssize_t k, bint skip_ambiguous=True):
    cdef dict counts = {}
    cdef object read
    cdef str sequence
    cdef str kmer
    cdef Py_ssize_t start
    cdef Py_ssize_t length

    if k < 1:
        raise ValueError("k must be >= 1")

    for read in reads:
        sequence = _normalize_sequence(read)
        length = len(sequence)
        if length < k:
            continue

        for start in range(length - k + 1):
            kmer = sequence[start:start + k]
            if skip_ambiguous and not _is_unambiguous_dna(kmer):
                continue
            counts[kmer] = counts.get(kmer, 0) + 1

    return sorted(counts.items())


cpdef tuple build_edges(
    object reads,
    Py_ssize_t node_k,
    Py_ssize_t min_abundance=1,
    bint skip_ambiguous=True,
):
    cdef Py_ssize_t edge_k
    cdef list counts
    cdef Py_ssize_t raw_edge_count = 0
    cdef list edges = []
    cdef str kmer
    cdef Py_ssize_t count

    if node_k < 1:
        raise ValueError("node_k must be >= 1")
    if min_abundance < 1:
        raise ValueError("min_abundance must be >= 1")

    edge_k = node_k + 1
    counts = count_kmers(reads, edge_k, skip_ambiguous)

    for kmer, count in counts:
        raw_edge_count += count
        if count < min_abundance:
            continue
        edges.append((kmer[:node_k], kmer[1:], kmer, count))

    return raw_edge_count, edges


def compact_contigs(Py_ssize_t node_k, object edge_rows, Py_ssize_t min_length=0):
    cdef dict out_edges = {}
    cdef dict in_edges = {}
    cdef set nodes = set()
    cdef set visited = set()
    cdef list contigs = []
    cdef object edge
    cdef object node
    cdef object row
    cdef Py_ssize_t discovery_index = 0

    if node_k < 1:
        raise ValueError("node_k must be >= 1")
    if min_length < 0:
        raise ValueError("min_length must be >= 0")

    for row in edge_rows:
        out_edges.setdefault(row[0], []).append(row)
        in_edges.setdefault(row[1], []).append(row)
        nodes.add(row[0])
        nodes.add(row[1])

    def edge_key(edge):
        return edge[0], edge[1]

    def is_one_in_one_out(node):
        return len(in_edges.get(node, [])) == 1 and len(out_edges.get(node, [])) == 1

    def walk(start_edge):
        visited.add(edge_key(start_edge))
        sequence = start_edge[2]
        counts = [start_edge[3]]
        current = start_edge[1]
        edge_count = 1

        while is_one_in_one_out(current):
            next_edges = out_edges.get(current, [])
            if not next_edges:
                break
            next_edge = next_edges[0]
            if edge_key(next_edge) in visited:
                break
            visited.add(edge_key(next_edge))
            sequence += next_edge[2][len(next_edge[2]) - 1]
            counts.append(next_edge[3])
            current = next_edge[1]
            edge_count += 1

        return sequence, sum(counts) / len(counts), edge_count

    for node in sorted(nodes):
        if is_one_in_one_out(node):
            continue
        for edge in out_edges.get(node, []):
            if edge_key(edge) in visited:
                continue
            sequence, mean_abundance, edge_count = walk(edge)
            if len(sequence) >= min_length:
                discovery_index += 1
                contigs.append((f"contig_{discovery_index}", sequence, mean_abundance, edge_count))

    for edge in edge_rows:
        if edge_key(edge) in visited:
            continue
        sequence, mean_abundance, edge_count = walk(edge)
        if len(sequence) >= min_length:
            discovery_index += 1
            contigs.append((f"contig_{discovery_index}", sequence, mean_abundance, edge_count))

    contigs.sort(key=lambda item: (-len(item[1]), item[0]))
    return [(sequence, mean_abundance, edge_count) for _, sequence, mean_abundance, edge_count in contigs]
