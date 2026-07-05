use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

mod mdbg;
#[path = "minimizers.rs"]
mod sketch;
mod union_find;

use union_find::UnionFind;

fn normalize_sequence(sequence: &str) -> Vec<u8> {
    sequence
        .bytes()
        .filter(|byte| !byte.is_ascii_whitespace())
        .map(|byte| byte.to_ascii_uppercase())
        .collect()
}

fn is_unambiguous_dna(sequence: &[u8]) -> bool {
    sequence
        .iter()
        .all(|byte| matches!(byte, b'A' | b'C' | b'G' | b'T'))
}

fn count_reads_chunk(reads: &[String], k: usize, skip_ambiguous: bool) -> HashMap<Vec<u8>, usize> {
    let mut counts: HashMap<Vec<u8>, usize> = HashMap::new();
    for read in reads {
        let sequence = normalize_sequence(read);
        if sequence.len() < k {
            continue;
        }
        for start in 0..=(sequence.len() - k) {
            let window = &sequence[start..start + k];
            if skip_ambiguous && !is_unambiguous_dna(window) {
                continue;
            }
            *counts.entry(window.to_vec()).or_insert(0) += 1;
        }
    }
    counts
}

fn merge_counts(
    mut into: HashMap<Vec<u8>, usize>,
    other: HashMap<Vec<u8>, usize>,
) -> HashMap<Vec<u8>, usize> {
    for (kmer, count) in other {
        *into.entry(kmer).or_insert(0) += count;
    }
    into
}

pub fn count_kmers_impl(
    reads: &[String],
    k: usize,
    skip_ambiguous: bool,
    threads: usize,
) -> Result<Vec<(String, usize)>, String> {
    if k == 0 {
        return Err("k must be >= 1".to_string());
    }
    if threads == 0 {
        return Err("threads must be >= 1".to_string());
    }

    // Integer merge is associative and commutative, so any thread count yields
    // the same counts; the final sort makes the row order deterministic too.
    let counts = if threads == 1 || reads.len() <= 1 {
        count_reads_chunk(reads, k, skip_ambiguous)
    } else {
        let chunk_size = reads.len().div_ceil(threads);
        let pool = ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .map_err(|err| err.to_string())?;
        pool.install(|| {
            reads
                .par_chunks(chunk_size)
                .map(|chunk| count_reads_chunk(chunk, k, skip_ambiguous))
                .reduce(HashMap::new, merge_counts)
        })
    };

    let mut rows: Vec<(String, usize)> = Vec::with_capacity(counts.len());
    for (kmer, count) in counts {
        let kmer = String::from_utf8(kmer)
            .map_err(|_| "sequence contained invalid UTF-8 after normalization".to_string())?;
        rows.push((kmer, count));
    }
    rows.sort_by(|left, right| left.0.cmp(&right.0));
    Ok(rows)
}

pub type EdgeRow = (String, String, String, usize);
pub type ContigRow = (String, f64, usize);

pub fn build_edges_impl(
    reads: &[String],
    node_k: usize,
    min_abundance: usize,
    skip_ambiguous: bool,
    threads: usize,
) -> Result<(usize, Vec<EdgeRow>), String> {
    if node_k == 0 {
        return Err("node_k must be >= 1".to_string());
    }
    if min_abundance == 0 {
        return Err("min_abundance must be >= 1".to_string());
    }

    let edge_k = node_k + 1;
    let counts = count_kmers_impl(reads, edge_k, skip_ambiguous, threads)?;
    let raw_edge_count = counts.iter().map(|(_, count)| *count).sum();
    let mut edges = Vec::new();

    for (kmer, count) in counts {
        if count < min_abundance {
            continue;
        }
        let prefix = kmer[..node_k].to_string();
        let suffix = kmer[1..].to_string();
        edges.push((prefix, suffix, kmer, count));
    }

    Ok((raw_edge_count, edges))
}

fn edge_key(edge: &EdgeRow) -> (String, String) {
    (edge.0.clone(), edge.1.clone())
}

fn is_one_in_one_out(
    node: &str,
    in_edges: &BTreeMap<String, Vec<EdgeRow>>,
    out_edges: &BTreeMap<String, Vec<EdgeRow>>,
) -> bool {
    in_edges.get(node).map_or(0, Vec::len) == 1 && out_edges.get(node).map_or(0, Vec::len) == 1
}

fn walk_contig(
    start_edge: EdgeRow,
    in_edges: &BTreeMap<String, Vec<EdgeRow>>,
    out_edges: &BTreeMap<String, Vec<EdgeRow>>,
    visited: &mut BTreeSet<(String, String)>,
) -> ContigRow {
    visited.insert(edge_key(&start_edge));
    let mut sequence = start_edge.2.clone();
    let mut count_sum = start_edge.3;
    let mut current = start_edge.1.clone();
    let mut edge_count = 1;

    while is_one_in_one_out(&current, in_edges, out_edges) {
        let Some(next_edges) = out_edges.get(&current) else {
            break;
        };
        let Some(next_edge) = next_edges.first() else {
            break;
        };
        let key = edge_key(next_edge);
        if visited.contains(&key) {
            break;
        }

        visited.insert(key);
        if let Some(last_base) = next_edge.2.chars().last() {
            sequence.push(last_base);
        }
        count_sum += next_edge.3;
        current = next_edge.1.clone();
        edge_count += 1;
    }

    (sequence, count_sum as f64 / edge_count as f64, edge_count)
}

pub fn compact_contigs_impl(
    node_k: usize,
    edge_rows: &[EdgeRow],
    min_length: usize,
) -> Result<Vec<ContigRow>, String> {
    if node_k == 0 {
        return Err("node_k must be >= 1".to_string());
    }

    let mut out_edges: BTreeMap<String, Vec<EdgeRow>> = BTreeMap::new();
    let mut in_edges: BTreeMap<String, Vec<EdgeRow>> = BTreeMap::new();
    let mut nodes: BTreeSet<String> = BTreeSet::new();
    let mut visited: BTreeSet<(String, String)> = BTreeSet::new();
    let mut contigs: Vec<(String, String, f64, usize)> = Vec::new();
    let mut discovery_index = 0usize;

    for edge in edge_rows {
        out_edges
            .entry(edge.0.clone())
            .or_default()
            .push(edge.clone());
        in_edges
            .entry(edge.1.clone())
            .or_default()
            .push(edge.clone());
        nodes.insert(edge.0.clone());
        nodes.insert(edge.1.clone());
    }

    for node in &nodes {
        if is_one_in_one_out(node, &in_edges, &out_edges) {
            continue;
        }
        let Some(edges) = out_edges.get(node) else {
            continue;
        };
        for edge in edges.clone() {
            if visited.contains(&edge_key(&edge)) {
                continue;
            }
            let (sequence, mean_abundance, edge_count) =
                walk_contig(edge, &in_edges, &out_edges, &mut visited);
            if sequence.len() >= min_length {
                discovery_index += 1;
                contigs.push((
                    format!("contig_{discovery_index}"),
                    sequence,
                    mean_abundance,
                    edge_count,
                ));
            }
        }
    }

    for edge in edge_rows {
        if visited.contains(&edge_key(edge)) {
            continue;
        }
        let (sequence, mean_abundance, edge_count) =
            walk_contig(edge.clone(), &in_edges, &out_edges, &mut visited);
        if sequence.len() >= min_length {
            discovery_index += 1;
            contigs.push((
                format!("contig_{discovery_index}"),
                sequence,
                mean_abundance,
                edge_count,
            ));
        }
    }

    contigs.sort_by(|left, right| {
        right
            .1
            .len()
            .cmp(&left.1.len())
            .then_with(|| left.1.cmp(&right.1))
    });

    Ok(contigs
        .into_iter()
        .map(|(_, sequence, mean_abundance, edge_count)| (sequence, mean_abundance, edge_count))
        .collect())
}

/// Reconstruct a single unitig (one union-find component) into a contig row.
///
/// The component's edges form a simple path or cycle. Paths start at the edge
/// whose prefix node has no in-edge inside the component; cycles start at the
/// lexicographically smallest edge so the rotation is deterministic and matches
/// the sequential walker.
fn reconstruct_component(component: &[usize], edge_rows: &[EdgeRow]) -> ContigRow {
    let mut next_by_prefix: HashMap<&str, usize> = HashMap::new();
    let mut has_incoming: HashSet<&str> = HashSet::new();
    for &index in component {
        next_by_prefix.insert(edge_rows[index].0.as_str(), index);
        has_incoming.insert(edge_rows[index].1.as_str());
    }

    let start = component
        .iter()
        .filter(|&&index| !has_incoming.contains(edge_rows[index].0.as_str()))
        .min_by(|&&a, &&b| edge_rows[a].2.cmp(&edge_rows[b].2))
        .copied()
        .unwrap_or_else(|| {
            *component
                .iter()
                .min_by(|&&a, &&b| edge_rows[a].2.cmp(&edge_rows[b].2))
                .expect("component is never empty")
        });

    let mut sequence = edge_rows[start].2.clone();
    let mut count_sum = edge_rows[start].3;
    let mut edge_count = 1usize;
    let mut visited: HashSet<usize> = HashSet::new();
    visited.insert(start);
    let mut current = edge_rows[start].1.as_str();

    while let Some(&next) = next_by_prefix.get(current) {
        if visited.contains(&next) {
            break;
        }
        visited.insert(next);
        if let Some(last) = edge_rows[next].2.chars().last() {
            sequence.push(last);
        }
        count_sum += edge_rows[next].3;
        current = edge_rows[next].1.as_str();
        edge_count += 1;
    }

    (sequence, count_sum as f64 / edge_count as f64, edge_count)
}

/// BCALM2-style parallel compaction: glue edges into unitigs with a union-find,
/// then reconstruct each component in parallel. Output is identical to the
/// sequential walker (same unitigs, same deterministic ordering).
pub fn compact_contigs_uf_impl(
    node_k: usize,
    edge_rows: &[EdgeRow],
    min_length: usize,
    threads: usize,
) -> Result<Vec<ContigRow>, String> {
    if node_k == 0 {
        return Err("node_k must be >= 1".to_string());
    }
    if threads == 0 {
        return Err("threads must be >= 1".to_string());
    }
    if edge_rows.is_empty() {
        return Ok(Vec::new());
    }

    let mut out_edges: BTreeMap<&str, Vec<usize>> = BTreeMap::new();
    let mut in_edges: BTreeMap<&str, Vec<usize>> = BTreeMap::new();
    for (index, edge) in edge_rows.iter().enumerate() {
        out_edges.entry(edge.0.as_str()).or_default().push(index);
        in_edges.entry(edge.1.as_str()).or_default().push(index);
    }

    let is_simple = |node: &str| -> bool {
        in_edges.get(node).map_or(0, Vec::len) == 1 && out_edges.get(node).map_or(0, Vec::len) == 1
    };

    // Glue the single in-edge and single out-edge at every one-in-one-out node.
    let mut uf = UnionFind::new(edge_rows.len());
    let simple_nodes: Vec<&str> = out_edges
        .keys()
        .copied()
        .filter(|node| is_simple(node))
        .collect();
    for node in simple_nodes {
        uf.union(in_edges[node][0], out_edges[node][0]);
    }

    let mut components: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for index in 0..edge_rows.len() {
        components.entry(uf.find(index)).or_default().push(index);
    }
    let component_lists: Vec<Vec<usize>> = components.into_values().collect();

    let reconstruct_all = || -> Vec<ContigRow> {
        component_lists
            .par_iter()
            .map(|component| reconstruct_component(component, edge_rows))
            .collect()
    };
    let rows: Vec<ContigRow> = if threads == 1 || component_lists.len() <= 1 {
        component_lists
            .iter()
            .map(|component| reconstruct_component(component, edge_rows))
            .collect()
    } else {
        let pool = ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .map_err(|err| err.to_string())?;
        pool.install(reconstruct_all)
    };

    let mut contigs: Vec<ContigRow> = rows
        .into_iter()
        .filter(|(sequence, _, _)| sequence.len() >= min_length)
        .collect();
    contigs.sort_by(|left, right| {
        right
            .0
            .len()
            .cmp(&left.0.len())
            .then_with(|| left.0.cmp(&right.0))
    });
    Ok(contigs)
}

#[pyfunction]
#[pyo3(signature = (reads, k, skip_ambiguous = true, threads = 1))]
fn count_kmers(
    reads: Vec<String>,
    k: usize,
    skip_ambiguous: bool,
    threads: usize,
) -> PyResult<Vec<(String, usize)>> {
    count_kmers_impl(&reads, k, skip_ambiguous, threads).map_err(PyValueError::new_err)
}

#[pyfunction]
#[pyo3(signature = (reads, node_k, min_abundance = 1, skip_ambiguous = true, threads = 1))]
fn build_edges(
    reads: Vec<String>,
    node_k: usize,
    min_abundance: usize,
    skip_ambiguous: bool,
    threads: usize,
) -> PyResult<(usize, Vec<EdgeRow>)> {
    build_edges_impl(&reads, node_k, min_abundance, skip_ambiguous, threads)
        .map_err(PyValueError::new_err)
}

#[pyfunction]
#[pyo3(signature = (node_k, edge_rows, min_length = 0, threads = 1))]
fn compact_contigs(
    node_k: usize,
    edge_rows: Vec<EdgeRow>,
    min_length: usize,
    threads: usize,
) -> PyResult<Vec<ContigRow>> {
    // Union-find compaction parallelizes across unitigs; the sequential walker
    // stays the single-thread path and the correctness reference.
    let result = if threads > 1 {
        compact_contigs_uf_impl(node_k, &edge_rows, min_length, threads)
    } else {
        compact_contigs_impl(node_k, &edge_rows, min_length)
    };
    result.map_err(PyValueError::new_err)
}

#[pyfunction]
#[pyo3(signature = (sequence, w, m))]
fn minimizers(sequence: &str, w: usize, m: usize) -> PyResult<Vec<(usize, String)>> {
    let seq = normalize_sequence(sequence);
    let out = sketch::minimizers(&seq, w, m).map_err(PyValueError::new_err)?;
    Ok(out
        .into_iter()
        .map(|(pos, bytes)| (pos, String::from_utf8_lossy(&bytes).into_owned()))
        .collect())
}

#[pyfunction]
#[pyo3(signature = (sequence, k, s, t = 0))]
fn syncmers(sequence: &str, k: usize, s: usize, t: usize) -> PyResult<Vec<(usize, String)>> {
    let seq = normalize_sequence(sequence);
    let out = sketch::syncmers(&seq, k, s, t).map_err(PyValueError::new_err)?;
    Ok(out
        .into_iter()
        .map(|(pos, bytes)| (pos, String::from_utf8_lossy(&bytes).into_owned()))
        .collect())
}

#[pyfunction]
#[pyo3(signature = (kmer, m, num_buckets))]
fn minimizer_bucket(kmer: &str, m: usize, num_buckets: usize) -> PyResult<usize> {
    let seq = normalize_sequence(kmer);
    sketch::minimizer_bucket(&seq, m, num_buckets).map_err(PyValueError::new_err)
}

#[pyfunction]
#[pyo3(signature = (reads, w, m, k, min_length = 0))]
fn mdbg_assemble(
    reads: Vec<String>,
    w: usize,
    m: usize,
    k: usize,
    min_length: usize,
) -> PyResult<Vec<(String, f64)>> {
    mdbg::mdbg_assemble_impl(&reads, w, m, k, min_length).map_err(PyValueError::new_err)
}

#[pymodule]
fn genome_assembly_native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(count_kmers, m)?)?;
    m.add_function(wrap_pyfunction!(build_edges, m)?)?;
    m.add_function(wrap_pyfunction!(compact_contigs, m)?)?;
    m.add_function(wrap_pyfunction!(minimizers, m)?)?;
    m.add_function(wrap_pyfunction!(syncmers, m)?)?;
    m.add_function(wrap_pyfunction!(minimizer_bucket, m)?)?;
    m.add_function(wrap_pyfunction!(mdbg_assemble, m)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts_kmers_deterministically() {
        let reads = vec!["ACGTAC".to_string()];
        let result = count_kmers_impl(&reads, 3, true, 1).unwrap();
        assert_eq!(
            result,
            vec![
                ("ACG".to_string(), 1),
                ("CGT".to_string(), 1),
                ("GTA".to_string(), 1),
                ("TAC".to_string(), 1),
            ]
        );
    }

    #[test]
    fn counts_repeated_kmers() {
        let reads = vec!["AAAAA".to_string()];
        let result = count_kmers_impl(&reads, 3, true, 1).unwrap();
        assert_eq!(result, vec![("AAA".to_string(), 3)]);
    }

    #[test]
    fn skips_ambiguous_bases() {
        let reads = vec!["ACNTG".to_string()];
        let result = count_kmers_impl(&reads, 3, true, 1).unwrap();
        assert!(result.is_empty());
    }

    #[test]
    fn can_keep_ambiguous_bases() {
        let reads = vec!["ACNTG".to_string()];
        let result = count_kmers_impl(&reads, 3, false, 1).unwrap();
        assert_eq!(
            result,
            vec![
                ("ACN".to_string(), 1),
                ("CNT".to_string(), 1),
                ("NTG".to_string(), 1),
            ]
        );
    }

    #[test]
    fn rejects_zero_k() {
        let reads = vec!["ACGT".to_string()];
        let error = count_kmers_impl(&reads, 0, true, 1).unwrap_err();
        assert_eq!(error, "k must be >= 1");
    }

    #[test]
    fn builds_edge_table() {
        let reads = vec!["ACGTTA".to_string(), "GTTACC".to_string()];
        let (raw_edge_count, edges) = build_edges_impl(&reads, 3, 1, true, 1).unwrap();
        assert_eq!(raw_edge_count, 6);
        assert_eq!(
            edges,
            vec![
                ("ACG".to_string(), "CGT".to_string(), "ACGT".to_string(), 1),
                ("CGT".to_string(), "GTT".to_string(), "CGTT".to_string(), 1),
                ("GTT".to_string(), "TTA".to_string(), "GTTA".to_string(), 2),
                ("TAC".to_string(), "ACC".to_string(), "TACC".to_string(), 1),
                ("TTA".to_string(), "TAC".to_string(), "TTAC".to_string(), 1),
            ]
        );
    }

    #[test]
    fn filters_low_abundance_edges() {
        let reads = vec!["ACGTTA".to_string(), "GTTACC".to_string()];
        let (raw_edge_count, edges) = build_edges_impl(&reads, 3, 2, true, 1).unwrap();
        assert_eq!(raw_edge_count, 6);
        assert_eq!(
            edges,
            vec![("GTT".to_string(), "TTA".to_string(), "GTTA".to_string(), 2)]
        );
    }

    #[test]
    fn rejects_zero_min_abundance() {
        let reads = vec!["ACGT".to_string()];
        let error = build_edges_impl(&reads, 3, 0, true, 1).unwrap_err();
        assert_eq!(error, "min_abundance must be >= 1");
    }

    #[test]
    fn compacts_linear_contig() {
        let reads = vec!["ACGTTA".to_string(), "GTTACC".to_string()];
        let (_, edges) = build_edges_impl(&reads, 3, 1, true, 1).unwrap();
        let contigs = compact_contigs_impl(3, &edges, 0).unwrap();
        assert_eq!(contigs, vec![("ACGTTACC".to_string(), 1.2, 5)]);
    }

    #[test]
    fn compacts_circular_contig() {
        let edges = vec![
            ("AT".to_string(), "TG".to_string(), "ATG".to_string(), 1),
            ("TG".to_string(), "GA".to_string(), "TGA".to_string(), 1),
            ("GA".to_string(), "AT".to_string(), "GAT".to_string(), 1),
        ];
        let contigs = compact_contigs_impl(2, &edges, 0).unwrap();
        assert_eq!(contigs, vec![("ATGAT".to_string(), 1.0, 3)]);
    }

    #[test]
    fn filters_short_contigs() {
        let reads = vec!["ACGTTA".to_string(), "GTTACC".to_string()];
        let (_, edges) = build_edges_impl(&reads, 3, 1, true, 1).unwrap();
        let contigs = compact_contigs_impl(3, &edges, 9).unwrap();
        assert!(contigs.is_empty());
    }

    #[test]
    fn rejects_zero_node_k_for_compaction() {
        let error = compact_contigs_impl(0, &[], 0).unwrap_err();
        assert_eq!(error, "node_k must be >= 1");
    }

    #[test]
    fn parallel_counts_match_single_thread() {
        let bases = [b'A', b'C', b'G', b'T'];
        let reads: Vec<String> = (0..64)
            .map(|i| {
                (0..24)
                    .map(|j| bases[(i * 7 + j * 3) % 4] as char)
                    .collect::<String>()
            })
            .collect();
        let single = count_kmers_impl(&reads, 5, true, 1).unwrap();
        let parallel = count_kmers_impl(&reads, 5, true, 4).unwrap();
        assert_eq!(single, parallel);
        assert!(single.iter().any(|(_, count)| *count > 1));
    }

    #[test]
    fn rejects_zero_threads() {
        let reads = vec!["ACGT".to_string()];
        let error = count_kmers_impl(&reads, 3, true, 0).unwrap_err();
        assert_eq!(error, "threads must be >= 1");
    }

    #[test]
    fn uf_compaction_matches_sequential_and_is_thread_stable() {
        let bases = [b'A', b'C', b'G', b'T'];
        let reads: Vec<String> = (0..40)
            .map(|i| {
                (0..30)
                    .map(|j| bases[(i * 5 + j * 3) % 4] as char)
                    .collect::<String>()
            })
            .collect();
        let (_, edges) = build_edges_impl(&reads, 7, 1, true, 1).unwrap();
        let sequential = compact_contigs_impl(7, &edges, 0).unwrap();
        let uf_single = compact_contigs_uf_impl(7, &edges, 0, 1).unwrap();
        let uf_multi = compact_contigs_uf_impl(7, &edges, 0, 4).unwrap();
        assert!(!sequential.is_empty());
        assert_eq!(sequential, uf_single);
        assert_eq!(uf_single, uf_multi);
    }

    #[test]
    fn uf_compacts_circular_contig() {
        let edges = vec![
            ("AT".to_string(), "TG".to_string(), "ATG".to_string(), 1),
            ("TG".to_string(), "GA".to_string(), "TGA".to_string(), 1),
            ("GA".to_string(), "AT".to_string(), "GAT".to_string(), 1),
        ];
        let contigs = compact_contigs_uf_impl(2, &edges, 0, 2).unwrap();
        assert_eq!(contigs, vec![("ATGAT".to_string(), 1.0, 3)]);
    }
}
