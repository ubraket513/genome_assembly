use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::collections::{BTreeMap, BTreeSet};

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

pub fn count_kmers_impl(
    reads: &[String],
    k: usize,
    skip_ambiguous: bool,
) -> Result<Vec<(String, usize)>, String> {
    if k == 0 {
        return Err("k must be >= 1".to_string());
    }

    let mut counts: BTreeMap<String, usize> = BTreeMap::new();
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
            let kmer = String::from_utf8(window.to_vec())
                .map_err(|_| "sequence contained invalid UTF-8 after normalization".to_string())?;
            *counts.entry(kmer).or_insert(0) += 1;
        }
    }

    Ok(counts.into_iter().collect())
}

pub type EdgeRow = (String, String, String, usize);
pub type ContigRow = (String, f64, usize);

pub fn build_edges_impl(
    reads: &[String],
    node_k: usize,
    min_abundance: usize,
    skip_ambiguous: bool,
) -> Result<(usize, Vec<EdgeRow>), String> {
    if node_k == 0 {
        return Err("node_k must be >= 1".to_string());
    }
    if min_abundance == 0 {
        return Err("min_abundance must be >= 1".to_string());
    }

    let edge_k = node_k + 1;
    let counts = count_kmers_impl(reads, edge_k, skip_ambiguous)?;
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
            .then_with(|| left.0.cmp(&right.0))
    });

    Ok(contigs
        .into_iter()
        .map(|(_, sequence, mean_abundance, edge_count)| (sequence, mean_abundance, edge_count))
        .collect())
}

#[pyfunction]
#[pyo3(signature = (reads, k, skip_ambiguous = true))]
fn count_kmers(
    reads: Vec<String>,
    k: usize,
    skip_ambiguous: bool,
) -> PyResult<Vec<(String, usize)>> {
    count_kmers_impl(&reads, k, skip_ambiguous).map_err(PyValueError::new_err)
}

#[pyfunction]
#[pyo3(signature = (reads, node_k, min_abundance = 1, skip_ambiguous = true))]
fn build_edges(
    reads: Vec<String>,
    node_k: usize,
    min_abundance: usize,
    skip_ambiguous: bool,
) -> PyResult<(usize, Vec<EdgeRow>)> {
    build_edges_impl(&reads, node_k, min_abundance, skip_ambiguous).map_err(PyValueError::new_err)
}

#[pyfunction]
#[pyo3(signature = (node_k, edge_rows, min_length = 0))]
fn compact_contigs(
    node_k: usize,
    edge_rows: Vec<EdgeRow>,
    min_length: usize,
) -> PyResult<Vec<ContigRow>> {
    compact_contigs_impl(node_k, &edge_rows, min_length).map_err(PyValueError::new_err)
}

#[pymodule]
fn genome_assembly_native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(count_kmers, m)?)?;
    m.add_function(wrap_pyfunction!(build_edges, m)?)?;
    m.add_function(wrap_pyfunction!(compact_contigs, m)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts_kmers_deterministically() {
        let reads = vec!["ACGTAC".to_string()];
        let result = count_kmers_impl(&reads, 3, true).unwrap();
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
        let result = count_kmers_impl(&reads, 3, true).unwrap();
        assert_eq!(result, vec![("AAA".to_string(), 3)]);
    }

    #[test]
    fn skips_ambiguous_bases() {
        let reads = vec!["ACNTG".to_string()];
        let result = count_kmers_impl(&reads, 3, true).unwrap();
        assert!(result.is_empty());
    }

    #[test]
    fn can_keep_ambiguous_bases() {
        let reads = vec!["ACNTG".to_string()];
        let result = count_kmers_impl(&reads, 3, false).unwrap();
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
        let error = count_kmers_impl(&reads, 0, true).unwrap_err();
        assert_eq!(error, "k must be >= 1");
    }

    #[test]
    fn builds_edge_table() {
        let reads = vec!["ACGTTA".to_string(), "GTTACC".to_string()];
        let (raw_edge_count, edges) = build_edges_impl(&reads, 3, 1, true).unwrap();
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
        let (raw_edge_count, edges) = build_edges_impl(&reads, 3, 2, true).unwrap();
        assert_eq!(raw_edge_count, 6);
        assert_eq!(
            edges,
            vec![("GTT".to_string(), "TTA".to_string(), "GTTA".to_string(), 2)]
        );
    }

    #[test]
    fn rejects_zero_min_abundance() {
        let reads = vec!["ACGT".to_string()];
        let error = build_edges_impl(&reads, 3, 0, true).unwrap_err();
        assert_eq!(error, "min_abundance must be >= 1");
    }

    #[test]
    fn compacts_linear_contig() {
        let reads = vec!["ACGTTA".to_string(), "GTTACC".to_string()];
        let (_, edges) = build_edges_impl(&reads, 3, 1, true).unwrap();
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
        let (_, edges) = build_edges_impl(&reads, 3, 1, true).unwrap();
        let contigs = compact_contigs_impl(3, &edges, 9).unwrap();
        assert!(contigs.is_empty());
    }

    #[test]
    fn rejects_zero_node_k_for_compaction() {
        let error = compact_contigs_impl(0, &[], 0).unwrap_err();
        assert_eq!(error, "node_k must be >= 1");
    }
}
