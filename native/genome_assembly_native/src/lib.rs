use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::collections::BTreeMap;

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

#[pymodule]
fn genome_assembly_native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(count_kmers, m)?)?;
    m.add_function(wrap_pyfunction!(build_edges, m)?)?;
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
}
