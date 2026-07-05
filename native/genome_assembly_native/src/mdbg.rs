//! Minimizer-space de Bruijn graph (mdBG) assembly.
//!
//! Reads are reduced to ordered lists of minimizer tokens; the de Bruijn graph
//! is built over k-min-mers (k consecutive minimizers) exactly as a base-space
//! dBG is built over (k+1)-mers. Each edge carries the base span it covers, so
//! contigs are reconstructed back into nucleotide sequence. This is the design
//! behind rust-mdbg: 10-100x smaller graphs for long, accurate reads.

use crate::sketch;
use std::collections::BTreeMap;

const SEP: char = '\u{1}'; // token separator, never appears in DNA

pub type MdbgContig = (String, f64); // (sequence, mean_abundance)

struct MdbgEdge {
    span: String,
    prefix_span_len: usize, // base length of the (k-1)-min-mer prefix node
    count: usize,
}

fn join(tokens: &[&str]) -> String {
    tokens.join(&SEP.to_string())
}

/// Build the minimizer-space edge table from reads.
///
/// Returns edges keyed by (prefix-node, suffix-node), where nodes are
/// (k-1)-min-mers rendered as separator-joined minimizer strings.
fn build_mdbg_edges(
    reads: &[String],
    w: usize,
    m: usize,
    k: usize,
) -> Result<BTreeMap<(String, String), MdbgEdge>, String> {
    if k < 2 {
        return Err("k (min-mers per k-min-mer) must be >= 2".to_string());
    }

    let mut edges: BTreeMap<(String, String), MdbgEdge> = BTreeMap::new();

    for read in reads {
        let sequence: Vec<u8> = read
            .bytes()
            .filter(|b| !b.is_ascii_whitespace())
            .map(|b| b.to_ascii_uppercase())
            .collect();
        let mins = sketch::minimizers(&sequence, w, m)?;
        if mins.len() < k {
            continue;
        }
        let positions: Vec<usize> = mins.iter().map(|(pos, _)| *pos).collect();
        let tokens: Vec<String> = mins
            .iter()
            .map(|(_, bytes)| String::from_utf8_lossy(bytes).into_owned())
            .collect();
        let token_refs: Vec<&str> = tokens.iter().map(String::as_str).collect();

        // Windows of k tokens = k-min-mers (edges); nodes are k-1 tokens.
        for i in 0..=(tokens.len() - k) {
            let prefix = join(&token_refs[i..i + k - 1]);
            let suffix = join(&token_refs[i + 1..i + k]);
            let span_start = positions[i];
            let span_end = positions[i + k - 1] + m;
            let span = String::from_utf8_lossy(&sequence[span_start..span_end]).into_owned();
            let prefix_span_len = (positions[i + k - 2] + m) - positions[i];

            edges
                .entry((prefix, suffix))
                .and_modify(|edge| edge.count += 1)
                .or_insert(MdbgEdge {
                    span,
                    prefix_span_len,
                    count: 1,
                });
        }
    }

    Ok(edges)
}

/// Assemble reads in minimizer space and reconstruct nucleotide contigs.
pub fn mdbg_assemble_impl(
    reads: &[String],
    w: usize,
    m: usize,
    k: usize,
    min_length: usize,
) -> Result<Vec<MdbgContig>, String> {
    let edges = build_mdbg_edges(reads, w, m, k)?;
    if edges.is_empty() {
        return Ok(Vec::new());
    }

    // Adjacency over minimizer-space nodes.
    let mut out_edges: BTreeMap<String, Vec<(String, String)>> = BTreeMap::new();
    let mut in_degree: BTreeMap<String, usize> = BTreeMap::new();
    let mut out_degree: BTreeMap<String, usize> = BTreeMap::new();
    for (prefix, suffix) in edges.keys() {
        out_edges
            .entry(prefix.clone())
            .or_default()
            .push((prefix.clone(), suffix.clone()));
        *out_degree.entry(prefix.clone()).or_insert(0) += 1;
        *in_degree.entry(suffix.clone()).or_insert(0) += 1;
    }

    let is_simple = |node: &str| -> bool {
        in_degree.get(node).copied().unwrap_or(0) == 1
            && out_degree.get(node).copied().unwrap_or(0) == 1
    };

    let mut visited: std::collections::BTreeSet<(String, String)> =
        std::collections::BTreeSet::new();
    let mut contigs: Vec<MdbgContig> = Vec::new();

    let walk = |start: &(String, String),
                visited: &mut std::collections::BTreeSet<(String, String)>|
     -> MdbgContig {
        visited.insert(start.clone());
        let start_edge = &edges[start];
        let mut sequence = start_edge.span.clone();
        let mut count_sum = start_edge.count;
        let mut edge_count = 1usize;
        let mut current = start.1.clone();

        while is_simple(&current) {
            let Some(nexts) = out_edges.get(&current) else {
                break;
            };
            let next_key = nexts[0].clone();
            if visited.contains(&next_key) {
                break;
            }
            visited.insert(next_key.clone());
            let next_edge = &edges[&next_key];
            // Append only the bases beyond the shared (k-1)-min-mer overlap.
            sequence.push_str(&next_edge.span[next_edge.prefix_span_len..]);
            count_sum += next_edge.count;
            current = next_key.1.clone();
            edge_count += 1;
        }

        (sequence, count_sum as f64 / edge_count as f64)
    };

    // Start unitigs at non-simple prefix nodes, then sweep remaining cycles.
    let mut start_keys: Vec<(String, String)> = edges.keys().cloned().collect();
    start_keys.sort();
    for key in &start_keys {
        if is_simple(&key.0) || visited.contains(key) {
            continue;
        }
        let contig = walk(key, &mut visited);
        if contig.0.len() >= min_length {
            contigs.push(contig);
        }
    }
    for key in &start_keys {
        if visited.contains(key) {
            continue;
        }
        let contig = walk(key, &mut visited);
        if contig.0.len() >= min_length {
            contigs.push(contig);
        }
    }

    contigs.sort_by(|left, right| {
        right
            .0
            .len()
            .cmp(&left.0.len())
            .then_with(|| left.0.cmp(&right.0))
    });
    Ok(contigs)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn random_genome(length: usize) -> String {
        // xorshift64 for a deterministic, reference-free genome.
        let bases = [b'A', b'C', b'G', b'T'];
        let mut x: u64 = 0x9E37_79B9_7F4A_7C15;
        (0..length)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                bases[(x & 3) as usize] as char
            })
            .collect()
    }

    #[test]
    fn contigs_are_substrings_of_source_genome() {
        // Correctness: every reconstructed contig must be an exact substring of
        // the genome the reads came from. A reconstruction bug would splice bases
        // that never appear contiguously in the source.
        let genome = random_genome(800);
        let reads: Vec<String> = (0..genome.len().saturating_sub(120))
            .step_by(40)
            .map(|start| genome[start..start + 120].to_string())
            .collect();

        let contigs = mdbg_assemble_impl(&reads, 5, 8, 3, 0).unwrap();
        assert!(!contigs.is_empty());
        for (sequence, _) in &contigs {
            assert!(
                genome.contains(sequence.as_str()),
                "contig is not a substring of the source genome"
            );
        }
    }

    #[test]
    fn assembly_is_deterministic() {
        let genome = random_genome(500);
        let reads: Vec<String> = (0..genome.len().saturating_sub(100))
            .step_by(30)
            .map(|start| genome[start..start + 100].to_string())
            .collect();
        let first = mdbg_assemble_impl(&reads, 5, 8, 3, 0).unwrap();
        let second = mdbg_assemble_impl(&reads, 5, 8, 3, 0).unwrap();
        assert_eq!(first, second);
    }

    #[test]
    fn rejects_small_k() {
        let err =
            mdbg_assemble_impl(std::slice::from_ref(&"ACGT".to_string()), 3, 2, 1, 0).unwrap_err();
        assert!(err.contains("k"));
    }
}
