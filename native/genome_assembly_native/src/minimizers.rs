//! Minimizer and syncmer sketching, mirroring `genome_assembly/sketch.py`.
//!
//! The Python module is the correctness oracle; these must produce identical
//! output for the same inputs (same FNV-1a hash, same tie-breaking).

const FNV64_OFFSET: u64 = 0xCBF2_9CE4_8422_2325;
const FNV64_PRIME: u64 = 0x0000_0100_0000_01B3;

pub fn hash64(bytes: &[u8]) -> u64 {
    let mut h = FNV64_OFFSET;
    for &byte in bytes {
        h = (h ^ byte as u64).wrapping_mul(FNV64_PRIME);
    }
    h
}

/// Windowed (w, m)-minimizers as `(position, m-mer bytes)` pairs.
pub fn minimizers(sequence: &[u8], w: usize, m: usize) -> Result<Vec<(usize, Vec<u8>)>, String> {
    if w == 0 {
        return Err("w must be >= 1".to_string());
    }
    if m == 0 {
        return Err("m must be >= 1".to_string());
    }
    if sequence.len() < m {
        return Ok(Vec::new());
    }

    let mmer_count = sequence.len() - m + 1;
    let hashes: Vec<u64> = (0..mmer_count)
        .map(|i| hash64(&sequence[i..i + m]))
        .collect();

    let (window_count, w_eff) = if mmer_count < w {
        (1, mmer_count)
    } else {
        (mmer_count - w + 1, w)
    };

    let mut selected = Vec::new();
    let mut last_position: Option<usize> = None;
    for start in 0..window_count {
        let mut best = start;
        for offset in (start + 1)..(start + w_eff) {
            if hashes[offset] < hashes[best] {
                best = offset;
            }
        }
        if last_position != Some(best) {
            selected.push((best, sequence[best..best + m].to_vec()));
            last_position = Some(best);
        }
    }
    Ok(selected)
}

/// Open (k, s)-syncmers as `(position, k-mer bytes)` pairs.
pub fn syncmers(
    sequence: &[u8],
    k: usize,
    s: usize,
    t: usize,
) -> Result<Vec<(usize, Vec<u8>)>, String> {
    if k == 0 {
        return Err("k must be >= 1".to_string());
    }
    if s == 0 || s > k {
        return Err("s must satisfy 1 <= s <= k".to_string());
    }
    let inner = k - s + 1;
    if t >= inner {
        return Err("t must satisfy 0 <= t < k - s + 1".to_string());
    }
    if sequence.len() < k {
        return Ok(Vec::new());
    }

    let mut result = Vec::new();
    for start in 0..=(sequence.len() - k) {
        let kmer = &sequence[start..start + k];
        let mut best_offset = 0usize;
        let mut best_hash = hash64(&kmer[0..s]);
        for offset in 1..inner {
            let candidate = hash64(&kmer[offset..offset + s]);
            if candidate < best_hash {
                best_hash = candidate;
                best_offset = offset;
            }
        }
        if best_offset == t {
            result.push((start, kmer.to_vec()));
        }
    }
    Ok(result)
}

/// Route a k-mer to one of `num_buckets` independent streams by minimizer hash.
pub fn minimizer_bucket(kmer: &[u8], m: usize, num_buckets: usize) -> Result<usize, String> {
    if num_buckets == 0 {
        return Err("num_buckets must be >= 1".to_string());
    }
    if m == 0 || m > kmer.len() {
        return Err("m must satisfy 1 <= m <= len(kmer)".to_string());
    }
    let best = (0..=(kmer.len() - m))
        .map(|i| hash64(&kmer[i..i + m]))
        .min()
        .expect("at least one m-mer exists");
    Ok((best % num_buckets as u64) as usize)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hash_is_deterministic() {
        assert_eq!(hash64(b"ACGT"), hash64(b"ACGT"));
        assert_ne!(hash64(b"ACGT"), hash64(b"ACGA"));
    }

    #[test]
    fn minimizers_dedup_consecutive_positions() {
        let seq = b"ACGTACGTACGT";
        let mins = minimizers(seq, 4, 3).unwrap();
        // Positions must be strictly increasing (deduped).
        for pair in mins.windows(2) {
            assert!(pair[0].0 < pair[1].0);
        }
        assert!(!mins.is_empty());
    }

    #[test]
    fn short_sequence_yields_global_minimizer() {
        let seq = b"ACGTA"; // 3 3-mers, window 10 > count
        let mins = minimizers(seq, 10, 3).unwrap();
        assert_eq!(mins.len(), 1);
    }

    #[test]
    fn syncmers_offset_zero() {
        let seq = b"ACGTACGTACGT";
        let syncs = syncmers(seq, 5, 2, 0).unwrap();
        for (pos, kmer) in &syncs {
            assert_eq!(kmer.len(), 5);
            assert!(*pos + 5 <= seq.len());
        }
    }

    #[test]
    fn bucket_is_stable_and_in_range() {
        let a = minimizer_bucket(b"ACGTACGT", 3, 8).unwrap();
        let b = minimizer_bucket(b"ACGTACGT", 3, 8).unwrap();
        assert_eq!(a, b);
        assert!(a < 8);
    }

    #[test]
    fn rejects_bad_params() {
        assert!(minimizers(b"ACGT", 0, 3).is_err());
        assert!(syncmers(b"ACGT", 3, 4, 0).is_err());
        assert!(minimizer_bucket(b"ACGT", 0, 4).is_err());
    }
}
