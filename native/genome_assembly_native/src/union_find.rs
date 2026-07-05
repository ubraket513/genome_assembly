//! Disjoint-set (Union-Find) with path compression and union by rank.
//!
//! Used to glue de Bruijn edges into unitigs, BCALM2-style: unioning the in-
//! and out-edge at every one-in-one-out node collapses each maximal
//! non-branching path into a single component.

pub struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<u32>,
}

impl UnionFind {
    pub fn new(n: usize) -> Self {
        UnionFind {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }

    pub fn find(&mut self, x: usize) -> usize {
        let mut root = x;
        while self.parent[root] != root {
            root = self.parent[root];
        }
        // Path compression: point every node on the path straight at the root.
        let mut current = x;
        while self.parent[current] != root {
            let next = self.parent[current];
            self.parent[current] = root;
            current = next;
        }
        root
    }

    pub fn union(&mut self, a: usize, b: usize) {
        let ra = self.find(a);
        let rb = self.find(b);
        if ra == rb {
            return;
        }
        match self.rank[ra].cmp(&self.rank[rb]) {
            std::cmp::Ordering::Less => self.parent[ra] = rb,
            std::cmp::Ordering::Greater => self.parent[rb] = ra,
            std::cmp::Ordering::Equal => {
                self.parent[rb] = ra;
                self.rank[ra] += 1;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn merges_components() {
        let mut uf = UnionFind::new(5);
        uf.union(0, 1);
        uf.union(1, 2);
        assert_eq!(uf.find(0), uf.find(2));
        assert_ne!(uf.find(0), uf.find(3));
        uf.union(3, 4);
        uf.union(2, 4);
        assert_eq!(uf.find(0), uf.find(3));
    }

    #[test]
    fn singletons_are_distinct() {
        let mut uf = UnionFind::new(3);
        assert_ne!(uf.find(0), uf.find(1));
        assert_ne!(uf.find(1), uf.find(2));
    }
}
