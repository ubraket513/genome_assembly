1. Generate random reads of coverage 3x from given coronavirus FASTA file
2. Store this collection of reads into some pkl file
3. Refer to this pkl file and construct de bruijn graph

## More on Step 3
- How to construct a compact de bruijn graph?
- Let's not think about parallelization.
- Consider the case when we only have one read.

### Procedural approach to de bruijn graph construction
If our de bruijn graph has vertices with k mers, then we start by...
1. Generate all (k+1)-mers
2. For each generate (k+1)-mer, we take its leftmost k-mer and rightmost k-mer.
3. We try to add these two k-mers as vertices to the de bruijn graph.
- If it already exists, then we do not add it.
- If it does not exist, then we add it.
4. After adding them, connect two vertices by a directed edge.
5. we continue this step for every (k+1)-mer we have.
6. We fully constructed the de bruijn graph.

### Generalization to handle multiple reads
1. Suppose that we have de bruijn graphs for each read.
2. We want to find their union.
3. 
