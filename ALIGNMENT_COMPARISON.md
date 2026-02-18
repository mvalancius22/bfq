# Alignment Strategies: Full Smith-Waterman vs Seed-and-Extend

## Overview

This project implements **two alignment approaches** to demonstrate the evolution from academic algorithms to production tools.

## Full Smith-Waterman (bfq_align)

### What It Is
Traditional dynamic programming that computes the optimal alignment by filling a complete scoring matrix.

### Algorithm Complexity
- **Time**: O(n × m) where n = query length, m = reference length
- **Space**: O(n) with 2-row optimization (originally O(n × m))

### When It's Used in Production
- **Short exact alignment** (< 100bp against < 5kb)
- **Validation and testing** (gold standard for accuracy)
- **Final polish step** after seed-and-extend
- **Specific domains** where accuracy is critical (antibody design, etc.)

### Limitations
```
For a 150bp read vs 1Mb reference:
- Matrix cells to compute: 150 × 1,000,000 = 150 million cells
- GPU time at 10 GCUPS: ~15ms per read
- For 1M reads: ~4 hours on GPU

This is too slow for real genomics workflows!
```

### Example Use Case
```bash
# Aligning to a short gene region (< 5kb)
./bfq_align reads.bfq -r gene.fasta --shared

# Validating seed-and-extend results
./bfq_align subset.bfq -s ACGTACGT...
```

---

## Seed-and-Extend (bfq_seed_align)

### What It Is
Modern heuristic approach used by **BWA-MEM, Bowtie2, minimap2** - the tools that actually power genomics pipelines.

### Algorithm Complexity
- **Indexing**: O(m) to build k-mer hash table (one-time cost)
- **Seeding**: O(n/k × log(|kmers|)) to find seeds
- **Extension**: O(n × w) where w = band width (typically 32-64)
- **Total**: Much faster than O(n × m) when m >> w

### Two-Phase Strategy

#### Phase 1: Seeding (Fast)
```
Goal: Find regions where query and reference might align

1. Extract k-mers from query (k=15 typical)
   Read: ACGTACGTACGTACGTACGT
   K-mers: ACGTACGTACGTACG, CGTACGTACGTACGT, ...
   
2. Look up each k-mer in reference hash table
   - Hash table built once from reference
   - O(1) lookup per k-mer
   - Returns all positions where k-mer occurs
   
3. Output: List of "seed hits"
   Seed: (query_pos=0, ref_pos=12345, diagonal=12345)
   Seed: (query_pos=5, ref_pos=12350, diagonal=12345)
   ...
```

**Why this is fast:**
- Only looking up small k-mers (15bp)
- No alignment computation yet
- Skips most of the reference (only checks k-mer matches)

#### Phase 2: Extension (Focused)
```
Goal: Verify and extend promising seed hits

1. For each seed, run BANDED Smith-Waterman
   - Only compute cells near the diagonal
   - Band width = ±32bp typically
   - Handles indels and mismatches
   
2. Matrix cells computed:
   Full SW: 150 × 1,000,000 = 150M cells
   Banded:  150 × 64 = 9,600 cells per seed
   
3. If you have 10 seeds: 10 × 9,600 = 96K cells
   Speedup: 1562x fewer cells!
```

### Performance Comparison

| Scenario | Full SW | Seed-and-Extend | Speedup |
|----------|---------|-----------------|---------|
| 150bp vs 5kb | 0.75M cells | ~10K cells (5 seeds) | **75x** |
| 150bp vs 100kb | 15M cells | ~100K cells (10 seeds) | **150x** |
| 150bp vs 1Mb | 150M cells | ~100K cells (10 seeds) | **1500x** |
| 150bp vs 3Gb (human) | 450B cells | ~200K cells (20 seeds) | **2.25Mx** |

### Why Production Tools Use This

1. **Scales to large references** (human genome = 3 billion bases)
2. **Near-optimal accuracy** (misses < 1% of alignments)
3. **Configurable speed/sensitivity tradeoff**
4. **Handles real-world data** (sequencing errors, variants)

### Example Use Cases

```bash
# Standard alignment to chromosome
./bfq_seed_align reads.bfq -r chr1.fasta -k 15 --band 32

# Large reference (full genome)
./bfq_seed_align reads.bfq -r hg38.fasta -k 17 --band 64

# Fast mode (sacrifice sensitivity for speed)
./bfq_seed_align reads.bfq -r ref.fasta -k 20 --interval 5 --band 32

# Sensitive mode (more seeds, slower)
./bfq_seed_align reads.bfq -r ref.fasta -k 13 --interval 1 --band 64
```

---

## Parameter Tuning Guide

### K-mer Size (-k)

| Value | Effect | Use Case |
|-------|--------|----------|
| k=11-13 | More seeds, sensitive | High error rate data |
| **k=15** | **Balanced (default)** | **Standard Illumina** |
| k=17-20 | Fewer seeds, faster | Long reads, low errors |

**Rule of thumb:** 
- Smaller k = more sensitive, slower
- Larger k = faster, may miss divergent sequences
- Expected hits per k-mer ≈ 4^k / genome_size

### Band Width (--band)

| Value | Effect | Use Case |
|-------|--------|----------|
| 16-24 | Narrow, fast | Low indel rate |
| **32** | **Balanced (default)** | **Standard genomics** |
| 48-64 | Wide, slower | High indel rate, SVs |

**Rule of thumb:**
- Band width ≥ 2 × expected indel size
- Wider band = can handle larger indels, but slower

### Seed Interval (--interval)

| Value | Effect |
|-------|--------|
| 1 | Sample every base (default, sensitive) |
| 2-5 | Sample every N bases (faster, less sensitive) |
| 10+ | Sparse sampling (very fast, may miss alignments) |

---

## When to Use Each Approach

### Use Full Smith-Waterman When:
- ✅ Reference is short (< 5kb)
- ✅ Need guaranteed optimal alignment
- ✅ Validating results
- ✅ Teaching/demonstrating DP algorithms
- ✅ Have time and compute resources

### Use Seed-and-Extend When:
- ✅ **Reference is large** (> 100kb, especially > 1Mb)
- ✅ **Production genomics pipeline**
- ✅ Need to process millions of reads
- ✅ Want results in reasonable time
- ✅ Near-optimal is acceptable

---

## Real-World Example

Aligning 1 million reads (150bp) to human chromosome 1 (249Mb):

### Full Smith-Waterman
```
Cells per read: 150 × 249,000,000 = 37.4 billion
Total cells: 37.4 trillion
Time at 10 GCUPS: ~1040 hours = 43 days
```

### Seed-and-Extend
```
Seeds per read: ~20 (with k=15)
Cells per seed: 150 × 64 = 9,600
Total cells per read: 192,000
Total cells (1M reads): 192 billion
Time at 10 GCUPS: ~19 seconds for extension + seeding overhead
Realistic total: ~30-60 seconds
```

**Result: ~25,000x speedup** 🚀

---

## Production Aligner Comparison

| Tool | Strategy | Speed (1M reads) |
|------|----------|------------------|
| BWA-MEM | Seed + extend + chain | ~5 min |
| Bowtie2 | FM-index + extend | ~3 min |
| minimap2 | Minimizer + chain | ~1 min |
| **This project** | K-mer + banded SW | ~1 min (estimated) |
| Naive SW | Full DP | ~43 days |

---

## Implementation Notes

### What This Project Demonstrates

1. **Educational Value**
   - Shows why full SW isn't practical for large-scale genomics
   - Demonstrates the algorithmic evolution
   - Clear comparison of O(n×m) vs O(n×w)

2. **Production Techniques**
   - K-mer indexing (hash table on GPU)
   - Banded alignment (focus computation where it matters)
   - Seed filtering (ignore repetitive k-mers)
   - Parallel processing (one thread per read)

3. **GPU Optimization**
   - Memory coalescing (threads access adjacent memory)
   - Efficient data structures (packed 2-bit format)
   - Minimal data transfer (pre-build index)

### What's Still Missing (for full production readiness)

- **Seed chaining**: Link seeds on same diagonal
- **CIGAR strings**: Detailed alignment description
- **Secondary alignments**: Report multiple mapping locations
- **Paired-end support**: Handle read pairs
- **Quality-aware scoring**: Use Phred scores
- **BWA-style FM-index**: Even faster seed finding

---

## Conclusion

**Full Smith-Waterman** is the gold standard for alignment accuracy, but it's computationally intractable for real genomics workflows with large references.

**Seed-and-extend** is the practical solution used by all modern production aligners. It sacrifices < 1% accuracy for 100-1000x speedup, making large-scale genomics possible.

This project demonstrates **both approaches** to show:
1. Why the naive approach doesn't scale
2. How production tools actually work
3. The tradeoffs between optimality and speed

For real work: **Use seed-and-extend** (bfq_seed_align) 🎯
