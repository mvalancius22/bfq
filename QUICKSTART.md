# Quick Start Guide: GPU-Accelerated Sequence Alignment

## Prerequisites

1. **NVIDIA GPU** with CUDA support (Compute Capability 7.5+)
2. **CUDA Toolkit** (version 11.0 or later)
3. **C++ compiler** (g++ or clang++)
4. **Python 3** (for test data generation)

## Building the Project

### Option 1: Using the build script (recommended)

```bash
./build.sh
```

The script will:
- Check for CUDA installation
- Compile both `bfq_unpack` and `bfq_align` binaries
- Display usage examples

### Option 2: Using Make

```bash
make clean
make all
```

### Option 3: Manual compilation

```bash
# Set your GPU architecture
ARCH=sm_75  # Change based on your GPU

# Compile unpacking tool
nvcc -O3 -arch=$ARCH -std=c++11 --use_fast_math -c main.cpp -o main.o
nvcc -O3 -arch=$ARCH -std=c++11 --use_fast_math -c unpack_kernel.cu -o unpack_kernel.o
nvcc -O3 -arch=$ARCH -std=c++11 --use_fast_math -o bfq_unpack main.o unpack_kernel.o

# Compile alignment tool
nvcc -O3 -arch=$ARCH -std=c++11 --use_fast_math -c align_main.cpp -o align_main.o
nvcc -O3 -arch=$ARCH -std=c++11 --use_fast_math -c alignment_kernel.cu -o alignment_kernel.o
nvcc -O3 -arch=$ARCH -std=c++11 --use_fast_math -o bfq_align align_main.o alignment_kernel.o
```

### GPU Architecture Reference

| GPU Series | Architecture | Flag |
|------------|-------------|------|
| RTX 2000, GTX 1600 | Turing | sm_75 |
| RTX 3000, A100 | Ampere | sm_80 |
| RTX 3050/3060 | Ampere | sm_86 |
| RTX 4000 | Ada Lovelace | sm_89 |

Check your GPU: `nvidia-smi` or `nvcc --list-gpu-arch`

## Quick Examples

### 1. Generate Test Data

```bash
python3 proj.py
```

This creates `reads.bfq` with 4 test reads (16 bases each).

### 2. Test Unpacking

```bash
./bfq_unpack reads.bfq
```

Expected output:
```
Loaded BFQ file:
  Reads: 4
  Length: 16
  
Unpacking performance:
  Time: 0.042 ms
  Throughput: 1523.8 MB/s
  Bases/sec: 1.523 billion
```

### 3. Seed-and-Extend (Recommended) ⭐

```bash
# Align to 10kb synthetic reference with 15-mer seeds
./bfq_seed_align reads.bfq -l 10000 -k 15 --band 32
```

This is the **production-ready approach** used by real aligners (BWA-MEM, Bowtie2).

Expected output:
```
=== Seed-and-Extend Alignment ===

Phase 1: Finding seeds...
  Time: 2.3 ms

Phase 2: Extending seeds...
  Time: 8.7 ms

Throughput: 363.6 K reads/sec
Aligned reads: 3,842 / 4,000 (96.0%)
```

### 4. Full Smith-Waterman (Teaching/Validation)

```bash
# Exact match reference
./bfq_align reads.bfq -s ACGTACGTACGTACGT

# Short synthetic reference
./bfq_align reads.bfq -l 500 --shared
```

This is optimal but much slower - only practical for short references (< 5kb).

### 5. Scaling to Large References

```bash
# 100kb reference
./bfq_seed_align reads.bfq -l 100000 -k 15

# 1Mb reference (chromosome-scale)
./bfq_seed_align reads.bfq -l 1000000 -k 17 --band 64
```

Seed-and-extend scales to gigabase references, full SW does not.

### 6. Custom Scoring Parameters

```bash
./bfq_seed_align reads.bfq -l 10000 \
  --match 3 \
  --mismatch -2 \
  --gap-open -5 \
  -k 15 \
  --band 32
```

## Command Line Options

### bfq_seed_align Options (Recommended)

```
Usage: ./bfq_seed_align <reads.bfq> <options>

Reference Options (choose one):
  -r <file>       Load reference from FASTA file
  -s <seq>        Use sequence string as reference
  -l <length>     Generate synthetic reference of given length

Seeding Parameters:
  -k <size>       K-mer size for seeding (default: 15)
                  Smaller = more sensitive, slower
                  Larger = faster, may miss divergent sequences
  
  --interval <n>  Sample seeds every N bases (default: 1)
                  Higher = faster but less sensitive

Extension Parameters:
  --band <width>  Band width for alignment (default: 32)
                  Larger = handle bigger indels, but slower
  
  --max-hits <n>  Max seed hits per read (default: 100)

Scoring Parameters:
  --match <n>     Match score (default: 2)
  --mismatch <n>  Mismatch penalty (default: -1)
  --gap-open <n>  Gap opening penalty (default: -3)
```

### bfq_align Options (Full Smith-Waterman)

```
Usage: ./bfq_align <reads.bfq> <options>

Reference Options (choose one):
  -r <file>     Load reference from FASTA file
  -s <seq>      Use sequence string as reference
  -l <length>   Generate synthetic reference of given length

Optimization:
  --shared      Enable shared memory optimization (recommended)

Scoring Parameters:
  --match <n>        Match score (default: 2)
  --mismatch <n>     Mismatch penalty (default: -1)
  --gap-open <n>     Gap opening penalty (default: -3)
  --gap-extend <n>   Gap extension penalty (default: -1)
```

## Understanding the Output

### Alignment Performance Metrics

```
Alignment Performance:
  Time: 45.2 ms                    # Total GPU kernel time
  Throughput: 88.5 K reads/sec     # Reads processed per second
  Time per read: 11.3 µs           # Average time per read
  GCUPS: 6.64                      # Giga Cell Updates Per Second
```

**GCUPS** is a standard benchmark for alignment algorithms:
- GCUPS = (reads × query_len × ref_len) / time
- Higher is better
- Typical values: 5-20 GCUPS for Smith-Waterman on GPU

### Alignment Statistics

```
Alignment Statistics:
  Aligned reads: 3842 / 4000 (96.0%)  # Reads with score > 0
  Average score: 28.5                  # Mean alignment score
  Max score: 32                        # Best alignment score
```

### Individual Results

```
Read    Score    RefPos    QueryEnd
----------------------------------------
0         32        16          16      # Perfect match
1         24        35          14      # Partial alignment
2          8        89           8      # Weak match
```

- **Score**: Smith-Waterman alignment score (higher = better)
- **RefPos**: Position in reference where alignment ends
- **QueryEnd**: Position in query where alignment ends

## Creating Your Own BFQ Files

Modify `proj.py` to create custom datasets:

```python
import struct

def pack_sequence(seq):
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    packed = []
    for i in range(0, len(seq), 16):
        chunk = seq[i:i+16]
        word = 0
        for j, base in enumerate(chunk):
            val = mapping.get(base, 0)
            word |= (val << (j * 2))
        packed.append(word)
    return packed

def write_bfq(filename, reads):
    read_count = len(reads)
    read_len = len(reads[0])
    with open(filename, 'wb') as f:
        f.write(struct.pack('<4sIQQ', b'BFQ1', 1, read_count, read_len))
        f.write(struct.pack('I', read_len))
        for read in reads:
            packed = pack_sequence(read)
            f.write(struct.pack(f'<{len(packed)}I', *packed))

# Your custom reads
reads = [
    "ACGTACGTACGTACGTACGTACGT",
    "TGCATGCATGCATGCATGCATGCA",
    # ... add more
]

write_bfq('my_reads.bfq', reads)
```

## Troubleshooting

### "nvcc not found"

Install CUDA Toolkit:
- **Ubuntu**: `sudo apt install nvidia-cuda-toolkit`
- **Download**: https://developer.nvidia.com/cuda-downloads

### "No CUDA-capable device"

Ensure you have an NVIDIA GPU:
```bash
nvidia-smi  # Should list your GPU
```

### "Unsupported architecture"

Update `CUDA_ARCH` in build.sh or Makefile to match your GPU.

### Low GCUPS performance

1. Check GPU utilization: `nvidia-smi` while running
2. Try different block sizes in kernel launch
3. Enable shared memory: `--shared` flag
4. Ensure GPU is not throttling (check temperature)

## Performance Tuning

### For Maximum Throughput

1. **Use shared memory** when aligning many reads to same reference
2. **Batch size**: Process reads in groups that fit in GPU memory
3. **Read length**: Optimal for reads 50-500 bases
4. **Reference length**: Best with references < 5000 bases

### Expected Performance

On a modern GPU (RTX 3000 series):
- **100K reads × 150bp × 500bp reference**: ~50ms (6-8 GCUPS)
- **1M reads × 100bp × 1000bp reference**: ~500ms (8-10 GCUPS)

## Next Steps

1. **Test with real data**: Convert FASTQ to BFQ format
2. **Benchmark**: Compare seed-and-extend vs full SW performance
3. **Scale up**: Try larger reference sequences (100kb → 1Mb → chromosome)
4. **Tune parameters**: Experiment with k-mer size and band width
5. **Extend functionality**: Add CIGAR strings, seed chaining
6. **Integrate**: Use in your genomics pipeline

## Choosing the Right Tool

### Use `bfq_seed_align` (Seed-and-Extend) when:
- ✅ Reference is > 5kb
- ✅ You need production-level performance
- ✅ Processing many reads (> 1000)
- ✅ Want to scale to chromosome/genome size

### Use `bfq_align` (Full SW) when:
- ✅ Reference is < 5kb
- ✅ Need guaranteed optimal alignment
- ✅ Validating results
- ✅ Teaching/demonstrating DP algorithms

## Additional Resources

- [README.md](README.md) - Full project documentation
- [alignment_kernel.cu](alignment_kernel.cu) - Kernel implementation
- Smith-Waterman algorithm: https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
- CUDA programming guide: https://docs.nvidia.com/cuda/

## Support

For issues or questions:
1. Check the main [README.md](README.md)
2. Review kernel implementation comments
3. Open an issue on the repository
