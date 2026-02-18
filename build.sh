#!/bin/bash
# Build script for BFQ GPU alignment tools

set -e

# Configuration
CUDA_ARCH=${CUDA_ARCH:-sm_75}  # Default to Turing, adjust as needed
# sm_75 = Turing (RTX 2000, GTX 1600)
# sm_80 = Ampere (RTX 3000, A100)
# sm_86 = Ampere (RTX 3050/3060)
# sm_89 = Ada Lovelace (RTX 4000)

echo "Building BFQ tools with CUDA architecture: $CUDA_ARCH"

# Check for nvcc
if ! command -v nvcc &> /dev/null; then
    echo "Error: nvcc not found. Please install CUDA toolkit."
    exit 1
fi

# Build flags
NVCC_FLAGS="-O3 -arch=$CUDA_ARCH -std=c++11 --use_fast_math"

# Clean previous builds
echo "Cleaning previous builds..."
rm -f bfq_unpack bfq_align bfq_seed_align *.o

# Build unpacking tool
echo "Building bfq_unpack..."
nvcc $NVCC_FLAGS -c main.cpp -o main.o
nvcc $NVCC_FLAGS -c unpack_kernel.cu -o unpack_kernel.o
nvcc $NVCC_FLAGS -o bfq_unpack main.o unpack_kernel.o

# Build alignment tool
echo "Building bfq_align..."
nvcc $NVCC_FLAGS -c align_main.cpp -o align_main.o
nvcc $NVCC_FLAGS -c alignment_kernel.cu -o alignment_kernel.o
nvcc $NVCC_FLAGS -o bfq_align align_main.o alignment_kernel.o

# Build seed-and-extend alignment tool
echo "Building bfq_seed_align..."
nvcc $NVCC_FLAGS -c seed_align_main.cpp -o seed_align_main.o
nvcc $NVCC_FLAGS -c seed_extend_kernel.cu -o seed_extend_kernel.o
nvcc $NVCC_FLAGS -o bfq_seed_align seed_align_main.o seed_extend_kernel.o

echo "Build complete!"
echo ""
echo "Binaries created:"
echo "  - bfq_unpack:      Unpack and benchmark BFQ format"
echo "  - bfq_align:       Full Smith-Waterman alignment (slow, optimal)"
echo "  - bfq_seed_align:  Seed-and-extend alignment (fast, like BWA-MEM)"
echo ""
echo "Usage examples:"
echo "  ./bfq_unpack reads.bfq"
echo "  ./bfq_align reads.bfq -s ACGTACGTACGTACGT"
echo "  ./bfq_seed_align reads.bfq -l 10000 -k 15 --band 32"
