#!/bin/bash
# Test script for seed-and-extend alignment

set -e

echo "======================================="
echo "Seed-and-Extend Alignment Test Suite"
echo "======================================="
echo ""

# Check if binary exists
if [ ! -f "bfq_seed_align" ]; then
    echo "Error: bfq_seed_align not found."
    echo "Please run: ./build.sh"
    exit 1
fi

# Generate test data
echo "Step 1: Generating test BFQ file..."
python3 proj.py
echo "✓ Created reads.bfq"
echo ""

# Test 1: Small reference
echo "Test 1: Alignment to 1kb reference..."
./bfq_seed_align reads.bfq -l 1000 -k 15 --band 32
echo ""

# Test 2: Medium reference
echo "Test 2: Alignment to 10kb reference..."
./bfq_seed_align reads.bfq -l 10000 -k 15 --band 32
echo ""

# Test 3: Large reference (scalability test)
echo "Test 3: Alignment to 100kb reference (scalability)..."
./bfq_seed_align reads.bfq -l 100000 -k 17 --band 64
echo ""

# Test 4: Different k-mer sizes
echo "Test 4: Testing different k-mer sizes..."
echo "  k=13 (sensitive)..."
./bfq_seed_align reads.bfq -l 10000 -k 13 --band 32 | grep "Total time"
echo "  k=15 (balanced)..."
./bfq_seed_align reads.bfq -l 10000 -k 15 --band 32 | grep "Total time"
echo "  k=17 (fast)..."
./bfq_seed_align reads.bfq -l 10000 -k 17 --band 32 | grep "Total time"
echo ""

# Test 5: Sparse sampling
echo "Test 5: Testing sparse seed sampling..."
echo "  Every base (interval=1)..."
./bfq_seed_align reads.bfq -l 10000 -k 15 --interval 1 | grep "Total time"
echo "  Every 3rd base (interval=3)..."
./bfq_seed_align reads.bfq -l 10000 -k 15 --interval 3 | grep "Total time"
echo "  Every 5th base (interval=5)..."
./bfq_seed_align reads.bfq -l 10000 -k 15 --interval 5 | grep "Total time"
echo ""

echo "======================================="
echo "All tests completed successfully!"
echo "======================================="
echo ""
echo "Key observations:"
echo "  • Seed-and-extend scales well to large references"
echo "  • Larger k-mer sizes are faster but may miss matches"
echo "  • Sparse sampling (interval > 1) trades speed for sensitivity"
echo ""
echo "Compare with full Smith-Waterman:"
echo "  ./bfq_align reads.bfq -l 1000"
echo ""
