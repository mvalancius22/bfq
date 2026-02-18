#!/bin/bash
# Test script for GPU alignment functionality

set -e

echo "====================================="
echo "BFQ GPU Alignment Test Suite"
echo "====================================="
echo ""

# Check if binaries exist
if [ ! -f "bfq_align" ]; then
    echo "Error: bfq_align binary not found."
    echo "Please run: ./build.sh"
    exit 1
fi

# Generate test data
echo "Step 1: Generating test BFQ file..."
python3 proj.py
echo "✓ Created reads.bfq"
echo ""

# Test 1: Simple alignment to exact match
echo "Step 2: Testing alignment to exact reference..."
./bfq_align reads.bfq -s ACGTACGTACGTACGT
echo ""

# Test 2: Alignment with synthetic reference
echo "Step 3: Testing alignment to synthetic reference (500bp)..."
./bfq_align reads.bfq -l 500
echo ""

# Test 3: Alignment with shared memory
echo "Step 4: Testing shared memory optimization..."
./bfq_align reads.bfq -l 500 --shared
echo ""

# Test 4: Custom scoring parameters
echo "Step 5: Testing custom scoring parameters..."
./bfq_align reads.bfq -l 500 --match 3 --mismatch -2 --gap-open -5
echo ""

echo "====================================="
echo "All tests completed successfully!"
echo "====================================="
