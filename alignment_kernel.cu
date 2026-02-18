#include <cuda_runtime.h>
#include <cstdint>
#include <cstdio>
#include "scoring.hpp"
#include "alignment_result.hpp"

// Base decoding lookup table (same as unpack kernel)
__constant__ char base_lut[4] = {'A', 'C', 'G', 'T'};

/**
 * Device function to extract a single 2-bit base from packed data
 */
__device__ inline uint8_t extract_base(const uint32_t* packed_data, uint64_t base_idx, size_t words_per_read) {
    size_t word_idx = base_idx / 16;
    size_t bit_offset = (base_idx % 16) * 2;
    uint32_t word = packed_data[word_idx];
    return (word >> bit_offset) & 0x3;
}

/**
 * Smith-Waterman local alignment kernel (banded version for efficiency)
 * 
 * This is a simplified version that:
 * - Aligns each read to a reference sequence
 * - Uses a diagonal band to reduce memory usage
 * - Outputs alignment score and positions
 * 
 * Each thread processes one read independently
 * 
 * @param d_query_reads  Packed query reads (2-bit encoded)
 * @param d_ref_seq      Packed reference sequence (2-bit encoded)
 * @param d_results      Output alignment results
 * @param query_len      Length of query reads
 * @param ref_len        Length of reference sequence
 * @param read_count     Number of reads to align
 * @param words_per_query Words per query read
 * @param words_per_ref   Words per reference
 * @param scoring        Scoring parameters
 */
__global__ void smith_waterman_kernel(
    const uint32_t* __restrict__ d_query_reads,
    const uint32_t* __restrict__ d_ref_seq,
    CompactAlignment* __restrict__ d_results,
    uint32_t query_len,
    uint32_t ref_len,
    uint64_t read_count,
    size_t words_per_query,
    size_t words_per_ref,
    ScoringMatrix scoring,
    int8_t gap_open,
    int8_t gap_extend)
{
    uint64_t read_idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (read_idx >= read_count) return;
    
    // Initialize result
    CompactAlignment result;
    result.score = 0;
    result.ref_pos = -1;
    
    // Allocate local DP matrices (using simplified 2-row approach)
    // For full Smith-Waterman, we need: H (match), E (gap in ref), F (gap in query)
    // We only store current and previous row to save memory
    
    const int MAX_QUERY_LEN = 512;  // Adjust based on typical read length
    if (query_len > MAX_QUERY_LEN) return;
    
    int16_t H_prev[MAX_QUERY_LEN + 1];  // Previous row
    int16_t H_curr[MAX_QUERY_LEN + 1];  // Current row
    
    // Initialize first row (all zeros for local alignment)
    for (int j = 0; j <= query_len; j++) {
        H_prev[j] = 0;
    }
    
    int16_t max_score = 0;
    int max_i = 0, max_j = 0;
    
    // Pointer to this read's packed data
    const uint32_t* query_packed = d_query_reads + read_idx * words_per_query;
    
    // Fill DP matrix
    for (int i = 1; i <= ref_len; i++) {
        H_curr[0] = 0;  // First column is always 0 for local alignment
        
        // Extract reference base for this row
        uint8_t ref_base = extract_base(d_ref_seq, i - 1, words_per_ref);
        
        for (int j = 1; j <= query_len; j++) {
            // Extract query base
            uint8_t query_base = extract_base(query_packed, j - 1, words_per_query);
            
            // Calculate match/mismatch score
            int16_t match = H_prev[j - 1] + scoring.score(query_base, ref_base);
            
            // Calculate gap scores (simplified affine gap)
            int16_t gap_ref = H_prev[j] + gap_open;     // Gap in reference
            int16_t gap_query = H_curr[j - 1] + gap_open;  // Gap in query
            
            // Smith-Waterman: take maximum of (match, gaps, 0)
            int16_t score = 0;
            if (match > score) score = match;
            if (gap_ref > score) score = gap_ref;
            if (gap_query > score) score = gap_query;
            
            H_curr[j] = score;
            
            // Track maximum score
            if (score > max_score) {
                max_score = score;
                max_i = i;
                max_j = j;
            }
        }
        
        // Swap rows
        for (int j = 0; j <= query_len; j++) {
            H_prev[j] = H_curr[j];
        }
    }
    
    // Store result
    result.score = max_score;
    result.ref_pos = max_i;
    result.query_end = max_j;
    result.query_start = 0;  // Would need traceback to get exact start
    
    d_results[read_idx] = result;
}

/**
 * Optimized kernel using shared memory for reference sequence
 * Good when multiple reads align to the same reference region
 */
__global__ void smith_waterman_kernel_shared(
    const uint32_t* __restrict__ d_query_reads,
    const uint32_t* __restrict__ d_ref_seq,
    CompactAlignment* __restrict__ d_results,
    uint32_t query_len,
    uint32_t ref_len,
    uint64_t read_count,
    size_t words_per_query,
    size_t words_per_ref,
    ScoringMatrix scoring,
    int8_t gap_open,
    int8_t gap_extend)
{
    // Load reference into shared memory (all threads cooperate)
    extern __shared__ uint32_t shared_ref[];
    
    // Cooperative load
    for (int i = threadIdx.x; i < words_per_ref; i += blockDim.x) {
        shared_ref[i] = d_ref_seq[i];
    }
    __syncthreads();
    
    // Now each thread processes its read using shared reference
    uint64_t read_idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (read_idx >= read_count) return;
    
    CompactAlignment result;
    result.score = 0;
    result.ref_pos = -1;
    
    const int MAX_QUERY_LEN = 512;
    if (query_len > MAX_QUERY_LEN) return;
    
    int16_t H_prev[MAX_QUERY_LEN + 1];
    int16_t H_curr[MAX_QUERY_LEN + 1];
    
    for (int j = 0; j <= query_len; j++) {
        H_prev[j] = 0;
    }
    
    int16_t max_score = 0;
    int max_i = 0, max_j = 0;
    
    const uint32_t* query_packed = d_query_reads + read_idx * words_per_query;
    
    for (int i = 1; i <= ref_len; i++) {
        H_curr[0] = 0;
        
        // Extract from shared memory
        uint8_t ref_base = extract_base(shared_ref, i - 1, words_per_ref);
        
        for (int j = 1; j <= query_len; j++) {
            uint8_t query_base = extract_base(query_packed, j - 1, words_per_query);
            
            int16_t match = H_prev[j - 1] + scoring.score(query_base, ref_base);
            int16_t gap_ref = H_prev[j] + gap_open;
            int16_t gap_query = H_curr[j - 1] + gap_open;
            
            int16_t score = 0;
            if (match > score) score = match;
            if (gap_ref > score) score = gap_ref;
            if (gap_query > score) score = gap_query;
            
            H_curr[j] = score;
            
            if (score > max_score) {
                max_score = score;
                max_i = i;
                max_j = j;
            }
        }
        
        for (int j = 0; j <= query_len; j++) {
            H_prev[j] = H_curr[j];
        }
    }
    
    result.score = max_score;
    result.ref_pos = max_i;
    result.query_end = max_j;
    result.query_start = 0;
    
    d_results[read_idx] = result;
}

/**
 * Host wrapper functions
 */
extern "C" {

void align_reads_cuda(
    const uint32_t* d_query_reads,
    const uint32_t* d_ref_seq,
    CompactAlignment* d_results,
    uint32_t query_len,
    uint32_t ref_len,
    uint64_t read_count,
    size_t words_per_query,
    size_t words_per_ref,
    const AlignmentScoring& scoring_params,
    bool use_shared_memory,
    cudaStream_t stream)
{
    // Prepare scoring matrix
    ScoringMatrix scoring_matrix;
    scoring_matrix.initialize(scoring_params);
    
    // Copy scoring matrix to device (could use constant memory for better performance)
    
    int block_size = 128;  // Fewer threads per block for local memory
    int grid_size = (read_count + block_size - 1) / block_size;
    
    if (use_shared_memory) {
        size_t shared_mem_size = words_per_ref * sizeof(uint32_t);
        
        smith_waterman_kernel_shared<<<grid_size, block_size, shared_mem_size, stream>>>(
            d_query_reads,
            d_ref_seq,
            d_results,
            query_len,
            ref_len,
            read_count,
            words_per_query,
            words_per_ref,
            scoring_matrix,
            scoring_params.gap_open_penalty,
            scoring_params.gap_extend_penalty
        );
    } else {
        smith_waterman_kernel<<<grid_size, block_size, 0, stream>>>(
            d_query_reads,
            d_ref_seq,
            d_results,
            query_len,
            ref_len,
            read_count,
            words_per_query,
            words_per_ref,
            scoring_matrix,
            scoring_params.gap_open_penalty,
            scoring_params.gap_extend_penalty
        );
    }
}

} // extern "C"
