#include <cuda_runtime.h>
#include <cstdint>
#include <cstdio>
#include "scoring.hpp"
#include "alignment_result.hpp"
#include "seed_index.hpp"

/**
 * Device function to extract a single base from packed data
 */
__device__ inline uint8_t extract_base_device(
    const uint32_t* packed_data, 
    uint64_t base_idx, 
    size_t words_per_seq)
{
    size_t word_idx = base_idx / 16;
    size_t bit_offset = (base_idx % 16) * 2;
    uint32_t word = packed_data[word_idx];
    return (word >> bit_offset) & 0x3;
}

/**
 * Device function to extract k-mer hash from packed sequence
 */
__device__ inline uint32_t extract_kmer_hash_device(
    const uint32_t* packed_data,
    uint64_t start_pos,
    int k,
    size_t words_per_seq)
{
    uint32_t hash = 0;
    for (int i = 0; i < k; i++) {
        uint8_t base = extract_base_device(packed_data, start_pos + i, words_per_seq);
        hash = (hash << 2) | base;
    }
    return hash;
}

/**
 * Kernel to find seed matches between query reads and reference
 * 
 * Each thread processes one query read and finds all matching k-mers
 * 
 * @param d_query_reads     Packed query reads
 * @param d_kmer_hashes     Array of unique k-mer hashes in index
 * @param d_position_starts Start indices in positions array
 * @param d_position_counts Number of positions per k-mer
 * @param d_positions       Flattened array of reference positions
 * @param d_seed_hits       Output: seed hits for each read
 * @param d_hit_counts      Output: number of hits per read
 * @param query_len         Length of query reads
 * @param read_count        Number of reads
 * @param words_per_query   Words per query read
 * @param k                 K-mer length
 * @param kmer_count        Number of unique k-mers in index
 * @param max_hits_per_read Maximum hits to store per read
 * @param seed_interval     Sample seeds every N bases
 */
__global__ void find_seeds_kernel(
    const uint32_t* __restrict__ d_query_reads,
    const uint32_t* __restrict__ d_kmer_hashes,
    const uint32_t* __restrict__ d_position_starts,
    const uint32_t* __restrict__ d_position_counts,
    const uint32_t* __restrict__ d_positions,
    SeedHit* __restrict__ d_seed_hits,
    uint32_t* __restrict__ d_hit_counts,
    uint32_t query_len,
    uint64_t read_count,
    size_t words_per_query,
    int k,
    uint32_t kmer_count,
    uint32_t max_hits_per_read,
    int seed_interval)
{
    uint64_t read_idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (read_idx >= read_count) return;
    
    const uint32_t* query_packed = d_query_reads + read_idx * words_per_query;
    SeedHit* read_hits = d_seed_hits + read_idx * max_hits_per_read;
    uint32_t hit_count = 0;
    
    // Extract k-mers from query and look them up in index
    for (uint32_t qpos = 0; qpos <= query_len - k && hit_count < max_hits_per_read; qpos += seed_interval) {
        uint32_t query_kmer = extract_kmer_hash_device(query_packed, qpos, k, words_per_query);
        
        // Binary search for k-mer in sorted hash array
        // For simplicity, using linear search (could optimize with binary search)
        for (uint32_t idx = 0; idx < kmer_count; idx++) {
            if (d_kmer_hashes[idx] == query_kmer) {
                // Found matching k-mer, get all its positions
                uint32_t start = d_position_starts[idx];
                uint32_t count = d_position_counts[idx];
                
                // Add all positions as seed hits
                for (uint32_t i = 0; i < count && hit_count < max_hits_per_read; i++) {
                    uint32_t ref_pos = d_positions[start + i];
                    read_hits[hit_count++] = SeedHit(qpos, ref_pos, query_kmer);
                }
                break;
            }
        }
    }
    
    d_hit_counts[read_idx] = hit_count;
}

/**
 * Banded Smith-Waterman alignment kernel
 * 
 * Extends seed hits using a diagonal band around the seed position
 * Much faster than full Smith-Waterman (O(n*w) instead of O(n*m))
 * 
 * @param d_query_reads  Packed query reads
 * @param d_ref_seq      Packed reference sequence
 * @param d_seed_hits    Seed hits to extend
 * @param d_hit_counts   Number of hits per read
 * @param d_results      Output alignment results
 * @param query_len      Length of query reads
 * @param ref_len        Length of reference
 * @param read_count     Number of reads
 * @param words_per_query Words per query
 * @param words_per_ref   Words per reference
 * @param scoring         Scoring matrix
 * @param gap_penalty     Gap penalty
 * @param band_width      Width of diagonal band
 * @param max_hits_per_read Maximum hits per read
 */
__global__ void banded_extend_kernel(
    const uint32_t* __restrict__ d_query_reads,
    const uint32_t* __restrict__ d_ref_seq,
    const SeedHit* __restrict__ d_seed_hits,
    const uint32_t* __restrict__ d_hit_counts,
    CompactAlignment* __restrict__ d_results,
    uint32_t query_len,
    uint32_t ref_len,
    uint64_t read_count,
    size_t words_per_query,
    size_t words_per_ref,
    ScoringMatrix scoring,
    int8_t gap_penalty,
    int band_width,
    uint32_t max_hits_per_read)
{
    uint64_t read_idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (read_idx >= read_count) return;
    
    const uint32_t* query_packed = d_query_reads + read_idx * words_per_query;
    const SeedHit* read_hits = d_seed_hits + read_idx * max_hits_per_read;
    uint32_t hit_count = d_hit_counts[read_idx];
    
    CompactAlignment best_result;
    best_result.score = 0;
    
    // Extend each seed hit
    for (uint32_t h = 0; h < hit_count; h++) {
        SeedHit seed = read_hits[h];
        
        // Define extension region around seed
        int32_t ref_start = max(0, (int32_t)seed.ref_pos - band_width);
        int32_t ref_end = min((int32_t)ref_len, (int32_t)seed.ref_pos + band_width);
        int32_t query_start = 0;  // Extend from beginning of query
        int32_t query_end = query_len;
        
        // Simple banded alignment (simplified for demo)
        // In production, would use proper banded DP with optimized band tracking
        
        int16_t max_score = 0;
        int16_t current_score = 0;
        int best_ref_pos = seed.ref_pos;
        int best_query_pos = seed.query_pos;
        
        // Extend forward from seed
        int qpos = seed.query_pos;
        int rpos = seed.ref_pos;
        
        while (qpos < query_end && rpos < ref_end) {
            uint8_t qbase = extract_base_device(query_packed, qpos, words_per_query);
            uint8_t rbase = extract_base_device(d_ref_seq, rpos, words_per_ref);
            
            if (qbase == rbase) {
                current_score += scoring.score(qbase, rbase);
            } else {
                current_score += scoring.score(qbase, rbase);
                if (current_score < 0) {
                    current_score = 0;
                }
            }
            
            if (current_score > max_score) {
                max_score = current_score;
                best_ref_pos = rpos;
                best_query_pos = qpos;
            }
            
            qpos++;
            rpos++;
        }
        
        // Update best result if this seed extension is better
        if (max_score > best_result.score) {
            best_result.score = max_score;
            best_result.ref_pos = best_ref_pos;
            best_result.query_end = best_query_pos + 1;
            best_result.query_start = seed.query_pos;
        }
    }
    
    d_results[read_idx] = best_result;
}

/**
 * Host wrapper functions
 */
extern "C" {

void find_seeds_cuda(
    const uint32_t* d_query_reads,
    const uint32_t* d_kmer_hashes,
    const uint32_t* d_position_starts,
    const uint32_t* d_position_counts,
    const uint32_t* d_positions,
    SeedHit* d_seed_hits,
    uint32_t* d_hit_counts,
    uint32_t query_len,
    uint64_t read_count,
    size_t words_per_query,
    int k,
    uint32_t kmer_count,
    uint32_t max_hits_per_read,
    int seed_interval,
    cudaStream_t stream)
{
    int block_size = 256;
    int grid_size = (read_count + block_size - 1) / block_size;
    
    find_seeds_kernel<<<grid_size, block_size, 0, stream>>>(
        d_query_reads,
        d_kmer_hashes,
        d_position_starts,
        d_position_counts,
        d_positions,
        d_seed_hits,
        d_hit_counts,
        query_len,
        read_count,
        words_per_query,
        k,
        kmer_count,
        max_hits_per_read,
        seed_interval
    );
}

void banded_extend_cuda(
    const uint32_t* d_query_reads,
    const uint32_t* d_ref_seq,
    const SeedHit* d_seed_hits,
    const uint32_t* d_hit_counts,
    CompactAlignment* d_results,
    uint32_t query_len,
    uint32_t ref_len,
    uint64_t read_count,
    size_t words_per_query,
    size_t words_per_ref,
    const AlignmentScoring& scoring_params,
    int band_width,
    uint32_t max_hits_per_read,
    cudaStream_t stream)
{
    ScoringMatrix scoring_matrix;
    scoring_matrix.initialize(scoring_params);
    
    int block_size = 128;
    int grid_size = (read_count + block_size - 1) / block_size;
    
    banded_extend_kernel<<<grid_size, block_size, 0, stream>>>(
        d_query_reads,
        d_ref_seq,
        d_seed_hits,
        d_hit_counts,
        d_results,
        query_len,
        ref_len,
        read_count,
        words_per_query,
        words_per_ref,
        scoring_matrix,
        scoring_params.gap_open_penalty,
        band_width,
        max_hits_per_read
    );
}

} // extern "C"
