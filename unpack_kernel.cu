#include <cuda_runtime.h>
#include <cstdint>
#include <cstdio>

// Base mapping: 0=A, 1=C, 2=G, 3=T
__constant__ char base_lut[4] = {'A', 'C', 'G', 'T'};

/**
 * CUDA kernel for unpacking 2-bit encoded bases from uint32 words
 * 
 * Memory Layout (coalesced access):
 * - Each uint32 contains 16 bases (2 bits each)
 * - Threads in a warp access consecutive uint32 words
 * - Each thread extracts all 16 bases from its word
 * 
 * @param packed_reads  Input: packed uint32 array
 * @param unpacked_bases Output: ASCII base array
 * @param read_count    Number of reads
 * @param read_len      Length of each read in bases
 * @param words_per_read Number of uint32 words per read
 */
__global__ void unpack_bases_kernel(
    const uint32_t* __restrict__ packed_reads,
    char* __restrict__ unpacked_bases,
    uint64_t read_count,
    uint64_t read_len,
    size_t words_per_read)
{
    // Global thread index
    uint64_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Calculate which read and which word within that read
    uint64_t read_idx = tid / words_per_read;
    uint64_t word_idx = tid % words_per_read;
    
    if (read_idx >= read_count) return;
    
    // Coalesced read: consecutive threads read consecutive uint32s
    uint32_t packed_word = packed_reads[tid];
    
    // Calculate output position for this word's bases
    uint64_t base_offset = read_idx * read_len + word_idx * 16;
    
    // Extract all 16 bases from this uint32 word
    #pragma unroll
    for (int i = 0; i < 16; i++) {
        uint64_t base_pos = base_offset + i;
        if (base_pos < read_idx * read_len + read_len) {
            // Extract 2 bits
            uint32_t base_code = (packed_word >> (i * 2)) & 0x3;
            // Convert to ASCII using LUT
            unpacked_bases[base_pos] = base_lut[base_code];
        }
    }
}

/**
 * Optimized kernel using shared memory for better bank access patterns
 * Good for when you need to process bases within a block
 */
__global__ void unpack_bases_kernel_shared(
    const uint32_t* __restrict__ packed_reads,
    char* __restrict__ unpacked_bases,
    uint64_t read_count,
    uint64_t read_len,
    size_t words_per_read)
{
    extern __shared__ uint32_t shared_packed[];
    
    uint64_t block_start = blockIdx.x * blockDim.x;
    uint64_t tid = threadIdx.x;
    uint64_t global_tid = block_start + tid;
    
    // Coalesced load into shared memory
    if (global_tid < read_count * words_per_read) {
        shared_packed[tid] = packed_reads[global_tid];
    }
    __syncthreads();
    
    if (global_tid >= read_count * words_per_read) return;
    
    uint64_t read_idx = global_tid / words_per_read;
    uint64_t word_idx = global_tid % words_per_read;
    
    uint32_t packed_word = shared_packed[tid];
    uint64_t base_offset = read_idx * read_len + word_idx * 16;
    
    #pragma unroll
    for (int i = 0; i < 16; i++) {
        uint64_t base_pos = base_offset + i;
        if (base_pos < read_idx * read_len + read_len) {
            uint32_t base_code = (packed_word >> (i * 2)) & 0x3;
            unpacked_bases[base_pos] = base_lut[base_code];
        }
    }
}

/**
 * Ultra-optimized kernel for extracting a single base from each read
 * Perfect for column-wise operations (e.g., getting base at position N from all reads)
 * 
 * @param packed_reads   Input: packed uint32 array
 * @param base_column    Output: single base per read
 * @param base_position  Which base position to extract (0-indexed)
 * @param read_count     Number of reads
 * @param words_per_read Number of uint32 words per read
 */
__global__ void extract_single_base_kernel(
    const uint32_t* __restrict__ packed_reads,
    char* __restrict__ base_column,
    uint64_t base_position,
    uint64_t read_count,
    size_t words_per_read)
{
    uint64_t read_idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (read_idx >= read_count) return;
    
    // Calculate which word and which position within word
    size_t word_idx = base_position / 16;
    size_t bit_offset = (base_position % 16) * 2;
    
    // Coalesced access: each thread reads different read, same word position
    uint32_t packed_word = packed_reads[read_idx * words_per_read + word_idx];
    
    // Extract the specific base
    uint32_t base_code = (packed_word >> bit_offset) & 0x3;
    base_column[read_idx] = base_lut[base_code];
}

/**
 * Host wrapper function for unpacking
 */
extern "C" {
    
void unpack_bases_cuda(
    const uint32_t* d_packed_reads,
    char* d_unpacked_bases,
    uint64_t read_count,
    uint64_t read_len,
    size_t words_per_read,
    cudaStream_t stream = 0)
{
    size_t total_words = read_count * words_per_read;
    
    // Configure launch parameters
    int block_size = 256;
    int grid_size = (total_words + block_size - 1) / block_size;
    
    unpack_bases_kernel<<<grid_size, block_size, 0, stream>>>(
        d_packed_reads,
        d_unpacked_bases,
        read_count,
        read_len,
        words_per_read
    );
}

void unpack_bases_cuda_shared(
    const uint32_t* d_packed_reads,
    char* d_unpacked_bases,
    uint64_t read_count,
    uint64_t read_len,
    size_t words_per_read,
    cudaStream_t stream = 0)
{
    size_t total_words = read_count * words_per_read;
    
    int block_size = 256;
    int grid_size = (total_words + block_size - 1) / block_size;
    size_t shared_mem_size = block_size * sizeof(uint32_t);
    
    unpack_bases_kernel_shared<<<grid_size, block_size, shared_mem_size, stream>>>(
        d_packed_reads,
        d_unpacked_bases,
        read_count,
        read_len,
        words_per_read
    );
}

void extract_base_column_cuda(
    const uint32_t* d_packed_reads,
    char* d_base_column,
    uint64_t base_position,
    uint64_t read_count,
    size_t words_per_read,
    cudaStream_t stream = 0)
{
    int block_size = 256;
    int grid_size = (read_count + block_size - 1) / block_size;
    
    extract_single_base_kernel<<<grid_size, block_size, 0, stream>>>(
        d_packed_reads,
        d_base_column,
        base_position,
        read_count,
        words_per_read
    );
}

} // extern "C"
