#include "bfq_reader.hpp"
#include <cuda_runtime.h>
#include <iostream>
#include <cstring>

// Forward declarations of CUDA functions
extern "C" {
    void unpack_bases_cuda(
        const uint32_t* d_packed_reads,
        char* d_unpacked_bases,
        uint64_t read_count,
        uint64_t read_len,
        size_t words_per_read,
        cudaStream_t stream);
    
    void extract_base_column_cuda(
        const uint32_t* d_packed_reads,
        char* d_base_column,
        uint64_t base_position,
        uint64_t read_count,
        size_t words_per_read,
        cudaStream_t stream);
}

#define CUDA_CHECK(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            std::cerr << "CUDA error in " << __FILE__ << ":" << __LINE__ \
                      << " - " << cudaGetErrorString(err) << std::endl; \
            exit(EXIT_FAILURE); \
        } \
    } while(0)

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <bfq_file>" << std::endl;
        return 1;
    }
    
    // Load BFQ file
    BFQReader reader;
    if (!reader.load(argv[1])) {
        std::cerr << "Failed to load " << argv[1] << std::endl;
        return 1;
    }
    
    std::cout << "Loaded BFQ file:" << std::endl;
    std::cout << "  Reads: " << reader.get_read_count() << std::endl;
    std::cout << "  Length: " << reader.get_read_length() << std::endl;
    std::cout << "  Words per read: " << reader.get_words_per_read() << std::endl;
    
    uint64_t read_count = reader.get_read_count();
    uint64_t read_len = reader.get_read_length();
    size_t words_per_read = reader.get_words_per_read();
    
    // Allocate device memory
    uint32_t* d_packed;
    char* d_unpacked;
    size_t packed_size = read_count * words_per_read * sizeof(uint32_t);
    size_t unpacked_size = read_count * read_len * sizeof(char);
    
    CUDA_CHECK(cudaMalloc(&d_packed, packed_size));
    CUDA_CHECK(cudaMalloc(&d_unpacked, unpacked_size));
    
    // Copy packed data to GPU
    CUDA_CHECK(cudaMemcpy(d_packed, reader.get_data(), packed_size, 
                          cudaMemcpyHostToDevice));
    
    // Benchmark: Create events for timing
    cudaEvent_t start, stop;
    CUDA_CHECK(cudaEventCreate(&start));
    CUDA_CHECK(cudaEventCreate(&stop));
    
    // Time the unpacking kernel
    CUDA_CHECK(cudaEventRecord(start));
    unpack_bases_cuda(d_packed, d_unpacked, read_count, read_len, 
                      words_per_read, 0);
    CUDA_CHECK(cudaEventRecord(stop));
    CUDA_CHECK(cudaEventSynchronize(stop));
    
    float milliseconds = 0;
    CUDA_CHECK(cudaEventElapsedTime(&milliseconds, start, stop));
    
    std::cout << "\nUnpacking performance:" << std::endl;
    std::cout << "  Time: " << milliseconds << " ms" << std::endl;
    std::cout << "  Throughput: " 
              << (unpacked_size / 1e6) / (milliseconds / 1000.0) 
              << " MB/s" << std::endl;
    std::cout << "  Bases/sec: " 
              << (read_count * read_len / 1e9) / (milliseconds / 1000.0)
              << " billion" << std::endl;
    
    // Copy results back and display first few reads
    char* unpacked = new char[unpacked_size];
    CUDA_CHECK(cudaMemcpy(unpacked, d_unpacked, unpacked_size, 
                          cudaMemcpyDeviceToHost));
    
    std::cout << "\nFirst 4 unpacked reads:" << std::endl;
    for (int i = 0; i < std::min(4UL, read_count); i++) {
        std::cout << "Read " << i << ": ";
        for (uint64_t j = 0; j < read_len; j++) {
            std::cout << unpacked[i * read_len + j];
        }
        std::cout << std::endl;
    }
    
    // Test column extraction
    char* d_column;
    CUDA_CHECK(cudaMalloc(&d_column, read_count));
    
    std::cout << "\nExtracting base at position 5 from all reads:" << std::endl;
    extract_base_column_cuda(d_packed, d_column, 5, read_count, 
                            words_per_read, 0);
    
    char* column = new char[read_count];
    CUDA_CHECK(cudaMemcpy(column, d_column, read_count, cudaMemcpyDeviceToHost));
    
    for (size_t i = 0; i < std::min(10UL, read_count); i++) {
        std::cout << "Read " << i << " base[5]: " << column[i] << std::endl;
    }
    
    // Cleanup
    delete[] unpacked;
    delete[] column;
    CUDA_CHECK(cudaFree(d_packed));
    CUDA_CHECK(cudaFree(d_unpacked));
    CUDA_CHECK(cudaFree(d_column));
    CUDA_CHECK(cudaEventDestroy(start));
    CUDA_CHECK(cudaEventDestroy(stop));
    
    return 0;
}
