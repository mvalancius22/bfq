#include "bfq_reader.hpp"
#include "reference_loader.hpp"
#include "scoring.hpp"
#include "alignment_result.hpp"
#include <cuda_runtime.h>
#include <iostream>
#include <iomanip>
#include <cstring>

// Forward declaration of CUDA alignment function
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

void print_usage(const char* prog_name) {
    std::cerr << "GPU-Accelerated Smith-Waterman Sequence Alignment\n\n";
    std::cerr << "Usage: " << prog_name << " <reads.bfq> <options>\n\n";
    std::cerr << "Options:\n";
    std::cerr << "  -r <file>    Reference FASTA file\n";
    std::cerr << "  -s <seq>     Reference sequence string\n";
    std::cerr << "  -l <length>  Generate synthetic reference of given length\n";
    std::cerr << "  --shared     Use shared memory optimization\n";
    std::cerr << "  --match <n>  Match score (default: 2)\n";
    std::cerr << "  --mismatch <n> Mismatch penalty (default: -1)\n";
    std::cerr << "  --gap-open <n> Gap open penalty (default: -3)\n";
    std::cerr << "  --gap-extend <n> Gap extend penalty (default: -1)\n";
    std::cerr << "\nExample:\n";
    std::cerr << "  " << prog_name << " reads.bfq -s ACGTACGTACGTACGT\n";
    std::cerr << "  " << prog_name << " reads.bfq -l 1000 --shared\n";
}

int main(int argc, char** argv) {
    if (argc < 3) {
        print_usage(argv[0]);
        return 1;
    }
    
    // Parse arguments
    std::string bfq_file = argv[1];
    std::string ref_fasta;
    std::string ref_sequence;
    uint64_t synthetic_length = 0;
    bool use_shared = false;
    
    AlignmentScoring scoring;  // Default values
    
    for (int i = 2; i < argc; i++) {
        std::string arg = argv[i];
        
        if (arg == "-r" && i + 1 < argc) {
            ref_fasta = argv[++i];
        } else if (arg == "-s" && i + 1 < argc) {
            ref_sequence = argv[++i];
        } else if (arg == "-l" && i + 1 < argc) {
            synthetic_length = std::stoull(argv[++i]);
        } else if (arg == "--shared") {
            use_shared = true;
        } else if (arg == "--match" && i + 1 < argc) {
            scoring.match_score = std::stoi(argv[++i]);
        } else if (arg == "--mismatch" && i + 1 < argc) {
            scoring.mismatch_penalty = std::stoi(argv[++i]);
        } else if (arg == "--gap-open" && i + 1 < argc) {
            scoring.gap_open_penalty = std::stoi(argv[++i]);
        } else if (arg == "--gap-extend" && i + 1 < argc) {
            scoring.gap_extend_penalty = std::stoi(argv[++i]);
        }
    }
    
    // Load query reads
    std::cout << "Loading query reads from: " << bfq_file << std::endl;
    BFQReader reader;
    if (!reader.load(bfq_file)) {
        std::cerr << "Failed to load " << bfq_file << std::endl;
        return 1;
    }
    
    uint64_t read_count = reader.get_read_count();
    uint64_t query_len = reader.get_read_length();
    size_t words_per_query = reader.get_words_per_read();
    
    std::cout << "  Reads: " << read_count << std::endl;
    std::cout << "  Length: " << query_len << " bases" << std::endl;
    
    // Load reference
    std::cout << "\nLoading reference sequence..." << std::endl;
    ReferenceLoader ref_loader;
    
    if (!ref_fasta.empty()) {
        if (!ref_loader.load_fasta(ref_fasta)) {
            std::cerr << "Failed to load reference from: " << ref_fasta << std::endl;
            return 1;
        }
        std::cout << "  Loaded from FASTA: " << ref_fasta << std::endl;
    } else if (!ref_sequence.empty()) {
        if (!ref_loader.load_string(ref_sequence)) {
            std::cerr << "Failed to load reference sequence" << std::endl;
            return 1;
        }
        std::cout << "  Loaded from string" << std::endl;
    } else if (synthetic_length > 0) {
        ref_loader.create_synthetic(synthetic_length);
        std::cout << "  Generated synthetic sequence" << std::endl;
    } else {
        std::cerr << "No reference provided. Use -r, -s, or -l option." << std::endl;
        return 1;
    }
    
    uint64_t ref_len = ref_loader.get_length();
    size_t words_per_ref = ref_loader.get_words_per_sequence();
    
    std::cout << "  Reference length: " << ref_len << " bases" << std::endl;
    std::cout << "  Words per ref: " << words_per_ref << std::endl;
    
    // Display scoring parameters
    std::cout << "\nScoring parameters:" << std::endl;
    std::cout << "  Match: +" << (int)scoring.match_score << std::endl;
    std::cout << "  Mismatch: " << (int)scoring.mismatch_penalty << std::endl;
    std::cout << "  Gap open: " << (int)scoring.gap_open_penalty << std::endl;
    std::cout << "  Gap extend: " << (int)scoring.gap_extend_penalty << std::endl;
    std::cout << "  Shared memory: " << (use_shared ? "enabled" : "disabled") << std::endl;
    
    // Allocate device memory
    uint32_t* d_query_reads;
    uint32_t* d_ref_seq;
    CompactAlignment* d_results;
    
    size_t query_size = read_count * words_per_query * sizeof(uint32_t);
    size_t ref_size = words_per_ref * sizeof(uint32_t);
    size_t results_size = read_count * sizeof(CompactAlignment);
    
    std::cout << "\nAllocating GPU memory..." << std::endl;
    std::cout << "  Query reads: " << query_size / (1024.0 * 1024.0) << " MB" << std::endl;
    std::cout << "  Reference: " << ref_size / 1024.0 << " KB" << std::endl;
    std::cout << "  Results: " << results_size / 1024.0 << " KB" << std::endl;
    
    CUDA_CHECK(cudaMalloc(&d_query_reads, query_size));
    CUDA_CHECK(cudaMalloc(&d_ref_seq, ref_size));
    CUDA_CHECK(cudaMalloc(&d_results, results_size));
    
    // Copy data to GPU
    std::cout << "\nTransferring data to GPU..." << std::endl;
    CUDA_CHECK(cudaMemcpy(d_query_reads, reader.get_data(), query_size, 
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_ref_seq, ref_loader.get_data(), ref_size, 
                          cudaMemcpyHostToDevice));
    
    // Run alignment
    std::cout << "\nRunning Smith-Waterman alignment on GPU..." << std::endl;
    
    cudaEvent_t start, stop;
    CUDA_CHECK(cudaEventCreate(&start));
    CUDA_CHECK(cudaEventCreate(&stop));
    
    CUDA_CHECK(cudaEventRecord(start));
    
    align_reads_cuda(
        d_query_reads,
        d_ref_seq,
        d_results,
        query_len,
        ref_len,
        read_count,
        words_per_query,
        words_per_ref,
        scoring,
        use_shared,
        0  // default stream
    );
    
    CUDA_CHECK(cudaEventRecord(stop));
    CUDA_CHECK(cudaEventSynchronize(stop));
    
    float milliseconds = 0;
    CUDA_CHECK(cudaEventElapsedTime(&milliseconds, start, stop));
    
    // Copy results back
    std::vector<CompactAlignment> results(read_count);
    CUDA_CHECK(cudaMemcpy(results.data(), d_results, results_size, 
                          cudaMemcpyDeviceToHost));
    
    // Analyze results
    std::cout << "\nAlignment Performance:" << std::endl;
    std::cout << "  Time: " << milliseconds << " ms" << std::endl;
    std::cout << "  Throughput: " << (read_count / 1e3) / (milliseconds / 1000.0) 
              << " K reads/sec" << std::endl;
    std::cout << "  Time per read: " << (milliseconds * 1000.0) / read_count 
              << " µs" << std::endl;
    
    uint64_t total_cells = (uint64_t)read_count * query_len * ref_len;
    double gcups = (total_cells / 1e9) / (milliseconds / 1000.0);
    std::cout << "  GCUPS: " << gcups << " (billion cell updates per second)" << std::endl;
    
    // Display alignment statistics
    int64_t total_score = 0;
    int aligned_count = 0;
    int32_t max_score = 0;
    
    for (size_t i = 0; i < read_count; i++) {
        total_score += results[i].score;
        if (results[i].score > 0) {
            aligned_count++;
        }
        if (results[i].score > max_score) {
            max_score = results[i].score;
        }
    }
    
    std::cout << "\nAlignment Statistics:" << std::endl;
    std::cout << "  Aligned reads: " << aligned_count << " / " << read_count 
              << " (" << std::fixed << std::setprecision(1) 
              << (100.0 * aligned_count / read_count) << "%)" << std::endl;
    std::cout << "  Average score: " << (double)total_score / read_count << std::endl;
    std::cout << "  Max score: " << max_score << std::endl;
    
    // Display first few alignments
    std::cout << "\nFirst 10 alignment results:" << std::endl;
    std::cout << std::setw(8) << "Read" 
              << std::setw(10) << "Score" 
              << std::setw(12) << "RefPos"
              << std::setw(12) << "QueryEnd" << std::endl;
    std::cout << std::string(42, '-') << std::endl;
    
    for (size_t i = 0; i < std::min(10UL, read_count); i++) {
        std::cout << std::setw(8) << i 
                  << std::setw(10) << results[i].score 
                  << std::setw(12) << results[i].ref_pos
                  << std::setw(12) << results[i].query_end << std::endl;
    }
    
    // Cleanup
    CUDA_CHECK(cudaFree(d_query_reads));
    CUDA_CHECK(cudaFree(d_ref_seq));
    CUDA_CHECK(cudaFree(d_results));
    CUDA_CHECK(cudaEventDestroy(start));
    CUDA_CHECK(cudaEventDestroy(stop));
    
    std::cout << "\nAlignment complete!" << std::endl;
    
    return 0;
}
