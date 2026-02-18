#include "bfq_reader.hpp"
#include "reference_loader.hpp"
#include "scoring.hpp"
#include "alignment_result.hpp"
#include "seed_index.hpp"
#include <cuda_runtime.h>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <chrono>

// Forward declarations of CUDA functions
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
        cudaStream_t stream);
    
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
    std::cerr << "GPU-Accelerated Seed-and-Extend Sequence Alignment\n\n";
    std::cerr << "Usage: " << prog_name << " <reads.bfq> <options>\n\n";
    std::cerr << "Options:\n";
    std::cerr << "  -r <file>       Reference FASTA file\n";
    std::cerr << "  -s <seq>        Reference sequence string\n";
    std::cerr << "  -l <length>     Generate synthetic reference of given length\n";
    std::cerr << "  -k <size>       K-mer size for seeding (default: 15)\n";
    std::cerr << "  --band <width>  Band width for extension (default: 32)\n";
    std::cerr << "  --max-hits <n>  Max seed hits per read (default: 100)\n";
    std::cerr << "  --interval <n>  Sample seeds every N bases (default: 1)\n";
    std::cerr << "  --match <n>     Match score (default: 2)\n";
    std::cerr << "  --mismatch <n>  Mismatch penalty (default: -1)\n";
    std::cerr << "  --gap-open <n>  Gap open penalty (default: -3)\n";
    std::cerr << "\nThis uses seed-and-extend strategy like BWA-MEM:\n";
    std::cerr << "  1. Find exact k-mer matches (seeds)\n";
    std::cerr << "  2. Extend seeds with banded alignment\n";
    std::cerr << "  3. Report best alignment\n";
    std::cerr << "\nExample:\n";
    std::cerr << "  " << prog_name << " reads.bfq -l 10000 -k 15 --band 32\n";
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
    
    SeedExtendParams params;
    AlignmentScoring scoring;
    
    for (int i = 2; i < argc; i++) {
        std::string arg = argv[i];
        
        if (arg == "-r" && i + 1 < argc) {
            ref_fasta = argv[++i];
        } else if (arg == "-s" && i + 1 < argc) {
            ref_sequence = argv[++i];
        } else if (arg == "-l" && i + 1 < argc) {
            synthetic_length = std::stoull(argv[++i]);
        } else if (arg == "-k" && i + 1 < argc) {
            params.k = std::stoi(argv[++i]);
        } else if (arg == "--band" && i + 1 < argc) {
            params.band_width = std::stoi(argv[++i]);
        } else if (arg == "--max-hits" && i + 1 < argc) {
            params.max_seed_hits = std::stoi(argv[++i]);
        } else if (arg == "--interval" && i + 1 < argc) {
            params.seed_interval = std::stoi(argv[++i]);
        } else if (arg == "--match" && i + 1 < argc) {
            scoring.match_score = std::stoi(argv[++i]);
        } else if (arg == "--mismatch" && i + 1 < argc) {
            scoring.mismatch_penalty = std::stoi(argv[++i]);
        } else if (arg == "--gap-open" && i + 1 < argc) {
            scoring.gap_open_penalty = std::stoi(argv[++i]);
        }
    }
    
    // Load query reads
    std::cout << "=== Seed-and-Extend Alignment ===" << std::endl;
    std::cout << "\nLoading query reads: " << bfq_file << std::endl;
    
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
            std::cerr << "Failed to load: " << ref_fasta << std::endl;
            return 1;
        }
    } else if (!ref_sequence.empty()) {
        if (!ref_loader.load_string(ref_sequence)) {
            std::cerr << "Failed to load reference" << std::endl;
            return 1;
        }
    } else if (synthetic_length > 0) {
        ref_loader.create_synthetic(synthetic_length);
    } else {
        std::cerr << "No reference provided" << std::endl;
        return 1;
    }
    
    uint64_t ref_len = ref_loader.get_length();
    size_t words_per_ref = ref_loader.get_words_per_sequence();
    
    std::cout << "  Length: " << ref_len << " bases" << std::endl;
    
    // Build seed index
    std::cout << "\nBuilding seed index..." << std::endl;
    auto index_start = std::chrono::high_resolution_clock::now();
    
    SeedIndex index(params.k);
    index.build(ref_loader.get_data(), ref_len, words_per_ref);
    
    auto index_end = std::chrono::high_resolution_clock::now();
    auto index_time = std::chrono::duration_cast<std::chrono::milliseconds>(index_end - index_start).count();
    
    std::cout << "  K-mer size: " << params.k << std::endl;
    std::cout << "  Unique k-mers: " << index.get_kmer_count() << std::endl;
    std::cout << "  Total positions: " << index.get_total_positions() << std::endl;
    std::cout << "  Index time: " << index_time << " ms" << std::endl;
    
    // Display parameters
    std::cout << "\nAlignment parameters:" << std::endl;
    std::cout << "  Match: +" << (int)scoring.match_score << std::endl;
    std::cout << "  Mismatch: " << (int)scoring.mismatch_penalty << std::endl;
    std::cout << "  Gap open: " << (int)scoring.gap_open_penalty << std::endl;
    std::cout << "  Band width: " << params.band_width << std::endl;
    std::cout << "  Max hits/read: " << params.max_seed_hits << std::endl;
    std::cout << "  Seed interval: " << params.seed_interval << std::endl;
    
    // Allocate device memory
    std::cout << "\nAllocating GPU memory..." << std::endl;
    
    uint32_t* d_query_reads;
    uint32_t* d_ref_seq;
    uint32_t* d_kmer_hashes;
    uint32_t* d_position_starts;
    uint32_t* d_position_counts;
    uint32_t* d_positions;
    SeedHit* d_seed_hits;
    uint32_t* d_hit_counts;
    CompactAlignment* d_results;
    
    size_t query_size = read_count * words_per_query * sizeof(uint32_t);
    size_t ref_size = words_per_ref * sizeof(uint32_t);
    size_t index_hashes_size = index.get_kmer_count() * sizeof(uint32_t);
    size_t seed_hits_size = read_count * params.max_seed_hits * sizeof(SeedHit);
    size_t results_size = read_count * sizeof(CompactAlignment);
    
    CUDA_CHECK(cudaMalloc(&d_query_reads, query_size));
    CUDA_CHECK(cudaMalloc(&d_ref_seq, ref_size));
    CUDA_CHECK(cudaMalloc(&d_kmer_hashes, index_hashes_size));
    CUDA_CHECK(cudaMalloc(&d_position_starts, index_hashes_size));
    CUDA_CHECK(cudaMalloc(&d_position_counts, index_hashes_size));
    CUDA_CHECK(cudaMalloc(&d_positions, index.get_total_positions() * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_seed_hits, seed_hits_size));
    CUDA_CHECK(cudaMalloc(&d_hit_counts, read_count * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_results, results_size));
    
    std::cout << "  Query: " << query_size / (1024.0 * 1024.0) << " MB" << std::endl;
    std::cout << "  Seed index: " << (index_hashes_size * 3 + index.get_total_positions() * sizeof(uint32_t)) / (1024.0 * 1024.0) << " MB" << std::endl;
    std::cout << "  Seed hits buffer: " << seed_hits_size / (1024.0 * 1024.0) << " MB" << std::endl;
    
    // Transfer data to GPU
    std::cout << "\nTransferring data to GPU..." << std::endl;
    CUDA_CHECK(cudaMemcpy(d_query_reads, reader.get_data(), query_size, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_ref_seq, ref_loader.get_data(), ref_size, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_kmer_hashes, index.get_kmer_hashes(), index_hashes_size, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_position_starts, index.get_position_starts(), index_hashes_size, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_position_counts, index.get_position_counts(), index_hashes_size, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_positions, index.get_positions(), index.get_total_positions() * sizeof(uint32_t), cudaMemcpyHostToDevice));
    
    // Run seed-and-extend alignment
    cudaEvent_t start, stop;
    CUDA_CHECK(cudaEventCreate(&start));
    CUDA_CHECK(cudaEventCreate(&stop));
    
    std::cout << "\n=== Running Seed-and-Extend ===" << std::endl;
    
    // Phase 1: Find seeds
    std::cout << "\nPhase 1: Finding seeds..." << std::endl;
    CUDA_CHECK(cudaEventRecord(start));
    
    find_seeds_cuda(
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
        params.k,
        index.get_kmer_count(),
        params.max_seed_hits,
        params.seed_interval,
        0
    );
    
    CUDA_CHECK(cudaEventRecord(stop));
    CUDA_CHECK(cudaEventSynchronize(stop));
    
    float seed_time = 0;
    CUDA_CHECK(cudaEventElapsedTime(&seed_time, start, stop));
    std::cout << "  Time: " << seed_time << " ms" << std::endl;
    
    // Phase 2: Extend seeds
    std::cout << "\nPhase 2: Extending seeds with banded alignment..." << std::endl;
    CUDA_CHECK(cudaEventRecord(start));
    
    banded_extend_cuda(
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
        scoring,
        params.band_width,
        params.max_seed_hits,
        0
    );
    
    CUDA_CHECK(cudaEventRecord(stop));
    CUDA_CHECK(cudaEventSynchronize(stop));
    
    float extend_time = 0;
    CUDA_CHECK(cudaEventElapsedTime(&extend_time, start, stop));
    std::cout << "  Time: " << extend_time << " ms" << std::endl;
    
    float total_time = seed_time + extend_time;
    
    // Copy results back
    std::vector<CompactAlignment> results(read_count);
    std::vector<uint32_t> hit_counts(read_count);
    CUDA_CHECK(cudaMemcpy(results.data(), d_results, results_size, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(hit_counts.data(), d_hit_counts, read_count * sizeof(uint32_t), cudaMemcpyDeviceToHost));
    
    // Analyze results
    std::cout << "\n=== Performance Summary ===" << std::endl;
    std::cout << "  Seed finding: " << seed_time << " ms" << std::endl;
    std::cout << "  Seed extension: " << extend_time << " ms" << std::endl;
    std::cout << "  Total time: " << total_time << " ms" << std::endl;
    std::cout << "  Throughput: " << (read_count / 1e3) / (total_time / 1000.0) << " K reads/sec" << std::endl;
    std::cout << "  Time per read: " << (total_time * 1000.0) / read_count << " µs" << std::endl;
    
    // Seed statistics
    uint64_t total_seeds = 0;
    uint32_t max_seeds = 0;
    for (uint32_t count : hit_counts) {
        total_seeds += count;
        if (count > max_seeds) max_seeds = count;
    }
    
    std::cout << "\n=== Seed Statistics ===" << std::endl;
    std::cout << "  Total seeds found: " << total_seeds << std::endl;
    std::cout << "  Avg seeds/read: " << (double)total_seeds / read_count << std::endl;
    std::cout << "  Max seeds/read: " << max_seeds << std::endl;
    
    // Alignment statistics
    int64_t total_score = 0;
    int aligned_count = 0;
    int32_t max_score = 0;
    
    for (const auto& result : results) {
        total_score += result.score;
        if (result.score > 0) aligned_count++;
        if (result.score > max_score) max_score = result.score;
    }
    
    std::cout << "\n=== Alignment Results ===" << std::endl;
    std::cout << "  Aligned reads: " << aligned_count << " / " << read_count 
              << " (" << std::fixed << std::setprecision(1) 
              << (100.0 * aligned_count / read_count) << "%)" << std::endl;
    std::cout << "  Average score: " << std::setprecision(2) << (double)total_score / read_count << std::endl;
    std::cout << "  Max score: " << max_score << std::endl;
    
    // Display sample results
    std::cout << "\nFirst 10 alignments:" << std::endl;
    std::cout << std::setw(6) << "Read" 
              << std::setw(8) << "Seeds"
              << std::setw(10) << "Score" 
              << std::setw(12) << "RefPos"
              << std::setw(10) << "QrySpan" << std::endl;
    std::cout << std::string(46, '-') << std::endl;
    
    for (size_t i = 0; i < std::min(10UL, read_count); i++) {
        std::cout << std::setw(6) << i 
                  << std::setw(8) << hit_counts[i]
                  << std::setw(10) << results[i].score 
                  << std::setw(12) << results[i].ref_pos
                  << std::setw(10) << (results[i].query_end - results[i].query_start) << std::endl;
    }
    
    // Cleanup
    CUDA_CHECK(cudaFree(d_query_reads));
    CUDA_CHECK(cudaFree(d_ref_seq));
    CUDA_CHECK(cudaFree(d_kmer_hashes));
    CUDA_CHECK(cudaFree(d_position_starts));
    CUDA_CHECK(cudaFree(d_position_counts));
    CUDA_CHECK(cudaFree(d_positions));
    CUDA_CHECK(cudaFree(d_seed_hits));
    CUDA_CHECK(cudaFree(d_hit_counts));
    CUDA_CHECK(cudaFree(d_results));
    CUDA_CHECK(cudaEventDestroy(start));
    CUDA_CHECK(cudaEventDestroy(stop));
    
    std::cout << "\nSeed-and-extend alignment complete!" << std::endl;
    
    return 0;
}
