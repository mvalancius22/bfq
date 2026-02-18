#ifndef ALIGNMENT_RESULT_HPP
#define ALIGNMENT_RESULT_HPP

#include <cstdint>
#include <string>
#include <vector>

/**
 * Represents alignment result for a single read
 */
struct AlignmentResult {
    uint64_t read_id;           // Index of the read
    int32_t ref_pos;            // Starting position in reference (0-based)
    int32_t alignment_score;    // Smith-Waterman alignment score
    uint32_t query_start;       // Start position in query (0-based)
    uint32_t query_end;         // End position in query (exclusive)
    uint32_t ref_start;         // Start position in reference (0-based)
    uint32_t ref_end;           // End position in reference (exclusive)
    bool aligned;               // Whether alignment succeeded
    
    AlignmentResult() 
        : read_id(0)
        , ref_pos(-1)
        , alignment_score(0)
        , query_start(0)
        , query_end(0)
        , ref_start(0)
        , ref_end(0)
        , aligned(false)
    {}
};

/**
 * Compact GPU-friendly alignment result (no CIGAR for now)
 */
struct CompactAlignment {
    int32_t score;        // Alignment score
    int32_t ref_pos;      // Reference start position
    uint16_t query_start; // Query start position
    uint16_t query_end;   // Query end position
    
    __host__ __device__ CompactAlignment() 
        : score(0), ref_pos(-1), query_start(0), query_end(0) {}
};

/**
 * Batch alignment results container
 */
struct BatchAlignmentResults {
    std::vector<AlignmentResult> results;
    size_t aligned_count;
    size_t total_count;
    
    BatchAlignmentResults(size_t count) 
        : results(count)
        , aligned_count(0)
        , total_count(count)
    {}
    
    void add_result(const AlignmentResult& result) {
        if (result.aligned) {
            aligned_count++;
        }
    }
    
    double get_alignment_rate() const {
        return total_count > 0 ? (double)aligned_count / total_count : 0.0;
    }
};

#endif // ALIGNMENT_RESULT_HPP
