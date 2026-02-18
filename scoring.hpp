#ifndef SCORING_HPP
#define SCORING_HPP

#include <cstdint>

/**
 * Scoring parameters for sequence alignment
 */
struct AlignmentScoring {
    int8_t match_score;        // Score for matching bases (e.g., +2)
    int8_t mismatch_penalty;   // Penalty for mismatches (e.g., -1)
    int8_t gap_open_penalty;   // Penalty for opening a gap (e.g., -3)
    int8_t gap_extend_penalty; // Penalty for extending a gap (e.g., -1)
    
    // Default constructor with common values
    AlignmentScoring() 
        : match_score(2)
        , mismatch_penalty(-1)
        , gap_open_penalty(-3)
        , gap_extend_penalty(-1)
    {}
    
    AlignmentScoring(int8_t match, int8_t mismatch, int8_t gap_open, int8_t gap_extend)
        : match_score(match)
        , mismatch_penalty(mismatch)
        , gap_open_penalty(gap_open)
        , gap_extend_penalty(gap_extend)
    {}
};

/**
 * Precomputed scoring matrix for 2-bit encoded bases
 * Matrix[query_base][ref_base] = score
 * 0=A, 1=C, 2=G, 3=T
 */
struct ScoringMatrix {
    int8_t matrix[4][4];
    
    void initialize(const AlignmentScoring& scoring) {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                matrix[i][j] = (i == j) ? scoring.match_score : scoring.mismatch_penalty;
            }
        }
    }
    
    __host__ __device__ inline int8_t score(uint8_t base_a, uint8_t base_b) const {
        return matrix[base_a & 0x3][base_b & 0x3];
    }
};

#endif // SCORING_HPP
