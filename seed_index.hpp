#ifndef SEED_INDEX_HPP
#define SEED_INDEX_HPP

#include <cstdint>
#include <vector>
#include <unordered_map>
#include <string>

/**
 * Seed position in reference
 */
struct SeedPosition {
    uint32_t position;    // Position in reference
    uint32_t kmer_hash;   // K-mer hash value
    
    SeedPosition() : position(0), kmer_hash(0) {}
    SeedPosition(uint32_t pos, uint32_t hash) : position(pos), kmer_hash(hash) {}
};

/**
 * Seed hit (match between query and reference)
 */
struct SeedHit {
    uint32_t query_pos;    // Position in query read
    uint32_t ref_pos;      // Position in reference
    uint32_t kmer_hash;    // K-mer hash (for validation)
    int32_t diagonal;      // Diagonal = ref_pos - query_pos (for chaining)
    
    __host__ __device__ SeedHit() 
        : query_pos(0), ref_pos(0), kmer_hash(0), diagonal(0) {}
    
    __host__ __device__ SeedHit(uint32_t qpos, uint32_t rpos, uint32_t hash)
        : query_pos(qpos), ref_pos(rpos), kmer_hash(hash)
        , diagonal(rpos - qpos) {}
};

/**
 * K-mer index for reference sequence
 * Uses hash table for fast seed lookup
 */
class SeedIndex {
public:
    static constexpr int DEFAULT_K = 15;  // 15-mer seeds (common in aligners)
    static constexpr int MAX_POSITIONS_PER_KMER = 1000;  // Ignore repetitive k-mers
    
    int k;  // K-mer length
    uint64_t ref_length;
    
    // Hash table: kmer_hash -> list of positions
    std::unordered_map<uint32_t, std::vector<uint32_t>> kmer_positions;
    
    // Flattened array for GPU (kmer_hash -> start_idx in positions array)
    std::vector<uint32_t> kmer_hashes;      // Unique k-mer hashes
    std::vector<uint32_t> position_starts;  // Start index in positions array
    std::vector<uint32_t> position_counts;  // Number of positions for this k-mer
    std::vector<uint32_t> positions;        // Flattened positions array
    
    SeedIndex(int kmer_size = DEFAULT_K) : k(kmer_size), ref_length(0) {}
    
    /**
     * Build index from packed reference sequence
     */
    void build(const uint32_t* packed_ref, uint64_t length, size_t words_per_ref) {
        ref_length = length;
        kmer_positions.clear();
        
        if (length < k) return;
        
        // Extract all k-mers from reference
        for (uint64_t pos = 0; pos <= length - k; pos++) {
            uint32_t kmer_hash = extract_kmer_hash(packed_ref, pos, words_per_ref);
            
            // Skip if this k-mer is too repetitive
            if (kmer_positions[kmer_hash].size() < MAX_POSITIONS_PER_KMER) {
                kmer_positions[kmer_hash].push_back(pos);
            }
        }
        
        // Flatten for GPU
        flatten_index();
    }
    
    /**
     * Get number of unique k-mers in index
     */
    size_t get_kmer_count() const {
        return kmer_hashes.size();
    }
    
    /**
     * Get total number of seed positions
     */
    size_t get_total_positions() const {
        return positions.size();
    }
    
    /**
     * Get GPU-friendly data
     */
    const uint32_t* get_kmer_hashes() const { return kmer_hashes.data(); }
    const uint32_t* get_position_starts() const { return position_starts.data(); }
    const uint32_t* get_position_counts() const { return position_counts.data(); }
    const uint32_t* get_positions() const { return positions.data(); }
    
private:
    /**
     * Extract k-mer starting at position and compute hash
     */
    uint32_t extract_kmer_hash(const uint32_t* packed, uint64_t start_pos, size_t words_per_ref) const {
        uint32_t hash = 0;
        
        for (int i = 0; i < k; i++) {
            uint8_t base = extract_base(packed, start_pos + i, words_per_ref);
            // Simple rolling hash: mix bases into 32-bit value
            hash = (hash << 2) | base;
        }
        
        return hash;
    }
    
    /**
     * Extract single base from packed data
     */
    uint8_t extract_base(const uint32_t* packed, uint64_t pos, size_t words_per_ref) const {
        size_t word_idx = pos / 16;
        size_t bit_offset = (pos % 16) * 2;
        return (packed[word_idx] >> bit_offset) & 0x3;
    }
    
    /**
     * Flatten hash table into arrays for GPU
     */
    void flatten_index() {
        kmer_hashes.clear();
        position_starts.clear();
        position_counts.clear();
        positions.clear();
        
        uint32_t current_start = 0;
        
        for (const auto& entry : kmer_positions) {
            uint32_t hash = entry.first;
            const auto& pos_list = entry.second;
            
            // Skip repetitive k-mers
            if (pos_list.size() > MAX_POSITIONS_PER_KMER) continue;
            
            kmer_hashes.push_back(hash);
            position_starts.push_back(current_start);
            position_counts.push_back(pos_list.size());
            
            // Add all positions
            for (uint32_t pos : pos_list) {
                positions.push_back(pos);
            }
            
            current_start += pos_list.size();
        }
    }
};

/**
 * Parameters for seed-and-extend alignment
 */
struct SeedExtendParams {
    int k;                    // K-mer length for seeding
    int max_seed_hits;        // Maximum hits per read to extend
    int band_width;           // Width of band for banded alignment
    int min_seed_score;       // Minimum score to keep a seed extension
    int seed_interval;        // Sample seeds every N bases (stride)
    
    SeedExtendParams()
        : k(15)
        , max_seed_hits(100)
        , band_width(32)
        , min_seed_score(20)
        , seed_interval(1)  // Sample every base by default
    {}
};

#endif // SEED_INDEX_HPP
