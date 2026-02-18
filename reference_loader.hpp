#ifndef REFERENCE_LOADER_HPP
#define REFERENCE_LOADER_HPP

#include <cstdint>
#include <string>
#include <vector>
#include <fstream>
#include <cctype>

/**
 * Simple reference sequence loader that packs sequences into 2-bit format
 * Supports loading from FASTA format
 */
class ReferenceLoader {
public:
    std::vector<uint32_t> packed_data;
    uint64_t sequence_length;
    
    ReferenceLoader() : sequence_length(0) {}
    
    /**
     * Load reference from FASTA file
     */
    bool load_fasta(const std::string& filename) {
        std::ifstream file(filename);
        if (!file) {
            return false;
        }
        
        std::string sequence;
        std::string line;
        
        while (std::getline(file, line)) {
            if (line.empty()) continue;
            
            // Skip header lines
            if (line[0] == '>') {
                continue;
            }
            
            // Accumulate sequence
            for (char c : line) {
                if (std::isspace(c)) continue;
                sequence += std::toupper(c);
            }
        }
        
        if (sequence.empty()) {
            return false;
        }
        
        sequence_length = sequence.length();
        pack_sequence(sequence);
        return true;
    }
    
    /**
     * Load reference from raw sequence string
     */
    bool load_string(const std::string& sequence) {
        if (sequence.empty()) {
            return false;
        }
        
        sequence_length = sequence.length();
        pack_sequence(sequence);
        return true;
    }
    
    /**
     * Create a synthetic reference (for testing)
     */
    void create_synthetic(uint64_t length, uint32_t seed = 42) {
        std::string sequence;
        sequence.reserve(length);
        
        const char bases[] = {'A', 'C', 'G', 'T'};
        uint32_t rng = seed;
        
        for (uint64_t i = 0; i < length; i++) {
            // Simple LCG random number generator
            rng = rng * 1103515245 + 12345;
            sequence += bases[(rng >> 16) & 0x3];
        }
        
        sequence_length = length;
        pack_sequence(sequence);
    }
    
    size_t get_length() const { 
        return sequence_length; 
    }
    
    size_t get_words_per_sequence() const { 
        return (sequence_length + 15) / 16; 
    }
    
    const uint32_t* get_data() const { 
        return packed_data.data(); 
    }
    
private:
    /**
     * Pack ASCII sequence into 2-bit format
     * Same encoding as BFQ: A=00, C=01, G=10, T=11
     */
    void pack_sequence(const std::string& sequence) {
        size_t words_needed = (sequence.length() + 15) / 16;
        packed_data.resize(words_needed);
        
        for (size_t word_idx = 0; word_idx < words_needed; word_idx++) {
            uint32_t word = 0;
            
            for (int bit_idx = 0; bit_idx < 16; bit_idx++) {
                size_t base_idx = word_idx * 16 + bit_idx;
                
                if (base_idx >= sequence.length()) {
                    break;
                }
                
                uint32_t base_val = encode_base(sequence[base_idx]);
                word |= (base_val << (bit_idx * 2));
            }
            
            packed_data[word_idx] = word;
        }
    }
    
    /**
     * Encode single base to 2-bit value
     */
    uint32_t encode_base(char base) const {
        switch (base) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': return 3;
            default: return 0;  // Unknown bases default to 'A'
        }
    }
};

#endif // REFERENCE_LOADER_HPP
