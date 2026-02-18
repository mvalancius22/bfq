#ifndef BFQ_READER_HPP
#define BFQ_READER_HPP

#include <cstdint>
#include <string>
#include <vector>
#include <fstream>
#include <stdexcept>

struct BFQHeader {
    char magic[4];      // "BFQ1"
    uint32_t version;   // Format version
    uint64_t read_count;
    uint64_t read_len;
    uint32_t read_len_dup;  // Duplicate field from Python code
};

class BFQReader {
public:
    BFQHeader header;
    std::vector<uint32_t> packed_data;  // All packed reads
    
    bool load(const std::string& filename) {
        std::ifstream file(filename, std::ios::binary);
        if (!file) {
            return false;
        }
        
        // Read header
        file.read(reinterpret_cast<char*>(&header), sizeof(BFQHeader));
        
        if (std::string(header.magic, 4) != "BFQ1") {
            throw std::runtime_error("Invalid BFQ magic number");
        }
        
        // Calculate packed words per read (16 bases per uint32)
        size_t words_per_read = (header.read_len + 15) / 16;
        size_t total_words = header.read_count * words_per_read;
        
        // Read all packed data
        packed_data.resize(total_words);
        file.read(reinterpret_cast<char*>(packed_data.data()), 
                  total_words * sizeof(uint32_t));
        
        return file.good();
    }
    
    size_t get_read_count() const { return header.read_count; }
    size_t get_read_length() const { return header.read_len; }
    size_t get_words_per_read() const { return (header.read_len + 15) / 16; }
    const uint32_t* get_data() const { return packed_data.data(); }
};

#endif // BFQ_READER_HPP
