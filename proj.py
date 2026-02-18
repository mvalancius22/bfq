import struct


def pack_qualities(qual_string):
    """Packs ASCII Prhed scores into 4-bit nibbles."""
    packed = bytearray()
    # Process 2 bases at a time
    for i in range(0, len(qual_string), 2):
        # Convert ASCII to 0-scale integer, clamped 0-15
        q1 = min(ord(qual_string[i]) - 33, 15)
        q2 = min(ord(qual_string[i+1]) - 33, 15) if i + 1 < len(qual_string) else 0
        # Pack into on byte: (q2 in high 4 bits, q1 in low 4 bits)
        combined = (q2 << 4) | q1
        packed.append(combined)
    return packed


def pack_sequence(seq):
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    packed = []
    for i in range(0, len(seq), 16):
        chunk = seq[i:i+16]
        word = 0
        for j, base in enumerate(chunk):
            val = mapping.get(base, 0)  # Default to 'A' if unknown base
            word |= (val << (j * 2))
        packed.append(word)
    return packed

def write_bfq(filename, reads):
    read_count = len(reads)
    read_len = len(reads[0])
    with open(filename, 'wb') as f:
        f.write(struct.pack('<4sIQQ', b'BFQ1', 1, read_count, read_len))
        f.write(struct.pack('I', read_len))
        for read in reads:
            packed_read = pack_sequence(read)
            f.write(struct.pack(f'<{len(packed_read)}I', *packed_read))
        
    
if __name__ == "__main__":
    reads = [
        "ACGTACGTACGTACGT",
        "TGCATGCATGCATGCA",
        "GGGGGGGGGGGGGGGG",
        "CCCCCCCCCCCCCCCC"
    ]
    write_bfq('reads.bfq', reads)
    print("BFQ file 'reads.bfq' has been created with the provided reads.") 
    