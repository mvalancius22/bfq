NVCC = nvcc
CXX = g++
CUDA_ARCH = -arch=sm_75  # Adjust for your GPU (sm_75 = Turing, sm_80 = Ampere, sm_89 = Ada)

NVCC_FLAGS = -O3 $(CUDA_ARCH) -std=c++11 --use_fast_math
CXX_FLAGS = -O3 -std=c++11

# Targets
TARGET_UNPACK = bfq_unpack
TARGET_ALIGN = bfq_align
TARGET_SEED = bfq_seed_align

OBJECTS_UNPACK = main.o unpack_kernel.o
OBJECTS_ALIGN = align_main.o alignment_kernel.o
OBJECTS_SEED = seed_align_main.o seed_extend_kernel.o

all: $(TARGET_UNPACK) $(TARGET_ALIGN) $(TARGET_SEED)

# Unpack target
$(TARGET_UNPACK): $(OBJECTS_UNPACK)
	$(NVCC) $(NVCC_FLAGS) -o $@ $^

main.o: main.cpp bfq_reader.hpp
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@

unpack_kernel.o: unpack_kernel.cu
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@

# Alignment target
$(TARGET_ALIGN): $(OBJECTS_ALIGN)
	$(NVCC) $(NVCC_FLAGS) -o $@ $^

align_main.o: align_main.cpp bfq_reader.hpp reference_loader.hpp scoring.hpp alignment_result.hpp
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@

alignment_kernel.o: alignment_kernel.cu scoring.hpp alignment_result.hpp
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@

# Seed-and-extend target
$(TARGET_SEED): $(OBJECTS_SEED)
	$(NVCC) $(NVCC_FLAGS) -o $@ $^

seed_align_main.o: seed_align_main.cpp bfq_reader.hpp reference_loader.hpp scoring.hpp alignment_result.hpp seed_index.hpp
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@

seed_extend_kernel.o: seed_extend_kernel.cu scoring.hpp alignment_result.hpp seed_index.hpp
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@

clean:
	rm -f $(TARGET_UNPACK) $(TARGET_ALIGN) $(TARGET_SEED) $(OBJECTS_UNPACK) $(OBJECTS_ALIGN) $(OBJECTS_SEED) *.bfq

test: $(TARGET_UNPACK)
	python3 proj.py
	./$(TARGET_UNPACK) reads.bfq

test_align: $(TARGET_ALIGN)
	python3 proj.py
	./$(TARGET_ALIGN) reads.bfq -s ACGTACGTACGTACGT

test_seed: $(TARGET_SEED)
	python3 proj.py
	./$(TARGET_SEED) reads.bfq -l 1000 -k 15

.PHONY: all clean test test_align test_seed
