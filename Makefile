NVCC = nvcc
CXX = g++
CUDA_ARCH = -arch=sm_75  # Adjust for your GPU (sm_75 = Turing, sm_80 = Ampere, sm_89 = Ada)

NVCC_FLAGS = -O3 $(CUDA_ARCH) -std=c++11 --use_fast_math
CXX_FLAGS = -O3 -std=c++11

TARGET = bfq_unpack
OBJECTS = main.o unpack_kernel.o

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(NVCC) $(NVCC_FLAGS) -o $@ $^

main.o: main.cpp bfq_reader.hpp
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@

unpack_kernel.o: unpack_kernel.cu
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJECTS) *.bfq

test: $(TARGET)
	python3 proj.py
	./$(TARGET) reads.bfq

.PHONY: all clean test
