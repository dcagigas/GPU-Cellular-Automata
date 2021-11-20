CC=g++
NVCC=nvcc
CXXFLAGS= -O3 -Wextra -std=c++11 -lm 
CUDAFLAGS= --gpu-architecture=compute_75 -std=c++11 -c   
LIBS= -lcudart

EXE=main

all: $(EXE)

ForestFire_AN5D_512_host.o : ForestFire_AN5D_512_host.cu
	$(NVCC) $(CUDAFLAGS) ForestFire_AN5D_512_host.cu 

ForestFire_AN5D_512_kernel.o: ForestFire_AN5D_512_kernel.cu 
	$(NVCC) $(CUDAFLAGS) ForestFire_AN5D_512_kernel.cu  

main: ForestFire_AN5D_512_kernel.o ForestFire_AN5D_512_host.o
	$(CC) -o ForestFire_an5d ForestFire_AN5D_512_kernel.o ForestFire_AN5D_512_host.o $(CXXFLAGS) $(LIBS)

clean:
	rm -rf ForestFire_an5d *.o


# TO TAKE INTO ACCOUNT:
# 1) The real architecture must be an implementation of the virtual architecture
# 2) The virtual architecture should always be chosen as low as possible, thereby maximizing the actual GPUs to run on. 
# 3) The real architecture should be chosen as high as possible 
#
# GPU REAL Architecture Feature List:
# https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html
# 
# 	                        Basic features
# sm_30 and sm_32           + Kepler support
#                           + Unified memory programming
# sm_35	                    + Dynamic parallelism support
# sm_50, sm_52, and sm_53	+ Maxwell support
# sm_60, sm_61, and sm_62	+ Pascal support
# sm_70 and sm_72	        + Volta support
# sm_75	                    + Turing support
#
#
# GPU VIRTUAL Architecture Feature List:
# https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html
#
# 	                                        Basic features
# compute_30 and compute_32                 + Kepler support
#                                           + Unified memory programming
# compute_35	                            + Dynamic parallelism support
# compute_50, compute_52, and compute_53	+ Maxwell support
# compute_60, compute_61, and compute_62	+ Pascal support
# compute_70 and compute_72	                + Volta support
# compute_75	                            + Turing support

