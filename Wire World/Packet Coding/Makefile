
all: clean WireWorld_cuda_pc8 


WireWorld_cuda_pc8: WireWorld_cuda_pc8.cu
	nvcc -arch=sm_75 -o WireWorld_cuda_pc8 WireWorld_cuda_pc8.cu

clean:
	rm -rf *.o WireWorld_cuda_pc8 
