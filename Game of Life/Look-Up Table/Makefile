
all: clean GOL_cu 


GOL_cu: GOL.cu
	nvcc --gpu-architecture=compute_75 -std=c++11 -o GOL_cu GOL.cu

clean:
	rm -rf *.o GOL_cu 
