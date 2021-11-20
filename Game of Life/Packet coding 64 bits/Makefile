
all: clean GOL_custom_int64 


GOL_custom_int64: GOL_custom_int64.cu
	nvcc -arch=sm_75 -o GOL_custom_int64 GOL_custom_int64.cu

clean:
	rm -rf *.o GOL_custom_int64 
