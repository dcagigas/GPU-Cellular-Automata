all: clean GOL_cu GOL_shared_cu GOL_texture_cu

GOL_cu: GOL.cu
	nvcc -arch=sm_70 -o GOL_cu.bin GOL.cu

GOL_shared_cu: GOL_shared.cu
	nvcc -arch=sm_70 -o GOL_shared_cu.bin GOL_shared.cu

GOL_texture_cu: GOL_texture.cu
	nvcc -arch=sm_70 -o GOL_texture_cu.bin GOL_texture.cu


clean:
	rm -rf *.bin
