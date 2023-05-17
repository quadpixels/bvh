all: bvh bvh_vis

NVCC=/usr/local/cuda-11.2/bin/nvcc
NVCCFLAGS=-arch=sm_50 -O2
CXXFLAGS=-O2 -DUSE_CUDA

bvh: main.cpp utils.hpp bvh.cpp utils.cpp cudastuff.cu reduction_kernel.cu
	$(NVCC) reduction_kernel.cu -std=c++11 -arch=sm_50 -O3 -c -o reduction_kernel.o
	$(NVCC) cudastuff.cu -std=c++11 $(NVCCFLAGS) -c -o cudastuff.o
	$(NVCC) -ccbin g++ -std=c++11 $(CXXFLAGS) main.cpp bvh.cpp utils.cpp cudastuff.o reduction_kernel.o -o $@ -lcuda

bvh_vis: main.cpp utils.hpp openglstuff.hpp bvh.cpp utils.cpp cudastuff.cu
	$(NVCC) cudastuff.cu -std=c++11 $(NVCCFLAGS) -c -o cudastuff.o
	$(NVCC) -std=c++11 $(CXXFLAGS) -D VISUALIZE main.cpp bvh.cpp utils.cpp cudastuff.o reduction_kernel.o -o $@ -lGL -lGLEW -lGLU -lglut

clean:
	@for x in *.o bvh bvh_vis bvh_cpuonly bvh_vis_cpuonly; do \
	  if [ -f $$x ]; then echo "Removing $$x"; rm -v $$x ; fi ; \
	  done
