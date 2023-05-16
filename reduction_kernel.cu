// Turning on -O2 seems to cause deadlock
// Must use -G

#include "bvh.hpp"
#include "utils.hpp"

// To prevent deadlock some instruction must be added after
//   the common end of the 2 if-else branches
__device__ __noinline__ void do_ProcessOneNode(int* d_lock, Vec3* d_min_ext, Vec3* d_max_ext, int par, int idx) {
	bool done = false;
	while (!done) {
		int* the_lock = &d_lock[par];
		if (atomicCAS(the_lock, 0, 1) == 0) {
			d_min_ext[par].x = min(d_min_ext[par].x, d_min_ext[idx].x);
			d_min_ext[par].y = min(d_min_ext[par].y, d_min_ext[idx].y);
			d_min_ext[par].z = min(d_min_ext[par].z, d_min_ext[idx].z);
			d_max_ext[par].x = max(d_max_ext[par].x, d_max_ext[idx].x);
			d_max_ext[par].y = max(d_max_ext[par].y, d_max_ext[idx].y);
			d_max_ext[par].z = max(d_max_ext[par].z, d_max_ext[idx].z);
			done = true;
			*the_lock = 0;
		} else {
		}
//		if (threadIdx.x == -999) printf("Ho!\n"); // works on NVS 5200M
		volatile bool enod = done;                // works on NVS 5200M
//		__asm("pmevent 0;");                      // works on NVS 5200M
//		volatile int dummy = *the_lock;           // works on NVS 5200M
//		int x = *the_lock;                        // does NOT work on NVS 5200M
	}
}

__global__ void BVHReduce_Kernel(int* d_left, int* d_right, int* d_parent, int* d_node_lb, int* d_node_ub,
	Vec3* d_min_ext, Vec3* d_max_ext, Vec3* d_points_sorted, int* d_mc, int* d_reduce_counter, int* d_lock,
	const int N) {
	const int nthd = blockDim.x * gridDim.x,
	          tidx = threadIdx.x + blockIdx.x * blockDim.x;

	for (int leaf_id = tidx, iter = 0; leaf_id < N; leaf_id += nthd) {
		int idx = N-1+leaf_id;
		iter = 0;
		while (idx != -999) {
			if (iter == 0) {
				int primitive_id = d_node_lb[idx];
				d_min_ext[idx] = d_max_ext[idx] = d_points_sorted[primitive_id];
			}
			int par = d_parent[idx];
			bool skip_to_next_leaf = false;
			bool done = false;
			if (par == -999) {
				skip_to_next_leaf = true;
				idx = -999;
				done = true;
			}

			int c = atomicAdd(&d_reduce_counter[par], 1);
			skip_to_next_leaf = (c == 0);

			if (!done) do_ProcessOneNode(d_lock, d_min_ext, d_max_ext, par, idx);

			if (skip_to_next_leaf) {
				idx = -999;
			} else {
				idx = par;
				iter ++;
			}
		}
	}
}
