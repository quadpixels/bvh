// 2017-10-08 ok, to run visualization, type 
//   `USE_CUDA=1 LD_LIBRARY_PATH=/usr/lib/nvidia-375 optirun ./bvh_vis`

#include "bvh.hpp"
#include "utils.hpp"
#include <stdlib.h>
#include <stdio.h>
#include <unordered_set>
#include <unordered_map>

#define CE(call) {\
	call; \
	cudaError_t err = cudaGetLastError(); \
	if(err != cudaSuccess) { \
		printf("%s\n", cudaGetErrorString(err)); \
		assert(0); \
	} \
}

int g_num_blk = 64, g_num_thd = 64;

__global__ void BVHReduce_Kernel(int* d_left, int* d_right, int* d_parent, int* d_node_lb, int* d_node_ub,
	Vec3* d_min_ext, Vec3* d_max_ext, Vec3* d_points_sorted, int* d_mc, int* d_reduce_counter, int* d_lock,
	const int N);

__global__ void BVHGen_Kernel(int* d_left, int* d_right, int* d_parent, int* d_node_lb, int* d_node_ub,
	Vec3* d_min_ext, Vec3* d_max_ext, Vec3* d_points_sorted, int* d_mc, int* d_reduce_counter, int* d_lock, const int N) {
	const int nthd = blockDim.x * gridDim.x,
	          tidx = threadIdx.x + blockIdx.x * blockDim.x;
	for (int i=tidx; i<=N-2; i += nthd) {
		int delta_left  = (i-1<0)   ? -1 : __clz(d_mc[i-1] ^ d_mc[i]),
		    delta_right = (i+1>N-1) ? -1 : __clz(d_mc[i] ^ d_mc[i+1]),
			d = (delta_right > delta_left) ? 1 : -1,
			delta_min = -999;
		if (i-d >= 0 && i-d <= N-1) {
			delta_min = __clz(d_mc[i] ^ d_mc[i-d]);
		} else delta_min = -1;

		// Find the other end
		int j=i, dj=1*d;

		// Binary Search
		while (j+dj >= 0 && j+dj <= N-1) {
			const int d1 = __clz(d_mc[i] ^ d_mc[j+dj]);
			if (d1 >= delta_min) dj *= 2; else break;
		}
		int dj_tot = 0;
		while (true) {
			dj /= 2;
			if (j+dj_tot+dj >= 0 && j+dj_tot+dj <= N-1) {
				const int d1 = __clz(d_mc[i] ^ d_mc[j+dj_tot+dj]);
				if (d1 >= delta_min) {
					dj_tot += dj;
				}
			}
			if (dj <= 1 && dj >= -1) break;
		}
		j = j+dj_tot;

		int lb = min(i, j), ub = max(i, j), 
		    delta_min_ij = __clz(d_mc[lb]^d_mc[ub]);

		// Find Split Point
		int lb1 = lb, ub1 = ub, split;
		while (lb1+1 < ub1) {
			int mid = (lb1 + ub1) / 2;
			if (__clz(d_mc[mid] ^ d_mc[lb]) > delta_min_ij) lb1 = mid;
			else ub1 = mid;
		}
		split = ub1;
		while (split-1 >= lb &&
			__clz(d_mc[split] ^ d_mc[lb]) <= delta_min_ij) {
			split --;
		}

		int left = split, right = split + 1;
		d_node_lb[i] = lb;
		d_node_ub[i] = ub;
		if (lb == ub) {
			*((int*)0x00000000) = 0xBAADCAFE;
		} else {
			if (left == lb) {
				d_left  [i]          = left + N-1;
				d_parent[left + N-1] = i;
			} else {
				d_left  [i]    = left;
				d_parent[left] = i;
			}

			if (right == ub) {
				d_right [i]           = right + N-1;
				d_parent[right + N-1] = i;
			} else {
				d_right [i]     = right;
				d_parent[right] = i;
			}
		}

		if (i==0) {
			printf("i=%d left=%d right=%d\n", i, d_left[i], d_right[i]);
		}
	}
}

BVH BVHGen_GPU(const std::vector<Vec3>& points,
               const std::vector<Vec3>& points_transformed) {
	// 1. compute morton code and Sort
	std::vector<unsigned int> mc;
	std::unordered_set<int> occ;
	std::unordered_map<unsigned, std::vector<int> > mc2point;
	for (int i=0; i<int(points_transformed.size()); i++) {
		Vec3 p = points_transformed[i];
		unsigned int x = std::max(std::min(1.0f, p.x), 0.0f) * 1024,
		             y = std::max(std::min(1.0f, p.y), 0.0f) * 1024,
		             z = std::max(std::min(1.0f, p.z), 0.0f) * 1024;
		int code = ExpandBits(x) * 4 + ExpandBits(y) * 2 + ExpandBits(z);
		if (occ.find(code) != occ.end()) {
		} else { mc.push_back(code); }
		mc2point[code].push_back(i);
		occ.insert(code);
	}

	printf("%lu unique Morton codes for %lu points\n",
		mc.size(), points.size());

	std::vector<Vec3> points_sorted;
	std::vector<std::pair<int, int> > mc2sortedidx;  // an MC corresponds to which
	                                                 //   sorted indices
	int* mapping;
	radixsort<unsigned int>(mc.data(), int(mc.size()), &mapping);
	for (int i=0; i<int(mc.size()); i++) {
		unsigned code = mc[i];
		int idx_min = int(points_sorted.size());
		for (const int orig_idx : mc2point.at(code)) {
			points_sorted.push_back(points.at(orig_idx));
		}
		int idx_max = int(points_sorted.size())-1;
		assert (idx_max >= idx_min);
		mc2sortedidx.push_back(std::make_pair(idx_min, idx_max));
	}
	assert (points_sorted.size() == points.size());

	// 2. transfer to gpu & alloc space

	int *d_mc, *d_left, *d_right, *d_parent,
	    *d_node_lb, *d_node_ub, *d_reduce_counter,
		*d_lock;
	Vec3* d_min_ext, *d_max_ext, *d_points_sorted;
	const int N = int(mc.size());

	CE(cudaMalloc(&d_mc, sizeof(int)*N));
	CE(cudaMemcpy(d_mc, mc.data(), sizeof(int)*N, cudaMemcpyHostToDevice));
	
	CE(cudaMalloc(&d_left,           sizeof(int) *(2*N-1)));
	CE(cudaMalloc(&d_right,          sizeof(int) *(2*N-1)));

	CE(cudaMalloc(&d_parent,         sizeof(int) *(2*N-1)));
	int *tmp = new int[2*N-1];
	for (int i=0; i<2*N-1; i++) tmp[i] = -999;
	CE(cudaMemcpy(d_parent, tmp, sizeof(int)*(2*N-1), cudaMemcpyHostToDevice));
	CE(cudaMemcpy(d_left  , tmp, sizeof(int)*(2*N-1), cudaMemcpyHostToDevice));
	CE(cudaMemcpy(d_right , tmp, sizeof(int)*(2*N-1), cudaMemcpyHostToDevice));

	CE(cudaMalloc(&d_node_lb,        sizeof(int) *(2*N-1)));
	CE(cudaMalloc(&d_node_ub,        sizeof(int) *(2*N-1)));

	CE(cudaMalloc(&d_min_ext,        sizeof(Vec3)*(2*N-1)));
	CE(cudaMalloc(&d_max_ext,        sizeof(Vec3)*(2*N-1)));
	Vec3* h_exts = new Vec3[2*N-1];
	for (int i=0; i<2*N-1; i++) h_exts[i] = Vec3( 1e20, 1e20, 1e20);
	CE(cudaMemcpy(d_min_ext, h_exts, sizeof(Vec3)*(2*N-1), cudaMemcpyHostToDevice));
	for (int i=0; i<2*N-1; i++) h_exts[i] = Vec3(-1e20,-1e20,-1e20);
	CE(cudaMemcpy(d_max_ext, h_exts, sizeof(Vec3)*(2*N-1), cudaMemcpyHostToDevice));
	delete[] h_exts;

	CE(cudaMalloc(&d_lock,           sizeof(int)*(2*N-1)));
	CE(cudaMemset(d_lock, 0x00,      sizeof(int)*(2*N-1)));

	CE(cudaMalloc(&d_points_sorted,  sizeof(Vec3)*points.size()));
	CE(cudaMemcpy(d_points_sorted, points_sorted.data(),
		sizeof(Vec3)*points_sorted.size(),
		cudaMemcpyHostToDevice));

	CE(cudaMalloc(&d_reduce_counter, sizeof(int)*(2*N-1)));
	CE(cudaMemset(d_reduce_counter, 0x00, sizeof(int)*(2*N-1)));

	cudaEvent_t t0, t1;
	CE(cudaEventCreate(&t0));
	CE(cudaEventCreate(&t1));
	CE(cudaEventRecord(t0));
	BVHGen_Kernel<<<120, 192>>>(d_left, d_right, d_parent, d_node_lb, d_node_ub,
		d_min_ext, d_max_ext, d_points_sorted, d_mc, d_reduce_counter, d_lock, N);
	CE(cudaEventRecord(t1));
	CE(cudaEventSynchronize(t0));
	CE(cudaEventSynchronize(t1));
	float ms;
	CE(cudaEventElapsedTime(&ms, t0, t1));
	printf("GPU Hierarchy Generation Kernel Time: %fms\n", ms);

	int *h_left  = new int[2*N-1],
	    *h_right = new int[2*N-1],
		*h_parent = new int[2*N-1],
		*h_node_lb = new int[2*N-1],
		*h_node_ub = new int[2*N-1];
	CE(cudaMemcpy(h_parent,  d_parent,  sizeof(int)*(2*N-1), cudaMemcpyDeviceToHost));

	// The lookup may need 2 be performed on CPU
	// Write the leaf nodes then xfer to GPU
	CE(cudaMemcpy(h_node_lb, d_node_lb, sizeof(int)*(2*N-1), cudaMemcpyDeviceToHost));
	CE(cudaMemcpy(h_node_ub, d_node_ub, sizeof(int)*(2*N-1), cudaMemcpyDeviceToHost));

	for (int i=0; i<N-1; i++) {
		h_node_lb[i] = mc2sortedidx.at(h_node_lb[i]).first;
		h_node_ub[i] = mc2sortedidx.at(h_node_ub[i]).second;
	}
	for (int i=N-1; i<2*N-1; i++) {
		h_node_lb[i] = mc2sortedidx.at(i-(N-1)).first;
		h_node_ub[i] = mc2sortedidx.at(i-(N-1)).second;
	}
	CE(cudaMemcpy(d_node_lb, h_node_lb, sizeof(int)*(2*N-1), cudaMemcpyHostToDevice));
	CE(cudaMemcpy(d_node_ub, h_node_ub, sizeof(int)*(2*N-1), cudaMemcpyHostToDevice));

	CE(cudaEventCreate(&t0));
	CE(cudaEventCreate(&t1));
	CE(cudaEventRecord(t0));
	BVHReduce_Kernel<<<g_num_blk, g_num_thd>>>(d_left, d_right, d_parent, d_node_lb, d_node_ub,
		d_min_ext, d_max_ext, d_points_sorted, d_mc, d_reduce_counter, d_lock, N);
	CE(cudaEventRecord(t1));
	CE(cudaEventSynchronize(t0));
	CE(cudaEventSynchronize(t1));
	CE(cudaEventElapsedTime(&ms, t0, t1));
	printf("GPU Parallel Reduction Kernel Time: %fms\n", ms);

	Vec3 *h_min_ext = new Vec3[2*N-1],
	     *h_max_ext = new Vec3[2*N-1];

	CE(cudaMemcpy(h_left,    d_left,    sizeof(int)*(2*N-1), cudaMemcpyDeviceToHost));
	CE(cudaMemcpy(h_right,   d_right,   sizeof(int)*(2*N-1), cudaMemcpyDeviceToHost));
	CE(cudaMemcpy(h_parent,  d_parent,  sizeof(int)*(2*N-1), cudaMemcpyDeviceToHost));
	CE(cudaMemcpy(h_node_lb, d_node_lb, sizeof(int)*(2*N-1), cudaMemcpyDeviceToHost));
	CE(cudaMemcpy(h_node_ub, d_node_ub, sizeof(int)*(2*N-1), cudaMemcpyDeviceToHost));
	CE(cudaMemcpy(h_min_ext, d_min_ext, sizeof(Vec3)*(2*N-1), cudaMemcpyDeviceToHost));
	CE(cudaMemcpy(h_max_ext, d_max_ext, sizeof(Vec3)*(2*N-1), cudaMemcpyDeviceToHost));
	BVH h_bvh;

	h_bvh.nodes.resize(2*N-1);
	for (int i=0; i<2*N-1; i++) {
		BVH::BVHNode* n = &(h_bvh.nodes[i]);
		n->left    = h_left[i];
		n->right   = h_right[i];
		n->parent  = h_parent[i];
		n->node_lb = h_node_lb[i];
		n->node_ub = h_node_ub[i];
		n->min_ext = h_min_ext[i];
		n->max_ext = h_max_ext[i];
	}

	printf("mapping=%p\n", mapping);
	h_bvh.mapping = mapping;
	h_bvh.PrintSummary();
	h_bvh.mc = mc;
	h_bvh.mc2point = mc2point;
	return h_bvh;
}
