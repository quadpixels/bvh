#include <stdio.h>
#include <math.h>
#include <signal.h>
#include <sys/time.h>

#include <vector>
#include <unordered_set>
#include "bvh.hpp"
#include "utils.hpp"

void BVH::Build() {
	ComputeMortonCodes();
	int N = int(mc.size());
	nodes.resize(N*2 - 1);
	internal_node_idx = 0;
	leaf_node_idx = N-1;

	struct timeval t0, t1;
	gettimeofday(&t0, NULL);
	radixsort<unsigned int>(mc.data(), int(mc.size()), &mapping);
	do_build(-999, 0, int(mc.size())-1);
	gettimeofday(&t1, NULL);
	float ms = (t1.tv_sec - t0.tv_sec) * 1000.0f +
	           (t1.tv_usec- t0.tv_usec)/ 1000.0f;
	printf("Build time: %f ms\n", ms);

	gettimeofday(&t0, NULL);

	std::vector<Vec3> points_sorted;
	std::vector<std::pair<int, int> > mc2sortedidx;  // an MC corresponds to which
	                                                 //   sorted indices
	for (int i=0, idx=0; i<int(mc.size()); i++) {
		unsigned code = mc[i];
		int idx_min = int(points_sorted.size());
		for (const int orig_idx : mc2point.at(code)) {
			points_sorted.push_back(points.at(orig_idx));
		}
		int idx_max = int(points_sorted.size())-1;
		assert (idx_max >= idx_min);
		mc2sortedidx.push_back(std::make_pair(idx_min, idx_max));
	}

	std::unordered_set<int> n2v;
	for (int i=0; i<N; i++) { 
		const int x = i+N-1;
		n2v.insert(x); 
		nodes[x].min_ext = nodes[x].max_ext = points_sorted[mc2sortedidx.at(i).first];
	}

	for (int i=0; i<N*2-1; i++) {
		// If we have triangles then we need to expand triangles here
		nodes.at(i).node_lb = mc2sortedidx.at(nodes.at(i).node_lb).first;
		nodes.at(i).node_ub = mc2sortedidx.at(nodes.at(i).node_ub).second;
	}

	while (!n2v.empty()) {
		std::unordered_set<int> n2v_next;
		for (const int x : n2v) {
			int par = nodes[x].parent;
			if (par != -999) {
				nodes[par].min_ext = Vec3Min(nodes[par].min_ext, nodes[x].min_ext);
				nodes[par].max_ext = Vec3Max(nodes[par].max_ext, nodes[x].max_ext);
				Vec3 m0 = nodes[par].min_ext,
				          m1 = nodes[par].max_ext;
				if (fabs(Vec3::dot(m0, m0)) > 1e20) {
					printf("X: (%f,%f,%f)->(%f,%f,%f)\n",
						nodes[x].min_ext.x,
						nodes[x].min_ext.y,
						nodes[x].min_ext.z,
						nodes[x].max_ext.x,
						nodes[x].max_ext.y,
						nodes[x].max_ext.z);
					assert(0);
				}
				if (fabs(Vec3::dot(m1, m1)) > 1e20) {
					printf("X: (%g,%g,%g)->(%g,%g,%g)\n",
						nodes[x].min_ext.x,
						nodes[x].min_ext.y,
						nodes[x].min_ext.z,
						nodes[x].max_ext.x,
						nodes[x].max_ext.y,
						nodes[x].max_ext.z);
					assert(0);
				}
				n2v_next.insert(par);
			}
		}
		n2v = n2v_next;
	}

	gettimeofday(&t1, NULL);
	ms = (t1.tv_sec - t0.tv_sec) * 1000.0f +
	     (t1.tv_usec- t0.tv_usec)/ 1000.0f;
	printf("Remap and propagate: %f ms\n", ms);

	int count = 0;
	do_count_node(-999, 0, 0, &count);
	printf("Node Count = %d\n", count);
}

int BVH::do_build(int parent, int lb, int ub) {
	BVHNode* n; // = &(nodes[node_idx]);
	int node_idx = -999;
	if (ub <= lb) {
		node_idx = AllocLeafNode();
		n = &(nodes[node_idx]);
		n->node_lb = n->node_ub = lb;
		n->parent = parent;
		return node_idx;
	} else {
		node_idx = AllocInternalNode();
		n = &(nodes[node_idx]);
		n->parent = parent;
	}

	assert (lb >= 0 && ub <= int(mc.size())-1);
	int split = FindSplit(mc.data(), lb, ub);
		/*
	int s1 = FindSplit_ref(mc.data(), lb, ub);
	if (s1 != split) {
		printf("%d vs %d\n", s1, split);
		printf("clz(lb,ub)=%d\n", clz(mc[lb]^mc[ub]));
		printf("lb=%d ub=%d [lb]: %X, [ub]: %X, [s1]:%X(clz=%d), [split]:%X(clz=%d)\n",
			lb, ub, 
			mc[lb], mc[ub], mc[s1], clz(mc[s1]^mc[lb]), mc[split], clz(mc[split]^
				mc[lb]));
		printf("\n");
		assert (0);
	}*/
	n->node_lb = lb;
	n->node_ub = ub;
	n->left   = do_build(node_idx, lb,      split);
	n->right  = do_build(node_idx, split+1, ub);
	n->parent = parent;

	if (false && (node_idx == 99997 || node_idx == 99998)) {
		printf("i=%d [LB,UB]=[%d,%d] split=%d left=%d right=%d\n",
			node_idx, lb, ub, split, n->left, n->right);
	}
	return node_idx;
}

void BVH::do_count_node(int parent, int level, int node_idx, int* count) {
	if (node_idx >= int(nodes.size())) return;
	BVHNode* n = &(nodes[node_idx]);
	assert (n->parent == parent);

	if (parent != -999) {
		BVHNode* pn = &(nodes[parent]);
		assert (pn->node_lb <= n->node_lb && pn->node_ub >= n->node_ub);
	}
	assert (fabs(n->node_lb) < 1e19);
	assert (fabs(n->node_ub) < 1e19);

	int n999 = 0;
	if (n->left == -999) n999++;
	if (n->right == -999) n999++;
	assert (n999 == 0 || n999 == 2);

	if (n->node_lb != -999 || n->node_ub != -999) {
		(*count) += 1;
	}

	if (n->left != -999)
		do_count_node(node_idx, level+1, n->left, count);
	if (n->right != -999)
		do_count_node(node_idx, level+1, n->right, count);

}

void BVH::ComputeMortonCodes() {
	std::unordered_set<int> occ; // 
	mc.clear();
	const int N = int(points_transformed.size());
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
	printf("%lu points, %lu unique Morton Codes\n", points_transformed.size(), mc.size());
}

void BVH::Build_Karras() {
	ComputeMortonCodes();
	const int N = int(mc.size());
	nodes.resize(N*2 - 1);
	internal_node_idx = 0;
	leaf_node_idx = N-1;

	struct timeval t0, t1;
	gettimeofday(&t0, NULL);
	radixsort<unsigned int>(mc.data(), int(mc.size()), &mapping);
	gettimeofday(&t1, NULL);
	float ms = (t1.tv_sec - t0.tv_sec)*1000.0f +
	           (t1.tv_usec- t0.tv_usec)/1000.0f;
	printf("Radix Sort: %f ms\n", ms);

	assert (N >= 2);

	gettimeofday(&t0, NULL);

	std::unordered_set<int> test;

	for (int i=0; i<=N-2; i++) {
		int delta_left  = (i-1<0) ?  -1 : clz(mc.at(i-1) ^ mc.at(i)),
		    delta_right = (i+1>N-1)? -1 : clz(mc.at(i) ^ mc.at(i+1)),
			d = (delta_right > delta_left) ? 1 : -1,
			delta_min = -999;
		if (i-d >= 0 && i-d <= N-1) {
			delta_min = clz(mc.at(i) ^ mc.at(i-d));
		} else { delta_min = -1; }

		// Find the other end
		int j = i;
		int dj = 1*d;

		// Binary Search idiom
		while (j+dj >= 0 && j+dj <= N-1) {
			const int d1 = clz(mc.at(i) ^ mc.at(j+dj));
			if (d1 >= delta_min) dj *= 2; else break;
		}
		int dj_tot = 0;
		while (true) {
			dj /= 2;
			if (j+dj_tot+dj >= 0 && j+dj_tot+dj <= N-1) {
				const int d1 = clz(mc.at(i) ^ mc.at(j+dj_tot+dj));
				if (d1 >= delta_min) { 
					dj_tot += dj;
				}
			}
			if (dj <= 1 && dj >= -1) break;
		}
		j = j+dj_tot;

//		while (mc.at(j-d) == mc.at(j) && nodes.at(j+N-1).parent!=-999) j-=d;

		#if 0
			int j_ref = i;
			while (j_ref+d >= 0 && j_ref+d <= N-1) {
				const int d1 = clz(mc.at(i) ^ mc.at(j_ref+d));
				if (d1 >= delta_min) { j_ref = j_ref+d; }
				else { break; }
			}
			if (not (j == j_ref)) {
				printf("j=%d j_ref=%d\n", j, j_ref);
				assert(0);
			}
			j = j_ref;
		#endif

		// Find the Split Point
		int lb = std::min(i, j), ub = std::max(i, j), split = lb;
		int delta_min_ij = clz(mc.at(lb) ^ mc.at(ub));
		#if 0
		for (int x = lb; x <= ub; x ++) {
			int delta = clz(mc.at(x) ^ mc.at(lb));
			if (clz(mc.at(x) ^ mc.at(lb)) > delta_min_ij) { split = x; }
			else { break; }
		}
		#endif

		{
			int lb1 = lb, ub1 = ub, split1;
			while (lb1+1 < ub1) {
				int mid = (lb1 + ub1) / 2;
				if (clz(mc.at(mid) ^ mc.at(lb)) > delta_min_ij) lb1 = mid;
				else ub1 = mid;
			}
			split1 = ub1;
			while (split1-1 >= lb &&
				clz(mc.at(split1) ^ mc.at(lb)) <= delta_min_ij) {
				split1 --;
			}
			#if 0
			if (not (split1 == split)) {
				printf("%d vs %d; lb=%d, ub=%d, lb1=%d, ub1=%d\n", split1, split,
					lb, ub, lb1, ub1);
				assert(0);
			}
			#endif
			split = split1;
		}

		BVHNode* node = &(nodes[i]);
		int left = split, right = split + 1;
		assert (split >= lb && split <= ub);

		node->node_lb = lb; node->node_ub = ub;

		if (lb == ub) { // Internal node must have at least 2 nodes !!
			printf("i=%d [LB,UB]=%d,%d split=%d\n", i, lb, ub, split);
			assert(0);
		} else {
			if (left == lb) {
				node->left = left + N-1;
				nodes.at(left + N-1).parent = i;
				test.insert(left);
			} else {
				node->left = left;
				nodes.at(left).parent = i;
			}
			
			if (right == ub) {
				node->right = right + N-1;
				nodes.at(right + N-1).parent = i;
				test.insert(right);
			} else {
				node->right = right;
				nodes.at(right).parent = i;
			}
		}
	}

	printf("%lu leaves touched\n", test.size());
	//assert (test.size() == N);

	gettimeofday(&t1, NULL);
	ms = (t1.tv_sec - t0.tv_sec)*1000.0f +
	     (t1.tv_usec - t0.tv_usec)/1000.0f;
	printf("BVH Hierarchy Generation: %f ms\n", ms);

	// Compute Bounding Boxes

	// bottom-up
	gettimeofday(&t0, NULL);

	// Remap
	std::vector<Vec3> points_sorted;
	std::vector<std::pair<int, int> > mc2sortedidx;  // an MC corresponds to which
	                                                 //   sorted indices
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

	std::unordered_set<int> n2v;
	for (int i=0; i<N; i++) {
		const int x = i+N-1;
		n2v.insert(x); 
		// If we have triangles then we need to expand triangles here
		nodes[x].min_ext = nodes[x].max_ext = points_sorted[mc2sortedidx.at(i).first];
	}

	for (int i=0; i<N-1; i++) {
		nodes.at(i).node_lb = mc2sortedidx.at(nodes.at(i).node_lb).first;
		nodes.at(i).node_ub = mc2sortedidx.at(nodes.at(i).node_ub).second;
	}
	for (int i=N-1; i<2*N-1; i++) {
		nodes.at(i).node_lb = mc2sortedidx.at(i-(N-1)).first;
		nodes.at(i).node_ub = mc2sortedidx.at(i-(N-1)).second;
	}

	while (!n2v.empty()) {
		std::unordered_set<int> n2v_next;
		for (const int x : n2v) {
			BVHNode* node = &(nodes[x]);
			const int par = node->parent;
			if (par != -999) {
				BVHNode* parent = &(nodes[par]);
				parent->min_ext = Vec3Min(parent->min_ext, node->min_ext);
				parent->max_ext = Vec3Max(parent->max_ext, node->max_ext);
				Vec3 m0 = nodes[par].min_ext,
				     m1 = nodes[par].max_ext;
				if (fabs(Vec3::dot(m0, m0)) > 1e20) {
					printf("X: (%f,%f,%f)->(%f,%f,%f)\n",
						nodes[x].min_ext.x,
						nodes[x].min_ext.y,
						nodes[x].min_ext.z,
						nodes[x].max_ext.x,
						nodes[x].max_ext.y,
						nodes[x].max_ext.z);
					assert(0);
				}
				if (fabs(Vec3::dot(m1, m1)) > 1e20) {
					printf("X: (%g,%g,%g)->(%g,%g,%g)\n",
						nodes[x].min_ext.x,
						nodes[x].min_ext.y,
						nodes[x].min_ext.z,
						nodes[x].max_ext.x,
						nodes[x].max_ext.y,
						nodes[x].max_ext.z);
					assert(0);
				}
				if (par == 0)
				printf("%d -> %d. BB:(%g,%g,%g)-->(%g,%g,%g) \n", x, par,
					nodes[par].min_ext.x,
					nodes[par].min_ext.y,
					nodes[par].min_ext.z,
					nodes[par].max_ext.x,
					nodes[par].max_ext.y,
					nodes[par].max_ext.z
				);
				assert (x!=par);
				n2v_next.insert(par);
			} else {
			}
		}
		n2v = n2v_next;
	}
	gettimeofday(&t1, NULL);
	ms = (t1.tv_sec - t0.tv_sec)*1000.0f +
	     (t1.tv_usec - t0.tv_usec)/1000.0f;
	printf("Remap and BB Computation: %f ms\n", ms);

	int count = 0;
	do_count_node(-999, 0, 0, &count);
	printf("Node Count = %d\n", count);
}

void BVH::PrintSummary() {
	printf("[BVH Summary]\n");
	int count = 0;
	do_count_node(-999, 0, 0, &count);
	printf("Node Count = %d\n", count);
	
}
