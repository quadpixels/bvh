#ifndef __BVH_H
#define __BVH_H

#include <vector>
#include <unordered_map>

#ifdef __CUDACC__
	#define __CUDA_CALLABLE__ __host__ __device
#else
	#define __CUDA_CALLABLE__
#endif

class Vec3 {
public:
	float x, y, z;
	Vec3(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}
	Vec3 operator+(const Vec3& other) const {
		return Vec3(x+other.x, y+other.y, z+other.z);
	}
	Vec3 operator*(const float a) const {
		return Vec3(x*a, y*a, z*a);
	}
	Vec3() : x(0), y(0), z(0) { }
	static float dot(const Vec3& a, const Vec3& b) {
		return a.x*b.x + a.y*b.y + a.z*b.z;
	}
	float operator[](int idx) const {
		switch (idx) {
			case 0: return x;
			case 1: return y;
			case 2: return z;
		}
		return 1.0f/0.0f;
	}
};


class BVH {
public:
	class BVHNode {
		public:
		int left,   right,    parent;
		int node_lb, node_ub;
		Vec3 min_ext, max_ext;
		BVHNode() { 
			left = right = parent = -999;
			node_lb = 1e9;
			node_ub = -1e9;
			min_ext = Vec3( 1e20,  1e20,  1e20);
			max_ext = Vec3(-1e20, -1e20, -1e20);
		}
	};
	int* mapping;
	std::vector<BVHNode>   nodes;
	std::vector<unsigned>  mc; // Morton Codes
	std::unordered_map<unsigned, std::vector<int> > mc2point;
	std::vector<Vec3> points, points_transformed;
	void Build();
	void Build_Karras(); // Using Karras's Method
	void ComputeMortonCodes();
	void PrintSummary();
	BVH(const std::vector<Vec3>& _p) :
		points(_p) {
		// Must be determined after computing Morton Codes
//		int N = int(_p.size());
//		nodes.resize(N*2 - 1);
//		internal_node_idx = 0;
//		leaf_node_idx = N-1;
	}
	BVH() {}
private:
	int internal_node_idx, leaf_node_idx;
	int do_build(int node_idx, int lb, int ub);
	void do_count_node(int parent, int level, int node_idx, int* count);
	int AllocInternalNode() {
		return internal_node_idx ++;
	}
	int AllocLeafNode() {
		return leaf_node_idx ++;
	}
};

#endif
