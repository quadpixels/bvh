#include "utils.hpp"

int clz(unsigned x) {
	int ret = 0;
	unsigned mask = 1U << 31;
	for (int i=0; i<32; i++, mask >>= 1) {
		if (!(x & mask)) { ret ++; } else break;
	}
	return ret;
}

// Expands a 10-bit integer into 30 bits
// by inserting 2 zeros after each bit.
unsigned int ExpandBits(unsigned int v) {
	v = (v * 0x00010001u) & 0xFF0000FFu;
	v = (v * 0x00000101u) & 0x0F00F00Fu;
	v = (v * 0x00000011u) & 0xC30C30C3u;
	v = (v * 0x00000005u) & 0x49249249u;
	return v;
}

int FindSplit_naive(unsigned int* sorted_mc, int lb, int ub) {
	int first = sorted_mc[lb], last = sorted_mc[ub];
	int common_prefix = clz(first ^ last);
	int ret = (lb + ub) >> 1;
	for (int i = lb; i <= ub; i++) {
		unsigned mid = sorted_mc[i];
		long cnt = clz(first ^ mid);
		if (cnt > common_prefix) { ret = i; }
	}
	return ret;
}

int FindSplit_ref(unsigned int* sorted_mc, int lb, int ub) {
	int first = sorted_mc[lb], last = sorted_mc[ub];
	if (first == last) return (lb + ub) >> 1;
	int common_prefix = clz(first ^ last);
	int split = lb, step = ub - lb;
	do {
		step = (step + 1) >> 1;
		int mid = split + step;
		if (mid < ub) {
			unsigned code = sorted_mc[mid],
			            c = clz(first ^ code);
			if (c > common_prefix) split = mid;
		}
	} while (step > 1);
	return split;
}

int FindSplit(unsigned int* sorted_mc, int lb, int ub) {
	int first = sorted_mc[lb], last = sorted_mc[ub];
	if (first == last) return (lb + ub) >> 1;
	int common_prefix = clz(first ^ last);
	int ll = lb, uu = ub;
	while (ll + 1 < uu) {
		int mm = (ll + uu) / 2;
		int mid = sorted_mc[mm],
		    c   = clz(first ^ mid);
		if (c   > common_prefix) { ll = mm; }
		else { uu = mm; }
	}
	int ret = ll, clz_ret = clz(first ^ sorted_mc[ret]);
	if (clz(first ^ sorted_mc[ll]) > common_prefix) ret = ll; else ret = uu;
	return ret;
}

Vec3 Vec3Min(const Vec3& x, const Vec3& y) {
	return Vec3(std::min(x.x, y.x),
	            std::min(x.y, y.y),
				std::min(x.z, y.z));
}

Vec3 Vec3Max(const Vec3& x, const Vec3& y) {
	return Vec3(std::max(x.x, y.x),
	            std::max(x.y, y.y),
				std::max(x.z, y.z));
}

template<>
void radixsort<unsigned int>(unsigned int* array, int size, int** mapping) {
	unsigned* temp = new unsigned[size];
	int* mapping1 = NULL;
	if (mapping) { 
		*mapping = new int[size]; 
		mapping1 = new int[size];
		for (int i=0; i<size; i++) (*mapping)[i] = i;
	}
	bool swapped = false;
	for (int bid=0; bid<sizeof(unsigned)*8; bid++) {
		unsigned bmask = (unsigned)(1 << bid);
		int nzeros = 0;
		for (int i=0; i<size; i++) {
			if ((array[i] & bmask) != bmask) nzeros ++;
		}
		int i0 = 0, i1 = nzeros;
		for (int i=0; i<size; i++) {
			if ((array[i] & bmask) != bmask) {
				temp[i0] = array[i];
				if (mapping) mapping1[i0] = (*mapping)[i];
				i0 ++;
			} else {
				temp[i1] = array[i];
				if (mapping) mapping1[i1] = (*mapping)[i];
				i1 ++;
			}
		}
		std::swap(array, temp);
		if (mapping) std::swap(*mapping, mapping1);
		swapped = !swapped;
	}
	if (swapped) {
		memcpy(array, temp, sizeof(unsigned)*size);
		if (mapping) memcpy(*mapping, mapping1, sizeof(int)*size);
	}
	delete[] temp;
	delete[] mapping1;
}
