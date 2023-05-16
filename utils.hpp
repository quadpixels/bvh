// Count Leading Zero
#include "bvh.hpp"

#ifndef _UTILS_H
#define _UTILS_H

#include <assert.h>
#include <string.h>
#include <stdlib.h>

int clz(unsigned);
unsigned int ExpandBits(unsigned int);
int FindSplit_naive(unsigned int*, int, int);
int FindSplit_ref(unsigned int*, int, int);
int FindSplit(unsigned int*, int, int);
Vec3 Vec3Min(const Vec3&, const Vec3&);
Vec3 Vec3Max(const Vec3&, const Vec3&);
template<typename T> void radixsort(T*, int, int**);
template<> void radixsort<unsigned int>(unsigned int*, int, int**);

#endif
