/* ParSGCT Code source file.
Copyright (c) 2015 Peter Edward Strazdins. All rights reserved.
Licensed under the terms of the BSD License as described in the LICENSE_SGCT file.
This comment must be retained in any redistributions of this source file.
*/

// Simple 3D (x,y)  array with a halo class. 
// Storage is contiguous in the x co-ordinate 
// written by Peter Strazdins, May 13
//   updated for block arrays, Apr 14
//   updated for 3D, Jun 14

#ifndef HALOARRAY3D_INCLUDED
#define HALOARRAY3D_INCLUDED

#include <stdio.h>
#include <assert.h>
#include <cmath>  //std::abs, log2 
#include <string> //std::string
#ifdef _OPENMP
#include <omp.h>
#endif

#include "Vec3D.h"

class HaloArray3D {
public:
  double *u;
  Vec3D<int> l, s; // local size, storage size (= local + halo) 
  Vec3D<int> halo;
  int B;           // size of individual elements

  HaloArray3D(Vec3D<int> l_, Vec3D<int> h, int blk=1) {
    l = l_; B = blk;
    l.x = l.x*B; 
    assert (h.prod() >= 0);
    halo = h;
    s = l + Vec3D<int>(B, 1, 1) * 2*halo;
    // check for potential overflow; should change to Vec3D<long> in future 
    if (s.x > 0 && s.y > 0 && s.z > 0) // this process owns part of the grid
      assert(log2(s.x) + log2(s.y) + log2(s.z) < sizeof(int)*8 - 1);
    if (s.prod() > 0)
      u = new double[s.prod()];
    else
      u = 0;
  }

  ~HaloArray3D() {
    if (u != 0)
      delete[] u;
  }

  // to take into account B>1, i = i'B + ib, where i' is the logical 
  // i-index of block and iB is offset within block, 0 <= ib < B
  inline double *ix(int i, int j, int k) {
#if 0
    assert (0 <= i && i < s.x  &&  0 <= j && j < s.y  &&  0 <= k && k < s.z);
#endif
    return(&u[i + s.x * (j + s.y*k)]);
  }

  inline double *ix_h(int i, int j, int k) {
    return ix(i + halo.x*B, j + halo.y, k + halo.z);
  }

  void zero() {
#pragma omp parallel for default(shared)
    for (int kj = 0; kj < l.z*l.y; kj++) {
      int k = kj / l.y, j = kj % l.y;  
      for (int i=0; i < l.x; i++) {
	*ix_h(i, j, k) = 0.0;
      }
    }
  }

  double norm1() {
    double norm = 0.0;
    for (int k=0; k < l.z; k++)
      for (int j=0; j < l.y; j++) 
	for (int i=0; i < l.x; i++) {
	  norm += std::abs(*ix_h(i, j, k));
	}
    return (norm);
  }

  double *pack(Vec3D<int> i0, Vec3D<int> m) {
    double* buf = new double [m.prod()*B];
#pragma omp parallel for default(shared)
    for (int kj = 0; kj < m.z*m.y; kj++) {
      int k = kj / m.y, j = kj % m.y;
      double *b = &buf[kj * m.x*B];
      for (int i=0; i < m.x*B; i++) {
	*b = *ix_h(i+i0.x*B, j+i0.y, k+i0.z);
	b++;
      }
    }
    return (buf); 
  }

  void unpack(double* buf, Vec3D<int> i0, Vec3D<int> m) {
#pragma omp parallel for default(shared)
    for (int kj = 0; kj < m.z*m.y; kj++) {
      int k = kj / m.y, j = kj % m.y;
      double *b = &buf[kj * m.x*B];
      for (int i=0; i < m.x*B; i++) {
	*ix_h(i+i0.x*B, j+i0.y, k+i0.z) = *b;
	b++;
      }
    }
  }

  void interpolate(double coeff, HaloArray3D *v, Vec3D<int> rV,
		   Vec3D<int> i0, Vec3D<int> n) {
#if 0
    printf("v.h=%d,%d,%d v.s=%d,%d,%d i0=%d,%d,%d n=%d,%d,%d\n", 
	   v->halo.x, v->halo.y, v->halo.z, v->s.x, v->s.y, v->s.z,
	   i0.x, i0.y, i0.z, n.x, n.y, n.z);
#endif
#pragma omp parallel for default(shared)
    for (int kj = 0; kj < n.z*n.y; kj++) {
      int k = kj / n.y, j = kj % n.y;
      double z = (1.0 * k) / rV.z;
      int iz = std::min((int) z, n.z-1);
      double rz = z - iz, rzc = 1.0 - rz;
      double y = (1.0 * j) / rV.y;
      int iy = std::min((int) y, n.y-1);
      double ry = y - iy, ryc = 1.0 - ry;
      for (int i=0; i < n.x; i++) {
	double x = (1.0 * i) / rV.x;
	int ix = std::min((int) x, n.x-1);
	double rx = x - ix, rxc = 1.0 - rx;
	int iB = (i+i0.x)*B, ixB = ix*B;
	for (int ib=0; ib < B; ib++) {
	  // avoid potentially undefined iy+1, iz+1 elts when ry, rz == 0
	  // if (ixB+B+ib >= 65) printf("oob on v B=%d\n", B);
	  double interpolant =
	    rxc*ryc*rzc * *(v->ix_h(ixB+ib,   iy,   iz))   + 
	    rx *ryc*rzc * *(v->ix_h(ixB+B+ib, iy,   iz));
	  if (rV.y > 1) 
	    interpolant +=
	      rxc*ry *rzc * *(v->ix_h(ixB+ib,   iy+1, iz))   + 
	      rx *ry *rzc * *(v->ix_h(ixB+B+ib, iy+1, iz));
	  if (rV.z > 1) { 
	    interpolant +=
	      rxc*ryc*rz  * *(v->ix_h(ixB+ib,   iy,   iz+1)) + 
	      rx *ryc*rz  * *(v->ix_h(ixB+B+ib, iy,   iz+1));
	    if (rV.y > 1)
	      interpolant +=
		rxc*ry *rz  * *(v->ix_h(ixB+ib,   iy+1, iz+1)) +
		rx *ry *rz  * *(v->ix_h(ixB+B+ib, iy+1, iz+1));
	  }
	  *ix_h(iB+ib, j+i0.y, k+i0.z) += coeff * interpolant;
	}
      }
    }
  } //interpolate()
 
  double *sample(Vec3D<int> i0, Vec3D<int> rV, Vec3D<int> n) {
    double *buf = new double [n.prod()*B];
#pragma omp parallel for default(shared)
    for (int kj = 0; kj < n.z*n.y; kj++) {
      int k = kj / n.y, j = kj % n.y;
      double *v = &buf[kj * n.x * B];
      for (int i=0; i < n.x; i++) {
	int iB = (i0.x + i*rV.x)*B;
	for (int ib=0; ib < B; ib++) {
	  *v = *ix_h(iB + ib, i0.y + j*rV.y, i0.z + k*rV.z);
	  v++;
	}
      }
    }
    return buf;
  } //sample()
 
  void print(int rank, std::string label) {
    printf("l.z=%d s.z=%d\n", l.z, s.z);
    for (int k=0; k < l.z; k++) {
      if (label.c_str()[0] != 0)
	printf("%d: %s[%d]:\n", rank, label.c_str(), k);    
      for (int j=0; j < l.y; j++) {
	if (rank >= 0)
	  printf("%d: ", rank);
	for (int i=0; i < l.x; i++) 
//	  printf("%0.2e ", *ix_h(i, j, k));
	  printf("(%d,%d)%0.2e ", i,j,*ix_h(i, j, k));
	printf("\n"); 
      }
      printf("\n"); 
    }
  }

  void printh(int rank, std::string label) {
    for (int k=0; k < l.z; k++) {
      if (label.c_str()[0] != 0)
	printf("%d: %s[%d]:\n", rank, label.c_str(), k);
      for (int j=0; j < s.y; j++) {
	if (rank >= 0)
	  printf("%d: ", rank);
	for (int i=0; i < s.x; i++) 
	  printf("%0.2e ", *ix(i, j, k));
	printf("\n"); 
      }
      printf("\n"); 
    }
  }
};

// macros for access elements for logical array without the halo
#define V(u, i, j, k) (*((u)->ix(i, j, k)))
// macros for accessing elements taking into account the halo
#define Vh(u, i, j, k) (*((u)->ix_h(i, j, k)))

#endif /*HALOARRAY3D_INCLUDED*/
