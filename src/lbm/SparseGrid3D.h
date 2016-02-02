/* ParSGCT Code source file.
Copyright (c) 2015 Peter Edward Strazdins. All rights reserved.
Licensed under the terms of the BSD License as described in the LICENSE_SGCT file.
This comment must be retained in any redistributions of this source file.
*/

// 3D (x,y,z)  sparse grid data structure 
// Storage is contiguous in the x co-ordinate 
// written by Peter Strazdins, June 15

/* This data structure is for a distributed, (partially) filled sparse grid.
   This arises from the truncated SGCT formula. Let g be the grid index
   and l be the level, with g >= l. The filled sparse grid will
   have a factor of 2^(g-l) elements filled in over the classic sparse grid
   (in the places where this is possible). 
 */

#ifndef SPARSEGRID3D_INCLUDED
#define SPARSEGRID3D_INCLUDED

// #define NOOPT2D  /* remove 2D B=1 optimization on interpolate() */

#include <stdio.h>
#include <assert.h>
#include <cmath>  //std::abs, log2 
#include <string> //std::string
#ifdef _OPENMP
#include <omp.h>
#endif

#include "Vec3D.h"
#include "ProcGrid3D.h"

class SparseGrid3D {
 public:
  int level;          // truncated SGCT level
  bool is2d;
  Vec3D<int> g;       // global s.g. index
  ProcGrid3D* pg;     // s.g. process grid vectors
  int B;              // size of individual elements
  bool useFullGrid;   // if true use uh, else use other fields
  HaloArray3D *uh;    // full grid (dense) data structure
  bool allocedUh;     // records if uh was allocated by constructor
                      // sparse grid data structure
  double *u;          // storage vector
  int nz;             // local s.g. length in z-dimension (=1 for 2D case)
  int *ny;            // local s.g. lengths (per z) in y-dimension
  int *cs;            // col (y) stride vector
  int **rx, **rs;     // row index vector (as in CSR format) and stride per z

 private:
  int sz_u;           // total number of elements in s.g. for this process 

  static inline bool isPower2(int v) {
    return (v & (v-1)) == 0;
  }

  static inline int gridSz(int g) {
    return (g < 0)? 0: ((1 << g) + 1);
  }

  inline Vec3D<int> gridSz(Vec3D<int> g) {
    return Vec3D<int>(gridSz(g.x), gridSz(g.y), (is2d == 0) + (1 << g.z));
  }

  static inline int numRightZeroes(int jG) {
    assert (jG > 0);
    int nz = 0;
    while ((jG & 1) == 0) 
      jG >>= 1, nz++;
    return(nz);
  }

 public:
  SparseGrid3D(HaloArray3D *uh_, bool is2dim, 
	       Vec3D<int> g_,  ProcGrid3D *pg_, int blk=1) { 
    allocedUh = false;
    useFullGrid = true; is2d = is2dim;  g = g_; pg = pg_;  B = blk;
    rx = 0; rs = 0; u = 0; // to trap bad references by clients
    ny = 0; cs = 0;
    uh = uh_;
    sz_u = uh->l.prod();
    nz = uh->l.z;
  } //SparseGrid3D()

  SparseGrid3D() {}

  SparseGrid3D(bool useFG, int lv, bool is2dim,
	       Vec3D<int> g_,  ProcGrid3D *pg_, int blk=1){
    allocedUh =false;
    useFullGrid = useFG; level = lv; is2d = is2dim; g = g_; pg = pg_; B = blk;
    rx = 0; rs = 0; u = 0; uh = 0; // to trap bad references by clients
    ny = 0; cs = 0;
    if (useFullGrid) {
      allocedUh = true;
      uh = new HaloArray3D(pg->G2L(gridSz(g)), Vec3D<int>(0), B);
      sz_u = uh->l.prod();
      nz = uh->l.z;
      return;
    }
    sz_u = 0; nz = 0;
    if (pg->myrank < 0) // this process does not hold part of s.g
      return;

    assert (isPower2(pg->P.x) && (is2d || isPower2(pg->P.y))); 
	      // needed to calc. nx, ny
    int Nz = gridSz(g.z);
    nz = is2d? 1:  pg->G2L(2, Nz);
    ny = new int[nz]; cs = new int [nz]; 
    rx = new int* [nz]; rs = new int* [nz];

    int kG = (g.z==0)? 0: pg->L2G0(2, Nz);
    for (int k=0; k < nz; k++,kG++) {
      int lk = numRightZeroes((2*kG) | (1 << level));
      assert (0 < lk  &&  lk <= level);
      assert (g.z > 0 || lk == level); //check for 2D case
      int Ny = gridSz(g.y - level + lk); //global number of s.g. elements at kG
      // calculate local length (y) for the sparse block distribution 
      ny[k] = (Ny >= pg->P.y)? pg->G2L(1, Ny):
	(pg->id.y % (pg->P.y / (Ny-1)) == 0);
      cs[k] = 1 << (level - lk); // =1 for the 2D case
      // this is needed for interopolate() in order to calc. sparse offsets.
      // should be true if isPower2(pg->P.y) is
      assert (Ny < pg->P.y  ||  pg->L2G0(1, gridSz(g.y)) % cs[k] == 0);
      rs[k] = new int[ny[k]];
      rx[k] = new int[ny[k]+1];
      rx[k][0] = k==0? 0: rx[k-1][ny[k-1]];
      int jG = pg->L2G0(1, Ny);
      for (int j=0; j < ny[k]; j++,jG++) {
	int lj = numRightZeroes((2*jG) | (1 << lk));
	assert (0 < lj  &&  lj <= lk);
	rs[k][j] = 1 << (level-lj);
	int Nx = gridSz(g.x - level + lj); // global # s.g. elements at kG,jG
	// calculate local length corresp. Nx for the sparse block distribution
	int nx = (Nx >= pg->P.x)? pg->G2L(0, Nx): 
	  (pg->id.x % (pg->P.x / (Nx-1)) == 0); 
	// this is needed for interopolate() in order to calc. sparse offsets
	// should be true if isPower2(pg->P.x) is
	assert (Nx < pg->P.x  ||  pg->L2G0(0, gridSz(g.x)) % rs[k][j] == 0);
	rx[k][j+1] = rx[k][j] + nx*B;
	sz_u += nx*B;
      } //for(j...)
    } //for(k...)
    if (sz_u > 0)
      u = new double[sz_u];
  } //SparseGrid3D()

  ~SparseGrid3D() {
    if (useFullGrid && allocedUh) {
      delete uh;
      return;
    }
    if (u != 0)
      delete[] u;
    if (rs != 0) {
      for (int k=0; k < nz; k++)
	delete[] rs[k];
      delete[] rs;
    }
    if (rx != 0) {
      for (int k=0; k < nz; k++)
	delete rx[k];
      delete[] rx;
    }
    if (cs != 0)
      delete[] cs;
    if (ny != 0)
      delete ny;
  } //~~SparseGrid3D

  int numElts() {
    return (sz_u);
  }

  void zero() {
    if (useFullGrid) {
      uh->zero();
      return;
    }
#pragma omp parallel for default(shared)    
    for (int i = 0; i < sz_u; i++) {
      u[i] = 0.0;
    }
  }

  /* in the following, we take advantage of 2 facts:
      1. every point in v must be represented in the s.g.
      2. v is aligned to the s.g. at offset i0
   */

  // simplified (hence optimized) case for 2D B=1
  void interpolate2DB1(double coeff, HaloArray3D *v, Vec3D<int> rV,
                   Vec3D<int> i0, Vec3D<int> n) {
#pragma omp parallel for default(shared)
    for (int j = 0; j < n.y; j++) {
      int jS = j+i0.y, sx = rs[0][jS];
      if (rx[0][jS+1] == rx[0][jS]) { // has no elements as rs[jS] > pg->P.x
        assert (sx > pg->P.x);
        continue;
      }
      // scale f.g. to s.g. indices; assumes L2G0(0, Nx) % rs[jS] == 0
      int xr0 = i0.x%sx, xOffs = xr0==0? 0: sx-xr0; // perform alignment
      int nx = (n.x - xOffs + sx-1) / sx;
      int i0x = (i0.x + sx-1) / sx;
      assert (nx == 0  ||  i0x + nx-1 <= rx[0][jS+1] - rx[0][jS]);

      double *uv = &u[rx[0][jS] + i0x];
      double y = (1.0 * j) / rV.y;  int iy = (int) y;
      double ry = y - iy, ryc = 1.0 - ry;
      for (int i=0; i < nx; i++) { // only interpolate on points in the s.g. 
        double x = (1.0 * (xOffs + i * sx)) / rV.x ;
        int ix = (int) x;
        double rx = x - ix, rxc = 1.0 - rx;
         double interpolant =
	   rxc*ryc * *(v->ix_h(ix,   iy,   0))   +
	   rx *ryc * *(v->ix_h(ix+1, iy,   0))   +
	   rxc*ry  * *(v->ix_h(ix,   iy+1, 0))   +
	   rx *ry  * *(v->ix_h(ix+1, iy+1, 0));
	 uv[i] += coeff * interpolant;
      } //for(i...)   
    } //for(j...)     
  } //interpolate2DB1()


  // interpolate coeff * v, where resolution(v)*rV = resolution(s.g.),
  // into local s.g. elements i0..i0+dn-1
  void interpolate(double coeff, HaloArray3D *v, Vec3D<int> rV,
		   Vec3D<int> i0, Vec3D<int> n) {
    if (useFullGrid) {
      uh->interpolate(coeff, v, rV, i0, n);
      return;
    }
#ifndef NOOPT2D
    if (is2d && B==1) {
      interpolate2DB1(coeff, v, rV, i0, n);
      return;
    }
#endif
    for (int k = 0; k < n.z; k++) {
      int kS = k+i0.z;
      assert (kS < nz);
      int sy = cs[kS];
      // scale f.g. to s.g. indices; assumes L2G0(1, Ny) % sy == 0
      int yr0 = i0.y % sy, yOffs = yr0==0? 0: sy-yr0; // perform alignment
      int ny_ = (n.y - yOffs + sy-1) / sy; 
      int i0y = (i0.y + sy-1) / sy;
      assert (ny == 0 ||  i0y + ny_-1 < ny[kS]); 

      double z = (1.0 * k) / rV.z;
      int iz = (int) z;
      double rz = z - iz, rzc = 1.0 - rz;
#pragma omp parallel for default(shared)
      for (int j = 0; j < ny_; j++) {
	int jS = j+i0y;
	if (rx[kS][jS+1] == rx[kS][jS]) { // has no elts. as rs[kS][jS]>pg->P.x
	  assert (rs[kS][jS] > pg->P.x);
	  continue;
	}
	int sx = rs[kS][jS];
	// scale f.g. to s.g. indices; assumes L2G0(0, Nx) % sx == 0
	int xr0 = i0.x % sx, xOffs = xr0==0? 0: sx-xr0; // perform alignment
	int nx = (n.x - xOffs + sx-1) / sx; 
	int i0x = (i0.x + sx-1) / sx;
	assert (nx == 0  ||  i0x + nx-1 < (rx[kS][jS+1] - rx[kS][jS])/B);

	double *uv = &u[rx[kS][jS] + i0x*B];
	double y = (1.0 * (yOffs + j * sy)) / rV.y;
	int iy = (int) y;
	double ry = y - iy, ryc = 1.0 - ry;
	for (int i=0; i < nx; i++) { // only interpolate on points in the s.g.
	  double x = (1.0 * (xOffs + i * sx)) / rV.x ;
	  int ix = (int) x;
	  double rx = x - ix, rxc = 1.0 - rx;
	  int iB = i*B, ixB = ix*B;
	  for (int ib=0; ib < B; ib++) {
	    double interpolant =
	      rxc*ryc*rzc *  *(v->ix_h(ixB+ib,   iy,   iz))   + 
	      rx *ryc*rzc *  *(v->ix_h(ixB+ib+B, iy,   iz));
	    if (rV.y > 1) 
	      interpolant +=
		rxc*ry*rzc * *(v->ix_h(ixB+ib,   iy+1, iz))  + 
		rx *ry*rzc * *(v->ix_h(ixB+ib+B, iy+1, iz));
	    if (rV.z > 1) {
	      interpolant +=
		rxc*ryc*rz  * *(v->ix_h(ixB+ib,   iy,   iz+1)) +
		rx *ryc*rz  * *(v->ix_h(ixB+B+ib, iy,   iz+1));
	      if (rV.y > 1)
		interpolant +=
                rxc*ry *rz  * *(v->ix_h(ixB+ib,   iy+1, iz+1)) +
                rx *ry *rz  * *(v->ix_h(ixB+B+ib, iy+1, iz+1));
	    }   
	    uv[iB+ib] += coeff * interpolant;
	    //	    printf("kSjS=%d,%d sy, rc=%.0f v=%.0f  x,y,z=%.1f,%.1f,%.1f \n", 
	    //	   kS, jS, rxc*ryc*rzc, *(v->ix_h(ixB+ib,   iy,   iz)), x, y, z);
	  }
	} //for(i...)
      } //for(j...)
    } //for(k...)
  } //interpolate() 

  // sample s.g. elements i0:(i0+dn-1)*rV:rV into v, 
  // where resolution(v)*rV = resolution(s.g.)
  double *sample(Vec3D<int> i0, Vec3D<int> rV, Vec3D<int> n) {
    if (useFullGrid) 
      return (uh->sample(i0, rV, n));
    double *buf = new double [n.prod()*B];
    assert(i0.z + (n.z-1)*rV.z < nz);
    for (int k=0; k < n.z; k++) {
      int kS = i0.z + k*rV.z;
      int sy = cs[kS];
      int jScale = rV.y / sy;   assert(jScale > 0);
      int jS = (i0.y + sy-1) / sy;
      assert(jS + (n.y-1)*jScale < ny[kS]);
#pragma omp parallel for default(shared)
      for (int j=0; j < n.y; j++) {
	double *v = &buf[(k*n.y + j) * n.x * B];
	int sx = rs[kS][jS];
	int iScale = (rV.x / sx) * B;   assert(iScale > 0);
	int iS = ((i0.x + sx-1) / sx) * B;
	double *uv = &u[rx[kS][jS]];
	assert (iS + iScale*(n.x-1) < rx[kS][jS+1] - rx[kS][jS]);
	for (int i=0;  i < n.x;  i++, iS+=iScale) {
	  for (int ib = 0; ib < B; ib++) {
	    *v = uv[iS + ib]; 
	    v++;
	  }
	}
	jS += jScale;
      } //for(j...)
    } //for(k...) 
    return buf;
  } //sample()
 
  void print(int rank, std::string label) {
    if (useFullGrid) {
      uh->print(rank, label);
      return;
    }
    if (rank < 0)
      return;
    for (int k=0; k < nz; k++) {
      if (label.c_str()[0] != 0)
	printf("%d: %s[%d]:\n", rank, label.c_str(), k);
      int Nx = gridSz(g.x), nx = pg->G2L(0, Nx);
      for (int j=0; j < ny[k]; j++) {
	printf("%d: %d", rank, rx[k][j]);
	double *v = &u[rx[k][j]];
	int iG = pg->L2G0(0, Nx);
	for (int i=0;  i < nx;  i++, iG++) 
	  for (int ib=0; ib < B; ib++) {
	    if (iG % rs[k][j] == 0) {
	      assert (v < &u[rx[k][j+1]]);
	      const char *pad = *v < 10? "  ": *v < 100? " ": "";
	      printf("%s%3.1f ", pad, *v);
	      v++; 
	    } else
	      printf("      ");
	  }
	printf("\n"); 
      }
      printf("\n"); 
    }
} //print()

}; //SparseGrid3D

#endif /*SPARSEGRID3D_INCLUDED*/
