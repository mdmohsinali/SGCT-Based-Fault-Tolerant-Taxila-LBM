/* SGCT-Based-Fault-Tolerant-Taxila-LBM Code source file.
   Copyright (c) 2015, Md Mohsin Ali. All rights reserved.
   Licensed under the terms of the BSD License as described in the LICENSE_FT_CODE file.
   This comment must be retained in any redistributions of this source file.
*/

// 3D (x,y,z)  sparse grid data structure
// Storage is contiguous in the x co-ordinate
// written by Mohsin Ali, June 15 (base class is written by Peter Strazdins, June 15)

/* This data structure is for a distributed, (partially) filled sparse grid.
   This arises from the truncated SGCT formula. Let g be the grid index
   and l be the level, with g >= l. The filled sparse grid will
   have a factor of 2^(g-l) elements filled in over the classic sparse grid
   (in the places where this is possible).
 */

#ifndef FTSPARSEGRID3D_INCLUDED
#define FTSPARSEGRID3D_INCLUDED

// #define NOOPT2D  /* remove 2D B=1 optimization on interpolate() */

//#include "SparseGrid3D.h"
#include <stdio.h>
#include <assert.h>
#include <cmath>  //std::abs, log2 
#include <string> //std::string
#ifdef _OPENMP
#include <omp.h>
#endif

#include "Vec3D.h"
#include "ProcGrid3D.h"
#include "FTHaloArray3D.h"
#include "SparseGrid3D.h"

class FTSparseGrid3D : public SparseGrid3D {
  public:
  MPI_Comm myCommWorld;

  private:
  int sz_u;           // total number of elements in s.g. for this process 

  bool isPower2(int v) {
    return (v & (v-1)) == 0;
  }
  
  Vec3D<int> gridSz(Vec3D<int> g) {
    return Vec3D<int>(gridSz(g.x), gridSz(g.y), (is2d == 0) + (1 << g.z));
  }
  
  int numRightZeroes(int jG) {
    assert (jG > 0);
    int nz = 0;
    while ((jG & 1) == 0) 
      jG >>= 1, nz++;
    return(nz);
  }

  public:
  
  public:
  int gridSz(int g) {
    return (g < 0)? 0: ((1 << g) + 1);
  }
  
  // Constructors
  // First constructor
  FTSparseGrid3D(MPI_Comm comm, FTHaloArray3D *uh_, bool is2dim, 
	       Vec3D<int> g_,  ProcGrid3D *pg_, int blk=1) : SparseGrid3D(uh_, is2dim, 
	       g_, pg_, blk) { 
    myCommWorld = comm;
    allocedUh = false;
    useFullGrid = true; is2d = is2dim;  g = g_; pg = pg_;  B = blk;
    rx = 0; rs = 0; u = 0; // to trap bad references by clients
    ny = 0; cs = 0;
    uh = uh_;
    sz_u = uh->l.prod();
    nz = uh->l.z;
  } //FTSparseGrid3D()

  FTSparseGrid3D() {}

  // Second constructor
  FTSparseGrid3D(MPI_Comm comm, Vec3D<int> l_, int nSpecies, bool useFG, int lv, bool is2dim,
	       Vec3D<int> g_,  ProcGrid3D *pg_, int blk=1) : SparseGrid3D(useFG, lv, is2dim,
	       g_, pg_, blk) {
    myCommWorld = comm;
    allocedUh =false;
    useFullGrid = useFG; level = lv; is2d = is2dim; g = g_; pg = pg_; B = blk;
    rx = 0; rs = 0; u = 0; uh = 0; // to trap bad references by clients
    ny = 0; cs = 0;
    if (useFullGrid) {
      allocedUh = true;
      uh = new HaloArray3D(Vec3D<int>(l_.x, l_.y, l_.z), Vec3D<int>(0), B);
      sz_u = uh->l.prod();
      nz = uh->l.z;
      return;
    }
    int Ny, Nz;
    sz_u = 0; nz = 0;
    if (pg->myrank < 0) // this process does not hold part of s.g
      return;

    assert (isPower2(pg->P.x) && (is2d || isPower2(pg->P.y))); 
	      // needed to calc. nx, ny
    if (is2d) {
      Nz = gridSz(g.z);
    }
    else {
      Nz = gridSz(g.z) * nSpecies - (nSpecies-1);
    }
    nz = is2d? 1:  pg->G2L(2, Nz);
    ny = new int[nz]; cs = new int [nz]; 
    rx = new int* [nz]; rs = new int* [nz];

    int kG = (g.z==0)? 0: pg->L2G0(2, Nz);
    for (int k=0; k < nz; k++,kG++) {
      int lk = numRightZeroes((2*kG) | (1 << level));
      assert (0 < lk  &&  lk <= level);
      assert (g.z > 0 || lk == level); //check for 2D case
      // calculate local length (y) for the sparse block distribution 
      if (is2d) { // 2D
        Ny = gridSz(g.y - level + lk) * nSpecies - (nSpecies-1); //global number of s.g. elements at kG
        ny[k] = (Ny >= pg->P.y)? pg->G2L(1, Ny) :
	  (pg->id.y % (pg->P.y / (Ny-1)) == 0);
      }
      else { // 3D
        Ny = gridSz(g.y - level + lk); //global number of s.g. elements at kG
        ny[k] = (Ny >= pg->P.y)? pg->G2L(1, Ny):
	  (pg->id.y % (pg->P.y / (Ny-1)) == 0);
      }
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
  } //FTSparseGrid3D()

  // Destructor
  ~FTSparseGrid3D(){}  

  // Private declaration of sz_u in base class causes this repetition
  int numElts() {
    return (sz_u);
  }

  // Private declaration of sz_u in base class causes this repetition
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

}; //FTSparseGrid3D

#endif /*FTSPARSEGRID3D_INCLUDED*/
