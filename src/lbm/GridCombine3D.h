/* ParSGCT Code source file.
Copyright (c) 2015 Peter Edward Strazdins. All rights reserved.
Licensed under the terms of the BSD License as described in the LICENSE_SGCT file.
This comment must be retained in any redistributions of this source file.
*/

// Parallel 3D Sparse Grid Combination class
// written by Peter Strazdins, Jun 14

#ifndef GRIDCOMBINE3D_INCLUDED
#define GRIDCOMBINE3D_INCLUDED

#include "Vec3D.h"
#include "Timer.h"
#include "ProcGrid3D.h"
#include "HaloArray3D.h"
#include "SparseGrid3D.h"
#include <algorithm>    /* sort() */
#include <cmath>        /* log2() */
#include <stdlib.h>     /* abs() */
#include <mpi.h>        /* MPI_Comm */

//#define USE_BSEND  //uncomment if want to use MPI_Bsend version

class GridCombine3D {
  public:

  // return size corresponding to co-ordinate gix
  int gridSz(int gix);
    
  // return size vector corresponding to co-ordinates gix & whether 2D
  Vec3D<int> gridSz(Vec3D<int> gix);

  static Vec3D<int> gridSz_(bool is2d, Vec3D<int> gix) {
    Vec3D<int> g;
    g.x = 1 + (1 << gix.x);
    g.y = 1 + (1 << gix.y);
    g.z = !is2d + (1 << gix.z);
    return (g);
  }

  // variant ignoring redundant points on bottom right boundaries
  // This corresponds to the scaling of grid points
  static Vec3D<int> gridSz1(Vec3D<int> gix) {
    Vec3D<int> g;
    g.x = (1 << gix.x);
    g.y = (1 << gix.y);
    g.z = (1 << gix.z);
    return (g);
  }

  // return default value of number of processes on next lower diagonal 
  // (to achieve load balance)
  static int getPD1(int pD0, bool fixedProcs) {
    return ((pD0+1)/2);
  }

  // for time being, assume fixed number of processes along hyperplanes
  // (fixed timestep)
  static int nProcs(bool is2D, int level, int pD[], bool fixedProcs) {
    assert (fixedProcs); 
    int nP = 0;
    for (int i=0; i < (is2D? 2: 3); i++) {
      int l = level - i;
      nP += pD[i]*(is2D? l: (l+1)*l/2) ;
    }
    return (nP);
  }

  // return an appropriate process grid for grid index g
  static Vec3D<int> getProcConfig(int nprocs, bool is2d, Vec3D<int> g) {
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int nDim = is2d ? 2: 3; // 2D case if former
    struct pairCmp { // container for nested functions
      static bool compRv( std::pair<int,int> i, std::pair<int,int> j) {
	return (i.second > j.second);
      }
      static bool compAo( std::pair<int,int> i,  std::pair<int,int> j) {
	return (i.first < j.first);
      }
    }; 
    Vec3D<int> P = Vec3D<int>(0, 0, 0);
    // C++ could not use  Vec3D<  std::pair<int,int> > due to its use of union
    std::pair<int,int> gp[3] = {{0, g.v[0]}, {1, g.v[1]}, {2, g.v[2]}};    
    // sort largest grid dim first
    std::sort(std::begin(gp), std::end(gp), pairCmp::compRv);
    if ((nprocs & (nprocs-1)) != 0) { // non-power of 2
      MPI_Dims_create(nprocs, nDim, P.v);
      if (nDim < 3) P.v[2] = 1;
    } else {
      int p = std::log2(nprocs); 
      // assign nprocs to largest grid dimensions first
      int dp = std::min(p, gp[0].second - gp[1].second);
      P.v[0] += dp; p -= dp;
      if (2*(gp[1].second - gp[2].second) <= p) 
	dp = 2*(gp[1].second - gp[2].second);
      else 
	dp = std::min(p, gp[1].second - gp[2].second);
      P.v[0] += dp/2, P.v[1] += dp/2;
      p -= dp;
      if ((p%nDim == 1) /*only P.v[0] will benefit from remainder below*/
	  || (P.v[0] > P.v[1])) /* o.w. they benefit equally; try to balance*/
	P.v[1] += dp%2;
      else
	P.v[0] += dp%2;
      // assign remainder evenly
      for (int i=0; i < nDim; i++) 
	P.v[i] += p/nDim + (i < p%nDim);
      // P should now be in proportion to gp[*].second, but may be too large
      if (nDim == 3 && P.v[2] > gp[2].second) { // redistribute to bigger 2
	int del = P.v[2] - gp[2].second;
	P.v[0] += del/2;
	P.v[1] += del/2 + del%2;
	P.v[2] = gp[2].second;
      }
      if (P.v[1] > gp[1].second) // redistribute to biggest
	P.v[0] += P.v[1] - gp[1].second, P.v[1] = gp[1].second;
      assert (P.v[0] <= gp[0].second);
      for (int i=0; i < 3; i++)   
	P.v[i] = (1 << P.v[i]);
    }
    for (int i=0; i < 3; i++)
      gp[i].second = P.v[i];
    // restore to original order of g
    std::sort(std::begin(gp), std::end(gp), pairCmp::compAo); 
    for (int i=0; i < 3; i++)
      P.v[i] = gp[i].second;
    assert (P.prod() == nprocs);  
    assert (nDim > 2 || P.v[2] == 1);
    return P;      
  } // getProcConfig()

  // return SG combination coefficient for hyperplane dRank
  double getCombCoeff(int dRank);    
  static const int GID_ALL = -1; //signifies all grid to be used in combination
  // only use grid gid in the combination (use GID_ALL to reset to all grids)
  void singleGridCombine(int gid);
  // override default combination coefficient with coeff (for non-zero vals)
  // useful for debugging, or when expanding a single grid to the full grid 
  void overrideCombCoeff(double coeff);
  void restoreCombCoeff(); 

  // return number of grids
  int nGrids();

  // return which grid id global rank is in
  int getGid(int rank);
    
  // return rank within grid for process with global rank
  int getGrank(int rank);

  // return co-ordinates of grid gid
  Vec3D<int> gridIx(int gid);

  // return hyperplane grid gid is on
  int diagRank(int gid);

  ProcGrid3D** pgs; // array of process grids for each component grid
                    // HierachDecomp3D needs to modify this
 protected:
  bool is2d;        // gxU.z==0 does not necessarily mean a 2D case  
  int pD[3], verbosity; 
  MPI_Comm pgSComm; // for component & sparse grid doing the combination 
  Timer *timer;
  int nprocs;       // total number of processes in pgs
  int myrank;       // global rank of this process; useful for debugging
  bool fixedProcs;  // whether number of process is fixed across diagonals
                    // in the grid index space. Otherwise it doubles as you
                    // go outwards from the middle
  int *gRanks, *gIds; //caches each process' rank within grid and grid index
  int *gRank0;  //each grid's rank in gSComm of 0th process
  Vec3D<int> *gIxs; //caches grid indices

  // if set, getCombCoeff(g)  returns overrideCoeff
  bool isOverrideCoeff; double overrideCoeff; 
  int gidSingle; // if !=GID_ALL, getCombCoeff(g) returns 0 if g!=gidSingle
    
 public:
  int level; 
  Vec3D<int> gxU; // this process' component grid's co-ordinates and index
  int gidU;
  ProcGrid3D* pgU; // this process' component grid's process grid

  Vec3D<int> gxS;  // co-ordinates for sparse grid
  ProcGrid3D* pgS; // process grid for sparse grid

  GridCombine3D(){}

  // create object for performing combination grids at level level
  // with pD[] determining the number of processes on each diagonal
  // It assumes all MPI processes in comm call this together.
  GridCombine3D(int level, int pD[], bool fixedProcs, int sgProcs, 
		bool is2dim, Vec3D<int> gridS, Timer *timer, MPI_Comm comm, 
		int verbosity=0);
  ~GridCombine3D();

  bool is2D(); // convenience method to determine if its a 2D grid

  // perform the whole operation
  void gatherScatter(HaloArray3D* u, SparseGrid3D* usg);

  // each stage, with options for just a single grid
  void gather(HaloArray3D* u, SparseGrid3D* usg);
  void scatter(HaloArray3D* u, SparseGrid3D* usg);

  // set up msg buffers for gather/scatter Send(), B is HaloArray3D block size
  void sendInit(bool gather, int B);

  // wait for buffers to be used and delete them
  void sendWait(bool gather);

  // current process sends its part of component grid (U) to the respective 
  // sparse grid (S) process
  void gatherSend(HaloArray3D* u);

  // and receive it back from the scatter
  void scatterRecv(HaloArray3D* u) ; 

  // for each component grid in pgs, gather contributions from 
  // respective processes to the sparse grid (S)
  void gatherRecv(SparseGrid3D* uS);

  // and scatter them back out
  void scatterSend(SparseGrid3D* uS);

private:
  double *selfBuf;   // used to buffer data sent to self on these functions
#ifdef USE_BSEND
  char *sendBuf;
  long int sendBufSz; 
  int *sendRanks; // array of ranks this process has sent data to
#else
  MPI_Request *msgSendRequest;
  double** msgSendBuff;
#endif
  int nSendMsgs, maxSendMsgs; 

  // helper functions to reduce common code duplication
  void gatherSendScatterRecv(bool send, HaloArray3D* u);
  void gatherRecvScatterSend(bool recv, SparseGrid3D* uS);

  // hierarchical decomposition area
 protected:
  bool hMode; // if hMode, in hierarchical decomp. mode for hierarchy hg
  Vec3D<int> hx; // index in grid space of (maximal, if coalesced) hierarchy
  Vec3D<int> gxUsave, gxSsave; // to saves non-hier. mode values
  ProcGrid3D* pgSsave; 
 public:
  // set hierarchical decomp. mode: perform combination on 
  // (possibly coalesced) hierarchy hx of corresponding grid size hgx
  // returns storage to perform the appropriate SGCT
  HaloArray3D * setHierComb(Vec3D<int> hx, Vec3D<int> hgx);
  // restore normal mode
  void resetHierComb();  

}; //GridCombine3D

#endif /*GRIDCOMBINE3D_INCLUDED*/
