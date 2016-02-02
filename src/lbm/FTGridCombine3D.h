/* SGCT-Based-Fault-Tolerant-Taxila-LBM Code source file.
   Copyright (c) 2015, Md Mohsin Ali. All rights reserved.
   Licensed under the terms of the BSD License as described in the LICENSE_FT_CODE file.
   This comment must be retained in any redistributions of this source file.
*/

#include "Vec3D.h"
#include "Timer.h"
#include "ProcGrid3D.h"
#include "FTHaloArray3D.h"
#include "GridCombine3D.h"
#include "FTSparseGrid3D.h"
#include <algorithm>    /* sort() */
#include <cmath>        /* log2() */
#include <stdlib.h>     /* abs() */

class FTGridCombine3D : public GridCombine3D {
  public:

  // for time being, assume fixed numbner of processes along hyperplanes
  // (fixed timestep)
  static int nProcs(bool is2D, int level, int pD[4], bool fixedProcs, bool haveExtraGrids = false) {
    assert (fixedProcs); 
    int nP = 0;
    for (int i=0; i < (haveExtraGrids? 4: is2D? 2: 3); ++i) {
      int l = level - i;
      nP += pD[i]*(is2D? l: (l+1)*l/2) ;
    }
    return (nP);
  }

  // return SG combination coeffient for hyperplane dRank
  double getCombCoeff(int dRank);    

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

  public:
  int pD[4]; bool debug;
  bool haveExtraGrids;           // to identify whether extra grids are used for fault-tolerant or not
  double * listCoeffs;           // caches list of coefficients of grids
  int degOfParall;               // degree of parallelization over non-SGCT dimension

  // Constructor
  FTGridCombine3D(int level, int pD[4], bool fixedProcs, int sgProcs, 
		bool is2dim, Vec3D<int> gridS, Timer *timer, int verbosity = 0, 
		bool dbg = false, int degOfParall = 1, MPI_Comm comm = MPI_COMM_WORLD, 
                bool haveExtraGrids_ = false);

  // Destructor
  ~FTGridCombine3D();

  // Set grid id (gId) corresponding to 3D (or 2D) coordinate of grid of grids
  void setGridId(int x, int y, int z, int gId);

  // Return grid id corresponding to 3D (or 2D) coordinate of grid of grids
  int getGridId(int x, int y, int z);

  // Setting grid coordinates and grid ids
  void setGridCoordinatesGridIds();

  // procedure for alternate combination coefficients
  void alternateCombCoeffs(int* listFails, int numOfFails);

  // Procedure for the selection and update of combination 
  // coefficients in the presence of grid failure and with one extra layer
  void coeffSelectUpdate3D(int dimension, int * failedList = NULL, int numFailed = 0); 

  // Procedure for selecting and updating combination coefficients for 2D
  void coeffSelectUpdate2D(int * failedGridList, int numFailedGrid = 0);

  // perform the whole operation
  void gatherScatter(HaloArray3D* u, FTSparseGrid3D* usg);

  // each stage, with options for just a single grid
  void gather(HaloArray3D* u, FTSparseGrid3D* usg);
  void scatter(HaloArray3D* u, FTSparseGrid3D* usg);

  // set up msg buffers for gather/scatter Send(), B is HaloArray3D block size
  void sendInit(bool gather, int B);

  // wait for buffers to be used and delete them
  void sendWait(bool gather);

  // current process sends its part of component grid (U) to the respective 
  // sparse grid (S) process
  void gatherSend(HaloArray3D* u);

  // and receive it back from the scatter
  void scatterRecv(HaloArray3D* u); 

  // for each component grid in pgs, gather contributions from 
  // respective processes to the sprse grid (S)
  void gatherRecv(FTSparseGrid3D* uS);

  // and scatter them back out
  void scatterSend(FTSparseGrid3D* uS);

  // Procedure for converting global rank to rank of a combination communicator
  int localCombRank(int rank); 

  private:
  Vec3D<int> * gxG;         // caches 3D coordinates of grid of grids
  int * gxGtoGridId;        // caches grid ids of 3D coordinates of grid of grids

  double *selfBuf;   // used to buffer data sent to self on these functions
#ifdef USE_BSEND
  char *sendBuf;
  long int sendBufSz; 
#else
  MPI_Request *msgSendRequest;
  double** msgSendBuff;
  int nSendMsgs, maxSendMsgs; 
#endif

  // helper functions to reduce common code duplication
  void gatherSendScatterRecv(bool send, HaloArray3D* u);
  void gatherRecvScatterSend(bool recv, FTSparseGrid3D* uS);
};

