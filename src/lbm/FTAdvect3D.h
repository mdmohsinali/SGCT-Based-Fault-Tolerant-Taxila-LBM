/* SGCT-Based-Fault-Tolerant-Taxila-LBM Code source file.
   Copyright (c) 2015, Md Mohsin Ali. All rights reserved.
   Licensed under the terms of the BSD License as described in the LICENSE_FT_CODE file.
   This comment must be retained in any redistributions of this source file.
*/

// Parallel 3D Advection Class
// based on 2D advection class

#include <stdio.h>
#include <complex>
#include <cmath> // M_PI
#include "FTHaloArray3D.h"
#include "ProcGrid3D.h"
#include "Advect3D.h"
#include "Timer.h"
#include "FTSparseGrid3D.h"

class FTAdvect3D : public Advect3D {
  public:
  MPI_Comm myCommWorld;

  // Constructor
  FTAdvect3D(Vec3D<int> gridSz, Vec3D<double> v, double timef, double CFL, 
	   int meth, int verb, Timer *timer, MPI_Comm myCW = MPI_COMM_WORLD, 
           int blk = 1);

  // Destructor
  ~FTAdvect3D(){}

  void initGridLBM(bool useFullGrid, FTSparseGrid3D *uvsg, HaloArray3D *u, double * doubleLbmF1, int LbmSize,
                   MPI_Comm comm, int spec, int specIndex = 0, 
                   bool componentGrid = true, bool is2D = false);

  void initGridLBM(bool useFullGrid, FTSparseGrid3D *uvsg, FTHaloArray3D *u, double * doubleLbmF1, int LbmSize,
                   MPI_Comm comm, int spec, int specIndex = 0, 
                   bool componentGrid = true, bool is2D = false);

  double checkErrorLBM(HaloArray3D *u, FTHaloArray3D *v, bool is2D = false, 
                        bool isReadFromDisk = false);
  double checkErrorLBM(FTSparseGrid3D *u, FTHaloArray3D *v, bool is2D = false, 
                        bool isReadFromDisk = false);

  void updateBoundary(HaloArray3D *u, ProcGrid3D *g);

};


