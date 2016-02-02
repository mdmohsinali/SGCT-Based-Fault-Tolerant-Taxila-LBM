/* ParSGCT Code source file.
Copyright (c) 2015 Peter Edward Strazdins. All rights reserved.
Licensed under the terms of the BSD License as described in the LICENSE_SGCT file.
This comment must be retained in any redistributions of this source file.
*/

// Parallel 3D Advection Class
// based on 2D advection class
// written by Peter Strazdins, June 14

#include <stdio.h>
#include <cmath> // M_PI
#include "HaloArray3D.h"
#include "ProcGrid3D.h"
#include "Timer.h"

class Advect3D {
  public:
  double tf /*final time*/, 
         dt /*timestep*/;
  Vec3D<double> delta /*grid spacing*/, V /*advection velocity*/;
  Vec3D<int> gridSize /* number of points in grid*/;
  int method; // 1st or 2nd order method}
  int verbosity; // for debugging output
  Timer *timer;
  int B;           // size of individual elements
  MPI_Comm myCommWorld;

  Advect3D(){}

  Advect3D(Vec3D<int> gridSz, Vec3D<double> v, double timef, double CFL, 
	   int meth, int verb, Timer *timer, MPI_Comm myCW = MPI_COMM_WORLD, int blk = 1);

  virtual ~Advect3D(){}

  // update u for a timestep by various schemes
  void updateGodunov(HaloArray3D *u);
  void updateLW(HaloArray3D *u);  
  void updateMacCormack(HaloArray3D *u);

  void updateGodunov2D(HaloArray3D *u);
  void updateLW2D(HaloArray3D *u);
  void updateMacCormack2D(HaloArray3D *u);

  static int bestMethod(bool is2Dsolver) {
    return(is2Dsolver? 2 /*McCormack*/: 1 /*LW*/);
  }

  inline bool is2D() {
    return (gridSize.z == 1);
  }

  inline double initialCondition(double x, double y, double z, double t=0.0,
				 double vx=1.0, double vy=1.0, double vz=1.0) {
    x = x - vx*t;
    y = y - vy*t;
    z = z - vz*t;
    double u = std::sin(4.0*M_PI*x) * std::sin(2.0*M_PI*y) ;
    if (!is2D())
	u *= std::sin(6.0*M_PI*z);
    return u;
  }

  void initGrid(HaloArray3D *u, ProcGrid3D *g);

  double checkError(double t, HaloArray3D *u, ProcGrid3D *g);

  void updateBoundary(HaloArray3D *u, ProcGrid3D *g);

  double simulateAdvection(HaloArray3D *u, ProcGrid3D *g, double dtA);

};


