/* ParSGCT Code source file.
Copyright (c) 2015 Peter Edward Strazdins. All rights reserved.
Licensed under the terms of the BSD License as described in the LICENSE_SGCT file.
This comment must be retained in any redistributions of this source file.
*/

// Linear Grid class; used for testing grid combination algorithn
// written by Peter Strazdins, May 13

#include  <cmath>
#include "HaloArray3D.h"
#include "ProcGrid3D.h"

class LinGrid3D {
  public:
  Vec3D<int> gridSize; /* number of points in grid */
  Vec3D<int> scaleR; 
  int B;

  LinGrid3D(Vec3D<int> gSize, Vec3D<int> sR, int b=1) {
    gridSize = gSize;
    scaleR = sR;
    B = b;
  }

  ~LinGrid3D(){}

  void initGrid(HaloArray3D *u, ProcGrid3D *g) {
    bool is3D = (u->l.z < u->s.z);
    for (int k=0; k < u->l.z+is3D; k++) {
      double z = scaleR.z * (k + g->L2G0(2, gridSize.z));
      for (int j=0; j < u->l.y+1; j++) {
	double y = scaleR.y * (j + g->L2G0(1, gridSize.y));
	for (int i=0; i < u->l.x+B; i++) {
	  double x = scaleR.x * (i/B + g->L2G0(0, gridSize.x));
	  Vh(u, i, j, k) = x + y + z;
	}
      }
    }
  }

  double checkError(double scale, HaloArray3D *u, ProcGrid3D *g) {
    double err = 0.0;
    for (int k=0; k < u->l.z; k++) {
      double z = scaleR.z * (k + g->L2G0(2, gridSize.z));
      for (int j=0; j < u->l.y; j++) {
	double y = scaleR.y * (j + g->L2G0(1, gridSize.y));
	for (int i=0; i < u->l.x; i++) {
	  double x = scaleR.x  * (i/B + g->L2G0(0, gridSize.x));
	  err += std::abs(Vh(u, i, j, k) - scale*(x + y + z));
	}
      }      
    }
    return (err);
  }
};


