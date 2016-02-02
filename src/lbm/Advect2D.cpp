/* ParSGCT Code source file.
Copyright (c) 2015 Peter Edward Strazdins. All rights reserved.
Licensed under the terms of the BSD License as described in the LICENSE_SGCT file.
This comment must be retained in any redistributions of this source file.
*/

// Parallel 2D Advection Solvers
// based on codes by Brendan Harding
// written by Peter Strazdins, May 13 

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "Advect3D.h"

//redefine these macros from HaloArray3D.h for 2D
#undef V
#undef Vh
#define V(u, i, j) (*((u)->ix(i, j, 0)))
#define Vh(u, i, j) (*((u)->ix_h(i, j, 0)))

void Advect3D::updateGodunov2D(HaloArray3D *u) {
  HaloArray3D * f_x = new HaloArray3D(u->l, Vec3D<int>(0), B);  
  HaloArray3D * f_y = new HaloArray3D(u->l, Vec3D<int>(0), B);  
  int fdx = (V.x >= 0.0)? -1: +1,
      fdy = (V.y >= 0.0)? -1: +1;
#pragma omp parallel for default(shared)
  for (int j=0; j < u->l.y; j++)     
    for (int i=0; i < u->l.x; i++) {
      V(f_x, i, j) = V.x * (Vh(u, i, j) - Vh(u, i+fdx, j));
      V(f_y, i, j) = V.y * (Vh(u, i, j) - Vh(u, i, j+fdy));
    }
  double dtdx = -fdx * dt / delta.x, dtdy = -fdy * dt / delta.y;
#pragma omp parallel for default(shared)
  for (int j=0; j < u->l.y; j++)     
    for (int i=0; i < u->l.x; i++) {
      Vh(u, i, j) += - dtdx * V(f_x, i, j) - dtdy * V(f_y, i, j);
    }
  delete f_x;
  delete f_y;
} //updateGodunov2D()

void Advect3D::updateLW2D(HaloArray3D *u) {
  HaloArray3D * uh = new HaloArray3D(Vec3D<int>(u->s.x-1, u->s.y-1, 1),
				     Vec3D<int>(0), B);  
  double sx = 0.5 * V.x / delta.x, sy = 0.5 * V.y / delta.y;
#pragma omp parallel for default(shared)
  for (int j=0; j < uh->l.y; j++)     
    for (int i=0; i < uh->l.x; i++) {
      V(uh,i,j) = 0.25*(Vh(u,i,j) + Vh(u,i-1,j) + Vh(u,i,j-1) + Vh(u,i-1,j-1))
	 -0.5*dt*(sx*(Vh(u,i,j) + Vh(u,i,j-1) - Vh(u,i-1,j) - Vh(u,i-1,j-1)) +
		  sy*(Vh(u,i,j) + Vh(u,i-1,j) - Vh(u,i,j-1) - Vh(u,i-1,j-1)));
    }
 
  double dtdx = 0.5 * dt / delta.x, dtdy = 0.5 * dt / delta.y;
#pragma omp parallel for default(shared)
  for (int j=0; j < u->l.y; j++)     
    for (int i=0; i < u->l.x; i++) {
      Vh(u, i, j) +=  
	- dtdx * (V(uh,i+1,j+1) + V(uh,i+1,j) - V(uh,i,j) - V(uh,i,j+1))
	- dtdy * (V(uh,i+1,j+1) + V(uh,i,j+1) - V(uh,i,j) - V(uh,i+1,j));
    }

  delete uh;
} //updateLW2D()

void Advect3D::updateMacCormack2D(HaloArray3D *u) {
  HaloArray3D * up = new HaloArray3D(Vec3D<int>(u->s.x-1*B, u->s.y-1, 1),
				     Vec3D<int>(0), B);  
  double sx = V.x * dt / delta.x, sy = V.y * dt / delta.y;  
#pragma omp parallel for default(shared)
  for (int j=0; j < up->l.y; j++) { 
    for (int i=0; i < up->l.x; i++) {
      V(up,i,j) = Vh(u,i-1,j-1) - sx * (Vh(u,i,j-1) - Vh(u,i-1,j-1)) 
	- sy * (Vh(u,i-1,j) - Vh(u,i-1,j-1));	
    } 
  }
  sx = 0.5 * sx; sy = 0.5 * sy; 
#pragma omp parallel for default(shared)
  for (int j=0; j < u->l.y; j++)     
    for (int i=0; i < u->l.x; i++) {
      Vh(u, i, j) = 0.5 * (Vh(u, i, j) + V(up, i+1, j+1))
			   - sx * (V(up, i+1, j+1) - V(up, i, j+1))
			   - sy * (V(up, i+1, j+1) - V(up, i+1, j));
    }
  delete up;
} //updateMacCormack2D()

