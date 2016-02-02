/* ParSGCT Code source file.
Copyright (c) 2015 Peter Edward Strazdins. All rights reserved.
Licensed under the terms of the BSD License as described in the LICENSE_SGCT file.
This comment must be retained in any redistributions of this source file.
*/

// Parallel 3D Advection Class
//   based on codes by Brendan Harding
// written by Peter Strazdins, May 13

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "Advect3D.h"

#include "3dadv_omp_cg.cpp" 

Advect3D::Advect3D(Vec3D<int> gridSz, Vec3D<double> v, double timef, double CFL, 
        int meth, int verb, Timer *t, MPI_Comm myCW, int blk) {
  B = blk;  
  myCommWorld = myCW; 
  gridSize = gridSz;
  V = v; 
  delta.x = 1.0 / (gridSize.x - 1);  // implicitly on unit square    
  delta.y = 1.0 / (gridSize.y - 1); 
  delta.z = is2D()? 1.0: 1.0 / (gridSize.z - 1); 
  tf = timef;
  dt = CFL * std::min<double>(delta.x, std::min<double>(delta.y, delta.z));
  method = meth;
  verbosity = verb;
  timer = t;
}

  
void Advect3D::initGrid(HaloArray3D *u, ProcGrid3D *g) {
#pragma omp parallel for default(shared)
  for (int kj = 0; kj < u->l.z*u->l.y; kj++) {
    int k = kj / u->l.y, j = kj % u->l.y;
    double z = delta.z * (k + g->L2G0(2, gridSize.z));
    double y = delta.y * (j + g->L2G0(1, gridSize.y));
    for (int i=0; i < u->l.x; i++) {
      double x = delta.x  * (i + g->L2G0(0, gridSize.x)*B);
      Vh(u, i, j, k) = initialCondition(x, y, z, 0.0, V.x, V.y, V.z);
    } 
  }
}

double Advect3D::checkError(double t, HaloArray3D *u, ProcGrid3D *g) {
  double err = 0.0;
#pragma omp parallel for default(shared) reduction(+:err)
  for (int kj = 0; kj < u->l.z*u->l.y; kj++) {
    int k = kj / u->l.y, j = kj % u->l.y;
    double z = delta.z * (k + g->L2G0(2, gridSize.z));
    double y = delta.y  * (j + g->L2G0(1, gridSize.y));
    for (int i=0; i < u->l.x; i++) {
      double x = delta.x  * (i + g->L2G0(0, gridSize.x)*B);
      err += std::abs(Vh(u, i, j, k) - initialCondition(x, y, z, t, 
							V.x, V.y, V.z));
    } 
  }
  return (err);
}

void Advect3D::updateLW(HaloArray3D *u) {
  timer->start("updateLW", u->l.prod(), 1);
  if (!is2D()) 
      compute_LW_scheme(dt, u->u, u->l.x, u->l.y, u->l.z, V.x, V.y, V.z, 
			delta.x, delta.y, delta.z);
  else
    updateLW2D(u);
  timer->stop("updateLW");
} //updateLW() 

void Advect3D::updateGodunov(HaloArray3D *u) {
  timer->start("updateGodunov", u->l.prod(), 1);
  if (!is2D())
    compute_godunov_scheme(dt, u->u, u->l.x, u->l.y, u->l.z, V.x, V.y, V.z, 
			   delta.x, delta.y, delta.z);
  else
    updateGodunov2D(u);
  timer->stop("updateGodunov");
} //updateLW() 

void Advect3D::updateMacCormack(HaloArray3D *u) {
  timer->start("updateMacCormack", u->l.prod(), 1);
  if (!is2D())
    compute_MacCormack_scheme(dt, u->u, u->l.x, u->l.y, u->l.z, V.x, V.y, V.z, 
			      delta.x, delta.y, delta.z);
  else
    updateMacCormack2D(u);
  timer->stop("updateMacCormack");
} //updateLW() 

// values at grid rows/columns 0 and N-1 are indentical; therefore  
// boundaries must get their value the next innermost opposite row/column
void Advect3D::updateBoundary(HaloArray3D *u, ProcGrid3D *g) {
  assert(u->halo.x == 1  &&  u->halo.y == 1 && u->halo.z <= 1) ; 
  int lx = u->l.x, ly = u->l.y, lz = u->l.z, sx = u->s.x, sy = u->s.y;
  int hz = u->halo.z; 
  timer->start("updateBoundary", hz*lx*ly + ly*lz + lz*lx, 1);
  double *bufS, *bufR, *b; int i, j, k;
  MPI_Request req; MPI_Status stat;

  if (g->P.x == 1) {  
    for (k=hz; k < lz+hz; k++) {
      for (j=1; j < ly+1; j++) {
         for (int ib = 0; ib < B; ib++) { 
	    V(u, 0+ib, j, k)    = V(u, lx-B+ib, j, k);
         }  
      }
      for (j=1; j < ly+1; j++) {
         for (int ib = 0; ib < B; ib++) {  
            V(u, lx+B+ib, j, k) = V(u, 2*B+ib, j, k);
         }
      }
    }
  } else {
    bufS = new double[ly*lz*B]; bufR = new double[ly*lz*B];
 
    int xOffs = (g->id.x == g->P.x-1) ? u->l.x-B: u->l.x;
    for (k=hz, b=bufS; k < lz+hz; k++) {
      for (j=1; j < ly+1; j++) {
         for (int ib = 0; ib < B; ib++, b++) {
	    *b = V(u, xOffs+ib, j, k);
         }
      }
    }
    //printf("%d: comm right boundary %dx%d to %d\n",g->myrank, ly, lz, g->neighbour(+1,0));
    MPI_Isend(bufS, ly*lz*B, MPI_DOUBLE, g->neighbour(+1, 0), 0, g->comm, &req);
    MPI_Recv(bufR, ly*lz*B, MPI_DOUBLE, g->neighbour(-1, 0), 0, g->comm, &stat);
    for (k=hz, b=bufR; k < lz+hz; k++) {
      for (j=1; j < ly+1; j++) {
         for (int ib = 0; ib < B; ib++, b++) { 
	    V(u, 0+ib, j, k) = *b;
         }
      }
    }
    MPI_Wait(&req, &stat);
        
    xOffs = (g->id.x == 0) ? 2*B: 1*B;
    for (k=hz, b=bufS; k < lz+hz; k++) {
      for (j=1; j < ly+1; j++) {
         for (int ib = 0; ib < B; ib++, b++) {  
	    *b = V(u, xOffs+ib, j, k);
         }
      }
    }
    // printf("%d: comm left boundary %dx%d to %d\n",g->myrank, ly, lz, g->neighbour(-1,0));
     MPI_Isend(bufS, ly*lz*B, MPI_DOUBLE, g->neighbour(-1, 0), 0, g->comm, &req);
     MPI_Recv(bufR, ly*lz*B, MPI_DOUBLE, g->neighbour(+1, 0), 0, g->comm, &stat);
     for (k=hz, b=bufR; k < lz+hz; k++) {
       for (j=1; j < ly+1; j++) {
          for (int ib = 0; ib < B; ib++, b++) { 
	     V(u, lx+B+ib, j, k) = *b;
          }
       }
     }
     MPI_Wait(&req, &stat);

     delete[] bufS; delete[] bufR;
  }
   
  if (g->P.y == 1) {
    for (k=hz; k < lz+hz; k++) {
      for (i=0; i < sx; i++)
	V(u, i, 0, k)    = V(u, i, ly-1, k);
      for (i=0; i < sx; i++)
	V(u, i, ly+1, k) = V(u, i, 2, k);
    }
  } else {
    bufS = new double[sx*lz]; bufR = new double[sx*lz]; 

    int yOffs = (g->id.y == g->P.y-1) ? ly-1: ly;
    for (k=hz, b=bufS; k < lz+hz; k++)
      for (i=0; i < sx; i++, b++) 
	*b = V(u, i, yOffs, k);
    //printf("%d: comm bottom boundary to %d\n",g->myrank,g->neighbour(+1,1));
     MPI_Isend(bufS, sx*lz, MPI_DOUBLE, g->neighbour(+1, 1), 0, g->comm, &req);
     MPI_Recv(bufR, sx*lz, MPI_DOUBLE, g->neighbour(-1, 1), 0, g->comm, &stat);
     for (k=hz, b=bufR; k < lz+hz; k++)
       for (i=0; i < sx; i++, b++)
	 V(u, i, 0, k) = *b;
     MPI_Wait(&req, &stat);
        
     yOffs = (g->id.y == 0) ? 2: 1;
     for (k=hz, b=bufS; k < lz+hz; k++)
       for (i=0; i < sx; i++, b++) 
	 *b = V(u, i, yOffs, k);
    //printf("%d: comm top boundary to %d\n", g->myrank, g->neighbour(-1, 1));
     MPI_Isend(bufS, sx*lz, MPI_DOUBLE, g->neighbour(-1, 1), 0, g->comm, &req);
     MPI_Recv(bufR, sx*lz, MPI_DOUBLE, g->neighbour(+1, 1), 0, g->comm, &stat);
     for (k=hz, b=bufR; k < lz+hz; k++)
       for (i=0; i < sx; i++, b++)
	 V(u, i, ly+1, k) = *b;
     MPI_Wait(&req, &stat);

     delete[] bufS; delete[] bufR;
  }

  if (hz == 0) { // no halo in z dimension
    timer->stop("updateBoundary");
    return;
  }
  if (g->P.z == 1) {
    for (j=0; j < sy; j++) {
      for (i=0; i < sx; i++)
	V(u, i, j, 0)    = V(u, i, j, lz-1);
      for (i=0; i < sx; i++)
	V(u, i, j, lz+1) = V(u, i, j, 2);
    }
  } else {
    bufS = new double[sx*sy]; bufR = new double[sx*sy]; 

    int zOffs = (g->id.z == g->P.z-1) ? lz-1: lz;
    for (j=0, b=bufS; j < sy; j++)
      for (i=0; i < sx; i++, b++) 
	*b = V(u, i, j, zOffs);
    //printf("%d: comm bottom boundary %dx%d to %d\n",g->myrank,sx, sy, g->neighbour(+1,1));
     MPI_Isend(bufS, sx*sy, MPI_DOUBLE, g->neighbour(+1, 2), 0, g->comm, &req);
     MPI_Recv(bufR, sx*sy, MPI_DOUBLE, g->neighbour(-1, 2), 0, g->comm, &stat);
     for (j=0, b=bufR; j < sy; j++)
       for (i=0; i < sx; i++, b++)
	 V(u, i, j, 0) = *b;
     MPI_Wait(&req, &stat);
        
     zOffs = (g->id.z == 0) ? 2: 1;
     for (j=0, b=bufS; j < sy; j++)
       for (i=0; i < sx; i++, b++) 
	 *b = V(u, i, j, zOffs);
     //printf("%d: comm top boundary %dx%d to %d\n", g->myrank, sx, sy, g->neighbour(-1, 1));
     MPI_Isend(bufS, sx*sy, MPI_DOUBLE, g->neighbour(-1, 2), 0, g->comm, &req);
     MPI_Recv(bufR, sx*sy, MPI_DOUBLE, g->neighbour(+1, 2), 0, g->comm, &stat);
     for (j=0, b=bufR; j < sy; j++)
       for (i=0; i < sx; i++, b++)
	 V(u, i, j, lz+1) = *b;
     MPI_Wait(&req, &stat);

     delete[] bufS; delete[] bufR;
  }
  timer->stop("updateBoundary");
} //updateBoundary()

double Advect3D::simulateAdvection(HaloArray3D *u, ProcGrid3D *g, double dtA) {
  timer->start("simulateAdvection", u->l.prod(), 0);
  double t = 0.0;  
  int s = 0;
  if (is2D()) // check that u, g have been set up appropriately
    assert (u->halo.z == 0  &&  u->l.z <= 1  &&  g->P.z == 1);
  while (t < dtA) {
    if (method == 0)
      updateGodunov(u);
    else if (method == 1)
      updateLW(u);
    else
      updateMacCormack(u);
    updateBoundary(u, g);

    t += dt; s++;
    if (verbosity > 3) {
      char s[64];
      sprintf(s, "after time %.4e, field is:\n", t);
      u->print(g->myrank, s);
    }
  }
  timer->stop("simulateAdvection");
  return t;
} //simulateAdvection()


