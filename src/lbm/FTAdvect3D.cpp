/* SGCT-Based-Fault-Tolerant-Taxila-LBM Code source file.
   Copyright (c) 2015, Md Mohsin Ali. All rights reserved.
   Licensed under the terms of the BSD License as described in the LICENSE_FT_CODE file.
   This comment must be retained in any redistributions of this source file.
*/

// Parallel 3D Advection Class
//   based on codes by Brendan Harding

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "FTAdvect3D.h"
#include <iostream>
#include<sys/stat.h>
#include<sys/types.h>

FTAdvect3D::FTAdvect3D(Vec3D<int> gridSz, Vec3D<double> v, double timef, double CFL, 
	   int meth, int verb, Timer *timer, MPI_Comm myCW, int blk) 
           : Advect3D(gridSz, v, timef, CFL, meth, verb, timer, myCW, blk) {
  myCommWorld = myCW;
}


void FTAdvect3D::initGridLBM(bool useFullGrid, FTSparseGrid3D * uvsg, HaloArray3D *u, double * doubleLbmFi, int lbmFiSize, 
                  MPI_Comm comm, int spec, int specIndex, bool isComponentGrid,
  bool is2D){
  int rank;
  MPI_Comm_rank(comm, &rank);
  int doubleLbmFiSize, doubleFiSizeWithoutSpec, startDoubleFiIndex, endDoubleFiIndex;
  int adv2DIndexX = 0, adv2DIndexY = 0,
      index2DPart = 0, adv3DIndexX = 0, adv3DIndexY = 0, adv3DIndexZ = 0;  
  int lz, ly, lx;
  if (is2D) { // 2D
     lz = u->l.z;
     if (u->l.z != 1)
        printf("***** Bad value. l.z = %d (expect 1)\n", u->l.z);
     assert(u->l.z == 1);
  }
  else // 3D
     lz = ((u->l.z % 2) == 0)? u->l.z: u->l.z - 1;  
  ly = ((u->l.y % 2) == 0)? u->l.y: u->l.y - 1;
  lx = (((u->l.x/B) % 2) == 0)? u->l.x: u->l.x - B; 
 
  if(lx == 0 || ly == 0){
     printf("[%d] lx or ly should not be zero.\n", rank);
     exit(1);
  } 

  // Count of double type elements in LBM's field array fi
  doubleLbmFiSize = lbmFiSize;

  if (isComponentGrid) { // component grid
     doubleFiSizeWithoutSpec = doubleLbmFiSize/spec;
     startDoubleFiIndex = doubleFiSizeWithoutSpec * specIndex;
     endDoubleFiIndex = startDoubleFiIndex + doubleFiSizeWithoutSpec - 1;
  }
  else { // full grid
     startDoubleFiIndex = 0;
     endDoubleFiIndex = doubleLbmFiSize - 1;
  }

  if (is2D) { // 2D
     #pragma omp parallel for default(shared)
     for(int i = startDoubleFiIndex; i <= endDoubleFiIndex; i++){
        adv2DIndexX = (i-startDoubleFiIndex) % lx;
        adv2DIndexY = (i-startDoubleFiIndex) / lx;
        Vh(u, adv2DIndexX, adv2DIndexY, 0) = doubleLbmFi[i];       
     }
  }
  else { // 3D
     #pragma omp parallel for default(shared)
     for(int i = startDoubleFiIndex; i <= endDoubleFiIndex; i++){
        adv3DIndexZ = (i-startDoubleFiIndex) / (lx*ly);
        index2DPart = (i-startDoubleFiIndex) % (lx*ly);
        adv3DIndexX = index2DPart % lx;
        adv3DIndexY = index2DPart / lx;
        Vh(u, adv3DIndexX, adv3DIndexY, adv3DIndexZ) = doubleLbmFi[i];     
     }
  }

  // Padding with 0's
  if(lx != u->l.x){
     #pragma omp parallel for default(shared)
     for (int k = 0; k < u->l.z; k++) {
        for (int j = 0; j < u->l.y; j++) {
           for (int ib = 0; ib < B; ib++)
              Vh(u, lx+ib, j, k) = 0.0;          
        }
     }
  }
  
  // Padding with 0's
  if(ly != u->l.y){
     #pragma omp parallel for default(shared)
     for (int k = 0; k < u->l.z; k++) {
        for (int i = 0; i < u->l.x; i++)
           Vh(u, i, ly, k) = 0.0;        
     }
  }

  // Padding with 0's
  if(lz != u->l.z){
     #pragma omp parallel for default(shared)
     for (int j = 0; j < u->l.y; j++) {
        for (int i = 0; i < u->l.x; i++)
           Vh(u, i, j, lz) = 0.0;       
     }
  }
  /*
  // Writing full grid results to files
  if (not isComponentGrid) {
     // Create "full_grid_results" directory
     mkdir("full_grid_results", S_IRWXU|S_IRGRP|S_IXGRP);

     FILE *fp;
     int rank;

     MPI_Comm_rank(myCommWorld, &rank);
     char ckptFileName[128] = "full_grid_results/saved_proc_", procRank[128] = " ";
     sprintf(procRank,"%d.dat", rank);
     strcat(ckptFileName, procRank);

     // Open the file for writing
     fp = fopen(ckptFileName, "w");
     if (fp == NULL) {
        printf("Could not open %s file for writing.\n", ckptFileName);
        exit(0);
     }

     if (useFullGrid) {
        if (is2D)
           assert(u->l.z == 1);
#pragma omp parallel for default(shared)        
        for (int k = 0; k < u->l.z; ++k) {
           for (int j = 0; j < u->l.y; ++j) {
              for (int i = 0; i < u->l.x; ++i) {
                 fprintf(fp, "%le ", Vh(u, i, j, k));
              }
           }
        }
     }
     else {
        int Nx = uvsg->gridSz(uvsg->g.x);
        int xIndex = uvsg->pg->G2L(0, Nx);
        int yIndex;
        int zIndex = uvsg->nz;
        int iG = uvsg->pg->L2G0(0, Nx);
        int k3;

        for (int k = 0; k < zIndex; ++k) {
           yIndex = uvsg->ny[k];
           for (int j = 0; j < yIndex; ++j) {
              double *u1 = &(uvsg->u[uvsg->rx[k][j]]);
              k3 = 0;
              for (int i = 0; i < xIndex*B; ++i, ++iG) {
	         if (iG % uvsg->rs[k][j] == 0 && (k*yIndex*xIndex*B + j*xIndex*B + k3) < uvsg->numElts()) {
                    fprintf(fp, "%le ", Vh(u, k3, j, k));
                    ++k3;
                 }
              }
           }
        }
     }

     // Close the file
     fclose(fp);
  }
  */
} //initGridLBM (with HaloArray)

void FTAdvect3D::initGridLBM(bool useFullGrid, FTSparseGrid3D * uvsg, FTHaloArray3D *u, double * doubleLbmFi, int lbmFiSize, 
                  MPI_Comm comm, int spec, int specIndex, bool isComponentGrid,
  bool is2D){
  int rank;
  MPI_Comm_rank(comm, &rank);
  int doubleLbmFiSize, doubleFiSizeWithoutSpec, startDoubleFiIndex, endDoubleFiIndex;
  int adv2DIndexX = 0, adv2DIndexY = 0,
      index2DPart = 0, adv3DIndexX = 0, adv3DIndexY = 0, adv3DIndexZ = 0;  
  int lz, ly, lx;
  if (is2D) { // 2D
     lz = u->l.z;
     if (u->l.z != 1)
        printf("***** Bad value. l.z = %d (expect 1)\n", u->l.z);
     assert(u->l.z == 1);
  }
  else // 3D
     lz = ((u->l.z % 2) == 0)? u->l.z: u->l.z - 1;  
  ly = ((u->l.y % 2) == 0)? u->l.y: u->l.y - 1;
  lx = (((u->l.x/B) % 2) == 0)? u->l.x: u->l.x - B; 
 
  if(lx == 0 || ly == 0){
     printf("[%d] lx or ly should not be zero.\n", rank);
     exit(1);
  } 

  // Count of double type elements in LBM's field array fi
  doubleLbmFiSize = lbmFiSize;

  if (isComponentGrid) { // component grid
     doubleFiSizeWithoutSpec = doubleLbmFiSize/spec;
     startDoubleFiIndex = doubleFiSizeWithoutSpec * specIndex;
     endDoubleFiIndex = startDoubleFiIndex + doubleFiSizeWithoutSpec - 1;
  }
  else { // full grid
     startDoubleFiIndex = 0;
     endDoubleFiIndex = doubleLbmFiSize - 1;
  }

  if (is2D) { // 2D
     #pragma omp parallel for default(shared)
     for(int i = startDoubleFiIndex; i <= endDoubleFiIndex; i++){
        adv2DIndexX = (i-startDoubleFiIndex) % lx;
        adv2DIndexY = (i-startDoubleFiIndex) / lx;
        Vh(u, adv2DIndexX, adv2DIndexY, 0) = doubleLbmFi[i];       
     }
  }
  else { // 3D
     #pragma omp parallel for default(shared)
     for(int i = startDoubleFiIndex; i <= endDoubleFiIndex; i++){
        adv3DIndexZ = (i-startDoubleFiIndex) / (lx*ly);
        index2DPart = (i-startDoubleFiIndex) % (lx*ly);
        adv3DIndexX = index2DPart % lx;
        adv3DIndexY = index2DPart / lx;
        Vh(u, adv3DIndexX, adv3DIndexY, adv3DIndexZ) = doubleLbmFi[i];     
     }
  }

  // Padding with 0's
  if(lx != u->l.x){
     #pragma omp parallel for default(shared)
     for (int k = 0; k < u->l.z; k++) {
        for (int j = 0; j < u->l.y; j++) {
           for (int ib = 0; ib < B; ib++)
              Vh(u, lx+ib, j, k) = 0.0;          
        }
     }
  }
  
  // Padding with 0's
  if(ly != u->l.y){
     #pragma omp parallel for default(shared)
     for (int k = 0; k < u->l.z; k++) {
        for (int i = 0; i < u->l.x; i++)
           Vh(u, i, ly, k) = 0.0;        
     }
  }

  // Padding with 0's
  if(lz != u->l.z){
     #pragma omp parallel for default(shared)
     for (int j = 0; j < u->l.y; j++) {
        for (int i = 0; i < u->l.x; i++)
           Vh(u, i, j, lz) = 0.0;       
     }
  }
  /*
  // Writing full grid results to files
  if (not isComponentGrid) {
     // Create "full_grid_results" directory
     mkdir("full_grid_results", S_IRWXU|S_IRGRP|S_IXGRP);

     FILE *fp;
     int rank;

     MPI_Comm_rank(myCommWorld, &rank);
     char ckptFileName[128] = "full_grid_results/saved_proc_", procRank[128] = " ";
     sprintf(procRank,"%d.dat", rank);
     strcat(ckptFileName, procRank);

     // Open the file for writing
     fp = fopen(ckptFileName, "w");
     if (fp == NULL) {
        printf("Could not open %s file for writing.\n", ckptFileName);
        exit(0);
     }

     if (useFullGrid) {
        if (is2D)
           assert(u->l.z == 1);
#pragma omp parallel for default(shared)        
        for (int k = 0; k < u->l.z; ++k) {
           for (int j = 0; j < u->l.y; ++j) {
              for (int i = 0; i < u->l.x; ++i) {
                 fprintf(fp, "%le ", Vh(u, i, j, k));
              }
           }
        }
     }
     else {
        int Nx = uvsg->gridSz(uvsg->g.x);
        int xIndex = uvsg->pg->G2L(0, Nx);
        int yIndex;
        int zIndex = uvsg->nz;
        int iG = uvsg->pg->L2G0(0, Nx);
        int k3;

        for (int k = 0; k < zIndex; ++k) {
           yIndex = uvsg->ny[k];
           for (int j = 0; j < yIndex; ++j) {
              double *u1 = &(uvsg->u[uvsg->rx[k][j]]);
              k3 = 0;
              for (int i = 0; i < xIndex*B; ++i, ++iG) {
	         if (iG % uvsg->rs[k][j] == 0 && (k*yIndex*xIndex*B + j*xIndex*B + k3) < uvsg->numElts()) {
                    fprintf(fp, "%le ", Vh(u, k3, j, k));
                    ++k3;
                 }
              }
           }
        }
     }

     // Close the file
     fclose(fp);
  }
  */
} //initGridLBM (with FTHaloArray)


double FTAdvect3D::checkErrorLBM(HaloArray3D *u, FTHaloArray3D *v, bool is2D, bool isReadFromDisk) {
  double err = 0.0;
  int i, j, k; // originally they are declared inside the loop
  int lx = v->l.x, ly = v->l.y, lz = v->l.z;
  lx = ((lx/B)%2 == 0)?lx: lx-B;
  ly = (ly%2 == 0)?ly: ly-1;
  if (is2D)
     assert(lz == 1);
  else
     lz = (lz%2 == 0)?lz: lz-1;
  if (not isReadFromDisk) {
     #pragma omp parallel for default(shared) reduction(+:err)
     for (k = 0; k < lz; k++) {
        for (j = 0; j < ly; j++) {
           for (i = 0; i < lx; i++)
              err += std::abs(Vh(u, i, j, k) - Vh(v, i, j, k));           
        }
     }
  }
  else { // isReadFromDisk
     double fullGridVal = 0.0;
     FILE *fp;
     int rank, ret;
     MPI_Comm_rank(myCommWorld, &rank);
     char ckptFileName[128] = "full_grid_results/saved_proc_", procRank[128] = " ";
     sprintf(procRank,"%d.dat", rank);
     strcat(ckptFileName, procRank);

     // Open the file for reading
     fp = fopen(ckptFileName, "r");
     if (fp == NULL) {
        printf("Could not open %s file for reading.\n", ckptFileName);
        exit(0);
     }
     #pragma omp parallel for default(shared) reduction(+:err)
     for (int k = 0; k < v->l.z; k++) {
        for (int j = 0; j < v->l.y; j++) {
           for (int i = 0; i < v->l.x; i++) {
              ret = fscanf(fp, "%lf ", &fullGridVal);
              err += std::abs(Vh(u, i, j, k) - fullGridVal);
           }
        }
     }

     // Close the file
     fclose(fp);
  }
  return (err);
} //checkErrorLBM (interpolation on full grid)

double FTAdvect3D::checkErrorLBM(FTSparseGrid3D *u, FTHaloArray3D *v, bool is2D, bool isReadFromDisk) {
  double err = 0.0;

  if (u->useFullGrid) {
     return (checkErrorLBM(u->uh, v, is2D, isReadFromDisk)); 
  }
  int Nx = u->gridSz(u->g.x); 
  int xIndex = u->pg->G2L(0, Nx);
  int yIndex;
  int zIndex = u->nz;
  int iG = u->pg->L2G0(0, Nx);
  int k3;  
  double *u1;

  if (!is2D) {
     zIndex = (zIndex % 2 == 1)? zIndex-1: zIndex;
  }

  if (not isReadFromDisk) {  
#pragma omp parallel for default(shared) reduction(+:err)
     for (int k = 0; k < zIndex; ++k) {
        yIndex = u->ny[k];
        if (is2D) {
           assert(zIndex == 1);                 
           yIndex = (yIndex % 2 == 1)? yIndex-1: yIndex;
        }
        for (int j = 0; j < yIndex; ++j) {
           u1 = &(u->u[u->rx[k][j]]);
           k3 = 0;
           for (int i = 0; i < xIndex*B; ++i, ++iG) {
	      if (iG % u->rs[k][j] == 0 && (k*yIndex*xIndex*B + j*xIndex*B + k3) < u->numElts()) {
                 err += std::abs(u1[k3] - Vh(v, k3, j, k));
                 ++k3;
              }
           }
        }
     }
  }
  else { // isReadFromDisk
     double sparseGridVal = 0.0;
     FILE *fp;
     int rank, ret;
     MPI_Comm_rank(myCommWorld, &rank);
     char ckptFileName[128] = "full_grid_results/saved_proc_", procRank[128] = " ";
     sprintf(procRank,"%d.dat", rank);
     strcat(ckptFileName, procRank);

     // Open the file for reading
     fp = fopen(ckptFileName, "r");
     if (fp == NULL) {
        printf("Could not open %s file for reading.\n", ckptFileName);
        exit(0);
     }  
#pragma omp parallel for default(shared) reduction(+:err)
     for (int k = 0; k < zIndex; ++k) {
        yIndex = u->ny[k];
        if (is2D) {
           assert(zIndex == 1);                 
           yIndex = (yIndex % 2 == 1)? yIndex-1: yIndex;
        }
        for (int j = 0; j < yIndex; ++j) {
           double *u1 = &(u->u[u->rx[k][j]]);
           k3 = 0;
           for (int i = 0; i < xIndex*B; ++i, ++iG) {
	      if (iG % u->rs[k][j] == 0 && (k*yIndex*xIndex*B + j*xIndex*B + k3) < u->numElts()) {
                 ret = fscanf(fp, "%le ", &sparseGridVal);
                 err += std::abs(u1[i] - sparseGridVal);
                 ++k3;
              }
           }
        }
     }

     // Close the file
     fclose(fp);
  }
  
  return (err);
} //checkErrorLBM (interpolation on sparse grid)

void FTAdvect3D::updateBoundary(HaloArray3D *u, ProcGrid3D *g) {

  assert(u->halo.x == 1  &&  u->halo.y == 1 && u->halo.z <= 1) ; 
  int lx = u->l.x, ly = u->l.y, lz = u->l.z, sx = u->s.x, sy = u->s.y;
  int hz = u->halo.z; 
  double *bufS, *bufR, *b; int i, j, k;
  MPI_Status stat;

  // Call the base method
     Advect3D::updateBoundary(u, g);

  // On x-direction
  if (g->P.x == 1) {
    for (k=hz; k < lz+hz; k++) {
      for (j=1; j < ly+1; j++) {
         for (int ib = 0; ib < B; ib++) {
            V(u, lx+ib, j, k) = V(u, B+ib, j, k);
         }
      }
    }
  }
  else {
     bufS = new double[ly*lz*B]; bufR = new double[ly*lz*B];
     if (g->id.x == 0) {
        int xOffs = B;
        for (k=hz, b=bufS; k < lz+hz; k++) {
           for (j=1; j < ly+1; j++) {
              for (int ib = 0; ib < B; ib++, b++) {
                 *b = V(u, xOffs+ib, j, k);
              }
           }
        }
        // printf("%d: comm left boundary %dx%d to %d\n",g->myrank, ly, lz, g->neighbour(-1,0));
        MPI_Send(bufS, ly*lz*B, MPI_DOUBLE, g->neighbour(-1, 0), 0, g->comm);
     }
     if (g->id.x == g->P.x-1) {
        MPI_Recv(bufR, ly*lz*B, MPI_DOUBLE, g->neighbour(+1, 0), 0, g->comm, &stat);
        for (k=hz, b=bufR; k < lz+hz; k++) {
           for (j=1; j < ly+1; j++) {
              for (int ib = 0; ib < B; ib++, b++) {
                 V(u, lx+ib, j, k) = *b;
              }
           }
        }
     }      
     delete[] bufS; delete[] bufR;
  }

  // On y-direction
  if (g->P.y == 1) {
    for (k=hz; k < lz+hz; k++) {
      for (i=0; i < sx; i++)
	V(u, i, ly, k) = V(u, i, 1, k);
    }
  } else {
     bufS = new double[sx*lz]; bufR = new double[sx*lz]; 
     if (g->id.y == 0) {
        int yOffs = 1;
        for (k=hz, b=bufS; k < lz+hz; k++)
          for (i=0; i < sx; i++, b++) 
	    *b = V(u, i, yOffs, k);
        //printf("%d: comm top boundary to %d\n", g->myrank, g->neighbour(-1, 1));
        MPI_Send(bufS, sx*lz, MPI_DOUBLE, g->neighbour(-1, 1), 0, g->comm);
     }
     if (g->id.y == g->P.y-1) {
        MPI_Recv(bufR, sx*lz, MPI_DOUBLE, g->neighbour(+1, 1), 0, g->comm, &stat);
        for (k=hz, b=bufR; k < lz+hz; k++)
          for (i=0; i < sx; i++, b++)
	    V(u, i, ly, k) = *b;
     }
     delete[] bufS; delete[] bufR;
  }

  // On z-direction  
  if (hz == 0) { // no halo in z dimension
    return;
  }
  if (g->P.z == 1) {
    for (j=0; j < sy; j++) {
      for (i=0; i < sx; i++)
	V(u, i, j, lz) = V(u, i, j, 1);
    }
  } else {
    bufS = new double[sx*sy]; bufR = new double[sx*sy];     
    if (g->id.z == 0) {        
       int zOffs = 1;
       for (j=0, b=bufS; j < sy; j++)
         for (i=0; i < sx; i++, b++) 
	   *b = V(u, i, j, zOffs);
       //printf("%d: comm top boundary %dx%d to %d\n", g->myrank, sx, sy, g->neighbour(-1, 1));
       MPI_Send(bufS, sx*sy, MPI_DOUBLE, g->neighbour(-1, 2), 0, g->comm);
    }
    if (g->id.z == g->P.z-1) {
       MPI_Recv(bufR, sx*sy, MPI_DOUBLE, g->neighbour(+1, 2), 0, g->comm, &stat);
       for (j=0, b=bufR; j < sy; j++)
         for (i=0; i < sx; i++, b++)
	   V(u, i, j, lz) = *b;
    }
    delete[] bufS; delete[] bufR;
  }

} //updateBoundary()

