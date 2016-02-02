/* SGCT-Based-Fault-Tolerant-Taxila-LBM Code source file.
   Copyright (c) 2015, Md Mohsin Ali. All rights reserved.
   Licensed under the terms of the BSD License as described in the LICENSE_FT_CODE file.
   This comment must be retained in any redistributions of this source file.
*/

// Parallel 3D Advection Class
// based on 2D advection class

#ifndef FTHALOARRAY3D_INCLUDED
#define FTHALOARRAY3D_INCLUDED

#include "HaloArray3D.h"
#include <stdio.h>
#include <assert.h>
#include <cmath>  //std::abs, log2 
#include <string> //std::string
#ifdef _OPENMP
#include <omp.h>
#endif

#include "Vec3D.h"

class FTHaloArray3D : public HaloArray3D {
  public:
  MPI_Comm myCommWorld;

  // Constructor
  FTHaloArray3D(MPI_Comm comm, Vec3D<int> l_, Vec3D<int> h, int blk=1):HaloArray3D(l_, h, blk) {
     myCommWorld = comm;
  }

  // Destructor
  ~FTHaloArray3D(){}

  double norm1(bool isReadFromDisk = false) {
     if (not isReadFromDisk) { // call base class function if no saved data
        return (HaloArray3D::norm1());
     }
     else { // isReadFromDisk
        double norm = 0.0, fullGridVal = 0.0;
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

        for (int k=0; k < l.z; k++)
           for (int j=0; j < l.y; j++) 
              for (int i=0; i < l.x; i++) {
                 ret = fscanf(fp, "%lf ", &fullGridVal);
	         norm += std::abs(fullGridVal);
              }

        // Close the file
        fclose(fp);

        return (norm);
     }
  } // end of "norm1()"
  
};

#endif /*FTHALOARRAY3D_INCLUDED*/
