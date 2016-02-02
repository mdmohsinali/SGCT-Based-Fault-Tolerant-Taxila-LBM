/*
   Copyright (c) 2014, Brendan Harding. All rights reserved.
*/

/*

Brendan Harding
October 2014

An installation on GNU GLPK is required.
This code has been tested with glpk-4.55 on both gnu and intel compilers.

Compiles and runs on my mac with: 
g++ -O2 -lglpk GcpGeneral.cpp
or
g++ GcpGeneral.cpp -O2 -lglpk
./a.out

Conclusions from 3d tests: 
	it can handle 969 grids in a 1/2 a second (level 16)
	2024 grids (level 21) takes about 20 seconds 
		(and still completes with success)
	2300 grids (level 22) takes about 29 seconds
	For most if not all practical purposes 1000 grids is about as high as one probably needs so it is probably okay

Original code by Brendan Harding: October 2014
Modified code by Mohsin Ali: November 2014

*/

#include <mpi.h>
#include "GcpGeneral.h"
#include "Vec3D.h"

// function prototypes
int coeffSelectUpdateAlgorithm3D(Vec3D<int> * v, int d, int n, int * failedGridList, 
	int failedGridCounter, double * listCoeffs, int verbose = 0);
bool isFailedGrid(int gridId, int * listGridFails, int numGridFails);

// function definitions
int coeffSelectUpdateAlgorithm3D(Vec3D<int> * v, int d, int n, int * failedGridList, 
	int failedGridCounter, double * listCoeffs, int verbose) {
	// n: level (definition varies)
	// d: dimension ,e.g., 3D

        int globalRank, ng = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &globalRank);

	if (d == 2) { //2D
	    ng = 4*(n+1) - 6;
	}
	else if(d == 3) { //3D
	    ng = (2*(n+1)*(n+1) - 4*(n+1) + 4);
	}
	else {
	    if(globalRank == 0)
		printf("Only supports 2D and 3D at the moment. Abort.\n");
	    exit(1);
	}

	int* grids = new int[d*ng];
	int ind = 0;
	#pragma omp parallel for default(shared)
	for (int i = 0; i < ng; ++i) {
	    grids[ind++] = v[i].x;
	    grids[ind++] = v[i].y;
	    grids[ind++] = v[i].z;
	}

	// optionally print the grids
	if (verbose > 0 && globalRank == 0) {
	    printf("Grids: (%d)\n", ng);
	    #pragma omp parallel for default(shared)
	    for (int i = 0; i < ng; ++i) {
		printf("\t");
		for (int k = 0; k < d; ++k) 
		    printf("%d ",grids[i*d+k]);
		printf("\n");
	    }
	}

	// initialise coefficient and status array
	double* coef = new double[ng];
	int* stat = new int[ng];
	
	for (int i = 0; i < ng; ++i) {
	    if (not isFailedGrid(i, failedGridList, failedGridCounter))
		stat[i] = 0; // grid is not failed
	    else
		stat[i] = 1; // grid is failed
	}
       
	// call wrapper to compute coefficients
	computeCoefficients(grids, ng, d, stat, coef);

	// copy the coefficients to return to the calle function  
	for (int i = 0; i < ng; ++i) {
	    listCoeffs[i] = coef[i];
	}

	// cleanup
	delete[] grids;
	delete[] coef;
	delete[] stat;

	return 0;
}


bool isFailedGrid(int gridId, int * listGridFails, int numGridFails) {

   int i;
   for(i = 0; i < numGridFails; ++i){
      if(gridId == listGridFails[i])
         return true;
   }

   if(i == numGridFails)
      return false;

   printf("Error in isFailedGrid() in GcpGeneral.cpp file.\n");
   exit(1);      
}//isFailedGrid()


