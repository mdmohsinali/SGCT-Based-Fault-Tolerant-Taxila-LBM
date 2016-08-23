/* SGCT-Based-Fault-Tolerant-Taxila-LBM Code source file.
   Copyright (c) 2015, Md Mohsin Ali. All rights reserved.
   Licensed under the terms of the BSD License as described in the LICENSE_FT_CODE file.
   This comment must be retained in any redistributions of this source file.
*/

/************************************************************************************ 
 * File       : FailureRecovery.h                                                   *
 * Description: Contains function prototypes for reconstructing faulty              *
 *              communicator by performing in-order failed process                  *
 *              replacement, and other supporting functions to achieve this         *
 * Author     : Md Mohsin Ali (mohsin.ali<AT>anu.edu.au)                            *
 * Created    : August 2013                                                         *
 * Updated    : December 2015                                                       *
 * Help       : See README file                                                     *
 ************************************************************************************/

#ifndef FAILURERECOVERY_INCLUDED
#define FAILURERECOVERY_INCLUDED

// Header files
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "mpi.h"
#include <unistd.h> //getopt(), gethostname()
#ifndef NON_FT
#include "mpi-ext.h"
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

// Function prototypes
MPI_Comm communicatorReconstruct(MPI_Comm myCommWorld, int childFlag,
                int * listFails, int * numFails, int * numNodeFails,
                int sumPrevNumNodeFails, int argc, char ** argv, int verbosity);
int numProcsFails(MPI_Comm comm);
void mpiErrorHandler(MPI_Comm * comm, int *errorCode, ...);
void repairComm(MPI_Comm * broken, MPI_Comm * repaired, int iteration,
                int * listFails, int * numFails, int * numNodeFails,
                int sumPrevNumNodeFails, int argc, char ** argv, int verbosity);
int rankIsNotOnFailedList(int rank, int * failedList, int numFails);
char * getHostToLaunch(int hostfileLineIndex);
int getHostfileLastLineIndex(void);
int getSlots(void);
#endif /*FAILURERECOVERY_INCLUDED*/
