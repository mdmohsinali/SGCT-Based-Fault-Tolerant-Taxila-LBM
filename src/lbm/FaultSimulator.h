/* SGCT-Based-Fault-Tolerant-Taxila-LBM Code source file.
   Copyright (c) 2015, Md Mohsin Ali. All rights reserved.
   Licensed under the terms of the BSD License as described in the LICENSE_FT_CODE file.
   This comment must be retained in any redistributions of this source file.
*/

//  Simulator for fault-tolerance
// 
//  One or more processes (numProcsFails) fail during computation except process 0
// 
//  Input constraints: numProcsFails < numProcs
//  failureType = 0: non-real process failure
//  failureType = 1: real process failure
// 
//  Original Author: Knut Imar Hagen, June 2007
//  Modified by: Mohsin Ali, August 2013, February 2014
 

#ifndef FAULTSIMULATOR_INCLUDED
#define FAULTSIMULATOR_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include "mpi.h"
#include <signal.h>


int faultSimulate(int failureType, int * failedList, int numProcsFails, int myRank, int numProcs, int startProcs){

   if(numProcsFails >= numProcs){
      printf("[%d] FaultSimulator ERROR: numProcs should be > numProcsFails\n\n", myRank);

      MPI_Finalize();
      exit(0);
   }
   else if(numProcsFails > (numProcs-startProcs)){
      printf("[%d] FaultSimulator ERROR: numProcsFails should be <= (numProcs-startProcs)\n\n", myRank);

      MPI_Finalize();
      exit(0);
   }   

   // If the simulator shall fail some processes, let myRank 0 make a random
   // list of the failing processes and send the list to every process    
   int * failingProcs = NULL;
   int i = 0;
   
   //process 0 is not allowed to fail
   if(startProcs == 0){
       startProcs = 1;
   }

   if(numProcsFails > 0){
      failingProcs = (int*) calloc(numProcsFails, sizeof(int));
      if(myRank == 0){
         srand(time(NULL));
         int k;
         for(k = 0; k < numProcsFails; k++){
            int tmp = 0;
            int ok = 1;
            do{
               ok = 1;
               //tmp = rand()%numProcs; // between 0 to (numProcs-1) which is not used
               //tmp = rand()%(numProcs-1) + 1; // between 1 to (numProcs-1) which is supported by the following
               tmp = rand()%(numProcs-startProcs) + startProcs; // between startProcs to (numProcs-1)               
	       int l;
               for(l = 0; l < k; l++){
                  if(failingProcs[l] == tmp)
                     ok = 0;
               }
            }while(!ok);
            failingProcs[k] = tmp;
         }// end of for
      }// end of if(myRank == 0)
      MPI_Bcast(failingProcs, numProcsFails, MPI_INT, 0, MPI_COMM_WORLD);
   }// end of if(numProcsFails > 0)


   // If the number of loops is set higher than 0, the computation will run
   // the given number of times, and the processes listed to be failed
   // will fail in a sequentially order in a regular manner
   int countFail = 0;
   int upcountFail = 0;
   if(numProcsFails == 0){
      numProcsFails = 1;
   }  
   else{
      upcountFail = numProcs/numProcsFails;
   }

   if(failureType == 0){//non-real process failure   
       //failingProcs[0] = 33; failingProcs[1] = 37; //failingProcs[2] = 41; failingProcs[3] = 45; failingProcs[4] = 47; failingProcs[5] = 48;
       if(myRank == 0){ 
          printf("\n");
       }
       for(i = 0; i < numProcsFails; i++){
          failedList[i] = failingProcs[i]; 
          if(myRank == 0){             
             printf("[%d] ===== *NON_REAL* (SIMULATED) FAULT INJECT: process %d failed =====\n", myRank, failedList[i]);
          }
       }
       if(myRank == 0){ 
          printf("\n");
       }
   }
   else if(failureType == 1){//real process failure
       if(myRank == 0){ 
          printf("\n");
       }
       for(i = 0; i < numProcs; i++){
          if(numProcs/numProcsFails == upcountFail){
             if(countFail < numProcsFails && failingProcs[countFail++] == myRank){

                printf("[%d] ===== *REAL* FAULT INJECT: process %d failed =====\n", myRank, myRank);
                fflush(stdout);
                kill(getpid(), SIGKILL);
                //exit(0);
             }
             upcountFail=0;
          }
          upcountFail++;
       }
       if(myRank == 0){ 
          printf("\n");    
       }
   }

   // release memory
   free(failingProcs);

   /*
   /////////////////////
   if(myRank == 44){
      printf("\n===== *REAL* FAULT INJECT: process %d failed =====\n\n", myRank);
      fflush(stdout);
      kill(getpid(), SIGKILL);
      //exit(0);
   }
   /////////////////////
   */

   return 0;
}

#endif /*FAULTSIMULATOR_INCLUDED*/
