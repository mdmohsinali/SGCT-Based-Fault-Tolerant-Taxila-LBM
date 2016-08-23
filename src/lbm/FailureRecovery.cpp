/* SGCT-Based-Fault-Tolerant-Taxila-LBM Code source file.
   Copyright (c) 2015, Md Mohsin Ali. All rights reserved.
   Licensed under the terms of the BSD License as described in the LICENSE_FT_CODE file.
   This comment must be retained in any redistributions of this source file.
*/

/*********************************************************************************** 
 * File       : FailureRecovery.cpp                                                *
 * Description: Contains functions for reconstructing faulty communicator by       *
 *              performing in-order failed process replacement, and other          *
 *              supporting functions to achieve this                               *
 * Author     : Md Mohsin Ali (mohsin.ali<AT>anu.edu.au)                           *
 * Created    : August 2013                                                        *
 * Updated    : December 2015                                                      *
 * Help       : See README file                                                    *
 ***********************************************************************************/

// Header files
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "mpi.h"
#ifndef NON_FT
#include "mpi-ext.h"
#endif
#include "FailureRecovery.h"
#include <unistd.h> //getopt(), gethostname()
#ifdef _OPENMP
#include <omp.h>
#endif

// Defined values
//#define SHRUNKEN_RECOVERY      // return shrunken comm as a reconstructed comm
#define HANG_ON_REMOVE           // defining this remove fault-tolerant mpi hang on
#define SLOTS                    getSlots()
//#define RECOV_ON_SPARE_NODES   // defining this causes processes will be spawned on
                                 // spare nodes (handling node failures),
                                 // otherwise, spawned on the same node where it was
                                 // before the failure (handling process failures)
//#define RUN_ON_COMPUTE_NODES   // defining this causes run on compute nodes,
                                 // otherwise, test run on head node/personal machine
#define GLOBAL_DETECTION         // defining this causes global detection, otherwise,
                                 // local detection
#define MERGE_TAG             20000
#define MERGE_TAG2            20001
#define NON_FAILED_FAILED_TAG 20002
#define PROCS_GRID_NF_TAG     20003
#define PROCS_GRID_NF_TAG2    20004

#ifndef NON_FT
///////////////////////////////////////////////////////////////////////////////////
MPI_Comm communicatorReconstruct(MPI_Comm myCommWorld, int childFlag, int * listFails,
        int * numFails, int * numNodeFails, int sumPrevNumNodeFails, int argc,
        char ** argv, int verbosity) {
	int i, ret, rank, nprocs, totFails = 0, * failedList, flag,
            iterCounter = 0, failure = 0;
	MPI_Comm parent, mcw, dupComm, tempIntracomm;
        MPI_Errhandler newEh;
        double startTime = 0.0, endTime;
#ifndef SHRUNKEN_RECOVERY
        int oldRank = 0, recvVal[2], length;
        MPI_Status mpiStatus;
        MPI_Comm unorderIntracomm;
        char hostName[MPI_MAX_PROCESSOR_NAME];
#endif //SHRUNKEN_RECOVERY

        // Error handler
        MPI_Comm_create_errhandler(mpiErrorHandler, &newEh);

	MPI_Comm_get_parent(&parent);
        MPI_Comm_rank(myCommWorld, &rank);
        startTime = MPI_Wtime();

        do{
           failure = 0;
           ret = MPI_SUCCESS;
           // Parent part
	   if(MPI_COMM_NULL == parent){
              if(iterCounter == 0)
                 mcw = myCommWorld;
              // Set error handler for communicator
              MPI_Comm_set_errhandler(mcw, newEh);
              // World information
              MPI_Comm_rank(mcw, &rank);
	      MPI_Comm_size(mcw, &nprocs);
              // Synchronize. Sometimes hangs on without this
#ifdef HANG_ON_REMOVE
              OMPI_Comm_agree(mcw, &flag); // to fix hang on problem
#endif             
              // Target function
              if(MPI_SUCCESS != (ret = MPI_Comm_dup(mcw, &dupComm))) {
              //if(MPI_SUCCESS != (ret = MPI_Barrier(mcw))) { // MPI_Comm_dup or
                                                              // MPI_Barrier
                 if(verbosity > 0 && rank == 0)
                    printf("[????? Process %d (nprocs %d)] MPI_Comm_dup (parent): "
                           "Unsuccessful (due to process failure) OK\n", rank, nprocs);
                 // Revoke the communicator
	         if(MPI_SUCCESS != (OMPI_Comm_revoke(mcw))){
                    if(rank == 0)
                       printf("[Process %d (nprocs %d)] Iteration %d: OMPI_Comm_revoke "
                              "(parent): Error!\n", rank, nprocs,  iterCounter);
                 }
                 else{
                    if(verbosity > 1 && rank == 0)
                       printf("[Process %d (nprocs %d)] Iteration %d: OMPI_Comm_revoke "
                              "(parent): SUCCESS\n", rank, nprocs, iterCounter);
                 }
                 // Call repair with splitted world
                 totFails = numProcsFails(mcw);
                 failedList = (int *) malloc(totFails*sizeof(int));
                 repairComm(&mcw, &tempIntracomm, iterCounter, failedList, numFails,
                            numNodeFails, sumPrevNumNodeFails, argc, argv, verbosity);
                 // Assign list of failed processes
#pragma omp parallel for default(shared)
                 for(i = 0; i < *numFails; ++i)
                    listFails[i] = failedList[i];
                 // Free memory
                 free(failedList);
                 // Operation failed: retry
                 failure = 1;
              } // end of "if MPI_Comm_dup fails"
              else{
                 if(verbosity > 0 && rank == 0)
                    printf("[..... Process %d (nprocs %d)] Iteration %d: MPI_Comm_dup "
                           "(parent): SUCCESS\n", rank, nprocs, iterCounter);
                 // Operation success: breaking iteration
                 failure = 0;
              }
	   } // end of "parent"
#ifndef SHRUNKEN_RECOVERY
           // Child part
	   else{
              MPI_Comm_set_errhandler(parent, newEh);
              // Synchronize. Sometimes hangs on without this
              // Position of code and intercommunicator, parent, (not intra) is
              // important
#ifdef HANG_ON_REMOVE
              OMPI_Comm_agree(parent, &flag);// to fix hang on problem
#endif              
              MPI_Comm_rank(parent, &rank);
              MPI_Comm_size(parent, &nprocs);

              if(verbosity > 0 && rank == 0){
                  MPI_Get_processor_name(hostName, &length);
                  printf("[Process %d, nprocs = %d] created on host %s (child)\n", 
                          rank, nprocs, hostName);
              }
              if(MPI_SUCCESS != (MPI_Intercomm_merge(parent, true, &unorderIntracomm))){
                 if(rank == 0)
                    printf("[Process %d] Iteration %d: MPI_Intercomm_merge (child): "
                           "Error!\n", rank, iterCounter);
              }
              else{
                 if(verbosity > 1 && rank == 0)
                    printf("[Process %d] Iteration %d: MPI_Intercomm_merge (child): "
                           "SUCCESS\n", rank, iterCounter);
              }
              // Receive failed ranks and number of fails from process 0 of parent
              if(MPI_SUCCESS != (MPI_Recv(&recvVal, 2, MPI_INT, 0, MERGE_TAG, 
                                unorderIntracomm, &mpiStatus))){
                 if(rank == 0)
                    printf("[Process %d] Iteration %d: MPI_Recv1 (child): Error!\n", 
                           rank, iterCounter);
              }
              else{
                 if(verbosity > 1 && rank == 0)
                    printf("[Process %d] Iteration %d: MPI_Recv1 (child): SUCCESS\n", 
                           rank, iterCounter);
                 oldRank = recvVal[0]; *numFails = recvVal[1];
              }              
              // Split the communicator to order the ranks.
              // No order is maintaining here. Actual ordering is done on parent side
              // This is a support only to parent side
              if(MPI_SUCCESS != (MPI_Comm_split(unorderIntracomm, 0, oldRank,
                 &tempIntracomm))) {
                 if(rank == 0)
                    printf("[Process %d] Iteration %d: MPI_Comm_split (child): " 
                           "Error!\n", rank, iterCounter);
              }
              else{
                 if(verbosity > 1 && rank == 0)
                    printf("[Process %d] Iteration %d: MPI_Comm_split (child): " 
                           "SUCCESS\n", rank, iterCounter);
              }
              // Operation on parent failed: retry
              ret = (!MPI_SUCCESS);
              failure = 1;
              // Free memory
              MPI_Comm_free(&unorderIntracomm);
              MPI_Comm_free(&parent);
	   }// end of "child"
#endif //SHRUNKEN_RECOVERY
           // Reset comm world
           if(ret != MPI_SUCCESS)
              mcw = tempIntracomm;
           // Reset parent value for parent
           if(parent == MPI_COMM_NULL && ret != MPI_SUCCESS)
              parent = mcw;
           // Reset parent value of child and make the operation collective
           if(MPI_SUCCESS != ret && MPI_COMM_NULL != parent)
              parent = MPI_COMM_NULL;
           ++iterCounter;
        }while(failure > 1);// replace 'failure > 1' with 'failure' if want fault
                            // tolerant recovery

        if(MPI_COMM_NULL == parent && childFlag == 0 && rank == 0){
           endTime = MPI_Wtime();
#ifndef SHRUNKEN_RECOVERY
           printf("[%d]----- Reconstructing failed communicator (including failed "
                  "list creation) takes %0.6f Sec (MPI_Wtime) -----\n", rank,
                  endTime - startTime);
#else
           printf("[%d]----- Reconstructing failed communicator with shrinking the "
                  "communicator (including failed list creation) takes %0.6f Sec "
                  "(MPI_Wtime) -----\n", rank, endTime - startTime);
#endif //SHRUNKEN_RECOVERY
        }
        // Memory release
        MPI_Errhandler_free(&newEh);

        return mcw;
}//communicatorReconstruct()

///////////////////////////////////////////////////////////////////////////////////
int numProcsFails(MPI_Comm mcw){
	int rank, ret, numFailures = 0, flag;
        MPI_Group fGroup;
        MPI_Errhandler newEh;
        MPI_Comm dupComm;

        // Error handler
        MPI_Comm_create_errhandler(mpiErrorHandler, &newEh);

        MPI_Comm_rank(mcw, &rank);

        // Set error handler for communicator
        MPI_Comm_set_errhandler(mcw, newEh);

        // Target function
        if(MPI_SUCCESS != (ret = MPI_Comm_dup(mcw, &dupComm))) {
        //if(MPI_SUCCESS != (ret = MPI_Barrier(mcw))) { // MPI_Comm_dup or MPI_Barrier
           OMPI_Comm_failure_ack(mcw);
           OMPI_Comm_failure_get_acked(mcw, &fGroup);
           // Get the number of failures
           MPI_Group_size(fGroup, &numFailures);
        }// end of "MPI_Comm_dup failure"

        OMPI_Comm_agree(mcw, &flag);
        // Memory release
	if(numFailures > 0)
           MPI_Group_free(&fGroup);
        MPI_Errhandler_free(&newEh);

        return numFailures;
}//numProcsFails()

///////////////////////////////////////////////////////////////////////////////////
void mpiErrorHandler(MPI_Comm * comm, int *errorCode, ...){
    MPI_Group failedGroup;

    OMPI_Comm_failure_ack(*comm);
    OMPI_Comm_failure_get_acked(*comm, &failedGroup);

    // *errorCode == MPI_ERR_PROC_FAILED classify failure type as "process" failure

    // Failed processes will NOT be synchronized Without delay
    // This delay MUST be through error handler (otherwise, problematic)
    usleep(10000); // 10 milliseconds delay (MPI_Comm_revoke is failed without
                   // this for a large number of processes)
    MPI_Group_free(&failedGroup);

    return;
}//mpiErrorHandler()

///////////////////////////////////////////////////////////////////////////////////
void repairComm(MPI_Comm * broken, MPI_Comm * repaired, int iteration, int * listFails,
        int * numFails, int * numNodeFails, int sumPrevNumNodeFails, int argc,
        char ** argv, int verbosity) {
	MPI_Comm tempShrink;
	int i, ret, result, procsNeeded = 0, oldRank, oldGroupSize, flag,
            * tempRanks, * failedRanks, * errCodes, rank, * procsNeededToLaunch;
	MPI_Group oldGroup, failedGroup, shrinkGroup;
        double startTime = 0.0, endTime, shrinkTime;
        char hostName[128], ** appToLaunch, *** argvToLaunch, ** hostNameToLaunch,
             * hostNameFailed;
        MPI_Info * hostInfoToLaunch;
#ifdef SHRUNKEN_RECOVERY
#undef RECOV_ON_SPARE_NODES // in case not disabled
#endif //SHRUNKEN_RECOVERY

#ifndef SHRUNKEN_RECOVERY
        MPI_Comm unorderIntracomm, tempIntercomm;
        int newRank, rankKey = 0, nprocs, j, * shrinkMergeList;
        double spawnTime, mergeTime, agreeTime;
#ifdef RECOV_ON_SPARE_NODES
#define RUN_ON_COMPUTE_NODES // in case this is not enabled
       int hostfileLastLineIndex, * failedNodeList = NULL, * nodeList = NULL,
           totNodeFailed = 0, failedNodeCounter;
#endif //RECOV_ON_SPARE_NODES
#endif //SHRUNKEN_RECOVERY

#ifdef RUN_ON_COMPUTE_NODES
        int hostfileLineIndex, tempLineIndex;
#endif //RUN_ON_COMPUTE_NODES

        gethostname(hostName, sizeof(hostName));
        MPI_Comm_rank(*broken, &rank);
        if(rank == 0)
            startTime = MPI_Wtime();

#ifndef GLOBAL_DETECTION
	MPI_Comm_size(*broken, &oldGroupSize);
	MPI_Comm_group(*broken, &oldGroup);
	MPI_Comm_rank(*broken, &oldRank);
	OMPI_Comm_failure_ack(*broken);
	OMPI_Comm_failure_get_acked(*broken, &failedGroup);
	MPI_Group_size(failedGroup, &procsNeeded);
	errCodes = (int *) malloc(sizeof(int) * procsNeeded);
	// Figure out ranks of the processes which had failed
	tempRanks = (int *) malloc(sizeof(int) * oldGroupSize);
	failedRanks = (int *) malloc(sizeof(int) * oldGroupSize);
#pragma omp parallel for default(shared)
	for(i = 0; i < oldGroupSize; ++i) 
	   tempRanks[i] = i;
	MPI_Group_translate_ranks(failedGroup, procsNeeded, tempRanks, oldGroup,
        	failedRanks);
#endif        
        shrinkTime = MPI_Wtime();
        // Shrink the broken communicator to remove failed procs
	if(MPI_SUCCESS != (ret = OMPI_Comm_shrink(*broken, &tempShrink)))
           printf("Iteration %d: OMPI_Comm_shrink (parent): ERROR!\n", iteration);
        else{
           if(verbosity > 1 )
              printf("Iteration %d: OMPI_Comm_shrink (parent): SUCCESS\n", iteration);
        }
        if (verbosity > 0 && rank == 0)
           printf("OMPI_Comm_shrink takes %0.6f Sec\n", MPI_Wtime() - shrinkTime);
#ifdef GLOBAL_DETECTION
	MPI_Comm_group(*broken, &oldGroup);
	MPI_Comm_group(tempShrink, &shrinkGroup);
	MPI_Comm_size(*broken, &oldGroupSize);
	MPI_Group_compare(oldGroup, shrinkGroup, &result);

	if(result != MPI_IDENT)
	   MPI_Group_difference(oldGroup, shrinkGroup, &failedGroup);
	MPI_Comm_rank(*broken, &oldRank);
	MPI_Group_size(failedGroup, &procsNeeded);
	errCodes = (int *) malloc(sizeof(int)*procsNeeded);
	// Figure out ranks of the processes which had failed
	tempRanks = (int*)malloc(sizeof(int)*oldGroupSize);
	failedRanks = (int*)malloc(sizeof(int)*oldGroupSize);
#pragma omp parallel for default(shared)
	for(i = 0; i < oldGroupSize; ++i)
	   tempRanks[i] = i;
	MPI_Group_translate_ranks(failedGroup, procsNeeded, tempRanks, oldGroup,
        	failedRanks);
	MPI_Group_free(&shrinkGroup);
#endif        
        // Assign number of failed processes
        *numFails = procsNeeded;
        hostNameToLaunch = (char **) malloc(procsNeeded * sizeof(char *));

        if(verbosity > 0 && rank == 0)
	       printf("*** Iteration %d: Application: Number of process(es) failed "
                      "in the corresponding communicator is %d ***\n", iteration,
                      procsNeeded);
        if(rank == 0){
            endTime = MPI_Wtime();
            printf("[%d]----- Creating failed process list takes %0.6f Sec (MPI_Wtime) "
                   "-----\n", rank, endTime - startTime);
        }
#ifndef SHRUNKEN_RECOVERY
#ifdef RUN_ON_COMPUTE_NODES
#ifdef RECOV_ON_SPARE_NODES
	// Determining total number of node failed, and a list of them
        hostfileLastLineIndex = getHostfileLastLineIndex(); // started from 0
	nodeList = (int *) malloc((hostfileLastLineIndex+1) * sizeof(int));
	memset(nodeList, 0, (hostfileLastLineIndex+1)*sizeof(int)); // initialize
                                                                    // nodeList with 0's	      
	for(int i = 0; i < procsNeeded; ++i){
	   tempLineIndex = failedRanks[i]/SLOTS; // started from 0
	   nodeList[tempLineIndex] = 1;
	}
	for(int nodeCounter = 0; nodeCounter < (hostfileLastLineIndex+1); ++nodeCounter)
	   totNodeFailed += nodeList[nodeCounter];
        *numNodeFails = totNodeFailed;
	// Check if there is sufficient spare node available for recovery
        if((hostfileLastLineIndex - totNodeFailed - sumPrevNumNodeFails) < 
           (oldGroupSize-1)/SLOTS) {
           if(rank == 0)
              printf("[%d] There is no sufficient spare node available for recovery.\n",
                     rank);
           exit(0);
        }

	failedNodeList = (int *) malloc(totNodeFailed * sizeof(int));
	memset(failedNodeList, 0, totNodeFailed * sizeof(int)); // initialize
                                                                // failedNodeList with 0's
	failedNodeCounter = 0;
	for(int nodeCounter = 0; nodeCounter < (hostfileLastLineIndex+1); ++nodeCounter){
	   if(nodeList[nodeCounter] == 1)
	      failedNodeList[failedNodeCounter++] = nodeCounter;
	}
#endif //RECOV_ON_SPARE_NODES
#endif //RUN_ON_COMPUTE_NODES
#endif //SHRUNKEN_RECOVERY
        hostNameFailed = NULL;
#pragma omp parallel for default(shared)
	for(i = 0; i < procsNeeded; ++i){
           // Assign list of processes failed
           listFails[i] = failedRanks[i];
#ifdef RUN_ON_COMPUTE_NODES
	   tempLineIndex = failedRanks[i]/SLOTS; // started from 0
#ifdef RECOV_ON_SPARE_NODES
           for(int j = 0; j < totNodeFailed; ++j){
              if(failedNodeList[j] == tempLineIndex)
                 hostfileLineIndex = hostfileLastLineIndex - j - sumPrevNumNodeFails;
           }
#else // Recovery on the same node (no node failure, only process failure)
	   hostfileLineIndex = tempLineIndex;	      
#endif //RECOV_ON_SPARE_NODES
	   hostNameToLaunch[i] = getHostToLaunch(hostfileLineIndex);
           hostNameFailed = getHostToLaunch(tempLineIndex);           
#else // Run on head node or personal machine
           hostNameToLaunch[i] = (char *)hostName;
           hostNameFailed = (char *)hostName;
#endif //RUN_ON_COMPUTE_NODES          
           if(verbosity > 0 && rank == 0)
	      printf("--- Iteration %d: Application: Process %d on node %s is failed! "
                     "---\n", iteration, failedRanks[i], hostNameFailed);
        }
#ifdef RUN_ON_COMPUTE_NODES
        // Release memory of hostNameFailed
        free(hostNameFailed);
#endif
        appToLaunch = (char **) malloc(procsNeeded * sizeof(char *));
        argvToLaunch = (char ***) malloc(procsNeeded * sizeof(char **));
        procsNeededToLaunch = (int *) malloc(procsNeeded * sizeof(int));
        hostInfoToLaunch = (MPI_Info *) malloc(procsNeeded * sizeof(MPI_Info));
        argv[argc] = NULL;
#pragma omp parallel for default(shared)
        for(i = 0; i < procsNeeded; ++i){
            appToLaunch[i] = (char *)argv[0];
            argvToLaunch[i] = (char **)argv;
            procsNeededToLaunch[i] = 1;
            // Host information where to spawn the processes
            MPI_Info_create(&hostInfoToLaunch[i]);
            MPI_Info_set(hostInfoToLaunch[i], (char *)"host", hostNameToLaunch[i]);
        }

#ifdef HANG_ON_REMOVE
        OMPI_Comm_agree(tempShrink, &flag);
#endif
#ifndef SHRUNKEN_RECOVERY
        spawnTime = MPI_Wtime();
	// Spawn the new process(es)
	if(MPI_SUCCESS != (ret = MPI_Comm_spawn_multiple(procsNeeded, appToLaunch,
	   argvToLaunch, procsNeededToLaunch, hostInfoToLaunch, 0, tempShrink,
           &tempIntercomm, MPI_ERRCODES_IGNORE))) {
	   free(tempRanks);
	   free(failedRanks);
	   free(errCodes);
	   if(MPI_ERR_COMM  == ret)
	      printf("Iteration %d: MPI_Comm_spawn_multiple: Invalid communicator "
                     "(parent)\n", iteration);
	   if(MPI_ERR_ARG  == ret)
	      printf("Iteration %d: MPI_Comm_spawn_multiple: Invalid argument "
                     "(parent)\n", iteration);
	   if(MPI_ERR_INFO  == ret)
	       printf("Iteration %d: MPI_Comm_spawn_multiple: Invalid info (parent)\n",
                      iteration);
	   if((MPI_ERR_PROC_FAILED == ret) || (MPI_ERR_REVOKED == ret)){
	      OMPI_Comm_revoke(tempShrink);
              return repairComm(broken, repaired, iteration, listFails, numFails,
                           numNodeFails, sumPrevNumNodeFails, argc, argv, verbosity);
	   }
	   else{
	      fprintf(stderr, "Iteration %d: Unknown error with MPI_Comm_spawn_multiple "
                              "(parent): %d\n", iteration, ret);
	      exit(1);
	   }
	}
	else{
	   if(verbosity > 0 && rank == 0){
	      for(i = 0; i < procsNeeded; ++i)
		 printf("Iteration %d: MPI_Comm_spawn_multiple (parent) [spawning "
                        "failed process %d on node %s]: SUCCESS\n", iteration,
                        failedRanks[i], hostNameToLaunch[i]);
	   }
	   // Memory release. Moving the last two to the end of the function causes
           // segmentation faults for 4 processes failure
	}
	if (verbosity > 0 && rank == 0)
	   printf("MPI_Comm_spawn_multiple takes %0.6f Sec\n", MPI_Wtime() - spawnTime);
	mergeTime = MPI_Wtime();
	// Merge the new processes into a new communicator
	if(MPI_SUCCESS != (ret = MPI_Intercomm_merge(tempIntercomm, false,
           &unorderIntracomm))) {
	   free(tempRanks);
	   free(failedRanks);
	   if((MPI_ERR_PROC_FAILED == ret) || (MPI_ERR_REVOKED == ret)){
	      // Start the recovery over again if there is a failure
	      OMPI_Comm_revoke(tempIntercomm);
              return repairComm(broken, repaired, iteration, listFails, numFails,
                                numNodeFails, sumPrevNumNodeFails, argc, argv, verbosity);
	   }
	   else if(MPI_ERR_COMM == ret){
	      fprintf(stderr, "Iteration %d: Invalid communicator in "
                              "MPI_Intercomm_merge (parent) %d\n", 
                              iteration, ret);
	      exit(1);
	   }
	   else if(MPI_ERR_INTERN == ret){
	      fprintf(stderr, "Iteration %d: Acquaring memory error in "
                              "MPI_Intercomm_merge () %d\n", iteration, ret);
	      exit(1);
	   }
	   else{
	      fprintf(stderr, "Iteration %d: Unknown error with MPI_Intercomm_merge: "
                              "%d\n", iteration, ret);
	      exit(1);
	   }
	}
	else{
	   if(verbosity > 1 )
	      printf("Iteration %d: MPI_Intercomm_merge (parent): SUCCESS\n", iteration);
	}
	if (verbosity > 0 && rank == 0)
	   printf("MPI_Intercomm_merge takes %0.6f Sec\n", MPI_Wtime() - mergeTime);

        agreeTime = MPI_Wtime();
	// Synchronize. sometimes hangs in without this
	// position of code and intercommunicator (not intra) is important
#ifdef HANG_ON_REMOVE
        OMPI_Comm_agree(tempIntercomm, &flag);// to fix hang on problem
#endif
	if (verbosity > 0 && rank == 0)
	   printf("OMPI_Comm_agree takes %0.6f Sec\n", MPI_Wtime() - agreeTime);
	// Sending failed ranks and number of processes failed to the the newly
        // created ranks.
	// oldGroupSize is the size of communicator before failure.
	// procsNeeded is the number of processes that are failed.
	int * child = (int *) malloc(procsNeeded*sizeof(int));
#pragma omp parallel for default(shared)
	for(i = 0; i < procsNeeded; ++i)
	   child[i] = oldGroupSize - procsNeeded + i;
	MPI_Comm_rank(unorderIntracomm, &newRank);
	if(newRank == 0){
	   int send_val[2];
	   for(i = 0; i < procsNeeded; ++i){
	      send_val[0] = failedRanks[i]; send_val[1] = procsNeeded;
	      if(MPI_SUCCESS != (ret = MPI_Send(&send_val, 2, MPI_INT, child[i],
                 MERGE_TAG, unorderIntracomm))) {
		 free(tempRanks);
		 free(failedRanks);
		 if((MPI_ERR_PROC_FAILED == ret) || (MPI_ERR_REVOKED == ret)){
		    // Start the recovery over again if there is a failure
		    OMPI_Comm_revoke(unorderIntracomm);
                    return repairComm(broken, repaired, iteration, listFails, numFails,
                                      numNodeFails, sumPrevNumNodeFails, argc, argv,
                                      verbosity);
		 }
		 else{
		    fprintf(stderr, "Iteration %d: Unknown error with MPI_Send1 "
                                    "(parent): %d\n", iteration, ret);
		    exit(1);
		 }
	      }
	      else{
		 if(verbosity > 1 )
		    printf("Iteration %d: MPI_Send1 (parent): SUCCESS\n", iteration);
	      }
	   }
	}
	// Split the current world (splitted from original) to order the ranks.
	MPI_Comm_rank(unorderIntracomm, &newRank);
	MPI_Comm_size(unorderIntracomm, &nprocs);
	// For one or more process failure (ordering)
	shrinkMergeList = (int *) malloc(nprocs*sizeof(int));

	j = 0;
	for(i = 0; i < nprocs; ++i){
	   if(rankIsNotOnFailedList(i, failedRanks, procsNeeded))
	      shrinkMergeList[j++] = i;
	}
	for(i = j; i < nprocs; ++i)
	   shrinkMergeList[i] = failedRanks[i-j];
	for(i = 0; i < (nprocs - procsNeeded); ++i){
	   if(newRank == i)
	      rankKey = shrinkMergeList[i];
	}
	if(MPI_SUCCESS != (MPI_Comm_split(unorderIntracomm, 0, rankKey, repaired))){
	   if((MPI_ERR_PROC_FAILED == ret) || (MPI_ERR_REVOKED == ret)){
	      // Start the recovery over again if there is a failure
	      OMPI_Comm_revoke(unorderIntracomm);
              return repairComm(broken, repaired, iteration, listFails, numFails,
                                numNodeFails, sumPrevNumNodeFails, argc, argv,
                                verbosity);
	   }
	   else{
	      fprintf(stderr, "Iteration %d: Unknown error with MPI_Comm_split "
                              "(parent): %d\n", iteration, ret);
	      exit(1);
	   }
	}
	else{
	   if(verbosity > 1 )
	      printf("Iteration %d: MPI_Comm_split (parent): SUCCESS\n", iteration);
	}
#endif //SHRUNKEN_RECOVERY

#ifdef SHRUNKEN_RECOVERY
	*repaired = tempShrink;
#endif //SHRUNKEN_RECOVERY

	// Release memory
	free(appToLaunch);
	free(argvToLaunch);
	free(procsNeededToLaunch);
	free(hostInfoToLaunch);
	free(hostNameToLaunch);
#ifndef SHRUNKEN_RECOVERY
	free(shrinkMergeList);
	MPI_Comm_free(&tempShrink);
	free(child);
	MPI_Comm_free(&tempIntercomm);
        MPI_Comm_free(&unorderIntracomm);
#endif //SHRUNKEN_RECOVERY
	free(errCodes);
	free(tempRanks);
	free(failedRanks);
	MPI_Group_free(&failedGroup);
	MPI_Group_free(&oldGroup);
#ifdef RECOV_ON_SPARE_NODES
        if(failedNodeList != NULL)
           free(failedNodeList);
        if(nodeList != NULL)
           free(nodeList);
#endif      
}//repairComm()
#endif

///////////////////////////////////////////////////////////////////////////////////
int rankIsNotOnFailedList(int rank, int * failedList, int numFails){
        int i;

        for(i = 0; i < numFails; ++i){
           if(rank == failedList[i])
              return 0;
        }
        if(i == numFails)
           return 1;
        printf("Error in rankIsNotOnFailedList().\n");
        exit(1);
}//rankIsNotOnFailedList()

///////////////////////////////////////////////////////////////////////////////////
int getHostfileLastLineIndex(void){
     // hostfile line index is started from 0
     int lastLineIndex = 0;
     char *buffer = (char *) malloc(1024);
     FILE *fPointer;

     if(NULL == (fPointer = fopen("hostfile", "r"))){
	 printf("Error in opening file. Exit\n");
	 exit(0);
     }
     fseek(fPointer, 0, SEEK_SET); // make sure start from 0
     while(!feof(fPointer)){
        memset(buffer, 0x00, 1024); // clean buffer
        int ret = fscanf(fPointer, "%[^\n]\n", buffer); // read file
        ++lastLineIndex;
        if(ret == EOF)
           break;
     }
     // Memory deallocation
     fclose(fPointer);
     if(lastLineIndex - 1 == -1) {
        printf("hotfile is empty...\n");
        fflush(NULL);
        exit(0);
     }
     else 
        return (lastLineIndex - 1);
}//getHostfileLastLineIndex()

///////////////////////////////////////////////////////////////////////////////////
char * getHostToLaunch(int hostfileLineIndex){
        FILE *fPointer;
        char lineRead[256];
        char * buffer;
        int curLineIndex = 0;

        if(NULL == (fPointer = fopen("hostfile", "r"))){
            printf("Error in opening file. Exit\n");
            exit(0);
        }
        // Extract a line of string from hostfile "hostfile" (with hostname followed
        // by slots) with line index "hostfileLineIndex" (started from 0)
        while(curLineIndex <= hostfileLineIndex)
        {
            if(NULL == (fgets(lineRead, 256, fPointer))){
                printf("fgets() encountered NULL\nProblem in hostfile\nExit\n");
                exit(0);
            }
            ++curLineIndex;
        }
        buffer = (char *) malloc(strlen(lineRead)+1);
        strcpy(buffer, lineRead);
        // Extract only hostname without slots information
        buffer = strtok (buffer, " ");
        // Memory deallocation
        fclose(fPointer);
        // buffer is deallocated from caller side

        return buffer;
}//getHostToLaunch()

///////////////////////////////////////////////////////////////////////////////////
int getSlots(void) {
   FILE *f;
   char string[1000], seps[] = " \n,( )!=";
   char *p;
   f = fopen("hostfile","r");
   if(!f) {
      printf("Probably executing without a \"hostfile\" file. Returning slots = 1.\n");
      return 1; // return slots = 1
   }
   while(fgets(string, sizeof(string)-1, f) != NULL) {
      // Break into tokens
      p = string;
      // Find first token
      p = strtok(string, seps); 
      while(p != NULL){
         //printf("Token: %s\n", p);
         p = strtok(NULL, seps); // find next token

         if (strncmp("slots", p, strlen("slots")) == 0) {
            // "slots" token is found
            // our "desired token" (a number) is after "="
            p = strtok(NULL, seps); 
            fclose(f); // close the opened file
            return atoi(p); // convert string into int and return
         }
      }
   }
   fclose(f); // close the opened file
   return 0;  // some compilers give warning without this
} //getSlots
