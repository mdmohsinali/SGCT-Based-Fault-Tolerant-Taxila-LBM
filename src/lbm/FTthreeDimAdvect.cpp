/* SGCT-Based-Fault-Tolerant-Taxila-LBM Code source file.
   Copyright (c) 2015, Md Mohsin Ali. All rights reserved.
   Licensed under the terms of the BSD License as described in the LICENSE_FT_CODE file.
   This comment must be retained in any redistributions of this source file.
*/

// Parallel 3D Advection program
// written by Peter Strazdins, Jun 14
// updated for FT-LBM by Mohsin Ali, November 14

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> //getopt(), gethostname()
#include <sys/stat.h>
#include <cstdlib>
#include <string.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include <assert.h>
#include <signal.h>
#ifndef NON_FT
#include "mpi-ext.h"
#endif
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef TAU_PROF
#include <Profile/Profiler.h>
#endif

#include "FTAdvect3D.h"
#include "LinGrid3D.h"
#include "FTGridCombine3D.h"
#include "Vec3D.h"
#include "FaultSimulator.h"
#include "FailureRecovery.h"
#include "LbmInterface.h"
#include "FTSparseGrid3D.h"

//#define USE_TIMER_CLASS             // use Timer class for time measurement
//#define ON_REAL_ENVIRONMENT         // defining this causes real process failure by killing selected 
                                    // process(es) manipulated by `ps' or `top'.
                                    // otherwise, simulated real or non-real.
#define CONSIST_TAG           10000
#define PROCS_LIST_FAILS_TAG  10001
#define B_TAG                 10002

#define USAGE   "\n \
        ./run3dAdvectRaijin [-d] [-2] [-h] [-i iFail] [-j jFail] [-k kFail] [-x xFail] \n \
        \t[-y yFail] [-z zFail] [-v v] [-V vT] [-m method] [-l level] \n \
        \t[-p pD0] [-q pD1] [-r pD2] [-P P] [-s sgP] [-T tF] [-S steps] \n \
        \t[-C nC] [-f] [-U] [-F][-n] [g.x [g.y [g.z]]] \n\t "

#define CONSTRAINTS "\n \
        l <= g.x, g.y, g.z; \n \
        2^sgP <= nprocs; \n \
        sgP<=g.x+g.y+g.z-dim*(l-1), where dim = 2 or 3 \n \
        pD0 => degP; \n \t"

#define DEFAULTS "\n \
        level = 4; \n \
        v = vT = 0; \n \
        iFail, jFail, kFail, xFail, yFail, zFail = -1; \n \
        m = 2 (2nd order); \n \
        g.x = g.y = g.z = 6; \n \
        degP = 1;\n \
        2^sgP ~= nprocsUg; \n \
        pD0 = degP; \n \
        pD1 = (pD0+1)/2; \n \
        pD2 = (pD1+1)/2; \n \
        pD3 = (pD2+1)/2 (on extra layer); \n \
        P = 1; \n \
        tF = 0.25; \n \
        nC = 1; \n \t"

#define NOTES "\n \
        -h: for help; \n \
        -2 for 2D (forces gz = 0) \n \
        -f: fixes dt across grids; \n \
        -U: read full grid saved data from disk; \n \
	\t(`U' option will be applicable for the repeated runs if uses \n \
	\tthe same number of MPI processes with the same degree of \n \
	\tparallelism on each run); \n \
        -F: full combined grid (must set for visualization); \n \
        -d: debug grid combination; \n \
        v: verbosity; \n \
        iFail, jFail, kFail, xFail, yFail, zFail: mpi rank which is simulate \n \
        \tto be failed; \n \
        vT: verbosity of timer output; \n \
        level = n: sparse grid combination for level n; \n \
        pD0 = n: n processes for each grid on first layer (upper); \n \
        pD1 = n: n processes for each grid on second layer; \n \
        pD2 = n: n processes for each grid on third layer; \n \
        pD3 = n: n processes for each grid on forth layer \n \
        \t(lower, that is, single extra layer); \n \
        P = n: set n OpenMP threads per process; \n \
        -n: print required nprocs only; \n \
        tF = n: final time n seconds; \n \
        nC = n: total n combinations; \n \
        steps = n: total n steps; \n \
        sgP = n: n processes on sparse grid; \n \
        g.x, g.y, g.z: grid size (g.x) X (g.y) X (g.z); \n \t"

#define OPTCHARS "i:j:k:x:y:z:v:m:l:s:dp:q:r:P:T:S:C:V:hfUFn2"

#define NOT_SET                        -1
#define NO_FAILURE_ON_RERUNNING_LBM    0
#define INDEX_ZERO                     0
#define BEFORE_FAILURE                 true
#define AFTER_FAILURE                  false
#define SETTING_PARAMETER_REQUIRED     true  // creating new LBM parameters files 
#define SETTING_PARAMETER_NOT_REQUIRED false // not creating new LBM parameters files
#define COMPONENT_GRID                 true  // this is component grid
#define FULL_GRID                      false // this is full grid
#define REPEATED_COMBI_FLAG_SET        true  // flag identifying repeated combination is set
#define REPEATED_COMBI_FLAG_NOT_SET    false // flag identifying repeated combination is not set
#define COMP_GRID_STRING_NOT_PASSED    NULL
#define DOUBLE_ARRAY_MOT_PASSED        NULL
#define REAL_CALL                      true
#define DUMMY_CALL                     false

/*  0: top-level messages from rank 0, 1: 1-off messages from all ranks,
    2: per-iteration messages, 3: dump data: one-off, 4: dump data, per itn.
 */
static bool haveExtraGrids = true;// specifying that application is fault-tolerant     
static int verbosity = 0;        // v above    
static bool isHelp = false;
static int verbTim  = 0;         // vT above. Verbosity (lint maxLevelevel) of timer output 
static bool is2D = false;        // if set, perform 2D simulation (force g.z==0)
static bool isDbgGridComb = false; // set to test/debug combination alg
static int method = NOT_SET;  
static int level = 4;            // use combination algorithm for l>1 
static Vec3D<int> grid;
static int sgProcs = NOT_SET /*unset*/; //2^sgProcs used for sparse grid collection
static double tF = 0.25,         // final time; must it be a multiple of 0.25?
	      CFL = 0.25;        // CFL condition number 
static int rank,                 // global MPI rank
           nprocs;               // global MPI process count
static int nprocsUg;             // process count for sparse grid combination technique
static int pD[4];                // number of processes for outer/inner/... planes
static int P = 1;                // number of OpenMP threads per process
static int steps = NOT_SET;      // #timesteps for sparse grid (dept. on grid & TF)
static int nC = 1;               // apply combination algorithm nC times over TF
static bool isFixedDt = true;    // if set, use a fixed dt across grids
                                 // default is true, otherwise too hard for 3D (at the momemt)
static bool isReadFromDisk = false;// if set, use full grid data saved on disk
static int totFails = 0;         // a total of totFails number of processes (non-real or real) is failed
static bool isFixedProcs;        // signifies if #processes isfixed across diagonal
                                 // Otherwise pD refers to the central points and 
                                 // the number doubles as you go outwards.
                                 // If set, only gets good load balance if isFixedDt
static int s;                    // counter of repeated combination

static int iFail = -1, jFail = -1, kFail = -1, xFail = -1, yFail = -1, zFail = -1;
static int simFailCounter = 0;

// Global variable declarations and initializations
static MPI_Errhandler newEh;     // error handler object
static Timer * timer = NULL;     // Timer object
static int grank;                // rank on a local component grid communicator
static int sumPrevNumNodeFails = 0,// count of sum of previous nodes fail 
           numNodeFails = 0;     // count of node fails
static MPI_Status procsListFailsStatus;// status object
static char ** argvChild = NULL; // argument of child
static int keyVal = 0,           // id of a splitted communicator
           keyCol = 0;           // local rank of a splitted communicator
static int isNewlySpawned = 0;   // determine whether newly spawned or not ('1' for newly spawned, '0' for not)
static int isParent = 1;         // used to determine parent ('1' for parent, '0' for child)
static FTGridCombine3D * gc = NULL; // an object of grid combination

static MPI_Comm myGridComm;      // communicator for a component grid
static MPI_Comm myFullGridComm;  // communicator for the full grid
static MPI_Comm myCommWorld = MPI_COMM_WORLD;// the global communicator
static MPI_Comm combinationComm; // communicator for each combination
static int randCombCount = 0;    // hold random number to determine when process failure will happen
static bool useFullGrid = false; // if 1, use a dense grid data struct for s.g. (must set for visualization)

static ProcGrid3D *g = NULL,     // process grid configuration of a component grid
                  *sg = NULL;    // process grid configuration of the sparse grid
static Vec3D<int> myGrid;        // grid configuration of a component grid
static HaloArray3D **u;        // contains field of a component grid
static FTSparseGrid3D **usg;      // contains field of the sparse grid for a single species
static FTHaloArray3D *usgFullGrid = NULL;// contains field of the full grid
static FTSparseGrid3D **usgPerGrid;// contains field when only a single component grid is used in 
                                 // combination for a single species
static FTSparseGrid3D *uvsg = NULL;// ontains field of the sparse grid for all the species'
static FTSparseGrid3D *uvsgPerGrid = NULL;// contains field when only a single component grid is used in 
                                 // combination for all the species'
static FTAdvect3D *adv = NULL;   // an object for advection

static int numFails = 0,         // count of the application processes failed
           *listFails = NULL;    // an array of processes that failed
static int totRealFails = 0;     // a copy of count of the application processes failed

static int combCommSize;         // size of a component grid communicator where multiple degree
                                 // of parallelism may occur             
static int fullCommPartitionSize;// size of a full grid communicator where multiple degree
                                 // of parallelism may occur    
static int degP = 1;             // degree of parallelization
static double * LbmFi = NULL;    // array (1D double) to cache LBM's field for component grid
static double * LbmFiSG = NULL;  // array (1D double) to cache LBM's sparse grid field
static double * LbmFiFull = NULL;// array (1D double) to cache LBM's field for full grid

static int xs, ys, zs,           // starting grid indices for component grid, like c (no ghost points)      
	   xm, ym, zm,           // widths of local grid for component grid, like c (no ghost points)
	   nDim = 3;             // grid dimension for component grid
static int xsFull, ysFull, zsFull,// starting grid indices for full grid, like c (no ghost points)      
           xmFull, ymFull, zmFull,// widths of local grid for full grid, like c (no ghost points)
           nDimFull = 3;         // grid dimension for full grid
static int LbmSize = 1;          // count of double type elements in LBM's component grid
static int LbmSizeFull = 1;      // count of double type elements in LBM's full grid
static int nSpecies = 1;         // number of species of application
static int B = 1;                // count of blocks of doubles per logical element 
                                 // (default is a single double)
MPI_Fint fortranCommWorld;                                 

static double fullGridStartTime = 0.0, 
              fullGridTotalTime =0.0;
static double gatherScatterStartTime = 0.0, 
              gatherScatterTotalTime = 0.0;
static double simulateLBMStartTime = 0.0, 
              simulateLBMTotalTime = 0.0;

// Function prototypes
void usage(std::string msg);
void getArgs(int argc, char *argv[], int rank, int nprocs);
int ranksOnFailedGrid(int rank, int numFails, int * listFails, FTGridCombine3D * gc);
int gridIsFailedGrid(int gridId, int numFails, int * listFails, FTGridCombine3D * gc);
static void printLocGlobAvgsRelatives(std::string name, double total, double fullNorm1, 
                   int nlVals, Vec3D<int> gix, MPI_Comm comm, int gridId = NOT_SET);
void createGridsBuildCommsConfigProcsGrids(void);
void allocateMemory(int s = NOT_SET, int nFails = NO_FAILURE_ON_RERUNNING_LBM);
void attemptFailureRecovery(MPI_Comm parent, int argc, char** argv, int s);
void simulateNonRealFailure(void);
int sumOfPrevGridsProcs(int gid, int level, int pd[]);
void combineFieldsOfSeveralSpecies(void);
void createConsistentValuesAfterFailure(void);
void callFullGridLBM(void);
void printErrorsComponentFields(void);
void printErrorsCombinedFields(void);
void printElapsedTimes(void);
void visualizeSparseGrid(void);
#ifdef __cplusplus
extern "C"{
#endif
void callLBMSubroutines(ProcGrid3D *g, int gridSizeX, int gridSizeY, int gridSizeZ, 
                   MPI_Comm myCComm, char* subFullIden, bool isSetParamRequired, 
                   double* u, int xs, int ys, int zs, int xm, int ym, int zm, int nDim, 
                   bool isComponentGrid, bool isRepeatComb, int nDoubles, bool isRealCall);
void runLBM(bool isSetParamRequired);
void reRunLBM(bool isRepeatedComb);
#ifdef __cplusplus
}
#endif

#define POWER_OF_TWO(x) (((x)&(x-1)) == 0)
#define IS_POWER2(n) (((n) & (n-1)) == 0)

/****************************************************************************************
 ****************************************************************************************
 **********                        START of Main Function                      **********
 ****************************************************************************************
 ****************************************************************************************/
int main(int argc, char** argv) {

#ifdef TAU_PROF
        TAU_PROFILE("int main(int, char **)", " ", TAU_DEFAULT);
        TAU_INIT(&argc, &argv); 
#ifndef TAU_MPI
        TAU_PROFILE_SET_NODE(0);
#endif // TAU_MPI
#endif	

        MPI_Init(&argc, &argv);

#ifndef NON_FT
	// Error handler
	MPI_Comm_create_errhandler(mpiErrorHandler, &newEh);
	MPI_Comm_set_errhandler(myCommWorld, newEh);
#endif

	// Getting parent
	MPI_Comm parent;
	MPI_Comm_get_parent(&parent);

	// Control starts from here
	if (MPI_COMM_NULL != parent) { // this is child (newly spawned)
            isNewlySpawned = 1; // this is newly spawned
            totFails = 0;
            isDbgGridComb = false;
	} // end of child (newly spawned)

	if (MPI_COMM_NULL == parent) { // this is parent
            isNewlySpawned = 0; // this is NOT newly spawned
            totFails = 0;

            MPI_Comm_rank(myCommWorld, &rank);
            MPI_Comm_size(myCommWorld, &nprocs);
            // Getting arguments
            getArgs(argc, argv, rank, nprocs);
            // Building communicators, creating grids, configuring process grids
            createGridsBuildCommsConfigProcsGrids();
            // Running LBM
            runLBM(SETTING_PARAMETER_REQUIRED);
            // Allocating necessary memory
            allocateMemory(NOT_SET, NO_FAILURE_ON_RERUNNING_LBM);
	    // Only for sub-grids
	    if (not isDbgGridComb) {
               for (int hi = 0; hi < nSpecies; hi++) {
	          // Initialize grids with LBM provided results
                  adv->initGridLBM(useFullGrid, uvsg, u[hi], LbmFi, LbmSize, myCommWorld, 
                                    nSpecies, hi, COMPONENT_GRID, gc->is2D());
                  // Update boundary
                  adv->updateBoundary(u[hi], g);
               }
	    }
	    if (verbosity > 2) {
               for (int hi = 0; hi < nSpecies; hi++) 
   	          u[hi]->print(rank, "initial field");
            }
	} // end of parent

        // Iterations for performing a number of combinations
	for (s = 0; s < nC; s++) {
            MPI_Comm_rank(myCommWorld, &rank);
            MPI_Comm_size(myCommWorld, &nprocs);
            isParent = 1;

            if (s > 0) { // multiple combinations are enabled
               // Rerun LBM with scatterred fields (after combination)
               reRunLBM(REPEATED_COMBI_FLAG_SET);
               // Initialize grid with re-runned LBM results and update the boundary 
               for (int hi = 0; hi < nSpecies; hi++) {
                  // Initialize grids with LBM provided results
                  adv->initGridLBM(useFullGrid, uvsg, u[hi], LbmFi, LbmSize, myCommWorld, 
                                    nSpecies, hi, COMPONENT_GRID, gc->is2D());
                  // Update boundary
                  adv->updateBoundary(u[hi], g);
               }
            }

#ifdef ON_REAL_ENVIRONMENT
            // Attempting process failure detection, failed processes identification, and
            // reconstructing the broken communicator
            attemptFailureRecovery(parent, argc, argv, s);
#else
            if (s == 0) {
               if (rank == 0) {
                  srand(time(NULL));
                  randCombCount = rand()%nC; // randCombCount is between 0 to (nC-1)
                  if (randCombCount == 2) randCombCount = 1; // we are going to do only 1, 2, or 4 combinations
                                                             // (randCombCount = 2 means 3 combinations, we are skipping it)
                  printf("***** Failure before %d%s combination ******\n", randCombCount+1, (randCombCount+1 == 1)? "st": 
                        (randCombCount+1 == 2)? "nd": "th");
               }
               // Broadcast this count to know everyone else about this
               MPI_Bcast(&randCombCount, 1, MPI_INT, 0, myCommWorld);
            }

            if (s == randCombCount) // failure will simulate on before any of the combination 
                                    // selected randomly
               simulateNonRealFailure();
#endif

            // Calculating alternate combination coefficients before combination
            gc->alternateCombCoeffs(listFails, numFails);
            if(s == 0 && verbosity > 0) {
               if (not isReadFromDisk) {
                  // Call LBM on full grid and initialize usgFull with LBM results
                  callFullGridLBM();
                  if (sg->myrank >= 0) { // this process holds part of the sparse grid
                     usgFullGrid->zero();
                     adv->initGridLBM(useFullGrid, uvsg, usgFullGrid, LbmFiFull, LbmSizeFull, myCommWorld, 
                                       nSpecies, INDEX_ZERO, FULL_GRID, gc->is2D());
                  }
               }
            }
#if 0
            if (verbosity > 0)
               // Print errors of fields before combination
               printErrorsComponentFields();            
#endif
            // Sunchronize everything before combination
            MPI_Barrier(myCommWorld);
            // Perform combination
            gatherScatterStartTime = MPI_Wtime();
            MPI_Pcontrol(1, "gatherScatter");
            for (int hi = 0; hi < nSpecies; hi++) {
               usg[hi]->zero();
               gc->gatherScatter(u[hi], usg[hi]);
            }
            MPI_Pcontrol(-1, "gatherScatter");
            gatherScatterTotalTime += (MPI_Wtime() - gatherScatterStartTime);

            // Reset numFails and totRealFails
            numFails = 0; totRealFails = 0;
	} // for (s...)

        if (sg->myrank >= 0) // this process holds part of the sparse grid
           // Copy combined fields for several species (if any) into a single array
           combineFieldsOfSeveralSpecies();       

        // Print errors of field after combination
        if(verbosity > 0) 
           printErrorsCombinedFields();

        // Printing elapsed time (time for full grid requires verbosity > 0)
        printElapsedTimes();
#if 1
#ifdef USE_TIMER_CLASS
        // Use Timer class for time measurement
        // Enable this for time measurement and speedup calculation
        timer->dump(myCommWorld, verbTim);
#endif

        // Visualize sparse grid field
        if ((sg->myrank >= 0) && useFullGrid) // this process holds part of the sparse grid
           visualizeSparseGrid();
#endif
	// Memory release
	delete gc;
	delete adv;
        delete timer;
        for (int hi = 0; hi < nSpecies; hi++) {
           delete u[hi];
           delete usg[hi];
        }
        delete[] u;
        delete[] usg;
        delete[] LbmFi;
        delete[] LbmFiFull;
        if (verbosity > 0) {
           for (int hi = 0; hi < nSpecies; hi++) {
              delete usgPerGrid[hi];
           }
           delete[] usgPerGrid;
           delete uvsgPerGrid;
           delete usgFullGrid;
        }
        delete uvsg;
	if (listFails != NULL) 
           free(listFails); // listFails available to both spawned and non-spawned processes
        //sleep(60);
        MPI_Barrier(myCommWorld);

#ifndef NON_FT
        MPI_Errhandler_free(&newEh);
	MPI_Comm_free(&myGridComm);
        MPI_Comm_free(&myFullGridComm);
        MPI_Comm_free(&myCommWorld);
        MPI_Comm_free(&combinationComm);
#endif
	MPI_Finalize();
	return 0;
} // main()

/****************************************************************************************
 ****************************************************************************************
 **********                         END of Main Function                       **********
 ****************************************************************************************
 ****************************************************************************************/

// print a usage message for this program and exit with a status of 1
void usage(std::string msg) {
    if (rank == 0) {
        printf("\nlbm_raijin_cluster: %s\n", msg.c_str());
        printf("\nUSAGE:%s\nCONSTRAINTS:%s\nDEFAULT VALUES:%s\nNOTES:%s",
        USAGE, CONSTRAINTS, DEFAULTS, NOTES);
    }
    exit(1);
}


/////////////////////////////////////////////////////////////////////////////////////////
void getArgs(int argc, char *argv[], int rank, int nprocs) {
	extern char *optarg;    // points to option argument (for -p option)
	extern int optind;      // index of last option parsed by getopt()
	extern int opterr;
	char optchar;           // option character returned my getopt()
	opterr = 0;             // suppress getopt() error message for invalid option
        pD[0] = degP; 
        pD[1] = 0;
        pD[2] = 0;              // for third plane (3D), or first extra layer (FT 2D)
        pD[3] = 0;              // for extra plane (FT 3D), or second extra layer (FT 2D)
	bool isOptNgiven = false; // records  if -n was given

	while ((optchar = getopt(argc, argv, OPTCHARS)) != NOT_SET) {
	   // extract next option from the command line
	   switch (optchar) {
	   case 'v':
	      if (sscanf(optarg, "%d", &verbosity) != 1)   // invalid integer
		  usage("bad value for verbose");
	      break;
           case 'i':
              if (sscanf(optarg, "%d", &iFail) != 1)   // invalid integer
                  usage("bad value for iFail");
              break;
           case 'j':
              if (sscanf(optarg, "%d", &jFail) != 1)   // invalid integer
                  usage("bad value for jFail");
              break;
           case 'k':
              if (sscanf(optarg, "%d", &kFail) != 1)   // invalid integer
                  usage("bad value for kFail");
              break;
           case 'x':
              if (sscanf(optarg, "%d", &xFail) != 1)   // invalid integer
                  usage("bad value for xFail");
              break;
           case 'y':
              if (sscanf(optarg, "%d", &yFail) != 1)   // invalid integer
                  usage("bad value for yFail");
              break;
           case 'z':
              if (sscanf(optarg, "%d", &zFail) != 1)   // invalid integer
                  usage("bad value for zFail");
              break;
	   case 'm':
	      if (sscanf(optarg, "%d", &method) != 1)      // invalid integer
		 usage("bad value for method");
	      break;
	   case 'l':
	      if (sscanf(optarg, "%d", &level) != 1)       // invalid integer
		 usage("bad value for level");
	      break;
	   case 's':
	      if (sscanf(optarg, "%d", &sgProcs) != 1)
		 usage("bad value for sgProcs");
	      break;
	   case 'p':
	      if (sscanf(optarg, "%d", &pD[0]) != 1)
		 usage("bad value for pD0");
	      break;
           case 'q':
	      if (sscanf(optarg, "%d", &pD[1]) != 1)
		 usage("bad value for pD1");
	      break;
           case 'r':
              if (sscanf(optarg, "%d", &pD[2]) != 1) 
                 usage("bad value for pD2");
              break;
	   case 'P':
	      if (sscanf(optarg, "%d", &P) != 1)
		 usage("bad value for P");
	      break;
	   case 'd':
	      isDbgGridComb = true;
	      break;
           case '2':
              is2D = true;
              break;
	   case 'S':
	      if (sscanf(optarg, "%d", &steps) != 1)
		 usage("bad value for steps");
	      break;
	   case 'C':
	      if (sscanf(optarg, "%d", &nC) != 1)
		 usage("bad value for nC");
	      break;
	   case 'T':
	      if (sscanf(optarg, "%lf", &tF) != 1)
		 usage("bad value for tF");
	      break;
	   case 'V':
	      if (sscanf(optarg, "%d", &verbTim) != 1)
		 usage("bad value for vT");
	      break;
	   case 'h':
	      isHelp = true;
	      break;
	   case 'f':
	      isFixedDt = true;
	      break;
	   case 'U':
	      isReadFromDisk = true;
	      break;
	   case 'F':
	      useFullGrid = true; // Must set for visualization
	      break;
	   case 'n':
	      isOptNgiven = true;
	      break;
	   default:
              usage("unknown option");              
	      break;
	   } //switch
	} //while

	grid.x = grid.y = grid.z = 6;

	if (optind < argc)
	   if (sscanf(argv[optind], "%d", &grid.x) != 1)
			usage("bad value g.x");
	if (optind + 1 < argc)
	   if (sscanf(argv[optind + 1], "%d", &grid.y) != 1)
	       usage("bad value g.y");
        if (optind+2 < argc)
           if (sscanf(argv[optind+2], "%d", &grid.z) != 1)
              usage("bad value g.z");
        if (grid.z == 0)
           is2D = true;
        else if (is2D)
           grid.z = 0;

        if (method == NOT_SET)
           method = FTAdvect3D::bestMethod(is2D);

	isFixedProcs = isFixedDt;
        if (pD[0] == 0)
	   pD[0] = degP;                                         // for first plane (3D), or first layer (2D)
	   // pD[0] = FTGridCombine3D::getPD0(nprocs, level, isFixedProcs);
	if (pD[1] == 0)
	   pD[1] = FTGridCombine3D::getPD1(pD[0], isFixedProcs); // for second plane (3D), or second layer (2D)
        if (pD[2] == 0)
           pD[2] = FTGridCombine3D::getPD1(pD[1], isFixedProcs); // for third plane (3D), or first extra layer (FT-2D)
        if (pD[3] == 0)
           pD[3] = FTGridCombine3D::getPD1(pD[2], isFixedProcs); // for extra plane (FT-3D), or second extra layer (FT-2D)

	if (sgProcs == NOT_SET) { // use as many process as possible for the sparse grid
	   sgProcs = (int) log2((double) FTGridCombine3D::nProcs(is2D, level, pD, isFixedProcs, haveExtraGrids));
           int maxSgProcs = grid.x+grid.y+grid.z - (is2D? 2: 3)*(level-1);
           if (sgProcs > maxSgProcs)
              sgProcs = maxSgProcs;
        }

	nprocsUg = FTGridCombine3D::nProcs(is2D, level, pD, isFixedProcs, haveExtraGrids);
	if (isOptNgiven) {
	   printf("%d MPI processes are needed for given command line parameters.\n", 
                FTGridCombine3D::nProcs(is2D, level, pD, isFixedProcs, haveExtraGrids));
	   exit(0);
	}

        int dx = 1 << std::max(grid.x,  std::max(grid.y, grid.z));
        if (steps == NOT_SET)
           steps = (int) (tF / CFL * dx);
        else
           tF = CFL * steps / (1.0 * dx);
	if (isHelp)
	   usage("Help");

#ifdef _OPENMP
	omp_set_num_threads(P);
#endif     
        if (level > std::min(grid.x,  std::min(grid.y, is2D? grid.y: grid.z)))
           usage("level must be <= min(g)");
        if ((1 << sgProcs) > nprocsUg)
           usage("sgProcs too large for required nprocs");
        if (sgProcs > (grid.x+grid.y+grid.z - (is2D? 2: 3)*(level-1)))
           usage("sgProcs too large for given g and level");
        if (nprocsUg != FTGridCombine3D::nProcs(is2D, level, pD, isFixedProcs, haveExtraGrids))
           usage("nprocsUg fails its constraint");
        //if (pD[0] < degP*8)
        //   usage("pD0 should be at least degP*8");

        if (rank == 0) {
           printf("Level %d %dD combination alg to be applied %d times with %d,%d,%d,%d "
                  "procs/grid %s, combine on 2^%d procs (%d MPI procs)\n", 
 	           level, is2D? 2: 3, nC, pD[0], pD[1], pD[2], pD[3],
 	           isFixedProcs? "": "(doubling outwards from center)", 
 	           sgProcs, nprocs); 
#ifdef _OPENMP
	   printf("[%d] ===== %d threads per process =====\n", rank, omp_get_max_threads());

#pragma omp parallel
	   {
	      printf("%d: hello from thread\n", omp_get_thread_num());
	   }
#endif
	}
} //getArgs()


/////////////////////////////////////////////////////////////////////////////////////////
int ranksOnFailedGrid(int rank, int numFails, int * listFails, FTGridCombine3D * gc) {

   int i;
   for(i = 0; i < numFails; i++){
      if(gc->getGid(rank) == gc->getGid(listFails[i]))
         return 1;
   }
   if(i == numFails)
      return 0;
   printf("Error in ranksOnFailedGrid() in FTthreeDimAdvect.cpp file.\n");
   exit(1);      
}//ranksOnFailedGrid()


/////////////////////////////////////////////////////////////////////////////////////////
int gridIsFailedGrid(int gridId, int numFails, int * listFails, FTGridCombine3D * gc) {

   int i;
   for(i = 0; i < numFails; i++){
      if(gridId == gc->getGid(listFails[i]))
         return 1;
   }
   if(i == numFails)
      return 0;
   printf("Error in gridIsFailedGrid() in FTthreeDimAdvect.cpp file.\n");
   exit(1);   
}//gridIsFailedGrid()


/////////////////////////////////////////////////////////////////////////////////////////
#ifdef __cplusplus
extern "C"{
#endif
void c_get_LBM_field(int * n_dims, double * fi_2d, double * fi_3d,
      int * x_start, int * y_start, int * z_start,
      int * x_width, int * y_width, int * z_width,
      bool * isComponentGrid, bool * isRepeatComb) {

      if (*isComponentGrid) { // calling by component grid comms which supposed to be 
                              // called maybe for multiple times                              
         if (not (*isRepeatComb)) {
            nDim = *n_dims;
            if (*n_dims == 2) { // 2D
               xs = *x_start; ys = *y_start;
               xm = *x_width; ym = *y_width;
               LbmSize = xm * ym;
               LbmFi = new double[LbmSize];
               for(int i = 0; i < LbmSize; i++) 
                  LbmFi[i] = fi_2d[i];
            }
            else if (*n_dims == 3) { // 3D
               xs = *x_start; ys = *y_start; zs = *z_start;
               xm = *x_width; ym = *y_width; zm = *z_width;
               LbmSize = xm * ym * zm;
               LbmFi = new double[LbmSize];
               for(int i = 0; i < LbmSize; i++) 
                  LbmFi[i] = fi_3d[i];
            }
            else {
               printf("========= n_dims out-of-range in c_get_LBM_field() \
                       of FTthreeDimAdvect.cpp file for component grid. Abort.\n");
               exit(1);
            }
         }
      }
      else { // calling by full grid comm which supposed to be called only single time
         if (*n_dims == 2) { // 2D
            xsFull = *x_start; ysFull = *y_start;
            xmFull = *x_width; ymFull = *y_width;
            LbmSizeFull = xmFull * ymFull;                   
            LbmFiFull = new double[LbmSizeFull];
            for(int i = 0; i < LbmSizeFull; i++)
               LbmFiFull[i] = fi_2d[i];
         }
         else if (*n_dims == 3) { // 3D
            xsFull = *x_start; ysFull = *y_start; zsFull = *z_start;
            xmFull = *x_width; ymFull = *y_width; zmFull = *z_width;
            LbmSizeFull = xmFull * ymFull * zmFull;
            LbmFiFull = new double[LbmSizeFull];
            for(int i = 0; i < LbmSizeFull; i++)
               LbmFiFull[i] = fi_3d[i];
         }
         else {
            printf("========= n_dims out-of-range in c_get_LBM_field() \
                    of FTthreeDimAdvect.cpp file for full grid. Abort.\n");
            exit(1);
         }
      }      
   }
#ifdef __cplusplus
}
#endif
// c_get_LBM_field()

/////////////////////////////////////////////////////////////////////////////////////////
#ifdef __cplusplus
extern "C"{
#endif
void callLBMSubroutines(ProcGrid3D *g, int gridSizeX, int gridSizeY, int gridSizeZ,
                         MPI_Comm myCComm , char* subFullIden, bool isSetParamRequired,
                         double* u, int xs, int ys, int zs, int xm, int ym, int zm, 
                         int nDim, bool isCompGrid, bool isRepeatComb, int nDoubles, bool isRealCall) {

        // Start wallclock time
        double startClock = MPI_Wtime();

        char parInDir[128] = "param_grid_", gridIdStr[128] = " ";
        int localRank;
        int procGridX, procGridY, procGridZ;
        MPI_Comm_rank(myCComm, &localRank);

        // Concatenate grid_id with parInDir
        if (subFullIden == NULL) { // chose sub-grid parameter dir
           sprintf(gridIdStr,"%d", gc->getGid(rank));
           strcat(parInDir, gridIdStr);
        }
        else // chose full grid parameter dir
           strcpy(parInDir, subFullIden);

        // Set process grid information
        procGridX = g->P.x;
        procGridY = g->P.y;
        procGridZ = g->P.z;

        if (gc->is2D())
           assert(g->P.z == 1);
           
        // Set parameters by local rank 0
        if (localRank == 0 && isRealCall) {
           if (isSetParamRequired)
              setParameters(parInDir, procGridX, procGridY, procGridZ,
                            degP*procGridX*procGridY*procGridZ, gridSizeX, gridSizeY, 
                            gridSizeZ, degP, nC, gc->is2D(), isCompGrid);
        }

        // Synchronize everything of each local communicator in a common point
        MPI_Barrier(myCComm);
 
        // Convert C communicator handle to Fortran handle
        fortranCommWorld = MPI_Comm_c2f(myCComm);

        if(localRank == 0 && isRealCall) {
           printf("\n*******************************************************\n");
           printf("*******             This is TaxilaLBM             *******\n");
           printf("*******       Calling by C++ routine on %s%s with px = %d, py = %d, pz = %d ******\n",
                 (subFullIden == NULL)? "grid ": "", (subFullIden == NULL)? gridIdStr: "full grid", 
                 procGridX, procGridY, procGridZ);
           printf("*******************************************************\n\n");
           fflush(stdout);
        }

        // Only one LBM simulation
        if (isRealCall) {
           c_runLBM(parInDir, &procGridX, &procGridY, &procGridZ, &fortranCommWorld, 
                 &isCompGrid, &isRepeatComb, &nDoubles, &gridSizeX, &gridSizeY, 
                 &gridSizeZ, &nC, &isRealCall);
        }
        else {
           int dummyNProcs;
           MPI_Comm_size(myCComm, &dummyNProcs);
           procGridX = dummyNProcs; 
           procGridY = 1; 
           procGridZ = 1;
           c_runLBM(parInDir, &procGridX, &procGridY, &procGridZ, &fortranCommWorld,
                 &isCompGrid, &isRepeatComb, &nDoubles, &gridSizeX, &gridSizeY,
                 &gridSizeZ, &nC, &isRealCall);
        }
        
        // Stop wallclock and measure elapsed time
        if(localRank == 0 && isRealCall){
           printf("\nTotal wallclock time for LBM (%s%s): %0.2lf sec\n", 
                 (subFullIden == NULL)? "grid ": "", 
                 (subFullIden == NULL)? gridIdStr: "full grid", MPI_Wtime() - startClock);
           fflush(stdout);
        }

        // Both B and nSpecies are 1 for LBM, at the moment
        if (isCompGrid) {
           if (gc->is2D()) {
              B = 1;
              nSpecies = 1;
           }
           else {              
              B = 1;
              nSpecies = 1;
           }
        }
} //callLBMSubroutines()
#ifdef __cplusplus
}
#endif

/////////////////////////////////////////////////////////////////////////////////////////
void createGridsBuildCommsConfigProcsGrids(void) {
	timer = new Timer;
        combCommSize = nprocs/degP;
        for (int dp = 0; dp < degP; dp++) {
           if (rank >= (dp*combCommSize) && rank <= ((dp+1)*combCommSize - 1)) {
              keyVal = dp; 
              keyCol = rank - (dp*combCommSize);
           }
        }
        if (!MPI_SUCCESS == MPI_Comm_split(myCommWorld, keyVal, keyCol,
                &combinationComm)) {
           printf("Error in MPI_Comm_split in createGridsBuildCommsConfigProcsGrids() "
                  "in FTthreeDimAdvect.cpp file. Abort\n");
           exit(1);
        }

        int myPd[4];
        myPd[0] = pD[0]/degP; myPd[1] = pD[1]/degP; myPd[2] = pD[2]/degP; myPd[3] = pD[3]/degP;
        int mySgProcs = sgProcs - (int) log2((double) degP);
        gc = new FTGridCombine3D(level, myPd, isFixedProcs, mySgProcs, is2D, grid, timer,
                 verbosity, isDbgGridComb, degP, combinationComm, haveExtraGrids);

        // Caches grid id corresponding to 3D (or 2D) coordinates, and vice-versa       
        gc->setGridCoordinatesGridIds();

        // Creating communicators for each grid (3D or 2D)
        int * procGridSize = new int [gc->nGrids()];
        for (int ng = 0; ng < (gc->nGrids()); ng++)
           procGridSize[ng] = NOT_SET;        
        if (gc->is2D())
           assert((gc->pgU)->P.z == 1);
        procGridSize[gc->getGid(rank)] = (gc->pgU)->P.x*(gc->pgU)->P.y*(gc->pgU)->P.z;
        combCommSize = nprocs/degP;

        // 2D:
        //   Use v and mu for combination, x, y, and z as block, and s as species
        //   Non-SGCT dimension used for parallelization: z
        // 3D:
        //   Use z, v, and mu for combination, x and y as block, and s as species
        //   Non-SGCT dimension used for parallelization: y
        for (int ng = 0; ng < (gc->nGrids()); ng++) {
           for (int dp = 0; dp < degP; dp++) {
              if (procGridSize[ng] != NOT_SET && rank >= sumOfPrevGridsProcs(ng, level, pD) + 
                  dp*combCommSize && rank < sumOfPrevGridsProcs(ng, level, pD) + 
                  dp*combCommSize + procGridSize[ng]) {
                 keyVal = ng;
                 keyCol = (rank - sumOfPrevGridsProcs(ng, level, pD) -
                          dp*combCommSize)*degP + dp;
              }
           }
        }

        if (!MPI_SUCCESS == MPI_Comm_split(myCommWorld, keyVal, keyCol,
                   &myGridComm)) {
           printf("Error in MPI_Comm_split in createGridsBuildCommsConfigProcsGrids() "
                  "in FTthreeDimAdvect.cpp file. Abort\n");
           exit(1);
        }
        delete[] procGridSize;

	MPI_Comm_rank(myGridComm, &grank);
} //createGridsBuildCommsConfigProcsGrids()


//////////////////////////////////////////////////////////////////////////////////////
#ifdef __cplusplus
extern "C"{
#endif
void runLBM(bool isSetParamRequired) {
#ifdef TAU_PROF
        TAU_PROFILE("void runLBM(bool)", " ", TAU_USER);
#endif

        int nDoubles = LbmSize;
        g = gc->pgU;
	sg = gc->pgS;
	myGrid = gc->gxU;

        if (rank == 0)
           int ret = system("rm -r param_grid_*"); // remove param_grid_* directories
        MPI_Barrier(myCommWorld); // Synchronize everything in this point

        timer->start("simulateLBM", 0, 0);
        MPI_Pcontrol(1, "simulateLBM");
        simulateLBMStartTime = MPI_Wtime();
        // LBM subroutines are called here (for sub-grids)
        callLBMSubroutines(g, (1 << myGrid.x), (1 << myGrid.y), (1 << myGrid.z), 
                            myGridComm, COMP_GRID_STRING_NOT_PASSED, isSetParamRequired, 
                            DOUBLE_ARRAY_MOT_PASSED, xs, ys, zs, xm, ym, zm, nDim, COMPONENT_GRID, 
                            REPEATED_COMBI_FLAG_NOT_SET, nDoubles, REAL_CALL);
        simulateLBMTotalTime += (MPI_Wtime() - simulateLBMStartTime);
        MPI_Barrier(myCommWorld); // Synchronize everything in this point
        MPI_Pcontrol(-1, "simulateLBM");
        timer->stop("simulateLBM");
} //runLBM()
#ifdef __cplusplus
}
#endif

/////////////////////////////////////////////////////////////////////////////////////////
#ifdef __cplusplus
extern "C"{
#endif
void reRunLBM(bool isRepeatedComb) {
        // Rerun LBM with scatterred fields which are getting after gatherScatter of combination
        int nDoubles = LbmSize;
        int reRunIndex = 0, lx, ly, lz;

        for (int hi = 0; hi < nSpecies; hi++) {
           lz = u[hi]->l.z;
           ly = u[hi]->l.y;
           lx = u[hi]->l.x;
           if (gc->is2D())
              assert(lz == 1);               
           else
              lz = (lz % 2 == 0)? lz: lz - 1;               
           ly = (ly % 2 == 0)? ly: ly - 1;
           lx = ((lx/B) % 2 == 0)? lx: lx - B; 
           for (int k = 0; k < lz; k++) {
              for (int j = 0; j < ly; j++) {
                 for (int i = 0; i < lx; i += 2) {
                    LbmFi[reRunIndex] = Vh(u[hi], i, j, k);
                    reRunIndex++;
                 }
              }
           }
        }

        timer->start("simulateLBM", 0, 0);
        MPI_Pcontrol(1, "simulateLBM");
        simulateLBMStartTime = MPI_Wtime();
        // LBM subroutines are called here (for sub-grids)
        callLBMSubroutines(g, (1 << myGrid.x), (1 << myGrid.y), (1 << myGrid.z), 
                            myGridComm, COMP_GRID_STRING_NOT_PASSED, SETTING_PARAMETER_NOT_REQUIRED, LbmFi, 
                            xs, ys, zs, xm, ym, zm, nDim, COMPONENT_GRID, isRepeatedComb, nDoubles, REAL_CALL);
        simulateLBMTotalTime += (MPI_Wtime() - simulateLBMStartTime);
        MPI_Pcontrol(-1, "simulateLBM");
        timer->stop("simulateLBM");
} //reRunLBM()
#ifdef __cplusplus
}
#endif

/////////////////////////////////////////////////////////////////////////////////////////
void allocateMemory(int s, int nFails) {
        // Process failures during re-running the LBM
        // LbmFi memory is not allocated in c_get_LBM_field() for this scenario
        int nDoubles = LbmSize;
        if (s != NOT_SET && nFails > 0)
           LbmFi = new double[nDoubles];
     
        u = new HaloArray3D* [nSpecies];
        usg = new FTSparseGrid3D* [nSpecies];
        for (int hi = 0; hi < nSpecies; hi++) {
	   u[hi] = new HaloArray3D(g->G2L(gc->gridSz_(is2D, myGrid)), Vec3D<int>(1, 1, !gc->is2D()), B);
	   usg[hi] = new FTSparseGrid3D(myCommWorld, sg->G2L(gc->gridSz_(is2D, grid)), 1, useFullGrid, level, is2D, grid, sg, B);
        }

        if (verbosity > 0) {
           if (gc->is2D())
              usgFullGrid = new FTHaloArray3D(myCommWorld, sg->G2L(gc->gridSz_(is2D, grid)*
                             Vec3D<int>(1, nSpecies, 1)-Vec3D<int>(0, nSpecies - 1, 0)), Vec3D<int>(0), B);           
           else
              usgFullGrid = new FTHaloArray3D(myCommWorld, sg->G2L(gc->gridSz_(is2D, grid)*
                             Vec3D<int>(1, 1, nSpecies)-Vec3D<int>(0, 0, nSpecies - 1)), Vec3D<int>(0), B);           

           usgPerGrid = new FTSparseGrid3D* [nSpecies];
           for (int hi = 0; hi < nSpecies; hi++) 
              usgPerGrid[hi] = new FTSparseGrid3D(myCommWorld, sg->G2L(gc->gridSz_(is2D, grid)), 1, useFullGrid, level, is2D, grid, sg, B);

           if (gc->is2D())
              uvsgPerGrid = new FTSparseGrid3D(myCommWorld, sg->G2L(gc->gridSz_(is2D, grid)*
                                Vec3D<int>(1, nSpecies, 1)-Vec3D<int>(0, nSpecies - 1, 0)), nSpecies, useFullGrid, level, is2D, grid, sg, B);
           else
              uvsgPerGrid = new FTSparseGrid3D(myCommWorld, sg->G2L(gc->gridSz_(is2D, grid)*
                                Vec3D<int>(1, 1, nSpecies)-Vec3D<int>(0, 0, nSpecies - 1)), nSpecies, useFullGrid, level, is2D, grid, sg, B);
        }
        if (gc->is2D())
           uvsg = new FTSparseGrid3D(myCommWorld, sg->G2L(gc->gridSz_(is2D, grid)*
                      Vec3D<int>(1, nSpecies, 1)-Vec3D<int>(0, nSpecies - 1, 0)), nSpecies, useFullGrid, level, is2D, grid, sg, B);
        else
           uvsg = new FTSparseGrid3D(myCommWorld, sg->G2L(gc->gridSz_(is2D, grid)*
                      Vec3D<int>(1, 1, nSpecies)-Vec3D<int>(0, 0, nSpecies - 1)), nSpecies, useFullGrid, level, is2D, grid, sg, B);
	adv = new FTAdvect3D(gc->gridSz_(is2D, myGrid), Vec3D<double>(1.0, 1.0, 1.0), tF, CFL,
			method, verbosity, timer, myCommWorld, B);
	if (verbosity > 0) {
	   char hostName[128];
	   gethostname(hostName, sizeof(hostName));
	   // only for sub-grids
	   printf("%d: process (%d,%d,%d) on grid (%d,%d,%d) with %dx%dx%d processes",
		     rank, g->id.x, g->id.y, g->id.z, myGrid.x, myGrid.y, myGrid.z, 
                     g->P.x, g->P.y, g->P.z);
           for (int hi = 0; hi < nSpecies; hi++) 
	      printf(" computes %dx%dx%d points from (%d,%d,%d) on host %s\n", u[hi]->l.x,
		     u[hi]->l.y, u[hi]->l.z, g->L2G0(0, adv->gridSize.x)*B,
		     g->L2G0(1, adv->gridSize.y), g->L2G0(2, adv->gridSize.z), hostName);
	}

	if (rank == 0) 
	   printf("Compute %ssparse grid (%d,%d,%d) with %dx%dx%d procs\n",
		"in-sub-grids ", grid.x, grid.y, grid.z,
		sg->P.x, sg->P.y, sg->P.z);
/*
        for (int hi = 0; hi < nSpecies; hi++) {
           u[hi]->zero(); // Initialize u
           usg[hi]->zero(); // in case gatherRecv() is not called to zero it
        }
        uvsg->zero();
*/
} //allocateMemory()


/////////////////////////////////////////////////////////////////////////////////////////
void attemptFailureRecovery(MPI_Comm parent, int argc, char** argv, int s) {
        // Set value to recognize child
        if (parent != MPI_COMM_NULL) 
            isParent = 0;
      
#ifndef NON_FT                 
        // Call communicator reconstruct function
        if (isNewlySpawned == 1 || isParent == 0) // newly spawned
           // NULL should be replaced by listFails, but totFails is not available on spawned.
           // totFails can be passed as command-line argument and then allocate memory for listFails
           // later
           myCommWorld = communicatorReconstruct(myCommWorld, isNewlySpawned, NULL, &numFails,
                                   &numNodeFails, sumPrevNumNodeFails, argc, argv, verbosity);
#endif
        // Call communicator reconstruct function
        // Make child to parent which become child automatically in communicatorReconstruct
        if (parent != MPI_COMM_NULL) 
           parent = MPI_COMM_NULL;

#ifndef NON_FT
        if (isNewlySpawned != 1) {// not newly spawned
           totRealFails = numProcsFails(myCommWorld);
	   if (totRealFails > 0) {
#ifdef ON_REAL_ENVIRONMENT
	      listFails = (int *) malloc(totRealFails * sizeof(int)); // released later
#endif
              myCommWorld = communicatorReconstruct(myCommWorld, isNewlySpawned, listFails, &numFails,
                                &numNodeFails, sumPrevNumNodeFails, argc, argv, verbosity);
	   }
        }
#endif
        // Bring child to its original position
        if (isParent == 0) 
           parent = myCommWorld;

#ifndef NON_FT
        // Set error handler for new communicator
        MPI_Comm_create_errhandler(mpiErrorHandler, &newEh);
        MPI_Comm_set_errhandler(myCommWorld, newEh);
#endif

        // Make child to parent
        if (MPI_COMM_NULL != parent) {// child
           isParent = 0;
           parent = MPI_COMM_NULL;
        } 
        else  // parent
           isParent = 1;
        isNewlySpawned = 0; // not newly spawned at the moment. maybe spawned before.

        if (numFails > 0) { 
           // numFails is available to all processes, including recovered
	   if (isParent == 0) { // this is child
	      // Passing command-line arguments to re-constructed process(es)
	      MPI_Comm_rank(myCommWorld, &rank);
	      MPI_Comm_size(myCommWorld, &nprocs);
	      argvChild = (char **) malloc((argc - 1) * sizeof(char *)); // released later
	      for (int i = 0; i < argc - 1; i++) 
	         argvChild[i] = strdup(argv[i + 1]);
	      getArgs(argc - 1, argvChild, rank, nprocs);
	   } // end of child

           // Accumulate number of node failures to track the spare node's index in hostfile              
           sumPrevNumNodeFails += numNodeFails;
		
	   // Make B, nSpecies, numFails, listFails, s, and sumPrevNumNodeFails 
           // consistent across all processes by process 0.
	   createConsistentValuesAfterFailure();         
	   createGridsBuildCommsConfigProcsGrids();                            
           if (isParent == 0 || (isParent != 0 && ranksOnFailedGrid(rank, numFails, listFails, gc))) {
              // This is child or child's same-grid process.
              // Initialize re-constructed process and re-initialize
              // processes with the same grid as failed process.
              g = gc->pgU;
	      sg = gc->pgS;
	      myGrid = gc->gxU;

              allocateMemory(s, numFails);

              MPI_Comm_rank(myCommWorld, &rank);
              MPI_Comm_size(myCommWorld, &nprocs);
              printf("[%d] (currently %s): re-initialization, recovered or processes with same gridId\n",
                      rank, (MPI_COMM_NULL == parent) ? "parent" : "child");
           } // end of child or child's same-grid process

           // Free allocated memory
           if (isParent == 0) // this is child
              free(argvChild);
        }// end of "if (numFails > 0)"
} //attemptFailureRecovery()

/////////////////////////////////////////////////////////////////////////////////////////
void simulateNonRealFailure(void) {
            if (iFail != -1)
               simFailCounter++;
            if (jFail != -1)
               simFailCounter++;
            if (kFail != -1)
               simFailCounter++;
            if (xFail != -1)
               simFailCounter++;
            if (yFail != -1)
               simFailCounter++;
            if (zFail != -1)
               simFailCounter++;
            numFails = simFailCounter;

            listFails = (int *) malloc(numFails * sizeof(int)); // released later

            simFailCounter = 0;
            if (iFail != -1)
               listFails[simFailCounter++] = iFail;
            if (jFail != -1)
               listFails[simFailCounter++] = jFail;
            if (kFail != -1)
               listFails[simFailCounter++] = kFail;
            if (xFail != -1)
               listFails[simFailCounter++] = xFail;
            if (yFail != -1)
               listFails[simFailCounter++] = yFail;
            if (zFail != -1)
               listFails[simFailCounter++] = zFail;
} //simulateNonRealFailure()

/////////////////////////////////////////////////////////////////////////////////////////
int sumOfPrevGridsProcs(int gid, int level, int pD[4]) {
        int pd[4];
        pd[0] = pD[0]/degP; 
        pd[1] = pD[1]/degP;
        pd[2] = pD[2]/degP; 
        pd[3] = pD[3]/degP;

        if (gc->is2D()) {
           if (gid <= level) 
              return (gid*pd[0]);
           else if (gid > level && gid <= (2*level - 1)) 
              return (level*pd[0] + (gid-level)*pd[1]);
           else if (gid > (2*level - 1) && gid <= 3*level - 3) 
              return (level*pd[0] + (level-1)*pd[1] + (gid - (2*level-1))*pd[1]/2);
           else if (gid > 3*level - 3) 
              return (level*pd[0] + (level-1)*pd[1] + (level-2)*pd[1]/2 + 
                     (gid - (3*level-3))*pd[1]/4);
        }
        else {
           int nGridsFirstLayer = level*(level+1)/2;
           int nGridsSecondLayer = level*(level-1)/2;
           int nGridsThirdLayer = (level-1)*(level-2)/2;
           int nGridsFirstSecondLayer = level*level;
           int nGridsFirstSecondThirdLayer = nGridsFirstSecondLayer +
                                          (level-2)*(level-1)/2;

           if (gid <= nGridsFirstLayer) 
              return (gid*pd[0]);
           else if (gid > nGridsFirstLayer && gid <= nGridsFirstSecondLayer) 
              return (nGridsFirstLayer*pd[0] + (gid-nGridsFirstLayer)*pd[1]);
           else if (gid > nGridsFirstSecondLayer && gid <= nGridsFirstSecondThirdLayer) 
              return (nGridsFirstLayer*pd[0] + nGridsSecondLayer*pd[1] + 
                     (gid-nGridsFirstSecondLayer)*pd[2]);
           else if (gid > nGridsFirstSecondThirdLayer) 
              return (nGridsFirstLayer*pd[0] + nGridsSecondLayer*pd[1] + 
                     nGridsThirdLayer*pd[2] + (gid-nGridsFirstSecondThirdLayer)*pd[3]);
        }
        return 0; // some compilers give warning without this
} //sumOfPrevGridsProcs()


/////////////////////////////////////////////////////////////////////////////////////////
void createConsistentValuesAfterFailure(void) {
        // Make B, numFails, listFails, s, and sumPrevNumNodeFails 
        // consistent across all processes by process 0.
        int errRet;
        int advConsis[5];
        MPI_Status advConsisStatus;

        MPI_Comm_rank(myCommWorld, &rank);
        MPI_Comm_size(myCommWorld, &nprocs);

        if (rank == 0) {
            advConsis[0] = B;
            advConsis[1] = nSpecies;
            advConsis[2] = numFails;
            advConsis[3] = sumPrevNumNodeFails;
            advConsis[4] = s;
            for (int i = 0; i < numFails; i++) {
                MPI_Send(advConsis, 5, MPI_INT, listFails[i], CONSIST_TAG, myCommWorld);

                // Sending failed ranks information to the failed ranks
                MPI_Send(listFails, numFails, MPI_INT, listFails[i], PROCS_LIST_FAILS_TAG,
                         myCommWorld);
            }
        }

        if (isParent == 0) { // child which was/were failed
            if (MPI_SUCCESS != (errRet = MPI_Recv(advConsis, 5, MPI_INT, 0,
                CONSIST_TAG, myCommWorld, &advConsisStatus))){
                printf("[%d] ERROR: MPI_Recv for advConsis in createConsistentValuesAfterFailure() "
                       "in FTthreeDimAdvect.cpp file.\n", rank);
                exit(1);
            }
            B = advConsis[0];
            nSpecies = advConsis[1];
            numFails = advConsis[2];
            sumPrevNumNodeFails = advConsis[3];
            s = advConsis[4];

            // Receiving failed ranks information to the failed ranks from process 0
            int * tempFails = (int *) malloc(numFails * sizeof(int));
            listFails = (int *) malloc(numFails * sizeof(int));
            if (MPI_SUCCESS != (errRet = MPI_Recv(tempFails, numFails, MPI_INT, 0,
                PROCS_LIST_FAILS_TAG, myCommWorld, &procsListFailsStatus))){
                printf("[%d] ERROR: MPI_Recv for tempFails in createConsistentValuesAfterFailure() "
                       "in FTthreeDimAdvect.cpp file.\n", rank);
                exit(1);
            }

            for (int i = 0; i < numFails; i++) 
                listFails[i] = tempFails[i];
            free(tempFails);
        }
} //createConsistentValuesAfterFailure()


/////////////////////////////////////////////////////////////////////////////////////////
void callFullGridLBM(void) {
        char fullGridStr[128] = "param_grid_full";
        int nDoubles = LbmSizeFull;
        // Creating communicator for full grid (3D or 2D)
        combCommSize = nprocs/degP;
        fullCommPartitionSize = pow(2, sgProcs-(int) log2(degP));

        // 2D:
        //   Use v and mu for combination, x, y, and z as block, and s as species
        //   Non-SGCT dimension used for parallelization: z
        // 3D:
        //   Use z, v, and mu for combination, x and y as block, and s as species
        //   Non-SGCT dimension used for parallelization: y
        //
        // The above information is applicable for GENE, not for LBM
        for (int dp = 0; dp < degP; dp++) {
           if (rank >= (dp*combCommSize) && rank <= (dp*combCommSize + fullCommPartitionSize - 1)) {
              keyVal = 0; 
              keyCol = (rank - (dp*combCommSize))*degP + dp;
           }
           else if (rank >= (dp*combCommSize + fullCommPartitionSize) &&
                 rank <= ((dp+1)*combCommSize - 1)) {
              keyVal = 1; 
              keyCol = rank - (dp*combCommSize + fullCommPartitionSize) +
                       (dp*(combCommSize - fullCommPartitionSize));
           }
        }
        if (!MPI_SUCCESS == MPI_Comm_split(myCommWorld, keyVal, keyCol, 
              &myFullGridComm)) {
           printf("Error in MPI_Comm_split in callFullGridLBM() "
                  "in FTthreeDimAdvect.cpp file. Abort\n");
           exit(1);
        }

        fullGridStartTime = MPI_Wtime();
        // LBM subroutines are called here (for full grid)               
        if (sg->myrank >= 0)
           callLBMSubroutines(sg, (1 << grid.x), (1 << grid.y), (1 << grid.z), 
                               myFullGridComm, fullGridStr, SETTING_PARAMETER_REQUIRED, 
                               DOUBLE_ARRAY_MOT_PASSED, xsFull, ysFull, zsFull, 
			       xmFull, ymFull, zmFull, nDimFull, FULL_GRID, 
                               REPEATED_COMBI_FLAG_NOT_SET, nDoubles, REAL_CALL);        
        else {
           callLBMSubroutines(sg, (1 << grid.x), (1 << grid.y), (1 << grid.z),
                               myFullGridComm, fullGridStr, SETTING_PARAMETER_REQUIRED,
                               DOUBLE_ARRAY_MOT_PASSED, xsFull, ysFull, zsFull,
                               xmFull, ymFull, zmFull, nDimFull, FULL_GRID,
                               REPEATED_COMBI_FLAG_NOT_SET, nDoubles, DUMMY_CALL);
        }
        fullGridTotalTime = (MPI_Wtime() - fullGridStartTime);
} //callFullGridLBM()


/////////////////////////////////////////////////////////////////////////////////////////
void combineFieldsOfSeveralSpecies(void) {
        int zIndex, yIndex, xIndex;
        uvsg->zero();
        if (useFullGrid) { // use full grid for interpolation
           int k2 = 0;
#pragma omp parallel for default(shared)
           for (int hi = 0; hi < nSpecies; ++hi) {
              zIndex = usg[hi]->uh->l.z;
              yIndex = usg[hi]->uh->l.y;
              xIndex = usg[hi]->uh->l.x;
              if (gc->is2D()) {
                 assert(zIndex == 1);                 
                 yIndex = (yIndex % 2 == 1)? yIndex-1: yIndex;
                 for (int k = 0; k < zIndex; ++k) {
                    for (int j = 0; j < yIndex; ++j, ++k2) {
		       for (int i = 0; i < xIndex; ++i) {
		          Vh(uvsg->uh, i, k2, k) = Vh(usg[hi]->uh, i, j, k);
		       }
                    }
                 }
              }
              else {
                 zIndex = (zIndex % 2 == 1)? zIndex-1: zIndex;
                 for (int k = 0; k < zIndex; ++k, ++k2) {
		    for (int j = 0; j < yIndex; ++j) {
		       for (int i = 0; i < xIndex; ++i) {
			  Vh(uvsg->uh, i, j, k2) = Vh(usg[hi]->uh, i, j, k);
		       }
		    }
                 }
              }
           } // end of for (int hi = 0;...)
        }
        else { // use sparse grid for interpolation
           int Nx, iG, k2 = 0, k3;
           double *u1, *v1;
#pragma omp parallel for default(shared)
           for (int hi = 0; hi < nSpecies; ++hi) {
              zIndex = usg[hi]->nz;
              Nx = usg[hi]->gridSz(usg[hi]->g.x);
              xIndex = usg[hi]->pg->G2L(0, Nx);
              iG = usg[hi]->pg->L2G0(0, Nx);
              if (gc->is2D()) {
                 assert(zIndex == 1);
                 for (int k = 0; k < zIndex; ++k) {
                    yIndex = usg[hi]->ny[k];
                    yIndex = (yIndex % 2 == 1)? yIndex-1: yIndex;
                    for (int j = 0; j < yIndex; ++j, ++k2) {
		       u1 = &uvsg->u[uvsg->rx[k][k2]];
		       v1 = &usg[hi]->u[usg[hi]->rx[k][j]];
		       k3 = 0;
		       for (int i = 0; i < xIndex*B; ++i, ++iG) {
			  if (iG % usg[hi]->rs[k][j] == 0 && (k*yIndex*xIndex*B + j*xIndex*B + k3) < usg[hi]->numElts()) {
			     u1[k3] = v1[k3];
			     ++k3;
			  }
		       }
                    }
                 }
              }
              else {
                 zIndex = (zIndex % 2 == 1)? zIndex-1: zIndex;
                 for (int k = 0; k < zIndex; ++k, ++k2) {
		    yIndex = usg[hi]->ny[k];
		    for (int j = 0; j < yIndex; ++j) {
		       u1 = &uvsg->u[uvsg->rx[k2][j]];
		       v1 = &usg[hi]->u[usg[hi]->rx[k][j]];
		       k3 = 0;
		       for (int i = 0; i < xIndex*B; ++i, ++iG) {
			  if (iG % usg[hi]->rs[k][j] == 0 && (k*yIndex*xIndex*B + j*xIndex*B + k3) < usg[hi]->numElts()) {
			     u1[k3] = v1[k3];
			     ++k3;
			  }
		       }
		    }
                 }
              }
           } // end of for (int hi = 0;...)
        }
} //combineFieldsOfSeveralSpecies()


/////////////////////////////////////////////////////////////////////////////////////////
static void printLocGlobAvgsRelatives(std::string name, double total, double fullNorm1, int nlVals,
                Vec3D<int> gix, MPI_Comm comm, int gridId) {
        int grank;
        MPI_Comm_rank(comm, &grank);
        Vec3D<int> fullGridSize;
        if (gc->is2D())
           fullGridSize = gc->gridSz_(is2D, gix) * Vec3D<int>(1, 2, 1) * Vec3D<int>(B, 1, 1);        
        else
           fullGridSize = gc->gridSz_(is2D, gix) * Vec3D<int>(1, 1, 2) * Vec3D<int>(B, 1, 1);        
        int ngVals = fullGridSize.prod();
        char buffer[128] = " ";

        sprintf(buffer, "%d", gridId);

        double v[1];
        double r[1];
        if (total != 0.0) {
           if (verbosity > 0){
              printf("%d: grid id %s: grid (%d,%d,%d): local avg %s is %.2e\n",
                   rank, (gridId != NOT_SET)? buffer: "not provided", gix.x, gix.y,
                   gix.z, name.c_str(), nlVals == 0 ? 0.0 : total / nlVals);

              printf("%d: grid id %s: grid (%d,%d,%d): local relative %s is %.2e, sgError %0.6lf\n",
                   rank, (gridId != NOT_SET)? buffer: "not provided", gix.x, gix.y,
                   gix.z, name.c_str(), fullNorm1 == 0.0 ? 0.0 : total / fullNorm1, total);
           }
        }

        MPI_Reduce(&total, v, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        MPI_Reduce(&fullNorm1, r, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (grank == 0 && v[0] != 0.0){
           printf("\n%d: GRID ID %s: grid (%d,%d,%d): AVERAGE %s %.2e\n",
                   rank, (gridId != NOT_SET)? buffer: "not provided", gix.x, gix.y,
                   gix.z, (char *) name.c_str(), v[0] / ngVals);

           printf("%d: GRID ID %s: grid (%d,%d,%d): RELATIVE %s %.2e, sum %0.6lf, fullNorm1 %0.6lf\n\n",
                   rank, (gridId != NOT_SET)? buffer: "not provided", gix.x, gix.y,
                   gix.z, (char *) name.c_str(), v[0] / r[0], v[0], r[0]);

        }
} //printLocGlobAvgsRelatives()

/////////////////////////////////////////////////////////////////////////////////////////
void printErrorsComponentFields(void) {
        // For calculating error before combination
        int zIndex, yIndex, xIndex;
        double sgError = 0.0, norm1Full = 0.0;
        for (int gId = 0; gId < (gc->nGrids()); ++gId) {
           uvsgPerGrid->zero();
           gc->singleGridCombine(gId);
#pragma omp parallel for default(shared)
           for (int hi = 0; hi < nSpecies; ++hi) {
              usgPerGrid[hi]->zero();
              gc->gather(u[hi], usgPerGrid[hi]);
           }
           gc->singleGridCombine(FTGridCombine3D::GID_ALL);
           // Copy combined fields for several species (if any) into a single array
           if (sg->myrank >= 0) { // this process holds part of the sparse grid              
              if (useFullGrid) { // use full grid for interpolation
                 int k2 = 0;
#pragma omp parallel for default(shared)
                 for (int hi = 0; hi < nSpecies; ++hi) {
                    zIndex = usgPerGrid[hi]->uh->l.z;
                    yIndex = usgPerGrid[hi]->uh->l.y;
                    xIndex = usgPerGrid[hi]->uh->l.x;

                    if (gc->is2D()) {
                       assert(zIndex == 1);
                       yIndex = (yIndex % 2 == 1)? yIndex - 1: yIndex;
                       for (int k = 0; k < zIndex; ++k) {
                          for (int j = 0; j < yIndex; ++j, ++k2) {
			     for (int i = 0; i < xIndex; ++i) {
				Vh(uvsgPerGrid->uh, i, k2, k) = Vh(usgPerGrid[hi]->uh, i, j, k);
			     }
                          }
                       }
                    }
                    else {
                       zIndex = (zIndex % 2 == 1)? zIndex - 1: zIndex;
                       for (int k = 0; k < zIndex; ++k, ++k2) {
			  for (int j = 0; j < yIndex; ++j) {
			     for (int i = 0; i < xIndex; ++i) {
				Vh(uvsgPerGrid->uh, i, j, k2) = Vh(usgPerGrid[hi]->uh, i, j, k);
			     }
			  }
                       }
                    }
                 }
              }
              else { // use sparse grid for interpolation 
                 int Nx, iG, k2 = 0, k3;
                 double *u1, *v1;
#pragma omp parallel for default(shared)
                 for (int hi = 0; hi < nSpecies; ++hi) {
                    zIndex = usgPerGrid[hi]->nz;
                    Nx = usgPerGrid[hi]->gridSz(usgPerGrid[hi]->g.x);
                    xIndex = usgPerGrid[hi]->pg->G2L(0, Nx);
                    iG = usgPerGrid[hi]->pg->L2G0(0, Nx);
                    if (gc->is2D()) {
                       assert(zIndex == 1);
                       for (int k = 0; k < zIndex; ++k) {
                          yIndex = usgPerGrid[hi]->ny[k];
                          yIndex = (yIndex % 2 == 1)? yIndex - 1: yIndex;
                          for (int j = 0; j < yIndex; ++j, ++k2) {
			     u1 = &uvsgPerGrid->u[uvsgPerGrid->rx[k][k2]];
			     v1 = &usgPerGrid[hi]->u[usgPerGrid[hi]->rx[k][j]];
			     k3 = 0;
			     for (int i = 0; i < xIndex*B; ++i, ++iG) {
				if (iG % usgPerGrid[hi]->rs[k][j] == 0 && (k*yIndex*xIndex*B + j*xIndex*B + k3) < usgPerGrid[hi]->numElts()) {
				   u1[k3] = v1[k3];
				   ++k3;
				}
			     }
                          }
                       }
                    }
                    else {
                       zIndex = (zIndex % 2 == 1)? zIndex - 1: zIndex;
                       for (int k = 0; k < zIndex; ++k, ++k2) {
			  for (int j = 0; j < yIndex; ++j) {
			     u1 = &uvsgPerGrid->u[uvsgPerGrid->rx[k2][j]];
			     v1 = &usgPerGrid[hi]->u[usgPerGrid[hi]->rx[k][j]];
			     k3 = 0;
			     for (int i = 0; i < xIndex*B; ++i, ++iG) {
				if (iG % usgPerGrid[hi]->rs[k][j] == 0 && (k*yIndex*xIndex*B + j*xIndex*B + k3) < usgPerGrid[hi]->numElts()) {
				   u1[k3] = v1[k3];
				   ++k3;
				}
			     }
			  }
                       }
                    }
                 }
              }
           }
           if (nC > 0 && !gridIsFailedGrid(gId, numFails, listFails, gc)) {
              sgError = 0.0; norm1Full = 0.0;
              if (sg->myrank >= 0) { // this process holds part of the sparse grid
                 sgError = adv->checkErrorLBM(uvsgPerGrid, usgFullGrid, gc->is2D(), 
                                               isReadFromDisk);
                 norm1Full = usgFullGrid->norm1(isReadFromDisk);
              }

              printLocGlobAvgsRelatives("l1 error of component field before combination", 
                       sgError, norm1Full, uvsgPerGrid->numElts(), grid, myCommWorld, gId);
           }
        }
} //printErrorsComponentFields()


/////////////////////////////////////////////////////////////////////////////////////////
void printErrorsCombinedFields(void) {
        if (nC > 0) {
            double sgError = 0.0, norm1Full = 0.0;
            if (sg->myrank >= 0) { // this process holds part of the sparse grid
               sgError = adv->checkErrorLBM(uvsg, usgFullGrid, gc->is2D(),
                                             isReadFromDisk);
               norm1Full = usgFullGrid->norm1(isReadFromDisk);
            }
            printLocGlobAvgsRelatives("l1 error of combined field", 
                                     sgError, norm1Full, uvsg->numElts(), grid, myCommWorld);

            if (verbosity > 2 && sg->myrank >= 0)
               uvsg->print(rank, "combined field");
        }

        if (verbosity > 2 && sg->myrank >= 0) {
           usgFullGrid->print(rank, "full grid field");
           fflush(stdout);
        }

        if (verbosity > 2) {
           for (int hi = 0; hi < nSpecies; hi++) 
              u[hi]->print(rank, "final field");
        }
} //printErrorsCombinedFields()


void printElapsedTimes(void) {
	// Measuring elapsed times
	// Max execution time of the application (excluding full grid computation)
	double simulateLBMMaxTotalTime;
	MPI_Reduce(&simulateLBMTotalTime, &simulateLBMMaxTotalTime, 1, MPI_DOUBLE, MPI_MAX, 0, myCommWorld);

	// Max execution time of gatherScatter
	double gatherScatterMaxTotalTime;
	MPI_Reduce(&gatherScatterTotalTime, &gatherScatterMaxTotalTime, 1, MPI_DOUBLE, MPI_MAX, 0, myCommWorld);

	// Max execution time of simulateLBM
	if (rank == 0)
	   printf("\nMaximum execution time (MPI_Wtime) of simulateLBM is: %0.6f sec\n\n",
			     simulateLBMMaxTotalTime);

	// Max execution time of gatherScatter
	if (rank == 0)
	   printf("\nMaximum execution time (MPI_Wtime) of gatherScatter is: %0.6f sec\n\n",
			     gatherScatterMaxTotalTime);

	// Max execution time of the application (excluding full grid computation)
	if (rank == 0)
	   printf("\nMaximum execution time (MPI_Wtime) of the application "
		  "(excluding the full grid computation, and communicator "
		  "reconstruction, if any) is: %0.6f sec\n\n",
		  simulateLBMMaxTotalTime + gatherScatterMaxTotalTime);

	// Max execution time for computing full grid
	double fullGridMaxTotalTime;
	MPI_Reduce(&fullGridTotalTime, &fullGridMaxTotalTime, 1, MPI_DOUBLE, MPI_MAX, 0, myCommWorld);
	if (rank == 0 && verbosity > 0)
	   printf("\nMaximum execution time (MPI_Wtime) only of the full grid computation is: %0.6f sec\n\n",
		   fullGridMaxTotalTime);
} //printElapsedTimes() 

void visualizeSparseGrid (void) {
        int gridSizeX = (1 << grid.x);
        int gridSizeY = (1 << grid.y);
        int gridSizeZ = (1 << grid.z);
        int fieldSize = gridSizeX * gridSizeY * gridSizeZ;

        int sgIndex = 0, lx, ly, lz;

        if (gc->is2D()) {
           assert(gridSizeZ == 1);
        }

        LbmFiSG = new double[fieldSize];

        lz = uvsg->uh->l.z;
        ly = uvsg->uh->l.y;
        lx = uvsg->uh->l.x;
        if (gc->is2D())
           assert(lz == 1);
        else
           lz = (lz % 2 == 0)? lz: lz - 1;
        ly = (ly % 2 == 0)? ly: ly - 1;
        lx = ((lx/B) % 2 == 0)? lx: lx - B;
        // For 2D only
        for (int k = 0; k < lz; k++) {
           for (int j = 0; j < ly; j++) {
              for (int i = 0; i < lx; i++) {
                 LbmFiSG[sgIndex] = Vh(uvsg->uh, i, j, k);
                 sgIndex++;
              }
           }
        }

        fortranCommWorld = MPI_Comm_c2f(combinationComm);
        visualizeField(&(sg->P.x), &(sg->P.y), &(sg->P.z), &fortranCommWorld,
                 &gridSizeX, &gridSizeY, &gridSizeZ, &fieldSize, LbmFiSG);

        delete [] LbmFiSG;

        // Call python script for visualization
        if (sg->myrank == 0) {
           int ret = system("python ../../src/testing/rho_visualize_lbm.py"); // May not work on PBS (in that case, manually do this)
        }
} //visualizeSparseGrid()        

