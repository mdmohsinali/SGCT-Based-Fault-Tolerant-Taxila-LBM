/* SGCT-Based-Fault-Tolerant-Taxila-LBM Code source file.
   Copyright (c) 2015, Md Mohsin Ali. All rights reserved.
   Licensed under the terms of the BSD License as described in the LICENSE_FT_CODE file.
   This comment must be retained in any redistributions of this source file.
*/

/* 
 * File       : LbmInterface.h
 * Description: Contains function prototypes to wrapper routines to be called 
 *              from C or Fortran, and a routine which creates individual
 *              directories where changed parameters are copied from a provided
 *              base parameter
 * Author     : Mohsin Ali
 * Created    : April 2014
 * Updated    : November 2014
 */

#ifndef LBMINTERFACE_INCLUDED
#define	LBMINTERFACE_INCLUDED

// Header files
#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>
#include<cstring>
#include<iostream>
#include<sys/stat.h>
#include<sys/types.h>

// Function prototypes
// Prototypes to wrapper routines to be called from C or Fortran
#ifdef __cplusplus
extern "C"{
#endif
   void c_runLBM(char *, int *, int *, int *, MPI_Fint *, bool *, bool *, int *,
            int * nx, int * ny, int * nz, int * n_comb, bool * isRealCall);
   void visualizeField(int * proc_x, int * proc_y, int * proc_z, MPI_Fint *,
            int * nx, int * ny, int * nz,
            int * lbm_size, double * lbm_field);
   void c_get_LBM_field(int * ndims, double * fi_2d, double * fi_3d,
            int * x_start, int * y_start, int * z_start, 
            int * x_width, int * y_width, int * z_width,
            bool * isComponentGrid, bool * isRepeatComb);
#ifdef __cplusplus
}
#endif
// Other prototypes
void setParameters(char * parInDir, int n_procs_x_dim, int n_procs_y_dim, int n_procs_z_dim,
                   int n_procs_sim, int n_x_dim, int n_y_dim, int n_z_dim, 
                   int degOfParall, int nC, bool is2D, bool isCompGrid);
int getTimeStepping(char * fName, char * iStr);
char * trim(char * str);
char * trimLeadingSpaces(char * str);

/////////////////////////////////////////////////////////////////////////////////////////
// Function definitions
void setParameters(char * parInDir, int n_procs_x_dim, int n_procs_y_dim, int n_procs_z_dim, 
        int n_procs_sim, int n_x_dim, int n_y_dim, int n_z_dim, int degOfParall, int nC, 
        bool is2D, bool isCompGrid) {
    char fileName[128] = " ";
    int MAX_PATH_LENGTH = 80;
    char path[MAX_PATH_LENGTH];
    char lineStringOfFile[128] = " ";
    char * lineStringOfFilePointer = NULL;
    char replacedLineStringOfFile[128] = " ";
    FILE *inputFilePointer, *outputFilePointer;
    char * cwdPointer = getcwd(path, MAX_PATH_LENGTH);

    // Set input file name
    char * inputFileName = (char *)"input_data";

    if(cwdPointer){
       sprintf(fileName, "%s/%s", path, inputFileName);
    }
    else{
       printf("getcwd() of fileName failed\n\n");
       exit(1);
    }

    if((inputFilePointer = fopen(fileName,"r")) == NULL) {
       printf("inputFilePointer: file opening error!\n");
    }

    // Creating "path/param_grid_*" directory
    mkdir(parInDir, S_IRWXU|S_IRGRP|S_IXGRP);

    if(cwdPointer){
       //sprintf(fileName, "%s/checkpoint_of_process_%d", path, rank);
       sprintf(fileName, "%s/%s/%s", path, parInDir, inputFileName);
    }
    else{
       printf("getcwd() of fileName failed\n\n");
       exit(1);
    }

    if((outputFilePointer = fopen(fileName,"w")) == NULL) {
       printf("outputFilePointer: file opening error!\n");
    }

    // Confirming z dimension parallelization for 2D
    if (is2D)
       assert(n_procs_z_dim == 1);

    // Running full grid LBM for whole timesteps
    // There is no division of timesteps for multiple combinations for full grid
    if (not isCompGrid)
       nC = 1;

    // Assuming there is NO space at the beginning of each entry (line) on the file opened
    while (fgets(lineStringOfFile, sizeof(lineStringOfFile), inputFilePointer) != NULL) {
        lineStringOfFilePointer = trimLeadingSpaces(lineStringOfFile);
	if (strncmp("-da_processors_x", lineStringOfFilePointer, strlen("-da_processors_x")) == 0){
	   sprintf(replacedLineStringOfFile,"#-da_processors_x %d\n", n_procs_x_dim); // Not taking this input from */input_data file
	   fputs(replacedLineStringOfFile, outputFilePointer);
	}
	else if (strncmp("-da_processors_y", lineStringOfFilePointer, strlen("-da_processors_y")) == 0){
	   sprintf(replacedLineStringOfFile,"#-da_processors_y %d\n", n_procs_y_dim); // Not taking this input from */input_data file
	   fputs(replacedLineStringOfFile, outputFilePointer);
	}
	else if (strncmp("-da_processors_z", lineStringOfFilePointer, strlen("-da_processors_z")) == 0) {
	   sprintf(replacedLineStringOfFile,"#-da_processors_z %d\n", n_procs_z_dim); // Not taking this input from */input_data file
	   fputs(replacedLineStringOfFile, outputFilePointer);
	}
	else if (strncmp("-NX", lineStringOfFilePointer, strlen("-NX")) == 0) {
	   sprintf(replacedLineStringOfFile,"#-NX %d\n", n_x_dim); // Not taking this input from */input_data file
	   fputs(replacedLineStringOfFile, outputFilePointer);
	}
	else if (strncmp("-NY", lineStringOfFilePointer, strlen("-NY")) == 0) {
	   sprintf(replacedLineStringOfFile,"#-NY %d\n", n_y_dim); // Not taking this input from */input_data file
	   fputs(replacedLineStringOfFile, outputFilePointer);
	}
	else if (strncmp("-NZ", lineStringOfFilePointer, strlen("-NZ")) == 0) {
	   sprintf(replacedLineStringOfFile,"#-NZ %d\n", n_z_dim); // Not taking this input from */input_data file
	   fputs(replacedLineStringOfFile, outputFilePointer);
	}
	else if (strncmp("-npasses", lineStringOfFilePointer, strlen("-npasses")) == 0) {
	   sprintf(replacedLineStringOfFile,"-npasses %d\n",
		   getTimeStepping(inputFileName, (char *)"-npasses")/nC);
	   fputs(replacedLineStringOfFile, outputFilePointer);
	}
	else if (strncmp("-kwrite", lineStringOfFilePointer, strlen("-kwrite")) == 0) {
	   sprintf(replacedLineStringOfFile,"-kwrite %d\n",
		   getTimeStepping(inputFileName, (char *)"-kwrite")/nC);
	   fputs(replacedLineStringOfFile, outputFilePointer);
	}
	else {
	   fputs(lineStringOfFilePointer, outputFilePointer);
	}
    }

    fclose(inputFilePointer);
    fclose(outputFilePointer);
}//setParameters()

/////////////////////////////////////////////////////////////////////////////////////////
int getTimeStepping(char * fName, char * iStr) {
   FILE *f;
   char string[1000], seps[] = " \n,( )!=";
   char *p;
   f = fopen(fName,"r");
   if(!f) {
      printf("Probably executing without an input file. Returning %s = 1.\n", iStr);
      return 1; // return ntimesteps = 1
   }

   while(fgets(string, sizeof(string)-1, f) != NULL) {
      // Break into tokens
      p = trim(string);

      // Find first token
      // It is obvious from the input file that iStr
      // is not the first token
      p = strtok(string, seps); 

      while(p != NULL){
         //printf("Token: %s\n", p);
         p = strtok(NULL, seps); // find next token

         if (strncmp(iStr, string, strlen(iStr)) == 0) {
            // iStr token is found
            // our "desired token" (a number) is just after this
            fclose(f); // close the opened file
            return atoi(p); // convert string into int and return
         }
      }
   }
   fclose(f); // close the opened file
   return 0; // some compilers give warning without this
} //getTimeStepping


/////////////////////////////////////////////////////////////////////////////////////////
char * trim(char * str) {
    size_t len = 0;
    char *frontp = str - 1;
    char *endp = NULL;

    if( str == NULL )
       return NULL;

    if(str[0] == '\0')
       return str;

    len = strlen(str);
    endp = str + len;

    // Move the front and back pointers to address
    // the first non-whitespace characters from
    // each end.
     
    while(isspace(*(++frontp)));
    while(isspace(*(--endp)) && endp != frontp);

    if( str + len - 1 != endp )
       *(endp + 1) = '\0';
    else if(frontp != str &&  endp == frontp)
       *str = '\0';

    // Shift the string so that it starts at str so
    // that if it's dynamically allocated, we can
    // still free it on the returned pointer.  Note
    // the reuse of endp to mean the front of the
    // string buffer now.     
    endp = str;
    if( frontp != str ) {
       while(*frontp) *endp++ = *frontp++;
          *endp = '\0';
    }

    return str;
} //trim()


/////////////////////////////////////////////////////////////////////////////////////////
char * trimLeadingSpaces(char * testString) {
    char * tmp = testString;
    int i;
 
    while (*tmp == ' ')
       *tmp++;

    i = 0;
    while (*tmp != '\0')
       testString[i++] = *tmp++;

    testString[i] = '\0';

    return testString;
} //trimLeadingSpaces()

#endif	/* LBMINTERFACE_INCLUDED */


