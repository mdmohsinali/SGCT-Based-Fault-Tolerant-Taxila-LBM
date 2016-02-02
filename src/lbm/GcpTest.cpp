/*
   Copyright (c) 2014, Brendan Harding. All rights reserved.
*/

/*

Brendan Harding
October 2014

This is a test for the routines in gcp_general.h
it also serves as an example of usage

An installation on GNU GLPK is required.
This code has been tested with glpk-4.55 on both gnu and intel compilers.

Compiles and runs on my mac with: 
g++ -lglpk -O2 gcp_test.cpp 
./a.out

*/

#include "gcp_general.h"

int quick_test(int d=2){
	// set d to 2,3,4 to test in 2,3,4 dimensions respectively
	printf("\n######################\nQuick_test: d=%d\n######################\n",d);
	// define some grids and print
	int grids[24]={6,4,2,5,3,6,4,8,6,2,5,4,2,5,6,4,1,9,5,6,4,6,5,3};
	int ng=24/d;
	printf("Grids: (%d)\n",ng);
	for (int i=0;i<ng;i++) {
		printf("\t");
		for (int k=0;k<d;k++) printf("%d ",grids[i*d+k]);
		printf("\n");
	}
	// get the (local) downset from the given grids
	int nng;
	int* newgrids=generate_downset(grids,ng,d,nng); // newgrids allocated here
	printf("Downset Grids: (%d)\n",nng);
	for (int i=0;i<nng;i++) {
		printf("\t");
		for (int k=0;k<d;k++) printf("%d ",newgrids[i*d+k]);
		printf("\n");
	}
	// initialise coefficient and status arrays
	int* coef=new int[nng];
	int* stat=new int[nng];
	for (int i=0;i<nng;i++) stat[i]=1; 
	// want to set stat to 0 for grids in the original list
	for (int i=0;i<nng;i++) {
		for (int j=0;j<ng;j++) {
			int s=0;
			for (int k=0;k<d;k++) {
				if (grids[j*d+k]==newgrids[i*d+k]) s++;
				else break;
			}
			if (s==d) {stat[i]=0;break;};
		}
	}
	// now call the gcp function to find the coefficients
	int ret=gcp(newgrids,stat,d,nng,coef);
	// check for success and print results
	if (ret==0) printf("success!\n");
	else printf("mip error\n");
	// now check that the coefficients sum to 1
	int s=0;
	for (int i=0;i<nng;i++) s+=coef[i];
	if (s!=1) printf("sum error: %d!=1\n",s);
	// check that coef is 0 for grids with status 1
	for (int i=0;i<nng;i++) {
		if (stat[i]==1 && coef[i]!=0) {
			printf("status error: stat[%d]==%d and coef[%d]=%d!=0\n",i,stat[i],i,coef[i]);
			break;
		}
	}
	// optional print the coefficients and status'
	/*
	printf("coeff (and stat): \n");
	for (int i=0;i<nng;i++) {
		printf("\t%d\t(%d)\n",coef[i],stat[i]);
	}
	*/
	// cleanup
	delete[] newgrids;
	delete[] coef;
	delete[] stat;
	return 0;
}


int comprehensive2d_test(int n=3){
	// set n for the level
	printf("\n######################\nComprehensive 2d test: n=%d\n######################\n",n);
	// some parameters, dimension, level
	int d=2;
	// generate list of grids in downset up to given level
	int ng=0;
	for (int l=0;l<=n;l++) ng+=(l+1);
	int* grids=new int[d*ng];
	int ind=0;
	for (int l=0;l<=n;l++) {
		for (int i=0;i<=l;i++) {
			grids[ind++]=i;
			grids[ind++]=l-i;
		}
	}
	// print the grids
	printf("Grids: (%d)\n",ng);
	for (int i=0;i<ng;i++) {
		printf("\t");
		for (int k=0;k<d;k++) printf("%d ",grids[i*d+k]);
		printf("\n");
	}
	// as we have a downset there is no need to call generate downset
	// initialise coefficient and status array
	int* coef=new int[ng];
	int* stat=new int[ng];
	for (int i=0;i<ng;i++) stat[i]=0;
	// now begin loop over all possible status'
	int its=(1<<ng)-1; // note the -1 is to skip the case with all status' 1
	int mip_errors=0;
	int sum_errors=0;
	int stat_errors=0;
	for (int it=0;it<its;it++) {
		// optionally print the current stat being solver
		//for (int i=0;i<ng;i++) printf("%d ",stat[i]);
		//printf("\n");
		// call gcp to solve current problem
		int ret=gcp(grids,stat,d,ng,coef,0);
		if (ret==0) {
			//printf("success!\n");
			// now check that the coefficients sum to 1
			int s=0;
			for (int i=0;i<ng;i++) s+=coef[i];
			if (s!=1) {
				sum_errors++;
			} else {
				// check that coef is 0 for grids with status 1
				for (int i=0;i<ng;i++) {
					if (stat[i]==1 && coef[i]!=0) {
						stat_errors++;
						break;
					}
				}
			}
		} else {
			//printf("error\n");
			mip_errors++;
		}
		// update stat
		stat[0]++;
		for (int i=0;i<ng-1;i++) if (stat[i]==2) {stat[i]=0;stat[i+1]++;};
	}
	// print total number of errors
	printf("total errors: %d %d %d\n",mip_errors,sum_errors,stat_errors);
	// cleanup
	delete[] grids;
	delete[] coef;
	delete[] stat;
	return 0;
}


int comprehensive3d_test(int n=2){
	// set n for the level
	printf("\n######################\nComprehensive 3d test: %d\n######################\n",n);
	// some parameters, dimension
	int d=3;
	// generate list of grids in downset up to given level
	int ng=0;
	for (int l=0;l<=n;l++) ng+=(l+2)*(l+1)/2;
	int* grids=new int[d*ng];
	int ind=0;
	for (int l=0;l<=n;l++) {
		for (int k=0;k<=l;k++) {
			for (int i=0;i<=l-k;i++) {
				grids[ind++]=i;
				grids[ind++]=k;
				grids[ind++]=l-i-k;
			}
		}
	}
	// print the grids
	printf("Grids: (%d)\n",ng);
	for (int i=0;i<ng;i++) {
		printf("\t");
		for (int k=0;k<d;k++) printf("%d ",grids[i*d+k]);
		printf("\n");
	}
	// as we have a downset there is no need to call generate downset
	// initialise coefficient and status array
	int* coef=new int[ng];
	int* stat=new int[ng];
	for (int i=0;i<ng;i++) stat[i]=0;
	// now begin loop over all possible status'
	int its=(1<<ng)-1; // note the -1 is to skip the case with all status' 1
	int mip_errors=0;
	int sum_errors=0;
	int stat_errors=0;
	for (int it=0;it<its;it++) {
		// call gcp to solve current problem
		int ret=gcp(grids,stat,d,ng,coef,0);
		if (ret==0) {
			// now check that the coefficients sum to 1
			int s=0;
			for (int i=0;i<ng;i++) s+=coef[i];
			if (s!=1) {
				sum_errors++;
			} else {
				// check that coef is 0 for grids with status 1
				for (int i=0;i<ng;i++) {
					if (stat[i]==1 && coef[i]!=0) {
						stat_errors++;
						break;
					}
				}
			}
		} else {
			mip_errors++;
		}
		// update stat
		stat[0]++;
		for (int i=0;i<ng-1;i++) if (stat[i]==2) {stat[i]=0;stat[i+1]++;};
	}
	// print total number of errors
	printf("total errors: %d %d %d\n",mip_errors,sum_errors,stat_errors);
	// cleanup
	delete[] grids;
	delete[] coef;
	delete[] stat;
	return 0;
}


int main() {

	quick_test(2);
	quick_test(3);
	comprehensive2d_test();
	comprehensive3d_test();

	return 0;
	
};

