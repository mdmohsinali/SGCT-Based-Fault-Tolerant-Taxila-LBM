/*
   Copyright (c) 2014, Brendan Harding. All rights reserved.
*/

/*

Brendan Harding
October 2014

Routines for the general sparse grid coefficient problem,
gcp(...) finds combination coefficients for an arbitrary
collection of grids in arbitrary dimensions.

An installation on GNU GLPK is required.
This code has been tested with glpk-4.55 on both gnu and intel compilers.

Compiles and runs on my mac with: 
g++ -O2 -lglpk GcpGeneral.cpp
or
g++ GcpGeneral.cpp -O2 -lglpk
./a.out

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <glpk.h>



// gcp: solves the general coefficient problem
int gcp(int* grids, int *status, int ndims, int ngrids, int* coeff, int verbosity = 1) {
	/*
	gcp: solves the general coefficient problem
	
	Note: it is expected that the given grids form a (local) downset
		it is unlikely to give sensible results if this is not the case
		(use the generateDownset() function to create a downset)
	grids is assumed to have length ndims*ngrids, it has the format 
		grid0[0],...,grid0[d-1],grid1[0],...,grid1[d-1],grid2[0],...
	status is assumed to have length ngrids, it gives the status of each grid
		status[i]=1 if the i'th grid is NOT to be used in the combination
	coef is assumed to have length ngrids and will be filled with the
		resulting combination coefficients (order consistent with grids)
	returns 0 if successful and 1 otherwise
	*/
	// initialise MIP problem
	glp_prob *mip;
	mip = glp_create_prob();
	glp_set_obj_dir(mip, GLP_MIN);
	glp_set_prob_name(mip, "gcp");

	if (verbosity == 0) 
		glp_term_out(GLP_OFF); // set to off for no output

	glp_add_cols(mip, ngrids);

	int dp2 = 1<<ndims;

	// calculate number of equality constraints
	int num_eq_cons = 0;
	for (int i = 0; i < ngrids; ++i) 
		num_eq_cons += status[i];

	// create matrix going from hierarchical to combination coefficients
	int *M = new int[ngrids*ngrids];
	for (int i = 0; i < ngrids; ++i) {
		for (int j = 0; j < ngrids; ++j) {
			M[i*ngrids+j] = 0;
			int c = 0;
			bool neighbour = true;
			for (int d = 0; d < ndims; ++d) {
				if (grids[j*ndims+d] == grids[i*ndims+d] || grids[j*ndims+d] == 1+grids[i*ndims+d]) {
					c += grids[j*ndims+d] - grids[i*ndims+d];
				} else {
					neighbour = false;
					break;
				}
			}
			if (neighbour) 
				M[i*ngrids+j] = pow(-1, c);
		}
	}

	// Now set up the function to minimise (and set variables as binary)
	for (int i = 0; i < ngrids; ++i) {
		int s = 0;
		for (int d = 0; d < ndims; ++d) 
			s += grids[i*ndims+d];
		glp_set_obj_coef(mip, i+1, -pow(0.25,s));
		glp_set_col_kind(mip, i+1, GLP_BV); // set parameter as binary variable
	}

	// initialise vectors for the problem matrix
	// number of constraints for the downset condition is no more than ndims*ngrids with 2 entries for each constraint
	// number of constraints for the status conditions is num_eq_cons with at most dp2 entries per constraint
	int *ia = new int[1+num_eq_cons*dp2+ngrids*ndims*2];
	int *ja = new int[1+num_eq_cons*dp2+ngrids*ndims*2];
	double* va = new double[1+num_eq_cons*dp2+ngrids*ndims*2];
	glp_add_rows(mip, num_eq_cons+ngrids*ndims); 
	int ind=1;
	int row=1;

	// Now set up the equality constraints
	for (int i = 0; i < ngrids; ++i) {
		if (status[i]==1) {
			for (int j = 0; j < ngrids; ++j) {
				if (M[i*ngrids+j] != 0) {
					ia[ind] = row;
					ja[ind] = j+1;
					va[ind] = M[i*ngrids+j];
					++ind;
				}
			}
			glp_set_row_bnds(mip,row,GLP_FX,0.0,0.0); // FX for fixed (equality constraint)
			++row;
		}
	}
	//Now add the downset constraints to this
	for (int i = 0; i < ngrids; ++i) {
		for (int j = 0; j < ngrids; ++j) {
			int s = 0;
			for (int d = 0; d < ndims; ++d) {
				if (grids[i*ndims+d] <= grids[j*ndims+d] && s <= 1) {
					s += grids[j*ndims+d] - grids[i*ndims+d];
				} else {
					s=0;
					break;
				}
			}
			if (s==1) {
				ia[ind] = row;
				ja[ind] = j+1;
				va[ind] = 1.0;
				++ind;
				ia[ind] = row;
				ja[ind] = i+1;
				va[ind] = -1.0;
				++ind;
				glp_set_row_bnds(mip, row, GLP_UP, 0.0, 0.0); // UP for upper bound
				++row;
			}
		}
	}
	ind--; // now ind is the total (non-zero) entries in the problem matrix
	row--; // now row is the total rows in the problem matrix

	// delete the unnecessary rows (not strictly necessary)
	const int ndr = num_eq_cons + ngrids*ndims-row;
	if (ndr>0) {
		int num[1+ndr];
		for (int i = 0; i < ndr; ++i) 
			num[i+1] = row+i+1;
		glp_del_rows(mip,ndr,num);
	};

	// now load the problem matrix
	glp_load_matrix(mip, ind, ia, ja, va); 

	// set some solver parameters, or just do a simplex presolve?
	glp_adv_basis(mip, 0);
	glp_simplex(mip, NULL);

	// call the integer programming solver
	glp_intopt(mip, NULL);
	int ret=glp_mip_status(mip);
	if (ret==GLP_OPT) {
		//double z=glp_mip_obj_val(mip); // the minimum of the result
		// get solution
		int h[ngrids];
		for (int i = 0; i < ngrids; ++i) {
			h[i] = glp_mip_col_val(mip, i+1);
		}
		// convert to combination coefficients
		for (int i = 0; i < ngrids; ++i) {
			int c = 0;
			for (int j = 0; j < ngrids; ++j) {
				c += M[i*ngrids+j]*h[j];
			}
			coeff[i] = c;
		}
		glp_delete_prob(mip);
		delete[] M;
		delete[] ia;
		delete[] ja;
		delete[] va;
		return 0;
	} else {
		printf("Error: gcp: GLPK could not solve problem\n");
		glp_delete_prob(mip);
		delete[] M;
		delete[] ia;
		delete[] ja;
		delete[] va;
		return 1;
	}
}




// downset_recursive: recursive function to count/generate a (local) downset
// it is far from optimal but will do the job for now
// Note: users should not call this directly!
//       call generateDownset instead
int downset_recursive(int* grids, int ngrids, int ndims, int* min, int* max, int* current, int dim, int *newgrids = NULL, int ind = 0) {
	/*
	downset_recursive: recursive function to count/generate a (local) downset
	
	grids is assumed to have length ndims*ngrids, it has the format 
		grid0[0],...,grid0[d-1],grid1[0],...,grid1[d-1],...
	min,max,current are of length ndims giving the min,max,current index in each dimension
	dim is the current dim to loop over (from 0 to ndim-1)
	newgrids is potentially where new grids are added with ind being the current index
	ind keeps track of the current index of the grid in the downset
	*/
	// loop over current dimension
	for (int i = min[dim]; i <= max[dim]; ++i) {
		current[dim] = i;
		if (dim > 0) 
			ind = downset_recursive(grids, ngrids, ndims, min, max, current, dim-1, newgrids, ind);
		else {
			for (int j = 0; j < ngrids; ++j) {
				int s = 0;
				for (int k = 0; k < ndims; ++k) 
					if (current[k] <= grids[j*ndims+k]) 
						++s;
				if (s == ndims) {
					if (newgrids != NULL) {
						for (int k = 0; k < ndims; ++k) 
							newgrids[ind*ndims+k] = current[k];
					}
					++ind;
					break;
				}
			}
		}
	}
	return ind;
}



// generateDownset: generates a (local) downset from a list of grids
int* generateDownset(int* grids, int ngrids, int ndims, int &new_len) {
	/*
	generateDownset: generates a (local) downset from a list of grids
	
	grids is of length ndims*ngrids, it has the format 
		grid0[0],...,grid0[d-1],grid1[0],...,grid1[d-1],grid2[0],...
	new_len returns with the number of grids in the (local) downset
	returns pointer to the array of new grids (total length: ndims*new_len)
	the caller is responsible for deallocation of the returned array
	returns 1 on failure and 0 with success
	*/
	// get min and max indices in each dimension
	int min[ndims];
	int max[ndims];
	for (int k = 0; k < ndims; ++k) {
		min[k] = 100;
		max[k] = 0;
	};
	for (int i = 0; i < ngrids; ++i) {
		for (int k = 0; k < ndims; ++k) {
			if (grids[i*ndims+k] < min[k]) 
				min[k] = grids[i*ndims+k];
			if (grids[i*ndims+k] > max[k]) 
				max[k] = grids[i*ndims+k];
		}
	}

	// now calculate the number of grids in the (local) downset
	int current[ndims];
	new_len = downset_recursive(grids,ngrids,ndims,min,max,current,ndims-1);

	// now allocate data and generate the list of grids in the downset
	int* newgrids = new int[ndims*new_len];
	downset_recursive(grids, ngrids, ndims, min, max, current, ndims-1, newgrids);
	return newgrids;
}



// a generic wrapper for setting up the gcp problem for a given list of grids
int computeCoefficients(int* grids, int ngrids, int ndims, int* status, double* coeff) {
	/*
	compute_coefficients: a wrapper forsetting up and solving the general coefficient problem given a list of grids
	
	ngrids is the number of grids
	ndims is the number of dimensions
	grids contains the grid levels and is of length ndims*ngrids, 
		it has the format grids[grid_number*ndims+dimension_index],
		there should NOT be any duplicate grids in the list
	status has length ngrids, it status[i]=1 means the i'th grid is considered as failed, status[i]=0 otherwise
	coeff should be preallocated with length ngrids and is where the calculated coefficients will be stored
	*/

	// get the (local) downset from the given grids
	int nng;
	int* newgrids = generateDownset(grids, ngrids, ndims, nng); // newgrids is allocated here
	if (nng > 1000) 
		printf("Warning: compute_coefficients: problem size is large an may take some time to solve: %d\n",nng);

	// initialise enlarged coefficient and status arrays
	int* temp_coef = new int[nng];
	int* temp_stat = new int[nng];
	int* mapping = new int[ngrids];
	for (int i = 0; i < nng; ++i) 
		temp_stat[i]=1;

	// want to set stat to 0 for grids in the original list
	for (int i = 0; i < ngrids; ++i) {
		for (int j = 0; j < nng; ++j) {
			int s = 0;
			for (int k = 0; k < ndims; ++k) {
				if (grids[i*ndims+k] == newgrids[j*ndims+k]) 
					++s;
				else 
					break;
			}
			if (s == ndims) {
				temp_stat[j] = status[i];
				mapping[i] = j;
				break;
			};
		}
	}

	// now call the gcp function to find the coefficients
	int ret = gcp(newgrids, temp_stat, ndims, nng, temp_coef, 0);

        // some status variables
        bool statusReturnSuccess = false, 
	     statusSumSuccess = false, 
	     statusCoeff0Status1Success = false;

	// check for success and print results
	if (ret != 0) {
		statusReturnSuccess = false;
		printf("Error: compute_coefficients: GLPK mip error\n");
		delete[] newgrids;
		delete[] temp_coef;
		delete[] temp_stat;
		delete[] mapping;
		return 1;                
	}
        else {
		statusReturnSuccess = true;
        }

	// perform a quick check that the coefficients are sensible
	// first check that the coefficients sum to 1
	int s = 0;
	for (int i = 0; i < nng; ++i) 
		s+=temp_coef[i];
	if (s != 1) {
		statusSumSuccess = false;
		printf("Error: compute_coeficients: sum error: %d != 1\n",s);
		delete[] newgrids;
		delete[] temp_coef;
		delete[] temp_stat;
		delete[] mapping;
		return 1;
	}
        else {
		statusSumSuccess = true;
        }
	// now check that coef is 0 for grids with status 1
	for (int i = 0; i < nng; ++i) {
		if (temp_stat[i] == 1 && temp_coef[i] != 0) {
			statusCoeff0Status1Success = false;
			printf("Error: compute_coefficients: status error: stat[%d]==%d and coef[%d]=%d!=0\n",i,temp_stat[i],i,temp_coef[i]);
			delete[] newgrids;
			delete[] temp_coef;
			delete[] temp_stat;
			delete[] mapping;
			return 1;
		}
		else {
			statusCoeff0Status1Success = true;
		}
	}
	// if everything is good, then copy coefficients from temp_coef to coeff
	if (statusReturnSuccess && statusSumSuccess && statusCoeff0Status1Success) {
		for (int i = 0; i < ngrids; ++i) 
	        	coeff[i] = (double) temp_coef[mapping[i]];
        }

	// cleanup
	delete[] newgrids;
	delete[] temp_coef;
	delete[] temp_stat;
	delete[] mapping;
	return 0;
}


