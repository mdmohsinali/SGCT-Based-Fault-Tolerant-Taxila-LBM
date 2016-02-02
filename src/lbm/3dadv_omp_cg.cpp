/* ParSGCT Code source file.
Copyright (c) 2015 Peter Edward Strazdins. All rights reserved.
Licensed under the terms of the BSD License as described in the LICENSE_SGCT file.
This comment must be retained in any redistributions of this source file.
*/

// adopted from Brendan Harding's 3dadv codes 08/13 to be callable
// for a distributed memory solver. This required:
//    - dx,dy,dz, to be passed as paraters (as set by global sizes)
//    - evolve only for a single time step (as a dist. memory update boundary
//      is required)

extern "C" {

void compute_godunov_scheme(double DT,double *u,int lsx,int lsy,int lsz,double Vx,double Vy,double Vz, double dx, double dy, double dz) {
#if 0
	int asx=lsx+2;
	int asy=lsy+2;
	// int asz=lsz+2;
	//Godunov fluxes for advection
	double *f_x=new double[lsx*lsy*lsz];
	double *f_y=new double[lsx*lsy*lsz];
	double *f_z=new double[lsx*lsy*lsz];
	int ix,iy,iz,kf,ku;
		if (Vx>=0.0) {		  
			for (iz=0;iz<lsz;iz++) {
				for (iy=0;iy<lsy;iy++) {
					kf=(iz*lsy+iy)*lsx;
					ku=((iz+1)*asy+(iy+1))*asx+1;
					for (ix=0;ix<lsx;ix++) {
						f_x[kf+ix]=Vx*(u[ku+ix]-u[ku+ix-1]);
					}
				}
			}
		} else {
			for (iz=0;iz<lsz;iz++) {
				for (iy=0;iy<lsy;iy++) {
					kf=(iz*lsy+iy)*lsx;
					ku=((iz+1)*asy+(iy+1))*asx+1;
					for (ix=0;ix<lsx;ix++) {
						f_x[kf+ix]=Vx*(u[ku+ix+1]-u[ku+ix]);
					}
				}
			}
		}
		if (Vy>=0.0) {
			for (iz=0;iz<lsz;iz++) {
				for (iy=0;iy<lsy;iy++) {
					kf=(iz*lsy+iy)*lsx;
					ku=((iz+1)*asy+(iy+1))*asx+1;
					for (ix=0;ix<lsx;ix++) {
						f_x[kf+ix]=Vx*(u[ku+ix]-u[ku+ix-asx]);
					}
				}
			}
		} else {
			for (iz=0;iz<lsz;iz++) {
				for (iy=0;iy<lsy;iy++) {
					kf=(iz*lsy+iy)*lsx;
					ku=((iz+1)*asy+(iy+1))*asx+1;
					for (ix=0;ix<lsx;ix++) {
						f_x[kf+ix]=Vx*(u[ku+ix+asx]-u[ku+ix]);
					}
				}
			}
		}
		//A diffusion part will possibly be added at a later stage
		//Now update u_loc
		for (iz=0;iz<lsz;iz++) {
			for (iy=0;iy<lsy;iy++) {
				kf=(iz*lsy+iy)*lsx;
				ku=((iz+1)*asy+(iy+1))*asx+1;
				for (ix=0;ix<lsx;ix++) {
					u[ku+ix]+=-DT/dx*f_x[kf+ix]-DT/dy*f_y[kf+ix]-DT/dz*f_z[kf+ix];
				}
			}
		}
		//Need to update the boundaries before movng on to the next time step
	delete f_x;
	delete f_y;
	delete f_z;
#else
	int asx=lsx+2;
	int asy=lsy+2;
	double *f_x=new double[lsx*lsy*lsz];
	double *f_y=new double[lsx*lsy*lsz];
	double *f_z=new double[lsx*lsy*lsz];
	int i,j,l,kf,ku;
 		if (Vx>=0.0) {
			for (l=0;l<lsz;l++) {
				kf=l*lsy*lsx-lsx;
				ku=(l+1)*asy*asx+1;
				for (i=0;i<lsy;i++) {
					kf+=lsx;
					ku+=asx;
					for (j=0;j<lsx;j++) {
						f_x[kf+j]=Vx*(u[ku+j]-u[ku+j-1]);
					}
				}
			}
		} else {
			for (l=0;l<lsz;l++) {
				kf=l*lsy*lsx-lsx;
				ku=(l+1)*asy*asx+1;
				for (i=0;i<lsy;i++) {
					kf+=lsx;
					ku+=asx;
					for (j=0;j<lsx;j++) {
						f_x[kf+j]=Vx*(u[ku+j+1]-u[ku+j]);
					}
				}
			}
		}
		if (Vy>=0.0) {
			for (l=0;l<lsz;l++) {
				kf=l*lsy*lsx-lsx;
				ku=(l+1)*asy*asx+1;
				for (i=0;i<lsy;i++) {
					kf+=lsx;
					ku+=asx;
					for (j=0;j<lsx;j++) {
						f_y[kf+j]=Vy*(u[ku+j]-u[ku+j-asx]);
					}
				}
			}
		} else {
			for (l=0;l<lsz;l++) {
				kf=l*lsy*lsx-lsx;
				ku=(l+1)*asy*asx+1;
				for (i=0;i<lsy;i++) {
					kf+=lsx;
					ku+=asx;
					for (j=0;j<lsx;j++) {
						f_y[kf+j]=Vy*(u[ku+j+asx]-u[ku+j]);
					}
				}
			}
		}
		if (Vz>=0.0) {
			for (l=0;l<lsz;l++) {
				kf=l*lsy*lsx-lsx;
				ku=(l+1)*asy*asx+1;
				for (i=0;i<lsy;i++) {
					kf+=lsx;
					ku+=asx;
					for (j=0;j<lsx;j++) {
						f_z[kf+j]=Vz*(u[ku+j]-u[ku+j-asx*asy]);
					}
				}
			}
		} else {
			for (l=0;l<lsz;l++) {
				kf=l*lsy*lsx-lsx;
				ku=(l+1)*asy*asx+1;
				for (i=0;i<lsy;i++) {
					kf+=lsx;
					ku+=asx;
					for (j=0;j<lsx;j++) {
						f_z[kf+j]=Vz*(u[ku+j+asx*asy]-u[ku+j]);
					}
				}
			}
		}
		//A diffusion part will possibly be added at a later stage
		//Now update u_loc
		for (l=0;l<lsz;l++) {
			kf=l*lsy*lsx-lsx;
			ku=(l+1)*asy*asx+1;
			for (i=0;i<lsy;i++) {
				kf+=lsx;
				ku+=asx;
				for (j=0;j<lsx;j++) {
					u[ku+j]+=-DT/dx*f_x[kf+j]-DT/dy*f_y[kf+j]-DT/dz*f_z[kf+j];
				}
			};
		}

	delete f_x;
	delete f_y;
	delete f_z;
#endif
}


// Note: the use of const is important for the openmp blocks.
// Without const, openmp has a lot of overhead in synchronising variables
void compute_LW_scheme(const double DT,double *u,const int lsx,const int lsy,const int lsz,const double Vx,const double Vy,const double Vz,
const double dx, const double dy, const double dz) {
	const int asx=2+lsx;
	const int asy=2+lsy;
	const int asz=2+lsz;
	const int asxy=asx*asy;
	const int asxm1ym1=(asx-1)*(asy-1);
	//Lax Wendroff approach
	double *u_half=new double[(asx-1)*(asy-1)*(asz-1)];
	// We wrap the whole time stepping in a parallel loop to save
	// repeated opening and closing of omp threads.
	// This should reduce overhead and also gets around a bug in Apple's gcc.
	// In this version we use a "coarse-grained" omp approach.
	// It should have less overhead/better scalability
	#pragma omp parallel
	{
	int ix,iy,iz,ku,kh;
	double fx,fy,fz;

	const int myOMPRank = omp_get_thread_num();
	const int numOMPThreads = omp_get_num_threads();
	const int iStart1 = (myOMPRank*(asz-1))/numOMPThreads;
	const int iEnd1 = ((myOMPRank+1)*(asz-1))/numOMPThreads;
	const int iStart2 = (myOMPRank*lsz)/numOMPThreads;
	const int iEnd2 = ((myOMPRank+1)*lsz)/numOMPThreads;
		//first the half step
		for (iz=iStart1;iz<iEnd1;iz++) {
			for (iy=0;iy<asy-1;iy++) {
				kh=(iz*(asy-1)+iy)*(asx-1);
				ku=(iz*asy+iy)*asx;
				for (ix=0;ix<asx-1;ix++) {
					fx=Vx/dx*0.25*(u[ku+ix+1]+u[ku+ix+asx+1]+u[ku+ix+1+asxy]+u[ku+ix+asx+1+asxy]-u[ku+ix]-u[ku+ix+asx]-u[ku+ix+asxy]-u[ku+ix+asx+asxy]);
					fy=Vy/dy*0.25*(u[ku+ix+asx]+u[ku+ix+asx+1]+u[ku+ix+asx+asxy]+u[ku+ix+asx+1+asxy]-u[ku+ix]-u[ku+ix+1]-u[ku+ix+asxy]-u[ku+ix+1+asxy]);
					fz=Vz/dz*0.25*(u[ku+ix+asxy]+u[ku+ix+asx+asxy]+u[ku+ix+1+asxy]+u[ku+ix+asx+1+asxy]-u[ku+ix]-u[ku+ix+asx]-u[ku+ix+1]-u[ku+ix+asx+1]);
					u_half[kh+ix]=0.125*(u[ku+ix]+u[ku+ix+1]+u[ku+ix+asx]+u[ku+ix+asx+1]+u[ku+ix+asxy]+u[ku+ix+1+asxy]+u[ku+ix+asx+asxy]+u[ku+ix+asx+1+asxy])-0.5*DT*(fx+fy+fz);
				}
			}
		} 
		// All threads must finish with u_half before writing to u in next loop
		// hence we place a barrier
		#pragma omp barrier
		//A diffusion part will possibly be added at a later stage
		//Now the full step
		for (iz=iStart2;iz<iEnd2;iz++) {
			for (iy=0;iy<lsy;iy++) {
				kh=(iz*(asy-1)+iy)*(asx-1);
				ku=((iz+1)*asy+iy+1)*asx+1;
				for (ix=0;ix<lsx;ix++) {
					fx=0.25*Vx/dx*(u_half[kh+ix+1]+u_half[kh+ix+asx]+u_half[kh+ix+1+asxm1ym1]+u_half[kh+ix+asx+asxm1ym1]-u_half[kh+ix]-u_half[kh+ix+asx-1]-u_half[kh+ix+asxm1ym1]-u_half[kh+ix+asx-1+asxm1ym1]);
					fy=0.25*Vy/dy*(u_half[kh+ix+asx-1]+u_half[kh+ix+asx]+u_half[kh+ix+asx-1+asxm1ym1]+u_half[kh+ix+asx+asxm1ym1]-u_half[kh+ix]-u_half[kh+ix+1]-u_half[kh+ix+asxm1ym1]-u_half[kh+ix+1+asxm1ym1]);
					fz=0.25*Vz/dz*(u_half[kh+ix+asx-1+asxm1ym1]+u_half[kh+ix+asx+asxm1ym1]+u_half[kh+ix+asxm1ym1]+u_half[kh+ix+1+asxm1ym1]-u_half[kh+ix]-u_half[kh+ix+1]-u_half[kh+ix+asx-1]-u_half[kh+ix+asx]);
					u[ku+ix]+=-DT*(fx+fy+fz);
				}
			}
		} 
	}// end omp parallel
	delete u_half;
};
	
	
	
// In this version we wrap the outer two loops into a single loop
// This should improve openmp scalability
// It may also speed up serial code a little??? (need to test)
void compute_LW_scheme2(const double DT,double *u,const int lsx,const int lsy,const int lsz,const double Vx,const double Vy,const double Vz,
const double dx, const double dy, const double dz) {
#if 0
	const int asx=2+lsx;
	const int asy=2+lsy;
	const int asz=2+lsz;
	const int asxy=asx*asy;
	const int asxm1ym1=(asx-1)*(asy-1);
	//Lax Wendroff approach
	double *u_half=new double[(asx-1)*(asy-1)*(asz-1)];
	// We wrap the whole time stepping in a parallel loop to save
	// repeated opening and closing of omp threads.
	// This should reduce overhead and also gets around a bug in Apple's gcc.
	// In this version we use a "coarse-grained" omp approach.
	// It should have less overhead/better scalability
	#pragma omp parallel
	{
		int ix,iz,iyz,ku,kh;
		double fx,fy,fz;
		const int myOMPRank = omp_get_thread_num();
		const int numOMPThreads = omp_get_num_threads();
		const int iStart1 = (myOMPRank*(asz-1)*(asy-1))/numOMPThreads;
		const int iEnd1 = ((myOMPRank+1)*(asz-1)*(asy-1))/numOMPThreads;
		const int iStart2 = (myOMPRank*lsz*lsy)/numOMPThreads;
		const int iEnd2 = ((myOMPRank+1)*lsz*lsy)/numOMPThreads;
			//first the half step
			for (iyz=iStart1;iyz<iEnd1;iyz++) {
				iz=iyz/(asy-1);
				kh=iyz*(asx-1);
				ku=(iyz+iz)*asx;
				for (ix=0;ix<asx-1;ix++) {
					fx=Vx/dx*0.25*(u[ku+ix+1]+u[ku+ix+asx+1]+u[ku+ix+1+asxy]+u[ku+ix+asx+1+asxy]-u[ku+ix]-u[ku+ix+asx]-u[ku+ix+asxy]-u[ku+ix+asx+asxy]);
					fy=Vy/dy*0.25*(u[ku+ix+asx]+u[ku+ix+asx+1]+u[ku+ix+asx+asxy]+u[ku+ix+asx+1+asxy]-u[ku+ix]-u[ku+ix+1]-u[ku+ix+asxy]-u[ku+ix+1+asxy]);
					fz=Vz/dz*0.25*(u[ku+ix+asxy]+u[ku+ix+asx+asxy]+u[ku+ix+1+asxy]+u[ku+ix+asx+1+asxy]-u[ku+ix]-u[ku+ix+asx]-u[ku+ix+1]-u[ku+ix+asx+1]);
					u_half[kh+ix]=0.125*(u[ku+ix]+u[ku+ix+1]+u[ku+ix+asx]+u[ku+ix+asx+1]+u[ku+ix+asxy]+u[ku+ix+1+asxy]+u[ku+ix+asx+asxy]+u[ku+ix+asx+1+asxy])-0.5*DT*(fx+fy+fz);
				}
			}
			// All threads must finish with u_half before writing to u in next loop
			// hence we place a barrier
			#pragma omp barrier
			//A diffusion part will possibly be added at a later stage
			//Now the full step
			for (iyz=iStart2;iyz<iEnd2;iyz++) {
				iz=iyz/lsy;
				kh=(iyz+iz)*(asx-1);
				ku=(iyz+2*iz+asy+1)*asx+1;
				for (ix=0;ix<lsx;ix++) {
					fx=0.25*Vx/dx*(u_half[kh+ix+1]+u_half[kh+ix+asx]+u_half[kh+ix+1+asxm1ym1]+u_half[kh+ix+asx+asxm1ym1]-u_half[kh+ix]-u_half[kh+ix+asx-1]-u_half[kh+ix+asxm1ym1]-u_half[kh+ix+asx-1+asxm1ym1]);
					fy=0.25*Vy/dy*(u_half[kh+ix+asx-1]+u_half[kh+ix+asx]+u_half[kh+ix+asx-1+asxm1ym1]+u_half[kh+ix+asx+asxm1ym1]-u_half[kh+ix]-u_half[kh+ix+1]-u_half[kh+ix+asxm1ym1]-u_half[kh+ix+1+asxm1ym1]);
					fz=0.25*Vz/dz*(u_half[kh+ix+asx-1+asxm1ym1]+u_half[kh+ix+asx+asxm1ym1]+u_half[kh+ix+asxm1ym1]+u_half[kh+ix+1+asxm1ym1]-u_half[kh+ix]-u_half[kh+ix+1]-u_half[kh+ix+asx-1]-u_half[kh+ix+asx]);
					u[ku+ix]+=-DT*(fx+fy+fz);
				}
			}
	}// end omp parallel
	delete u_half;
#else
	const int asx=2+lsx;
	const int asy=2+lsy;
	const int asz=2+lsz;
	double *u_temp=new double[asx*asy*asz];
	int i,j,l,ku;
		//first the x-direction
		for (l=0;l<asz;l++) {
			ku=l*asy*asx;
			for (i=0;i<asy;i++) {
				for (j=1;j<asx-1;j++) {
					u_temp[ku+j]=u[ku+j]-0.5*DT*Vx/dx*(u[ku+j+1]-u[ku+j-1])+0.5*DT*Vx/dx*DT*Vx/dx*(u[ku+j-1]-2.0*u[ku+j]+u[ku+j+1]);
				};
				ku+=asx;
			}
		}
		// now the y-direction
		for (l=0;l<asz;l++) {
			ku=l*asy*asx;
			for (i=1;i<asy-1;i++) {
				ku+=asx;
				for (j=1;j<asx-1;j++) {
					u[ku+j]=u_temp[ku+j]-0.5*DT*Vy/dy*(u_temp[ku+j+asx]-u_temp[ku+j-asx])+0.5*DT*Vy/dy*DT*Vy/dy*(u_temp[ku+j-asx]-2.0*u_temp[ku+j]+u_temp[ku+j+asx]);
				};
			}
		}
		// now the z-direction
		for (l=1;l<asz-1;l++) {
			ku=l*asy*asx;
			for (i=1;i<asy-1;i++) {
				ku+=asx;
				for (j=1;j<asx-1;j++) {
					u_temp[ku+j]=u[ku+j]-0.5*DT*Vz/dz*(u[ku+j+asy*asx]-u[ku+j-asy*asx])+0.5*DT*Vz/dz*DT*Vz/dz*(u[ku+j-asy*asx]-2.0*u[ku+j]+u[ku+j+asy*asx]);
				};
			}
		}
		//Now copy solution back to u (alternatively, could swap the pointers)
		for (l=1;l<asz-1;l++) {
			ku=l*asy*asx;
			for (i=1;i<asy-1;i++) {
				ku+=asx;
				for (j=1;j<asx-1;j++) {
					u[ku+j]=u_temp[ku+j];
				};
			}
		}
	delete u_temp;
	  
#endif
};


void compute_MacCormack_scheme(double DT,double *u,int lsx,int lsy,int lsz,double Vx,double Vy,double Vz, const double dx, const double dy, const double dz) {
	int asx=lsx+2;
	int asy=lsy+2;
	int asz=/*lsz-1*/lsz+2;
	int asxy=asx*asy;
	int asxm1ym1=(asx-1)*(asy-1);
	//MacCormack approach
	double *u_pred=new double[(asx-1)*(asy-1)*(asz-1)];
	int ix,iy,iz,ku,kp;
		//first the predictor step
		for (iz=0;iz<asz-1;iz++) {
			for (iy=0;iy<asy-1;iy++) {
				kp=(iz*(asy-1)+iy)*(asx-1);
				ku=(iz*asy+iy)*asx;
				for (ix=0;ix<asx-1;ix++) {
					u_pred[kp+ix]=u[ku+ix]-Vx*DT/dx*(u[ku+ix+1]-u[ku+ix])-Vy*DT/dy*(u[ku+ix+asx]-u[ku+ix])-Vz*DT/dz*(u[ku+ix+asxy]-u[ku+ix]);
				}
			}
		}
		//A diffusion part will possibly be added at a later stage
		//Now the corrector step
		for (iz=0;iz<asz-1;iz++) {
			for (iy=0;iy<lsy;iy++) {
				kp=(iz*(asy-1)+iy)*(asx-1);
				ku=((iz+1)*asy+iy+1)*asx+1;
				for (ix=0;ix<lsx;ix++) {
					u[ku+ix]=0.5*(u[ku+ix]+u_pred[kp+ix+asx+asxm1ym1])-0.5*Vx*DT/dx*(u_pred[kp+ix+asx+asxm1ym1]-u_pred[kp+ix+asx+asxm1ym1-1])-0.5*Vy*DT/dy*(u_pred[kp+ix+asx+asxm1ym1]-u_pred[kp+ix+asxm1ym1+1])-0.5*Vz*DT/dz*(u_pred[kp+ix+asx+asxm1ym1]-u_pred[kp+ix+asx]);
				}
			}
		}
		//Need to update the boundaries before moving on to the next time step
	delete u_pred;
}

}
