# list of CPP objects
CPP_OBJS = Advect2D.o \
	Advect3D.o \
	FTAdvect3D.o \
	GridCombine3D.o \
	FailureRecovery.o \
	FTGridCombine3D.o \
	FTthreeDimAdvect.o	

# switch to control TAU profile generation and running on non-ft version of mpi
TAU_PROFILE=no
NON_FT=yes
