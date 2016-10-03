ifndef LBM_DIR
  LBM_DIR:=$(shell pwd)/../..
endif

APP_NAME = lbm_raijin_cluster

LBM_OPT = yes

ifndef COMP
  HOST:= $(shell hostname | cut -d"-" -f1)
  ifeq ($(HOST),jaguar)
    COMP:=ftn
  else
    COMP:= $(shell basename $(shell mpif90 -show | cut -d" " -f 1))
  endif
endif

LBM_SRC_DIR= $(LBM_DIR)/src/lbm
LBM_PROBLEM_DIR= $(LBM_DIR)/src/problem_specs
LBM_INCLUDE_DIR= $(LBM_DIR)/include
LBM_LIB_DIR= $(LBM_DIR)/lib
LBM_MOD_DIR= $(LBM_LIB_DIR)/modules
LBM_LIB_NAME= lbm.a
PAPI_DIR= /apps/papi/5.1.1

# include CPP object files in CPP_OBJS variable, and set some 
# switches (TAU_PROFILE, NON_FT)
include $(LBM_SRC_DIR)/cpp_objs_switches.mk

# linking libraries
LIBS = -L$(PETSC_DIR)/lib -lglpk -lmpi_cxx -lstdc++ -lgomp

# libraries and headers related to TAU
TAU_INCPATHS = -I$(PAPI_DIR)/src -I$(PAPI_DIR)/include  \
        -I$(PETSC_DIR)/include
TAU_LIBS = -L$(PETSC_DIR)/x86_64/lib -L$(PAPI_DIR)/lib/ \
        -DPROFILING_ON -DTAU_DOT_H_LESS_HEADERS -ltau-papi-mpi-pdt -lpapi -lbfd -lrt

FFLAGS = -I$(LBM_INCLUDE_DIR) 

# module location
FFLAGS += -I$(LBM_MOD_DIR)

# optimization/debugging flags
ifndef LBM_OPT
  ifeq ($(COMP),gfortran)
    FFLAGS += -fbounds-check
  else ifeq ($(COMP),pgf90)
    FFLAGS += -Mbounds
  else ifeq ($(COMP),pgf95)
    FFLAGS += -Mbounds
  else ifeq ($(COMP),pgfortran)
    FFLAGS += -Mbounds
  else ifeq ($(COMP),ifort)
    FFLAGS += -C # all checks, includes bounds
  endif
else
  FFLAGS += -O3
endif

MAIN_WRAP = LBMWrap.o

ifeq ($(strip $(wildcard ft_initialize_state.F90)),)
  FORTRANOBJS = $(LBM_PROBLEM_DIR)/default_initialize_state.o
else
  FORTRANOBJS = ft_initialize_state.o
endif
ifeq ($(strip $(wildcard initialize_walls.F90)),)
  FORTRANOBJS += $(LBM_PROBLEM_DIR)/default_initialize_walls.o
else
  FORTRANOBJS += initialize_walls.o
endif
ifeq ($(strip $(wildcard initialize_bcs.F90)),)
  FORTRANOBJS += $(LBM_PROBLEM_DIR)/default_initialize_bcs.o
else
  FORTRANOBJS += initialize_bcs.o
endif
ifeq ($(strip $(wildcard initialize_bcs_transport.F90)),)
  FORTRANOBJS += $(LBM_PROBLEM_DIR)/default_initialize_bcs_transport.o
else
  FORTRANOBJS += initialize_bcs_transport.o
endif
ifeq ($(strip $(wildcard initialize_state_transport.F90)),)
  FORTRANOBJS += $(LBM_PROBLEM_DIR)/default_initialize_state_transport.o
else
  FORTRANOBJS += initialize_state_transport.o
endif

# add linking libraries and headers for TAU
ifeq ($(TAU_PROFILE),yes)
  INCPATHS += $(TAU_INCPATHS)
  LIBS += $(TAU_LIBS)
endif

fortran: copy_advectionGridCombine3D_files BEGIN source $(FORTRANOBJS) chkopts
	@echo ""
	@echo "Linking executable."
	@echo "-------------------"
	-${FLINKER} -o $(APP_NAME) $(FORTRANOBJS) \
	$(LBM_SRC_DIR)/$(MAIN_WRAP) \
	$(foreach OBJ, $(CPP_OBJS), $(LBM_SRC_DIR)/$(OBJ)) \
	$(INCPATHS) $(LIBS) \
	$(LBM_LIB_DIR)/$(LBM_LIB_NAME) ${PETSC_DM_LIB}

copy_advectionGridCombine3D_files:
	(pwd;\
	cd ../../../advectionGridCombine3D/;\
	find *.cpp ! -name 'threeDimAdvect.cpp' -print | xargs -I {} cp -p -u {} $(LBM_SRC_DIR);\
	cd $(LBM_DIR)/tests/bubble_3D/;\
	cd ../../../advectionGridCombine3D/;\
	find *.h ! -name 'threeDimAdvect.cpp' -print | xargs -I {} cp -p -u {} $(LBM_SRC_DIR);\
	cd $(LBM_DIR)/tests/bubble_3D/;\
	pwd)
	@echo "All *.cpp (except threeDimAdvect.cpp) and *.h files from \"../../../advectionGridCombine3D\" folder are copied to \"src/lbm\" directory"

BEGIN:
	@echo "Building lbm simulation in $(LBM_DIR) with $(COMP)" 

allclean: 
	@cd $(LBM_SRC_DIR); make allclean
	@cd $(LBM_PROBLEM_DIR); make allclean
	rm -f ./*.o
	rm -f $(APP_NAME)

source: 
	@cd $(LBM_SRC_DIR); make fortran

run2d:
	scp input_data.2d input_data
	mpirun -np 10 ./$(APP_NAME) -v 1 -2 -p 1 -q 1 -r 1 -l 4 7 7 0 #2D version, -F option needs for visualization

run3d:
	scp input_data.3d input_data
	mpirun -np 20 ./$(APP_NAME) -v 1 -p 1 -q 1 -r 1 -l 4 6 6 6 #3D version

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
