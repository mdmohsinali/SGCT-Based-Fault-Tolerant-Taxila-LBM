ifndef LBM_DIR
  LBM_DIR:=$(shell pwd)/../..
endif

ifndef COMP
  HOST:= $(shell hostname | cut -d"-" -f1)
  ifeq ($(HOST),jaguar)
    COMP:=ftn
  else
    COMP:= $(shell basename $(shell mpif90 -show | cut -d" " -f 1))
  endif
endif

LBM_SRC_DIR= $(LBM_DIR)/src/lbm
LBM_UTIL_DIR= $(LBM_DIR)/src/utils
LBM_INCLUDE_DIR= $(LBM_DIR)/include
LBM_LIB_DIR= $(LBM_DIR)/lib
LBM_MOD_DIR= $(LBM_LIB_DIR)/modules
LBM_LIB_NAME= lbm.a

FFLAGS += -I$(LBM_INCLUDE_DIR) 

# module location
ifeq ($(COMP),gfortran)
  FFLAGS += -J$(LBM_MOD_DIR) -I$(LBM_MOD_DIR)
else ifeq ($(COMP),pgf90)
  FFLAGS += -module $(LBM_MOD_DIR) -I$(LBM_MOD_DIR)
else ifeq ($(COMP),pgf95)
  FFLAGS += -module $(LBM_MOD_DIR) -I$(LBM_MOD_DIR)
else ifeq ($(COMP),pgfortran)
  FFLAGS += -module $(LBM_MOD_DIR) -I$(LBM_MOD_DIR)
else ifeq ($(COMP),ifort) # intel compiler
  FFLAGS += -module $(LBM_MOD_DIR)
endif

# optimization flags
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

allclean: 
	rm -f ./*.mod
	rm -f ./*.o

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
