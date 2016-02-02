! * SGCT-Based-Fault-Tolerant-Taxila-LBM Code source file.
! * Copyright (c) 2015, Md Mohsin Ali. All rights reserved.
! * Licensed under the terms of the BSD License as described in the LICENSE_FT_CODE file.
! * This comment must be retained in any redistributions of this source file.

! * 
! * File       : LBMWrap.F90
! * Description: contains wrapper routines that are able to communicate from
! *              C side to Fortran side
! * Author     : Mohsin Ali
! * Created    : 16 November 2014
! * Updated    : 
! *
MODULE list_module
       !USE mpi_ext
       use petsc
       use LBM_Options_module
       use LBM_BC_module
       use LBM_Logging_module
       use LBM_module      
       use c_interface_module
#ifndef NON_FT
       USE mpi_ext
#endif
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscdmdef.h"

#include "lbm_definitions.h"

       INTEGER  :: new_eh
       INTEGER  :: i_error
END MODULE list_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE mpiErrorHandler(comm, errorCode)

!       INTEGER                        :: comm
!       INTEGER                        :: errorCode
!       INTEGER                        :: failedGroup
!       INTEGER                        :: ierr
!       INTEGER                        :: rank
!       INTEGER                        :: nprocs

!       INTEGER                        :: tempShrink

!       CALL OMPI_Comm_failure_ack(comm, ierr)
!       CALL OMPI_Comm_failure_get_acked(comm, failedGroup, ierr)

!       CALL MPI_COMM_RANK(comm, rank)
!       CALL MPI_COMM_SIZE(comm, nprocs) 

       !if (errorCode .EQ. MPI_ERR_PROC_FAILED) then
       !   WRITE (*, '("(Fortran Error Handler: MPI_ERROR_PROC_FAILED
       !   Detected.) Process ", i2, " of ", i2)') rank, nprocs
       !else
       !   WRITE (*, '("(Fortran Error Handler: Other Failure Detected.)
       !   Process ", i2, " of ", i2)') rank, nprocs
       !endif

!       call sleep(1)       !1 second

!       CALL MPI_Group_free(failedGroup)

       ! Kill all the processes of this sub-grid
!       CALL MPI_ABORT(comm, errorCode, ierr)

END SUBROUTINE mpiErrorHandler


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
SUBROUTINE c_runLBM(in_dir, proc_x, proc_y, proc_z, f_comm, comp_g, rep_comb, n, &
                    nx, ny, nz, n_comb, real_call) BIND(C, NAME='c_runLBM')

       USE list_module
       USE, INTRINSIC :: iso_c_binding, ONLY : c_char, c_null_char, c_int, c_bool
       implicit none

       ! wrapper variables
       CHARACTER (KIND = c_char, LEN = 1), DIMENSION (128), INTENT (IN) :: in_dir
       INTEGER (c_int)                                                  :: proc_x
       INTEGER (c_int)                                                  :: proc_y
       INTEGER (c_int)                                                  :: proc_z
       INTEGER (c_int)                                                  :: f_comm
       INTEGER (c_int)                                                  :: n
       INTEGER (c_int)                                                  :: nx
       INTEGER (c_int)                                                  :: ny
       INTEGER (c_int)                                                  :: nz
       INTEGER (c_int)                                                  :: n_comb
       CHARACTER (LEN = 128)                                            :: f_par_in_dir
       CHARACTER (LEN = 128)                                            :: f_infile
       LOGICAL (c_bool)                                                 :: comp_g
       LOGICAL (c_bool)                                                 :: rep_comb
       LOGICAL (c_bool)                                                 :: real_call
       INTEGER                                                          :: c_i
#ifdef PETSC_USE_64BIT_INDICES 
       INTEGER (KIND = 8) :: proc_x64, proc_y64, proc_z64, n64, nx64, ny64, nz64, n_comb64
#endif

       ! TaxilaLBM variables
       PetscInt istep
       PetscInt ntimes, npasses
       PetscInt kwrite, kprint
       PetscErrorCode ierr
       character(len=MAXWORDLENGTH) prefix    
       ! Mohsin added this
       PetscInt select_field, TWO_DIMS, THREE_DIMS   

       ! TaxilaLBM subroutines
       external initialize_bcs
       external initialize_bcs_transport
       external initialize_state
       external initialize_state_transport
       external initialize_walls
       type(lbm_type),pointer:: lbm
       type(options_type),pointer:: options

       ! mpi error handler subroutine
       !EXTERNAL :: mpiErrorHandler

       f_par_in_dir = " "

       ! http://stackoverflow.com/questions/8207997/calling-a-fortran-subroutine-from-c/8208960
       loop_par_in_dir: DO c_i = 1, 128
          if(in_dir(c_i) == c_null_char) then
                exit loop_par_in_dir
          else
                f_par_in_dir(c_i:c_i) = in_dir(c_i)
          end if
       end do loop_par_in_dir

       ! disable the following (with trim) if create problem somewhere
       ! related to this
       f_par_in_dir = trim(f_par_in_dir)
       ! create and set mpi error handler for f_comm
       !CALL MPI_COMM_CREATE_ERRHANDLER(mpiErrorHandler, new_eh, i_error)
       !CALL MPI_COMM_SET_ERRHANDLER(f_comm, new_eh, i_error)

       !
       !
       ! call the original Petsc fortran routines of TaxilaLBM
       ! --- setup environment
       !call getarg(1, f_infile) !works only in Fortran MAIN_ program
       f_infile = trim(f_par_in_dir) // trim('/input_data')
       f_infile = trim(f_infile)
       call PetscInitialize(f_infile, ierr) !f_infile is coming from C side
       ! Mohsin updated this
       ! parameters are coming from C side       
#ifdef PETSC_USE_64BIT_INDICES 
       proc_x64 = INT(proc_x, 8)
       proc_y64 = INT(proc_y, 8)
       proc_z64 = INT(proc_z, 8)
       n64 = INT(n, 8)
       nx64 = INT(nx, 8)
       ny64 = INT(ny, 8)
       nz64 = INT(nz, 8)
       n_comb64 = INT(n_comb, 8)
       lbm => LBMCreate(f_comm, proc_x64, proc_y64, proc_z64, comp_g, rep_comb, n64, nx64, ny64, nz64, n_comb64)
#else
       lbm => LBMCreate(f_comm, proc_x, proc_y, proc_z, comp_g, rep_comb, n, nx, ny, nz, n_comb)
#endif
       options => lbm%options
       call LoggerCreate()
       call PetscLogStagePush(logger%stage(INIT_STAGE), ierr)

       ! initialize options and constants
       prefix = ''
       call PetscLogEventBegin(logger%event_init_options,ierr)
       call OptionsSetPrefix(options, prefix)
       call OptionsSetUp(options)
       call PetscLogEventEnd(logger%event_init_options,ierr)
       call LBMSetFromOptions(lbm, options, ierr);CHKERRQ(ierr)
       call LBMSetUp(lbm)
       ! set boundary conditions
       call PetscLogEventBegin(logger%event_init_bcsetup,ierr)
       call BCSetValues(lbm%flow%bc, lbm%flow%distribution, options, initialize_bcs)
       if (associated(lbm%transport)) then
           call BCSetValues(lbm%transport%bc, lbm%transport%distribution, &
               options, initialize_bcs_transport)
       end if
       call PetscLogEventEnd(logger%event_init_bcsetup,ierr)

       ! set initial conditions
       call PetscLogEventBegin(logger%event_init_icsetup,ierr)
       if (options%restart) then
           call LBMInitializeStateRestarted(lbm, options%istep, options%kwrite)
           istep = options%istep
       else if (options%ic_from_file) then
           call LBMInitializeStateFromFile(lbm)
           istep = 0
       else
           if (associated(lbm%transport)) then
               call LBMInitializeState(lbm, initialize_state, initialize_state_transport)
           else
               call LBMInitializeState(lbm, initialize_state)
           end if

           if (options%flow_at_steadystate_hasfile) then
               call LBMLoadSteadyStateFlow(lbm, options%flow_at_steadystate_flow_file)
           end if
           istep=0
       endif
       call PetscLogEventEnd(logger%event_init_icsetup,ierr)

       ! start lbm
       if (lbm%grid%info%rank.eq.0) then
           write(*,*) 'calling lbm from inital step', istep, 'to final step', &
               options%ntimes*options%npasses
       end if

       call LBMInit(lbm, istep, options%supress_ic_output)
       call LBMRun(lbm, istep, options%ntimes*options%npasses)

       ! Mohsin added this
       select_field = 2 ! 1 for fi and 2 for rho (now manually)
       TWO_DIMS = 2
       THREE_DIMS = 3
       if (real_call .and. lbm%flow%ndims.eq.TWO_DIMS) then
          select case(select_field)
          case(1)
             call access_LBM_field_fi_2d(lbm, lbm%flow%distribution%fi_a)
          case(2)
             call access_LBM_field_rho_2d(lbm, lbm%flow%distribution%rho_a)
          end select
       else if (real_call .and. lbm%flow%ndims.eq.THREE_DIMS) then
          select case(select_field)
          case(1)
             call access_LBM_field_fi_3d(lbm, lbm%flow%distribution%fi_a)
          case(2)
             call access_LBM_field_rho_3d(lbm, lbm%flow%distribution%rho_a)
          end select
       end if

       call PetscLogStagePop(ierr)
       call PetscLogStagePush(logger%stage(DESTROY_STAGE), ierr)
       call LBMDestroy(lbm, ierr)
       call PetscLogStagePop(ierr)
       call LoggerDestroy()
       ! release error handler for f_comm
       !CALL MPI_ERRHANDLER_FREE(new_eh, i_error)
       call PetscFinalize(ierr)
       !stop !works only in Fortran MAIN_ program
       !
       !

END SUBROUTINE c_runLBM

