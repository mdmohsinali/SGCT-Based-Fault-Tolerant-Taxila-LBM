! * SGCT-Based-Fault-Tolerant-Taxila-LBM Code source file.
! * Copyright (c) 2015, Md Mohsin Ali. All rights reserved.
! * Licensed under the terms of the BSD License as described in the LICENSE_FT_CODEfile.
! * This comment must be retained in any redistributions of this source file.

! * 
! * File       : LBMWrapPassField.F90
! * Description: contains wrapper routines that are able to communicate from
! *              C side to Fortran side
! * Author     : Mohsin Ali
! * Created    : 16 November 2014
! * Updated    : 
! *
MODULE list_module_2
       use petsc
       use LBM_Options_module
       use LBM_BC_module
       use LBM_Logging_module
       use LBM_module
       use LBM_IO_module
#ifndef NON_FT
       USE mpi_ext
#endif

#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscdmdef.h"

!#include <finclude/petscsys.h>
!#include <finclude/petscvec.h>
!#include <finclude/petscviewer.h>
!#include <finclude/petscviewer.h90>
!#include <finclude/petsclog.h>

#include "lbm_definitions.h"

       INTEGER  :: new_eh
       INTEGER  :: i_error

END MODULE list_module_2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
SUBROUTINE visualizeField(proc_x, proc_y, proc_z, f_comm, &
                    nx, ny, nz, l_dim, f_double_u) BIND(C, NAME='visualizeField')

       USE list_module_2
       USE, INTRINSIC :: iso_c_binding, ONLY : c_int, c_double
       implicit none

       ! wrapper variables
       INTEGER (c_int)                    :: proc_x
       INTEGER (c_int)                    :: proc_y
       INTEGER (c_int)                    :: proc_z
       INTEGER (c_int)                    :: f_comm
       INTEGER (c_int)                    :: nx
       INTEGER (c_int)                    :: ny
       INTEGER (c_int)                    :: nz
       INTEGER (c_int)                    :: l_dim
       REAL (c_double), Dimension(l_dim)  :: f_double_u

      ! Declare variables
#if defined(PETSC_USE_FORTRAN_DATATYPES)
      Type(Vec)  global, local
      Type(DM)   da
#else
      Vec global, local
      DM  da
#endif
     
      PetscInt                M, N, P, ne, nc, dof, s, px, py, pz
      PetscInt, Pointer    :: e(:)
      PetscScalar             value
      PetscInt                i,j, ixs, ixm, iys, iym, arr_index
      PetscScalar, Pointer :: global_ptr(:,:)
      PetscErrorCode ierr
      MPI_Comm comm
      PetscViewer viewer

      call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

      ! Initialize variables
      comm = f_comm
      px = proc_x
      py = proc_y
      pz = proc_z
      M = nx
      N = ny
      P = nz
      dof = 1
      s = 1

      ! Draw will be active until the user click/prompt
      call PetscOptionsSetValue('-draw_pause','-1',ierr)
      ! Creating viewer object
      call PetscViewerDrawOpen(comm,PETSC_NULL_CHARACTER,&
               PETSC_NULL_CHARACTER,0,0,600,600,viewer,ierr)
      ! Creating DMDA object
      call DMDACreate2d(comm,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,&
           M,N,px,py,dof,s,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,da,ierr)
      ! Getting global array
      call DMGetGlobalVector(da,global,ierr)
      ! Getting corner of global array
      call DMDAGetCorners(da, ixs, iys, PETSC_NULL_INTEGER, ixm, iym, PETSC_NULL_INTEGER,ierr)
      ! Getting global array
      call DMDAVecGetArrayF90(da,global,global_ptr,ierr)
      ! Copy input array to global array
      arr_index = 1
      do j=iys, iys+iym-1
         do i=ixs, ixs+ixm-1
            global_ptr(i, j)=f_double_u(arr_index)
            arr_index=arr_index+1
         end do
      end do
      ! Restore global array to global vector
      call DMDAVecRestoreArrayF90(da,global,global_ptr,ierr)
      ! Draw distribution of vector and print in terminal
      !call DMView(da,PETSC_VIEWER_DRAW_WORLD,ierr)
      !call VecView(global,PETSC_VIEWER_STDOUT_WORLD,ierr)

      ! Draw distribution of vector and perform contour plot
      !call DMView(da, viewer, ierr)
      !call VecView(global, viewer, ierr)

      ! Write in binary format to disk
      call PetscViewerCreate(comm, viewer, ierr)
      call PetscViewerSetType(viewer, PETSCVIEWERBINARY, ierr)
      !call PetscViewerBinarySetMPIIO(viewer, ierr)
      call PetscViewerFileSetMode(viewer, FILE_MODE_WRITE, ierr)
      call PetscViewerFileSetName(viewer, 'test_solution/vector.dat', ierr)
      call VecView(global, viewer, ierr)

      ! Clear memory and objects
      call PetscViewerDestroy(viewer, ierr)
      call DMRestoreGlobalVector(da,global,ierr)
      call DMDestroy(da, ierr)

      call PetscFinalize(ierr)

END SUBROUTINE visualizeField

