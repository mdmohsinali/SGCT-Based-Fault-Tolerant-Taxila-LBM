! * SGCT-Based-Fault-Tolerant-Taxila-LBM Code source file.
! * Copyright (c) 2015, Md Mohsin Ali. All rights reserved.
! * Licensed under the terms of the BSD License as described in the LICENSE_FT_CODE file.
! * This comment must be retained in any redistributions of this source file.

!!!==================================================================
!!! Fortran-file
!!!    author         : Mohsin Ali
!!!    filename       : c_interface.F90
!!!    version        :
!!!    created        : 18 November 2014
!!!    updated        : 17 February 2015
!!!    last modified  : 
!!!    URL            : http://users.cecs.anu.edu.au/~mohsin/
!!!    email          : md _dot_ ali _at_ anu.edu.au
!!!
!!!==================================================================

#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscdmdef.h"

module c_interface_module
    use petsc
    use LBM_Options_module
    use LBM_BC_module
    use LBM_Logging_module
    use LBM_module
    implicit none

    private

#include "lbm_definitions.h"

    type(lbm_type),pointer:: lbm
    PetscErrorCode ierr
    PetscInt  i, j, k, i1, j1, k1, &
              x_start, y_start, z_start, &
              x_width, y_width, z_width
    PetscScalar, Pointer ::  fi_2d(:,:), fi_3d(:,:,:)
    !x_start, y_start, z_start   - starting grid indices, like c (no ghost points)
    !x_width, y_width, z_width   - widths of local grid, like c (no ghost points)

    LOGICAL (KIND = 1) :: f_comp_grid, f_logi_rep_comb
    INTEGER (KIND = 4) :: ndims32, x_start32, y_start32, z_start32, x_width32, y_width32, z_width32

    public :: access_LBM_field_fi_2d, &
              access_LBM_field_rho_2d, &
              access_LBM_field_fi_3d, &
              access_LBM_field_rho_3d
    PetscInt select_component, select_dimension

contains

    subroutine access_LBM_field_fi_2d(lbm, fi)
       INTERFACE
          SUBROUTINE c_get_LBM_field(dimen, f_2d, f_3d, xs, ys, zs, xm, ym, zm, comp_g, rep_comb) BIND (C, NAME = 'c_get_LBM_field')
             USE, INTRINSIC :: iso_c_binding, ONLY : c_int, c_double, c_bool, c_null_char
             LOGICAL (c_bool)                       :: comp_g
             LOGICAL (c_bool)                       :: rep_comb
             INTEGER (c_int)                        :: xs, ys, zs
             INTEGER (c_int)                        :: xm, ym, zm
             INTEGER (c_int)                        :: dimen
             REAL (c_double), Dimension(xm, ym)     :: f_2d
             REAL (c_double), Dimension(xm, ym, zm) :: f_3d
          END SUBROUTINE c_get_LBM_field
       END INTERFACE

       type(lbm_type) lbm

       PetscScalar,dimension(lbm%flow%ncomponents, &
            0:lbm%flow%disc%b, &
            lbm%flow%distribution%info%gxs:lbm%flow%distribution%info%gxe, &
            lbm%flow%distribution%info%gys:lbm%flow%distribution%info%gye) :: fi

       !! User to choose from
       select_component = 1 ! n for selecting nth component out of ncomponents values
       select_dimension = 0 ! n for selecting nth dimension from (0 <= dimensions < 9 for D2Q9)

       call DMDAGetCorners(lbm%flow%grid%da(NCOMPONENTDOF), x_start, y_start, &
            z_start, x_width, y_width, z_width, ierr)

       ! Allocate memory
       ALLOCATE(fi_2d(x_width,y_width))

       ! Copy values return to C
       i1 = lbm%flow%distribution%info%gxs
       j1 = lbm%flow%distribution%info%gys
       do j = 1, y_width
           do i = 1, x_width
               fi_2d(i, j) = fi(select_component, select_dimension, i1+i, j1+j)
           enddo
       enddo

       f_comp_grid = lbm%grid%info%is_comp_grid
       f_logi_rep_comb = lbm%grid%info%is_rep_comb

       ndims32 = INT(lbm%flow%ndims, 4)
       x_start32 = INT(x_start, 4)
       y_start32 = INT(y_start, 4)
       z_start32 = INT(z_start, 4)
       x_width32 = INT(x_width, 4)
       y_width32 = INT(y_width, 4)
       z_width32 = INT(z_width, 4)

       call c_get_LBM_field(ndims32, fi_2d, fi_3d, &
            x_start32, y_start32, z_start32, &
            x_width32, y_width32, z_width32, &
            f_comp_grid, f_logi_rep_comb)

       ! Deallocate memory and nullify
       DEALLOCATE(fi_2d)
       NULLIFY(fi_2d)
       NULLIFY(fi_3d)
    end subroutine access_LBM_field_fi_2d

    subroutine access_LBM_field_rho_2d(lbm, rho)
       INTERFACE
          SUBROUTINE c_get_LBM_field(dimen, f_2d, f_3d, xs, ys, zs, xm, ym, zm, comp_g, rep_comb) BIND (C, NAME = 'c_get_LBM_field')
             USE, INTRINSIC :: iso_c_binding, ONLY : c_int, c_double, c_bool, c_null_char
             LOGICAL (c_bool)                       :: comp_g
             LOGICAL (c_bool)                       :: rep_comb
             INTEGER (c_int)                        :: xs, ys, zs
             INTEGER (c_int)                        :: xm, ym, zm
             INTEGER (c_int)                        :: dimen
             REAL (c_double), Dimension(xm, ym)     :: f_2d
             REAL (c_double), Dimension(xm, ym, zm) :: f_3d
          END SUBROUTINE c_get_LBM_field
       END INTERFACE

       type(lbm_type) lbm

       PetscScalar,dimension(lbm%flow%ncomponents, &
            lbm%flow%distribution%info%rgxs:lbm%flow%distribution%info%rgxe, &
            lbm%flow%distribution%info%rgys:lbm%flow%distribution%info%rgye) :: rho

       !! User to choose from
       select_component = 1 ! n for selecting nth component out of ncomponents values

       call DMDAGetCorners(lbm%flow%grid%da(NCOMPONENTDOF), x_start, y_start, &
            z_start, x_width, y_width, z_width, ierr)

       ! Allocate memory
       ALLOCATE(fi_2d(x_width,y_width))

       ! Copy values return to C
       i1 = lbm%flow%distribution%info%rgxs
       j1 = lbm%flow%distribution%info%rgys
       do j = 1, y_width
           do i = 1, x_width
               fi_2d(i, j) = rho(select_component, i1+i, j1+j)
               !print*, 'Fortran value = ', rho(select_component, i1+i, j1+j) 
           enddo
       enddo

       f_comp_grid = lbm%grid%info%is_comp_grid
       f_logi_rep_comb = lbm%grid%info%is_rep_comb

       ndims32 = INT(lbm%flow%ndims, 4)
       x_start32 = INT(x_start, 4)
       y_start32 = INT(y_start, 4)
       z_start32 = INT(z_start, 4)
       x_width32 = INT(x_width, 4)
       y_width32 = INT(y_width, 4)
       z_width32 = INT(z_width, 4)

       call c_get_LBM_field(ndims32, fi_2d, fi_3d, &
            x_start32, y_start32, z_start32, &
            x_width32, y_width32, z_width32, &
            f_comp_grid, f_logi_rep_comb)

       ! Deallocate memory and nullify
       DEALLOCATE(fi_2d)
       NULLIFY(fi_2d)
       NULLIFY(fi_3d)
    end subroutine access_LBM_field_rho_2d

    subroutine access_LBM_field_fi_3d(lbm, fi)
       INTERFACE
          SUBROUTINE c_get_LBM_field(dimen, f_2d, f_3d, xs, ys, zs, xm, ym, zm, comp_g, rep_comb) BIND (C, NAME = 'c_get_LBM_field')
             USE, INTRINSIC :: iso_c_binding, ONLY : c_int, c_double, c_bool, c_null_char
             LOGICAL (c_bool)                       :: comp_g
             LOGICAL (c_bool)                       :: rep_comb
             INTEGER (c_int)                        :: xs, ys, zs
             INTEGER (c_int)                        :: xm, ym, zm
             INTEGER (c_int)                        :: dimen
             REAL (c_double), Dimension(xm, ym)     :: f_2d
             REAL (c_double), Dimension(xm, ym, zm) :: f_3d
          END SUBROUTINE c_get_LBM_field
       END INTERFACE

       type(lbm_type) lbm

       PetscScalar,dimension(lbm%flow%ncomponents, &
            0:lbm%flow%disc%b, &
            lbm%flow%distribution%info%gxs:lbm%flow%distribution%info%gxe, &
            lbm%flow%distribution%info%gys:lbm%flow%distribution%info%gye, &
            lbm%flow%distribution%info%gzs:lbm%flow%distribution%info%gze) :: fi

       !! User to choose from
       select_component = 1 ! n for selecting nth component out of ncomponents values
       select_dimension = 0 ! n for selecting nth dimension from (0 <= dimensions < 19 for D3Q19)

       call DMDAGetCorners(lbm%flow%grid%da(NCOMPONENTDOF), x_start, y_start, &
            z_start, x_width, y_width, z_width, ierr)

       ! Allocate memory
       ALLOCATE(fi_3d(x_width,y_width,z_width))

       ! Copy values return to C
       i1 = lbm%flow%distribution%info%gxs
       j1 = lbm%flow%distribution%info%gys
       k1 = lbm%flow%distribution%info%gzs
       do k = 1, z_width
           do j = 1, y_width
               do i = 1, x_width
                   fi_3d(i, j, k) = fi(select_component, select_dimension, i1+i, j1+j, k1+k)
               enddo
           enddo
       enddo

       f_comp_grid = lbm%grid%info%is_comp_grid
       f_logi_rep_comb = lbm%grid%info%is_rep_comb

       ndims32 = INT(lbm%flow%ndims, 4)
       x_start32 = INT(x_start, 4)
       y_start32 = INT(y_start, 4)
       z_start32 = INT(z_start, 4)
       x_width32 = INT(x_width, 4)
       y_width32 = INT(y_width, 4)
       z_width32 = INT(z_width, 4)

       call c_get_LBM_field(ndims32, fi_2d, fi_3d, &
            x_start32, y_start32, z_start32, &
            x_width32, y_width32, z_width32, &
            f_comp_grid, f_logi_rep_comb)

       ! Deallocate memory and nullify
       DEALLOCATE(fi_3d)
       NULLIFY(fi_3d)
       NULLIFY(fi_2d)
    end subroutine access_LBM_field_fi_3d

    subroutine access_LBM_field_rho_3d(lbm, rho)
       INTERFACE
          SUBROUTINE c_get_LBM_field(dimen, f_2d, f_3d, xs, ys, zs, xm, ym, zm, comp_g, rep_comb) BIND (C, NAME = 'c_get_LBM_field')
             USE, INTRINSIC :: iso_c_binding, ONLY : c_int, c_double, c_bool, c_null_char
             LOGICAL (c_bool)                       :: comp_g
             LOGICAL (c_bool)                       :: rep_comb
             INTEGER (c_int)                        :: xs, ys, zs
             INTEGER (c_int)                        :: xm, ym, zm
             INTEGER (c_int)                        :: dimen
             REAL (c_double), Dimension(xm, ym)     :: f_2d
             REAL (c_double), Dimension(xm, ym, zm) :: f_3d
          END SUBROUTINE c_get_LBM_field
       END INTERFACE

       type(lbm_type) lbm

       PetscScalar,dimension(lbm%flow%ncomponents, &
            lbm%flow%distribution%info%rgxs:lbm%flow%distribution%info%rgxe, &
            lbm%flow%distribution%info%rgys:lbm%flow%distribution%info%rgye, &
            lbm%flow%distribution%info%rgzs:lbm%flow%distribution%info%rgze) :: rho

       !! User to choose from
       select_component = 1 ! n for selecting nth component out of ncomponents values

       call DMDAGetCorners(lbm%flow%grid%da(NCOMPONENTDOF), x_start, y_start, &
            z_start, x_width, y_width, z_width, ierr)

       ! Allocate memory
       ALLOCATE(fi_3d(x_width,y_width,z_width))

       ! Copy values return to C
       i1 = lbm%flow%distribution%info%rgxs
       j1 = lbm%flow%distribution%info%rgys
       k1 = lbm%flow%distribution%info%rgzs
       do k = 1, z_width
           do j = 1, y_width
               do i = 1, x_width
                   fi_3d(i, j, k) = rho(select_component, i1+i, j1+j, k1+k)
               enddo
           enddo
       enddo

       f_comp_grid = lbm%grid%info%is_comp_grid
       f_logi_rep_comb = lbm%grid%info%is_rep_comb

       ndims32 = INT(lbm%flow%ndims, 4)
       x_start32 = INT(x_start, 4)
       y_start32 = INT(y_start, 4)
       z_start32 = INT(z_start, 4)
       x_width32 = INT(x_width, 4)
       y_width32 = INT(y_width, 4)
       z_width32 = INT(z_width, 4)

       call c_get_LBM_field(ndims32, fi_2d, fi_3d, &
            x_start32, y_start32, z_start32, &
            x_width32, y_width32, z_width32, &
            f_comp_grid, f_logi_rep_comb)

       ! Deallocate memory and nullify
       DEALLOCATE(fi_3d)
       NULLIFY(fi_3d)
       NULLIFY(fi_2d)
    end subroutine access_LBM_field_rho_3d

end module c_interface_module

