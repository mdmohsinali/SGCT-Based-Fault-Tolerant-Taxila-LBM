!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_distribution_function.F90
!!!     version:         
!!!     created:         28 March 2011
!!!       on:            14:06:07 MDT
!!!     last modified:   07 November 2011
!!!       at:            11:35:30 MST
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!     Licensed under the terms of the BSD License as described in the LICENSE_FT_CODE, LICENSE_TAXILA and COPYRIGHT_TAXILA files.
!!!  
!!! ====================================================================

#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscdmdef.h"

module LBM_Distribution_Function_module
  use petsc
  use LBM_Error_module
  use LBM_Distribution_Function_type_module
  use LBM_Info_module
  use LBM_Discretization_Type_module
  implicit none

  private
#include "lbm_definitions.h"

  interface DistributionCalcFlux
    module procedure DistributionCalcFlux1
    module procedure DistributionCalcFlux2
  end interface

  interface DistributionCalcDensity
    module procedure DistributionCalcDensity1
    module procedure DistributionCalcDensity2
  end interface

  public:: DistributionCreate, &
       DistributionDestroy, &
       DistributionSetName, &
       DistributionSetInfo, &
       DistributionSetDiscretization, &
       DistributionSetSizes, &
       DistributionSetDAs, &
       DistributionSetTrackOld, &
       DistributionSetUp, &
       DistributionGetArrays, &
       DistributionRestoreArrays, &
       DistributionCommunicateAll, &
       DistributionCommunicateFi, &
       DistributionCommunicateFiBegin, &
       DistributionCommunicateFiEnd, &
       DistributionCommunicateDensity, &
       DistributionCommunicateDensityBegin, &
       DistributionCommunicateDensityEnd, &
       DistributionCalcDensity, &
       DistributionCalcFlux, &
       DistributionStream, &
       DistributionBounceback, &
       DistributionGatherValueToDirection, &
       DistributionCalcDeltaNorm

contains
  function DistributionCreate(comm) result(distribution)
    MPI_Comm comm
    type(distribution_type),pointer:: distribution
    allocate(distribution)
    distribution%comm = comm
    distribution%name = ''
    distribution%s = -1
    distribution%b = -1
    nullify(distribution%info)
    nullify(distribution%disc)
    nullify(distribution%da_fi)
    nullify(distribution%da_rho)

    distribution%fi = 0
    distribution%rho = 0
    distribution%fi_g = 0
    distribution%rho_g = 0
    distribution%fi_g_old = 0
    distribution%rho_g_old = 0
    nullify(distribution%fi_a)
    nullify(distribution%rho_a)
    nullify(distribution%flux)
    distribution%flux_required = PETSC_TRUE
    distribution%track_old_fi = PETSC_FALSE
    distribution%track_old_rho = PETSC_FALSE
  end function DistributionCreate

  subroutine DistributionDestroy(distribution, ierr)
    type(distribution_type) distribution
    PetscErrorCode ierr
    if (distribution%fi /= 0) call VecDestroy(distribution%fi,ierr)
    if (distribution%rho /= 0) call VecDestroy(distribution%rho,ierr)
    if (distribution%fi_g /= 0) call VecDestroy(distribution%fi_g,ierr)
    if (distribution%rho_g /= 0) call VecDestroy(distribution%rho_g,ierr)
    if (associated(distribution%flux)) deallocate(distribution%flux)
  end subroutine DistributionDestroy

  subroutine DistributionSetName(distribution, name)
    type(distribution_type) distribution
    character(len=MAXWORDLENGTH):: name
    distribution%name = name
  end subroutine DistributionSetName

  subroutine DistributionSetInfo(distribution, info)
    type(distribution_type) distribution
    type(info_type),pointer:: info
    distribution%info => info
  end subroutine DistributionSetInfo

  subroutine DistributionSetDiscretization(distribution, disc)
    type(distribution_type) distribution
    type(discretization_type),pointer:: disc
    distribution%disc => disc
    distribution%b = disc%b
  end subroutine DistributionSetDiscretization

  subroutine DistributionSetSizes(distribution, s)
    type(distribution_type) distribution
    PetscInt s
    distribution%s = s
  end subroutine DistributionSetSizes

  subroutine DistributionSetDAs(distribution, da_fi, da_rho)
    type(distribution_type) distribution
    DM,target:: da_fi, da_rho
    distribution%da_fi => da_fi
    distribution%da_rho => da_rho
  end subroutine DistributionSetDAs

  subroutine DistributionSetTrackOld(distribution, track_old_rho, track_old_fi)
    type(distribution_type) distribution
    PetscBool track_old_fi, track_old_rho
    distribution%track_old_fi = track_old_fi
    distribution%track_old_rho = track_old_rho
  end subroutine DistributionSetTrackOld

  subroutine DistributionSetUp(distribution)
    type(distribution_type) distribution
    PetscErrorCode ierr
    PetscScalar zero
    zero = 0.

    if (distribution%flux_required) then
       allocate(distribution%flux(1:distribution%s, 1:distribution%info%ndims, &
            1:distribution%info%gxyzl))
       distribution%flux = zero
    end if

    call DMCreateLocalVector(distribution%da_fi, distribution%fi, ierr)
    call VecSet(distribution%fi, zero, ierr)

    call DMCreateGlobalVector(distribution%da_fi, distribution%fi_g, ierr)
    call VecSet(distribution%fi_g, zero, ierr)
    call PetscObjectSetName(distribution%fi_g, trim(distribution%name)//'fi', ierr)
    call VecDuplicate(distribution%fi_g, distribution%fi_g_old, ierr)
    call VecSet(distribution%fi_g_old, zero, ierr)
    call PetscObjectSetName(distribution%fi_g_old, trim(distribution%name)//'fi_old', ierr)

    call DMCreateLocalVector(distribution%da_rho, distribution%rho, ierr)
    call VecSet(distribution%rho, zero, ierr)

    call DMCreateGlobalVector(distribution%da_rho, distribution%rho_g, ierr)
    call VecSet(distribution%rho_g, zero, ierr)
    call PetscObjectSetName(distribution%rho_g, trim(distribution%name)//'rho', ierr)
    call VecDuplicate(distribution%rho_g, distribution%rho_g_old, ierr)
    call VecSet(distribution%rho_g_old, zero, ierr)
    call PetscObjectSetName(distribution%rho_g_old, trim(distribution%name)//'rho_old', ierr)
  end subroutine DistributionSetUp

  subroutine DistributionGetArrays(distribution, ierr)
    type(distribution_type) distribution
    PetscErrorCode ierr
    call DMDAVecGetArrayF90(distribution%da_rho, distribution%rho, distribution%rho_a, ierr)
    call DMDAVecGetArrayF90(distribution%da_fi, distribution%fi, distribution%fi_a, ierr)
  end subroutine DistributionGetArrays

  subroutine DistributionRestoreArrays(distribution, ierr)
    type(distribution_type) distribution
    PetscErrorCode ierr
    call DMDAVecRestoreArrayF90(distribution%da_rho,distribution%rho,distribution%rho_a,ierr)
    call DMDAVecRestoreArrayF90(distribution%da_fi,distribution%fi,distribution%fi_a,ierr)
  end subroutine DistributionRestoreArrays

  subroutine DistributionCommunicateAll(distribution)
    type(distribution_type) distribution
    PetscErrorCode ierr
    call DistributionRestoreArrays(distribution, ierr)
    call DMLocalToLocalBegin(distribution%da_fi, distribution%fi, INSERT_VALUES, &
         distribution%fi, ierr)
    call DMLocalToLocalEnd(distribution%da_fi, distribution%fi, INSERT_VALUES, &
         distribution%fi, ierr)
    call DMLocalToLocalBegin(distribution%da_rho, distribution%rho, INSERT_VALUES, &
         distribution%rho, ierr)
    call DMLocalToLocalEnd(distribution%da_rho, distribution%rho, INSERT_VALUES, &
         distribution%rho, ierr)
    call DistributionGetArrays(distribution, ierr)
  end subroutine DistributionCommunicateAll

  subroutine DistributionCommunicateFi(distribution)
    type(distribution_type) distribution
    PetscErrorCode ierr
    call DistributionCommunicateFiBegin(distribution)
    call DistributionCommunicateFiEnd(distribution)
  end subroutine DistributionCommunicateFi

  subroutine DistributionCommunicateFiBegin(distribution)
    type(distribution_type) distribution
    PetscErrorCode ierr

    call DMDAVecRestoreArrayF90(distribution%da_fi, distribution%fi, &
         distribution%fi_a, ierr)
    call DMLocalToLocalBegin(distribution%da_fi, distribution%fi, INSERT_VALUES, &
         distribution%fi, ierr)

  end subroutine DistributionCommunicateFiBegin

  subroutine DistributionCommunicateFiEnd(distribution)
    type(distribution_type) distribution
    PetscErrorCode ierr

    call DMLocalToLocalEnd(distribution%da_fi, distribution%fi, INSERT_VALUES, &
         distribution%fi, ierr)
    call DMDAVecGetArrayF90(distribution%da_fi, distribution%fi, distribution%fi_a, ierr)
  end subroutine DistributionCommunicateFiEnd

  subroutine DistributionCommunicateDensity(distribution)
    type(distribution_type) distribution
    call DistributionCommunicateDensityBegin(distribution)
    call DistributionCommunicateDensityEnd(distribution)
  end subroutine DistributionCommunicateDensity

  subroutine DistributionCommunicateDensityBegin(distribution)
    type(distribution_type) distribution
    PetscErrorCode ierr
    call DMDAVecRestoreArrayF90(distribution%da_rho, distribution%rho, &
         distribution%rho_a, ierr)
    call DMLocalToLocalBegin(distribution%da_rho, distribution%rho, INSERT_VALUES, &
         distribution%rho, ierr)
  end subroutine DistributionCommunicateDensityBegin

  subroutine DistributionCommunicateDensityEnd(distribution)
    type(distribution_type) distribution
    PetscErrorCode ierr
    call DMLocalToLocalEnd(distribution%da_rho, distribution%rho, INSERT_VALUES, &
         distribution%rho, ierr)
    call DMDAVecGetArrayF90(distribution%da_rho, distribution%rho, distribution%rho_a, ierr)
  end subroutine DistributionCommunicateDensityEnd

  subroutine DistributionCalcDensity1(distribution, walls)
    type(distribution_type) distribution
    PetscScalar,dimension(distribution%info%rgxyzl):: walls

    call DistributionCalcDensity2(distribution,walls,distribution%rho_a)
  end subroutine DistributionCalcDensity1

  subroutine DistributionCalcDensity2(distribution, walls, rho)
    type(distribution_type) distribution
    PetscScalar,dimension(1:distribution%info%rgxyzl):: walls
    PetscScalar,dimension(distribution%s,distribution%info%rgxyzl):: rho

    select case(distribution%info%ndims)
    case(2)
       call DistributionCalcDensityD2(distribution, distribution%fi_a, walls, rho)
    case(3)
       call DistributionCalcDensityD3(distribution, distribution%fi_a, walls, rho)
    end select
  end  subroutine DistributionCalcDensity2

  subroutine DistributionCalcDensityD2(distribution, fi, walls, rho)
    type(distribution_type) distribution
    PetscScalar,dimension(1:distribution%s, 0:distribution%b, &
         distribution%info%gxs:distribution%info%gxe, &
         distribution%info%gys:distribution%info%gye):: fi
    PetscScalar,dimension(distribution%info%rgxs:distribution%info%rgxe, &
         distribution%info%rgys:distribution%info%rgye):: walls
    PetscScalar,dimension(distribution%s, &
         distribution%info%rgxs:distribution%info%rgxe, &
         distribution%info%rgys:distribution%info%rgye):: rho

    PetscInt i,j,m
    do j=distribution%info%ys,distribution%info%ye
    do i=distribution%info%xs,distribution%info%xe
       if (walls(i,j).eq.0) then
          rho(:,i,j) = sum(fi(:,:,i,j),2)
       else
          rho(:,i,j) = 0.
       end if
    end do
    end do
  end  subroutine DistributionCalcDensityD2

  subroutine DistributionCalcDensityD3(distribution, fi, walls, rho)
    type(distribution_type) distribution
    PetscScalar,dimension(1:distribution%s, 0:distribution%b, &
         distribution%info%gxs:distribution%info%gxe, &
         distribution%info%gys:distribution%info%gye, &
         distribution%info%gzs:distribution%info%gze):: fi
    PetscScalar,dimension(distribution%info%rgxs:distribution%info%rgxe, &
         distribution%info%rgys:distribution%info%rgye, &
         distribution%info%rgzs:distribution%info%rgze):: walls
    PetscScalar,dimension(distribution%s, &
         distribution%info%rgxs:distribution%info%rgxe, &
         distribution%info%rgys:distribution%info%rgye, &
         distribution%info%rgzs:distribution%info%rgze):: rho

    PetscInt i,j,k,m
    do k=distribution%info%zs,distribution%info%ze
    do j=distribution%info%ys,distribution%info%ye
    do i=distribution%info%xs,distribution%info%xe
       if (walls(i,j,k).eq.0) then
          rho(:,i,j,k) = sum(fi(:,:,i,j,k),2)
       else
          rho(:,i,j,k) = 0.
       end if
    end do
    end do
    end do
  end  subroutine DistributionCalcDensityD3

  subroutine DistributionCalcFlux1(distribution, walls)
    type(distribution_type) distribution
    PetscScalar,dimension(distribution%info%rgxyzl):: walls

    call DistributionCalcFlux2(distribution,walls,distribution%flux)
  end subroutine DistributionCalcFlux1

  subroutine DistributionCalcFlux2(distribution, walls, flux)
    type(distribution_type) distribution
    PetscScalar,dimension(distribution%info%rgxyzl):: walls
    PetscScalar,dimension(distribution%s, distribution%info%ndims, &
         distribution%info%gxyzl):: flux

    select case(distribution%info%ndims)
    case(2)
       call DistributionCalcFluxD2(distribution,distribution%fi_a,walls,flux)
    case(3)
       call DistributionCalcFluxD3(distribution,distribution%fi_a,walls,flux)
    end select
  end  subroutine DistributionCalcFlux2

  subroutine DistributionCalcFluxD3(distribution, fi, walls, u)
    type(distribution_type) distribution
    PetscScalar,dimension(distribution%s, 0:distribution%b, &
         distribution%info%gxs:distribution%info%gxe, &
         distribution%info%gys:distribution%info%gye, &
         distribution%info%gzs:distribution%info%gze):: fi
    PetscScalar,dimension(distribution%info%rgxs:distribution%info%rgxe, &
         distribution%info%rgys:distribution%info%rgye, &
         distribution%info%rgzs:distribution%info%rgze):: walls
    PetscScalar,dimension(distribution%s, distribution%info%ndims, &
         distribution%info%gxs:distribution%info%gxe, &
         distribution%info%gys:distribution%info%gye, &
         distribution%info%gzs:distribution%info%gze):: u

    PetscInt i,j,k,d,m
    do k=distribution%info%zs,distribution%info%ze
    do j=distribution%info%ys,distribution%info%ye
    do i=distribution%info%xs,distribution%info%xe
       if (walls(i,j,k).eq.0) then
          do m=1,distribution%s
          do d=1,distribution%info%ndims
             u(m,d,i,j,k) = sum(fi(m,:,i,j,k)*distribution%disc%ci(:,d),1)
          end do
          end do
       else
          u(:,:,i,j,k) = 0.
       end if
    end do
    end do
    end do
  end subroutine DistributionCalcFluxD3

  subroutine DistributionCalcFluxD2(distribution, fi, walls, u)
    type(distribution_type) distribution
    PetscScalar,dimension(distribution%s, 0:distribution%b, &
         distribution%info%gxs:distribution%info%gxe, &
         distribution%info%gys:distribution%info%gye):: fi
    PetscScalar,dimension(distribution%info%rgxs:distribution%info%rgxe, &
         distribution%info%rgys:distribution%info%rgye):: walls
    PetscScalar,dimension(distribution%s, distribution%info%ndims, &
         distribution%info%gxs:distribution%info%gxe, &
         distribution%info%gys:distribution%info%gye):: u

    PetscInt i,j,d,m
    do j=distribution%info%ys,distribution%info%ye
    do i=distribution%info%xs,distribution%info%xe
       if (walls(i,j).eq.0) then
          do m=1,distribution%s
          do d=1,distribution%info%ndims
             u(m,d,i,j) = sum(dble(distribution%disc%ci(:,d))*fi(m,:,i,j),1)
          end do
          end do
       else
          u(:,:,i,j) = 0.
       end if
    end do
    end do
  end subroutine DistributionCalcFluxD2

  subroutine DistributionGatherValueToDirection(distribution, val, out)
    type(distribution_type) distribution
    PetscScalar,intent(in),dimension(1:distribution%info%rgxyzl):: val
    PetscScalar,intent(out),dimension(0:distribution%b, 1:distribution%info%rgxyzl):: out
    PetscErrorCode ierr
    ! Mohsin added this
    PetscInt code_val

    if (distribution%info%ndims.eq.2) then
       call DistributionGatherValueToDirectionD2(distribution, val, out)
    else if (distribution%info%ndims.eq.3) then
       call DistributionGatherValueToDirectionD3(distribution, val, out)
    else
       ! Mohsin updated the followings
       code_val = 1
       call LBMError(PETSC_COMM_SELF, code_val, 'invalid ndims in LBM', ierr)
       !call LBMError(PETSC_COMM_SELF, 1, 'invalid ndims in LBM', ierr)
    end if
  end subroutine DistributionGatherValueToDirection

  subroutine DistributionGatherValueToDirectionD2(distribution, val, out)
    type(distribution_type) distribution
    PetscScalar,intent(in),dimension(distribution%info%rgxs:distribution%info%rgxe, &
         distribution%info%rgys:distribution%info%rgye):: val
    PetscScalar,intent(out),dimension(0:distribution%b, &
         distribution%info%rgxs:distribution%info%rgxe, &
         distribution%info%rgys:distribution%info%rgye):: out

    PetscInt n
    do n=0,distribution%b
       out(n,distribution%info%xs:distribution%info%xe, &
            distribution%info%ys:distribution%info%ye) = &
              val(distribution%info%xs + distribution%disc%ci(n,X_DIRECTION): &
                  distribution%info%xe + distribution%disc%ci(n,X_DIRECTION), &
                  distribution%info%ys + distribution%disc%ci(n,Y_DIRECTION): &
                  distribution%info%ye + distribution%disc%ci(n,Y_DIRECTION))
    end do
  end subroutine DistributionGatherValueToDirectionD2

  ! --- stream the fi
  subroutine DistributionStream(dist)
    type(distribution_type) dist
    PetscErrorCode ierr
    ! Mohsin added this
    PetscInt code_val

    select case(dist%info%ndims)
    case(3)
      call DistributionStreamD3(dist%fi_a, dist)
    case(2)
      call DistributionStreamD2(dist%fi_a, dist)
    case DEFAULT
       ! Mohsin updated the followings
       code_val = 1
       call LBMError(PETSC_COMM_SELF, code_val, 'invalid discretization in LBM', ierr)
       !call LBMError(PETSC_COMM_SELF, 1, 'invalid discretization in LBM', ierr)
    end select
  end subroutine DistributionStream

  ! --- stream each specie individually
  subroutine DistributionStreamD3(fi, dist)
    ! input variables
    type(distribution_type) dist
    PetscScalar,dimension(dist%s, 0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: fi

    ! local
    PetscScalar,dimension(dist%s, 0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: tmp
    PetscInt n
    PetscInt xs,xe,ys,ye,zs,ze          ! destination indices
    PetscInt sxs,sxe,sys,sye,szs,sze    ! source indices
    type(info_type),pointer:: info
    info => dist%info

    do n=0,dist%b
      xs = info%xs; sxs = info%xs
      xe = info%xe; sxe = info%xe
      ys = info%ys; sys = info%ys
      ye = info%ye; sye = info%ye
      zs = info%zs; szs = info%zs
      ze = info%ze; sze = info%ze

      if (dist%disc%ci(n,X_DIRECTION) < 0) then
        xs = info%xs + dist%disc%ci(n,X_DIRECTION)
        sxe = info%xe - dist%disc%ci(n,X_DIRECTION)
      else if (dist%disc%ci(n,X_DIRECTION) > 0) then
        xe = info%xe + dist%disc%ci(n,X_DIRECTION)
        sxs = info%xs - dist%disc%ci(n,X_DIRECTION)
      end if

      if (dist%disc%ci(n,Y_DIRECTION) < 0) then
        ys = info%ys + dist%disc%ci(n,Y_DIRECTION)
        sye = info%ye - dist%disc%ci(n,Y_DIRECTION)
      else if (dist%disc%ci(n,Y_DIRECTION) > 0) then
        ye = info%ye + dist%disc%ci(n,Y_DIRECTION)
        sys = info%ys - dist%disc%ci(n,Y_DIRECTION)
      end if

      if (dist%disc%ci(n,Z_DIRECTION) < 0) then
        zs = info%zs + dist%disc%ci(n,Z_DIRECTION)
        sze = info%ze - dist%disc%ci(n,Z_DIRECTION)
      else if (dist%disc%ci(n,Z_DIRECTION) > 0) then
        ze = info%ze + dist%disc%ci(n,Z_DIRECTION)
        szs = info%zs - dist%disc%ci(n,Z_DIRECTION)
      end if

      tmp(:,n,xs:xe,ys:ye,zs:ze) = fi(:,n,sxs:sxe,sys:sye,szs:sze)
    end do
    fi = tmp
  end subroutine DistributionStreamD3

  subroutine DistributionStreamD2(fi, dist)
    ! input variables
    type(distribution_type) dist
    PetscScalar,dimension(dist%s, 0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: fi

    ! local
    PetscScalar,dimension(dist%s, 0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: tmp
    PetscInt n
    PetscInt xs,xe,ys,ye
    PetscInt sxs,sxe,sys,sye
    type(info_type),pointer:: info
    info => dist%info

    do n=0,dist%b
      xs = info%xs; sxs = info%xs
      xe = info%xe; sxe = info%xe
      ys = info%ys; sys = info%ys
      ye = info%ye; sye = info%ye

      if (dist%disc%ci(n,X_DIRECTION) < 0) then
        xs = info%xs + dist%disc%ci(n,X_DIRECTION)
        sxe = info%xe - dist%disc%ci(n,X_DIRECTION)
      else if (dist%disc%ci(n,X_DIRECTION) > 0) then
        xe = info%xe + dist%disc%ci(n,X_DIRECTION)
        sxs = info%xs - dist%disc%ci(n,X_DIRECTION)
      end if

      if (dist%disc%ci(n,Y_DIRECTION) < 0) then
        ys = info%ys + dist%disc%ci(n,Y_DIRECTION)
        sye = info%ye - dist%disc%ci(n,Y_DIRECTION)
      else if (dist%disc%ci(n,Y_DIRECTION) > 0) then
        ye = info%ye + dist%disc%ci(n,Y_DIRECTION)
        sys = info%ys - dist%disc%ci(n,Y_DIRECTION)
      end if

      tmp(:,n,xs:xe,ys:ye) = fi(:,n,sxs:sxe,sys:sye)
    end do
    fi = tmp
  end subroutine DistributionStreamD2

  subroutine DistributionBounceback(dist, walls)
    type(distribution_type) dist
    PetscScalar,dimension(dist%info%rgxyzl):: walls
    PetscErrorCode ierr
    ! Mohsin added this
    PetscInt code_val

    select case(dist%info%ndims)
    case(3)
       call DistributionBouncebackD3(dist, dist%fi_a, walls)
    case(2)
       call DistributionBouncebackD2(dist, dist%fi_a, walls)
    case DEFAULT
       ! Mohsin updated the followings
       code_val = 1
       call LBMError(PETSC_COMM_SELF, code_val, 'invalid discretization in LBM', ierr)
       !call LBMError(PETSC_COMM_SELF, 1, 'invalid discretization in LBM', ierr)
    end select
  end subroutine DistributionBounceback

  subroutine DistributionBouncebackD3(dist, fi, walls)
    type(distribution_type) dist
    PetscScalar,dimension(dist%s,0:dist%b, &
         dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, &
         dist%info%gzs:dist%info%gze):: fi
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye, &
         dist%info%rgzs:dist%info%rgze):: walls

    PetscInt i,j,k,n
    PetscInt tmp(dist%s, 0:dist%b)
    PetscInt ni,nj,nk,nn ! new, bounced back values

    ! loop over ghost cells, turning around and streaming back
    do k=dist%info%gzs,dist%info%gze
    do j=dist%info%gys,dist%info%gye
    do i=dist%info%gxs,dist%info%gxe
      if (walls(i,j,k).eq.WALL_NORMAL_X) then
        do n=0,dist%b
          ni = i - dist%disc%ci(n,X_DIRECTION)
          if (ni <= dist%info%xe .and. ni >= dist%info%xs) then
            nj = j; nk = k
            nn = dist%disc%reflect_x(n)
            fi(:,nn,ni,nj,nk) = fi(:,n,i,j,k)
          endif
        enddo
        fi(:,:,i,j,k) = 0.d0
      else if (walls(i,j,k).eq.WALL_NORMAL_Y) then
        do n=0,dist%b
          nj = j - dist%disc%ci(n,Y_DIRECTION)
          if (nj <= dist%info%ye .and. nj >= dist%info%ys) then
            ni = i; nk = k
            nn = dist%disc%reflect_y(n)
            fi(:,nn,ni,nj,nk) = fi(:,n,i,j,k)
          endif
        enddo
        fi(:,:,i,j,k) = 0.d0
      else if (walls(i,j,k).eq.WALL_NORMAL_Z) then
        do n=0,dist%b
          nk = k - dist%disc%ci(n,Z_DIRECTION)
          if (nk <= dist%info%ze .and. nk >= dist%info%zs) then
            ni = i; nj = j
            nn = dist%disc%reflect_z(n)
            fi(:,nn,ni,nj,nk) = fi(:,n,i,j,k)
          endif
        enddo
        fi(:,:,i,j,k) = 0.d0
      else if (walls(i,j,k) > 0) then
        do n=0,dist%b
          ni = i - dist%disc%ci(n,X_DIRECTION)
          nj = j - dist%disc%ci(n,Y_DIRECTION)
          nk = k - dist%disc%ci(n,Z_DIRECTION)
          if ((ni <= dist%info%xe).and.(ni >= dist%info%xs).and. &
               (nj <= dist%info%ye).and.(nj >= dist%info%ys).and. &
               (nk <= dist%info%ze).and.(nk >= dist%info%zs)) then
            nn = dist%disc%opposites(n)
            fi(:,nn,ni,nj,nk) = fi(:,n,i,j,k)
          end if
        end do
        fi(:,:,i,j,k) = 0.d0
      end if
    end do
    end do
    end do
  end subroutine DistributionBouncebackD3

  subroutine DistributionBouncebackD2(dist, fi, walls)
    type(distribution_type) dist
    PetscScalar,dimension(dist%s,0:dist%b,&
         dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: fi
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: walls

    PetscInt i,j,n
    PetscInt tmp(dist%s, 0:dist%b)
    PetscInt ni,nj,nn ! new, bounced back values

    do j=dist%info%gys,dist%info%gye
    do i=dist%info%gxs,dist%info%gxe
       if (walls(i,j).eq.WALL_NORMAL_X) then
        do n=0,dist%b
          ni = i - dist%disc%ci(n,X_DIRECTION)
          if (ni <= dist%info%xe .and. ni >= dist%info%xs) then
            nj = j
            nn = dist%disc%reflect_x(n)
            fi(:,nn,ni,nj) = fi(:,n,i,j)
          endif
        enddo
        fi(:,:,i,j) = 0.d0
      else if (walls(i,j).eq.WALL_NORMAL_Y) then
        do n=0,dist%b
          nj = j - dist%disc%ci(n,Y_DIRECTION)
          if (nj <= dist%info%ye .and. nj >= dist%info%ys) then
            ni = i
            nn = dist%disc%reflect_y(n)
            fi(:,nn,ni,nj) = fi(:,n,i,j)
          endif
        enddo
        fi(:,:,i,j) = 0.d0
      else if (walls(i,j) > 0) then
        do n=0,dist%b
          ni = i - dist%disc%ci(n,X_DIRECTION)
          nj = j - dist%disc%ci(n,Y_DIRECTION)
          if ((ni <= dist%info%xe).and.(ni >= dist%info%xs).and. &
               (nj <= dist%info%ye).and.(nj >= dist%info%ys)) then
            nn = dist%disc%opposites(n)
            fi(:,nn,ni,nj) = fi(:,n,i,j)
          end if
        end do
        fi(:,:,i,j) = 0.d0
      end if
    end do
    end do
  end subroutine DistributionBouncebackD2

  subroutine DistributionGatherValueToDirectionD3(distribution, val, out)
    type(distribution_type) distribution
    PetscScalar,intent(in),dimension(distribution%info%rgxs:distribution%info%rgxe, &
         distribution%info%rgys:distribution%info%rgye, distribution%info%rgzs:distribution%info%rgze):: val
    PetscScalar,intent(out),dimension(0:distribution%b, &
         distribution%info%rgxs:distribution%info%rgxe, &
         distribution%info%rgys:distribution%info%rgye, &
         distribution%info%rgzs:distribution%info%rgze):: out

    PetscInt n
    do n=0,distribution%b
       out(n,distribution%info%xs:distribution%info%xe, &
            distribution%info%ys:distribution%info%ye, &
            distribution%info%zs:distribution%info%ze) = &
              val(distribution%info%xs + distribution%disc%ci(n,X_DIRECTION): &
                  distribution%info%xe + distribution%disc%ci(n,X_DIRECTION), &
                  distribution%info%ys + distribution%disc%ci(n,Y_DIRECTION): &
                  distribution%info%ye + distribution%disc%ci(n,Y_DIRECTION), &
                  distribution%info%zs + distribution%disc%ci(n,Z_DIRECTION): &
                  distribution%info%ze + distribution%disc%ci(n,Z_DIRECTION))
    end do
  end subroutine DistributionGatherValueToDirectionD3

  subroutine DistributionCalcDeltaNorm(distribution, norm)
    type(distribution_type) distribution
    PetscScalar norm
    PetscErrorCode ierr

    PetscScalar denom

    norm =1.d99
    call DistributionRestoreArrays(distribution, ierr)
    if (distribution%track_old_fi) then
      call DMLocalToGlobalBegin(distribution%da_fi,distribution%fi, &
           INSERT_VALUES, distribution%fi_g, ierr)
      call DMLocalToGlobalEnd(distribution%da_fi,distribution%fi, &
           INSERT_VALUES, distribution%fi_g, ierr)
      call VecAXPY(distribution%fi_g_old, -1.d0, distribution%fi_g, ierr)
      call VecPointwiseDivide(distribution%fi_g_old, distribution%fi_g_old, distribution%fi_g, ierr)
      call VecNorm(distribution%fi_g_old, NORM_INFINITY, norm, ierr)
      call VecCopy(distribution%fi_g, distribution%fi_g_old, ierr)
    else if (distribution%track_old_rho) then
      call DMLocalToGlobalBegin(distribution%da_rho,distribution%rho, &
           INSERT_VALUES, distribution%rho_g, ierr)
      call DMLocalToGlobalEnd(distribution%da_rho,distribution%rho, &
           INSERT_VALUES, distribution%rho_g, ierr)
      call VecAXPY(distribution%rho_g_old, -1.d0, distribution%rho_g, ierr)
      call VecPointwiseDivide(distribution%rho_g_old, distribution%rho_g_old, distribution%rho_g, ierr)
      call VecNorm(distribution%rho_g_old, NORM_INFINITY, norm, ierr)
      call VecCopy(distribution%rho_g, distribution%rho_g_old, ierr)
    endif
    call DistributionGetArrays(distribution, ierr)
    CHKERRQ(ierr)
  end subroutine DistributionCalcDeltaNorm


end module LBM_Distribution_Function_module
