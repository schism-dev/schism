#include "wwm_functions.h"
!> \file wwm_petsc_controller.F90
!> Controls which PETSC Module to use
!> \author    Thomas Huxhorn
!> \date      2012
#ifdef PETSC
    !> controls at compile time which PETSC module will be use
    Module PETSC_CONTROLLER
      implicit none
#include "finclude/petscsysdef.h"

      logical :: init = .false.
      contains

      !> calls the specific PETSC_INIT_* subroutine
      subroutine PETSC_INIT
        use datapool, only: AMETHOD, DBG
        use petscpool, only: petscpoolInit, stageInit, stageSolve, petscErr
#ifdef MPI_PARALL_GRID
        use petsc_parallel, only: PETSC_INIT_PARALLEL
        use petsc_block, only: PETSC_INIT_BLOCK
#endif
        use petsc_seriell, only: PETSC_INIT_SERIELL
        use petscsys
        implicit none

        init = .true.
        ! petscpoolInit must be called BEFORE PETSC_INIT_*
        call petscpoolInit
        call PetscLogStagePush(stageInit, petscErr);CHKERRQ(petscErr)
#ifdef MPI_PARALL_GRID
        if(AMETHOD == 4) then
          call PETSC_INIT_PARALLEL
        else if(AMETHOD == 5) then
          call PETSC_INIT_BLOCK
        else
          write(DBG%FHNDL,*) "PETSC_INIT() you can use AMETHOD=4 "
          write(DBG%FHNDL,*) "for petsc parallel or AMETHOD=5 for petsc block."
          write(DBG%FHNDL,*) "Other AMETHOD numbers are for other solvers."
          write(DBG%FHNDL,*) "AMETHOD =" ,  AMETHOD
        endif
#else
          call PETSC_INIT_SERIELL
#endif
          call PetscLogStagePop(stageInit, petscErr);CHKERRQ(petscErr)
      end subroutine

      !> calls the specific EIMPS_PETSC_{PARALLEL, SERIELL} subroutine
      !> \param ISS frequency. From 1 to MSC
      !> \param IDD direction. From 1 to MDC
      subroutine EIMPS_PETSC(ISS, IDD)
#ifdef MPI_PARALL_GRID
        use petsc_parallel, only: EIMPS_PETSC_PARALLEL
        use petsc_block,    only: EIMPS_PETSC_BLOCK
#endif
        use petsc_seriell,  only: EIMPS_PETSC_SERIELL
        implicit none
        integer, intent(in) :: ISS, IDD

#ifdef MPI_PARALL_GRID
        call EIMPS_PETSC_PARALLEL(ISS, IDD)
#else
        call EIMPS_PETSC_SERIELL(ISS, IDD)
#endif
      end subroutine


      !> clean up memory and calls PetscFinalize()
      subroutine PETSC_FINALIZE()
        use datapool, only: AMETHOD, DBG
        use petscpool,      only: petscErr, petscpoolFinalize, stageFin
#ifdef MPI_PARALL_GRID
        use petsc_parallel, only: PETSC_FINALIZE_PARALLEL
        use petsc_block,    only: PETSC_FINALIZE_BLOCK
#endif
        use petsc_seriell,  only: PETSC_FINALIZE_SERIELL
        use petscsys
        implicit none

        ! if petsc* was not initialized  we must not finalize it
        if(init.eqv..false.) return
          call PetscLogStagePush(stageFin, petscErr);CHKERRQ(petscErr)

          ! petscpoolFinalize must be called BEFORE PETSC_FINALIZE_*
          call petscpoolFinalize()
#ifdef MPI_PARALL_GRID
          if(AMETHOD == 4) then
            call PETSC_FINALIZE_PARALLEL
          else if(AMETHOD == 5) then
            call PETSC_FINALIZE_BLOCK
          else
            write(DBG%FHNDL,*) "PETSC_FINALIZE() you can use AMETHOD=4"
            write(DBG%FHNDL,*) "for petsc parallel or AMETHOD=5 for petsc block."
            write(DBG%FHNDL,*) "Other AMETHOD numbers are for other solvers."
            write(DBG%FHNDL,*) "AMETHOD =" , AMETHOD
          endif
# else
          call PETSC_FINALIZE_SERIELL
# endif
          call PetscFinalize(petscErr)
          
      end subroutine
    end Module
#endif
