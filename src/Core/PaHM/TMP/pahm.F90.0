!----------------------------------------------------------------
!               P R O G R A M   P A H M
!----------------------------------------------------------------
!> @file pahm.F90
!>
!>
!> @brief
!>   Main PaHM program, calls Init, Run and Finalize procedures.
!>
!> @details
!>   1) Initialize PaHM by establishing the logging facilities and calling the subroutine
!>      "GetProgramCmdlArgs" to get possible command line arguments and set the defaults.
!>      During the initialization stage, PaHM reads the mandatory input control file
!>      (defaults to pahm_control.in) to read in the definitions of different variables
!>      used in PaHM.
!>      At this stage we read the mesh/grid of the domain or the generic mesh/grid input file
!>      and the list of best track files supplied by the user.
!>
!>   2) Start the PaHM run (timestepping).
!>
!>   3) Finalize the PaHM run and exit the program.
!>
!> @author Panagiotis Velissariou <panagiotis.velissariou@noaa.gov>
!----------------------------------------------------------------

PROGRAM PaHM

  USE PaHM_DriverMod, ONLY : PaHM_Init, PaHM_Run, PaHM_Finalize

  IMPLICIT NONE

  CALL PaHM_Init()

  CALL PaHM_Run()

  CALL PaHM_Finalize()

END PROGRAM PaHM
