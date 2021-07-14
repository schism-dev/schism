!=======================================================================
!
! This module defines a variety of physical and numerical constants
! used throughout the ice model 
!
! author Elizabeth C. Hunke, LANL

      module icedrv_constants

      use icedrv_kinds

      implicit none

      !-----------------------------------------------------------------
      ! file units
      !-----------------------------------------------------------------

      integer (kind=int_kind), parameter, public :: &
         ice_stdin  =  5,   & ! reserved unit for standard input
         ice_stdout =  88,   & ! reserved unit for standard output
         ice_stderr =  87,   & ! reserved unit for standard error
         nu_nml     = 10,   &  ! unit for namelist
         nu_restart = 12,   &  ! unit for restart file
         nu_dump    = 13,   &  ! unit for dump file
         nu_forcing = 14,   &  ! unit for forcing file
         nu_open_clos = 15, &  ! unit for SHEBA forcing file
         nu_diag      = 86, &  ! unit for diagnostic output
         nu_diag_out  = 103

      !-----------------------------------------------------------------
      ! numerical constants
      !-----------------------------------------------------------------

      real (kind=dbl_kind), parameter, public :: &
         c0   = 0.0_dbl_kind, &
         c1   = 1.0_dbl_kind, &
         c1p5 = 1.5_dbl_kind, &
         c2   = 2.0_dbl_kind, &
         c3   = 3.0_dbl_kind, &
         c4   = 4.0_dbl_kind, &
         c5   = 5.0_dbl_kind, &
         c6   = 6.0_dbl_kind, &
         c8   = 8.0_dbl_kind, &
         c10  = 10.0_dbl_kind, &
         c15  = 15.0_dbl_kind, &
         c16  = 16.0_dbl_kind, &
         c20  = 20.0_dbl_kind, &
         c24  = 24.0_dbl_kind, &
         c25  = 25.0_dbl_kind, &
         c100 = 100.0_dbl_kind, &
         c1000= 1000.0_dbl_kind, &
         p001 = 0.001_dbl_kind, &
         p01  = 0.01_dbl_kind, &
         p1   = 0.1_dbl_kind, &
         p2   = 0.2_dbl_kind, &
         p4   = 0.4_dbl_kind, &
         p5   = 0.5_dbl_kind, &
         p6   = 0.6_dbl_kind, &
         p05  = 0.05_dbl_kind, &
         p15  = 0.15_dbl_kind, &
         p25  = 0.25_dbl_kind, &
         p75  = 0.75_dbl_kind, &
         p333 = c1/c3, &
         p666 = c2/c3
         
      !-----------------------------------------------------------------
      ! physical constants
      !-----------------------------------------------------------------

      real (kind=dbl_kind), parameter, public :: &
         omega     = 7.292e-5_dbl_kind   ! angular velocity of earth(rad/sec)

      !-----------------------------------------------------------------
      ! numbers used outside the column package
      !-----------------------------------------------------------------

      real (kind=dbl_kind), parameter, public :: &
         c9   = 9.0_dbl_kind, &
         c12  = 12.0_dbl_kind, &
         c30  = 30.0_dbl_kind, &
         c180 = 180.0_dbl_kind, &
         c360 = 360.0_dbl_kind, &
         c365 = 365.0_dbl_kind, &
         c400 = 400.0_dbl_kind, &
         c3600= 3600.0_dbl_kind, &
         p025 = 0.025_dbl_kind, &
         p166 = c1/c6, &
         p111 = c1/c9, &
         p055 = p111*p5, &
         p027 = p055*p5, &
         p222 = c2/c9, &
         eps13  = 1.0e-13_dbl_kind, &
         eps16  = 1.0e-16_dbl_kind

      !-----------------------------------------------------------------
      ! conversion factors
      !-----------------------------------------------------------------

      real (kind=dbl_kind), parameter, public :: &
         mps_to_cmpdy  = 8.64e6_dbl_kind   ! m per s to cm per day 

!=======================================================================

      end module icedrv_constants

!=======================================================================
