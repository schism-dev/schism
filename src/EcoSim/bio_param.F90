!   Copyright 2014 College of William and Mary
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

      MODULE bio_param

!     Routines: 
!        initialize_param      
!        initialize_scalars
!!======================================================================
!! April/May, 2007 - Original code                                     ! 
!!======================================================Marta Rodrigues=
!!                                                                     !
!! This module is adapted from ROMS (mod_param.F + mod_scalars.F):     !
!!                                                                     !
!! February, 2009 - Extended for oxygen cycle                          ! 
!!======================================================================   
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!                                                                      !
!  Tracer parameters:                                                  !
!                                                                      !
!                                                                      !
!  NBIT        Number of biological tracer type variables.              !
!                                                                      !
!  NBands     Number of spectral irradiance bands in EcoSim.           !
!  Nbac       Number of bacteria constituents in EcoSim.               !
!  Ndom       Number of dissolved matter constituents in EcoSim.       !
!  Nfec       Number of fecal matter constituents in EcoSim.           !
!  Nphy       Number of phytoplankton constituents in EcoSim.          !
!  Npig       Number of pigment constituents in EcoSim.                !
!  PHY        EcoSim indices of phytoplankton species considered.      !
!  PIG        EcoSim phytoplankton-pigment matrix.                     !
!                                                                      !
! Marta Rodrigues                                                      !
!  Nzoo       Number of zooplankton constituents in EcoSim.            !
!                                                                      !
!=======================================================================
!

        IMPLICIT NONE
	SAVE

#ifdef USE_SINGLE
        integer, parameter :: r8 = 4  
#else
        integer, parameter :: r8 = 8  
#endif
              
        integer, parameter :: NBands = 60
        integer, parameter :: NAT = 2
! Total number of tracer variables
!        integer :: NT

!
!  Biological tracers parameters.
!
        integer, parameter :: Nbac = 1  
        integer, parameter :: Ndom = 1  !MFR
        integer, parameter :: Nfec = 1  !MFR
        integer, parameter :: Nphy = 1  !MFR
        integer, parameter :: Npig = 7  
! Marta Rodrigues
        integer, parameter :: Nzoo = 1	
	
!
!  Determine number of EcoSim biological tracer. Currently, there is a
!  maximum of seven phytoplankton species and seven different pigments:
!
! [1] small diatom           [1] chlorophyll-a
! [2] large diatom           [2] chlorophyll-b
! [3] small dinoflagellate   [3] chlorophyll-c
! [4] large dinoflagellate   [4] photosythetic carotenoids
! [5] synechococcus          [5] photoprotective carotenoids
! [6] small prochlorococcus  [6] low  urobilin phycoeurythin carotenoids
! [7] large prochlorococcus  [7] high urobilin phycoeurythin carotenoids
!
!  The phytoplankton/pigment matrix is as follows:
!
!               P h y t o p l a n k t o n
!              [1]   [2]   [3]   [4]   [5]   [6]   [7]
!
!       t [7]   0     0     0     0     1     0     0
!       n [6]   0     0     0     0     0     0     0
!       e [5]   1     1     1     1     1     1     1
!       m [4]   1     1     1     1     0     0     0
!       g [3]   1     1     1     1     0     0     0
!       i [2]   0     0     0     0     0     1     1
!       P [1]   1     1     1     1     1     1     1
!
!        integer, parameter, dimension(7,7) :: PIG = reshape (           &
!     &                                      (/ 1, 1, 1, 1, 1, 1, 1,     &
!     &                                         0, 0, 0, 0, 0, 1, 1,     &
!     &                                         1, 1, 1, 1, 0, 0, 0,     &
!     &                                         1, 1, 1, 1, 0, 0, 0,     &
!     &                                         1, 1, 1, 1, 1, 1, 1,     &
!     &                                         0, 0, 0, 0, 0, 0, 0,     &
!     &                                         0, 0, 0, 0, 1, 0, 0 /),  &
!     &                                      (/ 7, 7 /) )

! MFR - Consider only chlorophyll a

        integer, parameter, dimension(7,7) :: PIG = reshape (           &
     &                                      (/ 1, 1, 1, 1, 1, 1, 1,     &
     &                                         0, 0, 0, 0, 0, 0, 0,     &
     &                                         0, 0, 0, 0, 0, 0, 0,     &
     &                                         0, 0, 0, 0, 0, 0, 0,     &
     &                                         0, 0, 0, 0, 0, 0, 0,     &
     &                                         0, 0, 0, 0, 0, 0, 0,     &
     &                                         0, 0, 0, 0, 0, 0, 0 /),  &
     &                                      (/ 7, 7 /) )

!
!  Set phytoplankton species to consider (see above classification):
!
        integer, parameter, dimension(Nphy) :: PHY = (/ 1 /) !MFR

!  MFR - Iron and CDOC flags (0-no calcultions; 1-calculations)
!        Note: If CDOC=0, then RtUVR_flag in ecosim.in must be set to F 
!              If CDOC=1, Ndom must be 2 

     integer, parameter :: IRON=0
     integer, parameter :: CDOC=0  
	
!-----------------------------------------------------------------------
!  Tracer identification indices.
!-----------------------------------------------------------------------
!

! Marta Rodrigues
! Temperature and salinity from SELFE
        integer :: itemp              ! Potential temperature
        integer :: isalt              ! Salinity

        integer, pointer :: idbio(:)  ! Biological tracers
        integer :: iBacC(Nbac)        ! Bacteria, Carbon group
        integer :: iBacN(Nbac)        ! Bacteria, Nitrogen group
        integer :: iBacP(Nbac)        ! Bacteria, Phosphorous group
        integer :: iBacF(Nbac)        ! Bacteria, Iron group
        integer :: iCDMC(2)           ! Color degradational matter
        integer :: iDOMC(2)           ! DOM, Carbon group
        integer :: iDOMN(2)           ! DOM, Nitrogen group
        integer :: iDOMP(2)           ! DOM, Phosphorous group
        integer :: iFecC(2)           ! Fecal matter, Carbon group
        integer :: iFecN(2)           ! Fecal matter, Nitrogen group
        integer :: iFecP(2)           ! Fecal matter, Phosphorous group
        integer :: iFecF(2)           ! Fecal matter, Iron group
        integer :: iFecS(2)           ! Fecal matter, Silica group
        integer :: iPhyC(Nphy)        ! Phytoplankton, Carbon group
        integer :: iPhyN(Nphy)        ! Phytoplankton, Nitrogen group
        integer :: iPhyP(Nphy)        ! Phytoplankton, Phosphorous group
        integer :: iPhyF(Nphy)        ! Phytoplankton, Iron group
        integer :: iPhyS(Nphy)        ! Phytoplankton, Silica group
        integer :: iPigs(Nphy,Npig)   ! Phytoplankton, pigment group
        integer :: iNO3_              ! Nitrate concentration
        integer :: iNH4_              ! Ammonium concentration
        integer :: iPO4_              ! Phosphate concentration
        integer :: iFeO_              ! Iron concentration
        integer :: iSiO_              ! Silica concentration
        integer :: iDIC_              ! Dissolved inorganic Carbon
        integer :: FirstPig           ! Index of first tracer pigment

! Marta Rodrigues
        integer :: iZooC(Nzoo)        ! Zooplankton, Carbon group
	integer :: iZooN(Nzoo)        ! Zooplankton, Nitrogen group
	integer :: iZooP(Nzoo)        ! Zooplankton, Phosphorous group	

        integer :: iDO_               ! Dissolved oxygen concentration
        integer :: iCOD_              ! Chemical oxygen demand concentration
	
!
!  EcoSim group names used on standard output.
!
        character (len=16), dimension(Nbac) :: BacName
        character (len=11), dimension(Ndom) :: DomName
        character (len=13), dimension(Nfec) :: FecName
        character (len=21), dimension(Nphy) :: PhyName
	
! Marta Rodrigues
        character (len=21), dimension(Nzoo) :: ZooName	

        character (len=39), dimension(7) :: PigName =                   &
     &            (/ 'chlorophyll-a                          ',         &
     &               'chlorophyll-b                          ',         &
     &               'chlorophyll-c                          ',         &
     &               'photosythetic carotenoids              ',         &
     &               'photoprotective carotenoids            ',         &
     &               'low urobilin phycoeurythin carotenoids ',         &
     &               'high urobilin phycoeurythin carotenoids' /)

!
! MFR


! Marta Rodrigues
! Other variables (for date and time)

       real(r8) :: day, month, year, hour, minutes, seconds, yday
 
!-----------------------------------------------------------------------
! Parameters
!-----------------------------------------------------------------------
!       real(r8), parameter :: pi = 3.14159265
       real(r8), parameter :: deg2rad = 3.14159265/180._r8
       real(r8), parameter :: rad2deg = 180._r8/3.14159265 
   

!-----------------------------------------------------------------------
!  Physical constants.   
!-----------------------------------------------------------------------
!
!    Cp            Specific heat for seawater (Joules/Kg/degC).
!    Csolar        Solar irradiantion constant (W/m2).
!    Eradius       Earth equatorial radius (m).
!    StefBo        Stefan-Boltzmann constant (W/m2/K4).
!    emmiss        Infrared emmissivity.
!    g             Acceleration due to gravity (m/s2).
!    gorho0        gravity divided by mean density anomaly.
!    rhow          fresh water density (kg/m3).
!    vonKar        von Karman constant.
!
        real(r8) :: Cp = 3985.0_r8              ! Joules/kg/degC
        real(r8) :: Csolar = 700.0_r8           ! 1360-1380 W/m2
        real(r8) :: Eradius = 6371315.0_r8      ! m
        real(r8) :: StefBo = 5.67E-8_r8         ! Watts/m2/K4
        real(r8) :: emmiss = 0.97_r8            ! non_dimensional
        real(r8) :: rhow = 1000.0_r8            ! kg/m3         
        real(r8) :: g = 9.81_r8                 ! m/s2
        real(r8) :: gorho0                      ! m4/s2/kg
        real(r8) :: vonKar = 0.41_r8            ! non-dimensional
!

!  MFR - verificar se alguns destes são usados!!!
!-----------------------------------------------------------------------
!  Water clarity parameters.
!-----------------------------------------------------------------------
!
!    lmd_mu1       Reciprocal of the absorption coefficient for solar
!                    wavelength band 1 as a function of the Jerlov
!                    water type.
!    lmd_mu2       Reciprocal of the absorption coefficient for solar
!                    wavelength band 2 as a function of the Jerlov
!                    water type.
!    lmd_r1        Fraction of total radiance for wavelength band 1 as
!                    a function of the Jerlov water type.
!
        real(r8), dimension(5) :: lmd_mu1 =                             &
     &            (/ 0.35_r8, 0.6_r8, 1.0_r8, 1.5_r8, 1.4_r8 /)

        real(r8), dimension(5) :: lmd_mu2 =                             &
     &            (/ 23.0_r8, 20.0_r8, 17.0_r8, 14.0_r8, 7.9_r8 /)

        real(r8), dimension(5) :: lmd_r1 =                              &
     &            (/ 0.58_r8, 0.62_r8, 0.67_r8, 0.77_r8, 0.78_r8 /)	
	
!
!  The following parameters are derived or only known during execution.
!
        common /param/ NBIT, NBIT2
	
        !integer,save :: NBIT,NBIT2
        integer :: NBIT,NBIT2

	CONTAINS
!=======================================================================
        SUBROUTINE initialize_param
!                                                                      !
!  This routine initializes several parameters in module "mod_param"   !
!  for all nested grids.                                               ! 
!                                                                      !
!=======================================================================
!Marta Rodrigues                                                       !
!  Kept only the ones with respect to EcoSim                           !
!                                                                      !
!=======================================================================
!
!  Local variable declarations
!

        integer :: i, j
!
!-----------------------------------------------------------------------
!  Determine number of EcoSim total bio-optical constituents:
!-----------------------------------------------------------------------
!
! MFR - Changed to allow that no calculations of CDOC (all DOC would
!       be in DOC pool) and Fe
!
!       Nutrients: NO3, NH4, PO4, SiO, DIC  (5 + FeO)
!        Bacteria: C, N, P                  (Nbac*3 + Fe)
!             DOM: C, N, P                  (Ndom*3 + CDM)
!           Fecal: C, N, P, Si              (Nfec*4 + Fe)
!    Phytoplakton: C, N, P                  (Nphy*3 + Si + Fe)
!        Pigments: look table
!
! Marta Rodrigues
!        Zooplankton: C, N, P               (Nzoo*3)
!        Dissolved oxygen: DO               (1)  
!        Chemical oxygen demand: COD        (1) 


        NBIT=5+(Nbac*3)+(Ndom*3)+(Nfec*4)+(Nphy*3)+(Nzoo*3)+2
!
!  Add phytoplankton silica constituents.
!
        DO i=1,Nphy
          IF (PHY(i).le.2) NBIT=NBIT+1
        END DO
!
!  Add pigments.  Check phytoplankton-pigment table for values greater
!  than zero.
!
        DO j=1,Npig
          DO i=1,Nphy
            IF (PIG(PHY(i),j).eq.1) NBIT=NBIT+1
          END DO
        END DO
!  MFR - Add iron (if iron calculations are made...)

        IF (IRON==1) NBIT=NBIT+1+Nbac*1+Nfec*1+Nphy*1

! MFR - Add CDOC (if colored and non-colored DOC are differentiate)

        IF (CDOC==1) NBIT=NBIT+Ndom*1

! Marta Rodrigues
! Total number of tracer variables + S + T
!      NT = NBIT + NAT   
       NBIT2=NBIT+2 !including T,S

      RETURN	
      END SUBROUTINE
!=======================================================================
             
!=======================================================================
     SUBROUTINE initialize_scalars
!                                                                      !
!  This routine initializes several variables in module "mod_scalars"  !
!  for all nested grids.                                               !
!                                                                      !
!=======================================================================
!Marta Rodrigues                                                       !
!  Kept only the ones with respect to EcoSim                           !
!                                                                      !
!=======================================================================

!
!  Local variable declarations.
!
      integer :: i, ic, j

      character (len=21), dimension(7) :: PhyGroups =                   &
     &                                 (/ 'small diatom         ',      &
     &                                    'large diatom         ',      &
     &                                    'small dinoflagellate ',      &
     &                                    'large dinoflagellate ',      &
     &                                    'synechococcus        ',      &
     &                                    'small prochlorococcus',      &
     &                                    'large prochlorococcus' /)

      allocate ( idbio(NBIT) )
!
!---------------------------------------------------------------------
!  Set tracer identification indices.
!---------------------------------------------------------------------
!
!      itemp=1
!      isalt=2
!      ic=NAT !MFR - NAT - number of active tracer variables (usually = 2
             !            for temperature and salinity)

! Marta Rodrigues
! In SELFE T and S are in diferent arrays from the tracers array...
       
       ic=0

!
!  Set bio-optical tracer indices.
!
      DO i=1,NBIT
        idbio(i)=ic+i
      END DO

      iDIC_=ic+1

      IF(IRON==1)THEN
       iFeO_=ic+2
       iNH4_=ic+3
       iNO3_=ic+4
       iPO4_=ic+5
       iSiO_=ic+6
       ic=ic+6
      ELSE
       iNH4_=ic+2
       iNO3_=ic+3
       iPO4_=ic+4
       iSiO_=ic+5
       ic=ic+5
      END IF

      DO i=1,Nbac
       IF(IRON==1)THEN 
        iBacC(i)=ic+1
        iBacF(i)=ic+2
        iBacN(i)=ic+3
        iBacP(i)=ic+4
        ic=ic+4
       ELSE
        iBacC(i)=ic+1
        iBacN(i)=ic+2
        iBacP(i)=ic+3
        ic=ic+3
       END IF
      END DO

      IF(CDOC==1.AND.Ndom==2)THEN 
        iCDMC(1)=ic+1
        iDOMC(1)=ic+2
        iDOMN(1)=ic+3
        iDOMP(1)=ic+4
        iCDMC(2)=ic+5
        iDOMC(2)=ic+6
        iDOMN(2)=ic+7
        ic=ic+7
!      ELSEIF(CDOC==1.AND.Ndom==1)THEN
!        iCDMC(1)=ic+1
!        iDOMC(1)=ic+2
!        iDOMN(1)=ic+3
!        iDOMP(1)=ic+4
!        ic=ic+4
      ELSEIF (CDOC==0.AND.Ndom==2)THEN
        iDOMC(1)=ic+1
        iDOMN(1)=ic+2
        iDOMP(1)=ic+3
        iDOMC(2)=ic+4
        iDOMN(2)=ic+5
        ic=ic+5
      ELSE
        iDOMC(1)=ic+1
        iDOMN(1)=ic+2
        iDOMP(1)=ic+3
        ic=ic+3
      ENDIF

      DO i=1,Nfec
       IF(IRON==1)THEN
        iFecC(i)=ic+1
        iFecF(i)=ic+2
        iFecN(i)=ic+3
        iFecP(i)=ic+4
        iFecS(i)=ic+5
        ic=ic+5
       ELSE
        iFecC(i)=ic+1
        iFecN(i)=ic+2
        iFecP(i)=ic+3
        iFecS(i)=ic+4
        ic=ic+4
       END IF
      END DO

      DO i=1,Nphy
       IF(IRON==1)THEN
        iPhyC(i)=ic+1
        iPhyF(i)=ic+2
        iPhyN(i)=ic+3
        iPhyP(i)=ic+4
        ic=ic+4
       ELSE
        iPhyC(i)=ic+1
        iPhyN(i)=ic+2
        iPhyP(i)=ic+3
        ic=ic+3
       ENDIF
      END DO

      DO i=1,Nphy
        IF (PHY(i).le.2) THEN
          ic=ic+1
          iPhyS(i)=ic
        ELSE
          iPhyS(i)=0
        END IF
      END DO
      FirstPig=ic+1
      DO j=1,Npig
        DO i=1,Nphy
          iPigs(i,j)=0
          IF (PIG(PHY(i),j).eq.1) THEN
            ic=ic+1
            iPigs(i,j)=ic
          END IF
        END DO
      END DO
!
! Marta Rodrigues
      DO i=1,Nzoo
        iZooC(i)=ic+1
	iZooN(i)=ic+2
	iZooP(i)=ic+3
	ic=ic+3
      END DO

      iDO_=ic+1
      iCOD_=ic+2

! Marta Rodrigues
! For salinity and temperature dimensions in Bio(mnv,mve,NT) array
      itemp = ic+3
      isalt = ic+4
	  

!  Set EcoSim group variable names.
!
      DO i=1,Nphy
        PhyName(i)=PhyGroups(PHY(i))
      END DO
      DO i=1,Nbac
        WRITE (BacName(i),'(a,1x,i1)') 'Bacteria Group', i
      END DO
      DO i=1,Ndom
        WRITE (DomName(i),'(a,1x,i1)') 'DOM Group', i
      END DO
      DO i=1,Nfec
        WRITE (FecName(i),'(a,1x,i1)') 'Fecal Group', i
      END DO

! Marta Rodrigues
      DO i=1,Nzoo
        WRITE (ZooName(i),'(a,1x,i1)') 'Zooplankton Group', i
      END DO	    
      
      RETURN
      END SUBROUTINE initialize_scalars    
!=====================================================================      

      END MODULE bio_param
