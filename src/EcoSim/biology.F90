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


  MODULE biology

!  Routines
!    initialize_biology

!!======================================================================
!! April/May, 2007 - Original code                                     ! 
!!======================================================Marta Rodrigues=
!!                                                                     !
!! This module is adapted from ROMS:                                   !
!!                                                                     ! 
!! February, 2009 - Extended for oxygen cycle                          !
!!======================================================================

!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group       W. Paul Bissett   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  EcoSim Model Phytoplaknton Parameters:                              !
!                                                                      !
!                                                                      !
!  HsNO3          Half-saturation for phytoplankton NO3 uptake         !
!                   (micromole_NO3/liter).                             !
!  HsNH4          Half-saturation for phytoplankton NH4 uptake         !
!                   (micromole_NH4/liter).                             !
!  HsSiO          Half-saturation for phytoplankton SiO uptake         !
!                   (micromole_SiO/liter).                             !
!  HsPO4          Half-saturation for phytoplankton PO4 uptake         !
!                   (micromole_PO4/liter).                             !
!  HsFe           Half-saturation for phytoplankton Fe uptake          !
!                  (micromole_Fe/liter).                               !
!  GtALG_max      Maximum 24 hour growth rate (1/day).                 !
!  PhyTbase       Phytoplankton temperature base for exponential       !
!                   response to temperature (Celsius).                 !
!  PhyTfac        Phytoplankton exponential temperature factor         !
!                   (1/Celsius).                                       !
!  BET_           Nitrate uptake inhibition for NH4 (l/micromole).     !
!  maxC2nALG      Maximum phytoplankton C:N ratio                      !
!                   (micromole_C/micromole_N).                         !
!  minC2nALG      Balanced phytoplankton C:N ratio                     !
!                   (micromole_C/micromole_N).                         !
!  C2nALGminABS   Absolute minimum phytoplankton C:N ratio             !
!                   (micromole_C/micromole_N).                         !
!  maxC2SiALG     Maximum phytoplankton C:Si ratio                     !
!                   (micromole_C/micromole_Si).                        !
!  minC2SiALG     Balanced phytoplankton C:Si ratio                    !
!                   (micromole_C/micromole_Si).                        !
!  C2SiALGminABS  Absolute minimum phytoplankton C:Si ratio            !
!                  (micromole_C/micromole_Si).                         !
!  maxC2pALG      Maximum phytoplankton C:P ratio                      !
!                   (micromole_C/micromole_P).                         !
!  minC2pALG      Balanced phytoplankton C:P ratio                     !
!                   (micromole_C/micromole_P).                         !
!  C2pALGminABS   Absolute minimum phytoplankton C:P ratio             !
!                   (micromole_C/micromole_P).                         !
!  maxC2FeALG     Maximum phytoplankton C:Fe ratio                     !
!                   (micromole_C/micromole_Fe).                        !
!  minC2FeALG     Balanced phytoplankton C:Fe ratio                    !
!                   (micromole_C/micromole_Fe).                        !
!  C2FeALGminABS  Absolute minimum phytoplankton C:Fe ratio            !
!                   (micromole_C/micromole_Fe).                        !
!  qu_yld         Maximum quantum yield                                !
!                   (micromole_C/micromole_quanta).                    !
!  E0_comp        Compensation light level (micromole_quanta).         !
!  E0_inhib       Light level for onset of photoinhibition             !
!                   (micromole_quanta).                                !
!  hib_fac      Exponential decay factor for light limited growth    !
!                   (1/micromole_quanta).                              !
!  C2Chl_max      Maximum lighted limited (nutrient replete) C:Chl     !
!                   ratio (microgram_C/microgram_Chl).                 !
!  mxC2Cl         Rate of change in the lighted limited C:Chl ratio    !
!                   (microgram_C/microgram_Chl/micromole_quanta).      !
!  b_C2Cl         Minimum lighted limited (nutrient replete) C:Chl     !
!                   ratio (microgram_C/microgram_Chl).                 !
!  mxC2Cn         Rate of change in the nutrient limited C:Chl ratio   !
!                   [(ug_C/ug_Chl)/(umole_C/umole_N)].                 !
!  b_C2Cn         Minimum nutrient limited C:Chl ratio                 !
!                   (microgram_C/microgram_Chl).                       !
!  mxPacEff       Rate of change in package effect                     !
!                   [1/(microgram_C/microgram_Chl)].                   !
!  b_PacEff       Maximum package effect                               !
!                   [1/(microgram_C/microgram_Chl)].                   !
!  mxChlB         Rate of change in the Chl_b:Chl_a ratio              !
!                   [(ug_Chl_b/ug_Chl_a)/(ug_C/ug_Chl_a)].             !
!  b_ChlB         Maximum Chl_b:Chl_a ratio                            !
!                   (microgram_Chl_b/microgram_Chl_a).                 !
!  mxChlC         Rate of change in the Chl_c:Chl_a ratio              !             
!                   [(ug_Chl_c/ug_Chl_a)/(ug_C/ug_Chl_a)].             !
!  b_ChlC         Maximum Chl_c:Chl_a ratio                            !
!                   (microgram_Chl_c/microgram_Chl_a).                 !
!  mxPSC          Rate of change in the PSC:Chl_a ratio                !
!                   [(ug_PSC/ug_Chl_a)/(ug_C/ug_Chl_a)].               !
!  b_PSC          Maximum PSC:Chl_a ratio                              !
!                  (microgram_Chl_c/microgram_Chl_a).                  !
!  mxPPC          Rate of change in the PPC:Chl_a ratio                !
!                   [(ug_PPC/ug_Chl_a)/(ug_C/ug_Chl_a)].               !
!  b_PPC          Maximum PPC:Chl_a ratio                              !
!                  (microgram_Chl_c/microgram_Chl_a).                  !
!  mxLPUb         Rate of change in the LPUb:Chl_a ratio               !
!                   [(ug_LPUb/ug_Chl_a)/(ug_C/ug_Chl_a)].              !
!  b_LPUb         Maximum LPUb:Chl_a ratio                             !
!                   (microgram_HPUb/microgram_Chl_a).                  !
!  mxHPUb         Rate of change in the HPUb:Chl_a ratio               !
!                   [(ug_HPUb/ug_Chl_a)/(ug_C/ug_Chl_a)].              !
!  b_HPUb         Maximum HPUb:Chl_a ratio                             !
!                   (microgram_HPUb/microgram_Chl_a).                  !
!  FecDOC         Proportion of grazing stress which is apportioned    !
!                   to DOM (nondimensional).                           !
!  FecPEL         Proportion of grazing stress which is apportioned    !
!                   to fecal pellets (nondimesional).                  !
!  FecCYC         Proportion of grazing stress which is apportioned    !
!                   to direct remineralization (nondimensional).       !
!  ExALG          Proportion of daily production that is lost to       !
!                   excretion (nondimensional).                        !
!  WS             Phytoplankton sinking speed (meters/day).            !
!----------------------------------------------------------------------!
! Marta Rodrigues                                                      !
!  HsGRZ          Phytoplankton celular lises parameter                !
!                   (nondimensional).                                  !
!  basalPhy       Phytoplankton basal specific respiration rate (/d)   !
!  QPhy           Temperature coefficient factor for phytoplankton     !
!                   respiration (nondimensional)                       !
!  gamaPhy        Fraction of phytoplankton assimilated production     !
!                   (nondimensional)                                   !
!----------------------------------------------------------------------!
!  MinRefuge      Refuge Phytoplankton population (micromole_C/liter). !
!  RefugeDep      Maximum Refuge Phytoplankton depth (meters).         !
!  Norm_Vol       Normalized Volume factor (nondimensional).           !
!  Norm_Surf      Normalized surface area factor (nondimensional).     !
!  HsDOP          Half Saturation Constant for DOP uptake              !
!                   (micromole_DOP/liter).                             !
!  C2pALKPHOS     C:P ratio where DOP uptake begins                    !
!                   (micromole_C/micromole_P).                         !
!  HsDON          Half Saturation Constant for DON uptake              !
!                   (micromole_DON/liter).                             !
!  C2nNupDON      C:N ratio where DON uptake begins                    !
!                   (micromole_C/micromole_N).                         !
!                                                                      !
! Bacteria Parameters:                                                 !
!                                                                      !
!  HsDOC_ba       Half saturation constant for bacteria DOC uptake     !
!                   (micromole_DOC/liter).                             !
!  GtBAC_max      Maximum 24 hour bacterial growth rate (1/day).       !
!  BacTbase       Phytoplankton temperature base for exponential       !
!                   response to temperature, (Celsius).                !
!  BacTfac        Phytoplankton exponential temperature factor         !
!                   (1/Celsius).                                       !
!  C2nBAC         Carbon to Nitrogen ratio of Bacteria                 !
!                   (micromole_C/micromole_N).                         !
!  C2pBAC         Carbon to Phosphorus ratio of Bacteria               !
!                   (micromole_C/micromole_P).                         !
!  C2FeBAC        Carbon to Iron ratio of Bacteria                     !
!                   (micromole_C/micromole_Fe)                         !
!  BacDOC         Proportion of bacterial grazing stress which is      !
!                   apportioned to DOM (nondimensional).               !
!  BacPEL         Proportion of bacterial grazing stress which is      !
!                   apportioned to fecal pellets (nondimensional).     !
!  BacCYC         Proportion of bacterial grazing stress which is      !
!                   apportioned to direct remineralization             !
!                   (nondimensional).                                  !
!  ExBAC_c        Bacterial recalcitrant carbon excretion as a         !
!                   proportion of uptake (nondimensional)              !
!  ExBacC2N       Bacterial recalcitrant excretion carbon to nitrogen  !
!                   ratio (micromole_C/micromole_N).                   !
!  Bac_Ceff       Bacterial gross growth carbon efficiency             !
!                   (nondimensional).                                  !
! Marta Rodrigues comment (an alternative formulation for nitrification!
! was added for nitrification - Vichi et al., 2007)                    !
!  RtNIT          Maximum bacterial nitrification rate (1/day).        !
!  HsNIT          Half saturation constant for bacterial nitrification !
!                   (micromole NH4/liter or micromole O2/liter, if     !
!                    Vichi et al., 2007 is used)                       !
!  QN             Nitrification temperature coefficient                !
!                   (nondimensionsal)                                  !
!----------------------------------------------------------------------!
! Marta Rodrigues:                                                     !
! For oxygen respiration and an alternative formulation for            !
!   nitrification                                                      !
!                                                                      !
!  basalBac       Bacterioplankton basal specific respiration rate (/d)!
!  QBac           Temperature coefficient factor for bacterioplankton  !
!                   respiration (nondimensional)                       !
!  GEE0C          Decrease in growth efficiency under anoxic           !
!                   conditions (nondimensional)                        !
!  HsBacO         Half-saturation constant for oxygen limitation       !
!                   (micromole O2/liter)                               !
!----------------------------------------------------------------------!
!                                                                      ! 
! Dissolved Organic Matter Parameters:                                 !
!                                                                      !
!  cDOCfrac_c     Colored fraction of DOC production from              !
!                   phytoplankton and bacterial losses                 !
!                   (nondimensional).                                  !
!  RtUVR_DIC      UV degradation of DOC into DIC at 410 nanometers     !
!                   (micromole/meter/liter/hour).                      !
!  RtUVR_DIC      UV degradation of DOC into colorless labile DOC at   !
!                   410 nanometers (micromole/meter/liter/hour).       !
!                                                                      !
! Fecal and detritus Parameters:                                       !
!                                                                      !
!  WF             Fecal sinking flux (meters/day).                     !
!  RegTbase       Fecal regeneration temperature base for exponential  !
!                   response to temperature (Celsius).                 !
!  RegTfac        Fecal regeneration exponential temperature factor    !
!                   (1/Celsius).                                       !
!  RegCR          Fecal carbon regeneration rate (1/day).              !
!  RegNR          Fecal nitrogen regeneration rate (1/day).            !
!  RegSR          Fecal silica regeneration rate (1/day).              !
!  RegPR          Fecal phosphorus regeneration rate (1/day).          !
!  RegFR          Fecal iron regeneration rate (1/day).                !
!                                                                      !
!======================================================================!
!
!====================================================Marta Rodrigues===!
! Zooplankton parameters:                                              !
!                                                                      !
!  ZooDOC         Proportion of zooplankton which is apportioned       !
!                   to DOM (nondimensional).                           !
!  ZooPEL         Proportion of zooplankton which is apportioned       !
!                   to fecal pellets (nondimensional).                 ! 
!  ZooCYC         Proportion of zooplankton which is apportioned       !
!                   to direct remineralization (nondimensional).       !
!  DeltaZoo       Availability of prey to predator (nondimensional)    !
!  EfcCap         Capture efficiency of zooplankton (nondimensional)   !
!  HsZoo          Half-saturation constant for total food ingestion    !
!                   (micromole_C/liter)                                ! 
!  EfcPrd         Assimilation efficiency of zooplankton's predators   !
!                   (nondimensional)                                   !
!  ExZOO          Zooplankton excretion rate (1/day)                   !
!  GZ             Zooplankton losses due to natural mortality          !
!                   and predation (1/day)                              !
!  basalZoo       Zooplankton basal specific respiration rate (/d)     !
!  QZoo           Temperature coefficient factor for zooplankton       !
!                   respiration (nondimensional)                       !
!  etaZoo         Zooplankton assimilation efficiency (nondimensional) !       
!                                                                      !
! Oxygen parameters:                                                   !
!                                                                      !
!  omegaO2C       Stoichiometric coefficient to convert carbon units   !
!                   to oxygen units (nondimensional)                   !
!  omegaO2N       Stoichiometric coefficient to convert nitrogen units !
!                   to oxygen units (nondimensional)                   !
!  omegaS2O       Stoichiometric coefficient to convert sulphur units  !
!                   to oxygen units (nondimensional)                   !
!  omegaO2NDenit  Stoichiometric coefficient to denitrification        !
!                   reaction (mmol O2/mmol N)                          !
!  QWind          Dimensionless coefficient for reaeration             !
!  psiWind        Dimensionless coefficient for OW (Hull et al., 2008) !
!  alfaWind       Dimensionless coefficient for OW (Hull et al., 2008) !
!  RtDENIT        Specific denitrification rate (/d)                   !
!  MDenit         Reference anoxic mineralization rate (micromol O2/   !
!                   /liter.d)                                          !
!  reox_COD       Specific reoxidation rate (/d)                       !
!  QCOD           Temperature coefficient for COD (nondimensional)     !
!  HsCOD          Half-saturation constant for COD (micromole O2/liter)!
! DIC reaeration                                                       !
!  pCO2a          CO2 partial pressure in the air (atm)                !
!  pH             pH in the water                                      !
!=======================================================================
!
! Marta Rodrigues - Nov/2015 ===========================================
! Atmospheric variables for specif.F90                                 ! 
!                                                                      ! 
! cloud  Cloud cover @ elements                                        !
! flag_cloud Flag to choose the input for cloud cover (currently       ! 
!            only 1 is available                                       !
! hair    Relative humidity @ elements                                 !
! tair    Air temperature @ elements                                   !
! pair    Atm. pressures @ elements                                    !
! uwind   U-wind @ elements                                            !
! vwind   V-wind @ elements                                            !
!=======================================================================
 

      USE bio_param
      USE eclight
      USE schism_glbl, only:wsett,irange_tr

      IMPLICIT NONE
      SAVE
!

!    Ngrids
! Marta Rodrigues
! Only one grid
! All real variables are initialized as -99.d0

!-----------------------------------------------------------------------
!  Standard input parameters.
!-----------------------------------------------------------------------
!
!  Number of biological iterations.
!
!MFR,Apr/2017
!      integer :: BioIter 

! Marta Rodrigues
! Zooplankton flag

      integer, dimension(Nzoo) :: zoo_sp
!
!  Control flags.
!
      logical :: RtUVR_flag
      logical :: NFIX_flag
      logical :: Regen_flag

      integer :: NIT_flag, DENIT_flag, REAER_flag
      integer :: flag_cloud = 1
!
!  Phytoplankton parameters.
!
      real(r8), dimension(Nphy) :: HsNO3 = -99.d0
      real(r8), dimension(Nphy) :: HsNH4 = -99.d0 
      real(r8), dimension(Nphy) :: HsSiO = -99.d0
      real(r8), dimension(Nphy) :: HsPO4 = -99.d0
      real(r8), dimension(Nphy) :: HsFe = -99.d0
      real(r8), dimension(Nphy) :: GtALG_max = -99.d0
      real(r8), dimension(Nphy) :: PhyTbase = -99.d0
      real(r8), dimension(Nphy) :: PhyTfac = -99.d0
      real(r8), dimension(Nphy) :: BET_ = -99.d0
      real(r8), dimension(Nphy) :: maxC2nALG = -99.d0
      real(r8), dimension(Nphy) :: minC2nALG = -99.d0
      real(r8), dimension(Nphy) :: C2nALGminABS = -99.d0
      real(r8), dimension(Nphy) :: maxC2SiALG = -99.d0
      real(r8), dimension(Nphy) :: minC2SiALG = -99.d0
      real(r8), dimension(Nphy) :: C2SiALGminABS = -99.d0
      real(r8), dimension(Nphy) :: maxC2pALG = -99.d0
      real(r8), dimension(Nphy) :: minC2pALG = -99.d0
      real(r8), dimension(Nphy) :: C2pALGminABS = -99.d0
      real(r8), dimension(Nphy) :: maxC2FeALG = -99.d0
      real(r8), dimension(Nphy) :: minC2FeALG = -99.d0
      real(r8), dimension(Nphy) :: C2FeALGminABS = -99.d0
      real(r8), dimension(Nphy) :: qu_yld = -99.d0
      real(r8), dimension(Nphy) :: E0_comp = -99.d0
      real(r8), dimension(Nphy) :: E0_inhib = -99.d0
      real(r8), dimension(Nphy) :: inhib_fac = -99.d0
      real(r8), dimension(Nphy) :: C2CHL_max = -99.d0
      real(r8), dimension(Nphy) :: mxC2Cl = -99.d0
      real(r8), dimension(Nphy) :: b_C2Cl = -99.d0
      real(r8), dimension(Nphy) :: mxC2Cn = -99.d0
      real(r8), dimension(Nphy) :: b_C2Cn = -99.d0
      real(r8), dimension(Nphy) :: mxPacEff = -99.d0
      real(r8), dimension(Nphy) :: b_PacEff = -99.d0
      real(r8), dimension(Nphy) :: mxChlB = -99.d0
      real(r8), dimension(Nphy) :: b_ChlB = -99.d0
      real(r8), dimension(Nphy) :: mxChlC = -99.d0
      real(r8), dimension(Nphy) :: b_ChlC = -99.d0
      real(r8), dimension(Nphy) :: mxPSC = -99.d0
      real(r8), dimension(Nphy) :: b_PSC = -99.d0
      real(r8), dimension(Nphy) :: mxPPC = -99.d0
      real(r8), dimension(Nphy) :: b_PPC = -99.d0
      real(r8), dimension(Nphy) :: mxLPUb = -99.d0
      real(r8), dimension(Nphy) :: b_LPUb = -99.d0
      real(r8), dimension(Nphy) :: mxHPUb = -99.d0
      real(r8), dimension(Nphy) :: b_HPUb = -99.d0
      real(r8), dimension(Nphy) :: FecDOC = -99.d0
      real(r8), dimension(Nphy,2) :: FecPEL = -99.d0
      real(r8), dimension(Nphy) :: FecCYC = -99.d0
      real(r8), dimension(Nphy) :: ExALG = -99.d0
      real(r8), dimension(Nphy) :: WS = -99.d0
      real(r8), dimension(Nphy) :: HsGRZ = -99.d0
      real(r8), dimension(Nphy) :: basalPhy = -99.d0
      real(r8), dimension(Nphy) :: QPhy = -99.d0
      real(r8), dimension(NPhy) :: gamaPhy = -99.d0
      real(r8), dimension(Nphy) :: MinRefuge = -99.d0
      real(r8), dimension(Nphy) :: RefugeDep = -99.d0
      real(r8), dimension(Nphy) :: Norm_Vol = -99.d0
      real(r8), dimension(Nphy) :: Norm_Surf = -99.d0
      real(r8), dimension(Nphy) :: HsDOP = -99.d0
      real(r8), dimension(Nphy) :: C2pALKPHOS = -99.d0
      real(r8), dimension(Nphy) :: HsDON = -99.d0
      real(r8), dimension(Nphy) :: C2nNupDON = -99.d0
!
!  Bacteria parameters.
!
      real(r8), dimension(Nbac) :: HsDOC_ba = -99.d0
      real(r8), dimension(Nbac) :: GtBAC_max = -99.d0
      real(r8), dimension(Nbac) :: BacTbase = -99.d0
      real(r8), dimension(Nbac) :: BacTfac = -99.d0
      real(r8), dimension(Nbac) :: basalBac = -99.d0
      real(r8), dimension(Nbac) :: QBac = -99.d0
      real(r8), dimension(Nbac) :: GEE0C = -99.d0
      real(r8), dimension(Nbac) :: HsBacO = -99.d0
      real(r8) :: C2nBAC = -99.d0
      real(r8) :: C2pBAC = -99.d0
      real(r8) :: C2FeBAC = -99.d0
      real(r8) :: BacDOC = -99.d0
      real(r8) :: BacPEL = -99.d0
      real(r8) :: BacCYC = -99.d0
      real(r8) :: ExBAC_c = -99.d0
      real(r8) :: ExBacC2N = -99.d0
      real(r8) :: Bac_Ceff = -99.d0
      real(r8) :: RtNIT = -99.d0
      real(r8) :: HsNIT   = -99.d0
      real(r8) :: QN = -99.d0
!
!  DOM parameters.
!
      real(r8), dimension(2) :: cDOCfrac_c = -99.d0
      real(r8) :: RtUVR_DIC = -99.d0
      real(r8) :: RtUVR_DOC = -99.d0
!
!  Fecal parameters.
!
      real(r8), dimension(Nfec) :: WF = -99.d0
      real(r8), dimension(Nfec) :: RegTbase = -99.d0
      real(r8), dimension(Nfec) :: RegTfac = -99.d0
      real(r8), dimension(Nfec) :: RegCR = -99.d0
      real(r8), dimension(Nfec) :: RegNR = -99.d0
      real(r8), dimension(Nfec) :: RegSR = -99.d0
      real(r8), dimension(Nfec) :: RegPR = -99.d0
      real(r8), dimension(Nfec) :: RegFR = -99.d0
!
! Marta Rodrigues
!  Zooplankton parameters.
      
      real(r8), dimension(Nzoo) :: ZooDOC = -99.d0
      real(r8), dimension(Nzoo,2) :: ZooPEL = -99.d0
      real(r8), dimension(Nzoo) :: ZooCYC = -99.d0
      real(r8), dimension(Nzoo,Nphy) :: DeltaZoo = -99.d0
      real(r8), dimension(Nzoo,Nphy) :: EfcCap = -99.d0
      real(r8), dimension(Nzoo) :: HsZoo = -99.d0
      real(r8), dimension(Nzoo) :: EfcPrd = -99.d0
      real(r8), dimension(Nzoo) :: ExZOO = -99.d0
      real(r8), dimension(Nzoo) :: GZ = -99.d0
      real(r8), dimension(Nzoo) :: basalZoo = -99.d0
      real(r8), dimension(Nzoo) :: QZoo = -99.d0
      real(r8), dimension(Nzoo) :: etaZoo = -99.d0

!  Oxygen parameters.

      real(r8) :: omegaO2C = -99.d0
      real(r8) :: omegaO2N = -99.d0
      real(r8) :: omegaO2NDenit  = -99.d0
      real(r8) :: omegaS2O  = -99.d0
      real(r8) :: QWind = -99.d0
      real(r8) :: psiWind = -99.d0
      real(r8) :: alfaWind = -99.d0
      real(r8) :: RtDenit = -99.d0
      real(r8) :: MDenit = -99.d0
      real(r8) :: reox_COD = -99.d0
      real(r8) :: QCOD = -99.d0
      real(r8) :: HsCOD = -99.d0
      real(r8) :: pCO2a = -99.d0
      real(r8) :: pH = -99.d0   

      real(r8), allocatable :: specir(:,:)
      real(r8), allocatable :: avcos(:,:)

! MFR - Nov/2015 - Atmospheric variable for specir.F90 
      real(r8), allocatable :: cloud(:)  ! cloud cover @ elements
      real(r8), allocatable :: hair(:)   ! relative humidity @ elements
      real(r8), allocatable :: tair(:)   ! air temperature @ elements
      real(r8), allocatable :: pair(:)   ! atm. pressures @ elements
      real(r8), allocatable :: uwind(:)  ! U-wind @ elements
      real(r8), allocatable :: vwind(:)  ! V-wind @ elements
     

!-----------------------------------------------------------------------
!  Internal parameters.
!-----------------------------------------------------------------------
!    
!  Spectral band width used in light calculations.

      real(r8), parameter :: DLAM  = 5.0_r8
!
!  Flags used for testing purposes.
!
      real(r8), parameter :: SMALL  = 1.0e-6_r8
      real(r8), parameter :: VSMALL = 1.0e-14_r8
      real(r8), parameter :: LARGE  = 1.0e+10_r8
      real(r8), parameter :: VLARGE = 1.0e+20_r8
!
!  Array indexes for frequently used constituents.
!
      integer, parameter :: ilab=1    ! labile index for DOC.
      integer, parameter :: irct=2    ! relict index for DOC.
      integer, parameter :: ichl=1    ! pigment index for chlorophyll-a
      integer, parameter :: isfc=1    ! slow fecal pellet index
      integer, parameter :: iffc=2    ! fast fecal pellet index
!
!  Phytoplankton calculated paramters.
!
      real(r8), dimension(Nphy) :: ImaxC2nALG = -99.d0   ! inverse C2nALG
      real(r8), dimension(Nphy) :: ImaxC2SiALG = -99.d0  ! inverse C2SiALG
      real(r8), dimension(Nphy) :: ImaxC2pALG = -99.d0   ! inverse C2pALG
      real(r8), dimension(Nphy) :: ImaxC2FeALG = -99.d0  ! inverse C2FeALG
!
!  Bacteria calculated parameters.
!
      real(r8) :: N2cBAC = -99.d0
      real(r8) :: P2cBAC = -99.d0
      real(r8) :: Fe2cBAC = -99.d0
      real(r8), dimension(Nbac) :: HsNH4_ba = -99.d0
      real(r8), dimension(Nbac) :: HsPO4_ba = -99.d0
      real(r8), dimension(Nbac) :: HsFe_ba = -99.d0
      real(r8) :: R_ExBAC_c = -99.d0
      real(r8) :: ExBAC_n = -99.d0
      real(r8) :: Frac_ExBAC_n = -99.d0
      real(r8) :: I_Bac_Ceff = -99.d0
!
!  Absorption parameters.
!
      real(r8), dimension(NBands) :: wavedp = -99.d0   ! a and b factor
      real(r8), dimension(2) :: aDOC410 = -99.d0    ! CDM absorption at 410
      real(r8), dimension(2) :: aDOC300 = -99.d0   ! CDM absorption at 300
      
! Marta Rodrigues
!  Zooplankton calculated parameters.
!  Keep the same name
!
!
! Other parameters

      real(r8), parameter :: sec2day=1.0_r8/86400.0_r8       

      CONTAINS
!=======================================================================
      SUBROUTINE initialize_biology
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!================================================== W. Paul Bissett ====
!                                                                      !
!  This routine initializes several parameters in module "mod_biology" !
!  for all nested grids.                                               !
!                                                                      !
! Marta Rodrigues                                                      !
!  "mod_biology" - called "biology"                                    !
!=======================================================================
      USE bio_param
      USE eclight
      USE schism_glbl, only:wsett,irange_tr
 
!
!  Local variable declarations
!
      integer :: ibac, iband, ifec, iphy, ng
! Marta Rodrigues
      integer :: izoo
      integer :: tmp00, tmp1 
!
!-----------------------------------------------------------------------
!  Derived parameters.
!-----------------------------------------------------------------------
   ibac=0

!  Convert rates from day-1 to second-1.
!
        tmp00 = irange_tr(1,6)-1

! Marta Rodrigues - DO ng=1,Ngrids
        DO iphy=1,Nphy
          GtALG_max(iphy)=GtALG_max(iphy)*sec2day
          ExALG(iphy)=ExALG(iphy)*sec2day
          HsGRZ(iphy)=HsGRZ(iphy)*sec2day
          
          !WS(iphy)=WS(iphy)*sec2day
          ! MFR, Nov/2015 - Set wsett for phytoplankton tracers
          tmp1 = tmp00+iPhyC(iphy)
          wsett(tmp1,:,:)=WS(iphy)*sec2day
          tmp1 = tmp00+iPhyN(iphy)
          wsett(tmp1,:,:)=WS(iphy)*sec2day
          tmp1 = tmp00+iPhyP(iphy)
          wsett(tmp1,:,:)=WS(iphy)*sec2day
          IF(IRON==1)THEN
             tmp1 = tmp00+iPhyF(iphy)
             wsett(tmp1,:,:)=WS(iphy)*sec2day
          ENDIF
          IF (iPhyS(iphy).gt.0) THEN
              tmp1 = tmp00+iPhyS(iphy)
              wsett(tmp1,:,:)=WS(iphy)*sec2day
          ENDIF 
  
          tmp1 = tmp00+iPigs(iphy,1)
          wsett(tmp1,:,:)=WS(iphy)*sec2day
          IF (iPigs(iphy,2).gt.0) THEN
             tmp1 = tmp00+iPigs(iphy,2)
             wsett(tmp1,:,:)=WS(iphy)*sec2day
          ENDIF
          IF (iPigs(iphy,3).gt.0) THEN
             tmp1 = tmp00+iPigs(iphy,3)
             wsett(tmp1,:,:)=WS(iphy)*sec2day
          ENDIF
          IF (iPigs(iphy,4).gt.0) THEN
             tmp1 = tmp00+iPigs(iphy,4)
             wsett(tmp1,:,:)=WS(iphy)*sec2day
          ENDIF
          IF (iPigs(iphy,5).gt.0) THEN
             tmp1 = tmp00+iPigs(iphy,5)
             wsett(tmp1,:,:)=WS(iphy)*sec2day
          ENDIF
          IF (iPigs(iphy,6).gt.0) THEN
             tmp1 = tmp00+iPigs(iphy,6)
             wsett(tmp1,:,:)=WS(iphy)*sec2day
          ENDIF
          IF (iPigs(iphy,7).gt.0) THEN
             tmp1 = tmp00+iPigs(iphy,7)
             wsett(tmp1,:,:)=WS(iphy)*sec2day
          ENDIF  
          ! End wsett for phytoplankton tracers

          basalPhy(iphy)=basalPhy(iphy)*sec2day 
        END DO
        DO ibac=1,Nbac
          GtBAC_max(ibac)=GtBAC_max(ibac)*sec2day
          basalBac(ibac)=basalBac(ibac)*sec2day
        END DO
        DO ifec=1,Nfec
          ! WF(ifec)=WF(ifec)*sec2day
          ! MFR, Nov/2015 - Set wsett for fecal/detritus tracers
          tmp1 = tmp00+iFecC(ifec)
          wsett(tmp1,:,:)=WF(ifec)*sec2day
          tmp1 = tmp00+iFecN(ifec)
          wsett(tmp1,:,:)=WF(ifec)*sec2day
          tmp1 = tmp00+iFecP(ifec)
          wsett(tmp1,:,:)=WF(ifec)*sec2day
          IF(IRON==1) THEN
            tmp1 = tmp00+iFecF(ifec)
            wsett(tmp1,:,:)=WF(ifec)*sec2day
          ENDIF 
          tmp1 = tmp00+iFecS(ifec)
          wsett(tmp1,:,:)=WF(ifec)*sec2day
          ! End wsett for fecal/detritus tracers
        END DO
        RtNIT=RtNIT*sec2day
        RtDENIT=RtDENIT*sec2day
! Marta Rodrigues	
	DO izoo=1,Nzoo
	  ExZOO(izoo)= ExZOO(izoo)*sec2day
	  GZ(izoo)=GZ(izoo)*sec2day
          basalZoo(izoo)=basalZoo(izoo)*sec2day
	 END DO 
!      END DO
       reox_COD=reox_COD*sec2day

!  Calculated reciprocal phytoplankton parameters.
!
! Marta Rodrigues - DO ng=1,Ngrids
        DO iphy=1,Nphy
          IF (maxC2nALG(iphy).gt.SMALL) THEN
            ImaxC2nALG(iphy)=1.0_r8/maxC2nALG(iphy)
          ELSE
            ImaxC2nALG(iphy)=0.0_r8
          END IF
          IF (maxC2SiALG(iphy).gt.SMALL) THEN
            ImaxC2SiALG(iphy)=1.0_r8/maxC2SiALG(iphy)
          ELSE
            ImaxC2SiALG(iphy)=0.0_r8
          END IF
          IF (maxC2pALG(iphy).gt.SMALL) THEN
            ImaxC2pALG(iphy)=1.0_r8/maxC2pALG(iphy)
          ELSE
            ImaxC2pALG(iphy)=0.0_r8
          END IF
          IF (maxC2FeALG(iphy).gt.SMALL) THEN
            ImaxC2FeALG(iphy)=1.0_r8/maxC2FeALG(iphy)
          ELSE
            ImaxC2FeALG(iphy)=0.0_r8
          END IF
        END DO
!     END DO
!
!  Calculated bacterial parameters.
!
! Marta Rodrigues - DO ng=1,Ngrids
          DO ibac=1,Nbac
            HsNH4_ba(ibac)=HsDOC_ba(ibac)/C2nBAC
            HsPO4_ba(ibac)=HsDOC_ba(ibac)/C2pBAC
            HsFe_ba (ibac)=HsDOC_ba(ibac)/C2FeBAC
          END DO
!     END DO
!
!  Inverse parameters for computational efficiency.
!
! Marta Rodrigues - DO ng=1,Ngrids
        N2cBAC=1.0_r8/C2nBAC
        P2cBAC=1.0_r8/C2pBAC
        Fe2cBAC=1.0_r8/C2FeBAC
        I_Bac_Ceff=1.0_r8/Bac_Ceff 
!     END DO
!
!  Reciprocal of non baterial recalcitran carbon excretion.
!
! Marta Rodrigues- DO ng=1,Ngrids
        R_ExBAC_c=1.0_r8/(1.0_r8-ExBAC_c)
!     END DO
!
!  Bacterial recalcitrant nitrogen excretion as a function of uptake.
!
! Marta Rodrigues - DO ng=1,Ngrids
        ExBAC_n=ExBAC_c*C2nBAC/ExBacC2N
        Frac_ExBAC_n=1.0_r8-ExBAC_n
!     END DO
!
!  Scale UV degradation parameters.
!
! Marta Rodrigues - DO ng=1,Ngrids
        RtUVR_DIC=RtUVR_DIC/3600.0_r8
        RtUVR_DOC=RtUVR_DOC/3600.0_r8 
!     END DO
!
!  If applicable, zero-out fecal regeneration parameters.
!
! Marta Rodrigues - DO ng=1,Ngrids
        IF (Regen_flag) THEN
          DO ifec=1,Nfec
            RegCR(ifec)=RegCR(ifec)*sec2day
            RegNR(ifec)=RegNR(ifec)*sec2day
            RegPR(ifec)=RegPR(ifec)*sec2day
            RegFR(ifec)=RegFR(ifec)*sec2day
            RegSR(ifec)=RegSR(ifec)*sec2day
          END DO
        ELSE
          DO ifec=1,Nfec
            RegCR(ifec)=0.0_r8
            RegNR(ifec)=0.0_r8
            RegPR(ifec)=0.0_r8
            RegFR(ifec)=0.0_r8
            RegSR(ifec)=0.0_r8
          END DO
        END IF
!     END DO
!
!  Spectral dependency for scattering and backscattering.
!
      DO iband=1,NBands
        wavedp(iband)=(550.0_r8/(397.0_r8+REAL(iband,r8)*DLAM))
      END DO
!
!  Calculated IOP parameter values.
!     
      aDOC =reshape( aDOCall, (/ 2, NBands /), ORDER=(/ 2, 1 /) )
      aDOC410(ilab)=aDOC(ilab,1)*EXP(0.014_r8*(ec_wave_ab(1)-410.0_r8))
      IF(Ndom==2) aDOC410(irct)=aDOC(irct,1)*EXP(0.025_r8*(ec_wave_ab(1)-410.0_r8))
      aDOC300(ilab)=EXP(0.0145_r8*(410.0_r8-300.0_r8))
      IF(Ndom==2) aDOC300(irct)=EXP(0.0145_r8*(410.0_r8-300.0_r8))

      RETURN
      END SUBROUTINE initialize_biology

      END MODULE biology
