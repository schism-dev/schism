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

module icm_sed_mod
!-------------------------------------------------------------------------------
!parameter definition for sediment flux model 
!-------------------------------------------------------------------------------
  use schism_glbl,only: iwp,nea
  implicit none

  !Hydro
  integer, save :: SED_IWC        !SED_IWC=id
  real(kind=iwp), save,allocatable,dimension(:) :: SFA
  real(kind=iwp), save,allocatable,dimension(:) :: SED_BL
  real(kind=iwp), save,allocatable,dimension(:) :: ZD

  !bottom layer concentration from water column
  real(kind=iwp), save,allocatable,dimension(:,:) :: SED_B !3 phyto. species
  real(kind=iwp), save,allocatable,dimension(:) :: SED_LPOP,SED_RPOP,SED_LPON,SED_RPON,SED_LPOC,SED_RPOC,SED_TSS
  real(kind=iwp), save,allocatable,dimension(:) :: SED_SU,SED_PO4,SED_NH4,SED_NO3,SED_SA,SED_DO,SED_COD,SED_SALT,SED_T,SSI
  
  !steady state
  integer, save :: iSteady
  real(kind=iwp), save :: TINTIM
  real(kind=iwp), save,allocatable,dimension(:) :: AG3CFL,AG3NFL,AG3PFL,ASDTMP

  !General Parameters
  real(kind=iwp),save :: HSEDALL,DIFFT,SALTSW,SALTND
  real(kind=iwp),save,dimension(3,3) :: FRPPH,FRNPH,FRCPH 
  real(kind=iwp),save,dimension(3) :: FRPPHB,FRNPHB,FRCPHB
  real(kind=iwp),save,dimension(3) :: KPDIAG,KNDIAG,KCDIAG,DPTHTA,DNTHTA,DCTHTA
  real(kind=iwp),save :: W2,H2,m1,m2,KSI,THTASI,THTADP,THTADD
  !ncai !sav related
  real(kind=iwp),save,dimension(3) :: frnsav,frpsav,frcsav


  !nutrients, parameters
  real(kind=iwp),save :: KAPPNH4F,KAPPNH4S,PIENH4,THTANH4,KMNH4,KMNH4O2 !!nitrification
  real(kind=iwp),save :: KAPPNO3F,KAPPNO3S,K2NO3,THTANO3 !denitrification
  real(kind=iwp),save :: KAPPD1,KAPPP1,PIE1S,PIE2S,THTAPD1,KMHSO2 !H2S
  real(kind=iwp),save :: CSISAT,DPIE1SI,PIE2SI,KMPSI,O2CRITSI,JSIDETR  !Si
  real(kind=iwp),save :: DPIE1PO4F,DPIE1PO4S,PIE2PO4,O2CRIT,KMO2DP  !PO4
  real(kind=iwp),save :: TEMPBEN,KBENSTR,KLBNTH,DPMIN  !benthic stress
  real(kind=iwp),save :: KAPPCH4,THTACH4,KMCH4O2,KMSO4 !CH4 reaction

  !initial concentration 
  real(kind=iwp),save :: CTEMPI,BBMI,CPOSI,PO4T2I,NH4T2I,NO3T2I,HST2I,CH4T2I,CH41TI,SO4T2I,SIT2I,BENSTI
  real(kind=iwp),save,dimension(3) :: CPOPI,CPONI,CPOCI

  !Sediment thickness, burial and mixing rates
  real(kind=iwp),save,allocatable,dimension(:) :: HSED,VSED,VPMIX,VDMIX
  real(kind=iwp),save :: W12,W12MIN,KL12

  !splits of refractory POM from Water Column into G2,G3 class POM in sediment
  real(kind=iwp), save,allocatable,dimension(:,:) :: FRPOP,FRPON,FRPOC !(1:nea,3), leave option for mapping
  
  !POM fluxes !unit:mg/m^2
  real(kind=iwp), save,allocatable,dimension(:,:) :: flxpop,flxpon,flxpoc
  real(kind=iwp), save,allocatable,dimension(:) :: flxpos

  !sediment concentration !unit:mg/m^3
  real(kind=iwp), save,allocatable,dimension(:) :: CTEMP,CPIP,CNO3,CNH4,CCH4,CSO4,CPOS,CH2S
  real(kind=iwp), save,allocatable,dimension(:,:) :: CPOP,CPON,CPOC
  real(kind=iwp), save,allocatable,dimension(:) :: CH4T2TM1S,CH41TM1S,SO4T2TM1S,BENSTR1S,BFORMAXS,ISWBENS
  real(kind=iwp), save,allocatable,dimension(:) :: POP1TM1S,POP2TM1S,POP3TM1S,PON1TM1S,PON2TM1S,PON3TM1S,POC1TM1S,POC2TM1S,POC3TM1S,PSITM1S
  real(kind=iwp), save,allocatable,dimension(:) :: NH41TM1S,NO31TM1S,HS1TM1S,SI1TM1S,PO41TM1S,NH4T2TM1S,NO3T2TM1S,HST2TM1S,SIT2TM1S,PO4T2TM1S
  real(kind=iwp), save,allocatable,dimension(:) :: DFEEDM1S !deposit feeder

  real(kind=iwp),save :: TEMPD,PON1,PON2,PON3,POC2,POC3,POP1,POP2,POP3,PSI
  real(kind=iwp),save :: NH41TM1,NO31TM1,HS1TM1,SI1TM1,PO41TM1,NH4T2TM1,NO3T2TM1,HST2TM1,SIT2TM1,PO4T2TM1,CH4T2TM1,CH41TM1,SO4T2TM1
  real(kind=iwp),save :: PON1TM1,PON2TM1,PON3TM1,POC1TM1,POC1,POC2TM1,POC3TM1,POP1TM1,POP2TM1,POP3TM1,PSITM1
  real(kind=iwp),save :: ROOTDO  !SAV
  real(kind=iwp),save :: DFEED,DFEEDM1 !deposit feeder 
  real(kind=iwp),save :: BENSTR1 !benthic stress

  real(kind=iwp),save :: SI1,SI2,SIT1,SIT2,PO41,PO42,PO4T1,PO4T2,NH41,NH42,NH4T1,NH4T2
  real(kind=iwp),save :: NO31,NO32,NO3T1,NO3T2,HS1,HS2,HST1,HST2,CH41,CH42,CH4T1,CH4T2,SO41,SO42,SO4T1,SO4T2

  !diagenesis fluxes
  real(kind=iwp),save :: XJN,XJC,XJP
  
  !sediment fluxes
  real(kind=iwp),save :: JSI,JPO4,JNH4,JNO3,JHS,JCH4,JCH4AQ,JCH4G,JN2GAS,JGAS

  !nutrient concentration in water column
  real(kind=iwp),save :: PO40,NH40,NO30,SI0,O20,HS0,SAL0,SO40MG
  real(kind=iwp),save :: CH4SAT

  !reaction rate (temp vars)
  real(kind=iwp),save :: TEMP5,TEMP20,TEMP202
  real(kind=iwp),save :: ZHTANH4F,ZHTANH4S,ZHTAD1,ZHTAP1,ZHTANO3F,ZHTANO3S,ZHTA2NO3,ZL12NOM,ZW12NOM,ZHTAPON1,ZHTAPON2,ZHTAPON3
  real(kind=iwp),save :: ZHTAPOC1,ZHTAPOC2,ZHTAPOC3,ZHTAPOP1,ZHTAPOP2,ZHTAPOP3,ZHTASI,ZHTACH4,ZHTAI0,ZHTAR,ZHTABETA

  !benthic stress
  real(kind=iwp),save :: BFORMAX,ISWBEN,BFOR,BENSTR

  !SOD calculation
  real(kind=iwp),save :: SOD,stc

  !diffusion under hypoxia
  real(kind=iwp),save :: O2CRITdif,stc0
  real(kind=iwp),save :: thtaTdif,alphaTdif

  !sediment fluxes
  real(kind=iwp),save,allocatable,dimension(:) :: SED_BENDO,SED_BENCOD,SED_BENNH4,SED_BENNO3,SED_BENPO4,SED_BENDOC,SED_BENSA

  !erosion fluxes
  real(kind=iwp),save,allocatable,dimension(:) :: SED_EROH2S,SED_EROLPOC,SED_ERORPOC !nea
  real(kind=iwp),save,allocatable,dimension(:) :: tau_c_elem !nea
  real(kind=iwp),save :: eroporo,erorate,erofrac,erodiso !0.9; 0.01kg/m^2/s; 80% in mud, 20% in sand
  real(kind=iwp),save :: depofracR,depofracL,depoWSR,depoWSL
  integer, save :: iERO,iDEPO

  !bottom Light (nea)
  real(kind=iwp),save,allocatable,dimension(:) :: sbLight
  
  !deposit feeder
  integer,save :: idf,ihypox
  real(kind=iwp),save :: XKMI0,ING0,THTAI0,R,THTAR,BETA,THBETA
  real(kind=iwp),save :: AMCN,AMCP,AA1,AA2,XKMG1,XKMG2
  real(kind=iwp),save :: XKBO2,TDD,DOLOW,DFDOH,DFDOQ,RDD,RMORT
  real(kind=iwp),save :: XKI0,XKR,XKBETA
  real(kind=iwp),save,allocatable,dimension(:) :: SFLUXP,SF_RPOP,SFLUXN,SF_RPON,SFLUXC,SF_RPOC,JSUSF,SF_SU
  
  !benthic algae
  integer, save :: iBalg
  real(kind=iwp),save :: PO4AVL, NH4AVL, NO3AVL
  real(kind=iwp),save :: PMB,ANCB,APCB,KTGB1,KTGB2,TMB
  real(kind=iwp),save :: ALPHB,CCHLB,KESED,KEBALG,KHNB,KHPB,KHRB
  real(kind=iwp),save :: BMRB,BPRB,KTBB,TRB,BALGMIN
  real(kind=iwp),save :: FNIB,FPIB
  real(kind=iwp),save :: FTB,PRNB
  real(kind=iwp),save,allocatable,dimension(:) :: BBM
  real(kind=iwp),save :: FIB,BLITE,NLB,PLB,BMB,PB,NPPB,PRB
  real(kind=iwp),save :: BAPOC,BAPON,BAPOP

end module icm_sed_mod
