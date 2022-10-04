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

module icm_mod
  use schism_glbl,only: rkind,nvrt,nea
  implicit none

  !-------------------------------------------------------------------------------
  !constants: molar weight for C,Ca,CaCo3,N
  !C1_PAR: srad to PAR; C2_PAR: PAR(W/m2) to PAR(E/m2/day) 
  !-------------------------------------------------------------------------------
  integer,parameter :: nPB=3,nZB=2
  real(rkind),parameter :: mC=12.011,mCACO3=100.086,mN=14.007
  real(rkind),parameter :: C1_PAR=0.47, C2_PAR=0.397

  !-------------------------------------------------------------------------------
  !global switch and variables
  !-------------------------------------------------------------------------------
  integer,save,target :: nsub,iKe,iLight,iPR,iLimit,iSed,iBA,iRad,isflux,ibflux
  integer,save,target :: iSilica,iZB,iPh,iCBP,isav_icm,iveg_icm,idry_icm
  real(rkind),save,target :: KeC,KeS,KeSalt,Ke0,tss2c,PRR(3),wqc0(29),WSP(29),WSPn(29)
  real(rkind),save,target,dimension(3) :: alpha
  integer,save,pointer :: jdry,jsav,jveg,ised_icm,iBA_icm

  integer,parameter :: nout_sav=4, nout_veg=12, nout_sed=26, nout_ba=1
  integer,save,target :: ntrs_icm,itrs(2,9),nout_icm,nout_d2d,nout_d3d,n2d(9),n3d(9),i2d(9),i3d(9)
  integer,save,pointer :: itrs_icm(:,:),elem_in(:,:)
  integer,save :: iPB1,iPB2,iPB3,iRPOC,iLPOC,iDOC,iRPON,iLPON,iDON,iNH4,iNO3,iRPOP,iLPOP, &
        & iDOP,iPO4,iCOD,iDOX,iSU,iSA,iZB1,iZB2,iTIC,iALK,iCA,iCACO3,iSRPOC,iSRPON,iSRPOP,iPIP
  character(len=6),save :: name_icm(100),name_d2d(100),name_d3d(100)
  integer,save,target :: ncid_icm(3),npt_icm(3)
  real(rkind),target,save :: time_icm(2,3),dt_icm(3)
  real(rkind),target,save,allocatable :: rad_in(:,:),sflux_in(:,:,:),bflux_in(:,:,:),wqc_d2d(:,:),wqc_d3d(:,:,:)

  !declear temporary variables to increase code readability (can be put in main loop)
  real(rkind),save,pointer,dimension(:,:) :: wqc,ZBS,PBS 
  real(rkind),save,pointer,dimension(:) :: temp,salt,ZB1,ZB2,PB1,PB2,PB3,RPOC,LPOC,DOC,RPON,LPON,DON,NH4, &
                              & NO3,RPOP,LPOP,DOP,PO4,SU,SA,COD,DOX,TIC,ALK,CA,CACO3,SRPOC,SRPON,SRPOP,PIP
  real(rkind),save,target,allocatable :: DIN(:),dwqc(:,:),zdwqc(:,:),sdwqc(:,:),vdwqc(:,:),gdwqc(:,:) 
  real(rkind),save,pointer,dimension(:,:) :: zdPBS,zdC,zdN,zdP,zdS
  real(rkind),save,pointer,dimension(:) :: zdDOX 
 
  !-------------------------------------------------------------------------------
  !ICM parameters and variables
  !-------------------------------------------------------------------------------
  real(rkind),save,target,dimension(3) :: GPM,TGP,MTB,TMT,KTMT,MTR
  real(rkind),save,target :: KTGP(3,2)
  real(rkind),save,target :: FCP(3,4),FNP(3,5),FPP(3,5),FCM(3,4),FNM(3,5),FPM(3,5)
  real(rkind),save,target :: Nit,TNit,KTNit(2),KhDOn,KhNH4n,KhDOox,KhNO3dn
  real(rkind),save,target,dimension(3) :: KC0,KN0,KP0,KCalg,KNalg,KPalg,TRM,KTRM,KSR0,TRSR,KTRSR
  real(rkind),save,target :: KCD,TRCOD,KTRCOD,KhCOD,KPIP
  real(rkind),save,target,dimension(3) :: KhN,KhP,KhSal,c2chl,n2c,p2c,KhDO,PBmin
  real(rkind),save,target :: o2c,o2n,dn2c,an2c,KPO4p,WRea,dz_flux(2)

  real(rkind),save :: dtw     !ICM time step (day)
  real(rkind),save:: time_ph  !time stamp for WQinput
  real(rkind),save :: rIa,rIavg

  !-------------------------------------------------------------------------------
  !silica parameters and variables
  !-------------------------------------------------------------------------------
  real(rkind),save,target :: FSP(2),FSM(2),KS,TRS,KTRS,KhS(3),s2c(3),KSAp

  !-------------------------------------------------------------------------------
  !zooplankton parameters and variables
  !-------------------------------------------------------------------------------
  real(rkind),save,target :: zGPM(8,2),zKhG(8,2),zTGP(2),zKTGP(2,2),zAG,zRG
  real(rkind),save,target,dimension(2) :: zMRT,zMTB,zTMT,zKTMT
  real(rkind),save,target :: zFCP(3),zFNP(4),zFPP(4),zFSP(2)
  real(rkind),save,target :: zFCM(2),zFNM(2,4),zFPM(2,4),zFSM(2,2)
  real(rkind),save,target :: zs2c(2),zn2c(2),zp2c(2),zKhDO(2),z2pr(2),p2pr

  !-------------------------------------------------------------------------------
  !pH parameters and variables
  !-------------------------------------------------------------------------------
  integer,save,target :: inu_ph
  real(rkind),save,target :: ppatch0,pKCACO3,pKCA,pRea

  integer, save :: irec_ph
  integer,save,allocatable :: ppatch(:)
  real(rkind),save,allocatable :: ph_nudge(:),ph_nudge_nd(:)
  real(rkind),save,allocatable,dimension(:,:) :: TIC_el,ALK_el

  !-------------------------------------------------------------------------------
  !SAV parameters and variables
  !-------------------------------------------------------------------------------
  real(rkind),save,target :: spatch0,stleaf0,ststem0,stroot0
  real(rkind),save,target :: sGPM,sTGP,sKTGP(2),sFAM,sFCP(3) !growth related coefficients
  real(rkind),save,target :: sMTB(3),sTMT(3),sKTMT(3)        !meta. coefficients (rate, temp, temp dependence)
  real(rkind),save,target :: sFCM(4),sFNM(4),sFPM(4)         !metabolism to (RPOM,RLOM,DOM,DIM)
  real(rkind),save,target :: sKhNw,sKhNs,sKhNH4,sKhPw,sKhPs  !half-saturation conc. of N&P
  real(rkind),save,target :: salpha,sKe,shtm(2),s2ht(3)      !(P-I curve, light attenu., canopy height)
  real(rkind),save,target :: sc2dw,s2den,sn2c,sp2c,so2c      !convert ratios

  integer,save,allocatable :: spatch(:)               !sav region
  real(rkind),save,allocatable,dimension(:) :: stleaf,ststem,stroot,sht
  real(rkind),save,allocatable,dimension(:,:) :: sleaf,sstem,sroot !(nvrt,nea), unit: g/m^2
  real(rkind),save,allocatable,dimension(:) :: sroot_POC,sroot_PON,sroot_POP,sroot_DOX !(nea), unit: g/m^2/day
  real(rkind),save,allocatable,dimension(:) :: sleaf_NH4,sleaf_PO4  !(nea), unit: g/m^2/day

  !-------------------------------------------------------------------------------
  !VEG parameters and variables
  !-------------------------------------------------------------------------------
  real(rkind),save,target,dimension(3) :: vtleaf0,vtstem0,vtroot0  !init conc.
  real(rkind),save,target :: vpatch0,vGPM(3),vFAM(3),vTGP(3),vKTGP(3,2),vFCP(3,3) !growth related coefficients
  real(rkind),save,target,dimension(3,3) :: vMTB,vTMT,vKTMT    !meta. coefficients (rate,temp,temp dependence)
  real(rkind),save,target,dimension(3,4) :: vFNM,vFPM,vFCM     !metabolism to (RPOM,RLOM,DOM,DIM)
  real(rkind),save,target,dimension(3) :: vKhNs,vKhPs,vScr,vSopt,vInun,valpha,vKe !growth limit(nutrent,light,salinity,inundation)
  real(rkind),save,target,dimension(3,2) :: vTMR,vKTMR,vMR0,vMRcr    !mortality coeffs
  real(rkind),save,target :: v2ht(3,2),vht0(3),vcrit(3)    !computing canopy height
  real(rkind),save,target,dimension(3) :: vc2dw,v2den,vn2c,vp2c,vo2c!convert ratios
  integer,save,target :: ivNc,ivPc,ivNs,ivPs,ivMRT              !flags for (N,P) limit, recycled (N,P) dest., mortality

  real(rkind),save :: mtemp !todo add function to read mtemp
  integer,save,allocatable :: vpatch(:)                     !reg region
  real(rkind),save,allocatable,dimension(:,:) :: vht !,ztcveg !(nea,3)
  real(rkind),save,allocatable,dimension(:,:) :: vtleaf,vtstem,vtroot !(nea,3)
  real(rkind),save,allocatable,dimension(:,:) :: vroot_POC,vroot_PON,vroot_POP,vroot_DOX !(nea,3)
  real(rkind),save,allocatable,dimension(:,:) :: vleaf_NH4,vleaf_PO4 !(nea,3)

  !-------------------------------------------------------------------------------
  !sediment flux model (SFM) parameters and variables
  !-------------------------------------------------------------------------------
  real(rkind),save,target :: bdz,bVb,bdiff,bsaltc,bsaltp,bsaltn,bsolid(2)
  real(rkind),save,target :: bKTVp,bKTVd,bVp,bVd,bTR
  real(rkind),save,target :: btemp0,bstc0,bSTR0,bThp0,bTox0,bNH40,bNO30,bPO40,bH2S0
  real(rkind),save,target :: bCH40,bPOS0,bSA0,bPOP0(3),bPON0(3),bPOC0(3)
  real(rkind),save,target,dimension(3) :: bKC,bKN,bKP,bKTC,bKTN,bKTP
  real(rkind),save,target :: bKS,bKTS
  real(rkind),save,target,dimension(3,3) :: bFPP,bFNP,bFCP,bFCv,bFNv,bFPv !(G1:G3,veg/PB)
  real(rkind),save,target,dimension(3) :: bFCs,bFNs,bFPs,bFCM,bFNM,bFPM !(G1:G3)
  real(rkind),save,target :: bKNH4f,bKNH4s,bpieNH4,bKTNH4,bKhNH4,bKhDO_NH4 !!nitrification
  real(rkind),save,target :: bKNO3f,bKNO3s,bKNO3,bKTNO3 !denitrification
  real(rkind),save,target :: bKH2Sd,bKH2Sp,bpieH2Ss,bpieH2Sb,bKTH2S,bKhDO_H2S !H2S oxidation
  real(rkind),save,target :: bSIsat,bKOSI,bpieSI,bKhPOS,bDOc_SI,bJPOSa  !Silica dissolution
  real(rkind),save,target :: bKOPO4f,bKOPO4s,bpiePO4,bDOc_PO4  !PO4
  real(rkind),save,target :: bKST,bSTmax,bp2d,bVpmin,bKhDO_Vp,bDOc_ST,banoxic,boxic !benthic stress
  real(rkind),save,target :: bKCH4,bKTCH4,bKhDO_CH4,bo2n !CH4 reaction

  !sediment concentrations and fluxes
  real(rkind),save,target,allocatable,dimension(:) :: bLight,bThp,bTox,btemp,bCH4,bSTR,bPOS
  real(rkind),save,target,allocatable,dimension(:) :: bNH4s,bNH4,bNO3,bH2S,bSA,bPO4,bstc
  real(rkind),save,target,allocatable,dimension(:,:) :: bPOC,bPON,bPOP
  real(rkind),save,target,allocatable,dimension(:) :: SOD,JNH4,JNO3,JPO4,JCOD,JSA

  !-------------------------------------------------------------------------------
  !Benthic Algea model (BA) parameters and variables
  !-------------------------------------------------------------------------------
  real(rkind),save,target :: gpatch0,BA0,gGPM,gTGP,gKTGP(2),gMTB,gPRR,gTR,gKTR,galpha
  real(rkind),save,target :: gKSED,gKBA,gKhN,gKhP,gp2c,gn2c,go2c,gFCP(3),gFNP(3),gFPP(3)

  integer,save,allocatable,dimension(:) :: gpatch
  real(rkind),save,target,allocatable,dimension(:) :: BA,gPR

  !-------------------------------------------------------------------------------
  !benthic erosion (ERO) parameters and variables
  !-------------------------------------------------------------------------------
  integer,save,target :: ierosion !0.9; 0.01kg/m^2/s; 80% in mud, 20% in sand 
  real(rkind),save,target :: erosion,etau,eporo,efrac,ediso,dfrac(2),dWS_POC(2)

  real(rkind),save,allocatable,dimension(:) :: eH2S,eLPOC,eRPOC !erosion flux

  !---------------------------------------------------------------------------
  !arguments for zbrent function
  !---------------------------------------------------------------------------
  type :: brent_var
    integer :: imed=0  !0: SFM;  1: pH
    integer :: id,ierr=0
    real(rkind) :: vmin,vmax
    real(rkind),pointer :: data 
    real(rkind) :: tdep,wsalt,Kd,Kp,wNH4,wNO3,wCOD,wDOX,XJC,XJN,SODrt,stc,SOD !SFM
    real(rkind) :: ph,K1,K2,Kw,Kb,Ct,Ca,Bt,rH !pH
  end type brent_var

  !---------------------------------------------------------------------------
  !spatially varying parameters: for different dimensions 
  !---------------------------------------------------------------------------
  type :: icm_spatial_param
    character(len=30) :: varname  !parameter name
    integer :: ndim=0  !parameter dimension
    integer :: dims(2) !dimension info
    real(rkind),dimension(30) :: data0 !oirginal value of data
    integer,allocatable,dimension(:,:) :: istat
    real(rkind),pointer :: p=>null()                !param. of single value
    real(rkind),pointer,dimension(:) :: p1=>null()   !param. of 1D array
    real(rkind),pointer,dimension(:,:) :: p2=>null() !param of 2D array
    real(rkind),allocatable,dimension(:,:,:) :: data
  end type icm_spatial_param
  type(icm_spatial_param),save,target,allocatable,dimension(:) :: sp

end module icm_mod

module icm_interface
  interface brent_def
    subroutine brent(bv)
      use icm_mod, only: brent_var
      type(brent_var),target,intent(inout) :: bv
    end subroutine brent
  end interface
end module icm_interface
