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
  integer,save,target :: nsub,iKe,iLight,iPR,iLimit,iSed,iBA,iRad,isflux,ibflux,iout_icm,nspool_icm
  integer,save,target :: iSilica,iZB,iPh,iCBP,isav_icm,imarsh_icm,nmarsh,idry_icm,iClam,nclam
  real(rkind),save,target :: KeC,KeS,KeSalt,Ke0,tss2c,PRR(3),wqc0(29),WSP(29),WSPn(29)
  real(rkind),save,target,dimension(3) :: alpha
  integer,save,pointer :: jdry,jsav,jmarsh,ised_icm,iBA_icm

  integer,save,target :: ntrs_icm,nout,nouts(10),iout(2,10),nout_icm_3d(2),nhot_icm
  integer,save,pointer :: itrs_icm(:,:),elem_in(:,:),nout_icm
  integer,save :: iPB1,iPB2,iPB3,iRPOC,iLPOC,iDOC,iRPON,iLPON,iDON,iNH4,iNO3,iRPOP,iLPOP, &
        & iDOP,iPO4,iCOD,iDOX,iSU,iSA,iZB1,iZB2,iTIC,iALK,iCA,iCACO3,iSRPOC,iSRPON,iSRPOP,iPIP
  integer,save,target :: iyear,imonth,iday,idoy
  character(len=6),save :: name_icm(100),name_d2d(100),name_d3d(100)
  integer,save,target :: ncid_icm(3),npt_icm(3)
  real(rkind),target,save :: time_icm(2,3),dt_icm(3)
  real(rkind),target,save,allocatable :: rad_in(:,:),sflux_in(:,:,:),bflux_in(:,:,:)

  !declear temporary variables to increase code readability (can be put in main loop)
  real(rkind),save,pointer,dimension(:,:) :: wqc,ZBS,PBS 
  real(rkind),save,pointer,dimension(:) :: salt,ZB1,ZB2,PB1,PB2,PB3,RPOC,LPOC,DOC,RPON,LPON,DON,NH4, &
                              & NO3,RPOP,LPOP,DOP,PO4,SU,SA,COD,DOX,TIC,ALK,CA,CACO3,SRPOC,SRPON,SRPOP,PIP
  real(rkind),save,target,allocatable :: temp(:),DIN(:),PO4d(:),TSS(:),dwqc(:,:),zdwqc(:,:),sdwqc(:,:), &
                              & vdwqc(:,:),gdwqc(:,:),cdwqc(:,:)
  real(rkind),save,pointer,dimension(:,:) :: zdPBS,zdC,zdN,zdP,zdS
  real(rkind),save,pointer,dimension(:) :: zdDOX 

  !debug variables
  real(rkind),save,pointer,dimension(:,:) :: dbCHLA
  real(rkind),save,pointer,dimension(:) :: dbTN,dbTP
 
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
  real(rkind),save,target :: spatch0,sav0(3)  !sav patch, and init conc.
  real(rkind),save,target :: sGPM,sTGP,sKTGP(2),sFAM,sFCP(3) !growth related coefficients
  real(rkind),save,target :: sMTB(3),sTMT(3),sKTMT(3)        !meta. coefficients (rate, temp, temp dependence)
  real(rkind),save,target :: sFCM(4),sFNM(4),sFPM(4)         !metabolism of leaf/stem to (RPOM,RLOM,DOM,DIM)
  real(rkind),save,target :: sFCMb(4),sFNMb(4),sFPMb(4)      !metabolism of root to (POM(G1-G3),dissolved nutrients)
  real(rkind),save,target :: sKhN(2),sKhP(2)                 !half-saturation conc. of N, P
  real(rkind),save,target :: salpha,sKe,shtm(2),s2ht(3)      !(P-I curve, light attenu., canopy height)
  real(rkind),save,target :: sc2dw,sn2c,sp2c,so2c      !convert ratios

  integer,save,allocatable :: spatch(:)               !sav region
  real(rkind),save,target,allocatable :: sht(:),sav(:,:),sFPOC(:,:),sFPON(:,:),sFPOP(:,:),sbNH4(:),sbPO4(:),sSOD(:)

  !-------------------------------------------------------------------------------
  !marsh parameters and variables
  !-------------------------------------------------------------------------------
  integer,save,target :: iNmarsh !swtich for N/P kentics
  real(rkind),save,target :: vpatch0 !init marsh region
  real(rkind),save,target,allocatable :: vmarsh0(:,:),vGPM(:),vFAM(:),vTGP(:),vKTGP(:,:),vFCP(:,:) !growth
  real(rkind),save,target,allocatable,dimension(:,:) :: vMTB,vTMT,vKTMT,vFCM,vFNM,vFPM !metabolsim
  !metabolism partition, growth limit of nutrient,light,salinity,inundation
  real(rkind),save,target,allocatable,dimension(:) :: vFW,vKhN,vKhP,valpha,vKe,vSopt,vKs,vInun
  real(rkind),save,target,allocatable :: vht0(:),vcrit(:),v2ht(:,:),vc2dw(:),vn2c(:),vp2c(:),vo2c(:) !misc

  integer,save,allocatable :: vpatch(:)  !marsh regions
  real(rkind),save,target,allocatable :: vmarsh(:,:,:),vht(:,:) !marsh biomass
  real(rkind),save,target,allocatable :: vFPOC(:,:),vFPON(:,:),vFPOP(:,:),vbNH4(:),vbPO4(:),vSOD(:) !sediment

  !-------------------------------------------------------------------------------
  !sediment flux model (SFM) parameters and variables
  !-------------------------------------------------------------------------------
  real(rkind),save,target :: bdz,bVb,bdiff,bsaltc,bsaltp,bsaltn,bsolid(2)
  real(rkind),save,target :: bKTVp,bKTVd,bVp,bVd,bTR
  real(rkind),save,target :: btemp0,bstc0,bSTR0,bThp0,bTox0,bNH40,bNO30,bPO40,bH2S0
  real(rkind),save,target :: bCH40,bPOS0,bSA0,bPOP0(3),bPON0(3),bPOC0(3)
  real(rkind),save,target,dimension(3) :: bKC,bKN,bKP,bKTC,bKTN,bKTP
  real(rkind),save,target :: bKS,bKTS
  real(rkind),save,target,dimension(3,3) :: bFPP,bFNP,bFCP!(G1:G3,PB)
  real(rkind),save,target,dimension(3) :: bFCM,bFNM,bFPM !(G1:G3)
  real(rkind),save,target :: bKNH4f,bKNH4s,bpieNH4,bKTNH4,bKhNH4,bKhDO_NH4 !!nitrification
  real(rkind),save,target :: bKNO3f,bKNO3s,bKNO3,bKTNO3 !denitrification
  real(rkind),save,target :: bKH2Sd,bKH2Sp,bpieH2Ss,bpieH2Sb,bKTH2S,bKhDO_H2S !H2S oxidation
  real(rkind),save,target :: bSIsat,bKOSI,bpieSI,bKhPOS,bDOc_SI,bJPOSa  !Silica dissolution
  real(rkind),save,target :: bKOPO4f,bKOPO4s,bpiePO4,bDOc_PO4  !PO4
  real(rkind),save,target :: bKST,bSTmax,bp2d,bVpmin,bKhDO_Vp,bDOc_ST,banoxic,boxic !benthic stress
  real(rkind),save,target :: bKCH4,bKTCH4,bKhDO_CH4,bo2n !CH4 reaction

  !sediment concentrations and fluxes
  real(rkind),save,target :: btau !sediment shear stress
  real(rkind),save,target,allocatable,dimension(:) :: bLight,bThp,bTox,btemp,bCH4,bSTR,bPOS
  real(rkind),save,target,allocatable,dimension(:) :: bNH4s,bNH4,bNO3,bH2S,bSA,bPO4,bstc
  real(rkind),save,target,allocatable,dimension(:,:) :: bPOC,bPON,bPOP
  real(rkind),save,target,allocatable,dimension(:) :: SOD,JNH4,JNO3,JPO4,JCOD,JSA

  !-------------------------------------------------------------------------------
  !Benthic Algea model (BA) parameters and variables
  !-------------------------------------------------------------------------------
  real(rkind),save,target :: gpatch0,gBA0,gGPM,gTGP,gKTGP(2),gMTB,gPRR,gTR,gKTR,galpha
  real(rkind),save,target :: gKSED,gKBA,gKhN,gKhP,gp2c,gn2c,go2c,gFCP(3),gFNP(3),gFPP(3)

  integer,save,allocatable,dimension(:) :: gpatch
  real(rkind),save,target,allocatable,dimension(:) :: gBA,gGP,gMT,gPR

  !-------------------------------------------------------------------------------
  !Clam model (CLAM) parameters and variables
  !-------------------------------------------------------------------------------
  real(rkind),save,target :: cpatch0
  real(rkind),save,target,allocatable,dimension(:) :: clam0,cfrmax,cTFR,csalt,cKDO,cDOh,cfTSSm,cRF, &
                                                    & cIFmax,cMTB,cTMT,cKTMT,cMRT,cn2c,cp2c
  real(rkind),save,target,allocatable,dimension(:,:) :: cKTFR,cKTSS,cTSS,calpha

  integer,save,allocatable,dimension(:) :: cpatch
  real(rkind),save,target,allocatable,dimension(:,:) :: CLAM,cFPOC,cFPON,cFPOP
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
  !pointer array defined for 1). spatially varying param,  2). output array
  !---------------------------------------------------------------------------
  type :: icm_pointer
    character(len=30) :: name  !parameter name
    integer :: ndim=0  !parameter dimension
    integer :: dims(3) !dimension info
    integer :: id=0    !id number
    integer :: itype=0 !used to identify data type
    real(rkind),dimension(50) :: data0 !oirginal value of data
    integer,allocatable,dimension(:,:) :: istat
    real(rkind),pointer :: p=>null()                   !param of single value
    real(rkind),pointer,dimension(:) :: p1=>null()     !param of 1D array
    real(rkind),pointer,dimension(:,:) :: p2=>null()   !param of 2D array
    real(rkind),pointer,dimension(:,:,:) :: p3=>null() !param of 2D array
    real(rkind),allocatable,dimension(:,:,:) :: data
  end type icm_pointer
  type(icm_pointer),save,target,allocatable,dimension(:) :: sp,wqout,wqhot

  !---------------------------------------------------------------------------
  !station outputs
  !---------------------------------------------------------------------------
  type :: station_var
    integer :: varid,vlen      !netcdf varid, var length
    character(len=30) :: varname=''  !variable name
    real(rkind),allocatable,dimension(:,:) :: data
  end type

  type :: station_data
    integer :: istat=0,it=0,nsta=0,nvar=0
    integer :: ncid,id_time
    real(rkind) :: time
    integer,allocatable :: ista(:)  !station index
    integer,allocatable :: iep(:)   !elem. in subdomain
    real(rkind),allocatable :: x(:),y(:),z(:)
    type(station_var) :: vars(200)
  end type
  type(station_data),save :: dg

end module icm_mod

module icm_interface
  interface brent_def
    subroutine brent(bv)
      use icm_mod, only: brent_var
      type(brent_var),target,intent(inout) :: bv
    end subroutine brent
  end interface
end module icm_interface
