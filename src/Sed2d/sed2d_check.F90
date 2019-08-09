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

SUBROUTINE sed2d_check
!--------------------------------------------------------------------
! This subroutine displays values computed by sed2d routines only for 
! the input parameters declared in this parameter list. It is only
! active if IPRE = 1 in sed2d.in.
! The main purpose of this routine is to check that there is
! no bug in the coded empirical formulae                            
!                                                                   
! Author: guillaume dodet (gdodet01@univ-lr.fr, gdodet@lnec.pt); 
! Date:   12/03/2013                                                 
!
! History:
!--------------------------------------------------------------------

  USE schism_glbl, ONLY : grav,rkind
  USE schism_msgp, ONLY : myrank,parallel_abort
  USE sed2d_mod, ONLY : idsed2d
  USE sed2d_friction
  USE sed2d_transport 

  IMPLICIT NONE

!- Local variable --------------------------------------------------- 
  REAL(rkind) :: beta,Cd,D_star,tau0,tau_cr,s,theta_cr,theta0,Ucr_c, &
                 Ucr_w,udir,z0cr,z0sw,z0wr
  REAL(rkind), DIMENSION(2) :: qb,qs,qtot,slope
  CHARACTER(LEN=*),PARAMETER :: FMT1='(A16,E11.4)'
!- User-defined parameters ------------------------------------------
  REAL(rkind), PARAMETER :: rho0 = 1027.d0  !Water density (kg/m3)
  REAL(rkind), PARAMETER :: rhos = 2650.d0  !Sediment density (kg/m3)
  REAL(rkind), PARAMETER :: nu   = 1.36e-6  !Water kine. visco.(m2/s)
  REAL(rkind), PARAMETER :: poro = 0.4d0    !Porosity (-)
  REAL(rkind), PARAMETER :: d50  = 0.0005d0 !d50 grain diameter (m) 
  REAL(rkind)            :: d90  = 0.000d0  !d90 grain diameter (m) 
  REAL(rkind), PARAMETER :: h    = 1.d0     !Water depth (m)
  REAL(rkind), PARAMETER :: u    = 2.d0     !Current velocity X (m/s)
  REAL(rkind), PARAMETER :: v    = 0.d0     !Current velocity Y (m/s)
  REAL(rkind)            :: z0   = 0.001d0  !Roughness length (m)  
  REAL(rkind), PARAMETER :: dpdx = 0.d0     !Bottom slope X (m/m)
  REAL(rkind), PARAMETER :: dpdy = 0.d0     !Bottom slope Y (m/m)
  REAL(rkind), PARAMETER :: urms = 1.d0     !Wave orbital vel. (m/s)
  REAL(rkind), PARAMETER :: hs   = 1.d0     !Sign. wave height (m)
  REAL(rkind), PARAMETER :: tp   = 10.d0    !Peak period (s)
  REAL(rkind), PARAMETER :: wl   = 50.d0     !Wave length (m)
!--------------------------------------------------------------------
! The following values can be specified in order to compare results
! with values from given in Soulsby, 1997:
!
! - Threshold of motion (p97)
! with d50 = 0.2 mm, d90 = 0.3 mm, h = 5 m, rho0 = 1027 kg/m*, 
! rhos = 2650 kg/m3, nu = 1.36e-6 m2/s:
!   * Ucr_c    = 0.39 m/s
!   * Ucr_w    = 0.18 m/s
!   * D*       = 4.06
!   * theta_cr = 0.0553
!   * tau_cr   = 0.176 N/m2
!
!  ...
!
!--------------------------------------------------------------------- 

  IF(myrank>0) CALL parallel_abort('sed2d_check only works on 1 proc')

!- Display user-defined parameters ----------------------------------
  WRITE(idsed2d,*)'----------------- CHECKS ------------------------'
  WRITE(idsed2d,*)'*** User-defined parameters ***'
  WRITE(idsed2d,*)'  rho0      = ',rho0
  WRITE(idsed2d,*)'  rhos      = ',rhos
  WRITE(idsed2d,*)'  nu        = ',nu
  WRITE(idsed2d,*)'  poro      = ',poro
  WRITE(idsed2d,*)'  d50       = ',d50
  WRITE(idsed2d,*)'  d90       = ',d90
  WRITE(idsed2d,*)'  depth     = ',h
  WRITE(idsed2d,*)'  u         = ',u
  WRITE(idsed2d,*)'  v         = ',v
  WRITE(idsed2d,*)'  dpdx      = ',dpdx
  WRITE(idsed2d,*)'  dpdy      = ',dpdy
  WRITE(idsed2d,*)'  urms      = ',urms
  WRITE(idsed2d,*)'  hs        = ',hs 
  WRITE(idsed2d,*)'  tp        = ',tp
  WRITE(idsed2d,*)'  wl        = ',wl
  WRITE(idsed2d,*)''

!- Display thresholds parameters ------------------------------------
  IF(d90<d50) d90 = 2.5*d50
  WRITE(idsed2d,*)'*** Threshold parameters ***'
  CALL compute_thresh(d50,d90,h,tp,D_star,Ucr_c,Ucr_w,theta_cr,tau_cr)
  WRITE(idsed2d,*)'  D*          = ',D_star
  WRITE(idsed2d,*)'  Ucr_c       = ',Ucr_c
  WRITE(idsed2d,*)'  Ucr_w       = ',Ucr_w
  WRITE(idsed2d,*)'  theta_cr    = ',theta_cr
  WRITE(idsed2d,*)'  tau_cr      = ',tau_cr
  WRITE(idsed2d,*)''

!- Display bedforms roughness ---------------------------------------
  WRITE(idsed2d,*)'*** Bedforms associated roughness ***'
  z0 = 0.d0
  z0cr = 0.d0
  z0wr = 0.d0
  z0sw = 0.d0
  CALL bedform_predictor_s97(d50,h,u,v,urms,tp,tau_cr,theta_cr,z0cr, &
                             z0wr,z0sw,z0)
  WRITE(idsed2d,*)'  z0 current-ripples = ',z0cr
  WRITE(idsed2d,*)'  z0 wave-ripples    = ',z0wr
  WRITE(idsed2d,*)'  z0 sandwave        = ',z0sw
  WRITE(idsed2d,*)'  z0 total           = ',z0
  WRITE(idsed2d,*)''

!- Display friction parameters --------------------------------------
  WRITE(idsed2d,*)'*** Friction parameters ***'
  WRITE(idsed2d,*)'  z0 skin     = ',z0
  WRITE(idsed2d,*)'  z0 ripple   = ',z0
  WRITE(idsed2d,*)'  z0 sandwave = ',z0

  CALL compute_drag(d50,h,z0,u,v,2,Cd)
  WRITE(idsed2d,*)'  Cd (z0s)    = ',Cd
  s = rhos/rho0
  tau0 = rho0*Cd*u*u
  theta0 = tau0/(grav*rho0*(s-1.d0)*d50)
  WRITE(idsed2d,*)'  tau0s       = ',tau0
  WRITE(idsed2d,*)'  theta0s       = ',theta0
  CALL compute_drag(d50,h,z0,u,v,2,Cd)
  WRITE(idsed2d,*)'  Cd (z0r)    = ',Cd
  CALL compute_drag(d50,h,z0,u,v,2,Cd)
  WRITE(idsed2d,*)'  Cd (z0sw)   = ',Cd
  CALL compute_drag(d50,h,z0,u,v,3,Cd)
  WRITE(idsed2d,*)'  Cd alluvial = ',Cd
  WRITE(idsed2d,*)''

!- Display transport rate -------------------------------------------
  udir = DATAN2(v,u)
  beta = -(dpdx*DCOS(udir)+dpdy*DSIN(udir))

  WRITE(idsed2d,*)'*** Transport formulae ***'
  CALL compute_drag(d50,h,z0,u,v,2,Cd)
  WRITE(idsed2d,*)'- Engelund and Hansen (1967)'
  CALL eh67(Cd,d50,u,v,beta,qtot)
  WRITE(idsed2d,*)'  Q total     = ',SQRT(qtot(1)**2+qtot(2)**2)
  WRITE(idsed2d,*)''

  WRITE(idsed2d,*)'- Ackers and White (1973)'
  CALL aw73(Cd,d50,u,v,h,beta,D_star,qtot)
  WRITE(idsed2d,*)'  Q total     = ',SQRT(qtot(1)**2+qtot(2)**2)
  WRITE(idsed2d,*)''

  WRITE(idsed2d,*)'- Soulsby and Van-Rijn (1997) '
  slope(1) = dpdx
  slope(2) = dpdy
  CALL svr97_bedl(Cd,d50,u,v,h,slope,Ucr_c,tau_cr,urms,qb)
  WRITE(idsed2d,*)'  Q bedload   = ',SQRT(qb(1)**2+qb(2)**2)
  CALL svr97_susp(Cd,d50,u,v,D_star,Ucr_c,urms,qs)
  WRITE(idsed2d,*)'  Q suspended = ',SQRT(qs(1)**2+qs(2)**2)
  qtot(:) = (qs(:)+qb(:))*(1.d0-1.6d0*beta)
  WRITE(idsed2d,*)'  Q total     = ',SQRT(qtot(1)**2+qtot(2)**2)
  WRITE(idsed2d,*)''

  WRITE(idsed2d,*)'- Van-Rijn (2007) '
  CALL vr07_bedl(Cd,d50,u,v,h,slope,urms,Ucr_c,Ucr_w,tau_cr,qb)
  WRITE(idsed2d,*)'  Q bedload   = ',SQRT(qb(1)**2+qb(2)**2)
  CALL vr07_susp(d50,u,v,urms,D_star,Ucr_c,Ucr_w,qs)
  WRITE(idsed2d,*)'  Q suspended = ',SQRT(qs(1)**2+qs(2)**2)
  qtot(:) = (qs(:)+qb(:))*(1.d0-1.6d0*beta)
  WRITE(idsed2d,*)'  Q total     = ',SQRT(qtot(1)**2+qtot(2)**2)
  WRITE(idsed2d,*)''
  WRITE(idsed2d,*)'-------------------------------------------------'
  WRITE(idsed2d,*)' '

END SUBROUTINE sed2d_check
