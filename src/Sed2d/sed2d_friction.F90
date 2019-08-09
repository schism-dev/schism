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

MODULE sed2d_friction
!--------------------------------------------------------------------
! This module contains several subroutines to compute friction -  
! related parameters:
! - compute_thresh
! - compute_drag
! - bedform_predictor_s97
! - bedform_predictor_v07
!                                                                  
! Author: guillaume dodet (gdodet01@univ-lr.fr, gdodet@lnec.pt)    
! Date:   20/02/2013
!
! History:                                               
! 02/2013 - G.Dodet: Removed compute_drag and passed irough as 
!                    argument in compute_rough
! 03/2013 - G.Dodet: - Added orbital threshold velocity in 
!                      compute_thresh;
!                    - Added alluvial friction method in compute_drag
! 06/2013 - G.Dodet: - Added two bedforms predictors;
!                    - Removed compute_rough (now included in 
!                      sed2d_main);
!--------------------------------------------------------------------

  CONTAINS

!--------------------------------------------------------------------
  SUBROUTINE compute_thresh(d50,d90,h,tp,D_star,Ucr_c,Ucr_w,theta_cr,&
                            tau_cr)
!--------------------------------------------------------------------
! This subroutine computes D* and threshold parameters
!
! Author: guillaume dodet (gdodet01@univ-lr.fr, gdodet@lnec.pt)
! Date:   20/02/2013
!
! History:
! 03/2013 - G.Dodet: Added orbital threshold velocity (Komar and
!                    Miller, 1974) 
!--------------------------------------------------------------------

  USE schism_glbl, ONLY : grav,rho0,rhosed,rkind
  USE schism_msgp, ONLY : parallel_abort
  USE sed2d_mod, ONLY : s,wvisco

  IMPLICIT NONE

!- Arguments --------------------------------------------------------
  REAL(rkind), INTENT(IN) :: d50,d90,h,tp
  REAL(rkind), INTENT(OUT) :: D_star,tau_cr,theta_cr,Ucr_c,Ucr_w 
!--------------------------------------------------------------------

!- Compute dimensionless grain size (Soulsby, 1997, p.41) -----------
  D_star = (grav*(s-1.d0)/wvisco**2.d0)**(1.d0/3.d0)*d50

!- Compute threshold Shields parameter ------------------------------
! Based on Soulsby and Whitehouse, 1997 (Soulsby, 1997, p.106)
!--------------------------------------------------------------------
  theta_cr = (0.3d0/(1.d0+1.2d0*D_star))+0.055d0*(1.d0-EXP(-0.02d0*D_star))

!- Compute threshold bed shear-stress (Soulsby, 1997, p.104) --------
  tau_cr = theta_cr*grav*d50*(rhosed-rho0)

!- Compute threshold velocities -------------------------------------
!- for currents (Soulsby, 1997, p.176)
!- for waves (Komar and Miller, 2004 in Soulsby, 1997, p.102)
!--------------------------------------------------------------------
!  IF((d50.GE.0.00005d0) .AND. (d50.LE.0.0005d0)) THEN
  IF(d50>=5d-5.AND.d50<=5d-4) THEN
    Ucr_c = 0.19d0*d50**0.1d0*LOG10(4.d0*h/d90)
    Ucr_w = (0.24d0*((s-1.d0)*grav)**0.66d0)*(d50**0.33d0)*(tp**0.33d0)
  ELSE IF(d50>5.d-4.AND.d50<=0.002d0) then
    Ucr_c = 8.5d0*d50**0.6d0*LOG10(4.d0*h/d90)
    Ucr_w = (0.95d0*((s-1.d0)*grav)**0.57d0)*(d50**0.43d0)*(tp**0.14d0)
  ELSE
    CALL parallel_abort('sed2d_friction: d50 not valid')
  ENDIF

  END SUBROUTINE compute_thresh

!--------------------------------------------------------------------
  SUBROUTINE compute_drag(d50,h,z0,u,v,idrag,Cd)
!--------------------------------------------------------------------
! This subroutine computes the drag coefficient
!
! Author: guillaume dodet (gdodet01@univ-lr.fr, gdodet@lnec.pt)
! Date:   12/03/2013
!
! History:
!--------------------------------------------------------------------
  USE schism_glbl, ONLY : grav,rho0,rkind
  USE schism_msgp, ONLY : parallel_abort
  USE sed2d_mod, ONLY : s

  IMPLICIT NONE

!- Arguments --------------------------------------------------------
  INTEGER, INTENT(IN) :: idrag
  REAL(rkind), INTENT(IN) :: d50,h,u,v,z0
  REAL(rkind), INTENT(OUT) :: Cd
!- Local variables --------------------------------------------------
  REAL(rkind) :: tau0,tau0s,theta,thetas,unorm,zrel
!- User-defined paramter --------------------------------------------
  REAL(rkind), PARAMETER :: zrel_max = 0.01d0 !0.05d0
!--------------------------------------------------------------------

  SELECT CASE(IABS(idrag))
    CASE(1) !Read Cd from main
    CASE(2) !Log-law formula
      zrel = MIN(z0/h,zrel_max)
      Cd = (0.4d0/(1+LOG(zrel)))**2.d0 !Cd<=0.0123
    CASE(3) !Alluvial friction method (Engelund, 1960, in S97, p.125)
      !Cd = (0.4d0/(1+LOG(z0/h)))**2.d0
      zrel = MIN(z0/h,zrel_max)
      Cd = (0.4d0/(1+LOG(zrel)))**2.d0 !Cd<=0.0123
      unorm = SQRT(u*u+v*v)
      tau0s = rho0*Cd*unorm*unorm
      if(s==1) call parallel_abort('SED_FRIC: (1)')
      thetas = tau0s/(grav*rho0*(s-1.d0)*d50)
      if(thetas<0.06) then
        Cd=0.0025
      else
        theta = 2.5d0*SQRT(thetas-0.06d0)
        tau0 = theta*grav*rho0*(s-1.d0)*d50
        IF(unorm == 0.d0) THEN
          CALL parallel_abort('sed2d_friction: alluvial method needs u')
        ELSE
          Cd = tau0/(rho0*unorm*unorm)
        ENDIF
      endif
  END SELECT
 
  END SUBROUTINE compute_drag

!--------------------------------------------------------------------
  SUBROUTINE bedform_predictor_vr07(d50,h,u,v,uorb,z0cr,z0mr,z0sw,z0)
!--------------------------------------------------------------------
! This subroutine computes the roughness associated to 
! current-induced bedforms and wave-induced bedforms based on 
! Van-Rijn(2007).

! Author: guillaume dodet (gdodet01@univ-lr.fr, gdodet@lnec.pt)
! Date:   17/05/2013
!
! History:
! 03/2014 - T.Guerin: Fixed error: z0 moved from local variable to 
!                     output argument
!--------------------------------------------------------------------
  USE schism_glbl, ONLY : grav,rkind
  USE schism_msgp, ONLY : parallel_abort
  USE sed2d_mod, ONLY : s

  IMPLICIT NONE

!- Arguments --------------------------------------------------------
  REAL(rkind), INTENT(IN) :: d50,h,u,v,uorb
  REAL(rkind), INTENT(OUT) :: z0,z0cr,z0mr,z0sw
!- Local variables --------------------------------------------------
  REAL(rkind) :: Uwc,phi,fcs,ffs,ks,kscr,ksmr,kss,kssw,unorm
!- User-defined parameters ------------------------------------------
  REAL(rkind), PARAMETER :: dgravel = 0.002, dsand = 0.000062
!--------------------------------------------------------------------
  unorm = SQRT(u*u+v*v)
  Uwc = SQRT(uorb*uorb+unorm*unorm)
  phi = Uwc*Uwc/((s-1)*grav*d50)

!- Skin roughness
  kss = 2.5d0*d50

!- Ripples
  IF(d50.LE.0.25d0*dgravel) THEN
    fcs = 1.d0
  ELSE
    fcs = (0.25d0*dgravel/d50)**1.5d0
  ENDIF
  kscr = MAX(fcs*d50*(85.d0-65.d0*TANH(0.015d0*(phi-150.d0))),0.d0)
  z0cr = kscr/30.d0

!- Mega-ripples
  IF(d50.GE.1.5d0*dsand) THEN
    ffs = 1.d0
  ELSE
    ffs = d50/(1.5*dsand)
  ENDIF
  ksmr = MAX(0.00002d0*ffs*h*(1.d0-DEXP(-0.05*phi))*(550-phi),0.d0)
  z0mr = ksmr/30.d0

!- Dunes
  kssw = MAX(0.00008d0*ffs*h*(1-DEXP(-0.002d0*phi))*(600-phi),0.d0) 
  z0sw = kssw/30.d0

  ks = kss+SQRT(kscr*kscr+ksmr*ksmr+kssw*kssw)
  z0 = ks/30.d0

  END SUBROUTINE bedform_predictor_vr07

!--------------------------------------------------------------------
  SUBROUTINE bedform_predictor_s97(d50,h,u,v,uorb,tp,tau_cr,theta_cr,&
                                  &z0cr,z0wr,z0sw,z0)
!--------------------------------------------------------------------
! This subroutine computes the roughness associated to 
! current-induced bedforms and wave-induced bedforms based on 
! Soulsby (1997).

! Author: guillaume dodet (gdodet01@univ-lr.fr, gdodet@lnec.pt)
! Date:   03/06/2013
!
! History:
!--------------------------------------------------------------------
  USE schism_glbl, ONLY : grav,pi,rho0,rhosed,rkind
  USE schism_msgp, ONLY : parallel_abort
  USE sed2d_mod, ONLY : s

  IMPLICIT NONE

!- Arguments --------------------------------------------------------
  REAL(rkind), INTENT(IN) :: d50,h,u,v,tau_cr,theta_cr,tp,uorb
  REAL(rkind), INTENT(INOUT) :: z0,z0cr,z0wr,z0sw
!- Local variables --------------------------------------------------
  REAL(rkind) :: A,fc,fwr,Hr_c,Hr_sw,Hr_w,Lr_c,Lr_sw,Lr_w,psi,r,     &
                 tau_c,tau_w,tau_wash,theta_w,Ts,unorm,z0s,z0cr0,    &
                 z0wr0,z0sw0
!- User-defined parameters ------------------------------------------
!--------------------------------------------------------------------
  unorm = SQRT(u*u+v*v)

!- Bedforms history
  z0cr0 = z0cr
  z0sw0 = z0sw
  z0wr0 = z0wr

!- Current-induced bed shear stress (skin friction only) 
  z0s = d50/12.d0
  if(LOG(z0s/h)==-1) call parallel_abort('SED_FRIC: (3)')
  fc = (0.4d0/(1.d0+LOG(z0s/h)))**2.d0
  tau_c = rho0*fc*unorm*unorm

  IF(tau_c > tau_cr) THEN !Initiation of motion
!--------------------------------------------------------------------
!- Dimensions of current ripples - Van Rijn, 1997 (eq. 81a-b)
!--------------------------------------------------------------------
!- Wash-out bed shear stress - Soulsby, 1997 (eq. 85b)
    tau_wash = 0.8d0*grav*rho0*(s-1.0d0)*d50
    
    IF(tau_c < tau_wash) THEN
      Lr_c = 1000.d0*d50
      Hr_c = Lr_c/7.d0
!      z0cr = MAX(1.d0*Hr_c*Hr_c/Lr_c,z0cr0)  
      z0cr = 1.d0*Hr_c*Hr_c/Lr_c
    ELSE !wash-out
      z0cr = 0.d0
    ENDIF

!--------------------------------------------------------------------
!- Dimensions of current sand-waves - Van Rijn, 1984 (eq. 83b)
!--------------------------------------------------------------------
!- Wash-out condition based on shear stress - Soulsby, 1997 (eq. 83b)
    tau_wash = 26.d0*tau_cr
    IF(tau_c < tau_wash) THEN
      if(h==0.or.tau_cr==0) call parallel_abort('SED_FRIC: (4)')
      Ts = (tau_c-tau_cr)/tau_cr
      Hr_sw = MAX(0.d0,0.11d0*h*(d50/h)**0.3d0*(1d0-DEXP(-0.5d0*Ts))*&
              (25.d0-Ts))
      Lr_sw = 7.3d0*h
!      z0sw = MAX(1.d0*Hr_sw*Hr_sw/Lr_sw,z0sw0) 
      z0sw = 1.d0*Hr_sw*Hr_sw/Lr_sw
    ELSE !wash-out
      z0sw = 0.d0
    ENDIF
  ELSE !No initiation of motion
    z0cr = z0cr0
    z0sw = z0sw0
  ENDIF

!-------------------------------------------------------------------- 
!- Dimensions of wave ripples - Nielsen, 1992 (eq. 89a-d)
!--------------------------------------------------------------------
!- Wave-induced friction coefficient - Swart, 1974 (eq. 60a-b)
  A = uorb*tp/(2.d0*pi)
  r = A/(2.5d0*d50)
  IF(r.GT.0.d0) THEN
    IF(r <= 1.57D0) THEN
      fwr = 0.3D0
    ELSE
      fwr = 0.00251D0*DEXP(5.21D0*r**(-0.19D0))
    ENDIF 
    tau_w = 0.5D0*rho0*fwr*uorb*uorb
  ELSE
    tau_w = 0.d0
  ENDIF
  theta_w = tau_w/(grav*(rhosed-rho0)*d50)

  IF(theta_w .GT. theta_cr) THEN !Initiation of motion
!- Wash-out condition based on  mobility parameter
    psi = uorb*uorb/((s-1)*grav*d50) 
    IF(psi < 156.d0 .AND. theta_w < 0.831) THEN
      Hr_w = (0.275d0-0.022d0*psi**0.5d0)*A
      Lr_w = Hr_w/(0.182d0-0.24*theta_w**1.5d0)
!      z0wr = MAX(0.267d0*Hr_w*Hr_w/Lr_w,z0wr0)
      z0wr = 0.267d0*Hr_w*Hr_w/Lr_w
    ELSE
      Lr_w = 1d0
      Hr_w = 0d0
      z0wr = 0.d0
    ENDIF
  ELSE !No initiation of motion
    z0wr = z0wr0
  ENDIF
  z0 = MAX((MAX(z0cr,z0wr)+z0sw),z0s)

  END SUBROUTINE bedform_predictor_s97
!--------------------------------------------------------------------
END MODULE sed2d_friction
