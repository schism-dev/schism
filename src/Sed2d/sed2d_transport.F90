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

MODULE sed2d_transport
!--------------------------------------------------------------------
! This module contains several subroutines/formulae to compute the 
! sediment transport rate at one location:
!
! - eh67       Total transport from Engelund and Hansen (1967)        
! - aw73       Total transport from Ackers and White (1973)
! - svr97_bedl Bedload transport from Soulsby and Van-Rijn (1997)
! - svr97_susp Suspended transport from Soulsby and Van-Rijn (1997)
! - vr07_bedl  Bedload transport from Soulsby (2007)
! - vr07_susp  Suspended transport from Soulsby (2007)
! - cl11_tot   Total transport (bedload + suspended) from Camenen
!              and Larson (2011)
! - wl14_tot   Total transport (bedload + suspended) from Wu and Lin (2014)
! - tr04_tot   Total transport (bedload + suspended from TRANSPOR2004 formula
!              of van Rijn et al. (2004) (final report Z3748, Delft Hydraulics)
!
! and
!
! - lesser04_slope Coefficients for slope effect on bed load from 
!                  Lesser et al. (2004)
! - wave_asymmetry_Elfrink Compute wave orbital bottom velocities
!                          for asymmetric waves
!                                                                
! Author: guillaume dodet (gdodet01@univ-lr.fr, gdodet@lnec.pt)    
! Date:   20/02/2013
!
! History:                                                
! 03/2013 - G.Dodet: - Decomposed svr97 into a bedload and a 
!                      suspended transport routine;
!                    - Added Van-Rijn (2007);
!                    - Put computation of coefficients for slope
!                      effects on bedload (Lesser et al. 2004) in a 
!                      distinct subroutine;
! 07/2013 - T.Guerin: - Added Camenen and Larson (2011) and wave 
!                       asymmetry computation based on Elfrink (2006)
! 05/2017 - T.Guerin: - Added Wu and Lin (2014) formula
!                     - Added TR2004 (van Rijn et al., 2004) formula
!--------------------------------------------------------------------

  IMPLICIT NONE

CONTAINS

!--------------------------------------------------------------------
  SUBROUTINE eh67(Cd,d50,u,v,beta,qtot)
!--------------------------------------------------------------------
! This subroutine computes the total sediment transport rate based on 
! Engelund-Hansen (1967) see S97, p.175 
!
! NB: Cd should be determined by alluvial friction method (Engelund,
!     1966, see S97,p175)
!
! Author: guillaume dodet (gdodet01@univ-lr.fr, gdodet@lnec.pt)
! Date: 20/02/2013
!
! History:
! 03/2013 - G.Dodet: Added slope effect on total transport
!--------------------------------------------------------------------

  USE schism_glbl, ONLY : grav,rkind
  USE schism_msgp, ONLY: parallel_abort
  USE sed2d_mod, ONLY : s

  IMPLICIT NONE

!- Arguments --------------------------------------------------------
  REAL(rkind), INTENT(IN) :: beta,Cd,d50,u,v
  REAL(rkind), DIMENSION(2), INTENT(OUT) :: qtot
!- Local variables --------------------------------------------------
  REAL(rkind) :: C,unorm
!--------------------------------------------------------------------

  unorm = DSQRT(u*u+v*v)
  if(Cd<=0) call parallel_abort('eh67: Cd<=0')
  C = DSQRT(grav/Cd)

  qtot(1) = (1.d0-1.6d0*beta)*(0.05d0/((s-1.d0)**2.d0*d50*C**3.d0*   &
            DSQRT(grav)))*(unorm**4.d0)*u
  qtot(2) = (1.d0-1.6d0*beta)*(0.05d0/((s-1.d0)**2.d0*d50*C**3.d0*   &
            DSQRT(grav)))*(unorm**4.d0)*v

  END SUBROUTINE eh67

!--------------------------------------------------------------------
  SUBROUTINE aw73(Cd,d50,u,v,h,beta,D_star,qtot)
!--------------------------------------------------------------------
! This subroutine computes the total sediment transport rate based on
! Ackers and White (1973), S97, p.175
!
! NB: - with this formula, d50 in sed2d.in should be replaced by
!       the equivalent d35
!     - Cd should be determined by White et al.(1980) friction
!       alluvial method. see S97, p.175
!
! Author: guillaume dodet (gdodet01@univ-lr.fr, gdodet@lnec.pt)
! Date: 20/02/2013
!
! History:
! 03/2013 - G.Dodet: Added slope effect on total transport
!--------------------------------------------------------------------

  USE schism_glbl, ONLY : grav,rkind
  USE schism_msgp, ONLY : parallel_abort
  USE sed2d_mod, ONLY : s

  IMPLICIT NONE

!- Arguments --------------------------------------------------------
  REAL(rkind), INTENT(IN) :: beta,Cd,D_star,d50,h,u,v
  REAL(rkind), DIMENSION(2), INTENT(OUT) :: qtot
!- Local variables --------------------------------------------------
  REAL(rkind) :: Aaw,Caw,Faw,lslope,m,n,udir,unorm,u_star
!--------------------------------------------------------------------

  unorm = DSQRT(u*u+v*v)
  udir = DATAN2(v,u)
  qtot = 0.d0

  IF((D_star.GE.1.d0).AND.(D_star.LE.60.d0)) THEN
    n = 1.d0-0.243d0*DLOG(D_star)
    Aaw = 0.23d0/DSQRT(D_star)+0.14d0
    m = 6.83d0/D_star+1.67d0
    Caw = DEXP(2.79d0*DLOG(D_star)-0.426d0*(DLOG(D_star))**2.d0-7.97d0)
  ELSEIF(D_star.GT.60) THEN
    n = 0.d0
    Aaw = 0.17d0
    m = 1.78d0
    Caw = 0.025d0
  ELSE
    CALL parallel_abort('sed2d: D_star out of range in AW73')
  ENDIF

  IF(unorm.NE.0.d0) THEN
    u_star = DSQRT(Cd)*unorm
    Faw = (u_star**n/DSQRT(grav*(s-1.d0)*d50))*(unorm/(2.46d0*DLOG   &
          (10.d0*h/d50)))**(1.d0-n)
    IF(Faw.GE.Aaw) THEN
      qtot(1) = (1.d0-1.6d0*beta)*((Caw*unorm*d50*(unorm/u_star)**n)*&
                ((Faw-Aaw)/Aaw)**m)*DCOS(udir)
      qtot(2) = (1.d0-1.6d0*beta)*((Caw*unorm*d50*(unorm/u_star)**n)*&
                ((Faw-Aaw)/Aaw)**m)*DSIN(udir)
    ENDIF
  ENDIF

  END SUBROUTINE aw73

!--------------------------------------------------------------------
  SUBROUTINE svr97_bedl(Cd,d50,u,v,h,dpdxy,U_cr,tau_cr,urms,qb)
!--------------------------------------------------------------------
! This subroutine computes the bed-load transport rate based on 
! Soulsby - Van Rijn (S97,p. 183)
!
! NB: This formula was computed for rippled bed. A value of z0 = 6 mm
!     should be used (S97, p.184)
! Cd = (0.4d0/(1.d0+DLOG(0.006d0/htot)))**2.d0
!
! Author: guillaume dodet (gdodet01@univ-lr.fr, gdodet@lnec.pt)
! Date: 11/03/2013
!
! History:
!--------------------------------------------------------------------

  USE schism_glbl, ONLY : grav,rkind
  USE schism_msgp, ONLY : parallel_abort
  USE sed2d_mod, ONLY : s

  IMPLICIT NONE

!- Arguments --------------------------------------------------------
  REAL(rkind), INTENT(IN) :: Cd,d50,h,tau_cr,u,U_cr,urms,v
  REAL(rkind), DIMENSION(2), INTENT(IN) :: dpdxy
  REAL(rkind), DIMENSION(2), INTENT(OUT) :: qb
!- Local variables --------------------------------------------------
  REAL(rkind) :: alpha_n,alpha_s,Asb,aux,cff,udir,unorm
!--------------------------------------------------------------------

!- Slope effect on bedload transport as in Lesser et al (2004)
  CALL lesser04_slope(Cd,d50,u,v,dpdxy,tau_cr,alpha_s,alpha_n)

!- Transport rate
  if(Cd<=0) call parallel_abort('svr97_bedl: Cd<=0')
  qb = 0.d0
  unorm = SQRT(u*u+v*v)
  udir = ATAN2(v,u)
  Asb = (0.005d0*h*(d50/h)**1.2d0)/(((s-1.d0)*grav*d50)**1.2d0)
  aux = SQRT(unorm**2.d0+(0.018d0/Cd)*urms**2.d0)
  IF(((aux-U_cr).GE.0.d0)) THEN
    cff = Asb*unorm*(aux-U_cr)**2.4d0
    qb(1) = alpha_s*cff*COS(udir)-alpha_n*alpha_s*cff*SIN(udir)
    qb(2) = alpha_s*cff*SIN(udir)-alpha_n*alpha_s*cff*COS(udir)
  ENDIF

  END SUBROUTINE svr97_bedl

!--------------------------------------------------------------------
  SUBROUTINE svr97_susp(Cd,d50,u,v,D_star,U_cr,urms,qs)
!--------------------------------------------------------------------
! This subroutine computes the suspended and total sediment transport
! rate based on Soulsby - Van Rijn (S97,p. 183)
!
! NB: This formula was computed for rippled bed. A value of z0 = 6 mm
!     should be used (S97, p.184)
! Cd = (0.4d0/(1.d0+DLOG(0.006d0/htot)))**2.d0
!
! Author: guillaume dodet (gdodet01@univ-lr.fr, gdodet@lnec.pt)
! Date: 11/03/2013
!
! History:
!--------------------------------------------------------------------

  USE schism_glbl, ONLY : grav,rkind
  USE schism_msgp, ONLY : parallel_abort
  USE sed2d_mod, ONLY : s

  IMPLICIT NONE

!- Arguments --------------------------------------------------------
  REAL(rkind), INTENT(IN) :: Cd,d50,u,v,D_star,U_cr,urms
  REAL(rkind), DIMENSION(2), INTENT(OUT) :: qs
!- Local variables --------------------------------------------------
  REAL(rkind) :: Ass,aux,udir,unorm
!--------------------------------------------------------------------

  if(Cd<=0.or.D_star==0) call parallel_abort('svr97_susp: ')
  qs = 0.d0
  unorm = SQRT(u*u+v*v)
  udir = ATAN2(v,u)
  Ass = (0.012d0*d50*D_star**(-0.6d0))/(((s-1.d0)*grav*d50)**1.2d0)
  aux = SQRT(unorm**2.d0+(0.018d0/Cd)*urms**2.d0)

!- Transport rate
  IF(((aux-U_cr).GE.0.d0)) THEN
    qs(1) = (Ass*unorm*(aux-U_cr)**2.4d0)*COS(udir)
    qs(2) = (Ass*unorm*(aux-U_cr)**2.4d0)*SIN(udir)
  ENDIF

  END SUBROUTINE svr97_susp

!--------------------------------------------------------------------
  SUBROUTINE vr07_bedl(Cd,d50,u,v,h,dpdxy,uorb,Ucr_c,Ucr_w,tau_cr,qb)
!--------------------------------------------------------------------
! This subroutine computes the bed-load transport rate based on
! Van Rijn (2007)
!
! Author: guillaume dodet (gdodet01@univ-lr.fr, gdodet@lnec.pt)
! Date: 11/03/2013
!
! History:
! 04/2013 - G.Dodet: Added a maximum velocity value based on Van Rijn
!                    (2007) to compute qb.
! 09/2014 - T.Guerin: Added condition for unorm+uorb=0 case
!--------------------------------------------------------------------

  USE schism_glbl, ONLY : grav,rkind
  USE schism_msgp, ONLY : parallel_abort
  USE sed2d_mod, ONLY : s  

  IMPLICIT NONE

!- Arguments --------------------------------------------------------
  REAL(rkind), INTENT(IN) :: Cd,d50,h,tau_cr,u,Ucr_c,Ucr_w,uorb,v
  REAL(rkind), DIMENSION(2), INTENT(IN) :: dpdxy
  REAL(rkind), DIMENSION(2), INTENT(OUT) :: qb
!- Local variables --------------------------------------------------
  REAL(rkind) :: alpha_n,alpha_s,beta,cff,k,Me,Ucr,ueff,udir,unorm  
!- User-defined parameters ------------------------------------------
  REAL(rkind), PARAMETER :: gama = 0.4d0       !0.8 for regular waves 
!--------------------------------------------------------------------

!- Slope effect on bedload transport as in Lesser et al (2004)
  CALL lesser04_slope(Cd,d50,u,v,dpdxy,tau_cr,alpha_s,alpha_n)

!- Transport rate
  qb = 0.d0
  unorm = MIN(SQRT(u*u+v*v),1.8d0) ! Range of validity for the formula
  udir = ATAN2(v,u)
  ueff = unorm+gama*uorb

  IF (unorm+uorb==0.d0) THEN
    qb(:) = 0.d0
  ELSE
    beta = unorm/(unorm+uorb)
    Ucr = beta*Ucr_c+(1.d0-beta)*Ucr_w
    Me = (ueff-Ucr)/SQRT((s-1.d0)*grav*d50)

    IF(ueff-Ucr.GE.0.d0) THEN
      cff = 0.015d0*unorm*h*((d50/h)**1.2d0)*(Me**1.5d0)
      qb(1) = alpha_s*cff*COS(udir)-alpha_n*alpha_s*cff*SIN(udir)
      qb(2) = alpha_s*cff*SIN(udir)-alpha_n*alpha_s*cff*COS(udir)
    ENDIF
  ENDIF

  END SUBROUTINE vr07_bedl

!--------------------------------------------------------------------
  SUBROUTINE vr07_susp(d50,u,v,uorb,D_star,Ucr_c,Ucr_w,qs)
!--------------------------------------------------------------------
! This subroutine computes the suspended transport rate based on
! Van Rijn (2007)
!
! Author: guillaume dodet (gdodet01@univ-lr.fr, gdodet@lnec.pt)
! Date: 11/03/2013
!
! History:
! 04/2013 - G.Dodet: Added a maximum velocity value based on Van Rijn
!                    (2007) to compute qs.
! 09/2014 - T.Guerin: Added condition for unorm+uorb=0 case
!--------------------------------------------------------------------

  USE schism_glbl, ONLY : grav,rkind
  USE schism_msgp, ONLY : parallel_abort
  USE sed2d_mod, ONLY : s

  IMPLICIT NONE

!- Arguments --------------------------------------------------------
  REAL(rkind), INTENT(IN) :: D_star,d50,u,Ucr_c,Ucr_w,uorb,v
  REAL(rkind), DIMENSION(2), INTENT(OUT) :: qs
!- Local variables --------------------------------------------------
  REAL(rkind) :: beta,cff,k,Me,Ucr,udir,ueff,unorm
!- User-defined parameters ------------------------------------------
  REAL(rkind), PARAMETER :: gama = 0.4d0       !0.8 for regular waves
!--------------------------------------------------------------------

!- Transport rate
  qs = 0.d0
  unorm = MIN(SQRT(u*u+v*v),1.8d0) !Range of validity for the formula
  udir = DATAN2(v,u)
  ueff = unorm+gama*uorb

  IF (unorm+uorb==0.d0) THEN
    qs(:) = 0.d0
  ELSE
    beta = unorm/(unorm+uorb)
    Ucr = beta*Ucr_c+(1.d0-beta)*Ucr_w
    Me = (ueff-Ucr)/SQRT((s-1.d0)*grav*d50)

    IF(ueff-Ucr.GE.0.d0) THEN
      if(D_star<=0) call parallel_abort('vr07_susp: D_star<=0')
      cff = 0.012d0*unorm*d50*Me**2.4d0*D_star**(-0.6d0)
      qs(1) = cff*COS(udir)
      qs(2) = cff*SIN(udir)
    ENDIF
  ENDIF

  END SUBROUTINE vr07_susp

!--------------------------------------------------------------------
  SUBROUTINE cl11_tot(d50,D_star,dpdxy,hs,htot,theta_cr,tp,u,v,wdir, &
                     &wlpeak,qb,qs)
!--------------------------------------------------------------------
! This subroutine computes total sediment transport rate (m^2/s) 
! (bedload + suspended load) from Camenen and Larson (2011)
!
! Author: thomas guerin (thomas.guerin@univ-lr.fr)    
! Date: 24/04/2013
!
! History:
!--------------------------------------------------------------------

  USE schism_glbl, ONLY: grav,pi,rho0,rkind
  USE schism_msgp, ONLY : parallel_abort
  USE sed2d_mod, ONLY: iasym,s,wvisco

#ifdef USE_WWM
  USE DATAPOOL, ONLY: gammab =>BRHD
#endif

  IMPLICIT NONE

!- Arguments --------------------------------------------------------  
  REAL(rkind), INTENT(IN) :: d50,D_star,hs,htot,theta_cr,tp,u,v,wdir,&
                             wlpeak
  REAL(rkind), DIMENSION(2), INTENT(IN) :: dpdxy
  REAL(rkind), DIMENSION(2), INTENT(OUT) :: qb,qs
!- User-defined parameters ------------------------------------------  
  REAL(rkind), PARAMETER :: ksd_coef = 2.5d0 ! = 2 or 2.5
!- Constants --------------------------------------------------------  
  REAL(rkind), PARAMETER :: kappa = 0.4d0
  REAL(rkind), PARAMETER :: an = 12.d0
  REAL(rkind), PARAMETER :: b = 4.5d0
  REAL(rkind), PARAMETER :: pwr_cr = 100.d0
  REAL(rkind), PARAMETER :: Ad = 0.7d0
  REAL(rkind), PARAMETER :: d_ks_c = 0.0001d0
  INTEGER, PARAMETER :: ech = 50 !sampling of wave orbital velocity uw(t)
!- Local variables --------------------------------------------------
  REAL(rkind) :: a_b,A_epsi,A_w,AcR,alpha_off,alpha_on,alpha_plb,    &
                 alpha_pls,aux,aw,cR,D_b,D_c,D_w,delta_w,dpdxy_tot,  &
                 dpdxy_wdir,Dpl_off,Dpl_on,epsi,fc,fcw,fw,Hr,Hr_c,   &
                 Hr_w,Hrms,irribaren,k_b,k_c,k_cw_c,k_cw_w,          &
                 kpeak,ksd,ks_c,ks_w,ksf_c,ksf_w,kss_w,Lr_c,Lr_w,    &
                 omega,omega_off,omega_on,phi,psi_w,pwr,qb_n,qb_w,   &
                 qs_n,qs_w,R,rw,sigma_c,sigma_cw,sigma_w,slope_dir,  &
                 slope_wdir_angle,T_crest,T_trough,tau_c,tau_w,      &
                 theta_c,theta_cn,theta_cw,theta_cw_on,theta_cw_off, &
                 theta_cwm,theta_cx,theta_cy,theta_net,theta_w,      &
                 theta_wm,U_crest,U_trough,Uc,Uc_n,Uc_net,Uc_w,      &
                 Uc_x,Uc_y,Ucw_off,Ucw_on,Ucw2_off,Ucw2_on,udir,Uorb,&
                 ustar_c,ustar_w,Uw_crsf,v_off,v_on,wl0,Ws,Xt,Xv,Y, &
                 tmp,tmp2
  REAL(rkind), DIMENSION(ech-1) :: Uorb_Elfrink
  INTEGER :: i
!--------------------------------------------------------------------

  !- constant
  aux = 2.d0*(s-1.d0)*grav*d50 !>0

  !- settling velocity (Soulsby 1997)
  Ws = wvisco*((10.36d0**2+1.049d0*D_star**3)**0.5d0 - 10.36d0) /d50 !>0

  !- current amplitude and direction
  Uc = DSQRT(u*u+v*v)
  udir = DATAN2(v,u)

  IF (hs>0.AND.tp>0.AND.wlpeak>0) THEN
    !- wave angular frequency
    omega = 2.d0*pi/tp !>0

    !- peak wave number
    kpeak = 2.d0*pi/wlpeak !>0

    !- angle between current and waves directions
    phi = udir - wdir

    !- current velocity components in the wave frame
    Uc_w = Uc*DCOS(phi)
    Uc_n = Uc*DSIN(phi)
  ELSE
    Uc_x = Uc*DCOS(udir)
    Uc_y = Uc*DSIN(udir)
  ENDIF

  IF (hs>0.AND.tp>0.AND.wlpeak>0) THEN
    IF (iasym==0) THEN
      !- Airy waves
      tmp=DSINH(kpeak*htot)
      if(tmp==0) call parallel_abort('cl11_tot: tmp=0')
      Uorb = pi*hs/tp/tmp !>0
      U_crest = Uorb !>0
      T_crest = tp/2.d0 !>0
      T_trough = T_crest !>0

    ELSEIF (iasym==1) THEN
      !- Wave asymmetry treatment : compute wave orbital velocity
      !  (horizontal component) from Elfrink (2006)
      CALL wave_asymmetry_Elfrink(hs,tp,wdir,htot,-dpdxy(1),-dpdxy(2),&
           ech,Uorb_Elfrink,U_crest,U_trough,T_crest,T_trough)

      IF (U_crest<=0.OR.U_trough<=0.OR.T_crest<=0.OR.T_trough<=0) THEN
        !- Airy waves
        tmp=DSINH(kpeak*htot)
        if(tmp==0) call parallel_abort('cl11_tot: tmp=0')
        Uorb = pi*hs/tp/tmp !>0
        U_crest = Uorb !>0
        T_crest = tp/2.d0 !>0
        T_trough = T_crest !>0
      ELSE
        Uorb = (U_crest+U_trough)/2.d0 !>0
      ENDIF

    ENDIF !iasym
  ENDIF

  !- total bed roughness and friction factors -----------------------

  !- grain-related roughness
  ksd = ksd_coef*d50

  !--- total bed roughness for current alone
  Lr_c = 1000.d0*d50              !ripple length
  Hr_c = Lr_c/7.d0                !ripple height
  ksf_c = 7.5d0*Hr_c*Hr_c/Lr_c   !form-drag roughness

  !- iterative method
  if(htot<=0) call parallel_abort('cl11_tot: htot<=0')
  ks_c = d_ks_c
  DO WHILE ((-ks_c+ksd+ksf_c + 5.d0*kappa*kappa*Uc*Uc/               &
            &(grav*(s-1.d0)*(1.d0+DLOG(ks_c/30.d0/htot))**2.d0))      &
            &.GT. 0.d0)
    ks_c = ks_c + d_ks_c !total bed roughness for current alone       
    IF (ks_c .GE. (0.1d0*htot)) THEN
      ks_c = 0.1d0*htot !maximum value of ks_c is set to htot/10
      EXIT
    ENDIF
  ENDDO

  !- friction factor for current alone
  fc = 2.d0*(kappa/(1.d0+DLOG(ks_c/30.d0/htot)))**2.d0 !>0

  IF (hs>0.AND.tp>0.AND.wlpeak>0) THEN
    !--- total bed roughness for waves alone
    Aw = Uorb*tp/2.d0/pi
    psi_w = 2.d0*Uorb*Uorb/aux
    IF (psi_w .LT. 10.d0) THEN
      Hr_w = 0.22d0*Aw
      Lr_w = 1.25d0*Aw
    ELSEIF ((psi_w .GE. 10.d0) .AND. (psi_w .LT. 250.d0)) THEN
      Hr_w = 2.8d0*1.0e-13*(250d0-psi_w)**5d0*Aw
      Lr_w = 1.4d0*1.0e-6*(250d0-psi_w)**2.5d0*Aw
    ELSEIF (psi_w .GE. 250.d0) THEN
      Hr_w = 0.d0
      Lr_w = 0.d0
    ENDIF
    ksf_w = 7.5d0*Hr_w*Hr_w/(Lr_w + 1.0e-6)   !form-drag roughness

    !- calculation of friction factor for waves (Swart 1974)
    R = Aw/ksd
    fw = DEXP(5.21d0*R**(-0.19d0) - 6.d0)
    kss_w = 2.5d0*fw*Uorb*Uorb/grav   !sediment-related roughness
    ks_w = ksd + ksf_w + kss_w
    IF ((Aw/ks_w) .LE. 1.57d0) THEN
      fw = 0.3d0
    ENDIF

    !--- friction factor for waves and current
    Xv = Uc/(Uc+Uorb+1.0e-6)
    fcw = Xv*fc + (1.d0-Xv)*fw

    !- bed shear-stress for waves only
    tau_w = 0.5d0*rho0*fw*Uorb*Uorb !>=0
  ENDIF

  !- bed shear-stress for current only
  tau_c = 0.5d0*rho0*fc*Uc*Uc

  !--- Shields number
  !- current Shield nb
  theta_c = tau_c / (rho0*0.5d0*aux) !>=0

  IF (hs>0.AND.tp>0.AND.wlpeak>0) THEN
    !- current Shield nb relative to normal current component in the wave frame
    theta_cn = fc*Uc_n*Uc_n / aux !>=0

    !- maximum and mean wave Shield nb
    theta_w = tau_w / (rho0*0.5d0*aux) !>=0
    theta_wm = theta_w/2.d0   !valid for a sinusoidal wave profile (!)

    !- maximum and mean wave Shield nb due to wave-current interaction
    tmp=theta_c*theta_c + theta_w*theta_w+ 2.d0*theta_w*theta_c*DCOS(phi)
    tmp2=theta_c*theta_c + theta_wm*theta_wm+ 2.d0*theta_wm*theta_c*DCOS(phi)
    if(tmp<=0.or.tmp2<0) call parallel_abort('cl11_tot: tmp<=0.or.tmp2<0')

    theta_cw=sqrt(tmp)
    theta_cwm=sqrt(tmp2)
!  theta_cw = DSQRT(theta_c*theta_c + theta_w*theta_w                 &
!                   + 2.d0*theta_w*theta_c*DCOS(phi))
!  theta_cwm = DSQRT(theta_c*theta_c + theta_wm*theta_wm              &
!                    + 2.d0*theta_wm*theta_c*DCOS(phi))
  ELSE
    theta_cx = fc*Uc_x*Uc_x / aux !>=0
    theta_cy = fc*Uc_y*Uc_y / aux !>=0
  ENDIF

  !--- Bed load transport -------------------------------------------

  IF (hs>0.AND.tp>0.AND.wlpeak>0) THEN
    !- mean values of the instantaneous Shield nb over the wave half periods
    if(omega==0.or.T_trough<=0.or.T_crest<=0) then
       write(12,*) 'omega=',omega
       write(12,*) 'T_trough=',T_trough
       write(12,*) 'T_crest=',T_crest
       call parallel_abort('cl11_tot: omega==0')
    endif
    Ucw2_on = (Uc_w*Uc_w*T_crest                                       &
               + 2.d0/omega*Uc_w*Uorb*(1.d0-DCOS(omega*T_crest))       &
               + Uorb*Uorb/2.d0*(T_crest-DSIN(2.d0*omega*T_crest)      &
               /2.d0/omega)) /T_crest
    Ucw2_off = (Uc_w*Uc_w*T_trough + 2.d0/omega*Uc_w*Uorb              &
                *(DCOS(omega*T_crest)-DCOS(omega*tp))                  &
                + Uorb*Uorb/2.d0*(T_trough+(DSIN(2.d0*omega*T_crest)   &
                -DSIN(2.d0*omega*tp))/2.d0/omega)) /T_trough
    theta_cw_on = fcw*Ucw2_on / aux
    theta_cw_off = -fcw*Ucw2_off / aux

    !- calculation of the phase-lag effects coefficient (Camenen and Larson 2006)
    delta_w = DSQRT(wvisco*tp/pi)
    rw = U_crest/Uorb - 1.d0   !wave asymmetry parameter
    Uw_crsf = 8.35d0*DSQRT((s-1.d0)*grav*DSQRT(d50*delta_w))*(1.d0+rw)
    if(Ucw2_on<=0.or.Ucw2_off<=0) call parallel_abort('cl11_tot: Ucw2_on<0')
    Ucw_on = DSQRT(Ucw2_on)
    Ucw_off = DSQRT(Ucw2_off)
    alpha_on = wvisco**0.25d0*Ucw_on**0.5d0                            &
               *DEXP(-(Uw_crsf/Ucw_on)**2.d0) /(Ws*T_crest**0.75d0)
    alpha_off = wvisco**0.25d0*Ucw_off**0.5d0                          &
                *DEXP(-(Uw_crsf/Ucw_off)**2.d0) /(Ws*T_trough**0.75d0)
    alpha_plb = alpha_on-alpha_off

    !- net sediment transporting Shields nb
    theta_net = (1.d0-alpha_plb)*theta_cw_on                           &
                + (1.d0+alpha_plb)*theta_cw_off

    !--- bed load transport components 
    !- parallel and normal to the wave direction

    if(theta_c+theta_w==0) call parallel_abort('cl11_tot: theta_c+theta_w==0')
    Xt = theta_c/(theta_c+theta_w)
    aw = 6.d0+6.d0*Xt

    IF (theta_net .GE. 0.d0) THEN
        qb_w = DSQRT((s-1.d0)*grav*d50*d50*d50)*aw*DSQRT(theta_net)    &
               *theta_cwm*DEXP(-b*theta_cr/theta_cw)
    ELSE
        qb_w = -DSQRT((s-1.d0)*grav*d50*d50*d50)*aw*DSQRT(-theta_net)  &
               *theta_cwm*DEXP(-b*theta_cr/theta_cw)
    ENDIF
    IF (Uc_n .GE. 0.d0) THEN
        qb_n = DSQRT((s-1.d0)*grav*d50*d50*d50)*an*DSQRT(theta_cn)     &
               *theta_cwm*DEXP(-b*theta_cr/theta_cw)
    ELSE
        qb_n = -DSQRT((s-1.d0)*grav*d50*d50*d50)*an*DSQRT(theta_cn)    &
               *theta_cwm*DEXP(-b*theta_cr/theta_cw)
    ENDIF

    !- x,y components
    qb(1) = qb_w*DCOS(wdir) - qb_n*DSIN(wdir)
    qb(2) = qb_w*DSIN(wdir) + qb_n*DCOS(wdir)

  ELSE
    aw = 12.d0

    IF (theta_c .GT. 0.d0) THEN
      IF (Uc_x .GE. 0.d0) THEN
        qb(1) = DSQRT((s-1.d0)*grav*d50*d50*d50)*aw*DSQRT(theta_cx)    &
                 *theta_c*DEXP(-b*theta_cr/theta_c)
      ELSE
        qb(1) = -DSQRT((s-1.d0)*grav*d50*d50*d50)*aw*DSQRT(theta_cx)   &
                 *theta_c*DEXP(-b*theta_cr/theta_c)
      ENDIF
      IF (Uc_y .GE. 0.d0) THEN
        qb(2) = DSQRT((s-1.d0)*grav*d50*d50*d50)*an*DSQRT(theta_cy)    &
                 *theta_c*DEXP(-b*theta_cr/theta_c)
      ELSE
        qb(2) = -DSQRT((s-1.d0)*grav*d50*d50*d50)*an*DSQRT(theta_cy)   &
                 *theta_c*DEXP(-b*theta_cr/theta_c)
      ENDIF
    ENDIF
  ENDIF

  !- put zero if NaN value for qb (temporary solution)
!  IF ((qb(1) .NE. qb(1)) .OR. (qb(2) .NE. qb(2))) THEN
!    qb(:) = 0.d0
!  ENDIF

  !--- Suspended load transport -------------------------------------

  IF (hs>0.AND.tp>0.AND.wlpeak>0) THEN
    !- net mean current after a wave period
    Uc_net = Ucw_on - Ucw_off
  ENDIF

  !- bed reference sediment concentration
  AcR = 0.0035d0*DEXP(-0.3d0*D_star)
  IF (hs>0.AND.tp>0.AND.wlpeak>0) THEN
    IF (theta_cw .GT. 0.d0) THEN
      cR = AcR*theta_cwm*DEXP(-4.5d0*theta_cr/theta_cw)
    ELSE
      cR = 0.d0
    ENDIF
  ELSE
    IF (theta_c .GT. 0.d0) THEN
      cR = AcR*theta_c*DEXP(-4.5d0*theta_cr/theta_c)
    ELSE
      cR = 0.d0
    ENDIF
  ENDIF

  !--- calcul of sediment diffusivity

  !-- sediment diffusivity due to current alone
  ustar_c = DSQRT(tau_c/rho0)
  D_c = tau_c*ustar_c  !energy dissipation from bottom friction due to current
  !- Schmidt nb calculation
  IF (ustar_c.GT.0) THEN
    IF ((Ws/ustar_c) .LE. 1.d0) THEN
      sigma_c = 0.4d0 + 3.5d0*(DSIN(pi*Ws/2.d0/ustar_c))**2.d0
    ELSE
      sigma_c = 1.d0 + 2.9d0*(DSIN(pi*ustar_c/2.d0/Ws))**2.d0
      ! sigma_c = 1.d0 + 2.9d0*(DSIN(pi*Ws/2.d0/ustar_c))**2.d0
      ! sigma_c = 1.d0 + 2.9d0*(DSIN(pi*ustar_c/2.d0/Ws))**2.5d0
    ENDIF
    k_c = kappa*sigma_c/6.d0
  ELSE
    k_c = 0.d0
  ENDIF

  IF (hs>0.AND.tp>0.AND.wlpeak>0) THEN
    !-- sediment diffusivity due to waves alone
    ustar_w = DSQRT(tau_w/rho0)
    D_w = tau_w*ustar_w   !energy dissipation from bottom friction due to waves                    
    !- Schmidt nb calculation
    IF ((Ws/ustar_w) .LE. 1.d0) THEN
      if(ustar_w==0) call parallel_abort('cl11_tot: ustar_w==0')
      sigma_w = 0.15d0 + 1.5d0*(DSIN(pi*Ws/2.d0/ustar_w))**2.d0
    ELSE
      sigma_w = 1.d0 + 0.65d0*(DSIN(pi*ustar_w/2.d0/Ws))**2.d0
      ! sigma_w = 1.d0 + 0.65d0*(DSIN(pi*Ws/2.d0/ustar_w))**2.d0
      ! sigma_w = 1.d0 + 0.65d0*(DSIN(pi*ustar_w/2.d0/Ws))**2.5d0
    ENDIF

    !- sediment diffusivity due to waves and current
    if(theta_c+theta_w==0) call parallel_abort('cl11_tot: theta_c+theta_w==0 (2)')
    Y = theta_c/(theta_c+theta_w)
    sigma_cw = Y*sigma_c + (1.d0-Y)*sigma_w
    k_cw_c = kappa*sigma_cw/6.d0
    k_cw_w = kappa*sigma_cw/3.d0/pi

    !-- sediment diffusivity due to wave breaking
    !- calculation of the slope along the local wave direction
    dpdxy_tot = DSQRT(dpdxy(1)*dpdxy(1) + dpdxy(2)*dpdxy(2))
    slope_dir = DATAN2(-dpdxy(2),-dpdxy(1))
    slope_wdir_angle = DABS(wdir-slope_dir)
    IF (slope_wdir_angle .GT. pi) THEN
      slope_wdir_angle = 2.d0*pi - slope_wdir_angle
    ENDIF
    dpdxy_wdir = DCOS(slope_wdir_angle)*dpdxy_tot

    !- irribaren nb
    if(hs<0.or.wlpeak<=0) call parallel_abort('cl11_tot: hs<0')
    IF (dpdxy_wdir .GT. 0.d0) THEN
      irribaren = dpdxy_wdir/DSQRT(hs/wlpeak)
    ELSE
      irribaren = 0.001d0/DSQRT(hs/wlpeak)
    ENDIF

    !- calculation of energy dissipation from wave breaking
    Hrms = hs/DSQRT(2.d0)
    if(Hrms==0) call parallel_abort('cl11_tot: Hrms=0')
#ifdef USE_WWM
    a_b = DEXP(-(gammab*htot/Hrms)**2.d0)
#else
    a_b = 1.d0
#endif
    A_epsi = 2.d0*DTANH(5.d0*irribaren)
    tmp=tp*(4.d0*htot*htot-hs*hs)
    if(tmp==0) call parallel_abort('cl11_tot: tmp=0 (2)')
    D_b = a_b*A_epsi*rho0*grav*htot*hs**3.d0/tmp
!        / (tp*(4.d0*htot*htot-hs*hs))

    !- parametrization of the efficiency parameter for wave breaking (see Camenen 2007, p. 132)
    k_b = 0.01d0
    ! k_b = 0.012d0*(1.d0+DTANH(irribaren))
    ! k_b = 0.062d0*(1.d0-0.9d0*DTANH(0.25d0*ustar_w/Ws))

    !-- total sediment diffusivity
    epsi = htot*((D_b*k_b**3.d0 + D_c*k_cw_c**3.d0 + D_w*k_cw_w**3.d0) &
           /rho0)**(1.d0/3.d0)

    !--- suspended load transport components
    !- parallel and normal to the wave direction
    qs_w = Uc_net*cR*epsi/Ws * (1.d0-DEXP(-Ws*htot/epsi))
    qs_n = Uc_n*cR*epsi/Ws * (1.d0-DEXP(-Ws*htot/epsi))

    !- x,y components
    qs(1) = qs_w*DCOS(wdir) - qs_n*DSIN(wdir)
    qs(2) = qs_w*DSIN(wdir) + qs_n*DCOS(wdir)

  ELSE
    !-- total sediment diffusivity
    epsi = htot * ((D_c*k_c**3.d0)/rho0)**(1.d0/3.d0)

    !--- suspended load transport components
    IF (epsi.GT.0.d0) THEN
      qs(1) = Uc_x*cR*epsi/Ws * (1.d0-DEXP(-Ws*htot/epsi))
      qs(2) = Uc_y*cR*epsi/Ws * (1.d0-DEXP(-Ws*htot/epsi))
    ELSE
      qs(1) = 0.d0
      qs(2) = 0.d0
    ENDIF
  ENDIF

  !- put zero if NaN value for qs (temporary solution)
!  IF ((qs(1) .NE. qs(1)) .OR. (qs(2) .NE. qs(2))) THEN
!    qs(:) = 0.d0
!  ENDIF

  END SUBROUTINE cl11_tot

!--------------------------------------------------------------------
  SUBROUTINE wl14_tot(d,dclass,d50,d90,D_star,F_class,htot,dpdxy,    &
                     &u,v,hs,tp,wdir,wlpeak,qb,qs)
!--------------------------------------------------------------------
! This subroutine computes bedload and suspended load sediment
! transport rate (m^2/s) from Wu and Lin (2014) which is a formulation
! created for nonuniform (i.e. multi-class) sediment transport
!
! Author: thomas guerin (thomas.guerin@univ-lr.fr)    
! Date: 12/2014
!
! History:
!--------------------------------------------------------------------

  USE schism_glbl, ONLY: grav,pi,rho0,rhosed,rkind
  USE schism_msgp, ONLY: parallel_abort
  USE sed2d_mod, ONLY: iasym,nb_class,s,wvisco

  IMPLICIT NONE

!- Arguments --------------------------------------------------------  
  REAL(rkind), INTENT(IN) :: d,d50,d90,D_star,htot,u,v,hs,tp,wdir,   &
                            &wlpeak
  REAL(rkind), DIMENSION(2), INTENT(IN) :: dpdxy
  REAL(rkind), DIMENSION(nb_class), INTENT(IN) :: dclass,F_class
  REAL(rkind), DIMENSION(2), INTENT(OUT) :: qb,qs
!- Constants --------------------------------------------------------  
  REAL(rkind), PARAMETER :: Ar = 12.d0
  REAL(rkind), PARAMETER :: shields_cr = 0.03d0
  REAL(rkind), PARAMETER :: m = 0.6d0
  INTEGER, PARAMETER :: ech = 50
!- Local variables --------------------------------------------------
  REAL(rkind) :: ac,at,Aw,Cd,Deltar_c,Deltar_c_mm,Deltar_w,f_grain_c,&
                &f_grain_cw,f_grain_w,fw,kpeak,ks_c,ks_form_c,       &
                &ks_form_w,ks_grain,ks_w,Lambdar_c,Lambdar_c_mm,     &
                &Lambdar_w,n,n_grain,omega,pe,ph,phi,qb_off,qb_on,   &
                &qs_c,rw,T_crest,T_trough,tau_b,tau_b_c,tau_b_wm,    &
                &tau_cr,tau_grain_b_off,tau_grain_b_on,              &
                &tau_grain_b_wm_off,tau_grain_b_wm_on,theta_off,     &
                &theta_on,tmp,U_crest,U_trough,Uc,udir,Uw,Uwm_off,   &
                &Uwm_on,Ws,Xu
  REAL(rkind), DIMENSION(ech) :: uorb
  INTEGER :: i
!--------------------------------------------------------------------

!- current amplitude and direction
  Uc = dsqrt(u*u+v*v)
  udir = datan2(v,u)

!- angle between current and waves directions
  phi = wdir - udir

!- Wave parameters
  IF (hs>0.d0.and.tp>0.d0.and.wlpeak>0.d0) THEN
    !- wave angular frequency
    omega = 2.d0*pi/tp !>0

    !- peak wave number
    kpeak = 2.d0*pi/wlpeak !>0
  ENDIF

!- Wave orbital velocity calculation
  IF (hs>0.d0.and.tp>0.d0.and.wlpeak>0.d0) THEN
    IF (iasym==0) THEN
      !- Airy waves
      tmp=dsinh(kpeak*htot)
      if(tmp==0.d0) call parallel_abort('wl14_tot(1)')
      Uw = dmax1(pi*hs/tp/tmp,0.001d0)
      U_crest = Uw !>0
      T_crest = tp/2.d0 !>0
      T_trough = T_crest !>0

    ELSEIF (iasym==1) THEN
      !- Wave asymmetry treatment : compute wave orbital velocity
      !  (horizontal component) from Elfrink (2006)
      CALL wave_asymmetry_Elfrink(hs,tp,wdir,htot,-dpdxy(1),-dpdxy(2),&
                          &ech,uorb,U_crest,U_trough,T_crest,T_trough)

      Uw = (U_crest+U_trough)/2.d0 ! >0

      IF (Uw/=Uw.or.U_crest<=0.or.U_trough<=0.or.T_crest<=0.or.T_trough<=0) THEN
        write(12,*)'Airy instead of Elfrink'
        write(12,*)'hs',hs,'tp',tp,'wdir',wdir,'wlpeak',wlpeak,'htot',htot
!        call parallel_abort('wl14_tot(2)')
        !- Airy waves
        tmp=dsinh(kpeak*htot)
        if(tmp==0.d0) call parallel_abort('wl14_tot(2)')
        Uw = pi*hs/tp/tmp !>0
        U_crest = Uw !>0
        T_crest = tp/2.d0 !>0
        T_trough = T_crest !>0
      ENDIF

    ENDIF !iasym
  ELSE !no waves
    Uw = 0.d0
  ENDIF

!- Bed-load transport

  !exposed probability of particles
  pe = 0.d0
  do i=1,nb_class
     pe = pe + F_class(i)*d/(d+dclass(i))
  enddo

  !hidden probability of particles
  ph = 1.d0 - pe

  !critical shear stress for the incipient motion of corresponding 
  !size class
  tau_cr = grav*(rhosed-rho0)*d*shields_cr*(pe/ph)**(-m) ! >0
  !NB: grav forgotten in Wu and Lin 2014 (Eq. 24)

  !equivalent roughness height calculation
  ks_grain = 1.5d0*d90 !grain roughness height
  Lambdar_c = 1000.d0*d50 !ripple wavelength (m) according to Soulsby (97)
  Deltar_c = Lambdar_c/7.d0 !ripple height (m) according to Soulsby (97)
!  Lambdar_c_mm = 245.d0*(d50*1000.d0)**0.35d0 !ripple wavelength (mm)
!  Lambdar_c = Lambdar_c_mm/1000.d0 !ripple wavelength (m)
!  Deltar_c_mm = 0.074d0*Lambdar_c_mm*(d50*1000.d0)**(-0.253d0) !ripple height (mm)
!  Deltar_c = Deltar_c_mm/1000.d0 !ripple height (m)
  ks_form_c = Ar*Deltar_c**2.d0/Lambdar_c ! >0
  ks_c = ks_grain + ks_form_c !equivalent roughness height

  !grain friction coefficient under current only
  n_grain = 1.d0/20.d0*d50**(1.d0/6.d0) !Manning coef of grain roughness
  n = htot**(1.d0/6.d0)/(18.d0*dlog(12.d0*htot/ks_c)) !Manning coef of bed roughness
  f_grain_c = 2.d0*(n_grain/n)**(3.d0/2.d0)*grav*n**2.d0             &
             &/(htot**(1.d0/3.d0))

  !grain friction coefficient under waves only
  Aw = Uw*tp/2.d0/pi !wave excursion
  if (Aw>0.d0) then
    f_grain_w = 0.237d0*(Aw/ks_grain)**(-0.52d0)
  else
    f_grain_w = 0.d0
  endif

  !grain friction coefficient under combined waves and current
  if (Uc>0.d0) then
    Xu = Uc**2.d0/(Uc**2.d0 + 0.5d0*Uw**2.d0)
  else
    Xu = 0.d0
  endif
  f_grain_cw = Xu*f_grain_c + (1.d0-Xu)*f_grain_w

  !bed grain shear stress averaged over the two half cycles
  if (hs>0.d0.and.tp>0.d0) then
    rw = U_crest/Uw - 1.d0 !wave asymmetry coef
    ac = pi*T_crest/tp
    at = pi*T_trough/tp

    tmp = 0.5d0*rho0*f_grain_w*Uw**2.d0/2.d0

    tau_grain_b_wm_on = tmp*(1.d0 + rw**2.d0 +                       &
      &13.d0/6.d0*rw*dsin(ac)/ac + 1.d0/6.d0*dsin(2.d0*ac)/(2.d0*ac))

    tau_grain_b_wm_off = tmp*(1.d0 + rw**2.d0 -                      &
      &13.d0/6.d0*rw*dsin(at)/at + 1.d0/6.d0*dsin(2.d0*at)/(2.d0*at))
  else
    tau_grain_b_wm_on = 0.d0
    tau_grain_b_wm_off = 0.d0
  endif

  !root-mean-square values of wave orbital velocity over the two half cycles
  if (hs>0.d0.and.tp>0.d0) then
    tmp = 0.5d0*rho0*f_grain_w ! >0
    if (tmp==0.d0) then
      write(12,*)'f_grain_w',f_grain_w,'Aw',Aw,'Uw',Uw,'tp',tp,'ks_grain',ks_grain
      call parallel_abort('wl14_tot(3)')
    endif
    Uwm_on = dsqrt(tau_grain_b_wm_on/tmp)
    Uwm_off = dsqrt(tau_grain_b_wm_off/tmp)
  else
    Uwm_on = 0.d0
    Uwm_off = 0.d0
  endif

  !bed grain shear stress due to combined waves and current
  tmp = 0.5d0*rho0*f_grain_cw
  tau_grain_b_on = tmp*(Uc**2.d0 + Uwm_on**2.d0                      &
                       & + 2.d0*Uc*Uwm_on*dcos(phi))
  tau_grain_b_off = tmp*(Uc**2.d0 + Uwm_off**2.d0                    &
                       & + 2.d0*Uc*Uwm_off*dcos(pi-phi))

  !resultant bed-load transport components during the two half cycles
  tmp = 0.0053d0*dsqrt((s-1.d0)*grav*d**3.d0)
  if (tau_grain_b_on>tau_cr) then
    qb_on = tmp*(tau_grain_b_on/tau_cr - 1.d0)**2.2d0
  else
    qb_on = 0.d0
  endif
  if (tau_grain_b_off>tau_cr) then
    qb_off = tmp*(tau_grain_b_off/tau_cr - 1.d0)**2.2d0
  else
    qb_off = 0.d0
  endif

  !angle of the resultant components
  theta_on = datan2(Uc*dsin(udir)+Uwm_on*dsin(udir+phi),             &
                   &Uc*dcos(udir)+Uwm_on*dcos(udir+phi))

  theta_off = datan2(Uc*dsin(udir)+Uwm_off*dsin(udir-pi+phi),        &
                    &Uc*dcos(udir)+Uwm_off*dcos(udir-pi+phi))

  !final x and y components of bed-load transport
  if (hs>0.d0.and.tp>0.d0) then
    qb(1) = T_crest/tp*dcos(theta_on)*qb_on +                        &
           &T_trough/tp*dcos(theta_off)*qb_off

    qb(2) = T_crest/tp*dsin(theta_on)*qb_on +                        &
           &T_trough/tp*dsin(theta_off)*qb_off
  else !only current
    qb(1) = dcos(udir)*qb_on
    qb(2) = dsin(udir)*qb_on
  endif

!- Suspended-load transport

  !settling velocity
  Ws = wvisco/d*(dsqrt(10.36d0**2.d0+1.049d0*D_star**3.d0)-10.36d0)

  !bed shear stress due to current
  Cd = grav*n**2.d0/(htot**(1.d0/3.d0))
  tau_b_c = rho0*Cd*Uc**2.d0

  !equivalent roughness height calculation
  if (hs>0.d0.and.tp>0.d0) then
    Lambdar_w = Aw/(1.d0 + 0.00187d0*Aw/d50*                         &
               &(1.d0-dexp(-(0.0002d0*Aw/d50)**1.5d0))) !ripple wavelength (m)
    Deltar_w = 0.15d0*Lambdar_w*(1.d0 - dexp(-(5000.d0*d50/Aw)**3.5d0)) !ripple height (m)
    ks_form_w = Ar*Deltar_w**2.d0/Lambdar_w !form roughness height
  else
    ks_form_w = 0.d0
  endif
  ks_w = ks_grain + ks_form_w !equivalent roughness height

  !mean bed wave stress averaged over a wave cycle
  if (Aw>0.d0) then
    fw = 0.237d0*(Aw/ks_w)**(-0.52d0)
  else
    fw = 0.d0
  endif
  tau_b_wm = 1.d0/4.d0*rho0*fw*Uw**2.d0

  !total bed shear stress
  tau_b = dsqrt(tau_b_c**2.d0 + tau_b_wm**2.d0 +                     &
               &2.d0*tau_b_c*tau_b_wm*dcos(phi))

  !transport along the current direction
  if (tau_b>tau_cr) then
    qs_c = 0.0000262*dsqrt((s-1.d0)*grav*d**3.d0)*                   &
                      &((tau_b/tau_cr-1.d0)*Uc/Ws)**1.74d0
  else
    qs_c = 0.d0
  endif

  !final x and y components of suspended-load transport
  qs(1) = dcos(udir)*qs_c
  qs(2) = dsin(udir)*qs_c

  END SUBROUTINE wl14_tot

!--------------------------------------------------------------------
  SUBROUTINE tr04_tot(d_class,d50,d90,htot,dpdxy,u,v,hs,tp,   &
                     &wdir,wlpeak,qb_class,qs_class)
!--------------------------------------------------------------------
! This subroutine computes bedload and suspended load sediment
! transport rate (m^2/s) from TRANSPOR2004 model (Final report Z3748,
! Delft Hydraulics) developed for multi-sized sediment transport 
! under combined current and waves
!
! Author: thomas guerin (thomas.guerin@univ-lr.fr)    
! Date: 02/2016
!
! History:
!--------------------------------------------------------------------

  USE schism_glbl, ONLY: grav,pi,rho0,rhosed,rkind
  USE schism_msgp, ONLY: parallel_abort
  USE sed2d_mod, ONLY: iasym,s,wvisco

  IMPLICIT NONE

!- Arguments --------------------------------------------------------  
  REAL(rkind), INTENT(IN) :: d_class,d50,d90,htot,u,v,hs,tp,wdir,wlpeak
  REAL(rkind), DIMENSION(2), INTENT(IN) :: dpdxy
  REAL(rkind), DIMENSION(2), INTENT(OUT) :: qb_class,qs_class
!- Constants --------------------------------------------------------  
  REAL(rkind), PARAMETER :: kappa = 0.4d0
  REAL(rkind), PARAMETER :: dz = 0.02d0        !m
  REAL(rkind), PARAMETER :: hs_min = 0.05d0    !m
  REAL(rkind), PARAMETER :: tp_min = 2.0d0     !sec
  REAL(rkind), PARAMETER :: wlpeak_min = 1.d0  !m
  REAL(rkind), PARAMETER :: p_mud = 0.d0       !mud fraction: not taken into account
  INTEGER, PARAMETER :: ech = 25
  INTEGER, PARAMETER :: ndz = 50
!- Local variables --------------------------------------------------
  REAL(rkind) :: Ad,alpha_cw,alpha_v,b,b_bis,b_max,B1,B2,beta,       &
                &beta_c_class,beta_d_class,beta_f,    &
                &beta_w_class,c_a_class,c_z_dw_class,Chezy,cphase,   &
                &D_star_class,dc_dz_class,dm,ds,dt,dw, &
                &epsi_s_w_bed_class,epsi_w_max,F,f_c,f_grain_c,      &
                &f_grain_cw,f_grain_w,f_w,gama,gama_br,k_a,kpeak,    &
                &ks_c,ks_c_d,ks_c_mr,ks_c_r,ks_w_r,lambda_class,mob, &
                &mu_c,mu_w_class,omega,phi,qs_c_class,qs_c_class_x,  &
                &qs_c_class_y,qs_w_class,qs_w_class_x,qs_w_class_y,  &
                &T_for,T_back,tau_adim_cw_class,tau_b_c,tau_b_cw,    &
                &tau_b_w,tau_cr,tau_cr_1,tau_grain_b_cw_class,       &
                &theta_cr,tmp,u_d,U_d_back,U_d_for,U_d_w_r,u_star_c, &
                &u_star_cw,u_star_w,Uasym,Uc,Uc_a,Uc_a_x,Uc_a_y,     &
                &Uc_dm,Uc_x,Uc_y,udir,Uw,Uwc,Ws_class,xi_class,z_a
  REAL(rkind), DIMENSION(ech) :: beta_slope_t,f_slope1_t,f_slope2_t, &
                                &uw_t,qb_t_class,qb_x_t_class,       &
                                &qb_y_t_class,t,tau_adim_cw_t_class, &
                                &tau_grain_b_cw_t,u_grain_star_cw_t, &
                                &Uvec_a_x_t,Uvec_a_t,Uvec_a_y_t,     &
                                &uw_x_t,uw_y_t
  REAL(rkind), DIMENSION(ndz) :: c_z_class,epsi_s_c_class,           &
                                &epsi_s_cw_class,epsi_s_w_class,Uc_z,z                                 
  INTEGER :: i
!--------------------------------------------------------------------

!- current amplitude and direction
  Uc = dsqrt(u*u+v*v)
  udir = datan2(v,u)

  if(Uc/=Uc) call parallel_abort('Uc = NaN')
  if(udir/=udir) call parallel_abort('udir = NaN')

! angle between wave and current directions
  phi = wdir-udir

  if(phi/=phi) call parallel_abort('phi = NaN')

! dimensionless diameter
  D_star_class = d_class*(grav*(s-1.d0)/wvisco**2.d0)**(1.d0/3.d0) !>0

  if(D_star_class/=D_star_class) call parallel_abort('D_star_class = NaN')

! initiation of motion (critical shields)
  if (D_star_class<=1.d0) then
      call parallel_abort('tr04_tot_multi_size: D_star_class<=1')
      
  elseif (D_star_class>1.d0.and.D_star_class<=4.d0) then
      theta_cr = 0.115d0*D_star_class**(-0.5d0)
    
  elseif (D_star_class>4.d0.and.D_star_class<=10.d0) then
      theta_cr = 0.14d0*D_star_class**(-0.64d0)
    
  elseif (D_star_class>10.d0.and.D_star_class<=20.d0) then
      theta_cr = 0.04d0*D_star_class**(-0.1d0)
    
  elseif (D_star_class>20.d0.and.D_star_class<=150.d0) then
      theta_cr = 0.013d0*D_star_class**0.29d0
    
  elseif (D_star_class>150.d0) then
      theta_cr = 0.055d0
  endif

  if(theta_cr/=theta_cr) call parallel_abort('theta_cr = NaN')

! critical bed shear stress
  tau_cr = (rhosed-rho0)*grav*d50*theta_cr
  tau_cr_1 = (1.d0+p_mud)**3.d0*tau_cr ! NB: (1+p_mud) or (1-p_mud) ?

  if(tau_cr/=tau_cr) call parallel_abort('tau_cr = NaN')
  if(tau_cr_1/=tau_cr_1) call parallel_abort('tau_cr_1 = NaN')

! time vector corresponding to a wave period
  t(:) = 0.d0
  if (hs>hs_min.and.tp>tp_min.and.wlpeak>wlpeak_min) then
    dt = tp/ech
    do i=1,ech
       t(i) = dt*(i-1)
    enddo
  endif

! wave orbital velocity
  uw_t(:) = 0.d0
  if (hs>hs_min.and.tp>tp_min.and.wlpeak>wlpeak_min) then
    
    ! wave frequency
    omega = 2.d0*pi/tp  !>0
    
    ! peak wave number
    kpeak = 2.d0*pi/wlpeak !>0
    
    ! phase velocity
    cphase = omega/kpeak
    
    if (iasym==0) then ! Airy waves
        
        tmp=dsinh(kpeak*htot)
        if(tmp==0.d0) call parallel_abort('tr04_tot_multi_size: sinh(kpeak*htot) = 0')
        
        Uw = pi*hs/tp/tmp
        U_d_for = Uw
        U_d_back = Uw ! NB: U_d_back > 0
        
        do i=1,ech
           uw_t(i) = Uw*dsin(omega*t(i))
        enddo
        
    elseif (iasym==1) then ! asymmetric waves (Elfrink, 2006)
        
        if (dabs(dpdxy(1))>0.d0) then
            CALL wave_asymmetry_Elfrink(hs,tp,wdir,htot,-dpdxy(1),-dpdxy(2),&
                          &ech+1,uw_t,U_d_for,U_d_back,T_for,T_back)
        else
            CALL wave_asymmetry_Elfrink(hs,tp,wdir,htot,0.001d0,-dpdxy(2),&
                          &ech+1,uw_t,U_d_for,U_d_back,T_for,T_back)
        endif
        Uw = 0.5d0*(U_d_for+U_d_back)
    endif
  else
    omega = 0.d0
    kpeak = 0.d0
    cphase = 0.d0
    Uw = 0.d0
    U_d_for = 0.d0
    U_d_back = 0.d0
  endif

  if(U_d_for/=U_d_for) call parallel_abort('U_d_for = NaN')
  if(U_d_back/=U_d_back) call parallel_abort('U_d_back = NaN')

! near-bed peak orbital excursion (linear waves)
  tmp = dsinh(kpeak*htot)
  if (hs>hs_min.and.tmp>0.d0) then
    Ad = hs/(2.d0*tmp)
  else
    Ad = 0.d0
  endif

  if(Ad/=Ad) call parallel_abort('Ad = NaN')

! representative peak orbital velocity
  U_d_w_r = (0.5d0*U_d_for**3.d0 + 0.5d0*U_d_back**3.d0)**(1.d0/3.d0)

  if(U_d_w_r/=U_d_w_r) call parallel_abort('U_d_w_r = NaN')

! mobility parameter
  Uwc = dsqrt(U_d_w_r**2.d0 + Uc**2.d0 + 2.d0*U_d_w_r*Uc*dcos(phi)) ! method 1
!  Uwc = dsqrt(U_d_w_r**2.d0 + Uc_a**2.d0 + 2.d0*U_d_w_r*Uc_a*dcos(phi)) ! method 2
  mob = Uwc**2.d0/((s-1.d0)*grav*d50)

  if(Uwc/=Uwc) call parallel_abort('Uwc = NaN')
  if(mob/=mob) call parallel_abort('mob = NaN')

! sediment fall velocity
  Ws_class = wvisco/d_class*(dsqrt(10.36d0**2.d0+1.049d0*D_star_class**3.d0)-10.36d0)
!  if (d_class>0.000065d0.and.d_class<=0.0001d0) then
!    Ws_class = (s-1.d0)*grav*d_class**2.d0/(18.d0*wvisco)
!  elseif (d_class>0.0001d0.and.d_class<=0.001d0) then
!    Ws_class = 10.d0*wvisco/d_class*        &
!      &(dsqrt(1.d0+0.01d0*(s-1.d0)*grav*d_class**3.d0/wvisco**2.d0)-1.d0)
!  elseif (d_class>0.001d0) then
!    Ws_class = 1.1d0*dsqrt((s-1.d0)*grav*d_class)
!  endif

  if(Ws_class/=Ws_class) call parallel_abort('Ws_class = NaN')

! current and wave-related bed roughness due to ripples
  if (mob<=50.d0) then
    ks_c_r = 150.d0*d50
  elseif (mob>=250.d0) then
    ks_c_r = 20.d0*d50
  elseif (mob>50.d0.and.mob<250.d0) then
    ks_c_r = (182.5d0-0.65d0*mob)*d50
  endif
  ks_c_r = min(ks_c_r,0.075d0)
  ks_w_r = ks_c_r ! NB: ks_c_r and ks_w_r are the same in TR2004 (make sense?)

  if(ks_c_r/=ks_c_r) call parallel_abort('ks_c_r = NaN')

! current-related bed roughness due to mega-ripples
  if (htot>1.d0) then
    if (mob<=50.d0) then
        ks_c_mr = 0.0002d0*mob*htot
    elseif (mob>=550.d0) then
        ks_c_mr = 0.d0
    elseif (mob>50.d0.and.mob<550.d0) then
        ks_c_mr = (0.011d0-0.00002d0*mob)*htot
    endif
    ks_c_mr = min(ks_c_mr,0.2d0)
  else ! htot<=1
    ks_c_mr = 0.d0
  endif

  if(ks_c_mr/=ks_c_mr) call parallel_abort('ks_c_mr = NaN')

! current-related bed roughness due to dunes (only rivers)
  if (hs<hs_min.and.htot>1.d0) then
    if (mob<=100.d0) then
        ks_c_d = 0.0004d0*mob*htot
    elseif (mob>=600.d0) then
        ks_c_d = 0.d0
    elseif (mob>100.d0.and.mob<600.d0) then
        ks_c_d = (0.048d0-0.00008d0*mob)*htot
    endif
    ks_c_d = min(ks_c_d,1.d0)
  else
    ks_c_d = 0.d0
  endif

  if(ks_c_d/=ks_c_d) call parallel_abort('ks_c_d = NaN')

! overall current-related bed roughness
  ks_c = dsqrt(ks_c_r**2.d0 + ks_c_mr**2.d0 + ks_c_d**2.d0)
  if(ks_c==0.d0) call parallel_abort('tr04_tot_multi_size: ks_c = 0')

  if(ks_c/=ks_c) call parallel_abort('ks_c = NaN')

! apparent bed roughness
  if (dabs(phi)<=pi) then
    beta = dabs(phi)
  else
    beta = 2.d0*pi-dabs(phi)
  endif
  gama = 0.8d0 + beta - 0.3d0*beta**2.d0
  if (Uc>0.d0) then
    k_a = ks_c*dexp(gama*U_d_w_r/Uc) ! method 1
  else
    k_a = ks_c
  endif
!  if (Uc_a>0.d0) then
!      k_a = ks_c*dexp(gama*U_d_w_r/Uc_a) ! method 2
!  endif
  k_a = min(k_a,10.d0*ks_c)
  if (k_a<ks_c) then
    k_a = ks_c
  endif

  if(k_a/=k_a) call parallel_abort('k_a = NaN')

! thickness of wave boundary layer (m)
  if (ks_w_r>0.d0.and.Ad>0.d0) then
    dw = 0.36d0*Ad*(Ad/ks_w_r)**(-0.25d0)
  else
    dw = 0.d0
  endif  

! reference level (m)
  z_a = max(0.5d0*ks_c_r,0.5d0*ks_w_r)
  z_a = max(z_a,0.01d0)
  z_a = min(z_a,0.2d0*htot)

! thickness of effective fluid mixing layer (m)
  dm = 2.d0*dw
  if (dm<0.05d0) then
    dm = 0.05d0
  elseif (dm>0.2d0) then
    dm = 0.2d0
  endif

! current velocity at the top of effective fluid mixing layer
  Uc_dm = Uc*dlog(30.d0*dm/k_a)/(-1.d0+dlog(30.d0*htot/k_a))
  
!  if(Uc_dm/=Uc_dm) call parallel_abort('Uc_dm = NaN')
  if(Uc_dm/=Uc_dm) write(12,*) 'Uc=',Uc,'dm=',dm,'k_a=',k_a,'htot=',htot,'dw',dw,'Ad',Ad,'ks_w_r',ks_w_r

! near bed current velocity
  if (z_a>=dm) then ! make sense ?
    Uc_a = Uc*dlog(30.d0*z_a/k_a)/(-1.d0+dlog(30.d0*htot/k_a))
  else ! make sense ?
    Uc_a = Uc_dm*dlog(30.d0*z_a/ks_c)/(dlog(30.d0*dm/ks_c))
  endif
  if(Uc_a<0.d0) Uc_a = 0.d0

  if(Uc_a/=Uc_a) call parallel_abort('Uc_a = NaN')

! x and y current velocity
  Uc_x = Uc*dcos(udir)
  Uc_y = Uc*dsin(udir)

! x and y near-bed current velocity
  ! method 1
  Uc_a_x = Uc_a*dcos(udir)
  Uc_a_y = Uc_a*dsin(udir)
  ! method 2
!  Uc_a_x = Uc_dm*dcos(udir)
!  Uc_a_y = Uc_dm*dsin(udir)

! x and y instantaneous near-bed orbital velocity
  uw_x_t(:) = 0.d0
  uw_y_t(:) = 0.d0
  do i=1,ech
     uw_x_t(i) = uw_t(i)*dcos(wdir)
     uw_y_t(i) = uw_t(i)*dsin(wdir)
  enddo

! instantaneous near-bed global velocity
  Uvec_a_x_t(:) = 0.d0
  Uvec_a_y_t(:) = 0.d0
  Uvec_a_t(:) = 0.d0
  do i=1,ech
     Uvec_a_x_t(i) = Uc_a_x + uw_x_t(i)
     Uvec_a_y_t(i) = Uc_a_y + uw_y_t(i)
     Uvec_a_t(i) = dsqrt(Uvec_a_x_t(i)**2.d0 + Uvec_a_y_t(i)**2.d0)
  enddo  

! instantaneous slope effect
  beta_slope_t(:) = 0.d0
  f_slope1_t(:) = 0.d0
  do i=1,ech
     if (Uvec_a_t(i)>0.d0) then
       beta_slope_t(i) = -dpdxy(1)*(Uvec_a_x_t(i)/Uvec_a_t(i)) + &
                        &-dpdxy(2)*(Uvec_a_y_t(i)/Uvec_a_t(i))
     else
       beta_slope_t(i) = 0.d0
     endif

     tmp = 1.d0 + beta_slope_t(i)/0.6d0
     if (tmp>0.d0) then
       f_slope1_t(i) = tmp**(-1.d0)
     else
       f_slope1_t(i) = 1.d0
     endif
  enddo

! grain-related friction factor coefficient for current and waves
  tmp = dlog10(12.d0*htot/d90)
  if(tmp<=0.d0) call parallel_abort('pb f_grain_c')
  f_grain_c = 0.24d0*(dlog10(12.d0*htot/d90))**(-2.d0)

  tmp = 5.2d0*(Ad/d90)
  if (tmp>0.d0) then
    f_grain_w = dexp(-6.d0+5.2d0*(Ad/d90)**(-0.19d0))
    f_grain_w = min(f_grain_w,0.05d0)
  else
    f_grain_w = 0.05d0
  endif

! grain-related friction factor coefficient for combined current and waves
  tmp = U_d_for+Uc_a
  if (tmp>0.d0) then
    alpha_v = Uc_a/tmp
  else
    alpha_v = 0.d0
  endif
  beta_f = 0.25d0*(-1.d0+dlog(30.d0*htot/ks_c))**2.d0/(dlog(30.d0*z_a/ks_c))**2.d0
  f_grain_cw = dsqrt(alpha_v)*beta_f*f_grain_c + (1.d0-dsqrt(alpha_v))*f_grain_w

  if(beta_f/=beta_f) call parallel_abort('beta_f = NaN')
  if(f_grain_cw/=f_grain_cw) call parallel_abort('f_grain_cw = NaN')

! instantaneous grain-related bed-shear stress and velocity
  tau_grain_b_cw_t(:) = 0.d0
  u_grain_star_cw_t(:) = 0.d0
  do i=1,ech
     tau_grain_b_cw_t(i) = 0.5d0*rho0*f_grain_cw*Uvec_a_t(i)**2.d0
     u_grain_star_cw_t(i) = dsqrt(tau_grain_b_cw_t(i)/rho0)
  enddo

! correction factor of excess bed-shear stress related to grain roughness effects
  lambda_class = dsqrt(d_class/d50)

! hiding factor of Egiazaroff
  xi_class = (dlog10(19.d0)/dlog10(19.d0*d_class/d50))**2.d0

  if(xi_class/=xi_class) call parallel_abort('xi_class = NaN')

! instantaneous bed-shear stress
  f_slope2_t(:) = 0.d0
  tau_adim_cw_t_class(:) = 0.d0
  do i=1,ech
     f_slope2_t(i) = dsin(datan(0.6d0)+datan(beta_slope_t(i)))/dsin(datan(0.6d0))
     if(f_slope2_t(i)/=f_slope2_t(i)) call parallel_abort('f_slope2_t(i) = NaN')

     tmp = (tau_grain_b_cw_t(i) - f_slope2_t(i)*d_class/d50*xi_class*tau_cr_1)
     if (tmp>0.d0) then
        tau_adim_cw_t_class(i) = lambda_class*tmp/(tau_cr*d_class/d50)
     else
        tau_adim_cw_t_class(i) = 0.d0
     endif
  enddo

! instantaneous bed-load transport (kg.m-1.s-1)
  qb_t_class(:) = 0.d0
  qb_x_t_class(:) = 0.d0
  qb_y_t_class(:) = 0.d0
  do i=1,ech
     qb_t_class(i) = 0.5d0*rhosed*f_slope1_t(i)*d_class*   &
                    &u_grain_star_cw_t(i)*tau_adim_cw_t_class(i)/D_star_class**0.3d0

     if(qb_t_class(i)<0) qb_t_class(i)=0
     if(qb_t_class(i)/=qb_t_class(i)) call parallel_abort('qb_t_class(i) = NaN')
  
     if (dabs(Uvec_a_t(i))>0.d0) then
        qb_x_t_class(i) = (Uvec_a_x_t(i)/Uvec_a_t(i))*qb_t_class(i)
        qb_y_t_class(i) = (Uvec_a_y_t(i)/Uvec_a_t(i))*qb_t_class(i)
     else
        qb_x_t_class(i) = 0.d0
        qb_y_t_class(i) = 0.d0
     endif
  enddo

! mean bed-load transport (m2.s-1)
  qb_class(1) = sum(qb_x_t_class(:))/ech/rhosed
  qb_class(2) = sum(qb_y_t_class(:))/ech/rhosed

  if(qb_class(1)/=qb_class(1)) call parallel_abort('qb_class(1) = NaN')
  if(qb_class(2)/=qb_class(2)) call parallel_abort('qb_class(2) = NaN')

! vertical dicretization
  z(:) = 0.d0
  do i=1,ndz-1
     z(i) = z_a*(htot/z_a)**(dfloat(i-1)/dfloat(ndz-1))
  enddo
  z(ndz) = htot

! vertical profile of current velocity
  Uc_z(:) = 0.d0
  do i=1,ndz
     if (z(i)>=dm) then
        Uc_z(i) = Uc*dlog(30.d0*z(i)/k_a)/(-1.d0+dlog(30.d0*htot/k_a))
     else
        Uc_z(i) = Uc_dm*dlog(30.d0*z(i)/ks_c)/(dlog(30.d0*dm/ks_c))
     endif
     if(Uc_z(i)<0.d0) Uc_z(i) = 0.d0
  enddo

! current-related friction coefficient
  if (htot>0.d0.and.ks_c>0.d0) then
    f_c = 0.24d0*(dlog10(12.d0*htot/ks_c))**(-2.d0)
  else
    f_c = f_grain_c
  endif

  if(f_c/=f_c) call parallel_abort('f_c = NaN')

! wave-related friction coefficient
  if (ks_w_r>0.d0) then
    tmp = 5.2d0*(Ad/ks_w_r)
    if (tmp>0.d0) then
      f_w = dexp(-6.d0+tmp**(-0.19d0))
      f_w = min(f_w,0.3d0)
    else
      f_w = f_grain_w
    endif
  else
    f_w = f_grain_w
  endif  

  if(f_w/=f_w) call parallel_abort('f_w = NaN')

! current and wave-related efficiency factors
  if (f_c>0.d0) then
    mu_c = f_grain_c/f_c
  else
    call parallel_abort('tr04_tot_multi_size: f_c = 0')
  endif

  mu_w_class = 0.7d0/D_star_class
  if (D_star_class>=5.d0) then
    mu_w_class = max(mu_w_class,0.14d0)
  elseif (D_star_class<=2.d0) then
    mu_w_class = min(mu_w_class,0.35d0)
  endif

! wave-current interaction coefficient
  alpha_cw = (dlog(30.d0*dm/k_a)/dlog(30.d0*dm/ks_c))**2.d0*&
   &((-1.d0+dlog(30.d0*htot/ks_c))/(-1.d0+dlog(30.d0*htot/k_a)))**2.d0
  alpha_cw = min(alpha_cw,1.d0)

  if(alpha_cw/=alpha_cw) call parallel_abort('alpha_cw = NaN')

! current and wave-related bed-shear stress
  tau_b_c = 0.125d0*rho0*f_c*Uc**2.d0
  tau_b_w = 0.25d0*rho0*f_w*U_d_w_r**2.d0

! combined current and wave-related bed shear stress
  tau_b_cw = alpha_cw*tau_b_c + tau_b_w

! combined current and wave-related friction velocity
  u_star_cw = dsqrt(tau_b_cw/rho0)

! grain-related bed-shear stress for combined current and waves
  ! method 1
  tau_grain_b_cw_class = alpha_cw*mu_c*tau_b_c + mu_w_class*tau_b_w
  ! method 2
!  tau_grain_b_cw_class_2 = sum(tau_grain_b_cw_t(:))/ech

  if(tau_grain_b_cw_class/=tau_grain_b_cw_class) call parallel_abort('tau_grain_b_cw_class = NaN')

! dimensionless bed-shear stress for combined current and waves
  tmp = tau_grain_b_cw_class-tau_cr_1*d_class/d50*xi_class
  if (tmp>0.d0) then
    tau_adim_cw_class = lambda_class*tmp/(tau_cr*d_class/d50)
  else
    tau_adim_cw_class = 0.d0
  endif

  if(tau_adim_cw_class/=tau_adim_cw_class) call parallel_abort('tau_adim_cw_class = NaN')

! sediment mixing coef for current
  if (ks_c>0.d0) then
    Chezy = 18.d0*dlog10(12.d0*htot/ks_c)
  else
    call parallel_abort('pb Chezy coef: ks_c <= 0')
  endif
  if (Chezy>0.d0) then
    u_star_c = dsqrt(grav)/Chezy*Uc
  else
    call parallel_abort('Chezy coef <= 0')
  endif
!  u_star_c = dsqrt(tau_b_c/rho0) !TEST
  if (u_star_c>0.d0) then
    beta_c_class = 1.d0+2.d0*(Ws_class/u_star_c)**2.d0
    beta_c_class = min(beta_c_class,1.5d0)
  else
    beta_c_class = 1.d0
  endif

  epsi_s_c_class(:) = 0.d0
  do i=1,ndz
     if (z(i)<0.5d0*htot) then
       epsi_s_c_class(i) = kappa*beta_c_class*u_star_c*z(i)*(1.d0-z(i)/htot)
     elseif (z(i)>=0.5d0*htot) then
       epsi_s_c_class(i) = 0.25d0*kappa*beta_c_class*u_star_c*htot
     endif
     if(epsi_s_c_class(i)/=epsi_s_c_class(i)) call parallel_abort('epsi_s_c_class(i) = NaN')
  enddo

! sediment mixing coef for waves
  epsi_s_w_class(:) = 0.d0
  if (hs>hs_min.and.tp>tp_min) then
    if (hs/htot<0.4d0) then
      gama_br = 1.d0
    else
      gama_br = 1.d0+dsqrt(hs/htot-0.4d0)
    endif

    ds = 2.d0*gama_br*dw
    if(ds<0.05d0) ds=0.05d0
    if(ds>0.2d0) ds=0.2d0

    u_star_w = dsqrt(tau_b_w/rho0)
!    if (u_star_w>0.d0) then
!      beta_w_class = 1.d0+2.d0*(Ws_class/u_star_w)**2.d0
!      beta_w_class = min(beta_w_class,1.5d0)
!    else
!      call parallel_abort('u_star_w <= 0')
!    endif
    if (u_star_w>0.d0) then
      beta_w_class = 1.d0+2.d0*(Ws_class/u_star_w)**2.d0
      beta_w_class = min(beta_w_class,1.5d0)
    elseif (u_star_w==0.d0) then
      beta_w_class = 1.5d0
    else
      call parallel_abort('u_star_w < 0')
    endif

    epsi_s_w_bed_class = 0.018d0*beta_w_class*gama_br*ds*U_d_w_r
    epsi_w_max = 0.035d0*gama_br*htot*hs/tp

    do i=1,ndz
       if (z(i)<=ds) then
         epsi_s_w_class(i) = epsi_s_w_bed_class
       elseif (z(i)>ds.and.z(i)<0.5d0*htot) then
         epsi_s_w_class(i) = epsi_s_w_bed_class + (epsi_w_max-epsi_s_w_bed_class) &
               &*((z(i)-ds)/(0.5d0*htot-ds))
       elseif (z(i)>=0.5d0*htot) then
         epsi_s_w_class(i) = epsi_w_max
       endif
       if(epsi_s_w_class(i)/=epsi_s_w_class(i)) call parallel_abort('epsi_s_w_class(i) = NaN')
    enddo
  endif

! sediment mixing coef for current and waves
  epsi_s_cw_class(:) = 0.d0
  do i=1,ndz
     epsi_s_cw_class(i) = dsqrt(epsi_s_c_class(i)**2.d0+epsi_s_w_class(i)**2.d0)
  enddo

! near-bed sediment concentration (dimensionless)
  if(z_a<=0.d0) call parallel_abort('z_a <= 0')
  c_a_class = 0.015d0*(d_class/z_a)*tau_adim_cw_class**1.5d0/D_star_class**0.3d0
  c_a_class = min(c_a_class,0.05d0)

  if(c_a_class/=c_a_class) call parallel_abort('c_a_class = NaN')

! sediment concentration vertical profile (dimensionless)
  c_z_class(:) = 0.d0
  if (Uc>0.d0) then
    c_z_class(1) = c_a_class
    do i=2,ndz
       beta_d_class = 1.d0+(c_z_class(i-1)/0.65d0)**0.8d0 - 2.d0*(c_z_class(i-1)/0.65d0)**0.4d0
       if (beta_d_class>0.d0.and.epsi_s_cw_class(i-1)>0.d0) then
         dc_dz_class = -(1.d0-c_z_class(i-1))**5.d0*c_z_class(i-1)*Ws_class &
                       &/(beta_d_class*epsi_s_cw_class(i-1))
       else
         write(12,*)'beta_d_class',beta_d_class,'epsi_s_cw_class(i-1)',epsi_s_cw_class(i-1)
         call parallel_abort('pb dc_dz_class')
       endif

       c_z_class(i) = (z(i)-z(i-1))*dc_dz_class + c_z_class(i-1)

       if(c_z_class(i)<0.d0) c_z_class(i) = 0.d0

       if (c_z_class(i)/=c_z_class(i)) then
         write(12,*)'c_z_class(i)',c_z_class(i)
         call parallel_abort('c_z_class(i)')
       endif

    enddo
  endif

#ifdef INCLUDE_TIMING
      timer_ns(13) = timer_ns(13)+mpi_wtime()-time_tmp !end of timer 
      time_tmp = mpi_wtime()                         !start of timer
#endif

! current-related suspended load transport (kg.m-1.s-1)
! vertical integration
  qs_c_class = z_a*Uc_a*c_a_class
  do i=2,ndz
     qs_c_class = qs_c_class + (z(i)-z(i-1))  &
          &*(Uc_z(i)*c_z_class(i)+Uc_z(i-1)*c_z_class(i-1))/2.d0
  enddo
  qs_c_class = qs_c_class*rhosed

  if(qs_c_class<0.d0) qs_c_class = 0.d0
  if(qs_c_class/=qs_c_class) call parallel_abort('qs_c_class = NaN')

! streaming velocity at edge of wave boundary layer
  if (hs>hs_min.and.tp>tp_min.and.wlpeak>wlpeak_min) then

    if(ks_w_r==0.d0) call parallel_abort('tr04_tot_multi_size: ks_w_r = 0')
    tmp = Ad/ks_w_r

    if(cphase==0.d0) call parallel_abort('tr04_tot_multi_size: cphase = 0')

    if (tmp>1.d0.and.tmp<100.d0) then
      u_d = (-1.d0+0.875d0*dlog10(tmp))*(Uw**2.d0/cphase)
    elseif (tmp>=100.d0) then
      u_d = 0.75d0*(Uw**2.d0/cphase)
    elseif (tmp<=1.d0) then
      u_d = -(Uw**2.d0/cphase)
    endif
  else
    u_d = 0.d0
  endif

! wave-related suspended load transport (kg.m-1.s-1)
  if (hs>hs_min.and.tp>tp_min.and.wlpeak>wlpeak_min) then
    tmp = U_d_for**3.d0 + U_d_back**3.d0
    if (dabs(tmp)>0.d0) then
      Uasym = (U_d_for**4.d0 - U_d_back**4.d0)/tmp
    else
      Uasym = 0.d0
    endif
    F = 0.1d0*(Uasym + u_d)

    !vertical integration
    qs_w_class = z_a*c_a_class
    i=2
    do while(z(i)<=0.5d0.and.i<=ndz)
       qs_w_class = qs_w_class + (z(i)-z(i-1))*(c_z_class(i)+c_z_class(i-1))/2.d0
       i=i+1
       if (z(i-1)>0.5d0) then
         write(12,*)'z(i-1)',z(i-1),'z(i)',z(i)
         call parallel_abort('pb while loop')
       endif
    enddo
    qs_w_class = qs_w_class*rhosed*F

    if(qs_w_class<0.d0) qs_w_class = 0.d0

  else
    qs_w_class = 0.d0
  endif

  if(qs_w_class/=qs_w_class) call parallel_abort('qs_w_class = NaN')

! suspended-load transport (m2.s-1)
! current
  qs_c_class_x = qs_c_class*dcos(udir)/rhosed
  qs_c_class_y = qs_c_class*dsin(udir)/rhosed
! waves
  qs_w_class_x = qs_w_class*dcos(wdir)/rhosed
  qs_w_class_y = qs_w_class*dsin(wdir)/rhosed
! current and waves
  qs_class(1) = qs_c_class_x + qs_w_class_x
  qs_class(2) = qs_c_class_y + qs_w_class_y

  if(qs_class(1)/=qs_class(1)) call parallel_abort('qs_class(1) = NaN')
  if(qs_class(2)/=qs_class(2)) call parallel_abort('qs_class(2) = NaN')

  END SUBROUTINE tr04_tot

!--------------------------------------------------------------------
  SUBROUTINE lesser04_slope(Cd,d50,u,v,dpdxy,tau_cr,alpha_s,alpha_n)
!--------------------------------------------------------------------
! This subroutine computes the coefficients of Lesser et al. (2004) 
! to take into account the slope effect on bed-load transport
!
! Author: guillaume dodet (gdodet01@univ-lr.fr, gdodet@lnec.pt)
! Date: 11/03/2013
!
! History:
!--------------------------------------------------------------------

  USE schism_glbl, ONLY : pi,rho0,rkind
  USE schism_msgp, ONLY : parallel_abort
  USE sed2d_mod, ONLY : islope

  IMPLICIT NONE

!- Arguments --------------------------------------------------------
  REAL(rkind), INTENT(IN) :: Cd,d50,tau_cr,u,v
  REAL(rkind), DIMENSION(2), INTENT(IN) :: dpdxy
  REAL(rkind), INTENT(OUT) :: alpha_n,alpha_s
!- Local variables --------------------------------------------------
  REAL(rkind) :: cff,lslope,sed_angle,tau,tslope,udir,unorm
!--------------------------------------------------------------------

  unorm = SQRT(u*u+v*v)
  udir = ATAN2(v,u)
  lslope = dpdxy(1)*COS(udir)+dpdxy(2)*SIN(udir) !longitudinal slope
  tslope = dpdxy(2)*COS(udir)-dpdxy(1)*SIN(udir) !transverse slope

  IF(islope == 2) THEN
    sed_angle = TAN(30.d0*pi/180.d0)
    cff = MIN(ABS(lslope),0.9d0*sed_angle)*SIGN(1.d0,lslope)
    alpha_s = 1.d0+1.d0*((sed_angle/(COS(ATAN(cff))*(sed_angle-cff)))-1.d0)
    tau = rho0*Cd*unorm**2.d0
    IF(tau == 0.d0) THEN
       alpha_n = 0.d0
    ELSE
       alpha_n = 1.5d0*SQRT(tau_cr/tau)*tslope
    ENDIF
  ELSE
    alpha_n = 0.d0
    alpha_s = 1.d0
  ENDIF

  END SUBROUTINE lesser04_slope

!--------------------------------------------------------------------
  SUBROUTINE wave_asymmetry_Elfrink(H,wave_per,w_dir,depth,dhxi,dhyi,&
                     ech,Uorbi,Ucrest,Utrough,T_crest,T_trough)
!--------------------------------------------------------------------
! This subroutine computes wave asymmetry based on Elfrink et al. 
! (2006, Coastal Engineering)
!
! NB: offshore wave height and wave period are taken into account to
! compute irribaren number and offshore wave length
!
! Author: thomas guerin (thomas.guerin@univ-lr.fr)    
! Date: 26/04/2013
!--------------------------------------------------------------------
  
  USE schism_msgp, ONLY : parallel_abort
  IMPLICIT NONE

!- Arguments --------------------------------------------------------  
  REAL(8), INTENT(IN) :: H,wave_per,w_dir,depth,dhxi,dhyi
  INTEGER, INTENT(IN) :: ech
  REAL(8), INTENT(OUT) :: Ucrest,Utrough,T_crest,T_trough
  REAL(8), DIMENSION(ech-1), INTENT(OUT) :: Uorbi
!- Constants --------------------------------------------------------  
  REAL(8), PARAMETER :: g = 9.80665d0
  REAL(8), PARAMETER :: pi = 3.141592653589793d0 !DACOS(-1.d0)
  REAL(8), PARAMETER :: a1 = 0.38989d0
  REAL(8), PARAMETER :: a2 = -0.0145d0
  REAL(8), PARAMETER :: a3 = -0.0005d0
  REAL(8), PARAMETER :: a4 = 0.5028d0
  REAL(8), PARAMETER :: a5 = 0.9209d0
  REAL(8), PARAMETER :: b1 = 0.5366d0
  REAL(8), PARAMETER :: b2 = 1.16d0
  REAL(8), PARAMETER :: b3 = -0.2615d0
  REAL(8), PARAMETER :: b4 = 0.0958d0
  REAL(8), PARAMETER :: b5 = -0.5623d0
!- Local variables --------------------------------------------------
  REAL(8) :: a1bis,C1,C2,C3,C4,C5,D1,D2,D3,D4,D5,E1,E2,E3,F1,F2,F3,  &
             F4,G1,G2,G3,G4,G5,G6,G7,G8,Hadim,kh,L0,Ladim,P1,P2,P3,  &
             P4,P5,psi,Slope,T0,T1,T2,tv,U0,U1,U2,Uairy,uorb,Ur,     &
             Ustar,Zeta,tmp
  INTEGER :: t
!--------------------------------------------------------------------

!- Compute offshore wave length, local adimensional wave length, and
!- orbital velocity according to linear wave theory
  if(depth<=0) call parallel_abort('SED_TRANS: (1)')
  psi = 4.d0*pi**2.d0*depth/(g*wave_per**2.d0) !>0
  !kh>0
  IF (psi .LE. 1.d0) THEN
      kh = DSQRT(psi)*(1.d0 + 0.2d0*psi)
  ELSE
      kh = psi*(1.d0 + 0.2d0*DEXP(2.d0-2.d0*psi))
  ENDIF
  Ladim = 2.d0*pi/kh !>0
  Hadim = H/depth
  uorb = pi*Hadim*depth/wave_per/DSINH(kh)

!- Compute bed slope in the direction of wave propagation, irribaren 
!- number, and Ursell number
  Slope = -DSQRT(dhxi**2.d0+dhyi**2.d0)*DCOS(w_dir+DATAN2(dhyi,dhxi))
  if(Hadim<0.or.Ladim<=0) call parallel_abort('SED_TRANS: (2)')
  Zeta = DTAN(Slope)/DSQRT(Hadim/Ladim)
  Ur = Hadim*Ladim**2.d0

!- Compute the normalized maximal orbital velocity
  C1 = Ladim-10.d0
  C2 = DABS(C1-(Hadim-DABS(Zeta)))
  C3 = Zeta*(1.d0-C1)
  C4 = DTANH(DABS(C3-C2)/Ur)
  tmp=DABS(Zeta)+DTANH(C4)
  if(Hadim<=0.or.tmp<0) call parallel_abort('SED_TRANS: (3)')
  !Ur>0
  C5 = DSQRT(tmp) !DABS(Zeta)+DTANH(C4))
  P1 = DSQRT(Hadim)-C5*Hadim
  U1 = b1*P1+a1

!- Compute the velocity asymmetry parameter
  D1 = 3.d0*Zeta+2.d0*Ladim/Ur
  D2 = DSQRT(Ladim)-DTANH(DABS(D1))
  D3 = (2.d0*Zeta+DSQRT(Ladim/Ur))**2.d0
  D4 = Ur+Ladim/D3/Ur
!  if(D2/D4<0) call parallel_abort('SED_TRANS: (4)')
  IF (D2/D4<=0) THEN
    D5 = 0.d0
  ELSE
    D5 = DSQRT(D2/D4)
  ENDIF
  P2 = 1.2001d0*D5+0.4758d0
  U2 = b2*P2+a2

!- Compute the normalized phase of wave crest
  E1 = Hadim*Ladim*Zeta
  E2 = E1*(-9.8496d0*Zeta*Hadim)**2.d0
  E3 = DTANH(E2)+DTANH(E1)+Ladim-1.d0
  P3 = DTANH(-9.3852d0/E3)
  T1 = b3*P3+a3

!- Compute the normalized phase of zero down crossing
  F1 = 0.0113d0*Zeta*Ladim**2.d0
  F2 = 0.00035667d0*Zeta*Ladim**4.d0
  F3 = 0.1206d0*Ladim*DTANH(DTANH(Zeta))
  F4 = Hadim*DTANH(F2)/DTANH(F3)
  P4 = Hadim*DTANH(0.02899d0*Ladim*F1)-DTANH(F4)
  T0 = b4*P4+a4

!- Compute the normalized phase of wave trough
  G1 = Zeta+0.9206d0
  G2 = Ladim-DSQRT(Ur)+DSQRT(2.5185d0/Ladim)-4.6505d0
  G3 = DSQRT(DABS(G2/Hadim))
  G4 = DABS(Zeta+Ladim)-4.4995d0+Zeta
  G5 = DABS(G4+DABS(Zeta)-5.3981d0)
  G6 = DABS(Ladim+DSQRT(3.0176d0/Hadim)-5.2868d0+Hadim)
  G7 = DABS(Zeta+0.1950d0*(G6+Zeta))
  G8 = DABS(Zeta)+Ladim
  P5 = 4.1958d0/(G1+G3+G5+G7+G8)
  T2 = b5*P5+a5

!- Compute orbital velocity parameters and correct U0
  Uairy = uorb
  Ustar = 2.d0*U2*Uairy
  Ucrest = U1*Ustar
  Utrough = Ustar-Ucrest
  U0 = (Ucrest*T0-Utrough*(1.d0-T0))/(T0-T1)
  IF (U0 .GT. (0.25d0*Ucrest)) THEN
      U0 = 0.25d0*Ucrest
  ENDIF

  tmp=T0*(U0-Ucrest-Utrough)
  if(tmp==0) call parallel_abort('SED_TRANS: (5)')
  a1bis = (-Utrough+T1*U0)/tmp !(T0*(U0-Ucrest-Utrough))
  IF (a1bis .LT. 0.99d0) THEN
    T0 = 0.99d0*T0
    if(T0==1) call parallel_abort('SED_TRANS: (6)')
    Utrough = (-(T0-T1)*U0+T0*Ucrest)/(1.d0-T0)
  ELSE
    T0 = a1bis*T0
  ENDIF

!- Compute wave half-periods
  T_crest = T0*wave_per
  T_trough = wave_per-T_crest

!- Rebuilt of orbital velocity over one wave period
  if(ech==1) call parallel_abort('SED_TRANS: (7)')
  DO t=1,ech-1
    tv = (t-1.d0)/(ech-1.d0)
    IF ((tv .GE. 0.d0) .AND. (tv .LE. T1)) THEN
      if(T1==0) call parallel_abort('SED_TRANS: (8)')
      Uorbi(t) = Ucrest*DSIN(pi/2.d0*tv/T1)
    ELSEIF ((tv .GT. T1) .AND. (tv .LE. T0)) THEN
      if(T1==T0) call parallel_abort('SED_TRANS: (9)')
      Uorbi(t) = Ucrest*DCOS(pi/2.d0*(tv-T1)/(T0-T1))                &
                 - U0*DSIN(pi*(tv-T1)/(T0-T1))
    ELSEIF ((tv .GT. T0) .AND. (tv .LE. T2)) THEN
      if(T2==T0) call parallel_abort('SED_TRANS: (10)')
      Uorbi(t) = -Utrough*DSIN(pi/2.d0*(tv-T0)/(T2-T0))
    ELSEIF ((tv .GT. T2) .AND. (tv .LE. 1.d0)) THEN
      if(T2==1) call parallel_abort('SED_TRANS: (11)')
      Uorbi(t) = -Utrough*DCOS(pi/2.d0*(tv-T2)/(1.d0-T2))
    ENDIF
  ENDDO

  END SUBROUTINE wave_asymmetry_Elfrink
!--------------------------------------------------------------------
END MODULE sed2d_transport
