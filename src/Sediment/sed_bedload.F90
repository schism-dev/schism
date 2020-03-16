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

!=====================================================================
!=====================================================================
! MORSELFE BEDLOAD SUBROUTINES
!
! subroutine sed_bedload_sd
! subroutine sed_bedload_vr
! subroutine sed_bedload_mpm
! subroutine sed_bedload_wl14
! subroutine bedchange_bedload
! subroutine wave_asymmetry_Elfrink
! subroutine bedload_flux_limiter
! subroutine weno_scheme
! subroutine acceleration_bedload_hoe2003
! subroutine acceleration_bedload_dub2015
! subroutine bedslope_effects_lesser2004
! 
!=====================================================================
!=====================================================================
     SUBROUTINE sed_bedload_sd(ised,inea)
!--------------------------------------------------------------------!
! This routine computes bedload according to Soulsby and Damgaard (2005)         !
!                                                                    !
! Author: Knut Kramer                                                !
! Date: 27/11/2012                                                   !
!                                                                    !
! History: 2012/12 - F.Ganthy : form homogenisation of sediments     !
!          routines                                                  !
!          2013/03 - F.Ganthy : implement wave effects on bedload    !
!          2017/06 - A. de Bakker: couple to schism                  !
!          2020/02 - B.Mengual: > use of ustress/vstress, wdir       !
!                       > wave asymmetry effect: rewritten and       !
!                         introduction of a limiter w_asym_max       !
!                       > bed slope effects: now under option with   !
!                         an external call to a new subroutine       !
!                                                                    !
!--------------------------------------------------------------------!

      USE sed_mod
      USE schism_glbl, ONLY: rkind,errmsg,dpe,nea,dt,kbe,i34,elnode,eta2,pi, & 
                                uu2,vv2,ielg
      USE schism_msgp, ONLY: myrank,parallel_abort

      IMPLICIT NONE
      SAVE

!- Local variables --------------------------------------------------!

      INTEGER, INTENT(IN)  :: ised,inea
      INTEGER :: j,inode
      ! arbitrary coefficients
      REAL(rkind) :: cff, cff1, cff2, cff3, cff4
      ! Hydrodynamic characteristiques
      REAL(rkind) :: htot, udir, phi, wave_h, wave_per, wave_dir
      ! Param for theta_cr
      REAL(rkind),PARAMETER :: nu=1.36d-6 ! Cinematic viscosity
      ! Shields parameter etc.
      REAL(rkind) :: theta_m, theta_w, smg, osmgd, smgdr  
      ! Wave asymmetry    
      REAL(rkind) :: w_asym 
      ! Initiation of motion
      REAL(rkind) :: ratio, dstar, theta_cr, theta_max1, theta_max2, &
                    & theta_max
      ! Transport nondimensional
      REAL(rkind) :: phi_x1, phi_x2, phi_x, phi_y
      ! Transport dimensional
      REAL(rkind) :: bedld_x, bedld_y
      ! Element bottom vertical indices +1 
      REAL(rkind) :: kb1


!- Start Statement --------------------------------------------------!

      cff1  = 0.d0
      cff2  = 0.d0
      cff3  = 0.d0
      cff4  = 0.d0
      bedld_x = 0.d0
      bedld_y = 0.d0

!--------------------------------------------------------------------
!- Water depth and wave conditions                                   
!--------------------------------------------------------------------
      htot = dpe(inea)+sum(eta2(elnode(1:i34(inea),inea)))/i34(inea)
      wave_per = tp(inea)
      wave_h = hs(inea) 
      wave_dir = wdir(inea)

!--------------------------------------------------------------------
!- Bedload characteristics                                           
!--------------------------------------------------------------------
      smg = (Srho(ised)/rhom-1.0d0)*g
!      smgd = smg * Sd50(ised)
      osmgd = 1.0d0/smgd
      smgdr = sqrt(smgd)*Sd50(ised)
 
!--------------------------------------------------------------------
!- Current amplitude and direction                                   
!-------------------------------------------------------------------- 

      ! Reference code
      !kb1 = kbe(inea)+1
      !cff1 = sum(uu2(kb1,elnode(1:i34(inea),inea)))/i34(inea)
      !cff2 = sum(vv2(kb1,elnode(1:i34(inea),inea)))/i34(inea)
      !udir = DATAN2(cff2,cff1)

      ! BM: see ustress, vstress estimates in sediment.F90
      udir = ATAN2(vstress(inea),ustress(inea))

!--------------------------------------------------------------------
!- Angle between current and wave directions                         
!--------------------------------------------------------------------
#ifdef USE_WWM
      phi = wave_dir - udir
#else
      phi = 0.d0
#endif

!--------------------------------------------------------------------
!- Wave asymmetry treatment : compute wave orbital velocity          
!  (horizontal component) from Elfrink (2006)                        
!--------------------------------------------------------------------

      !BM: useless CALL, already done in sediment.F90
!      CALL sed_wavecurrent_stress() !!! compute bottom shear stress for wave current

#ifdef USE_WWM
      IF ( U_crest(inea)+U_trough(inea) > 0) THEN
        ! BM: limiting w_asym to w_asym_max according to Soulsby and
        ! Damgaard (2005), who recommend w_asym_max=0.2
        ! Symmetric wave U_crest=U_trough -> w_asym=0
        w_asym = MAX(MIN(U_crest(inea)/((U_crest(inea)+U_trough(inea))/2.0d0)-1.0d0,&
                        &  w_asym_max),&
                        &  0.d0)
      ELSE 
        w_asym = 0.d0
      ENDIF

#else
      w_asym = 0.d0
#endif


!      WRITE(16,*),'wave_h',wave_h,'htot',htot,'w_asym',w_asym
      
!--------------------------------------------------------------------
!- Compute non dimensional stresses                                  
!--------------------------------------------------------------------
#ifdef USE_WWM
      theta_w = tau_w(inea)*osmgd ! wave component
#else   
      theta_w = 0.d0
#endif
         
      theta_m = tau_wc(inea)*osmgd ! component due to the waves+current 
      cff1 = theta_w*(1.0d0+w_asym)
      cff2 = theta_w*(1.0d0-w_asym)

!-------------------------------------------------------------------- 
!- Equations 34 and 35, Soulsby and Damgaard 2005                    
!--------------------------------------------------------------------
      theta_max1 = SQRT((theta_m+cff1*cos(phi))**2+(cff1*sin(phi))**2)
      theta_max2 = SQRT((theta_m+cff2*cos(phi+pi))**2+(cff2*sin(phi+pi))**2)
      theta_max = max(theta_max1,theta_max2)

!--------------------------------------------------------------------
!- Motion initiation factor                                          
!--------------------------------------------------------------------
!      theta_cr = 0.04 ! based on S&D page 678

      ratio    = Srho(ised)/rhom
      dstar    = Sd50(ised) * (g*(ratio-1.d0)/nu**2.d0)**(1.d0/3d0)
      IF(1.d0+1.2d0*dstar==0) THEN
         WRITE(errmsg,*)'SED3D sed_bedload: dev. by 0'
         CALL parallel_abort(errmsg)
      ENDIF 
      theta_cr = 0.3d0/(1.d0+1.2d0*dstar) + 0.055d0 * (1.d0-EXP(-0.02d0*dstar))
      cff3 = 0.5d0*(1.0d0+sign(1.0d0,theta_max/theta_cr-1.0d0))

!-------------------------------------------------------------------------
!- Calculate bedload in direction of currents and perpendicular direction 
!- Equation 31 and 32 SD05                                                
!-------------------------------------------------------------------------
      phi_x1 = 12.0d0 * sqrt(theta_m)* (theta_m - theta_cr)
      phi_x2 = 12.0d0 * (0.9534d0 + 0.1907d0*cos(2.d0*phi))*sqrt(theta_w)*theta_m + 12.0d0*(0.229d0*w_asym*theta_w**1.5d0*cos(phi))

      IF (ABS(phi_x2).gt.phi_x1) THEN
         phi_x = phi_x2
      ELSE
         phi_x = phi_x1
      ENDIF
      bedld_x = phi_x*smgdr*cff3
      
!--------------------------------------------------------------------
!- Equation 33 SD05                                                 
!--------------------------------------------------------------------
      cff4 = theta_w**1.5d0+1.5d0*(theta_m**1.5d0)
      phi_y=12.0d0*0.1907d0*theta_w**2.0d0*(theta_m*sin(2.0d0*phi) + 1.2d0*w_asym*theta_w*sin(phi))/cff4
      IF (phi_y /= phi_y) THEN
        bedld_y = 0.d0
      ELSE
        bedld_y = phi_y*smgdr*cff3
      ENDIF

!--------------------------------------------------------------------
!- Fluxes                                                            
!--------------------------------------------------------------------
!      FX_r(inea) = sqrt(bedld_x**2.d0+bedld_y**2.d0)*dcos(udir)*dt
!      FY_r(inea) = sqrt(bedld_x**2.d0+bedld_y**2.d0)*dsin(udir)*dt
      FX_r(inea) = (bedld_x*cos(udir)-bedld_y*sin(udir))*dt
      FY_r(inea) = (bedld_y*cos(udir)+bedld_x*sin(udir))*dt


!--------------------------------------------------------------------
!- Bed slope effects (Lesser et al. 2004)                                                            
!--------------------------------------------------------------------
      IF (slope_formulation==4) THEN
        CALL bedslope_effects_lesser2004(inea,ised,FX_r(inea),FY_r(inea))
      END IF

!--------------------------------------------------------------------
!- Sanity check                                                      
!--------------------------------------------------------------------
      IF(FX_r(inea)/=FX_r(inea)) THEN
        WRITE(errmsg,*)'FX_r0 is NaN',myrank,inea,htot,wave_h,FX_r(inea),bedld_x,  &
        &              phi_x,phi
        CALL parallel_abort(errmsg)
      ENDIF

      IF(FY_r(inea)/=FY_r(inea)) THEN
        WRITE(errmsg,*)'FY_r0 is NaN',myrank,inea,htot,wave_h,FY_r(inea),bedld_y,  &
        &              phi_y,phi
        CALL parallel_abort(errmsg)
      endif

!---------------------------------------------------------------------
      END SUBROUTINE sed_bedload_sd

!=====================================================================
!=====================================================================

      SUBROUTINE sed_bedload_vr(ised,inea,dave)
!--------------------------------------------------------------------!
! This routine computes bedload according to van Rijn (2007)         !
!                                                                    !
! Author: Knut Kramer                                                !
! Date: 27/11/2012                                                   !
!                                                                    !
! History: 2012/12 - F.Ganthy : form homogenisation of sediments     !
!          routines                                                  !
!          2013/03 - F.Ganthy : implement wave effects on bedload    !
!          2020/02 - B.Mengual: bed slope effects: now under option  !
!                        with an external call to a new subroutine   !
!                                                                    !
!--------------------------------------------------------------------!

      USE sed_mod
      USE schism_glbl, ONLY: rkind,errmsg,dpe,nea,dt,i34,elnode,eta2
      USE schism_msgp, ONLY: myrank,parallel_abort

      IMPLICIT NONE
      SAVE

!- Local variables --------------------------------------------------!

      INTEGER, INTENT(IN)  :: ised,inea
      ! depth averaged elementwise hvel
      REAL(rkind), INTENT(IN) :: dave(nea)

      ! arbitrary coefficients
      REAL(rkind) :: cff, cff1, cff2, smg
      REAL(rkind) :: a_slopex, a_slopey
      REAL(rkind) :: bedld, me
      ! Total water depth
      REAL(rkind) :: htot
      ! Critical velocities
      REAL(rkind) :: ucrc ! Current alone
      REAL(rkind) :: ucrw ! Wave alone
      REAL(rkind) :: ucrt ! Current-wave
      REAL(rkind) :: ueff ! effective velocity
      REAL(rkind) :: beta
      REAL(rkind), PARAMETER :: gama = 0.4d0 !0.8 for regular waves

!- Start Statement --------------------------------------------------!

      ucrc  = 0.d0
      ucrw  = 0.d0
      ucrt  = 0.d0
      ueff  = 0.d0
      cff   = 0.d0
      cff1  = 0.d0
      cff2  = 0.d0
      me    = 0.d0
      bedld = 0.d0

      smg = (Srho(ised)/rhom-1.0d0)*g

!---------------------------------------------------------------------
! - Critical velocity for currents based on Shields (initiation of 
! motion), (van Rijn 2007a)
! - For current:
! Ucrc = 0.19*d50**0.1*log10(4*h/12*d90) if 5e-5 < d50 < 5e-4
! Ucrc = 8.5*d50**0.6*log10(4*h/12*d90) if 5e-4<= d50 <2e-3
! We do not have d90 in current implementation
! - For waves:
! Ucrw = 0.24*((s-1)*g)**0.66*d50**0.33*Tp**0.33 if 5e-5 < d50 < 5e-4
! Ucrw = 0.95*((s-1)*g)**0.57*d50**0.43*Tp**0.14 if 5e-4<= d50 <2e-3
!---------------------------------------------------------------------

      !dpe is the min of nodes
      htot = dpe(inea)+sum(eta2(elnode(1:i34(inea),inea)))/i34(inea) 

      IF (Sd50(ised)>5.0d-5.AND.Sd50(ised)<5.0d-4) THEN

        ucrc = 0.19d0*Sd50(ised)**0.1d0*LOG10(4.0d0*htot/Sd50(ised))
        ucrw = 0.24d0*smg**0.66d0*Sd50(ised)**0.33d0*tp(inea)**0.33d0

      ELSEIF (Sd50(ised)>=5.0d-4.AND.Sd50(ised)<2.0d-3)THEN

        ucrc = 8.5d0*Sd50(ised)**0.6d0*LOG10(4.0d0*htot/Sd50(ised))
        ucrw = 0.95d0*smg**0.57d0*Sd50(ised)**0.43d0*tp(inea)**0.14d0

      ELSE

        WRITE(errmsg,*)'Sediment diameter out of range:',ised,       &
        &              Sd50(ised)
        CALL parallel_abort(errmsg)

      ENDIF

!---------------------------------------------------------------------
! - Effective velocity and wave-current (total) critical velocity
! Currently Uorb RMS is used instead of Uorb peak, to prevent
! instabilities due to linear wave theory applied in breaking zone
! to compute Uorb peak
!---------------------------------------------------------------------
      ueff = dave(inea) + gama*uorb(inea)
      beta = dave(inea)/(dave(inea)+uorb(inea)) !denom /=0
      ucrt = beta*ucrc + (1.0d0-beta)*ucrw

!---------------------------------------------------------------------
! - Mobility parameter (van Rijn 2007a)
!---------------------------------------------------------------------

      IF (ueff.GT.ucrt) THEN
        me = (ueff-ucrt)/SQRT(smgd) !smgd checked in sediment.F90
      ENDIF

!---------------------------------------------------------------------
! - Simplified bed load transport formula for steady flow 
! (van Rijn, 2007)
! rho_s (van Rijn) missing at this point
! instead of [kg/m/s] (van Rijn) here [m2/s]
!
!jl. Looks to me like the units are [m2/s] 
!    In ROMS the bedld is integrated in time(s) and space(1-d, size of
!    RHO-grid spacing to end up with [kg]. Here it looks like rho is 
!    not used yet so the JCG can be used to calculate the change 
!    layer thickness due bedload flux.
!
!    rho is eventually considered at 762
! 
!---------------------------------------------------------------------
      IF (ueff.GT.ucrt) THEN
        bedld = 0.015d0*dave(inea)*htot*(Sd50(ised)/htot)**1.2d0*    &
        &       me**1.5d0 ![m^2/s]
      ENDIF

!---------------------------------------------------------------------
! - Partition bedld into x and y directions, at the center of each
! element (rho points), and integrate in time.
! FX_r and FY_r have dimensions of [m2]
!---------------------------------------------------------------------

      FX_r(inea) = bedld*angleu*dt
      FY_r(inea) = bedld*anglev*dt

      IF(FX_r(inea)/=FX_r(inea)) THEN
        WRITE(errmsg,*)'FX_r0 is NaN',myrank,inea,FX_r(inea),bedld,  &
        &              angleu,dt
        CALL parallel_abort(errmsg)
      ENDIF

      IF(FY_r(inea)/=FY_r(inea)) THEN
        WRITE(errmsg,*)'FY_r0 is NaN',myrank,inea,FY_r(inea),bedld,  &
        &              anglev,dt
        CALL parallel_abort(errmsg)
      endif


!--------------------------------------------------------------------
!- Bed slope effects (Lesser et al. 2004)
!--------------------------------------------------------------------
      IF (slope_formulation==4) THEN
        CALL bedslope_effects_lesser2004(inea,ised,FX_r(inea),FY_r(inea))
      END IF


!---------------------------------------------------------------------
! - Sanity check
!---------------------------------------------------------------------

      IF (FX_r(inea)/=FX_r(inea)) THEN
        WRITE(errmsg,'(A,I5,A,E15.8,A,E15.8,A,E15.8,A,E15.8,A,E15.8)')&
        &     'FX_r NaN nea:',inea,' bustr:',bustr(inea),             &
        &' bvstr:',bvstr(inea),' cff:',cff,' cff1:',cff1,' cff2:',cff2
        CALL parallel_abort(errmsg)
      ENDIF

      IF (FY_r(inea)/=FY_r(inea)) THEN
        WRITE(errmsg,*)'FY_r1 is NaN',myrank,inea,FY_r(inea),        &
        &              a_slopey,a_slopex
        CALL parallel_abort(errmsg)
      ENDIF

!---------------------------------------------------------------------
      END SUBROUTINE sed_bedload_vr

!=====================================================================
!=====================================================================

      SUBROUTINE sed_bedload_mpm()

! this whole section is currently unused
! commented for better overview
!... Compute critical stress for horizontal bed
!... For Meyer-Peter Muller formulation tauc0= 0.047.
!... Compute bottom stress and tauc. Effect of bottom slope (Carmo,1995).
!            tauc0 =tau_ce(ised)
!            alphas=datan(-(dzdx*angleu+dzdy*anglev))

!... Magnitude of bed load at rho points. Meyer-Peter Muller formulation.
!... bedld has dimensions  (m2 s-1)
!!!#if defined CARMO
!!!          if(slope_formulation==sf_carmo)
!!!            call stress(angleu,anglev,dzdx,dzdy,alphas,sed_angle,tauc0,tauc,i)
!!!            Eq. 46 q_{b,q} = 8*(tau_k
!!!            bedld=8.0_r8*(MAX((tau_w(i)-tauc)/smgd,0.0_r8)**1.5_r8)*smgdr
!!!#elif defined SOULSBY
!!!          elseif(slope_formulation==sf_soul)
!!!            call stress_soulsby(angleu,anglev,dzdx,dzdy,alphas,sed_angle,tauc0,tauc,i)
!!!            bedld=8.0_r8*(MAX((tau_w(i)-tauc)/smgd,0.0_r8)**1.5_r8)*smgdr
!!!#elif defined DELFT
!!!          elseif(slope_formulation==sf_del)
!!!            bedld=8.0_r8*(MAX((tau_w(i)*osmgd-0.047_r8),0.0_r8)**1.5_r8)*smgdr!!!          endif

!!!#elif defined DAMGAARD
!!!            if (abs(alphas)>datan(sed_angle))alphas=datan(sed_angle)*SIGN(1.0_r8,alphas)
!!!            tauc=tauc0*(sin(datan(sed_angle)+(alphas))/sin(datan(sed_angle)))
!!!            bedld=8.0_r8*(MAX((tau_w(i)-tauc)/smgd,0.0_r8)**1.5_r8)*smgdr

!!!            if (alphas>0)then
!!!                 cff=1
!!!             else if (alphas<=0)then
!!!                 cff=1+0.8*((tauc0/tau_w(i))**0.2)*(1-(tauc/tauc0))**(1.5+(tau_w(i)/tauc0))
!!!             endif
!!!             bedld=bedld*cff
!!!             if (time>1500) then
!!!!ZYL
!!!               if (isnan(cff)==.true.) call parallel_abort('cff is NaN "i" 1')
!!!               if (isnan(bedld)==.true.) call parallel_abort('bedld is NaN "i"')
!!!!"
!!!             endif
!!!#endif  ! bedld

!!!!-------------------------------------------------------------------
!!!!... Partition bedld into x  and y directions, at the center of each 
!!!!... element (rho points), and integrate in time.
!!!!... FX_r and FY_r have dimensions of m2
!!!!-------------------------------------------------------------------

!!!#if defined CARMO || defined SOULSBY
!!!          if(slope_formulation==3) ! || soulsby
!!!            FX_r(i)=bedld*dcos(alphas)*angleu*dt
!!!            FY_r(i)=bedld*dcos(alphas)*anglev*dt
!!!#elif defined DAMGAARD || defined DELFT
!!!          else
!!!            FX_r(i)=bedld*angleu*dt
!!!            FY_r(i)=bedld*anglev*dt
!!!          endif
!!!#endif !CARMO||SOULSBY

!!!#ifdef DELFT
!!!          if(slope_formulation== 
!!!!... Bed_slope effects
!!!!... longitudinal bed slope
!!!!... limit slope to 0.9*(sed_angle)

!!!            cff=(dzdx*angleu+dzdy*anglev)
!!!            cff1=min(abs(cff),0.9*sed_angle)*sign(1.0_r8,cff)

!!!!... add contribution of longitudinal bed slope to bed load

!!!            cff2=datan(cff1)
!!!            a_slopex=1+1*((sed_angle/(cos(cff2)*(sed_angle-cff1)))-1)

!!!            FX_r(i)=FX_r(i)*a_slopex
!!!            FY_r(i)=FY_r(i)*a_slopex

!!!!... Transverse bed slope

!!!            cff=(-(dzdx*anglev)+dzdy*angleu)
!!!            a_slopey=1.5*sqrt(tauc0/(abs(bustr(i))+abs(bvstr(i))))*cff
!!!            FX_r(i)=FX_r(i)-(FY_r(i)*a_slopey)
!!!            FY_r(i)=FY_r(i)+(FX_r(i)*a_slopey)
!!!#endif !DELFT
      END SUBROUTINE sed_bedload_mpm

!=====================================================================


!=====================================================================
!--------------------------------------------------------------------
      SUBROUTINE sed_bedload_wl14(ised,inea)
!--------------------------------------------------------------------
! This subroutine computes bedload load sediment
! transport rate (m^2/s) from Wu and Lin (2014) which is a formulation
! created for nonuniform (i.e. multi-class) sediment transport.
! The output flux is integrated over the time step (m^2).
!
! Author: thomas guerin (thomas.guerin@univ-lr.fr)
! Date: 12/2014
!
! History: - 2019/02 - B.Mengual : Implementation of the original code
!                                  of T.Guerin in Sed2d
!          - 2020/02 - B.Mengual: > use of ustress/vstress, wdir       
!                         > wave asymmetry effect: rewritten and       
!                           introduction of a limiter w_asym_max       
!                         > bed slope effects: now under option with   
!                           an external call to a new subroutine       
!--------------------------------------------------------------------

      USE sed_mod
      USE schism_glbl, ONLY: rkind,dt,kbe,pi,i34,elnode,eta2,dpe,uu2,vv2,&
                             out_wwm,h0
      USE schism_msgp, ONLY: parallel_abort
 
      IMPLICIT NONE
   
!- Arguments --------------------------------------------------------
      INTEGER, INTENT(IN)  :: ised,inea
!- Constants --------------------------------------------------------
      REAL(rkind), PARAMETER :: Ar = 12.d0
      REAL(rkind), PARAMETER :: shields_cr = 0.03d0
      REAL(rkind), PARAMETER :: m = 0.6d0
      REAL(rkind), PARAMETER :: nu=1.36d-6   ! Cinematic viscosity

      INTEGER,     PARAMETER :: top = 1      ! Top layer of bed
!- Local variables --------------------------------------------------
      REAL(rkind) :: d,d50,d90,D_star,htot,hs_e,tp_e,wdir_e,&
                     &wlp_e,kb1,s,cff,cff1,cff2,&
                     &a_slopex,a_slopey,beta
      REAL(rkind) :: ac,at,Aw,Cd,Deltar_c,Deltar_c_mm,Deltar_w,f_grain_c,&
                    &f_grain_cw,f_grain_w,fw,kpeak,ks_c,ks_form_c,       &
                    &ks_form_w,ks_grain,ks_w,Lambdar_c,Lambdar_c_mm,     &
                    &Lambdar_w,n,n_grain,omega,pe,ph,phi,qb_off,qb_on,   &
                    &qs_c,rw,tau_b,tau_b_c,tau_b_wm,    &
                    &tau_cr,tau_grain_b_off,tau_grain_b_on,              &
                    &tau_grain_b_wm_off,tau_grain_b_wm_on,theta_off,     &
                    &theta_on,tmp,Uc,udir,Uw,Uwm_off,   &
                    &Uwm_on,Ws,Xu
      REAL(rkind), DIMENSION(ntr_l) :: dclass,F_class
      INTEGER :: i,j,inode
      
!--------------------------------------------------------------------
   
!--------------------------------------------------------------------
!- Sediment parameters
!--------------------------------------------------------------------
      d=Sd50(ised)
      dclass=Sd50
      d50=bottom(inea,isd50)
      d90=2.5d0*d50
      s=Srho(ised)/rhom
      D_star=Sd50(ised) * (g*(s-1.d0)/nu**2.d0)**(1.d0/3d0)
      F_class=bed_frac(top,inea,ised)
   
!--------------------------------------------------------------------
!- Water depth and wave conditions
!--------------------------------------------------------------------  
      htot = dpe(inea)+sum(eta2(elnode(1:i34(inea),inea)))/i34(inea)
      tp_e = tp(inea)
      hs_e = hs(inea)
      wlp_e = wlpeak(inea)
      wdir_e = wdir(inea)
   
!--------------------------------------------------------------------
!- Current amplitude and direction
!--------------------------------------------------------------------

      ! Reference code
      !kb1 = kbe(inea)+1
      !cff1 = sum(uu2(kb1,elnode(1:i34(inea),inea)))/i34(inea)
      !cff2 = sum(vv2(kb1,elnode(1:i34(inea),inea)))/i34(inea)

      ! BM: see ustress, vstress estimates in sediment.F90
      cff1 = ustress(inea)
      cff2 = vstress(inea)

      Uc = dsqrt(cff1*cff1+cff2*cff2)
      udir = DATAN2(cff2,cff1)
   
!--------------------------------------------------------------------
!- Angle between current and wave directions
!--------------------------------------------------------------------
#ifdef USE_WWM
      phi = wdir_e - udir
#else
      phi = 0.d0
#endif
   
   
!- Wave parameters
      IF (hs_e>0.01d0.and.tp_e>0.01d0.and.wlp_e>0.01d0) THEN
        !- wave angular frequency
        omega = 2.d0*pi/tp_e !>0
    
        !- peak wave number
        kpeak = 2.d0*pi/wlp_e !>0
      ENDIF
   
!- Wave orbital velocity calculation
      IF (hs_e>0.01d0.and.tp_e>0.01d0.and.wlp_e>0.01d0.and.htot.GT.h0) THEN
        IF (iasym==0) THEN ! Airy waves
          Uw=U_crest(inea)
        ELSE IF (iasym==1) THEN ! Elfrink et al. (2006)
          Uw = (U_crest(inea)+U_trough(inea))/2.d0 ! >0
        END IF
      ELSE !no waves
        Uw = 0.d0
      ENDIF

!- Bed-load transport
   
      !exposed probability of particles
      pe = 0.d0
      do i=1,ntr_l
         pe = pe + F_class(i)*d/(d+dclass(i))
      enddo
    
      !hidden probability of particles
      ph = 1.d0 - pe
    
      !critical shear stress for the incipient motion of corresponding
      !size class
      tau_cr = g*(Srho(ised)-rhom)*d*shields_cr*(pe/ph)**(-m) ! >0
      !NB: grav forgotten in Wu and Lin 2014 (Eq. 24)
   
      !equivalent roughness height calculation
      ks_grain = 1.5d0*d90 !grain roughness height
      Lambdar_c = 1000.d0*d50 !ripple wavelength (m) according to Soulsby (97)
      Deltar_c = Lambdar_c/7.d0 !ripple height (m) according to Soulsby (97)
      ! Lambdar_c_mm = 245.d0*(d50*1000.d0)**0.35d0 !ripple wavelength (mm)
      !Lambdar_c = Lambdar_c_mm/1000.d0 !ripple wavelength (m)
      !Deltar_c_mm = 0.074d0*Lambdar_c_mm*(d50*1000.d0)**(-0.253d0) !ripple
      !height (mm)
      !Deltar_c = Deltar_c_mm/1000.d0 !ripple height (m)
      ks_form_c = Ar*Deltar_c**2.d0/Lambdar_c ! >0
      ks_c = ks_grain + ks_form_c !equivalent roughness height
    
      !grain friction coefficient under current only
      n_grain = 1.d0/20.d0*d50**(1.d0/6.d0) !Manning coef of grain roughness
      n = htot**(1.d0/6.d0)/(18.d0*dlog(12.d0*htot/ks_c)) !Manning coef of bed roughness
      f_grain_c = 2.d0*(n_grain/n)**(3.d0/2.d0)*g*n**2.d0             &
                 &/(htot**(1.d0/3.d0))
    
      !grain friction coefficient under waves only
      Aw = Uw*tp_e/2.d0/pi !wave excursion
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
      if (hs_e>0.01d0.and.tp_e>0.01d0) then
        rw = U_crest(inea)/Uw - 1.d0 !wave asymmetry coef
        ac = pi*T_crest(inea)/tp_e
        at = pi*T_trough(inea)/tp_e
    
        tmp = 0.5d0*rhom*f_grain_w*Uw**2.d0/2.d0
    
        tau_grain_b_wm_on = tmp*(1.d0 + rw**2.d0 +                       &
          &13.d0/6.d0*rw*dsin(ac)/ac + 1.d0/6.d0*dsin(2.d0*ac)/(2.d0*ac))
    
        tau_grain_b_wm_off = tmp*(1.d0 + rw**2.d0 -                      &
          &13.d0/6.d0*rw*dsin(at)/at + 1.d0/6.d0*dsin(2.d0*at)/(2.d0*at))
      else
        tau_grain_b_wm_on = 0.d0
        tau_grain_b_wm_off = 0.d0
      endif
   
      !root-mean-square values of wave orbital velocity over the two half
      !cycles
      !if (hs_e>0.d0.and.tp_e>0.d0) then
      if (f_grain_w>0.0d0) then !BM
        tmp = 0.5d0*rhom*f_grain_w ! >0
        !if (tmp==0.d0) then
        !  write(12,*)'f_grain_w',f_grain_w,'Aw',Aw,'Uw',Uw,'tp_e',tp_e
        !  write(12,*)'ks_grain',ks_grain
        !  call parallel_abort('sed_bedload_wl14(3)')
        !endif
        Uwm_on = dsqrt(tau_grain_b_wm_on/tmp)
        Uwm_off = dsqrt(tau_grain_b_wm_off/tmp)
      else
        Uwm_on = 0.d0
        Uwm_off = 0.d0
      endif
   
      !bed grain shear stress due to combined waves and current
      tmp = 0.5d0*rhom*f_grain_cw
      tau_grain_b_on = tmp*(Uc**2.d0 + Uwm_on**2.d0                      &
                           & + 2.d0*Uc*Uwm_on*dcos(phi))
      tau_grain_b_off = tmp*(Uc**2.d0 + Uwm_off**2.d0                    &
                           & + 2.d0*Uc*Uwm_off*dcos(pi-phi))
   
      !resultant bed-load transport components during the two half cycles
      tmp = 0.0053d0*dsqrt((s-1.d0)*g*d**3.d0)
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
   
      !final x and y components of bed-load transport [m2/s]
      if (hs_e>0.01d0.and.tp_e>0.01d0) then
        FX_r(inea) = T_crest(inea)/tp_e*dcos(theta_on)*qb_on + &
               &T_trough(inea)/tp_e*dcos(theta_off)*qb_off
    
        FY_r(inea) = T_crest(inea)/tp_e*dsin(theta_on)*qb_on + &
               &T_trough(inea)/tp_e*dsin(theta_off)*qb_off
      else !only current
        FX_r(inea) = dcos(udir)*qb_on
        FY_r(inea) = dsin(udir)*qb_on
      endif
    
      !Integration over time step [m2]
      FX_r(inea) = FX_r(inea)*dt
      FY_r(inea) = FY_r(inea)*dt


!--------------------------------------------------------------------
!- Bed slope effects (Lesser et al. 2004)
!--------------------------------------------------------------------
      IF (slope_formulation==4) THEN
        CALL bedslope_effects_lesser2004(inea,ised,FX_r(inea),FY_r(inea))
      END IF


!--------------------------------------------------------------------
!- Sanity check
!--------------------------------------------------------------------
      IF(FX_r(inea)/=FX_r(inea) .OR. FY_r(inea)/=FY_r(inea)) THEN
        WRITE(12,*)'FX_r or FY_r is NaN'
        WRITE(12,*)'ised,inea',ised,inea
        WRITE(12,*)'d,dclass,d50',d,dclass,d50
        WRITE(12,*)'D_star,F_class',D_star,F_class
        WRITE(12,*)'htot,tp_e,hs_e,wlp_e',htot,tp_e,hs_e,wlp_e
        WRITE(12,*)'cff1,cff2,Uc,udir,wdir,phi',cff1,cff2,Uc,udir,wdir,phi
        WRITE(12,*)'omega,kpeak',omega,kpeak
        WRITE(12,*)'Uw,U_crest,T_crest,T_trough',Uw,U_crest(inea),T_crest(inea),T_trough(inea)
        WRITE(12,*)'pe,ph,tau_cr',pe,ph,tau_cr
        WRITE(12,*)'ks_form_c,ks_c',ks_form_c,ks_c
        WRITE(12,*)'f_grain_c,f_grain_w,Xu,f_grain_cw',f_grain_c,f_grain_w,Xu,f_grain_cw
        WRITE(12,*)'rw,ac,at',rw,ac,at
        WRITE(12,*)'tau_grain_b_wm_on,tau_grain_b_wm_off',tau_grain_b_wm_on,tau_grain_b_wm_off
        WRITE(12,*)'Uwm_on,Uwm_off',Uwm_on,Uwm_off
        WRITE(12,*)'tau_grain_b_on,tau_grain_b_off',tau_grain_b_on,tau_grain_b_off
        WRITE(12,*)'qb_on,qb_off',qb_on,qb_off
        WRITE(12,*)'theta_on,theta_off',theta_on,theta_off
        WRITE(12,*)'FX_r,FY_r',FX_r,FY_r
      ENDIF


      END SUBROUTINE sed_bedload_wl14
!=====================================================================


!=====================================================================
!--------------------------------------------------------------------
      SUBROUTINE bedload_flux_limiter
!--------------------------------------------------------------------
! This subroutine limits the bedload rate components
! FX_r/FY_r (in m^2) according to the sediment mass available within
! the active layer.
! The unit of the output flux is unchanged (m^2; i.e. integrated over
! the time step).
!
! Author: Baptiste Mengual (baptiste.mengual@univ-lr.fr)
! Date: 02/2020
!--------------------------------------------------------------------

      USE sed_mod
      USE schism_glbl, ONLY : rkind,i34,nea,idry_e,xnd,ynd,area

      IMPLICIT NONE

!- Local variables --------------------------------------------------!
      INTEGER     :: i,j
      REAL(rkind) :: dist_tmp,dx,dy,FN_r,Fmax

      DO i=1,nea
        IF(idry_e(i)==1) CYCLE
        dist_tmp=0.0d0
        DO j=1,i34(i)-1
          dx=xnd(j+1)-xnd(j)
          dy=ynd(j+1)-ynd(j)
          dist_tmp=dist_tmp+(sqrt(dx*dx+dy*dy)/i34(i))
        END DO
        FN_r=sqrt(FX_r(i)*FX_r(i) + FY_r(i)*FY_r(i))
        ! Fmax in m2
        Fmax=(1.0d0/dist_tmp)*(1.0d0-bed(1,i,iporo))*area(i)*bottom(i,iactv)
        IF (FN_r > Fmax) THEN
          FX_r(i)=(FX_r(i)/FN_r)*Fmax
          FY_r(i)=(FY_r(i)/FN_r)*Fmax
        END IF
      END DO

      END SUBROUTINE bedload_flux_limiter
!=====================================================================




!=====================================================================

      SUBROUTINE bedchange_bedload(ised,it,moitn,mxitn,rtol,qsan,    &
      &                            hbed,hbed_ised)    
!--------------------------------------------------------------------!
! This computes changes in bed characteristics and thickness induced !
! by bedload transport                                               !
!                                                                    !
! Author: Knut Kraemer                                               !
! Date: 27/11/2012                                                   !
!                                                                    !
! History: 2012/12 - F.Ganthy : form homogenisation of sediments     !
!          routines
!          2020/02 - A.de Bakker : Implementation of Thomas Guerin   !
!                    code to solve the Exner equation with a WENO    !
!                    scheme (imeth_bed_evol=2)                       !
!          2020/02 - B.Mengual : > morph_fac applied to hbed_ised    !
!               instead of FX_r/FY_r (CFL issue for large morph_fac) !
!               > FX_r/FY_r not multiplied by bed_frac(1,i,ised),    !
!                 already done in sediment.F90                       !
!--------------------------------------------------------------------!

      USE sed_mod
      USE schism_glbl, ONLY : rkind,i34,elnode,nea,npa,nne,idry,idry_e,xctr,   &
     &yctr,np,indel,errmsg,elside,iself,nxq,xcj,ycj,ics,xel,yel,mnei_p
      USE schism_msgp
      
      IMPLICIT NONE

!- Local variables --------------------------------------------------!

      INTEGER, INTENT(IN) :: ised,it ! sediment class and time step
      INTEGER, INTENT(IN) :: moitn,mxitn ! JCG solver iterations
      REAL(rkind),INTENT(IN) :: rtol     ! Relative tolerance
      REAL(rkind),DIMENSION(npa),INTENT(OUT) :: qsan,hbed,hbed_ised

      INTEGER     :: i,j,nm1,nm2,nm3,ks,k,ie,id,id1,id2,id3,isd1,isd2
      REAL(rkind) :: yp,xp,flux,cff,cff1,cff2,cff3,cff4

      REAL(rkind),DIMENSION(np) :: bed_poro 
      
!- Start Statement --------------------------------------------------!
      
!      IF(myrank.EQ.0) WRITE(16,*)'SED: Entering bedchange_bedload'
      
!---------------------------------------------------------------------
! -  Apply morphology factor to bedload transport
!   This routine is called inside a ised-loop, and F[XY]_r will be
!   overwritten for each ised-loop
!---------------------------------------------------------------------

      !BM morph_fac will be applied to hbed_ised and FX_r/FY_r are 
      ! already multiplied by bed_frac(1,i,ised) in sediment.F90
      ! -> entire loop commented
      !do i = 1,nea
        !IF(idry_e(i)==1) CYCLE
        !FX_r(i) = FX_r(i)*morph_fac(ised)*bed_frac(1,i,ised) !m^2
        !FY_r(i) = FY_r(i)*morph_fac(ised)*bed_frac(1,i,ised)
      !enddo !i

!---------------------------------------------------------------------
! -  qsan=\int q_n*dt d\Gamma (although dimensioned up to npa, only 1:np
! are used)
! qsaxy is the sand flux at the element center integrated in time. Now 
! compute the line integral of qsaxy*normal along the control volume. 
! Add for each node in qsan. The unit normal is directed outward.
!---------------------------------------------------------------------

     ! Anouk; implemented Thomas' code of the weno scheme --> imeth_bed_evol=2
     IF (imeth_bed_evol == 1) THEN  
       qsan=0.0d0 ![m^3]
       do i=1,np !resident only required
         do j=1,nne(i)        
           ie=indel(j,i)
           id=iself(j,i)
           if(idry_e(ie)==1) cycle
 
           !Wet elem.
           if(ics==1) then
             isd1=elside(nxq(i34(ie)-2,id,i34(ie)),ie)
             isd2=elside(nxq(i34(ie)-1,id,i34(ie)),ie)
             qsan(i)=qsan(i)+FX_r(ie)*(ycj(isd1)-ycj(isd2))+FY_r(ie)*(xcj(isd2)-xcj(isd1)) !m^3
           else !ll
             id1=nxq(1,id,i34(ie))
             id3=nxq(i34(ie)-1,id,i34(ie))
             cff1=(xel(id,ie)+xel(id1,ie))/2 !xcj(isd2)-> between i and id1   
             cff2=(yel(id,ie)+yel(id1,ie))/2 !ycj(isd2)-> between i and id1   
             cff3=(xel(id,ie)+xel(id3,ie))/2 !xcj(isd1)-> between i and id3   
             cff4=(yel(id,ie)+yel(id3,ie))/2 !ycj(isd1)-> between i and id3   
             qsan(i)=qsan(i)+FX_r(ie)*(cff4-cff2)+FY_r(ie)*(cff1-cff3) !m^3
           endif !ics
         enddo !j
       enddo !i=1,np
     ELSEIF (imeth_bed_evol == 2) THEN
       CALL weno_scheme(it,qsan)
     ELSE
       WRITE (errmsg,*)'SED: Scheme not existing'
     ENDIF


!---------------------------------------------------------------------
! - Compute erosion rates
! change bed(1,mne,iporo) from elements to nodes 
! initalize before adding
!---------------------------------------------------------------------
      bed_poro = 0.0d0
      DO i=1,np
        IF(idry(i)==1) CYCLE

        ks=0 !Number of wet neighbor elements
        DO j=1,nne(i)
          k = indel(j,i)
          IF(idry_e(k)==1) CYCLE
          ks = ks+1
          bed_poro(i) = bed_poro(i)+bed(1,k,iporo) 
        ENDDO ! End loop nne
        
        IF(ks==0) CALL parallel_abort('SEDIMENT: (2)')
        bed_poro(i) = bed_poro(i)/ks
      ENDDO ! End loop np

      ! Exchange ghosts
      !call exchange_p2d(bed_poro(:))

!---------------------------------------------------------------------
! - Take porosity into account and adjust qsan accordingly
!jl. mcoefd      [m2] ----  aux1/aux2 are related to Exner equation 
!    hdbed_ised  [m]
!    qsan        [m3]
!    A      x          = B
!    mcoefd hdbed_ised = qsan
!    [m2]   [m]        = [m3]
!---------------------------------------------------------------------
!jl. Why is qsan divided by (1-porosity)?  Doesn't this magically
!   dd volume? I think the idea is to scale the value by (1-porosity).
!FG. qsan is divided by (1-porosity) to take account for empty space
!    between grain (i.e. porosity) on bed level change

      qsan(1:np) = qsan(1:np)/(1-bed_poro(1:np))

      ! Use JCG solver to calc hbed_ised
      hbed_ised=0.0d0 !initial guess
      !Right now zero-flux b.c. is imposed
      CALL solve_jcg(mnei_p,np,npa,it,moitn,mxitn,rtol,mcoefd,hbed_ised,qsan,bc_sed,lbc_sed)

!---------------------------------------------------------------------
! - Bed/bottom level change due to bedload transport in [m]
!---------------------------------------------------------------------
      !aggregate over all sed. classes (this routine is called inside a sed class loop)
      !hbed(:) = hbed(:)+hbed_ised(:)
      !BM morphological factor is now applied
      hbed(:) = hbed(:)+(hbed_ised(:)*morph_fac) 

      ! Consistency check
      DO i=1,np
        IF (hbed(i)/=hbed(i)) THEN
          WRITE(errmsg,*)'hbed(i) is NaN',myrank,i,hbed(i),qsan(i),bed_poro(i)
          CALL parallel_abort(errmsg)
        ENDIF
      ENDDO

!---------------------------------------------------------------------
! - Evaluate depth at the elements and changes in bed properties
!---------------------------------------------------------------------
      DO i=1, nea
        IF(idry_e(i)==1) CYCLE

        ! Average bed change in element [m]
        cff=sum(hbed_ised(elnode(1:i34(i),i)))/i34(i)
        ! Average bed change in element [kg/m2]
        cff1=cff*Srho(ised)

!---------------------------------------------------------------------
! - Update bed mass [kg/m2] according to bed change
!jl. You can lose a maximum of bed_mass per sediment class per time 
!    step based on changes to hbed_ised.  Is there a lower limit on 
!    hbed_ised?...
!    No lower limit on hbed_ised, but bed_mass and bed(1, nea, ithick)     
!    have a lower bound of 0.
!---------------------------------------------------------------------
!YJZ: should be -cff1?
        !bed_mass(1,i,nnew,ised)=MAX(bed_mass(1,i,nstp,ised)+cff1,0.0d0)
        bed_mass(1,i,nnew,ised)=MAX(bed_mass(1,i,nstp,ised)-cff1,0.0d0)

       !Jan what is this for?
       !Save bed mass for next step
        IF(suspended_load==1) THEN
          DO k=2,Nbed
            bed_mass(k,i,nnew,ised) = bed_mass(k,i,nstp,ised)
          ENDDO
        ENDIF

!---------------------------------------------------------------------
! - Update top layer thickness according to bed change in [m]
!---------------------------------------------------------------------

        !bed(1,i,ithck) = MAX((bed(1,i,ithck)+cff),0.0d0)
        bed(1,i,ithck) = MAX(bed(1,i,ithck)-cff,0.0d0)
      ENDDO !i=1,nea


!--------------------------------------------------------------------!
      END SUBROUTINE bedchange_bedload


!--------------------------------------------------------------------
      SUBROUTINE wave_asymmetry_Elfrink(H,wave_per,w_dir,depth,dhxi,dhyi,&
                       & Ucrest,Utrough,Tcrest,Ttrough,Uorbi)
!--------------------------------------------------------------------
! This subroutine computes wave asymmetry based on Elfrink et al.
! (2006, Coastal Engineering)
!
! NB: offshore wave height and wave period are taken into account to
! compute irribaren number and offshore wave length
!
! Author: thomas guerin (thomas.guerin@univ-lr.fr)
! Date: 26/04/2013
!
! Updates: 18/12/2014 - Slope calculation corrected
!                     - Checks added
!--------------------------------------------------------------------
 
      USE sed_mod, ONLY : ech_uorb  
      USE schism_glbl, ONLY : pi,errmsg,ielg
      USE schism_msgp, ONLY : parallel_abort
   
      IMPLICIT NONE
   
!- Arguments --------------------------------------------------------
      REAL(8), INTENT(IN) :: H,wave_per,w_dir,depth,dhxi,dhyi
      REAL(8), INTENT(OUT) :: Ucrest,Utrough,Tcrest,Ttrough
      REAL(8), DIMENSION(ech_uorb+1), INTENT(OUT) :: Uorbi
!- Constants --------------------------------------------------------
      REAL(8), PARAMETER :: g = 9.80665d0
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
      IF (depth<=0.d0.or.wave_per<=0.d0) THEN
      WRITE (errmsg,*)'depth',depth,'wave_per',wave_per 
      CALL parallel_abort(errmsg)
      ENDIF
      psi = 4.d0*pi**2.d0*depth/(g*wave_per**2.d0) 
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
      if (dabs(dcos(w_dir))>0.d0) then
        Slope = dcos(w_dir)*(dhxi + dhyi*dtan(w_dir))
      else
        if (dsin(w_dir)==1.d0) then
          Slope = dhyi
        elseif (dsin(w_dir)==-1.d0) then
          Slope = -dhyi
        endif
      endif
!  Slope = -DSQRT(dhxi**2.d0+dhyi**2.d0)*DCOS(w_dir+DATAN2(dhyi,dhxi))
      IF (Hadim<=0.d0.or.Ladim<=0.d0) THEN
        CALL parallel_abort('SED_TRANS: (2)')
      END IF
      Zeta = DTAN(Slope)/DSQRT(Hadim/Ladim)
      Ur = Hadim*Ladim**2.d0 ! >0
   
!- Compute the normalized maximal orbital velocity
      C1 = Ladim-10.d0
      C2 = DABS(C1-(Hadim-DABS(Zeta)))
      C3 = Zeta*(1.d0-C1)
      C4 = DTANH(DABS(C3-C2)/Ur)
      tmp=DABS(Zeta)+DTANH(C4)
      if(Hadim<=0.d0.or.tmp<0.d0) call parallel_abort('SED_TRANS: (3)')
!Ur>0
      C5 = DSQRT(tmp) !DABS(Zeta)+DTANH(C4))
      P1 = DSQRT(Hadim)-C5*Hadim
      U1 = b1*P1+a1
   
!- Compute the velocity asymmetry parameter
      D1 = 3.d0*Zeta+2.d0*Ladim/Ur
      D2 = DSQRT(Ladim)-DTANH(DABS(D1))
      if (D2<0.d0) then
        D2 = 0.d0
      endif
      D3 = (2.d0*Zeta+DSQRT(Ladim/Ur))**2.d0 ! >0
      D4 = Ur+Ladim/D3/Ur ! >0
      D5 = DSQRT(D2/D4)
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
   
!WRITE(16,*)'sed_bedload','T0',T0,'T1',T1,'Ucrest',Ucrest,'Utrough',Utrough

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
!WRITE(16,*)'sed_bedload','Ustar',Ustar,'Uairy',Uairy,'U2',U2,'U1',U1,'P1',P1,'Hadim',Hadim
   
      Ucrest = U1*Ustar
      Utrough = Ustar-Ucrest
   
!WRITE(16,*)'sed_bedload','Utrough',Utrough,'Ucrest',Ucrest,'T0',T0,'T1',T1
      U0 = (Ucrest*T0-Utrough*(1.d0-T0))/(T0-T1)
   
! WRITE(16,*)'sed_bedload','U0',U0,'Ucrest',Ucrest
         
      IF (U0 .GT. (0.25d0*Ucrest)) THEN
          U0 = 0.25d0*Ucrest
      ENDIF
   
      tmp=T0*(U0-Ucrest-Utrough)

      IF (tmp==0.d0) THEN

        WRITE (errmsg,*)'T0',T0,'U0',U0,'Ucrest',Ucrest,'Utrough',Utrough
        WRITE(16,*)'sed_bedload','U0',U0,'Ucrest',Ucrest
        WRITE(16,*)'Ladim Hadim Zeta ',Ladim,Hadim,Zeta

        Ucrest=0.0d0
        Utrough=0.0d0
        Tcrest=0.0d0
        Ttrough=0.0d0
        Uorbi(:)=0.0d0

      ELSE

        a1bis = (-Utrough+T1*U0)/tmp !(T0*(U0-Ucrest-Utrough))
        IF (a1bis .LT. 0.99d0) THEN
          T0 = 0.99d0*T0
          if(T0==1) call parallel_abort('SED_TRANS: (6)')
          Utrough = (-(T0-T1)*U0+T0*Ucrest)/(1.d0-T0)
        ELSE
          T0 = a1bis*T0
        ENDIF
  
        !- Compute wave half-periods
        Tcrest = T0*wave_per
        Ttrough = wave_per-Tcrest
  
        !- Rebuilt of orbital velocity over one wave period
        if(ech_uorb==1) call parallel_abort('SED_TRANS: (7)')
        DO t=0,ech_uorb
          tv=DBLE(t)/DBLE(ech_uorb)
          IF ((tv .GE. 0.d0) .AND. (tv .LE. T1)) THEN
            if(T1==0) call parallel_abort('SED_TRANS: (8)')
            Uorbi(t) = Ucrest*DSIN(pi/2.d0*tv/T1)
          ELSEIF ((tv .GT. T1) .AND. (tv .LE. T0)) THEN
            if(T1==T0) call parallel_abort('SED_TRANS: (9)')
            Uorbi(t) = Ucrest*DCOS(pi/2.d0*(tv-T1)/(T0-T1))  &
                       - U0*DSIN(pi*(tv-T1)/(T0-T1))
          ELSEIF ((tv .GT. T0) .AND. (tv .LE. T2)) THEN
            if(T2==T0) call parallel_abort('SED_TRANS: (10)')
            Uorbi(t) = -Utrough*DSIN(pi/2.d0*(tv-T0)/(T2-T0))
          ELSEIF ((tv .GT. T2) .AND. (tv .LE. 1.d0)) THEN
            if(T2==1) call parallel_abort('SED_TRANS: (11)')
            Uorbi(t) = -Utrough*DCOS(pi/2.d0*(tv-T2)/(1.d0-T2))
          ENDIF
        ENDDO

      END IF

   
      END SUBROUTINE wave_asymmetry_Elfrink
!--------------------------------------------------------------------

  SUBROUTINE weno_scheme(it,qsan)
!--------------------------------------------------------------------
! This subroutine computes bottom evolution using a WENO-base scheme
! to solve the Exner equation
!
! Authors: Thomas Guerin   (thomas.guerin@univ-lr.fr)
! Date: 11/2016
!--------------------------------------------------------------------

  USE schism_glbl, ONLY : area,dp,dt,elnode,eta2,idry,idry_e,indel,moitn0,mxitn0,nea,nne,np,npa,rkind,rtol0,xctr,xnd,yctr,ynd
  USE schism_msgp, ONLY : exchange_p2d
  USE sed_mod, ONLY : FX_r,FY_r,vc_area

  IMPLICIT NONE

!- Arguments --------------------------------------------------------
  INTEGER, INTENT(IN) :: it
  REAL(rkind), DIMENSION(npa), INTENT(OUT) :: qsan
!- Parameters -------------------------------------------------------
  INTEGER, PARAMETER :: nb_comb = 8
  REAL(rkind), PARAMETER :: ONETHIRD = 1.d0/3.d0
  REAL(rkind), PARAMETER :: c = 0.5d0 + DSQRT(3.d0)/6.d0
  REAL(rkind), PARAMETER :: epsi = 1.e-10
  REAL(rkind), PARAMETER :: rp = 1.d0
  REAL(rkind), PARAMETER :: beta_comp = 2.d0
!- Local variables --------------------------------------------------
  INTEGER :: el1,el2,el3,i,j,k,m,nd1,nd2,nd3,p
  REAL(rkind) :: alpha,flux,h1,h2,h3,htot,OItmp,OIx,OIxtmp,OIy,OIytmp,          &
                 P1,P1G1,P1G2,P1nd,P2,P2G1,P2G2,P2nd,Px1G1,Px1G2,Px1nd,Px2G1,Px2G2,Px2nd,  &
                 Py1G1,Py1G2,Py1nd,Py2G1,Py2G2,Py2nd,r_flux,w_sum,xG1,xG2,xi,xnd1,xnd2, &
                 xnd3,xp,yG1,yG2,yi,ynd1,ynd2,ynd3,yp
  REAL(rkind), DIMENSION(nea,2) :: qdt_e
  REAL(rkind), DIMENSION(npa,56,3) :: Cx,Cy
  REAL(rkind), DIMENSION(3) :: Cxtmp,Cytmp
  REAL(rkind), DIMENSION(56) :: OI
  REAL(rkind), DIMENSION(npa,nb_comb) :: w
!--------------------------------------------------------------------


!- Initialization
  qdt_e(:,1) = FX_r
  qdt_e(:,2) = FY_r
  Cx(:,:,:) = 0.d0; Cy(:,:,:) = 0.d0
  w(:,:) = 0.d0
  DO i=1,np
     IF(idry(i)==1) CYCLE

     IF(nne(i)<3) CYCLE

     !Compute linear polynom for all stencils
     m=1
     OI(:) = 0.d0

     DO j=1,nne(i)
        el1 = indel(j,i)
        IF(el1<1) CYCLE

        IF (j<=nne(i)-2) THEN
          el2 = indel(j+1,i)
          el3 = indel(j+2,i)

        ELSEIF (j==nne(i)-1) THEN
          el2 = indel(j+1,i)
          el3 = indel(1,i)

        ELSEIF (j==nne(i)) THEN
          el2 = indel(1,i)
          el3 = indel(2,i)
        ENDIF

        IF((el2<1).or.(el3<1)) CYCLE

        CALL lin_polyn(xctr(el1),yctr(el1),xctr(el2),yctr(el2),&
            &xctr(el3),yctr(el3),qdt_e(el1,:),qdt_e(el2,:),  &
            &qdt_e(el3,:),Cx(i,m,:),Cy(i,m,:))

        !Compute oscillation indicator
        OIx = DSQRT(vc_area(i)*(Cx(i,m,1)**2.d0+Cx(i,m,2)**2.d0)&
                   &/((area(el1)+area(el2)+area(el3))/3.d0))

        OIy = DSQRT(vc_area(i)*(Cy(i,m,1)**2.d0+Cy(i,m,2)**2.d0)&
                   &/((area(el1)+area(el2)+area(el3))/3.d0))

        OI(m) = OIx+OIy

        !Sorting oscillation indicator
        IF (m>1) THEN
          DO p=1,m-1
             IF (OI(m)<OI(p)) THEN
               OItmp = OI(p)
               Cxtmp(:) = Cx(i,p,:); Cytmp(:) = Cy(i,p,:)

               OI(p) = OI(m)
               Cx(i,p,:) = Cx(i,m,:); Cy(i,p,:) = Cy(i,m,:)

               OI(m) = OItmp
               Cx(i,m,:) = Cxtmp(:); Cy(i,m,:) = Cytmp(:)
             ENDIF
          ENDDO
        ENDIF

        m=m+1

     ENDDO!j

     w_sum = 0.d0
     IF (nne(i)<=nb_comb) THEN
       DO j=1,nne(i)
          w_sum = w_sum + 1.d0/((epsi+OI(j))**rp)
       ENDDO
     ELSE
       DO j=1,nb_comb
          w_sum = w_sum + 1.d0/((epsi+OI(j))**rp)
       ENDDO
     ENDIF

     IF (nne(i)<=nb_comb) THEN
       DO j=1,nne(i)
          w(i,j) = 1.d0/((epsi+OI(j))**rp)/w_sum
       ENDDO
     ELSE
       DO j=1,nb_comb
          w(i,j) = 1.d0/((epsi+OI(j))**rp)/w_sum
       ENDDO
     ENDIF

  ENDDO !i=1,np

  DO j=1,nb_comb
     CALL exchange_p2d(w(:,j))
     DO k=1,3
        CALL exchange_p2d(Cx(:,j,k)); CALL exchange_p2d(Cy(:,j,k))
     ENDDO
  ENDDO

  qsan = 0.d0
  DO i = 1,nea
     xi = xctr(i); yi = yctr(i)
     nd1 = elnode(1,i); nd2 = elnode(2,i); nd3 = elnode(3,i)
     xnd1 = xnd(nd1); xnd2 = xnd(nd2); xnd3 = xnd(nd3)
     ynd1 = ynd(nd1); ynd2 = ynd(nd2); ynd3 = ynd(nd3)

     !Compute total water depth at the 3 nodes
     h1 = dp(nd1)+eta2(nd1)
     h2 = dp(nd2)+eta2(nd2)
     h3 = dp(nd3)+eta2(nd3)
     htot = (h1+h2+h3)*ONETHIRD

     IF((idry_e(i)==1)) CYCLE

     !1st segment
     xp = 0.5d0*(xnd1+xnd2)
     yp = 0.5d0*(ynd1+ynd2)
     xG1 = c*xp+(1.d0-c)*xi; yG1 = c*yp+(1.d0-c)*yi
     xG2 = c*xi+(1.d0-c)*xp; yG2 = c*yi+(1.d0-c)*yp

     Px1G1 = 0.d0; Py1G1 = 0.d0; Px2G1 = 0.d0; Py2G1 = 0.d0
     Px1G2 = 0.d0; Py1G2 = 0.d0; Px2G2 = 0.d0; Py2G2 = 0.d0
     Px1nd = 0.d0; Py1nd = 0.d0; Px2nd = 0.d0; Py2nd = 0.d0
     IF (nne(nd1)>3) THEN
       DO j=1,nb_comb
          Px1G1 = Px1G1 + w(nd1,j)*(Cx(nd1,j,1)*xG1+Cx(nd1,j,2)*yG1+  &
                                    Cx(nd1,j,3))
          Py1G1 = Py1G1 + w(nd1,j)*(Cy(nd1,j,1)*xG1+Cy(nd1,j,2)*yG1+  &
                                    Cy(nd1,j,3))
          Px1G2 = Px1G2 + w(nd1,j)*(Cx(nd1,j,1)*xG2+Cx(nd1,j,2)*yG2+  &
                                    Cx(nd1,j,3))
          Py1G2 = Py1G2 + w(nd1,j)*(Cy(nd1,j,1)*xG2+Cy(nd1,j,2)*yG2+  &
                                    Cy(nd1,j,3))
          Px1nd = Px1nd + w(nd1,j)*(Cx(nd1,j,1)*xnd1+Cx(nd1,j,2)*ynd1+&
                                    Cx(nd1,j,3))
          Py1nd = Py1nd + w(nd1,j)*(Cy(nd1,j,1)*xnd1+Cy(nd1,j,2)*ynd1+&
                                    Cy(nd1,j,3))
       ENDDO
     ELSEIF (nne(nd1)==3) THEN
       Px1G1 = Cx(nd1,1,1)*xG1+Cx(nd1,1,2)*yG1+Cx(nd1,1,3)
       Py1G1 = Cy(nd1,1,1)*xG1+Cy(nd1,1,2)*yG1+Cy(nd1,1,3)
       Px1G2 = Cx(nd1,1,1)*xG2+Cx(nd1,1,2)*yG2+Cx(nd1,1,3)
       Py1G2 = Cy(nd1,1,1)*xG2+Cy(nd1,1,2)*yG2+Cy(nd1,1,3)
       Px1nd = Cx(nd1,1,1)*xnd1+Cx(nd1,1,2)*ynd1+Cx(nd1,1,3)
       Py1nd = Cy(nd1,1,1)*xnd1+Cy(nd1,1,2)*ynd1+Cy(nd1,1,3)
     ELSEIF (nne(nd1)<3) THEN
       Px1G1 = qdt_e(i,1); Py1G1 = qdt_e(i,2)
       Px1G2 = qdt_e(i,1); Py1G2 = qdt_e(i,2)
       Px1nd = qdt_e(i,1); Py1nd = qdt_e(i,2)
     ENDIF

     IF (nne(nd2)>3) THEN
       DO j=1,nb_comb
          Px2G1 = Px2G1 + w(nd2,j)*(Cx(nd2,j,1)*xG1+Cx(nd2,j,2)*yG1+  &
                                    Cx(nd2,j,3))
          Py2G1 = Py2G1 + w(nd2,j)*(Cy(nd2,j,1)*xG1+Cy(nd2,j,2)*yG1+  &
                                    Cy(nd2,j,3))
          Px2G2 = Px2G2 + w(nd2,j)*(Cx(nd2,j,1)*xG2+Cx(nd2,j,2)*yG2+  &
                                    Cx(nd2,j,3))
          Py2G2 = Py2G2 + w(nd2,j)*(Cy(nd2,j,1)*xG2+Cy(nd2,j,2)*yG2+  &
                                    Cy(nd2,j,3))
          Px2nd = Px2nd + w(nd2,j)*(Cx(nd2,j,1)*xnd2+Cx(nd2,j,2)*ynd2+&
                                    Cx(nd2,j,3))
          Py2nd = Py2nd + w(nd2,j)*(Cy(nd2,j,1)*xnd2+Cy(nd2,j,2)*ynd2+&
                                    Cy(nd2,j,3))
       ENDDO
     ELSEIF (nne(nd2)==3) THEN
       Px2G1 = Cx(nd2,1,1)*xG1+Cx(nd2,1,2)*yG1+Cx(nd2,1,3)
       Py2G1 = Cy(nd2,1,1)*xG1+Cy(nd2,1,2)*yG1+Cy(nd2,1,3)
       Px2G2 = Cx(nd2,1,1)*xG2+Cx(nd2,1,2)*yG2+Cx(nd2,1,3)
       Py2G2 = Cy(nd2,1,1)*xG2+Cy(nd2,1,2)*yG2+Cy(nd2,1,3)
       Px2nd = Cx(nd2,1,1)*xnd2+Cx(nd2,1,2)*ynd2+Cx(nd2,1,3)
       Py2nd = Cy(nd2,1,1)*xnd2+Cy(nd2,1,2)*ynd2+Cy(nd2,1,3)
     ELSEIF (nne(nd2)<3) THEN
       Px2G1 = qdt_e(i,1); Py2G1 = qdt_e(i,2)
       Px2G2 = qdt_e(i,1); Py2G2 = qdt_e(i,2)
       Px2nd = qdt_e(i,1); Py2nd = qdt_e(i,2)
     ENDIF

     P1G1 = Px1G1*(yi-yp)-Py1G1*(xi-xp)
     P1G2 = Px1G2*(yi-yp)-Py1G2*(xi-xp)
     P1 = 0.5d0*(P1G1+P1G2)
     P2G1 = Px2G1*(yi-yp)-Py2G1*(xi-xp)
     P2G2 = Px2G2*(yi-yp)-Py2G2*(xi-xp)
     P2 = 0.5d0*(P2G1+P2G2)
     P1nd = Px1nd*(yi-yp)-Py1nd*(xi-xp)
     P2nd = Px2nd*(yi-yp)-Py2nd*(xi-xp)

     !compute the r_flux parameter (quantifying the local upwinding applied)
     IF (0.5d0*(h1+h2)>0.d0) THEN
       r_flux = abs(h2-h1)/(0.5d0*(h1+h2))
     ELSE
       r_flux = 0.d0
     ENDIF

     !compute weight parameter for numerical flux
     alpha = MAX(0.d0,MIN(r_flux,beta_comp))

     !compute the numerical flux
     IF (h2.LE.h1) THEN
       IF (P1nd<=P2nd) THEN
         flux = P1+0.5d0*alpha*(P1nd-P1)
       ELSE
         flux = P2+0.5d0*alpha*(P2nd-P2)
       ENDIF
     ELSE
       IF (P1nd>=P2nd) THEN
         flux = P1+0.5d0*alpha*(P1nd-P1)
       ELSE
         flux = P2+0.5d0*alpha*(P2nd-P2)
       ENDIF
     ENDIF

     qsan(nd1) = qsan(nd1) + flux
     qsan(nd2) = qsan(nd2) - flux

     !2nd segment
     xp = 0.5d0*(xnd(nd2)+xnd(nd3))
     yp = 0.5d0*(ynd(nd2)+ynd(nd3))
     xG1 = c*xp+(1.d0-c)*xi; yG1 = c*yp+(1.d0-c)*yi
     xG2 = c*xi+(1.d0-c)*xp; yG2 = c*yi+(1.d0-c)*yp

     Px1G1 = 0.d0; Py1G1 = 0.d0; Px2G1 = 0.d0; Py2G1 = 0.d0
     Px1G2 = 0.d0; Py1G2 = 0.d0; Px2G2 = 0.d0; Py2G2 = 0.d0
     Px1nd = 0.d0; Py1nd = 0.d0; Px2nd = 0.d0; Py2nd = 0.d0
     IF (nne(nd2)>3) THEN
       DO j=1,nb_comb
          Px1G1 = Px1G1 + w(nd2,j)*(Cx(nd2,j,1)*xG1+Cx(nd2,j,2)*yG1+  &
                                    Cx(nd2,j,3))
          Py1G1 = Py1G1 + w(nd2,j)*(Cy(nd2,j,1)*xG1+Cy(nd2,j,2)*yG1+  &
                                    Cy(nd2,j,3))
          Px1G2 = Px1G2 + w(nd2,j)*(Cx(nd2,j,1)*xG2+Cx(nd2,j,2)*yG2+  &
                                    Cx(nd2,j,3))
          Py1G2 = Py1G2 + w(nd2,j)*(Cy(nd2,j,1)*xG2+Cy(nd2,j,2)*yG2+  &
                                    Cy(nd2,j,3))
          Px1nd = Px1nd + w(nd2,j)*(Cx(nd2,j,1)*xnd2+Cx(nd2,j,2)*ynd2+&
                                    Cx(nd2,j,3))
          Py1nd = Py1nd + w(nd2,j)*(Cy(nd2,j,1)*xnd2+Cy(nd2,j,2)*ynd2+&
                                    Cy(nd2,j,3))
       ENDDO
     ELSEIF (nne(nd2)==3) THEN
       Px1G1 = Cx(nd2,1,1)*xG1+Cx(nd2,1,2)*yG1+Cx(nd2,1,3)
       Py1G1 = Cy(nd2,1,1)*xG1+Cy(nd2,1,2)*yG1+Cy(nd2,1,3)
       Px1G2 = Cx(nd2,1,1)*xG2+Cx(nd2,1,2)*yG2+Cx(nd2,1,3)
       Py1G2 = Cy(nd2,1,1)*xG2+Cy(nd2,1,2)*yG2+Cy(nd2,1,3)
       Px1nd = Cx(nd2,1,1)*xnd2+Cx(nd2,1,2)*ynd2+Cx(nd2,1,3)
       Py1nd = Cy(nd2,1,1)*xnd2+Cy(nd2,1,2)*ynd2+Cy(nd2,1,3)
     ELSEIF (nne(nd2)<3) THEN
       Px1G1 = qdt_e(i,1); Py1G1 = qdt_e(i,2)
       Px1G2 = qdt_e(i,1); Py1G2 = qdt_e(i,2)
       Px1nd = qdt_e(i,1); Py1nd = qdt_e(i,2)
     ENDIF

     IF (nne(nd3)>3) THEN
       DO j=1,nb_comb
          Px2G1 = Px2G1 + w(nd3,j)*(Cx(nd3,j,1)*xG1+Cx(nd3,j,2)*yG1+  &
                                     Cx(nd3,j,3))
          Py2G1 = Py2G1 + w(nd3,j)*(Cy(nd3,j,1)*xG1+Cy(nd3,j,2)*yG1+  &
                                    Cy(nd3,j,3))
          Px2G2 = Px2G2 + w(nd3,j)*(Cx(nd3,j,1)*xG2+Cx(nd3,j,2)*yG2+  &
                                    Cx(nd3,j,3))
          Py2G2 = Py2G2 + w(nd3,j)*(Cy(nd3,j,1)*xG2+Cy(nd3,j,2)*yG2+  &
                                    Cy(nd3,j,3))
          Px2nd = Px2nd + w(nd3,j)*(Cx(nd3,j,1)*xnd3+Cx(nd3,j,2)*ynd3+&
                                    Cx(nd3,j,3))
          Py2nd = Py2nd + w(nd3,j)*(Cy(nd3,j,1)*xnd3+Cy(nd3,j,2)*ynd3+&
                                    Cy(nd3,j,3))
       ENDDO
     ELSEIF (nne(nd3)==3) THEN
       Px2G1 = Cx(nd3,1,1)*xG1+Cx(nd3,1,2)*yG1+Cx(nd3,1,3)
       Py2G1 = Cy(nd3,1,1)*xG1+Cy(nd3,1,2)*yG1+Cy(nd3,1,3)
       Px2G2 = Cx(nd3,1,1)*xG2+Cx(nd3,1,2)*yG2+Cx(nd3,1,3)
       Py2G2 = Cy(nd3,1,1)*xG2+Cy(nd3,1,2)*yG2+Cy(nd3,1,3)
       Px2nd = Cx(nd3,1,1)*xnd3+Cx(nd3,1,2)*ynd3+Cx(nd3,1,3)
       Py2nd = Cy(nd3,1,1)*xnd3+Cy(nd3,1,2)*ynd3+Cy(nd3,1,3)
     ELSEIF (nne(nd3)<3) THEN
       Px2G1 = qdt_e(i,1); Py2G1 = qdt_e(i,2)
       Px2G2 = qdt_e(i,1); Py2G2 = qdt_e(i,2)
       Px2nd = qdt_e(i,1); Py2nd = qdt_e(i,2)
     ENDIF

     P1G1 = Px1G1*(yi-yp)-Py1G1*(xi-xp)
     P1G2 = Px1G2*(yi-yp)-Py1G2*(xi-xp)
     P1 = 0.5d0*(P1G1+P1G2)
     P2G1 = Px2G1*(yi-yp)-Py2G1*(xi-xp)
     P2G2 = Px2G2*(yi-yp)-Py2G2*(xi-xp)
     P2 = 0.5d0*(P2G1+P2G2)
     P1nd = Px1nd*(yi-yp)-Py1nd*(xi-xp)
     P2nd = Px2nd*(yi-yp)-Py2nd*(xi-xp)

     !compute the r_flux parameter (quantifying the local upwinding applied)
     IF (0.5d0*(h2+h3)>0.d0) THEN
       r_flux = abs(h3-h2)/(0.5d0*(h2+h3))
     ELSE
       r_flux = 0.d0
     ENDIF

     !compute weight parameter for numerical flux
     alpha = MAX(0.d0,MIN(r_flux,beta_comp))

     !compute the numerical flux
     IF (h3.LE.h2) THEN
       IF (P1nd<=P2nd) THEN
         flux = P1+0.5d0*alpha*(P1nd-P1)
       ELSE
         flux = P2+0.5d0*alpha*(P2nd-P2)
       ENDIF
     ELSE
       IF (P1nd>=P2nd) THEN
         flux = P1+0.5d0*alpha*(P1nd-P1)
       ELSE
         flux = P2+0.5d0*alpha*(P2nd-P2)
       ENDIF
     ENDIF

     qsan(nd2) = qsan(nd2) + flux
     qsan(nd3) = qsan(nd3) - flux

     !3rd segment
     xp = 0.5d0*(xnd(nd3)+xnd(nd1))
     yp = 0.5d0*(ynd(nd3)+ynd(nd1))
     xG1 = c*xp+(1.d0-c)*xi; yG1 = c*yp+(1.d0-c)*yi
     xG2 = c*xi+(1.d0-c)*xp; yG2 = c*yi+(1.d0-c)*yp

     Px1G1 = 0.d0; Py1G1 = 0.d0; Px2G1 = 0.d0; Py2G1 = 0.d0
     Px1G2 = 0.d0; Py1G2 = 0.d0; Px2G2 = 0.d0; Py2G2 = 0.d0
     Px1nd = 0.d0; Py1nd = 0.d0; Px2nd = 0.d0; Py2nd = 0.d0
     IF (nne(nd3)>3) THEN
       DO j=1,nb_comb
          Px1G1 = Px1G1 + w(nd3,j)*(Cx(nd3,j,1)*xG1+Cx(nd3,j,2)*yG1+  &
                                    Cx(nd3,j,3))
          Py1G1 = Py1G1 + w(nd3,j)*(Cy(nd3,j,1)*xG1+Cy(nd3,j,2)*yG1+  &
                                    Cy(nd3,j,3))
          Px1G2 = Px1G2 + w(nd3,j)*(Cx(nd3,j,1)*xG2+Cx(nd3,j,2)*yG2+  &
                                    Cx(nd3,j,3))
          Py1G2 = Py1G2 + w(nd3,j)*(Cy(nd3,j,1)*xG2+Cy(nd3,j,2)*yG2+  &
                                    Cy(nd3,j,3))
          Px1nd = Px1nd + w(nd3,j)*(Cx(nd3,j,1)*xnd3+Cx(nd3,j,2)*ynd3+&
                                    Cx(nd3,j,3))
          Py1nd = Py1nd + w(nd3,j)*(Cy(nd3,j,1)*xnd3+Cy(nd3,j,2)*ynd3+&
                                    Cy(nd3,j,3))
       ENDDO
     ELSEIF (nne(nd3)==3) THEN
       Px1G1 = Cx(nd3,1,1)*xG1+Cx(nd3,1,2)*yG1+Cx(nd3,1,3)
       Py1G1 = Cy(nd3,1,1)*xG1+Cy(nd3,1,2)*yG1+Cy(nd3,1,3)
       Px1G2 = Cx(nd3,1,1)*xG2+Cx(nd3,1,2)*yG2+Cx(nd3,1,3)
       Py1G2 = Cy(nd3,1,1)*xG2+Cy(nd3,1,2)*yG2+Cy(nd3,1,3)
       Px1nd = Cx(nd3,1,1)*xnd3+Cx(nd3,1,2)*ynd3+Cx(nd3,1,3)
       Py1nd = Cy(nd3,1,1)*xnd3+Cy(nd3,1,2)*ynd3+Cy(nd3,1,3)
     ELSEIF (nne(nd3)<3) THEN
       Px1G1 = qdt_e(i,1); Py1G1 = qdt_e(i,2)
       Px1G2 = qdt_e(i,1); Py1G2 = qdt_e(i,2)
       Px1nd = qdt_e(i,1); Py1nd = qdt_e(i,2)
     ENDIF

     IF (nne(nd1)>3) THEN
       DO j=1,nb_comb
          Px2G1 = Px2G1 + w(nd1,j)*(Cx(nd1,j,1)*xG1+Cx(nd1,j,2)*yG1+  &
                                    Cx(nd1,j,3))
          Py2G1 = Py2G1 + w(nd1,j)*(Cy(nd1,j,1)*xG1+Cy(nd1,j,2)*yG1+  &
                                    Cy(nd1,j,3))
          Px2G2 = Px2G2 + w(nd1,j)*(Cx(nd1,j,1)*xG2+Cx(nd1,j,2)*yG2+  &
                                    Cx(nd1,j,3))
          Py2G2 = Py2G2 + w(nd1,j)*(Cy(nd1,j,1)*xG2+Cy(nd1,j,2)*yG2+  &
                                    Cy(nd1,j,3))
          Px2nd = Px2nd + w(nd1,j)*(Cx(nd1,j,1)*xnd1+Cx(nd1,j,2)*ynd1+&
                                    Cx(nd1,j,3))
          Py2nd = Py2nd + w(nd1,j)*(Cy(nd1,j,1)*xnd1+Cy(nd1,j,2)*ynd1+&
                                    Cy(nd1,j,3))
       ENDDO
     ELSEIF (nne(nd1)==3) THEN
       Px2G1 = Cx(nd1,1,1)*xG1+Cx(nd1,1,2)*yG1+Cx(nd1,1,3)
       Py2G1 = Cy(nd1,1,1)*xG1+Cy(nd1,1,2)*yG1+Cy(nd1,1,3)
       Px2G2 = Cx(nd1,1,1)*xG2+Cx(nd1,1,2)*yG2+Cx(nd1,1,3)
       Py2G2 = Cy(nd1,1,1)*xG2+Cy(nd1,1,2)*yG2+Cy(nd1,1,3)
       Px2nd = Cx(nd1,1,1)*xnd1+Cx(nd1,1,2)*ynd1+Cx(nd1,1,3)
       Py2nd = Cy(nd1,1,1)*xnd1+Cy(nd1,1,2)*ynd1+Cy(nd1,1,3)
     ELSEIF (nne(nd1)<3) THEN
       Px2G1 = qdt_e(i,1); Py2G1 = qdt_e(i,2)
       Px2G2 = qdt_e(i,1); Py2G2 = qdt_e(i,2)
       Px2nd = qdt_e(i,1); Py2nd = qdt_e(i,2)
     ENDIF

     P1G1 = Px1G1*(yi-yp)-Py1G1*(xi-xp)
     P1G2 = Px1G2*(yi-yp)-Py1G2*(xi-xp)
     P1 = 0.5d0*(P1G1+P1G2)
     P2G1 = Px2G1*(yi-yp)-Py2G1*(xi-xp)
     P2G2 = Px2G2*(yi-yp)-Py2G2*(xi-xp)
     P2 = 0.5d0*(P2G1+P2G2)
     P1nd = Px1nd*(yi-yp)-Py1nd*(xi-xp)
     P2nd = Px2nd*(yi-yp)-Py2nd*(xi-xp)

     !compute the r_flux parameter (quantifying the local upwinding applied)
     IF (0.5d0*(h3+h1)>0.d0) THEN
       r_flux = abs(h1-h3)/(0.5d0*(h3+h1))
     ELSE
       r_flux = 0.d0
     ENDIF

     !compute weight parameter for numerical flux
     alpha = MAX(0.d0,MIN(r_flux,beta_comp))

     !compute the numerical flux
     IF (h1.LE.h3) THEN
       IF (P1nd<=P2nd) THEN
         flux = P1+0.5d0*alpha*(P1nd-P1)
       ELSE
         flux = P2+0.5d0*alpha*(P2nd-P2)
       ENDIF
     ELSE
       IF (P1nd>=P2nd) THEN
         flux = P1+0.5d0*alpha*(P1nd-P1)
       ELSE
         flux = P2+0.5d0*alpha*(P2nd-P2)
       ENDIF
     ENDIF

     qsan(nd3) = qsan(nd3) + flux
     qsan(nd1) = qsan(nd1) - flux

  ENDDO !nea

  CALL exchange_p2d(qsan)

  END SUBROUTINE weno_scheme

  SUBROUTINE lin_polyn(x1,y1,x2,y2,x3,y3,q1,q2,q3,Cx,Cy)
!--------------------------------------------------------------------
! Compute coefficients for bilinear interpolation of vectors
!--------------------------------------------------------------------
  USE schism_glbl, only : rkind
  USE schism_msgp, ONLY : parallel_abort

  IMPLICIT NONE

!- Arguments --------------------------------------------------------
  REAL(rkind), INTENT(IN) :: x1,y1,x2,y2,x3,y3
  REAL(rkind), DIMENSION(2), INTENT(IN) :: q1,q2,q3
  REAL(rkind), DIMENSION(3), INTENT(OUT) :: Cx,Cy
!- Local variables --------------------------------------------------
  REAL(rkind) :: detA
  REAL(rkind), DIMENSION(2) :: det1,det2,det3
!--------------------------------------------------------------------

  CALL det_3x3(x1,y1,1.d0,x2,y2,1.d0,x3,y3,1.d0,detA)

  CALL det_3x3(q1(1),y1,1.d0,q2(1),y2,1.d0,q3(1),y3,1.d0,det1(1))
  CALL det_3x3(q1(2),y1,1.d0,q2(2),y2,1.d0,q3(2),y3,1.d0,det1(2))

  CALL det_3x3(x1,q1(1),1.d0,x2,q2(1),1.d0,x3,q3(1),1.d0,det2(1))
  CALL det_3x3(x1,q1(2),1.d0,x2,q2(2),1.d0,x3,q3(2),1.d0,det2(2))

  CALL det_3x3(x1,y1,q1(1),x2,y2,q2(1),x3,y3,q3(1),det3(1))
  CALL det_3x3(x1,y1,q1(2),x2,y2,q2(2),x3,y3,q3(2),det3(2))

  IF (detA /= 0.d0) THEN
    Cx(1) = det1(1)/detA
    Cy(1) = det1(2)/detA
    Cx(2) = det2(1)/detA
    Cy(2) = det2(2)/detA
    Cx(3) = det3(1)/detA
    Cy(3) = det3(2)/detA
  ELSE
    CALL parallel_abort('detA=0 in lin_polyn subroutine')
  ENDIF

  END SUBROUTINE lin_polyn

  SUBROUTINE det_3x3(a,b,c,d,e,f,g,h,i,det)
!--------------------------------------------------------------------
! Compute determinant of a 3x3 matrice
!--------------------------------------------------------------------
  USE schism_glbl, only : rkind

  IMPLICIT NONE

!- Arguments --------------------------------------------------------
  REAL(rkind), INTENT(IN) :: a,b,c,d,e,f,g,h,i
  REAL(rkind), INTENT(OUT) :: det
!--------------------------------------------------------------------

  det = a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)

  END SUBROUTINE det_3x3


!--------------------------------------------------------------------
      SUBROUTINE acceleration_bedload_hoe2003(inea,wave_per,&
                                            wave_dir,Uc,Ut,Uorbi,&
                                            Qaccx,Qaccy)

!--------------------------------------------------------------------
! This subroutine computes bedload transport caused by wave
! acceleration, Qacc, following Hoefel and Elgar, Science, 2003.
! Need the reconstitution of Uorbi(t) thanks to the theory of
! Elfrink et al. (2006), i.e. iasym=1
!
! Author: baptiste mengual (baptiste.mengual@univ-lr.fr)
! Date: 18/03/2019
!
!--------------------------------------------------------------------

      USE sed_mod
      USE schism_glbl, ONLY : pi,ielg,errmsg,rkind,dt
      USE schism_msgp, ONLY : parallel_abort

      IMPLICIT NONE

!- Arguments --------------------------------------------------------
      INTEGER, INTENT(IN) :: inea
      REAL(rkind), INTENT(IN) :: wave_per,wave_dir, Uc, Ut
      REAL(rkind), DIMENSION(ech_uorb+1), INTENT(IN) :: Uorbi
      REAL(rkind), INTENT(OUT) :: Qaccx,Qaccy
!- Local variables --------------------------------------------------
      INTEGER :: t
      REAL(rkind), PARAMETER :: nu=1.36d-6 ! Cinematic viscosity
      REAL(rkind) :: cff1,cff2,acc3,acc2,w_asym,aspike,Qacc
      REAL(rkind) :: abskb,fw,osmgd,ks,d50,tau_loc,theta_loc,&
                     ratio,dstar,theta_cr
      REAL(rkind), DIMENSION(ech_uorb) :: acc
!--------------------------------------------------------------------

      d50=bottom(inea,isd50) ! Median grain size (m)
      ks =2.5d0*d50

      cff1=wave_per/DBLE(ech_uorb)
      cff2=0.0d0

      ! Compute acceleration along a wave period
      DO t=0,ech_uorb-1
        cff2=cff2+1.0d0
        acc(t)=(Uorbi(t+1)-Uorbi(t))/cff1
      END DO

      acc3=sum(acc(:)**3.0d0)/cff2
      acc2=sum(acc(:)**2.0d0)/cff2

      ! Compute wave asymmetry coefficient
      IF ( Uc+Ut > 0.0d0) THEN
        w_asym =MAX(Uc/((Uc+Ut)/2.0d0)-1.0d0,0.0d0)
      ELSE
        w_asym = 0.d0
      ENDIF

      ! Compute aspike (impulses that transfer
      ! momentum to the near-bed fluid and 
      ! sediment)
      IF (    acc2 .GT. 0.0d0 .AND. &
          & w_asym .GT. 0.0d0 .AND. &
          & w_asym .LT. w_asym_max) THEN
        ! Only considered if the average acceleration
        ! is higher than 0 and waves are not too 
        ! asymmetric (w_asym_max criterion)
        aspike=acc3/acc2
      ELSE
        aspike=0.0d0
      END IF


      ! Compute acceleration-skewness induced transport
      ! Qacc: [m2.s-1]
      ! kacc_hoe: [m.s]
      ! aspike: [m.s-2]
      IF (thresh_acc_opt==0) THEN
        !!! No condition
        Qacc=kacc_hoe*aspike
      ELSE IF (thresh_acc_opt==1) THEN
        !!! Criterion based on Uc
        ! Compute a mobility parameter based on U_crest
        ! -> theta_loc
        abskb = Uc*wave_per/2.0d0/pi/ks
        IF (abskb.LE.1.57d0) THEN
          fw = 0.3d0
        ELSE
          fw = 0.00251d0*EXP(5.21*abskb**(-0.19d0))
        ENDIF
        tau_loc = 0.5d0*fw*Uc*Uc ! m2.s-2
        osmgd = 1.0d0/smgd
        theta_loc = tau_loc*osmgd
        ! Compute a Shields mobility criterion
        ! -> theta_cr
        ratio=bottom(inea,idens)/rhom
        dstar=bottom(inea,isd50)*(g*(ratio-1.d0)/nu**2.d0)**(1.d0/3d0)
        if(1.d0+1.2d0*dstar==0) call parallel_abort('hoe2003, div. by 0')
        theta_cr = 0.3d0/(1.d0+1.2d0*dstar) + 0.055d0 *                &
        &         (1.d0-EXP(-0.02d0*dstar))
        ! Qacc considered only if theta_loc > theta_cr
        IF (theta_loc .GT. theta_cr) THEN
          Qacc=kacc_hoe*aspike
        ELSE
          Qacc=0.0d0
        END IF
      ELSE IF (thresh_acc_opt==2) THEN
        !!! Critical acceleration acrit
        IF (DABS(aspike) .GE. acrit) THEN
          Qacc=kacc_hoe*(aspike-SIGN(1.0d0,aspike)*acrit)
        ELSE
          Qacc=0.0d0
        END IF
      ELSE
        call parallel_abort('hoe2003, invalid thresh_acc_opt')
      END IF

      Qacc=MAX(Qacc,0.0d0) ! m2/s

      ! Projection according to wave direction
      Qaccx=Qacc*DCOS(wave_dir)*dt !m2
      Qaccy=Qacc*DSIN(wave_dir)*dt

      ! Check
      IF(Qacc/=Qacc) THEN
        WRITE(errmsg,*)'Qacc is NaN : inea,wave_per,wave_dir',&
                    & 'Uc,Ut,w_asym,acc3,acc2,aspike,Qacc,Qaccx,Qaccy :',&
                    & inea,wave_per,wave_dir, &
                    & Uc,Ut,w_asym,acc3,acc2,aspike,Qacc,Qaccx,Qaccy
      ENDIF


      END SUBROUTINE acceleration_bedload_hoe2003
!--------------------------------------------------------------------



!--------------------------------------------------------------------
      SUBROUTINE acceleration_bedload_dub2015(inea,H,wave_per,&
                                        & wave_dir,depth,Qaccx,Qaccy)

!--------------------------------------------------------------------
! This subroutine computes bedload transport caused by wave
! acceleration, Qacc, following Dubarbier et al. (2015, Coastal 
! Engineering).
! Velocity asymmetry Au is computed according to Ruessink et al. (2012,
! Coastal Engineering).
!
! Author: baptiste mengual (baptiste.mengual@univ-lr.fr)
! Date: 26/03/2019
!
!--------------------------------------------------------------------

      USE sed_mod
      USE schism_glbl, ONLY : pi,ielg,errmsg,rkind,dt
      USE schism_msgp, ONLY : parallel_abort

      IMPLICIT NONE

!- Arguments --------------------------------------------------------
      INTEGER, INTENT(IN) :: inea
      REAL(rkind), INTENT(IN) :: H,wave_per,wave_dir,depth
      REAL(rkind), INTENT(OUT) :: Qaccx,Qaccy
!- Local variables --------------------------------------------------
      REAL(rkind), PARAMETER :: p1=0.0d0
      REAL(rkind), PARAMETER :: p2=0.857d0
      REAL(rkind), PARAMETER :: p3=-0.471d0
      REAL(rkind), PARAMETER :: p4=0.297d0
      REAL(rkind), PARAMETER :: p5=0.815
      REAL(rkind), PARAMETER :: p6=0.672
      REAL(rkind), PARAMETER :: nu=1.36d-6 ! Cinematic viscosity
      REAL(rkind) :: Hs2,L,k,Ur,psi,B,Au,w,Hrms,Uw,Aw,&
                     ratio,dstar,theta_cr,Qacc
      REAL(rkind) :: abskb,fw,ks,d50,osmgd,tau_loc,theta_loc
!--------------------------------------------------------------------

      ! Compute velocity asymmetry Au (Ruessink et al., 2012)
      Hs2=0.5d0*H
      L=wave_per*DSQRT(g*depth)
      k=2.0d0*pi/L
      Ur=0.75d0*Hs2*k/((k*depth)**3.0d0)
      psi=-pi/2.0d0 + ((pi/2.0d0)*DTANH(p5/(Ur**p6)))
      B=p1+((p2-p1)/(1.0d0+EXP((p3-LOG(Ur))/p4)))
      Au=B*DSIN(psi)

      ! Compute near-bed acceleration amplitude
      w=2.0d0*pi/wave_per
      Hrms=H/DSQRT(2.0d0)
      Uw=(pi*Hrms)/(wave_per*DSINH(k*depth))
      Aw=w*Uw

      ! Compute acceleration-skewness induced transport
      ! [m2.s-1]=[m.s]*[-]*[m.s-2]
      ! Then multiplied by dt -> [m2]

      IF (thresh_acc_opt==0) THEN
        !! No condition
        Qacc=DABS(-kacc_dub*Au*Aw)
      ELSE IF (thresh_acc_opt==1) THEN
        !!! Criterion based on Uc
        ! Compute a mobility parameter based on Uw
        ! -> theta_loc
        abskb  = 0.0d0
        fw     = 0.0d0
        d50 = bottom(inea,isd50)      ! Median grain size (m)
        ks  = 2.5d0*d50
        abskb = Uw*wave_per/2.0d0/pi/ks ! * Ratio Ab/kb
        IF (abskb.LE.1.57d0) THEN
          fw = 0.3d0
        ELSE
          fw = 0.00251d0*EXP(5.21*abskb**(-0.19d0))
        ENDIF !
        tau_loc = 0.5d0*fw*Uw*Uw ! m2.s-2
        osmgd = 1.0d0/smgd
        theta_loc = tau_loc*osmgd
        ! Compute a Shields mobility criterion
        ! -> theta_cr
        ratio=bottom(inea,idens)/rhom
        dstar=bottom(inea,isd50)*(g*(ratio-1.d0)/nu**2.d0)**(1.d0/3d0)
        if(1.d0+1.2d0*dstar==0) call parallel_abort('dub2015, div. by 0')
        theta_cr = 0.3d0/(1.d0+1.2d0*dstar) + 0.055d0 *                &
        &         (1.d0-EXP(-0.02d0*dstar))
        ! Qacc considered only if theta_loc > theta_cr
        IF (theta_loc .GT. theta_cr) THEN
          Qacc=DABS(-kacc_dub*Au*Aw)
        ELSE
          Qacc=0.0d0
        END IF
      ELSE IF (thresh_acc_opt==2) THEN
        !!! Critical acceleration acrit
        IF (Aw .GT. acrit) THEN
          Qacc=DABS(-kacc_dub*Au*Aw)
        ELSE
          Qacc=0.0d0
        END IF
      ELSE
        call parallel_abort('dub2015, invalid thresh_acc_opt')
      END IF

      Qacc=MAX(Qacc,0.0d0) ! m2/s

      ! Projection according to wave direction
      Qaccx=Qacc*DCOS(wave_dir)*dt !m2
      Qaccy=Qacc*DSIN(wave_dir)*dt


      ! Check
      IF(Qacc/=Qacc) THEN
        WRITE(errmsg,*)'Qacc is NaN : inea,H,wave_per',&
                      'wave_dir,depth,L,k,Ur,psi,B,Au,w,Hrms',    &
                      'Uw,Aw,Qacc,Qaccx,Qaccy',                   &
                      inea,&
                    & H,wave_per,wave_dir,depth,L,    &
                    & k,Ur,psi,B,Au,w,Hrms,Uw,Aw,Qacc,Qaccx,Qaccy
      ENDIF


      END SUBROUTINE acceleration_bedload_dub2015
!--------------------------------------------------------------------


!--------------------------------------------------------------------
      SUBROUTINE bedslope_effects_lesser2004(inea,ised,Fx,Fy)

      USE sed_mod
      USE schism_glbl, ONLY : pi,ielg,errmsg,rkind
      USE schism_msgp, ONLY : parallel_abort

      IMPLICIT NONE

!- Arguments --------------------------------------------------------
      INTEGER,     INTENT(IN)    :: inea,ised
      REAL(rkind), INTENT(INOUT) :: Fx,Fy
!- Local variables --------------------------------------------------
      REAL(rkind) :: a_slopex, a_slopey, beta, &
                     cff,cff1,cff2
!--------------------------------------------------------------------

      !---------------------------------------------------------------------
      ! longitudinal bed slope \beta_s
      ! limit slope to 0.9*(sed_angle) (repose angle)
      !---------------------------------------------------------------------

      cff  = dzdx*angleu+dzdy*anglev !tan(\beta_s)
      cff1 = MIN(ABS(cff),0.9d0*sed_angle)*SIGN(1.0d0,cff)
      cff2 = ATAN(cff1) !\beta_s
      if(COS(cff2)==0.or.sed_angle-cff1==0) then
        CALL parallel_abort('SED3D error in bedslope_effects_lesser2004')
      endif
      !\alpha_s
      a_slopex = 1.0d0+alpha_bs*((sed_angle/(COS(cff2)*(sed_angle-cff1)))-1.0d0)

      !---------------------------------------------------------------------
      ! - Modify bedload transport due to longitudinal bed slope
      !---------------------------------------------------------------------

      Fx = Fx*a_slopex
      Fy = Fy*a_slopex

      !---------------------------------------------------------------------
      ! - Transverse bed slope
      !---------------------------------------------------------------------

      cff = -dzdx*anglev+dzdy*angleu !tan(\beta_n)
      cff1 = ABS(bustr(inea))+ABS(bvstr(inea))

      ! - Test used to prevent Inf & NaNs with very small bustr/bvstr
      !   May still produce unrealistic values
      IF(cff1<1d-10) THEN
        a_slopey = 0.d0
      ELSE
        cff2 = SQRT(tau_ce(ised)/cff1)
        a_slopey=alpha_bn*cff2*cff !\alpha_n
      ENDIF

      !---------------------------------------------------------------------
      ! - Add contribution of transverse to bed load
      !---------------------------------------------------------------------

      beta=Fx !temp. save
      Fx = Fx-Fy*a_slopey
 
      END SUBROUTINE bedslope_effects_lesser2004
!--------------------------------------------------------------------

