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
! MORSELFE BOTTOM FRICTION SUBROUTINES
!
! subroutine sed_current_stress
! subroutine sed_wavecurrent_stress
! subroutine sed_roughness (called inside schism_step)
! subroutine stress_soulsby (not used now)
! subroutine stress (not used now)
!
!=====================================================================
!=====================================================================

      SUBROUTINE sed_current_stress()
!--------------------------------------------------------------------!
! This routine computes the bottom shear stress for currents alone   !
! This subroutine is adapted from the former subroutine set_vbc(),   !
! adapted from a ROMS routine                                        !
! Copyright (c) 2002-2007 The ROMS/TOMS Group                        !
!   Licensed under a MIT/X style license                             !
!   See License_ROMS.txt                                             !
!                                                                    !
! Author: Ligia Pinto                                                !
! Date: xx/08/2007                                                   !
!                                                                    !
! History: 2012/12 - F.Ganthy : form homogenisation of sediments     !
!                               routines                             !
!          2013/03 - F.Ganthy : cleaning and some update related to  ! 
!                               implementation of wave-current bottom!
!                               stress and bedforms predictor        !
!          2013/05 - F.Ganthy : Further cleaning                     !
!          2020/02 - X.Bertin, B.Mengual : shear stress estimates    !
!                               based on the skin roughness          !
!                  - B.Mengual : > introduction of tau_option        !
!                                > tau_c limited by tau_max          !
!                                                                    !
!--------------------------------------------------------------------!

      USE sed_mod,   ONLY: Cdb_min,Cdb_max,Zob,vonKar,bustr,bvstr,   &
                           tau_c,drag_formulation,isd50,bottom, &
                           ustress,vstress,zstress,tau_max,tau_option
      USE schism_glbl, ONLY: rkind,nvrt,nea,dfv,idry_e,kbe,i34,elnode, &
     &uu2,vv2,ze,dpe,isbnd,ifltype,errmsg,rho0,eta2
      USE schism_msgp, ONLY: myrank,parallel_abort,exchange_e2d

      IMPLICIT NONE
      SAVE

!- Local variables --------------------------------------------------!

      INTEGER     :: i,j,n1,n2,n3,itmp,ibnd
      INTEGER     :: kb0,kb1,kb2
      REAL(rkind) :: cff1, cff2, cff3, cff4, cff5, hh
      REAL(rkind) :: wrk
      REAL(rkind) :: htot,d50,z0s

!- Start Statement --------------------------------------------------!

!---------------------------------------------------------------------
! - Set kinematic bottom momentum flux (m2/s2).
!---------------------------------------------------------------------

      bustr=0; bvstr=0; tau_c=0 !init.

      IF (drag_formulation == 1) THEN
      ! - Set logarithmic bottom stress.
        DO i=1,nea
          IF(idry_e(i)==1) CYCLE
          !Also remove elem. that has an open bnd node where vel. is
          !imposed, because a uniform profile is imposed there leading
          !to large stress
!          itmp=0 !flag
!          do j=1,i34(i)
!            ibnd=isbnd(1,elnode(j,i)) !global
!            if(ibnd>0) then; if(ifltype(ibnd)/=0) then !0324
!              itmp=1; exit
!            endif; endif
!          enddo !j
!          if(itmp==1) cycle

          kb0 = kbe(i)
          kb1 = kbe(i)+1
          htot = dpe(i)+sum(eta2(elnode(1:i34(i),i)))/i34(i)

          IF (tau_option .EQ. 1) THEN
            hh  = ze(kb1,i)-ze(kb0,i)
          ELSE IF (tau_option .EQ. 2) THEN
            hh  = zstress
          ELSE IF (tau_option .EQ. 3) THEN
            hh  = MAX(ze(kb1,i)-ze(kb0,i),zstress)
          END IF

          d50      = bottom(i,isd50)      ! Median grain size (m)         
          z0s    = d50/12.0d0             ! Nikuradse roughness  (m)

          IF(hh>z0s) THEN
          !IF(hh>Zob(i)) THEN
            !IF(hh<=0.d0.OR.Zob(i)<=0.d0) THEN
            IF(hh<=0.OR.z0s<=0) THEN
              CALL parallel_abort('SED-current_stress: Cd failed')
            ENDIF
            cff1 = 1.0d0/LOG(hh/z0s)
            !cff1 = 1.0d0/LOG(hh/Zob(i))
            cff2 = vonKar*vonKar*cff1*cff1
            wrk  = MIN(Cdb_max,MAX(Cdb_min,cff2))
          else
            wrk=Cdb_max
          endif !hh

          cff3 = vstress(i)
          cff4 = ustress(i)

          !Limit vel. in shallows - need better approach
          if(htot<=0.5d0) then
            cff3=min(cff3,1.d0)
            cff4=min(cff4,1.d0)
          endif

          cff5 = SQRT(cff4*cff4+cff3*cff3)
          ! Bottom stress [m^2/s/s]
          bustr(i) = wrk*cff4*cff5
          bvstr(i) = wrk*cff3*cff5
!          ELSE !approx.
!            kb2 = min(nvrt,kbe(i)+2)
!            cff3 = sum(vv2(kb2,elnode(1:i34(i),i)))/i34(i)- &
!                  &sum(vv2(kb1,elnode(1:i34(i),i)))/i34(i)
!            cff4 = sum(uu2(kb2,elnode(1:i34(i),i)))/i34(i)- &
!                  &sum(uu2(kb1,elnode(1:i34(i),i)))/i34(i)
!            cff5 = sum(dfv(kb1,elnode(1:i34(i),i))+dfv(kb2,elnode(1:i34(i),i)))/2/i34(i)
!            hh = ze(kb2,i)-ze(kb1,i)
!            IF(hh<=0) CALL parallel_abort('SED-current_stress: div. by 0')
!            !' Bottom stress
!            bustr(i) = cff5*cff4/hh
!            bvstr(i) = cff5*cff3/hh
!          ENDIF ! End test on Zob>ze

          ! Bottom stress mag. [m^2/s/s]
          tau_c(i) = SQRT(bustr(i)*bustr(i)+bvstr(i)*bvstr(i))

          ! BM: Limitation of current-induced shear stress
          tau_c(i) = MIN(tau_c(i),tau_max/rho0)

          IF (tau_c(i)/=tau_c(i)) THEN
            WRITE(errmsg,*)'Sed current bot. stress is NaN',myrank,i,&
            &              tau_c(i),bustr(i),bvstr(i),cff3,cff4
            CALL parallel_abort(errmsg)
          ENDIF

        END DO !i

!      ELSEIF (drag_formulation == 2) THEN
!      ! Set quadratic bottom stress.
!        DO i=1,nea
!          if(idry_e(i)==1) cycle
!          n1=elnode(1,i)
!          n2=elnode(2,i)
!          n3=elnode(3,i)
!          cff1=(vv2(kbe(i)+1,n1)+vv2(kbe(i)+1,n2)+vv2(kbe(i)+1,n3))/3
!          cff2=(uu2(kbe(i)+1,n1)+uu2(kbe(i)+1,n2)+uu2(kbe(i)+1,n3))/3
!          cff3=SQRT(cff2*cff2+cff1*cff1)
!!Error: rdrg2 not defined
!!LLP user defined, if defined valor must be read from ...
!         bustr(i)=rdrg2*cff2*cff3
!         bvstr(i)=rdrg2*cff1*cff3
!         ! Bottom stress norm
!         tau_c(i) = SQRT(bustr(i)*bustr(i)+bvstr(i)*bvstr(i))
!       END DO
!       else if (drag_formulation == 3) then
!!#elif defined UV_LDRAG
!!  Set linear bottom stress.
!       DO i=1,nea
!         if(idry_e(i)==1) cycle
!         n1=elnode(1,i)
!         n2=elnode(2,i)
!         n3=elnode(3,i)
!         cff1=(vv2(kbe(i)+1,n1)+vv2(kbe(i)+1,n2)+vv2(kbe(i)+1,n3)/3
!         cff2=(uu2(kbe(i)+1,n1)+uu2(kbe(i)+1,n2)+uu2(kbe(i)+1,n3)/3
!!Error: rdrg not defined
!!LLP user defined 
!         bustr(i)=rdrg*cff2
!         bvstr(i)=rdrg*cff1
!         ! Bottom stress norm
!         tau_c(i) = SQRT(bustr(i)*bustr(i)+bvstr(i)*bvstr(i))
!       END DO
!!#endif
      ENDIF ! End test on drag_formulation 

      ! Exchange "ghost" elements 
      ! * FG : Check if necessary
!      CALL exchange_e2d(bustr)
!      CALL exchange_e2d(bvstr)
!      CALL exchange_e2d(tau_c)

!---------------------------------------------------------------------
       END SUBROUTINE sed_current_stress

!=====================================================================
!=====================================================================

      SUBROUTINE sed_wavecurrent_stress()
!--------------------------------------------------------------------!
! This routine computes the bottom shear stress for waves-currents   !
! - Bottom stress for currents alone is firstly computed             !
! - Bottom stress for wave alone is then computed using friction     !
!   factor from Warner et al. (2008) following the formulation of    !
!   Madsen (1994)                                                    !
! - Wave-current mean botom stress is finally computed following     !
!   Soulsby (1997) eq. (69)                                          !
!                                                                    !
! Author: florian ganthy (fganthy@lnec.pt ; florian.ganthy@gmail.com)!
! Date: 2013/03/27                                                   !
!                                                                    !
! History:                                                           !
!          2020/02 - X.Bertin, B.Mengual : fw based on Swart (1974)  !
!                    using the skin friction                         !
!                  - B.Mengual : tau_w limited by tau_max            !
!                                                                    !
!--------------------------------------------------------------------!

      USE sed_mod,   ONLY: uorb,tp,Zob,bustr,bvstr,tau_c,      &
     &                      tau_w,tau_wc,isd50,bottom,tau_max 
      USE schism_glbl, ONLY: rkind,pi,nea,idry_e,errmsg,rho0
      USE schism_msgp, ONLY: myrank,parallel_abort,exchange_e2d

      IMPLICIT NONE
      SAVE

!- Local variables --------------------------------------------------!
      INTEGER     :: i
      REAL(rkind) :: abskb,fw,ks,d50
!- Start Statement --------------------------------------------------!

!---------------------------------------------------------------------
! - Compute bottom stress for currents alone
!---------------------------------------------------------------------

      CALL sed_current_stress()

!---------------------------------------------------------------------
! - Compute wave-currents bottom stress
!---------------------------------------------------------------------
#ifdef USE_WWM
      DO i = 1,nea
        abskb  = 0.0d0
        fw     = 0.0d0
        IF (idry_e(i).EQ.1) CYCLE

!---------------------------------------------------------------------
! - First computes wave alone bottom stress using the friction factor 
!   from Warner et al. (2008) following the formulation of Madsen 
!   (1994):
!           (1) Computes ratio of the wave-orbital excursion amplitude
!               to the bottom roughness (Ab/kb) where:
!                                        - Ab = Uorb*Tp/(2*pi)
!                                        - kb = 30.0*z0
!           (2) Computes the friction factor
!           (3) Computes the wave bottom stress
!---------------------------------------------------------------------
        ! * Ratio Ab/kb
        !if(Zob(i)==0) call parallel_abort('sed_current_stress: Zob=0')
        !abskb = uorb(i)*tp(i)/(30.0d0*Zob(i)*2.0d0*pi)

        !! * Friction factor
        !IF (abskb.LE.0.2d0) THEN
        !  fw = 0.3d0
        !ELSEIF ((abskb.GT.0.2d0).AND.(abskb.LE.100.0d0)) THEN
        !  fw = EXP(-8.82d0 + 7.02d0*abskb**(-0.078d0))
        !ELSE
        !  fw = EXP(-7.30d0 + 5.61d0*abskb**(-0.109d0))
        !ENDIF ! End test on abskb value

! XB: compute fw based on Swart (1974) using the skin friction
        d50 = bottom(i,isd50)      ! Median grain size (m)
        ks  = 2.5d0*d50
        abskb = uorb(i)*tp(i)/2.0d0/pi/ks ! * Ratio Ab/kb
        IF (abskb.LE.1.57d0) THEN
          fw = 0.3d0
        ELSE
          fw = 0.00251d0*EXP(5.21*abskb**(-0.19d0))
        ENDIF !
! end XB

        ! * Wave bottom stress (in m2.s-2)
        tau_w(i) = 0.5d0*fw*uorb(i)*uorb(i)

        !BM: same limitation than for tau_c
        tau_w(i)=MIN(tau_w(i),tau_max/rho0)

        ! * Consistency check
        IF (tau_w(i)/=tau_w(i)) THEN
          WRITE(errmsg,*)'Sed wave bot. stress is NaN',myrank,i,     &
          &              tau_w(i),fw,uorb(i)
          CALL parallel_abort(errmsg)
        ENDIF

!---------------------------------------------------------------------
! - Computes mean wave-current bottom stress (Soulsby, 1997, eq. 69)
!---------------------------------------------------------------------

        if(tau_c(i)+tau_w(i)==0) then
          tau_wc(i)=0
        else
          tau_wc(i) = tau_c(i)*(1.0d0+1.2d0*(tau_w(i)/(tau_c(i)+tau_w(i)))**3.2d0)
        endif

      ENDDO !i = 1,nea

#else 
      ! No wave model
      DO i = 1,nea
        tau_wc(i)  = tau_c(i)
      ENDDO ! End loop nea
#endif
      ! Exchange "ghost" elements 
!      CALL exchange_e2d(tau_c) ! * FG : Check if necessary

!---------------------------------------------------------------------
      END SUBROUTINE sed_wavecurrent_stress

!=====================================================================
!=====================================================================

      SUBROUTINE sed_roughness
!--------------------------------------------------------------------!
! This routine computes ripples and sand-waves parameters, and then  !
! update the roughness length for hydrodynamics (rough_p) and/or for !
! sediment model (Zob) accordingly.                                  !
! By setting bedforms_rough to 1, 2, 3 or 4 within sediment.in, you  !
! can choose to update both hydro and sediment, only hydro or only   !
! sediment related roughness length                                  !
!                                                                    !
! Currently the procedure is:                                        !
!    1- computation of current alone skin stress                     !
!    2- computation of current-induced ripples (Soulsby, 1997)       !
!    3- computation of current-induced sand-waves (Van Rijn 1984, in !
!       Soulsby, 1997)                                               !
!    4- if WWM is activated, computes wave-induced ripples following !
!       Grant and Madsen (1982) or Nielsen (1992)                    !
!    5- Computes the total roughness length                          !
!    6- Depending on choosen option within sediment.in hydrodynamic, !
!       sediment or both roughness length are updated according with !
!                                                                    !
! Author: florian ganthy (fganthy@lnec.pt ; florian.ganthy@gmail.com)!
! Date:   2013/01/02                                                 !
!                                                                    !
! History: 2013/05 - F.Ganthy : Some update (computation at element  !
!                               centers)                             !
!          2013/05 - F.Ganthy : Implementation of wave-induced ripple!
!          2013/05 - F.Ganthy : Updates to the ripple predictor:     !
!                               - Changes on the total bedform       !
!                                 roughness computation              !
!                               - Add wave-ripple computation from   !
!                                 Nielsen (1992)                     !
!          2020/02 - B.Mengual : Correction of the Hsw expression    !
!                                for current sand waves              !
!                                                                    !
!--------------------------------------------------------------------!

      USE schism_msgp, ONLY: myrank,parallel_abort,exchange_e2d,exchange_p2d
      USE schism_glbl, ONLY: rkind,nvrt,nea,i34,elnode,dpe,eta2,idry_e,dfv,    &
     &                     kbe,uu2,vv2,ze,pi,grav,rough_p,errmsg,    &
     &                     rho0,dzb_min,np,nne,indel,area,ielg
      USE sed_mod,   ONLY: Cdb_min,Cdb_max,bottom,vonKar,Zob,uorb,   &
     &                     tp,bedforms_rough,iwave_ripple,           &
     &                     irough_bdld,bottom,isd50,idens,itauc,     &
     &                     izdef,izNik,izcr,izsw,izwr,izbld,izapp,   &
     &                     bed_rough,sed_debug


      IMPLICIT NONE


!- Local variables --------------------------------------------------!

      INTEGER     :: i,n1,n2,n3,ie,j
      INTEGER     :: kb0,kb1,kb2

      REAL(rkind) :: tau_cr,rhosed,d50,s,theta_cr,theta_b,dstar
      REAL(rkind) :: Hcr,Lcr,z0cr0,z0cr
      REAL(rkind) :: Hsw,Lsw,z0sw0,z0sw
      REAL(rkind) :: Hwr,Lwr,z0wr0,z0wr
      REAL(rkind) :: z0s,z0bld,z0
      REAL(rkind) :: hwat,hh,tau_c2,tau_wash,Ts
      REAL(rkind) :: cff1,cff2,cff3,cff4,cff5,wrk
      REAL(rkind) :: A,r,fwr,tau_w,theta_w,psi



!- User-defined parameters ------------------------------------------!

      REAL(rkind),PARAMETER :: nu=1.36d-6 ! Cinematic viscosity

!- Start Statement --------------------------------------------------!

      DO i = 1,nea
!---------------------------------------------------------------------
! ** Some preliminar initializations 
!---------------------------------------------------------------------
        ! Main bed sediment characteristics (previous time-step)
        tau_cr   = bottom(i,itauc)*rho0 ! Critical shear stress (N.m-2)
        rhosed   = bottom(i,idens)      ! Sediment density (kg.m-3)
        d50      = bottom(i,isd50)      ! Median grain size (m)
        s        = rhosed/rho0

        if(d50<=0.or.s<=1) then
          write(errmsg,*)'sed_roughness: d50<=0,',d50,s,rhosed,rho0,bottom(i,:),ielg(i)
          call parallel_abort(errmsg)
        endif

        dstar    = d50*(grav*(s-1.0d0)/(nu*nu))**(1.0d0/3.0d0)
        theta_cr = tau_cr/(grav*d50*(rhosed-rho0))
        if(theta_cr<=0) call parallel_abort('SED3D, sed_roughness; theta_cr<=0')
!'
        ! Bed forms history
        z0cr0  = bottom(i,izcr)       ! Current ripple roughness length (m)
        z0sw0  = bottom(i,izsw)       ! Sand waves roughness length (m)
        z0wr0  = bottom(i,izwr)       ! Wave ripple roughness length (m)
        ! Nikuradse roughness  (m)
        z0s    = d50/12.0d0
        bottom(i,izNik) = z0s

        IF(idry_e(i).EQ.1) THEN
          ! Dry element
          z0cr  = z0cr0
          z0sw  = z0sw0
          z0wr  = z0wr0
          z0bld = 0.0d0
        ELSE !wet
          ! Water depth
          hwat = dpe(i)+sum(eta2(elnode(1:i34(i),i)))/i34(i)

! ** Current-induced bed shear stress (skin friction only) 
!Error: limit shallows
          kb0 = kbe(i)
          kb1 = kbe(i)+1
          hh  = ze(kb1,i)-ze(kb0,i)                       
          IF(hh.GT.z0s) THEN
            IF((hh.LE.0.0d0).OR.(z0s.LE.0.0d0)) THEN
              CALL parallel_abort('SED-Roughness: Cd failed')
            ENDIF
            cff1 = 1.0d0/DLOG(hh/z0s)
            cff2 = vonKar*vonKar*cff1*cff1 !Cd
            wrk  = MIN(Cdb_max,MAX(Cdb_min,cff2)) !Cd
            cff3 = sum(vv2(kb1,elnode(1:i34(i),i)))/i34(i)
            cff4 = sum(uu2(kb1,elnode(1:i34(i),i)))/i34(i)
            !cff5 = DSQRT(cff4*cff4+cff3*cff3)
            ! Bottom stress [m^2/s/s]
            !tau_c2 = wrk*DSQRT((cff4*cff5)**2.0d0+(cff3*cff5)**2.0d0)
            tau_c2 = wrk*(cff3*cff3+cff4*cff4)
          ELSE !use approx.
!            kb2 = kbe(i)+2
            cff3 = sum(vv2(kb1,elnode(1:i34(i),i)))/i34(i) - &
            &      sum(vv2(kb0,elnode(1:i34(i),i)))/i34(i)
            cff4 = sum(uu2(kb1,elnode(1:i34(i),i)))/i34(i) - &
            &      sum(uu2(kb0,elnode(1:i34(i),i)))/i34(i) 
            cff5 = sum(dfv(kb1,elnode(1:i34(i),i))+dfv(kb0,elnode(1:i34(i),i)))/2/i34(i)
            !hh = ze(kb1,i)-ze(kb0,i)
            IF(hh<=0) CALL parallel_abort('SED-Roughness: div. by 0')
            ! Bottom stress
            tau_c2 = cff5*SQRT((cff4/hh)**2.0d0+(cff3/hh)**2.0d0)
          ENDIF ! End test on hh>z0s
          tau_c2 = tau_c2*rho0 ! Tau in Pa

! ** Current Ripples (Soulsby, 1997)
          IF(tau_c2.GT.tau_cr) THEN ! Initiation of motion
            ! Wash-out shear stress (Soulsby, 1997, eq. 85b)
            tau_wash = 0.8d0*grav*rho0*(s-1.0d0)*d50
            IF(tau_c2.LT.tau_wash)THEN
              ! Current ripples dimensions (Soulsby, 1997, eq. 81a-b)
              Lcr  = 1000.0d0*d50
              Hcr  = Lcr/7.0d0
              z0cr = MAX(1.d0*Hcr*Hcr/Lcr , z0cr0)
            ELSE !Wash-out
              Hcr  = 0.0d0
              Lcr  = 0.0d0 
              z0cr = 0.0d0
            ENDIF
          ELSE ! Keep roughness from previous time step
            z0cr = z0cr0
          ENDIF ! End test on initiation of motion

! ** Current Sand-waves (Van Rijn, 1984, In: Soulsby, 1997)
          IF(tau_c2.GT.tau_cr) THEN ! Initiation of motion
            ! Wash-out shear stress (Soulsby, 1997, eq. 83b)
            tau_wash = 26.0d0*tau_cr
            IF(tau_c2.LT.tau_wash)THEN
              ! Current ripples dimensions (Soulsby, 1997, eq. 83a-d)
              if(tau_cr==0) call parallel_abort('SED3D, sed_roughness; tau_cr==0')
!'
              Ts   = (tau_c2-tau_cr)/tau_cr
              Lsw  = 7.3d0*hwat
              !Hsw  = MAX(0.0d0 , 0.11d0*hwat*(d50*hwat)**3.0d0*    &
              !&      (1.0d0-DEXP(-0.5d0*Ts))*(25.0d0-Ts))
              !BM
              Hsw  = MAX(0.0d0 , 0.11d0*hwat*(d50/hwat)**0.3d0*    &
              &      (1.0d0-DEXP(-0.5d0*Ts))*(25.0d0-Ts)) 
              z0sw = MAX(1.d0*Hsw*Hsw/Lsw , z0sw0)
            ELSE !Wash-out
              Hsw  = 0.0d0
              Lsw  = 0.0d0 
              z0sw = 0.0d0
            ENDIF
          ELSE ! Keep roughness from previous time step
            z0sw = z0sw0
          ENDIF ! End test on initiation of motion

#ifdef USE_WWM
! ** Wave ripples (Grant and Madsen, 1982 or Nielsen, 1992)
          ! - Wave-induced friction coefficient - Swart, 1974 (eq. 60a-b)
          A = uorb(i)*tp(i)/(2.0d0/pi) !>=0
          r = A/(2.5d0*d50)
          IF(r.GT.0.0d0) THEN
            IF(r.LE.1.57d0)THEN
              fwr = 0.3d0
            ELSE
              fwr = 0.00251d0*DEXP(5.21d0*r**(-0.19d0))
            ENDIF
            tau_w = 0.5d0*rho0*fwr*uorb(i)*uorb(i)
          ELSE
            tau_w = 0.0d0
          ENDIF ! End test on r
          theta_w = tau_w/(grav*(rhosed-rho0)*d50)
!          if(s==1) call parallel_abort('sed_roughness; s==1')
          psi=uorb(i)*uorb(i)/((s-1.0d0)*grav*d50)

          ! - Wave ripples dimensions
          IF (iwave_ripple.EQ.0) THEN
            ! Grant and Madsen 1982, In: Soulsby, 1997, eq. 88a-f
            theta_b  = 1.8d0*theta_cr*(dstar**1.5d0/4.0d0)**0.6d0
            IF(theta_w.GT.theta_cr)THEN
              IF(theta_w.LE.theta_b)THEN
                !Initiation of motion
                Hwr  = 0.22d0*A*(theta_w/theta_cr)**(-0.16d0)
                Lwr  = Hwr/(0.16d0*(theta_w/theta_cr)**(-0.04d0))
              ELSE
                ! Washing condition
                Hwr = 0.48d0*A*((dstar**1.5d0/4.0d0)**0.8d0)*      &
                &     (theta_w/theta_cr)**(-1.5d0)
                Lwr = Hwr/(0.28d0*((dstar**1.5d0/4.0d0)**0.6d0)*   &
                &     (theta_w/theta_cr)**(-1.0d0))
                !z0wr = MAX(0.923d0*Hwr*Hwr/Lwr , z0wr0)
              ENDIF

              if(Lwr==0) then
                z0wr =z0wr0
              else
                z0wr = MAX(0.923d0*Hwr*Hwr/Lwr , z0wr0)
              endif

              ! Sediment transport component (Soulsby, 1997 eq. 91)
              IF(irough_bdld.EQ.1)THEN
                z0bld = 5.33d0*(s+0.5d0)*d50*theta_cr*             &
                &       ((theta_w/theta_cr)**0.5d0-0.7d0)**2.0d0
              ELSE
                z0bld = 0.0d0
              ENDIF
            ELSE
              ! No sediment motion
              Hwr   = 0.0d0
              Lwr   = 1.0d0
              z0wr  = z0wr0 !0.0d0
              z0bld = 0.0d0
            ENDIF

          ELSEIF (iwave_ripple.EQ.1) THEN

            ! Nielsen 1992, In: Soulsby, 1997, eq. 89a-d
            IF(theta_w.GT.theta_cr)THEN
              !Initiation of motion
              IF((psi.LT.156.0d0).AND.(theta_w.LT.0.831d0))THEN
                Hwr = (0.275d0-0.022d0*psi**0.5d0)*A
                if(0.182d0-0.24d0*theta_w**1.5d0==0) call parallel_abort('SED3D,sed_roughness; div. by 0 (11)')
!'
                Lwr = Hwr/(0.182d0-0.24d0*theta_w**1.5d0)
                if(Lwr==0) then
                  z0wr=z0wr0
                else
                  z0wr=max(0.267d0*Hwr*Hwr/Lwr,z0wr0)
                endif
              ELSE
                ! Washing condition
                Hwr  = 0.0d0
                Lwr  = 1.0d0
                z0wr = 0.0d0
              ENDIF

              ! Sediment transport component (Soulsby, 1997 eq. 92)
              IF(irough_bdld.EQ.1)THEN
                if(theta_w<0.05) then
                  !BM: add test on sed_debug to prevent too big nonfatal
                  !output files ...
                  IF (sed_debug .EQ. 1) THEN
                    write(12,*)'sed_roughness: theta_w<0.05:', &
       &theta_w,theta_cr,psi,Hwr,Lwr,d50,hwat,tau_c2,uorb(i),tp(i),tau_w
                  END IF
                  !call parallel_abort(errmsg)
                  z0bld = 0.0d0
                else
                  z0bld = 5.67d0*d50*sqrt(theta_w-0.05d0)
                endif
              ELSE
                z0bld = 0.0d0
              ENDIF
            ELSE
              ! No sediment motion
              Hwr   = 0.0d0
              Lwr   = 1.0d0
              z0wr  = z0wr0 !0.0d0
              z0bld = 0.0d0
            ENDIF
              
          ENDIF ! End test of wave-ripple formulation
#else
! WWM is not used
          z0wr  = 0.0d0
          z0bld = 0.0d0
#endif
        ENDIF ! idry_e

!---------------------------------------------------------------------
! ** Compute total roughness length
!---------------------------------------------------------------------
        z0 = MAX(z0s , MAX(z0cr , z0wr+z0bld)+z0sw)

!---------------------------------------------------------------------
! ** Consistency checks
!---------------------------------------------------------------------
        IF(z0.LE.0.d0) THEN
          WRITE(errmsg,*)'SED ROUGHNESS: roughness <= 0',ielg(i),z0
          CALL parallel_abort(errmsg)
        ENDIF
        !IF(z0.GE.dzb_min) THEN
        !  WRITE(errmsg,*)'SED ROUGHNESS: To high roughness!', &
        !  &              i,z0,z0cr,z0sw,z0wr,z0bld
        !  CALL parallel_abort(errmsg)
        !ENDIF
        !Limit z0
        z0=min(z0,dzb_min*1.d-1)

        bottom(i,izcr)  = z0cr  ! Current ripple roughness length (m)
        bottom(i,izsw)  = z0sw  ! Sand waves roughness length (m)
        bottom(i,izwr)  = z0wr  ! Wave ripple roughness length (m)
        bottom(i,izbld) = z0bld ! Sed Transport roughness length(m)
        bottom(i,izapp) = z0    ! Total roughness length

!---------------------------------------------------------------------
! ** Applying total roughness length to element for sediment transport
!---------------------------------------------------------------------
        IF(bedforms_rough==1) THEN
          ! Application of nikuradse roughness length
          Zob(i) = z0s
        ELSE IF(bedforms_rough==2)THEN
          ! Application of total bedforms roughness length
          Zob(i) = z0
        ENDIF

!---------------------------------------------------------------------
      ENDDO ! i=1,nea

      ! Exchange "ghost" elements 
      CALL exchange_e2d(Zob)             ! * FG : Check if necessary
      CALL exchange_e2d(bottom(:,izNik)) ! * FG : Check if necessary
      CALL exchange_e2d(bottom(:,izcr))  ! * FG : Check if necessary
      CALL exchange_e2d(bottom(:,izsw))  ! * FG : Check if necessary
      CALL exchange_e2d(bottom(:,izwr))  ! * FG : Check if necessary
      CALL exchange_e2d(bottom(:,izbld)) ! * FG : Check if necessary
      CALL exchange_e2d(bottom(:,izapp)) ! * FG : Check if necessary

!---------------------------------------------------------------------
! ** Convert roughness from elements to node for outputs
! ** And applying total roughness length to hydrodynamics roughness
!---------------------------------------------------------------------
      bed_rough   = 0.0d0
      DO i = 1,np
        cff1 = 0.0d0
        DO j = 1,nne(i)
          ie  = indel(j,i)
          cff1 = cff1 + area(ie)
          bed_rough(i)   = bed_rough(i)+bottom(ie,izapp)*area(ie)
        ENDDO ! End loop nne
        IF(cff1.EQ.0.0d0) THEN
          CALL parallel_abort('SED_ROUGHNESS: elem2nod: div. by 0')
        ELSE
          bed_rough(i) = bed_rough(i)/cff1
          IF(bedforms_rough.GE.1)THEN
            rough_p(i) = bed_rough(i)
          ENDIF
        ENDIF
      ENDDO ! End loop np
     
      call exchange_p2d(bed_rough)
      call exchange_p2d(rough_p)

!---------------------------------------------------------------------
!      IF(myrank.EQ.0) WRITE(16,*)'Leaving sed_roughness'
!---------------------------------------------------------------------
      END SUBROUTINE sed_roughness   

!=====================================================================
!=====================================================================

!n      SUBROUTINE stress_soulsby(cosv,sinv,dhx,dhy,alphas,sed_angle,  &
!n                 &              tauc0,tauc,i)
!n
!n!---------------------------------------------------------------------
!n!
!n!       AUTHOR:         Andre Fortunato - 96-10-22 modified by
!n!                       Ligia Pinto - 09-09-15
!n!       PURPOSE:
!n!
!n!---------------------------------------------------------------------
!n!
!n      USE schism_glbl, ONLY: errmsg,rkind
!n      USE schism_msgp, ONLY: parallel_abort 
!n
!n      IMPLICIT NONE
!n      SAVE     
!n  
!- Local variables --------------------------------------------------!
!n
!n      INTEGER      :: i
!n      REAL(rkind)  :: dhx,dhy
!n      REAL(rkind)  :: cosv,sinv,alphas,tauc
!n      REAL(rkind)  :: cosas,sinas,cosat,tan_alphat,aux,aux1
!n      REAL(rkind)  :: tauc0,sed_angle,sinat,aux2
!n
!n!- Start Statement --------------------------------------------------!
!n
!n!      alphas = DATAN(-(dhx*cosv+dhy*sinv)) !angulo alpha(s)
!n! Limitar o inclinacao a 0.9*(sed_angle) ou seja o angulo a atan(sed_angle*0.9)
!n!      aux=min(ABS(alphas),DATAN(sed_angle*0.9))
!n!      alphas = sign(aux,alphas)
!n!      IF (ABS(alphas)>DATAN(sed_angle))alphas=DATAN(sed_angle)*SIGN(1.0,alphas)
!n      IF (ABS(alphas)>DATAN(sed_angle))alphas=DATAN(sed_angle)*SIGN(real(1.0,8),alphas)
!n      cosas = DCOS(alphas)
!n      sinas = DSIN(alphas)
!n      tan_alphat = cosas*(dhx*sinv-dhy*cosv)
!n      cosat = DCOS(DATAN(tan_alphat))
!n      sinat = DSIN(DATAN(tan_alphat))
!n      tauc = tauc0*((sinas*cosat+DSQRT((cosas*cosas*sed_angle*sed_angle)-(sinat*sinat*sinas*sinas)))/(sed_angle))
!n!      IF (isnan(tauc)==-1) then
!n      IF (isnan(tauc).EQV..TRUE.) THEN
!n        WRITE(errmsg,*)'tauc is NaN "i"',tauc
!n        CALL parallel_abort(errmsg)
!n      ENDIF
!n
!n      RETURN
!n      END SUBROUTINE stress_soulsby

!=====================================================================
!=====================================================================

!n      SUBROUTINE stress(cosv,sinv,dhx,dhy,alphas,sed_angle,tauc0,tauc,i)
!n
!n!---------------------------------------------------------------------
!n!
!n!      AUTHOR:        Andre Fortunato - 96-10-22
!n!
!n!      PURPOSE:
!n!
!n!---------------------------------------------------------------------
!n
!n      USE schism_glbl, ONLY: errmsg,rkind
!n      USE schism_msgp, ONLY: parallel_abort
!n
!n      IMPLICIT NONE
!n      SAVE
!n
!n!- Local variables --------------------------------------------------!
!n
!n      INTEGER :: i
!n      REAL(rkind) :: dhx,dhy
!n      REAL(rkind) :: cosv,sinv,alphas,tauc
!n      REAL(rkind) :: cosas,sinas,cosat,tan_alphat,aux
!n      REAL(rkind) :: tauc0,sed_angle
!n      
!n
!n!- Start Statement --------------------------------------------------!
!n
!n!      alphas = DATAN(-(dhx*cosv+dhy*sinv)) !alpha(s)
!n!      IF (ABS(alphas).GT.DATAN(sed_angle))alphas=DATAN(sed_angle)*SIGN(1.0,alphas)
!n      IF (ABS(alphas).GT.DATAN(sed_angle))alphas=DATAN(sed_angle)*SIGN(REAL(1.0,rkind),alphas)
!n      cosas = DCOS(alphas)
!n      sinas = DSIN(alphas)
!n      tan_alphat = cosas*(dhx*sinv-dhy*cosv)
!n      cosat = DCOS(DATAN(tan_alphat))
!n      aux = sed_angle*sed_angle-tan_alphat*tan_alphat
!n      IF (aux.GT.0.d0) THEN
!n         tauc = tauc0*(cosas*cosat*DSQRT(1-((tan_alphat*tan_alphat)/(sed_angle*sed_angle)))+(sinas/sed_angle))
!n      ELSE
!n         tauc=tauc0*(sinas/sed_angle)
!n      ENDIF
!n
!n!ZYL
!n!       IF (isnan(tauc)==-1) then
!n       IF (isnan(tauc).EQV..TRUE.) THEN
!n            WRITE(errmsg,*)'tauc is NaN "i"',tauc,i,tauc0,sed_angle,cosas,cosat,tan_alphat,sinas
!n            CALL parallel_abort(errmsg)
!n       ENDIF
!n
!n
!n      RETURN
!n      END SUBROUTINE stress

!=====================================================================
!=====================================================================

