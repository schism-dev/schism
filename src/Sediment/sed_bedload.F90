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
! subroutine sed_bedload_vr
! subroutine sed_bedload_mpm
! subroutine bedchange_bedload
!
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
      ! Param for long bed slope effects
      REAL(rkind), PARAMETER :: alpha_bs = 1.d0
      ! Param for orth bed slope effects
      REAL(rkind), PARAMETER :: alpha_bn = 1.5d0

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
      htot = dpe(inea)+sum(eta2(elnode(1:i34(inea),inea)))/dble(i34(inea))

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

!---------------------------------------------------------------------
! - Bed_slope effects (Lesser et al. 2004)
! longitudinal bed slope \beta_s
! limit slope to 0.9*(sed_angle) (repose angle)
!---------------------------------------------------------------------

      cff  = dzdx*angleu+dzdy*anglev !tan(\beta_s)
      !sed_angle=TAN(\phi)
      cff1 = MIN(ABS(cff),0.9d0*sed_angle)*SIGN(1.0d0,cff)
      cff2 = ATAN(cff1) !\beta_s
      if(COS(cff2)==0.d0.or.sed_angle-cff1==0.d0) call parallel_abort('SED3D,sed_bedload_vr; div. by 0 (12)')
!'
      a_slopex = 1.0d0+alpha_bs*((sed_angle/(COS(cff2)*(sed_angle-cff1)))-1.0d0) !\alpha_s

!---------------------------------------------------------------------
! - Modify bedload transport due to longitudinal bed slope 
!---------------------------------------------------------------------

      FX_r(inea) = FX_r(inea)*a_slopex
      FY_r(inea) = FY_r(inea)*a_slopex

!---------------------------------------------------------------------
! - Transverse bed slope
!---------------------------------------------------------------------

      cff = -dzdx*anglev+dzdy*angleu !tan(\beta_n)
      cff1 = ABS(bustr(inea))+ABS(bvstr(inea))

      ! - Test used to prevent Inf & NaNs with very small bustr/bvstr 
      !   May still produce unrealistic values 
      IF(cff1<1.d-10) THEN
        a_slopey = 0.d0
      ELSE
        cff2 = SQRT(tau_ce(ised)/cff1)
        a_slopey=alpha_bn*cff2*cff !\alpha_n
      ENDIF

!---------------------------------------------------------------------
! - Add contribution of transverse to bed load 
!---------------------------------------------------------------------

      beta=FX_r(inea) !temp. save
      FX_r(inea) = FX_r(inea)-FY_r(inea)*a_slopey
      !FY_r(inea) = FY_r(inea)+FX_r(inea)*a_slopey
      FY_r(inea) = FY_r(inea)+a_slopey*beta

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
!          routines                                                  !
!                                                                    !
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
      do i = 1,nea
        IF(idry_e(i)==1) CYCLE
        FX_r(i) = FX_r(i)*morph_fac(ised)*bed_frac(1,i,ised) !m^2
        FY_r(i) = FY_r(i)*morph_fac(ised)*bed_frac(1,i,ised)
      enddo !i

!---------------------------------------------------------------------
! -  qsan=\int q_n*dt d\Gamma (although dimensioned up to npa, only 1:np
! are used)
! qsaxy is the sand flux at the element center integrated in time. Now 
! compute the line integral of qsaxy*normal along the control volume. 
! Add for each node in qsan. The unit normal is directed outward.
!---------------------------------------------------------------------
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
            cff1=(xel(id,ie)+xel(id1,ie))/2.d0 !xcj(isd2)-> between i and id1   
            cff2=(yel(id,ie)+yel(id1,ie))/2.d0 !ycj(isd2)-> between i and id1   
            cff3=(xel(id,ie)+xel(id3,ie))/2.d0 !xcj(isd1)-> between i and id3   
            cff4=(yel(id,ie)+yel(id3,ie))/2.d0 !ycj(isd1)-> between i and id3   
            qsan(i)=qsan(i)+FX_r(ie)*(cff4-cff2)+FY_r(ie)*(cff1-cff3) !m^3
          endif !ics
        enddo !j
      enddo !i=1,np

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
        bed_poro(i) = bed_poro(i)/dble(ks)
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
!   add volume? I think the idea is to scale the value by (1-porosity).
!FG. qsan is divided by (1-porosity) to take account for empty space
!    between grain (i.e. porosity) on bed level change

      qsan(1:np) = qsan(1:np)/(1.d0-bed_poro(1:np))

      ! Use JCG solver to calc hbed_ised
      hbed_ised=0.0d0 !initial guess
      !Right now zero-flux b.c. is imposed
      CALL solve_jcg(mnei_p,np,npa,it,moitn,mxitn,rtol,mcoefd,hbed_ised,qsan,bc_sed,lbc_sed)

!---------------------------------------------------------------------
! - Bed/bottom level change due to bedload transport in [m]
!---------------------------------------------------------------------
      !aggregate over all sed. classes (this routine is called inside a sed class loop)
      hbed(:) = hbed(:)+hbed_ised(:) 

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
        cff=sum(hbed_ised(elnode(1:i34(i),i)))/dble(i34(i))
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

