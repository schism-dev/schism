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
!--------------------------------------------------------------------

SUBROUTINE sed2d_main(it)

!--------------------------------------------------------------------
! Main sed2d routine
!
! Authors: Guillaume Dodet (guillaume.dodet@univ-brest.fr)
!          Thomas Guerin   (thomas.guerin@univ-lr.fr)
!
! History:
! 03/2013 - G.Dodet: - Corrected bugs in node centered FE-method;
!                    - Modified filtering section;
!                    - Added space-variable d50 input and Cdsed
!                      output;
! 04/2013 - G.Dodet: - Added nskip parameter to skip first iterations
!                      before updating bed level;
!                    - Split in 2 routines for each numerical method;
!                    - Added diffusive filter on total transport;
! 06/2013 - G.Dodet:  - Added filter on velocity before computing
!                       sediment fluxes to prevent development of
!                       oscillations along the coastline induced by
!                       wave-force/barotropic pressure gradient
!                       inconsistency (temporary solution);
!                     - Added diffusive filter on sediment fluxes
!                     - Added bedform predictors and new roughness;
!                     - Added timers;
! 07/2013 - G.Dodet:  - Converted output transport rate in kg/m/s;
!         - G.Dodet:  Added iterative shapiro filter for currents
!         - T.Guerin: Added Camenen and Larson (2011)
! 04/2014 - T.Guerin: - Removed part for getting offshore wave
!                       parameters (needed for Camemen and Larson
!                       formula and wave asymmetry calculation) and
!                       replaced it by local parameters computation
!                     - Added initialization of Cd_e, z0_e
!                     - Removed morphological time step (dtsed2d) to
!                       be consistent with multi-class mode (will be
!                       adapted in the future)
! 10/2016 - T.Guerin: Merging single-class and multi-class routines
!                     (one single routine is also kept for the
!                     different numerical schemes)
! 11/2016 - T.Guerin: Adding WENO scheme to solve the Exner equation
! 04/2017 - T.Guerin: Added morphological ramp factor
! 05/2017 - T.Guerin: - Added Wu an Lin (2014) formula
!                     - Added TR2004 formula (van Rijn et al., 2004)
!--------------------------------------------------------------------

  USE schism_glbl, ONLY : Cdp,dldxy,dp,dt,eta2,idry,idry_e,idry_s,ielg,iegl,&
                         &ipgl,indel,iplg,isdel,isidenode,isidenei2,        &
                         &elside,nea,nne,i34,elnode,np,npa,ns,nsa,pi,       &
                         &out_wwm,su2,sv2,timer_ns,rhosed,rkind,kbs,zs,nvrt,&
                         &ihydraulics,isblock_sd,errmsg,time_stamp,area, &
                         &in_dir,out_dir,len_in_dir,len_out_dir
  USE schism_msgp, ONLY : comm,exchange_e2d,exchange_p2d,exchange_s2d,      &
                         &ierr,myrank,parallel_abort,rtype,parallel_finalize
  USE hydraulic_structures, ONLY : nhtblocks
  USE sed2d_mod, ONLY : Cdsed,cflsed,d_class,diffac,dpdxy,dpdxy_e,       &
                       &dtsed2d,h0_sed,idrag,idsed2d,ifilt,imeth,        &
                       &imorpho,irough,iskip,islope,itrans,nskip,        &
                       &poro,qav,qb,qb_e,qdt_e,qfilter,qramp,qs,qs_e,    &
                       &qtot,qtot_e,ufilter,vc_area,                     &
                       &transfac,z0_e,z0cr_e,z0sw_e,z0wr_e,u2_tmp,v2_tmp,&
                       &bed_delta,d50,d90,F_class,h_inf,                 &
                       &h_lim_max,h_lim_min,h_top,h2,nb_class,           &
                       &nb_layer,morfac,ised_dump,imnp
  USE sed2d_friction
  USE sed2d_transport
  USE sed2d_filter
  USE sed2d_morpho

  IMPLICIT NONE

  INCLUDE 'mpif.h'

!- Arguments --------------------------------------------------------
  INTEGER, INTENT(IN) :: it
!- Local variables --------------------------------------------------
  INTEGER :: i,id,iel,inode,iside,j,k,k2,n,neigh,rank_tmp
  REAL(rkind) :: beta,D_star,D_star_class,d50_e,d90_e,hs,htot,      &
                &kpeak,rampfac,sum_F_class,suru,surv,tau_cr,        &
                &tau_cr_class,theta_cr,theta_cr_class,tmp,tp,Ucr_c, &
                &Ucr_c_class,Ucr_w,Ucr_w_class,u_star,udir,ue,uorb, &
                &uorbp,utmp,vtmp,ve,wdir,wdirc,wdirs,wlpeak,wtmp,z0
  REAL(rkind), DIMENSION(npa) :: dp_filt,dp_tmp,dp1,neigh2,sum_dh
  REAL(rkind), DIMENSION(nea) :: Cd_e,tmp_e
  REAL(rkind), DIMENSION(nsa) :: su2_tmp,su2_tmp1,su2_tmp2,sv2_tmp, &
                                &sv2_tmp1,sv2_tmp2,tmp_s
  REAL(rkind), DIMENSION(npa,2) :: qdt
  REAL(rkind), DIMENSION(nb_class) :: F_class_e
  REAL(rkind), DIMENSION(npa,nb_class) :: dh_class
  REAL(rkind), DIMENSION(nea,2,nb_class) :: qb_class_e,qs_class_e,  &
                                           &qtot_class_e
!- Dumping-related variables ---------------------------------------
  INTEGER,save :: ne_dump
  INTEGER, allocatable :: ie_dump(:)
  REAL(rkind),save :: t_dump  !time in dumping option
  REAL(rkind), allocatable :: vol_dump(:)
  logical, save :: first_call=.true.
!--------------------------------------------------------------------

!- Compute ramp factor
  IF(it<iskip) THEN
    rampfac = 0.d0
  ELSE
    IF(qramp>0) THEN
      rampfac = TANH(2.d0*(it-iskip)*dt/qramp)
    ELSE
      rampfac = 1.d0
    ENDIF
  ENDIF

#ifdef INCLUDE_TIMING
  wtmp = mpi_wtime() !start of timer
#endif

!- Side velocities filter before computing transport ----------------
!- Final vel used is s[u,v]2_tmp
  su2_tmp=0
  sv2_tmp=0
  do i=1,nsa
    if(idry_s(i)==1) cycle

    do k=kbs(i),nvrt-1
      su2_tmp(i)=su2_tmp(i)+(su2(k,i)+su2(k+1,i))/2*(zs(k+1,i)-zs(k,i))
      sv2_tmp(i)=sv2_tmp(i)+(sv2(k,i)+sv2(k+1,i))/2*(zs(k+1,i)-zs(k,i))
    enddo !k
    htot=zs(nvrt,i)-zs(kbs(i),i)
    if(htot<=0) CALL parallel_abort('SED2D: htot<=0')
    su2_tmp(i)=su2_tmp(i)/htot
    sv2_tmp(i)=sv2_tmp(i)/htot
  enddo !i 

  IF(ufilter /= 0) THEN
    su2_tmp1 = su2_tmp !su2(1,:)
    sv2_tmp1 = sv2_tmp !sv2(1,:)
    su2_tmp2 = su2_tmp !su2(1,:)
    sv2_tmp2 = sv2_tmp !sv2(1,:)
!- Skip if the element has a side at the wet-dry interface 
    DO i = 1,nsa
      IF(isdel(1,i)>0 .AND. isdel(2,i)>0) THEN !then resident and internal
        IF(idry_e(isdel(1,i))+idry_e(isdel(2,i))==1) THEN !wet/dry interface
          su2_tmp1(i) = 0.d0 
          sv2_tmp1(i) = 0.d0
        ELSE
          su2_tmp1(i) = su2_tmp(i) !su2(1,i) 
          sv2_tmp1(i) = sv2_tmp(i) !sv2(1,i)
        ENDIF
      ENDIF
    ENDDO !nsa
    CALL exchange_s2d(su2_tmp1)
    CALL exchange_s2d(sv2_tmp1)

!- Iterative shapiro filter applied on the current field    
!JZ: lon/lat; 
    DO k = 1,ufilter
      DO i = 1,ns
        if(isdel(2,i)==0.or.idry_s(i)==1) CYCLE
        if(ihydraulics/=0.and.nhtblocks>0) then
          if(isblock_sd(1,i)/=0) cycle
        endif

        suru = 0.d0
        surv = 0.d0
        DO j = 1,4
          id = isidenei2(j,i)
          suru = suru+su2_tmp1(id)
          surv = surv+sv2_tmp1(id)
        ENDDO
        su2_tmp2(i) = su2_tmp1(i)+0.25d0*(suru-4.d0*su2_tmp1(i))
        sv2_tmp2(i) = sv2_tmp1(i)+0.25d0*(surv-4.d0*sv2_tmp1(i))
      ENDDO !ns
      DO i = 1,ns
        if(isdel(2,i)==0.or.idry_s(i)==1) CYCLE
        if(ihydraulics/=0.and.nhtblocks>0) then
          if(isblock_sd(1,i)/=0) cycle
        endif

        su2_tmp(i) = su2_tmp2(i)
        sv2_tmp(i) = sv2_tmp2(i)
      ENDDO
      CALL exchange_s2d(su2_tmp)
      CALL exchange_s2d(sv2_tmp)
      su2_tmp1 = su2_tmp
      sv2_tmp1 = sv2_tmp
    ENDDO !k=1,ufilter
    IF(myrank==0) WRITE(16,*)'done filtering velocities (sed2d)'
  ENDIF !ufilter/=0

#ifdef INCLUDE_TIMING
      timer_ns(11) = timer_ns(11)+mpi_wtime()-wtmp !end of timer 
#endif 

!- Dumping/dredging
  if(ised_dump/=0) then
    !For 1st call (including hot), init. read
    if(first_call) then
      !Time stamps in this file must be one of the time steps
      open(18,file=in_dir(1:len_in_dir)//'sed_dump.in',status='old')
      read(18,*)
      do 
        read(18,*,iostat=k)t_dump,ne_dump !time in sec
        if(k/=0) then
          ised_dump=0 !reset
          exit
        endif

        if(t_dump>=time_stamp) exit

        read(18,*) !vol
      enddo
    endif !first_call

    if(ised_dump/=0.and.abs(t_dump-time_stamp)<1.e-4) then !in case end of file
      allocate(ie_dump(ne_dump),vol_dump(ne_dump),stat=k)
      if(k/=0) call parallel_abort('SED2D: alloc (9)')
      read(18,*)(ie_dump(k),vol_dump(k),k=1,ne_dump)
      if(myrank==0) write(16,*)'SED2D,start dumping at time:',real(t_dump),ne_dump

      !Modify depth but not bottom() or bed_frac
      do k=1,ne_dump
        iel=ie_dump(k) !global index
        if(iegl(iel)%rank==myrank) then
          i=iegl(iel)%id !local index
          tmp=vol_dump(k)/area(i) !m
          do j=1,i34(i)
            dp(elnode(j,i))=dp(elnode(j,i))-tmp
          enddo !j
        endif !iegl
      enddo !k

      deallocate(ie_dump,vol_dump)

      !Prep for next record
      read(18,*,iostat=k)t_dump,ne_dump 
      if(k/=0) ised_dump=0 !reset
    endif !ised_dump/=0
  endif !ised_dump/=0

!--------------------------------------------------------------------
!- Compute sediment transport at element center for each grain size
!- class in the top layer of the sediment column
!--------------------------------------------------------------------

!- Initialization
  Cd_e(:)             = 0.d0 !Drag coefficient (-)
  dpdxy_e(:,:)        = 0.d0 !Bed slope (m/m)
  z0_e(:)             = 0.d0 !Roughness length (m)
  qb_class_e(:,:,:)   = 0.d0 !Bed load transport for each grain size class (m3/s/m)
  qs_class_e(:,:,:)   = 0.d0 !Suspended load transport for each grain size class (m3/s/m)
  qtot_class_e(:,:,:) = 0.d0 !Total transport for each grain size class (m3/s/m)
  qb_e(:,:)           = 0.d0 !Total bed load transport (sum over grain size classes) (m3/s/m)
  qs_e(:,:)           = 0.d0 !Total suspended load transport (sum over grain size classes) (m3/s/m)
  qtot_e(:,:)         = 0.d0 !Total total transport (sum over grain size classes) (m3/s/m)

!- Start loop over elements
  DO i = 1,nea

#ifdef INCLUDE_TIMING
     wtmp = mpi_wtime() !start of timer
#endif

!- Initialization
     d50_e     = 0.d0 !d50 (m)
     d90_e     = 0.d0 !d90 (m)
     htot      = 0.d0 !Total depth (m)
     ue        = 0.d0 !X-velocity (m/s)
     ve        = 0.d0 !Y-velocity (m/s)
     uorb      = 0.d0 !Init orbital velocity computed by wwm (m/s)
     hs        = 0.d0 !Significant wave height computed by wwm(m)
     tp        = 0.d0 !Peak period computed by wwm (s)
     wlpeak    = 0.d0 !Wave length associated to peak period by wwm(m)
     uorbp     = 0.d0 !Peak orbital velocity (m/s)
     wdir      = 0.d0 !Mean average energy transport direction (rad)
     wdirc     = 0.d0 !cosine value needed for calculating mean wave direction in each element
     wdirs     = 0.d0 !sine value needed for calculating mean wave direction in each element

     dp1 = dp !Store depth before bed update for output purpose

     DO j = 1,i34(i)
        inode = elnode(j,i)
        htot = htot+(eta2(inode)+dp(inode))/i34(i) !Total water depth
     ENDDO

     IF(idry_e(i)==1.OR.htot<=h0_sed) CYCLE

     DO j = 1,i34(i)
        inode = elnode(j,i)
        iside = elside(j,i)

!- Velocity 
        ue = ue+su2_tmp(iside)/i34(i)
        ve = ve+sv2_tmp(iside)/i34(i)

!- d50 and d90 corresponding to the top layer
        d50_e = d50_e+d50(inode,1)/i34(i)
        d90_e = d90_e+d90(inode,1)/i34(i)

!- bed slope  
        dpdxy_e(i,1) = dpdxy_e(i,1)+dp(inode)*dldxy(j,1,i)
        dpdxy_e(i,2) = dpdxy_e(i,2)+dp(inode)*dldxy(j,2,i)

!- Wave parameters 
#ifdef USE_WWM
        uorb   = uorb+out_wwm(inode,22)/i34(i)
        hs     = hs+out_wwm(inode,1)/i34(i)
        tp     = tp+out_wwm(inode,12)/i34(i)
        wlpeak = wlpeak+out_wwm(inode,17)/i34(i)
        wdirc  = wdirc+(COS((270.d0-out_wwm(inode,9))*pi/180.d0))/i34(i)
        wdirs  = wdirs+(SIN((270.d0-out_wwm(inode,9))*pi/180.d0))/i34(i)
#endif

!- Drag coefficient
        IF(idrag==1)THEN
          Cd_e(i) = Cd_e(i)+Cdp(inode)/i34(i)
        ENDIF
     ENDDO! j = 1,i34

#ifdef USE_WWM
     wdir = ATAN2(wdirs,wdirc)
#endif

#ifdef INCLUDE_TIMING
     timer_ns(4) = timer_ns(4)+mpi_wtime()-wtmp !end of timer 
     wtmp = mpi_wtime()                         !start of timer
#endif 

!- Compute bed slope in streamwise direction used for slope effect on
!  total transport (Soulsby, 1997). Positive if flows runs uphill.
     IF(islope == 1) THEN
       udir = ATAN2(ve,ue)
       beta = -(dpdxy_e(i,1)*COS(udir)+dpdxy_e(i,2)*SIN(udir))
       beta = MIN(beta,0.6d0) !May lead to negative transport otherwise
     ELSE
       beta = 0.d0
     ENDIF

!- Compute threshold parameters, roughness and drag coefficients
!  related to the global sediment characteristics (i.e. d50 and d90)
     CALL compute_thresh(d50_e,d90_e,htot,tp,D_star,Ucr_c,Ucr_w,  &
                         &theta_cr,tau_cr)

     SELECT CASE(irough)
       CASE(1) !Skin friction only
         z0_e(i) = d50_e/12.d0
       CASE(2) !Bedform associated roughness (Soulsby, 1997)
         CALL bedform_predictor_s97(d50_e,htot,ue,ve,uorb,tp,     &
           &tau_cr,theta_cr,z0cr_e(i),z0wr_e(i),z0sw_e(i),z0_e(i))
       CASE(3) !Bedform associated roughness (Van-Rijn, 2007)
         CALL bedform_predictor_vr07(d50_e,htot,ue,ve,uorb,       &
           &z0cr_e(i),z0wr_e(i),z0sw_e(i),z0_e(i))
     END SELECT

     IF(IABS(idrag).GT.1) CALL compute_drag(d50_e,htot,z0_e(i),ue,ve,idrag,Cd_e(i))

#ifdef INCLUDE_TIMING
     timer_ns(5) = timer_ns(5)+mpi_wtime()-wtmp !end of timer 
     wtmp = mpi_wtime()                         !start of timer
#endif 

#ifdef USE_DEBUG
!- Test if parameters are physically correct
     if(htot<0.d0 .or. SQRT(ue*ue+ve*ve)>7.5d0 .or. d50_e>0.05d0.or. &
       &d90_e>0.10d0 .or. SQRT(dpdxy_e(i,1)**2.d0+dpdxy_e(i,2)**2.d0)&
       &>5.d0 .or. uorb>10.d0 .or. hs>15.d0 .or. tp>30.d0 .or.       &
       &wlpeak>1500.d0 .or. Cd_e(i)>10.d0 .or. z0>0.5d0) then
       write(12,*)'Warning (1) for it:',it,' and el: ', ielg(i)
       write(12,*)'Check following values: '
       write(12,*)'total depth = ',htot
       write(12,*)'current velocity = ',SQRT(ue*ue+ve*ve)
       write(12,*)'d50 = ',d50_e
       write(12,*)'d90 = ',d90_e
       write(12,*)'slope = ',SQRT(dpdxy_e(i,1)**2.d0+dpdxy_e(i,2)**2.d0)
       write(12,*)'orbital velocity = ',uorb
       write(12,*)'hs = ',hs
       write(12,*)'tp = ',tp
       write(12,*)'wlpeak = ',wlpeak
       write(12,*)'Cd = ',Cd_e(i)
       write(12,*)'z0 = ',z0_e(i)
!       CALL parallel_abort('Sed2d: chech value in outputs/nonfatal')
     endif
#endif

!- Start loop over grain size classes
     DO k=1,nb_class

!- Compute sediment class fractions at element center
        if (k==1) then
           F_class_e(:) = 0.d0
           do j=1,i34(i)
              inode = elnode(j,i)
              do k2 = 1,nb_class
                 F_class_e(k2) = F_class_e(k2)+F_class(inode,k2,1)/i34(i)
              enddo
           enddo !j=1,i34(i)
        endif

!- Compute threshold parameters related to grain size class k
        CALL compute_thresh(d_class(k),d90_e,htot,tp,D_star_class,    &
               &Ucr_c_class,Ucr_w_class,theta_cr_class,tau_cr_class)

!- Compute sediment transport for grain size class k
!  NB: we always keep Cd computed from the d50

        SELECT CASE(itrans)
         CASE(1) !Engelund-Hansen (1967)
           CALL eh67(Cd_e(i),d_class(k),ue,ve,beta,qtot_class_e(i,:,k))

         CASE(2) !Ackers and White (1973)
           CALL aw73(Cd_e(i),d_class(k),ue,ve,htot,beta,D_star_class, &
                  &qtot_class_e(i,:,k))

         CASE(3) !Soulsby - Van Rijn (Soulsby, 1997)
           CALL svr97_bedl(Cd_e(i),d_class(k),ue,ve,htot,dpdxy_e(i,:),&
                  &Ucr_c_class,tau_cr_class,uorb,qb_class_e(i,:,k))

           CALL svr97_susp(Cd_e(i),d_class(k),ue,ve,D_star_class,     &
                  &Ucr_c_class,uorb,qs_class_e(i,:,k))

         CASE(4) !Van-Rijn (2007)
           CALL vr07_bedl(Cd_e(i),d_class(k),ue,ve,htot,dpdxy_e(i,:),uorb,&
                  &Ucr_c_class,Ucr_w_class,tau_cr_class,qb_class_e(i,:,k))

           CALL vr07_susp(d_class(k),ue,ve,uorb,D_star_class,Ucr_c_class, &
                  &Ucr_w_class,qs_class_e(i,:,k))

         CASE(5) !Camenen and Larson (2011)
           CALL cl11_tot(d_class(k),D_star_class,dpdxy_e(i,:),hs,htot,&
                  &theta_cr_class,tp,ue,ve,wdir,wlpeak,               &
                  &qb_class_e(i,:,k),qs_class_e(i,:,k))

         CASE(6) !Wu and Lin (2014)
           CALL wl14_tot(d_class(k),d_class(:),d50_e,d90_e,D_star_class,  &
                  &F_class_e(:),htot,dpdxy_e(i,:),ue,ve,hs,tp,wdir,wlpeak,&
                  &qb_class_e(i,:,k),qs_class_e(i,:,k))

         CASE(7) !TRANSPOR2004 (van Rijn et al., 2004)
           CALL tr04_tot(d_class(k),d50_e,d90_e,htot,dpdxy_e(i,:),ue,ve,hs,&
                  &tp,wdir,wlpeak,qb_class_e(i,:,k),qs_class_e(i,:,k))
        END SELECT

!- Apply slope effect and compute total transport for grain size class k
        if (itrans==3.or.itrans==4.or.itrans==5.or.itrans==6) then
          qb_class_e(i,:,k)   = qb_class_e(i,:,k)*(1.d0-1.6d0*beta)
          qs_class_e(i,:,k)   = qs_class_e(i,:,k)*(1.d0-1.6d0*beta)
          qtot_class_e(i,:,k) = qb_class_e(i,:,k)+qs_class_e(i,:,k)
        elseif (itrans==7) then
          !NB: no slope effect with tr04 formula because already done in tr04_tot
          qtot_class_e(i,:,k) = qb_class_e(i,:,k)+qs_class_e(i,:,k)
        endif

!- Multiply the transport by the corresponding grain size fraction
        qb_class_e(i,:,k)   = F_class_e(k)*qb_class_e(i,:,k)
        qs_class_e(i,:,k)   = F_class_e(k)*qs_class_e(i,:,k)
        qtot_class_e(i,:,k) = F_class_e(k)*qtot_class_e(i,:,k)

!- Apply ramp, correction, and diffusion factors if required
        qtot_class_e(i,1,k) = rampfac*transfac*(qtot_class_e(i,1,k)+(1.d0-poro)*diffac* &
                             &abs(qtot_class_e(i,1,k))*dpdxy_e(i,1))
        qtot_class_e(i,2,k) = rampfac*transfac*(qtot_class_e(i,2,k)+(1.d0-poro)*diffac* &
                             &abs(qtot_class_e(i,2,k))*dpdxy_e(i,2))

!- Sum over grain size classes
        qb_e(i,:)   = qb_e(i,:)   + qb_class_e(i,:,k)
        qs_e(i,:)   = qs_e(i,:)   + qs_class_e(i,:,k)
        qtot_e(i,:) = qtot_e(i,:) + qtot_class_e(i,:,k)

#ifdef INCLUDE_TIMING
        timer_ns(6) = timer_ns(6)+mpi_wtime()-wtmp !end of timer
        wtmp = mpi_wtime()                         !start of timer
#endif

     ENDDO !k=1,nb_class
  ENDDO !i=1,nea

#ifdef INCLUDE_TIMING
    timer_ns(7) = timer_ns(7)+mpi_wtime()-wtmp !end of timer
    wtmp = mpi_wtime()                         !start of timer
#endif


!--------------------------------------------------------------------
!- Compute bed change at nodes for each grain size class
!--------------------------------------------------------------------

  DO k=1,nb_class

!- Apply diffusive filter on total transport if required
     IF(qfilter == 1) THEN
       CALL sed2d_filter_diffu(qtot_class_e(:,1,k),tmp_e,nea)
       qtot_class_e(:,1,k) = tmp_e
       CALL sed2d_filter_diffu(qtot_class_e(:,2,k),tmp_e,nea)
       qtot_class_e(:,2,k) = tmp_e
       CALL exchange_e2d(qtot_class_e(:,1,k))
       CALL exchange_e2d(qtot_class_e(:,2,k))
       IF(myrank==0) WRITE(16,*)'done filtering transport rate (sed2d) for grain size class:',k
       !Recompute total transport (sum over grain size classes)
       if(k==1) qtot_e(:,:) = 0.d0
       qtot_e(:,1) = qtot_e(:,1) + qtot_class_e(:,1,k)
       qtot_e(:,2) = qtot_e(:,2) + qtot_class_e(:,2,k)
     ENDIF

!- Time-integrated sediment transport (m^2) for grain size class k
!  NB: qdt_e is related to the transport of grain size class k
!      but not total transport over all grain size classes
     qdt_e(:,:) = 0.d0 !initialization
     qdt_e(:,1) = qtot_class_e(:,1,k)*dt
     qdt_e(:,2) = qtot_class_e(:,2,k)*dt

!- Apply morphological factor
     qdt_e(:,1) = qdt_e(:,1)*morfac
     qdt_e(:,2) = qdt_e(:,2)*morfac

     IF (imorpho==1 .AND. it>iskip) THEN
       !Compute bottom evolution for grain size class k
       if (imeth==1) then
         CALL base_scheme(it)
       elseif (imeth==2) then
         CALL weno_scheme(it)
       endif
       !Save bottom evolution for grain size class k
       dh_class(:,k) = bed_delta(:)
     ENDIF
  ENDDO !k=1,nb_class


!--------------------------------------------------------------------
!-                    Update bathymetry                             -
!--------------------------------------------------------------------

  IF(imorpho == 1 .AND. it>iskip) THEN
    sum_dh(:) = 0.d0
    DO i=1,npa
       DO k=1,nb_class
          sum_dh(i) = sum_dh(i) + dh_class(i,k)
       ENDDO
       sum_dh(i) = sum_dh(i)*imnp(i) !apply morphological ramp factor
       dp(i) = dp(i) + sum_dh(i)     !update bathy
    ENDDO
  ENDIF


!--------------------------------------------------------------------
!-                    Update grain size fractions                   -
!--------------------------------------------------------------------

  IF (nb_class>1.and.imorpho==1.and.it>iskip) THEN
    DO i=1,npa

       IF (sum(dh_class(i,:))==0.d0) CYCLE

!- Compute new sediment class fractions
      DO k=1,nb_class
        IF (sum_dh(i)>0.d0 .or. &
           &(sum_dh(i)==0.d0.and.dh_class(i,k)>0.d0)) THEN !EROSION case

!YJZ: this check may be too restrictive?
          IF(sum_dh(i)>=h_lim_min) THEN
             write(errmsg,*)'negative thickness for layer #2: &
                 &need to increase h_lim_min, or decrease time step, or decrease dtsed2d:',sum_dh(i),h_lim_min
             CALL parallel_abort(errmsg)
          ENDIF

          !New fractions for surface layer; draw from layer #2 so top layer
          !thickness stays same
          F_class(i,k,1)=(h_top*F_class(i,k,1)+sum_dh(i)*F_class(i,k,2)-dh_class(i,k))/h_top

          IF(F_class(i,k,1)<0.d0) F_class(i,k,1)=0.d0


        ELSEIF (sum_dh(i)<0.d0 .or. &
               &(sum_dh(i)==0.d0.and.dh_class(i,k)<0.d0)) THEN !DEPOSITION case

!YJZ: this check may be too restrictive?
          IF(ABS(sum_dh(i))>=h_lim_max) THEN
            write(errmsg,*)'layer #2 thickness too large: &
              &need to increase h_lim_max, or decrease time step, or decrease dtsed2d:',sum_dh(i),h_lim_max
            CALL parallel_abort(errmsg)
          ENDIF

          !New fractions for surface layer
          !Deposit into layer #1 (-dh), and swap out same amount to layer #2
          F_class(i,k,1)=((h_top+sum_dh(i))*F_class(i,k,1)-dh_class(i,k))/h_top

          IF(F_class(i,k,1)<0.d0) F_class(i,k,1)=0.d0

          !New sub-surface layer fractions
          F_class(i,k,2)=(h2(i)*F_class(i,k,2)-sum_dh(i)*F_class(i,k,1))/(h2(i)-sum_dh(i)) !>=0
        ENDIF !sum_dh(i) sign
      ENDDO !k=1,nb_class

!-    Update layer #2
      IF (sum_dh(i)>0.d0) THEN !EROSION case
        !New sub-surface layer thickness
        h2(i) = h2(i) - sum_dh(i) !should >0 b/cos sum_dh(i)<=h_lim_min

        IF (h2(i)<h_lim_min) THEN !Merging sub-surface layer with layer #3
          !New sub-surface layer fractions and thickness
          if(h2(i)+h_inf==0) call parallel_abort('SED2D: (3)')
          F_class(i,:,2)=(h_inf*F_class(i,:,3)+h2(i)*F_class(i,:,2))/(h2(i)+h_inf) !>0
          h2(i) = h2(i) + h_inf

          !Adding new layer => updating fractions - implies infinite supply
          !of sediment
          DO n=3,nb_layer-1
            F_class(i,:,n) = F_class(i,:,n+1)
          ENDDO
        ENDIF !h2(i)

      ELSEIF (sum_dh(i)<0.d0) THEN !DEPOSITION
        !New sub-surface layer thickness
        h2(i) = h2(i) - sum_dh(i) !sum_dh(i)<0

        IF (h2(i)>h_lim_max) THEN !Splitting sub-surface layer into 2 layers
          !New sub-surface layer thickness
          h2(i) = h2(i) - h_inf 
          if(h2(i)<h_lim_min) then
            write(errmsg,*)'SED2D: check h_lim_max etc: ',h2(i),h_lim_min,h_lim_max
            call parallel_abort(errmsg)      
          endif

          !Removing lowest layer => updating fractions
          DO n=1,nb_layer-2
            F_class(i,:,nb_layer-n+1) = F_class(i,:,nb_layer-n)
          ENDDO
        ENDIF !h2(i)>h_lim_max
      ENDIF !sum_dh(i)

!- Compute sum of fractions and correct fractions if sum /= 1
      DO n=1,nb_layer
        sum_F_class = 0.d0
        DO k=1,nb_class
          sum_F_class = sum_F_class + F_class(i,k,n)
        ENDDO

        IF (sum_F_class/=1.d0) THEN
          if(sum_F_class==0) then
            if (n==1 .and. imnp(i)==0) then     !non-erodable surface layer case
              F_class(i,:,n) = 0.d0
              F_class(i,nb_class,n) = 1.d0
            else
              write(errmsg,*)'SED2D: all eroded; ',n,iplg(i)
              call parallel_abort(errmsg)
            endif
          else
            F_class(i,:,n) = F_class(i,:,n)/sum_F_class
          endif
        ENDIF
      ENDDO !nb_layer

!- Compute new d50 (weighted geometric mean)
      DO n=1,nb_layer
        d50(i,n) = d_class(1)**F_class(i,1,n)
        DO k=2,nb_class
          d50(i,n) = d50(i,n) * d_class(k)**F_class(i,k,n)
        ENDDO
      ENDDO !nb_layer

!- Compute new d90 (to be improved)
      d90(i,:) = 2.5d0*d50(i,:)
    ENDDO !i=1,npa

  ENDIF !imorpho=1 .and. it>iskip

#ifdef INCLUDE_TIMING
  timer_ns(9) = timer_ns(9)+mpi_wtime()-wtmp !end of timer 
  wtmp = mpi_wtime()                         !start of timer
#endif 

!--------------------------------------------------------------------
!-                      Apply user-defined filter                   -
!--------------------------------------------------------------------
  IF (imorpho==1.AND.it>iskip.AND.ifilt>0) THEN
    SELECT CASE(ifilt)
     CASE(1)
       CALL sed2d_filter_extrema(dp,dp_filt)
     CASE(2)
       CALL sed2d_filter_slope(dp,dp_filt)
     CASE(3)
       CALL sed2d_filter_slope(dp,dp_tmp)
       CALL sed2d_filter_extrema(dp_tmp,dp_filt)
    END SELECT
! prevent filters from changing depths where imnp=0
    do i = 1, npa
      if (imnp(i) .eq. 1) dp(i) = dp_filt(i)
    end do
    CALL exchange_p2d(dp)
    IF(myrank==0)WRITE(16,*)'done filtering new bathy. (sed2d)'
  ENDIF

!--------------------------------------------------------------------
!-                     Check if new bathymetry /= NaN               -
!--------------------------------------------------------------------
  DO i = 1,npa
    IF(dp(i)/=dp(i)) THEN
      WRITE(errmsg,*)'Sed2d: NaN in depth at step:',it,' and node: ', iplg(i),';dp = ',dp(i)
      CALL parallel_abort(errmsg)
    ENDIF
  ENDDO

#ifdef INCLUDE_TIMING
  timer_ns(10) = timer_ns(10)+mpi_wtime()-wtmp !end of timer 
  wtmp = mpi_wtime()                         !start of timer
#endif 

!--------------------------------------------------------------------
!-       Interpolate variables at node for output purpose and also for 
!        passing back Cdsed to schism_step 
!--------------------------------------------------------------------
  Cdsed = 0.d0
  cflsed = 0.d0
  dpdxy = 0.d0
  qav = 0.d0
  qdt = 0.d0
  qb = 0.d0
  qs = 0.d0
  qtot = 0.d0

  DO i = 1,np
    IF(idry(i) == 1) CYCLE

    neigh = 0
    DO j = 1,nne(i)
      iel = indel(j,i)
      Cdsed(i)   = Cdsed(i)+Cd_e(iel)
      dpdxy(i,:) = dpdxy(i,:)+dpdxy_e(iel,:)
      qav(i,:)   = qav(i,:)+qtot_e(iel,:)*dt
      qdt(i,:)   = qdt(i,:)+qtot_e(iel,:)*dt
      qb(i,:)    = qb(i,:)+qb_e(iel,:)
      qs(i,:)    = qs(i,:)+qs_e(iel,:)
      qtot(i,:)  = qtot(i,:)+qtot_e(iel,:)
      neigh = neigh+1
    ENDDO !j

    if(neigh==0.or.eta2(i)+dp1(i)==0) call parallel_abort('SED2D: div. by 0 (12)')
    Cdsed(i)   = Cdsed(i)/neigh
    dpdxy(i,1) = dpdxy(i,1)/neigh
    dpdxy(i,2) = dpdxy(i,2)/neigh
    qav(i,:)   = rhosed*qav(i,:)/(neigh*dt*dtsed2d) !(kg/m/s)
    qdt(i,:)   = qdt(i,:)/neigh                     !(m3/m)
    qb(i,:)    = rhosed*qb(i,:)/neigh               !(kg/m/s)
    qs(i,:)    = rhosed*qs(i,:)/neigh               !(kg/m/s)
    qtot(i,:)  = rhosed*qtot(i,:)/neigh             !(kg/m/s)
    cflsed(i) = 3.d0*SQRT(qdt(i,1)*qdt(i,1)+qdt(i,2)*qdt(i,2))/  &
     &((eta2(i)+dp1(i))*SQRT(vc_area(i)))
  ENDDO !np

  CALL exchange_p2d(Cdsed(:))
  DO i=1,2
    CALL exchange_p2d(dpdxy(:,i))
    CALL exchange_p2d(qav(:,i))
    CALL exchange_p2d(qdt(:,i))
    CALL exchange_p2d(qb(:,i))
    CALL exchange_p2d(qs(:,i))
    CALL exchange_p2d(qtot(:,i))
  ENDDO
  CALL exchange_p2d(cflsed(:))

  u2_tmp = 0.d0
  v2_tmp = 0.d0
  neigh2 = 0.d0
  DO i = 1,nsa
    DO j = 1,2
      inode = isidenode(j,i)
      IF(idry(inode) == 1) CYCLE
        u2_tmp(inode) = u2_tmp(inode)+su2_tmp(i)
        v2_tmp(inode) = v2_tmp(inode)+sv2_tmp(i)
        neigh2(inode) = neigh2(inode)+1
    ENDDO !j=1,2
  ENDDO !nsa
  CALL exchange_p2d(u2_tmp(:)) !check if necessary
  CALL exchange_p2d(v2_tmp(:)) !check if necessary
  CALL exchange_p2d(neigh2(:)) !check if necessary

  DO i = 1,npa
    IF(idry(i) == 1) CYCLE
    if(neigh2(i)==0) call parallel_abort('sed2d_main_node: impossible (1)')
    u2_tmp(i) = u2_tmp(i)/neigh2(i)
    v2_tmp(i) = v2_tmp(i)/neigh2(i)   
  ENDDO !npa
 
  tmp=sum(dp)
  WRITE(12,*)'SED2D, it=',it,' ; sum of depth=',tmp
  IF(myrank==0) WRITE(16,*)'done computing total sediment load (sed2d) and &
    &exiting sed2d_main; check nonfatal* for sum of depths'
  if(tmp/=tmp) call parallel_abort('sed2d_main_node: NaN in sum')

  first_call=.false.

#ifdef INCLUDE_TIMING
  timer_ns(8) = timer_ns(8)+mpi_wtime()-wtmp !end of timer 
#endif 

  END SUBROUTINE sed2d_main
