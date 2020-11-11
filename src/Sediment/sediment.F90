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

      subroutine sediment(it,moitn,mxitn,rtol,dave,tot_bedmass)
!--------------------------------------------------------------------!
!  This routine computes the sediment sources and sinks and adds     !
!  then the global sediment tracer fields. Currently, it includes    !
!  the following:                                                    !
!                                                                    !
!  * Vertical settling of sediment in the water column.              !
!  * Erosive and depositional flux interactions of sediment          !
!    between water column and the bed.                               !
!  * Transport of multiple grain sizes.                              !
!  * Bed layer stratigraphy.                                         !
!  * Bed morphology.                                                 !
!  * Bedload based on Meyer Peter Mueller or Van Rijn, 2007, Journal !
!    of Hydraulic Engineering                                        !
!  * Bedload slope term options: Damgaard et al, 1997, Journal       !
!    of Hydraulic Engineerig v123, p 1130-1138; Antunes do Carmo,    !
!    1995 PhD Thesis; Lesser et al, 2004, Coastal Engineering, v 51, !
!    p 883-915.                                                      !
!                                                                    !
!  * Seawater/sediment vertical level distribution:                  !
!                                                                    !
!         W-level  RHO-level                                         !
!                                                                    !
!            N     _________                                         !
!                 |         |                                        !
!                 |    N    |                                        !
!          N-1    |_________|  S                                     !
!                 |         |  E                                     !
!                 |   N-1   |  A                                     !
!            3    |_________|  W                                     !
!                 |         |  A                                     !
!                 |    3    |  T                                     !
!            2    |_________|  E                                     !
!                 |         |  R                                     !
!                 |    2    |                                        !
!            1    |_________|_____ bathymetry                        !
!                 |/////////|                                        !
!                 |    1    |                                        !
!            1    |_________|  S                                     !
!                 |         |  E                                     !
!                 |    2    |  D                                     !
!            2    |_________|  I                                     !
!                 |         |  M                                     !
!                 |  Nbed-1 |  E                                     !
!        Nbed-1   |_________|  N                                     !
!                 |         |  T                                     !
!                 |  Nbed   |                                        !
!         Nbed    |_________|                                        !
!                                                                    !
!                                                                    !
! This subroutine is adapted from ROMS routines                      !
! Copyright (c) 2002-2007 The ROMS/TOMS Group                        !
!   Licensed under a MIT/X style license                             !
!   See License_ROMS.txt                                             !
!                                                                    !
! Author: Ligia Pinto                                                !
! Date: xx/08/2007                                                   !
!                                                                    !
! History: 2015 - Joseph Zhang: overhauled some parts
!          2012/12 - F.Ganthy : form homogenisation of sediments     !
!          routines                                                  !
!          2012/12 - F.Ganthy : modifications for Bottom Composition !
!                               Generation (BCG) purpose (added      !
!                               output files - bed median grain size,!
!                               bed fraction - and modification of   !
!                               morphodynamics management            !
!          2013/01 - F.Ganthy : Implementation of roughness predictor!
!          2013/01 - F.Ganthy : Implementation of avalanching        !
!          2013/03 - F.Ganthy : Implementation of bedmass filter     !
!          2013/03 - F.Ganthy : Implementation wave-induced bedload  !
!                               transport                            !
!          2013/04 - F.Ganthy : Implementation of wave-current bottom!
!                               stress                               !
!          2013/05 - F.Ganthy : Add different sediment behavior:     !
!                                - MUD-like or SAND-like             !
!          2013/06 - F.Ganthy : Modification for BCG related to the  !
!                               multiple bed layer model             !
!
!          2020/02 - B.Mengual                                       !
!                  * Implementation of a method to compute and update!
!                    a porosity representative of non-cohesive       !
!                    matrix                                          !
!                  * Implementation of several methods to compute    !
!                    the current-induced bottom shear stress         !
!                    (tau_option)                                    !
!                  * Wave asymmetry computations: restructuration    !
!                    and introduction of filtering possibilities     !
!                  * Implementation of bedload flux caused by        !
!                    wave acceleration skewness                      !
!                  * Implementation of bedload transport formulations!
!                    - Soulsby and Damgaard (2005) initially         !
!                      introduced by Anouk de Bakker                 !
!                    - Wu and Lin (2014) adapted from SED2D          !
!                      (original code: Thomas Guerin)                !
!                  * Add the possibility to filter bedload fluxes    !
!                  * Add the possibility to use the WENO scheme      !
!                    developed by Guerin et al (2016) to solve the   !
!                    Exner equation in bedchange_bedload             !
!                    (firstly implemented by Anouk de Bakker)        !
!                  * Implementation of a limiter on bedload fluxes   !
!                    based on the active layer criterion             !
!                  * Suspended transport: correction of depo_mss     !
!                    for ised_bc_bot=1                               !
!                  * Introduction of a maximum active layer thickness!
!                  * Morphodynamics:                                 !
!                    - common morph_fac over all sediment classes    !
!                    - implementation of morphological ramp value    !
!                      imnp, initially defined in the imorphogrid.gr3!
!                      input file (same procedure than in SED2D)     !                 
!                  * New outputs at nodes: erosion flux, deposition  ! 
!                    flux, top layer porosity, bedload fluxes due to !
!                    wave acceleration                               !
!                                                                    !
!--------------------------------------------------------------------!

      USE sed_mod
      USE schism_glbl, ONLY : rkind,nvrt,ne,nea,npa,np,irange_tr,ntrs,idry_e,   &
     & idry,area,znl,dt,i34,elnode,xctr,yctr,   &
     & kbe,ze,pi,nne,indel,tr_el,flx_bt,    &
     & errmsg,ielg,iplg,nond_global,iond_global,&
     & ipgl,nope_global,np_global,dp,h0,dpe,    &
     & iegl,out_wwm,pi,eta2,dp00,dldxy,we,time_stamp, &
     & itur,Phai,in_dir,out_dir,len_in_dir,len_out_dir, &
     & uu2,vv2,dav
                            !Tsinghua group:+alphd,im_pick_up,Two_phase_mix
                            !phai_m !1120:-alphd,im_pick_up,Two_phase_mix,phai_m  +itur,Phai
      USE schism_msgp
      use misc_modules

      IMPLICIT NONE
      include 'mpif.h'

!      SAVE

!- Local variables --------------------------------------------------!

      INTEGER,INTENT(IN)     :: it          ! time step
      INTEGER,INTENT(IN)     :: moitn,mxitn ! JCG solver int and max 
                                            ! iterations
      REAL(rkind),INTENT(IN) :: rtol        ! Relative tolerance 
      REAL(rkind),INTENT(IN) :: dave(nea)   ! Depth-averaged vel. at centroids
      REAL(rkind),INTENT(OUT) :: tot_bedmass !total bed mass [kg]


      INTEGER, PARAMETER :: top = 1      ! Top layer of bed
      REAL(rkind), PARAMETER :: eps = 1.0d-14

!      INTEGER :: nwild(nea+12)
!      INTEGER :: ibnd,isd,isd00,ind1,ind2,ndo,ndf(npa)
      INTEGER,save :: ie,nd,ip,ifl
      INTEGER,save :: Ksed,i,indx,ised,j,k,ks,l
      INTEGER,save :: nm1,nm2,nm3 !bnew

      REAL(rkind),save :: time,ta,tmp
      REAL(rkind),save :: cff, cff1, cff2, cff3, cffL, cffR, dltL, dltR
      REAL(rkind),save :: cu, cff4, cff6, aref, cff7, cff8, cff9
      REAL(rkind),save :: thck_avail,thck_to_add,eros_mss,depo_mss,flux_eros,flux_depo
 
      ! - For suspended sediment
!      INTEGER, dimension(nvrt,nea)     :: ksource
!      REAL(rkind)                      :: Hz_inv3
!      REAL(rkind), DIMENSION(nvrt)     :: Hz_inv
!      REAL(rkind), DIMENSION(nvrt)     :: Hz_inv2
!      REAL(rkind), DIMENSION(nvrt,nea) :: FC
!      REAL(rkind), DIMENSION(nvrt,nea) :: qc
!      REAL(rkind), DIMENSION(nvrt,nea) :: qR
!      REAL(rkind), DIMENSION(nvrt,nea) :: qL
!      REAL(rkind), DIMENSION(nvrt,nea) :: WR
!      REAL(rkind), DIMENSION(nvrt,nea) :: WL
      ! - For bed load
      REAL(rkind),save :: cff5
      REAL(rkind),save :: bedld_mass
      REAL(rkind),save :: smgdr, osmgd, Umag
      REAL(rkind),save :: tauc0
      REAL(rkind),save :: derx1,derx2,derx3,dery1,dery2,dery3
      REAL(rkind),save :: yp,xp,flux      
      REAL(rkind), allocatable :: dep_mass(:,:)

      ! - For MPM bed load
      REAL(rkind),save :: alphas,tauc

      ! - For VR bed load
      REAL(rkind),save :: tsta,dpar

      ! - For morphology
      INTEGER,save     :: kbed

      REAL(rkind),save, allocatable :: hdep(:),hbed(:) !depth change (in dt) due to suspended load and bedload
      REAL(rkind),save, allocatable :: hbed_ised(:),hdep_nd(:)
      REAL(rkind),save, allocatable :: qsan(:) 
      REAL(rkind),save, allocatable :: dhnd(:) !total depth change in dt (=sus+bedload)

      ! - For waves
      INTEGER          :: t
      REAL(rkind),save :: htot
      REAL(rkind),save :: kpeak
      REAL(rkind)      :: wdirc,wdirs,T1,T0,T2,U0,tv

      !Dumping
      INTEGER,save :: ne_dump
      INTEGER, save, allocatable :: ie_dump(:)
      REAL(rkind),save :: t_dump  !time in dumping option
      REAL(rkind), save, allocatable :: vol_dump(:)
    
      logical, save :: first_call=.true.

      ! Variables associated to tau_option
      INTEGER :: kk,kb1,kb0
      REAL(rkind) :: wkk,wkkm1,dzb

      REAL(rkind),DIMENSION(nea) :: tmp_e ! temporary array used in filters


!- Start Statement --------------------------------------------------!
      allocate(dep_mass(nea,ntr_l),stat=i)
      if(i/=0) call parallel_abort('SED: alloc failed')

      if(.not.allocated(hdep)) allocate(hdep(nea))
      if(.not.allocated(hbed)) allocate(hbed(npa))
      if(.not.allocated(hbed_ised)) allocate(hbed_ised(npa))
      if(.not.allocated(hdep_nd)) allocate(hdep_nd(npa))
      if(.not.allocated(qsan)) allocate(qsan(npa))
      if(.not.allocated(dhnd)) allocate(dhnd(npa))

      ! Used to laternately store bed_mass from 2 steps
      nstp = 1+MOD(it-1,2)
      nnew = 3-nstp
      !bnew = nnew

      ! Initialize  variables
      ks       = 0
      !ksource  = 0
      Hz       = 0.d0
!      qL       = 0.d0
!      qR       = 0.d0
!      qc       = 0.d0
      cffR     = 0.d0
      cffL     = 0.d0
!      FC       = 0.d0
      cff      = 0.d0
      dep_mass = 0.d0
      hdep_nd  = 0.d0
      hdep     = 0.d0
      hbed     = 0.d0
      angleu   = 0.d0
      anglev   = 0.d0
      bedldu   = 0.d0
      bedldv   = 0.d0
      sedcaty  = 0.d0 !Tsinghua group

      DO i=1,nea
        bustr(i) = 0.d0
        bvstr(i) = 0.d0
      ENDDO

      ! Bedload
      FX_r = 0.0d0
      FY_r = 0.0d0
      Qaccu = 0.0d0
      Qaccv = 0.0d0
      r_accu = 0.0d0
      r_accv = 0.0d0
      Qaccun = 0.0d0
      Qaccvn = 0.0d0

      ! Wave parameters initialisation
      hs     = 0.0d0
      tp     = 0.0d0
      wlpeak = 0.0d0
      uorb   = 0.0d0
      uorbp  = 0.0d0
      dirpeak = 0.0d0
      wdir = 0.0d0

      ! Wave asymmetry (Elfrink)
      U_crest = 0.d0
      U_trough = 0.d0
      T_crest = 0.d0
      T_trough = 0.d0
      Uorbi_elfrink = 0.0d0

      ! Bottom shear stress
      ustress=0.0d0
      vstress=0.0d0

      ! RUNTIME
      time=dt*it

!---------------------------------------------------------------------
!     Dumping/dredging
!---------------------------------------------------------------------
      if(ised_dump/=0) then
        !For 1st call (including hot), init. read
        if(first_call) then
          !!Time stamps in this file must be one of the time steps
          open(18,file=in_dir(1:len_in_dir)//'sed_dump.in',status='old')
          read(18,*)
          do 
            read(18,*,iostat=k)t_dump,ne_dump !time in sec
            if(k/=0) then
              ised_dump=0 !reset
              exit
            endif

            if(t_dump>=time) exit

            read(18,*) !vol
          enddo
        endif !first_call
 
        if(ised_dump/=0.and.abs(t_dump-time)<1.e-4) then !in case end of file
          allocate(ie_dump(ne_dump),vol_dump(ne_dump),stat=l)
          if(l/=0) call parallel_abort('SED: alloc (9)')
          read(18,*)(ie_dump(l),vol_dump(l),l=1,ne_dump)
          if(myrank==0) write(16,*)'start dumping at time:',real(t_dump),ne_dump

          !Modify depth, bed(), but not bottom() or bed_frac
          do l=1,ne_dump
            ie=ie_dump(l) !global index
            if(iegl(ie)%rank==myrank) then
              i=iegl(ie)%id !local index
              cff=vol_dump(l)/area(i) !m
              tmp=bed(top,i,ithck)+cff
              if(tmp>0) then !enough sed on top
                bed(top,i,ithck)=tmp
              else !re-init.
                tmp=max(0.d0,sum(bed(:,i,ithck))+cff)
                bed(:,i,ithck)=tmp/Nbed
              endif
              do ised=1,ntr_l
                bed_mass(:,i,1,ised)=bed(:,i,ithck)*Srho(ised)*(1.0d0-bed(:,i,iporo))*bed_frac(:,i,ised)
                bed_mass(:,i,2,ised)=bed_mass(:,i,1,ised)
              enddo !ised

              do j=1,i34(i)
                dp(elnode(j,i))=dp(elnode(j,i))-cff
              enddo !j
            endif !iegl
          enddo !l

          deallocate(ie_dump,vol_dump)

          !Prep for next record
          read(18,*,iostat=k)t_dump,ne_dump 
          if(k/=0) ised_dump=0 !reset
        endif !ised_dump/=0
      endif !ised_dump/=0

!---------------------------------------------------------------------
! - Get wave parameters (defined at nodes) and converte to element
!   centres
!---------------------------------------------------------------------
#ifdef USE_WWM
      DO i = 1,nea
        IF (idry_e(i).EQ.1) CYCLE
!        htot  = 0.0d0
        kpeak = 0.0d0
        wdirc = 0.0d0
        wdirs = 0.0d0
        DO j = 1,i34(i)
!          htot      = htot + (dp(elnode(j,i))+eta2(elnode(j,i)))/i34(i)
          hs(i)     = hs(i) + out_wwm(elnode(j,i),1)/i34(i)
          tp(i)     = tp(i) + out_wwm(elnode(j,i),12)/i34(i)
          wlpeak(i) = wlpeak(i) + out_wwm(elnode(j,i),17)/i34(i) !Peak wave length
          uorb(i)   = uorb(i) + out_wwm(elnode(j,i),22)/i34(i) !orbital vel.
          dirpeak(i) = dirpeak(i) + out_wwm(elnode(j,i),18)/i34(i) ! direction of peak Anouk
          wdirc  = wdirc+(COS((270.d0-out_wwm(elnode(j,i),9))*pi/180.d0))/i34(i) !BM
          wdirs  = wdirs+(SIN((270.d0-out_wwm(elnode(j,i),9))*pi/180.d0))/i34(i) !BM
        ENDDO !j 
        wdir(i) = ATAN2(wdirs,wdirc) !BM
        ! * FG - Uorbp unused now
        !kpeak    = 2.0d0*pi/wlpeak(i)
        !uorbp(i) = pi*hs(i)/(tp(i)*DSINH(kpeak*htot))
      ENDDO ! End loop nea
#endif

!---------------------------------------------------------------------
! - Compute bottom stress
!   Test of whether WWM is used or not is performed within the subroutine
!---------------------------------------------------------------------

!BM: Several methods to compute current-induced bottom shear stress 
!    (tau_option, see sediment.in).
!    Current components (ustress, vstress) are then used in sed_current_stress 
!    to compute BSS, or in sed_bedload routines to derive current direction.
  
      DO i=1,nea
        IF (idry_e(i)==1) CYCLE  
        kb0 = kbe(i)
        kb1 = kbe(i)+1
        dzb=ze(kb1,i)-ze(kb0,i)

        IF ((tau_option .EQ. 1) .OR. (tau_option .EQ. 3 .AND. dzb .GE. zstress)) THEN
            ustress(i) = sum(uu2(kb1,elnode(1:i34(i),i)))/i34(i)
            vstress(i) = sum(vv2(kb1,elnode(1:i34(i),i)))/i34(i)

        ELSE IF ((tau_option .EQ. 2) .OR. (tau_option .EQ. 3 .AND. dzb .LT. zstress)) THEN
          DO j=1,i34(i)
            nd=elnode(j,i)
            kk=nvrt
            DO k=kb1,nvrt
              IF (znl(k,nd)-znl(kb0,nd) .GE. zstress) THEN
                kk=k
                EXIT
              END IF
            END DO
  
            IF (kk .LT. nvrt) THEN
              wkk=(zstress+znl(kb0,nd)-znl(kk-1,nd))/(znl(kk,nd)-znl(kk-1,nd))
              wkkm1=(znl(kk,nd)-zstress-znl(kb0,nd))/(znl(kk,nd)-znl(kk-1,nd))
              ustress(i)=ustress(i)+(wkkm1*uu2(kk-1,nd) + wkk*uu2(kk,nd))
              vstress(i)=vstress(i)+(wkkm1*vv2(kk-1,nd) + wkk*vv2(kk,nd))
   
            ELSE
              ustress(i)=ustress(i)+dav(1,nd)
              vstress(i)=vstress(i)+dav(2,nd)
  
            END IF
          END DO
          ustress(i) = ustress(i)/i34(i)
          vstress(i) = vstress(i)/i34(i)
  
        ELSE
          call parallel_abort('SED: unknown tau_option')
        END IF ! test on tau_option and dzb
      END DO
  
      CALL exchange_e2d(ustress)
      CALL exchange_e2d(vstress)
 
      ! Compute current/wave-induced and total shear stresses
      CALL sed_wavecurrent_stress()

!---------------------------------------------------------------------
! - Sediment debug output
!---------------------------------------------------------------------

      IF(sed_debug==1) CALL sed_write_debug(it)  

!---------------------------------------------------------------------
!           *** Compute bedload sediment transport. ***
!                Adapted from [ROMS sed_bedload.F]
!---------------------------------------------------------------------
      IF(myrank==0) WRITE(16,*)'SED: start bedload...'

      ! Compute some constant bed slope parameters.
      ! Friction angle = 33º
      ! sed_angle in rad van rijn 35�
      sed_angle = TAN(30.0_r8*pi/180.0_r8)


      ! Compute bedload boundary conditions
      IF (sed_morph.GE.1) THEN

        ! Identify open boundaries and impose flux of 0 for JCG
        bc_sed = -9999 !flags
        DO i=1,nope_global
          DO j=1,nond_global(i)
            nd = iond_global(i,j) !global
            IF(ipgl(nd)%rank==myrank) THEN
              ip = ipgl(nd)%id
              bc_sed(ip) = 0
            ENDIF !ipgl(nd)%rank==myrank
          ENDDO !End loop nond_global
        ENDDO !End loop nope_global

        ! Compute b.c. flag for all nodes for the matrix
        DO i=1,npa
          IF(bc_sed(i)>-9998) THEN
            lbc_sed(i)=.TRUE.
          ELSE
            lbc_sed(i)=.FALSE.
          ENDIF
        ENDDO ! End loop npa

      ENDIF ! sed_morph.GE.1

!---------------------------------------------------------------------
!                 **   Wave effects **
!---------------------------------------------------------------------
#ifdef USE_WWM
!-------------------
! - Wave asymmetry
!-------------------
      IF (bedload > 0) THEN
        DO i=1,nea
          IF (idry_e(i)==1) CYCLE
  
          dzdx=dot_product(dp(elnode(1:i34(i),i)),dldxy(1:i34(i),1,i))
          dzdy=dot_product(dp(elnode(1:i34(i),i)),dldxy(1:i34(i),2,i))
  
          ! Compute vector components of bottom stress induced by
          ! currents.
          ! Here it is assumed that effects of waves on current
          ! direction are already taken into account.
          IF (tau_c(i).EQ.0.d0) THEN
            angleu = 0.d0
            anglev = 0.d0
          ELSE
            angleu = bustr(i)/tau_c(i) !direction of stress
            anglev = bvstr(i)/tau_c(i) !
          ENDIF
  
          ! Total water depth
          htot=dpe(i)+sum(eta2(elnode(1:i34(i),i)))/i34(i) 
  
  
          ! Compute wave asymmetry features if bedload is taken into
          ! account and if iasym is equal to 1
          IF (hs(i)>0.01d0 .AND. tp(i)>0.01d0 &
              & .AND. wlpeak(i)>0.01d0 .AND. htot > h0) THEN

              IF (iasym ==1) THEN

                CALL wave_asymmetry_Elfrink(hs(i),tp(i),wdir(i),htot,dzdx,dzdy,&  ! inputs
                       & U_crest(i),U_trough(i), & ! Outputs
                       & T_crest(i),T_trough(i), &
                       & Uorbi_elfrink(i,:))
 
              END IF

              IF (iasym == 0                      & ! Airy waves
                & .OR. U_crest(i)  <= 0.0d0       & ! or issue with
                & .OR. U_trough(i) <= 0.0d0       & ! iasym = 1
                & .OR. T_crest(i)  <= 0.0d0       &
                & .OR. T_trough(i) <= 0.0d0       &
                & .OR. U_crest(i)  /= U_crest(i)  &
                & .OR. U_trough(i) /= U_trough(i) &
                & .OR. T_crest(i)  /= T_crest(i)  &
                & .OR. T_trough(i) /= T_trough(i) ) THEN

                tmp=dsinh((2.d0*pi/wlpeak(i))*htot)
                IF (tmp==0.d0) call parallel_abort('wave asym issue (1)')
                U_crest(i) = dmax1(pi*hs(i)/tp(i)/tmp,0.001d0) !>0
                U_trough(i) = U_crest(i)
                T_crest(i) = tp(i)/2.d0 !>0
                T_trough(i) = T_crest(i) !>0

                T1=tp(i)/4.0d0
                T0=2.0d0*T1
                T2=3.0d0*T1

                IF (T1==0) call parallel_abort('wave asym issue (2)')

                U0 = (U_crest(i)*T0-U_trough(i)*(1.d0-T0))/(T0-T1)

                IF (U0 .GT. (0.25d0*U_crest(i))) THEN
                  U0 = 0.25d0*U_crest(i)
                ENDIF


                DO t=0,ech_uorb
                  tv=DBLE(t)/DBLE(ech_uorb)

                  IF ((tv .GE. 0.d0) .AND. (tv .LE. T1)) THEN
                    Uorbi_elfrink(i,t) = U_crest(i)*DSIN(pi/2.d0*tv/T1)

                  ELSEIF ((tv .GT. T1) .AND. (tv .LE. T0)) THEN
                    Uorbi_elfrink(i,t) = U_crest(i)*DCOS(pi/2.d0*(tv-T1)/(T0-T1)) &
                               - U0*DSIN(pi*(tv-T1)/(T0-T1))

                  ELSEIF ((tv .GT. T0) .AND. (tv .LE. T2)) THEN
                    Uorbi_elfrink(i,t) = -U_trough(i)*DSIN(pi/2.d0*(tv-T0)/(T2-T0))

                  ELSEIF ((tv .GT. T2) .AND. (tv .LE. 1.d0)) THEN
                    Uorbi_elfrink(i,t) = -U_trough(i)*DCOS(pi/2.d0*(tv-T2)/(1.d0-T2))

                  ENDIF
                ENDDO
              END IF

              IF (iasym<0 .OR. iasym>1) THEN
                CALL parallel_abort('incorrect iasym parameter (O or 1)')
              END IF

          END IF ! if waves

        END DO ! nea loop

        IF (iasym == 1 .AND. elfrink_filter == 1) THEN
          ! Filtering of Elfrink outputs
          ! -> diffusive filter from sed2d (Dodet, 2013)

          CALL sed2d_filter_diffu(U_crest(:),tmp_e,nea)
          U_crest(:) = tmp_e
          CALL sed2d_filter_diffu(U_trough(:),tmp_e,nea)
          U_trough(:) = tmp_e
          CALL sed2d_filter_diffu(T_crest(:),tmp_e,nea)
          T_crest(:) = tmp_e
          CALL sed2d_filter_diffu(T_trough(:),tmp_e,nea)
          T_trough(:) = tmp_e

          CALL exchange_e2d(U_crest(:))
          CALL exchange_e2d(U_trough(:))
          CALL exchange_e2d(T_crest(:))
          CALL exchange_e2d(T_trough(:))

        END IF

      END IF ! if bedload>0


!---------------------------------------------------------------------
! - Compute acceleration-skewness induced transport: Qaccu/Qaccv [m2]
!---------------------------------------------------------------------

      IF (bedload_acc > 0) THEN
        DO i=1,nea
          IF (idry_e(i)==1) CYCLE
  
          htot=dpe(i)+sum(eta2(elnode(1:i34(i),i)))/i34(i)

          IF (hs(i)>0.01d0 .AND. tp(i)>0.01d0 &
            & .AND. wlpeak(i)>0.01d0 .AND. htot > h0) THEN
  
            IF (bedload_acc == 1) THEN
              CALL acceleration_bedload_hoe2003(i,                    &
                                             & tp(i),wdir(i),         &
                                             & U_crest(i),U_trough(i),&
                                             & Uorbi_elfrink(i,:),    &
                                             & Qaccu(i),              &
                                             & Qaccv(i))
  
            ELSE IF (bedload_acc == 2) THEN
              CALL acceleration_bedload_dub2015(i,                    &
                                             & hs(i),tp(i),wdir(i),   &
                                             & htot,                  &
                                             & Qaccu(i),              &
                                             & Qaccv(i))

             ELSE
              CALL parallel_abort('incorrect bedload_acc parameter')
             END IF
          ELSE ! no waves
            Qaccu(i)=0.0d0
            Qaccv(i)=0.0d0
          END IF
        END DO ! nea
  
        IF (bedload_acc_filter == 1) THEN
          ! Filtering of Qaccu/Qaccv
          ! -> diffusive filter from sed2d (Dodet, 2013)
          CALL sed2d_filter_diffu(Qaccu(:),tmp_e,nea)
          Qaccu(:) = tmp_e
          CALL exchange_e2d(Qaccu(:))

          CALL sed2d_filter_diffu(Qaccv(:),tmp_e,nea)
          Qaccv(:) = tmp_e
          CALL exchange_e2d(Qaccv(:))
        END IF
      END IF ! bedload_acc > 0
#endif /*USE_WWM*/


!---------------------------------------------------------------------
! - Loop over all non-cohesive sediment classes
!---------------------------------------------------------------------
      DO ised=1,ntr_l 

        !Use in some routines
        smgd  = (Srho(ised)/rhom-1.0d0)*g*Sd50(ised)
        if(smgd<=0) call parallel_abort('SED3D: smgd<=0')
        osmgd = 1.0d0/smgd
        smgdr = SQRT(smgd)*Sd50(ised)

!---------------------------------------------------------------------
! - Loop over elements
!---------------------------------------------------------------------
        DO i=1,nea

          IF (idry_e(i)==1) CYCLE

!---------------------------------------------------------------------
! - Compute bottom slopes
!   Jan check x-y-assignment
!jl. Derivatives of shape functions
!    An equivalent implementation is found in schism_main for dl 
!    Dphi/Dx, Dphi/Dy
!---------------------------------------------------------------------
          dzdx=dot_product(dp(elnode(1:i34(i),i)),dldxy(1:i34(i),1,i))
          dzdy=dot_product(dp(elnode(1:i34(i),i)),dldxy(1:i34(i),2,i))

          ! Compute vector components of bottom stress induced by curents
          ! Here it is assumed that effects of waves on current direction
          ! are already taken into account
          IF (tau_c(i).EQ.0.d0) THEN
            angleu = 0.d0
            anglev = 0.d0
          ELSE
            angleu = bustr(i)/tau_c(i) !direction of stress
            anglev = bvstr(i)/tau_c(i) !
          ENDIF

!---------------------------------------------------------------------
! - Computation of bedload
!   returns FX_r and FY_r for class ised - each ised loop will
!   overrwrite F[XY]_r (common via sed_mod)
!---------------------------------------------------------------------

          IF (bedload == 0) THEN
            ! No bedload
            FX_r(i)=0.0d0
            FY_r(i)=0.0d0

          ELSEIF (bedload == 1) THEN

            IF(iSedtype(ised).EQ.0) THEN
              !* MUD-like sediment type, no bedload
              FX_r(i) = 0.0d0
              FY_r(i) = 0.0d0
            ELSEIF(iSedtype(ised).EQ.1) THEN
              !* SAND-like sediment (0.05 < D50 < 2.0 mm) 
              ! van Rijn bedload can be applied
              CALL sed_bedload_vr(ised,i,dave)
            ELSE !IF(iSedtype(ised).EQ.2) THEN
              !* GRAVEL-like sediment (>= 2.0 mm)
              ! Need a specific bedload transport subroutine
              ! NOT AVAILABLE YET
              WRITE(errmsg,*)'SED: GRAVEL type not available yet'
              CALL parallel_abort(errmsg)
            ENDIF

          ELSEIF (bedload == 2) THEN
            ! bedload mpm UNFINISHED
            !CALL sed_bedload_mpm()
            WRITE(errmsg,*)'SED: MPM bedload not ready yet'
            CALL parallel_abort(errmsg)
              
          ! Anouk Soulsby and Damgaard (2005) 
          ELSEIF (bedload == 3) THEN ! Anouk
            IF(iSedtype(ised).EQ.0) THEN
              !* MUD-like sediment type, no bedload
              FX_r(i) = 0.0d0
              FY_r(i) = 0.0d0
            ELSEIF(iSedtype(ised).EQ.1) THEN 
              !* SAND-like sediment 
              CALL sed_bedload_sd(ised,i) 
            ENDIF    

          ! BM Wu and Lin (2014), adapted from Sed2d
          ELSEIF (bedload == 4) THEN
            IF(iSedtype(ised).EQ.0) THEN
              !* MUD-like sediment type, no bedload
              FX_r(i) = 0.0d0
              FY_r(i) = 0.0d0
            ELSEIF(iSedtype(ised).EQ.1) THEN
              !* SAND-like sediment 
              CALL sed_bedload_wl14(ised,i)
            ENDIF

          ELSE
            CALL parallel_abort('SED3D: incorrect bedload option')
          ENDIF

!---------------------------------------------------------------------
! - Apply bedload transport rate coefficient (nondimensional).
! jl. Added bedload transport as implemented in ROMS 
!---------------------------------------------------------------------

          FX_r(i) = FX_r(i)*bedload_coeff ![m^2]
          FY_r(i) = FY_r(i)*bedload_coeff

!---------------------------------------------------------------------
! - Diffusive fluxes in flow direction (Fortunato et al., 2009; eq 2)
!---------------------------------------------------------------------

          FX_r(i) = FX_r(i)+bdldiffu*(1.0-bed(top,i,iporo))*DABS(FX_r(i))*dzdx
          FY_r(i) = FY_r(i)+bdldiffu*(1.0-bed(top,i,iporo))*DABS(FY_r(i))*dzdy


!---------------------------------------------------------------------
! - Consistency check
!---------------------------------------------------------------------
          
          IF (FX_r(i)/=FX_r(i)) THEN
            WRITE(errmsg,*)'FX_r is NaN',myrank,i,FX_r(i)
            CALL parallel_abort(errmsg)
          ENDIF
          IF (FY_r(i)/=FY_r(i)) THEN
            WRITE(errmsg,*)'FY_r is NaN',myrank,i,FY_r(i)
            CALL parallel_abort(errmsg)
          ENDIF

!---------------------------------------------------------------------
! - END Loop over elements
!---------------------------------------------------------------------
        ENDDO !End loop i=1,nea

        IF (bedload_filter == 1) THEN
          !BM test bed filter from sed2d on bedload fluxes
          ! X
          CALL sed2d_filter_diffu(FX_r(:),tmp_e,nea)
          FX_r(:) = tmp_e
          CALL exchange_e2d(FX_r(:))
          ! Y
          CALL sed2d_filter_diffu(FY_r(:),tmp_e,nea)
          FY_r(:) = tmp_e
          CALL exchange_e2d(FY_r(:))
        END IF

!---------------------------------------------------------------------
! - Add bedload fluxes caused by acceleration skewness of waves
!---------------------------------------------------------------------
#ifdef USE_WWM
        IF (bedload_acc > 0) THEN
          FX_r(:) = FX_r(:) + Qaccu(:) ! [m2]
          FY_r(:) = FY_r(:) + Qaccv(:)

          ! Save the part of bedload fluxes caused by wave acceleration
          ! before the flux limitation linked to active layer
          r_accu(:) = Qaccu(:)/FX_r(:)
          r_accv(:) = Qaccv(:)/FY_r(:)
        END IF
#endif


!---------------------------------------------------------------------
! - Apply a limiter on bedload fluxes based on the active layer
!   criterion --> modification of FX_r and FY_r according to the 
!   sediment mass available in the active layer
!---------------------------------------------------------------------
        IF (bedload_limiter == 1) THEN
          ! Modify FX_r/FY_r (same unit: m2)
          CALL bedload_flux_limiter
        END IF
   
!---------------------------------------------------------------------
! - Adjust FX_r/FY_r of sediment class ised according to its mass
!   fraction within the bed mixture.
!   Previously done in subroutine bedchange_bedload, but inconsistent
!   in case of sed_morph==0.
!   Note that in case of bedload_acc>0, bed_frac is also applied to 
!   to Qaccu/Qaccv, already included in FX_r/FY_r.
!---------------------------------------------------------------------

        FX_r(:) = FX_r(:)*bed_frac(top,:,ised)
        FY_r(:) = FY_r(:)*bed_frac(top,:,ised)


!---------------------------------------------------------------------
! - Compute morphology/bed change (characteristics and dh) due to bed load
!---------------------------------------------------------------------

        IF (sed_morph.GE.1) THEN
          IF(myrank.EQ.0) WRITE(16,*)'SED: Entering bedchange_bedload:',ised,it
          CALL bedchange_bedload(ised,it,moitn,mxitn,rtol,qsan,      &
          &                      hbed,hbed_ised)
          IF(myrank.EQ.0) WRITE(16,*)'SED: leaving bedchange_bedload:',ised,it
        ELSE !0326a
          bed_mass(:,:,nnew,:)=bed_mass(:,:,nstp,:)
        ENDIF

!-----------------------------------------------------------------------
!  Output bedload fluxes.
!-----------------------------------------------------------------------
        DO i=1,np
          IF(idry(i)==1) CYCLE
          ks = 0 !total count
          cff1=0.0d0
          cff2=0.0d0
          DO j=1,nne(i)
            k = indel(j,i)
            IF(idry_e(k)==1) CYCLE
            ks = ks+1
            ! compute bedload transport rate in [kg/m/s]
            bedldu(i,ised) = bedldu(i,ised) + FX_r(k)*Srho(ised)/dt
            bedldv(i,ised) = bedldv(i,ised) + FY_r(k)*Srho(ised)/dt
            ! bedload flux caused by wave acceleration [kg/m/s]
            cff1 = cff1 + r_accu(k)*FX_r(k)*Srho(ised)/dt
            cff2 = cff2 + r_accv(k)*FY_r(k)*Srho(ised)/dt
          ENDDO ! End loop nne

          if(ks==0) call parallel_abort('SED3D: (4)')
          bedldu(i,ised) = bedldu(i,ised)/ks
          bedldv(i,ised) = bedldv(i,ised)/ks
          ! Integration over all sediment classes for Qacc
          Qaccun(i) = Qaccun(i)+cff1/ks
          Qaccvn(i) = Qaccvn(i)+cff2/ks

        ENDDO !End loop np

        ! Exchange ghosts
        CALL exchange_p2d(bedldu(:,ised))
        CALL exchange_p2d(bedldv(:,ised))

        CALL exchange_p2d(Qaccun(:))
        CALL exchange_p2d(Qaccvn(:))

!---------------------------------------------------------------------
! - END Loop over all non-cohesive sediment classes
!---------------------------------------------------------------------
      ENDDO ! ised=1,ntr_l


!---------------------------------------------------------------------
! - Update mean bed surface properties.
! Sd50 must be positive definite, due to BBL routines.
! Srho must be >1000, due to (s-1) in BBL routines.
!---------------------------------------------------------------------

      DO i=1,nea
! * FG. Here, I removed testing for dry element, because bed
!       characteristics have to be updated over the whole domain
!       IF (idry_e(i)==1) CYCLE

        ! update sediment fractions
        cff3 = 0.0d0
        DO ised=1,ntr_l
          cff3 = cff3+bed_mass(top,i,nnew,ised)
        ENDDO
        ! Test to prevent div by 0
        cff3=max(cff3,eps)

        DO ised=1,ntr_l
          bed_frac(top,i,ised) = bed_mass(top,i,nnew,ised)/cff3 !\in [0,1]
        ENDDO

        ! BM : update of the top layer porosity according to changes in
        ! sediment fractions
        IF (poro_option .EQ. 2) THEN
          CALL sed_comp_poro_noncoh(bed_frac(top,i,:),bed(top,i,iporo))
        END IF

        ! Weighted geometric mean 
        cff1 = 1.0d0
        cff2 = 1.0d0
        cff3 = 1.0d0
        cff4 = 1.0d0
        cff5 = 0.0d0
        DO ised=1,ntr_l
          cff1 = cff1*tau_ce(ised)**bed_frac(top,i,ised)
          cff2 = cff2*Sd50(ised)**bed_frac(top,i,ised)
          cff3 = cff3*(Wsed(ised)+eps)**bed_frac(top,i,ised)
          cff4 = cff4*Srho(ised)**bed_frac(top,i,ised)
          cff5 = cff5+bed_frac(top,i,ised)
        ENDDO !End loop ntracer

        ! Update mean bottom properties
        if(cff5<0) then
          WRITE(errmsg,*)'SED3D: cff5<0 (2); ',cff5
          call parallel_abort(errmsg)
        else if(cff5==0) then
          !WRITE(12,*)'SED3D: all eroded at elem. ',ielg(i),it
          !Care takers
          bottom(i,itauc) = tau_ce(1)
          bottom(i,isd50) = Sd50(1)
          bottom(i,iwsed) = Wsed(1)+eps
          bottom(i,idens) = Srho(1)
        else !cff5>0
          bottom(i,itauc) = cff1**(1.0d0/cff5)
!          bottom(i,isd50) = MIN(cff2,Zob(i))
          bottom(i,isd50) = MAX(MIN(cff2**(1.0d0/cff5),MAXVAL(Sd50(:))), &
     &MINVAL(Sd50(:)))
          bottom(i,iwsed) = cff3**(1.0d0/cff5)
          bottom(i,idens) = MAX(cff4**(1.0d0/cff5),1050.0d0)
        endif !cff5

! * FG. I don't understand why the grain diameter is the minimum
!       between the grain diameter and the roughness length.
!       Roughness length depends on the grain size, but
!       the grain size does not depends on the roughness...
!       However, total D50 should not be lower than the lowest
!       D50(ised) and higher than the highest D50(ised)
! * FG. I added the sum of weights (bed_frac) for the computation
!       of geometric mean, because the rounding of bed fractions
!       leads to total bed fraction slightly different than 1.0,
!       inducing significant errors for bed properties

      ENDDO !i=1,nea
!---------------------------------------------------------------------
!  End computing bedload sediment transport.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!           *** Compute suspended sediment transport. ***
!
!---------------------------------------------------------------------
      IF (suspended_load == 1) THEN
        IF(myrank.EQ.0) WRITE(16,*)'SED: Entering suspended load...'

        IF(ised_bc_bot==2) call sed_pickup(im_pick_up) !Tsinghua group
        SED_LOOP: DO ised=1,ntr_l
          indx=isand(ised) !into 1:ntracers
          DO i=1,nea
            IF (idry_e(i)==1) CYCLE

            select case(ised_bc_bot)
              case(1) !Warner
                !Depositional mass (kg/m/m) in a time step
!                cff=ze(nvrt,i)-ze(kbe(i),i) !total depth
!                cff1=ze(kbe(i)+1,i)-ze(kbe(i),i) !bottom cell thickness
!
!                if(cff1.LE.0.d0) then
!                  WRITE(errmsg,*)'SED, wrong bottom layer0:',cff,ze(kbe(i)+1,i),ze(kbe(i),i)
!                  CALL parallel_abort(errmsg)
!                endif
!
!                !Limit ratio between reference depth and bottom depth
!                !ta=min(0.5d0,relath*cff/cff1) !usually the relath is 0.01; for natural river, it should be even smaller
!
!                !aref=ta*we(kbe(i)+1,i) !w-vel. at ref. height                
!                !Estimate the starting pt
!                cff8=(ze(kbe(i)+1,i)+ze(kbe(i),i))/2 !bottom half cell -starting pt of ELM
!                cff9=cff8+Wsed(ised)*dt !>cff8; foot of char.
!                if(cff9<=cff8) call parallel_abort('SED: (9)')
!                Ksed=nvrt+1 !init. for abnormal case
!                do k=kbe(i)+2,nvrt
!                  if(cff9<=(ze(k,i)+ze(k-1,i))/2) then
!                    Ksed=k
!                    exit
!                  endif
!                enddo !k
!
!                depo_mss=0 !(min(cff9,cff8)-ze(kbe(i),i))*tr_el(indx,kbe(i)+1,i)
!                do k=kbe(i)+2,min(nvrt,Ksed)
!                  cff7=(ze(k-2,i)+ze(k-1,i))/2
!                  cff8=(ze(k,i)+ze(k-1,i))/2
!                  if(cff9<cff7) then
!                    WRITE(errmsg,*)'SED, wrong cff9:',cff9,cff7,ielg(i),ised,k
!                    CALL parallel_abort(errmsg)
!                  endif
!
!                  if(cff9<cff8) then
!                    cff5=(cff9-cff7)/(cff8-cff7) !ratio
!                    if(cff5<0.or.cff5>1) then
!                      WRITE(errmsg,*)'SED, wrong ratio:',cff5,ielg(i),ised,k 
!                      CALL parallel_abort(errmsg)
!                    endif
!                    cff6=(1-cff5)*tr_el(indx,k-1,i)+cff5*tr_el(indx,k,i) !conc @upper layer
!                  else
!                    cff6=tr_el(indx,k,i)
!                  endif
!
!                  depo_mss=depo_mss+(cff6+tr_el(indx,k-1,i))/2*(min(cff9,cff8)-cff7) !>=0
!                  !depo_mss=depo_mss+cff6*(min(cff9,cff8)-cff7) !>=0 (lower depos. flux)
!                enddo !k
!
!                !Deal with above F.S. case
!                if(Ksed==nvrt+1) then
!                  cff7=(ze(nvrt,i)+ze(nvrt-1,i))/2
!                  depo_mss=depo_mss+(min(ze(nvrt,i),cff9)-cff7)*tr_el(indx,nvrt,i)
!                endif

                !BM: Correction of depo_mss, otherwise conservativity
                !    issues occur 
                depo_mss=dt*Wsed(ised)*tr_el(indx,kbe(i)+1,i) !s * m/s * kg/m3

                !BM: output
                depflxel(i,ised)=depo_mss/dt ! en kg/m2/s

! - Compute erosion, eros_mss (kg/m/m) following 
!  (original erosion flux is in kg/m/m/s; note dt below)
                cff1=(1.0d0-bed(top,i,iporo))*bed_frac(top,i,ised)
                if(ierosion==0) then 
                  !Ariathurai & Arulanandan (1978)
                  eros_mss=MAX(0.0d0,dt*Erate(ised)*cff1*(tau_wc(i)/tau_ce(ised)-1.0d0)) !kg/m/m
                else if(ierosion==1) then 
                  !Winterwerp et al. (JGR, vol 117, 2012, C10020); eq. (15)
                  cff2=tau_wc(i)/tau_ce(ised) !tau_b/tau_cr [-]
                  if(cff2<0.52) then
                    cff4=0
                  else if(cff2<=1.7) then
                    cff4=-0.144*cff2**3+0.904*cff2*cff2-0.823*cff2+0.204 ![-]
                  else
                    cff4=cff2-1 ![-]
                  endif !cff2
                  cff3=tau_ce(ised)*rhom ![Pa]; tau_ce: critical shear stress in m^2/s/s
                  eros_mss=Erate(ised)*cff3*cff4 !kg/m/m/s; Erate (M_E in original paper) in [s/m]
                  eros_mss=MAX(0.0d0,cff1*dt*eros_mss) !kg/m/m
                else
                  CALL parallel_abort('SED3D: unknown erosion formula')
                endif !ierosion
                eros_mss=MIN(eros_mss,MIN(Srho(ised)*cff1*bottom(i,iactv),bed_mass(top,i,nnew,ised))+depo_mss) !>=0

!...            Save erosion and depo. fluxes [kg/m/m/s] for b.c. of transport eq.
!               This should not be scaled by morph_fac, b/cos both mass
!               and dt are multiplied by same factor
                flux_eros=eros_mss/dt !>=0
                !flux_depo=depo_mss/dt

                !BM: output
                eroflxel(i,ised)=flux_eros ! en kg/m2/s            



!...            Update flx_bt for transport solver
                flx_bt(indx,i)=-flux_eros ![kg/m/m/s]

                IF (sed_morph>=1) THEN
! - Apply morphology factor to flux and settling...
                  eros_mss=eros_mss*morph_fac
                  depo_mss=depo_mss*morph_fac

! - Depth change due to erosion/deposition of suspended sediment
                  if(bed(top,i,iporo)==1) then
                    WRITE(errmsg,*)'SED3D: bed(top,i,iporo)==1; ',top,i,iporo
                    CALL parallel_abort(errmsg)
                  endif
                  hdep(i)=hdep(i)+(eros_mss-depo_mss)/Srho(ised)/(1.0d0-bed(top,i,iporo))
                ENDIF !sed_morph

! - If first time step of deposit, then store deposit material in
! temporary array, dep_mass.
                IF (eros_mss<depo_mss) THEN
                  !IF(time>bed(top,i,iaged)+1.1d0*dt.and.bed(top,i,ithck)>newlayer_thick) THEN
                  IF(bed(top,i,ithck)>newlayer_thick) THEN !171027
                    dep_mass(i,ised)=depo_mss-eros_mss !>0 !kg/m/m
                  ENDIF
                  bed(top,i,iaged)=time
                ENDIF

! - Update bed mass arrays.
!jl.  The whole nnew/bnew,nstp is way more confusing than need be
                bed_mass(top,i,nnew,ised)=MAX(bed_mass(top,i,nnew,ised)-eros_mss+depo_mss,0.0d0)
                DO k=2,Nbed
                  bed_mass(k,i,nnew,ised)=bed_mass(k,i,nstp,ised)
                  if(bed_mass(k,i,nnew,ised)<0) then
                    WRITE(errmsg,*)'SED3D: bed_m<0',k,bed_mass(k,i,nnew,ised)
                    CALL parallel_abort(errmsg)
                  endif
                ENDDO
                !bed_mass(top,i,nnew,ised)=MAX(bed_mass(top,i,nnew,ised),0.0d0) 

              case(2) !Tsinghua Univ. group
                !Xiaonan's addition here
                !Calculate flx_bt=D-E-w_s*T_{kbe+1}, and several other variables (e.g. bed_mass)
                ! - Compute erosion, eros_mss (kg/m/m) following Ariathurai and Arulanandan (1978)
                cff1=(1.0d0-bed(top,i,iporo))*bed_frac(top,i,ised)
                eros_mss=MAX(0.0d0,dt*sedcaty(i,ised)*Srho(ised)*cff1) !kg/m/m  !sedcaty(i,ised) m/s

!Tsinghua group---------------------------------------------------
                if (itur==5) then !1120:itur==5
!                  depo_mss=alphd*dt*Wsed(ised)*tr_el(indx,kbe(i)+1,i)*phai_m(kbe(i)+1,indx,i) !1120:close
                  depo_mss=alphd*dt*Wsed(ised)*tr_el(indx,kbe(i)+1,i)*sum(Phai(kbe(i)+1,ised,elnode(1:i34(i),i)))/i34(i)
                else
                  depo_mss=alphd*dt*Wsed(ised)*tr_el(indx,kbe(i)+1,i)    !dt*Wsed(ised)*tr_el(indx,kbe(i)+1,i)
                endif !itur
!Tsinghua group---------------------------------------------------
                eros_mss=MIN(eros_mss,MIN(Srho(ised)*cff1*bottom(i,iactv),bed_mass(top,i,nnew,ised))+depo_mss) !>=0

                !Save erosion and depo. fluxes [kg/m/m/s] for b.c. of transport eq.
                flux_eros= eros_mss/dt !>=0
                !Update flx_bt for transport solver
                flx_bt(indx,i)=-flux_eros ![kg/m/m/s]

                if(flx_bt(indx,i)/=flx_bt(indx,i)) then
                  WRITE(errmsg,*)'flx_bt: has nan; ',indx,ielg(i),eros_mss,depo_mss,sedcaty(i,ised)
                  CALL parallel_abort(errmsg)
                endif

                IF (sed_morph>=1) THEN
                  ! - Apply morphology factor to flux and settling...
                  eros_mss=eros_mss*morph_fac
                  depo_mss=depo_mss*morph_fac

                  ! - Depth change due to erosion/deposition of suspended sediment
                  if(bed(top,i,iporo)==1) then
                    WRITE(errmsg,*)'SED3D: bed(top,i,iporo)==1; ',top,i,iporo
                    CALL parallel_abort(errmsg)
                  endif
                  hdep(i)=hdep(i)+(eros_mss-depo_mss)/Srho(ised)/(1.0d0-bed(top,i,iporo))
                ENDIF !sed_morph

                ! - If first time step of deposit, then store deposit material in
                ! temporary array, dep_mass.
                IF (eros_mss<depo_mss) THEN
                  IF(time>bed(top,i,iaged)+1.1d0*dt.and.bed(top,i,ithck)>newlayer_thick) THEN
                    dep_mass(i,ised)=depo_mss-eros_mss !>0 !kg/m/m
                  ENDIF
                  bed(top,i,iaged)=time
                ENDIF

                ! - Update bed mass arrays.
                !jl.  The whole nnew/bnew,nstp is way more confusing than need be
                bed_mass(top,i,nnew,ised)=MAX(bed_mass(top,i,nnew,ised)-eros_mss+depo_mss,0.0d0)
                DO k=2,Nbed
                  bed_mass(k,i,nnew,ised)=bed_mass(k,i,nstp,ised)
                  if(bed_mass(k,i,nnew,ised)<0) then
                    WRITE(errmsg,*)'SED3D: bed_m<0',k,bed_mass(k,i,nnew,ised)
                    CALL parallel_abort(errmsg)
                  endif
                ENDDO
                !bed_mass(top,i,nnew,ised)=MAX(bed_mass(top,i,nnew,ised),0.0d0) 

              case default
                call parallel_abort('SED: unknown ised_bc_bot')
            end select
          ENDDO !i=1,nea

!          IF (sed_morph.GE.1) THEN
!            CALL exchange_e2d(hdep(:))
!          ENDIF

        ENDDO SED_LOOP !ised=1,ntr_l

        IF(myrank==0) WRITE(16,*)'SED: out of sus. load loop...'

!        call schism_output_custom(indx,5,1,222,'hdep',1,nea,hdep)

! - If first time step of deposit, create new layer and combine bottom
! two bed layers.
!YJZ: the splitting/combining of layers may cause noise. Using a large
!newlayer_thick would avoid this part.

        DO i=1,nea
          IF (idry_e(i)==1) CYCLE

          IF(Nbed>1) THEN
            cff = 0.0d0
            DO ised=1,ntr_l
              cff = cff+dep_mass(i,ised)
            ENDDO

            IF (cff>0.0d0) THEN !deposition ocurred

              ! Combine bottom layers.
              ! BM : poro updated after
              !bed(Nbed,i,iporo) = 0.5d0*(bed(Nbed-1,i,iporo)+        &
              !&                          bed(Nbed,i,iporo))
              bed(Nbed,i,iaged) = 0.5d0*(bed(Nbed-1,i,iaged)+        &
              &                          bed(Nbed,i,iaged))

              DO ised=1,ntr_l
                bed_mass(Nbed,i,nnew,ised) =                         &
                &                   bed_mass(Nbed-1,i,nnew,ised)+    &
                &                   bed_mass(Nbed,i,nnew,ised)
              ENDDO

              ! Push layers down.
              DO k=Nbed-1,2,-1
                !bed(k,i,iporo) = bed(k-1,i,iporo) ! BM : idem
                bed(k,i,iaged) = bed(k-1,i,iaged)
                DO ised =1,ntr_l
                  bed_mass(k,i,nnew,ised) = bed_mass(k-1,i,nnew,ised)
                ENDDO
              ENDDO

              ! Set new top layer parameters.
              DO ised=1,ntr_l
                if(dep_mass(i,ised)<0) then
                  WRITE(errmsg,*)'SED3D: dep_mass(i,ised)<0; ',dep_mass(i,ised),i,ised
                  CALL parallel_abort(errmsg)
                endif
                bed_mass(2,i,nnew,ised)=MAX(bed_mass(2,i,nnew,ised)-dep_mass(i,ised),0.0d0) !>=0
                bed_mass(top,i,nnew,ised)=dep_mass(i,ised) !>=0
              ENDDO !ised

            ENDIF ! Deposition occurred (cff>0)
          ENDIF !Nbed>1

! If Nbed = 1
! - Recalculate thickness and fractions for all layers.
!   bed(ithck) = bed_mass/Srho
!     [m]      = [m.m-2] / [kg.m-2] --> Do not have to integrate over
!     area
          DO k=1,Nbed
            cff3=0.0d0 !total mass
            DO ised=1,ntr_l
              if(bed_mass(k,i,nnew,ised)<0) then
                WRITE(errmsg,*)'SED3D, mass<0(3):',bed_mass(k,i,nnew,ised),k,ielg(i),ised 
                CALL parallel_abort(errmsg)
              endif
              cff3 = cff3+bed_mass(k,i,nnew,ised)
            ENDDO

            cff3=max(cff3,eps)


!!! BM : update of bed_frac --> compute new porosity --> deduce new thickness

!            bed(k,i,ithck) = 0.0d0
!            if(bed(k,i,iporo)==1) call parallel_abort('SED3D: div. by 0 (6)')
!!'
!            DO ised=1,ntr_l
!              bed_frac(k,i,ised)=bed_mass(k,i,nnew,ised)/cff3
!              bed(k,i,ithck)=bed(k,i,ithck)+bed_mass(k,i,nnew,ised)/(Srho(ised)*(1.0d0-bed(k,i,iporo)))
!            ENDDO !ised

            DO ised=1,ntr_l
              bed_frac(k,i,ised)=bed_mass(k,i,nnew,ised)/cff3
            ENDDO !ised
            IF (poro_option .EQ. 2) THEN
              CALL sed_comp_poro_noncoh(bed_frac(k,i,:),bed(k,i,iporo))
            END IF
            if(bed(k,i,iporo)==1) call parallel_abort('SED3D: div. by 0(6)')
!'
            bed(k,i,ithck) = 0.0d0
            DO ised=1,ntr_l
              bed(k,i,ithck)=bed(k,i,ithck)+bed_mass(k,i,nnew,ised)/(Srho(ised)*(1.0d0-bed(k,i,iporo)))
            ENDDO !ised
          ENDDO !k=1,Nbed

          ! BM: bottom(i,isd50) and  bottom(i,tauc) need to be updated 
          ! in the top layer to compute an active layer thickness bottom(i,iactv) 
          ! representative of new bed features, resulting from previous bedload 
          ! and suspended load dynamics.
          cff1 = 1.0d0
          cff2 = 1.0d0
          cff5 = 0.0d0
          DO ised=1,ntr_l
            cff1 = cff1*tau_ce(ised)**bed_frac(1,i,ised)
            cff2 = cff2*Sd50(ised)**bed_frac(1,i,ised)
            cff5 = cff5+bed_frac(top,i,ised)
          ENDDO !ntr_l
          if(cff5<0) then
            WRITE(errmsg,*)'SED3D: cff5<0 (1a); ',cff5
            call parallel_abort(errmsg)
          else if(cff5==0) then
            bottom(i,itauc) = tau_ce(1)
            bottom(i,isd50) = Sd50(1)
          else !cff5>0
            bottom(i,itauc) = cff1**(1.0d0/cff5)
            bottom(i,isd50) = MAX(MIN(cff2**(1.0d0/cff5),MAXVAL(Sd50(:))),MINVAL(Sd50(:)))
          endif !cff5

        ENDDO !i=1,nea

        IF(myrank==0) WRITE(16,*)'SED: End of suspended sediment'
      ENDIF !suspended_load==1
!---------------------------------------------------------------------
!         *** End of Suspended Sediment only section ***
!---------------------------------------------------------------------

!---------------------------------------------------------------------
! * Application of the bed_mass filter - not quite working so turn off
!---------------------------------------------------------------------
!     IF(bedmass_filter.GE.1) THEN
!       IF(myrank==0) WRITE(16,*)'SED: doing sed_bedmass_filter '
!       CALL sed_bedmass_filter(hdep)
!       IF(myrank==0) WRITE(16,*)'SED: done sed_bedmass_filter'
!     ENDIF

     IF(myrank==0) WRITE(16,*)'SED: adjusting bed layers...'

!---------------------------------------------------------------------
! - Ensure top bed layer thickness is greater or equal than active 
! layer thickness. If need to add sed to top layer, then entrain from 
! lower levels. Create new layers at bottom to maintain Nbed.
!---------------------------------------------------------------------

      DO i=1,nea
        IF (idry_e(i)==1) CYCLE

        ! - Computation of the active layer thickness, bottom(i,j,iactv),
        !   based on the relation of Harris and Wiberg (1997)
        bottom(i,iactv)=MAX(0.0d0,7.0d-3*(tau_wc(i)-bottom(i,itauc))*rhom)+6.0d0*bottom(i,isd50)
        ! - BM: the active layer thickness cannot exceed a used-defined
        !       one, actv_max
        bottom(i,iactv)=MIN(actv_max,bottom(i,iactv))
!>0



! - Apply morphology factor
!jl. The application of morph_fac is arbitrary here.  This is not in  
!    need of immediate attention, but if morph_facs differ between  
!    sediment classes there should be some attention, but if  
!    morph_facs differ between sediment classes there should be some 
!    kind of averaging, or other method, to determine the morph_fac to
!    use.

!BM. I don't think that it makes sense to have different morph_fac between
!    sediment classes ... morph_fac is now set as a parameter in
!    sediment.in

        IF(sed_morph>=1) bottom(i,iactv)=MAX(bottom(i,iactv)*morph_fac,bottom(i,iactv))

        IF(bed(top,i,ithck)<bottom(i,iactv)) THEN

          IF(Nbed==1) THEN
            bottom(i,iactv) = bed(top,i,ithck) !possibly 0
          ELSE ! Nbded>1
            !Increase top bed layer - entrain from layers below
            thck_to_add = bottom(i,iactv)-bed(top,i,ithck) !>0
            thck_avail = 0.0d0
            Ksed=Nbed !initialize to include all
            DO k=2,Nbed
              thck_avail=thck_avail+bed(k,i,ithck)
              IF (thck_avail>=thck_to_add) THEN
                Ksed=k; exit
              ENDIF
            ENDDO !k

            IF(thck_avail<0) THEN
              write(errmsg,*)'SED3D: thck_avail<0:',ielg(i),thck_avail
              CALL parallel_abort(errmsg)
            ELSE IF(thck_avail==0) THEN
              IF (sed_debug .EQ. 1) THEN
                write(12,*)'SED3D: not enough sed; likely all eroded:', &
       &ielg(i),thck_avail,thck_to_add,bed(:,i,ithck),it
              END IF
              bottom(i,iactv)=bed(top,i,ithck) !0326
              bed(2:Nbed,i,ithck)=0
              bed_frac(2:Nbed,i,:)=0
              bed_mass(2:Nbed,i,nnew,:)=0

!              bed(:,i,ithck)=0
!              bed_frac(:,i,:)=0
!              bed_mass(:,i,nnew,:)=0
            ELSE !thck_avail>0
! - Catch here if there was not enough bed material
              IF(thck_avail<thck_to_add) THEN
                bottom(i,iactv) = bed(top,i,ithck)+thck_avail
                thck_to_add = thck_avail
              ENDIF

              if(thck_avail<thck_to_add) then
                write(errmsg,*)'SED3D, thck_avail<thck_to_add:',ielg(i),thck_avail,thck_to_add
                CALL parallel_abort(errmsg)
              endif

! - Update bed mass of top layer and fraction of layer Ksed
              cff2=(thck_avail-thck_to_add)/MAX(bed(Ksed,i,ithck),eps) !\in [0,1]; fraction of Ksed left
              DO ised=1,ntr_l
                cff1=0.0d0
                DO k=1,Ksed
                  cff1 = cff1+bed_mass(k,i,nnew,ised)
                ENDDO
                cff3 = cff2*bed_mass(Ksed,i,nnew,ised) !what's left there
                if(cff1-cff3<0) then
                  write(errmsg,*)'SED3D, cff1-cff3<0:',ielg(i),cff1-cff3
                  CALL parallel_abort(errmsg)
                endif
                bed_mass(top,i,nnew,ised)=cff1-cff3 !>=0
                bed_mass(Ksed,i,nnew,ised)=cff3 !>=0
              ENDDO !ised

! - Update thickness of fractional layer Ksed
              bed(Ksed,i,ithck)=thck_avail-thck_to_add !>=0

! - Upate bed fraction of top layer
              cff3=0.0d0
              DO ised=1,ntr_l
                cff3=cff3+bed_mass(top,i,nnew,ised)
              ENDDO
              cff3=max(cff3,eps)
              DO ised=1,ntr_l
                bed_frac(top,i,ised)=bed_mass(top,i,nnew,ised)/cff3 !>=0
              ENDDO

! Upate bed thickness of top layer
              bed(top,i,ithck)=bottom(i,iactv)
              !!! >BM : new porosity --> Hyp : Top layer can be different
              !!!      from theoritical active layer thickness
              IF (poro_option .EQ. 2) THEN
                CALL sed_comp_poro_noncoh(bed_frac(top,i,:),bed(top,i,iporo))
              END IF
              IF (bed(top,i,iporo)==1) call parallel_abort('SED3D: div. by 0(6)')
!'
              bed(top,i,ithck) = 0.0d0
              DO ised=1,ntr_l
                bed(top,i,ithck)=bed(top,i,ithck)+bed_mass(top,i,nnew,ised)/(Srho(ised)*(1.0d0-bed(top,i,iporo)))
              ENDDO !ised
              !!! <BM


            
! Pull all layers closer to the surface 
              ks=Ksed-2 !reduction of # of layers; >=0
              DO k=Ksed,Nbed
                bed(k-ks,i,ithck) = bed(k,i,ithck)
                bed(k-ks,i,iporo) = bed(k,i,iporo)
                bed(k-ks,i,iaged) = bed(k,i,iaged)
                DO ised=1,ntr_l
                  bed_frac(k-ks,i,ised) = bed_frac(k,i,ised)
                  bed_mass(k-ks,i,nnew,ised) = bed_mass(k,i,nnew,ised)
                ENDDO
              ENDDO !k

! By now there are only Nbed-ks layers left
! Split bottom-most layer to make Nbed layers in total and
! conserve total mass
              if(ks+1==0) call parallel_abort('SED3D: ks+1==0')
              cff=1.0d0/(ks+1) !fraction
              DO k=Nbed,Nbed-ks,-1
                bed(k,i,ithck) = bed(Nbed-ks,i,ithck)*cff
                bed(k,i,iaged) = bed(Nbed-ks,i,iaged)
                !!! BM porosity 
                bed(k,i,iporo) = bed(Nbed-ks,i,iporo)
                DO ised=1,ntr_l
                  bed_frac(k,i,ised) = bed_frac(Nbed-ks,i,ised)
                  bed_mass(k,i,nnew,ised)=bed_mass(Nbed-ks,i,nnew,ised)*cff
                ENDDO
              ENDDO ! k
            ENDIF !thck_avail

          ENDIF !Nbed > 1
        ENDIF  ! End test increase top bed layer
      ENDDO !i=1,nea

!---------------------------------------------------------------------
! Update mean surface properties.
! write(*,*)'update mean surface prop' 
! Sd50 must be positive definite, due to BBL routines.
! Srho must be >1000, due to (s-1) in BBL routines
!---------------------------------------------------------------------

! * FG. Here, I removed testing for dry element, because bed
!       characteristics have to be updated over the whole domain
      DO i=1,nea
        cff1 = 1.0d0
        cff2 = 1.0d0
        cff3 = 1.0d0
        cff4 = 1.0d0
        cff5 = 0.0d0
        ! weighted geometric mean
        DO ised=1,ntr_l
          cff1 = cff1*tau_ce(ised)**bed_frac(1,i,ised)
          cff2 = cff2*Sd50(ised)**bed_frac(1,i,ised)
          cff3 = cff3*(Wsed(ised)+eps)**bed_frac(1,i,ised)
          cff4 = cff4*Srho(ised)**bed_frac(1,i,ised)
          cff5 = cff5+bed_frac(top,i,ised)
        ENDDO !ntr_l

        if(cff5<0) then
          WRITE(errmsg,*)'SED3D: cff5<0 (1); ',cff5
          call parallel_abort(errmsg)
        else if(cff5==0) then
          !WRITE(12,*)'SED3D: all eroded at elem. (2):',ielg(i),it
          !Caretakers
          bottom(i,itauc) = tau_ce(1)
          bottom(i,isd50) = Sd50(1)
          bottom(i,iwsed) = Wsed(1)+eps
          bottom(i,idens) = Srho(1)
        else !cff5>0
          bottom(i,itauc) = cff1**(1.0d0/cff5)
!          bottom(i,isd50) = MIN(cff2,Zob(i))
          bottom(i,isd50) = MAX(MIN(cff2**(1.0d0/cff5),MAXVAL(Sd50(:))),MINVAL(Sd50(:)))
          bottom(i,iwsed) = cff3**(1.0d0/cff5)
          bottom(i,idens) = MAX(cff4**(1.0d0/cff5),1050.0d0)
        endif !cff5

! * FG. I don't understand why the grain diameter is the minimum
!       between the grain diameter and the roughness length.
!       Roughness length depends on the grain size, but
!       the grain size does not depends on the roughness...
!       However, total D50 should not be lower than the lowest
!       D50(ised) and higher than the highest D50(ised)
! * FG. I added the sum of weights (bed_frac) for the computation
!       of geometric mean, because the rounding of bed fractions
!       leads to total bed fraction slightly different than 1.0,
!       inducing significant errors for bed properties
      ENDDO !i=1,nea

!---------------------------------------------------------------------
! Convert bed sediment properties from elements to node for outputs
!---------------------------------------------------------------------
      IF(myrank==0) WRITE(16,*)'SED: converting arrays to nodes...'
      
      !First, properties which cannot be equal to zero even if dry
      bed_d50n   = 0.0d0
      bed_fracn  = 0.0d0

      !BM: new outputs
      eroflxn = 0.0d0
      depflxn = 0.0d0
      poron = 0.0d0

      DO i = 1,np
        ta = 0.0d0
        DO j = 1,nne(i)
          ie = indel(j,i)
          ta = ta + area(ie)
          bed_d50n(i)=bed_d50n(i)+bottom(ie,isd50)*area(ie)
          poron(i)=poron(i)+bed(top,ie,iporo)*area(ie)   
          DO ised=1,ntr_l
            bed_fracn(i,ised)=bed_fracn(i,ised)+bed_frac(1,ie,ised)*area(ie)
            eroflxn(i)=eroflxn(i)+eroflxel(ie,ised)*area(ie)
            depflxn(i)=depflxn(i)+depflxel(ie,ised)*area(ie)
            if(bed_frac(1,ie,ised)<0) then
              WRITE(errmsg,*)'SED3D, frac<0 (1):',bed_frac(1,ie,ised),ielg(ie)
              CALL parallel_abort(errmsg)
            endif
          ENDDO !ised
        ENDDO !j
        IF(ta.EQ.0) THEN
          CALL parallel_abort('SEDIMENT: elem2nod (1)')
        ELSE
          bed_d50n(i)=bed_d50n(i)/ta
          poron(i)=poron(i)/ta
          eroflxn(i)=eroflxn(i)/ta
          depflxn(i)=depflxn(i)/ta

          DO ised=1,ntr_l
            bed_fracn(i,ised)=bed_fracn(i,ised)/ta
          ENDDO ! END loop ntr_l
        ENDIF
      ENDDO ! END loop i=1,np
      CALL exchange_p2d(bed_d50n(:))
      CALL exchange_p2d(poron(:)) 
      CALL exchange_p2d(eroflxn(:))
      CALL exchange_p2d(depflxn(:))

      DO  ised=1,ntr_l
        CALL exchange_p2d(bed_fracn(:,ised))
      ENDDO !ised

      !bottom shear stress
      bed_taun  = 0.0d0
      DO i=1,np
        IF(idry(i)==1) CYCLE
        ta=0.0d0
        DO j=1,nne(i)
          ie=indel(j,i)
          IF(idry_e(ie)==0)THEN
            ta=ta+area(ie)
            bed_taun(i)=bed_taun(i)+tau_wc(ie)*area(ie)
          ENDIF
        ENDDO !j
        IF(ta==0)THEN
          CALL parallel_abort('SEDIMENT: elem2nod (2)')
        ELSE
          bed_taun(i)=bed_taun(i)/ta
        ENDIF
      ENDDO !i
      CALL exchange_p2d(bed_taun(:))

!---------------------------------------------------------------------
! - BCG: If only one bed layer, re-initialize bed thickness to initial
! thickness and update bed_mass accordingly
!---------------------------------------------------------------------
      IF(sed_morph==2.AND.Nbed==1)THEN
        DO i=1,nea
          tmp=sum(bedthick_overall(elnode(1:i34(i),i)))/i34(i)
          bed(1,i,ithck)=tmp !bedthick_overall
          DO ised=1,ntr_l
            bed_mass(1,i,nnew,ised)=tmp*Srho(ised)*(1.0d0-bed(1,i,iporo))*bed_frac(1,i,ised)
            if(bed_mass(1,i,nnew,ised)<0) then
              WRITE(errmsg,*)'SED3D: mass<0; ',bed_mass(1,i,nnew,ised),i,nnew,ised
              CALL parallel_abort(errmsg)
            endif
          ENDDO
        ENDDO
      ENDIF

!---------------------------------------------------------------------
! - Store old bed thickness.
!---------------------------------------------------------------------
!jl.  Do we not need to save the total bed thickness???
!     I didn't see where it is used in ROMS other than here although
!     it gets passed around a lot.
!
!      IF (sed_morph.GE.1) THEN
!          do i=1,nea
!            bed_thick(i,nnew)=0.0_r8
!            do kbed=1,Nbed
!              bed_thick(i,nnew)=bed_thick(i,nnew)+bed(kbed,i,ithck)
!            end do
!          end do


!     Changing depth variation due to erosion/deposition from 
!     elements to nodes
      IF(sed_morph>=1) THEN
        DO i=1,np
          IF(idry(i)==1) CYCLE
  
          tmp = 0
          ta = 0
          DO j=1,nne(i)
            ie=indel(j,i)
            IF(idry_e(ie)==0) THEN
              ta=ta+area(ie)
              tmp=tmp+hdep(ie)*area(ie)
            ENDIF
          ENDDO !j

          IF(ta==0) THEN
            CALL parallel_abort('SEDIMENT: (4)')
          ELSE
            hdep_nd(i)=tmp/ta
          ENDIF
        ENDDO !i

        ! Exchange
        CALL exchange_p2d(hdep_nd)
  
!---------------------------------------------------------------------
! Loop over all nodes to collect depth change due to bedload
! and suspended load:  dhnd=hdep_nd+hbed
!---------------------------------------------------------------------

        dhnd=0.0d0
        DO i=1,np
          IF(suspended_load==1) dhnd(i)=dhnd(i)+hdep_nd(i)
          dhnd(i)=dhnd(i)+hbed(i) !from bedload
        ENDDO
        CALL exchange_p2d(dhnd)

!---------------------------------------------------------------------
! Compute slope avalanching
!---------------------------------------------------------------------
        IF (slope_avalanching.EQ.1)THEN
          IF(myrank.EQ.0) WRITE(16,*)'start sed_avalanching'
          CALL sed_avalanching(dhnd)
          IF(myrank.EQ.0) WRITE(16,*)'done sed_avalanching'
        ENDIF

!---------------------------------------------------------------------
! Update depth according to bedload and suspended fluxes and 
! avalanching ;  update dp()
! Only if sed_morph==1, else, no application of depth changes:
! no morphology or BCG
!---------------------------------------------------------------------
        IF(sed_morph==1.and.time_stamp>=sed_morph_time*86400) THEN
          DO i=1,np
            !dp(i)=dp(i)+dhnd(i)
            dp(i)=dp(i)+dhnd(i)*imnp(i) !BM: add morphological ramp value (-)
!           Impose bare rock limit
            dp(i)=min(dp(i),dp00(i)+bedthick_overall(i))
          ENDDO !i 
!---------------------------------------------------------------------
! Exchange depths of nodes before going back to circulation code 
!---------------------------------------------------------------------
          CALL exchange_p2d(dp)
        ENDIF ! End test sed_morph == 1 (Not for BCG)
      ENDIF ! End sed_morph >=1 (Morpho or BCG)

!     Output to main routine the total bed mass (kg)
      cff=0
      do i=1,ne
        do k=1,Nbed
          do ised=1,ntr_l
            cff=cff+bed_mass(k,i,2,ised)*area(i) !kg
          enddo !ised
        enddo !k
      enddo !i

      call mpi_allreduce(cff,tot_bedmass,1,rtype,MPI_SUM,comm,ierr)

!     Check interface arrays to main routine
!     Additonal arrays: rough_p (sed_friction.F90)
      cff=sum(dp)+sum(flx_bt)+tot_bedmass
      if(cff/=cff) then
        WRITE(errmsg,*)'SED: has nan; ',sum(dp),sum(flx_bt),tot_bedmass
        CALL parallel_abort(errmsg)
      endif

      deallocate(dep_mass)
      first_call=.false.

      end SUBROUTINE sediment
    
!**************************************************************************** tsinghua group
! Numerical integration in the calculation percentage of sediment-carring capacity
!****************************************************************************
      FUNCTION ERF(x)
      use schism_glbl, only : rkind,pi
      implicit none
   
      real(rkind) :: ERF
      real(rkind), intent(in) ::x
      ERF=2.d0/SQRT(pi)*(x-x**3.0/3.d0+x**5.0/10.d0-x**7.d0/42.d0+x**9.d0/216.d0)
      END FUNCTION ERF


!*****************************************************************************
!pick-up flux-Tsinghua
!*****************************************************************************
      subroutine sed_pickup(flag)

      USE sed_mod
      USE schism_glbl, ONLY : rkind,nea,nvrt,kbe,ze,pi,&
                               xnd,ynd,h0,errmsg,ielg,idry_e,dave,&
                               rho0,grav,q2,i34,elnode,tr_el
                               !Tsinghua group:+tr_el,refht,Tbp !1120:-refht,Tbp 
      USE schism_msgp
      
      IMPLICIT NONE
      include 'mpif.h'

      SAVE

      !- Local variables !
      integer, intent(in) :: flag
      INTEGER :: ised,k,i,indx,j,nd
      INTEGER, PARAMETER :: top = 1      ! Top layer of bed
      INTEGER, PARAMETER :: mirror = 16  ! FD for mirror.out
      REAL(rkind),PARAMETER :: nuf = 1.36d-6
      REAL(rkind),PARAMETER :: kf = 0.4d0
      REAL(rkind),PARAMETER :: muf  = 0.3d0
      REAL(rkind),PARAMETER :: Spmax = 0.6d0
      REAL(rkind),PARAMETER :: Bstar = 0.06d0
      REAL(rkind),PARAMETER :: Eta0  = 0.5d0
      REAL(rkind),PARAMETER :: Clift = 0.1d0

      REAL(rkind) :: cff0,cff1,cff2,cff3,cff4,cff5
      REAL(rkind) :: tmp,FAI,s,downstmp,upstmp,thetal,q2tmp,q2ha0_m
      REAL(rkind) :: sc1,sc2,sc3,dz1,dz2,dz3
      REAL(rkind) :: sum_sc,sum_dz,sum_scdz,sum_dz2,sp_limit

      REAL(rkind),allocatable :: theta(:),Rep(:)
      REAL(rkind),allocatable :: Scdp(:),tau_p(:),tau_f(:)
      REAL(rkind),allocatable :: theta_cr(:),Dstar(:)
      REAL(rkind),allocatable :: Dpm(:),Keci(:),Prob(:)

      REAL(rkind) :: ERF,thetacr

      IF(myrank==0) WRITE(mirror,*)'SED: Start the pick-up function'
      !- Start Statement --------------------------------------------------!
      allocate(theta(ntr_l),Rep(ntr_l),Scdp(ntr_l), &
                &tau_p(ntr_l),tau_f(ntr_l),Prob(ntr_l), &
                &Dstar(ntr_l),Dpm(ntr_l),Keci(ntr_l), &
                &theta_cr(ntr_l),stat=i)
      if(i/=0) call parallel_abort('Pick-up function: fail to allocate')
      
      DO i=1,nea
          
        cff0=0.d0; cff1=0.d0; cff2=0.d0; cff3=0.d0; cff4=0.d0; cff5=0.d0
        tmp=0.d0; downstmp=0.d0; upstmp=0.d0;FAI=0.d0;thetal=0.d0
        theta=0.d0;Rep=0.d0;Scdp=0.d0;tau_p=0.d0;tau_f=0.d0;q2tmp=0.d0
        theta_cr=0.d0;Dstar=0.d0;Dpm=0.d0;Keci=0.d0;Prob=0.d0;q2ha0_m=0.d0

        IF(idry_e(i)==1) CYCLE

!       cff0=(bustr(i)*bustr(i)+bvstr(i)*bvstr(i))**0.25d0                               !friction velocity  
        cff0=(abs(bustr(i))+abs(bvstr(i)))**0.5d0          
!************************************************************************** 
        do ised=1,ntr_l
          indx=isand(ised) !into 1:ntracers
          if(cff0.LE.Thero_ustar) cff0=Thero_ustar
          s=Srho(ised)/rho0
          cff5=((s-1)*grav*Sd50(ised))**0.5d0
          theta(ised)=min(10.d0,max(1.d-14,cff0*cff0/(grav*Sd50(ised)*(s-1))))            !Shields number
          thetal=4.d0/(3.d0*Clift)
          Rep(ised)=cff0*Sd50(ised)/nuf                                !Sediment reynolds number
          if (Rep(ised).LE.1.d-14) THEN
            WRITE(errmsg,*)'sediment reynolds number is zero',ised,ielg(i),Rep(ised)
            CALL parallel_abort(errmsg)
          end if
          Scdp(ised)=((32.d0/Rep(ised))**0.67d0+1.d0)**1.5d0                           !Cd
          Dstar(ised)=Sd50(ised)*(grav*(s-1.0d0)/(nuf*nuf))**(1.0d0/3.0d0)             !D*
          theta_cr(ised)=thetacr(Dstar(ised))
          Dpm(ised)=1.2d0*(1.d0-exp(-0.095d0*Dstar(ised)))                             !Dp
!********************************************************************* above is the 
          tau_p(ised)=s*Wsed(ised)/((s-1)*grav)                      
          tau_f(ised)=0.4d0*kf*refht*Sd50(ised)/muf/cff0
          cff1=tau_f(ised)/tau_p(ised)*(1+2.d0*s)                                       !beta
          cff2=(3.d0+cff1)/(1.d0+cff1+2.d0*s)                                           !Ct
          Keci(ised)=2.d0*cff2*cff2*Dpm(ised)/(3.d0*muf)
          if (Keci(ised).LE.1.d-10) THEN
            WRITE(errmsg,*)'Keci is zero',ised,ielg(i),Keci(ised)
            CALL parallel_abort(errmsg)
          end if
!*********************************************************************  
          tmp=Bstar/theta(ised)-1.d0/Eta0
          if(tmp.GE.0.d0) then
            Prob(ised)=max(0.d0,0.5d0-0.5d0*ERF(tmp)) 
          else
            Prob(ised)=min(1.d0,0.5d0+0.5d0*ERF(abs(tmp))) 
          endif

          select case(flag)
            case(0) !Zhong
              upstmp=Spmax*Prob(ised)*sqrt(Keci(ised))*cff0
              downstmp=sqrt(2.d0*pi)*(1+3.d0*Scdp(ised)*refht/(2.d0*s))
              FAI=exp(-1.d0*(1.d0/theta(ised)-1.d0/thetal)*refht/(s*Keci(ised)))
              sedcaty(i,ised)=Erate(ised)*upstmp*FAI/downstmp
            case(1) !Van
              cff3=cff5*(theta_cr(ised))**0.5d0
              upstmp=max(0.d0,(cff0*cff0/cff3/cff3-1.d0))
              sedcaty(i,ised)=Erate(ised)*0.00033*Dstar(ised)**0.3d0*upstmp**1.5d0*cff5
            case(2) !Cao
              cff3=cff5*(theta_cr(ised))**0.5d0
              upstmp=max(0.d0,(cff0*cff0/cff3/cff3-1.d0))
              downstmp=(0.02d0*Spmax/nuf/Tbp)*((s-1.d0)*grav)**0.5d0
              sedcaty(i,ised)=Erate(ised)*downstmp*Sd50(ised)**1.5d0*upstmp*theta(ised)*cff5
            case(3) !Cheng
              do k=kbe(i)+1,nvrt
                 do j=1,i34(i)
                    nd=elnode(j,i)
                    q2ha0_m=q2ha0_m+(q2(k,nd)+q2(k-1,nd))/2.d0  
                 end do !j
                q2ha0_m=q2ha0_m/i34(i)
                q2tmp=q2ha0_m*(ze(k,i)-ze(k-1,i))+q2tmp 
              end do
              q2tmp=q2tmp/(ze(nvrt,i)-ze(kbe(i),i))     !depth-average kinetic energy
              upstmp=(dave(i)/cff5)**8.d0
              downstmp=(q2tmp/(cff5*cff5))**(-0.8d0)
              sedcaty(i,ised)=Erate(ised)*2.81d0*1.d-13*upstmp*downstmp*Dstar(ised)**2.4d0*cff5
!Tsinghua group----------------------------
!            case(4) !Zhou
!              sedcaty(i,ised)=Erate(ised)*(tau_wc(i)/tau_ce(ised)-1.0d0)/Srho(ised)
             case(4) !cheng-new
                upstmp=dave(i)/cff5
                if (abs(dave(i)).LE.1.d-5) then
                    downstmp=0.d0
                else
                    downstmp=exp(-40.0/upstmp)
                endif
                sedcaty(i,ised)=Erate(ised)*1.d-4*upstmp*downstmp*Dstar(ised)**2.5d0*cff5
!Tsinghua group----------------------------
            case default
              call parallel_abort('SED: unknown im_pick_up')
          end select

!          if(ielg(i)==1) write(200,*),sedcaty(i,ised),FAI,upstmp,downstmp
!          if(ielg(i)==1) write(300,*),theta(ised),Rep(ised),Scdp(ised),theta_cr(ised),Dpm(ised)
!          if(ielg(i)==1) write(400,*),cff0,cff1,cff2,cff5,tmp
!          if(ielg(i)==1) write(500,*),tau_p(ised),tau_f(ised),Keci(ised),Prob(ised),Erate(ised)

          if (sedcaty(i,ised)/=sedcaty(i,ised)) THEN
            WRITE(errmsg,*)'pick function is NaN',ised,ielg(i),sedcaty(i,ised)
            CALL parallel_abort(errmsg)
          end if
        enddo  !ised=1,ntr_l
      ENDDO !i=1,nea

      deallocate(theta,Rep,Scdp,tau_p,tau_f,Dstar,Dpm,Keci,Prob,theta_cr)
      IF(myrank==0) WRITE(mirror,*)'SED: End the pick-up function'
      end subroutine sed_pickup  !End the calculation sediment carring-capacity. 

      function thetacr(D)
      use schism_glbl, only : rkind
      implicit none
   
      real(rkind) :: thetacr
      real(rkind), intent(in) ::D
      if(D.LE.4.0) thetacr=0.24d0/D
      if((D.GT.4.0).AND.(D.LE.10.0)) thetacr=0.14d0/(D**0.64d0)
      if((D.GT.10.0).AND.(D.LE.20.0)) thetacr=0.04d0/(D**0.10d0)
      if((D.GT.20.0).AND.(D.LE.150.0)) thetacr=0.013d0/(D**0.29d0)
      if(D.GT.150.0) thetacr=0.055d0
      end function thetacr 
