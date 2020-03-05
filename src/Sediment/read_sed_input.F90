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

      SUBROUTINE read_sed_input
!--------------------------------------------------------------------!
! This subroutine reads sediment model inputs                        !
!                                                                    !
! Author: Ligia Pinto                                                !
! Date: xx/09/2007                                                   !
!                                                                    !
! History: 2012/12 - F.Ganthy : form homogenisation of sediments     !
!          routines                                                  !
!          2013/01 - F.Ganthy : Introduction of switchs for roughness!
!          predictor                                                 !
!          2013/01 - F.Ganthy : Implementation of avalanching        !
!          2013/05 - F.Ganthy : Updates related to ripple predictor  !
!          2013/05 - F.Ganthy : Updates to the ripple predictor:     !
!                               - Changes on the total bedform       !
!                                 roughness computation              !
!                               - Add wave-ripple computation from   !
!                                 Nielsen (1992)                     !
!          2013/05 - F.Ganthy : Add different sediment behavior:     !
!                                - MUD-like or SAND-like             !
!          2013/05 - F.Ganthy : Re-introduction of the number of bed !
!                               sediment layers within sediment.in   !
!          2020/02 - B.Mengual                                       !
!                  * Implementation of a method to compute and update!
!                    a porosity representative of non-cohesive       !
!                    matrix (poro_option,poro_cst,Awooster,Bwooster) !
!                  * Implementation of several methods to compute    !
!                    the current-induced bottom shear stress         !
!                    (tau_option,tau_max,zstress)                    !
!                  * Wave asymmetry computations: restructuration    !
!                    and introduction of filtering possibilities     !
!                    (iasym,w_asym_max,elfrink_filter,ech_uorb)      !
!                  * Implementation of bedload flux caused by        !
!                    wave acceleration skewness                      !
!                    (bedload_acc,kacc_hoe,                          !
!                    kacc_dub,thresh_acc_opt,acrit)                  !
!                  * Add the possibility to filter and limit         !
!                    bedload fluxes (bedload_filter,bedload_limiter, !
!                    bedload_acc_filter)
!                  * Add the possibility to use different methods    !
!                    to solve the Exner equation in bedchange_bedload!
!                    (IMETH_BED_EVOL)                                !
!                  * Introduction of a maximum active layer thickness!
!                    (actv_max)                                      !
!                  * Morphodynamics: common morph_fac over all       !
!                    sediment classes (morph_fac)                    !
!                  * Initialisation of the vertical discretisation of!
!                    sediment layers (sedlay_ini_opt,toplay_inithick)!
!                  * Bed slope effects for bedload fluxes            !
!                    (slope_formulation=4 with alpha_bs,alpha_bn)    !
!--------------------------------------------------------------------!

      USE schism_glbl, ONLY: rkind,ntrs,iwsett,wsett,irange_tr,ddensed, &
   &grav,rho0,errmsg,in_dir,out_dir,len_in_dir,len_out_dir
      USE schism_msgp, ONLY: myrank,parallel_abort
      USE sed_mod      

      IMPLICIT NONE
      SAVE  

!- Local variables --------------------------------------------------!

      CHARACTER(len=100) :: var1, var2
      INTEGER :: i, j, ierror

      INTEGER :: line, len_str, loc, loc2, lstr_tmp
      CHARACTER(len=90) :: line_str,str_tmp,str_tmp2

!- Start Statement --------------------------------------------------!
      !defined # of classes (in sed_mod) for all routines
      ntr_l=ntrs(5) 

      IF(myrank==0) write(16,*)'Entering read_sed_input routine'

      ALLOCATE(Sd50(ntr_l),Srho(ntr_l),Wsed(ntr_l),Erate(ntr_l),tau_ce(ntr_l), &
     &poros(ntr_l),iSedtype(ntr_l),stat=i)
      IF(i/=0) CALL parallel_abort('main: Sd50 allocation failure')
       
!--------------------------------------------------------------------
! - Scan sediment.in.
!--------------------------------------------------------------------

      !Take some old parameters out of sediment.in
      g=grav
      rhom=rho0
      vonKar=0.4

      !Init. new parameters to make sure they are read in
      ised_bc_bot=-100
      !depo_scale=-100
      ised_dump=-100
      ierosion=-100

      OPEN(5,FILE=in_dir(1:len_in_dir)//'sediment.in',STATUS='old')
      IF(myrank==0) WRITE(16,*)'reading sediment.in'
      IF(myrank==0) WRITE(16,'(A,A,I2,A)')'##### Number of Tracers',  &
      &             ' / Sediment Classes required in sediment.in: ', &
      &             ntr_l,' #####'

      REWIND(5) !go to beginning of file
      line=0

      ! - PS: start reading
      DO
        ! - PS: read line an break on '!'
        line=line+1
        READ(5,'(a)',END=99)line_str
        line_str=ADJUSTL(line_str)
        len_str=LEN_TRIM(line_str)
        IF(len_str==0.OR.line_str(1:1)=='!') CYCLE

        ! - PS: locate '==' and '!'
        loc=INDEX(line_str,'==')
        loc2=INDEX(line_str,'!')
        IF(loc2/=0.AND.loc2-1<loc+1)THEN
          CALL parallel_abort('READ_PARAM: ! before =')
        ENDIF

        ! - PS: get name of the variable
        str_tmp=''
        str_tmp(1:loc-1)=line_str(1:loc-1)
        str_tmp=TRIM(str_tmp)
        lstr_tmp=LEN_TRIM(str_tmp)

        ! PS: get the value
        !if(loc2/=0) then
        !  str_tmp2=line_str(loc+2:loc2-1)
        !else
        !  str_tmp2=line_str(loc+2:len_str)
        !endif

        !str_tmp2=adjustl(str_tmp2)
        !str_tmp2=trim(str_tmp2)

        ! - PS: switch between variables and set values
        SELECT CASE(str_tmp)

          CASE('NEWLAYER_THICK')
            READ(line_str,*) var1, var2, newlayer_thick

          CASE('BEDLOAD_COEFF')
            READ(line_str,*) var1, var2, bedload_coeff

          CASE('SAND_SD50')
            READ(line_str,*) var1, var2, (Sd50(i), i=1,ntr_l)

          CASE('SAND_SRHO')
            READ(line_str,*) var1, var2, (Srho(i), i=1,ntr_l)

          CASE('SAND_WSED')
            READ(line_str,*) var1, var2, (Wsed(i), i=1,ntr_l)

          CASE('SAND_ERATE')
            READ(line_str,*) var1, var2, (Erate(i), i=1,ntr_l)

          CASE('SAND_TAU_CE')
            READ(line_str,*) var1, var2, (tau_ce(i), i=1,ntr_l)

          CASE('SED_TYPE')
            READ(line_str,*) var1, var2, (iSedtype(i), i=1,ntr_l)

          CASE('sed_debug')
            READ(line_str,*) var1, var2, sed_debug

!          CASE('g')
!            READ(line_str,*) var1, var2, g

!          CASE('vonKar')
!            READ(line_str,*) var1, var2, vonKar

          CASE('Cdb_min')
            READ(line_str,*) var1, var2, Cdb_min

          CASE('Cdb_max')
            READ(line_str,*) var1, var2, Cdb_max

          CASE('actv_max')
            READ(line_str,*) var1, var2, actv_max

!          CASE('rhom')
!            READ(line_str,*) var1, var2, rhom

          CASE('poro_option')
            READ(line_str,*) var1, var2, poro_option

          CASE('poro_cst')
            READ(line_str,*) var1, var2, porosity

          CASE('Awooster')
            READ(line_str,*) var1, var2, Awooster

          CASE('Bwooster')
            READ(line_str,*) var1, var2, Bwooster

          CASE('sedlay_ini_opt')
            READ(line_str,*) var1, var2, sedlay_ini_opt

          CASE('toplay_inithick')
            READ(line_str,*) var1, var2, toplay_inithick

          CASE('bdldiffu')
            READ(line_str,*) var1, var2, bdldiffu

          CASE('Nbed')
            READ(line_str,*) var1, var2, Nbed

          CASE('bedload')
            READ(line_str,*) var1, var2, bedload

          CASE('bedload_filter')
            READ(line_str,*) var1, var2, bedload_filter

          CASE('bedload_limiter')
            READ(line_str,*) var1, var2, bedload_limiter

          CASE('suspended_load')
            READ(line_str,*) var1, var2, suspended_load

          CASE('slope_formulation')
            READ(line_str,*) var1, var2, slope_formulation

          CASE('alpha_bs')
            READ(line_str,*) var1, var2, alpha_bs

          CASE('alpha_bn')
            READ(line_str,*) var1, var2, alpha_bn

          CASE('ised_bc_bot')
            READ(line_str,*) var1, var2, ised_bc_bot

          CASE('sed_morph')
            READ(line_str,*) var1, var2, sed_morph

          !CASE('depo_scale')
          !  READ(line_str,*) var1, var2, depo_scale

          CASE('sed_morph_time')
            READ(line_str,*) var1, var2, sed_morph_time !in days

          CASE('sed_morph_fac')
            READ(line_str,*) var1, var2, morph_fac ! morphological factor [-]

          CASE('drag_formulation')
            READ(line_str,*) var1, var2, drag_formulation

          CASE('ddensed')
            READ(line_str,*) var1, var2, ddensed

          CASE('comp_ws')
            READ(line_str,*) var1, var2, comp_ws

          CASE('comp_tauce')
            READ(line_str,*) var1, var2, comp_tauce

          CASE('bedforms_rough')
            READ(line_str,*) var1, var2, bedforms_rough

          CASE('iwave_ripple')
            READ(line_str,*) var1, var2, iwave_ripple

          CASE('irough_bdld')
            READ(line_str,*) var1, var2, irough_bdld

          CASE('slope_avalanching')
            READ(line_str,*) var1, var2, slope_avalanching

          CASE('dry_slope_cr')
            READ(line_str,*) var1, var2, dry_slope_cr
 
          CASE('wet_slope_cr')
            READ(line_str,*) var1, var2, wet_slope_cr

          CASE('bedmass_filter')
            READ(line_str,*) var1, var2, bedmass_filter

          CASE('bedmass_threshold')
            READ(line_str,*) var1, var2, bedmass_threshold  

!          CASE('relath')
!            READ(line_str,*) var1, var2, relath

          CASE('ised_dump')
            READ(line_str,*) var1, var2, ised_dump

          CASE('ierosion')
            READ(line_str,*) var1, var2, ierosion

          CASE('alphd') !1120:+alphd,refht,Tbp,im_pick_up
            READ(line_str,*) var1, var2, alphd

          CASE('refht')
            READ(line_str,*) var1, var2, refht

          CASE('Tbp')
            READ(line_str,*) var1, var2, Tbp

          CASE('im_pick_up')
            READ(line_str,*) var1, var2, im_pick_up

          CASE('IMETH_BED_EVOL')
            READ(line_str,*) var1, var2, imeth_bed_evol

          CASE('iasym')
            READ(line_str,*) var1, var2, iasym

          CASE('w_asym_max')
            READ(line_str,*) var1, var2, w_asym_max

          CASE('elfrink_filter')
            READ(line_str,*) var1, var2, elfrink_filter

          CASE('ech_uorb')
            READ(line_str,*) var1, var2, ech_uorb

          CASE('bedload_acc')
            READ(line_str,*) var1, var2, bedload_acc

          CASE('bedload_acc_filter')
            READ(line_str,*) var1, var2, bedload_acc_filter

          CASE('kacc_hoe')
            READ(line_str,*) var1, var2, kacc_hoe

          CASE('kacc_dub')
            READ(line_str,*) var1, var2, kacc_dub

          CASE('thresh_acc_opt')
            READ(line_str,*) var1, var2, thresh_acc_opt

          CASE('acrit')
            READ(line_str,*) var1, var2, acrit

          CASE('tau_option')
            READ(line_str,*) var1, var2, tau_option

          CASE('tau_max')
            READ(line_str,*) var1, var2, tau_max

          CASE('zstress')
            READ(line_str,*) var1, var2, zstress

        END SELECT

      ENDDO !scan sediment.in
99    CLOSE(5)


!---------------------
!      g = g*1.0d0
!      vonKar = vonKar*1.0d0
!      Cdb_min = Cdb_min*1.0d0
!      Cdb_max = Cdb_max*1.0d0
!      rhom = rhom*1.0d0

      ! Conversion of threshold from mm to m
      bedmass_threshold = bedmass_threshold/1000.0d0
!---------------------

!     Catch unread parameters (for new additions only)
      if(ised_bc_bot<0.or.ised_dump<0.or.ierosion<0) then
        write(errmsg,*)'READ_SED, some pars. not read in:',ised_bc_bot,ised_dump,ierosion
        call parallel_abort(errmsg)
      endif 

! BM: check
      IF ((tau_option .EQ. 0) .OR. (tau_option .GT. 3)) THEN
        WRITE(errmsg,*)'READ_SED: Invalid tau_option (from 1 to 3)'
        CALL parallel_abort(errmsg)
      END IF



!---------------------------------------------------------------------
! - Scale input parameters
!
! Particel settling velocity (Wsed) input is in [mm/s]
! Median sediment grain diameter (Sd50) input is in [mm]
! Critical shear for erosion and deposition (Tau_ce) input is in [N/m2]
!---------------------------------------------------------------------

      DO i=1,ntr_l !ntracers
        Sd50(i) = Sd50(i)*0.001d0 !to [m]
        IF(comp_ws.EQ.1.AND.iSedtype(i).EQ.1) THEN
          ! Sed type is SAND
          CALL sed_settleveloc(i)
        ELSE
          ! Sed type is other (MUD or GRAVEL)
          Wsed(i) = Wsed(i)*0.001d0 !to m/s
        ENDIF
        IF(comp_tauce.EQ.1.AND.iSedtype(i).EQ.1)THEN
          ! Sed type is SAND
          CALL sed_taucrit(i)
        ENDIF

        ! tau_ce already read in for other Sed type (MUD or GRAVEL); scale to
        ! get m^2/s/s
        tau_ce(i) = tau_ce(i)/rhom
      ENDDO !i

      !Update global array for settling vel.
      do i=1,ntr_l
        wsett(i-1+irange_tr(1,5),:,:)=Wsed(i)
        iwsett(i-1+irange_tr(1,5))=1
      enddo !i

!---------------------------------------------------------------------
! - Screen output of sediment-model parameters
!---------------------------------------------------------------------

       IF(myrank==0) THEN
         WRITE(16,*)
         WRITE(16,*) 'Sediment Parameters:'
         WRITE(16,*) '--------------------'
         WRITE(16,*) 'sed_debug: ', sed_debug
         WRITE(16,*) 'Nbed: ',Nbed
         IF (poro_option .EQ. 1) THEN
           WRITE(16,*) 'Constant porosity: ',porosity
         ELSE
           WRITE(16,*) 'Variable porosity poro_option > 1'
         END IF
         IF (sed_morph .EQ. 1) THEN
           WRITE(16,*) 'Accounting for morphodynamics'
           WRITE(16,*) ' > sed_morph_time (ramp in days)=',sed_morph_time
           WRITE(16,*) ' > morphological factor: ',morph_fac
         END IF
!         WRITE(16,*) 'bedthick_overall(1): ',bedthick_overall(1)
         WRITE(16,*) 'New layer thickness',newlayer_thick
         WRITE(16,*) '----------------------------------------------'
         WRITE(16,*) 'Sed type  ;  Sd50 [m]  ;  Srho [kg/m3]  ;',     &
         &          '  Wsed [m/s]  ;  Erate ;  tau_ce', &
         &          ' [Pa]  ;  ierosion '
         DO i = 1,ntr_l !ntracers
           WRITE(16,*) iSedtype(i),Sd50(i),Srho(i),Wsed(i),Erate(i),   &
           &           tau_ce(i)*rhom,ierosion
         ENDDO
         WRITE(16,*)
       ENDIF


       IF(myrank==0) WRITE(16,*)'Leaving read_sed_input routine'
       RETURN
!--------------------------------------------------------------------!
       END SUBROUTINE read_sed_input
