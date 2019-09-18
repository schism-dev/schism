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
! MORSELFE INITIALIZATION SUBROUTINES
!
! subroutine sed_alloc
! subroutine sed_init_phase[12]
!
!=====================================================================
!=====================================================================

      SUBROUTINE sed_alloc()
!--------------------------------------------------------------------!
! This subroutine allocates and pre-initialize sediment model arrays !
!                                                                    !
! Adapted from former subroutines initialize_scalars and             !
! initialize_ocean (former init_sed.F90), other allocation formerly  !
! done within schism_init.F90 were also moved within this subroutine. !
!                                                                    !
! Author: florian ganthy (fganthy@lnec.pt ; florian.ganthy@gmail.com)!
! Date: 2013/06/07                                                   !
!                                                                    !
! History:                                                           !
!            ***  Previous history (former subroutines) ***          !
!                                                                    !
!   2012/12 - F.Ganthy : form homogenisation of sediments  routines  !
!   2012/12 - F.Ganthy : modifications for Bottom Composition        !
!                        Generation (BCG) purpose (added output      !
!                        files - bed median grain size, bed fraction)!
!   2013/01 - F.Ganthy : Implementation of roughness predictor       !
!   2013/03 - F.Ganthy : Implementation of wave-induced bedload      !
!                        transport                                   !
!   2013/04 - F.Ganthy : Implementation of wave-current bottom stress!
!   2013/05 - F.Ganthy : Updates related to ripple predictor         !
!   2013/05 - F.Ganthy : Added node-centered  volume control area    !
!   2013/05 - F.Ganthy : Updates to the ripple predictor:            !
!                         - Changes on the total bedform             !
!                           roughness computation                    !
!                         - Add wave-ripple computation from         !
!                           Nielsen (1992)                           !
!                                                                    !
!            ***  Current history (former subroutines) ***           !
!                                                                    !
!                                                                    !
!--------------------------------------------------------------------!

      USE sed_mod

      USE schism_glbl, ONLY: rkind,dav,dave,nea,npa,nvrt,mnei_p, &
     &ntrs,irange_tr,ipgl,ielg,i34,elnode,np_global,  &
     & ifile_char,ifile_len,area,np,nne,indel,iself,nnp,indnd,nxq,     &
     & isbnd,rough_p,errmsg,xnd,ynd,xcj,ycj,xctr,yctr,elside,ics,xel,yel, &
     &in_dir,out_dir,len_in_dir,len_out_dir
      USE schism_msgp, ONLY: myrank,parallel_abort,exchange_p2d

      IMPLICIT NONE

!- Local variables --------------------------------------------------!

      REAL(rkind), PARAMETER :: IniVal = 0.0d0
      CHARACTER(len=48) :: inputfile
      INTEGER :: i,j,k,ie,jj,id,id1,id2,id3,nd,indx,m,nwild(3),nwild2(3)
      INTEGER :: ised,ic,itmp,istat

      REAL(rkind) :: aux1,aux2,xtmp,ytmp,tmp1,cff1,cff2,cff3,cff4,cff5,cff6
      REAL(rkind) :: bed_frac_sum,ar1,ar2
      real(rkind) :: signa

      !swild98 used for exchange (deallocate immediately afterwards)
      REAL(rkind),ALLOCATABLE :: swild98(:,:,:)

!- Start Statement --------------------------------------------------!

      IF(myrank==0) write(16,*)'Entering sed_alloc'

!--------------------------------------------------------------------!
!* ARRAYS ALLOCATION (a few have been done in read_sed_input)
!--------------------------------------------------------------------!

      !--------------------------------------------------------------!
      !* 1D arrays defined on elements
      !--------------------------------------------------------------!
      ALLOCATE(Zob(nea),stat=i)
        IF(i/=0) CALL parallel_abort('Sed: Zob allocation failure')
      ALLOCATE(dave(nea),stat=i)
        IF(i/=0) CALL parallel_abort('Sed: dave allocation failure')
      ALLOCATE(FX_r(nea), stat=i )
        IF(i/=0) CALL parallel_abort('Sed: FX_r allocation failure')
      ALLOCATE(FY_r(nea), stat=i )
        IF(i/=0) CALL parallel_abort('Sed: FY_r allocation failure')
      ALLOCATE(bustr(nea), stat=i )
        IF(i/=0) CALL parallel_abort('Sed: bustr allocation failure')
      ALLOCATE(bvstr(nea), stat=i )
        IF(i/=0) CALL parallel_abort('Sed: bvstr allocation failure')
      ALLOCATE(bed_thick(nea),stat=i)
        IF(i/=0) CALL parallel_abort('Sed: bed_thick allocation failure')
      ALLOCATE(hs(nea),stat=i)
        IF(i/=0) CALL parallel_abort('Sed: hs allocation failure')
      ALLOCATE(tp(nea),stat=i)
        IF(i/=0) CALL parallel_abort('Sed: tp allocation failure')
      ALLOCATE(wlpeak(nea),stat=i)
        IF(i/=0) CALL parallel_abort('Sed: wlpeak allocation failure')
      ALLOCATE(uorb(nea),stat=i)
        IF(i/=0) CALL parallel_abort('Sed: uorb allocation failure')
      ALLOCATE(uorbp(nea),stat=i)
        IF(i/=0) CALL parallel_abort('Sed: uorbp allocation failure')
      ALLOCATE(tau_c(nea),stat=i)
        IF(i/=0) CALL parallel_abort('Sed: tau_c allocation failure')
      ALLOCATE(tau_w(nea),stat=i)
        IF(i/=0) CALL parallel_abort('Sed: tau_w allocation failure')
      ALLOCATE(tau_wc(nea),stat=i)
        IF(i/=0) CALL parallel_abort('Sed: tau_wc allocation failure')

      !--------------------------------------------------------------!
      !* 2D arrays defined on elements
      !--------------------------------------------------------------!
      ALLOCATE(mcoefd(0:mnei_p,np),stat=i)
        IF (i/=0) CALL parallel_abort('Sed: mcoefd allocation failure')
      ALLOCATE(Hz(nvrt,nea),stat=i)
        IF(i/=0) CALL parallel_abort('Sed: Hz allocation failure')
      ALLOCATE(bottom(nea,MBOTP),stat=i)
        IF(i/=0) CALL parallel_abort('Sed: bottom allocation failure')
      ALLOCATE(sedcaty(nea,ntr_l),stat=i)                         !Tsinghua group:pick up sediment flux
      IF(i/=0) CALL parallel_abort('Sed: sedcaty allocation failure')

      !--------------------------------------------------------------!
      !* 3D arrays defined on elements
      !--------------------------------------------------------------!
      ALLOCATE(bed(Nbed,nea,MBEDP),stat=i)
      IF(i/=0) CALL parallel_abort('Sed: bed allocation failure')
      ALLOCATE(bed_frac(Nbed,nea,ntr_l),stat=i)
      IF(i/=0) CALL parallel_abort('Sed: bed_frac allocation failure')

      !--------------------------------------------------------------!
      !* 4D arrays defined on elements
      !--------------------------------------------------------------!
      ALLOCATE(bed_mass(Nbed,nea,2,ntr_l),stat=i)
      IF(i/=0) CALL parallel_abort('Sed: bed_mass allocation failure')

      !--------------------------------------------------------------!
      !* 1D arrays defined on nodes
      !--------------------------------------------------------------!
      ALLOCATE(vc_area(npa),stat=i)
      IF(i/=0) CALL parallel_abort('sed_alloc: vc_area allocation failure')
      ALLOCATE(bedthick_overall(npa),bed_d50n(npa),stat=i)
      IF(i/=0) CALL parallel_abort('sed_alloc: bed_d50n allocation failure')
      ALLOCATE(bed_taun(npa),stat=i)
      IF(i/=0) CALL parallel_abort('sed_alloc: bed_taun allocation failure')
      ALLOCATE(bed_rough(npa),stat=i)
      IF(i/=0) CALL parallel_abort('sed_alloc: bed_rough allocation failure')
      ALLOCATE(lbc_sed(npa),stat=i)
      IF(i/=0) CALL parallel_abort('sed_alloc: lbc_sed allocation failure')
      ALLOCATE(bc_sed(npa),stat=i)
      IF(i/=0) CALL parallel_abort('sed_alloc: bc_sed allocation failure')

      !--------------------------------------------------------------!
      !* 2D arrays defined on nodes
      !--------------------------------------------------------------!
      ALLOCATE(bedldu(npa,ntr_l),stat=i)
      IF(i/=0) CALL parallel_abort('sed_alloc: bedldu alloc failure')
      ALLOCATE(bedldv(npa,ntr_l),stat=i)
      IF(i/=0) CALL parallel_abort('sed_alloc: bedldv alloc failure')
      ALLOCATE(bed_fracn(npa,ntr_l),stat=i)
      IF(i/=0) CALL parallel_abort('sed_alloc: bed_fracn alloc failure')

      !--------------------------------------------------------------!
      !* Other 1D arrays
      !--------------------------------------------------------------!
      ALLOCATE ( isand(ntr_l), stat=i )
      IF(i/=0) CALL parallel_abort('Main: allocation failure')

!--------------------------------------------------------------------!
!* Pre-initialization
!--------------------------------------------------------------------!

      !--------------------------------------------------------------!
      !* 1D arrays defined on elements
      !--------------------------------------------------------------!
      Zob(:)       = IniVal
      dave(:)      = IniVal
      FX_r(:)      = IniVal
      FY_r(:)      = IniVal
      bustr(:)     = IniVal
      bvstr(:)     = IniVal
      bed_thick(:) = IniVal
      hs(:)        = IniVal
      tp(:)        = IniVal
      wlpeak(:)    = IniVal
      uorb(:)      = IniVal
      uorbp(:)     = IniVal
      tau_c(:)     = IniVal
      tau_w(:)     = IniVal
      tau_wc(:)    = IniVal

      !--------------------------------------------------------------!
      !* 2D arrays defined on elements
      !--------------------------------------------------------------!
      mcoefd(:,:) = IniVal
      Hz(:,:)     = IniVal
      bottom(:,:) = IniVal
      sedcaty(:,:)= IniVal                                           !Tsinghua group

      !--------------------------------------------------------------!
      !* 3D arrays defined on elements
      !--------------------------------------------------------------!
      bed(:,:,:)      = IniVal
      bed_frac(:,:,:) = IniVal

      !--------------------------------------------------------------!
      !* 3D arrays defined on elements
      !--------------------------------------------------------------!
      bed_mass(:,:,:,:) = IniVal

      !--------------------------------------------------------------!
      !* 1D arrays defined on nodes
      !--------------------------------------------------------------!
      vc_area(:)   = IniVal
      bed_d50n(:)  = IniVal
      bed_taun(:)  = IniVal
      bed_rough(:) = IniVal

      !Special initializations
      bc_sed(:)    = -9999
      lbc_sed(:)   = .FALSE.

      !--------------------------------------------------------------!
      !* 2D arrays defined on nodes
      !--------------------------------------------------------------!
      bedldu(:,:)    = IniVal
      bedldv(:,:)    = IniVal
      bed_fracn(:,:) = IniVal

      !--------------------------------------------------------------!
      !* Other 1D arrays
      !--------------------------------------------------------------!
      isand(:) = 1 ! init.

!--------------------------------------------------------------------!
! - Set tracer indices into 1:ntracers
! In SCHISM T and S are in different arrays from the tracers array...
!--------------------------------------------------------------------!
      ic = irange_tr(1,5)
      DO i = 1,ntr_l
        isand(i) = ic
        ic = ic+1
      END DO

!--------------------------------------------------------------------!
! - Computes matrix coefficients for the JCG solver
! Used for the computation of depth variation induced by bedload
!--------------------------------------------------------------------!
      mcoefd = 0
      aux1=22.0d0/108.0d0
      aux2=7.0d0/108.0d0
      DO i=1,np !residents
        DO j=1,nne(i)
          ie=indel(j,i)
          id=iself(j,i)
          if(i34(ie)==3) then
            mcoefd(0,i) = mcoefd(0,i)+area(ie)*aux1 !diagonal

            !Other 2 nodes
            do jj=1,2 !other 2 nodes
              nd=elnode(nxq(jj,id,i34(ie)),ie)
              indx=0
              do m=1,nnp(i)
                if(indnd(m,i)==nd) then
                  indx=m; exit
                endif
              enddo !m
              if(indx==0) call parallel_abort('SED_INIT: failed to find')
            
              mcoefd(indx,i)=mcoefd(indx,i)+area(ie)*aux2
            enddo !jj
!            IF(isbnd(1,i)==0.and.j==nne(i)) THEN !internal ball
!               mcoefd(1,i) = mcoefd(1,i)+area(ie)*aux2
!            ELSE
!               mcoefd(j+1,i) = mcoefd(j+1,i)+area(ie)*aux2
!            ENDIF
          else !quad
            id1=nxq(1,id,i34(ie))
            id2=nxq(2,id,i34(ie))
            id3=nxq(3,id,i34(ie))
            if(ics==1) then
              ar1=signa(xnd(i),xcj(elside(id3,ie)),xctr(ie),ynd(i),ycj(elside(id3,ie)),yctr(ie)) 
              ar2=signa(xnd(i),xctr(ie),xcj(elside(id2,ie)),ynd(i),yctr(ie),ycj(elside(id2,ie)))
            else !ll
              cff1=(xel(id,ie)+xel(id1,ie))/2 !xcj(elside(id3,ie))
              cff2=(yel(id,ie)+yel(id1,ie))/2 !ycj(elside(id3,ie))
              cff3=sum(xel(1:4,ie))/4 !xctr
              cff4=sum(yel(1:4,ie))/4 !yctr
              cff5=(xel(id,ie)+xel(id3,ie))/2 !xcj(elside(id2,ie))
              cff6=(yel(id,ie)+yel(id3,ie))/2 !ycj(elside(id2,ie))
              ar1=signa(xel(id,ie),cff1,cff3,yel(id,ie),cff2,cff4)
              ar2=signa(xel(id,ie),cff3,cff5,yel(id,ie),cff4,cff6)
            endif !ics
            if(ar1<=0.or.ar2<=0) call parallel_abort('SED_INIT:area<=0')
            mcoefd(0,i)=mcoefd(0,i)+(ar1+ar2)*7./12 !diagonal

            !Find indices
            do jj=1,3
              nd=elnode(nxq(jj,id,i34(ie)),ie)
              indx=0
              do m=1,nnp(i)
                if(indnd(m,i)==nd) then
                  indx=m; exit
                endif
              enddo !m
              if(indx==0) call parallel_abort('SED_INIT: faile to find2')
              nwild(jj)=indx
            enddo !jj

            mcoefd(nwild(1),i)=mcoefd(nwild(1),i)+ar1/4+ar2/12
            mcoefd(nwild(3),i)=mcoefd(nwild(3),i)+ar1/12+ar2/4
            mcoefd(nwild(2),i)=mcoefd(nwild(2),i)+ar1/12+ar2/12
          endif !i34
        ENDDO ! END loop nne
      ENDDO ! END loop np

!--------------------------------------------------------------------!
! - Control volume at each node used in filter
!   Split quads into 2 tri's
!--------------------------------------------------------------------!
      vc_area = 0.0d0
      nwild2(1:3)=(/1,3,4/) !prep. indices for 2nd tri of quad
      DO i=1,nea
        if(i34(i)==4) then !2 areas
          !nwild(1:3)=elnode(1:3,i)
          !ar1=signa(xnd(nwild(1)),xnd(nwild(2)),xnd(nwild(3)),ynd(nwild(1)),ynd(nwild(2)),ynd(nwild(3)))
          ar1=signa(xel(1,i),xel(2,i),xel(3,i),yel(1,i),yel(2,i),yel(3,i))
          !nwild(2)=elnode(3,i)
          !nwild(3)=elnode(4,i)
          !ar2=signa(xnd(nwild(1)),xnd(nwild(2)),xnd(nwild(3)),ynd(nwild(1)),ynd(nwild(2)),ynd(nwild(3)))
          ar2=signa(xel(1,i),xel(3,i),xel(4,i),yel(1,i),yel(3,i),yel(4,i))
          if(ar1<=0.or.ar2<=0) call parallel_abort('SED_INIT:area2<=0')
        endif

        DO j=1,3
          if(i34(i)==3) then
            nd=elnode(j,i)
            vc_area(nd)=vc_area(nd)+area(i)/3
          else !quad
            !1st tri
            nd=elnode(j,i)
            vc_area(nd)=vc_area(nd)+ar1/3

            !2nd tri
            nd=elnode(nwild2(j),i)
            vc_area(nd)=vc_area(nd)+ar2/3
          endif !i34(i)
        ENDDO !j
      ENDDO ! i
      CALL exchange_p2d(vc_area)

!     Read in total bed thickness at nodes
      OPEN(10,FILE=in_dir(1:len_in_dir)//'bedthick.ic',STATUS='OLD')
      READ(10,*); READ(10,*)
      DO i = 1,np_global
        READ(10,*) itmp,xtmp,ytmp,tmp1
        IF(tmp1<0) CALL parallel_abort('Sed_init: bed_thick<0!')
        IF(ipgl(i)%rank==myrank) bedthick_overall(ipgl(i)%id)=tmp1
      ENDDO !i=1,np_global
      CLOSE(10)

!     For cold start only
      !if(ihot==0) then
!========================================================
!--------------------------------------------------------------------!
! - Reading bed_frac files and convert bed_fraction from nodes to 
! elements.
! For instance, bed_fraction is applied to all bed layers, i.e. no vertical variation
!--------------------------------------------------------------------!
        ALLOCATE(swild98(Nbed,npa,ntr_l),stat=istat)
        ! * READING bed_frac_x.ic
        DO ised = 1,ntr_l
          WRITE(ifile_char,'(i03)')ised
          ifile_char=ADJUSTL(ifile_char); ifile_len=LEN_TRIM(ifile_char)
          inputfile='bed_frac_'//ifile_char(1:ifile_len)//'.ic'
          OPEN(10,FILE=in_dir(1:len_in_dir)//inputfile,STATUS='OLD')
          READ(10,*) !read in first line, no need to store it
          READ(10,*) !read in second line, no need to store it
          DO i = 1,np_global
            READ(10,*) itmp,xtmp,ytmp,tmp1
            IF(tmp1<0.or.tmp1>1) CALL parallel_abort('Sed: bed_frac wrong!')
            IF(ipgl(i)%rank==myrank) swild98(:,ipgl(i)%id,ised) = tmp1
          ENDDO !i=1,np_global
          CLOSE(10)
        ENDDO !ised=1,ntr_l

!'-------------------------------------------------------------------!
! - Mapping bed fraction (conversion from node to elements)
!--------------------------------------------------------------------!
        DO ised = 1,ntr_l
          DO k = 1,Nbed
            DO i = 1,nea
              bed_frac(k,i,ised) = sum(swild98(k,elnode(1:i34(i),i),ised))/i34(i)
            ENDDO ! END loop nea
          ENDDO ! END loop Nbed
        ENDDO ! END loop ntr_l
        DEALLOCATE(swild98,stat=istat)

!--------------------------------------------------------------------!
! - Check the sum of bed fractions 
!   Sum should be equal to 1 but not exactly possible due to rounding
!       - Sum must not exceed 1.01
!       - Sum must not be lower than 0.99
!--------------------------------------------------------------------!
        DO k=1,Nbed
          DO i=1,nea
            bed_frac_sum = 0.0d0
            DO ised=1,ntr_l
              bed_frac_sum = bed_frac_sum+bed_frac(k,i,ised)
            ENDDO !End loop ntr_l
            IF (bed_frac_sum>1.01d0.or.bed_frac_sum<0.99d0) THEN
              WRITE(errmsg,*)'SED: sum of bed_frac/=1 at elem ',ielg(i),bed_frac_sum
              CALL parallel_abort(errmsg)
            ENDIF
          ENDDO !i
        ENDDO !k

!--------------------------------------------------------------------!
! - Initialize the bed model (layer thickness, age and porosity)
!--------------------------------------------------------------------!
        DO i=1,Nbed
          DO j=1,nea
            bed(i,j,ithck) = sum(bedthick_overall(elnode(1:i34(j),j)))/i34(j)/real(Nbed) !>=0
            bed(i,j,iaged) = 0.0d0
            bed(i,j,iporo) = porosity
          ENDDO ! End loop Nbed
        ENDDO !End loop nea
!========================================================
!      endif !ireset

      IF(myrank==0) write(16,*)'Leaving sed_alloc'
!--------------------------------------------------------------------!
      END SUBROUTINE sed_alloc      
      
      
!=====================================================================
!=====================================================================

      SUBROUTINE sed_init
!--------------------------------------------------------------------!
! This subroutine initialize variables and arrays for the sediment   ! 
! model:                                                             !
!        - computes coefficients for the jcg solver                  !
!        - computes control volume area at nodes                     !
!        - read bed fraction file (bed_frac_x.ic)                    !
!        - mapping of bed fraction                                   !
!        - checking sum of bed fraction                              !
!        - initialize the bed model (layer thickness, age, porosity) !
!        - initialize the bed mass                                   !
!        - initialize exposed sediment layer properties              !
!        - initialize sediment roughness length                      !
!        - initialize total bed thickness                            !
!        - initialization of arrays defined at nodes                 !
!        - if debuging: write sediment initialization                !
!                                                                    !
! Adapted from former subroutines initialize_scalars,                !
! initialize_ocean (former init_sed.F90), and sed_init (former       !
! sed_init.F90) and other initializations formerly done within       !
!schism_init.F90 were also moved within this subroutine.              !
!                                                                    !
! The former subroutine sed_init.F90 was adapted from ROMS routine   !
! ana_sediment.h                                                     !
! Copyright (c) 2002-2007 The ROMS/TOMS Group                        !
!   Licensed under a MIT/X style license                             !
!   See License_ROMS.txt                                             !
!                                                                    !
! Author: florian ganthy (fganthy@lnec.pt ; florian.ganthy@gmail.com)!
! Date: 2013/06/07                                                   !
!                                                                    !
! History:                                                           !
!            ***  Previous history (former subroutines) ***          !
!                                                                    !
!   2012/12 - F.Ganthy : form homogenisation of sediments  routines  !
!   2012/12 - F.Ganthy : modifications for Bottom Composition        !
!                        Generation (BCG) purpose (added output      !
!                        files - bed median grain size, bed fraction)!
!   2013/01 - F.Ganthy : Implementation of roughness predictor       !
!   2013/03 - F.Ganthy : Implementation of wave-induced bedload      !
!                        transport                                   !
!   2013/04 - F.Ganthy : Implementation of wave-current bottom stress!
!   2013/05 - F.Ganthy : Updates related to ripple predictor         !
!   2013/05 - F.Ganthy : Added node-centered  volume control area    !
!   2013/05 - F.Ganthy : Updates to the ripple predictor:            !
!                         - Changes on the total bedform             !
!                           roughness computation                    !
!                         - Add wave-ripple computation from         !
!                           Nielsen (1992)                           !
!                                                                    !
!            ***  Current history (former subroutines) ***           !
!   2013/06 - F.Ganthy : Also removed old initializations currently  !
!                        which can be found in former svn revision   !
!
!                                                                    !
!                                                                    !
!--------------------------------------------------------------------!      

      USE sed_mod

      USE schism_glbl, ONLY: nea,npa,mnei_p,ntrs,irange_tr,ipgl,ielg,i34,elnode,np_global,  &
     & ifile_char,ifile_len,area,np,nne,indel,iself,nnp,indnd,nxq,     &
     & isbnd,rough_p,errmsg,xnd,ynd,xcj,ycj,xctr,yctr,elside,ics,xel,yel, &
     &in_dir,out_dir,len_in_dir,len_out_dir
      USE schism_msgp, ONLY: myrank,parallel_abort,exchange_p2d

      IMPLICIT NONE

      real(rkind) :: signa

!- Local variables --------------------------------------------------!

      INTEGER :: i,j,k,ie,jj,id,id1,id2,id3,nd,indx,m,nwild(3),nwild2(3)
      INTEGER :: ised,ic,itmp,istat

      REAL(rkind) :: aux1,aux2,xtmp,ytmp,tmp1,cff1,cff2,cff3,cff4,cff5,cff6
      REAL(rkind) :: bed_frac_sum,ar1,ar2

      REAL(rkind),DIMENSION(npa) :: bdfc

      CHARACTER(len=48) :: inputfile

!- Start Statement --------------------------------------------------!

      IF(myrank==0) WRITE(16,*)'Entering sed_init'

!--------------------------------------------------------------------!
! - Initialize the bed mass ([kg/m2]=[m]*[kg.m-3]*[-])
!--------------------------------------------------------------------!
      DO k=1,Nbed
        DO j=1,nea
          DO ised=1,ntr_l
             bed_mass(k,j,:,ised) = bed(k,j,ithck)*                  &
             &                      Srho(ised)*                      &
             &                      (1.0d0-bed(k,j,iporo))*          &
             &                      bed_frac(k,j,ised)
          ENDDO ! End loop ntr_l
        ENDDO ! End loop nea
      END DO ! End loop Nbed

!--------------------------------------------------------------------!
! - Initialize exposed sediment layer properties:
!      - total D50 (m)
!      - total density (kg.m-3)
!      - total settling velocity (m.s-1)
!      - total critical shear stress (m2.s-2)
!--------------------------------------------------------------------!
      DO j=1,nea
        cff1 = 1.0d0
        cff2 = 1.0d0
        cff3 = 1.0d0
        cff4 = 1.0d0
        cff5 = 0.0d0
        ! Weighted geometric mean
        DO ised=1,ntr_l
          cff1 = cff1*Sd50(ised)**bed_frac(1,j,ised)
          cff2 = cff2*Srho(ised)**bed_frac(1,j,ised)
          cff3 = cff3*Wsed(ised)**bed_frac(1,j,ised)
          cff4 = cff4*tau_ce(ised)**bed_frac(1,j,ised)
          cff5 = cff5+bed_frac(1,j,ised)
        ENDDO ! End loop ntr_l

        if(cff5<0) then
          call parallel_abort('SED_INIT: cff5<0 (1)')
        else if(cff5==0) then !care-takers for all-eroded case
          WRITE(12,*)'SED_INIT: all eroded at elem. ',ielg(j)
          bottom(j,isd50) = Sd50(1)
          bottom(j,idens) = Srho(1)
          bottom(j,iwsed) = Wsed(1)
          bottom(j,itauc) = tau_ce(1)
        else
          bottom(j,isd50) = cff1**(1.0d0/cff5)
          bottom(j,idens) = cff2**(1.0d0/cff5)
          bottom(j,iwsed) = cff3**(1.0d0/cff5)
          bottom(j,itauc) = cff4**(1.0d0/cff5)
        endif !cff5
      ENDDO ! j

      !Init. active layer thickness
      bottom(:,iactv)=bed(1,:,ithck)
!--------------------------------------------------------------------!
! - Initialize sediment main roughness length and bedform predictor
! related roughness length
! Here, the initialization is the same whether or not the bedfrom
! predictor is used
!--------------------------------------------------------------------!
      ! Current ripple roughness length (Soulsby, 1997)
      bottom(:,izcr)  = 0.0d0
      ! Sand waves roughness length (Van Rijn, 1984)
      bottom(:,izsw)  = 0.0d0
      ! Wave ripples roughness length (Grant and Madsen, 1982 OR Nielsen, 1992)
      bottom(:,izwr)  = 0.0d0
      ! Bed load bottom roughness length (Grant and Madsen, 1982 OR Nielsen, 1992)
      bottom(:,izbld) = 0.0d0

      DO i=1,nea
        ! Nikurasde roughness length
        bottom(i,izNik) = bottom(i,isd50)/12.0d0
        ! Default roughness
        bottom(i,izdef) = sum(rough_p(elnode(1:i34(i),i)))/i34(i) 
        ! Apparent initial roughness
        bottom(i,izapp) = bottom(i,izdef)
        ! Roughness length effectively used (even if bedform predictor is not used)
        Zob(i) = bottom(i,izdef)
        IF(Zob(i).LE.0.0d0)THEN
          WRITE(errmsg,*) 'SED_INIT: Zob <=0 at elem, rank:', Zob(i),i,myrank
          CALL parallel_abort(errmsg)
        ENDIF
      ENDDO ! i

!--------------------------------------------------------------------!
! - Initialize total bed thickness by summing over all layers
!   This array is not actively updated after this
!--------------------------------------------------------------------!
      bed_thick(:)=0.0d0
      DO i=1,Nbed
        DO j=1,nea
          bed_thick(j) = bed_thick(j)+bed(i,j,ithck)
        ENDDO !j
      ENDDO !i 

!--------------------------------------------------------------------!
! - Initialization of variables defined at nodes
!--------------------------------------------------------------------!
      bed_d50n(:)    = 0.0d0
      bed_taun(:)    = 0.0d0
      bed_fracn(:,:) = 0.0d0
      bdfc(:)        = 0.0d0
      DO i=1,nea
        DO j=1,i34(i)
          DO ised=1,ntr_l
            bed_fracn(elnode(j,i),ised) = bed_fracn(elnode(j,i),ised)+       &
            &                         bed_frac(1,i,ised)
          ENDDO !ised
          bed_d50n(elnode(j,i))  = bed_d50n(elnode(j,i))+bottom(i,isd50)
          bdfc(elnode(j,i))      = bdfc(elnode(j,i))+1.0d0
        ENDDO !j
      ENDDO !i
      DO i=1,np
        if(bdfc(i)==0) call parallel_abort('SED_INIT: bdfc(i)==0')
        DO ised=1,ntr_l
          bed_fracn(i,ised) = bed_fracn(i,ised)/bdfc(i)
        ENDDO !ised
        bed_d50n(i)    = bed_d50n(i)/bdfc(i)
        bed_rough(i)   = rough_p(i) !roughness length
      ENDDO !i
      DO ised=1,ntr_l
        CALL exchange_p2d(bed_fracn(:,ised))
      ENDDO !END loop ntr_l
      CALL exchange_p2d(bed_d50n)
      CALL exchange_p2d(bed_rough)

!--------------------------------------------------------------------!
! - Sediment debug outputs
!--------------------------------------------------------------------!
!      IF(sed_debug==1) THEN
      DO i=1,Nbed
        DO j=1,nea
          DO k=1,ntr_l
            WRITE(12,*)'Nbed:',i,' nea:',j,' ntr_l:',k,       &
              &                    ' bed_frac.ic:',real(bed_frac(i,j,k))
          ENDDO !k
        ENDDO !j
      ENDDO !i

      DO i=1,Nbed
        DO j=1,nea
          WRITE(12,*)'Nbed:',i,' nea:',j,' bed_thick:',real(bed_thick(j)),  &
            &     ' bed(ithck):',real(bed(i,j,ithck)),' bed(iaged):',    &
            &     real(bed(i,j,iaged)),' bed(iporo):',real(bed(i,j,iporo))
        ENDDO !j
      ENDDO !i

      DO i=1,nea
        WRITE(12,*)'nea:',i,' bottom(itauc):',real(bottom(i,itauc)),        &
          &     ' bottom(isd50):',real(bottom(i,isd50)),                 &
          &     ' bottom(iwsed):',real(bottom(i,iwsed)),                 &
          &     ' bottom(idens):',real(bottom(i,idens))
      ENDDO !i
!      ENDIF !End sed_debug

!--------------------------------------------------------------------!
      IF(myrank==0) write(16,*)'Leaving sed_init'
!--------------------------------------------------------------------!
      END SUBROUTINE sed_init

