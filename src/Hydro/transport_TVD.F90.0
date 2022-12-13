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

!===============================================================================
!===============================================================================
! SCHISM transport models
!
!  subroutine do_transport_tvd
!
!===============================================================================
!===============================================================================
!

!     Do upwind and TVD transport
      subroutine do_transport_tvd(it,ltvd,ntr,difnum_max_l) !,nvrt1,npa1,dfh1)

!#ifdef USE_MPIMODULE
!      use mpi
!#endif
      use schism_glbl
      use schism_msgp
      use misc_modules

!#ifdef USE_TIMOR
!      USE flmud_pool, only: wsink !wsink([],nvrt,npa)>=0 (positive down)
!#endif /*USE_TIMOR*/
      implicit none
!#ifndef USE_MPIMODULE
      include 'mpif.h'
!#endif

      integer, intent(in) :: it !time stepping #; info only
      logical, intent(in) :: ltvd !true if TVD is used (must be for all tracers)
!      character(len=2), intent(in) :: flimiter
      integer, intent(in) :: ntr !# of tracers (for dimensioning)
!      integer, intent(in) :: nvrt1 !,npa1 !for dimensioning
!      real(rkind), intent(in) :: dfh1(nvrt1,npa1) 
      real(rkind), intent(out) :: difnum_max_l !max. horizontal diffusion number reached by this process (check stability)

!     Functions used
      real(rkind) :: flux_lim


!     Working temporary arrays in this routine
      real(rkind) :: iupwind_e(nea) !to mark upwind prisms when TVD is used
      real(rkind), allocatable :: trel_tmp(:,:,:) !tracer @ elements and half levels
      real(rkind), allocatable :: flux_adv_hface(:,:) ! original horizontal flux (the local x-driection) 
      real(rkind), allocatable :: flux_mod_hface(:,:,:) !limited advective fluxes on horizontal faces
      real(rkind), allocatable :: flux_mod_vface(:,:,:) !limited advective fluxes on vertical faces
      real(rkind), allocatable :: up_rat_hface(:,:,:) !upwind ratios for horizontal faces
      real(rkind), allocatable :: up_rat_vface(:,:,:) !upwind ratios for vertical faces
      real(rkind) :: buf(2,1),buf2(2,1)

!#ifdef DEBUG
!      real(rkind) :: dtbe(ne)
!#endif

      real(rkind) :: psumtr(ntr),delta_tr(ntr),adv_tr(ntr),dtb_min3(ne), &
     &alow(nvrt),bdia(nvrt),cupp(nvrt),rrhs(ntr,nvrt),soln(ntr,nvrt),gam(nvrt), &
     &swild(max(3,nvrt)),swild4(3,2),trel_tmp_outside(ntr),surf_vol(ntr)
      integer :: nwild(2)

      integer :: istat,i,j,k,l,khh2,ie,n1,n2,n3,isd,isd0,isd1,isd2,isd3,j0, &
                 &nd,it_sub,ntot_v,ntot_vgb,ntot_h,ntot_hgb,kup,kdo,jsj,kb, &
                 &kb1,iup,ido,ie01,lev01,in_st,jj,ll,lll,ndim,kin,iel,ibnd, &
                 &ndo,ind1,ind2,nd1,nd2,ibio
      real(rkind) :: vnor1,vnor2,xcon,ycon,zcon,dot1,sum1,tmp,cwtmp,toth, &
                     &time_r,psum,rat,dtbl,dtb,vj,bigv,av_df,av_dz,hdif_tmp, &
                     &av_h,difnum,cwtmp2,dtb_by_bigv,vol

      real(rkind) :: ref_flux
      logical     :: same_sign, is_land
!      logical, save :: first_call
      

#ifdef INCLUDE_TIMING
      cwtmp2=mpi_wtime()
#endif

      allocate(trel_tmp(ntr,nvrt,nea), &
              &flux_adv_hface(nvrt,nsa), &
              &flux_mod_hface(ntr,nvrt,ns),flux_mod_vface(ntr,nvrt,ne), &
              &up_rat_hface(ntr,nvrt,nsa),up_rat_vface(ntr,nvrt,nea),stat=istat)
      if(istat/=0) call parallel_abort('Transport: fail to allocate')

!     For TVD, prepare some arrays for 2-tier ghosts
#ifdef INCLUDE_TIMING
      cwtmp=mpi_wtime()
      timer_ns(1)=timer_ns(1)+cwtmp-cwtmp2
#endif
      if(ltvd) then
        idry_e_2t(1:ne)=idry_e(1:ne)
        call exchange_e2di_2t(idry_e_2t) !now has values up to nea2
        call exchange_e3d_2t_tr(tr_el)
      endif !ltvd
#ifdef INCLUDE_TIMING
      wtimer(9,2)=wtimer(9,2)+mpi_wtime()-cwtmp
#endif

!$OMP parallel default(shared) private(j,is_land,k,vnor1,vnor2,i,nd,toth,kup,kdo,psum,psumtr, &
!$OMP jsj,ie,tmp,iup,ido,ind1,delta_tr,l,rat,jj,ref_flux,same_sign,vj,ndim,kin,alow,bdia,cupp, &
!$OMP bigv,dtb_by_bigv,av_df,av_dz,adv_tr,iel,trel_tmp_outside,ibnd,nwild,ll,ndo,lll,ind2,rrhs, &
!$OMP nd1,nd2,hdif_tmp,av_h,difnum,soln,gam)

!'    Modify here 3D velocity for transport (for whatever reason) - make sure volume conservation is not violated
!     For rewetted elements, tr_el takes the value from last wet step

!     Compute (pre-limiting) fluxes at all faces 
!$OMP workshare
      flux_adv_hface=-1.d34 !flags
!$OMP end workshare
!      flux_adv_vface=-1.d34 !flags

!     Horizontal fluxes
!$OMP do
      do j=1,ns !resident side
        if(idry_s(j)==1) cycle
        is_land=(isdel(2,j)==0.and.isbs(j)<=0)

        do k=kbs(j)+1,nvrt
          if(is_land) then !land
            flux_adv_hface(k,j)=0.d0
          else            
!            if(ics==1) then
            vnor1=su2(k,j)*snx(j)+sv2(k,j)*sny(j)
            vnor2=su2(k-1,j)*snx(j)+sv2(k-1,j)*sny(j)
!            else !lat/lon
!              vnor1=su2(k,j)
!              vnor2=su2(k-1,j)
!            endif !ics
            flux_adv_hface(k,j)=(zs(k,j)-zs(k-1,j))*distj(j)*(vnor1+vnor2)/2 !normal * area = flux (in local x-direction)

!         Debug
!         if(it==46.and.i==58422) write(99,*)j,k,vnor1,vnor2,flux_adv_hface(k,jsj)
          endif !is_land
        enddo !k=kbs(i)+1,nvrt
      enddo !j=1,ns
!$OMP end do

!     Compute vertical fluxes - already done in schism_step

!$OMP master
!     Exchange flux_adv
#ifdef INCLUDE_TIMING
      cwtmp=mpi_wtime()
      timer_ns(1)=timer_ns(1)+cwtmp-cwtmp2
#endif
      call exchange_s3dw(flux_adv_hface)
#ifdef INCLUDE_TIMING
      wtimer(9,2)=wtimer(9,2)+mpi_wtime()-cwtmp
#endif
 
#ifdef INCLUDE_TIMING
      cwtmp2=mpi_wtime()
#endif
!$OMP end master
!$OMP barrier

!     Mark upwind prisms for efficiency
      if(ltvd) then
!$OMP   workshare
        iupwind_e=0
!$OMP   end workshare

!$OMP   do 
        do i=1,nea
          if(itvd_e(i)==0) then
            iupwind_e(i)=1 
          else !itvd_e=1
            do j=1,i34(i)
              nd=elnode(j,i)
              toth=eta2(nd)+dp(nd)
              if(toth<h_tvd) then
                iupwind_e(i)=1; exit
              endif
            enddo !j
          endif !itvd_e
        enddo !i=1,ne
!$OMP   end do
      endif !ltvd

      do i=1,ntr
!$OMP   workshare
        flux_mod_hface(i,1:nvrt,1:ns)=flux_adv_hface(1:nvrt,1:ns)
        !flux_adv_vface from step routine. This routine cannot handle
        !settling vel. and assumes flux_adv_vface is same across all tracers
        flux_mod_vface(i,1:nvrt,1:ne)=flux_adv_vface(1:nvrt,1,1:ne)
!$OMP   end workshare
      enddo !i

!     Debug
!      do i=1,ne
!        if(idry_e(i)==1) cycle
!        do k=kbe(i)+1,nvrt
!          if(flux_mod_vert(1,k,i)<-1.d33) then
!            write(errmsg,*)'Vertical flux: out of bound',ielg(i),k,flux_mod(1,k,2,i),flux_adv_vface(k,1,i)
!            call parallel_abort(errmsg)
!          endif
!        enddo !k
!      enddo !i

!#ifdef DEBUG
!      dtbe=dt !min (over all subcycles and all levels) time step allowed at each element
!#endif

!$OMP single
      it_sub=0
      time_r=dt !time remaining
      difnum_max_l=0 !max. diffusion number reached by this process (check stability)
      surf_vol(:)=0 !sum of surface fluxes for conservation check
!$OMP end single
      loop11: do
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!$OMP single
      it_sub=it_sub+1
!$OMP end single

!     Compute flux limiters and modify fluxes
      if(ltvd) then !TVD is used for all tracers
!$OMP   workshare
        up_rat_hface=-1.d34 !flags
        up_rat_vface=-1.d34 !flags
!$OMP   end workshare

!       Vertical limiters
!$OMP   do
        do i=1,ne
          if(idry_e(i)==1) cycle

          up_rat_vface(:,:,i)=-1.d0 !initialize upwind ratio for abnormal cases; \phi=0 when r=-1
          do k=kbe(i)+1,nvrt-1 !bottom and surface flux unchanged at -1
            if(flux_adv_vface(k,1,i)<-1.d33) then
              write(errmsg,*)'Transport: Left out vertical flux (3):',i,k
              call parallel_abort(errmsg)
            endif
            if(flux_adv_vface(k,1,i)>0) then
              kup=k !upwind prism
              kdo=k+1 !downwind prism
            else
              kup=k+1 
              kdo=k
            endif

            psum=0 !sum of original fluxes
            psumtr(1:ntr)=0 !sum of products (|Q|*(T-T))
#ifdef DEBUG
            if(flux_adv_vface(kup,1,i)<-1.d33.or.flux_adv_vface(kup-1,1,i)<-1.d33) then
              write(errmsg,*)'Left out vertical flux (4):',i,kup
              call parallel_abort(errmsg)
            endif
#endif
            if(flux_adv_vface(kup,1,i)<0.and.kup/=nvrt) then
              psum=psum+abs(flux_adv_vface(kup,1,i))
              psumtr(1:ntr)=psumtr(1:ntr)+abs(flux_adv_vface(kup,1,i))*(tr_el(1:ntr,kup+1,i)-tr_el(1:ntr,kup,i))
            endif
            if(flux_adv_vface(kup-1,1,i)>0.and.kup/=kbe(i)+1) then
              psum=psum+abs(flux_adv_vface(kup-1,1,i))
              psumtr(1:ntr)=psumtr(1:ntr)+abs(flux_adv_vface(kup-1,1,i))*(tr_el(1:ntr,kup-1,i)-tr_el(1:ntr,kup,i))
            endif
            do j=1,i34(i)
              jsj=elside(j,i)
              ie=ic3(j,i)
              if(ie/=0) then; if(idry_e(ie)==0.and.kup>=kbs(jsj)+1.and.ssign(j,i)*flux_adv_hface(kup,jsj)<0) then
#ifdef DEBUG
                if(flux_adv_hface(kup,jsj)<-1.d33) then
                  write(errmsg,*)'Left out horizontal flux (5):',jsj,kup
                  call parallel_abort(errmsg)
                endif
#endif
                psum=psum+abs(flux_adv_hface(kup,jsj))
                psumtr(1:ntr)=psumtr(1:ntr)+abs(flux_adv_hface(kup,jsj))*(tr_el(1:ntr,kup,ie)-tr_el(1:ntr,kup,i))
              endif; endif
            enddo !j

! This is the calculation of the TVD stability/variation. Selection is a performance killer.
            do j=1,ntr
              tmp=(tr_el(j,kup,i)-tr_el(j,kdo,i))*abs(flux_adv_vface(k,1,i))
              if(abs(tmp)>1.e-20) up_rat_vface(j,k,i)=psumtr(j)/tmp !otherwise it remains at -1
            enddo !j

!#ifdef DEBUG
!            if( flux_lim( up_rat_vface(1,k,i))> 0.1) ntot_v=ntot_v+1
!#endif            
          enddo !k=kbe(i)+1,nvrt-1
        enddo !i=1,ne
!$OMP   end do

!       Horizontal limiters
!$OMP   do
        do i=1,ns !residents
          if(idry_s(i)==1) cycle

!         At least one element is wet
          up_rat_hface(:,:,i)=-1.d0 !initialize (for below bottom and abnormal cases)
          if(isdel(2,i)==0.or.(isdel(2,i)/=0.and.idry_e(max(1,isdel(2,i)))==1).or.idry_e(isdel(1,i))==1) cycle

!         Not bnd face; 2 elements are wet
!          kb1=min(kbe(isdel(1,i)),kbe(isdel(2,i)))
!          kb=max(kbe(isdel(1,i)),kbe(isdel(2,i)))
!          do k=kb1+1,kb-1
!            if(flux_adv_hface(k,i)/=0) then
!              write(errmsg,*)'Pls zero out the excess layers:',flux_adv_hface(k,i),i,isdel(1,i),isdel(2,i),k,kb1,kb
!              call parallel_abort(errmsg)
!            endif
!          enddo !k
 
!         Leave k=kbs unchanged
          do k=kbs(i)+1,nvrt !faces
            if(flux_adv_hface(k,i)<-1.d33) then
              write(errmsg,*)'Left out horizontal flux (3):',i,k
              call parallel_abort(errmsg)
            endif
            if(flux_adv_hface(k,i)>0) then
              iup=isdel(1,i); ido=isdel(2,i) !up/downwind prisms
            else
              iup=isdel(2,i); ido=isdel(1,i)
            endif

            psum=0
            psumtr(1:ntr)=0
            if(flux_adv_vface(k,1,iup)<-1.d33.or.flux_adv_vface(k-1,1,iup)<-1.d33) then
              write(errmsg,*)'Left out vertical flux (6):',iup,k,flux_adv_vface(k,1,iup), flux_adv_vface(k-1,1,iup)
              call parallel_abort(errmsg)
            endif
            if(flux_adv_vface(k,1,iup)<0.and.k/=nvrt) then
              psum=psum+abs(flux_adv_vface(k,1,iup))
              psumtr(1:ntr)=psumtr(1:ntr)+abs(flux_adv_vface(k,1,iup))*(tr_el(1:ntr,k+1,iup)-tr_el(1:ntr,k,iup))
            endif
            if(flux_adv_vface(k-1,1,iup)>0.and.k>kbe(iup)+1) then
              psum=psum+abs(flux_adv_vface(k-1,1,iup))
              psumtr(1:ntr)=psumtr(1:ntr)+abs(flux_adv_vface(k-1,1,iup))*(tr_el(1:ntr,k-1,iup)-tr_el(1:ntr,k,iup))
            endif

            do j=1,i34(iup)
              jsj=elside(j,iup)
              ie=ic3(j,iup) !must be inside 2-tier aug. domain
#ifdef DEBUG
              if(ie>0) then !inside 1-tier aug. domain
                !Check consistency between iegl and iegl2 etc
                if(ielg(ie)/=ielg2(ie)) call parallel_abort('TRANS:2.1')
                ind1=ielg(ie)
                if(iegl2(1,ind1)/=myrank) call parallel_abort('TRANS:2.3')
                if(iegl(ind1)%id/=iegl2(2,ind1)) call parallel_abort('TRANS:2.2')
                if(idry_e_2t(ie)/=idry_e(ie)) call parallel_abort('TRANS:2.4')
!'
              endif
#endif
              if(ie<0) then !outside 1-tier aug. domain
                ie=iabs(ie) !global elem.
                write(12,*)'inside 2-tier ghost zone:',ie
!Error: eventually into DEBUG or assert mode
                if(iegl2(1,ie)/=myrank) then
                  write(errmsg,*)'TVD: element outside:',ie
                  call parallel_abort(errmsg)
                endif
                ind1=iegl2(2,ie) !local elem. index in 2-tier aug. domain
                if(ind1<=nea.or.ind1>nea2) then
                  write(errmsg,*)'TVD: element wrong:',ind1,nea,nea2
                  call parallel_abort(errmsg)
                endif
                ie=ind1
              endif !ie<0

              !idry_e_2t, tr_el are valid up to 2-tier aug.
              if(ie/=0) then; if(idry_e_2t(ie)==0.and.k>=kbs(jsj)+1.and.ssign(j,iup)*flux_adv_hface(k,jsj)<0) then
#ifdef DEBUG
                if(flux_adv_hface(k,jsj)<-1.d33) then
                  write(errmsg,*)'Left out horizontal flux (6):',jsj,k
                  call parallel_abort(errmsg)
                endif
#endif
                psum=psum+abs(flux_adv_hface(k,jsj))
                psumtr(1:ntr)=psumtr(1:ntr)+abs(flux_adv_hface(k,jsj))*(tr_el(1:ntr,k,ie)-tr_el(1:ntr,k,iup))
              endif; endif
            enddo !j

            do j=1,ntr
              tmp=(tr_el(j,k,iup)-tr_el(j,k,ido))*abs(flux_adv_hface(k,i))
              if(abs(tmp)>1.e-20) up_rat_hface(j,k,i)=psumtr(j)/tmp
            enddo !j
!#ifdef DEBUG
!            if(flux_lim( up_rat_hface(1,k,i))>0.1) ntot_h=ntot_h+1
!#endif
          enddo !k=kbs(i)+1,nvrt
        enddo !i=1,ns
!$OMP   end do

!       Debug
!        if(it==1.and.it_sub==1) then
!          do i=1,ne
!            do j=1,3
!              jsj=elside(j,i)
!              write(99,*)isdel(1,jsj),isdel(2,jsj),up_rat()
!            enddo !j
!          enddo !i
!          stop
!        endif

!       Reset upwind ratios for upwind prism faces
!$OMP   do
        do i=1,ns
          do j=1,2
            ie=isdel(j,i)
            if(ie>0) then; if(iupwind_e(ie)/=0) then
              up_rat_hface(:,:,i)=0; exit
            endif; endif
          enddo !j
        enddo !i
!$OMP   end do

!$OMP   master
!       Exchange up_rat
#ifdef INCLUDE_TIMING
        cwtmp=mpi_wtime()
#endif
        call exchange_s3d_tr2(up_rat_hface)
        call exchange_e3d_tr2(up_rat_vface)
#ifdef INCLUDE_TIMING
        wtimer(9,2)=wtimer(9,2)+mpi_wtime()-cwtmp
#endif

#ifdef INCLUDE_TIMING
        cwtmp2=mpi_wtime()
#endif
!$OMP   end master

!       Reset upwind ratios 
!$OMP   do
        do i=1,nea
          if(iupwind_e(i)/=0) then
            up_rat_vface(:,:,i)=0
!            do j=1,i34(i) !sides
!              up_rat_hface(:,:,elside(j,i))=0
!            enddo !j
          endif
        enddo !i=1,ne
!$OMP   end do

!       Modifed fluxes flux_mod (their signs do not change) 
!       Vertical fluxes
!$OMP   do
        do i=1,ne !residents
          if(idry_e(i)==1) cycle

          do k=kbe(i)+1,nvrt-1 !leave out the bnd
!           Compute \delta_i
            if(flux_adv_vface(k,1,i)>0) then
              kup=k !upwind prism
            else
              kup=k+1
            endif

            delta_tr(1:ntr)=0.d0
            do l=0,1 !two vertical faces of upwind prism
              if(flux_adv_vface(kup-l,1,i)*(1-2*l)>0) then !outflow
                do j=1,ntr
                  rat=up_rat_vface(j,kup-l,i)
#ifdef DEBUG
                  if(rat<-1.d33) then
                    write(errmsg,*)'Left out (1):',i,kup-l,rat,it_sub,j
                    call parallel_abort(errmsg)
                  endif
#endif
                  if(abs(rat)>1.e-5) then
                    tmp=flux_lim(rat)/rat/2.d0
#ifdef DEBUG
                    if(tmp<0.or.tmp>1) then
                      write(errmsg,*)'Flux limiting failed (1):',tmp,rat,flux_adv_vface(kup-l,1,i),l,kup
                      call parallel_abort(errmsg)
                    endif 
#endif
                    delta_tr(j)=delta_tr(j)+tmp
                  endif
                enddo !j=1,ntr
              endif !outflow face
            enddo !l=0,1

            do j=1,i34(i)
              jsj=elside(j,i) !residents
              !ie=ic3(j,i)
              if(kup>=kbs(jsj)+1.and.ssign(j,i)*flux_adv_hface(kup,jsj)>0) then
                do jj=1,ntr
                  rat=up_rat_hface(jj,kup,jsj)
#ifdef DEBUG
                  if(rat<-1.d33) then
                    write(errmsg,*)'Left out (3):',i,j,kup,rat,jj
                    call parallel_abort(errmsg)
                  endif 
#endif
                  if(abs(rat)>1.d-5) then
                    tmp=flux_lim(rat)/rat/2.d0
#ifdef DEBUG
                    if(tmp<0.or.tmp>1) then
                      write(errmsg,*)'Flux limiting failed (3):',tmp,rat,jj
                      call parallel_abort(errmsg)
                    endif 
#endif
                    delta_tr(jj)=delta_tr(jj)+tmp
                  endif
                enddo !jj=1,ntr
              endif
            enddo !j

            do j=1,ntr
              flux_mod_vface(j,k,i)=flux_adv_vface(k,1,i)*(1.d0- &
     &flux_lim(up_rat_vface(j,k,i))/2.d0+delta_tr(j))
            enddo !j
          enddo !k=kbe(i)+1,nvrt-1  
        enddo !i=1,ne
!$OMP   end do

!       Horizontal fluxes
!$OMP   do
        do i=1,ns
          if(idry_s(i)==1.or.isdel(2,i)==0.or.idry_e(isdel(1,i))==1) cycle
          if(idry_e(isdel(2,i))==1) cycle

!         Both elements are wet
          do k=kbs(i)+1,nvrt
            if(flux_adv_hface(k,i)>0) then
              iup=isdel(1,i)
            else
              iup=isdel(2,i)
            endif
 
            delta_tr(1:ntr)=0
            do l=0,1 !two vertical faces of upwind prism
              if(flux_adv_vface(k-l,1,iup)*(1-2*l)>0) then !outflow
                do j=1,ntr
                  rat=up_rat_vface(j,k-l,iup)
#ifdef DEBUG
                  if(rat<-1.d33) then
                    write(errmsg,*)'Left out (5):',iup,k-l,rat,j
                    call parallel_abort(errmsg)
                  endif
#endif
                  if(abs(rat)>1.d-5) then
                    tmp=flux_lim(rat)/rat/2.d0
#ifdef DEBUG
                    if(tmp<0.or.tmp>1) then
                      write(errmsg,*)'Flux limiting failed (5):',tmp,rat,j
                      call parallel_abort(errmsg)
                    endif
#endif
                    delta_tr(j)=delta_tr(j)+tmp
                  endif
                enddo !j=1,ntr
              endif !outflow face
            enddo !l=0,1

            do j=1,i34(iup)
              jsj=elside(j,iup) !inside aug. domain
!              ie=ic3(j,iup) !not really used
              if(k>=kbs(jsj)+1.and.ssign(j,iup)*flux_adv_hface(k,jsj)>0) then !outflow
                do jj=1,ntr
                  rat=up_rat_hface(jj,k,jsj)
#ifdef DEBUG
                  if(rat<-1.d33) then
                    write(errmsg,*)'Left out (7):',iup,ielg(ie),k,rat,jj
                    call parallel_abort(errmsg)
                  endif
#endif
                  if(abs(rat)>1.e-5) then
                    tmp=flux_lim(rat)/rat/2.d0
#ifdef DEBUG
                    if(tmp<0.or.tmp>1) then
                      write(errmsg,*)'Flux limiting failed (7):',tmp,rat,jj
                      call parallel_abort(errmsg)
                    endif
#endif
                    delta_tr(jj)=delta_tr(jj)+tmp
                  endif
                enddo !jj=1,ntr
              endif !outflow
            enddo !j

            do j=1,ntr
              flux_mod_hface(j,k,i)=flux_adv_hface(k,i)*(1.d0- &
     &flux_lim(up_rat_hface(j,k,i))/2.d0+delta_tr(j)) 
            enddo !j
          enddo !k=kbs(i)+1,nvrt
        enddo !i=1,ns
!$OMP   end do

      endif !ltvd; flux limiter

!     Compute sub time step
!     Strike out \hat{S}^- (including all horizontal and vertical bnds, and where ic3(j,i) is dry)
!     Caution: \hat{S}^- conditions must be consistent later in the advective flux part!!!!!!
!     Implicit vertical flux for upwind; explicit for TVD

      if(ltvd.or.it_sub==1) then !for upwind, only compute dtb for the first step
!!$OMP   single
!        dtbl=time_r
!        ie01=0 !element # where the exteme is attained (local)
!        lev01=0 !level #
!        in_st=0 !tracer #
!!$OMP   end single

!$OMP   workshare
        dtb_min3(:)=time_r !init
!$OMP   end workshare

!$OMP   do
        do i=1,ne
          if(idry_e(i)==1) cycle

          do k=kbe(i)+1,nvrt !prism
            psumtr(1:ntr)=0.d0 !sum of modified fluxes for all inflow bnds
   
            if(ltvd.and.iupwind_e(i)==0) then !TVD for all tracers
              if(k/=nvrt.and.flux_mod_vface(1,k,i)<0) then !flux_mod and flux_adv same sign
                psumtr(1:ntr)=psumtr(1:ntr)+abs(flux_mod_vface(1:ntr,k,i))
!               Debug
!                  if(it==46.and.it_sub==1.and.i==58422) write(99,*)k,flux_adv_vface(k,1,i)
              endif
              if(k-1/=kbe(i).and.flux_mod_vface(1,k-1,i)>0) then
                psumtr(1:ntr)=psumtr(1:ntr)+abs(flux_mod_vface(1:ntr,k-1,i))
!               Debug
!                  if(it==46.and.it_sub==1.and.i==58422) write(99,*)k,flux_adv_vface(k-1,1,i)
              endif
            endif !TVD

            do j=1,i34(i)
              jsj=elside(j,i) !resident side
              ie=ic3(j,i)
              if(k>=kbs(jsj)+1) then
                ref_flux = flux_mod_hface(1,k,jsj) !flux_mod(:) same sign as flux_adv
                same_sign = (ssign(j,i)*ref_flux)<0
!!DIR$ IVDEP 
                if((ie/=0.and.idry_e(max(1,ie))==0.or.ie==0.and.isbs(jsj)>0).and.same_sign) then
                  do jj=1,ntr
#ifdef DEBUG
                    if(flux_mod_hface(jj,k,jsj)<-1.d33) then
                      write(errmsg,*)'Left out horizontal flux (10):',i,k,j,jj
                      call parallel_abort(errmsg)
                    endif
#endif

                    psumtr(jj)=psumtr(jj)+abs(flux_mod_hface(jj,k,jsj))
                  enddo !jj
!                     Debug
!                     if(it==46.and.it_sub==1.and.i==58422) write(99,*)j,k,flux_adv_hface(k,jsj)
!                   if(jj==1.and.ssign(j,i)*flux_adv_hface(k,jsj)>0) nplus=nplus+1
                endif !ie
              endif !k>=kbs
            enddo !j

            vj=area(i)*(ze(k,i)-ze(k-1,i))

!               Debug
!                if(it==46.and.it_sub==1.and.i==58422) write(99,*)k,nplus,vj

            do jj=1,ntr
              if(psumtr(jj)/=0) then
                tmp=vj/psumtr(jj)*(1-1.e-6) !safety factor included
                if(tmp<dtb_min3(i)) then
                  !dtbl=tmp 
                  dtb_min3(i)=tmp
                  !ie01=i; lev01=k; in_st=jj
                endif
!#ifdef DEBUG
!                if(tmp<dtbe(i)) dtbe(i)=tmp
!#endif
              endif
            enddo !jj

!            if(qj/=0) dtb_altl=min(dtb_altl,vj/(1+nplus)/qj*(1-1.e-10)) !safety factor included
          enddo !k=kbe(i)+1,nvrt
        enddo !i=1,ne
!$OMP   end do

!$OMP   workshare
        dtbl=minval(dtb_min3)
!$OMP   end workshare

!$OMP   master
#ifdef INCLUDE_TIMING
        cwtmp=mpi_wtime()
        timer_ns(1)=timer_ns(1)+cwtmp-cwtmp2
#endif
        buf(1,1)=dtbl; buf(2,1)=myrank
        call mpi_allreduce(buf,buf2,1,MPI_2DOUBLE_PRECISION,MPI_MINLOC,comm,ierr)
        dtb=buf2(1,1)
#ifdef INCLUDE_TIMING
        cwtmp2=mpi_wtime()
        wtimer(9,2)=wtimer(9,2)+cwtmp2-cwtmp
#endif

#ifdef DEBUG
        if(dtb<=0.or.dtb>time_r) then
          write(errmsg,*)'Transport: Illegal sub step:',dtb,time_r
          call parallel_abort(errmsg)
        endif
#endif

!       Output time step
!        if(myrank==int(buf2(2,1)).and.ie01>0) &
!     &write(12,'(a20,5(1x,i10),1x,f14.3,1x,e22.10)') &
!     &'TVD-upwind dtb info:',it,it_sub,ielg(ie01),lev01,in_st,dtb,it*dt !,dtb_alt 
!$OMP   end master

      endif !ltvd.or.it_sub==1; compute dtb

!$OMP master
      dtb=min(dtb,time_r) !for upwind
      time_r=time_r-dtb
!$OMP end master
!$OMP barrier

!     Store last step's S,T
!$OMP workshare
      trel_tmp(1:ntr,:,:)=tr_el(1:ntr,:,1:nea)
!$OMP end workshare

!$OMP do reduction(max: difnum_max_l)
      do i=1,ne
        if(idry_e(i)==1) cycle

!       Wet elements with wet nodes
!       Matrix
        ndim=nvrt-kbe(i)
        do k=kbe(i)+1,nvrt
          kin=k-kbe(i) 
          alow(kin)=0
          cupp(kin)=0
          bigv=area(i)*(ze(k,i)-ze(k-1,i)) !volume
          dtb_by_bigv = dtb/bigv
#ifdef DEBUG
          if(bigv<=0) then
            write(errmsg,*)'Negative volume (5): ',bigv,i,k
            call parallel_abort(errmsg)
          endif
#endif
          bdia(kin)=1
          if(k<nvrt) then
            av_df=sum(dfh(k,elnode(1:i34(i),i)))/i34(i)
            av_dz=(ze(k+1,i)-ze(k-1,i))/2.d0
#ifdef DEBUG
            if(av_dz<=0) then
              write(errmsg,*)'Impossible 111'
              call parallel_abort(errmsg)
            endif
#endif
            tmp=area(i)*dtb_by_bigv*av_df/av_dz
            cupp(kin)=cupp(kin)-tmp
            bdia(kin)=bdia(kin)+tmp
          endif !k<nvrt

          if(k>kbe(i)+1) then
            av_df=sum(dfh(k-1,elnode(1:i34(i),i)))/i34(i)
            av_dz=(ze(k,i)-ze(k-2,i))/2.d0
#ifdef DEBUG
            if(av_dz<=0) then
              write(errmsg,*)'Impossible 112'
              call parallel_abort(errmsg)
            endif
#endif
            tmp=area(i)*dtb_by_bigv*av_df/av_dz
            alow(kin)=alow(kin)-tmp
            bdia(kin)=bdia(kin)+tmp
          endif !k>kbe(i)+1

!         Advective flux
!         Strike out \hat{S}^- (see above)
!          psumtr(1:ntr)=0 !sum of modified fluxes at all inflow bnds 
!          delta_tr(1:ntr)=0 !sum of tracer fluxes at all inflow bnds
!         Alternative mass conservative form for the advection part (Eq. C32); contribute to rrhs
          adv_tr(1:ntr)=trel_tmp(1:ntr,k,i) 

!         2 vertical faces
#ifdef DEBUG
          if(ntr>1) then; if(flux_mod_vface(1,k,i)*flux_mod_vface(2,k,i)<0) then
            write(errmsg,*)'Left out vertical flux (0):',i,k,flux_mod_vface(1:2,k,i)
            call parallel_abort(errmsg)
          endif; endif
          do jj=1,ntr
            if(flux_mod_vface(jj,k,i)<-1.d33) then
              write(errmsg,*)'Left out vertical flux:',i,k,flux_mod_vface(jj,k,i),jj
              call parallel_abort(errmsg)
            endif
          enddo !jj
#endif

          !Upper face
          if(flux_mod_vface(1,k,i)<0) then !all flux_mod(:) same sign; inflow
            if(ltvd.and.iupwind_e(i)==0) then !TVD for all tracers
              do jj=1,ntr
!                psumtr(jj)=psumtr(jj)+abs(flux_mod_vface(jj,k,i))
                adv_tr(jj)=adv_tr(jj)-dtb_by_bigv*flux_adv_vface(k,1,i)*trel_tmp(jj,min(nvrt,k+1),i)
              enddo !jj
            else !upwind
              tmp=flux_mod_vface(1,k,i)*dtb_by_bigv !flux_mod(:) all same for upwind
              if(k==nvrt) then !downwind
                bdia(kin)=bdia(kin)+tmp
              else !upwind
                cupp(kin)=cupp(kin)+tmp
              endif
            endif
          else !outflow
            if(ltvd.and.iupwind_e(i)==0) then !TVD for all tracers
              do jj=1,ntr
                adv_tr(jj)=adv_tr(jj)-dtb_by_bigv*flux_adv_vface(k,1,i)*trel_tmp(jj,k,i)
              enddo !jj
            else !upwind
              tmp=flux_mod_vface(1,k,i)*dtb_by_bigv !flux_mod(:) all same for upwind
              bdia(kin)=bdia(kin)+tmp
            endif
          endif !flux_mod_vface

          !Lower face: note that Q_j=-flux_mod_vface!
          if(k-1/=kbe(i)) then
            if(flux_mod_vface(1,k-1,i)>0) then !inflow
              if(ltvd.and.iupwind_e(i)==0) then !TVD for all tracers
                do jj=1,ntr
!                  psumtr(jj)=psumtr(jj)+abs(flux_mod_vface(jj,k-1,i))
                  adv_tr(jj)=adv_tr(jj)+dtb_by_bigv*flux_adv_vface(k-1,1,i)*trel_tmp(jj,k-1,i)
                enddo !jj
              else !upwind
                tmp=flux_mod_vface(1,k-1,i)*dtb_by_bigv
                alow(kin)=alow(kin)-tmp
              endif
            else !outflow
              if(ltvd.and.iupwind_e(i)==0) then !TVD for all tracers
                do jj=1,ntr
                  adv_tr(jj)=adv_tr(jj)+dtb_by_bigv*flux_adv_vface(k-1,1,i)*trel_tmp(jj,k,i)
                enddo !jj
              else !upwind
                tmp=flux_mod_vface(1,k-1,i)*dtb_by_bigv
                bdia(kin)=bdia(kin)-tmp
              endif
            endif !in/out
          endif !k-1/=kbe(i)

!         Additional terms in adv_tr (Eq. C32)
          if(ltvd) then !for upwind prism, up_rat_vface=0
            if(k/=nvrt) then
              do jj=1,ntr
                adv_tr(jj)=adv_tr(jj)+dtb_by_bigv*abs(flux_adv_vface(k,1,i))*(trel_tmp(jj,k,i)- &
     &trel_tmp(jj,k+1,i))*flux_lim(up_rat_vface(jj,k,i))/2.d0
              enddo !jj
            endif
            if(k-1/=kbe(i)) then
              do jj=1,ntr
                adv_tr(jj)=adv_tr(jj)+dtb_by_bigv*abs(flux_adv_vface(k-1,1,i))*(trel_tmp(jj,k,i)- &
     &trel_tmp(jj,k-1,i))*flux_lim(up_rat_vface(jj,k-1,i))/2.d0
              enddo !jj
            endif
          endif !TVD

!         Horizontal faces
          do j=1,i34(i)
            jsj=elside(j,i) !resident side
            iel=ic3(j,i)

            if(iel/=0) then
              if(idry_e(iel)==1) cycle
              trel_tmp_outside(:)=trel_tmp(:,k,iel)
            else !bnd side
              if(isbs(jsj)<=0.or.k>=kbs(jsj)+1.and.ssign(j,i)*flux_mod_hface(1,k,jsj)>=0) cycle
       
              !Open bnd side with _inflow_; compute trel_tmp from outside and save it as trel_tmp_outside(1:ntr)
              ibnd=isbs(jsj) !global bnd #
              !Find node indices on bnd segment for the 2 nodes (for type 4 b.c.)
              nwild(1:2)=0
              do ll=1,2 !nodes
                ndo=isidenode(ll,jsj)
                do lll=1,2 !2 possible bnds
                  if(isbnd(lll,ndo)==ibnd) then
                    nwild(ll)=isbnd(-lll,ndo) !global index
                    exit
                  endif
                enddo !lll
              enddo !ll
              ind1=nwild(1); ind2=nwild(2);
#ifdef DEBUG
              if(ind1==0.or.ind2==0) then
                write(errmsg,*)'Cannot find a local index'
                call parallel_abort(errmsg)
              endif
#endif

              do jj=1,natrm
                if(ntrs(jj)<=0) cycle

                do ll=irange_tr(1,jj),irange_tr(2,jj)
                  if(itrtype(jj,ibnd)==0) then !set to be same as interior (so cancel out below)
                    trel_tmp_outside(ll)=trel_tmp(ll,k,i)
                  else if(itrtype(jj,ibnd)==1.or.itrtype(jj,ibnd)==2) then
                    trel_tmp_outside(ll)=trobc(jj,ibnd)*trth(ll,1,1,ibnd)+(1-trobc(jj,ibnd))*trel_tmp(ll,k,i)
                  else if(itrtype(jj,ibnd)==3) then
                    tmp=sum(tr_nd0(ll,k,elnode(1:i34(i),i))+tr_nd0(ll,k-1,elnode(1:i34(i),i)))/2/i34(i)
                    trel_tmp_outside(ll)=trobc(jj,ibnd)*tmp+(1-trobc(jj,ibnd))*trel_tmp(ll,k,i)
                    !trel_tmp_outside(ll)=trobc(jj,ibnd)*trel0(ll,k,i)+(1-trobc(jj,ibnd))*trel_tmp(ll,k,i)
                  else if(itrtype(jj,ibnd)==4) then
                    trel_tmp_outside(ll)=trobc(jj,ibnd)* &
     &(trth(ll,k,ind1,ibnd)+trth(ll,k,ind2,ibnd)+trth(ll,k-1,ind1,ibnd)+trth(ll,k-1,ind2,ibnd))/4+ &
     &(1-trobc(jj,ibnd))*trel_tmp(ll,k,i)
                  else
                    write(errmsg,*)'TRASNPORT: INVALID VALUE FOR ITRTYPE:',jj,ibnd
!'
                    call parallel_abort(errmsg)
                  endif !itrtype
                enddo !ll  
              enddo !jj  
            endif !iel

            if(k>=kbs(jsj)+1.and.ssign(j,i)*flux_mod_hface(1,k,jsj)<0) then !inflow
              do jj=1,ntr
#ifdef DEBUG
                if(flux_mod_hface(jj,k,jsj)<-1.d33) then
                  write(errmsg,*)'Left out horizontal flux:',i,j,k,flux_mod_hface(jj,k,jsj),jj
                  call parallel_abort(errmsg)
                endif
#endif

                !psumtr(jj)=psumtr(jj)+abs(flux_mod_hface(jj,k,jsj))
                adv_tr(jj)=adv_tr(jj)-dtb_by_bigv*ssign(j,i)*flux_adv_hface(k,jsj)*trel_tmp_outside(jj)
              enddo !jj
            else !outflow
              do jj=1,ntr
                adv_tr(jj)=adv_tr(jj)-dtb_by_bigv*ssign(j,i)*flux_adv_hface(k,jsj)*trel_tmp(jj,k,i)
              enddo !jj
            endif !in/out

            if(ltvd.and.k>=kbs(jsj)+1) then !for upwind prism, up_rat_hface=0
              do jj=1,ntr
                adv_tr(jj)=adv_tr(jj)+dtb_by_bigv*abs(flux_adv_hface(k,jsj))*(trel_tmp(jj,k,i)-trel_tmp_outside(jj))* &
     &flux_lim(up_rat_hface(jj,k,jsj))/2.d0
              enddo !jj
            endif
          enddo !j

!         Check Courant number
!new11
!          do jj=1,ntr
!            if(1-dtb_by_bigv*psumtr(jj)<0) then
!              write(errmsg,*)'Courant # condition violated:',i,k,1-dtb_by_bigv*psumtr(jj),jj
!              call parallel_abort(errmsg)
!            endif
!          enddo !jj
!new11

          rrhs(1:ntr,kin)=adv_tr(1:ntr)

!         Check consistency between 2 formulations in TVD
!            if(ltvd) then 
!              if(abs(adv_t-rrhs(1,kin))>1.e-4.or.abs(adv_s-rrhs(2,kin))>1.e-4) then
!                write(11,*)'Inconsistency between 2 TVD schemes:',i,k,adv_t,rrhs(1,kin),adv_s,rrhs(2,kin)
!                stop
!              endif
!            endif !TVD

!         Body source
          rrhs(1:ntr,kin)=rrhs(1:ntr,kin)+dtb*bdy_frc(1:ntr,k,i)

!         Horizontal diffusion
          if(ihdif/=0) then
            do j=1,i34(i) !sides
              jsj=elside(j,i) !residents
              iel=ic3(j,i)
              if(iel==0.or.idry_e(max(1,iel))==1) cycle

              nd1=isidenode(1,jsj)
              nd2=isidenode(2,jsj)
              hdif_tmp=(hdif(k,nd1)+hdif(k,nd2)+hdif(k-1,nd1)+hdif(k-1,nd2))/4
              if(k>=kbs(jsj)+1) then
                !av_h=(znl(k,nd1)-znl(k-1,nd1)+znl(k,nd2)-znl(k-1,nd2))/2.d0 !average height
                av_h=zs(k,jsj)-zs(k-1,jsj)
                if(av_h<=0) call parallel_abort('TRANSPORT: Height<=0')
                !Check diffusion number; write warning message
                difnum=dtb_by_bigv*hdif_tmp/delj(jsj)*av_h*distj(jsj)
!                if(difnum>difnum_max_l) difnum_max_l=difnum
                difnum_max_l=max(difnum_max_l,difnum)
                rrhs(1:ntr,kin)=rrhs(1:ntr,kin)+difnum*(trel_tmp(1:ntr,k,iel)-trel_tmp(1:ntr,k,i))
              endif !k>=
            enddo !j    
          endif !ihdif/=0

!         b.c.
          if(k==nvrt) rrhs(1:ntr,kin)=rrhs(1:ntr,kin)+area(i)*dtb_by_bigv*flx_sf(1:ntr,i)
          if(k==kbe(i)+1) rrhs(1:ntr,kin)=rrhs(1:ntr,kin)-area(i)*dtb_by_bigv*flx_bt(1:ntr,i)
        enddo !k=kbe(i)+1,nvrt

        call tridag(nvrt,ntr,ndim,ntr,alow,bdia,cupp,rrhs,soln,gam)
        do k=kbe(i)+1,nvrt
          kin=k-kbe(i)
          tr_el(:,k,i)=soln(:,kin)
          if(ihconsv/=0) tr_el(1,k,i)=max(tempmin,min(tempmax,soln(1,kin)))
          if(isconsv/=0) tr_el(2,k,i)=max(saltmin,min(saltmax,soln(2,kin)))

          !tr_el(1:ntr,k,i)=soln(1:ntr,kin)

#ifdef USE_SED
          do j=irange_tr(1,3),irange_tr(2,3) !1,ntr
            if(tr_el(j,k,i).lt. -1.0E-20) then
!             write(12,*)'negative sediment',i,k,tr_el(j,k,i)
              tr_el(j,k,i)=0.d0
            endif
          enddo
#endif /*USE_SED*/
        enddo !k

!       Extend
        do k=1,kbe(i)
          tr_el(:,k,i)=tr_el(:,kbe(i)+1,i)
        enddo !k
      enddo !i=1,ne
!$OMP end do

!$OMP master
!     Update ghosts
#ifdef INCLUDE_TIMING
      cwtmp=mpi_wtime()
      timer_ns(1)=timer_ns(1)+cwtmp-cwtmp2
#endif
      if(ltvd) then !extend to 2-tier aug.
        call exchange_e3d_2t_tr(tr_el)
      else !pure upwind
        call exchange_e3d_tr2(tr_el)
      endif
#ifdef INCLUDE_TIMING
      cwtmp2=mpi_wtime()
      wtimer(9,2)=wtimer(9,2)+cwtmp2-cwtmp
#endif      
!$OMP end master
!$OMP barrier

#ifdef DEBUG
!surface contrinution to conservation
      do i=1,ne
        do j=1,ntr
          if(iupwind_e(i)/=0) then !upwind; implicit
            tmp=tr_el(j,nvrt,i)
          else
            tmp=trel_tmp(j,nvrt,i)
          endif
          surf_vol(j)=surf_vol(j)+dtb*flux_adv_vface(nvrt,1,i)*tmp
        enddo !j
      enddo !i
#endif

      if(time_r<1.e-8) exit loop11
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       end do loop11

!$OMP end parallel

!     Output warning for diffusion number
      if(difnum_max_l>0.5) write(12,*)'Transport: diffusion # exceeds 0.5:',it,difnum_max_l
!'

!#ifdef DEBUG
!!     Output _estimated_ # of divisions etc.
!      if(ltvd) then 
!#ifdef INCLUDE_TIMING
!        cwtmp=mpi_wtime()
!#endif
!        call mpi_allreduce(ntot_h,ntot_hgb,1,itype,MPI_SUM,comm,ierr)
!        call mpi_allreduce(ntot_v,ntot_vgb,1,itype,MPI_SUM,comm,ierr)
!#ifdef INCLUDE_TIMING
!        wtimer(9,2)=wtimer(9,2)+mpi_wtime()-cwtmp
!#endif
!        if(myrank==0) &
!          write(16,*)'Total # of vertical and S faces limited = ',ntot_hgb,ntot_vgb
!      endif !ltvd
!#endif /*DEBUG*/
      if(myrank==0) write(17,*)it,it_sub
      
#ifdef DEBUG
!debug conservation 
      psumtr=surf_vol !0
      do i=1,ne
        if(idry_e(i)==1) cycle

        do k=kbe(i)+1,nvrt
          vol=(ze(k,i)-ze(k-1,i))*area(i)
          psumtr(1:ntr)=psumtr(1:ntr)+vol*tr_el(1:ntr,k,i)
        enddo !k
      enddo !i=1,ne
      call mpi_allreduce(psumtr,adv_tr,ntr,rtype,MPI_SUM,comm,ierr)
      if(myrank==0) write(25,*)real(time_stamp/86400),adv_tr(1:ntr), &
     &' ; after itr=2'
#endif

!     Deallocate
      deallocate(trel_tmp,flux_adv_hface,flux_mod_hface,flux_mod_vface,&
     &           up_rat_hface, up_rat_vface)

!     Debug output of time steps allowed at each element
#ifdef DEBUG
!      call schism_output_custom(istat,5,1,205,'dtbe',1,ne,dtb_min3)
!      if(myrank==0.and.istat==1) write(16,*)'done outputting dtbe.66'
#endif

      end subroutine do_transport_tvd

