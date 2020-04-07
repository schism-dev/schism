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

!
!===============================================================================
!===============================================================================
! SCHISM transport models using implicit TVD in the vertical, explicit TVD
! in the horizontal.
!
!  subroutine do_transport_tvd_imp
!  function flux_lim
!
!===============================================================================
!===============================================================================
!

!     Do upwind and TVD transport
      subroutine do_transport_tvd_imp(it,ltvd,ntr,difnum_max_l) !,nvrt1,npa1,dfh1)

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
      logical, intent(in) :: ltvd !true if TVD is used (must be for all tracers) - always true in this routine
!      character(len=2), intent(in) :: flimiter
      integer, intent(in) :: ntr !# of tracers (=ntracers)
!      integer, intent(in) :: nvrt1 !,npa1 !for dimensioning
!      real(rkind), intent(in) :: dfh1(nvrt1,npa1) 
      real(rkind), intent(out) :: difnum_max_l !max. horizontal diffusion number reached by this process (check stability)


!     Functions used
      real(rkind) :: flux_lim

      integer, parameter :: iter_back=50 !maximum iterations in implicit part (upwind is used after that)
!      real(rkind), parameter :: eps1=1e-4 !1e-9 !convergence criteria (implicit)
!      real(rkind), parameter :: eps2=1e-14 !convergence criteria (implicit)

!     Working temporary arrays in this routine
      real(rkind), allocatable :: trel_tmp(:,:,:) !tracer @ elements and half levels
      real(rkind), allocatable :: flux_adv_hface(:,:) ! original horizontal flux (the local x-driection) 
      real(rkind), allocatable :: flux_mod_hface(:,:,:) !limited advective fluxes on horizontal faces
      !weno>
      real(rkind), allocatable :: trel_tmp0(:,:,:) !tracer @ elements and half levels
      real(rkind), allocatable :: trsd_tmp(:,:,:) !tmp concentration on horizontal faces
      real(rkind), allocatable :: trsd_tmp0(:,:,:) !tmp concentration on horizontal faces, only used for message-passing
      !<weno
      real(rkind), allocatable :: up_rat_hface(:,:,:) !upwind ratios for horizontal faces
!      real(rkind), allocatable :: psum2(:,:,:)

      real(rkind) :: iupwind_e(nea) !to mark upwind prisms when TVD is used
      real(rkind) :: dtb_min3(ne),buf(2,1),buf2(2,1)
      real(rkind) :: flux_mod_v1(nvrt) !coefficient of limited advective fluxes on vertical faces (space)
      real(rkind) :: flux_mod_v2(nvrt) !coefficient of limited advective fluxes on vertical faces (time)
      real(rkind) :: rrat(nvrt) !upwind ratios for vertical faces (spatial limiter)
      real(rkind) :: srat(nvrt) !s-ratio for vertical faces (temporal limiter)
      real(rkind) :: phi(nvrt) !spatial limiter 
      real(rkind) :: bigv_m(nvrt) !prism volume
      real(rkind) :: psi1(nvrt) !time limiter 
!      real(rkind) :: psi2(nvrt) !time limiter from the current iteration
!      real(rkind) :: vdf_c1(nvrt) !coefficients related to vertical diffusive flux
!      real(rkind) :: vdf_c2(nvrt) !coefficients related to vertical diffusive flux
      real(rkind) :: r_s(nvrt),r_s0(nvrt) !local Courant number

      real(rkind) :: psumtr(ntr),delta_tr(ntr),adv_tr(ntr),tmass(ntr),h_mass_in(ntr), &
     &alow(nvrt),bdia(nvrt),cupp(nvrt),rrhs(1,nvrt),soln(1,nvrt),gam(nvrt), &
     &swild(max(3,nvrt)),swild4(3,2),trel_tmp_outside(ntr),swild5(3)
      integer :: nwild(2),ielem_elm(ne)

      integer :: istat,i,j,k,kk,l,m,khh2,ie,n1,n2,n3,n4,isd,isd0,isd1,isd2,isd3,j0,je, &
                 &nd,it_sub,ntot_v,ntot_vgb,kup,kdo,jsj,jsj0,kb, &
                 &kb1,iup,ido,ie01,lev01,in_st,ie02,lev02,in_st2,jj,ll,lll,ndim,kin,iel,ibnd, &
                 &ndo,ind1,ind2,nd1,nd2,ibio,iterK,iele_max,iterK_MAX,it_sum1,it_sum2
      real(rkind) :: vnor1,vnor2,xcon,ycon,zcon,dot1,sum1,tmp,cwtmp,toth, &
                     &time_r,psum,rat,dtbl,dtbl2,dtb,vj,av_df,av_dz,hdif_tmp, &
                     &av_h,difnum,cwtmp2,bigv,dt_by_bigv,dtb_by_bigv, &
                     &term1,term2,term6,strat1,strat2,denom, &
                     &b1,b2,b3,b4,b5,vol

      real(rkind) :: ref_flux
      logical     :: same_sign, is_land
!      logical, save :: first_call
      !weno>
      real(rkind) :: trsd(2),rctr,ang1,ui,vi,t_out,rk_coef(3) 
      real(rkind),allocatable :: wm1(:),wm2(:),wm(:) !weight
      real(rkind),allocatable :: tr_min_max(:,:) 
      logical :: iweno
      integer :: itd_weno       

      !character(72) :: ftest  ! Name of debugging file
      !integer :: lftest       ! Length of debugging file name
      !integer :: iorder !temporary identifier for order of accuracy

!     Total mass (for conservation)
      if(max_iadjust_mass_consv>0) then
        psumtr=0
        do i=1,ne
          if(idry_e(i)==1) cycle

          do k=kbe(i)+1,nvrt
            vol=(ze(k,i)-ze(k-1,i))*area(i)
            psumtr(1:ntr)=psumtr(1:ntr)+vol*tr_el(1:ntr,k,i)
          enddo !k
        enddo !i=1,ne
        call mpi_allreduce(psumtr,tmass,ntr,rtype,MPI_SUM,comm,ierr)
!        if(myrank==0) write(25,*)'mass entering transport:',real(time_stamp/86400),adv_tr(1:ntr)
      endif !max_iadjust_mass_consv

!#define weno_debug
      !<weno
      
#ifdef INCLUDE_TIMING
      cwtmp2=mpi_wtime()
#endif

      allocate(trel_tmp(ntr,nvrt,0:nea),flux_adv_hface(nvrt,nsa),tr_min_max(2,ntr), &
              &flux_mod_hface(ntr,nvrt,ns),up_rat_hface(ntr,nvrt,nsa),stat=istat) 
      trel_tmp=0.d0
!      allocate(psum2(ntr,nvrt,ne))
      if(istat/=0) call parallel_abort('Transport: fail to allocate')

!     For TVD, prepare some arrays for 2-tier ghosts
#ifdef INCLUDE_TIMING
      cwtmp=mpi_wtime()
      timer_ns(1)=timer_ns(1)+cwtmp-cwtmp2
#endif
!      if(ltvd) then
      idry_e_2t(1:ne)=idry_e(1:ne)
      call exchange_e2di_2t(idry_e_2t) !now has values up to nea2
      call exchange_e3d_2t_tr(tr_el)
!      endif !ltvd
#ifdef INCLUDE_TIMING
      wtimer(9,2)=wtimer(9,2)+mpi_wtime()-cwtmp
#endif

!'    Modify here 3D velocity for transport (for whatever reason) - make sure volume conservation is not violated
!     For rewetted elements, tr_el takes the value from last wet step

!$OMP parallel default(shared) private(j,is_land,k,vnor1,vnor2,i,nd,toth)
!     Compute (pre-limited) horizontal fluxes at all faces (vertical
!     fluxes done in schism_step)
!$OMP workshare
      flux_adv_hface=-1.d34 !flags
      iupwind_e=0
!$OMP end workshare
      
!#ifdef DEBUG
!      !-------constant velocity filed used in a few benchmarks------
!      !RUN01a,RUN03*
!      !su2=0.1d0; sv2=0.0d0
!      !RUN03g1*
!      su2=10d0; sv2=0.0d0
!      !RUN07
!      !su2=1.0d0; sv2=0.0d0
!
!      !-------constant rotational velocity filed------
!      !RUN02a
!      !do j=1,ns
!      !  n1=isidenode(1,j)
!      !  n2=isidenode(2,j)
!      !  rctr=sqrt((ynd(n2)+ynd(n1))**2+(xnd(n2)+xnd(n1))**2)/2.0d0
!      !  ang1=datan2(ynd(n2)+ynd(n1),xnd(n2)+xnd(n1))
!      !  ui=-0.002094395102393d0*rctr*sin(ang1)
!      !  vi=0.002094395102393d0*rctr*cos(ang1)
!      !  su2(:,j)=ui
!      !  sv2(:,j)=vi
!      !enddo !j
!
!      ftest='trelm_xxxx'
!      lftest=len_trim(ftest)
!      write(ftest(lftest-3:lftest),'(i4.4)') myrank
!      open(95,file=out_dir(1:len_out_dir)//ftest,status='replace')
!      ftest='order_xxxx'
!      lftest=len_trim(ftest)
!      if (it==1) then
!        write(ftest(lftest-3:lftest),'(i4.4)') myrank
!        open(96,file=out_dir(1:len_out_dir)//ftest,status='replace')
!      endif
!#endif

!     Horizontal fluxes
!$OMP do 
      do j=1,ns !resident side
        if(idry_s(j)==1) cycle
        is_land=(isdel(2,j)==0.and.isbs(j)<=0)

        do k=kbs(j)+1,nvrt
          if(is_land) then !land
            flux_adv_hface(k,j)=0.d0
          else            
            vnor1=su2(k,j)*snx(j)+sv2(k,j)*sny(j)
            vnor2=su2(k-1,j)*snx(j)+sv2(k-1,j)*sny(j)
            flux_adv_hface(k,j)=(zs(k,j)-zs(k-1,j))*distj(j)*(vnor1+vnor2)/2 !normal * area = flux (in local x-direction)

!           Debug
!           if(it==46.and.i==58422) write(99,*)j,k,vnor1,vnor2,flux_adv_hface(k,jsj)
          endif !is_land
        enddo !k=kbs(i)+1,nvrt
      enddo !j=1,ns
!$OMP end do

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
!barrier below 

!     Mark upwind prisms for efficiency
!$OMP do 
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
      enddo !i=1,nea
!$OMP end do

#ifdef USE_DVD
      !Init rkai_num. Use it to first store the value of RHS
!$OMP workshare
      rkai_num(1:ntrs(12),:,1:ne)=tr_el(irange_tr(1,12):irange_tr(2,12),:,1:ne) 
!$OMP end workshare
#endif
!$OMP end parallel


!     Init horizontal mass influx (m^3*concentration) for conservation 
      h_mass_in(1:ntr)=0 !positive in

!weno>
      if(itr_met==4) then 
!-------------------------------------------------------------------------------------
      !WENO in horizontal
!-------------------------------------------------------------------------------------

      allocate(trsd_tmp(ntr,nvrt,ns),trsd_tmp0(ntr,nvrt,ns),stat=istat)
      if(istat/=0) call parallel_abort('failed in alloc. trsd_tmp and trsd_tmp0') 
      allocate(wm(max(mnweno1,mnweno2)),wm1(max(mnweno1,mnweno2)),wm2(max(mnweno1,mnweno2)), stat=istat)
      if(istat/=0) call parallel_abort('failed in alloc. wm') 

      tr_min_max(1,1)=tempmin; tr_min_max(2,1)=tempmax;
      tr_min_max(1,2)=saltmin; tr_min_max(2,2)=saltmax;
      do i=3,ntr !todo: other tracers to be added
        tr_min_max(1,i)=0.0d0
        tr_min_max(2,i)=1.0d8
      enddo

      if (ntd_weno>1) then !additional arrays for Runge-Kutta time-stepping
        allocate(trel_tmp0(ntr,nvrt,0:nea),stat=istat) 
        if(istat/=0) call parallel_abort('Transport: fail to allocate')
        trel_tmp0=0.d0
      endif

      !RK coefficients (Shu and Osher, 1988)
      rk_coef(1)=1.d0
      if (ntd_weno==2) then
        !not implemented yet
      elseif (ntd_weno==3) then
        rk_coef(2)=0.25d0; rk_coef(3)=2.d0/3.d0
      elseif (ntd_weno==4) then
        !not implemented yet
      endif

!$OMP parallel default(shared) private(i,dtbl2,k,j,ie02,lev02,in_st2,psumtr,jsj,ie,ref_flux, &
!$OMP same_sign,vj,tmp,m,wm,sum1,wm1,wm2,b1,b2,b3,b4,b5,n1,n2,trsd,kk,bigv,dtb_by_bigv, &
!$OMP iweno,iel,ibnd,nwild,ll,ndo,lll,ind1,ind2,jj,adv_tr,trel_tmp_outside)
!'
#ifdef USE_ANALYSIS
!$OMP workshare
      dtbe=dt !min (over all subcycles and all levels) time step allowed at each element
!$OMP end workshare
#endif

!$OMP single
      it_sub=0
      time_r=dt !time remaining
!$OMP end single



      loop12: do
!$OMP   single
        it_sub=it_sub+1
!$OMP   end single


!$OMP   workshare
        dtb_min3(:)=time_r !init 
!$OMP   end workshare

!$OMP     workshare
        ielem_elm=0
!$OMP     end workshare

        if (it_sub==1) then !only compute dtb for the first step, not related to scalar field
!$OMP     do 
          do i=1,ne
            if(idry_e(i)==1) cycle
            !if(idry_e(i)==1) then 
            !  dtb_min3=huge(1.d0)
            !  cycle
            !endif
            dtbl2=huge(1.d0) !init local

            ie02=0 !element # where the exteme is attained (local)
            lev02=0 !level #
            in_st2=0
            do k=kbe(i)+1,nvrt !prism
              psumtr(1)=0.d0 !sum of horizontal fluxes for all inflow bnds
     
              do j=1,i34(i)
                jsj=elside(j,i) !resident side
                ie=ic3(j,i)

                if(k>=kbs(jsj)+1) then
                  ref_flux = flux_adv_hface(k,jsj)
                  same_sign = (ssign(j,i)*ref_flux)<0 !inflow

                  if((ie/=0.and.idry_e(max(1,ie))==0.or.ie==0.and.isbs(jsj)>0).and.same_sign) then 
#ifdef DEBUG
                      if(flux_adv_hface(k,jsj)<-1.d33) then
                        write(errmsg,*)'Left out horizontal flux (10):',i,k,j !,jj
                        call parallel_abort(errmsg)
                      endif
#endif

                    psumtr(1)=psumtr(1)+abs(flux_adv_hface(k,jsj))
                  endif !ie
                endif !k>=kbs

              enddo !j


              vj=area(i)*(ze(k,i)-ze(k-1,i))

              if(psumtr(1)/=0.d0) then
                if (iupwind_e(i)==1) then !upwind
                  tmp=vj/psumtr(1)*(1.d0-1.d-6) !safety factor included
                else !weno
                  tmp=vj/psumtr(1)*courant_weno*(1.d0-1.d-6) !safety factor 1.e-6 included
                endif

                if(tmp<dtbl2) then
                  dtbl2=tmp
                  ie02=i; lev02=k !; in_st2=jj
                endif
              endif !psumtr

            enddo !k=kbe(i)+1,nvrt

            !Hybrid ELM
            if(ielm_transport/=0) then
              if(dtbl2<dtb_min_transport) then
                !write(99,*) iplg(i)
                ielem_elm(i)=1
                dtbl2=dt !unlimited
              endif
            endif

            dtb_min3(i)=dtbl2


#ifdef USE_ANALYSIS
            dtbe(i)=dtb_min3(i) !only calculate during 1st iteration
#endif

!!WARNING: 'critical' has to be outside if; otherwise some threads are
!modifying while others may be comapring a transient 'dtbl'!!
!!!$OMP       critical
!            if(dtbl2<dtbl) then
!              dtbl=dtbl2
!              ie01=ie02; lev01=lev02; in_st=in_st2
!            endif
!!!$OMP       end critical
          enddo !i=1,ne
!$OMP     end do

!$OMP     workshare
          dtbl=minval(dtb_min3)
!$OMP     end workshare

!$OMP     master
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

!         Output time step
!          if(myrank==int(buf2(2,1)).and.ie01>0) &
!     &write(12,'(a20,5(1x,i10),1x,f14.3,1x,e22.10)') &
!     &'TVD-upwind dtb info:',it,it_sub,ielg(ie01),lev01,in_st,dtb,it*dt !,dtb_alt 
!$OMP     end master
!$OMP     barrier

        endif !if it_sub==1; compute dtb, todo: see if it's necessary to compute dtb for each sub_step

!debug>
!#ifdef DEBUG
!        write(*,*) 'rank ',myrank,': ',it, dtb,min(dtb,time_r) 
!#endif
!<debug

!$OMP   single
        dtb=min(dtb,time_r) !for upwind
        time_r=time_r-dtb
!$OMP   end  single

        if (ntd_weno>1) then
!$OMP     workshare
          trel_tmp0(1:ntr,:,1:nea)=tr_el(1:ntr,:,1:nea) !Runge-Kutta, u0
!$OMP     end workshare
        endif


        do itd_weno=1,ntd_weno !temporal discretization steps (reduces to Euler when ntd_weno=1)
          !Store last step's tracer concentrations
!$OMP     workshare
          trel_tmp(1:ntr,:,1:nea)=tr_el(1:ntr,:,1:nea)
!$OMP     end workshare

          !Note:
          !tracer concentrations are re-constructed at sides, saved as trsd_tmp(1:ntr,1:nvrt,1:ns).
          !If an inflow side is also an interface side between two ranks, it will keep the
          !initial value (0) inside the current rank; whereas the same side will be 
          !an outflow side in the neighboring rank, where
          !its concentration will be re-constructed and shared.
          !Since the flow direction determines inflow/outflow, the side communication table 
          !may be changing. The message-passing procedure implemented here
          !circumvents this dynamic table by ASSUMING the interface side are resident in 2 ranks.
          !For each interface side, its value is reconstructed in one and only one rank,
          !whereas its value in the other rank is strictly 0.0d0

!$OMP     workshare
          trsd_tmp=0.0d0 
!$OMP     end workshare

!$OMP     do     
          do ie=1,ne

            if(ip_weno==2 .and. nweno2(ie)>0 .and. isten_qual2(ie) .and. isbe(ie)==0) then !p2 weno method
              do m=1,ntr; do k=kbe(ie)+1,nvrt !!!possible optimization
                !calculate each polynomial's weight 
                wm=0.d0; sum1=0.d0; 
                wm1=0.d0; wm2=0.d0; 
                do i=1,nweno2(ie)
                  b1=dot_product(wts2(:,1,i,ie),tr_el(m,k,isten2(1:6,i,ie)))
                  b2=dot_product(wts2(:,2,i,ie),tr_el(m,k,isten2(1:6,i,ie)))
                  b3=dot_product(wts2(:,3,i,ie),tr_el(m,k,isten2(1:6,i,ie)))
                  b4=dot_product(wts2(:,4,i,ie),tr_el(m,k,isten2(1:6,i,ie)))
                  b5=dot_product(wts2(:,5,i,ie),tr_el(m,k,isten2(1:6,i,ie)))
                  b1=b1*b1+b2*b2+(4.d0*b1*b3+2*b2*b4)*fwts2(1,ie)+(4.d0*b2*b5+2*b1*b4)*fwts2(2,ie)
                  b2=(4.d0*b3*b3+b4*b4)*fwts2(3,ie)+4.d0*b4*(b3+b5)*fwts2(4,ie)+(4.d0*b5*b5+b4*b4)*fwts2(5,ie)
                  b3=(4.d0*b3*b3+b4*b4+4.d0*b5*b5)*area(ie)
                  !-----------------Hu and Shu (1999)'s method----------------
                  b4=b1+b2+b3

!Error: add PRODUCTION CPP
#ifdef DEBUG
                  if(b4<0.d0.and.abs(b4)>1.0d-50) then
                    write(errmsg,*)'b4<0',b1,b2,b3,b4
                    call parallel_abort(errmsg)
                  endif
#endif

                  wm(i)=1.0d0/((epsilon2+b4)*(epsilon2+b4))
                  !----test-------------- 
                  !wm1(i)=1.0/((epsilon2+b1+b2)*(epsilon2+b1+b2)); ! sum1=sum1+wm1(i)
                  !wm2(i)=1.0/((epsilon3+b3)*(epsilon3+b3)); 
                  !wm(i)=wm1(i)*wm2(i)
                  !---------------------- 
                  sum1=sum1+wm(i)
                enddo !i
                wm=wm/sum1

                do j=1,i34(ie) ! side
                  jsj=elside(j,ie)
                  n1=isidenode(1,jsj);   n2=isidenode(2,jsj)

                  if (ssign(j,ie)*flux_adv_hface(k,jsj)>=0.d0) then !outflow face or bnd face
                    trsd=0.d0
                    do kk=1,nquad !1 or 2 quadrature points
                      do i=1,nweno2(ie)
                        trsd(kk)=trsd(kk)+wm(i)*dot_product(wmat2(1:6,i,kk,j,ie),tr_el(m,k,isten2(1:6,i,ie)))
                      enddo
                    enddo !kk
                    trsd_tmp(m,k,jsj)=sum(trsd(1:nquad))/real(nquad,rkind) !mean value of the two quadrature points
                  endif !outward flux
                enddo !j=1,i34(ie) ! side

              enddo ; enddo ! do k=kbe(ie)+1,nvrt; do m=1,ntr
              

            elseif((ip_weno==1.or.ip_weno==2).and.nweno1(ie)>0.and.isbe(ie)==0) then !p1 weno method
              !calculate each polynomial's weight 
              do m=1,ntr; do k=kbe(ie)+1,nvrt!!!possible optimization
                wm=0.d0; sum1=0.d0
                do i=1,nweno1(ie) !for each polynomial
                  b1=dot_product(wts1(:,1,i,ie),tr_el(m,k,isten1(1:3,i,ie)))
                  b2=dot_product(wts1(:,2,i,ie),tr_el(m,k,isten1(1:3,i,ie)))
                  b1=b1*b1+b2*b2 !smoothness indicator
                  b2=1.d0/((epsilon1+b1)*(epsilon1+b1))
                  wm(i)=b2
                  sum1=sum1+b2
                  !wm(i)=1.d0
                  !sum1=sum1+1.d0

!Error: PRODUCTION
#ifdef DEBUG
                  if(.not.(wm(i)>=0.d0.or.wm(i)<0.d0)) then
                    write(errmsg,*)'wm(i)=nan',wts1(:,:,i,ie),tr_el(m,k,isten1(1:3,i,ie))
                    call parallel_abort(errmsg)
                  endif
#endif

                enddo !i
                wm=wm/sum1 !final weight

                do j=1,i34(ie) ! side
                  jsj=elside(j,ie)
                  n1=isidenode(1,jsj);   n2=isidenode(2,jsj)

                  if ((ssign(j,ie)*flux_adv_hface(k,jsj)>=0.d0) .or. isbs(jsj).ne.0) then !outflow face or bnd face
                    trsd=0.d0
                    do kk=1,nquad !1 or 2 quadrature points
                      do i=1,nweno1(ie)
                        trsd(kk)=trsd(kk)+wm(i)*dot_product(wmat1(1:3,i,kk,j,ie),tr_el(m,k,isten1(1:3,i,ie)))
                      enddo
                    enddo !kk
                    trsd_tmp(m,k,jsj)=sum(trsd(1:nquad))/real(nquad,rkind) !mean value of the two quadrature points
                  endif !outward flux
                enddo !j=1,i34(ie) ! side
              enddo ; enddo !nvrt; ntr 

            else !simple upwind method
              do j=1,i34(ie)
                jsj=elside(j,ie)
                n1=isidenode(1,jsj);   n2=isidenode(2,jsj)
                do k=kbe(ie)+1,nvrt !!!possible optimization
                  if(ssign(j,ie)*flux_adv_hface(k,jsj)>=0.d0) then !outflow face
                    trsd_tmp(:,k,jsj)=tr_el(:,k,ie) 
                  endif
                enddo  !vertical layers
              enddo !j=1,i34(ie) ! side
            endif !order of accuracy

          enddo !ie
!$OMP     end do

!$OMP     master
          !message passing
          if (nproc>1) then
            !make a copy of side concentrations
            trsd_tmp0=trsd_tmp 

            !this exchange routine swaps interface concentration between two partitions
            !trsd_tmp0 receives neighboring info, while trsd_tmp remains the same
            call exchange_s3d_tr3(trsd_tmp0,trsd_tmp) 

            !Two sets of concentrations: before message-passing (trsd_tmp) and after (trsd_tmp0).
            !(1) trsd_tmp0 is the same as trsd_tmp everywhere except on interface sides;
            !(2) For an interface side, non-0 concentration can not be present in both trsd_tmp and trsd_tmp0;


            !This obtains final values for interface sides
            trsd_tmp(:,:,iside_table)= trsd_tmp(:,:,iside_table)+trsd_tmp0(:,:,iside_table)
          endif !nproc>1
          !message passing end
!$OMP     end master
!$OMP     barrier


!         Do advection; conc @ dry elem will not be changed
!         Follows the style of explicit TVD, but using WENO correction instead of TVD correction
!$OMP     do reduction(+: h_mass_in)
          do i=1,ne
            if(idry_e(i)==1) cycle


!           Wet elements with 3 wet nodes
            do k=kbe(i)+1,nvrt


              bigv=area(i)*(ze(k,i)-ze(k-1,i)) !volume
              dtb_by_bigv = dtb/bigv
    
!             Advective flux
              if (ntd_weno==1 .or. itd_weno==1) then
                adv_tr(1:ntr)=trel_tmp(1:ntr,k,i) 
              else
                adv_tr(1:ntr)=(1.d0-rk_coef(itd_weno))*trel_tmp0(1:ntr,k,i)+rk_coef(itd_weno)*trel_tmp(1:ntr,k,i) 
              endif

              !psum=0.d0 !sum of modified fluxes for all inflow bnds
              do j=1,i34(i)
                jsj=elside(j,i) !side
                iweno=.true. !use weno reconstruction as default
                iel=ic3(j,i) !neighboring elem

                if(iel/=0) then
                  if(idry_e(iel)==1) cycle
                  trel_tmp_outside(:)=trel_tmp(:,k,iel)
                  if(iupwind_e(i)+iupwind_e(iel)>0) then !reset to upwind
                    iweno=.false.
                  endif
                else !land or open bnd side
                  !skip land bnd and outflowing open bnd, 
                  !in this case iweno is not changed from initial value (.true.), but won't be used
                  if(isbs(jsj)<=0.or.k>=kbs(jsj)+1.and.ssign(j,i)*flux_adv_hface(k,jsj)>=0.d0) then
                    !For outflow open bnd side, estimate mass in (open
                    !side cannot be interface side)
                    if(max_iadjust_mass_consv>0.and.isbs(jsj)>0) then !outflow @open bnd
                      h_mass_in(:)=h_mass_in(:)-trel_tmp(1:ntr,k,i)*flux_adv_hface(k,jsj)*dtb
                    endif

                    cycle
                  endif

                  !Open bnd side with inflow (must be wet elem) or k<kbs(jsj)+1
                  iweno=.false.  !reset to upwind
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

                  do jj=1,natrm
                    if(ntrs(jj)<=0) cycle

                    if(itrtype(jj,ibnd)==0) then !set to be same as interior (so cancel out below)
                      do ll=irange_tr(1,jj),irange_tr(2,jj)
                        trel_tmp_outside(ll)=trel_tmp(ll,k,i)
                      enddo !ll
                    else if(itrtype(jj,ibnd)==1.or.itrtype(jj,ibnd)==2) then
                      do ll=irange_tr(1,jj),irange_tr(2,jj)
                        trel_tmp_outside(ll)=trobc(jj,ibnd)*trth(ll,1,1,ibnd)+(1-trobc(jj,ibnd))*trel_tmp(ll,k,i)
                      enddo !ll
                    else if(itrtype(jj,ibnd)==3) then
                      do ll=irange_tr(1,jj),irange_tr(2,jj)
                        tmp=sum(tr_nd0(ll,k,elnode(1:i34(i),i))+tr_nd0(ll,k-1,elnode(1:i34(i),i)))/2.d0/real(i34(i),rkind)
                        trel_tmp_outside(ll)=trobc(jj,ibnd)*tmp+(1.d0-trobc(jj,ibnd))*trel_tmp(ll,k,i)
                      enddo !ll
                    else if(itrtype(jj,ibnd)==4) then
                      do ll=irange_tr(1,jj),irange_tr(2,jj)
                        trel_tmp_outside(ll)=trobc(jj,ibnd)* &
       &(trth(ll,k,ind1,ibnd)+trth(ll,k,ind2,ibnd)+trth(ll,k-1,ind1,ibnd)+trth(ll,k-1,ind2,ibnd))/4.d0+ &
       &(1.d0-trobc(jj,ibnd))*trel_tmp(ll,k,i)
                      enddo !ll
                    else
                      write(errmsg,*)'TRASNPORT: INVALID VALUE FOR ITRTYPE:',jj,ibnd
  !'
                      call parallel_abort(errmsg)
                    endif !itrtype

                  enddo !jj

                  !Tally horizontal mass in for inflow case
                  if(max_iadjust_mass_consv>0) h_mass_in(:)=h_mass_in(:)-trel_tmp_outside(1:ntr)*flux_adv_hface(k,jsj)*dtb
                endif !iel

                !!upwind part
                !if(k>=kbs(jsj)+1.and.ssign(j,i)*flux_adv_hface(k,jsj)<0) then !inflow
                !  do jj=1,ntr
                !    adv_tr(jj)=adv_tr(jj)+dtb_by_bigv*abs(flux_adv_hface(k,jsj))*(trel_tmp_outside(jj)-trel_tmp(jj,k,i))
                !  enddo
                !endif
                !!weno correction
                !if (iweno) then
                !  if(k>=kbs(jsj)+1.and.ssign(j,i)*flux_adv_hface(k,jsj)<0) then !inflow
                !    do jj=1,ntr
                !      adv_tr(jj)=adv_tr(jj)+dtb_by_bigv*abs(flux_adv_hface(k,jsj))*(trsd_tmp(jj,k,jsj)-trel_tmp_outside(jj))
                !    enddo
                !  else if(k>=kbs(jsj)+1.and.ssign(j,i)*flux_adv_hface(k,jsj)>=0) then !outflow
                !    do jj=1,ntr
                !      adv_tr(jj)=adv_tr(jj)+dtb_by_bigv*abs(flux_adv_hface(k,jsj))*(trel_tmp(jj,k,i)-trsd_tmp(jj,k,jsj) )
                !    enddo
                !  endif !inflow or outflow
                !endif !if use weno

                !alternative formulation
                if (k>=kbs(jsj)+1) then
                  if (iweno) then !weno
                    do jj=1,ntr
                      tmp=rk_coef(itd_weno)*dtb_by_bigv*( ssign(j,i)*(flux_adv_hface(k,jsj))*(trel_tmp(jj,k,i)-trsd_tmp(jj,k,jsj)) )
                      adv_tr(jj)=adv_tr(jj)+tmp
                      if (adv_tr(jj)<tr_min_max(1,jj) .or. adv_tr(jj)>tr_min_max(2,jj) ) then !reset to upwind
                        adv_tr(jj)=adv_tr(jj)-tmp !reset to previous value
                        if(ssign(j,i)*flux_adv_hface(k,jsj)<0.d0) then  !inflow 
                          adv_tr(jj)=adv_tr(jj)+dtb_by_bigv*abs(flux_adv_hface(k,jsj))*(trel_tmp_outside(jj)-trel_tmp(jj,k,i))
                        endif
                      endif
                      !when ntd_weno=1, reduces to: adv_tr(jj)=adv_tr(jj)+dtb_by_bigv*( ssign(j,i)*(flux_adv_hface(k,jsj))*(trel_tmp(jj,k,i)-trsd_tmp(jj,k,jsj)) )
                    enddo
                  else !upwind
                    if(ssign(j,i)*flux_adv_hface(k,jsj)<0) then  !inflow 
                      do jj=1,ntr
                        adv_tr(jj)=adv_tr(jj)+dtb_by_bigv*abs(flux_adv_hface(k,jsj))*(trel_tmp_outside(jj)-trel_tmp(jj,k,i))
                      enddo
                    endif
                  endif !weno or upwind
                endif

              enddo !j

              tr_el(1:ntr,k,i)=adv_tr(1:ntr) 
            enddo !k=kbe(i)+1,nvrt

            !Overwrite with ELM value if enabled.
            if(ielem_elm(i)/=0) then
              rat=time_r/dt !time ratio
              rat=max(0.d0,min(1.d0,rat))
              do k=kbe(i)+1,nvrt
                do jj=1,ntr
                  psumtr(jj)=0.d0
                  do j=1,i34(i)
                    jsj=elside(j,i)
                    n1=isidenode(1,jsj); n2=isidenode(2,jsj)
                    swild4(1,1)=0.5d0*(tr_nd(jj,k,n1)+tr_nd(jj,k,n2))
                    swild4(2,1)=0.5d0*(tr_nd(jj,k-1,n1)+tr_nd(jj,k-1,n2))
  
                    swild4(1,1)=swild4(1,1)*rat+sdbt(2+jj,max(k,kbs(jsj)),jsj)*(1-rat)
                    swild4(2,1)=swild4(2,1)*rat+sdbt(2+jj,max(k-1,kbs(jsj)),jsj)*(1-rat)
                    psumtr(jj)=psumtr(jj)+0.5d0*(swild4(1,1)+swild4(2,1)) !(sdbt(3:2+ntr,max(k,kbs(jsj)),jsj)+sdbt(3:2+ntr,max(k-1,kbs(jsj)),jsj))*0.5d0
                  enddo !j
                enddo !jj
                tr_el(1:ntr,k,i)=psumtr(:)/dble(i34(i))
              enddo !k
            endif !ielem_elm(i)/=0


!           Extend
            do k=1,kbe(i)
              tr_el(1:ntr,k,i)=tr_el(1:ntr,kbe(i)+1,i)
            enddo !k
          enddo !i=1,ne
!$OMP     end do 

!$OMP     master
!         Update ghosts
#ifdef INCLUDE_TIMING
          cwtmp=mpi_wtime()
          timer_ns(1)=timer_ns(1)+cwtmp-cwtmp2
#endif
          call exchange_e3d_2t_tr(tr_el)

#ifdef INCLUDE_TIMING
          cwtmp2=mpi_wtime()
          wtimer(9,2)=wtimer(9,2)+cwtmp2-cwtmp
#endif
!$OMP     end master
!$OMP     barrier
        enddo !itd_weno


#ifdef weno_debug
!$OMP     master
!        write(95,'(f8.1,x,360000(f15.8,x))') dt-time_r,tr_el(1,2,1:ne)
!        flush(95)
        if (myrank==0) then
          write(*,*) 'time_r=',time_r
        endif
        flush(16)
!$OMP     end master
!$OMP     barrier
#endif

        if(time_r<1.d-8) then 
          exit loop12
        endif

      end do loop12
      !write(errmsg,*)'force stop debugging'
      !call parallel_abort(errmsg)
      !if(myrank==0) write(17,*)it,it_sub

!$OMP end parallel

!     Deallocate temp. arrays
      deallocate(trsd_tmp,trsd_tmp0,wm1,wm2,wm)
!     write number of sub steps
      if(myrank==0) then
        write(17,*)it,it_sub
        flush(17)
      endif

!<weno

!-------------------------------------------------------------------------------------
      else if(itr_met==3) then !TVD in horizontal
!-------------------------------------------------------------------------------------

!$OMP parallel default(shared) private(i,k,iup,ido,psum,psumtr,j,jsj,ie,ind1,ind2,tmp, &
!$OMP delta_tr,jj,rat,ref_flux,same_sign,vj,bigv,dtb_by_bigv,adv_tr,iel,trel_tmp_outside, &
!$OMP ibnd,nwild,ll,ndo,lll,dtbl2,ie02,lev02,in_st2,n1,n2,swild4)

#ifdef USE_ANALYSIS
!$OMP workshare
      dtbe=dt !min (over all subcycles and all levels) time step allowed at each element
!$OMP end workshare
#endif

      do i=1,ntr
!$OMP   workshare
        flux_mod_hface(i,1:nvrt,1:ns)=flux_adv_hface(1:nvrt,1:ns)
!$OMP   end workshare
      enddo !i

!$OMP single
      it_sub=0
      time_r=dt !time remaining
!$OMP end single

      loop11: do
!$OMP   single
        it_sub=it_sub+1
!$OMP   end single

!       Compute flux limiters and modify fluxes
        !Use h_tvd as a flag to bypass this and other parts for efficiency
        if(h_tvd<1.d5) then
!$OMP     workshare
          up_rat_hface=-1.d34 !flags
!$OMP     end workshare

!         Horizontal limiters
!#ifdef DEBUG
!          ntot_h=0 !total # of horizontal faces that have large limiters (for 1st tracer)
!#endif

!$OMP     do 
          do i=1,ns !residents
            if(idry_s(i)==1) cycle

!           At least one element is wet
            up_rat_hface(:,:,i)=0.d0 !-1.d0 !initialize (for below bottom and abnormal cases)
            if(isdel(2,i)==0.or.(isdel(2,i)/=0.and.idry_e(max(1,isdel(2,i)))==1).or.idry_e(isdel(1,i))==1) cycle

!           Both elem r wet
            !Bypass sides with at least 1 'upwind' elem
            if(iupwind_e(isdel(1,i))/=0.or.iupwind_e(isdel(2,i))/=0) cycle

!           Leave k=kbs unchanged
            do k=kbs(i)+1,nvrt !faces
              if(flux_adv_hface(k,i)<-1.d33) then
                write(errmsg,*)'Left out horizontal flux (3):',i,k
                call parallel_abort(errmsg)
              endif
              if(flux_adv_hface(k,i)>0.d0) then
                iup=isdel(1,i); ido=isdel(2,i) !up/downwind prisms
              else
                iup=isdel(2,i); ido=isdel(1,i)
              endif

              psum=0.d0 !!sum of original fluxes
              psumtr(1:ntr)=0.d0 !sum of products (|Q|*(T-T))
              do j=1,i34(iup)
                jsj=elside(j,iup)
                ie=ic3(j,iup)
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

!Error: add PRODUCTION CPP for the following 2 checks?
#ifdef DEBUG
                  if(iegl2(1,ie)/=myrank) then
                    write(errmsg,*)'TVD: element outside:',ie
                    call parallel_abort(errmsg)
                  endif
#endif

                  ind1=iegl2(2,ie) !local elem. index in 2-tier aug. domain

#ifdef DEBUG
                  if(ind1<=nea.or.ind1>nea2) then
                    write(errmsg,*)'TVD: element wrong:',ind1,nea,nea2
                    call parallel_abort(errmsg)
                  endif
#endif

                  ie=ind1
                endif !ie<0

                !idry_e_2t, tr_el are valid up to 2-tier aug.
                if(ie/=0) then; if(idry_e_2t(ie)==0.and.k>=kbs(jsj)+1.and.ssign(j,iup)*flux_adv_hface(k,jsj)<0.d0) then
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
                if(abs(tmp)>1.d-20) up_rat_hface(j,k,i)=psumtr(j)/tmp
              enddo !j

!#ifdef DEBUG
!              if(flux_lim( up_rat_hface(1,k,i))>0.1) ntot_h=ntot_h+1
!#endif
            enddo !k=kbs(i)+1,nvrt
          enddo !i=1,ns
!$OMP     end do

!         Debug
!          if(it==1.and.it_sub==1) then
!            do i=1,ne
!              do j=1,i34(i)
!                jsj=elside(j,i)
!                write(99,*)isdel(1,jsj),isdel(2,jsj),up_rat()
!              enddo !j
!            enddo !i
!            stop
!          endif

!         Reset upwind ratios and flux_mod for upwind prism faces
!          do i=1,ne
!            if(iupwind_e(i)/=0) then
!              do j=1,i34(i) !sides
!                up_rat_hface(:,:,elside(j,i))=0
!              enddo !j
!            endif
!          enddo !i=1,ne

!$OMP     master
#ifdef INCLUDE_TIMING
          timer_ns(1)=timer_ns(1)+mpi_wtime()-cwtmp2
#endif

!         Exchange up_rat
#ifdef INCLUDE_TIMING
          cwtmp=mpi_wtime()
#endif
          call exchange_s3d_tr2(up_rat_hface)
#ifdef INCLUDE_TIMING
          wtimer(9,2)=wtimer(9,2)+mpi_wtime()-cwtmp
#endif

#ifdef INCLUDE_TIMING
          cwtmp2=mpi_wtime()
#endif
!$OMP     end master
!$OMP     barrier

!         Modified horizontal fluxes
!$OMP     do 
          do i=1,ns
            if(idry_s(i)==1.or.isdel(2,i)==0.or.idry_e(isdel(1,i))==1) cycle
            if(idry_e(isdel(2,i))==1) cycle
            !Bypass if both elem r upwind (then up_rat_hface=0 on all sides)
            if(iupwind_e(isdel(1,i))/=0.and.iupwind_e(isdel(2,i))/=0) cycle

!           Both elements are wet
            do k=kbs(i)+1,nvrt
              if(flux_adv_hface(k,i)>0) then
                iup=isdel(1,i)
              else
                iup=isdel(2,i)
              endif
 
              delta_tr(1:ntr)=0.d0
              do j=1,i34(iup)
                jsj=elside(j,iup) !inside aug. domain
!                ie=ic3(j,iup) !not really used
                if(k>=kbs(jsj)+1.and.ssign(j,iup)*flux_adv_hface(k,jsj)>0.d0) then !outflow
                  do jj=1,ntr
                    rat=up_rat_hface(jj,k,jsj)
#ifdef DEBUG
                    if(rat<-1.d33) then
                      write(errmsg,*)'Left out (7):',iup,ielg(ie),k,rat,jj
                      call parallel_abort(errmsg)
                    endif
#endif
                    if(abs(rat)>1.d-5) then
                      tmp=flux_lim(rat)/(rat*2.d0)
#ifdef DEBUG
                      if(tmp<0.d0.or.tmp>1.d0) then
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
     &flux_lim(up_rat_hface(j,k,i))/2.d0+ delta_tr(j)) 
              enddo !j
            enddo !k=kbs(i)+1,nvrt
          enddo !i=1,ns
!$OMP     end do

        endif !flux limiter; h_tvd

!       Compute sub time step
!       Strike out \hat{S}^- (including all horizontal and vertical bnds, and where ic3(j,i) is dry)
!       Caution: \hat{S}^- conditions must be consistent later in the advective flux part!!!!!!
!       Implicit vertical flux for upwind; explicit for TVD

!        if(ltvd.or.it_sub==1) then !for upwind, only compute dtb for the first step
        if(h_tvd<1.d5.or.it_sub==1) then !for upwind in entire domain, only compute dtb for the first step
!$OMP     single
          dtbl=time_r !init
          ie01=0 !element # where the exteme is attained (local)
          lev01=0 !level #
          in_st=0 !tracer #
!$OMP     end single

!$OMP     workshare
          dtb_min3(:)=time_r !init
!$OMP     end workshare

!!!$OMP     workshare
!          psum2=-1.e34
!!!$OMP     end workshare

!$OMP     workshare
          ielem_elm=0
!$OMP     end workshare

!$OMP     do 
          do i=1,ne
            if(idry_e(i)==1) cycle

            dtbl2=huge(1.d0) !init local
            ie02=0 !element # where the exteme is attained (local)
            lev02=0 !level #
            in_st2=0
            do k=kbe(i)+1,nvrt !prism
              psumtr(1:ntr)=0.d0 !sum of modified fluxes for all inflow bnds
     
              do j=1,i34(i)
                jsj=elside(j,i) !resident side
                ie=ic3(j,i)
  
                if(k>=kbs(jsj)+1) then
                  ref_flux = flux_mod_hface(1,k,jsj)
                  same_sign = (ssign(j,i)*ref_flux)<0.d0
!!DIR$ IVDEP 
                  if((ie/=0.and.idry_e(max(1,ie))==0.or.ie==0.and.isbs(jsj)>0).and.same_sign) then !flux_mod(:) same sign as flux_adv
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

!              psum2(1:ntr,k,i)=psumtr(1:ntr)

              vj=area(i)*(ze(k,i)-ze(k-1,i))

!               Debug
!                if(it==46.and.it_sub==1.and.i==58422) write(99,*)k,nplus,vj

              do jj=1,ntr
                if(psumtr(jj)/=0.d0) then
                  tmp=vj/psumtr(jj)*(1.d0-1.d-6) !safety factor included
!                  if(tmp<dtbl) then
!                    dtbl=tmp 
!                    ie01=i; lev01=k; in_st=jj
!                  endif

                  if(tmp<dtbl2) then
                    dtbl2=tmp
                    ie02=i; lev02=k; in_st2=jj

!                    if(dtbl2/vj*psumtr(jj)>1) then
!                      write(errmsg,*)'WOW9:',ielg(i),k,jj,dtbl2/vj*psumtr(jj)
!                      call parallel_abort(errmsg)
!                    endif
                  endif
                endif !psumtr
              enddo !jj

!            if(qj/=0) dtb_altl=min(dtb_altl,vj/(1+nplus)/qj*(1-1.e-10)) !safety factor included
            enddo !k=kbe(i)+1,nvrt

            if(ielm_transport/=0) then
              if(dtbl2<dtb_min_transport) then
                ielem_elm(i)=1
                dtbl2=dt !unlimited
              endif
            endif
            dtb_min3(i)=dtbl2

#ifdef USE_ANALYSIS
            if(dtbl2<dtbe(i)) dtbe(i)=dtbl2
#endif

!!WARNING: 'critical' has to be outside if; otherwise some threads are
!modifying while others may be comparing a transient 'dtbl'!!
!!!$OMP       critical
!            if(dtbl2<dtbl) then
!              dtbl=dtbl2
!              ie01=ie02; lev01=lev02; in_st=in_st2
!            endif
!!!$OMP       end critical
          enddo !i=1,ne
!$OMP     end do

!$OMP     workshare
          dtbl=minval(dtb_min3)
!$OMP     end workshare

!$OMP     master
#ifdef INCLUDE_TIMING
          cwtmp=mpi_wtime()
          timer_ns(1)=timer_ns(1)+cwtmp-cwtmp2
#endif
          buf(1,1)=dtbl; buf(2,1)=real(myrank,rkind)
          call mpi_allreduce(buf,buf2,1,MPI_2DOUBLE_PRECISION,MPI_MINLOC,comm,ierr)
          dtb=buf2(1,1)
#ifdef INCLUDE_TIMING
          cwtmp2=mpi_wtime()
          wtimer(9,2)=wtimer(9,2)+cwtmp2-cwtmp
#endif

#ifdef DEBUG
          if(dtb<=0.d0.or.dtb>time_r) then
            write(errmsg,*)'Transport: Illegal sub step:',dtb,time_r
            call parallel_abort(errmsg)
          endif
#endif

!         Output time step
!          if(myrank==int(buf2(2,1)).and.ie01>0) &
!     &write(12,'(a20,5(1x,i10),1x,f14.3,1x,e22.10)') &
!     &'TVD-upwind dtb info:',it,it_sub,ielg(ie01),lev01,in_st,dtb,it*dt !,dtb_alt 
!$OMP     end master

        endif !h_tvd.or.it_sub==1; compute dtb

!$OMP   master
        dtb=min(dtb,time_r) !for upwind
        time_r=time_r-dtb
!$OMP   end master
!$OMP   barrier

!       Store last iteration's S,T
!$OMP   workshare
        trel_tmp(1:ntr,:,1:nea)=tr_el(1:ntr,:,1:nea)
!$OMP   end workshare

!$OMP   do reduction(+: h_mass_in)
        do i=1,ne
          if(idry_e(i)==1) cycle

!         Wet elements with wet nodes
          do k=kbe(i)+1,nvrt
            bigv=area(i)*(ze(k,i)-ze(k-1,i)) !volume
            dtb_by_bigv = dtb/bigv
  
!           Advective flux
!           Strike out \hat{S}^- (see above)
            psumtr(1:ntr)=0.d0 !sum of modified fluxes at all inflow bnds 
!           Alternative mass conservative form for the advection part (Eq. C32); contribute to rrhs
            adv_tr(1:ntr)=trel_tmp(1:ntr,k,i) 

!           Horizontal faces
            do j=1,i34(i)
              jsj=elside(j,i) !resident side
              iel=ic3(j,i)

              if(iel/=0) then
                if(idry_e(iel)==1) cycle
                trel_tmp_outside(:)=trel_tmp(:,k,iel)
              else !bnd side

                if(isbs(jsj)<=0.or.k>=kbs(jsj)+1.and.ssign(j,i)*flux_mod_hface(1,k,jsj)>=0.d0) then
                  !For outflow open bnd side, estimate mass in (open
                  !side cannot be interface side)
                  if(max_iadjust_mass_consv>0.and.isbs(jsj)>0) then !outflow @open bnd
                    h_mass_in(:)=h_mass_in(:)-trel_tmp(1:ntr,k,i)*flux_adv_hface(k,jsj)*dtb
                  endif

                  cycle
                endif
       
                !Open bnd side with _inflow_ or k<kbs(jsj)+1; compute trel_tmp from outside and save it as trel_tmp_outside(1:ntr)
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
     !@         if(ind1==0.or.ind2==0) then
     !@           write(errmsg,*)'Cannot find a local index'
     !@           call parallel_abort(errmsg)
     !@        endif

                do jj=1,natrm
                  if(ntrs(jj)<=0) cycle

                  if(itrtype(jj,ibnd)==0) then !set to be same as interior (so cancel out below)
                    do ll=irange_tr(1,jj),irange_tr(2,jj)
                      trel_tmp_outside(ll)=trel_tmp(ll,k,i)
                    enddo !ll
                  else if(itrtype(jj,ibnd)==1.or.itrtype(jj,ibnd)==2) then
                    do ll=irange_tr(1,jj),irange_tr(2,jj)
                      trel_tmp_outside(ll)=trobc(jj,ibnd)*trth(ll,1,1,ibnd)+(1.d0-trobc(jj,ibnd))*trel_tmp(ll,k,i)
                    enddo
                  else if(itrtype(jj,ibnd)==3) then
                    do ll=irange_tr(1,jj),irange_tr(2,jj)
                      tmp=sum(tr_nd0(ll,k,elnode(1:i34(i),i))+tr_nd0(ll,k-1,elnode(1:i34(i),i)))/2.d0/real(i34(i),rkind)
                      trel_tmp_outside(ll)=trobc(jj,ibnd)*tmp+(1.d0-trobc(jj,ibnd))*trel_tmp(ll,k,i)
                    enddo
                  else if(itrtype(jj,ibnd)==4) then
                    do ll=irange_tr(1,jj),irange_tr(2,jj)
                      trel_tmp_outside(ll)=trobc(jj,ibnd)* &
     &(trth(ll,k,ind1,ibnd)+trth(ll,k,ind2,ibnd)+trth(ll,k-1,ind1,ibnd)+trth(ll,k-1,ind2,ibnd))/4.d0+ &
     &(1.d0-trobc(jj,ibnd))*trel_tmp(ll,k,i)
                    enddo
                  else
                    write(errmsg,*)'TRASNPORT: INVALID VALUE FOR ITRTYPE:',jj,ibnd
!'
                    call parallel_abort(errmsg)
                  endif !itrtype
                enddo !jj

                !Tally horizontal mass in fro inflow case
                if(max_iadjust_mass_consv>0) h_mass_in(:)=h_mass_in(:)-trel_tmp_outside(1:ntr)*flux_adv_hface(k,jsj)*dtb
              endif !iel

              if(k>=kbs(jsj)+1.and.ssign(j,i)*flux_mod_hface(1,k,jsj)<0.d0) then !inflow
                do jj=1,ntr
#ifdef DEBUG
                  if(flux_mod_hface(jj,k,jsj)<-1.d33) then
                    write(errmsg,*)'Left out horizontal flux:',i,j,k,flux_mod_hface(jj,k,jsj),jj
                    call parallel_abort(errmsg)
                  endif
                  if(flux_mod_hface(1,k,jsj)*flux_mod_hface(2,k,jsj)<0.d0) then
                    write(errmsg,*)'Left out horizontal flux (0):',i,j,k,flux_mod_hface(1:2,k,jsj)
                    call parallel_abort(errmsg)
                  endif
#endif
                  psumtr(jj)=psumtr(jj)+abs(flux_mod_hface(jj,k,jsj))
                  adv_tr(jj)=adv_tr(jj)+dtb_by_bigv*abs(flux_adv_hface(k,jsj))*(trel_tmp_outside(jj)-trel_tmp(jj,k,i))
                enddo !jj

#ifdef USE_DVD
                rkai_num(1:ntrs(12),k,i)=rkai_num(1:ntrs(12),k,i)+dtb_by_bigv*flux_mod_hface(irange_tr(1,12):irange_tr(2,12),k,jsj)* &
     &(trel_tmp_outside(irange_tr(1,12):irange_tr(2,12))-trel_tmp(irange_tr(1,12):irange_tr(2,12),k,i))
#endif
              endif !inflow

              if(h_tvd<1.d5.and.k>=kbs(jsj)+1) then
                do jj=1,ntr
                  adv_tr(jj)=adv_tr(jj)+dtb_by_bigv*abs(flux_adv_hface(k,jsj))*(trel_tmp(jj,k,i)-trel_tmp_outside(jj))* &
     &flux_lim( up_rat_hface(jj,k,jsj))/2.d0
                enddo !jj
              endif
            enddo !j=1,i34(i)

!           Check Courant number
            do jj=1,ntr
!              if(abs(psumtr(jj)-psum2(jj,k,i))>1.e-12) then
!                write(errmsg,*)'WOW8:',ielg(i),k,jj,psumtr(jj),psum2(jj,k,i)
!                call parallel_abort(errmsg)
!              endif

!Error: PRODUCTION
#ifdef DEBUG
              if(1.d0-dtb_by_bigv*psumtr(jj)<0.d0) then
                write(errmsg,*)'Courant # condition violated:',ielg(i),k,1-dtb_by_bigv*psumtr(jj),jj,dtb,bigv
                call parallel_abort(errmsg)
              endif
#endif
            enddo !jj

            tr_el(1:ntr,k,i)=adv_tr(1:ntr) 
          enddo !k=kbe(i)+1,nvrt

!         Check consistency between 2 formulations in TVD
!            if(ltvd) then 
!              if(abs(adv_t-rrhs(1,kin))>1.e-4.or.abs(adv_s-rrhs(2,kin))>1.e-4) then
!                write(11,*)'Inconsistency between 2 TVD schemes:',i,k,adv_t,rrhs(1,kin),adv_s,rrhs(2,kin)
!                stop
!              endif
!            endif !TVD

          !Overwrite with ELM value if enabled. The previous sections
          !are necessary for diagnostics like DVD etc
          if(ielem_elm(i)/=0) then
            rat=time_r/dt !time ratio
            rat=max(0.d0,min(1.d0,rat))
            do k=kbe(i)+1,nvrt
              do jj=1,ntr
                psumtr(jj)=0.d0
                do j=1,i34(i)
                  jsj=elside(j,i)
                  n1=isidenode(1,jsj); n2=isidenode(2,jsj)
                  swild4(1,1)=0.5d0*(tr_nd(jj,k,n1)+tr_nd(jj,k,n2))
                  swild4(2,1)=0.5d0*(tr_nd(jj,k-1,n1)+tr_nd(jj,k-1,n2))

                  swild4(1,1)=swild4(1,1)*rat+sdbt(2+jj,max(k,kbs(jsj)),jsj)*(1-rat)
                  swild4(2,1)=swild4(2,1)*rat+sdbt(2+jj,max(k-1,kbs(jsj)),jsj)*(1-rat)
                  psumtr(jj)=psumtr(jj)+0.5d0*(swild4(1,1)+swild4(2,1)) !(sdbt(3:2+ntr,max(k,kbs(jsj)),jsj)+sdbt(3:2+ntr,max(k-1,kbs(jsj)),jsj))*0.5d0
                enddo !j
              enddo !jj
              tr_el(1:ntr,k,i)=psumtr(:)/dble(i34(i))
            enddo !k
          endif !ielem_elm(i)/=0

!         Extend
          do k=1,kbe(i)
            tr_el(1:ntr,k,i)=tr_el(1:ntr,kbe(i)+1,i)
          enddo !k
        enddo !i=1,ne
!$OMP   end do

!$OMP   master
!       Update ghosts
#ifdef INCLUDE_TIMING
        cwtmp=mpi_wtime()
        timer_ns(1)=timer_ns(1)+cwtmp-cwtmp2
#endif
!        if(ltvd) then !extend to 2-tier aug.
        call exchange_e3d_2t_tr(tr_el)
!      else !pure upwind
!        call exchange_e3d_tr2(tr_el)
!      endif

#ifdef INCLUDE_TIMING
        cwtmp2=mpi_wtime()
        wtimer(9,2)=wtimer(9,2)+cwtmp2-cwtmp
#endif      
!$OMP   end master
!$OMP   barrier

        if(time_r<1.d-8) exit loop11

      end do loop11

!$OMP end parallel

      if(myrank==0) then
        write(17,*)it,it_sub
        flush(17)
      endif
!-------------------------------------------------------------------------------------
      endif !itr_met=3,4
      
!'    Save the final array from horizontal part as trel_tmp
!$OMP parallel default(shared) private(i,k,bigv_m,r_s,r_s0,m,iterK,rrat,phi,kup,kdo,psumtr, &
!$OMP tmp,flux_mod_v1,flux_mod_v2,psum,l,srat,psi1,kin,ndim,alow,bdia,cupp,dt_by_bigv,denom, &
!$OMP rrhs,soln,gam,term1,term6,strat1,strat2,bigv,av_df,av_dz,swild,j,jsj,iel,nd1,nd2,hdif_tmp,av_h,difnum)

!$OMP workshare
      trel_tmp(1:ntr,:,1:nea)=tr_el(1:ntr,:,1:nea)
!$OMP end workshare

!...  2nd step, vertical advection: Fei's addition

!      adv_tr(:)=0
!$OMP do 
      do i=1,nea 
        if(idry_e(i)==1) cycle

!       Wet elements

!--------------------Parameters that do not change through iterations------------------

        do k=kbe(i)+1,nvrt !prism
          bigv_m(k)=area(i)*(ze(k,i)-ze(k-1,i)) !volume
        enddo !k

!       Coefficients related to vertical diffusivity
!        do k=kbe(i)+1,nvrt
!          if(k<nvrt) then
!            av_df=sum(dfh(k,elnode(1:i34(i),i)))/i34(i) !diffusivity
!            av_dz=(ze(k+1,i)-ze(k-1,i))/2.d0
!            vdf_c2(k)=area(i)*dt/bigv_m(k)*av_df/av_dz !coeff. of T_{k+1}-T_k
!          endif
!
!          if(k>kbe(i)+1) then
!            av_df=sum(dfh(k-1,elnode(1:i34(i),i)))/i34(i)
!            av_dz=(ze(k,i)-ze(k-2,i))/2.d0
!            vdf_c1(k)=area(i)*dt/bigv_m(k)*av_df/av_dz !coeff. of T_{k}-T_{k-1}
!          endif
!        enddo !k

#ifdef DEBUG
        !r_s (local Courant number), only used for determining temporal limiter psi
        do k=kbe(i)+1,nvrt-1
          !flux_adv_vface(k,1,i)+flux_adv_vface(k-1,1,i) )
          r_s(k)=dt/bigv_m(k)*max(abs(flux_adv_vface(k+1,1,i)),abs(flux_adv_vface(k,1,i)),abs(flux_adv_vface(k-1,1,i)))+epsilon(1.0)
        enddo !k
        r_s(kbe(i))=r_s(kbe(i)+1)
        r_s(nvrt)=r_s(nvrt-1)
        r_s0(kbe(i)+1:nvrt)=abs(r_s(kbe(i)+1:nvrt))
        r_s0(kbe(i))=maxval(r_s0(kbe(i)+1:nvrt))
        write(12,*)'Courant #:',real(xctr(i)), real(yctr(i)),real(r_s0)
#endif
        !--------------------End: parameters that do not change through iterations------------------ 

        do m=1,ntr !cycle through tracers.
          !time limiter set to 0 before the first iteration, which will be updated after each iteration
          iterK=0
          do !iterations
            iterK=iterK+1

            !space limiter
            rrat(:)=0.d0 !init for F.S., bottom etc
            phi(:)=0.d0 !init for F.S., bottom etc (also 2D prism)
            do k=kbe(i)+1,nvrt-1 !intermediate levels (excelude bnds)
              if(flux_adv_vface(k,m,i)>=0.d0) then
                kup=k; kdo=k+1
              else !if(flux_adv_vface(k,m,i)<0) then
                kup=k+1; kdo=k
              endif
      
              psumtr(m)=0.d0 !sum of products (|Q|*(T-T))
#ifdef DEBUG
              if(flux_adv_vface(kup,m,i)<-1.d33.or.flux_adv_vface(kup-1,m,i)<-1.d33) then
                write(errmsg,*)'Left out vertical flux (4):',i,kup
                call parallel_abort(errmsg)
              endif
#endif
              if(flux_adv_vface(kup,m,i)<0.d0.and.kup/=nvrt) then !inflow at upper face (for up-upwind)
                psumtr(m)=psumtr(m)+abs(flux_adv_vface(kup,m,i))*(tr_el(m,kup+1,i)-tr_el(m,kup,i))
              endif
              if(flux_adv_vface(kup-1,m,i)>0.d0.and.kup/=kbe(i)+1) then !inflow at lower face
                psumtr(m)=psumtr(m)+abs(flux_adv_vface(kup-1,m,i))*(tr_el(m,kup-1,i)-tr_el(m,kup,i))
              endif
     
              tmp=(tr_el(m,kup,i)-tr_el(m,kdo,i))*abs(flux_adv_vface(k,m,i))
              if(abs(tmp)>1.d-20) rrat(k)=psumtr(m)/tmp !otherwise it remains at 0

              !phi(k)=max(0.d0,min(1.d0,rrat(k))) !MM
              phi(k)=max(0.d0,min(1.d0,2.0d0*rrat(k)),min(2.d0,rrat(k))) !SB
              !phi(k)=(rrat(k)+abs(rrat(k)))/(1.0d0+abs(rrat(k))) !VL
              !phi(k)=0 !upwind
            enddo !k

            !reset to upwind for upwind elem. abnormal cases
            if(iupwind_e(i)==1.or.iterK==iter_back) phi(:)=0.d0

            !The _coefficient_ of modified flux (space) at intermediate levels
            flux_mod_v1(:)=1.d0 !init
            do k=kbe(i)+1,nvrt-1 !intermediate levels (exclude bnds); k='j' in notes
              !Find downwind prism 'i'
              if(flux_adv_vface(k,m,i)<=0.d0) then
                kdo=k
              else !if(flux_adv_vface(k,m,i)>0) then
                kdo=k+1
              endif
      
              psum=0.d0
              do l=0,1 !two vertical faces of downwind prism
                if(flux_adv_vface(kdo-l,m,i)*real(1-2*l,rkind)>0.d0) then !outflow
                  if(abs(rrat(kdo-l))>1.d-6) psum=psum+phi(kdo-l)/rrat(kdo-l)
                endif !outflow
              enddo !l

              tmp=1.d0+0.5d0*(psum-phi(k))
              if(tmp<0.d0) then
                write(errmsg,*)'TRANS_IMP: mod. flux<0:',ielg(i),k,tmp
                call parallel_abort(errmsg)
              endif !tmp
              flux_mod_v1(k)=tmp
            enddo !k

            !Debug
            !write(12,*)'flux_adv_vface:',it,iterK,i,ielg(i),m,flux_adv_vface(:,m,i)
            !write(12,*)'flux_mod_v1:',it,iterK,i,ielg(i),m,flux_mod_v1

            !TVD temporal modification
            !s-ratios, defined at all levels
            srat=0.d0 !init for bottom & F.S.
            psi1=0.d0 !init for all first
            do k=kbe(i)+1,nvrt-1 !intermediate levels (excelude bnds)
              if(flux_adv_vface(k,m,i)>=0) then
                kup=k; kdo=k+1 !prisms
              else !if(flux_adv_vface(k,m,i)<0) then
                kup=k+1; kdo=k
              endif

              psum=0.d0 !sum of |Q|*(T-T)
              if(flux_adv_vface(kdo,m,i)>0.d0) then !outflow at upper face 
                psum=psum+abs(flux_adv_vface(kdo,m,i))
              endif
              if(flux_adv_vface(kdo-1,m,i)<0.d0.and.kdo/=kbe(i)+1) then !outflow at lower face
                psum=psum+abs(flux_adv_vface(kdo-1,m,i))
              endif
              psum=psum*(tr_el(m,kdo,i)-trel_tmp(m,kdo,i))

              tmp=(tr_el(m,kup,i)-trel_tmp(m,kup,i))*abs(flux_adv_vface(k,m,i))
              if(abs(tmp)>1.d-20) srat(k)=psum/tmp !otherwise it remains at 0
          
              !Prep undetermined faces
              psi1(k)=max(0.d0,min(1.d0,srat(k))) !MM
              !psi1(k)=max(0.d0,min(1.d0,2.0d0*srat(k)),min(2.d0,srat(k))) !SB
            enddo !k

            !For all levels that have a uni-directional upwind prism and (s_rat>0 or bnd), redo psi1
            do k=kbe(i),nvrt !including bnd now
              kup=0 !init
              if(flux_adv_vface(k,m,i)>=0.d0.and.k/=kbe(i)) then
                kup=k !prism
              else if(flux_adv_vface(k,m,i)<0.d0.and.k/=nvrt) then
                kup=k+1
              endif

              if(kup/=0) then; if(flux_adv_vface(kup,m,i)*flux_adv_vface(kup-1,m,i)>0.d0.and. &
     &(srat(k)>0.d0.or.k==nvrt.or.k==kbe(i))) then
                !flux_adv_vface(k,m,i)/=0 as k is one of kup | kup-1
                tmp=2.d0*(1.d0-1.d-4)*bigv_m(kup)/dt/abs(flux_adv_vface(k,m,i)) !>0
                if(tmp<=0.d0) call parallel_abort('TRANS_IMP: tmp<=0(1)')
                psi1(k)=min(1.d0,tmp) !>0
              endif; endif !uni-directional
            enddo !k
            
            !Modified flux (time) - at all levels that have a
            !uni-directional upwind prism
            flux_mod_v2(:)=-1.d20 !init as flag
            do k=kbe(i),nvrt !include bnds
              kup=0 !init
              if(flux_adv_vface(k,m,i)>=0.d0.and.k/=kbe(i)) then
                kup=k
              else if(flux_adv_vface(k,m,i)<0.d0.and.k/=nvrt) then
                kup=k+1
              endif

              if(kup/=0) then; if(flux_adv_vface(kup,m,i)*flux_adv_vface(kup-1,m,i)>0.d0) then !uni-directional
                !There is exactly 1 inflow face - k is outflow face
                kin=2*kup-1-k !inflow face
                if(abs(srat(kin))>1.d-10) then
                  psum=psi1(kin)/srat(kin) !should be >=0
                else
                  psum=0.d0
                endif
                flux_mod_v2(k)=0.5d0*(psum-psi1(k))*abs(flux_adv_vface(k,m,i))
              endif; endif !uni-directional
            enddo !k=kbe(i)+1,nvrt-1

            !Matrix
            ndim=nvrt-kbe(i) !# of eqs/unknowns
            alow=0.d0; bdia=0.d0; cupp=0.d0
            do k=kbe(i)+1,nvrt !prism
              kin=k-kbe(i) !eq. #
              dt_by_bigv = dt/bigv_m(k)

              !Diffusivity
!              if(k<nvrt) then
!                cupp(kin)=cupp(kin)-vdf_c2(k)
!                bdia(kin)=bdia(kin)+vdf_c2(k)
!              endif !k<nvrt
!              if(k>kbe(i)+1) then
!                alow(kin)=alow(kin)-vdf_c1(k)
!                bdia(kin)=bdia(kin)+vdf_c1(k)
!              endif

              !Advection part
              denom=1.d0 !denom. of Eq. (5)
              if(flux_adv_vface(k,m,i)*flux_adv_vface(k-1,m,i)>0.d0) then !uni-directional
                if(flux_adv_vface(k,m,i)>0.d0) then !outflow at upper face (including rising F.S.)
                  if(flux_mod_v2(k)<-1.d10) then !check if flux_mod_v2 has valid values
                    write(errmsg,*)'TRANS_IMP, flux_mod_v2(1):',it,iterK,m,ielg(i)
                    call parallel_abort(errmsg)
                  endif

                  denom=denom+dt_by_bigv*flux_mod_v2(k)
                endif !outflow at upper face
                !if(k-1/=kbe(i).and.flux_adv_vface(k-1,m,i)<0) then !outflow at lower
                if(flux_adv_vface(k-1,m,i)<0.d0) then !outflow at lower (including sinking bottom)
                  if(flux_mod_v2(k-1)<-1.d10) then
                    write(errmsg,*)'TRANS_IMP, flux_mod_v2(2):',it,iterK,m,ielg(i)
                    call parallel_abort(errmsg)
                  endif
                  denom=denom+dt_by_bigv*flux_mod_v2(k-1)
                endif !outflow at lower
              endif !uni-directional

              if(denom<=0.d0) then
                write(errmsg,*)'TRANS_IMP, mod.  flux<=0:',it,iterK,m,ielg(i),k,denom
                call parallel_abort(errmsg)
              endif !denom

              !Reset to upwind for upwind elem. or abnormal case
              if(iupwind_e(i)==1.or.iterK==iter_back) then
                denom=1.d0
                psi1(:)=0.d0 !reset for conservation check
              endif

              !DEBUG; new23
              !denom=1.d0

              bdia(kin)=bdia(kin)+denom
              rrhs(1,kin)=trel_tmp(m,k,i)*denom !# of columns=1 because tracer loop is outer

              if(k/=nvrt.and.flux_adv_vface(k,m,i)<0.d0) then !inflow at upper face
                tmp=dt_by_bigv*abs(flux_adv_vface(k,m,i)*flux_mod_v1(k)) !flux_mod_v1>=0
                cupp(kin)=cupp(kin)-tmp
                bdia(kin)=bdia(kin)+tmp
              endif !inflow at upper face

              if(k-1/=kbe(i).and.flux_adv_vface(k-1,m,i)>0.d0) then !inflow at lower
                tmp=dt_by_bigv*abs(flux_adv_vface(k-1,m,i)*flux_mod_v1(k-1))
                alow(kin)=alow(kin)-tmp
                bdia(kin)=bdia(kin)+tmp
              endif !inflow at lower

            enddo !k=kbe(i)+1,nvrt

            !Other RHS
            !Debug
            !write(12,*)'RHS:',it,iterK,i,ielg(i),m,rrhs(1,1:ndim)
            !write(12,*)'alow:',alow
            !write(12,*)'bdia:',bdia
            !write(12,*)'cupp:',cupp

            call tridag(nvrt,1,ndim,1,alow,bdia,cupp,rrhs,soln,gam)

            !check convergence, based on increment
            term1=sqrt(sum((soln(1,1:ndim)-tr_el(m,kbe(i)+1:nvrt,i))**2.d0))
            term6=sqrt(sum(tr_el(m,kbe(i)+1:nvrt,i)**2.d0))
            !update concentration
            tr_el(m,kbe(i)+1:nvrt,i)=soln(1,1:ndim)

            !Debug
            !write(12,*)'soln:',it,iterK,i,ielg(i),m,tr_el(m,kbe(i)+1:nvrt,i)
            !write(12,*)'diff:',term1,eps1_tvd_imp*term6+eps2_tvd_imp,term6

            !Calculate strat in prep for exit
            if(iterK==iter_back-1) strat1=maxval(soln(1,1:ndim))-minval(soln(1,1:ndim))

            !Done upwind for abnormal cases and exit
            if(iterK==iter_back.or.term1<=eps1_tvd_imp*term6+eps2_tvd_imp) then
              !strat2=maxval(soln(1,1:ndim))-minval(soln(1,1:ndim))
              !DEBUG
              !write(12,*)'TRANS_IMP, strat loss:',real(strat1),real(strat2),real(strat2-strat1),ielg(i),m,it
#ifdef USE_DVD
!Error: did not add time limiter yet due to complication @F.S.
              if(i<=ne.and.m>=irange_tr(1,12).and.m<=irange_tr(2,12)) then
                l=m-irange_tr(1,12)+1
                do k=kbe(i)+1,nvrt !prism
                  if(k/=nvrt.and.flux_adv_vface(k,m,i)<0.d0) then !inflow at upper face
                    rkai_num(l,k,i)=rkai_num(l,k,i)+dt_by_bigv* &
     &abs(flux_adv_vface(k,m,i))*(tr_el(m,k+1,i)-tr_el(m,k,i))*(1.d0-0.5d0*phi(k))
                  else if(k/=nvrt.and.flux_adv_vface(k,m,i)>0.d0) then !outflow at upper face
                    rkai_num(l,k,i)=rkai_num(l,k,i)-0.5d0*dt_by_bigv* &
     &abs(flux_adv_vface(k,m,i))*(tr_el(m,k+1,i)-tr_el(m,k,i))*phi(k)
                  endif

                  if(k-1/=kbe(i).and.flux_adv_vface(k-1,m,i)>0.d0) then !inflow at lower face
                    rkai_num(l,k,i)=rkai_num(l,k,i)+dt_by_bigv* &
     &abs(flux_adv_vface(k-1,m,i))*(tr_el(m,k-1,i)-tr_el(m,k,i))*(1.d0-0.5d0*phi(k-1))
                  else if(k-1/=kbe(i).and.flux_adv_vface(k-1,m,i)<0.d0) then !outflow at lower
                    rkai_num(l,k,i)=rkai_num(l,k,i)-0.5d0*dt_by_bigv* &
     &abs(flux_adv_vface(k-1,m,i))*(tr_el(m,k-1,i)-tr_el(m,k,i))*phi(k-1)
                  endif
               enddo !k
             endif !i<=ne etc
#endif /*USE_DVD*/

              exit
            endif !iterK

!            if(term1<=eps1_tvd_imp*term6+eps2_tvd_imp) then
!              !DEBUG
!              !write(12,*) "converged in ", iterK,i,ielg(i),m,it
!              exit
!            endif   
          enddo !nonlinear iteration

!          if(iterK>=iterK_MAX) then
!            iterK_MAX=iterK
!            iele_max=ielg(i)
!          endif   
          !it_sum1=it_sum1+iterK

          !Tally mass flux @ F.S.
!          if(i<=ne) adv_tr(m)=adv_tr(m)+dt*flux_adv_vface(nvrt,m,i)*(tr_el(m,nvrt,i)* &
!     &(1-psi1(nvrt)/2)+psi1(nvrt)/2*trel_tmp(m,nvrt,i))

        enddo !m: tracers

!       Extend
        do k=1,kbe(i)
          tr_el(:,k,i)=tr_el(:,kbe(i)+1,i)
        enddo !k

      enddo !i=1,nea
!$OMP end do

!!!$OMP master
!!      call mpi_reduce(iterK_MAX,jj,1,itype,MPI_MAX,0,comm,ierr)
!      call mpi_reduce(it_sum1,it_sum2,1,itype,MPI_SUM,0,comm,ierr)
!      if(myrank==0) write(20,*)it,real(it_sum2)/ne_global/ntr !,jj
!!!$OMP end master

!     Finalize DVD
#ifdef USE_DVD
!$OMP workshare
      !Dim= [C^2]/sec
      rkai_num(1:ntrs(12),:,1:ne)=(rkai_num(1:ntrs(12),:,1:ne)-tr_el(irange_tr(1,12):irange_tr(2,12),:,1:ne))/dt
!$OMP end workshare
#endif /*USE_DVD*/

!$OMP master
!     conservation
      if(max_iadjust_mass_consv>0) then
        psumtr=0.d0 
        do i=1,ne
          if(idry_e(i)==1) cycle

          do k=kbe(i)+1,nvrt
            vol=(ze(k,i)-ze(k-1,i))*area(i)
            psumtr(1:ntr)=psumtr(1:ntr)+vol*tr_el(1:ntr,k,i)
          enddo !k
        enddo !i=1,ne
        call mpi_allreduce(psumtr,adv_tr,ntr,rtype,MPI_SUM,comm,ierr)
        !Mass 'error' after advection step
        total_mass_error=adv_tr-tmass-h_mass_in 
!        if(myrank==0) write(25,*)'mass after 2 steps with correction:',real(time_stamp/86400),adv_tr(1:ntr)+total_mass_error(1:ntr)
      endif
!$OMP end master

!     3rd step: non-advection terms 
!     Save the final array from previous step as trel_tmp
!$OMP workshare
      trel_tmp(1:ntr,:,1:nea)=tr_el(1:ntr,:,1:nea)
!$OMP end workshare

!$OMP single
      difnum_max_l=0.d0 !max. diffusion number reached by this process (check stability)
!$OMP end single

!$OMP do reduction(max: difnum_max_l) 
      do i=1,ne
        if(idry_e(i)==1) cycle

!       Wet elements with 3 wet nodes
        ndim=nvrt-kbe(i) !# of eqs/unknowns
        do m=1,ntr !cycle through tracers
          ! Vertical movement of POM (by Richard Hofmeister)
          ! for vertically varying velocities of vertical movement.
          ! If iwsett=0 (default), set wsett(nvrt|kbe)=0 (actually bypassed below with if) and use wsett(kbe(i)+1:nvrt-1) for
          ! vertical velocities at whole levels
          ! Limit wsett to avoid char line out of boundary 
          do k=kbe(i),nvrt 
            tmp=0.9*(ze(nvrt,i)-ze(kbe(i),i))/dt !>0; safety factor added
!Error: to assert
            if(tmp<=0) call parallel_abort('Transport:tmp<=0')
            swild(k)=max(-tmp,min(tmp,wsett(m,k,i)))
          enddo !k

          !Matrix
          alow=0.d0; bdia=1.d0; cupp=0.d0
          do k=kbe(i)+1,nvrt !prism
            kin=k-kbe(i) !eq. #
            bigv=area(i)*(ze(k,i)-ze(k-1,i)) !volume
            dt_by_bigv = dt/bigv
  
            if(k<nvrt) then
              if(itur==5.and.m>=irange_tr(1,5).and.m<=irange_tr(1,5)) then !1018:itur==5
                av_df=sum(dfhm(k,m-irange_tr(1,5)+1,elnode(1:i34(i),i)))/real(i34(i),rkind) !1007
              else
                av_df=sum(dfh(k,elnode(1:i34(i),i)))/real(i34(i),rkind) !diffusivity
              endif              
              av_dz=(ze(k+1,i)-ze(k-1,i))/2.d0
              tmp=area(i)*dt_by_bigv*av_df/av_dz
              cupp(kin)=cupp(kin)-tmp
              bdia(kin)=bdia(kin)+tmp

              !if(wsett(m,k,i)<=0.0d0) then !upwinding for conc
              if(swild(k)<=0.0d0) then !upwinding for conc
                bdia(kin) = bdia(kin) - area(i)*dt_by_bigv*swild(k) !wsett(m,k,i)
              else
                cupp(kin) = cupp(kin) - area(i)*dt_by_bigv*swild(k) !wsett(m,k,i)
              endif
            endif !k<nvrt

            if(k>kbe(i)+1) then
              if(itur==5.and.m>=irange_tr(1,5).and.m<=irange_tr(1,5)) then !1018:itur==5
                av_df=sum(dfhm(k-1,m-irange_tr(1,5)+1,elnode(1:i34(i),i)))/real(i34(i),rkind) !1007
              else
                av_df=sum(dfh(k-1,elnode(1:i34(i),i)))/real(i34(i),rkind)
              endif
              av_dz=(ze(k,i)-ze(k-2,i))/2.d0
              tmp=area(i)*dt_by_bigv*av_df/av_dz
              alow(kin)=alow(kin)-tmp
              bdia(kin)=bdia(kin)+tmp

              !if(wsett(m,k-1,i)<=0.0d0) then
              if(swild(k-1)<=0.0d0) then
                alow(kin) = alow(kin) + area(i)*dt_by_bigv*swild(k-1) !wsett(m,k-1,i)
              else
                bdia(kin) = bdia(kin) + area(i)*dt_by_bigv*swild(k-1) !wsett(m,k-1,i)
              endif
            endif !k>

            !Extra terms for sediment at bottom
            if(iwsett(m)==1.and.k==kbe(i)+1) then
              !if(wsett(m,k-1,i)<0.0d0) then
              if(swild(k-1)<0.0d0) then
                write(errmsg,*)'TRAN_IMP: wsett<0,',m,k,ielg(i),swild(k-1) !wsett(m,k-1,i)
                call parallel_abort(errmsg)   
              endif
              bdia(kin)=bdia(kin)+area(i)*dt_by_bigv*swild(k-1) !wsett(m,k-1,i)
            endif !k==kbe(i)+1
         
            !# of column=1 as tracer loop is outside
            rrhs(1,kin)=trel_tmp(m,k,i)
            !b.c.
            if(k==nvrt) then
              rrhs(1,kin)=rrhs(1,kin)+area(i)*dt_by_bigv*flx_sf(m,i)
            endif !k==nvrt

            if(k==kbe(i)+1) then
              !NOTE: with settling vel., flx_bt=D-E-w_s*T_{kbe+1}, since
              !in well-formulated b.c., D \approx w_s*T_{kbe+1}. D&E are
              !deposi. & erosional fluxes respectively
              rrhs(1,kin)=rrhs(1,kin)-area(i)*dt_by_bigv*flx_bt(m,i)
            endif !k==kbe(i)+1

            !Body source
            rrhs(1,kin)=rrhs(1,kin)+dt*bdy_frc(m,k,i)

            !Horizontal diffusion
            if(ihdif/=0) then
              do j=1,i34(i) !sides
                jsj=elside(j,i) !residents
                iel=ic3(j,i)
                if(iel==0.or.idry_e(max(1,iel))==1) cycle
 
                nd1=isidenode(1,jsj)
                nd2=isidenode(2,jsj)
                hdif_tmp=(hdif(k,nd1)+hdif(k,nd2)+hdif(k-1,nd1)+hdif(k-1,nd2))/4.d0
                if(k>=kbs(jsj)+1) then
                  !av_h=(znl(k,nd1)-znl(k-1,nd1)+znl(k,nd2)-znl(k-1,nd2))/2.d0
                  !!average height
                  av_h=zs(k,jsj)-zs(k-1,jsj)
                  if(av_h<=0.d0) call parallel_abort('TRAN_IMP: av_h<=0')
                  !Check diffusion number; write warning message
                  difnum=dt_by_bigv*hdif_tmp/delj(jsj)*av_h*distj(jsj)
!                  if(difnum>difnum_max_l) difnum_max_l=difnum
                  difnum_max_l=max(difnum_max_l,difnum)
                  rrhs(1,kin)=rrhs(1,kin)+difnum*(trel_tmp(m,k,iel)-trel_tmp(m,k,i))
                endif !k>=
              enddo !j
            endif !ihdif/=0

          enddo !k=kbe(i)+1,nvrt

          call tridag(nvrt,1,ndim,1,alow,bdia,cupp,rrhs,soln,gam)

          do k=kbe(i)+1,nvrt
            kin=k-kbe(i)
            tr_el(m,k,i)=soln(1,kin)
          enddo !k

!          !171217 cal. depo_mss
!          if(m>=irange_tr(1,5).and.m<=irange_tr(2,5)) then
!            if(iwsett(m)==2.and.wsett(m,kbe(i),i).GE.0.0d0) then 
!              depo_tvd(m-irange_tr(1,5)+1,i)=dt/area(i)*area(i)*wsett(m,kbe(i),i)*tr_el(m,kbe(i)+1,i)
!            endif
!          endif
        enddo !m: tracers

        !Post-proc
        do k=kbe(i)+1,nvrt
          if(ihconsv/=0) tr_el(1,k,i)=max(tempmin,min(tempmax,tr_el(1,k,i)))
          if(isconsv/=0) tr_el(2,k,i)=max(saltmin,min(saltmax,tr_el(2,k,i)))

!#ifdef USE_SED
!          do j=irange_tr(1,5),irange_tr(2,5) !1,ntr
!            if(tr_el(j,k,i).lt.-1.0E-20) then
!!             write(12,*)'negative sediment',i,k,tr_el(j,k,i)
!              tr_el(j,k,i)=0.d0
!            endif
!          enddo
!#endif /*USE_SED*/

        enddo !k

!       Extend
        do k=1,kbe(i)
          tr_el(:,k,i)=tr_el(:,kbe(i)+1,i)
        enddo !k

      enddo !i=1,ne
!$OMP end do

!$OMP end parallel

!     Update ghosts
#ifdef INCLUDE_TIMING
      cwtmp=mpi_wtime()
      timer_ns(1)=timer_ns(1)+cwtmp-cwtmp2
#endif
      call exchange_e3d_2t_tr(tr_el)
#ifdef INCLUDE_TIMING
      cwtmp2=mpi_wtime()
      wtimer(9,2)=wtimer(9,2)+cwtmp2-cwtmp
#endif

!     Output warning for diffusion number
      if(difnum_max_l>0.5d0) write(12,*)'TRAN_IMP, diffusion # exceeds 0.5:',it,difnum_max_l
!'

#ifdef DEBUG
!debug conservation
      psumtr=0.d0
      do i=1,ne
        if(idry_e(i)==1) cycle

        do k=kbe(i)+1,nvrt
          vol=(ze(k,i)-ze(k-1,i))*area(i)
          psumtr(1:ntr)=psumtr(1:ntr)+vol*tr_el(1:ntr,k,i)
        enddo !k
      enddo !i=1,ne
      call mpi_allreduce(psumtr,adv_tr,ntr,rtype,MPI_SUM,comm,ierr)
      if(myrank==0) write(25,*)'mass after all steps:',real(time_stamp/86400.d0),adv_tr(1:ntr)
#endif

!     Deallocate temp. arrays
      deallocate(trel_tmp,flux_adv_hface,flux_mod_hface,up_rat_hface)
!      deallocate(psum2)

      end subroutine do_transport_tvd_imp

!===============================================================================
!     Flux limiter functions used in TVD schemes
!===============================================================================
      function flux_lim(ss)
      use schism_glbl, only : rkind,errmsg
      use schism_msgp, only : parallel_abort
      implicit none
   
      real(rkind) :: flux_lim
      real(rkind), intent(in) :: ss !upwind ratio
!      character(len=2), intent(in) :: flimiter

#ifdef TVD_SB
      !Superbee
      flux_lim=max(0.d0,min(1.d0,2.0d0*ss),min(2.d0,ss))
#elif TVD_VL
      !Van Leer
      flux_lim=(ss+abs(ss))/(1.0d0+abs(ss))
#elif TVD_MM
      !MINMOD
      flux_lim=max(0.d0,min(1.d0,ss))
#elif TVD_OS
      !OSHER
      flux_lim=max(0.d0,min(2.d0,ss))
#else
      write(errmsg,*)'TVD limiter not defined'
      call parallel_abort(errmsg)
#endif

      end function flux_lim

