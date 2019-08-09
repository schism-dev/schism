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
! This file collect subroutines implementing FE-FCT
! advection scheme by Loehner et al.
! There is a tunable paremeter gamma in ts_solve_low_order and fem_fct.
! Increasing it leads to positivity preserving solution.

! Driving routine is fct_ice_solve. It calles other routines
! that do low-order and figh order solutions and then combine them in a flux
! corrected way. Taylor-Galerkin type of advection is used.
!===============================================================================
!===============================================================================
! subroutine  ice_fct

  subroutine ice_fct
  use schism_glbl, only: rkind,nea,np,npa,elnode,nnp,indnd,time_stamp,rnday, &
 &fdb,lfdb,xnd,ynd,iplg
  use schism_msgp, only: nproc,myrank,parallel_abort,exchange_p2d
  use ice_module

  implicit none
 
  integer :: i,j,m,n,nd,nd2,icoef(3,3)
  real(rkind) :: um,vm,diff,sum1,flux,ae
  real(rkind) :: dx(3),dy(3),entries(3),icefluxes(3,nea),tmax(np),tmin(np), &
 &icepplus(npa),icepminus(npa),swild(npa) !,swild2(nea)
  real(rkind),allocatable :: fct_rhs(:,:),ice_tr_lo(:,:),d_ice_tr(:,:)

  allocate(fct_rhs(ntr_ice,npa),ice_tr_lo(ntr_ice,npa),d_ice_tr(ntr_ice,npa),stat=i)
  if(i/=0) call parallel_abort('ice_fct: alloc(1)')

! RHS of Taylor Galerkin
  ! Computes the rhs in a Taylor-Galerkin way (with upwind type of
  ! correction for the advection operator)
  fct_rhs=0 !-A_{jk}*a_k^n in notes; [m^2]
  !fct_rhs invalid at ghost nodes!!!
  do i=1,nea !assembling rhs over elements 
    !derivatives of shape functions
    dx=bafux(:,i)
    dy=bafuy(:,i)
    um=sum(u_ice(elnode(1:3,i)))
    vm=sum(v_ice(elnode(1:3,i)))
    !diffusivity: diff. coeff. is scaled as (area/scalevol)^(1/2) - not used
    do j=1,3
      nd=elnode(j,i)
      do m=1,3
        nd2=elnode(m,i)
        ![entries]=m^2
        entries(m)=voltriangle(i)*dt_ice*((dx(j)*(um+u_ice(nd2))+dy(j)*(vm+v_ice(nd2)))/12.0- &
  &0.5*dt_ice*(um*dx(j)+vm*dy(j))*(um*dx(m)+vm*dy(m))/9.0)
  !&diff*(dx(j)*dx(m)+dy(j)*dy(m))-0.5*dt_ice*(um*dx(j)+vm*dy(j))*(um*dx(m)+vm*dy(m))/9.0)
      enddo !m

      do m=1,ntr_ice
        fct_rhs(m,nd)=fct_rhs(m,nd)+sum(entries*ice_tr(m,elnode(1:3,i)))
      enddo !m
    enddo !j
  enddo !i=1,nea

  !Debug
!  do m=1,ntr_ice
!    swild=fct_rhs(m,:)
!    call exchange_p2d(swild)
!    fct_rhs(m,:)=swild
!  enddo !m
!  do i=1,np
!    write(12,*)'RHS:',it_main,iplg(i),xnd(i),ynd(i),fct_rhs(1:2,i)
!  enddo !i

!----------------------------------------------------------------------
! High-order solve 1st as it uses low-order arrays as temp storage
  !Taylor-Galerkin solution
  !the first approximation
  do i=1,np
    d_ice_tr(:,i)=fct_rhs(:,i)/lump_ice_matrix(i)
  enddo
  do m=1,ntr_ice
    swild=d_ice_tr(m,:)
    call exchange_p2d(swild)
    d_ice_tr(m,:)=swild
  enddo !m

  !Iterate. Use ice_tr_lo as temp
  do n=1,niter_fct-1
    do i=1,np
      do m=1,ntr_ice
        sum1=ice_matrix(0,i)*d_ice_tr(m,i)
        do j=1,nnp(i)
          nd=indnd(j,i)
          sum1=sum1+ice_matrix(j,i)*d_ice_tr(m,nd)
        enddo !j
        !Update d_ice_tr
        !lump_ice_matrix/=0
        ice_tr_lo(m,i)=d_ice_tr(m,i)+(fct_rhs(m,i)-sum1)/lump_ice_matrix(i)
      enddo !m
    enddo !i

    d_ice_tr=ice_tr_lo
    do m=1,ntr_ice
      swild=d_ice_tr(m,:)
      call exchange_p2d(swild)
      d_ice_tr(m,:)=swild
    enddo !m
  enddo !n

!----------------------------------------------------------------------
 ! Low-order solution
 ! It is assumed that m_ice, a_ice and m_snow from the previous time step
 ! are known at 1:myDim_nod2D+eDim_nod2D.
 ! One adds diffusive contribution to the rhs. It is realized as
 ! difference between the  consistent and lumped mass matrices
 ! acting on the field from the previous time step. The mass matrix on the
 ! lhs is replaced with lumped one.
  do i=1,np
    do m=1,ntr_ice
      sum1=ice_matrix(0,i)*ice_tr(m,i)
      do j=1,nnp(i)
        nd=indnd(j,i)
        sum1=sum1+ice_matrix(j,i)*ice_tr(m,nd)
      enddo !j
      ice_tr_lo(m,i)=(fct_rhs(m,i)+ice_gamma_fct*sum1)/lump_ice_matrix(i)+(1-ice_gamma_fct)*ice_tr(m,i)
    enddo !m
  enddo !i=1,np

  do m=1,ntr_ice
    swild=ice_tr_lo(m,:)
    call exchange_p2d(swild)
    ice_tr_lo(m,:)=swild
  enddo !m

  !Debug
!  do i=1,npa
!    write(12,*)'LOw:',iplg(i),ice_tr_lo(1:2,i)
!  enddo !i

!----------------------------------------------------------------------
! Flux corrected transport algorithm for tracer advection
!
! It is based on Loehner et al. (Finite-element flux-corrected
! transport (FEM-FCT) for the Euler and Navier-Stokes equation,
! Int. J. Numer. Meth. Fluids, 7 (1987), 1093--1109) as described by Kuzmin and
! Turek. (kuzmin@math.uni-dortmund.de)
  ! Auxiliary elemental operator (mass matrix- lumped mass matrix)
  icoef=1
  do i=1,3   ! three upper nodes
    ! Cycle over rows =elnodes(i)
    icoef(i,i)=-2
  enddo
 
  do m=1,ntr_ice
!+++++++++++++++++++++++++++++++++++++++++++++++++
    do i=1,nea !augmented
      do j=1,3 
        icefluxes(j,i)=-sum(icoef(:,j)*(ice_gamma_fct*ice_tr(m,elnode(1:3,i))+ &
  &d_ice_tr(m,elnode(1:3,i))))*voltriangle(i)/lump_ice_matrix(elnode(j,i))/12.0
      enddo !j
    enddo !i

    !Debug
!    do j=1,3
!      swild2=icefluxes(j,:)
!      call exchange_e2d(swild2)
!      icefluxes(j,:)=swild2
!    enddo !j
     
    !==========================   
    ! Screening the low-order solution
    !==========================
    ! Screening means comparing low-order solutions with the
    ! solution on the previous time step and using whichever 
    ! is greater/smaller in computations of max/min below
  
    !==========================
    ! Cluster min/max
    !==========================
    do i=1,np 
      tmax(i)=maxval(ice_tr_lo(m,indnd(1:nnp(i),i)))  !nghbr_nod2D(row)%addresses(1:n)))
      tmax(i)=max(tmax(i),ice_tr_lo(m,i))
      tmin(i)=minval(ice_tr_lo(m,indnd(1:nnp(i),i)))
      tmin(i)=min(tmin(i),ice_tr_lo(m,i))
      ! Admissible increments
      tmax(i)=tmax(i)-ice_tr_lo(m,i) !>=0
      tmin(i)=tmin(i)-ice_tr_lo(m,i) !<=0
    enddo !i

    !=========================
    ! Sums of positive/negative fluxes to node
    !=========================
    icepplus=0.
    icepminus=0.
    do i=1,nea !icepplus invalid @ ghost nodes after this loop
      do j=1,3
        nd=elnode(j,i) 
        flux=icefluxes(j,i)
        if(flux>0) then
          icepplus(nd)=icepplus(nd)+flux !>0
        else
          icepminus(nd)=icepminus(nd)+flux !<=0
        endif
      enddo  
    enddo !i

    !========================
    ! The least upper bound for the correction factors
    !========================
    do i=1,np
      flux=icepplus(i)
      if(flux/=0) then
        icepplus(i)=min(1.d0,tmax(i)/flux) !>=0
      else
        icepplus(i)=0.
      endif

      flux=icepminus(i)
      if(flux/=0) then
        icepminus(i)=min(1.d0,tmin(i)/flux) !>=0
      else
        icepminus(i)=0.
      endif
    enddo !i
    call exchange_p2d(icepminus)
    call exchange_p2d(icepplus) 

    !========================	 
    ! Limiting
    !========================	 
    do i=1,nea !augmented
      ae=1.0
      do j=1,3
        nd=elnode(j,i)
        flux=icefluxes(j,i)
        if(flux>=0.) ae=min(ae,icepplus(nd)) !\in [0,1]
        if(flux<0.) ae=min(ae,icepminus(nd)) !\in [0,1]
      enddo !j
      !ae \in [0,1]
      icefluxes(:,i)=ae*icefluxes(:,i) !valid @ ghost
      !if (ae.le.0.0) write (*,*) 'ae is too large', ae 
    enddo !i
  
    !Debug
!    do j=1,3
!      swild2=icefluxes(j,:)
!      call exchange_e2d(swild2)
!      icefluxes(j,:)=swild2
!    enddo !j

    !==========================
    ! Update the solution 
    !==========================
    ice_tr=ice_tr_lo
    do i=1,nea !ice_tr invalid @ ghost nodes
      do j=1,3
        nd=elnode(j,i)
        ice_tr(m,nd)=ice_tr(m,nd)+icefluxes(j,i)
      enddo !j
    enddo !i  

    swild=ice_tr(m,:)
    call exchange_p2d(swild)
    ice_tr(m,:)=swild

    !Check NaN
    do i=1,npa
      if(ice_tr(m,i)/=ice_tr(m,i)) call parallel_abort('ice_fct: NaN')
    enddo !i
!+++++++++++++++++++++++++++++++++++++++++++++++++
  enddo !m: tracers

!----------------------------------------------------------------------
! Impose bounds
  do i=1,npa
    ice_tr(2,i)=min(ice_tr(2,i),1.d0)
    if(ice_tr(2,i)<1.d-9) ice_tr(2,i)=0
    if(ice_tr(1,i)<1.d-9) ice_tr(1,i)=0
  enddo !i

!Debug
!  if(abs(time_stamp-rnday*86400)<0.1) then
!    fdb='iceha_0000'
!    lfdb=len_trim(fdb)
!    write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
!    open(10,file='outputs/'//fdb,status='replace')
!    write(10,*)np,nproc
!    do i=1,np
!      write(10,'(i11,11(1x,e20.12))')iplg(i),xnd(i),ynd(i),ice_tr(:,i)
!    enddo !
!    close(10)
!  endif

! Deallocate
  deallocate(fct_rhs,ice_tr_lo,d_ice_tr)

!  call parallel_finalize
!  stop

  end subroutine ice_fct

!===============================================================================
