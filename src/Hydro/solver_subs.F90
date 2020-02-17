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
! ELCIRC LINEAR SOLVER ROUTINES
!
! subroutine solve_jcg
! subroutine solve_jcg_qnon (non-hydrostatic)
! subroutine tridag
!
!===============================================================================
!===============================================================================

subroutine solve_jcg(mnei_p1,np1,npa1,itime,moitn,mxitn,rtol,s,x,b,bc,lbc)
!-------------------------------------------------------------------------------
! Jacobi Preconditioned Conjugate Gradient for elevation-like variables.
!-------------------------------------------------------------------------------
!#ifdef USE_MPIMODULE
!  use mpi
!#endif
  use schism_glbl, only : rkind,np,npa,wtimer,iplg,ipgl,nnp,indnd,errmsg
  use schism_msgp
  implicit none
!#ifndef USE_MPIMODULE
  include 'mpif.h'
!#endif
  integer,intent(in) :: mnei_p1,np1,npa1 !for dimensioning only
  integer,intent(in) :: itime !for outputting/debug only
  integer,intent(in) :: moitn,mxitn    !output interval and max iterations
  real(rkind),intent(in) :: rtol       !relative tolerance
  real(rkind),intent(in) :: s(0:mnei_p1,np1)  !sparse matrix
  real(rkind),intent(inout) :: x(npa1)  !eta2 -- with initial guess
  real(rkind),intent(in) :: b(np1)      !qel
  real(rkind),intent(in) :: bc(npa1)    !b.c. (elbc)
  logical,intent(in) :: lbc(npa1)       !b.c. flag

  integer :: itn,ip,jp,j,omp_get_num_threads,omp_get_thread_num
  real(rkind) :: rdotrl,rdotr,rdotzl,rdotz,old_rdotz,beta,alphal,alpha,cwtmp
  real(rkind) :: rtol2,rdotr0
  real(rkind) :: z(np1),r(np1),p(npa1),sp(np1)
  integer :: inz(mnei_p1,np1),nnz(np1),icolor_intf_node(np1)
  real(rkind) :: snz(0:mnei_p1,np1),bb(np1)
!-------------------------------------------------------------------------------

  itn=0
!$OMP parallel default(shared) private(ip,j,jp)

! Load local storage for non-zeros of sparse matrix and include essential BCs
! Since the calcualtion is done for resident nodes only, entire row of
! the global matrix will be contained in the local row of the matrix.
!$OMP workshare
  nnz=0 !# of non-zero entries for each node
  !1: not interface node btw MPI processes; 0: interface node
  icolor_intf_node=0
!$OMP end workshare

!$OMP do
  do ip=1,np
    if(lbc(ip)) then !nnz(ip)=0: no non-zero entries besides diagonal
      if(bc(ip)<-9998.d0) call parallel_abort('JCG: wrong b.c.')
      snz(0,ip)=1 !modified sparsem
      x(ip)=bc(ip)
      bb(ip)=bc(ip) !modified rrhs (qel)
    else !not essential b.c.
      snz(0,ip)=s(0,ip)
      bb(ip)=b(ip)
      do j=1,nnp(ip)
        jp=indnd(j,ip)
        if(lbc(jp)) then
          if(bc(jp)<-9998.d0) call parallel_abort('JCG: wrong b.c. (2)')
          bb(ip)=bb(ip)-s(j,ip)*bc(jp)
        else
          nnz(ip)=nnz(ip)+1
          inz(nnz(ip),ip)=jp !local node #; can be ghost
          snz(nnz(ip),ip)=s(j,ip)
        endif
      enddo
    endif
  enddo !i=1,np
!$OMP end do

! Initialization -- assume x is initial guess
!#ifdef INCLUDE_TIMING
!  cwtmp=mpi_wtime()
!#endif
!  call exchange_p2d(x) !update x ghost values
!#ifdef INCLUDE_TIMING
!  wtimer(7,2)=wtimer(7,2)+mpi_wtime()-cwtmp
!#endif

! Residual
!$OMP single
  rdotrl=0
!$OMP end single
!$OMP do 
  do ip=1,np
    if(snz(0,ip)==0.d0) call parallel_abort('JCG: zero diagonal')
    sp(ip)=snz(0,ip)*x(ip)+sum(snz(1:nnz(ip),ip)*x(inz(1:nnz(ip),ip)))
!    do j=1,nnz(ip) 
!      sp(ip)=sp(ip)+snz(j,ip)*x(inz(j,ip))
!    enddo
    r(ip)=bb(ip)-sp(ip) !residual
    if(associated(ipgl(iplg(ip))%next)) then !interface node
      if(ipgl(iplg(ip))%next%rank<myrank) cycle !already in the sum so skip
    endif
    icolor_intf_node(ip)=1
!    rdotrl=rdotrl+r(ip)*r(ip)
  enddo !ip
!$OMP end do

!$OMP workshare
  rdotrl=sum(r(:)*r(:)*real(icolor_intf_node(:),rkind))
!$OMP end workshare

!$OMP master
#ifdef INCLUDE_TIMING
  cwtmp=mpi_wtime()
#endif
  call mpi_allreduce(rdotrl,rdotr,1,rtype,MPI_SUM,comm,ierr)
#ifdef INCLUDE_TIMING
  wtimer(7,2)=wtimer(7,2)+mpi_wtime()-cwtmp
#endif

! Convergence based on square of 2norm
  rtol2=rtol*rtol
  rdotr0=rdotr !initial 2norm of residual
  if(rdotr0<0.d0) then
    write(errmsg,*)'JCG: 0 initial error:',rdotr0
    call parallel_abort(errmsg)
  endif
!$OMP end master
!$OMP barrier

!  if(rdotr0==0) return
!  Illegal to branch in/out of parallel region
  if(rdotr0/=0.d0) then
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!$OMP master
  if(myrank==0) then
    write(33,'(//a,i8)') '********CG Solve at timestep ',itime
    write(33,'(a,i6,2e14.6)') &
    'Itn, 2Norm, Rnorm: ',itn,sqrt(rdotr),sqrt(rdotr/rdotr0)
  endif
!$OMP end master

! CG Iteration
  do
    !Need to be checked by all threads to ensure synchronization
    if(rdotr<=rtol2*rdotr0) then
!$OMP master
      if(myrank==0) write(33,*)'JCG converged in ',itn,' iterations'
!$OMP end master
      exit
    endif
    if(itn>=mxitn) then
!$OMP master
      if(myrank==0) write(33,*)'JCG did not converge in ',mxitn,' iterations' 
!$OMP end master
      exit
    endif

!$OMP do
    do ip=1,np 
      !snz(0,ip)/=0 checked
      z(ip)=r(ip)/snz(0,ip) !jacobi
    enddo
!$OMP end do
!#ifdef INCLUDE_TIMING
!    cwtmp=mpi_wtime()
!#endif
!    call exchange_p2d(z) !jacobi
!#ifdef INCLUDE_TIMING
!    wtimer(7,2)=wtimer(7,2)+mpi_wtime()-cwtmp
!#endif

!$OMP single
    rdotzl=0.d0
!$OMP end single

!!!$OMP do reduction(+: rdotzl)
!    do ip=1,np
!      if(associated(ipgl(iplg(ip))%next)) then !interface node
!        if(ipgl(iplg(ip))%next%rank<myrank) cycle !already in the sum so skip
!      endif
!      rdotzl=rdotzl+r(ip)*z(ip)
!    enddo
!!!$OMP end do

!$OMP workshare
    rdotzl=sum(r(:)*z(:)*real(icolor_intf_node(:),rkind))
!$OMP end workshare

!$OMP master
#ifdef INCLUDE_TIMING
    cwtmp=mpi_wtime()
#endif
    call mpi_allreduce(rdotzl,rdotz,1,rtype,MPI_SUM,comm,ierr)
#ifdef INCLUDE_TIMING
    wtimer(7,2)=wtimer(7,2)+mpi_wtime()-cwtmp
#endif
    itn=itn+1
!$OMP end master
!$OMP barrier

    if(itn==1) then 
!$OMP workshare
      p(1:np)=z(1:np)
!$OMP end workshare
    else
!$OMP single
      if(old_rdotz==0.d0) call parallel_abort('JCG: old_rdotz=0')
      beta=rdotz/old_rdotz 
!$OMP end single

!$OMP workshare
      p(1:np)=z(1:np)+beta*p(1:np)
!$OMP end workshare
    endif

!$OMP master
#ifdef INCLUDE_TIMING
    cwtmp=mpi_wtime()
#endif
    call exchange_p2d(p) 
#ifdef INCLUDE_TIMING
    wtimer(7,2)=wtimer(7,2)+mpi_wtime()-cwtmp
#endif

    alphal=0.d0 !temporarily used for the denominator
!$OMP end master
!$OMP barrier

!$OMP do 
    do ip=1,np
      sp(ip)=snz(0,ip)*p(ip)+sum(snz(1:nnz(ip),ip)*p(inz(1:nnz(ip),ip)))
!      do j=1,nnz(ip) 
!        sp(ip)=sp(ip)+snz(j,ip)*p(inz(j,ip)) 
!      enddo
!      if(associated(ipgl(iplg(ip))%next)) then !interface node
!        if(ipgl(iplg(ip))%next%rank<myrank) cycle !already in the sum so skip
!      endif
!      alphal=alphal+p(ip)*sp(ip)
    enddo !ip
!$OMP end do

!$OMP workshare
    alphal=sum(p(1:np)*sp(:)*real(icolor_intf_node(:),rkind))
!$OMP end workshare

!$OMP master
#ifdef INCLUDE_TIMING
    cwtmp=mpi_wtime()
#endif
!    call exchange_p2d(sp) !update sp ghost values
    call mpi_allreduce(alphal,alpha,1,rtype,MPI_SUM,comm,ierr)
#ifdef INCLUDE_TIMING
    wtimer(7,2)=wtimer(7,2)+mpi_wtime()-cwtmp
#endif

    if(alpha==0.d0) call parallel_abort('JCG: division by zero')
    alpha=rdotz/alpha
    old_rdotz=rdotz 
    rdotrl=0.d0
!$OMP end master
!$OMP barrier

!$OMP workshare
    x=x+alpha*p !augmented update
    r=r-alpha*sp !non-augmented update
!$OMP end workshare

!!$OMP do reduction(+: rdotrl)
!    do ip=1,np 
!      if(associated(ipgl(iplg(ip))%next)) then !interface node
!        if(ipgl(iplg(ip))%next%rank<myrank) cycle !already in the sum so skip
!      endif
!      rdotrl=rdotrl+r(ip)*r(ip)
!    enddo
!!$OMP end do

!$OMP workshare
    rdotrl=sum(r(:)*r(:)*real(icolor_intf_node(:),rkind))
!$OMP end workshare

!$OMP master
#ifdef INCLUDE_TIMING
    cwtmp=mpi_wtime()
#endif
    call mpi_allreduce(rdotrl,rdotr,1,rtype,MPI_SUM,comm,ierr)
#ifdef INCLUDE_TIMING
    wtimer(7,2)=wtimer(7,2)+mpi_wtime()-cwtmp
#endif
    if(mod(itn,moitn)==0.and.myrank==0) write(33,'(a,i6,2e14.6)') &
    'Itn, 2Norm, Rnorm: ',itn,sqrt(rdotr),sqrt(rdotr/rdotr0)
!$OMP end master
!$OMP barrier
  enddo !CG Iteration

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  endif !rdotr0/=0
!$OMP end parallel

end subroutine solve_jcg

!===============================================================================
!===============================================================================

subroutine solve_jcg_qnon(itime,moitn,mxitn,rtol,nvrt1,mnei1,np1,npa1,ihydro2,qmatr,qhat,qir)
!-------------------------------------------------------------------------------
! Jacobi Preconditioned Conjugate Gradient for non-hydrostatic model
!-------------------------------------------------------------------------------
!#ifdef USE_MPIMODULE
!  use mpi
!#endif
  use schism_glbl !, only : rkind,np,npa,wtimer,iplg,ipgl,mnei,nnp,indnd,errmsg
  use schism_msgp
  implicit none
!#ifndef USE_MPIMODULE
  include 'mpif.h'
!#endif
  integer,intent(in) :: itime !for info only
  integer,intent(in) :: moitn,mxitn    !output interval and max iterations
  integer,intent(in) :: nvrt1,mnei1,np1,npa1 !used in dimensioning arrays (=nvrt,mnei,np,npa)
  integer,intent(in) :: ihydro2(npa1) !hydrostatic region flag 
  real(rkind),intent(in) :: rtol       !relative tolerance
  real(rkind),intent(in) :: qmatr(nvrt1,-1:1,0:(mnei1+1),np1)  !sparse matrix
  real(rkind),intent(out) :: qhat(nvrt1,npa1)  !final solution
  real(rkind),intent(in) :: qir(nvrt1,np1)     !RHS
  integer :: itn,ip,jp,j,k,kin,l,nd,ndim
  real(rkind) :: rdotrl,rdotr,rdotzl,rdotz,old_rdotz,beta,alphal,alpha,cwtmp,tmp,threshold_rat
  real(rkind) :: rtol2,rdotr0,dx_min2,dz_max2,dist2,rat_max2,rat_max2_gb
  real(rkind) :: zz(nvrt1,np1),rr(nvrt1,np1),pp(nvrt1,npa1),sp(nvrt1,np1)
  real(rkind) :: alow(nvrt1),bdia(nvrt1),cupp(nvrt1),gam(nvrt1),soln(nvrt1,nvrt1),rrhs(nvrt1,nvrt1)
  real(rkind) :: blockj(nvrt1,nvrt1,np)
  logical :: lhbc(npa1) !horizontal essential b.c. nodes
  logical :: large_rat(np1) !flag indicating if theshold grid aspect ratio is exceeded
!-------------------------------------------------------------------------------

! Horizontal essential b.c. nodes (where qhat=0)
  do ip=1,npa
    lhbc(ip)=idry(ip)==1.or.ihydro2(ip)==1 !.or.isbnd(1,ip)>0
  enddo !ip

! Block-Jacobi pre-conditioner
  rat_max2=-1.d0 !max. horizontal-vertical ratio
  threshold_rat=0.4d0 !threshold for grid aspect ratio (dz/dx)
  large_rat(:)=.false.
  do ip=1,np
    if(lhbc(ip)) cycle
    
    !Compute horizontal-vertical ratio
    dx_min2=1.d25
    do j=1,nnp(ip)
      nd=indnd(j,ip)
      dist2=(xnd(ip)-xnd(nd))**2.d0+(ynd(ip)-ynd(nd))**2.d0
      dx_min2=min(dx_min2,dist2)
    enddo !j
    dz_max2=-1
    do k=kbp(ip),nvrt-1
      dist2=(znl(k+1,ip)-znl(k,ip))**2.d0
      dz_max2=max(dz_max2,dist2)
    enddo !k
    if(dx_min2<=0.d0) call parallel_abort('CG2: dx_min2<=0')
    tmp=dz_max2/dx_min2
    rat_max2=max(rat_max2,tmp)
    if(sqrt(tmp)>=threshold_rat) large_rat(ip)=.true.
    if(large_rat(ip)) cycle

    !Invert matrix only for small-aspect-ratio nodes
    ndim=nvrt-kbp_e(ip)
    do k=kbp_e(ip),nvrt-1 !no F.S.
      kin=k-kbp_e(ip)+1
      if(k/=kbp_e(ip)) alow(kin)=qmatr(k,-1,0,ip)
      bdia(kin)=qmatr(k,0,0,ip)
      cupp(kin)=qmatr(k,1,0,ip)
    enddo !k
    !RHS
    rrhs=0.d0
    do l=1,ndim
      rrhs(l,l)=1.d0
    enddo !l

    call tridag(nvrt,nvrt,ndim,ndim,alow,bdia,cupp,rrhs,soln,gam)
    !indice order of blockj: (row #, column #, node)
    blockj(kbp_e(ip):(nvrt-1),kbp_e(ip):(nvrt-1),ip)=transpose(soln(1:ndim,1:ndim))
    
    do k=kbp_e(ip),nvrt-1
      !Check symmetry
      do l=kbp_e(ip),nvrt-1
        if(abs(blockj(k,l,ip)-blockj(l,k,ip))>1.d-4) then
          write(errmsg,*)'Not symmetric:',iplg(ip),k,l,blockj(k,l,ip)-blockj(l,k,ip)
          call parallel_abort(errmsg)
        endif
      enddo !l
      !write(12,*)1/qmatr(k,0,0,ip),ip,k,blockj(k,kbp_e(ip):(nvrt-1),ip) 
    enddo !k
  enddo !ip=1,np

#ifdef INCLUDE_TIMING
  cwtmp=mpi_wtime()
#endif
  call mpi_allreduce(rat_max2,rat_max2_gb,1,rtype,MPI_MAX,comm,ierr)
  if(myrank==0) then
    write(29,'(//a,i8)') '********CG2 Solve at timestep ',itime
    write(29,*)'done pre-conditioner'
    if(rat_max2_gb>0) write(16,*)'Max. vertical/horizontal ratio and threshold=',sqrt(rat_max2_gb),threshold_rat
  endif
  !if(rat_max2_gb>0.16) call parallel_abort('JCG2: grid aspect ratio exceeded')
#ifdef INCLUDE_TIMING
  wtimer(7,2)=wtimer(7,2)+mpi_wtime()-cwtmp
#endif
  
! Residual
  qhat=0.d0 !initial guess
  rdotrl=0.d0
  rr=0.d0 !for b.c. pts etc
  do ip=1,np
    if(lhbc(ip)) cycle

    do k=kbp_e(ip),nvrt-1 !no F.S.
      if(qmatr(k,0,0,ip)<=0.d0) call parallel_abort('JCG2: zero diagonal')
      rr(k,ip)=qir(k,ip) !residual since qhat=0
      if(associated(ipgl(iplg(ip))%next)) then !interface node
        if(ipgl(iplg(ip))%next%rank<myrank) cycle !already in the sum so skip
      endif
      rdotrl=rdotrl+rr(k,ip)*rr(k,ip)
    enddo !k
  enddo !ip=1,np

#ifdef INCLUDE_TIMING
  cwtmp=mpi_wtime()
#endif
  call mpi_allreduce(rdotrl,rdotr,1,rtype,MPI_SUM,comm,ierr)
#ifdef INCLUDE_TIMING
  wtimer(7,2)=wtimer(7,2)+mpi_wtime()-cwtmp
#endif

! Convergence based on square of 2norm
  rtol2=rtol*rtol
  rdotr0=rdotr !initial 2norm of residual
  if(rdotr0<0.d0) call parallel_abort('JCG2: 0 initial error')
  if(rdotr0==0.d0) return

  itn=0
  if(myrank==0) then
    write(29,'(a,i6,3e14.6)')'Itn, 2Norm, Rnorm, tol2: ',itn,rdotr,rdotr/rdotr0,rtol2
  endif

! CG Iteration
  do
    if(rdotr<=rtol2*rdotr0) then
      if(myrank==0) write(29,*)'JCG2 converged in ',itn,' iterations'
      exit
    endif
    if(itn>=mxitn) then
      if(myrank==0) write(29,*)'JCG2 did not converge in ',mxitn,' iterations'
      exit
    endif

    zz=0.d0 !for b.c. pts
    do ip=1,np 
      if(lhbc(ip)) cycle

      !Block or simple Jacobian
      do k=kbp_e(ip),nvrt-1
        if(large_rat(ip)) then !simple
          zz(k,ip)=rr(k,ip)/qmatr(k,0,0,ip) !diagonal checked
        else !block
          do l=kbp_e(ip),nvrt-1
            zz(k,ip)=zz(k,ip)+blockj(k,l,ip)*rr(l,ip)
          enddo !l
        endif !large_rat
      enddo !k
    enddo !ip

    rdotzl=0.d0
    do ip=1,np
      if(lhbc(ip)) cycle
      if(associated(ipgl(iplg(ip))%next)) then !interface node
        if(ipgl(iplg(ip))%next%rank<myrank) cycle !already in the sum so skip
      endif
      do k=kbp_e(ip),nvrt-1
        rdotzl=rdotzl+rr(k,ip)*zz(k,ip)
      enddo !k
    enddo !ip
#ifdef INCLUDE_TIMING
    cwtmp=mpi_wtime()
#endif
    call mpi_allreduce(rdotzl,rdotz,1,rtype,MPI_SUM,comm,ierr)
#ifdef INCLUDE_TIMING
    wtimer(7,2)=wtimer(7,2)+mpi_wtime()-cwtmp
#endif
    itn=itn+1
    if(itn==1) then 
      pp(:,1:np)=zz(:,1:np)
    else
      if(old_rdotz==0.d0) call parallel_abort('JCG: old_rdotz=0')
      beta=rdotz/old_rdotz 
      pp(:,1:np)=zz(:,1:np)+beta*pp(:,1:np)
    endif

#ifdef INCLUDE_TIMING
    cwtmp=mpi_wtime()
#endif
    call exchange_p3dw(pp) 
#ifdef INCLUDE_TIMING
    wtimer(7,2)=wtimer(7,2)+mpi_wtime()-cwtmp
#endif

    alphal=0.d0 !temporarily used for the denominator
    sp=0.d0 !for b.c. pts and initialize
    do ip=1,np
      if(lhbc(ip)) cycle

      do k=kbp_e(ip),nvrt-1 !no F.S.
        !sp(k,ip)=0
        do l=-1,1
          if(k+l<kbp_e(ip).or.k+l>=nvrt) cycle !repeat F.S. essential b.c. (no harm)
          do j=0,nnp(ip)
            if(j==0) then
              nd=ip
            else
              nd=indnd(j,ip)
            endif
            if(lhbc(nd).or.k+l==nvrt) cycle !essential b.c.
            sp(k,ip)=sp(k,ip)+qmatr(k,l,j,ip)*pp(k+l,nd)
          enddo !j
        enddo !l
        if(associated(ipgl(iplg(ip))%next)) then !interface node
          if(ipgl(iplg(ip))%next%rank<myrank) cycle !already in the sum so skip
        endif
        alphal=alphal+pp(k,ip)*sp(k,ip)
      enddo !k=kbp_e(ip),nvrt
    enddo !ip

#ifdef INCLUDE_TIMING
    cwtmp=mpi_wtime()
#endif
    call mpi_allreduce(alphal,alpha,1,rtype,MPI_SUM,comm,ierr)
#ifdef INCLUDE_TIMING
    wtimer(7,2)=wtimer(7,2)+mpi_wtime()-cwtmp
#endif

    if(alpha==0.d0) call parallel_abort('JCG2: division by zero')
    alpha=rdotz/alpha

    qhat=qhat+alpha*pp !augmented update
    rr=rr-alpha*sp !non-augmented update
    old_rdotz=rdotz 
    rdotrl=0.d0
    do ip=1,np 
      !if(idry(ip)==1.or.isbnd(1,ip)>0) cycle
      if(lhbc(ip)) cycle
      if(associated(ipgl(iplg(ip))%next)) then !interface node
        if(ipgl(iplg(ip))%next%rank<myrank) cycle !already in the sum so skip
      endif
      do k=kbp_e(ip),nvrt-1
        rdotrl=rdotrl+rr(k,ip)*rr(k,ip)
      enddo !k
    enddo
#ifdef INCLUDE_TIMING
    cwtmp=mpi_wtime()
#endif
    call mpi_allreduce(rdotrl,rdotr,1,rtype,MPI_SUM,comm,ierr)
#ifdef INCLUDE_TIMING
    wtimer(7,2)=wtimer(7,2)+mpi_wtime()-cwtmp
#endif
    if(mod(itn,moitn)==0.and.myrank==0) write(29,'(a,i6,3e14.6)') &
    'Itn, 2Norm, Rnorm, tol2: ',itn,rdotr,rdotr/rdotr0,rtol2
  enddo !CG Iteration

end subroutine solve_jcg_qnon

!===============================================================================
!===============================================================================
subroutine tridag(nmax,nvec,n,nc,a,b,c,r,u,gam)
!-------------------------------------------------------------------------------
! This program solves a tridiagonal system. It was adapted from "Numerical 
! Recipes in FORTRAN (pp.43 ).
!
! a,b,c,r: input vectors and are not modified by this program.
! b is the main diagonal, a below and c above. a(1) and c(n) are not used.
! r is the r.h.s.
! (nvec,nmax) is the dimension of r() _and_ u() in the driving routine.
! n: actual rank of the system.
! nc: input; actual # of columns of rhs (<= nvec). For efficiency, this is # of
! rows in u() and r();
! u: output with nc columns (transposed to rows)
! gam: a working array.
!-------------------------------------------------------------------------------
  use schism_glbl, only : rkind,errmsg
  use schism_msgp, only : parallel_abort
  implicit none

  integer, intent(in) :: nmax,nvec,n,nc
  real(rkind), dimension(nmax), intent(in) :: a,b,c
  real(rkind), dimension(nvec,nmax), intent(in) :: r
  real(rkind), dimension(nmax), intent(out) :: gam
  real(rkind), dimension(nvec,nmax), intent(out) :: u

  integer :: i,j,ifl
  real(rkind) :: bet

  if(n<1) call parallel_abort('TRIDAG: n must be >= 1')
  if(nc>nvec) call parallel_abort('TRIDAG: Increase # of columns')
  if(b(1)==0.d0) call parallel_abort('TRIDAG:  b(1)=0')

  bet=b(1)
  u(1:nc,1)=r(1:nc,1)/bet

  ifl=0 !flag for abort (for better vectorization)
  do j=2,n
    gam(j)=c(j-1)/bet
    bet=b(j)-a(j)*gam(j)
    if(bet==0.d0) ifl=1
    u(1:nc,j)=(r(1:nc,j)-a(j)*u(1:nc,j-1))/bet
  enddo !j

  if(ifl==1) call parallel_abort('TRIDAG: failed')

! Backsubstitution
  do j=n-1,1,-1
    u(1:nc,j)=u(1:nc,j)-gam(j+1)*u(1:nc,j+1)
  enddo
  
end subroutine tridag

