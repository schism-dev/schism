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
! SCHISM GRID SUBROUTINES used by QSim
!
! subroutine aquire_vgrid
!
!===============================================================================
!===============================================================================

subroutine aquire_vgrid
!-------------------------------------------------------------------------------
! Aquire vertical grid data from vgrid.in
!-------------------------------------------------------------------------------
  use schism_glbl
  use schism_msgp
  implicit none
  integer :: i,j,k,l,jki,stat,kin,m
  real(rkind) :: buf1(100),hmod2,zz

  !ivcor: types of vertical coord.; surface must all be nvrt (for sflux routines)
  if(myrank==0) then
    open(19,file=in_dir(1:len_in_dir)//'vgrid.in',status='old',iostat=stat)
    if(stat/=0) call parallel_abort('AQUIRE_VGIRD: open(19) failure')
    read(19,*)ivcor
  endif !myrank
  call mpi_bcast(ivcor,1,itype,0,comm,stat)

  if(ivcor==2) then !SZ coordinates
    if(myrank==0) then
      read(19,*) nvrt,kz,h_s !kz>=1
      if(nvrt<2) call parallel_abort('nvrt<2')
      if(kz<1) then !.or.kz>nvrt-2) then
        write(errmsg,*)'Wrong kz:',kz
        call parallel_abort(errmsg)
      endif
      if(h_s<6.d0) then
        write(errmsg,*)'h_s needs to be larger:',h_s
        call parallel_abort(errmsg)
      endif
    endif !myrank
    call mpi_bcast(nvrt,1,itype,0,comm,stat)
    call mpi_bcast(kz,1,itype,0,comm,stat)
    call mpi_bcast(h_s,1,rtype,0,comm,stat)

    ! Allocate vertical layers arrays
    allocate(ztot(nvrt),sigma(nvrt),cs(nvrt),dcs(nvrt),stat=stat)
    if(stat/=0) call parallel_abort('AQUIRE_VGIRD: ztot allocation failure')
    nsig=nvrt-kz+1 !# of S levels (including "bottom" & f.s.)

    if(myrank==0) then
      ! # of z-levels excluding "bottom" at h_s
      read(19,*) !for adding comment "Z levels"
      do k=1,kz-1
        read(19,*)j,ztot(k)
        if(ztot(k)>=-h_s) then
          write(errmsg,*)'Illegal Z level:',k
          call parallel_abort(errmsg)
        endif
        if(k>1) then; if(ztot(k)<=ztot(k-1)) then
          write(errmsg,*)'z-level inverted:',k
          call parallel_abort(errmsg)
        endif; endif
      enddo !k
      read(19,*) !level kz       
      ! In case kz=1, there is only 1 ztot(1)=-h_s
      ztot(kz)=-h_s

      read(19,*) !for adding comment "S levels"
      read(19,*)h_c,theta_b,theta_f
      if(h_c<5._rkind.or.h_c>=h_s) then !large h_c to avoid 2nd type abnormaty
        write(errmsg,*)'h_c needs to be larger:',h_c
        call parallel_abort(errmsg)
      endif
      if(theta_b<0._rkind.or.theta_b>1._rkind) then
        write(errmsg,*)'Wrong theta_b:',theta_b
        call parallel_abort(errmsg)
      endif
      if(theta_f<=0._rkind) then
        write(errmsg,*)'Wrong theta_f:',theta_f
        call parallel_abort(errmsg)
      endif

      sigma(1)=-1._rkind !bottom
      sigma(nsig)=0._rkind !surface
      read(19,*) !level kz
      do k=kz+1,nvrt-1
        kin=k-kz+1
        read(19,*) j,sigma(kin)
        if(sigma(kin)<=sigma(kin-1).or.sigma(kin)>=0._rkind) then
          write(errmsg,*)'Check sigma levels at:',k,sigma(kin),sigma(kin-1)
          call parallel_abort(errmsg)
        endif
      enddo !k
      read(19,*) !level nvrt
      close(19)
    endif !myrank
    call mpi_bcast(ztot,nvrt,rtype,0,comm,stat)
    call mpi_bcast(h_c,1,rtype,0,comm,stat)
    call mpi_bcast(theta_b,1,rtype,0,comm,stat)
    call mpi_bcast(theta_f,1,rtype,0,comm,stat)
    call mpi_bcast(sigma,nvrt,rtype,0,comm,stat)

    !Pre-compute constants
    s_con1=sinh(theta_f)
  else if(ivcor==1) then !localized sigma
    if(myrank==0) read(19,*)nvrt !needs hgrid to read the rest
    call mpi_bcast(nvrt,1,itype,0,comm,stat)
    close(19)
    allocate(ztot(nvrt),sigma(nvrt),stat=stat)
    if(stat/=0) call parallel_abort('AQUIRE_VGIRD: ztot allocation failure (2)')
    !for output only - remove later
    ztot=0._rkind; sigma=0._rkind 
    kz=1; h_s=0._rkind; h_c=0._rkind; theta_b=0._rkind; theta_f=0._rkind
  else
    call parallel_abort('GRID_SUBS: Unknown ivcor')
  endif !ivcor
!  endif !lm2d

  if(ivcor==2) then
    ! Compute C(s) and C'(s)
    do k=1,nsig
      cs(k)=(1._rkind-theta_b)*sinh(theta_f*sigma(k))/sinh(theta_f)+ &
       &theta_b*(tanh(theta_f*(sigma(k)+0.5_rkind))-tanh(theta_f*0.5_rkind))/2._rkind/tanh(theta_f*0.5_rkind)
      dcs(k)=(1._rkind-theta_b)*theta_f*cosh(theta_f*sigma(k))/sinh(theta_f)+ &
       &theta_b*theta_f/2._rkind/tanh(theta_f*0.5_rkind)/cosh(theta_f*(sigma(k)+0.5_rkind))**2._rkind
    enddo !k
  endif !ivcor==2

  ! Output some sample z-coordinates
!  if(myrank==0) then
!    open(10,file='sample_z.out',status='replace')
!    write(10,*)'Sample z coordinates'
!    buf1(1)=h_s; buf1(2)=h_c; buf1(2:11)=(/(10*(i-1),i=2,11)/); buf1(12:28)=(/(200+50*(i-12), i=12,28)/)
!    write(10,*)'h_c= ',h_c,', h_s=',h_s
!    do i=1,28
!      write(10,*)'Depth= ',buf1(i)
!      do k=kz,nvrt
!        kin=k-kz+1
!        hmod2=min(buf1(i),h_s)
!        if(hmod2<=h_c) then
!          zz=sigma(kin)*hmod2
!        else
!          zz=h_c*sigma(kin)+(hmod2-h_c)*cs(kin)
!        endif
!        write(10,*)k,zz
!      enddo !k
!    enddo !i
!    close(10)
!  endif

end subroutine aquire_vgrid
