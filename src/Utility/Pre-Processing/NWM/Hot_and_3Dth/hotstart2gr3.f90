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

! Modify hotstart.nc based on available data within certain regions

! Inputs:
! hotstart.nc

! Outputs:
! salt_bot.gr3, salt_surf.gr3
! temp_surf.gr3, temp_surf.gr3
! elev.gr3


! ifort -O2 -mcmodel=medium -assume byterecl -o hotstart2gr3 compute_zcor.f90 hotstart2gr3.f90 -I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff

  program gen_hot
  use netcdf
  use compute_zcor

  integer :: status ! netcdf local status variable
  integer :: sid, ncids(100) ! Netcdf file IDs
  integer :: latvid, lonvid, zmvid, hvid ! positional variables
  integer :: xdid, xvid ! longitude index
  integer :: ydid, yvid ! latitude index
  integer :: ldid, lvid ! vertical level, 1 is top
  integer :: svid, tvid ! salt & temp variable IDs
  integer :: uvid, vvid ! vel variable IDs
  integer :: evid ! SSH variable IDs
  integer, dimension(nf90_max_var_dims) :: dids

  real(4), allocatable :: eout(:),eout1(:),xl(:),yl(:),dp(:), &
    & sigma_lcl(:,:),z(:,:),ztot(:),sigma(:),x0(:),y0(:),val(:,:),z0(:,:),temp1(:), &
    & temp2(:),saltout(:,:),tempout(:,:),tr_nd(:,:,:),tr_nd1(:,:,:),tsel(:,:,:)
  integer, allocatable :: elnode(:,:),i34(:),kbp2(:),nDepth(:),nn(:)
    
  character(Len = 1000),allocatable :: stName(:)
  
! -----------------Read grids------------------------
  open(16,file='hgrid.ll',status='old')
  open(14,file='hgrid.gr3',status='old') !only need depth info and connectivity
  open(19,file='vgrid.in',status='old')
  open(11,file='fort.11',status='replace')
  read(14,*)
  read(14,*)ne,np

  allocate(xl(np),yl(np),dp(np),stat=ierr)
  allocate(elnode(4,ne),i34(ne),stat=ierr)
  if(ierr/=0) stop 'Allocation failed (1)'

  write(*,*) "Reading hgrid"
  read(16,*)
  read(16,*)
  do i=1,np
    read(14,*)j,xtmp,ytmp,dp(i)
    read(16,*)j,xl(i),yl(i) 
  enddo !i
  do i=1,ne
    read(14,*)j,i34(i),(elnode(l,i),l=1,i34(i))
  enddo !i

  close(14)
  close(16)

  print*, 'xl(1): ',xl(1)
  print*, 'xl(np): ',xl(np)
  print*, 'min/max xl: ',minval(xl),maxval(xl)
  print*, 'dp(1): ',dp(1)
  print*, 'dp(np): ',dp(np)
  print*, 'min/max dp: ',minval(dp),maxval(dp)


  write(*,*) "Reading vgrid.in"

!     V-grid
  read(19,*)ivcor
  read(19,*)nvrt
  rewind(19)
  allocate(sigma_lcl(nvrt,np),z(np,nvrt),kbp2(np),stat=ierr)
  allocate(ztot(0:nvrt),sigma(nvrt),stat=ierr)
  if(ierr/=0) stop 'Allocation failed (2)'
  call get_vgrid_single('vgrid.in',np,nvrt,ivcor,kz,h_s,h_c,theta_b,theta_f,ztot,sigma,sigma_lcl,kbp2)

!     Compute z-coord.
  do i=1,np
    if(ivcor==2) then
      call zcor_SZ_single(max(0.11,dp(i)),0.,0.1,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot,sigma,z(i,:),idry,kbp2(i))
    else if(ivcor==1) then
      z(i,kbp2(i):nvrt)=max(0.11,dp(i))*sigma_lcl(kbp2(i):nvrt,i)
    else
      write(11,*)'Unknown ivcor:',ivcor
      stop
    endif

    !Extend below bottom for interpolation later
    z(i,1:kbp2(i)-1)=z(i,kbp2(i))
  enddo !i

  print*, 'z(1,1): ',z(1,1)
  print*, 'z(np,1): ',z(np,1)
  print*, 'min/max z: ',minval(z),maxval(z)


! -----------------Open old hotstart.nc ------------------------
  status = nf90_open('hotstart.nc', nf90_nowrite, sid)
  call check(status)

  status = nf90_inq_varid(sid, "eta2", evid); call check(status)
  print*, 'Done reading eout variable ID'
  status = nf90_Inquire_Variable(sid, evid, dimids = dids); call check(status)
  print*, dids(1),dids(2),dids(3),dids(4)

  status = nf90_Inquire_Dimension(sid, dids(1), len = np0); call check(status)
  if (np0.ne.np) then
    print*, 'inconsistent np'
    stop
  endif
  print*, 'np: ',np

  allocate(eout(np),stat=ier)
  if(ierr/=0) stop 'Allocation failed (3)'
  status = nf90_get_var(sid, evid, eout); call check(status)
  print*, 'Done reading eout'
  print*, 'eout(1): ',eout(1)
  print*, 'eout(np): ',eout(np)
  print*, 'min/max eout: ',minval(eout),maxval(eout)


  status = nf90_inq_varid(sid, "tr_nd", svid); call check(status)
  print*, 'Done reading tr_nd variable ID'
  status = nf90_Inquire_Variable(sid, svid, dimids = dids); call check(status)
  print*, dids(1),dids(2),dids(3),dids(4)

  allocate(tr_nd(2,nvrt,np),stat=ierr)
  if(ierr/=0) stop 'Allocation failed (4)'
  status = nf90_get_var(sid, svid, tr_nd); call check(status)
  print*, 'Done reading tr_nd'
  print*, 'tr_nd(1,1,1): ',tr_nd(1,1,1)
  print*, 'tr_nd(2,nvrt,np): ',tr_nd(2,nvrt,np)
  print*, 'min/max tem: ',minval(tr_nd(1,:,:)),maxval(tr_nd(1,:,:))
  print*, 'min/max sal: ',minval(tr_nd(2,:,:)),maxval(tr_nd(2,:,:))


  iret=nf90_close(sid)


!--------------diagnostic outputs---------------------
  open(27,file='temp_surf.gr3',status='replace')
  open(28,file='temp_bot.gr3',status='replace')
  open(37,file='salt_surf.gr3',status='replace')
  open(38,file='salt_bot.gr3',status='replace')
  open(29,file='elev.gr3',status='replace')

  write(27,*); write(28,*); write(37,*); write(38,*); write(29,*)
  write(27,*) ne,np; write(28,*) ne,np; write(37,*) ne,np; write(38,*) ne,np; write(29,*) ne,np
  do i=1,np
    write(27,*)i,xl(i),yl(i),tr_nd(1,nvrt,i)
    write(28,*)i,xl(i),yl(i),tr_nd(1,1,i)
    write(37,*)j,xl(i),yl(i),tr_nd(2,nvrt,i)
    write(38,*)j,xl(i),yl(i),tr_nd(2,1,i)
    write(29,*)j,xl(i),yl(i),eout(i)
  enddo !i
  do i=1,ne
    write(27,*)j,i34(i),(elnode(l,i),l=1,i34(i))
    write(28,*)j,i34(i),(elnode(l,i),l=1,i34(i))
    write(37,*)j,i34(i),(elnode(l,i),l=1,i34(i))
    write(38,*)j,i34(i),(elnode(l,i),l=1,i34(i))
    write(29,*)j,i34(i),(elnode(l,i),l=1,i34(i))
  enddo !i

  close(27); close(28); close(37); close(38); close(29)


!---------------subroutines--------------------
  contains
!     Internal subroutine - checks error status after each netcdf, 
!     prints out text message each time an error code is returned. 
  subroutine check(status)
  integer, intent ( in) :: status

  if(status /= nf90_noerr) then 
    print *, trim(nf90_strerror(status))
    print *, 'failed to open nc files'
    stop
  end if
  end subroutine check  


  end program 
