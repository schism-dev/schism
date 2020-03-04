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
! hostart.nc.1


! ifort -O2 -mcmodel=medium -assume byterecl -o replace_hot_in_region compute_zcor.f90 replace_hot_in_region.f90 -I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff

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
  integer, allocatable :: iest(:),elnode(:,:),i34(:),kbp2(:),nDepth(:),nn(:)
    
  character(Len = 1000),allocatable :: stName(:)
  
! -----------------Read grids------------------------
  open(17,file='estuary.gr3',status='old')
  open(16,file='hgrid.ll',status='old')
  open(14,file='hgrid.gr3',status='old') !only need depth info and connectivity
  open(19,file='vgrid.in',status='old')
  open(11,file='fort.11',status='replace')
  read(14,*)
  read(14,*)ne,np

  allocate(iest(np),xl(np),yl(np),dp(np),stat=ierr)
  allocate(elnode(4,ne),i34(ne),stat=ierr)
  if(ierr/=0) stop 'Allocation failed (1)'

  read(16,*)
  read(16,*)
  read(17,*)
  read(17,*)
  do i=1,np
    read(14,*)j,xtmp,ytmp,dp(i)
    read(16,*)j,xl(i),yl(i) 
    read(17,*)j,xtmp,ytmp,tmp
    iest(i)=nint(tmp)
    if(iest(i)<0.and.iest(i)>2) then
      write(11,*)'Estuary flag wrong:',i,iest(i)
      stop
    endif
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


! -----------------Open old hotstart.nc (#0, to be replaced)------------------------
  status = nf90_open('hotstart.nc.0', nf90_nowrite, sid)
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


  allocate(tempout(nvrt,np),saltout(nvrt,np),stat=ierr)
  if(ierr/=0) stop 'Allocation failed (5)'
  tempout=tr_nd(1,:,:); saltout=tr_nd(2,:,:)

  iret=nf90_close(sid)

! -----------------Open another hotstart.nc (#1: to replace #0 with)------------------------
  status = nf90_open('hotstart.nc.1', nf90_nowrite, sid)
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

  allocate(eout1(np),stat=ier)
  if(ierr/=0) stop 'Allocation failed (3)'
  status = nf90_get_var(sid, evid, eout1); call check(status)
  print*, 'Done reading eout'
  print*, 'eout(1): ',eout1(1)
  print*, 'eout(np): ',eout1(np)
  print*, 'min/max eout: ',minval(eout1),maxval(eout1)


  status = nf90_inq_varid(sid, "tr_nd", svid); call check(status)
  print*, 'Done reading tr_nd variable ID'
  status = nf90_Inquire_Variable(sid, svid, dimids = dids); call check(status)
  print*, dids(1),dids(2),dids(3),dids(4)

  allocate(tr_nd1(2,nvrt,np),stat=ierr)
  if(ierr/=0) stop 'Allocation failed (4)'
  status = nf90_get_var(sid, svid, tr_nd1); call check(status)
  print*, 'Done reading tr_nd1'
  print*, 'tr_nd1(1,1,1): ',tr_nd1(1,1,1)
  print*, 'tr_nd1(2,nvrt,np): ',tr_nd1(2,nvrt,np)
  print*, 'min/max tem: ',minval(tr_nd1(1,:,:)),maxval(tr_nd1(1,:,:))
  print*, 'min/max sal: ',minval(tr_nd1(2,:,:)),maxval(tr_nd1(2,:,:))

  iret=nf90_close(sid)


!--------------diagnostic outputs---------------------
  open(27,file='temp_surf.0.gr3',status='replace')
  open(28,file='salt_bot.0.gr3',status='replace')

  write(27,*)
  write(28,*)
  write(27,*) ne,np
  write(28,*) ne,np
  do i=1,np
    write(27,*)i,xl(i),yl(i),tempout(nvrt,i)
    write(28,*)j,xl(i),yl(i),saltout(1,i)
  enddo !i
  do i=1,ne
    write(27,*)j,i34(i),(elnode(l,i),l=1,i34(i))
    write(28,*)j,i34(i),(elnode(l,i),l=1,i34(i))
  enddo !i

  close(27)
  close(28)


!-------------------modify sal/tem------------------------------
  !!DB
  iTest=0
  do i=1,np
    if (iest(i)==2) then !in DB
      !eout(i)=eout1(i)
      saltout(:,i)=tr_nd(2,:,i)
      tempout(:,i)=tr_nd(1,:,i)
      if (iTest==0.and.dp(i)>5) then
        iTest=1
        !print *,'eout in DB at: ',i,xl(i),yl(i),eout(i)
        print *,'saltout in DB: ',i,saltout(:,i)
        print *,'tempout in DB: ',i,tempout(:,i)
        !pause
      endif
    endif
  enddo !i

  !node to element 
  allocate(tsel(2,nvrt,ne),stat=ierr)
  if(ierr/=0) stop 'Allocation failed (8)'
  do i=1,ne
    do k=2,nvrt
      tsel(1,k,i)=(sum(tempout(k,elnode(1:i34(i),i)))+sum(tempout(k-1,elnode(1:i34(i),i))))/2/i34(i) 
      tsel(2,k,i)=(sum(saltout(k,elnode(1:i34(i),i)))+sum(saltout(k-1,elnode(1:i34(i),i))))/2/i34(i)
    enddo !k
    tsel(1,1,i)=tsel(1,2,i) !mainly for hotstart format
    tsel(2,1,i)=tsel(2,2,i)
  enddo !i


!------------write new values for eta2, tr_nd, tr_nd0, tr_el-----------------
  status = nf90_open('hotstart.nc',nf90_write, sid); call check(status)

  status = nf90_inq_varid(sid, "eta2", evid); call check(status)
  iret=nf90_put_var(sid,evid,dble(eout(1:np)),(/1/),(/np/))

  status = nf90_inq_varid(sid, "tr_nd", svid); call check(status)
  iret=nf90_put_var(sid,svid,dble(tempout(1:nvrt,1:np)),(/1,1,1/),(/1,nvrt,np/))
  iret=nf90_put_var(sid,svid,dble(saltout(1:nvrt,1:np)),(/2,1,1/),(/1,nvrt,np/))
  status = nf90_inq_varid(sid, "tr_nd0", svid); call check(status)
  iret=nf90_put_var(sid,svid,dble(tempout(1:nvrt,1:np)),(/1,1,1/),(/1,nvrt,np/))
  iret=nf90_put_var(sid,svid,dble(saltout(1:nvrt,1:np)),(/2,1,1/),(/1,nvrt,np/))

  status = nf90_inq_varid(sid, "tr_el", svid); call check(status)
  iret=nf90_put_var(sid,svid,dble(tsel(1:2,1:nvrt,1:ne)))

  iret=nf90_close(sid)

!--------------diagnostic outputs---------------------
  open(27,file='temp_surf.gr3',status='replace')
  open(28,file='salt_bot.gr3',status='replace')

  write(27,*)
  write(28,*)
  write(27,*) ne,np
  write(28,*) ne,np
  do i=1,np
    write(27,*)i,xl(i),yl(i),tempout(nvrt,i)
    write(28,*)j,xl(i),yl(i),saltout(1,i)
  enddo !i
  do i=1,ne
    write(27,*)j,i34(i),(elnode(l,i),l=1,i34(i))
    write(28,*)j,i34(i),(elnode(l,i),l=1,i34(i))
  enddo !i

  close(27)
  close(28)

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
