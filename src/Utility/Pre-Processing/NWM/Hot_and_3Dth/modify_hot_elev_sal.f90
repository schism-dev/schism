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
! hgrid.gr3
! vgrid.in

! Outputs:
! hotstart.nc (overwrite overland elev, raise initial elev to 0.1 m below
! ground)


! ifort -O2 -mcmodel=medium -assume byterecl -o modify_hot_elev_sal compute_zcor.f90 modify_hot_elev_sal.f90 -I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff

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

  real(4), allocatable :: eout(:),xl(:),yl(:),dp(:), &
    & sigma_lcl(:,:),z(:,:),ztot(:),sigma(:),x0(:),y0(:),val(:,:),z0(:,:),temp1(:), &
    & temp2(:),saltout(:,:),tempout(:,:),tr_nd(:,:,:),tsel(:,:,:)
  integer, allocatable :: i_ic_tem(:),i_ic_sal(:),elnode(:,:),i34(:),kbp2(:),nDepth(:),nn(:), &
    & iTest(:)
  character(Len = 1000),allocatable :: stName(:)

  h0=1e-1
  
! -----------------Read grids------------------------
  open(14,file='hgrid.gr3',status='old') 
  read(14,*)
  read(14,*)ne,np

  allocate(xl(np),yl(np),dp(np),stat=ierr)
  allocate(elnode(4,ne),i34(ne),stat=ierr)
  if(ierr/=0) stop 'Allocation failed (1)'

  do i=1,np
    read(14,*)j,xl(i),yl(i),dp(i)
  enddo !i
  do i=1,ne
    read(14,*)j,i34(i),(elnode(l,i),l=1,i34(i))
  enddo !i
  close(14);

  print*, 'dp(1): ',dp(1)
  print*, 'dp(np): ',dp(np)
  print*, 'min/max dp: ',minval(dp),maxval(dp)


  open(19,file='vgrid.in',status='old')
  read(19,*)ivcor
  read(19,*)nvrt
  close(19)
! -----------------Open old hotstart.nc------------------------
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

  allocate(tempout(nvrt,np),saltout(nvrt,np),stat=ierr)
  if(ierr/=0) stop 'Allocation failed (5)'
  tempout=tr_nd(1,:,:); saltout=tr_nd(2,:,:)

  iret=nf90_close(sid)

  ! coastal regions: 
  ! raise elev to just below-ground
  do i=1,np
    if (dp(i)<0.0) then ! above vertical datum
      eout(i)=max(0.d0,-dp(i)-h0)
    endif
  enddo
  !set salinity above 3m NAVD to be 0, linear interpolation between original value
  !and 0 psu if -3 <= dp <= 0
  do i=1,np
    if (dp(i)<-3)then
      saltout(:,i)=0.0
    elseif (dp(i)<0)then
      rat=max(min(1.0, (dp(i)+3.0)/3.0), 0.0)
      saltout(:,i)=saltout(:,i) * rat
    endif
  enddo!i

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
  open(26,file='elev.gr3',status='replace')
  open(27,file='temp_surf.gr3',status='replace')
  open(28,file='salt_bot.gr3',status='replace')

  write(26,*)
  write(27,*)
  write(28,*)
  write(26,*) ne,np
  write(27,*) ne,np
  write(28,*) ne,np
  do i=1,np
    write(26,*)i,xl(i),yl(i),eout(i)
    write(27,*)i,xl(i),yl(i),tempout(nvrt,i)
    write(28,*)j,xl(i),yl(i),saltout(1,i)
  enddo !i
  do i=1,ne
    write(26,*)j,i34(i),(elnode(l,i),l=1,i34(i))
    write(27,*)j,i34(i),(elnode(l,i),l=1,i34(i))
    write(28,*)j,i34(i),(elnode(l,i),l=1,i34(i))
  enddo !i

  close(26)
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
