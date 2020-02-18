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


! ifort -O2 -mcmodel=medium -assume byterecl -o modify_hot compute_zcor.f90 modify_hot.f90 -I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff

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
  integer, allocatable :: iest(:),elnode(:,:),i34(:),kbp2(:),nDepth(:),nn(:), &
    & iTest(:)
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

!---------------subroutines--------------------

  iret=nf90_close(sid)

!-------------------modify CB's sal/tem------------------------------
  open(71,file='ChesBay_Data/hot.in',status='old')
  read (71,*) nSt, nSt_J

  allocate(temp1(4),temp2(nvrt),nn(2),iTest(nSt+nSt_J),stat=ierr)
  if(ierr/=0) stop 'Allocation failed (6)'
  allocate(x0(nSt+nSt_J),y0(nSt+nSt_J),nDepth(nSt+nSt_J),val(nSt+nSt_J,200),stat=ierr)
  allocate(z0(nSt+nSt_J,200),stName(nSt+nSt_J),stat=ierr)
  if(ierr/=0) stop 'Allocation failed (7)'

  do k=1,nSt+nSt_J
      read (71,*) x0(k), y0(k), stName(k)
  enddo
  close(71)


  do iter=1,2  !1: salinity; 2: temperature
    iTest=0
    !Read CBP casts
    do k=1,nSt+nSt_J
        if (iter==1) then
            open(72,file='ChesBay_Data/'//trim(adjustl(stName(k)))//'.sal.hot',status='old')
        else
            open(72,file='ChesBay_Data/'//trim(adjustl(stName(k)))//'.tem.hot',status='old')
        endif
        read (72,*) nDepth(k)
        print *, "station: ", k, stName(k)
        do m=1,nDepth(k)
            read(72,*) z0(k,m), val(k,m)
            !write(*,*) z0(k,m), val(k,m)
        enddo
    enddo
    close(72)

    !inside ChesBay, interpolate from CBP casts, FY
    do i=1,np
      if(dp(i)<=0) then
          temp2(:)=0; !stop 'Dry'
          cycle
      endif
      ! find lat between two CB*.* stations
      if(iest(i).ne.1) then !outside Ches Bay
        cycle
      else !main stem
          if (yl(i)>=y0(1)) then
              nn(1)=1; nn(2)=1
          elseif (yl(i)<y0(nSt)) then
              nn(1)=nSt; nn(2)=nSt
          else
              do k=1,nSt-1
                  if (yl(i)<y0(k) .and. yl(i)>=y0(k+1)) then
                      nn(1)=k; nn(2)=k+1
                  endif
              enddo
          endif
      endif

      !find z-level
      do k=1,nvrt
          thisS=abs(z(i,k)/z(i,1))
          do j=1,2
              n0=nn(j)
              if (thisS<=z0(n0,1)/z0(n0,nDepth(n0))) then
                  temp1(j)=val(n0,1)
              else
                  do m=1,nDepth(n0)-1
                      S1=z0(n0,m)/z0(n0,nDepth(n0)) !slayer of the lower
                      S2=z0(n0,m+1)/z0(n0,nDepth(n0)) !slayer of the upper
                      if (thisS>=S1 .and. thisS<=S2) then
                          rat1=abs(S2-thisS)/(S2-S1)
                          rat2=abs(S1-thisS)/(S2-S1)
                          temp1(j)=rat1*val(n0,m)+rat2*val(n0,m+1)
                          exit
                      endif
                  enddo
              endif
          enddo !j
          !horizontal
          if (nn(1)==nn(2)) then
              rat1=.5; rat2=.5
          else
              rat1=abs(y0(nn(2))-yl(i))/abs(y0(nn(2))-y0(nn(1)))
              rat2=abs(y0(nn(1))-yl(i))/abs(y0(nn(2))-y0(nn(1)))
          endif
          temp2(k)=rat1*temp1(1)+rat2*temp1(2)
          if (iter==1) then
              saltout(k,i)=rat1*temp1(1)+rat2*temp1(2)
          else
              tempout(k,i)=rat1*temp1(1)+rat2*temp1(2)
          endif
      enddo !k

      if (iTest(nn(1))==0) then
          print *, nn(1), nn(2)
          print *,  temp2(:)
          iTest(nn(1))=1
      endif

    enddo !i=1,np
  enddo ! iter sal tem

  !!DB
  iTest=0
  open(27,file='DB_elev_ic.gr3',status='old')
  open(28,file='DB_surf_S_ic.gr3',status='old')
  read(27,*); read(28,*)
  read(27,*); read(28,*)
  do i=1,np
    read(27,*) j,xtmp,ytmp,tmp1
    read(28,*) j,xtmp,ytmp,tmp2
    if (iest(i)==2) then !in DB
      eout(i)=tmp1
      saltout(1,i)=tmp2
      tempout(1,i)=21.0 !set 21 oC in DB
      do k=2,nvrt
        tmp=abs( (z(i,k)-z(i,1))/(z(i,nvrt)-z(i,1))  )
        saltout(k,i)=max(0.0,-2.5*tmp+saltout(1,i))
        tempout(k,i)=max(0.0,2.0*tmp+tempout(1,i))
        !if (iTest(1)==0.and.dp(i)>15) then
        !  print *,'k,tmp,z: ',k,tmp,z(i,k),z(i,1),z(i,nvrt)
        !endif
      enddo
      if (iTest(1)==0.and.dp(i)>15) then
        iTest(1)=1
        print *,'eout in DB at: ',i,xl(i),yl(i),eout(i)
        print *,'saltout in DB: ',i,saltout(:,i)
        print *,'tempout in DB: ',i,tempout(:,i)
        !pause
      endif
    else !outside DB
      if (iTest(2)==0.and.dp(i)>150) then
        iTest(2)=1
        print *,'eout outside DB at: ',i,xl(i),yl(i),eout(i)
        print *,'saltout outside DB: ',i,saltout(:,i)
        print *,'tempout outside DB: ',i,tempout(:,i)
      endif
    endif
  enddo !i
  close(27); close(28)

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
