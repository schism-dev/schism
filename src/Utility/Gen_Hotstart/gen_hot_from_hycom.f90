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

! Generate hotstart.nc only (no *.th.nc) from gridded HYCOM data (nc file); works for global grid.
! Changed algorithm from gen_hot_3Dth_from_hycom: no longer do
! horizontal extension but simply fill invalid pts with T,S values from a valid
! location (0 for SSH, U,V).

! The section on nc read needs to be modified as appropriate- search for
! 'start11' and 'end11'
! Beware the order of vertical levels in the nc file!!!
! Assume elev=0 in interpolation of 3D variables
! The code will do all it can to extrapolate: below bottom/above surface.
! If a pt in hgrid.ll is outside the background nc grid, const. values will be filled (from gen_hot_from_nc.in) 
!
! Dan Yu: added checking vertical layer for each 3D var from HYCOM nc, and some easy remedy for nc file error  
! Wet-dry points are defined by ssh.
! For U,V, if junk values in the middle of water, simply fill with 0 and record in fort.20
! For S,T, if junk values in the middle of water, simply fill with bottom value, and record in fort.20
! If SSH shows wet/dry in time, search nearby points to fill. If none found, fill 0, and record in fort.20
! Consider checking HYCOM nc files with ncview or other tools.

!   Input: 
!     (1) hgrid.gr3;
!     (2) hgrid.ll;
!     (3) vgrid.in (SCHISM R1703 and up);
!     (4) estuary.gr3 (flags for extrapolating S,T, vel.): depth=0: outside; =1: inside
!     (5) gen_hot_from_nc.in: 
!                     1st line: 1: include vel and elev. in hotstart.nc; 0: only T,S
!                     2nd line: T,S values for estuary points defined in estuary.gr3
!                     3rd line: T,S values for pts outside background grid in nc
!                     4th line: time step in .nc in sec
!     (6) HYCOM files: [SSH,TS,UV]_[1,2,..nfiles].nc (beware scaling etc)
!   Output: hotstart.nc
!   Debug outputs: fort.11 (fatal errors); fort.20 (warning); fort.2[1-9], fort.9[5-9], fort.100; backup.out

! ifort -O2 -mcmodel=medium -assume byterecl -CB -o gen_hot_from_hycom.exe ../UtilLib/schism_geometry.f90 \
! ../UtilLib/extract_mod.f90 ../UtilLib/compute_zcor.f90 ../UtilLib/pt_in_poly_test.f90 gen_hot_from_hycom.f90 \
!-I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff

      program gen_hot
      use netcdf
      use schism_geometry_mod
      use compute_zcor
      use pt_in_poly_test, only: signa_single

!      implicit none

      integer, parameter :: nbyte=4
      real, parameter :: small1=1.e-2 !used to check area ratios
  
!     netcdf related variables
      integer :: sid, ncids(100) ! Netcdf file IDs
      integer :: latvid, lonvid, zmvid, hvid ! positional variables
      integer :: xdid, xvid ! longitude index
      integer :: ydid, yvid ! latitude index
      integer :: ldid, lvid ! vertical level, 1 is top
      integer :: svid, tvid ! salt & temp variable IDs
      integer :: uvid, vvid ! vel variable IDs
      integer :: evid ! SSH variable IDs
      integer :: one_dim,ivarid(100),var1d_dims(1),var2d_dims(2),var3d_dims(3), &
     &var4d_dims(4),itime_id(4)
      integer, dimension(nf90_max_var_dims) :: dids
             
!     Local variables for data
      real (kind = 4), allocatable :: xind(:), yind(:), lind(:) 
!     Lat, lon, bathymetry
      real (kind = 4), allocatable :: lat(:,:), lon(:,:), hnc(:)
!     Vertical postion, salinity, and temperature
      real (kind = 4), allocatable :: zm(:,:,:),salt(:,:,:),temp(:,:,:), &
     &uvel(:,:,:),vvel(:,:,:),ssh(:,:),salt0(:,:,:),temp0(:,:,:),uvel0(:,:,:),vvel0(:,:,:),ssh0(:,:)
      integer, allocatable :: kbp(:,:),ihope(:,:)
!     File names for netcdf files
      character(len=1024) :: ncfile1
!     Command line arguments
      character(len=1024) :: s, yr, md
      character(len=4) :: iyear_char
      character(len=1) :: char1,char2
      character(len=20) :: char3,char4
!     external function for number of command line arguments
!      integer :: iargc

      integer :: status ! netcdf local status variable
      integer :: ier ! allocate error return.
      integer :: ixlen, iylen, ilen ! sampled lengths in each coordinate direction
      integer :: ixlen1, iylen1,ixlen2, iylen2 !reduced indices for CORIE grid to speed up interpolation
      integer :: klev(4)
      real :: wild(100),wild2(100,2)

      integer, allocatable :: elnode(:,:),elside(:,:),isdel(:,:),idry_s(:),i34(:),iest(:), &
     &ixy(:,:),ic3(:,:),isidenode(:,:),nond(:),iond(:,:),iob(:),iond2(:)
      real, allocatable :: xl(:),yl(:),dp(:),ztot(:),sigma(:),arco(:,:), &
     &tempout(:,:),saltout(:,:),uout(:,:),vout(:,:),eout(:),su2(:,:),sv2(:,:), &
     &eout_tmp(:),tsel(:,:,:),zeros(:,:),xcj(:),ycj(:),ts_lat(:,:,:)
      allocatable :: z(:,:),sigma_lcl(:,:),kbp2(:),iparen_of_dry(:,:)
      real*8 :: aa1(1)

!     First statement
!     Currently we assume rectangular grid in HYCOM
!     (interp_mode=0). interp_mode=1 uses generic UG search (splitting
!     quads) and is kept for more generic cases
      interp_mode=0

      open(10,file='gen_hot_from_nc.in',status='old')
      read(10,*) iuv !1: include vel and elev. in hotstart.in; 0: only T,S
      read(10,*) tem_es,sal_es !T,S values for estuary points defined in estuary.gr3
      read(10,*) tem_outside,sal_outside !T,S values for pts outside bg grid in nc
!      read(10,*) dtout !time step in .nc [sec]
!      read(10,*) nndays !# of days needed in output
!      read(10,*) nfiles !# of stacks of HYCOM files
      close(10)
      if(tem_es<0.or.sal_es<0.or.tem_outside<0.or.sal_outside<0) &
     &stop 'Invalid T,S constants'

!     Read in hgrid and vgrid
      open(17,file='estuary.gr3',status='old')
      open(16,file='hgrid.ll',status='old')
      open(14,file='hgrid.gr3',status='old') !only need depth info and connectivity
      open(19,file='vgrid.in',status='old')
      open(11,file='fort.11',status='replace')
      read(14,*)
      read(14,*)ne,np
      allocate(xl(np),yl(np),dp(np),i34(ne),elnode(4,ne),iest(np),ixy(np,3),stat=istat)
      if(istat/=0) stop 'Failed to alloc. (0)'
      read(16,*)
      read(16,*)
      read(17,*)
      read(17,*)
      do i=1,np
        read(14,*)j,xtmp,ytmp,dp(i)
        read(16,*)j,xl(i),yl(i) 
        read(17,*)j,xtmp,ytmp,tmp
        iest(i)=nint(tmp)
        if(iest(i)/=0.and.iest(i)/=1) then
          write(11,*)'Estuary flag wrong:',i,iest(i)
          stop
        endif
      enddo !i
      do i=1,ne
        read(14,*)j,i34(i),(elnode(l,i),l=1,i34(i))
      enddo !i
      close(14)
      close(16)
      close(17)

!     V-grid
      read(19,*)ivcor
      read(19,*)nvrt
      rewind(19)
      allocate(sigma_lcl(nvrt,np),z(nvrt,np),kbp2(np),ztot(0:nvrt),sigma(nvrt),arco(4,np), &
     &tempout(nvrt,np),saltout(nvrt,np),uout(nvrt,np),vout(nvrt,np),eout(np), &
     &eout_tmp(np),tsel(2,nvrt,ne),zeros(nvrt,max(ne,np)),stat=istat)
      if(istat/=0) stop 'Failed to alloc. (1)'
      call get_vgrid_single('vgrid.in',np,nvrt,ivcor,kz,h_s,h_c,theta_b,theta_f,ztot,sigma,sigma_lcl,kbp2)

!     Compute z-coord.
      do i=1,np
        if(ivcor==2) then
          call zcor_SZ_single(max(0.11,dp(i)),0.,0.1,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot,sigma,z(:,i),idry,kbp2(i))
        else if(ivcor==1) then
          z(kbp2(i):nvrt,i)=max(0.11,dp(i))*sigma_lcl(kbp2(i):nvrt,i)
        else
          write(11,*)'Unknown ivcor:',ivcor
          stop
        endif

        !Extend below bottom for interpolation later
        z(1:kbp2(i)-1,i)=z(kbp2(i),i)
      enddo !i

!     Compute geometry
      call compute_nside(np,ne,i34,elnode,ns)
      allocate(ic3(4,ne),elside(4,ne),isdel(2,ns),isidenode(2,ns),xcj(ns),ycj(ns), &
     &idry_s(ns),su2(nvrt,ns),sv2(nvrt,ns),stat=istat)
      if(istat/=0) stop 'Failed to alloc. (2)'

      call schism_geometry_single(np,ne,ns,xl,yl,i34,elnode,ic3,elside,isdel,isidenode,xcj,ycj)

      if(ns<ne.or.ns<np) then
        write(11,*)'Weird grid with ns < ne or ns < np',np,ne,ns
        stop
      endif
      print*, 'Done computing geometry...'

      open(12,file='backup.out',status='replace')

!start11
!     Define limits for variables for sanity checks
      tempmin=-10 
      tempmax=50
      saltmin=0
      saltmax=60
      vmag_max=10 !max. |u| or |v|
      ssh_max=5 !max. of |SSH|

!     Assuming elev, u,v, S,T have same dimensions (and time step) and grids do not change over time
!      timeout=0
      do ifile=1,1 !nfiles
!--------------------------------------------------------------------
      write(char3,'(i20)')ifile
      char3=adjustl(char3); len_char3=len_trim(char3)
!     Open nc file 
      status = nf90_open('TS_'//char3(1:len_char3)//'.nc', nf90_nowrite, sid)
      call check(status)
      status = nf90_open('SSH_'//char3(1:len_char3)//'.nc', nf90_nowrite, ncids(1))
      status = nf90_open('UV_'//char3(1:len_char3)//'.nc', nf90_nowrite, ncids(2))

      if(ifile==1) then
!       Get dimnesions from S
        status = nf90_inq_varid(sid, "salinity", svid)
        print*, 'Done reading variable ID'
        print*, sid,svid,ncids(1:2),'TS_'//char3(1:len_char3)//'.nc'

!       Get dimensions from the 1st file
!       Assumed same as for all variables
!       WARNING: indices reversed from ncdump!
        status = nf90_Inquire_Variable(sid, svid, dimids = dids)
        status = nf90_Inquire_Dimension(sid, dids(1), len = ixlen)
        status = nf90_Inquire_Dimension(sid, dids(2), len = iylen)
        status = nf90_Inquire_Dimension(sid, dids(3), len = ilen)
        status = nf90_Inquire_Dimension(sid, dids(4), len = ntime)
        print*, dids(1),dids(2),dids(3),dids(4)
        print*, 'ixlen,iylen,ilen,ntime= ',ixlen,iylen,ilen,ntime

!       allocate arrays
        allocate(xind(ixlen),stat=ier)
        allocate(yind(iylen),stat=ier)
        allocate(lind(ilen),stat=ier)
        allocate(lat(ixlen,iylen))
        allocate(lon(ixlen,iylen))
        allocate(zm(ixlen, iylen, ilen))
!        allocate(hnc(ixlen,iylen))
        allocate(hnc(ilen))
        allocate(kbp(ixlen,iylen))
        allocate(ihope(ixlen,iylen))
        allocate(uvel(ixlen,iylen,ilen),stat=ier)
        allocate(vvel(ixlen,iylen,ilen),stat=ier)
        allocate(salt(ixlen,iylen,ilen),stat=ier)
        allocate(ts_lat(iylen,ilen,2),stat=ier)
        allocate(temp(ixlen,iylen,ilen),stat=ier)
        allocate(ssh(ixlen,iylen),stat=ier)
        uvel=0; vvel=0; ssh=0
  
!       get static info (lat/lon grids etc) 
        status = nf90_inq_varid(ncids(1), "xlon", xvid)
        status = nf90_get_var(ncids(1), xvid, xind)
        status = nf90_inq_varid(ncids(1), "ylat", yvid)
        status = nf90_get_var(ncids(1), yvid, yind)
        !lind may be sigma coord.
!        status = nf90_inq_varid(sid, "sigma", lvid)
!        status = nf90_get_var(sid, lvid, lind)
        status = nf90_inq_varid(sid, "depth", hvid)
        status = nf90_get_var(sid, hvid, hnc)

!       processing static info
        do i=1,ixlen
          lon(i,:)=xind(i)
          if(i<ixlen) then; if(xind(i)>=xind(i+1)) then
            write(11,*)'Lon must be increasing:',i,xind(i),xind(i+1)
            stop
          endif; endif;
        enddo !i
        do j=1,iylen
          lat(:,j)=yind(j)
          if(j<iylen) then; if(yind(j)>=yind(j+1)) then
            write(11,*)'Lat must be increasing:',j,yind(j),yind(j+1)
            stop
          endif; endif;
        enddo !j
!        lon=lon-360 !convert to our long.

!       Compute z-coord. (assuming eta=0)
!       WARNING: In zm(), 1 is bottom; ilen is surface (SCHISM convention)
        do i=1,ixlen
          do j=1,iylen
            do k=1,ilen
              zm(i,j,k)=-hnc(1+ilen-k) !lind(k)*hnc(i,j)
            enddo !k
          enddo !j
        enddo !i

!       Get rest of varid's
        status = nf90_inq_varid(sid, "temperature", tvid)
        status = nf90_inq_varid(ncids(1), "surf_el", evid)
        status = nf90_inq_varid(ncids(2), "water_u", uvid)
        status = nf90_inq_varid(ncids(2), "water_v", vvid)

!       Arrays no longer used after this: hnc,lind
        deallocate(hnc,lind)
      endif !ifile==1

!     Vertical convention follows SCHISM from now on; i.e., 1 is at bottom
!     Compute bottom indices
!      if(ifile==1) then
!        ilo=it0
!      else
!        ilo=1
!      endif !ifile
      ilo=1

      do it2=1,1 !ntime
!        print*, 'Time out (days)=',timeout/86400,it2

!       Read T,S,u,v,SSH
!       WARNING! Make sure the order of vertical indices is 1 
!                at bottom, ilen at surface; revert if necessary!
        status = nf90_get_var(sid,svid,salt(:,:,ilen:1:-1),start=(/1,1,1,it2/),count=(/ixlen,iylen,ilen,1/))
        status = nf90_get_var(sid,tvid,temp(:,:,ilen:1:-1),start=(/1,1,1,it2/),count=(/ixlen,iylen,ilen,1/))
        status = nf90_get_var(ncids(2),uvid,uvel(:,:,ilen:1:-1),start=(/1,1,1,it2/),count=(/ixlen,iylen,ilen,1/))
        status = nf90_get_var(ncids(2),vvid,vvel(:,:,ilen:1:-1),start=(/1,1,1,it2/),count=(/ixlen,iylen,ilen,1/))
        status = nf90_get_var(ncids(1),evid,ssh,start=(/1,1,it2/),count=(/ixlen,iylen,1/))

        !Scaling etc
        salt=salt*1.e-3+20
        temp=temp*1.e-3+20
        uvel=uvel*1.e-3
        vvel=vvel*1.e-3
        ssh=ssh*1.e-3
        !Define junk value for sid; the test is salt<rjunk+0.1
        rjunk=-3.e4*1.e-3+20

!       Assume these won't change over time iteration!!
!       Find bottom index and extend
        ndrypt=0 !# of dry nodes in nc
        do i=1,ixlen
          do j=1,iylen
!             if(salt(i,j,ilen)<rjunk+0.1) then
            if(ssh(i,j)<rjunk+0.1) then !use ssh to judge
              kbp(i,j)=-1 !dry
              ndrypt=ndrypt+1
            else !wet
              !Extend near bottom
              klev0=-1 !flag
              klev=-1 !flag
              do k=1,ilen !uvel
                if(uvel(i,j,k)>rjunk) then
                  klev(1)=k; exit
                endif
              enddo !k
              do k=1,ilen !vvel
                if(vvel(i,j,k)>rjunk) then
                  klev(2)=k; exit
                endif
              enddo !k
              do k=1,ilen !temp
                if(temp(i,j,k)>rjunk) then
                  klev(3)=k; exit
                endif
              enddo !k
              do k=1,ilen !salt
                if(salt(i,j,k)>rjunk) then
                  klev(4)=k; exit
                endif
              enddo !k
              klev0=maxval(klev)
              
!               if(klev0<=0) then
              if(minval(klev)<=0) then
                write(11,*)'Impossible (1):',i,j,klev!salt(i,j,ilen)
                stop
              endif !klev0

!               Fill junk uvel, vvel with zero to eliminate HYCOM nc error
              do k=klev0,ilen
                if (uvel(i,j,k)<rjunk) then
                  write(20,*) 'Warn! Uvel Junk in the middle:',it2,i,j,k
                  uvel(i,j,k)=0.
                end if
                if (vvel(i,j,k)<rjunk) then
                  write(20,*) 'Warn! Vvel Junk in the middle:',it2,i,j,k
                  vvel(i,j,k)=0.
                 end if
              end do !k
!               Fill junk salt ,temp with bottom value
              do k=klev0,ilen
                if (temp(i,j,k)<rjunk) then
                  if(k==klev0) then
                    write(11,*)'Bottom T is junk:',it2,i,j,k,temp(i,j,k) 
                    stop
                  endif
                  write(20,*) 'Warn! Temp Junk in the middle:',it2,i,j,k
                  temp(i,j,k)=temp(i,j,klev0)
                end if
                if (salt(i,j,k)<rjunk) then
                  if(k==klev0) then
                    write(11,*)'Bottom S is junk:',it2,i,j,k,salt(i,j,k) 
                    stop
                  endif
                  write(20,*) 'Warn! Salt Junk in the middle:',it2,i,j,k
                  salt(i,j,k)=salt(i,j,klev0)
                end if
              end do !k

              salt(i,j,1:klev0-1)=salt(i,j,klev0)
              temp(i,j,1:klev0-1)=temp(i,j,klev0)
              uvel(i,j,1:klev0-1)=uvel(i,j,klev0)
              vvel(i,j,1:klev0-1)=vvel(i,j,klev0)             
              kbp(i,j)=klev0 !>0; <=ilen

              !Check
              do k=1,ilen
                if(salt(i,j,k)<saltmin.or.salt(i,j,k)>saltmax.or. &
                  &temp(i,j,k)<tempmin.or.temp(i,j,k)>tempmax.or. &
                  &abs(uvel(i,j,k))>vmag_max.or.abs(vvel(i,j,k))>vmag_max.or. &
                  &abs(ssh(i,j))>ssh_max) then
                  write(11,*)'Fatal: no valid values:',it2,i,j,k,salt(i,j,k), &
   &temp(i,j,k),uvel(i,j,k),vvel(i,j,k),ssh(i,j)
                  stop
                endif
              enddo !k
            endif !salt
          enddo !j
        enddo !i  
!       print*, 'ndrypt=',ndrypt

!       Find a valid T,S at each lat point for extrapolation later
        ts_lat=-9999 !flag
        do j=1,iylen
          ifl=0 !flag
          do i=1,ixlen
            if(salt(i,j,ilen)>rjunk) then
              ifl=1
              ts_lat(j,:,1)=temp(i,j,:)
              ts_lat(j,:,2)=salt(i,j,:)
              exit
            endif
          enddo !i
          if(ifl==0) write(20,*)'Did not find a valid T,S at lat pt:',j,lat(1,j) !,salt(:,j,ilen)
        enddo !j

!       Further extend to fill rest of lat pts
        do j=1,iylen
          if(ts_lat(j,1,1)<rjunk+0.1) then
            ifl=0
            do jj=j+1,iylen
              if(ts_lat(jj,1,1)>rjunk) then
                ts_lat(j,:,1:2)=ts_lat(jj,:,1:2) 
                ifl=1; exit 
              endif   
            enddo !jj
            if(ifl==0) then
              write(11,*)'Failed to find a valid T,S at lat pt:',j,lat(1,j)
              stop
            endif
          endif !ts_lat
        enddo !j

!       Compute S,T etc@ invalid pts based on nearest neighbor
!       Search around neighborhood of a pt
        allocate(iparen_of_dry(2,max(1,ndrypt)))
        call cpu_time(tt0)
      
        icount=0
        do i=1,ixlen
          do j=1,iylen
            if(kbp(i,j)==-1) then !invalid pts  
              salt(i,j,1:ilen)=ts_lat(j,:,2)
              temp(i,j,1:ilen)=ts_lat(j,:,1)
              uvel(i,j,1:ilen)=0.
              vvel(i,j,1:ilen)=0.
              ssh(i,j)=0.


!                icount=icount+1
!                if(icount>ndrypt) stop 'overflow'
!                !Compute max possible tier # for neighborhood
!                mmax=max(i-1,ixlen-i,j-1,iylen-j)
!
!                m=0 !tier #
!                loop6: do
!                  m=m+1
!                  do ii=max(-m,1-i),min(m,ixlen-i)
!                    i3=max(1,min(ixlen,i+ii))
!                    do jj=max(-m,1-j),min(m,iylen-j)
!                      j3=max(1,min(iylen,j+jj))
!                      if(kbp(i3,j3)>0) then !found
!                        i1=i3; j1=j3
!                        exit loop6   
!                      endif
!                    enddo !jj
!                  enddo !ii
!
!                  if(m==mmax) then
!                    write(11,*)'Max. exhausted:',i,j,mmax
!                    write(11,*)'kbp'
!                    do ii=1,ixlen
!                      do jj=1,iylen
!                        write(11,*)ii,jj,kbp(ii,jj)
!                      enddo !jj
!                    enddo !ii
!                    stop
!                  endif
!                end do loop6
!
!                salt(i,j,1:ilen)=salt(i1,j1,1:ilen)
!                temp(i,j,1:ilen)=temp(i1,j1,1:ilen)
!                uvel(i,j,1:ilen)=uvel(i1,j1,1:ilen)
!                vvel(i,j,1:ilen)=vvel(i1,j1,1:ilen)
!                ssh(i,j)=ssh(i1,j1)
!
!                !Save for other steps
!                iparen_of_dry(1,icount)=i1
!                iparen_of_dry(2,icount)=j1

              !Check 
              do k=1,ilen
                if(salt(i,j,k)<saltmin.or.salt(i,j,k)>saltmax.or. &
                   temp(i,j,k)<tempmin.or.temp(i,j,k)>tempmax.or. &
                   abs(uvel(i,j,k))>vmag_max.or.abs(vvel(i,j,k))>vmag_max.or. &
                   abs(ssh(i,j))>ssh_max) then
                  write(11,*)'Fatal: no valid values after searching:',it2,i,j,k,salt(i,j,k), &
     &temp(i,j,k),uvel(i,j,k),vvel(i,j,k),ssh(i,j)
                  stop
                endif
              enddo !k

!                write(12,*)icount,i,j,iparen_of_dry(1:2,icount)
            endif !kbp(i,j)==-1
          enddo !j=iylen1,iylen2
        enddo !i=ixlen1,ixlen2
        call cpu_time(tt1)
!          if(icount/=ndrypt) stop 'mismatch(7)'
!          write(20,*)'extending took (sec):',tt1-tt0,ndrypt
!          call flush(20)

!         Test outputs
        icount=0
        do i=1,ixlen
          do j=1,iylen
            icount=icount+1
            write(95,*)icount,lon(i,j),lat(i,j),ssh(i,j)
            write(96,*)icount,lon(i,j),lat(i,j),uvel(i,j,ilen)
            write(97,*)icount,lon(i,j),lat(i,j),vvel(i,j,ilen)
            write(98,*)icount,lon(i,j),lat(i,j),salt(i,j,1) !,i,j
            write(99,*)icount,lon(i,j),lat(i,j),temp(i,j,ilen)
            write(100,*)icount,lon(i,j),lat(i,j),temp(i,j,1)
          enddo !j
        enddo !i
        print*, 'done outputting test outputs for nc'
   
!       Find parent elements for hgrid.ll
        call cpu_time(tt0)
        loop4: do i=1,np
          ixy(i,1:2)=0

          if(interp_mode==0) then !SG search
            do ix=1,ixlen-1
              if(xl(i)>=xind(ix).and.xl(i)<=xind(ix+1)) then
                !Lower left corner index
                ixy(i,1)=ix
                xrat=(xl(i)-xind(ix))/(xind(ix+1)-xind(ix))
                exit
              endif
            enddo !ix
            do iy=1,iylen-1
              if(yl(i)>=yind(iy).and.yl(i)<=yind(iy+1)) then
                !Lower left corner index
                ixy(i,2)=iy
                yrat=(yl(i)-yind(iy))/(yind(iy+1)-yind(iy))
                exit
              endif
            enddo !ix

            if(ixy(i,1)==0.or.ixy(i,2)==0) then
              write(11,*)'Did not find parent:',i,ixy(i,1:2),xl(i),yl(i)
              stop
            endif
            if(xrat<0.or.xrat>1.or.yrat<0.or.yrat>1) then
              write(11,*)'Ratio out of bound:',i,xl(i),yl(i),xrat,yrat
              stop
            endif

            !Bilinear shape function
            arco(1,i)=(1-xrat)*(1-yrat)
            arco(2,i)=xrat*(1-yrat)
            arco(4,i)=(1-xrat)*yrat
            arco(3,i)=xrat*yrat
          else !interp_mode=1; generic search with UG
            do ix=1,ixlen-1 
              do iy=1,iylen-1 
                x1=lon(ix,iy); x2=lon(ix+1,iy); x3=lon(ix+1,iy+1); x4=lon(ix,iy+1)
                y1=lat(ix,iy); y2=lat(ix+1,iy); y3=lat(ix+1,iy+1); y4=lat(ix,iy+1)
                a1=abs(signa_single(xl(i),x1,x2,yl(i),y1,y2))
                a2=abs(signa_single(xl(i),x2,x3,yl(i),y2,y3))
                a3=abs(signa_single(xl(i),x3,x4,yl(i),y3,y4))
                a4=abs(signa_single(xl(i),x4,x1,yl(i),y4,y1))
                b1=abs(signa_single(x1,x2,x3,y1,y2,y3))
                b2=abs(signa_single(x1,x3,x4,y1,y3,y4))
                rat=abs(a1+a2+a3+a4-b1-b2)/(b1+b2)
                if(rat<small1) then
                  ixy(i,1)=ix; ixy(i,2)=iy
!                 Find a triangle
                  in=0 !flag
                  do l=1,2
                    ap=abs(signa_single(xl(i),x1,x3,yl(i),y1,y3))
                    if(l==1) then !nodes 1,2,3
                      bb=abs(signa_single(x1,x2,x3,y1,y2,y3))
                      wild(l)=abs(a1+a2+ap-bb)/bb
                      if(wild(l)<small1*5) then
                        in=1
                        arco(1,i)=max(0.,min(1.,a2/bb))
                        arco(2,i)=max(0.,min(1.,ap/bb))
                        arco(3,i)=max(0.,min(1.,1-arco(1,i)-arco(2,i)))
                        arco(4,i)=0.
                        exit
                      endif
                    else !nodes 1,3,4
                      bb=abs(signa_single(x1,x3,x4,y1,y3,y4))
                      wild(l)=abs(a3+a4+ap-bb)/bb
                      if(wild(l)<small1*5) then
                        in=2
                        arco(1,i)=max(0.,min(1.,a3/bb))
                        arco(3,i)=max(0.,min(1.,a4/bb))
                        arco(4,i)=max(0.,min(1.,1-arco(1,i)-arco(3,i)))
                        arco(2,i)=0.
                        exit
                      endif
                    endif
                  enddo !l=1,2
                  if(in==0) then
                    write(11,*)'Cannot find a triangle:',(wild(l),l=1,2)
                    stop
                  endif
                  !ixy(i,3)=in
                  cycle loop4
                endif !rat<small1
              enddo !iy=iylen1,iylen2-1
            enddo !ix=ixlen1,ixlen2-1
          endif !interp_mode
        end do loop4 !i=1,np

        call cpu_time(tt1)
        write(20,*)'weights took (sec):',tt1-tt0
        call flush(20)

    
!       Do interpolation
        tempout=-99; saltout=-99
        do i=1,np
          if(ixy(i,1)==0.or.ixy(i,2)==0) then
            write(20,*)'Cannot find a parent element:',i
            tempout(:,i)=tem_outside
            saltout(:,i)=sal_outside
            uout(:,i)=0
            vout(:,i)=0
            eout(i)=0
          else !found parent
            ix=ixy(i,1); iy=ixy(i,2) !; in=ixy(i,3)
            !Find vertical level
            do k=1,nvrt
              if(kbp(ix,iy)==-1) then
                lev=ilen-1; vrat=1
              else if(z(k,i)<=zm(ix,iy,kbp(ix,iy))) then
                !lev=kbp(ix,iy); vrat=0
                lev=1; vrat=0
              else if(z(k,i)>=zm(ix,iy,ilen)) then !above f.s.
                lev=ilen-1; vrat=1
              else
                lev=-99 !flag
                do kk=1,ilen-1
                  if(z(k,i)>=zm(ix,iy,kk).and.z(k,i)<=zm(ix,iy,kk+1)) then
                    lev=kk
                    vrat=(zm(ix,iy,kk)-z(k,i))/(zm(ix,iy,kk)-zm(ix,iy,kk+1))
                    exit
                  endif
                enddo !kk
                if(lev==-99) then
                  write(11,*)'Cannot find a level:',i,k,z(k,i),(zm(ix,iy,l),l=1,kbp(ix,iy))
                  stop
                endif
              endif
          
!              write(18,*)i,k,ix,iy,lev,vrat,kbp(ix,iy)
              if(lev>=ilen) then
                write(11,*)'lev:',lev,ix,iy,k,i,ilen,vrat,kbp(ix,iy),z(k,i),zm(ix,iy,1:ilen)
                stop
              endif
              wild2(1,1)=temp(ix,iy,lev)*(1-vrat)+temp(ix,iy,lev+1)*vrat
              wild2(1,2)=salt(ix,iy,lev)*(1-vrat)+salt(ix,iy,lev+1)*vrat
              wild2(2,1)=temp(ix+1,iy,lev)*(1-vrat)+temp(ix+1,iy,lev+1)*vrat
              wild2(2,2)=salt(ix+1,iy,lev)*(1-vrat)+salt(ix+1,iy,lev+1)*vrat
              wild2(3,1)=temp(ix+1,iy+1,lev)*(1-vrat)+temp(ix+1,iy+1,lev+1)*vrat
              wild2(3,2)=salt(ix+1,iy+1,lev)*(1-vrat)+salt(ix+1,iy+1,lev+1)*vrat
              wild2(4,1)=temp(ix,iy+1,lev)*(1-vrat)+temp(ix,iy+1,lev+1)*vrat
              wild2(4,2)=salt(ix,iy+1,lev)*(1-vrat)+salt(ix,iy+1,lev+1)*vrat

              wild2(5,1)=uvel(ix,iy,lev)*(1-vrat)+uvel(ix,iy,lev+1)*vrat
              wild2(5,2)=vvel(ix,iy,lev)*(1-vrat)+vvel(ix,iy,lev+1)*vrat
              wild2(6,1)=uvel(ix+1,iy,lev)*(1-vrat)+uvel(ix+1,iy,lev+1)*vrat
              wild2(6,2)=vvel(ix+1,iy,lev)*(1-vrat)+vvel(ix+1,iy,lev+1)*vrat
              wild2(7,1)=uvel(ix+1,iy+1,lev)*(1-vrat)+uvel(ix+1,iy+1,lev+1)*vrat
              wild2(7,2)=vvel(ix+1,iy+1,lev)*(1-vrat)+vvel(ix+1,iy+1,lev+1)*vrat
              wild2(8,1)=uvel(ix,iy+1,lev)*(1-vrat)+uvel(ix,iy+1,lev+1)*vrat
              wild2(8,2)=vvel(ix,iy+1,lev)*(1-vrat)+vvel(ix,iy+1,lev+1)*vrat

              tempout(k,i)=dot_product(wild2(1:4,1),arco(1:4,i))
              saltout(k,i)=dot_product(wild2(1:4,2),arco(1:4,i))
              uout(k,i)=dot_product(wild2(5:8,1),arco(1:4,i))
              vout(k,i)=dot_product(wild2(5:8,2),arco(1:4,i))
              eout(i)=ssh(ix,iy)*arco(1,i)+ssh(ix+1,iy)*arco(2,i)+ssh(ix+1,iy+1)*arco(3,i)+ssh(ix,iy+1)*arco(4,i)

              !Check
              if(tempout(k,i)<tempmin.or.tempout(k,i)>tempmax.or. &
                saltout(k,i)<saltmin.or.saltout(k,i)>saltmax.or. &
                abs(uout(k,i))>vmag_max.or.abs(vout(k,i))>vmag_max.or.abs(eout(i))>ssh_max) then
                write(11,*)'Interpolated values invalid:',i,k,tempout(k,i),saltout(k,i),uout(k,i), &
     &vout(k,i),eout(i),temp(ix,iy,lev)
                stop
              endif

!             Enforce lower bound for temp. for eqstate
              tempout(k,i)=max(-10.,tempout(k,i))
              saltout(k,i)=min(50.,saltout(k,i))

!             Enforce lower bound for salt (this is the only occurence in the code)
!             if(z(k,i)<=-ht) saltout(k,i)=max(saltout(k,i),smin)
            enddo !k=1,nvrt

            !Estuary pts
            if(iest(i)==1) then
              tempout(:,i)=tem_es
              saltout(:,i)=sal_es
            endif
          endif !ixy(i,1)==0.or.
        enddo !i

!       hotstart.nc
        if(ifile==1.and.it2==ilo) then
          write(20,*)'outputting hot...'
          call flush(20)

          do i=1,ns
            n1=isidenode(1,i)
            n2=isidenode(2,i)
            do k=1,nvrt
!              tsd(i,k)=(tempout(n1,k)+tempout(n2,k))/2
!              ssd(i,k)=(saltout(n1,k)+saltout(n2,k))/2
              if(iuv==0) then
                su2(k,i)=0; sv2(k,i)=0
              else
                su2(k,i)=(uout(k,n1)+uout(k,n2))/2
                sv2(k,i)=(vout(k,n1)+vout(k,n2))/2
              endif !iuv
            enddo !k
!           write(88,*)i,xcj(i),ycj(i),ssd(i,1),ssd(i,nvrt)
          enddo !i

          do i=1,ne
            do k=2,nvrt
              tsel(1,k,i)=(sum(tempout(k,elnode(1:i34(i),i)))+sum(tempout(k-1,elnode(1:i34(i),i))))/2/i34(i) 
              tsel(2,k,i)=(sum(saltout(k,elnode(1:i34(i),i)))+sum(saltout(k-1,elnode(1:i34(i),i))))/2/i34(i)
            enddo !k
            tsel(1,1,i)=tsel(1,2,i) !mainly for hotstart format
            tsel(2,1,i)=tsel(2,2,i)
          enddo !i

          if(iuv==0) then
            eout_tmp=0
          else
            eout_tmp=eout
          endif

!         Debug
          do i=1,np
            write(26,*)i,xl(i),yl(i),saltout(nvrt,i)
            write(30,*)i,xl(i),yl(i),tempout(nvrt,i)
            write(21,*)i,xl(i),yl(i),tempout(1,i)
            write(24,*)xl(i),yl(i),uout(1,i),vout(1,i)
            write(27,*)xl(i),yl(i),uout(nvrt,i),vout(nvrt,i)
            write(25,*)i,xl(i),yl(i),eout(i)
!            write(29,*)'T profile at node ',i,dp(i)
!            do k=1,nvrt
!              write(29,*)k,z(k,i),tempout(k,i)
!            enddo !k
          enddo !i
          do i=1,ne
            write(22,*)i,tsel(1,nvrt,i) !T
            write(23,*)i,tsel(1,1,i)
          enddo !i

          iret=nf90_create('hotstart.nc',OR(NF90_NETCDF4,NF90_CLOBBER),ncids(4))
          iret=nf90_def_dim(ncids(4),'node',np,node_dim)
          iret=nf90_def_dim(ncids(4),'elem',ne,nele_dim)
          iret=nf90_def_dim(ncids(4),'side',ns,nedge_dim)
          iret=nf90_def_dim(ncids(4),'nVert',nvrt,nv_dim)
          iret=nf90_def_dim(ncids(4),'ntracers',2,ntr_dim)
          iret=nf90_def_dim(ncids(4),'one',1,one_dim)
          var1d_dims(1)=one_dim
          iret=nf90_def_var(ncids(4),'time',NF90_DOUBLE,var1d_dims,ivarid(1))
          iret=nf90_def_var(ncids(4),'iths',NF90_INT,var1d_dims,ivarid(2))
          iret=nf90_def_var(ncids(4),'ifile',NF90_INT,var1d_dims,ivarid(3))
          var1d_dims(1)=nele_dim
          iret=nf90_def_var(ncids(4),'idry_e',NF90_INT,var1d_dims,ivarid(4))
          var1d_dims(1)=nedge_dim
          iret=nf90_def_var(ncids(4),'idry_s',NF90_INT,var1d_dims,ivarid(5))
          var1d_dims(1)=node_dim
          iret=nf90_def_var(ncids(4),'idry',NF90_INT,var1d_dims,ivarid(6))
          iret=nf90_def_var(ncids(4),'eta2',NF90_DOUBLE,var1d_dims,ivarid(7))
          var2d_dims(1)=nv_dim; var2d_dims(2)=nele_dim
          iret=nf90_def_var(ncids(4),'we',NF90_DOUBLE,var2d_dims,ivarid(8))
          var3d_dims(1)=ntr_dim; var3d_dims(2)=nv_dim; var3d_dims(3)=nele_dim
          iret=nf90_def_var(ncids(4),'tr_el',NF90_DOUBLE,var3d_dims,ivarid(9))
          var3d_dims(1)=ntr_dim; var3d_dims(2)=nv_dim; var3d_dims(3)=node_dim
          iret=nf90_def_var(ncids(4),'tr_nd',NF90_DOUBLE,var3d_dims,ivarid(10))
          iret=nf90_def_var(ncids(4),'tr_nd0',NF90_DOUBLE,var3d_dims,ivarid(11))
          var2d_dims(1)=nv_dim; var2d_dims(2)=nedge_dim
          iret=nf90_def_var(ncids(4),'su2',NF90_DOUBLE,var2d_dims,ivarid(12))
          iret=nf90_def_var(ncids(4),'sv2',NF90_DOUBLE,var2d_dims,ivarid(13))
          var2d_dims(1)=nv_dim; var2d_dims(2)=node_dim
          iret=nf90_def_var(ncids(4),'q2',NF90_DOUBLE,var2d_dims,ivarid(14))
          iret=nf90_def_var(ncids(4),'xl',NF90_DOUBLE,var2d_dims,ivarid(15))
          iret=nf90_def_var(ncids(4),'dfv',NF90_DOUBLE,var2d_dims,ivarid(16))
          iret=nf90_def_var(ncids(4),'dfh',NF90_DOUBLE,var2d_dims,ivarid(17))
          iret=nf90_def_var(ncids(4),'dfq1',NF90_DOUBLE,var2d_dims,ivarid(18))
          iret=nf90_def_var(ncids(4),'dfq2',NF90_DOUBLE,var2d_dims,ivarid(19))
          iret=nf90_enddef(ncids(4))
   
          idry_s=0; zeros=0
          iret=nf90_put_var(ncids(4),ivarid(1),dble(0.))
          iret=nf90_put_var(ncids(4),ivarid(2),0)
          iret=nf90_put_var(ncids(4),ivarid(3),1)
          iret=nf90_put_var(ncids(4),ivarid(4),idry_s(1:ne),(/1/),(/ne/))
          iret=nf90_put_var(ncids(4),ivarid(5),idry_s(1:ns),(/1/),(/ns/))
          iret=nf90_put_var(ncids(4),ivarid(6),idry_s(1:np),(/1/),(/np/))
          iret=nf90_put_var(ncids(4),ivarid(7),dble(eout_tmp(1:np)),(/1/),(/np/))
          iret=nf90_put_var(ncids(4),ivarid(8),dble(zeros(1:nvrt,1:ne)),(/1,1/),(/nvrt,ne/))
          iret=nf90_put_var(ncids(4),ivarid(9),dble(tsel(1:2,1:nvrt,1:ne)))
          iret=nf90_put_var(ncids(4),ivarid(10),dble(tempout(1:nvrt,1:np)),(/1,1,1/),(/1,nvrt,np/))
          iret=nf90_put_var(ncids(4),ivarid(10),dble(saltout(1:nvrt,1:np)),(/2,1,1/),(/1,nvrt,np/))
          iret=nf90_put_var(ncids(4),ivarid(11),dble(tempout(1:nvrt,1:np)),(/1,1,1/),(/1,nvrt,np/))
          iret=nf90_put_var(ncids(4),ivarid(11),dble(saltout(1:nvrt,1:np)),(/2,1,1/),(/1,nvrt,np/))
          iret=nf90_put_var(ncids(4),ivarid(12),dble(su2(1:nvrt,1:ns)))
          iret=nf90_put_var(ncids(4),ivarid(13),dble(sv2(1:nvrt,1:ns)))
          iret=nf90_put_var(ncids(4),ivarid(14),dble(zeros(1:nvrt,1:np)))
          iret=nf90_put_var(ncids(4),ivarid(15),dble(zeros(1:nvrt,1:np)))
          iret=nf90_put_var(ncids(4),ivarid(16),dble(zeros(1:nvrt,1:np)))
          iret=nf90_put_var(ncids(4),ivarid(17),dble(zeros(1:nvrt,1:np)))
          iret=nf90_put_var(ncids(4),ivarid(18),dble(zeros(1:nvrt,1:np)))
          iret=nf90_put_var(ncids(4),ivarid(19),dble(zeros(1:nvrt,1:np)))

          iret=nf90_close(ncids(4))
        endif !ifile; hotstart

!        timeout=timeout+dtout !sec

!        write(20,*)'done time:',timeout/86400
!        call flush(20)

!        if(timeout/86400>nndays) stop 'finished'

      enddo !it2=1,1

      status = nf90_close(sid)
      status = nf90_close(ncids(1))
      status = nf90_close(ncids(2))
!end11
      print*, 'done reading nc for file ',ifile
      print*, 'Finished!'
      stop
!--------------------------------------------------------------------
      enddo !ifile=1,

!      deallocate(lat,lon,zm,h,kbp,ihope,xind,yind,lind,salt,temp)

!      print*, 'Finished'
!     End of main
!     Subroutines
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

      end program gen_hot

