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

! Generate boundary conditions (*[23]D.th.nc) from gridded HYCOM data (nc file)
! The section on nc read needs to be modified as appropriate- search for
! 'start11' and 'end11'
! Beware the order of vertical levels in the nc file!!!
! Assume elev=0 in interpolation of 3D variables
! The code will do all it can to extrapolate: below bottom/above surface.
! If a bnd pt is outside the background nc grid, const. values will be filled (from gen_3Dth_from_nc.in) 

! Dan Yu: added checking vertical layer for each 3D var from HYCOM nc, and some easy remedy for nc file error  
! Wet-dry points are defined by ssh.
! For U,V, if junk values in the middle of water, simply fill with 0 and record in fort.20
! For S,T, if junk values in the middle of water, simply fill with bottom value, and record in fort.20
! If SSH shows wet/dry in time, search nearby points to fill. If none found, fill 0, and record in fort.20
! Consider checking HYCOM nc files with ncview or other tools.

!   Tip: to speed up, reduce the extent of HYCOM inputs to only cover
!   the open bnd.

!   Input: 
!     (1) hgrid.gr3;
!     (2) hgrid.ll;
!     (3) vgrid.in (SCHISM R1703 and up);
!     (4) gen_3Dth_from_nc.in: 
!                     1st line: T,S values for pts outside bg grid in nc
!                     2nd line: time step in .nc in sec
!                     3rd line: nob, iob(1:nob) - # of open bnd seg's that need
!                               *[23D].th; list of seg IDs. All *[23D].th must share same set of bnd seg's
!                     4th line: # of days needed
!                     5th line: # of HYCOM file stacks
!     (5) HYCOM files: [SSH,TS,UV]_[1,2,..nfiles].nc (beware scaling etc)
!   Output: *[23D].th.nc
!   Debug outputs: fort.11 (fatal errors); fort.20 (warning); fort.2[1-9], fort.9[5-9], fort.100; backup.out

! ifort -O2 -mcmodel=medium -assume byterecl -CB -o gen_3Dth_from_hycom.exe ../UtilLib/schism_geometry.f90 \
! ../UtilLib/extract_mod.f90 ../UtilLib/compute_zcor.f90 ../UtilLib/pt_in_poly_test.f90 gen_3Dth_from_hycom.f90 \
!-I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff

      program gen_bc
      use netcdf
      use schism_geometry_mod
      use compute_zcor
      use pt_in_poly_test, only: signa_single

!      implicit none

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
      integer :: lldim ! check lat/lon dimension
             
!     Local variables for data
      real (kind = 4), allocatable :: xind(:), yind(:), lind(:), xind2(:,:), yind2(:,:)
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

      integer, allocatable :: elnode(:,:),elside(:,:),isdel(:,:),idry_s(:),i34(:), &
     &ixy(:,:),ic3(:,:),isidenode(:,:),nond(:),iond(:,:),iob(:),iond2(:)
      !integer :: indel(mnei,mnp),idry_s(mns)
      !integer :: idry_s(mns)
      real, allocatable :: xl(:),yl(:),dp(:),ztot(:),sigma(:),arco(:,:), &
     &tempout(:,:),saltout(:,:),uout(:,:),vout(:,:),eout(:),su2(:,:),sv2(:,:), &
     &eout_tmp(:),tsel(:,:,:),zeros(:,:),xcj(:),ycj(:)
      allocatable :: z(:,:),sigma_lcl(:,:),kbp2(:),iparen_of_dry(:,:)
      real*8 :: aa1(1)

!     First statement
!     Currently we assume rectangular grid in HYCOM
!     (interp_mode=0). interp_mode=1 uses generic UG search (splitting
!     quads) and is kept for more generic cases
      interp_mode=0

      open(10,file='gen_3Dth_from_nc.in',status='old')
      read(10,*) tem_outside,sal_outside !T,S values for pts outside bg grid in nc
      read(10,*) dtout !time step in .nc [sec]
      read(10,*) nob !,iob(1:nob) !# of open bnds that need *3D.th; list of IDs
      read(10,*) nndays !# of days needed in output
      read(10,*) nfiles !# of stacks of HYCOM files

      allocate(iob(nob))
      rewind(10)
      do i=1,2; read(10,*); enddo
      read(10,*) nob,iob(1:nob)
      close(10)
      if(tem_outside<0.or.sal_outside<0) &
     &stop 'Invalid T,S constants'

!     Read in hgrid and vgrid
      open(16,file='hgrid.ll',status='old')
      open(14,file='hgrid.gr3',status='old') !only need depth info and connectivity
      open(19,file='vgrid.in',status='old')
      open(11,file='fort.11',status='replace')
      read(14,*)
      read(14,*)ne,np
      allocate(xl(np),yl(np),dp(np),i34(ne),elnode(4,ne),ixy(np,3),stat=istat)
      if(istat/=0) stop 'Failed to alloc. (0)'
      read(16,*); read(16,*)
      do i=1,np
        read(14,*)j,xtmp,ytmp,dp(i)
        read(16,*)j,xl(i),yl(i) 
      enddo !i
      do i=1,ne
        read(14,*)j,i34(i),(elnode(l,i),l=1,i34(i))
      enddo !i
!     Open bnds
      read(14,*) nope
      read(14,*) neta
      ntot=0
      allocate(nond(max(1,nope)))
      do k=1,nope
        read(14,*) nond(k)
        do i=1,nond(k)
          read(14,*) !iond(k,i)
        enddo
      enddo

      rewind(14)
      mnope=maxval(nond)
      if(mnope<=0) stop 'no open bnd found'
      allocate(iond(mnope,nope),iond2(mnope*nope))
      do i=1,2+np+ne+2; read(14,*); enddo
      do k=1,nope
        read(14,*) !nond(k)
        do i=1,nond(k)
          read(14,*)iond(i,k)
        enddo
      enddo
      close(14)

      nond0=0
      do i=1,nob
        ibnd=iob(i)
        do j=1,nond(ibnd)
          nond0=nond0+1
          iond2(nond0)=iond(j,ibnd)
        enddo !j
      enddo !i

      close(16)

!     V-grid
      read(19,*)ivcor
      read(19,*)nvrt
!      if(nvrt>mnv) stop 'increase mnv'
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

      iret=nf90_create('elev2D.th.nc',OR(NF90_NETCDF4,NF90_CLOBBER),ncids(5))
      iret=nf90_def_dim(ncids(5),'nOpenBndNodes',nond0,nond0_dim)
      iret=nf90_def_dim(ncids(5),'one',1,one_dim)
      iret=nf90_def_dim(ncids(5),'time', NF90_UNLIMITED,itime_dim)
      iret=nf90_def_dim(ncids(5),'nLevels',1,nvrt_dim)
      iret=nf90_def_dim(ncids(5),'nComponents',1,ivs_dim)
      var1d_dims(1)=one_dim
      iret=nf90_def_var(ncids(5),'time_step',NF90_FLOAT,var1d_dims,idt)
      var1d_dims(1)=itime_dim
      iret=nf90_def_var(ncids(5),'time',NF90_DOUBLE,var1d_dims,itime_id(1))
      var4d_dims(1)=ivs_dim; var4d_dims(2)=nvrt_dim; var4d_dims(3)=nond0_dim
      var4d_dims(4)=itime_dim
      iret=nf90_def_var(ncids(5),'time_series',NF90_FLOAT,var4d_dims,ivarid(51))
      iret=nf90_enddef(ncids(5))
      iret=nf90_put_var(ncids(5),idt,dtout)

      iret=nf90_create('uv3D.th.nc',OR(NF90_NETCDF4,NF90_CLOBBER),ncids(6))
      iret=nf90_def_dim(ncids(6),'nOpenBndNodes',nond0,nond0_dim)
      iret=nf90_def_dim(ncids(6),'one',1,one_dim)
      iret=nf90_def_dim(ncids(6),'time', NF90_UNLIMITED,itime_dim)
      iret=nf90_def_dim(ncids(6),'nLevels',nvrt,nvrt_dim)
      iret=nf90_def_dim(ncids(6),'nComponents',2,ivs_dim)
      var1d_dims(1)=one_dim
      iret=nf90_def_var(ncids(6),'time_step',NF90_FLOAT,var1d_dims,idt)
      var1d_dims(1)=itime_dim
      iret=nf90_def_var(ncids(6),'time',NF90_DOUBLE,var1d_dims,itime_id(2))
      var4d_dims(1)=ivs_dim; var4d_dims(2)=nvrt_dim; var4d_dims(3)=nond0_dim
      var4d_dims(4)=itime_dim
      iret=nf90_def_var(ncids(6),'time_series',NF90_FLOAT,var4d_dims,ivarid(52))
      iret=nf90_enddef(ncids(6))
      iret=nf90_put_var(ncids(6),idt,dtout)

      iret=nf90_create('TEM_3D.th.nc',OR(NF90_NETCDF4,NF90_CLOBBER),ncids(7))
      iret=nf90_def_dim(ncids(7),'nOpenBndNodes',nond0,nond0_dim)
      iret=nf90_def_dim(ncids(7),'one',1,one_dim)
      iret=nf90_def_dim(ncids(7),'time', NF90_UNLIMITED,itime_dim)
      iret=nf90_def_dim(ncids(7),'nLevels',nvrt,nvrt_dim)
      iret=nf90_def_dim(ncids(7),'nComponents',1,ivs_dim)
      var1d_dims(1)=one_dim
      iret=nf90_def_var(ncids(7),'time_step',NF90_FLOAT,var1d_dims,idt)
      var1d_dims(1)=itime_dim
      iret=nf90_def_var(ncids(7),'time',NF90_DOUBLE,var1d_dims,itime_id(3))
      var4d_dims(1)=ivs_dim; var4d_dims(2)=nvrt_dim; var4d_dims(3)=nond0_dim
      var4d_dims(4)=itime_dim
      iret=nf90_def_var(ncids(7),'time_series',NF90_FLOAT,var4d_dims,ivarid(53))
      iret=nf90_enddef(ncids(7))
      iret=nf90_put_var(ncids(7),idt,dtout)

      iret=nf90_create('SAL_3D.th.nc',OR(NF90_NETCDF4,NF90_CLOBBER),ncids(8))
      iret=nf90_def_dim(ncids(8),'nOpenBndNodes',nond0,nond0_dim)
      iret=nf90_def_dim(ncids(8),'one',1,one_dim)
      iret=nf90_def_dim(ncids(8),'time', NF90_UNLIMITED,itime_dim)
      iret=nf90_def_dim(ncids(8),'nLevels',nvrt,nvrt_dim)
      iret=nf90_def_dim(ncids(8),'nComponents',1,ivs_dim)
      var1d_dims(1)=one_dim
      iret=nf90_def_var(ncids(8),'time_step',NF90_FLOAT,var1d_dims,idt)
      var1d_dims(1)=itime_dim
      iret=nf90_def_var(ncids(8),'time',NF90_DOUBLE,var1d_dims,itime_id(4))
      var4d_dims(1)=ivs_dim; var4d_dims(2)=nvrt_dim; var4d_dims(3)=nond0_dim
      var4d_dims(4)=itime_dim
      iret=nf90_def_var(ncids(8),'time_series',NF90_FLOAT,var4d_dims,ivarid(54))
      iret=nf90_enddef(ncids(8))
      iret=nf90_put_var(ncids(8),idt,dtout)

      open(12,file='backup.out',status='replace')

!start11
!     Set time step # (must be in 1st stack) corresponding to t=0
      it0=1

!     Define limits for variables for sanity checks
      tempmin=-10 
      tempmax=45
      saltmin=0
      saltmax=50
      vmag_max=10 !max. |u| or |v|
      ssh_max=6 !max. of |SSH|

!     Assuming elev, u,v, S,T have same dimensions (and time step) and grids do not change over time
      timeout=0
      irecout=1
      do ifile=1,nfiles
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

!       Check static info (lat/lon) dimension & allocate
        status = nf90_inq_varid(ncids(1), "xlon", xvid)
        status = nf90_Inquire_Variable(ncids(1), xvid,ndims = lldim)
        print*, 'xlon, ylat is ', lldim ,' dimension.'
        if (lldim.eq.1) then
           allocate(xind(ixlen),stat=ier)
           allocate(yind(iylen),stat=ier)
        else if (lldim.eq.2) then
           allocate(xind2(ixlen,iylen),stat=ier)
           allocate(yind2(ixlen,iylen),stat=ier)
           interp_mode=1
        else
           print*, 'Error dimension in xlon,ylat!'
           stop
        end if

!       allocate arrays
!       allocate(xind(ixlen),stat=ier)
!       allocate(yind(iylen),stat=ier)
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
        allocate(temp(ixlen,iylen,ilen),stat=ier)
        allocate(ssh(ixlen,iylen),stat=ier)
        uvel=0; vvel=0; ssh=0
  
!       get static info (lat/lon grids etc) 
        if (lldim.eq.1) then
           status = nf90_inq_varid(ncids(1), "xlon", xvid)
           status = nf90_get_var(ncids(1), xvid, xind)
           status = nf90_inq_varid(ncids(1), "ylat", yvid)
           status = nf90_get_var(ncids(1), yvid, yind)
        elseif (lldim.eq.2) then
           status = nf90_inq_varid(ncids(1), "xlon", xvid)
           status = nf90_get_var(ncids(1), xvid, xind2)
           status = nf90_inq_varid(ncids(1), "ylat", yvid)
           status = nf90_get_var(ncids(1), yvid, yind2)
        end if

        !lind may be sigma coord.
!        status = nf90_inq_varid(sid, "sigma", lvid)
!        status = nf90_get_var(sid, lvid, lind)
        status = nf90_inq_varid(sid, "depth", hvid)
        status = nf90_get_var(sid, hvid, hnc)

!       processing static info
        if (lldim.eq.1) then
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
        elseif (lldim.eq.2) then
           lon=xind2
           lat=yind2
           do j=1,iylen
            do i=1,ixlen
             if(i<ixlen) then; if(lon(i,j)>=lon(i+1,j)) then
               write(11,*)'Lon must be increasing:',i,lon(i,j),lon(i+1,j)
               stop
             endif; endif;
            enddo !i
           enddo !j
           do i=1,ixlen
            do j=1,iylen
             if(j<iylen) then; if(lat(i,j)>=lat(i,j+1)) then
               write(11,*)'Lat must be increasing:',i,lat(i,j),lat(i,j+1)
               stop
             endif; endif;
            enddo !j
           enddo !i
        end if !lldim

!       lon=lon-360 !convert to our long.

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
      if(ifile==1) then
        ilo=it0
      else
        ilo=1
      endif !ifile

      do it2=ilo,ntime
        print*, 'Time out (days)=',timeout/86400,it2

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
        !Define junk value for sid; the test is var <rjunk+0.1
        rjunk=-3.e4*1.e-3+20

!       Do sth for 1st step
        if(ifile==1.and.it2==ilo) then
!         Assume these won't change over time iteration!!
!         Find bottom index and extend
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
                  if(uvel(i,j,k)>=rjunk+0.1) then
                    klev(1)=k; exit
                  endif
                enddo !k
                do k=1,ilen !vvel
                  if(vvel(i,j,k)>=rjunk+0.1) then
                    klev(2)=k; exit
                  endif
                enddo !k
                do k=1,ilen !temp
                  if(temp(i,j,k)>=rjunk+0.1) then
                    klev(3)=k; exit
                  endif
                enddo !k
                do k=1,ilen !salt
                  if(salt(i,j,k)>=rjunk+0.1) then
                    klev(4)=k; exit
                  endif
                enddo !k
                klev0=maxval(klev(1:4))
                
!               if(klev0<=0) then
                if(minval(klev)<=0) then
                  write(11,*)'Impossible (1):',i,j,klev !salt(i,j,ilen)
                  stop
                endif !klev0

!               Fill junk uvel, vvel with zero to eliminate HYCOM nc error
                do k=klev0,ilen
                  if (uvel(i,j,k)<rjunk+0.1) then
                    write(20,*) 'Warn! Uvel Junk in the middle:',it2,i,j,k
                    uvel(i,j,k)=0.
                  end if
                  if (vvel(i,j,k)<rjunk+0.1) then
                    write(20,*) 'Warn! Vvel Junk in the middle:',it2,i,j,k
                    vvel(i,j,k)=0.
                   end if
                end do !k

!               Fill junk salt ,temp with bottom value
                do k=klev0,ilen
                  if (temp(i,j,k)<rjunk+0.1) then
                    if(k==klev0) then
                      write(11,*)'Bottom T is junk:',it2,i,j,k,temp(i,j,k) 
                      stop
                    endif
                    write(20,*) 'Warn! Temp Junk in the middle:',it2,i,j,k
                    temp(i,j,k)=temp(i,j,klev0)
                  end if
                  if (salt(i,j,k)<rjunk+0.1) then
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
!          print*, 'ndrypt=',ndrypt

!         Compute S,T etc@ invalid pts based on nearest neighbor
!         Search around neighborhood of a pt
          allocate(iparen_of_dry(2,max(1,ndrypt)))
          call cpu_time(tt0)
        
          icount=0
          do i=1,ixlen
            do j=1,iylen
              if(kbp(i,j)==-1) then !invalid pts  
                icount=icount+1
                if(icount>ndrypt) stop 'overflow'
                !Compute max possible tier # for neighborhood
                mmax=max(i-1,ixlen-i,j-1,iylen-j)

                m=0 !tier #
                loop6: do
                  m=m+1
                  do ii=max(-m,1-i),min(m,ixlen-i)
                    i3=max(1,min(ixlen,i+ii))
                    do jj=max(-m,1-j),min(m,iylen-j)
                      j3=max(1,min(iylen,j+jj))
                      if(kbp(i3,j3)>0) then !found
                        i1=i3; j1=j3
                        exit loop6   
                      endif
                    enddo !jj
                  enddo !ii

                  if(m==mmax) then
                    write(11,*)'Max. exhausted:',i,j,mmax
                    write(11,*)'kbp'
                    do ii=1,ixlen
                      do jj=1,iylen
                        write(11,*)ii,jj,kbp(ii,jj)
                      enddo !jj
                    enddo !ii
                    stop
                  endif
                end do loop6

                salt(i,j,1:ilen)=salt(i1,j1,1:ilen)
                temp(i,j,1:ilen)=temp(i1,j1,1:ilen)
                uvel(i,j,1:ilen)=uvel(i1,j1,1:ilen)
                vvel(i,j,1:ilen)=vvel(i1,j1,1:ilen)
                ssh(i,j)=ssh(i1,j1)

                !Save for other steps
                iparen_of_dry(1,icount)=i1
                iparen_of_dry(2,icount)=j1

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

                write(12,*)icount,i,j,iparen_of_dry(1:2,icount)
              endif !kbp(i,j)==-1
            enddo !j=iylen1,iylen2
          enddo !i=ixlen1,ixlen2
          call cpu_time(tt1)
          if(icount/=ndrypt) stop 'mismatch(7)'
          write(20,*)'extending took (sec):',tt1-tt0,ndrypt
          call flush(20)

          !Save for abnormal cases later
!          salt0=salt
!          temp0=temp         
!          ssh0=ssh
!          uvel0=uvel
!          vvel0=vvel

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
     
!         Find parent elements for bnd pts
          call cpu_time(tt0)
          loop4: do ii=1,nond0 !np
            i=iond2(ii) !global node #

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

              if(ixy(i,1)/=0.and.ixy(i,2)/=0) then !found
!                write(11,*)'Did not find parent:',i,ixy(i,1:2)
!                stop
                if(xrat<0.or.xrat>1.or.yrat<0.or.yrat>1) then
                  write(11,*)'Ratio out of bound:',i,xrat,yrat
                  stop
                endif

                !Bilinear shape function
                arco(1,i)=(1-xrat)*(1-yrat)
                arco(2,i)=xrat*(1-yrat)
                arco(4,i)=(1-xrat)*yrat
                arco(3,i)=xrat*yrat
              endif !ixy

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
!                   Find a triangle
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

        else !skip to extend
          do i=1,ixlen
            do j=1,iylen
              if(kbp(i,j)>0) then !valid
!               Check ssh consistency by kbp
!               If not consistent, seek nearby points to replace
!               Force ssh = 0 if none are found
                if (ssh(i,j)<rjunk+0.1) then
                  do ii=1,2
                    do jj=1,2
                      if(i+ii<=ixlen.and.j+jj<=iylen) then; if (ssh(i+ii,j+jj)>=rjunk+0.1) then
                        ssh(i,j)=ssh(i+ii,j+jj); exit
                      endif; endif
                    end do
                  end do !ii
                  if (ssh(i,j)<rjunk+0.1) then
                    write(20,*) 'Warning! SSH Junk, fill with 0: ',it2,i,j
                    ssh(i,j)=0. !fill with 0 if no finding
                  end if
                end if

                !Extend bottom (kbp changes over time)
                klev0=-1 !flag
                klev=-1 !flag
                do k=1,ilen !uvel
                  if(uvel(i,j,k)>=rjunk+0.1) then
                    klev(1)=k; exit
                  endif
                enddo !k
                do k=1,ilen !vvel
                  if(vvel(i,j,k)>=rjunk+0.1) then
                    klev(2)=k; exit
                  endif
                enddo !k
                do k=1,ilen !temp
                  if(temp(i,j,k)>=rjunk+0.1) then
                    klev(3)=k; exit
                  endif
                enddo !k
                do k=1,ilen !salt
                  if(salt(i,j,k)>=rjunk+0.1) then
                    klev(4)=k; exit
                  endif
                enddo !k
                klev0=maxval(klev)

!               if(klev0<=0) then
                if(minval(klev)<=0) then
                  write(11,*)'Impossible (4):',i,j,klev !salt(i,j,ilen)
                  stop
                endif !klev0

!               Fill junk uvel, vvel with zero to easy eliminate HYCOM nc error
                do k=klev0,ilen
                  if (uvel(i,j,k)<rjunk+0.1) then
                    write(20,*) 'Warn! Uvel Junk in the middle:',it2,i,j,k
                    uvel(i,j,k)=0.
                  endif
                  if (vvel(i,j,k)<rjunk+0.1) then
                    write(20,*) 'Warn! Vvel Junk in the middle:',it2,i,j,k
                    vvel(i,j,k)=0.
                  endif
                end do !k
!               Fill junk salt ,temp with bottom value
                do k=klev0,ilen
                  if (temp(i,j,k)<rjunk+0.1) then
                    if(k==klev0) then
                      write(11,*)'Bottom T is junk(2):',it2,i,j,k,temp(i,j,k)
                      stop
                    endif
                    write(20,*) 'Warn! Temp Junk in the middle:',it2,i,j,k
                    temp(i,j,k)=temp(i,j,klev0)
                  endif
                  if (salt(i,j,k)<rjunk+0.1) then
                    if(k==klev0) then
                      write(11,*)'Bottom S is junk(2):',it2,i,j,k,salt(i,j,k)
                      stop
                    endif

                    write(20,*) 'Warn! Salt Junk in the middle:',it2,i,j,k
                    salt(i,j,k)=salt(i,j,klev0)
                  endif
                end do !k

                salt(i,j,1:klev0-1)=salt(i,j,klev0)
                temp(i,j,1:klev0-1)=temp(i,j,klev0)
                uvel(i,j,1:klev0-1)=uvel(i,j,klev0)
                vvel(i,j,1:klev0-1)=vvel(i,j,klev0)

                !Check
                do k=1,ilen
                  if(salt(i,j,k)<saltmin.or.salt(i,j,k)>saltmax.or. &
                    &temp(i,j,k)<tempmin.or.temp(i,j,k)>tempmax.or. &
                    &abs(uvel(i,j,k))>vmag_max.or.abs(vvel(i,j,k))>vmag_max.or. &
                    &abs(ssh(i,j))>ssh_max) then
                    write(11,*)'no valid values after searching(1):',it2,i,j,k,salt(i,j,:), &
     &temp(i,j,:),uvel(i,j,:),vvel(i,j,:),'; SSH:',ssh(i,j)
!'                   salt(i,j,k)=salt0(i,j,k)
!                    temp(i,j,k)=temp0(i,j,k)
!                    uvel(i,j,k)=uvel0(i,j,k)
!                    vvel(i,j,k)=vvel0(i,j,k)
!                    ssh(i,j)=ssh0(i,j)
                    stop
                  endif
                enddo !k
              endif
            enddo !j
          enddo !i

          icount=0
          do i=1,ixlen
            do j=1,iylen
              if(kbp(i,j)==-1) then !dry pt
                icount=icount+1
                !Use extended parents
                i1=iparen_of_dry(1,icount); j1=iparen_of_dry(2,icount)
                salt(i,j,1:ilen)=salt(i1,j1,1:ilen)
                temp(i,j,1:ilen)=temp(i1,j1,1:ilen)
                uvel(i,j,1:ilen)=uvel(i1,j1,1:ilen)
                vvel(i,j,1:ilen)=vvel(i1,j1,1:ilen)
                ssh(i,j)=ssh(i1,j1)

                !Check
                do k=1,ilen
                  if(salt(i,j,k)<saltmin.or.salt(i,j,k)>saltmax.or. &
                    &temp(i,j,k)<tempmin.or.temp(i,j,k)>tempmax.or. &
                    &abs(uvel(i,j,k))>vmag_max.or.abs(vvel(i,j,k))>vmag_max.or. &
                    &abs(ssh(i,j))>ssh_max) then
                    write(11,*)'Warning: no valid values after searching(2):',it2,i,j,k,salt(i,j,k), &
     &temp(i,j,k),uvel(i,j,k),vvel(i,j,k),ssh(i,j),i1,j1,'; ',salt(i1,j1,:),temp(i1,j1,:),uvel(i1,j1,:)
!'
                    stop
                  endif
                enddo !k
              endif !kbp
            enddo !j
          enddo !i
!          if(icount/=ndrypt) stop 'mismatch(8)'

          !Save for abnormal cases later
!          salt0=salt
!          temp0=temp
!          ssh0=ssh
!          uvel0=uvel
!          vvel0=vvel

          write(20,*)'done prep for step:',it2
          call flush(20)
        endif !ifile==1.and.it2==
    
!       Do interpolation
        tempout=-99; saltout=-99
        do ii=1,nond0
          i=iond2(ii)
        
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

              !Impose bounds for odd cases
              lev2=lev+1
              lev=max(1,min(ilen,lev))
              lev2=max(1,min(ilen,lev2))

              wild2(1,1)=temp(ix,iy,lev)*(1-vrat)+temp(ix,iy,lev2)*vrat
              wild2(1,2)=salt(ix,iy,lev)*(1-vrat)+salt(ix,iy,lev2)*vrat
              wild2(2,1)=temp(ix+1,iy,lev)*(1-vrat)+temp(ix+1,iy,lev2)*vrat
              wild2(2,2)=salt(ix+1,iy,lev)*(1-vrat)+salt(ix+1,iy,lev2)*vrat
              wild2(3,1)=temp(ix+1,iy+1,lev)*(1-vrat)+temp(ix+1,iy+1,lev2)*vrat
              wild2(3,2)=salt(ix+1,iy+1,lev)*(1-vrat)+salt(ix+1,iy+1,lev2)*vrat
              wild2(4,1)=temp(ix,iy+1,lev)*(1-vrat)+temp(ix,iy+1,lev2)*vrat
              wild2(4,2)=salt(ix,iy+1,lev)*(1-vrat)+salt(ix,iy+1,lev2)*vrat

              wild2(5,1)=uvel(ix,iy,lev)*(1-vrat)+uvel(ix,iy,lev2)*vrat
              wild2(5,2)=vvel(ix,iy,lev)*(1-vrat)+vvel(ix,iy,lev2)*vrat
              wild2(6,1)=uvel(ix+1,iy,lev)*(1-vrat)+uvel(ix+1,iy,lev2)*vrat
              wild2(6,2)=vvel(ix+1,iy,lev)*(1-vrat)+vvel(ix+1,iy,lev2)*vrat
              wild2(7,1)=uvel(ix+1,iy+1,lev)*(1-vrat)+uvel(ix+1,iy+1,lev2)*vrat
              wild2(7,2)=vvel(ix+1,iy+1,lev)*(1-vrat)+vvel(ix+1,iy+1,lev2)*vrat
              wild2(8,1)=uvel(ix,iy+1,lev)*(1-vrat)+uvel(ix,iy+1,lev2)*vrat
              wild2(8,2)=vvel(ix,iy+1,lev)*(1-vrat)+vvel(ix,iy+1,lev2)*vrat

              tempout(k,i)=dot_product(wild2(1:4,1),arco(1:4,i))
              saltout(k,i)=dot_product(wild2(1:4,2),arco(1:4,i))
              uout(k,i)=dot_product(wild2(5:8,1),arco(1:4,i))
              vout(k,i)=dot_product(wild2(5:8,2),arco(1:4,i))
              eout(i)=ssh(ix,iy)*arco(1,i)+ssh(ix+1,iy)*arco(2,i)+ssh(ix+1,iy+1)*arco(3,i)+ssh(ix,iy+1)*arco(4,i)


!              if(in==1) then
!                tempout(k,i)=wild2(1,1)*arco(1,i)+wild2(2,1)*arco(2,i)+wild2(3,1)*arco(3,i)
!                saltout(k,i)=wild2(1,2)*arco(1,i)+wild2(2,2)*arco(2,i)+wild2(3,2)*arco(3,i)
!                uout(k,i)=wild2(5,1)*arco(1,i)+wild2(6,1)*arco(2,i)+wild2(7,1)*arco(3,i)
!                vout(k,i)=wild2(5,2)*arco(1,i)+wild2(6,2)*arco(2,i)+wild2(7,2)*arco(3,i)
!                eout(i)=ssh(ix,iy)*arco(1,i)+ssh(ix+1,iy)*arco(2,i)+ssh(ix+1,iy+1)*arco(3,i)
!              else
!                tempout(k,i)=wild2(1,1)*arco(1,i)+wild2(3,1)*arco(2,i)+wild2(4,1)*arco(3,i)
!                saltout(k,i)=wild2(1,2)*arco(1,i)+wild2(3,2)*arco(2,i)+wild2(4,2)*arco(3,i)
!                uout(k,i)=wild2(5,1)*arco(1,i)+wild2(7,1)*arco(2,i)+wild2(8,1)*arco(3,i)
!                vout(k,i)=wild2(5,2)*arco(1,i)+wild2(7,2)*arco(2,i)+wild2(8,2)*arco(3,i)
!                eout(i)=ssh(ix,iy)*arco(1,i)+ssh(ix+1,iy+1)*arco(2,i)+ssh(ix,iy+1)*arco(3,i)
!              endif

              !Check
              if(tempout(k,i)<tempmin.or.tempout(k,i)>tempmax.or. &
                saltout(k,i)<saltmin.or.saltout(k,i)>saltmax.or. &
                abs(uout(k,i))>vmag_max.or.abs(vout(k,i))>vmag_max.or.abs(eout(i))>ssh_max) then
                write(11,*)'Interpolated values invalid:',i,k,tempout(k,i),saltout(k,i),uout(k,i), &
     &vout(k,i),eout(i),temp(ix,iy,lev)
                stop
              endif

!             Enforce lower bound for temp. for eqstate
              tempout(k,i)=max(0.,tempout(k,i))
              saltout(k,i)=min(40.,saltout(k,i))

!             Enforce lower bound for salt (this is the only occurence in the code)
!             if(z(k,i)<=-ht) saltout(k,i)=max(saltout(k,i),smin)
            enddo !k=1,nvrt

            !Estuary pts
!            if(iest(i)==1) then
!              tempout(:,i)=tem_es
!              saltout(:,i)=sal_es
!            endif
          endif !ixy(i,1)==0.or.
        enddo !ii

!       *[23D].th
        aa1(1)=timeout
        iret=nf90_put_var(ncids(5),itime_id(1),aa1,(/irecout/),(/1/))
        iret=nf90_put_var(ncids(6),itime_id(2),aa1,(/irecout/),(/1/))
        iret=nf90_put_var(ncids(7),itime_id(3),aa1,(/irecout/),(/1/))
        iret=nf90_put_var(ncids(8),itime_id(4),aa1,(/irecout/),(/1/))

        iret=nf90_put_var(ncids(5),ivarid(51),eout(iond2(1:nond0)),(/1,1,1,irecout/),(/1,1,nond0,1/))
        iret=nf90_put_var(ncids(6),ivarid(52),uout(1:nvrt,iond2(1:nond0)), &
     &(/1,1,1,irecout/),(/1,nvrt,nond0,1/))
        iret=nf90_put_var(ncids(6),ivarid(52),vout(1:nvrt,iond2(1:nond0)), &
     &(/2,1,1,irecout/),(/1,nvrt,nond0,1/))
        iret=nf90_put_var(ncids(7),ivarid(53),tempout(1:nvrt,iond2(1:nond0)), &
     &(/1,1,1,irecout/),(/1,nvrt,nond0,1/))
        iret=nf90_put_var(ncids(8),ivarid(54),saltout(1:nvrt,iond2(1:nond0)), &
     &(/1,1,1,irecout/),(/1,nvrt,nond0,1/))

        irecout=irecout+1
        timeout=timeout+dtout !sec

!        write(20,*)'done time:',timeout/86400
!        call flush(20)

        if(timeout/86400>nndays) stop 'finished'

      enddo !it2=ilo,ntime

      status = nf90_close(sid)
      status = nf90_close(ncids(1))
      status = nf90_close(ncids(2))
!end11
      print*, 'done reading nc for file ',ifile
!--------------------------------------------------------------------
      enddo !ifile

      iret=nf90_close(ncids(5))
      iret=nf90_close(ncids(6))
      iret=nf90_close(ncids(7))
      iret=nf90_close(ncids(8))

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

      end program gen_bc

