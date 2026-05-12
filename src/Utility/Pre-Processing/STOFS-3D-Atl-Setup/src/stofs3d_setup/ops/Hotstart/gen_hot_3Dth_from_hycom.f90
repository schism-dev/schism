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

! Generate hotstart.nc, and *[23]D.th.nc from gridded HYCOM data (nc file)
! The section on nc read needs to be modified as appropriate- search for
! 'start11' and 'end11'
! Beware the order of vertical levels in the nc file!!!
! Assume elev=0 in interpolation of 3D variables
! The code will do all it can to extrapolate: below bottom/above surface.
! If a pt in hgrid.ll is outside the background nc grid, const. values will be filled (from gen_hot_3Dth_from_nc.in) 
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
!     (5) gen_hot_3Dth_from_nc.in: 
!                     1st line: 1: include vel and elev. in hotstart.nc (and *[23D].th will start from non-0 values); 0: only T,S
!                     2nd line: T,S values for estuary points defined in estuary.gr3
!                     3rd line: T,S values for pts outside bg grid in nc
!                     4th line: time step in .nc in sec
!                     5th line: nob, iob(1:nob) - # of open bnd seg's that need
!                               *[23D].th; list of seg IDs. All *[23D].th must share same set of bnd seg's
!                     6th line: # of days needed
!                     7th line: # of HYCOM file stacks
!     (6) HYCOM files: [SSH,TS,UV]_[1,2,..nfiles].nc (beware scaling etc)
!   Output: hotstart.nc; *[23D].th.nc
!   Debug outputs: fort.11 (fatal errors); fort.20 (warning); fort.2[1-9], fort.9[5-9], fort.100; backup.out
!   Use note: if the domain is large and you wish to download HYCOM
!   fast, you may want to run this script twice: first to generate
!   hotstart only (so you need only to download a few HYCOM files); the
!   .th.nc are junk in this step.
!   Second, download HYCOM to only cover the open boundary segments to
!   generate .th.nc (hotstart would be junk).

! ifort -O2 -mcmodel=medium -assume byterecl -CB -o gen_hot_3Dth_from_hycom.exe ../UtilLib/schism_geometry.f90 \
! ../UtilLib/extract_mod.f90 ../UtilLib/compute_zcor.f90 ../UtilLib/pt_in_poly_test.f90 gen_hot_3Dth_from_hycom.f90 \
!-I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff

      program gen_hot
      use netcdf
      use schism_geometry_mod
      use compute_zcor
      use pt_in_poly_test, only: signa_single

!      implicit none

!      integer, parameter :: debug=1
!      integer, parameter :: mnp=480000*6
!      integer, parameter :: mne=960000*6
!      integer, parameter :: mns=1440000*6
!      integer, parameter :: mnv=50
!      integer, parameter :: mnope=10 !max # of open bnd segments
!      integer, parameter :: mnond=5000 !max # of open bnd nodes in each segment
!      integer, parameter :: mnei=20 !neighbor
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
      !integer :: indel(mnei,mnp),idry_s(mns)
      !integer :: idry_s(mns)
      real, allocatable :: xl(:),yl(:),dp(:),ztot(:),sigma(:),arco(:,:), &
     &tempout(:,:),saltout(:,:),uout(:,:),vout(:,:),eout(:),su2(:,:),sv2(:,:), &
     &eout_tmp(:),tsel(:,:,:),zeros(:,:),xcj(:),ycj(:)
!      dimension ztot(0:mnv),sigma(mnv),cs(mnv),arco(4,mnp)
!      dimension tempout(mnv,mnp), saltout(mnv,mnp) !,month_day(12)
!      dimension uout(mnv,mnp),vout(mnv,mnp),eout(mnp),su2(mns,mnv),sv2(mns,mnv),eout_tmp(mnp)
!      dimension tsel(2,mnv,mne),zeros(mnv,max(mne,mnp))
!      dimension xcj(mns),ycj(mns) !,nond(mnope),iond(mnope,mnond),iob(mnope),iond2(mnope*mnond)
      allocatable :: z(:,:),sigma_lcl(:,:),kbp2(:),iparen_of_dry(:,:)
      real*8 :: aa1(1)

!     First statement
!     Currently we assume rectangular grid in HYCOM
!     (interp_mode=0). interp_mode=1 uses generic UG search (splitting
!     quads) and is kept for more generic cases
      interp_mode=0

      open(10,file='gen_hot_3Dth_from_nc.in',status='old')
      read(10,*) iuv !1: include vel and elev. in hotstart.in; 0: only T,S
      read(10,*) tem_es,sal_es !T,S values for estuary points defined in estuary.gr3
      read(10,*) tem_outside,sal_outside !T,S values for pts outside bg grid in nc
      read(10,*) dtout !time step in .nc [sec]
      read(10,*) nob !,iob(1:nob) !# of open bnds that need *3D.th; list of IDs
      read(10,*) nndays !# of days needed in output
      read(10,*) nfiles !# of stacks of HYCOM files

      allocate(iob(nob))
      rewind(10)
      do i=1,4; read(10,*); enddo
      read(10,*) nob,iob(1:nob)
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
!     Open bnds
      read(14,*) nope
      read(14,*) neta
      ntot=0
      allocate(nond(max(1,nope)))
!      if(nope>mnope) stop 'Increase mnope (2)'
      do k=1,nope
        read(14,*) nond(k)
        do i=1,nond(k)
          read(14,*) !iond(k,i)
        enddo
      enddo

      rewind(14)
      mnope=maxval(nond)
      if(mnope<=0) stop 'need open bnd'
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

!      do k=3,4
!        do i=1,k
!          do j=1,k-1
!            nx(k,i,j)=i+j
!            if(nx(k,i,j)>k) nx(k,i,j)=nx(k,i,j)-k
!            if(nx(k,i,j)<1.or.nx(k,i,j)>k) then
!              write(*,*)'nx wrong',i,j,k,nx(k,i,j)
!              stop
!            endif
!          enddo !j
!        enddo !i
!      enddo !k
!
!      do i=1,np
!        nne(i)=0
!      enddo
!
!      do i=1,ne
!        do j=1,i34(i)
!          nd=elnode(j,i)
!          nne(nd)=nne(nd)+1
!          if(nne(nd)>mnei) then
!            write(11,*)'Too many neighbors',nd
!            stop
!          endif
!          indel(nne(nd),nd)=i
!        enddo
!      enddo
!
!!     Compute ball info; this won't be affected by re-arrangement below
!      do i=1,ne
!        do j=1,i34(i)
!          ic3(j,i)=0 !index for bnd sides
!          nd1=elnode(nx(i34(i),j,1),i)
!          nd2=elnode(nx(i34(i),j,2),i)
!          do k=1,nne(nd1)
!            ie=indel(k,nd1)
!            if(ie/=i.and.(elnode(1,ie)==nd2.or.elnode(2,ie)==nd2.or.elnode(3,ie)==nd2.or.(i34(ie)==4.and.elnode(4,ie)==nd2))) ic3(j,i)=ie
!          enddo !k
!        enddo !j
!      enddo !i
!
!      ns=0 !# of sides
!      do i=1,ne
!        do j=1,i34(i)
!          nd1=elnode(nx(i34(i),j,1),i)
!          nd2=elnode(nx(i34(i),j,2),i)
!          if(ic3(j,i)==0.or.i<ic3(j,i)) then !new sides
!            ns=ns+1
!            if(ns>mns) then
!              write(11,*)'Too many sides'
!              stop
!            endif
!            elside(j,i)=ns
!            isdel(1,ns)=i
!            isidenode(1,ns)=nd1
!            isidenode(2,ns)=nd2
!            xcj(ns)=(xl(nd1)+xl(nd2))/2
!            ycj(ns)=(yl(nd1)+yl(nd2))/2
!!            dps(ns)=(dp(nd1)+dp(nd2))/2
!!            distj(ns)=dsqrt((x(nd2)-x(nd1))**2+(y(nd2)-y(nd1))**2)
!!            if(distj(ns)==0) then
!!              write(11,*)'Zero side',ns
!!              stop
!!            endif
!!            thetan=datan2(x(nd1)-x(nd2),y(nd2)-y(nd1))
!!            snx(ns)=dcos(thetan)
!!            sny(ns)=dsin(thetan)
!
!            isdel(2,ns)=ic3(j,i) !bnd element => bnd side
!!           Corresponding side in element ic3(j,i)
!            if(ic3(j,i)/=0) then !old internal side
!              iel=ic3(j,i)
!              index=0
!              do k=1,i34(iel)
!                if(ic3(k,iel)==i) then
!                  index=k
!                  exit
!                endif
!              enddo !k
!              if(index==0) then
!                write(11,*)'Wrong ball info',i,j
!                stop
!              endif
!              elside(index,iel)=ns
!            endif !ic3(j,i).ne.0
!          endif !ic3(j,i)==0.or.i<ic3(j,i)
!        enddo !j=1,i34
!      enddo !i=1,ne
!
      if(ns<ne.or.ns<np) then
        write(11,*)'Weird grid with ns < ne or ns < np',np,ne,ns
        stop
      endif
     print*, 'Done computing geometry...'

!     open(54,file='elev2D.th',access='direct',recl=nbyte*(1+nond0),status='replace')
!     open(55,file='uv3D.th',access='direct',recl=nbyte*(1+nond0*nvrt*2),status='replace')
!     open(56,file='TEM_3D.th',access='direct',recl=nbyte*(1+nond0*nvrt),status='replace')
!     open(57,file='SAL_3D.th',access='direct',recl=nbyte*(1+nond0*nvrt),status='replace')

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
      tempmax=40
      saltmin=0
      saltmax=45
      vmag_max=10 !max. |u| or |v|
      ssh_max=5 !max. of |SSH|

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
        allocate(temp(ixlen,iylen,ilen),stat=ier)
        allocate(ssh(ixlen,iylen),stat=ier)
!        allocate(uvel0(ixlen,iylen,ilen),stat=ier)
!        allocate(vvel0(ixlen,iylen,ilen),stat=ier)
!        allocate(ssh0(ixlen,iylen),stat=ier)
!        allocate(salt0(ixlen,iylen,ilen),stat=ier)
!        allocate(temp0(ixlen,iylen,ilen),stat=ier)
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
        !Define junk value for sid; the test is salt<rjunk+0.1
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
     
!         Find parent elements for hgrid.ll
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
                      if (ssh(i+ii,j+jj)>rjunk) then
                        ssh(i,j)=ssh(i+ii,j+jj);exit
                      end if
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
                  write(11,*)'Impossible (4):',i,j,klev !salt(i,j,ilen)
                  stop
                endif !klev0

!               Fill junk uvel, vvel with zero to easy eliminate HYCOM nc error
                do k=klev0,ilen
                  if (uvel(i,j,k)<rjunk) then
                    write(20,*) 'Warn! Uvel Junk in the middle:',it2,i,j,k
                    uvel(i,j,k)=0.
                  end if
                  if (vvel(i,j,k)<rjunk) then
                    write(20,*) 'Warn! Vvel Junk in the middle:',it2,i,j,k
                    vvel(i,j,k)=0.
                  end if
                end do
!               Fill junk salt ,temp with bottom value
                do k=klev0,ilen
                  if (temp(i,j,k)<rjunk) then
                    if(k==klev0) then
                      write(11,*)'Bottom T is junk(2):',it2,i,j,k,temp(i,j,k)
                      stop
                    endif
                    write(20,*) 'Warn! Temp Junk in the middle:',it2,i,j,k
                    temp(i,j,k)=temp(i,j,klev0)
                  end if
                  if (salt(i,j,k)<rjunk) then
                    if(k==klev0) then
                      write(11,*)'Bottom S is junk(2):',it2,i,j,k,salt(i,j,k)
                      stop
                    endif

                    write(20,*) 'Warn! Salt Junk in the middle:',it2,i,j,k
                    salt(i,j,k)=salt(i,j,klev0)
                  end if
                end do

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
          if(icount/=ndrypt) stop 'mismatch(8)'

          !Save for abnormal cases later
!          salt0=salt
!          temp0=temp
!          ssh0=ssh
!          uvel0=uvel
!          vvel0=vvel

          write(20,*)'done prep for step:',it2
          call flush(20)
        endif !ifile==1.and.it2==
    
!       Do interpolation: all at 1st  step, bnd only for the rest
        if(ifile==1.and.it2==ilo) then
          iup=np
        else
          iup=nond0
        endif

        tempout=-99; saltout=-99
        do ii=1,iup
          if(ifile==1.and.it2==ilo) then
            i=ii
          else
            i=iond2(ii)
          endif
        
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
            if(iest(i)==1) then
              tempout(:,i)=tem_es
              saltout(:,i)=sal_es
            endif
          endif !ixy(i,1)==0.or.
        enddo !ii

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
            write(29,*)'T profile at node ',i,dp(i)
            do k=1,nvrt
              write(29,*)k,z(k,i),tempout(k,i)
            enddo !k
          enddo !i
          do i=1,ne
            write(22,*)i,tsel(1,nvrt,i) !T
            write(23,*)i,tsel(1,1,i)
          enddo !i

!         Output hotstart 
!          open(36,file='hotstart.in',form='unformatted',status='replace')
!          write(36) 0.d0,0,1
!          do i=1,ne
!            write(36) i,0,(0.d0,dble(tsel(1:2,j,i)),j=1,nvrt)
!          enddo !i
!          do i=1,ns
!            write(36) i,0,(dble(su2(j,i)),dble(sv2(j,i)),dble(tsd(i,j)),dble(ssd(i,j)),j=1,nvrt)
!          enddo !i
!          do i=1,np
!            write(36) i,dble(eout_tmp(i)),0,(dble(tempout(j,i)),dble(saltout(j,i)), &
!                      dble(tempout(j,i)),dble(saltout(j,i)),0.d0,0.d0, &
!                      0.d0,0.d0,0.d0,0.d0,0.d0,j=1,nvrt)
!          enddo !i
!          close(36)

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
          iret=nf90_def_var(ncids(4),'nsteps_from_cold',NF90_INT,var1d_dims,ivarid(20))

          var1d_dims(1)=nele_dim
          iret=nf90_def_var(ncids(4),'idry_e',NF90_INT,var1d_dims,ivarid(4))
          var1d_dims(1)=nedge_dim
          iret=nf90_def_var(ncids(4),'idry_s',NF90_INT,var1d_dims,ivarid(5))
          var1d_dims(1)=node_dim
          iret=nf90_def_var(ncids(4),'idry',NF90_INT,var1d_dims,ivarid(6))
          iret=nf90_def_var(ncids(4),'eta2',NF90_DOUBLE,var1d_dims,ivarid(7))
          iret=nf90_def_var(ncids(4),'cumsum_eta',NF90_DOUBLE,var1d_dims,ivarid(21))

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
          iret=nf90_put_var(ncids(4),ivarid(20),0)

          iret=nf90_put_var(ncids(4),ivarid(4),idry_s(1:ne),(/1/),(/ne/))
          iret=nf90_put_var(ncids(4),ivarid(5),idry_s(1:ns),(/1/),(/ns/))
          iret=nf90_put_var(ncids(4),ivarid(6),idry_s(1:np),(/1/),(/np/))
          iret=nf90_put_var(ncids(4),ivarid(7),dble(eout_tmp(1:np)),(/1/),(/np/))
          iret=nf90_put_var(ncids(4),ivarid(21),dble(eout_tmp(1:np)),(/1/),(/np/)) !cumsum

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

!          stop

        endif !ifile; hotstart

!       *[23D].th
!        write(54,rec=irecout)timeout,eout(iond2(1:nond0))
!        write(55,rec=irecout)timeout,((uout(k,iond2(j)),vout(k,iond2(j)),k=1,nvrt),j=1,nond0)
!        write(56,rec=irecout)timeout,((tempout(k,iond2(j)),k=1,nvrt),j=1,nond0)
!        write(57,rec=irecout)timeout,((saltout(k,iond2(j)),k=1,nvrt),j=1,nond0)

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
      enddo !ifile=1

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

      end program gen_hot

