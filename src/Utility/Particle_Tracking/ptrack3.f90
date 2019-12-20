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

!										
!	SCHISM Particle tracking code for nc outputs. Works for mixed
!	tri/quads, but may not work if the bottom is moving (as the fill values in
!	nc outputs are based on initial bottom). Euler tracking only.
!	Routines adpated from SCHISM:						
!	cpp, quicksearch, intersect2, signa, area_coord, levels			
!       Warning: indices in 2D arrays are not swapped (i.e., (np,nv)).	
!                 Also interpolation is along S-coord in tracking; no
!                 bilinear interp for quads.
!                 Assume the quads are not split in the nc outputs.
!
!       Behavior when particles hit horizontal bnd or dry interface:
!       reflect off like LTRAN
!										
!	Inputs: 
!          a) hgrid.ll (if ics=2 in particle.bp) or hgrid.gr3 (if ics=1 in
!          particle.bp), with b.c. part
!          b) vgrid.in
!          c) particle.bp (see format below; also sample in this dir) 
!          d) schout*.nc (must have 3D vel, elevation; also diffusivity and wind for oil spill)
!          e) (for advanced users of oil spill etc) some constants specified near the
!          beginning of the main routine (e.g. 'di' near line 262)

!	Input particle.bp:							
!	  (1) description;						
!	  (2) nscreen;								
!         (3) mod_part: model #. 0-passive; 1: oil spill (Dr. Jung)
!	  (4) ibf: forward (=1) or backward (=-1) tracking.			
!         (5) istiff: stiff (fixed distance frm f.s.; 1) or not (0).		
!	  (6) ics,slam0,sfea0: coordinate system and center for CPP projection
!             (same as in param.in). If ics=2, inputs/outputs are in
!             lon/lat;
!	  (7) h0,rnday,dtm,nspool,ihfskip,ndeltp: min. depth, # of days, 	
!	      time step used in the original run, nspool and ihfskip used 	
!             in the run (see param.in); # of sub-division in the tracking;
!             !!!WARNING: ihfskip must be equal to the actual # of steps 	
!                         contained in a file, in the case there is only 1 file.
!             Also: rnday may be smaller than the original run.
!	  (8) nparticle: # of particles;					
!	  (9) idp(i),st_p(i),xpar(i),ypar(i),zpar0(i): particle id, start time (sec),
!		starting x,y, and z relative to the instant f.s. (<=0).		
!         (10) additional parameters for oil spill
!										
!	Output: particle.pth, particle.pth.more (more info), fort.11 (fatal errors).		
!										
! ifort -mcmodel=medium -assume byterecl -O2 -o ptrack3.exe ../UtilLib/compute_zcor.f90 ../UtilLib/schism_geometry.f90 ptrack3.f90 -I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff
! pgf90 -mcmodel=medium -O2 -o ptrack3.exe ../UtilLib/compute_zcor.f90 ../UtilLib/schism_geometry.f90 ptrack3.f90 -I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff
! gfortran -O2 -ffree-line-length-none  -o ptrack3.exe ../UtilLib/compute_zcor.f90 ../UtilLib/schism_geometry.f90 ptrack3.f90 -I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff

!...  Data type consts
      module kind_par
        implicit none
        integer, parameter :: sng_kind1=4
        integer, parameter :: dbl_kind1=8
        real(kind=dbl_kind1), parameter :: small1=1.e-6 !small non-negative number; must be identical to that in global
      end module kind_par


!...  definition of variables
!...
      module global
        implicit none
        public

        integer, parameter :: sng_kind=4
        integer, parameter :: dbl_kind=8

!...  	Dimensioning parameters
        integer, parameter :: nbyte=4
      	real(kind=dbl_kind), parameter :: zero=1.e-5 !small postive number in lieu of 0; usually used to check areas 
        real(kind=dbl_kind), parameter :: small1=1.e-6 !small non-negative number; must be identical to that in kind_par
        real(kind=dbl_kind), parameter :: pi=3.1415926d0 

!...  	Important variables
        integer, save :: np,ne,ns,nvrt,mnei,mod_part,ibf,istiff,ivcor,kz,nsig
      	real(kind=dbl_kind), save :: h0,rho0,dt
        real(kind=dbl_kind), save :: h_c,theta_b,theta_f,h_s !s_con1

!...    Output handles
        character(len=48), save :: start_time,version
        character(len=12), save :: ifile_char
        integer,save :: ihot,nrec,nspool,igmp,noutgm,ifile,noutput,ifort12(100)
        
        integer,save, allocatable :: nne(:),indel(:,:),idry(:),idry_e(:),idry_e0(:)
        integer,save, allocatable :: kbp(:),kbs(:),kbe(:),kbp00(:),isbnd(:)

        real(kind=dbl_kind),save, allocatable :: x(:),y(:),dp(:),eta1(:),eta2(:),eta3(:) !,hmod(:)
        real(kind=dbl_kind),save, allocatable :: area(:),xctr(:),yctr(:)
        real(kind=dbl_kind),save, allocatable :: snx(:),sny(:),distj(:),dps(:),dldxy(:,:,:)

        real(kind=dbl_kind),save, allocatable :: zpar0(:)
        !For interface with util routines
        real(kind=dbl_kind),save, allocatable :: ztot(:),sigma(:),xcj(:),ycj(:),sigma_lcl(:,:)

        integer,save, allocatable :: i34(:),elnode(:,:),ic3(:,:),elside(:,:),isdel(:,:),isidenode(:,:)
        integer,save, allocatable :: icum1(:,:),icum2(:,:,:)
        integer,save :: nxq(3,4,4),nodel(3)

        real(kind=dbl_kind),save, allocatable :: z(:,:)
        real(kind=dbl_kind),save, allocatable :: uu1(:,:),vv1(:,:),ww1(:,:),uu2(:,:),vv2(:,:),ww2(:,:)
        real*8,save, allocatable :: wnx1(:),wnx2(:),wny1(:),wny2(:),hf1(:,:),vf1(:,:),hf2(:,:),vf2(:,:)
        real*8,save, allocatable :: hvis_e(:,:)
      end module global

!...  Main program
      program ptrack
      use global
      use netcdf
      use compute_zcor
      use schism_geometry_mod

      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)

      character(len=30) :: file63
      real(kind=sng_kind) :: floatout,floatout2
      real*8, allocatable :: xpar(:),ypar(:),zpar(:),st_p(:),upar(:),vpar(:),wpar(:)
      real*8, allocatable :: ztmp(:),ztmp2(:),dhfx(:),dhfy(:),dhfz(:),grdx(:),grdy(:), &
     &grdz(:),amas(:),wndx(:),wndy(:),timeout(:)
      integer, allocatable :: ielpar(:),levpar(:),iabnorm(:),ist(:),inbr(:)
      character(len=25), allocatable :: idp(:)
      real*8 :: vxl(4,2),vyl(4,2),vzl(4,2),vxn(4),vyn(4),vzn(4),arco(3), &
     &dx(10),dy(10),dz(10),val(4,2),vbl(4,2),vcl(4,2),vdl(4,2),van(4),vcn(4),vdn(4),vwx(4),vwy(4)
      real, allocatable :: real_ar(:,:)
      integer :: nodel2(3)
      integer :: varid1,varid2,dimids(3),istat,nvtx,iret

      !Random seed used only for oil spill model
      iseed=5
      !Ekman effects
      rotate_angle=3.0*(pi/180.0)   !angle between wind and current directions
      drag_c=0.03 !reducing wind speed to get surface current
      !Evaporation const.
      T_half=86400 ! half-life [sec]
      remain_ratio=0.6  ! remain ratio after a long time

      !Cyclical index
      do k=3,4 !elem. type
        do i=1,k  !local index
          do j=1,k-1 !offset
            nxq(j,i,k)=i+j
            if(nxq(j,i,k)>k) nxq(j,i,k)=nxq(j,i,k)-k
            if(nxq(j,i,k)<1.or.nxq(j,i,k)>k) then
              write(*,*)'nx wrong',i,j,k,nxq(j,i,k)
              stop
            endif
          enddo !j
        enddo !i
      enddo !k

      ifort12=0 !init
      open(11,file='fort.11',status='replace')
 
!...  Read in particles
      open(95,file='particle.bp',status='old')
      read(95,*)
      read(95,*) nscreen
!     Model #: 0-passive; 1:oil spill
      read(95,*) mod_part
      if(mod_part<0.or.mod_part>1) stop 'Unknown model'
      read(95,*) ibf
      if(iabs(ibf)/=1) then
        write(*,*)'Wrong ibf',ibf
        stop
      endif
      if(mod_part==1.and.ibf/=1) stop 'Oil spill must have ibf=1'

      read(95,*) istiff !1: fixed distance from F.S.
      if(istiff/=0.and.istiff/=1) then
        write(*,*)'Wrong istiff',istiff
        stop
      endif
      read(95,*) ics,slam0,sfea0
      slam0=slam0/180*pi
      sfea0=sfea0/180*pi
      read(95,*) h0,rnday,dtm,nspool,ihfskip,ndeltp !# of sub-divisions
      if(mod(ihfskip,nspool)/=0) then
        write(*,*)'ihfskip must be a multiple of nspool'
        stop
      endif
      nrec=ihfskip/nspool !# of records (steps) per stack

      read(95,*) nparticle
      allocate(zpar0(nparticle),xpar(nparticle),ypar(nparticle),zpar(nparticle), &
     &st_p(nparticle),idp(nparticle),ielpar(nparticle),levpar(nparticle),upar(nparticle), &
     &vpar(nparticle),wpar(nparticle),iabnorm(nparticle),ist(nparticle),inbr(nparticle), &
     &dhfx(nparticle),dhfy(nparticle),dhfz(nparticle),grdx(nparticle),grdy(nparticle), &
     &grdz(nparticle),amas(nparticle),wndx(nparticle),wndy(nparticle),stat=istat)
      if(istat/=0) stop 'Failed to alloc (1)'

      levpar=-99 !vertical level
      iabnorm=0 !abnormal tracking exit flag

      dt=dtm*nspool !output time step
      st_m=(ibf+1)/2*rnday*86400 !initialize min. or max. starting time
      do i=1,nparticle
        if(ics.eq.1) then
!	  zpar0: relative to f.s.
          read(95,*)idp(i),st_p(i),xpar(i),ypar(i),zpar0(i)
        else
          read(95,*)idp(i),st_p(i),xparl,yparl,zpar0(i)
          xparl=xparl/180*pi
          yparl=yparl/180*pi
          call cpp(xpar(i),ypar(i),xparl,yparl,slam0,sfea0)
        endif
        if(st_p(i)<0.or.st_p(i)>rnday*86400) then
          write(11,*)'Starting time for particle ',i,' out of range:',st_p(i)
          stop
        endif
        if(zpar0(i)>0) then
          write(11,*)'Starting z-coord above f.s.',i
          stop
        endif
        if(ibf*st_p(i)<ibf*st_m) st_m=st_p(i)
      enddo !i

!...  Additional parameters for oil spill
! ... Description of parameters
!     ihdf  : turn on Smagorinsky algorithm - off(0), on(1)
!     ibuoy : turn buoyancy of particle  - off(0), on(1)
!     iwind : turn wind effect - off(0), on(1)
!     pbeach  : set percentage of stranding on shore
! ........................................................
      if(mod_part==1) then
        read(95,*) !comment line
        read(95,*) ihdf,hdc,horcon
        read(95,*) ibuoy,iwind
        read(95,*) pbeach
      endif !mod_part=1

      close(95)

!...  Init. ist etc; some of these are only used in certain models
!     ist(1:nparticles): 0 - inactive b4 release
!                        -2 - stranded at wet/dry interface, after some
!                        time has elapsed (only in this case is inbr (ngbring elem) defined)
!                        -1 - permanently stranded at land bnd (hit
!                        there and some time has elapsed)
!                        2 - hit the open bnd
!                        1 - active but not -1 or 2
      ist=0      ! particle status flag
      inbr=0     !neighboring elem
!      amas=1.0   ! mass of particles

! ... treatment for buoyant oil particles
      if(mod_part==1) then
        if(ibuoy==1) then
          !compute the rising velocity(m/s) based on Proctor et al., 1994
          gr=9.8               ! m/s^2
          rho_o=900.0d3        ! kg/m^3 (oil)
          rho_w=1025.0d3       ! kg/m^3
          di=500.0d-6          ! m
          smu=1.05d-6          ! m^2/s
! ... critical diameter, dc
          dc=9.52*smu**(2./3.)/(gr**(1./3.)*(1.-rho_o/rho_w)**(1./3.)) !m
!... compute the rising velocity, m/s
          if(di>=dc) then
            rsl=sqrt(8./3.*gr*di*(1.-rho_o/rho_w))  ! large droplet
          else
            rsl=gr*di**2*(1.-rho_o/rho_w)/(18.*smu) ! small droplet
          endif
        else if(ibuoy==0) then
          rsl=0.0
        else
          write(11,*)'ibuoy is incorrect'
          stop
        endif !ibuoy
        print*, 'Rising vel=',rsl
      endif !mod_part=1

!...  Read in header (if quads r split in binary, the hgrid read in here
!     has different conn table from hgrid.gr3)
      ntime=nint(rnday*86400/dt) !total # of records
      ifile=1 !for st_m=0; starting stack #
      do i=1,ntime/nrec+1
        if(st_m/dt>(i-1)*nrec.and.st_m/dt<=i*nrec) then
          ifile=i
          exit
        endif
      enddo
      print*, 'Starting stack =',ifile,st_m,ntime
      if(ibf==1) then
        iths=(ifile-1)*nrec+1 !starting _total_ record # 
      else !ibf=-1
        iths=ifile*nrec !starting _total_ record #
      endif
      print*, 'Starting step # =',iths

      write(ifile_char,'(i12)') ifile
      ifile_char=adjustl(ifile_char); len_char=len_trim(ifile_char)
      file63='schout_'//ifile_char(1:len_char)//'.nc'
      iret=nf90_open(trim(adjustl(file63)),OR(NF90_NETCDF4,NF90_NOWRITE),ncid)
      if(iret/=nf90_NoErr) stop '1st stack not opened'
      iret=nf90_inq_varid(ncid,'elev',ielev_id)
      if(iret/=nf90_NoErr) stop 'elev not found'
      iret=nf90_inq_varid(ncid,'hvel',luv)
      if(iret/=nf90_NoErr) stop 'hvel not found'
      iret=nf90_inq_varid(ncid,'vertical_velocity',lw)
      if(iret/=nf90_NoErr) stop 'w not found'
      if(mod_part==1) then
        iret=nf90_inq_varid(ncid,'wind_speed',lwind)
        if(iret/=nf90_NoErr) stop 'wind not found'
        iret=nf90_inq_varid(ncid,'diffusivity',ltdff)
        if(iret/=nf90_NoErr) stop 'tdiff not found'
      endif !mod_part

!      open(60,file=trim(ifile_char)//'_elev.61',access='direct',recl=nbyte)
!      open(61,file=trim(ifile_char)//'_hvel.64',access='direct',recl=nbyte)
!      open(62,file=trim(ifile_char)//'_vert.63',access='direct',recl=nbyte)
!      if(mod_part==1) then
!        open(63,file=trim(ifile_char)//'_wind.62',access='direct',recl=nbyte)
!        open(64,file=trim(ifile_char)//'_tdff.63',access='direct',recl=nbyte)
!      endif !mod_part

      iret=nf90_inq_dimid(ncid,'nSCHISM_vgrid_layers',i)
      if(iret/=nf90_NoErr) stop 'nSCHISM_vgrid_layers not found'
      iret=nf90_Inquire_Dimension(ncid,i,len=nvrt)
      if(iret/=nf90_NoErr) stop 'read error(1)'
      iret=nf90_inq_varid(ncid,'SCHISM_hgrid_face_nodes',varid1)
      if(iret/=nf90_NoErr) stop 'read error(2)'
      iret=nf90_Inquire_Variable(ncid,varid1,dimids=dimids(1:2))
      if(iret/=nf90_NoErr) stop 'read error(3)'
      iret=nf90_Inquire_Dimension(ncid,dimids(1),len=nvtx)
      if(iret/=nf90_NoErr) stop 'read error(4)'
      iret=nf90_Inquire_Dimension(ncid,dimids(2),len=ne)
      if(iret/=nf90_NoErr) stop 'read error(5)'
      if(nvtx/=4) stop 'vtx/=4'
      iret=nf90_inq_varid(ncid,'SCHISM_hgrid_node_x',varid2)
      if(iret/=nf90_NoErr) stop 'read error(6)'
      iret=nf90_Inquire_Variable(ncid,varid2,dimids=dimids)
      if(iret/=nf90_NoErr) stop 'read error(7)'
      iret=nf90_Inquire_Dimension(ncid,dimids(1),len=np)
      if(iret/=nf90_NoErr) stop 'read error(8)'
      iret=nf90_inq_varid(ncid,'time',itime_id)
      if(iret/=nf90_NoErr) stop 'read error(9)'
      iret=nf90_Inquire_Variable(ncid,itime_id,dimids=dimids)
      if(iret/=nf90_NoErr) stop 'read error(10)'
      iret=nf90_Inquire_Dimension(ncid,dimids(1),len=nrec)
      if(iret.ne.NF90_NOERR) then
        print*, nf90_strerror(iret) 
        stop 'readheader: error reading header'
      endif

      allocate(x(np),y(np),dp(np),kbp00(np),i34(ne),elnode(4,ne),timeout(nrec), &
     &area(ne),nne(np),idry(np),idry_e(ne),idry_e0(ne),kbp(np),kbe(ne), &
     &eta1(np),eta2(np),eta3(np),xctr(ne),yctr(ne),isbnd(np),wnx1(np), &
     &wnx2(np),wny1(np),wny2(np),dldxy(3,2,ne),real_ar(nvrt,np),stat=istat)
      if(istat/=0) stop 'failed to allocate (3)'

      iret=nf90_get_var(ncid,varid1,elnode)
      if(iret/=nf90_NoErr) stop 'read error(11)'
      iret=nf90_get_var(ncid,varid2,x)
      if(iret/=nf90_NoErr) stop 'read error(12)'
      iret=nf90_inq_varid(ncid,'SCHISM_hgrid_node_y',varid1)
      if(iret/=nf90_NoErr) stop 'read error(13)'
      iret=nf90_get_var(ncid,varid1,y)
      if(iret/=nf90_NoErr) stop 'read error(14)'
      iret=nf90_inq_varid(ncid,'depth',varid1)
      if(iret/=nf90_NoErr) stop 'read error(15)'
      iret=nf90_get_var(ncid,varid1,dp)
      if(iret/=nf90_NoErr) stop 'read error(16)'
      iret=nf90_inq_varid(ncid,'node_bottom_index',varid1)
      if(iret/=nf90_NoErr) stop 'read error(17)'
      iret=nf90_get_var(ncid,varid1,kbp00)
      if(iret/=nf90_NoErr) stop 'kbp00 not found'

      !Leave it open as this is the 1st stack to read from
!      iret=nf90_close(ncid)

      !print*, 'nc dim:',nvrt,np,ne,nrec,start_time

      !Calc i34
      i34=4 !init
      do i=1,ne
        if(elnode(4,i)<0) i34(i)=3
      enddo !i


!      irec00=5*48/nbyte+5 !+7+nvrt 
!      read(60,rec=irec00+1) nvrt
!      irec00=irec00+7+nvrt
!      read(60,rec=irec00+1) np
!      read(60,rec=irec00+2) ne
!      print*, 'nvrt,np,ne=',nvrt,np,ne

!      irec00=irec00+2
!      do m=1,np
!        !read x,y from hgrid for mroe precision
!        read(60,rec=irec00+4)kbp00(m)
!        irec00=irec00+4
!      enddo !m=1,np
!      do m=1,ne
!        read(60,rec=irec00+1)i34(m)
!        irec00=irec00+1
!        do mm=1,i34(m)
!          read(60,rec=irec00+1)elnode(mm,m)
!          irec00=irec00+1
!        enddo !mm
!      enddo !m

!      irec00=irec00+4*ne

!     Initialize kbp for levels()
!     if(ivcor==2) kbp=kbp00

!...  Read in h- and v-grid and compute geometry
!...  Since binary may split quads, do not read in conn table
      if(ics==1) then
        open(14,file='hgrid.gr3',status='old')
      else
        open(14,file='hgrid.ll', status='old')
      endif
      read(14,*) 
      read(14,*) ne2,np2
      if(np/=np2) stop 'mismatch (3)'
      do i=1,np
        if(ics==1) then
          read(14,*) j,x(i),y(i),dp(i)
        else if(ics==2) then
          read(14,*) j,xlon,ylat,dp(i)
          ylat=ylat/180*pi
          xlon=xlon/180*pi
          call cpp(x(i),y(i),xlon,ylat,slam0,sfea0)
        endif
        !hmod(i)=min(dp(i),dble(h_s))
      enddo !i=1,np

      do i=1,ne2
        !Use conn table from binary
        read(14,*) !j,i34(i),elnode(1:i34(i),i)
      enddo !i

      do i=1,ne
        n1=elnode(1,i)
        n2=elnode(2,i)
        n3=elnode(3,i)
        area(i)=signa(x(n1),x(n2),x(n3),y(n1),y(n2),y(n3))

        !Derivative of shape function for triangles only
        !For quds, use nodes 1-3
        do j=1,3
          id1=j+1
          id2=j+2
          if(id1>3) id1=id1-3
          if(id2>3) id2=id2-3
          dldxy(j,1,i)=(y(elnode(id1,i))-y(elnode(id2,i)))/2/area(i)
          dldxy(j,2,i)=(x(elnode(id2,i))-x(elnode(id1,i)))/2/area(i)
        enddo !j

        if(i34(i)==4) then
          n4=elnode(4,i)
          area(i)=area(i)+signa(x(n1),x(n3),x(n4),y(n1),y(n3),y(n4))
        endif

        if(area(i)<=0) then
          write(11,*)'Negative area at',i
          stop
        endif
      enddo !i=1,ne

!     Open bnds (not affecetd by splitting of quads)
      isbnd=0
      read(14,*) nope
      read(14,*) neta
      ntot=0
      do k=1,nope
        read(14,*) nond
        do i=1,nond
          read(14,*) iond
          isbnd(iond)=k
        enddo
        ntot=ntot+nond
      enddo

      if(neta/=ntot) then
        write(11,*)'neta /= total # of open bnd nodes',neta,ntot
        stop
      endif

!     Land bnds
      read(14,*) nland
      read(14,*) nvel
      do k=1,nland
        read(14,*) nlnd
        do i=1,nlnd
          read(14,*) ilnd
          if(isbnd(ilnd)==0) isbnd(ilnd)=-1 !overlap of open bnd
        enddo
      enddo !k=1,nland

      close(14)
!...  End fort.14

!     vgrid
      open(19,file='vgrid.in',status='old')
      read(19,*)ivcor
      read(19,*)nvrt
      rewind(19)
      allocate(sigma_lcl(nvrt,np),z(np,nvrt),icum1(np,nvrt),icum2(np,nvrt,2), &
     &uu1(np,nvrt),vv1(np,nvrt),ww1(np,nvrt),uu2(np,nvrt),vv2(np,nvrt),ww2(np,nvrt), &
     &ztmp(nvrt),ztmp2(nvrt),sigma(nvrt),ztot(nvrt),hf1(np,nvrt), &
     &hf2(np,nvrt),vf1(np,nvrt),vf2(np,nvrt),hvis_e(ne,nvrt),stat=istat)
      if(istat/=0) stop 'Failed to alloc (3)'
      call get_vgrid_double('vgrid.in',np,nvrt,ivcor,kz,h_s,h_c,theta_b,theta_f,ztot,sigma,sigma_lcl,kbp)
      !kbp has been assigned for ivcor=1

!     Init some arrays (for below bottom etc)
      uu2=0; vv2=0; ww2=0; vf2=0; hf2=0; wnx2=0; wny2=0

!******************************************************************************
!                                                                             *
!     			Compute geometry 				      *
!                                                                             *
!******************************************************************************
!                                                                             *
!                                                                             *
!     We also need elem ball
      nne=0
      do i=1,ne
        do j=1,i34(i)
          nd=elnode(j,i)
          nne(nd)=nne(nd)+1
        enddo
      enddo
      mnei=maxval(nne)

      allocate(indel(mnei,np),stat=istat)
      if(istat/=0) stop 'Failed to alloc. indel'
      nne=0
      do i=1,ne
        do j=1,i34(i)
          nd=elnode(j,i)
          nne(nd)=nne(nd)+1
          if(nne(nd)>mnei) then
            write(11,*)'Too many neighbors',nd
            stop
          endif
          indel(nne(nd),nd)=i
        enddo
      enddo !i

      call compute_nside(np,ne,i34,elnode,ns)
      if(nscreen.eq.1) write(*,*) 'There are',ns,' sides in the grid...'

!     Allocate side-related arrays
      allocate(ic3(4,ne),elside(4,ne),isdel(2,ns),isidenode(2,ns),xcj(ns),ycj(ns), &
     &snx(ns),sny(ns),distj(ns),kbs(ns),dps(ns),stat=istat)
      if(istat/=0) stop 'Failed to alloc (4)'
!     Then compute the rest of side related arrays with additional
!     inputs (xnd,ynd) (x,y coordinates of each node)
      call schism_geometry_double(np,ne,ns,x,y,i34,elnode,ic3,elside,isdel,isidenode,xcj,ycj)

      !Remaining side arrays
      do i=1,ns
        nd1=isidenode(1,i)
        nd2=isidenode(2,i)
        dps(i)=(dp(nd1)+dp(nd2))/2
        distj(i)=dsqrt((x(nd2)-x(nd1))**2+(y(nd2)-y(nd1))**2)
        if(distj(i)==0) then
          write(11,*)'Zero side',i
          stop
        endif

        thetan=atan2(x(nd1)-x(nd2),y(nd2)-y(nd1))
        snx(i)=cos(thetan) !from elem 1->2
        sny(i)=sin(thetan)
      enddo !i

!...  compute elem centers 
!...
      do i=1,ne
        xctr(i)=0
        yctr(i)=0
        do j=1,i34(i)
          xctr(i)=xctr(i)+x(elnode(j,i))/i34(i)
          yctr(i)=yctr(i)+y(elnode(j,i))/i34(i)
        enddo !j
      enddo !i=1,ne

      if(nscreen.eq.1) write(*,*)'done computing geometry...'

!...  Compute record # offset for a node and level for 3D outputs
!...
!      icount1=0
!      icount2=0
!      do i=1,np
!        do k=max0(1,kbp00(i)),nvrt
!          do m=1,2
!            icount2=icount2+1
!            icum2(i,k,m)=icount2
!          enddo !m
!          icount1=icount1+1
!          icum1(i,k)=icount1
!        enddo !k
!      enddo !i=1,np

!      if(ibf==1) then
!        irec01= !starting time record # in the 1st stack
!        irec02=irec00 !hvel.64
!        irec03=irec00 !vert.63
!        irec04=irec00 !wind.62
!        irec05=irec00 !tdff.63
!      else !ibf=-1; no oil spill
!        irec01=irec00+(2+2*np)*nrec+1 !last record in ifile plus 1 (for backward reading) 
!        irec02=irec00+(2+np+icum2(np,nvrt,2))*nrec+1
!        irec03=irec00+(2+np+icum1(np,nvrt))*nrec+1
!        irec04=0 !not used
!        irec05=0 !not used
!      endif
!      irec1=irec01
!      irec2=irec02
!      irec3=irec03
!      irec4=irec04 !junk if ibf/=1
!      irec5=irec05 !junk if ibf/=1

!...  Compute initial elements for particle tracking
!...
      lp1: do i=1,nparticle
        do k=1,ne
          call pt_in_poly2(i34(k),x(elnode(1:i34(k),k)),y(elnode(1:i34(k),k)),xpar(i),ypar(i),inside)
          if(inside/=0) then
            ielpar(i)=k
            cycle lp1
          endif
        enddo !k=1,ne
        write(11,*)'Cannot find init. element for particle',i
        stop
      end do lp1 !i=1,nparticle

      open(95,file='particle.pth',status='replace')
      open(97,file='particle.pth.more',status='replace')
      write(95,*)'Drogues'
      if(ibf==1) then
        write(95,*) ntime-iths+1
      else
        write(95,*) iths
      endif

      if(nscreen.eq.1) write(*,*)'done initialization...'


!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!									!
!	       Time iteration						!
!									!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Initialize for output before particle moving
      upar=0; vpar=0; wpar=0
      zpar=zpar0

      !Starting record # in the starting stack
      if(ibf==1) then
        it2=ntime !last record (total)
        irec1=1
      else
        it2=1
        irec1=nrec
      endif
      do it=iths,it2,ibf !it is total time record #
!--------------------------------------------------------------------------
      time=it*dt
      
!...  Read in elevation and vel. info
      if((ibf==1.and.it>nrec*ifile).or.(ibf==-1.and.it<=nrec*(ifile-1))) then
        !Open next stack
        iret=nf90_close(ncid)
        ifile=ifile+ibf !stack #
        write(ifile_char,'(i12)') ifile
        ifile_char=adjustl(ifile_char); len_char=len_trim(ifile_char)
        file63='schout_'//ifile_char(1:len_char)//'.nc'
        iret=nf90_open(trim(adjustl(file63)),OR(NF90_NETCDF4,NF90_NOWRITE),ncid)
        if(iret/=nf90_NoErr) stop 'schout not opened (2)'
        !time is double
        iret=nf90_inq_varid(ncid,'time',itime_id)
        if(iret/=nf90_NoErr) stop 'itime_id'
        iret=nf90_get_var(ncid,itime_id,timeout,(/1/),(/nrec/))
        if(iret/=nf90_NoErr) stop 'time not read'

        !Reset record #
        if(ibf==1) then
          irec1=1
        else
          irec1=nrec
        endif

!        close(60)
!        close(61)
!        close(62)
!        open(60,file=trim(ifile_char)//'_elev.61',access='direct',recl=nbyte)
!        open(61,file=trim(ifile_char)//'_hvel.64',access='direct',recl=nbyte)
!        open(62,file=trim(ifile_char)//'_vert.63',access='direct',recl=nbyte)
!        irec2=irec02
!        irec3=irec03

!        if(mod_part==1) then
!          close(63)
!          close(64)
!          open(63,file=trim(ifile_char)//'_wind.62',access='direct',recl=nbyte)
!          open(64,file=trim(ifile_char)//'_tdff.63',access='direct',recl=nbyte)
!          irec4=irec04
!          irec5=irec05
!        endif !mod_part
      endif !ibf

      if(irec1<1.or.irec1>nrec) stop 'record out of bound'

      iret=nf90_get_var(ncid,ielev_id,real_ar(1,1:np),(/1,irec1/),(/np,1/))
      if(iret/=nf90_NoErr) stop 'elev not read'
      eta2=real_ar(1,1:np)
      iret=nf90_get_var(ncid,luv,real_ar(1:nvrt,1:np),(/1,1,1,irec1/),(/1,nvrt,np,1/))
      if(iret/=nf90_NoErr) stop 'uu2 not read'
      uu2(:,:)=transpose(real_ar(1:nvrt,1:np))
      iret=nf90_get_var(ncid,luv,real_ar(1:nvrt,1:np),(/2,1,1,irec1/),(/1,nvrt,np,1/))
      if(iret/=nf90_NoErr) stop 'vv2 not read'
      vv2(:,:)=transpose(real_ar(1:nvrt,1:np))
      iret=nf90_get_var(ncid,lw,real_ar(1:nvrt,1:np),(/1,1,irec1/),(/nvrt,np,1/))
      if(iret/=nf90_NoErr) stop 'ww2 not read'
      ww2(:,:)=transpose(real_ar(1:nvrt,1:np))
      if(mod_part==1) then
        iret=nf90_get_var(ncid,lwind,real_ar(1:2,1:np),(/1,1,irec1/),(/2,np,1/))
        if(iret/=nf90_NoErr) stop 'wind not read'
        wnx2(:)=real_ar(1,1:np)
        wny2(:)=real_ar(2,1:np)
        iret=nf90_get_var(ncid,ltdff,real_ar(1:nvrt,1:np),(/1,1,irec1/),(/nvrt,np,1/))
        if(iret/=nf90_NoErr) stop 'tdiff not read'
        vf2(:,:)=transpose(real_ar(1:nvrt,1:np))
      endif !mod_part

      irec1=irec1+ibf

!...  Store info for first step
      if(it==iths) eta1=eta2

!...  Compute elevation eta3
      do i=1,np
        if(ibf==1) then
          eta3(i)=eta1(i)
        else
          eta3(i)=eta2(i)
        endif
      enddo !i

!...  Compute z-cor
      call levels

!...  Deal with junks
      do i=1,np
        if(idry(i)==1) then
          uu2(i,:)=0; vv2(i,:)=0; ww2(i,:)=0; vf2(i,:)=0
        else !note that the fill values in nc are based on init bottom, which may be different from kbp
          do k=1,nvrt
            if(k<=kbp(i)-1.or.abs(uu2(i,k))>1.e8) then
              uu2(i,k)=0
              vv2(i,k)=0
              ww2(i,k)=0
              vf2(i,k)=0
            endif
          enddo !k
        endif 
      enddo !i
      
!...  Store info for first step
      if(it==iths) then
        uu1=uu2; vv1=vv2; ww1=ww2
        vf1=vf2; wnx1=wnx2; wny1=wny2
      endif 

!...  compute hvis_e & hvis_e based on Smagorinsky Algorithm
      if(mod_part==1) then
        if(ihdf==0) then
          hf2=hdc     ! constant diffusivity [m^2/s]
        else !Smag.
          hvis_e=0
          do i=1,ne
            if(idry_e(i)==1) cycle

            do k=kbe(i),nvrt
              dudx_e=dot_product(uu2(elnode(1:3,i),k),dldxy(1:3,1,i)) 
              dudy_e=dot_product(uu2(elnode(1:3,i),k),dldxy(1:3,2,i))
              dvdx_e=dot_product(vv2(elnode(1:3,i),k),dldxy(1:3,1,i))
              dvdy_e=dot_product(vv2(elnode(1:3,i),k),dldxy(1:3,2,i))
              hvis_e(i,k)=horcon*area(i)*sqrt(dudx_e**2+dvdy_e**2+0.5*(dudy_e+dvdx_e)**2) !m^2/s
            enddo !k
          enddo !i=1,ne

!...      Convert to nodes
          hf2=0
          do i=1,np
            if(idry(i)==1) cycle

            do k=kbp(i),nvrt
              spm=0;slm=0
              do j=1,nne(i)
                nel=indel(j,i)
                rl=sqrt((x(i)-xctr(nel))**2+(y(i)-yctr(nel))**2)
                if(rl==0) stop 'rl=0'
                slm=slm+1./rl
                spm=spm+hvis_e(nel,k)/rl
              enddo
              hf2(i,k)=spm/slm !slm>0
            enddo !k
          enddo !i=1,np
        endif !ihdf

!...    Store info for diffusion term in first step
        if(it==iths) then
          hf1=hf2; grdx=0.0; grdy=0.0; grdz=0.0
          dhfx=hdc; dhfy=hdc; dhfz=3.0d-4
        endif
      endif !mod_part=1

      if(nscreen.eq.1) write(*,*)'begin ptrack...'

!...  Particle tracking
      write(95,*) time,nparticle
      write(97,*) time,nparticle
      do i=1,nparticle
        eta_p=0; dp_p=0 !for output before moving
        if((ibf==1.and.time<=st_p(i)).or.(ibf==-1.and.time>st_p(i)-dt)) go to 449 !output directly

        if(mod_part==1) then
          if(ist(i)==0) ist(i)=1 !particle activated
          if(ist(i)==-2) then
            if(inbr(i)>0) then; if(idry_e(inbr(i))==0) then
              ist(i)=1 !re-activated
            endif; endif
          endif
          if(ist(i)/=1) go to 449
        endif !mod_part

        pt=dt !tracking time step
!...    Initialize starting level 
        if(levpar(i)==-99) then !just started
          pt=(time-st_p(i))*ibf
          if(pt<=0) then
            write(*,*)'Tracking step negative:',pt
            stop
          endif
          iel=ielpar(i)
          if(idry_e(iel)==1) then !dry
            levpar(i)=-1
          else !wet
            call pt_in_poly3(i34(iel),x(elnode(1:i34(iel),iel)),y(elnode(1:i34(iel),iel)), &
     &xpar(i),ypar(i),arco,nodel)

            do k=kbe(iel),nvrt
              ztmp2(k)=0
              do j=1,3
                nd=elnode(nodel(j),iel)
                ztmp2(k)=ztmp2(k)+z(nd,max(k,kbp(nd)))*arco(j)
              enddo !j
            enddo !k
            zpar(i)=max(zpar0(i)+ztmp2(nvrt),ztmp2(kbe(iel))) !zpar0<=0
            jlev=0
            do k=kbe(iel),nvrt-1
              if(zpar(i)>=ztmp2(k).and.zpar(i)<=ztmp2(k+1)) then
                jlev=k+1
                exit
              endif
            enddo !k
            if(jlev==0) then
              write(11,*)'Cannot find an init. level:',i,zpar(i),(ztmp2(k),k=kbe(iel),nvrt)
              stop
            endif
            levpar(i)=jlev
          
            upar(i)=0; vpar(i)=0; wpar(i)=0; wndx(i)=0; wndy(i)=0
            do j=1,3
              nd=elnode(nodel(j),iel)
              upar(i)=upar(i)+uu2(nd,jlev)*arco(j)
              vpar(i)=vpar(i)+vv2(nd,jlev)*arco(j)
              wpar(i)=wpar(i)+ww2(nd,jlev)*arco(j)
              if(mod_part==1) then
                wndx(i)=wndx(i)+wnx2(nd)*arco(j)
                wndy(i)=wndy(i)+wny2(nd)*arco(j)
              endif !mod_part
            enddo !j

          endif !wet
        endif !levpar=-99

!	Wetting/drying
        if(idry_e(ielpar(i))==1) then
          levpar(i)=-1
          go to 449
        endif

        nnel=ielpar(i)
        jlev=levpar(i)
!       Rewetted elements
        if(jlev==-1) then !element nnel wet
          jlev=nvrt
          zpar(i)=sum(eta3(elnode(1:i34(nnel),nnel)))/i34(nnel)
        endif
  
!	Tracking
        x0=xpar(i)
        y0=ypar(i)
        z0=zpar(i)
        nnel0=nnel
        jlev0=jlev
        dtb=pt/ndeltp
        do idt=1,ndeltp
          if(ibf==1) then
            trat=real(idt)/ndeltp
          else
            trat=real(idt-1)/ndeltp
          endif

          if(mod_part==1) then !oil spill
! ...       wind rotation & apply to surface current
            dir = atan2(wndy(i),wndx(i))
            speed = sqrt(wndx(i)**2+wndy(i)**2)
            dir = dir+rotate_angle
            if(iwind==0) then
              cur_x=0; cur_y=0
            else
              cur_x = speed*sin(dir)*drag_c !approx surface current
              cur_y = speed*cos(dir)*drag_c
            endif  
            if(z0<-0.1) then  ! sub_surface particles are not influenced by wind
              cur_x=0; cur_y=0
            endif 

! ...      generating random number
            do k=1,3
              dx(k)=ran1(iseed)
              dy(k)=ran2(iseed)
              dz(k)=ran3(iseed)
            enddo

            rndx=dx(1)
            rndy=dx(2)
            rndz=dx(3)
            xadv=(upar(i)+cur_x+grdx(i))*dtb
            yadv=(vpar(i)+cur_y+grdy(i))*dtb
            zadv=(wpar(i)+rsl+grdz(i))*dtb        
            xdif=(2*rndx-1)*sqrt(6*dhfx(i)*dtb) !random walk
            ydif=(2*rndy-1)*sqrt(6*dhfy(i)*dtb)
            zdif=(2*rndz-1)*sqrt(6*dhfz(i)*dtb)   
!          xadv=0;yadv=0;zadv=0;zdif=0
!          if(i.eq.1) then
!            iel=nnel
!            k=jlev
!            nd1=elnode(1,iel)
!            nd2=elnode(2,iel)
!            nd3=elnode(3,iel)
!            print*,iel,idry_e(iel),hvis_e(iel,k)
!            print*,hvis_p(nd1,k),hvis_p(nd2,k),hvis_p(nd3,k)
!            print*,grdx(i),grdy(i),grdz(i)
!            print*,dhfx(i),dhfy(i),dhfz(i)
!            pause  
!          endif

            xt=x0+xadv+xdif
            yt=y0+yadv+ydif
            zt=z0+zadv+zdif    

            !rnds - random number used for stranding (oil spill only)
            rnds=dy(1)          
          else !not oil spill
            xt=x0+ibf*dtb*upar(i)
            yt=y0+ibf*dtb*vpar(i)
            zt=z0+ibf*dtb*wpar(i)
            rnds=0 !not used
          endif !mod_part

          call quicksearch(1,idt,i,nnel0,jlev0,dtb,x0,y0,z0,xt,yt,zt,nnel,jlev, &
     &nodel2,arco,zrat,nfl,eta_p,dp_p,ztmp,kbpl,ist(i),inbr(i),rnds,pbeach)

!	  nnel not dry
!	  Interpolate in time
          do j=1,i34(nnel)
            nd=elnode(j,nnel)
            if(mod_part==1) then
              vwx(j)=wnx1(nd)*(1-trat)+wnx2(nd)*trat  ! windx
              vwy(j)=wny1(nd)*(1-trat)+wny2(nd)*trat  ! windy
            endif !mod_part
            do l=1,2
              lev=jlev+l-2
              vxl(j,l)=uu1(nd,lev)*(1-trat)+uu2(nd,lev)*trat
              vyl(j,l)=vv1(nd,lev)*(1-trat)+vv2(nd,lev)*trat
              vzl(j,l)=ww1(nd,lev)*(1-trat)+ww2(nd,lev)*trat
              if(mod_part==1) then
                val(j,l)=hf1(nd,lev)*(1-trat)+hf2(nd,lev)*trat !viscosity [m^2/s]
                vbl(j,l)=vf1(nd,lev)*(1-trat)+vf2(nd,lev)*trat !diffusivity
                vcl(j,l)=hf1(nd,lev)*(1-trat/2)+hf2(nd,lev)*trat/2 !viscosity @ half-step
                vdl(j,l)=vf1(nd,lev)*(1-trat/2)+vf2(nd,lev)*trat/2 !diffusivity @ half-step
              endif !mod_part
            enddo !l
          enddo !j

!	  Interpolate in vertical 
          do j=1,i34(nnel)
            vxn(j)=vxl(j,2)*(1-zrat)+vxl(j,1)*zrat
            vyn(j)=vyl(j,2)*(1-zrat)+vyl(j,1)*zrat
            vzn(j)=vzl(j,2)*(1-zrat)+vzl(j,1)*zrat
            if(mod_part==1) then
              van(j)=val(j,2)*(1-zrat)+val(j,1)*zrat !viscosity
              vcn(j)=vcl(j,2)*(1-zrat)+vcl(j,1)*zrat !viscosity @ half-step
              vdn(j)=vdl(j,2)*(1-zrat)+vdl(j,1)*zrat !diffusivity @ half-step
            endif !mod_part
          enddo !j

!	  Interpolate in horizontal
          upar(i)=0; vpar(i)=0; wpar(i)=0
          wndx(i)=0;wndy(i)=0; dhfx(i)=0; dhfz(i)=0
          do j=1,3 !1st tri for quads
            id=nodel2(j)
            upar(i)=upar(i)+vxn(id)*arco(j)
            vpar(i)=vpar(i)+vyn(id)*arco(j)
            wpar(i)=wpar(i)+vzn(id)*arco(j)
            if(mod_part==1) then
              wndx(i)=wndx(i)+vwx(id)*arco(j)  !windx
              wndy(i)=wndy(i)+vwy(id)*arco(j)  !windy
              dhfx(i)=dhfx(i)+vcn(id)*arco(j)  !viscosity @half-step
              dhfz(i)=dhfz(i)+vdn(id)*arco(j)  !diffusivity @half-step
            endif !mod_part
          enddo !j
          dhfy(i)=dhfx(i)

!         Compute the vertical gradient of diffusivity
          if(mod_part==1) then
            az=0
            do j=1,3
              id=nodel2(j)
              az=az+(vbl(id,2)-vbl(id,1))*arco(j)
            enddo !j
            tmp=ztmp(jlev)-ztmp(jlev-1)
            if(tmp==0) then
              grdz(i)=0
            else
              grdz(i)=az/tmp
            endif
            !Gradient of viscosity
            grdx(i)=dot_product(van(nodel2(1:3)),dldxy(1:3,1,nnel))
            grdy(i)=dot_product(van(nodel2(1:3)),dldxy(1:3,2,nnel))
          endif !mod_part

          if(nfl==1) then
            iabnorm(i)=1
            exit
          endif

          x0=xt
          y0=yt
          z0=zt
          nnel0=nnel
          jlev0=jlev
        enddo !idt=1,ndeltp

        xpar(i)=xt
        ypar(i)=yt
        zpar(i)=zt
        ielpar(i)=nnel
        levpar(i)=jlev

449     continue

        if(mod_part==1) then
          !Calculate remaining mass
          if(ist(i)/=0) then
            y0=1 ! initial mass
            yc=y0*remain_ratio
            arg=-log(2.0)/T_half*(time-st_p(i))
            amas(i)=yc+(y0-yc)*exp(max(arg,-20.d0))
          endif
        endif !mod_part        

        if(ics==2) then
          call cppinverse(xout, yout, xpar(i), ypar(i), slam0, sfea0)
          xout = xout * 180.0 / pi
          yout = yout * 180.0 / pi
        else
          xout = xpar(i)
          yout = ypar(i)
        endif

        !drogue format for xmvis6s; no extra lines after this
        write(95,'(i12,2(1x,e22.14),1x,f12.3)')i,xout,yout,zpar(i)-eta_p
!       write(95,*) i,ist(i),amas(i),xout,yout,real(zpar(i)-eta_p)
        write(97,*)i,real(xout),real(yout),real(zpar(i)),ielpar(i),levpar(i), &
     &real(eta_p),real(dp_p),iabnorm(i),real(upar(i)),real(vpar(i)),real(wpar(i))
        if(levpar(i)>0) then
          ie4=ielpar(i)
          write(97,*)'wet:',i34(ie4),real(arco(1:3)),real(uu2(elnode(1:i34(ie4),ie4),levpar(i))), &
     &real(vv2(elnode(1:i34(ie4),ie4),levpar(i))),real(eta3(elnode(1:i34(ie4),ie4)))
        endif !levpar

!!        write(95,'(2e14.4)')time,ztmp2(nvrt)
!!        write(*,'(2e14.4)')time,zpar(i)-eta3(ielpar(i))
      enddo !i=1,nparticle

!...  Store info for next step
      uu1=uu2; vv1=vv2; ww1=ww2; eta1=eta2
      hf1=hf2; vf1=vf2; wnx1=wnx2; wny1=wny2

      if(nscreen.eq.1) write(*,*)'Time (days)=',time/86400

!--------------------------------------------------------------------------
      enddo !it

      if(nscreen.eq.1) write(*,*)'Completed'

      stop
      end


! ****************************************************      
! RANDOM NUMBER GENERATORS from numerical recipes       
! ran[1-3] are from 3 methods
! ****************************************************
  FUNCTION ran1(idum)
   Integer idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,&
      & IR2,NTAB,NDIV
   Real*8 ran1,AM,EPS,RNMX
   Parameter(IM1=2147483563,IM2=2147483399,AM=1./&
      & IM1,IMM1=IM1-1,IA1=40014,IA2=40692,IQ1=53668,&
      & IQ2=52774,IR1=12211,IR2=3791,NTAB=32,NDIV=1+&
      & IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
   Integer idum2,j,k,iv(NTAB),iy
   SAVE iv,iy,idum2
   DATA idum2/123456789/,iv/NTAB*0/,iy/0/
! ....................................................
   if(idum.le.0) then
     idum=max(-idum,1)
     idum2=idum
     do j=NTAB+8,1,-1
        k=idum/IQ1
        idum=IA1*(idum-k*IQ1)-k*IR1
        if(idum.lt.0) idum=idum+IM1
        if(j.le.NTAB) iv(j)=idum
     enddo
     iy=iv(1) 
   endif
   k=idum/IQ1
   idum=IA1*(idum-k*IQ1)-k*IR1
   if(idum.lt.0) idum=idum+IM1
     k=idum2/IQ2
     idum2=IA2*(idum2-k*IQ2)-k*IR2
   if(idum2.lt.0) idum2=idum2+IM2
     j=1+iy/NDIV
     iy=iv(j)-idum2
     iv(j)=idum
   if(iy.lt.1) iy=iy+IMM1
     ran1=min(AM*iy,RNMX)  
   return
  end
      
  FUNCTION ran2(idum)
   Integer idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,&
      & IR2,NTAB,NDIV
   Real*8 ran2,AM,EPS,RNMX
   Parameter(IM1=2147483563,IM2=2147483399,AM=1./&
      & IM1,IMM1=IM1-1,IA1=40014,IA2=40692,IQ1=53668,&
      & IQ2=52774,IR1=12211,IR2=3791,NTAB=32,NDIV=1+&
      & IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
   Integer idum2,j,k,iv(NTAB),iy
   SAVE iv,iy,idum2
   DATA idum2/123456789/,iv/NTAB*0/,iy/0/
! ....................................................
   if(idum.le.0) then
     idum=max(-idum,1)
     idum2=idum
     do j=NTAB+8,1,-1
        k=idum/IQ1
        idum=IA1*(idum-k*IQ1)-k*IR1
        if(idum.lt.0) idum=idum+IM1
        if(j.le.NTAB) iv(j)=idum
     enddo
     iy=iv(1) 
   endif
   k=idum/IQ1
   idum=IA1*(idum-k*IQ1)-k*IR1
   if(idum.lt.0) idum=idum+IM1
     k=idum2/IQ2
     idum2=IA2*(idum2-k*IQ2)-k*IR2
   if(idum2.lt.0) idum2=idum2+IM2
     j=1+iy/NDIV
     iy=iv(j)-idum2
     iv(j)=idum
   if(iy.lt.1) iy=iy+IMM1
     ran2=min(AM*iy,RNMX)  
   return
  end

  FUNCTION ran3(idum)
   Integer idum
   Integer MBIG,MSEED,MZ
   Real*8 ran3,FAC 
   Parameter(MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
   Integer i,iff,ii,inext,inextp,k
   Integer mj,mk,ma(55) 
   SAVE iff,inext,inextp,ma
   DATA iff /0/
! ........................................................... 
   if(idum.lt.0.or.iff.eq.0) then
     iff=1
     mj=abs(MSEED-abs(idum))
     mj=mod(mj,MBIG)
     ma(55)=mj
     mk=1
     do i=1,54
        ii=mod(21*i,55)
        ma(ii)=mk
        mk=mj-mk
        if(mk.lt.MZ) mk=mk+MBIG
          mj=ma(ii)
     enddo   
     do k=1,4
       do i=1,55
          ma(i)=ma(i)-ma(1+mod(i+30,55))
          if(ma(i).lt.MZ) ma(i)=ma(i)+MBIG
       enddo
     enddo
     inext=0
     inextp=31
     idum=1
   endif  
     inext=inext+1
     if(inext.eq.56) inext=1
       inextp=inextp+1
     if(inextp.eq.56) inextp=1
       mj=ma(inext)-ma(inextp)
     if(mj.lt.MZ) mj=mj+MBIG
       ma(inext)=mj
       ran3=mj*FAC
   return
  end

!******************************************************************************
!                                                                             *
!    Transform from lon,lat (rlambda,phi) coordinates into CPP coordinates.     *
!    Lon,Lat must be in radians.                                              *
!                                                                             *
!******************************************************************************

      subroutine cpp(x,y,rlambda,phi,rlambda0,phi0)
      use kind_par
      implicit real(kind=dbl_kind1)(a-h,o-z), integer(i-n)

      r=6378206.4
      x=r*(rlambda-rlambda0)*cos(phi0)
      y=phi*r

      return
      end

      subroutine cppinverse(rlambda,phi,x,y,rlambda0,phi0)
      use kind_par
      implicit real(kind=dbl_kind1)(a-h,o-z), integer(i-n)

      r=6378206.4
      rlambda=x / (r * cos(phi0)) + rlambda0
      phi=y/r

      return
      end


!
!********************************************************************************
!										*
!     Program to detect if two segments (1,2) and (3,4) have common pts   	*
!     Assumption: the 4 pts are distinctive.					*
!     The eqs. for the 2 lines are: X=X1+(X2-X1)*tt1 and X=X3+(X4-X3)*tt2.	*
!     Output: iflag: 0: no intersection or colinear; 1: exactly 1 intersection.	*
!     If iflag=1, (xin,yin) is the intersection.				*
!										*
!********************************************************************************
!

      subroutine intersect2(x1,x2,x3,x4,y1,y2,y3,y4,iflag,xin,yin,tt1,tt2)
      use kind_par
      implicit real(kind=dbl_kind1)(a-h,o-z), integer(i-n)
      real(kind=dbl_kind1), parameter :: zero1=0.0 !small positive number or 0

      real(kind=dbl_kind1), intent(in) :: x1,x2,x3,x4,y1,y2,y3,y4
      integer, intent(out) :: iflag
      real(kind=dbl_kind1), intent(out) :: xin,yin,tt1,tt2

      tt1=-1000
      tt2=-1000
      iflag=0
      delta=(x2-x1)*(y3-y4)-(y2-y1)*(x3-x4)
      delta1=(x3-x1)*(y3-y4)-(y3-y1)*(x3-x4)
      delta2=(x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)

      if(delta.ne.0.0d0) then
        tt1=delta1/delta
        tt2=delta2/delta
        if(tt1.ge.-zero1.and.tt1.le.1+zero1.and.tt2.ge.-zero1.and.tt2.le.1+zero1) then
          iflag=1
          xin=x1+(x2-x1)*tt1
          yin=y1+(y2-y1)*tt1
        endif
      endif

      return
      end



       function signa(x1,x2,x3,y1,y2,y3)
       use kind_par
       implicit real(kind=dbl_kind1)(a-h,o-z),integer(i-n)

       signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2

       return
       end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                       !
!       Compute area coordinates of pt (xt,yt), which must be inside element nnel.      !
!       Impose bounds for area coordinates.                                             !
!                                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine area_coord(nnel,xt,yt,arco)
      use global, only : dbl_kind,elnode,area,x,y
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)
      integer, intent(in) :: nnel
      real(kind=dbl_kind), intent(in) :: xt,yt
      real(kind=dbl_kind), intent(out) :: arco(3)

      n1=elnode(1,nnel)
      n2=elnode(2,nnel)
      n3=elnode(3,nnel)
      arco(1)=signa(xt,x(n2),x(n3),yt,y(n2),y(n3))/area(nnel)
      arco(2)=signa(x(n1),xt,x(n3),y(n1),yt,y(n3))/area(nnel)
      arco(1)=max(0.0d0,min(1.0d0,arco(1)))
      arco(2)=max(0.0d0,min(1.0d0,arco(2)))
      if(arco(1)+arco(2)>1) then
        arco(3)=0
        arco(1)=min(1.d0,max(0.d0,arco(1)))
        arco(2)=1-arco(1)
      else
        arco(3)=1-arco(1)-arco(2)
      endif

      return
      end


!
!********************************************************************
!								    *
!	Routine to update z-coordinates and wetting and drying      *
!								    *
!********************************************************************
!
      subroutine levels
      use global
      use compute_zcor
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)

      dimension idry_new(np) !,out2(12+mnv)
      real(dbl_kind) :: zlcl(nvrt)

!...  z-coor. for nodes
!...  
      do i=1,np
        if(dp(i)+eta3(i)<=h0) then !dry
          idry_new(i)=1 
          if(ivcor==2) then; if(dp(i)>=h_s) then
            write(11,*)'Deep depth dry:',i,h_s,dp(i)
            stop
          endif; endif
!          kbp(i)=0
        else !wet
          idry_new(i)=0

          if(ivcor==2) then
            call zcor_SZ_double(dp(i),eta3(i),h0,h_s,h_c,theta_b, &
     &theta_f,kz,nvrt,ztot,sigma,zlcl,idry_tmp,kbpl)
            if(idry_tmp==1) then
              write(11,*)'Impossible dry (7):',i,idry_tmp,dp(i),eta1(i),eta2(i),eta3(i),kbpl
              stop
            endif
            !Cannot use kbp00 b/cos wet/dry
            kbp(i)=kbpl
            z(i,kbpl:nvrt)=zlcl(kbp(i):nvrt)
          else if(ivcor==1) then
            z(i,kbp(i):nvrt)=(eta3(i)+dp(i))*sigma_lcl(kbp(i):nvrt,i)+eta3(i)
          else
            write(11,*)'Unknown ivcor:',ivcor
            stop
          endif

          do k=kbp(i)+1,nvrt
            if(z(i,k)-z(i,k-1)<=0) then
              write(11,*)'Inverted z-levels at:',i,k,z(i,k)-z(i,k-1),eta3(i)
              stop
            endif
          enddo !k
        endif !wet ot dry
      enddo !i=1,np

!...  Set wet/dry flags for element; element is "dry" if one of nodes is dry; conversely, 
!...  an element is wet if all nodes are wet (and all sides are wet as well)
!...  Weed out fake we nodes; a node is wet if and only if at least one surrounding element is wet
!...
!      idry_e0=idry_e !save
      idry=1 !dry unless wet
      kbe=0
      do i=1,ne
        idry_e(i)=maxval(idry_new(elnode(1:i34(i),i)))
        if(idry_e(i)==0) then
          idry(elnode(1:i34(i),i))=0
          kbe(i)=minval(kbp(elnode(1:i34(i),i)))
        endif
      enddo !i

      return
      end

!
!********************************************************************************
!	
!     Straightline search algorithm. Initially nnel0 is an element that 	
!     encompasses (x0,y0). iloc=0: do not nudge initial pt; iloc=1: nudge.	 
!     Input: iloc,idt,ipar,nnel0,x0,y0,z0,xt,yt,zt,jlev0, time, and uu2,vv2,ww2 for 	
!	abnormal cases;								
!     Output the updated end pt (xt,yt,zt) (if so), nnel1, jlev1, area          
!       coordinates, vertical ratio, a flag nfl, and local elevation and depth. 
!       I/O for oil spill: ist2,inbr2,rnds,pbeach
!     nfl=1 if a bnd or dry element is hit and vel. there is small,		 
!	or death trap is reached.						
!										
!********************************************************************************
!
      subroutine quicksearch(iloc,idt,ipar,nnel0,jlev0,time,x0,y0,z0,xt,yt,zt,nnel1,jlev1, &
     &nodel2,arco,zrat,nfl,etal,dp_p,ztmp,kbpl,ist2,inbr2,rnds,pbeach)

      use global
      use compute_zcor
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)

      integer, intent(in) :: iloc,idt,ipar,nnel0,jlev0
      real(kind=dbl_kind), intent(in) :: time,x0,y0,z0,rnds,pbeach
      integer, intent(out) :: nnel1,jlev1,nodel2(3),nfl,kbpl,ist2,inbr2
      real(kind=dbl_kind), intent(inout) :: xt,yt,zt
      real(kind=dbl_kind), intent(out) :: arco(3),zrat,etal,dp_p,ztmp(nvrt)

      !Local
      real(dbl_kind) :: zlcl(nvrt)
      logical ::  ltmp1,ltmp2

      if(iloc>1) then
        write(11,*)'iloc > 1'
        stop
      endif
      if(idry_e(nnel0)==1) then
        write(11,*)'Starting element is dry'
        stop
      endif

      nfl=0
      trm=time !time remaining

!     Starting element nel
      nel=nnel0

!     An interior pt close to (x0,y0) to prevent underflow for iloc >=1.
      if(iloc==0) then
        xcg=x0
        ycg=y0
      else if(iloc==1) then
        xcg=(1-1.0d-4)*x0+1.0d-4*xctr(nel)
        ycg=(1-1.0d-4)*y0+1.0d-4*yctr(nel)
      endif

      call pt_in_poly2(i34(nel),x(elnode(1:i34(nel),nel)),y(elnode(1:i34(nel),nel)),xcg,ycg,inside)
!      aa=0
!      aa1=0
!      do i=1,3
!        n1=elnode(i,nel)
!        n2=elnode(nxq(1,i,i34(nel)),nel)
!        aa=aa+dabs(signa(x(n1),x(n2),x0,y(n1),y(n2),y0))
!        aa1=aa1+dabs(signa(x(n1),x(n2),xt,y(n1),y(n2),yt))
!      enddo !i
!      ae=dabs(aa-area(nel))/area(nel)
!      if(ae>small1) then

      if(inside==0) then
        write(11,*)'(x0,y0) not in nnel0 initially',nnel0,xcg,ycg
        stop
      endif

      call pt_in_poly2(i34(nel),x(elnode(1:i34(nel),nel)),y(elnode(1:i34(nel),nel)),xt,yt,inside)
!      ae=dabs(aa1-area(nel))/area(nel)
!      if(ae<small1) then
      if(inside/=0) then
        nnel1=nel
        go to 400
      endif

!     (xt,yt) not in nel, and thus (xcg,ycg) and (xt,yt) are distinctive
      pathl=dsqrt((xt-xcg)**2+(yt-ycg)**2)
      if(pathl==0) then
        write(11,*)'Zero path',x0,y0,xt,yt,xcg,ycg
        stop
      endif

!     Starting edge nel_j
      nel_j=0
      do j=1,i34(nel)
        jd1=elnode(nxq(1,j,i34(nel)),nel)
        jd2=elnode(nxq(2,j,i34(nel)),nel)
        call intersect2(xcg,xt,x(jd1),x(jd2),ycg,yt,y(jd1),y(jd2),iflag,xin,yin,tt1,tt2)
        if(iflag==1) then
          nel_j=j
          exit
        endif
      enddo !j=1,3
      if(nel_j==0) then
        write(11,*)'Found no intersecting edges I:',nel,xcg,ycg,xt,yt
        stop
      endif

      zin=z0 !intialize
      it=0
      loop4: do
!----------------------------------------------------------------------------------------
      it=it+1
      if(it>1000) then
        if(ifort12(3)==0) then
          ifort12(3)=1
          write(12,*)'Death trap reached',idt
        endif
        nfl=1
        xt=xin
        yt=yin
        zt=zin
        nnel1=nel
        exit loop4
      endif
      md1=elnode(nxq(1,nel_j,i34(nel)),nel)
      md2=elnode(nxq(2,nel_j,i34(nel)),nel)
      
!     Compute z position 
      dist=sqrt((xin-xt)**2+(yin-yt)**2)
!      if(dist/pathl.gt.1+1.0d-4) then
!        write(11,*)'Path overshot'
!        stop
!      endif
      tmp=min(1.d0,dist/pathl)
      zin=zt-tmp*(zt-zin)
      trm=trm*dist/pathl !time remaining
      
      pathl=sqrt((xin-xt)**2+(yin-yt)**2)
      if(pathl==0.or.trm==0) then
        write(11,*)'Target reached'
        stop
      endif

      lit=0 !flag
!     For horizontal exit and dry elements
      if(ic3(nel_j,nel)==0.or.idry_e(max(1,ic3(nel_j,nel)))==1) then
        lit=1
        isd=elside(nel_j,nel)
        if(isidenode(1,isd)+isidenode(2,isd)/=md1+md2) then
          write(11,*)'Wrong side'
          stop
        endif

        !Oil spill
        if(mod_part==1) then
          !Set status flag
          if(ic3(nel_j,nel)==0) then
            ltmp1=isbnd(md1)>0.or.isbnd(md2)>0 !open bnd
            ltmp2=isbnd(md1)==-1.and.isbnd(md2)==-1.and.rnds>=pbeach !permanently stranded @ land bnd
            if(ltmp1.or.ltmp2) then
              if(ltmp1) then
                ist2=2
              else
                ist2=-1 !permanently stranded
              endif
              nfl=1
              xt=xin
              yt=yin
              zt=zin
              nnel1=nel
              exit loop4
            endif
          else !internal side with a dry elem.
            if(rnds>=pbeach) then !% exceeded
              ist2=-2  !stranded onshore
              inbr2=ic3(nel_j,nel)
              nfl=1
              xt=xin
              yt=yin
              zt=zin
              nnel1=nel
              exit loop4
            endif
          endif !ic3
        endif !mod_part

!       Reflect off
!        eps=1.e-2
!        xin=(1-eps)*xin+eps*xctr(nel)
!        yin=(1-eps)*yin+eps*yctr(nel)
        xcg=xin
        ycg=yin

        !Original vel
        uvel0=(xt-xin)/trm
        vvel0=(yt-yin)/trm
        vnorm=uvel0*snx(isd)+vvel0*sny(isd)
        vtan=-uvel0*sny(isd)+vvel0*snx(isd)
        !vtan=-(uu2(md1,jlev0)+uu2(md2,jlev0))/2*sny(isd)+(vv2(md1,jlev0)+vv2(md2,jlev0))/2*snx(isd)
        !Reverse normal vel
        vnorm=-vnorm

        !tmp=max(abs(vtan),1.d-2) !to prevent getting stuck
        !vtan=tmp*sign(1.d0,tmp)
        xvel=vnorm*snx(isd)-vtan*sny(isd)
        yvel=vnorm*sny(isd)+vtan*snx(isd)
        zvel=(ww2(md1,jlev0)+ww2(md2,jlev0))/2
        xt=xin+xvel*trm
        yt=yin+yvel*trm
        zt=zin+zvel*trm
!        hvel=dsqrt(xvel**2+yvel**2)
!        if(hvel<1.e-4) then
!          write(11,*)'Impossible (5):',hvel
!          nfl=1
!          xt=xin
!          yt=yin
!          zt=zin
!          nnel1=nel
!          exit loop4
!        endif
        !pathl unchanged since hvel is unchanged
!        pathl=hvel*trm
      endif !abnormal cases

!     Search for nel's neighbor with edge nel_j, or in abnormal cases, the same element
      if(lit==0) nel=ic3(nel_j,nel) !next front element

!      aa=0
!      do i=1,3
!        k1=elnode(i,nel)
!        k2=elnode(nxq(1,i,i34(nel)),nel)
!        aa=aa+dabs(signa(x(k1),x(k2),xt,y(k1),y(k2),yt))
!      enddo !i
!      ae=dabs(aa-area(nel))/area(nel)
!      if(ae<small1) then

      call pt_in_poly2(i34(nel),x(elnode(1:i34(nel),nel)),y(elnode(1:i34(nel),nel)),xt,yt,inside)
      if(inside/=0) then
        nnel1=nel
        exit loop4
      endif

!     Next intersecting edge
      do j=1,i34(nel)
         jd1=elnode(nxq(1,j,i34(nel)),nel)
         jd2=elnode(nxq(2,j,i34(nel)),nel)
!        For abnormal case, same side (border side) cannot be hit again
         if(jd1==md1.and.jd2==md2.or.jd2==md1.and.jd1==md2) cycle
         call intersect2(xcg,xt,x(jd1),x(jd2),ycg,yt,y(jd1),y(jd2),iflag,xin,yin,tt1,tt2)
         if(iflag==1) then
           nel_j=j !next front edge          
           cycle loop4
         endif
      enddo !j
      write(11,*)'Failed to find next edge I:',lit,xin,yin,xt,yt,nel,md1,md2,idt,rnds
      stop
!----------------------------------------------------------------------------------------
      end do loop4 

400   continue
!     No vertical exit from domain
      if(idry_e(nnel1)==1) then
        write(11,*)'Ending element is dry'
        stop
      endif

!     Compute area & sigma coord.
      !call area_coord(nnel1,xt,yt,arco)
      call pt_in_poly3(i34(nnel1),x(elnode(1:i34(nnel1),nnel1)),y(elnode(1:i34(nnel1),nnel1)), &
     &xt,yt,arco,nodel2)
      n1=elnode(nodel2(1),nnel1)
      n2=elnode(nodel2(2),nnel1)
      n3=elnode(nodel2(3),nnel1)
      etal=eta3(n1)*arco(1)+eta3(n2)*arco(2)+eta3(n3)*arco(3)
      dep=dp(n1)*arco(1)+dp(n2)*arco(2)+dp(n3)*arco(3)
      dp_p=dep
      if(etal+dep<h0) then
        write(11,*)'Weird wet element in quicksearch:',nnel1,eta3(n1),eta3(n2),eta3(n3),h0
        stop
      endif

      if(istiff==1) zt=etal+zpar0(ipar)

!     Compute z-levels
      if(ivcor==2) then
        call zcor_SZ_double(dep,etal,h0,h_s,h_c,theta_b, &
     &theta_f,kz,nvrt,ztot,sigma,zlcl,idry_tmp,kbpl)
        ztmp(kbpl:nvrt)=zlcl(kbpl:nvrt)
      else if(ivcor==1) then
        kbpl=nvrt+1 !local bottom index (maybe degenerate)
        do j=1,3
          nd=elnode(nodel2(j),nnel1)
          if(kbp(nd)<kbpl) kbpl=kbp(nd)
        enddo !j
        !ztmp(kbpl)=-dep 
        !ztmp(nvrt)=etal
        do k=kbpl,nvrt
          ztmp(k)=0
          do j=1,3
            nd=elnode(nodel2(j),nnel1)
            tmp=(eta3(nd)+dp(nd))*sigma_lcl(max(kbp(nd),k),nd)+eta3(nd)
            ztmp(k)=ztmp(k)+tmp*arco(j)
          enddo !j
        enddo !k
      else
        write(11,*)'Unknown ivcor:',ivcor
        stop
      endif

      do k=kbpl+1,nvrt
        if(ztmp(k)-ztmp(k-1)<0) then !can be 0 for deg. case
          write(11,*)'Inverted z-level in quicksearch:',nnel1,etal,dep,k,ztmp(k)-ztmp(k-1),ztmp(kbpl:nvrt)
          stop
        endif
      enddo !k

      if(zt<=ztmp(kbpl)) then
        !Avoid getting stuck at bottom
        zt=ztmp(kbpl+1)
        zrat=0
        jlev1=kbpl+1
      else if(zt>=ztmp(nvrt)) then
        zt=ztmp(nvrt)
        zrat=0
        jlev1=nvrt
      else
        jlev1=0
        do k=kbpl,nvrt-1
          if(zt>=ztmp(k).and.zt<=ztmp(k+1)) then 
            jlev1=k+1
            exit
          endif
        enddo !k
        if(jlev1==0) then
          write(11,*)'Cannot find a vert. level:',zt,etal,dep
          write(11,*)(ztmp(k),k=kbpl,nvrt)
          stop
        endif
        zrat=(ztmp(jlev1)-zt)/(ztmp(jlev1)-ztmp(jlev1-1))
      endif

      if(zrat<0.or.zrat>1) then
        write(11,*)'Sigma coord. wrong (4):',jlev1,zrat
        stop
      endif

!      if(kbpl==kz) then !in pure S region
!        ss=(1-zrat)*sigma(jlev1-kz+1)+zrat*sigma(jlev1-kz)
!      else
!        ss=-99
!      endif

      end subroutine quicksearch

!======================================================================
      subroutine pt_in_poly2(i34,x,y,xp,yp,inside)
!     (Double-precision) Routine to perform point-in-polygon
!     (triangle/quads) test.
!     Inputs:
!            i34: 3 or 4 (type of elem)
!            x(i34),y(i34): coord. of polygon/elem. (counter-clockwise)
!            xp,yp: point to be tested
!     Outputs:
!            inside: 1, inside
      use global, only : small1
      implicit real*8(a-h,o-z)
      integer, intent(in) :: i34
      real*8, intent(in) :: x(i34),y(i34),xp,yp
      integer, intent(out) :: inside

      real*8 :: swild(i34)

      inside=0
      do j=1,i34
        j1=j+1
        if(j1>i34) j1=j1-i34
        swild(j)=signa(x(j),x(j1),xp,y(j),y(j1),yp)
      enddo !j
      ae=minval(swild(1:i34))
      if(ae>-small1) inside=1

      end subroutine pt_in_poly2
     
!======================================================================
      subroutine pt_in_poly3(i34,x,y,xp,yp,arco,nodel)
!     (Double-precision) Routine to perform point-in-polygon
!     (triangle/quads) test with assumption that it's inside, and calculate the area coord.
!     (for quad, split it into 2 triangles and return the 3 nodes and
!     area coord.)
!     Inputs:
!            i34: 3 or 4 (type of elem)
!            x(i34),y(i34): coord. of polygon/elem. (counter-clockwise)
!            xp,yp: point to be tested
!     Outputs:
!            arco(3), nodel(3) : area coord. and 3 local node indices (valid only if inside)
      implicit real*8(a-h,o-z)
      integer, intent(in) :: i34
      real*8, intent(in) :: x(i34),y(i34),xp,yp
      integer, intent(out) :: nodel(3)
      real*8, intent(out) :: arco(3)

      !Local
      integer :: list(3)
      real*8 :: ar(2),swild(2,3)

      !Areas of up to 2 triangles
      ar(1)=signa(x(1),x(2),x(3),y(1),y(2),y(3))
      ar(2)=0 !init
      if(i34==4) ar(2)=signa(x(1),x(3),x(4),y(1),y(3),y(4))
      if(ar(1)<=0.or.i34==4.and.ar(2)<=0) then
        print*, 'Negative area:',i34,ar,x,y
        stop
      endif

      ae_min=huge(1.0d0)
      do m=1,i34-2 !# of triangles
        if(m==1) then
          list(1:3)=(/1,2,3/) !local indices
        else !quads
          list(1:3)=(/1,3,4/)
        endif !m
        aa=0
        do j=1,3
          j1=j+1
          j2=j+2
          if(j1>3) j1=j1-3
          if(j2>3) j2=j2-3
          swild(m,j)=signa(x(list(j1)),x(list(j2)),xp,y(list(j1)),y(list(j2)),yp) !temporary storage
          aa=aa+abs(swild(m,j))
        enddo !j=1,3

        ae=abs(aa-ar(m))/ar(m)
        if(ae<=ae_min) then
          ae=ae_min
          nodel(1:3)=list(1:3)
          arco(1:3)=swild(m,1:3)/ar(m)
          arco(1)=max(0.d0,min(1.d0,arco(1)))
          arco(2)=max(0.d0,min(1.d0,arco(2)))
          if(arco(1)+arco(2)>1) then
            arco(3)=0
            arco(2)=1-arco(1)
          else
            arco(3)=1-arco(1)-arco(2)
          endif
        endif
      enddo !m

      end subroutine pt_in_poly3

!======================================================================
