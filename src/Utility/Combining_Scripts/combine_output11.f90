!===============================================================================
! Read in netcdf outputs (rank-specific) from parallel code and combine them into
! one global outputs in classic nc format. 
! Global-local mappings are read in from separate files.
! Run this program inside the directory outputs/, where some of the input files below
! can be found.

! Usage: ./combine_output11 -h for help
! Inputs:
!        rank-specific nc files (from SCHISM outputs/); 
!        local_to_global_* (from SCHISM outputs/);
! Output: combined nc file
!

!  ifort -cpp -O2 -assume byterecl -o combine_output11.exe ../UtilLib/argparse.f90 ../UtilLib/schism_geometry.f90 combine_output11.f90 -I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff

!  History: 
!          2018-1      Richard Hofmeister changed the combine method (all vars read in at a
!                      time) in order to speed up combine on many cores. This
!                      version requires more memory.
!          2017-9      revamped whole thing (all-in-1 pre- and post-combine
!          files)
!          2015-12-04  changed mesh, dim and var names, added support of LSC2 and mixed triangle and quad mesh
!                      used standard and long name as CF convention

!           2013-12-09 drf: NF_FLOATs, add SZ metadata, netcdf error checking towards ugrid v0.9
!              Add sigma attributes and sigma_theta_f, sigma_theta_b, sigma_h_c' variables for CF 
!              Add mesh variable for ugrid compliance
!              Parse binary file time (from bctides.in header) into udunits2 compatible units  
!   
!           changed non-standard outputs which can be viz'ed with vis6 now (vertical
!           coordinates not dependable).

!           From combine_output9: netcdf
subroutine combine_output11(ibgn,iend,iwetdry,to_be_combined,output_prefix)
!-------------------------------------------------------------------------------
  use netcdf
  use schism_geometry_mod
  implicit real(4)(a-h,o-z),integer(i-n)
  
!  integer                              :: nfile  !< number of file base names
!  character(len=30),dimension(nfile)   :: files  !< base names (e.g. elev.61) of files to combine
  integer, intent(in)                  :: ibgn   !< first output file index to process
  integer, intent(in)                  :: iend   !< last output file index to process
  integer, intent(in)                  :: iwetdry !dry option (0: use last-wet value; 1: junk)
  
  character(len=30) :: file63
  character(len=12) :: it_char
  character(len=48) :: start_time,version,variable_dim
  character(len=48) :: data_format
  character(len=72) :: fgb,fgb2,fdb  ! Processor specific global output file name
  integer :: start_year,start_month,start_day
  real*8 :: start_hour,utc_start
  integer :: lfgb,lfdb       ! Length of processor specific global output file name
  character(len=4) :: a_4
  integer,allocatable :: elnode(:,:)
  integer,allocatable :: elside(:,:)
  integer,allocatable :: isdel(:,:),vlen(:),ne(:),np(:),ns(:),ihot_len(:), &
 &i34(:),nm2(:,:),kbp00(:),iplg(:,:),ielg(:,:),islg(:,:),kbs(:),kbe(:), &
 &ic3(:,:),isidenode(:,:),i23d(:),ivs(:),ndims(:),iu_id(:),idry(:),idry_e(:),idry_s(:)
  real,allocatable :: ztot(:),sigma(:),cs(:),eta2(:),outeb(:,:,:), &
 &outsb(:,:,:),worka(:,:,:),xctr(:),yctr(:),dpe(:),x(:),y(:),dp(:),xcj(:),ycj(:), &
 &dps(:),eta2s(:),eta2e(:)
  character(len=48),allocatable :: variable_nm(:)
  character(len=1024)           :: to_be_combined
  character(len=1024)           :: default_variables='time,wetdry_elem,depth'
  character(len=1024)           :: output_prefix
  logical                       :: check_vars=.false.
  logical,allocatable           :: skip_var(:)

  type :: type_outd
    real,pointer :: data(:,:,:)=>null()
  end type

  type(type_outd),allocatable :: outd(:)

  integer :: invalid_index
  
  !netcdf variables
  integer :: coords_type ! lat/long or project coords
  character(len=48) :: lat_coord_standard_name
  character(len=48) :: lon_coord_standard_name
  character(len=48) :: x_units,y_units
  integer ::  lat_str_len,lon_str_len ! length of coord name
  character(len=50) :: fname
  character(len=256) :: cbuffer, stdname,units,longname
  character(len=8) :: varname,vname,uname
  integer :: time_dims(1),ele_dims(2),x_dims(1),y_dims(1),z_dims(1),sigma_dims(1), &
            &var1d_dims(1),var2d_dims(2),var3d_dims(3),var4d_dims(4), &
            &data_start_1d(1),data_start_2d(2),data_start_3d(3),data_start_4d(4), &
            &data_count_1d(1),data_count_2d(2),data_count_3d(3),data_count_4d(4), &
            & int_buffer(4),dummy_dim(1),ihgrid_id, tempint_array(1),&
            & chunks(3),one_dim
  real  :: hc_array(1),hs_array(1),thetab_array(1),thetaf_array(1),real_buffer(4) 
  
!  character(len=long_name_len) :: netcdf_var_long_name(2) 
!  character(len=long_name_len) :: netcdf_var_standard_name(2) 
!  character(len=dataset_name_len) :: netcdf_var_name
  character(len=4) :: netcdf_out_location
  character(len=5) :: netcdf_level_location
  integer          :: varid,dimid
!  integer :: netcdf_var_dim
!  logical :: found_netcdf_var

  invalid_index = -99999

! Read local_to_global_0000 for global info
  open(10,file='local_to_global_0000',status='old')
  read(10,*)ns_global,ne_global,np_global,nvrt,nproc !,ntracers
  close(10)
  allocate(x(np_global),y(np_global),dp(np_global),kbp00(np_global),kbe(ne_global),i34(ne_global), &
  &np(0:nproc-1),ns(0:nproc-1),ne(0:nproc-1),elnode(4,ne_global),nm2(4,ne_global), &
  &ztot(nvrt),sigma(nvrt),cs(nvrt),ihot_len(0:nproc-1), &
  &dpe(ne_global),xctr(ne_global),yctr(ne_global), &
  &idry(np_global),idry_e(ne_global),idry_s(ns_global),stat=istat)
  if(istat/=0) stop 'Allocation error: x,y'

  elnode=invalid_index
  !-------------------------------------------------------------------------------
  ! Read rank-specific local_to_global*
  !-------------------------------------------------------------------------------
  ! Read in local-global mappings from all ranks
  fdb='local_to_global_0000'
  lfdb=len_trim(fdb)

  !Find max. for dimensioning
  do irank=0,nproc-1
    write(fdb(lfdb-3:lfdb),'(i4.4)') irank
    open(10,file=fdb,status='old')
    read(10,*) !global info
    read(10,*) !info
    read(10,*)ne(irank)
    do i=1,ne(irank)
      read(10,*)!j,ielg(irank,i)
    enddo !i
    read(10,*)np(irank)
    do i=1,np(irank)
      read(10,*)
    enddo !i
    read(10,*)ns(irank)
    close(10)
  enddo !irank
  np_max=maxval(np(:))
  ns_max=maxval(ns(:))
  ne_max=maxval(ne(:))

  allocate(iplg(0:nproc-1,np_max),ielg(0:nproc-1,ne_max),islg(0:nproc-1,ns_max),stat=istat)
  if(istat/=0) stop 'Allocation error (2)'

  !Re-read
  !ns_global=0
  do irank=0,nproc-1
    write(fdb(lfdb-3:lfdb),'(i4.4)') irank
    open(10,file=fdb,status='old')
    read(10,*) !global info
    read(10,*) !info
    read(10,*)ne(irank)
    do i=1,ne(irank)
      read(10,*)j,ielg(irank,i)
    enddo !i
    read(10,*)np(irank)
    do i=1,np(irank)
      read(10,*)j,iplg(irank,i)
    enddo
    read(10,*)ns(irank) !sides
    do i=1,ns(irank)
      read(10,*)j,islg(irank,i)
      !if(ns_global<islg(irank,i)) ns_global=islg(irank,i)
    enddo

    read(10,*) !'Header:'
    read(10,*)start_year,start_month,start_day,start_hour,utc_start 
    version='v10'
!    read(10,'(a)')version  ! hardcoded as 'description' by schism_init.F90
!    read(10,'(a)')start_time
    read(10,*)nrec,dtout,nspool,nvrt,kz,h0,h_s,h_c,theta_b,theta_f,ics
    read(10,*)(ztot(k),k=1,kz-1),(sigma(k),k=1,nvrt-kz+1)
    read(10,*)np(irank),ne(irank),(x(iplg(irank,m)),y(iplg(irank,m)), &
  &dp(iplg(irank,m)),kbp00(iplg(irank,m)),m=1,np(irank)), &
  &(i34(ielg(irank,m)),(nm2(mm,m),mm=1,i34(ielg(irank,m))),m=1,ne(irank))

    close(10)

    !Debug
    !write(98,*)irank,(i34(ielg(irank,m)),m=1,ne(irank))

    !Compute C(s) for output
    do klev=kz,nvrt
      k=klev-kz+1
      cs(k)=(1-theta_b)*sinh(theta_f*sigma(k))/sinh(theta_f)+ &
    &theta_b*(tanh(theta_f*(sigma(k)+0.5))-tanh(theta_f*0.5))/2/tanh(theta_f*0.5)
    enddo !klev

    !Compute kbp00 (to avoid mismatch of indices) - larger rank prevails
!    do m=1,np(irank)
!      ipgb=iplg(irank,m)
!      kbp00(ipgb)=kbp01(irank,m)
!    enddo !m
 
    !   Reconstruct connectivity table
    do m=1,ne(irank)
      iegb=ielg(irank,m)
      if(iegb>ne_global) stop 'Overflow!'
      do mm=1,i34(iegb)
        itmp=nm2(mm,m)
        if(itmp>np(irank).or.itmp<=0) then
          write(*,*)'Overflow:',m,mm,itmp
          stop
        endif
        elnode(mm,iegb)=iplg(irank,itmp)
      enddo !mm
    enddo !m
  enddo !irank=0,nproc-1

  ! Compute geometry
  call compute_nside(np_global,ne_global,i34,elnode(1:4,1:ne_global),ns2)
  allocate(ic3(4,ne_global),elside(4,ne_global),isdel(2,ns2),isidenode(2,ns2),xcj(ns2),ycj(ns2), &
  &worka(2,nvrt,ns2),stat=istat)
  if(istat/=0) stop 'Allocation error: side(0)'
  call schism_geometry_single(np_global,ne_global,ns2,x,y,i34,elnode(1:4,1:ne_global),ic3(1:4,1:ne_global), &
  &elside(1:4,1:ne_global),isdel,isidenode,xcj,ycj)

  if(ns2/=ns_global) then
    write(*,*)'Mismatch in side:',ns2,ns_global
    stop
  endif
  !For dimensioning purpose
  if(np_global>ns_global.or.ne_global>ns_global) stop 'ns is not largest'

  ! Allocate side arrays
  allocate(dps(ns_global),kbs(ns_global),stat=istat)
  if(istat/=0) stop 'Allocation error: side'

  ! Compute side/element bottom index
  do i=1,ne_global
    kbe(i)=minval(kbp00(elnode(1:i34(i),i))) 
    dpe(i)=sum(dp(elnode(1:i34(i),i)))/i34(i)
    xctr(i)=sum(x(elnode(1:i34(i),i)))/i34(i)
    yctr(i)=sum(y(elnode(1:i34(i),i)))/i34(i)
  enddo !i
  do i=1,ns_global
    kbs(i)=minval(kbp00(isidenode(1:2,i)))
    dps(i)=sum(dp(isidenode(1:2,i)))/2
  enddo !i

  !-------------------------------------------------------------------------------
  ! Time iteration -- select "node" data
  !-------------------------------------------------------------------------------

  ! Loop over input files (stacks)
  do iinput=ibgn,iend
    write(it_char,'(i12)')iinput
    it_char=adjustl(it_char)  !place blanks at end
    it_len=len_trim(it_char)  !length without trailing blanks
    fname=trim(output_prefix)//'_'//it_char(1:it_len)//'.nc'
    !fname=it_char(1:it_len)//'_hvel.nc'
    fname=adjustl(fname)

!   Some compiles do not like 3 arguemnts in OR()
    iret = nf90_create(trim(fname), OR(NF90_CLOBBER,OR(NF90_NETCDF4,NF90_CLASSIC_MODEL)), ncid)
    !define dimensions
    iret = nf90_def_dim(ncid, 'nSCHISM_hgrid_node',np_global, node_dim)
    iret = nf90_def_dim(ncid, 'nSCHISM_hgrid_face',ne_global, nele_dim)
    iret = nf90_def_dim(ncid, 'nSCHISM_hgrid_edge',ns_global, nedge_dim)
    iret = nf90_def_dim(ncid, 'nMaxSCHISM_hgrid_face_nodes',4,nfour_dim)
    iret = nf90_def_dim(ncid, 'nSCHISM_vgrid_layers',nvrt, nv_dim)
    iret = nf90_def_dim(ncid, 'one',1, one_dim)
    iret = nf90_def_dim(ncid, 'two',2, ntwo_dim)
    iret = nf90_def_dim(ncid, 'sigma',nvrt-kz+1, nsigma_dim)
    if(kz/=1) iret = nf90_def_dim(ncid, 'nz',kz-1, nz_dim)
    iret = nf90_def_dim(ncid, 'time', NF90_UNLIMITED, ntime_dim)
      
!    coords_type =1 ! default project coords
      
    if(ics==1) then
      lat_coord_standard_name = "projection_y_coordinate"
      lon_coord_standard_name = "projection_x_coordinate"
      x_units = "m"
      y_units = "m"
      lat_str_len = 23
      lon_str_len = 23
    else
      lat_coord_standard_name = "latitude"
      lon_coord_standard_name = "longitude"
      y_units = "degrees_north"
      x_units = "degrees_east"
      lat_str_len = 8
      lon_str_len = 9
    endif
      
    !define variables
    time_dims(1) = ntime_dim
    iret=nf90_def_var(ncid,'time',NF90_DOUBLE,time_dims,itime_id)
    if(iret.ne.NF90_NOERR) then
      print*, nf90_strerror(iret); stop
    endif

    write(start_time,'(i5,2(1x,i2),2(1x,f10.2))')start_year,start_month,start_day,start_hour,utc_start

    iret=nf90_put_att(ncid,itime_id,'long_name','Time')
!  The time coordinate variable needs a units attribute compliant with udunits2 to be CF compliant and
!  for most of the NetCDF tools to be able to interpret the time coordinate. 
    write(cbuffer,20) start_year,start_month,start_day,int(start_hour),-int(utc_start*100)
20  FORMAT('seconds since ',I0.4,'-',I0.2,'-',I0.2,' ',I0.2,':00:00 ',SP,I0.4)

    iret=nf90_put_att(ncid,itime_id,'units',cbuffer)
    iret=nf90_put_att(ncid,itime_id,'base_date',start_time) !len(start_time),start_time)
    iret=nf90_put_att(ncid,itime_id,'standard_name','time')
    !     write time stamps later

    !define mesh
    !Error: not sure what dummy_dim is?
    iret=nf90_def_var(ncid,'SCHISM_hgrid',NF90_INT,one_dim,ihgrid_id) !0,dummy_dim,ihgrid_id)
    iret=nf90_put_att(ncid,ihgrid_id,'long_name','Topology data of 2d unstructured mesh')
!    tempint_array(1)=2
    iret=nf90_put_att(ncid,ihgrid_id,'topology_dimension',2) !tempint_array)
    iret=nf90_put_att(ncid,ihgrid_id,'cf_role','mesh_topology')
    iret=nf90_put_att(ncid,ihgrid_id,'node_coordinates','SCHISM_hgrid_node_x SCHISM_hgrid_node_y')
    iret=nf90_put_att(ncid,ihgrid_id,'face_node_connectivity','SCHISM_hgrid_face_nodes')
    iret=nf90_put_att(ncid,ihgrid_id,'edge_coordinates','SCHISM_hgrid_edge_x SCHISM_hgrid_edge_y')
    iret=nf90_put_att(ncid,ihgrid_id,'face_coordinates','SCHISM_hgrid_face_x SCHISM_hgrid_face_y')
    iret=nf90_put_att(ncid,ihgrid_id,'edge_node_connectivity','SCHISM_hgrid_edge_nodes')

!    iret=nf90_def_var(ncid,'Mesh3D',NF90_INT,one_dim,i3Dhgrid_id) !0,dummy_dim,i3Dhgrid_id)
!    iret=nf90_put_att(ncid,i3Dhgrid_id,'long_name','Topology data of 3d unstructured mesh')
!!    tempint_array(1)=3
!    iret=nf90_put_att(ncid,i3Dhgrid_id,'topology_dimension',3) !tempint_array)
!    iret=nf90_put_att(ncid,i3Dhgrid_id,'cf_role','mesh_topology')
     
    ele_dims(2)=nele_dim; ele_dims(1)=nfour_dim
    iret=nf90_def_var(ncid,'SCHISM_hgrid_face_nodes',NF90_INT,ele_dims,iele_id)
    iret=nf90_put_att(ncid,iele_id,'long_name','Horizontal Element Table')
!   iret=nf90_put_att(ncid,iele_id,'units','non-dimensional')
    iret=nf90_put_att(ncid,iele_id,'cf_role','face_node_connectivity')
!    int_buffer(1) = 1  ! fortran indexing starts at 1
    iret=nf90_put_att(ncid,iele_id,'start_index',1) !int_buffer)
!    int_buffer(1) = invalid_index
    iret=nf90_put_att(ncid,iele_id,'_FillValue',invalid_index) !int_buffer)
      
!    iret = nf90_enddef(ncid) 
!    if(iret.ne.NF90_NOERR) then
!      print*, 'WOW2:',nf90_strerror(iret); stop
!    endif
!    stop

    ele_dims(2)=nedge_dim; ele_dims(1)=ntwo_dim
    iret=nf90_def_var(ncid,'SCHISM_hgrid_edge_nodes',NF90_INT,ele_dims,iedge_id)
    iret=nf90_put_att(ncid,iedge_id,'long_name','Map every edge to the two nodes that it connects')
!   iret=nf90_put_att(ncid,iedge_id,'units','non-dimensional')
    iret=nf90_put_att(ncid,iedge_id,'cf_role','edge_node_connectivity')
!    int_buffer(1) = 1  ! fortran indexing starts at 1
    iret=nf90_put_att(ncid,iedge_id,'start_index',1) !int_buffer)
      
    x_dims(1)=node_dim
    iret=nf90_def_var(ncid,'SCHISM_hgrid_node_x',NF90_FLOAT,x_dims,ix_id)
    iret=nf90_put_att(ncid,ix_id,'long_name','node x-coordinate')
    iret=nf90_put_att(ncid,ix_id,'standard_name',lon_coord_standard_name)
    iret=nf90_put_att(ncid,ix_id,'units',x_units)
    iret=nf90_put_att(ncid,ix_id,'mesh','SCHISM_hgrid')
      
    iret=nf90_def_var(ncid,'SCHISM_hgrid_node_y',NF90_FLOAT,x_dims,iy_id)
    iret=nf90_put_att(ncid,iy_id,'long_name','node y-coordinate')
    iret=nf90_put_att(ncid,iy_id,'standard_name',lat_coord_standard_name)
    iret=nf90_put_att(ncid,iy_id,'units',y_units)
    iret=nf90_put_att(ncid,iy_id,'mesh','SCHISM_hgrid')
      
    iret=nf90_def_var(ncid,'node_bottom_index',NF90_INT,x_dims,inode_bottom_id)
    iret=nf90_put_att(ncid,inode_bottom_id,'long_name','bottom level index at each node')
    iret=nf90_put_att(ncid,inode_bottom_id,'units','non-dimensional')
    iret=nf90_put_att(ncid,inode_bottom_id,'mesh','SCHISM_hgrid')
    iret=nf90_put_att(ncid,inode_bottom_id,'location','node')
!    int_buffer(1) = 1  
    iret=nf90_put_att(ncid,inode_bottom_id,'start_index',1) !int_buffer)
      
    x_dims(1)=nele_dim
    iret=nf90_def_var(ncid,'SCHISM_hgrid_face_x',NF90_FLOAT,x_dims,ix_face_id)
    iret=nf90_put_att(ncid,ix_face_id,'long_name','x_coordinate of 2D mesh face')
    iret=nf90_put_att(ncid,ix_face_id,'standard_name',lon_coord_standard_name)
    iret=nf90_put_att(ncid,ix_face_id,'units',x_units)
    iret=nf90_put_att(ncid,ix_face_id,'mesh','SCHISM_hgrid')

    iret=nf90_def_var(ncid,'SCHISM_hgrid_face_y',NF90_FLOAT,x_dims,iy_face_id)
    iret=nf90_put_att(ncid,iy_face_id,'long_name','y_coordinate of 2D mesh face')
    iret=nf90_put_att(ncid,iy_face_id,'standard_name',lat_coord_standard_name)
    iret=nf90_put_att(ncid,iy_face_id,'units',y_units)
    iret=nf90_put_att(ncid,iy_face_id,'mesh','SCHISM_hgrid')
      
    iret=nf90_def_var(ncid,'ele_bottom_index',NF90_INT,x_dims,iele_bottom_id)
    iret=nf90_put_att(ncid,iele_bottom_id,'long_name','bottom level index at each element')
    iret=nf90_put_att(ncid,iele_bottom_id,'units','non-dimensional')
    iret=nf90_put_att(ncid,iele_bottom_id,'mesh','SCHISM_hgrid')
    iret=nf90_put_att(ncid,iele_bottom_id,'location','elem')
!    int_buffer(1) = 1  
    iret=nf90_put_att(ncid,iele_bottom_id,'start_index',1) !int_buffer)
      
    x_dims(1)=nedge_dim
    iret=nf90_def_var(ncid,'SCHISM_hgrid_edge_x',NF90_FLOAT,x_dims,ix_edge_id)
    iret=nf90_put_att(ncid,ix_edge_id,'long_name','x_coordinate of 2D mesh edge')
    iret=nf90_put_att(ncid,ix_edge_id,'standard_name',lon_coord_standard_name)
    iret=nf90_put_att(ncid,ix_edge_id,'units',x_units)
    iret=nf90_put_att(ncid,ix_edge_id,'mesh','SCHISM_hgrid')

    iret=nf90_def_var(ncid,'SCHISM_hgrid_edge_y',NF90_FLOAT,x_dims,iy_edge_id)
    iret=nf90_put_att(ncid,iy_edge_id,'long_name','y_coordinate of 2D mesh edge')
    iret=nf90_put_att(ncid,iy_edge_id,'standard_name',lat_coord_standard_name)
    iret=nf90_put_att(ncid,iy_edge_id,'units',y_units)
    iret=nf90_put_att(ncid,iy_edge_id,'mesh','SCHISM_hgrid')
    
    iret=nf90_def_var(ncid,'edge_bottom_index',NF90_INT,x_dims,iedge_bottom_id)
    iret=nf90_put_att(ncid,iedge_bottom_id,'long_name','bottom level index at each edge')
    iret=nf90_put_att(ncid,iedge_bottom_id,'units','non-dimensional')
    iret=nf90_put_att(ncid,iedge_bottom_id,'mesh','SCHISM_hgrid')
    iret=nf90_put_att(ncid,iedge_bottom_id,'location','edge')
!    int_buffer(1) = 1  
    iret=nf90_put_att(ncid,iedge_bottom_id,'start_index',1) !int_buffer)
      
    x_dims(1)=node_dim
    iret=nf90_def_var(ncid,'depth',NF90_FLOAT,x_dims,idepth_id)
    iret=nf90_put_att(ncid,idepth_id,'long_name','Bathymetry')
    iret=nf90_put_att(ncid,idepth_id,'units','meters')
    iret=nf90_put_att(ncid,idepth_id,'positive','down')
    iret=nf90_put_att(ncid,idepth_id,'mesh','SCHISM_hgrid')
    iret=nf90_put_att(ncid,idepth_id,'location','node')
      
    sigma_dims(1)=nsigma_dim
    ! See http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.7-draft1/apd.html 
    ! section "Ocean s-coordinate"  for the CF calculation
    ! See trunk/src/Utility/Post-Processing-Fortran/compute_zcor.f90 for the SCHISM calc.
    ! CF:SCHISM corresponence: depth_c:h_c s:sigma C(k):cs(k) a:theta_f b:theta_b 
    hs_array(1)=h_s
    hc_array(1)=h_c
    thetab_array(1)=theta_b
    thetaf_array(1)=theta_f
    iret=nf90_def_var(ncid,'sigma',NF90_FLOAT,sigma_dims,isigma_id)
    iret=nf90_put_att(ncid,isigma_id,'long_name','S coordinates at whole levels')
    iret=nf90_put_att(ncid,isigma_id,'units','1')
    cbuffer='ocean_s_coordinate'
    iret=nf90_put_att(ncid,isigma_id,'standard_name',cbuffer)
    iret=nf90_put_att(ncid,isigma_id,'positive','up')
    iret=nf90_put_att(ncid,isigma_id,'h_s',h_s) !hs_array)
    iret=nf90_put_att(ncid,isigma_id,'h_c',h_c) !hc_array)
    iret=nf90_put_att(ncid,isigma_id,'theta_b',theta_b) !thetab_array)
    iret=nf90_put_att(ncid,isigma_id,'theta_f',theta_f) !thetaf_array)
    cbuffer='s: sigma eta: elev depth: depth a: sigma_theta_f b: sigma_theta_b depth_c: sigma_h_c'
    iret=nf90_put_att(ncid,isigma_id,'formula_terms',cbuffer)

    var1d_dims(1)=one_dim
    iret=nf90_def_var(ncid,'dry_value_flag',NF90_INT,var1d_dims,iwetdry_id)
    iret=nf90_put_att(ncid,iwetdry_id,'values','0: use last-wet value; 1: use junk')
    iret=nf90_def_var(ncid,'coordinate_system_flag',NF90_INT,var1d_dims,icoord)
    iret=nf90_def_var(ncid,'minimum_depth',NF90_FLOAT,var1d_dims,ih0)
    iret=nf90_def_var(ncid,'sigma_h_c',NF90_FLOAT,var1d_dims,ihc_id)
    if(iret.ne.NF90_NOERR) then
      print*, nf90_strerror(iret); stop
    endif
    cbuffer='ocean_s_coordinate h_c constant'
    iret=nf90_put_att(ncid,ihc_id,'long_name',cbuffer)
    iret=nf90_put_att(ncid,ihc_id,'units','meters')
    iret=nf90_put_att(ncid,ihc_id,'positive','down')

    var1d_dims(1)=one_dim
    iret=nf90_def_var(ncid,'sigma_theta_b',NF90_FLOAT,var1d_dims,itheta_b_id)
    cbuffer='ocean_s_coordinate theta_b constant'
    iret=nf90_put_att(ncid,itheta_b_id,'long_name',cbuffer)

    var1d_dims(1)=one_dim
    iret=nf90_def_var(ncid,'sigma_theta_f',NF90_FLOAT,var1d_dims,itheta_f_id)
    cbuffer='ocean_s_coordinate theta_f constant'
    iret=nf90_put_att(ncid,itheta_f_id,'long_name',cbuffer)

    var1d_dims(1)=one_dim
    iret=nf90_def_var(ncid,'sigma_maxdepth',NF90_FLOAT,var1d_dims,ihs_id)
    cbuffer='ocean_s_coordinate maximum depth cutoff (mixed s over z boundary)'
    iret=nf90_put_att(ncid,ihs_id,'long_name',cbuffer)
    iret=nf90_put_att(ncid,ihs_id,'units','meters')
    iret=nf90_put_att(ncid,ihs_id,'positive','down')
       
    iret=nf90_def_var(ncid,'Cs',NF90_FLOAT,sigma_dims,ics_id)
    iret=nf90_put_att(ncid,ics_id,'long_name','Function C(s) at whole levels')
!   iret=nf90_put_att(ncid,ics_id,'units','non-dimensional')
    iret=nf90_put_att(ncid,ics_id,'positive','up')

    if(kz/=1) then !(kz-1) z levels 
      z_dims(1)=nz_dim
      iret=nf90_def_var(ncid,'z',NF90_FLOAT,z_dims,iz_id)
      iret=nf90_put_att(ncid,iz_id,'long_name','Z coordinates at whole levels')
      iret=nf90_put_att(ncid,iz_id,'units','meters')
      iret=nf90_put_att(ncid,iz_id,'positive','up')
    endif
      
!          call get_netcdf_var_names(file63(1:lfile63-3),&
!                                    &file63((lfile63-1):lfile63),&
!                                    &netcdf_var_name,&
!                                    &netcdf_var_long_name,&
!                                    &netcdf_var_standard_name,&
!                                    &netcdf_out_location,&
!                                    &netcdf_level_location,&
!                                    &netcdf_var_dim,&
!                                    &found_netcdf_var)
!          if(.not.found_netcdf_var) then
!              print *, "unsupported data type :", file63(1:lfile63-3)
!              exit
!          endif

    !Find and define variables
    file63='schout_0000_'//it_char(1:it_len)//'.nc'
    file63=adjustl(file63)
    iret=nf90_open(trim(file63),OR(NF90_NETCDF4,NF90_NOWRITE),ncid2)
    !iret=nf_inq_nvars(ncid2,nvars)
    iret=nf90_inquire(ncid2,nVariables=nvars)
!    write(99,*)'nvars=',nvars,file63

    if(iinput==ibgn) then 
      allocate(i23d(nvars),ivs(nvars),variable_nm(nvars),vlen(nvars),ndims(nvars), &
    &iu_id(nvars))
      allocate(outd(nvars),stat=istat)
      allocate(skip_var(nvars),stat=istat)
      skip_var(:)=.false. ! by default do not skip variables
      check_vars = len_trim(to_be_combined)>0
    endif

    !Debug
    !var2d_dims(1)=node_dim; var2d_dims(2)=ntime_dim
    !iret=nf_def_var(ncid,'elev',NF_FLOAT,2,var2d_dims,iu_id(2))
    !write(98,*)'Before:',ncid,var2d_dims,iu_id(2),nvars,iinput,ibgn

    do m=1,nvars
      iret=nf90_inquire_variable(ncid2,m,name=variable_nm(m)) !,itype,ndims(m),int_buffer,natts)
      variable_nm(m)=trim(adjustl(variable_nm(m))); vlen(m)=len_trim(variable_nm(m))
      if (check_vars) then
        if ((index(to_be_combined,trim(variable_nm(m)))>0).or. &
           (index(default_variables,trim(variable_nm(m)))>0)) then
          skip_var(m)=.false.
        else
          skip_var(m)=.true.
          cycle ! continue with next variable in loop
        end if
        if (m==nvars) check_vars=.false. ! continue with skip_var array
      else
        if (skip_var(m)) cycle
      end if

      iret=nf90_get_att(ncid2,m,'i23d',i23d(m))

!      write(99,*)'i23d:',m,variable_nm(m),i23d(m)
         
      if(i23d(m)>0) then !not time array
        iret=nf90_get_att(ncid2,m,'ivs',ivs(m))
        if(i23d(m)<=3) then 
          netcdf_out_location="node"
        else if(i23d(m)<=6) then 
          netcdf_out_location="elem"
        else if(i23d(m)<=9) then 
          netcdf_out_location="side"
        else
          stop 'Unknown i23d'
        endif !i23d
        netcdf_level_location="full"
        if(mod(i23d(m),3)==0) netcdf_level_location="half"

        if(i23d(m)<=3) then
          var2d_dims(1)=node_dim
          var3d_dims(2)=node_dim
          var4d_dims(3)=node_dim
          !chunks(1)=np_global
          npse=np_global
        else if(i23d(m)<=6) then 
          var2d_dims(1)=nele_dim
          var3d_dims(2)=nele_dim
          var4d_dims(3)=nele_dim
          !chunks(1)=ne_global  
          npse=ne_global  
        else
          var2d_dims(1)=nedge_dim
          var3d_dims(2)=nedge_dim
          var4d_dims(3)=nedge_dim
          !chunks(1)=ns_global
          npse=ns_global
        endif

        ! allocate output arrays
        ! todo: use ivs to control array size here and
        !       for copying into global arrays/masking later
        if (iinput==ibgn) then
          select case (i23d(m))
          case (1)
            allocate(outd(m)%data(2,1,np_global))
          case (2)
            allocate(outd(m)%data(2,nvrt,np_global))
          case (3)
            allocate(outd(m)%data(2,nvrt,np_global))
          case (4)
            allocate(outd(m)%data(2,1,ne_global))
          case (5)
            allocate(outd(m)%data(2,nvrt,ne_global))
          case (6)
            allocate(outd(m)%data(2,nvrt,ne_global))
          case (7)
            allocate(outd(m)%data(2,1,ns2))
          case (8)
            allocate(outd(m)%data(2,nvrt,ns2))
          case (9)
            allocate(outd(m)%data(2,nvrt,ns2))
          end select
          if (associated(outd(m)%data)) outd(m)%data = -9999.0d0
        end if

        if(mod(i23d(m)-1,3)==0) then !2D
          if(ivs(m)==1) then
            var2d_dims(2)=ntime_dim
            iret=nf90_def_var(ncid,variable_nm(m),NF90_FLOAT,var2d_dims,iu_id(m))
          else !2D vector
            var3d_dims(1)=ntwo_dim; var3d_dims(3)=ntime_dim
            iret=nf90_def_var(ncid,variable_nm(m),NF90_FLOAT,var3d_dims,iu_id(m))
          endif !ivs
        else !3D
          if(ivs(m)==1) then
            var3d_dims(1)=nv_dim; var3d_dims(3)=ntime_dim
            chunks(1)=nvrt; chunks(2)=npse; chunks(3)=1
            iret=nf90_def_var(ncid,variable_nm(m),NF90_FLOAT,var3d_dims,iu_id(m))
!#ifdef NETCDF_4
            iret=nf90_def_var_chunking(ncid,iu_id(m),NF90_CHUNKED,chunks)
!function nf90_def_var_deflate(ncid, varid, shuffle, deflate, deflate_level)
!where deflate_level\in[0,9] with 9 being most compression
            iret=nf90_def_var_deflate(ncid,iu_id(m),0,1,4)
!#endif
          else !3D vector
            var4d_dims(2)=nv_dim; var4d_dims(1)=ntwo_dim; var4d_dims(4)=ntime_dim
            chunks(2)=nvrt; chunks(1)=2; chunks(3)=npse
            iret=nf90_def_var(ncid,variable_nm(m),NF90_FLOAT,var4d_dims,iu_id(m))
!#ifdef NETCDF_4
            iret=nf90_def_var_chunking(ncid,iu_id(m),NF90_CHUNKED, chunks)
            iret=nf90_def_var_deflate(ncid,iu_id(m),0,1,4)
!#endif
          endif !ivs
        endif !i23d: 2|3D
! debug header information in stderr:
!write(0,*) m,trim(variable_nm(m)),chunks,iret

!       iret=nf_put_att_text(ncid,iu_id,'long_name',len_trim(netcdf_var_long_name(1)),netcdf_var_long_name(1))
!       iret=nf_put_att_text(ncid,iu_id,'standard_name',len_trim(netcdf_var_standard_name(1)),netcdf_var_standard_name(1))
        iret=nf90_put_att(ncid,iu_id(m),'missing_value',NF90_FILL_FLOAT)
        iret=nf90_put_att(ncid,iu_id(m),'mesh','SCHISM_hgrid')
        iret=nf90_put_att(ncid,iu_id(m),'data_horizontal_center',netcdf_out_location)
        iret=nf90_put_att(ncid,iu_id(m),'data_vertical_center',netcdf_level_location)
        iret=nf90_put_att(ncid,iu_id(m),'i23d',i23d(m)) 
        iret=nf90_put_att(ncid,iu_id(m),'ivs',ivs(m)) 
      endif !i23d(m)>0
    enddo !m=1,nvars

!   global attributes per
!   http://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html#_attributes
    iret = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions','CF-1.0, UGRID-1.0')
    iret = nf90_put_att(ncid, NF90_GLOBAL, 'title','SCHISM Model output')
    iret = nf90_put_att(ncid, NF90_GLOBAL, 'institution','SCHISM Model output')
    write(cbuffer,30) version
30  FORMAT('SCHISM model output version ',A48)
    iret = nf90_put_att(ncid, NF90_GLOBAL, 'source',cbuffer)
    iret = nf90_put_att(ncid, NF90_GLOBAL, 'references','http://ccrm.vims.edu/schismweb/')
    iret = nf90_put_att(ncid, NF90_GLOBAL, 'history','created by combine_output11')
    iret = nf90_put_att(ncid, NF90_GLOBAL, 'comment','SCHISM Model output')

! Extra global attributes
    iret = nf90_put_att(ncid, NF90_GLOBAL, 'type','SCHISM Model output')
    iret = nf90_put_att(ncid, NF90_GLOBAL, 'VisIT_plugin','https://schism.water.ca.gov/library/-/document_library/view/3476283') 

    !leave define mode
    iret = nf90_enddef(ncid)

    !get nrec from time dimension
    !write(0,*) 'Number of time records from local_to_global: ',nrec
    iret = nf90_inquire(ncid2, unlimitedDimId=dimid)
    iret = nf90_inquire_dimension(ncid2, dimid, len=nrec)
    !write(0,*) 'Number of time records from netcdf: ',nrec

    !done with reading header
    iret=nf90_close(ncid2)

!    write(99,*)'out of def mode'

    !Write mode (static part only)
    data_start_2d(1:2)=1
    data_count_2d(1)=4; data_count_2d(2)=ne_global
    iret=nf90_put_var(ncid,iele_id,elnode,data_start_2d,data_count_2d)
    if(iret.ne.NF90_NOERR) then
      print*, nf90_strerror(iret); stop
    endif
      
    data_start_2d(1:2)=1
    data_count_2d(1)=2; data_count_2d(2)=ns_global
    iret=nf90_put_var(ncid,iedge_id,isidenode,data_start_2d,data_count_2d)
    if(iret.ne.NF90_NOERR) then
      print*, nf90_strerror(iret); stop
    endif
 
    ! fill bottom index
!    data_start_1d(1)=1; data_count_1d(1)=np_global
    iret=nf90_put_var(ncid,inode_bottom_id,kbp00,(/1/),(/np_global/))
!    data_count_1d(1)=ns_global
    iret=nf90_put_var(ncid,iedge_bottom_id,kbs,(/1/),(/ns_global/))
!    data_count_1d(1)=ne_global
    iret=nf90_put_var(ncid,iele_bottom_id,kbe,(/1/),(/ne_global/))
  
    ! fill node,side center coords
    iret=nf90_put_var(ncid,ix_id,x,(/1/),(/np_global/))
    iret=nf90_put_var(ncid,iy_id,y,(/1/),(/np_global/))
    iret=nf90_put_var(ncid,ix_edge_id,xcj,(/1/),(/ns_global/))
    iret=nf90_put_var(ncid,iy_edge_id,ycj,(/1/),(/ns_global/))
    iret=nf90_put_var(ncid,ix_face_id,xctr,(/1/),(/ne_global/))
    iret=nf90_put_var(ncid,iy_face_id,yctr,(/1/),(/ne_global/))
          
    iret=nf90_put_var(ncid,idepth_id,dp,(/1/),(/np_global/))
    iret=nf90_put_var(ncid,isigma_id,sigma,(/1/),(/nvrt-kz+1/))
    iret=nf90_put_var(ncid,ics_id,cs,(/1/),(/nvrt-kz+1/))
    iret=nf90_put_var(ncid,iwetdry_id,iwetdry)
    iret=nf90_put_var(ncid,icoord,ics) !,(/1/),(/1/))
    iret=nf90_put_var(ncid,ih0,h0) !,(/1/),(/1/))
    iret=nf90_put_var(ncid,ihc_id,h_c) !,(/1/),(/1/))
    iret=nf90_put_var(ncid,itheta_b_id,theta_b) !,(/1/),(/1/))
    iret=nf90_put_var(ncid,itheta_f_id,theta_f) !,(/1/),(/1/))
    iret=nf90_put_var(ncid,ihs_id,h_s) !,(/1/),(/1/))
    if(kz/=1) iret=nf90_put_var(ncid,iz_id,ztot,(/1/),(/kz-1/))
      
    !print*, 'Last element:',elnode(1:3,ne_global)
    !end output header

    ! Loop over records in file
    do ispool=1,nrec
      !Gather all ranks
      do irank=0,nproc-1
        !Open input file
        write(a_4,'(i4.4)') irank
        file63='schout_'//a_4//'_'//it_char(1:it_len)//'.nc'
        file63=adjustl(file63)
        iret=nf90_open(trim(file63),OR(NF90_NETCDF4,NF90_NOWRITE),ncid2)
        !write(99,*)'nvars=',nvars,file63

        do m=1,nvars
          if (skip_var(m)) cycle 
          ! get variable id (maybe not necessary)
          iret = nf90_inq_varid(ncid2,variable_nm(m),varid)
          !if (iret/=NF90_NOERR) then
          !  write(0,*) trim(file63),trim(nf90_strerror(iret))
          !  stop
          !endif
          !! alternatively use same id as in first file:
          !varid = m

          if(i23d(m)<=3) then !node
            data_count_2d(1)=np(irank)
            data_count_3d(2)=np(irank)
            data_count_4d(3)=np(irank)
          else if(i23d(m)<=6) then
            data_count_2d(1)=ne(irank)
            data_count_3d(2)=ne(irank)
            data_count_4d(3)=ne(irank)
          else
            data_count_2d(1)=ns(irank)
            data_count_3d(2)=ns(irank)
            data_count_4d(3)=ns(irank)
          endif

          if(i23d(m)==0) then !time
            data_start_1d(1)=ispool; data_count_1d(1)=1
            iret=nf90_get_var(ncid2,varid,hc_array,data_start_1d,data_count_1d) 
            time=hc_array(1)
            !write(99,*)'time=',time
          else if(mod(i23d(m)-1,3)==0) then !2D
            if(ivs(m)==1) then
              data_start_2d(1)=1; data_start_2d(2)=ispool 
              data_count_2d(2)=1
              iret=nf90_get_var(ncid2,varid,worka(1,1,1:data_count_2d(1)),data_start_2d,data_count_2d)
            else !vector
              data_start_3d(1:2)=1; data_start_3d(3)=ispool 
              data_count_3d(1)=2; data_count_3d(3)=1
              iret=nf90_get_var(ncid2,varid,worka(1:2,1,1:data_count_3d(2)),data_start_3d,data_count_3d)
            endif !ivs
          else !3D
            if(ivs(m)==1) then
              data_start_3d(1:2)=1; data_start_3d(3)=ispool
              data_count_3d(1)=nvrt; data_count_3d(3)=1
              iret=nf90_get_var(ncid2,varid,worka(1,:,1:data_count_3d(2)),data_start_3d,data_count_3d)
            else !vector
              data_start_4d(1:3)=1; data_start_4d(4)=ispool
              data_count_4d(1)=2; data_count_4d(2)=nvrt; data_count_4d(4)=1
              iret=nf90_get_var(ncid2,varid,worka(1:2,:,1:data_count_4d(3)),data_start_4d,data_count_4d)
            endif !ivs
          endif !i23d

          !Put into global index
          if(i23d(m)==0) then !do nothing
          else if(i23d(m)<=3) then !node
            do i=1,np(irank)
              nd=iplg(irank,i)
              outd(m)%data(1:2,:,nd)=worka(1:2,1:ubound(outd(m)%data,2),i)
            enddo !i
          else if(i23d(m)<=6) then
            do i=1,ne(irank)
              ie=ielg(irank,i)
              outd(m)%data(1:2,:,ie)=worka(1:2,1:ubound(outd(m)%data,2),i)
            enddo !i
          else
            do i=1,ns(irank)
              isd=islg(irank,i)
              outd(m)%data(1:2,:,isd)=worka(1:2,1:ubound(outd(m)%data,2),i)
            enddo !i
          endif !i23d
        enddo !m=1,nvars
        iret=nf90_close(ncid2)
      enddo !irank

      do m=1,nvars
        if (skip_var(m)) cycle
        !Compute wet/dry flags
        if(m==2) then !idry_e
          idry_e=nint(outd(m)%data(1,1,1:ne_global))
          idry=1 !init as dry
          idry_s=1
          do i=1,ne_global
            do j=1,i34(i)
              nd=elnode(j,i)
              isd=elside(j,i)
              if(idry_e(i)==0) then
                idry(nd)=0
                idry_s(isd)=0
              endif
            enddo !j
          enddo !i
        endif !m=2

!        !Compute bottom indices based on zcor
!        if(m==3) then !zcor
!          do i=1,np_global
!            kbp00(i)=nvrt+1 !init for dry
!            do k=1,nvrt
!              if(outd(m)%data(1,k,i)>-1.e19) then
!                kbp00(i)=k; exit
!              endif
!            enddo !k
!          enddo !i
!
!          do i=1,ne_global
!            kmax=maxval(kbp00(elnode(1:i34(i),i)))
!            if(kmax==nvrt+1) then !dry
!              kbe(i)=nvrt+1
!            else !wet
!              kbe(i)=minval(kbp00(elnode(1:i34(i),i)))
!            endif
!          enddo !i
!
!          do i=1,ns_global
!            kmax=maxval(kbp00(isidenode(1:2,i)))
!            if(kmax==nvrt+1) then !dry
!              kbs(i)=nvrt+1
!            else !wet
!              kbs(i)=minval(kbp00(isidenode(1:2,i)))
!            endif
!          enddo !i
!        endif !zcor
!
        !Fill below-bottom for 3D vars
        if(mod(i23d(m)-1,3)/=0.and.i23d(m)>0) then !3D
          if(i23d(m)<=3) then !node
            do i=1,np_global
              outd(m)%data(1:2,1:kbp00(i)-1,i)=NF90_FILL_FLOAT
              if(idry(i)==1.and.iwetdry/=0) outd(m)%data(1:2,:,i)=NF90_FILL_FLOAT
            enddo !i
          else if(i23d(m)<=6) then
            do i=1,ne_global
              outd(m)%data(1:2,1:kbe(i)-1,i)=NF90_FILL_FLOAT
              if(idry_e(i)==1.and.iwetdry/=0) outd(m)%data(1:2,:,i)=NF90_FILL_FLOAT
            enddo !i
          else
            do i=1,ns_global
              outd(m)%data(1:2,1:kbs(i)-1,i)=NF90_FILL_FLOAT
              if(idry_s(i)==1.and.iwetdry/=0) outd(m)%data(1:2,:,i)=NF90_FILL_FLOAT
            enddo !i
          endif !i23d
        endif !3D

        !Debug
!        if(i23d(m)>0) then
!          write(98,*)'ispool=',ispool,iinput,m,i23d(m),ivs(m),iu_id(m),variable_nm(m)
!          do i=1,np_global
!            if(mod(i23d(m)-1,3)==0) then
!              write(98,*)i,(1:ivs(m),1,i)
!            else
!              write(98,*)i,(1:ivs(m),nvrt,i)
!            endif
!          enddo !i
!        endif !i23d

        !write combined arrays
        if(i23d(m)<=3) then
          data_count_2d(1)=np_global
          data_count_3d(2)=np_global
          data_count_4d(3)=np_global
        else if(i23d(m)<=6) then 
          data_count_2d(1)=ne_global
          data_count_3d(2)=ne_global
          data_count_4d(3)=ne_global
        else
          data_count_2d(1)=ns_global
          data_count_3d(2)=ns_global
          data_count_4d(3)=ns_global
        endif !i23d

        if(i23d(m)==0) then !time
          hc_array(1)=time
          iret=nf90_put_var(ncid,itime_id,dble(hc_array),(/ispool/),(/1/))
        else if(mod(i23d(m)-1,3)==0) then !2D
          if(ivs(m)==1) then
            data_start_2d(1)=1; data_start_2d(2)=ispool
            data_count_2d(2)=1
            iret=nf90_put_var(ncid,iu_id(m),outd(m)%data(1,1,1:data_count_2d(1)),data_start_2d,data_count_2d)
          else
            data_start_3d(1:2)=1; data_start_3d(3)=ispool
            data_count_3d(1)=2; data_count_3d(3)=1
            iret=nf90_put_var(ncid,iu_id(m),outd(m)%data(1:2,1,1:data_count_3d(2)),data_start_3d,data_count_3d)
          endif !ivs
        else !3D
          if(ivs(m)==1) then
            data_start_3d(1:2)=1; data_start_3d(3)=ispool
            data_count_3d(1)=nvrt; data_count_3d(3)=1
            iret=nf90_put_var(ncid,iu_id(m),outd(m)%data(1,:,1:data_count_3d(2)),data_start_3d,data_count_3d)
          else
            data_start_4d(1:3)=1; data_start_4d(4)=ispool
            data_count_4d(1)=2; data_count_4d(2)=nvrt; data_count_4d(4)=1
            iret=nf90_put_var(ncid,iu_id(m),outd(m)%data(:,:,1:data_count_4d(3)),data_start_4d,data_count_4d)
          endif
        endif !i23d
      enddo !m=1,nvars
      iret = nf90_sync(ncid)
    enddo !ispool=1,nrec
    ! Close output file
    iret = nf90_close(ncid)
  enddo ! do iinput=ibgn,iend
  
  deallocate(x,y,dp,kbp00,kbe,np,ns,ne,elnode,nm2,ztot,sigma,cs,ihot_len, &
  &dpe,xctr,yctr,iplg,ielg,islg,ic3,elside,isdel,isidenode, &
  &xcj,ycj,dps,kbs,idry,idry_s,idry_e)

  ! deallocate output arrays
  do m=1,size(outd)
    if (skip_var(m)) cycle 
    if (associated(outd(m)%data)) deallocate(outd(m)%data)
    nullify(outd(m)%data)
  end do
  deallocate(outd)
  
end subroutine combine_output11

subroutine combine_output11_input(ibgn,iend,iwetdry,to_be_combined,output_prefix)
use argparse
implicit none
integer, parameter                   :: nfilemax = 20
!integer                              :: nfile   !< number of file base names
!character(len=30),dimension(nfilemax) :: files   !< base names (e.g. elev.61) of files to combine
integer,intent(out)                              :: ibgn    !< first output file index to process
integer,intent(out)                              :: iend    !< last output file index to process
integer,intent(out)                              :: iwetdry !dry value option
!integer                              :: inetcdf !< netcdf flag
character(len=80) :: infile
character(len=30) :: cfile = ""
character(len=1024) :: varlist = ""
character(len=1024) :: output_prefix
character(len=1024) :: to_be_combined ! list of variables to be combined

! local
integer :: comcount
integer :: ifile

!infile = ""
cfile = ""
!files=""

cmd_name = "combine_output11"
call cla_init(cmd_name,"Combine time blocked per-processor binary outputs (e.g. 'schout_0000_1.nc') into time blocked global outputs ('schout_1.nc')")

!call cla_register('-i','--in', 'input file (e.g. combine_input.in) containing options (overridden by command line specs)', cla_char,'')
call cla_register('-b','--begin', 'start day', cla_int,'-1')
call cla_register('-e','--end','end day', cla_int  ,'-1')
call cla_register('-w','--wetdry','dry option (0: last wet value; 1: junk value)', cla_int  ,'0')
!call cla_register('-n','--nc','combine to NetCDF format (1) or ordinary binary (0)', cla_int  ,'0')
!call cla_register('-f','--file','base file name like elev.61',  cla_char,'') 
call cla_register('-v','--vars','comma separated list of variables',  cla_char,'') 
call cla_register('-o','--output','output file prefix',  cla_char,'schout') 
call cla_validate
    
comcount = command_argument_count()
if(comcount == 0) then
  stop './"combine_output11" -h for help'
  !infile = "combine_output.in"
!else
!  call cla_get("--in",infile)
end if

! possibly read files from input file if there is one
!if (len_trim(infile)>1) then
!  ! read inputs from combine_output.in or other provided file
!  open(10,file=infile,status='old')
!  read(10,'(a12)') cfile !e.g. 'hvel.64'
!  read(10,*) ibgn,iend  
!  read(10,*) inetcdf       ! netcdf option
!  close(10)
!! command line takes precedence
!  if (nfile == 0) then
!    nfile = 1
!    files(1) = trim(cfile)
!  end if 
!  if (cla_key_present("--begin"))call cla_get("--begin",ibgn)
!  if (cla_key_present("--end"))  call cla_get("--end",iend)
!  if (cla_key_present("--nc"))   call cla_get("--nc",inetcdf)
!  if (cla_key_present("--file")) call cla_get("--file",cfile)

!else
  call cla_get("--begin",ibgn)
  call cla_get("--end",iend)
  call cla_get("--wetdry",iwetdry)
  call cla_get("--vars", varlist)
  call cla_get("--output", output_prefix)
!  call cla_get("--nc",inetcdf)
!  call cla_get("--file", cfile)
!end if

!if (len_trim(cfile) > 1) then
!  nfile = 1
!  files(1) = cfile
!else
!  print*, "No base files specified."
!end if

to_be_combined=varlist

! validate and prioritize
if (ibgn > iend) then
  stop ("Beginning index after end. Check input file or -b and -e options")
end if
print '("Begin: ",i4,", End: ",i4,", dry flag: ", i4)',ibgn,iend,iwetdry 

if (len_trim(to_be_combined)==0) then
  print *,'combine all variables'
else
  print *,'combine variables: ',trim(to_be_combined)
end if
!print*,"File:"
!do ifile = 1,nfile
!  print*, files(ifile)
!end do

end subroutine


!===============================================================================
program combine_output
implicit none
integer, parameter                   :: nfilemax = 20
integer                              :: nfile  !< number of file base names
character(len=30),dimension(nfilemax) :: files  !< base names (e.g. elev.61) of files to combine
integer                              :: ibgn   !< first output file index to process
integer                              :: iend   !< last output file index to process
integer                              :: iwetdry !dry option
integer                              :: inetcdf !< netcdf flag
character(len=1024)                  :: to_be_combined ! list of variables to be combined
character(len=1024)                  :: output_prefix ! output file
call combine_output11_input(ibgn,iend,iwetdry,to_be_combined,output_prefix)
call combine_output11(ibgn,iend,iwetdry,to_be_combined,output_prefix)
end program



