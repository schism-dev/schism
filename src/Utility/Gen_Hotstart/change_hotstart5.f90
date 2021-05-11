! Modify hotstart.nc; use this as a template

! Inputs:
!        hotstart.nc, hgrid.gr3, vgrid.in; some const below
! Output: hotstart.nc.new 
!
!  ifort -O2 -CB -g -traceback -mcmodel=medium -assume byterecl -o change_hotstart5 change_hotstart5.f90 ../UtilLib/schism_geometry.f90 ../UtilLib/compute_zcor.f90 -I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff

!===============================================================================

program combine_hotstart1
  use netcdf
  use schism_geometry_mod
  use compute_zcor

  implicit real(8)(a-h,o-z),integer(i-n)
!  character(12) :: it_char
  allocatable idry_e(:),we(:,:),tsel(:,:,:),idry_s(:),su2(:,:),sv2(:,:)
  allocatable idry(:),eta2(:),tnd(:,:),snd(:,:)
  allocatable tem0(:,:),sal0(:,:),q2(:,:),xl(:,:),dfv(:,:),dfh(:,:)
  allocatable dfq1(:,:),dfq2(:,:),tr_nd(:,:,:),tr_el(:,:,:),tr_nd0(:,:,:)
  allocatable intv(:),zrat(:),swild(:,:)
  real*8, allocatable :: xnd(:),ynd(:),dp(:),xcj(:),ycj(:),ztot(:),sigma(:),znl(:,:),sigma_lcl(:,:)
  integer, allocatable :: elnode(:,:),i34(:),ic3(:,:),elside(:,:),isdel(:,:),isidenode(:,:), &
 &kbp(:)
  real*8 :: time(1)
  integer :: iths(1),ifile(1)

  integer :: elem_dim,side_dim,one_dim,three_dim,nwild(10000),var1d_dim(1), &
 &var2d_dim(2),var3d_dim(3)
      
! Consts
  ntracers=2

! hgrid
  open(14,file='hgrid.gr3',status='old')
  read(14,*); read(14,*)ne,np
  allocate(xnd(np),ynd(np),dp(np),i34(ne),elnode(4,ne),stat=istat)
  do i=1,np
    read(14,*)j,xnd(i),ynd(i),dp(i)
  enddo !i
  xmin=minval(xnd); ymin=minval(ynd)
  xmax=maxval(xnd); ymax=maxval(ynd)
  do i=1,ne
    read(14,*)j,i34(i),elnode(1:i34(i),i)
  enddo !i
  close(14)

! Gemoetry
  call compute_nside(np,ne,i34,elnode(1:4,1:ne),ns)
  print*, '# of sides=',ns
  allocate(ic3(4,ne),elside(4,ne),isdel(2,ns),isidenode(2,ns),xcj(ns),ycj(ns),stat=istat)
  if(istat/=0) stop 'Allocation error: side(0)'
  call schism_geometry_double(np,ne,ns,xnd,ynd,i34,elnode(1:4,1:ne),ic3(1:4,1:ne), &
  &elside(1:4,1:ne),isdel,isidenode,xcj,ycj)

! Vgrid (SZ)
  open(19,file='vgrid.in',status='old')
  read(19,*); read(19,*)nvrt 
  close(19)

  allocate(ztot(nvrt),sigma(nvrt),sigma_lcl(nvrt,np),znl(nvrt,np),kbp(np))
  call get_vgrid_double('vgrid.in',np,nvrt,ivcor,kz,h_s,h_c,theta_b, &
 &theta_f,ztot,sigma,sigma_lcl,kbp)
!  allocate(ztmp(nvrt),ztmp2(nvrt,3),stat=istat)

!  open(19,file='vgrid.in',status='old')
!  read(19,*)ivcor
!  read(19,*) nvrt,kz,h_s !kz>=1
!  if(ivcor/=2) stop 'ivcor/=2'

!  if(nvrt<2) stop 'nvrt<2'
!  if(kz<1) then !.or.kz>nvrt-2) then
!    write(*,*)'Wrong kz:',kz
!    stop
!  endif
!  if(h_s<6) then
!    write(*,*)'h_s needs to be larger:',h_s
!    stop
!  endif

!  ! # of z-levels excluding "bottom" at h_s
!  read(19,*) !for adding comment "Z levels"
!  do k=1,kz-1
!    read(19,*)j,ztot(k)
!    if(ztot(k)>=-h_s) then
!      print*, 'Illegal Z level:',k
!      stop
!    endif
!    if(k>1) then; if(ztot(k)<=ztot(k-1)) then
!      print*, 'z-level inverted:',k
!      stop
!    endif; endif
!  enddo !k
!  read(19,*) !level kz       
!  ! In case kz=1, there is only 1 ztot(1)=-h_s
!  ztot(kz)=-h_s

!  nsig=nvrt-kz+1 !# of S levels (including "bottom" & f.s.)
!  read(19,*) !for adding comment "S levels"
!  read(19,*)h_c,theta_b,theta_f
!  if(h_c<5.or.h_c>=h_s) then !large h_c to avoid 2nd type abnormality
!    print*, 'h_c needs to be larger avoid 2nd type abnormality; &
!  &do u want to continue? '
!    stop
!  endif
!  if(theta_b<0.or.theta_b>1) then
!    write(*,*)'Wrong theta_b:',theta_b
!    stop
!  endif
!  if(theta_f<=0) then
!    write(*,*)'Wrong theta_f:',theta_f
!    stop
!  endif

!  sigma(1)=-1 !bottom
!  sigma(nsig)=0 !surface
!  read(19,*) !level kz
!  do k=kz+1,nvrt-1
!    kin=k-kz+1
!    read(19,*) j,sigma(kin)
!    if(sigma(kin)<=sigma(kin-1).or.sigma(kin)>=0) then
!      write(*,*)'Check sigma levels at:',k,sigma(kin),sigma(kin-1)
!      stop
!    endif

!    write(98,*)'sigma=',k,sigma(kin)
!  enddo !k
!  read(19,*) !level nvrt
!  close(19)

  allocate(idry_e(ne),we(nvrt,ne), &
           idry_s(ns),su2(nvrt,ns),sv2(nvrt,ns), &
           idry(np),eta2(np),tr_nd(ntracers,nvrt,np), &
           tr_nd0(ntracers,nvrt,np),tr_el(ntracers,nvrt,ne),q2(nvrt,np), &
           xl(nvrt,np),dfv(nvrt,np),dfh(nvrt,np), &
           dfq1(nvrt,np),dfq2(nvrt,np), &
           swild(nvrt,7+2*ntracers),zrat(nvrt),stat=istat)
  if(istat/=0) stop 'Allocation error (2)'

! Compute zcor (assuming elev=0)
  do i=1,np
    if(ivcor==1) then
      !znl junks if dry
      do k=kbp(i)+1,nvrt-1
        znl(k,i)=dp(i)*sigma_lcl(k,i) !+eta2(nd,irec)
      enddo !k
      znl(kbp(i),i)=-dp(i) !to avoid underflow
      znl(nvrt,i)=0. !eta2(nd,irec) !to avoid underflow
    else
      call zcor_SZ_double(dp(i),0.d0,0.01d0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot,sigma,znl(:,i),idry(i),kbp(i))
    endif
  enddo !i

! Read in old hotstart
  iret=nf90_open('hotstart.nc',OR(NF90_NETCDF4,NF90_NOWRITE),ncid2)
  iret=nf90_inq_varid(ncid2, "time",i);
  iret=nf90_get_var(ncid2,i,time);
  iret=nf90_inq_varid(ncid2, "iths",i);
  iret=nf90_get_var(ncid2,i,iths);
  iret=nf90_inq_varid(ncid2, "ifile",i);
  iret=nf90_get_var(ncid2,i,ifile);
  print*, 'static info:',time(1),iths(1),ifile(1)

  iret=nf90_inq_varid(ncid2, "idry",i);
  iret=nf90_get_var(ncid2,i,idry);
  iret=nf90_inq_varid(ncid2, "idry_e",i);
  iret=nf90_get_var(ncid2,i,idry_e);
  iret=nf90_inq_varid(ncid2, "idry_s",i);
  iret=nf90_get_var(ncid2,i,idry_s);
  iret=nf90_inq_varid(ncid2, "eta2",i);
  iret=nf90_get_var(ncid2,i,eta2);
  iret=nf90_inq_varid(ncid2, "we",i);
  iret=nf90_get_var(ncid2,i,we);
  iret=nf90_inq_varid(ncid2, "tr_el",i);
  iret=nf90_get_var(ncid2,i,tr_el);
  iret=nf90_inq_varid(ncid2, "su2",i);
  iret=nf90_get_var(ncid2,i,su2);
  iret=nf90_inq_varid(ncid2, "sv2",i);
  iret=nf90_get_var(ncid2,i,sv2);
  iret=nf90_inq_varid(ncid2, "tr_nd",i);
  iret=nf90_get_var(ncid2,i,tr_nd);
  iret=nf90_inq_varid(ncid2, "tr_nd0",i);
  iret=nf90_get_var(ncid2,i,tr_nd0);
  iret=nf90_inq_varid(ncid2, "q2",i);
  iret=nf90_get_var(ncid2,i,q2);
  iret=nf90_inq_varid(ncid2, "xl",i);
  iret=nf90_get_var(ncid2,i,xl);
  iret=nf90_inq_varid(ncid2, "dfv",i);
  iret=nf90_get_var(ncid2,i,dfv);
  iret=nf90_inq_varid(ncid2, "dfh",i);
  iret=nf90_get_var(ncid2,i,dfh);
  iret=nf90_inq_varid(ncid2, "dfq1",i);
  iret=nf90_get_var(ncid2,i,dfq1);
  iret=nf90_inq_varid(ncid2, "dfq2",i);
  iret=nf90_get_var(ncid2,i,dfq2);
  iret=nf90_close(ncid2)
  print*, 'done reading hotstart'

  !Debug
!  write(99,*)'T:'
!  write(99,*)np
!  do i=1,np
!    write(99,*)i,real(xnd(i)),real(ynd(i)),real(tr_nd(1,nvrt,i))
!  enddo !i

! Modify i.c. here
  do i=1,np
    do k=1,nvrt
      !tr_nd(1,k,i)=temp_0-delta_T*tmp+hat_T
      tr_nd(1,k,i)=tr_nd(1,k,i)+2
    enddo !k

!    write(96,'(i12,3(1x,e20.12))')i,xnd(i),ynd(i),tr_nd(1,nvrt,i)
!    write(97,'(i12,3(1x,e20.12))')i,xnd(i),ynd(i),tr_nd(1,1,i)
  enddo !i=1,np
  tr_nd0=tr_nd

  do i=1,ne
    do k=1,nvrt
!      tr_el(1,k,i)=sum(tr_nd(1,k,elnode(1:i34(i),i))+tr_nd(1,k-1,elnode(1:i34(i),i)))/i34(i)*0.5d0
      tr_el(1,k,i)=tr_el(1,k,i)+2
    enddo !k
!    tr_el(1,1,i)=tr_el(1,2,i)
  enddo !i

  !write
  j=nf90_create('hotstart.nc.new',OR(NF90_CLOBBER,NF90_NETCDF4),ncid_hot)
  j=nf90_def_dim(ncid_hot,'node',np,node_dim)
  j=nf90_def_dim(ncid_hot,'elem',ne,elem_dim)
  j=nf90_def_dim(ncid_hot,'side',ns,side_dim)
  j=nf90_def_dim(ncid_hot,'nVert',nvrt,nvrt_dim)
  j=nf90_def_dim(ncid_hot,'ntracers',ntracers,ntracers_dim)
  j=nf90_def_dim(ncid_hot,'one',1,one_dim)
  j=nf90_def_dim(ncid_hot,'three',3,three_dim)
  var1d_dim(1)=one_dim
  j=nf90_def_var(ncid_hot,'time',NF90_DOUBLE,var1d_dim,nwild(1))
  j=nf90_def_var(ncid_hot,'iths',NF90_INT,var1d_dim,nwild(2))
  j=nf90_def_var(ncid_hot,'ifile',NF90_INT,var1d_dim,nwild(3))

  var1d_dim(1)=elem_dim
  j=nf90_def_var(ncid_hot,'idry_e',NF90_INT,var1d_dim,nwild(4))
  var1d_dim(1)=side_dim
  j=nf90_def_var(ncid_hot,'idry_s',NF90_INT,var1d_dim,nwild(5))
  var1d_dim(1)=node_dim
  j=nf90_def_var(ncid_hot,'idry',NF90_INT,var1d_dim,nwild(6))
  j=nf90_def_var(ncid_hot,'eta2',NF90_DOUBLE,var1d_dim,nwild(7))

  var2d_dim(1)=nvrt_dim; var2d_dim(2)=elem_dim
  j=nf90_def_var(ncid_hot,'we',NF90_DOUBLE,var2d_dim,nwild(8))
  var3d_dim(1)=ntracers_dim; var3d_dim(2)=nvrt_dim; var3d_dim(3)=elem_dim
  j=nf90_def_var(ncid_hot,'tr_el',NF90_DOUBLE,var3d_dim,nwild(9))
  var2d_dim(1)=nvrt_dim; var2d_dim(2)=side_dim
  j=nf90_def_var(ncid_hot,'su2',NF90_DOUBLE,var2d_dim,nwild(10))
  j=nf90_def_var(ncid_hot,'sv2',NF90_DOUBLE,var2d_dim,nwild(11))
  var3d_dim(1)=ntracers_dim; var3d_dim(2)=nvrt_dim; var3d_dim(3)=node_dim
  j=nf90_def_var(ncid_hot,'tr_nd',NF90_DOUBLE,var3d_dim,nwild(12))
  j=nf90_def_var(ncid_hot,'tr_nd0',NF90_DOUBLE,var3d_dim,nwild(13))
  var2d_dim(1)=nvrt_dim; var2d_dim(2)=node_dim
  j=nf90_def_var(ncid_hot,'q2',NF90_DOUBLE,var2d_dim,nwild(14))
  j=nf90_def_var(ncid_hot,'xl',NF90_DOUBLE,var2d_dim,nwild(15))
  j=nf90_def_var(ncid_hot,'dfv',NF90_DOUBLE,var2d_dim,nwild(16))
  j=nf90_def_var(ncid_hot,'dfh',NF90_DOUBLE,var2d_dim,nwild(17))
  j=nf90_def_var(ncid_hot,'dfq1',NF90_DOUBLE,var2d_dim,nwild(18))
  j=nf90_def_var(ncid_hot,'dfq2',NF90_DOUBLE,var2d_dim,nwild(19))
  j=nf90_enddef(ncid_hot)

  !Write
  j=nf90_put_var(ncid_hot,nwild(1),time) 
  j=nf90_put_var(ncid_hot,nwild(2),it) 
  j=nf90_put_var(ncid_hot,nwild(3),ifile) 
  j=nf90_put_var(ncid_hot,nwild(4),idry_e,(/1/),(/ne/))
  j=nf90_put_var(ncid_hot,nwild(5),idry_s,(/1/),(/ns/))
  j=nf90_put_var(ncid_hot,nwild(6),idry,(/1/),(/np/))
  j=nf90_put_var(ncid_hot,nwild(7),eta2,(/1/),(/np/))
  j=nf90_put_var(ncid_hot,nwild(8),we(:,1:ne),(/1,1/),(/nvrt,ne/))
  j=nf90_put_var(ncid_hot,nwild(9),tr_el(:,:,1:ne),(/1,1,1/),(/ntracers,nvrt,ne/))
  j=nf90_put_var(ncid_hot,nwild(10),su2(:,1:ns),(/1,1/),(/nvrt,ns/))
  j=nf90_put_var(ncid_hot,nwild(11),sv2(:,1:ns),(/1,1/),(/nvrt,ns/))
  j=nf90_put_var(ncid_hot,nwild(12),tr_nd(:,:,1:np),(/1,1,1/),(/ntracers,nvrt,np/))
  j=nf90_put_var(ncid_hot,nwild(13),tr_nd0(:,:,1:np),(/1,1,1/),(/ntracers,nvrt,np/))
  j=nf90_put_var(ncid_hot,nwild(14),q2(:,1:np),(/1,1/),(/nvrt,np/))
  j=nf90_put_var(ncid_hot,nwild(15),xl(:,1:np),(/1,1/),(/nvrt,np/))
  j=nf90_put_var(ncid_hot,nwild(16),dfv(:,1:np),(/1,1/),(/nvrt,np/))
  j=nf90_put_var(ncid_hot,nwild(17),dfh(:,1:np),(/1,1/),(/nvrt,np/))
  j=nf90_put_var(ncid_hot,nwild(18),dfq1(:,1:np),(/1,1/),(/nvrt,np/))
  j=nf90_put_var(ncid_hot,nwild(19),dfq2(:,1:np),(/1,1/),(/nvrt,np/))
  
  j=nf90_close(ncid_hot)

  print*, 'Finished'
end program combine_hotstart1
