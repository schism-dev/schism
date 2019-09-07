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

!===============================================================================
! Convert hotstart.in to hotstart.nc (no module data)

! Inputs:
!        hotstart.in (unformatted binary); screen inputs.
! Output: hotstart.nc 
!
!  ifort -O2 -CB -mcmodel=medium -assume byterecl -o convert_hotstart_nc.exe convert_hotstart_nc.f90 -I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff

!===============================================================================

program combine_hotstart1
  use netcdf
  implicit real(8)(a-h,o-z),integer(i-n)
  parameter(nbyte=4)
!  character(12) :: it_char
!  character(72) :: fgb,fgb2,fdb  ! Processor specific global output file name
!  integer :: lfgb,lfdb       ! Length of processor specific global output file name
!  allocatable ner(:),npr(:),nsr(:)
!  allocatable ielg(:,:),iplg(:,:),islg(:,:)
  allocatable idry_e(:),we(:,:),tsel(:,:,:),idry_s(:),su2(:,:),sv2(:,:)
  allocatable tsd(:,:),ssd(:,:),idry(:),eta2(:),tnd(:,:),snd(:,:)
  allocatable tem0(:,:),sal0(:,:),q2(:,:),xl(:,:),dfv(:,:),dfh(:,:)
  allocatable dfq1(:,:),dfq2(:,:),qnon(:,:),trnd(:,:,:),trel(:,:,:),trnd0(:,:,:)
  allocatable intv(:),zrat(:),swild(:,:)

  integer :: elem_dim,side_dim,one_dim,three_dim,nwild(10000),var1d_dim(1), &
 &var2d_dim(2),var3d_dim(3)
      
  print*, 'Input ne, np, ns, nvrt, and ntracers (including modules):'
  read*, ne,np,ns,nvrt,ntracers
  if(ntracers<2) stop 'wrong ntracers'

  allocate(idry_e(ne),we(nvrt,ne), &
           idry_s(ns),su2(nvrt,ns),sv2(nvrt,ns), &
           tsd(nvrt,ns),ssd(nvrt,ns), &
           idry(np),eta2(np),trnd(ntracers,nvrt,np), &
           trnd0(ntracers,nvrt,np),trel(ntracers,nvrt,ne),q2(nvrt,np), &
           xl(nvrt,np),dfv(nvrt,np),dfh(nvrt,np), &
           dfq1(nvrt,np),dfq2(nvrt,np),qnon(nvrt,np), &
           swild(nvrt1,7+2*ntracers),intv(nvrt1),zrat(nvrt1),stat=istat)
  if(istat/=0) stop 'Allocation error (2)'

!-------------------------------------------------------------------------------
! Read hotstart files
!-------------------------------------------------------------------------------
  open(36,file='hotstart.in',form='unformatted',status='old')
  read(36) time,it,ifile
  do i=1,ne
    read(36) j,idry_e(i),(we(j,i),(trel(l,j,i),l=1,ntracers),j=1,nvrt)
  enddo !i
  print*, 'done reading elem...'
  do i=1,ns
    read(36) j,idry_s(i),(su2(j,i),sv2(j,i),tsd(j,i),ssd(j,i),j=1,nvrt)
  enddo !i
  print*, 'done reading side...'
  do i=1,np
    read(36) j,eta2(i),idry(i),(trnd(:,j,i),trnd0(:,j,i),q2(j,i),xl(j,i), &
             dfv(j,i),dfh(j,i),dfq1(j,i),dfq2(j,i),qnon(j,i),j=1,nvrt)
  enddo !i
  print*, 'done reading node...'

! Other modules (SED etc)

  close(36)

  !write
  j=nf90_create('hotstart.nc',OR(NF90_CLOBBER,NF90_NETCDF4),ncid_hot)
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

  !Note the order of multi-dim arrays not reversed here!
  !As long as the write is consistent with def it's fine 
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

  !var1d_dim(1)=side_dim
  !j=nf90_def_var(ncid_hot,'xcj',NF90_DOUBLE,var1d_dim,nwild(20))
  !j=nf90_def_var(ncid_hot,'ycj',NF90_DOUBLE,var1d_dim,nwild(21))
  !var1d_dim(1)=node_dim
  !j=nf90_def_var(ncid_hot,'xnd',NF90_DOUBLE,var1d_dim,nwild(22))
  !j=nf90_def_var(ncid_hot,'ynd',NF90_DOUBLE,var1d_dim,nwild(23))
  !var2d_dim(1)=nvrt_dim; var2d_dim(2)=node_dim
  !j=nf90_def_var(ncid_hot,'uu2',NF90_DOUBLE,var2d_dim,nwild(24))
  !j=nf90_def_var(ncid_hot,'vv2',NF90_DOUBLE,var2d_dim,nwild(25))

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
  j=nf90_put_var(ncid_hot,nwild(9),trel(:,:,1:ne),(/1,1,1/),(/ntracers,nvrt,ne/))
  j=nf90_put_var(ncid_hot,nwild(10),su2(:,1:ns),(/1,1/),(/nvrt,ns/))
  j=nf90_put_var(ncid_hot,nwild(11),sv2(:,1:ns),(/1,1/),(/nvrt,ns/))
  j=nf90_put_var(ncid_hot,nwild(12),trnd(:,:,1:np),(/1,1,1/),(/ntracers,nvrt,np/))
  j=nf90_put_var(ncid_hot,nwild(13),trnd0(:,:,1:np),(/1,1,1/),(/ntracers,nvrt,np/))
  j=nf90_put_var(ncid_hot,nwild(14),q2(:,1:np),(/1,1/),(/nvrt,np/))
  j=nf90_put_var(ncid_hot,nwild(15),xl(:,1:np),(/1,1/),(/nvrt,np/))
  j=nf90_put_var(ncid_hot,nwild(16),dfv(:,1:np),(/1,1/),(/nvrt,np/))
  j=nf90_put_var(ncid_hot,nwild(17),dfh(:,1:np),(/1,1/),(/nvrt,np/))
  j=nf90_put_var(ncid_hot,nwild(18),dfq1(:,1:np),(/1,1/),(/nvrt,np/))
  j=nf90_put_var(ncid_hot,nwild(19),dfq2(:,1:np),(/1,1/),(/nvrt,np/))
  
  j=nf90_close(ncid_hot)
end program combine_hotstart1
