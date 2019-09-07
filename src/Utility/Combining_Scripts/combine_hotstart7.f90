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
! Read in rank-specific hotstart outputs and combine them into hotstart.in.
! Gobal-local mappings are read in from separate files.

! Usage: ./combine_hotstart6 -h for help, inside outputs/
! Inputs:
!        rank-specific hotstart files; 
!        local_to_global_*; 
!        screen: ntracers; it_char (iteration #)
! Output: hotstart_it=[time step].nc 
!
!  ifort -O2 -cpp -CB -mcmodel=medium -assume byterecl -g -traceback -o combine_hotstart7.exe ../UtilLib/argparse.f90 combine_hotstart7.f90 -I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff

! Revisions: v6 with nc
!================================================================================
program combine_hotstart
integer :: istep
integer :: nproc
integer :: ntracer
call combine_hotstart7_input(istep)
call combine_hotstart7(istep)
end program combine_hotstart

subroutine combine_hotstart7_input(istep)
use argparse
implicit none
integer       :: istep   !< global iteration/step 
!integer       :: nproc   !< number of processors in original simulation
!integer       :: ntracer !< number of tracers 

cmd_name = "combine_hotstart7"
call cla_init(cmd_name)
call cla_register('-i','--iteration', 'global iteration to combine (before _00*_hotstart', cla_int,'')
!call cla_register('-p','--nproc', 'number of procs to combine', cla_int,'1')
!call cla_register('-t','--ntracer','number of tracers (including T,S)', cla_int  ,'2')
call cla_validate
        
call cla_get("--iteration",istep)
!call cla_get("--nproc",nproc)
!call cla_get("--ntracer",ntracer)
end subroutine

subroutine combine_hotstart7(istep)
!-------------------------------------------------------------------------------
  use netcdf
  implicit real(8)(a-h,o-z),integer(i-n)
  parameter(nbyte=4)
  character(len=12) :: it_char
  character(len=72) :: fgb,fgb2,fdb,variable_nm
  character(len=72),allocatable :: dim_name(:)
  integer :: lfgb,lfdb       ! Length of processor specific global output file name
  integer :: one_dim,xtype,var1d_dims(1),var2d_dims(2),var3d_dims(3), &
 &data_start_1d(1),data_start_2d(2),data_start_3d(3), &
 &data_count_1d(1),data_count_2d(2),data_count_3d(3)
  integer,allocatable :: dimids(:),idims(:),iwork1(:),iwork2(:), &
 &ner(:),npr(:),nsr(:), & !resident (no ghosts)
 &ielg(:,:),iplg(:,:),islg(:,:)
  real(kind=8),allocatable :: work1(:,:,:),work2(:,:,:)
!-------------------------------------------------------------------------------
      
!-------------------------------------------------------------------------------
! Aquire user inputs
!-------------------------------------------------------------------------------

  !Debug
  print*, 'istep=',istep

! Read mapping info
  open(10,file='local_to_global_0000',status='old')
  read(10,*)ns_global,ne_global,np_global,nvrt,nproc,ntracers
  close(10)
  !For dim purpose
  if(np_global>ns_global.or.ne_global>ns_global) stop 'ns_global not max'
  allocate(ner(0:nproc-1),npr(0:nproc-1),nsr(0:nproc-1))

  fdb='local_to_global_0000'; fdb=adjustl(fdb)
  lfdb=len_trim(fdb)
  mxner=0; mxnpr=0; mxnsr=0
  do irank=0,nproc-1
    write(fdb(lfdb-3:lfdb),'(i4.4)') irank
    open(10,file=fdb,status='old')
    read(10,*)!ns_global,ne_global,np_global,nvrt,nproc,ntracers

    read(10,*) 
    read(10,*)ner(irank); do i=1,ner(irank); read(10,*); enddo;
    read(10,*)npr(irank); do i=1,npr(irank); read(10,*); enddo;
    read(10,*)nsr(irank); do i=1,nsr(irank); read(10,*); enddo;
    close(10)
  enddo !irank
  mxner=maxval(ner)
  mxnpr=maxval(npr)
  mxnsr=maxval(nsr)

  allocate(ielg(0:nproc-1,mxner),iplg(0:nproc-1,mxnpr),islg(0:nproc-1,mxnsr), &
     &stat=istat)
  if(istat/=0) stop 'Allocation error (3)'

  do irank=0,nproc-1
    write(fdb(lfdb-3:lfdb),'(i4.4)') irank
    open(10,file=fdb,status='old')
    read(10,*); read(10,*)
    read(10,*)ner(irank)
    do i=1,ner(irank)
      read(10,*)j,ielg(irank,i)
    enddo
    read(10,*)npr(irank)
    do i=1,npr(irank)
      read(10,*)j,iplg(irank,i)
    enddo
    read(10,*)nsr(irank)
    do i=1,nsr(irank)
      read(10,*)j,islg(irank,i)
    enddo
    close(10)
  enddo !irank

  print*, 'Global quantities:',ne_global,np_global,ns_global

!-------------------------------------------------------------------------------
! Read/write hotstart files
! See schism_step.F90 for rules
!-------------------------------------------------------------------------------
  ! Open file
  write(it_char,'(i12)') istep
  it_char=adjustl(it_char)  !place blanks at end
  it_len=len_trim(it_char)  !length without trailing blanks
  fgb='0000_'//it_char(1:it_len); 
  fgb=adjustl(fgb); lfgb=len_trim(fgb);
  !print*, 'suffix is:',fgb(1:lfgb)

  !Query # of vars
  iret=nf90_open('hotstart_0000_'//it_char(1:it_len)//'.nc',OR(NF90_NETCDF4,NF90_NOWRITE),ncid2)
  iret=nf90_inquire(ncid2,nDimensions=ndimensions,nVariables=nvars)
  write(*,*)'nvars=',nvars,ndimensions
  allocate(idims(ndimensions),dim_name(ndimensions),dimids(ndimensions))
  do i=1,ndimensions
    iret=nf90_inquire_dimension(ncid2,i,name=dim_name(i),len=idims(i))
    !write(99,*)'dim=',i,idims(i),dim_name(i)
  enddo !i
  !Get type II arrays (including from modules)
  iret=nf90_inq_varid(ncid2, "time",i);
  iret=nf90_get_var(ncid2,i,time);
  iret=nf90_inq_varid(ncid2, "it",i);
  iret=nf90_get_var(ncid2,i,iths);
  iret=nf90_inq_varid(ncid2, "ifile",i);
  iret=nf90_get_var(ncid2,i,ifile);
  !From modules
  ice_free_flag=-99 !flag
  iret=nf90_inq_varid(ncid2, "ice_free_flag",i);
  if(iret==NF90_NOERR) then
    iret=nf90_get_var(ncid2,i,ice_free_flag);
    if(iret/=NF90_NOERR) stop 'cannot get ice_free_flag'
    print*, 'Found ice_free_flag:',ice_free_flag
  endif
  print*, 'static info:',time,iths,ifile
  iret=nf90_close(ncid2)

  !Open output file
  iret=nf90_create('hotstart_it='//it_char(1:it_len)//'.nc',OR(NF90_NETCDF4,NF90_CLOBBER),ncid)
  iret=nf90_def_dim(ncid,'node',np_global,node_dim)
  iret=nf90_def_dim(ncid,'elem',ne_global,nele_dim)
  iret=nf90_def_dim(ncid,'side',ns_global,nedge_dim)
  do i=4,ndimensions !skip local node/elem/side #s
    iret=nf90_def_dim(ncid,trim(adjustl(dim_name(i))),idims(i),itmp)
    !so that we can easily refer to the dim id in writes
    if(itmp/=i) stop 'pre- and post-comb conflict' 
  enddo !i
  !Extra dims
  iret=nf90_def_dim(ncid,'one_new',1,one_dim)
  iret=nf90_enddef(ncid)
  
  loop1: do m=1,nvars
    do irank=0,nproc-1
      fgb2=fgb
      fgb2=adjustl(fgb2)
      write(fgb2(1:4),'(i4.4)') irank
      iret=nf90_open('hotstart_'//fgb2(1:lfgb)//'.nc',OR(NF90_NETCDF4,NF90_NOWRITE),ncid2)
      if(iret.ne.NF90_NOERR) then
        print*, nf90_strerror(iret); stop
      endif
!      write(99,*)'file=',m,irank,'hotstart_'//fgb2(1:lfgb)//'.nc',ncid2

      iret=nf90_inquire_variable(ncid2,m,name=variable_nm,xtype=xtype,ndims=ndims)
      if(ndims<=0.or.ndims>ndimensions) stop '0 dim encountered'
      variable_nm=adjustl(variable_nm); len_var=len_trim(variable_nm)
!      write(99,*)'var dim:',m,irank,variable_nm(1:len_var),ndims

      if(ndims>ndimensions) stop 'dim overflow'
      iret=nf90_inquire_variable(ncid2,m,dimids=dimids(1:ndims))
      if(iret.ne.NF90_NOERR) then
        print*, nf90_strerror(iret); stop
      endif
      !Re-read dims (rank specific)
      do i=1,ndimensions
        iret=nf90_inquire_dimension(ncid2,i,len=idims(i))
      enddo !i
!      write(99,*)'var dim(2):',m,irank,dimids(1:ndims),(idims(dimids(i)),i=1,ndims)

      if(dimids(ndims)>3) then
        iret=nf90_close(ncid2)
        cycle loop1 !skip static arrays
      endif

      !Define output dims
      if(dimids(ndims)==1) then !node array
        npse=np_global
        npse_dim=node_dim
        npse_lcl=npr(irank)
      else if(dimids(ndims)==2) then !elem array
        npse=ne_global
        npse_dim=nele_dim
        npse_lcl=ner(irank)
      else if(dimids(ndims)==3) then !side array
        npse=ns_global
        npse_dim=nedge_dim
        npse_lcl=nsr(irank)
      else
        stop 'last dim is not node/elem/side'
      endif !dimids

      if(idims(dimids(ndims))/=npse_lcl) then
        print*, 'last dim wrong:',m,irank,npse_lcl,idims(dimids(ndims))
        stop
      endif

      !alloc
      if(ndims==1) then !1D
        itmp1=idims(dimids(1))
        allocate(iwork1(npse_lcl),work1(1,1,itmp1),stat=istat)
        if(irank==0) allocate(iwork2(npse),work2(1,1,npse),stat=istat)
      else if(ndims==2) then !2D array
        itmp1=idims(dimids(1)); itmp2=idims(dimids(2))
        allocate(iwork1(npse_lcl),work1(1,itmp1,itmp2),stat=istat)
        if(irank==0) allocate(iwork2(npse),work2(1,itmp1,npse),stat=istat)
      else if(ndims==3) then
        itmp1=idims(dimids(1))
        itmp2=idims(dimids(2)); itmp3=idims(dimids(3))
        allocate(iwork1(npse_lcl),work1(itmp1,itmp2,itmp3),stat=istat)
        if(irank==0) allocate(iwork2(npse),work2(itmp1,itmp2,npse),stat=istat)
      else
       stop '>3D array'
      endif !ndims
      if(istat/=0) stop 'Allocation error (7)'

      !Check if the var is already defined
      iret2=nf90_inq_varid(ncid,variable_nm(1:len_var),i)

      if(xtype==NF90_INT) then !1D arrays only
        data_start_1d(1)=1; data_count_1d(1)=npse_lcl !npr(irank)
        iret=nf90_get_var(ncid2,m,iwork1(1:npse_lcl),data_start_1d,data_count_1d)

        !Define output var
        if(iret2/=NF90_NOERR) then !not defined yet
          iret=nf90_redef(ncid)
          var1d_dims(1)=npse_dim
          iret=nf90_def_var(ncid,variable_nm(1:len_var),NF90_INT,var1d_dims,ivarid)
          iret=nf90_enddef(ncid)
        endif

        !write(99,*)'Im in:',m,iwork1,ivarid
      else !doubles
        if(ndims==1) then
          data_start_1d(1)=1; data_count_1d(1)=npse_lcl !npr(irank)
          iret=nf90_get_var(ncid2,m,work1(1,1,1:npse_lcl),data_start_1d,data_count_1d)

          !Define output var
          if(iret2/=NF90_NOERR) then !not defined yet
            iret=nf90_redef(ncid)
            var1d_dims(1)=npse_dim
            iret=nf90_def_var(ncid,variable_nm(1:len_var),NF90_DOUBLE,var1d_dims,ivarid)
            iret=nf90_enddef(ncid)
          endif
!          write(99,*)'Im in:',m,irank,ivarid,work1(1,1,1:npse_lcl)

        else if(ndims==2) then !2D
          data_start_2d(1:2)=1
          data_count_2d(1)=idims(dimids(1)); data_count_2d(2)=npse_lcl !npr(irank)
          iret=nf90_get_var(ncid2,m,work1(1,1:data_count_2d(1),1:npse_lcl),data_start_2d,data_count_2d)
 
          !Define output var
          if(iret2/=NF90_NOERR) then !not defined yet
            iret=nf90_redef(ncid)
            !Pre- and post-comb nc share same dims of ID>3
            var2d_dims(1)=dimids(1); var2d_dims(2)=npse_dim
            iret=nf90_def_var(ncid,variable_nm(1:len_var),NF90_DOUBLE,var2d_dims,ivarid)
            iret=nf90_enddef(ncid)
          endif

        else if(ndims==3) then !3D array
          data_start_3d(1:3)=1
          data_count_3d(1:2)=idims(dimids(1:2)); data_count_3d(3)=npse_lcl !npr(irank)
          iret=nf90_get_var(ncid2,m,work1(1:data_count_3d(1),1:data_count_3d(2),1:npse_lcl), &
     &data_start_3d,data_count_3d)
 
          !Define output var
          if(iret2/=NF90_NOERR) then !not defined yet
            iret=nf90_redef(ncid)
            var3d_dims(1:2)=dimids(1:2); var3d_dims(3)=npse_dim
            iret=nf90_def_var(ncid,variable_nm(1:len_var),NF90_DOUBLE,var3d_dims,ivarid)
            iret=nf90_enddef(ncid)
          endif
        else
          stop 'unknown array(2)'
        endif !ndims
      endif !xtype

      !Put into global array
      do i=1,npse_lcl !npr(irank)
        if(dimids(ndims)==1) then !node
          ip=iplg(irank,i)
        else if(dimids(ndims)==2) then
          ip=ielg(irank,i)
        else !side
          ip=islg(irank,i)
        endif
        iwork2(ip)=iwork1(i)
        work2(:,:,ip)=work1(:,:,i)
      enddo !i

      iret=nf90_close(ncid2)
      deallocate(iwork1,work1)
    enddo !irank

    !write
    if(xtype==NF90_INT) then !1D arrays only
      data_start_1d(1)=1; data_count_1d(1)=npse
      iret=nf90_put_var(ncid,ivarid,iwork2,data_start_1d,data_count_1d)
    else !doubles
      if(ndims==1) then
        data_start_1d(1)=1; data_count_1d(1)=npse
        iret=nf90_put_var(ncid,ivarid,work2(1,1,1:npse),data_start_1d,data_count_1d)
      else if(ndims==2) then
        data_start_2d(1:2)=1
        !dimids(1) same across all ranks
        data_count_2d(1)=idims(dimids(1)); data_count_2d(2)=npse
        iret=nf90_put_var(ncid,ivarid,work2(1,1:data_count_2d(1),1:npse),data_start_2d,data_count_2d)
      else if(ndims==3) then
        data_start_3d(1:3)=1
        data_count_3d(1:2)=idims(dimids(1:2)); data_count_3d(3)=npse
        iret=nf90_put_var(ncid,ivarid,work2(1:data_count_3d(1),1:data_count_3d(2),1:npse), &
     &data_start_3d,data_count_3d)
      endif !ndims
    endif !xtype

    !Debug
!    if(variable_nm(1:len_var).eq.'eta2') then
!      write(99,*)'Check output:',m,variable_nm(1:len_var)
!      do i=1,npse
!        write(99,*)i,work2(1,1,i)
!      enddo !i

      !iret=nf90_close(ncid)
      !stop
!    endif

    deallocate(iwork2,work2)
  enddo loop1 !m: vars

  !Append type II arrays
  iret=nf90_redef(ncid)
  var1d_dims(1)=one_dim
  iret=nf90_def_var(ncid,'time',NF90_DOUBLE,var1d_dims,i1)
  iret=nf90_def_var(ncid,'iths',NF90_INT,var1d_dims,i2)
  iret=nf90_def_var(ncid,'ifile',NF90_INT,var1d_dims,i3)
  if(ice_free_flag>=0) iret=nf90_def_var(ncid,'ice_free_flag',NF90_INT,var1d_dims,i4)
  iret=nf90_enddef(ncid)
  iret=nf90_put_var(ncid,i1,time)
  iret=nf90_put_var(ncid,i2,iths)
  iret=nf90_put_var(ncid,i3,ifile)
  if(ice_free_flag>=0) iret=nf90_put_var(ncid,i4,ice_free_flag)

  iret=nf90_close(ncid)

end subroutine combine_hotstart7

