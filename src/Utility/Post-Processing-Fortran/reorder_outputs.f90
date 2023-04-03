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
!****************************************************************************************

!	Read nc outputs from scribe I/O versions and re-order according to
!	sub-domain for efficient subsetting.
!       Works for mixed tri/quad outputs on NODE based vars.
!       Inputs: screen; nc file, out2d*.nc, and local_to_global*
!       Outputs: *_ro.nc (re-ordered)
!****************************************************************************************
!     ifort -O2 -assume byterecl -o reorder_outputs.exe reorder_outputs.f90 -I$NETCDF/include  -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff

      program read_out
      use netcdf
!      use extract_mod2
!      use schism_geometry_mod
!      use compute_zcor

!      parameter(nbyte=4)
      character(len=30) :: file63,varname,varname2,outname(3),file62,file64,file_char,file61
      character(len=12) :: it_char
!      character*48 data_format
      logical::lexist,lexist2
      integer :: nx(4,4,3),dimids(100),idims(100)
      allocatable :: sigma(:),cs(:),ztot(:)
      allocatable :: outvar(:,:,:),out(:,:,:),icum(:,:,:),eta2(:),ztmp(:,:)
      allocatable :: xcj(:),ycj(:),dp(:),dps(:),kbs(:),xctr(:),yctr(:),dpe(:),kbe(:)
      allocatable :: zs(:,:),ze(:,:),idry(:),outs(:,:,:),oute(:,:,:),rstat2d(:,:),rviolator(:)
      allocatable :: ic3(:,:),isdel(:,:),isidenode(:,:),kbp(:),sigma_lcl(:,:),kbp00(:)
      allocatable :: dpout(:),etaout(:),out3d(:,:)
      integer, allocatable :: elside(:,:),idry_e(:),nne(:),indel(:,:),iviolator(:)
      integer,allocatable :: i34(:),elnode(:,:),np_lcl(:),ne_lcl(:),ns_lcl(:),nm(:,:), &
     &iplg(:,:),ielg(:,:),islg(:,:),iegl_rank(:),nm2(:,:)
      
      integer :: one_dim,time_dim,char_len,start_1d(1),start_2d(2),start_3d(3),start_4d(4), &
     &count_1d(1),count_2d(2),count_3d(3),count_4d(4),var1d_dim(1),var2d_dim(2),var3d_dim(3)
      real*8 :: dble2
      real*8,allocatable :: timeout(:),xnd(:),ynd(:),xout(:),yout(:)

      integer :: icomp_stats(3)

      print*, 'Input NODE-based output name (e.g. out2d):'
      read(*,'(a30)')file63
      print*, '<<<<<file name: ',file63

!      print*, 'Is the var node-based or ele-based? 0: node based; 1: element based'
!      read(*,*)inode_ele
!      print*, '<<<<<inode_ele: ',inode_ele
!      inode_ele=0
      if(trim(adjustl(file63)).eq."out2d") then !out2d
        i23d=1
      else !3D
        i23d=2 
      endif

!      varname=adjustl(varname); len_var=len_trim(varname)
      
      print*, 'Input start and end stack # to read:'
      read(*,*) iday1,iday2
      print*, '<<<<<start and end stack #: ',iday1,iday2

!...  Get geometry from local_to_global*
      open(10,file='local_to_global_000000',status='old')
      read(10,*)ns,ne,np,nvrt,nproc
      close(10)
      allocate(xnd(np),ynd(np),dp(np),kbp00(np),i34(ne),elnode(4,ne), &
     &np_lcl(0:nproc-1),ns_lcl(0:nproc-1),ne_lcl(0:nproc-1),nm2(4,ne), &
     &ztot(nvrt),sigma(nvrt),stat=istat)
      if(istat/=0) stop 'Failed to allocate'
    
      elnode=-99999 !init
      !-------------------------------------------------------------------------------
      ! Read rank-specific local_to_global*
      !-------------------------------------------------------------------------------
      ! Read in local-global mappings from all ranks
      file_char='local_to_global_000000'
      char_len=len_trim(file_char)
    
      !Find max. for dimensioning
      do irank=0,nproc-1
        write(file_char(char_len-5:char_len),'(i6.6)') irank
        open(10,file=file_char,status='old')
        read(10,*) !global info
        read(10,*) !info
        read(10,*)ne_lcl(irank) !resident only
        do i=1,ne_lcl(irank)
          read(10,*)!j,ielg(irank,i)
        enddo !i
        read(10,*)np_lcl(irank) !resident only (no halo)
        do i=1,np_lcl(irank)
          read(10,*)
        enddo !i
        read(10,*)ns_lcl(irank)
        close(10)
      enddo !irank
      np_max=maxval(np_lcl(:))
      ns_max=maxval(ns_lcl(:))
      ne_max=maxval(ne_lcl(:))
    
      if(allocated(iplg)) deallocate(iplg)
      if(allocated(ielg)) deallocate(ielg)
      if(allocated(islg)) deallocate(islg)
      if(allocated(iegl_rank)) deallocate(iegl_rank)
      allocate(iplg(0:nproc-1,np_max),ielg(0:nproc-1,ne_max),islg(0:nproc-1,ns_max), &
     &iegl_rank(ne),stat=istat)
      if(istat/=0) stop 'Allocation error (2)'
      iplg=-1 !init as junk
    
      !Re-read
      do irank=0,nproc-1
        write(file_char(char_len-5:char_len),'(i6.6)') irank
        open(10,file=file_char,status='old')
        read(10,*) !global info
        read(10,*) !info
        read(10,*)ne_lcl(irank)
        do i=1,ne_lcl(irank)
          read(10,*)j,ielg(irank,i)
          iegl_rank(ielg(irank,i))=irank
        enddo !i
        read(10,*)np_lcl(irank)
        do i=1,np_lcl(irank)
          read(10,*)j,iplg(irank,i)
        enddo
        read(10,*)ns_lcl(irank) !sides
        do i=1,ns_lcl(irank)
          read(10,*)j,islg(irank,i)
        enddo
    
        read(10,*) !'Header:'
        read(10,*) itmp,itmp,itmp,dble2,dble2 !start_year,start_month,start_day,start_hour,utc_start 
        read(10,*)nrec,dtout,nspool,nvrt,kz,h0,h_s,h_c,theta_b,theta_f,itmp !ics
        do k=1,kz-1
          read(10,*)ztot(k)
        enddo !k
        do k=1,nvrt-kz+1
          read(10,*)sigma(k)
        enddo

        read(10,*)np_lcl(irank),ne_lcl(irank),(xnd(iplg(irank,m)),ynd(iplg(irank,m)), &
      &dp(iplg(irank,m)),kbp00(iplg(irank,m)),m=1,np_lcl(irank)), &
      &(i34(ielg(irank,m)),(nm2(mm,m),mm=1,i34(ielg(irank,m))),m=1,ne_lcl(irank))
      !nm2 is local node index
    
        close(10)
    
        !Debug
        !write(98,*)irank,(i34(ielg(irank,m)),m=1,ne(irank))
    
        !Reconstruct connectivity table
        do m=1,ne_lcl(irank)
          iegb=ielg(irank,m)
          if(iegb>ne) stop 'get_global_geo: Overflow!'
          do mm=1,i34(iegb)
            itmp=nm2(mm,m)
            if(itmp>np_lcl(irank).or.itmp<=0) then
              write(*,*)'get_global_geo: Overflow,',m,mm,itmp
              stop
            endif
            elnode(mm,iegb)=iplg(irank,itmp)
          enddo !mm
        enddo !m
      enddo !irank=0,nproc-1

      deallocate(nm2,ztot,sigma)

      !Returned vars: ne,np,ns,nrec,[xnd ynd dp](np),
      !elnode,i34,nvrt,h0,dtout,kbp
!      call get_dims(iday1,np,ne,ns,nvrt,h0)
!      allocate(xnd(np),ynd(np),dp(np),kbp(np),i34(ne),elnode(4,ne),stat=istat)
!      if(istat/=0) stop 'alloc (1)'
!      call readheader(iday1,np,ne,ns,kbp,i34,elnode,nrec,xnd,ynd,dp,dtout)

      print*
      print*, 'After header:',ne,np,nrec,i34(ne),elnode(1:i34(ne),ne),nvrt,h0,xnd(np),ynd(np),dp(np) !,start_time

      !For dimensioning purpose
      if(np>ns.or.ne>ns) stop 'ns is not largest'

!      if (inode_ele==0) then
      last_dim=np
      ivs=1
!      elseif (inode_ele==1) then
!        last_dim=ne
!      endif

      allocate(outvar(nvrt,last_dim,ivs),eta2(np),idry(np),ztmp(nvrt,np),timeout(nrec), &
     &stat=istat)
      if(istat/=0) stop 'Allocation error (3)'

!     Prep output node dim: each rank outputs in turn
      np_out=sum(np_lcl)
      print*, '# of node=',np_out,np
      allocate(xout(np_out),yout(np_out),dpout(np_out),etaout(np_out),out3d(nvrt,np_out),stat=istat)
      if(istat/=0) stop 'Alloc (5)'
      out3d=-huge(1.0) !test mem
      icount=0
      do irank=0,nproc-1
        do i=1,np_lcl(irank)
          nd=iplg(irank,i)
          icount=icount+1
          if(icount>np_out) then
            print*, 'Oveflow(0):',irank,i,nd,icount,np_out
            stop
          endif
          xout(icount)=xnd(nd)
          yout(icount)=ynd(nd)
          dpout(icount)=dp(nd)
        enddo !i
      enddo !irank
      if(icount/=np_out) stop 'mismatch(1)'


!     Read vgrid.in
!      call get_vgrid_single('vgrid.in',np,nvrt,ivcor,kz,h_s,h_c,theta_b,theta_f,ztot,sigma,sigma_lcl,kbp)

!     Calculate kbp00
!      if(ivcor==1) then
!        kbp00=kbp
!      else
!        do i=1,np
!          call zcor_SZ_single(dp(i),1.e8,h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot,sigma, &
!     &ztmp(:,i),idry2,kbp00(i))
!        enddo !i
!      endif !ivcor
      
!...  Time iteration
!...
      outvar=-99 !init.
      ztmp=-99

      do iday=iday1,iday2 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      write(it_char,'(i12)')iday
      it_char=adjustl(it_char)
      leng=len_trim(it_char)

      file_char=trim(adjustl(file63))//'_'//it_char(1:leng)//'.nc'
      iret=nf90_open(trim(adjustl(file_char)),OR(NF90_NETCDF4,NF90_NOWRITE),ncid)
      if(iret/=nf90_NoErr) stop 'Failed to open file63'
      iret=nf90_inq_varid(ncid,'time',itime_id)
      iret=nf90_get_var(ncid,itime_id,timeout,(/1/),(/nrec/))
      print*, 'time=',timeout !,trim(adjustl(file63))

      if(i23d==1) then !out2d
        !Get IDs
        iret=nf90_inq_varid(ncid,'elevation',ieta_id)
      else !3D var
        iret=nf90_inq_varid(ncid,trim(adjustl(file63)),ivarid1)
        if(iret/=nf90_NoErr) stop 'Var not found'
        iret=nf90_Inquire_Variable(ncid,ivarid1,ndims=ndims,dimids=dimids)
        if(ndims>100) stop 'increase dimension of dimids & idims'
        do i=1,ndims
         iret=nf90_Inquire_Dimension(ncid,dimids(i),len=idims(i))
        enddo !i
        npes=idims(ndims-1) !np|ne|ns
        if(npes/=np) stop 'can only handle node-based'
!'
        if(idims(ndims)/=nrec) stop 'last dim is not time'

!      if(ivs==2) then !vector
!        iret=nf90_open(trim(adjustl(file64)),OR(NF90_NETCDF4,NF90_NOWRITE),ncid2)
!        if(iret/=nf90_NoErr) stop 'Failed to open file64'
!        iret=nf90_inq_varid(ncid2,varname2(1:len_var),ivarid2)
!        if(iret/=nf90_NoErr) stop 'Var2 not found'
!      endif !ivs

!        iret=nf90_get_att(ncid,ivarid1,'i23d',i23d)
!        if(i23d<=0.or.i23d>6) stop 'wrong i23d'
!        if(i23d>3.and.ics==2) stop 'ics=2 with elem-based var'
!        iret=nf90_get_att(ncid,ivarid1,'ivs',ivs)
!        !print*, 'i23d:',i23d,ivs,idims(1:ndims)
!
      endif !out2d or 3D

      !Open output
      file61=trim(adjustl(file63))//'_ro_'//it_char(1:leng)//'.nc'
      j=nf90_create(trim(adjustl(file61)),OR(NF90_NETCDF4,NF90_CLOBBER),ncid_out)
      if(j/=nf90_NoErr) stop 'Failed to open output'
      j=nf90_def_dim(ncid_out,'np_out',np_out,node_dim)
      j=nf90_def_dim(ncid_out,'nVert',nvrt,nvrt_dim)
      j=nf90_def_dim(ncid_out,'one',1,one_dim)
      j=nf90_def_dim(ncid_out,'numProssesors',nproc,nproc_dim)
      j=nf90_def_dim(ncid_out,'maxNodesPerRank',np_max,np_max_dim)
      j=nf90_def_dim(ncid_out,'time',nrec,time_dim)
      var1d_dim(1)=time_dim
      j=nf90_def_var(ncid_out,'time',NF90_DOUBLE,var1d_dim,itime_id2)
      if(i23d==1) then
        var1d_dim(1)=node_dim
        j=nf90_def_var(ncid_out,'SCHISM_hgrid_node_x',NF90_DOUBLE,var1d_dim,ix_id2)
        j=nf90_def_var(ncid_out,'SCHISM_hgrid_node_y',NF90_DOUBLE,var1d_dim,iy_id2)
        j=nf90_def_var(ncid_out,'depth',NF90_FLOAT,var1d_dim,ih_id2)
        var1d_dim(1)=nproc_dim
        j=nf90_def_var(ncid_out,'numNodesPerRank',NF90_INT,var1d_dim,inp_lcl_id2)
 
        var2d_dim(1)=nproc_dim; var2d_dim(2)=np_max_dim
        j=nf90_def_var(ncid_out,'node_local_to_global',NF90_INT,var2d_dim,iplg_id2)
        var2d_dim(1)=node_dim; var2d_dim(2)=time_dim
        j=nf90_def_var(ncid_out,'elevation',NF90_FLOAT,var2d_dim,ieta_id2)
        j=nf90_enddef(ncid_out)
        !Static
        start_1d(1)=1; count_1d(1)=np_out
        j=nf90_put_var(ncid_out,ix_id2,xout,start_1d,count_1d)
        j=nf90_put_var(ncid_out,iy_id2,yout,start_1d,count_1d)
        j=nf90_put_var(ncid_out,ih_id2,dpout,start_1d,count_1d)
        !Warning: the proc ID implicitly changed to 1:nproc
        start_1d(1)=1; count_1d(1)=nproc
        j=nf90_put_var(ncid_out,inp_lcl_id2,np_lcl(0:nproc-1),start_1d,count_1d)
  
        start_2d(1)=1; start_2d(2)=1 
        count_2d(1)=nproc; count_2d(2)=np_max
        j=nf90_put_var(ncid_out,iplg_id2,iplg(0:nproc-1,1:np_max),start_2d,count_2d)
      else !3D
        var3d_dim(1)=nvrt_dim; var3d_dim(2)=node_dim; var3d_dim(3)=time_dim
        j=nf90_def_var(ncid_out,trim(adjustl(file63)),NF90_FLOAT,var3d_dim,iout3d_id2)
        j=nf90_enddef(ncid_out)
      endif !i23d
      !put time
      j=nf90_put_var(ncid_out,itime_id2,timeout,(/1/),(/nrec/))
 
      do irec=1,nrec
        if(i23d==1) then !2D
          start_2d(1)=1; start_2d(2)=irec
          count_2d(1)=np; count_2d(2)=1
          iret=nf90_get_var(ncid,ieta_id,eta2,start_2d,count_2d)

!          start_2d(1)=1; start_2d(2)=irec
!          count_2d(1)=np; count_2d(2)=1
!          iret=nf90_get_var(ncid,,outvar(1,1:npes,1),start_2d,count_2d)
!          if(ivs==2) iret=nf90_get_var(ncid2,ivarid2,outvar(1,1:npes,2),start_2d,count_2d)

          !Output
          icount=0
          do irank=0,nproc-1
            do i=1,np_lcl(irank)
              nd=iplg(irank,i) 
              icount=icount+1
              if(icount>np_out) then
                print*, 'Oveflow(1):',irank,i,nd,icount,np_out
                stop 
              endif
              etaout(icount)=eta2(nd)
            enddo !i
          enddo !irank
          if(icount/=np_out) stop 'mismatch(1)'
          start_2d(1)=1; start_2d(2)=irec
          count_2d(1)=np_out; count_2d(2)=1
          j=nf90_put_var(ncid_out,ieta_id2,etaout,start_2d,count_2d)
        else !3D
          start_3d(1:2)=1; start_3d(3)=irec
          count_3d(1)=nvrt; count_3d(2)=np; count_3d(3)=1
          iret=nf90_get_var(ncid,ivarid1,outvar(:,1:npes,1),start_3d,count_3d)
!          if(ivs==2) iret=nf90_get_var(ncid2,ivarid2,outvar(:,1:npes,2),start_3d,count_3d)

          !Output
          icount=0
          do irank=0,nproc-1
            do i=1,np_lcl(irank)
              nd=iplg(irank,i)
              icount=icount+1
              if(icount>np_out) then
                print*, 'Oveflow(2):',irank,i,nd,icount,np_out
                stop
              endif
              out3d(:,icount)=outvar(:,nd,1)
            enddo !i
          enddo !irank
          if(icount/=np_out) stop 'mismatch(2)'
          start_3d(1:2)=1; start_3d(3)=irec
          count_3d(1)=nvrt; count_3d(2)=np_out; count_3d(3)=1
          j=nf90_put_var(ncid_out,iout3d_id2,out3d,start_3d,count_3d)
        endif !i23d
      enddo !irec

      iret=nf90_close(ncid)
      j=nf90_close(ncid_out)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      enddo !iday

      

      print*, 'Finished!'

      stop
      end
