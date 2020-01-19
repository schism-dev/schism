!  Basic routines for reading SCHISM netcdf outputs
!  WARNING: beware conflicts with global vars defined in this module, esp
!  scalars like start_time,nrec,h0,dtout!
!  Author: Joseph Zhang
!  Date: Sept 2017

!  Routines in this module
!  subroutine readheader (combined nc)
!  subroutine get_timeout (combined or uncombined nc)
!  subroutine get_elev (combined nc)
!  subroutine get_outvar (combined or uncombined nc)
!  subroutine get_global_geo (uncombined nc)
!  subroutine get_outvar_multirecord (combined or uncombined nc)

    module extract_mod
    use netcdf
    implicit none
    public

    character(len=48), save :: start_time 
    integer, save :: nrec,nspool,nvrt,kz,np,ne,ns,nproc
    real, save :: h0,fill_in,dtout
    real, save, allocatable :: x(:),y(:),dp(:)
    integer, save, allocatable :: kbp00(:),i34(:),elnode(:,:), &
  &iplg(:,:),ielg(:,:),islg(:,:),np_lcl(:),ne_lcl(:),ns_lcl(:),iegl_rank(:)

    integer :: char_len,start_2d(2),start_3d(3),start_4d(4), &
  &count_2d(2),count_3d(3),count_4d(4)
    real*8,allocatable :: timeout2(:)
    character(len=12) :: stack_char
    character(len=30) :: file_char

      contains

!================================================================
!     Read in static info from _combined_ nc
!     Returned vars: ne,np,ns,nrec,start_time,[x y dp](np),elnode,i34,nvrt,
!                    h0,dtout
!================================================================
      subroutine readheader(istack)
!      character(len=*),intent(in) :: fname !schout_*.nc (combined)
      integer, intent(in) :: istack
      integer :: i,varid1,varid2,dimids(3),istat,nvtx,iret,ncid2,itime_id,ielev_id

      write(stack_char,'(i12)')istack
      stack_char=adjustl(stack_char)
      char_len=len_trim(stack_char)
      file_char='schout_'//stack_char(1:char_len)//'.nc'

      iret=nf90_open(trim(adjustl(file_char)),OR(NF90_NETCDF4,NF90_NOWRITE),ncid2)
      iret=nf90_inq_dimid(ncid2,'nSCHISM_hgrid_edge',i)
      iret=nf90_Inquire_Dimension(ncid2,i,len=ns)
      iret=nf90_inq_dimid(ncid2,'nSCHISM_vgrid_layers',i)
      iret=nf90_Inquire_Dimension(ncid2,i,len=nvrt)
      iret=nf90_inq_varid(ncid2,'minimum_depth',varid1)
      iret=nf90_get_var(ncid2,varid1,h0)
      iret=nf90_inq_varid(ncid2,'SCHISM_hgrid_face_nodes',varid1)
      iret=nf90_Inquire_Variable(ncid2,varid1,dimids=dimids(1:2))
      iret=nf90_Inquire_Dimension(ncid2,dimids(1),len=nvtx)
      iret=nf90_Inquire_Dimension(ncid2,dimids(2),len=ne)
      if(nvtx/=4) stop 'readheader: vtx/=4'
      iret=nf90_inq_varid(ncid2,'SCHISM_hgrid_node_x',varid2)
      iret=nf90_Inquire_Variable(ncid2,varid2,dimids=dimids)
      iret=nf90_Inquire_Dimension(ncid2,dimids(1),len=np)
      iret=nf90_inq_varid(ncid2,'time',itime_id)
      iret=nf90_Inquire_Variable(ncid2,itime_id,dimids=dimids)
      iret=nf90_Inquire_Dimension(ncid2,dimids(1),len=nrec)
      iret=nf90_get_att(ncid2,itime_id,'base_date',start_time)
      iret=nf90_inq_varid(ncid2,'elev',ielev_id)
      if(iret.ne.NF90_NOERR) then
        print*, nf90_strerror(iret); stop 'readheader: error reading header'
      endif

      if(allocated(x)) deallocate(x)
      if(allocated(y)) deallocate(y)
      if(allocated(dp)) deallocate(dp)
      if(allocated(kbp00)) deallocate(kbp00)
      if(allocated(elnode)) deallocate(elnode)
      if(allocated(i34)) deallocate(i34)
      if(allocated(timeout2)) deallocate(timeout2)
      allocate(x(np),y(np),dp(np),kbp00(np),i34(ne),elnode(4,ne),timeout2(nrec),stat=istat)
      if(istat/=0) stop 'readheader: failed to allocate (3)'
      iret=nf90_get_var(ncid2,varid1,elnode)
      iret=nf90_get_var(ncid2,varid2,x)
      iret=nf90_get_var(ncid2,itime_id,timeout2,(/1/),(/nrec/))
      dtout=timeout2(2)-timeout2(1)

      iret=nf90_inq_varid(ncid2,'SCHISM_hgrid_node_y',varid1)
      iret=nf90_get_var(ncid2,varid1,y)
      iret=nf90_inq_varid(ncid2,'depth',varid1)
      iret=nf90_get_var(ncid2,varid1,dp)
      !iret=nf90_inq_varid(ncid2,'node_bottom_index',varid1)
      !iret=nf90_get_var(ncid2,varid1,kbp00)

      iret=nf90_close(ncid2)

      !print*, 'nc dim:',nvrt,np,ne,nrec,start_time

      !Calc i34
      i34=4 !init
      do i=1,ne
        if(elnode(4,i)<0) i34(i)=3
      enddo !i

      end subroutine readheader

!================================================================
!     Read in timeout from combined or uncombined nc
!     Returned vars: timeout()
!================================================================
      subroutine get_timeout(istack,nrec2,timeout,icomb)
      integer, intent(in) :: istack,nrec2  
      integer, optional, intent(in) :: icomb
      real*8, intent(out) :: timeout(nrec2) 
      integer :: iret,ncid2,itime_id
      real :: timeout5(nrec2)

      write(stack_char,'(i12)')istack
      stack_char=adjustl(stack_char)
      char_len=len_trim(stack_char)
      if(present(icomb)) then !uncombined
        file_char='schout_0000_'//stack_char(1:char_len)//'.nc'
      else
        file_char='schout_'//stack_char(1:char_len)//'.nc'
      endif
      iret=nf90_open(trim(adjustl(file_char)),OR(NF90_NETCDF4,NF90_NOWRITE),ncid2)
      !time is double for combined but single for uncombined!
      iret=nf90_inq_varid(ncid2,'time',itime_id)
      if(present(icomb)) then
        iret=nf90_get_var(ncid2,itime_id,timeout5,(/1/),(/nrec2/))
        timeout=timeout5
      else
        iret=nf90_get_var(ncid2,itime_id,timeout,(/1/),(/nrec2/))
      endif
!      print*, 'time=',timeout,trim(adjustl(file63))
      iret=nf90_close(ncid2)

      end subroutine get_timeout


!================================================================
!     Read in elev from _combined_ nc
!     Returned vars: eta2()
!================================================================
      subroutine get_elev(istack,irec,np2,eta2)
      integer, intent(in) :: istack,irec,np2
      real, intent(out) :: eta2(np2)
      integer :: iret,ncid2,ielev_id

      write(stack_char,'(i12)')istack
      stack_char=adjustl(stack_char)
      char_len=len_trim(stack_char)
      file_char='schout_'//stack_char(1:char_len)//'.nc'
      iret=nf90_open(trim(adjustl(file_char)),OR(NF90_NETCDF4,NF90_NOWRITE),ncid2)
      iret=nf90_inq_varid(ncid2,'elev',ielev_id)

      start_2d(1)=1; start_2d(2)=irec
      count_2d(1)=np2; count_2d(2)=1
      iret=nf90_get_var(ncid2,ielev_id,eta2,start_2d,count_2d)
      iret=nf90_close(ncid2)

      end subroutine get_elev

!================================================================
!     Read in the desired var (node/elem based) from combined or uncombined nc,
!     1 time record at a time.
!     Returned vars: outvar(2,nvrt,np|ne),i23d,ivs and eta2(np) (on global
!     indices). For uncombined nc, only subdomain part of arrays will be filled with values
!================================================================
      subroutine get_outvar(ics,istack,irec,varname,np2,last_dim,nvrt2,outvar,i23d,ivs,eta2,irank)
      !use ics=1 if linear interp
      integer, intent(in) :: ics,istack,irec,np2,last_dim,nvrt2
      character(len=*),intent(inout) :: varname
      real, intent(out) :: outvar(2,nvrt2,last_dim),eta2(np2)
      integer, intent(out) :: i23d,ivs
      integer, optional, intent(in) :: irank

      integer :: iret,i,npes,len_var,ivarid1,ndims,dimids(100),idims(100),ncid2,ielev_id
      character(len=4) :: a_4
      real, allocatable :: worka(:,:,:),workb(:)

      if(allocated(worka)) deallocate(worka) 
      if(allocated(workb)) deallocate(workb)

      varname=adjustl(varname); len_var=len_trim(varname)

      write(stack_char,'(i12)')istack
      stack_char=adjustl(stack_char)
      char_len=len_trim(stack_char)
      if(present(irank)) then !uncombined
        write(a_4,'(i4.4)') irank
        file_char='schout_'//a_4//'_'//stack_char(1:char_len)//'.nc'
      else
        file_char='schout_'//stack_char(1:char_len)//'.nc'
      endif
      print*, 'File=',trim(adjustl(file_char)),ics,istack,irec,varname
      iret=nf90_open(trim(adjustl(file_char)),OR(NF90_NETCDF4,NF90_NOWRITE),ncid2)
      if(iret/=nf90_NoErr) stop 'get_outvar: cannot open'
      iret=nf90_inq_varid(ncid2,'elev',ielev_id)
      if(iret/=nf90_NoErr) stop 'get_outvar: cannot find elev'

      start_2d(1)=1; start_2d(2)=irec
      count_2d(2)=1
      if(present(irank)) then
        count_2d(1)=np_lcl(irank) 
        allocate(workb(np_lcl(irank)))
        iret=nf90_get_var(ncid2,ielev_id,workb,start_2d,count_2d)
        eta2(iplg(irank,1:np_lcl(irank)))=workb
        deallocate(workb)
      else
        count_2d(1)=np2 
        iret=nf90_get_var(ncid2,ielev_id,eta2,start_2d,count_2d)
      endif

      iret=nf90_inq_varid(ncid2,varname(1:len_var),ivarid1)
      if(iret/=nf90_NoErr) stop 'get_outvar: Var not found'
      iret=nf90_Inquire_Variable(ncid2,ivarid1,ndims=ndims,dimids=dimids)
      if(ndims>100) stop 'get_outvar:increase dimension of dimids & idims'
      do i=1,ndims
        iret=nf90_Inquire_Dimension(ncid2,dimids(i),len=idims(i))
      enddo !i
      if(nvrt2/=nvrt) stop 'get_outvar: nvrt2/=nvrt'
      if(idims(ndims)/=nrec) stop 'get_outvar: last dim is not time'

      npes=idims(ndims-1) !np|ne|ns
      if(present(irank)) then
        if(npes/=np_lcl(irank).and.npes/=ne_lcl(irank)) stop 'get_outvar: can only handle node- or elem-based'
      else
        if(npes/=np.and.npes/=ne) stop 'get_outvar: can only handle node- or elem-based'
      endif

      allocate(worka(2,nvrt,npes))      

      iret=nf90_get_att(ncid2,ivarid1,'i23d',i23d)
      if(i23d<=0.or.i23d>6) stop 'get_outvar: wrong i23d'
      if(i23d>3.and.ics==2) stop 'get_outvar: ics=2 with elem-based var'
      iret=nf90_get_att(ncid2,ivarid1,'ivs',ivs)
      !print*, 'i23d:',i23d,ivs,idims(1:ndims)

      if(ivs==1) then !scalar
        if(mod(i23d-1,3)==0) then !2D
          start_2d(1)=1; start_2d(2)=irec
          count_2d(1)=npes; count_2d(2)=1
          iret=nf90_get_var(ncid2,ivarid1,worka(1,1,1:npes),start_2d,count_2d)
        else !3D
          start_3d(1:2)=1; start_3d(3)=irec
          count_3d(2)=npes; count_3d(1)=nvrt; count_3d(3)=1
          iret=nf90_get_var(ncid2,ivarid1,worka(1,:,1:npes),start_3d,count_3d)
        endif 
      else !vector
        if(mod(i23d-1,3)==0) then !2D
          start_3d(1:2)=1; start_3d(3)=irec
          count_3d(2)=npes; count_3d(1)=2; count_3d(3)=1
          iret=nf90_get_var(ncid2,ivarid1,worka(1:2,1,1:npes),start_3d,count_3d)
        else if(ndims-1==3) then !3D vector
          start_4d(1:3)=1; start_4d(4)=irec
          count_4d(3)=npes; count_4d(2)=nvrt; count_4d(1)=2; count_4d(4)=1
          iret=nf90_get_var(ncid2,ivarid1,worka(:,:,1:npes),start_4d,count_4d)
        else
          stop 'get_outvar: Unknown type(2)'
        endif
      endif !ivs

      iret=nf90_close(ncid2)

      if(present(irank)) then !uncombined; not all may be filled
        if(i23d<=3) then !node
          outvar(:,:,iplg(irank,1:np_lcl(irank)))=worka
        else if(i23d<=6) then !elem
          outvar(:,:,ielg(irank,1:ne_lcl(irank)))=worka
        else
          stop 'get_outvar: cannot be side based'
        endif
      else
        outvar(:,:,1:npes)=worka
      endif

      deallocate(worka)

      end subroutine get_outvar
      
!================================================================
!     Gather global geometry from uncombined nc outputs
!     local_to_global*
!     Returned vars: ne,np,ns,nrec,[x y dp](np),elnode,i34,nvrt,h0,dtout, and
!                    iegl_rank
!================================================================
      subroutine get_global_geo

      integer :: istat,i,j,k,irank,iegb,m,mm,itmp,np_max,ne_max,ns_max
      integer, allocatable :: nm2(:,:)
      real :: h_s,h_c,theta_b,theta_f
      real*8 :: dble2
      real, allocatable :: ztot(:),sigma(:)

      !Read local_to_global_0000 for global info
      open(10,file='local_to_global_0000',status='old')
      read(10,*)ns,ne,np,nvrt,nproc
      close(10)
      if(allocated(x)) deallocate(x)
      if(allocated(y)) deallocate(y)
      if(allocated(dp)) deallocate(dp)
      if(allocated(kbp00)) deallocate(kbp00)
      if(allocated(elnode)) deallocate(elnode)
      if(allocated(i34)) deallocate(i34)
      if(allocated(np_lcl)) deallocate(np_lcl)
      if(allocated(ne_lcl)) deallocate(ne_lcl)
      if(allocated(ns_lcl)) deallocate(ns_lcl)
      if(allocated(nm2)) deallocate(nm2)
      if(allocated(ztot)) deallocate(ztot)
      if(allocated(sigma)) deallocate(sigma)
    
      allocate(x(np),y(np),dp(np),kbp00(np),i34(ne),elnode(4,ne), &
  &np_lcl(0:nproc-1),ns_lcl(0:nproc-1),ne_lcl(0:nproc-1),nm2(4,ne), &
  &ztot(nvrt),sigma(nvrt),stat=istat)
      if(istat/=0) stop 'get_global_geo: failed to allocate'
    
      elnode=-99999 !init
      !-------------------------------------------------------------------------------
      ! Read rank-specific local_to_global*
      !-------------------------------------------------------------------------------
      ! Read in local-global mappings from all ranks
      file_char='local_to_global_0000'
      char_len=len_trim(file_char)
    
      !Find max. for dimensioning
      do irank=0,nproc-1
        write(file_char(char_len-3:char_len),'(i4.4)') irank
        open(10,file=file_char,status='old')
        read(10,*) !global info
        read(10,*) !info
        read(10,*)ne_lcl(irank)
        do i=1,ne_lcl(irank)
          read(10,*)!j,ielg(irank,i)
        enddo !i
        read(10,*)np_lcl(irank)
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
      if(istat/=0) stop 'get_global_geo: Allocation error (2)'
    
      !Re-read
      do irank=0,nproc-1
        write(file_char(char_len-3:char_len),'(i4.4)') irank
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
        read(10,*)(ztot(k),k=1,kz-1),(sigma(k),k=1,nvrt-kz+1)
        read(10,*)np_lcl(irank),ne_lcl(irank),(x(iplg(irank,m)),y(iplg(irank,m)), &
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

      end subroutine get_global_geo

!================================================================
!     Read in the desired var (node/elem based) from combined or uncombined nc
!     for multiple time records [irec1,irec2]. This routine requires storage of
!     data in some memory-intensive arrays, so the driving routine will ask for
!     either max array size (usually used for uncombined) or max # of records (
!     usually for combined) allowed. 
!     Returned vars: outvar(2,nvrt,np|ne,nrec3),i23d,ivs and eta2(np,nrec3) (on global
!     indices). For uncombined nc, only subdomain part of arrays will be filled with values
!================================================================
      subroutine get_outvar_multirecord(ics,istack,varname,irec1,irec2,np2,last_dim,nvrt2,nrec3,outvar,i23d,ivs,eta2,irank)
      !use ics=1 if linear interp
      integer, intent(in) :: ics,istack,irec1,irec2,np2,last_dim,nvrt2,nrec3
      character(len=*),intent(inout) :: varname
      real, intent(out) :: outvar(2,nvrt2,last_dim,nrec3),eta2(np2,nrec3)
      integer, intent(out) :: i23d,ivs
      integer, optional, intent(in) :: irank

      integer :: iret,i,npes,len_var,ivarid1,ndims,dimids(100),idims(100),ncid2,ielev_id,nrec4
      character(len=4) :: a_4
      real, allocatable :: worka(:,:,:,:),workb(:,:)

      if(allocated(worka)) deallocate(worka)
      if(allocated(workb)) deallocate(workb)

      nrec4=irec2-irec1+1
      if(nrec4<=0.or.nrec4>nrec3) stop 'get_outvar_multirecord: nrec4<=0'
      if(np2/=np.or.nvrt2/=nvrt) stop 'get_outvar_multirecord: check inputs'
      varname=adjustl(varname); len_var=len_trim(varname)

      write(stack_char,'(i12)')istack
      stack_char=adjustl(stack_char)
      char_len=len_trim(stack_char)
      if(present(irank)) then !uncombined
        write(a_4,'(i4.4)') irank
        file_char='schout_'//a_4//'_'//stack_char(1:char_len)//'.nc'
      else
        file_char='schout_'//stack_char(1:char_len)//'.nc'
      endif
      iret=nf90_open(trim(adjustl(file_char)),OR(NF90_NETCDF4,NF90_NOWRITE),ncid2)
      iret=nf90_inq_varid(ncid2,'elev',ielev_id)

      start_2d(1)=1; start_2d(2)=irec1
      count_2d(2)=nrec4
      if(present(irank)) then
        count_2d(1)=np_lcl(irank) 
        allocate(workb(np_lcl(irank),nrec4))
        iret=nf90_get_var(ncid2,ielev_id,workb,start_2d,count_2d)
        eta2(iplg(irank,1:np_lcl(irank)),1:nrec4)=workb
        deallocate(workb)
      else
        count_2d(1)=np2 
        iret=nf90_get_var(ncid2,ielev_id,eta2(:,1:nrec4),start_2d,count_2d)
      endif

      iret=nf90_inq_varid(ncid2,varname(1:len_var),ivarid1)
      if(iret/=nf90_NoErr) stop 'get_outvar_multirecord: Var not found'
      iret=nf90_Inquire_Variable(ncid2,ivarid1,ndims=ndims,dimids=dimids)
      if(ndims>100) stop 'get_outvar_multirecord:increase dimension of dimids & idims'
      do i=1,ndims
        iret=nf90_Inquire_Dimension(ncid2,dimids(i),len=idims(i))
      enddo !i
      if(idims(ndims)/=nrec) stop 'get_outvar_multirecord: last dim is not time'

      npes=idims(ndims-1) !np|ne|ns
      if(present(irank)) then
        if(npes/=np_lcl(irank).and.npes/=ne_lcl(irank)) stop 'get_outvar_multirecord: can only handle node- or elem-based'
      else
        if(npes/=np.and.npes/=ne) stop 'get_outvar_multirecord: can only handle node- or elem-based'
      endif

      allocate(worka(2,nvrt,npes,nrec4))
      worka(2,nvrt,npes,nrec4)=0 !test mem

      iret=nf90_get_att(ncid2,ivarid1,'i23d',i23d)
      if(i23d<=0.or.i23d>6) stop 'get_outvar_multirecord: wrong i23d'
      if(i23d>3.and.ics==2) stop 'get_outvar_multirecord: ics=2 with elem-based var'
      iret=nf90_get_att(ncid2,ivarid1,'ivs',ivs)
      !print*, 'i23d:',i23d,ivs,idims(1:ndims)

      if(ivs==1) then !scalar
        if(mod(i23d-1,3)==0) then !2D
          start_2d(1)=1; start_2d(2)=irec1
          count_2d(1)=npes; count_2d(2)=nrec4
          iret=nf90_get_var(ncid2,ivarid1,worka(1,1,1:npes,:),start_2d,count_2d)
        else !3D
          start_3d(1:2)=1; start_3d(3)=irec1
          count_3d(1)=nvrt; count_3d(2)=npes; count_3d(3)=nrec4
          iret=nf90_get_var(ncid2,ivarid1,worka(1,:,1:npes,:),start_3d,count_3d)
        endif 
      else !vector
        if(mod(i23d-1,3)==0) then !2D
          start_3d(1:2)=1; start_3d(3)=irec1
          count_3d(1)=2; count_3d(2)=npes; count_3d(3)=nrec4
          iret=nf90_get_var(ncid2,ivarid1,worka(1:2,1,1:npes,:),start_3d,count_3d)
        else if(ndims-1==3) then !3D vector
          start_4d(1:3)=1; start_4d(4)=irec1
          count_4d(1)=2; count_4d(2)=nvrt; count_4d(3)=npes; count_4d(4)=nrec4
          iret=nf90_get_var(ncid2,ivarid1,worka(:,:,1:npes,:),start_4d,count_4d)
        else
          stop 'get_outvar_multirecord: Unknown type(2)'
        endif
      endif !ivs

      iret=nf90_close(ncid2)

      if(present(irank)) then !uncombined; not all may be filled
        if(i23d<=3) then !node
          outvar(:,:,iplg(irank,1:np_lcl(irank)),1:nrec4)=worka
        else if(i23d<=6) then !elem
          outvar(:,:,ielg(irank,1:ne_lcl(irank)),1:nrec4)=worka
        else
          stop 'get_outvar_multirecord: cannot be side based'
        endif
      else
        outvar(:,:,1:npes,1:nrec4)=worka
      endif

      deallocate(worka)

      end subroutine get_outvar_multirecord
      
!================================================================
!================================================================
    end module extract_mod
