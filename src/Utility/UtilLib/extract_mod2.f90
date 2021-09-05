!  Basic routines for reading SCHISM netcdf outputs
!  Author: Joseph Zhang
!  Date: Sept 2017

!  Routines in this module, for scribe I/O
!  subroutine readheader 

    module extract_mod2
    use netcdf
    implicit none
    private

!    character(len=48), save :: start_time 
!    integer, save :: nrec,nspool,nvrt,kz,np,ne,ns,nproc
    real, save :: fill_in
    !coordinates in double
!    real*8, save, allocatable :: xnd(:),ynd(:)
!    real, save, allocatable :: dp(:)
!    integer, save, allocatable :: kbp00(:),i34(:),elnode(:,:), &
!  &iplg(:,:),ielg(:,:),islg(:,:),np_lcl(:),ne_lcl(:),ns_lcl(:),iegl_rank(:)

    integer :: char_len,start_2d(2),start_3d(3),start_4d(4), &
  &count_2d(2),count_3d(3),count_4d(4)
!    real*8,allocatable :: timeout2(:)
    character(len=12) :: stack_char
    character(len=30) :: file_char

    public :: get_dims
    public :: readheader

    contains

!================================================================
!     Get dimensions from out2d_*.nc
      subroutine get_dims(istack,np,ne,ns,nvrt,h0)
      integer, intent(in) :: istack
      integer, intent(out) :: np,ne,ns,nvrt
      real, intent(out) :: h0

      integer :: i,varid1,varid2,dimids(3),istat,nvtx,iret,ncid2,itime_id,ielev_id
      real*8 :: tmp

      write(stack_char,'(i12)')istack
      stack_char=adjustl(stack_char)
      char_len=len_trim(stack_char)
      file_char='out2d_'//stack_char(1:char_len)//'.nc'

      iret=nf90_open(trim(adjustl(file_char)),OR(NF90_NETCDF4,NF90_NOWRITE),ncid2)
      if(iret.ne.NF90_NOERR) then
        print*, nf90_strerror(iret); stop 'extract_mod2: (1)'
      endif
      iret=nf90_inq_dimid(ncid2,'nSCHISM_hgrid_node',i)
      iret=nf90_Inquire_Dimension(ncid2,i,len=np)
      iret=nf90_inq_dimid(ncid2,'nSCHISM_hgrid_face',i)
      iret=nf90_Inquire_Dimension(ncid2,i,len=ne)
      iret=nf90_inq_dimid(ncid2,'nSCHISM_hgrid_edge',i)
      iret=nf90_Inquire_Dimension(ncid2,i,len=ns)
      iret=nf90_inq_dimid(ncid2,'nSCHISM_vgrid_layers',i)
      iret=nf90_Inquire_Dimension(ncid2,i,len=nvrt)
      iret=nf90_inq_varid(ncid2,'minimum_depth',varid1)
      iret=nf90_get_var(ncid2,varid1,tmp)
      h0=tmp
      iret=nf90_close(ncid2)
      
      end subroutine get_dims

!================================================================
!     Read in static info from out2d*nc
!================================================================
      subroutine readheader(istack,np,ne,ns,kbp00,i34,elnode,nrec,xnd,ynd,dp,dtout)
      integer, intent(in) :: istack,np,ne,ns
      integer, intent(out) :: kbp00(np),i34(ne),elnode(4,ne),nrec
      real*8, intent(out) :: xnd(np),ynd(np)
      real, intent(out) :: dp(np),dtout

      integer :: i,varid1,varid2,dimids(3),istat,nvtx,iret,ncid2,itime_id,ielev_id
      real*8,allocatable :: timeout2(:)

      write(stack_char,'(i12)')istack
      stack_char=adjustl(stack_char)
      char_len=len_trim(stack_char)
      file_char='out2d_'//stack_char(1:char_len)//'.nc'

      iret=nf90_open(trim(adjustl(file_char)),OR(NF90_NETCDF4,NF90_NOWRITE),ncid2)
      if(iret.ne.NF90_NOERR) then
        print*, nf90_strerror(iret); stop 'extract_mod2: (2)'
      endif
      iret=nf90_inq_varid(ncid2,'SCHISM_hgrid_face_nodes',varid1)
      iret=nf90_Inquire_Variable(ncid2,varid1,dimids=dimids(1:2))
      iret=nf90_Inquire_Dimension(ncid2,dimids(1),len=nvtx)
      if(nvtx/=4) then
        print*, 'readheader:',nvtx
        stop 'readheader: vtx/=4'
      endif
      iret=nf90_inq_varid(ncid2,'time',itime_id)
      iret=nf90_Inquire_Variable(ncid2,itime_id,dimids=dimids)
      iret=nf90_Inquire_Dimension(ncid2,dimids(1),len=nrec)
!      iret=nf90_get_att(ncid2,itime_id,'base_date',start_time)
!      iret=nf90_inq_varid(ncid2,'elev',ielev_id)
!      if(iret.ne.NF90_NOERR) then
!        print*, nf90_strerror(iret); stop 'readheader: error reading header'
!      endif

!      if(allocated(xnd)) deallocate(xnd)
!      if(allocated(ynd)) deallocate(ynd)
!      if(allocated(dp)) deallocate(dp)
!      if(allocated(kbp00)) deallocate(kbp00)
!      if(allocated(elnode)) deallocate(elnode)
!      if(allocated(i34)) deallocate(i34)
      allocate(timeout2(nrec),stat=istat)
      if(istat/=0) stop 'readheader: failed to allocate (3)'
      iret=nf90_get_var(ncid2,varid1,elnode)
      iret=nf90_get_var(ncid2,itime_id,timeout2,(/1/),(/nrec/))
      dtout=timeout2(2)-timeout2(1)

      iret=nf90_inq_varid(ncid2,'SCHISM_hgrid_node_x',varid1)
      iret=nf90_get_var(ncid2,varid1,xnd)
      iret=nf90_inq_varid(ncid2,'SCHISM_hgrid_node_y',varid1)
      iret=nf90_get_var(ncid2,varid1,ynd)
      iret=nf90_inq_varid(ncid2,'depth',varid1)
      iret=nf90_get_var(ncid2,varid1,dp)
      iret=nf90_inq_varid(ncid2,'bottom_index_node',varid1)
      iret=nf90_get_var(ncid2,varid1,kbp00)
      iret=nf90_close(ncid2)

      !print*, 'nc dim:',nvrt,np,ne,nrec,start_time

      !Calc i34
      i34=4 !init
      do i=1,ne
        if(elnode(4,i)<0) i34(i)=3
      enddo !i

      deallocate(timeout2)

      end subroutine readheader
      
!================================================================
!================================================================
    end module extract_mod2
