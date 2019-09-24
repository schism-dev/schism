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
!===============================================================================
! SCHISM FILE I/O SUBROUTINES
!
! subroutine write_obe
! subroutine report_timers
! subroutine writeout_nc
! subroutine fill_nc_header

!===============================================================================
!===============================================================================

    module schism_io
    use schism_glbl
    use schism_msgp
    use netcdf
    implicit none
!    include 'netcdf.inc'
    private

    integer,save :: ncid,node_dim,nele_dim,nedge_dim,four_dim,nv_dim, &
    &one_dim,two_dim,time_dim,time_dims(1),itime_id,ele_dims(2),x_dims(1), &
    &y_dims(1),z_dims(1),var2d_dims(2),var3d_dims(3),var4d_dims(4),dummy_dim(1), &
    &data_start_1d(1),data_start_2d(2),data_start_3d(3),data_start_4d(4), &
    &data_count_1d(1),data_count_2d(2),data_count_3d(3),data_count_4d(4)

    public :: write_obe
    public :: report_timers
    public :: writeout_nc
    public :: fill_nc_header

    contains

      subroutine write_obe
!-------------------------------------------------------------------------------
! Output centers.bp and sidecenters.bp
! NOTE: Valid for single processor only!
!-------------------------------------------------------------------------------
!      use schism_glbl
!      use schism_msgp
      implicit none
      integer :: i,j
      real(rkind) ::  tmp1,tmp2

      !Output sidecenters.bp
      open(32,file='sidecenters.bp')
      write(32,*) 'Sidegrid'
      write(32,*) ns
      do i=1,ns
        if(ics==1) then
          write(32,*) i,real(xcj(i)),real(ycj(i)),real(dps(i))
        else
          tmp1=sum(xlon(isidenode(1:2,i)))/2
          tmp2=sum(ylat(isidenode(1:2,i)))/2
          write(32,*) i,real(tmp1),real(tmp2),real(dps(i))
        endif
      enddo !ics
      close(32)

      !Output centers.bp
      open(32,file='centers.bp')
      write(32,*) 'centers pts'
      write(32,*) ne
      do i=1,ne
        if(ics==1) then
          write(32,'(i12,2(1x,e22.14),e12.4)') i,xctr(i),yctr(i),dpe(i)
        else
          tmp1=sum(xlon(elnode(1:i34(i),i)))/i34(i)
          tmp2=sum(ylat(elnode(1:i34(i),i)))/i34(i)
          write(32,*) i,real(tmp1),real(tmp2),real(dpe(i))
        endif !ics
      enddo !i
      close(32)

      end subroutine write_obe


!===============================================================================
!===============================================================================

      subroutine report_timers
!-------------------------------------------------------------------------------
! Write timing data for all tasks to timer.out file
!-------------------------------------------------------------------------------
!      use schism_glbl, only : rkind,mxtimer,wtimer
!      use schism_msgp
      implicit none
      include 'mpif.h'
      integer :: i
      real(rkind) :: wavg1(0:mxtimer),wavg2(0:mxtimer)
      real(rkind) :: wbuf(2,0:mxtimer)
      real(rkind) :: wmin1(2,0:mxtimer),wmin2(2,0:mxtimer)
      real(rkind) :: wmax1(2,0:mxtimer),wmax2(2,0:mxtimer)
!-------------------------------------------------------------------------------

      ! Sum communication time for timestepping section
      do i=3,13; wtimer(2,2)=wtimer(2,2)+wtimer(i,2); enddo;

      ! Total communication time
      wtimer(0,2)=wtimer(1,2)+wtimer(2,2)

      ! Make computation time excluding communication time
      wtimer(:,1)=wtimer(:,1)-wtimer(:,2)

      ! Compute average time over all tasks
      call mpi_allreduce(wtimer(0,1),wavg1,mxtimer+1,rtype,MPI_SUM,comm,ierr)
      wavg1=wavg1/dble(nproc)
      call mpi_allreduce(wtimer(0,2),wavg2,mxtimer+1,rtype,MPI_SUM,comm,ierr)
      wavg2=wavg2/dble(nproc)

      ! Compute min & max computation time over all tasks
      wbuf(1,:)=wtimer(:,1); wbuf(2,:)=myrank;
      call mpi_allreduce(wbuf,wmin1,mxtimer+1,MPI_2DOUBLE_PRECISION,MPI_MINLOC,comm,ierr)
      call mpi_allreduce(wbuf,wmax1,mxtimer+1,MPI_2DOUBLE_PRECISION,MPI_MAXLOC,comm,ierr)

      ! Compute min & max communication time over all tasks
      wbuf(1,:)=wtimer(:,2); wbuf(2,:)=myrank;
      call mpi_allreduce(wbuf,wmin2,mxtimer+1,MPI_2DOUBLE_PRECISION,MPI_MINLOC,comm,ierr)
      call mpi_allreduce(wbuf,wmax2,mxtimer+1,MPI_2DOUBLE_PRECISION,MPI_MAXLOC,comm,ierr)

      ! Rank 0 create new file and write header and avg time
      if(myrank==0) then
        ! Open new file
        open(36,file='timer.out',form='formatted',status='replace')

        ! Write ledger
        write(36,'(2a)') '# ','********** Timer Index Mapping **********'
        write(36,'(2a)') '# ','00 -- Total'
        write(36,'(2a)') '# ','01 -- Init Section'
        write(36,'(2a)') '# ','02 -- Timestepping Section'
        write(36,'(2a)') '# ','03 -- Forcings & Prep Section'
        write(36,'(2a)') '# ','04 -- Backtracking Section'
        write(36,'(2a)') '# ','05 -- Turbulence Closure Section'
        write(36,'(2a)') '# ','06 -- Matrix Preparation Section'
        write(36,'(2a)') '# ','07 -- Wave-Cont. Solver Section'
        write(36,'(2a)') '# ','08 -- Momentum Eqs. Solve Section'
        write(36,'(2a)') '# ','09 -- Transport Section'
        write(36,'(2a)') '# ','10 -- Recomputing Levels Section'
        write(36,'(2a)') '# ','11 -- Conservation Check Section'
        write(36,'(2a)') '# ','12 -- Global Output Section'
        write(36,'(2a)') '# ','13 -- Hotstart Section'
!'

        ! Write average, min & max times
        write(36,'(/)')
        write(36,'(2a)') '# ','********** Average, Min & Max Times in secs **********'
!'
        write(36,'(11a)') 'ID', &
     '        CompT','     MinCompT',' RankMinCompT','     MaxCompT',' RankMaxCompT', &
     '        CommT','     MinCommT',' RankMinCommT','     MaxCommT',' RankMaxCommT'
        do i=0,13
          write(36,'(i2.2,2(e13.6,2(e13.6,i13)))') i, &
          wavg1(i),wmin1(1,i),int(wmin1(2,i)),wmax1(1,i),int(wmax1(2,i)), &
          wavg2(i),wmin2(1,i),int(wmin2(2,i)),wmax2(1,i),int(wmax2(2,i))
        enddo

        ! Close file
        if(nproc>1) close(36)

      endif !myrank=0

      ! Initiate round-robin synchronization
      call parallel_rrsync(1)

      ! Open file to append
      if(nproc>1) &
        open(36,file='timer.out',form='formatted',status='old',position='append')

      ! Task 0 write next header
      if(myrank==0) then
        write(36,'(/)')
        write(36,'(a)') '# ********** Computation Times (sec) For Each MPI Task **********'
        write(36,'(a)') '# ********** Rows = Ranks; Columns = Timers      **********'
        write(36,'(a,20i13)') '# Rank',(i,i=0,13)
      endif

      ! Each task write its own timer data
      write(36,'(a,i4.4,20e13.6)') '# ',myrank,(wtimer(i,1),i=0,13)

      ! Close file
      if(nproc>1) close(36)

      ! Round-robin synchronization step
      call parallel_rrsync(2)

      ! Open file to append
      if(nproc>1) &
        open(36,file='timer.out',form='formatted',status='old',position='append')

      ! Task 0 write next header
      if(myrank==0) then
        write(36,'(/)')
        write(36,'(a)') '# ********** Communication Times For Each MPI Task **********'
        write(36,'(a)') '# ********** Rows = Ranks; Columns = Timers        **********'
        write(36,'(a,20i13)') '# Rank',(i,i=0,13)
      endif

      ! Each task write its own timer data
      write(36,'(a,i4.4,20e13.6)') '# ',myrank,(wtimer(i,2),i=0,13)

      ! Close file
      if(nproc>1) close(36)

      ! Round-robin synchronization final
      call parallel_rrsync(-2)

      end subroutine report_timers

!===============================================================================
!===============================================================================
      subroutine writeout_nc(varid,var_nm,i23d,idim1,idim2,outvar1,outvar2)
!-------------------------------------------------------------------------------
!     Netcdf outputs for global arrays. Can be called from any routine, but make sure that
!     the calling routine is called inside the main time loop 
!     exactly ONCE per step! 
!
!     Inputs:
!            var_nm: name of the output variable (to appear in nc file). 
!            i23d: indicates location where outputs are defined. 1:3 - node 2D/3D whole/3D half level
!                  4:6 - elem 2D/3D whole/half levels; 7:9 - side 2D/3D whole/half levels
!            idim1,idim2: dimensions of output array(s) in the driving routine. 
!                         For 2D variables (e.g., bottom
!                         stress), idim1 must be 1; for 3D variables, idim1 must be nvrt.
!                         idim2 must be consistent with the type of output as given by
!                         i23d (e.g., idim2=nea or ne for i23d=4);
!            outvar[1,2](idim1,idim2): output array. outvar2 is optional [for vectors]
!     In/out: varid: 1st call will generate variable ID, which is used later
!-------------------------------------------------------------------------------
      implicit none
   
      character(len=*),intent(in) :: var_nm
      integer,intent(in) :: i23d,idim1,idim2
      real(rkind),intent(in) :: outvar1(idim1,idim2)
      real(rkind),optional,intent(in) :: outvar2(idim1,idim2)
      integer,intent(inout) :: varid
 
      !character(len=3) :: sfix
      character(len=1000) :: var_nm2
      logical :: lex1,lex2
      integer :: i,k,iret,irec,len_var,idim2p,iret2,ivs
      real*4 :: a1d(1) 
      
!     Return if not output step
      if(mod(it_main,nspool)/=0) return

      ivs=1
      if(present(outvar2)) ivs=2

      irec=(it_main-(ifile-1)*ihfskip)/nspool !time recod #
      if(irec<=0) call parallel_abort('writeout_nc: irec<=0')
      var_nm2=var_nm
      var_nm2=adjustl(var_nm2); len_var=len_trim(var_nm2)

      !Define dim/vars
      !nf90_put_var(ncid,varid,values,start,count,stride)
      !values can be of any type, (optional) start, count, stride are of same dim
      !as values. e.g., to write to 1st entry in an array, start=count=1

      !Dump time
      !Note: using scalar directly won't work; must use array
      a1d(1)=real(time_stamp)
      data_start_1d(1)=irec; data_count_1d(1)=1
      iret=nf90_put_var(ncid,itime_id,a1d,data_start_1d,data_count_1d)

      !Use original dim order in nc
      if(i23d<=3) then !node
        var2d_dims(1)=node_dim
        var3d_dims(2)=node_dim
        var4d_dims(3)=node_dim
        idim2p=np !final output dim
      else if(i23d<=6) then !elem
        var2d_dims(1)=nele_dim
        var3d_dims(2)=nele_dim
        var4d_dims(3)=nele_dim
        idim2p=ne
      else if(i23d<=9) then !side
        var2d_dims(1)=nedge_dim
        var3d_dims(2)=nedge_dim
        var4d_dims(3)=nedge_dim
        idim2p=ns
      else
        call parallel_abort('writeout_nc: unknown i23d')       
      endif

      iret2=nf90_inq_varid(ncid,var_nm2(1:len_var),i)

      if(mod(i23d-1,3)==0) then !2D var (2D array in nc that has time dim)
        if(iret2/=NF90_NOERR) then !not defined yet
          iret=nf90_redef(ncid)
          if(ivs==1) then
            var2d_dims(2)=time_dim
            iret=nf90_def_var(ncid,var_nm2(1:len_var),NF90_FLOAT,var2d_dims,varid)
          else
            var3d_dims(1)=two_dim; var3d_dims(3)=time_dim
            iret=nf90_def_var(ncid,var_nm2(1:len_var),NF90_FLOAT,var3d_dims,varid)
          endif !ivs
          iret=nf90_put_att(ncid,varid,'i23d',i23d)
          iret=nf90_put_att(ncid,varid,'ivs',ivs)
          iret=nf90_enddef(ncid)
        endif !iret

        if(ivs==1) then
          data_start_2d(1)=1; data_start_2d(2)=irec
          data_count_2d(1)=idim2p; data_count_2d(2)=1
          iret=nf90_put_var(ncid,varid,real(outvar1(1,1:idim2p)),data_start_2d,data_count_2d)
        else !vector
          data_start_3d(1)=1; data_start_3d(2)=1; data_start_3d(3)=irec
          data_count_3d(1)=1; data_count_3d(2)=idim2p; data_count_3d(3)=1
          iret=nf90_put_var(ncid,varid,real(outvar1(1,1:idim2p)),data_start_3d,data_count_3d)
          data_start_3d(1)=2
          iret=nf90_put_var(ncid,varid,real(outvar2(1,1:idim2p)),data_start_3d,data_count_3d)
        endif !ivs
        !write(12,*)'2D:',it_main,varid,var_nm2(1:len_var),iret2
      else !3D
        if(iret2/=NF90_NOERR) then !not defined yet
          iret=nf90_redef(ncid)
          if(ivs==1) then
            var3d_dims(1)=nv_dim; var3d_dims(3)=time_dim
            iret=nf90_def_var(ncid,var_nm2(1:len_var),NF90_FLOAT,var3d_dims,varid)
          else
            var4d_dims(1)=two_dim; var4d_dims(2)=nv_dim; var4d_dims(4)=time_dim
            iret=nf90_def_var(ncid,var_nm2(1:len_var),NF90_FLOAT,var4d_dims,varid)
          endif !ivs
          !write(12,*)'3D def:',var3d_dims,varid,var_nm2(1:len_var),iret
!Add chunking option as well?
          iret=nf90_def_var_deflate(ncid,varid,0,1,4)
          iret=nf90_put_att(ncid,varid,'i23d',i23d)
          iret=nf90_put_att(ncid,varid,'ivs',ivs)
          iret=nf90_enddef(ncid)
        endif !iret

        if(ivs==1) then
          data_start_3d(1)=1; data_start_3d(2)=1; data_start_3d(3)=irec
          data_count_3d(1)=nvrt; data_count_3d(2)=idim2p; data_count_3d(3)=1
          iret=nf90_put_var(ncid,varid,real(outvar1(:,1:idim2p)),data_start_3d,data_count_3d)
        else !vector
          data_start_4d(1:3)=1; data_start_4d(4)=irec
          data_count_4d(1)=1; data_count_4d(2)=nvrt; data_count_4d(3)=idim2p; data_count_4d(4)=1
          iret=nf90_put_var(ncid,varid,real(outvar1(:,1:idim2p)),data_start_4d,data_count_4d)
          data_start_4d(1)=2
          iret=nf90_put_var(ncid,varid,real(outvar2(:,1:idim2p)),data_start_4d,data_count_4d)
        endif !ivs
        !write(12,*)'3D:',it_main,varid,var_nm2(1:len_var),iret2,NF90_NOERR
      endif !2/3D
 
      end subroutine writeout_nc

!===============================================================================
      subroutine fill_nc_header(iopen)
!-------------------------------------------------------------------------------
!     Create nc file and define dimension and static info for netcdf output
!     Input: iopen=1: close nc file handle. 0: do not close
!-------------------------------------------------------------------------------
      implicit none

      integer, intent(in) :: iopen
      character(len=140) :: fname
      character(len=4) :: fgb

      integer :: iret

      write(ifile_char,'(i12)') ifile !convert ifile to a string
      ifile_char=adjustl(ifile_char)  !place blanks at end
      ifile_len=len_trim(ifile_char)  !length without trailing blanks
      fgb='0000' 
      write(fgb,'(i4.4)') myrank
      fname=out_dir(1:len_out_dir)//('schout_'//fgb//'_'//ifile_char(1:ifile_len)//'.nc')
!'

      if(iopen==1) iret=nf90_close(ncid)
      iret=nf90_create(trim(adjustl(fname)),OR(NF90_NETCDF4,NF90_CLOBBER),ncid)
      iret=nf90_def_dim(ncid,'nSCHISM_hgrid_node',np,node_dim)
      iret=nf90_def_dim(ncid,'nSCHISM_hgrid_face',ne,nele_dim)
      iret=nf90_def_dim(ncid,'nSCHISM_hgrid_edge',ns,nedge_dim)
      iret=nf90_def_dim(ncid,'nMaxSCHISM_hgrid_face_nodes',4, four_dim)
      iret=nf90_def_dim(ncid,'nSCHISM_vgrid_layers',nvrt,nv_dim)
      iret=nf90_def_dim(ncid,'one',1,one_dim)
      iret=nf90_def_dim(ncid,'two',2,two_dim)
      iret=nf90_def_dim(ncid,'time', NF90_UNLIMITED,time_dim)

      time_dims(1)=time_dim
      iret=nf90_def_var(ncid,'time',NF90_FLOAT,time_dims,itime_id)
      if(iret.ne.NF90_NOERR) call parallel_abort('fill_nc_header: time dim')
!'
      iret=nf90_put_att(ncid,itime_id,'i23d',0) !set i23d flag
      iret=nf90_enddef(ncid)

      end subroutine fill_nc_header
!===============================================================================
!===============================================================================
! END FILE I/O module
!===============================================================================
!===============================================================================
    end module schism_io
