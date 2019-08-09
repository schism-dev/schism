!===============================================================================
!===============================================================================
! ELCIRC FILE I/O SUBROUTINES
!
! subroutine write_obe
! subroutine report_timers
! subroutine get_param

!===============================================================================
!===============================================================================

      subroutine write_obe
!-------------------------------------------------------------------------------
! Output centers.bp and sidecenters.bp
! NOTE: Valid for single processor only!
!-------------------------------------------------------------------------------
      use elfe_glbl
      use elfe_msgp
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
          write(32,*) i,real(xctr(i)),real(yctr(i)),real(dpe(i))
        else
          tmp1=sum(xlon(elnode(1:3,i)))/3
          tmp2=sum(ylat(elnode(1:3,i)))/3
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
!#ifdef USE_MPIMODULE
!      use mpi
!#endif
      use elfe_glbl, only : rkind,mxtimer,wtimer
      use elfe_msgp
      implicit none
!#ifndef USE_MPIMODULE
      include 'mpif.h'
!#endif
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
! Routine to read in param.in-like inputs
!===============================================================================
!===============================================================================
      subroutine get_param(fname,varname,vartype,ivarvalue,varvalue1,varvalue2)
! Get a parameter from param.in
! Inputs:
!        fname: the file name (e.g. 'param.in')
!        varname: parameter name (string no longer than 90)
!        vartype: parameter value type (0: 2-char string; 1: integer; 2: float)
! Outputs:
!        ivarvalue: integer output;
!        varvalue1: float output;
!        varvalue2: 2-char string output.
! Format rules for param.in:
! (1) Lines beginning with "!" are comments; blank lines are ignored;
! (2) one line for each parameter in the format: keywords= value;
!     keywords are case sensitive; spaces allowed between keywords and "=" and value;
!     comments starting with "!"  allowed after value;
! (3) value is an integer, double, or 2-char string; for double, any of the format is acceptable:
!     40 40. 4.e1
!     Use of decimal point in integers is OK but discouraged.
      use elfe_glbl, only : rkind,errmsg
      use elfe_msgp, only : parallel_abort,myrank
      implicit none

      character(*),intent(in) :: fname 
      character(*),intent(in) :: varname
      integer,intent(in) :: vartype
      integer,intent(out) :: ivarvalue
      real(rkind),intent(out) :: varvalue1
      character(len=2),intent(out) :: varvalue2

      character(len=90) :: line_str,str_tmp,str_tmp2
      integer :: lstr_tmp,lstr_tmp2,line,len_str,loc,loc2

      str_tmp2=adjustl(varname)
      lstr_tmp2=len_trim(str_tmp2)
!      print*, varname !,str_tmp2(1:lstr_tmp2)

      ! Scan param.in
      open(15,file=fname,status='old')
      !open(15,file='param.in',status='old')

      rewind(15)
      line=0
      do
        line=line+1
        read(15,'(a)',end=99)line_str
        line_str=adjustl(line_str) !place blanks at end
        len_str=len_trim(line_str)
        if(len_str==0.or.line_str(1:1)=='!') cycle

        loc=index(line_str,'=')
        loc2=index(line_str,'!')
        if(loc2/=0.and.loc2-1<loc+1) call parallel_abort('READ_PARAM: comments before =')
!'

        str_tmp=''
        str_tmp(1:loc-1)=line_str(1:loc-1) !keyword
        str_tmp=trim(str_tmp)
        lstr_tmp=len_trim(str_tmp)
    
        if(str_tmp(1:lstr_tmp)==str_tmp2(1:lstr_tmp2)) then
          if(loc2/=0) then
            str_tmp2=line_str(loc+1:loc2-1)
          else
            str_tmp2=line_str(loc+1:len_str)
          endif
          str_tmp2=adjustl(str_tmp2)
          str_tmp2=trim(str_tmp2)
          if(vartype==0) then !string
            varvalue2=str_tmp2(1:2)
#ifdef DEBUG
            if(myrank==0) write(86,*)varname,' = ',varvalue2
#endif
          else if(vartype==1) then !integer
            read(str_tmp2,*)ivarvalue
#ifdef DEBUG
            if(myrank==0) write(86,*)varname,' = ',ivarvalue
#endif
          else if(vartype==2) then !float
            read(str_tmp2,*)varvalue1
#ifdef DEBUG
            if(myrank==0) write(86,*)varname,' = ',real(varvalue1)
#endif
          else
            write(errmsg,*)'read_param: unknown type:',vartype
            call parallel_abort(errmsg)
          endif
          exit
        endif
      enddo !scan param.in
  
!     print*, 'Found it on line: ',line
      close(15)
      return

99    close(15)
      write(errmsg,*)'Failed to find parameter:',varname
      call parallel_abort(errmsg)

      end subroutine get_param

!===============================================================================
!===============================================================================
! END ELFE FILE I/O SUBROUTINES
!===============================================================================
!===============================================================================
