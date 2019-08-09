!-----------------------------------------------------------------------
! cjd_bin - calculates the Chronological Julian Date
!
! Mike Zulauf, January 2007

      program cjd

! modules, implicit none, etc
#ifdef NEED_F2KCLI
      use f2kcli  ! if compiler doesn't implement f2003 CLI capability
#endif
      implicit none

! variables needed for command line options
      integer, parameter :: max_args = 20
      character(len=40)  :: arg(max_args)
      integer            :: num_args, i_arg
      character(len=2)   :: reverse_option = '-r'

! time variables
      integer yyyy, mm, dd, jd, jd_in
      logical reverse
      character(len=2) :: m2, d2

! determine how many command-line arguments there are
      num_args = command_argument_count()

! store command line arguments
      if (num_args .ge. 1) then
        do i_arg = 1, num_args
          call get_command_argument(i_arg, arg(i_arg))
        end do
      endif

! determine if this is to be a "forward" or "reverse" operation
      reverse = (trim(arg(1)) .eq. reverse_option)

      if ( (.not. reverse .and. (num_args .lt. 3)) .or.
     &     (reverse       .and. (num_args .lt. 2)) ) then
        write(*,*)
        write(*,*) 'you need to supply: YYYY MM DD'
        write(*,*) 'or reverse option:  -r NNNNNNN'
        stop
      endif

! forward
      if (.not. reverse) then

! extract inputs from command line arguments
        read(arg(1),*) yyyy
        read(arg(2),*) mm
        read(arg(3),*) dd

! output the Chronological Julian Date
50      format(i7)
        write(*,50) jd(yyyy, mm, dd)

! reverse
      else

! extract inputs from command line arguments
        read(arg(2),*) jd_in

! get calendar date
        call cdate(jd_in, yyyy, mm, dd)

! convert month and day to 2 character format
10      format ('0', i1)
20      format (i2)

        if (mm .le. 9) then
          write(m2,10) mm
        else
          write(m2,20) mm
        endif

        if (dd .le. 9) then
          write(d2,10) dd
        else
          write(d2,20) dd
        endif

! output the year month date
60      format(i4, 1x, a2, 1x, a2)
        write(*,60) yyyy, m2, d2
      
      endif

      end
!-----------------------------------------------------------------------
      FUNCTION JD(YYYY,MM,DD)
      INTEGER YYYY,MM,DD
!              DATE ROUTINE JD(YYYY,MM,DD) CONVERTS CALENDER DATE TO
!              JULIAN DATE.  SEE CACM 1968 11(10):657, LETTER TO THE
!              EDITOR BY HENRY F. FLIEGEL AND THOMAS C. VAN FLANDERN.
!    EXAMPLE JD(1970,1,1)=2440588
      JD=DD-32075+1461*(YYYY+4800+(MM-14)/12)/4
     &         +367*(MM-2-((MM-14)/12)*12)/12-3*
     &         ((YYYY+4900+(MM-14)/12)/100)/4
      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE CDATE(JD,YYYY,MM,DD)
!=======GIVEN A JULIAN DAY NUMBER, NNNNNNNN, YYYY,MM,DD ARE RETURNED AS
!              AS THE CALENDAR DATE. JD=NNNNNNNN IS THE JULIAN DATE
!              FROM AN EPOCK IN THE VERY DISTANT PAST.  SEE CACM
!              1968 11(10):657, LETTER TO THE EDITOR BY FLIEGEL AND
!              VAN FLANDERN.
!    EXAMPLE CALL CDATE(2440588,YYYY,MM,DD) RETURNS 1970 1 1 .
!
      INTEGER JD,YYYY,MM,DD,L,N
      L=JD+68569
      N=4*L/146097
      L=L-(146097*N + 3)/4
      YYYY=4000*(L+1)/1461001
      L=L-1461*YYYY/4+31
      MM=80*L/2447
      DD=L-2447*MM/80
      L=MM/11
      MM=MM + 2 - 12*L
      YYYY=100*(N-49) + YYYY + L
      RETURN
      END
!-----------------------------------------------------------------------
