
module icepack_warnings

      use icepack_kinds

      implicit none

      private

      ! warning messages
      character(len=char_len_long), dimension(:), allocatable :: warnings
      integer :: nWarnings = 0
      integer, parameter :: nWarningsBuffer = 10 ! incremental number of messages

      ! abort flag, accessed via icepack_warnings_setabort and icepack_warnings_aborted
      logical :: warning_abort = .false.

      ! public string for all subroutines to use
      character(len=char_len_long), public :: warnstr

      public :: &
        icepack_warnings_clear,    &
        icepack_warnings_print,    &
        icepack_warnings_flush,    &
        icepack_warnings_aborted,  &
        icepack_warnings_add,      &
        icepack_warnings_setabort

      private :: &
        icepack_warnings_getall,   &
        icepack_warnings_getone

!=======================================================================

contains

!=======================================================================
!autodocument_start icepack_warnings_aborted
! turn on the abort flag in the icepack warnings package
! pass in an optional error message

      logical function icepack_warnings_aborted(instring)

        character(len=*),intent(in), optional :: instring

!autodocument_end

        character(len=*),parameter :: subname='(icepack_warnings_aborted)'

        icepack_warnings_aborted = warning_abort
        if (warning_abort .and. present(instring)) then
           call icepack_warnings_add(subname//' ... '//trim(instring))
        endif

      end function icepack_warnings_aborted

!=======================================================================

      subroutine icepack_warnings_setabort(abortflag,file,line)

        logical, intent(in) :: abortflag
        character(len=*), intent(in), optional :: file
        integer, intent(in), optional :: line

        character(len=*),parameter :: subname='(icepack_warnings_setabort)'

        ! try to capture just the first setabort call

        if (abortflag) then
          write(warnstr,*) subname,abortflag
          if (present(file)) write(warnstr,*) trim(warnstr)//' :file '//trim(file)
          if (present(line)) write(warnstr,*) trim(warnstr)//' :line ',line
          call icepack_warnings_add(warnstr)
        endif

        warning_abort = abortflag

      end subroutine icepack_warnings_setabort

!=======================================================================
!autodocument_start icepack_warnings_clear
! clear all warning messages from the icepack warning buffer

      subroutine icepack_warnings_clear()

!autodocument_end

        character(len=*),parameter :: subname='(icepack_warnings_clear)'

        nWarnings = 0

      end subroutine icepack_warnings_clear

!=======================================================================
      
      subroutine icepack_warnings_getall(warningsOut)

        character(len=char_len_long), dimension(:), allocatable, intent(out) :: &
             warningsOut
 
        integer :: iWarning
        character(len=*),parameter :: subname='(icepack_warnings_getall)'

        if (allocated(warningsOut)) deallocate(warningsOut)
        allocate(warningsOut(nWarnings))

        do iWarning = 1, nWarnings
           warningsOut(iWarning) = trim(icepack_warnings_getone(iWarning))
        enddo

      end subroutine icepack_warnings_getall

!=======================================================================
!autodocument_start icepack_warnings_print
! print all warning messages from the icepack warning buffer

      subroutine icepack_warnings_print(iounit)

        integer, intent(in) :: iounit

!autodocument_end

        integer :: iWarning
        character(len=*),parameter :: subname='(icepack_warnings_print)'

! tcraig
! this code intermittenly aborts on recursive IO errors with intel
! not sure if it's OMP or something else causing this
!$OMP MASTER
        do iWarning = 1, nWarnings
          write(iounit,*) trim(icepack_warnings_getone(iWarning))
        enddo
!$OMP END MASTER

      end subroutine icepack_warnings_print

!=======================================================================
!autodocument_start icepack_warnings_flush
! print and clear all warning messages from the icepack warning buffer

      subroutine icepack_warnings_flush(iounit)

        integer, intent(in) :: iounit

!autodocument_end

        character(len=*),parameter :: subname='(icepack_warnings_flush)'

        if (nWarnings > 0) then
          call icepack_warnings_print(iounit)
        endif
        call icepack_warnings_clear()

      end subroutine icepack_warnings_flush

!=======================================================================

      subroutine icepack_warnings_add(warning)

        character(len=*), intent(in) :: warning ! warning to add to array of warnings

        ! local 

        character(len=char_len_long), dimension(:), allocatable :: warningsTmp
        integer :: &
             nWarningsArray, & ! size of warnings array at start
             iWarning ! warning index
        character(len=*),parameter :: subname='(icepack_warnings_add)'

!$OMP CRITICAL (omp_warnings_add)
        ! check if warnings array is not allocated
        if (.not. allocated(warnings)) then

           ! allocate warning array with number of buffer elements
           allocate(warnings(nWarningsBuffer))

           ! set initial number of nWarnings
           nWarnings = 0

        ! already allocated
        else

           ! find the size of the warnings array at the start
           nWarningsArray = size(warnings)
           
           ! check to see if need more space in warnings array
           if (nWarnings + 1 > nWarningsArray) then
           
              ! allocate the temporary warning storage
              allocate(warningsTmp(nWarningsArray))

              ! copy the warnings to temporary storage
              do iWarning = 1, nWarningsArray
                 warningsTmp(iWarning) = trim(warnings(iWarning))
              enddo ! iWarning

              ! increase the size of the warning array by the buffer size
              deallocate(warnings)
              allocate(warnings(nWarningsArray + nWarningsBuffer))

              ! copy back the temporary stored warnings
              do iWarning = 1, nWarningsArray
                 warnings(iWarning) = trim(warningsTmp(iWarning))
              enddo ! iWarning

              ! deallocate the temporary storage
              deallocate(warningsTmp)

           endif
              
        endif

        ! increase warning number
        nWarnings = nWarnings + 1
!$OMP END CRITICAL (omp_warnings_add)

        ! add the new warning
        warnings(nWarnings) = trim(warning)

      end subroutine icepack_warnings_add

!=======================================================================

      function icepack_warnings_getone(iWarning) result(warning)

        integer, intent(in) :: iWarning

        character(len=char_len_long) :: warning

        character(len=*),parameter :: subname='(icepack_warnings_getone)'

        if (iWarning <= nWarnings) then
           warning = warnings(iWarning)
        else
           warning = ""
        endif

      end function icepack_warnings_getone

!=======================================================================

end module icepack_warnings
