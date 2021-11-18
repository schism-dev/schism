!----------------------------------------------------------------
!               M O D U L E   C S V _ U T I L I T I E S
!----------------------------------------------------------------
!> @file csv_utilities.F90
!>
!> @brief
!>   Utility routines.
!>
!> @details
!>   
!>
!> @author Jacob Williams
!> @copyright License BSD
!----------------------------------------------------------------

    module csv_utilities

    USE PaHM_Sizes, ONLY : WP, IP
    use csv_parameters

    private

    integer,parameter :: max_size_for_insertion_sort = 20 !> max size for using insertion sort.

    public :: unique
    public :: expand_vector
    public :: sort_ascending

    contains
!*******************************************************************************

  !----------------------------------------------------------------
  ! S U B R O U T I N E   E X P A N D _ V E C T O R
  !----------------------------------------------------------------
  !> @brief
  !>   Add elements to the integer vector in chunks.
  !>
  !> @details
  !>   
  !>
  !> @param
  !>   vec          The input integer vector (input/output)
  !> @param
  !>   n            Counter for last element added to `vec`; must be initialized to `size(vec)` \n
  !>                (or 0 if not allocated) before first call  (input/output)
  !> @param
  !>   chunk_size   Allocate `vec` in blocks of this size (>0)
  !> @param
  !>   val          The value to add to `vec` (optional)
  !> @param
  !>   finished     Set to true to return `vec` as its correct size (`n`) (optional)
  !>
  !----------------------------------------------------------------
    pure subroutine expand_vector(vec, n, chunk_size, val, finished)

    implicit none

    integer,dimension(:),allocatable,intent(inout) :: vec
    integer,intent(inout)       :: n           ! counter for last element added to `vec`.
                                               ! must be initialized to `size(vec)`
                                               ! (or 0 if not allocated) before first call
    integer,intent(in)          :: chunk_size  ! allocate `vec` in blocks of this size (>0)
    integer,intent(in),optional :: val         ! the value to add to `vec`
    logical,intent(in),optional :: finished    ! set to true to return `vec`
                                               ! as its correct size (`n`)

    integer,dimension(:),allocatable :: tmp  ! temporary array

    if (present(val)) then
        if (allocated(vec)) then
            if (n==size(vec)) then
                ! have to add another chunk:
                allocate(tmp(size(vec)+chunk_size))
                tmp(1:size(vec)) = vec
                call move_alloc(tmp,vec)
            end if
            n = n + 1
        else
            ! the first element:
            allocate(vec(chunk_size))
            n = 1
        end if
        vec(n) = val
    end if

    if (present(finished)) then
        if (finished) then
            ! set vec to actual size (n):
            if (allocated(tmp)) deallocate(tmp)
            allocate(tmp(n))
            tmp = vec(1:n)
            call move_alloc(tmp,vec)
        end if
    end if

    end subroutine expand_vector
!*******************************************************************************

!*******************************************************************************
!
!> Returns only the unique elements of the vector.

    function unique(vec,chunk_size) result(ivec_unique)

    implicit none

    integer,dimension(:),intent(in)    :: vec        !> a vector of integers
    integer,intent(in)                 :: chunk_size  !> chunk size for adding to arrays
    integer,dimension(:),allocatable   :: ivec_unique !> unique elements of `ivec`

    integer,dimension(size(vec)) :: ivec !> temp copy of vec
    integer :: i !> counter
    integer :: n !> number of unique elements

    ! first we sort it:
    ivec = vec ! make a copy
    call sort_ascending(ivec)

    ! add the first element:
    n = 1
    ivec_unique = [ivec(1)]

    ! walk through array and get the unique ones:
    if (size(ivec)>1) then
        do i = 2, size(ivec)
            if (ivec(i)/=ivec(i-1)) then
                call expand_vector(ivec_unique,n,chunk_size,val=ivec(i))
            end if
        end do
        call expand_vector(ivec_unique,n,chunk_size,finished=.true.)
    end if

    end function unique
!*******************************************************************************

!*******************************************************************************
!>
!> Sorts an integer array `ivec` in increasing order.
!> Uses a basic recursive quicksort
!> (with insertion sort for partitions with \f$ \le \f$ 20 elements).

    subroutine sort_ascending(ivec)

    implicit none

    integer,dimension(:),intent(inout) :: ivec

    call quicksort(1,size(ivec))

    contains

        recursive subroutine quicksort(ilow,ihigh)

        !> Sort the array

        implicit none

        integer,intent(in) :: ilow
        integer,intent(in) :: ihigh

        integer :: ipivot !> pivot element
        integer :: i      !> counter
        integer :: j      !> counter

        if ( ihigh-ilow<=max_size_for_insertion_sort .and. ihigh>ilow ) then

            ! do insertion sort:
            do i = ilow + 1,ihigh
                do j = i,ilow + 1,-1
                    if ( ivec(j) < ivec(j-1) ) then
                        call swap(ivec(j),ivec(j-1))
                    else
                        exit
                    end if
                end do
            end do

        else if ( ihigh-ilow>max_size_for_insertion_sort ) then

            ! do the normal quicksort:
            call partition(ilow,ihigh,ipivot)
            call quicksort(ilow,ipivot - 1)
            call quicksort(ipivot + 1,ihigh)

        end if

        end subroutine quicksort

        subroutine partition(ilow,ihigh,ipivot)

        !> Partition the array, based on the
        !> lexical ivecing comparison.

        implicit none

        integer,intent(in)  :: ilow
        integer,intent(in)  :: ihigh
        integer,intent(out) :: ipivot

        integer :: i,ip

        call swap(ivec(ilow),ivec((ilow+ihigh)/2))
        ip = ilow
        do i = ilow + 1, ihigh
            if ( ivec(i) < ivec(ilow) ) then
                ip = ip + 1
                call swap(ivec(ip),ivec(i))
            end if
        end do
        call swap(ivec(ilow),ivec(ip))
        ipivot = ip

        end subroutine partition

    end subroutine sort_ascending
!*******************************************************************************

!*******************************************************************************
!>
!> Swap two integer values.

    pure elemental subroutine swap(i1,i2)

    implicit none

    integer,intent(inout) :: i1
    integer,intent(inout) :: i2

    integer :: tmp

    tmp = i1
    i1  = i2
    i2  = tmp

    end subroutine swap
!*******************************************************************************

!*******************************************************************************
    end module csv_utilities
!*******************************************************************************
