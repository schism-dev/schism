!> Command line parsing module
!! This module extends the CLA module by Ed Zaron
!! Method of use:
!! replace cla_fatal and cla_message
!! call cla_init to perform some initialization
!! call cla_register('-f','--foo', 'an integer', cla_int ,'72')
module argparse

#ifdef f2003
  use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit
#else
#define stdin  5
#define stdout 6
#define stderr 0
#endif
  
  implicit none

  ! Command Line Arguments
  ! A key-value parser for the command line

  ! Precision constants
  integer(kind=4), public, parameter :: int_kind = 4
  integer(kind=4), public, parameter :: real_kind = 4
  integer(kind=4), public, parameter :: ptr_kind = 8
  
  integer(kind=4), public, parameter :: STRLEN = 80
  integer(kind=4), public, parameter :: XSTRLEN = 256
  
  
  integer(kind=int_kind), parameter :: &
       cla_int   = 1, &
       cla_float = 2, &
       cla_char  = 3, &
       cla_xchar = 4, & ! NOT IMPLEMENTED
       cla_logical=5, &
       cla_flag  = 6

  character(len=STRLEN), dimension(6) :: cla_kindstr
  character(len=STRLEN), private :: cla_empty
  character(len=STRLEN), dimension(6) :: cla_true_str
 
  ! For key word arguments
  type, private :: cla_t
     character(len=2)  :: key
     character(len=STRLEN)  :: longkey 
     character(len=XSTRLEN) :: description
     integer(kind=int_kind) :: kind
     character(len=STRLEN)  :: default
  end type cla_t

  ! For positional arguments
  type, private :: posarg_t
     character(len=STRLEN)  :: key
     character(len=XSTRLEN) :: description
     integer(kind=int_kind) :: kind
     character(len=STRLEN)  :: default     
  end type posarg_t
  
  
  type(cla_t), private, dimension(:), pointer :: cla_registry
  type(posarg_t), private, dimension(:), pointer :: pcla_registry
  
  
  integer(kind=int_kind), private :: cla_num
  integer(kind=int_kind), private :: pcla_num  

  character(len=STRLEN) :: cmd_name        !< name of command for help
  character(len=2048) :: cmd_description !< description for help

  interface cla_get
     module procedure &
          cla_get_float_r4, &
          cla_get_float_r8, &
          cla_get_int_i4, &
          cla_get_int_i8, &
          cla_get_char, &
          cla_get_logical
  end interface

  contains

    subroutine cla_message(message)
    character(LEN=*) message
    ! Default message without print statements. 
    ! May need to be modified for, e.g. MPI codes 
    write(stdout,*)message
    end subroutine  
  
    subroutine cla_fatal(message)
    character(LEN=*) message
    ! Default stop with message without print statements. 
    ! May need to be modified for, e.g. MPI codes 
    write(stderr,*)message
    stop 6
    end subroutine
  
    subroutine cla_init(cmd,description)
      character(len=*) :: cmd
      character(len=*),optional :: description
      cmd_name = trim(cmd)
      if (present(description))then
          cmd_description = trim(description)
      else
          cmd_description = ''
      end if
      ! Allocate a zero size registry, just so that it gets
      ! associated.
      cla_num = 0
      allocate(cla_registry(0))
      allocate(pcla_registry(0))      
      cla_kindstr(cla_int)     = 'integer'
      cla_kindstr(cla_float)   = 'float'
      cla_kindstr(cla_char)    = 'character'
      cla_kindstr(cla_xchar)   = 'xcharacter' !NOT IMPLEMENTED
      cla_kindstr(cla_logical) = 'logical'
      cla_kindstr(cla_flag)    = 'flag'
      cla_empty='THIS_IS_THE_EMPTY_STRING'
      cla_true_str(1)='true'
      cla_true_str(2)='on'
      cla_true_str(3)='1'
      cla_true_str(4)='t'
      cla_true_str(5)='T'
      cla_true_str(6)='.true.'
    end subroutine cla_init

    subroutine pcla_register(key,description,kkind,default)
      character(len=*)  :: key
      character(len=*) :: description
      integer(kind=int_kind) :: kkind
      character(len=*)  :: default
      type(posarg_t), dimension(:), pointer :: pcla_registry_tmp
      integer(kind=int_kind) :: i


      ! This is a dumb way to increase the size of the
      ! registry of command line arguments, but there
      ! should not be so many arguments that either speed
      ! or memory is an issue.
      allocate(pcla_registry_tmp(pcla_num+1))
      do i=1,pcla_num
         pcla_registry_tmp(i)%key         = pcla_registry(i)%key
         pcla_registry_tmp(i)%description = pcla_registry(i)%description
         pcla_registry_tmp(i)%kind        = pcla_registry(i)%kind
         pcla_registry_tmp(i)%default     = pcla_registry(i)%default
         if (index(trim(key),' ') /= 0) then
            call cla_fatal('Error: pcla key contains a space character.')
         end if
         if (cla_str_eq(trim(pcla_registry(i)%key),trim(key))) then
            call cla_fatal('pcla key already registered: ' // &
                           pcla_registry(i)%key)
         end if
      end do
      pcla_num = pcla_num + 1
      deallocate(pcla_registry)
      allocate(pcla_registry(pcla_num))
      do i=1,pcla_num-1
         pcla_registry(i)%key         = pcla_registry_tmp(i)%key
         pcla_registry(i)%description = pcla_registry_tmp(i)%description
         pcla_registry(i)%kind        = pcla_registry_tmp(i)%kind
         pcla_registry(i)%default     = pcla_registry_tmp(i)%default
      end do
      i = pcla_num
      pcla_registry(i)%key         = key
      pcla_registry(i)%description = description
      pcla_registry(i)%description = description
      pcla_registry(i)%kind        = kkind
      pcla_registry(i)%default     = default
      deallocate(pcla_registry_tmp)
    end subroutine

    
    subroutine cla_register(key,longkey,description,kkind,default)
      character(len=2)  :: key
      character(len=*) :: longkey
      character(len=*) :: description
      integer(kind=int_kind) :: kkind
      character(len=*)  :: default
      type(cla_t), dimension(:), pointer :: cla_registry_tmp
      integer(kind=int_kind) :: i

      if (len(key) > 2)then
        call cla_fatal("The short key should be a dash plus one character (-e)")
      end if
      
      ! This is a dumb way to increase the size of the
      ! registry of command line arguments, but there
      ! should not be so many arguments that either speed
      ! or memory is an issue.
      allocate(cla_registry_tmp(cla_num+1))
      do i=1,cla_num
         cla_registry_tmp(i)%key         = cla_registry(i)%key
         cla_registry_tmp(i)%longkey         = cla_registry(i)%longkey
         cla_registry_tmp(i)%description = cla_registry(i)%description
         cla_registry_tmp(i)%kind        = cla_registry(i)%kind
         cla_registry_tmp(i)%default     = cla_registry(i)%default
         if (index(trim(key),' ') /= 0) then
            call cla_fatal('Attempt to register cla key containing a space.')
         end if
         if (index(trim(longkey),' ') /= 0) then
            call cla_fatal('Attempt to register long key containing space')
         end if         
         if (cla_str_eq(trim(cla_registry(i)%key),trim(key))) then
            call cla_fatal('cla key already registered:'// &
                           trim(key))
         end if
      end do
      cla_num = cla_num + 1
      deallocate(cla_registry)
      allocate(cla_registry(cla_num))
      do i=1,cla_num-1
         cla_registry(i)%key         = cla_registry_tmp(i)%key
         cla_registry(i)%longkey     = cla_registry_tmp(i)%longkey
         cla_registry(i)%description = cla_registry_tmp(i)%description
         cla_registry(i)%kind        = cla_registry_tmp(i)%kind
         cla_registry(i)%default     = cla_registry_tmp(i)%default
      end do
      i = cla_num
      cla_registry(i)%key         = key
      cla_registry(i)%longkey     = longkey
      cla_registry(i)%description = description
      cla_registry(i)%kind        = kkind
      cla_registry(i)%default     = default
      deallocate(cla_registry_tmp)
    end subroutine

    subroutine cla_show
      integer(kind=int_kind) :: i
      call cla_message('General usage:')
      call cla_message('  command -[key] [value] --[longkey] ARG1')
      call cla_message('  The key/value pairs must be matched if they appear.')
      call cla_message('  Key/value pairs and flags may be in any order.')
      call cla_message(' ')
      call cla_message('The following command line arguments and switches are expected:')
      do i=1,cla_num
         !call cla_message('---------- i: '// i)
         call cla_message('         key: '// trim(cla_registry(i)%key))
         call cla_message('     longkey: '// trim(cla_registry(i)%longkey))
         call cla_message(' description: '// trim(cla_registry(i)%description))
         call cla_message('        kind: '// &
                          trim(cla_kindstr(cla_registry(i)%kind)))
         call cla_message('     default: '// trim(cla_registry(i)%default))
      end do
      call cla_message('The following positional (non-keyword) arguments are expected')
      do i=1,pcla_num
         !call cla_message('---------- i: '// i)
         call cla_message('         key: '// trim(pcla_registry(i)%key))
         call cla_message(' description: '// trim(pcla_registry(i)%description))
         call cla_message('        kind: '// &
                          trim(cla_kindstr(pcla_registry(i)%kind)))
         call cla_message('     default: '// trim(pcla_registry(i)%default))
      end do
      
      call cla_message(' ')
      call cla_message('Also, -?, -h, -H, -help, --help, and --usage are recognized.')
      call cla_message(' ')
    end subroutine cla_show

    subroutine cla_help
      implicit none
      character(len=STRLEN) :: usage_line
      integer(kind=int_kind) :: i
      usage_line = '  '//trim(cmd_name)// ' [--key] [value] | [-flag] '
      do i = 1,pcla_num
         usage_line = trim(usage_line) // " " // pcla_registry(i)%key
      end do
      write(stdout,*)'Usage:'
      write(stdout,*)trim(usage_line)
      !write(stdout,*)'  ',trim(cmd_name),' [--key] [value] | -[flag] arg1 arg2 ...'
      write(stdout,*)' '
      if (len_trim(cmd_description) > 0)then
        write(stdout,*)trim(cmd_description)
        write(stdout,*)' '
      end if
      write(stdout,*)'Options and flags:'
      if (cla_num == 0) write(stdout,*)"None"
      do i=1,cla_num
         if (cla_registry(i)%kind == cla_flag) then
            write(stdout,'(1x,a,1x,a24,":",4x,a)')trim(cla_registry(i)%key), &
                                       trim(cla_registry(i)%longkey), &
                                       trim(cla_registry(i)%description)
         else
            write(stdout,'(1x,a,1x,a24,":",4x,a,2x,"{",a,"}")')trim(cla_registry(i)%key), &
                                 trim(cla_registry(i)%longkey), &   
                                 trim(cla_registry(i)%description), & 
                                 trim(cla_registry(i)%default) 
         endif
      end do
      write(stdout,*)' '

      write(stdout,*)'Positional arguments:'
      if (pcla_num == 0) write(stdout,*)"None"
      do i=1,pcla_num
        write(stdout,'(1x,a,":",1x,a,4x,a)')trim(pcla_registry(i)%key), &
                                trim(pcla_registry(i)%description)
      end do
      
      write(stdout,*)' '
      write(stdout,*)'Also, -?, -h, -H, -help, --help, and --usage are recognized.'
      write(stdout,*)' '
    end subroutine cla_help

    integer function cla_eq(str1,str2)
      implicit none
      character(*) :: str1, str2
      cla_eq = index(trim(str1),trim(str2))*index(trim(str2),trim(str1))
    end function cla_eq

    logical function key_arg_match(key,longkey,arg)
      implicit none
      ! do a match that includes two alternate keys and possibility of = in arg
      integer :: iequal
      character(*) :: key,longkey,arg
      key_arg_match = .false.
      key_arg_match = cla_str_eq(trim(key),trim(arg)) .or. &
                      cla_str_eq(trim(longkey),trim(arg))
      if (key_arg_match) return
      iequal = index(arg,"=")
      if (iequal > 1) &
         key_arg_match = cla_str_eq(trim(key),arg(1:(iequal-1))) .or. &
                         cla_str_eq(trim(longkey),arg(1:(iequal-1)))      
    end function key_arg_match   
    

    logical function cla_str_eq(str1,str2)
      implicit none
      character(*) :: str1, str2
      integer :: str_test
      str_test = index(trim(str1),trim(str2))*index(trim(str2),trim(str1))
      cla_str_eq = .false.
      if (str_test /= 0) cla_str_eq = .true.
    end function cla_str_eq

    subroutine cla_validate
      implicit none
      character(len=STRLEN) :: arg
      character(len=STRLEN)  :: value, key
      integer(kind=int_kind) :: ncla, k, kk, iequal
      ncla = command_argument_count()
      if (ncla == 0) return
      
      ! First check for -?, -h, -H, -help, or --help flags.
      call get_command_argument(1,arg)
      key = trim(arg)
      if (cla_str_eq(trim(key),'-h')      .or. &
          cla_str_eq(trim(key),'-?')      .or. &
          cla_str_eq(trim(key),'/?')      .or. &
          cla_str_eq(trim(key),'-H')      .or. &
          cla_str_eq(trim(key),'-help')   .or. &
          cla_str_eq(trim(key),'--help')  .or. &
          cla_str_eq(trim(key),'--usage')      &
          ) then
         call cla_help
         stop
      endif
    end subroutine cla_validate
    
    logical function cla_key_present(key)
      implicit none
      character(len=STRLEN) :: arg
      character(len=*)  :: key
      character(len=STRLEN) :: longkey
      character(len=2) :: shortkey
      character(len=STRLEN)  :: value
      
      integer(kind=int_kind) :: ncla, k, kk
!      integer :: command_argument_count
!      external command_argument_count

      !     Loop over the command line arguments to assign to
      !     value.
      !     Note that no error is reported if the key was NOT
      !     registered, but it is present on the command line.
 
      cla_key_present = .false.
      

!      print *,'Calling cla_key_present with key = ',trim(key)
      value = trim(cla_empty)
      do kk=1,cla_num
         ! must test for exact match, not just substring
         if (key_arg_match(cla_registry(kk)%key, &
                           cla_registry(kk)%longkey, &
                           key))then
            value = trim(cla_registry(kk)%default)
            longkey = cla_registry(kk)%longkey
            shortkey = cla_registry(kk)%key
            exit
         end if
      end do
      
      if (index(trim(value),trim(cla_empty)) /= 0) then
         call cla_help
         call cla_fatal('Unknown command line argument: '//trim(key))
      endif
      
      ncla = command_argument_count()
      if (ncla == 0) return
      
      do k=1,ncla
         call get_command_argument(k,arg)
         ! test for exact match
         if (key_arg_match(shortkey,longkey,arg))then
            cla_key_present = .true.
            return
         endif
      enddo
      
    end function cla_key_present

    subroutine cla_get_char(key,value)
      implicit none
      character(len=STRLEN) :: arg
      character(len=*)  :: key
      character(len=2)  :: shortkey
      character(len=STRLEN)  :: longkey
      character(len=STRLEN)  :: pvalue
      character(len=*)  :: value
      character(len=STRLEN) :: kkey
      integer(kind=int_kind) :: ncla, k, kkind, iequal
      integer :: kk, kmatch, ordinal
      logical :: just_matched, prev_matched, is_match, carryover
!      integer :: command_argument_count
!      external command_argument_count

      !     Loop over the command line arguments to assign to
      !     value.
      !     Note that no error is reported if the key was NOT
      !     registered, but it is present on the command line.
      if (index(key,"-") /= 1)then
         ! assume positional argument, confirm by matching name
         ordinal = -1
         do k=1,pcla_num
             ! must test for exact match, not just substring
             if (cla_str_eq(trim(pcla_registry(k)%key),trim(key))) then
                pvalue = trim(pcla_registry(k)%default)
                ordinal = k
             end if
          end do
      
          if (index(trim(value),trim(cla_empty)) /= 0) then
             print *,'Unknown command line argument: ',trim(key)
             call cla_help
             stop 5
          endif
          
          if (ordinal > 0) then
              ncla = command_argument_count()
              if (ncla == 0) return
              kmatch = 0
              prev_matched = .False.
              
              do k=1,ncla
                 call get_command_argument(k,arg)
                 ! test for exact match among key args
                 just_matched = .False.    
                 do kk = 1, cla_num
                     kkey = cla_registry(kk)%key
                     kkind = cla_registry(kk)%kind
                     is_match = key_arg_match(kkey,cla_registry(kk)%longkey,arg)
                     carryover = kkind/=cla_flag
                     iequal = index(arg,"=")
                     if (is_match .and. iequal > 1)then
                        carryover = .False.
                     end if
                     just_matched = just_matched .or. is_match
                     if (just_matched) exit  ! preserve kkey and kkind
                 end do
                 if (just_matched .or. prev_matched )then
                     ! current arg part of keyword arg constructal
                     ! so this is not positional
                     prev_matched = just_matched .and. carryover
                     carryover = .False.
                     cycle
                 end if
                 kmatch = kmatch + 1
                     ! increment # of positionals
                 if(kmatch == ordinal) then
                     pvalue = trim(arg)
                     value=pvalue
                     return
                 end if

              end do
          end if
          value = pvalue
          return
      end if
      
      
      ! keyword
      do k=1,cla_num
         ! must test for exact match, not just substring
         if (key_arg_match(cla_registry(k)%key, &
                           cla_registry(k)%longkey, &
                           key)) then
            shortkey = cla_registry(k)%key
            longkey = cla_registry(k)%longkey
            
            value = trim(cla_registry(k)%default)
            kkind = cla_registry(k)%kind
         end if
      end do
      
      if (index(trim(value),trim(cla_empty)) /= 0) then
         print *,'Unknown command line argument: ',trim(key)
         call cla_help
         stop 5
      endif
      
      ncla = command_argument_count()
      if (ncla == 0) return
      
      do k=1,ncla
         call get_command_argument(k,arg)
         ! test for exact match
         if (key_arg_match(shortkey,longkey,trim(arg))) then
            if (kkind == cla_flag) then
               value = 't'
               return
            else
               iequal = index(arg,"=")
               if (iequal < 1)then
                 call get_command_argument(k+1,arg)       
                 value = trim(arg)
                 return
               else
                 value=arg(iequal+1:len_trim(arg))
                 return
               end if                             
            endif
        end if
      enddo

    end subroutine cla_get_char
    

  subroutine cla_get_float_r4(key,float_value)
    implicit none
    character(len=*)       :: key
    character(len=STRLEN)  :: value
    real(kind=4)           :: float_value

    
    call cla_get_char(key,value)
    if (index(trim(value),trim(cla_empty)) == 0) read(value,*)float_value
  end subroutine cla_get_float_r4

  subroutine cla_get_float_r8(key,float_value)
    implicit none
    character(len=*)       :: key
    character(len=STRLEN)  :: value
    real(kind=8)           :: float_value

    call cla_get_char(key,value)
    if (index(trim(value),trim(cla_empty)) == 0) &
       read(value,*,ERR=100)float_value
    return
 100  call cla_fatal("Input value not correct type: "//value)
  end subroutine cla_get_float_r8


  subroutine cla_get_int_i4(key,int_value)
    implicit none
    character(len=*)       :: key
    character(len=STRLEN)  :: value
    integer(kind=4)        :: int_value
    
    call cla_get_char(key,value)
    if (index(trim(value),trim(cla_empty)) == 0) read(value,*,ERR=100)int_value
    return
 100  call cla_fatal("Input value not correct type: "//value)    
    
  end subroutine cla_get_int_i4

  subroutine cla_get_int_i8(key,int_value)
    implicit none
    character(len=*)       :: key
    character(len=STRLEN)  :: value
    integer(kind=8)        :: int_value
    
    call cla_get_char(key,value)
    if (index(trim(value),trim(cla_empty)) == 0) read(value,*,ERR=101)int_value
    return
 101  call cla_fatal("Input value not correct type: "//value)        
  end subroutine cla_get_int_i8

  subroutine cla_get_logical(key,logical_value)
    implicit none
    character(len=*)  :: key
    character(len=STRLEN)  :: value
    logical :: logical_value
    integer(kind=int_kind) :: k
    
    logical_value = .false.

    call cla_get_char(key,value)
    if (index(trim(value),trim(cla_empty)) == 0) then
       do k=1,6
          if (index(trim(value),trim(cla_true_str(k))) /= 0) then
             logical_value = .true.
          endif
       end do
    end if
  end subroutine cla_get_logical

  subroutine cla_get_flag(key,logical_value)
    implicit none
    character(len=*)  :: key
    character(len=STRLEN)  :: value
    logical :: logical_value
    integer(kind=int_kind) :: k
    
    logical_value = .false.

    call cla_get_char(key,value)
    if (index(trim(value),trim(cla_empty)) == 0) then
       do k=1,6
          if (index(trim(value),trim(cla_true_str(k))) /= 0) then
             logical_value = .true.
          endif
       end do
    end if
  end subroutine cla_get_flag
  

end module argparse
