module gen_modules_clock
  !combining RT and Lars version
  !Used in USE_MICE
  use schism_glbl
  use schism_msgp, only : myrank
  implicit none
  save
  real(rkind)              :: timeold, timenew     !time in a day, unit: sec
  integer                  :: dayold, daynew       !day in a year
  integer                  :: yearold, yearnew     !year before and after time step
  integer                  :: month_mice, day_in_month  !month and day in a month
  integer                  :: fleapyear            !1 fleapyear, 0 not 
  integer                  :: ndpyr                !number of days in yearnew 
  integer                  :: num_day_in_month(0:1,12)
  character(4)             :: cyearold, cyearnew   !year as character string      
  data num_day_in_month(0,:) /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
  data num_day_in_month(1,:) /31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/


contains
  !
  !--------------------------------------------------------------------------------
  !
  subroutine clock

    implicit none
    integer          :: i
    real(rkind)    :: aux1, aux2
    !
    timeold=timenew 
    dayold=daynew
    yearold=yearnew

    ! update time
    timenew=timenew+dt          
     
    ! update day
    if (timenew>86400.d0) then  !assumed that time step is less than one day!
       daynew=daynew+1
       timenew=timenew-86400.d0
    endif

    ! update year
    if (daynew>ndpyr) then
       daynew=1
       yearnew=yearnew+1
       call check_fleapyr(yearnew, fleapyear)
       ndpyr=365+fleapyear
       write(cyearold,'(i4)') yearold
       write(cyearnew,'(i4)') yearnew
    endif

    ! find month and dayinmonth at new time step
    aux1=0
    do i=1,12
       aux2=aux1+num_day_in_month(fleapyear,i)
       if(daynew>aux1 .and. daynew<=aux2) then
          month_mice=i
          day_in_month=daynew-aux1
          exit
       end if
       aux1=aux2
    end do
       
  end subroutine clock
  !
  !--------------------------------------------------------------------------------
  !
  subroutine clock_init(time)
    
    !use g_parsup
    !use g_config
    implicit none

    real(rkind), intent(in) :: time
    integer          :: i,j, daystart, yearstart
    real(rkind)    :: aux1, aux2, timestart
 
    ! the model initialized at
    timenew=start_hour*3600
    yearnew=start_year
    call check_fleapyr(yearnew, fleapyear)
    ndpyr=365+fleapyear
    aux1=0
      do i=1,start_month
         aux1=aux1+num_day_in_month(fleapyear,i)
      end do
    aux1=aux1+start_day-num_day_in_month(fleapyear,start_month)
    daynew=aux1
   
   j=time/dt

   do i=1,j
      call clock_newyear 
      call clock
   enddo

    timestart=timenew
    daystart=daynew
    yearstart=yearnew

    ! init clock for this run
    !open(99,file=trim(ResultPath)//trim(runid)//'.clock',status='old')
    !read(99,*) timeold, dayold, yearold
    !read(99,*) timenew, daynew, yearnew
    !close(99)
    !if(daynew==0) daynew=1
    
    ! check if this is a restart or not
    !if(yearnew==yearstart .and. daynew==daystart .and. timenew==timestart) then
    !  r_restart=.false.
    !  yearold=yearnew-1 !required for checking if create new output files
    !else
    !   r_restart=.true.
    !end if

    ! year as character string 
    write(cyearold,'(i4)') yearold
    write(cyearnew,'(i4)') yearnew

    ! if restart model at beginning of a day, set timenew to be zero
    if (timenew==86400.d0) then  
       timenew=0.0
       daynew=daynew+1
    endif

    ! set timeold to be timenew, ready for initializing forcing fields,
    ! yearold should not be updated here, which is requird to open input files.
    ! timeold=timenew 
    ! dayold=daynew
    
    ! check fleap year
    
    !ndpyr=365+fleapyear

    ! find month and dayinmonth at the new time step
    aux1=0
    do i=1,12
       aux2=aux1+num_day_in_month(fleapyear,i)
       if(daynew>aux1 .and. daynew<=aux2) then
          month_mice=i
          day_in_month=daynew-aux1
          exit
       end if
       aux1=aux2
    end do

    !if(myrank==0) then
        !if(r_restart) then
           ! write(*,*)
          !  print *, achar(27)//'[31m'    //'____________________________________________________________'//achar(27)//'[0m'
         !   print *, achar(27)//'[5;7;31m'//' --> THIS IS A RESTART RUN !!!                              '//achar(27)//'[0m'
        !    write(*,"(A, F5.2, I4, I5)") '     > clock restarted at time:', timenew, daynew, yearnew
       !     write(*,*)
      !  else
     !       write(*,*)
    !        print *, achar(27)//'[32m'  //'____________________________________________________________'//achar(27)//'[0m'
    !        print *, achar(27)//'[7;32m'//' --> THIS IS A INITIALISATION RUN !!!                       '//achar(27)//'[0m'
    !        write(*,"(A, F5.2, I4, I5)")'     > clock initialized at time:', timenew, daynew, yearnew
    !        write(*,*)
    !    end if
    !end if
  
  end subroutine clock_init
  !
  !-------------------------------------------------------------------------------
  !
  subroutine clock_finish
    implicit none
    !
    real(rkind)            :: dum_timenew     !time in a day, unit: sec
    integer                  :: dum_daynew       !day in a year
    integer                  :: dum_yearnew     !year before and after time step
    
    dum_timenew = timenew
    dum_daynew  = daynew
    dum_yearnew = yearnew
    if ((dum_daynew==ndpyr) .and. (dum_timenew==86400.d0)) then
       dum_timenew=0.0
       dum_daynew=1
       dum_yearnew=yearold+1
    endif

    !open(99,file=trim(ResultPath)//trim(runid)//'.clock',status='unknown')
    !write(99,*) timeold, dayold, yearold
    !write(99,*) dum_timenew, dum_daynew, dum_yearnew
    !close(99)
  end subroutine clock_finish
  !
  !----------------------------------------------------------------------------
  !
  subroutine clock_newyear
    implicit none
    !
    if ((daynew>=ndpyr).and.(timenew==86400.d0)) then
       timenew=0.0
       daynew=1
       yearnew=yearold+1
       write(cyearnew,'(i4)') yearnew
    endif
  end subroutine clock_newyear
  !
  !----------------------------------------------------------------------------
  !
  subroutine check_fleapyr(year, flag)
    implicit none
    integer, intent(in) :: year      
    integer, intent(out):: flag

    flag=0

    !if(.not.include_fleapyear) return

    if ((mod(year,4)==0.and.mod(year,100)/=0) .or. mod(year,400)==0) then
       flag=1
    endif
  end subroutine check_fleapyr
  !
  !----------------------------------------------------------------------------
  !
end module gen_modules_clock
