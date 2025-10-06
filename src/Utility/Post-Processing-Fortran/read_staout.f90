! Read in 3D station outputs and output profile (at wet instances) 
! Inputs: screen, vgrid.in,station.in, outputs/staout_*
! Outputs: staout_*_prof.dat
! ifx -O2 -o read_staout read_staout.f90

  implicit real*8(a-h,o-z)
  character(len=10) :: mychar
  real*8, allocatable :: wild(:,:),zcor(:,:)

  print*, 'Input the variable ID in station.in (e.g, 5 for T):'
  read*, idvar
  if(idvar<=4) stop '2D output has no profile'
  write(mychar,'(i10)')idvar

  open(8,file='vgrid.in',status='old')
  read(8,*)ivcor
  read(8,*)nvrt
  close(8)

  open(9,file='station.in',status='old')
  open(10,file='outputs/staout_'//trim(adjustl(mychar)),status='old')
  open(12,file='staout_'//trim(adjustl(mychar))//'_prof.dat',status='replace')
  read(9,*); read(9,*)nsta
  allocate(wild(nvrt,nsta),zcor(nvrt,nsta))

  nt0=0
  do
    nt0=nt0+1

    !Odd line is station output, even line is 3D profile
    if(mod(nt0,2)==1) then
      read(10,*,end=99,err=99) !time,wild(1,:)
    else !profile
      read(10,*,end=99,err=99)time,wild(:,:),zcor(:,:)

      !Process dry etc. Raw output uses -1.e7 for dry or below bottom; +1.e7 if the station is outside mesh
      write(12,*)'Time (days)= ',time/86400.
      do i=1,nsta
        write(12,*)'Profile @station #= ',i
        do k=1,nvrt
          if(abs(wild(k,i))<5.d6) then
            if(abs(wild(k,i))>1.d2) print*, 'Abnormal values:',zcor(k,i),wild(k,i)
            write(12,*)k,real(zcor(k,i)),real(wild(k,i))
          endif !abs()
        enddo !k
      enddo !i
    endif !mod
  enddo
99 continue
  nt0=nt0-1
  print*, nt0,' lines read;'

  end
