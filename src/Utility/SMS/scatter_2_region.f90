! Convert SMS scatter formar .xy to .rgn (gredit), assuming it's single
! connected and in proper order
! To generate .xy in SMS: create a map and polygon (mesh generator), and edit, clean/build. Then convert 
!  map to Scatter (Arc end points and vertices), and highlight Scatter set and save as 
!  'Scatter pt .xy'.

! Inputs: scatter.xy (SMS)
! Outputs: scatter.rgn

! ifort -O2 -mcmodel=medium -CB -Bstatic -o scatter_2_region scatter_2_region.f90

  implicit real*8(a-h,o-z)
  character(len=100) :: str1,str2
  real*8, allocatable :: xy(:,:)

  open(10,file='scatter.xy',status='old')
  open(12,file='scatter.rgn',status='replace')
  line=0
  line_start=-1
  do 
    read(10,*,err=99)str1
    line=line+1
    str2=adjustl(str1)
!    len1=len_trim(str2)
    if(str2(1:3).eq."IXY") then
      line_start=line
      exit
    endif
  enddo 

99 if(line_start==-1) stop 'did not find starting line'
  print*, 'Starting line found:',line_start

  write(12,*)'Region written by ACE/gredit'
  write(12,*)1

  rewind(10)
  do i=1,line_start-1; read(10,*); enddo;
  read(10,*)str1,npts 
  write(12,*)npts,1
  allocate(xy(2,npts))
  print*, '# of points in scatter=',npts
  do i=1,npts
    read(10,*)j,xy(1:2,i)
    write(12,*)xy(1:2,i)
  enddo !i
 
  !Repeat 1st pt to close region: no need to do this as .rgn can be either way
!  write(12,*)xy(1:2,1)
  close(12)

  stop
  end
