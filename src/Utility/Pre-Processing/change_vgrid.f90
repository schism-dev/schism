! Change ivcor=1 format of vgrid.in for newer versions from 11 Oct 2021
! Inputs: hgrid.gr3, vgrid.in.old
! Outputs: vgrid.in.new
! ifort -O2 -mcmodel=medium -CB -Bstatic -o change_vgrid change_vgrid.f90 

  implicit real*8(a-h,o-z)
  integer, allocatable :: kbp(:)
  real*8, allocatable :: sigma(:,:)

  open(14,file='hgrid.gr3',status='old')
  open(19,file='vgrid.in.old',status='old')
  open(18,file='vgrid.in.new',status='replace')
  read(14,*); read(14,*)ne,np
  close(14)

  read(19,*)ivcor
  if(ivcor/=1) stop 'ivcor must be 1'
  read(19,*)nvrt
  allocate(kbp(np),sigma(nvrt,np))
  sigma=-9. !init as any junk value outside [-1,0]
  do i=1,np
    read(19,*)k,kbp(i),sigma(kbp(i):nvrt,i)
  enddo !i
  close(19)

  write(18,*)ivcor; write(18,*)nvrt
  if(np>10000000) stop 'Please increase write length'
  write(18,'(10000000(1x,i10))')kbp(:)
  do k=1,nvrt
    write(18,'(i10,10000000(1x,f14.6))')k,sigma(k,:)
  enddo !k

  stop
  end 
