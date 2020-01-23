! Augument fixed.gr3 (input to grid_spring.f90) with regions surrounding
! stations

! Inputs: fixed.gr3.0, station.gr3 (depth is radius desired; connectivity table
! not used)
! Outputs: fixed_new.gr3

! ifort -O2 -mcmodel=medium -CB -Bstatic -o augment_fixed_grid augment_fixed_grid.f90

  implicit real*8(a-h,o-z)
  integer :: nm(4)
  integer, allocatable :: i34(:),elnode(:,:),ifixed(:)
  real*8, allocatable :: x(:),y(:),xsta(:),ysta(:),radius(:)

  open(10,file='fixed.gr3.0',status='old')
  open(11,file='station.gr3',status='old')
  read(10,*); read(10,*)ne,np
  allocate(x(np),y(np),ifixed(np),i34(ne),elnode(4,ne))
  do i=1,np
    read(10,*) j,x(i),y(i),tmp
    ifixed(i)=nint(tmp)
  enddo !i=1,np
  do i=1,ne
    read(10,*)j,i34(i),elnode(1:i34(i),i)
  enddo !i
  close(10)
  
  read(11,*); read(11,*)ntmp,nsta
  allocate(xsta(nsta),ysta(nsta),radius(nsta))
  do i=1,nsta
    read(11,*)j,xsta(i),ysta(i),radius(i)
  enddo !i
  close(11)

  !Add to fixed nodes
  do i=1,np
    if(ifixed(i)>0) cycle

    do j=1,nsta
      rl=(x(i)-xsta(j))**2+(y(i)-ysta(j))**2
      if(rl<radius(j)*radius(j)) then
        ifixed(i)=1; exit
      endif
    enddo !j
  enddo !i

  !Output
  open(12,file='fixed_new.gr3',status='replace')
  write(12,*); write(12,*)ne,np
  do i=1,np
    write(12,*)i,x(i),y(i),ifixed(i)
  enddo !i=1,np
  do i=1,ne
    write(12,*)i,i34(i),elnode(1:i34(i),i)
  enddo !i

  stop
  end
