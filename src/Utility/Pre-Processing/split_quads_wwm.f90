! Split all quads for WWM (compatible with main SCHISM code)
! Inputs: hgrid.gr3; output: hgrid_WWM.gr3 (same if pure tri)

! ifort -O2 -mcmodel=medium -Bstatic -o split_quads_wwm split_quads_wwm.f90

  implicit real*8(a-h,o-z)

  integer :: nwild(10)
  integer, allocatable :: i34(:),elnode(:,:),elnode_wwm(:,:)
  real*8, allocatable :: xnd(:),ynd(:),dp(:)  

  open(14,file='hgrid.gr3',status='old') 
  open(12,file='hgrid_WWM.gr3',status='replace') 
  read(14,*); read(14,*)ne,np
  allocate(i34(ne),elnode(4,ne),xnd(np),ynd(np),dp(np))
  do i=1,np
    read(14,*)j,xnd(i),ynd(i),dp(i)
  enddo !i

  new=0
  do i=1,ne
    read(14,*)j,i34(i),elnode(1:i34(i),i)
    if(i34(i)==4) new=new+1
  enddo !i
  close(14)
  print*, 'Adding ',new,' new elem'

  ne_wwm=ne+new
  allocate(elnode_wwm(3,ne_wwm))

  elnode_wwm(1:3,1:ne)=elnode(1:3,1:ne)
  nwild(1:3)=(/1,3,4/)
  new=0
  do i=1,ne
    if(i34(i)==4) then
      new=new+1
      elnode_wwm(1:3,ne+new)=elnode(nwild(1:3),i)
    endif
  enddo !i
  if(ne+new/=ne_wwm) stop 'mismatch'

! Output
  write(12,*)'WWM pure tri grid'
  write(12,*)ne_wwm,np
  do i=1,np
    write(12,'(i12,3(1x,e22.14))')i,xnd(i),ynd(i),dp(i)
  enddo !i
  do i=1,ne_wwm
    write(12,*)i,3,elnode_wwm(:,i)
  enddo !i

  stop
  end
