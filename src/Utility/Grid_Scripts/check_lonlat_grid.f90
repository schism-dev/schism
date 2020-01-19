! Check quality of lon/lat grid (including global with or w/o poles), b/cos
! gredit has issues (e.g. hi distortion near poles), using SCHISM's 3D coordinates on an ellipsoid.
! Fix negative elements and output list of skew elem.
!
! Inputs: screen; hgrid.ll (tri only)
! Outputs: hgrid.ll.new; fort.99 (diagnostics)

! ifort -O2 -mcmodel=medium -CB -Bstatic -o check_lonlat_grid check_lonlat_grid.f90

  implicit real*8(a-h,o-z)
  real*8, parameter :: pi=3.1415926d0
  real*8, parameter :: deg2rad=pi/180.d0
  real*8, parameter :: rad2deg=180.d0/pi
  real*8, parameter :: rearth_eq=6378206.4 ![m]
  real*8, parameter :: rearth_pole=6378206.4

  real*8 :: xx(4),yy(4),swild(4)
  integer, allocatable :: i34(:),elnode(:,:),i34_new(:),elnode_new(:,:)
  real*8, allocatable :: xlon(:),ylat(:),dp(:),area(:),quad_loc(:,:), &
 &xctr(:),yctr(:),xnd(:),ynd(:),znd(:),pframe(:,:,:)

  print*, 'Input max allowable skewness for triangles:'
  read*, skew_max
  
  open(14,file='hgrid.ll',status='old')
  read(14,*); read(14,*)ne,np

  allocate(xlon(np),ylat(np),dp(np),area(ne),i34(ne),elnode(4,ne), &
 &i34_new(2*ne),elnode_new(4,2*ne),quad_loc(2,ne),xctr(ne),yctr(ne),xnd(np),ynd(np),znd(np), &
 &pframe(3,3,np))

  do i=1,np
    read(14,*) j,xlon(i),ylat(i),dp(i)
    xtmp=xlon(i)*deg2rad !to rad
    ytmp=ylat(i)*deg2rad !to rad
    xnd(i)=rearth_eq*cos(ytmp)*cos(xtmp)
    ynd(i)=rearth_eq*cos(ytmp)*sin(xtmp)
    znd(i)=rearth_pole*sin(ytmp)
  enddo !i

  !Compute pframe
  do i=1,np
    xtmp=xlon(i)*deg2rad !to rad
    ytmp=ylat(i)*deg2rad
    pframe(1,1,i)=-sin(xtmp) !zonal dir.
    pframe(2,1,i)=cos(xtmp)
    pframe(3,1,i)=0
    pframe(1,2,i)=-cos(xtmp)*sin(ytmp) !meri. dir.
    pframe(2,2,i)=-sin(xtmp)*sin(ytmp)
    pframe(3,2,i)=cos(ytmp)
    call cross_product(pframe(1,1,i),pframe(2,1,i),pframe(3,1,i), &
   &pframe(1,2,i),pframe(2,2,i),pframe(3,2,i),pframe(1,3,i),pframe(2,3,i),pframe(3,3,i))
  enddo !i=1,npa

  do i=1,ne
    read(14,*) j,i34(i),elnode(1:i34(i),i)
    if(i34(i)/=3) stop 'No quads plz'
  enddo !i
  close(14)  
  
  i34_new(1:ne)=i34
  elnode_new(:,1:ne)=elnode

  do i=1,ne
    !Project to a local plane @ node 1
    nd1=elnode(1,i)
    do j=1,i34(i)
      nd=elnode(j,i)
      call project_pt('g2l',xnd(nd),ynd(nd),znd(nd),(/xnd(nd1),ynd(nd1),znd(nd1)/),pframe(:,:,nd1),xx(j),yy(j),tmp)
    enddo !j

    do j=1,i34(i)
      j1=j+1
      if(j1>i34(i)) j1=j1-i34(i)
      swild(j)=sqrt((xx(j)-xx(j1))**2+(yy(j)-yy(j1))**2) !side length
    enddo !j

    n1=elnode(1,i)
    n2=elnode(2,i)
    n3=elnode(3,i)
    area(i)=signa(xx(1),xx(2),xx(3),yy(1),yy(2),yy(3))

    if(i34(i)==3) then !tria
      ifl=0
      if(area(i)<0) then
        write(*,*)'Fixing negative area at elem:',i,area(i)
        write(99,*)'Fixing negative area at elem:',i,area(i)
        elnode_new(1,i)=elnode(2,i)
        elnode_new(2,i)=elnode(1,i)
        area(i)=-area(i)
      else if(area(i)==0) then
        ifl=1
        write(99,*)'Elem_has_0_area:',i,xlon(i),ylat(i)
      endif

      !Skewness
      if(ifl==0) then
        skew=maxval(swild(1:3))/sqrt(area(i)/pi)
        if(skew>skew_max) write(99,*)'Skew_elem:',i,xlon(i),ylat(i)
      endif !ifl
    !else if(i34(i)==4) then !quad
      
    else 
      stop 'unknown elem type'
    endif !i34(i)
  enddo !i=1,ne

! Output
  open(12,file='hgrid.ll.new',status='replace')
  write(12,*)
  write(12,*)ne,np
  do i=1,np
    write(12,*)i,real(xlon(i)),real(ylat(i)),real(dp(i))
  enddo !i
  do i=1,ne
    write(12,*)i,i34_new(i),elnode_new(1:i34_new(i),i)
  enddo !i
  close(12)

  print*, 'finished'

  stop
  end


      function signa(x1,x2,x3,y1,y2,y3)
!...  Compute signed area formed by pts 1,2,3
      implicit real*8(a-h,o-z)

      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2

      return
      end

      subroutine project_pt(dir,xi,yi,zi,origin0,frame0,xo,yo,zo)
      implicit real*8(a-h,o-z)

      character(len=3), intent(in) :: dir
      real(8), intent(in) :: xi,yi,zi,origin0(3),frame0(3,3)
      real(8), intent(out) :: xo,yo,zo

      !Local
      real(8) :: wild(3)

      if(dir.eq.'g2l') then
        wild(1:3)=(xi-origin0(1))*frame0(1,1:3)+(yi-origin0(2))*frame0(2,1:3)+ &
                 &(zi-origin0(3))*frame0(3,1:3)
      else if(dir.eq.'l2g') then
        wild(1:3)=origin0(1:3)+xi*frame0(1:3,1)+yi*frame0(1:3,2)+ &
     &zi*frame0(1:3,3)
      else
        stop 'PROJECT_PT: unknown flag'
      endif
      xo=wild(1)
      yo=wild(2)
      zo=wild(3)

      end subroutine project_pt

      subroutine cross_product(x1,y1,z1,x2,y2,z2,x3,y3,z3)
      implicit real*8(a-h,o-z)
      real(8),intent(in) :: x1,y1,z1,x2,y2,z2
      real(8),intent(out) :: x3,y3,z3

      x3=y1*z2-y2*z1
      y3=x2*z1-x1*z2
      z3=x1*y2-x2*y1

      end subroutine cross_product

