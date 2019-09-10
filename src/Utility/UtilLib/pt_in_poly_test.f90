! Modules for point in polgon tests

!  signa_[single,double]
!  pt_in_poly_[single,double]: point-in-polygon (tri/quad) test
!  pt_in_poly_ray_method_[single,double]: point-in-polygon test using ray tracing method (double prevision). 
!====================================================================
  module pt_in_poly_test
   implicit none

   public

   contains

      function signa_single(x1,x2,x3,y1,y2,y3)
!...  Compute signed area formed by pts 1,2,3
      real :: signa_single,x1,x2,x3,y1,y2,y3

      signa_single=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2

      end function signa_single

      function signa_double(x1,x2,x3,y1,y2,y3)
!...  Compute signed area formed by pts 1,2,3
      real(kind=8) :: signa_double,x1,x2,x3,y1,y2,y3

      signa_double=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2

      end function signa_double

      subroutine pt_in_poly_single(i34,x,y,xp,yp,inside,arco,nodel)
!     (Single-precision) Routine to perform point-in-polygon
!     (triangle/quads) test and if it's inside, calculate the area coord.
!     (for quad, split it into 2 triangles and return the 3 nodes and
!     area coord.)
!     Inputs:
!            i34: 3 or 4 (type of elem)
!            x(i34),y(i34): coord. of polygon/elem. (counter-clockwise)
!            xp,yp: point to be tested
!     Outputs:
!            inside: 0, outside; 1, inside
!            arco(3), nodel(3) : area coord. and 3 local node indices (valid only if inside)
      implicit real(a-h,o-z), integer(i-n)
      integer, intent(in) :: i34
      real, intent(in) :: x(i34),y(i34),xp,yp
      integer, intent(out) :: inside,nodel(3)
      real, intent(out) :: arco(3)

      !Local
      integer :: list(3)
      real :: ar(2),swild(2,3)

      !Areas
      ar(1)=signa_single(x(1),x(2),x(3),y(1),y(2),y(3))
      ar(2)=0 !init
      if(i34==4) ar(2)=signa_single(x(1),x(3),x(4),y(1),y(3),y(4))
      if(ar(1)<=0.or.i34==4.and.ar(2)<=0) then
        print*, 'Negative area:',i34,ar,x,y
        stop
      endif

      zero=0; one=1
      inside=0
      do m=1,i34-2 !# of triangles
        if(m==1) then
          list(1:3)=(/1,2,3/) !local indices
        else !quads
          list(1:3)=(/1,3,4/)
        endif !m
        aa=0
        do j=1,3
          j1=j+1
          j2=j+2
          if(j1>3) j1=j1-3
          if(j2>3) j2=j2-3
          swild(m,j)=signa_single(x(list(j1)),x(list(j2)),xp,y(list(j1)),y(list(j2)),yp) !temporary storage
          aa=aa+abs(swild(m,j))
        enddo !j=1,3

        ae=abs(aa-ar(m))/ar(m)
        if(ae<=1.e-5) then
          inside=1
          nodel(1:3)=list(1:3)
          arco(1:3)=swild(m,1:3)/ar(m)
          arco(1)=max(zero,min(one,arco(1)))
          arco(2)=max(zero,min(one,arco(2)))
          if(arco(1)+arco(2)>1) then
            arco(3)=0
            arco(2)=1-arco(1)
          else
            arco(3)=1-arco(1)-arco(2)
          endif
          exit
        endif
      enddo !m

      end subroutine pt_in_poly_single

      subroutine pt_in_poly_double(i34,x,y,xp,yp,inside,arco,nodel)
!     (Double-precision) Routine to perform point-in-polygon
!     (triangle/quads) test and if it's inside, calculate the area coord.
!     (for quad, split it into 2 triangles and return the 3 nodes and
!     area coord.)
!     Inputs:
!            i34: 3 or 4 (type of elem)
!            x(i34),y(i34): coord. of polygon/elem. (counter-clockwise)
!            xp,yp: point to be tested
!     Outputs:
!            inside: 0, outside; 1, inside
!            arco(3), nodel(3) : area coord. and 3 local node indices (valid only if inside)
      implicit real*8(a-h,o-z), integer(i-n)
      integer, intent(in) :: i34
      real*8, intent(in) :: x(i34),y(i34),xp,yp
      integer, intent(out) :: inside,nodel(3)
      real*8, intent(out) :: arco(3)

      !Local
      integer :: list(3)
      real*8 :: ar(2),swild(2,3)

      !Areas
      ar(1)=signa_double(x(1),x(2),x(3),y(1),y(2),y(3))
      ar(2)=0 !init
      if(i34==4) ar(2)=signa_double(x(1),x(3),x(4),y(1),y(3),y(4))
      if(ar(1)<=0.or.i34==4.and.ar(2)<=0) then
        print*, 'Negative area:',i34,ar,x,y
        stop
      endif

      zero=0; one=1
      inside=0
      do m=1,i34-2 !# of triangles
        if(m==1) then
          list(1:3)=(/1,2,3/) !local indices
        else !quads
          list(1:3)=(/1,3,4/)
        endif !m
        aa=0
        do j=1,3
          j1=j+1
          j2=j+2
          if(j1>3) j1=j1-3
          if(j2>3) j2=j2-3
          swild(m,j)=signa_double(x(list(j1)),x(list(j2)),xp,y(list(j1)),y(list(j2)),yp) !temporary storage
          aa=aa+abs(swild(m,j))
        enddo !j=1,3

        ae=abs(aa-ar(m))/ar(m)
        if(ae<=1.e-5) then
          inside=1
          nodel(1:3)=list(1:3)
          arco(1:3)=swild(m,1:3)/ar(m)
          arco(1)=max(zero,min(one,arco(1)))
          arco(2)=max(zero,min(one,arco(2)))
          if(arco(1)+arco(2)>1) then
            arco(3)=0
            arco(2)=1-arco(1)
          else
            arco(3)=1-arco(1)-arco(2)
          endif
          exit
        endif
      enddo !m

      end subroutine pt_in_poly_double

!====================================================================
!====================================================================

!     Routines to perform point-in-polygon
!     test. Polygon may be multi-connected (and concave), and
!     the list of vertices can be either clock- or counter-clockwise.
!     The 1st vertex of each poly must be repeated, and 1 sub-polygon must
!     completely enclose all others, and sub-poly's do not cross each other.
!     Inputs:
!            nvertices: # of polygonal vertices. 1st vertex of each sub-poly
!            must be repeated.
!            small1,ray_angle: tolerance and search angle (see below)
!            xpoly(nvertices),ypoly(nvertices): list of vertices
!            xtest,ytest: test point
!     Outputs:
!            in_out: -1 if outside, 0 if on the polygon, 1 if inside
!            inters: # of intersections (no accurate if in_out=0). If 
!            outside bounding box of poly, inters=-1.
!            npoly: # of sub-poly's found
!     Error: fort.11 (fatal)
!     Method: ray tracing method by testing each side of polygons and a
!     ray from the test point to infty along an angle. The point is
!     outside if the ray intersects polygonal sides even # of times,
!     inside if odd # of times. Check separately if on the side based
!     on a tolerance (small1).
!     Assume that the ray is not parallel to any side and also the
!     intersection does not coincide with any vertices if the test pt
!     is not one of them. If fatal error occurs, adjust
!     ray_angle or small1.

      subroutine pt_in_poly_ray_method_double(nvertices,small1,ray_angle,xpoly,ypoly,xtest,ytest,in_out,inters,npoly)
      implicit real*8(a-h,o-z), integer(i-n)
      integer, intent(in) :: nvertices
      real*8, intent(in) :: small1,ray_angle,xpoly(nvertices),ypoly(nvertices),xtest,ytest
      integer, intent(out) :: in_out,inters,npoly

      !Local
      real*8, parameter :: pi=3.1415926d0
      !real*8, parameter :: small1=1.e-5 !tolerance
      !real*8 :: ray_angle=3.13192 !degr - use a real number to enhance randomness
      integer :: list(3),firstv_poly(nvertices),endv_poly(nvertices)

      include 'pt_in_poly_ray_method.txt'

      end subroutine pt_in_poly_ray_method_double

      subroutine pt_in_poly_ray_method_single(nvertices,small1,ray_angle,xpoly,ypoly,xtest,ytest,in_out,inters,npoly)
      implicit real(a-h,o-z), integer(i-n)
      integer, intent(in) :: nvertices
      real, intent(in) :: small1,ray_angle,xpoly(nvertices),ypoly(nvertices),xtest,ytest
      integer, intent(out) :: in_out,inters,npoly

      !Local
      real, parameter :: pi=3.1415926
      !real, parameter :: small1=1.e-5 !tolerance
      !real :: ray_angle=3.13192 !degr - use a real number to enhance randomness
      integer :: list(3),firstv_poly(nvertices),endv_poly(nvertices)

      include 'pt_in_poly_ray_method.txt'

      end subroutine pt_in_poly_ray_method_single

!====================================================================
  end module pt_in_poly_test

!     Test driver
!     ifort -Bstatic -O2 -mcmodel=medium -o 
!      program test_pt_in_poly
!      implicit real*8(a-h,o-z)
!      real*8, allocatable :: xpoly(:),ypoly(:)
!
!      !Use gredit/regions to generate a multi-connected poly
!      !Make sure to repeat 1st vertex of each sub-poly
!      open(12,file='multi.xy',status='old')
!      open(10,file='test_pts.bp',status='old')
!      read(12,*)nvertices
!      allocate(xpoly(nvertices),ypoly(nvertices))
!      do i=1,nvertices
!        read(12,*)xpoly(i),ypoly(i)
!      enddo !i
!      close(12)
!      read(10,*); read(10,*)npts
!      do i=1,npts
!        read(10,*)j,xtest,ytest
!
!        !Debug
!        write(98,*)'Test pt #',i
!
!        call pt_in_poly_ray_method(nvertices,xpoly,ypoly,xtest,ytest,in_out,inters)
!        print*, 'Point ',i,in_out,inters
!      enddo !i
!
!      !Test a pt on poly
!      xtest=0.53*xpoly(nvertices)+(1-0.53)*xpoly(nvertices-1)
!      ytest=0.53*ypoly(nvertices)+(1-0.53)*ypoly(nvertices-1)
!      call pt_in_poly_ray_method(nvertices,xpoly,ypoly,xtest,ytest,in_out,inters,npoly)
!      print*, 'Extra pt on poly:',in_out,inters
!
!      end program test_pt_in_poly
