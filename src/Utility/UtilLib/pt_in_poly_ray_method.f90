!  pt_in_poly_ray_method: point-in-polygon test using ray tracing method (double prevision). 
!====================================================================
!====================================================================

      subroutine pt_in_poly_ray_method(nvertices,small1,ray_angle,xpoly,ypoly,xtest,ytest,in_out,inters,npoly)
!     (Double-precision) Routine to perform point-in-polygon
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
      implicit real*8(a-h,o-z)
      integer, intent(in) :: nvertices
      real*8, intent(in) :: xpoly(nvertices),ypoly(nvertices),xtest,ytest
      integer, intent(out) :: in_out,inters,npoly

      !Local
      real*8, parameter :: pi=3.1415926d0
      !real*8, parameter :: small1=1.e-5 !tolerance
      !real*8 :: ray_angle=3.13192 !degr - use a real number to enhance randomness
      integer :: list(3),firstv_poly(nvertices),endv_poly(nvertices)

      !Find # of poly's
      endv_poly=-99 !flag
      npoly=1
      firstv_poly(npoly)=1 !1st vertex 
      i=2 !starting vertex for searching
      do 
        j=firstv_poly(npoly)
        rl=(xpoly(i)-xpoly(j))**2+(ypoly(i)-ypoly(j))**2
        if(rl<small1*1.d-3) then
          endv_poly(npoly)=i
          if(i+1==nvertices) then
            write(11,*)'pt_in_poly: last poly has only 1 vertex'
            stop
          else if(i+1<nvertices) then
            npoly=npoly+1
            firstv_poly(npoly)=i+1
            i=i+1 !skip 1st vertex in the next iteration
          endif
        endif

        i=i+1
        if(i>nvertices) exit
      enddo !i

      !print*, npoly,' polygons found'
 
      do i=1,npoly
        if(endv_poly(i)<=0) then
          write(11,*)'pt_in_poly: end not closed, ',i,firstv_poly(i)
          stop
        endif
        if(endv_poly(i)<firstv_poly(i)+3) then
          print*, 'pt_in_poly: improper poly, ',i,firstv_poly(i),endv_poly(i)
          stop
        endif

        !Debug
        !write(40+i,*)
        !write(40+i,*)endv_poly(i)-firstv_poly(i)+1
        !do j=firstv_poly(i),endv_poly(i)
        !  write(40+i,*)j,real(xpoly(j)),real(ypoly(j)),real(i)
        !enddo !j
      enddo !i

      in_out=-1 !init as outside; change sign each time the ray crosses a side
      inters=0
      !Bounding box
      xmin=minval(xpoly); xmax=maxval(xpoly)
      ymin=minval(ypoly); ymax=maxval(ypoly)
      if(xtest<xmin.or.xtest>xmax.or.ytest<ymin.or.ytest>ymax) then
        !outside
        in_out=-1; inters=-1
        return 
      endif

      an2=ray_angle/180*pi
      cosa=cos(an2); sina=sin(an2)
      loop1: do i=1,npoly
        do j=firstv_poly(i),endv_poly(i)-1
          j2=j+1
          !Calc intersection btw 2 lines:
          !L1: x=x_0+(cos,sin)*t1 (t1>=0) (x_0 is the test pt)
          !L2: x=x_1+(x_2-x_1)*t2 (0<=t2<=1) (x_{1,2} are side pts)
          x0=xtest; y0=ytest
          x1=xpoly(j); x2=xpoly(j2)
          y1=ypoly(j); y2=ypoly(j2)
          delta=(y1-y2)*cosa-(x1-x2)*sina
          if(delta==0) then
            write(11,*)'pt_in_poly, delta=0:',i,j,delta
            stop
          endif
          t1=((x1-x0)*(y1-y2)-(y1-y0)*(x1-x2))/delta
          t2=((y1-y0)*cosa-(x1-x0)*sina)/delta
          if(t2>=0.and.t2<=1) then
            if(abs(t1)<small1) then !on the side
              in_out=0; exit loop1
            else if(t1>0) then
              if(t2==0.or.t2==1) then !ambiguous intersection (with at least 2 sides)
                write(11,*)'pt_in_poly: plz adjust ray_angle:',i,j,t1,t2
                stop
              endif

              !Exactly 1 intersection
              in_out=-in_out
              inters=inters+1

              !Debug
              !xin=x0+cosa*t1
              !yin=y0+sina*t1
              !write(98,*)inters,real(xin),real(yin),i,j,j2

            endif
          endif !t2
        enddo !j
      end do loop1 !i=1,npoly

      end subroutine pt_in_poly_ray_method

!====================================================================

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
