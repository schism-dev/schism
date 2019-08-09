!  pt_in_poly: point-in-polygon (tri/quad) test. For general polygon,
!              see pt_in_poly2.f90 (single-precision)
!  signa (single-precision)
!====================================================================
!====================================================================

      subroutine pt_in_poly(i34,x,y,xp,yp,inside,arco,nodel)
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
!      implicit real*8(a-h,o-z)
      integer, intent(in) :: i34
      real, intent(in) :: x(i34),y(i34),xp,yp
      integer, intent(out) :: inside,nodel(3)
      real, intent(out) :: arco(3)

      !Local
      integer :: list(3)
      real :: ar(2),swild(2,3)

      !Areas
      ar(1)=signa(x(1),x(2),x(3),y(1),y(2),y(3))
      ar(2)=0 !init
      if(i34==4) ar(2)=signa(x(1),x(3),x(4),y(1),y(3),y(4))
      if(ar(1)<=0.or.i34==4.and.ar(2)<=0) then
        print*, 'Negative area:',i34,ar,x,y
        stop
      endif

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
          swild(m,j)=signa(x(list(j1)),x(list(j2)),xp,y(list(j1)),y(list(j2)),yp) !temporary storage
          aa=aa+abs(swild(m,j))
        enddo !j=1,3

        ae=abs(aa-ar(m))/ar(m)
        if(ae<=1.e-5) then
          inside=1
          nodel(1:3)=list(1:3)
          arco(1:3)=swild(m,1:3)/ar(m)
          arco(1)=max(0.,min(1.,arco(1)))
          arco(2)=max(0.,min(1.,arco(2)))
          if(arco(1)+arco(2)>1) then
            arco(3)=0
            arco(2)=1-arco(1)
          else
            arco(3)=1-arco(1)-arco(2)
          endif
          exit
        endif
      enddo !m

      end subroutine pt_in_poly

!====================================================================
      function signa(x1,x2,x3,y1,y2,y3)
!...  Compute signed area formed by pts 1,2,3

      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2

      return
      end
!====================================================================

