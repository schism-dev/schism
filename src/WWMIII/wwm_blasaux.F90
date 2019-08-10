#include "wwm_functions.h"
!-----------------------------------------------------------------------
!     subroutine from blas1.f90
!-----------------------------------------------------------------------
      function dnrm2 ( n, dx, incx)
!Mathieu here zero was not defined and implicit none not specified ...
!        major bug!
      use datapool, only : rkind, THR, zero, one
      real(rkind) :: dnrm2
      integer, intent(in) :: incx
      integer          next
      real(rkind)   dx(*), cutlo, cuthi, hitest, sumd, xmax
      integer n, nn, i, j
!
!     euclidean norm of the n-vector stored in dx() with storage
!     increment incx .
!     if    n .le. 0 return with result = 0.
!     if n .ge. 1 then incx must be .ge. 1
!
!           c.l.lawson, 1978 jan 08
!
!     four phase method     using two built-in constants that are
!     hopefully applicable to all machines.
!         cutlo = maximum of  dsqrt(u/eps)  over all known machines.
!         cuthi = minimum of  dsqrt(v)      over all known machines.
!     where
!         eps = smallest no. such that eps + 1. .gt. 1.
!         u   = smallest positive no.   (underflow limit)
!         v   = largest  no.            (overflow  limit)
!
!     brief outline of algorithm..
!
!     phase 1    scans zero components.
!     move to phase 2 when a component is nonzero and .le. cutlo
!     move to phase 3 when a component is .gt. cutlo
!     move to phase 4 when a component is .ge. cuthi/m
!     where m = n for x() real and m = 2*n for complex.
!
!     values for cutlo and cuthi..
!     from the environmental parameters listed in the imsl converter
!     document the limiting values are as follows..
!     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
!                   univac and dec at 2**(-103)
!                   thus cutlo = 2**(-51) = 4.44089e-16
!     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.
!                   thus cuthi = 2**(63.5) = 1.30438e19
!     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.
!                   thus cutlo = 2**(-33.5) = 8.23181d-11
!     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19
!     data cutlo, cuthi / 8.232d-11,  1.304d19 /
!     data cutlo, cuthi / 4.441e-16,  1.304e19 /
      data cutlo, cuthi / 8.232d-11,  1.304d19 /
!
      if(n .gt. 0) go to 10
         dnrm2  = zero
         go to 300
!
   10 assign 30 to next
      sumd = zero
      nn = n * incx
!                                                 begin main loop
      i = 1
   20    go to next,(30, 50, 70, 110)
   30 if( MyABS(dx(i)) .gt. cutlo) go to 85
      assign 50 to next
      xmax = zero
!
!                        phase 1.  sumd is zero
!
   50 if( abs(dx(i)) .lt. THR) go to 200
      if( MyABS(dx(i)) .gt. cutlo) go to 85
!
!                                prepare for phase 2.
      assign 70 to next
      go to 105
!
!                                prepare for phase 4.
!
  100 i = j
      assign 110 to next
      sumd = (sumd / dx(i)) / dx(i)
  105 xmax = MyABS(dx(i))
      go to 115
!
!                   phase 2.  sumd is small.
!                             scale to avoid destructive underflow.
!
   70 if( MyABS(dx(i)) .gt. cutlo ) go to 75
!
!                     common code for phases 2 and 4.
!                     in phase 4 sumd is large.  scale to avoid overflow.
!
  110 if( MyABS(dx(i)) .le. xmax ) go to 115
         sumd = one + sumd * (xmax / dx(i))**2
         xmax = MyABS(dx(i))
         go to 200
!
  115 sumd = sumd + (dx(i)/xmax)**2
      go to 200
!
!
!                  prepare for phase 3.
!
   75 sumd = (sumd * xmax) * xmax
!
!
!     for real or d.p. set hitest = cuthi/n
!     for complex      set hitest = cuthi/(2*n)
!
   85 hitest = cuthi/float( n )
!
!                   phase 3.  sumd is mid-range.  no scaling.
!
      do 95 j =i,nn,incx
      if(MyABS(dx(j)) .ge. hitest) go to 100
   95    sumd = sumd + dx(j)**2
      dnrm2 = MySQRT( sumd )
      go to 300
!
  200 continue
      i = i + incx
      if ( i .le. nn ) go to 20
!
!              end of main loop.
!
!              compute square root and adjust for scaling.
!
      dnrm2 = xmax * MySQRT(sumd)
  300 continue
      return
      end
!-------------------------------------------------------------------------
#ifndef SCHISM
      real(rkind) function ddot(n,dx,incx,dy,incy)
      use datapool, only : rkind, ZERO
!
!     forms the dot product of two vectors.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!
      real(rkind) dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
!
      ddot = ZERO
      dtemp = ZERO
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) + &
     &   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end
#endif
!----------------------------------------------------------------------
      subroutine daxpy(n,da,dx,incx,dy,incy)
      use datapool, only : rkind, THR
!
!     constant times a vector plus a vector.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!
      real(rkind) dx(1),dy(1),da
      integer i,incx,incy,ix,iy,m,mp1,n
!
      if(n.le.0)return
      if (abs(da) .lt. THR) return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end
