
!
      subroutine avstrp(x,nx,lstrip,xf,nxf)
!
!   Copyright Goverment of Canada
!   Department of Fisheries and Oceans
!   Institute of Ocean Sciences
!   Sidney, B.C.
!
!   Name: avstrp
!
!   Purpose: To filter using the simple moving average technique.
!   
!   Input:    x      an equally spaced time series.
!             nx     the number of points in X
!             lstrip number of consecutive x values in each average (<=nx)
!             
!   Output:   xf     the lowpassed x values.
!             nxf    number of points in xf.
!
!
!  Written by: MEDS Ottawa
!              Changed to T&C Apollo system by Anne Ballantyne.
!
!	Apollo references deleted by Mike Foreman
!
!      include '/tctree/sources/dsee/tcf/tcf_header_defs.h'
!     Modified by Joseph Zhang Oct 2006
!     Bad values of time series are > 9998 (the filtered results are 9999 in this case)
!     The first point of xf() corresponds to the average of first lstrip points
!     The last point of xf() corresponds to the average of last lstrip points
!     The time stamp of 1st point of xf() is t0+(lstrip-1)*dt/2 (t0 is the starting time)

      integer, intent(in) :: nx,lstrip
      integer, intent(out) :: nxf
      real, intent(in) :: x(nx)
      real, intent(out) :: xf(nx)

      integer mxf,nbad,i,inew,lsm1
      real    big, sum, xold, xnew

      if(nx<lstrip) then
        write(*,*)'Make sure nx >= lstrip in filter avstrp'
        stop
      endif

      lsm1 = lstrip-1
      mxf = nx - lstrip +1
!      big = SHORT_PAD*1.0
      big = 9998.
      nbad = 0
      sum = 0.0
      do i = 1, lsm1
       if(x(i).gt.big) then
         nbad = nbad + 1
       else
         sum = sum + x(i)
       endif
      end do

      xold = 0.0
      do i = 1, mxf
        inew = i+lstrip-1
        xnew = x(inew)
        if(xnew.gt.big) then
          nbad = nbad + 1
        else
          sum = sum + xnew
        endif
        if(xold.gt.big) then
          nbad = nbad - 1
        else
          sum = sum - xold
        endif
        xold = x(i)
        if (nbad.gt.0) then
          xf(i) = big+1
        else
          xf(i) = sum/lstrip
        endif
!	write(6,*) ' i,nbad,xnew,xf(i)=',i,nbad,xnew,xf(i)
      end do

      nxf = mxf

      return
      end
