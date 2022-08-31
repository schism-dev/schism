!   Copyright 2014 College of William and Mary
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

!     Interpolate depths from a structured grid to unstructured grid.
!     Inputs: struc.grd (structured grid in *.asc (ASCII Raster Grid File with,
!             Arc Header format - more restrictive than Mansour's format; see interpolate_depth_structured1)
!             hgrid.old (unstructured grid, mixed tri and quads)
!     Output: hgrid.new (for pts outside the structured grid or the depth is junk there, 
!                        the original depths are preserved).
!     ifort -Bstatic -O3 -o interpolate_depth_structured2 interpolate_depth_structured2.f90
      implicit real*8(a-h,o-z)
!      parameter(mnx=9073) !max. # of pts in x
!      parameter(mny=12961)
      character*5 cha1
      character*9 cha2
      character*12 cha3
      dimension nm(4)
      allocatable :: dp1(:,:)
      
      print*, 'Reverse the sign of the depth? (1: no; -1: yes; say yes)'
!'
      read*, ih
      print*, 'Add vertical const. to outputs (i.e. change of vdatum):'
      read*, vshift
!      print*, 'Adjust 1/2 cell (for corner based .asc)? 1: yes'
!      read*, iadjust_corner

      open(62,file='struc.grd',status='old')
      read(62,*) cha1,nx !# of nodes in x
      read(62,*) cha1,ny !# of nodes in y
      read(62,*) cha2,xmin0
      cha2=adjustl(cha2)
      if(cha2(7:7).eq."n".or.cha2(7:7).eq."N") then !lower-left is corner based
        iadjust_corner=1
      else !center based
        iadjust_corner=0
      endif

      read(62,*) cha2,ymin0
      read(62,*) cha2,dxy
      read(62,*) cha3,fill_value
      dx=dxy
      dy=dxy

!     Define min/max for both vertex and center. '0' means corner/vertex location; otherwise
!     (xmin/xmax etc) are cell center locations.
      if(iadjust_corner/=0) then !corner based
        xmin=xmin0+dx/2
        ymin=ymin0+dy/2
      else !cell center based
        xmin = xmin0
        ymin = ymin0
        xmin0=xmin0-dx/2 !redefine vertex location
        ymin0=ymin0-dy/2
      endif

      allocate(dp1(nx,ny),stat=istat)
      if(istat/=0) stop 'Failed to allocate (1)'

!     Max
      ymax=ymin+(ny-1)*dy !center location
      xmax=xmin+(nx-1)*dx
      xmax0=xmax+dx/2 ! right edge of raster
      ymax0=ymax+dy/2 ! upper edge of raster

!     Read starts from upper left corner and goes along x
      do iy=1,ny
        read(62,*)(dp1(ix,ny-iy+1),ix=1,nx)
!        write(99,*)'line read in:',iy+6
      enddo !iy
      close(62)
   
      open(14,file='hgrid.old',status='old')
      open(13,file='hgrid.new')
      read(14,*)
      read(14,*)ne,np
      write(13,*)'Bathymetry loaded grid'
      write(13,*)ne,np
      do i=1,np
        read(14,*)j,x,y,dp

!       Interpolate
        if(x>xmax0.or.x<xmin0.or.y>ymax0.or.y<ymin0) then
          write(13,101)j,x,y,dp
        else !inside structured grid
          !1/2 cell shift case: extrap to cover lower&left
!          if(iadjust_corner/=0) then
!            x=max(x,xmin)
!            y=max(y,ymin)
!          endif
!          x2=x 
!          y2=y 

          x2=min(xmax,max(x,xmin))
          y2=min(ymax,max(y,ymin))

          ix=(x2-xmin)/dx+1 !i-index of the lower corner of the parent box 
          iy=(y2-ymin)/dy+1
          if(ix.lt.1.or.ix.gt.nx.or.iy.lt.1.or.iy.gt.ny) then
            print*, 'Impossible:',i,ix,iy
            stop
          endif

          if(ix.eq.nx) then !for pts right on the right bnd
            ix=nx-1
            xrat=1
          else
            xrat=(x2-xmin)/dx-ix+1
          endif
          if(iy.eq.ny) then !for pts right on the upper bnd
            iy=ny-1
            yrat=1
          else
            yrat=(y2-ymin)/dy-iy+1
          endif
          if(xrat.lt.0.or.xrat.gt.1.or.yrat.lt.0.or.yrat.gt.1) then
            print*, 'ratios out of bound:',i,xrat,yrat
            stop
          endif
 
          if(abs(dp1(ix,iy)-fill_value)<1.e-2.or.abs(dp1(ix+1,iy)-fill_value)<1.e-2.or. &
     &abs(dp1(ix,iy+1)-fill_value)<1.e-2.or.abs(dp1(ix+1,iy+1)-fill_value)<1.e-2) then
            write(13,101)j,x,y,dp
          else !all valid
            hy1=dp1(ix,iy)*(1-xrat)+xrat*dp1(ix+1,iy)
            hy2=dp1(ix,iy+1)*(1-xrat)+xrat*dp1(ix+1,iy+1)
            h=hy1*(1-yrat)+hy2*yrat
            write(13,101)j,x,y,h*ih+vshift
          endif !junk

        endif
      enddo !i
101   format(i9,2(1x,e24.16),1x,f13.6)

      do i=1,ne
        read(14,*)j,k,(nm(l),l=1,k)
        write(13,*)j,k,(nm(l),l=1,k)
        if(i/=j) then
          write(*,*)'Wrong elem:',i,j,k,nm(1:k)
          stop
        endif
      enddo !i

      stop
      end
