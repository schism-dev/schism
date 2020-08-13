!     Interpolate depths from a structured grid to unstructured grid.
!     Inputs: struc.grd (structured grid in *.asc (ASCII Raster Grid File with,
!             Arc Header format - more restrictive than Mansour's format; see interpolate_depth_structured1)
!             hgrid.old (unstructured grid)
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
      logical :: exist
      
      print*, 'Reverse the sign of the depth? (1: no; -1: yes; say yes)'
!'
      read*, ih
      print*, 'Add vertical const. to outputs (i.e. change of vdatum):'
      read*, vshift

      open(62,file='struc.grd',status='old')
      read(62,*) cha1,nx !# of nodes in x
      read(62,*) cha1,ny !# of nodes in y
      read(62,*) cha2,xmin
      read(62,*) cha2,ymin
      read(62,*) cha2,dxy
      read(62,*) cha3,fill_value
      dx=dxy
      dy=dxy
      tol= 0.5 * dxy

!      if(nx.gt.mnx.or.ny.gt.mny) then
!        print*, 'Increase mnx,y to ',nx,ny
!        stop
!      endif
      allocate(dp1(nx,ny),stat=istat)
      if(istat/=0) stop 'Failed to allocate (1)'

!     Coordinates for upper left corner (the starting point for *.asc)
      ymax=ymin+(ny-1)*dy
!     xmax
      xmax=xmin+(nx-1)*dx

!      write(12,*)1,xmin,ymin
!      write(12,*)2,xmax,ymax
!      stop

      do iy=1,ny
        read(62,*)(dp1(ix,ny-iy+1),ix=1,nx)
      enddo !iy
      close(62)

      inquire(file="test.txt", exist=exist)
      if (exist) then
        open(19, file="test.txt", status="old", position="append", action="write")
      else
        open(19, file="test.txt", status="new", action="write")
      endif
   
      open(14,file='hgrid.old',status='old')
      open(13,file='hgrid.new')
      read(14,*)
      read(14,*)ne,np
      write(13,*)'Bathymetry loaded grid'
      write(13,*)ne,np
      do i=1,np
        read(14,*)j,x,y,dp

!       Interpolate
        if(x.gt.(xmax+tol).or.x.lt.(xmin-tol).or.y.gt.(ymax+tol).or.y.lt.(ymin-tol)) then
          write(13,101)j,x,y,dp
        else !inside structured grid
          x2=x 
          y2=y 
          ix=(x2-xmin)/dx+1 !i-index of the lower corner of the parent box 
          iy=(y2-ymin)/dy+1
          if(ix.lt.1.or.ix.gt.nx.or.iy.lt.1.or.iy.gt.ny) then
            write(19,*) 'Interface padded:i,ix,iy,x2,y2,xmin,ymin,xmax,ymax: ',i,ix,iy,x2,y2,xmin,ymin,xmax,ymax
            x2=min(max(xmin,x2),xmax)
            y2=min(max(ymin,y2),ymax)
            ix=(x2-xmin)/dx+1 !i-index of the lower corner of the parent box 
            iy=(y2-ymin)/dy+1
          endif
          if(ix.lt.1.or.ix.gt.nx.or.iy.lt.1.or.iy.gt.ny) then
            write(19,*)  'impossible:',i,ix,iy,x2,y2,xmin,ymin,xmax,ymax
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
            write(19,*) 'ratios out of bound:',i,xrat,yrat
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
