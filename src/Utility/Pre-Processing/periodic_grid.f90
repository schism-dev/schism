!     Generate a one-way slant rectangular grid, and add wrap-around elem. for sphere.
!     This is the way to use periodic b.c. along the x-direction in
!     SCHISM (pay attention to Coriolis as well).
!     Inputs: constants below
!     Outputs: hgrid.ll
!     ifort -Bstatic -O3 -CB -o periodic_grid periodic_grid.f90
      implicit real*8(a-h,o-z)
      !# of divisions in latitude, min/max latitudes (degrees)
      parameter(ny=3,ymin=-0.15,ymax=0.15,maxdim=9000000)
      dimension x(maxdim),y(maxdim),ic1(maxdim,3),ic2(maxdim,20)
      dimension nsur(maxdim),x3(maxdim),y3(maxdim),z3(maxdim)
      dimension x3_0(maxdim),y3_0(maxdim),z3_0(maxdim)

      !Const
      !radii at poles and equator [m]. Adjust these to get right
      !distance in the actual grid. These 2 constants should be
      !consistent with those in param.in (rearth_pole,rearth_eq)
      r_pole=6378206.4d0 
      r_eq=r_pole/2 
      !Resolution in longitude [degrees] - cannot be coarser than 0.5 deg(?)
      dx=0.1 
      xmin=-180 
      xmax=180-dx !leave a gap at dateline to be filled later
      nx=(xmax-xmin)/dx+0.1
      pi=acos(-1.d0)

      print*, 'nx=',nx

      if((nx+1)*(ny+1).gt.maxdim.or.nx*ny*2.gt.maxdim) then
        write(*,*)'Increase maxdim'
        stop
      endif

      open(12,file='hgrid.ll')
      do ix=1,nx+1
      do iy=1,ny+1
        nd=(ny+1)*(ix-1)+iy
        x(nd)=xmin+(xmax-xmin)/nx*(ix-1)
        y(nd)=ymin+(ymax-ymin)/ny*(iy-1)
      enddo
      enddo
      nn=(nx+1)*(ny+1)

      do ix=1,nx
      do iy=1,ny
         ne=ny*(ix-1)+iy
         nd1=(ny+1)*(ix-1)+iy
         nd2=(ny+1)*ix+iy
         ic1(2*ne-1,1)=nd1
         ic1(2*ne-1,2)=nd2
         ic1(2*ne-1,3)=nd2+1
         ic1(2*ne,1)=nd1
         ic1(2*ne,2)=nd2+1
         ic1(2*ne,3)=nd1+1
      enddo
      enddo
      nt=nx*ny*2
      nt0=nt

!     Add wrap-around elements
      do iy=1,ny
        nt=nt+2
        nd1=nn-ny+iy-1
        nd2=iy
        ic1(nt-1,1)=nd1
        ic1(nt-1,2)=nd2
        ic1(nt-1,3)=nd2+1
        ic1(nt,1)=nd1
        ic1(nt,2)=nd2+1
        ic1(nt,3)=nd1+1
      enddo !iy

      write(12,*)'hgrid.ll'
      write(12,*)nt,nn
      do i=1,nn
        write(12,100)i,x(i),y(i),50.

        !global coord.
        x3(i)=r_eq*cos(y(i)/180*pi)*cos(x(i)/180*pi)
        y3(i)=r_eq*cos(y(i)/180*pi)*sin(x(i)/180*pi)
        z3(i)=r_pole*sin(y(i)/180*pi)
        write(98,*)i,x3(i),y3(i),z3(i)
      enddo !i
      write(12,101)(j,3,(ic1(j,i),i=1,3),j=1,nt)

      !Calc distances
      rlmin=huge(1.d0)
      rlmax=-huge(1.d0)
      do i=1,nt0
        do j=1,3
          j1=j+1
          if(j1>3) j1=j1-3
          n1=ic1(i,j); n2=ic1(i,j1)
          rl=sqrt((x3(n1)-x3(n2))**2+(y3(n1)-y3(n2))**2+(z3(n1)-z3(n2))**2)
          rlmin=min(rl,rlmin)
          rlmax=max(rl,rlmax)
          write(99,*)i,j,real(rl)
        enddo !j
      enddo !i
      print*, 'Min/max side length (m) =',rlmin,rlmax

      !Bnd
      write(12,*)0,' = open' !# of open bnd
      write(12,*)0 !total # of open bnd  nodes
      write(12,*)4 !# of land bnd
      write(12,*)(nx+3)*2 !total # of land bnd  nodes
      write(12,*)nx+1,0 !1st land bnd
      do i=1,nx+1
        write(12,*)1+(i-1)*(ny+1)
      enddo !i
      write(12,*)2,0 !2nd land bnd, across dateline
      write(12,*)nn-ny
      write(12,*)1
      write(12,*)2,0 !3rd land bnd
      write(12,*)ny+1
      write(12,*)nn
      write(12,*)nx+1,0 !4th bnd
      do i=1,nx+1
        write(12,*)nn-(i-1)*(ny+1)
      enddo !i

      close(12)

 100  format(i10,3(1x,e16.8))
 101  format(i10,1x,i2,1x,3i10)

      stop
      end
  
