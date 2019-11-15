!     Given a rectangular cartesian grid, generate a simple rectangular lon/lat grid, and add wrap-around elem. for sphere
!     (a virtual 'earth') so the grid can be used for periodic b.c.
!     After this is done, use ics=2 (and a constant Coriolis if
!     desired), and set rearth_pole and rearth_eq to be r_pole,r_eq below. 
!     The grid actually has no boundary at dateline (\pm 180
!     degrees), so no b.c. is needed there (i.e. a periodic construction).
!     The script basically drapes the original Cartesian grid onto a
!     sphere, with lat=[-0.5,0.5] in order to have minimal distortion to
!     the distance along the meridional

!     Inputs: constants below
!     Outputs: hgrid.ll; debug outputs fort.9[89]
!     ifort -Bstatic -O3 -o periodic_grid periodic_grid.f90
      implicit real*8(a-h,o-z)
      parameter(maxdim=9000000) !for dimensioning
      dimension x(maxdim),y(maxdim),ic1(maxdim,3),ic2(maxdim,20)
      dimension nsur(maxdim),x3(maxdim),y3(maxdim),z3(maxdim)
      dimension x3_0(maxdim),y3_0(maxdim),z3_0(maxdim)

      pi=acos(-1.d0)
      !Lengths in x,y (m)
      rlx=100e3; rly=10e3
      !Resolution in x,y (m)
      !# of divisions in x must >=360
      dx=rlx/361; dy=1e3

      !Min/max lat (degr) - use a small meridional band to minimize
      !distorion
      !The global (x,y,z) are calculated via:
      !x=R_e*cos(\phi)*cos(\lambda)
      !y=R_e*cos(\phi)*sin(\lambda)
      !z=R_p*sin(\phi)
      ymin=-0.5; ymax=0.5 

      !Radius @ pole to match rly
      r_pole=rly/(ymax-ymin)/pi*180
      
      nx=rlx/dx+0.01 !# of divisions in x
      ny=rly/dy+0.01 !# of divisions in y
      if(nx<360) stop 'Make sure # of divisions in x >=360!'

      !Back calculate radius of 'earth' @ equator to match length in x
      r_eq=rlx/2/pi 

      !Resolutions in lon/lat (degr)
      dx_lon=360./(nx+1)
      dy_lat=(ymax-ymin)/ny
      !Min/max lon (degr)
      xmin=-180
      xmax=180-dx_lon !leave a gap 

      print*, 'Radii @ pole and equator=',r_pole,r_eq
      print*, 'nx,ny=',nx,ny
      print*, 'Resolutions in meter=',dx,dy
      print*, 'Resolutions in lon/lat (degr)=',dx_lon,dy_lat

      open(12,file='hgrid.ll',status='replace')
      if((nx+1)*(ny+1).gt.maxdim.or.nx*ny*2.gt.maxdim) then
        write(*,*)'Increase maxdim'
        stop
      endif

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
         !Split into pairs of tri's
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
      write(10,*)
      write(10,*)nt,nn
      write(9,*)
      write(9,*)nt,nn
      do i=1,nn
        write(12,100)i,x(i),y(i),50.

        !global coord. (SCHISM)
        x3(i)=r_eq*cos(y(i)/180*pi)*cos(x(i)/180*pi)
        y3(i)=r_eq*cos(y(i)/180*pi)*sin(x(i)/180*pi)
        z3(i)=r_pole*sin(y(i)/180*pi)

        write(98,*)i,x3(i),y3(i),z3(i)
      enddo !i
      write(12,101)(j,3,(ic1(j,i),i=1,3),j=1,nt)

      !Calc side lengths 
      rl_min=huge(1.d0)
      rl_max=-rl_min
      do i=1,nt
        do j=1,3
          j1=j+1
          if(j1>3) j1=j1-3
          n1=ic1(i,j); n2=ic1(i,j1)
          rl=sqrt((x3(n1)-x3(n2))**2+(y3(n1)-y3(n2))**2+(z3(n1)-z3(n2))**2)
          rl_min=min(rl_min,rl)
          rl_max=max(rl_max,rl)
          write(99,*)i,j,real(rl)
        enddo !j
      enddo !i
      print*, 'Max/min side length=',rl_max,rl_min,sqrt(dx*dx+dy*dy)

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
 101  format(i10,1x,i2,1x,3i9)

      stop
      end
  
