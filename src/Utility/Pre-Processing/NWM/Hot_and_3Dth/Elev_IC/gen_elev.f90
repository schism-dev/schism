!     This script was initially intended for
!     making a simple source_sink.in at all elem., and also .th
!     
!     ifort -O2 -mcmodel=medium -CB -Bstatic -o gen_elev gen_elev.f90
      implicit real*8(a-h,o-z)
      allocatable :: x(:),y(:),dp(:),area(:),vso(:),r_rough(:),nlayers(:)
      integer, allocatable :: i34(:),elnode(:,:),i_rain(:),isource(:),iDB(:)
      real(8) :: slope(4)

!     Rain rate
      rain_rate= 0.001 !m/hour
      h0=1e-6

      open(10,file='include.gr3',status='old')
      read(10,*); read(10,*)

      open(14,file='hgrid.gr3')
      read(14,*)
      read(14,*)ne,np
      allocate(x(np),y(np),dp(np),r_rough(np),nlayers(np),iDB(np),i_rain(np),isource(ne),i34(ne),elnode(4,ne),area(ne),vso(ne))
      do i=1,np
        read(14,*)j,x(i),y(i),dp(i)
        read(10,*)dummy,dummy,dummy,iDB(i)
      enddo !i
      close(10)
      vol0=0 !inti volume
      do i=1,ne
        read(14,*)j,i34(i),elnode(1:i34(i),i)
        n1=elnode(1,i)
        n2=elnode(2,i)
        n3=elnode(3,i)
        area(i)=signa(x(n1),x(n2),x(n3),y(n1),y(n2),y(n3))
        if(area(i)<=0) then
          print*, 'area<=0 (0):',i,area(i),i34(i)
          stop
        endif

        if(i34(i)==3) then
        else if(i34(i)==4) then
          n4=elnode(4,i)
          area(i)=area(i)+signa(x(n1),x(n3),x(n4),y(n1),y(n3),y(n4))
          if(area(i)<=0) then
            print*, 'area<=0:',i,area(i)
            stop
          endif
        else
          stop 'Unknown elem'
        endif

        tmp=minval(dp(elnode(1:i34(i),i)))
        if(tmp>0) vol0=vol0+area(i)*sum(dp(elnode(1:i34(i),i)))/i34(i)
      enddo !i
      close(14)

      open(15,file='vgrid.in')
      read(15,*) ivcor
      read(15,*) nvrt
      if (ivcor==2) then !s-layer
        nlayers=nvrt-1
      elseif (ivcor==1) then
        do i=1,np
          read(15,*) j,kbp
          nlayers(i)=nvrt-kbp
        enddo
      endif
      close(15)

      total_area=sum(area)
      print*, 'Total area=',total_area !sum(area)
      print*, 'Init volume=',vol0

      open(8,file='rough.gr3',status='replace')
      open(9,file='elev.ic',status='replace')
      write(8,*); write(8,*)ne,np
      write(9,*); write(9,*)ne,np


      r_rough_background=2e-4
      r_rough_land=0.01
      slope_thres=999999.; depth2=-3.0; depth1=-1.0
      !init
      r_rough = r_rough_background

      !based on bathymetry
      do i=1,np
        if (iDB(i) ==0 ) cycle

        if (y(i) > 4452200 ) then !upstream river, y value in utm espg:26918
          r_rough(i)=r_rough_land
        else !DB
          if(dp(i)<depth2) then !upland
            r_rough(i)=r_rough_land
          elseif(dp(i)<=depth1) then !transitional zone
            r_rough(i)=(  (r_rough_land-r_rough_background) * abs(depth1-dp(i))/(depth1-depth2) + r_rough_background ) 
          endif
        endif ! DE R. or DB
        if (r_rough(i) < -1000) then
          print *, i,dp(i),r_rough(i)
        endif
      enddo

      !multi-layers use small r_rough
      do i=1,ne
        if (maxval(nlayers(elnode(1:i34(i),i))) >1 ) then
          do j=1,i34(i)
            n1=elnode(j,i)
            r_rough(n1) = r_rough_background
          enddo
        endif
      enddo

      !!based on slope and number of layers
      !do i=1,ne
      !  !low land uses small r_rough
      !  if (maxval(dp(elnode(1:i34(i),i))) >= depth_thres) cycle
      !  !multi-layers use small r_rough
      !  if (maxval(nlayers(elnode(1:i34(i),i))) >1 ) cycle

      !  slope = 0.0
      !  do j=1,i34(i)
      !    n1=elnode(j,i);
      !    if ((j+1) <= i34(i)) then 
      !      n2=elnode(j+1,i)
      !    else
      !      n2=elnode(1,i)
      !    endif
      !    distance = sqrt((x(n1)-x(n2))**2+(y(n1)-y(n2))**2)
      !    slope(j)=abs((dp(n1)-dp(n2))/distance)
      !  enddo
      !  if (maxval(slope) > slope_thres) then
      !    do j=1,i34(i)-1
      !      n1=elnode(j,i)
      !      r_rough(n1)=-r_CD_manning !negative as drag  
      !    enddo
      !  else
      !    do j=1,i34(i)-1
      !      n1=elnode(j,i)
      !      r_rough(n1)=-r_CD_manning !negative as drag
      !      r_rough(n1)=-(  (-r_rough(n1)-drag_min) * maxval(slope)/slope_thres + drag_min ) !linearly varying, slope=0: drag_min; slope=slope_thres: largest drag
      !    enddo
      !  endif
      !enddo

      do i=1,np
        write(8,'(i8,3(1x,f15.6))')i,x(i),y(i),r_rough(i)
        if (iDB(i) ==0 ) then
          write(9,'(i8,3(1x,f15.6))')i,x(i),y(i),max(0.d0,-dp(i))
        else
          write(9,'(i8,3(1x,f15.6))')i,x(i),y(i),max(0.d0,-dp(i)-1e-6)
        endif
      enddo 
      do i=1,ne
        write(8,*)i,i34(i),elnode(1:i34(i),i)
        write(9,*)i,i34(i),elnode(1:i34(i),i)
      enddo
      close(8)
      close(9)

      nsource = 0
      do i=1,ne
         !only in selected regions
         if ( i==64871 .or. i==497175 )  then
           nsource = nsource +1
         endif
      enddo !i

      open(10,file='source_sink.in',status='replace')
      write(10,*)nsource
      vso=0.0
      isource =0
      nsource = 0
      do i=1,ne
         !only in selected regions
         if ( i==64871 .or. i==497175 )  then
           write(10,*)i
           nsource = nsource +1
           isource(nsource)=i
           if (i==497175) then
             vso(i)=1500.0
             print *, 'setting Susquehanna'
           elseif (i==64871) then
             vso(i)=400.0
             print *, 'setting Delaware'
           endif
         endif
      enddo !i
      write(10,*)
      write(10,*)0 !# of sinks
      close(10)

      open(10,file='vsource.th',status='replace')
      write(10,'(9000000(1x,e16.5))')0.,vso(isource(1:nsource))
      write(10,'(9000000(1x,e16.5))')8640000.,vso(isource(1:nsource))
      close(10)

      open(12,file='msource.th',status='replace')
      write(12,'(9000000(1x,f3.0))')0.,(24.,j=1,nsource),(0.,j=1,nsource)
      write(12,'(f14.1,1x,9000000(1x,f3.0))')8640000.,(4.,j=1,nsource),(0.,j=1,nsource)
      close(12)

!     Analytical change in volume from init volume (m^3)
!     When comparing to model, remember to remove the init volume in model
      dt=100 !sec
      do it=1,5.*86400/dt
        time=it*dt/86400 !days
        if(time<=0.5) then
          dV=total_area*rain_rate*24*time
        else if(time<=1) then
          dV=total_area*rain_rate*24*(0.75-(time-1)**2)
        else
          dV=total_area*rain_rate*24*0.75
        endif
        write(7,*)time,dV
      enddo !it

      stop
      end

function signa(x1,x2,x3,y1,y2,y3)
!-------------------------------------------------------------------------------
! Compute signed area formed by pts 1,2,3
!-------------------------------------------------------------------------------
  implicit real*8(a-h,o-z)
  real*8,intent(in) :: x1,x2,x3,y1,y2,y3

  signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2d0

end function signa

