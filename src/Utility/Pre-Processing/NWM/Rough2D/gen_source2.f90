!     Generate rough.gr3
!     Inputs: hgrid.gr3 (lon/lat), include.gr3 and vgrid.in, the input file include.gr3 (1
!     inside DB) 
!     Outputs: rough.gr3
!     ifort -O2 -mcmodel=medium -CB -Bstatic -o gen_source2 gen_source2.f90
      implicit real*8(a-h,o-z)
      allocatable :: x0(:),y0(:),x(:),y(:),dp(:),area(:),vso(:),r_rough(:),nlayers(:)
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
      allocate(x0(np),y0(np),x(np),y(np),dp(np),r_rough(np),nlayers(np),iDB(np),i_rain(np),isource(ne),i34(ne),elnode(4,ne),area(ne),vso(ne))
      do i=1,np
        read(14,*)j,x(i),y(i),dp(i)
        read(10,*)dummy,dummy,dummy,tmp
        iDB(i)=nint(tmp)
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
      !if(ivcor/=1) stop 'ivcor/=1'
      read(15,*) nvrt
      !do i=1,np
      !  read(15,*) j,kbp
      !  nlayers(i)=nvrt-kbp
      !enddo
      close(15)

      total_area=sum(area)
      print*, 'Total area=',total_area !sum(area)
      print*, 'Init volume=',vol0

      open(8,file='rough.gr3',status='replace')
      open(9,file='elev.ic',status='replace')
      write(8,*); write(8,*)ne,np
      write(9,*); write(9,*)ne,np

      r_rough_ocean=0.019
      r_rough_bay=0.019
      r_rough_upperbay=0.019
      r_rough_land=0.025

      x_min=-75.21
      y_min=39.0 !start transition

      !rotate
      x0=x; y0=y
      x=(x0-x_min)*sqrt(2.0)+(y0-y_min)*sqrt(2.0)
      y=-(x0-x_min)*sqrt(2.0)+(y0-y_min)*sqrt(2.0)
      y_min=0.0

      y_innerbay=0.09 !max y where r_rough_bay starts (between Ship John and Brandy)
      y_upperbay=0.92 !Phily
      slope_thres=999999.; depth2=-3.0; depth1=-1.0
      !init
      r_rough = r_rough_background

      !based on bathymetry
      do i=1,np
        if (iDB(i) ==0 ) then !ocean
          r_rough(i)=r_rough_ocean
!          write(99,*)i,r_rough(i)
        else !DB
          !Transition channel rough
          if(y(i)<y_upperbay) then
            r_chan=r_rough_ocean+(r_rough_bay-r_rough_ocean)*(y(i)-y_min)/(y_innerbay-y_min)
            r_chan=max(r_rough_bay,min(r_rough_ocean,r_chan))
          else !upper Bay; quick transition
            r_chan=r_rough_bay+(r_rough_upperbay-r_rough_bay)*(y(i)-y_upperbay)/0.015
            r_chan=max(r_rough_bay,min(r_rough_upperbay,r_chan))
          endif

          r_rough(i)=r_chan+(r_rough_land-r_chan)*(depth1-dp(i))/(depth1-depth2)
          r_rough(i)=max(r_chan,min(r_rough_land,r_rough(i)))
        endif
      enddo !i

      !multi-layers use small r_rough
      do i=1,ne
        itmp=minval(iDB(elnode(1:i34(i),i)))
        if (itmp/=0.and.maxval(nlayers(elnode(1:i34(i),i))) >1 ) then
          do j=1,i34(i)
            n1=elnode(j,i)
            if(abs(r_rough(n1)-r_rough_land)<1.d-7) then
!              write(99,*)i,r_rough(i)
              r_rough(n1) = r_rough_bay
            endif
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
        write(8,'(i8,3(1x,f15.6))')i,x0(i),y0(i),r_rough(i)
        if (iDB(i) ==0 ) then
          write(9,'(i8,3(1x,f15.6))')i,x0(i),y0(i),max(0.d0,-dp(i))
        else
          write(9,'(i8,3(1x,f15.6))')i,x0(i),y0(i),max(0.d0,-dp(i))+0.02
        endif
      enddo 
      do i=1,ne
        write(8,*)i,i34(i),elnode(1:i34(i),i)
        write(9,*)i,i34(i),elnode(1:i34(i),i)
      enddo
      close(8)
      close(9)

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

