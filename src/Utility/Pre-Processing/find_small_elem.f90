! Output all elem's that have area <= a given threshold and all small sides
! Inputs: screen; old.gr3 (any proj)
! Outputs: small_[area, side].bp (add header yourself)
! ifort -Bstatic -O2 -mcmodel=medium -o find_small_elem find_small_elem.f90

  implicit real*8(a-h,o-z)
  allocatable :: xnd(:),ynd(:),area(:),ielnode(:,:),i34(:)
   
  print*, 'Input a threshold for area:'
  read*, small1

  print*, 'Input a threshold for side length:'
  read*, small2

  open(14,file='old.gr3',status='old')
  open(12,file='small_area.bp',status='replace')
  open(13,file='small_side.bp',status='replace')
  read(14,*); read(14,*)ne,np
  allocate(xnd(np),ynd(np),area(ne),ielnode(4,ne),i34(ne))
  do i=1,np
    read(14,*)j,xnd(i),ynd(i)
  enddo !i
  icount=0
  icount2=0
  write(12,*)'Threshold=',small1
  write(13,*)'Threshold=',small2
  do i=1,ne
    read(14,*)j,i34(i),ielnode(1:i34(i),i)

    n1=ielnode(1,i); n2=ielnode(2,i); n3=ielnode(3,i);
    ar1=signa(xnd(n1),xnd(n2),xnd(n3),ynd(n1),ynd(n2),ynd(n3))
    if(ar1<=0) then
      print*, 'area<=0:',i,ar1
    endif
    area(i)=ar1
    if(i34(i)==3) then
    else if(i34(i)==4) then
      n4=ielnode(4,i)
      ar2=signa(xnd(n1),xnd(n3),xnd(n4),ynd(n1),ynd(n3),ynd(n4))
      if(ar2<=0) then
        print*, 'area<=0(2):',i,ar2
      endif
      area(i)=area(i)+ar2
    else
      stop 'Unknown elem type'
    endif

    if(area(i)<=small1) then
      icount=icount+1
      xctr=sum(xnd(ielnode(1:i34(i),i)))/i34(i)
      yctr=sum(ynd(ielnode(1:i34(i),i)))/i34(i)
      write(12,'(i12,2(1x,e20.12),1x,e12.5,1x,i12)')icount,xctr,yctr,real(area(i)),i
    endif

    !Side length
    do j=1,i34(i)
      j2=j+1
      if(j2>i34(i)) j2=j2-i34(i)
      nd1=ielnode(j,i)
      nd2=ielnode(j2,i)
      rl=sqrt((xnd(nd1)-xnd(nd2))**2+(ynd(nd1)-ynd(nd2))**2)
      if(rl<=small2) then
        icount2=icount2+1
        xj=(xnd(nd1)+xnd(nd2))/2
        yj=(ynd(nd1)+ynd(nd2))/2
        write(13,'(i12,2(1x,e20.12),1x,e12.5,2(1x,i12))')icount2,xj,yj,rl,nd1,nd2
      endif
    enddo !j
  enddo !i=1,ne
  close(14)
  close(13)
  close(12)
 
  print*, icount,' small elem found; done'
  print*, icount2,' small side found; done'

  stop
  end

      function signa(x1,x2,x3,y1,y2,y3)
!...  Compute signed area formed by pts 1,2,3
      implicit real*8(a-h,o-z)

      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2

      return
      end

