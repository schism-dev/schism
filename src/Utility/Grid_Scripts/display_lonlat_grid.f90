!     Remove wrap-around elements so gredit can display spherical grid
!     The output grid has hanging nodes
!     Input: screen
!     Output: nowrap_[]

!     ifort -O2 -CB -o display_lonlat_grid display_lonlat_grid.f90

      implicit real*8(a-h,o-z)
      character*80 fn1
      integer :: nwild(3)
      integer, allocatable :: i34(:),elnode(:,:),nne(:),indel(:,:),nnp(:), &
     &indnd(:,:),isbnd(:),ifixed(:),iremove(:)
      real*8, allocatable :: xnd(:),ynd(:),dp(:),area(:)

      print*, 'Input full names of .gr3:'
      read*, fn1
      print*, 'Do u want to shift discontinbuity to prime meridian? (1: yes)'
      read*, ishift
!'
    
      fn1=adjustl(fn1); len1=len_trim(fn1)
      open(14,file=fn1(1:len1),status='old')
      read(14,*)
      read(14,*) ne,np
      allocate(xnd(np),ynd(np),dp(np),area(ne),i34(ne),elnode(4,ne),iremove(ne))

      do i=1,np
        read(14,*)j,xnd(i),ynd(i),dp(i)
        if(ishift/=0.and.xnd(i)<0.d0) xnd(i)=xnd(i)+360
      enddo !i

      iremove=0
      do i=1,ne
        read(14,*) j,i34(i),elnode(1:i34(i),i)

        do j=1,i34(i)
          n1=elnode(j,i)
          j2=j+1
          if(j==i34(i)) j2=j2-i34(i)
          n2=elnode(j2,i)
          if(abs(xnd(n1)-xnd(n2))>180.d0) then
            iremove(i)=1
            exit
          endif
        enddo !j
      enddo !i=1,ne      

      ne_remove=sum(iremove)
      ne2=ne-ne_remove
      print*, ne_remove,' elements removed'

!     Output
      open(12,file='nowrap_'//fn1(1:len1),status='replace')
      write(12,*)'Removed wrap-around:',ishift
      write(12,*)ne2,np
      do i=1,np
        write(12,'(i12,2(1x,e22.14),1x,f13.3)')i,xnd(i),ynd(i),dp(i)
      enddo !i
      icount=0
      do i=1,ne
        if(iremove(i)==0) then
          icount=icount+1
          write(12,*)icount,i34(i),elnode(1:i34(i),i)
        endif
      enddo !i
      if(icount/=ne2) then
        print*, 'Mismatch:',ne2,icount
        stop
      endif
      close(12)

      stop
      end

