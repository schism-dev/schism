!     Generate tvd.prop for hybrid WENO/ELM based on depth, and remove
!     'isolated' WENO elements etc.
!     Typically, you need to run this program twice: (1) account for
!     eddying regime by setting iaugment=1;
!     (2) account for nearshore regime by setting iaugment=0 (also
!     possibly a different hmin). Then
!     create a nearshore region and adjust flags inside using (2) (e.g.
!     using combine_TVD.f90). You may do additional fine adjustments also.
!     Inputs: screen; 
!             hgrid.gr3 (in any projection or lon/lat; b.c. part not needed)
!     Output: tvd.prop.0

!     ifort -O2 -o gen_tvd_WENO UtilLib/schism_geometry.f90 gen_tvd_WENO.f90

      use schism_geometry_mod
      implicit real*8(a-h,o-z)
      integer :: nwild(3)
      integer, allocatable :: i34(:),elnode(:,:),nne(:),indel(:,:),nnp(:), &
     &indnd(:,:),isbnd(:),itvd(:),itvd0(:)
      integer, allocatable :: ic3(:,:),elside(:,:),isdel(:,:),isidenode(:,:)
      real*8, allocatable :: xnd(:),ynd(:),dp(:),area(:),xpoly(:),ypoly(:),xcj(:),ycj(:)

      print*, 'Input min depth in meters:'
      read*, hmin

      print*, 'Augment upwind by 1 extra layer? (1: yes); 0: no'
      read*, iaugment

      open(14,file='hgrid.gr3',status='old')
      read(14,*)
      read(14,*) ne,np
      allocate(xnd(np),ynd(np),dp(np),area(ne),i34(ne),elnode(4,ne),nne(np), &
     &isbnd(np),itvd(ne),itvd0(ne))

      do i=1,np
        read(14,*) j,xnd(i),ynd(i),dp(i)
      enddo !i

      do i=1,ne
        read(14,*) j,i34(i),elnode(1:i34(i),i)

        n1=elnode(1,i)
        n2=elnode(2,i)
        n3=elnode(3,i)
        area(i)=signa(xnd(n1),xnd(n2),xnd(n3),ynd(n1),ynd(n2),ynd(n3))
        if(area(i)<=0) then
          write(*,*)'Negative area at elem:',i
          stop
        endif
         
        if(i34(i)==4) then
          n4=elnode(4,i)
          tmp=signa(xnd(n1),xnd(n3),xnd(n4),ynd(n1),ynd(n3),ynd(n4))
          if(tmp<=0) then
            write(*,*)'Negative area at elem:',i
            stop
          endif
          area(i)=area(i)+tmp
        endif
      enddo !i=1,ne      
      close(14)

!     Neighborhood
      nne=0
      do i=1,ne
        do j=1,i34(i)
          nd=elnode(j,i)
          nne(nd)=nne(nd)+1
        enddo
      enddo !i
      mnei=maxval(nne)

      allocate(indel(mnei,np))
      nne=0
!      nnp=0 !estimate 1st
      do i=1,ne
        do j=1,i34(i)
          nd=elnode(j,i)
          nne(nd)=nne(nd)+1
          if(nne(nd)>mnei) stop 'impossible'
          indel(nne(nd),nd)=i
!          nnp(nd)=nnp(nd)+i34(i)-1
        enddo
      enddo !i

      call compute_nside(np,ne,i34,elnode,ns)
!     Allocate side-related arrays
      allocate(ic3(4,ne),elside(4,ne),isdel(2,ns),isidenode(2,ns),xcj(ns),ycj(ns))
!     Then compute the rest of side related arrays with additional
!     inputs (xnd,ynd) (x,y coordinates of each node)
      call schism_geometry_double(np,ne,ns,xnd,ynd,i34,elnode,ic3,elside,isdel,isidenode,xcj,ycj)
      
!     Set TVD flag
      itvd(:)=0 !init
      do i=1,ne
        hmin2=minval(dp(elnode(1:i34(i),i)))
!        xctr=sum(xnd(elnode(1:i34(i),i)))/i34(i)
!        yctr=sum(ynd(elnode(1:i34(i),i)))/i34(i)
        if(hmin2>=hmin) then
          itvd(i)=1
        endif 
      enddo !i

!     Augment upwind zone by 1 extra layer
      if(iaugment/=0) then
        itvd0=itvd
        do i=1,ne
          if(itvd0(i)==0) then
            do j=1,i34(i)
              nd=elnode(j,i)
              do m=1,nne(nd)
                ie=indel(m,nd)
                itvd(ie)=0
              enddo !m
            enddo !j
          endif !in_out
        enddo !i
      endif !iaugment/

!     Remove 'isolated' WENO elem
      loop1: do
        itouched=0 !# of elem changed in this iteration
        itvd0=itvd
        loop2: do i=1,ne
          if(itvd0(i)==0) cycle

          icount=0
          do j=1,i34(i)
            ie=ic3(j,i)
            if(ie/=0) then 
              if(itvd0(ie)/=0) icount=icount+1
            endif !ie
          enddo !j 

          if(icount<=1) then
            itvd(i)=0
            itouched=itouched+1
          endif
        end do loop2 !i=1,ne

        print*, '# of elem flipped=',itouched
        if(itouched==0) exit loop1
      end do loop1

!     Output
      open(12,file='tvd.prop.0',status='replace')
      do i=1,ne
        write(12,*)i,itvd(i)
      enddo !i
      write(12,*)'hmin=',hmin
      close(12)

      stop
      end

      function signa(x1,x2,x3,y1,y2,y3)
!...  Compute signed area formed by pts 1,2,3
      implicit real*8(a-h,o-z)

      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2

      return
      end

