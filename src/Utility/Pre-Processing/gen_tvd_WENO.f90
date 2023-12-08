!     Generate tvd.prop for hybrid WENO/ELM and cross-scale applications, based on depth and remove
!     'isolated' WENO elements in shallow water, etc. The goal is to
!     minimize the direct 'contact' between WENO and ELM cells, which
!     can cause dispersion.
!     Inputs: (1) screen; 
!             (2) hgrid.gr3 (in any projection or lon/lat; b.c. part not needed)
!             (3) nearshore.gr3: '0' for eddying regime; non-0 for nearshore. Nearshore region 
!                                 should be slightly offshore of isobath 'hmin_of'
!     Output: tvd.prop.0 (may be further edited, e.g. deep channels/fjords)

!     ifort -O2 -o gen_tvd_WENO ../UtilLib/schism_geometry.f90 gen_tvd_WENO.f90

      use schism_geometry_mod
      implicit real*8(a-h,o-z)
      integer :: nwild(3)
      integer, allocatable :: i34(:),elnode(:,:),nne(:),indel(:,:),nnp(:), &
     &indnd(:,:),isbnd(:),itvd(:),itvd0(:),iest(:),iest_e(:)
      integer, allocatable :: ic3(:,:),elside(:,:),isdel(:,:),isidenode(:,:)
      real*8, allocatable :: xnd(:),ynd(:),dp(:),area(:),xpoly(:),ypoly(:),xcj(:),ycj(:)

      print*, 'Input cut-off depth (meters) for offshore (e.g., 30m));'
      print*, 'Upwind will be used shallower than this depth in offshore:'
      read*, hmin_of

      print*, 'Input cut-off depth (meters) for nearshore.gr3 (usually <hmin_of; e.g. same as h_tvd).'
      print*, 'Upwind will be used shallower than this depth in nearshore region:'
      read*, hmin_near
!'

      open(14,file='hgrid.gr3',status='old')
      open(13,file='nearshore.gr3',status='old')
      read(14,*); read(14,*) ne,np
      read(13,*); read(13,*) 
      allocate(xnd(np),ynd(np),dp(np),area(ne),i34(ne),elnode(4,ne),nne(np), &
     &isbnd(np),itvd(ne),itvd0(ne),iest(np),iest_e(ne))

      do i=1,np
        read(14,*) j,xnd(i),ynd(i),dp(i)
        read(13,*) j,tmp,tmp,tmp2
        iest(i)=nint(tmp2)
      enddo !i

      do i=1,ne
        read(14,*) j,i34(i),elnode(1:i34(i),i)

        iest_e(i)=maxval(iest(elnode(1:i34(i),i)))

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
        if(iest_e(i)==0) then !offshore
          if(hmin2>=hmin_of) itvd(i)=1
        else !nearshore
          if(hmin2>=hmin_near) itvd(i)=1
        endif !iest_e
      enddo !i

!     Augment upwind zone by 1 extra layer for offshore
      itvd0=itvd
      do i=1,ne
        if(iest_e(i)==0.and.itvd0(i)==0) then
          do j=1,i34(i)
            nd=elnode(j,i)
            do m=1,nne(nd)
              ie=indel(m,nd)
              itvd(ie)=0
            enddo !m
          enddo !j
        endif !ifl
      enddo !i

!     Remove 'isolated' WENO elem nearshore
      loop1: do
        itouched=0 !# of elem changed in this iteration
        itvd0=itvd
        loop2: do i=1,ne
          if(iest_e(i)==0.or.itvd0(i)==0) cycle

          !WENO elem nearshore 
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
      write(12,*)'hmin=',hmin_of,hmin_near
      close(12)

      stop
      end

      function signa(x1,x2,x3,y1,y2,y3)
!...  Compute signed area formed by pts 1,2,3
      implicit real*8(a-h,o-z)

      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2

      return
      end

