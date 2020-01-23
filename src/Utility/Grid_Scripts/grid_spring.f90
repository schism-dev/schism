!     Use spring to improve grid quality. All quad nodes are preserved,
!     so are nodes on bnd or with depth>=threshold or specified inside fixed.gr3.
!     Inputs: screen; hgrid.gr3 (in any projection or lon/lat; with b.c. part); fixed.gr3 (optional)
!     Output: hgrid.spring (no depths or bndinfo); fixed_out.gr3 (debug)

!     ifort -O2 -CB -g -traceback -o grid_spring grid_spring.f90

      implicit real*8(a-h,o-z)
      integer :: nwild(3)
      integer, allocatable :: i34(:),elnode(:,:),nne(:),indel(:,:),nnp(:), &
     &indnd(:,:),isbnd(:),ifixed(:)
      real*8, allocatable :: xnd(:),ynd(:),dp(:),area(:)

      print*, 'Input # of spring iteration wanted:'
      read*, niter
      print*, 'Input min area allowed:'
      read*, area_min_allowed
      if(area_min_allowed<0) stop 'area_min_allowed<0'
      print*, 'Input depth (m) threshold (e.g. channel):'
      read*, h_chan
      print*, 'Do you want (1) or not (0) to specify fixed.gr3?'
      read*, ifixed_gr3

      open(14,file='hgrid.gr3',status='old')
      read(14,*)
      read(14,*) ne,np
      allocate(xnd(np),ynd(np),dp(np),area(ne),i34(ne),elnode(4,ne),nne(np), &
     &nnp(np),isbnd(np),ifixed(np))

      ifixed=0 !init
      do i=1,np
        read(14,*) j,xnd(i),ynd(i),dp(i)
        if(dp(i)>=h_chan) ifixed(i)=1
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

        if(i34(i)==4) ifixed(elnode(1:i34(i),i))=1
      enddo !i=1,ne      

!     Bnd
      isbnd=0
      read(14,*) nope
      read(14,*) neta
      ntot=0
      do k=1,nope
        read(14,*) nond
        do i=1,nond
          read(14,*) nd !iond(k,i)
          isbnd(nd)=k
        enddo !i
      enddo !k

!     Land bnds
      read(14,*) nland
      read(14,*) nvel !total #
      do k=1,nland
        read(14,*)nlnd
        do i=1,nlnd
          read(14,*)ilnd
          if(isbnd(ilnd)==0) isbnd(ilnd)=-1
        enddo !i
      enddo !k
      close(14)

      do i=1,np
        if(isbnd(i)/=0) ifixed(i)=1
      enddo !i

!     fixed.gr3
      if(ifixed_gr3/=0) then
        open(12,file='fixed.gr3',status='old')
        read(12,*); read(12,*)
        do i=1,np
          read(12,*)j,tmp,tmp,tmp2
          if(nint(tmp2)>0) ifixed(i)=1
        enddo !i
        close(12)
      endif !ifixed_gr3/=0

      !Output final fixed grid
      open(12,file='fixed_out.gr3',status='replace')
      write(12,*)'Final fixed nodes'
      write(12,*)ne,np
      do i=1,np
        write(12,'(i12,2(1x,e22.14),1x,i4)')i,xnd(i),ynd(i),ifixed(i)
      enddo !i
      do i=1,ne
        write(12,*)i,i34(i),elnode(1:i34(i),i)
      enddo !i
      close(12)

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
      nnp=0 !estimate 1st
      do i=1,ne
        do j=1,i34(i)
          nd=elnode(j,i)
          nne(nd)=nne(nd)+1
          if(nne(nd)>mnei) stop 'impossible'
          indel(nne(nd),nd)=i
          nnp(nd)=nnp(nd)+i34(i)-1
        enddo
      enddo !i
      mnei_p=maxval(nnp)
      allocate(indnd(mnei_p,np),stat=istat)
      if(istat/=0) stop 'failed to alloc'
      
!     Node-node
      nnp=0
      do i=1,np
        do j=1,nne(i)
          ie=indel(j,i)
          do m=1,i34(ie)
            nd=elnode(m,ie)
            if(nd==i) cycle
            iexist=0 !flag
            do k=1,nnp(i)
              if(nd==indnd(k,i)) then
                iexist=1; exit
              endif
            enddo !k 
            if(iexist==0) then
               nnp(i)=nnp(i)+1
               if(nnp(i)>mnei_p) then
                 write(*,*)'Too many ngbr nodes:',i,mnei_p
                 stop
               endif
               indnd(nnp(i),i)=nd
            endif !iexist
          enddo !m
        enddo !j

        !Debug
!        write(99,*)'Node ngbr:',i,nnp(i),indnd(1:nnp(i),i)
      enddo !i=1,np

!     Spring
      do it=1,niter
        do i=1,np
          if(ifixed(i)==1) cycle

          xold=xnd(i); yold=ynd(i) !save
          xnd(i)=sum(xnd(indnd(1:nnp(i),i)))/nnp(i) 
          ynd(i)=sum(ynd(indnd(1:nnp(i),i)))/nnp(i)
          !Check if valid
          iabort=0
          do j=1,nne(i)
            ie=indel(j,i)
            
            n1=elnode(1,ie)
            n2=elnode(2,ie)
            n3=elnode(3,ie)
            tmp=signa(xnd(n1),xnd(n2),xnd(n3),ynd(n1),ynd(n2),ynd(n3))
            if(tmp<=area_min_allowed) then
              iabort=1; exit
            endif

            if(i34(ie)==4) then
              n4=elnode(4,ie)
              tmp=signa(xnd(n1),xnd(n3),xnd(n4),ynd(n1),ynd(n3),ynd(n4))
              if(tmp<=area_min_allowed) then
                iabort=1; exit
              endif
            endif
          enddo !j
          if(iabort==1) then !restore original position
            xnd(i)=xold; ynd(i)=yold
          endif
        enddo !i=1,np
      enddo !it=1,niter

!     Output
      open(12,file='hgrid.spring',status='replace')
      write(12,*)'After spring:',niter
      write(12,*)ne,np
      do i=1,np
        write(12,'(i12,2(1x,e22.14),1x,f13.3)')i,xnd(i),ynd(i),dp(i) !depth not reliable
      enddo !i
      do i=1,ne
        write(12,*)i,i34(i),elnode(1:i34(i),i)
      enddo !i
      close(12)

      stop
      end

      function signa(x1,x2,x3,y1,y2,y3)
!...  Compute signed area formed by pts 1,2,3
      implicit real*8(a-h,o-z)

      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2

      return
      end

