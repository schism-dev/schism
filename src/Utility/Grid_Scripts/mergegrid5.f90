!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Program to merge 2 hybrid grids (triangles and quads) at their common bnd, 
!     determined by a small distance between corresponding bnd nodes from the 2 grids.
!     It won't check negative elem (to accommodate global lon/lat grids).

!     Input files: screen
!
!     Output files: merged.gr3 (with original depths; when depths mismatch at common bnd, grid 1 prevails)

!     ifort -mcmodel=medium -Bstatic -CB -O2 -o mergegrid5 ~/git/schism/src/Utility/UtilLib/schism_geometry.f90 mergegrid5.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      program merge
      use schism_geometry_mod

      implicit real*8(a-h,o-z)
      character*26 grid1,grid2
      character*30 fgrid1,fgrid2,bnd1,bnd2
      allocatable :: x1(:),y1(:),x2(:),y2(:),ibd1(:),ibd2(:),nm1(:,:),nm2(:,:), &
     &i34_1(:),i34_2(:),ired(:),map(:),x(:),y(:),i34(:),iwork(:),dp1(:), &
     &dp2(:),dp(:)
      integer, allocatable :: elnode(:,:),ic3_1(:,:),elside_1(:,:),isidenode_1(:,:),isdel_1(:,:), &
                             &ic3_2(:,:),elside_2(:,:),isidenode_2(:,:),isdel_2(:,:)
      real(kind=8), allocatable :: xcj_1(:),ycj_1(:),xcj_2(:),ycj_2(:)

      print*, 'Input 2 grid names (larger grid first):'
      read(*,'(a26)')grid1
      read(*,'(a26)')grid2
    
      print*, 'Input tolerance distance:'
      read*, eps

!      do len=1,30
!        fgrid1(len:len)=' '
!        fgrid2(len:len)=' '
!        bnd1(len:len)=' '
!        bnd2(len:len)=' '
!      enddo
!
!      icount1=0
!      icount2=0
!      do len=1,26
!        if(grid1(len:len).ne.' ') then
!          icount1=icount1+1
!          fgrid1(icount1:icount1)=grid1(len:len)
!        endif
!        if(grid2(len:len).ne.' ') then
!          icount2=icount2+1
!          fgrid2(icount2:icount2)=grid2(len:len)
!        endif
!      enddo
!
!      bnd1(1:icount1)=fgrid1(1:icount1)
!      bnd2(1:icount2)=fgrid2(1:icount2)
!      fgrid1(icount1+1:icount1+4)='.gr3'
!      fgrid2(icount2+1:icount2+4)='.gr3'
!      bnd1(icount1+1:icount1+4)='.ebn'
!      bnd2(icount2+1:icount2+4)='.ebn'

      open(10,file=trim(adjustl(grid1)),status='old')
      open(12,file=trim(adjustl(grid2)),status='old')
      open(14,file='merged.gr3',status='replace')

      read(10,*)
      read(10,*) ne1,np1
      allocate(x1(np1),y1(np1),dp1(np1),nm1(4,ne1),i34_1(ne1),stat=istat)
      if(istat/=0) stop 'Failed to alloc. (1)'

      do i=1,np1
        read(10,*)j,x1(i),y1(i),dp1(i)
      enddo
      do i=1,ne1
        read(10,*)j,i34_1(i),(nm1(l,i),l=1,i34_1(i))
      enddo !i
      close(10)

      read(12,*)
      read(12,*) ne2,np2
      allocate(x2(np2),y2(np2),dp2(np2),nm2(4,ne2),i34_2(ne2),stat=istat)
      if(istat/=0) stop 'Failed to alloc. (2)'

      do i=1,np2
        read(12,*)j,x2(i),y2(i),dp2(i)
      enddo
      do i=1,ne2
        read(12,*)j,i34_2(i),(nm2(l,i),l=1,i34_2(i))
      enddo
      close(12)

!     Compute geometry for 2 grids
      call compute_nside(np1,ne1,i34_1,nm1,ns1)
      print*, 'ns1=',ns1
      !Allocate side-related arrays
      allocate(ic3_1(4,ne1),elside_1(4,ne1),isdel_1(2,ns1),isidenode_1(2,ns1),xcj_1(ns1),ycj_1(ns1),stat=istat)
      if(istat/=0) stop 'Failed to alloc. side1'
      !Then compute the rest of side related arrays with additional
      !inputs (xnd,ynd) (x,y coordinates of each node)
      call schism_geometry_double(np1,ne1,ns1,x1,y1,i34_1,nm1,ic3_1,elside_1,isdel_1,isidenode_1,xcj_1,ycj_1)

      call compute_nside(np2,ne2,i34_2,nm2,ns2)
      print*, 'ns2=',ns2
      !Allocate side-related arrays
      allocate(ic3_2(4,ne2),elside_2(4,ne2),isdel_2(2,ns2),isidenode_2(2,ns2),xcj_2(ns2),ycj_2(ns2),stat=istat)
      if(istat/=0) stop 'Failed to alloc. side2'
      call schism_geometry_double(np2,ne2,ns2,x2,y2,i34_2,nm2,ic3_2,elside_2,isdel_2,isidenode_2,xcj_2,ycj_2)

!     Compute bnd
      allocate(map(max(np1,np2)),ibd1(np1),ibd2(np2),stat=istat)
      if(istat/=0) stop 'Failed to alloc. bnd'

      map=0
      do i=1,ns1
        if(isdel_1(2,i)==0) then
          map(isidenode_1(1:2,i))=1
        endif
      enddo !i
      nbd1=0 !# of bnd nodes
      do i=1,np1
        if(map(i)==1) then
          nbd1=nbd1+1
          ibd1(nbd1)=i
        endif
      enddo !i

      map=0
      do i=1,ns2
        if(isdel_2(2,i)==0) then
          map(isidenode_2(1:2,i))=1
        endif
      enddo !i
      nbd2=0 !# of bnd nodes
      do i=1,np2
        if(map(i)==1) then
          nbd2=nbd2+1
          ibd2(nbd2)=i
        endif
      enddo !i
      deallocate(map)
      if(nbd1==0.or.nbd2==0) then
        print*, 'No bndy:',nbd1,nbd2
        stop
      endif
      print*, '# of bndy nodes in 2 grids:',nbd1,nbd2

!     Compute bounding box for grid2 to speed up
      xmin2=minval(x2(ibd2(1:nbd2)))
      xmax2=maxval(x2(ibd2(1:nbd2)))
      ymin2=minval(y2(ibd2(1:nbd2)))
      ymax2=maxval(y2(ibd2(1:nbd2)))

!     Add grid 1 first to combined grid
      np=np1+np2
      ne=ne1+ne2
      allocate(x(np),y(np),dp(np),elnode(4,ne),i34(ne),stat=istat)
      if(istat/=0) stop 'Failed to alloc. (5)'

      x(1:np1)=x1(:)
      y(1:np1)=y1(:)
      dp(1:np1)=dp1(:)

      do i=1,ne1
        i34(i)=i34_1(i)
        do l=1,i34_1(i)
          elnode(l,i)=nm1(l,i)
        enddo !l
      enddo !i

!     Add grid2 to grid1
      do i=1,np2
        x(i+np1)=x2(i) 
        y(i+np1)=y2(i) 
        dp(i+np1)=dp2(i) 
      enddo
      do i=1,ne2
        i34(i+ne1)=i34_2(i)
        do j=1,i34_2(i)
          elnode(j,i+ne1)=nm2(j,i)+np1
        enddo
      enddo

!     Done with some arrays 
      deallocate(dp1,nm1,dp2,nm2)
      if(allocated(i34_1)) deallocate(i34_1)
      if(allocated(i34_2)) deallocate(i34_2)


!     Read bnd info
!      read(11,*)
!      read(11,*) nop1
!      nbd1=0
!      do i=1,nop1
!        read(11,*) ndiv1
!        nbd1=nbd1+ndiv1
!        do j=1,ndiv1
!          read(11,*) 
!        enddo !j
!      enddo !i
!
!      read(13,*)
!      read(13,*) nop2
!      nbd2=0
!      do i=1,nop2
!        read(13,*) ndiv2
!        nbd2=nbd2+ndiv2
!        do j=1,ndiv2
!          read(13,*) 
!        enddo !j
!      enddo !i
!      rewind(11)
!      rewind(13)
!
!      allocate(ibd1(nbd1),ibd2(nbd2),stat=istat)
!      if(istat/=0) stop 'Failed to alloc. (3)'
!
!      read(11,*); read(11,*);
!      nbd1=0
!      do i=1,nop1
!        read(11,*) ndiv1
!        do j=1,ndiv1
!          nbd1=nbd1+1
!          read(11,*) ibd1(nbd1)
!        enddo !j
!      enddo !i
!      close(11)
!
!      read(13,*); read(13,*)
!      nbd2=0
!      do i=1,nop2
!        read(13,*) ndiv2
!        do j=1,ndiv2
!          nbd2=nbd2+1
!          read(13,*) ibd2(nbd2)
!        enddo !j
!      enddo !i
!      close(13)

!     Create mapping
      allocate(map(np),stat=istat)
      if(istat/=0) stop 'Failed to alloc. (6)'

      map=0
      do i=1,nbd1
        nd1=ibd1(i)
        if(x1(nd1)<xmin2-eps*2.or.x1(nd1)>xmax2+eps*2.or.y1(nd1)<ymin2-eps*2.or.y1(nd1)>ymax2+eps*2) cycle

        do j=1,nbd2
          nd2=ibd2(j)
          dist=(x1(nd1)-x2(nd2))**2+(y1(nd1)-y2(nd2))**2
          if(dist<eps*eps) then
            if(map(np1+nd2)/=0) then
              write(*,*)'Duplicate red node in grid2',nd2
              stop
            endif
            map(np1+nd2)=nd1
          endif
        enddo !j
      enddo !i

!     Done with some arrays 
      deallocate(x1,y1,x2,y2,ibd1,ibd2)

!     Order all redundant nodes to be removed
      nred=0 !# of red nodes
      do i=np1+1,np
        if(map(i)/=0) then
          nred=nred+1
!          ired(nred)=i
        endif
      enddo !i

      if(nred==0) then
        write(*,*)'2 grids have no common bnd pts'
        stop
      endif
      write(*,*)nred,' common bnd pts found'
      
      allocate(ired(0:nred+1),iwork(np+1),stat=istat)
      if(istat/=0) stop 'Failed to alloc. (4)'

      nred=0 !# of red nodes
      do i=np1+1,np
        if(map(i)/=0) then
          nred=nred+1
          ired(nred)=i
        endif
      enddo !i
      ired(nred+1)=np+1 !for convenience

!     Remove red nodes
      do i=ne1+1,ne
        do j=1,i34(i)
          if(map(elnode(j,i))/=0) elnode(j,i)=map(elnode(j,i))
        enddo !j
      enddo !i

      do i=1,nred
        do k=ired(i)+1,ired(i+1)-1
          x(k-i)=x(k)
          y(k-i)=y(k)
          dp(k-i)=dp(k)
        enddo !k
      enddo

!     Generate a search table to speed up next section
      ired(0)=0
      do i=0,nred
        do k=ired(i)+1,ired(i+1)-1
          iwork(k)=i
        enddo
        if(i/=0) iwork(ired(i))=-1 !Flag
      enddo !i

      do i=ne1+1,ne
        do j=1,i34(i)
          if(iwork(elnode(j,i))==-1) then
            write(*,*)'Fake hanging node',elnode(j,i)
            stop
          endif
          elnode(j,i)=elnode(j,i)-iwork(elnode(j,i))
        enddo !j
      enddo !i
      np=np-nred !no change for ne
 
!     Output and check
!     Re-initialize map to catch hanging nodes
      map=0

      write(14,*)'Merged grid'
      write(14,*) ne,np
      do i=1,np
        write(14,99)i,x(i),y(i),dp(i)
      enddo
      do i=1,ne
        write(14,*)i,i34(i),(elnode(l,i),l=1,i34(i))
        do l=1,i34(i)
          if(elnode(l,i)>np) then
            write(*,*)'Out of bound',i,elnode(l,i)
            stop
          endif
          map(elnode(l,i))=1 !color
        enddo !l
      enddo
      close(14)

      do i=1,np
        if(map(i)==0) write(*,*)'Hanging node',i
      enddo

99    format(i8,2(1x,e23.16),1x,e16.7)

      stop
      end

