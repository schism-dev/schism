!     Generate fg.bp for FES etc; works for mixed tri/quad
!     Input: 
!           hgrid.gr3 (for FES, use hgrid.ll with boundary info); 
!           gen_fg.in: 
!                    1st line: nob - # of open boundary segments that need *3D.th;
!                    2nd line: iob(1:nob) - list of segment IDs in hgrid.gr3
!     Output: fg.bp (list of open bnd nodes along the specified segments)

!     ifort -Bstatic -O3 -o gen_fg gen_fg.f90
!     pgf90 -O2 -mcmodel=medium  -Bstatic -o  gen_fg gen_fg.f90

      implicit real*8(a-h,o-z)
!      integer,allocatable :: elnode(:,:)
      allocatable :: x(:),y(:),dp(:),nond(:),iond(:,:),iob(:)

      open(13,file='gen_fg.in',status='old')
      read(13,*)nob
      allocate(iob(nob),stat=istat)
      if(istat/=0) stop 'Failed to allocate (0)'

      read(13,*)iob(:)
      close(13)

      open(14,file='hgrid.gr3',status='old')
      open(12,file='fg.bp')
      read(14,*)
      read(14,*) ne,np
      allocate(x(np),y(np),dp(np),stat=istat)
      if(istat/=0) stop 'Failed to allocate (1)'     
 
      do i=1,np
        read(14,*) j,x(i),y(i),dp(i)
      enddo !i
      do i=1,ne
        read(14,*) !j,l,(elnode(k,i),k=1,3)
      enddo !ii

!     Open bnds
      read(14,*) nope
      allocate(nond(nope),stat=istat)
      if(istat/=0) stop 'Failed to allocate (2)'

      read(14,*) neta
      ntot=0
      mnond=0
      do k=1,nope
        read(14,*) nond(k)
        if(nond(k)>mnond) mnond=nond(k)
        do i=1,nond(k)
          read(14,*) !iond(k,i)
        enddo
        ntot=ntot+nond(k)
      enddo !k

      if(neta/=ntot) then
        write(11,*)'neta /= total # of open bnd nodes',neta,ntot
        stop
      endif

      allocate(iond(nope,mnond),stat=istat)
      if(istat/=0) stop 'Failed to allocate (3)'
      rewind(14)

      do i=1,2+np+ne+2
        read(14,*)
      enddo !i

      do k=1,nope
        read(14,*) nond(k)
        do i=1,nond(k)
          read(14,*) iond(k,i)
        enddo
        if(iond(k,1)==iond(k,nond(k))) then
          write(11,*)'Looped open bnd:',k
          stop
        endif
      enddo !k
      close(14)

!     Output
      write(12,*)'fg.bp'
      write(12,*)sum(nond(iob(1:nob)))
      do k=1,nob
        ibnd=iob(k)
        if(ibnd>nope) then
          write(11,*)'Boundary segment # exceeds max:',ibnd,nope
          stop
        endif
!        print*, k,ibnd,nond(ibnd)
        do i=1,nond(ibnd)
          nd=iond(ibnd,i)
          write(12,'(i10,3e17.8)')nd,x(nd),y(nd),dp(nd)
        enddo !i
      enddo !k
 
      stop
      end
