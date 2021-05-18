!     Create wwmbnd.gr3 based on SELFE open bnd segments
!     Inputs:
!            hgrid.gr3, gen_wwmbnd.in (see sample below)
!     Output: wwmbnd.gr3

!     ifort -Bstatic -o gen_wwmbnd gen_wwmbnd.f90

!     Sample gen_wwmbnd.in
!     2 !nope - total # of open bnd segments
!     1 0 !seg #; bc. flag
!     2 2 !seg #; bc. flag

      allocatable :: ifl_wwm(:),xnd(:),ynd(:),ibnd(:),nm(:,:)

      open(14,file='hgrid.gr3',status='old')
      open(10,file='gen_wwmbnd.in',status='old')
      open(12,file='wwmbnd.gr3',status='replace')

      read(10,*)nope2
      allocate(ifl_wwm(nope2))
      do k=1,nope2
        read(10,*)j,ifl_wwm(k)
      enddo !k
      close(10)

      read(14,*); read(14,*)ne,np
      allocate(xnd(np),ynd(np),ibnd(np),nm(ne,3))
      ibnd=0 !detault
      do i=1,np
        read(14,*)j,xnd(i),ynd(i),tmp
      enddo !i
      do i=1,ne
        read(14,*)j,k,nm(i,:)
      enddo !i
      read(14,*)nope
      if(nope/=nope2) stop 'nope/=nope2'
      read(14,*)neta
      do k=1,nope
        read(14,*) nond
        do i=1,nond
          read(14,*) iond
          if(iond>np.or.iond<=0) stop 'iond>np'
          ibnd(iond)=ifl_wwm(k)
        enddo
      enddo
      close(14)

!     Output
      write(12,*); write(12,*)ne,np
      do i=1,np
        write(12,*)i,xnd(i),ynd(i),real(ibnd(i))
      enddo !i
      do i=1,ne
        write(12,*)i,3,nm(i,:)
      enddo !i

      stop
      end
