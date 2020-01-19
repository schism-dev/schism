!     Compute geometric arrays used in SCHISM (side info etc)
!     In the driver routine, the following 2 routines should be called
!     in the following fashion (all float arguments are real*4):
!
!     First get # of sides 'ns' with inputs: np (# of nodes),ne (# of elements), i34, and elnode (connectivity table)
!     call compute_nside(np,ne,i34,elnode,ns) 
!     Allocate side-related arrays
!     allocate(ic3(4,ne),elside(4,ne),isdel(2,ns),isidenode(2,ns),xcj(ns),ycj(ns))
!     Then compute the rest of side related arrays with additional inputs (xnd,ynd) (x,y coordinates of each node)
!     call schism_geometry_[single,double](np,ne,ns,xnd,ynd,i34,elnode,ic3,elside,isdel,isidenode,xcj,ycj)
!     The outputs:
!                ns: # of sides
!                ic3(4,ne): neighboring elements of an element
!                elside(4,ne): 4 sides of an element
!                isdel(2,ns): 2 adjacent elements of a side
!                isidenode(2,ns): 2 end nodes of a side
!                xcj(ns),ycj(ns): x,y of each side center

    module schism_geometry_mod
!#ifdef USE_DOUBLE
!      integer,parameter,private :: RKIND=8
!#else
!      integer,parameter,private :: RKIND=4
!#endif

      integer,save,private :: nx(4,4,3)
      integer,save,allocatable,private :: nne(:),indel(:,:),ic3(:,:)
      public

      contains

      subroutine compute_nside(np,ne,i34,elnode,ns)
      !This routine deals only with int; no real
      implicit real(8)(a-h,o-z),integer(i-n)
      integer, intent(in) :: np,ne,i34(ne),elnode(4,ne)
      integer, intent(out) :: ns !,ic3(3,ne)

      if(allocated(nne)) deallocate(nne)
      if(allocated(ic3)) deallocate(ic3)
      if(allocated(indel)) deallocate(indel)

      allocate(nne(np),ic3(4,ne),stat=istat)
      if(istat/=0) stop 'Failed to alloc. nne'

!     Compute geometry
      do k=3,4
        do i=1,k
          do j=1,k-1
            nx(k,i,j)=i+j
            if(nx(k,i,j)>k) nx(k,i,j)=nx(k,i,j)-k
            if(nx(k,i,j)<1.or.nx(k,i,j)>k) then
              write(*,*)'nx wrong',i,j,k,nx(k,i,j)
              stop
            endif
          enddo !j
        enddo !i
      enddo !k

      nne=0
      do i=1,ne
        do j=1,i34(i)
          nd=elnode(j,i)
          nne(nd)=nne(nd)+1
!          indel(nne(nd),nd)=i
        enddo
      enddo
      mnei=maxval(nne)

      allocate(indel(mnei,np),stat=istat)
      if(istat/=0) stop 'Failed to alloc. indel'
      nne=0
      do i=1,ne
        do j=1,i34(i)
          nd=elnode(j,i)
          nne(nd)=nne(nd)+1
          if(nne(nd)>mnei) then
            write(*,*)'Too many neighbors',nd
            stop
          endif
          indel(nne(nd),nd)=i
        enddo
      enddo !i

!     Compute ball info; this won't be affected by re-arrangement below
      do i=1,ne
        do j=1,i34(i) 
          ic3(j,i)=0 !index for bnd sides
          nd1=elnode(nx(i34(i),j,1),i)
          nd2=elnode(nx(i34(i),j,2),i)
          do k=1,nne(nd1)
            ie=indel(k,nd1)
            if(ie/=i.and.(elnode(1,ie)==nd2.or.elnode(2,ie)==nd2.or.elnode(3,ie)==nd2.or. &
     &(i34(ie)==4.and.elnode(4,ie)==nd2))) ic3(j,i)=ie
          enddo !k
        enddo !j
      enddo !i

      ns=0
      do ie=1,ne
        do j=1,i34(ie) !visit each side associated with element ie
          if(ic3(j,ie)==0.or.ie<ic3(j,ie)) then !new side
            ns=ns+1
          endif
        enddo !j
      enddo !ie

      end subroutine compute_nside

      !Double precision version
      subroutine schism_geometry_double(np,ne,ns0,xnd,ynd,i34,elnode,ic3_out,&
     &elside,isdel,isidenode,xcj,ycj)
      implicit real(8)(a-h,o-z),integer(i-n)
      integer, intent(in) :: np,ne,ns0,i34(ne),elnode(4,ne)
      real(8), intent(in) :: xnd(np),ynd(np)
      integer, intent(out) :: ic3_out(4,ne),elside(4,ne),isdel(2,ns0),isidenode(2,ns0)
      real(8), intent(out) :: xcj(ns0),ycj(ns0)
      
      include 'schism_geometry.txt'

      end subroutine schism_geometry_double

      !Single precision version
      subroutine schism_geometry_single(np,ne,ns0,xnd,ynd,i34,elnode,ic3_out,&
     &elside,isdel,isidenode,xcj,ycj)
      implicit real(4)(a-h,o-z),integer(i-n)
      integer, intent(in) :: np,ne,ns0,i34(ne),elnode(4,ne)
      real(4), intent(in) :: xnd(np),ynd(np)
      integer, intent(out) :: ic3_out(4,ne),elside(4,ne),isdel(2,ns0),isidenode(2,ns0)
      real(4), intent(out) :: xcj(ns0),ycj(ns0)
      
      include 'schism_geometry.txt'

      end subroutine schism_geometry_single

    end module schism_geometry_mod
