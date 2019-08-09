! Generate bnd list for GNOME: nbe(ne,3) (elem ball); bnd(:,4) - 2 end nodes of a bnd seg,
! followed by 2 flags. First flag indicates if the seg is on exterior
! bnd (0) or island (1,2,...). The second flag indicates if it's an open ocean 
! bnd (0: no (inclduing river bnd); 1: yes).
! GNOME only works for pure triangular grids.
!
! ifort -O2 -Bstatic -CB -g -traceback -o gen_GNOME_bnd gen_GNOME_bnd.f90

!   Input: 
!     (1) hgrid.gr3 with bnd info, generated with gredit. No open bnd's
!     on islands; (2) screen input
!   Output: nbe,out; GNOME_bnd_ext.out (exterior bnd) and GNOME_bnd_island.out (island); debug: fort.98
!           Then: cat GNOME_bnd_ext.out GNOME_bnd_island.out > GNOME_bnd.out (i.e. bnd(:,4))
!   Use gen_GNOME_current.m to generate nc input to GNOME

!      implicit none

!      integer, parameter :: debug=1
      integer, parameter :: mnp=500000
      integer, parameter :: mne=1000000
      integer, parameter :: mns=1500000
      integer, parameter :: mnope=10 !max # of open bnd segments
      integer, parameter :: mnond=1000 !max # of open bnd nodes in each segment
      integer, parameter :: mnei=20 !neighbor
      integer, parameter :: nbyte=4
      real, parameter :: small1=1.e-3 !used to check area ratios
  
!     Vertical postion, salinity, and temperature
      integer, dimension(:,:), allocatable :: kbp

      integer :: ier ! allocate error return.

      dimension xnd(mnp),ynd(mnp),nm(mne,4),dp(mnp),i34(mne)
      dimension iest(mnp),ixy(mnp,2),arco(3)
      dimension wild(100),wild2(100,2)
      dimension nne(mnp),ine(mnp,mnei),ic3(mne,4),nx(4,4,3),js(mne,4),is(mns,2),isidenode(mns,2)
      dimension xcj(mns),ycj(mns),nond(mnope),iond(mnope,mnond),isbnd(mnp)
      dimension ie_land(mne),icolor(mnp),list_ext(mnp)

      print*, 'Input open bnd segment # corresponding to open ocean bnd:'
!'
      read*, iocean

!     Read in hgrid and vgrid
      open(13,file='GNOME_bnd_ext.out',status='replace')
      open(15,file='GNOME_bnd_island.out',status='replace')
      open(17,file='nbe.out',status='replace')
      open(14,file='hgrid.gr3',status='old') 
      read(14,*)
      read(14,*)ne,np
      if(np.gt.mnp.or.ne.gt.mne) then
        write(*,*)'Increase mnp/mne'
        stop
      endif
      do i=1,np
        read(14,*)j,xnd(i),ynd(i),dp(i)
      enddo !i
      do i=1,ne
        read(14,*)j,i34(i),(nm(i,l),l=1,i34(i))
        n1=nm(i,1)
        n2=nm(i,2)
        n3=nm(i,3)
        if(i34(i)==3) then
          area=signa(xnd(n1),xnd(n2),xnd(n3),ynd(n1),ynd(n2),ynd(n3))
        else !quad
          n4=nm(i,4)
!         Check convexity
          ar1=signa(xnd(n1),xnd(n2),xnd(n3),ynd(n1),ynd(n2),ynd(n3))
          ar2=signa(xnd(n1),xnd(n3),xnd(n4),ynd(n1),ynd(n3),ynd(n4))
          ar3=signa(xnd(n1),xnd(n2),xnd(n4),ynd(n1),ynd(n2),ynd(n4))
          ar4=signa(xnd(n2),xnd(n3),xnd(n4),ynd(n2),ynd(n3),ynd(n4))
          if(ar1.le.0.or.ar2.le.0.or.ar3.le.0.or.ar4.le.0) then
            write(*,*)'Concave quadrangle',i,ar1,ar2,ar3,ar4
            stop
          endif

          area=ar1+ar2
        endif
        if(area<=0) then
          write(*,*)'Negative area at',i
          stop
        endif
        write(12,*)i,area
      enddo !i=1,ne
      close(12)

!     Open bnds: isbnd()
!     0: internal; >0: open bnd segmt #; -1: exterior land bnd (not on
!     any open bnd); -2: island bnd nodes
      isbnd=0 
      read(14,*) nope
      read(14,*) neta
      ntot=0
      if(nope>mnope) stop 'Increase mnope (2)'
      do k=1,nope
        read(14,*) nond(k)
        if(nond(k)>mnond) stop 'Increase mnond'
        do i=1,nond(k)
          read(14,*) iond(k,i)
          isbnd(iond(k,i))=k
        enddo !i
      enddo !k

      !Init list of exterior bnd nodes
      nlist_ext=0
      if(nond(1)>0) then
        nlist_ext=nond(1) !# of exterior bnd nodes in the list now
        list_ext(1:nlist_ext)=iond(1,1:nlist_ext)
      endif 

!     Land bnds
      read(14,*) nland
      read(14,*) nvel !total #
      nseg=0 !final # of bnd segments (exterior counted as a single bnd)
      do k=1,nland
        read(14,*)nlnd,ifl
        do i=1,nlnd
          read(14,*)ilnd
          icolor(i)=ilnd !temp save

          if(isbnd(ilnd)==0) then
            if(ifl==0) then !exterior
              isbnd(ilnd)=-1
            else !island
              isbnd(ilnd)=-2
            endif
          endif
        enddo !i 

        if(nlist_ext==0.and.ifl==0) then
          nlist_ext=nlnd
          list_ext(1:nlist_ext)=icolor(1:nlnd)     
        endif

        !Output islands
        if(ifl==1) then
          icolor(nlnd+1)=icolor(1) !cyclic
          nseg=nseg+1
          do i=1,nlnd
            write(15,*)icolor(i),icolor(i+1),nseg,0
          enddo !i
        endif !ifl
      enddo !k
      close(14)

      print*, 'List of exterior bnd init:',nlist_ext,list_ext(1:nlist_ext)

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

      do i=1,np
        nne(i)=0
      enddo

      do i=1,ne
        do j=1,i34(i)
          nd=nm(i,j)
          nne(nd)=nne(nd)+1
          if(nne(nd)>mnei) then
            write(*,*)'Too many neighbors',nd
            stop
          endif
          ine(nd,nne(nd))=i
        enddo
      enddo

!     Compute ball info; this won't be affected by re-arrangement below
      do i=1,ne
        do j=1,i34(i)
          ic3(i,j)=0 !index for bnd sides
          nd1=nm(i,nx(i34(i),j,1))
          nd2=nm(i,nx(i34(i),j,2))
          do k=1,nne(nd1)
            ie=ine(nd1,k)
            if(ie/=i.and.(nm(ie,1)==nd2.or.nm(ie,2)==nd2.or.nm(ie,3)==nd2.or.(i34(ie)==4.and.nm(ie,4)==nd2))) ic3(i,j)=ie
          enddo !k
        enddo !j

        write(17,*)i,ic3(i,1:3)
      enddo !i

      ns=0 !# of sides
      do i=1,ne
        do j=1,i34(i)
          nd1=nm(i,nx(i34(i),j,1))
          nd2=nm(i,nx(i34(i),j,2))
          if(ic3(i,j)==0.or.i<ic3(i,j)) then !new sides
            ns=ns+1
            if(ns>mns) then
              write(*,*)'Too many sides'
              stop
            endif
            js(i,j)=ns
            is(ns,1)=i
            isidenode(ns,1)=nd1
            isidenode(ns,2)=nd2
            xcj(ns)=(xnd(nd1)+xnd(nd2))/2
            ycj(ns)=(ynd(nd1)+ynd(nd2))/2

            is(ns,2)=ic3(i,j) !bnd element => bnd side
!           Corresponding side in element ic3(i,j)
            if(ic3(i,j)/=0) then !old internal side
              iel=ic3(i,j)
              index=0
              do k=1,i34(iel)
                if(ic3(iel,k)==i) then
                  index=k
                  exit
                endif
              enddo !k
              if(index==0) then
                write(*,*)'Wrong ball info',i,j
                stop
              endif
              js(iel,index)=ns
            endif !ic3(i,j).ne.0
          endif !ic3(i,j)==0.or.i<ic3(i,j)
        enddo !j=1,i34
      enddo !i=1,ne

      if(ns<ne.or.ns<np) then
        write(*,*)'Weird grid with ns < ne or ns < np',np,ne,ns
        stop
      endif

!     Construct the exterior bnd
      icolor=0 !flag
      icolor(list_ext(1:nlist_ext))=1
      ifront=list_ext(nlist_ext)
      loop1: do
        do m=1,nne(ifront)
          lflound=0 !flag
          ie=ine(ifront,m)
          do j=1,i34(ie) !sides
            isd=js(ie,j)
            if(is(isd,2)/=0) cycle
            if(isidenode(isd,1)/=ifront.and.isidenode(isd,2)/=ifront) cycle
            nd=isidenode(isd,1)+isidenode(isd,2)-ifront
            !write(99,*)ifront,nd,icolor(nd),isbnd(nd)

            if(nd==list_ext(1)) then !closed the loop
              nlist_ext=nlist_ext+1
              if(nlist_ext>np) stop 'overflow(0)'
              list_ext(nlist_ext)=nd
              exit loop1
            endif
            if(icolor(nd)==0.and.(isbnd(nd)>0.or.isbnd(nd)==-1)) then !new front on exterior
              lflound=1
              icolor(nd)=1
              nlist_ext=nlist_ext+1
              if(nlist_ext>np) stop 'overflow(1)'
              list_ext(nlist_ext)=nd
              ifront=nd
              write(98,*)'new front node:',nd
              cycle loop1
            endif
          enddo !j
        enddo !m
        if(lflound==0) then
          write(*,*)'Unable to find next front node:',ifront
          stop
        endif
      enddo loop1

      write(98,*)nlist_ext,' nodes on exterior bnd'
      write(98,*)list_ext(1:nlist_ext)
     
      do i=1,nlist_ext-1
        n1=list_ext(i); n2=list_ext(i+1)
        if(isbnd(n1)==iocean.and.isbnd(n2)==iocean) then
          itmp=1
        else
          itmp=0
        endif
        write(13,*)n1,n2,0,itmp
      enddo !i
       
      stop
      end

      function signa(x1,x2,x3,y1,y2,y3)
!...  Compute signed area formed by pts 1,2,3
!      implicit real*8(a-h,o-z)

      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2

      return
      end

