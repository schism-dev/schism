! Analyze global grid (w/o south pole; no land bnd near north pole), and the only exterior
! land bnd is in the Southern Ocean. This is b/cos gredit has large distortions
! in hi latitudes.
! Inputs: node index below; hgrid.ll (tri only)
! Outputs: bnd.out (b.c. part of hgrid.ll)
! ifort -Bstatic -O2 -mcmodel=medium -o gen_boundary_global gen_boundary_global.f90

  implicit real*8(a-h,o-z)
  integer :: nx(4,4,3)
  real*8, allocatable :: xnd(:),ynd(:),dp(:),area(:)
  integer, allocatable :: elnode(:,:),i34(:),nne(:),indel(:,:),ic3(:,:), &
  &icolor(:),icolor2(:),ibnd(:),isdel(:,:),elside(:,:),isidenode(:,:),ilnd(:)
   
  !Input a node on Antarctica bnd:
  node_ant=70665

  open(14,file='hgrid.ll',status='old')
  read(14,*); read(14,*)ne,np
  allocate(xnd(np),ynd(np),dp(np),area(ne),i34(ne),elnode(3,ne),nne(np),ic3(3,ne), &
 &icolor(np),icolor2(np),ibnd(np),ilnd(np))
  do i=1,np
    read(14,*)j,xnd(i),ynd(i),dp(i)
  enddo !i
  i34=3
  do i=1,ne
    read(14,*)j,k,elnode(:,i)

!    n1=elnode(1,i); n2=elnode(2,i); n3=elnode(3,i);
!    ar1=signa(xnd(n1),xnd(n2),xnd(n3),ynd(n1),ynd(n2),ynd(n3))
!    if(ar1<=0) then
!      print*, 'area<=0:',i,ar1
!    endif
!    area(i)=ar1
!
!    if(area(i)<=0) then
!      write(*,*)i,area(i)
!      stop
!    endif
  enddo !i=1,ne
  close(14)
 
! Compute geometry
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
!      indel(nne(nd),nd)=i
    enddo
  enddo
  mnei=maxval(nne)

  allocate(indel(mnei,np),stat=istat)
  if(istat/=0) stop 'Failed to alloc. indel'
  nne=0
  do i=1,ne
    do j=1,3
      nd=elnode(j,i)
      nne(nd)=nne(nd)+1
      if(nne(nd)>mnei) then
        write(*,*)'Too many neighbors',nd
        stop
      endif
      indel(nne(nd),nd)=i
    enddo
  enddo !i

! Compute ball info; this won't be affected by re-arrangement below
  do i=1,ne
    do j=1,3
      ic3(j,i)=0 !index for bnd sides
      nd1=elnode(nx(i34(i),j,1),i)
      nd2=elnode(nx(i34(i),j,2),i)
      do k=1,nne(nd1)
        ie=indel(k,nd1)
        if(ie/=i.and.(elnode(1,ie)==nd2.or.elnode(2,ie)==nd2.or.elnode(3,ie)==nd2)) ic3(j,i)=ie
      enddo !k
    enddo !j
  enddo !i

! Sides
  ns0=0
  do ie=1,ne
    do j=1,3 !visit each side associated with element ie
      nd1=elnode(nx(i34(ie),j,1),ie)
      nd2=elnode(nx(i34(ie),j,2),ie)

      if(ic3(j,ie)==0.or.ie<ic3(j,ie)) ns0=ns0+1
    enddo !j
  enddo !i

  allocate(elside(3,ne),isdel(2,ns0),isidenode(2,ns0))

  ns=0 !# of sides
  do i=1,ne
    do j=1,i34(i)
      nd1=elnode(nx(i34(i),j,1),i)
      nd2=elnode(nx(i34(i),j,2),i)
      if(ic3(j,i)==0.or.i<ic3(j,i)) then !new sides
        ns=ns+1
        if(ns>ns0) then
          write(*,*)'Too many sides'
          stop
        endif
        elside(j,i)=ns
        isdel(1,ns)=i
        isidenode(1,ns)=nd1
        isidenode(2,ns)=nd2
        !xcj(ns)=(xnd(nd1)+xnd(nd2))/2
        !ycj(ns)=(ynd(nd1)+ynd(nd2))/2

        isdel(2,ns)=ic3(j,i) !bnd element => bnd side
!       Corresponding side in element ic3(j,i)
        if(ic3(j,i)/=0) then !old internal side
          iel=ic3(j,i)
          index=0
          do k=1,i34(iel)
            if(ic3(k,iel)==i) then
              index=k
              exit
            endif
          enddo !k
          if(index==0) then
            write(*,*)'Wrong ball info',i,j
            stop
          endif
          elside(index,iel)=ns
        endif !ic3(j,i).ne.0
      endif !ic3(j,i)==0.or.i<ic3(j,i)
    enddo !j
  enddo !i=1,ne

! Mark all bnd nodes
  icolor=0
  ibnd=0
  do ie=1,ne
    do j=1,3 !visit each side associated with element ie
      nd1=elnode(nx(i34(ie),j,1),ie)
      nd2=elnode(nx(i34(ie),j,2),ie)

      if(ic3(j,ie)==0) then
        ibnd(nd1)=1; ibnd(nd2)=1
        !Find bnd near North pole
        ytmp=maxval(ynd(elnode(:,ie)))
        if(ytmp>84) then
          write(*,*)'Plz no land bnd near North pole:',nd1,nd2
          stop 
        endif

        !Antarctica
        ytmp=minval(ynd(elnode(:,ie)))
        if(ytmp<-63.06) then
          icolor(nd1)=1
          icolor(nd2)=1
        endif
      endif !ic3(j,ie)==0
      
    enddo !j
  enddo !ie

! Sort Antarctica bnd
  open(13,file='bnd.tmp',status='replace')
  if(icolor(node_ant)/=1) stop 'Input node is not on southern bnd'
  next=node_ant
  write(13,*)next,0,0 !last 2 values not important (for re-reading only)
  icolor2=0
  icolor2(next)=1 !in the list
  nlnd=1
  loop2: do 
    next0=next !save
    next=-1 !flag
    loop1: do j=1,nne(next0)
      ie=indel(j,next0)
      do m=1,3 !sides
        isd=elside(m,ie)
        if(isdel(2,isd)/=0) cycle
        if(isidenode(1,isd)/=next0.and.isidenode(2,isd)/=next0) cycle

        nd=isidenode(1,isd)+isidenode(2,isd)-next0
        if(nd/=next0.and.icolor(nd)==1) then
          if(nlnd>2.and.nd==node_ant) then !closed the loop (assuming at least 3 nodes)
            exit loop2
          else if(icolor2(nd)==0) then !not in list; new node
            next=nd
            last_node=nd !save last node
            icolor2(nd)=1 !put in the list
            write(13,*)next,0,0
            nlnd=nlnd+1
            exit loop1
          endif !nd
        endif
      enddo !m
    enddo loop1 !j
    if(next<0) then
      write(*,*)'Failed to find next node:',i,next0
      stop
    endif
 
  enddo loop2
  
  print*, nlnd,' land bnd nodes found in the south'

  write(13,*)-1,nlnd,0 !'-1' as a flag for re-reading
  !Close the bnd
  write(13,*)last_node,0,0
  write(13,*)node_ant,0,0
  write(13,*)-1,2,0
  nlnd=nlnd+2

! Deal with islands
  islands=0 !# of islands
  do i=1,np
    if(ibnd(i)==0.or.icolor2(i)==1) cycle

    !bnd node
    write(13,*)i,0,0
    i00=i !save
    next=i
    icolor2(next)=1 !in the list
    nlnd=nlnd+1
    nlnd2=1 !# of nodes on this island
    islands=islands+1
    print*, 'starting node of new island:',i00
    loop4: do 
      next0=next !save
      next=-1 !flag
      loop3: do j=1,nne(next0)
        ie=indel(j,next0)
        do m=1,3 !side
          isd=elside(m,ie)
          if(isdel(2,isd)/=0) cycle
          if(isidenode(1,isd)/=next0.and.isidenode(2,isd)/=next0) cycle

          nd=isidenode(1,isd)+isidenode(2,isd)-next0
          if(nd/=next0.and.ibnd(nd)==1) then
            if(nlnd2>2.and.nd==i00) then !closed the loop (assuming at least 3 nodes)
              exit loop4
            else if(icolor2(nd)==0) then !not in list; new node
              next=nd
              icolor2(nd)=1 !put in the list
              write(13,*)next,0,0
              nlnd2=nlnd2+1
              nlnd=nlnd+1
              exit loop3
            endif !nd
          endif
        enddo !m
      enddo loop3 !j
      if(next<0) then
        write(*,*)'Failed to find next node:',i,next0
        stop
      endif
    enddo loop4
    write(13,*)-1,nlnd2,1
    print*, nlnd2,' nodes on island # ',islands
  enddo !i=1,np
  print*, 'Total # of land nodes=',nlnd,sum(ibnd) !2 nodes repeated on outer bnd

  close(13)
  open(13,file='bnd.tmp',status='old')
  open(10,file='bnd.out',status='replace')
  write(10,*)0, '=open'; write(10,*)0;
  write(10,*)2+islands !number of land boundaries
  write(10,*)nlnd !Total number of open boundary nodes
  nlines=0
  indx=0 !array index
  ilnd=-99 !flag
  do 
    read(13,*,end=100,err=100)i1,i2,i3
    nlines=nlines+1
    if(i1<0) then !reset index and output this seg
      write(10,*)i2,i3
      do i=1,indx
        write(10,*)ilnd(i)
      enddo !i
      ilnd=-99 !flag
      indx=0
    else
      indx=indx+1
      ilnd(indx)=i1
    endif
  enddo
100 print*, nlines,'  lines read from bnd output'
  

  stop
  end

      function signa(x1,x2,x3,y1,y2,y3)
!...  Compute signed area formed by pts 1,2,3
      implicit real*8(a-h,o-z)

      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2

      return
      end

