! Simple prep script for METIS

! Inputs: hgrid.gr3 (with b.c.); vgrid.in
! Outputs: graphinfo (required by METIS)

! ifort -mcmodel=medium  -O2 -CB -g -traceback -o metis_prep metis_prep.f90 

! METIS usage:
! ./gpmetis graphinfo <nproc> -ufactor=1.01 -seed=15
! Afterward, use awk to generate partition.prop:
! awk '{print NR,$0}' graphinfo.part.88 > partition.prop      (using nproc=88 cores)
! and add 2 extra lines at the end of partition.prop:
! nproc
! <nproc value>

  program metis_prep
  implicit none
!  include 'metis.h'

  integer,parameter :: rkind = 8
  logical :: found
  integer :: i,j,k,l,ie,je,ne,np,ip,ip0,mnei,n1,n2,nope,neta,mnond,nt,stat,nn,icount,ii, &
 &new,ivcor,nvrt,kz,nsig,kin,mxnedge,nedge,ntedge,wgtflag,numflag,ncon,ndims,kbetmp, &
 &nproc,nxq(3,4,4),nland,mnlnd,nvel,objval
!  integer :: options(METIS_NOPTIONS)
  real(rkind) :: h_c,h_s,theta_f,theta_b,dtmp,etmp,ptmp,stmp

  integer, allocatable :: i34(:),elnode(:,:),nne(:),indel(:,:),ic3(:,:), &
 &isbnd(:),nond(:),iond(:,:),kbp(:),adjncy(:),xadj(:),nlev(:),vwgt(:),adjwgt(:), & 
 &vtxdist(:),nlnd(:),ilnd(:,:),part(:)
  real(rkind), allocatable :: ztot(:),sigma(:),dp(:)
  real(4),allocatable :: tpwgts(:),ubvec(:)
  

!  print*, 'Input # of MPI processes (compute):'
!  read*, nproc

  do k=3,4 !elem. type
    do i=1,k  !local index
      do j=1,k-1 !offset
        nxq(j,i,k)=i+j
        if(nxq(j,i,k)>k) nxq(j,i,k)=nxq(j,i,k)-k
        if(nxq(j,i,k)<1.or.nxq(j,i,k)>k) then
          write(*,*)'INIT: nx wrong',i,j,k,nxq(j,i,k)
          stop
        endif
      enddo !j
    enddo !i
  enddo !k

  !-----------------------------------------------------------------------------
  ! Aquire and construct global data tables
  !-----------------------------------------------------------------------------

  ! Aquire global grid size
  open(14,file='hgrid.gr3',status='old')
  read(14,*); read(14,*) ne,np

  ! Aquire global element-node tables from hgrid.gr3
  allocate(i34(ne),elnode(4,ne),dp(np),kbp(np),stat=stat)

  do i=1,np 
    read(14,*)j,dtmp,dtmp,dp(i)
  enddo;
  do i=1,ne
    read(14,*) ie,i34(ie),(elnode(k,ie),k=1,i34(ie))
    if(i34(ie)/=3.and.i34(ie)/=4) then
      write(*,*) 'AQUIRE_HGRID: Unknown type of element',ie,i34(ie)
      stop
    endif
  enddo !i

  ! Count number of elements connected to each node (global)
  allocate(nne(np),stat=stat)
  nne=0
  do ie=1,ne
    do k=1,i34(ie)
      ip=elnode(k,ie)
      nne(ip)=nne(ip)+1
    enddo !k
  enddo !ie

  ! Check hanging nodes
  found=.false.
  do ip=1,np
    if(nne(ip)==0) then
      found=.true.
      write(*,*) 'Hanging node:',ip
    endif
  enddo
  if(found) then
    write(*,*)'check hanging nodes'
    stop
  endif

  ! Maximum number of elements connected to a node
  mnei=0
  do ip=1,np
    mnei=max(mnei,nne(ip))
  enddo

  ! Build global node-element and self-reference table table
  allocate(indel(mnei,np),stat=stat)
  nne=0
  do ie=1,ne
    do k=1,i34(ie)
      ip=elnode(k,ie)
      nne(ip)=nne(ip)+1
      indel(nne(ip),ip)=ie
    enddo !k
  enddo !ie


  !Compute mnei_p (max. # of nodes around a node)
!  mnei_p=0
!  do i=1,np
!    icount=0
!    do j=1,nne(i)
!      icount=icount+i34(ie)-2 !# of new surrounding nodes
!    enddo !j
!    mnei_p=max(mnei_p,icount+1) !account for bnd ball
!  enddo !i

!  if(myrank==0.and..not.full_aquire) write(16,*)'mnei, mnei_p = ',mnei,mnei_p

  ! Build global element-side-element table; this won't be affected by re-arrangement below
  allocate(ic3(4,ne),stat=stat)
  do ie=1,ne
    do k=1,i34(ie)
      ic3(k,ie)=0 !index for boundary sides
      n1=elnode(nxq(1,k,i34(ie)),ie)
      n2=elnode(nxq(2,k,i34(ie)),ie)
      do l=1,nne(n1)
        je=indel(l,n1)
        if(je/=ie.and.(elnode(1,je)==n2.or.elnode(2,je)==n2.or. &
     &elnode(3,je)==n2.or.(i34(je)==4.and.elnode(4,je)==n2))) ic3(k,ie)=je
      enddo !l
      je=ic3(k,ie)
      if(je/=0) then
        do l=1,i34(je)
          if(elnode(nxq(1,l,i34(je)),je)==n1.and.elnode(nxq(2,l,i34(je)),je)==n2) then
            write(*,*) 'Elem ', ie, ' and ', je, ' have opposite orientation'
            stop
          endif
        end do  !l
      endif
    enddo !k
  enddo !ie

  !new39
!  ip0=718
!  do i=1,nne(ip0)
!    ie=indel(i,ip0)
!    write(99,*)'Init ball:',nne(ip0),ie
!    write(99,*)'ic3:',ic3(1:i34(ie),ie)
!  enddo !i  

  ! ine to be re-arranged in counter-clockwise fashion after boundary info is read in
  ! Count global number of sides and build global element-side index table
!  if(allocated(js)) deallocate(js); allocate(js(4,ne),stat=stat);
!  if(stat/=0) call parallel_abort('AQUIRE_HGRID: js allocation failure')
!  ns=0
!  do ie=1,ne
!    do j=1,i34(ie) !visit each side associated with element ie
!      if(ic3(j,ie)==0.or.ie<ic3(j,ie)) then !new global side
!        ns=ns+1
!        js(j,ie)=ns
!        if(ic3(j,ie)/=0) then !old internal side
!          je=ic3(j,ie)
!          l=0
!          do k=1,i34(je)
!            if(ic3(k,je)==ie) then
!              l=k
!              exit
!            endif
!          enddo !k
!          if(l==0) then
!            write(errmsg,'(a,10i6)') 'AQUIRE_HGRID: Wrong ball info',ie,j,ns
!            call parallel_abort(errmsg)
!          endif
!          js(l,je)=ns
!        endif !ic3(j,ie)/=0
!      endif !ic3(j,ie)==0.or.ie<ic3(j,ie)
!    enddo !j
!  enddo !ie
!  if(ns.lt.ne.or.ns.lt.np) then
!    write(errmsg,*)'AQUIRE_HGRID: weird grid with ns < ne or ns < np', &
!    &np,ne,ns
!    call parallel_abort(errmsg)
!  endif

  !-----------------------------------------------------------------------------
  ! Aquire global open boundary segments from hgrid.gr3
  !-----------------------------------------------------------------------------

  ! Allocate and assign global node-to-open-boundary-segment flags
  allocate(isbnd(np),stat=stat)
  isbnd=0

  ! Global number of open boundary segments and nodes
  rewind(14); read(14,*); read(14,*);
  do i=1,np; read(14,*); enddo;
  do i=1,ne; read(14,*); enddo;
  read(14,*) nope
  read(14,*) neta

  ! Scan segments to count number of open boundary segments and nodes
  mnond=0 !global max number of nodes per segment
  nt=0    !global total node count
  do k=1,nope
    read(14,*) nn
    mnond=max(mnond,nn);
    nt=nt+nn
    do i=1,nn; read(14,*); enddo;
  enddo !k
  if(neta/=nt) then
    write(*,*) 'neta /= total # of open bnd nodes',neta,nt
    stop
  endif

  ! Allocate arrays for global open boundary segments
  allocate(nond(nope),stat=stat)
  allocate(iond(nope,mnond),stat=stat)

  ! Aquire global open boundary segments and nodes
    rewind(14); read(14,*); read(14,*);
    do i=1,np; read(14,*); enddo;
    do i=1,ne; read(14,*); enddo;
    read(14,*); read(14,*);
    nond=0; iond=0;
    do k=1,nope
      read(14,*) nn
      do i=1,nn
        read(14,*) ip
        nond(k)=nond(k)+1
        iond(k,nond(k))=ip
        isbnd(ip)=k
      enddo !i
      if(iond(k,1)==iond(k,nond(k))) then
        write(*,*) 'Looped open bnd:',k
        stop
      endif
    enddo !k

  !-----------------------------------------------------------------------------
  ! Aquire global land boundary segments from hgrid.gr3
  !-----------------------------------------------------------------------------

  ! Global total number of land boundary segments and nodes
    rewind(14); read(14,*); read(14,*);
    do i=1,np; read(14,*); enddo;
    do i=1,ne; read(14,*); enddo;
    read(14,*); read(14,*);
    do k=1,nope; read(14,*) nn; do i=1,nn; read(14,*); enddo; enddo;
    read(14,*) nland
    read(14,*) nvel

    ! Scan segments to count number of land boundary segments and nodes
    mnlnd=0 !global max number of nodes per segment
    nt=0    !global total node count
    do k=1,nland
      read(14,*) nn
      mnlnd=max(mnlnd,nn)
      nt=nt+nn
      do i=1,nn; read(14,*); enddo;
    enddo !k
    if(nvel/=nt) then
      write(*,*) 'AQUIRE_HGRID: nvel /= total # of land bnd nodes', &
                    &nvel,nt
      stop
    endif

  ! Allocate arrays for global land boundary segments
  allocate(nlnd(nland),stat=stat)
  allocate(ilnd(nland,mnlnd),stat=stat)

  ! Aquire global land boundary segments and nodes
    rewind(14); read(14,*); read(14,*);
    do i=1,np; read(14,*); enddo;
    do i=1,ne; read(14,*); enddo;
    read(14,*); read(14,*);
    do k=1,nope; read(14,*) nn; do i=1,nn; read(14,*); enddo; enddo;
    read(14,*); read(14,*);
    nlnd=0; ilnd=0;
    do k=1,nland
      read(14,*) nn
      do i=1,nn
        read(14,*) ip
        nlnd(k)=nlnd(k)+1
        ilnd(k,nlnd(k))=ip
        if(isbnd(ip)==0) isbnd(ip)=-1 !overlap of open bnd
      enddo !i
    enddo !k

    !-----------------------------------------------------------------------------
    ! Done with global grid -- close grid file
    !-----------------------------------------------------------------------------
    close(14)

! Re-arrange in counter-clockwise fashion
  do i=1,np
    if(isbnd(i)/=0) then !bnd ball
!     Look for starting bnd element
      icount=0
      do j=1,nne(i)
        ie=indel(j,i)
        ii=0 !local index
        do l=1,i34(ie)
          if(elnode(l,ie)==i) then
            ii=l; exit
          endif
        enddo !l
        if(ii==0) stop 'AQUIRE_HGRID: bomb (1)'

        if(ic3(nxq(i34(ie)-1,ii,i34(ie)),ie)==0) then
          icount=icount+1
          indel(1,i)=ie
        endif
      enddo !j=1,nne(i)
      if(icount/=1) then
        write(*,*)'Illegal bnd node',i,isbnd(i),icount
        stop
      endif
    endif !bnd ball

!   For internal balls, starting elem. is not altered
!   Sequential search for the rest of elements
!    nnp(i)=2
!    inp(i,1)=elnode(nx(iself(1,i),1),indel(i,1))
!    inp(i,2)=elnode(nx(iself(1,i),2),indel(i,1))
    do j=2,nne(i)
      ie=indel(j-1,i)
      ii=0 !local index
      do l=1,i34(ie)
        if(elnode(l,ie)==i) then
          ii=l; exit
        endif
      enddo !l
      if(ii==0) stop 'AQUIRE_HGRID: bomb (2)'

      new=ic3(nxq(i34(ie)-2,ii,i34(ie)),ie)
      if(new==0) then
        write(*,*)'Incomplete ball',i,j,indel(1:j-1,i)
        stop
      endif
      indel(j,i)=new
    enddo !j=2,nne(i)
  enddo !i=1,np

  !new39
!  do i=1,np
!    write(99,*)'Node ball:',i,nne(i),indel(1:nne(i),i)
!  enddo !i

  !Vgrid
  open(19,file='vgrid.in',status='old')
  read(19,*)ivcor

  if(ivcor==2) then !SZ coordinates
      read(19,*) nvrt,kz,h_s !kz>=1
      if(nvrt<2) stop 'nvrt<2'
      if(kz<1) then !.or.kz>nvrt-2) then
        write(*,*)'Wrong kz:',kz
        stop
      endif
      if(h_s<6.d0) then
        write(*,*)'h_s needs to be larger:',h_s
        stop
      endif

    ! Allocate vertical layers arrays
    allocate(ztot(nvrt),sigma(nvrt),stat=stat)
    nsig=nvrt-kz+1 !# of S levels (including "bottom" & f.s.)

      ! # of z-levels excluding "bottom" at h_s
      read(19,*) !for adding comment "Z levels"
      do k=1,kz-1
        read(19,*)j,ztot(k)
        if(ztot(k)>=-h_s) then
          write(*,*)'Illegal Z level:',k
          stop
        endif
        if(k>1) then; if(ztot(k)<=ztot(k-1)) then
          write(*,*)'z-level inverted:',k
          stop
        endif; endif
      enddo !k
      read(19,*) !level kz       
      ! In case kz=1, there is only 1 ztot(1)=-h_s
      ztot(kz)=-h_s

      read(19,*) !for adding comment "S levels"
      read(19,*)h_c,theta_b,theta_f
      if(h_c<5._rkind.or.h_c>=h_s) then !large h_c to avoid 2nd type abnormaty
        write(*,*)'h_c needs to be larger:',h_c
        stop
      endif
      if(theta_b<0._rkind.or.theta_b>1._rkind) then
        write(*,*)'Wrong theta_b:',theta_b
        stop
      endif
      if(theta_f<=0._rkind) then
        write(*,*)'Wrong theta_f:',theta_f
        stop
      endif

      sigma(1)=-1._rkind !bottom
      sigma(nsig)=0._rkind !surface
      read(19,*) !level kz
      do k=kz+1,nvrt-1
        kin=k-kz+1
        read(19,*) j,sigma(kin)
        if(sigma(kin)<=sigma(kin-1).or.sigma(kin)>=0._rkind) then
          write(*,*)'Check sigma levels at:',k,sigma(kin),sigma(kin-1)
          stop
        endif
      enddo !k
      read(19,*) !level nvrt
      close(19)

  else if(ivcor==1) then !localized sigma
    read(19,*)nvrt 
    read(19,*)kbp(1:np)
    close(19)

    allocate(ztot(nvrt),sigma(nvrt),stat=stat)
    !for output only - remove later
    ztot=0._rkind; sigma=0._rkind 
    kz=1; h_s=0._rkind; h_c=0._rkind; theta_b=0._rkind; theta_f=0._rkind
  else
    stop 'GRID_SUBS: Unknown ivcor'
  endif !ivcor


  ! Count number of edges in dual graph (element as 'vertex')
  allocate(adjncy(1000),stat=stat) !for single element
  ntedge=0 !total # of edges in the grid
  mxnedge=0 !max. # of local edges
  do ie=1,ne
    nedge=0 !# of local edges
    adjncy=0 !list of global element indices
    do j=1,i34(ie)
      ip=elnode(j,ie)
      do k=1,nne(ip)
        je=indel(k,ip)
        if(je/=ie) then
          found=.false.
          do l=1,nedge
            if(adjncy(l)==je) then
              found=.true.
              exit
            endif
          enddo
          if(.not.found) then !new edge
            nedge=nedge+1
            if(nedge>1000) stop 'PARTITION: bound (1)'
            adjncy(nedge)=je
          endif
        endif
      enddo !k
    enddo !j
    ntedge=ntedge+nedge
    mxnedge=max(mxnedge,nedge)
  enddo !ie
  deallocate(adjncy)

  ! Use vertex (elem) and/or edge weights
  wgtflag = 3   ! 0: none; 1: edges; 2: vertices; 3: vertices & edges

  ! Number of weights associated with each vertex, used to optimize the
  ! partition
  ncon = 4

  ! Allocate storage for dual graph
  if(allocated(xadj)) deallocate(xadj); allocate(xadj(ne+1),stat=stat);
  if(allocated(adjncy)) deallocate(adjncy); allocate(adjncy(ntedge),stat=stat)
!  if(allocated(xyz)) deallocate(xyz); allocate(xyz(2*ne),stat=stat)
  if(allocated(nlev)) deallocate(nlev); allocate(nlev(ne),stat=stat) !estimate of number of active levels
  if(wgtflag==2.or.wgtflag==3) then
    if(allocated(vwgt)) deallocate(vwgt); allocate(vwgt(ne*ncon),stat=stat)
  endif
  if(wgtflag==1.or.wgtflag==3) then
    if(allocated(adjwgt)) deallocate(adjwgt); allocate(adjwgt(ntedge),stat=stat)
  endif

  ! Assign vertex coordinates & weights
  ! Weights based on estimated number of active levels
  do ie=1,ne
!    xtmp=0._rkind !xctr
!    ytmp=0._rkind
    dtmp=real(5.d10,rkind) !min. depth
    do j=1,i34(ie)
      ip=elnode(j,ie)
!      xtmp=xtmp+real(xproj(ip),rkind)/real(i34(ie),rkind)
!      ytmp=ytmp+real(yproj(ip),rkind)/real(i34(ie),rkind)
      if(dp(ip)<dtmp) dtmp=dp(ip)
    enddo !j
!    xyz(2*(ie-1)+1)=xtmp
!    xyz(2*(ie-1)+2)=ytmp
    if(wgtflag==2.or.wgtflag==3) then
      if(ivcor==1) then
        kbetmp=minval(kbp(elnode(1:i34(ie),ie)))
      else if(ivcor==2) then !SZ (including 2D)
        if(dtmp<=0._rkind) then
          kbetmp=nvrt !only for estimating nlev
        else 
          if(dtmp<=h_s) then
            kbetmp=kz
          else !>h_s
            kbetmp=0 !element bottom index; also works for 2D model
            do j=1,kz-1
              if(-dtmp>=ztot(j).and.-dtmp<ztot(j+1)) then
                kbetmp=j
                exit
              endif
            enddo !j
          endif
        endif !dtmp
      endif !ivcor

      nlev(ie)=max(1,nvrt-kbetmp) !estimate of number of active levels (excluding bottom as most procedures do not use it)
      etmp=1._rkind; ptmp=0._rkind; stmp=0._rkind
      do j=1,i34(ie)
        ip=elnode(j,ie)
        ptmp=ptmp+1._rkind/real(nne(ip),rkind) !each node contributes 1/nne
        if(ic3(j,ie)/=0) then
          stmp=stmp+0.5_rkind !each side contributes 1/2
        else
          stmp=stmp+1.0_rkind !bndry sides contribute 1
        endif
      enddo
      etmp=etmp*real(nlev(ie),rkind)
      ptmp=ptmp*real(nlev(ie),rkind)
      stmp=stmp*real(nlev(ie),rkind)
      vwgt(ncon*(ie-1)+1)=nint(etmp)   !weight 1: 3D element ctr
      vwgt(ncon*(ie-1)+2)=nint(ptmp)   !weight 2: 3D vertical edge; less weight for nodes with more neighbors (so more of those nodes will be put into a sub-domain to prevent edge cuts)
      vwgt(ncon*(ie-1)+3)=nint(stmp)   !weight 3: 3D vertical face (put more internal sides in one sub-domain)
      vwgt(ncon*(ie-1)+4)=1            !weight 4: 2D element
    endif !(wgtflag==2.or.wgtflag==3)
  enddo !ie=1,ne

  ! Build edge list for dual graph and assign edge weights 
  adjncy=0
  xadj(1)=1
  !In the non-aug. domain?
  do ie=1,ne
    nedge=0
    do j=1,i34(ie) !node
      ip=elnode(j,ie)
      do k=1,nne(ip)
        je=indel(k,ip)
        if(je/=ie) then
          found=.false.
          do l=xadj(ie),xadj(ie)+nedge-1
            if(adjncy(l)==je) then
              !side sharing
              !contribute i34-2 nodes and i34-1 sides to comm
              !Reduce weight adjwgt()?
              if(wgtflag==1.or.wgtflag==3) adjwgt(l)=nlev(je)*(2*i34(je)-3)
              found=.true.
              exit
            endif
          enddo !l
          if(.not.found) then
            adjncy(xadj(ie)+nedge)=je
            !node sharing
            !contribute i34-1 nodes and i34 sides to comm
            if(wgtflag==1.or.wgtflag==3) adjwgt(xadj(ie)+nedge)=nlev(je)*(2*i34(je)-1)
            nedge=nedge+1
          endif
        endif !je/=ie
      enddo !k=1,nne(ip)
    enddo !j=1,
    xadj(ie+1)=xadj(ie)+nedge
  enddo !ie=1,ne

  !Output
  open(10,file='graphinfo',status='replace')
  !Edge cuts in METIS count only once for 2 of the same edge btw 2 elements
  !Vertex/edge counts start from 1
  write(10,*)ne,ntedge/2,'011',ncon
  do ie=1,ne
    !METIS 5.1 manual is wrong: no vertex size if '011'
!    write(10,'(i2,4(1x,i5),100000(1x,i11,1x,i5))')1,vwgt(ncon*(ie-1)+1:ncon*(ie-1)+4),(adjncy(j),adjwgt(j),j=xadj(ie),xadj(ie+1)-1)
    write(10,'(4(1x,i5),100000(1x,i11,1x,i5))')vwgt(ncon*(ie-1)+1:ncon*(ie-1)+4),(adjncy(j),adjwgt(j),j=xadj(ie),xadj(ie+1)-1)
  enddo !ie=1,ne
  close(10)

  ! ParMeTiS vertex distribution array (starts at 1)
!  allocate(vtxdist(nproc+1),stat=stat)
!  vtxdist(nproc+1)=ne+1
!
!  ! Fortran-style numbering that starts from 1
!  numflag = 1
!
!  ! Number of dimensions of the space in which the graph is embedded
!  ndims = 2
!
!  ! Vertex weight fraction
!  allocate(tpwgts(ncon*nproc),stat=stat)
!  tpwgts=1.0/real(nproc)
!
!  ! Imbalance tolerance (1: perfect balance; nproc: perfect imbalance)
!  allocate(ubvec(ncon),stat=stat)
!  ubvec=1.01

  ! MeTiS options
!  options(0)=1   ! 0: default options; 1: user options
!  options(1)=0   ! Level of information returned: see defs.h in ParMETIS-Lib dir
!  options(2)=15  ! Random number seed

!  call METIS_SetDefaultOptions(options)
!  options[METIS_OPTION_NUMBERING]=1 !FORTRAN style
!  options[METIS_OPTION_SEED]=15 !Random number seed

  ! Partition array returned from ParMeTiS
!  allocate(part(ne),stat=stat)

!int METIS PartGraphKway(idx t *nvtxs, idx t *ncon, idx t *xadj, idx t *adjncy,
!  idx t *vwgt, idx t *vsize, idx t *adjwgt, idx t *nparts, real t *tpwgts,
!  real t ubvec, idx t *options, idx t *objval, idx t *part)
!  call PartGraphKway(ne,ncon,xadj,adjncy,vwgt,NULL,adjwgt,nproc,tpwgts,ubvec,options,objval,part)

!  part=part-1 !partition numbering starts from 0
!  do i=1,ne
!    write(90,*)i,part(i)
!  enddo !i

  !Dealloc

  stop
  end
