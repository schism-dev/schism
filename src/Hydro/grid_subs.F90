!   Copyright 2014 College of William and Mary
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.


!===============================================================================
!===============================================================================
! SCHISM GRID SUBROUTINES
!
! subroutine aquire_vgrid
! subroutine partition_hgrid
! subroutine aquire_hgrid
! function signa
! subroutine dump_hgrid
!
!===============================================================================
!===============================================================================

subroutine aquire_vgrid
!-------------------------------------------------------------------------------
! Aquire vertical grid data from vgrid.in
!-------------------------------------------------------------------------------
  use schism_glbl
  use schism_msgp
  implicit none
  integer :: i,j,k,l,jki,stat,kin,m
  real(rkind) :: buf1(100),hmod2,zz

  !ivcor: types of vertical coord.; surface must all be nvrt (for sflux routines)
  !ivcor=2 !SZ coordinates

!  if(lm2d) then
!    !2D
!    ivcor=2; nvrt=2; kz=1; nsig=2; h_s=1.e6; h_c=h_s !5
!    theta_b=0; theta_f=1.e-4; s_con1=sinh(theta_f)
!    allocate(ztot(nvrt),sigma(nvrt),cs(nvrt),dcs(nvrt),stat=stat) !for outputs only
!    if(stat/=0) call parallel_abort('VGRID: ztot allocation failure')
!    ztot(kz)=-h_s
!    sigma(1)=-1 !bottom
!    sigma(nsig)=0 !surface
!  else !3D
  open(19,file=in_dir(1:len_in_dir)//'vgrid.in',status='old',iostat=stat)
  if(stat/=0) call parallel_abort('AQUIRE_VGIRD: open(19) failure')
  read(19,*)ivcor
  if(ivcor==2) then !SZ coordinates
    read(19,*) nvrt,kz,h_s !kz>=1
    if(nvrt<2) call parallel_abort('nvrt<2')
    if(kz<1) then !.or.kz>nvrt-2) then
      write(errmsg,*)'Wrong kz:',kz
      call parallel_abort(errmsg)
    endif
    if(h_s<6) then
      write(errmsg,*)'h_s needs to be larger:',h_s
      call parallel_abort(errmsg)
    endif

    ! Allocate vertical layers arrays
    allocate(ztot(nvrt),sigma(nvrt),cs(nvrt),dcs(nvrt),stat=stat)
    if(stat/=0) call parallel_abort('AQUIRE_VGIRD: ztot allocation failure')

    ! # of z-levels excluding "bottom" at h_s
    read(19,*) !for adding comment "Z levels"
    do k=1,kz-1
      read(19,*)j,ztot(k)
      if(ztot(k)>=-h_s) then
        write(errmsg,*)'Illegal Z level:',k
        call parallel_abort(errmsg)
      endif
      if(k>1) then; if(ztot(k)<=ztot(k-1)) then
        write(errmsg,*)'z-level inverted:',k
        call parallel_abort(errmsg)
      endif; endif
    enddo !k
    read(19,*) !level kz       
    ! In case kz=1, there is only 1 ztot(1)=-h_s
    ztot(kz)=-h_s

    nsig=nvrt-kz+1 !# of S levels (including "bottom" & f.s.)
    read(19,*) !for adding comment "S levels"
    read(19,*)h_c,theta_b,theta_f
    if(h_c<5._rkind.or.h_c>=h_s) then !large h_c to avoid 2nd type abnormaty
      write(errmsg,*)'h_c needs to be larger:',h_c
      call parallel_abort(errmsg)
    endif
    if(theta_b<0._rkind.or.theta_b>1._rkind) then
      write(errmsg,*)'Wrong theta_b:',theta_b
      call parallel_abort(errmsg)
    endif
    if(theta_f<=0._rkind) then
      write(errmsg,*)'Wrong theta_f:',theta_f
      call parallel_abort(errmsg)
    endif
    !Pre-compute constants
    s_con1=sinh(theta_f)

    sigma(1)=-1._rkind !bottom
    sigma(nsig)=0._rkind !surface
    read(19,*) !level kz
    do k=kz+1,nvrt-1
      kin=k-kz+1
      read(19,*) j,sigma(kin)
      if(sigma(kin)<=sigma(kin-1).or.sigma(kin)>=0._rkind) then
        write(errmsg,*)'Check sigma levels at:',k,sigma(kin),sigma(kin-1)
        call parallel_abort(errmsg)
      endif
    enddo !k
    read(19,*) !level nvrt
    close(19)
  else if(ivcor==1) then !localized sigma
    read(19,*)nvrt !needs hgrid to read the rest
    close(19)
    allocate(ztot(nvrt),sigma(nvrt),stat=stat)
    if(stat/=0) call parallel_abort('AQUIRE_VGIRD: ztot allocation failure (2)')
    !for output only - remove later
    ztot=0._rkind; sigma=0._rkind 
    kz=1; h_s=0._rkind; h_c=0._rkind; theta_b=0._rkind; theta_f=0._rkind
  else
    call parallel_abort('VGRID: Unknown ivcor')
  endif !ivcor
!  endif !lm2d

  if(ivcor==2) then
    ! Compute C(s) and C'(s)
    do k=1,nsig
      cs(k)=(1._rkind-theta_b)*sinh(theta_f*sigma(k))/sinh(theta_f)+ &
       &theta_b*(tanh(theta_f*(sigma(k)+0.5_rkind))-tanh(theta_f*0.5_rkind))/2._rkind/tanh(theta_f*0.5_rkind)
      dcs(k)=(1._rkind-theta_b)*theta_f*cosh(theta_f*sigma(k))/sinh(theta_f)+ &
       &theta_b*theta_f/2._rkind/tanh(theta_f*0.5_rkind)/cosh(theta_f*(sigma(k)+0.5_rkind))**2._rkind
    enddo !k
  endif !ivcor==2

  ! Output some sample z-coordinates
!  if(myrank==0) then
!    open(10,file='sample_z.out',status='replace')
!    write(10,*)'Sample z coordinates'
!    buf1(1)=h_s; buf1(2)=h_c; buf1(2:11)=(/(10*(i-1),i=2,11)/); buf1(12:28)=(/(200+50*(i-12), i=12,28)/)
!    write(10,*)'h_c= ',h_c,', h_s=',h_s
!    do i=1,28
!      write(10,*)'Depth= ',buf1(i)
!      do k=kz,nvrt
!        kin=k-kz+1
!        hmod2=min(buf1(i),h_s)
!        if(hmod2<=h_c) then
!          zz=sigma(kin)*hmod2
!        else
!          zz=h_c*sigma(kin)+(hmod2-h_c)*cs(kin)
!        endif
!        write(10,*)k,zz
!      enddo !k
!    enddo !i
!    close(10)
!  endif

end subroutine aquire_vgrid

!===============================================================================
!===============================================================================

subroutine partition_hgrid
!-------------------------------------------------------------------------------
! Compute a load-balanced partition of the horizontal grid using ParMeTiS
! graph partitioning routine applied to dual-graph constructed from elements.
! > vertical grid must be aquired prior to calling
!-------------------------------------------------------------------------------
!#ifdef USE_MPIMODULE
!  use mpi
!#endif
  use schism_glbl
  use schism_msgp
  implicit none
!#ifndef USE_MPIMODULE
  include 'mpif.h'
!#endif

  ! ParMeTiS
  integer :: wgtflag,numflag,ncon,ndims,options(0:2)
  integer :: mxnedge,nedge,ntedge,edgecut
  integer,allocatable :: vtxdist(:),xadj(:),adjncy(:),part(:)
  real(4),allocatable :: xyz(:),xproj(:),yproj(:)
  integer,allocatable :: nlev(:),vwgt(:),adjwgt(:)
  real(4),allocatable :: tpwgts(:),ubvec(:)

  ! Miscellaneous
  integer :: i,j,k,l,ip,ie,je,iegb,jegb,stat,kbetmp
  real(rkind) :: xtmp,ytmp,dtmp,tmp,etmp,ptmp,stmp
  integer,allocatable :: eprocv(:),neproc(:),neprocsum(:)
  logical :: found
!-------------------------------------------------------------------------------

  ! Open global grid file and read global grid size
  open(14,file=in_dir(1:len_in_dir)//'hgrid.gr3',status='old')
  read(14,*); read(14,*) ne_global,np_global
  close(14)

  ! Allocate global element to resident processor vector (in global module)
  if(allocated(iegrpv)) deallocate(iegrpv)
  allocate(iegrpv(ne_global),stat=stat)
  if(stat/=0) call parallel_abort('partition: iegrpv allocation failure')

  ! Handle single processor case
  if(nproc==1) then
    ne=ne_global
    iegrpv=0
    return
  endif

  ! Setup initial partition
  ! Equal # of elements in each processor (except for the last one).
  allocate(neproc(0:nproc-1),stat=stat)
  if(stat/=0) call parallel_abort('partition: neproc allocation failure')
  allocate(neprocsum(0:nproc-1),stat=stat)
  if(stat/=0) call parallel_abort('partition: neprocsum allocation failure')
  ne=ne_global/nproc
  neprocsum(0)=0
  do l=0,nproc-2
    neproc(l)=ne !# of resident elements
    neprocsum(l+1)=neprocsum(l)+ne
  enddo
  neproc(nproc-1)=ne_global-neprocsum(nproc-1)
  do i=1,ne_global
    iegrpv(i)=min((i-1)/ne,nproc-1)
  enddo
  ne=neproc(myrank)

  ! Aquire horizontal grid based on initial partition
  call aquire_hgrid(.false.)

  !The following needs info 
  !from aquire_hgrid: i34, ielg; elnode; ic3; nne; indel; ne; nea; npa; xnd, ynd, znd, dp
  !from aquire_vgrid: nvrt; ztot; kz; h_s.
  !from schism_init: ics

  !Do map projection for lat/lon
  if(allocated(xproj)) deallocate(xproj)
  if(allocated(yproj)) deallocate(yproj)
  allocate(xproj(npa),yproj(npa),stat=stat)
  if(stat/=0) call parallel_abort('partition: xproj allocation failure')
  if(ics==1) then
    xproj=xnd; yproj=ynd
  else !lat/lon
    do i=1,npa
      !Stereographic projection for partition only
      !Center is at south pole and plane is at north pole; must remove south pole
      if(rearth_pole+znd(i)<=0._rkind) then
        write(errmsg,*)'PARTITION: remove south pole:',i,rearth_pole+znd(i)
        call parallel_abort(errmsg)
      endif
      xproj(i)=2._rkind*rearth_eq*xnd(i)/(rearth_pole+znd(i))
      yproj(i)=2._rkind*rearth_eq*ynd(i)/(rearth_pole+znd(i))
    enddo !i=1,npa
  endif !ics

  !Read in kbp from vgrid.in if ivcor=1 for partioning
  !This will be temporary as it'll be recomputed once the final partioning is
  !done
  if(ivcor==1) then
    allocate(kbp(npa))
    open(19,file=in_dir(1:len_in_dir)//'vgrid.in',status='old')
    read(19,*); read(19,*)nvrt
    do i=1,np_global
      read(19,*)j,kbetmp
      if(ipgl(i)%rank==myrank) kbp(ipgl(i)%id)=kbetmp
    enddo !i
    close(19)
  endif !ivcor==1

  ! Count number of edges in dual graph
  allocate(adjncy(1000),stat=stat) !for single element
  if(stat/=0) call parallel_abort('partition: adjncy(1000) allocation failure')
  ntedge=0 !total # of edges in the grid
  mxnedge=0 !max. # of local edges
  do ie=1,ne
    iegb=ielg(ie)
    nedge=0 !# of local edges
    adjncy=0 !list of global element indices
    do j=1,i34(ie)
      ip=elnode(j,ie)
      do k=1,nne(ip)
        jegb=ielg(indel(k,ip)) !indel(k,ip) inside aug. domain
        if(jegb/=iegb) then
          found=.false.
          do l=1,nedge
            if(adjncy(l)==jegb) then
              found=.true.
              exit
            endif
          enddo
          !if(l==nedge+1) then
          if(.not.found) then
            nedge=nedge+1
            if(nedge>1000) call parallel_abort('PARTITION: bound (1)')
            adjncy(nedge)=jegb
          endif
        endif
      enddo !k
    enddo !j
    ntedge=ntedge+nedge
    mxnedge=max(mxnedge,nedge)
  enddo !ie
  deallocate(adjncy)

  ! Use vertex and/or edge weights
  wgtflag = 3   ! 0: none; 1: edges; 2: vertices; 3: vertices & edges

  ! Number of weights associated with each vertex, used to optimize the
  ! partition
  ncon = 4

  ! Allocate storage for dual graph
  if(allocated(xadj)) deallocate(xadj); allocate(xadj(nea+1),stat=stat);
  if(stat/=0) call parallel_abort('partition: xadj allocation failure')
  if(allocated(adjncy)) deallocate(adjncy); allocate(adjncy(ntedge),stat=stat)
  if(stat/=0) call parallel_abort('partition: adjncy allocation failure')
  if(allocated(xyz)) deallocate(xyz); allocate(xyz(2*nea),stat=stat)
  if(stat/=0) call parallel_abort('partition: xyz allocation failure')
  if(allocated(nlev)) deallocate(nlev); allocate(nlev(nea),stat=stat) !estimate of number of active levels
  if(stat/=0) call parallel_abort('partition: nlev allocation failure')
  if(wgtflag==2.or.wgtflag==3) then
    if(allocated(vwgt)) deallocate(vwgt); allocate(vwgt(nea*ncon),stat=stat)
    if(stat/=0) call parallel_abort('partition: vwgt allocation failure')
  endif
  if(wgtflag==1.or.wgtflag==3) then
    if(allocated(adjwgt)) deallocate(adjwgt); allocate(adjwgt(ntedge),stat=stat)
    if(stat/=0) call parallel_abort('partition: adjwgt allocation failure')
  endif

  ! Assign vertex coordinates & weights
  ! Weights based on estimated number of active levels
  do ie=1,nea
    xtmp=0._rkind !xctr
    ytmp=0._rkind
    dtmp=real(5.d10,rkind) !min. depth
    do j=1,i34(ie)
      ip=elnode(j,ie)
      xtmp=xtmp+real(xproj(ip),rkind)/real(i34(ie),rkind)
      ytmp=ytmp+real(yproj(ip),rkind)/real(i34(ie),rkind)
      if(dp(ip)<dtmp) dtmp=dp(ip)
    enddo !j
    xyz(2*(ie-1)+1)=xtmp
    xyz(2*(ie-1)+2)=ytmp
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
  enddo !ie=1,nea

  ! Build edge list for dual graph and assign edge weights
  adjncy=0
  xadj(1)=1
  !In the non-aug. domain?
  do ie=1,ne
    iegb=ielg(ie)
    nedge=0
    do j=1,i34(ie) !node
      ip=elnode(j,ie)
      do k=1,nne(ip)
        je=indel(k,ip)
        jegb=ielg(je)
        if(jegb/=iegb) then
          found=.false.
          do l=xadj(ie),xadj(ie)+nedge-1
            if(adjncy(l)==jegb) then
              !side sharing
              !contribute i34-2 nodes and i34-1 sides to comm
              !Reduce weight adjwgt()?
              if(wgtflag==1.or.wgtflag==3) adjwgt(l)=nlev(je)*(2*i34(je)-3)
              found=.true.
              exit
            endif
          enddo !l
          if(.not.found) then
            adjncy(xadj(ie)+nedge)=jegb
            !node sharing
            !contribute i34-1 nodes and i34 sides to comm
            if(wgtflag==1.or.wgtflag==3) adjwgt(xadj(ie)+nedge)=nlev(je)*(2*i34(je)-1)
            nedge=nedge+1
          endif
        endif !jegb/=iegb
      enddo !k=1,nne(ip)
    enddo !j=1,
    xadj(ie+1)=xadj(ie)+nedge
  enddo !ie=1,ne

  ! ParMeTiS vertex distribution array (starts at 1)
  allocate(vtxdist(nproc+1),stat=stat)
  if(stat/=0) call parallel_abort('partition: vtxdist allocation failure')
  call mpi_allgather(neprocsum(myrank)+1,1,itype,vtxdist,1,itype,comm,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('partition: mpi_allgather',ierr)
  vtxdist(nproc+1)=ne_global+1

  ! Fortran-style numbering that starts from 1
  numflag = 1

  ! Number of dimensions of the space in which the graph is embedded
  ndims = 2

  ! Vertex weight fraction
  allocate(tpwgts(ncon*nproc),stat=stat)
  if(stat/=0) call parallel_abort('partition: tpwgts allocation failure')
  tpwgts=1.0/real(nproc)

  ! Imbalance tolerance
  allocate(ubvec(ncon),stat=stat)
  if(stat/=0) call parallel_abort('partition: ubvec allocation failure')
  ubvec=1.01

  ! ParMeTiS options
!#ifdef DEBUG
!  options(0)=1   ! 0: default options; 1: user options
!  options(1)=15  ! Level of information returned: see defs.h in ParMETIS-Lib dir
!  options(2)=15  ! Random number seed
!  if(myrank==0) write(*,'(/a)') 'ParMETIS Partitioning:'
!#else
  options(0)=1   ! 0: default options; 1: user options
  options(1)=0   ! Level of information returned: see defs.h in ParMETIS-Lib dir
  options(2)=15  ! Random number seed
!#endif

  ! Partition array returned from ParMeTiS
  allocate(part(ne),stat=stat)
  if(stat/=0) call parallel_abort('partition: part allocation failure')

  ! Partition dual graph
!  My notes from manual:
!  p: # of processors;
!  n: total # of vertices (local) in graph sense;
!  m: total # of neighboring vertices ("edges"); double counted between neighboring vertice u and v.
!  ncon: # of weights for each vertex;
!  int(in) vtxdist(p+1): Processor j stores vertices vtxdist(j):vtxdist(j+1)-1
!  int (in) xadj(n+1), adjncy(m): locally, vertex j's neighboring vertices are adjncy(xadj(j):xadj(j+1)-1). adjncy points to global index;
!  int(in) vwgt(ncon*n), adjwgt(m): weights at vertices and "edges". Format of adjwgt follows adjncy;
!  int(in) wgtflag: 0: none (vwgt and adjwgt are NULL); 1: edges (vwgt is NULL); 2: vertices (adjwgt is NULL); 3: both vertices & edges;
!  int(in) numflag: 0: C-style numbering from 0; 1: FORTRAN style from 1;
!  int(in) ndims: 2 or 3 (D);
!  float(in) xyz(ndims*n): coordinate for vertex j is xyz(j*ndims:(j+1)*ndims-1);
!  int(in) nparts: # of desired sub-domains (usually nproc);
!  float(in) tpwgts(ncon*nparts): =1/nparts if sub-domains are to be of same size for each vertex weight; useful for heterogeneous cluster
!  float(in) ubvec(ncon): imbalance tolerance for each weight;
!  int(in) options: additonal parameters for the routine (see above);
!  int(out) edgecut: # of edges that are cut by the partitioning;
!  int(out) part(): array size = # of local vertices. It stores indices of local vertices.
 
  call ParMETIS_V3_PartGeomKway(vtxdist,xadj,adjncy,vwgt,adjwgt,wgtflag, &
                            &numflag,ndims,xyz,ncon,nproc,tpwgts,ubvec,options, &
                            &edgecut,part,comm)

  ! Change partition numbering to start from 0
  part=part-1

  ! Construct global element-to-resident-processor assignment vector
  call mpi_allgatherv(part,ne,itype,iegrpv,neproc,neprocsum,itype,comm,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('partition: mpi_allgatherv',ierr)

  ! Deallocate neproc arrays
  deallocate(neproc,neprocsum)

  ! Deallocate ParMeTiS arrays
  deallocate(vtxdist,xadj,adjncy,part)
  deallocate(xyz,nlev,tpwgts,ubvec,xproj,yproj)
  if(wgtflag==2.or.wgtflag==3) deallocate(vwgt)
  if(wgtflag==1.or.wgtflag==3) deallocate(adjwgt)
  if(ivcor==1) deallocate(kbp)

end subroutine partition_hgrid

!===============================================================================
!===============================================================================

subroutine aquire_hgrid(full_aquire)
!-------------------------------------------------------------------------------
! Based on a partition specified by iegrpv(i) aquire horizontal grid for
! augmented subdomain (resident + ghost elements/node/sides) from hgrid.gr3
! > construct global-to-local and local-to-global element, node and side index tables
! > construct geometry for augmented subdomain
! > aquire boundary segments and construct boundary segment mappings
!-------------------------------------------------------------------------------
!#ifdef USE_MPIMODULE
!  use mpi
!#endif
  use schism_glbl
  use schism_msgp
  use schism_io
  implicit none
!#ifndef USE_MPIMODULE
  include 'mpif.h'
!#endif
  logical,intent(in) :: full_aquire  ! Aquire and construct all tables
  integer :: i,j,k,l,ii,jj,irank,ip,jp,ie,je,ic,iegb,jegb,ngb1,ngb2,ipgb,isgb,icount,new,id,isd
  integer,allocatable :: ibuf(:),isbuf(:),irbuf(:)
  real(rkind),allocatable :: dbuf1(:),dbuf2(:)
  type(llist_type),pointer :: llp,node,nodep,side,sidep
  logical :: found,local1,local2,found1,found2
  real(rkind),parameter :: deg2rad=pi/180._rkind
  integer :: n1,n2,n3,n4,jsj,nt,nn,nscnt,nsgcnt
  real(rkind) :: ar1,ar2,ar3,ar4,signa !slam0,sfea0,
  real(rkind) :: xtmp,ytmp,dptmp,thetan,egb1,egb2,egb
!  integer :: ipre 

  integer,allocatable :: neproc(:)   ! Number of resident elements on each processor
  integer,allocatable :: npproc(:)   ! Number of resident nodes on each processor
  integer,allocatable :: nsproc(:)   ! Number of resident sides on each processor

  integer,allocatable :: i34gb(:),nmgb(:,:)  ! Global element-node tables
  integer,allocatable :: ic3gb(:,:)          ! Global element-side-element table
  integer,allocatable :: jsgb(:,:)           ! Global element-side table
  integer,allocatable :: nnegb(:),inegb(:,:) ! Global node-element tables
  integer,allocatable :: isbnd_global(:) ! Node to open bndry segment flags (global)
!  integer,allocatable :: nnpgb(:),inpgb(:,:) ! Global node-node tables
!  integer,allocatable :: iselfgb(:,:)        ! Global node-element self-reference table

  integer :: npi                           ! Number of interface nodes
  integer,allocatable :: ieg(:)            ! List of ghost elements (global index)
  integer,allocatable :: ipg(:)            ! List of ghost nodes (global index)
  integer,allocatable :: isg(:)            ! List of ghost sides (global index)
  integer,allocatable :: ieg2(:)            ! List of 2-tier ghost elements (global index)
  integer,allocatable :: ipg2(:)            ! List of 2-tier ghost nodes (global index)
  integer,allocatable :: isg2(:)            ! List of 2-tier ghost sides (global index)

  integer :: mnopep     ! Max number of partitions of local open boundary segments

  integer :: stat

  integer :: intvalue
  real(rkind) :: realvalue
  character(len=2) :: stringvalue
!-------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! Initialize
  !-----------------------------------------------------------------------------

  !mntr=max(ntracers,2)
  !ntracers2=ntracers+2

  !-----------------------------------------------------------------------------
  ! Open global grid file
  !-----------------------------------------------------------------------------
  open(14,file=in_dir(1:len_in_dir)//'hgrid.gr3',status='old',iostat=stat)
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: open(14) failure')

  !-----------------------------------------------------------------------------
  ! Aquire and construct global data tables
  !-----------------------------------------------------------------------------

  ! Aquire global grid size
  read(14,*); read(14,*) ne_global,np_global

  ! Aquire global element-node tables from hgrid.gr3
  if(allocated(i34gb)) deallocate(i34gb); allocate(i34gb(ne_global),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: i34gb allocation failure')
  if(allocated(nmgb)) deallocate(nmgb); allocate(nmgb(4,ne_global),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: nmgb allocation failure')
  do i=1,np_global; read(14,*); enddo;
  do i=1,ne_global
    read(14,*) iegb,i34gb(iegb),(nmgb(k,iegb),k=1,i34gb(iegb))
    if(i34gb(iegb)/=3.and.i34gb(iegb)/=4) then
      write(errmsg,*) 'AQUIRE_HGRID: Unknown type of element',iegb,i34gb(iegb)
      call parallel_abort(errmsg)
    endif
    if(i34gb(iegb)==4) lhas_quad=.true.
  enddo
  if(myrank==0) write(16,*)'lhas_quad=',lhas_quad

  ! Count number of elements connected to each node (global)
  if(allocated(nnegb)) deallocate(nnegb); allocate(nnegb(np_global),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: nnegb allocation failure')
  nnegb=0
  do iegb=1,ne_global
    do k=1,i34gb(iegb)
      ipgb=nmgb(k,iegb)
      nnegb(ipgb)=nnegb(ipgb)+1
    enddo !k
  enddo !iegb

  ! Check hanging nodes
  if(myrank==0) then
    found=.false.
    do ipgb=1,np_global
      if(nnegb(ipgb)==0) then
        found=.true.
        write(11,*) 'Hanging node:',ipgb
      endif
    enddo
    if(found) call parallel_abort('AQUIRE_HGRID: check fort.11 for hanging nodes')
  endif

  ! Maximum number of elements connected to a node
  mnei=0
  do ipgb=1,np_global
    mnei=max(mnei,nnegb(ipgb))
  enddo

  ! Build global node-element and self-reference table table
  if(allocated(inegb)) deallocate(inegb); allocate(inegb(np_global,mnei),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: inegb allocation failure')
!  if(allocated(iselfgb)) deallocate(iselfgb); allocate(iselfgb(mnei,np_global),stat=stat);
!  if(stat/=0) call parallel_abort('AQUIRE_HGRID: iselfgb allocation failure')
  nnegb=0
  do iegb=1,ne_global
    do k=1,i34gb(iegb)
      ipgb=nmgb(k,iegb)
      nnegb(ipgb)=nnegb(ipgb)+1
      inegb(ipgb,nnegb(ipgb))=iegb
!      iselfgb(nnegb(ipgb),ipgb)=k
    enddo !k
  enddo !iegb

  !Compute mnei_p (max. # of nodes around a node)
  mnei_p=0
  do i=1,np_global
    icount=0
    do j=1,nnegb(i)
      iegb=inegb(i,j)
      icount=icount+i34gb(iegb)-2 !# of new surrounding nodes
    enddo !j
    mnei_p=max(mnei_p,icount+1) !account for bnd ball
  enddo !i

  if(myrank==0.and..not.full_aquire) write(16,*)'mnei, mnei_p = ',mnei,mnei_p

  ! Build global element-side-element table; this won't be affected by re-arrangement below
  if(allocated(ic3gb)) deallocate(ic3gb); allocate(ic3gb(4,ne_global),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: ic3gb allocation failure')
  do iegb=1,ne_global
    do k=1,i34gb(iegb)
      ic3gb(k,iegb)=0 !index for boundary sides
      ngb1=nmgb(nxq(1,k,i34gb(iegb)),iegb)
      ngb2=nmgb(nxq(2,k,i34gb(iegb)),iegb)
      do l=1,nnegb(ngb1)
        jegb=inegb(ngb1,l)
        if(jegb/=iegb.and.(nmgb(1,jegb)==ngb2.or.nmgb(2,jegb)==ngb2.or. &
     &nmgb(3,jegb)==ngb2.or.(i34gb(jegb)==4.and.nmgb(4,jegb)==ngb2))) ic3gb(k,iegb)=jegb
      enddo !l
      jegb=ic3gb(k,iegb)
      if(jegb/=0) then
        do l=1,i34gb(jegb)
          if(nmgb(nxq(1,l,i34gb(jegb)),jegb)==ngb1.and.nmgb(nxq(2,l,i34gb(jegb)),jegb)==ngb2) then
            write(errmsg,*) 'Elem ', iegb, ' and ', jegb, ' have opposite orientation'
            call parallel_abort(errmsg)
          endif
        end do  !l
      endif
    enddo !k
  enddo !iegb

  ! inegb to be re-arrange in counter-clockwise fashion after boundary info is read in
  ! Count global number of sides and build global element-side index table
  if(allocated(jsgb)) deallocate(jsgb); allocate(jsgb(4,ne_global),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: jsgb allocation failure')
  ns_global=0
  do iegb=1,ne_global
    do j=1,i34gb(iegb) !visit each side associated with element iegb
      if(ic3gb(j,iegb)==0.or.iegb<ic3gb(j,iegb)) then !new global side
        ns_global=ns_global+1
        jsgb(j,iegb)=ns_global
        if(ic3gb(j,iegb)/=0) then !old internal side
          jegb=ic3gb(j,iegb)
          l=0
          do k=1,i34gb(jegb)
            if(ic3gb(k,jegb)==iegb) then
              l=k
              exit
            endif
          enddo !k
          if(l==0) then
            write(errmsg,'(a,10i6)') 'AQUIRE_HGRID: Wrong ball info',iegb,j,ns_global
            call parallel_abort(errmsg)
          endif
          jsgb(l,jegb)=ns_global
        endif !ic3gb(j,iegb)/=0
      endif !ic3gb(j,iegb)==0.or.iegb<ic3gb(j,iegb)
    enddo !j
  enddo !iegb
  if(ns_global.lt.ne_global.or.ns_global.lt.np_global) then
    write(errmsg,*)'AQUIRE_HGRID: weird grid with ns_global < ne_global or ns_global < np_global', &
    &np_global,ne_global,ns_global
    call parallel_abort(errmsg)
  endif
  if(full_aquire) then
    if(myrank==0) then
      write(16,'(/a,4i10)')'Global Grid Size (ne,np,ns,nvrt): ',ne_global,np_global,ns_global,nvrt
    endif
    call parallel_barrier
  endif

  !-----------------------------------------------------------------------------
  ! Aquire global open boundary segments from hgrid.gr3
  !-----------------------------------------------------------------------------

  ! Allocate and assign global node-to-open-boundary-segment flags
  !  isbnd_global = 0 (internal); >0 (open bnd segment #); -1 (land bnd)
  if(allocated(isbnd_global)) deallocate(isbnd_global);
  allocate(isbnd_global(np_global),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: isbnd_global allocation failure')
  isbnd_global=0;

  ! Global number of open boundary segments and nodes
  rewind(14); read(14,*); read(14,*);
  do i=1,np_global; read(14,*); enddo;
  do i=1,ne_global; read(14,*); enddo;
  read(14,*) nope_global
  read(14,*) neta_global

  ! Scan segments to count number of open boundary segments and nodes
  mnond_global=0 !global max number of nodes per segment
  nt=0    !global total node count
  do k=1,nope_global
    read(14,*) nn
    mnond_global=max(mnond_global,nn);
    nt=nt+nn
    do i=1,nn; read(14,*); enddo;
  enddo !k
  if(neta_global/=nt) then
    write(errmsg,*) 'neta_global /= total # of open bnd nodes',neta_global,nt
    call parallel_abort(errmsg)
  endif

  ! Allocate arrays for global open boundary segments
  if(allocated(nond_global)) deallocate(nond_global);
  allocate(nond_global(nope_global),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: nond_global allocation failure')
  if(allocated(iond_global)) deallocate(iond_global);
  allocate(iond_global(nope_global,mnond_global),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: iond_global allocation failure')

  ! Aquire global open boundary segments and nodes
  rewind(14); read(14,*); read(14,*);
  do i=1,np_global; read(14,*); enddo;
  do i=1,ne_global; read(14,*); enddo;
  read(14,*); read(14,*);
  nond_global=0; iond_global=0;
  do k=1,nope_global
    read(14,*) nn
    do i=1,nn
      read(14,*) ipgb
      nond_global(k)=nond_global(k)+1
      iond_global(k,nond_global(k))=ipgb
      isbnd_global(ipgb)=k
    enddo !i
    if(iond_global(k,1)==iond_global(k,nond_global(k))) then
      write(errmsg,*) 'Looped open bnd:',k
      call parallel_abort(errmsg)
    endif
  enddo !k

  !-----------------------------------------------------------------------------
  ! Aquire global land boundary segments from hgrid.gr3
  !-----------------------------------------------------------------------------

  ! Global total number of land boundary segments and nodes
  rewind(14); read(14,*); read(14,*);
  do i=1,np_global; read(14,*); enddo;
  do i=1,ne_global; read(14,*); enddo;
  read(14,*); read(14,*);
  do k=1,nope_global; read(14,*) nn; do i=1,nn; read(14,*); enddo; enddo;
  read(14,*) nland_global
  read(14,*) nvel_global

  ! Scan segments to count number of land boundary segments and nodes
  mnlnd_global=0 !global max number of nodes per segment
  nt=0    !global total node count
  do k=1,nland_global
    read(14,*) nn
    mnlnd_global=max(mnlnd_global,nn)
    nt=nt+nn
    do i=1,nn; read(14,*); enddo;
  enddo !k
  if(nvel_global/=nt) then
    write(errmsg,*) 'AQUIRE_HGRID: nvel_global /= total # of land bnd nodes', &
                    &nvel_global,nt
    call parallel_abort(errmsg)
  endif

  ! Allocate arrays for global land boundary segments
  if(allocated(nlnd_global)) deallocate(nlnd_global);
  allocate(nlnd_global(nland_global),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: nlnd_global allocation failure')
  if(allocated(ilnd_global)) deallocate(ilnd_global);
  allocate(ilnd_global(nland_global,mnlnd_global),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: ilnd_global allocation failure')

  ! Aquire global land boundary segments and nodes
  rewind(14); read(14,*); read(14,*);
  do i=1,np_global; read(14,*); enddo;
  do i=1,ne_global; read(14,*); enddo;
  read(14,*); read(14,*);
  do k=1,nope_global; read(14,*) nn; do i=1,nn; read(14,*); enddo; enddo;
  read(14,*); read(14,*);
  nlnd_global=0; ilnd_global=0;
  do k=1,nland_global
    read(14,*) nn
    do i=1,nn
      read(14,*) ipgb
      nlnd_global(k)=nlnd_global(k)+1
      ilnd_global(k,nlnd_global(k))=ipgb
      if(isbnd_global(ipgb)==0) isbnd_global(ipgb)=-1 !overlap of open bnd
    enddo !i
  enddo !k

! Re-arrange in counter-clockwise fashion
!  if(allocated(nnpgb)) deallocate(nnpgb);
!  allocate(nnpgb(np_global),stat=stat);
!  if(stat/=0) call parallel_abort('AQUIRE_HGRID: nnpgb allocation failure')
!  if(allocated(inpgb)) deallocate(inpgb);
!  allocate(inpgb(np_global,mnei+1),stat=stat);
!  if(stat/=0) call parallel_abort('AQUIRE_HGRID: inpgb allocation failure')

  do i=1,np_global
    if(isbnd_global(i)/=0) then !bnd ball
!     Look for starting bnd element
      icount=0
      do j=1,nnegb(i)
        ie=inegb(i,j)
!        ii=iselfgb(j,i)
        ii=0 !local index
        do l=1,i34gb(ie)
          if(nmgb(l,ie)==i) then
            ii=l; exit
          endif
        enddo !l
        if(ii==0) call parallel_abort('AQUIRE_HGRID: bomb (1)')

        if(ic3gb(nxq(i34gb(ie)-1,ii,i34gb(ie)),ie)==0) then
          icount=icount+1
          inegb(i,1)=ie
        endif
      enddo !j=1,nnegb(i)
      if(icount/=1) then
        write(errmsg,*)'Illegal bnd node',i,isbnd_global(i),icount
        call parallel_abort(errmsg)
      endif
    endif !bnd ball

!   For internal balls, starting elem. is not altered
!   Sequential search for the rest of elements
!    nnpgb(i)=2
!    inpgb(i,1)=nmgb(nx(iselfgb(1,i),1),inegb(i,1))
!    inpgb(i,2)=nmgb(nx(iselfgb(1,i),2),inegb(i,1))
    do j=2,nnegb(i)
      ie=inegb(i,j-1)
      ii=0 !local index
      do l=1,i34gb(ie)
        if(nmgb(l,ie)==i) then
          ii=l; exit
        endif
      enddo !l
      if(ii==0) call parallel_abort('AQUIRE_HGRID: bomb (2)')

      new=ic3gb(nxq(i34gb(ie)-2,ii,i34gb(ie)),ie)
      if(new==0) then
        write(errmsg,*)'Incomplete ball',i
        call parallel_abort(errmsg)
      endif
      inegb(i,j)=new
!      ii=0
!      do l=1,3
!        if(nmgb(l,new)==i) ii=l
!      enddo !l
!      if(ii==0) then
!        write(errmsg,*)'Failed to find local index:',i,new
!        call parallel_abort(errmsg)
!      endif
!      iselfgb(j,i)=ii
!
!      if(isbnd_global(i)==0.and.j==nnegb(i)) then !complete internal ball
!!	Check completeness
!        if(nmgb(nx(ii,2),new)/=inpgb(i,1)) then
!          write(errmsg,*)'Broken ball:',i
!          call parallel_abort(errmsg)
!        endif
!      else !one more node
!        nnpgb(i)=nnpgb(i)+1
!        if(nnpgb(i)>mnei+1) then
!          write(errmsg,*)'Too many neighbor nodes:',i,mnei
!          call parallel_abort(errmsg)
!        endif
!        inpgb(i,nnpgb(i))=nmgb(nx(ii,2),new)
!      endif
    enddo !j=2,nnegb(i)
  enddo !i=1,np_global
  !-----------------------------------------------------------------------------
  ! Done with global grid -- close grid file
  !-----------------------------------------------------------------------------
  close(14)

  !-----------------------------------------------------------------------------
  ! Count number of resident elements, nodes and sides for each processor.
  ! Build global-to-local index tables for resident element, nodes and sides.
  !-----------------------------------------------------------------------------

  ! Allocate and initialize resident elements/nodes/sides count arrays
  allocate(neproc(0:nproc-1),stat=stat)
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: neproc allocation failure')
  neproc=0
  allocate(npproc(0:nproc-1),stat=stat)
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: npproc allocation failure')
  npproc=0

  ! Allocate global-to-local element index table
  if(associated(iegl)) call release_gl(ne_global,iegl)
  allocate(iegl(ne_global),stat=stat)
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: iegl allocation failure')

  ! Allocate global-to-local node index table
  if(associated(ipgl)) call release_gl(np_global,ipgl)
  allocate(ipgl(np_global),stat=stat)
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: ipgl allocation failure')

  allocate(nsproc(0:nproc-1),stat=stat)
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: nsproc allocation failure')
  nsproc=0
  ! Allocate global-to-local side index table
  if(associated(isgl)) call release_gl(ns_global,isgl)
  allocate(isgl(ns_global),stat=stat)
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: isgl allocation failure')

  ! Allocate 2-tier arrays for later
  if(allocated(iegl2)) deallocate(iegl2); allocate(iegl2(2,ne_global),stat=stat)
  if(allocated(ipgl2)) deallocate(ipgl2); allocate(ipgl2(2,np_global),stat=stat)
  if(allocated(isgl2)) deallocate(isgl2); allocate(isgl2(2,ns_global),stat=stat)
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: iegl2 allocation failure')

  ! Build global-to-local element, node and side index tables
  if(nproc>1) then
    do iegb=1,ne_global
      irank=iegrpv(iegb)
      neproc(irank)=neproc(irank)+1
      ie=neproc(irank)
      iegl(iegb)%rank=irank
      iegl(iegb)%id=ie
      do k=1,i34gb(iegb) !visit each node/side attached to element iegb
        !Node table
        ipgb=nmgb(k,iegb) !ipgb is global node index
        if(ipgl(ipgb)%id==0) then !node not yet visited -- add to list
          npproc(irank)=npproc(irank)+1
          ipgl(ipgb)%rank=irank
          ipgl(ipgb)%id=npproc(irank) !local node index
        else !node visited at least once
          found=.false.
          node=>ipgl(ipgb)
          loopn: do !check if in linked-list
            if(node%rank==irank) then !already resident in processor irank (which may not be current processor myrank!)
              found=.true.
              exit loopn
            endif
            nodep=>node !save last entry in the list
            node=>node%next !search next entry in the list
            if(.not.associated(node)) exit loopn
          enddo loopn
          if(.not.found) then !add interface node to the end of linked-list
! nodep%next=>null()
            npproc(irank)=npproc(irank)+1
            allocate(node,stat=stat)
            if(stat/=0) call parallel_abort('AQUIRE_HGRID: node allocation failure')
            node=llist_type(irank,npproc(irank),null())
            nodep%next=>node
          endif !.not.found
        endif !ipgl(ipgb)%id==0

        !Side table
        isgb=jsgb(k,iegb) !isgb is global side index
        if(isgl(isgb)%id==0) then !side not yet visited -- add to list
          nsproc(irank)=nsproc(irank)+1
          isgl(isgb)%rank=irank
          isgl(isgb)%id=nsproc(irank) !local side index
        else !side visited at least once
          found=.false.
          side=>isgl(isgb)
          loops: do !check if in linked-list
            if(side%rank==irank) then
              found=.true.
              exit loops
            endif
            sidep=>side
            side=>side%next
            if(.not.associated(side)) exit loops
          enddo loops
          if(.not.found) then !add interface side to linked-list
            nsproc(irank)=nsproc(irank)+1
            allocate(side,stat=stat)
            if(stat/=0) call parallel_abort('AQUIRE_HGRID: side allocation failure')
            side=llist_type(irank,nsproc(irank),null())
            sidep%next=>side
          endif !.not.found
        endif !isgl(isgb)%id==0
      enddo !k=1,3
    enddo !iegb
  else !nproc==1
    do iegb=1,ne_global
      neproc(0)=neproc(0)+1
      iegl(iegb)%rank=0
      iegl(iegb)%id=iegb
    enddo !iegb
    do ipgb=1,np_global
      npproc(0)=npproc(0)+1
      ipgl(ipgb)%rank=0
      ipgl(ipgb)%id=ipgb
    enddo !iegb

    do isgb=1,ns_global
      nsproc(0)=nsproc(0)+1
      isgl(isgb)%rank=0
      isgl(isgb)%id=isgb
    enddo !iegb
  endif !nproc

  ! Place resident node/sides at front of node/side link-lists in ipgl/isgl
  if(nproc>1) then
    call swap_llrank(np_global,ipgl)
    call swap_llrank(ns_global,isgl)
  endif

  ! Sort rest of each node/side link-list in ascending order according to rank
  ! Hereafter, ipgl%rank=myrank (if resident in myrank), ipgl%id=local index in myrank, and ipgl%next%next%next... is the list
  if(nproc>1) then
    call sort_llrank(np_global,ipgl)
    call sort_llrank(ns_global,isgl)
  endif

  ! Number of resident elements, nodes and sides
  ne=neproc(myrank)
  np=npproc(myrank)
  ns=nsproc(myrank)

  ! Deallocated resident count arrays
  if(allocated(neproc)) deallocate(neproc)
  if(allocated(npproc)) deallocate(npproc)
  if(allocated(nsproc)) deallocate(nsproc)

  !-----------------------------------------------------------------------------
  ! Build lists of ghost elements, nodes and sides
  !-----------------------------------------------------------------------------

  ! Only build ghost tables if nproc>1
  if(nproc>1) then
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Count number of interface nodes (for allocation)
    npi=0
    do ipgb=1,np_global
      if(ipgl(ipgb)%rank==myrank.and.associated(ipgl(ipgb)%next)) npi=npi+1
    enddo !ipgb

    ! Build sorted list of ghost elements
    if(allocated(ieg)) deallocate(ieg); allocate(ieg(npi*mnei),stat=stat);
    if(stat/=0) call parallel_abort('AQUIRE_HGRID: ieg allocation failure')
    neg=0; ieg=0;
    iegloop: do iegb=1,ne_global
      if(iegrpv(iegb)==myrank) cycle iegloop !skip resident elements
      do j=1,i34gb(iegb)
        ipgb=nmgb(j,iegb)
        if(ipgl(ipgb)%rank==myrank) then !local interface node => ghost element
          neg=neg+1
!         check bound
          if(neg>npi*mnei) then
            write(errmsg,*)'Overflow in ieg:',npi*mnei
            call parallel_abort(errmsg)
          endif
          ieg(neg)=iegb
          cycle iegloop
        endif
      enddo !j
    enddo iegloop

    ! Build ordered list of ghost nodes/sides
    if(allocated(ipg)) deallocate(ipg); allocate(ipg(neg*4),stat=stat);
    if(stat/=0) call parallel_abort('AQUIRE_HGRID: ipg allocation failure')
    if(allocated(isg)) deallocate(isg); allocate(isg(neg*4),stat=stat);
    if(stat/=0) call parallel_abort('AQUIRE_HGRID: isg allocation failure')

    npg=0; ipg=0;
    nsg=0; isg=0;
    do i=1,neg
      iegb=ieg(i)
      do j=1,i34gb(iegb)
        !Add to node list if not local interface node
        ipgb=nmgb(j,iegb)
        if(ipgl(ipgb)%rank/=myrank) then
          found=.false.
          kn: do k=1,npg !check if already counted
            if(ipg(k)==ipgb) then
              found=.true.
              exit kn
            endif
          enddo kn
          if(.not.found) then !add to list
            npg=npg+1
            if(npg>neg*4) then
              write(errmsg,*)'Overflow in ipg:',neg*4
              call parallel_abort(errmsg)
            endif
            ipg(npg)=ipgb
          endif !.not.found
        endif 

        !Add to side list if not local interface side
        isgb=jsgb(j,iegb)
        if(isgl(isgb)%rank/=myrank) then
          found=.false.
          ks: do k=1,nsg !check if already counted
            if(isg(k)==isgb) then
              found=.true.
              exit ks
            endif
          enddo ks
          if(.not.found) then !add to list
            nsg=nsg+1
            if(nsg>neg*4) then
              write(errmsg,*)'Overflow in isg:',neg*4
              call parallel_abort(errmsg)
            endif
            isg(nsg)=isgb
          endif !.not.found
        endif !not in myrank
      enddo !j=1,
    enddo !i=1,neg

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Handle case of nproc==1
    else

      neg=0
      npg=0
      nsg=0

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  endif !nproc >1

  !-----------------------------------------------------------------------------
  ! Build local-to-global element, node and side index tables for augmented
  ! subdomain. Adjust global-to-local element, node and side index tables for
  ! augmented subdomain.
  !-----------------------------------------------------------------------------

  ! Set size of augmented subdomain (elements & nodes)
  nea=ne+neg
  npa=np+npg
  nsa=ns+nsg

  ! Allocate and build augmented local-to-global element index tables
  if(allocated(ielg)) deallocate(ielg); allocate(ielg(nea),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: ielg allocation failure')
  do iegb=1,ne_global
    if(iegl(iegb)%rank==myrank) ielg(iegl(iegb)%id)=iegb
  enddo
  if(nproc>1) ielg(ne+1:nea)=ieg(1:neg)

  ! Adjust global-to-local element index table to account for augmented subdomain
! iegl%rank is myrank if inside the aug. domain.
  do i=1,neg
    iegb=ieg(i)
    allocate(llp,stat=stat)
    if(stat/=0) call parallel_abort('AQUIRE_HGRID: aug-element allocation failure')
    llp=llist_type(iegl(iegb)%rank,iegl(iegb)%id,iegl(iegb)%next)
!  add to the front of the list
    iegl(iegb)%rank=myrank
    iegl(iegb)%id=ne+i
    iegl(iegb)%next=>llp
  enddo !i

  ! Allocate and build augmented local-to-global node index tables
  if(allocated(iplg)) deallocate(iplg); allocate(iplg(npa),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: iplg allocation failure')
  do ipgb=1,np_global
    if(ipgl(ipgb)%rank==myrank) iplg(ipgl(ipgb)%id)=ipgb
  enddo
  if(nproc>1) iplg(np+1:npa)=ipg(1:npg)

  ! Adjust global-to-local node index table to account for augmented subdomain
!  ipgl%rank is still myrank (if in the aug. domain), but the rest may not be in ascending order!
  do i=1,npg
    ipgb=ipg(i)
    allocate(llp,stat=stat)
    if(stat/=0) call parallel_abort('AQUIRE_HGRID: aug-node allocation failure')
    llp=llist_type(ipgl(ipgb)%rank,ipgl(ipgb)%id,ipgl(ipgb)%next)
    ipgl(ipgb)%rank=myrank
    ipgl(ipgb)%id=np+i
    ipgl(ipgb)%next=>llp
  enddo !i

  ! Allocate and build augmented local-to-global side index tables
  if(allocated(islg)) deallocate(islg); allocate(islg(nsa),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: islg allocation failure')
  do isgb=1,ns_global
    if(isgl(isgb)%rank==myrank) islg(isgl(isgb)%id)=isgb
  enddo
  if(nproc>1) islg(ns+1:nsa)=isg(1:nsg)

  ! Adjust global-to-local side index table to account for augmented subdomain
  do i=1,nsg
    isgb=isg(i)
    allocate(llp,stat=stat)
    if(stat/=0) call parallel_abort('AQUIRE_HGRID: aug-side allocation failure')
    llp=llist_type(isgl(isgb)%rank,isgl(isgb)%id,isgl(isgb)%next)
    isgl(isgb)%rank=myrank
    isgl(isgb)%id=ns+i
    isgl(isgb)%next=>llp
  enddo !i

  ! Deallocated ghost list arrays
  if(allocated(ieg)) deallocate(ieg)
  if(allocated(ipg)) deallocate(ipg)
  if(allocated(isg)) deallocate(isg)

  !-----------------------------------------------------------------------------
  ! Build geometry tables for augmented subdomain using local indices.
  !-----------------------------------------------------------------------------

  ! Allocate and build augmented element-node, element-side and
  ! element-side-element tables using local indices.
  ! NOTE: The negative global element index is used for elements that
  ! lie outside the augmented subdomain in the element-side-element table.
  if(allocated(i34)) deallocate(i34); allocate(i34(nea),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: i34 allocation failure')
  if(allocated(elnode)) deallocate(elnode); allocate(elnode(4,nea),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: elnode allocation failure')
  if(allocated(ic3)) deallocate(ic3); allocate(ic3(4,nea),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: ic3 allocation failure')
  if(allocated(elside)) deallocate(elside); allocate(elside(4,nea),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: elside (element-side table) allocation failure')

  do ie=1,nea
    iegb=ielg(ie)
    i34(ie)=i34gb(iegb)
    do j=1,i34(ie)
      elnode(j,ie)=ipgl(nmgb(j,iegb))%id
      elside(j,ie)=isgl(jsgb(j,iegb))%id
      jegb=ic3gb(j,iegb) !>=0
      if(jegb/=0.and.iegl(max(1,jegb))%rank==myrank) then
        ic3(j,ie)=iegl(jegb)%id
      else
        ic3(j,ie)=-jegb !including bnd side
      endif
    enddo !j
  enddo !ie

  ! Allocate and build augmented node-element  & node-node tables using local indices.
  ! NOTE: The negative of the global element/node index is used for elements/nodes that lie
  ! outside the augmented subdomain in the table.
  if(allocated(nne)) deallocate(nne); allocate(nne(npa),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: nne allocation failure')
  if(allocated(indel)) deallocate(indel); allocate(indel(mnei,npa),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: indel allocation failure')
  if(allocated(iself)) deallocate(iself); allocate(iself(mnei,npa),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: iself allocation failure')
  do ip=1,npa
    ipgb=iplg(ip)
    nne(ip)=nnegb(ipgb)
    do k=1,nnegb(ipgb)
      iegb=inegb(ipgb,k)
      if(iegl(iegb)%rank==myrank) then
        indel(k,ip)=iegl(iegb)%id
      else
        indel(k,ip)=-iegb
      endif
!      iself(k,ip)=iselfgb(k,ipgb)
      iself(k,ip)=0 !always follow global
      do l=1,i34gb(iegb)
        if(nmgb(l,iegb)==ipgb) then
          iself(k,ip)=l; exit
        endif
      enddo !l
      if(iself(k,ip)==0) call parallel_abort('AQUIRE_HGRID: bomb (3)')
    enddo !k
  enddo !ip

  if(allocated(nnp)) deallocate(nnp); allocate(nnp(npa),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: nnp allocation failure')
  if(allocated(indnd)) deallocate(indnd); allocate(indnd(mnei_p,npa),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: indnd allocation failure')

  nnp=0 
  do i=1,npa
    ipgb=iplg(i)
    do j=1,nnegb(ipgb)
      ie=inegb(ipgb,j)
      ii=iself(j,i) !self index

      do k=1,i34gb(ie)-2 !# of new nodes (counter-clockwise order)
        nnp(i)=nnp(i)+1
        if(nnp(i)>mnei_p) then
          write(errmsg,*)'Too many neighbor nodes (0):',i,j,k,ie,mnei_p
          call parallel_abort(errmsg)
        endif
        !Use indnd to temporarily store global node #
        indnd(nnp(i),i)=nmgb(nxq(k,ii,i34gb(ie)),ie)
      enddo !k=1,i34gb(ie)-2

!     Check last element
      if(j==nnegb(ipgb)) then 
        intvalue=nmgb(nxq(i34gb(ie)-1,ii,i34gb(ie)),ie) !last node (global) to check
        if(isbnd_global(ipgb)==0) then !complete internal ball
!         Check completeness
          if(intvalue/=indnd(1,i)) then
            write(errmsg,*)'Broken ball:',i
            call parallel_abort(errmsg)
          endif
        else !one more node for bnd ball
          nnp(i)=nnp(i)+1
          if(nnp(i)>mnei_p) then
            write(errmsg,*)'Too many neighbor nodes:',i,mnei_p
            call parallel_abort(errmsg)
          endif
          indnd(nnp(i),i)=intvalue !global
        endif !isbnd_global
      endif !j==nnegb(ipgb)
    enddo !j=1,nnegb(i)

    !Convert indnd to local node #
    do j=1,nnp(i)
      new=indnd(j,i)
      if(ipgl(new)%rank==myrank) then
        indnd(j,i)=ipgl(new)%id
      else
        indnd(j,i)=-new
      endif
    enddo !j
  enddo !i=1,npa

  !-----------------------------------------------------------------------------
  ! Build 2-tier ghost and augmented domain for higher-order schemes
  ! Note that not all arrays as in 1-tier ghost are defined
  ! Divide into different cases in which different definitions of 2-tier ghost
  ! are used in order to maximize efficiency.
  ! The 2-tier augmented domain is constructed by appending 2-tier ghosts
  ! (elem/node/side) to the 1-tier augmented domain, and thus the original local order
  ! is not altered.
  !-----------------------------------------------------------------------------
  i2tier_case=2
  if(full_aquire) then
    if(i2tier_case==1) then
      !-----------------------------------------------------------------------------
      ! Usual definition of 2-tier ghosts: an elem. that (1) surrounds a (1-tier) ghost
      ! node and (2) is not in the (1-tier) aug. domain. All of its nodes/sides that are
      ! not in the (1-tier) aug. domain are 2-tier ghosts.
      ! The global-to-local arrays (iegl2 etc) are not linked list but simple 2D arrays.
      !-----------------------------------------------------------------------------
      if(nproc>1) then
        if(allocated(ieg2)) deallocate(ieg2); allocate(ieg2(npg*mnei),stat=stat)
        if(allocated(ipg2)) deallocate(ipg2); allocate(ipg2(npg*mnei*4),stat=stat)
        if(allocated(isg2)) deallocate(isg2); allocate(isg2(npg*mnei*4),stat=stat)
        if(stat/=0) call parallel_abort('AQUIRE_HGRID: ieg2 allocation failure')

        !Count # of 2-tier ghost elem. (neg2)
        neg2=0; npg2=0; nsg2=0
        do ip=np+1,npa !ghost nodes
          do l=1,nne(ip) 
            k=indel(l,ip)
            if(k==0) call parallel_abort('aquire_hgrid: bomb (1)')
            if(k<0) then !outside aug.
              k=iabs(k) !global elem. #
              !Check against list
              found=.false.
              do j=1,neg2
                if(k==ieg2(j)) then
                  found=.true.; exit
                endif
              enddo !j

              if(.not.found) then
                neg2=neg2+1
                if(neg2>npg*mnei) call parallel_abort('AQUIRE_HGRID: overflow (9)')
                ieg2(neg2)=k !global

                do jp=1,i34gb(k) !node/side
                  ii=nmgb(jp,k) !global node
                  if(ipgl(ii)%rank/=myrank) then
                    found1=.false.
                    do j=1,npg2
                      if(ii==ipg2(j)) then
                        found1=.true.; exit
                      endif
                    enddo !j
                    if(.not.found1) then
                      npg2=npg2+1
                      if(npg2>4*npg*mnei) call parallel_abort('AQUIRE_HGRID: overflow (8)')
                      ipg2(npg2)=ii
                    endif
                  endif !ipgl

                  ii=jsgb(jp,k) !global side
                  if(isgl(ii)%rank/=myrank) then
                    found2=.false.
                    do j=1,nsg2
                      if(ii==isg2(j)) then
                        found2=.true.; exit
                      endif
                    enddo !j
                    if(.not.found2) then
                      nsg2=nsg2+1
                      if(nsg2>4*npg*mnei) call parallel_abort('AQUIRE_HGRID: overflow (5)')
                      isg2(nsg2)=ii !global
                    endif !new
                  endif !ipgl

                enddo !jp
              endif !not found
            endif !k
          enddo !l=1,nne(ip)
        enddo !ip=np+1,npa

      else !nproc==1
        neg2=0; npg2=0; nsg2=0
      endif !nproc

    else if(i2tier_case==2) then
      !-----------------------------------------------------------------------------
      ! TVD way of defining 2-tier ghosts (greedy): starting from a ghost elem,
      ! and if at least 1 of its sides is on the interface, then check its
      ! surrounding elem. and add those outside aug. domain as 2-tier ghost
      ! elem. No 2-tier nodes/sides defined. 
      !-----------------------------------------------------------------------------
      if(nproc>1) then
        if(allocated(ieg2)) deallocate(ieg2); allocate(ieg2(npg*mnei),stat=stat)
        if(allocated(ipg2)) deallocate(ipg2); allocate(ipg2(npg*mnei*4),stat=stat)
        if(allocated(isg2)) deallocate(isg2); allocate(isg2(npg*mnei*4),stat=stat)
        if(stat/=0) call parallel_abort('AQUIRE_HGRID: ieg2 allocation failure(2)')

        !Count # of 2-tier ghost elem. (neg2)
        neg2=0; npg2=0; nsg2=0
        do i=ne+1,nea
          found=.false.
          do j=1,i34(i)
            if(elside(j,i)<=ns) then !must be interface
              found=.true.; exit
            endif
          enddo !j
          if(found) then
            do j=1,i34(i)
              ie=ic3(j,i)
              if(ie<0) then !outside; check against list
                found=.false.
                ie=iabs(ie) !global index
                do l=1,neg2
                  if(ie==ieg2(l)) then
                    found=.true.; exit
                  endif
                enddo !l
                if(.not.found) then
                  neg2=neg2+1
                  if(neg2>npg*mnei) call parallel_abort('AQUIRE_HGRID: overflow (9.2)')
                  ieg2(neg2)=ie !global
                endif
              endif
            enddo !j
          endif !found
        enddo !i

        !Check pure tri
        if(.not.lhas_quad.and.neg2/=0) then
          write(errmsg,*)'Pure tri should not have 2-tier:',neg2,ieg2(1:neg2)
          call parallel_abort(errmsg)
        endif
      else !nproc==1
        neg2=0; npg2=0; nsg2=0
      endif !nproc
    else
      call parallel_abort('AQUIRE_HGRID: Unknown i2tier_case')
    endif !i2tier_case

    nea2=nea+neg2
    npa2=npa+npg2
    nsa2=nsa+nsg2
    if(allocated(ielg2)) deallocate(ielg2); allocate(ielg2(nea2),stat=stat)
    if(allocated(iplg2)) deallocate(iplg2); allocate(iplg2(npa2),stat=stat)
    if(allocated(islg2)) deallocate(islg2); allocate(islg2(nsa2),stat=stat)
    if(stat/=0) call parallel_abort('AQUIRE_HGRID: ielg2 allocation failure')

    !Update local_to_global and global_to_local arrays (these are not linked lists)
    ielg2(1:nea)=ielg(1:nea) !init
    iplg2(1:npa)=iplg(1:npa)
    islg2(1:nsa)=islg(1:nsa)

    do i=1,ne_global
      iegl2(1,i)=iegl(i)%rank !init
      iegl2(2,i)=iegl(i)%id !init
    enddo !i
    do i=1,np_global
      ipgl2(1,i)=ipgl(i)%rank
      ipgl2(2,i)=ipgl(i)%id
    enddo !i
    do i=1,ns_global
      isgl2(1,i)=isgl(i)%rank
      isgl2(2,i)=isgl(i)%id 
    enddo !i
    
    if(nproc>1) then
      ielg2(nea+1:nea2)=ieg2(1:neg2)
      iplg2(npa+1:npa2)=ipg2(1:npg2)
      islg2(nsa+1:nsa2)=isg2(1:nsg2)
      do ii=1,neg2
        iegb=ieg2(ii)
        iegl2(1,iegb)=myrank
        iegl2(2,iegb)=nea+ii !local index
      enddo !iegb

      do ii=1,npg2
        ipgb=ipg2(ii)
        ipgl2(1,ipgb)=myrank
        ipgl2(2,ipgb)=npa+ii !local index
      enddo !iegb

      do ii=1,nsg2
        isgb=isg2(ii)
        isgl2(1,isgb)=myrank
        isgl2(2,isgb)=nsa+ii !local index
      enddo !iegb

      !Done with temp arrays
      if(allocated(ieg2)) deallocate(ieg2)
      if(allocated(ipg2)) deallocate(ipg2)
      if(allocated(isg2)) deallocate(isg2)
    endif !nproc>1

    ! Output augmented subdomain size to mirror and/or screen
    allocate(isbuf(15)); allocate(irbuf(15*nproc));
    isbuf(1)=nea; isbuf(2)=ne; isbuf(3)=neg; isbuf(4)=nea2; isbuf(5)=neg2;
    isbuf(6)=npa; isbuf(7)=np; isbuf(8)=npg; isbuf(9)=npa2; isbuf(10)=npg2;
    isbuf(11)=nsa; isbuf(12)=ns; isbuf(13)=nsg; isbuf(14)=nsa2; isbuf(15)=nsg2;
    call mpi_gather(isbuf,15,itype,irbuf,15,itype,0,comm,ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('AQUIRE_HGRID: gather subdomain size',ierr)
    if(myrank==0) then
      write(16,'(/a)') '**********Augmented Subdomain Sizes**********'
      write(16,'(16a)') ' rank', &
      &'     nea','      ne','     neg','     nea2','     neg2', &
      &'     npa','      np','     npg','     npa2','     npg2', &
      &'     nsa','      ns','     nsg','     nsa2','     nsg2'
      do i=0,nproc-1
        write(16,'(i5,15i8)') i, &
        &irbuf(15*i+1),irbuf(15*i+2),irbuf(15*i+3),irbuf(15*i+4),irbuf(15*i+5), &
        &irbuf(15*i+6),irbuf(15*i+7),irbuf(15*i+8),irbuf(15*i+9),irbuf(15*i+10), &
        &irbuf(15*i+11),irbuf(15*i+12),irbuf(15*i+13),irbuf(15*i+14),irbuf(15*i+15)
      enddo
    endif !myrank==0
    call parallel_barrier
    deallocate(isbuf,irbuf)

  endif !full_aquire

  !-----------------------------------------------------------------------------
  ! Deallocate global grid geometry arrays
  !-----------------------------------------------------------------------------
  if(allocated(i34gb)) deallocate(i34gb)
  if(allocated(nmgb)) deallocate(nmgb)
  if(allocated(ic3gb)) deallocate(ic3gb)
  if(allocated(jsgb)) deallocate(jsgb)
  if(allocated(nnegb)) deallocate(nnegb)
  if(allocated(inegb)) deallocate(inegb)
  if(allocated(isbnd_global)) deallocate(isbnd_global)
!  if(allocated(nnpgb)) deallocate(nnpgb)
!  if(allocated(inpgb)) deallocate(inpgb)
!  if(allocated(iselfgb)) deallocate(iselfgb)

  !-----------------------------------------------------------------------------
  ! Allocate and read node coordinates and depths from hgrid.gr3 for augmented
  ! subdomain
  !-----------------------------------------------------------------------------

  if(allocated(xnd)) deallocate(xnd); allocate(xnd(npa),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: x allocation failure')
  if(allocated(ynd)) deallocate(ynd); allocate(ynd(npa),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: y allocation failure')
  if(allocated(znd)) deallocate(znd); allocate(znd(npa),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: z allocation failure')
  znd=0 !for ics=1
  if(allocated(dp)) deallocate(dp); allocate(dp(npa),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: dp allocation failure')
  if(allocated(xlon)) deallocate(xlon); allocate(xlon(npa),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: xlon allocation failure')
  if(allocated(ylat)) deallocate(ylat); allocate(ylat(npa),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: ylat allocation failure')

  open(14,file=in_dir(1:len_in_dir)//'hgrid.gr3',status='old',iostat=stat)
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: open(14) failure')
  read(14,*); read(14,*);
  do i=1,np_global
    read(14,*) ipgb,xtmp,ytmp,dptmp
    node=>ipgl(ipgb)
    if(node%rank==myrank) then
      ii=node%id
      if(ics==1) then
        xnd(ii)=xtmp
        ynd(ii)=ytmp
        dp(ii)=dptmp
      else !lat/lon
        xlon(ii)=xtmp*deg2rad
        ylat(ii)=ytmp*deg2rad
        dp(ii)=dptmp
        lreadll=.true. !flag to indicate lat/lon already read in
        !global coordi.
        xnd(ii)=rearth_eq*cos(ylat(ii))*cos(xlon(ii))
        ynd(ii)=rearth_eq*cos(ylat(ii))*sin(xlon(ii))
        znd(ii)=rearth_pole*sin(ylat(ii))
      endif !ics
!      if(dp(node%id).le.0) idrynode(node%id)=1
    endif !node%rank
  enddo !i
  close(14)

  !-----------------------------------------------------------------------------
  ! Only do the rest if full_aquire
  !-----------------------------------------------------------------------------
  if(full_aquire) then
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Compute element properties: areas, equivalent radii etc
  ! Compute transformation tensor for element frame eframe(i,j,ie) for ics=2
  ! where j is the axis id, i is the component id, ie is the local element id
  ! Compute xel, yel - coord. in element frame for ics=2 (for ics=1 they 
  ! are copied from xnd,ynd)
  ! WARINING: be careful when using average to compute lat/lon at center around dateline
  if(allocated(xctr)) deallocate(xctr); allocate(xctr(nea),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: xctr allocation failure')
  if(allocated(yctr)) deallocate(yctr); allocate(yctr(nea),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: yctr allocation failure')
  if(allocated(zctr)) deallocate(zctr); allocate(zctr(nea),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: zctr allocation failure')
  if(allocated(area)) deallocate(area); allocate(area(nea),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: area allocation failure')
  if(allocated(radiel)) deallocate(radiel); allocate(radiel(nea),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: radiel allocation failure')
  zctr=0._rkind !for ics=1
  if(allocated(dpe)) deallocate(dpe); allocate(dpe(nea),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: dpe allocation failure')
  if(allocated(eframe)) deallocate(eframe); allocate(eframe(3,3,nea),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: eframe allocation failure')
  if(allocated(xel)) deallocate(xel); allocate(xel(4,nea),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: xel allocation failure')
  if(allocated(yel)) deallocate(yel); allocate(yel(4,nea),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: yel allocation failure')

  eframe=0._rkind !for ics=1
  thetan=-1.d10 !max. dot product for checking only
  found=.false. !check quad areas
  do ie=1,nea
    xctr(ie)=0._rkind
    yctr(ie)=0._rkind
    dpe(ie)=real(1.d10,rkind)
    do j=1,i34(ie)
      xctr(ie)=xctr(ie)+xnd(elnode(j,ie))/real(i34(ie),rkind)
      yctr(ie)=yctr(ie)+ynd(elnode(j,ie))/real(i34(ie),rkind)
      !zctr initialized
      if(ics==2) zctr(ie)=zctr(ie)+znd(elnode(j,ie))/real(i34(ie),rkind)
      if(dp(elnode(j,ie))<dpe(ie)) dpe(ie)=dp(elnode(j,ie))
    enddo !j

    if(ics==1) then
      xel(1:i34(ie),ie)=xnd(elnode(1:i34(ie),ie))
      yel(1:i34(ie),ie)=ynd(elnode(1:i34(ie),ie))
    else !lat/lon
      call compute_ll(xctr(ie),yctr(ie),zctr(ie),ar1,ar2)
      !local ll frame
      !1: zonal axis; 2: meridional axis; 3: outward radial
      eframe(1,1,ie)=-sin(ar1)
      eframe(2,1,ie)=cos(ar1)
      eframe(3,1,ie)=0._rkind
      eframe(1,2,ie)=-cos(ar1)*sin(ar2)
      eframe(2,2,ie)=-sin(ar1)*sin(ar2)
      eframe(3,2,ie)=cos(ar2)
      call cross_product(eframe(1,1,ie),eframe(2,1,ie),eframe(3,1,ie), &
     &                   eframe(1,2,ie),eframe(2,2,ie),eframe(3,2,ie), &
     &                   eframe(1,3,ie),eframe(2,3,ie),eframe(3,3,ie))

      egb1=eframe(1,3,ie)*xctr(ie)+eframe(2,3,ie)*yctr(ie)+eframe(3,3,ie)*zctr(ie) !dot product
      if(egb1<=0) then
        write(errmsg,*)'AQUIRE_HGRID: orientation wrong:',ielg(ie),egb1,&
     &xlon(elnode(1:i34(ie),ie)),ylat(elnode(1:i34(ie),ie)),&
     &xnd(elnode(1:i34(ie),ie)),ynd(elnode(1:i34(ie),ie)),znd(elnode(1:i34(ie),ie))
        call parallel_abort(errmsg)
      endif

      !Check
      xtmp=dot_product(eframe(1:3,1,ie),eframe(1:3,3,ie))
      ytmp=dot_product(eframe(1:3,1,ie),eframe(1:3,2,ie))
      dptmp=dot_product(eframe(1:3,3,ie),eframe(1:3,2,ie))
      ar1=max(abs(xtmp),abs(ytmp),abs(dptmp))
      if(ar1>real(1.e-7,rkind)) then
        write(errmsg,*)'AQUIRE_HGRID: axes wrong',ielg(ie),xtmp,ytmp,dptmp
        call parallel_abort(errmsg)
      endif
      if(ar1>thetan) thetan=ar1

      !Compute local x,y coord.
      do j=1,i34(ie) !nodes
        nn=elnode(j,ie) 
        xel(j,ie)=(xnd(nn)-xctr(ie))*eframe(1,1,ie)+(ynd(nn)-yctr(ie))*eframe(2,1,ie)+ &
     &(znd(nn)-zctr(ie))*eframe(3,1,ie)
        yel(j,ie)=(xnd(nn)-xctr(ie))*eframe(1,2,ie)+(ynd(nn)-yctr(ie))*eframe(2,2,ie)+ &
     &(znd(nn)-zctr(ie))*eframe(3,2,ie)
      enddo !j
    endif !ics

    area(ie)=signa(xel(1,ie),xel(2,ie),xel(3,ie),yel(1,ie),yel(2,ie),yel(3,ie))
    if(area(ie)<=0._rkind) then
      found=.true.
      write(12,'(a,2i8)') 'AQUIRE_HGRID: negative area at',ielg(ie)
    endif 

    if(i34(ie)==4) then
      ar1=signa(xel(1,ie),xel(3,ie),xel(4,ie),yel(1,ie),yel(3,ie),yel(4,ie))
      area(ie)=area(ie)+ar1
      !Also check other diagonal
      ar2=signa(xel(1,ie),xel(2,ie),xel(4,ie),yel(1,ie),yel(2,ie),yel(4,ie))
      ar3=signa(xel(2,ie),xel(3,ie),xel(4,ie),yel(2,ie),yel(3,ie),yel(4,ie))
      if(min(ar1,ar2,ar3)<=0._rkind) then
        found=.true.
        write(12,*) 'AQUIRE_HGRID: concave quad at ',ielg(ie)
      endif
    endif !i34

    if(.not.found) radiel(ie)=sqrt(area(ie)/pi) !Equivalent radius
  enddo !ie=1,nea

  !Report fatal area errors
  if(found) call parallel_abort('Negative elem. areas; check nonfatal_*')

  if(ics==2) then
    call mpi_reduce(thetan,dptmp,1,rtype,MPI_MAX,0,comm,ierr)
    if(myrank==0) then
      write(16,*)'Max. dot product of 3 axes=',real(dptmp) !thetan
    endif
  endif

  ! Allocate side data arrays for augmented subdomain
  if(allocated(isdel)) deallocate(isdel); allocate(isdel(2,nsa),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: is allocation failure')
  if(allocated(isidenode)) deallocate(isidenode); allocate(isidenode(2,nsa),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: isidenode allocation failure')
  if(allocated(xcj)) deallocate(xcj); allocate(xcj(nsa),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: xcj allocation failure')
  if(allocated(ycj)) deallocate(ycj); allocate(ycj(nsa),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: ycj allocation failure')
  if(allocated(zcj)) deallocate(zcj); allocate(zcj(nsa),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: zcj allocation failure')
  zcj=0._rkind !for ics=1
  if(allocated(dps)) deallocate(dps); allocate(dps(nsa),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: dps allocation failure')
  if(allocated(distj)) deallocate(distj); allocate(distj(nsa),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: distj allocation failure')

  ! Build side data for augmented subdomain
  do ie=1,nea
    iegb=ielg(ie)
    do j=1,i34(ie)
      jsj=elside(j,ie)
      n1=elnode(nxq(1,j,i34(ie)),ie)
      n2=elnode(nxq(2,j,i34(ie)),ie)
      if(ic3(j,ie)==0.or.(ic3(j,ie)>0.and.iegb<ielg(max(1,ic3(j,ie)))) &
      &.or.(ic3(j,ie)<0.and.iegb<iabs(ic3(j,ie)))) then !new local side
        isdel(1,jsj)=ie
        isdel(2,jsj)=ic3(j,ie)
        isidenode(1,jsj)=n1
        isidenode(2,jsj)=n2
      else !old augmented subdomain boundary side (reverse orientation)
        isdel(1,jsj)=ic3(j,ie)
        isdel(2,jsj)=ie
        isidenode(1,jsj)=n2
        isidenode(2,jsj)=n1
      endif !ic3(j,ie)==0....
      xcj(jsj)=(xnd(n1)+xnd(n2))/2._rkind
      ycj(jsj)=(ynd(n1)+ynd(n2))/2._rkind
      if(ics==2) zcj(jsj)=(znd(n1)+znd(n2))/2._rkind
      dps(jsj)=(dp(n1)+dp(n2))/2._rkind
      distj(jsj)=sqrt((xnd(n2)-xnd(n1))**2._rkind+(ynd(n2)-ynd(n1))**2._rkind+(znd(n2)-znd(n1))**2._rkind)
      if(distj(jsj)==0._rkind) then
        write(errmsg,*) 'AQUIRE_HGRID: Zero side',jsj
        call parallel_abort(errmsg)
      endif
    enddo !j=1,
  enddo !ie=1,nea

  ! Allocate and compute sign (for outer normal vector) associated
  ! with each side of an element
  if(allocated(ssign)) deallocate(ssign); allocate(ssign(4,nea),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: ssign allocation failure')
  do ie=1,nea
    iegb=ielg(ie)
    do j=1,i34(ie)
      jsj=elside(j,ie)
      je=isdel(1,jsj)
      if(je==0) then !impossible
        write(errmsg,*)'First element empty:',ie,iegb,j
        call parallel_abort(errmsg)
      else if(je>0) then !resident
        jegb=ielg(je)
      else !je<0
        jegb=-je
      endif
      if(iegb==jegb) then
        ssign(j,ie)=1._rkind
      else
        ssign(j,ie)=-1._rkind
      endif    
    enddo !j
  enddo !ie

  ! Allocate and compute side frame tensor for ics=1 or 2
  ! sframe(i,j,isd): where j is the axis id, i is the component id, isd is the local side id
  ! For ics=1, only sframe(1:2,1:2,isd) are used
  if(allocated(sframe)) deallocate(sframe); allocate(sframe(3,3,nsa),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: sframe allocation failure')

  thetan=-1.d10 !max. deviation between ze and zs axes
  realvalue=-1.d10 !max. dot product of zs and ys axes
  sframe=0._rkind !for ics=1
  do j=1,nsa
    n1=isidenode(1,j)
    n2=isidenode(2,j)
    if(ics==1) then
      thetan=atan2(xnd(n1)-xnd(n2),ynd(n2)-ynd(n1))
      sframe(1,1,j)=cos(thetan) 
      sframe(2,1,j)=sin(thetan)
      sframe(1,2,j)=-sframe(2,1,j)
      sframe(2,2,j)=sframe(1,1,j)
    else !lat/lon
      !ys axis
      ar1=xnd(n2)-xnd(n1)
      ar2=ynd(n2)-ynd(n1)
      ar3=znd(n2)-znd(n1)
      ar4=sqrt(ar1*ar1+ar2*ar2+ar3*ar3)
      if(ar4==0._rkind) then
        write(errmsg,*)'AQUIRE_HGRID: 0 ys-vector',iplg(isidenode(1:2,j))
        call parallel_abort(errmsg)
      endif
      sframe(1,2,j)=ar1/ar4
      sframe(2,2,j)=ar2/ar4
      sframe(3,2,j)=ar3/ar4

      !zs axis
      ar4=sqrt(xcj(j)**2+ycj(j)**2+zcj(j)**2)
      if(ar4==0._rkind) then
        write(errmsg,*)'AQUIRE_HGRID: 0 zs-vector',iplg(isidenode(1:2,j))
        call parallel_abort(errmsg)
      endif
      sframe(1,3,j)=xcj(j)/ar4
      sframe(2,3,j)=ycj(j)/ar4
      sframe(3,3,j)=zcj(j)/ar4

      !Orthogonality between zs and ys
      egb1=abs(dot_product(sframe(1:3,2,j),sframe(1:3,3,j)))
      if(egb1>realvalue) realvalue=egb1

      !xs axis
      call cross_product(sframe(1,2,j),sframe(2,2,j),sframe(3,2,j), &
     &                   sframe(1,3,j),sframe(2,3,j),sframe(3,3,j),ar1,ar2,ar3)
      ar4=sqrt(ar1*ar1+ar2*ar2+ar3*ar3)
      if(ar4==0._rkind) then
        write(errmsg,*)'AQUIRE_HGRID: 0 xs-vector',iplg(isidenode(1:2,j))
        call parallel_abort(errmsg)
      endif
      sframe(1,1,j)=ar1/ar4
      sframe(2,1,j)=ar2/ar4
      sframe(3,1,j)=ar3/ar4

      !Check zs and ze axes (from isdel(1,j))
      if(j<=ns) then !resident
        ie=isdel(1,j)
        egb1=dot_product(sframe(1:3,3,j),eframe(1:3,3,ie))-1
        if(abs(egb1)>thetan) thetan=abs(egb1)
      endif !j<=ns

      !Debug
      !if(islg(j)==1.or.islg(j)==ns_global.or.islg(j)==1000) then
      !  xtmp=(xlon(n1)+xlon(n2))/2/pi*180
      !  ytmp=(ylat(n1)+ylat(n2))/2/pi*180
      !  write(12,*)'sample sframe:',iplg(n1),iplg(n2),xtmp,ytmp,sframe(:,:,j)
      !endif
    endif !ics
  enddo !j=1,nsa
  if(ics==2) then
    call mpi_reduce(thetan,xtmp,1,rtype,MPI_MAX,0,comm,ierr)
    call mpi_reduce(realvalue,ytmp,1,rtype,MPI_MAX,0,comm,ierr)
    if(myrank==0) then
      write(16,*)'Max. deviation between ze and zs axes=',real(xtmp) !thetan
      write(16,*)'Max. dot prod. between ys and zs axes=',real(ytmp) !realvalue
    endif
  endif !ics

  !-----------------------------------------------------------------------------
  ! Aquire open boundary segments from global for _augmented_ subdomain
  !-----------------------------------------------------------------------------

  ! Scan segments to count number of local open boundary segments and nodes
  !iopegl(0,nope_global) is the # of local fragmentations ('1' is not used yet)
  if(allocated(iopegl)) deallocate(iopegl); allocate(iopegl(0:1,nope_global),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: iopegl allocation failure')
  iopegl=0  !count of partitions of local segments
  nope=0  !local segment count
  neta=0  !local total node count
  mnond=0 !local max number of nodes per segment
  do k=1,nope_global
    j=0 !# of nodes in each local segment
    nn=nond_global(k)
    n1=iond_global(k,1)
    local1=(myrank==ipgl(n1)%rank)
    if(local1) then !start of segment is local
      iopegl(0,k)=iopegl(0,k)+1 !fragmentation of segment k in local domain
      nope=nope+1
      neta=neta+1
      j=j+1
      mnond=max(mnond,j);
    endif !local1
    do i=2,nn
      n2=iond_global(k,i)
      local2=(myrank==ipgl(n2)%rank)
      if(.not.local1.and.local2) then !segment is locally partitioned; starting new fragment
        iopegl(0,k)=iopegl(0,k)+1
        nope=nope+1; j=0;
      endif
      if(local2) then !segment continues locally
        j=j+1
        neta=neta+1
        mnond=max(mnond,j);
      endif !local2
      n1=n2
      local1=local2
    enddo !i=2,nn
  enddo !k

  ! Find maximum number partitions of local segments; i.e., global segment k is fragmented into iopegl(0,k) pieces
  mnopep=0; do k=1,nope_global; mnopep=max(mnopep,iopegl(0,k)); enddo;

  ! Allocate arrays for local open boundary segments
  if(allocated(iopegl)) deallocate(iopegl); allocate(iopegl(0:mnopep,nope_global),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: iopegl allocation failure')
  if(allocated(iopelg)) deallocate(iopelg); allocate(iopelg(nope),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: iopelg allocation failure')
  if(allocated(nond)) deallocate(nond); allocate(nond(nope),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: nond allocation failure')
  if(allocated(iond)) deallocate(iond); allocate(iond(nope,mnond),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: iond allocation failure')

  ! Initialize some for all cases
!  iopegl(0,k): # of local fragmentations of global segment k.
!  iopegl(j,k): (1<=j<=iopegl(0,k)) local open bnd segment # of jth fragmentation of  global segment k.
  nond=0; iopegl=0;
  ! Aquire local open boundary segments and nodes
  jj=nope !temporary storage for checking dimension
  nope=0
  do k=1,nope_global
    nn=nond_global(k)
    n1=iond_global(k,1)
    local1=(myrank==ipgl(n1)%rank)
    if(local1) then !start of segment is local
      nope=nope+1
      if(nope>jj) call parallel_abort('AQUIRE_HGRID: nope>jj')
      iopelg(nope)=k
      iopegl(0,k)=iopegl(0,k)+1
      if(iopegl(0,k)>mnopep) call parallel_abort('AQUIRE_HGRID: iopegl(0,k)>mnopep')
      iopegl(iopegl(0,k),k)=nope
      nond(nope)=nond(nope)+1
      if(nond(nope)>mnond) then
        write(errmsg,*)'AQUIRE_HGRID: nond(nope)>mnond',nond(nope),mnond
        call parallel_abort(errmsg)
      endif
      iond(nope,nond(nope))=ipgl(n1)%id
    endif !local1
    do i=2,nn
      n2=iond_global(k,i)
      local2=(myrank==ipgl(n2)%rank)
      if(.not.local1.and.local2) then !segment is locally partitioned
        nope=nope+1
        if(nope>jj) call parallel_abort('AQUIRE_HGRID: nope>jj (2)')
        iopelg(nope)=k
        iopegl(0,k)=iopegl(0,k)+1
        if(iopegl(0,k)>mnopep) call parallel_abort('AQUIRE_HGRID: iopegl(0,k)>mnopep (2)')
        iopegl(iopegl(0,k),k)=nope
      endif
      if(local2) then !segment continues locally
        nond(nope)=nond(nope)+1
        if(nond(nope)>mnond) call parallel_abort('AQUIRE_HGRID: nond(nope)>mnond (2)')
        iond(nope,nond(nope))=ipgl(n2)%id
      endif !local2
      n1=n2
      local1=local2
    enddo !i=2,nn
  enddo !k: global segment

  !-----------------------------------------------------------------------------
  ! Aquire land boundary segments from global for augmented subdomain
  !-----------------------------------------------------------------------------

  ! Scan segments to count number of local land boundary segments and nodes
  nland=0 !local segment count
  nvel=0  !local total node count
  mnlnd=0 !local max number of nodes per segment
  do k=1,nland_global
    j=0 !# of local nodes in each local segment
    nn=nlnd_global(k)
    n1=ilnd_global(k,1)
    local1=(myrank==ipgl(n1)%rank)
    if(local1) then
      nland=nland+1
      nvel=nvel+1
      j=j+1
      mnlnd=max(mnlnd,j);
    endif !local1
    do i=2,nn
      n2=ilnd_global(k,i)
      local2=(myrank==ipgl(n2)%rank)
      if(.not.local1.and.local2) then; nland=nland+1; j=0; endif;
      if(local2) then
        j=j+1
        nvel=nvel+1
        mnlnd=max(mnlnd,j);
      endif !local2
      n1=n2
      local1=local2
    enddo !i
  enddo !k

  ! Allocate arrays for local land boundary segments
  if(allocated(nlnd)) deallocate(nlnd); allocate(nlnd(nland),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: nlnd allocation failure')
  if(allocated(ilnd)) deallocate(ilnd); allocate(ilnd(nland,mnlnd),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: ilnd allocation failure')

  ! Initialize some for all cases
  nlnd=0; ilnd=0;

  ! Aquire local land boundary segments and nodes
  nland=0
  do k=1,nland_global
    nn=nlnd_global(k)
    n1=ilnd_global(k,1)
    local1=(myrank==ipgl(n1)%rank)
    if(local1) then
      nland=nland+1
      nlnd(nland)=nlnd(nland)+1
      ilnd(nland,nlnd(nland))=ipgl(n1)%id
    endif !local1
    do i=2,nn
      n2=ilnd_global(k,i)
      local2=(myrank==ipgl(n2)%rank)
      if(.not.local1.and.local2) nland=nland+1
      if(local2) then
        nlnd(nland)=nlnd(nland)+1
        ilnd(nland,nlnd(nland))=ipgl(n2)%id
      endif !local2
      n1=n2
      local1=local2
    enddo !i
  enddo !k

  !-----------------------------------------------------------------------------
  ! Setup various mappings and data for open boundary segment elements,
  ! nodes and sides.
  !-----------------------------------------------------------------------------

  ! Allocate and assign node-to-open-boundary-segment flags
  ! _global_ open bnd segment # if isbnd(1,ip)>0 (in this case isbnd(2,ip) may also be positive,
  !     even though isbnd(2,ip) may be outside the aug. domain); in this case, isbnd(-2:-1,ip) 
  !     are global index for the open bnd node;
  ! _global_ land bnd if isbnd(1,ip)=-1 (not on any open bnd);
  ! isbnd(1,ip)=0 if ip is internal node
  if(allocated(isbnd)) deallocate(isbnd); allocate(isbnd(-2:2,npa),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: isbnd allocation failure')
  isbnd=0
  do k=1,nope_global
    do j=1,nond_global(k)
      ipgb=iond_global(k,j)
      if(ipgl(ipgb)%rank==myrank) then
        ip=ipgl(ipgb)%id
        if(isbnd(1,ip)==0) then
          isbnd(1,ip)=k !point to global segment
          isbnd(-1,ip)=j !global index
        else if(isbnd(2,ip)==0) then
          isbnd(2,ip)=k !point to global segment
          isbnd(-2,ip)=j !global index
        else
          write(errmsg,*)'agquire_hgrid: node on more than 2 open bnds:',ipgb
          call parallel_abort(errmsg)
        endif
      endif !resident
    enddo !j
  enddo !k

  do k=1,nland_global
    do j=1,nlnd_global(k)
      ipgb=ilnd_global(k,j)
      if(ipgl(ipgb)%rank==myrank) then
        ip=ipgl(ipgb)%id
        if(isbnd(1,ip)==0) isbnd(1,ip)=-1 !overlap with open bnd
      endif
    enddo !j
  enddo !k

  ! Allocate and classify boundary sides
  ! isbs >0 if on open bnd (points to global segment #); =-1 if land bnd; =0 if internal
  if(allocated(isbs)) deallocate(isbs); allocate(isbs(nsa),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: isbs allocation failure')
  isbs=0
  do i=1,nope_global
    do j=1,nond_global(i)-1
      n1=iond_global(i,j) !global
      n2=iond_global(i,j+1)
      if(ipgl(n1)%rank/=myrank.or.ipgl(n2)%rank/=myrank) cycle

!     Both nodes are in aug. domain. However, side (n3,n4) may not be in this
!     rank -e.g., ghost zone includes n3 and n4 on 2 separate fronts
      n3=ipgl(n1)%id !local node index
      n4=ipgl(n2)%id
      do ii=1,nne(n3)
        ie=indel(ii,n3)
        if(ie>0) then !inside aug.
          k=0 !flag
          do jj=1,i34(ie)
            isd=elside(jj,ie)
            if((isidenode(1,isd)==n3.or.isidenode(2,isd)==n3).and. &
              &(isidenode(1,isd)==n4.or.isidenode(2,isd)==n4)) then
              k=isd; exit
            endif
          enddo !jj

          if(k>0) then
            if(isdel(2,k)/=0) then
              write(errmsg,*)'aquire_hgrid: impossible (1)',n1,n2,isdel(1:2,k),ielg(isdel(1:2,k)),ielg(ie),iplg(isidenode(1:2,isd))
              call parallel_abort(errmsg)
            endif
            isbs(k)=i; exit
          endif !k
        endif !ie>0
      enddo !ii=1,nne(n3)
    enddo !j
  enddo !i=1,nope_global

  ! Land bnd
  do i=1,nsa
    if(isdel(2,i)==0.and.isbs(i)==0) isbs(i)=-1
  enddo !i 

  ! Release needed buffer
  if(allocated(ibuf)) deallocate(ibuf)

  ! Output boundary info to mirror and/or screen
  if(myrank==0) then
    write(16,'(/a)') '**********Global Boundary Sizes**********'
    write(16,'(4a)') '    nope','    neta','   nland','    nvel'
    write(16,'(4i8)') nope_global,neta_global,nland_global,nvel_global
  endif

  allocate(isbuf(4),irbuf(4*nproc),stat=stat);
  if(stat/=0) call parallel_abort('AQUIRE_HGRID: isbuf allocation failure')

  isbuf(1)=nope; isbuf(2)=neta; isbuf(3)=nland; isbuf(4)=nvel;
  call mpi_gather(isbuf,4,itype,irbuf,4,itype,0,comm,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('AQUIRE_HGRID: gather subdomain bnd size',ierr)
  if(myrank==0) then
    write(16,'(/a)') '**********Augmented Subdomain Boundary Sizes**********'
    write(16,'(5a)') '    rank','    nope','    neta','   nland','    nvel'
    do i=0,nproc-1
      write(16,'(5i8)') i,irbuf(4*i+1),irbuf(4*i+2),irbuf(4*i+3),irbuf(4*i+4)
    enddo
    write(16,*)
  endif !myrank==0
  deallocate(isbuf,irbuf)
  call parallel_barrier

  ! Output centers.bp and sidecenters.bp if ipre/=0 and nproc==1
  if(nproc==1.and.ipre/=0) call write_obe

  !-----------------------------------------------------------------------------
  ! End of full_aquire if block
  !-----------------------------------------------------------------------------
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  endif !full_aquire

!-------------------------------------------------------------------------------
! Internal subroutines for aquire_hgrid
!-------------------------------------------------------------------------------
contains


subroutine swap_llrank(n,llarray)
!-------------------------------------------------------------------------------
! Place resident linked-list entries at top of list
!-------------------------------------------------------------------------------
  implicit none
  integer,intent(in) :: n
  type(llist_type),pointer :: llarray(:)
  type(llist_type),pointer :: llp
  type(llist_type) :: lltmp
  integer :: i
  do i=1,n
    if(llarray(i)%rank==myrank) cycle
    llp=>llarray(i)%next
    do
      if(.not.associated(llp)) exit
      if(llp%rank==myrank) then
        lltmp=llarray(i)
        llarray(i)%rank=llp%rank
        llarray(i)%id=llp%id
        llp%rank=lltmp%rank
        llp%id=lltmp%id
        exit
      endif
      llp=>llp%next
    enddo
  enddo
end subroutine swap_llrank


subroutine sort_llrank(n,llarray)
!-------------------------------------------------------------------------------
! Sort remaining linked-list entries in ascending order according to rank
!-------------------------------------------------------------------------------
  implicit none
  integer,intent(in) :: n
  type(llist_type),pointer :: llarray(:)
  type(llist_type),pointer :: llp
  integer :: i,j,k,t1(100),t2(100)
  do i=1,n
    if(.not.associated(llarray(i)%next)) cycle
    if(.not.associated(llarray(i)%next%next)) cycle
    k=0
    llp=>llarray(i)%next
    do
      k=k+1
      t1(k)=llp%rank
      t2(k)=llp%id
      llp=>llp%next
      if(.not.associated(llp)) exit
    enddo
    call sort(k,t1,t2)
    llp=>llarray(i)%next
    do j=1,k
      llp%rank=t1(j)
      llp%id=t2(j)
      llp=>llp%next
    enddo
  enddo
end subroutine sort_llrank


subroutine sort(n,ra,rb)
!---------------------------------------------------------------------------        
!  sorts array ra of length n into ascending order using heapsort algorithm.
!  n is input; ra is replaced on its output by its sorted rearrangement.
!  if second array is present then sort also
!  ref: numerical recipes
!--------------------------------------------------------------------------- 
  implicit none
  integer :: n, l, ir, rra, rrb, i, j
  integer :: ra(n)
  integer,optional :: rb(n)
   l = n/2 + 1
   ir = n
10 continue
   if (l.gt.1)then
     l=l-1
     rra = ra(l)
     if(present(rb)) rrb=rb(l)
   else
     rra=ra(ir)
     ra(ir)=ra(1)
     if(present(rb)) then
       rrb=rb(ir)
       rb(ir)=rb(1)
     endif
     ir=ir-1
     if (ir.eq.1) then
       ra(1)=rra
       if(present(rb)) rb(1)=rrb
       return
     endif
   endif
   i=l
   j=l+l
20 if (j.le.ir) then
     if (j.lt.ir) then
       if(ra(j).lt.ra(j+1)) j=j+1
     endif
     if (rra.lt.ra(j)) then
       ra(i)=ra(j)
       if(present(rb)) rb(i)=rb(j)
       i=j
       j=j+j
     else
       j=ir+1
     endif
     go to 20
   endif
   ra(i)=rra
   if(present(rb)) rb(i)=rrb
   go to 10
end subroutine sort


end subroutine aquire_hgrid

!===============================================================================
!===============================================================================

!dir$ attributes forceinline :: signa
function signa(x1,x2,x3,y1,y2,y3)
!-------------------------------------------------------------------------------
! Compute signed area formed by pts 1,2,3 (positive counter-clockwise)
!-------------------------------------------------------------------------------
  use schism_glbl, only : rkind,errmsg
  implicit none
  real(rkind) :: signa
  real(rkind),intent(in) :: x1,x2,x3,y1,y2,y3

  signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2._rkind
  
end function signa

!===============================================================================
!===============================================================================

subroutine dump_hgrid
!-------------------------------------------------------------------------------
! Dump horizontal grid data to processor specific formatted files
! Write local-global mapping info
!-------------------------------------------------------------------------------
  use schism_glbl
  use schism_msgp
  implicit none
  integer, parameter :: maxbuf=max(100,3)
  integer :: ie,ip,i,j,k,ngb1,ngb2,isd,isdgb,iegb1,iegb2
  integer :: ibuf1(maxbuf),ibuf2(maxbuf),ibuf3(maxbuf)
  type(llist_type),pointer :: llp
!-------------------------------------------------------------------------------

#ifdef DEBUG
  ! Dump elements
  fdb='helem_0000'
  lfdb=len_trim(fdb)
  write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
  open(10,file=out_dir(1:len_out_dir)//fdb,status='unknown')
  write(10,'(a,4i10)') '#',nea,ne,neg
  do ie=1,nea
    j=0
    llp=>iegl(ielg(ie))
    do
      j=j+1
! Check bound
      if(j>maxbuf) call parallel_abort('Increase buffer size in dump_hgrid (1)')
      ibuf1(j)=llp%rank
      ibuf2(j)=llp%id
      llp=>llp%next
      if(.not.associated(llp)) exit
    enddo !j
    if(ie<=ne) then
      write(10,'(a,2i8,4e14.6)') 'Element ',ie,ielg(ie),xctr(ie),yctr(ie),zctr(ie),dpe(ie)
    else
      write(10,'(a,2i8,4e14.6)') '# Element ',ie,ielg(ie),xctr(ie),yctr(ie),zctr(ie),dpe(ie)
    endif
    write(10,'(a,4i8)') '####NODE:  ',(iplg(elnode(k,ie)),k=1,i34(ie))
    do k=1,i34(ie)
      if(ic3(k,ie)>0) then
        ibuf3(k)=ielg(ic3(k,ie))
      elseif(ic3(k,ie)<0) then
        if(ie<=ne) then !resident must have valid ic3
          write(errmsg,*)'Resident element having wrong nbr:',ie,ielg(ie),myrank
          call parallel_abort(errmsg)
        endif
        ibuf3(k)=ic3(k,ie)
      else
        ibuf3(k)=0
      endif

!     Check validity of elnode etc
      if(elnode(k,ie)<=0.or.elside(k,ie)<=0) then
        write(errmsg,*)'Check elnode or elside:',ielg(ie),(elnode(ip,ie),elside(ip,ie),ip=1,i34(ie))
        call parallel_abort(errmsg)
      endif

    enddo !k
    write(10,'(a,4i8)') '####IC3:   ',(ibuf3(k),k=1,i34(ie))
    write(10,'(a,4i8)') '####JS:    ',(islg(elside(k,ie)),k=1,i34(ie))
    write(10,'(a,4(1x,f10.3))') '####SSIGN: ',(ssign(k,ie),k=1,i34(ie))
    write(10,'(a,1000i8)') '####PList:',(ibuf1(k),ibuf2(k),k=1,j)
  enddo !ie=1,nea

  !2-tier ghost elem
  do ie=nea+1,nea2
    iegb1=ielg2(ie)
    if(iegl(iegb1)%rank==myrank) call parallel_abort('DUMP_HGRID: 2-tier elem (0)')
    if(iegl2(1,iegb1)/=myrank) call parallel_abort('DUMP_HGRID: 2-tier elem')
    if(iegb1<1.or.iegb1>ne_global) call parallel_abort('DUMP_HGRID: 2-tier elem (2)')
    write(10,'(a,3i8)') '2-t element ',ie-nea,ie,iegb1
  enddo !ie
  close(10)

  ! Dump nodes
  fdb='hnode_0000'
  lfdb=len_trim(fdb)
  write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
  open(10,file=out_dir(1:len_out_dir)//fdb,status='unknown')
  write(10,'(a,4i10)') '#',npa,np,npg
  do ip=1,npa
    j=0
    llp=>ipgl(iplg(ip))
    do
      j=j+1
! Check bound
      if(j>maxbuf) call parallel_abort('Increase buffer size in dump_hgrid: node (1)')
      ibuf1(j)=llp%rank
      ibuf2(j)=llp%id
      llp=>llp%next
      if(.not.associated(llp)) exit
    enddo !j

   if(nnp(ip)>maxbuf.or.nne(ip)>maxbuf) call parallel_abort('Increase buffer size in dump_hgrid: node (2)')
    do k=1,nne(ip)
      if(indel(k,ip)>0) then
        ibuf3(k)=ielg(indel(k,ip))
      elseif(indel(k,ip)<0) then
        if(ip<=np) then
          write(errmsg,*)'Surrounding element outside:',indel(k,ip),iplg(ip),k
          call parallel_abort(errmsg)
        endif

        ibuf3(k)=indel(k,ip)
      else
!  Fatal error
!        ibuf3(k)=0
        write(errmsg,*)'Surrounding element not exist:',indel(k,ip),iplg(ip),k
        call parallel_abort(errmsg) 
      endif
    enddo !k
    if(ip<=np) then
      write(10,'(a,2i8,4e14.6,2i4,1000(i8,i4))') 'Node ',ip,iplg(ip),xnd(ip),ynd(ip),znd(ip),dp(ip), &
      &isbnd(-2:2,ip),nne(ip),(ibuf3(k),iself(k,ip),k=1,nne(ip))
    else
      write(10,'(a,2i8,4e14.6,2i4,1000(i8,i4))') '# Node ',ip,iplg(ip),xnd(ip),ynd(ip),znd(ip),dp(ip), &
      &isbnd(-2:2,ip),nne(ip),(ibuf3(k),iself(k,ip),k=1,nne(ip))
    endif
    write(10,'(a,1000i8)') '####PList:',(ibuf1(k),ibuf2(k),k=1,j)

    do k=1,nnp(ip)
      if(indnd(k,ip)>0) then
        ibuf3(k)=iplg(indnd(k,ip))
      elseif(indnd(k,ip)<0) then
        if(ip<=np) then !resident
          write(errmsg,*)'Surrounding node outside:',indnd(k,ip),iplg(ip),k
          call parallel_abort(errmsg)
        endif

        ibuf3(k)=indnd(k,ip)
      else
        write(errmsg,*)'Surrounding node not exist:',indnd(k,ip),iplg(ip),k
        call parallel_abort(errmsg) 
      endif
    enddo !k
    write(10,'(a,1000i8)')'Nbr nodes:',(ibuf3(k),k=1,nnp(ip))
    
  enddo !ip=1,npa

  !2-tier ghost 
  do ip=npa+1,npa2
    ngb1=iplg2(ip)
    if(ipgl(ngb1)%rank==myrank) call parallel_abort('DUMP_HGRID: 2-tier node (0)')
    if(ipgl2(1,ngb1)/=myrank) call parallel_abort('DUMP_HGRID: 2-tier node')
    if(ngb1<1.or.ngb1>np_global) call parallel_abort('DUMP_HGRID: 2-tier node (2)')
    write(10,'(a,3i8)') '2-t node:',ip-npa,ip,ngb1
  enddo !ip
  close(10)

  ! Dump sides
  fdb='hside_0000'
  lfdb=len_trim(fdb)
  write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
  open(10,file=out_dir(1:len_out_dir)//fdb,status='unknown')
  write(10,'(a,4i10)') '#',nsa,ns,nsg
  do isd=1,nsa
    isdgb=islg(isd)
    ngb1=iplg(isidenode(1,isd)); ngb2=iplg(isidenode(2,isd));
    if(isdel(1,isd)>0) then
      iegb1=ielg(isdel(1,isd)) 
    else if(isdel(1,isd)<0) then
      iegb1=isdel(1,isd) 
    else !=0
      write(errmsg,*)'isdel(1,:) =0:',ngb1,ngb2
      call parallel_abort(errmsg)
    endif
    if(isdel(2,isd)>0) then; iegb2=ielg(isdel(2,isd)); else; iegb2=isdel(2,isd); endif;
    j=0
    llp=>isgl(islg(isd))
    do
      j=j+1
      ibuf1(j)=llp%rank
      ibuf2(j)=llp%id
      llp=>llp%next
      if(.not.associated(llp)) exit
    enddo
    if(isd<=ns) then
      write(10,'(a,6i8,7e14.6,i4)') 'Side ',isd,isdgb,ngb1,ngb2,iegb1,iegb2, &
      &xcj(isd),ycj(isd),zcj(isd),dps(isd),distj(isd),sframe(1,1,isd),sframe(2,1,isd),isbs(isd)
    else
      write(10,'(a,6i8,7e14.6,i4)') '# Side', isd,isdgb,ngb1,ngb2,iegb1,iegb2, &
      &xcj(isd),ycj(isd),zcj(isd),dps(isd),distj(isd),sframe(1,1,isd),sframe(2,1,isd),isbs(isd)
    endif
    write(10,'(a,1000i8)') '####PList:',(ibuf1(k),ibuf2(k),k=1,j)
  enddo !isd

  !2-tier ghost 
  do isd=nsa+1,nsa2
    isdgb=islg2(isd)
    if(isgl(isdgb)%rank==myrank) call parallel_abort('DUMP_HGRID: 2-tier side (0)')
    if(isgl2(1,isdgb)/=myrank) call parallel_abort('DUMP_HGRID: 2-tier side')
    if(isdgb<1.or.isdgb>ns_global) call parallel_abort('DUMP_HGRID: 2-tier side (2)')
    write(10,'(a,3i8)') '2-t side:',isd-nsa,isd,isdgb
  enddo !ip
  close(10)

  ! Dump local bnd info
  fdb='bndinfo_0000'
  lfdb=len_trim(fdb)
  write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
  open(10,file=out_dir(1:len_out_dir)//fdb,status='unknown')
  write(10,'(a,i10)') 'Open bnd:',nope
  do i=1,nope
    write(10,*)'open bnd #',i,iopelg(i),(iplg(iond(i,j)),j=1,nond(i))
  enddo !i
  write(10,'(a,i10)') 'Land bnd:',nland
  do i=1,nland
    write(10,*)'land bnd #',i,(iplg(ilnd(i,j)),j=1,nlnd(i))
  enddo !i
  close(10)

#endif /*DEBUG*/

  ! Rank 0 writes global to local element info
  if(myrank==0) then
    open(32,file=out_dir(1:len_out_dir)//'global_to_local.prop',status='unknown')
    write(32,'(i8,1x,i4)')(ie,iegrpv(ie),ie=1,ne_global)
    close(32)
  endif

end subroutine dump_hgrid

!===============================================================================
!===============================================================================
! END GRID SUBROUTINES
!===============================================================================
!===============================================================================
