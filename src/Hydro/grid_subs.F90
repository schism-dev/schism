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
! SCHISM GRID SUBROUTINES not used by QSim
!
! include subroutine aquire_vgrid
! subroutine partition_hgrid
! include subroutine aquire_hgrid
! function signa
! subroutine dump_hgrid
!
!===============================================================================
!===============================================================================

   include 'aquire_vgrid.F90'

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
  if(myrank==0) then
    open(14,file=in_dir(1:len_in_dir)//'hgrid.gr3',status='old')
    read(14,*); read(14,*) ne_global,np_global
    close(14)
  endif
  call mpi_bcast(ne_global,1,itype,0,comm,i)
  call mpi_bcast(np_global,1,itype,0,comm,i)

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

  ! Setup initial naive partition
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
  !from aquire_vgrid: nvrt; ztot; kz; h_s
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
    if(allocated(nlev)) deallocate(nlev); allocate(nlev(np_global),stat=stat) 
    if(myrank==0) then
      open(19,file=in_dir(1:len_in_dir)//'vgrid.in',status='old')
      read(19,*); read(19,*)nvrt !not needed
!      do i=1,np_global
      read(19,*)nlev(1:np_global) !kbetmp
!      enddo !i
      close(19)
    endif !myrank
    call mpi_bcast(nlev,np_global,itype,0,comm,stat)

    do i=1,np_global
      if(ipgl(i)%rank==myrank) kbp(ipgl(i)%id)=nlev(i) !kbetmp
    enddo !i
    deallocate(nlev)
  endif !ivcor==1

#ifdef NO_PARMETIS
    !Offline paritition by reading from input similar to global_to_local.prop
    if(myrank==0) then
      open(10,file=in_dir(1:len_in_dir)//'partition.prop',status='old')
      do i=1,ne_global 
        read(10,*)j,iegrpv(i)
      enddo
!      read(10,*); read(10,*)k
      close(10)

      k=maxval(iegrpv); l=minval(iegrpv)
      if(k/=nproc-1.or.l/=0) then
        write(errmsg,*)'Offline partition: different nproc,',k,l
        call parallel_abort(errmsg)
      endif
    endif !myrank
    call mpi_bcast(iegrpv,ne_global,itype,0,comm,stat)

#else
!Use ParMETIS
  
  ! Count number of edges in dual graph (element as 'vertex')
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

  ! Use vertex (elem) and/or edge weights
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

  ! Imbalance tolerance (1: perfect balance; nproc: perfect imbalance)
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

!  ParMETIS_V3_PartGeomKway (idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, 
!    idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, idx_t *ndims, real_t *xyz, 
!    idx_t *ncon, idx_t *nparts, real_t *tpwgts, real_t *ubvec, idx_t *options, 
!    idx_t *edgecut, idx_t *part, MPI Comm *comm)

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

  ! Construct global element-to-resident-processor assignment vector (iegrpv)
  call mpi_allgatherv(part,ne,itype,iegrpv,neproc,neprocsum,itype,comm,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('partition: mpi_allgatherv',ierr)

  ! Deallocate ParMeTiS arrays
  deallocate(vtxdist,xadj,adjncy,part)
  deallocate(xyz,tpwgts,ubvec)
  if(wgtflag==2.or.wgtflag==3) deallocate(vwgt)
  if(wgtflag==1.or.wgtflag==3) deallocate(adjwgt)
#endif /*NO_PARMETIS*/

  ! Deallocate arrays
  deallocate(neproc,neprocsum)
!  if(ivcor==1) deallocate(kbp)
  if(allocated(nlev)) deallocate(nlev)
  if(allocated(kbp)) deallocate(kbp)
  deallocate(xproj,yproj)

end subroutine partition_hgrid

!===============================================================================
!===============================================================================
   include 'aquire_hgrid.F90'

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
  fdb='helem_000000'
  lfdb=len_trim(fdb)
  write(fdb(lfdb-5:lfdb),'(i6.6)') myrank
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
  fdb='hnode_000000'
  lfdb=len_trim(fdb)
  write(fdb(lfdb-5:lfdb),'(i6.6)') myrank
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
  fdb='hside_000000'
  lfdb=len_trim(fdb)
  write(fdb(lfdb-5:lfdb),'(i6.6)') myrank
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
      &xcj(isd),ycj(isd),zcj(isd),dps(isd),distj(isd),snx(isd),sny(isd),isbs(isd)
    else
      write(10,'(a,6i8,7e14.6,i4)') '# Side', isd,isdgb,ngb1,ngb2,iegb1,iegb2, &
      &xcj(isd),ycj(isd),zcj(isd),dps(isd),distj(isd),snx(isd),sny(isd),isbs(isd)
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
  fdb='bndinfo_000000'
  lfdb=len_trim(fdb)
  write(fdb(lfdb-5:lfdb),'(i6.6)') myrank
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
#ifdef USE_QSIM
    !write(32,*)ne_global
    do ie=1,ne_global
      write(32,'(i8,1x,i6,1x,i6,1x,i6)')ie,iegrpv(ie),iegl2(1,ie),iegl2(2,ie)
    enddo !ie
    !write(32,*)np_global
    do ie=1,np_global
      write(32,'(i8,1x,i6,1x,i6)')ie, ipgl(ie)%rank, ipgl(ie)%id
    enddo !ie
#else
    write(32,'(i8,1x,i6)')(ie,iegrpv(ie),ie=1,ne_global)
#endif
    !Add more info
!    write(32,*)
!    write(32,*)nproc
    close(32)
  endif

end subroutine dump_hgrid

!===============================================================================
!===============================================================================
! END GRID SUBROUTINES
!===============================================================================
!===============================================================================
