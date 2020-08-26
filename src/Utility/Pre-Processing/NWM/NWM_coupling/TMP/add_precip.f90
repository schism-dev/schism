!This script is used to post-process the NWM coupling files.
!It combines the sinks with adjacent sources if there is any.
!Two options are available: distance-based and neighboring element based.

!Inputs files: 
!      rain rate: (m/hour)
!
!Input files
!      hgrid.cpp        : hgrid in cpp projection
!      source_sink.in   : contains the element ID for each
!                         intersection of the NWM streams and the land boundary.
!      msource.th       : contains the salinity and temprature
!                         of the source element along the land boundary.
!                         Salinity is set to be 0, temp = -9999.
!      vsource.th       : input of the stream flows of source elements.
!      vsink.th         : input of the stream flows of sink elements.
!
!Output files
!      source_sink.in.1, msource.th.1, vsource.th.1, vsink.th.1
!
!serial:
!ifort -CB -O2 -o add_precip add_precip.f90 
!
!openmp:
!ifort -CB -O2 -qopenmp -o add_precip add_precip.F90 
!pgf90 -O2 -mp -o add_precip add_precip.F90 
!export OMP_NUM_TREADS=8
!setenv OMP_NUM_TREADS 8

      program coupling_nwm

       use omp_lib
       implicit none
      

       real(8) :: start_time,end_time
       integer,parameter::max_nei=1000
       real(8) :: tmp,distance,dist,max_dist,rain_rate,signa
       real(8), allocatable :: temp(:),xx(:),yy(:),xel(:),yel(:)
       integer, allocatable :: elnode(:,:),i34(:),indel(:,:),n_sink_source(:), &
         & i_sink_source(:,:),ntier1(:),ele_source(:),ele_sink(:),imap_ele2source(:), &
         & tier1(:,:),tier_n(:,:),ntier_n(:),nne(:),nxq(:,:,:),ic3(:,:),isbnd(:), &
         & ncount(:),nlbnd(:),lbnd(:,:),nobnd(:),obnd(:,:),i_island(:)
       integer :: mnei,inbr,i,j,k,nsource,nsink,itmp,nt,istat,nd,ntracer,n_tier_dist
       integer :: ne,np,id,id1,ie,ie_nbr,isource,kk,k0,nope,nland,nn,n1,n2,n3,n4,l, &
         & ip,je,new,ii,icount,tid,nthreads,ntmp
       real(8), allocatable :: vsource(:,:),vsink(:,:),msource(:,:,:), &
         & time_stamp(:),area(:),vsource_precip(:),vsource_no_precip(:),max_source_precip(:)

       !$omp parallel
       !tid = omp_get_thread_num()
       !if (tid==0) then
       !  nthreads= omp_get_num_threads()
       !  print*,'Number 0f threads=',nthreads
       !  allocate(ncount(0:nthreads-1),stat=istat)
       !  if (istat/=0) stop 'Failed to alloc. ncount'
       !endif
       !$omp end parallel

       !Read in cmd line inputs
       print*, 'rain rate (m/h)'
       read*, rain_rate
       !distance=1
  
       !read hgrid
       !call cpu_time(start_time)
       !start_time=omp_get_wtime();
       print*, 'reading inputs ...'

       open(14,file='hgrid.cpp',status='old')!lambert projection, NWM has shorter precision
       read(14,*)
       read(14,*)ne,np
       write(*,*)'# of elements',ne,'# of nodes',np
       allocate(xx(np),yy(np),nne(np),stat=istat)
       if(istat/=0) stop 'Failed to alloc. xx, yy, nne'
       allocate(i34(ne),elnode(4,ne),stat=istat)
       if(istat/=0) stop 'Failed to alloc. i34, elnode'
       do i=1,np
         read(14,*)j,xx(i),yy(i),tmp
       enddo     
       do i=1,ne
         read(14,*) j,i34(i),elnode(1:i34(i),i)
       enddo

       if (inbr==1) then
         allocate(isbnd(np),stat=istat)
         if (istat/=0) stop 'Failed to alloc. isbnd'
         isbnd=0
         read(14,*) nope
         read(14,*) ntmp
         allocate(obnd(ntmp,nope),nobnd(nope),stat=istat)
         if (istat/=0) stop 'Failed to alloc. obnd, nobnd'
         do k=1,nope
           read(14,*) nobnd(k)
           do i=1,nobnd(k)
             read(14,*) ip
             obnd(i,k)=ip
             isbnd(ip)=1 !open bnd
           enddo !i
         enddo !k

         read(14,*) nland
         read(14,*) ntmp
         allocate(lbnd(ntmp,nland),nlbnd(nland),i_island(nland),stat=istat)
         if (istat/=0) stop 'Failed to alloc. lbnd, nlbnd, i_island'
         do k=1,nland
           read(14,*) nlbnd(k), i_island(k)
           do i=1,nlbnd(k)
             read(14,*) ip
             lbnd(i,k)=ip
             if (isbnd(ip)==1) then !open bnd
               isbnd(ip)=-1  !open and land bnd
             else
               isbnd(ip)=1  !land bnd
             endif
           enddo !i
         enddo !k
       endif

       close(14)


!      read existing source/sink files
       open(13,file='source_sink.in',status='old') 
       read(13,*) nsource
       allocate(ele_source(nsource),stat=istat)
       if (istat/=0) stop 'Failed to alloc. ele_source'
       do i=1,nsource
          read(13,*)ele_source(i)
       enddo
       read(13,*) 
       read(13,*) nsink
       allocate(ele_sink(nsink),stat=istat)
       if (istat/=0) stop 'Failed to alloc. ele_sink'
       do i=1,nsink
          read(13,*)ele_sink(i)
       enddo
       close(13)

       nt=0
       open(15,file='vsource.th',status='old') 
       do while (.true.)
         read (15,*,end=999) 
         nt=nt+1
       enddo
       999 continue
       close(15)

       allocate(vsource(nsource,nt),time_stamp(nt),stat=istat)
       if (istat/=0) stop 'Failed to alloc. vsource'
       open(15,file='vsource.th',status='old') 
       do i=1,nt
         read (15,*) time_stamp(i),vsource(:,i)
       enddo
       close(15)


!     add precipitation to all elements
      allocate(vsource_precip(ne),vsource_no_precip(ne),max_source_precip(ne),stat=istat)
      if (istat/=0) stop 'Failed to alloc. vsource_precip'
      vsource_no_precip=0d0
      allocate(area(ne),stat=istat)
      if (istat/=0) stop 'Failed to alloc. area'
      do i=1,ne
        n1=elnode(1,i)
        n2=elnode(2,i)
        n3=elnode(3,i)
        area(i)=signa(xx(n1),xx(n2),xx(n3),yy(n1),yy(n2),yy(n3))
        if(area(i)<=0) then
          print*, 'area<=0 (0):',i,area(i),i34(i)
          stop
        endif

        if(i34(i)==3) then
        else if(i34(i)==4) then
          n4=elnode(4,i)
          area(i)=area(i)+signa(xx(n1),xx(n2),xx(n3),yy(n1),yy(n2),yy(n3))
          if(area(i)<=0) then
            print*, 'area<=0:',i,area(i)
            stop
          endif
        else
          stop 'Unknown elem'
        endif

        vsource_precip(i)=rain_rate/3600.*area(i) 
        
      enddo


!     add max previous vsource for each source
      do i=1,nsource
        ie=ele_source(i)
        vsource_precip(ie)=vsource_precip(ie)+maxval(vsource(i,:))
        vsource_no_precip(ie)=maxval(vsource(i,:))
      enddo
!      print*,vsource_no_precip(ie)

!     write source_sink.in
      open(16,file='source_sink.in.1',status='replace')
      write(16,*) ne
      do i=1,ne
        write(16,*) i
      end do

      write(16,*) 
      write(16,*) 0
      close(16)
!     find max for each source
!      do i=1,ne
!        max_source_precip(i)=maxval(vsource_precip(i,:))+vsource_precip(i)
!      enddo

!     write vsource.th file
      open(16,file='vsource.th.1',status='replace')
      write(16,'(10000000F15.3)')0.0,(vsource_precip(k),k=1,ne)
      write(16,'(10000000F15.3)')43200.0,(vsource_precip(k),k=1,ne)
      write(16,'(10000000F15.3)')86400.0,(vsource_no_precip(k),k=1,ne)
      write(16,'(10000000F15.3)')1000000.0,(vsource_no_precip(k),k=1,ne)
      close(16)

!     write msource.th file
      open(16,file='msource.th.1',status='replace')
      write(16,'(10000000F15.3)')0.0,(-9999d0,k=1,ne),(0d0,k=1,ne)
      write(16,'(10000000F15.3)')43200.0,(-9999d0,k=1,ne),(0d0,k=1,ne)
      write(16,'(10000000F15.3)')86400.0,(-9999d0,k=1,ne),(0d0,k=1,ne)
      write(16,'(10000000F15.3)')1000000.0,(-9999d0,k=1,ne),(0d0,k=1,ne)
      close(16)

!!     write vsink.th
!      open(17,file='vsink.th.1')
!      do i=1,nt
!        write(17,'(10000000F15.3)')time_stamp(i),(vsink(k,i),k=1,nsink)
!      end do



     contains
     recursive subroutine mark_bnd_neighbors(isink,ie_ctr,i_depth,n_tier_dist)
       implicit none

       integer, intent(in) :: isink, ie_ctr, i_depth,n_tier_dist

       integer :: iter,ie000,kk,ie1,ie2,nd,nd1,nd2,ie3,ie_new(2),nn,i
       logical :: inew

       if (i_depth >= n_tier_dist) then
         return
       endif

       !find tier 1 neighbors, at most 2, at least 1 (when self is quad)
       !the following searching procedure didn't stop at non-land bnd elements.
       !This is okay since non-land bnd elements do not have sources anyway
       nn=0
       do i=1,i34(ie_ctr)
         nd = elnode(i,ie_ctr)
         if (isbnd(nd)==0) cycle !internal nodes
         ie1=indel(1,nd)
         ie2=indel(nne(nd),nd)
         if (ie1==ie_ctr .and. ie2==ie_ctr) then 
           cycle !tip of a small stream
         elseif (ie1==ie_ctr .and. ie2.ne.ie_ctr) then 
           nn=nn+1
           if (nn>2) then
             print*, 'tier 1 neighbors > 2: ', isink,ie_ctr,i_depth,ie1,ie2
           endif
           ie_new(nn)=ie2
         elseif (ie1.ne.ie_ctr .and. ie2==ie_ctr) then 
           nn=nn+1
           if (nn>2) then
             print*, 'tier 1 neighbors > 2: ', isink,ie_ctr,i_depth,ie1,ie2
           endif
           ie_new(nn)=ie1
         else
           !could happen in narrow 2DV channel, where the other node is also a bnd 
           !node but on the other side of the channel
           cycle
         endif
       enddo

       !push new elements in tier n into record
       loop1: do iter=1,nn
         ie000=ie_new(iter)
         inew=.true.
         do kk=0,ntier_n(isink)
           if (ie000==tier_n(kk,isink)) then
             inew=.false. !avoid duplicated points, but its neighbor may be new
             exit
           endif
         enddo

         !register new neighbor
         if (inew) then
           ntier_n(isink)=ntier_n(isink)+1
           tier_n(ntier_n(isink),isink)=ie000

           !debug
           !write(*,*) 'ie_ctr,i_depth,ntier_n(isink),ie000:',ie_ctr,i_depth,ntier_n(isink),ie000

           if (imap_ele2source(ie000)>0) then
             n_sink_source(isink)=n_sink_source(isink)+1
             i_sink_source(n_sink_source(isink),isink)=imap_ele2source(ie000)
           endif
         endif

         call mark_bnd_neighbors(isink,ie000,i_depth+1,n_tier_dist)

       enddo loop1
     end subroutine mark_bnd_neighbors  


     recursive subroutine mark_neighbors(isink,ie_ctr,i_depth,n_tier_dist)
       implicit none

       integer, intent(in) :: isink, ie_ctr, i_depth,n_tier_dist

       integer :: iter,ie000,kk
       logical :: inew

       if (i_depth >= n_tier_dist) then
         return
       endif

       loop1: do iter=1,ntier1(ie_ctr)
         ie000=tier1(iter,ie_ctr)

         inew=.true.
         do kk=0,ntier_n(isink)
           if (ie000==tier_n(kk,isink)) then
             inew=.false. !avoid duplicated points, but its neighbor may be new
             exit
           endif
         enddo

         !register new neighbor
         if (inew) then
           ntier_n(isink)=ntier_n(isink)+1
           tier_n(ntier_n(isink),isink)=ie000

           !debug
           !write(*,*) 'ie_ctr,i_depth,ntier_n(isink),ie000:',ie_ctr,i_depth,ntier_n(isink),ie000

           if (imap_ele2source(ie000)>0) then
             n_sink_source(isink)=n_sink_source(isink)+1
             i_sink_source(n_sink_source(isink),isink)=imap_ele2source(ie000)
           endif
         endif

         call mark_neighbors(isink,ie000,i_depth+1,n_tier_dist)

       enddo loop1
     end subroutine mark_neighbors  



end

function signa(x1,x2,x3,y1,y2,y3)
!-------------------------------------------------------------------------------
! Compute signed area formed by pts 1,2,3
!-------------------------------------------------------------------------------
real*8 :: signa
real*8,intent(in) :: x1,x2,x3,y1,y2,y3

signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2d0

end function signa
