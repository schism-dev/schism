!This script is used to post-process the NWM coupling files.
!It combines the sinks with adjacent sources if there is any.
!Two options are available: distance-based and neighboring element based.

!Inputs files: 
!      inbr: 0 (distance-based); 1 (neighboring element based)
!      dist: distance if inbr; number of neighboring tiers
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
!ifort -CB -O2 -qopenmp -o post_couple_omp post_couple_omp.f90 
!pgf90 -O2 -mp -o post_couple_omp post_couple_omp.f90 
!export OMP_NUM_TREADS=8
!setenv OMP_NUM_TREADS 8

      program coupling_nwm

       use omp_lib
       implicit none
      

       real(8) :: start_time,end_time
       integer,parameter::max_nei=1000
       real(8) :: tmp,distance,dist,max_dist
       real(8), allocatable :: temp(:),xx(:),yy(:),xel(:),yel(:)
       integer, allocatable :: elnode(:,:),i34(:),indel(:,:),n_sink_source(:), &
         & i_sink_source(:,:),ntier1(:),ele_source(:),ele_sink(:),is_source_ele(:), &
         & tier1(:,:),tier_n(:,:),ntier_n(:),nne(:),nxq(:,:,:),ic3(:,:),isbnd(:), &
         & ncount(:)
       integer :: mnei,inbr,i,j,k,nsource,nsink,itmp,nt,istat,nd,ntracer,n_tier_dist
       integer :: ne,np,id,id1,ie,ie_nbr,isource,kk,k0,nope,nland,nn,n1,n2,l, &
         & ip,je,new,ii,icount,tid,nthreads
       real(8), allocatable :: vsource(:,:),vsink(:,:),msource(:,:,:), &
         & time_stamp(:)

       !$omp parallel
       tid = omp_get_thread_num()
       if (tid==0) then
         nthreads= omp_get_num_threads()
         print*,'Number 0f threads=',nthreads
         allocate(ncount(0:nthreads-1))
       endif
       !$omp end parallel

       !Read in cmd line inputs
       print*, 'Input search option (0: distance-based; 1: neighboring element based):'
       read*, inbr
       !inbr=1
       print*, 'Input search radius (if inbr=0: distance; if inbr=1: neighboring tiers):'
       if (inbr==0) then
         read*, distance
       elseif (inbr==1) then
         read*, n_tier_dist
       else
         print*, 'wrong inbr ', inbr
         stop
       endif
       !distance=1
  
       !read hgrid
       !call cpu_time(start_time)
       start_time=omp_get_wtime();
       print*, 'reading inputs ...'

       open(14,file='hgrid.cpp',status='old')!lambert projection, NWM has shorter precision
       read(14,*)
       read(14,*)ne,np
       write(*,*)'# of elements',ne,'# of nodes',np
       allocate(xx(np),yy(np),nne(np))
       allocate(i34(ne),elnode(4,ne))
       do i=1,np
         read(14,*)j,xx(i),yy(i),tmp
       enddo     
       do i=1,ne
         read(14,*) j,i34(i),elnode(1:i34(i),i)
       enddo

       if (inbr==1) then
         allocate(isbnd(np)); isbnd=0
         read(14,*) nope
         read(14,*)
         do k=1,nope
           read(14,*) nn
           do i=1,nn
             read(14,*) ip
             isbnd(ip)=1
           enddo !i
         enddo !k

         read(14,*) nland
         read(14,*)
         do k=1,nland
           read(14,*) nn
           do i=1,nn
             read(14,*) ip
             isbnd(ip)=1
           enddo !i
         enddo !k
       endif

       close(14)

!      read existing source/sink files
       open(13,file='source_sink.in',status='old') 
       read(13,*) nsource
       allocate(ele_source(nsource))
       do i=1,nsource
          read(13,*)ele_source(i)
       enddo
       read(13,*) 
       read(13,*) nsink
       allocate(ele_sink(nsink))
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

       allocate(vsource(nsource,nt),time_stamp(nt))
       open(15,file='vsource.th',status='old') 
       do i=1,nt
         read (15,*) time_stamp(i),vsource(:,i)
       enddo
       close(15)

       ntracer=2
       allocate(msource(nsource,ntracer,nt))
       open(15,file='msource.th',status='old') 
       do i=1,nt
         read (15,*) tmp,msource(:,:,i)
       enddo
       close(15)

       allocate(vsink(nsink,nt))
       open(15,file='vsink.th',status='old') 
       do i=1,nt
         read (15,*) tmp,vsink(:,i)
       enddo
       close(15)

       !call cpu_time(end_time)
       end_time=omp_get_wtime();
       print*, 'reading inputs takes ',(end_time-start_time)/60d0,' minutes'
!------done reading inputs-----------------------------------


       !call cpu_time(start_time)
       start_time=omp_get_wtime();
       print*, 'calculating geometry ...'
       if (inbr==1) then
  !      elem ball
         nne=0
         do i=1,ne
           do j=1,i34(i)
             nd=elnode(j,i)
             nne(nd)=nne(nd)+1
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

         allocate(nxq(3,4,4),ic3(4,ne))
         !re-order indel(:,nd) in counter-clockwise order
         !setup nxq (cyclic node index)
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

         !ele-side-ele table
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
               if(ii==0) stop('AQUIRE_HGRID: bomb (1)')

               if(ic3(nxq(i34(ie)-1,ii,i34(ie)),ie)==0) then
                 icount=icount+1
                 indel(1,i)=ie
               endif
             enddo !j=1,nne
             if(icount/=1) then
               write(*,*)'Illegal bnd node',i,isbnd(i),icount
               stop
             endif
           endif !bnd ball

       !   For internal balls, starting elem. is not altered
       !   Sequential search for the rest of elements
           do j=2,nne(i)
             ie=indel(j-1,i)
             ii=0 !local index
             do l=1,i34(ie)
               if(elnode(l,ie)==i) then
                 ii=l; exit
               endif
             enddo !l
             if(ii==0) stop('AQUIRE_HGRID: bomb (2)')

             new=ic3(nxq(i34(ie)-2,ii,i34(ie)),ie)
             if(new==0) then
               write(*,*)'Incomplete ball',i
               stop
             endif
             indel(j,i)=new
           enddo !j=2,nne(i)
         enddo !i=1,np

       endif

!---------------------------------------------------------------
!      find neighboring sources for each sink
!---------------------------------------------------------------
       allocate(n_sink_source(nsink),i_sink_source(max_nei,nsink))
       n_sink_source=0; i_sink_source=0
!---------------------------------------------------------------
       if (inbr==0) then !distance-based
!---------------------------------------------------------------
  !      elem coordinates
         allocate(xel(ne),yel(ne))
         xel=0d0; yel=0d0
         do i=1,ne
           xel(i)=sum(xx(elnode(1:i34(i),i)))/dble(i34(i))
           yel(i)=sum(yy(elnode(1:i34(i),i)))/dble(i34(i))
         enddo

         max_dist=distance**2

         do i=1,nsink
           ie=ele_sink(i)
           do j=1,nsource
             ie_nbr=ele_source(j)
             dist=(xel(ie)-xel(ie_nbr))**2+(yel(ie)-yel(ie_nbr))**2
             if (dist<max_dist) then
               n_sink_source(i)=n_sink_source(i)+1
               i_sink_source(n_sink_source(i),i)=j
             endif
           enddo
         enddo
!---------------------------------------------------------------
       else !neighboring element tier 
!---------------------------------------------------------------
         allocate(tier1(0:mnei*4,ne),ntier1(ne),stat=istat) !tier1(0) should be 0 at all time
         if(istat/=0) then
           write(*,*)'failed in alloc. tier1'
           stop
         endif
         allocate(tier_n(0:mnei*400,nsink),ntier_n(nsink),stat=istat) 
         if(istat/=0) then
           write(*,*)'failed in alloc. tier_n'
           stop
         endif

         ncount=0
!$omp parallel do private (ie,j,k,n1,k0,kk,tid)
         do ie=1,ne
           tid = omp_get_thread_num()
           ncount(tid)=ncount(tid)+1
           ntier1(ie)=0 !number of tier 1 elements
           do j=1,i34(ie)
             n1=elnode(j,ie)
             !find position of ie in the nodal ball
             do k=1,nne(n1)
               if (indel(k,n1).eq.ie) then
                 k0=k
                 exit
               endif
             enddo !k
             do k=1,nne(n1)-1
               !counter-clockwise: from the first non-ie element to the last non-ie element
               !thus, the possible overlapping of tier1 elements (from next node of ie) can only 
               !occur at the first or last element in the current tier1
               kk=k+k0; if (kk>nne(n1)) kk=kk-nne(n1)
               if (indel(kk,n1).ne.0 .and. indel(kk,n1).ne.tier1(1,ie) .and. indel(kk,n1).ne.tier1(ntier1(ie),ie)) then
                 ntier1(ie)=ntier1(ie)+1
                 if (ntier1(ie)<=mnei*4) then
                   tier1(ntier1(ie),ie)=indel(kk,n1)
                 else
                   write(*,*) 'tier 1 > mnei*4'
                   stop
                 endif
               endif
             enddo !k
           enddo
         enddo !ie
!$omp end parallel do

         print*, ncount

         !mark sources
         allocate(is_source_ele(ne))
         is_source_ele=0
         do i=1,nsource
           is_source_ele(ele_source(i))=i
         enddo

         tier_n=0; ntier_n=0
!$omp parallel do private (i,ie)
         do i=1,nsink
           ie=ele_sink(i)
           tier_n(0,i)=ie !self
           call mark_neighbors(i,ie,0,n_tier_dist)
         enddo
!$omp end parallel do

       !call cpu_time(end_time)
       end_time=omp_get_wtime();
       print*, 'calculating geometry takes ',(end_time-start_time)/60d0,' minutes'

!---------------------------------------------------------------
       endif !inbr
!---------------------------------------------------------------

       !call cpu_time(start_time)
       start_time=omp_get_wtime();
       print*, 'redistributing sinks ...'

!      redistribute vsink to neighboring vsources
!$omp parallel do private (i,j,k,isource)
       do i=1,nsink
         if (n_sink_source(i)>0) then
           do k=1,nt
             do j=1,n_sink_source(i)
               isource=i_sink_source(j,i) 
               if (vsource(isource,k)>abs(vsink(i,k))) then
                 vsource(isource,k)=vsource(isource,k)+vsink(i,k)
                 vsink(i,k)=0d0
                 exit !no need for finding next neighboring source
               else
                 vsink(i,k)=vsink(i,k)+vsource(isource,k)
                 vsource(isource,k)=0.0d0
               endif
             enddo ! j=1,nsource
           enddo !k=1,nt
         endif
       enddo !i=1,nsink
!$omp end parallel do

       !call cpu_time(end_time)
       end_time=omp_get_wtime();
       print*, 'redistributing sinks takes ',(end_time-start_time)/60d0,' minutes'

      
!!write source_sink.in file
!      open(13,file='source_sink.in.1') 
!      if (n_source.ne.0) then
!        write(13,*)nsource
!        do i=1,nsource
!           write(13,*)eleso_uni(i)
!        enddo
!      else
!        write(13,*)0
!      endif
!      write(13,*)
!      if (n_sink.ne.0) then
!        write(13,*)nsi
!        do i=1,nsi
!           write(13,*)elesi_uni(i)
!        enddo
!      else
!        write(13,*)0
!      endif
!
!
!!write msource.th file
!      allocate(salt(nsource),temp(nsource))
!      do i=1,nsource
!         salt(i)=0
!         temp(i)=-9999
!      enddo
!      open(15,file='msource.th.1')
!      do i=1,ntime
!      write(15,'(10000F15.3)')tstep(i),(temp(k),k=1,nsource),(salt(k),k=1,nsource)
!     
!      enddo
! 
!     write vsource.th file
      open(16,file='vsource.th.1',status='replace')
      do i=1,nt
        write(16,'(10000F15.3)')time_stamp(i),(vsource(k,i),k=1,nsource)
      end do

!     write vsink.th
      open(17,file='vsink.th.1')
      do i=1,nt
        write(17,'(10000F15.3)')time_stamp(i),(vsink(k,i),k=1,nsink)
      end do



     contains
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

           if (is_source_ele(ie000)>0) then
             n_sink_source(isink)=n_sink_source(isink)+1
             i_sink_source(n_sink_source(isink),isink)=is_source_ele(ie000)
           endif
         endif

         call mark_neighbors(isink,ie000,i_depth+1,n_tier_dist)

       enddo loop1
     end subroutine mark_neighbors  

end
