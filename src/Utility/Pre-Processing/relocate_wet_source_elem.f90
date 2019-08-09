!     Generate new source_sink.in to mv source from dry elem to nearest wet elem
!     
!     Input: 
!           (1) hgrid.gr3 (with bathy, 
!               bnd part is not necessary for this script, but needed for former one to build source_sink.in)
!           (2) source_sink.in
!           (3) screen
!     Output: source_sink.in.new

!     ifort -mcmodel=medium -O2 -o relocate_wet_source_elem relocate_wet_source_elem.f90

      implicit real*8(a-h,o-z)
      parameter(nbyte=4)
      real, allocatable :: xnd(:),ynd(:),dp(:)
      real, allocatable :: xctr(:),yctr(:)
      integer, allocatable ::nne(:),ndelem(:,:)
      integer, allocatable ::elnode(:,:),i34(:),elld(:),i34ld(:),el_s(:),el_sn(:),el_tar(:)
      integer, allocatable ::dryel(:),dryld(:),elld_new(:),elldnode(:,:)
      integer, allocatable :: icolor(:),idup(:)
      integer :: ne,nd
      integer :: el_c !tmp el id for each loop
      real :: dmin

      print*, 'Input depth for wet/dry determination in m:'
      read*, dmin 


      !read in hgrid.gr3
      !nodes: xnd(nd),ynd(nd),dp(nd)
      !elem: elnode(4,ne),i34(ne),xcrt(ne),yctr(ne)
      open(14,file='hgrid.gr3',status='old')
      read(14,*); read(14,*)ne,nd
      allocate(xnd(nd),ynd(nd),dp(nd),i34(ne),elnode(4,ne),xctr(ne),yctr(ne),icolor(ne),idup(ne),stat=istat)
      if(istat/=0) stop 'Failed to alloc: hgrid.gr3'
      elnode=0!init
      icolor=0
      idup=0
      do i=1,nd
        read(14,*)j,xnd(i),ynd(i),dp(i)
      enddo !i
      do i=1,ne
        read(14,*)j,k,elnode(1:k,i)
        i34(i)=k
        xctr(i)=sum(xnd(elnode(1:k,i)))/k
        yctr(i)=sum(ynd(elnode(1:k,i)))/k
      enddo !i
      close(14)

      bathmax=maxval(dp)
      if(bathmax<dmin) stop 'check bathymetry'


      !nodes: nne(nd),indel(maxval(nne),nd)
      allocate(nne(nd),stat=istat)
      if(istat/=0) stop 'Failed to alloc: nne'
      nne=0!init counting
      do i=1,ne
        do j=1,i34(i)
          nd=elnode(j,i)
          nne(nd)=nne(nd)+1
        enddo
      enddo
      nne_d=maxval(nne)

      allocate(ndelem(nne_d,nd),stat=istat)
      if(istat/=0) stop 'Failed to alloc: ndelem'
      ndelem=0;nne=0
      do i=1,ne
        do j=1,i34(i)
          nd=elnode(j,i)
          nne(nd)=nne(nd)+1
          ndelem(nne(nd),nd)=i
        enddo
      enddo


      !read in source_sink.in
      open(13,file='source_sink.in',status='old')
      read(13,*)nld
      allocate(elld(nld),dryld(nld),elld_new(nld),stat=istat)
      if(istat/=0) stop 'Failed to alloc: source_sink.in'
      elld_new=0!init 
      do i=1,nld
        read(13,*)elld(i)
      enddo!i
      close(13)


      !pick up dry elem from lding elem
      !icolor(ne),mark checked&&dry elem
      idry=0 !counting dry elem
      dryld=0 !init
      do i=1,nld
        el_c=elld(i)!elem id
!        i34_c=i34(el_c)!i34 
!        do j=1,i34_c
!          nd_c=elnode(j,el_c)!node id
!          if(dmin>dp(nd_c))then
!            idry=idry+1
!            dryld(i)=1 !dry flag
!            exit !j
!          endif!dmin
!        enddo!j

        !print*, 'el_c=',el_c
        call ifdry(el_c,ne,elnode,i34,nd,dp,dmin,isdry)
        !print*, 'isdry=',isdry
        dryld(i)=isdry
        if(isdry==1)then
          idry=idry+1
          icolor(el_c)=1!icolor(ne)
        endif!isdry
      enddo!i

      if(idry==0)then
        print*, 'All wet elem, take original source_sink.in'
        stop
      endif!idry


!      !set dry/wet prop table for each elem of ori domain
!      allocate(dryel(ne),stat=istat)
!      do i=1,ne !elem id
!        do j=1,i34(i)
!          nd_c=elnode(j,i)!node id
!          if(dmin>dp(nd_c))then
!            dryel(i)=1 !dry flag
!            exit !j
!          endif!dmin
!        enddo!j
!      enddo!i


      !renew elld to elld_new
      !el_s(ne):center elem to search 
      !el_tar(ne):target wet elem to check distance
      allocate(el_s(ne),el_sn(ne),el_tar(ne),stat=istat)
      if(istat/=0) stop 'Failed to alloc: el_s,el_sn,el_tar'
      do i=1,nld
        if(dryld(i)==1)then!dry
          !go through local search
          !init with local dry elem
          el_ori=elld(i)
          print*,'start searching round for elem',el_ori
          el_s(1)=elld(i)!elem id
          nel_s=1
          el_sn=0;nel_sn=0
          nel_tar=0
          do !unlimited rounds of surrounding search
            idup=0 !renew for each round, flag for checked elem
            do j=1,nel_s !go through all the center elem
              i34_c=i34(el_s(j))
              do k=1,i34_c !go throgh all the nodes of the center elem
                nd_c=elnode(k,el_s(j))
                nne_c=nne(nd_c)
                do p=1,nne_c !go through all the surrounding elem to this node
                  el_c=ndelem(p,nd_c)
                  !print*,'check elem, icolor, idup',el_c,icolor(el_c),idup(el_c)
                  !saving center elem list for next round
                  if(idup(el_c)==0)then!have not been checked this round
                    nel_sn=nel_sn+1
                    el_sn(nel_sn)=el_c       
                  endif!idup
                  if(icolor(el_c)==0.and.idup(el_c)==0)then!not dry, have not been checked this round
                    call ifdry(el_c,ne,elnode,i34,nd,dp,dmin,isdry)
                    !print*,'check elem_wet/dry',el_c,isdry
                    if(isdry==0)then!wet elem, go into el_tar
                      nel_tar=nel_tar+1
                      el_tar(nel_tar)=el_c   
                      idup(el_c)=1
                    else
                      icolor(el_c)=1
                      idup(el_c)=1
                    endif!isdry
                  endif!icolor&idup
                enddo!p,each surounding elem
              enddo!k,each node for certer elem
            enddo !j, center elem

            !check nel_tar for this round
            print*,'nel_tar=',nel_tar
            if(nel_tar>0)then
              exit !unlimited rounds
            else
              !renew el_s and nel_s for next round
              el_s=0;el_s(1:nel_sn)=el_sn(1:nel_sn)
              nel_s=nel_sn
            endif!nel_tar

            !extreme case, no wet elem
            if(minval(icolor)>0) stop 'check bathymetry (2)' 

          enddo !unlimited rounds of surrounding search

          !find el_pk, which has the smallest distance amoung el_tar
          distm2=0
          do j=1,nel_tar
            el_c=el_tar(j)
            distm_tmp=(xctr(el_c)-xctr(el_ori))**2+(yctr(el_c)-yctr(el_ori))**2
            if(distm2==0)then
              distm2=distm_tmp
              el_pk=el_c
            elseif(distm2>distm_tmp)then
              distm2=distm_tmp
              el_pk=el_c
            endif!distm2
          enddo!j=1:nel_tar
          elld_new(i)=el_pk
        else!keep wet elem
          elld_new(i)=elld(i)
        endif!dryld(i)

      enddo !i=1,nld


      !write out new source_sink.in.new
      open(12,file='source_sink.in.new',status='replace')
      write(12,*)nld,"#source"
      do i=1,nld
        write(12,*)elld_new(i)
      enddo!i

      contains
      subroutine ifdry(el_id,ne,elnode,i34,nd,dp,dmin,isdry)
        implicit real*8(a-h,o-z)
        integer, intent(in):: el_id 
        integer, intent(in):: ne,nd
        integer, allocatable :: elnode(:,:),i34(:),nd_ar(:)
        real,allocatable :: dp(:)
        integer :: isdry,i34_c,nd_c 
        real, intent(in) :: dmin

        allocate(dp(nd),elnode(4,ne),i34(ne),stat=istat)
        !print*,'el_c=',el_id
        isdry=0
        i34_c=i34(el_id)
!        do j=1,i34_c
!          nd_c=elnode(j,el_id)!node id
!          !print*,'nd_c=',nd_c
!          !print*,'dp=',dp(nd_c)
!          if(dmin>dp(nd_c))then
!            isdry=1 !dry flag
!            return !j
!          endif!dmin
!        enddo!j

        allocate(nd_ar(i34_c),stat=istat)
        nd_ar=elnode(:,el_id)
        if(dmin>minval(dp(nd_ar(:))))then
          isdry=1
        endif!dmin

      endsubroutine ifdry

      end 

