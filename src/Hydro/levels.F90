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
! SCHISM MISCELLANEOUS SUBROUTINES
! subroutine zcoor
! subroutine levels1
! subroutine levels0
! function lindex
!===============================================================================
!===============================================================================
      subroutine zcoor(itag,inode,kbpl,ztmp)
!-------------------------------------------------------------------------------
!     Calculate z-coord. at a _wet_ node
!     Search for 'ivcor' for other changes
!-------------------------------------------------------------------------------
!#ifdef USE_MPIMODULE
!      use mpi
!#endif
      use schism_glbl, only: rkind,errmsg,ivcor,eta2,dp,kbp,nvrt,kz,h0,h_s, &
     &h_c,theta_b,theta_f,s_con1,sigma,ztot,cs,sigma_lcl,iplg
      use schism_msgp, only: parallel_abort
      implicit none
!#ifndef USE_MPIMODULE
      include 'mpif.h'
!#endif
      integer, intent(in) :: itag,inode !tag to indicate where this routine is called from
!      real(rkind), intent(in) :: dep,etal
      integer, intent(out) :: kbpl
      real(rkind), intent(out) :: ztmp(nvrt)

!     Local
      integer :: k,kin,m
      real(rkind) :: hmod2,z0,z_1,sp,tmp,z_pws(nvrt),z_sigma(nvrt)

      !Make sure it's wet
      if(dp(inode)+eta2(inode)<=h0) then
        write(errmsg,*)'ZCOOR: dry location:',dp(inode),eta2(inode),itag
        call parallel_abort(errmsg)
      endif

!     WARNING: explicitly specify bottom/surface to avoid underflow
      if(ivcor==2) then !SZ
        hmod2=min(dp(inode),h_s)
        ztmp(kz)=-hmod2 !to avoid underflow
        ztmp(nvrt)=eta2(inode)

        if(hmod2<=h_c) then
          do k=kz+1,nvrt-1
            kin=k-kz+1
            ztmp(k)=sigma(kin)*(hmod2+eta2(inode))+eta2(inode)
          enddo !k
        else if(eta2(inode)<=-h_c-(dp(inode)-h_c)*theta_f/s_con1) then
          write(errmsg,*)'ZCOOR: Pls choose a larger h_c:',eta2(inode),h_c,itag
          call parallel_abort(errmsg)
        else
          do k=kz+1,nvrt-1
            kin=k-kz+1
            ztmp(k)=eta2(inode)*(1+sigma(kin))+h_c*sigma(kin)+(hmod2-h_c)*cs(kin)
          enddo !k
        endif

        if(dp(inode)<=h_s) then
          kbpl=kz
        else !z levels
!         Find bottom index
          kbpl=0
          do k=1,kz-1
            if(-dp(inode)>=ztot(k).and.-dp(inode)<ztot(k+1)) then
              kbpl=k
              exit
            endif
          enddo !k
          !todo: assert
          if(kbpl==0) then
            write(errmsg,*)'ZCOOR: Cannot find a bottom level:',dp(inode),itag
            call parallel_abort(errmsg)
          endif
          ztmp(kbpl)=-dp(inode)
          do k=kbpl+1,kz-1
            ztmp(k)=ztot(k)
          enddo !k
        endif !dep<=h_s

      else if(ivcor==1) then !localized simga
!        if(eta<=-hsm(m_pws)) then
!          write(errmsg,*)'ZCOOR: elev<hsm:',eta,itag
!          call parallel_abort(errmsg)
!        endif

        kbpl=kbp(inode)
        do k=kbpl+1,nvrt-1
          ztmp(k)=(eta2(inode)+dp(inode))*sigma_lcl(k,inode)+eta2(inode)
        enddo !k

        ztmp(kbpl)=-dp(inode) !to avoid underflow
        ztmp(nvrt)=eta2(inode) !to avoid underflow
      else
        call parallel_abort('ZCOOR: unknown z-coor.')
      endif !ivcor

!#ifdef DEBUG
      do k=kbpl+1,nvrt
        !todo: assert
        if(ztmp(k)-ztmp(k-1)<=0._rkind) then
          write(12,*)'ZCOOR: Inverted z-level:',itag,ivcor,k,kbpl,iplg(inode),eta2(inode),dp(inode),ztmp(k),ztmp(k-1),sigma_lcl(kbpl:nvrt,inode)
          write(errmsg,*)'ZCOOR: Inverted z-level:',itag,ivcor,k,kbpl,iplg(inode),eta2(inode),dp(inode),ztmp(k),ztmp(k-1)
          call parallel_abort(errmsg)
        endif
      enddo !k
!#endif

      end subroutine zcoor
      
!===============================================================================

      ! Compute local index of a node (0 if not a local node)
      function lindex(node,ie)
      use schism_glbl
      use schism_msgp, only : parallel_abort
      implicit none
      integer :: lindex
      integer,intent(in) :: node,ie
      integer :: j

      lindex=0 !error flag
      do j=1,i34(ie)
        if(node==elnode(j,ie)) lindex=j
      enddo
!     if(lindex.eq.0) then
!       write(errmsg,*)'LINDEX: ',node,' is not in element ',ie
!       call parallel_abort(errmsg)
!     endif

      end function lindex
!===============================================================================


      subroutine levels1(iths,it)
!-------------------------------------------------------------------------------
! Routine to update level indices and wetting and drying.
! Used when resolution is fine enough.
! ONLY WORKS WITH PURE TRI's
!-------------------------------------------------------------------------------
!#ifdef USE_MPIMODULE
!      use mpi
!#endif
      use schism_glbl
      use schism_msgp
      implicit none
!#ifndef USE_MPIMODULE
      include 'mpif.h'
!#endif
      integer, intent(in) :: iths,it

!     Function
      integer :: lindex

!     Local
      integer :: i,j,k,l,m,nd,istop,itr,nsdf,nsdf_gb,isd,isd2,ie,ie2, &
                 &n1,n2,n3,n4,nodeA,inun,id,id1,l0,istat,iwet,icount,jj,kin
    
      real(rkind) :: cwtmp,tmp,flux_t,etm,dot11,dot12,dot21,dot22,stmp,ttmp

      integer :: idry2(npa),idry_s2(nsa),idry_e2(nea),isdf(nsa),inew(nsa), &
                 &icolor(npa),icolor2(nsa)
      real(rkind) :: out2(12+nvrt),sutmp(nvrt),svtmp(nvrt),swild2(2,nvrt)

      real(rkind),allocatable :: swild(:,:,:)
      logical :: srwt_xchng(1),prwt_xchng(1),ltmp
      logical :: srwt_xchng_gb(1),prwt_xchng_gb(1)
      logical :: cwtime
!-------------------------------------------------------------------------------

!     Flag for comm timing
      cwtime=it/=iths

!...  An element is wet if and only if depths at all nodes >h0 
!...  A node is wet if and only if at least one surrounding element is wet
!...  A side is wet if and only if at least one surrounding element is wet
!     Initialize element flags for first step

!$OMP parallel default(shared) private(i,j,nd)

      if(it==iths) then
!$OMP   workshare
        idry_e=0
!$OMP   end workshare

!$OMP   do
        do i=1,nea
          do j=1,i34(i)
            nd=elnode(j,i)
            if(eta2(nd)+dp(nd)<=h0) then
              idry_e(i)=1
              exit
            endif
          enddo !j
        enddo !i
!$OMP   end do
      endif !it

!      if(it/=iths) idry_e0=idry_e !save only for upwindtrack()

!...  Wetting/drying algorithm
!$OMP workshare
      idry_e2=idry_e !starting from step n's indices
!$OMP end workshare
!$OMP end parallel 

      if(it/=iths) then

!       Make dry first (to speed up iteration)
!        do i=1,np
!          if(dp(i)+eta2(i)<=h0) idry_e2(indel(1:nne(i),i))=1
!        enddo !i

!$OMP parallel do default(shared) private(i,j,nd)
        do i=1,ne
          do j=1,i34(i)
            nd=elnode(j,i)
            if(eta2(nd)+dp(nd)<=h0) then
              idry_e2(i)=1
              exit
            endif
          enddo !j
        enddo !i
!$OMP end parallel do

        call exchange_e2di(idry_e2)

!Debug
!        write(12,*)'it=',it  
!        if(it==321) then
!          fdb='tmp_0000'
!          lfdb=len_trim(fdb)
!          write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
!          open(10,file=out_dir(1:len_out_dir)//fdb,status='replace')
!          write(10,*)np
!          do i=1,np
!            write(10,*)iplg(i),real(eta2(i))
!          enddo !i
!          write(10,*)ns
!          do i=1,ns
!            write(10,*)i,iplg(isidenode(1:2,i)),real(su2(nvrt,i)),real(sv2(nvrt,i))
!          enddo !i
!          close(10)
!        endif

        !istop: 1- ready for final extrap. stage; 2- ready for final
        !checks and exit loop15
        istop=0 
        itr=0
        loop15: do
          itr=itr+1
          if(itr>100) call parallel_abort('LEVELS1: Too many iterations in wet/dry')
!'

!$OMP parallel default(shared) private(i,j,ie,id,m,isd)

!         Interface (shoreline) sides
!$OMP     workshare
!          icolor=0 !nodes on the interface sides (not needed)
          icolor2=0 !interface sides
!$OMP     end workshare

!$OMP     do
          do i=1,ns
            if(isdel(2,i)/=0) then; if(idry_e2(isdel(1,i))+idry_e2(isdel(2,i))==1) then
              icolor2(i)=1
            endif; endif
          enddo !i
!$OMP     end do

!!$OMP     do
!          loopinun: do i=1,np
!            do j=1,nne(i)
!              ie=indel(j,i)
!              id=iself(j,i)
!              do m=1,2 !2 neighboring sides
!                isd=elside(nxq(m+i34(ie)-3,id,i34(ie)),ie)
!                if(icolor2(isd)==1) then
!                  icolor(i)=1
!                  cycle loopinun
!                endif
!              enddo !m
!            enddo !j
!          end do loopinun !i
!!$OMP     end do
!$OMP end parallel

!          call exchange_p2di(icolor)
          call exchange_s2di(icolor2)
          
!         Aug. shoreline sides (must be internal sides)
          nsdf=0
          do i=1,nsa
            if(icolor2(i)==1) then
              nsdf=nsdf+1
              isdf(nsdf)=i
            endif
          enddo !i

          call mpi_allreduce(nsdf,nsdf_gb,1,itype,MPI_SUM,comm,ierr)
          if(nsdf_gb==0) exit loop15 !all wet

!         Final extrapolation
          srwt_xchng(1)=istop==1
          call mpi_allreduce(srwt_xchng,srwt_xchng_gb,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
          if(srwt_xchng_gb(1)) then !all ranks ready
            if(myrank==0) write(16,*)'doing final extrapolation in levels1...'
!'
            icolor=0 !frontier nodes for extrapolation
            inew=0 !for initializing and counting su2 sv2
            do i=1,nsdf !aug.
              isd=isdf(i)
!              if(isdel(1,isd)<0.or.isdel(2,isd)<0) cycle
              if(isdel(1,isd)==0.or.isdel(2,isd)==0) then
                write(errmsg,*)'LEVELS1: bnd side (2):',isdel(:,isd),iplg(isidenode(1:2,isd))
                call parallel_abort(errmsg)
              endif
!              if(idry_e2(isdel(1,isd))+idry_e2(isdel(2,isd))/=1) cycle

              !Try to find a dry elem (to take care of some odd cases where
              !nodeA is interface btw sub-domains)
!              if(idry_e2(isdel(1,isd))==1) then
!                ie=isdel(1,isd)
!              else 
!                ie=isdel(2,isd)
!              endif
              ie=0
              do m=1,2
                if(isdel(m,isd)>0) then; if(idry_e2(isdel(m,isd))==1) then
                  ie=isdel(m,isd); exit
                endif; endif
              enddo !m
              if(ie==0) cycle

              n1=isidenode(1,isd)
              n2=isidenode(2,isd)
              nodeA=elnode(1,ie)+elnode(2,ie)+elnode(3,ie)-n1-n2

              if(icolor(nodeA)==1) cycle !this node is done

              icolor(nodeA)=1 !this node will be done
              if(nodeA>np) cycle
!             nodeA is resident

              inun=0 !inundation flag
              do j=1,nne(nodeA)
                ie2=indel(j,nodeA)
                id=iself(j,nodeA)
                isd2=elside(id,ie2)
                if(icolor2(isd2)==1) then
!                  if(ics==1) then
                  tmp=su2(nvrt,isd2)*snx(isd2)+sv2(nvrt,isd2)*sny(isd2)
!                  else !ics=2
!                    tmp=su2(nvrt,isd2)
!                  endif !ics
                  flux_t=-tmp*ssign(id,ie2) !inward normal
                  if(flux_t>0._rkind) then
                    n1=isidenode(1,isd2)
                    n2=isidenode(2,isd2)
!                    avh=(eta2(n1)+dp(n1)+eta2(n2)+dp(n2))/2
!                    vol=flux_t*dt*avh*distj(isd2) !inflow volume in one step
!                    avh3=(eta2(n1)+dp(n1)+eta2(n2)+dp(n2))/3 !assume total depth at nodeA=0
!                    volmin=avh3*area(ie2)
                    etm=max(eta2(n1),eta2(n2))
                    if(etm+dp(nodeA)>h0) then
                      inun=1
                      exit
                    endif
                  endif !flux_t>0
                endif !icolor2(isd2)==1
              enddo !j

              if(inun==1) then
                eta2(nodeA)=max(eta2(nodeA),-dp(nodeA)+2._rkind*h0)
                do j=1,nne(nodeA)
                  ie2=indel(j,nodeA)
                  id=iself(j,nodeA)
                  isd2=elside(id,ie2)
                  if(icolor2(isd2)==1) then
                    do l=1,3
                      nd=elnode(l,ie2)
                      if(eta2(nd)+dp(nd)<=h0) then 
                        write(errmsg,*)'LEVELS1: Failed to wet element:',ielg(ie2),iplg(nodeA)
                        call parallel_abort(errmsg)
                      endif
                    enddo !l=1,3
                    idry_e2(ie2)=0
                    do l=1,2 !sides sharing nodeA
                      id1=elside(nx(id,l),ie2)
!                      if(ics==1) then
                      swild2(1,1:nvrt)=su2(1:nvrt,isd2)
                      swild2(2,1:nvrt)=sv2(1:nvrt,isd2)
!                      else !ics=2
!                        !Assuming plane rotation
!                        dot11=dot_product(sframe(1:3,1,isd2),sframe(1:3,1,id1))
!                        dot21=dot_product(sframe(1:3,2,isd2),sframe(1:3,1,id1))
!                        swild2(1,1:nvrt)=su2(1:nvrt,isd2)*dot11+sv2(1:nvrt,isd2)*dot21
!                        dot12=dot_product(sframe(1:3,1,isd2),sframe(1:3,2,id1))
!                        dot22=dot_product(sframe(1:3,2,isd2),sframe(1:3,2,id1))
!                        swild2(2,1:nvrt)=su2(1:nvrt,isd2)*dot12+sv2(1:nvrt,isd2)*dot22
!                      endif !ics
                      if(inew(id1)==0) then
                        su2(1:nvrt,id1)=swild2(1,1:nvrt)
                        sv2(1:nvrt,id1)=swild2(2,1:nvrt)
                        inew(id1)=1
                      else
                        su2(1:nvrt,id1)=su2(1:nvrt,id1)+swild2(1,1:nvrt)
                        sv2(1:nvrt,id1)=sv2(1:nvrt,id1)+swild2(2,1:nvrt)
                        inew(id1)=inew(id1)+1
                      endif
                    enddo !l=1,2
                  endif !icolor2(isd2)==1
                enddo !j=1,nne(nodeA)
              endif !inun==1
            enddo !i=1,nsdf

            call exchange_e2di(idry_e2)
            call exchange_p2d(eta2)

!            srwt_xchng(1)=.false. !flag for wetting occurring
            ltmp=.false. !flag for wetting occurring
!$OMP parallel do default(shared) private(i) reduction(.or.: ltmp)
            do i=1,ns
              ltmp=ltmp.or.inew(i)/=0
              if(inew(i)/=0) then
!                srwt_xchng(1)=.true.
                su2(1:nvrt,i)=su2(1:nvrt,i)/dble(inew(i))
                sv2(1:nvrt,i)=sv2(1:nvrt,i)/dble(inew(i))
              endif
            enddo !i
!$OMP end parallel do 
            srwt_xchng(1)=ltmp

            istop=2
!            go to 991
          endif !srwt_xchng_gb; final extrapolation

          if(istop/=2) then
!=========
            istop=1 !stop iteration and go to extrapolation stage; initialize first
            do i=1,nsdf !aug.
              isd=isdf(i)
              do j=1,2
                nd=isidenode(j,isd)
                if(eta2(nd)+dp(nd)<=h0) then
!Debug
!                write(12,*)'Make dry:',itr,iplg(nd)

                  istop=0
                  do l=1,nne(nd)
                    ie=indel(l,nd)
                    if(ie>0) idry_e2(ie)=1
                  enddo !l
                endif
              enddo !j=1,2 nodes
            enddo !i=1,nsdf
            call exchange_e2di(idry_e2)

!           Wetting
            inew=0 !for initializing and counting su2 sv2
            srwt_xchng(1)=.false. !flag for wetting occurring
            do i=1,nsdf !aug. domain for updating vel. at interfacial sides (between 2 sub-domains)
              isd=isdf(i) !must be internal side
              if(isdel(1,isd)<0.or.isdel(2,isd)<0) cycle !neither element can have interfacial sides
              if(isdel(1,isd)==0.or.isdel(2,isd)==0) then
                write(errmsg,*)'LEVELS1: bnd side:',isdel(:,isd),iplg(isidenode(1:2,isd))
                call parallel_abort(errmsg)
              endif
              if(idry_e2(isdel(1,isd))+idry_e2(isdel(2,isd))/=1) cycle
!             2 end nodes have total depths > h0

              if(idry_e2(isdel(1,isd))==1) then
                ie=isdel(1,isd) !>0
              else
                ie=isdel(2,isd) !>0
              endif
              n1=isidenode(1,isd)
              n2=isidenode(2,isd)
              nodeA=elnode(1,ie)+elnode(2,ie)+elnode(3,ie)-n1-n2   ! eli: is the 2,ie one right?
              l0=lindex(nodeA,ie)
!            if(l0==0.or.icolor(nodeA)==1.or.nodeA==n1.or.nodeA==n2) then
              if(l0==0.or.nodeA==n1.or.nodeA==n2) then
                write(errmsg,*)'Frontier node outside, or on the interface:', &
       &l0,iplg(nodeA),iplg(n1),iplg(n2),itr,it,iths !icolor(nodeA)
!'
                write(12,*)'LEVELS1: fatal error message'
                do l=1,ns
                  if(icolor2(l)==1) then
                    write(12,*)l,iplg(isidenode(1:2,l))
                    write(12,*)l,ielg(isdel(1:2,l)),idry_e2(isdel(1:2,l)),idry_e(isdel(1:2,l))
                  endif
                enddo !l
                do l=1,nea
                  write(12,*)l,idry_e2(l),idry_e(l)
                enddo !l
                call parallel_abort(errmsg)
              endif !end fatal

              if(eta2(nodeA)+dp(nodeA)>h0) then !all 3 nodes have depths > h0
!               Check
                do j=1,3
                  nd=elnode(j,ie)
                  if(eta2(nd)+dp(nd)<=h0) then
                    write(errmsg,*)'Failed to wet element (13):',ielg(ie),iplg(nd),iplg(nodeA)
                    call parallel_abort(errmsg)
                  endif
                enddo !j

!Debug
!              write(12,*)'Make wet:',itr,iplg(nodeA),ielg(ie)

                srwt_xchng(1)=.true.
                istop=0
                idry_e2(ie)=0

                do j=1,2 !sides sharing nodeA
                  id1=elside(nx(l0,j),ie)
                  if(icolor2(id1)==0) then

!                  if(ics==1) then
                    swild2(1,1:nvrt)=su2(1:nvrt,isd)
                    swild2(2,1:nvrt)=sv2(1:nvrt,isd)
!                  else !ics=2
!                    !Assuming plane rotation
!                    dot11=dot_product(sframe(1:3,1,isd),sframe(1:3,1,id1))
!                    dot21=dot_product(sframe(1:3,2,isd),sframe(1:3,1,id1))
!                    swild2(1,1:nvrt)=su2(1:nvrt,isd)*dot11+sv2(1:nvrt,isd)*dot21
!                    dot12=dot_product(sframe(1:3,1,isd),sframe(1:3,2,id1))
!                    dot22=dot_product(sframe(1:3,2,isd),sframe(1:3,2,id1))
!                    swild2(2,1:nvrt)=su2(1:nvrt,isd)*dot12+sv2(1:nvrt,isd)*dot22
!                  endif !ics

                    if(inew(id1)==0) then
                      !vel. only accurate in resident domain
                      su2(1:nvrt,id1)=swild2(1,1:nvrt) !su2(1:nvrt,isd)
                      sv2(1:nvrt,id1)=swild2(2,1:nvrt) !sv2(1:nvrt,isd)
                      inew(id1)=1
                    else
                      su2(1:nvrt,id1)=su2(1:nvrt,id1)+swild2(1,1:nvrt)
                      sv2(1:nvrt,id1)=sv2(1:nvrt,id1)+swild2(2,1:nvrt)
                      inew(id1)=inew(id1)+1
                    endif
                  endif !icolor2(id)==0
                enddo !j=1,2
              endif !eta2(nodeA)+dp(nodeA)>h0
            enddo !i=1,nsdf; shoreline sides

!         Compute average vel. for rewetted sides
!$OMP parallel do default(shared) private(i)
            do i=1,ns
              if(inew(i)/=0) then
                su2(1:nvrt,i)=su2(1:nvrt,i)/real(inew(i),rkind)
                sv2(1:nvrt,i)=sv2(1:nvrt,i)/real(inew(i),rkind)
              endif !inew(i)/=0
            enddo !i=1,ns
!$OMP end parallel do

!991       continue
!=========
          endif !istop/=2

          call mpi_allreduce(srwt_xchng,srwt_xchng_gb,1,MPI_LOGICAL,MPI_LOR,comm,ierr)
          if(srwt_xchng_gb(1)) then
            call exchange_e2di(idry_e2)
            allocate(swild(2,nvrt,nsa),stat=istat)
            if(istat/=0) call parallel_abort('Levels1: fail to allocate (9)')
!'
            swild(1,:,:)=su2(:,:)
            swild(2,:,:)=sv2(:,:)
!#ifdef INCLUDE_TIMING
            if(cwtime) cwtmp=mpi_wtime()
!#endif
            call exchange_s3d_2(swild)
!#ifdef INCLUDE_TIMING
            if(cwtime) wtimer(10,2)=wtimer(10,2)+mpi_wtime()-cwtmp
!#endif
            su2(:,:)=swild(1,:,:)
            sv2(:,:)=swild(2,:,:)
            deallocate(swild)
          endif !srwt_xchng_gb

          ltmp=.false. !for vel. exchange
!$OMP parallel default(shared) private(i,j,iwet,ie,sutmp,svtmp,icount,m,jj,isd2)

!         Enforce wet/dry flag consistency between nodes and elements due to added wet elements
!$OMP     workshare
          idry2=1
!$OMP     end workshare
!          do i=1,nea
!            if(idry_e2(i)==0) idry2(elnode(1:3,i))=0
!          enddo !i

!$OMP     do
          do i=1,np
            do j=1,nne(i)
              if(idry_e2(indel(j,i))==0) then
                idry2(i)=0; exit
              endif
            enddo !j
          enddo !i
!$OMP     end do

!$OMP     master
          call exchange_p2di(idry2)
!$OMP     end master
!$OMP     barrier

!         Compute su2 sv2 for newly wetted sides (due to reasons other than the wetting above)
!$OMP     do
          do i=1,nea
            inew(i)=0 !use for temp. storage of new element wet/dry flags
            do j=1,3
              if(idry2(elnode(j,i))==1) inew(i)=1
            enddo !j
          enddo !i=1,nea
!$OMP     end do

!$OMP     do reduction(.or.: ltmp)
!          srwt_xchng(1)=.false. !for vel. exchange
          do i=1,ns
            if(.not.(idry_e2(isdel(1,i))==1.and.(isdel(2,i)==0.or.isdel(2,i)>0.and.idry_e2(max(1,isdel(2,i)))==1))) cycle
!           Dry side that may need new vel.

            iwet=0 !flag
            do j=1,2
              ie=isdel(j,i)
              if(ie>0.and.idry_e2(max(1,ie))==1.and.inew(max(1,ie))==0) iwet=1
            enddo !j

            if(iwet==1) then !vel. as average
              sutmp=0._rkind; svtmp=0._rkind; icount=0
              do m=1,2 !2 elements
                ie=isdel(m,i)
                if(ie<=0) cycle

                do jj=1,3 !3 sides
                  !Find wet side
                  isd2=elside(jj,ie)
                  if(isdel(1,isd2)>0.and.idry_e2(max(1,isdel(1,isd2)))==0.or. &
     &isdel(2,isd2)>0.and.idry_e2(max(1,isdel(2,isd2)))==0) then !at least one wet element
                    icount=icount+1

!                    swild2(1,1:nvrt)=su2(1:nvrt,isd2)
!                    swild2(2,1:nvrt)=sv2(1:nvrt,isd2)
                    sutmp(1:nvrt)=sutmp(1:nvrt)+su2(1:nvrt,isd2)
                    svtmp(1:nvrt)=svtmp(1:nvrt)+sv2(1:nvrt,isd2)
                  endif
                enddo !jj
              enddo !m=1,2; 2 elements

              ltmp=ltmp.or.icount/=0
              if(icount/=0) then
!                srwt_xchng(1)=.true.
                su2(1:nvrt,i)=sutmp(1:nvrt)/real(icount,rkind)
                sv2(1:nvrt,i)=svtmp(1:nvrt)/real(icount,rkind)
              endif
            endif !iwet
          enddo !i=1,ns
!$OMP     end do

!$OMP     workshare
          idry_e2(1:nea)=inew(1:nea)
!$OMP     end workshare
!$OMP end parallel
          srwt_xchng(1)=ltmp

          call mpi_allreduce(srwt_xchng,srwt_xchng_gb,1,MPI_LOGICAL,MPI_LOR,comm,ierr)
          if(srwt_xchng_gb(1)) then
            allocate(swild(2,nvrt,nsa),stat=istat)
            if(istat/=0) call parallel_abort('Levels1: fail to allocate (8)')
!'
            swild(1,:,:)=su2(:,:)
            swild(2,:,:)=sv2(:,:)
!#ifdef INCLUDE_TIMING
            if(cwtime) cwtmp=mpi_wtime()
!#endif
            call exchange_s3d_2(swild)
!#ifdef INCLUDE_TIMING
            if(cwtime) wtimer(10,2)=wtimer(10,2)+mpi_wtime()-cwtmp
!#endif
            su2(:,:)=swild(1,:,:)
            sv2(:,:)=swild(2,:,:)
            deallocate(swild)
          endif !srwt_xchng_gb

!         Sync
          call parallel_barrier

          if(istop==2) exit loop15

        end do loop15

        if(myrank==0) then
          write(16,*)'see fort.7 for # of iterations used in LEVELS1...'
          write(7,*)it,itr
        endif
      endif !it/=iths

!$OMP parallel default(shared) private(i,j,nd,n1,n2,n3,k,stmp,ttmp,icount)

!...  Isolated dry nodes (do nothing for isolated wet)
!      do i=1,np
!        if(dp(i)+eta2(i)<=h0) idry_e2(indel(1:nne(i),i))=1
!      enddo !i

!$OMP do
      do i=1,ne
        do j=1,i34(i)
          nd=elnode(j,i)
          if(eta2(nd)+dp(nd)<=h0) then
            idry_e2(i)=1
            exit
          endif
        enddo !j
      enddo !i
!$OMP end do

!$OMP master
      call exchange_e2di(idry_e2)
!$OMP end master
!$OMP barrier

!...  Wet/dry flags for nodes/sides
!$OMP workshare
      idry2=1; idry_s2=1
!$OMP end workshare
!      do i=1,nea
!        if(idry_e2(i)==0) then
!          idry2(elnode(1:3,i))=0
!          idry_s2(elside(1:3,i))=0
!        endif
!      enddo !i

!$OMP do
      do i=1,np
        do j=1,nne(i)
          if(idry_e2(indel(j,i))==0) then
            idry2(i)=0; exit
          endif
        enddo !j
      enddo !i
!$OMP end do

!$OMP do
      do i=1,ns
        do j=1,2
          if(isdel(j,i)>0) then; if(idry_e2(isdel(j,i))==0) then
            idry_s2(i)=0; exit
          endif; endif
        enddo !j
      enddo !i
!$OMP end do

!$OMP master
      call exchange_p2di(idry2)
      call exchange_s2di(idry_s2)
!$OMP end master
!$OMP barrier

!...  Reset vel. at dry sides
!$OMP do
      do i=1,nsa
        if(idry_s2(i)==1) then
          su2(1:nvrt,i)=0._rkind
          sv2(1:nvrt,i)=0._rkind
        endif
      enddo !i
!$OMP end do

!...  Limit elevation at dry nodes
!$OMP do
      do i=1,npa
        if(idry2(i)==1) then
          !eta2(i)=min(0.d0,-dp(i))
          eta2(i)=min(eta2(i),-dp(i))
        endif
      enddo !i
!$OMP end do

!...  z-coor. for nodes
!...  
!$OMP do
      do i=1,npa
        if(ivcor==2) then; if(eta2(i)<=h0-h_s) then
          write(errmsg,*)'Deep depth dry:',iplg(i)
          call parallel_abort(errmsg)
        endif; endif

        if(idry2(i)==1) then
          if(ivcor/=1) kbp(i)=0
        else !wet
          call zcoor(1,i,kbp(i),znl(:,i))
        endif !wet ot dry
      enddo !i=1,npa
!$OMP end do

!     Debug
!      fdb='dry_0000'
!      lfdb=len_trim(fdb)
!      write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
!      open(10,file='outputs/'//fdb,status='unknown')
!      rewind(10)
!      write(10,*)'Time step=',it
!      write(10,*)'Node'
!      do i=1,npa
!        write(10,*)i,iplg(i),dp(i),eta2(i)
!      enddo !i

!     Compute element bottom index
!$OMP workshare
      kbe=0
!$OMP end workshare

!$OMP do
      do i=1,nea
        if(idry_e2(i)/=0) cycle

!       Wet
        n1=elnode(1,i); n2=elnode(2,i); n3=elnode(3,i)
        if(idry2(n1)/=0.or.idry2(n2)/=0.or.idry2(n3)/=0) then
          write(errmsg,*)'level1: Element-node inconsistency (0):',ielg(i),idry_e(i), &
     &iplg(elnode(1:3,i)),idry2(elnode(1:3,i))
          call parallel_abort(errmsg)
        endif
        kbe(i)=min(kbp(n1),kbp(n2),kbp(n3))
        do k=kbe(i),nvrt
          ze(k,i)=(znl(max(k,kbp(n1)),n1)+znl(max(k,kbp(n2)),n2)+znl(max(k,kbp(n3)),n3))/3._rkind
          if(k>=kbe(i)+1) then; if(ze(k,i)-ze(k-1,i)<=0._rkind) then
            write(errmsg,*)'Weird element (1):',k,i,ze(k,i),ze(k-1,i)
            call parallel_abort(errmsg)
          endif; endif
        enddo !k
      enddo !i
!$OMP end do

!     Compute side bottom index. For wet side and its wet adjacent element,
!     kbs>=kbe
!$OMP do
      do i=1,nsa
        kbs(i)=0 !dry
        if(idry_s2(i)==0) then !wet side with 2 wet nodes
          n1=isidenode(1,i)
          n2=isidenode(2,i)
          if(idry2(n1)/=0.or.idry2(n2)/=0) then
            write(errmsg,*)'Side-node inconsistency:',it,islg(i),'node:',iplg(n1),iplg(n2), &
     &eta2(n1),eta2(n2),idry2(n1),idry2(n2),';element:', &
     &(isdel(j,i),ielg(isdel(j,i)),idry_e2(isdel(j,i)),j=1,2)
            call parallel_abort(errmsg)
          endif
          if(dps(i)+(eta2(n1)+eta2(n2))/2._rkind<=h0) then
            write(errmsg,*)'Weird side (0):',islg(i),iplg(n1),iplg(n2),eta2(n1),eta2(n2)
            call parallel_abort(errmsg)
          endif
          kbs(i)=min(kbp(n1),kbp(n2))
          do k=kbs(i),nvrt
            zs(k,i)=(znl(max(k,kbp(n1)),n1)+znl(max(k,kbp(n2)),n2))/2
            if(k>=kbs(i)+1) then; if(zs(k,i)-zs(k-1,i)<=0._rkind) then
              write(errmsg,*)'Weird side (1):',k,iplg(n1),iplg(n2),znl(max(k,kbp(n1)),n1), &
     &znl(max(k,kbp(n2)),n2),znl(max(k-1,kbp(n1)),n1),znl(max(k-1,kbp(n2)),n2)
              call parallel_abort(errmsg)
            endif; endif
          enddo !k
        endif !wet side
      enddo !i=1,nsa
!$OMP end do

!     Compute vel., S,T for re-wetted nodes (q2 and xl are fine)
      if(it/=iths) then
!$OMP   do
        do i=1,np 
          if(idry(i)==1.and.idry2(i)==0) then
            do k=1,nvrt
              uu2(k,i)=0._rkind
              vv2(k,i)=0._rkind
              ttmp=0._rkind
              stmp=0._rkind
              icount=0
              do j=1,nnp(i)
                nd=indnd(j,i) !must be inside the aug. domain
!               Wet nbrs not affected by this part and so each sub-domain should use same values
                if(idry(nd)==0) then !all indices extended
                  icount=icount+1
                  uu2(k,i)=uu2(k,i)+uu2(k,nd)
                  vv2(k,i)=vv2(k,i)+vv2(k,nd)
                  ttmp=ttmp+tr_nd(1,k,nd) !tnd(k,nd)
                  stmp=stmp+tr_nd(2,k,nd) !snd(k,nd)
                endif
              enddo !j
              if(icount==0) then
                !Use last wet value
              else
                uu2(k,i)=uu2(k,i)/real(icount,rkind)
                vv2(k,i)=vv2(k,i)/real(icount,rkind)
                tr_nd(1,k,i)=ttmp/real(icount,rkind)
                tr_nd(2,k,i)=stmp/real(icount,rkind)
              endif
            enddo !k=1,nvrt
          endif !rewetted
        enddo !i=1,np
!$OMP   end do
      endif !it/=iths

!$OMP end parallel

!     Check wet/dry in ghost zone
      prwt_xchng(1)=.false.
      if(it/=iths) then
        do i=np+1,npa !check ghosts wet/dry
          if(idry(i)==1.and.idry2(i)==0) then
            prwt_xchng(1)=.true. !ghost rewetted; need exchange
            exit
          endif

!          if(idry(i)==1.and.idry2(i)==0) then
!            if(.not.prwt_xchng(1).and.i>np) prwt_xchng(1)=.true. !ghost
!            rewetted; need exchange
!            if(i>np) cycle !do rest for residents
        enddo !i
      endif !it/=iths

      if(nproc>1) then
!#ifdef INCLUDE_TIMING
        if(cwtime) cwtmp=mpi_wtime()
!#endif
!       See if the node exchange is needed
        call mpi_allreduce(prwt_xchng,prwt_xchng_gb,1,MPI_LOGICAL,MPI_LOR,comm,ierr)
        if(ierr/=MPI_SUCCESS) call parallel_abort('levels1: allreduce prwt_xchng_gb',ierr)
!'
!#ifdef INCLUDE_TIMING
        if(cwtime) wtimer(10,2)=wtimer(10,2)+mpi_wtime()-cwtmp
!#endif

!       update ghost nodes
        if(prwt_xchng_gb(1)) then
          allocate(swild(4,nvrt,nsa),stat=istat)
          if(istat/=0) call parallel_abort('Levels0: fail to allocate swild')
!'
          swild(1,:,1:npa)=uu2(:,:)
          swild(2,:,1:npa)=vv2(:,:)
          swild(3,:,1:npa)=tr_nd(1,:,:) !tnd(:,:)
          swild(4,:,1:npa)=tr_nd(2,:,:) !snd(:,:)
!#ifdef INCLUDE_TIMING
          if(cwtime) cwtmp=mpi_wtime()
!#endif
          call exchange_p3d_4(swild)
!#ifdef INCLUDE_TIMING
          if(cwtime) wtimer(10,2)=wtimer(10,2)+mpi_wtime()-cwtmp
!#endif
          uu2(:,:)=swild(1,:,1:npa)
          vv2(:,:)=swild(2,:,1:npa)
          tr_nd(1,:,:)=swild(3,:,1:npa)
          tr_nd(2,:,:)=swild(4,:,1:npa)
          deallocate(swild)
        endif !prwt_xchng_gb
      endif !nproc>1

!      close(10)

!...  Update wet/dry flags
      idry=idry2
      idry_s=idry_s2
      idry_e=idry_e2

      end subroutine levels1

!===============================================================================
!===============================================================================

      subroutine levels0(iths,it)
!-------------------------------------------------------------------------------
! Routine to update level indices and wetting and drying.
! Use levels1() for better inundation if resolution is fine enough.
!-------------------------------------------------------------------------------
!#ifdef USE_MPIMODULE
!      use mpi
!#endif
      use schism_glbl
      use schism_msgp
      implicit none
!#ifndef USE_MPIMODULE
      include 'mpif.h'
!#endif
      integer, intent(in) :: iths,it

!     Local
      integer :: i,j,k,kin,ie,ifl,n1,n2,n3,n4,icount,nd,isd,jj,istat
      real(rkind) :: cwtmp,utmp,vtmp,stmp,ttmp,dot11,dot12,dot21,dot22

      integer :: idry2(npa),idry_s2(nsa),idry_e2(nea)
      real(rkind) :: swild2(2)
      real(rkind),allocatable :: swild(:,:,:)
      logical :: srwt_xchng(1),prwt_xchng(1)
      logical :: srwt_xchng_gb(1),prwt_xchng_gb(1)
      logical :: cwtime
!-------------------------------------------------------------------------------

! Flag for comm timing
      cwtime=it/=iths
!$OMP parallel default(shared) private(i,j,ie,n1,n2,n3,n4,k,utmp,vtmp,ttmp,stmp,icount,nd,jj,isd)

!...  z-coor. for nodes
!...  
!$OMP do
      do i=1,npa
        if(dp(i)+eta2(i)<=h0) then !dry
          idry2(i)=1 
          if(ivcor==2) then; if(dp(i)>=h_s) then
            write(errmsg,*)'Deep depth dry:',iplg(i)
            call parallel_abort(errmsg)
          endif; endif
          if(ivcor/=1) kbp(i)=0
        else !wet
          idry2(i)=0
          call zcoor(0,i,kbp(i),znl(:,i))
        endif !wet ot dry
      enddo !i=1,npa
!$OMP end do

!     Debug
!      fdb='dry_0000'
!      lfdb=len_trim(fdb)
!      write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
!      open(10,file='outputs/'//fdb,status='unknown')
!      rewind(10)
!      write(10,*)'Time step=',it
!      write(10,*)'Node'
!      do i=1,npa
!        write(10,*)i,iplg(i),dp(i),eta2(i),idry2(i)
!      enddo !i

!...  Set wet/dry flags for element; element is "dry" if one of nodes is dry; conversely, 
!...  an element is wet if all nodes are wet (and all sides are wet as well)
!...  Weed out fake wet nodes; a node is wet if and only if at least one surrounding element is wet
!...
!      if(it/=iths) idry_e0=idry_e !save only for upwindtrack()

!$OMP do
      do i=1,nea
        idry_e2(i)=maxval(idry2(elnode(1:i34(i),i)))
      enddo !i
!$OMP end do

!      write(10,*)'Element'
!      do i=1,nea
!        write(10,*)i,ielg(i),idry_e2(i)
!      enddo !i

!$OMP workshare
      idry2=1 !dry unless wet
!$OMP end workshare

!$OMP do
      do i=1,np
        do j=1,nne(i)
          ie=indel(j,i)
          if(idry_e2(ie)==0) then
            idry2(i)=0; exit
          endif
        enddo !j
      enddo !i
!$OMP end do

!$OMP master
!#ifdef INCLUDE_TIMING
      if(cwtime) cwtmp=mpi_wtime()
!#endif
      call exchange_p2di(idry2) !update ghost values
!#ifdef INCLUDE_TIMING
      if(cwtime) wtimer(10,2)=wtimer(10,2)+mpi_wtime()-cwtmp
!#endif
!$OMP end master
!$OMP barrier

!      write(10,*)'nodes'
!      do i=1,npa
!        write(10,*)i,iplg(i),idry2(i),np
!      enddo !i

!     Consistency check
!#ifdef DEBUG
!      do i=1,npa
!        if(idry2(i)==1) cycle
! 
!        if(eta2(i)+dp(i)<=h0) then
!          write(errmsg,*)'levels0: weird wet node:',iplg(i),eta2(i),dp(i),idry2(i)
!          call parallel_abort(errmsg)
!        endif
!
!        if(i>np) cycle !do rest for residents only
!        ifl=0
!        do j=1,nne(i)
!          ie=indel(j,i)
!          if(idry_e2(ie)==0) then
!            ifl=1; exit
!          endif 
!        enddo !j
!        if(ifl==0) then
!          write(errmsg,*)'Node-element inconsistency:',iplg(i),idry2(i),(idry_e2(indel(j,i)),j=1,nne(i))
!          call parallel_abort(errmsg)
!        endif
!      enddo !i=1,npa
!#endif

!     Compute element bottom index
!$OMP workshare
      kbe=0
!$OMP end workshare

!$OMP do
      do i=1,nea
        if(idry_e2(i)/=0) cycle

!       Wet
        n1=elnode(1,i); n2=elnode(2,i); n3=elnode(3,i)
        if(maxval(idry2(elnode(1:i34(i),i)))/=0) then
          write(errmsg,*)'level0: Element-node inconsistency (0):',ielg(i),idry_e2(i), &
     &iplg(elnode(1:i34(i),i)),idry2(elnode(1:i34(i),i)),idry(elnode(1:i34(i),i))
          call parallel_abort(errmsg)
        endif
        kbe(i)=minval(kbp(elnode(1:i34(i),i)))
        do k=kbe(i),nvrt
          ze(k,i)=znl(max(k,kbp(n1)),n1)+znl(max(k,kbp(n2)),n2)+znl(max(k,kbp(n3)),n3)
          if(i34(i)==4) then
            n4=elnode(4,i)
            ze(k,i)=ze(k,i)+znl(max(k,kbp(n4)),n4)
          endif
          ze(k,i)=ze(k,i)/real(i34(i),rkind)
          if(k>=kbe(i)+1) then; if(ze(k,i)-ze(k-1,i)<=0._rkind) then
            write(errmsg,*)'Weird element (2):',k,i,ze(k,i),ze(k-1,i)
            call parallel_abort(errmsg)
          endif; endif
        enddo !k
      enddo !i
!$OMP end do

!     Compute vel., S,T for re-wetted nodes (q2 and xl are fine)
      if(it/=iths) then
!$OMP   do
        do i=1,np
          if(idry(i)==1.and.idry2(i)==0) then !rewetted
            do k=1,nvrt
              !uu2(k,i)=0
              !vv2(k,i)=0
              utmp=0._rkind
              vtmp=0._rkind
              ttmp=0._rkind
              stmp=0._rkind
              icount=0
              do j=1,nnp(i)
                nd=indnd(j,i) !must be inside the aug. domain
!               Wet nbrs not affected by this part and so each sub-domain should use same values
                if(idry(nd)==0) then !all indices extended
                  icount=icount+1
                  !Assume small element size in wet/dry zone so pframes are close to each other
                  utmp=utmp+uu2(k,nd)
                  vtmp=vtmp+vv2(k,nd)
                  ttmp=ttmp+tr_nd(1,k,nd) !tnd(k,nd)
                  stmp=stmp+tr_nd(2,k,nd) !snd(k,nd)
                endif
              enddo !j
              if(icount==0) then
!                if(ifort12(7)==0) then
!                  ifort12(7)=1
!                  write(12,*)'Isolated rewetted node:',iplg(i)
!                endif
!                tr_nd(1,k,i)=(k,i)
!                tr_nd(2,k,i)=(k,i)
              else
                uu2(k,i)=utmp/real(icount,rkind)
                vv2(k,i)=vtmp/real(icount,rkind)
                tr_nd(1,k,i)=ttmp/real(icount,rkind)
                tr_nd(2,k,i)=stmp/real(icount,rkind)
              endif
            enddo !k=1,nvrt
          endif !rewetted
        enddo !i=1,npa
!$OMP   end do
      endif !it/=iths

!...  z-coor. for sides
!...  A side is wet if and only if at least one of its elements is wet
!$OMP workshare
      idry_s2=1 !reinitialize to wipe out previous temp. storage
!$OMP end workshare

!$OMP do
      do i=1,ns
        do j=1,2 !elements
          ie=isdel(j,i)
          if(ie/=0.and.idry_e2(max(1,ie))==0) idry_s2(i)=0
        enddo !j
      enddo !i
!$OMP end do

!$OMP master
!#ifdef INCLUDE_TIMING
      if(cwtime) cwtmp=mpi_wtime()
!#endif
      call exchange_s2di(idry_s2) !update ghost values
!#ifdef INCLUDE_TIMING
      if(cwtime) wtimer(10,2)=wtimer(10,2)+mpi_wtime()-cwtmp
!#endif
!$OMP end master
!$OMP barrier

!      write(10,*)'Side'
!      do i=1,nsa
!        write(10,*)i,islg(i),idry_s2(i),ns
!      enddo !i

!     Consistency checks
!#ifdef DEBUG
!      do i=1,nea
!        if(idry_e2(i)/=0) cycle
!!       Wet
!        do j=1,3
!          isd=elside(j,i)
!          if(idry_s2(isd)/=0) then
!            write(errmsg,*)'Element-side inconsistency:',ielg(i),islg(isd),idry_s2(isd)
!            call parallel_abort(errmsg)
!          endif
!        enddo !j
!      enddo !i
!
!      do i=1,ns
!        if(idry_s2(i)==1) cycle
!
!        ifl=0
!        do j=1,2
!          ie=isdel(j,i)
!          if(ie/=0.and.idry_e2(max(1,ie))==0) then
!            ifl=1; exit
!          endif
!        enddo !j
!        if(ifl==0) then
!          write(errmsg,*)'Side-element inconsistency:',islg(i),idry_s2(i), &
!                         (isdel(j,i),idry_e2(isdel(j,i)),j=1,2)
!          call parallel_abort(errmsg)
!        endif
!      enddo !i
!#endif

!     Compute side bottom index
!$OMP do
      do i=1,nsa
        n1=isidenode(1,i)
        n2=isidenode(2,i)
        kbs(i)=0 !dry
        if(idry_s2(i)==0) then !wet side with 2 wet nodes
          if(idry2(n1)/=0.or.idry2(n2)/=0) then
            write(errmsg,*)'Side-node inconsistency (1):',it,islg(i),'node:',iplg(n1),iplg(n2), &
!'
             &eta2(n1),eta2(n2),idry2(n1),idry2(n2),';element:', &
             &(isdel(j,i),ielg(isdel(j,i)),idry_e2(isdel(j,i)),j=1,2)
            call parallel_abort(errmsg)
          endif
          if(dps(i)+(eta2(n1)+eta2(n2))/2._rkind<=h0) then
            write(errmsg,*)'Weird side (2):',islg(i),iplg(n1),iplg(n2),eta2(n1),eta2(n2)
            call parallel_abort(errmsg)
          endif
          kbs(i)=min(kbp(n1),kbp(n2))
          do k=kbs(i),nvrt
            zs(k,i)=(znl(max(k,kbp(n1)),n1)+znl(max(k,kbp(n2)),n2))/2._rkind
            if(k>=kbs(i)+1) then; if(zs(k,i)-zs(k-1,i)<=0._rkind) then
              write(errmsg,*)'Weird side (3):',k,iplg(n1),iplg(n2),znl(max(k,kbp(n1)),n1), &
     &znl(max(k,kbp(n2)),n2),znl(max(k-1,kbp(n1)),n1),znl(max(k-1,kbp(n2)),n2)
              call parallel_abort(errmsg)
            endif; endif
          enddo !k
        endif !wet side
      enddo !i=1,nsa
!$OMP end do

!     Compute vel., S,T for re-wetted sides 
      if(it/=iths) then
!$OMP   do
        do i=1,ns
          if(idry_s(i)==1.and.idry_s2(i)==0) then
            n1=isidenode(1,i)
            n2=isidenode(2,i)
            do k=1,nvrt
              utmp=0._rkind
              vtmp=0._rkind
              !ttmp=0
              !stmp=0
              icount=0
              do j=1,2
                ie=isdel(j,i)
                if(ie/=0) then
                  if(ie<0) call parallel_abort('levels0: ghost element')
                  do jj=1,i34(ie) !side; in the aug. domain
                    isd=elside(jj,ie)
                    if(idry_s(isd)==0) then
                      icount=icount+1

!                      if(ics==1) then
!                        swild2(1)=su2(k,isd)
!                        swild2(2)=sv2(k,isd)
!                      else !ics=2
!                        !Assuming plane rotation
!                        dot11=dot_product(sframe(1:3,1,isd),sframe(1:3,1,i))
!                        dot21=dot_product(sframe(1:3,2,isd),sframe(1:3,1,i))
!                        swild2(1)=su2(k,isd)*dot11+sv2(k,isd)*dot21
!                        dot12=dot_product(sframe(1:3,1,isd),sframe(1:3,2,i))
!                        dot22=dot_product(sframe(1:3,2,isd),sframe(1:3,2,i))
!                        swild2(2)=su2(k,isd)*dot12+sv2(k,isd)*dot22
!                      endif !ics

                      utmp=utmp+su2(k,isd)
                      vtmp=vtmp+sv2(k,isd)
                    endif
                  enddo !jj
                endif !ie/=0
              enddo !j
              if(icount==0) then
              else
                su2(k,i)=utmp/real(icount,rkind)
                sv2(k,i)=vtmp/real(icount,rkind)
              endif
            enddo !k
          endif !rewetted
        enddo !i=1,ns
!$OMP   end do
      endif !it/=iths

!$OMP end parallel

!     Check wet/dry in ghost zone
      prwt_xchng(1)=.false. !node
      srwt_xchng(1)=.false. !side

      if(it/=iths) then
        do i=np+1,npa !check ghosts wet/dry
          if(idry(i)==1.and.idry2(i)==0) then
            prwt_xchng(1)=.true. !ghost rewetted; need exchange
            exit
          endif

!          if(idry(i)==1.and.idry2(i)==0) then
!            if(.not.prwt_xchng(1).and.i>np) prwt_xchng(1)=.true. !ghost rewetted; need exchange
!            if(i>np) cycle !do rest for residents
        enddo !i
  
        do i=ns+1,nsa !ghost
          if(idry_s(i)==1.and.idry_s2(i)==0) then
            srwt_xchng(1)=.true.
            exit
          endif
 
!          if(idry_s(i)==1.and.idry_s2(i)==0) then
!            if(.not.srwt_xchng(1).and.i>ns) srwt_xchng(1)=.true. !rewetted ghost side; needs exchange
!            if(i>ns) cycle !do the rest only for residents
        enddo !i
      endif !it/

      if(nproc>1) then
!#ifdef INCLUDE_TIMING
        if(cwtime) cwtmp=mpi_wtime()
!#endif
!       See if the node/side exchange is needed
        call mpi_allreduce(prwt_xchng,prwt_xchng_gb,1,MPI_LOGICAL,MPI_LOR,comm,ierr)
        if(ierr/=MPI_SUCCESS) call parallel_abort('levels0: allreduce prwt_xchng_gb',ierr)
        call mpi_allreduce(srwt_xchng,srwt_xchng_gb,1,MPI_LOGICAL,MPI_LOR,comm,ierr)
        if(ierr/=MPI_SUCCESS) call parallel_abort('levels0: allreduce srwt_xchng_gb',ierr)
!'
!#ifdef INCLUDE_TIMING
        if(cwtime) wtimer(10,2)=wtimer(10,2)+mpi_wtime()-cwtmp
!#endif

!       Allocate temporary array
        if(prwt_xchng_gb(1).or.srwt_xchng_gb(1)) then
          allocate(swild(4,nvrt,nsa),stat=istat)
          if(istat/=0) call parallel_abort('Levels0: fail to allocate swild')
!'
        endif

!       update ghost nodes
        if(prwt_xchng_gb(1)) then
          swild(1,:,1:npa)=uu2(:,:)
          swild(2,:,1:npa)=vv2(:,:)
          swild(3,:,1:npa)=tr_nd(1,:,:) !tnd(:,:)
          swild(4,:,1:npa)=tr_nd(2,:,:) !snd(:,:)
!#ifdef INCLUDE_TIMING
          if(cwtime) cwtmp=mpi_wtime()
!#endif
          call exchange_p3d_4(swild)
!#ifdef INCLUDE_TIMING
          if(cwtime) wtimer(10,2)=wtimer(10,2)+mpi_wtime()-cwtmp
!#endif
          uu2(:,:)=swild(1,:,1:npa)
          vv2(:,:)=swild(2,:,1:npa)
          tr_nd(1,:,:)=swild(3,:,1:npa)
          tr_nd(2,:,:)=swild(4,:,1:npa)
        endif

!       update ghost sides
        if(srwt_xchng_gb(1)) then
          swild(1,:,:)=su2(:,:)
          swild(2,:,:)=sv2(:,:)
          swild(3,:,:)=0 !tsd(:,:) - not used
          swild(4,:,:)=0 !ssd(:,:)
!#ifdef INCLUDE_TIMING
          if(cwtime) cwtmp=mpi_wtime()
!#endif
          call exchange_s3d_4(swild)
!#ifdef INCLUDE_TIMING
          if(cwtime) wtimer(10,2)=wtimer(10,2)+mpi_wtime()-cwtmp
!#endif
          su2(:,:)=swild(1,:,:)
          sv2(:,:)=swild(2,:,:)
          !tsd(:,:)=swild(3,:,:)
          !ssd(:,:)=swild(4,:,:)
        endif

        if(prwt_xchng_gb(1).or.srwt_xchng_gb(1)) deallocate(swild)
      endif !nproc>1

!      close(10)

!     Update flags
      idry=idry2
      idry_s=idry_s2
      idry_e=idry_e2

      end subroutine levels0

!===============================================================================
!===============================================================================
