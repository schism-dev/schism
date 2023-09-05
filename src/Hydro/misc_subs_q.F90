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
! SCHISM MISCELLANEOUS SUBROUTINES used by QSim
! subroutine levels0
! subroutine compute_ll
! subroutine cross_product
! function signa2

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
#ifdef INCLUDE_TIMING
      if(cwtime) cwtmp=mpi_wtime()
#endif
      call exchange_p2di(idry2) !update ghost values
#ifdef INCLUDE_TIMING
      if(cwtime) wtimer(10,2)=wtimer(10,2)+mpi_wtime()-cwtmp
#endif
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
#ifdef INCLUDE_TIMING
      if(cwtime) cwtmp=mpi_wtime()
#endif
      call exchange_s2di(idry_s2) !update ghost values
#ifdef INCLUDE_TIMING
      if(cwtime) wtimer(10,2)=wtimer(10,2)+mpi_wtime()-cwtmp
#endif
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
#ifdef INCLUDE_TIMING
        if(cwtime) cwtmp=mpi_wtime()
#endif
!       See if the node/side exchange is needed
        call mpi_allreduce(prwt_xchng,prwt_xchng_gb,1,MPI_LOGICAL,MPI_LOR,comm,ierr)
        if(ierr/=MPI_SUCCESS) call parallel_abort('levels0: allreduce prwt_xchng_gb',ierr)
        call mpi_allreduce(srwt_xchng,srwt_xchng_gb,1,MPI_LOGICAL,MPI_LOR,comm,ierr)
        if(ierr/=MPI_SUCCESS) call parallel_abort('levels0: allreduce srwt_xchng_gb',ierr)
!'
#ifdef INCLUDE_TIMING
        if(cwtime) wtimer(10,2)=wtimer(10,2)+mpi_wtime()-cwtmp
#endif

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
#ifdef INCLUDE_TIMING
          if(cwtime) cwtmp=mpi_wtime()
#endif
          call exchange_p3d_4(swild)
#ifdef INCLUDE_TIMING
          if(cwtime) wtimer(10,2)=wtimer(10,2)+mpi_wtime()-cwtmp
#endif
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
#ifdef INCLUDE_TIMING
          if(cwtime) cwtmp=mpi_wtime()
#endif
          call exchange_s3d_4(swild)
#ifdef INCLUDE_TIMING
          if(cwtime) wtimer(10,2)=wtimer(10,2)+mpi_wtime()-cwtmp
#endif
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
!     Given global coord. (may not be on surface of earth), find lat/lon in radian
!===============================================================================
      subroutine compute_ll(xg,yg,zg,rlon,rlat)
      use schism_glbl, only : rkind,pi,errmsg,rearth_pole,rearth_eq
      use schism_msgp, only : parallel_abort
      implicit none
      real(rkind),intent(in) :: xg,yg,zg
      real(rkind),intent(out) :: rlon,rlat
      real(rkind) :: rad

      rad=sqrt(xg*xg+yg*yg+zg*zg)
      if(rad==0._rkind.or.abs(zg)>rad) then
        write(errmsg,*)'COMPUTE_LL: rad=0:',xg,yg,zg,rad
        call parallel_abort(errmsg)
      endif

      rlon=atan2(yg,xg) !(-pi,pi]
      if(abs(rearth_pole-rearth_eq)<1.d-2) then !for backward compatibility
        rlat=asin(zg/rad)
      else
        rlat=asin(zg/rearth_pole)
      endif
 
      end subroutine compute_ll
      
!===============================================================================
!     Cross-product of two vectors: (x1,y1,z1) x (x2,y2,z2) = (x3,y3,z3)
!===============================================================================
      subroutine cross_product(x1,y1,z1,x2,y2,z2,x3,y3,z3)
      use schism_glbl, only : rkind
      implicit none
      real(rkind),intent(in) :: x1,y1,z1,x2,y2,z2
      real(rkind),intent(out) :: x3,y3,z3

      x3=y1*z2-y2*z1
      y3=x2*z1-x1*z2
      z3=x1*y2-x2*y1

      end subroutine cross_product
      

      !dir$ attributes forceinline :: signa2
      function signa2(x1,x2,x3,y1,y2,y3)
!-------------------------------------------------------------------------------
! Compute signed area formed by pts 1,2,3 (positive counter-clockwise)
!-------------------------------------------------------------------------------
      use schism_glbl, only : rkind,errmsg
      implicit none
      real(rkind) :: signa2
      real(rkind),intent(in) :: x1,x2,x3,y1,y2,y3

      signa2=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2._rkind
  
      end function signa2

!===============================================================================
!===============================================================================
