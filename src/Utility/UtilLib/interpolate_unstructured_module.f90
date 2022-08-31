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


!   Module to find weights from a background unstructured grid
!   (tri-quad). _Double precision_ version.  
!   Search for parent elements along x- or y- strips (use small mne_bin for efficiency)

!   Calling sequence
!   (1) use interpolate_unstructured_module
!   (2) get npfg, npbg,nebg from fore- and back-ground grids (fg
!      actually can be list of ponts); 
!   (3) allocate(xfg(npfg),yfg(npfg),i34bg(nebg),nmbg(4,nebg),iparen(npfg), &
!     &xbg(npbg),ybg(npbg),dpbg(npbg),arco(4,npfg))
!   (4) Obtain xfg(),yfg() from fg;
!   (5) call interpolate_unstructured_weights(bgname,small1,is_xy,nbin,mne_bin,npfg,npbg,nebg,xfg,yfg, &
!      &i34bg,nmbg,xbg,ybg,dpbg,iparen,arco)

!     Input arguments:
!       bgname: background grid name;
!       small1: small positive number used to check if a pt is inside a polygon (e.g. 1.d-4);
!       is_xy: search bins varies along x (1) or y (2) axis;
!       nbin,mne_bin: # of bins & max # of elements in each bin (e.g. 20000 34000); 
!       npfg: # of points in foreground
!       npbg,nebg: # of nodes/elements in bg grid (bgname)
!       xfg(),yfg(): coord for fg pts;

!     Output arguments:  
!       i34bg(),nmbg(:,:): elem type and connectivity table of bg grid;
!       xbg(),ybg(),dpbg(): coord and depth at bg nodes;
!       iparen(npfg): parent elem #; <0 if no parent elements are found (abs value is the not-so-precise
!                     nearest elem); 
!       arco(4,npfg): area coord if a  parent elem is found (iparen>0). 
!   (6) The interp'ed value at fg pt i can be calculated as:
!         ie=iparen(i)
!         if(ie>0) then
!           final=0
!           do j=1,i34bg(ie)
!             nd=nmbg(j,ie)
!             final=final+var(nd)*arco(j,i)
!           enddo !j
!         endif !ie
         
!     Additional outputs in this module:
!              fort.11: fatal errors.

!     ifort -Bstatic -O3 -c interpolate_unstructured_module.f90

      module interpolate_unstructured_module 
      implicit real*8(a-h,o-z)

      private

      parameter(pi=3.141592653589793)
!      parameter(small1=1.e-4) !small number used to check if a pt is inside a polygon

!        real*8,save,allocatable :: xbg(:),ybg(:),dpbg(:),areabg(:),arbg(:),xybin(:)
!        integer,save,allocatable :: nmbg(:,:),i34bg(:),ne_bin(:),ie_bin(:,:)

      public :: interpolate_unstructured_weights

      contains

      subroutine interpolate_unstructured_weights(bgname,small1,is_xy,nbin,mne_bin,npfg,npbg,nebg,xfg,yfg, &
     &i34bg,nmbg,xbg,ybg,dpbg,iparen,arco)
      character(*), intent(in) :: bgname
      integer, intent(in) :: is_xy,nbin,mne_bin,npfg,npbg,nebg
      real*8, intent(in) :: small1,xfg(npfg),yfg(npfg)

      integer, intent(out) :: i34bg(nebg),nmbg(4,nebg),iparen(npfg)
      real*8, intent(out) :: xbg(npbg),ybg(npbg),dpbg(npbg),arco(4,npfg)

      !Local
      integer :: ne_bin(nbin),ie_bin(mne_bin,nbin),nind(3)
      real*8 :: binwid,xy_min,xy_max
      real*8 :: areabg(nebg),xybin(nbin+1)

!     Check inputs
      if(is_xy/=1.and.is_xy/=2) stop 'interpolate_unstructured_module: Check is_xy'
      if(small1<=0.d0) stop 'interpolate_unstructured_module: small1'
      if(nbin<2) stop 'interpolate_unstructured_module: nbin'
!'
      ie_bin(mne_bin,nbin)=0 !test memory leak

!      open(9,file='interpolate_unstructured.in',status='old')
!!     Input x- or y- search (1: x; 2: y)
!      read(9,*) is_xy
!!     Input # of bins (larger the better; e.g. 20000), 
!!     and max # of elements in each bin (smaller the better; e.g. # of elements/20000):'
!      read(9,*) nbin,mne_bin
!
!     Set threshold for aspect ratio in background grid (suggest 8) so those triangular elements 
!     will be discarded in interpolation
!     Aspect ratio of a triangle is defined as AR=max(Li)/R, where Li are 3 sides, and R=sqrt(A/pi)
!     is equivalent radius (A is area). Note that AR>=2*sqrt(pi/sqrt(3))=2.6935
!      read(9,*) ar_bgmax
!      if(mne_bin<=0.or.ar_bgmax<2.6935) stop 'Check mne_bin or ar_bgmax'

!'    Offset in depth in m for output (h_off is added to depths at output if interpolation is done)
!      read(9,*) h_off
!      close(9)

      open(10,file=trim(adjustl(bgname)),status='old')
      read(10,*)
      read(10,*)nebg2,npbg2
      if(nebg2/=nebg.or.npbg2/=npbg) stop 'interpolate_unstructured_module: mismatch'
!'
      do i=1,npbg
        read(10,*)j,xbg(i),ybg(i),dpbg(i)
!        if(i==1) then
!          xmin_bg=xbg(1); xmax_bg=xbg(1)
!          ymin_bg=ybg(1); ymax_bg=ybg(1)
!        else
!          if(xbg(i)<xmin_bg) xmin_bg=xbg(i)
!          if(xbg(i)>xmax_bg) xmax_bg=xbg(i)
!          if(ybg(i)<ymin_bg) ymin_bg=ybg(i)
!          if(ybg(i)>ymax_bg) ymax_bg=ybg(i)
!        endif
      enddo !i=1,npbg

!     Nudge min/max values to avoid underflow
      xmin_bg=minval(xbg)-1.e-5
      xmax_bg=maxval(xbg)+1.e-5
      ymin_bg=minval(ybg)-1.e-5
      ymax_bg=maxval(ybg)+1.e-5

!      arbg=0 !A.R. for quad
      do i=1,nebg
        read(10,*)j,i34bg(i),(nmbg(k,i),k=1,i34bg(i))
        if(i34bg(i)/=3.and.i34bg(i)/=4) then
          write(*,*)'interpolate_unstructured_module: only triangles/quad allowed:',i
!'
          stop
        endif
        n1=nmbg(1,i)
        n2=nmbg(2,i)
        n3=nmbg(3,i)
        areabg(i)=signa(xbg(n1),xbg(n2),xbg(n3),ybg(n1),ybg(n2),ybg(n3))
        if(areabg(i)<=0) then
          write(*,*)'interpolate_unstructured_module:Concave element:',i
          stop
        endif
        if(i34bg(i)==4) then
          n4=nmbg(4,i)
          tmp=signa(xbg(n1),xbg(n3),xbg(n4),ybg(n1),ybg(n3),ybg(n4))
          if(tmp<=0) then
            write(*,*)'interpolate_unstructured_module:Concave element (2):',i
!'
            stop
          endif
          areabg(i)=areabg(i)+tmp
        endif
      enddo !i=1,nebg
      close(10)

!     Bucket strip sorting
!     If an element i is in ie_bin(:,l), it's 'physically' in it (i.e. at least one internal point
!     is inside the bin l; 1<=l<=nbin).
!     If a given pt is inside bin l, search ie_bin(:,l); in addition, also search neighboring 
!     (if any) bins if it's right on the border
      ne_bin=0
!      write(98,*)'Max/min:',xmin_bg,xmax_bg,ymin_bg,ymax_bg
      do i=1,nbin+1
        xbin=xmin_bg+(xmax_bg-xmin_bg)/nbin*(i-1)
        ybin=ymin_bg+(ymax_bg-ymin_bg)/nbin*(i-1)
        if(is_xy==1) then
          xybin(i)=xbin
        else
          xybin(i)=ybin
        endif
!        write(98,*)i,xybin(i)
      enddo !i
      binwid=(xybin(nbin+1)-xybin(1))/nbin !bin width
      if(binwid<=0) stop 'interpolate_unstructured_module:negative bin width'

!'    Set limits
      if(is_xy==1) then
        xy_min=xmin_bg; xy_max=xmax_bg
      else
        xy_min=ymin_bg; xy_max=ymax_bg
      endif

      iabort=0 !flag for success in binning elements
      do i=1,nebg
!        if(arbg(i)>=ar_bgmax) cycle

!        do j=1,i34bg(i)
!          nd=nmbg(i,j)
!          if(j==1) then
!            xe_min=xbg(nd); xe_max=xbg(nd)
!            ye_min=ybg(nd); ye_max=ybg(nd)
!          else
!            if(xbg(nd)<xe_min) xe_min=xbg(nd)
!            if(xbg(nd)>xe_max) xe_max=xbg(nd)
!            if(ybg(nd)<ye_min) ye_min=ybg(nd)
!            if(ybg(nd)>ye_max) ye_max=ybg(nd)
!          endif
!        enddo !j

        xe_min=minval(xbg(nmbg(1:i34bg(i),i)))
        ye_min=minval(ybg(nmbg(1:i34bg(i),i)))
        xe_max=maxval(xbg(nmbg(1:i34bg(i),i)))
        ye_max=maxval(ybg(nmbg(1:i34bg(i),i)))
        if(is_xy==1) then
          bmin=xe_min; bmax=xe_max
        else
          bmin=ye_min; bmax=ye_max
        endif
        ibin_min=min(nbin,int((bmin-xybin(1))/binwid+1)) !estimate

        ibin1=0 !start bin #
        ibin2=0 !end bin #
        do l=ibin_min,nbin !min0(ibin_max+1,nbin)
          if(bmin>=xybin(l).and.bmin<xybin(l+1)) ibin1=l
          if(bmax>=xybin(l).and.bmax<=xybin(l+1)) ibin2=l
          if(ibin1/=0.and.ibin2/=0) exit
        enddo !l

        if(ibin1==0.or.ibin2==0) then
          write(11,*)'Cannot find a bin:',i,ibin1,ibin2,bmin,bmax
          do l=1,nbin+1
            write(11,*)xybin(l)
          enddo !l
          stop
        endif
        if(ibin1>ibin2) then
          write(*,*)'interpolate_unstructured_module: Reversed bin'
          stop
        endif
        do l=ibin1,ibin2
          ne_bin(l)=ne_bin(l)+1
          if(ne_bin(l)>mne_bin) then
            if(iabort==0) print*, 'Aborting...'
            iabort=1
          endif
          if(iabort==0) ie_bin(ne_bin(l),l)=i
        enddo !l
      enddo !i=1,nebg

      mne_bin2=maxval(ne_bin)
      if(iabort==0) then
        print*, 'done bucket sorting; actual max. elements in a bin = ',mne_bin2
      else !failed
        write(*,*)'Falied in bucket sort; max. elements in a bin needs to be ',mne_bin2
!'
        stop
      endif
!     end bucket sort

!     Find parent elem
      iparen=0
      do i=1,npfg
        x0=xfg(i)
        y0=yfg(i)
        if(is_xy==1) then
          xy=xfg(i)
        else
          xy=yfg(i)
        endif
       
        if(xy<xy_min) then
          ibin1=1; ibin2=1
        else if(xy>xy_max) then
          ibin1=nbin; ibin2=nbin
        else !in [xy_min,xy_max]
          l=min(nbin,int((xy-xybin(1))/binwid+1))
          if(xy==xybin(l)) then
            ibin1=max(l-1,1); ibin2=l
          else if(xy==xybin(l+1)) then
            ibin1=l; ibin2=min(l+1,nbin)
          else if(xy>xybin(l).and.xy<xybin(l+1)) then
            ibin1=l; ibin2=l
          else
            write(11,*)'Cannot find a bin (2):',node_num,x0,y0,xybin(l),xybin(l+1),binwid,l
            write(11,*)(m,xybin(m),m=1,nbin+1)
            stop
          endif
        endif !xy

        rmin_ac=huge(1.d0) !min area coord sum -1 to find nearest elem
        loop1: do l=ibin1,ibin2
          do k=1,ne_bin(l)
            ie=ie_bin(k,l)
            suma=0
            do j=1,3 !compute 3 area coordinates of tri (1,2,3)
              j_1=j+1
              j_2=j+2
              if(j_1>3) j_1=j_1-3
              if(j_2>3) j_2=j_2-3
              n1=nmbg(j_1,ie)
              n2=nmbg(j_2,ie)
              n3=nmbg(j,ie)
              if(j==1) area2=abs(signa(xbg(n1),xbg(n2),xbg(n3),ybg(n1),ybg(n2),ybg(n3)))
              arco(j,i)=abs(signa(xbg(n1),xbg(n2),x0,ybg(n1),ybg(n2),y0))/area2
              suma=suma+arco(j,i)
            enddo !j
            tmp=abs(suma-1.d0)
            if(tmp<small1) then
              iparen(i)=ie
              arco(4,i)=0
              exit loop1
            else if(tmp<rmin_ac) then
              rmin_ac=tmp
              iparen(i)=-ie
            endif

            if(i34bg(ie)==4) then 
              nind(1)=1; nind(2)=3; nind(3)=4
              suma=0
              do j=1,3 !compute 3 area coordinates of tri (1,3,4)
                j_1=j+1
                j_2=j+2
                if(j_1>3) j_1=j_1-3
                if(j_2>3) j_2=j_2-3
                n1=nmbg(nind(j_1),ie)
                n2=nmbg(nind(j_2),ie)
                n3=nmbg(nind(j),ie)
                if(j==1) area2=abs(signa(xbg(n1),xbg(n2),xbg(n3),ybg(n1),ybg(n2),ybg(n3)))
                arco(nind(j),i)=abs(signa(xbg(n1),xbg(n2),x0,ybg(n1),ybg(n2),y0))/area2
                suma=suma+arco(nind(j),i)
              enddo !j

              tmp=abs(suma-1.d0)
              if(tmp<small1) then
                iparen(i)=ie
                arco(2,i)=0
                exit loop1
              else if(tmp<rmin_ac) then
                rmin_ac=tmp
                iparen(i)=-ie
              endif
            endif !quad
          enddo !k=1,ne_bin(l)
        end do loop1 !l
      enddo !i=1,npfg

      end subroutine interpolate_unstructured_weights

      function signa(x1,x2,x3,y1,y2,y3)
!...  Compute signed area formed by pts 1,2,3

      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2
      
      return
      end

      end module interpolate_unstructured_module
