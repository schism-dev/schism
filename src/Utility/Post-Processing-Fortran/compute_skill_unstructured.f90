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

!   Compute error metrics using unstructured as background grid (tri or quad; e.g. from model), by 
!   interpolating random points (e.g. from obs) onto the UG. 
!   Main routine for outside driver is compute_skill_unstructured2() (see below)

      module compute_skill_unstructured
      implicit real*8(a-h,o-z)
      private

      parameter(pi=3.141592653589793d0)
      parameter(small1=1.d-3) !small number used to check if a pt is inside a polygon
!      integer,save :: mne_bin !max # of elements in each bin 
!      integer,save :: is_xy,nbin,nebg,npbg
      real*8,save :: binwid,ar_bgmax,xy_min,xy_max

      real*8,save,allocatable :: areabg(:),arbg(:),xybin(:)
      integer,save,allocatable :: ne_bin(:),ie_bin(:,:)

      public :: compute_skill_unstructured2

      contains 

!     Main routine
!     Search for parent elements along x- or y- strips (use small mne_bin for efficiency)
!     Input arguments:
!     is_xy: search bins varies along x (1) or y (2) axis;
!     nbin,mne_bin: # of bins & max # of elements in each bin; 
!     nebg,npbg,xbg,ybg,dpbg,nmbg,i34bg: info from bg grid (dpbg is 'modeled' result)
!     npfg,xfg,yfg,vfg: fg points and value (vfg; e.g. obs)
!     v_range(1:2): min/max of valid range (to check both model and obs)
!     Output arguments: scores (1: bias; 2: RMSE; 3: MAE); exception is large negative #
!     Additional output file: fort.11 (fatal errors)
      subroutine compute_skill_unstructured2(is_xy,nbin,mne_bin,nebg,npbg,xbg,ybg,dpbg,nmbg,i34bg, &
     &npfg,xfg,yfg,vfg,v_range,scores)
      integer, intent(in) :: is_xy,nbin,mne_bin,nebg,npbg,nmbg(nebg,4),i34bg(nebg),npfg
      real*8, intent(in) :: xbg(npbg),ybg(npbg),dpbg(npbg),xfg(npfg),yfg(npfg),vfg(npfg), &
     &v_range(2)
      real*8, intent(out) :: scores(10)

      real*8 :: arco(4),dpfg(npfg)
      integer :: indx(npfg)

      !init
      scores=-huge(1.d0)

!      open(9,file='interpolate_unstructured.in',status='old')
!!     Input x- or y- search (1: x; 2: y)
!      read(9,*) is_xy
      if(is_xy/=1.and.is_xy/=2) stop 'Check is_xy'
!     Input # of bins (larger the better; e.g. 20000), 
!     and max # of elements in each bin (smaller the better; e.g. # of elements/20000):'
!      read(9,*) nbin,mne_bin
!
!     Set threshold for aspect ratio in background grid (suggest 8) so those triangular elements 
!     will be discarded in interpolation
!     Aspect ratio of a triangle is defined as AR=max(Li)/R, where Li are 3 sides, and R=sqrt(A/pi)
!     is equivalent radius (A is area). Note that AR>=2*sqrt(pi/sqrt(3))=2.6935
!      read(9,*) ar_bgmax
!      if(mne_bin<=0.or.ar_bgmax<2.6935) stop 'Check mne_bin or ar_bgmax'
!
!     Offset in depth in m for output (h_off is added to depths at output if interpolation is done)
!      read(9,*) h_off
!      close(9)

!      open(10,file='bg.gr3',status='old')
!      read(10,*)
!      read(10,*)nebg,npbg
      allocate(areabg(nebg),arbg(nebg),xybin(nbin+1),ne_bin(nbin),ie_bin(nbin,mne_bin),stat=istat)
      ie_bin(nbin,mne_bin)=0 !test memory leak
      if(istat/=0) stop 'Inside score routine: allocation failed (1)'
      do i=1,npbg
!        read(10,*)j,xbg(i),ybg(i),dpbg(i)
        if(i==1) then
          xmin_bg=xbg(1); xmax_bg=xbg(1)
          ymin_bg=ybg(1); ymax_bg=ybg(1)
        else
          if(xbg(i)<xmin_bg) xmin_bg=xbg(i)
          if(xbg(i)>xmax_bg) xmax_bg=xbg(i)
          if(ybg(i)<ymin_bg) ymin_bg=ybg(i)
          if(ybg(i)>ymax_bg) ymax_bg=ybg(i)
        endif
      enddo !i=1,npbg

!     Nudge min/max values to avoid underflow
      xmin_bg=xmin_bg-1.e-2
      xmax_bg=xmax_bg+1.e-2
      ymin_bg=ymin_bg-1.e-2
      ymax_bg=ymax_bg+1.e-2

      arbg=0 !A.R. for quad
      do i=1,nebg
!        read(10,*)j,i34bg(i),(nmbg(i,k),k=1,i34bg(i))
        if(i34bg(i)/=3.and.i34bg(i)/=4) then
          write(*,*)'Only triangles/quad allowed:',i
          stop
        endif
        n1=nmbg(i,1)
        n2=nmbg(i,2)
        n3=nmbg(i,3)
        areabg(i)=signa(xbg(n1),xbg(n2),xbg(n3),ybg(n1),ybg(n2),ybg(n3))
        if(areabg(i)<=0) then
          write(*,*)'Concave element:',i,areabg(i),nmbg(i,1:3)
          stop
        endif
        if(i34bg(i)==4) then
          n4=nmbg(i,4)
          tmp=signa(xbg(n1),xbg(n3),xbg(n4),ybg(n1),ybg(n3),ybg(n4))
          if(tmp<=0) then
            write(*,*)'Concave element (2):',i,tmp,nmbg(i,:)
            stop
          endif
          areabg(i)=areabg(i)+tmp
        endif
       
        !Aspect ratio
        if(i34bg(i)==3) then
          rl1=sqrt((xbg(n2)-xbg(n3))**2+(ybg(n2)-ybg(n3))**2)
          rl2=sqrt((xbg(n1)-xbg(n3))**2+(ybg(n1)-ybg(n3))**2)
          rl3=sqrt((xbg(n2)-xbg(n1))**2+(ybg(n2)-ybg(n1))**2)
          arbg(i)=max(rl1,rl2,rl3)/sqrt(areabg(i)/pi)
        endif !i34bg(i)==3
      enddo !i=1,nebg
!      close(10)

!     Bucket strip sorting
!     If an element i is in ie_bin(l,:), it's 'physically' in it (i.e. at least one internal point
!     is inside the bin l; 1<=l<=nbin).
!     If a given pt is inside bin l, search ie_bin(l,:); in addition, also search neighboring 
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
      if(binwid<=0) stop 'negative bin width'

!     Set limits
      if(is_xy==1) then
        xy_min=xmin_bg; xy_max=xmax_bg
      else
        xy_min=ymin_bg; xy_max=ymax_bg
      endif

      iabort=0 !flag for success in binning elements
      do i=1,nebg
!        if(arbg(i)>=ar_bgmax) cycle

        do j=1,i34bg(i)
          nd=nmbg(i,j)
          if(j==1) then
            xe_min=xbg(nd); xe_max=xbg(nd)
            ye_min=ybg(nd); ye_max=ybg(nd)
          else
            if(xbg(nd)<xe_min) xe_min=xbg(nd)
            if(xbg(nd)>xe_max) xe_max=xbg(nd)
            if(ybg(nd)<ye_min) ye_min=ybg(nd)
            if(ybg(nd)>ye_max) ye_max=ybg(nd)
          endif
        enddo !j
        if(is_xy==1) then
          bmin=xe_min; bmax=xe_max
        else
          bmin=ye_min; bmax=ye_max
        endif
        ibin_min=min(nbin,int((bmin-xybin(1))/binwid+1)) !estimate
!        ibin_max=min(nbin,int((bmax-xybin(1))/binwid+1))

        ibin1=0 !start bin #
        ibin2=0 !end bin #
        do l=ibin_min,nbin !min0(ibin_max+1,nbin)
          if(bmin>=xybin(l).and.bmin<xybin(l+1)) ibin1=l
          if(bmax>=xybin(l).and.bmax<=xybin(l+1)) ibin2=l
          if(ibin1/=0.and.ibin2/=0) exit
        enddo !l

        if(ibin1==0.or.ibin2==0) then
          write(*,*)'Cannot find a bin; see fort.11'
          write(11,*)'Cannot find a bin:',i,ibin1,ibin2,bmin,bmax
          do l=1,nbin+1
            write(11,*)xybin(l)
          enddo !l
          stop
        endif
        if(ibin1>ibin2) then
          write(*,*)'Reversed bin'
          stop
        endif
        do l=ibin1,ibin2
          ne_bin(l)=ne_bin(l)+1
          if(ne_bin(l)>mne_bin) then
            if(iabort==0) print*, 'Aborting...'
            iabort=1
          endif
          if(iabort==0) ie_bin(l,ne_bin(l))=i
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

!     Foreground build point file
!      open(12,file='fg.gr3',status='old')
!      open(14,file='include.gr3',status='old')
!      open(13,file='fg.new',status='replace')
!      read(12,*)
!      read(12,*)nefg,npfg
!      read(14,*); read(14,*)
!      write(13,'(a30,2i6,2(1x,f9.3))')'nbin,mne_bin,ar_bgmax,h_off=',nbin,mne_bin,ar_bgmax,h_off
!      write(13,*)nefg,npfg
!     Interpolate
      dpfg=-1.e20 !init
      n_valid=0 !# of valid pts
      do i=1,npfg
        !Check obs range 
        if(vfg(i)<v_range(1).or.vfg(i)>v_range(2)) cycle

        call binsearch(i,nbin,xfg(i),yfg(i),is_xy,nebg,npbg,nmbg,i34bg,xbg,ybg,iparen,arco)
!        write(19,*)i,iparen,arco(1:4)
        if(iparen/=0) then
          dpfg(i)=0.
          do j=1,i34bg(iparen)
            nd=nmbg(iparen,j)
            dpfg(i)=dpfg(i)+dpbg(nd)*arco(j) !e.g., model value @ obs pt
          enddo !j

          !Check range
          if(dpfg(i)>=v_range(1).and.dpfg(i)<=v_range(2)) then
            n_valid=n_valid+1
            if(n_valid>npfg) stop 'Overflow(4)'
            indx(n_valid)=i
          endif
        endif !iparen/
      enddo !i=1,npfg

      !Output for debug
!      do i=1,npfg
!        write(99,'(i10,2(1x,e20.12),1x,e20.12)')i,xfg(i),yfg(i),dpfg(i)
!      enddo !i

!     Scores
      if(n_valid>0) then
        !bias
        scores(1)=sum(dpfg(indx(1:n_valid))-vfg(indx(1:n_valid)))/n_valid
        !RMSE
        scores(2)=sum((dpfg(indx(1:n_valid))-vfg(indx(1:n_valid)))**2)/n_valid
        scores(2)=sqrt(scores(2))
        !MAE
        scores(3)=sum(abs(dpfg(indx(1:n_valid))-vfg(indx(1:n_valid))))/n_valid
      endif

      deallocate(areabg,arbg,xybin,ne_bin,ie_bin)

      end subroutine compute_skill_unstructured2

!     Search for parent element of (x0,y0) and compute area coordinates arco(4)
!     If the parent element cannot be found, iparen=0
!     Sum of arco may not be exactly 1
      subroutine binsearch(node_num,nbin,x0,y0,is_xy,nebg,npbg,nmbg,i34bg,xbg,ybg,iparen,arco)
      integer, intent(in) :: node_num !info only
      integer, intent(in) :: nbin,is_xy,nebg,npbg,nmbg(nebg,4),i34bg(nebg)
      real*8, intent(in) :: x0,y0,xbg(npbg),ybg(npbg)
      integer, intent(out) :: iparen
      real*8, intent(out) :: arco(4)

      integer :: nind(3)

      if(is_xy==1) then
        xy=x0
      else
        xy=y0
      endif
     
      iparen=0
      if(xy<xy_min.or.xy>xy_max) return

      l=min(nbin,int((xy-xybin(1))/binwid+1))
      if(xy==xybin(l)) then
        ibin1=max(l-1,1); ibin2=l
      else if(xy==xybin(l+1)) then
        ibin1=l; ibin2=min(l+1,nbin)
      else if(xy>xybin(l).and.xy<xybin(l+1)) then
        ibin1=l; ibin2=l
      else
        write(*,*)'Cannot find a bin (2); see fort.11:',node_num,x0,y0,xybin(l),xybin(l+1),binwid,l
        write(11,*)(i,xybin(i),i=1,nbin+1)
        stop
      endif

      loop1: do l=ibin1,ibin2
        do k=1,ne_bin(l)
          ie=ie_bin(l,k)
          suma=0
          do j=1,3 !compute 3 area coordinates of tri (1,2,3)
            j_1=j+1
            j_2=j+2
            if(j_1>3) j_1=j_1-3
            if(j_2>3) j_2=j_2-3
            n1=nmbg(ie,j_1)
            n2=nmbg(ie,j_2)
            n3=nmbg(ie,j)
            if(j==1) area2=dabs(signa(xbg(n1),xbg(n2),xbg(n3),ybg(n1),ybg(n2),ybg(n3)))
            arco(j)=dabs(signa(xbg(n1),xbg(n2),x0,ybg(n1),ybg(n2),y0))/area2
            suma=suma+arco(j)
          enddo !j
          if(dabs(suma-1)<small1) then
            iparen=ie
            arco(4)=0
            exit loop1
          endif

          if(i34bg(ie)==4) then 
            nind(1)=1; nind(2)=3; nind(3)=4
            suma=0
            do j=1,3 !compute 3 area coordinates of tri (1,3,4)
              j_1=j+1
              j_2=j+2
              if(j_1>3) j_1=j_1-3
              if(j_2>3) j_2=j_2-3
              n1=nmbg(ie,nind(j_1))
              n2=nmbg(ie,nind(j_2))
              n3=nmbg(ie,nind(j))
              if(j==1) area2=dabs(signa(xbg(n1),xbg(n2),xbg(n3),ybg(n1),ybg(n2),ybg(n3)))
              arco(nind(j))=dabs(signa(xbg(n1),xbg(n2),x0,ybg(n1),ybg(n2),y0))/area2
              suma=suma+arco(nind(j))
            enddo !j
            if(dabs(suma-1)<small1) then
              iparen=ie
              arco(2)=0
              exit loop1
            endif
          endif !quad
        enddo !k=1,ne_bin(l)
      end do loop1 !l

!      if(iparen==0) then
!        write(*,*)'Cannot find a parent:',x0,y0
!        stop
!      endif

      end subroutine binsearch

      function signa(x1,x2,x3,y1,y2,y3)
!...  Compute signed area formed by pts 1,2,3
!      implicit real*8(a-h,o-z)

      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2
      
      end function signa

      end module compute_skill_unstructured
