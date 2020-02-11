!     Given a list of pts, find parent elements in hgrid.gr3 (tri or quad), and generate source_sink.in etc 
!     Search for parent elements along x- or y- strips (use small mne_bin for efficiency)
!     Works for mixed tri/quad
!     Inputs: ss_pt.bp (list of pts, with z related to volume rate (positive for sinks); see line 217 onwards)
!             hgrid.gr3 
!             gen_source.in (see sample below): 
!                                         (1) is_xy: search along x (1) or y (2);
!                                         (2) nbin,mne_bin: # of bins & max # of elements in each bin; 
!     Output:  source_sink.in, vsource.th, vsink.th, msource.th.
!              If no parent elements are found, the pts are excluded. If >1 pt are in same elem, 
!              the results are combined.
!              fort.11: fatal errors.
!
!     ifort -O2 -mcmodel=medium -CB -Bstatic -o gen_source gen_source.f90

!     Sample interpolate_unstructured.in
!     1 !is_xy
!     20000 34000 !nbin,mne_bin

      module global
        implicit real*8(a-h,o-z)

        parameter(pi=3.141592653589793)
        parameter(small1=1.e-4) !small number used to check if a pt is inside a polygon
        integer,save :: mne_bin !max # of elements in each bin 

        integer,save :: is_xy,nbin,nebg,npbg
        real*8,save :: binwid,ar_bgmax,xy_min,xy_max

        real*8,save,allocatable :: xbg(:),ybg(:),dpbg(:),areabg(:),arbg(:),xybin(:)
        integer,save,allocatable :: nmbg(:,:),i34bg(:),ne_bin(:),ie_bin(:,:)
      end module global

      program cross
      use global
      implicit real*8(a-h,o-z)
      dimension arco(4),nmfg(4)
      integer,allocatable :: iparen(:),imap_source(:),imap_sink(:),ifl_el(:)
      real*8,allocatable :: dpfg(:),vrate(:),vsource(:),vsink(:),msource(:),total_r(:)

      open(9,file='gen_source.in',status='old')
!     Input x- or y- search (1: x; 2: y)
      read(9,*) is_xy
      if(is_xy/=1.and.is_xy/=2) stop 'Check is_xy'
!     Input # of bins (larger the better; e.g. 20000), 
!     and max # of elements in each bin (smaller the better; e.g. # of elements/20000):'
      read(9,*) nbin,mne_bin

!     Set threshold for aspect ratio in background grid (suggest 8) so those triangular elements 
!     will be discarded in interpolation
!     Aspect ratio of a triangle is defined as AR=max(Li)/R, where Li are 3 sides, and R=sqrt(A/pi)
!     is equivalent radius (A is area). Note that AR>=2*sqrt(pi/sqrt(3))=2.6935
!      read(9,*) ar_bgmax
!      if(mne_bin<=0.or.ar_bgmax<2.6935) stop 'Check mne_bin or ar_bgmax'
!     Offset in depth in m for output (h_off is added to depths at output if interpolation is done)
!      read(9,*) h_off
      close(9)

      open(10,file='hgrid.gr3',status='old')
      read(10,*)
      read(10,*)nebg,npbg
      allocate(xbg(npbg),ybg(npbg),dpbg(npbg),nmbg(nebg,4),i34bg(nebg),areabg(nebg),arbg(nebg), &
     &xybin(nbin+1),ne_bin(nbin),ie_bin(nbin,mne_bin),stat=istat)
      ie_bin(nbin,mne_bin)=0 !test memory leak
      if(istat/=0) stop 'Allocation failed (1)'
      do i=1,npbg
        read(10,*)j,xbg(i),ybg(i),dpbg(i)
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

      !arbg=0 !init. A.R. for quad
      do i=1,nebg
        read(10,*)j,i34bg(i),(nmbg(i,k),k=1,i34bg(i))
        if(i34bg(i)/=3.and.i34bg(i)/=4) then
          write(*,*)'Only triangles/quad allowed:',i
          stop
        endif
        n1=nmbg(i,1)
        n2=nmbg(i,2)
        n3=nmbg(i,3)
        areabg(i)=signa(xbg(n1),xbg(n2),xbg(n3),ybg(n1),ybg(n2),ybg(n3))
        if(areabg(i)<=0) then
          write(*,*)'Concave element:',i
          stop
        endif
        if(i34bg(i)==4) then
          n4=nmbg(i,4)
          tmp=signa(xbg(n1),xbg(n3),xbg(n4),ybg(n1),ybg(n3),ybg(n4))
          if(tmp<=0) then
            write(*,*)'Concave element (2):',i
            stop
          endif
          areabg(i)=areabg(i)+tmp
        endif
       
        !Aspect ratio
        !if(i34bg(i)==3) then
        !  rl1=sqrt((xbg(n2)-xbg(n3))**2+(ybg(n2)-ybg(n3))**2)
        !  rl2=sqrt((xbg(n1)-xbg(n3))**2+(ybg(n1)-ybg(n3))**2)
        !  rl3=sqrt((xbg(n2)-xbg(n1))**2+(ybg(n2)-ybg(n1))**2)
        !  arbg(i)=max(rl1,rl2,rl3)/sqrt(areabg(i)/pi)
        !endif !i34bg(i)==3
      enddo !i=1,nebg
      close(10)

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
        !if(arbg(i)>=ar_bgmax) cycle

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
      open(12,file='ss_pt.bp',status='old')
      read(12,*)
      read(12,*)npfg
      allocate(iparen(npfg),dpfg(npfg),vrate(npfg),ifl_el(nebg),total_r(nebg), &
     &imap_sink(nebg),imap_source(nebg),vsource(nebg),vsink(nebg),msource(nebg))
!     Interpolate and output
      total_r=0
      ifl_el=0 !flag for each bg elem.
      do i=1,npfg
        read(12,*)j,xfg,yfg,vrate(i) !vrate: m/year (positive ->sink)
        if(vrate(i)==0) cycle

        call binsearch(i,xfg,yfg,iparen(i),arco)

        if(iparen(i)/=0) then
          ie=iparen(i)
          ifl_el(ie)=1
          dpfg(i)=0
          do j=1,i34bg(ie)
            nd=nmbg(ie,j)
            dpfg(i)=dpfg(i)+dpbg(nd)*arco(j)
          enddo !j
          total_r(ie)=total_r(ie)+vrate(i)
        endif
!        write(13,*)i,real(xfg),real(yfg),real(dpfg),iparen(i)
      enddo !i=1,npfg
      close(12)

      nsource=0
      nsink=0
      do i=1,nebg
        if(total_r(i)>0) then !if positive values are source, change here
          nsink=nsink+1
          imap_sink(nsink)=i
        endif
        if(total_r(i)<0) then
          nsource=nsource+1
          imap_source(nsource)=i
        endif
      enddo !i

!     Process source outputs
      open(12,file='source_sink.in',status='replace')
      write(12,*)nsource
      do i=1,nsource
        write(12,*)imap_source(i)
      enddo !i
      write(12,*)
      write(12,*)nsink
      do i=1,nsink
        write(12,*)imap_sink(i)
      enddo !i
      close(12)

!     Const. in time for .th
      if(nsource>0) then
        open(13,file='vsource.th',status='replace')
        do i=1,nsource
          ie=imap_source(i)
          dp=max(1.d0,sum(dpbg(nmbg(ie,1:i34bg(ie))))/i34bg(ie))
          rl=100 !100m as a cell? !sqrt(areabg(ie)/pi)
          vsource(i)=abs(total_r(ie))*dp*rl/365/86400 !m^3/s
        enddo !i
        print*, 'Max source=',maxval(vsource(1:nsource))
        write(13,'(1000000(1x,e12.4))')0.,vsource(1:nsource)
        write(13,'(1000000(1x,e12.4))')86400*1000.,vsource(1:nsource)
        close(13)

        open(16,file='msource.th',status='replace')
        ntracers=3
        write(16,'(1000000(1x,e12.4))')0.,(10.,0.,2650.,2650.,2650.,i=1,nsource)
        write(16,'(1000000(1x,e12.4))')86400*1000.,(10.,0.,2650.,2650.,2650.,i=1,nsource)
        close(16)
      endif !nsource

      if(nsink>0) then
        open(13,file='vsink.th',status='replace')
        do i=1,nsink
          ie=imap_sink(i)
          dp=max(1.0d0,sum(dpbg(nmbg(ie,1:i34bg(ie))))/i34bg(ie))
          rl=100 !sqrt(areabg(ie)/pi)
          vsink(i)=-abs(total_r(ie))*dp*rl/365/86400 !m^3/s
        enddo !i
        print*, 'Max sink=',minval(vsink(1:nsink))

        write(13,'(1000000(1x,e12.4))')0.,vsink(1:nsink)
        write(13,'(1000000(1x,e12.4))')86400*1000.,vsink(1:nsink)
        close(13)
      endif !nsink

      stop
      end

!     Search for parent element of (x0,y0) and compute area coordinates arco(4)
!     If the parent element cannot be found, iparen=0
!     Sum of arco may not be exactly 1
      subroutine binsearch(node_num,x0,y0,iparen,arco)
      use global
      implicit real*8(a-h,o-z)
      integer, intent(in) :: node_num !info only
      real*8, intent(in) :: x0,y0
      integer, intent(out) :: iparen
      real*8, intent(out) :: arco(4)
      dimension nind(3)

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
        write(*,*)'Cannot find a bin (2):',node_num,x0,y0,xybin(l),xybin(l+1),binwid,l
        write(99,*)(i,xybin(i),i=1,nbin+1)
        stop
      endif

      loop1: do l=ibin1,ibin2
        do k=1,ne_bin(l)
          ie=ie_bin(l,k)
          suma=0
          do j=1,3
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
            do j=1,3
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

      return
      end

      function signa(x1,x2,x3,y1,y2,y3)
!...  Compute signed area formed by pts 1,2,3
      implicit real*8(a-h,o-z)

      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2
      
      return
      end
