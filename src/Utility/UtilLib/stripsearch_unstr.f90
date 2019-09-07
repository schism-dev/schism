!     Sort each unstru. grid element into strips along x or y to speed up interpolation later
!     In:
!            is_xy - 1: search along x; 2: along y
!            nbin - # of bins along x/y
!            mne_bin - max. # of elem. per bin
!            nebg,npbg - # of elem. and nodes for the unstr. grid (background grid)
!            xbg,ybg,i34,nmbg - x,y, and connectivity table for unstr. grid 
!     Out:
!            ne_bin(1:nbin) - # of elem. for each bin
!            ie_bin(1:nbin,1:ne_bin) - list of elem. for each bin
!            xybin(1:nbin+1) - x or y coord. of each bin line

      subroutine stripsearch_unstr(is_xy,nbin,mne_bin,nebg,npbg,xbg,ybg,i34,nmbg, &
     &ne_bin,ie_bin,xybin,binwid)
      implicit real(8)(a-h,o-z)

      integer, intent(in) :: is_xy,nbin,mne_bin,nebg,npbg,i34(nebg),nmbg(nebg,4)
      real(8), intent(in) :: xbg(npbg),ybg(npbg)
      integer, intent(out) :: ne_bin(nbin),ie_bin(nbin,mne_bin)
      real(8), intent(out) :: xybin(nbin+1),binwid
      
      if(is_xy/=1.and.is_xy/=2.or.mne_bin<=0.or.nbin<=0) stop 'STRIPSEARCH: Check is_xy etc'
!'
      ie_bin(nbin,mne_bin)=0 !test memory leak
      xmin_bg=minval(xbg)-1.e-2
      xmax_bg=maxval(xbg)+1.e-2
      ymin_bg=minval(ybg)-1.e-2
      ymax_bg=maxval(ybg)+1.e-2

!     Bucket strip sorting
!     If an element i is in ie_bin(l,:), it's 'physically' in it (i.e. at least one internal point
!     is inside the bin l; 1<=l<=nbin).
!     If a given pt is inside bin l, search ie_bin(l,:); in addition, also search neighboring 
!     (if any) bins if it's right on the border
      ne_bin=0
      write(95,*)'Max/min:',xmin_bg,xmax_bg,ymin_bg,ymax_bg,'; coord. of bin lines'
      do i=1,nbin+1
        xbin=xmin_bg+(xmax_bg-xmin_bg)/nbin*(i-1)
        ybin=ymin_bg+(ymax_bg-ymin_bg)/nbin*(i-1)
        if(is_xy==1) then
          xybin(i)=xbin
        else
          xybin(i)=ybin
        endif
        write(95,*)i,xybin(i)
      enddo !i
      binwid=(xybin(nbin+1)-xybin(1))/nbin !bin width
      if(binwid<=0) stop 'STRIPSEARCH: negative bin width'

!      print*, 'min/max=',xmin_bg,xmax_bg,ymin_bg,ymax_bg,xybin(1),xybin(nbin+1),binwid

!     Set limits
!      if(is_xy==1) then
!        xy_min=xmin_bg; xy_max=xmax_bg
!      else
!        xy_min=ymin_bg; xy_max=ymax_bg
!      endif

      iabort=0 !flag for success in binning elements
      do i=1,nebg
        if(is_xy==1) then
          bmin=minval(xbg(nmbg(i,:))); bmax=maxval(xbg(nmbg(i,1:i34(i))))
        else
          bmin=minval(ybg(nmbg(i,:))); bmax=maxval(ybg(nmbg(i,1:i34(i))))
        endif
        ibin_min=(bmin-xybin(1))/binwid+1 !estimate
        ibin_min=min(nbin,max(1,ibin_min))
!        ibin_max=min(nbin,int((bmax-xybin(1))/binwid+1))

        ibin1=0 !start bin #
        ibin2=0 !end bin #
        do l=ibin_min,nbin !min0(ibin_max+1,nbin)
          if(bmin>=xybin(l).and.bmin<xybin(l+1)) ibin1=l
          if(bmax>=xybin(l).and.bmax<=xybin(l+1)) ibin2=l
          if(ibin1/=0.and.ibin2/=0) exit
        enddo !l

        if(ibin1==0.or.ibin2==0) then
          write(11,*)'STRIPSEARCH: Cannot find a bin:',i,ibin1,ibin2,bmin,bmax,ibin_min,binwid
          do l=1,nbin+1
            write(11,*)l,xybin(l)
          enddo !l
          stop
        endif
        if(ibin1>ibin2) then
          write(*,*)'STRIPSEARCH: Reversed bin'
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
        print*, 'STRIPSEARCH: done bucket sorting; actual max. elements in a bin = ',mne_bin2
      else !failed
        write(*,*)'STRIPSEARCH: Failed in bucket sort; max. elements in a bin needs to be ',mne_bin2
!'
        stop
      endif

      end subroutine stripsearch_unstr

