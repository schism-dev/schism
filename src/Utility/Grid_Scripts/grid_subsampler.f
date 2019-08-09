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

!     Sub-sample nodes based on a minimum distance (from Mike Foreman)
!     Input: hgrid.gr3; sample_region.gr3 (1: do sampling); min. distance
!     Output: nodeflags.bp; nodesample.bp (only selected nodes)
!     ifort -Bstatic -O3 -o grid_subsampler grid_subsampler.f
      program sample
      implicit real*8(a-h,o-z)
      parameter(mnp=2200000)
!      parameter(mne=100000)
      dimension x(mnp),y(mnp),iregion(mnp),map(mnp),icolor(mnp)
      
!     Min. distance
      print*, 'Input min. distance (m):'
      read*, dis0
!      dis0=1.e3      

      open(14,file='hgrid.gr3',status='old')
      open(13,file='sample_region.gr3',status='old')
      read(14,*)
      read(14,*)ne,np
      read(13,*)
      read(13,*)
      if(np.gt.mnp) then !.or.ne.gt.mne) then
        write(*,*)'Increase mnp/mne'
        stop
      endif

      icount=0
      do i=1,np
        read(14,*)j,x(i),y(i) !,dp(i)
        read(13,*)j,xtmp,ytmp,tmp
        iregion(i)=tmp
        if(iregion(i).ne.1.and.iregion(i).ne.0) then
          write(*,*)'Unknown region flag:',i
          stop
        endif
        icolor(i)=0
        if(icount.eq.0.and.iregion(i).eq.1) then
          icount=icount+1
          map(1)=i
          icolor(i)=1
        endif
      enddo !i
      close(13)
      close(14)
      if(icount.ne.1) then
        write(*,*)'Initialization falied'
        stop
      endif
      print*, '1st node=',map(1)

      do 11 i=1,np
        if(i.eq.map(1).or.iregion(i).eq.0) go to 11
        do j=1,icount
          dis2=(x(i)-x(map(j)))**2+(y(i)-y(map(j)))**2
          if(dis2.lt.dis0*dis0) go to 11
        enddo !j
        icount=icount+1
        map(icount)=i
        icolor(i)=1
11    continue  !i=1,np 
      print*, icount,' nodes found'
      
!     Output
      open(15,file='nodeflags.bp')
      open(17,file='nodesample.bp')
      write(15,*)
      write(15,*)np
      write(17,*)
      write(17,*)icount
      do i=1,np
        write(15,101)i,x(i),y(i),icolor(i)
        if(icolor(i)==1) write(17,101)i,x(i),y(i),1
      enddo !i
101   format(i11,2(1x,e24.15),1x,i4)

      stop
      end
