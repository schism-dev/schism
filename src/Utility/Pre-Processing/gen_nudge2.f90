!     Generate [MOD]_nudge.gr3 based on a distance from open bnd
!     segments
!     Works for mixed grids
!     Input: 
!           (1) hgrid.gr3 (with bnd part)
!           (2) screen
!     Output: nudge.gr3 (edit out unwanted portions in gredit)

!     ifort -mcmodel=medium -O2 -o gen_nudge2 gen_nudge2.f90

      implicit real*8(a-h,o-z)
      parameter(nbyte=4)
      allocatable :: xnd(:),ynd(:),dp(:),nond(:),iond(:,:)
      allocatable :: i34(:),iobnd(:),list(:),xctr(:),yctr(:)
      integer, allocatable :: elnode(:,:)
      allocatable :: icolor(:)

      print*, 'Input # of open bnd segments used in nudging:'
      read*, nobnd
      allocate(iobnd(nobnd))
      print*, 'Input segment indices for those segments:'
      read*, iobnd(:)
      print*, nobnd,iobnd
      print*, 'Input max relax distance in m or degr:'
      !Can also adjust this for individual segments
      read*, rlmax
      print*, 'Input max relax strength in days:'
      read*, rnu_day
  
      open(14,file='hgrid.gr3',status='old')
      open(12,file='nudge.gr3',status='replace')
      read(14,*); read(14,*)ne,np
      allocate(xnd(np),ynd(np),dp(np),i34(ne),elnode(4,ne),xctr(ne),yctr(ne), &
     &icolor(np),stat=istat)
      if(istat/=0) stop 'Failed to alloc (1)'
      do i=1,np
        read(14,*)j,xnd(i),ynd(i),dp(i)
      enddo !i
      do i=1,ne
        read(14,*)j,k,elnode(1:k,i)
        i34(i)=k
        xctr(i)=sum(xnd(elnode(1:k,i)))/k
        yctr(i)=sum(ynd(elnode(1:k,i)))/k
      enddo !i
      read(14,*)nope
      read(14,*)neta
      nlines_b4=2+np+ne+2 !for rewind and re-read the open bnd part
      allocate(nond(nope))
      do i=1,nope
        read(14,*)nond(i)
        do j=1,nond(i)
          read(14,*) !iond(j)
        enddo !j
      enddo !i
      mnond=maxval(nond(:))
      nond_all=sum(nond(:))
      nlist=sum(nond(iobnd(:)))
      allocate(iond(nope,mnond),list(nlist))
      rewind(14)
      do l=1,nlines_b4
        read(14,*)
      enddo !l

      icount=0
      do i=1,nope
        read(14,*) !nond(i)
        do j=1,nond(i)
          read(14,*) iond(i,j)
          icount=icount+1
          !write(12,*)icount,xnd(iond(i,j)),ynd(iond(i,j)),i
        enddo !j
      enddo !i
      close(14)

      if(icount/=nond_all) stop 'Mismatch' 

      !Populate list(): list of all open bnd nodes used in nudging
      icount=0
      icolor=0 !flag
      do i2=1,nobnd
        i=iobnd(i2) !segment index
        do j=1,nond(i)
          icount=icount+1
          nd0=iond(i,j)
          list(icount)=nd0
          icolor(nd0)=1
        enddo !j
      enddo !i2

      if(icount/=nlist) then
        print*, 'Mismatch:',icount,nlist
        stop
      endif
!      print*, 'Final list of open bnd nodes:',list

      !Generate nudging zones; adjust here to use different values for
      !each seg
      write(12,*)rlmax,rnu_day
      write(12,*)ne,np
      rnu_max=1./rnu_day/86400.

      do id=1,np
        if(icolor(id)==1) then
          rnu=rnu_max
          distmin=0; in0=id 
        else !not on open seg's
          distmin=huge(1.d0) !min distance
          in0=0 !node index for min distance
          do i2=1,nobnd
            i=iobnd(i2) !segment index
            do j=1,nond(i)
              nd0=iond(i,j)
              rl2=sqrt((xnd(id)-xnd(nd0))**2+(ynd(id)-ynd(nd0))**2)
              !if(id==1.and.i2==2) write(99,*)j,nd0,rl2
              if(rl2<distmin) then
                distmin=rl2; in0=nd0 
              endif !rl2
            enddo !j
          enddo !i2

          rnu=0
          if(distmin<=rlmax) rnu=(1-distmin/rlmax)*rnu_max !\in [0,1]
        endif !icolor

        !if(mod(id,100)==0) print*, 'done node ',id
        write(12,*)id,real(xnd(id)),real(ynd(id)),real(rnu),in0 !real(distmin),in0
      enddo !id

      do i=1,ne
        write(12,*)i,i34(i),elnode(1:i34(i),i) 
      enddo !i

      stop
      end 
