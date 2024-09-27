!     Redistribute sources to limit flow so we can set S=0 at those points
!     and the output S is not too fresh. Must set lev_tr_source(2)=0 in param.nml 
!     (inject at all depths).
!     If the original source cell does not have at least 1 land bnd node, the code 
!     will not try to redistribute.
!     Works for mixed tri/quads, lon/lat or proj. The code will look for neigboring bnd cells 
!     (with at least 1 node on land bnd) to spread flow.
!     Inputs: redistribute_source.in, constants below, hgrid.gr3, source_sink.in.0, vsource.th.0
!     Output: source_sink.in,vsource.th
!             Warning messages in warning.out; fatal on screen

!     ifx -CB -O2 -g -traceback -o redistribute_source redistribute_source.f90

      implicit real*8(a-h,o-z)
      parameter(rearth_eq=6378206.4d0)
      parameter(rearth_pole=6378206.4d0)
      parameter(pi=3.1415926d0)

      real*8 :: lframe(3,3),xloc(4),yloc(4)
      real*8, allocatable :: xctr(:),yctr(:),dp(:),xnd(:),ynd(:),flow(:),av_flow(:),area(:), &
     &tot_perc(:),perc(:),x(:),y(:),znd(:),flow2(:),bbox_ll(:,:),bbox_ur(:,:)
      integer,allocatable :: i34(:),elnode(:,:),isbnd(:),ie_src(:),icolor(:), &
     &nne(:),indel(:,:),icolor2(:),nd_search(:),ie_src2(:,:),ie_sink(:)

!     Inputs
      max_iter=100 !max # of iterations for searching neighbors
!     Target input and output S (on mean flow)
      s_in=33. !if h<=100m
      s_out=30.
      !ratio s_in/s_out>1
      s_rat_shallow=s_in/s_out !if h<100m
      s_rat_deep=1./0.99 !if h>100m
      if(min(s_rat_shallow,s_rat_deep)<1.d0) stop 'S ratio<1'

!     Read in redistribute_source.in
      open(9,file='redistribute_source.in',status='old')      
      read(9,*)lonlat !>0 for lon/lat
      read(9,*)dt !Estimated max time step in .nml
      read(9,*)nbbox !# of bounding boxes in which to exclude sources
      allocate(bbox_ll(2,max(nbbox,1)),bbox_ur(2,max(nbbox,1)))
      do i=1,nbbox
        read(9,*)bbox_ll(:,i),bbox_ur(:,i) !lower-left and upper right coordinates
      enddo !i
      close(9)

!      print*, 'hgrid in lon/lat? (1: yes)'
!      read*, lonlat
!      print*, 'Estimated max time step in .nml:'
!      read*, dt

      open(14,file='hgrid.gr3',status='old')
      open(10,file='source_sink.in.0',status='old')
      open(12,file='vsource.th.0',status='old')
      open(13,file='warning.out',status='replace')
      open(15,file='source_sink.in',status='replace')
      open(17,file='vsource.th',status='replace')
      read(14,*); read(14,*)ne,np
      allocate(xctr(ne),yctr(ne),dp(np),xnd(np),ynd(np),i34(ne),elnode(4,ne),isbnd(np), &
     &icolor(ne),icolor2(np),area(ne),nne(np),tot_perc(ne),perc(ne),nd_search(np),x(np), &
     &y(np),znd(np))
      do i=1,np
        read(14,*)j,x(i),y(i),dp(i)
        if(lonlat>0) then
!          xnd(i)=xnd(i)/180.*pi
!          ynd(i)=ynd(i)/180.*pi
!          call cpp(xcpp,ycpp,xnd(i),ynd(i),cpp_lon,cpp_lat) 
!          xnd(i)=xcpp; ynd(i)=ycpp
          !global coordi.
          xtmp=x(i)/180*pi
          ytmp=y(i)/180*pi
          xnd(i)=rearth_eq*cos(ytmp)*cos(xtmp)
          ynd(i)=rearth_eq*cos(ytmp)*sin(xtmp)
          znd(i)=rearth_pole*sin(ytmp)
        endif
      enddo !i

      do i=1,ne
        read(14,*)j,k,elnode(1:k,i)
        i34(i)=k
        !Area
        n1=elnode(1,i)
        n2=elnode(2,i)
        n3=elnode(3,i)
        if(lonlat==0) then
          area(i)=signa(x(n1),x(n2),x(n3),y(n1),y(n2),y(n3))
          if(area(i)<=0.d0) then
            print*, 'Area<=0:',i,area(i)
            stop
          endif

          if(i34(i)==4) then
            n4=elnode(4,i)
            tmp=signa(x(n1),x(n3),x(n4),y(n1),y(n3),y(n4))
            if(tmp<=0.d0) then
              print*, 'Quad area<=0:',i,tmp
              stop
            endif
            area(i)=area(i)+tmp
          endif
        else !lon/lat
          !Construct local frame with 1,2 as local x-axis, and
          lframe(1,1)=xnd(n2)-xnd(n1)
          lframe(2,1)=ynd(n2)-ynd(n1)
          lframe(3,1)=znd(n2)-znd(n1)
          !(2,1)x(3,1) as local z.
          !Compute local z-axis 1st. [lframe(1:3,1:3): 1st index is
          !component]
          call cross_product(lframe(1,1),lframe(2,1),lframe(3,1), &
                              &xnd(n3)-xnd(n1),ynd(n3)-ynd(n1),znd(n3)-znd(n1), &
                              &lframe(1,3),lframe(2,3),lframe(3,3))


            !y= z \cross x
          call cross_product(lframe(1,3),lframe(2,3),lframe(3,3), &
     &lframe(1,1),lframe(2,1),lframe(3,1),lframe(1,2),lframe(2,2),lframe(3,2))

          !Make unit vectors
          !Compute local coords
          xloc(1)=0; yloc(1:2)=0
          do j=1,3
            rnorm=sqrt(lframe(1,j)**2+lframe(2,j)**2+lframe(3,j)**2)
            if(rnorm==0) stop '0 vector'
            lframe(:,j)=lframe(:,j)/rnorm
            if(j==1) xloc(2)=rnorm
          enddo !j
 
          xloc(3)=(xnd(n3)-xnd(n1))*lframe(1,1)+(ynd(n3)-ynd(n1))*lframe(2,1)+ &
     &(znd(n3)-znd(n1))*lframe(3,1)
          yloc(3)=(xnd(n3)-xnd(n1))*lframe(1,2)+(ynd(n3)-ynd(n1))*lframe(2,2)+ &
     &(znd(n3)-znd(n1))*lframe(3,2)

          area(i)=signa(xloc(1),xloc(2),xloc(3),yloc(1),yloc(2),yloc(3))
          if(area(i)<=0.d0) then
            print*, 'Area<=0(2):',i,area(i)
            stop
          endif
    
          if(i34(i)==4) then
            n4=elnode(4,i)
            xloc(4)=(xnd(n4)-xnd(n1))*lframe(1,1)+(ynd(n4)-ynd(n1))*lframe(2,1)+ &
     &(znd(n4)-znd(n1))*lframe(3,1)
            yloc(4)=(xnd(n4)-xnd(n1))*lframe(1,2)+(ynd(n4)-ynd(n1))*lframe(2,2)+ &
     &(znd(n4)-znd(n1))*lframe(3,2)
            tmp=signa(xloc(1),xloc(3),xloc(4),yloc(1),yloc(3),yloc(4))
            if(tmp<=0.d0) then
              print*, 'Quad area<=0(2):',i,tmp
              stop
            endif
            area(i)=area(i)+tmp
          endif !i34
        endif !lon/lat
      enddo !i=1,ne
      read(14,*)nope
      read(14,*)neta
      isbnd=0 !init
      do i=1,nope
        read(14,*)nond
        do j=1,nond
          read(14,*) itmp
          isbnd(itmp)=i
        enddo !j
      enddo !i
      read(14,*) nland
      read(14,*) nvel
      do k=1,nland
        read(14,*) nn
        do i=1,nn 
          read(14,*) itmp
          if(isbnd(itmp)==0) isbnd(itmp)=-1 !land bnd
        enddo;
      enddo !k
      close(14)

      !Neighborhood
      nne=0
      do i=1,ne
        do j=1,i34(i)
          nd=elnode(j,i)
          nne(nd)=nne(nd)+1
        enddo
      enddo
      mnei=maxval(nne)
      allocate(indel(mnei,np))

      nne=0
      do i=1,ne
        do j=1,i34(i)
          nd=elnode(j,i)
          nne(nd)=nne(nd)+1
          if(nne(nd)>mnei) then
            write(*,*)'Too many neighbors',nd,mnei,nne(nd)
            stop
          endif
          indel(nne(nd),nd)=i
        enddo
      enddo

      !Read source flows
      read(10,*)nsource
      if(nsource<=0.or.nsource>ne) stop 'No sources'
      allocate(ie_src(nsource),flow(nsource),av_flow(nsource),ie_src2(ne,2))
      icolor=0 !init
      do i=1,nsource
        read(10,*)ie_src(i)
        icolor(ie_src(i))=1
      enddo !i
 
      read(10,*)
      read(10,*)nsink
      allocate(ie_sink(nsink))
      do i=1,nsink
        read(10,*)ie_sink(i)
      enddo !i
      close(10)
      nsource2=nsource !init output # of sources
      ie_src2(1:nsource,1)=ie_src
      ie_src2(1:nsource,2)=ie_src

      nline=0
      av_flow=0 !mean flow @ each source
      do
        read(12,*,end=99)time,flow(1:nsource)
        av_flow=av_flow+flow
        nline=nline+1
      end do 
99    if(nline==0) stop 'vsource.th empty'
      av_flow=av_flow/nline 

      !Iteratively assess flow at each source cell
      perc(:)=1. !init ratio of flow received at each source 
      lp3: do i=1,nsource 
        ie=ie_src(i)
        nd1=elnode(1,ie) !use 1st node to check bbox

          !new11
          if(i==21) write(13,*)'CORIE (1):',ie

        do m=1,nbbox
          if(x(nd1)>=bbox_ll(1,m).and.x(nd1)<=bbox_ur(1,m).and. &
            &y(nd1)>=bbox_ll(2,m).and.y(nd1)<=bbox_ur(2,m)) then
            cycle lp3
          endif
        enddo !m

        av_h=sum(dp(elnode(1:i34(ie),ie)))/i34(ie)
        if(av_h<=0.d0) then
          write(13,*)'Warning, dry source elem in original list:',ie,av_h
          cycle lp3
        endif
        if(av_h>100.d0) then
          s_rat=s_rat_deep
        else
          s_rat=s_rat_shallow    
        endif
        vol=area(ie)*av_h !>0
        !Max flow allowed for S
        flowmax=vol/dt*(s_rat-1) !>0

          !new11
          if(i==21) write(13,*)'CORIE b4:',ie

        if(flowmax<av_flow(i)) then !look for neighboring bnd cells
 
          !new11
          if(i==21) write(13,*)'CORIE in'

          !Prep list of land bnd nodes to search         
          icolor2=0 !flag
          !Find first <=2 starting bnd nodes
          icount=0
          do j=1,i34(ie)
            nd=elnode(j,ie)
            if(isbnd(nd)==-1.and.icolor2(nd)==0) then
              icolor2(nd)=1
              icount=icount+1
              if(icount>np) then
                print*, 'Too many search nodes:',ie,icount,nd_search(1:icount-1)
                stop
              endif
              nd_search(icount)=nd
            endif
          enddo !j
          if(icount==0) then
            write(13,*)'Starting land nodes not found:',ie
            cycle lp3
          endif

          !Search along both fronts to add bnd nodes to the list (up to max iterations)
          m=0 !init
          lp1: do 
            m=m+1
            if(m>max_iter.or.m>icount) exit lp1 !exhausted
            nd=nd_search(m)
            do j=1,nne(nd)
              ie2=indel(j,nd)
              do jj=1,i34(ie2)
                nd2=elnode(jj,ie2)
                if(icolor2(nd2)==0.and.isbnd(nd2)==-1) then
                  icolor2(nd2)=1
                  icount=icount+1
                  if(icount>np) then
                    print*, 'Too many search nodes(2):',ie,icount,nd_search(1:icount-1)
                    stop
                  endif
                  nd_search(icount)=nd2
                  cycle !lp1
                endif
              enddo !jj
            enddo !j
          end do lp1 !m: fronts
          !Debug
          write(99,*)'List of land bnd nodes for searching:',ie,nd_search(1:icount)

          perc(i)=flowmax/av_flow(i)
          flow_left=av_flow(i)-flowmax !amount to spread
          lp2: do m=1,icount
            nd=nd_search(m)
            do j=1,nne(nd)
              ie2=indel(j,nd)
              if(icolor(ie2)==0) then !new candidate cell
                av_h=sum(dp(elnode(1:i34(ie2),ie2)))/i34(ie2)
                if(av_h<=0.d0) cycle

                icolor(ie2)=1
                nsource2=nsource2+1
                if(nsource2>ne) stop 'Overflow'
                ie_src2(nsource2,1)=ie2
                ie_src2(nsource2,2)=ie !origin/parent source cell

                vol=area(ie2)*av_h !>0
                if(av_h>100.d0) then
                  s_rat=s_rat_deep
                else
                  s_rat=s_rat_shallow
                endif
                !Max flow allowed for S
                flowmax=vol/dt*(s_rat-1) !>0
                if(flowmax>=flow_left) then !done; append the new sources to the end of list
                  !% of flow received in this cell
                  perc(nsource2)=flow_left/av_flow(i)
                  exit lp2
                else !continue search
                  perc(nsource2)=flowmax/av_flow(i)
                  flow_left=flow_left-flowmax
                endif
              endif !icolor
            enddo !j=1,nne(nd)
          end do lp2 !m=1,icount
        endif !flowmax<av_flow
      end do lp3 !i=1,nsource
      print*, nsource2-nsource,' new sources added'

      !Check total %
      tot_perc(:)=0
      do i=nsource+1,nsource2
        ie0=ie_src2(i,2)
        ie=ie_src2(i,1)
        tot_perc(ie0)=tot_perc(ie0)+perc(i)
      enddo !i
      do i=1,nsource
        ie=ie_src(i)
        tmp=tot_perc(ie)+perc(i)
        write(98,*)'Total %=',i,ie,tmp
        if(abs(tmp-1.d0)>1.d-5) then
          print*, '% not right; check fort.98:',i,ie,tmp
          stop 
        endif
      enddo !i

      !Outputs
      write(15,*)nsource2
      do i=1,nsource2
        write(15,*)ie_src2(i,1)
      enddo !i

      write(15,*)
      write(15,*)nsink
      do i=1,nsink
        write(15,*)ie_sink(i)
      enddo !i
      close(15)

      allocate(flow2(nsource2))
      nline=0
      rewind(12)
      do
        read(12,*,end=98)time,flow(1:nsource)
        flow2(1:nsource)=flow(1:nsource)*perc(i)
        do i=nsource+1,nsource2
          ie0=ie_src2(i,2) !origin cell
          indx=0
          do j=1,nsource
            if(ie_src(j)==ie0) then
              indx=j
              exit
            endif
          enddo !j
          if(indx==0) then
            print*, 'Not found:',ie0
            stop 
          endif
          flow2(i)=flow(indx)*perc(i)
        enddo !i
        write(17,'(e22.12,1000000(1x,e14.6))')time,flow2(1:nsource2)
        nline=nline+1
      end do
98    close(12)

      print*, 'Finished!'

      stop
      end

      function signa(x1,x2,x3,y1,y2,y3)
!...  Compute signed area formed by pts 1,2,3
      implicit real*8(a-h,o-z)

      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2

      return
      end

!===============================================================================
!     Cross-product of two vectors: (x1,y1,z1) x (x2,y2,z2) = (x3,y3,z3)
!===============================================================================
      subroutine cross_product(x1,y1,z1,x2,y2,z2,x3,y3,z3)
      implicit real*8(a-h,o-z)
      real(8),intent(in) :: x1,y1,z1,x2,y2,z2
      real(8),intent(out) :: x3,y3,z3

      x3=y1*z2-y2*z1
      y3=x2*z1-x1*z2
      z3=x1*y2-x2*y1

      end subroutine cross_product

