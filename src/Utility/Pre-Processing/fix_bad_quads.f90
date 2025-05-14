!     Fix bad quality quads based on angle ratio criterion (like
!     xmgredit5).
!     Inputs: screen; hgrid.gr3 (projection or lon/lat; quads can be
!     concave)
!     Output: hgrid.gr3.new (node order not changed); split_loc.bp (location of quads that are
!     split)

!     ifx -O2 -CB -g -traceback -o fix_bad_quads fix_bad_quads.f90

      implicit real*8(a-h,o-z)
      integer :: nwild(3),nwild2(4)
      real*8 :: skew(2,2)
      integer, allocatable :: i34(:),elnode(:,:),i34_new(:),elnode_new(:,:),isplit(:)
      real*8, allocatable :: xnd(:),ynd(:),dp(:),area(:),quad_loc(:,:), &
     &xctr(:),yctr(:)

      pi=acos(-1.d0)
      !Quality is measured by the ratio of smallest to largest interior
      !angle of quads
      print*, 'Input threshold for angle ratio (e.g. 0.5):'
      read*, rat_angle

      print*, 'Input threshold for quad aspect ratio (e.g. 5):'
      read*, aspect_rat

      open(14,file='hgrid.gr3',status='old')
      read(14,*)
      read(14,*) ne,np
      allocate(xnd(np),ynd(np),dp(np),i34(ne),elnode(4,ne),isplit(ne), &
     &i34_new(2*ne),elnode_new(4,2*ne),quad_loc(2,ne),xctr(ne),yctr(ne))

      do i=1,np
        read(14,*) j,xnd(i),ynd(i),dp(i)
      enddo !i

      !Init isplit: local index (=1 or 2) where splitting occurs; quad is split into tri
      !(isplit,isplit+1,isplit+2)
      isplit(:)=0 
      do i=1,ne
        read(14,*) j,i34(i),elnode(1:i34(i),i)

        n1=elnode(1,i)
        n2=elnode(2,i)
        n3=elnode(3,i)
        area2=signa(xnd(n1),xnd(n2),xnd(n3),ynd(n1),ynd(n2),ynd(n3))
        !Quads can be concave
        if(i34(i)==3) then
          if(area2<=0) then
            write(*,*)'Negative area at elem:',i,i34(i)
            stop
          endif
        else !quad
          n4=elnode(4,i)
          tmp=signa(xnd(n1),xnd(n3),xnd(n4),ynd(n1),ynd(n3),ynd(n4))
          if(area2<=0.or.tmp<=0) isplit(i)=2
!            write(*,*)'Negative area at elem:',i
!            stop
!          endif
!          area(i)=area(i)+tmp
        endif

        xctr(i)=sum(xnd(elnode(1:i34(i),i)))/i34(i)
        yctr(i)=sum(ynd(elnode(1:i34(i),i)))/i34(i)
      enddo !i=1,ne      
      i34_new(1:ne)=i34
      elnode_new(:,1:ne)=elnode(:,1:ne)

      ne_extra=0
      do i=1,ne
        if(i34(i)/=4) cycle

        !Quads: calc min/max side for aspect ratio later
        sdmin=huge(1.d0); sdmax=-armin
        do j=1,4
          n1=elnode(j,i) 
          j1=j+1
          if(j1>4) j1=j1-4
          n2=elnode(j1,i)
          distl=sqrt((xnd(n1)-xnd(n2))**2+(ynd(n1)-ynd(n2))**2)
          sdmin=min(sdmin,distl)
          sdmax=max(sdmax,distl)
        enddo !j
        if(min(sdmax,sdmin)<=0.) then
          write(*,*)'Wrong aspect ratio::',i,sdmax,sdmin
          stop
        endif

        !Quads; interior angles
        ang_max=-1; ang_min=pi*1.1
        do j=1,4
          n1=elnode(j,i) !angle centered @ n1
          j1=j+1; j3=j+3
          if(j1>4) j1=j1-4
          if(j3>4) j3=j3-4
          n2=elnode(j1,i)
          n3=elnode(j3,i)

          ar2=signa(xnd(n1),xnd(n2),xnd(n3),ynd(n1),ynd(n2),ynd(n3))
          if(ar2<=0) then !concave quad; split from this node
            isplit(i)=j; exit
          endif

          !Convex @ n1; interior angle
          dot=(xnd(n2)-xnd(n1))*(xnd(n3)-xnd(n1))+(ynd(n2)-ynd(n1))*(ynd(n3)-ynd(n1))
          rl1=sqrt((xnd(n2)-xnd(n1))**2+(ynd(n2)-ynd(n1))**2)
          rl3=sqrt((xnd(n3)-xnd(n1))**2+(ynd(n3)-ynd(n1))**2)
          if(rl1==0.or.rl3==0) stop '0 side'
          costmp=dot/rl1/rl3
          if(abs(costmp)>1) then
            print*, 'Impossible cosine:',i,j,costmp
            stop
          endif
          ang=acos(costmp)
          if(ang>ang_max) ang_max=ang
          if(ang<ang_min) ang_min=ang
        enddo !j=1,4

        if(isplit(i)==0) then !calc angle ratio for convex quad
          if(ang_max==0.or.ang_max<ang_min) stop 'failed to find max/min'
!'
          if(ang_min/ang_max<rat_angle.or.sdmax/sdmin>aspect_rat) then
            !Decide which node (1 or 2) to split from, based on best aspect ratios
            !calc 2 A.R.'s
            do j=1,2 !2 options of splitting
              if(j==1) then
                nwild2(1:4)=(/1,2,3,4/)
              else
                nwild2(1:4)=(/2,3,4,1/)
              endif

              do k=1,2 !pair of tri's
                if(k==1) then
                  nwild(1:3)=(/1,2,3/)
                else
                  nwild(1:3)=(/1,3,4/)
                endif
                n1=elnode(nwild2(nwild(1)),i)
                n2=elnode(nwild2(nwild(2)),i)
                n3=elnode(nwild2(nwild(3)),i)
                ar3=signa(xnd(n1),xnd(n2),xnd(n3),ynd(n1),ynd(n2),ynd(n3))
                if(ar3<=0) then
!                  print*, 'ar3<=0:',i,j,k,n1,n2,n3
!                  stop
                  skew(k,j)=-1. !flag
                else
                  rl1=sqrt((xnd(n1)-xnd(n2))**2+(ynd(n1)-ynd(n2))**2)
                  rl2=sqrt((xnd(n1)-xnd(n3))**2+(ynd(n1)-ynd(n3))**2)
                  rl3=sqrt((xnd(n2)-xnd(n3))**2+(ynd(n2)-ynd(n3))**2)
                  rlmax=max(rl1,rl2,rl3)
                  skew(k,j)=rlmax/sqrt(ar3/pi) !skewness
                endif
              enddo !k
            enddo !j=1,2
            
            !Change here if you want to use different criteria
            if(minval(skew(:,1))<0.and.minval(skew(:,2))<0) then
              print*, 'Cannot find a split point at elem:',i
              stop
            else if(minval(skew(:,1))<0) then
              isplit(i)=2
            else if(minval(skew(:,2))<0) then
              isplit(i)=1
            else
              arat1=sum(skew(:,1))
              arat2=sum(skew(:,2))
              if(arat1<arat2) then
                isplit(i)=1 !split at node 1
              else
                isplit(i)=2 !split at node 2
              endif
            endif
          endif !ang_min/
        endif !isplit==0

        !Do splitting of quad
        if(isplit(i)/=0) then !add extra tri to the end of conn table
          ne_extra=ne_extra+1
          if(ne_extra>ne) stop 'overflow(1)'
          i34_new(i)=3; i34_new(ne+ne_extra)=3 
          i1=isplit(i); i2=isplit(i)+1; i3=isplit(i)+2; i4=isplit(i)+3
          if(i2>4) i2=i2-4
          if(i3>4) i3=i3-4
          if(i4>4) i4=i4-4
  
          !Check
          n1=elnode(i1,i); n2=elnode(i2,i); n3=elnode(i3,i); n4=elnode(i4,i)
          ar3=signa(xnd(n1),xnd(n2),xnd(n3),ynd(n1),ynd(n2),ynd(n3))
          ar4=signa(xnd(n1),xnd(n3),xnd(n4),ynd(n1),ynd(n3),ynd(n4))
          if(ar3<=0.or.ar4<=0) then
            print*, 'Final orientation wrong:',i,ar3,ar4,n1,n2,n3,n4
            stop
          endif

          nwild(1:3)=(/i1,i2,i3/)
          elnode_new(1:3,i)=elnode(nwild(1:3),i)
          nwild(1:3)=(/i1,i3,i4/)
          elnode_new(1:3,ne+ne_extra)=elnode(nwild(1:3),i)

          quad_loc(1,ne_extra)=xctr(i)
          quad_loc(2,ne_extra)=yctr(i)
        endif !isplit/=0
      enddo !i=1,ne

      print*, ne_extra,' quads split; see split_loc.bp'
      ne=ne+ne_extra

!     Output
      open(13,file='split_loc.bp',status='replace')
      write(13,*); write(13,*)ne_extra
      do i=1,ne_extra
        write(13,'(i12,2(1x,e22.14),1x,f8.1)')i,quad_loc(1:2,i),0.
      enddo !i
      close(13)

      open(12,file='hgrid.gr3.new',status='replace')
      write(12,*)'Threshold of ang ratio=',rat_angle
      write(12,*)ne,np
      do i=1,np
        write(12,'(i12,2(1x,e22.14),1x,e16.5)')i,xnd(i),ynd(i),dp(i) 
      enddo !i
      nquads=0
      do i=1,ne
        write(12,*)i,i34_new(i),elnode_new(1:i34_new(i),i)
        if(i34_new(i)==4) nquads=nquads+1
      enddo !i
      close(12)
      print*, '# of quads remaining=',nquads

      stop
      end

      function signa(x1,x2,x3,y1,y2,y3)
!...  Compute signed area formed by pts 1,2,3
      implicit real*8(a-h,o-z)

      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2

      return
      end

