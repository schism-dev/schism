!     Combine 2 tvd.prop.[1,2] based a specified region (final flag is
!     from '2' inside region2.rgn). The 2 TVD files may be generated
!     from gen_tvd_WENO.f90.
!     Inputs: 
!             hgrid.gr3 (in any projection or lon/lat; b.c. part not needed)
!             tvd.prop.[1,2]
!             region2.rgn (gredit region format)
!     Output: tvd.prop

!     ifort -O2 -o combine_TVD UtilLib/pt_in_poly_test.f90 combine_TVD.f90

      use pt_in_poly_test
      implicit real*8(a-h,o-z)
      integer :: nwild(3)
      integer, allocatable :: i34(:),elnode(:,:),nne(:),indel(:,:),nnp(:), &
     &indnd(:,:),isbnd(:),itvd(:),itvd1(:),itvd2(:)
      real*8, allocatable :: xnd(:),ynd(:),dp(:),area(:),xpoly(:),ypoly(:)

      open(14,file='hgrid.gr3',status='old')
      read(14,*)
      read(14,*) ne,np
      allocate(xnd(np),ynd(np),dp(np),area(ne),i34(ne),elnode(4,ne),nne(np), &
     &isbnd(np),itvd(ne),itvd1(ne),itvd2(ne))

      do i=1,np
        read(14,*) j,xnd(i),ynd(i),dp(i)
      enddo !i

      do i=1,ne
        read(14,*) j,i34(i),elnode(1:i34(i),i)

        n1=elnode(1,i)
        n2=elnode(2,i)
        n3=elnode(3,i)
        area(i)=signa(xnd(n1),xnd(n2),xnd(n3),ynd(n1),ynd(n2),ynd(n3))
        if(area(i)<=0) then
          write(*,*)'Negative area at elem:',i
          stop
        endif
         
        if(i34(i)==4) then
          n4=elnode(4,i)
          tmp=signa(xnd(n1),xnd(n3),xnd(n4),ynd(n1),ynd(n3),ynd(n4))
          if(tmp<=0) then
            write(*,*)'Negative area at elem:',i
            stop
          endif
          area(i)=area(i)+tmp
        endif
      enddo !i=1,ne      
      close(14)

      !Read in .rgn
      open(10,file='region2.rgn',status='old')
      read(10,*); read(10,*)npoly
      if(npoly/=1) stop 'can only have 1 poly in .rgn'
      read(10,*)nvertices
      nvertices=nvertices+1 !repeat 1st vertex
      allocate(xpoly(nvertices),ypoly(nvertices))
      do i=1,nvertices-1
        read(10,*)xpoly(i),ypoly(i)
      enddo !i
      xpoly(nvertices)=xpoly(1)
      ypoly(nvertices)=ypoly(1)
      close(10)

!     Neighborhood
      nne=0
      do i=1,ne
        do j=1,i34(i)
          nd=elnode(j,i)
          nne(nd)=nne(nd)+1
        enddo
      enddo !i
      mnei=maxval(nne)

      allocate(indel(mnei,np))
      nne=0
!      nnp=0 !estimate 1st
      do i=1,ne
        do j=1,i34(i)
          nd=elnode(j,i)
          nne(nd)=nne(nd)+1
          if(nne(nd)>mnei) stop 'impossible'
          indel(nne(nd),nd)=i
!          nnp(nd)=nnp(nd)+i34(i)-1
        enddo
      enddo !i
      
!     Output
      open(8,file='tvd.prop.1',status='old')
      open(9,file='tvd.prop.2',status='old')
      do i=1,ne
        read(8,*)j,tmp1
        read(9,*)j,tmp2
        itvd1(i)=nint(tmp1)
        itvd2(i)=nint(tmp2)
      enddo !i

      open(12,file='tvd.prop',status='replace')
      small1=1.d-4
      ray_angle=3.13192d0
      itvd(:)=1 !init
      do i=1,ne
        !Beware lon jumps
        xctr=sum(xnd(elnode(1:i34(i),i)))/i34(i)
        yctr=sum(ynd(elnode(1:i34(i),i)))/i34(i)
        call pt_in_poly_ray_method_double(nvertices,small1,ray_angle,xpoly,ypoly,xctr,yctr,in_out,inters,npoly)
        if(in_out==1) then
          itvd(i)=itvd2(i)
        else
          itvd(i)=itvd1(i)
        endif !in_out
      enddo !i

      do i=1,ne
        write(12,*)i,itvd(i)
      enddo !i
      close(12)

      stop
      end

      function signa(x1,x2,x3,y1,y2,y3)
!...  Compute signed area formed by pts 1,2,3
      implicit real*8(a-h,o-z)

      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2

      return
      end

