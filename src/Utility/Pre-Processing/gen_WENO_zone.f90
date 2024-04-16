!     Generate tvd.prop for hybrid WENO/ELM and cross-scale applications, based on depth and remove
!     'isolated' WENO elements in shallow water, etc. The goal is to
!     minimize the direct 'contact' between WENO and ELM cells, which can cause dispersion.
!     Inputs: 
!             (1) hgrid.gr3 (in any projection or lon/lat; b.c. part not needed)
!             (2) gen_tvd_WENO.in:
!                   hmin_estu: cut-off depth for estuaries (e.g. same as h_tvd)
!                   nreg: # of regions (nearshore + estuaries; >=1)
!                   list of region names: 1st region is coastal region (e.g. using 20m isobath); the rest
!                        are optional estuaries. Inside coastal region,
!                        upwind is used except in estuaries. Make sure estuaries
!                        overlap with nearshore (so no gap in TVD/WENO zone). Outside coastal region,
!                        upwind/WENO is controlled by h_tvd only.
!     Output: tvd.prop.0 (may be further edited, e.g. look for very skew elem nearby,
!              and connectivity)

!     ifort -O2 -o gen_tvd_WENO ../UtilLib/schism_geometry.f90 ../UtilLib/pt_in_poly_test.f90 gen_tvd_WENO.f90

      use schism_geometry_mod
      use pt_in_poly_test
      implicit real*8(a-h,o-z)
      real*8, parameter :: small1=1.d-4
      real*8, parameter :: ray_angle=3.13192d0

      character(len=100) :: file_reg(200)
      integer :: nwild(3)
      integer, allocatable :: i34(:),elnode(:,:),nne(:),indel(:,:),nnp(:), &
     &indnd(:,:),isbnd(:),itvd(:),itvd0(:),iest(:),iest_e(:),inear(:),inear_e(:)
      integer, allocatable :: ic3(:,:),elside(:,:),isdel(:,:),isidenode(:,:),nvertices(:)
      real*8, allocatable :: xnd(:),ynd(:),dp(:),area(:),xpoly(:),ypoly(:),xcj(:),ycj(:)

      type :: poly_vertex
        real*8,pointer :: xy(:,:)=>null()
      end type
      type(poly_vertex),allocatable :: poly(:)

      !Read in control input
      open(9,file='gen_tvd_WENO.in',status='old')
      read(9,*) hmin_estu !cut-off depth for estuaries (e.g. h_tvd)
      read(9,*) nreg
      if(nreg<1.or.nreg>200) stop 'nreg<1 or nreg>200'
      do i=1,nreg
        read(9,*)file_reg(i)
      enddo !i
      close(9)

!      print*, 'Input # of estuary regions (can be 0):'
!      read*, nestuaries
!      if(nestuaries<0) stop 'nestuaries<0'
!      print*, 'Input cut-off depth (meters) for estuaries (e.g. h_tvd):'
!      read*, hmin_estu

      !Read in nearshore (1st file) and (optional) estuary regions
      nestuaries=nreg-1 !>=0
      allocate(nvertices(nreg),poly(nreg)) !first is the nearshore.rgn
!      open(10,file='nearshore.rgn',status='old')
!      read(10,*); read(10,*)npoly
!      if(npoly/=1) stop 'can only have 1 poly in .rgn'
!      read(10,*)nvertices(1)
!      allocate(poly(1)%xy(nvertices(1)+1,2)) !extra pt to close the region
!      poly(1)%xy(nvertices(1)+1,:)=0. !init
!      do i=1,nvertices(1)
!        read(10,*)poly(1)%xy(i,1:2) !xpoly(i,1),ypoly(i,1)
!      enddo !i
!      close(10)
!      
!      !Check if the region is closed
!      rl=(poly(1)%xy(1,1)-poly(1)%xy(nvertices(1),1))**2+(poly(1)%xy(1,2)-poly(1)%xy(nvertices(1),2))**2
!      if(rl>1.d-14) then
!        nvertices(1)=nvertices(1)+1
!        poly(1)%xy(nvertices(1),:)=poly(1)%xy(1,:)
!      endif

      do irgn=1,nreg
        open(10,file=trim(adjustl(file_reg(irgn))),status='old')
        read(10,*); read(10,*)npoly
        if(npoly/=1) then
          print*, 'Can only have 1 poly in .rgn:',irgn
          stop 
        endif
        read(10,*)nvertices(irgn)
        allocate(poly(irgn)%xy(nvertices(irgn)+1,2)) !extra pt to close the region
        poly(irgn)%xy(nvertices(irgn)+1,:)=0.
        do i=1,nvertices(irgn)
          read(10,*)poly(irgn)%xy(i,1:2) 
        enddo !i
        close(10)
      
        !Check if the region is closed
        rl=(poly(irgn)%xy(1,1)-poly(irgn)%xy(nvertices(irgn),1))**2+ &
         &+(poly(irgn)%xy(1,2)-poly(irgn)%xy(nvertices(irgn),2))**2
        if(rl>1.d-14) then
          nvertices(irgn)=nvertices(irgn)+1
          poly(irgn)%xy(nvertices(irgn),:)=poly(irgn)%xy(1,:)
        endif
      enddo !irgn

      open(14,file='hgrid.gr3',status='old')
      read(14,*); read(14,*) ne,np
      allocate(xnd(np),ynd(np),dp(np),area(ne),i34(ne),elnode(4,ne),nne(np), &
     &isbnd(np),itvd(ne),itvd0(ne),iest(np),iest_e(ne),inear(np),inear_e(ne))

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

      call compute_nside(np,ne,i34,elnode,ns)
!     Allocate side-related arrays
      allocate(ic3(4,ne),elside(4,ne),isdel(2,ns),isidenode(2,ns),xcj(ns),ycj(ns))
!     Then compute the rest of side related arrays with additional
!     inputs (xnd,ynd) (x,y coordinates of each node)
      call schism_geometry_double(np,ne,ns,xnd,ynd,i34,elnode,ic3,elside,isdel,isidenode,xcj,ycj)
      
!     Mark nearshore & estuary nodes
      iest=0 !init as outside any estuary of interest
      do i=1,np
        call pt_in_poly_ray_method_double(nvertices(1),small1,ray_angle, &
     &poly(1)%xy(1:nvertices(1),1),poly(1)%xy(1:nvertices(1),2),xnd(i),ynd(i),in_out,inters,npoly)
        if(in_out==1) then 
          inear(i)=1 !nearshore
        else
          inear(i)=0 !offshore
        endif

        !Estuaries
        do irgn=1,nestuaries
          call pt_in_poly_ray_method_double(nvertices(irgn+1),small1,ray_angle, &
     &poly(irgn+1)%xy(1:nvertices(irgn+1),1),poly(irgn+1)%xy(1:nvertices(irgn+1),2), &
     &xnd(i),ynd(i),in_out,inters,npoly)
          
          if(in_out==1) iest(i)=1 !inside any estuaries
        enddo !irgn
      enddo !i=1,np

      do i=1,ne
        inear_e(i)=maxval(inear(elnode(1:i34(i),i)))
        iest_e(i)=minval(iest(elnode(1:i34(i),i)))
      enddo !i

!     Set TVD flag
      itvd(:)=0 !init
      do i=1,ne
        hmin2=minval(dp(elnode(1:i34(i),i)))
        if((inear_e(i)==0.or.iest_e(i)/=0).and.hmin2>=hmin_estu) itvd(i)=1
      enddo !i

!     Augment upwind zone by 1 extra layer for offshore
      itvd0=itvd
      do i=1,ne
        if(inear_e(i)==0.and.itvd0(i)==0) then
          do j=1,i34(i)
            nd=elnode(j,i)
            do m=1,nne(nd)
              ie=indel(m,nd)
              itvd(ie)=0
            enddo !m
          enddo !j
        endif !ifl
      enddo !i

!     Remove 'isolated' WENO elem 
      loop1: do
        itouched=0 !# of elem changed in this iteration
        itvd0=itvd
        loop2: do i=1,ne
          if(itvd0(i)==0) cycle

          !WENO elem 
          icount=0
          do j=1,i34(i)
            ie=ic3(j,i)
            if(ie/=0) then 
              if(itvd0(ie)/=0) icount=icount+1
            endif !ie
          enddo !j 

          if(icount<=1) then
            itvd(i)=0
            itouched=itouched+1
          endif
        end do loop2 !i=1,ne

        print*, '# of elem flipped=',itouched
        if(itouched==0) exit loop1
      end do loop1

!     Output
      open(12,file='tvd.prop.0',status='replace')
      do i=1,ne
        write(12,*)i,itvd(i)
      enddo !i
      write(12,*)'hmin=',hmin_estu
      close(12)

      ! deallocate pointer arrays
      do m=1,size(poly)
        if (associated(poly(m)%xy)) deallocate(poly(m)%xy)
        nullify(poly(m)%xy)
      enddo
      deallocate(poly)

      stop
      end

      function signa(x1,x2,x3,y1,y2,y3)
!...  Compute signed area formed by pts 1,2,3
      implicit real*8(a-h,o-z)

      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2

      return
      end

