!     Generate bottom_friction.gr3
!     Inputs: hgrid.gr3 (lon/lat), and vgrid.in, the input file
!     Outputs: bottom_friction.gr3
!     ifort -O2 -mcmodel=medium -CB -Bstatic -o gen_bottom_fric gen_bottom_fric.f90
      implicit real*8(a-h,o-z)
      allocatable :: x(:),y(:),dp(:),area(:),vso(:),r_bfric(:),nlayers(:)
      integer, allocatable :: i34(:),elnode(:,:),i_rain(:),isource(:),iDB(:)
      real(8) :: slope(4)

      print*, '--------------Input two hgrid depths to distinguish river, land, and the transition zone------------:'
      read(*,*) depth1,depth2
      print*,'<<<<<river: lower than ',depth1,' m;'
      print*,'<<<<<transition zone: ',depth1,'~',depth2,' m'
      print*,'<<<<<watershed: higher than ',depth2,' m;'

      print*, '-----------Input bottom friction in the river and on the land respectively:--------------'
      read(*,*) r_bfric_river,r_bfric_land
      print*,'<<<<<r_bfric_river',r_bfric_river,'; r_bfric_land: ',r_bfric_land

      print*, '-----------Force small friction on multilayer nodes:--------------'
      read(*,*) i_multi_layer
      print*,'<<<<<i_multi_layer',i_multi_layer
      if (i_multi_layer<0.and.i_multi_layer>1) then
        print *, 'wrong i_multi_layer', i_multi_layer
      endif

      open(14,file='hgrid.gr3')
      read(14,*)
      read(14,*)ne,np
      allocate(x(np),y(np),dp(np),r_bfric(np),nlayers(np),iDB(np),i_rain(np),isource(ne),i34(ne),elnode(4,ne),area(ne),vso(ne))
      do i=1,np
        read(14,*)j,x(i),y(i),dp(i)
      enddo !i
      do i=1,ne
        read(14,*)j,i34(i),elnode(1:i34(i),i)
      enddo
      close(14)

      open(15,file='vgrid.in')
      read(15,*) ivcor
      read(15,*) nvrt
      if(ivcor==1) then
        do i=1,np
          read(15,*) j,kbp
          nlayers(i)=nvrt-kbp
        enddo
      elseif (ivcor==2) then
        nlayers=nvrt
      endif
      close(15)


      open(8,file='bottom_friction.gr3',status='replace')
      write(8,*); write(8,*)ne,np

      !init
      r_bfric = r_bfric_river

      !based on bathymetry
      do i=1,np
        r_bfric(i)=r_bfric_river+(r_bfric_land-r_bfric_river)*(depth1-dp(i))/(depth1-depth2)
        r_bfric(i)=max(r_bfric_river,min(r_bfric_land,r_bfric(i)))
      enddo !i

      !multi-layers use small r_bfric
      if (i_multi_layer==1) then
        do i=1,ne
          if (maxval(nlayers(elnode(1:i34(i),i))) >1 ) then
            do j=1,i34(i)
              n1=elnode(j,i)
              r_bfric(n1) = r_bfric_river
            enddo
          endif
        enddo
      endif

      !write
      do i=1,np
        write(8,'(i8,3(1x,f15.6))')i,x(i),y(i),r_bfric(i)
      enddo 
      do i=1,ne
        write(8,*)i,i34(i),elnode(1:i34(i),i)
      enddo
      close(8)

      stop
      end
