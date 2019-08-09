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

!
!****************************************************************************************
!											*
!	Read binary format v5.0 (hybrid S-Z models) convert into netcdf (structured grid CF-1.0) format
!       readable into ArcMap. Only works for 2D scalar at the moment.
!
!       Inputs: from screen; binary must be in CPP projection (not lat/lon)
!       Outputs: 
!               *.nc (all stacks combined into 1 file; in lat/lon); fort.12 (non fatal); fort.11 (fatal; also screen)
!       History: 
!       WARNING: CF convention is very picky; (x,y) must use lat/lon (otherwise cannot overlay on maps)
!****************************************************************************************
!     Tsuanmi:
!     ifort -Bstatic -assume byterecl -O2 -o binary_2_ncstruc binary_2_ncstruc.f90 ../Grid_Scripts/stripsearch_unstr.f90  -Vaxlib -I/share/apps/netcdf/include/ -L/share/apps/netcdf/lib/ -lnetcdf

!     INDIAN (not Pacific):
!     module load compiler64/pgi-11.5
!     pgf90 -Bstatic -O2 -o binary_2_ncstruc binary_2_ncstruc.f90 ~/Scripts/stripsearch_unstr.f90 -I/usr/vimssw/pgi-11.5-linux86-64/netcdf-3.6.2/include -L/usr/vimssw/pgi-11.5-linux86-64/netcdf-3.6.2/lib -lnetcdf

      program binary_2_ncstruc
      include 'netcdf.inc'

      integer,parameter :: nbyte=4
      character*30 file63
      character*12 it_char
      character*48 start_time,version,variable_nm,variable_dim
      character*48 data_format

!     Double precision for search
      real(8), allocatable :: x(:),y(:),arco(:,:,:),xybin(:),xlon(:),ylat(:),area(:), &
                              &timeout(:),xproj(:,:),yproj(:,:)
      real(8) :: xy,suma,binwid,signa,dx,dy,xmin,xmax,ymin,ymax,xlmin,xlmax,ylmin,ylmax,dxl,dyl,rlambda0,phi0
      integer,allocatable :: elnode(:,:)
      integer,allocatable :: elside(:,:)
      integer,allocatable :: isdel(:,:)
      allocatable :: sigma(:),cs(:),ztot(:),dp(:),kbp00(:),kfp(:)
      allocatable :: out(:,:,:),icum(:,:,:),eta2(:),ztmp(:,:)
      allocatable :: xcj(:),ycj(:),dps(:),kbs(:),xctr(:),yctr(:),dpe(:),kbe(:)
      allocatable :: zs(:,:),ze(:,:),idry(:),outnc(:,:),outs(:,:,:),oute(:,:,:)
      allocatable :: ic3(:,:),isidenode(:,:)
      allocatable :: ne_bin(:),ie_bin(:,:),iparen(:,:),nne(:),indel(:,:)
      integer :: iwest_south(4,2)

!     netcdf variables
      character(len=50) :: fname
      integer :: stat,time_dims(1),ele_dims(2),xy_dims(1),var2d_dims(2),var3d_dims(3), &
                &data_start_2d(2),data_start_3d(3),data_count_2d(2),data_count_3d(3),z_dims(1), &
                &intval(1)
      real :: realval(1)
      
      print*, 'Input file to read from (without *_; e.g. salt.63):'
      read(*,'(a30)')file63
      file63=adjustl(file63)
      lfile63=len_trim(file63)
      
      print*, 'Input # of files (stacks)  to read:'
      read(*,*) ndays

!      print*, 'Input (xmin,ymin) in CPP projection:'
!      read(*,*) xmin,ymin

      print*, 'Input (xmin,ymin) in lon/lat:'
      read(*,*) xlmin,ylmin

!      print*, 'Input (xmax,ymax) in projection:'
!      read(*,*) xmax,ymax

      print*, 'Input (xmax,ymax) in lon/lat:'
      read(*,*) xlmax,ylmax

      print*, 'Input CPP projection center (degrees):'
      read(*,*) rlambda0,phi0

      print*, 'Input output dx,dy in meters:'
      read(*,*) dx,dy

      print*, 'Input search direction (1: x; 2: y):'
      read(*,*) is_xy

      print*, 'Input # of bins used in search (suggest 50000):'
      read(*,*) nbin

      print*, 'Input max. # of elem. per bin (suggest 30000):'
      read(*,*) mne_bin

      allocate(xybin(nbin+1),ne_bin(nbin),ie_bin(nbin,mne_bin),stat=istat)
      if(istat/=0) stop 'Falied to allocate (-1)'

      call cpp(1,xlmin,ylmin,xmin,ymin,rlambda0,phi0)
      call cpp(1,xlmax,ylmax,xmax,ymax,rlambda0,phi0)

      nx=(xmax-xmin)/dx+1
      ny=(ymax-ymin)/dy+1
      print*, 'xmin,ymin=',xmin,ymin
      print*, 'xmax,ymax=',xmax,ymax
      print*, 'Output nx,ny=',nx,ny

!     Compute x,y
      allocate(xproj(nx,ny),yproj(nx,ny),xlon(nx),ylat(ny),iparen(nx,ny),arco(3,nx,ny),outnc(nx,ny),stat=istat)
      if(istat/=0) stop 'Falied to allocate (-2)'
!      do i=1,nx
!        xproj(i,:)=xmin+dx*(i-1) !in meters
!      enddo !i
!      do j=1,ny
!        yproj(:,j)=ymin+dy*(j-1)
!      enddo !j

!     lat/lon for output nc
      dxl=(xlmax-xlmin)/(nx-1)
      dyl=(ylmax-ylmin)/(ny-1)

      print*, 'dx,dy in degrees=',dxl,dyl
      do i=1,nx
        xlon(i)=xlmin+dxl*(i-1)
!        write(95,*)i,xlon(i)
      enddo !i
      do j=1,ny
        ylat(j)=ylmin+dyl*(j-1)
!        write(95,*)j,ylat(j)
      enddo !i

!     CPP
      do i=1,nx
        do j=1,ny
          call cpp(1,xlon(i),ylat(j),xproj(i,j),yproj(i,j),rlambda0,phi0)
          write(99,*)i,j,real(xproj(i,j)),real(yproj(i,j))
        enddo !i
      enddo !i

!...  Header
!...
      open(63,file='1_'//file63,status='old',access='direct',recl=nbyte)
      irec=0
      do m=1,48/nbyte
        read(63,rec=irec+m) data_format(nbyte*(m-1)+1:nbyte*m)
      enddo
      if(data_format.ne.'DataFormat v5.0') then
        print*, 'This code reads only v5.0:  ',data_format
        stop
      endif
      irec=irec+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) version(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) start_time(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) variable_nm(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) variable_dim(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte

      !Overwrite start_time
      start_time='1900-1-1'

      write(*,'(a48)')data_format
      write(*,'(a48)')version
      write(*,'(a48)')start_time
      write(*,'(a48)')variable_nm
      write(*,'(a48)')variable_dim

      read(63,rec=irec+1) nrec
      read(63,rec=irec+2) dtout
      read(63,rec=irec+3) nspool
      read(63,rec=irec+4) ivs
      read(63,rec=irec+5) i23d
      irec=irec+5

      print*, 'i23d=',i23d,' nrec= ',nrec
      if(i23d/=2) stop '3D not working at the moment'

!     Vertical grid
      read(63,rec=irec+1) nvrt
      read(63,rec=irec+2) kz
      read(63,rec=irec+3) h0
      read(63,rec=irec+4) h_s
      read(63,rec=irec+5) h_c
      read(63,rec=irec+6) theta_b
      read(63,rec=irec+7) theta_f
      irec=irec+7
      allocate(ztot(nvrt),sigma(nvrt),cs(nvrt),timeout(ndays*nrec),stat=istat)
      if(istat/=0) stop 'Falied to allocate (1)'

      do k=1,kz-1
        read(63,rec=irec+k) ztot(k)
      enddo
      do k=kz,nvrt
        kin=k-kz+1
        read(63,rec=irec+k) sigma(kin)
        cs(kin)=(1-theta_b)*sinh(theta_f*sigma(kin))/sinh(theta_f)+ &
     &theta_b*(tanh(theta_f*(sigma(kin)+0.5))-tanh(theta_f*0.5))/2/tanh(theta_f*0.5)
      enddo
      irec=irec+nvrt

!     Horizontal grid
      read(63,rec=irec+1) np !could be ns,ne also
      read(63,rec=irec+2) ne
      irec=irec+2
      allocate(x(np),y(np),dp(np),kbp00(np),kfp(np),&
     &elnode(3,ne),out(np,nvrt,2),idry(np), &
     !&icum(np,nvrt,2),eta2(np),stat=istat)
     &eta2(np),xctr(ne),yctr(ne),dpe(ne),kbe(ne),ztmp(nvrt,np),ze(nvrt,ne), &
     &area(ne),nne(np),stat=istat)
      if(istat/=0) stop 'Falied to allocate (2)'

      do m=1,np
        read(63,rec=irec+1)tmp 
        x(m)=tmp !real*8
        read(63,rec=irec+2)tmp
        y(m)=tmp
        read(63,rec=irec+3)dp(m)
        read(63,rec=irec+4)kbp00(m)
        irec=irec+4
      enddo !m=1,np

!     Additional for non-standard outputs
!      if(i23d==4) then !3D side @ whole levels
!        read(63,rec=irec+1) ns
!        irec=irec+1
!        allocate(xcj(ns),ycj(ns),dps(ns),kbs(ns),zs(nvrt,ns),outs(2,nvrt,ns),stat=istat)
!        if(istat/=0) stop 'Falied to allocate (5)'
!        do m=1,ns
!          read(63,rec=irec+1)xcj(m)
!          read(63,rec=irec+2)ycj(m)
!          read(63,rec=irec+3)dps(m)
!          read(63,rec=irec+4)kbs(m)
!          irec=irec+4
!        enddo !m
!      else if (i23d>4) then !3D elem. @ whole/half
!        read(63,rec=irec+1) netmp
!        irec=irec+1
!        do m=1,ne
!          read(63,rec=irec+1)xctr(m)
!          read(63,rec=irec+2)yctr(m)
!          read(63,rec=irec+3)dpe(m)
!          read(63,rec=irec+4)kbe(m)
!          irec=irec+4
!        enddo !m
!      endif !i23d

      do m=1,ne
        read(63,rec=irec+1)i34
        if(i34/=3) stop 'No quads pls'
        irec=irec+1
        do mm=1,i34
          read(63,rec=irec+1)elnode(mm,m)
          irec=irec+1
        enddo !mm

        !Area
        n1=elnode(1,m); n2=elnode(2,m); n3=elnode(3,m)
        area(m)=dabs(signa(x(n1),x(n2),x(n3),y(n1),y(n2),y(n3)))
        if(area(m)<=0) stop 'Zero area'
      enddo !m
      irec0=irec

!     Compute geometry
!      call compute_nside(np,ne,elnode,ns2)
!      if(.not.allocated(xcj)) allocate(xcj(ns2))
!      if(.not.allocated(ycj)) allocate(ycj(ns2))
!      allocate(ic3(ne,3),elside(3,ne),isdel(2,ns2),isidenode(2,ns2),stat=istat)
!      if(istat/=0) stop 'Allocation error: side(0)'
!      call schism_geometry(np,ne,ns2,x,y,elnode,ic3,elside,isdel,isidenode,xcj,ycj)

      print*, 'last element',(elnode(j,ne),j=1,3)

!...  Compute ball
      nne=0
      do i=1,ne
        do j=1,3
          nd=elnode(j,i)
          nne(nd)=nne(nd)+1
        enddo
      enddo
      mnei=maxval(nne)

      allocate(indel(mnei,np),stat=istat)
      if(istat/=0) stop 'Failed to alloc. indel'
      nne=0
      do i=1,ne
        do j=1,3
          nd=elnode(j,i)
          nne(nd)=nne(nd)+1
          if(nne(nd)>mnei) then
            write(*,*)'Too many neighbors',nd
            stop
          endif
          indel(nne(nd),nd)=i
        enddo
      enddo !i

!...  Compute relative record # for a node and level for 3D outputs
!...
!      icount=0
!      do i=1,np
!        do k=max0(1,kbp00(i)),nvrt
!          do m=1,ivs
!            icount=icount+1
!            icum(i,k,m)=icount
!          enddo !m
!        enddo !k
!      enddo !i=1,np

!...  Do strip sort
      call stripsearch_unstr(is_xy,nbin,mne_bin,ne,np,x,y,elnode, &
     &ne_bin,ie_bin,xybin,binwid)

!...  Find parent elements and weights
      iparen=0
      do i=1,nx
        do j=1,ny
          if(mod(i,1000)==0.and.mod(j,1000)==0) print*, 'finding parent for (i,j):',i,j

          !First try neighbors of previously found elem. to speed up
          iwest_south(1,1)=i-1; iwest_south(1,2)=j-1 !(i,j) pair
          iwest_south(2,1)=i-1; iwest_south(2,2)=j
          iwest_south(3,1)=i-1; iwest_south(3,2)=j+1
          iwest_south(4,1)=i; iwest_south(4,2)=j-1
          if(i>1.or.j>1) then !search all west and south neighbors
            loop2: do i2=1,4 !4 pairs
              if(iwest_south(i2,1)>0.and.iwest_south(i2,2)>0) then
                !if(j==1) then !use i-1
                !  i0=i-1; j0=j
                !else !j>1
                !  i0=i; j0=j-1
                !endif !j
                ie0=iparen(iwest_south(i2,1),iwest_south(i2,2)) !previous parent
                if(ie0>0) then
                  do m=1,3
                    nd=elnode(m,ie0)
                    do k=1,nne(nd)
                      ie=indel(k,nd)
                      do jj=1,3
                        j_1=jj+1
                        j_2=jj+2
                        if(j_1>3) j_1=j_1-3
                        if(j_2>3) j_2=j_2-3
                        n1=elnode(j_1,ie)
                        n2=elnode(j_2,ie)
                        arco(jj,i,j)=dabs(signa(x(n1),x(n2),xproj(i,j),y(n1),y(n2),yproj(i,j)))/area(ie)
                      enddo !jj
                      if(dabs(sum(arco(1:3,i,j))-1)<1.e-4) then
                        iparen(i,j)=ie
                        exit loop2
                      endif
                    enddo !k
                  enddo !m
                endif !ie0>0 
              endif !iwest_south(i2,1)
            end do loop2 !i2=1,4
          endif !i>1.or.j>1
   
          if(iparen(i,j)==0) then !full blown search
            if(is_xy==1) then
              xy=xproj(i,j)
            else
              xy=yproj(i,j)
            endif
     
!            if(xy<xy_min.or.xy>xy_max) return
            l=(xy-xybin(1))/binwid+1
            if(l<=0.or.l>nbin+1) then !outside
              write(12,*)'Out of bg grid:',i,j,xproj(i,j),yproj(i,j),l
              cycle
            endif

            l=min(nbin,l)
            if(xy==xybin(l)) then
              ibin1=max(l-1,1); ibin2=l
            else if(xy==xybin(l+1)) then
              ibin1=l; ibin2=min(l+1,nbin)
            else if(xy>xybin(l).and.xy<xybin(l+1)) then
              ibin1=l; ibin2=l
            else
              write(*,*)'Cannot find a bin (2):',i,j,xproj(i,j),yproj(i,j),xybin(l),xybin(l+1),binwid,l
              write(11,*)xybin
              stop
            endif
  
!          print*, 'ibin1,ibin2:',ibin1,ibin2
!          do l=ibin1,ibin2
!            do k=1,ne_bin(l)
!              ie=ie_bin(l,k)
!              write(95,*)l,k,ie
!            enddo !k
!          enddo !l

            loop1: do l=ibin1,ibin2
              do k=1,ne_bin(l)
                ie=ie_bin(l,k)
                do jj=1,3
                  j_1=jj+1
                  j_2=jj+2
                  if(j_1>3) j_1=j_1-3
                  if(j_2>3) j_2=j_2-3
                  n1=elnode(j_1,ie)
                  n2=elnode(j_2,ie)
                  arco(jj,i,j)=dabs(signa(x(n1),x(n2),xproj(i,j),y(n1),y(n2),yproj(i,j)))/area(ie)
                enddo !jj
                if(dabs(sum(arco(1:3,i,j))-1)<1.e-4) then
                  iparen(i,j)=ie
                  exit loop1
                endif
              enddo !k=1,ne_bin(l)
            end do loop1 !l
          endif !iparen(i,j)==0
        enddo !j=1,ny
      enddo !i=1,nx

!     Check parents
      do i=1,nx
        do j=1,ny
          if(iparen(i,j)==0) then
            write(12,*)'Failed to find a parent elem:',i,j,xproj(i,j),yproj(i,j)
!          else
!            write(55,*)'Parent elem=',iparen(i,j),i,j,xproj(i,j),yproj(i,j)
          endif
        enddo !j=1,ny
      enddo !i=1,nx

      print*,'done finding parents...'

!     Write nc header
!     enter define mode
      fname=file63(1:lfile63-3)//'.nc'
      iret = nf_create(trim(fname), NF_CLOBBER, ncid)
!     define dimensions
      iret = nf_def_dim(ncid, 'lon',nx, nx_dim)
      iret = nf_def_dim(ncid, 'lat',ny, ny_dim)
!      iret = nf_def_dim(ncid, 'time', NF_UNLIMITED, ntime_dim)
      iret = nf_def_dim(ncid, 'time', nrec*ndays, ntime_dim)

!     define variables
      xy_dims(1)=nx_dim
      iret=nf_def_var(ncid,'lon',NF_DOUBLE,1,xy_dims,ix_id)

!      var2d_dims(1)=nbnds_dim; var2d_dims(2)=nx_dim
!      iret=nf_def_var(ncid,'lon_bnds',NF_DOUBLE,nbnds,var2d_dims,ixbnd_id)

      xy_dims(1)=ny_dim
      iret=nf_def_var(ncid,'lat',NF_DOUBLE,1,xy_dims,iy_id)

      time_dims(1) = ntime_dim
      iret=nf_def_var(ncid,'time',NF_DOUBLE,1,time_dims,itime_id)

      var3d_dims(1)=nx_dim; var3d_dims(2)=ny_dim; var3d_dims(3)=ntime_dim
      iret=nf_def_var(ncid,'elev',NF_REAL,3,var3d_dims,ie_id)

!     Global attributes
      stat = nf_put_att_text(ncid, NF_GLOBAL, 'Conventions', 6, 'CF-1.0')

!     assign per-variable attributes
      iret=nf_put_att_text(ncid,ix_id,'standard_name',9,'longitude')
      iret=nf_put_att_text(ncid,ix_id,'long_name',9,'longitude')
      iret=nf_put_att_text(ncid,ix_id,'units',12,'degrees_east')
      iret=nf_put_att_text(ncid,ix_id,'axis',1,'X')
      iret=nf_put_att_text(ncid,ix_id,'original_units',12,'degrees_east')

      iret=nf_put_att_text(ncid,iy_id,'standard_name',8,'latitude')
      iret=nf_put_att_text(ncid,iy_id,'long_name',8,'latitude')
      iret=nf_put_att_text(ncid,iy_id,'units',13,'degrees_north')
      iret=nf_put_att_text(ncid,iy_id,'axis',1,'Y')
      iret=nf_put_att_text(ncid,iy_id,'original_units',13,'degrees_north')
!'

      iret=nf_put_att_text(ncid,itime_id,'standard_name',4,'time')
      iret=nf_put_att_text(ncid,itime_id,'long_name',4,'time')
      iret=nf_put_att_text(ncid,itime_id,'units',22,'seconds since 0000-1-1')
!'
      iret=nf_put_att_text(ncid,itime_id,'axis',1,'T')
      stat = nf_put_att_text(ncid,itime_id,'calendar',7,'365_day')
      stat = nf_put_att_text(ncid,itime_id,'original_units',22,'seconds since 0000-1-1')
!'

      iret=nf_put_att_text(ncid,ie_id,'standard_name',13,'Elevation (m)')
      stat = nf_put_att_text(ncid,ie_id,'long_name',19,'Water surface elev.')
      iret=nf_put_att_text(ncid,ie_id,'units',1,'m')
      realval(1)=-9999
      stat = nf_put_att_real(ncid,ie_id,'_FillValue',nf_float,1, realval)
      stat = nf_put_att_real(ncid,ie_id,'missing_value',nf_float,1,realval)

!     leave define mode
      iret = nf_enddef(ncid)

!     Write mode (header part only)
!     write time stamps 
      timeout=(/(i*dtout,i=1,ndays*nrec)/)

      !Scale x,y to fake as lon/lat
      iret=nf_put_vara_double(ncid,ix_id,1,nx,xlon(:))
      iret=nf_put_vara_double(ncid,iy_id,1,ny,ylat(:))
      iret=nf_put_vara_double(ncid,itime_id,1,ndays*nrec,timeout)
!     End nc header part

!...  Time iteration
!...
      it_tot=0
      out=-99 !init.
      ztmp=-99
      do iday=1,ndays
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      write(it_char,'(i12)')iday
      it_char=adjustl(it_char)  !place blanks at end
      it_len=len_trim(it_char)
      open(63,file=it_char(1:it_len)//'_'//file63,status='old',access='direct',recl=nbyte)
!'
      irec=irec0

      do it1=1,nrec
        read(63,rec=irec+1) time
        read(63,rec=irec+2) it
        irec=irec+2
        it_tot=it_tot+1
        time=it_tot*dtout

        print*, 'time in seconds =',time

        do i=1,np
          read(63,rec=irec+i) eta2(i)
        enddo !i
        irec=irec+np

        if(1==2) then
!---------------------------------------------------------
!       Compute z coordinates
        do i=1,np
          if(eta2(i)+dp(i)<h0) then !node dry
            idry(i)=1
          else !wet
            idry(i)=0
!           Compute z-coordinates
            do k=kz,nvrt
              kin=k-kz+1
              hmod2=min(dp(i),h_s)
              if(hmod2<=h_c) then
                ztmp(k,i)=sigma(kin)*(hmod2+eta2(i))+eta2(i)
              else if(eta2(i)<=-h_c-(hmod2-h_c)*theta_f/sinh(theta_f)) then
                write(*,*)'Pls choose a larger h_c (2):',eta2(i),h_c
                stop
              else
                ztmp(k,i)=eta2(i)*(1+sigma(kin))+h_c*sigma(kin)+(hmod2-h_c)*cs(kin)
              endif

!             Following to prevent underflow
              if(k==kz) ztmp(k,i)=-hmod2
              if(k==nvrt) ztmp(k,i)=eta2(i)
            enddo !k

            if(dp(i)<=h_s) then
              kbpl=kz
            else !z levels
!             Find bottom index
              kbpl=0
              do k=1,kz-1
                if(-dp(i)>=ztot(k).and.-dp(i)<ztot(k+1)) then
                  kbpl=k
                  exit
                endif
              enddo !k
              if(kbpl==0) then
                write(*,*)'Cannot find a bottom level:',dp(i)
                stop
              endif
              ztmp(kbpl,i)=-dp(i)
              do k=kbpl+1,kz-1
                ztmp(k,i)=ztot(k)
              enddo !k
            endif
            if(kbpl/=kbp00(i)) stop 'Bottom index wrong'

            do k=kbpl+1,nvrt
              if(ztmp(k,i)-ztmp(k-1,i)<=0) then
                write(*,*)'Inverted z-level:',eta2(i),dp(i),ztmp(k,i)-ztmp(k-1,i)
                stop
              endif
            enddo !k
          endif !dry/wet
        enddo !i=1,np
!---------------------------------------------------------
        endif !bypass
          
        if(i23d==2) then !2D
          do i=1,np
            do m=1,ivs
              read(63,rec=irec+1) out(i,1,m)
              irec=irec+1
            enddo !m
          enddo !i

          !Interpolate
          outnc=-9999 !init.
          do i=1,nx
            do j=1,ny
              if(iparen(i,j)==0) cycle
              ie=iparen(i,j)
              outnc(i,j)=dot_product(out(elnode(1:3,ie),1,1),arco(:,i,j))
            enddo !j
          enddo !i

          !Debug
!          do i=1,nx
!            do j=1,ny
!              write(97,*)i,real(xlon(i)),real(ylat(j)),outnc(i,j)
!            enddo !j
!          enddo !i

          !Output nc
          data_start_3d(1:2)=1; data_start_3d(3)=(iday-1)*nrec+it1
          data_count_3d(1)=nx; data_count_3d(2)=ny; data_count_3d(3)=1
          iret=nf_put_vara_real(ncid,ie_id,data_start_3d,data_count_3d,outnc(:,:))

        else !if(i23d==3) then !3D 
          do i=1,np
            do k=max0(1,kbp00(i)),nvrt
              do m=1,ivs
                read(63,rec=irec+1) out(i,k,m)
                irec=irec+1
              enddo !m
            enddo !k
            do k=max0(1,kbp00(i)),nvrt
!             Output: time, node #, level #, z-coordinate, 3D variable 
              if(idry(i)==1) then
!                write(65,*)time/86400,i,k,-99.,(out(i,k,m),m=1,ivs)
              else
!                write(65,*)time/86400,i,k,ztmp(k,i),(out(i,k,m),m=1,ivs)
              endif
            enddo !k

            !Debug
            !write(65,*)x(i),y(i),out(i,1,1:ivs) 

          enddo !i
!        else if(i23d==4) then !3D side @ whole level
!          do i=1,ns
!            !Compute zcor
!            n1=isidenode(1,i); n2=isidenode(2,i)
!            zs(:,i)=(ztmp(:,n1)+ztmp(:,n2))/2
!            do k=max0(1,kbs(i)),nvrt
!              do m=1,ivs
!                read(63,rec=irec+1) outs(m,k,i)
!                irec=irec+1
!              enddo !m
!          
!              !Output
!              write(65,*)time/86400,i,k,zs(k,i),outs(1:ivs,k,i)
!            enddo !k
!          enddo !i
!        else !3D elem. @whole/half level
!          do i=1,ne
!            !Compute zcor
!            n1=elnode(1,i); n2=elnode(2,i); n3=elnode(3,i)
!            ze(:,i)=(ztmp(:,n1)+ztmp(:,n2)+ztmp(:,n3))/3
!            do k=max0(1,kbe(i)),nvrt
!              do m=1,ivs
!                read(63,rec=irec+1) oute(m,k,i)
!                irec=irec+1
!              enddo !m
!
!              !Output
!              write(65,*)time/86400,i,k,ze(k,i),oute(1:ivs,k,i)
!            enddo !k
!          enddo !i
        endif !i23d
      enddo !it1=1,nrec
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      enddo !iday=1,ndays

      iret = nf_close(ncid)

      stop
      end

      subroutine check_err(stat)
      integer stat
      include 'netcdf.inc'
      if (stat .ne. NF_NOERR) then
      print *, nf_strerror(stat)
      stop
      endif
      end

!     rlambda0,phi0: in degrees
      subroutine cpp(ifl,xin,yin,xout,yout,rlambda0,phi0)
      implicit real*8(a-h,o-z)
      parameter(r=6378206.4)
      parameter(pi=3.1415926)
      integer, intent(in) :: ifl
      real*8, intent(in) :: xin,yin,rlambda0,phi0
      real*8, intent(out) :: xout,yout

      rlambda=rlambda0*pi/180
      phi=phi0*pi/180
      if(ifl.eq.1) then !lat/lon to CPP
        xin1=xin*pi/180
        yin1=yin*pi/180
        xout=r*(xin1-rlambda)*cos(phi)
        yout=yin1*r
      else !ifl=-1
        xout=rlambda+xin/r/cos(phi)
        yout=yin/r
        xout=xout/pi*180
        yout=yout/pi*180
      endif

      return
      end

