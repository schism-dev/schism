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
!       For scribe I/O
!	Read in (x,y,z) from station.bp or station.sta (sta format; x,y, with z), 
!       where z is either distance from F.S. or z-coord. If below bottom or above F.S., 
!       const. extrapolation is used, except if z=1.e10, in which case
!       depth averaged value will be calculated (for 3D vars).
!       Output time series for 3D variables (surface values for 2D variables), DEFINED AT NODES OR
!       ELEMENTS.

!       Inputs: 
!              (1) screen; 
!              (2) station.bp or station.sta
!              (3) vgrid.in: in this dir or ../
!              (4) out2d*.nc
!              (4) nc outputs for that variable(tri-quad)

!       Outputs: fort.1[89]; fort.21 (magnitude), fort.22 (dir in deg in math convention); fort.20 - local depth for each pt.
!       For ics=2 (e.g. for lon/lat), use nearest node for output
!											
!   ifort -mcmodel=medium -assume byterecl -CB -O2 -o read_output10_xyz.exe ../UtilLib/extract_mod2.f90  ../UtilLib/compute_zcor.f90 ../UtilLib/pt_in_poly_test.f90 read_output10_xyz.f90 -I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff
!****************************************************************************************
!
      program read_out
      use netcdf
      use extract_mod2
      use compute_zcor
      use pt_in_poly_test

      character(len=30) :: file63,file62,varname,varname2,file64
      character(len=12) :: it_char
!      character(len=48) :: version,variable_nm,variable_dim
      logical :: lexist,lexist2
      dimension swild(3)
      integer, allocatable :: i34(:),elnode(:,:),node3(:,:),kbp(:),iep(:),irank_read(:),kbp00(:)
      real*8,allocatable :: timeout(:),xnd(:),ynd(:)
      real,allocatable :: dp(:),ztot(:),sigma(:),sigma_lcl(:,:),outvar(:,:,:), &
     &out2(:,:,:),eta2(:),arco(:,:),ztmp(:),x00(:),y00(:), &
     &out3(:,:),z00(:),rl2min(:),dep(:),ztmp2(:,:),sum1(:)
      integer :: nodel(3),dimids(100),idims(100)
      integer :: char_len,start_2d(2),start_3d(3),start_4d(4), &
     &count_2d(2),count_3d(3),count_4d(4)
!      real*8 :: h0
      
      pi=3.1515926

      print*, 'Input extraction pts format (1: .bp; 2:.sta):'
      read(*,*)ibp
      if(ibp/=1.and.ibp/=2) stop 'Unknown format'

      print*, 'Input ics (1-linear interp; 2-nearest neighbor interp. 2 for node-based variables only! 2 is suggested for sub-meter resolution!):'
      read(*,*)ics

      print*, 'Input variable name to read from nc (e.g. elevation):'
      print*, '(if vector, input X component)'
      read(*,'(a30)')varname
      varname=adjustl(varname); len_var=len_trim(varname)
      
      print*, 'Is the var scalar (1) or vector (2)?'
      read(*,*) ivs
      if(ivs/=1.and.ivs/=2) stop 'check ivs'

      print*, 'Is the var node (1) or elem (2) based?'
      read(*,*) inode_elem

      print*, 'Input start and end file # to read:'
      read(*,*) iday1,iday2

      print*, 'Is the z-coord. in station.* relative to surface (1) or a fixed level (0)?'
      read(*,*) ifs
!'
      if(ibp==1) then !.bp format
        open(10,file='station.bp',status='old')
        read(10,*) 
        read(10,*) nxy
      else !.sta format
        open(10,file='station.sta',status='old')
        read(10,*) nxy
      endif !ibp

      allocate(x00(nxy),y00(nxy),z00(nxy),rl2min(nxy),dep(nxy),stat=istat)
      if(istat/=0) stop 'Failed to allocate (1)'

      do i=1,nxy
        if(ibp==1) then !.bp format
          read(10,*)j,x00(i),y00(i),z00(i)
        else !.sta format
          read(10,*)
          read(10,*)x00(i),y00(i),z00(i) 
        endif !ibp
      enddo !i
      close(10)

!...  Get basic info from out2d*.nc
      !Returned vars: ne,np,ns,nrec,[xnd ynd dp](np),
      !elnode,i34,nvrt,h0,dtout,kbp
      call get_dims(iday1,np,ne,ns,nvrt,h0)
      allocate(xnd(np),ynd(np),dp(np),kbp(np),i34(ne),elnode(4,ne),stat=istat)
      if(istat/=0) stop 'alloc (1)'
      call readheader(iday1,np,ne,ns,kbp,i34,elnode,nrec,xnd,ynd,dp,dtout)

      print*, 'After header:',ne,np,ns,nvrt,nrec,i34(ne), &
     &elnode(1:i34(ne),ne),h0,xnd(np),ynd(np),dp(np),dtout !,start_time

!...  Read in time records in segments for mem
      if(inode_elem==1) then !node based
        last_dim=np
      else !elem
        last_dim=ne
      endif

      allocate(timeout(nrec),ztot(nvrt),sigma(nvrt),sigma_lcl(nvrt,np),kbp00(np), &
     &outvar(nvrt,last_dim,ivs),eta2(np),out2(nxy,nvrt,ivs),out3(nxy,ivs),sum1(ivs), &
     &node3(nxy,3),arco(nxy,3),iep(nxy),ztmp(nvrt),ztmp2(nvrt,3))
      outvar=-huge(1.0) !test mem

!     Read in vgrid
      call get_vgrid_single('vgrid.in',np,nvrt,ivcor,kz,h_s,h_c,theta_b, &
     &theta_f,ztot,sigma,sigma_lcl,kbp)

!     Calculate kbp00 
      if(ivcor==1) then
        kbp00=kbp
      else
        do i=1,np
          !Use large eta to get true bottom
          call zcor_SZ_single(dp(i),1.e8,h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot,sigma, &
     &ztmp(:),idry2,kbp00(i))
        enddo !i
      endif !ivcor

!...  Find parent element for (x00,y00) 
      iep=0
      if(ics==1) then !Cartesian
        do i=1,ne
          do l=1,nxy
            if(iep(l)/=0) cycle
            call pt_in_poly_single(i34(i),real(xnd(elnode(1:i34(i),i))), &
     &real(ynd(elnode(1:i34(i),i))),x00(l),y00(l),inside,arco(l,1:3),nodel)
            if(inside==1) then
              iep(l)=i
              !print*, 'Found:',l,arco(l,1:3),nodel
              node3(l,1:3)=elnode(nodel(1:3),i)
            endif !inside
          enddo !l; build pts

          ifl=0 !flag
          do l=1,nxy
            if(iep(l)==0) then
              ifl=1
              exit
            endif
          enddo !l
          if(ifl==0) exit
        enddo !i=1,ne
      else !lat/lon; needed for node-based only
        rl2min=1.e25 !min distance^2
        do ie=1,ne
          do j=1,i34(ie)
            i=elnode(j,ie)
            do l=1,nxy
              rl2=(xnd(i)-x00(l))**2+(ynd(i)-y00(l))**2
              if(rl2<rl2min(l)) then
                rl2min(l)=rl2
                iep(l)=ie
                node3(l,1:3)=i
                arco(l,1:3)=1./3
              endif
            enddo !l=1,nxy
          enddo !j
        enddo !i=1,np
      endif !ics

      do j=1,nxy
        if(iep(j)<=0) then
          print*, 'Cannot find a parent for pt:',j,x00(j),y00(j)
          stop
        endif
      enddo !j

!...  Time iteration
!...
      do iday=iday1,iday2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      write(it_char,'(i12)')iday
      it_char=adjustl(it_char)
      leng=len_trim(it_char)
      file62='out2d_'//it_char(1:leng)//'.nc'
      iret=nf90_open(trim(adjustl(file62)),OR(NF90_NETCDF4,NF90_NOWRITE),ncid4)
      !time is double
      iret=nf90_inq_varid(ncid4,'time',itime_id)
      iret=nf90_get_var(ncid4,itime_id,timeout,(/1/),(/nrec/))
      print*, 'time=',timeout !,trim(adjustl(file63))
 
      !Find nc file
      file63=varname(1:len_var)//'_'//it_char(1:leng)//'.nc'
      if(ivs==2) varname2=varname(1:len_var-1)//'Y'
      inquire(file=file63,exist=lexist)
      if(lexist) then !3D var
        i23d=2 
        if(ivs==2) then
          file64=varname2(1:len_var)//'_'//it_char(1:leng)//'.nc'
          inquire(file=file64,exist=lexist2)
          if(.not.lexist2) then
            print*, 'Missing y-component:',file64
            print*, file63
            stop
          endif
        endif
      else
        i23d=1 !2D
        file63=file62
        file64=file62
      endif

      iret=nf90_open(trim(adjustl(file63)),OR(NF90_NETCDF4,NF90_NOWRITE),ncid)
      iret=nf90_inq_varid(ncid,varname(1:len_var),ivarid1)
      if(iret/=nf90_NoErr) stop 'Var not found'
      iret=nf90_Inquire_Variable(ncid,ivarid1,ndims=ndims,dimids=dimids)
      if(ndims>100) stop 'increase dimension of dimids & idims'
      do i=1,ndims
       iret=nf90_Inquire_Dimension(ncid,dimids(i),len=idims(i))
      enddo !i
      npes=idims(ndims-1) !np|ne|ns
      if(npes/=np.and.npes/=ne) stop 'can only handle node- or elem-based'
!'
      if(idims(ndims)/=nrec) stop 'last dim is not time'

      if(ivs==2) then !vector
        iret=nf90_open(trim(adjustl(file64)),OR(NF90_NETCDF4,NF90_NOWRITE),ncid2)
        if(iret/=nf90_NoErr) stop 'Failed to open file64'
        iret=nf90_inq_varid(ncid2,varname2(1:len_var),ivarid2)
        if(iret/=nf90_NoErr) stop 'Var2 not found'
      endif !ivs

      do irec=1,nrec
        !Get elev
        iret=nf90_inq_varid(ncid4,'elevation',itmp)
        start_2d(1)=1; start_2d(2)=irec
        count_2d(1)=np; count_2d(2)=1
        iret=nf90_get_var(ncid4,itmp,eta2,start_2d,count_2d)

        if(i23d==1) then !2D
          start_2d(1)=1; start_2d(2)=irec
          count_2d(1)=npes; count_2d(2)=1
          iret=nf90_get_var(ncid,ivarid1,outvar(1,1:npes,1),start_2d,count_2d)
          if(ivs==2) iret=nf90_get_var(ncid2,ivarid2,outvar(1,1:npes,2),start_2d,count_2d)
        else !3D
          start_3d(1:2)=1; start_3d(3)=irec
          count_3d(1)=nvrt; count_3d(2)=npes; count_3d(3)=1
          iret=nf90_get_var(ncid,ivarid1,outvar(:,1:npes,1),start_3d,count_3d)
          if(ivs==2) iret=nf90_get_var(ncid2,ivarid2,outvar(:,1:npes,2),start_3d,count_3d)
        endif 

        !Available now: outvar(nvrt,np|ne,ivs), eta2(np)
        out2=0
        out3=0
        if(i23d==1) then !2D
          do i=1,nxy
            dep(i)=0
            do j=1,3 !nodes
              nd=node3(i,j)
              !Compute local depth
              dep(i)=dep(i)+arco(i,j)*dp(nd)
              if(inode_elem==1) then !node
                out2(i,1,:)=out2(i,1,:)+arco(i,j)*outvar(1,nd,:)
              else if (inode_elem==2) then !elem
                if(iep(i)<=0) stop 'iep(i)<=0'
                out2(i,1,:)=outvar(1,iep(i),:)
              endif
            enddo !j
          enddo !i
          write(18,'(e16.8,20000(1x,e14.6))')timeout(irec)/86400,(out2(i,1,1),i=1,nxy)
          if(ivs==2) then
            write(19,'(e16.8,20000(1x,e14.6))')timeout(irec)/86400,(out2(i,1,2),i=1,nxy)
            write(21,'(e16.8,20000(1x,f14.6))')timeout(irec)/86400,(sqrt(out2(i,1,1)**2+out2(i,1,2)**2),i=1,nxy)
            write(22,'(e16.8,20000(1x,f14.6))')timeout(irec)/86400,(atan2(out2(i,1,2),out2(i,1,1))/pi*180,i=1,nxy)
          endif
        else !3D
!         Do interpolation
          do i=1,nxy
            etal=0; dep(i)=0; idry=0
            do j=1,3
              nd=node3(i,j)
              if(eta2(nd)+dp(nd)<h0) idry=1
              etal=etal+arco(i,j)*eta2(nd)
              dep(i)=dep(i)+arco(i,j)*dp(nd)
!             Debug
!              write(11,*)i,j,nd,dp(nd),arco(i,j)
            enddo !j
            if(idry==1) then
              out3(i,:)=-99
!              write(65,*)'Dry'
            else !element wet
              !Compute z-coordinates
              if(ivcor==1) then !localized
                !Strictly speaking for elem-based vars, we need to use i34 nodes
                do j=1,3
                  nd=node3(i,j)
                  do k=kbp(nd)+1,nvrt-1
                    ztmp2(k,j)=(eta2(nd)+dp(nd))*sigma_lcl(k,nd)+eta2(nd)
                  enddo !k
                  ztmp2(kbp(nd),j)=-dp(nd) !to avoid underflow
                  ztmp2(nvrt,j)=eta2(nd) !to avoid underflow
                enddo !j

                ztmp=0
                kbpl=minval(kbp(node3(i,1:3)))
                do k=kbpl,nvrt
                  do j=1,3
                    nd=node3(i,j)
                    ztmp(k)=ztmp(k)+arco(i,j)*ztmp2(max(k,kbp(nd)),j)
                  enddo !j
                enddo !k
              else if(ivcor==2) then !SZ
                call zcor_SZ_single(dep(i),etal,h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot,sigma, &
     &ztmp(:),idry2,kbpl)
              endif

!             Horizontal interpolation
              if(inode_elem==1) then !node based
                do k=kbpl,nvrt
                  do j=1,3
                    nd=node3(i,j)
                    kin=max(k,kbp00(nd))
                    out2(i,k,:)=out2(i,k,:)+arco(i,j)*outvar(kin,nd,:)
                  enddo !j
                enddo !k
              endif !inode_elem

              if(abs(z00(i)-1.e10)<0.1) then !depth average
                total_dp=ztmp(nvrt)-ztmp(kbpl)
                if(total_dp<=0) then
                  write(*,*)'depth<=0:',total_dp,i
                  stop
                endif
                sum1=0
                do k=kbpl,nvrt-1
                  if(inode_elem==1) then !node
                    sum1(:)=sum1(:)+(ztmp(k+1)-ztmp(k))*(out2(i,k,:)+out2(i,k+1,:))/2
                  else if(inode_elem==2) then !elem
                    sum1(:)=sum1(:)+(ztmp(k+1)-ztmp(k))*outvar(k+1,iep(i),:)
                  endif    
                enddo !k
                out3(i,:)=sum1(:)/total_dp
              else !Interplate in vertical
                if(ifs==0) then !relative to MSL
                  z2=z00(i)
                else
                  z2=ztmp(nvrt)-z00(i)
                endif
                if(z2>=ztmp(nvrt)) then !above F.S.
                  k0=nvrt-1; rat=1
                else if(z2<=ztmp(kbpl)) then !below bottom; extrapolate
                  k0=kbpl; rat=0
                else !above bottom; cannot be above F.S.
                  k0=0
                  do k=kbpl,nvrt-1
                    if(z2>=ztmp(k).and.z2<=ztmp(k+1)) then
                      k0=k
                      rat=(z2-ztmp(k))/(ztmp(k+1)-ztmp(k))
                      exit
                    endif
                  enddo !k
                endif !ztmp

                if(k0==0) then
                  write(*,*)'read_output7b_xyz: failed to find a vertical level:',irec_real,i,ifs,z2,ztmp(:)
!'
                  stop
                else
                  if(inode_elem==1) then !node
                    out3(i,:)=out2(i,k0,:)*(1-rat)+out2(i,k0+1,:)*rat
                  else if(inode_elem==2) then !elem
                    if(iep(i)<=0) stop 'iep(i)<=0(2)'
                    out3(i,:)=outvar(k0+1,iep(i),:) !FV
                  endif
                endif
              endif !depth average or not
            endif !dry/wet
          enddo !i=1,nxy
          write(18,'(e16.8,20000(1x,f14.6))')timeout(irec)/86400,(out3(i,1),i=1,nxy)
          if(ivs==2) then
            write(19,'(e16.8,20000(1x,f14.6))')timeout(irec)/86400,(out3(i,2),i=1,nxy)
            write(21,'(e16.8,20000(1x,f14.6))')timeout(irec)/86400,(sqrt(out3(i,1)**2+out3(i,2)**2),i=1,nxy)
            write(22,'(e16.8,20000(1x,f14.6))')timeout(irec)/86400,(atan2(out3(i,2),out3(i,1))/pi*180,i=1,nxy)
          endif
         
        endif !i23d
      enddo !irec
      iret=nf90_close(ncid)
      iret=nf90_close(ncid4)
      if(ivs==2) iret=nf90_close(ncid2)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      enddo !iday

!     Output local depths info
      do i=1,nxy
        write(20,*)i,dep(i)
      enddo !i

      print*, 'Finished!'

      end program read_out
