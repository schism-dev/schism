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
!	(x,y) read in from station.bp (build pts), and output
!       results along a transect (defined below) for 3D variables. 
!       Need to manually modify ntran etc. 
!       Works for mixed tri/quad outputs and node/elem based vars from
!       scribe I/O versions.
!       Will extrapolate above surface but not below bottom.

!       Inputs: screen; vgrid (in this dir or ../); station.bp (build pts; depths not used); relevant nc and out2d*.nc
!       Outputs: transect.out & transect_grd.[zr]0 (ascii on struc'ed grid; 
!                use plot_transect.m, but also consider using SCHISM_TRANSECT.m)
!                average_transect.out (averaged; in .bp with header) 

! ifort -mcmodel=medium -assume byterecl -CB -O2 -o read_output10_transect.exe ../UtilLib/extract_mod2.f90 ../UtilLib/compute_zcor.f90 ../UtilLib/pt_in_poly_test.f90 read_output10_transect.f90 -I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -L$NETCDF/lib -lnetcdf -lnetcdff
!											*
!****************************************************************************************
!
      program read_out
      use netcdf
      use extract_mod2
      use compute_zcor
      use pt_in_poly_test

!      parameter(nbyte=4)
      character(len=30) :: file63,varname,file62
      character(len=12) :: it_char
      logical::lexist
      allocatable :: sigma(:),cs(:),ztot(:)
      allocatable:: outvar(:,:),out(:,:,:),out2(:,:),eta2(:),node3(:,:),arco(:,:)
      allocatable :: ztmp(:),x00(:),y00(:),iep(:),out3(:,:),out4(:),av_out3(:,:)
      allocatable :: nmxz(:,:),z0(:),r0(:),r00(:),z00(:),dp(:)
      allocatable :: sigma_lcl(:,:),kbp(:),kbp00(:),ztmp2(:,:),irank_read(:)
      integer :: nodel(3),dimids(100),idims(100)
      integer :: char_len,start_2d(2),start_3d(3),start_4d(4), &
     &count_2d(2),count_3d(3),count_4d(4)
      integer, allocatable :: i34(:),elnode(:,:)

      real :: swild(3)
      real*8,allocatable :: timeout(:),xnd(:),ynd(:)
      
      print*, 'Input variable name to read from nc (e.g. salinity):'
      read(*,'(a30)')varname
      varname=adjustl(varname); len_var=len_trim(varname)

      print*, 'Is the var node (1) or elem (2) based?'
      read(*,*) inode_elem
      
      print*, 'Input start and end file # to read:'
      read(*,*) iday1,iday2

!     Input transect depths
      ntran=201
      allocate(z0(ntran))
      z0(1:ntran)=(/(-200+i, i=0,ntran-1) /)

      open(10,file='station.bp',status='old')
      read(10,*) 
      read(10,*) nxy
      nxz=nxy*ntran
      nxz_e=(nxy-1)*(ntran-1)*2
      allocate(x00(nxy),y00(nxy),r0(nxy),r00(nxz),z00(nxz),nmxz(nxz_e,4),stat=istat)
      if(istat/=0) stop 'Falied to allocate (1)'

      do i=1,nxy
        read(10,*)j,x00(i),y00(i)
      enddo !i
      !Along transect distance
      r0(1)=0
      do i=2,nxy
        r0(i)=r0(i-1)+sqrt((x00(i)-x00(i-1))**2+(y00(i)-y00(i-1))**2)
      enddo !i
      close(10)
 
!     Output transect grid
      open(22,file='transect_grd.z0',status='replace')
      open(23,file='transect_grd.r0',status='replace')
      write(22,'(f12.3)')z0
      write(23,'(e15.6)')r0
      close(22) 
      close(23)

!     Compute connectivity
      do i=1,nxy
        do j=1,ntran
          nd=ntran*(i-1)+j
          r00(nd)=r0(i)
          z00(nd)=z0(j)
        enddo !j
      enddo !i

      do i=1,nxy-1
        do j=1,ntran-1
          ie=(ntran-1)*(i-1)+j
          n1=ntran*(i-1)+j
          n2=ntran*i+j
          nmxz(2*ie-1,1)=n1
          nmxz(2*ie-1,2)=n2
          nmxz(2*ie-1,3)=n2+1
          nmxz(2*ie,1)=n1
          nmxz(2*ie,2)=n2+1
          nmxz(2*ie,3)=n1+1
        enddo !j
      enddo !i

      open(21,file='transect.out',status='replace')
      
!...  Get basic info from out2d*.nc
      !Returned vars: ne,np,ns,nrec,[xnd ynd dp](np),
      !elnode,i34,nvrt,h0,dtout,kbp
      call get_dims(iday1,np,ne,ns,nvrt,h0)
      allocate(xnd(np),ynd(np),dp(np),kbp(np),i34(ne),elnode(4,ne),stat=istat)
      if(istat/=0) stop 'alloc (1)'
      call readheader(iday1,np,ne,ns,kbp,i34,elnode,nrec,xnd,ynd,dp,dtout)

      print*, 'After header:',ne,np,ns,nvrt,nrec,i34(ne), &
     &elnode(1:i34(ne),ne),h0,xnd(np),ynd(np),dp(np),dtout !,start_time

!      if(iskip<1) then
!        print*, 'Output interval too small; resetting to dtout=',dtout
!        iskip=1
!      endif
      !nstep=(iday2-iday1+1)*nrec/iskip
!      if(nxz*ivs>6000) stop 'Increase output statement below!'
 
      allocate(out(nxy,3,nvrt),out2(nxy,nvrt),node3(nxy,3),arco(nxy,3),iep(nxy), &
     &out3(ntran,nxy),av_out3(ntran,nxy),out4(ntran),stat=istat)
      if(istat/=0) stop 'Falied to allocate (5)'

      if(inode_elem==1) then !node based
        last_dim=np
      else !elem
        last_dim=ne
      endif

      allocate(ztot(nvrt),sigma(nvrt),sigma_lcl(nvrt,np),kbp00(np), &
     &timeout(nrec),outvar(nvrt,last_dim),eta2(np))
      outvar=-huge(1.0) !init

!     Read in vgrid.in
      call get_vgrid_single('vgrid.in',np,nvrt,ivcor,kz,h_s,h_c,theta_b,theta_f,ztot,sigma,sigma_lcl,kbp)

      allocate(ztmp(nvrt),ztmp2(nvrt,3))

!     Calculate kbp00
      if(ivcor==1) then
        kbp00=kbp
      else
        do i=1,np
          call zcor_SZ_single(dp(i),1.e8,h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot,sigma, &
     &ztmp(:),idry2,kbp00(i))
        enddo !i
      endif !ivcor

!...  Find parent element for (x00,y00)
      iep=0
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

      iabort=0
      do j=1,nxy
        if(iep(j)<=0) then
          print*, 'Cannot find a parent for pt:',j,x00(j),y00(j)
          iabort=1
        endif
      enddo !j
      if(iabort==1) stop 'check station points'

      if(varname(1:len_var).eq.'hvel') then
        rjunk=0 !invalid values
      else
        rjunk=-9999
      endif

!...  Time iteration
!...
      av_out3=0
      nstep=0
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
      inquire(file=file63,exist=lexist)
      if(lexist) then
        i23d=2 !3D var
      else
        i23d=1 !2D
        file63=file62
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

!        iret=nf90_get_att(ncid,ivarid1,'i23d',i23d)
!        if(i23d<=0.or.i23d>6) stop 'wrong i23d'
!        if(i23d>3.and.ics==2) stop 'ics=2 with elem-based var'
!        iret=nf90_get_att(ncid,ivarid1,'ivs',ivs)
!        !print*, 'i23d:',i23d,ivs,idims(1:ndims)
!
      do irec=1,nrec
        !Get elev
        iret=nf90_inq_varid(ncid4,'elevation',itmp)
        start_2d(1)=1; start_2d(2)=irec
        count_2d(1)=np; count_2d(2)=1
        iret=nf90_get_var(ncid4,itmp,eta2,start_2d,count_2d)

        if(i23d==1) then !2D
          start_2d(1)=1; start_2d(2)=irec
          count_2d(1)=npes; count_2d(2)=1
          iret=nf90_get_var(ncid,ivarid1,outvar(1,1:npes),start_2d,count_2d)
        else !3D
          start_3d(1:2)=1; start_3d(3)=irec
          count_3d(1)=nvrt; count_3d(2)=npes; count_3d(3)=1
          iret=nf90_get_var(ncid,ivarid1,outvar(:,1:npes),start_3d,count_3d)
        endif 

        !Available now: outvar(nvrt,np|ne), eta2(np)
        out2=0
        if(i23d==1) then !2D
          stop 'No 2D vars plz'
        else !3D
          if(inode_elem==1) then !node
            do i=1,nxy
              do j=1,3 !nodes
                nd=node3(i,j)
                do k=max0(1,kbp00(nd)),nvrt
                  out(i,j,k)=outvar(k,nd)
                enddo !k
              enddo !j
            enddo !i
          endif !node

!         Do interpolation
          do i=1,nxy
            etal=0; dep=0; idry=0
            do j=1,3
              nd=node3(i,j)
              if(eta2(nd)+dp(nd)<h0) idry=1
              etal=etal+arco(i,j)*eta2(nd)
              dep=dep+arco(i,j)*dp(nd)
      
!             Debug
!              write(11,*)i,j,nd,dp(nd),arco(i,j)

            enddo !j
            if(idry==1) then
              out4(:)=rjunk
            else !element wet
              !Compute z-coordinates
              if(ivcor==1) then !localized
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
                call zcor_SZ_single(dep,etal,h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot, &
     &sigma,ztmp(:),idry2,kbpl)
              endif

              if(inode_elem==1) then !node based
                do k=kbpl,nvrt
                  do j=1,3
                    nd=node3(i,j)
                    kin=max(k,kbp00(nd))
                    out2(i,k)=out2(i,k)+arco(i,j)*out(i,j,kin)
                  enddo !j
!                  write(65,*)i,k,ztmp(k),out2(i,k)
                enddo !k
              endif

!             Interplate in vertical
              do kk=1,ntran
                if(z0(kk)>=ztmp(nvrt)) then !extrap above F.S.
                  k0=nvrt-1; rat=1
                else !no extrap below bottom
                  k0=0
                  do k=kbpl,nvrt-1
                    if(z0(kk)>=ztmp(k).and.z0(kk)<=ztmp(k+1)) then
                      k0=k
                      rat=(z0(kk)-ztmp(k))/(ztmp(k+1)-ztmp(k))
                      exit
                    endif
                  enddo !k
                endif 

                if(k0==0) then !no extrap below bottom
!                  write(12,*)'Warning: failed to find a vertical level:',it,i
                  out4(kk)=rjunk
                else
                  if(i23d<=3) then !node based
                    out4(kk)=out2(i,k0)*(1-rat)+out2(i,k0+1)*rat
                  else if(i23d<=6) then !elem
                    out4(kk)=outvar(k0+1,iep(i)) !FV
                  endif
        
                  !Debug
                  !if(mod(its,iskip)==0) write(99,*)time,i,kk,out4(kk,m)
                endif
              enddo !kk
            endif !dry/wet

            do kk=1,ntran
              out3(kk,i)=out4(kk)
            enddo !kk
          enddo !i=1,nxy

          !Along all vertical levels of each build pt
          do i=1,nxy
            do kk=1,ntran
              write(21,'(e22.12,2(1x,f10.3))')timeout(irec),out3(kk,i)
            enddo !kk
          enddo !i

          av_out3=av_out3+out3
          nstep=nstep+1
        endif !i23d
!----------------------------------------------------------------------------
      enddo !irec
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      enddo !iday
 
!     Output average
      open(24,file='average_transect.out',status='replace')
      av_out3=av_out3/nstep
      do i=1,nxy
        do kk=1,ntran
          j=(i-1)*ntran+kk
          write(24,'(i13,4(1x,e22.10))')j,r0(i),z0(kk),av_out3(kk,i)
        enddo !kk
      enddo !i
      close(24)

      print*, 'Finished!'

      stop
      end

