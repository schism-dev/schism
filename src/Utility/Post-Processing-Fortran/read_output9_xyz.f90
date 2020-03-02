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
!											
!	Read in (x,y,z) from station.bp or station.sta (sta format; x,y, with z), 
!       where z is either distance from F.S. or z-coord. If below bottom or above F.S., 
!       const. extrapolation is used, except if z=1.e10, in which case
!       depth averaged value will be calculated (for 3D vars).
!       Output time series for 3D variables (surface values for 2D variables), DEFINED AT NODES OR
!       ELEMENTS.
!       Works with combined or uncombined nc outputs.

!       Inputs: 
!              (1) screen; 
!              (2) station.bp or station.sta
!              (3) vgrid.in: in this dir or ../
!              (4) combined or uncombined nc outputs (schout*.nc; tri-quad)

!              max_array_size - this constant sets max
!              array size and may need to be adjusted depending on your
!              systems memory.
!       Outputs: fort.1[89]; ; fort.20 - local depth for each pt.
!       For ics=2 (e.g. for lon/lat), use nearest node for output
!											
!   ifort -mcmodel=medium -assume byterecl -CB -O2 -o read_output9_xyz.exe ../UtilLib/extract_mod.f90 \
! ../UtilLib/compute_zcor.f90 ../UtilLib/pt_in_poly_test.f90 read_output9_xyz.f90 -I$NETCDF/include \
!-I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff
!****************************************************************************************
!
      program read_out
      use netcdf
      use extract_mod
      use compute_zcor
      use pt_in_poly_test

      character(len=30) :: file63,varname
      character(len=12) :: it_char
!      character(len=48) :: version,variable_nm,variable_dim
      logical :: lexist
      dimension swild(3)
      integer, allocatable :: node3(:,:),kbp(:),iep(:),irank_read(:)
      real*8,allocatable :: timeout(:)
      real,allocatable :: ztot(:),sigma(:),sigma_lcl(:,:),outvar(:,:,:,:), &
     &out(:,:,:,:),out2(:,:,:),eta2(:,:),arco(:,:),ztmp(:),x00(:),y00(:), &
     &out3(:,:),z00(:),rl2min(:),dep(:),ztmp2(:,:)
      integer :: nodel(3)
      
      print*, 'Do you work on uncombined (0) or combined (1) nc?'
      read(*,*)icomb
      if(icomb/=0.and.icomb/=1) stop 'Unknown icomb'

!...  Set max array size for system memory
      print*, 'Recommendation: for uncombined nc, specify max array size (e.g., <=2.e9);'
      print*, 'for combined nc, specify # of records to read each time.'
      print*, 'Do you want to specify max array size (1) or # of records (2)'
      read(*,*)ispec_max
      if(ispec_max==1) then
        print*, 'Input max array size (e.g., <=2.e9):'
        read(*,*)max_array_size
        print*, 'max_array_size read in=',max_array_size
      else
        print*, 'Input # of records:'
        read(*,*)nrec3
      endif

      print*, 'Input extraction pts format (1: .bp; 2:.sta):'
      read(*,*)ibp
      if(ibp/=1.and.ibp/=2) stop 'Unknown format'

      print*, 'Input ics (1-linear interp; 2-nearest neighbor interp. 2 for node-based variables only! 2 is suggested for sub-meter resolution!):'
      read(*,*)ics

      print*, 'Input variable name to read from nc (e.g. elev):'
      read(*,'(a30)')varname
      varname=adjustl(varname); len_var=len_trim(varname)
      
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

!...  Header
      !Returned vars: ne,np,ns,nrec,[x y dp](np),
      !elnode,i34,nvrt,h0,dtout
      !If icomb=0, additonal vars: nproc,iegl_rank,iplg,ielg,islg,np_lcl(:),ne_lcl(:),ns_lcl(:)
      if(icomb==0) then !uncombined
        call get_global_geo
      else
        call readheader(iday1)
      endif

      print*, 'After header:',ne,np,ns,nrec,i34(ne), &
     &elnode(1:i34(ne),ne),nvrt,h0,x(np),y(np),dp(np) !,start_time

!...  Read in time records in segments for mem
      if(inode_elem==1) then !node based
        last_dim=np
      else !elem
        last_dim=ne
      endif

      !nrec specified if ispec_max=2
      if(ispec_max==1) nrec3=max_array_size/(2*nvrt*last_dim)-1
      nrec3=min(nrec,max(nrec3,1))

      print*, '# of records in each stack is: ',nrec, &
     &'; actual # of records that can be read in each time is:',nrec3
 
      allocate(timeout(nrec),ztot(nvrt),sigma(nvrt),sigma_lcl(nvrt,np), &
     &kbp(np),outvar(2,nvrt,last_dim,nrec3),eta2(np,nrec3),out2(nxy,nvrt,2),out3(nxy,2), &
     &out(nxy,3,nvrt,2),node3(nxy,3),arco(nxy,3),iep(nxy))
      call get_vgrid_single('vgrid.in',np,nvrt,ivcor,kz,h_s,h_c,theta_b, &
     &theta_f,ztot,sigma,sigma_lcl,kbp)
      allocate(ztmp(nvrt),ztmp2(nvrt,3),stat=istat)
      outvar=-huge(1.0) !test mem

!     Read in vgrid.in 
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

!...  Find parent element for (x00,y00) (global indices)
      iep=0
      if(ics==1) then !Cartesian
        do i=1,ne
          do l=1,nxy
            if(iep(l)/=0) cycle
            call pt_in_poly_single(i34(i),x(elnode(1:i34(i),i)),y(elnode(1:i34(i),i)),x00(l),y00(l),inside,arco(l,1:3),nodel)
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
              rl2=(x(i)-x00(l))**2+(y(i)-y00(l))**2
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

      !Mark ranks that need to be read in
      if(icomb==0) then
        allocate(irank_read(0:nproc-1))
        irank_read=0
        do j=1,nxy
          irank_read(iegl_rank(iep(j)))=1 
          write(99,*)'reading from rank #:',iegl_rank(iep(j))
        enddo
         write(99,*)'Need to read from ',sum(irank_read),' ranks'
      endif !icomb

!...  Time iteration
!...
      do iday=iday1,iday2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!      write(it_char,'(i12)')iday
!      it_char=adjustl(it_char)
!      leng=len_trim(it_char)
!      file63='schout_'//it_char(1:leng)//'.nc'
!      iret=nf90_open(trim(adjustl(file63)),OR(NF90_NETCDF4,NF90_NOWRITE),ncid)
!      !time is double
!      iret=nf90_inq_varid(ncid,'time',itime_id)
!      iret=nf90_get_var(ncid,itime_id,timeout,(/1/),(/nrec/))
      if(icomb==0) then !uncombined
        call get_timeout(iday,nrec,timeout,icomb)
      else
        call get_timeout(iday,nrec,timeout)
      endif
!      print*, 'time=',timeout !,trim(adjustl(file63))

!      print*, 'done reading at time:',timeout(:)/86400

!        iret=nf90_inq_varid(ncid,varname(1:len_var),ivarid1)
!        if(iret/=nf90_NoErr) stop 'Var not found'
!        iret=nf90_Inquire_Variable(ncid,ivarid1,ndims=ndims,dimids=dimids)
!        if(ndims>100) stop 'increase dimension of dimids & idims'
!        do i=1,ndims
!          iret=nf90_Inquire_Dimension(ncid,dimids(i),len=idims(i))
!        enddo !i
!        npes=idims(ndims-1) !np|ne|ns
!        if(npes/=np.and.npes/=ne) stop 'can only handle node- or elem-based'
!!'
!        if(idims(ndims)/=nrec) stop 'last dim is not time'
!
!        iret=nf90_get_att(ncid,ivarid1,'i23d',i23d)
!        if(i23d<=0.or.i23d>6) stop 'wrong i23d'
!        if(i23d>3.and.ics==2) stop 'ics=2 with elem-based var'
!        iret=nf90_get_att(ncid,ivarid1,'ivs',ivs)
!        !print*, 'i23d:',i23d,ivs,idims(1:ndims)
!
!        if(ivs==1) then !scalar
!          if(mod(i23d-1,3)==0) then !2D
!            start_2d(1)=1; start_2d(2)=irec
!            count_2d(1)=npes; count_2d(2)=1
!            iret=nf90_get_var(ncid,ivarid1,outvar(1,1,1:npes),start_2d,count_2d)
!          else !3D
!            start_3d(1:2)=1; start_3d(3)=irec
!            count_3d(2)=npes; count_3d(1)=nvrt; count_3d(3)=1
!            iret=nf90_get_var(ncid,ivarid1,outvar(1,:,1:npes),start_3d,count_3d)
!          endif 
!        else !vector
!          if(mod(i23d-1,3)==0) then !2D
!            start_3d(1:2)=1; start_3d(3)=irec
!            count_3d(2)=npes; count_3d(1)=2; count_3d(3)=1
!            iret=nf90_get_var(ncid,ivarid1,outvar(1:2,1,1:npes),start_3d,count_3d)
!          else if(ndims-1==3) then !3D vector
!            start_4d(1:3)=1; start_4d(4)=irec
!            count_4d(3)=npes; count_4d(2)=nvrt; count_4d(1)=2; count_4d(4)=1
!            iret=nf90_get_var(ncid,ivarid1,outvar(:,:,1:npes),start_4d,count_4d)
!          else
!            stop 'Unknown type(2)'
!          endif
!        endif !ivs


      !do irec=1,nrec
      irec1=1 !start record
      loop1: do
        irec2=min(nrec,irec1+nrec3-1)

        if(icomb==0) then !uncombined
          do irank=0,nproc-1
            if(irank_read(irank)>0) then
              call get_outvar_multirecord(ics,iday,varname,irec1,irec2,np,last_dim,nvrt,nrec3,outvar,i23d,ivs,eta2,irank)
            endif
          enddo !irank
        else
          call get_outvar_multirecord(ics,iday,varname,irec1,irec2,np,last_dim,nvrt,nrec3,outvar,i23d,ivs,eta2)
        endif
        if(inode_elem==1) then !node based
          if(i23d>3) stop 'U said it is node based'
        else !elem
          if(i23d<=3) stop 'U said it is elem based'
        endif

        !Available now: outvar(2,nvrt,np|ne,irec2-irec1+1),i23d,ivs,eta2(np,irec2-irec1+1)
        !However, for uncombined nc, values in untouched ranks are junk

        do irec=1,irec2-irec1+1 !offeset record #
!----------------------------------------------------------------------------
        irec_real=irec1+irec-1 !actual record #
        out2=0
        out3=0
        if(mod(i23d-1,3)==0) then !2D
          do i=1,nxy
            dep(i)=0
            do j=1,3 !nodes
              nd=node3(i,j)
              !Compute local depth
              dep(i)=dep(i)+arco(i,j)*dp(nd)
              do m=1,ivs
                if(i23d<=3) then !node
                  out2(i,1,m)=out2(i,1,m)+arco(i,j)*outvar(m,1,nd,irec)
                else if (i23d<=6) then !elem
                  if(iep(i)<=0) stop 'iep(i)<=0'
                  out2(i,1,m)=outvar(m,1,iep(i),irec)
                endif
              enddo !m
            enddo !j
          enddo !i
          write(18,'(e16.8,20000(1x,e14.6))')timeout(irec_real)/86400,(out2(i,1,1),i=1,nxy)
          if(ivs==2) write(19,'(e16.8,20000(1x,e14.6))')timeout(irec_real)/86400,(out2(i,1,2),i=1,nxy)
        else !3D
          if(i23d<=3) then !node
            do i=1,nxy
              do j=1,3 !nodes
                nd=node3(i,j)
                do k=max0(1,kbp00(nd)),nvrt
                  do m=1,ivs
                    out(i,j,k,m)=outvar(m,k,nd,irec)
                  enddo !m
                enddo !k
              enddo !j
            enddo !i
          endif !node

!         Do interpolation
          do i=1,nxy
            etal=0; dep(i)=0; idry=0
            do j=1,3
              nd=node3(i,j)
              if(eta2(nd,irec)+dp(nd)<h0) idry=1
              etal=etal+arco(i,j)*eta2(nd,irec)
              dep(i)=dep(i)+arco(i,j)*dp(nd)
      
!             Debug
!              write(11,*)i,j,nd,dp(nd),arco(i,j)

            enddo !j
            if(idry==1) then
              if(ivs==2) then
                out3(i,1:2)=0
              else
                out3(i,1:2)=-99
              endif
!              write(65,*)'Dry'
            else !element wet
              !Compute z-coordinates
              if(ivcor==1) then !localized
                !Strictly speaking for elem-based vars, we need to use
                !i34 nodes
                do j=1,3
                  nd=node3(i,j)
                  do k=kbp(nd)+1,nvrt-1
                    ztmp2(k,j)=(eta2(nd,irec)+dp(nd))*sigma_lcl(k,nd)+eta2(nd,irec)
                  enddo !k
                  ztmp2(kbp(nd),j)=-dp(nd) !to avoid underflow
                  ztmp2(nvrt,j)=eta2(nd,irec) !to avoid underflow
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
              if(i23d<=3) then !node based
                do k=kbpl,nvrt
                  do m=1,ivs
                    do j=1,3
                      nd=node3(i,j)
                      kin=max(k,kbp00(nd))
                      out2(i,k,m)=out2(i,k,m)+arco(i,j)*out(i,j,kin,m)
                    enddo !j
                  enddo !m
!                  write(65,*)i,k,ztmp(k),(out2(i,k,m),m=1,ivs)     
                enddo !k
              endif !i23d<=3

              if(abs(z00(i)-1.e10)<0.1) then !depth average
                total_dp=ztmp(nvrt)-ztmp(kbpl)
                if(total_dp<=0) then
                  write(*,*)'depth<=0:',total_dp,i
                  stop
                endif
                do m=1,ivs
                  sum1=0
                  do k=kbpl,nvrt-1
                    if(i23d<=3) then !node
                      sum1=sum1+(ztmp(k+1)-ztmp(k))*(out2(i,k,m)+out2(i,k+1,m))/2
                    else if(i23d<=6) then !elem
                      sum1=sum1+(ztmp(k+1)-ztmp(k))*outvar(m,k+1,iep(i),irec)
                    endif    
                  enddo !k
                  out3(i,m)=sum1/total_dp
                enddo !m
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
                  do m=1,ivs
                    if(i23d<=3) then !node
                      out3(i,m)=out2(i,k0,m)*(1-rat)+out2(i,k0+1,m)*rat
                    else if(i23d<=6) then !elem
                      if(iep(i)<=0) stop 'iep(i)<=0(2)'
                      out3(i,m)=outvar(m,k0+1,iep(i),irec) !FV
                    endif
                  enddo !m
                endif
              endif !depth average or not
            endif !dry/wet
          enddo !i=1,nxy
          write(18,'(e16.8,20000(1x,f14.6))')timeout(irec_real)/86400,(out3(i,1),i=1,nxy)
          if(ivs==2) write(19,'(e16.8,20000(1x,f14.6))')timeout(irec_real)/86400,(out3(i,2),i=1,nxy)
         
        endif !i23d
!----------------------------------------------------------------------------
        enddo !irec

        if(irec2==nrec) exit loop1
        irec1=irec1+nrec3
      end do loop1
      !enddo !irec=1,nrec
      !iret=nf90_close(ncid)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      enddo !iday

!     Output local depths info
      do i=1,nxy
        write(20,*)i,dep(i)
      enddo !i

      print*, 'Finished!'

      end program read_out
