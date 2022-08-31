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
!	Read in (x,y,z,time) from station.xyzt (z>=0 is distance from F.S.; time in sec) 
!         for 3D variables (surface values for 2D variables) DEFINED @ nodes or elem. Interpolation in time.
!         Not working for lon/lat.
!         Works for mixed tri/quad outputs, combined or uncombined nc outputs.
!       Inputs: (1) nc files;
!               (2) station.xyzt: make sure all times are after 1st record (to ensure interpolation in time); 
!                                 pad extra days before and after if necessary.
!                                 z>=0 from surface.
!               (3) screen inputs: varname; invalid value (for out of domain, dry etc)
!               (4) vgrid.in: in this dir or ../ 
!       Outputs: fort.1[89]; fort.11 (fatal errors); fort.12: nonfatal errors.
!											
! ifort -mcmodel=medium -CB -O2 -o read_output9_xyzt.exe ../UtilLib/extract_mod.f90 \
! ../UtilLib/compute_zcor.f90 ../UtilLib/pt_in_poly_test.f90 read_output9_xyzt.f90 \
!-I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff
!****************************************************************************************
!
      program read_out
      use netcdf
      use extract_mod
      use compute_zcor
      use pt_in_poly_test

      character(len=30) :: file63,varname
      character(len=12) :: it_char
      integer,allocatable :: kbp(:),iday(:,:),irecord(:,:),node3(:,:),iep(:),irank_read(:)
      real,allocatable ::  sigma(:),cs(:),ztot(:),times(:,:),out(:,:,:),out2(:,:,:), &
    &eta2(:),arco(:,:),ztmp(:),x00(:),y00(:),z00(:),t00(:),sigma_lcl(:,:),ztmp2(:,:), &
    &outvar(:,:,:)
      real*8,allocatable :: timeout(:)

      integer :: nodel(3) !,dimids(100),idims(100), &
!     &start_2d(2),start_3d(3),start_4d(4), &
!     &count_2d(2),count_3d(3),count_4d(4)

      dimension swild(3),out3(2,2),out4(2)
      !long int for large files
      !integer(kind=8) :: irec
      
      print*, 'Do you work on uncombined (0) or combined (1) nc?'
      read(*,*)icomb
      if(icomb/=0.and.icomb/=1) stop 'Unknown icomb'

      print*, 'Input variable name to read from (e.g. elev):'
      read(*,'(a30)')varname
      varname=adjustl(varname); len_var=len_trim(varname)
      
!     Invliad number used for 3D variables: below bottom; dry spot; no parents
      print*, 'Input values to be used for invalid place:'
      read(*,*)rjunk

      open(10,file='station.xyzt',status='old')
      read(10,*) 
      read(10,*) nxy
      allocate(x00(nxy),y00(nxy),z00(nxy),t00(nxy),stat=istat)
      if(istat/=0) stop 'Falied to allocate (1)'

      do i=1,nxy
        read(10,*)j,x00(i),y00(i),z00(i),t00(i) !z00>=0 from F.S.; t00 in sec
        if(z00(i)<0) then
          write(*,*)'Invalid z value:',i; stop
        endif
      enddo !i
      close(10)

!...  Header
      !Returned vars: ne,np,ns,nrec,[x y dp](np),
      !elnode,i34,nvrt,h0,dtout
      !If icomb=0, additonal vars:
      !nproc,iegl_rank,iplg,ielg,islg,np_lcl(:),ne_lcl(:),ns_lcl(:)
      if(icomb==0) then !uncombined
        call get_global_geo
      else
        call readheader(1)
      endif

      print*, 'After header:',ne,np,ns,nrec,i34(ne), &
     &elnode(1:i34(ne),ne),nvrt,h0,x(np),y(np),dp(np) !,start_time

      last_dim=max(np,ne,ns)
      allocate(timeout(nrec),out(3,nvrt,2),out2(2,nvrt,2),eta2(np),node3(nxy,3),arco(nxy,3), &
    &iep(nxy),iday(2,nxy),irecord(2,nxy),times(2,nxy),outvar(2,nvrt,last_dim),stat=istat)
      outvar=-huge(1.0)
      if(istat/=0) stop 'Failed to allocate (3)'

!     Read in vgrid.in 
      allocate(ztot(nvrt),sigma(nvrt),sigma_lcl(nvrt,np),kbp(np))
      call get_vgrid_single('vgrid.in',np,nvrt,ivcor,kz,h_s,h_c,theta_b,theta_f,ztot,sigma,sigma_lcl,kbp)

      allocate(ztmp(nvrt),ztmp2(nvrt,3))

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
!      arco=1./3 !initialize for pts without parents
!      do l=1,nxy
!        node3(l,1:3)=elnode(1:3,1) !initialize for pts without parents
!      enddo !l

      do i=1,ne
        do l=1,nxy
          if(iep(l)/=0) cycle

          call pt_in_poly_single(i34(i),real(x(elnode(1:i34(i),i))), &
     &real(y(elnode(1:i34(i),i))),x00(l),y00(l),inside,arco(l,1:3),nodel)
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
          write(11,*)'Cannot find a parent for pt:',j,x00(j),y00(j)
          iabort=1
        endif
      enddo !j
      if(iabort==1) stop 'check fort.11 for pts outside'

      !Mark ranks that need to be read in
      if(icomb==0) then
        allocate(irank_read(0:nproc-1))
        irank_read=0
        do j=1,nxy
          irank_read(iegl_rank(iep(j)))=1
          write(99,*)'reading from rank #:',iegl_rank(iep(j))
        enddo
      endif !icomb

!...  Compute stack and record # for each pt
      do i=1,nxy
!       Check if time is before first record
        if(t00(i)<dtout) then
          write(11,*)'Time before first record:',i,t00(i)
          stop
        endif

!       Lower and upper bound stacks and record #s for t00
        iday(1,i)=(t00(i)-dtout)/nrec/dtout+1
        if(iday(1,i)<1) then
          write(11,*)'Impossible'; stop
        else
          irecord(1,i)=(t00(i)-(iday(1,i)-1)*nrec*dtout)/dtout
          !Bounding record time just b4 the cast time, corresponding to record
          !irecord(1,i) in stack iday(1,i)
          times(1,i)=((iday(1,i)-1)*nrec+irecord(1,i))*dtout
          iday(2,i)=t00(i)/nrec/dtout+1
          irecord(2,i)=(t00(i)-(iday(2,i)-1)*nrec*dtout)/dtout+1
          !Bounding record time just after the cast time, corresponding to
          !record irecord(2,i) in stack iday(2,i). Note that irecord(2,i)
          !may<irecord(1,i) (e.g. t00 before 1st record of iday(2,i))
          times(2,i)=((iday(2,i)-1)*nrec+irecord(2,i))*dtout
        endif

        if(irecord(1,i)>nrec.or.irecord(2,i)>nrec) then
          write(11,*)'Record # overflow: ',i,irecord(:,i)
          stop
        endif
        if(t00(i)<times(1,i).or.t00(i)>times(2,i)) then
          write(11,*)'Wrong time bounds:',i,t00(i),times(:,i),iday(:,i),irecord(:,i)
          stop
        endif
      enddo !i=1,nxy

!...  Time iteration
      do i=1,nxy
        loop1: do l=1,2 !2 times
!          write(it_char,'(i12)')iday(l,i)
!          it_char=adjustl(it_char); leng=len_trim(it_char)
!          file63='schout_'//it_char(1:leng)//'.nc'
!          iret=nf90_open(trim(adjustl(file63)),OR(NF90_NETCDF4,NF90_NOWRITE),ncid)
!          !time is double
!          iret=nf90_inq_varid(ncid,'time',itime_id)
!          iret=nf90_get_var(ncid,itime_id,timeout,(/1/),(/nrec/))
          if(icomb==0) then !uncombined
            call get_timeout(iday(l,i),nrec,timeout,icomb)
          else
            call get_timeout(iday(l,i),nrec,timeout)
          endif

!           print*, 'time=',timeout,trim(adjustl(file63))

          irec=irecord(l,i)
          if(icomb==0) then !uncombined
            do irank=0,nproc-1
              if(irank_read(irank)>0) then
                call get_outvar(iday(l,i),irec,varname,np,last_dim,nvrt,outvar,i23d,ivs,eta2,irank)
              endif
            enddo !irank
          else
            call get_outvar(iday(l,i),irec,varname,np,last_dim,nvrt,outvar,i23d,ivs,eta2)
          endif

          !Get elev
!          start_2d(1)=1; start_2d(2)=irec
!          count_2d(1)=np; count_2d(2)=1
!          iret=nf90_get_var(ncid,ielev_id,eta2,start_2d,count_2d)
!
!          iret=nf90_inq_varid(ncid,varname(1:len_var),ivarid1)
!          if(iret/=nf90_NoErr) stop 'Var not found'
!          iret=nf90_Inquire_Variable(ncid,ivarid1,ndims=ndims,dimids=dimids)
!          if(ndims>100) stop 'increase dimension of dimids & idims'
!          do ii=1,ndims
!            iret=nf90_Inquire_Dimension(ncid,dimids(ii),len=idims(ii))
!          enddo !i
!          npes=idims(ndims-1) !np|ne|ns
!          if(npes/=np.and.npes/=ne) stop 'can only handle node- or elem-based'
!          if(idims(ndims)/=nrec) stop 'last dim is not time'
!
!          iret=nf90_get_att(ncid,ivarid1,'i23d',i23d)
!          if(i23d<=0.or.i23d>6) stop 'wrong i23d'
!          iret=nf90_get_att(ncid,ivarid1,'ivs',ivs)
!          !print*, 'i23d:',i23d,ivs,idims(1:ndims)
!          if(ivs==1) then !scalar
!            if(mod(i23d-1,3)==0) then !2D
!              start_2d(1)=1; start_2d(2)=irec
!              count_2d(1)=npes; count_2d(2)=1
!              iret=nf90_get_var(ncid,ivarid1,outvar(1,1,1:npes),start_2d,count_2d)
!            else !3D
!              start_3d(1:2)=1; start_3d(3)=irec
!              count_3d(2)=npes; count_3d(1)=nvrt; count_3d(3)=1
!              iret=nf90_get_var(ncid,ivarid1,outvar(1,:,1:npes),start_3d,count_3d)
!            endif
!          else !vector
!            if(mod(i23d-1,3)==0) then !2D
!              start_3d(1:2)=1; start_3d(3)=irec
!              count_3d(2)=npes; count_3d(1)=2; count_3d(3)=1
!              iret=nf90_get_var(ncid,ivarid1,outvar(1:2,1,1:npes),start_3d,count_3d)
!            else if(ndims-1==3) then !3D vector
!              start_4d(1:3)=1; start_4d(4)=irec
!              count_4d(3)=npes; count_4d(2)=nvrt; count_4d(1)=2; count_4d(4)=1
!              iret=nf90_get_var(ncid,ivarid1,outvar(:,:,1:npes),start_4d,count_4d)
!            else
!              stop 'Unknown type(2)'
!            endif
!          endif !ivs

          !Available now: outvar(2,nvrt,np|ne),eta2(np)
          !However, for uncombined nc, values in untouched ranks are junk
    
          !Debug
!          write(98,*)i,l
!          do j=1,np 
!            write(98,*)j,kbp00(j),outvar(1,:,j)
!          enddo !j
!          stop

          out2(l,:,:)=0
          out3(l,:)=0
          if(mod(i23d-1,3)==0) then !2D
            do j=1,3 !nodes
              nd=node3(i,j)
              do m=1,ivs
                if(i23d<=3) then !node
                  out2(l,1,m)=out2(l,1,m)+arco(i,j)*outvar(m,1,nd)
                else if (i23d<=6) then !elem
                  out2(l,1,m)=outvar(m,1,iep(i))
                endif
              enddo !m
            enddo !j
          else !3D
!           Do interpolation
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
              out3(:,:)=rjunk
              exit loop1
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
       
              do k=kbpl,nvrt
                do m=1,ivs
                  do j=1,3
                    nd=node3(i,j)
                    kin=max(k,kbp00(nd))
                    if(i23d<=3) then !node
                      out2(l,k,m)=out2(l,k,m)+arco(i,j)*outvar(m,kin,nd)
                    else if (i23d<=6) then !elem
                      out2(l,k,m)=outvar(m,k,iep(i)) 
                    endif
                  enddo !j
                enddo !m
              enddo !k

!             Interplate in vertical
              k0=0
              do k=kbpl,nvrt-1
                if(ztmp(nvrt)-z00(i)>=ztmp(k).and.ztmp(nvrt)-z00(i)<=ztmp(k+1)) then
                  k0=k
                  rat=(ztmp(nvrt)-z00(i)-ztmp(k))/(ztmp(k+1)-ztmp(k))
                  exit
                endif
              enddo !k
              if(k0==0) then
                out3(:,:)=rjunk
                exit loop1
!               write(12,*)'Warning: failed to find a vertical level:',it,i
              else
                do m=1,ivs
                  if(i23d<=3) then !node
                    out3(l,m)=out2(l,k0,m)*(1-rat)+out2(l,k0+1,m)*rat
                  else if (i23d<=6) then !elem
                    out3(l,m)=out2(l,k0+1,m) !FV
                  endif
                enddo !m
              endif
            endif !dry/wet
          endif !i23d
        enddo loop1 !l=1,2; 2 times

!       Interpolate in time
        trat=(t00(i)-times(1,i))/(times(2,i)-times(1,i)) !must be [0,1]
        if(mod(i23d-1,3)==0) then
          if(iep(i)==0) then !no parents
            out4(1:ivs)=rjunk
          else
            out4(1:ivs)=out2(1,1,1:ivs)*(1-trat)+out2(2,1,1:ivs)*trat
          endif
          write(18,'(e16.8,1000(1x,f15.3))')t00(i)/86400,out4(1:ivs),x00(i),y00(i)
        else !3D
          if(iep(i)==0) then !no parents
            out4(1:ivs)=rjunk
          else
            out4(1:ivs)=out3(1,1:ivs)*(1-trat)+out3(2,1:ivs)*trat
          endif
          write(18,*)t00(i)/86400,out4(1),x00(i),y00(i),z00(i)
          !write(18,'(e16.8,1000(1x,f15.3))')t00(i)/86400,out4(1),x00(i),y00(i),z00(i)
          if(ivs==2) write(19,'(e16.8,1000(1x,f15.3))')t00(i)/86400,out4(2),x00(i),y00(i),z00(i)
        endif
      enddo !i=1,nxy

      print*, 'Finished!'

      stop
      end

