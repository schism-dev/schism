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
!       Inputs: (1) nc files and out2d*.nc;
!               (2) station.xyzt: make sure all times are after 1st record (to ensure interpolation in time); 
!                                 pad extra days before and after if necessary.
!                                 z>=0 from surface.
!               (3) screen inputs
!               (4) vgrid.in: in this dir or ../ 
!       Outputs: fort.1[89]; fort.11 (fatal errors); fort.12: nonfatal errors.
!											
! ifort -mcmodel=medium -CB -O2 -o read_output10_xyzt.exe ../UtilLib/extract_mod2.f90 ../UtilLib/compute_zcor.f90 ../UtilLib/pt_in_poly_test.f90 read_output10_xyzt.f90 -I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff
!****************************************************************************************
!
      program read_out
      use netcdf
      use extract_mod2
      use compute_zcor
      use pt_in_poly_test

      character(len=30) :: file63,varname,file62,varname2,file64
      character(len=12) :: it_char
      logical :: lexist,lexist2
      integer,allocatable :: kbp(:),kbp00(:),iday(:,:),irecord(:,:),node3(:,:),iep(:),irank_read(:), &
     &i34(:),elnode(:,:)
      real,allocatable ::  sigma(:),cs(:),ztot(:),times(:,:),out(:,:,:),out2(:,:,:), &
    &eta2(:),arco(:,:),ztmp(:),x00(:),y00(:),z00(:),t00(:),sigma_lcl(:,:),ztmp2(:,:),dp(:), &
    &outvar(:,:,:),out3(:,:),out4(:)
      real*8,allocatable :: timeout(:),xnd(:),ynd(:)
      integer :: char_len,start_2d(2),start_3d(3),start_4d(4), &
     &count_2d(2),count_3d(3),count_4d(4),dimids(100),idims(100)

      integer :: nodel(3) !,dimids(100),idims(100), &
!     &start_2d(2),start_3d(3),start_4d(4), &
!     &count_2d(2),count_3d(3),count_4d(4)

      dimension swild(3)
      !long int for large files
      !integer(kind=8) :: irec
      
      print*, 'Input variable name to read from (e.g. elevation):'
      print*, '(if vector, input X component)'
      read(*,'(a30)')varname
      varname=adjustl(varname); len_var=len_trim(varname)
      
!     Invliad number used for 3D variables: below bottom; dry spot; no parents
      print*, 'Input values to be used for invalid place:'
      read(*,*)rjunk

      print*, 'Variable is node (1) or elem (2) based?'
      read(*,*)inode_elem

      print*, 'Variable is scalar (1) or vector (2)?'
      read(*,*)ivs
      if(ivs/=1.and.ivs/=2) stop 'check ivs'

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

!...  Get basic info from out2d*.nc
      !Returned vars: ne,np,ns,nrec,[xnd ynd dp](np),
      !elnode,i34,nvrt,h0,dtout,kbp
      call get_dims(1,np,ne,ns,nvrt,h0)
      allocate(xnd(np),ynd(np),dp(np),kbp(np),i34(ne),elnode(4,ne),stat=istat)
      if(istat/=0) stop 'alloc (1)'
      call readheader(1,np,ne,ns,kbp,i34,elnode,nrec,xnd,ynd,dp,dtout)

      print*, 'After header:',ne,np,ns,nvrt,nrec,i34(ne), &
     &elnode(1:i34(ne),ne),h0,xnd(np),ynd(np),dp(np),dtout !,start_time

      last_dim=max(np,ne,ns)
      allocate(timeout(nrec),out(3,nvrt,2),out2(2,nvrt,ivs),eta2(np),node3(nxy,3),arco(nxy,3), &
    &iep(nxy),iday(2,nxy),irecord(2,nxy),times(2,nxy),outvar(nvrt,last_dim,ivs),out3(2,ivs),out4(ivs),stat=istat)
      outvar=-huge(1.0)
      if(istat/=0) stop 'Failed to allocate (3)'

!     Read in vgrid.in 
      allocate(ztot(nvrt),sigma(nvrt),sigma_lcl(nvrt,np),kbp00(np))
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
          write(11,*)'Cannot find a parent for pt:',j,x00(j),y00(j)
          iabort=1
        endif
      enddo !j
      if(iabort==1) stop 'check fort.11 for pts outside'

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
          print*, 'reading stack ',iday(l,i),' for point ',i
          write(it_char,'(i12)')iday(l,i)
          it_char=adjustl(it_char)
          leng=len_trim(it_char)
          file62='out2d_'//it_char(1:leng)//'.nc'
          iret=nf90_open(trim(adjustl(file62)),OR(NF90_NETCDF4,NF90_NOWRITE),ncid4)
          !time is double
          iret=nf90_inq_varid(ncid4,'time',itime_id)
          iret=nf90_get_var(ncid4,itime_id,timeout,(/1/),(/nrec/))
!          print*, 'time=',timeout !,trim(adjustl(file63))
     
          !Find nc file
          file63=varname(1:len_var)//'_'//it_char(1:leng)//'.nc'
          if(ivs==2) then
            if(varname(len_var:len_var).ne."X") stop 'Vector name must end with X'
            varname2=varname(1:len_var-1)//'Y'
          endif
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
            endif !ivs
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
          do ii=1,ndims
           iret=nf90_Inquire_Dimension(ncid,dimids(ii),len=idims(ii))
          enddo !ii
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

          irec=irecord(l,i)

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
    
          !Debug
!          write(98,*)i,l
!          do j=1,np 
!            write(98,*)j,kbp00(j),outvar(1,:,j)
!          enddo !j
!          stop

          out2(l,:,:)=0
          out3(l,:)=0
          if(i23d==1) then !2D
            do j=1,3 !nodes
              nd=node3(i,j)
              if(inode_elem==1) then !node based
                out2(l,1,:)=out2(l,1,:)+arco(i,j)*outvar(1,nd,:)
              else !elem
                out2(l,1,:)=outvar(1,iep(i),:)
              endif
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
                do j=1,3
                  nd=node3(i,j)
                  kin=max(k,kbp00(nd))
                  if(inode_elem==1) then !node
                    out2(l,k,:)=out2(l,k,:)+arco(i,j)*outvar(kin,nd,:)
                  else !elem
                    out2(l,k,:)=outvar(k,iep(i),:) 
                  endif
                enddo !j
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
                if(inode_elem==1) then !node
                  out3(l,:)=out2(l,k0,:)*(1-rat)+out2(l,k0+1,:)*rat
                else if (i23d<=6) then !elem
                  out3(l,:)=out2(l,k0+1,:) !FV
                endif
              endif
            endif !dry/wet
          endif !i23d

          iret=nf90_close(ncid)
          iret=nf90_close(ncid4)
          if(ivs==2) iret=nf90_close(ncid2)

        enddo loop1 !l=1,2; 2 times

!       Interpolate in time
        trat=(t00(i)-times(1,i))/(times(2,i)-times(1,i)) !must be [0,1]
        if(i23d==1) then !2D
          if(iep(i)==0) then !no parents
            out4(:)=rjunk
          else
            out4(:)=out2(1,1,:)*(1-trat)+out2(2,1,:)*trat
          endif
          write(18,'(e16.8,1000(1x,f15.3))')t00(i)/86400,out4(1),x00(i),y00(i)
          if(ivs==2) write(19,'(e16.8,1000(1x,f15.3))')t00(i)/86400,out4(2),x00(i),y00(i)
        else !3D
          if(iep(i)==0) then !no parents
            out4(:)=rjunk
          else
            out4(:)=out3(1,:)*(1-trat)+out3(2,:)*trat
          endif
          write(18,*)t00(i)/86400,out4(1),x00(i),y00(i),z00(i)
          if(ivs==2) write(19,*)t00(i)/86400,out4(2),x00(i),y00(i),z00(i)
        endif
      enddo !i=1,nxy

      print*, 'Finished!'

      stop
      end

