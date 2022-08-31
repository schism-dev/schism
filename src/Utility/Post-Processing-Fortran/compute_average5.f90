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
!********************************************************************************
!										
!	Read nc outputs and compute average field from scribed I/O
!       at a particular Z-level (use above surface/below bottom to get surface/bottom).
!       Skip dry times for 3D variables.
!       Works for mixed quad/tri on NODE based variables only.
!       Input: out2d* and corresponding nc file if the variable is 3D; vgrid.in; screen
!       Output: average.out (gredit fromat)
!										
!       ifort -O2 -mcmodel=medium -assume byterecl -CB -o compute_average5.exe ../UtilLib/extract_mod2.f90 ../UtilLib/compute_zcor.f90 compute_average5.f90 -I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff
!********************************************************************************
!
      program read_out
      use netcdf
      use extract_mod2
      use compute_zcor
      parameter(nbyte=4)
      character(len=30) :: file62,file63,varname
      character(len=12) :: it_char
      character(len=48) :: data_format
  
      integer,allocatable :: i34(:),elnode(:,:)
      integer :: iday1, iday2, iskipst
      allocatable :: sigma(:),cs(:),ztot(:)
      allocatable :: outvar(:,:),icum(:,:,:),eta2(:),ztmp(:,:)
      allocatable :: xcj(:),ycj(:),dp(:),dps(:),kbs(:),xctr(:),yctr(:),dpe(:),kbe(:)
      allocatable :: idry(:),outs(:,:,:),residual(:,:),icounter(:)
      allocatable :: ic3(:,:),isdel(:,:),isidenode(:,:),kbp(:),sigma_lcl(:,:)
      real*8,allocatable :: timeout(:),xnd(:),ynd(:)
      integer :: char_len,start_2d(2),start_3d(3),start_4d(4), &
     &count_2d(2),count_3d(3),count_4d(4),dimids(100),idims(100)

      pi=3.1415926
      ivs=1 !for vectors, use own script to combine components

      print*, 'Input variable name (e.g. salinity):'
      read(*,'(a30)')varname
      varname=adjustl(varname); len_var=len_trim(varname)
      
      print*, 'Input start and end stack #s to read:'
      read(*,*) iday1,iday2

      print*, 'Input stride in stacks (1 to include all):'
      read(*,*) iskipst

      print*, 'Input start and end record #s in start|end stack respectively:'
!'
      read(*,*) irec_start,irec_end

      print*, 'Input z-coord. (<=0 below MSL):'
      read(*,*) z00

!      if(mod(iday2-iday1,iskipst)/=0) then
!        write(*,*)'should be n skips over stack1 and stack2:',iday1,iday2,iskipst
!        stop
!      endif

!      file63=adjustl(file63)
!      len_file63=len_trim(file63)

!...  Header
!...  Get basic info from out2d*.nc
      !Returned vars: ne,np,ns,nrec,[xnd ynd dp](np),
      !elnode,i34,nvrt,h0,dtout,kbp
      call get_dims(iday1,np,ne,ns,nvrt,h0)
      allocate(xnd(np),ynd(np),dp(np),kbp(np),i34(ne),elnode(4,ne),stat=istat)
      if(istat/=0) stop 'alloc (1)'
      call readheader(iday1,np,ne,ns,kbp,i34,elnode,nrec,xnd,ynd,dp,dtout)

      print*, 'After header:',ne,np,nrec,i34(ne),elnode(1:i34(ne),ne),nvrt,h0,xnd(np),ynd(np),dp(np) !,start_time

!...  Read in time records in segments for mem
      last_dim=np
      allocate(ztot(nvrt),sigma(nvrt),sigma_lcl(nvrt,np),timeout(nrec), &
     &outvar(nvrt,last_dim),eta2(np),ztmp(nvrt,np),residual(np,2),icounter(np),idry(np))
      outvar=-huge(1.0) !test mem

      call get_vgrid_single('vgrid.in',np,nvrt,ivcor,kz,h_s,h_c,theta_b,theta_f,ztot,sigma,sigma_lcl,kbp)

!...  Time iteration
!...
      outvar=-99
      ztmp=-99
      residual=0
      icounter=0 !counter for each node (wet/dry)
      do iday=iday1,iday2,iskipst
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
      if(npes/=np) stop 'can only handle node-based'
!'
      if(idims(ndims)/=nrec) stop 'last dim is not time'

      do irec=1,nrec
!----------------------------------------------------------------------------
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

        if(i23d==1) then !2D
          !print*,'irec=',irec,iday1,irec_start,iday2,irec_end
          if(.not.(iday==iday1.and.irec<irec_start.or.iday==iday2.and.irec>irec_end)) then
            do i=1,np
              icounter(i)=icounter(i)+1
              residual(i,1)=residual(i,1)+outvar(1,i)
            enddo !i
          endif
        else !3D
          !Compute z coordinates
          do i=1,np
            if(ivcor==1) then !localized
              if(dp(i)+eta2(i)<=h0) then
                idry(i)=1
              else !wet
                idry(i)=0
                do k=kbp(i),nvrt
                  ztmp(k,i)=(eta2(i)+dp(i))*sigma_lcl(k,i)+eta2(i)
                enddo !k
              endif !dp
            else if(ivcor==2) then !SZ
              call zcor_SZ_single(dp(i),eta2(i),h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot, &
     &sigma,ztmp(:,i),idry(i),kbpl)
              kbp(i)=kbpl
            endif !ivcor

            if(idry(i)==0) then !wet
              !Interplate in vertical
              if(z00>=ztmp(nvrt,i)) then !above F.S.
                k0=nvrt-1; rat=1
              else if(z00<=ztmp(kbp(i),i)) then !below bottom; extrapolate
                k0=kbp(i); rat=0
              else !above bottom; cannot be above F.S.
                k0=0
                do k=kbp(i),nvrt-1
                  if(z00>=ztmp(k,i).and.z00<=ztmp(k+1,i)) then
                    k0=k
                    rat=(z00-ztmp(k,i))/(ztmp(k+1,i)-ztmp(k,i))
                    exit
                  endif
                enddo !k
              endif !ztmp

              if(k0==0) then
                write(*,*)'failed to find a vertical level:',irec_real,i,ifs,z2,ztmp(:,i)
                stop
              endif
              if(.not.(iday==iday1.and.irec<irec_start.or.iday==iday2.and.irec>irec_end)) then
                icounter(i)=icounter(i)+1
                do m=1,ivs
                  tmp=outvar(k0,i)*(1-rat)+outvar(k0+1,i)*rat
                  residual(i,m)=residual(i,m)+tmp
                enddo !m
              endif !not
            endif !idry
          enddo !i=1,np
        endif !2/3D
!----------------------------------------------------------------------------
      enddo !irec=1,nrec
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      enddo !iday=iday1,iday2,iskipst

!     Output
      do i=1,np
        if(icounter(i)==0) then
          residual(i,:)=0
        else
          residual(i,:)=residual(i,:)/icounter(i)
        endif !icounter
      enddo !i

      open(65,file='average.out')

      !For wave dir
!      if(varname(1:len_var).eq.'meanWaveDirection') then
!        !Read in sub-samples
!        open(20,file='nodeflags.bp',status='old')
!        read(20,*); read(20,*)
!
!        eta2(:)=residual(:,1) !temp save
!        do i=1,np
!          read(20,*)j,tmp,tmp,tmp2
!          residual(i,1)=-sin(eta2(i)/180*pi)
!          residual(i,2)=-cos(eta2(i)/180*pi)
!          if(nint(tmp2)==0) residual(i,1:2)=0
!        enddo !i
!        ivs=2 !xyuv format
!      endif !file63

      if(ivs==1) then
        write(65,*)iday1,iday2
        write(65,*)ne,np
        do i=1,np
          write(65,*)i,xnd(i),ynd(i),residual(i,1)
        enddo !i
        do i=1,ne
          write(65,*)i,i34(i),elnode(1:i34(i),i)
        enddo !i
!      else !vectors
!        do i=1,np
!          write(65,*)x(i),y(i),residual(i,1:2)
!        enddo !i
      endif !ivs
      close(65)

      print*, 'Finished!'

      stop
      end
