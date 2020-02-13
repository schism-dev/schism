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
!	Read (combined or uncombined) nc outputs and compute average field 
!       at a particular Z-level (use above surface/below bottom to get surface/bottom).
!       Skip dry times for 3D variables.
!       Works for mixed quad/tri on NODE based variables only.
!       Input: schout*.nc (combined or uncombined); vgrid.in; screen
!       Output: average.out (gredit for scalar or xmgr5 format for vector)
!										
!       ifort -O2 -mcmodel=medium -assume byterecl -CB -o compute_average3.exe ../UtilLib/extract_mod.f90 ../UtilLib/compute_zcor.f90 compute_average3.f90 -I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff
!********************************************************************************
!
      program read_out
      use netcdf
      use extract_mod
      use compute_zcor
      parameter(nbyte=4)
      character(len=30) :: file63,varname
      character(len=12) :: it_char
      character(len=48) :: data_format
  
!      integer,allocatable :: i34(:),elnode(:,:)
      integer :: iday1, iday2, iskipst
      allocatable :: sigma(:),cs(:),ztot(:)
      allocatable :: outvar(:,:,:,:),icum(:,:,:),eta2(:,:),ztmp(:,:)
      allocatable :: xcj(:),ycj(:),dps(:),kbs(:),xctr(:),yctr(:),dpe(:),kbe(:)
      allocatable :: idry(:),outs(:,:,:),residual(:,:),icounter(:)
      allocatable :: ic3(:,:),isdel(:,:),isidenode(:,:),kbp(:),sigma_lcl(:,:)
      real*8,allocatable :: timeout(:)

      pi=3.1415926

!...  Set max array size for system memory
      print*, 'Input max array size (e.g., <=2.e9):'
      read(*,*)max_array_size
      print*, 'max_array_size read in=',max_array_size

      print*, 'Do you work on uncombined (0) or combined (1) nc?'
      read(*,*)icomb
      if(icomb/=0.and.icomb/=1) stop 'Unknown icomb'

      print*, 'Input variable name (e.g. salt):'
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

      file63=adjustl(file63)
      len_file63=len_trim(file63)

!...  Header
      !Returned vars: ne,np,ns,nrec,[x y dp](np),
      !elnode,i34,nvrt,h0,dtout
      !If icomb=0, additonal vars:
      !nproc,iegl_rank,iplg,ielg,islg,np_lcl(:),ne_lcl(:),ns_lcl(:)
      if(icomb==0) then !uncombined
        call get_global_geo
      else
        call readheader(iday1)
      endif

      print*, 'After header:',ne,np,nrec,i34(ne), &
     &elnode(1:i34(ne),ne),nvrt,h0,x(np),y(np),dp(np) !,start_time
      
!...  Read in time records in segments for mem
      last_dim=np
      nrec3=max_array_size/(2*nvrt*last_dim)-1
      nrec3=min(nrec,max(nrec3,1))
      print*, '# of records in each stack is: ',nrec, &
     &'; actual # of records that can be read in each time is:',nrec3

      allocate(timeout(nrec),ztot(nvrt),sigma(nvrt),sigma_lcl(nvrt,np),kbp(np), &
     &outvar(2,nvrt,last_dim,nrec3),eta2(np,nrec3),ztmp(nvrt,np),residual(np,2),icounter(np),idry(np))
      outvar=-huge(1.0) !test mem

      call get_vgrid_single('vgrid.in',np,nvrt,ivcor,kz,h_s,h_c,theta_b,theta_f,ztot,sigma,sigma_lcl,kbp)

      print*, 'last element',elnode(1:i34(ne),ne)

!...  Time iteration
!...
      out=-99 !init.
      ztmp=-99
      residual=0
!      ndays=iday2-iday1+1
      icounter=0 !counter for each node (wet/dry)
      do iday=iday1,iday2,iskipst
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!      write(it_char,'(i12)')iday
!      it_char=adjustl(it_char)
!      leng=len_trim(it_char)
!      file63='schout_'//it_char(1:leng)//'.nc'
!      iret=nf90_open(trim(adjustl(file63)),NF90_NOWRITE,ncid)
!      !time is double
!      iret=nf90_inq_varid(ncid,'time',itime_id)
!      iret=nf90_get_var(ncid,itime_id,timeout,(/1/),(/nrec/))
      if(icomb==0) then !uncombined
        call get_timeout(iday,nrec,timeout,icomb)
      else
        call get_timeout(iday,nrec,timeout)
      endif
!      print*, 'time=',timeout

      irec1=1 !start record
      loop1: do
        irec2=min(nrec,irec1+nrec3-1)

        if(icomb==0) then !uncombined
          do irank=0,nproc-1
            call get_outvar_multirecord(1,iday,varname,irec1,irec2,np,last_dim,nvrt,nrec3,outvar,i23d,ivs,eta2,irank)
          enddo !irank
        else
          call get_outvar_multirecord(1,iday,varname,irec1,irec2,np,last_dim,nvrt,nrec3,outvar,i23d,ivs,eta2)
        endif
!        print*, 'done reading records:',irec1,irec2

        !Available now: outvar(2,nvrt,np|ne,irec2-irec1+1),i23d,ivs,eta2(np,irec2-irec1+1)

        do irec=1,irec2-irec1+1 !offeset record #
!----------------------------------------------------------------------------
        irec_real=irec1+irec-1 !actual record #
        if(mod(i23d-1,3)==0) then !2D
          !print*,'irec=',irec,iday1,irec_start,iday2,irec_end
          if(.not.(iday==iday1.and.irec_real<irec_start.or.iday==iday2.and.irec_real>irec_end)) then
            do i=1,np
              icounter(i)=icounter(i)+1
              residual(i,1:ivs)=residual(i,1:ivs)+outvar(1:ivs,1,i,irec)
            enddo !i
          endif
        else !3D
          !Compute z coordinates
          do i=1,np
            if(ivcor==1) then !localized
              if(dp(i)+eta2(i,irec)<=h0) then
                idry(i)=1
              else !wet
                idry(i)=0
                do k=kbp(i),nvrt
                  ztmp(k,i)=(eta2(i,irec)+dp(i))*sigma_lcl(k,i)+eta2(i,irec)
                enddo !k
              endif !dp
            else if(ivcor==2) then !SZ
              call zcor_SZ_single(dp(i),eta2(i,irec),h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot, &
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
              if(.not.(iday==iday1.and.irec_real<irec_start.or.iday==iday2.and.irec_real>irec_end)) then
                icounter(i)=icounter(i)+1
                do m=1,ivs
                  tmp=outvar(m,k0,i,irec)*(1-rat)+outvar(m,k0+1,i,irec)*rat
                  residual(i,m)=residual(i,m)+tmp
                enddo !m
              endif !not
            endif !idry
          enddo !i=1,np
        endif !2/3D
!----------------------------------------------------------------------------
        enddo !irec

        if(irec2==nrec) exit loop1
        irec1=irec1+nrec3
      end do loop1

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      enddo !iday

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
      if(varname(1:len_var).eq.'WWM_18') then
        !Read in sub-samples
        open(20,file='nodeflags.bp',status='old')
        read(20,*); read(20,*)

        eta2(:,1)=residual(:,1) !temp save
        do i=1,np
          read(20,*)j,tmp,tmp,tmp2
          residual(i,1)=-sin(eta2(i,1)/180*pi)
          residual(i,2)=-cos(eta2(i,1)/180*pi)
          if(nint(tmp2)==0) residual(i,1:2)=0
        enddo !i
        ivs=2 !xyuv format
      endif !file63

      if(ivs==1) then
        write(65,*)iday1,iday2
        write(65,*)ne,np
        do i=1,np
          write(65,*)i,x(i),y(i),residual(i,1)
        enddo !i
        do i=1,ne
          write(65,*)i,i34(i),elnode(1:i34(i),i)
        enddo !i
      else !vectors
        do i=1,np
          write(65,*)x(i),y(i),residual(i,1:2)
        enddo !i
      endif
      close(65)

      print*, 'Finished!'

      stop
      end
