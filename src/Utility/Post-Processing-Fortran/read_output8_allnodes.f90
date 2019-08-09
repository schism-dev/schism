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

!	Read (combined or uncombined) nc outputs for multiple files at all nodes 
!       Works for mixed tri/quad outputs on NODE based vars only.
!       Inputs: screen; combined or uncombined nc file; vgrid.in (in this dir or ../)
!       Outputs: extract.out (ascii)
!       History: (1) added non-standard outputs (April 2012) - transparent to most scripts
!               as format is same; (2) added ivcor=1 (Dec 2013); (3)
!               added quads (Nov. 2014) (4) changed to nc outputs (Sept
!               2017); (5) added uncombined option (Feb 2019)
!****************************************************************************************
!     ifort -O2 -assume byterecl -o read_output8_allnodes.WW ../UtilLib/extract_mod.f90 read_output8_allnodes.f90 ../UtilLib/compute_zcor.f90 -I$NETCDF/include  -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff
      program read_out
      use netcdf
      use extract_mod

!      parameter(nbyte=4)
      character(len=30) :: file63,varname
      character(len=12) :: it_char
!      character*48 data_format
  
!      integer,allocatable :: i34(:),elnode(:,:)
      allocatable :: sigma(:),cs(:),ztot(:)
      allocatable :: outvar(:,:,:,:),out(:,:,:),icum(:,:,:),eta2(:,:),ztmp(:,:)
      allocatable :: xcj(:),ycj(:),dps(:),kbs(:),xctr(:),yctr(:),dpe(:),kbe(:)
      allocatable :: zs(:,:),ze(:,:),idry(:),outs(:,:,:),oute(:,:,:)
      allocatable :: ic3(:,:),isdel(:,:),isidenode(:,:),kbp(:),sigma_lcl(:,:)
      real*8,allocatable :: timeout(:)

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

      print*, 'Input NODE-based variable name to read from nc (e.g. elev):'
!'
      read(*,'(a30)')varname
      varname=adjustl(varname); len_var=len_trim(varname)
      
      print*, 'Input start and end stack # to read:'
      read(*,*) iday1,iday2

      open(65,file='extract.out')
      
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

      last_dim=np
!     nrec specified if ispec_max=2
      if(ispec_max==1) nrec3=max_array_size/(2*nvrt*last_dim)-1
      nrec3=min(nrec,max(nrec3,1))

      allocate(ztot(nvrt),sigma(nvrt),sigma_lcl(nvrt,np),kbp(np), &
     &outvar(2,nvrt,last_dim,nrec3),eta2(np,nrec3),idry(np),ztmp(nvrt,np),timeout(nrec))
      outvar=-huge(1.0) !test mem

!     Read vgrid.in
      call get_vgrid('vgrid.in',np,nvrt,ivcor,kz,h_s,h_c,theta_b,theta_f,ztot,sigma,sigma_lcl,kbp)

!     Calculate kbp00
      if(ivcor==1) then
        kbp00=kbp
      else
        do i=1,np
          call zcor_SZ(dp(i),1.e8,h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot,sigma, &
     &ztmp(:,i),idry2,kbp00(i))
        enddo !i
      endif !ivcor
      
!...  Time iteration
!...
      outvar=-99 !init.
      ztmp=-99
      do iday=iday1,iday2 !1,ndays
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

!      print*, 'time=',timeout

!      do irec=1,nrec
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

        !Available now:
        !outvar(2,nvrt,np|ne,irec2-irec1+1),i23d,ivs,eta2(np,irec2-irec1+1)

        do irec=1,irec2-irec1+1 !offeset record #
!----------------------------------------------------------------------------
        irec_real=irec1+irec-1 !actual record #

        if(mod(i23d-1,3)==0) then !2D
!         Output: time, 2D variable at all nodes
          write(65,'(e14.6,1000000(1x,e14.4))')time/86400,((outvar(m,1,i,irec),m=1,ivs),i=1,np)
        else !if(i23d==3) then !3D 
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
              endif !wet/dry
            else if(ivcor==2) then !SZ
              call zcor_SZ(dp(i),eta2(i,irec),h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot, &
     &sigma,ztmp(:,i),idry(i),kbpl)
            endif

            do k=max0(1,kbp00(i)),nvrt
!             Output: time, node #, level #, z-coordinate, 3D variable 
              if(idry(i)==1) then
                write(65,*)timeout(irec_real)/86400,i,k,-1.e6,(outvar(m,k,i,irec),m=1,ivs)
              else
                write(65,*)timeout(irec_real)/86400,i,k,ztmp(k,i),(outvar(m,k,i,irec),m=1,ivs)
              endif
            enddo !k

            !Debug
            !write(65,*)x(i),y(i),out(i,1,1:ivs) 

          enddo !i
        endif !i23d
      !enddo !irec=1,nrec
!----------------------------------------------------------------------------
        enddo !irec

        if(irec2==nrec) exit loop1
        irec1=irec1+nrec3
      end do loop1

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      enddo !iday

      print*, 'Finished!'

      stop
      end
