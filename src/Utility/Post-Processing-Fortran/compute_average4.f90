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
!       This script also can be used as a simple template for other analyses, as
!       it shows how to read in and make sense of outputs.

!	Read (combined or uncombined) nc outputs and compute average field.
!       Skip dry times for 3D variables.
!       Works for mixed quad/tri on node/elem/side based variables
!       @whole levels (for half levels, change how the depth
!       integration is done in the code)
!       Input: schout*.nc (combined or uncombined); vgrid.in; screen
!       Output: average.out (gredit for scalar or xmgr5 format for vector for nodes)
!										
!       ifort -O2 -mcmodel=medium -assume byterecl -CB -o compute_average4.exe ../UtilLib/extract_mod.f90 ../UtilLib/schism_geometry.f90 ../UtilLib/compute_zcor.f90 compute_average4.f90 -I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff
!********************************************************************************
!
      program read_out
      use netcdf
      use extract_mod
      use compute_zcor
      use schism_geometry_mod
      parameter(nbyte=4)
      character(len=30) :: file63,varname
      character(len=12) :: it_char
      character(len=48) :: data_format
  
!      integer,allocatable :: i34(:),elnode(:,:)
      integer :: iday1, iday2, iskipst
      allocatable :: sigma(:),cs(:),ztot(:)
      allocatable :: outvar(:,:,:),icum(:,:,:),eta2(:),ztmp(:,:)
      allocatable :: xcj(:),ycj(:),dps(:),kbs(:),xctr(:),yctr(:),dpe(:),kbe(:)
      allocatable :: idry(:),outs(:,:,:),residual(:,:),icounter(:)
      allocatable :: ic3(:,:),isdel(:,:),isidenode(:,:),kbp(:),sigma_lcl(:,:)
      integer, allocatable :: elside(:,:)
      real :: wild(100)
      real*8,allocatable :: timeout(:)

      pi=3.1415926

!...  Set max array size for system memory
      print*, 'Do you work on uncombined (0) or combined (1) nc?'
      read(*,*)icomb
      if(icomb/=0.and.icomb/=1) stop 'Unknown icomb'

      print*, 'Input variable name (e.g. salt):'
      read(*,'(a30)')varname
      varname=adjustl(varname); len_var=len_trim(varname)
      
      print*, 'Is the var node(1)/elem(2)/side(3) based?'
      read(*,*)ibase

      print*, 'Input start and end stack #s to read:'
      read(*,*) iday1,iday2

!      print*, 'Input stride in stacks (1 to include all):'
!      read(*,*) iskipst
!      print*, 'Input start and end record #s in start|end stack respectively:'
!!'
!      read(*,*) irec_start,irec_end

!      print*, 'Input z-coord. (<=0 below MSL; not used at the moment):'
!      read(*,*) z00

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

      print*, 'After header:',ne,np,ns,nrec,i34(ne), &
     &elnode(1:i34(ne),ne),nvrt,h0,x(np),y(np),dp(np) !,start_time

!     Geometry
      call compute_nside(np,ne,i34,elnode(1:4,1:ne),ns2)
      if(ns2/=ns) then
        write(*,*)'Mismatch in side:',ns2,ns
        stop
      endif
      allocate(ic3(4,ne),elside(4,ne),isdel(2,ns2),isidenode(2,ns2),xcj(ns2),ycj(ns2),stat=istat)
      if(istat/=0) stop 'Allocation error: side(0)'
      call schism_geometry_single(np,ne,ns2,real(x),real(y),i34,elnode(1:4,1:ne),ic3(1:4,1:ne), &
     &elside(1:4,1:ne),isdel,isidenode,xcj,ycj)

      if(ibase==1) then !node based
        last_dim=np
      else if(ibase==2) then
        last_dim=ne
      else
        last_dim=ns
      endif 
      print*, '# of records in each stack is: ',nrec

      allocate(timeout(nrec),ztot(nvrt),sigma(nvrt),sigma_lcl(nvrt,np),kbp(np), &
     &outvar(2,nvrt,last_dim),eta2(np),ztmp(nvrt,np),residual(2,last_dim), &
     &icounter(last_dim),idry(np),xctr(ne),yctr(ne))
      outvar=-huge(1.0) !test mem

      !Center
      do i=1,ne
        xctr(i)=sum(x(elnode(1:i34(i),i)))/i34(i)
        yctr(i)=sum(y(elnode(1:i34(i),i)))/i34(i)
      enddo !i

      call get_vgrid_single('vgrid.in',np,nvrt,ivcor,kz,h_s,h_c,theta_b,theta_f,ztot,sigma,sigma_lcl,kbp)

      print*, 'last element',elnode(1:i34(ne),ne)

!...  Time iteration
!...
      out=-99 !init.
      ztmp=-99
      residual=0
      icounter=0 !counter for each node (wet/dry)
      do iday=iday1,iday2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if(icomb==0) then !uncombined
        call get_timeout(iday,nrec,timeout,icomb)
      else
        call get_timeout(iday,nrec,timeout)
      endif
!      print*, 'time=',timeout

      do irec=1,nrec
!----------------------------------------------------------------------------
        if(icomb==0) then !uncombined
          do irank=0,nproc-1
            call get_outvar(iday,irec,varname,np,last_dim,nvrt,outvar,i23d,ivs,eta2,irank)

            !You can also read in additional var here
            !call get_outvar(iday,irec,varname2,np,last_dim,nvrt,outvar2,i23d,ivs,eta2,irank)
          enddo !irank
        else
          call get_outvar(iday,irec,varname,np,last_dim,nvrt,outvar,i23d,ivs,eta2)
          !You can also read in additional var here; similar to uncombined above
        endif

        !Check ibase
        if((i23d-1)/3+1 /= ibase) then
          print*, 'Incorrect input ibase:',i23d,ibase
          stop
        endif

        !Available now: outvar(2,nvrt,np|ne|ns),i23d,ivs,eta2(np)

        if(mod(i23d-1,3)==0) then !2D
          do i=1,last_dim
            icounter(i)=icounter(i)+1
            residual(1:ivs,i)=residual(1:ivs,i)+outvar(1:ivs,1,i)
          enddo !i
        else !3D
          !Compute z coordinates, which may be used for node-based vars
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

            !Following is code for vertical interpolation
!            if(idry(i)==0) then !wet
!              !Interplate in vertical
!              if(z00>=ztmp(nvrt,i)) then !above F.S.
!                k0=nvrt-1; rat=1
!              else if(z00<=ztmp(kbp(i),i)) then !below bottom; extrapolate
!                k0=kbp(i); rat=0
!              else !above bottom; cannot be above F.S.
!                k0=0
!                do k=kbp(i),nvrt-1
!                  if(z00>=ztmp(k,i).and.z00<=ztmp(k+1,i)) then
!                    k0=k
!                    rat=(z00-ztmp(k,i))/(ztmp(k+1,i)-ztmp(k,i))
!                    exit
!                  endif
!                enddo !k
!              endif !ztmp
!
!              if(k0==0) then
!                write(*,*)'failed to find a vertical level:',irec,i,ifs,z2,ztmp(:,i)
!                stop
!              endif
!            endif !idry
             !Vertical interp
!            tmp=outvar(1,k0,i)*(1-rat)+outvar(1,k0+1,i)*rat
          enddo !i=1,np

          !Depth-averged var as an e.g., assuming @ whole level (if @
          !half level, simply change below by assuming constant in each prism as in F.V.)
          do i=1,last_dim
            wild(:)=0
            htot=0
            if(ibase==1) then !node
              if(idry(i)==0) then
                htot=ztmp(nvrt,i)-ztmp(kbp(i),i)
                do k=kbp(i),nvrt-1
                  wild(1:ivs)=wild(1:ivs)+(ztmp(k+1,i)-ztmp(k,i))*(outvar(1:ivs,k,i)+outvar(1:ivs,k+1,i))*0.5
                enddo !k
              endif !idry
            else if(ibase==2) then !elem
              if(maxval(idry(elnode(1:i34(i),i)))==0) then ! all wet
                htot=0
                ktmp=maxval(kbp(elnode(1:i34(i),i)))
                do k=ktmp,nvrt-1
                  tmp=sum(ztmp(k+1,elnode(1:i34(i),i))-ztmp(k,elnode(1:i34(i),i)))/i34(i) !layer thickness
                  htot=htot+tmp
                  wild(1:ivs)=wild(1:ivs)+(outvar(1:ivs,k,i)+outvar(1:ivs,k+1,i))*0.5*tmp    
                enddo !k
              endif !all wet
            else !side
              if(maxval(idry(isidenode(1:2,i)))==0) then ! all wet
                htot=0
                ktmp=maxval(kbp(isidenode(1:2,i)))
                do k=ktmp,nvrt-1
                  tmp=sum(ztmp(k+1,isidenode(1:2,i))-ztmp(k,isidenode(1:2,i)))*0.5
                  htot=htot+tmp
                  wild(1:ivs)=wild(1:ivs)+(outvar(1:ivs,k,i)+outvar(1:ivs,k+1,i))*0.5*tmp
                enddo !k
              endif !all wet
            endif !ibase

            if(htot/=0) then
              residual(1:ivs,i)=residual(1:ivs,i)+wild(1:ivs)/htot
              icounter(i)=icounter(i)+1
            endif
          enddo !i
        endif !2/3D
!----------------------------------------------------------------------------
      enddo !irec
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      enddo !iday

!     Output
      do i=1,last_dim
        if(icounter(i)==0) then
          residual(:,i)=0
        else
          residual(:,i)=residual(:,i)/icounter(i)
        endif !icounter
      enddo !i

      open(65,file='average.out')
      if(ibase==1) then !nodes
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
            write(65,*)x(i),y(i),residual(1:2,i)
          enddo !i
        endif
        close(65)
      else if(ibase==2) then !elem
        write(65,*)iday1,iday2
        write(65,*)ne
        do i=1,ne
          write(65,*)i,xctr(i),yctr(i),residual(1:ivs,i)
        enddo !i
      else !side
        write(65,*)iday1,iday2
        write(65,*)ns
        do i=1,ns
          write(65,*)i,xcj(i),ycj(i),residual(1:ivs,i) 
        enddo !i
      endif !ibase

      print*, 'Finished!'

      stop
      end
