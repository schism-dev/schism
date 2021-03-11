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

!       Simplified read_output8_allnodes.f90 for 2D vars only, with
!       flexible inputs to facilitate parallelization.
!	Read (combined or uncombined) nc outputs for multiple files at all nodes 
!       Works for mixed tri/quad outputs on NODE based vars.
!       Inputs: screen; combined or uncombined nc file; filter_flag (for
!       filtering outputs), output file name
!       Outputs: time series (ascii); 
!       History: (1) added non-standard outputs (April 2012) - transparent to most scripts
!               as format is same; (2) added ivcor=1 (Dec 2013); (3)
!               added quads (Nov. 2014) (4) changed to nc outputs (Sept
!               2017); (5) added uncombined option (Feb 2019); (6) added special
!               treatment on max_elev and other min/max/avg for other vars (Mar,
!               2020)
!               (7) simplified the code
!****************************************************************************************
!     ifort -O2 -assume byterecl -o read_output8_allnodes_simple ../UtilLib/extract_mod.f90 ../UtilLib/schism_geometry.f90 ../UtilLib/compute_zcor.f90 read_output8_allnodes_simple.f90 -I$NETCDF/include  -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff

      program read_out
      use netcdf
      use extract_mod
      use schism_geometry_mod
      use compute_zcor

!      parameter(nbyte=4)
      character(len=30) :: file63,varname,outname(3)
      character(len=100) :: filter_flag
      character(len=12) :: it_char
!      character*48 data_format
  
!      integer,allocatable :: i34(:),elnode(:,:)
      integer :: nx(4,4,3)
      allocatable :: sigma(:),cs(:),ztot(:)
      allocatable :: outvar(:,:,:,:),outvar2(:,:),icum(:,:,:),eta2(:,:),ztmp(:,:)
      allocatable :: xcj(:),ycj(:),dps(:),kbs(:),xctr(:),yctr(:),dpe(:),kbe(:)
      allocatable :: zs(:,:),ze(:,:),idry(:),outs(:,:,:),oute(:,:,:),rstat2d(:,:),rviolator(:)
      allocatable :: ic3(:,:),isdel(:,:),isidenode(:,:),kbp(:),sigma_lcl(:,:)
      integer, allocatable :: elside(:,:),idry_e(:),nne(:),indel(:,:),iviolator(:),include2(:)
      real*8,allocatable :: timeout(:)

      integer :: icomp_stats(3)

      print*, 'Do you work on uncombined (0) or combined (1) nc?'
      read(*,*)icomb
      if(icomb==0) then
        print*,'<<<<<uncombined'
      elseif(icomb==1) then
        print*,'<<<<<combined'
      endif

!...  Set max array size for system memory
!...
      print*, 'Recommendation: for uncombined nc, specify max array size (e.g., <=2.e9);'
      print*, 'for combined nc, specify # of records to read each time.'
      print*, 'Do you want to specify max array size (1) or # of records (2)'
      read(*,*)ispec_max
      if(ispec_max==1) then
        print*, 'Input max array size (e.g., <=2.e9):'
        read(*,*)max_array_size
        print*, '<<<<<max_array_size read in=',max_array_size
      else
        print*, 'Input # of records:'
        read(*,*)nrec3
        print*, '<<<<<# of records=',nrec3
      endif

      print*, 'Input NODE-based variable name to read from nc (e.g. elev):'
      read(*,'(a30)')varname
!!'
      print*, '<<<<<var name: ',varname

!      print*, 'Is the var node-based or ele-based? 0: node based; 1: element based'
!      read(*,*)inode_ele
!!!'
!      print*, '<<<<<inode_ele: ',inode_ele
      inode_ele=0

      varname=adjustl(varname); len_var=len_trim(varname)
      
      print*, 'Input start and end stack # to read:'
      read(*,*) iday1,iday2
      print*, '<<<<<start and end stack #: ',iday1,iday2

      print*, 'name for output time series:'
      read(*,*) file63
      print*, 'Output name is:',file63

      print*, 'File name for filter outputs:'
      read(*,*) filter_flag
      print*, 'Filter input is:',filter_flag

!      print*, 'Do you want to compute stats for 2D var? (0/1:min; 0/1:max; 0/1:time avg)'
!!!'
!      read(*,*) icomp_stats(1),icomp_stats(2),icomp_stats(3)
!      print*, '<<<<<icomp_stats: ',icomp_stats
!      outname(1)='min';outname(2)='max';outname(3)='avg'
!
!      print*, 'Input a threshold: values below the threshold will be output as a bp file.'
!!!'
!      read(*,*) thres
!      print*, '<<<<<threshold: ',thres
!
!      print*, 'How do you want to initialize the variable values for min/max?'
!!!'
!      print*, '0: do nothing;  1: intialized to -dp (useful for maxelev)'
!!!'
!      read(*,*) iinitial
!!
!      print*, '<<<<<initialization flag: ',iinitial
!      if (iinitial.ne.0 .and. iinitial.ne.1) stop 'wrong initialization flag'

!!'
!      if(varname(1:len_var).eq.'elev' .and. iinitial==1) then
!        is_elev=1 
!        print*, '<<<<<special treatment will be implemented for elev'
!      else
!        is_elev=0 
!      endif
 
      open(65,file=trim(adjustl(file63)),status='replace')
      
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

      print*
      print*, 'After header:',ne,np,nrec,i34(ne),elnode(1:i34(ne),ne),nvrt,h0,x(np),y(np),dp(np) !,start_time

      !For dimensioning purpose
!      if(np>ns.or.ne>ns) stop 'ns is not largest'

      if (inode_ele==0) then
        last_dim=np
      elseif (inode_ele==1) then
        last_dim=ne
      endif

      print*, '<<<<<last_dim: ',last_dim

!     nrec specified if ispec_max=2
      if(ispec_max==1) nrec3=max_array_size/(2*nvrt*last_dim)-1
      nrec3=min(nrec,max(nrec3,1))

      allocate(idry_e(ne),rstat2d(3,ne),ztot(nvrt),sigma(nvrt),sigma_lcl(nvrt,np),kbp(np), &
     &outvar(2,nvrt,last_dim,nrec3),eta2(np,nrec3),idry(np),ztmp(nvrt,np),timeout(nrec), &
     &rviolator(ne),iviolator(ne),include2(np),outvar2(2,np))
      outvar=-huge(1.0) !test mem

!     Read in filtering flags for output
      open(60,file=trim(adjustl(filter_flag)),status='old')
      do i=1,np
        read(60,*)include2(i)
      enddo !i
      close(60)

!     Read vgrid.in
!      call get_vgrid_single('vgrid.in',np,nvrt,ivcor,kz,h_s,h_c,theta_b,theta_f,ztot,sigma,sigma_lcl,kbp)
!
!!     Calculate kbp00
!      if(ivcor==1) then
!        kbp00=kbp
!      else
!        do i=1,np
!          call zcor_SZ_single(dp(i),1.e8,h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot,sigma, &
!     &ztmp(:,i),idry2,kbp00(i))
!        enddo !i
!      endif !ivcor
      
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

      print*, 'doing stack # ',iday

      irec1=1 !start record
      loop1: do
        irec2=min(nrec,irec1+nrec3-1)

        if(icomb==0) then !uncombined
          do irank=0,nproc-1
            call get_outvar_multirecord(iday,varname,irec1,irec2,np,last_dim,nvrt,nrec3,outvar,i23d,ivs,eta2,irank)
          enddo !irank
        else
          call get_outvar_multirecord(iday,varname,irec1,irec2,np,last_dim,nvrt,nrec3,outvar,i23d,ivs,eta2)
        endif

        !Available now:
        !outvar(2,nvrt,np|ne,irec2-irec1+1),i23d,ivs,eta2(np,irec2-irec1+1)

        do irec=1,irec2-irec1+1 !offeset record #
          irec_real=irec1+irec-1 !actual record #

          if(mod(i23d-1,3)==0) then !2D
!           Output: time, 2D variable at selected nodes
            icount=0
            do i=1,np
              if(include2(i)/=0) then
                icount=icount+1
                outvar2(1:ivs,icount)=outvar(1:ivs,1,i,irec)
              endif
            enddo !i 
            !write(65,'(e14.6,10000000(1x,e14.4))')timeout(irec_real)/86400,((outvar(m,1,i,irec),m=1,ivs),i=1,np)
            write(65,'(e14.6,10000000(1x,e14.4))')timeout(irec_real)/86400,outvar2(1:ivs,1:icount)

          else !if(i23d==3) then !3D 
            stop 'Cannot be 3D var'
          endif !i23d
        enddo !irec

        if(irec2==nrec) exit loop1
        irec1=irec1+nrec3
      end do loop1

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      enddo !iday

      print*, 'Finished!'

      stop
      end
