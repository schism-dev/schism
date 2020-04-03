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
!       Works for mixed tri/quad outputs on NODE/ELEMENT based vars.
!       Inputs: screen; combined or uncombined nc file; vgrid.in (in this dir or ../)
!       Outputs: extract.out (ascii); optional: max/min/avg 2D output, with
        !       self-explanatory file names
!       History: (1) added non-standard outputs (April 2012) - transparent to most scripts
!               as format is same; (2) added ivcor=1 (Dec 2013); (3)
!               added quads (Nov. 2014) (4) changed to nc outputs (Sept
!               2017); (5) added uncombined option (Feb 2019); (6) added special
!               treatment on max_elev and other min/max/avg for other vars (Mar,
!               2020)
!****************************************************************************************
!     ifort -O2 -assume byterecl -o read_output8_allnodes.exe ../UtilLib/extract_mod.f90 ../UtilLib/schism_geometry.f90 ../UtilLib/compute_zcor.f90 read_output8_allnodes.f90 -I$NETCDF/include  -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff
!     ifort -g -assume byterecl -o read_output8_allnodes.exe ../UtilLib/extract_mod.f90 ../UtilLib/schism_geometry.f90 ../UtilLib/compute_zcor.f90 read_output8_allnodes.f90 -I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff

      program read_out
      use netcdf
      use extract_mod
      use schism_geometry_mod
      use compute_zcor

!      parameter(nbyte=4)
      character(len=30) :: file63,varname,outname(3)
      character(len=12) :: it_char
!      character*48 data_format
  
!      integer,allocatable :: i34(:),elnode(:,:)
      integer :: nx(4,4,3)
      allocatable :: sigma(:),cs(:),ztot(:)
      allocatable :: outvar(:,:,:,:),out(:,:,:),icum(:,:,:),eta2(:,:),ztmp(:,:)
      allocatable :: xcj(:),ycj(:),dps(:),kbs(:),xctr(:),yctr(:),dpe(:),kbe(:)
      allocatable :: zs(:,:),ze(:,:),idry(:),outs(:,:,:),oute(:,:,:),rstat2d(:,:),rviolator(:)
      allocatable :: ic3(:,:),isdel(:,:),isidenode(:,:),kbp(:),sigma_lcl(:,:)
      integer, allocatable :: elside(:,:),idry_e(:),nne(:),indel(:,:),iviolator(:)
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

      print*, 'Is the var node-based or ele-based? 0: node based; 1: element based'
      read(*,*)inode_ele
!!'
      print*, '<<<<<inode_ele: ',inode_ele

      varname=adjustl(varname); len_var=len_trim(varname)
      
      print*, 'Input start and end stack # to read:'
      read(*,*) iday1,iday2
      print*, '<<<<<start and end stack #: ',iday1,iday2

      print*, 'Do you want to compute stats for 2D var? (0/1:min; 0/1:max; 0/1:time avg)'
!!'
      read(*,*) icomp_stats(1),icomp_stats(2),icomp_stats(3)
      print*, '<<<<<icomp_stats: ',icomp_stats
      outname(1)='min';outname(2)='max';outname(3)='avg'

      print*, 'Input a threshold: values below the threshold will be output as a bp file.'
!!'
      read(*,*) thres
      print*, '<<<<<threshold: ',thres

      print*, 'How do you want to initialize the variable values for min/max?'
!!'
      print*, '0: do nothing;  1: intialized to -dp (useful for maxelev)'
!!'
      read(*,*) iinitial
!!
      print*, '<<<<<initialization flag: ',iinitial
      if (iinitial.ne.0 .and. iinitial.ne.1) stop 'wrong initialization flag'

!!'
      if(varname(1:len_var).eq.'elev' .and. iinitial==1) then
        is_elev=1 
        print*, '<<<<<special treatment will be implemented for elev'
      else
        is_elev=0 
      endif


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

      print*
      print*, 'After header:',ne,np,nrec,i34(ne),elnode(1:i34(ne),ne),nvrt,h0,x(np),y(np),dp(np) !,start_time

      ! Compute geometry
      allocate(nne(np))
      do k=3,4
        do i=1,k
          do j=1,k-1
            nx(k,i,j)=i+j
            if(nx(k,i,j)>k) nx(k,i,j)=nx(k,i,j)-k
            if(nx(k,i,j)<1.or.nx(k,i,j)>k) then
              write(*,*)'nx wrong',i,j,k,nx(k,i,j)
              stop
            endif
          enddo !j
        enddo !i
      enddo !k

      nne=0
      do i=1,ne
        do j=1,i34(i)
          nd=elnode(j,i)
          nne(nd)=nne(nd)+1
        enddo
      enddo
      mnei=maxval(nne)

      allocate(indel(mnei,np),stat=istat)
      if(istat/=0) stop 'Failed to alloc. indel'
      nne=0
      do i=1,ne
        do j=1,i34(i)
          nd=elnode(j,i)
          nne(nd)=nne(nd)+1
          if(nne(nd)>mnei) then
            write(*,*)'Too many neighbors',nd
            stop
          endif
          indel(nne(nd),nd)=i
        enddo
      enddo !i

      call compute_nside(np,ne,i34,elnode(1:4,1:ne),ns2)
      if(ns2/=ns) then
        write(*,*)'Mismatch in side:',ns2,ns
        stop
      endif
      allocate(ic3(4,ne),elside(4,ne),isdel(2,ns2),isidenode(2,ns2),xcj(ns2),ycj(ns2),stat=istat)
      if(istat/=0) stop 'Allocation error: side(0)'
      call schism_geometry_single(np,ne,ns2,x,y,i34,elnode(1:4,1:ne),ic3(1:4,1:ne), &
     &elside(1:4,1:ne),isdel,isidenode,xcj,ycj)

      !For dimensioning purpose
      if(np>ns.or.ne>ns) stop 'ns is not largest'

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
     &rviolator(ne),iviolator(ne))
      outvar=-huge(1.0) !test mem

!     Read vgrid.in
      call get_vgrid_single('vgrid.in',np,nvrt,ivcor,kz,h_s,h_c,theta_b,theta_f,ztot,sigma,sigma_lcl,kbp)

!     Calculate kbp00
      if(ivcor==1) then
        kbp00=kbp
      else
        do i=1,np
          call zcor_SZ_single(dp(i),1.e8,h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot,sigma, &
     &ztmp(:,i),idry2,kbp00(i))
        enddo !i
      endif !ivcor
      
!...  Time iteration
!...
      outvar=-99 !init.
      ztmp=-99
      if(iinitial==1) then
        print*,'setting initial elev = -dp for all (dry and wet) nodes'
        rstat2d(1,1:np)=-dp !for elev, min is -h
        rstat2d(2,1:np)=-dp !for elev, min is -h
        rstat2d(3,1:np)=-dp !for elev, min is -h
      else
        if (icomp_stats(1)==1) then  !min
            rstat2d(1,:)=huge(1.0)
        elseif (icomp_stats(2)==1) then !max
            rstat2d(2,:)=-huge(1.0)
        elseif (icomp_stats(3)==1) then !avg
            rstat2d(3,:)=0.0
        endif
      endif

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
            call get_outvar_multirecord(1,iday,varname,irec1,irec2,np,last_dim,nvrt,nrec3,outvar,i23d,ivs,eta2,irank)
          enddo !irank
        else
          call get_outvar_multirecord(1,iday,varname,irec1,irec2,np,last_dim,nvrt,nrec3,outvar,i23d,ivs,eta2)
        endif

        !Available now:
        !outvar(2,nvrt,np|ne,irec2-irec1+1),i23d,ivs,eta2(np,irec2-irec1+1)

        do irec=1,irec2-irec1+1 !offeset record #
          irec_real=irec1+irec-1 !actual record #

          if(mod(i23d-1,3)==0) then !2D
!           Output: time, 2D variable at all nodes
            write(65,'(e14.6,1000000(1x,e14.4))')time/86400,((outvar(m,1,i,irec),m=1,ivs),i=1,np)

            !Compute stats (magnitude for vectors)
            if(sum(icomp_stats(:))/=0) then
              if(is_elev==1 .and. icomp_stats(2)==1) then !maxelev

                !Mark wet elem
                do i=1,ne
                  if(minval(outvar(1,1,elnode(1:i34(i),i),irec)+dp(elnode(1:i34(i),i)))>0) then
                    idry_e(i)=0
                  else
                    idry_e(i)=1
                  endif 
                enddo !i
                do i=1,ne
                  if(idry_e(i)==0) then !wet
                    ifl=1 !isolated wet
                    loop2: do j=1,i34(i)
                      nd=elnode(j,i)
                      do m=1,nne(nd)
                        ie=indel(m,nd)
                        if(i/=ne.and.idry_e(ie)==0) then
                          ifl=0
                          exit loop2
                        endif
                      enddo !m
                    enddo loop2 !j

                    if(ifl==0) then ! not isolated wet
                      do j=1,i34(i)
                        nd=elnode(j,i)
                        tmp=outvar(1,1,nd,irec)
                        if(tmp+dp(nd)>0) rstat2d(2,nd)=max(rstat2d(2,nd),tmp)
                      enddo
                    endif 
                  endif !idry_e
                enddo !i=1,ne 
              else !other variables, or min_ele or avg_ele
                do i=1,last_dim
                  if(ivs==2) then
                    tmp=sqrt(outvar(1,1,i,irec)**2+outvar(2,1,i,irec)**2)
                  else
                    tmp=outvar(1,1,i,irec)
                  endif

                  if (icomp_stats(1)==1) then !min
                    rstat2d(1,i)=min(rstat2d(1,i),tmp)
                  endif
                  if(icomp_stats(2)==1) then !max
                    rstat2d(2,i)=max(rstat2d(2,i),tmp)
                  endif
                  if(icomp_stats(3)==1) then !min
                    rstat2d(3,i)=tmp+rstat2d(3,i)
                  endif

                enddo !i 
              endif !is_elev
            endif ! icomp_stats
          else !if(i23d==3) then !3D 
            !Compute z coordinates
            do i=1,last_dim
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
                call zcor_SZ_single(dp(i),eta2(i,irec),h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot, &
     &sigma,ztmp(:,i),idry(i),kbpl)
              endif

              do k=max0(1,kbp00(i)),nvrt
!               Output: time, node #, level #, z-coordinate, 3D variable 
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
        enddo !irec

        if(irec2==nrec) exit loop1
        irec1=irec1+nrec3
      end do loop1

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      enddo !iday

      
      nrec_total=nrec*(iday2-iday1+1)
      if(icomp_stats(3)==1) then
        rstat2d(3,:)=rstat2d(3,:)/real(nrec_total)
      endif
      print*, '<<<<<total rec: ', nrec_total

      do k=1,3
        if(icomp_stats(k)==0) cycle

        if (inode_ele==0) then !node-based, write gr3
          open(12,file=trim(adjustl(varname))//'_'//trim(adjustl(outname(k)))//'.gr3',status='replace')
          write(12,*)iday1,iday2
          write(12,*)ne,np
          do i=1,np
            write(12,*)i,x(i),y(i),rstat2d(k,i)
          enddo !i
          do i=1,ne
            write(12,*)i,i34(i),elnode(1:i34(i),i)
          enddo !i
          close(12)
        elseif(inode_ele==1) then !ele based, write prop
          open(12,file=trim(adjustl(varname))//'_'//trim(adjustl(outname(k)))//'.prop',status='replace')
          do i=1,ne
            write(12,*)i,rstat2d(k,i)
          enddo
        endif

        close(12)
      enddo !k=1,3

      !record violators
      do k=1,3
        nrec=0; rviolator=0.0; iviolator=0
        if (icomp_stats(k)==0) cycle
        do i=1,last_dim
          if (rstat2d(k,i)<thres) then
            nrec=nrec+1
            rviolator(nrec)=rstat2d(k,i) !val
            iviolator(nrec)=i !id
          endif
        enddo
        
        open(13,file=trim(adjustl(varname))//'_'//trim(adjustl(outname(k)))//'_violators.bp',status='replace')
        write(13,*) 'violators <', thres
        write(13,*) nrec

        if (inode_ele==0) then !node-based
          do i=1,nrec
            ip=iviolator(i)
            write(13,*) i,x(ip),y(ip),rviolator(i),ip
          enddo
        elseif (inode_ele==1) then !ele based
          do i=1,nrec
            ie=iviolator(i)
            write(13,*) i,sum(x(elnode(1:i34(ie),ie)))/i34(ie),&
              & sum(y(elnode(1:i34(ie),ie)))/i34(ie), rviolator(i),ie
          enddo
        endif
        close(13)
      enddo !k=1,3

      print*, 'Finished!'

      stop
      end
