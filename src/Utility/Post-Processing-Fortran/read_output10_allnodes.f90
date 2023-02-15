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

!	Read nc outputs from scribe I/O versions for multiple files at all nodes
!       Works for mixed tri/quad outputs on NODE based vars.
!       Inputs: screen; nc file and out2d*.nc; vgrid.in (in this dir or ../)
!       Outputs: extract.out (ascii); optional: max/min/avg 2D output (*_[max|min|avg].gr3 and *_violators.bp)
!****************************************************************************************
!     ifort -O2 -assume byterecl -o read_output10_allnodes.exe ../UtilLib/extract_mod2.f90 ../UtilLib/schism_geometry.f90 ../UtilLib/compute_zcor.f90 read_output10_allnodes.f90 -I$NETCDF/include  -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff

      program read_out
      use netcdf
      use extract_mod2
      use schism_geometry_mod
      use compute_zcor

!      parameter(nbyte=4)
      character(len=30) :: file63,varname,varname2,outname(3),file62,file64
      character(len=12) :: it_char
!      character*48 data_format
      logical::lexist,lexist2
      integer :: nx(4,4,3),dimids(100),idims(100)
      allocatable :: sigma(:),cs(:),ztot(:)
      allocatable :: outvar(:,:,:),out(:,:,:),icum(:,:,:),eta2(:),ztmp(:,:)
      allocatable :: xcj(:),ycj(:),dp(:),dps(:),kbs(:),xctr(:),yctr(:),dpe(:),kbe(:)
      allocatable :: zs(:,:),ze(:,:),idry(:),outs(:,:,:),oute(:,:,:),rstat2d(:,:),rviolator(:)
      allocatable :: ic3(:,:),isdel(:,:),isidenode(:,:),kbp(:),sigma_lcl(:,:),kbp00(:)
      integer, allocatable :: elside(:,:),idry_e(:),nne(:),indel(:,:),iviolator(:)
      integer,allocatable :: i34(:),elnode(:,:)
      integer :: char_len,start_2d(2),start_3d(3),start_4d(4), &
     &count_2d(2),count_3d(3),count_4d(4)
      real*8,allocatable :: timeout(:),xnd(:),ynd(:)

      integer :: icomp_stats(3)

      print*, 'Input NODE-based variable name to read from nc (e.g. elevation):'
      print*, '(if vector, input X component)'
      read(*,'(a30)')varname
!!'
      print*, '<<<<<var name: ',varname

      print*, 'Is variable scalar (1) or vector (2)?'
      read(*,*) ivs
      if(ivs/=1.and.ivs/=2) stop 'check ivs'     

!      print*, 'Is the var node-based or ele-based? 0: node based; 1: element based'
!      read(*,*)inode_ele
!      print*, '<<<<<inode_ele: ',inode_ele
      inode_ele=0

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
      if(varname(1:len_var).eq.'elevation' .and. iinitial==1) then
        is_elev=1 
        print*, '<<<<<special treatment will be implemented for elev'
      else
        is_elev=0 
      endif

      open(65,file='extract.out')
      
!...  Get basic info from out2d*.nc
      !Returned vars: ne,np,ns,nrec,[xnd ynd dp](np),
      !elnode,i34,nvrt,h0,dtout,kbp
      call get_dims(iday1,np,ne,ns,nvrt,h0)
      allocate(xnd(np),ynd(np),dp(np),kbp(np),i34(ne),elnode(4,ne),stat=istat)
      if(istat/=0) stop 'alloc (1)'
      call readheader(iday1,np,ne,ns,kbp,i34,elnode,nrec,xnd,ynd,dp,dtout)

      print*
      print*, 'After header:',ne,np,nrec,i34(ne),elnode(1:i34(ne),ne),nvrt,h0,xnd(np),ynd(np),dp(np) !,start_time

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
      call schism_geometry_single(np,ne,ns2,real(xnd),real(ynd),i34,elnode(1:4,1:ne),ic3(1:4,1:ne), &
     &elside(1:4,1:ne),isdel,isidenode,xcj,ycj)

      !For dimensioning purpose
      if(np>ns.or.ne>ns) stop 'ns is not largest'

!      if (inode_ele==0) then
      last_dim=np
!      elseif (inode_ele==1) then
!        last_dim=ne
!      endif

      print*, '<<<<<last_dim: ',last_dim

      allocate(idry_e(ne),rstat2d(3,ne),ztot(nvrt),sigma(nvrt),sigma_lcl(nvrt,np), &
     &outvar(nvrt,last_dim,ivs),eta2(np),idry(np),ztmp(nvrt,np),timeout(nrec), &
     &rviolator(ne),iviolator(ne),kbp00(np))
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
      write(it_char,'(i12)')iday
      it_char=adjustl(it_char)
      leng=len_trim(it_char)
      file62='out2d_'//it_char(1:leng)//'.nc'
      iret=nf90_open(trim(adjustl(file62)),OR(NF90_NETCDF4,NF90_NOWRITE),ncid4)
      if(iret/=nf90_NoErr) stop 'Failed to open file62'
      !time is double
      iret=nf90_inq_varid(ncid4,'time',itime_id)
      iret=nf90_get_var(ncid4,itime_id,timeout,(/1/),(/nrec/))
      print*, 'time=',timeout !,trim(adjustl(file63))
 
      !Find nc file
      file63=varname(1:len_var)//'_'//it_char(1:leng)//'.nc'
      if(ivs==2) varname2=varname(1:len_var-1)//'Y'
      inquire(file=file63,exist=lexist)
      if(lexist) then
        i23d=2 !3D var
        if(ivs==2) then
          file64=varname2(1:len_var)//'_'//it_char(1:leng)//'.nc'
          inquire(file=file64,exist=lexist2)
          if(.not.lexist2) then
            print*, 'Missing y-component:',file64
            print*, file63
            stop
          endif
        endif
      else
        i23d=1 !2D
        file63=file62
        file64=file62
      endif

      iret=nf90_open(trim(adjustl(file63)),OR(NF90_NETCDF4,NF90_NOWRITE),ncid)
      if(iret/=nf90_NoErr) stop 'Failed to open file63'
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

      if(ivs==2) then !vector
        iret=nf90_open(trim(adjustl(file64)),OR(NF90_NETCDF4,NF90_NOWRITE),ncid2)
        if(iret/=nf90_NoErr) stop 'Failed to open file64'
        iret=nf90_inq_varid(ncid2,varname2(1:len_var),ivarid2)
        if(iret/=nf90_NoErr) stop 'Var2 not found'
      endif !ivs

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
          iret=nf90_get_var(ncid,ivarid1,outvar(1,1:npes,1),start_2d,count_2d)
          if(ivs==2) iret=nf90_get_var(ncid2,ivarid2,outvar(1,1:npes,2),start_2d,count_2d)
        else !3D
          start_3d(1:2)=1; start_3d(3)=irec
          count_3d(1)=nvrt; count_3d(2)=npes; count_3d(3)=1
          iret=nf90_get_var(ncid,ivarid1,outvar(:,1:npes,1),start_3d,count_3d)
          if(ivs==2) iret=nf90_get_var(ncid2,ivarid2,outvar(:,1:npes,2),start_3d,count_3d)
        endif 

        !Available now: outvar(nvrt,np,ivs), eta2(np)
          if(i23d==1) then !2D
!           Output: time, 2D variable at all nodes/elem
            if(maxval(icomp_stats(1:3))/=0) &
     &write(65,'(e14.6,1000000(1x,e14.4))')timeout(irec)/86400,(outvar(1,i,1),i=1,np)

            !Compute stats (magnitude for vectors)
            if(sum(icomp_stats(:))/=0) then
              if(is_elev==1 .and. icomp_stats(2)==1) then !maxelev

                !Mark wet elem
                do i=1,ne
                  !if(minval(outvar(1,elnode(1:i34(i),i),1)+dp(elnode(1:i34(i),i)))>0) then
                  if(minval(eta2(elnode(1:i34(i),i))+dp(elnode(1:i34(i),i)))>0) then
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
                        tmp=outvar(1,nd,1)
                        if(tmp+dp(nd)>0) rstat2d(2,nd)=max(rstat2d(2,nd),tmp)
                      enddo
                    endif 
                  endif !idry_e
                enddo !i=1,ne 
              else !other variables, or min_ele or avg_ele
                do i=1,last_dim
                  tmp=outvar(1,i,1)
                  if (icomp_stats(1)==1) then !min
                    rstat2d(1,i)=min(rstat2d(1,i),tmp)
                  endif
                  if(icomp_stats(2)==1) then !max
                    rstat2d(2,i)=max(rstat2d(2,i),tmp)
                  endif
                  if(icomp_stats(3)==1) then !average
                    rstat2d(3,i)=tmp+rstat2d(3,i)
                  endif
                enddo !i 
              endif !is_elev
            endif ! icomp_stats
          else !3D 
            !Compute z coordinates
            do i=1,last_dim
              if(ivcor==1) then !localized
                if(dp(i)+eta2(i)<=h0) then
                  idry(i)=1
                else !wet
                  idry(i)=0
                  do k=kbp(i),nvrt
                    ztmp(k,i)=(eta2(i)+dp(i))*sigma_lcl(k,i)+eta2(i)
                  enddo !k
                endif !wet/dry
              else if(ivcor==2) then !SZ
                call zcor_SZ_single(dp(i),eta2(i),h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot, &
     &sigma,ztmp(:,i),idry(i),kbpl)
              endif

              do k=max0(1,kbp00(i)),nvrt
!               Output: time, node #, level #, z-coordinate, 3D variable 
                if(idry(i)==1) then
                  write(65,*)timeout(irec)/86400,i,k,-1.e6,outvar(k,i,1:ivs)
                else
                  write(65,*)timeout(irec)/86400,i,k,ztmp(k,i),outvar(k,i,1:ivs)
                endif
              enddo !k
 
              !Debug
              !write(65,*)x(i),y(i),out(i,1,1:ivs) 

            enddo !i
          endif !i23d
      enddo !irec

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      enddo !iday

      
      nrec_total=nrec*(iday2-iday1+1)
      if(icomp_stats(3)==1) then
        rstat2d(3,:)=rstat2d(3,:)/real(nrec_total)
      endif
      print*, '<<<<<total rec: ', nrec_total

      do k=1,3
        if(icomp_stats(k)==0) cycle

!        if (inode_ele==0) then !node-based, write gr3
          open(12,file=trim(adjustl(varname))//'_'//trim(adjustl(outname(k)))//'.gr3',status='replace')
          write(12,*)iday1,iday2
          write(12,*)ne,np
          do i=1,np
            write(12,*)i,xnd(i),ynd(i),rstat2d(k,i)
          enddo !i
          do i=1,ne
            write(12,*)i,i34(i),elnode(1:i34(i),i)
          enddo !i
          close(12)
!        elseif(inode_ele==1) then !ele based, write prop
!          open(12,file=trim(adjustl(varname))//'_'//trim(adjustl(outname(k)))//'.prop',status='replace')
!          do i=1,ne
!            write(12,*)i,rstat2d(k,i)
!          enddo
!        endif

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

!        if (inode_ele==0) then !node-based
        do i=1,nrec
          ip=iviolator(i)
          write(13,*) i,xnd(ip),ynd(ip),rviolator(i),ip
        enddo
!        elseif (inode_ele==1) then !ele based
!          do i=1,nrec
!            ie=iviolator(i)
!            write(13,*) i,sum(xnd(elnode(1:i34(ie),ie)))/i34(ie),&
!              & sum(ynd(elnode(1:i34(ie),ie)))/i34(ie), rviolator(i),ie
!          enddo
!        endif
        close(13)
      enddo !k=1,3

      print*, 'Finished!'

      stop
      end
