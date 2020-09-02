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
!       Outputs: fort.11 (fatal errors); fort.12: nonfatal errors.
!                out (diagnostic outputs);
!                f1[89]_????.txt (extracted values on local build points on each process;
!                since build points are distributed to each process in order, you can
!                do "cat f18_* > fort.18" ("cat f19_* > fort.19") to get results from all build points)
!											
! mpf90 -mcmodel=medium -CB -O2 -o read_output9_xyzt_mpi.exe ../UtilLib/extract_mod.f90 ../UtilLib/compute_zcor.f90 ../UtilLib/pt_in_poly_test.f90 read_output9_xyzt_mpi.F90 -I$NETCDF/include  -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf  -lnetcdff
!****************************************************************************************
!
      program read_out
      use netcdf
      use extract_mod
      use compute_zcor
      use pt_in_poly_test
      implicit none

      include 'mpif.h'


      integer :: myrank,myrank2,errcode,color,comm,mycomm,itmp,ierr,nm(4),i_bp_glb,nxy_glb,npr,nxy_init
      integer, allocatable :: imap(:), iep_glb(:),node3_glb(:,:)
      real, allocatable :: arco_glb(:,:)
      character(len=4) :: fgb
      character(len=140) :: fname

      character(len=30) :: file63,varname
      character(len=12) :: it_char
      integer,allocatable :: kbp(:),iday(:,:),irecord(:,:),node3(:,:),iep(:),irank_read(:)
      real,allocatable ::  sigma(:),cs(:),ztot(:),times(:,:),out(:,:,:),out2(:,:,:), &
    &eta2(:),arco(:,:),ztmp(:),x00(:),y00(:),z00(:),t00(:),sigma_lcl(:,:),ztmp2(:,:), &
    &outvar(:,:,:)
      real*4,allocatable :: out5(:),out_glb(:)
      real*8,allocatable :: timeout(:)

      integer :: nodel(3) !,dimids(100),idims(100), &
!     &start_2d(2),start_3d(3),start_4d(4), &
!     &count_2d(2),count_3d(3),count_4d(4)
      
      real :: etal,dep,rat,trat,h_s,h_c,theta_b,theta_f,len_var,rjunk
      integer :: idry,irec,irank,i23d,ivs,j,nd,k,kbpl,idry2,kin,k0,m,i,i_read_iep,n
      integer :: l,ntmp,inside,icomb,nxy,istat,iabort,last_dim,ivcor,ie,ibp,iglb

      real :: swild(3),out3(2,2),out4(2)
      !long int for large files
      !integer(kind=8) :: irec

      !MPI
      call MPI_INIT(errcode)
      call mpi_comm_dup(MPI_COMM_WORLD,comm,errcode)
      call mpi_comm_size(comm,npr,ierr)
      call MPI_COMM_RANK(comm, myrank, errcode)

      write(fgb,'(i4.4)') myrank
      fname='out_'//fgb//'.txt'
      open(99,file=fname,status='replace')
      fname='f18_'//fgb//'.txt'
      open(98,file=fname,status='replace')
      fname='f19_'//fgb//'.txt'
      open(97,file=fname,status='replace')


      if (myrank==0) then 
        print*, 'Do you work on uncombined (0) or combined (1) nc?'
        read(*,*)icomb
        if(icomb/=0.and.icomb/=1) stop 'Unknown icomb'
        !icomb=0

        print*, 'Input variable name to read from (e.g. elev):'
        read(*,'(a30)')varname
        !varname='temp'
        
  !     Invliad number used for 3D variables: below bottom; dry spot; no parents
        print*, 'Input values to be used for invalid place:'
        read(*,*)rjunk
        !rjunk=-9999

  !     If read iep
        !print*, 'Read enclosing elements from iep:'
        !read(*,*)i_read_iep
        i_read_iep=0
      endif ! Rank 0 reads cmd line inputs

  !     Distribute to other ranks
      call mpi_bcast(icomb,1,MPI_INTEGER,0,comm,istat)
      call mpi_bcast(varname,30,MPI_CHARACTER,0,comm,istat)
      call mpi_bcast(rjunk,1,MPI_REAL4,0,comm,istat)
      call mpi_bcast(i_read_iep,1,MPI_INTEGER,0,comm,istat)

      varname=adjustl(varname); len_var=len_trim(varname)
      write(99,*) 'icomb, varname, rjunk, i_read_iep: ', &
      &icomb, varname, rjunk, i_read_iep

      open(10,file='station.xyzt',status='old')
      read(10,*) 
      read(10,*) nxy_glb
      allocate(x00(nxy_glb),y00(nxy_glb),z00(nxy_glb),t00(nxy_glb),stat=istat)
      if(istat/=0) stop 'Falied to allocate (1)'

      do i=1,nxy_glb
        read(10,*)j,x00(i),y00(i),z00(i),t00(i) !z00>=0 from F.S.; t00 in sec
        if(z00(i)<0) then
          write(*,*)'Invalid z value:',i; stop
        endif
      enddo !i
      close(10)

      !Distribute bp points to each rank
      nxy_init=ceiling(float(nxy_glb)/npr)
      allocate(imap(nxy_init))
      nxy=0
      do i=1,nxy_init
        i_bp_glb=myrank*nxy_init+i
        if (i_bp_glb<=nxy_glb) then
          imap(i)=i_bp_glb
          nxy=nxy+1
        else
          exit
        endif
      enddo

      write(99,*) 'Local bp: ', nxy
      write(99,*) imap(1:nxy)

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
      allocate(timeout(nrec),out(3,nvrt,2),out2(2,nvrt,2),out5(2*5*nxy_init),eta2(np),node3(nxy,3),arco(nxy,3), &
    &iep(nxy),iday(2,nxy),irecord(2,nxy),times(2,nxy),outvar(2,nvrt,last_dim), &
    &iep_glb(nxy_glb),node3_glb(nxy_glb,3),arco_glb(nxy_glb,3),stat=istat)
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
      iep_glb=0
!      arco=1./3 !initialize for pts without parents
!      do l=1,nxy
!        node3(l,1:3)=elnode(1:3,1) !initialize for pts without parents
!      enddo !l
      iep=0

      if (i_read_iep==0) then
        do i=1,nxy
          iglb=imap(i)
          do ie=1,ne
            call pt_in_poly_single(i34(ie),x(elnode(1:i34(ie),ie)),y(elnode(1:i34(ie),ie)),x00(iglb),y00(iglb),inside,arco(i,1:3),nodel)
            if(inside==1) then
              iep(i)=ie
              print*, 'Found:',i,arco(i,1:3),nodel
              node3(i,1:3)=elnode(nodel(1:3),ie)
              exit
            endif !inside
          enddo !i=1,ne
        enddo !i; local build pts
      elseif (i_read_iep==1) then
  !     read enclosing element ids
        open(10,file='iep',status='old')
        read(10,*) ntmp
        if (myrank==0) then
          if (ntmp.ne.nxy_glb) then
            print*, 'nxy in iep and station.xyzt not consistent'
            stop
          else
            print*, 'nxy in iep and station.xyzt consistent'
          endif
        endif
        do i=1,nxy_glb
          read(10,*) iep_glb(i),arco_glb(i,1:3),node3_glb(i,1:3)
        enddo !i
        close(10)
  !     convert to local
        do i=1,nxy
          iep(i)=iep_glb(imap(i))
          arco(i,1:3)=arco_glb(imap(i),1:3)
          node3(i,1:3)=node3_glb(imap(i),1:3)
        enddo
      else
        write(11,*) 'wrong i_read_iep'
        stop
      endif

      do i=1,nxy
        write(99,'(I10,3F15.10,3I10)') iep(i),arco(i,1:3),node3(i,1:3)
      enddo

!      open(10,file='iep',status='replace')
!      write(10,*) nxy
!      do l=1,nxy
!        write(10,'(I10,3F15.10,3I10)') iep(l),arco(l,1:3),node3(l,1:3)
!      enddo
!      close(10)
!
      iabort=0
      do i=1,nxy
        iglb=imap(i)
        if(iep(i)<=0) then
          write(11,*)'Cannot find a parent for pt:',iglb,i,x00(iglb),y00(iglb)
          iabort=1
        endif
      enddo !i
      if(iabort==1) stop 'check fort.11 for pts outside'

      !Mark ranks that need to be read in
      if(icomb==0) then
        allocate(irank_read(0:nproc-1))
        irank_read=0
        do i=1,nxy
          ie=iep(i)
          irank_read(iegl_rank(ie))=1
          write(99,*)'reading from rank #:',iegl_rank(ie)
        enddo
      endif !icomb

!...  Compute stack and record # for each pt
      do i=1,nxy
        iglb=imap(i)
!       Check if time is before first record
        if(t00(iglb)<dtout) then
          write(11,*)'Time before first record:',iglb,i,t00(iglb)
          stop
        endif

!       Lower and upper bound stacks and record #s for t00
        iday(1,i)=(t00(iglb)-dtout)/nrec/dtout+1
        if(iday(1,i)<1) then
          write(11,*)'Impossible'; stop
        else
          irecord(1,i)=(t00(iglb)-(iday(1,i)-1)*nrec*dtout)/dtout
          !Bounding record time just b4 the cast time, corresponding to record
          !irecord(1,i) in stack iday(1,i)
          times(1,i)=((iday(1,i)-1)*nrec+irecord(1,i))*dtout
          iday(2,i)=t00(iglb)/nrec/dtout+1
          irecord(2,i)=(t00(iglb)-(iday(2,i)-1)*nrec*dtout)/dtout+1
          !Bounding record time just after the cast time, corresponding to
          !record irecord(2,i) in stack iday(2,i). Note that irecord(2,i)
          !may<irecord(1,i) (e.g. t00 before 1st record of iday(2,i))
          times(2,i)=((iday(2,i)-1)*nrec+irecord(2,i))*dtout
        endif

        if(irecord(1,i)>nrec.or.irecord(2,i)>nrec) then
          write(11,*)'Record # overflow: ',iglb,i,irecord(:,i)
          stop
        endif
        if(t00(iglb)<times(1,i).or.t00(iglb)>times(2,i)) then
          write(11,*)'Wrong time bounds:',iglb,i,t00(iglb),times(:,i),iday(:,i),irecord(:,i)
          stop
        endif
      enddo !i=1,nxy

!...  Time iteration
      do i=1,nxy
        iglb=imap(i)
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
!          call get_outvar(1,iday(l,i),irec,varname,np,last_dim,nvrt,outvar,i23d,ivs,eta2)          
          if(icomb==0) then !uncombined
            do irank=0,nproc-1
              if(irank_read(irank)>0) then
                call get_outvar(1,iday(l,i),irec,varname,np,last_dim,nvrt,outvar,i23d,ivs,eta2,irank)
              endif
            enddo !irank
          else
            call get_outvar(1,iday(l,i),irec,varname,np,last_dim,nvrt,outvar,i23d,ivs,eta2)
          endif

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
                if(ztmp(nvrt)-z00(iglb)>=ztmp(k).and.ztmp(nvrt)-z00(iglb)<=ztmp(k+1)) then
                  k0=k
                  rat=(ztmp(nvrt)-z00(iglb)-ztmp(k))/(ztmp(k+1)-ztmp(k))
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
        trat=(t00(iglb)-times(1,i))/(times(2,i)-times(1,i)) !must be [0,1]
        if(mod(i23d-1,3)==0) then
          if(iep(i)==0) then !no parents
            out4(1:ivs)=rjunk
          else
            out4(1:ivs)=out2(1,1,1:ivs)*(1-trat)+out2(2,1,1:ivs)*trat
          endif
          write(98,'(e16.8,1000(1x,f15.3))')t00(iglb)/86400,out4(1:ivs),x00(iglb),y00(iglb)
          out5((i-1)*10+1)=t00(iglb)/86400
          out5((i-1)*10+2)=out4(1)
          out5((i-1)*10+3)=x00(iglb)
          out5((i-1)*10+4)=y00(iglb)
          out5((i-1)*10+5)=z00(iglb);
        else !3D
          if(iep(i)==0) then !no parents
            out4(1:ivs)=rjunk
          else
            out4(1:ivs)=out3(1,1:ivs)*(1-trat)+out3(2,1:ivs)*trat
          endif
          write(98,'(f20.10,1000(f20.10))')t00(iglb)/86400,out4(1),x00(iglb),y00(iglb),z00(iglb)
          out5((i-1)*10+1)=t00(iglb)/86400
          out5((i-1)*10+2)=out4(1)
          out5((i-1)*10+3)=x00(iglb)
          out5((i-1)*10+4)=y00(iglb)
          out5((i-1)*10+5)=z00(iglb);
          if(ivs==2) then
            write(97,'(e16.8,1000(1x,f15.3))')t00(iglb)/86400,out4(2),x00(iglb),y00(iglb),z00(iglb)
            out5((i-1)*10+5+1)=t00(iglb)/86400
            out5((i-1)*10+5+2)=out4(1)
            out5((i-1)*10+5+3)=x00(iglb)
            out5((i-1)*10+5+4)=y00(iglb)
            out5((i-1)*10+5+5)=z00(iglb);
          endif
        endif
      enddo !i=1,nxy

      !if (myrank==0) then
      !  allocate(out_glb(2*5*nxy_glb))
      !endif
      !call mpi_gather(out5,nxy_init*10,MPI_REAL4, out_glb,nxy_init*10,MPI_REAL4, 0,MPI_COMM_WORLD)
      !if (myrank==0) then
      !  do i=1,nxy_glb
      !    write(28,*) out_glb((i-1)*10+1:(i-1)*10+5)
      !    if(ivs==2) write(29,*) out_glb((i-1)*10+5+1:(i-1)*10+5+5)
      !  enddo !i=1,nxy_glb
      !endif

      print*, 'Finished!'

      call MPI_FINALIZE(errcode)
      end program
