!********************************************************************************
!										
!	Interpolate 2 or 3D variables from SCHISM nc outputs to get *3D.th.nc 
!       Use generic method to find parents.
!       Changes from interpolate_variables_selfe2: abandoned bucket sort; vgrid.fg can be SZ now.
!       Changes from interpolate_variables_selfe3: *3D.th now binary; added ivcor=1.
!       Changes from interpolate_variables_selfe4: added 1 extra record at the beginning for new .th format
!       Changes from interpolate_variables_selfe5: add quads
!       Changes from interpolate_variables6: changed to nc
!										
!       Inputs: 
!          (1) bg.gr3 (for more precise x,y): hgrid.gr3 from large-domain run;
! 	   (2) fg.gr3: hgrid.gr3 from small-domain run; boundary segments 
!                      (specified in interpolate_variables.in) that need *3D.th 
!                      may be lumped; make sure that the parent elem. in the large-domain run
!                      never becomes dry!
!          (3) vgrid.bg: background vgrid.in
!	   (4) vgrid.fg: vgrid.in from the small-domain run
!          (5) interpolate_variables.in (see sample in this dir):
!              1st line: ifile rndays - ifile=1: generate elev2D.th; =2: SAL_3D.th 
!                        and TEM_3D.th; =3: uv3D.th); rndays is the # of days needed;
!              2nd line: total # of open bnd seg. (that need *3D.th in fg.gr3), list of segments IDs
!              3rd line: idebug - more outputs for debug if idebug=1
!          (6) combined nc outputs from large-domain run (must have elev or temp or salt or hvel depending on ifile)
!       Outputs: 
!          (1) *[23]D.th.nc, depending on the choice in interpolate_variables.in
!          (2) fort.11: fatal errors; fort.12: non-fatal errors.
!
! ifort -O2 -mcmodel=medium -CB -g -traceback -o interpolate_variables7.exe ../UtilLib/extract_mod.f90 \
!../UtilLib/compute_zcor.f90 ../UtilLib/pt_in_poly_test.f90 interpolate_variables7.f90 \
!-I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff
!********************************************************************************
!
      program interpolate
      use netcdf
      use extract_mod
      use compute_zcor
      use pt_in_poly_test

!      implicit real*8(a-h,o-z)
      character(len=30) :: file63,varname(2),varname2
      character(len=12) :: it_char
      !bg arrays: whole grid
      real,allocatable :: area(:),eta2(:),tnd(:,:),snd(:,:), &
     &ztot(:),sigma2(:),sigma_lcl(:,:),z(:,:),ze(:),x_fg(:),y_fg(:), &
     &dp_fg(:),ztot_fg(:),sigma_fg(:),sigma_lcl_fg(:,:),cs_fg(:)
      integer,allocatable :: kbp(:),kbp_fg(:)
      !fg arrays: open bnd part only
      real,allocatable :: xfg(:),yfg(:),dpfg(:),zfg(:,:),etafg(:),ratmin(:), &
     &tfg(:,:),sfg(:,:),arco(:,:),outvar(:,:,:),outvar0(:,:,:)
      integer,allocatable :: kbpfg(:),imap(:),iob(:),nond(:),iond(:,:),iparen(:),iparen_b(:)
      integer :: nwild(3),one_dim
      !long int for large files
!      integer(kind=8) :: irec,long_rec
      integer :: var1d_dims(1),var4d_dims(4) !dimids(100),idims(100), &
!     &start_2d(2),start_3d(3),start_4d(4), &
!     &var1d_dims(1),var4d_dims(4)
      real*8 :: aa1(1)
      real*8,allocatable :: timeout(:)

      open(21,file='interpolate_variables.in',status='old')
      read(21,*)ifile,rndays
      read(21,*)nob
      allocate(iob(nob))
      rewind(21)
      read(21,*)
      read(21,*)nob,iob(:) !points to fg grid
      read(21,*)idebug
      close(21)

!..   Read bg.gr3 for more precise x,y
      open(14,file='bg.gr3',status='old')
      read(14,*)
      read(14,*)ne,np
      allocate(x(np),y(np),dp(np),i34(ne),elnode(4,ne),kbp(np),kbp00(np),area(ne), &
     &eta2(np),stat=istat)
      if(istat/=0) stop 'Failed to alloc (1)'
      do i=1,np
        read(14,*)j,x(i),y(i),dp(i)
      enddo !i
      close(14)

!...  Read header 
      if(ifile==1) then
        varname(1)=adjustl('elev')
        n_vars=1
      else if(ifile==2) then
        varname(1)=adjustl('temp')
        varname(2)=adjustl('salt')
        n_vars=2
      else if(ifile==3) then
        varname(1)=adjustl('hvel')
        n_vars=1
      else 
        stop 'unknown ifile'
      endif
      len_var=4
!      step_nu=900

!...  Header
      call readheader(1) !'schout_1.nc')
      !Returned vars: ne,np,ns,nrec,[x y dp kbp00](np),
      !elnode,i34,nvrt,h0,dtout

      print*, 'After header:',ne,np,nrec,i34(ne), &
     &elnode(1:i34(ne),ne),nvrt,h0,x(np),y(np),dp(np),kbp00(np) !,start_time

!     Read in vgrid.in from bg grid
      allocate(timeout(nrec),ztot(nvrt),sigma2(nvrt),sigma_lcl(nvrt,np),z(np,nvrt), &
     &tnd(np,nvrt),snd(np,nvrt),ze(nvrt),stat=istat)
      if(istat/=0) stop 'Failed to alloc (2)'

      if(nvrt==2) then !2D
        ivcor=2; kz=1; nsig=2; h_s=1.e6; h_c=h_s 
        theta_b=0; theta_f=1.e-4;
        do i=1,np
          if(dp(i)<=h0) then
            kbp(i)=0
          else
            kbp(i)=1
            z(i,1)=-dp(i)
            z(i,nvrt)=0
          endif
        enddo !i
      else !3D
        call get_vgrid_single('vgrid.bg',np,nvrt,ivcor,kz,h_s,h_c,theta_b,theta_f,ztot,sigma2,sigma_lcl,kbp)
        do i=1,np
          if(ivcor==2) then
            call zcor_SZ_single(dp(i),0.,h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot,sigma2,z(i,:),idry,kbp(i))
          else if(ivcor==1) then
            if(dp(i)<=h0) then
              kbp(i)=0
            else
              z(i,kbp(i):nvrt)=dp(i)*sigma_lcl(kbp(i):nvrt,i)
            endif
          else
            write(11,*)'Unknown ivcor:',ivcor
            stop
          endif
        enddo !i
      endif !2/3D

      !Check
      do i=1,np
        if(kbp(i)/=0) then
          do k=kbp(i)+1,nvrt
            if(z(i,k)<=z(i,k-1)) then
              write(11,*)'Inverted z for large grid:',i,k,z(i,kbp(i):nvrt)
              stop
            endif
          enddo !k

          if(idebug==1) write(95,*)'zcor for bg:',i,kbp(i),z(i,kbp(i):nvrt)
        endif !kbp
      enddo !i

!     Read in fg.gr3
      open(13,file='fg.gr3',status='old')
      read(13,*)
      read(13,*)ne_fg,np_fg
      allocate(x_fg(np_fg),y_fg(np_fg),dp_fg(np_fg))
      do i=1,np_fg
        read(13,*)j,x_fg(i),y_fg(i),dp_fg(i)
!        if(dp_fg(i)<=0) then
!          write(*,*)'Fg pt has negative depth:',i
!          stop
!        endif
      enddo !i
      do i=1,ne_fg; read(13,*); enddo
      read(13,*)nope
      read(13,*)neta
      allocate(nond(nope))
      do i=1,nope
        read(13,*)nond(i)
        do j=1,nond(i)
          read(13,*) !iond
        enddo !j
      enddo !i
      rewind(13)
      mnope=maxval(nond)
      allocate(iond(nope,mnope))
      do i=1,np_fg+ne_fg+2+2; read(13,*); enddo
      do i=1,nope
        read(13,*) !nond(i)
        do j=1,nond(i)
          read(13,*) iond(i,j)
        enddo !j
      enddo !i
      close(13)

      npfg=sum(nond(iob(1:nob))) !# of b.c. nodes
      print*, '# of relevant open bnd nodes=',npfg

      allocate(imap(npfg),xfg(npfg),yfg(npfg),dpfg(npfg),etafg(npfg),ratmin(npfg), &
     &iparen(npfg),iparen_b(npfg),arco(npfg,4))

      npfg=0
      do i=1,nob
        ibnd=iob(i)
        do j=1,nond(ibnd)
          npfg=npfg+1
          nd=iond(ibnd,j)
          imap(npfg)=nd
          xfg(npfg)=x_fg(nd)
          yfg(npfg)=y_fg(nd)
          dpfg(npfg)=dp_fg(nd)

          if(idebug==1) write(95,*)'List of fg nodes:',npfg,imap(npfg)
        enddo !j 
      enddo !i

!     Read in vgrid.fg 
      open(19,file='vgrid.fg',status='old')
      read(19,*); read(19,*)nvrt_fg
      close(19)

      last_dim=np
      allocate(ztot_fg(nvrt_fg),sigma_fg(nvrt_fg),sigma_lcl_fg(nvrt_fg,np_fg),kbp_fg(np_fg), &
     &zfg(npfg,nvrt_fg),kbpfg(npfg),tfg(npfg,nvrt_fg),sfg(npfg,nvrt_fg),cs_fg(nvrt_fg), &
     &outvar0(2,nvrt,np),outvar(2,nvrt,np))

      if(nvrt_fg==2) then !2D
!        ivcor=2; kz_fg=1; nsig_fg=2; h_s_fg=1.e6; h_c_fg=h_s
!        theta_b_fg=0; theta_f_fg=1.e-4;
        do i=1,npfg
          !Use h0 from bg run
          if(dpfg(i)<=h0) then
            write(11,*)'Small-domain run depth too small:',imap(i),h0
            stop
            !kbpfg(i)=0
          else
            kbpfg(i)=1
            zfg(i,1)=-dpfg(i)
            zfg(i,nvrt_fg)=0
          endif
        enddo !i
      else !3D
        call get_vgrid_single('vgrid.fg',np_fg,nvrt_fg,ivcor_fg,kz_fg,h_s_fg,h_c_fg, &
     &theta_b_fg,theta_f_fg,ztot_fg,sigma_fg,sigma_lcl_fg,kbp_fg)

        !zcor
        do i=1,npfg
          if(dpfg(i)<=h0) then
            write(11,*)'Small-domain run depth too small:',imap(i),h0
            stop
          endif

          if(ivcor_fg==2) then !SZ
            !Use h0 from bg run
            call zcor_SZ_single(dpfg(i),0.,h0,h_s_fg,h_c_fg,theta_b_fg,theta_f_fg,kz_fg, &
     &nvrt_fg,ztot_fg,sigma_fg,zfg(i,:),idry2,kbpfg(i))
          else !=1
            nd=imap(i)
            kbpfg(i)=kbp_fg(nd)
            zfg(i,kbpfg(i):nvrt_fg)=dpfg(i)*sigma_lcl_fg(kbpfg(i):nvrt_fg,nd)
          endif !ivcor_fg
          !Extend 
          zfg(i,1:kbpfg(i)-1)=zfg(i,kbpfg(i))
        enddo !i=1,npfg
      endif !nvrt_fg

      !Check z-cor
      do i=1,npfg
        do k=kbpfg(i)+1,nvrt_fg
          if(zfg(i,k)<=zfg(i,k-1)) then
            write(11,*)'Inverted z-cor in small grid:',imap(i),k,zfg(i,1:nvrt_fg)
            stop
          endif
        enddo !k

        if(idebug==1) write(95,*)'fg zcor:',i,imap(i),kbpfg(i),zfg(i,:)
      enddo !i

!...  Find parent elements - split quads
      arco=0
      iparen=0 !flags
      iparen_b=0 !flags; back-up parent based on min ratio
      ratmin=1.e25 !min. area ratio for fg pt
      loop1: do ie=1,ne
        iexit=1 !flag
        do i=1,npfg
          if(iparen(i)/=0) cycle

          iexit=0

          do j=1,i34(ie)-2
            if(j==1) then
              nwild(1:3)=(/1,2,3/) !pt to local indices
            else !quads
              nwild(1:3)=(/1,3,4/)
            endif !j
            n1=elnode(nwild(1),ie); n2=elnode(nwild(2),ie); n3=elnode(nwild(3),ie)
            ar1=signa_single(xfg(i),x(n2),x(n3),yfg(i),y(n2),y(n3))  
            ar2=signa_single(x(n1),xfg(i),x(n3),y(n1),yfg(i),y(n3))
            ar3=signa_single(x(n1),x(n2),xfg(i),y(n1),y(n2),yfg(i))
            bb=abs(ar1)+abs(ar2)+abs(ar3)
            aa=abs(signa_single(x(n1),x(n2),x(n3),y(n1),y(n2),y(n3)))
            rat=abs(bb-aa)/aa
            if(rat<ratmin(i)) then
              ratmin(i)=rat
              iparen_b(i)=ie
            endif
            if(rat<1.e-4) then
              iparen(i)=ie
              arco(i,nwild(1))=max(0.,min(1.,ar1/aa))
              arco(i,nwild(2))=max(0.,min(1.,ar2/aa))
              if(arco(i,nwild(1))+arco(i,nwild(2))>1) then
                arco(i,nwild(3))=0
                arco(i,nwild(2))=1-arco(i,nwild(1))
              else
                arco(i,nwild(3))=1-arco(i,nwild(1))-arco(i,nwild(2))
              endif
          
              !1 more acro for quads
              if(i34(ie)==4) then
                itmp=1+2+3+4-sum(nwild(1:3))
                arco(i,itmp)=0
              endif
            endif
          enddo !j=1,i34(ie)-2
        enddo !i=1,npfg
        if(iexit==1) exit
      end do loop1 !ie
      
!     Back-up parents
      icount=0
      write(12,*)'Following pts have no parents:'
      do i=1,npfg
        if(iparen(i)==0) then
          icount=icount+1
          if(iparen_b(i)==0) then
            write(11,*) 'No back-up:',i
            stop
          endif
          write(12,*)i,iparen_b(i)
          iparen(i)=iparen_b(i)
          arco(i,1:i34(iparen(i)))=1.0/i34(iparen(i))
        endif

        if(idebug==1) write(95,*)'Parents elem:',i,imap(i),iparen(i),arco(i,1:i34(iparen(i)))
      enddo !i
      print*, icount,' pts have no immediate parent elements; see fort.12'
!'

!...  Time iteration
!...
!     Open outputs
      if(ifile==1) then
        iret=nf90_create('elev2D.th.nc',OR(NF90_NETCDF4,NF90_CLOBBER),ncid)
        nv=1; ivs=1
      else if(ifile==2) then !T,S
        iret=nf90_create('TEM_3D.th.nc',OR(NF90_NETCDF4,NF90_CLOBBER),ncid)
        iret=nf90_create('SAL_3D.th.nc',OR(NF90_NETCDF4,NF90_CLOBBER),ncid2)
        nv=nvrt_fg; ivs=1
      else !uv
        iret=nf90_create('uv3D.th.nc',OR(NF90_NETCDF4,NF90_CLOBBER),ncid)
        nv=nvrt_fg; ivs=2
      endif

      iret=nf90_def_dim(ncid,'nOpenBndNodes',npfg,nond0_dim)
      iret=nf90_def_dim(ncid,'nLevels',nv,nvrt_dim)
      iret=nf90_def_dim(ncid,'nComponents',ivs,ivs_dim)
      iret=nf90_def_dim(ncid,'one',1,one_dim)
      iret=nf90_def_dim(ncid,'time', NF90_UNLIMITED,itime_dim)

      var1d_dims(1)=one_dim
      iret=nf90_def_var(ncid,'time_step',NF90_FLOAT,var1d_dims,idt)
      var1d_dims(1)=itime_dim
      iret=nf90_def_var(ncid,'time',NF90_DOUBLE,var1d_dims,itimeout_id)
      var4d_dims(1)=ivs_dim; var4d_dims(2)=nvrt_dim; var4d_dims(3)=nond0_dim
      var4d_dims(4)=itime_dim
      iret=nf90_def_var(ncid,'time_series',NF90_FLOAT,var4d_dims,ivarid)
      iret=nf90_enddef(ncid)
      iret=nf90_put_var(ncid,idt,dtout)
 
      if(ifile==2) then
        iret=nf90_def_dim(ncid2,'nOpenBndNodes',npfg,nond0_dim)
        iret=nf90_def_dim(ncid2,'nLevels',nv,nvrt_dim)
        iret=nf90_def_dim(ncid2,'nComponents',ivs,ivs_dim)
        iret=nf90_def_dim(ncid2,'one',1,one_dim)
        iret=nf90_def_dim(ncid2,'time', NF90_UNLIMITED,itime_dim)

        var1d_dims(1)=one_dim
        iret=nf90_def_var(ncid2,'time_step',NF90_FLOAT,var1d_dims,idt)
        var1d_dims(1)=itime_dim
        iret=nf90_def_var(ncid2,'time',NF90_DOUBLE,var1d_dims,itimeout_id2)
        var4d_dims(1)=ivs_dim; var4d_dims(2)=nvrt_dim; var4d_dims(3)=nond0_dim
        var4d_dims(4)=itime_dim
        iret=nf90_def_var(ncid2,'time_series',NF90_FLOAT,var4d_dims,ivarid2)
        iret=nf90_enddef(ncid2)
        iret=nf90_put_var(ncid2,idt,dtout)
      endif !ifile

      irec_out=0
      it_tot=0
      e_max=-huge(1.0) !max of elev. over all nodes and time
      e_min=huge(1.0) !min of elev.
      u_max=-huge(1.0) !max of u|T
      v_max=-huge(1.0) !max of v|S
      u_min=huge(1.0) !min of u|T
      v_min=huge(1.0) !min of v|S
      do iday=1,max(1,int(rndays*86400/dtout/nrec+0.1)) !1+int(rndays*86400/dtout/nrec+0.1)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!      write(it_char,'(i12)')iday
!      it_char=adjustl(it_char)
!      leng=len_trim(it_char)
!      file63='schout_'//it_char(1:leng)//'.nc'
!      iret=nf90_open(trim(adjustl(file63)),OR(NF90_NETCDF4,NF90_NOWRITE),ncid3)
!      !time is double
!      iret=nf90_inq_varid(ncid3,'time',itime_id)
!      iret=nf90_get_var(ncid3,itime_id,timeout,(/1/),(/nrec/))
      call get_timeout(iday,nrec,timeout)
      !write(98,*)'time=',timeout,trim(adjustl(file63))

      do it1=1,nrec
        !Get elev
!        start_2d(1)=1; start_2d(2)=it1
!        count_2d(1)=np; count_2d(2)=1
!        iret=nf90_get_var(ncid3,ielev_id,eta2,start_2d,count_2d)
        call get_elev(iday,it1,np,eta2)

        if(ifile>1) then !get 3D vars
          do ii=1,n_vars
            if(ifile==2) then !TS
              call get_outvar(1,iday,it1,varname(ii),np,last_dim,nvrt,outvar0,i23d,ivs,eta2)
              outvar(ii,:,:)=outvar0(1,:,:)
            else if(ifile==3) then !uv
              call get_outvar(1,iday,it1,varname(ii),np,last_dim,nvrt,outvar,i23d,ivs,eta2)
            endif !

!            varname2=adjustl(varname(ii))
!            iret=nf90_inq_varid(ncid3,varname2(1:len_var),ivarid1)
!            if(iret/=nf90_NoErr) stop 'Var not found'
!            iret=nf90_Inquire_Variable(ncid3,ivarid1,ndims=ndims,dimids=dimids)
!            if(ndims>100) stop 'increase dimension of dimids & idims'
!            do i=1,ndims
!              iret=nf90_Inquire_Dimension(ncid3,dimids(i),len=idims(i))
!            enddo !i
!            if(idims(ndims-1)/=np) stop 'can only handle node-based'
!            if(idims(ndims)/=nrec) stop 'last dim is not time'
!
!            iret=nf90_get_att(ncid3,ivarid1,'i23d',i23d)
!            iret=nf90_get_att(ncid3,ivarid1,'ivs',ivs)
!            !print*, 'i23d:',i23d,ivs,idims(1:ndims)
!
!            if(ifile==2) then !TS
!              start_3d(1:2)=1; start_3d(3)=it1
!              count_3d(2)=np; count_3d(1)=nvrt; count_3d(3)=1
!              iret=nf90_get_var(ncid3,ivarid1,outvar(ii,:,:),start_3d,count_3d)
!            else if(ifile==3) then !uv
!              start_4d(1:3)=1; start_4d(4)=it1
!              count_4d(3)=np; count_4d(2)=nvrt; count_4d(1)=2; count_4d(4)=1
!              iret=nf90_get_var(ncid3,ivarid1,outvar(:,:,:),start_4d,count_4d)
!            endif !ivs
          enddo !ii

          !Available now: outvar(2,nvrt,np) (for ifile=2, 1:T; 2:S)
        endif !ifile>1

        !print*, 'processing time (days) = ',time/86400,' from stack:',iday

        if(ifile==1) then !elev
!         Do interpolation
          do i=1,npfg
            ie=iparen(i)
            iflag=0
            etafg(i)=0
            do j=1,i34(ie)
              nd=elnode(j,ie)

              !if(eta2(nd)+dp(nd)<=h0) iflag=1
              etafg(i)=etafg(i)+eta2(nd)*arco(i,j)
            enddo !j
            if(etafg(i)+dpfg(i)<=h0) iflag=1
!           Bg dry nodes
            if(iflag==1) then
              write(11,*)'Dry node:',i
              stop
            endif
          enddo !i

          if(idebug==1) write(23,'(10000(1x,e15.7))')time,etafg(1:npfg)
          
          !Add a record at t=0
          if(iday==1.and.it1==1) then
            irec_out=irec_out+1
            aa1(1)=0
            iret=nf90_put_var(ncid,itimeout_id,aa1,(/irec_out/),(/1/))
            iret=nf90_put_var(ncid,ivarid,etafg(1:npfg),(/1,1,1,irec_out/),(/1,1,npfg,1/))
!            write(16,rec=irec_out+1)0.,(etafg(i),i=1,npfg)
          endif
          irec_out=irec_out+1
          aa1(1)=timeout(it1)
          iret=nf90_put_var(ncid,itimeout_id,aa1,(/irec_out/),(/1/))
!          write(16,rec=irec_out+1)time,(etafg(i),i=1,npfg)
          iret=nf90_put_var(ncid,ivarid,etafg(1:npfg),(/1,1,1,irec_out/),(/1,1,npfg,1/))

          e_max=max(e_max,maxval(etafg(1:npfg)))
          e_min=min(e_min,minval(etafg(1:npfg)))
        else !T,S or uv
!         Do interpolation (no true vertical interpolation)
!          print*, 'Day ',iday
          do i=1,npfg
            ie=iparen(i)
            kbe=maxval(kbp(elnode(1:i34(ie),ie)))
            if(kbe==0) then
              write(11,*)'All dry'
              stop
            endif
            do k=kbe,nvrt
              ze(k)=-1.e6 !take the max. of z-cor from all nodes
              do j=1,i34(ie)
                nd=elnode(j,ie)
                if(kbp(nd)/=0.and.z(nd,k)>ze(k)) ze(k)=z(nd,k)
              enddo !j
            enddo !k

            !Debug
            !write(97,*)'Time=',time/86400,i,ie,kbe,ze(kbe:nvrt)

            do k=1,nvrt_fg
              if(zfg(i,k)<=ze(kbe)) then
                klev=kbe
              else if(zfg(i,k)>=ze(nvrt)) then
                klev=nvrt
              else
                klev=0
                do kk=kbe,nvrt-1
                  if(zfg(i,k)>=ze(kk).and.zfg(i,k)<=ze(kk+1)) then
                    klev=kk
                    exit
                  endif
                enddo !kk
                if(klev==0) then
                  !ze() may not be in proper order to cause this
                  write(11,*)'Cannot find a level:',i,k,zfg(i,k),ze(kbe:nvrt)
                  stop
                endif
              endif

              tfg(i,k)=0 !for u if ifile=3
              sfg(i,k)=0 !for v if ifile=3
              iflag=0
              tanc=-99; sanc=-99 !back-up values
              do j=1,i34(ie)
                nd=elnode(j,ie)
                 
                if(ifile==2.and.outvar(2,klev,nd)<0) then
                  iflag=1
                else
                  tanc=outvar(1,klev,nd) !tnd(nd,klev)
                  sanc=outvar(2,klev,nd) !snd(nd,klev)
                endif
                tfg(i,k)=tfg(i,k)+outvar(1,klev,nd)*arco(i,j)
                sfg(i,k)=sfg(i,k)+outvar(2,klev,nd)*arco(i,j)
              enddo !j
!             Bg dry nodes
              if(iflag==1) then !ifile=2
                if(sanc<0) then 
                  write(11,*)'Parent element dry:',i
                  stop
                endif
                tfg(i,k)=tanc
                sfg(i,k)=sanc
              endif

!             Debug
!              write(97,*)'bg:',k,klev,tnd(elnode(:,ie),klev),snd(elnode(:,ie),klev)
!              write(97,*)'fg:',arco(i,:),tfg(i,k),sfg(i,k)

            enddo !k=1,nvrt_fg
!            if(i==1) then
!              write(98,*)kbe,(ze(k),k=kbe,nvrt)
!              write(98,*)(zfg(i,k),k=1,nvrt_fg)
!            endif
          enddo !i=1,npfg

!         Output
          if(ifile==2) then !T,S
            !Add a record at t=0
            if(iday==1.and.it1==1) then
              irec_out=irec_out+1
              aa1(1)=0
              iret=nf90_put_var(ncid,itimeout_id,aa1,(/irec_out/),(/1/))
              iret=nf90_put_var(ncid,ivarid,transpose(tfg(1:npfg,1:nvrt_fg)), &
     &(/1,1,1,irec_out/),(/1,nvrt_fg,npfg,1/))
              iret=nf90_put_var(ncid2,itimeout_id2,aa1,(/irec_out/),(/1/))
              iret=nf90_put_var(ncid2,ivarid2,transpose(sfg(1:npfg,1:nvrt_fg)), &
     &(/1,1,1,irec_out/),(/1,nvrt_fg,npfg,1/))
!              write(17,rec=irec_out)0.,((tfg(i,k),k=1,nvrt_fg),i=1,npfg)
!              write(18,rec=irec_out)0.,((sfg(i,k),k=1,nvrt_fg),i=1,npfg)
            endif

            irec_out=irec_out+1
            aa1(1)=timeout(it1)
            iret=nf90_put_var(ncid,itimeout_id,aa1,(/irec_out/),(/1/))
            iret=nf90_put_var(ncid,ivarid,transpose(tfg(1:npfg,1:nvrt_fg)), &
     &(/1,1,1,irec_out/),(/1,nvrt_fg,npfg,1/))
            iret=nf90_put_var(ncid2,itimeout_id2,aa1,(/irec_out/),(/1/))
            iret=nf90_put_var(ncid2,ivarid2,transpose(sfg(1:npfg,1:nvrt_fg)), &
     &(/1,1,1,irec_out/),(/1,nvrt_fg,npfg,1/))
!            write(17,rec=irec_out+1)time,((tfg(i,k),k=1,nvrt_fg),i=1,npfg)
!            write(18,rec=irec_out+1)time,((sfg(i,k),k=1,nvrt_fg),i=1,npfg)
          else !uv
            !Add a record at t=0
            if(iday==1.and.it1==1) then
              irec_out=irec_out+1
              aa1(1)=0
              iret=nf90_put_var(ncid,itimeout_id,aa1,(/irec_out/),(/1/))
              iret=nf90_put_var(ncid,ivarid,transpose(tfg(1:npfg,1:nvrt_fg)), &
     &(/1,1,1,irec_out/),(/1,nvrt_fg,npfg,1/))
              iret=nf90_put_var(ncid,ivarid,transpose(sfg(1:npfg,1:nvrt_fg)), &
     &(/2,1,1,irec_out/),(/1,nvrt_fg,npfg,1/))
!              write(17,rec=irec_out+1)0.,((tfg(i,k),sfg(i,k),k=1,nvrt_fg),i=1,npfg)
            endif
!            write(17,rec=irec_out+1)time,((tfg(i,k),sfg(i,k),k=1,nvrt_fg),i=1,npfg)

            irec_out=irec_out+1
            aa1(1)=timeout(it1)
            iret=nf90_put_var(ncid,itimeout_id,aa1,(/irec_out/),(/1/))
            iret=nf90_put_var(ncid,ivarid,transpose(tfg(1:npfg,1:nvrt_fg)), &
     &(/1,1,1,irec_out/),(/1,nvrt_fg,npfg,1/))
            iret=nf90_put_var(ncid,ivarid,transpose(sfg(1:npfg,1:nvrt_fg)), &
     &(/2,1,1,irec_out/),(/1,nvrt_fg,npfg,1/))

          endif

          u_max=max(u_max,maxval(tfg))
          v_max=max(v_max,maxval(sfg))
          u_min=min(u_min,minval(tfg))
          v_min=min(v_min,minval(sfg))

          if(idebug==1) then !use check_3D.m to get quiver plot
            write(23,'(70000(1x,e16.8))')time,((tfg(i,k),k=1,nvrt_fg),i=1,npfg)
            write(24,'(70000(1x,e16.8))')time,((sfg(i,k),k=1,nvrt_fg),i=1,npfg)
          endif
        endif !i23d
!        print*, 'finished writing time (days) = ',time/86400
      enddo !it1=1,nrec
      iret=nf90_close(ncid3)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      enddo !iday=1,

      !Output x,y for debug plots
      if(idebug==1) write(25,'(2(1x,i4),10000(1x,e16.8))')ifile,nvrt_fg,xfg,yfg

!     Output extrema
      if(ifile==1) then
        print*, 'max/min of elev=',e_max,e_min
      else !2 or 3
        print*, 'max/min of vel or T,S=',u_max,u_min,v_max,v_min
      endif !ifile

      stop
      end

