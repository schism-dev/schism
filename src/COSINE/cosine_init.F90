!Routines & functions
!cosine_init: allocate and initilize variables
!read_cosine_param: read cosine parameters
!read_cosine_stainfo: read info. for station outputs

subroutine cosine_init
!---------------------------------------------------------------------------
!allocate COSINE arrays and initialize
!---------------------------------------------------------------------------
  use schism_glbl, only : nea,npa,nvrt
  use schism_msgp, only : parallel_abort
  use cosine_mod
  implicit none
  
  !local variables
  integer :: istat

  !allocate 
  allocate( NO3(nvrt),NH4(nvrt),SiO4(nvrt),S1(nvrt),S2(nvrt),Z1(nvrt),Z2(nvrt),&
          & DN(nvrt),DSi(nvrt),PO4(nvrt),DOX(nvrt),CO2(nvrt),ALK(nvrt),temp(nvrt),&
          & salt(nvrt),bgraze(nea),SPM(nvrt,nea),bcos(nvrt,ntrc),& 
          & mS2(ndelay,nvrt,nea),mDN(ndelay,nvrt,nea),mZ1(ndelay,nvrt,nea),mZ2(ndelay,nvrt,nea),&
          & sS2(nvrt,nea),sDN(nvrt,nea),sZ1(nvrt,nea),sZ2(nvrt,nea),nstep(nvrt,nea), &
          & nclam(nea), stat=istat) 
  if(istat/=0) call parallel_abort('failure in alloc. mS2')

  !initialize
  NO3=0.0;  NH4=0.0; SiO4=0.0; S1=0.0;   S2=0.0;   Z1=0.0;  Z2=0.0
  DN=0.0;   DSi=0.0; PO4=0.0;  DOX=0.0;  CO2=0.0;  ALK=0.0; temp=0.0
  salt=0.0; SPM=0.0; bcos=0.0  !;  bio0=0.0; qcos=0.0; sqcos=0.0 
  mS2=0.0;  mDN=0.0; mZ1=0.0;  mZ2=0.0; 
  sS2=0.0;  sDN=0.0; sZ1=0.0;  sZ2=0.0; nstep=0
  nclam=0;

  !read cosine parameters
  call read_cosine_param
   
end subroutine cosine_init

subroutine read_cosine_param
!---------------------------------------------------------------------------
!read parameters in cosine.in
!---------------------------------------------------------------------------
  use schism_glbl, only : rkind,npa,nea,ne_global,np_global,ipgl,iegl,elnode,i34, &
 &in_dir,out_dir,len_in_dir,len_out_dir,ihot
  use schism_msgp, only : myrank,parallel_abort
  use cosine_misc, only : read_gr3_prop
  use cosine_mod
  implicit none

  !local variables
  integer :: i,j,k,m,negb,npgb,nd,itmp,itmp1(1),itmp2(1,1),istat
  real(rkind) :: xtmp,ytmp,tSPM,tSPMs(npa),rtmp,rtmp1(1),rtmp2(1,1) 
  integer :: tnclam !,nclams(npa)
  character(len=2) :: stmp
  character(len=100) :: snum 

  !define namelist
  namelist /MARCO/ idelay,ndelay,ibgraze,idapt,alpha_corr,zeptic,iz2graze,&
          & iout_cosine,nspool_cosine,ico2s,ispm,spm0,ised 
  namelist /CORE/ gmaxs1,gmaxs2,pis1,pis2,kno3s1,knh4s1,kpo4s1,kco2s1,kno3s2,&
          & knh4s2,kpo4s2,kco2s2,ksio4s2,kns1,kns2,alpha1,alpha2,beta,ak1,ak2,&
          & ak3,gammas1,gammas2,beta1,beta2,kgz1,kgz2,rho1,rho2,rho3,gamma1,&
          & gamma2,gammaz,kex1,kex2,wss2,wsdn,wsdsi,si2n,p2n,o2no,o2nh,c2n,&
          & kox,ipo4,kmdn1,kmdn2,kmdsi1,kmdsi2,gamman,TR,pco2a
  namelist /MISC/ iws,NO3c,ws1,ws2,iclam,deltaZ,kcex,Nperclam,Wclam,Fclam,&
          & nclam0,fS2,fDN,fDSi,rkS2,rkDN,rkDSi,mkS2,mkDN,mkDSi

  !initialize parameter values
  idelay=0; ndelay=7; ibgraze=0; idapt=0; alpha_corr=1.25; zeptic=10.0; iz2graze=1
  iout_cosine=0; nspool_cosine=60; ico2s=0; ispm=0; spm0=20.0; ised=1; 
  gmaxs1=3.0; gmaxs2=2.5; pis1=1.5; pis2=1.5; kno3s1=1.0; 
  knh4s1=0.15; kpo4s1=0.1; kco2s1=50.0; kno3s2=3.0; knh4s2=0.45; kpo4s2=0.1;
  kco2s2=50.0; ksio4s2=4.5; kns1=0.0; kns2=0.0; alpha1=0.1; alpha2=0.1; beta=0.0; 
  ak1=0.75; ak2=0.03; ak3=0.066; gammas1=0.5; gammas2=0.3; beta1=0.75; beta2=0.5; 
  kgz1=0.5; kgz2=0.25; rho1=0.6; rho2=0.3; rho3=0.1; gamma1=0.75; gamma2=0.75; 
  gammaz=0.05; kex1=0.2; kex2=0.3; wss2=0.25; wsdn=0.5; wsdsi=0.5; si2n=1.2; p2n=0.0625; 
  o2no=8.625; o2nh=6.625; c2n=7.3; kox=30.0; ipo4=0; kmdn1=0.009; kmdn2=0.075; kmdsi1=0.0114; 
  kmdsi2=0.015; gamman=0.07; TR=20.0; pco2a=400.0; iws=0; NO3c=2.0; ws1=2.5; ws2=2.0
  iclam=0; deltaZ=1.0; kcex=0.002; Nperclam=0.39032; Wclam=5.45e-3; Fclam=40.0; 
  nclam0=2000; fS2=0.0; fDN=0.0; fDSi=0.0; rkS2=4e-3; rkDN=4e-3; rkDSi=4e-3;
  mkS2=0.1; mkDN=0.1; mkDSi=0.1

  !read parameter values
  open(31,file=in_dir(1:len_in_dir)//'cosine.nml',delim='apostrophe',status='old')
  read(31,nml=MARCO); read(31,nml=CORE); read(31,nml=MISC)
  close(31)
  
  !allocate sediment variables
  if(ised==1) then
    allocate(PS2(3,nea),PDN(3,nea),PDSi(3,nea),RS2(3,nea),RDN(3,nea),RDSi(3,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. fS2')
  endif

  if(myrank==0) then
    open(31,file=out_dir(1:len_out_dir)//'cosine.out.nml',status='replace')
    write(31,nml=MARCO); write(31,nml=CORE); write(31,nml=MISC)
    close(31)
  endif

  !read in station info. 
  if(iout_cosine==1) call read_cosine_stainfo

  !read bottom grazing information
  if(ibgraze==1) then
    bgraze=0.0
    call read_gr3_prop('bgraze',-9999.d0,bgraze,nea)
  elseif(ibgraze==2) then !temporally and spatially varying inputs
    open(455,file=in_dir(1:len_in_dir)//'bgraze.th',status='old') 
    bgraze=0.0
    time_cosine(2)=-999.0
  elseif(ibgraze/=0) then
    call parallel_abort('unknown ibgraze')
  endif

  !read dynamic clam number
  nclam=0;
  if(iclam==1) then
    nclam=nclam0  
  elseif (iclam==2) then
    call read_gr3_prop('nclam',-999.d0,nclam,nea)
  elseif (iclam==3) then
    open(456,file=in_dir(1:len_in_dir)//'nclam.th',status='old') 
    time_cosine(3)=-999.0
  endif 

  !read SPM information
  SPM=0.0
  if(ispm==0) then !constant
    SPM=spm0
  elseif(ispm==1) then !spatial varying
    call read_gr3_prop('SPM',-999.d0,SPM(1,:),nea)
    do k=2,nvrt; SPM(k,:)=SPM(1,:); enddo
  elseif(ispm==2) then !call SED3D model
#ifndef USE_SED
    call parallel_abort('ispm=2, need to turn on SED module')
#endif 
  elseif(ispm==3) then !use spatial and temporal varying SPM
    open(453,file=in_dir(1:len_in_dir)//'SPM.th',status='old')
    time_cosine(1)=-999.0
  else
    call parallel_abort
  endif

  !read sediment flux model parameters
  if(ised==1) then
    !initialize sediment variables
    !todo, include these variables in hotstart.nc; update to sediment flux model
    PS2=0.0; PDN=0.0; PDSi=0.0; RS2=0.0; RDN=0.0; RDSi=0.0
    if(ihot==0) then !temporary fix, need to update
      do i=1,3
        write(snum,*)i
        call read_gr3_prop('PS2_'//trim(adjustl(snum)),-999.d0,PS2(i,:),nea)
        call read_gr3_prop('RS2_'//trim(adjustl(snum)),-999.d0,RS2(i,:),nea)
      enddo

      do i=1,3
        write(snum,*)i
        call read_gr3_prop('PDN_'//trim(adjustl(snum)),-999.d0,PDN(i,:),nea)
        call read_gr3_prop('RDN_'//trim(adjustl(snum)),-999.d0,RDN(i,:),nea)
      enddo

      do i=1,3
        write(snum,*)i
        call read_gr3_prop('PDSi_'//trim(adjustl(snum)),-999.d0,PDSi(i,:),nea)
        call read_gr3_prop('RDSi_'//trim(adjustl(snum)),-999.d0,RDSi(i,:),nea)
      enddo
    endif! ihot=0
  endif !ised
 
  return
end subroutine read_cosine_param

subroutine read_cosine_stainfo
!---------------------------------------------------------------------------
!CoSiNE outputs
!---------------------------------------------------------------------------
  use schism_glbl, only : rkind,dt,ihot,ne,i34,xnd,ynd,elnode,ielg, &
      & in_dir,out_dir,len_in_dir,len_out_dir,ics,pi,rearth_eq,ne_global,np_global
  use schism_msgp, only : myrank,nproc,comm,parallel_abort,parallel_barrier, &
      & itype,rtype
  use cosine_misc, only : pt_in_poly
  use cosine_mod, only : nvar,dvar,nsta_lc,nsta,sie,sid,sdep,nstas,displ,sids,dvars
  implicit none
  include 'mpif.h'

  !local variables
  integer :: i,j,m,id,irank,istat,inside,nodel(3),ierr,negb,npgb
  real(rkind) :: rtmp,xtmp,ytmp,x(4),y(4),arco(3)
  integer,allocatable :: ista(:),iep(:),i34gb(:),elnodegb(:,:),sid_gb(:)
  real(rkind), allocatable :: sx(:),sy(:),sz(:),xgb(:),ygb(:)

  !read grid and station information on myrank=0
  if(myrank==0) then
    !read hgird information
    open(31,file=in_dir(1:len_in_dir)//'hgrid.gr3',status='old')
    read(31,*); read(31,*)negb,npgb
    if(negb/=ne_global.or.npgb/=np_global) call parallel_abort('Check: negb and npgb in hgrid.gr3') 
    allocate(i34gb(negb),elnodegb(4,negb),xgb(npgb),ygb(npgb),stat=istat)
    if(istat/=0) call parallel_abort('failed to alloc. i34gb')
    do i=1,npgb; read(31,*)j,xgb(i),ygb(i),rtmp; enddo
    do i=1,negb; read(31,*)j,i34gb(i),elnodegb(1:i34gb(i),i); enddo
    close(31)

    !read station info. 
    open(31,file=in_dir(1:len_in_dir)//'cstation.in',status='old')
    read(31,*); read(31,*)nsta
    allocate(sx(nsta),sy(nsta),sz(nsta),ista(nsta),iep(nsta),sids(nsta), &
           & nstas(nproc),displ(nproc),dvars(nvar),stat=istat)
    if(istat/=0) call parallel_abort('failed to alloc. sx')
    do i=1,nsta; read(31,*)j,sx(i),sy(i),sz(i); enddo
    close(31)
  
    !find parent element
    ista=0; iep=-1
    do m=1,nsta
      do i=1,ne_global
        x(1:i34gb(i))=xgb(elnodegb(1:i34gb(i),i))
        y(1:i34gb(i))=ygb(elnodegb(1:i34gb(i),i))
        call pt_in_poly(i34gb(i),x(1:i34gb(i)),y(1:i34gb(i)),sx(m),sy(m),inside,arco,nodel)
        if(inside==1) then
          ista(m)=1 
          iep(m)=i
        endif !if
        if(ista(m)==1) exit
      enddo !i
    enddo !m

    deallocate(ista,i34gb,elnodegb,xgb,ygb,sx,sy)
  endif !if(myrank==0)

  !boradcast station information
  call mpi_bcast(nsta,1,itype,0,comm,ierr) 
  if(.not.allocated(iep)) then
     allocate(iep(nsta),sz(nsta),stat=istat)
     if(istat/=0) call parallel_abort('failed to alloc. iep')
  endif  
  call mpi_bcast(iep,nsta,itype,0,comm,ierr) 
  call mpi_bcast(sz,nsta,rtype,0,comm,ierr) 

  !compute local nsta,sid,sie
  do m=1,2
    !initilize
    if(m==2) then
      allocate(sdep(nsta_lc),sie(nsta_lc),sid(nsta_lc),sid_gb(nsta_lc),dvar(nvar),stat=istat) 
      if(istat/=0) call parallel_abort('failed to alloc. sdep')
    endif

    !check each station pt in local domain
    nsta_lc=0
    do j=1,nsta
      do i=1,ne
        if(ielg(i)==iep(j)) then
          nsta_lc=nsta_lc+1
          if(m==1) exit
          sdep(nsta_lc)=sz(j)
          sie(nsta_lc)=i
          sid_gb(nsta_lc)=j
        endif 
      enddo !i
    enddo !j

    !find the local index for each station in local domain
    if(m==2) then
      sid=0; id=0
      do j=1,ne
        do i=1,nsta_lc
          if(sie(i)==j) then
            id=id+1
            sid(id)=i
          endif
        enddo
      enddo
    endif !if
  enddo!m

  !write(361+myrank,*) sie
  !call parallel_barrier
  !stop

  !collect nsta_lc
  call mpi_gather(nsta_lc,1,itype,nstas,1,itype,0,comm,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('failed in gather nsta_lc')
  if(myrank==0) then
    displ=0
    do i=2,nproc
      displ(i)=displ(i-1)+nstas(i-1)
    enddo
  endif
  call mpi_gatherv(sid_gb,nsta_lc,itype,sids,nstas,displ,itype,0,comm,ierr)

  deallocate(iep,sid_gb,sz)
end subroutine read_cosine_stainfo

subroutine cosine_output(vid,id,varname,ndata_lc,rarray,imode)
  use schism_glbl, only : rkind,errmsg,out_dir,len_out_dir
  use schism_msgp, only : myrank,nproc,comm,parallel_abort,parallel_barrier, &
      & itype,rtype
  use cosine_mod, only: nvar,dcosine,dvar,dvars,nsta_lc,nsta,nstas,ndata_gb,istat_cosine,& 
      & sids,displ
  use netcdf
  implicit none
  include 'mpif.h'
  
  integer :: stype=MPI_CHARACTER
  integer,intent(in) :: imode,vid,id
  character(*),intent(in) :: varname
  integer,intent(in) :: ndata_lc
  real(rkind),dimension(*),intent(in) :: rarray
 
  !local variables 
  integer :: i,j,m,n,istat,ierr
  integer :: iret,ncid,ndim,time_dim,nsta_dim
  integer,allocatable :: ndata(:),ndatas(:),dims(:),var_dims(:)
  real(rkind),allocatable :: swild(:),swild_gb(:) 
  character(len=30) :: stmp
  character(len=30), allocatable :: varname_gb(:)

  if(imode==0) then
    !allocate variables 
    if(dvar(vid)%istat==0) then
       dvar(vid)%istat=1

       dvar(vid)%ndata=ndata_lc
       dvar(vid)%varname=varname 
       allocate(dvar(vid)%data(ndata_lc,nsta_lc),stat=istat)
       if(istat/=0) call parallel_abort('failed in alloc. dvar(vid)%data')
    endif

    !save variable values 
    dvar(vid)%data(:,id)=rarray(1:ndata_lc)
  elseif(imode==1) then
    !------------------------------------------------------------------------------ 
    !compute total datasize for each station, and initilize dvars
    !------------------------------------------------------------------------------ 
    if(istat_cosine==0) then
      istat_cosine=1

      !compute dimensions
      allocate(ndata(nvar),ndatas(nvar),varname_gb(nproc),dims(nvar),var_dims(nvar), stat=istat)
      if(istat/=0) call parallel_abort('failed in alloc. ndata')
      ndata=0; do i=1,nvar; ndata(i)=dvar(i)%ndata; enddo
      call mpi_allreduce(ndata,ndatas,nvar,itype,MPI_MAX,comm,ierr)
      ndata_gb=sum(ndatas)

      !initialize dvars
      ndim=0; dims=0
      do i=1,nvar
        call mpi_gather(dvar(i)%varname,30,stype,varname_gb,30,stype,0,comm,ierr)

        if(myrank==0) then
          !determine varname
          do j=1,nproc
            stmp=varname_gb(j)
            if(len(trim(adjustl(stmp)))/=0) exit 
          enddo 

          !initilize dvars variables
          if(ndatas(i)>0) then
            dvars(i)%istat=1; dvars(i)%ndata=ndatas(i)
            dvars(i)%varname=stmp
            allocate(dvars(i)%data(ndatas(i),nsta),stat=istat)
            if(istat/=0) call parallel_abort('failed in alloc. dvars')
          endif

          !find variables with ndata>1
          if(ndatas(i)>1) then
            istat=0
            do j=1,ndim
              if(dims(j)==ndatas(i)) istat=1
            enddo
            
            !new dimension found
            if(istat==0) ndim=ndim+1; dims(ndim)=ndatas(i)
          endif
        endif !myrank
      enddo 

      !create station output
      if(myrank==0) then
        iret=nf90_create(trim(adjustl(out_dir(1:len_out_dir)//'cosine.nc')),OR(NF90_NETCDF4,NF90_CLOBBER),ncid)
        dcosine%ncid=ncid

        !define dimension
        iret=nf90_def_dim(ncid,'time',NF90_UNLIMITED,time_dim)
        iret=nf90_def_dim(ncid,'nstation',nsta,nsta_dim)
        do i=1,ndim
          write(stmp,*)dims(i)
          iret=nf90_def_dim(ncid,trim(adjustl(stmp)),dims(i),var_dims(i))
        enddo

        !define variables
        iret=nf90_def_var(ncid,'time',nf90_double,(/time_dim/),dcosine%varid)
        do i=1,nvar
          if(dvars(i)%ndata==1) then
            iret=nf90_def_var(ncid,trim(adjustl(dvars(i)%varname)),nf90_FLOAT,(/time_dim,nsta_dim/),dvars(i)%varid)
          elseif(dvars(i)%ndata>1) then
            do j=1,ndim; if(dims(j)==dvars(i)%ndata) exit; enddo
            iret=nf90_def_var(ncid,trim(adjustl(dvars(i)%varname)),nf90_FLOAT,(/time_dim,var_dims(j),nsta_dim/),dvars(i)%varid)
          elseif(dvars(i)%ndata/=0) then
            call parallel_abort('error in def cosine diagnostic variables')
          endif
        enddo

        iret=nf90_enddef(ncid)
      endif !myrank

      deallocate(ndata,ndatas,varname_gb) 
    endif !istat_cosine

    !------------------------------------------------------------------------------ 
    !allocate data for storing all diagnostic values 
    !------------------------------------------------------------------------------ 
    allocate(swild(ndata_gb*nsta_lc),stat=istat)
    if(istat/=0) call parallel_abort('failed in alloc. swild')
    if(myrank==0) then
      allocate(swild_gb(ndata_gb*nsta),stat=istat)
      if(istat/=0) call parallel_abort('failed in alloc. swild_gb')
    endif

    !assemble data 
    swild=0.0; m=0 
    do i=1,nsta_lc
      do j=1,nvar
        do n=1,dvar(j)%ndata
          m=m+1
          swild(m)=dvar(j)%data(n,i)
        enddo !n
      enddo !j=1,nvar
    enddo !i=1,nsta_lc

    !pass all data to myrank=0
    call mpi_gatherv(swild,nsta_lc*ndata_gb,rtype, swild_gb, nstas*ndata_gb,displ*ndata_gb,rtype,0,comm,ierr)
    if(myrank==0) then
      m=0
      do i=1,nsta
        do j=1,nvar
          do n=1,dvars(j)%ndata
            m=m+1; dvars(j)%data(n,sids(i))=swild_gb(m) 
          enddo
        enddo
      enddo

      !write station output
      iret=nf90_put_var(dcosine%ncid,dcosine%varid,(/dcosine%time/),start=(/dcosine%it/),count=(/1/))
      do i=1,nvar
        if(dvars(i)%ndata==1) then
          iret=nf90_put_var(dcosine%ncid,dvars(i)%varid,dvars(i)%data(:,1),start=(/dcosine%it,1/),count=(/1,nsta/))
        elseif(dvars(i)%ndata>1) then
          iret=nf90_put_var(dcosine%ncid,dvars(i)%varid,dvars(i)%data(:,:),start=(/dcosine%it,1,1/),count=(/1,dvars(i)%ndata,nsta/))
        endif
      enddo
      iret=nf90_sync(dcosine%ncid) 

    endif

    deallocate(swild)
    if(myrank==0) deallocate(swild_gb)
  endif !imode
    
end subroutine cosine_output

!subroutine read_cosine_stainfo
!!---------------------------------------------------------------------------
!!output cosine parameters to check
!!---------------------------------------------------------------------------
!  use cosine_mod, only : nsta,ista,depsta,stanum,nspool_cosine
!  use schism_glbl, only : rkind,dt,ihot,ne,i34,xnd,ynd,elnode, &
!      & in_dir,out_dir,len_in_dir,len_out_dir,ics,pi,rearth_eq
!  use schism_msgp, only : myrank,nproc,parallel_abort
!  use cosine_misc, only : pt_in_poly
!  implicit none
!
!  !local variables
!  integer,parameter :: maxsta=10000,maxl=100 !maximum station
!  integer :: i,j,istat,nstation,nodel(3),inside,id,iflag,mid,msta,nstai(ne),stanumi(maxl,ne)
!  real(rkind) :: slx(maxsta),sly(maxsta),sdep(maxsta),x(3),y(3),arco(3),depstai(maxl,ne)
!  real(rkind) :: xtmp,ytmp
!  character(len=4) :: fn
!  logical :: lexist
!
!  !read station info.
!  open(450,file=in_dir(1:len_in_dir)//'cstation.in',status='old')
!  read(450,*)
!  read(450,*)nstation
!  do i=1,nstation
!    read(450,*)j,slx(i),sly(i),sdep(i) 
!    if(ics==2) then
!      xtmp=slx(i)*pi/180.0; ytmp=sly(i)*pi/180.0
!      slx(i)=rearth_eq*cos(ytmp)*cos(xtmp)
!      sly(i)=rearth_eq*cos(ytmp)*sin(xtmp)
!    endif
!  enddo
!  close(450)
!
!  !alloc.
!  allocate(ista(ne),stat=istat) 
!  if(istat/=0) call parallel_abort('failure in alloc. ista')
!
!  !determine the elements with values to be checked
!  id=0; ista=0; nstai=0;
!  msta=-100; depstai=-9999
!  do i=1,ne
!    iflag=0
!    do j=1,nstation
!      x=xnd(elnode(1:i34(i),i))
!      y=ynd(elnode(1:i34(i),i))
!      call pt_in_poly(i34(i),x,y,slx(j),sly(j),inside,arco,nodel)
!      if(inside==1) then
!        if(ista(i)==0) then
!          id=id+1
!          ista(i)=id
!        endif
!        nstai(id)=nstai(id)+1
!        depstai(nstai(id),id)=sdep(j)
!        stanumi(nstai(id),id)=j
!        msta=max(msta,nstai(id))
!      endif 
!    enddo !j
!  enddo !i
!  mid=id !number of elements
!
!  if(mid==0) return 
!
!  !alloc.
!  allocate(nsta(mid),depsta(msta,mid),stanum(msta,mid),stat=istat) 
!  if(istat/=0) call parallel_abort('failure in alloc. nsta')
!
!  nsta=0; depsta=-9999
!  do i=1,mid
!    nsta(i)=nstai(i)
!    do j=1,nsta(i)
!      depsta(j,i)=depstai(j,i)
!      stanum(j,i)=stanumi(j,i)
!    enddo
!  enddo
!
!  !open a station output file
!  write(fn,'(i4.4)')myrank
!  inquire(file=out_dir(1:len_out_dir)//'cstation_'//fn//'.out',exist=lexist)
!  if(ihot<=1.or.(ihot==2.and.(.not.lexist))) then
!    open(451,file=out_dir(1:len_out_dir)//'cstation_'//fn//'.out',form='unformatted',status='replace')
!    write(451)sum(nsta),dt*nspool_cosine
!  elseif(ihot==2.and.lexist) then
!    open(451,file=out_dir(1:len_out_dir)//'cstation_'//fn//'.out',form='unformatted',access='append',status='old')
!  else
!    call parallel_abort('unknown ihot, CoSiNE')
!  endif
!
!  return
!end subroutine read_cosine_stainfo
