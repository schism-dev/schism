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
  logical :: lexist

  !define namelist
  namelist /MARCO/ idelay,ndelay,ibgraze,idapt,alpha_corr,zeptic,iz2graze,&
          & iout_cosine,nspool_cosine,ico2s,ispm,spm0,ised 
  namelist /CORE/ gmaxs,pis,kno3s,knh4s,kpo4s,kco2s,&
          & ksio4,kns,alphas,betas,aks,gammas,betaz,kgz,rhoz,alphaz,&
          & gammaz,kez,wss2,wsdn,wsdsi,si2n,p2n,o2no,o2nh,c2n,&
          & kox,ipo4,kmdn,kmdsi,gamman,TR,pco2a
  namelist /MISC/ iws,NO3c,ws1,ws2,iclam,deltaZ,kcex,Nperclam,Wclam,Fclam,&
          & nclam0,fS2,fDN,fDSi,rkS2,rkDN,rkDSi,mkS2,mkDN,mkDSi

  !initialize parameter values
  idelay=0; ndelay=7; ibgraze=0; idapt=0; alpha_corr=1.25; zeptic=10.0; iz2graze=1
  iout_cosine=0; nspool_cosine=60; ico2s=0; ispm=0; spm0=20.0; ised=1; 
  gmaxs=0; pis=0; kno3s=0; knh4s=0; kpo4s=0; kco2s=0; ksio4=0; kns=0; alphas=0; betas=0; 
  aks=0; gammas=0;  betaz=0; kgz=0; rhoz=0;  alphaz=0;
  gammaz=0; kez=0; wss2=0.25; wsdn=0.5; wsdsi=0.5; si2n=1.2; p2n=0.0625; 
  o2no=8.625; o2nh=6.625; c2n=7.3; kox=30.0; ipo4=0; kmdn=0; kmdsi=0; 
  gamman=0.07; TR=20.0; pco2a=400.0; iws=0; NO3c=2.0; ws1=2.5; ws2=2.0
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

  !variable names for outputs
  allocate(name_cos(ntrc),stat=istat)
  if(istat/=0) call parallel_abort('failed in alloc. name_cos')
  name_cos(1:13)=(/'NO3 ','SiO4','NH4 ','S1  ','S2  ','Z1  ','Z2  ', & 
                 & 'DN  ','DSi ','PO4 ','DOX ','CO2 ','ALK '/)

  !read in station info. 
  if(iout_cosine/=0) call read_cosine_stainfo

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
        inquire(file=in_dir(1:len_in_dir)//'PS2_'//trim(adjustl(snum))//'.gr3',exist=lexist)
        if(lexist) call read_gr3_prop('PS2_'//trim(adjustl(snum)),-999.d0,PS2(i,:),nea)

        inquire(file=in_dir(1:len_in_dir)//'RS2_'//trim(adjustl(snum))//'.gr3',exist=lexist)
        if(lexist) call read_gr3_prop('RS2_'//trim(adjustl(snum)),-999.d0,RS2(i,:),nea)
      enddo

      do i=1,3
        write(snum,*)i
        inquire(file=in_dir(1:len_in_dir)//'PDN_'//trim(adjustl(snum))//'.gr3',exist=lexist)
        if(lexist) call read_gr3_prop('PDN_'//trim(adjustl(snum)),-999.d0,PDN(i,:),nea)

        inquire(file=in_dir(1:len_in_dir)//'RDN_'//trim(adjustl(snum))//'.gr3',exist=lexist)
        if(lexist) call read_gr3_prop('RDN_'//trim(adjustl(snum)),-999.d0,RDN(i,:),nea)
      enddo

      do i=1,3
        write(snum,*)i
        inquire(file=in_dir(1:len_in_dir)//'PDSi_'//trim(adjustl(snum))//'.gr3',exist=lexist)
        if(lexist) call read_gr3_prop('PDSi_'//trim(adjustl(snum)),-999.d0,PDSi(i,:),nea)

        inquire(file=in_dir(1:len_in_dir)//'RDSi_'//trim(adjustl(snum))//'.gr3',exist=lexist)
        if(lexist) call read_gr3_prop('RDSi_'//trim(adjustl(snum)),-999.d0,RDSi(i,:),nea)
      enddo
    endif! ihot=0
  endif !ised

  !read spatially varying parameters
  call cosine_vars_init
 
end subroutine read_cosine_param

subroutine cosine_vars_init
  !--------------------------------------------------------------------------------
  !allocate cosine arrays and initialize
  !--------------------------------------------------------------------------------
  use schism_glbl, only : rkind,nea,npa,nvrt,ntrs,in_dir,len_in_dir,np_global, &
                        & ne_global,ielg,iplg,i34,elnode
  use schism_msgp, only : parallel_abort,myrank,comm,itype,rtype
  use netcdf
  use cosine_mod
  use misc_modules
  implicit none

  !local variables
  integer :: istat,i,j,k,ie,m,n,ip,ncid,varid,npt,nsp,ndim,dimid(3),dims(3)
  real(rkind) :: data0,swild(max(np_global,ne_global))
  character(len=15),allocatable :: pname(:)
  character(len=20) :: fname
  type(cosine_spatial_param),pointer :: p

  !---------------------------------------------------------------------------
  !to add a spatially varying parameter
  !  1). append parameter name in pname array
  !  2). make links by piointing p/p1/p2 for scalar/1D/2D variable
  !---------------------------------------------------------------------------
  !define spatial varying parameters 
  fname='COS_param.nc';  nsp=40
  allocate(pname(nsp),sp(nsp),stat=istat)
  if(istat/=0) call parallel_abort('Failed in alloc. pname')

  m=0
  !global and core modules
  pname(1:32)=(/'gmaxs ','gammas','pis   ','kno3s ','knh4s ', &
              & 'kpo4s ','kco2s ','ksio4 ','kns   ','alphas', &
              & 'betas ','aks   ','betaz ','alphaz','gammaz', &
              & 'kez   ','kgz   ','rhoz  ','TR    ','kox   ', &
              & 'wss2  ','wsdn  ','wsdsi ','si2n  ','p2n   ', &
              & 'o2no  ','o2nh  ','c2n   ','gamman','pco2a ', &
              & 'kmdn  ','kmdsi '/)

  sp(m+1)%p1=>gmaxs; sp(m+2)%p1=>gammas; sp(m+3)%p1=>pis;    sp(m+4)%p1=>kno3s; sp(m+5)%p1=>knh4s;   m=m+5
  sp(m+1)%p1=>kpo4s; sp(m+2)%p1=>kco2s;  sp(m+3)%p=>ksio4;   sp(m+4)%p1=>kns;   sp(m+5)%p1=>alphas;  m=m+5
  sp(m+1)%p1=>betas; sp(m+2)%p1=>aks;    sp(m+3)%p1=>betaz;  sp(m+4)%p1=>alphaz;sp(m+5)%p1=>gammaz;  m=m+5
  sp(m+1)%p1=>kez;   sp(m+2)%p1=>kgz;    sp(m+3)%p1=>rhoz;   sp(m+4)%p=>TR;     sp(m+5)%p=>kox;      m=m+5
  sp(m+1)%p=>wss2;   sp(m+2)%p=>wsdn;    sp(m+3)%p=>wsdsi;   sp(m+4)%p=>si2n;   sp(m+5)%p=>p2n;      m=m+5
  sp(m+1)%p=>o2no;   sp(m+2)%p=>o2nh;    sp(m+3)%p=>c2n;     sp(m+4)%p=>gamman; sp(m+5)%p=>pco2a;    m=m+5
  sp(m+1)%p1=>kmdn;  sp(m+2)%p1=>kmdsi;  m=m+2

  !read spatially varying parameters
  do m=1,nsp
    p=>sp(m)
    !get dimension info. about parameter
    p%dims=(/1,1/)
    if(associated(p%p)) then
      p%ndim=1; p%data0(1)=p%p
    elseif(associated(p%p1)) then
      p%ndim=2; p%data0(1:size(p%p1))=p%p1; p%dims(1)=size(p%p1)
    elseif(associated(p%p2)) then
      p%ndim=3; p%data0(1:size(p%p2))=reshape(p%p2,(/n/)); p%dims=shape(p%p2)
    else
      cycle
    endif
    p%varname=trim(adjustl(pname(m)))
    allocate(p%istat(p%dims(1),p%dims(2))); p%istat=0

    !read parameter data
    ip=0
    do i=1,p%dims(2)
      do k=1,p%dims(1)
        ip=ip+1; data0=p%data0(ip); npt=0
        if(abs(data0+999.d0)>1.d-6.and.abs(data0+9999.d0)>1.d-6) cycle
        if(.not.allocated(p%data)) allocate(p%data(nea,p%dims(1),p%dims(2)))
        p%istat(k,i)=1

        !read value on myrank=0, then bcast
        if(myrank==0) then
          j=nf90_open(in_dir(1:len_in_dir)//trim(adjustl(fname)),OR(NF90_NETCDF4,NF90_NOWRITE),ncid)
          if(j/=NF90_NOERR) call parallel_abort(trim(adjustl(fname))//': open')
          j=nf90_inq_varid(ncid,trim(adjustl(p%varname)),varid)
          if(j/=NF90_NOERR) call parallel_abort(trim(adjustl(p%varname))//': wrong varid' )
          j=nf90_inquire_variable(ncid,varid,ndims=ndim)
          if(j/=NF90_NOERR) call parallel_abort(trim(adjustl(p%varname))//': wrong ndim')
          j=nf90_inquire_variable(ncid,varid,dimids=dimid(1:ndim))
          if(j/=NF90_NOERR) call parallel_abort(trim(adjustl(p%varname))//': wrong dimid')
          j=nf90_inquire_dimension(ncid,dimid(1),len=npt)
          if(j/=NF90_NOERR) call parallel_abort(trim(adjustl(p%varname))//': wrong npt')
          if(npt/=np_global.and.npt/=ne_global) call parallel_abort(trim(adjustl(p%varname))//': npt/=ne,np' )
          if(p%ndim==1) j=nf90_get_var(ncid,varid,swild(1:npt), (/1/),(/npt/))
          if(p%ndim==2) j=nf90_get_var(ncid,varid,swild(1:npt), (/1,k/),(/npt,1/))
          if(p%ndim==3) j=nf90_get_var(ncid,varid,swild(1:npt), (/1,k,i/),(/npt,1,1/))
          if(j/=NF90_NOERR) call parallel_abort(trim(adjustl(p%varname))//': wrong in read value')
          j=nf90_close(ncid)
        endif
        call mpi_bcast(npt,1,itype,0,comm,istat)
        call mpi_bcast(swild(1:npt),npt,rtype,0,comm,istat)

        !get parameter value for each rank
        do ie=1,nea
          p%data(ie,k,i)=0.d0
          if(npt==ne_global) then
            p%data(ie,k,i)=swild(ielg(ie))
          elseif(npt==np_global) then
            do n=1,i34(ie); p%data(ie,k,i)=p%data(ie,k,i)+swild(iplg(elnode(n,ie)))/dble(i34(ie)); enddo
          else
            call parallel_abort(trim(adjustl(p%varname))//': wrong npt')
          endif
        enddo !ie

      enddo !k
    enddo !i
  enddo !m
end subroutine cosine_vars_init


subroutine update_cosine_vars(id)
!--------------------------------------------------------------------
!get 2D parameter value of element
!--------------------------------------------------------------------
  use schism_glbl, only : rkind,irange_tr,tr_el,nvrt,i34,elside,isdel, &
                        & elnode
  use cosine_mod
  implicit none
  integer, intent(in) :: id

  !local variable
  integer :: i,j,k,m
  type(cosine_spatial_param),pointer :: p

  !spatial varying parameters
  do m=1,size(sp)
    p=>sp(m)
    if(p%ndim==0) cycle
    do i=1,p%dims(2)
      do k=1,p%dims(1)
        if(p%istat(k,i)==0) cycle
        if(p%ndim==1) p%p=p%data(id,1,1)
        if(p%ndim==2) p%p1(k)=p%data(id,k,1)
        if(p%ndim==3) p%p2(k,i)=p%data(id,k,i)
      enddo !k
    enddo !i
  enddo !m

end subroutine update_cosine_vars


subroutine read_cosine_stainfo
!---------------------------------------------------------------------------
!read CoSiNE station 
!---------------------------------------------------------------------------
  use schism_glbl, only : rkind,dt,ihot,ne,i34,xnd,ynd,elnode,ielg, &
      & in_dir,out_dir,len_in_dir,len_out_dir,ne_global,np_global
  use schism_msgp, only : myrank,nproc,comm,parallel_abort,parallel_barrier, &
      & itype,rtype
  use cosine_misc, only : pt_in_poly
  use cosine_mod, only : dg,dl,mval
  implicit none
  include 'mpif.h'

  !local variables
  integer :: i,j,m,id,irank,istat,inside,nodel(3),ierr,negb,npgb,nsta
  real(rkind) :: rtmp,xtmp,ytmp,x(4),y(4),arco(3)
  integer,allocatable :: i34gb(:),elnodegb(:,:),ista(:),iep(:)
  real(rkind), allocatable :: xgb(:),ygb(:),z(:)

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
    allocate(dg%x(nsta),dg%y(nsta),dg%z(nsta),dg%nstas(nproc),dg%sids(nsta), &
           & dg%displ(nproc),dg%iep(nsta),ista(nsta),iep(nsta),z(nsta),stat=istat)
    if(istat/=0) call parallel_abort('failed to alloc. dcosine%x')
    do i=1,nsta; read(31,*)j,dg%x(i),dg%y(i),dg%z(i); enddo
    close(31)
  
    !find parent element
    ista=0; dg%iep=-1
    do m=1,nsta
      do i=1,ne_global
        x(1:i34gb(i))=xgb(elnodegb(1:i34gb(i),i))
        y(1:i34gb(i))=ygb(elnodegb(1:i34gb(i),i))
        call pt_in_poly(i34gb(i),x(1:i34gb(i)),y(1:i34gb(i)),dg%x(m),dg%y(m),inside,arco,nodel)
        if(inside==1) then
          ista(m)=1 
          dg%iep(m)=i
        endif !if
        if(ista(m)==1) exit
      enddo !i
    enddo !m

    deallocate(i34gb,elnodegb,xgb,ygb)
  endif !if(myrank==0)

  !boradcast global station information
  call mpi_bcast(nsta,1,itype,0,comm,ierr) 
  if(.not.allocated(dg%iep)) then
     allocate(dg%iep(nsta),dg%z(nsta),stat=istat)
     if(istat/=0) call parallel_abort('failed to alloc. iep')
  endif  
  call mpi_bcast(dg%iep,nsta,itype,0,comm,ierr) 
  call mpi_bcast(dg%z,nsta,rtype,0,comm,ierr) 
  dg%nsta=nsta
  call parallel_barrier 

  !compute local nsta,sid,iep
  do m=1,2
    !initilize
    if(m==2) then
      allocate(dl%z(nsta),dl%iep(nsta),stat=istat) 
      if(istat/=0) call parallel_abort('failed to alloc. sdep')
    endif

    !check each station pt in local domain
    nsta=0
    do j=1,dg%nsta
      do i=1,ne
        if(ielg(i)==dg%iep(j)) then
          nsta=nsta+1
          if(m==1) exit
          dl%z(nsta)=dg%z(j)
          dl%iep(nsta)=dg%iep(j)
        endif 
      enddo !i
    enddo !j
  enddo!m
  dl%nsta=nsta

  !collect local nsta and iep,z
  call mpi_gather(dl%nsta,1,itype,dg%nstas,1,itype,0,comm,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('failed in gather local nsta')
  if(myrank==0) then
    dg%displ=0
    do i=2,nproc
      dg%displ(i)=dg%displ(i-1)+dg%nstas(i-1)
    enddo
  endif
  call mpi_gatherv(dl%iep,dl%nsta,itype,iep,dg%nstas,dg%displ,itype,0,comm,ierr)
  call mpi_gatherv(dl%z,dl%nsta,rtype,z,dg%nstas,dg%displ,rtype,0,comm,ierr)

  !compute the indices of local iep
  if(myrank==0) then
    ista=0
    do i=1,dg%nsta
      do j=1,dg%nsta
        if(ista(j)==1) cycle 
        if(iep(i)==dg%iep(j) .and. abs(z(i)-dg%z(j))<mval) then
          dg%sids(i)=j; ista(j)=1; exit
        endif
      enddo
    enddo

    deallocate(ista,z,iep)
  endif

end subroutine read_cosine_stainfo

subroutine cosine_output(imode,id,varname,ndim,rarray)
!---------------------------------------------------------------------------
!CoSiNE station outputs
!imode=0: for local rank to initialize and store local diagnostic variables
!         id: the station index 
!         varname: variable name
!         ndim: number of variable dimensions
!         rarray: variable values
!imode=1: for myrank=0 to initialize and store all diagnostic variables
!---------------------------------------------------------------------------
  use schism_glbl, only : rkind,errmsg,out_dir,len_out_dir,ihot
  use schism_msgp, only : myrank,nproc,comm,parallel_abort,parallel_barrier, &
      & itype,rtype
  use cosine_mod, only: dl,dg,dlv,dgv,cosine_diagnostic_variable
  use netcdf
  implicit none
  include 'mpif.h'
  
  integer :: stype=MPI_CHARACTER
  integer,intent(in) :: imode,id
  character(*),intent(in) :: varname
  integer,intent(in) :: ndim
  real(rkind),dimension(*),intent(in) :: rarray
 
  !local variables 
  integer :: i,j,m,n,istat,iflag,ierr,id_x,id_y,id_z,id_iep
  integer :: ndim_lc,ndim_nc,time_dim,nsta_dim,ncid,iret
  integer,allocatable :: ndim_gb(:),dims_nc(:),var_dims_nc(:)
  character(len=30) :: varname_lc
  character(len=30), allocatable :: varname_gb(:)
  real(rkind),allocatable :: swild_lc(:),swild_gb(:) 
  type(cosine_diagnostic_variable), pointer :: dtv,drv
  logical :: lexist

  !1). initilize local diagnotic variable; 2) save variable values
  if(imode==0) then
    iflag=0
    if(associated(dlv)) then
      !dlv exists; search dlv for the variable
      dtv=>dlv
      do while(.True.)
        if(.not. associated(dtv)) exit
        if(trim(adjustl(dtv%varname)) .eq. trim(adjustl(varname))) then
          dtv%data(:,id)=rarray(1:ndim) !if variable found
          iflag=1; exit
        endif
        drv=>dtv; dtv=>dtv%next
      enddo
    endif

    !varable not found 
    if(iflag==0) then
      !initialize new variable
      allocate(dtv,stat=istat)
      if(istat/=0) call parallel_abort('failed in alloc. dtv')
      if(associated(dlv)) drv%next=>dtv
      if(.not. associated(dlv)) dlv=>dtv
      
      dtv%ndim=ndim
      dtv%varname=varname
      allocate(dtv%data(ndim,dl%nsta),stat=istat)
      if(istat/=0) call parallel_abort('failed in alloc. dtv%data')
      dtv%data(:,id)=rarray(1:ndim)

      dl%nvar=dl%nvar+1; dl%ndim=dl%ndim+ndim
      if(trim(adjustl(varname)).eq.'temp') dl%istat=1
    endif
  
  elseif(imode==1) then
    !------------------------------------------------------------------------------ 
    !compute total datasize for each station, and initilize dvars
    !------------------------------------------------------------------------------ 
    if(dg%istat==0) then
      !sync total nvar and ndim
      call mpi_allreduce(dl%nvar,dg%nvar,1,itype,MPI_MAX,comm,ierr)
      call mpi_allreduce(dl%ndim,dg%ndim,1,itype,MPI_MAX,comm,ierr)
      if(dg%nvar>0) dg%istat=1
      if(dg%istat==0) return

      if(myrank==0) then
        allocate(varname_gb(nproc),ndim_gb(nproc),dims_nc(dg%nvar),var_dims_nc(dg%nvar),stat=istat)
        if(istat/=0) call parallel_abort('failed in alloc. varname_gb')
        ndim_nc=0; dims_nc=0
      endif

      !initialize dgv
      if(associated(dlv)) dtv=>dlv
      do i=1,dg%nvar
        !get local ndim and varname
        if(associated(dtv)) then
          ndim_lc=dtv%ndim; varname_lc=dtv%varname
          dtv=>dtv%next
        else
          ndim_lc=0; varname_lc=''
        endif

        !initialize dgv
        call mpi_gather(ndim_lc,1,itype,ndim_gb,1,itype,0,comm,ierr)
        call mpi_gather(varname_lc,30,stype,varname_gb,30,stype,0,comm,ierr)
        if(myrank==0) then
          if(i==1) then
            allocate(drv,stat=istat)
            if(istat/=0) call parallel_abort('failed in alloc. drv')
            dgv=>drv
          else
            allocate(drv%next,stat=istat)
            if(istat/=0) call parallel_abort('failed in alloc. drv%next')
            drv=>drv%next
          endif

          !determine/check global ndim and varname
          ndim_lc=0; varname_lc=''
          do j=1,nproc
            if((varname_lc.eq.'').and.(varname_gb(j).ne.'')) varname_lc=varname_gb(j) 
            if(ndim_lc==0 .and. ndim_gb(j)/=0) ndim_lc=ndim_gb(j) 

            if((varname_lc.ne.'').and.(varname_gb(j).ne.'').and.(varname_lc.ne.varname_gb(j))) call parallel_abort('varname_lc: cosine') 
            if(ndim_lc/=0 .and. ndim_gb(j)/=0 .and. ndim_lc/=ndim_gb(j)) call parallel_abort('ndim_lc: cosine') 
          enddo
          drv%ndim=ndim_lc; drv%varname=varname_lc
          allocate(drv%data(ndim_lc,dg%nsta),stat=istat)
          if(istat/=0) call parallel_abort('failed in alloc. drv%data')

          !find variable with ndim>1
          if(ndim_lc>1) then
            istat=0
            do j=1,ndim_nc
              if(ndim_lc==dims_nc(j)) istat=1
            enddo
           
            !new dimension found
            if(istat==0) ndim_nc=ndim_nc+1; dims_nc(ndim_nc)=ndim_lc
          endif
        endif !if(myrank==0)

      enddo !i=1,dg%nvar

      !initialize station output file
      if(myrank==0) then
        inquire(file=trim(adjustl(out_dir(1:len_out_dir)//'cosine.nc')),exist=lexist)
        if(ihot==2 .and. lexist) then
          iret=nf90_open(trim(adjustl(out_dir(1:len_out_dir)//'cosine.nc')),NF90_WRITE,ncid)
          dg%ncid=ncid

          !get varid 
          iret=nf90_inq_varid(ncid,'time',dg%id_time)
          iret=nf90_inquire_dimension(ncid,dg%id_time,len=dg%it); dg%it=dg%it+1

          dtv=>dgv
          do i=1,dg%nvar
            iret=nf90_inq_varid(ncid,trim(adjustl(dtv%varname)),dtv%varid)
            dtv=>dtv%next
          enddo
        else
          iret=nf90_create(trim(adjustl(out_dir(1:len_out_dir)//'cosine.nc')),OR(NF90_NETCDF4,NF90_CLOBBER),ncid)
          dg%ncid=ncid

          !define dimension
          iret=nf90_def_dim(ncid,'time',NF90_UNLIMITED,time_dim)
          iret=nf90_def_dim(ncid,'nstation',dg%nsta,nsta_dim)
          do i=1,ndim_nc
            write(varname_lc,*)dims_nc(i)
            iret=nf90_def_dim(ncid,trim(adjustl(varname_lc)),dims_nc(i),var_dims_nc(i))
          enddo

          !define variables
          iret=nf90_def_var(ncid,'time',nf90_double,(/time_dim/),dg%id_time)
          iret=nf90_def_var(ncid,'x',nf90_double,(/nsta_dim/),id_x)
          iret=nf90_def_var(ncid,'y',nf90_double,(/nsta_dim/),id_y)
          iret=nf90_def_var(ncid,'z',nf90_double,(/nsta_dim/),id_z)
          iret=nf90_def_var(ncid,'ie',nf90_int,(/nsta_dim/),id_iep)

          dtv=>dgv
          do i=1,dg%nvar
            if(dtv%ndim==1) then
              iret=nf90_def_var(ncid,trim(adjustl(dtv%varname)),nf90_FLOAT,(/time_dim,nsta_dim/),dtv%varid)
            elseif(dtv%ndim>1) then
              do j=1,ndim_nc; if(dims_nc(j)==dtv%ndim) exit; enddo
              iret=nf90_def_var(ncid,trim(adjustl(dtv%varname)),nf90_FLOAT,(/time_dim,var_dims_nc(j),nsta_dim/),dtv%varid)
            endif
            dtv=>dtv%next
          enddo
          iret=nf90_enddef(ncid)

          !put x,y,z and iep
          iret=nf90_put_var(dg%ncid,id_iep,dg%iep,start=(/1/),count=(/dg%nsta/))
          iret=nf90_put_var(dg%ncid,id_x,dg%x,start=(/1/),count=(/dg%nsta/))
          iret=nf90_put_var(dg%ncid,id_y,dg%y,start=(/1/),count=(/dg%nsta/))
          iret=nf90_put_var(dg%ncid,id_z,dg%z,start=(/1/),count=(/dg%nsta/))
        endif!ihot

        deallocate(varname_gb,ndim_gb,dims_nc,var_dims_nc) 
      endif !myrank
    endif !dg%istat==0

    !------------------------------------------------------------------------------ 
    !allocate data for storing all diagnostic values 
    !------------------------------------------------------------------------------ 
    !check nvar and ndim
    !if(dl%nsta/=0) then !this can happen if dry_ie(i)==1
    !  if(dl%nvar/=dg%nvar) call parallel_abort('dl%nvar/=dg%nvar')
    !  if(dl%ndim/=dg%ndim) call parallel_abort('dl%ndim/=dg%ndim')
    !endif
  
    allocate(swild_lc(dg%ndim*dl%nsta),stat=istat)
    if(istat/=0) call parallel_abort('failed in alloc. swild_lc')
    if(myrank==0) then
      allocate(swild_gb(dg%ndim*dg%nsta),stat=istat)
      if(istat/=0) call parallel_abort('failed in alloc. swild_gb')
    endif

    !assemble data on each rank 
    swild_lc=-9999; m=0 
    do i=1,dl%nsta
      if(.not. associated(dlv)) cycle
      dtv=>dlv
      do j=1,dl%nvar
        do n=1,dtv%ndim
          m=m+1
          swild_lc(m)=dtv%data(n,i)
        enddo !n
        dtv=>dtv%next
      enddo !j=1,dl%nvar
    enddo !i=1,dl%nsta

    !pass all data to myrank=0
    call mpi_gatherv(swild_lc,dl%nsta*dg%ndim,rtype,swild_gb,dg%nstas*dg%ndim,dg%displ*dg%ndim,rtype,0,comm,ierr)
    if(myrank==0) then
      m=0
      do i=1,dg%nsta
        dtv=>dgv
        do j=1,dg%nvar
          do n=1,dtv%ndim
            m=m+1; dtv%data(n,dg%sids(i))=swild_gb(m) 
          enddo
          dtv=>dtv%next
        enddo !j1,dg%nvar
      enddo !i=1,dg%nsta

      !write station output
      iret=nf90_put_var(dg%ncid,dg%id_time,(/dg%time/),start=(/dg%it/),count=(/1/))
      dtv=>dgv
      do i=1,dg%nvar
        if(dtv%ndim==1) then
          iret=nf90_put_var(dg%ncid,dtv%varid,dtv%data(:,1),start=(/dg%it,1/),count=(/1,dg%nsta/))
        elseif(dtv%ndim>1) then
          iret=nf90_put_var(dg%ncid,dtv%varid,dtv%data(:,:),start=(/dg%it,1,1/),count=(/1,dtv%ndim,dg%nsta/))
        endif
        dtv=>dtv%next
      enddo
      iret=nf90_sync(dg%ncid) 
    endif !myrank

    deallocate(swild_lc)
    if(myrank==0) deallocate(swild_gb)
  endif !imode
    
end subroutine cosine_output
