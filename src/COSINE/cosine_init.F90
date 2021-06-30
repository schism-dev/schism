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
          & iout_cosine,nspool_cosine,ico2s,ispm,spm0,ised,nsedS2,nsedDN,nsedDSi
  namelist /CORE/ gmaxs1,gmaxs2,pis1,pis2,kno3s1,knh4s1,kpo4s1,kco2s1,kno3s2,&
          & knh4s2,kpo4s2,kco2s2,ksio4s2,kns1,kns2,alpha1,alpha2,beta,ak1,ak2,&
          & ak3,gammas1,gammas2,beta1,beta2,kgz1,kgz2,rho1,rho2,rho3,gamma1,&
          & gamma2,gammaz,kex1,kex2,wss2,wsdn,wsdsi,si2n,p2n,o2no,o2nh,c2n,&
          & kox,kmdn1,kmdn2,kmdsi1,kmdsi2,gamman,TR,pco2a
  namelist /MISC/ iws,NO3c,ws1,ws2,iclam,deltaZ,kcex,Nperclam,Wclam,Fclam,&
          & nclam0,nsedS2,psedS2,rsedS2,rsedS2m,nsedDN,psedDN,rsedDN,rsedDNm,&
          & nsedDSi,psedDSi,rsedDSi,rsedDSim

  !initialize parameter values
  idelay=0; ndelay=7; ibgraze=0; idapt=0; alpha_corr=1.25; zeptic=10.0; iz2graze=1
  iout_cosine=0; nspool_cosine=60; ico2s=0; ispm=0; spm0=20.0; ised=1; nsedS2=2; 
  nsedDN=2; nsedDSi=1; gmaxs1=3.0; gmaxs2=2.5; pis1=1.5; pis2=1.5; kno3s1=1.0; 
  knh4s1=0.15; kpo4s1=0.1; kco2s1=50.0; kno3s2=3.0; knh4s2=0.45; kpo4s2=0.1;
  kco2s2=50.0; ksio4s2=4.5; kns1=0.0; kns2=0.0; alpha1=0.1; alpha2=0.1; beta=0.0; 
  ak1=0.75; ak2=0.03; ak3=0.066; gammas1=0.5; gammas2=0.3; beta1=0.75; beta2=0.5; 
  kgz1=0.5; kgz2=0.25; rho1=0.6; rho2=0.3; rho3=0.1; gamma1=0.75; gamma2=0.75; 
  gammaz=0.05; kex1=0.2; kex2=0.3; wss2=0.25; wsdn=0.5; wsdsi=0.5; si2n=1.2; p2n=0.0625; 
  o2no=8.625; o2nh=6.625; c2n=7.3; kox=30.0; kmdn1=0.009; kmdn2=0.075; kmdsi1=0.0114; 
  kmdsi2=0.015; gamman=0.07; TR=20.0; pco2a=400.0; iws=0; NO3c=2.0; ws1=2.5; ws2=2.0
  iclam=0; deltaZ=1.0; kcex=0.002; Nperclam=0.39032; Wclam=5.45e-3; Fclam=40.0; 
  nclam0=2000; psedS2=0.0; rsedS2=0.0; rsedS2m=0.1; psedDN=0.0; rsedDN=0.0; rsedDNm=0.1 
  psedDSi=0.0; rsedDSi=4e-3; rsedDSim=0.1

  !read parameter values
  open(31,file=in_dir(1:len_in_dir)//'cosine.nml',delim='apostrophe',status='old')
  read(31,nml=MARCO); read(31,nml=CORE)
  
  !allocate sediment variables
  if(ised==1) then
    nsed=nsedS2+nsedDN+nsedDSi
    allocate(psedS2(nsedS2),rsedS2(nsedS2),rsedS2m(nsedS2),psedDN(nsedDN), &
          & rsedDN(nsedDN),rsedDNm(nsedDN),psedDSi(nsedDSi),rsedDSi(nsedDSi), &
          & rsedDSim(nsedDSi),rsed(nsed),rsedm(nsed),sedcon(nsed,nea),  &
          & sedrate(nsed,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. psedDN')
  endif

  read(31,nml=MISC)
  close(31)

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
    m=0
    do i=1,nsedS2
      m=m+1 
      rsed(m)=rsedS2(i)
      rsedm(m)=rsedS2m(i)
    enddo
    do i=1,nsedDN
      m=m+1 
      rsed(m)=rsedDN(i)
      rsedm(m)=rsedDNm(i)
    enddo
    do i=1,nsedDSi
      m=m+1 
      rsed(m)=rsedDSi(i)
      rsedm(m)=rsedDSim(i)
    enddo

    !initialize sediment variables
    !todo, include sedcon, sedrate in hotstart.in; update to sediment flux model
    sedcon=0.d0; sedrate=0.d0
    !sedcon(1,:)=1; sedcon(2,:)=100; sedcon(3,:)=1; sedcon(4,:)=100; sedcon(5,:)=10
    if(ihot==0) then !temporary fix, need to update
      do i=1,nsed
        write(snum,*)i
        call read_gr3_prop('sedcon_'//trim(adjustl(snum)),-999.d0,sedcon(i,:),nea)
        call read_gr3_prop('sedrate_'//trim(adjustl(snum)),-999.d0,sedrate(i,:),nea)
      enddo
    endif!
  endif !ised
 
  return
end subroutine read_cosine_param

subroutine read_cosine_stainfo
!---------------------------------------------------------------------------
!output cosine parameters to check
!---------------------------------------------------------------------------
  use cosine_mod, only : nsta,ista,depsta,stanum,nspool_cosine
  use schism_glbl, only : rkind,dt,ihot,ne,i34,xnd,ynd,elnode, &
      & in_dir,out_dir,len_in_dir,len_out_dir
  use schism_msgp, only : myrank,nproc,parallel_abort
  use cosine_misc, only : pt_in_poly
  implicit none

  !local variables
  integer,parameter :: maxsta=10000,maxl=100 !maximum station
  integer :: i,j,istat,nstation,nodel(3),inside,id,iflag,mid,msta,nstai(ne),stanumi(maxl,ne)
  real(rkind) :: slx(maxsta),sly(maxsta),sdep(maxsta),x(3),y(3),arco(3),depstai(maxl,ne)
  character(len=4) :: fn
  logical :: lexist

  !read station info.
  open(450,file=in_dir(1:len_in_dir)//'cstation.in',status='old')
  read(450,*)
  read(450,*)nstation
  do i=1,nstation
    read(450,*)j,slx(i),sly(i),sdep(i) 
  enddo
  close(450)

  !alloc.
  allocate(ista(ne),stat=istat) 
  if(istat/=0) call parallel_abort('failure in alloc. ista')

  !determine the elements with values to be checked
  id=0; ista=0; nstai=0;
  msta=-100; depstai=-9999
  do i=1,ne
    iflag=0
    do j=1,nstation
      x=xnd(elnode(1:i34(i),i))
      y=ynd(elnode(1:i34(i),i))
      call pt_in_poly(i34(i),x,y,slx(j),sly(j),inside,arco,nodel)
      if(inside==1) then
        if(ista(i)==0) then
          id=id+1
          ista(i)=id
        endif
        nstai(id)=nstai(id)+1
        depstai(nstai(id),id)=sdep(j)
        stanumi(nstai(id),id)=j
        msta=max(msta,nstai(id))
      endif 
    enddo !j
  enddo !i
  mid=id !number of elements

  if(mid==0) return 

  !alloc.
  allocate(nsta(mid),depsta(msta,mid),stanum(msta,mid),stat=istat) 
  if(istat/=0) call parallel_abort('failure in alloc. nsta')

  nsta=0; depsta=-9999
  do i=1,mid
    nsta(i)=nstai(i)
    do j=1,nsta(i)
      depsta(j,i)=depstai(j,i)
      stanum(j,i)=stanumi(j,i)
    enddo
  enddo

  !open a station output file
  write(fn,'(i4.4)')myrank
  inquire(file=out_dir(1:len_out_dir)//'cstation_'//fn//'.out',exist=lexist)
  if(ihot<=1.or.(ihot==2.and.(.not.lexist))) then
    open(451,file=out_dir(1:len_out_dir)//'cstation_'//fn//'.out',form='unformatted',status='replace')
    write(451)sum(nsta),dt*nspool_cosine
  elseif(ihot==2.and.lexist) then
    open(451,file=out_dir(1:len_out_dir)//'cstation_'//fn//'.out',form='unformatted',access='append',status='old')
  else
    call parallel_abort('unknown ihot, CoSiNE')
  endif

  return
end subroutine read_cosine_stainfo
