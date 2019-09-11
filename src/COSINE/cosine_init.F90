!Routines & functions
!cosine_init: allocate and initilize variables
!read_cosine_param: read cosine parameters
!read_cosine_stainfo: read info. for station outputs
!check_cosine_param: output cosine parameters for check
!pt_in_poly: determine if point is in a polygon

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
          & salt(nvrt),bgraze(nea),SPM(nvrt,nea),bio(nvrt,ntrc),bio0(nvrt,ntrc),qcos(ntrc),sqcos(nvrt,ntrc),&
          & mS2(ndelay,nvrt,nea),mDN(ndelay,nvrt,nea),mZ1(ndelay,nvrt,nea),mZ2(ndelay,nvrt,nea),&
          & sS2(nvrt,nea),sDN(nvrt,nea),sZ1(nvrt,nea),sZ2(nvrt,nea),nstep(nvrt,nea), &
          & nclam(nea), stat=istat) 
  if(istat/=0) call parallel_abort('failure in alloc. mS2')

  !initialize
  NO3=0.0;  NH4=0.0; SiO4=0.0; S1=0.0;   S2=0.0;   Z1=0.0;  Z2=0.0
  DN=0.0;   DSi=0.0; PO4=0.0;  DOX=0.0;  CO2=0.0;  ALK=0.0; temp=0.0
  salt=0.0; SPM=0.0; bio=0.0;  bio0=0.0; qcos=0.0; sqcos=0.0 
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
 &in_dir,out_dir,len_in_dir,len_out_dir
  use schism_msgp, only : myrank,parallel_abort
  use cosine_mod
  use misc_modules
  implicit none

  !local variables
  integer :: i,j,k,m,negb,npgb,nd,itmp,itmp1(1),itmp2(1,1),istat
  real(rkind) :: xtmp,ytmp,tSPM,tSPMs(npa),rtmp,rtmp1(1),rtmp2(1,1)
  integer :: tnclam,nclams(npa)
  character(len=2) :: stmp

  !read switches and macro parameters
  call get_param('cosine.in','niter',1,niter,rtmp,stmp)
  call get_param('cosine.in','idelay',1,idelay,rtmp,stmp)
  if(idelay==1) call get_param('cosine.in','ndelay',1,ndelay,rtmp,stmp)
  call get_param('cosine.in','ibgraze',1,ibgraze,rtmp,stmp)
  call get_param('cosine.in','idapt',1,idapt,rtmp,stmp)
  call get_param('cosine.in','iz2graze',1,iz2graze,rtmp,stmp)
  call get_param('cosine.in','iout_cosine',1,iout_cosine,rtmp,stmp)
  call get_param('cosine.in','nspool_cosine',1,nspool_cosine,rtmp,stmp)
  call get_param('cosine.in','ico2s',1,ico2s,rtmp,stmp)
  call get_param('cosine.in','ispm',1,ispm,rtmp,stmp)
  if(ispm==0) then
    call get_param('cosine.in','spm0',2,itmp,spm0,stmp)
  endif
  call get_param('cosine.in','icheck',1,icheck,rtmp,stmp)
  call get_param('cosine.in','ised',1,ised,rtmp,stmp)
  call get_param('cosine.in','iws',1,iws,rtmp,stmp)
  if(iws/=0) then
    call get_param('cosine.in','NO3c',2,itmp,NO3c,stmp)
    call get_param('cosine.in','ws1',2,itmp,ws1,stmp)
    call get_param('cosine.in','ws2',2,itmp,ws2,stmp)
  endif
  call get_param('cosine.in','iclam',1,iclam,rtmp,stmp)
  if(iclam/=0) then
    call get_param('cosine.in','deltaZ',2,itmp,deltaZ,stmp)
    call get_param('cosine.in','kcex',2,itmp,kcex,stmp)
    call get_param('cosine.in','Nperclam',2,itmp,Nperclam,stmp)
    call get_param('cosine.in','Wclam',2,itmp,Wclam,stmp)
    call get_param('cosine.in','Fclam',2,itmp,Fclam,stmp)
    if(iclam==1) then
      call get_param('cosine.in','nclam0',1,nclam0,rtmp,stmp)
    endif
  endif

  !read cosine kinetics parameters
  !phytoplankton
  call get_param('cosine.in','gmaxs1',2,itmp,gmaxs1,stmp)
  call get_param('cosine.in','alpha1',2,itmp,alpha1,stmp)
  call get_param('cosine.in','pis1',2,itmp,pis1,stmp)
  call get_param('cosine.in','kno3s1',2,itmp,kno3s1,stmp)
  call get_param('cosine.in','knh4s1',2,itmp,knh4s1,stmp)
  call get_param('cosine.in','kpo4s1',2,itmp,kpo4s1,stmp)
  call get_param('cosine.in','kco2s1',2,itmp,kco2s1,stmp)
  call get_param('cosine.in','kns1',2,itmp,kns1,stmp)
  call get_param('cosine.in','gammas1',2,itmp,gammas1,stmp)

  call get_param('cosine.in','gmaxs2',2,itmp,gmaxs2,stmp)
  call get_param('cosine.in','alpha2',2,itmp,alpha2,stmp)
  call get_param('cosine.in','pis2',2,itmp,pis2,stmp)
  call get_param('cosine.in','kno3s2',2,itmp,kno3s2,stmp)
  call get_param('cosine.in','knh4s2',2,itmp,knh4s2,stmp)
  call get_param('cosine.in','kpo4s2',2,itmp,kpo4s2,stmp)
  call get_param('cosine.in','kco2s2',2,itmp,kco2s2,stmp)
  call get_param('cosine.in','kns2',2,itmp,kns2,stmp)
  call get_param('cosine.in','gammas2',2,itmp,gammas2,stmp)
  call get_param('cosine.in','ksio4s2',2,itmp,ksio4s2,stmp)

  call get_param('cosine.in','ak1',2,itmp,ak1,stmp)
  call get_param('cosine.in','ak2',2,itmp,ak2,stmp)
  call get_param('cosine.in','ak3',2,itmp,ak3,stmp)
  call get_param('cosine.in','alpha_corr',2,itmp,alpha_corr,stmp)
  call get_param('cosine.in','zeptic',2,itmp,zeptic,stmp)
  call get_param('cosine.in','beta',2,itmp,beta,stmp)

  !zooplankton
  call get_param('cosine.in','kex1',2,itmp,kex1,stmp)
  call get_param('cosine.in','gamma1',2,itmp,gamma1,stmp)
  call get_param('cosine.in','kex2',2,itmp,kex2,stmp)
  call get_param('cosine.in','gamma2',2,itmp,gamma2,stmp)

  call get_param('cosine.in','beta1',2,itmp,beta1,stmp)
  call get_param('cosine.in','beta2',2,itmp,beta2,stmp)
  call get_param('cosine.in','kgz1',2,itmp,kgz1,stmp)
  call get_param('cosine.in','kgz2',2,itmp,kgz2,stmp)
  call get_param('cosine.in','rho1',2,itmp,rho1,stmp)
  call get_param('cosine.in','rho2',2,itmp,rho2,stmp)
  call get_param('cosine.in','rho3',2,itmp,rho3,stmp)

  call get_param('cosine.in','gammaz',2,itmp,gammaz,stmp)

  !other
  call get_param('cosine.in','kox',2,itmp,kox,stmp)
  call get_param('cosine.in','kbmdn',2,itmp,kbmdn,stmp)
  call get_param('cosine.in','kmdn1',2,itmp,kmdn1,stmp)
  call get_param('cosine.in','kmdn2',2,itmp,kmdn2,stmp)
  call get_param('cosine.in','kbmdsi',2,itmp,kbmdsi,stmp)
  call get_param('cosine.in','kmdsi1',2,itmp,kmdsi1,stmp)
  call get_param('cosine.in','kmdsi2',2,itmp,kmdsi2,stmp)
  call get_param('cosine.in','TR',2,itmp,TR,stmp)
  call get_param('cosine.in','gamman',2,itmp,gamman,stmp)
  call get_param('cosine.in','pco2a',2,itmp,pco2a,stmp)
  call get_param('cosine.in','wss2',2,itmp,wss2,stmp)
  call get_param('cosine.in','wsdn',2,itmp,wsdn,stmp)
  call get_param('cosine.in','wsdsi',2,itmp,wsdsi,stmp)
  call get_param('cosine.in','si2n',2,itmp,si2n,stmp)
  call get_param('cosine.in','p2n',2,itmp,p2n,stmp)
  call get_param('cosine.in','o2no',2,itmp,o2no,stmp)
  call get_param('cosine.in','o2nh',2,itmp,o2nh,stmp)
  call get_param('cosine.in','c2n',2,itmp,c2n,stmp)

  !read in station info. 
  if(iout_cosine==1) call read_cosine_stainfo

  !read bottom grazing information
  if(ibgraze==1) then
    bgraze=0.0
    open(454,file=in_dir(1:len_in_dir)//'bgraze.prop',status='old')
    do i=1,ne_global
      read(454,*)itmp,ytmp
      if(iegl(itmp)%rank==myrank) then
        bgraze(iegl(itmp)%id)=ytmp 
      endif
    enddo
    close(454)
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
    open(452,file=in_dir(1:len_in_dir)//'nclam.gr3',status='old')
    read(452,*);read(452,*)negb,npgb
    if(negb/=ne_global.or.npgb/=np_global) call parallel_abort('Check nclam.gr3')
    do i=1,np_global
      read(452,*)itmp,xtmp,ytmp,tnclam
      if(ipgl(itmp)%rank==myrank) then
        nclams(ipgl(itmp)%id)=tnclam
      endif
    enddo
    close(452)

    do i=1,nea
      do j=1,i34(i)
        nd=elnode(j,i)
        nclam(i)=nclam(i)+nclams(nd)
      enddo
      nclam(i)=nclam(i)/i34(i)
    enddo
  elseif (iclam==3) then
    open(456,file=in_dir(1:len_in_dir)//'nclam.th',status='old') 
    time_cosine(3)=-999.0
  endif 

  !read SPM information
  SPM=0.0
  if(ispm==0) then !constant
    SPM=spm0
  elseif(ispm==1) then !spatial varying
    open(452,file=in_dir(1:len_in_dir)//'SPM.gr3',status='old')
    read(452,*); read(452,*)negb,npgb
    if(negb/=ne_global.or.npgb/=np_global) call parallel_abort('Check SPM.gr3')
    do i=1,np_global
      read(452,*)itmp,xtmp,ytmp,tSPM
      if(ipgl(itmp)%rank==myrank) then
        tSPMs(ipgl(itmp)%id)=tSPM
      endif
    enddo !i
    close(452)

    do i=1,nea
      do j=1,i34(i)
        nd=elnode(j,i)
        SPM(:,i)=SPM(:,i)+tSPMs(nd)/i34(i)
      enddo !j
    enddo !i
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
    call get_param('cosine.in','nsedS2',1,nsedS2,rtmp,stmp)
    call get_param('cosine.in','nsedDN',1,nsedDN,rtmp,stmp)
    call get_param('cosine.in','nsedDSi',1,nsedDSi,rtmp,stmp)
    nsed=nsedS2+nsedDN+nsedDSi
    allocate(psedS2(nsedS2),rsedS2(nsedS2),rsedS2m(nsedS2),psedDN(nsedDN), &
          & rsedDN(nsedDN),rsedDNm(nsedDN),psedDSi(nsedDSi),rsedDSi(nsedDSi), &
          & rsedDSim(nsedDSi),rsed(nsed),rsedm(nsed),sedcon(nsed,nea),  &
          & sedrate(nsed,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. psedDN')
    call get_param_1D('cosine.in','psedS2',2,itmp1,psedS2,stmp,nsedS2)
    call get_param_1D('cosine.in','rsedS2',2,itmp1,rsedS2,stmp,nsedS2)
    call get_param_1D('cosine.in','rsedS2m',2,itmp1,rsedS2m,stmp,nsedS2)
    call get_param_1D('cosine.in','psedDN',2,itmp1,psedDN,stmp,nsedDN)
    call get_param_1D('cosine.in','rsedDN',2,itmp1,rsedDN,stmp,nsedDN)
    call get_param_1D('cosine.in','rsedDNm',2,itmp1,rsedDNm,stmp,nsedDN)
    call get_param_1D('cosine.in','psedDSi',2,itmp1,psedDSi,stmp,nsedDSi)
    call get_param_1D('cosine.in','rsedDSi',2,itmp1,rsedDSi,stmp,nsedDSi)
    call get_param_1D('cosine.in','rsedDSim',2,itmp1,rsedDSim,stmp,nsedDSi)

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

    !initialize
    !ZG, include sedcon, sedrate in hotstart.in
    sedcon=0.d0
    sedrate=0.d0
  endif !ised
 
  !output cosine parameter 
  if(icheck==1) call check_cosine_param

  return
end subroutine read_cosine_param

subroutine read_cosine_stainfo
!---------------------------------------------------------------------------
!output cosine parameters to check
!---------------------------------------------------------------------------
  use cosine_mod, only : nsta,ista,depsta,stanum,nspool_cosine
  use schism_glbl, only : rkind,dt,ihot,ne,i34,xnd,ynd,elnode, &
 &in_dir,out_dir,len_in_dir,len_out_dir
  use schism_msgp, only : myrank,nproc,parallel_abort
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

subroutine check_cosine_param
!---------------------------------------------------------------------------
!output cosine parameters to check
!---------------------------------------------------------------------------
  use cosine_mod
  use schism_msgp, only : myrank,parallel_abort
  use schism_glbl, only : in_dir,out_dir,len_in_dir,len_out_dir
  implicit none

  !local variables
  integer :: i,j
  if(myrank==0) then
    open(31,file=out_dir(1:len_out_dir)//'cosine_param.out',status='replace')
    write(31,*) 'Cosine Parameters used in the model'

    write(31,*)
    write(31,*)'!------switches and macro parameters-----'
    write(31,'(a10,i5)')'niter = ',niter
    write(31,'(a10,i5)')'idelay = ',idelay
    write(31,'(a10,i5)')'ibgraze = ',ibgraze
    write(31,'(a10,i5)')'idapt = ',idapt
    write(31,'(a10,i5)')'iz2graze = ',iz2graze
    write(31,'(a10,i5)')'iout_cosine = ',iout_cosine
    write(31,'(a10,i5)')'nspool_cosine = ',nspool_cosine
    write(31,'(a10,i5)')'ico2s = ',ico2s
    write(31,'(a10,i5)')'ispm = ',ispm
    if(ispm==0) then
      write(31,'(a10,f12.3)')'spm0 = ',spm0
    endif
    write(31,'(a10,i5)')'icheck = ',icheck
    write(31,'(a10,i5)')'ised = ',ised
    write(31,'(a10,i5)')'iws = ',iws
    if(iws/=0) then
      write(31,'(a10,f12.3)')'NO3c = ',NO3c
      write(31,'(a10,f12.3)')'ws1 = ',ws1
      write(31,'(a10,f12.3)')'ws2 = ',ws2
    endif
    write(31,'(a10,i5)')'iclam = ',iclam
    if(iclam/=0) then
      write(31,'(a10,f12.3)')'deltaZ = ',deltaZ
      write(31,'(a10,f12.3)')'kcex = ',kcex
      write(31,'(a10,f12.3)')'Nperclam = ',Nperclam
      write(31,'(a10,f12.3)')'Wclam = ',Wclam
      write(31,'(a10,f12.3)')'Fclam = ',Fclam
      if(iclam==1) then
        write(31,'(a10,i5)')'nclam0 = ',nclam0
      endif
    endif

    write(31,*)
    write(31,*)'!--------cosine kinetics parameters------'
    write(31,'(a10,i5)')'icheck = ',icheck
    !write(31,'(a10,i5)')' =',

    write(31,*)
    write(31,*)'!--------cosine kinetics parameters------'
    write(31,*)
    write(31,*)'!phytoplankton'
    write(31,'(a10,100(f8.5 x))')'gmaxs1 = ',gmaxs1
    write(31,'(a10,100(f8.5 x))')'alpha1 = ',alpha1
    write(31,'(a10,100(f8.5 x))')'pis1 = ',pis1
    write(31,'(a10,100(f8.5 x))')'kno3s1 = ',kno3s1
    write(31,'(a10,100(f8.5 x))')'knh4s1 = ',knh4s1
    write(31,'(a10,100(f8.5 x))')'kpo4s1 = ',kpo4s1
    write(31,'(a10,100(f18.9 x))')'kco2s1 = ',kco2s1
    write(31,'(a10,100(f18.9 x))')'kns1 = ',kns1
    write(31,'(a10,100(f8.5 x))')'gammas1 = ',gammas1

    write(31,'(a10,100(f8.5 x))')'gmaxs2 = ',gmaxs2
    write(31,'(a10,100(f8.5 x))')'alpha2 = ',alpha2
    write(31,'(a10,100(f8.5 x))')'pis2 = ',pis2
    write(31,'(a10,100(f8.5 x))')'kno3s2 = ',kno3s2
    write(31,'(a10,100(f8.5 x))')'knh4s2 = ',knh4s2
    write(31,'(a10,100(f8.5 x))')'kpo4s2 = ',kpo4s2
    write(31,'(a10,100(f18.9 x))')'kco2s2 = ',kco2s2
    write(31,'(a10,100(f18.9 x))')'kns2 = ',kns2
    write(31,'(a10,100(f8.5 x))')'gammas2 = ',gammas2
    write(31,'(a10,100(f8.5 x))')'ksio4s2 = ',ksio4s2

    write(31,'(a10,100(f8.5 x))')'ak1 = ',ak1
    write(31,'(a10,100(f8.5 x))')'ak2 = ',ak2
    write(31,'(a10,100(f8.5 x))')'ak3 = ',ak3
    write(31,'(a10,100(f8.5 x))')'alpha_corr = ',alpha_corr
    write(31,'(a10,100(f8.5 x))')'zeptic = ', zeptic
    write(31,'(a10,100(f8.5 x))')'beta = ',beta

    write(31,*)
    write(31,*)'!phytoplankton'
    write(31,'(a10,100(f8.5 x))')'kex1 = ',kex1
    write(31,'(a10,100(f8.5 x))')'gamma1 = ',gamma1
    write(31,'(a10,100(f8.5 x))')'kex2 = ',kex2
    write(31,'(a10,100(f8.5 x))')'gamma2 = ',gamma2

    write(31,'(a10,100(f8.5 x))')'beta1 = ',beta1
    write(31,'(a10,100(f8.5 x))')'beta2 = ',beta2
    write(31,'(a10,100(f8.5 x))')'kgz1 = ',kgz1
    write(31,'(a10,100(f8.5 x))')'kgz2 = ',kgz2
    write(31,'(a10,100(f8.5 x))')'rho1 = ',rho1
    write(31,'(a10,100(f8.5 x))')'rho2 = ',rho2
    write(31,'(a10,100(f8.5 x))')'rho3 = ',rho3
    write(31,'(a10,100(f8.5 x))')'gammaz = ',gammaz

    write(31,*)
    write(31,*)'!other'
    write(31,'(a10,100(f8.5 x))')'kox = ',kox
    write(31,'(a10,100(f8.5 x))')'kbmdn = ',kbmdn
    write(31,'(a10,100(f8.5 x))')'kmdn1 = ',kmdn1
    write(31,'(a10,100(f8.5 x))')'kmdn2 = ',kmdn2
    write(31,'(a10,100(f8.5 x))')'kbmdsi = ',kbmdsi
    write(31,'(a10,100(f8.5 x))')'kmdsi1 = ',kmdsi1
    write(31,'(a10,100(f8.5 x))')'kmdsi2 = ',kmdsi2
    write(31,'(a10,100(f8.5 x))')'TR = ',TR
    write(31,'(a10,100(f8.5 x))')'gamman = ',gamman
    write(31,'(a10,100(f18.9 x))')'pco2a = ',pco2a
    write(31,'(a10,100(f8.5 x))')'wss2 = ',wss2
    write(31,'(a10,100(f8.5 x))')'wsdn = ',wsdn
    write(31,'(a10,100(f8.5 x))')'wsdsi = ',wsdsi
    write(31,'(a10,100(f8.5 x))')'si2n = ',si2n
    write(31,'(a10,100(f8.5 x))')'p2n = ',p2n
    write(31,'(a10,100(f8.5 x))')'o2no = ',o2no
    write(31,'(a10,100(f8.5 x))')'o2nh = ',o2nh
    write(31,'(a10,100(f8.5 x))')'c2n = ',c2n

    !write(31,'(a10,100(f8.5 x))')' = ',
    close(31)
  endif !myrank
  
  return
end subroutine check_cosine_param

subroutine pt_in_poly(i34,x,y,xp,yp,inside,arco,nodel)
!---------------------------------------------------------------------------
!subroutine from Utility/UtilLib
!---------------------------------------------------------------------------
!     (Single-precision) Routine to perform point-in-polygon
!     (triangle/quads) test and if it's inside, calculate the area coord.
!     (for quad, split it into 2 triangles and return the 3 nodes and
!     area coord.)
!     Inputs:
!            i34: 3 or 4 (type of elem)
!            x(i34),y(i34): coord. of polygon/elem. (counter-clockwise)
!            xp,yp: point to be tested
!     Outputs:
!            inside: 0, outside; 1, inside
!            arco(3), nodel(3) : area coord. and 3 local node indices (valid only if inside)
      implicit real*8(a-h,o-z)
      integer, intent(in) :: i34
      real(kind=8), intent(in) :: x(i34),y(i34),xp,yp
      integer, intent(out) :: inside,nodel(3)
      real(kind=8), intent(out) :: arco(3)

      !Local
      integer :: list(3)
      real(kind=8) :: ar(2),swild(2,3)

      !Areas
      ar(1)=signa(x(1),x(2),x(3),y(1),y(2),y(3))
      ar(2)=0 !init
      if(i34==4) ar(2)=signa(x(1),x(3),x(4),y(1),y(3),y(4))
      if(ar(1)<=0.or.i34==4.and.ar(2)<=0) then
        print*, 'Negative area:',i34,ar,x,y
        stop
      endif

      inside=0
      do m=1,i34-2 !# of triangles
        if(m==1) then
          list(1:3)=(/1,2,3/) !local indices
        else !quads
          list(1:3)=(/1,3,4/)
        endif !m
        aa=0
        do j=1,3
          j1=j+1
          j2=j+2
          if(j1>3) j1=j1-3
          if(j2>3) j2=j2-3
          swild(m,j)=signa(x(list(j1)),x(list(j2)),xp,y(list(j1)),y(list(j2)),yp) !temporary storage
          aa=aa+abs(swild(m,j))
        enddo !j=1,3

        ae=abs(aa-ar(m))/ar(m)
        if(ae<=1.e-5) then
          inside=1
          nodel(1:3)=list(1:3)
          arco(1:3)=swild(m,1:3)/ar(m)
          arco(1)=max(0.,min(1.,arco(1)))
          arco(2)=max(0.,min(1.,arco(2)))
          if(arco(1)+arco(2)>1) then
            arco(3)=0
            arco(2)=1-arco(1)
          else
            arco(3)=1-arco(1)-arco(2)
          endif
          exit
        endif
      enddo !m

end subroutine pt_in_poly

subroutine get_param_1D(fname,varname,vartype,ivar,rvar,svar,idim1)
!--------------------------------------------------------------------
!Read a one-Dimensional CoSiNE parameter
!--------------------------------------------------------------------
  use schism_glbl, only : rkind,errmsg
  use schism_msgp, only : parallel_abort,myrank
  use misc_modules
  implicit none

  character(*),intent(in) :: fname
  character(*),intent(in) :: varname
  integer,intent(in) :: vartype
  integer,intent(in) :: idim1
  integer,intent(out) :: ivar(idim1)
  real(rkind),intent(out) :: rvar(idim1)
  character(len=2),intent(out) :: svar

  !local variables
  integer :: itmp,iarray(10000)
  real(rkind) :: rtmp,rarray(10000)
  character(len=2) :: stmp

  svar='  '
  if(vartype==1) then  !read 1D integer array
    call get_param(fname,varname,3,itmp,rtmp,stmp,ndim1=idim1,iarr1=iarray)
    ivar=iarray(1:idim1)
  elseif(vartype==2) then !read 1D float array
    call get_param(fname,varname,4,itmp,rtmp,stmp,ndim1=idim1,arr1=rarray)
    rvar=rarray(1:idim1)
  else
    write(errmsg,*)'unknown vartype :',varname
    call parallel_abort(errmsg)
  endif

  return
end subroutine get_param_1D

