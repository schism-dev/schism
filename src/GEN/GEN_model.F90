module GEN_mod
!------------------------------------------------------------
!define GEN model variables
!------------------------------------------------------------
  use schism_glbl,only: rkind,nvrt,nea
  implicit none

  !constant parameter
  real(rkind),save,target :: dz_flux,WRea,KT,TR

  !spatial varying parameters
  real(rkind),save,target,allocatable :: dr(:)

  !for reading parameters
  type :: gen_param
    character(len=30) :: varname  !parameter name
    integer :: nv=0               !number of values
    integer,allocatable,dimension(:,:) :: istat
    real(rkind),pointer :: p=>null()                 !param. of single value
    real(rkind),pointer,dimension(:) :: p1=>null()   !param. of 1D array
  end type gen_param
  type(gen_param),save,target,allocatable,dimension(:) :: sp

end module GEN_mod

subroutine GEN_model
!------------------------------------------------------------
!GEN model kinetics
!------------------------------------------------------------
  use GEN_mod
  use schism_glbl, only : rkind,dt,elnode,i34,nea,tr_el,nvrt,kbe,ze,windx,windy
  implicit none

  !local variables
  integer :: i,j,k,ie,m,kb
  real(rkind) :: z1,z2,s
  real(rkind) :: dtw,zid(nvrt),dz(nvrt),tdep,srat(nvrt),wspd,rKa,sflux,DOsat
  real(rkind),pointer,dimension(:) :: DOX,temp,salt

  dtw=dt/86400.0
  do ie=1,nea
    !link variable
    temp=>tr_el(1,:,ie); salt=>tr_el(2,:,ie); DOX=>tr_el(3,:,ie)

    !compute geometry data
    zid=0; dz=0; srat=0; kb=min(kbe(ie),nvrt-1)
    do k=1,nvrt; zid(k)=ze(max(k,kb),ie); enddo             !zid : zcoor of each level
    do k=kb+1,nvrt; dz(k)=max(zid(k)-zid(k-1),1.d-2); enddo !dz : depth of each layer
    tdep=sum(dz((kb+1):nvrt))                               !tdep: total water depth 
    do k=kb+1,nvrt !surface flux ratio: y=2*(s-x)/s**2
        s=min(tdep,dz_flux); m=nvrt+kb+1-k
        z1=min(zid(nvrt)-zid(m),s); z2=min(zid(nvrt)-zid(m-1),s)
        srat(m)=min(max((z2-z1)*(2.0-(z1+z2)/s)/s,0.d0),1.d0)
    enddo

    !compute sflux of DO
    wspd=sum(sqrt(windx(elnode(1:i34(ie),ie))**2.d0+windy(elnode(1:i34(ie),ie))**2.d0))/dble(i34(ie))
    DOsat=14.5532-0.38217*temp(nvrt)+5.4258e-3*temp(nvrt)*temp(nvrt)-salt(nvrt)*(1.665e-4-5.866e-6*temp(nvrt)+9.796e-8*temp(nvrt)*temp(nvrt))/1.80655
    rKa=WRea+0.157*(0.54+0.0233*temp(nvrt)-0.002*salt(nvrt))*wspd**1.5
    sflux=rKa*(DOsat-DOX(nvrt))

    !update DO 
    do k=kb+1,nvrt
      DOX(k)=DOX(k)+dtw*(srat(k)*sflux/dz(k)-dr(ie)*(KT**(temp(k)-TR)))
      DOX(k)=max(DOX(k),0.d0)
    enddo !k

  enddo !nea

end subroutine GEN_model

subroutine GEN_init 
!------------------------------------------------------------
!initilize GEN model
!------------------------------------------------------------
  use GEN_mod
  use schism_glbl, only : rkind,nea,tr_el,in_dir,len_in_dir,np_global,ne_global, &
                        & elnode,ielg,iplg,i34
  use schism_msgp, only : myrank,parallel_abort,comm,itype,rtype
  use netcdf
  implicit none
  include 'mpif.h'

  !local variables
  integer :: i,j,ie,m,npt,istat,ncid,varid,attrid,dimid(1)
  real(rkind) :: swild(max(np_global,ne_global))
  character(len=20) :: fname
  integer, parameter :: nsp=5
  character(len=15) :: pname(nsp)
  type(gen_param),pointer :: p

  !init parameter
  allocate(dr(nea),stat=istat)
  if(istat/=0) call parallel_abort('Failed in alloc. dr')

  !read parameters
  fname='GEN_param.nc'
  allocate(sp(nsp),stat=istat)
  if(istat/=0) call parallel_abort('Failed in alloc. pname') 
  pname(1:nsp)=(/'dz_flux','dr','WRea','KT','TR'/)
  sp(1)%p=>dz_flux; sp(2)%p1=>dr; sp(3)%p=>WRea; sp(4)%p=>KT; sp(5)%p=>TR

  !read constant parameters
  do m=1,nsp
    p=>sp(m)
    !get dimension info
    if(associated(p%p)) then
      p%nv=1
    else
      p%nv=size(p%p1)
    endif
    p%varname=trim(adjustl(pname(m)))
    
    !read value on myrank=0   
    if(myrank==0) then
      j=nf90_open(in_dir(1:len_in_dir)//trim(adjustl(fname)),OR(NF90_NETCDF4,NF90_NOWRITE),ncid)
      if(j/=NF90_NOERR) call parallel_abort(trim(adjustl(fname))//': open')
      if(p%nv==1) then !single value, saved as global attribute
         j=nf90_get_att(ncid, NF90_GLOBAL,trim(adjustl(p%varname)),p%p)
         if(j/=NF90_NOERR) call parallel_abort(trim(adjustl(p%varname))//': wrong att value')
      else  !data array
          j=nf90_inq_varid(ncid,trim(adjustl(p%varname)),varid)
          if(j/=NF90_NOERR) call parallel_abort(trim(adjustl(p%varname))//': wrong varid' )
          j=nf90_inquire_variable(ncid,varid,dimids=dimid(1:1))
          if(j/=NF90_NOERR) call parallel_abort(trim(adjustl(p%varname))//': wrong dimid')
          j=nf90_inquire_dimension(ncid,dimid(1),len=npt)
          if(j/=NF90_NOERR) call parallel_abort(trim(adjustl(p%varname))//': wrong npt')
          if(npt/=np_global.and.npt/=ne_global) call parallel_abort(trim(adjustl(p%varname))//': npt/=ne,np' ) 
          j=nf90_get_var(ncid,varid,swild(1:npt), (/1/),(/npt/))
      endif
    endif !myrank

    !bcast values
    if(p%nv==1) then
      call mpi_bcast(p%p,1,rtype,0,comm,istat)
    else
      call mpi_bcast(npt,1,itype,0,comm,istat)
      call mpi_bcast(swild,max(np_global,ne_global),rtype,0,comm,istat)
    endif

    !get value for each element
    if(p%nv/=1) then
      do ie=1,nea
        p%p1(ie)=0.d0
        if(npt==ne_global) then 
          p%p1(ie)=swild(ielg(ie))
        else
          do i=1,i34(ie); p%p1(ie)=p%p1(ie)+swild(iplg(elnode(i,ie)))/dble(i34(ie)); enddo 
        endif
      enddo
    endif

  enddo !m
end subroutine GEN_init 

