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

!  Adapted from FESOM's ice module. Special thanks to Dr. Sergey Danilov's group
!  for their generous help.
!  Standing alone sea ice Version 2, based on version 1, with several new features added
!  Most important are true VP solver and FCT advection  
!  Questions to S. Danilov (dynamics) and Q. Wang and R. Timmermann 
! (thermodynamics). Many updates and corrections
!  to version 1 are due to R. Timmermann
! ======================
!==============================================================================

subroutine ice_step
  use schism_glbl, only: rkind,pi,np,npa,nvrt,uu2,vv2,time_stamp,windx,windy,xnd,ynd, &
 &tau_oi,nws,ihconsv,isconsv,iplg,fresh_wa_flux,net_heat_flux,rho0,rnday,fdb,lfdb, &
 &lice_free_gb,lhas_ice
  use schism_msgp, only : myrank,nproc,parallel_abort,comm,ierr
  use ice_module
  use ice_therm_mod

  implicit none
  include 'mpif.h'

!  integer, intent(in) :: npa1 !for dim only
!  real(rkind), intent(inout) :: tau_oi(2,npa1) !ocean-ice stress (junk if no ice)

  integer :: i
  real(rkind) :: tmp1,uwind,vwind,umod
  logical :: lice_free

  !Set wind and ocean vel
  do i=1,npa
    if(ice_tests==0) then !not test
      u_ocean(i)=uu2(nvrt,i); v_ocean(i)=vv2(nvrt,i)
      uwind=windx(i); vwind=windy(i)
    else !box test
      u_ocean(i)=0.1*(2*(ynd(i)-ymin_ice)-rly_ice)/rly_ice
      v_ocean(i)=-0.1*(2*(xnd(i)-xmin_ice)-rlx_ice)/rlx_ice
      uwind=5+(sin(2*pi*time_stamp/4/86400)-3)*sin(2*pi*(xnd(i)-xmin_ice)/rlx_ice)*sin(pi*(ynd(i)-ymin_ice)/rly_ice) 
      vwind=5+(sin(2*pi*time_stamp/4/86400)-3)*sin(2*pi*(ynd(i)-ymin_ice)/rly_ice)*sin(pi*(xnd(i)-xmin_ice)/rlx_ice)
    endif !ice_tests

    tmp1=rhoair*cdwin*sqrt(uwind**2+vwind**2)
    stress_atm_ice(1,i)=tmp1*uwind !Pa
    stress_atm_ice(2,i)=tmp1*vwind

    !Debug
!    if(abs(time_stamp-21600)<1.e-2) then
!      write(98,*)real(xnd(i)),real(ynd(i)),real(u_ocean(i)),real(v_ocean(i))
!      write(99,*)real(xnd(i)),real(ynd(i)),real(uwind),real(vwind)
!    endif
  enddo !i

  if(.not.lice_free_gb) then
    !EVP dynamics
    if(ievp==1) then
      call ice_evp
      if(myrank==0) write(16,*)'done ice EVP dynamics'
    else if (ievp==2) then
      call ice_mevp
      if(myrank==0) write(16,*)'done ice mEVP dynamics'
    endif !ievp

    !Transport: operator splitting
    if(ice_advection/=0) then
      call ice_fct
      if(myrank==0) write(16,*)'done ice FCT advection'
    endif
  endif !not ice free

  if(ice_therm_on==1) then
    if(ice_tests==0.and.(nws/=2.or.ihconsv/=1.or.isconsv/=1)) &
  &call parallel_abort('ice_step: ice therm needs nws=2 etc')
    !Atmos variables are read in for thermodynamics

    call ice_thermodynamics
    if(myrank==0) write(16,*)'done ice thermodynamics'
  endif

  !Mark ice/no ice
  lice_free=.true. !over all aug domain
  do i=1,npa
    if(ice_tr(1,i)<=ice_cutoff.or.ice_tr(2,i)<=ice_cutoff) then
      lhas_ice(i)=.false.
    else
      lhas_ice(i)=.true.
      lice_free=.false.
    endif
  enddo !i
  call mpi_allreduce(lice_free,lice_free_gb,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
  if(myrank==0) write(16,*)'lice_free_gb=',lice_free_gb

  !Update tau_oi (for ocean)
  if(ice_tests==0) then !not test
    tau_oi=-huge(1.d0) !init as junk
    do i=1,npa
      if(lhas_ice(i)) then
        umod=sqrt((u_ice(i)-u_ocean(i))**2+(v_ice(i)-v_ocean(i))**2)
        tmp1=ice_tr(2,i)*cdwat*umod
        tau_oi(1,i)=tmp1*((u_ice(i)-u_ocean(i))*cos_io-(v_ice(i)-v_ocean(i))*sin_io) !m^2/s/s
        tau_oi(2,i)=tmp1*((u_ice(i)-u_ocean(i))*sin_io+(v_ice(i)-v_ocean(i))*cos_io)
      endif
    enddo !i
  else !test
    tau_oi=0
  endif !ice_tests

  !Debug
!  if(abs(time_stamp-rnday*86400)<0.1) then
!    fdb='icefw_0000'
!    lfdb=len_trim(fdb)
!    write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
!    open(10,file='outputs/'//fdb,status='replace')
!    write(10,*)np,nproc
!    do i=1,np
!      write(10,'(i11,3(1x,e20.12))')iplg(i),xnd(i),ynd(i),fresh_wa_flux(i),net_heat_flux(i),t_oi(i)
!    enddo !i
!    close(10)
!  endif

end subroutine ice_step
