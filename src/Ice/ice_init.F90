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
!====================================================================
! Init ice vars
subroutine ice_init
  use schism_glbl, only : rkind,pi,np,npa,ne,nea,mnei,mnei_p,nne,indel,xctr,yctr,area, &
 &nstep_ice,fresh_wa_flux,net_heat_flux,xlon,ylat,rearth_eq,elnode,nnp,indnd,iplg,dt, &
 &xnd,ynd,errmsg,lice_free_gb,in_dir,out_dir,len_in_dir,len_out_dir
  use schism_msgp, only : myrank,parallel_abort,parallel_finalize,exchange_p2d
  use ice_module
  use ice_therm_mod
  implicit none
  integer :: i,j,ie,istat,nd,nd2,m,mm,indx
  real(rkind) :: sum1,meancos,local_cart(2,3),jacobian2D(2,2),jacobian2D_inv(2,2), &
 &det,der_transp(3,2),derivative_stdbf(2,3) 
  namelist /ice_in/ice_tests,ice_advection,ice_therm_on,ievp,ice_cutoff,evp_rheol_steps,mevp_rheol_steps, &
 &delta_min,theta_io,mevp_alpha1,mevp_alpha2,pstar,ellipse,c_pressure,niter_fct, &
 &ice_gamma_fct,h_ml0,salt_ice,salt_water
  
  !Init parameters
  !integers
  ice_tests=-1e6; ice_advection=-1e6; ice_therm_on=-1e6; ievp=-1e6; evp_rheol_steps=-1e6;
  mevp_rheol_steps=-1e6; niter_fct=-1e6; 
  !Doubles
  ice_cutoff=-huge(1.d0); delta_min=-huge(1.d0); theta_io=-huge(1.d0);
  mevp_alpha1=-huge(1.d0); mevp_alpha2=-huge(1.d0); pstar=-huge(1.d0);
  ellipse=-huge(1.d0); c_pressure=-huge(1.d0); ice_gamma_fct=-huge(1.d0);
  h_ml0=-huge(1.d0); salt_ice=-huge(1.d0); salt_water=-huge(1.d0)

  open(10,file=in_dir(1:len_in_dir)//'ice.nml',status='old')
  read(10,nml=ice_in)
  close(10)

  !Check
  if(ice_tests/=0.and.ice_tests/=1) call parallel_abort('ice_init: ice_tests')
  if(ice_advection/=0.and.ice_advection/=1) call parallel_abort('ice_init: ice_advection')
  if(ice_therm_on/=0.and.ice_therm_on/=1) call parallel_abort('ice_init: ice_therm_on')
  if(ievp/=1.and.ievp/=2) call parallel_abort('ice_init: ievp')
  if(ice_cutoff<=0) call parallel_abort('ice_init: ice_cutoff')
  if(evp_rheol_steps<=0.or.mevp_rheol_steps<=0) call parallel_abort('ice_init: evp_rheol_steps')
  if(delta_min<=0) call parallel_abort('ice_init: delta_min')
  if(abs(theta_io)>360) call parallel_abort('ice_init: theta_io')
  if(mevp_alpha1<=0.or.mevp_alpha2<=0) call parallel_abort('ice_init: mevp_alpha1')
  if(pstar<=0) call parallel_abort('ice_init: pstar')
  if(ellipse<=0.or.c_pressure<=0) call parallel_abort('ice_init: ellipse')
  if(niter_fct<=0) call parallel_abort('ice_init: niter_fct')
  if(ice_gamma_fct<0) call parallel_abort('ice_init: ice_gamma_fct')
  if(h_ml0<=0) call parallel_abort('ice_init: h_ml0')
  if(salt_ice<0.or.salt_water<0) call parallel_abort('ice_init: salt_water')
  
  dt_ice=dt*nstep_ice
  cos_io=cos(theta_io/180*pi)
  sin_io=sin(theta_io/180*pi)

  allocate(u_ice(npa),v_ice(npa),ice_tr(ntr_ice,npa),stress_atm_ice(2,npa), &
     &sigma11(nea),sigma12(nea),sigma22(nea),weit_elem2node(mnei,np),u_ocean(npa), &
     &v_ocean(npa),area_median(np),voltriangle(nea),bafux(3,nea),bafuy(3,nea), &
     &ice_matrix(0:mnei_p,np),lump_ice_matrix(npa),delta_ice(nea),t_oi(npa),stat=istat)
  if(istat/=0) call parallel_abort('ice_init: alloc (1)')
!  if(ice_therm_on==1) then
!    allocate(t_oi(npa),stat=istat)
!    if(istat/=0) call parallel_abort('ice_init: alloc (2)')
!  endif

  t_oi=0 !init T @snow/ice surface in C
  u_ice=0; v_ice=0; sigma11=0; sigma12=0; sigma22=0
  fresh_wa_flux=0; net_heat_flux=0

  !Box test
  if(ice_tests==0) then !normal
    ice_tr=0
    lice_free_gb=.true. !for normal cases, always start from ice free
  else !box test
    xmin_ice=-33.34015; ymin_ice=3334060.
    xmax_ice=1222472.; ymax_ice=4556532.
    rlx_ice=xmax_ice-xmin_ice; rly_ice=ymax_ice-ymin_ice
    ice_tr(1,:)=2
    do i=1,npa
      ice_tr(2,i)=(xnd(i)-xmin_ice)/(xmax_ice-xmin_ice)
      ice_tr(2,i)=max(0.d0,min(1.d0,ice_tr(2,i)))
    enddo !i
    !'Net' ice volume [m]
    ice_tr(1,:)=ice_tr(1,:)*ice_tr(2,:)
    ice_tr(3,:)=0
    lice_free_gb=.false.
  endif !ice_tests

  !Calc weights for interpolating from elem to node (via the ball)
  !and area of median dual grid (sum of integrals)
  !Node must not be ghost
  !The ball averaging is done as:
  !dot_product(weit_elem2node(1:nne(i),i),h_ice(indel(1:nne(i),i)))
  !where i=1,np and h_ice(1:nea) is defined @ elem
  weit_elem2node=0
  area_median=0
  do i=1,np
    sum1=0
    do j=1,nne(i)
      ie=indel(j,i)
      sum1=sum1+1/area(ie)
!Error: tri
      area_median(i)=area_median(i)+area(ie)/3
    enddo !j
    do j=1,nne(i)
      ie=indel(j,i)
      weit_elem2node(j,i)=1/area(ie)/sum1
    enddo !j
  enddo !i

! Prep for FCT transport
  !Derivatives of shape function
  !stdbafu(1)=1-x-y, stdbafu(2)=x, stdbafu(3)=y
  derivative_stdbf=0.
  derivative_stdbf(:,1)=-1.
  derivative_stdbf(1,2)=1.
  derivative_stdbf(2,3)=1.

  bafux=0.
  bafuy=0.
  voltriangle=0.
  do i=1,nea
    meancos=sum(cos(ylat(elnode(1:3,i))))/3.d0

    do j=1,3 !nodes
      nd=elnode(j,i)
      local_cart(1,j)=xlon(nd) !radian
      local_cart(2,j)=ylat(nd)*rearth_eq
    enddo !j

    !Jacobian of transform
    do j=1,2 !nodes 2&3
      jacobian2D(:,j)=local_cart(:,j+1)-local_cart(:,1)
      !make sure \in [-pi,pi]
      !Error: distance not rite across dateline
      if(jacobian2D(1,j)>pi) jacobian2D(1,j)=jacobian2D(1,j)-2*pi
      if(jacobian2D(1,j)<-pi) jacobian2D(1,j)=jacobian2D(1,j)+2*pi
    enddo !j
    jacobian2D(1,:)=jacobian2D(1,:)*meancos*rearth_eq
    det=jacobian2D(1,1)*jacobian2D(2,2)-jacobian2D(1,2)*jacobian2D(2,1)
    if(det==0) then
      write(errmsg,*)'ice_init: det=0,',jacobian2D
      call parallel_abort(errmsg)     
    else
      jacobian2D_inv(1,1)=jacobian2D(2,2)/det
      jacobian2D_inv(1,2)=-jacobian2D(1,2)/det
      jacobian2D_inv(2,1)=-jacobian2D(2,1)/det
      jacobian2D_inv(2,2)=jacobian2D(1,1)/det
    endif !det

    der_transp=matmul(transpose(derivative_stdbf),jacobian2D_inv)
    !derivative_loczeta=transpose(der_transp)

    do j=1,3
      bafux(j,i)=der_transp(j,1)
      bafuy(j,i)=der_transp(j,2)
    enddo
    voltriangle(i)=abs(det)*0.5d0

    !Debug
!    write(12,*)'voltriangle:',ielg(i),voltriangle(i)
!    write(12,*)'bafux:',ielg(i),bafux(:,i)
!    write(12,*)'bafuy:',ielg(i),bafuy(:,i)
  enddo !i=1,nea

  !Mass matrix: ice_matrix(0:mnei_p,np)
  ice_matrix=0.0
  do i=1,ne
    do j=1,3
      nd=elnode(j,i)
      do m=1,3
        nd2=elnode(m,i)
        if(j==m) then !diagonal
          ice_matrix(0,nd)=ice_matrix(0,nd)+voltriangle(i)/6.0
        else
          indx=0
          do mm=1,nnp(nd)
            if(indnd(mm,nd)==nd2) then
              indx=mm; exit
            endif
          enddo !mm
          if(indx==0) call parallel_abort('ice_init: failed (9.1)')
          ice_matrix(indx,nd)=ice_matrix(indx,nd)+voltriangle(i)/12.0
        endif
      enddo !m
    enddo !j
  enddo !i=1,ne

  do i=1,np
    lump_ice_matrix(i)=sum(ice_matrix(0:nnp(i),i))
    if(lump_ice_matrix(i)==0) then
      write(errmsg,*)'ice_init: lump=0,',iplg(i)
      call parallel_abort(errmsg)
    endif

    !Debug
!    write(12,*)'mass:',iplg(i),lump_ice_matrix(i) !,area_median(i)
  enddo !i=1,np

  call exchange_p2d(lump_ice_matrix)

!  call parallel_finalize
!  stop

end subroutine ice_init
