! implicit vp solver of ice module
#ifdef USE_PETSC
subroutine ice_vp_imp(u_ice0,v_ice0)
  use schism_glbl,only: rkind,time_stamp,rnday,eta2,np,ne,nea, &
 &elnode,i34,dldxy,cori,grav,isbnd,indel,nne,area,iself,fdb,lfdb, &
 &xnd,ynd,iplg,ielg,elside,mnei,rho0,idry,errmsg,indnd,nnp,mnei_p,npa,&
 &lelbc
  use schism_msgp, only: myrank,nproc,parallel_abort,exchange_p2d,ierr,comm
  use ice_module
  use petsc_schism
  implicit none

  real(rkind), dimension(npa), intent(in) :: u_ice0,v_ice0

  integer :: isub,i,j,k,ie,n1,n2,icount,kk,ll,m,n,id,n_columns,&
            &petsc_its_u,petsc_its_v,mxitn0,nd0,istat

  real(rkind) :: dtevp,t_evp_inv,det1,det2,tmp1,tmp2,sum1,sum2, &
 &pp0,delta,delta_inv,rr1,rr2,rr3,sig1,sig2,x10,x20,y10,y20,rl10, &
 &rl20,sintheta,bb1,bb2,h_ice_el,a_ice_el,h_snow_el,dsig_1,dsig_2,mass, &
 &cori_nd,umod,gam1,rx,ry,elevx,elevy,eps11,eps12,eps22, &
 &zeta,eta,zeta_over_T,inv_m,vale,pressure

  integer :: iball(mnei)

  real(rkind) :: swild(2,3),swild2(nea),deta(2,nea),sparsem(0:mnei_p,np),rhsu_vp_local(np),&
  &rhsv_vp_local(np),rhs_a(npa),rhs_m(npa),entries(3)

  integer :: column_ix_vp(0:mnei_p)
  real(rkind) :: coeff_vals_vp(0:mnei_p),uice_npi(npi),vice_npi(npi),rhsu_vp(npi),rhsv_vp(npi),&
    &uice_guess_npi(npi),vice_guess_npi(npi)

  mxitn0=1500*2

!Fills rhs into rhsu_vp and rhsv_vp
!Fills the vp rheology stiffness matrix into coeff_vals_vp and solve it by PETSC
!Refer to Danilov et al., GMD, 2014

!Refer to eta2 in schism_hydro.F90
  do i=1,nea
    swild2(i)=sum(cori(elside(1:i34(i),i)))/i34(i)
  enddo

  do i=1,np !resident only
    do j=0,nnp(i)
      sparsem(j,i)=0.d0
    enddo

    rhsu_vp_local(i)=0.d0
    rhsv_vp_local(i)=0.d0

    mass=rhoice*ice_tr(1,i)+rhosnow*ice_tr(3,i)
    inv_m=max(mass,0.1_8)
    inv_m=ice_tr(2,i)/inv_m

!ice-ocean drag
    umod=sqrt((u_ice(i)-u_ocean(i))**2+(v_ice(i)-v_ocean(i))**2)
    gam1=inv_m*cdwat(i)*rho0*umod

    sparsem(0,i)=lump_ice_matrix(i)*(1.0/dt_ice+gam1*cos_io) !diagonal 

    !Coriolis @ node
    iball(1:nne(i))=indel(1:nne(i),i)
    cori_nd=dot_product(weit_elem2node(1:nne(i),i),swild2(iball(1:nne(i))))
    !RHS: mass matrix part+ocean/ice friction part + the Coriolis part.
    rhsu_vp_local(i) = lump_ice_matrix(i)*(u_ice0(i)/dt_ice+&
                       &gam1*cos_io*u_ocean(i)+cori_nd*v_ice(i)+&
                       stress_atm_ice(1,i)*inv_m-&
                       &sin_io*gam1*(v_ocean(i)-v_ice(i)))

    rhsv_vp_local(i) = lump_ice_matrix(i)*(v_ice0(i)/dt_ice+&
                       &gam1*cos_io*v_ocean(i)-cori_nd*u_ice(i)+&
                       stress_atm_ice(2,i)*inv_m+&
                       &sin_io*gam1*(u_ocean(i)-u_ice(i)))

    if(idry(i)/=0 .or. ice_tr(2,i)<0.01_rkind) then ! skip if dry or very thin ice
    ! thin ice velocity is expected to be zero, so we can skip the rheology part and set the rhs to zero.
    !if(idry(i)/=0) then ! skip if dry or very thin ice
      rhsu_vp_local(i)=0.d0
      rhsv_vp_local(i)=0.d0
    endif
  enddo

  !ice rheology part
  do i=1,ne

    delta=product(ice_tr(1,elnode(1:3,i)))*product(ice_tr(2,elnode(1:3,i)))

    if(maxval(idry(elnode(1:i34(i),i)))/=0) cycle ! skip if dry
    if(delta==0) cycle ! skip if no ice

    eps11=dot_product(u_ice(elnode(1:i34(i),i)),dldxy(1:i34(i),1,i)) !du_dx
    eps22=dot_product(v_ice(elnode(1:i34(i),i)),dldxy(1:i34(i),2,i)) !dv_dy
    tmp1=dot_product(u_ice(elnode(1:i34(i),i)),dldxy(1:i34(i),2,i)) !du_dy
    tmp2=dot_product(v_ice(elnode(1:i34(i),i)),dldxy(1:i34(i),1,i)) !dv_dx

    eps12=0.5*(tmp1+tmp2)

    vale=1.0_8/(ellipse**2)

    delta=(eps11**2+eps22**2)*(1.0_8+vale)+4.0_8*vale*eps12**2 + &
            2.0_8*eps11*eps22*(1.0_8-vale)

    delta=sqrt(delta)
    delta_ice(i)=delta

    h_ice_el=sum(ice_tr(1,elnode(1:3,i)))/3.0
    a_ice_el=sum(ice_tr(2,elnode(1:3,i)))/3.0
    h_snow_el=sum(ice_tr(3,elnode(1:3,i)))/3.0_rkind

    mass = rhoice*h_ice_el + rhosnow*h_snow_el

    pp0 = h_ice_el*pstar*exp(-c_pressure*(1-a_ice_el)) !P_0

    zeta=pp0*0.5/(delta+delta_min)

    pressure=zeta*delta
    eta=zeta*vale

    if(maxval(idry(elnode(1:i34(i),i)))/=0) then
      elevx = 0.d0
      elevy = 0.d0
    else
      elevx=-1*dot_product(eta2(elnode(1:i34(i),i)),dldxy(1:i34(i),1,i))*area(i)/3*grav
      elevy=-1*dot_product(eta2(elnode(1:i34(i),i)),dldxy(1:i34(i),2,i))*area(i)/3*grav
    endif

    do j=1,i34(i)
      k=elnode(j,i)

      if(ice_tr(2,k)<0.01_rkind) cycle ! skip if very thin ice，so there is no ice rheology effect
      if(isbnd(1,k)/=0) cycle ! skip if boundary node
      rhsu_vp_local(k) = rhsu_vp_local(k) + elevx
      rhsv_vp_local(k) = rhsv_vp_local(k) + elevy

      mass=rhoice*ice_tr(1,k)+rhosnow*ice_tr(3,k)
      !mass = rhoice*h_ice_el + rhosnow*h_snow_el

      inv_m=max(mass,0.1_8)
      inv_m=1.0/inv_m

      entries = 0.d0

      do m=1,i34(i)
        n=elnode(m,i)

        if(ice_tr(2,n)<0.01_rkind) cycle ! skip if very thin ice，so there is no ice rheology effect
        if(isbnd(1,n)/=0) cycle ! skip if boundary node

        entries(m)=(eta+zeta)*(dldxy(j,1,i)*dldxy(m,1,i)+dldxy(j,2,i)*dldxy(m,2,i))* &
                  voltriangle(i)*inv_m

        if(n==k) then
          sparsem(0,k)=sparsem(0,k)+entries(m)
        else
          id=0
          do kk=1,nnp(k)
            ll=indnd(kk,k)
            if(ll==n) then
              id=kk
              exit
            endif
          enddo
          if(id==0) call parallel_abort('STEP: icemap(0)')
          sparsem(id,k)=sparsem(id,k)+entries(m)
        endif
      enddo

      rhsu_vp_local(k) = rhsu_vp_local(k) + &
       ( -eta*(dldxy(j,1,i)*eps11+dldxy(j,2,i)*tmp2) &
       -(zeta-eta)*dldxy(j,1,i)*(eps11+eps22) &
       + dldxy(j,1,i)*pressure ) * voltriangle(i) * inv_m &
       + zeta/(eta+zeta)*dot_product(entries,u_ice(elnode(1:i34(i),i)))

      rhsv_vp_local(k) = rhsv_vp_local(k) + &
       ( -eta*(dldxy(j,1,i)*tmp1+dldxy(j,2,i)*eps22) &
       -(zeta-eta)*dldxy(j,2,i)*(eps11+eps22) &
       + dldxy(j,2,i)*pressure ) * voltriangle(i) * inv_m &
       + zeta/(eta+zeta)*dot_product(entries,v_ice(elnode(1:i34(i),i)))

    enddo
  enddo

  do i=1,np
    ! Skip interface nodes that do not belong to this rank
    if(npa2npi(i)==-999) cycle

    ! Apply essential BC - 0 out rows and columns.
    n_columns=1
    nd0=npa2npi(i)
    if(nd0<=0) call parallel_abort('STEP: icemap(1)')

    column_ix_vp(0)=npa2npia(i)-1
    coeff_vals_vp=0.d0
    coeff_vals_vp(0)=sparsem(0,i)

    if(isbnd(1,i)/=0) then
      ! skip if boundary node, set diagonal to 1 and rhs to 0, so velocity is fixed to 0.
      coeff_vals_vp(0)=1.d0
      rhsu_vp(nd0)=0.d0
      rhsv_vp(nd0)=0.d0
    else
      rhsu_vp(nd0)=rhsu_vp_local(i)
      rhsv_vp(nd0)=rhsv_vp_local(i)

      do j=1,nnp(i)
        k=indnd(j,i)
        if(isbnd(1,k)/=0) then ! no change to rhs if boundary node, since no external force from boundary node
          rhsu_vp(nd0) = rhsu_vp(nd0)
          rhsv_vp(nd0) = rhsv_vp(nd0)
        else
          n_columns=n_columns+1
          if(npa2npia(k)<=0) call parallel_abort('STEP: icemap(2)')
          column_ix_vp(n_columns-1)=npa2npia(k)-1
          coeff_vals_vp(n_columns-1)=sparsem(j,i)
        endif
      enddo
    endif

    call load_mat_row_ice(icevp_M,row_ix=npa2npia(i)-1,n_columns=n_columns, &
                          column_ix=column_ix_vp(0:n_columns-1), &
                          coeff_vals=coeff_vals_vp(0:n_columns-1))
  enddo

  if(myrank==0) write(16,*)'call ice VP petsc solver...'
  !Use the ice velocity from the previous time step or previous nonlinear iteration as the initial guess,
  !instead of starting from a zero velocity field.
  do i=1,npi
    uice_guess_npi(i)=u_ice(npi2np(i))
    vice_guess_npi(i)=v_ice(npi2np(i))
  enddo

  call petsc_solve_VP(npi,rhsu_vp,rhsv_vp,uice_guess_npi,vice_guess_npi,&
  &uice_npi,vice_npi,petsc_its_u,petsc_its_v)


  if(myrank==0) then
    write(33,'(//a,i8)') '********ice VP PetSc Solve'
    if((petsc_its_u+petsc_its_v)>=0.and.(petsc_its_u+petsc_its_v)<mxitn0) then
      write(33,*)'converged in ',petsc_its_u,petsc_its_v
    else
      write(33,*)'diverged:',petsc_its_u,petsc_its_v,mxitn0
    endif
  endif

  do i=1,npi
    !if(ice_tr(2,npi2np(i))<0.01_rkind) then
    !  u_ice(npi2np(i)) = u_ocean(npi2np(i))
    !  v_ice(npi2np(i)) = v_ocean(npi2np(i))
    !else
      u_ice(npi2np(i))=max(-5.d0,min(5.d0,uice_npi(i)))
      v_ice(npi2np(i))=max(-5.d0,min(5.d0,vice_npi(i)))
    !endif
  enddo

  if(myrank==0) write(16,*)'done ice VP petsc...'

  call exchange_p2d(u_ice)
  call exchange_p2d(v_ice)

end subroutine ice_vp_imp

subroutine ice_vp_imp_step
  use schism_glbl,only: rkind,npa
  use schism_msgp, only: myrank,exchange_p2d
  use ice_module

  implicit none

  integer :: isub

  real(rkind) :: u_ice0(npa),v_ice0(npa)

  u_ice0 = u_ice
  v_ice0 = v_ice

  if(myrank==0) write(16,*)'start ice VP petsc...'

  call ice_vp_imp(u_ice0,v_ice0)

  do isub=1,ice_VP_iter

    u_ice=0.5_rkind*(u_ice+u_ice0)
    v_ice=0.5_rkind*(v_ice+v_ice0)

    call exchange_p2d(u_ice)
    call exchange_p2d(v_ice)

    call ice_vp_imp(u_ice0,v_ice0)

  enddo

end subroutine ice_vp_imp_step
#endif /*USE_PETSC*/
