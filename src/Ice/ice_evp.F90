! EVP dynamics of ice module
subroutine ice_evp
  use schism_glbl,only: rkind,time_stamp,rnday,eta2,np,ne,nea, &
 &elnode,i34,dldxy,cori,grav,isbnd,indel,nne,area,iself,fdb,lfdb, &
 &xnd,ynd,iplg,ielg,elside,mnei,rho0,idry,errmsg
  use schism_msgp, only: myrank,nproc,parallel_abort,exchange_p2d
  use ice_module
  implicit none
  integer :: isub,i,j,ie,n1,n2,icount,m,id
  real(rkind) :: dtevp,t_evp_inv,det1,det2,tmp1,tmp2,sum1,sum2, &
 &pp0,delta,delta_inv,rr1,rr2,rr3,sig1,sig2,x10,x20,y10,y20,rl10, &
 &rl20,sintheta,bb1,bb2,h_ice_el,a_ice_el,h_snow_el,dsig_1,dsig_2,mass, &
 &cori_nd,umod,gam1,rx,ry,rxp,ryp,eps11,eps12,eps22, &
 &zeta,delta_nd

  integer :: iball(mnei)
  real(rkind) :: swild(2,3),swild2(nea),deta(2,nea)

  dtevp=dt_ice/evp_rheol_steps
  t_evp_inv=3/dt_ice !inverse of T
  det1=1./(1+0.5*dtevp*t_evp_inv) !used in solving EVP eqs
  det2=1./(1+0.5*ellipse*ellipse*dtevp*t_evp_inv) 

  do isub=1,evp_rheol_steps
    !Update stress @ elem
    do i=1,nea
!      if(h_ice(i)<=ice_cutoff.or.a_ice(i)<=ice_cutoff) then
!        sigma11(i)=0; sigma12(i)=0; sigma22(i)=0; delta_ice(i)=0
!        cycle
!      endif

      eps11=dot_product(u_ice(elnode(1:i34(i),i)),dldxy(1:i34(i),1,i)) !epsilon11=du_dx
      eps22=dot_product(v_ice(elnode(1:i34(i),i)),dldxy(1:i34(i),2,i)) !epsilon22=dv_dy
      tmp1=dot_product(u_ice(elnode(1:i34(i),i)),dldxy(1:i34(i),2,i)) !du_dy
      tmp2=dot_product(v_ice(elnode(1:i34(i),i)),dldxy(1:i34(i),1,i)) !dv_dx
      eps12=0.5*(tmp1+tmp2) !epsilon12
      !Deformation rate
      tmp1=eps11+eps22 !divergence
      tmp2=(eps11-eps22)**2+4*eps12*eps12 !shear strain rate squared
      delta_ice(i)=sqrt(tmp1*tmp1+tmp2/ellipse/ellipse)
      h_ice_el=sum(ice_tr(1,elnode(1:3,i)))/3.0
      a_ice_el=sum(ice_tr(2,elnode(1:3,i)))/3.0
      pp0=h_ice_el*pstar*exp(-c_pressure*(1-a_ice_el)) !P_0
      zeta=pp0*0.5/max(delta_ice(i),delta_min)

      rr1=zeta*t_evp_inv*(eps11+eps22-delta_ice(i)) !RHS for 1st eq
      rr2=zeta*t_evp_inv*(eps11-eps22)
      rr3=zeta*t_evp_inv*eps12
      sig1=sigma11(i)+sigma22(i) !from previous step
      sig2=sigma11(i)-sigma22(i)

      sig1=det1*(sig1+rr1*dtevp)
      sig2=det2*(sig2+rr2*dtevp)

      sigma12(i)=det2*(sigma12(i)+rr3*dtevp)
      sigma11(i)=0.5*(sig1+sig2)
      sigma22(i)=0.5*(sig1-sig2)
    enddo !i=1,nea

!    if(isub==evp_rheol_steps.and.it_main==1) then
!      do i=1,nea
!        write(96,*)i,real(sigma11(i))
!        write(95,*)i,real(sigma12(i))
!        write(94,*)i,real(sigma22(i))
!      enddo !i
!      close(94); close(95); close(96);
!    endif

    !RHS of mom eq
    !Derivatives of stress: projection/reconstruction method
!    dsigdxy=0 !init
!    do i=1,ne
!      dsigdxy(i,1,1)=dot_product(sigma11(elnode(1:i34(i),i)),dldxy(1:i34(i),1,i)) !d{sigma11}/dx
!      dsigdxy(i,2,3)=dot_product(sigma22(elnode(1:i34(i),i)),dldxy(1:i34(i),2,i)) !d{sigma22}/dy
!      dsigdxy(i,1,2)=dot_product(sigma12(elnode(1:i34(i),i)),dldxy(1:i34(i),1,i)) !d{sigma12}/dx
!      dsigdxy(i,2,2)=dot_product(sigma12(elnode(1:i34(i),i)),dldxy(1:i34(i),2,i)) !d{sigma12}/dy
!    enddo !i
!
!    do i=1,3
!      do j=1,2
!        swild2=dsigdxy(:,j,i)
!        call exchange_e2d(swild2)
!        dsigdxy(:,j,i)=swild2
!      enddo !j
!    enddo !i

    !Coriolis @ elem
    do i=1,nea
      swild2(i)=sum(cori(elside(1:i34(i),i)))/i34(i)
      !Pressure gradient @elem
      if(ice_tests==1) then
        deta(1,i)=swild2(i)*sum(v_ocean(elnode(1:i34(i),i)))/i34(i)/grav
        deta(2,i)=-swild2(i)*sum(u_ocean(elnode(1:i34(i),i)))/i34(i)/grav
      else
        if(maxval(idry(elnode(1:i34(i),i)))/=0) then !dry
          deta(1:2,i)=0
        else !wet
          deta(1,i)=dot_product(eta2(elnode(1:i34(i),i)),dldxy(1:i34(i),1,i))
          deta(2,i)=dot_product(eta2(elnode(1:i34(i),i)),dldxy(1:i34(i),2,i))
        endif !dry
      endif !
    enddo !i

    !Solve mom eq.
    do i=1,np !resident
      if(isbnd(1,i)/=0) then !b.c. (including open)
        u_ice(i)=0; v_ice(i)=0
        cycle
      endif

      iball(1:nne(i))=indel(1:nne(i),i)
      if(ice_tr(1,i)<=ice_cutoff.or.ice_tr(2,i)<=ice_cutoff) then !no ice
        u_ice(i)=0; v_ice(i)=0
        cycle 
      endif
   
      !Not bnd node; has ice
      mass=rhoice*ice_tr(1,i)+rhosnow*ice_tr(3,i) !>0
      !Error: limit mass>9?
      !Coriolis @ node
      cori_nd=dot_product(weit_elem2node(1:nne(i),i),swild2(iball(1:nne(i))))
      !Debug
      !if(isub==1.and.it_main==1) then
      !  write(93,*)i,real(xnd(i)),real(ynd(i)),real(),real(cori_nd) !,real(h_ice_nd),real(mass)
      !endif

      !RHS
!Error: tri only
      umod=sqrt((u_ice(i)-u_ocean(i))**2+(v_ice(i)-v_ocean(i))**2)
      gam1=ice_tr(2,i)/mass*cdwat*rho0*umod
      rxp=ice_tr(2,i)*stress_atm_ice(1,i)/mass
      ryp=ice_tr(2,i)*stress_atm_ice(2,i)/mass
      !Pressure gradient
      sum1=0; sum2=0
      do j=1,nne(i)
        ie=indel(j,i)
        h_ice_el=sum(ice_tr(1,elnode(1:3,ie)))/3.0
        h_snow_el=sum(ice_tr(3,elnode(1:3,ie)))/3.0
        tmp1=rhoice*h_ice_el+rhosnow*h_snow_el !mass @elem
        sum1=sum1+tmp1*deta(1,ie)*area(ie)/3
        sum2=sum2+tmp1*deta(2,ie)*area(ie)/3
      enddo !j
      rxp=rxp-grav*sum1/mass/area_median(i)
      ryp=ryp-grav*sum2/mass/area_median(i)
      !Stress
      sum1=0; sum2=0
      do j=1,nne(i)
        ie=indel(j,i)
        id=iself(j,i)
        sum1=sum1+area(ie)*(dldxy(id,1,ie)*sigma11(ie)+dldxy(id,2,ie)*sigma12(ie))
        sum2=sum2+area(ie)*(dldxy(id,1,ie)*sigma12(ie)+dldxy(id,2,ie)*sigma22(ie))
      enddo !j
      rxp=rxp-sum1/mass/area_median(i)
      ryp=ryp-sum2/mass/area_median(i)

      rx=u_ice(i)+rxp*dtevp+gam1*dtevp*(u_ocean(i)*cos_io-v_ocean(i)*sin_io)
      ry=v_ice(i)+ryp*dtevp+gam1*dtevp*(u_ocean(i)*sin_io+v_ocean(i)*cos_io)
     
      tmp1=1+gam1*dtevp*cos_io
      tmp2=dtevp*(cori_nd+gam1*sin_io)
      delta=tmp1*tmp1+tmp2*tmp2
      if(delta<=0) call parallel_abort('ice_evp: delta<=0')
      u_ice(i)=(rx*tmp1+ry*tmp2)/delta
      v_ice(i)=(-rx*tmp2+ry*tmp1)/delta

      !Debug
      !if(isub==evp_rheol_steps.and.it_main==1) then
      !  write(92,*)real(xnd(i)),real(ynd(i)),real(u_ice(i)),real(v_ice(i))
      !endif
    enddo !i=1,np
    call exchange_p2d(u_ice)
    call exchange_p2d(v_ice)
  enddo !isub: sub-cycling

  !Check NaN
!  do i=1,nea
!    if(sigma11(i)/=sigma11(i).or.sigma12(i)/=sigma12(i).or.sigma22(i)/=sigma22(i)) then
!      write(errmsg,*)'NaN in ice_evp (1):',ielg(i),sigma11(i),sigma12(i),sigma22(i)
!      call parallel_abort(errmsg)
!    endif
!  enddo !i
!  do i=1,npa
!    if(u_ice(i)/=u_ice(i).or.v_ice(i)/=v_ice(i)) then
!      write(errmsg,*)'NaN in ice_evp (2):',iplg(i),u_ice(i),v_ice(i)
!      call parallel_abort(errmsg)
!    endif
!  enddo !i

  !Debug
!  if(abs(time_stamp-rnday*86400)<0.1) then
!    fdb='iceuv_0000'
!    lfdb=len_trim(fdb)
!    write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
!    open(10,file='outputs/'//fdb,status='replace')
!    write(10,*)np,nproc
!    do i=1,np
!      write(10,'(i11,3(1x,e20.12))')iplg(i),xnd(i),ynd(i),u_ice(i),v_ice(i)
!    enddo !i
!    close(10)
!
!    fdb='icedel_0000'
!    lfdb=len_trim(fdb)
!    write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
!    open(10,file='outputs/'//fdb,status='replace')
!    write(10,*)ne
!    do i=1,ne
!      write(10,*)ielg(i),delta_ice(i)
!    enddo !i
!    close(10)
!  endif

end subroutine ice_evp
