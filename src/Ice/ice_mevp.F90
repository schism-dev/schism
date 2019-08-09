! Modified EVP dynamics of ice module (Bouillon 2013)
subroutine ice_mevp
  use schism_glbl,only: rkind,time_stamp,eta2,np,npa,ne,nea,dldxy,elnode,i34,cori, &
 &grav,isbnd,nne,indel,area,iself,time_stamp,rnday,fdb,lfdb,xnd,ynd,iplg,ielg, &
 &elside,mnei,rho0,idry,errmsg
  use schism_msgp, only: myrank,nproc,parallel_abort,exchange_p2d
  use ice_module
  implicit none
  integer :: isub,i,j,ie,n1,n2,icount,m,id
  real(rkind) :: tmp1,tmp2,pp0,delta,rr1,rr2,rr3,sig1,sig2,x10,x20,y10,y20,rl10, &
 &rl20,sintheta,bb1,bb2,h_ice_el,a_ice_el,h_snow_el,dsig_1,dsig_2,mass, &
 &cori_nd,umod,gam1,rx,ry,rxp,ryp,eps11,eps12,eps22, &
 &zeta,delta_nd,sum1,sum2,dt_by_mass

  integer :: iball(mnei)
  real(rkind) :: swild(2,3),deta(2,nea), &
 &swild2(nea),alow(4),bdia(4),rrhs(3,4),u_ice_0(npa),v_ice_0(npa)

  !Save u^n
  u_ice_0=u_ice; v_ice_0=v_ice
  do isub=1,mevp_rheol_steps !iterations
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
      zeta=pp0/max(delta_ice(i),delta_min) !actually 2*zeta

      rr1=zeta*(eps11+eps22-delta_ice(i)) !part of RHS for 1st eq
      rr2=zeta*(eps11-eps22)/ellipse/ellipse
      rr3=zeta*eps12/ellipse/ellipse
      sig1=sigma11(i)+sigma22(i) !from previous iteration
      sig2=sigma11(i)-sigma22(i)

      sig1=sig1+(rr1-sig1)/mevp_alpha1
      sig2=sig2+(rr2-sig2)/mevp_alpha1

      sigma12(i)=sigma12(i)+(rr3-sigma12(i))/mevp_alpha1
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

    !Coriolis @ elem
    do i=1,nea
      swild2(i)=sum(cori(elside(1:i34(i),i)))/i34(i)
      !Pressure gradient 
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
      endif
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
      !mass=max(mass,9.d0*ice_tr(2,i)) !limit m/a>=9
      !Coriolis @ node
      cori_nd=dot_product(weit_elem2node(1:nne(i),i),swild2(iball(1:nne(i))))
      !Debug
      !if(isub==1.and.it_main==1) then
      !  write(93,*)i,real(xnd(i)),real(ynd(i)),real(),real(cori_nd) !,real(h_ice_nd),real(mass)
      !endif

      !RHS
!Error: tri only
      umod=sqrt((u_ice(i)-u_ocean(i))**2+(v_ice(i)-v_ocean(i))**2)
      dt_by_mass=dt_ice/mass
      gam1=ice_tr(2,i)*dt_by_mass*cdwat*rho0*umod
      rx=mevp_alpha2*u_ice(i)+u_ice_0(i)+gam1*(u_ocean(i)*cos_io-v_ocean(i)*sin_io)+ &
    &dt_by_mass*ice_tr(2,i)*stress_atm_ice(1,i)
      ry=mevp_alpha2*v_ice(i)+v_ice_0(i)+gam1*(u_ocean(i)*sin_io+v_ocean(i)*cos_io)+ &
    &dt_by_mass*ice_tr(2,i)*stress_atm_ice(2,i)

      !Pressure gradient
      sum1=0; sum2=0
      do j=1,nne(i)
        ie=indel(j,i)
        h_ice_el=sum(ice_tr(1,elnode(1:3,ie)))/3.0
        h_snow_el=sum(ice_tr(3,elnode(1:3,ie)))/3.0
        tmp2=rhoice*h_ice_el+rhosnow*h_snow_el !mass @elem
        sum1=sum1+tmp2*deta(1,ie)*area(ie)/3
        sum2=sum2+tmp2*deta(2,ie)*area(ie)/3
      enddo !j      
      rx=rx-grav*dt_by_mass/area_median(i)*sum1
      ry=ry-grav*dt_by_mass/area_median(i)*sum2

      !Stress
      sum1=0; sum2=0
      do j=1,nne(i)
        ie=indel(j,i)
        id=iself(j,i)
        sum1=sum1+area(ie)*(dldxy(id,1,ie)*sigma11(ie)+dldxy(id,2,ie)*sigma12(ie))
        sum2=sum2+area(ie)*(dldxy(id,1,ie)*sigma12(ie)+dldxy(id,2,ie)*sigma22(ie))
      enddo !j
      rx=rx-dt_by_mass/area_median(i)*sum1
      ry=ry-dt_by_mass/area_median(i)*sum2
     
      tmp1=1+mevp_alpha2+gam1*cos_io
      tmp2=dt_ice*cori_nd+gam1*sin_io
      delta=tmp1*tmp1+tmp2*tmp2
      if(delta<=0) call parallel_abort('ice_mevp: delta<=0')
      u_ice(i)=(rx*tmp1+ry*tmp2)/delta
      v_ice(i)=(-rx*tmp2+ry*tmp1)/delta

      !Debug
      !if(isub==evp_rheol_steps.and.it_main==1) then
      !  write(92,*)real(xnd(i)),real(ynd(i)),real(u_ice(i)),real(v_ice(i))
      !endif
    enddo !i=1,np
    call exchange_p2d(u_ice)
    call exchange_p2d(v_ice)
  enddo !isub: iteration

  !Check NaN
  do i=1,nea
    if(sigma11(i)/=sigma11(i).or.sigma12(i)/=sigma12(i).or.sigma22(i)/=sigma22(i)) then
      write(12,*)'u_ice:',u_ice,v_ice
      write(errmsg,*)'NaN in ice_mevp (1):',ielg(i),sigma11(i),sigma12(i),sigma22(i)
      call parallel_abort(errmsg)
    endif
  enddo !i
  do i=1,npa
    if(u_ice(i)/=u_ice(i).or.v_ice(i)/=v_ice(i)) then
      write(errmsg,*)'NaN in ice_mevp (2):',iplg(i),u_ice(i),v_ice(i)
      call parallel_abort(errmsg)
    endif
  enddo !i

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

end subroutine ice_mevp
