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

!  Adapted from FESOM's ice module.
!  Changed by Joseph Zhang marked by 'YJZ'

!  subroutine ice_thermodynamics (driver)
!  subroutine therm_ice
!  subroutine ocean_budget
!  subroutine ice_budget
!  function TFrez
! ==================================================================
subroutine ice_thermodynamics
  ! For every surface node, this subroutine extracts the information
  ! needed for computation of thermodydnamics, calls the relevant
  ! subroutine, and returns the information to the vectors of prognostic
  ! variables.
  !
  ! Originally programmed by N. Yakovlev/S. Danilov.
  ! Adjusted for upgraded model physics (NCEP forcing data; parameterization
  ! of ice-ocean heat flux considering friction velocity) by Ralph Timmermann.
  ! Note that adjustments are not indivually labeled any more.
  use schism_glbl,only: rkind,npa,tr_nd,iplg,pr,fluxprc,rho0,windx,windy, &
  &nvrt,srad,albedo,hradd,airt2,shum2,errmsg,xlon,ylat,fresh_wa_flux,net_heat_flux
  use schism_msgp, only: myrank,nproc,parallel_abort,parallel_finalize
  use ice_module
  use ice_therm_mod

  implicit none
  real(rkind) :: h,hsn,A,t,fsh,flo,Ta,qa,preslev,acl,precrate
  real(rkind) :: ug,ustar,T_oc,S_oc,fw,ehf,srad2
  integer :: i,j

  do i=1,npa 
!    h=m_ice(i)
!    hsn=m_snow(i)
!    A=a_ice(i)
!    fsh=shortwave(i)
!    flo=longwave(i)
!    Ta=Tair(i)
!    td=Tdew(i)
!    qa=shum(i)  !specific humidity
!    acl     = cloudiness(i) !not used
!    t_oi(i)=0. !T@ ice-ocean interface [C]
!    fw=0
!    ehf=0.

    if(ice_tests==1) then
      preslev=1.e3 !hPa
      precrate=1.e-5 
      T_oc=-1
      S_oc=34 
      srad2=1
      hradd(i)=1
      airt2(i)=0 !C
      shum2(i)=1.e-3
!      ustar=0
      ug=5 !test
    else
      preslev=pr(i)/100 !air pressure in hPa      
      precrate=fluxprc(i)/rho0 !precipitation [m water/s]
      !In dry ocean case, use last wet values
      T_oc=tr_nd(1,nvrt,i) !T@ mixed layer - may want to average top layers??
      S_oc=tr_nd(2,nvrt,i) !S
      srad2=srad(i)/(1-albedo(i)) !(denom checked) solar without albedo
      ug=sqrt(windx(i)**2+windy(i)**2)
    endif !ice_tests

    ustar=sqrt(cdwat*((u_ice(i)-u_ocean(i))**2+(v_ice(i)-v_ocean(i))**2)) !ice-ocean frictional vel
!YJZ
!    h_ml=0.1 !Mixed layer (ML) depth [m]

!    write(99,*)'doing node:',i,ustar,ug,u_ice(i),v_ice(i),u_ocean(i),v_ocean(i),ice_tr(1:3,i),&
!    &dt_ice,srad(i),hradd(i),airt2(i),shum2(i),preslev,precrate,T_oc,S_oc,h_ml0

    call therm_ice(ice_tr(1,i),ice_tr(3,i),ice_tr(2,i),t_oi(i),srad2,hradd(i),airt2(i), &
  &shum2(i),preslev,precrate,ug,ustar,T_oc,S_oc,h_ml0,dt_ice,iplg(i), &
  &fresh_wa_flux(i),net_heat_flux(i))

!    write(99,*)'After therm_ice:',t_oi(i),fresh_wa_flux(i),net_heat_flux(i),ice_tr(1:3,i)

!    m_ice(i)=h
!    m_snow(i)=hsn
!    a_ice(i)=A
!    fresh_wa_flux(i)=fw
!    net_heat_flux(i)=ehf

    !Check outputs
    if(t_oi(i)/=t_oi(i).or.fresh_wa_flux(i)/=fresh_wa_flux(i).or. &
  &net_heat_flux(i)/=net_heat_flux(i).or.ice_tr(1,i)/=ice_tr(1,i).or. &
  &ice_tr(2,i)/=ice_tr(2,i).or.ice_tr(3,i)/=ice_tr(3,i)) then
      write(errmsg,*)'ice_thermodynamics: nan',iplg(i),t_oi(i), &
  &fresh_wa_flux(i),net_heat_flux(i),ice_tr(:,i)
      call parallel_abort(errmsg)
    endif
  enddo !i=1,npa

!  write(99,*)'done'
!  call parallel_finalize
!  stop

end subroutine ice_thermodynamics

!====================================================================
subroutine therm_ice(h,hsn,A,toi,fsh,flo,tair,qa,preslev,precrate, &
  &ug,ustar,T_oc,S_oc,h_ml,dt,inode,fw,ehf)

  !===========================================================
  !             Ice Thermodynamic growth model               !
  !===========================================================

  ! Input parameters:
  !------------------
  ! h - ice mass [m]
  ! hsn - snow mass [m]
  ! A - ice compactness [-]
  ! toi - temperature of snow/ice top surface [C] (T_sfc in P&W)
  ! fsh - shortwave radiation [W/m/m] (before albedo correction) (positive down)
  ! flo - longwave radiation [W/m/m] (positive down)
  ! tair - air temperature [C]
  ! qa - specific humidity [-]
  ! preslev - air pressure [hPa]
  ! precrate - precipitation rate [m water /s]
  ! ug - wind speed [m/s]
  ! T_oc, S_oc - ocean temperature and salinity beneath the ice (mixed layer) [C,PSU]
  !              In case of dry, use last wet values
  ! h_ml - mixed layer depth - should be specified [m]
  ! dt - ice time step [sec]
  ! ustar - friction velocity [m/s]
  ! inode - global node # (info only)
  ! acl - not used
  ! td - air dew point temp [C] - not needed
  
  ! Output parameters:
  !-------------------
  ! h - ice mass
  ! hsn - snow mass
  ! A - ice compactness 
  ! toi - temperature of snow/ice top surface [C]
  ! fw - freshwater flux due to ice melting [m water/dt]; >0 to freshen the
  ! ocean surface
  ! ehf - net heat flux at the ocean surface [W/m/m]; >0 to warm the ocean

  use schism_glbl,only: rkind
  use ice_therm_mod
  use ice_module, only: salt_ice,salt_water

  implicit none
  integer, intent(in) :: inode
  real(rkind), intent(inout) :: h,hsn,A,toi
  real(rkind), intent(in) :: fsh,flo,tair,qa,preslev,precrate,ug,ustar,T_oc,S_oc,h_ml,dt
  real(rkind), intent(out) :: fw,ehf

  real(rkind) :: dhgrowth,dhsngrowth,ahf,prec,fh,gh,hfsen,hfsenow,hfrad,hfradow,hflatow,hftotow, & 
 &hflat,hftot,tmp,tmp2,tmp3
  real(rkind) :: rhow,show,rhice,shice,sh,thick,thact,thsn,hdraft,hflood
  real(rkind) :: rh,ra,qhst,qfm,qtm,sn,snow,hsntmp,qhsto,qt,o2ihf

  integer :: k

  real(rkind),external :: TFrez  ! Sea water freeze temperature.

  !------------------------------------------------------------------
  ! determine h(i,j)/a(i,j) = actual ice thickness.
  ! if snow layer is present, add hsn weighted with quotient
  ! of conductivities of ice and snow, according to 0-layer approach
  ! of Semtner (1976).
  !------------------------------------------------------------------
  
  dhgrowth=h  ! Store ice thickness at start of growth routine
  thsn=hsn*(con/consn)/max(A,Armin)    ! Effective snow thickness
  thick=h/max(A,Armin)+thsn            ! Effective ice thickness
  !!!!!!!thick=max(thick,hmin)  ! Now there is IF block 28.10.2003.
  rhice=0.0                            ! Init? Reset average heat flux
  
  !write(*,*) T_oc vor ocean_budget,T_oc
  !No ice cover
  !rhow: Growth rate for thin ice (open ocean); >0: ice grows [m/s]
  call ocean_budget(preslev,qa,fsh,flo,T_oc,ug,tair,rhow)

!  write(99,*)'enter therm_ice:',dhgrowth,thsn,thick,rhow,h,A,Armin,hmin
  
  !Ice cover with no snow?
  if (thick>hmin) then
    ! K loop for all ice classes
    ! write(mype+30,*) 'iclasses',istep,inode,iclasses,thick
    do k=1,iclasses ! K loop
      thact = (2*k-1)*thick/dble(iclasses)  ! Thicknesses of actual ice class
      !  write(mype+30,*) 'in iclasses',&
      !   thact,hsn,t,tair,qa,preslev,fsh,flo,ug,S_oc,shice
      call ice_budget(thact,hsn,toi,tair,qa,preslev,fsh,flo,ug,S_oc,shice) ! Thick ice K-class growth rate
      rhice=rhice+shice/dble(iclasses)      ! Add to average heat flux
    end do ! k loop
  END IF
!  write(99,*) 'after iclasses',inode,iclasses,rhice
  
  !---------------------------------------------------------------
  ! Convert growth rates [m ice/sec] into growth per time step DT.
  !---------------------------------------------------------------
  rhow=rhow*dt !>0: ice grows  [m]
!  rA=rhow ! for historical reasons {remove}
  rhice=rhice*dt !>0: ice grows?  [m]
  
  !---------------------------------------------------------------
  ! Multiply ice growth of open water and ice-covered areas
  ! with the corresponding areal fractions of grid cell
  !---------------------------------------------------------------
  show=rhow*(1.0-A) !m
  shice=rhice*A
  sh=show+shice !m
  
  ahf=-cl*sh/dt    !Atmospheric heat flux, average over grid cell [W/m**2]
                   !ahf<0 (sh>0): ice grows

!  write(99,*)'STEP1:',rhow,rhice,show,shice,sh,ahf
  
  if(tair>=0.) then ! Precipitation going to ocean
    prec=precrate*dt ![m]
    snow=0.
  else
    prec=precrate*dt*(1.-A)
    snow=precrate*dt*A*rhowat/rhosnow ![m]
  endif
  
  hsn=hsn+snow     ! Add snow fall to temporary snow thickness [m]
  dhsngrowth=hsn   ! Store snow thickness after snow fall 

!  write(99,*)'STEP2:',prec,snow,hsn,dhsngrowth
  
  ! ---------------------------------------------------------------------
  ! If there is atmospheric melting, first melt any snow that is present.
  ! Atmospheric heat flux available for melting
  ! sh = MINUS atm. heat flux / specific latent heat of sea ice
  ! Note: (sh<0) for melting, (sh>0) for freezing
  !----------------------------------------------------------------------
  hsntmp=-min(sh,0.)*rhoice/rhosnow !>=0; melting amount; [m]
  
  !--------------------------------------------------------------------
  ! hsntmp is the decrease in snow thickness due to atmospheric melting
  ! [m/DT]. Do not melt more snow than available
  !--------------------------------------------------------------------
  hsntmp=min(hsntmp,hsn) !m
  hsn=hsn-hsntmp  ! Update snow thickness after atmospheric snow melt; >=0

!  write(99,*)'STEP3:',hsntmp,hsn
  
  !-------------------------------------------------------------------
  ! Negative atmospheric heat flux left after melting of snow
  ! Note: (sh<0) and (hsntmp>0) for melting conditions
  ! hsntmp=0 for no-snow-melting conditions
  ! RT-FR099: snow melt will be added to fw flux at the end of this routine
  !-------------------------------------------------------------------
  rh=sh+hsntmp*rhosnow/rhoice !>=0 [m]
  h=max(h,0.d0)
  
  !-------------------------------------------------------------------
  ! Compute heat content qhst of mixed layer - sea ice system
  !
  ! Total heat content is the sum of
  !	h	ice thickness after calculation of dynamic effects
  !	rh	change in ice thickness due to atmospheric forcing
  ! and heat available in mixed layer, with
  !	T_oc	temperature of ocean surface layer
  !	Tfrez	freezing point of sea water
  !	h_ml	thickness of uppermost layer
  !
  !RT:
  ! There are three possibilities to do this.
  ! 1.: Assume an instantaneous adjustment of mixed layer heat content.
  !     Any heat available is then instantaneously used to melt ice.
  !     (so-called ice-bath approach)
  !     This is what used to be used in the Lemke sea ice-mixed layer model.
  ! qhst=h+rh -(T_oc-TFrez(S_oc))*h_ml*cc/cl
  !
  ! 2.: Parameterize the ocean-to-ice heat flux (o2ihf)
  !     as a function of temperature difference. For a first step 
  !     we can assume a constant exchange coefficient gamma_t:
  ! o2ihf= (T_oc-TFrez(S_oc))*gamma_t*cc*A     &
  !        +(T_oc-Tfrez(S_oc))*h_ml/dt*cc*(1-A)  ! [W/m2]
  ! qhst=h+rh-o2ihf*dt/cl                        ! [m]
  !
  ! 3.  Parameterize the ocean-to-ice heat flux (o2ihf)
  !     as a function of temperature difference and the
  !     friction velocity:
  
  !Option 2
  o2ihf=(T_oc-TFrez(S_oc))*0.012*ustar*cc*A+  &
      &(T_oc-Tfrez(S_oc))*h_ml/dt*cc*(1-A)  ! [W/m2]
  qhst=h+rh-o2ihf*dt/cl  ! [m]

!  write(99,*)'STEP4:',rh,h,o2ihf,qhst

  !    
  !RT--
  !-------------------------------------------------------------------
  ! qhst > 0	temporary ice thickness [m]
  ! qhst < 0  	temporary heat stored in mixed layer (all ice melts) [m]
  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  ! New temporary ice thickness after changes due to atmospheric
  !  forcing and mixed layer heat content:
  !-------------------------------------------------------------------
  tmp=max(qhst,0.d0)
  !-------------------------------------------------------------------
  ! Melt snow if there is any ML heat content left (qhst<0).
  ! This may be the case if advection moves ice (with snow) to regions
  ! with a warm mixed layer.
  !-------------------------------------------------------------------
  sn=hsn+min(qhst,0.d0)*rhoice/rhosnow
  !-------------------------------------------------------------------
  ! When the temporary snow thickness is negative, there is some ML
  ! heat content left: 
  !-------------------------------------------------------------------
  tmp2=min(sn,0.d0)*rhosnow/rhoice !<=0 [m]
  !-------------------------------------------------------------------
  ! tmp2 <= 0 is the ice melt equivalent [m] of the mixed layer heat
  ! content left after melting of snow.
  ! New temporary snow thickness must not be negative:
  !-------------------------------------------------------------------
!  write(99,*)'STEP5:',tmp,sn,tmp2

  sn=max(sn,0.d0)
  !-------------------------------------------------------------------
  ! Add snow melt to precipitation
  !-------------------------------------------------------------------
  !RT FR099: this will be done at the end of this routine
  !RT prec=prec +(hsn-sn)*rhosnow/rhowat
  !-------------------------------------------------------------------
  ! "Add" snow melt to ML heat content (there will be no snow melt for qhst>0)
  !-------------------------------------------------------------------
  qhst=qhst+(hsn-sn)*rhosnow/rhoice
  !-------------------------------------------------------------------
  ! Additional fresh water flux from ice thickness changes
  ! due to atmospheric forcing and ML heat content [m water/DT]
  !-------------------------------------------------------------------
  qfm=(h-tmp)*rhoice/rhowat !m; can be <0?
  !rtsd  fw=-qfm  !++RT: Why minus??? prec is positive for downward fresh water
           !      water flux too. And besides, salinity of sea ice is
           !      about 5 psu, not zero. In BRIOS we use:
  !fw=qfm*(34.-5.)/34./dt !m/s
!YJZ
  if(salt_water==0) then
    fw=qfm/dt !m/s
  else
    fw=qfm*(salt_water-salt_ice)/salt_water/dt !m/s
  endif
           !++RT This is the fresh water flux that should be passed to
           !     the ocean for coupling.                           -RT--
  !-------------------------------------------------------------------
  ! Heat flux into ML =    !++RT: actually it is the heat flux FROM ML
  !  ML heat content left after melting of ice and snow
  !  minus ML heat content of previous time step
  ! Note: In coupled ocean-ice models, the heat content of the
  !	mixed layer should be calculated in the ocean part
  !qtm<0: ice melts (heat content goes down in ML), so qtm=sflux in SCHISM
  !-------------------------------------------------------------------
  qtm=-tmp2-(T_oc-TFrez(S_oc))*h_ml*cc/cl !m
!YJZ: this formulation is clearer
  ehf=qtm/dt*cl ![W/m/m]
  
  !++RT: This qtm can be used as the ocean surface heat loss.
  !++RT: However, I suggest to use ehf further down.
  
  !-------------------------------------------------------------------
  ! Update snow depth
  !-------------------------------------------------------------------
  hsn=sn !>=0

!  write(99,*)'STEP6:',sn,qhst,qfm,fw,qtm,hsn
   
  h=max(qhst,0.d0)             ! New ice thickness >=0
  if(h<1d-6) h=0.        ! Avoid very small ice thicknesses
  dhgrowth=h-dhgrowth        ! Change in ice thickness due to thermodynamic effects
  dhsngrowth=hsn-dhsngrowth  ! Change in snow thickness due to thermodynamic melting
                             ! (without snow fall).
                             ! This is a negative value (= MINUS snow melt)
  
  prec=prec/dt               ! Conversion: 'per time step' -> 'per second' [m/s]
  dhgrowth=dhgrowth/dt       ! Conversion: 'per time step' -> 'per second' [m/s]
  dhsngrowth=dhsngrowth/dt   ! Conversion: 'per time step' -> 'per second'
  
  !RT: This is what I suggest to use as the ocean surface boundary
  !RT: condition for coupling:
  !ehf>0: heat to warm the ocean
!YJZ: use the formulation above
!  ehf=ahf+cl*(dhgrowth+(rhosnow/rhoice)*dhsngrowth) !W/m/m

  !RT FR099:
  !fw<0: evaporation
!YJZ
  if(salt_water==0) then
    tmp3=1
  else
    tmp3=(salt_water-salt_ice)/salt_water
  endif
  !fw=-dhgrowth*rhoice/rhowat*(34.-5.)/34.-dhsngrowth*rhosnow/rhowat+prec !m/s
  fw=-dhgrowth*rhoice/rhowat*tmp3-dhsngrowth*rhosnow/rhowat+prec !m/s
  !RT-

!  write(99,*)'STEP7:',h,dhgrowth,dhsngrowth,prec,dhgrowth,ehf,fw
  
  !======================================================
  ! Changes in compactnesses (equation 16 of Hibler 1979)
  !======================================================
  
  rh=-min(h,-rh)   ! Make sure we do not try to melt more ice than is available
  rA= rhow -(T_oc-TFrez(S_oc))*h_ml*cc/cl !m
  
  A=A+Saterm*min(rh,0.d0)*A/max(h,hmin)+max(rA,0.d0)*(1.-A)/h0   ! Total change. [-]
  
  A=min(A,h*1.e6)     ! A -> 0 for h -> 0; impose h>=1.e-6m
  A=min(max(A,0.d0),1.d0) ! A >= 0, A <= 1

!  write(99,*)'STEP8:',rh,rA,A
  
!  call flooding(h,hsn)
  ! Flooding (snow to ice conversion)
  hdraft=(rhosnow*hsn+h*rhoice)/rhowat ! Archimedes: displaced water
  hflood=hdraft-min(hdraft,h)         ! Increase in mean ice thickness due to flooding; >=0
  h=h+hflood                          ! Add converted snow to ice volume; >=0
  hsn=hsn-hflood*rhoice/rhosnow       ! Subtract snow from snow layer
  !hsn>=0 ensured b/cos: (1) if hdraft<=h, hflood=0; (2) if hdraft>h,
  !hflood=hdraft-h, and
  !hsn_new=hsn-rhoice/rhosnow*(hdraft-h)=(hsn+h*rhoice/rhosnow)*(1-rhoice/rhow)>=0
  !sice rhoice<rhow

  !RT: This is what all AWI sea ice models do, but
  !I wonder whether it really is correct for the heat budget.
  !I suggest we initially keep it to allow for a comparison with BRIOS
  !results and rethink it at a later stage.

!  write(99,*)'After flooding:',h,hsn

end subroutine therm_ice

!==================================================================
subroutine ice_budget(hice,hsn,toi,tair,qa,preslev,fsh,flo,ug,S_oc,fh)
  !++RT This subroutine needs modification if it is supposed to work with 
  !RT   NCEP atmospheric forcing data - which it is.
  !RT   I am gonna do it in the final code once it is clear where in the code
  !RT   we read the atmospheric daily forcing data.
  !RT   I am also gonna do the final crosscheck then.
  !--RT
  
  ! Thick ice growth rate from air-ice exchange [m ice/sec]
  
  ! INPUT:
  ! hice - actual ice thickness [m]
  ! hsn - snow thickness, used for albedo parameterization [m]
  ! toi - temperature of snow/ice surface [C] (T_sfc in P&W)
  ! tair - air temperature [C]
  ! qa - specific humidity
  ! preslev - air pressure [hPa]
  ! fsh - shortwave radiation [W/m**2]
  ! flo - longwave radiation  [W/m**2]
  ! ug - wind speed [m/sec]
  ! S_oc - ocean salinity for the temperature of the ice base calculation [ppt]
  
  ! OUTPUT: 
  ! fh (growth rate); 
  ! toi: modified when thick ice occurs
  
  use schism_glbl,only: rkind
  use schism_msgp, only : parallel_abort
  use ice_therm_mod
  implicit none
  real(rkind), intent(in) :: hice,hsn,tair,qa,preslev,fsh,flo,ug,S_oc !hfsen,hfrad,hflat,hftot,S_oc
  real(rkind), intent(out) :: fh
  real(rkind), intent(inout) :: toi
  
  real(rkind) ::  alb             ! Albedo of sea ice
  real(rkind) ::  a1,a2,a3,b,c,hfsen,hfrad,hflat,hftot    
  integer :: iter, imax      ! Number of iterations
  
  real(rkind), external :: TFrez
  
  imax = 5
!  if ((qa == 1.e20).and.(td < 1.e19)) then
!   write(*,*)'compute qa from td'
!   qa = (0.622/preslev)*6.11*exp(qsi*td/(td+tqi))
!  endif
   
  !-----------------------------------------------------------------
  ! set albedo
  ! ice and snow, freezing and melting conditions are distinguished.
  !-----------------------------------------------------------------
  if(toi<0.0) then ! freezing condition    ----
    if(hsn>0.0) then !   snow cover present  -   i
      alb=albsn         !                        i  i
    else             !   no snow cover       -i  i
      alb=albi      !                        i  i
    endif !   snow cover?         -   
  else ! melting condition     ----i
    if(hsn>0.0) then !   snow cover present  -   i
      alb=albsnm !                        i  i
    else !   no snow cover       -i  i
      alb=albm !                        i  i
    endif !   snow cover?         -   i
  endif ! freezing or melting?  ----
  
  
  ! total incoming atmospheric heat flux (part of f(x0))
  a1=(1.0-alb)*fsh +flo +d1*ug*tair +d2i*ug*qa
  
  !
  ! NEWTON-RHAPSON TO GET TEMPERATURE AT THE TOP OF THE ICE LAYER (toi)
  ! Eq. (16) of Parkinson & Washington?
  do iter=1,imax
    if(toi+tqi==0) call parallel_abort('ice_budget: toi+tqi=0')
    !Changed 10^a to exp(a*log(10)) compared to P&W
    b=6.11*0.622/preslev*exp(qsi*toi/(toi+tqi)) ! specific humidity for ice
    a3=d2i*ug*b*qsi*tqi/((toi+tqi)**2)        !negative gradient coefficient for the humidity part
    a2=-d1*ug*toi-d2i*ug*b-d3*((toi+tmelt)**4)  ! remaining f(x0) from sensible, latent heat, upward long wave
    !hice>0 checked in the calling routine
    if(hice==0) call parallel_abort('ice_budget: hice=0')
    c=con/hice                              !negative gradient coefficient for heat conductivity part
    a3=a3+4.0*d3*((toi+tmelt)**3)+c+d1*ug   !negative gradient coefficient for sensible & latent heat 
    c=c*(TFrez(S_oc)-toi)                   !downward conductivity term (f(x0))
  
    if(a3==0) call parallel_abort('ice_budget: a3=0')
    toi=min(0.d0,toi+(a1+a2+c)/a3)  !Newton-Raphson
  enddo !iter
  
  !Heat fluxes [W/m**2]
  hfrad= (1.0-alb)*fsh & ! absorbed short wave radiation
                 &+flo & ! long wave radiation coming in
 &-d3*((toi+tmelt)**4) ! long wave radiation going out
  
  hfsen=d1*ug*(tair-toi)                   ! sensible heat
  hflat=d2i*ug*(qa-b)                     ! latent   heat
  hftot=hfrad+hfsen+hflat                 ! total (atm.)
  
  fh=-hftot/cl    ! growth rate [m ice/sec]; >0: ice grows
end subroutine ice_budget

!======================================================

subroutine ocean_budget(preslev,qa,fsh,flo,toc,ug,tair,fh)  
  ! Ice growth rate for open ocean [m ice/sec]
  
  ! INPUT:
  ! preslev - air pressure [hPa]
  ! qa   - specific humidity              !RT
  ! fsh - shortwave radiation
  ! flo - longwave radiation
  ! toc - temperature of open water ML [C]
  ! ug - wind speed [m/sec]
  ! tair - air temperature [C]
  
  ! OUTPUT: fh
  
  !RT: This routine can use either the dew point temperature or the specific humidity.
  !RT: Put 1.e20 as a dummy value into the field you do not want to provide.
  
  use schism_glbl,only: rkind
  use ice_therm_mod
  implicit none
  real(rkind),intent(in) :: preslev,qa,fsh,flo,toc,ug,tair
  real(rkind),intent(out) :: fh

  real(rkind) :: b,hfsenow,hfradow,hflatow,hftotow
  
  !RT specific humidity [kg/kg] 
  
  !write(*,*)'qaline',qa,td,toc,preslev,tqw,qsw
  
!  if((qa>=9.e19).and.(td< 9.e19)) then
!   write (*,*) 'compute specific humidity from dew point temperature'
!   qa= (0.622/preslev)*6.11*exp(qsw*td/(td+tqw))
!  endif
  
  !specific humidity @ surface
  !Since preslev is in hPa, Eq. (10) of Parkinson & Washington is deivided by
  !100
  b=6.11*0.622/preslev*exp(qsw*toc/(toc+tqw))
  
  !                                 Heat fluxes [W/m**2]
  hfradow= (1.0-albw)*fsh & ! absorbed short wave radiation
                    &+flo & ! long wave radiation coming in
        &-d3*((toc+tmelt)**4) ! long wave radiation going out
  
  hfsenow=d1 *ug*(tair-toc)                 ! sensible heat
  hflatow=d2w*ug*(qa-b)                   ! latent   heat
  hftotow=hfradow+hfsenow+hflatow ! total downward [W/m/m]
  
  fh= -hftotow/cl   ! growth rate [m ice/sec]
  !fh>0: ice grows
end subroutine ocean_budget

!==========================================================
function TFrez(S)
  ! Nonlinear correlation for the water freezing temperature.
  ! Millero (1978) - UNESCO. Reference - See A. Gill, 1982.
  use schism_msgp, only: parallel_abort
  implicit none
  real*8 :: S, TFrez

  if(S<0) call parallel_abort('TFrez: S<0')
  TFrez=-0.0575*S+1.7105e-3*sqrt(S**3)-2.155e-4*S*S ![C]

  !TFrez=-1.86

end function TFrez

!==========================================================

