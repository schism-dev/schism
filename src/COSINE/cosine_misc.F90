!--------------------------------------------------------------------
!routines and functions for BGC models
!--------------------------------------------------------------------
!get_param_1D: read 1D parameter array
!pt_in_poly:   check whether point-in-polygon (from Utility/UtilLib) 
!o2flux:       compute air-sea o2 exchange
!co2flux:      compute air-sea co2 exchange
!ceqstate:     compuate seawater density
!ph_zbrent:    Brent's method to compute ph value
!ph_f:         nonlinear equation value of PH

include 'misc_module.txt'

contains

#ifdef USE_COSINE

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

#endif

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
   !            arco(3), nodel(3) : area coord. and 3 local node indices (valid
   !            only if inside)
         implicit real*8(a-h,o-z)
         integer, intent(in) :: i34
         real(rkind), intent(in) :: x(i34),y(i34),xp,yp
         integer, intent(out) :: inside,nodel(3)
         real(rkind), intent(out) :: arco(3)
   
         !Local
         integer :: list(3)
         real(rkind) :: ar(2),swild(2,3)
   
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

   subroutine o2flux(exflux,Tc,S,DOX,Uw)
   !-----------------------------------------------------------------
   !calculate O2 flux
   !Input:
   !     1)Tc: temperature in [C]
   !     2)S: salinity in [PSU]
   !     3)DOX: dissolved oxygen concentration [mmol/m3]
   !     4)Uw: wind velocity speed [m/s]
   !Output:
   !     1)exflux: O2 flux from air to water [m*mol/day.m2]
   !-----------------------------------------------------------------
     implicit none
     real(rkind),intent(in) :: Tc,S,DOX,Uw
     real(rkind),intent(out) :: exflux 
   
     !local variables
     real(rkind) :: Ts,eC0,C0,rho_m,rho,DOs,Sc,KwO2
   
     !pre-calculation
     call ceqstate(rho_m,Tc,S,0.0d0); rho=rho_m/1000;  !sea water density [kg/L]
     !calculate O2 saturation concentration (Garcia and Gordon 1992, Limnology and
     !Oceanography)
     Ts=log((298.15d0-Tc)/(273.15+Tc))
     eC0=5.80818d0+3.20684*Ts+4.11890d0*Ts*Ts+4.93845d0*Ts**3+1.01567d0*Ts**4 &
        & +1.41575d0*Ts**5-S*(7.01211d-3+7.25958d-3*Ts+7.93334d-3*Ts*Ts  &
        & +5.54491d-3*Ts**3)-1.32412d-7*S*S
     C0=exp(eC0)  !unit in [umol/kg]
     DOs=C0*rho  !O2 saturation concentration in [mmol/m3]
   
     !gas exchange coefficient (Keeling 1998,Global Biogeochemical Cycles)
     Sc=1638.d0-81.83d0*Tc+1.483d0*Tc*Tc-8.004d-3*Tc**3
     KwO2=0.39d0*Uw*Uw*sqrt(660.d0/Sc)
   
     exflux=KwO2*(DOs-DOX)*0.24d0  !unit in [m.mmol/day.m3]
   
     return
   end subroutine o2flux
   
   subroutine co2flux(imed,ph,exflux,Tc,S,Ct_in,Sit_in,Pt_in,Alk_in,pco2,Uw)
   !-----------------------------------------------------------------
   !written by Zhengui Wang on Jun 20, 2017
   !1)calculate co2 flux, 2)calculate ph value
   !
   !Input: 
   !     1)Tc: temperature in [C]
   !     2)S: salinity in [PSU]
   !     3)Ct_in: dissolved inorganic carbon [mol/L]
   !     4)Sit_in: Silicate [mol/L]
   !     5)Pt_in: phosphate [mol/L] 
   !     6)Alk_in: Alkalinity [eq/L]
   !     7)pco2: atmospheric CO2 concentration [ppm]
   !     8)Uw: wind velocity speed [m/s]
   !     9)imed: Ph calculation method (imed=1: simple form; imed=2: complete form) 
   !Output:
   !     1)ph: ph value
   !     2)exflux: co2 flux from air to water [mmol/day.m2]
   !-----------------------------------------------------------------
     implicit none
     integer, intent(in) :: imed
     real(rkind),intent(in) :: Tc,S,Ct_in,Sit_in,Pt_in,Alk_in,pco2,Uw
     real(rkind),intent(out) :: ph,exflux
   
     !local variables
     integer :: ierr
     real(rkind) :: rho_m,rho,T,TI,TL,T2,S05,S15,S2,I,I05,I15,I2
     real(rkind) :: ek0,ek1,ek2,ekb,ek1p,ek2p,ek3p,eksi,ekw,eks,ekf
     real(rkind) :: Ct,Sit,Pt,Alk,Bt,St,Ft,k0,k1,k2,kb,k1p,k2p,k3p,ksi,kw,ks,kf
     real(rkind) :: h,a,f0,f1,f2,CO2a,CO2w,Sc,KwCO2
   
     !pre-calculation
     T=273.15+Tc   
     TI=1.d0/T
     TL=log(T)
     T2=T*T
   
     S05=sqrt(S)
     S15=S*S05
     S2=S*S
     
     I=19.924d0*S/(1.d3-1.005d0*S)
     I05=sqrt(I)
     I15=I*I05
     I2=I*I
   
     !concentration for borate, sulfate, and fluoride [mol/Kg]
     Bt=4.16d-4*S/35d0   !(uppstrom 1974, DOE 1994)
     St=0.02824d0*S/35d0 !(Morris and Riley 1966, DOE 1994)
     Ft=7d-5*S/35d0      !(Riley 1965, DOE 1994)
   
     !convert form mol/L to mol/Kg
     call ceqstate(rho_m,Tc,S,0.0d0) 
     rho=rho_m/1.d3 !sea water density [Kg/L]
     Ct=Ct_in/rho
     Sit=Sit_in/rho
     Pt=Pt_in/rho
     Alk=Alk_in/rho
   
     !----------------------------------------------
     !calcuate reaction constants for ph calculation
     !k0,k1,k2,kb,k1p,k2p,k3p,ksi,kw,ks,kf
     !----------------------------------------------
   
     !CO2 solubility (Weiss,1974)
     ek0=9345.17d0*TI-60.2409d0+23.3585d0*log(T/100)+ &
        & S*(0.023517d0-2.3656d-4*T+4.7036d-7*T2)
     k0=exp(ek0)
    
     !carbonate dissociation (Millero, 1995) 
     !k1=[H+][HCO3-]/[H2CO3]
     !k2=[H+][CO3-]/[HCO3]
     !(Millero,1995, 'Mehrbach')
     !pk1=3670.7d0*TI-62.008d0+9.7944d0*TL-1.118d-2*S+1.16d-4*S2
     !pk2=1394.7d0*TI+4.777d0-1.84d-2*S+1.18d-4*S2
     !(Roy et al. 1993)
     ek1=2.83655d0-2307.1266*TI-1.5529413*TL-(0.207608410d0+4.0484d0*TI)*S05+ &
        & 0.0846834d0*S-6.54208d-3*S15+log(1.d0-1.005d-3*S)
     ek2=-9.226508d0-3351.6106d0*TI-0.2005743d0*TL-(0.106901773d0+23.9722d0*TI)*S05+ &
        & 0.1130822d0*S-0.00846934d0*S15+log(1.d0-1.005d-3*S) 
     k1=exp(ek1); k2=exp(ek2)
   
     !boric acid (millero,1995 <= Dickson, 1990)
     ekb=(-8966.9d0-2890.53d0*S05-77.942d0*S+1.728d0*S15-9.96d-2*S2)*TI+148.0248d0+ &
        & 137.1942d0*S05+1.621142d0*S-(24.4344d0+25.085d0*S05+0.2474d0*S)*TL+0.053105*S05*T
     kb=exp(ekb)
   
     !Water (DOE 1994, Miller 1995)
     ekw=148.96502d0-13847.26d0*TI-23.6521d0*TL+(118.67d0*TI-5.977d0+1.0495d0*TL)*S05-0.01615d0*S
     kw=exp(ekw)
    
     !Bisulfate
     eks=-4276.1d0*TI+141.328d0-23.093d0*TL+(-13856d0*TI+324.57d0-47.986d0*TL)*I05+ & 
        & (35474d0*TI-771.54d0+114.723d0*TL)*I-2698d0*TI*I15+1776d0*TI*I2+log(1.d0-1.005d-3*S)
     ks=exp(eks)
   
     if(imed==2) then 
       !phosphoric acid (DOE,1994)
       !k1p=[H][H2PO4]/[H3PO4]
       !k2p=[H][HPO4]/[H2PO4]
       !k3p=[H][PO4]/[HPO4]
       ek1p=-4576.752d0*TI+115.525d0-18.453d0*TL+(-106.736d0*TI+0.69171d0)*S05+(-0.65643d0*TI-0.01844d0)*S
       ek2p=-8814.715d0*TI+172.0883d0-27.927d0*TL+(-160.34d0*TI+1.3566d0)*S05+(0.37335*TI-0.05778d0)*S
       ek3p=-3070.75d0*TI-18.141d0+(17.27039d0*TI+2.81197d0)*S05+(-44.99486*TI-0.09984d0)*S
       k1p=exp(ek1p); k2p=exp(ek2p); k3p=exp(ek3p)
   
       !silicic acid (DOE 1994, Millero 1995)
       eksi=-8904.2d0*TI+117.385d0-19.334d0*TL+(3.5913d0-458.79d0*TI)*I05+(188.74d0*TI-1.5998d0)*I+ &
           & (0.07871d0-12.1652d0*TI)*I2+log(1.d0-1.005d-3*S)
       ksi=exp(eksi)
   
       !Hydrogen fluoride (DOE 1994, Dickson and Riley 1979) 
       ekf=1590.2d0*TI-12.641d0+1.525d0*I05+log(1.d0-1.005d-3*S)+log(1.d0+St/ks)
       kf=exp(ekf)
     endif!imed
     !-----------------------------------------------------
   
     !calculate [H+] using Brent's method
     call ph_zbrent(ierr,imed,ph,k1,k2,kb,k1p,k2p,k3p,ksi,kw,ks,kf,Bt,St,Ft,Ct,Sit,Pt,Alk)
     !write(*,*) ierr,imed,ph,k1,k2,kb,k1p,k2p,k3p,ksi,kw,ks,kf,Bt,St,Ft,Ct,Sit,Pt,Alk,rho
   
     !-----------------------------------------------------
     !calculate CO2 flux
     !-----------------------------------------------------
     h=10.d0**(-ph)
     CO2a=k0*pco2*1.d-6 !saturation CO2 concentration [mol/kg]
     CO2w=Ct*h*h/(h*h+k1*h+k1*k2); !aqueous CO2 concentration [mol/kg]
     !Schmidt number of CO2 (Wanninkhof 1992, JGR)
     Sc=2073.1d0-125.62d0*Tc+3.6276d0*Tc*Tc-0.043219d0*Tc**3
     KwCO2=0.39d0*Uw*Uw*sqrt(660.d0/Sc) !unit in [cm/hr]
   
     exflux=KwCO2*(CO2a-CO2w)*(0.24d0*rho*1.d6);  !unit in [m.mmol/day.m3]
     
     return
   end subroutine co2flux
   
   subroutine ceqstate(rho,T,S,P)
   !--------------------------------------------------------------------
   !written by Zhengui Wang on Jun 19, 2017
   !calculate seawater density (Millero and Poisson 1981, Gill, 1982)
   !Input:
   !    1)T: temperature in [C]
   !    2)S: salinity in [PSU]
   !    3)P: pressure [bars] (P=0: at 1 atm. pressure)
   !Output:
   !    1)rho: seawater density in [Kg/m3]
   !--------------------------------------------------------------------
     implicit none
     real(rkind), intent(in) :: T,S,P
     real(rkind), intent(out) :: rho
   
     !local 
     real(rkind) :: T2,T3,T4,T5,S05,S15,S2,P2
     real(rkind) :: rho_pw,rho_st,A,B,C,K_pw,K_st,K_stp
   
     !pre_calculation
     T2=T*T
     T3=T**3
     T4=T**4
     T5=T**5
     S05=sqrt(S)
     S15=S*S05
     S2=S*S
     P2=P*P
   
     !pure water S=0,at 1 atm.  
     rho_pw=999.842594d0+6.793952d-2*T-9.095290d-3*T2+1.001685d-4*T3 &
           & -1.120083d-6*T4+6.536332d-9*T5
   
     !density with Salinity
     A=0.0; B=0.0; C=0.0
     if(S/=0.0) then
       A=8.24493d-1-4.0899d-3*T+7.6438d-5*T2-8.2467d-7*T3+5.3875d-9*T4 
       B=-5.72466d-3+1.0227d-4*T-1.6546d-6*T2
       C=4.8314d-4
     endif !S 
     rho_st=rho_pw+A*S+B*S15+C*S2
   
     !pressure not zero
     if(P/=0.0) then
       K_pw=19652.21d0+148.4206d0*T-2.327105d0*T2+1.360477d-2*T3-5.155288d-5*T4
       K_st=K_pw
       if(S/=0.0) then
         K_st=K_st+S*(54.6746d0-0.603459d0*T+1.09987d-2*T2-6.1670d-5*T3) &
             & +S15*(7.944d-2+1.6483d-2*T-5.3009d-4*T2);
       endif !S
       K_stp=K_st+P*(3.239908d0+1.43713d-3*T+1.16092d-4*T2-5.77905d-7*T3) &
            & +P*S*(2.2838d-3-1.0981d-5*T-1.6078d-6*T2)+1.91075d-4*P*S15 &
            & +P2*(8.50935d-5-6.12293d-6*T+5.2787d-8*T2) &
            & +P2*S*(-9.9348d-7+2.0816d-8*T+9.1697d-10*T2)
       rho=rho_st/(1.d0-P/K_stp)
     else
       rho=rho_st
     endif !P
     
     return
   end subroutine ceqstate

   subroutine ph_zbrent(ierr,imed,ph,k1,k2,kb,k1p,k2p,k3p,ksi,kw,ks,kf,Bt,St,Ft,Ct,Sit,Pt,Alk)
   !---------------------------------------------------------------------
   !Brent's method to find ph value
   !numerical recipes from William H. Press, 1992
   !
   !Input:
   !    1)k1,k2,kb,k1p,k2p,k3p,ksi,kw,ks,kf: these are reaction constant 
   !    2)Bt,St,Ft,Ct,Sit,Pt,Alk: borate,sulfate,fluoride,carbon,silicate,phosphate [mol/kg]
   !    3)imed: method (imed=1: simple form; imed=2: complete form)
   !Output:
   !    1)ph: ph value
   !    2)ierr: error flag 
   !---------------------------------------------------------------------
     implicit none
     integer, parameter :: nloop=100
     real(kind=rkind), parameter :: eps=3.0d-8, tol=1.d-6,phmin=3.0,phmax=13.0
     integer, intent(in) :: imed
     integer, intent(out) :: ierr
     real(kind=rkind),intent(in) :: k1,k2,kb,k1p,k2p,k3p,ksi,kw,ks,kf,Bt,St,Ft,Ct,Sit,Pt,Alk
     real(kind=rkind),intent(out) :: ph
   
     !local variables
     integer :: i
     real(kind=rkind) :: a,b,c,d,e,m1,m2,fa,fb,fc,p,q,r,s,tol1,xm
     real(kind=rkind) :: rtmp,h
   
     !initilize upper and lower limits
     ierr=0
     a=phmin
     b=phmax
   
     h=10.0**(-a); call ph_f(fa,imed,h,k1,k2,kb,k1p,k2p,k3p,ksi,kw,ks,kf,Bt,St,Ft,Ct,Sit,Pt,Alk)
     h=10.0**(-b); call ph_f(fb,imed,h,k1,k2,kb,k1p,k2p,k3p,ksi,kw,ks,kf,Bt,St,Ft,Ct,Sit,Pt,Alk)
   
     !root must be bracketed in brent
     if(fa*fb>0.d0) then
       ierr=1
       return
     endif
   
     fc=fb
     do i=1,nloop
       if(fb*fc>0.d0) then
         c=a
         fc=fa
         d=b-a
         e=d
       endif !fb*fc>0.d0
       if(abs(fc)<abs(fb)) then
         a=b
         b=c
         c=a
         fa=fb
         fb=fc
         fc=fa
       endif !abs(fc)
       tol1=2.d0*eps*abs(b)+0.5d0*tol !convergence check
       xm=0.5d0*(c-b)
       if(abs(xm)<=tol1.or.fb==0.d0) then
       !if(abs(xm)<=tol1.or.abs(fb)<=1.d-8) then
         ph=b
         return
       endif
       if(abs(e)>=tol1.and.abs(fa)>abs(fb)) then
         s=fb/fa
         if(a==c) then
           p=2.d0*xm*s
           q=1.d0-s
         else
           q=fa/fc
           r=fb/fc
           p=s*(2.d0*xm*q*(q-r)-(b-a)*(r-1.d0))
           q=(q-1.d0)*(r-1.d0)*(s-1.d0)
         endif !a==c
         if(p>0.d0) q=-q
         p=abs(p)
         m1=3.d0*xm*q-abs(tol1*q)
         m2=abs(e*q)
         if(2.d0*p<min(m1,m2)) then
           e=d
           d=p/q
         else
           d=xm
           e=d
         endif !2.d0*p<min
       else
         d=xm
         e=d
       endif !abs(e)
       a=b;
       fa=fb
       if(abs(d)>tol1) then
         b=b+d
       else
         b=b+sign(tol1,xm)
       endif !abs(d)
       h=10.0**(-b); 
       call ph_f(fb,imed,h,k1,k2,kb,k1p,k2p,k3p,ksi,kw,ks,kf,Bt,St,Ft,Ct,Sit,Pt,Alk)
     enddo !i
   
     ierr=2
     ph=b
   
     return
   
   end subroutine ph_zbrent
   
   subroutine ph_f(f,imed,h,k1,k2,kb,k1p,k2p,k3p,ksi,kw,ks,kf,Bt,St,Ft,Ct,Sit,Pt,Alk)
   !--------------------------------------------------------------------
   !calculate the nonlinear equation value of PH
   !--------------------------------------------------------------------
     implicit none
     !integer, parameter :: rkind=rk
     integer, intent(in) :: imed
     real(kind=rkind), intent(in) :: h,k1,k2,kb,k1p,k2p,k3p,ksi,kw,ks,kf,Bt,St,Ft,Ct,Sit,Pt,Alk
     real(kind=rkind), intent(out):: f
   
     if(imed==1) then !simple form!
       !f=[HCO3]+2[CO3]+[B(OH)4]+[OH]+-[H]_f-Alk
       f=(h+2.0*k2)*Ct*k1/(h*h+k1*h+k1*k2)+Bt*kb/(h+kb)+kw/h-h*ks/(St+ks)-Alk
     elseif(imed==2) then !only boric
       !f=[HCO3]+2[CO3]+[B(OH)4]+[OH]+[HPO4]+2[PO4]-[H3PO4]+[H3SiO4]-[H]_f-[HF]-[HSO4]-Alk
       f=(h+2.0*k2)*Ct*k1/(h*h+k1*h+k1*k2)+Bt*kb/(h+kb)+kw/h+ &
        & Pt*((h+2*k3p)*k1p*k2p-h**3)/(h**3+k1p*h*h+k1p*k2p*h+k1p*k2p*k3p)+ &
        & Sit*ksi/(h+ksi)-h*ks/(St+ks)-Ft*h/(h+kf)-St*h/(h+St+ks)-Alk
     else
       stop 'unknown imed in PH calculation'
     endif
   
     return
   end subroutine ph_f

end module 
