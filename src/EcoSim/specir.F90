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
    
    SUBROUTINE spec_ir(wtratio)

    
!!======================================================================
!! April/May, 2007                                                     ! 
!!======================================================Marta Rodrigues=
!!                                                                     !
!! This subroutine is from ROMS (where is called ana_specir.h):        !
!! Revision: Nov/2015 - Updated to use nws=2                           !
!!                                                                     ! 
!!======================================================================
!! Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!                                                                     !
!=======================================================================
!                                                                      !
!  This subroutine sets surface solar downwelling spectral irradiance  !
!  at just beneath the sea surface, Ed(lambda,0-), in micromol quanta  !
!  per meter squared per second.                                       !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Gregg, W.W. and K.L. Carder, 1990:  A simple spectral solar       !
!           irradiance model for cloudless maritime atmospheres,       !
!           Limmol. Oceanogr., 35(8), 1657-1675.                       !
!                                                                      !
!=======================================================================
!
      USE bio_param
      USE eclight
      USE biology
      USE schism_glbl, only : nea,pi,xlon_el,ylat_el,idry_e,&
                              windx1,windx2,windy1,windy2,&
                              pr1,pr2,&
                              airt1,airt2,&
                              shum1,shum2,&
                              elnode, i34 
!
      IMPLICIT NONE
      SAVE

!
!  Imported variable declarations.
!
!      integer, intent(in) :: ne
!
!      real(r8), intent(in) :: yday, hour
!      real(r8), intent(in) :: lonr(ne)
!      real(r8), intent(in) :: latr(ne)
!      real(r8), intent(in) :: cloud(nea)
!      real(r8), intent(in) :: Hair(nea)
!      real(r8), intent(in) :: Tair(nea)
!      real(r8), intent(in) :: Pair(nea)
!      real(r8), intent(in) :: Uwind(nea)
!      real(r8), intent(in) :: Vwind(nea)
!      real(r8), intent(inout) :: specir(nea,Nbands)
!      real(r8), intent(inout) :: avcos(nea,Nbands)
      real(r8), intent(in) :: wtratio

! Local variables
      real(r8) :: eT,eTH2O,eTice
      real(r8) :: u_wind1,u_wind2,v_wind1,v_wind2
      real(r8) :: p_1,p_2
      real(r8) :: t_air1,t_air2
      real(r8) :: s_hum1,s_hum2
!
!  Local constant declarations.

      real(r8) :: am = 1.0_r8        ! Aerosol type 1-10: ocean to land
      real(r8) :: betalam = 0.55_r8
      real(r8) :: p0 = 29.92_r8      ! Standard pressure (inches of Hg)
      real(r8) :: rex = -1.6364_r8
      real(r8) :: roair = 1200.0_r8  ! Density of air (g/m3)
      real(r8) :: rn = 1.341_r8      ! Index of refraction of pure seawater
      real(r8) :: vis = 15.0_r8      ! Visibility (km)
      real(r8) :: wv = 1.5_r8        ! Precipitable water (cm): water vapor
!
!  Local variable declarations.
!
      
      integer :: i, iband, ic, j, nc
      integer :: iday !, month, year

      real(r8) :: Dangle, Hangle, LatRad, LonRad
      real(r8) :: cff, cff1, cff2  !, hour, yday
      real(r8) :: alpha, beta, gamma, theta, rtheta, rthetar
      real(r8) :: atra, gtra, otra, rtra, wtra
      real(r8) :: alg, arg, asymp, cosunz, Fa
      real(r8) :: frh, rh, rlam, rlogc
      real(r8) :: rm, rmin, rmo, rmp, rod, rof
      real(r8) :: ros, rospd, rosps, rpls 
      real(r8) :: sumx, sumx2, sumxy, sumy
      real(r8) :: taa, tas, to3, wa, wspeed, zenith

      real(r8), dimension(NBands) :: Fo, Edir, Edif, Ed, qlam

      real(r8), dimension(3) :: a_arr, dndr
      real(r8), dimension(3) :: ro     = (/ 0.03_r8, 0.24_r8, 2.00_r8 /)
      real(r8), dimension(3) :: r_arr  = (/ 0.10_r8, 1.00_r8, 10.0_r8 /)
!          

!
!-----------------------------------------------------------------------
!  Compute spectral irradiance: Using RADTRAN formulations.
!-----------------------------------------------------------------------
!
!  Assume time is in modified Julian day.  Get hour and year day.
!
      
      !!!!MFR - Change this... only hour and year day are needed
!            CALL caldate (r_date, tdays(ng), year, yday, month, iday, hour)
      !!!!!

! 
!
!  Estimate solar declination angle (radians).
!
      Dangle=23.44_r8*COS((172.0_r8-yday)*2.0_r8*pi/365.25_r8)
      Dangle=Dangle*deg2rad
!
!  Compute hour angle (radians).
!
     Hangle=(12.0_r8-hour)*pi/12.0_r8
!
!  Conversion constant from E to micromol quanta.
!   1/(Plank*c*(avos#*1.0e6).
!
      cff=1.0E-9_r8/(6.6256E-34_r8*2.998E8_r8*6.023E17_r8)
      DO iband=1,NBands
        qlam(iband)=ec_wave_ab(iband)*cff
      END DO
!
!  Correct solar constant for Earth-Sun distance.
!
      cff=(1.0_r8+0.0167_r8*COS(2.0_r8*pi*(yday-3.0_r8)/365.0_r8))**2
      DO iband=1,NBands
        Fo(iband)=ec_Fobar(iband)*cff
      END DO
!
!  Compute spectral irradiance.
!

! Marta Rodrigues
! For each element
      DO i=1,nea
        if(idry_e(i)==1) cycle

	  LatRad=ylat_el(i)*deg2rad
          LonRad=xlon_el(i)*deg2rad

! MFR - Nov/2015 - Compute atmos variables from nws=2 and cloud cover
!                  @ elements

            u_wind1=sum(windx1(elnode(1:i34(i),i)))/i34(i)
            u_wind2=sum(windx2(elnode(1:i34(i),i)))/i34(i)
            v_wind1=sum(windy1(elnode(1:i34(i),i)))/i34(i)
            v_wind2=sum(windy2(elnode(1:i34(i),i)))/i34(i)
            p_1=sum(pr1(elnode(1:i34(i),i)))/i34(i)
            p_2=sum(pr2(elnode(1:i34(i),i)))/i34(i)
            t_air1=sum(airt1(elnode(1:i34(i),i)))/i34(i)
            t_air2=sum(airt2(elnode(1:i34(i),i)))/i34(i)
            s_hum1=sum(shum1(elnode(1:i34(i),i)))/i34(i)
            s_hum2=sum(shum2(elnode(1:i34(i),i)))/i34(i)

            uwind(i)=u_wind1+wtratio*(u_wind2-u_wind1)
            vwind(i)=v_wind1+wtratio*(v_wind2-v_wind1)
            pair(i)=p_1+wtratio*(p_2-p_1)
            tair(i)=t_air1+wtratio*(t_air2-t_air1)
            hair(i)=s_hum1+wtratio*(s_hum2-s_hum1)

            ! Compute relative humidity
            ! 1st compute water and ice saturation vapor pressure, in hPa
            ! (CIMO Guide,WMO,2008)
            eTH2O = 6.112*exp(17.62*tair(i)/(243.12+tair(i)))

            eTice = 6.112*exp(22.46*tair(i)/(272.62+tair(i)))

            eT = MIN(eTH2O,eTice)

            ! epsilon = Rdry/Rvap = 287.58/461.5 ~ 0.622
            ! HR = 100*w/ws; q~w; ws~0.622*eT/p; HR~q*p*100/(0.622*eT)
            ! pair in hPa (pair=pair*0.01)
            ! hair in %
            ! tair in ºC 
           
            pair(i)=pair(i)*0.01d0
            hair(i) = 100.d0*hair(i)*pair(i)*461.5d0/(287.58d0*eT)

            IF(flag_cloud==1)THEN
               ! Rough estimation of cloud cover from relative humidity (Sundqvist et al.,
               ! 1989)
               cloud(i)=1-dsqrt((1-hair(i)/100)/(1-0.75))
               IF(cloud(i)<0.d0) cloud(i)=0.0d0
            ELSE
               ! Cloud cover from file
               ! MFR - In implementation
            ENDIF
! -------------------------------------------------------------------------------


!
!  Compute Climatological Ozone.
!
          to3=(235.0_r8+(150.0_r8+40.0_r8*                              &
     &         SIN(0.9865_r8*(yday-30.0_r8)*deg2rad)+                   &
     &         20.0_r8*SIN(3.0_r8*LonRad))*                             &
     &         SIN(1.28_r8*LatRad)*SIN(1.28_r8*LatRad))*                &
     &        0.001_r8                                 ! sco3 conversion

!
!  Local daylight is a function of the declination (Dangle) and hour 
!  angle adjusted for the local meridian (Hangle-lonr(i,j)*deg2rad). 
          
          cosunz=SIN(LatRad)*SIN(Dangle)+                               &
     &           COS(LatRad)*COS(Dangle)*COS(Hangle-xlon_el(i)*deg2rad)
          zenith=ACOS(cosunz)
          theta=zenith*rad2deg
!
!  Model for atmospheric transmittance of solar irradiance through
!  a maritime atmosphere.  Computes direct and diffuse separately.
!  Includes water vapor and oxygen absorption.
!
!  Compute atmospheric path lengths (air mass); pressure-corrected
!
          IF ((theta.ge.0.0_r8).and.(theta.le.90.0_r8)) THEN
!
!  Modified March, 1994 according to Kasten and Young 1989.
!
            rm=1.0_r8/(cosunz+0.50572_r8*(96.07995_r8-theta)**rex)
            rmp=rm*(pair(i)*0.02952756_r8)/p0
            rmo=(1.0_r8+22.0_r8/6370.0_r8)/                             &
     &          SQRT(cosunz*cosunz+44.0_r8/6370.0_r8)

!
!  Computes aerosol parameters according to a simplified version
!  of the Navy marine aerosol model.
!
!  Compute wind speed (24 hour mean is equal to current wind).
!
            wspeed=SQRT(uwind(i)*uwind(i)+vwind(i)*vwind(i))
!
!  Relative humidity factor, frh.
!
            rh=hair(i)
            IF (rh.ge.100.0_r8) rh=99.9_r8
            frh=((2.0_r8-rh*0.01_r8)/                                   &
                 (6.0_r8*(1.0_r8-rh*0.01_r8)))**0.333_r8
!
!  Size distribution amplitude components.
!
            a_arr(1)=2000.0_r8*am*am
            a_arr(2)=5.866_r8*(wspeed-2.2_r8)
            a_arr(3)=0.01527_r8*(wspeed-2.2_r8)*0.05_r8   !from Hughes 1987
            IF (a_arr(2).lt.0.5_r8) a_arr(2)=0.5_r8
            IF (a_arr(3).lt.0.000014_r8) a_arr(3)=0.000014_r8
!
!  Compute size distribution at three selected radii according to
!  Navy method.
!
            cff=1.0_r8/frh
            DO nc=1,3
              dndr(nc)=0.0_r8
              DO ic=1,3
                arg=LOG(r_arr(nc)/(frh*ro(ic)))
                dndr(nc)=dndr(nc)+a_arr(ic)*EXP(-arg*arg)*cff
              END DO
            END DO
!
!  Least squares approximation
!
            sumx=0.0_r8
            sumy=0.0_r8
            sumxy=0.0_r8
            sumx2=0.0_r8
            DO ic=1,3
              cff1=LOG10(r_arr(ic))
              cff2=LOG10(dndr(ic))
              sumx=sumx+cff1 
              sumy=sumy+cff2 
              sumxy=sumxy+cff1*cff2
              sumx2=sumx2+cff1*cff1
            END DO
            gamma=sumxy/sumx2
            rlogc=sumy/3.0_r8-gamma*sumx/3.0_r8          ! no used
            alpha=-(gamma+3.0_r8)
!
!  Compute aerosol turbity coefficient, beta.
!
            beta=(3.91_r8/vis)*betalam**alpha
!
!  Compute asymmetry parameter -- a function of alpha.
!
            IF (alpha.gt.1.2_r8) THEN
              asymp=0.65_r8
            ELSE IF (alpha .lt. 0.0_r8) THEN
              asymp=0.82_r8
            ELSE
              asymp=-0.14167_r8*alpha+0.82_r8
            END IF
!
!  Single scattering albedo at 550; function of RH.
!
            wa=(-0.0032_r8*am+0.972_r8)*EXP(0.000306_r8*rh)
!
!  Forward scattering probability.
!
            alg=LOG(1.0_r8-asymp)
            Fa=1.0_r8-0.5_r8*                                           &
     &         EXP((alg*(1.459_r8+alg*(0.1595_r8+alg*0.4129_r8))+       &
     &              alg*(0.0783_r8+alg*(-0.3824_r8-alg*0.5874_r8))*     &
     &             cosunz)*cosunz)
!
!  Compute surface reflectance for direct (rod) and diffuse (ros)
!  components separately, as a function of theta, wind speed or
!  stress.
!
!  Foam and diffuse reflectance
!
            IF (wspeed.gt.4.0_r8) THEN
              IF (wspeed.le.7.0_r8) THEN
                rof=roair*(0.00062_r8+0.00156_r8/wspeed)*               &
     &              0.000022_r8*wspeed*wspeed-0.00040_r8
              ELSE
                rof=(roair*(0.00049_r8+0.000065_r8*wspeed)*             &
     &               0.000045_r8-0.000040_r8)*wspeed*wspeed
              END IF
              rosps=0.057_r8
            ELSE
              rof=0.0_r8
              rosps=0.066_r8
            END IF
!
!  Direct Fresnel reflectance for theta < 40, wspeed < 2 m/s.
!
            IF ((theta.lt.40.0_r8).or.(wspeed.lt.2.0_r8)) THEN
              IF (theta.eq.0.0_r8) THEN
                rospd=0.0211_r8
              ELSE
                rtheta=zenith
                rthetar=ASIN(SIN(rtheta)/rn)
                rmin=rtheta-rthetar
                rpls=rtheta+rthetar
                rospd=0.5_r8*((SIN(rmin)*SIN(rmin))/                    &
     &                        (SIN(rpls)*SIN(rpls))+                    &
     &                        (TAN(rmin)*TAN(rmin))/                    &
     &                        (TAN(rpls)*TAN(rpls)))
              END IF
!
!  Empirical fit otherwise.
!
            ELSE
              rospd=0.0253_r8*EXP((-0.000714_r8*wspeed+0.0618_r8)*      &
     &                            (theta-40.0_r8))
            END IF
!
!  Reflectance totals.
!
            rod=rospd+rof
            ros=rosps+rof
!
!  Compute spectral irradiance for selected spectral bands.
!
            DO iband=1,NBands
              rlam=ec_wave_ab(iband)*0.001_r8
!
!  Transmittance, Rayleigh, by the method of Bird.
!
              rtra=EXP(-rmp/(115.6406_r8*rlam**4-1.335_r8*rlam**2))
!
!  Ozone.
!
              otra=EXP(-ec_aoz(iband)*to3*rmo)

!
!  Aerosols.
!
              arg=beta*rm*rlam**(-alpha)
              atra=EXP(-arg)
              taa=EXP(-(1.0_r8-wa)*arg)
              tas=EXP(-wa*arg)
!
!  Oxygen/gases.
!
              gtra=EXP((-1.41_r8*ec_ag(iband)*rmp)/                     &
     &                 ((1.0_r8+118.3_r8*ec_ag(iband)*rmp)**0.45_r8))

!
!  Water Vapor.
!
              wtra=EXP((-0.2385_r8*ec_aw(iband)*wv*rm)/                 &
     &                 ((1.0_r8+20.07_r8*ec_aw(iband)*wv*rm)**0.45_r8))

!
!  Direct irradiance.
!
              Edir(iband)=Fo(iband)*cosunz*rtra*otra*atra*gtra*         &
     &                    wtra*(1.0_r8-rod)

!
!  Total diffuse irradiance.
!
              Edif(iband)=(1.0_r8-ros)*                                 &
     &                    Fo(iband)*cosunz*gtra*wtra*otra*              &
     &                    (taa*0.5_r8*(1.0_r8-rtra**0.95_r8)+           &
     &                     taa*Fa*(1.0_r8-tas)*rtra**1.5_r8)

!
!  Cloud effects approximations, Kasten and Czeplak (1980).
!  (See Hydrolight Technical Notes).
!
              IF (cloud(i).gt.0.25_r8) THEN
                Ed(iband)=(Edir(iband)+Edif(iband))*                    &
     &                    (1.0_r8-0.75_r8*cloud(i)**3.4_r8)
                Edif(iband)=Ed(iband)*                                  &
     &                      (0.3_r8+0.7_r8*cloud(i)**2.0_r8)
              ELSE
                Ed(iband)=Edir(iband)+Edif(iband)
              END IF
!
!  Convert from W/cm/um to micromole quanta/m2/s.
!
              specir(i,iband)=Ed(iband)*10.0_r8*qlam(iband)
!
!  Correction of zenith angle after crossing air/sea interface.
!
              cff1=COS(ASIN((SIN(zenith))/rn))
              cff2=Edif(iband)/Ed(iband)
	      
              avcos(i,iband)=cff1*(1.0_r8-cff2)+0.86_r8*cff2
            END DO
          ELSE
            DO iband=1,NBands
              specir(i,iband)=0.0_r8
              avcos(i,iband)=0.66564_r8
            END DO
          END IF
      END DO
   
      RETURN
      END SUBROUTINE spec_ir
      
      
