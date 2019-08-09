#include "wwm_functions.h"
#ifdef SCHISM
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE STOKES_STRESS_INTEGRAL_SCHISM
        USE DATAPOOL
        use schism_glbl, only: iplg,errmsg,hmin_radstress,kbp
        USE schism_msgp
        implicit none
        integer     :: IP, k, ID, IS, IL
        real(rkind) :: eF1, eF2, eDelta, TheInt, eDep, eHeight
        real(rkind) :: eFrac, eFracB, eQuot
        real(rkind) :: eQuot1, eScal
        real(rkind) :: eOmega, eMult, kD, eSinc
        real(rkind) :: USTOKESpart, VSTOKESpart, eJPress
        real(rkind) :: eWk, eSigma, eLoc, eSinhkd, eSinh2kd, eSinhkd2
        real(rkind) :: PPTAIL, CETAIL, CKTAIL
        logical     :: DoTail
        real(rkind) :: eWkReal
        real(rkind) :: SumHeight
        real(rkind) :: eJPress_loc, eProd, eUint, eVint
        real(rkind) :: eUSTOKES_loc(NVRT), eVSTOKES_loc(NVRT)
        real(rkind) :: ZZETA
        REAL(rkind) :: tmp,TMP_X,TMP_Y,SINT,COST
        REAL(rkind) :: HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2
        REAL(rkind) :: ACLOC(MSC,MDC) 
 
        DO IP=1,MNP
          tmp=max(DEP8(IP)+ETA2(IP),hmin_radstress)
          eDep=tmp
!!!          eDep=SHYFZETA(NLEV(IP),MNP)
          eUSTOKES_loc=0
          eVSTOKES_loc=0
          eJpress_loc=0
          DO IS=1,MSC
            eMult=SPSIG(IS)*DDIR*DS_INCR(IS)
            eWk=WK(IS,IP)
            kD=MIN(KDMAX, eWk*eDep)
            eWkReal=kD/eDep
            eSinh2kd=DSINH(2*kD)
            eSinhkd=DSINH(kD)
            eSinhkd2=eSinhkd**2
            eSigma=SPSIG(IS)
            eUint=0
            eVint=0
            DO ID=1,MDC
              eLoc=AC2(IS,ID,IP)*eMult
              eJPress=G9*(kD/eSinh2kd)*(1/eDep) * eLoc
              eJPress_loc=eJPress_loc + eJPress
              eUint=eUint + eLoc*COSTH(ID)
              eVint=eVint + eLoc*SINTH(ID)
            END DO
            DO IL = KBP(IP), NVRT
              ZZETA = ZETA(IL,IP)-ZETA(KBP(IP),IP) !from bottom; 'z+D'
              eFrac=ZZETA/eDep
! Need some better understanding of vertical levels in SCHISM  
! for putting those correction terms.
!              eHeight=z_w_loc(k)-z_w_loc(k-1)
!              eFracB=eHeight/eDep
!              eSinc=SINH(kD*eFracB)/(kD*eFracB)
!              eQuot1=eSinc*DCOSH(2*kD*eFrac)/eSinhkd2
              eQuot1=DCOSH(2*kD*eFrac)/eSinhkd2
              eProd=eSigma*eWkReal*eQuot1
!YJZ: error
              eUSTOKES_loc(IL)=eUSTOKES_loc(IL) + eUint*eProd
              eVSTOKES_loc(IL)=eVSTOKES_loc(IL) + eVint*eProd
            ENDDO
          END DO
          !STOKES_X(:,IP)=eUSTOKES_loc
          !STOKES_Y(:,IP)=eVSTOKES_loc
          STOKES_VEL(1,:,IP)=eUSTOKES_loc
          STOKES_VEL(2,:,IP)=eVSTOKES_loc
          JPRESS(IP)=eJPress_loc
        ENDDO ! IP

      END SUBROUTINE STOKES_STRESS_INTEGRAL_SCHISM
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE RADIATION_STRESS_SCHISM
        USE DATAPOOL
        use schism_glbl, only: iplg,errmsg,hmin_radstress,kbp,wwave_force,idry
        USE schism_msgp !, only : myrank,parallel_abort
        IMPLICIT NONE

        INTEGER :: IP,IL,IS,ID
        REAL(rkind)  :: ACLOC(MSC,MDC)
        REAL(rkind)  :: COSE2, SINE2, COSI2
        REAL(rkind)  :: EWK(MNP),EWS(MNP),EWN(MNP),ETOT(MNP),MDIR(MNP)
        REAL(rkind)  :: m0, m0d, tmp, EHFR, ELOC, EFTAIL, ZZETA, DVEC2RAD
        REAL(rkind)  :: DS, D, KW, KD, SINH2KD, SINHKW, COSH2KW, COSHKW, COSHKD, ETOTS, ETOTC, EWSIG, S11, S22
        REAL(rkind)  :: WNTMP,WKTMP,WCGTMP,WCTMP,WN,WKDEPTMP
        REAL(rkind)  :: WSTMP, DEPLOC

        INTEGER      :: ND1,ND2
        REAL(rkind)  :: SINHKD,FSS(NVRT,MNP),FCS(NVRT,MNP),FSC(NVRT,MNP),FCC(NVRT,MNP)
        REAL(rkind)  :: dr_dxy(2,NVRT,nsa),HTOT,SXX3D0(NVRT,MNP),SYY3D0(NVRT,MNP),SXY3D0(NVRT,MNP)
        REAL(rkind)  :: WILD1(NVRT,MNP),WILD2(NVRT,MNP),WILD3(2,NVRT,nsa),WILD4(3,NVRT,MNP),DSPX,DSPY, WILD5(10)

!GD: imet_dry allows to choose between 2 different methods to compute the
!derivative at the sides between wet and dry elements:
!
! imet_dry=1 : only the values at the 2 nodes of the side are used to
! compute the derivative (this older method showed to provide inconsistent
! wave force at the wet/dry interface).
!
! imet_dry=2 : a 4-point stencil (the 3 wet nodes and an artificial
! node at the center of the side) is used to compute the derivative.
! This method is similar to using shape functions to compute the
! derivative at the center of the element and assigning this value to the
! the side center.
     
        IMET_DRY = 2 

        SXX3D(:,:) = ZERO
        SYY3D(:,:) = ZERO
        SXY3D(:,:) = ZERO

        EFTAIL = ONE / (PTAIL(1) - ONE)

        ETOT = ZERO
        MDIR = ZERO

        IF (LETOT) THEN
!AR: Estimate zeroth moment m0, mean wave direction, dominant wave number, dominant sigma ...
          DO IP = 1, MNP
!            IF (ABS(IOBP(IP)) .GT. 0) CYCLE
            !IF (DEP(IP) .LT. DMIN) CYCLE
            IF (IDRY(IP)==1) CYCLE

            DEPLOC = DEP(IP)
            ACLOC = AC2(:,:,IP)
            m0    = ZERO
            EWSIG  = ZERO
            ETOTS  = ZERO
            ETOTC  = ZERO
            IF (MSC .GE. 2) THEN
              DO ID = 1, MDC
                m0d = ZERO
                DO IS = 2, MSC
                  tmp = 0.5_rkind*(SPSIG(IS)*ACLOC(IS,ID)+SPSIG(IS-1)*ACLOC(IS-1,ID))*DS_INCR(IS)*DDIR
                  m0 = m0 + tmp
                  EWSIG  = EWSIG  + SPSIG(IS) * tmp
                  m0d = m0d + tmp
                END DO
                IF (MSC > 3) THEN
                  EHFR = ACLOC(MSC,ID) * SPSIG(MSC)
                  m0 = m0 + DDIR * EHFR * SPSIG(MSC) * EFTAIL
                endif
                ETOTC  = ETOTC + m0d * COS(SPDIR(ID))
                ETOTS  = ETOTS + m0d * SIN(SPDIR(ID))
              END DO
            ELSE
              DS = SGHIGH - SGLOW
              DO ID = 1, MDC
                m0d = ACLOC(1,ID) * DS * DDIR
                m0 = m0 + m0d
              END DO
            END IF
            ETOT(IP) = m0
            IF (m0 .GT. small .and. .not. dep(ip) .lt. dmin) then
              EWS(IP) = EWSIG/m0
              WSTMP = EWSIG/m0
              CALL ALL_FROM_TABLE(WSTMP,DEPLOC,WKTMP,WCGTMP,WKDEPTMP,WNTMP,WCTMP)
              EWN(IP) = WNTMP
              EWK(IP) = WKTMP
              MDIR(IP) = DVEC2RAD (ETOTC, ETOTS)
            ELSE
              EWS(IP)  = ZERO 
              EWN(IP)  = ZERO 
              EWK(IP)  = 10. 
              MDIR(IP) = ZERO 
            END IF 
          END DO !IP
        END IF !LETOT

!AR: Here comes the whole story ... 
! Etot = 1/16 * Hs² = 1/8 * Hmono² => Hs² = 2 * Hmono² => Hs = sqrt(2) * Hmono => Hmono = Hs / SQRT(2) 
! Etot = 1/16 * Hs² = 1/16 * (4 * sqrt(m0))² = m0 
! Etot = 1/8 * Hmono² ... so the problem for the analytical solution evolved because we treat the Etot from Hs and Hmono there is a factor of 2 between this!
! Or in other words for the analytical solution we impose a Hs = X[m], we integrate m0 out of it and get Etot, since this Etot is a function of Hs and not Hmono^X^O
! it needs the factor of 2 between it! This should make now things clear forever. So the question is not how we calculate the total energy the question is 
! what is defined on the boundary that means we should always recalculate the boundary in terms of Hs =  SQRT(2) * Hmono !!!
! Or saying it again in other words our boundary conditions is wrong if we impose Hmono in wwminput.nml !!!

        SXX3D = ZERO
        SXY3D = ZERO
        SYY3D = ZERO
        WWAVE_FORCE=ZERO
        IF (RADFLAG .EQ. 'LON') THEN
          RSXX = ZERO
          RSXY = ZERO
          RSYY = ZERO
          DO IP = 1, MNP
!            IF (ABS(IOBP(IP)) .GT. 0) CYCLE
            !IF (DEP(IP) .LT. DMIN) CYCLE
            IF (IDRY(IP)==1) CYCLE

            IF (.NOT. LETOT) THEN
              ACLOC = AC2(:,:,IP)
              DO ID = 1, MDC
                DO IS = 2, MSC
                  ELOC  = 0.5_rkind * (SPSIG(IS)*ACLOC(IS,ID)+SPSIG(IS-1)*ACLOC(IS-1,ID))*DS_INCR(IS)*DDIR
                  COSE2 = COS(SPDIR(ID))**TWO
                  SINE2 = SIN(SPDIR(ID))**TWO
                  COSI2 = COS(SPDIR(ID)) * SIN(SPDIR(ID))
                  WN    = CG(IS,IP) / ( SPSIG(IS)/WK(IS,IP) )
                  RSXX(IP) = RSXX(IP) + ( WN * COSE2 + WN - 0.5_rkind) * ELOC   ! Units = [ 1/s + 1/s - 1/s ] * m²s = m²
                  RSXY(IP) = RSXY(IP) + ( WN * COSI2          ) * ELOC
                  RSYY(IP) = RSYY(IP) + ( WN * SINE2 + WN - 0.5_rkind) * ELOC
                ENDDO
              ENDDO
            ELSE IF (LETOT) THEN
              RSXX(IP) =  ETOT(IP) * (EWN(IP)*((EWK(IP)*SIN(MDIR(IP)))**TWO/EWK(IP)**TWO+ONE)-0.5_rkind)
              RSXY(IP) =  ETOT(IP) *  EWN(IP)* EWK(IP)*SIN(MDIR(IP))*EWK(IP)*COS(MDIR(IP))* ONE/EWK(IP)
              RSYY(IP) =  ETOT(IP) * (EWN(IP)*((EWK(IP)*COS(MDIR(IP)))**TWO/EWK(IP)**TWO+ONE)-0.5_rkind)
            END IF 
          END DO

          DO IP = 1, MNP
!            IF (ABS(IOBP(IP)) .GT. 0) CYCLE
            !IF (DEP(IP) .LT. DMIN)  THEN
            IF (IDRY(IP)==1) THEN
              SXX3D(:,IP) = ZERO
              SXY3D(:,IP) = ZERO
              SYY3D(:,IP) = ZERO
            ELSE
              SXX3D(:,IP) = RSXX(IP) / DEP(IP) * G9
              SXY3D(:,IP) = RSXY(IP) / DEP(IP) * G9
              SYY3D(:,IP) = RSYY(IP) / DEP(IP) * G9
            END IF
          END DO
          !Store as double for force later
          SXX3D0 = SXX3D
          SXY3D0 = SXY3D
          SYY3D0 = SYY3D

        ELSE
          call parallel_abort('R.S.: unknown R.S. model') 
        END IF !RADFLAG 

!       Integrate over depth for checking
        RSXX = ZERO
        DO IP = 1, MNP
          IF(IDRY(IP)==1) CYCLE
          DO IL = KBP(IP)+1, NVRT 
            RSXX(IP) = RSXX(IP) + 0.5_rkind*( SXX3D(IL,IP)+SXX3D(IL-1,IP) ) * ABS((ZETA(IL,IP) - ZETA(IL-1,IP)))/G9
          END DO !IL
        END DO !IP
!'
!       Computation in double precision here
!       SXX3D0() etc. should have dimension of m^2/s/s, defined at nodes and whole levels.
!       Use same arrays to temporarily store properly scaled Sxx etc
!       write(12,*)'Checking Sxx,Sxy,Syy:'
        do IP=1,MNP
          IF(IDRY(IP)==1) then
            SXX3D0(:,IP)=ZERO
            SYY3D0(:,IP)=ZERO
            SXY3D0(:,IP)=ZERO
            cycle
          endif

          do IL=KBP(IP),NVRT
!           D*(Sxx, Sxy, Syy)/rho in Xia et al. (2004) 
!           After this the dimension should be m^3/s/s
            tmp=max(DEP8(IP)+ETA2(IP),hmin_radstress)
            !SXX3D0(IL,IP)=(DEP8(IP)+ETA2(IP))*SXX3D0(IL,IP) !D*Sxx/rho
            !SXY3D0(IL,IP)=(DEP8(IP)+ETA2(IP))*SXY3D0(IL,IP) !D*Sxy/rho
            !SYY3D0(IL,IP)=(DEP8(IP)+ETA2(IP))*SYY3D0(IL,IP) !D*Syy/rho
            SXX3D0(IL,IP)=tmp*SXX3D0(IL,IP) !D*Sxx/rho
            SXY3D0(IL,IP)=tmp*SXY3D0(IL,IP) !D*Sxy/rho
            SYY3D0(IL,IP)=tmp*SYY3D0(IL,IP) !D*Syy/rho
          enddo !k
        enddo !IP

!       Compute radiation stress force 
!       wwave_force(:,1:nsa,1:2) = Rsx, Rsy in my notes (the terms in momen. eq.)
!       and has a dimension of m/s/s
        call hgrad_nodes(IMET_DRY,0,NVRT,MNP,nsa,SXX3D0,dr_dxy)
        call exchange_s3d_2(dr_dxy)

        do IS=1,nsa
          if(idry_s(IS)==0) then
            HTOT=(eta2(isidenode(1,IS))+eta2(isidenode(2,IS))+DEP8(isidenode(1,IS))+DEP8(isidenode(2,IS)))/2
!            write(*,*) sum(eta2), sum(dep8)
!            write(*,*) HTOT, eta2(isidenode(1,IS)), eta2(isidenode(2,IS)), DEP8(isidenode(1,IS)), DEP8(isidenode(1,IS)), isidenode(1,IS), isidenode(1,IS) 
!            if (isidenode(1,IS) == 150 .AND. isidenode(2,IS) == 149) stop
            if(HTOT<=0) call parallel_abort('RADIATION_STRESS: (999)')
            HTOT=max(HTOT,hmin_radstress)
            do il = 1, nvrt
              WWAVE_FORCE(1,il,IS)=WWAVE_FORCE(1,il,IS)-dr_dxy(1,il,IS)/HTOT
            end do
          endif
        enddo !IS

        call hgrad_nodes(IMET_DRY,0,NVRT,MNP,nsa,SYY3D0,dr_dxy)
        call exchange_s3d_2(dr_dxy)
        do IS=1,nsa
          if(idry_s(IS)==0) then
            HTOT=(eta2(isidenode(1,IS))+eta2(isidenode(2,IS))+DEP8(isidenode(1,IS))+DEP8(isidenode(2,IS)))/2
            if(HTOT<=0) call parallel_abort('RADIATION_STRESS: (998)')
            HTOT=max(HTOT,hmin_radstress)
            do il = 1, nvrt 
              WWAVE_FORCE(2,il,IS)=WWAVE_FORCE(2,il,IS)-dr_dxy(2,il,IS)/HTOT
            end do
          endif
        enddo !IS

        call hgrad_nodes(IMET_DRY,0,NVRT,MNP,nsa,SXY3D0,dr_dxy)
        call exchange_s3d_2(dr_dxy)
!
!          write(12,*)'Checking R.S.'
!
        do IS=1,nsa
          if(idry_s(IS)==0) then
            HTOT=(eta2(isidenode(1,IS))+eta2(isidenode(2,IS))+DEP8(isidenode(1,IS))+DEP8(isidenode(2,IS)))/2
            if(HTOT<=0) call parallel_abort('RADIATION_STRESS: (997)')
            HTOT=max(HTOT,hmin_radstress)
            WWAVE_FORCE(1,:,IS)=WWAVE_FORCE(1,:,IS)-dr_dxy(2,:,IS)/HTOT
            WWAVE_FORCE(2,:,IS)=WWAVE_FORCE(2,:,IS)-dr_dxy(1,:,IS)/HTOT
          endif
        enddo !IS

      END SUBROUTINE RADIATION_STRESS_SCHISM
#endif /*SCHISM*/
!**********************************************************************
!*                                                                    *
!**********************************************************************
