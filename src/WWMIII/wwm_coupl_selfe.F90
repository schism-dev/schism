!Note: most arrays in this file are from SCHISM directly (too account for
!quads)
#include "wwm_functions.h"
#ifdef SCHISM
!**********************************************************************
!*   This routine is for RADFLAG=LON (Longuet-Higgins)
!**********************************************************************
      SUBROUTINE RADIATION_STRESS_SCHISM
        USE DATAPOOL
        use schism_glbl, only: errmsg,hmin_radstress,npa,nsa,idry_s,isidenode 
        USE schism_msgp 
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
            IF (idry(IP)==1) CYCLE

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
            IF (idry(IP)==1) CYCLE

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
            IF (idry(IP)==1) THEN
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
          IF(idry(IP)==1) CYCLE
          DO IL = KBP(IP)+1, NVRT 
            RSXX(IP) = RSXX(IP) + 0.5_rkind*( SXX3D(IL,IP)+SXX3D(IL-1,IP) ) * ABS((ZETA(IL,IP) - ZETA(IL-1,IP)))/G9
          END DO !IL
        END DO !IP
!'
!       Computation in double precision here
!       SXX3D0() etc. should have dimension of m^3/s/s, defined at nodes and whole levels.
!       Use same arrays to temporarily store properly scaled Sxx etc
!       write(12,*)'Checking Sxx,Sxy,Syy:'
        do IP=1,MNP
          IF(idry(IP)==1) then
            SXX3D0(:,IP)=ZERO
            SYY3D0(:,IP)=ZERO
            SXY3D0(:,IP)=ZERO
            cycle
          endif

          do IL=KBP(IP),NVRT
!           D*(Sxx, Sxy, Syy)/rho in Xia et al. (2004) 
            tmp=max(DEP8(IP)+ETA2(IP),hmin_radstress)
            SXX3D0(IL,IP)=tmp*SXX3D0(IL,IP) !D*Sxx/rho
            SXY3D0(IL,IP)=tmp*SXY3D0(IL,IP) !D*Sxy/rho
            SYY3D0(IL,IP)=tmp*SYY3D0(IL,IP) !D*Syy/rho
          enddo !k
        enddo !IP

!       Compute radiation stress force 
!       wwave_force(:,1:MNS,1:2) = Rsx, Rsy in my notes (the terms in momen. eq.)
!       and has a dimension of m/s/s
        call hgrad_nodes(IMET_DRY,0,NVRT,npa,nsa,SXX3D0,dr_dxy) ![dr_dxy]=m^2/s/s
        call exchange_s3d_2(dr_dxy)

        do IS=1,nsa
          if(idry_s(IS)==0) then
            HTOT=(eta2(isidenode(1,IS))+eta2(isidenode(2,IS))+DEP8(isidenode(1,IS))+DEP8(isidenode(2,IS)))/2
!            write(*,*) sum(eta2), sum(dep8)
!            write(*,*) HTOT, eta2(isidenode(1,IS)), eta2(isidenode(2,IS)), DEP8(isidenode(1,IS)), DEP8(isidenode(1,IS)), isidenode(1,IS), isidenode(1,IS) 
!            if (isidenode(1,IS) == 150 .AND. isidenode(2,IS) == 149) stop
            if(HTOT<=0) call parallel_abort('RADIATION_STRESS: (999)')
            HTOT=max(HTOT,hmin_radstress)
            do il = 1, NVRT
              WWAVE_FORCE(1,il,IS)=WWAVE_FORCE(1,il,IS)-dr_dxy(1,il,IS)/HTOT !m/s/s
            end do
          endif
        enddo !IS

        call hgrad_nodes(IMET_DRY,0,NVRT,npa,nsa,SYY3D0,dr_dxy)
        call exchange_s3d_2(dr_dxy)
        do IS=1,nsa
          if(idry_s(IS)==0) then
            HTOT=(eta2(isidenode(1,IS))+eta2(isidenode(2,IS))+DEP8(isidenode(1,IS))+DEP8(isidenode(2,IS)))/2
            if(HTOT<=0) call parallel_abort('RADIATION_STRESS: (998)')
            HTOT=max(HTOT,hmin_radstress)
            do il = 1, NVRT 
              WWAVE_FORCE(2,il,IS)=WWAVE_FORCE(2,il,IS)-dr_dxy(2,il,IS)/HTOT
            end do
          endif
        enddo !IS

        call hgrad_nodes(IMET_DRY,0,NVRT,npa,nsa,SXY3D0,dr_dxy)
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

!**********************************************************************
!*  This routine is used with RADFLAG=VOR (3D vortex formulation, after Bennis et al., 2011)
!*  => Computation of the wave-induced pressure term.		 
!**********************************************************************
      SUBROUTINE STOKES_STRESS_INTEGRAL_SCHISM
        USE DATAPOOL
        use schism_glbl, only: errmsg,hmin_radstress 
        USE schism_msgp
        implicit none
        integer     :: IP, k, ID, IS, IL
        real(rkind) :: eDep, eHeight, eFrac, eFracB
        real(rkind) :: eQuot1, eMult, kD, eSinc, eWkReal
        real(rkind) :: eWk, eSigma, eLoc, eSinhkd, eSinh2kd, eSinhkd2
        real(rkind) :: eJPress, eJPress_loc, eProd, eUint, eVint
        real(rkind) :: eUSTOKES_loc(NVRT), eVSTOKES_loc(NVRT)
        real(rkind) :: ZZETA
 
        DO IP = 1,MNP
          if(idry(IP) == 1) CYCLE 

          ! Total water depth at the node
          eDep = max(DEP8(IP)+ETA2(IP),hmin_radstress)
          ! Initialization of the local Stokes drift and J variables
          eUSTOKES_loc = 0.d0
          eVSTOKES_loc = 0.d0
          eJPress_loc  = 0.d0

          ! Loop on the frequencies
          DO IS = 1,MSC
            eMult = SPSIG(IS)*DDIR*DS_INCR(IS)
            eWk = WK(IS,IP)
            kD = MIN(KDMAX, eWk*eDep)
            eWkReal = kD/eDep
            eSinh2kd = DSINH(2.d0*kD)
            eSinhkd = DSINH(kD)
            eSinhkd2 = eSinhkd**2.d0
            eSigma = SPSIG(IS)
            eUint = 0.d0
            eVint = 0.d0

            ! Loop on the directions
            DO ID = 1,MDC
              eLoc = AC2(IS,ID,IP)*eMult
              eJPress = G9*(kD/eSinh2kd)*(1.d0/eDep)*eLoc
              eJPress_loc = eJPress_loc + eJPress
              eUint = eUint + eLoc*COSTH(ID)
              eVint = eVint + eLoc*SINTH(ID)
            ENDDO !MDC

            ! Loop on the vertical nodes
            DO IL = KBP(IP), NVRT
              ZZETA = ZETA(IL,IP)-ZETA(KBP(IP),IP) !from bottom; 'z+D'
              eFrac = ZZETA/eDep
! Need some better understanding of vertical levels in SCHISM  
! for putting those correction terms.
!              eHeight = z_w_loc(k)-z_w_loc(k-1)
!              eFracB = eHeight/eDep
!              eSinc = SINH(kD*eFracB)/(kD*eFracB)
!              eQuot1 = eSinc*DCOSH(2*kD*eFrac)/eSinhkd2
              eQuot1 = DCOSH(2.d0*kD*eFrac)/eSinhkd2
              eProd = eSigma*eWkReal*eQuot1
              eUSTOKES_loc(IL) = eUSTOKES_loc(IL) + eUint*eProd
              eVSTOKES_loc(IL) = eVSTOKES_loc(IL) + eVint*eProd
            ENDDO !NVRT
          ENDDO !MSC

          ! Storing Stokes drift and J term variables
          STOKES_VEL(1,:,IP) = eUSTOKES_loc
          STOKES_VEL(2,:,IP) = eVSTOKES_loc
          JPRESS(IP) = eJPress_loc
        ENDDO !MNP

      END SUBROUTINE STOKES_STRESS_INTEGRAL_SCHISM

!**********************************************************************
!*  This routine is used with RADFLAG=VOR (3D vortex formulation, after Bennis et al., 2011)
!*  => Computation of the conservative terms A1 and B1 from Eq. (11) and (12) respectively
!**********************************************************************
      SUBROUTINE COMPUTE_CONSERVATIVE_VF_TERMS_SCHISM
        USE DATAPOOL
        use schism_glbl, only: errmsg,hmin_radstress,kbs,kbe,npa,nsa,ns,ne,nea,idry_e,idry_s, & 
                             & nne,isidenode,isdel,indel,elnode,elside,i34,area,dldxy,dr_dxy, &
                             &cori,zs,su2,sv2
        USE schism_msgp
        implicit none

        integer     :: IP, ID, IS, IE, isd, k, j, l, n1, n2, n3, icount
        real(rkind) :: tmp0, tmp1, tmp2, ztmp, ubar, vbar, dhdx, dhdy
        real(rkind) :: stokes_w(nvrt,nea), stokes_w_sd(nvrt,nsa), ws_tmp1(nvrt,nsa), ws_tmp2(nvrt,nsa)

!YJZ: init here for better readability
        WWAVE_FORCE=0.d0
!...    Pressure term (grad(J))
        do IS = 1,ns !resident
          if(idry_s(IS) == 1) cycle

          ! Wet side
          icount = 0 
          tmp1 = 0; tmp2 = 0
          do l = 1,2 !elements
            IE = isdel(l,IS)
            if(ie /= 0) then
              if(idry_e(IE) == 0) then
                icount = icount + 1
                tmp1 = tmp1 + dot_product(JPRESS(elnode(1:i34(IE),IE)),dldxy(1:i34(IE),1,IE)) !in eframe
                tmp2 = tmp2 + dot_product(JPRESS(elnode(1:i34(IE),IE)),dldxy(1:i34(IE),2,IE))
              endif
            endif 
          enddo !l

          ! Averaging the values from the two surrounding elements
          if(icount > 2) call parallel_abort('Pressure term: icount>2')
          if(icount == 2) then
            tmp1 = tmp1/2.d0
            tmp2 = tmp2/2.d0
          endif

          ! Updating the wave forces
          do k = kbs(IS),NVRT
            WWAVE_FORCE(1,k,IS) = WWAVE_FORCE(1,k,IS) - tmp1
            WWAVE_FORCE(2,k,IS) = WWAVE_FORCE(2,k,IS) - tmp2
          enddo
        enddo !ns

        ! Exchange between ghost regions
        call exchange_s3d_2(WWAVE_FORCE) 

!...    Stokes vel. at side (in pframe if ics=2)
        ! The average of the values from vertically adjacent nodes is taken
        STOKES_VEL_SD=0.d0
        do IS = 1,nsa
          if(idry_s(IS) == 0) then
            n1 = isidenode(1,IS); n2 = isidenode(2,IS)
            do k = kbs(IS),NVRT
              STOKES_VEL_SD(1,k,IS) = (stokes_vel(1,k,n1) + stokes_vel(1,k,n2))/2.d0
              STOKES_VEL_SD(2,k,IS) = (stokes_vel(2,k,n1) + stokes_vel(2,k,n2))/2.d0
            enddo
          endif
        enddo !MNS
 
!...    Contribution from the terms with Coriolis force and the spatial derivative of u
        call hgrad_nodes(2,0,NVRT,npa,nsa,uu2,dr_dxy)
        do IS = 1,ns
          if(idry_s(IS) == 0) then
            do k = kbs(IS),NVRT
               ! f*v_s-du/dy*v_s
               WWAVE_FORCE(1,k,IS) = WWAVE_FORCE(1,k,IS) + (cori(IS) - dr_dxy(2,k,IS))*STOKES_VEL_SD(2,k,IS)
               ! -f*u_s+du/dy*u_s
               WWAVE_FORCE(2,k,IS) = WWAVE_FORCE(2,k,IS) + (-cori(IS) + dr_dxy(2,k,IS))*STOKES_VEL_SD(1,k,IS)
            enddo
          endif
        enddo !ns

!...    Contribution from the terms with spatial derivative of v
        call hgrad_nodes(2,0,NVRT,npa,nsa,vv2,dr_dxy)
        do IS = 1,ns
          if(idry_s(IS) == 0) then
            do k = kbs(IS),NVRT
               ! dv/dx*v_s
               WWAVE_FORCE(1,k,IS) = WWAVE_FORCE(1,k,IS) + dr_dxy(1,k,IS)*STOKES_VEL_SD(2,k,IS)
!     & (STOKES_VEL(2,k,isidenode(1,IS)) + STOKES_VEL(2,k,isidenode(2,IS)))/2.d0
               ! -dv/dx*u_s
               WWAVE_FORCE(2,k,IS) = WWAVE_FORCE(2,k,IS) - dr_dxy(1,k,IS)*STOKES_VEL_SD(1,k,IS) 
!     & (STOKES_VEL(1,k,isidenode(1,IS)) + STOKES_VEL(1,k,isidenode(2,IS)))/2.d0
            enddo
          endif
        enddo !ns
        ! Exchange between ghost regions
        call exchange_s3d_2(WWAVE_FORCE)

!...    Compute _bottom_ Stokes z-vel. at elements
        stokes_w = 0.d0
        do IE = 1,nea
           if(idry_e(IE) == 1) cycle
!           n1 = elnode(1,IE)
!           n2 = elnode(2,IE)
!           n3 = elnode(3,IE)
           if(kbe(IE) == 0) then
             write(errmsg,*)'Error: kbe(i) == 0'
             call parallel_abort(errmsg)
           endif

           ubar=0; vbar=0
           do j=1,i34(IE)
             n1=elnode(j,IE)
             ubar=ubar+STOKES_VEL(1,kbp(n1),n1)/i34(IE)
             vbar=vbar+STOKES_VEL(2,kbp(n1),n1)/i34(IE)
           enddo !j
           dhdx=dot_product(DEP8(elnode(1:i34(ie),ie)),dldxy(1:i34(ie),1,ie))
           dhdy=dot_product(DEP8(elnode(1:i34(ie),ie)),dldxy(1:i34(ie),2,ie))

!           ubar = (STOKES_VEL(1,max(kbp(n1),kbe(IE)),n1) + STOKES_VEL(1,max(kbp(n2),kbe(IE)),n2) &
!               & + STOKES_VEL(1,max(kbp(n3),kbe(IE)),n3))/3.d0 !average bottom stokes-x-vel
!           vbar = (STOKES_VEL(2,max(kbp(n1),kbe(IE)),n1) + STOKES_VEL(2,max(kbp(n2),kbe(IE)),n2) &
!               & + STOKES_VEL(2,max(kbp(n3),kbe(IE)),n3))/3.d0 !average bottom stokes-y-vel
!           dhdx = DEP8(n1)*dldxy(1,1,IE) + DEP8(n2)*dldxy(2,1,IE) + DEP8(n3)*dldxy(3,1,IE) !eframe
!           dhdy = DEP8(n1)*dldxy(1,2,IE) + DEP8(n2)*dldxy(2,2,IE) + DEP8(n3)*dldxy(3,2,IE)
           stokes_w(kbe(IE),IE) = -dhdx*ubar - dhdy*vbar
        enddo !nea

!...    Compute bottom Stokes z-vel. at nodes
        STOKES_W_ND = 0.d0
        do IP = 1,np !residents only
           if(idry(IP) == 1) cycle

           !Bottom Stokes z-vel.
           tmp0 = 0.d0
           do j = 1,nne(IP)
              ie = indel(j,IP)
              if(idry_e(ie)==0) then
                STOKES_W_ND(kbp(IP),IP) = STOKES_W_ND(kbp(IP),IP) + stokes_w(kbe(ie),ie)*area(ie)
              endif
              tmp0 = tmp0 + area(ie)
           enddo !j
           STOKES_W_ND(kbp(IP),IP) = STOKES_W_ND(kbp(IP),IP)/tmp0
        enddo !np

!...    Compute horizontal gradient of Stokes x and y-vel. (to compute Stokes z-vel.)
        ws_tmp1 = 0.d0; ws_tmp2 = 0.d0
        call hgrad_nodes(2,0,NVRT,npa,nsa,STOKES_VEL(1,:,:),dr_dxy)
        ws_tmp1(:,:) = dr_dxy(1,:,:)
        call hgrad_nodes(2,0,NVRT,npa,nsa,STOKES_VEL(2,:,:),dr_dxy)
        ws_tmp2(:,:) = dr_dxy(2,:,:)

!...    Compute Stokes z-vel. at all levels: stokes_w_sd(NVRT,nsa)
        stokes_w_sd = 0.d0
        do IS = 1,ns !residents only
          if(idry_s(IS) == 1) cycle
          n1 = isidenode(1,IS)
          n2 = isidenode(2,IS)

          !Bottom Stokes z-vel.
          stokes_w_sd(kbs(IS),IS) = (STOKES_W_ND(max(kbs(IS),kbp(n1)),n1)+STOKES_W_ND(max(kbs(IS),kbp(n2)),n2))/2.d0

          !Stokes z-vel. at all levels
          do k = kbs(IS)+1,NVRT
            ztmp = zs(k,IS)-zs(k-1,IS)
            stokes_w_sd(k,IS) = stokes_w_sd(k-1,IS)-(ws_tmp1(k,IS)+ws_tmp1(k-1,IS))/2.d0*ztmp &
     & -(ws_tmp2(k,IS)+ws_tmp2(k-1,IS))/2.d0*ztmp
           enddo
        enddo !ns

!...    Compute Stokes z-vel. at elements and all levels
!       (this does not appear to be used)
        do IE = 1,ne
          if(idry_e(IE) == 1) cycle

          !Wet elem
          do k = kbe(IE)+1,NVRT
            stokes_w(k,IE) = 0.d0
            do j = 1,i34(IE)
              isd = elside(j,IE) !wet
              stokes_w(k,IE) = stokes_w(k,IE) + stokes_w_sd(max(k,kbs(isd)),isd)/i34(IE)
            enddo !j
          enddo !k
        enddo !ne


!...    Adding term -w_s*(du/dz,dv/dz) to the wave forces
        do IS = 1,ns
          if(idry_s(IS) == 1) cycle
          tmp1 = 0.d0
          tmp2 = 0.d0
          do k = kbs(IS),NVRT
            if (k == kbs(IS)) then
              ztmp = zs(k+1,IS) - zs(k,IS)
              tmp1 = su2(k+1,IS) - su2(k,IS)
              tmp2 = sv2(k+1,IS) - sv2(k,IS)
            elseif (k == NVRT) then
              ztmp = zs(k,IS)-zs(k-1,IS)
              tmp1 = su2(k,IS) - su2(k-1,IS)
              tmp2 = sv2(k,IS) - sv2(k-1,IS)
            else
              ztmp = zs(k+1,IS) - zs(k-1,IS)
              tmp1 = su2(k+1,IS) - su2(k-1,IS)
              tmp2 = sv2(k+1,IS) - sv2(k-1,IS)
            endif

            if(ztmp==0) call parallel_abort('wwm_coupl_selfe: ztmp=0')
            WWAVE_FORCE(1,k,IS) = WWAVE_FORCE(1,k,IS) - stokes_w_sd(k,IS)*tmp1/ztmp
            WWAVE_FORCE(2,k,IS) = WWAVE_FORCE(2,k,IS) - stokes_w_sd(k,IS)*tmp2/ztmp

          enddo !k
        enddo !i=1,ns
        ! Exchange between ghost regions
        call exchange_s3d_2(WWAVE_FORCE)

      END SUBROUTINE COMPUTE_CONSERVATIVE_VF_TERMS_SCHISM

!**********************************************************************
!*  This routine is used with RADFLAG=VOR (3D vortex formulation, after Bennis et al., 2011)
!*  => Computation of the non-conservative terms due to depth-induced breaking (term Fb from Eq. (11) and (12))
!**********************************************************************
      SUBROUTINE COMPUTE_BREAKING_VF_TERMS_SCHISM
        USE DATAPOOL
        USE schism_glbl, only: errmsg,hmin_radstress,kbs,kbe,nsa,ns,np,ne,nea,idry_e, & 
                             & idry_s,isidenode,isbs,i34,dps,dldxy,dr_dxy,h0,out_wwm,zs 
        USE schism_msgp 
        IMPLICIT NONE

        integer     :: IS, isd, k, j, l, n1, n2, n3, icount
        real(rkind) :: eta_tmp, tmp0, htot, sum1
        real(rkind) :: swild(NVRT)

        ! Compute sink of momentum due to wave breaking 
        do IS = 1,nsa
          ! Check if dry segment or open bnd segment
          if(idry_s(IS) == 1 .or. isbs(IS) > 0) cycle
          
          ! Water depth at side
          n1 = isidenode(1,IS); n2 = isidenode(2,IS)
          eta_tmp = (eta2(n1) + eta2(n2))/2.d0
          htot = max(h0,dps(IS)+eta_tmp,hmin_radstress) ! KM
          !htot = max(h0,dps(IS)+eta_tmp)
   
          if(kbs(IS)+1 == NVRT) then !2D

            ! Breaking acceleration: average between the two adjacent nodes
            WWAVE_FORCE(1,:,IS) = WWAVE_FORCE(1,:,IS) - (SBR(1,n1) + SBR(1,n2))/2.d0/htot
            WWAVE_FORCE(2,:,IS) = WWAVE_FORCE(2,:,IS) - (SBR(2,n1) + SBR(2,n2))/2.d0/htot

          else !3D

            ! Threshold on Hs
            tmp0 = (out_wwm(n1,1) + out_wwm(n2,1))/2.d0 !Hs
            if(tmp0 <= 0.05d0) cycle
            if(tmp0/htot < 0.1d0) cycle

            ! Vertical distribution function of qdm (due to wave breaking)
            swild = 0.d0
            do k = kbs(IS),NVRT
              ! Homogeneous vertical distribution
!              swild(k)=1.d0

              ! All qdm injected in the surface layer
!              if (k == NVRT) then
!                swild(k)=1.d0
!              else
!                swild(k)=0.d0
!              endif

              ! Profile 1
              swild(k) = cosh((dps(IS)+zs(k,IS))/(0.6d0*tmp0))

              ! Profile 2
!              swild(k) = 1.d0 - dtanh(((eta_tmp-zs(k,j))/(0.2d0*tmp0))**2.d0)

              ! Profile 3
!              swild(k) = 1.d0 - dtanh(((eta_tmp-zs(k,j))/(0.2d0*tmp0))**4.d0)
            enddo !k
 
            ! Integral of the vertical distribution function
            sum1 = 0
            do k = kbs(IS),NVRT-1
              sum1 = sum1 + (swild(k+1) + swild(k))/2.d0*(zs(k+1,IS) - zs(k,IS))
            enddo !NVRT-1
            if(sum1 == 0) call parallel_abort('wave_break: integral=0')

            do k = kbs(IS),NVRT
              ! Breaking acceleration
              WWAVE_FORCE(1,k,IS) = WWAVE_FORCE(1,k,IS) - swild(k)/sum1*(SBR(1,n1) + SBR(1,n2))/2.d0
              WWAVE_FORCE(2,k,IS) = WWAVE_FORCE(2,k,IS) - swild(k)/sum1*(SBR(2,n1) + SBR(2,n2))/2.d0
            enddo
          endif !2D/3D
        enddo !nsa

      END SUBROUTINE COMPUTE_BREAKING_VF_TERMS_SCHISM

!**********************************************************************
!*   This routine fixes the wave forces to the barotropic gradient at the numerical shoreline (boundary between dry and wet elements)			
!**********************************************************************
      SUBROUTINE SHORELINE_WAVE_FORCES
        USE DATAPOOL
        USE schism_glbl, only: errmsg,hmin_radstress,idry_e,idry_s,isidenode, & 
                             & isdel,elnode,i34,dldxy,grav,thetai,ns 
        USE schism_msgp 
        IMPLICIT NONE

        INTEGER      :: IP,IS,INODE_1,INODE_2,IELEM
        REAL(rkind)  :: TMP_X,TMP_Y,BPGR(2)

!... 	Loop on the resident sides
        do IS = 1,ns
          if(idry_s(IS) == 1) cycle

          ! Adjacent nodes index
          INODE_1 = isidenode(1,IS) ! Side node #1
          INODE_2 = isidenode(2,IS) ! Side node #2
          if(idry(INODE_1) == 1 .OR. idry(INODE_2) == 1) cycle  

          ! Sides we are not interested in
          if(isdel(1,IS) == 0 .OR. isdel(2,IS) == 0) cycle ! Boundaries
          if(idry_e(isdel(1,IS)) == 0 .AND. idry_e(isdel(2,IS)) == 0) cycle ! Case where both adjacent elements are wet    
          if(idry_e(isdel(1,IS)) == 1 .AND. idry_e(isdel(2,IS)) == 1) cycle ! Case where both adjacent elements are dry (should never occur anyway)
          !if(isbnd(1,INODE_1) == 1 .OR. isbnd(1,INODE_2) == 1) cycle ! Case where the side touches open boundaries
          if(isbnd(1,INODE_1)>0.OR.isbnd(1,INODE_2)>0) cycle ! Case where the side touches open boundaries

          ! Reinitializing the wave force
          WWAVE_FORCE(:,:,IS) = 0

          ! We are left with sides that belong to one dry element and one wet element
          ! We store the elements indexes for future use
          if(idry_e(isdel(1,IS)) == 0 .AND. idry_e(isdel(2,IS)) == 1) then 
            IELEM = isdel(1,IS)
          elseif(idry_e(isdel(2,IS)) == 0 .AND. idry_e(isdel(1,IS)) == 1) then 
            IELEM = isdel(2,IS)
          else
            cycle
          endif

          ! We compute the barotropic gradient at the shoreline (only the wet element is used)
          BPGR = 0
          do IP = 1,i34(IELEM)
            ! x and y-components of grad(eta)
            TMP_X = (1-thetai)*eta1(elnode(IP,IELEM))*dldxy(IP,1,IELEM) + thetai*eta2(elnode(IP,IELEM))*dldxy(IP,1,IELEM)
            TMP_Y = (1-thetai)*eta1(elnode(IP,IELEM))*dldxy(IP,2,IELEM) + thetai*eta2(elnode(IP,IELEM))*dldxy(IP,2,IELEM)
            ! Barotropic gradient = g*grad(eta)
            BPGR(1) = BPGR(1) + grav*TMP_X
            BPGR(2) = BPGR(2) + grav*TMP_Y
          enddo !IP

          ! Overwriting wave forces to balance out pressure gradient
          WWAVE_FORCE(1,:,IS) = BPGR(1)
          WWAVE_FORCE(2,:,IS) = BPGR(2)
        enddo !IS

        call exchange_s3d_2(WWAVE_FORCE)

      END SUBROUTINE SHORELINE_WAVE_FORCES
#endif /*SCHISM*/
!**********************************************************************
!*                                                                    *
!**********************************************************************
