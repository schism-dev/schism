#include "wwm_functions.h"
!#define LTESTWAMSOURCES
#undef LTESTWAMSOURCES
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE INIT_ARRAYS
       USE DATAPOOL
       IMPLICIT NONE

#ifdef MPI_PARALL_GRID
       INTEGER :: IE
#endif

       IF (DIMMODE .EQ. 1) THEN
         ALLOCATE( DX1(0:MNP+1), DX2(0:MNP+1), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 1')
         DX1 = zero
         DX2 = zero 
       ENDIF

       ALLOCATE( XP(MNP), YP(MNP), DEP(MNP), stat=istat)
       IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 2')
       XP  = zero
       YP  = zero
       DEP = zero

       ALLOCATE( INVSPHTRANS(MNP,2), stat=istat)
       IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 3')
       INVSPHTRANS = zero

       ALLOCATE( INE(3,MNE), IEN(6,MNE), TRIA(MNE), stat=istat)
       IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 4')
       INE = 0
       IEN = zero
       TRIA = zero
#ifdef MPI_PARALL_GRID
# ifdef PDLIB
       INE(:,:) = INETMP(:,:)
# else
       INE = INETMP(1:3,:)
# endif
#endif
!
! spectral grid - shared
!
       ALLOCATE(MSC_HF(MNP), stat=istat)
       IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 5')
       MSC_HF = MSC
!
! action densities and source terms - shared
!
       ALLOCATE (AC2(MSC,MDC,MNP), stat=istat)
       IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 8')
       AC2 = zero

       ALLOCATE (AC1(MSC,MDC,MNP), stat=istat)
       IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 9')
       AC1 = zero

       IF ((.NOT. BLOCK_GAUSS_SEIDEL).and.(AMETHOD .eq. 7)) THEN
         ALLOCATE (U_JACOBI(MSC,MDC,MNP), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 9a')
         U_JACOBI = zero
       END IF

       IF (ICOMP .GE. 2) THEN
         ALLOCATE (IMATRAA(MSC,MDC,MNP), IMATDAA(MSC,MDC,MNP), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 10')
         IMATRAA = zero
         IMATDAA = zero
       END IF

!       ALLOCATE(SBR(2,MNP),SBF(2,MNP), stat=istat)
!       IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 10a')
!       SBR = ZERO
!       SBF = ZERO

!       ALLOCATE(STOKES_X(NLVT,MNP), STOKES_Y(NLVT,MNP), JPRESS(MNP), stat=istat)
!       IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 10b')
!       STOKES_X = ZERO
!       STOKES_Y = ZERO
!       JPRESS = ZERO

       IF (LITERSPLIT) THEN
!         ALLOCATE (AC1(MNP,MSC,MDC)); AC1 = zero
         ALLOCATE (DAC_ADV(2,MNP,MSC,MDC), DAC_THE(2,MNP,MSC,MDC), DAC_SIG(2,MNP,MSC,MDC), DAC_SOU(2,MNP,MSC,MDC), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 11')
         DAC_ADV = zero
         DAC_THE = zero
         DAC_SIG = zero
         DAC_SOU = zero
       END IF

       IF ((ICOMP .eq. 3).and.(AMETHOD .eq. 7).AND.(ASPAR_LOCAL_LEVEL .eq. 0)) THEN
#ifdef WWM_SOLVER
         IF (REFRACTION_IMPL) THEN
           allocate(CAD_THE(MSC,MDC,NP_RES), stat=istat)
           IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 11.1')
         END IF
         IF (FREQ_SHIFT_IMPL) THEN
           allocate(CAS_SIG(MSC,MDC,NP_RES), stat=istat)
           IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 11.2')
         END IF
#else
         CALL WWM_ABORT('Needs WWM_SOLVER for JACOBI_ITERATION (AMETHOD 7)')
#endif
       END IF

#ifdef SHYFEM_COUPLING
       IF (LSHYFEM) THEN
         ALLOCATE(SHYFZETA(NLVT,MNP),NLEV(MNP), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 12')
         SHYFZETA = ZERO
         NLEV = 0
         ALLOCATE(STOKES_X(NLVT,MNP), STOKES_Y(NLVT,MNP), JPRESS(MNP), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 13')
         STOKES_X = ZERO; STOKES_Y = ZERO; JPRESS = ZERO
       END IF
#endif
!
! Boundary conditions - shared
!
       ALLOCATE( IOBDP(MNP), IOBPD(MDC,MNP), IOBP(MNP), IOBWB(MNP), stat=istat)
       IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 15')
       IOBDP  = 1
       IOBPD  = 0
       IOBP   = 0
       IOBWB  = 1
!
! phase velocity, wave number, group velocity, dwdh, kh
!
       ALLOCATE( WK(MSC,MNP), CG(MSC,MNP), WC(MNP,MSC), stat=istat)
       IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 16')
       WK = ZERO 
       CG = ZERO 
       WC = ZERO
!
! phase velocity, wave number, group velocity, dwdh, kh
!
       IF (MESTR .EQ. 7) THEN
         ALLOCATE( DWKDX(MNP,MSC), DCGDY(MNP,MSC), DWKDY(MNP,MSC), DCGDX(MNP,MSC), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 16')
         DWKDX = ZERO
         DCGDX = ZERO
         DWKDY = ZERO
         DCGDY = ZERO
       ENDIF
!
       ALLOCATE( TABK (0:IDISPTAB), TABCG(0:IDISPTAB), stat=istat)
       IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 17')
       TABK  = zero
       TABCG = zero
!
! diffraction parameter - shared
!
       IF (LDIFR) THEN
         ALLOCATE(DIFRM(MNP), DIFRX(MNP), DIFRY(MNP), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 18')
         DIFRM = zero
         DIFRX = zero
         DIFRY = zero
       END IF
!
! water level, currents and depths ...
!
       ALLOCATE(WINDXY(MNP,2), PRESSURE(MNP), stat=istat)
       IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 19')
       WINDXY = zero
       PRESSURE = zero
       ALLOCATE( DVWIND(MNP,2), CURTXY(MNP,2), DVCURT(MNP,2), stat=istat)
       IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 20')
       DVWIND = zero
       CURTXY = zero
       DVCURT = zero
       ALLOCATE( DDEP(MNP,2), DCUX(MNP,2), DCUY(MNP,2), stat=istat)
       IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 21')
       DDEP = zero
       DCUX = zero
       DCUY = zero
       ALLOCATE( WATLEV(MNP), WATLEVOLD(MNP), stat=istat)
       IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 22')
       WATLEV = zero
       WATLEVOLD = zero
       ALLOCATE( DVWALV(MNP), WLDEP(MNP), DEPDT(MNP), stat=istat)
       IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 23')
       DVWALV = zero
       WLDEP = zero
       DEPDT = zero
!
!  convergence analysis - shared
!
       IF (LCONV .OR. (LQSTEA .AND. LCHKCONV)) THEN
         ALLOCATE ( SUMACOLD(MNP), HSOLD(MNP), KHSOLD(MNP), TM02OLD(MNP), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 24')
         SUMACOLD = zero
         HSOLD  = zero
         TM02OLD = zero
         KHSOLD  = zero
       END IF
!
!  new stuff not ready  
!
       IF (LCONV) THEN ! more work needed ...
         ALLOCATE ( IP_IS_STEADY(MNP), IE_IS_STEADY(MNE), STAT2D(MSC,MDC), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 24')
         IP_IS_STEADY = 0
         IE_IS_STEADY = 0
         STAT2D = ZERO
       END IF
!
!  output - shared
!
       ALLOCATE( QBLOCAL(MNP), DISSIPATION(MNP), AIRMOMENTUM(MNP), stat=istat)
       IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 25')
       QBLOCAL = zero
       DISSIPATION = zero
       AIRMOMENTUM = zero

       ALLOCATE( UFRIC(MNP), ALPHA_CH(MNP), stat=istat)
       IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 26')
       UFRIC = zero
       ALPHA_CH = zero

       ALLOCATE( TAUW(MNP), TAUTOT(MNP), TAUWX(MNP), TAUWY(MNP), TAUHF(MNP), TAUHFT2(0:IUSTAR,0:IALPHA,0:ILEVTAIL), TAUHFT(0:IUSTAR,0:IALPHA,MSC), TAUT(0:ITAUMAX,0:JUMAX,JPLEVT), stat=istat)
       IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 27')
       TAUW = zero
       TAUTOT = zero
       TAUWX = zero
       TAUWY = zero
       TAUHF = zero
       TAUHFT2 = zero

       ALLOCATE( Z0(MNP), CD(MNP), USTDIR(MNP), stat=istat)
       IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 28')
       Z0 = zero
       CD = zero
       USTDIR = zero

       ALLOCATE( RSXX(MNP), RSXY(MNP), RSYY(MNP), FORCEXY(MNP,2), stat=istat)
       IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 29')
       RSXX = zero
       RSXY = zero
       RSYY = zero
       FORCEXY = zero

#ifdef SCHISM
       ALLOCATE( SXX3D(NVRT,MNP), SXY3D(NVRT,MNP), SYY3D(NVRT,MNP), stat=istat)
       IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 30')
       SXX3D = zero
       SXY3D = zero
       SYY3D = zero
#endif
#ifndef SCHISM
       IF (LCPL) THEN
         IF (LTIMOR.or.LSHYFEM) THEN
           ALLOCATE( SXX3D(NLVT,MNP), SXY3D(NLVT,MNP), SYY3D(NLVT,MNP), stat=istat)
           IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 31')
           SXX3D = zero
           SXY3D = zero
           SYY3D = zero
         END IF
       END IF
#endif
       ALLOCATE(HMAX(MNP), ISHALLOW(MNP), stat=istat)
       IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 32')
       HMAX = zero
       ISHALLOW = 0

       IF (LSOURCESWAM .OR. MESIN == 2) THEN
         ALLOCATE( FL(MNP,MDC,MSC), FL3(MNP,MDC,MSC), SL(MNP,MDC,MSC), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 32a')
         FL = ZERO
         FL3 = ZERO
         SL = ZERO
         ALLOCATE( FMEAN(MNP), EMEAN(MNP), FMEANWS(MNP), MIJ(MNP), ENH(MNP,MSC+4,1), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 32b')
         FMEAN = ZERO
         EMEAN = ZERO
         FMEANWS = ZERO
         MIJ = 0
         ENH = 1.d0
         ALLOCATE( USOLD(MNP,1), THWOLD(MNP,1), THWNEW(MNP), Z0OLD(MNP,1), Z0NEW(MNP), ROAIRO(MNP,1), ROAIRN(MNP), U10OLD(MNP,1), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 32c')
         U10OLD = ZERO
         USOLD = ZERO
         THWOLD = ZERO 
         THWNEW = ZERO
         Z0OLD = ZERO
         Z0NEW = ZERO 
         ROAIRO = 1.2250000238 
         ROAIRN = 1.2250000238 
         ALLOCATE( ZIDLOLD(MNP,1), ZIDLNEW(MNP), U10NEW(MNP), USNEW(MNP), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 32d')
         ZIDLOLD = ZERO
         ZIDLNEW = ZERO
         U10NEW = ZERO
         USNEW = ZERO
         ALLOCATE( FCONST(MNP,MSC), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 32e')
         FCONST = 1
       ENDIF
!
!      init source term parameter 
!      
      IF (IPHYS.EQ.0) THEN
!       ECMWF PHYSICS:
        BETAMAX = 1.20
        ZALP    = 0.008
        TAUWSHELTER=0.0
        IF(MSC.GT.30) THEN
          ALPHA   = 0.0060
        ELSE
          ALPHA   = 0.0075
        ENDIF
      ELSE IF (IPHYS.EQ.1) THEN
!       METEO FRANCE PHYSICS:
        BETAMAX = 1.52
        ZALP    = 0.0060
        TAUWSHELTER=0.6
        IF(MSC.GT.30) THEN
          ALPHA   = 0.0090
        ELSE
          ALPHA   = 0.0095
        ENDIF
      ELSE IF (IPHYS.EQ.2) THEN
!      COMBINED ECMWF/METEO FRANCE PHYSICS:
        BETAMAX = 1.52
        ZALP    = 0.0060
        TAUWSHELTER=0.0
        IF(MSC.GT.30) THEN
          ALPHA   = 0.003
        ELSE
          ALPHA   = 0.004
        ENDIF
      ELSE
         CALL WWM_ABORT('UKNOWN PHYSICS SELECTION') 
      ENDIF ! IPHYS

      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'LEAVING INIT_ARRAYS'
      FLUSH(STAT%FHNDL)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DEALLOC_ARRAYS
      USE DATAPOOL
      IMPLICIT NONE

      IF (DIMMODE == 1) THEN
        DEALLOCATE( DX1, DX2)
      END IF

      DEALLOCATE(XP, YP, DEP)
      DEALLOCATE(INVSPHTRANS)
      DEALLOCATE(INE, IEN, TRIA)
!
! spectral grid - shared
!
      DEALLOCATE( SPSIG, SPDIR, FR, COSTH, SINTH, COS2TH, SIN2TH)
      DEALLOCATE( SINCOSTH, SIGPOW, DS_BAND, DS_INCR)
      DEALLOCATE(MSC_HF)
!
! action densities and source terms - shared
!
      IF ((.NOT. BLOCK_GAUSS_SEIDEL).and.(AMETHOD .eq. 7)) THEN
        DEALLOCATE (U_JACOBI)
      END IF
      DEALLOCATE (AC2, AC1)
      IF (ICOMP .GE. 2) THEN
        DEALLOCATE (IMATRAA, IMATDAA)
      END IF

      IF (LITERSPLIT) THEN
        DEALLOCATE (DAC_ADV, DAC_THE, DAC_SIG, DAC_SOU)
      END IF

#ifdef WWM_SOLVER
      IF ((ICOMP .eq. 3).and.(AMETHOD .eq. 7).AND.(ASPAR_LOCAL_LEVEL .eq. 0)) THEN
         IF (REFRACTION_IMPL) THEN
           deallocate(CAD_THE)
         END IF
         IF (FREQ_SHIFT_IMPL) THEN
           deallocate(CAS_SIG)
         END IF
      END IF
#endif
#ifdef SHYFEM_COUPLING
      IF (LSHYFEM) THEN
        DEALLOCATE(SHYFZETA, NLEV)
      END IF
#endif
!
! WAM Cycle 4.5 - shared
!
      IF (MESIN .EQ. 2) THEN
        DEALLOCATE ( TAUHFT, TAUT)
      ENDIF
!
! Boundary conditions - shared
!
      DEALLOCATE( IOBPD, IOBP, IOBWB)
!
! phase velocity, wave number, group velocity, dwdh, kh
!
      DEALLOCATE( WK, CG, TABK, TABCG)
!
! diffraction parameter - shared
!
      IF (LDIFR) THEN
        DEALLOCATE ( DIFRM, DIFRX, DIFRY )
      END IF
!
! water level, currents and depths ...
!
      DEALLOCATE( WINDXY, PRESSURE, DVWIND, CURTXY, DVCURT)
      DEALLOCATE( DDEP, DCUX, DCUY, WATLEV, WATLEVOLD)
      DEALLOCATE( DVWALV, WLDEP, DEPDT)
!
!  convergence analysis - shared
!
      IF (LCONV .OR. (LQSTEA .AND. LCHKCONV)) THEN
        DEALLOCATE ( SUMACOLD, HSOLD, KHSOLD, TM02OLD)
      END IF
!
!  output - shared
!
      DEALLOCATE( QBLOCAL, DISSIPATION, AIRMOMENTUM, UFRIC, ALPHA_CH)
      DEALLOCATE( TAUW, TAUTOT, TAUWX, TAUWY, TAUHF)
      DEALLOCATE( Z0, CD, USTDIR)
      DEALLOCATE( RSXX, RSXY, RSYY)
#ifdef SCHISM
      DEALLOCATE( SXX3D, SXY3D, SYY3D)
#endif
#ifndef SCHISM
      IF (LCPL) THEN
        IF (LTIMOR.or.LSHYFEM) THEN
          DEALLOCATE( SXX3D, SXY3D, SYY3D)
        END IF
      END IF
#endif
      DEALLOCATE(HMAX, ISHALLOW)
#ifdef NCDF
      IF (GRIDWRITE) THEN
        DEALLOCATE(XPtotal, YPtotal, IOBPtotal, DEPtotal, INEtotal)
      END IF
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INITIALIZE_WWM
#if defined ROMS_WWM_PGMCL_COUPLING || defined MODEL_COUPLING_ATM_WAV || defined MODEL_COUPLING_OCN_WAV
      USE WWMaOCN_PGMCL
#endif
#if defined MODEL_COUPLING_ATM_WAV || defined MODEL_COUPLING_OCN_WAV
      USE pgmcl_lib_WWM, only : WAV_common_initialize
#endif
      USE DATAPOOL
#ifdef SCHISM
      use schism_glbl, only: xlon,ylat
#endif
#ifdef PDLIB
      USE yowpd, only : initFromGridDim
#endif
#if !defined PDLIB && defined MPI_PARALL_GRID
      USE schism_glbl, only : ics
#endif
#ifdef PETSC
      USE PETSC_CONTROLLER, ONLY : PETSC_INIT
#endif
#ifdef ST41
      USE W3SRC4MD_OLD
#endif
#ifdef ST42
      USE W3SRC4MD
#endif
      IMPLICIT NONE
!
      INTEGER        :: i, j, IP
      REAL(rkind)    :: TIME1, TIME2
      
#ifdef TIMINGS
      CALL WAV_MY_WTIME(TIME1)
#endif
     
      CALL SET_WWMINPULNML 

! variable nx1 should be initialized in SCHISM code, not here!
#if defined MPI_PARALL_GRID && !defined PDLIB
      do i=1,3
        do j=1,2
          nx1(i,j)=i+j
          if(nx1(i,j)>3) nx1(i,j)=nx1(i,j)-3
          if(nx1(i,j)<1.or.nx1(i,j)>3) then
            write(wwmerr,*)'MAIN: nx1 wrong',i,j,nx1(i,j)
            call wwm_abort(wwmerr)
          endif
        enddo
      enddo
#endif
      CALL INIT_FILE_HANDLES
      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'DONE SETTING FHNDL'
      FLUSH(STAT%FHNDL)

      CALL READ_WWMINPUT
      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'DONE READING NAMELIST'
      FLUSH(STAT%FHNDL)
      CALL READ_SPATIAL_GRID_TOTAL

#ifndef MPI_PARALL_GRID
      MNP=NP_TOTAL
      MNE=NE_TOTAL
      NP_RES=MNP
      CALL INIT_ARRAYS
      XP=XPtotal
      YP=YPtotal
      DEP=DEPtotal
      INE=INEtotal
#else
# ifdef PDLIB
      IF (IGRIDTYPE .eq. 2) THEN
        CALL WWM_ABORT('Not yet support for PDLIB and IGRIDTYPE=2')
      END IF
      write(DBG%FHNDL,*) 'sum(XPtotal)=', sum(XPtotal)
      write(DBG%FHNDL,*) 'sum(YPtotal)=', sum(YPtotal)
      write(DBG%FHNDL,*) 'sum(DEPtotal)=', sum(DEPtotal)
      write(DBG%FHNDL,*) 'sum(INEtotal)=', sum(INEtotal)
      write(DBG%FHNDL,*) NP_TOTAL, NE_TOTAL, MDC, MSC
      CALL initFromGridDim(NP_TOTAL, XPtotal, YPtotal, DEPtotal, NE_TOTAL, INEtotal, MDC, MSC, comm)
      call fillPublicVars
      CALL INIT_ARRAYS
      XP = XPTMP
      YP = YPTMP
      DEP=DEP8
      INETMP=INE
# else
#  ifndef SCHISM
      call partition_hgrid
      call aquire_hgrid(.true.)
      call msgp_tables
      call msgp_init
      call parallel_barrier
#  endif
      CALL INIT_ARRAYS

      DEP  = DEP8
      IF (ics .eq. 2) THEN
        XP = XLON*RADDEG
        YP = YLAT*RADDEG
      ELSE
        XP = XPTMP
        YP = YPTMP
      END IF
# endif
      CALL COLLECT_ALL_IPLG
      CALL SETUP_ONED_SCATTER_ARRAY
#endif
      WLDEP=DEP
      IF (CART2LATLON) THEN
        XP = XP / 111111.
        YP = YP / 111111.
      ELSE IF (LATLON2CART) THEN
        XP = XP * 111111.
        YP = YP * 111111. 
      ELSE IF (CART2LATLON .AND. LATLON2CART) THEN
        CALL  WWM_ABORT('CART2LATLON .AND. LATLON2CART cannot be T')
      ENDIF 
      CALL INIT_SPATIAL_GRID
      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'INIT SPATIAL GRID'
      FLUSH(STAT%FHNDL)
      !
      ! Main inits done, now the secondary ones.
      !
#ifdef VDISLIN
      CALL INIT_DISLIN()
      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'INIT DISLIN                '
      FLUSH(STAT%FHNDL)
#endif

      CALL CHECK_LOGICS
      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'CHECK LOGICS                '
      FLUSH(STAT%FHNDL)
      !
      IF (LADVTEST) THEN
        ALLOCATE(UTEST(MNP), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 33')
        UTEST = 0.
        CALL ADVTEST(UTEST)
        AC2(1,1,:) = UTEST
        CALL CHECKCONS(UTEST,SUMACt0)
!AR: check the advections test if it still works ...
        DEALLOCATE(UTEST)
      END IF
#ifdef MPI_PARALL_GRID
      CALL BUILD_WILD_ARRAY
#endif
      CALL BUILD_IPSTATUS
      CALL BUILD_TRIANGLE_CORRESPONDENCES
      CALL SET_IOBPD_BY_DEP
      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'SET DEPTH POINTER'
      FLUSH(STAT%FHNDL)

      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'INITIALIZE SPECTRAL GRID'
      FLUSH(STAT%FHNDL)
      CALL INIT_SPECTRAL_GRID

      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'INITIALIZE BOUNDARY POINTER 1/2'
      FLUSH(STAT%FHNDL)
#if defined SCHISM
!AR: let dmin free ...
!      DMIN = DMIN_SCHISM
#endif
      CALL SET_IOBP_NEXTGENERATION
      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'INITIALIZE BOUNDARY POINTER 2/2'
      FLUSH(STAT%FHNDL)
      CALL SET_IOBPD

      IF (DIMMODE .EQ. 2) THEN
        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'THE FLUCTUATION SPLITTING PREPROCESSOR HAS STARTED'
        FLUSH(STAT%FHNDL)
        CALL INIT_FLUCT_ARRAYS
        CALL INIT_FLUCT

        IF (AMETHOD .EQ. 4 .OR. AMETHOD .EQ. 5) THEN
#ifdef PETSC
          CALL PETSC_INIT
#endif
        END IF
        IF (AMETHOD .EQ. 6) THEN
#if defined WWM_SOLVER && defined MPI_PARALL_GRID
          CALL WWM_SOLVER_INIT
#endif
        END IF
        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'THE FLUCTUATION SPLITTING PREPROCESSOR HAS ENDED'
        FLUSH(STAT%FHNDL)
      END IF

      IF (LZETA_SETUP) THEN
#ifdef WWM_SETUP
        CALL INIT_WAVE_SETUP
#else
        CALL WWM_ABORT('Need WWM_SEZUP if LZETA is selected')
#endif
      END IF


      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'INITIALIZE WIND CURRENT WATERLEVEL'
      FLUSH(STAT%FHNDL)
#if !defined ROMS_WWM_PGMCL_COUPLING && !defined MODEL_COUPLING_ATM_WAV && !defined MODEL_COUPLING_OCN_WAV
      IF (LWINDFROMWWM) THEN
        CALL INIT_WIND_INPUT
      END IF
#endif
#ifdef SCHISM
      IF (.NOT. LWINDFROMWWM) THEN
        WINDXY(:,1) = WINDX0
        WINDXY(:,2) = WINDY0
      END IF
#endif
#if !defined ROMS_WWM_PGMCL_COUPLING && !defined MODEL_COUPLING_OCN_WAV
      IF (.NOT. LCPL) THEN
        CALL INIT_CURRENT_INPUT
        CALL INIT_WATLEV_INPUT
      END IF
#endif
     

      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'COMPUTE THE WAVE PARAMETER'
      FLUSH(STAT%FHNDL)
      CALL INITIATE_WAVE_PARAMETER
      CALL SETSHALLOW
      CALL SET_HMAX

      IF ( (MESIN .EQ. 1 .OR. MESDS .EQ. 1) .AND. SMETHOD .GT. 0 .AND. (LSOURCESWAM .OR. LSOURCESWWIII)) THEN
        CALL WWM_ABORT('DO NOT USE MESIN = 1 OR MESDS = 1 together with LSOURCESWAM .OR. LSOURCESWWIII')
      END IF

      IF ( (MESIN .EQ. 1 .OR. MESDS .EQ. 1) .AND. SMETHOD .GT. 0 .AND. .NOT. (LSOURCESWAM .OR. LSOURCESWWIII)) THEN
        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'INIT ARDHUIN et al.'
        FLUSH(STAT%FHNDL)
#ifdef ST41
        CALL PREPARE_ARDHUIN_OLD
#elif ST42
        CALL PREPARE_ARDHUIN
#endif
      ENDIF
      
      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'SET THE INITIAL WAVE BOUNDARY CONDITION'
      FLUSH(STAT%FHNDL)
      CALL INIT_WAVE_BOUNDARY_CONDITION
      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'SET THE INITIAL CONDITION'
      FLUSH(STAT%FHNDL)
      CALL INITIAL_CONDITION
!      CALL Print_SumAC2("After INITIAL_CONDITION")
      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'INIT STATION OUTPUT'
      FLUSH(STAT%FHNDL)
      CALL INIT_STATION_OUTPUT
      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'WRITING INITIAL TIME STEP'
      FLUSH(STAT%FHNDL)
#ifdef MPI_PARALL_GRID
      CALL EXCHANGE_P4D_WWM(AC2)
#endif
      CALL WWM_OUTPUT(ZERO,.TRUE.)
      IF (LWXFN) THEN
        CALL WRINPGRD_XFN
      ELSE IF(LWSHP) THEN
        CALL WRINPGRD_SHP
      END IF
#if !defined SCHISM && !defined ROMS_WWM_PGMCL_COUPLING && !defined MODEL_COUPLING_ATM_WAV && !defined MODEL_COUPLING_OCN_WAV
      IF (LCPL) THEN
        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'OPEN PIPES FOR COUPLING'
        FLUSH(STAT%FHNDL)
        IF (LTIMOR) THEN
          CALL INIT_PIPES_TIMOR()
#ifdef SHYFEM_COUPLING
        ELSE IF (LSHYFEM) THEN
          CALL INIT_PIPES_SHYFEM()
#endif
        ELSE IF (LROMS) THEN
          CALL INIT_PIPES_ROMS()
        END IF
      END IF
#endif
#ifdef ROMS_WWM_PGMCL_COUPLING
      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'Before ROMS_COUPL_INITIALIZE'
      FLUSH(STAT%FHNDL)
      CALL WWM_common_coupl_initialize
      CALL WWM_a_OCN_COUPL_INITIALIZE
      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'After ROMS_COUPL_INITIALIZE'
      FLUSH(STAT%FHNDL)
#endif
#if defined MODEL_COUPLING_ATM_WAV || defined MODEL_COUPLING_OCN_WAV
      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'Before WAV_common_initialize'
      FLUSH(STAT%FHNDL)
      CALL WAV_common_initialize
      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'After WAV_common_initialize'
      FLUSH(STAT%FHNDL)
#endif
#ifdef TIMINGS
      CALL WAV_MY_WTIME(TIME2)
#endif

#if defined SCHISM
      IF (MSC_SCHISM .NE. MSC .OR. MDC_SCHISM .NE. MDC) THEN
        WRITE(DBG%FHNDL,*) 'MSC_SCHISM', MSC_SCHISM
        WRITE(DBG%FHNDL,*) 'MSC', MSC
        WRITE(DBG%FHNDL,*) 'MDC_SCHISM', MDC_SCHISM
        WRITE(DBG%FHNDL,*) 'MDC', MDC
        FLUSH(DBG%FHNDL)
        CALL PARALLEL_ABORT('THERE IS AND ERROR IN MSC2 OR MDC2 IN PARAM.IN')
      END IF
#endif

      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'INITIALIZE_WWM'
#ifdef TIMINGS
      WRITE(STAT%FHNDL,'("+TRACE...",A,F15.4)') 'CPU Time for the preprocessing', TIME2-TIME1
#endif
      FLUSH(STAT%FHNDL)
      AC1 = AC2
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE TERMINATE_WWM
#ifdef ROMS_WWM_PGMCL_COUPLING
       USE WWMaOCN_PGMCL
#endif
#if defined MODEL_COUPLING_ATM_WAV || defined MODEL_COUPLING_OCN_WAV
       USE pgmcl_lib_WWM, only : WAV_all_deallocate
#endif
       USE DATAPOOL
#ifdef ST41
       USE W3SRC4MD_OLD
#endif
#ifdef ST42
       USE W3SRC4MD
#endif

#ifdef PETSC
       USE PETSC_CONTROLLER, ONLY : PETSC_FINALIZE
#endif

       CALL CLOSE_FILE_HANDLES
       CALL DEALLOC_ARRAYS

#ifdef MPI_PARALL_GRID
       CALL DEALLOC_WILD_ARRAY
#endif
       IF (DIMMODE .EQ. 2) THEN
         CALL DEALLOC_FLUCT_ARRAYS
         CALL DEALLOC_FLUCT

         IF (AMETHOD .EQ. 4 .OR. AMETHOD .EQ. 5) THEN
#ifdef PETSC
           CALL PETSC_FINALIZE
#endif
         END IF
       END IF
       IF (LZETA_SETUP) THEN
         CALL FINALIZE_WAVE_SETUP
       END IF
       CALL DEALLOC_SPECTRAL_GRID
       CALL CLOSE_IOBP
       CALL TERMINATE_STATION_OUTPUT

#if !defined SCHISM && !defined ROMS_WWM_PGMCL_COUPLING && !defined MODEL_COUPLING_ATM_WAV && !defined MODEL_COUPLING_OCN_WAV
       IF (LCPL) THEN
         IF (LTIMOR) THEN
           CALL TERMINATE_PIPES_TIMOR()
# ifdef SHYFEM_COUPLING
         ELSE IF (LSHYFEM) THEN
           CALL TERMINATE_PIPES_SHYFEM()
# endif
         ELSE IF (LROMS) THEN
           CALL TERMINATE_PIPES_ROMS()
         END IF
       END IF
#endif
#ifdef ROMS_WWM_PGMCL_COUPLING
       CALL WWM_a_OCN_COUPL_DEALLOCATE
#endif
#if defined MODEL_COUPLING_ATM_WAV || defined MODEL_COUPLING_OCN_WAV
       CALL WAV_all_deallocate
#endif
       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef MPI_PARALL_GRID 
       SUBROUTINE BUILD_WILD_ARRAY
       USE DATAPOOL, only : MNP, NP_RES, ONE, rkind
       USE DATAPOOL, only : nwild_gb, nwild_loc, nwild_loc_res
       USE datapool, only : iplg, np_global
       USE datapool, only : istatus, itype, comm, ierr, myrank, rtype, nproc
       use datapool, only : MPI_SUM
       IMPLICIT NONE
       integer, allocatable :: nwild_i(:), nwild_gbi(:)
       integer :: iProc, istat
       INTEGER :: IP
       real(rkind), allocatable :: nwild_gb_res(:)
       ALLOCATE(nwild_i(np_global), nwild_gbi(np_global), stat=istat)
       IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 33')
       !
       nwild_i=0
       nwild_gbi=0
       DO IP=1,MNP
         nwild_i(IPLG(IP))=1
       END DO
       call mpi_reduce(nwild_i,nwild_gbi,NP_GLOBAL,itype,MPI_SUM,0,comm,ierr)
       ALLOCATE(nwild_gb(NP_GLOBAL), stat=istat)
       IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 34')
       IF (myrank.eq.0) THEN
         DO IP=1,NP_GLOBAL
           nwild_gb(IP)=ONE/MyREAL(nwild_gbi(IP))
         END DO
         DO iProc=2,nproc
           CALL MPI_SEND(nwild_gb,np_global,rtype, iProc-1, 2037, comm, ierr)
         END DO
       ELSE
         CALL MPI_RECV(nwild_gb,np_global,rtype, 0, 2037, comm, istatus, ierr)
       END IF
       ALLOCATE(nwild_loc(MNP), stat=istat)
       IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 35')
       DO IP=1,MNP
         nwild_loc(IP)=nwild_gb(iplg(IP))
       END DO
       IF (myrank.ne.0) THEN
         deallocate(nwild_gb)
       END IF
       !
       nwild_i=0
       nwild_gbi=0
       DO IP=1,NP_RES
         nwild_i(IPLG(IP))=1
       END DO
       call mpi_reduce(nwild_i,nwild_gbi,NP_GLOBAL,itype,MPI_SUM,0,comm,ierr)
       ALLOCATE(nwild_gb_res(NP_GLOBAL), stat=istat)
       IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 36')
       IF (myrank.eq.0) THEN
         DO IP=1,NP_GLOBAL
           nwild_gb_res(IP)=ONE/MyREAL(nwild_gbi(IP))
         END DO
         DO iProc=2,nproc
           CALL MPI_SEND(nwild_gb_res,np_global,rtype, iProc-1, 277, comm, ierr)
         END DO
       ELSE
         CALL MPI_RECV(nwild_gb_res,np_global,rtype, 0, 277, comm, istatus, ierr)
       END IF
       ALLOCATE(nwild_loc_res(NP_RES), stat=istat)
       IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 37')
       DO IP=1,NP_RES
         nwild_loc_res(IP)=nwild_gb_res(iplg(IP))
       END DO
       deallocate(nwild_gb_res)
       !
       DEALLOCATE(nwild_i)
       DEALLOCATE(nwild_gbi)
       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE DEALLOC_WILD_ARRAY
!AR: what is this kind of shit ?
       USE DATAPOOL
       implicit none
       DEALLOCATE(nwild_loc)
       DEALLOCATE(nwild_loc_res)
       END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE BASIC_PARAMETER
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER :: IP
         REAL(rkind)    :: PPTAIL
!
!2do: Check with WAM parameterization and WW3
!
         PQUAD(1) = 0.25
         IF (MESIN .EQ. 1) THEN
           PQUAD(2) = 2.5E7
         ELSE
           PQUAD(2) = 2.78E7
         END IF
         PQUAD(3) = 5.5_rkind
         PQUAD(4) = 6._rkind/7._rkind
         PQUAD(5) = -1.25_rkind
         PQUAD(6) = zero
!
! High frequency tail integration factors as defined in the SWAN model.
!
         PTAIL(1) = 4._rkind
         PTAIL(2) = 2.5_rkind
         PTAIL(3) = PTAIL(1) + 1._rkind
         IF (MESIN .EQ. 2) THEN
           PTAIL(1) = 5._rkind
           PTAIL(3) = PTAIL(1) + 1._rkind
         ENDIF
         PTAIL(4) = 3._rkind
         DO IP = 0, 3
           PPTAIL = PTAIL(1) - MyREAL(IP)
           PTAIL(5+IP) = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
         ENDDO
!
! Output Variables ... LVARS, NVARS
!
         OUTT_VARNAMES(1)  = 'HS'
         OUTT_VARNAMES(2)  = 'TM01'
         OUTT_VARNAMES(3)  = 'TM02'
         OUTT_VARNAMES(4)  = 'TM10'
         OUTT_VARNAMES(5)  = 'KLM'
         OUTT_VARNAMES(6)  = 'WLM'
         OUTT_VARNAMES(7)  = 'ETOTS'
         OUTT_VARNAMES(8)  = 'ETOTC'
         OUTT_VARNAMES(9)  = 'DM'
         OUTT_VARNAMES(10)  = 'DSPR'
         OUTT_VARNAMES(11) = 'TPPD'
         OUTT_VARNAMES(12) = 'TPP'
         OUTT_VARNAMES(13) = 'CPP'
         OUTT_VARNAMES(14) = 'WNPP'
         OUTT_VARNAMES(15) = 'CGPP'
         OUTT_VARNAMES(16) = 'KPP'
         OUTT_VARNAMES(17) = 'LPP'
         OUTT_VARNAMES(18) = 'PEAKD'
         OUTT_VARNAMES(19) = 'PEAKDSPR'
         OUTT_VARNAMES(20) = 'DPEAK'
         OUTT_VARNAMES(21) = 'UBOT'
         OUTT_VARNAMES(22) = 'ORBITAL'
         OUTT_VARNAMES(23) = 'BOTEXPER'
         OUTT_VARNAMES(24) = 'TMBOT'
         OUTT_VARNAMES(25) = 'URSELL'
         OUTT_VARNAMES(26) = 'USTAR'
         OUTT_VARNAMES(27) = 'ALPHA'
         OUTT_VARNAMES(28) = 'Z0'
         OUTT_VARNAMES(29) = 'WIND-X'
         OUTT_VARNAMES(30) = 'WIND-Y'
         OUTT_VARNAMES(31) = 'CD'
         OUTT_VARNAMES(32) = 'CURR-X'
         OUTT_VARNAMES(33) = 'CURR-Y'
         OUTT_VARNAMES(34) = 'DEPTH'
         OUTT_VARNAMES(35) = 'ELEVATION'
       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE INITIATE_WAVE_PARAMETER
         USE DATAPOOL, ONLY: STAT, LSTCU, LSECU, MESNL, SPSIG, SPDIR, MSC, MDC, DELALP
         USE DATAPOOL, ONLY: G9, DEP, MNP, MESTR, LSOURCESWWIII, LSOURCESWAM, DELTAIL
         USE DATAPOOL, ONLY: LPRECOMP_EXIST, TAUHFT, TAUHFT2, TAUT, DELU, DELTAUW, DELUST
         USE DATAPOOL, ONLY: IPHYS, MESIN, SMETHOD
#ifdef MPI_PARALL_GRID
         USE DATAPOOL, ONLY: MYRANK,COMM
#endif
         USE M_CONSTANTS
         USE M_XNLDATA
         USE M_FILEIO

         IMPLICIT NONE
         INTEGER IQGRID, INODE, IERR

         WRITE(STAT%FHNDL,*) 'START WAVE PARAMETER'
         FLUSH(STAT%FHNDL)
         CALL GRADDEP
         WRITE(STAT%FHNDL,*) 'GRADDEP'
         FLUSH(STAT%FHNDL)
         IF (LSTCU .OR. LSECU) CALL GRADCURT
         WRITE(STAT%FHNDL,*) 'GRADCURT'
         FLUSH(STAT%FHNDL)
         CALL BASIC_PARAMETER
         WRITE(STAT%FHNDL,*) 'BASIC'
         FLUSH(STAT%FHNDL)
         CALL MAKE_WAVE_TABLE
         CALL WAVE_K_C_CG
         WRITE(STAT%FHNDL,*) 'WAVEKCG'
         FLUSH(STAT%FHNDL)

         IF (.NOT. (LSOURCESWAM .OR. LSOURCESWWIII) .AND. SMETHOD .GT. 0) THEN
           IF (MESNL .LT. 5) THEN
             CALL PARAMETER4SNL
           ELSE IF (MESNL .EQ. 5) THEN
             IQGRID = 3
             INODE  = 1
             CALL XNL_INIT(REAL(SPSIG),REAL(SPDIR),MSC,MDC,-4.0,REAL(G9),REAL(DEP),MNP,1,IQGRID,INODE,IERR)
             CALL INIT_CONSTANTS
             CALL XNL_INIT(REAL(SPSIG),REAL(SPDIR),MSC,MDC,-4.0,REAL(G9),REAL(DEP),MNP,1,IQGRID,INODE,IERR)
             IF (IERR .GT. 0) CALL WWM_ABORT('IERR XNL_INIT')
           ENDIF
         ELSE IF (LSOURCESWAM .AND. .NOT. LSOURCESWWIII) THEN
           WRITE(STAT%FHNDL,'("+TRACE...",A)')'COMPUTING NONLINEAR COEFFICIENTS' 
           CALL NLWEIGT
           WRITE(STAT%FHNDL,'("+TRACE...",A)')'COMPUTING NONLINEAR COEFFICIENTS'
           CALL INISNONLIN
           INQUIRE(FILE='fort.5011',EXIST=LPRECOMP_EXIST)
#ifdef MPI_PARALL_GRID
           CALL MPI_BARRIER(COMM, ierr)
#endif
           IF (LPRECOMP_EXIST) THEN
             WRITE(STAT%FHNDL,'("+TRACE...",A)')'READING STRESS TABLES'
             OPEN(5011, FILE='fort.5011', FORM='UNFORMATTED') 
             IF (IPHYS == 0) THEN
               READ(5011) DELU, DELTAUW
               READ(5011) TAUT
               READ(5011) DELALP, DELUST, DELTAIL
               READ(5011) TAUHFT
             ELSE
               READ(5011) DELU, DELTAUW
               READ(5011) TAUT
               READ(5011) DELALP, DELUST, DELTAIL
               READ(5011) TAUHFT, TAUHFT2
             ENDIF
           ELSE
             !WRITE(*,*) SIZE(TAUT)
             !WRITE(*,*) SIZE(TAUHFT)
             IF (MESIN .GT. 0) THEN 
               WRITE(STAT%FHNDL,'("+TRACE...",A)')'COMPUTING STRESS TABLES'
               CALL STRESS
               WRITE(STAT%FHNDL,'("+TRACE...",A)')'COMPUTING HF TABLES'
               CALL TAUHF_WAM(MSC)
             ENDIF
           ENDIF

           WRITE(STAT%FHNDL,'("+TRACE...",A)')'INITIALIZING STRESS ARRAYS'
           IF (MESIN .GT. 0) CALL BUILDSTRESS
         ELSE IF (LSOURCESWWIII .AND. .NOT. LSOURCESWAM) THEN
         ENDIF

         IF (MESTR == 6) CALL GRAD_CG_K 
       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE INITIAL_CONDITION
       USE WWM_HOTFILE_MOD
       USE DATAPOOL
#ifdef NCDF
       USE NETCDF
#endif
       IMPLICIT NONE
       INTEGER         :: IP, I, K, L, M, IS, ID
       REAL(rkind)     :: SPPAR(8)
       REAL(rkind)     :: MS
       REAL(rkind)     :: HS, TP, HSLESS, TPLESS, FDLESS
       REAL(rkind)     :: WIND10, WINDTH, VEC2DEG
       REAL(rkind)     :: WINDX, WINDY
       REAL(rkind)     :: ACLOC(MSC,MDC)
       REAL(rkind)     :: DEG
       REAL(rkind)     :: TMPPAR(8,MNP), SSBRL(MSC,MDC)
#ifdef NCDF
       CHARACTER(len=25) :: CALLFROM
#endif
       TMPPAR = 0.
       IF (.NOT. LHOTR .AND. LINID) THEN
         IF (INITSTYLE == 2 .AND. IBOUNDFORMAT == 3) THEN
#ifdef NCDF
           CALL READ_NETCDF_WW3_PARAM
#else
           CALL WWM_ABORT('compile with DNCDF PPFLAG')
#endif
           CALL INTER_STRUCT_DOMAIN(NDX_BND,NDY_BND,DX_BND,DY_BND,OFFSET_X_BND,OFFSET_Y_BND,TMPPAR)
!test the initial conditions ...
           WRITE(2001) 1.
           WRITE(2001) (TMPPAR(3,IP), TMPPAR(2,IP), TMPPAR(1,IP), IP = 1, MNP)
         END IF
         DO IP = 1, MNP
           !WRITE(*,*) 'wwm_initio.F90 l.1022', IP, SUM(AC2(:,:,IP))
           IF (ABS(IOBP(IP)) .GT. 0 .AND. .NOT. LSOUBOUND) THEN
             AC2(:,:,IP) = ZERO
             CYCLE
           ENDIF
           ACLOC = 0.
           WINDX  = WINDXY(IP,1)
           WINDY  = WINDXY(IP,2)
           WIND10 = SQRT(WINDX**2+WINDY**2)
           IF(DEP(IP) .GT. DMIN) THEN
             IF (DIMMODE .EQ. 1 .AND. IP .EQ. 1) CYCLE
             IF (INITSTYLE == 1) THEN
               IF (WIND10 .GT. 1.) THEN !AR: why one? 
                 WINDTH = VEC2DEG(WINDX,WINDY)
                 CALL DEG2NAUT(WINDTH, DEG, LNAUTIN)
                 FDLESS = G9*AVETL/WIND10**2
                 HSLESS = MIN(0.21_rkind, 0.00288_rkind*FDLESS**0.45_rkind)
                 TPLESS = MIN(ONE/0.13_rkind, 0.46_rkind*FDLESS**0.27_rkind)
                 HS = HSLESS*WIND10**2/G9
                 TP = TPLESS*WIND10/G9
                 MS = 2.0_rkind
                 SPPAR(1) = HS
                 SPPAR(2) = TP ! Check whether TP is inside of spectal representation !!!
                 SPPAR(3) = DEG
                 SPPAR(4) = MS
                 SPPAR(5) = 2.
                 SPPAR(6) = 2.
                 SPPAR(7) = 0.1
                 SPPAR(8) = 3.3
                 CALL SPECTRAL_SHAPE(SPPAR,ACLOC,.FALSE.,'INITIAL CONDITION PARA', .FALSE.)
                 AC2(:,:,IP) = ACLOC
               ELSE
                 ACLOC = 1.E-8
               END IF
             ELSE IF (INITSTYLE == 2 .AND. IBOUNDFORMAT == 3) THEN
               TMPPAR(5,IP) = 2.
               TMPPAR(6,IP) = 1.
               TMPPAR(7,IP) = 0.1
               TMPPAR(8,IP) = 3.3
               CALL SPECTRAL_SHAPE(TMPPAR(:,IP),ACLOC,.FALSE.,'INITIAL CONDITION WW3', .FALSE.)
               AC2(:,:,IP) = ACLOC
             ELSE IF (INITSTYLE == 3) THEN
               OPEN(1113,FILE='fort.10003',STATUS='OLD')
               DO ID=1,MDC
                 DO IS=1,MSC
                   READ(1113,*) K, M, AC2(IS,ID,IP)
                   AC2(IS,ID,IP) =  AC2(IS,ID,IP) / PI2 / SPSIG(IS)
                 ENDDO
               ENDDO
               REWIND(1113)
             END IF ! INITSTYLE
             IF (LMAXETOT .AND. ISHALLOW(IP) .EQ. 1) CALL BREAK_LIMIT(IP,ACLOC,SSBRL) ! Miche for initial cond.
           ELSE
             FDLESS      = 0.
             HS          = 0.
             TP          = 0.
             WINDTH      = 0.
             AC2(:,:,IP) = 0.
           END IF ! DEP(IP) .GT. DMIN .AND. WIND10 .GT. SMALL
         END DO ! IP
       ELSE IF (LHOTR .AND. .NOT. LINID) THEN
         CALL INPUT_HOTFILE
       END IF
       CALL SET_WAVE_BOUNDARY
       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE INIT_FILE_HANDLES()
         USE DATAPOOL
         IMPLICIT NONE
#ifdef MPI_PARALL_GRID
         CHARACTER (LEN = 30) :: FDB
         INTEGER              :: LFDB
#endif
            CHK%FNAME  = 'wwmcheck.nml'
           QSTEA%FNAME = 'qstea.out'
         WINDBG%FNAME  = 'winddbg.out'
#if defined DEBUG && defined IOBPDOUT
         IOBPOUT%FNAME = 'iobp.out'
        IOBPDOUT%FNAME = 'iobpd.out'
#endif
        SRCDBG%FNAME = 'srcdbg.out'
!
!2do ... dinstinguish between binary and ascii stuff ...
!
             BND%FHNDL  = STARTHNDL + 1
             WIN%FHNDL  = STARTHNDL + 2
             CUR%FHNDL  = STARTHNDL + 3
             WAT%FHNDL  = STARTHNDL + 4
             WAV%FHNDL  = STARTHNDL + 5
             CHK%FHNDL  = STARTHNDL + 6
           HOTIN%FHNDL  = STARTHNDL + 7
          HOTOUT%FHNDL  = STARTHNDL + 8
             INP%FHNDL  = STARTHNDL + 9
             GRD%FHNDL  = STARTHNDL + 10
          GRDCOR%FHNDL  = STARTHNDL + 11

           IF (LQSTEA) QSTEA%FHNDL  = STARTHNDL + 12

#if defined DEBUG && defined IOBPDOUT
         IOBPOUT%FHNDL  = STARTHNDL + 13
         IOBPDOUT%FHNDL = STARTHNDL + 14
#endif

         DBG%FHNDL      = STARTHNDL + 15 
         STAT%FHNDL     = STARTHNDL + 16 
         WINDBG%FHNDL   = STARTHNDL + 17 
         SRCDBG%FHNDL   = STARTHNDL + 18

         IU06           = STAT%FHNDL

#ifndef MPI_PARALL_GRID
         open(DBG%FHNDL,file='wwmdbg.out',status='unknown') !non-fatal errors
         open(STAT%FHNDL,file='wwmstat.out',status='unknown') !non-fatal errors
         open(WINDBG%FHNDL,file='windbg.out',status='unknown') !non-fatal errors
         open(SRCDBG%FHNDL,file='srcdbg.out',status='unknown') !non-fatal errors
#else
# ifdef SCHISM
         FDB  ='wwmdbg_0000'
         LFDB =len_trim(FDB)
         write(FDB(LFDB-3:LFDB),'(i4.4)') MYRANK
         open(DBG%FHNDL,file='outputs/'//fdb,status='replace') 
         FDB  ='wwmstat_0000'
         LFDB =len_trim(FDB)
         write(FDB(LFDB-3:LFDB),'(i4.4)') MYRANK
         open(STAT%FHNDL,file='outputs/'//fdb,status='replace') 
         FDB  ='windbg_0000'
         LFDB =len_trim(FDB)
         write(FDB(LFDB-3:LFDB),'(i4.4)') MYRANK
         open(WINDBG%FHNDL,file='outputs/'//fdb,status='replace') 
# else
         FDB  ='wwmdbg_0000'
         LFDB =len_trim(FDB)
         write(FDB(LFDB-3:LFDB),'(i4.4)') MYRANK
         open(DBG%FHNDL,file=fdb,status='replace') 
         FDB  ='wwmstat_0000'
         LFDB =len_trim(FDB)
         write(FDB(LFDB-3:LFDB),'(i4.4)') MYRANK
         open(STAT%FHNDL,file=fdb,status='replace') 
         FDB  ='windbg_0000'
         LFDB =len_trim(FDB)
         write(FDB(LFDB-3:LFDB),'(i4.4)') MYRANK
         open(WINDBG%FHNDL,file=fdb,status='replace')
# endif
#endif


         WRITE(DBG%FHNDL, *) 'THR=', THR
         WRITE(DBG%FHNDL, *) 'THR8=', THR8
         CALL TEST_FILE_EXIST_DIE("Missing input file : ", TRIM(INP%FNAME))
#ifdef MPIP_PARALL_GRID
         IF (myrank == 0) THEN
#endif
         WRITE(STAT%FHNDL,*) 'Input Filename   =', TRIM(INP%FNAME)
         WRITE(STAT%FHNDL,*) 'Check Filename   =', TRIM(CHK%FNAME)
         WRITE(STAT%FHNDL,*) 'Qstea Filename   =', TRIM(QSTEA%FNAME)
#if defined DEBUG && defined IOBPDOUT
         WRITE(STAT%FHNDL,*) 'Iobp Filename    =', TRIM(IOBPOUT%FNAME)
         WRITE(STAT%FHNDL,*) 'Iobpd Filename   =', TRIM(IOBPDOUT%FNAME)
#endif
         WRITE(STAT%FHNDL,*) 'WindDbg Filename =', TRIM(WINDBG%FNAME)
#ifdef MPIP_PARALL_GRID
         ENDIF
#endif
         CALL TEST_FILE_EXIST_DIE("Missing input file : ", INP%FNAME)
         OPEN( INP%FHNDL,      FILE = TRIM(INP%FNAME))
         OPEN( CHK%FHNDL,      FILE = TRIM(CHK%FNAME))
         IF (LQSTEA) OPEN( QSTEA%FHNDL,    FILE = TRIM(QSTEA%FNAME))
#if defined DEBUG && defined IOBPDOUT
         OPEN( IOBPOUT%FHNDL,  FILE = TRIM(IOBPOUT%FNAME))
         OPEN( IOBPDOUT%FHNDL, FILE = TRIM(IOBPDOUT%FNAME))
#endif

           OUT1D%FHNDL = STARTHNDL + 19 
            MISC%FHNDL = STARTHNDL + 20 
         OUTSP1D%FHNDL = STARTHNDL + 21 
         OUTPARM%FHNDL = STARTHNDL + 22 
         OUTSP2D%FHNDL = STARTHNDL + 23 

         OUT%FHNDL     = STARTHNDL + 24
         
       FHNDL_EXPORT_BOUC_WW3 = STARTHNDL + 25
       FHNDL_EXPORT_WIND_WW3 = STARTHNDL + 26
       FHNDL_EXPORT_CURR_WW3 = STARTHNDL + 27
       FHNDL_EXPORT_WALV_WW3 = STARTHNDL + 28
       FHNDL_EXPORT_GRID_WW3 = STARTHNDL + 29

         
       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE CLOSE_FILE_HANDLES
         USE DATAPOOL
         IMPLICIT NONE
         close(DBG%FHNDL)
         close(STAT%FHNDL)
         IF (LQSTEA) close( QSTEA%FHNDL)
#if defined DEBUG && defined IOBPDOUT
         close( IOBPOUT%FHNDL)
         close( IOBPDOUT%FHNDL)
#endif
         close( WINDBG%FHNDL)
       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_STATION_OUTPUT
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER           :: I, NI(3), IP, IS
      REAL(rkind)              :: XYTMP(2,MNP)
#ifdef MPI_PARALL_GRID
      integer :: iProc
      integer :: rbuf_int(1)
#endif
!
!    set the site output
!
      IF (LOUTS) THEN
        ALLOCATE (ACLOC_STATIONS(MSC,MDC,IOUTS), CDLOC_STATIONS(IOUTS), Z0LOC_STATIONS(IOUTS), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 38')

        ALLOCATE (ALPHALOC_STATIONS(IOUTS), WINDXLOC_STATIONS(IOUTS), WINDYLOC_STATIONS(IOUTS), USTARLOC_STATIONS(IOUTS), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 39')

        ALLOCATE (DEPLOC_STATIONS(IOUTS), WKLOC_STATIONS(IOUTS,MSC), CURTXYLOC_STATIONS(IOUTS,2), WATLEVLOC_STATIONS(IOUTS), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 40')
        ACLOC_STATIONS      = 0.
        DEPLOC_STATIONS     = 0.
        WKLOC_STATIONS      = 0.
        CURTXYLOC_STATIONS  = 0.
        USTARLOC_STATIONS   = 0.
        CDLOC_STATIONS      = 0.
        Z0LOC_STATIONS      = 0.
        WINDXLOC_STATIONS   = 0.
        WINDYLOC_STATIONS   = 0.
        WATLEVLOC_STATIONS  = 0.
      END IF
      IF (LOUTS .and. (DIMMODE .EQ. 2)) THEN
#ifdef MPI_PARALL_GRID
        XYTMP(1,:) = XP
        XYTMP(2,:) = YP
        WRITE(DBG%FHNDL,*) 'SEARCHING FOR STATION ACROSS RANKS', myrank
        DO I = 1, IOUTS
          STATION(I)%ELEMENT=0
          CALL FIND_ELE ( MNE,MNP,INE,XYTMP,STATION(I)%XCOORD, STATION(I)%YCOORD,STATION(I)%ELEMENT )
          IF (STATION(I)%ELEMENT .GT. 0) THEN
            STATION(I)%IFOUND  = 1
            NI                 = INE(:,STATION(I)%ELEMENT)
            STATION(I)%XELE(:) = XP(NI)
            STATION(I)%YELE(:) = YP(NI)
            CALL INTELEMENT_COEF(XP(NI),YP(NI), STATION(I)%XCOORD,STATION(I)%YCOORD, STATION(I)%WI)
            !WRITE(DBG%FHNDL,'(A10,I10,A20,I10,A15,2I10)') 'MYRANK', MYRANK, 'STATION =',I, 'IN ELEMENT =', IELG(STATION(I)%ELEMENT), STATION(I)%IFOUND
            WRITE(DBG%FHNDL,'(A10,I10,A20,I10,A15,2I10)') 'MYRANK', MYRANK, 'STATION =',I, 'IN ELEMENT =',STATION(I)%IFOUND
            FLUSH(DBG%FHNDL)
          ELSE
            STATION(I)%IFOUND  = 0
            NI                 = 0
            STATION(I)%XELE(:) = 0.
            STATION(I)%YELE(:) = 0.
            WRITE(DBG%FHNDL,'(A10,I10,A20,I10,A15,2I10)') 'MYRANK', MYRANK, 'STATION =',I, 'IN ELEMENT =', STATION(I)%ELEMENT, STATION(I)%IFOUND
            FLUSH(DBG%FHNDL)
          END IF
        END DO
        DO I = 1, IOUTS
          CALL MPI_REDUCE(STATION(I)%IFOUND,STATION(I)%ISUM,1, itype,MPI_SUM,0,COMM,IERR)
          IF (myrank == 0) THEN
            WRITE(DBG%FHNDL,'(A30,3I10)') 'SUM OF THE FOUND STATIONS MYRANK', MYRANK, I, STATION(I)%ISUM
            FLUSH(DBG%FHNDL)
            rbuf_int(1)=STATION(I)%ISUM
            DO iProc=2,nproc
              CALL MPI_SEND(rbuf_int,1,itype, iProc-1, 144, COMM, ierr)
            ENDDO
          ELSE
            CALL MPI_RECV(rbuf_int,1,itype, 0, 144, COMM, istatus, ierr)
            STATION(I)%ISUM=rbuf_int(1)
          END IF
        END DO
        IF (myrank == 0) THEN
          DO I = 1, IOUTS
            IF (STATION(I)%ISUM .EQ. 0) THEN
              WRITE(DBG%FHNDL,'(A20,I10,A10,2F15.8)') 'STATION NOT FOUND', I, STATION(I)%NAME, STATION(I)%XCOORD, STATION(I)%YCOORD
            ELSE
              WRITE(DBG%FHNDL,'(A25,I10,A10,2F15.8)') 'STATION FOUND    ', I, STATION(I)%NAME, STATION(I)%XCOORD, STATION(I)%YCOORD
            END IF
          END DO
        END IF
        ALLOCATE (DEPLOC_SUM(IOUTS), WKLOC_SUM(IOUTS,MSC), CURTXYLOC_SUM(IOUTS,2), ACLOC_SUM(MSC,MDC,IOUTS), USTAR_SUM(IOUTS), ALPHA_SUM(IOUTS), WINDY_SUM(IOUTS), WINDX_SUM(IOUTS), Z0_SUM(IOUTS), CD_SUM(IOUTS), WATLEVLOC_SUM(IOUTS), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 42')
        ACLOC_SUM           = 0.
        WKLOC_SUM           = 0.
        DEPLOC_SUM          = 0.
        CURTXYLOC_SUM       = 0.
        USTAR_SUM           = 0.
        CD_SUM              = 0.
        Z0_SUM              = 0.
        WINDX_SUM           = 0.
        WINDY_SUM           = 0.
        WATLEVLOC_SUM = 0.
#else
        XYTMP(1,:) = XP
        XYTMP(2,:) = YP
        IF (DIMMODE .EQ. 2) THEN
          WRITE(STAT%FHNDL,*) 'FINDING ELEMENT CONNECTED TO STATION'
          DO I = 1, IOUTS
            STATION(I)%ELEMENT=0
            CALL FIND_ELE ( MNE,MNP,INE,XYTMP,STATION(I)%XCOORD, STATION(I)%YCOORD,STATION(I)%ELEMENT )
            IF (STATION(I)%ELEMENT == 0) THEN
              STATION(I)%IFOUND = 0
              WRITE(STAT%FHNDL,*) STATION(I)%NAME, STATION(I)%XCOORD, STATION(I)%YCOORD, ' is out of mesh !'
            ELSE
              STATION(I)%IFOUND = 1
              NI = INE(:,STATION(I)%ELEMENT)
              STATION(I)%XELE(:) = XP(NI)
              STATION(I)%YELE(:) = YP(NI)
              CALL INTELEMENT_COEF(XP(NI),YP(NI), STATION(I)%XCOORD,STATION(I)%YCOORD, STATION(I)%WI)
              WRITE(STAT%FHNDL,*)'Site    ',STATION(I)%NAME, STATION(I)%XCOORD,STATION(I)%YCOORD,STATION(I)%IFOUND
            END IF
          END DO
        END IF
#endif
        STATION(:)%ISMAX = MSC
        IF (LSIGMAX) THEN
          DO IP = 1, IOUTS
            DO IS = 1, MSC
              IF (SPSIG(IS)/PI2 .GT. STATION(IP)%CUTOFF) THEN
                STATION(IP)%ISMAX = IS - 1
                EXIT
              END IF
            END DO
            WRITE(DBG%FHNDL,*) 'CUT-OFF FREQ. OF STATION =', IP, STATION(IP)%CUTOFF, 'RAD - IS =', STATION(IP)%ISMAX
          END DO
        END IF
      END IF
      WRITE(STAT%FHNDL,'("+TRACE...",A)')'FINISHED WITH INIT_STATION_OUTPUT'
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE TERMINATE_STATION_OUTPUT
      USE DATAPOOL
      IMPLICIT NONE
      IF (LOUTS) THEN
        DEALLOCATE (ACLOC_STATIONS)
        DEALLOCATE (CDLOC_STATIONS)
        DEALLOCATE (Z0LOC_STATIONS)
        DEALLOCATE (ALPHALOC_STATIONS)
        DEALLOCATE (WINDXLOC_STATIONS)
        DEALLOCATE (WINDYLOC_STATIONS)
        DEALLOCATE (USTARLOC_STATIONS)
        DEALLOCATE (DEPLOC_STATIONS)
        DEALLOCATE (WKLOC_STATIONS)
        DEALLOCATE (CURTXYLOC_STATIONS)
        DEALLOCATE (WATLEVLOC_STATIONS)
#ifdef MPI_PARALL_GRID
        DEALLOCATE (DEPLOC_SUM)
        DEALLOCATE (WKLOC_SUM)
        DEALLOCATE (CURTXYLOC_SUM)
        DEALLOCATE (ACLOC_SUM)
        DEALLOCATE (USTAR_SUM, ALPHA_SUM)
        DEALLOCATE (WINDY_SUM, WINDX_SUM)
        DEALLOCATE (Z0_SUM, CD_SUM, WATLEVLOC_SUM)
#endif
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READWAVEPARWWM
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER         :: IP
#ifdef MPI_PARALL_GRID
      INTEGER         :: IPP
      REAL(rkind)     :: RTMP
#endif

!     SPPARM(1): Hs, sign. wave height
!     SPPARM(2): Wave period given by user (either peak or mean)
!     SPPARM(3): average direction
!     SPPARM(4): directional spread
!     SPPARM(5): spectral shape (1-4),
!                (1 - Pierson-Moskowitz,
!                 2 - JONSWAP,
!                 3 - BIN,
!                 4 - Gauss)
!                     negative peak (+) 
!                     or mean frequency (-)
!     SPPARM(6): directional spreading in degree (1) or exponent (2)
!     SPPARM(7): gaussian width for the gauss spectrum 0.1
!     SPPARM(8): peak enhancement factor for the JONSWAP spectra 3.3

      IF (LINHOM) THEN
        READ(WAV%FHNDL,*)
      END IF

#ifdef MPI_PARALL_GRID
      IPP = 0
      IF (LINHOM) THEN
        DO IP = 1, IWBMNPGL
          IF(ipgl(IWBNDGL(IP))%rank == myrank) THEN ! if boundary nodes belong to local domain ...
            IPP = IPP + 1
            READ (WAV%FHNDL, *) SPPARM(:,IPP) ! ... read values into boundary array
          ELSE
            READ (WAV%FHNDL, *) RTMP, RTMP, RTMP, RTMP, RTMP, RTMP, RTMP, RTMP ! ... else ... throw them away
          ENDIF
        END DO
      ELSE
        READ (WAV%FHNDL, *) SPPARM(:,1)
        DO IP = 2, IWBMNPGL
          SPPARM(:,IP) = SPPARM(:,1)
        END DO
      END IF
#else 
      IF (LINHOM) THEN
        DO IP = 1, IWBMNP
          READ (WAV%FHNDL, *) SPPARM(:,IP)
        END DO
      ELSE
        READ (WAV%FHNDL, *) SPPARM(:,1)
      END IF
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_HMAX
        USE DATAPOOL

        HMAX = BRHD * DEP

        IF (LMONO_IN) HMAX = HMAX * SQRT(TWO)

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_WWMINPULNML
      USE DATAPOOL, only : INP
#ifdef ROMS_WWM_PGMCL_COUPLING
      USE mod_coupler, only : Iwaves, INPname
#endif
#if defined MODEL_COUPLING_ATM_WAV || defined MODEL_COUPLING_OCN_WAV
      USE coupling_var, only : WWM_InputFile
#endif
      IMPLICIT NONE
      INTEGER nbArg
#ifdef SCHISM
      INP%FNAME  = 'wwminput.nml'
#else
# if !defined ROMS_WWM_PGMCL_COUPLING && !defined MODEL_COUPLING_ATM_WAV && !defined MODEL_COUPLING_OCN_WAV
      nbArg=command_argument_count()
      IF (nbArg > 1) THEN
        CALL WWM_ABORT('Number of argument is 0 or 1')
      ENDIF
      IF (nbArg.eq.0) THEN
        INP%FNAME  = 'wwminput.nml'
      ELSE
        CALL GET_COMMAND_ARGUMENT(1, INP%FNAME)
      ENDIF
# else
#  ifdef ROMS_WWM_PGMCL_COUPLING
      INP%FNAME=INPname(Iwaves)
#  endif
#  if defined MODEL_COUPLING_ATM_WAV || defined MODEL_COUPLING_OCN_WAV
      INP%FNAME=WWM_InputFile
#  endif
# endif
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
