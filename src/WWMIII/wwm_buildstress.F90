#ifdef WAM_ECMWF
      SUBROUTINE BUILDSTRESS(U10OLD,THWOLD,USOLD,TAUW,Z0OLD, &
     &                       ROAIRO, ZIDLOLD, ICEMASK, &
     &                       IREAD)
#endif
      SUBROUTINE BUILDSTRESS

! ----------------------------------------------------------------------
!     J. BIDLOT    ECMWF   APRIL 1998 

!     J. BIDLOT    ECMWF   FEBRUARY 1999 TAUT --> SQRT(TAUT)

!     S. ABDALLA   ECMWF   OCTOBER 1999 MODIFICATION THE CALL TO GETWND
 
!     J. BIDLOT    ECMWF   AUGUST 2008 : MAKE IT MORE PARALLEL.

!*    PURPOSE.
!     --------
!     CREATES WIND AND STRESS FIELDS FROM GRIB WINDS AND CD.

!**   INTERFACE.
!     ----------
!     CALL *BUILDSTRESS*(U10OLD,THWOLD,USOLD,TAUW,Z0OLD,ROAIRO,
!    &                   ROAIRO, ZIDLOLD, ICEMASK,
!    &                   IREAD)*
!     *U10OLD*   WIND SPEED.
!     *THWOLD*   WIND DIRECTION (RADIANS).
!     *USOLD*    FRICTION VELOCITY.
!     *TAUW*     WAVE STRESS.
!     *Z0OLD*    ROUGHNESS LENGTH IN M.
!     *RAD0OLD*   AIR DENSITY IN KG/M3.
!     *RZIDL0OLD* Zi/L (Zi: INVERSION HEIGHT, L: MONIN-OBUKHOV LENGTH).
!     *ICEMASK*   SEA ICE MASK
!     *IREAD*     PROCESSOR WHICH WILL ACCESS THE FILE ON DISK

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------
!     *ABORT1*
!     *AIRSEA*
!     *GETWND*
!     *PBOPEN*
!     *PBREAD*
!     *PBCLOSE*
!     *READWGRIB*

!     REFERENCE.
!     ----------
!     NONE
! ----------------------------------------------------------------------

!      USE YOWCOUP  , ONLY : LWCOU    ,ALPHA    ,XKAPPA   ,XNLEV
!      USE YOWGRIBHD, ONLY : NKSEK1 
!      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL
!      USE YOWMPP   , ONLY : NINF     ,NSUP
!      USE YOWMESPAS, ONLY : LMESSPASS,LNOCDIN  ,LWAVEWIND
!      USE YOWPARAM , ONLY : NBLO     ,NIBLO    ,NGX      ,NGY
!      USE YOWPCONS , ONLY : G        ,ROAIR    ,EPSUS    ,EPSU10
!      USE YOWSTAT  , ONLY : CDATEA   ,CDTPRO   ,NWAM_BLKS
!      USE YOWTABL  , ONLY : ITAUMAX  ,JUMAX    ,JPLEVT   ,TAUT
!      USE YOWTEST  , ONLY : IU06     ,ITEST
!      USE YOWWIND  , ONLY : CDAWIFL  ,CDATEWO  ,CDATEFL  ,FIELDG

!      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE DATAPOOL, ONLY  : MNP, EPSU10, CD, U10OLD, USOLD, Z0OLD, ILEV, EPSUS, TAUW, WINDXY
      USE DATAPOOL, ONLY  : RKIND, ITEST, IU06
      IMPLICIT NONE

! ----------------------------------------------------------------------
      INTEGER     :: IREAD, IP
      REAL(rkind) :: CDINV
! ----------------------------------------------------------------------

!     1.3 INITIALISE CD USING THE FRICTION VELOCITY FOR TAUW=0.
!         ----------------------------------------------------

      DO IP = 1, MNP
        TAUW(IP)=0.
        U10OLD(IP,1) = SQRT(WINDXY(IP,1)**2.+WINDXY(IP,2)**2.) 
        CALL AIRSEA (U10OLD(IP,1),TAUW(IP),USOLD(IP,1), Z0OLD(IP,1), IP, IP, ILEV)
        CDINV = MAX(U10OLD(IP,1)**2,EPSU10)/MAX(USOLD(IP,1)**2,EPSUS)
        CDINV = MIN(CDINV,10000.0_rkind) 
        CD(IP) = 1./CDINV
        IF (ITEST.GT.0) WRITE (IU06,*) ' SUB. AIRSEA DONE AT 1'
      ENDDO

      IF (ITEST.GE.1) WRITE(IU06,*) ' SUB. BUILDSTRESS: INPUT OF RESTART FILES DONE'

      END SUBROUTINE BUILDSTRESS
! ----------------------------------------------------------------------
