#ifdef VISDISLIN 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_DISLIN()

        USE DISLIN
        IMPLICIT NONE

        CALL METAFL('XWIN')
!        CALL WINDOW (0,0,640,480)
!        CALL WINSIZ(600,600)
        CALL X11MOD('STORE')
        CALL SETPAG('DA4L')
        CALL DISINI()

        RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DISPLAY_GRAPH(X,Y,N,XMIN,XMAX,YMIN,YMAX,TIT,TIT_X_AXIX,TIT_Y_AXIX)

         USE DISLIN
         IMPLICIT NONE

         CHARACTER (LEN = 100), INTENT(IN) :: TIT
         CHARACTER (LEN = 100), INTENT(IN) :: TIT_X_AXIX
         CHARACTER (LEN = 100), INTENT(IN) :: TIT_Y_AXIX

         INTEGER, INTENT(IN) :: N

         REAL(rkind), INTENT(IN)    :: X(N)
         REAL(rkind), INTENT(IN)    :: Y(N)

         REAL(rkind), INTENT(INOUT)   :: XMIN,XMAX
         REAL(rkind), INTENT(INOUT)   :: YMIN,YMAX

         REAL(rkind) :: YSTEP,XSTEP

         YSTEP = (YMAX-YMIN)/10.
         XSTEP = (XMAX-XMIN)/10.

         CALL ERASE()
         CALL PAGERA()
         CALL CENTER()
         CALL CHNCRV('BOTH')
         CALL COMPLX()
         CALL TICKS(5,'X')
         CALL TICKS(5,'Y')
         CALL LABELS('EXP','X')
         CALL LABELS('EXP','Y')
         CALL NAME(TIT_X_AXIX,'X')
         CALL NAME(TIT_Y_AXIX,'Y')
         CALL TITLIN(TIT,1)
         CALL GRAF(XMIN,XMAX,XMIN,XSTEP,YMIN,YMAX,YMIN,YSTEP)
         CALL COLOR('RED')
         CALL CURVE(X,Y,N)
         CALL TITLE()
         CALL ENDGRF()

       RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PLOT_FEM_SURFACE(XCOOR,YCOOR,ZCOOR,NTOTAL,I1RAY,I2RAY,I3RAY,NELMNT)

      USE DISLIN
      IMPLICIT NONE

      INTEGER, PARAMETER    :: DX = 5
      INTEGER, PARAMETER    :: DY = 5
      INTEGER, PARAMETER    :: DZ = 10
      INTEGER               :: AX_X, AX_Y, AX_Z
      INTEGER               :: WINSIZE_X, WINSIZE_Y, NTOTAL, NELMNT
      INTEGER               :: I1RAY(NELMNT), I2RAY(NELMNT), I3RAY(NELMNT)
      REAL(rkind)           :: XCOOR(NTOTAL), YCOOR(NTOTAL), ZCOOR(NTOTAL)
      REAL(rkind)           :: XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, DXIN, DYIN, DZIN, XYRATIO

      CALL ERASE()
      CALL HWFONT()
      CALL TITLIN('3-D Colour Plot of the Function',2)
      CALL TITLIN('Nordsee',4)
      CALL NAME('X-axis','X')
      CALL NAME('Y-axis','Y')
      CALL NAME('Z-axis','Z')
      CALL SHDMOD('SMOOTH', 'SURFACE')
      CALL INTAX()
      CALL AUTRES(100,100)
      CALL AXSPOS(300,1850)

      XMIN = MINVAL(XCOOR(:))
      XMAX = MAXVAL(XCOOR(:))
      YMIN = MINVAL(YCOOR(:))
      YMAX = MAXVAL(YCOOR(:))
      ZMIN = MINVAL(ZCOOR(:))
      ZMAX = MAXVAL(ZCOOR(:))

      DXIN = (XMAX - XMIN) / DX
      DYIN = (XMAX - XMIN) / DY
      DZIN = (ZMAX - ZMIN) / DZ

      WRITE (*,*) 'X-AXIS   ',XMIN, XMAX, DXIN
      WRITE (*,*) 'Y-AXIS   ',YMIN, YMAX, DYIN
      WRITE (*,*) 'Z-AXIS   ',ZMIN, ZMAX, DZIN

      XYRATIO = (MAXVAL(XCOOR(:))-MINVAL(XCOOR(:)))/(MAXVAL(YCOOR(:))-MINVAL(YCOOR(:)))

      AX_X = INT(1500 * XYRATIO)
      AX_Y = INT(1500)
      AX_Z = INT(MAXVAL(ZCOOR(:))-MINVAL(ZCOOR(:)))

      CALL AX3LEN(AX_X, AX_Y, AX_Z)
      CALL GRAF3(XMIN,XMAX,XMIN,DXIN,YMIN,YMAX,YMIN,DYIN,ZMIN,ZMAX,ZMIN,DZIN)
      CALL CRVTRI(XCOOR,YCOOR,ZCOOR,NTOTAL,I1RAY,I2RAY,I3RAY,NELMNT)
      CALL HEIGHT(50)
      CALL PSFONT('Palatino-BoldItalic')
      CALL TITLE()

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PLOT_SHADED_CONTOUR_GRAF(XRAY,YRAY,M,N,ZMAT,NLV,NXSTEP,NYSTEP,STR1)
      USE DISLIN
      IMPLICIT NONE

      REAL(rkind), INTENT(IN)      :: XRAY(M),YRAY(N)
      REAL(rkind), INTENT(IN)      :: ZMAT(M,N)
      INTEGER, INTENT(IN)   :: N,M,NLV

      INTEGER, INTENT(IN)         :: NXSTEP, NYSTEP
      CHARACTER(len=5),INTENT(IN) :: STR1

      INTEGER               :: I

      REAL(rkind)                  :: ZLEV(NLV)
      REAL(rkind)                  :: XSTEP,YSTEP,ZSTEP
      REAL(rkind)                  :: XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX

      CHARACTER(len=10)     :: str2 = 'IP ='

      CALL ERASE()
      CALL PAGERA()
      CALL COMPLX()
      CALL MIXALF()
!      CALL TITLIN('Shaded Contour Plot',1)
!      CALL TITLIN(str1,1)
      CALL HEIGHT(30)
      CALL NAME('X-axis - IP ='//str1,'X')
      CALL NAME('Y-axis','Y')
      CALL AXENDS('NOENDS','X')
!     CALL AXSSCL('LOG','X')
      CALL AXSLEN(1500,1500)
      CALL AXSORG(1500,1050)

      XMIN = MINVAL(XRAY)
      XMAX = MAXVAL(XRAY)
      YMIN = MINVAL(YRAY)
      YMAX = MAXVAL(YRAY)
      ZMIN = MINVAL(ZMAT)
      ZMAX = MAXVAL(ZMAT)
      WRITE (*,*) 'X',XMIN, XMAX
      WRITE (*,*) 'Y',YMIN, YMAX
      WRITE (*,*) 'Z',ZMIN, ZMAX
      XSTEP = (XMAX-XMIN)/NXSTEP
      YSTEP = (YMAX-YMIN)/NYSTEP
      ZLEV(1) = ZMIN
      ZSTEP = (ZMAX-ZMIN)/NLV
      DO I = 2, NLV
        ZLEV(I) = ZLEV(I-1) + ZSTEP
      END DO
      WRITE (*,*) XSTEP, YSTEP, ZSTEP
      CALL SHDMOD('POLY','CONTUR')
      CALL GRAF(XMIN,XMAX,XMIN,XSTEP,YMIN,YMAX,YMIN,YSTEP)
      CALL CONSHD(XRAY, M, YRAY, N, ZMAT, ZLEV, NLV)
      CALL HEIGHT(28)
      CALL TITLE()
      CALL ENDGRF()
      RETURN
      END SUBROUTINE PLOT_SHADED_CONTOUR_GRAF
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PLOT_SHADED_CONTOUR_POLAR(XRAY,YRAY,M,N,ZMAT,NLV,NXSTEP,NYSTEP,STR1)
      USE DISLIN
      IMPLICIT NONE

      REAL(rkind), INTENT(IN)      :: XRAY(M),YRAY(N)
      REAL(rkind), INTENT(IN)      :: ZMAT(M,N)
      INTEGER, INTENT(IN)          :: N,M,NLV

      INTEGER, INTENT(IN)   :: NXSTEP, NYSTEP

      CHARACTER(len=5),INTENT(IN) :: STR1

      INTEGER               :: I

      REAL(rkind)                  :: ZLEV(NLV)
      REAL(rkind)                  :: XSTEP,YSTEP,ZSTEP
      REAL(rkind)                  :: XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX

      REAL(rkind), PARAMETER       :: PI       = 3.141593
      REAL(rkind), PARAMETER       :: PI2      = 2*PI
      REAL(rkind), PARAMETER       :: DEGRAD   = PI/180.0
      REAL(rkind), PARAMETER       :: RADDEG   = 180.0/PI

      CHARACTER(len=10)     :: str2 = 'IP ='

      CALL ERASE()
      CALL PAGERA()
      CALL COMPLX()
      CALL MIXALF()

      CALL HEIGHT(28)
      CALL TICKS(2,'X')
      CALL TICKS(10,'Y')
      CALL TITLIN('Shaded Contour Plot',1)
      CALL TITLIN(str1,2)
      CALL NAME('X-axis - IP ='//str1,'X')
      CALL NAME('Y-axis','Y')
      CALL AXENDS('NOENDS','X')

      CALL AXSLEN(1500,1500)
      CALL AXSORG(1500,1050)
 
      XMIN = MINVAL(XRAY)
      XMAX = MAXVAL(XRAY)
      YMIN = MINVAL(YRAY)
      YMAX = MAXVAL(YRAY)
      ZMIN = MINVAL(ZMAT)
      ZMAX = MAXVAL(ZMAT)

      XSTEP = (XMAX-XMIN)/NXSTEP
      YSTEP = (360.)/NYSTEP
      ZLEV(1) = ZMIN
      ZSTEP = (ZMAX-ZMIN)/NLV
      DO I = 2, NLV
        ZLEV(I) = ZLEV(I-1) + ZSTEP
      END DO
 
      WRITE(*,*) ZLEV 

# ifndef SCHISM
      WRITE(*,*) MAXVAL(ZLEV), MAXVAL(XRAY), MAXVAL(YRAY)
      WRITE(*,*) MINVAL(ZLEV), MINVAL(XRAY), MINVAL(YRAY)
# endif

      CALL HEIGHT(28)
      CALL SHDMOD('POLY','CONTUR')
      CALL POLAR (XMAX/PI2, 0., XSTEP/PI2, 0., YSTEP)
      CALL CONSHD(XRAY/PI2, M,  YRAY*DEGRAD, N, ZMAT, ZLEV, NLV)
      CALL TITLE()
      CALL ENDGRF()

      RETURN
      END SUBROUTINE PLOT_SHADED_CONTOUR_POLAR
!**********************************************************************
!*                                                                    *
!**********************************************************************
#endif
