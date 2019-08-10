      REAL FUNCTION TRANSF(XK,D)
!
!***  DETERMINE NARROW BAND LIMIT NONLINEAR TRANSFER FUNCTION            
!     BASED ON TECH MEMO 464 BY P. JANSSEN AND M. ONORATO
!     
!
!     AUTHOR:  P.A.E.M. JANSSEN ECMWF JUNE 2005
!     ------
!
!     VARIABLE       TYPE         PURPOSE
!     --------       ----         -------
!
!     XK             REAL         WAVE NUMBER
!     D              REAL         DEPTH
!
!----------------------------------------------------------------------

!      USE YOWPCONS , ONLY : G     ,DKMAX

       USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, FRINTF, COFRM4, CG, WK, DKMAX, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, RKIND, EMEAN, FMEAN, TH, ENH, DEP, AF11, &
     &                      IKP, IKP1, IKM, IKM1, K1W, K2W, K11W, K21W, FKLAP, FKLAP1, FKLAM, FKLAM1, FRH, &
     &                      DELTH => DDIR, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NANG => MDC, &
     &                      NFRE => MSC, &
     &                      INDEP => DEP

!----------------------------------------------------------------------
 
      IMPLICIT NONE

      REAL :: EPS,D,X,XK,T_0,OM,C_0,V_G,DV_G,XNL_1,XNL_2,XNL

      EPS=0.0001
!
!*    1. DETERMINE TRANSFER FUNCTION.
!     ------------------------------
!     
      IF(D.LT.999. .AND. D.GT.0.) THEN
        X   = XK*D
        IF ( X .GT. DKMAX) THEN
          TRANSF = 1. 
        ELSE
          T_0 = TANH(X)
          OM  = SQRT(G*XK*T_0)
          C_0 = OM/XK
          IF(X .LT. EPS) THEN
            V_G = 0.5*C_0
            V_G = C_0
          ELSE
            V_G = 0.5*C_0*(1.+2.*X/SINH(2.*X))
          ENDIF
          DV_G = (T_0-X*(1.-T_0**2))**2+4.*X**2*T_0**2*(1.-T_0**2)
       
          XNL_1 = (9.*T_0**4-10.*T_0**2+9.)/(8.*T_0**3)
          XNL_2 = ((2.*V_G-0.5*C_0)**2/(G*D-V_G**2)+1.)/X

          XNL = XNL_1-XNL_2
          TRANSF = XNL**2/(DV_G*T_0**8)
        ENDIF
      ELSE
        TRANSF = 1. 
      ENDIF
!
      RETURN
      END
