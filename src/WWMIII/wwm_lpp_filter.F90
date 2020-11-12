!**********************************************************************
!*  This routine is used with LPP_FILT_FLAG = T; Use LLP_FRAC in wwminput.nml to choose the pourcentage of wave length. 
!*  => Filter the non-conservative terms due to depth-induced breaking (term Fb from Eq. (11) and (12))
!**********************************************************************

      SUBROUTINE LPP_FILT(VAR)

        USE DATAPOOL,    ONLY : AC2,MNP,MSC,rkind,LPP_FRAC
        USE schism_glbl, ONLY : idry,xnd,ynd
        USE schism_msgp, ONLY : exchange_p2d 

        IMPLICIT NONE
  
        INTEGER                    :: IP,IPP,icount
        REAL(rkind),INTENT(INOUT)  :: VAR(MNP)
        REAL(rkind),DIMENSION(MNP) :: dist,dx,dy
        REAL(rkind)                :: var_tmp
        REAL(rkind)                :: FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,   &
                                    & PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD, &
                                    & CGPD,CPPD
        DO IP=1,MNP

          IF(idry(IP)==1) CYCLE
  
          dx=ABS(xnd(IP)-xnd(:))
          dy=ABS(ynd(IP)-ynd(:))
          dist=SQRT(dx**2.0d0 + dy **2.0d0)
  
          call PEAK_PARAMETER(IP,AC2(:,:,IP),MSC,FPP,TPP,CPP,WNPP,CGPP,&
                    & KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD)
  
          icount=0
          var_tmp=0.0d0       
          DO IPP=1,MNP
            IF (idry(IPP)==1) CYCLE
            ! LPP_FRAC: fraction of peak wave length read in wwminput.nml 
            IF (dist(IPP) .LT. LPP_FRAC*LPP) THEN
              var_tmp=var_tmp+VAR(IPP)
              icount=icount+1
            END IF
          END DO

          IF (icount > 0) THEN
            VAR(IP)=var_tmp/icount
          END IF

        END DO
 
        CALL exchange_p2d(VAR)       
		

      END SUBROUTINE LPP_FILT

