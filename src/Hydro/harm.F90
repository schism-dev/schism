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

!******************************************************************************
! PADCIRC VERSION 45.11 02/02/2006                                            *
!    last changes in this file prior to VERSION 44.15                         *
!                                                                             *
!                                                                             *
!   added local LOGICAL CHARMV declaration 04/21/2004                         *
!******************************************************************************
!                                                                      *
!   PADCIRC MODULE  ( HARM )                                           *
!                                                                      *
!   HA_SUBS.FOR     V3.01        11/9/95                               *
!                                                                      *
!   Least Square harmonic analysis of timeseries from ADCIRC2DDI_v27   *
!                                                                      *
!    Notes:                                                            *
!    1.)  Both the left hand side matrix and the right hand side       *
!         forcing vectors are continuously updated in time.  This      *
!         eliminates the need to store time series outputs for later   *
!         harmonic analysis.                                           *
!    2.)  The left hand side matrix and the right hand side forcing    *
!         vectors are output in the hotstart file and can be used to   *
!         perform harmonic analysis on an incomplete run.              *
!    3.)  Frequencies should be in rad/sec,times should be in sec.     *
!                                                                      *
!***********************************************************************
!                                                                      *
!    Program Written by:                                               *
!          R.A. Luettich, IMS UNC                                      *
!          J.J. Westerink, CE ND                                       *
!                                                                      *
!    Program Development History:                                      *
!    1.) lsq_stations_v004 by JJW                                      *
!    2.) LSQEX by RL used in 2D binary extr program                    *
!    3.) LSQRL by RL used in 1D test codes                             *
!    4.) LSQ2D v1.00-v2.26 by RL real time Harmonic Analysis for ADCIRC*
!    5.) HA_SUBS v3.01 by RL real time HA for ADCIRC separate          *
!        subroutines for elevation station, velocity station,          *
!        global elevation and global velocity harmonic analysis        *
!                                                                      *
!***********************************************************************
!                                                                      *
! SUBROUTINE LSQUPDLHS updates the LHS matrix                          *
! SUBROUTINE LSQUPDEG updates the RHS load vector for elevation global *
! SUBROUTINE LSQUPDVG updates the RHS load vector for velocity global  *
! SUBROUTINE FULSOL fills out, decomposes and solves the matricies     *
! SUBROUTINE LSQSOLEG solves & writes output for elevation global      *
! SUBROUTINE LSQSOLVG solves & writes output for velocity global       *
! SUBROUTINE HAHOUT writes HA parameters & LHS matrix to hotstart file *
! SUBROUTINE HAHOUTEG writes glob elev RHS load vector to hotstart file*
! SUBROUTINE HAHOUTVG writes glob vel RHS load vector to hotstart file *
! SUBROUTINE HACOLDS initializes HA param & LHS matrix for cold start  *
! SUBROUTINE HACOLDSEG initializes glob ele RHS load vec for cold start*
! SUBROUTINE HACOLDSVG initializes glob vel RHS load vec for cold start*
! SUBROUTINE HAHOTS initializes HA params & LHS matrix for a hot start *
! SUBROUTINE HAHOTSEG initializes glob elev RHS load vec for hot start *
! SUBROUTINE HAHOTSVG initializes glob vel RHS load vec for hot start  *
!                                                                      *
!***********************************************************************
!                                                                      *
!    INPUT FILES:                                                      *
!      - Frequency information is read in by ADCIRC from unit 15.      *
!        This information is passed in common block LSQFREQS.          *
!                                                                      *
!      - If the model is hot start, input is read from UNIT 67 or 68   *
!                                                                      *
!    OUTPUT FILES:                                                     *
!      UNIT 53 : HARMONIC CONSTITUENT ELEVATIONS AT ALL NODES (ASCII)  *
!      UNIT 54 : HARMONIC CONSTITUENT VELOCITIES AT ALL NODES (ASCII)  *
!      UNIT 55 : COMPARISON BETWEEN THE MEAN AND VARIANCE OF THE TIME  *
!                  SERIES GENERATED BY THE MODEL AND THE MEAN AND      *
!                  VARIANCE OF A TIME SERIES RESYNTHESIZED FROM THE    *
!                  COMPUTED HARMONIC CONSTITUENTS.  THIS GIVES AN      *
!                  INDICATION OF HOW COMPLETE THE HARMONIC ANALYSIS    *
!                  WAS. (ASCII)                                        *
!      UNIT 67 or 68 : HOT START FILES (BINARY)                        *
!                                                                      *
!***********************************************************************
!
      MODULE HARM
      USE schism_glbl, only: rkind,pi,MNHARF,CHARMV,ipgl,in_dir,out_dir, &
     &len_in_dir,len_out_dir
!
!     REAL(8),PRIVATE,PARAMETER :: PI=3.141592653589793D0
      INTEGER NFREQ
      CHARACTER*10,ALLOCATABLE ::  NAMEFR(:)
      REAL(rkind),    ALLOCATABLE ::  HAFREQ(:),HAFF(:),HAFACE(:)
!
      INTEGER, PRIVATE,SAVE :: NZ, NF, MM, ITUD, ICALL
      REAL(8), PRIVATE,SAVE :: TIMEUD
      REAL(rkind),PRIVATE,ALLOCATABLE ::  HA(:,:)
      REAL(rkind),PRIVATE,ALLOCATABLE ::  HAP(:),HAX(:)
      REAL(rkind),PRIVATE,ALLOCATABLE ::  GLOELV(:,:)
      REAL(rkind),PRIVATE,ALLOCATABLE ::  GLOULV(:,:),GLOVLV(:,:)
      REAL(rkind),PRIVATE,ALLOCATABLE ::  STAELV(:,:)
      REAL(rkind),PRIVATE,ALLOCATABLE ::  STAULV(:,:),STAVLV(:,:)


!-----------------END OF DECLARATIONS---------------------------------------

      CONTAINS


!
!***********************************************************************
!  Allocate arays used by LSQ_HARM.
!
!  vjp 8/99
!***********************************************************************
!
      SUBROUTINE ALLOC_HA(NP)
      integer, intent(in) :: NP

      ALLOCATE ( HAFREQ(MNHARF),HAFF(MNHARF),HAFACE(MNHARF) )
      ALLOCATE ( NAMEFR(MNHARF) )
!     
      ALLOCATE ( HA(2*MNHARF,2*MNHARF) )
      ALLOCATE ( HAP(2*MNHARF),HAX(2*MNHARF) )
      ALLOCATE ( GLOELV(2*MNHARF,NP) )
      ALLOCATE ( GLOULV(2*MNHARF,NP),GLOVLV(2*MNHARF,NP) )
!      ALLOCATE ( GLOELV(2*MNHARF,NP) )
!      ALLOCATE ( GLOULV(2*MNHARF,NP),GLOVLV(2*MNHARF,NP) )
!      ALLOCATE ( STAELV(2*MNHARF,MNSTAE) )
!      ALLOCATE ( STAULV(2*MNHARF,MNSTAV),STAVLV(2*MNHARF,MNSTAV) )
      
      RETURN
      END SUBROUTINE
      

!***********************************************************************
!   Subroutine to update the Left Hand Side Matrix                     *
!                                                                      *
!  TIME  - ABSOLUTE MODEL TIME (SEC)                                   *
!  IT    - MODEL TIME STEP                                             *
!  icall - number of times the subroutine has been called              *
!  a     - Left Hand Side Matrix                                       *
!                                                                      *
!                        RL 11/7/95                                    *
!***********************************************************************
!
      SUBROUTINE LSQUPDLHS(TIME,IT)
      IMPLICIT NONE
      INTEGER IT,I,J,I1,I2,J1,J2
      REAL(rkind) TF1,TF2
      REAL(8) TIME
!
      icall = icall + 1
!     
!***** Update the Left Hand Side Matrix
!     Note: this is a symmetric matrix and therefore only store the
!     upper triangular part.  The lower part will be filled out in
!     SUBROUTINE FULSOL prior to the matrix's decomposition

!     Take care of the steady constituent if included in the analysis

      if(nf.eq.1) then
         ha(1,1)=icall
         do j=1,nfreq
            tf1=hafreq(j+nf)*time
            ha(1,2*j)   = ha(1,2*j) + cos(tf1)
            ha(1,2*j+1) = ha(1,2*j+1) + sin(tf1)
         end do
      endif

!   Take care of the other constituents

      do i=1,nfreq
         do j=i,nfreq
            i1=2*i-(1-nf)
            i2=i1+1
            j1=2*j-(1-nf)
            j2=j1+1
            tf1=hafreq(i+nf)*time
            tf2=hafreq(j+nf)*time
            ha(i1,j1) = ha(i1,j1) + cos(tf1)*cos(tf2)
            ha(i1,j2) = ha(i1,j2) + cos(tf1)*sin(tf2)
            ha(i2,j2) = ha(i2,j2) + sin(tf1)*sin(tf2)
            if(i2.le.j1) ha(i2,j1) = ha(i2,j1) + sin(tf1)*cos(tf2)
         end do
      end do

!   Record update time and time step

      TIMEUD = TIME
      ITUD = IT
      
      return
      end subroutine

!***********************************************************************
!   Subroutine to update the Right Hand Side Load Vectors for the      *
!   global elevation harmonic analysis.                                *
!                                                                      *
!  GLOE  - GLOBAL ELEVATION VALUES USED TO UPDATE LOAD VECTORS         *
!  NP    - NUMBER OF POINTS IN GLOBAL GRID                             *
!                                                                      *
!  GLOELV - global elevation load vector                               *
!                                                                      *
!                        RL 11/8/95                                    *
!***********************************************************************
!
      SUBROUTINE LSQUPDEG(GLOE,NP)
      IMPLICIT NONE
      INTEGER I,J,NP,N,I1,I2,IR,IRE,K,JR
      REAL(rkind) TF1,CTF1,STF1
      REAL(rkind) GLOE(NP)
!      REAL(rkind) GLOE(NP)
!     
!*****Update the Right Hand Side Load Vectors
!     
!     Take care of the steady constituent if included in the analysis

      if(nz.eq.0) then
         do n=1,np
            GLOELV(1,N)=GLOELV(1,N)+GLOE(N)
         end do
      endif

!     Take care of the other constituents

      do i=1,nfreq
         i1=2*i-nz
         i2=i1+1
         tf1=hafreq(i+nf)*TIMEUD
         ctf1 = cos(tf1)
         stf1 = sin(tf1)
         do n=1,np
            GLOELV(I1,N)=GLOELV(I1,N)+GLOE(N)*CTF1
            GLOELV(I2,N)=GLOELV(I2,N)+GLOE(N)*STF1
         end do
      end do
!     
      return
      end subroutine


!***********************************************************************
!   Subroutine to update the Right Hand Side Load Vectors for the      *
!   global velocity harmonic analysis.                                 *
!                                                                      *
!  GLOU  - GLOBAL U VELOCITY VALUES USED TO UPDATE LOAD VECTORS        *
!  GLOV  - GLOBAL V VELOCITY VALUES USED TO UPDATE LOAD VECTORS        *
!  NP    - NUMBER OF POINTS IN GLOBAL GRID                             *
!                                                                      *
!  GLOULV - global u velocity load vector                              *
!  GLOVLV - global v velocity load vector                              *
!                                                                      *
!                        RL 11/8/95                                    *
!***********************************************************************
!
      SUBROUTINE LSQUPDVG(GLOU,GLOV,NP)
      IMPLICIT NONE
      INTEGER NP,NPI1,I1,I2,N,I,J,IR,IRE,K,JR
      REAL(rkind) TF1,CTF1,STF1
      REAL(rkind) GLOU(NP),GLOV(NP)
!      REAL(rkind) GLOU(NP),GLOV(NP)
!     
!*****Update the Right Hand Side Load Vectors
!     
!     Take care of the steady constituent if included in the analysis

      if(nz.eq.0) then
         do n=1,np
            GLOULV(1,N) = GLOULV(1,N) + GLOU(N)
            GLOVLV(1,N) = GLOVLV(1,N) + GLOV(N)
         end do
      endif

!     Take care of the other constituents

      do i=1,nfreq
         i1=2*i-nz
         i2=i1+1
         tf1=hafreq(i+nf)*TIMEUD
         ctf1 = cos(tf1)
         stf1 = sin(tf1)
         do n=1,np
            GLOULV(I1,N) = GLOULV(I1,N) + GLOU(N)*CTF1
            GLOVLV(I1,N) = GLOVLV(I1,N) + GLOV(N)*CTF1
            GLOULV(I2,N) = GLOULV(I2,N) + GLOU(N)*STF1
            GLOVLV(I2,N) = GLOVLV(I2,N) + GLOV(N)*STF1
         end do
      end do
!     
      return
      end subroutine

!***********************************************************************
!   Subroutine to fill out, decompose and solve the lsq system         *
!   Solves system a*x=b by l*d*l(tr) decomp in full storage mode       *
!                                                                      *
!   NOTE: This routine has been modified so that the filling out and   *
!         decomposition (and only those operations) are done if        *
!         idecom=0.                                                    *
!                                                                      *
!   mm  -  actual dimension of a matrix                                *
!                                                                      *
!                        rl 11/7/95                                    *
!***********************************************************************
!
      subroutine fulsol(idecom)
      implicit none
      integer idecom,i,j,ir,ire,k,jr
      real(rkind),allocatable ::  c(:),y(:)

!     
!**** If only want to fill out matrix and decomponse
!     
      if(idecom.eq.0) then
         
!     Set up the lower triangular part of the LHS a matrix
         
         do j=1,mm
            do i=j,mm
               ha(i,j)=ha(j,i)
            end do
         end do
         
!     Decomposition of matrix a

         do 100 ir=1,mm
            ire=ir+1
            do 20 j=ire,mm
 20         ha(ir,j)=ha(ir,j)/ha(ir,ir)
            if(ire.gt.mm) goto 100
            do 40 j=ire,mm
              do 40 k=ire,mm
 40           ha(k,j)=ha(k,j)-ha(k,ir)*ha(ir,j)
            do 50 j=ire,mm
 50           ha(j,ir)=0.0
 100     continue
         return
      endif

!...  solve for y by forward substitution for l*y=p

      allocate ( c(2*MNHARF),y(2*MNHARF) )
!     
      do 120 ir=1,mm
         y(ir)=hap(ir)
         do 110 jr=1,ir-1
 110        y(ir)=y(ir)-ha(jr,ir)*y(jr)
 120     continue

!...  calculate c=d**(-1)*y

         do 130 ir=1,mm
 130        c(ir)=y(ir)/ha(ir,ir)

!...  solve for x by back-substituting for l(tr)*x=c

            ir=mm
 140        continue
            hax(ir)=c(ir)
            do 150 jr=ir+1,mm
 150          hax(ir)=hax(ir)-ha(ir,jr)*hax(jr)
            ir=ir-1
            if(ir.ge.1) goto 140
      deallocate(c,y)
      return
      end subroutine

!***********************************************************************
!   Subroutine to solve the system and write output for elevation      *
!   globally.                                                          *
!                                                                      *
!   nf=0  if no steady constituent                                     *
!   nf=1  if steady constituent                                        *
!                                                                      *
!                        R.L. 11/8/95                                  *
!                                                                      *
!   added local LOGICAL CHARMV declaration 04/09/2004                  *
!   closed unit 53 08/23/05                                            *
!   Modified for SELFE - Andre Fortunato 2009/06/08                    *
!***********************************************************************
!
      SUBROUTINE LSQSOLEG(NP,ELAV,ELVA,CHARMV, &
     &                    myrank,TIMEBEG,DT,FMV,NTSTEPS,ITMV)
      IMPLICIT NONE
      LOGICAL CHARMV            !rl 04/09/2004
      integer J,NP,N,K,I,I1,I2,IT,IFR,NEAVMAX,NEAVMIN, &
     &  NEVAMAX,NEVAMIN,myrank
      REAL(8)  CONVRD 
      REAL(rkind) EAVMAX,EVAMAX,EAVMIN,EVAMIN,EMAGT
      REAL(rkind) PHASEDE,EAV,ESQ,TIME,RSE,FTIME,EAVDIF,EVADIF
      REAL(rkind) ELAV(NP),ELVA(NP)
!      REAL(rkind) ELAV(NP),ELVA(NP)
      REAL(rkind),ALLOCATABLE  ::  PHASEE(:),EMAG(:)
!
      INTEGER NTSTEPS,ITMV
      REAL(8) TIMEBEG
      REAL(rkind) DT,FMV
!      COMMON /MEANSQ/ TIMEBEG,DT,FMV,NTSTEPS,ITMV
      CHARACTER*4 DIRNAME
!
      convrd=180.d0/pi
!
      ALLOCATE ( PHASEE(MNHARF),EMAG(MNHARF) )
!
!**** Open velocity station harmonic output file and write header information
!
      write(DIRNAME(1:4),'(i4.4)') myrank
      open(53,file=out_dir(1:len_out_dir)//'harme.53'//DIRNAME(1:4))
      write(53,*) nfreq+nf
      do j=1,nfreq+nf
         write(53,3679) hafreq(j),HAFF(j),HAFACE(j),namefr(j)
      end do
 3679 format(1x,e20.10,1x,f10.7,1x,f12.8,1x,a10)
      write(53,*) NP
!
      if (CHARMV) then
         EAVMAX=-999.
         EVAMAX=-999.
         EAVMIN= 999.
         EVAMIN= 999.
      end if
!
!***** AT each node transfer each load vector to p, solve and write output
!
      DO N=1,NP
         do k=1,mm
            hap(k) = GLOELV(k,n)
         end do
         call fulsol(n)
!
!        Compute amplitude and phase for each frequency making sure that the
!        phase is between 0 and 360 deg.  Then write output.
!
         write(53,*) N
         do i=1,nfreq+nf
            if((nf.eq.1).and.(i.eq.1)) then
               emag(i)=hax(i)
               emagt=emag(i)/haff(i)
               phasee(i)=0.
            else
               i1=2*i-1-nf
               i2=i1+1
               emag(i)=sqrt(hax(i1)*hax(i1)+hax(i2)*hax(i2))
               emagt=emag(i)/haff(i)
               if((hax(i1).eq.0.).and.(hax(i2).eq.0.)) then
                  phasee(i)=0.
               else
                  phasee(i) = atan2(hax(i2),hax(i1))
               endif
            endif
            phasede=convrd*phasee(i)+haface(i)
            if(phasede.lt.0.) phasede=phasede+360.d0
            if(phasede.ge.360.d0) phasede=phasede-360.d0
            write(53,6635) emagt,phasede
 6635       format(2x,e16.8,1x,f11.4)
         end do
         
         if (CHARMV) then
            eav = 0.
            esq = 0.
            do it=1,ntsteps
               TIME=TIMEBEG+DT*IT
               rse=0.
               do ifr=1,nfreq+nf
                  ftime=hafreq(ifr)*time
                  rse=rse+emag(ifr)*cos(ftime-phasee(ifr))
               end do
               eav=eav+rse
               esq=esq+rse*rse
            end do
            
         eav=eav/ntsteps
         esq=esq/ntsteps-eav*eav
         if(elav(n).eq.0.) then
            if(eav.eq.0.) eavdif=1.0d0
            if(eav.ne.0.) eavdif=99d19
         else
            eavdif=eav/elav(n)
         endif
         if(elva(n).eq.0.) then
            if(esq.eq.0.) evadif=1.0d0
            if(esq.ne.0.) evadif=99e19
         else
            evadif=esq/elva(n)
         endif
         write(55,*) n
         write(55,7637) elav(n),eav,eavdif,elva(n),esq,evadif
 7637    format(2x,3(e16.8,1x),2x,3(e16.8,1x))
         
         IF(EAVDIF.GT.EAVMAX) THEN
            EAVMAX=EAVDIF
            NEAVMAX=n
         ENDIF
         IF(EAVDIF.LT.EAVMIN) THEN
            EAVMIN=EAVDIF
            NEAVMIN=n
         ENDIF
         IF(EVADIF.GT.EVAMAX) THEN
            EVAMAX=EVADIF
            NEVAMAX=n
         ENDIF
         IF(EVADIF.LT.EVAMIN) THEN
            EVAMIN=EVADIF
            NEVAMIN=n
         ENDIF
      endif                     ! charmv
      
      end do

      CLOSE(53)

!     Deallocate arrays
      DEALLOCATE (PHASEE,EMAG)
      
      if (charmv) then
!
      WRITE(16,7740)
 7740 FORMAT(///,5X,'THE LARGEST VALUES OF THE RATIO ', &
     &              'RESYNTHESIZED ELEV TIME SERIES/RAW TIME SERIES:',/)
      WRITE(16,7741) EAVMAX,NEAVMAX
      WRITE(16,7742) EVAMAX,NEVAMAX
      WRITE(16,7747)
 7747 FORMAT(/,5X,'THE LOWEST VALUES OF THE RATIO ', &
     &            'RESYNTHESIZED ELEV TIME SERIES/RAW TIME SERIES:',/)
      WRITE(16,7741) EAVMIN,NEAVMIN
      WRITE(16,7742) EVAMIN,NEVAMIN
 7741 FORMAT(9X,'  AVERAGE ELEVATION RATIO = ',E15.7,' AT NODE ',I8)
 7742 FORMAT(9X,' VARIANCE ELEVATION RATIO = ',E15.7,' AT NODE ',I8)
!     
      endif                     ! charmv
!
      return
      end subroutine

!***********************************************************************
!   Subroutine to solve the system and write output for velocity       *
!   globally.                                                          *
!                                                                      *
!   nf=0  if no steady constituent                                     *
!   nf=1  if steady constituent                                        *
!                                                                      *
!                        R.L. 11/10/95                                 *
!                                                                      *
!   added local LOGICAL CHARMV declaration 04/09/2004                  *
!   closed unit 54 08/23/05                                            *
!   Modified for SELFE - Andre Fortunato 2009/06/08                    *
!***********************************************************************
!
      SUBROUTINE LSQSOLVG(NP,XVELAV,YVELAV,XVELVA,YVELVA,CHARMV, &
     &                    myrank,TIMEBEG,DT,FMV,NTSTEPS,ITMV)
      IMPLICIT NONE
      LOGICAL CHARMV            !rl 04/09/2004      
      INTEGER NP,I,J,N,K,I1,I2,IT,IFR
      INTEGER NUAVMAX,NUAVMIN,NVAVMAX,NVAVMIN,NUVAMAX,NUVAMIN, &
     &  NVVAMAX,NVVAMIN,myrank
      REAL(rkind) UAV,VAV,USQ,VSQ,TIME,FTIME,RSU,RSV
      REAL(rkind) UAVMAX,VAVMAX,UVAMAX,VVAMAX,UAVMIN,VAVMIN, &
     & UVAMIN,VVAMIN,PHASEDU,PHASEDV,UAVDIF,VAVDIF,UMAGT,VMAGT, &
     & UVADIF,VVADIF
      REAL(8) CONVRD
!      REAL(rkind) XVELAV(NP),YVELAV(NP),XVELVA(NP),YVELVA(NP)
      REAL(rkind) XVELAV(NP),YVELAV(NP),XVELVA(NP),YVELVA(NP)
      REAL(rkind),ALLOCATABLE :: UMAG(:),VMAG(:),PHASEU(:),PHASEV(:)
      REAL(rkind),ALLOCATABLE :: Y(:)
      INTEGER NTSTEPS,ITMV
      REAL(8) TIMEBEG
      REAL(rkind) DT,FMV
!      COMMON /MEANSQ/ TIMEBEG,DT,FMV,NTSTEPS,ITMV
      CHARACTER*4 DIRNAME
!
      convrd=180.d0/pi
!
      ALLOCATE ( Y(2*MNHARF) )
      ALLOCATE ( UMAG(MNHARF),VMAG(MNHARF) )
      ALLOCATE ( PHASEU(MNHARF),PHASEV(MNHARF) )
!
!**** Open velocity station harmonic output file and write header information
!
      write(DIRNAME(1:4),'(i4.4)') myrank
      open(54,file=out_dir(1:len_out_dir)//'harmv.54'//DIRNAME(1:4))
      write(54,*) nfreq+nf
      do j=1,nfreq+nf
         write(54,3679) hafreq(j),HAFF(j),HAFACE(j),namefr(j)
      end do
 3679 format(1x,e20.10,1x,f10.7,1x,f12.8,1x,a10)
      write(54,*) NP
      
      if ( charmv ) then
         UAVMAX=-999.
         VAVMAX=-999.
         UVAMAX=-999.
         VVAMAX=-999.
         UAVMIN= 999.
         VAVMIN= 999.
         UVAMIN= 999.
         VVAMIN= 999.
      endif                     ! charmv
!
!***** AT each node transfer each load vector to p, solve and write output
!
      DO N=1,NP
         do k=1,mm
            hap(k) = GLOVLV(k,n)
         end do
         call fulsol(n)
         do k=1,mm
            y(k)=hax(k)
         end do
         do k=1,mm
            hap(k) = GLOULV(k,n)
         end do
         call fulsol(n)
         write(54,*) n
         do i=1,nfreq+nf
            if((nf.eq.1).and.(i.eq.1)) then
               umag(i)=hax(i)
               umagt=umag(i)/haff(i)
               vmag(i)=y(i)
               vmagt=vmag(i)/haff(i)
               phaseu(i)=0.
               phasev(i)=0.
            else
               i1=2*i-1-nf
               i2=i1+1
               umag(i)=sqrt(hax(i1)*hax(i1)+hax(i2)*hax(i2))
               umagt=umag(i)/haff(i)
               vmag(i)=sqrt(y(i1)*y(i1)+y(i2)*y(i2))
               vmagt=vmag(i)/haff(i)
               if((hax(i1).eq.0.).and.(hax(i2).eq.0.)) then
                  phaseu(i)=0.
               else
                  phaseu(i)=atan2(hax(i2),hax(i1))
               endif
               if((y(i1).eq.0.).and.(y(i2).eq.0.)) then
                  phasev(i)=0.
               else
                  phasev(i)=atan2(y(i2),y(i1))
               endif
            endif
            phasedu=convrd*phaseu(i)+haface(i)
            if(phasedu.lt.0.) phasedu=phasedu+360.d0
            if(phasedu.ge.360.d0) phasedu=phasedu-360.d0
            phasedv=convrd*phasev(i)+haface(i)
            if(phasedv.lt.0.) phasedv=phasedv+360.d0
            if(phasedv.ge.360.d0) phasedv=phasedv-360.d0

            write(54,6636) umagt,phasedu,vmagt,phasedv
 6636       format(2x,e16.8,1x,f11.4,2x,e16.8,1x,f11.4)
         end do

!HARMV...UNCOMMENT THE FOLLOWING LINES TO COMPUTE MEANS AND VARIANCES
!HARMV...FOR CHECKING THE HARMONIC ANALYSIS RESULTS.
!HARMV...Resynthesize the time series to compute the average and variances.
!HARMV...Compare resynthesized values with those computed during time stepping.
         if ( charmv ) then
            uav = 0.
            vav = 0.
            usq = 0.
            vsq = 0.
            do it=1,ntsteps
               TIME=TIMEBEG+DT*IT
               rsu=0.
               rsv=0.
               do ifr=1,nfreq+nf
                  ftime=hafreq(ifr)*time
                  rsu=rsu+umag(ifr)*cos(ftime-phaseu(ifr))
                  rsv=rsv+vmag(ifr)*cos(ftime-phasev(ifr))
               end do
               uav=uav+rsu
               vav=vav+rsv
               usq=usq+rsu*rsu
               vsq=vsq+rsv*rsv
            end do

            uav=uav/ntsteps
            vav=vav/ntsteps
            usq=usq/ntsteps-uav*uav
            vsq=vsq/ntsteps-vav*vav
            if(xvelav(n).eq.0.) then
               if(uav.eq.0.) uavdif=1.0d0
               if(uav.ne.0.) uavdif=99e19
            else
               uavdif=uav/xvelav(n)
            endif
            if(yvelav(n).eq.0.) then
               if(vav.eq.0.) vavdif=1.0d0
               if(vav.ne.0.) vavdif=99e19
            else
               vavdif=vav/yvelav(n)
            endif
            if(xvelva(n).eq.0.) then
               if(usq.eq.0.) uvadif=1.0d0
               if(usq.ne.0.) uvadif=99e19
            else
               uvadif=usq/xvelva(n)
            endif
            if(yvelva(n).eq.0.) then
               if(vsq.eq.0.) vvadif=1.0d0
               if(vsq.ne.0.) vvadif=99e19
            else
               vvadif=vsq/yvelva(n)
            endif
            write(55,*) n
            write(55,7637) xvelav(n),uav,uavdif,xvelva(n),usq,uvadif
            write(55,7637) yvelav(n),vav,vavdif,yvelva(n),vsq,vvadif
 7637       format(2x,3(e16.8,1x),2x,3(e16.8,1x))

            IF(UAVDIF.GT.UAVMAX) THEN
               UAVMAX=UAVDIF
               NUAVMAX=n
            ENDIF
            IF(UAVDIF.LT.UAVMIN) THEN
               UAVMIN=UAVDIF
               NUAVMIN=n
            ENDIF
            IF(VAVDIF.GT.VAVMAX) THEN
               VAVMAX=VAVDIF
               NVAVMAX=n
            ENDIF
            IF(VAVDIF.LT.VAVMIN) THEN
               VAVMIN=VAVDIF
               NVAVMIN=n
            ENDIF
            IF(UVADIF.GT.UVAMAX) THEN
               UVAMAX=UVADIF
               NUVAMAX=n
            ENDIF
            IF(UVADIF.LT.UVAMIN) THEN
               UVAMIN=UVADIF
               NUVAMIN=n
            ENDIF
            IF(VVADIF.GT.VVAMAX) THEN
               VVAMAX=VVADIF
               NVVAMAX=n
            ENDIF
            IF(VVADIF.LT.VVAMIN) THEN
               VVAMIN=VVADIF
               NVVAMIN=n
            ENDIF

         endif                  !  charmv

      end do

      CLOSE(54)

!  Deallocate arrays
      DEALLOCATE ( Y )
      DEALLOCATE ( UMAG,VMAG )
      DEALLOCATE ( PHASEU,PHASEV )

      if ( charmv ) then 
!
         WRITE(16,7740)
 7740    FORMAT(///,5X,'THE LARGEST VALUES OF THE RATIO ', &
     &              'RESYNTHESIZED VEL TIME SERIES/RAW TIME SERIES:',/)
         WRITE(16,7743) UAVMAX,NUAVMAX
         WRITE(16,7744) UVAMAX,NUVAMAX
         WRITE(16,7745) VAVMAX,NVAVMAX
         WRITE(16,7746) VVAMAX,NVVAMAX
         WRITE(16,7747)
 7747    FORMAT(//,5X,'THE LOWEST VALUES OF THE RATIO ', &
     &             'RESYNTHESIZED VEL TIME SERIES/RAW TIME SERIES:',/)
         WRITE(16,7743) UAVMIN,NUAVMIN
         WRITE(16,7744) UVAMIN,NUVAMIN
         WRITE(16,7745) VAVMIN,NVAVMIN
         WRITE(16,7746) VVAMIN,NVVAMIN
 7743    FORMAT(9X,' AVERAGE U VELOCITY RATIO = ',E15.7,' AT NODE ',I8)
 7744    FORMAT(9X,'VARIANCE U VELOCITY RATIO = ',E15.7,' AT NODE ',I8)
 7745    FORMAT(9X,' AVERAGE V VELOCITY RATIO = ',E15.7,' AT NODE ',I8)
 7746    FORMAT(9X,'VARIANCE V VELOCITY RATIO = ',E15.7,' AT NODE ',I8)
!     
      endif                     ! charmv
!
      return
      end subroutine

!
!***********************************************************************
!   Subroutine to initialize parameters for harmonic analysis with a   *
!   cold start.                                                        *
!                                                                      *
!                        R.L.  11/9/95                                 *
!***********************************************************************
!
      SUBROUTINE HACOLDS(HAFREQ)
      implicit none
      INTEGER I,J
      REAL(rkind) HAFREQ(MNHARF)
!
      if (hafreq(1).eq.0.0) then
         nz=0
         nf=1
      else
         nz=1
         nf=0
      endif
!
      nfreq=nfreq-nf
      mm=2*nfreq+nf
!
      do i=1,mm
         do j=1,mm
            ha(i,j)=0.
         end do
      end do
      icall=0
!     
      return
      end subroutine

!***********************************************************************
!   Subroutine to initialize global elevation load vectors for         *
!   harmonic analysis with a cold start.                               *
!                                                                      *
!                        R.L.  11/9/95                                 *
!***********************************************************************
!
      SUBROUTINE HACOLDSEG(NP)
      implicit none
      INTEGER NP,I,N
!
      do i=1,mm
         do N=1,NP
            GLOELV(I,N)=0.
         end do
      end do
!
      return
      end subroutine

!***********************************************************************
!   Subroutine to initialize global velocity load vectors for          *
!   harmonic analysis with a cold start.                               *
!                                                                      *
!                        R.L.  11/9/95                                 *
!***********************************************************************
!
      SUBROUTINE HACOLDSVG(NP)
      implicit none
      INTEGER NP,I,N
!
      do i=1,mm
         do N=1,NP
            GLOULV(I,N)=0.
            GLOVLV(I,N)=0.
         end do
      end do
!
      return
      end subroutine

!***********************************************************************
!   Subroutine to read in and initialize harmonic analysis for a hot   *
!   start.                                                             *
!                                                                      *
!   Checks are made to ensure agreement between values read in from    *
!   the hotstart file and values read in from the UNIT 15 file.        *
!                                                                      *
!                        R.L. 11/9/95                                  *
!   Adapted for SELFE: Andre Fortunato, 2009/06/08
!***********************************************************************
!
      SUBROUTINE HAHOTS(NSTAE,NSTAV,NP,ISTAE,ISTAV,IGLOE,IGLOV, &
     &   IHOT,myrank)
      IMPLICIT NONE
      INTEGER NSTAE,NSTAV,NP, myrank
      INTEGER ISTAE, ISTAV, IGLOE, IGLOV, IFLAG, I, J
!     
      INTEGER INFREQ, INSTAE, INSTAV, INP, INZ, INF
      INTEGER IISTAE, IISTAV, IIGLOE, IIGLOV,NSCREEN 
      INTEGER IHOT, IMM, IICALL
      REAL(rkind) FDIFF
!     
      REAL(rkind),ALLOCATABLE ::  IFREQ(:),IFF(:),IFACE(:)
      CHARACTER*10,ALLOCATABLE :: INAMEFR(:)
!     
      CHARACTER*16 FNAME
      CHARACTER*8 FNAM8(2)
      EQUIVALENCE (FNAM8(1),FNAME)
!     
!***** Compute parameter values for checking
!
      if (hafreq(1).eq.0.0) then
         nz=0
         nf=1
      else
         nz=1
         nf=0
      endif
      nfreq=nfreq-nf
      mm=2*nfreq+nf

      ALLOCATE ( IFREQ(MNHARF),IFF(MNHARF),IFACE(MNHARF) )
      ALLOCATE ( INAMEFR(MNHARF) )
!
!***** Read in and check various parameter values
!
      READ(IHOT) inz
      READ(IHOT) inf
      READ(IHOT) imm
      READ(IHOT) inp
      READ(IHOT) instae
      READ(IHOT) instav
      READ(IHOT) iistae
      READ(IHOT) iistav
      READ(IHOT) iigloe
      READ(IHOT) iiglov
      READ(IHOT) iicall
      READ(IHOT) infreq
!
      iflag=0
      if(nz.ne.inz) iflag=1
      if(nf.ne.inf) iflag=1
      if(mm.ne.imm) iflag=1
      if(np.ne.inp) iflag=1
      if(nstae.ne.instae) iflag=1
      if(nstav.ne.instav) iflag=1
      if(istae.ne.iistae) iflag=1
      if(istav.ne.iistav) iflag=1
      if(igloe.ne.iigloe) iflag=1
      if(iglov.ne.iiglov) iflag=1
      if(nfreq.ne.infreq) iflag=1
!
      do i=1,nfreq+nf
         READ(IHOT) FNAM8(1)
         READ(IHOT) FNAM8(2)
         INAMEFR(I) = FNAME
         read(IHOT) ifreq(i)
         read(IHOT) iff(i)
         read(IHOT) iface(i)

         if(namefr(i).ne.inamefr(i)) iflag=1
         if(abs(hafreq(i)+ifreq(i)).lt.1.0d-30) then
            fdiff=0.
         else
            fdiff=abs(hafreq(i)-ifreq(i))/abs(hafreq(i)+ifreq(i))
         endif
         if(fdiff.ge.1.d-6) iflag=1
         if(abs(HAFF(i)+iFF(i)).lt.1d-30) then
            fdiff=0.
         else
            fdiff=abs(HAFF(i)-iFF(i))/abs(HAFF(i)+iFF(i))
         endif
         if(fdiff.ge.1.d-6) iflag=1
         if(abs(HAFACE(i)+iFACE(i)).lt.1d-30) then
            fdiff=0.
         else
            fdiff=abs(HAFACE(i)-iFACE(i))/abs(HAFACE(i)+iFACE(i))
         endif
         if(fdiff.ge.1.d-6) iflag=1
      end do
      if(iflag.eq.1) goto 999
!
!***** Read in time of most recent H.A. update
!
      READ(IHOT) TIMEUD
      READ(IHOT) ITUD
!
!***** Read in RHS Matrix
!
      do i=1,mm
         do j=1,mm
            READ(IHOT) HA(I,J)
         end do
      end do

!
!***** FATAL Error Messages
!
 999  continue
      if(iflag.ne.0) then
         if(NSCREEN.EQ.1.AND.myrank.EQ.0) write(6,1000)
!        Use 12 for non-fatal messages; need parallel_abort() if using 11
         write(12,1000)
 1000    FORMAT(////,5x,'***** DISCREPANCY IN HARMONIC ANALYSIS HOT ', &
     &        'START FILE *****',/)
      endif

      if(iflag.eq.1) then
         if(nz.ne.inz) then
            if(NSCREEN.EQ.1.AND.myrank.EQ.0) write(6,2010) inz,nz
            write(12,2010) inz,nz
 2010       format(5x,'NZ COMPUTED FROM UNIT 14 INPUT = ',I2, &
     &           ', NZ READ IN FROM HOT START FILE = ',I2,/)
         endif
         if(nf.ne.inf) then
            if(NSCREEN.EQ.1.AND.myrank.EQ.0) write(6,2020) inf,nf
            write(12,2020) inf,nf
 2020       format(5x,'NF COMPUTED FROM UNIT 14 INPUT = ',I2, &
     &           ', NF READ IN FROM HOT START FILE = ',I2,/)
         endif
         if(mm.ne.imm) then
            if(NSCREEN.EQ.1.AND.myrank.EQ.0) write(6,2030) imm,mm
            write(12,2030) imm,mm
 2030       format(5x,'MM COMPUTED FROM UNIT 14 INPUT = ',I2, &
     &           ', MM READ IN FROM HOT START FILE = ',I2,/)
         endif
         if(np.ne.inp) then
            if(NSCREEN.EQ.1.AND.myrank.EQ.0) write(6,2040) inp,np
            write(12,2040) inp,np
 2040       format(5x,'NP READ IN FROM UNIT 15 = ',I2, &
     &           ', NP READ IN FROM HOT START FILE = ',I2,/)
         endif
         if(nstae.ne.instae) then
            if(NSCREEN.EQ.1.AND.myrank.EQ.0) write(6,2050) instae,nstae
            write(12,2050) instae,nstae
 2050       format(5x,'NSTAE READ IN FROM UNIT 15 = ',I2, &
     &           ', NSTAE READ IN FROM HOT START FILE = ',I2,/)
         endif
         if(nstav.ne.instav) then
            if(NSCREEN.EQ.1.AND.myrank.EQ.0) write(6,2060) instav,nstav
            write(12,2060) instav,nstav
 2060       format(5x,'NSTAV READ IN FROM UNIT 15 = ',I2, &
     &           ', NSTAV READ IN FROM HOT START FILE = ',I2,/)
         endif
         if(istae.ne.iistae) then
            if(NSCREEN.EQ.1.AND.myrank.EQ.0) write(6,2070) iistae,istae
            write(12,2070) iistae,istae
 2070       format(5x,'ISTAE READ IN FROM UNIT 15 = ',I2, &
     &           ', ISTAE READ IN FROM HOT START FILE = ',I2,/)
         endif
         if(istav.ne.iistav) then
            if(NSCREEN.EQ.1.AND.myrank.EQ.0) write(6,2080) iistav,istav
            write(12,2080) iistav,istav
 2080       format(5x,'ISTAV READ IN FROM UNIT 15 = ',I2, &
     &           ', ISTAV READ IN FROM HOT START FILE = ',I2,/)
         endif
         if(igloe.ne.iigloe) then
            if(NSCREEN.EQ.1.AND.myrank.EQ.0) write(6,2090) iigloe,igloe
            write(12,2090) iigloe,igloe
 2090       format(5x,'IGLOE READ IN FROM UNIT 15 = ',I2, &
     &           ', IGLOE READ IN FROM HOT START FILE = ',I2,/)
         endif
         if(iglov.ne.iiglov) then
            if(NSCREEN.EQ.1.AND.myrank.EQ.0) write(6,2100) iiglov,iglov
            write(12,2100) iiglov,iglov
 2100       format(5x,'IGLOV READ IN FROM UNIT 15 = ',I2, &
     &           ', IGLOV READ IN FROM HOT START FILE = ',I2,/)
         endif
         if(nfreq.ne.infreq) then
            if(NSCREEN.EQ.1.AND.myrank.EQ.0) write(6,2110) infreq,nfreq
            write(12,2110) infreq,nfreq
 2110       format(5x,'NFREQ COMPUTED FROM UNIT 15 INPUT = ',I2, &
     &           ', NFREQ READ IN FROM HOT START FILE = ',I2,/)
         endif
         do i=1,nfreq+nf
            if(namefr(i).ne.inamefr(i)) then
               if(NSCREEN.EQ.1.AND.myrank.EQ.0) &
     &              write(6,2120) i,inamefr(i),namefr(i)
               write(12,2120) i,namefr(i),namefr(i)
 2120          format(5x,'FOR CONSTITUENT # ',I3, &
     &              ', NAMEFR READ IN FROM UNIT 15 = ',A10, &
     &              ', NAMEFR READ IN FROM HOT START FILE = ',A10,/)
            endif
            if(hafreq(i).ne.ifreq(i)) then
               if(NSCREEN.EQ.1.AND.myrank.EQ.0)  &
     &              write(6,2130) i,ifreq(i),hafreq(i)
               write(12,2130) i,ifreq(i),hafreq(i)
 2130          format(5x,'FOR CONSTITUENT # ',I3, &
     &              ', FREQ READ IN FROM UNIT 15 = ',D20.10, &
     &              ', FREQ READ IN FROM HOT START FILE = ',D20.10,/)
            endif
            if(HAFF(i).ne.iFF(i)) then
               if(NSCREEN.EQ.1.AND.myrank.EQ.0) &
     &              write(6,2140) i,iff(i),haff(i)
               write(12,2140) i,iff(i),haff(i)
 2140          format(5x,'FOR CONSTITUENT # ',I3, &
     &              ', FF READ IN FROM UNIT 15 = ',F10.5, &
     &              ', FF READ IN FROM HOT START FILE = ',F10.5,/)
            endif
            if(HAFACE(i).ne.iFACE(i)) then
               if(NSCREEN.EQ.1.AND.myrank.EQ.0) &
     &              write(6,2150) i,iface(i),haface(i)
               write(12,2150) i,iface(i),haface(i)
 2150          format(5x,'FOR CONSTITUENT # ',I3, &
     &              ', FACE READ IN FROM UNIT 15 = ',F10.5, &
     &              ', FACE READ IN FROM HOT START FILE = ',F10.5,/)
            endif
         end do

         DEALLOCATE ( IFREQ,IFF,IFACE )
         DEALLOCATE ( INAMEFR )

         if(NSCREEN.EQ.1.AND.myrank.EQ.0) write(6,1010)
         write(12,1010)
 1010    FORMAT(//,5x,'********** RUN TERMINATED **********',/)
         stop
      endif

      return
      end subroutine


!***********************************************************************
!   Subroutine to read in and initialize the global elevation load     *
!   vector for harmonic analysis with a hot start.                     *
!                                                                      *
!                        R.L. 11/9/95                                  *
!   Adapted for SELFE: Andre Fortunato, 2009/06/08
!***********************************************************************
!
      SUBROUTINE HAHOTSEG(NP_GLOBAL,IHOT,myrank)
      implicit none
      integer np, ihot, n,i,np_global,myrank,nl
      real(rkind)  tmp
!
!***** Read in Global Elevation LHS load vector
!
      do n=1,np_global
         do i=1,mm
!            READ(IHOT) GLOELV(I,N)
            READ(IHOT) tmp
            if (ipgl(i)%rank==myrank) then
               nl = ipgl(n)%id
               GLOELV(i,nl)=tmp
            endif
         end do
      end do
!
      return
      end subroutine

!***********************************************************************
!   Subroutine to read in and initialize the global velocity load      *
!   vector for harmonic analysis with a hot start.                     *
!                                                                      *
!                        R.L. 11/9/95                                  *
!   Adapted for SELFE: Andre Fortunato, 2009/06/08
!***********************************************************************
!
      SUBROUTINE HAHOTSVG(NP_GLOBAL,IHOT,myrank)
      implicit none
      integer np, ihot, n, i,np_global,myrank,nl
      real(rkind)  tmp1,tmp2
!
!***** Read in Global Velocity LHS load vector
!
      do n=1,np_global
         do i=1,mm
!            READ(IHOT) GLOULV(I,N)
!            READ(IHOT) GLOVLV(I,N)
            READ(IHOT) tmp1
            READ(IHOT) tmp2
            if (ipgl(i)%rank==myrank) then
               nl = ipgl(n)%id
               GLOULV(i,nl)=tmp1
               GLOVLV(i,nl)=tmp2
            endif
         end do
      end do
!
      return
      end subroutine

!***********************************************************************
!   Subroutine to write out to the hotstart file (UNITS 67 and 68)     *
!   header information and the LHS matrix for the harmonic analysis    *
!                                                                      *
!                        R.L.  11/8/95                                 *
!     jgf45.07 name changes to variables ISTAE:IGLOV -> NHASE:NHAGV    *
!***********************************************************************
!
      SUBROUTINE HAHOUT(NP,NSTAE,NSTAV,NHASE,NHASV,NHAGE,NHAGV, &
     &  IOUNIT,IHOTSTP)
      implicit none
      INTEGER NP,NSTAE,NSTAV,NHASE,AE,NHASV
      INTEGER NHAGE,NHAGV,IOUNIT,IHOTSTP,I,J
      CHARACTER*16 FNAME
      CHARACTER*8 FNAM8(2)
      EQUIVALENCE (FNAM8(1),FNAME)

!
!***** Write Out various parameter values
!
      WRITE(IOUNIT,REC=IHOTSTP+1) NZ
      WRITE(IOUNIT,REC=IHOTSTP+2) NF
      WRITE(IOUNIT,REC=IHOTSTP+3) MM
      WRITE(IOUNIT,REC=IHOTSTP+4) NP
      WRITE(IOUNIT,REC=IHOTSTP+5) NSTAE
      WRITE(IOUNIT,REC=IHOTSTP+6) NSTAV
      WRITE(IOUNIT,REC=IHOTSTP+7) NHASE
      WRITE(IOUNIT,REC=IHOTSTP+8) NHASV
      WRITE(IOUNIT,REC=IHOTSTP+9) NHAGE
      WRITE(IOUNIT,REC=IHOTSTP+10) NHAGV
      WRITE(IOUNIT,REC=IHOTSTP+11) ICALL
      WRITE(IOUNIT,REC=IHOTSTP+12) NFREQ
      IHOTSTP = IHOTSTP+12

      do i=1,nfreq+nf
         FNAME=NAMEFR(I)
         WRITE(IOUNIT,REC=IHOTSTP+1) FNAM8(1)
         WRITE(IOUNIT,REC=IHOTSTP+2) FNAM8(2)
         IHOTSTP=IHOTSTP+2
         WRITE(IOUNIT,REC=IHOTSTP+1) hafreq(i)
         WRITE(IOUNIT,REC=IHOTSTP+2) HAFF(i)
         WRITE(IOUNIT,REC=IHOTSTP+3) HAFACE(i)
         IHOTSTP=IHOTSTP+3
      end do

!
!***** Write Out time of most recent H.A. update
!
      WRITE(IOUNIT,REC=IHOTSTP+1) TIMEUD
      WRITE(IOUNIT,REC=IHOTSTP+2) ITUD
      IHOTSTP=IHOTSTP+2
!
!***** Write Out LHS Matrix
!
      do i=1,mm
         do j=1,mm
            IHOTSTP = IHOTSTP + 1
            WRITE(IOUNIT,REC=IHOTSTP) HA(I,J)
         END DO
      END DO

      return
      end subroutine

!***********************************************************************
!   Subroutine to write global elevation harmonic analysis RHS load    *
!   vector to a hot start file (UNITS 67 and 68)                       *
!                                                                      *
!                        R.L.  11/8/95                                 *
!***********************************************************************
!
      SUBROUTINE HAHOUTEG(NP,IOUNIT,IHOTSTP)
      implicit none
      INTEGER IOUNIT
      INTEGER NP,IHOTSTP,N,I 
!
!***** Write Out Global Elevation RHS load vector
!
      do n=1,np
         do i=1,mm
            IHOTSTP=IHOTSTP+1
            WRITE(IOUNIT,REC=IHOTSTP) GLOELV(I,N)
         end do
      end do
      
      return
      end subroutine

!***********************************************************************
!   Subroutine to write global velocity harmonic analysis RHS load     *
!   vector to a hot start file (UNITS 67 and 68)                       *
!                                                                      *
!                        R.L.  11/8/95                                 *
!***********************************************************************
!
      SUBROUTINE HAHOUTVG(NP,IOUNIT,IHOTSTP)
      implicit none
      INTEGER NP,IOUNIT,IHOTSTP,N,I
!
!***** Write Out Global Velocity RHS load vector
!
      do n=1,np
         do i=1,mm
            IHOTSTP=IHOTSTP+1
            WRITE(IOUNIT,REC=IHOTSTP) GLOULV(I,N)
            IHOTSTP=IHOTSTP+1
            WRITE(IOUNIT,REC=IHOTSTP) GLOVLV(I,N)
         end do
      end do
      
      return
      end subroutine

      END MODULE HARM
