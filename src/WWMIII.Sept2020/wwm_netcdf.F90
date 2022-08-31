#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CF_EXTRACT_TIME(eStrUnitTime, ConvertToDay, eTimeStart)
      USE DATAPOOL
      IMPLICIT NONE
      character(len=100), intent(in) :: eStrUnitTime
      real(rkind), intent(out) :: ConvertToDay, eTimeStart
      character (len=100) :: Xname, Yname
      character (len=10) :: YnameYear, YnameMonth, YnameDay
      character (len=10) :: YnameHour, YnameMin, YnameSec
      character (len=50) :: YnameB, YnameD, YnameE
      character (len=50) :: YnameDate, YnameTime, YnameTimeP
      character (len=15) :: eStrTime
      integer alenB, alenC, alenD, alenE, alenTime, alenDate
      integer alen, posBlank
      integer lenHour, lenMin, lenSec, lenMonth, lenDay, posSepDateTime
      alen=LEN_TRIM(eStrUnitTime)
      posBlank=INDEX(eStrUnitTime(1:alen), ' ')
      Xname=eStrUnitTime(1:posBlank-1) ! should be days/hours/seconds
      IF (TRIM(Xname) .eq. 'days') THEN
        ConvertToDay=1
      ELSEIF (TRIM(Xname) .eq. 'hours') THEN
        ConvertToDay=1/24
      ELSEIF (TRIM(Xname) .eq. 'seconds') THEN
        ConvertToDay=1/86400
      ELSE
        CALL WWM_ABORT('Error in the code for conversion')
      END IF
      !
      Yname=eStrUnitTime(posBlank+1:alen)
      alenB=LEN_TRIM(Yname)
      posBlank=INDEX(Yname(1:alenB), ' ')
      YnameB=Yname(posBlank+1:alenB) ! should be 1990-01-01 0:0:0
      !
      alenC=LEN_TRIM(YnameB)
      posSepDateTime=INDEX(YnameB(1:alenC), ' ')
      IF (posSepDateTime .gt. 0) THEN
        YnameDate=YnameB(1:posSepDateTime-1) ! should be 1990-01-01
        YnameTimeP=YnameB(posSepDateTime+1:alenC) ! should be 0:0:0
        alenC=LEN_TRIM(YnameTimeP)
        posBlank=INDEX(YnameTimeP(1:alenC), ' ')
        IF (posBlank .eq. 0) THEN
          YnameTime=YnameTimeP
        ELSE
          YnameTime=YnameTimeP(1:posBlank-1)
        END IF
      ELSE
        YnameDate=YnameB
        eStrTime(10:10)='0'
        eStrTime(11:11)='0'
        eStrTime(12:12)='0'
        eStrTime(13:13)='0'
        eStrTime(14:14)='0'
        eStrTime(15:15)='0'
      END IF
      !
      alenDate=LEN_TRIM(YnameDate)
      posBlank=INDEX(YnameDate(1:alenDate), '-')
      YnameYear=YnameDate(1:posBlank-1) ! should be 1990
      YnameD=YnameDate(posBlank+1:alenDate)
      alenD=LEN_TRIM(YnameD)
      posBlank=INDEX(YnameD(1:alenD), '-')
      YnameMonth=YnameD(1:posBlank-1) ! should be 01
      YnameDay=YnameD(posBlank+1:alenD) ! should be 01
      !
      ! year
      eStrTime( 1: 1)=YnameYear( 1: 1)
      eStrTime( 2: 2)=YnameYear( 2: 2)
      eStrTime( 3: 3)=YnameYear( 3: 3)
      eStrTime( 4: 4)=YnameYear( 4: 4)
      !
      ! month
      WRITE(WINDBG%FHNDL,*) 'YnameMonth=', YnameMonth
      lenMonth=LEN_TRIM(YnameMonth)
      IF (lenMonth .eq. 2) THEN
        eStrTime( 5: 5)=YnameMonth( 1: 1)
        eStrTime( 6: 6)=YnameMonth( 2: 2)
      ELSE
        IF (lenMonth .eq. 1) THEN
          eStrTime( 5: 5)='0'
          eStrTime( 6: 6)=YnameMonth( 1: 1)
        ELSE
          CALL WWM_ABORT('DIE in trying to get the month')
        END IF
      END IF
      !
      ! day
      lenDay=LEN_TRIM(YnameDay)
      IF (lenDay .eq. 2) THEN
        eStrTime( 7: 7)=YnameDay( 1: 1)
        eStrTime( 8: 8)=YnameDay( 2: 2)
      ELSE
        IF (lenDay .eq. 1) THEN
          eStrTime( 7: 7)='0'
          eStrTime( 8: 8)=YnameDay( 1: 1)
        ELSE
          CALL WWM_ABORT('DIE in trying to get the day')
        END IF
      END IF
      !
      eStrTime( 9: 9)='.'
      !
      IF (posSepDateTime .gt. 0) THEN
        !
        alenTime=LEN_TRIM(YnameTime)
        posBlank=INDEX(YnameTime(1:alenTime), ':')
        YnameHour=YnameTime(1:posBlank-1) ! should be 0
        YnameE=YnameTime(posBlank+1:alenTime)
        alenE=LEN_TRIM(YnameE)
        posBlank=INDEX(YnameE(1:alenE), ':')
        YnameMin=YnameE(1:posBlank-1) ! should be 0
        YnameSec=YnameE(posBlank+1:alenE) ! should be 0
        !
        !
        ! Hour
        lenHour=LEN_TRIM(YnameHour)
        IF (lenHour .eq. 2) THEN
          eStrTime(10:10)=YnameHour( 1: 1)
          eStrTime(11:11)=YnameHour( 2: 2)
        ELSE
          IF (lenHour .eq. 1) THEN
            eStrTime(10:10)='0'
            eStrTime(11:11)=YnameHour( 1: 1)
          ELSE
            CALL WWM_ABORT('DIE in trying to get the hour')
          END IF
        END IF
        !
        ! Min
        lenMin=LEN_TRIM(YnameMin)
        IF (lenMin .eq. 2) THEN
          eStrTime(12:12)=YnameMin( 1: 1)
          eStrTime(13:13)=YnameMin( 2: 2)
        ELSE
          IF (lenMin .eq. 1) THEN
            eStrTime(12:12)='0'
            eStrTime(13:13)=YnameMin( 1: 1)
          ELSE
            CALL WWM_ABORT('DIE in trying to get the min')
          END IF
        END IF
        !
        ! Sec
        lenSec=LEN_TRIM(YnameSec)
        IF (lenSec .eq. 2) THEN
          eStrTime(14:14)=YnameSec( 1: 1)
          eStrTime(15:15)=YnameSec( 2: 2)
        ELSE
          IF (lenSec .eq. 1) THEN
            eStrTime(14:14)='0'
            eStrTime(15:15)=YnameSec( 1: 1)
          ELSE
            WRITE(WINDBG%FHNDL,*) 'YnameSec=', TRIM(Ynamesec)
            WRITE(WINDBG%FHNDL,*) 'lenSec=', lenSec
            FLUSH(WINDBG%FHNDL)
            CALL WWM_ABORT('DIE in trying to get the sec')
          END IF
        END IF
      END IF
      WRITE(WINDBG%FHNDL,*) 'eStrTime=', eStrTime
      CALL CT2MJD(eStrTime, eTimeStart)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DETERMINE_NEEDED_COMPUTATION(eVar)
      USE DATAPOOL, only : VAROUT, OUTVARS_COMPLETE
      implicit none
      type(VAROUT), intent(inout) :: eVar
      LOGICAL     ::   HS, TM01, TM02, TM10, KLM, WLM,                  &
     &   ETOTC, ETOTS, DM, DSPR,                                        &
     &   TPPD, CPPD, KPPD, CGPD,                                        &
     &   TPP, CPP, WNPP, CGPP, KPP, LPP, PEAKD, PEAKDSPR,               &
     &   DPEAK, UBOT, ORBITAL, BOTEXPER, TMBOT,                         &
     &   URSELL, UFRIC, Z0, ALPHA_CH, WINDX, WINDY, CD,                 &
     &   CURRTX, CURRTY, WATLEV, WATLEVOLD, DEPDT, DEP,                 &
     &   WINDMAG, TAUW, TAUWX, TAUWY, TAUHF, TAUTOT,                    &
     &   STOKESBOTTX, STOKESBOTTY,                                      &
     &   STOKESSURFX, STOKESSURFY, STOKESBAROX, STOKESBAROY,            &
     &   RSXX, RSXY, RSYY, CFL1, CFL2, CFL3, ZETA_SETUP
      LOGICAL :: ComputeMean, ComputeDirSpread, ComputePeak
      LOGICAL :: ComputeCurr, ComputeUrsell, ComputeStokes
      integer iVar, idx, nbOutVarEff
      integer istat
      HS           = eVar%LVAR( 1)
      TM01         = eVar%LVAR( 2)
      TM02         = eVar%LVAR( 3)
      TM10         = eVar%LVAR( 4)
      KLM          = eVar%LVAR( 5)
      WLM          = eVar%LVAR( 6)
      ETOTC        = eVar%LVAR( 7)
      ETOTS        = eVar%LVAR( 8)
      DM           = eVar%LVAR( 9)
      DSPR         = eVar%LVAR(10)
      TPPD         = eVar%LVAR(11)
      CPPD         = eVar%LVAR(12)
      KPPD         = eVar%LVAR(13)
      CGPD         = eVar%LVAR(14)
      TPP          = eVar%LVAR(15)
      CPP          = eVar%LVAR(16)
      WNPP         = eVar%LVAR(17)
      CGPP         = eVar%LVAR(18)
      KPP          = eVar%LVAR(19)
      LPP          = eVar%LVAR(20)
      PEAKD        = eVar%LVAR(21)
      PEAKDSPR     = eVar%LVAR(22)
      DPEAK        = eVar%LVAR(23)
      UBOT         = eVar%LVAR(24)
      ORBITAL      = eVar%LVAR(25)
      BOTEXPER     = eVar%LVAR(26)
      TMBOT        = eVar%LVAR(27)
      URSELL       = eVar%LVAR(28)
      UFRIC        = eVar%LVAR(29)
      Z0           = eVar%LVAR(30)
      ALPHA_CH     = eVar%LVAR(31)
      WINDX        = eVar%LVAR(32)
      WINDY        = eVar%LVAR(33)
      CD           = eVar%LVAR(34)
      CURRTX       = eVar%LVAR(35)
      CURRTY       = eVar%LVAR(36)
      WATLEV       = eVar%LVAR(37)
      WATLEVOLD    = eVar%LVAR(38)
      DEPDT        = eVar%LVAR(39)
      DEP          = eVar%LVAR(40)
      WINDMAG      = eVar%LVAR(41)
      TAUW         = eVar%LVAR(42)
      TAUWX        = eVar%LVAR(43)
      TAUWY        = eVar%LVAR(44)
      TAUHF        = eVar%LVAR(45)
      TAUTOT       = eVar%LVAR(46)
      STOKESBOTTX  = eVar%LVAR(47)
      STOKESBOTTY  = eVar%LVAR(48)
      STOKESSURFX  = eVar%LVAR(49)
      STOKESSURFY  = eVar%LVAR(50)
      STOKESBAROX  = eVar%LVAR(51)
      STOKESBAROY  = eVar%LVAR(52)
      RSXX         = eVar%LVAR(53)
      RSXY         = eVar%LVAR(54)
      RSYY         = eVar%LVAR(55)
      CFL1         = eVar%LVAR(56)
      CFL2         = eVar%LVAR(57)
      CFL3         = eVar%LVAR(58)
      ZETA_SETUP   = eVar%LVAR(59)
      ComputeMean=.FALSE.
      ComputeDirSpread=.FALSE.
      ComputePeak=.FALSE.
      ComputeCurr=.FALSE.
      ComputeUrsell=.FALSE.
      ComputeStokes=.FALSE.
      IF (HS .or. TM01 .or. TM02 .or. TM10 .or. KLM .or. WLM) THEN
        ComputeMean=.TRUE.
      END IF
      IF (ETOTC .or. ETOTS .or. DM .or. DSPR) THEN
        ComputeDirspread=.TRUE.
      END IF
      IF (TPPD .or. CPPD .or. KPPD .or. CGPD .or. TPP .or. CPP .or. WNPP .or. CGPP .or. KPP .or. LPP .or. PEAKD .or. PEAKDSPR .or. DPEAK) THEN
        ComputePeak=.TRUE.
      END IF
      IF (UBOT .or. ORBITAL .or. BOTEXPER .or. TMBOT) THEN
        ComputeCurr=.TRUE.
      END IF
      IF (URSELL) THEN
        ComputeUrsell=.TRUE.
        ComputeMean=.TRUE.
        ComputePeak=.TRUE.
      END IF
      IF (STOKESSURFX .or. STOKESSURFY .or. STOKESBOTTX .or. STOKESBOTTY .or. STOKESBAROX .or. STOKESBAROY) THEN
        ComputeStokes=.TRUE.
      END IF
      eVar%ComputeMean=ComputeMean
      eVar%ComputeDirspread=ComputeDirspread
      eVar%ComputePeak=ComputePeak
      eVar%ComputeCurr=ComputeCurr
      eVar%ComputeUrsell=ComputeUrsell
      eVar%ComputeStokes=ComputeStokes
      nbOutVarEff=0
      DO iVar=1,OUTVARS_COMPLETE
        IF (eVar%LVAR(iVar)) THEN
          nbOutVarEff=nbOutVarEff+1
        END IF
      END DO
      eVar%nbOutVarEff=nbOutVarEff
      allocate(eVar%ListIdxEff(nbOutVarEff), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_netcdf, allocate error 1')
      idx=0
      DO iVar=1,OUTVARS_COMPLETE
        IF (eVar%LVAR(iVar)) THEN
          idx=idx+1
          eVar%ListIdxEff(idx)=iVar
        END IF
      END DO
      END SUBROUTINE
#ifdef NCDF
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SERIAL_GET_BOUNDARY_NEXTGENERATION(TheBound)
      USE DATAPOOL
      IMPLICIT NONE
      type(BoundaryInfo), intent(inout) :: TheBound
      integer :: ContElements(np_total)
      integer :: ListDegWork(np_total)
      integer, allocatable :: ListAdjWithDupl(:,:)
      integer, allocatable :: IEcontain(:,:)
      integer, allocatable :: StatusAdj(:)
      integer IE, I, IP, INEXT, IPREV, IPadj
      integer IP_N, IP_P, eDeg, nb1, nb2, nb, J
      integer NumberAllTwo, NumberBoundary, NumberPathological
      integer SatMaxDeg, pos, MaxIEcont, eDegVert
      integer IsPointClassicalBoundary(np_total)
      integer IsPointPathologicalBoundary(np_total)
      integer NumberContainedEdges(np_total)
      integer ListDegVertBound(np_total)
      integer,allocatable :: ListDegEdgeBound(:)
      integer idxEdgeBound, idx
      integer iEdgeBound, MaxNbContEdge
      integer, allocatable :: MappingIP_iEdgeBound(:,:)
      integer jEdgeBound, eRealDeg, eRealDegBound
      integer, allocatable :: ListAdjVert(:), ListAdjVertBound(:), ListMatchVert(:)
      integer idxDeg, idxDegBound, nbContIE, idxIE, jVertBound
      real(rkind), allocatable :: ListAng(:)
      logical WeDoOperation, WeWork
      integer iVertBound, IP1, IP_C
      real(rkind) DeltaX, DeltaY, eAng, eAngLargest, eAngSmallest, eDiffAng
      integer idxLargest, idxSmallest, iEdgeBoundFound, IPwork1, IPwork2
      integer, allocatable :: ListPrevVertBound(:), ListNextVertBound(:)
      integer eDegEdgeBound, iEdge, jVertBoundFound, NbCycle
      integer iEdgeBoundFirst, iEdgeBoundWork, IPfound, len
      integer IP_prev, IP_next, IP2
      integer, allocatable :: ListIedgeBound(:), LenCycle(:)
      integer iEdgeBoundNext
      ContElements=0
      DO IE=1,NE_TOTAL
        DO I=1,3
          IP=INEtotal(I,IE)
          ContElements(IP)=ContElements(IP)+1
        END DO
      END DO
      MaxIEcont=maxval(ContElements)
      SatMaxDeg=2*MaxIEcont
      allocate(ListAdjWithDupl(SatMaxDeg,NP_TOTAL))
      allocate(IEcontain(MaxIEcont,NP_TOTAL))
      ListDegWork=0
      DO IE=1,NE_TOTAL
        DO I=1,3
          IF (I.eq.3) THEN
            INEXT=1
          ELSE
            INEXT=I+1
          END IF
          IF (I.eq.1) THEN
            IPREV=3
          ELSE
            IPREV=I-1
          END IF
          IP=INEtotal(I,IE)
          IP_N=INEtotal(INEXT,IE)
          IP_P=INEtotal(IPREV,IE)
          pos=ListDegWork(IP)
          ListAdjWithDupl(2*pos+1,IP)=IP_N
          ListAdjWithDupl(2*pos+2,IP)=IP_P
          IF ((IP.eq.IP_N).or.(IP.eq.IP_P)) THEN
            WRITE(DBG%FHNDL,*) 'IE=', IE
            WRITE(DBG%FHNDL,*) 'I=', I, 'IP=', IP
            WRITE(DBG%FHNDL,*) 'INEXT=', INEXT, ' IP_N=', IP_N
            WRITE(DBG%FHNDL,*) 'IPREV=', IPREV, ' IP_P=', IP_P
            CALL WWM_ABORT("logical error")
          END IF
          IEcontain(pos+1,IP)=IE
          ListDegWork(IP)=pos+1
        END DO
      END DO
      WRITE(STAT%FHNDL,*) 'Stage 1 finished'
      FLUSH(STAT%FHNDL)
      allocate(StatusAdj(SatMaxDeg))
      allocate(TheBound % IOBP(np_total))
      NumberAllTwo=0
      NumberBoundary=0
      NumberPathological=0
      TheBound % nbEdgeBound=0
      TheBound % nbVertBound=0
      IsPointClassicalBoundary=0
      IsPointPathologicalBoundary=0
      DO IP=1,NP_TOTAL
        eDeg=ListDegWork(IP)
        StatusAdj=0
        nb1=0
        nb2=0
        eDegVert=2*eDeg
        DO I=1,eDegVert
          IF (StatusAdj(I) .eq. 0) THEN
            IPadj=ListAdjWithDupl(I,IP)
            nb=0
            DO J=I,eDegVert
              IF (ListAdjWithDupl(J,IP) .eq. IPadj) THEN
                nb=nb+1
                StatusAdj(J)=1
              END IF
            END DO
            IF (nb .eq. 0) CALL WWM_ABORT("Clear bug in code")
            IF (nb .gt. 2) THEN
              WRITE(DBG%FHNDL,*) 'IP=', IP, 'IPadj=', IPadj
              DO J=1,eDeg
                IE=IEcontain(J,IP)
                WRITE(DBG%FHNDL,*) 'IE=', IE
                WRITE(DBG%FHNDL,*) 'INE=', INEtotal(1,IE), INEtotal(2,IE), INEtotal(3,IE)
              END DO
              CALL WWM_ABORT("Hopelessly pathological grid")
            END IF
            IF (nb .eq. 1) nb1=nb1+1
            IF (nb .eq. 2) nb2=nb2+1
            IF ((nb .eq. 1).and.(IP .gt. IPadj)) TheBound % nbEdgeBound = TheBound % nbEdgeBound + 1
          END IF
        END DO
        IF (nb1 .eq. 0) NumberAllTwo=NumberAllTwo + 1
        IF (nb1 .eq. 1) CALL WWM_ABORT("Number 1 should not happen")
        IF (nb1 .eq. 2) THEN
          NumberBoundary=NumberBoundary + 1
          IsPointClassicalBoundary(IP)=1
        END IF
        IF (nb1 .gt. 2) THEN
          NumberPathological=NumberPathological + 1
          IsPointPathologicalBoundary(IP)=1
        END IF
        IF (nb1 .gt. 0) THEN
          TheBound % IOBP(IP)=1
          TheBound % nbVertBound=TheBound % nbVertBound+1
        END IF
      END DO
      WRITE(STAT%FHNDL,*) 'NumberAllTwo      =', NumberAllTwo
      WRITE(STAT%FHNDL,*) 'NumberBoundary    =', NumberBoundary
      WRITE(STAT%FHNDL,*) 'NumberPathological=', NumberPathological
      FLUSH(STAT%FHNDL)
      WRITE(STAT%FHNDL,*) 'Stage 2 finished'
      FLUSH(STAT%FHNDL)
      allocate(TheBound % ListBoundEdge(2, TheBound % nbEdgeBound))
      idxEdgeBound=0
      TheBound % IOBP = 0
      DO IP=1,NP_TOTAL
        eDeg=ListDegWork(IP)
        StatusAdj=0
        eDegVert=2*eDeg
        DO I=1,eDegVert
          IF (StatusAdj(I) .eq. 0) THEN
            IPadj=ListAdjWithDupl(I,IP)
            nb=0
            DO J=I,eDegVert
              IF (ListAdjWithDupl(J,IP) .eq. IPadj) THEN
                nb=nb+1
                StatusAdj(J)=1
              END IF
            END DO
            IF ((nb .eq. 1).and.(IP .gt. IPadj)) THEN
              idxEdgeBound=idxEdgeBound+1
              TheBound % ListBoundEdge(1,idxEdgeBound)=IPadj
              TheBound % ListBoundEdge(2,idxEdgeBound)=IP
            END IF
          END IF
        END DO
      END DO
      WRITE(STAT%FHNDL,*) 'Stage 3 finished'
      FLUSH(STAT%FHNDL)
      allocate(TheBound % ListVertBound(TheBound % nbVertBound))
      idx=0
      DO IP=1,np_total
        IF (TheBound % IOBP(IP) .eq. 1) THEN
          idx=idx+1
          TheBound % ListVertBound(idx)=IP
        END IF
      END DO

      NumberContainedEdges=0
!      Print *, ' maxval(NumberContainedEdges)=', maxval(NumberContainedEdges)
      DO iEdgeBound=1,TheBound % nbEdgeBound
        DO I=1,2
          IP=TheBound % ListBoundEdge(I,iEdgeBound)
          NumberContainedEdges(IP)=NumberContainedEdges(IP)+1
        END DO
!        Print *, 'iEdgeBound=', iEdgeBound, ' max=', maxval(NumberContainedEdges)
      END DO
      MaxNbContEdge=maxval(NumberContainedEdges)
      ListDegVertBound=0
!      Print *, 'nbEdgeBound=', TheBound % nbEdgeBound
!      Print *, 'MaxNbContEdge=', MaxNbContEdge, ' np_total=', np_total
      allocate(MappingIP_iEdgeBound(MaxNbContEdge,np_total))
      DO iEdgeBound=1,TheBound % nbEdgeBound
        DO I=1,2
          IP=TheBound % ListBoundEdge(I,iEdgeBound)
          eDeg=ListDegVertBound(IP)
          MappingIP_iEdgeBound(eDeg+1,IP)=iEdgeBound
          ListDegVertBound(IP)=eDeg+1
        END DO
      END DO
      WRITE(STAT%FHNDL,*) 'Stage 4 finished'
      FLUSH(STAT%FHNDL)
      allocate(TheBound % AdjacencyEdgeBound(2,TheBound % nbEdgeBound))
      allocate(ListDegEdgeBound(TheBound % nbEdgeBound))
      ListDegEdgeBound=0
      DO IP=1,np_total
        IF (IsPointClassicalBoundary(IP) .eq. 1) THEN
          iEdgeBound=MappingIP_iEdgeBound(1,IP)
          jEdgeBound=MappingIP_iEdgeBound(2,IP)
          !
          eDeg=ListDegEdgeBound(iEdgeBound)
          TheBound % AdjacencyEdgeBound(eDeg+1,iEdgeBound)=jEdgeBound
          ListDegEdgeBound(iEdgeBound)=eDeg+1
          !
          eDeg=ListDegEdgeBound(jEdgeBound)
          TheBound % AdjacencyEdgeBound(eDeg+1,jEdgeBound)=iEdgeBound
          ListDegEdgeBound(jEdgeBound)=eDeg+1
        END IF
      END DO
      DO IP=1,np_total
        IF (IsPointPathologicalBoundary(IP) .eq. 1) THEN
          eDeg=ListDegWork(IP)
          eDegVert=2*eDeg
          StatusAdj=0
          eRealDeg=0
          eRealDegBound=0
          DO I=1,eDegVert
            IF (StatusAdj(I) .eq. 0) THEN
              IPadj=ListAdjWithDupl(I,IP)
              nb=0
              DO J=I,eDegVert
                IF (ListAdjWithDupl(J,IP) .eq. IPadj) THEN
                  nb=nb+1
                  StatusAdj(J)=1
                END IF
              END DO
              eRealDeg=eRealDeg+1
              IF (nb .eq. 1) eRealDegBound=eRealDegBound+1
            END IF
          END DO
          allocate(ListAdjVert(eRealDeg))
          allocate(ListAdjVertBound(eRealDegBound))
          StatusAdj=0
          idxDeg=0
          idxDegBound=0
          DO I=1,eDegVert
            IF (StatusAdj(I) .eq. 0) THEN
              IPadj=ListAdjWithDupl(I,IP)
              IF (IPadj .eq. 0) CALL WWM_ABORT("IPadj should be non-zero") 
              nb=0
              DO J=I,eDegVert
                IF (ListAdjWithDupl(J,IP) .eq. IPadj) THEN
                  nb=nb+1
                  StatusAdj(J)=1
                END IF
              END DO
              idxDeg=idxDeg+1
              ListAdjVert(idxDeg)=IPadj
              IF (nb .eq. 1) THEN
                idxDegBound=idxDegBound+1
                ListAdjVertBound(idxDegBound)=IPadj
              END IF
            END IF
          END DO
          if (idxDegBound .ne. eRealDegBound) CALL WWM_ABORT("Logical error on eRealDegBound")
          nbContIE=ListDegWork(IP)
          allocate(ListMatchVert(eRealDegBound))
          ListMatchVert=ListAdjVertBound
          DO 
            WeDoOperation=.false.
            DO iVertBound=1,eRealDegBound
              IP1=ListMatchVert(iVertBound)
              DO idxIE=1,nbContIE
                IE=IEcontain(idxIE,IP)
                DO I=1,3
                  IF (I.eq.3) THEN
                    INEXT=1
                  ELSE
                    INEXT=I+1
                  END IF
                  IF (I.eq.1) THEN
                    IPREV=3
                  ELSE
                    IPREV=I-1
                  END IF
                  IP_C=INEtotal(I,IE)
                  IP_N=INEtotal(INEXT,IE)
                  IP_P=INEtotal(IPREV,IE)
                  IF ((IP_C.eq.IP).and.(IP_N.eq.IP1)) THEN
                    WeDoOperation=.true.
                    ListMatchVert(iVertBound)=IP_P
                  END IF
                END DO
              END DO
            END DO
            IF (WeDoOperation .eqv. .false.) THEN
              EXIT
            END IF
          END DO
          allocate(ListAng(eRealDegBound))
          DO iVertBound=1,eRealDegBound
            IPadj=ListAdjVertBound(iVertBound)
            DeltaX=XPtotal(IPadj) - XPtotal(IP)
            DeltaY=YPtotal(IPadj) - YPtotal(IP)
            eAng=ATAN2(DeltaY, DeltaX)
            ListAng(iVertBound)=eAng
          END DO
          allocate(ListPrevVertBound(eRealDegBound), ListNextVertBound(eRealDegBound))
!          Print *, 'Before determination of next/prev s'
          DO iVertBound=1,eRealDegBound
!            Print *, 'iVertBound=', iVertBound
            idxLargest=0
            idxSmallest=0
            eAngLargest=-1.0
            eAngSmallest=400.0
            DO jVertBound=1,eRealDegBound
              IF (jVertBound .ne. iVertBound) THEN
                eDiffAng=ListAng(jVertBound) - ListAng(iVertBound)
                if (eDiffAng .lt. 0) eDiffAng = eDiffAng + PI2
                IF (eDiffAng .gt. eAngLargest) THEN
                  eAngLargest=eDiffAng
                  idxLargest=jVertBound
                END IF
                IF (eDiffAng .lt. eAngSmallest) THEN
                  eAngSmallest=eDiffAng
                  idxSmallest=jVertBound
                END IF
              END IF
            END DO
            ListPrevVertBound(iVertBound)=ListAdjVertBound(idxSmallest)
            ListNextVertBound(iVertBound)=ListAdjVertBound(idxLargest)
!            Print *, 'IP=', ListAdjVertBound(iVertBound), '   match=', ListMatchVert(iVertBound)
!            Print *, 'iV=', iVertBound, ' prev/next=', ListPrevVertBound(iVertBound), ListNextVertBound(iVertBound)
          END DO
          deallocate(ListAng)
          allocate(ListIedgeBound(eRealDegBound))
          DO iVertBound=1,eRealDegBound
            IP1=ListAdjVertBound(iVertBound)
            IF (IP .lt. IP1) THEN
              IPwork1=IP
              IPwork2=IP1
            ELSE
              IPwork1=IP1
              IPwork2=IP
            END IF
            iEdgeBoundFound=-1
            eDegEdgeBound=ListDegVertBound(IP)
            DO iEdge=1,eDegEdgeBound
              iEdgeBound=MappingIP_iEdgeBound(iEdge,IP)
              IF ((TheBound % ListBoundEdge(1,iEdgeBound) .eq. IPwork1).and.(TheBound % ListBoundEdge(2,iEdgeBound).eq.IPwork2)) THEN
                iEdgeBoundFound=iEdgeBound
              END IF
            END DO
            IF (iEdgeBoundFound .eq. -1) CALL WWM_ABORT("Error in finding iEdgeBoundFound")
            ListIedgeBound(iVertBound)=iEdgeBoundFound
          END DO
!          Print *, 'eRealDegBound=', eRealDegBound
          DO iVertBound=1,eRealDegBound
            IP1=ListAdjVertBound(iVertBound)
            IP2=ListMatchVert(iVertBound)
            IF (IP1 .ne. IP2) THEN
              IP_prev=ListPrevVertBound(iVertBound)
              IP_next=ListNextVertBound(iVertBound)
              IPadj=0
              IF (IP2 .eq. IP_prev) IPadj=IP_next
              IF (IP2 .eq. IP_next) IPadj=IP_prev
              IF (IPadj .eq. 0) CALL WWM_ABORT("Clear bug to solve")
              jVertBoundFound=-1
              DO jVertBound=1,eRealDegBound
                IF (ListAdjVertBound(jVertBound) .eq. IPadj) jVertBoundFound=jVertBound
              END DO
              IF (jVertBoundFound .eq. -1) CALL WWM_ABORT("Another error")
              jVertBound=jVertBoundFound
              iEdgeBound=ListIedgeBound(iVertBound)
              jEdgeBound=ListIedgeBound(jVertBound)
              !
              eDeg=ListDegEdgeBound(iEdgeBound)
              TheBound % AdjacencyEdgeBound(eDeg+1,iEdgeBound)=jEdgeBound
              ListDegEdgeBound(iEdgeBound)=eDeg+1
              !
              eDeg=ListDegEdgeBound(jEdgeBound)
              TheBound % AdjacencyEdgeBound(eDeg+1,jEdgeBound)=iEdgeBound
              ListDegEdgeBound(jEdgeBound)=eDeg+1
            END IF
          END DO
          deallocate(ListIedgeBound)
          deallocate(ListMatchVert)
          deallocate(ListPrevVertBound, ListNextVertBound)
          deallocate(ListAdjVert, ListAdjVertBound)
        END IF
      END DO
      WRITE(STAT%FHNDL,*) 'Stage 5 finished'
      FLUSH(STAT%FHNDL)
      DO iEdgeBound=1,TheBound % nbEdgeBound
        WRITE(DBG%FHNDL,*) 'iEdgeBound/eDeg=', iEdgeBound, ListDegEdgeBound(iEdgeBound)
      END DO
      DO iEdgeBound=1,TheBound % nbEdgeBound
        IF (ListDegEdgeBound(iEdgeBound) .ne. 2) THEN
          WRITE(DBG%FHNDL,*) 'iEdgeBound=', iEdgeBound
          WRITE(DBG%FHNDL,*) 'eDeg=', ListDegEdgeBound(iEdgeBound)
          CALL WWM_ABORT("Degree error")
        END IF
      END DO
      allocate(TheBound % NEIGHBORedge(TheBound % nbEdgeBound))
      allocate(TheBound % CorrespVertex(TheBound % nbEdgeBound))
      TheBound % NEIGHBORedge=0
      allocate(LenCycle(TheBound % nbEdgeBound))
      allocate(TheBound % TheCycleBelong(TheBound % nbEdgeBound))
      LenCycle=0
      NbCycle=0
      DO iEdgeBound=1,TheBound % nbEdgeBound
        IF (TheBound % NEIGHBORedge(iEdgeBound) .eq. 0) THEN
          NbCycle = NbCycle + 1
          iEdgeBoundFirst=iEdgeBound
          iEdgeBoundWork=iEdgeBound 
          iEdgeBoundNext=TheBound % AdjacencyEdgeBound(1,iEdgeBound)
          len=0
          DO
            TheBound % NEIGHBORedge(iEdgeBoundWork)=iEdgeBoundNext
            IPfound=-1
            DO I=1,2
              IP=TheBound % ListBoundEdge(I,iEdgeBoundWork)
              DO J=1,2
                IP2=TheBound % ListBoundEdge(J,iEdgeBoundNext)
                IF (IP2 .eq. IP) THEN
                  IPfound=IP
                END IF
              END DO
            END DO
            IF (IPfound .eq. -1) CALL WWM_ABORT("Error in looking for IPfound")
            TheBound % CorrespVertex(iEdgeBoundWork)=IPfound
            TheBound % TheCycleBelong(iEdgeBoundWork)=NbCycle
            len=len+1 
            WeWork=.false.
            DO I=1,2
              IF (Wework .eqv. .false.) THEN
                jEdgeBound=TheBound % AdjacencyEdgeBound(I,iEdgeBoundNext)
                IF (jEdgeBound .ne. iEdgeBoundWork) THEN
                  WeWork=.true.
                  iEdgeBoundWork=iEdgeBoundNext
                  iEdgeBoundNext=jEdgeBound
                END IF
              END IF
            END DO
            if (iEdgeBoundNext .eq. iEdgeBoundFirst) EXIT
          END DO
          LenCycle(NbCycle)=len
        END IF
      END DO
      TheBound % NbCycle=NbCycle
      allocate(TheBound % LenCycle(NbCycle))
      DO i=1,NbCycle
        TheBound % LenCycle(i)=LenCycle(i)
      END DO
      deallocate(StatusAdj, ListAdjWithDupl)
      deallocate(MappingIP_iEdgeBound, ListDegEdgeBound)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DeallocateBoundaryInfo(TheBound)
      USE DATAPOOL
      IMPLICIT NONE
      type(BoundaryInfo), intent(inout) :: TheBound
      deallocate(TheBound % ListVertBound)
      deallocate(TheBound % ListBoundEdge)
      deallocate(TheBound % AdjacencyEdgeBound)
      deallocate(TheBound % NEIGHBORedge)
      deallocate(TheBound % CorrespVertex)
      deallocate(TheBound % TheCycleBelong)
      deallocate(TheBound % LenCycle)
      deallocate(TheBound % IOBP)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SERIAL_WRITE_BOUNDARY(ncid, np_glob, ne_glob, INEglob, Oper)
      USE NETCDF
      USE DATAPOOL
      implicit none
      integer, intent(in) :: ncid, Oper
      integer, intent(in) :: np_glob, ne_glob
      integer, intent(in) :: INEglob(ne_glob,3)
      integer iret, var_id, two_dims
      integer nbVertBound_dims, nbEdgeBound_dims, NbCycle_dims
      integer idx
      character (len = *), parameter :: CallFct="SERIAL_WRITE_BOUNDARY"
      character (len = *), parameter :: FULLNAME = "full-name"
      type(BoundaryInfo) TheBound
!      Print *, 'Entering SERIAL_WRITE_BOUNDARY'
      CALL SERIAL_GET_BOUNDARY_NEXTGENERATION(TheBound)
!      Print *, 'After SERIAL_GET_BOUNDARY_NEXTGEN, Oper=', Oper
      IF ((Oper == 1).and.(TheBound % nbEdgeBound.gt.0)) THEN
        iret=nf90_inq_dimid(ncid, 'two', two_dims)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, iret)
        !
        iret = nf90_def_dim(ncid, 'nbEdgeBound', TheBound % nbEdgeBound, nbEdgeBound_dims)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, iret)
        iret = nf90_def_dim(ncid, 'nbVertBound', TheBound % nbVertBound, nbVertBound_dims)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, iret)
        iret = nf90_def_dim(ncid, 'NbCycle', TheBound % NbCycle, NbCycle_dims)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, iret)
        !
        iret=nf90_def_var(ncid,'ListVertBound',NF90_INT,(/ nbVertBound_dims/),var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, iret)
        iret=nf90_put_att(ncid,var_id,FULLNAME,'IP of boundary element')
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, iret)
        !
        iret=nf90_def_var(ncid,'ListBoundEdge',NF90_INT,(/ two_dims, nbEdgeBound_dims/),var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 7, iret)
        iret=nf90_put_att(ncid,var_id,FULLNAME,'boundary edges')
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 8, iret)
        !
        iret=nf90_def_var(ncid,'AdjacencyEdgeBound',NF90_INT,(/ two_dims, nbEdgeBound_dims/),var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 9, iret)
        iret=nf90_put_att(ncid,var_id,FULLNAME,'boundary edges')
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 10, iret)
        !
        iret=nf90_def_var(ncid,'NEIGHBORedge',NF90_INT,(/ nbEdgeBound_dims/),var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, iret)
        iret=nf90_put_att(ncid,var_id,FULLNAME,'next boundary edge in the cycle')
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 12, iret)
        !
        iret=nf90_def_var(ncid,'CorrespVertex',NF90_INT,(/ nbEdgeBound_dims/),var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 13, iret)
        iret=nf90_put_att(ncid,var_id,FULLNAME,'Corresponding vertex to the boundary edge')
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 14, iret)
        !
        iret=nf90_def_var(ncid,'CycleBelong',NF90_INT,(/ nbEdgeBound_dims/),var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 15, iret)
        iret=nf90_put_att(ncid,var_id,FULLNAME,'The cycle to which the boundary edge belongs')
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 16, iret)
        !
        iret=nf90_def_var(ncid,'LenCycle',NF90_INT,(/ NbCycle_dims/),var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 17, iret)
        iret=nf90_put_att(ncid,var_id,FULLNAME,'Length of cycles')
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 18, iret)
      END IF
!      Print *, 'After Oper=1'
!      Print *, 'TheBound % nbEdgeBound=', TheBound % nbEdgeBound
      IF ((Oper == 2).and.(TheBound % nbEdgeBound.gt.0)) THEN
        iret=nf90_inq_varid(ncid, "ListVertBound", var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 19, iret)
        iret=nf90_put_var(ncid,var_id,TheBound % ListVertBound, start = (/1/), count = (/ TheBound % nbVertBound/))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 20, iret)
        !
        iret=nf90_inq_varid(ncid, "ListBoundEdge", var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 21, iret)
        iret=nf90_put_var(ncid,var_id,TheBound % ListBoundEdge, start = (/1, 1/), count = (/ 2, TheBound % nbEdgeBound/))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 22, iret)
        !
        iret=nf90_inq_varid(ncid, "AdjacencyEdgeBound", var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 23, iret)
        iret=nf90_put_var(ncid,var_id,TheBound % AdjacencyEdgeBound, start = (/1, 1/), count = (/ 2, TheBound % nbEdgeBound/))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 24, iret)
        !
        iret=nf90_inq_varid(ncid, "NEIGHBORedge", var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 25, iret)
        iret=nf90_put_var(ncid,var_id,TheBound % NEIGHBORedge, start = (/1/), count = (/ TheBound % nbEdgeBound/))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 26, iret)
        !
        iret=nf90_inq_varid(ncid, "CorrespVertex", var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 27, iret)
        iret=nf90_put_var(ncid,var_id,TheBound % CorrespVertex, start = (/1/), count = (/ TheBound % nbEdgeBound/))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 28, iret)
        !
        iret=nf90_inq_varid(ncid, "CycleBelong", var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 29, iret)
        iret=nf90_put_var(ncid,var_id,TheBound % TheCycleBelong, start = (/1/), count = (/ TheBound % nbEdgeBound/))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 30, iret)
        !
        iret=nf90_inq_varid(ncid, "LenCycle", var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 31, iret)
        iret=nf90_put_var(ncid,var_id,TheBound % LenCycle, start = (/1/), count = (/ TheBound % NbCycle/))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 32, iret)
      END IF
!      Print *, 'After Oper=2'
      CALL DeallocateBoundaryInfo(TheBound)
!      Print *, 'After deallocate'
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE REPORT_ERROR_INQ(iret, str)
      INTEGER, intent(in) :: iret
      character(*) :: str
      character(256) :: ErrMsg
      IF (iret.ne.0) THEN
        WRITE (ErrMsg,10) 'ERROR while inquiring', str
  10    FORMAT (a,' ',a)
        CALL WWM_ABORT(TRIM(ErrMsg))
      ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE REPORT_ERROR_DEF(iret, str)
      INTEGER, intent(in) :: iret
      character(*) :: str
      character(256) :: ErrMsg
      IF (iret.ne.0) THEN
        WRITE (ErrMsg,10) 'ERROR while defining', str
  10    FORMAT (a,' ',a)
        CALL WWM_ABORT(TRIM(ErrMsg))
      ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GENERIC_NETCDF_ERROR_WWM_CLEAR(ncid, CallFct, idx, iret)
      USE NETCDF
      USE DATAPOOL, only : wwmerr
      implicit none
      integer, intent(in) :: iret, idx, ncid
      character(*), intent(in) :: CallFct
      character(len=500) :: CHRERR
      integer iretB
      IF (iret .NE. nf90_noerr) THEN
        CHRERR = nf90_strerror(iret)
        WRITE(wwmerr,*) TRIM(CallFct), ' -', idx, '-', CHRERR
        iretB=nf90_close(ncid)
        Print *, 'iretB=', iretB
        CALL WWM_ABORT(wwmerr)
      ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GENERIC_NETCDF_ERROR_WWM(CallFct, idx, iret)
      USE NETCDF
      USE DATAPOOL, only : wwmerr
      implicit none
      integer, intent(in) :: iret, idx
      character(*), intent(in) :: CallFct
      character(len=500) :: CHRERR
      IF (iret .NE. nf90_noerr) THEN
        CHRERR = nf90_strerror(iret)
        WRITE(wwmerr,*) 'NETCDF error in routine ', TRIM(CallFct), ' Error Message: ', TRIM(CHRERR), ' Position in the routine :', idx
        CALL WWM_ABORT(wwmerr)
      ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE NAMEVARIABLE(idx, eStr, eStrFullName, eStrUnit)
      INTEGER, intent(in) :: idx
      character(len=40), intent(out) :: eStr
      character(len=80), intent(out) :: eStrFullName
      character(len=40), intent(out) :: eStrUnit

      IF (IDX.eq.1) THEN
        eStr="HS"
        eStrFullName="Significant wave height"
        eStrUnit="meter"
      ELSE IF (IDX.eq.2) THEN
        eStr="TM01"
        eStrFullName="mean wave period"
        eStrUnit="second"
      ELSE IF (IDX.eq.3) THEN
        eStr="TM02"
        eStrFullName="Zero-crossing wave period"
        eStrUnit="second"
      ELSE IF (IDX.eq.4) THEN
        eStr="TM10"
        eStrFullName="Mean period of wave over topping/run-up"
        eStrUnit="second"
      ELSE IF (IDX.eq.5) THEN
        eStr="KLM"
        eStrFullName="mean wave number"
        eStrUnit="meter-1"
      ELSE IF (IDX.eq.6) THEN
        eStr="WLM"
        eStrFullName="Mean wave length"
        eStrUnit="meter"
      ELSE IF (IDX.eq.7) THEN
        eStr="ETOTC"
        eStrFullName="model variable"
        eStrUnit="unk"
      ELSE IF (IDX.eq.8) THEN
        eStr="ETOTS"
        eStrFullName="model variable"
        eStrUnit="unk"
      ELSE IF (IDX.eq.9) THEN
        eStr="DM"
        eStrFullName="Mean wave direction"
        eStrUnit="degree"
      ELSE IF (IDX.eq.10) THEN
        eStr="DSPR"
        eStrFullName="Directional spreading"
        eStrUnit="degree"
      ELSE IF (IDX.eq.11) THEN
        eStr="TPPD"
        eStrFullName="Discrete peak wave period"
        eStrUnit="second"
      ELSE IF (IDX.eq.12) THEN
        eStr="CPPD"
        eStrFullName="Discrete peak wave speed"
        eStrUnit="meter second -1"
      ELSE IF (IDX.eq.13) THEN
        eStr="KPPD"
        eStrFullName="discrete peak wave number"
        eStrUnit="meter-1"
      ELSE IF (IDX.eq.14) THEN
        eStr="CGPD"
        eStrFullName="discrete peak group speed"
        eStrUnit="meter second-1"
      ELSE IF (IDX.eq.15) THEN
        eStr="TPP"
        eStrFullName="peak wave period"
        eStrUnit="second"
      ELSE IF (IDX.eq.16) THEN
        eStr="CPP"
        eStrFullName="peak wave speed"
        eStrUnit="meter second-1"
      ELSE IF (IDX.eq.17) THEN
        eStr="WNPP"
        eStrFullName="unk"
        eStrUnit="unk"
      ELSE IF (IDX.eq.18) THEN
        eStr="CGPP"
        eStrFullName="peak group velocity"
        eStrUnit="meter second-1"
      ELSE IF (IDX.eq.19) THEN
        eStr="KPP"
        eStrFullName="peak wave number"
        eStrUnit="meter-1"
      ELSE IF (IDX.eq.20) THEN
        eStr="LPP"
        eStrFullName="Peak wave length"
        eStrUnit="meter"
      ELSE IF (IDX.eq.21) THEN
        eStr="PEAKD"
        eStrFullName="Peak wave direction"
        eStrUnit="degree"
      ELSE IF (IDX.eq.22) THEN
        eStr="PEAKDSPR"
        eStrFullName="Peak directional spreading"
        eStrUnit="degree"
      ELSE IF (IDX.eq.23) THEN
        eStr="DPEAK"
        eStrFullName="discrete peak direction"
        eStrUnit="degree"
      ELSE IF (IDX.eq.24) THEN
        eStr="UBOT"
        eStrFullName="wind-induced bottom orbital velocity"
        eStrUnit="meter second-1"
      ELSE IF (IDX.eq.25) THEN
        eStr="ORBITAL"
        eStrFullName="unk"
        eStrUnit="unk"
      ELSE IF (IDX.eq.26) THEN
        eStr="BOTEXPER"
        eStrFullName="unk"
        eStrUnit="unk"
      ELSE IF (IDX.eq.27) THEN
        eStr="TMBOT"
        eStrFullName="unk"
        eStrUnit="unk"
      ELSE IF (IDX.eq.28) THEN
        eStr="URSELL"
        eStrFullName="ursell number"
        eStrUnit="non-dimensional"
      ELSE IF (IDX.eq.29) THEN
        eStr="UFRIC"
        eStrFullName="air friction velocity"
        eStrUnit="meter second-1"
      ELSE IF (IDX.eq.30) THEN
        eStr="Z0"
        eStrFullName="air roughness length"
        eStrUnit="meter"
      ELSE IF (IDX.eq.31) THEN
        eStr="ALPHA_CH"
        eStrFullName="air Charnock coefficient"
        eStrUnit="non-dimensional"
      ELSE IF (IDX.eq.32) THEN
        eStr="Uwind"
        eStrFullName="wind in X direction"
        eStrUnit="meter second-1"
      ELSE IF (IDX.eq.33) THEN
        eStr="Vwind"
        eStrFullName="wind in Y direction"
        eStrUnit="meter second-1"
      ELSE IF (IDX.eq.34) THEN
        eStr="CD"
        eStrFullName="drag coefficient"
        eStrUnit="non-dimensional"
      ELSE IF (IDX.eq.35) THEN
        eStr="CURTX"
        eStrFullName="current in X direction"
        eStrUnit="meter second-1"
      ELSE IF (IDX.eq.36) THEN
        eStr="CURTY"
        eStrFullName="current in Y direction"
        eStrUnit="meter second-1"
      ELSE IF (IDX.eq.37) THEN
        eStr="WATLEV"
        eStrFullName="water level"
        eStrUnit="meter"
      ELSE IF (IDX.eq.38) THEN
        eStr="WATLEVOLD"
        eStrFullName="water level at previous time step"
        eStrUnit="meter"
      ELSE IF (IDX.eq.39) THEN
        eStr="DEPDT"
        eStrFullName="water level change"
        eStrUnit="meter second-1"
      ELSE IF (IDX.eq.40) THEN
        eStr="DEP"
        eStrFullName="bathymetry"
        eStrUnit="meter"
      ELSE IF (IDX.eq.41) THEN
        eStr="WINDMAG"
        eStrFullName="10-m wind magnitude"
        eStrUnit="meter second-1"
      ELSE IF (IDX.eq.42) THEN
        eStr="TAUW"
        eStrFullName="wave supported surface stress"
        eStrUnit="unk"
      ELSE IF (IDX.eq.43) THEN
        eStr="TAUWX"
        eStrFullName="wave supported surface X-stress"
        eStrUnit="unk"
      ELSE IF (IDX.eq.44) THEN
        eStr="TAUWY"
        eStrFullName="wave supported surface Y-stress"
        eStrUnit="unk"
      ELSE IF (IDX.eq.45) THEN
        eStr="TAUHF"
        eStrFullName="high frequency surface stress"
        eStrUnit="unk"
      ELSE IF (IDX.eq.46) THEN
        eStr="TAUTOT"
        eStrFullName="total surface stress"
        eStrUnit="unk"
      ELSE IF (IDX.eq.47) THEN
        eStr="STOKESBOTTX"
        eStrFullName="bottom Stokes velocity in X direction"
        eStrUnit="meter second -1"
      ELSE IF (IDX.eq.48) THEN
        eStr="STOKESBOTTY"
        eStrFullName="bottom Stokes velocity in Y direction"
        eStrUnit="meter second -1"
      ELSE IF (IDX.eq.49) THEN
        eStr="STOKESSURFX"
        eStrFullName="surface Stokes velocity in X direction"
        eStrUnit="meter second -1"
      ELSE IF (IDX.eq.50) THEN
        eStr="STOKESSURFY"
        eStrFullName="surface Stokes velocity in Y direction"
        eStrUnit="meter second -1"
      ELSE IF (IDX.eq.51) THEN
        eStr="STOKESBAROX"
        eStrFullName="barotropic Stokes velocity in X direction"
        eStrUnit="meter second -1"
      ELSE IF (IDX.eq.52) THEN
        eStr="STOKESBAROY"
        eStrFullName="barotropic Stokes velocity in Y direction"
        eStrUnit="meter second -1"
      ELSE IF (IDX.eq.53) THEN
        eStr="RSXX"
        eStrFullName="barotropic Stress potential Sxx"
        eStrUnit="meter2 second2"
      ELSE IF (IDX.eq.54) THEN
        eStr="RSXY"
        eStrFullName="barotropic Stress potential Sxy"
        eStrUnit="meter2 second2"
      ELSE IF (IDX.eq.55) THEN
        eStr="RSYY"
        eStrFullName="barotropic Stress potential Syy"
        eStrUnit="meter2 second2"
      ELSE IF (IDX.eq.56) THEN
        eStr="CFL1"
        eStrFullName="CFL number 1"
        eStrUnit="non-dimensional"
      ELSE IF (IDX.eq.57) THEN
        eStr="CFL2"
        eStrFullName="CFL number 2"
        eStrUnit="non-dimensional"
      ELSE IF (IDX.eq.58) THEN
        eStr="CFL3"
        eStrFullName="CFL number 3"
        eStrUnit="non-dimensional"
      ELSE IF (IDX.eq.59) THEN
        eStr="ZETA_SETUP"
        eStrFullName="Free-surface elevation induced setup"
        eStrUnit="m"
      ELSE
        CALL WWM_ABORT('Wrong Number')
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WRITE_NETCDF_TIME_HEADER(ncid, nbTime, ntime_dims)
      USE DATAPOOL
      USE NETCDF
      IMPLICIT NONE
      integer, intent(in) :: ncid, nbTime
      integer, intent(inout) :: ntime_dims
      character (len = *), parameter :: UNITS = "units"
      character (len = *), parameter :: CallFct="WRITE_NETCDF_TIME_HEADER"
      integer iret, fifteen_dims, var_id
      iret=nf90_inq_dimid(ncid, 'fifteen', fifteen_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, iret)
      IF (nbTime.gt.0) THEN
        iret = nf90_def_dim(ncid, 'ocean_time', nbTime, ntime_dims)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, iret)
      ELSE
        iret = nf90_def_dim(ncid, 'ocean_time', NF90_UNLIMITED, ntime_dims)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, iret)
      END IF
      iret=nf90_def_var(ncid,'ocean_time',NF90_RUNTYPE,(/ ntime_dims/), var_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, iret)
      iret=nf90_put_att(ncid,var_id,UNITS,'seconds since 1858-11-17 00:00:00')
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, iret)
      iret=nf90_put_att(ncid,var_id,"calendar",'gregorian')
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, iret)
      !
      iret=nf90_def_var(ncid,'ocean_time_day',NF90_RUNTYPE,(/ ntime_dims/), var_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 7, iret)
      iret=nf90_put_att(ncid,var_id,UNITS,'days since 1858-11-17 00:00:00')
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 8, iret)
      iret=nf90_put_att(ncid,var_id,"calendar",'gregorian')
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 9, iret)
      !
      iret=nf90_def_var(ncid,'ocean_time_str',NF90_CHAR,(/ fifteen_dims, ntime_dims/), var_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 10, iret)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WRITE_NETCDF_TIME(ncid, idx, eTimeDay)
      USE DATAPOOL, only : DAY2SEC,RKIND, wwmerr
      USE NETCDF
      IMPLICIT NONE
      integer, intent(in) :: ncid, idx
      REAL(rkind), intent(IN) :: eTimeDay
      character (len = *), parameter :: CallFct="WRITE_NETCDF_TIME"
      integer oceantimeday_id, oceantimestr_id, oceantime_id
      integer iret, I
      CHARACTER          :: eChar
      REAL(rkind) eTimeSec
      CHARACTER(LEN=15) :: eTimeStr
      !
      CALL MJD2CT(eTimeDay,eTimeStr)
      eTimeSec=eTimeDay*DAY2SEC
      iret=nf90_inq_varid(ncid, 'ocean_time', oceantime_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, iret)
      iret=nf90_put_var(ncid,oceantime_id,eTimeSec,start=(/idx/) )
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, iret)
      !
      iret=nf90_inq_varid(ncid, 'ocean_time_day', oceantimeday_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, iret)
      iret=nf90_put_var(ncid,oceantimeday_id,eTimeDay,start=(/idx/) )
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, iret)
      !
      iret=nf90_inq_varid(ncid, 'ocean_time_str', oceantimestr_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, iret)
      DO i=1,15
        eChar=eTimeStr(i:i)
        iret=nf90_put_var(ncid,oceantimestr_id,eChar,start=(/i, idx/) )
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, iret)
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GET_IOBPD_OUTPUT(IOBPDoutput, np_write)
      USE DATAPOOL, only : MDC, IOBPD, np_total, MULTIPLEOUT_HIS, MNP
# ifdef MPI_PARALL_GRID
      USE datapool, only : iplg, comm, nproc, istatus, ierr, myrank, itype
# endif
      IMPLICIT NONE
      INTEGER, intent(in)  :: np_write
      INTEGER, INTENT(OUT) :: IOBPDoutput(MDC, np_write)
# ifdef MPI_PARALL_GRID
      integer, allocatable :: rIOBPD(:,:), rStatus(:), Status(:)
      integer iProc, IP, istat
      IF (MULTIPLEOUT_HIS .eq. 1) THEN
        IOBPDoutput=IOBPD
      ELSE
        allocate(rIOBPD(MDC, np_total), rStatus(np_total), Status(np_total), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_netcdf, allocate error 7')
        IOBPDoutput=0
        Status=0
        DO IP=1,MNP
          IOBPDoutput(:, iplg(IP))=IOBPD(:,IP)
          Status(iplg(IP)) = 1
        END DO
        IF (myrank .eq. 0) THEN
          DO iProc=2,nproc
            CALL MPI_RECV(rIOBPD, MDC*np_total, itype, iProc-1, 193, comm, istatus, ierr)
            CALL MPI_RECV(rStatus, np_total, itype, iProc-1, 197, comm, istatus, ierr)
            DO IP=1,np_total
              IF (rStatus(IP) .eq. 1) THEN
                IOBPDoutput(:,IP)=rIOBPD(:,IP)
              END IF
            END DO
          END DO
          DO iProc=2,nproc
            CALL MPI_SEND(IOBPDoutput, MDC*np_total, itype, iProc-1, 199, comm, ierr)
          END DO
        ELSE
          CALL MPI_SEND(IOBPDoutput, MDC*np_total, itype, 0, 193, comm, ierr)
          CALL MPI_SEND(Status, np_total, itype, 0, 197, comm, ierr)
          CALL MPI_RECV(IOBPDoutput, MDC*np_total, itype, 0, 199, comm, istatus, ierr)
        ENDIF
        deallocate(rIOBPD, rStatus, Status)
      END IF
# else
      IOBPDoutput=IOBPD
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WRITE_NETCDF_HEADERS_1(ncid, nbTime, MULTIPLEOUT, GRIDWRITE_W, IOBPD_HISTORY_W, np_write, ne_write)
      USE DATAPOOL
      USE NETCDF
      implicit none
      integer, intent(in) :: ncid, nbTime, MULTIPLEOUT
      integer, intent(in) :: np_write, ne_write
      logical, intent(in) :: GRIDWRITE_W, IOBPD_HISTORY_W
      !
      character (len = *), parameter :: UNITS = "units"
      integer one_dims, two_dims, three_dims, fifteen_dims
      integer mnp_dims, mne_dims, nfreq_dims, ndir_dims
# ifdef MPI_PARALL_GRID
      integer np_global_dims, ne_global_dims
# endif
      integer iret, var_id
      integer ntime_dims
      integer p_dims, e_dims
      integer Oper
      character (len = *), parameter :: CallFct="WRITE_NETCDF_HEADERS_1"
      IF ((np_write.eq.0).or.(ne_write.eq.0)) THEN
        CALL WWM_ABORT('np_write=0 or ne_write=0, not allowed by any mean')
      ENDIF
      iret = nf90_def_dim(ncid, 'one', 1, one_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, iret)
      iret = nf90_def_dim(ncid, 'two', 2, two_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, iret)
      iret = nf90_def_dim(ncid, 'three', 3, three_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, iret)
      iret = nf90_def_dim(ncid, 'fifteen', 15, fifteen_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, iret)
      iret = nf90_def_dim(ncid, 'mnp', np_write, mnp_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, iret)
      iret = nf90_def_dim(ncid, 'mne', ne_write, mne_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, iret)
      iret = nf90_def_dim(ncid, 'nfreq', MSC, nfreq_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 7, iret)
      iret = nf90_def_dim(ncid, 'ndir', MDC, ndir_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 8, iret)
# ifdef MPI_PARALL_GRID
      IF (MULTIPLEOUT.eq.1) THEN
        iret = nf90_def_dim(ncid, 'np_global', np_global, np_global_dims)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 9, iret)
        iret = nf90_def_dim(ncid, 'ne_global', ne_global, ne_global_dims)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 9, iret)
        iret = nf90_def_var(ncid,'iplg',NF90_INT,(/ mnp_dims/),var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 10, iret)
        iret = nf90_put_att(ncid,var_id,'description','local to global indexes')
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, iret)
        iret=nf90_def_var(ncid,'nproc',NF90_INT,(/one_dims/),var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 12, iret)
        iret=nf90_put_att(ncid,var_id,UNITS,'integer')
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 13, iret)
        iret=nf90_put_att(ncid,var_id,'description','number of processors')
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 14, iret)
      END IF
      !
      iret=nf90_def_var(ncid,'MULTIPLEOUT',NF90_INT,(/one_dims/),var_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 15, iret)
      iret=nf90_put_att(ncid,var_id,UNITS,'integer')
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 16, iret)
      iret=nf90_put_att(ncid,var_id,'description','multiple status')
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 17, iret)
# endif
      !
      IF (PARAMWRITE_HIS) THEN
        CALL WRITE_PARAM_1(ncid, one_dims)
      ENDIF
      !
      CALL WRITE_NETCDF_TIME_HEADER(ncid, nbTime, ntime_dims)
      !
      IF (GRIDWRITE_W) THEN
# ifdef MPI_PARALL_GRID
        IF (MULTIPLEOUT.eq.1) THEN
          e_dims=ne_global_dims
          p_dims=np_global_dims
        ELSE
          e_dims=mne_dims
          p_dims=mnp_dims
        END IF
# else
        e_dims=mne_dims
        p_dims=mnp_dims
# endif
! element
        iret=nf90_def_var(ncid,'ele',NF90_INT,(/three_dims, e_dims/),var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 27, iret)
        iret=nf90_put_att(ncid,var_id,UNITS,'non-dimensional')
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 28, iret)
! lon
        IF (LSPHE) THEN
          iret=nf90_def_var(ncid,"lon",NF90_RUNTYPE,(/ p_dims/),var_id)
        ELSE
          iret=nf90_def_var(ncid,"x",NF90_RUNTYPE,(/ p_dims/),var_id)
        END IF
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 29, iret)
        IF (LSPHE) THEN
          iret=nf90_put_att(ncid,var_id,UNITS,'degree')
        ELSE
          iret=nf90_put_att(ncid,var_id,UNITS,'meter')
        END IF
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 30, iret)
! lat
        IF (LSPHE) THEN
          iret=nf90_def_var(ncid,"lat",NF90_RUNTYPE,(/ p_dims/),var_id)
        ELSE
          iret=nf90_def_var(ncid,"y",NF90_RUNTYPE,(/ p_dims/),var_id)
        END IF
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 31, iret)
        IF (LSPHE) THEN
          iret=nf90_put_att(ncid,var_id,UNITS,'degree')
        ELSE
          iret=nf90_put_att(ncid,var_id,UNITS,'meter')
        END IF
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 32, iret)
! IOBP
        iret=nf90_def_var(ncid,"IOBP",NF90_RUNTYPE,(/ p_dims/),var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 33, iret)
! depth
        iret=nf90_def_var(ncid,'depth',NF90_RUNTYPE,(/ p_dims/),var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 34, iret)
        iret=nf90_put_att(ncid,var_id,UNITS,'meters')
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 35, iret)
! boundary
        Oper=1
        CALL SERIAL_WRITE_BOUNDARY(ncid, np_total, ne_total, INEtotal, Oper)
        !
      END IF
      IF (IOBPD_HISTORY_W) THEN
# ifdef MPI_PARALL_GRID
        IF (MULTIPLEOUT.eq.1) THEN
          p_dims=np_global_dims
        ELSE
          p_dims=mnp_dims
        END IF
# else
        p_dims=mnp_dims
# endif
        iret=nf90_def_var(ncid,'IOBPD',NF90_INT,(/ ndir_dims, p_dims, ntime_dims/), var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 20, iret)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WRITE_NETCDF_HEADERS_2(FILERET, MULTIPLEOUT, WriteOutputProcess, GRIDWRITE_W, np_write, ne_write)
      USE DATAPOOL
      USE NETCDF
      implicit none
      character(len=140), intent(in) :: FILERET
      integer, intent(in) :: MULTIPLEOUT
      logical, intent(in) :: WriteOutputProcess, GRIDWRITE_W
      integer, intent(in) :: np_write, ne_write
      integer var_id, iret, ncid
      character (len = *), parameter :: CallFct="WRITE_NETCDF_HEADERS_2"
      integer Oper
#ifdef MPI_PARALL_GRID
      integer eInt(1)
#endif
      IF (WriteOutputProcess) THEN
        iret = nf90_open(FILERET, nf90_write, ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, iret)
      END IF
#ifdef MPI_PARALL_GRID
      IF (MULTIPLEOUT.eq.1) THEN
        iret=nf90_inq_varid(ncid, 'iplg', var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, iret)
        iret=nf90_put_var(ncid,var_id,iplg,start=(/1/), count = (/ np_write /))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, iret)
      END IF
      !
      IF (WriteOutputProcess) THEN
        IF (MULTIPLEOUT.eq.1) THEN
          iret=nf90_inq_varid(ncid,'nproc',var_id)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, iret)
          eInt(1)=nproc
          iret=nf90_put_var(ncid,var_id,eInt,start = (/1/), count=(/1/))
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, iret)
        END IF
        !
        iret=nf90_inq_varid(ncid,'MULTIPLEOUT',var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, iret)
        eInt(1)=MULTIPLEOUT
        iret=nf90_put_var(ncid,var_id,eInt,start = (/1/), count=(/1/))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 7, iret)
      END IF
#endif
      IF (PARAMWRITE_HIS.and.WriteOutputProcess) THEN
        CALL WRITE_PARAM_2(ncid)
      END IF
      IF (GRIDWRITE_W.and.WriteOutputProcess) THEN
        !
        iret=nf90_inq_varid(ncid, "ele", var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 8, iret)
        iret=nf90_put_var(ncid,var_id,INEtotal, start = (/1,1/), count = (/ 3, ne_total/))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 9, iret)
        !
        IF (LSPHE) THEN
          iret=nf90_inq_varid(ncid, "lon", var_id)
        ELSE
          iret=nf90_inq_varid(ncid, "x", var_id)
        ENDIF
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 10, iret)
        iret=nf90_put_var(ncid,var_id,XPtotal, start = (/1/), count = (/ np_total/))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, iret)
        !
        IF (LSPHE) THEN
          iret=nf90_inq_varid(ncid, "lat", var_id)
        ELSE
          iret=nf90_inq_varid(ncid, "y", var_id)
        ENDIF
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 12, iret)
        iret=nf90_put_var(ncid,var_id,YPtotal, start = (/1/), count = (/ np_total/))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 13, iret)
        !
        iret=nf90_inq_varid(ncid, "IOBP", var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 14, iret)
        iret=nf90_put_var(ncid,var_id,IOBPtotal, start = (/1/), count = (/ np_total/))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 15, iret)
        !
        iret=nf90_inq_varid(ncid, "depth", var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 16, iret)
        iret=nf90_put_var(ncid,var_id,DEPtotal, start = (/1/), count = (/ np_write/))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 17, iret)
        !
        Oper=2
        CALL SERIAL_WRITE_BOUNDARY(ncid, np_total, ne_total, INEtotal, Oper)
      ENDIF
      IF (WriteOutputProcess) THEN
        iret = nf90_close(ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 18, iret)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WRITE_PARAM_1(ncid, one_dims)
      USE NETCDF
      USE DATAPOOL, only : NF90_RUNTYPE
      implicit none
      integer, intent(in) :: ncid, one_dims
      integer :: iret, var_id
      character (len = *), parameter :: UNITS = "units"
      iret = nf90_def_var(ncid,'frlow', NF90_RUNTYPE,(/ one_dims/), var_id)
      CALL REPORT_ERROR_DEF(iret, 'frlow')
      iret = nf90_put_att(ncid,var_id,UNITS,'lower_frequency')
      !
      iret = nf90_def_var(ncid,'frhigh',NF90_RUNTYPE,(/ one_dims/),var_id)
      CALL REPORT_ERROR_DEF(iret, 'frhigh')
      iret = nf90_put_att(ncid,var_id,UNITS,'higher_frequency')
      !
      iret = nf90_def_var(ncid,'MESNL',NF90_INT,(/ one_dims/),var_id)
      CALL REPORT_ERROR_DEF(iret, 'MESNL')
      iret = nf90_put_att(ncid,var_id,'description','nonlinear interaction nl4')
      !
      iret = nf90_def_var(ncid,'MESIN',NF90_INT,(/ one_dims/),var_id)
      CALL REPORT_ERROR_DEF(iret, 'MESIN')
      iret = nf90_put_att(ncid,var_id,'description','wind input source term')
      !
      iret = nf90_def_var(ncid,'MESDS',NF90_INT,(/ one_dims/),var_id)
      CALL REPORT_ERROR_DEF(iret, 'MESDS')
      iret = nf90_put_att(ncid,var_id,'description','dissipation source term')
      !
      iret = nf90_def_var(ncid,'MESBF',NF90_INT,(/ one_dims/),var_id)
      CALL REPORT_ERROR_DEF(iret, 'MESBF')
      iret = nf90_put_att(ncid,var_id,'description','bottom friction')
      !
      iret = nf90_def_var(ncid,'ICOMP',NF90_INT,(/ one_dims/),var_id)
      CALL REPORT_ERROR_DEF(iret, 'ICOMP')
      iret = nf90_put_att(ncid,var_id,'description','implicitness')
      !
      iret = nf90_def_var(ncid,'AMETHOD',NF90_INT,(/ one_dims/),var_id)
      CALL REPORT_ERROR_DEF(iret, 'AMETHOD')
      iret = nf90_put_att(ncid,var_id,'description','advection method')
      !
      iret = nf90_def_var(ncid,'FMETHOD',NF90_INT,(/ one_dims/),var_id)
      CALL REPORT_ERROR_DEF(iret, 'FMETHOD')
      iret = nf90_put_att(ncid,var_id,'description','frequency shifting method')
      !
      iret = nf90_def_var(ncid,'DMETHOD',NF90_INT,(/ one_dims/),var_id)
      CALL REPORT_ERROR_DEF(iret, 'DMETHOD')
      iret = nf90_put_att(ncid,var_id,'description','directional shifting method(refraction)')
      !
      iret = nf90_def_var(ncid,'SMETHOD',NF90_INT,(/ one_dims/),var_id)
      CALL REPORT_ERROR_DEF(iret, 'SMETHOD')
      iret = nf90_put_att(ncid,var_id,'description','source term integration method')
      !
      iret = nf90_def_var(ncid,'LSPHE',NF90_INT,(/ one_dims/),var_id)
      CALL REPORT_ERROR_DEF(iret, 'LSPHE')
      iret = nf90_put_att(ncid,var_id,'description','spherical coordinates')
      END SUBROUTINE WRITE_PARAM_1
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WRITE_PARAM_2(ncid)
      USE DATAPOOL, only : FRLOW, FRHIGH, ICOMP, AMETHOD, FMETHOD,        &
     &    DMETHOD, SMETHOD, MESIN, MESBF, MESDS, MESNL, LSPHE
      USE NETCDF
      IMPLICIT NONE
      integer, intent(in) :: ncid
      integer iret, var_id
      integer LSPHE_INT
      character (len = *), parameter :: CallFct="WRITE_PARAM_2"
      iret=nf90_inq_varid(ncid, "frlow", var_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, iret)
      iret=nf90_put_var(ncid,var_id,FRLOW,start=(/1/) )
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, iret)
      !
      iret=nf90_inq_varid(ncid, "frhigh", var_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, iret)
      iret=nf90_put_var(ncid,var_id,FRHIGH, start=(/1/) )
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, iret)
      !
      iret=nf90_inq_varid(ncid, "MESIN", var_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 15, iret)
      iret = nf90_put_var(ncid,var_id,MESIN, start=(/1/) )
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 16, iret)
      !
      iret=nf90_inq_varid(ncid, "MESBF", var_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 17, iret)
      iret=nf90_put_var(ncid,var_id,MESBF, start=(/1/) )
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 18, iret)
      !
      iret=nf90_inq_varid(ncid, "MESDS", var_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 19, iret)
      iret=nf90_put_var(ncid,var_id,MESDS, start=(/1/) )
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 20, iret)
      !
      iret=nf90_inq_varid(ncid, "MESNL", var_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 21, iret)
      iret=nf90_put_var(ncid,var_id,MESNL, start=(/1/) )
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 22, iret)
      !
      iret=nf90_inq_varid(ncid, "ICOMP", var_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, iret)
      iret=nf90_put_var(ncid,var_id,ICOMP, start=(/1/) )
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, iret)
      !
      iret=nf90_inq_varid(ncid, "AMETHOD", var_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 7, iret)
      iret=nf90_put_var(ncid,var_id,AMETHOD, start=(/1/) )
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 8, iret)
      !
      iret=nf90_inq_varid(ncid, "FMETHOD", var_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 9, iret)
      iret=nf90_put_var(ncid,var_id,FMETHOD, start=(/1/) )
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 10, iret)
      !
      iret=nf90_inq_varid(ncid, "DMETHOD", var_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, iret)
      iret=nf90_put_var(ncid,var_id,DMETHOD, start=(/1/) )
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 12, iret)
      !
      iret=nf90_inq_varid(ncid, "SMETHOD", var_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 13, iret)
      iret=nf90_put_var(ncid,var_id,SMETHOD, start=(/1/) )
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 14, iret)
      !
      iret=nf90_inq_varid(ncid, "LSPHE", var_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 13, iret)
      IF (LSPHE) THEN
        LSPHE_INT=1
      ELSE
        LSPHE_INT=0
      END IF
      iret=nf90_put_var(ncid,var_id,LSPHE_INT, start=(/1/) )
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 14, iret)
      END SUBROUTINE WRITE_PARAM_2
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DEFINE_STATION_NC(FILE_NAME, MULTIPLEOUT)
      USE NETCDF
      USE DATAPOOL
      implicit none
      character(len=256), intent(in) :: FILE_NAME
      integer, intent(in) :: MULTIPLEOUT
      character (len = *), parameter :: CallFct="DEFINE_STATION_NC"
      character (len = *), parameter :: UNITS = "units"
      character (len = *), parameter :: FULLNAME = "full-name"
      character(len=40) :: eStr, eStrUnit
      character(len=80) :: eStrFullName
      integer iret, ncid, nbstat_dims, ntime_dims, nfreq_dims, ndir_dims
      integer one_dims, three_dims, var_id, I
      iret = nf90_create(TRIM(FILE_NAME), NF90_CLOBBER, ncid)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, iret)

      CALL WRITE_NETCDF_HEADERS_STAT_1(ncid, -1, MULTIPLEOUT)

      iret=nf90_inq_dimid(ncid, 'nbstation', nbstat_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, iret)

      iret=nf90_inq_dimid(ncid, 'ocean_time', ntime_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, iret)

      iret=nf90_inq_dimid(ncid, 'nfreq', nfreq_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, iret)

      iret=nf90_inq_dimid(ncid, 'ndir', ndir_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, iret)

      iret=nf90_inq_dimid(ncid, 'one', one_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, iret)

      iret=nf90_inq_dimid(ncid, 'three', three_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 7, iret)

      IF (VAROUT_STATION%AC) THEN
        iret=nf90_def_var(ncid,'AC',NF90_OUTTYPE_STAT,(/nbstat_dims, nfreq_dims, ndir_dims, ntime_dims /),var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, iret)

        iret=nf90_put_att(ncid,var_id,UNITS,'unk')
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 12, iret)

        iret=nf90_put_att(ncid,var_id,FULLNAME,'spectral energy density')
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 13, iret)
      END IF
      IF (VAROUT_STATION%WK) THEN
        iret=nf90_def_var(ncid,'WK',NF90_OUTTYPE_STAT,(/nbstat_dims, nfreq_dims, ntime_dims /),var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 14, iret)

        iret=nf90_put_att(ncid,var_id,UNITS,'unk')
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 15, iret)

        iret=nf90_put_att(ncid,var_id,FULLNAME,'wave number by frequency')
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 16, iret)
      END IF
      IF (VAROUT_STATION%ACOUT_1D) THEN
        iret=nf90_def_var(ncid,'ACOUT_1D',NF90_OUTTYPE_STAT,(/nbstat_dims, nfreq_dims, three_dims, ntime_dims /),var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 17, iret)

        iret=nf90_put_att(ncid,var_id,UNITS,'unk')
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 18, iret)

        iret=nf90_put_att(ncid,var_id,FULLNAME,'1-dimensional spectrum')
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 19, iret)
      END IF
      IF (VAROUT_STATION%ACOUT_2D) THEN
        iret=nf90_def_var(ncid,'ACOUT_2D',NF90_OUTTYPE_STAT,(/nbstat_dims, nfreq_dims, ndir_dims, ntime_dims /),var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 20, iret)

        iret=nf90_put_att(ncid,var_id,UNITS,'unk')
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 21, iret)

        iret=nf90_put_att(ncid,var_id,FULLNAME,'2-dimensional spectrum')
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 22, iret)
      END IF
      DO I=1,OUTVARS_COMPLETE
        IF (VAROUT_STATION%LVAR(I)) THEN
          CALL NAMEVARIABLE(I, eStr, eStrFullName, eStrUnit)
          iret=nf90_def_var(ncid,TRIM(eStr),NF90_OUTTYPE_STAT,(/ nbstat_dims, ntime_dims /),var_id)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 23, iret)

          iret=nf90_put_att(ncid,var_id,UNITS,TRIM(eStrUnit))
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 24, iret)

          iret=nf90_put_att(ncid,var_id,FULLNAME,TRIM(eStrFullName))
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 25, iret)
        END IF
      END DO
      iret=nf90_close(ncid)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 26, iret)

      iret=nf90_open(TRIM(FILE_NAME), NF90_WRITE, ncid)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 27, iret)

      CALL WRITE_NETCDF_HEADERS_STAT_2(ncid, MULTIPLEOUT_STAT)
      iret=nf90_close(ncid)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 28, iret)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WRITE_NETCDF_HEADERS_STAT_1(ncid, nbTime, MULTIPLEOUT)
      USE DATAPOOL
      USE NETCDF
      implicit none
      integer, intent(in) :: ncid, nbTime, MULTIPLEOUT
      character (len = *), parameter :: CallFct="WRITE_NETCDF_HEADERS_STAT_1"
      character (len = *), parameter :: UNITS = "units"
      integer one_dims, two_dims, three_dims, fifteen_dims
      integer nfreq_dims, ndir_dims
      integer nbstat_dims
      integer iret, var_id
      integer ntime_dims
      iret = nf90_def_dim(ncid, 'one', 1, one_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, iret)
      iret = nf90_def_dim(ncid, 'two', 2, two_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, iret)
      iret = nf90_def_dim(ncid, 'three', 3, three_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, iret)
      iret = nf90_def_dim(ncid, 'fifteen', 15, fifteen_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, iret)
      iret = nf90_def_dim(ncid, 'nbstation', IOUTS, nbstat_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, iret)
      iret = nf90_def_dim(ncid, 'nfreq', MSC, nfreq_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, iret)
      iret = nf90_def_dim(ncid, 'ndir', MDC, ndir_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 7, iret)
      !
# ifdef MPI_PARALL_GRID
      iret=nf90_def_var(ncid,'nproc',NF90_INT,(/one_dims/),var_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 8, iret)
      iret=nf90_put_att(ncid,var_id,UNITS,'integer')
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 9, iret)
      iret=nf90_put_att(ncid,var_id,'description','number of processors')
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 10, iret)
      !
      iret=nf90_def_var(ncid,'MULTIPLEOUT',NF90_INT,(/one_dims/),var_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, iret)
      iret=nf90_put_att(ncid,var_id,UNITS,'integer')
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 12, iret)
      iret=nf90_put_att(ncid,var_id,'description','multiple status')
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 13, iret)
# endif
      !
      IF (PARAMWRITE_STAT) THEN
        CALL WRITE_PARAM_1(ncid, one_dims)
      ENDIF
      !
      CALL WRITE_NETCDF_TIME_HEADER(ncid, nbTime, ntime_dims)
      IF (LSPHE) THEN
        iret=nf90_def_var(ncid,"lon",NF90_RUNTYPE,(/ nbstat_dims/),var_id)
      ELSE
        iret=nf90_def_var(ncid,"x",NF90_RUNTYPE,(/ nbstat_dims/),var_id)
      END IF
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 23, iret)
      IF (LSPHE) THEN
        iret=nf90_put_att(ncid,var_id,UNITS,'degree')
      ELSE
        iret=nf90_put_att(ncid,var_id,UNITS,'meter')
      END IF
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 24, iret)
! lat
      IF (LSPHE) THEN
        iret=nf90_def_var(ncid,"lat",NF90_RUNTYPE,(/ nbstat_dims/),var_id)
      ELSE
        iret=nf90_def_var(ncid,"y",NF90_RUNTYPE,(/ nbstat_dims/),var_id)
      END IF
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 25, iret)
      IF (LSPHE) THEN
        iret=nf90_put_att(ncid,var_id,UNITS,'degree')
      ELSE
        iret=nf90_put_att(ncid,var_id,UNITS,'meter')
      END IF
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 26, iret)
! cutoff frequency
      iret=nf90_def_var(ncid,'cutoff',NF90_RUNTYPE,(/ nbstat_dims/),var_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 27, iret)
      iret=nf90_put_att(ncid,var_id,UNITS,'Hz')
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 28, iret)
! ismax value
      iret=nf90_def_var(ncid,'ismax',NF90_INT,(/ nbstat_dims/),var_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 29, iret)
      iret=nf90_put_att(ncid,var_id,UNITS,'integer')
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 30, iret)
      IF (MULTIPLEOUT.gt.0) THEN
! Ifound
        iret=nf90_def_var(ncid,'ifound',NF90_INT,(/ nbstat_dims/),var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 31, iret)
        iret=nf90_put_att(ncid,var_id,UNITS,'integer')
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 32, iret)
      END IF
! Isum
      iret=nf90_def_var(ncid, 'isum', NF90_INT,(/ nbstat_dims/),var_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 33, iret)
      iret=nf90_put_att(ncid,var_id,UNITS,'integer')
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 34, iret)
! SPSIG
      iret=nf90_def_var(ncid, 'spsig', NF90_RUNTYPE,(/ nfreq_dims/),var_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 35, iret)
      iret=nf90_put_att(ncid,var_id,UNITS,'integer')
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 36, iret)
! SPDIR
      iret=nf90_def_var(ncid, 'spdir', NF90_RUNTYPE,(/ ndir_dims/),var_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 31, iret)
      iret=nf90_put_att(ncid,var_id,UNITS,'integer')
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 37, iret)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WRITE_NETCDF_HEADERS_STAT_2(ncid, MULTIPLEOUT)
      USE DATAPOOL
      USE NETCDF
      implicit none
      integer, intent(in) :: ncid, MULTIPLEOUT
      integer :: eWriteInt(1)
      real(rkind) :: eWriteReal(1)
      integer var_id, iret
      integer I
# ifdef MPI_PARALL_GRID
      integer eInt(1)
# endif
      character (len = *), parameter :: CallFct="WRITE_NETCDF_HEADERS_STAT_2"
      !
# ifdef MPI_PARALL_GRID
      iret=nf90_inq_varid(ncid,'nproc',var_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, iret)
      eInt(1)=nproc
      iret=nf90_put_var(ncid,var_id,eInt,start = (/1/), count=(/1/))
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, iret)
      !
      iret=nf90_inq_varid(ncid,'MULTIPLEOUT',var_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, iret)
      eInt(1)=MULTIPLEOUT
      iret=nf90_put_var(ncid,var_id,eInt,start = (/1/), count=(/1/))
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, iret)
# endif
      !
      IF (PARAMWRITE_STAT) THEN
        CALL WRITE_PARAM_2(ncid)
      ENDIF
      DO I=1,IOUTS 
        eWriteReal(1)=STATION(I) % XCOORD
        IF (LSPHE) THEN
          iret=nf90_inq_varid(ncid, "lon", var_id)
        ELSE
          iret=nf90_inq_varid(ncid, "x", var_id)
        END IF
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, iret)
        iret=nf90_put_var(ncid,var_id,eWriteReal, start=(/I/), count =(/1/) )
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, iret)
        !
        eWriteReal(1)=STATION(I) % YCOORD
        IF (LSPHE) THEN
          iret=nf90_inq_varid(ncid, "lat", var_id)
        ELSE
          iret=nf90_inq_varid(ncid, "y", var_id)
        END IF
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, iret)
        iret=nf90_put_var(ncid,var_id,eWriteReal, start=(/I/), count =(/1/) )
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, iret)
        !
        eWriteReal(1)=STATION(I) % CUTOFF
        iret=nf90_inq_varid(ncid, "cutoff", var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, iret)
        iret=nf90_put_var(ncid,var_id,eWriteReal, start=(/I/), count =(/1/) )
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, iret)
        !
        eWriteInt(1)=STATION(I) % ISMAX
        iret=nf90_inq_varid(ncid, "ismax", var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 7, iret)
        iret=nf90_put_var(ncid,var_id,eWriteInt, start=(/I/), count =(/1/) )
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 8, iret)
# ifdef MPI_PARALL_GRID
        IF (MULTIPLEOUT.eq.0) THEN
          eWriteInt(1)=STATION(I) % ISUM
          iret=nf90_inq_varid(ncid, "isum", var_id)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 9, iret)
          iret=nf90_put_var(ncid,var_id,eWriteInt, start=(/I/), count=(/1/) )
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 10, iret)
        ELSE
          eWriteInt(1)=STATION(I) % IFOUND
          iret=nf90_inq_varid(ncid, "ifound", var_id)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, iret)
          iret=nf90_put_var(ncid,var_id,eWriteInt, start=(/I/), count=(/1/) )
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 12, iret)
          !
          eWriteInt(1)=STATION(I) % ISUM
          iret=nf90_inq_varid(ncid, "isum", var_id)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 13, iret)
          iret=nf90_put_var(ncid,var_id,eWriteInt, start=(/I/), count=(/1/) )
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 14, iret)
        ENDIF
# else
        eWriteInt(1)=STATION(I) % IFOUND
        iret=nf90_inq_varid(ncid, "isum", var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 15, iret)
        iret=nf90_put_var(ncid,var_id,eWriteInt, start=(/I/), count =(/1/) )
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 16, iret)
# endif
        iret=nf90_inq_varid(ncid, "spsig", var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 17, iret)
        iret=nf90_put_var(ncid,var_id,SPSIG, start=(/1/), count =(/MSC/) )
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 18, iret)
        !
        iret=nf90_inq_varid(ncid, "spdir", var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 19, iret)
        iret=nf90_put_var(ncid,var_id,SPDIR, start=(/1/), count =(/MDC/) )
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 20, iret)
      ENDDO
      END SUBROUTINE
#endif
