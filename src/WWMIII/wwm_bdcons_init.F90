#include "wwm_functions.h"
#define DEBUG
#undef DEBUG
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_IOBPD
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER :: I
      REAL(rkind) :: DXP1, DXP2, DXP3, DYP1, DYP2, DYP3
      REAL(rkind) :: x1, y1, x2, y2
      INTEGER :: I1, I2, I3, IE, IP, ID
#ifdef MPI_PARALL_GRID
      INTEGER :: iwild(mnp)
#endif
      REAL(rkind) :: EVX, EVY
      REAL(rkind) :: eDet1, eDet2
      IOBPD = 0

      DO IE=1,MNE
        I1   =   INE(1,IE)
        I2   =   INE(2,IE)
        I3   =   INE(3,IE)
        DXP1 =   IEN(6,IE)
        DYP1 = - IEN(5,IE)
        DXP2 =   IEN(2,IE)
        DYP2 = - IEN(1,IE)
        DXP3 =   IEN(4,IE)
        DYP3 = - IEN(3,IE)
!AR: ... modifly wave direction by currents ...
        DO ID=1,MDC
          EVX=COSTH(ID)
          EVY=SINTH(ID)
          DO I=1,3
            IF (I.eq.1) THEN
              x1=   DXP1
              y1=   DYP1
              x2= - DXP3
              y2= - DYP3
              IP=   I1
            END IF
            IF (I.eq.2) THEN
              x1 =   DXP2
              y1 =   DYP2
              x2 = - DXP1
              y2 = - DYP1
              IP =   I2
            END IF
            IF (I.eq.3) THEN
              x1 =   DXP3
              y1 =   DYP3
              x2 = - DXP2
              y2 = - DYP2
              IP =   I3
            END IF
!AR: MDS please check if the new thr can pose a problem ...
            IF (abs(iobp(ip)) .eq. 1) THEN
              eDet1 = THR-x1*EVY+y1*EVX
              eDet2 = THR+x2*EVY-y2*EVX
              IF ((eDet1.gt.ZERO).and.(eDet2.gt.ZERO)) THEN
                IOBPD(ID,IP)=1
              ENDIF 
            ELSE ! land boundary ...
              IOBPD(ID,IP)=1
            END IF
          END DO
        END DO
      END DO
      DO IP = 1, MNP
        IF ( (LBCWA .OR. LBCSP) ) THEN
          IF ( IOBP(IP) == 2 .OR. IOBP(IP) == 4) THEN
            IOBWB(IP) = 0
            IOBPD(:,IP) = 1
          ENDIF
        END IF
        IF ( IOBP(IP) == 3 .OR. IOBP(IP) == 4) THEN ! If Neumann boundary condition is given set IOBP to 3
          IOBPD(:,IP) = 1 ! Update Neumann nodes ...
        END IF
      END DO
!2do: recode for mpi 
!        IF (LBCWA .OR. LBCSP) THEN
!          IF (.NOT. ANY(IOBP .EQ. 2)) THEN
!            CALL WWM_ABORT('YOU IMPOSED BOUNDARY CONDITIONS BUT IN THE BOUNDARY FILE ARE NO NODES WITH FLAG = 2')
!          ENDIF
!        ENDIF
#ifdef MPI_PARALL_GRID
      CALL exchange_p2di(IOBWB)
      DO ID = 1, MDC
        iwild = IOBPD(ID,:)
        CALL exchange_p2di(iwild)
        IOBPD(ID,:) = iwild
      ENDDO
#endif

#if defined DEBUG && defined IOBPDOUT
# ifdef MPI_PARALL_GRID
      IF (myrank == 0) THEN
# endif
        DO IP = 1, MNP
          WRITE(IOBPOUT%FHNDL,*) IP, ID, IOBWB(IP)
        END DO
        DO IP = 1, MNP
          DO ID = 1, MDC
            WRITE(IOBPDOUT%FHNDL,*) IP, ID, SPDIR(ID)*RADDEG, IOBPD(ID,IP), IOBP(IP)
          ENDDO
        END DO
# ifdef MPI_PARALL_GRID
      END IF
# endif
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_IOBPD_BY_DEP
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER              :: IP
      DO IP = 1, MNP
        IF (DEP(IP) .LT. DMIN) THEN 
          IOBDP(IP) = 0 
        ELSE 
          IOBDP(IP) = 1
        ENDIF
      END DO
      END SUBROUTINE
!**********************************************************************
!* Determine status of node points by purely combinatorial and fast   *
!* return values:                                                     *
!* ---  0: should not happen in return                                *
!* ---  1: interior point of the grid.                                *
!* --- -1: boundary node                                              *
!**********************************************************************
      SUBROUTINE GET_BOUNDARY_STATUS(STATUS)
      USE DATAPOOL
      implicit none
      integer, intent(out) :: STATUS(MNP)
      INTEGER :: COLLECTED(MNP), NEXTVERT(MNP), PREVVERT(MNP)
      INTEGER          :: ISFINISHED, INEXT, IPREV
      INTEGER          :: IPNEXT, IPPREV, ZNEXT, IP, I, IE
      integer nb0, nb1, nbM1
      STATUS(:) = 0
      DO IE=1,MNE
        DO I=1,3
          IF (I.EQ.1) THEN
            IPREV=3
          ELSE
            IPREV=I-1
          END IF
          IF (I.EQ.3) THEN
            INEXT=1
          ELSE
            INEXT=I+1
          END IF
          IP=INE(I,IE)
          IPNEXT=INE(INEXT,IE)
          IPPREV=INE(IPREV,IE)
          IF (STATUS(IP).EQ.0) THEN
            STATUS(IP)=1
            PREVVERT(IP)=IPPREV
            NEXTVERT(IP)=IPNEXT
          END IF
        END DO
      END DO
      STATUS(:)=0
      DO
        COLLECTED(:)=0
        DO IE=1,MNE
          DO I=1,3
            IF (I.EQ.1) THEN
              IPREV=3
            ELSE
              IPREV=I-1
            END IF
            IF (I.EQ.3) THEN
              INEXT=1
            ELSE
              INEXT=I+1
            END IF
            IP=INE(I,IE)
            IPNEXT=INE(INEXT,IE)
            IPPREV=INE(IPREV,IE)
            IF (STATUS(IP).eq.0) THEN
              ZNEXT=NEXTVERT(IP)
              IF (ZNEXT.eq.IPPREV) THEN
                COLLECTED(IP)=1
                NEXTVERT(IP)=IPNEXT
                IF (NEXTVERT(IP).eq.PREVVERT(IP)) THEN
                  STATUS(IP)=1
                END IF
              END IF
            END IF
          END DO
        END DO
        ISFINISHED=1
        DO IP=1,MNP
          IF ((COLLECTED(IP).eq.0).and.(STATUS(IP).eq.0)) THEN
            STATUS(IP)=-1
          END IF
          IF (STATUS(IP).eq.0) THEN
            ISFINISHED=0
          END IF
        END DO
        IF (ISFINISHED.eq.1) THEN
          EXIT
        END IF
      END DO
#ifdef MPI_PARALL_GRID
      CALL exchange_p2di(STATUS)
#endif
      nb0=0
      nb1=0
      nbM1=0
      DO IP=1,MNP
        IF (STATUS(IP) .eq. 0) THEN
          nb0=nb0+1
        END IF
        IF (STATUS(IP) .eq. 1) THEN
          nb1=nb1+1
        END IF
        IF (STATUS(IP) .eq. -1) THEN
          nbM1=nbM1+1
        END IF
      END DO
      WRITE(STAT%FHNDL,*) 'Number of  0 in STATUS=', nb0
      WRITE(STAT%FHNDL,*) 'Number of  1 in STATUS=', nb1
      WRITE(STAT%FHNDL,*) 'Number of -1 in STATUS=', nbM1
      FLUSH(STAT%FHNDL)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SINGLE_READ_IOBP_TOTAL(IOBPtotal, IGRIDTYPE, eBND, np_total)
      USE DATAPOOL, only : rkind, FILEDEF, istat, STAT, wwmerr
#ifdef NCDF
      USE NETCDF
#endif
      IMPLICIT NONE
      INTEGER, intent(out) :: IOBPtotal(np_total)
      INTEGER, intent(in) :: IGRIDTYPE
      type(FILEDEF), intent(in) :: eBND
      integer, intent(in) :: np_total
      !
      INTEGER I, IP, IFSTAT
      REAL(rkind)       :: ATMP, BTMP, BNDTMP
      INTEGER ITMP
      integer nb0, nb1, nb2, nb3, nb4
#ifdef NCDF
      INTEGER ncid, var_id
      character (len = *), parameter :: CallFct="SINGLE_READ_IOBP_TOTAL"
#endif
      CALL TEST_FILE_EXIST_DIE('Missing boundary file : ', TRIM(eBND%FNAME))
      IOBPtotal = 0
!
! Reading of raw boundary file
!
      IF (IGRIDTYPE.eq.1) THEN ! XFN 
        OPEN(eBND%FHNDL, FILE = eBND%FNAME, STATUS = 'OLD')
        DO I = 1, 2
          READ(eBND%FHNDL,*)
        END DO
        READ(eBND%FHNDL,*)
        READ(eBND%FHNDL,*)
        READ(eBND%FHNDL,*)
        DO I = 1, 7
          READ(eBND%FHNDL,*)
        END DO
        DO IP = 1, NP_TOTAL
          READ(eBND%FHNDL, *, IOSTAT = IFSTAT) ITMP, BNDTMP, BNDTMP, BNDTMP
          IF ( IFSTAT /= 0 ) THEN
            CALL WWM_ABORT('error in the bnd file 1')
          END IF
          ITMP=INT(BNDTMP)
          IOBPtotal(IP) = ITMP
        END DO
        CLOSE(eBND%FHNDL)
      ELSE IF (IGRIDTYPE.eq.2) THEN ! Periodic  
        OPEN(eBND%FHNDL, FILE = eBND%FNAME, STATUS = 'OLD')
        READ(eBND%FHNDL,*)
        READ(eBND%FHNDL,*)
        DO IP = 1, NP_TOTAL
          READ(eBND%FHNDL, *, IOSTAT = IFSTAT) ITMP, BNDTMP, BNDTMP, BNDTMP
          IF ( IFSTAT /= 0 ) THEN
            CALL WWM_ABORT('error in the bnd file 2')
          END IF
          ITMP=INT(BNDTMP)
          IOBPtotal(IP) = ITMP
        END DO
        CLOSE(eBND%FHNDL)
      ELSE IF (IGRIDTYPE.eq.3) THEN ! SCHISM 
        OPEN(eBND%FHNDL, FILE = eBND%FNAME, STATUS = 'OLD')
        READ(eBND%FHNDL,*)
        READ(eBND%FHNDL,*)
        DO IP = 1, NP_TOTAL
          READ(eBND%FHNDL, *, IOSTAT = IFSTAT) ITMP, BNDTMP, BNDTMP, BNDTMP
          IF ( IFSTAT /= 0 ) THEN
            CALL WWM_ABORT('error in the bnd file 3')
          END IF
          ITMP=INT(BNDTMP)
          IOBPtotal(IP) = ITMP
        END DO
        CLOSE(eBND%FHNDL)
      ELSE IF (IGRIDTYPE.eq.4) THEN ! WWMOLD 
        OPEN(eBND%FHNDL, FILE = eBND%FNAME, STATUS = 'OLD')
        READ(eBND%FHNDL,*)
        READ(eBND%FHNDL,*)
        DO IP = 1, NP_TOTAL
          READ(eBND%FHNDL, *, IOSTAT = IFSTAT) ATMP, BTMP, BNDTMP
          IF ( IFSTAT /= 0 ) THEN
            CALL WWM_ABORT('error in the bnd file 4')
          END IF
          ITMP=INT(BNDTMP)
          IOBPtotal(IP) = ITMP
        END DO
        CLOSE(eBND%FHNDL)
#ifdef NCDF
      ELSE IF (IGRIDTYPE.eq.5) THEN ! netcdf boundary format
        ISTAT = NF90_OPEN(eBND%FNAME, NF90_NOWRITE, ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, ISTAT)

        ISTAT = nf90_inq_varid(ncid, 'IOBP', var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, ISTAT)

        ISTAT = nf90_get_var(ncid, var_id, IOBPtotal)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, ISTAT)

        ISTAT = NF90_CLOSE(ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, ISTAT)
#endif
      END IF
      DO IP = 1, NP_TOTAL
        IF (IOBPtotal(IP) .GT. 4) THEN
          WRITE(wwmerr, *) 'NextGen: We need iobp<=2 but ip=', IP, ' iobp=', IOBPtotal(IP)
          CALL WWM_ABORT(wwmerr)
        ENDIF
      ENDDO
#if defined DEBUG && defined IOBPDOUT
# ifdef MPI_PARALL_GRID
      IF (myrank == 0) THEN
# endif
        write(16,*) 'WWM: BC,',IOBPOUT%FHNDL, IOBPOUT%FNAME
        DO IP = 1, NP_TOTAL
          WRITE(IOBPOUT%FHNDL,*) IP, IOBPtotal(IP), 'TEST'
        END DO
        FLUSH(IOBPOUT%FHNDL)
# ifdef MPI_PARALL_GRID
      ENDIF 
# endif
#endif
      nb0=0
      nb1=0
      nb2=0
      nb3=0
      nb4=0
      DO IP=1,NP_TOTAL
        IF (IOBPtotal(IP) .eq. 0) THEN
          nb0=nb0+1
        END IF
        IF (IOBPtotal(IP) .eq. 1) THEN
          nb1=nb1+1
        END IF
        IF (IOBPtotal(IP) .eq. 2) THEN
          nb2=nb2+1
        END IF
        IF (IOBPtotal(IP) .eq. 3) THEN
          nb3=nb3+1
        END IF
        IF (IOBPtotal(IP) .eq. 4) THEN
          nb4=nb4+1
        END IF
      END DO
      WRITE(STAT%FHNDL,*) 'Number of 0 in IOBPtotal=', nb0
      WRITE(STAT%FHNDL,*) 'Number of 1 in IOBPtotal=', nb1
      WRITE(STAT%FHNDL,*) 'Number of 2 in IOBPtotal=', nb2
      WRITE(STAT%FHNDL,*) 'Number of 3 in IOBPtotal=', nb3
      WRITE(STAT%FHNDL,*) 'Number of 4 in IOBPtotal=', nb4
      FLUSH(STAT%FHNDL)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_IOBP_TOTAL
      USE DATAPOOL
      IMPLICIT NONE
      integer iProc
      allocate(IOBPtotal(np_total), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 1')
#ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_BOUND) THEN
        CALL SINGLE_READ_IOBP_TOTAL(IOBPtotal, IGRIDTYPE, BND, np_total)
      ELSE
        IF (myrank .eq. 0) THEN
          CALL SINGLE_READ_IOBP_TOTAL(IOBPtotal, IGRIDTYPE, BND, np_total)
          DO iProc=2,nproc
            CALL MPI_SEND(IOBPtotal,np_total,itype, iProc-1, 30, comm, ierr)
          END DO
        ELSE
          CALL MPI_RECV(IOBPtotal,np_total,itype, 0, 30, comm, istatus, ierr)
        END IF
      END IF
#else
      CALL SINGLE_READ_IOBP_TOTAL(IOBPtotal, IGRIDTYPE, BND, np_total)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PRINT_STATISTICS_IOBP_TOTAL
      USE DATAPOOL
      IMPLICIT NONE
      integer :: ContElements(np_total)
      integer :: ListDegWork(np_total)
      integer, allocatable :: ListAdjWithDupl(:,:)
      integer, allocatable :: IEcontain(:,:)
      integer, allocatable :: StatusAdj(:)
      integer IE, I, IP, INEXT, IPREV, IPadj
      integer IP_N, IP_P, eDeg, nb1, nb2, nb, J
      integer NumberAllTwo, NumberBoundary, NumberPathological
      integer SatMaxDeg, pos, MaxIEcont, eDegVert
      ContElements=0
      DO IE=1,NE_TOTAL
        DO I=1,3
          IP=INEtotal(I,IE)
          ContElements(IP)=ContElements(IP)+1
        END DO
      END DO
      MaxIEcont=maxval(ContElements)
      SatMaxDeg=2*MaxIEcont
      allocate(ListAdjWithDupl(SatMaxDeg,NP_TOTAL), stat=istat)
      allocate(IEcontain(MaxIEcont,NP_TOTAL), stat=istat)
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
            WRITE(DBG%FHNDL, *) 'IE=', IE
            WRITE(DBG%FHNDL, *) 'I=', I, 'IP=', IP
            WRITE(DBG%FHNDL, *) 'INEXT=', INEXT, ' IP_N=', IP_N
            WRITE(DBG%FHNDL, *) 'IPREV=', IPREV, ' IP_P=', IP_P
            CALL WWM_ABORT("logical error")
          END IF
          IEcontain(pos+1,IP)=IE
          ListDegWork(IP)=pos+1
        END DO
      END DO
      allocate(StatusAdj(SatMaxDeg), stat=istat)
      NumberAllTwo=0
      NumberBoundary=0
      NumberPathological=0
      DO IP=1,NP_TOTAL
        eDeg=ListDegWork(IP)
!        Print *, 'IP=', IP, ' eDeg=', eDeg
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
!           Print *, '  nb=', nb
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
          END IF
        END DO
        IF (nb1 .eq. 0) NumberAllTwo=NumberAllTwo + 1
        IF (nb1 .eq. 1) CALL WWM_ABORT("Number 1 should not happen")
        IF (nb1 .eq. 2) NumberBoundary=NumberBoundary + 1
        IF (nb1 .gt. 2) NumberPathological=NumberPathological + 1
      END DO
      deallocate(StatusAdj, ListAdjWithDupl)
      WRITE(STAT%FHNDL,*) 'NumberAllTwo      =', NumberAllTwo
      WRITE(STAT%FHNDL,*) 'NumberBoundary    =', NumberBoundary
      WRITE(STAT%FHNDL,*) 'NumberPathological=', NumberPathological
      FLUSH(STAT%FHNDL)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_IOBP_NEXTGENERATION
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER     :: IP, IFSTAT, SPsize
      REAL(rkind) :: BNDTMP
      INTEGER     :: STATUS(MNP)
      INTEGER     :: idx, PosWBAC
      IOBPD   = 0
      CALL READ_IOBP_TOTAL
      CALL PRINT_STATISTICS_IOBP_TOTAL
      IOBP    = 0
      DO IP = 1, NP_TOTAL
#ifdef MPI_PARALL_GRID
# ifndef PDLIB
        IF (ipgl(ip)%rank == myrank) THEN
          IOBP(ipgl(ip)%id) = IOBPtotal(IP)
        END IF
# else
        IF (ipgl(ip)%id .gt. 0) THEN
          IOBP(ipgl(ip)%id) = IOBPtotal(IP)
        END IF
# endif
#else
        IOBP(IP) = IOBPtotal(IP)
#endif
      END DO
#ifdef SCHISM
      DO IP = 1, NP_RES ! reset boundary flag in the case that wave boundary are not used but defined in the boundary file
        IF (.NOT. LBCWA .AND. .NOT. LBCSP) THEN
          IF (IOBP(IP) .EQ. 2 .OR. IOBP(IP) .EQ. 4) IOBP(IP) = 1
        ENDIF
      ENDDO

      DO IP = 1, NP_RES
        IF (IOBP(IP) .ne. 2 .and. IOBP(IP) .ne. 3 .and. IOBP(IP) .ne. 4) THEN
          IF (abs(ibnd_ext_int(IP)) == 1) THEN
            IOBP(IP) = ibnd_ext_int(IP)
          ENDIF 
        END IF
      END DO
      CALL EXCHANGE_P2DI(IOBP)
#endif
!
! indexing boundary nodes ...
!
      ! Local  boundary nodes ...
      IWBMNP = 0
      DO IP = 1, MNP
        IF (IOBP(IP) == 2 .OR. IOBP(IP) == 4) IWBMNP = IWBMNP + 1
      END DO
      WRITE(STAT%FHNDL,*) 'IWBMNP=', IWBMNP
      FLUSH(STAT%FHNDL)
      ALLOCATE( IWBNDLC(IWBMNP), IWBNDLC_REV(MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 2')
      IWBNDLC_REV=0
      idx = 0
      DO IP = 1, MNP
        IF (IOBP(IP) == 2 .OR. IOBP(IP) == 4) THEN
          idx = idx + 1
          IWBNDLC(idx) = IP ! Stores local wave boundary index
          IF (LINHOM) THEN
            PosWBAC=idx
          ELSE
            PosWBAC=1
          END IF
          IWBNDLC_REV(IP) = PosWBAC
        END IF
      END DO
      ! Global boundary nodes
      IWBMNPGL = 0
      DO IP = 1, NP_TOTAL
        IF (IOBPtotal(IP) == 2 .OR. IOBPtotal(IP) == 4) IWBMNPGL = IWBMNPGL + 1
      END DO
      WRITE(STAT%FHNDL,*) 'IWBMNPGL=', IWBMNPGL
      FLUSH(STAT%FHNDL)
      ALLOCATE( IWBNDGL(IWBMNPGL), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 3')
      idx=0
      DO IP = 1, NP_TOTAL
        IF (IOBPtotal(IP) == 2 .OR. IOBPtotal(IP) == 4) THEN
          idx = idx + 1
          IWBNDGL(idx) = IP
        END IF
      END DO
!
! find islands and domain boundary ....
!
      CALL GET_BOUNDARY_STATUS(STATUS)
      DO IP=1,MNP
        IF (STATUS(IP).eq.-1 .AND. IOBP(IP) .EQ. 0) THEN
          IOBP(IP)=1
        END IF
      END DO
#ifdef MPI_PARALL_GRID
      CALL EXCHANGE_P2DI(IOBP)
#endif
      IF (HACK_HARD_SET_IOBP) THEN
        ! For hackish purposes, sometimes, we want to have an interior point
        ! set to Neumann condition 2 for example.
        ! Use that option only if you are sure of yourself.
        DO IP = 1, NP_TOTAL
#ifdef MPI_PARALL_GRID
# ifndef PDLIB
          IF (ipgl(ip)%rank == myrank) THEN
            IOBP(ipgl(ip)%id) = IOBPtotal(IP)
          END IF
# else
          IF (ipgl(ip)%id .gt. 0) THEN
            IOBP(ipgl(ip)%id) = IOBPtotal(IP)
          END IF
# endif
#else
          IOBP(IP) = IOBPtotal(IP)
#endif
        END DO
      END IF
!
! allocate wave boundary arrays ... 
!
      IF (LINHOM) THEN
        SPsize=IWBMNP
      ELSE
        SPsize=1
      ENDIF
      IF (LBCWA) THEN
        ALLOCATE( SPPARM(8,SPsize), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 4')
        SPPARM = 0.
      ENDIF
      IF (LBCWA .OR. LBCSP) THEN
        ALLOCATE( WBAC(MSC,MDC,SPsize), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 5')
        WBAC = 0.
        IF (LBINTER) THEN
          ALLOCATE( WBACOLD(MSC,MDC,SPsize), WBACNEW(MSC,MDC,SPsize), DSPEC(MSC,MDC,SPsize), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 6')
          WBACOLD = 0.
          WBACNEW = 0.
          DSPEC   = 0.
        ENDIF
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CLOSE_IOBP
      USE DATAPOOL
      IMPLICIT NONE
      DEALLOCATE( IWBNDLC)
#ifdef MPI_PARALL_GRID
      DEALLOCATE( IWBNDGL)
#endif
      IF (LBCWA .OR. LBCSP) THEN
        CLOSE(WAV%FHNDL)
      END IF
      IF (LBCWA) THEN
        DEALLOCATE( SPPARM)
      ENDIF
      IF (LBCWA .OR. LBCSP) THEN
        DEALLOCATE( WBAC)
        IF (LBINTER) DEALLOCATE( WBACOLD, WBACNEW, DSPEC)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
