MODULE SplitInterpo
      USE DATAPOOL
!**********************************************************************
!*  Design goal is for fast interpolation                             *
!*  For that we need to do the interpolation on blocks                *
!*  ---One level of blocks is the one of the parallelization          *
!*  ---Another level is by blocks of points (the quad)                *
!*  ---Thus we use the INE, MNP, MNE for everything                   *
!**********************************************************************
      TYPE QUADCOORD
        real(rkind) :: MinLon
        real(rkind) :: MaxLon
        real(rkind) :: MinLat
        real(rkind) :: MaxLat
      END TYPE QUADCOORD
      TYPE FULL_SPLIT_INFO
        integer :: nbQuad
        TYPE(QUADCOORD), dimension(:), pointer :: ListQuadReal
        integer, dimension(:), pointer :: ListStartIDXrange
        integer, dimension(:), pointer :: ListIE
      END TYPE FULL_SPLIT_INFO
      CONTAINS
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE IS_CONTAINED_IN_QUAD(eQuad, eLon, eLat, eAns)
      IMPLICIT NONE
      TYPE(QUADCOORD), intent(in) :: eQuad
      REAL(rkind), intent(in) :: eLon, elat
      LOGICAL, intent(out) :: eAns
      eAns=.TRUE.
      IF (eLon .lt. eQuad % MinLon) THEN
        eAns=.FALSE.
      END IF
      IF (eLon .gt. eQuad % MaxLon) THEN
        eAns=.FALSE.
      END IF
      IF (eLat .lt. eQuad % MinLat) THEN
        eAns=.FALSE.
      END IF
      IF (eLat .gt. eQuad % MaxLat) THEN
        eAns=.FALSE.
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_DECOMPOSITION(eFull, MaxSizeBlock)
      IMPLICIT NONE
      integer, intent(in) :: MaxSizeBlock
      TYPE(FULL_SPLIT_INFO), intent(inout) :: eFull
      integer NbQuadAlloc, NbQuadReal
      real(rkind), allocatable :: ListCentTrig(:,:)
      integer, allocatable :: ListStartIDX(:)
      integer, allocatable :: ListEndIDX(:)
      integer, allocatable :: ListIE(:)
      TYPE(QUADCOORD), allocatable :: ListQuadFormal(:)
      TYPE(QUADCOORD) :: Quad1, Quad2, Quad3, Quad4
      integer I, iPt, IE
      real(rkind) :: MinLon, MinLat, MaxLon, MaxLat, midLon, midLat
      real(rkind) :: eLon, eLat
      LOGICAL IsFirst
      LOGICAL eAns1, eAns2, eAns3, eAns4
      REAL(rkind) diff12, diff34
      integer idx, idxIE
      integer IP, iQuad, iQuadFind, lenFind, len
      integer iStart, iEnd, nb1, nb2, nb3, nb4
      CALL INITIAL_ALLOCATION()
      DO
        iQuadFind=-1
        lenFind=0
        DO iQuad=1,NbQuadReal
          len=ListEndIDX(iQuad) + 1 - ListStartIDX(iQuad)
          IF (len .gt. lenFind) THEN
            lenFind=len
            iQuadFind=iQuad
          END IF
        END DO
        IF (lenFind .lt. MaxSizeBlock) EXIT
        !
        iStart=ListStartIDX(iQuad)
        iEnd=ListEndIDX(iQuad)
        MinLon=ListQuadFormal(iQuadFind) % MinLon
        MinLat=ListQuadFormal(iQuadFind) % MinLat
        MaxLon=ListQuadFormal(iQuadFind) % MaxLon
        MaxLat=ListQuadFormal(iQuadFind) % MaxLat
        midLon = (MinLon + MaxLon)/2.0_rkind
        midLat = (MinLat + MaxLat)/2.0_rkind
        !
        Quad1 % MinLon = MinLon
        Quad1 % MaxLon = midLon !
        Quad1 % MinLat = MinLat
        Quad1 % MaxLat = MaxLat
        !
        Quad2 % MinLon = midLon !
        Quad2 % MaxLon = MaxLon
        Quad2 % MinLat = MinLat
        Quad2 % MaxLat = MaxLat
        !
        Quad3 % MinLon = MinLon
        Quad3 % MaxLon = MaxLon
        Quad3 % MinLat = midLat !
        Quad3 % MaxLat = MaxLat
        !
        Quad4 % MinLon = MinLon
        Quad4 % MaxLon = MaxLon
        Quad4 % MinLat = MinLat
        Quad4 % MaxLat = midLat !
        !
        nb1=0
        nb2=0
        nb3=0
        nb4=0
        DO idx=iStart,iEnd
          IE=ListIE(idx)
          eLon=ListCentTrig(1,IE)
          eLat=ListCentTrig(2,IE)
          CALL IS_CONTAINED_IN_QUAD(Quad1, eLon, eLat, eAns1)
          CALL IS_CONTAINED_IN_QUAD(Quad2, eLon, eLat, eAns2)
          CALL IS_CONTAINED_IN_QUAD(Quad3, eLon, eLat, eAns3)
          CALL IS_CONTAINED_IN_QUAD(Quad4, eLon, eLat, eAns4)
          IF (eAns1) nb1=nb1+1
          IF (eAns2) nb2=nb2+1
          IF (eAns3) nb3=nb3+1
          IF (eAns4) nb4=nb4+1
        END DO
        diff12=abs(nb1 - nb2)
        diff34=abs(nb3 - nb4)
        IF (diff12 .lt. diff34) THEN
          CALL SPLIT_ENTRY(iQuadFind, Quad1, Quad2)
        ELSE
          CALL SPLIT_ENTRY(iQuadFind, Quad3, Quad4)
        END IF
      END DO
      eFull % nbQuad = NbQuadReal
      allocate(eFull % ListStartIDXrange(NbQuadReal+1), eFull % ListQuadReal(NbQuadReal), eFull % ListIE(MNE), stat=istat)
      if (istat /= 0) CALL WWM_ABORT('allocate error 1')
      eFull % ListStartIDXrange(1)=1
      idxIE=0
      DO iQuad=1,NbQuadReal
        len=ListEndIDX(iQuad) + 1 - ListStartIDX(iQuad)
        eFull % ListStartIDXrange(iQuad+1) = eFull % ListStartIDXrange(iQuad) + len
        iStart=ListStartIDX(iQuad)
        IsFirst=.TRUE.
        DO iPt=1,len
          idxIE = idxIE + 1
          IE = ListIE(iStart + iPt - 1)
          eFull % ListIE(idxIE) = IE
          DO I=1,3
            IP=INE(I,IE)
            eLon=XP(IP)
            eLat=XP(IP)
            IF (IsFirst) THEN
              MinLON = eLon
              MinLAT = eLat
              MaxLON = eLon
              MaxLAT = eLat
              IsFirst=.FALSE.
            ELSE
              IF (eLon .gt. MaxLON) MaxLon=eLon
              IF (eLon .lt. MinLON) MinLon=eLon
              IF (eLat .gt. MaxLAT) MaxLat=eLat
              IF (eLat .lt. MinLAT) MinLat=eLat
            END IF
          END DO
        END DO
        eFull % ListQuadReal(iQuad) % MinLon = MinLon
        eFull % ListQuadReal(iQuad) % MinLat = MinLat
        eFull % ListQuadReal(iQuad) % MaxLon = MaxLon
        eFull % ListQuadReal(iQuad) % MaxLat = MaxLat
      END DO
      deallocate(ListStartIDX, ListEndIDX, ListIE)
      CONTAINS
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INITIAL_ALLOCATION()
      IMPLICIT NONE
      LOGICAL IsFirst
      integer IE, I
      real(rkind) :: sumLon, sumLAT, eLon, eLat, CentLon, CentLat
      allocate(ListIE(MNE), ListCentTrig(2,MNE), stat=istat)
      if (istat /= 0) CALL WWM_ABORT('allocate error 1')
      IsFirst=.TRUE.
      DO IE=1,MNE
        sumLON=0
        sumLAT=0
        DO I=1,3
          IP=INE(I,IE)
          eLon=XP(IP)
          eLat=YP(IP)
          IF (IsFirst) THEN
            MinLON=eLon
            MinLAT=eLat
            MaxLON=eLon
            MaxLAT=eLat
            IsFirst=.FALSE.
          ELSE
            IF (eLon .gt. MaxLON) MaxLon=eLon
            IF (eLon .lt. MinLON) MinLon=eLon
            IF (eLat .gt. MaxLAT) MaxLat=eLat
            IF (eLat .lt. MinLAT) MinLat=eLat
          END IF
          sumLON = sumLON + eLon
          sumLAT = sumLAT + eLat
        END DO
        CentLon = sumLON / 3.0_rkind
        CentLat = sumLAT / 3.0_rkind
        ListIE(IE)=IE
        ListCentTrig(1,IE)=CentLon
        ListCentTrig(2,IE)=CentLat
      END DO
      NbQuadReal=1
      NbQuadAlloc=1
      allocate(ListQuadFormal(NbQuadAlloc), ListStartIDX(NbQuadAlloc), ListEndIDX(NbQuadAlloc), stat=istat)
      if (istat /= 0) CALL WWM_ABORT('allocate error 1')
      ListStartIDX(1)=1
      ListEndIDX  (1)=MNE
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SPLIT_ENTRY(idxSplit, Quad1, Quad2)
      IMPLICIT NONE
      integer, intent(in) :: idxSplit
      TYPE(QUADCOORD), intent(in) :: Quad1, Quad2
      !
      TYPE(QUADCOORD), allocatable :: ListQuadFormal_copy(:)
      integer, allocatable :: ListStartIDX_copy(:), ListEndIDX_copy(:)
      integer iQuad, nb1, nb2, idx1, idx2
      integer, allocatable :: ListIEcopy(:)
      logical eAns1, eAns2
      integer nbError
      IF (NbQuadReal .eq. NbQuadAlloc) THEN
        allocate(ListQuadFormal_copy(NbQuadReal), ListStartIDX_copy(NbQuadReal), ListEndIDX_copy(NbQuadReal), stat=istat)
        if (istat /= 0) CALL WWM_ABORT('allocate error 1')
        DO iQuad=1,NbQuadReal
          ListQuadFormal_copy(iQuad)=ListQuadFormal(iQuad)
          ListStartIDX_copy(iQuad)=ListStartIDX(iQuad)
          ListEndIDX_copy(iQuad)=ListEndIDX(iQuad)
        END DO
        deallocate(ListQuadFormal, ListStartIDX, ListEndIDX)
        NbQuadAlloc=2*NbQuadAlloc
        allocate(ListQuadFormal(NbQuadAlloc), ListStartIDX(NbQuadAlloc), ListEndIDX(NbQuadAlloc), stat=istat)
        if (istat /= 0) CALL WWM_ABORT('allocate error 1')
        DO iQuad=1,NbQuadReal
          ListQuadFormal(iQuad)=ListQuadFormal_copy(iQuad)
          ListStartIDX(iQuad)=ListStartIDX_copy(iQuad)
          ListEndIDX(iQuad)=ListEndIDX_copy(iQuad)
        END DO
        deallocate(ListQuadFormal_copy, ListStartIDX_copy, ListEndIDX_copy)
      END IF
      NbQuadReal=NbQuadReal+1
      ListQuadFormal(idxSplit)=Quad1
      ListQuadFormal(NbQuadReal)=Quad2
      iStart=ListStartIDX(idxSplit)
      iEnd=ListStartIDX(idxSplit)
      nb1=0
      nb2=0
      nbError=0
      allocate(ListIEcopy(iStart:iEnd), stat=istat)
      if (istat /= 0) CALL WWM_ABORT('allocate error 1')
      len=iEnd+1-iStart
      DO idx=iStart,iEnd
        IE=ListIE(idx)
        ListIEcopy(idx)=IE
        eLon=ListCentTrig(1,IE)
        eLat=ListCentTrig(2,IE)
        CALL IS_CONTAINED_IN_QUAD(Quad1, eLon, eLat, eAns1)
        CALL IS_CONTAINED_IN_QUAD(Quad2, eLon, eLat, eAns2)
        IF (eAns1) nb1=nb1+1
        IF (eAns2) nb2=nb2+1
        IF (eAns1 .and. eAns2) nbError=nbError+1
        IF (.NOT.(eAns1) .and. .NOT.(eAns2)) nbError=nbError+1
      END DO
      IF (nbError .gt. 0) CALL WWM_ABORT('Error in group construction')
      ListStartIDX(idxSplit) = iStart
      ListEndIDX(idxSplit) = iStart+nb1-1
      ListStartIDX(NbQuadReal) = iStart+nb1
      ListEndIDX(NbQuadReal) = iEnd
      idx1=0
      idx2=0
      DO idx=iStart,iEnd
        IE=ListIEcopy(idx)
        CALL IS_CONTAINED_IN_QUAD(Quad1, eLon, eLat, eAns1)
        CALL IS_CONTAINED_IN_QUAD(Quad2, eLon, eLat, eAns2)
        IF (eAns1) THEN
          ListIE(iStart + idx1)=IE
          idx1 = idx1 + 1
        END IF
        IF (eAns2) THEN
          ListIE(iStart + nb1 + idx2)=IE
          idx2 = idx2+1
        END IF
      END DO
      deallocate(ListIEcopy)
      END SUBROUTINE
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_CONTAINING_TRIANGLE(eFull, eLon, eLat, IEcont)
      IMPLICIT NONE
      TYPE(FULL_SPLIT_INFO), intent(in) :: eFull
      REAL(rkind), intent(in) :: eLon, elat
      integer, intent(out) :: IEcont
      integer I, IP, idx, iQuad, IDXstart, IDXend, IE
      real(rkind) :: X(3), Y(3), WI(3)
      LOGICAL eAns
      DO iQuad=1,eFull % nbQuad
        CALL IS_CONTAINED_IN_QUAD(eFull % ListQuadReal(iQuad), eLon, eLat, eAns)
        IF (eAns) THEN
          IDXstart=eFull % ListStartIDXrange(iQuad)
          IDXend  =eFull % ListStartIDXrange(iQuad+1)-1
          DO idx=IDXstart,IDXend
            IE=eFull % ListIE(idx)
            DO I=1,3
              IP=INE(I,IE)
              X(I)=XP(IP)
              Y(I)=YP(IP)
            END DO
            CALL INTELEMENT_COEF(X,Y,eLon,eLat,WI)
            IF (minval(WI) .ge. -THR) THEN
              IEcont=IE
              RETURN
            END IF
          END DO
        END IF
      END DO
      IEcont=-1
      END SUBROUTINE
END MODULE
