!  Single precision function to calculate Julian date
      FUNCTION julian_date(YYYY,MM,DD)
!      IMPLICIT NONE
      INTEGER, INTENT(IN) :: YYYY,MM,DD
!              DATE ROUTINE julian_date(YYYY,MM,DD) CONVERTS CALENDER DATE TO
!              JULIAN DATE.  SEE CACM 1968 11(10):657, LETTER TO THE
!              EDITOR BY HENRY F. FLIEGEL AND THOMAS C. VAN FLANDERN.
!    EXAMPLE julian_date(1970,1,1)=2440588
      INTEGER :: julian_date

      julian_date=DD-32075+1461*(YYYY+4800+(MM-14)/12)/4 & !days
     &         +367*(MM-2-((MM-14)/12)*12)/12-3* &
     &         ((YYYY+4900+(MM-14)/12)/100)/4
      RETURN
      END

