!$Id: gotm_lib_version.F90,v 1.5 2005-11-15 11:35:02 lars Exp $
#include"cppdefs.h"
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Printing GOTM library version
!
! !INTERFACE:
      subroutine gotm_lib_version(unit)
!
! !DESCRIPTION:
!  Simply prints the version number of the GOTM turbulence library to unit.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: unit
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: gotm_lib_version.F90,v $
!  Revision 1.5  2005-11-15 11:35:02  lars
!  documentation finish for print
!
!  Revision 1.4  2005/06/27 13:44:07  kbk
!  modified + removed traling blanks
!
!  Revision 1.3  2003/03/28 09:20:35  kbk
!  added new copyright to files
!
!  Revision 1.2  2003/03/10 09:02:05  gotm
!  Added new Generic Turbulence Model + 
!  improved documentation and cleaned up code
!
!  Revision 1.1.1.1  2001/02/12 15:55:58  gotm
!  initial import into CVS
!
!EOP
!-------------------------------------------------------------------------
!BOC
   write(unit,*) 'GOTM library version: ',RELEASE

   return
   end
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
