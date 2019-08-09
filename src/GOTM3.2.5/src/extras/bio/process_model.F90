!$Id: process_model.F90,v 1.6 2004-07-30 09:22:20 hb Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Initialise the bio module
!
! !INTERFACE:
   subroutine process_model(first,numc,nlev,cc,pp,dd,h,t)
!
! !DESCRIPTION:
!
! !USES:
   use bio_var, only: bio_model
   use bio_template, only: do_bio_template
   use bio_npzd, only: do_bio_npzd
   use bio_iow, only: do_bio_iow
   use bio_sed, only: do_bio_sed
   use bio_fasham, only: do_bio_fasham
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   logical, intent(in)                 :: first
   integer, intent(in)                 :: numc,nlev
   REALTYPE, intent(in)                :: cc(1:numc,0:nlev)
   REALTYPE, intent(in)                :: h(0:nlev),t(0:nlev)
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE, intent(inout)             :: pp(1:numc,1:numc,0:nlev)
   REALTYPE, intent(inout)             :: dd(1:numc,1:numc,0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
   select case (bio_model)
      case (-1)
         call do_bio_template(numc,nlev)
      case (1)
         call do_bio_npzd(first,numc,nlev,cc,pp,dd)
      case (2)
         call do_bio_iow(first,numc,nlev,cc,pp,dd,h,t)
      case (3)
         call do_bio_sed()
      case (4)
         call do_bio_fasham(first,numc,nlev,cc,pp,dd)
      case default
         stop "bio: no valid biomodel specified in bio.inp !"
   end select
   return

   end subroutine process_model
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!----------------------------------------------------------------------- 
