!$Id: bio_template.F90,v 1.2 2004-07-30 09:22:20 hb Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bio_template --- TEMPLATE bio model \label{sec:bio_template}
!
! !INTERFACE:
   module bio_template
!
! !DESCRIPTION:
!  Remember this Hans
!
! !USES:
!  default: all is private.
   use bio_var
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_bio_template, init_var_template, var_info_template, &
          light_template, do_bio_template, end_bio_template
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: bio_template.F90,v $
!  Revision 1.2  2004-07-30 09:22:20  hb
!  use bio_var in specific bio models - simpliefied internal interface
!
!  Revision 1.1  2003/07/23 12:27:31  hb
!  more generic support for different bio models
!
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the template bio module
!
! !INTERFACE:
   subroutine init_bio_template(namlst,fname,unit)
!
! !DESCRIPTION:
!  Here, the bio namelist {\tt bio_template.inp} is read and memory is
!  allocated - and various variables are initialised.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,          intent(in)   :: namlst
   character(len=*), intent(in)   :: fname
   integer,          intent(in)   :: unit
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL2 'init_bio_template'

   numc=5

   LEVEL3 'TEMPLATE bio module initialised ...'

   return

   end subroutine init_bio_template
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the concentration variables
!
! !INTERFACE:
   subroutine init_var_template(nlev)
!
! !DESCRIPTION:
!  Here, the cc and ws varibles are filled with initial conditions
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: nlev
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding

! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC

   LEVEL3 'TEMPLATE bio model should be initialised here ...'

   return
   end subroutine init_var_template
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Providing info on variables
!
! !INTERFACE:
   subroutine var_info_template()
!
! !DESCRIPTION:
!  This subroutine provides information on the variables. To be used
!  when storing data in NetCDF files.
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   var_names(1) = 'name_1'
   var_units(1) = 'units_1'
   var_long(1)  = 'long_1'

   var_names(2) = 'name_2'
   var_units(2) = 'units_2'
   var_long(2)  = 'long_2'

   var_names(3) = 'name_3'
   var_units(3) = 'units_3'
   var_long(3)  = 'long_3'

   var_names(4) = 'name_4'
   var_units(4) = 'units_4'
   var_long(4)  = 'long_4'

   var_names(5) = 'name_5'
   var_units(5) = 'units_5'
   var_long(5)  = 'long_5'

   return
   end subroutine var_info_template
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Light properties for the template model
!
! !INTERFACE
   subroutine light_template(nlev,h,rad,bioshade_feedback,bioshade)
!
! !DESCRIPTION
!
! !USES
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer                              :: nlev
   REALTYPE, intent(in)                 :: h(0:nlev)
   REALTYPE, intent(in)                 :: rad(0:nlev)
   logical, intent(in)                  :: bioshade_feedback
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: bioshade(0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC

   LEVEL2 'calculating TEMPLATE light properties'

   end subroutine light_template
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of geobiochemical model
!
! !INTERFACE
   subroutine do_bio_template(numc,nlev)
!
! !DESCRIPTION
!
! !USES
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
  integer                              :: numc,nlev
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC

   STDERR 'Here the biological process model is defined.'
   STDERR 'Since this is only a template the program is now terminated.'
   stop 'do_bio_template()'
   return
   end subroutine do_bio_template
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finish the bio calculations
!
! !INTERFACE:
   subroutine end_bio_template
!
! !DESCRIPTION:
!  Nothing done yet --- supplied for completeness.
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC

   return
   end subroutine end_bio_template
!EOC

!-----------------------------------------------------------------------

   end module bio_template

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
