!$Id: bio_template.F90,v 1.5 2007-01-06 11:49:15 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bio_template --- template bio model \label{sec:biotemplate}
!
! !INTERFACE:
   module bio_template
!
! !DESCRIPTION:
! This is a template for including new biogeochemical models into
! GOTM. It has the full structural functionality of a GOTM
! biogeochemical model, but is terminated  in {\tt do\_bio\-template},
! where the right hand sides should be calculated.
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
!  Revision 1.5  2007-01-06 11:49:15  kbk
!  namelist file extension changed .inp --> .nml
!
!  Revision 1.4  2006-10-26 13:12:46  kbk
!  updated bio models to new ode_solver
!
!  Revision 1.3  2005-12-02 20:57:27  hb
!  Documentation updated and some bugs fixed
!
!  Revision 1.2  2004/07/30 09:22:20  hb
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
!  Here, the bio namelist {\tt bio\_template.nml} should be read and
!  various variables (rates and settling velocities)
!  should be transformed into SI units.
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
!  Here, the the initial conditions should be set and the settling 
!  velocities should be
!  transferred to all vertical levels. Non-negative concentrations 
!  should be declared  as non-negative variables, 
!  and it should be defined which variables would be
!  taken up by benthic filter feeders.
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
!  This subroutine provides information about the variable names as they
!  will be used when storing data in NetCDF files.
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
! !INTERFACE:
   subroutine light_template(nlev,bioshade_feedback)
!
! !DESCRIPTION:
! Here, the photosynthetically available radiation should be calculated
! by simply assuming that the short wave part of the total
! radiation is available for photosynthesis. The user should make
! sure that this is consistent with the light class given in the
! {\tt extinct} namelist of the {\tt obs.nml} file.
! The self-shading effect should also be calculated in this subroutine,
! which may be used to consider the effect of bio-turbidity also
! in the temperature equation (if {\tt bioshade\_feedback} is set
! to true in {\tt bio.nml}). For details, see section \ref{sec:do-bio}.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: nlev
   logical, intent(in)                 :: bioshade_feedback
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC

   STDERR 'TEMPLATE light properties are calculated here.'

   end subroutine light_template
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of geobiochemical model
!
! !INTERFACE:
   subroutine do_bio_template(numc,nlev)
!
! !DESCRIPTION:
! Here, the source and sink terms for the right hand sides need to be
! defined. Since this is a template file only, nothing is done here,
! and the execution of the program is terminated here. 
!
! !USES:
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
