!$Id: bio_sed.F90,v 1.3 2004-08-02 08:34:36 hb Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bio_sed --- simple sedimentation model \label{sec:bio_sed}
!
! !INTERFACE:
   module bio_sed
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
   public init_bio_sed, init_var_sed, var_info_sed, &
          do_bio_sed, end_bio_sed
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:!
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: bio_sed.F90,v $
!  Revision 1.3  2004-08-02 08:34:36  hb
!  updated init routines to reflect new internal bio interface
!
!  Revision 1.2  2004/07/30 09:22:20  hb
!  use bio_var in specific bio models - simpliefied internal interface
!
!  Revision 1.1  2003/10/28 10:22:45  hb
!  added support for sedimentation only 1 compartment bio model
!
!
! !LOCAL VARIABLES:
!  from a namelist
   REALTYPE                  :: C_initial=4.5
   REALTYPE                  :: w_C=-5.787037e-05
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the bio module
!
! !INTERFACE:
   subroutine init_bio_sed(namlst,fname,unit)
!
! !DESCRIPTION:
!  Here, the bio namelist {\tt bio_sed.inp} is read and memory is
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
! !LOCAL VARIABLES:
   namelist /bio_sed_nml/ numc,C_initial,w_C
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL2 'init_bio_sed'

   open(namlst,file=fname,action='read',status='old',err=98)
   read(namlst,nml=bio_sed_nml,err=99)
   close(namlst)

   LEVEL3 'Sedimentation bio module initialised ...'

   w_C = w_C/86400.

   return

98 LEVEL2 'I could not open bio_sed.inp'
   LEVEL2 'If thats not what you want you have to supply bio_sed.inp'
   LEVEL2 'See the bio example on www.gotm.net for a working bio_sed.inp'
   return
99 FATAL 'I could not read bio_sed.inp'
   stop 'init_bio_sed'
   end subroutine init_bio_sed
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the concentration variables
!
! !INTERFACE:
   subroutine init_var_sed(nlev)
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

!EOP
!-----------------------------------------------------------------------
!BOC
   cc(1,:)=C_initial

   ws(1,:) = w_C

   mussels_inhale(1) = .true.


   LEVEL3 'Sedimentation variables initialised ...'

   return

   end subroutine init_var_sed
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Providing info on variables
!
! !INTERFACE:
   subroutine var_info_sed()
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
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
   var_names(1) = 'conc'
   var_units(1) = 'fractions'
   var_long(1) = 'concentration'

   return
   end subroutine var_info_sed
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of geobiochemical model
!
! !INTERFACE
   subroutine do_bio_sed()
!
! !DESCRIPTION
!
! !USES
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC

!  no right hand sides necessary

   return
   end subroutine do_bio_sed
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finish the bio calculations
!
! !INTERFACE:
   subroutine end_bio_sed
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
   end subroutine end_bio_sed
!EOC

!-----------------------------------------------------------------------

   end module bio_sed

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
