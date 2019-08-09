!$Id: tridiagonal.F90,v 1.7 2006-11-24 15:13:41 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mtridiagonal --- solving the system\label{sec:tridiagonal}
!
! !INTERFACE:
   MODULE mtridiagonal
!
! !DESCRIPTION:
!
!  Solves a linear system of equations with a tridiagonal matrix
!  using Gaussian elimination.
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_tridiagonal, tridiagonal, clean_tridiagonal
!
! !PUBLIC DATA MEMBERS:
   REALTYPE, dimension(:), allocatable     :: au,bu,cu,du
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!  $Log: tridiagonal.F90,v $
!  Revision 1.7  2006-11-24 15:13:41  kbk
!  de-allocate memory and close open files
!
!  Revision 1.6  2005/06/27 13:44:07  kbk
!  modified + removed traling blanks
!
!  Revision 1.5  2004/08/17 15:33:47  lars
!  removed tabs
!
!  Revision 1.4  2003/03/28 09:20:36  kbk
!  added new copyright to files
!
!  Revision 1.3  2003/03/28 08:06:33  kbk
!  removed tabs
!
!  Revision 1.2  2003/03/10 08:54:16  gotm
!  Improved documentation and cleaned up code
!
!  Revision 1.1.1.1  2001/02/12 15:55:58  gotm
!  initial import into CVS
!
!EOP
!
!  private data members
   REALTYPE, private, dimension(:),allocatable  ::  ru,qu
!
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Allocate memory
!
! !INTERFACE:
   subroutine init_tridiagonal(N)
!
! !DESCRIPTION:
!  This routines allocates memory necessary to perform the Gaussian
!  elimination.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: N
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: rc
!
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'init_tridiagonal'
   allocate(au(0:N),stat=rc)
   if (rc /= 0) stop 'init_tridiagonal: Error allocating au)'
   au = 0.

   allocate(bu(0:N),stat=rc)
   if (rc /= 0) stop 'init_tridiagonal: Error allocating bu)'
   bu = 0.

   allocate(cu(0:N),stat=rc)
   if (rc /= 0) stop 'init_tridiagonal: Error allocating cu)'
   cu = 0.

   allocate(du(0:N),stat=rc)
   if (rc /= 0) stop 'init_tridiagonal: Error allocating du)'
   du = 0.

   allocate(ru(0:N),stat=rc)
   if (rc /= 0) stop 'init_tridiagonal: Error allocating ru)'

   allocate(qu(0:N),stat=rc)
   if (rc /= 0) stop 'init_tridiagonal: Error allocating qu)'

   return
   end subroutine init_tridiagonal
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Simplified Gaussian elimination
!
! !INTERFACE:
   subroutine tridiagonal(N,fi,lt,value)
!
! !DESCRIPTION:
! A linear equation with tridiagonal matrix structure is solved here. The main
! diagonal is stored on {\tt bu}, the upper diagonal on {\tt au}, and the
! lower diagonal on {\tt cu}, the right hand side is stored on {\tt du}.
! The method used here is the simplified Gauss elimination, also called
! \emph{Thomas algorithm}.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: N,fi,lt
!
! !OUTPUT PARAMETERS:
   REALTYPE                            :: value(0:N)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!  $Log: tridiagonal.F90,v $
!  Revision 1.7  2006-11-24 15:13:41  kbk
!  de-allocate memory and close open files
!
!  Revision 1.6  2005/06/27 13:44:07  kbk
!  modified + removed traling blanks
!
!  Revision 1.5  2004/08/17 15:33:47  lars
!  removed tabs
!
!  Revision 1.4  2003/03/28 09:20:36  kbk
!  added new copyright to files
!
!  Revision 1.3  2003/03/28 08:06:33  kbk
!  removed tabs
!
!  Revision 1.2  2003/03/10 08:54:16  gotm
!  Improved documentation and cleaned up code
!
!  Revision 1.1.1.1  2001/02/12 15:55:58  gotm
!  initial import into CVS
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: i
!
!-----------------------------------------------------------------------
!BOC
   ru(lt)=au(lt)/bu(lt)
   qu(lt)=du(lt)/bu(lt)

   do i=lt-1,fi+1,-1
      ru(i)=au(i)/(bu(i)-cu(i)*ru(i+1))
      qu(i)=(du(i)-cu(i)*qu(i+1))/(bu(i)-cu(i)*ru(i+1))
   end do

   qu(fi)=(du(fi)-cu(fi)*qu(fi+1))/(bu(fi)-cu(fi)*ru(fi+1))

   value(fi)=qu(fi)
   do i=fi+1,lt
      value(i)=qu(i)-ru(i)*value(i-1)
   end do


   return
   end subroutine tridiagonal
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: De-allocate memory
!
! !INTERFACE:
   subroutine clean_tridiagonal()
!
! !DESCRIPTION:
!  De-allocates memory allocated in init\_tridiagonal.
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'clean_tridiagonal'

   if (allocated(au)) deallocate(au)

   if (allocated(bu)) deallocate(bu)

   if (allocated(cu)) deallocate(cu)

   if (allocated(du)) deallocate(du)

   if (allocated(ru)) deallocate(ru)

   if (allocated(qu)) deallocate(qu)

   return
   end subroutine clean_tridiagonal
!EOC
!-----------------------------------------------------------------------

   end module mtridiagonal

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
