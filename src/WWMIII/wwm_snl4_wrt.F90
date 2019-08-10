!----------------------------------------------------------------------
subroutine wwmquad_wrt (aquad,sigma,dir,ndir,nsig,depth,iquad,xnl,diag,ierr)
!----------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant Ph. van Vledder
!   |   +---+
!   |   | +---+  Last update:  9 September 2002
!   +---+ |   |  Release: 1.0
!         +---+
!
!  +--------+ +-+   +-+ +------+  Unversity of Technology Darmstadt
!  |        | | |   | | |+--+  |  Institute for Hydraulic Engineering and Water Resources Managment   
!  +--|  |--+ | |   | | ||  |  |  
!     |  |    | |   | | ||  |  |  Aaron Roland
!     |  |    | +---+ | |+--+  |    
!     +--+    +-------+ +------+  Last update: 29.12.2005 
!
use serv_xnl4v5
use m_xnldata
!
implicit none
!-----------------------------------------------------------------------------------
!
!  0. Update history
!
!     29/12/2005 Implementation in the WindWaveModel, Aaron Roland, Vrdnik, Yugoslawija
!
!  1. Purpose:
!
!     interface with WWM model to compute nonlinear transfer for
!     given action density spectrum
!
!  2. Method
!
!     Resio/Tracy deep water geometric scaling
!     Rewritten by Gerbrant van Vledder
!
!  3. Parameter list:
!
! Type     I/O         Name      Description
!----------------------------------------------------------------------------------
integer, intent(in) :: ndir     ! number of directions
integer, intent(in) :: nsig     ! number of sigma values in array
integer, intent(in) :: iquad    ! method of computing nonlinear quadruplet interactions
!
real, intent(in)    :: aquad(nsig,ndir)        ! action density spectrum as a function of (sigma,dir)
real, intent(in)    :: sigma(nsig)             ! Intrinsic frequencies
real, intent(in)    :: dir(ndir)               ! directions in radians 
real, intent(in)    :: depth                   ! local depth 
!
!
integer, intent(out):: ierr                    ! Error indicator. If no errors are detected IERR=0
!
real, intent(out)   :: xnl(nsig,ndir)          ! transfer rate dA/dt(sigma,dir)
real, intent(out)   :: diag(nsig,ndir)         ! diagonal term (dXnl/dA)
!                                              ! a certain exact method (sigma,dir)
!--------------------------------------------------------------------------------------------------
!
!  4. Error messages
!
!     An error message is produced within the QUAD system.
!     If no errors are detected ierror = 0 
!     1, incorrect iquad
!
!  5. Called by
!
!     Program XNLTEST4
!
!  6. Subroutines used
!
!  7. Remarks
!
!     Old swan call routine ... rewritten by Aaron Roland
!
!  8. Structure
!
!  9. Switches
!
! 10. Source code
!------------------------------------------------------------------------------
!     Local parameters
!
integer isig,idir            ! counters
integer igrid                ! grid index
integer iproc                ! processor number for MPI and file naming
!
!------------------------------------------------------------------------------
!
! initialisations
!
ierr    = 0
diag    = 0.
!
!  check value of iquad
!
! iquad = 1   ! deep water transfer
! iquad = 2   ! deep water transfer with WAM depth scaling
! iquad = 3   ! finite depth transfer
!
if (iquad .lt. 1 .or. iquad .gt. 3) then
  ierr  = 1
end if
!
xnl  = 0.
diag = 0.
iproc = 0
!
!
!WRITE (*,*) 'ACLOC' 
!WRITE (*,*) aquad
!WRITE (*,*) 'SIGMA'
!WRITE (*,*) sigma 
!WRITE (*,*) 'SPDIR'
!WRITE (*,*) dir 
!WRITE (*,*) 'nsig & ndir'
!WRITE (*,*) nsig, ndir 
!WRITE (*,*) 'dep'
!WRITE (*,*) depth
!WRITE (*,*) 'iquad'
!WRITE (*,*)  iquad
!WRITE (*,*) 'iproc'
!WRITE (*,*) iproc 
!WRITE (*,*) 'ierr'
!WRITE (*,*) ierr 
!
call xnl_main(aquad,sigma,dir,nsig,ndir,depth,iquad,xnl,diag,iproc,ierr)
!
return
end subroutine
