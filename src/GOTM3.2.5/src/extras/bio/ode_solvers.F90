!$Id: ode_solvers.F90,v 1.5 2004-07-30 09:22:20 hb Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: General ODE solver
!
! !INTERFACE
   subroutine ode_solver(solver,numc,nlev,dt,h,cc,t)
!
! !DESCRIPTION
!
! !USES
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: solver,nlev,numc
   REALTYPE, intent(in)                :: dt
   REALTYPE, intent(in)                :: t(0:nlev)
   REALTYPE, intent(in)                :: h(0:nlev)
!
! !INPUT/OUTPUT PARAMETER:
   REALTYPE, intent(inout)             :: cc(1:numc,0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
   select case (solver)
      case (1)
         call euler_forward(dt,numc,nlev,cc,h,t)
      case (2)
         call runge_kutta_2(dt,numc,nlev,cc,h,t)
      case (3)
         call runge_kutta_4(dt,numc,nlev,cc,h,t)
      case (4)
         call patankar(dt,numc,nlev,cc,h,t)
      case (5)
         call patankar_runge_kutta_2(dt,numc,nlev,cc,h,t)
      case (6)
         call patankar_runge_kutta_4(dt,numc,nlev,cc,h,t)
      case (7)
         call modified_patankar(dt,numc,nlev,cc,h,t)
      case (8)
         call modified_patankar_2(dt,numc,nlev,cc,h,t)
      case (9)
         call modified_patankar_4(dt,numc,nlev,cc,h,t)
      case default
         stop "bio: no valid solver method specified in bio.inp !"
   end select

   return
   end subroutine ode_solver
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Euler-forward scheme for geobiochemical models
!
! !INTERFACE
   subroutine euler_forward(dt,numc,nlev,cc,h,t)
!
! !DESCRIPTION
!
! !USES
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,nlev
   REALTYPE, intent(in)                :: h(0:nlev)
   REALTYPE, intent(in)                :: t(0:nlev)
!
! !INPUT/OUTPUT PARAMETER:
  REALTYPE, intent(inout)              :: cc(1:numc,0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  logical  :: first
  REALTYPE :: rhs
  REALTYPE :: pp(1:numc,1:numc,0:nlev),dd(1:numc,1:numc,0:nlev)
  integer  :: i,j,ci
!EOP
!-----------------------------------------------------------------------
!BOC
   first=.true.

   call process_model(first,numc,nlev,cc,pp,dd,h,t)

   do ci=1,nlev
      do i=1,numc
         rhs=0.
         do j=1,numc
            rhs=rhs+pp(i,j,ci)-dd(i,j,ci)
         end do
         cc(i,ci)=cc(i,ci)+dt*rhs
      end do
   end do

   return
   end subroutine euler_forward
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Second-order Runge-Kutta scheme for geobiochemical models
!
! !INTERFACE
   subroutine runge_kutta_2(dt,numc,nlev,cc,h,t)
!
! !DESCRIPTION
!
! !USES
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,nlev
   REALTYPE, intent(in)                :: h(0:nlev)
   REALTYPE, intent(in)                :: t(0:nlev)
!
! !INPUT/OUTPUT PARAMETER:
   REALTYPE, intent(inout)             :: cc(1:numc,0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  logical  :: first
  REALTYPE :: rhs(1:numc,0:nlev),rhs1(1:numc)
  REALTYPE :: cc1(1:numc,0:nlev)
  REALTYPE :: pp(1:numc,1:numc,0:nlev),dd(1:numc,1:numc,0:nlev)
  integer  :: i,j,ci
!EOP
!-----------------------------------------------------------------------
!BOC
   first=.true.
   call process_model(first,numc,nlev,cc,pp,dd,h,t)

   do ci=1,nlev
      do i=1,numc
         rhs(i,ci)=0.
         do j=1,numc
            rhs(i,ci)=rhs(i,ci)+pp(i,j,ci)-dd(i,j,ci)
         end do
         cc1(i,ci)=cc(i,ci)+dt*rhs(i,ci)
      end do
   end do

   call process_model(first,numc,nlev,cc1,pp,dd,h,t)

   do ci=1,nlev
      do i=1,numc
         rhs1(i)=0.
         do j=1,numc
            rhs1(i)=rhs1(i)+pp(i,j,ci)-dd(i,j,ci)
         end do
         cc(i,ci)=cc(i,ci)+dt*0.5*(rhs(i,ci)+rhs1(i))
      end do
   end do

   return
   end subroutine runge_kutta_2
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Fourth-order Runge-Kutta scheme for geobiochemical models
!
! !INTERFACE
   subroutine runge_kutta_4(dt,numc,nlev,cc,h,t)
!
! !DESCRIPTION
!
! !USES
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,nlev
   REALTYPE, intent(in)                :: h(0:nlev)
   REALTYPE, intent(in)                :: t(0:nlev)
!
! !INPUT/OUTPUT PARAMETER:
   REALTYPE, intent(inout)             :: cc(1:numc,0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  logical  :: first
  REALTYPE :: rhs(1:numc,0:nlev),rhs1(1:numc,0:nlev)
  REALTYPE :: rhs2(1:numc,0:nlev),rhs3(1:numc,0:nlev)
  REALTYPE :: cc1(1:numc,0:nlev)
  REALTYPE :: pp(1:numc,1:numc,0:nlev),dd(1:numc,1:numc,0:nlev)
  integer  :: i,j,ci
!EOP
!-----------------------------------------------------------------------
!BOC
   first=.true.
   call process_model(first,numc,nlev,cc,pp,dd,h,t)

   do ci=1,nlev
      do i=1,numc
         rhs(i,ci)=0.
         do j=1,numc
            rhs(i,ci)=rhs(i,ci)+pp(i,j,ci)-dd(i,j,ci)
         end do
         cc1(i,ci)=cc(i,ci)+dt*rhs(i,ci)
      end do
   end do

   call process_model(first,numc,nlev,cc1,pp,dd,h,t)

   do ci=1,nlev
      do i=1,numc
         rhs1(i,ci)=0.
         do j=1,numc
            rhs1(i,ci)=rhs1(i,ci)+pp(i,j,ci)-dd(i,j,ci)
         end do
         cc1(i,ci)=cc(i,ci)+dt*rhs1(i,ci)
      end do
   end do

   call process_model(first,numc,nlev,cc1,pp,dd,h,t)

   do ci=1,nlev
      do i=1,numc
         rhs2(i,ci)=0.
         do j=1,numc
            rhs2(i,ci)=rhs2(i,ci)+pp(i,j,ci)-dd(i,j,ci)
         end do
         cc1(i,ci)=cc(i,ci)+dt*rhs2(i,ci)
      end do
   end do

   call process_model(first,numc,nlev,cc1,pp,dd,h,t)

   do ci=1,nlev
      do i=1,numc
         rhs3(i,ci)=0.
         do j=1,numc
            rhs3(i,ci)=rhs3(i,ci)+pp(i,j,ci)-dd(i,j,ci)
         end do
         cc(i,ci)=cc(i,ci)+dt*1./3.*(0.5*rhs(i,ci)+rhs1(i,ci)+rhs2(i,ci)+0.5*rhs3(i,ci))
      end do
   end do

   return
   end subroutine runge_kutta_4
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Patankar scheme for geobiochemical models
!
! !INTERFACE
   subroutine patankar(dt,numc,nlev,cc,h,t)
!
! !DESCRIPTION
!
! !USES
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,nlev
   REALTYPE, intent(in)                :: t(0:nlev)
   REALTYPE, intent(in)                :: h(0:nlev)
!
! !INPUT/OUTPUT PARAMETER:
   REALTYPE, intent(inout)             :: cc(1:numc,0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  logical  :: first
  REALTYPE :: ppsum,ddsum
  REALTYPE :: pp(1:numc,1:numc,0:nlev),dd(1:numc,1:numc,0:nlev)
  integer  :: i,j,ci
!EOP
!-----------------------------------------------------------------------
!BOC
   first=.true.
   call process_model(first,numc,nlev,cc,pp,dd,h,t)
   do ci=1,nlev
      do i=1,numc
         ppsum=0.
         ddsum=0.
         do j=1,numc
            ppsum=ppsum+pp(i,j,ci)
            ddsum=ddsum+dd(i,j,ci)
         end do
         cc(i,ci)=(cc(i,ci)+dt*ppsum)/(1.+dt*ddsum/cc(i,ci))
      end do
   end do

   return
   end subroutine patankar
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Patankar-Runge-Kutta (2nd-order) scheme for geobiochemical
!  models
!
! !INTERFACE
   subroutine patankar_runge_kutta_2(dt,numc,nlev,cc,h,t)
!
! !DESCRIPTION
!
! !USES
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,nlev
   REALTYPE, intent(in)                :: h(0:nlev)
   REALTYPE, intent(in)                :: t(0:nlev)
!
! !INPUT/OUTPUT PARAMETER:
   REALTYPE, intent(inout)             :: cc(1:numc,0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  logical  :: first
  REALTYPE :: ppsum(1:numc,0:nlev),ddsum(1:numc,0:nlev)
  REALTYPE :: pp(1:numc,1:numc,0:nlev),dd(1:numc,1:numc,0:nlev)
  REALTYPE :: cc1(1:numc,0:nlev)
  integer  :: i,j,ci
!EOP
!-----------------------------------------------------------------------
!BOC
   first=.true.
   call process_model(first,numc,nlev,cc,pp,dd,h,t)

   do ci=1,nlev
      do i=1,numc
         ppsum(i,ci)=0.
         ddsum(i,ci)=0.
         do j=1,numc
            ppsum(i,ci)=ppsum(i,ci)+pp(i,j,ci)
            ddsum(i,ci)=ddsum(i,ci)+dd(i,j,ci)
         end do
         cc1(i,ci)=(cc(i,ci)+dt*ppsum(i,ci))/(1.+dt*ddsum(i,ci)/cc(i,ci))
      end do
   end do

   call process_model(first,numc,nlev,cc1,pp,dd,h,t)

   do ci=1,nlev
      do i=1,numc
         do j=1,numc
            ppsum(i,ci)=ppsum(i,ci)+pp(i,j,ci)
            ddsum(i,ci)=ddsum(i,ci)+dd(i,j,ci)
         end do
         cc(i,ci)=(cc(i,ci)+0.5*dt*ppsum(i,ci))/(1.+0.5*dt*ddsum(i,ci)/cc1(i,ci))
      end do
   end do

   return
   end subroutine patankar_runge_kutta_2
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Patankar-Runge-Kutta (4th-order) scheme for geobiochemical
!  models
!
! !INTERFACE
   subroutine patankar_runge_kutta_4(dt,numc,nlev,cc,h,t)
!
! !DESCRIPTION
!
! !USES
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,nlev
   REALTYPE, intent(in)                :: h(0:nlev)
   REALTYPE, intent(in)                :: t(0:nlev)
!
! !INPUT/OUTPUT PARAMETER:
   REALTYPE, intent(inout)             :: cc(1:numc,0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  logical  :: first
  REALTYPE :: ppsum(1:numc,0:nlev),ddsum(1:numc,0:nlev)
  REALTYPE :: ppsum1(1:numc,0:nlev),ddsum1(1:numc,0:nlev)
  REALTYPE :: ppsum2(1:numc,0:nlev),ddsum2(1:numc,0:nlev)
  REALTYPE :: ppsum3(1:numc,0:nlev),ddsum3(1:numc,0:nlev)
  REALTYPE :: pp(1:numc,1:numc,0:nlev),dd(1:numc,1:numc,0:nlev)
  REALTYPE :: cc1(1:numc,0:nlev)
  integer  :: i,j,ci
!EOP
!-----------------------------------------------------------------------
!BOC
   first=.true.
   call process_model(first,numc,nlev,cc,pp,dd,h,t)

   do ci=1,nlev
      do i=1,numc
         ppsum(i,ci)=0.
         ddsum(i,ci)=0.
         do j=1,numc
            ppsum(i,ci)=ppsum(i,ci)+pp(i,j,ci)
            ddsum(i,ci)=ddsum(i,ci)+dd(i,j,ci)
         end do
         cc1(i,ci)=(cc(i,ci)+dt*ppsum(i,ci))/(1.+dt*ddsum(i,ci)/cc(i,ci))
      end do
   end do

   call process_model(first,numc,nlev,cc1,pp,dd,h,t)

   do ci=1,nlev
      do i=1,numc
         ppsum1(i,ci)=0.
         ddsum1(i,ci)=0.
         do j=1,numc
            ppsum1(i,ci)=ppsum1(i,ci)+pp(i,j,ci)
            ddsum1(i,ci)=ddsum1(i,ci)+dd(i,j,ci)
         end do
         cc1(i,ci)=(cc(i,ci)+dt*ppsum1(i,ci))/(1.+dt*ddsum1(i,ci)/cc1(i,ci))
      end do
   end do

   call process_model(first,numc,nlev,cc1,pp,dd,h,t)

   do ci=1,nlev
      do i=1,numc
         ppsum2(i,ci)=0.
         ddsum2(i,ci)=0.
         do j=1,numc
            ppsum2(i,ci)=ppsum2(i,ci)+pp(i,j,ci)
            ddsum2(i,ci)=ddsum2(i,ci)+dd(i,j,ci)
         end do
         cc1(i,ci)=(cc(i,ci)+dt*ppsum2(i,ci))/(1.+dt*ddsum2(i,ci)/cc1(i,ci))
      end do
   end do

   call process_model(first,numc,nlev,cc1,pp,dd,h,t)

   do ci=1,nlev
      do i=1,numc
         ppsum3(i,ci)=0.
         ddsum3(i,ci)=0.
         do j=1,numc
            ppsum3(i,ci)=ppsum3(i,ci)+pp(i,j,ci)
            ddsum3(i,ci)=ddsum3(i,ci)+dd(i,j,ci)
         end do
         ppsum(i,ci)=1./3.*(0.5*ppsum(i,ci)+ppsum1(i,ci)+ppsum2(i,ci)+0.5*ppsum3(i,ci))
         ddsum(i,ci)=1./3.*(0.5*ddsum(i,ci)+ddsum1(i,ci)+ddsum2(i,ci)+0.5*ddsum3(i,ci))
         cc(i,ci)=(cc(i,ci)+dt*ppsum(i,ci))/(1.+dt*ddsum(i,ci)/cc1(i,ci))
      end do
   end do

   return
   end subroutine patankar_runge_kutta_4
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Modified Patankar scheme for geobiochemical models
!
! !INTERFACE
   subroutine modified_patankar(dt,numc,nlev,cc,h,t)
!
! !DESCRIPTION
!
! !USES
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,nlev
   REALTYPE, intent(in)                :: h(0:nlev)
   REALTYPE, intent(in)                :: t(0:nlev)
!
! !INPUT/OUTPUT PARAMETER:
   REALTYPE, intent(inout)             :: cc(1:numc,0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  logical  :: first
  REALTYPE :: ppsum,ddsum
  REALTYPE :: pp(1:numc,1:numc,0:nlev),dd(1:numc,1:numc,0:nlev)
  REALTYPE :: a(1:numc,1:numc),r(1:numc)
  integer  :: i,j,ci
!EOP
!-----------------------------------------------------------------------
!BOC
   first=.true.
   call process_model(first,numc,nlev,cc,pp,dd,h,t)

   do ci=1,nlev
      do i=1,numc
         a(i,i)=0.
         do j=1,numc
            a(i,i)=a(i,i)+dd(i,j,ci)
            if (i.ne.j) a(i,j)=-dt*pp(i,j,ci)/cc(j,ci)
         end do
         a(i,i)=dt*a(i,i)/cc(i,ci)
         a(i,i)=1.+a(i,i)
         r(i)=cc(i,ci)+dt*pp(i,i,ci)
      end do
      call matrix(numc,a,r,cc(:,ci))
   end do

   return
   end subroutine modified_patankar
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Modified Patankar-Runge-Kutta (2nd-order) scheme for
!  geobiochemical models
!
! !INTERFACE
   subroutine modified_patankar_2(dt,numc,nlev,cc,h,t)
!
! !DESCRIPTION
!
! !USES
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,nlev
   REALTYPE, intent(in)                :: h(0:nlev)
   REALTYPE, intent(in)                :: t(0:nlev)
!
! !INPUT/OUTPUT PARAMETER:
   REALTYPE, intent(inout)             :: cc(1:numc,0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  logical  :: first
  REALTYPE :: ppsum,ddsum
  REALTYPE :: pp(1:numc,1:numc,0:nlev),dd(1:numc,1:numc,0:nlev)
  REALTYPE :: pp1(1:numc,1:numc,0:nlev),dd1(1:numc,1:numc,0:nlev)
  REALTYPE :: a(1:numc,1:numc),r(1:numc)
  REALTYPE :: cc1(1:numc,0:nlev)
  integer  :: i,j,ci
!EOP
!-----------------------------------------------------------------------
!BOC
   first=.true.
   call process_model(first,numc,nlev,cc,pp,dd,h,t)

   do ci=1,nlev
      do i=1,numc
         a(i,i)=0.
         do j=1,numc
            a(i,i)=a(i,i)+dd(i,j,ci)
            if (i.ne.j) a(i,j)=-dt*pp(i,j,ci)/cc(j,ci)
         end do
         a(i,i)=dt*a(i,i)/cc(i,ci)
         a(i,i)=1.+a(i,i)
         r(i)=cc(i,ci)+dt*pp(i,i,ci)
      end do
      call matrix(numc,a,r,cc1(:,ci))
   end do

   call process_model(first,numc,nlev,cc1,pp1,dd1,h,t)

   pp=0.5*(pp+pp1)
   dd=0.5*(dd+dd1)

   do ci=1,nlev
      do i=1,numc
         a(i,i)=0.
         do j=1,numc
            a(i,i)=a(i,i)+dd(i,j,ci)
            if (i.ne.j) a(i,j)=-dt*pp(i,j,ci)/cc1(j,ci)
         end do
         a(i,i)=dt*a(i,i)/cc1(i,ci)
         a(i,i)=1.+a(i,i)
         r(i)=cc(i,ci)+dt*pp(i,i,ci)
      end do
      call matrix(numc,a,r,cc(:,ci))
   end do

   return
   end subroutine modified_patankar_2
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Modified Patankar-Runge-Kutta (4th-order) scheme for
!  geobiochemical models
!
! !INTERFACE
   subroutine modified_patankar_4(dt,numc,nlev,cc,h,t)
!
! !DESCRIPTION
!
! !USES
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
   integer, intent(in)                 :: numc,nlev
   REALTYPE, intent(in)                :: h(0:nlev)
   REALTYPE, intent(in)                :: t(0:nlev)
!
! !INPUT/OUTPUT PARAMETER:
   REALTYPE, intent(inout)             :: cc(1:numc,0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  logical  :: first
  REALTYPE :: ppsum,ddsum
  REALTYPE :: pp(1:numc,1:numc,0:nlev),dd(1:numc,1:numc,0:nlev)
  REALTYPE :: pp1(1:numc,1:numc,0:nlev),dd1(1:numc,1:numc,0:nlev)
  REALTYPE :: pp2(1:numc,1:numc,0:nlev),dd2(1:numc,1:numc,0:nlev)
  REALTYPE :: pp3(1:numc,1:numc,0:nlev),dd3(1:numc,1:numc,0:nlev)
  REALTYPE :: a(1:numc,1:numc),r(1:numc)
  REALTYPE :: cc1(1:numc,0:nlev)
  integer  :: i,j,ci
!EOP
!-----------------------------------------------------------------------
!BOC
   first=.true.
   call process_model(first,numc,nlev,cc,pp,dd,h,t)

   do ci=1,nlev
      do i=1,numc
         a(i,i)=0.
         do j=1,numc
            a(i,i)=a(i,i)+dd(i,j,ci)
            if (i.ne.j) a(i,j)=-dt*pp(i,j,ci)/cc(j,ci)
         end do
         a(i,i)=dt*a(i,i)/cc(i,ci)
         a(i,i)=1.+a(i,i)
         r(i)=cc(i,ci)+dt*pp(i,i,ci)
      end do
      call matrix(numc,a,r,cc1(:,ci))
   end do

   call process_model(first,numc,nlev,cc1,pp1,dd1,h,t)

   do ci=1,nlev
      do i=1,numc
         a(i,i)=0.
         do j=1,numc
            a(i,i)=a(i,i)+dd1(i,j,ci)
            if (i.ne.j) a(i,j)=-dt*pp1(i,j,ci)/cc1(j,ci)
         end do
         a(i,i)=dt*a(i,i)/cc1(i,ci)
         a(i,i)=1.+a(i,i)
         r(i)=cc(i,ci)+dt*pp1(i,i,ci)
      end do
      call matrix(numc,a,r,cc1(:,ci))
   end do

   call process_model(first,numc,nlev,cc1,pp2,dd2,h,t)

   do ci=1,nlev
      do i=1,numc
         a(i,i)=0.
         do j=1,numc
            a(i,i)=a(i,i)+dd2(i,j,ci)
            if (i.ne.j) a(i,j)=-dt*pp2(i,j,ci)/cc1(j,ci)
         end do
         a(i,i)=dt*a(i,i)/cc1(i,ci)
         a(i,i)=1.+a(i,i)
         r(i)=cc(i,ci)+dt*pp2(i,i,ci)
      end do
      call matrix(numc,a,r,cc1(:,ci))
   end do

   call process_model(first,numc,nlev,cc1,pp3,dd3,h,t)

   pp=1./3.*(0.5*pp+pp1+pp2+0.5*pp3)
   dd=1./3.*(0.5*dd+dd1+dd2+0.5*dd3)

   do ci=1,nlev
      do i=1,numc
         a(i,i)=0.
         do j=1,numc
            a(i,i)=a(i,i)+dd(i,j,ci)
            if (i.ne.j) a(i,j)=-dt*pp(i,j,ci)/cc1(j,ci)
         end do
         a(i,i)=dt*a(i,i)/cc1(i,ci)
         a(i,i)=1.+a(i,i)
         r(i)=cc(i,ci)+dt*pp(i,i,ci)
      end do
      call matrix(numc,a,r,cc(:,ci))
   end do

   return
   end subroutine modified_patankar_4
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Matrix solver
!
! !INTERFACE
   subroutine matrix(n,a,r,c)
!
! !DESCRIPTION
!
! !USES
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: n
!
! INPUT/OUTPUT PARAMETERS:
  REALTYPE                             :: a(1:n,1:n),r(1:n)
!
! OUTPUT PARAMETERS:
  REALTYPE, intent(out)                :: c(1:n)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  integer  :: i,j,k
!EOP
!-----------------------------------------------------------------------
!BOC
   do i=1,n
      r(i)=r(i)/a(i,i)
      do j=n,i,-1
         a(i,j)=a(i,j)/a(i,i)
      end do
      do k=i+1,n
         r(k)=r(k)-a(k,i)*r(i)
         do j=i+1,n
            a(k,j)=a(k,j)-a(k,i)*a(i,j)
         end do
      end do
   end do

   do i=n,1,-1
      c(i)=r(i)
      do j=i+1,n
         c(i)=c(i)-a(i,j)*c(j)
      end do
   end do

   return
   end subroutine matrix
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
