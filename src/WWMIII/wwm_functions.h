#ifdef USE_SINGLE
# define MyREAL(xinp) REAL(xinp)
# define MySNGL(xinp) xinp
# define MySQRT(xinp) SQRT(xinp)
# define MyATAN2(xinp,yinp) ATAN2(xinp,yinp)
# define MyABS(xinp) ABS(xinp)
# define MyEXP(xinp) EXP(xinp)
# define MySINH(xinp) SINH(xinp)
# define MyCOSH(xinp) COSH(xinp)
# define MyTANH(xinp) TANH(xinp)
# define MySIN(xinp) SIN(xinp)
# define MyCOS(xinp) COS(xinp)
# define MyASIN(xinp) ASIN(xinp)
# define MyACOS(xinp) ACOS(xinp)
#else
# define MySNGL(xinp) SNGL(xinp)
# define MyREAL(xinp) DBLE(xinp)
# define MySQRT(xinp) DSQRT(xinp)
# define MyATAN2(xinp,yinp) DATAN2(xinp,yinp)
# define MyABS(xinp) DABS(xinp)
# define MyEXP(xinp) DEXP(xinp)
# define MySINH(xinp) DSINH(xinp)
# define MyCOSH(xinp) DCOSH(xinp)
# define MyTANH(xinp) DTANH(xinp)
# define MySIN(xinp) DSIN(xinp)
# define MyCOS(xinp) DCOS(xinp)
# define MyASIN(xinp) DASIN(xinp)
# define MyACOS(xinp) DACOS(xinp)
#endif
#if defined ST41 || defined ST42
# define ST_DEF
#endif
! The cppdefs creates lots of trouble
! since MPI is one symbol and it gets
! translated from use MPI ---> use 1
!
!#ifdef ROMS_WWM_PGMCL_COUPLING
!# include "cppdefs.h"
!#endif
#ifdef SCHISM
# define MPI_PARALL_GRID
#endif
#ifdef WWM_MPI
# define MPI_PARALL_GRID
#endif
