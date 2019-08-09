#include "wwm_functions.h"
#ifdef VDISLIN
!/****************************************************************/
!/**                       DISLIN.F90                           **/
!/**                                                            **/
!/** Module file for DISLIN Fortran 90.                         **/
!/**                                                            **/
!/** Date     :  29.05.2009                                     **/
!/** Routines :  694                                            **/
!/** Version  :  9.5 A / explicit-shape                         **/
!/****************************************************************/

module dislin
  interface
    subroutine abs3pt(x,y,z,xp,yp)
      implicit none
      real, intent (in)  :: x,y,z
      real, intent (out) :: xp,yp
    end subroutine abs3pt
 
    subroutine addlab(cstr,v,itic,cax)
      implicit none
      character (len = *), intent (in) :: cstr,cax
      real, intent (in) :: v
      integer, intent (in) :: itic
    end subroutine addlab

    subroutine angle(i)
      implicit none
      integer, intent (in) :: i
    end subroutine angle
 
    subroutine arcell(nx,ny,na,nb,alpha,beta,theta)
      implicit none
      integer, intent (in) :: nx,ny,na,nb
      real, intent (in)   :: alpha,beta,theta
    end subroutine arcell
 
    subroutine areaf(ix,iy,n)
      implicit none
      integer, intent (in) :: n
      integer, dimension (n), intent (in) :: ix,iy
    end subroutine areaf
 
    subroutine autres(i,j)
      implicit none
      integer, intent (in) :: i,j
    end subroutine autres
 
    subroutine ax2grf()
    end subroutine ax2grf
 
    subroutine ax3len(i,j,k)
      implicit none
      integer, intent (in) :: i,j,k
    end subroutine ax3len
 
    subroutine axclrs(n,copt,cax)
      implicit none
      integer, intent(in) :: n
      character (len = *) , intent (in) :: copt, cax
    end subroutine axclrs
 
    subroutine axends(copt,cax)
      implicit none
      character (len = *) , intent (in) :: copt, cax
    end subroutine axends
 
    subroutine axgit()
    end subroutine axgit
 
    subroutine axis3d(x,y,z)
      implicit none
      real, intent (in) :: x,y,z
    end subroutine axis3d
 
    subroutine axsbgd(n)
      implicit none
      integer, intent (in) :: n
    end subroutine axsbgd
 
    subroutine axslen(i,j)
      implicit none
      integer, intent (in) :: i,j
    end subroutine axslen
 
    subroutine axsorg(i,j)
      implicit none
      integer, intent (in) :: i,j
    end subroutine axsorg
 
    subroutine axspos(i,j)
      implicit none
      integer, intent (in) :: i,j
    end subroutine axspos
 
    subroutine axsscl(copt,cax)
      implicit none
      character (len = *), intent (in) :: copt,cax
    end subroutine axsscl
 
    subroutine axstyp(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine axstyp

    subroutine barbor(iclr)
      implicit none
      integer, intent (in) :: iclr
    end subroutine barbor

    subroutine barclr(ic1,ic2,ic3)
      implicit none
      integer, intent (in) :: ic1,ic2,ic3
    end subroutine barclr
 
    subroutine bargrp(n,xgap)
      implicit none
      integer, intent (in) :: n
      real, intent (in) :: xgap
    end subroutine bargrp

    subroutine barmod(cmode,copt)
      implicit none
      character (len = *), intent (in) :: cmode,copt
    end subroutine barmod

    subroutine baropt(x1,x2)
      implicit none
      real, intent (in) :: x1,x2
    end subroutine baropt
 
    subroutine barpos(cpos)
      implicit none
      character (len = *), intent (in) :: cpos
    end subroutine barpos
 
    subroutine bars(xray,y1ray,y2ray,n)
      implicit none
      integer, intent (in) :: n
      real, intent (in out), dimension (n) :: xray,y1ray,y2ray
    end subroutine bars
 
    subroutine bars3d(xray,yray,z1ray,z2ray,xwray,ywray,icray,n)
      implicit none
      integer, intent (in) :: n
      real, intent (in), dimension (n) :: xray,yray,z1ray,z2ray,xwray,ywray
      integer, intent (in), dimension (n) :: icray
    end subroutine bars3d

    subroutine bartyp(ctyp)
      implicit none
      character (len = *), intent (in) :: ctyp
    end subroutine bartyp
 
    subroutine barwth(fact)
      implicit none
      real, intent (in) :: fact
    end subroutine barwth
 
    subroutine basalf(calph)
      implicit none
      character (len = *), intent (in) :: calph
    end subroutine basalf
 
    subroutine basdat(id,im,iy)
      implicit none
      integer, intent (in) :: id, im, iy
    end subroutine basdat

    subroutine bezier(xray,yray,nray,x,y,n)
      implicit none
      integer, intent (in)  :: nray,n
      real, dimension (nray), intent (in)  :: xray, yray
      real, dimension (n), intent (out) :: x, y
    end subroutine bezier
 
    subroutine bfcclr(ic)
      implicit none
      integer, intent (in) :: ic
    end subroutine bfcclr

    subroutine bfcmsh(ic)
      implicit none
      integer, intent (in) :: ic
    end subroutine bfcmsh

    subroutine bitsi2(nbits,mher,iher,mhin,ihin,lob)
      implicit none
      integer, intent (in) :: nbits,iher,ihin,lob
      integer (kind=selected_int_kind(4)), intent (in) :: mher
      integer (kind=selected_int_kind(4)), intent (in out) :: mhin
    end subroutine bitsi2
 
    subroutine bitsi4(nbits,mher,iher,mhin,ihin,lob)
      implicit none
      integer, intent (in) :: nbits,mher,iher,ihin,lob
      integer, intent (in out) :: mhin
    end subroutine bitsi4
 
    subroutine bmpfnt(cfnt)
      implicit none
      character (len = *), intent (in) :: cfnt
    end subroutine bmpfnt

    subroutine bmpmod(n,cval,copt)
      implicit none
      character (len = *), intent (in) :: cval,copt
      integer, intent (in) :: n
    end subroutine bmpmod

    subroutine box2d()
    end subroutine box2d
 
    subroutine box3d()
    end subroutine box3d
 
    subroutine center()
    end subroutine center
 
    subroutine cgmbgd(xr,xg,xb)
      implicit none
      real, intent (in) :: xr,xg,xb
    end subroutine cgmbgd
 
    subroutine cgmpic(ct)
      implicit none
      character (len = *), intent (in) :: ct
    end subroutine cgmpic
 
    subroutine cgmver(n)
      implicit none
      integer, intent (in) :: n
    end subroutine cgmver
 
    subroutine chaang(x)
      implicit none
      real, intent (in) :: x
    end subroutine chaang
 
    subroutine chacod(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine chacod

    subroutine chaspc(x)
      implicit none
      real, intent (in) :: x
    end subroutine chaspc
 
    subroutine chawth(x)
      implicit none
      real, intent (in) :: x
    end subroutine chawth
 
    subroutine chnatt()
    end subroutine chnatt
 
    subroutine chncrv(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine chncrv
 
    subroutine chndot()
    end subroutine chndot
 
    subroutine chndsh()
    end subroutine chndsh

    subroutine chnbar(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine chnbar
 
    subroutine chnpie(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine chnpie

    subroutine circ3p(x1,y1,x2,y2,x3,y3,xm,ym,r)
      implicit none
      real, intent (in) :: x1,y1,x2,y2,x3,y3
      real, intent (out) :: xm,ym,r
    end subroutine circ3p
 
    subroutine circle(nx,ny,nr)
      implicit none
      integer, intent (in) :: nx,ny,nr
    end subroutine circle
 
    subroutine circsp(n)
      implicit none
      integer, intent (in) :: n
    end subroutine circsp
 
    subroutine clip3d(ctyp)
      implicit none
      character (len = *), intent (in) :: ctyp
    end subroutine clip3d
 
    subroutine closfl(nlu)
      implicit none
      integer, intent (in) :: nlu
    end subroutine closfl
 
    subroutine clpbor(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine clpbor
 
    subroutine clpmod(cmod)
      implicit none
      character (len = *), intent (in) :: cmod
    end subroutine clpmod
 
    subroutine clpwin(nx,ny,nw,nh)
      implicit none
      integer, intent (in) :: nx,ny,nw,nh
    end subroutine clpwin
 
    subroutine clrcyc(i,iclr)
      implicit none
      integer, intent (in) :: i,iclr
    end subroutine clrcyc
 
    subroutine clrmod(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine clrmod
 
    subroutine clswin(id)
      implicit none
      integer, intent (in) :: id
    end subroutine clswin
 
    subroutine color(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine color
 
    subroutine colran(i,j)
      implicit none
      integer, intent (in) :: i,j
    end subroutine colran
 
    subroutine colray(z,ncol,n)
      implicit none
      integer, intent (in) :: n
      real, dimension (n), intent (in) :: z
      integer, dimension (n), intent (out) :: ncol
    end subroutine colray
 
    subroutine complx()
    end subroutine complx

    subroutine conclr(iray,n)
      implicit none
      integer, intent (in) :: n
      integer, dimension (n), intent (in) :: iray
    end subroutine conclr
 
    subroutine concrv(x,y,n,zlev)
      implicit none
      integer, intent (in) :: n
      real, dimension (n), intent (in) :: x,y
      real, intent (in) :: zlev
    end subroutine concrv

    subroutine cone3d(x,y,z,r,h1,h2,nsk1,nsk2)
      implicit none
      real, intent (in) :: x,y,z,r,h1,h2
      integer, intent (in) :: nsk1,nsk2
    end subroutine cone3d

    subroutine confll(xray,yray,zray,n,i1ray,i2ray,i3ray,ntri,zlev,nlev)
      implicit none
      integer, intent (in) :: n,ntri,nlev
      real, dimension (n), intent (in) :: xray,yray,zray
      real, dimension (nlev), intent (in) :: zlev
      integer, dimension (ntri), intent (in) :: i1ray,i2ray,i3ray
    end subroutine confll
 
    subroutine congap(xgap)
      implicit none
      real, intent (in) :: xgap
    end subroutine congap
 
    subroutine conlab(cstr)
      implicit none
      character (len = *), intent (in) :: cstr
    end subroutine conlab
 
    subroutine conmat(zmat,n,m,zlev)
      implicit none
      integer, intent (in) :: n,m
      real, dimension (n,m), intent (in) :: zmat
      real, intent (in) :: zlev
    end subroutine conmat
 
    subroutine conmod (xf1,xf2)
      implicit none
      real, intent (in) :: xf1,xf2
    end subroutine conmod
 
    subroutine conn3d(x2,y2,z2)
      implicit none
      real, intent (in) :: x2,y2,z2
    end subroutine conn3d
 
    subroutine connpt(x,y)
      implicit none
      real, intent (in) :: x,y
    end subroutine connpt
 
    subroutine conpts(x,n,y,m,z,zlev,xpts,ypts,maxpts,nray,maxray,nlins)
      implicit none
      integer, intent (in) :: n,m,maxpts,maxray
      real, dimension (n), intent (in) :: x
      real, dimension (m), intent (in) :: y
      real, dimension (n,m), intent (in) :: z
      real, intent (in) :: zlev
      integer, dimension (maxray), intent (out) :: nray
      integer, intent (out) :: nlins
      real, dimension (maxpts), intent (out) :: xpts,ypts
    end subroutine conpts
 
    subroutine conshd(xray,n,yray,m,zmat,zlev,nlray)
      implicit none
      integer, intent (in) :: n,m,nlray
      real, dimension (n), intent (in) :: xray
      real, dimension (m), intent (in) :: yray
      real, dimension (nlray), intent (in) :: zlev
      real, dimension (n,m), intent (in) :: zmat
    end subroutine conshd
 
    subroutine contri(xray,yray,zray,n,i1ray,i2ray,i3ray,ntri,zlev)
      implicit none
      integer, intent (in) :: n,ntri
      real, dimension (n), intent (in) :: xray,yray,zray
      integer, dimension (ntri), intent (in) :: i1ray,i2ray,i3ray
      real, intent (in) :: zlev
    end subroutine contri

    subroutine contur(x,n,y,m,z,zlev)
      implicit none
      integer, intent (in) :: n,m
      real, dimension (n), intent (in) :: x
      real, dimension (m), intent (in) :: y
      real, dimension (n,m), intent (in) :: z
      real, intent (in) :: zlev
    end subroutine contur
 
    subroutine cross()
    end subroutine cross
 
    subroutine crvmat(zmat,ixdim,iydim,ixpts,iypts)
      implicit none
      integer, intent (in) :: ixdim,iydim,ixpts,iypts
      real, dimension (ixdim,iydim), intent (in) :: zmat
    end subroutine crvmat

    subroutine crvtri(xray,yray,zray,n,i1ray,i2ray,i3ray,ntri)
      implicit none
      integer, intent (in) :: n,ntri
      real, dimension (n), intent (in) :: xray,yray,zray
      integer, dimension (ntri), intent (in) :: i1ray,i2ray,i3ray
    end subroutine crvtri

    subroutine curv3d(x,y,z,n)
      implicit none
      integer, intent (in) :: n
      real, dimension (n), intent (in) :: x,y,z
    end subroutine curv3d

    subroutine csrkey(ik)
      implicit none
      integer, intent (out) :: ik
    end subroutine csrkey

    subroutine csrmod(cmod,ckey)
      implicit none
      character (len = *), intent (in) :: cmod,ckey
    end subroutine csrmod

    subroutine csrpos(ix,iy,ik)
      implicit none
      integer, intent (in out) :: ix,iy
      integer, intent (out) :: ik
    end subroutine csrpos

    subroutine csrpt1(ix,iy)
      implicit none
      integer, intent (out) :: ix,iy
    end subroutine csrpt1

    subroutine csrmov(ixray,iyray,nmax,n,iret)
      implicit none
      integer, intent (in) :: nmax
      integer, dimension (nmax), intent (out) :: ixray,iyray
      integer, intent (out) :: n, iret
    end subroutine csrmov

    subroutine csrpts(ixray,iyray,nmax,n,iret)
      implicit none
      integer, intent (in) :: nmax
      integer, dimension (nmax), intent (out) :: ixray,iyray
      integer, intent (out) :: n, iret
    end subroutine csrpts

    subroutine csrrec(ix1,iy1,ix2,iy2)
      implicit none
      integer, intent (out) :: ix1,iy1,ix2,iy2
    end subroutine csrrec

    subroutine csrtyp(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine csrtyp

    subroutine csruni(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine csruni

    subroutine curve(x,y,n)
      implicit none
      integer, intent (in) :: n
      real, dimension (n), intent (in) :: x,y
    end subroutine curve
 
    subroutine curve3(x,y,z,n)
      implicit none
      integer, intent (in) :: n
      real, dimension (n), intent (in) :: x,y,z
    end subroutine curve3
 
    subroutine curvmp(x,y,n)
      implicit none
      integer, intent (in) :: n
      real, dimension (n), intent (in) :: x,y
    end subroutine curvmp
 
    subroutine curvx3(x,y,z,n)
      implicit none
      integer, intent (in) :: n
      real, dimension (n), intent (in) :: x,z
      real, intent (in) :: y
    end subroutine curvx3
 
    subroutine curvy3(x,y,z,n)
      implicit none
      integer, intent (in) :: n
      real, dimension (n), intent (in) :: y,z
      real, intent (in) :: x
    end subroutine curvy3
 
    subroutine cyli3d(x,y,z,r,h,nsk1,nsk2)
      implicit none
      real, intent (in) :: x,y,z,r,h
      integer, intent (in) :: nsk1,nsk2
    end subroutine cyli3d

    subroutine dash()
    end subroutine dash
 
    subroutine dashl()
    end subroutine dashl
 
    subroutine dashm()
    end subroutine dashm
 
    subroutine dattim(cdat,ctim)
      implicit none
      character (len = *), intent (out) :: cdat,ctim
    end subroutine dattim

    subroutine dbffin()
    end subroutine dbffin
 
    subroutine dbfini(iret)
      implicit none
      integer, intent (out) :: iret
    end subroutine dbfini
 
    subroutine delglb()
      implicit none
    end subroutine delglb
 
    subroutine digits(i,cax)
      implicit none
      integer, intent (in) :: i
      character (len = *), intent (in) :: cax
    end subroutine digits
 
    subroutine disalf()
    end subroutine disalf
 
    subroutine disfin()
    end subroutine disfin
 
    subroutine disini()
    end subroutine disini
 
    subroutine disk3d(x,y,z,r1,r2,nsk1,nsk2)
      implicit none
      real, intent (in) :: x,y,z,r1,r2
      integer, intent (in) :: nsk1,nsk2
    end subroutine disk3d

    subroutine doevnt()
    end subroutine doevnt

    subroutine dot()
    end subroutine dot
 
    subroutine dotl()
    end subroutine dotl
 
    subroutine duplx()
    end subroutine duplx
 
    subroutine dwgbut(cstr,ival)
      implicit none
      character (len=*), intent (in) :: cstr
      integer, intent (in out) :: ival
    end subroutine dwgbut
 
    subroutine dwgfil(clab,cstr,cmask)
      implicit none
      character (len=*), intent (in) :: clab,cmask
      character (len=*), intent (in out) :: cstr
    end subroutine dwgfil
 
    subroutine dwglis(clab,clis,ilis)
      implicit none
      character (len=*), intent(in) :: clab,clis
      integer, intent (in out) :: ilis
    end subroutine dwglis
 
    subroutine dwgmsg(cstr)
      implicit none
      character (len=*), intent (in) :: cstr
    end subroutine dwgmsg
 
    subroutine dwgtxt(clab,cstr)
      implicit none
      character (len=*), intent(in) :: clab
      character (len=*), intent(in out) :: cstr
    end subroutine dwgtxt
 
    subroutine ellips(nx,ny,na,nb)
      implicit none
      integer, intent (in) :: nx,ny,na,nb
    end subroutine ellips
 
    subroutine endgrf()
    end subroutine endgrf
 
    subroutine erase()
    end subroutine erase
 
    subroutine errbar(x,y,err1,err2,n)
      implicit none
      integer, intent (in) :: n
      real, dimension (n), intent (in) :: x,y,err1,err2
    end subroutine errbar
 
    subroutine errdev(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine errdev
 
    subroutine errfil(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine errfil

    subroutine errmod(cstr,cmode)
      implicit none
      character (len = *), intent (in) :: cstr,cmode
    end subroutine errmod
 
    subroutine eushft(calph,csft)
      implicit none
      character (len = *), intent (in) :: calph,csft
    end subroutine eushft
 
    subroutine expzlb(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine expzlb
 
    subroutine fcha(x,ndez,nl,cstr)
      implicit none
      real, intent (in) :: x
      integer, intent (in) :: ndez
      integer, intent (out) :: nl
      character (len = *), intent (out) :: cstr
    end subroutine fcha
 
    subroutine field(xray,yray,uray,vray,n,ivec)
      implicit none
      integer, intent (in) :: n,ivec
      real, dimension (n), intent (in) :: xray,yray,uray,vray
    end subroutine field

    subroutine field3d(x1ray,y1ray,z1ray,x2ray,y2ray,z2ray,n,ivec)
      implicit none
      integer, intent (in) :: n,ivec
      real, dimension (n), intent (in) :: x1ray,y1ray,z1ray,x2ray,y2ray,z2ray
    end subroutine field3d
 
    subroutine filbox(nx,ny,nw,nh)
      implicit none
      integer, intent (in) :: nx,ny,nw,nh
    end subroutine filbox
 
    subroutine filclr(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine filclr
 
    subroutine filmod(cmod)
      implicit none
      character (len = *), intent (in) :: cmod
    end subroutine filmod
 
    subroutine filopt(copt,ckey)
      implicit none
      character (len = *), intent (in) :: copt,ckey
    end subroutine filopt

    subroutine fixspc(x)
      implicit none
      real, intent (in) :: x
    end subroutine fixspc
 
    subroutine flab3d()
    end subroutine flab3d
 
    subroutine flen(x,ndez,nx)
      implicit none
      real, intent (in) :: x
      integer, intent (in) :: ndez
      integer, intent (out) :: nx
    end subroutine flen
 
    subroutine frame(i)
      implicit none
      integer, intent (in):: i
    end subroutine frame

    subroutine frmclr(i)
      implicit none
      integer, intent (in):: i
    end subroutine frmclr
 
    subroutine frmess(i)
      implicit none
      integer, intent (in):: i
    end subroutine frmess
 
    subroutine gapcrv(x)
      implicit none
      real, intent (in) :: x
    end subroutine gapcrv

    subroutine gaxpar(a1,a2,copt,cax,a,b,or,stp,ndig)
      implicit none
      real, intent (in) :: a1,a2
      real, intent (out) :: a,b,or,stp
      integer, intent (out) :: ndig
      character (len=*), intent (in) :: copt, cax
    end subroutine gaxpar 
      
    subroutine getalf(cstr)
      implicit none
      character (len = *), intent (out) :: cstr
    end subroutine getalf
 
    subroutine getang(i)
      implicit none
      integer, intent (out) :: i
    end subroutine getang
 
    subroutine getbpp(i)
      implicit none
      integer, intent (out) :: i
    end subroutine getbpp
 
    subroutine getclp(nx,ny,nw,nh)
      implicit none
      integer, intent (out) :: nx,ny,nw,nh
    end subroutine getclp
 
    subroutine getclr(i)
      implicit none
      integer, intent (out) :: i
    end subroutine getclr
 
    subroutine getdig(i,j,k)
      implicit none
      integer, intent (out) :: i,j,k
    end subroutine getdig
 
    subroutine getdsp(cdsp)
      implicit none
      character (len = *), intent (out) :: cdsp
    end subroutine getdsp
 
    subroutine getfil(cstr)
      implicit none
      character (len = *), intent (out) :: cstr
    end subroutine getfil
 
    subroutine getgrf(a,e,or,step,copt)
      implicit none
      real, intent (out) :: a,e,or,step
      character (len = *), intent (in) :: copt
    end subroutine getgrf
 
    subroutine gethgt(i)
      implicit none
      integer, intent (out) :: i
    end subroutine gethgt
 
    subroutine gethnm(i)
      implicit none
      integer, intent (out) :: i
    end subroutine gethnm

    subroutine getind(i,xr,xg,xb)
      implicit none
      integer, intent (in) :: i
      real, intent (out) :: xr,xg,xb
    end subroutine getind
 
    subroutine getlab(c1,c2,c3)
      implicit none
      character (len = *), intent (out) :: c1,c2,c3
    end subroutine getlab
 
    subroutine getlen(i,j,k)
      implicit none
      integer, intent (out) :: i,j,k
    end subroutine getlen
 
    subroutine getlev(i)
      implicit none
      integer, intent (out) :: i
    end subroutine getlev
 
    subroutine getlin(i)
      implicit none
      integer, intent (out) :: i
    end subroutine getlin

    subroutine getlit(xp,yp,zp,xn,yn,zn,i)
      implicit none
      integer, intent (in) :: xp,yp,zp,xn,yn,zn
      integer, intent (out) :: i
    end subroutine getlit
 
    subroutine getmat(x,y,z,n,zmat,nx,ny,zval,imat,wmat)
      implicit none
      integer, intent (in) :: n,nx,ny
      real, dimension (n), intent (in) :: x,y,z
      real, dimension (nx,ny), intent (out) :: zmat
      real, intent (in) :: zval
      integer, dimension (nx,ny), intent (in out) :: imat
      real, dimension (nx,ny), intent (in out) :: wmat
    end subroutine getmat
 
    subroutine getmfl(cstr)
      implicit none
      character (len = *), intent (out) :: cstr
    end subroutine getmfl
 
    subroutine getmix(c,cstr)
      implicit none
      character (len = *), intent (out) :: c
      character (len = *), intent (in) :: cstr
    end subroutine getmix
 
    subroutine getor(i,j)
      implicit none
      integer, intent (out) :: i,j
    end subroutine getor
 
    subroutine getpag(i,j)
      implicit none
      integer, intent (out) :: i,j
    end subroutine getpag
 
    subroutine getpat(i)
      implicit none
      integer, intent (out) :: i
    end subroutine getpat
 
    subroutine getplv(i)
      implicit none
      integer, intent (out) :: i
    end subroutine getplv
 
    subroutine getpos(i,j)
      implicit none
      integer, intent (out) :: i,j
    end subroutine getpos
 
    subroutine getran(i,j)
      implicit none
      integer, intent (out) :: i,j
    end subroutine getran
 
    subroutine getres(i,j)
      implicit none
      integer, intent (out) :: i,j
    end subroutine getres
 
    subroutine getrgb(xr,xg,xb)
      implicit none
      real, intent (out) :: xr,xg,xb
    end subroutine getrgb
 
    subroutine getscl(i,j,k)
      implicit none
      integer, intent (out) :: i,j,k
    end subroutine getscl

    subroutine getscr(i,j)
      implicit none
      integer, intent (out) :: i,j
    end subroutine getscr
 
    subroutine getshf(cstr,c)
      implicit none
      character (len = *), intent (out) :: c
      character (len = *), intent (in) :: cstr
    end subroutine getshf
 
    subroutine getsp1(i,j,k)
      implicit none
      integer, intent (out) :: i,j,k
    end subroutine getsp1
 
    subroutine getsp2(i,j,k)
      implicit none
      integer, intent (out) :: i,j,k
    end subroutine getsp2
 
    subroutine getsym(i,j)
      implicit none
      integer, intent (out) :: i,j
    end subroutine getsym
 
    subroutine gettcl(i,j)
      implicit none
      integer, intent (out) :: i,j
    end subroutine gettcl
 
    subroutine gettic(i,j,k)
      implicit none
      integer, intent (out) :: i,j,k
    end subroutine gettic
 
    subroutine gettyp(i)
      implicit none
      integer, intent (out) :: i
    end subroutine gettyp
 
    subroutine getuni(i)
      implicit none
      integer, intent (out) :: i
    end subroutine getuni
 
    subroutine getver(xver)
      implicit none
      real, intent (out) :: xver
    end subroutine getver
 
    subroutine getvk(i,j,k)
      implicit none
      integer, intent (out) :: i,j,k
    end subroutine getvk
 
    subroutine getvlt(ctab)
      implicit none
      character (len = *), intent (out) :: ctab
    end subroutine getvlt
 
    subroutine getwid(i)
      implicit none
      integer, intent (out) :: i
    end subroutine getwid
 
    subroutine getwin(ix,iy,nw,nh)
      implicit none
      integer, intent (out) :: ix,iy,nw,nh
    end subroutine getwin
 
    subroutine getxid (ival, copt)
      implicit none
      integer, intent (out) :: ival
      character (len = *), intent (in) :: copt
    end subroutine getxid
 
    subroutine gifmod(cmod,ckey)
      implicit none
      character (len = *), intent (in) :: cmod,ckey
    end subroutine gifmod

    subroutine gmxalf(calph,ca,cb,n)
      implicit none
      character (len = *), intent (in) :: calph
      character (len = *), intent (out) :: ca,cb
      integer, intent (out) :: n
    end subroutine gmxalf
 
    subroutine gothic()
    end subroutine gothic
 
    subroutine grace(i)
      implicit none
      integer, intent (in) :: i
    end subroutine grace
 
    subroutine graf(ax,ex,orx,stepx,ay,ey,ory,stepy)
      implicit none
      real, intent (in) :: ax,ex,orx,stepx,ay,ey,ory,stepy
    end subroutine graf
 
    subroutine graf3(ax,ex,orx,stepx,ay,ey,ory,stepy,az,ez,orz,stepz)
      implicit none
      real, intent (in) :: ax,ex,orx,stepx,ay,ey,ory,stepy,az,ez,orz,stepz
    end subroutine graf3
 
    subroutine graf3d(ax,ex,orx,stepx,ay,ey,ory,stepy,az,ez,orz,stepz)
      implicit none
      real, intent (in) :: ax,ex,orx,stepx,ay,ey,ory,stepy,az,ez,orz,stepz
    end subroutine graf3d
 
    subroutine grafmp(ax,ex,orx,stepx,ay,ey,ory,stepy)
      implicit none
      real, intent (in) :: ax,ex,orx,stepx,ay,ey,ory,stepy
    end subroutine grafmp

    subroutine grafp(ex,orx,stepx,ory,stepy)
      implicit none
      real, intent (in) :: ex,orx,stepx,ory,stepy
    end subroutine grafp
 
    subroutine grdpol(igrd,jgrd)
      implicit none
      integer, intent (in) :: igrd,jgrd
    end subroutine grdpol
 
    subroutine grffin()
      implicit none
    end subroutine grffin
 
    subroutine grfini(x1,y1,z1,x2,y2,z2,x3,y3,z3)
      implicit none
      real, intent (in) :: x1,y1,z1,x2,y2,z2,x3,y3,z3
    end subroutine grfini
 
    subroutine grid(i,j)
      implicit none
      integer, intent (in) :: i,j
    end subroutine grid
 
    subroutine grid3d(igrid,jgrid,copt)
      implicit none
      integer, intent (in) :: igrid,jgrid
      character (len = *), intent (in) :: copt
    end subroutine grid3d
 
    subroutine gridmp(i,j)
      implicit none
      integer, intent (in) :: i,j
    end subroutine gridmp
 
    subroutine gwgatt (id,ival,copt)
      implicit none
      integer, intent (in) :: id
      integer, intent (out) :: ival
      character (len=*), intent (in) :: copt
    end subroutine gwgatt
 
    subroutine gwgbox(id,ival)
      implicit none
      integer, intent (in) :: id
      integer, intent (out) :: ival
    end subroutine gwgbox
 
    subroutine gwgbut(id,ival)
      implicit none
      integer, intent (in) :: id
      integer, intent (out) :: ival
    end subroutine gwgbut
 
    subroutine gwgfil(id,cstr)
      implicit none
      integer, intent (in) :: id
      character (len=*), intent (out) :: cstr
    end subroutine gwgfil

    subroutine gwgflt(id,xv)
      implicit none
      integer, intent (in) :: id
      real, intent (out) :: xv
    end subroutine gwgflt

    subroutine gwgint(id,iv)
      implicit none
      integer, intent (in) :: id
      integer, intent (out) :: iv
    end subroutine gwgint
 
    subroutine gwglis(id,ival)
      implicit none
      integer, intent (in) :: id
      integer, intent (out) :: ival
    end subroutine gwglis
 
    subroutine gwgscl(id,xval)
      implicit none
      integer, intent (in) :: id
      real, intent (out) :: xval
    end subroutine gwgscl

    subroutine gwgtbf(id,i,j,xv)
      implicit none
      integer, intent (in) :: id,i,j
      real, intent (out) :: xv
    end subroutine gwgtbf

    subroutine gwgtbi(id,i,j,iv)
      implicit none
      integer, intent (in) :: id,i,j
      integer, intent (out) :: iv
    end subroutine gwgtbi

    subroutine gwgtbl(id,xray,n,idx,copt)
      implicit none
      integer, intent (in) :: id,n,idx
      real, dimension (n), intent (out) :: xray
      character (len=*), intent (in) :: copt
    end subroutine gwgtbl

    subroutine gwgtbs(id,i,j,cstr)
      implicit none
      integer, intent (in) :: id,i,j
      character (len=*), intent (out) :: cstr
    end subroutine gwgtbs
 
    subroutine gwgtxt(id,cstr)
      implicit none
      integer, intent (in) :: id
      character (len=*), intent (out) :: cstr
    end subroutine gwgtxt

    subroutine gwgxid(id,ival)
      implicit none
      integer, intent (in) :: id
      integer, intent (out) :: ival
    end subroutine gwgxid
 
    subroutine height(i)
      implicit none
      integer, intent (in) :: i
    end subroutine height
 
    subroutine helve()
    end subroutine helve
 
    subroutine helves()
    end subroutine helves
 
    subroutine helvet()
    end subroutine helvet
 
    subroutine histog(xray,n,x,y,m)
      implicit none
      integer, intent (in) :: n
      real, dimension (n), intent (in) :: xray
      integer, intent (out) :: m
      real, dimension (n), intent (out) :: x,y
    end subroutine histog
 
    subroutine hname(i)
      implicit none
      integer, intent (in) :: i
    end subroutine hname
 
    subroutine hpgmod(cmod,ckey)
      implicit none
      character (len = *), intent (in) :: cmod,ckey
    end subroutine hpgmod

    subroutine hsvrgb(xh,xs,xv,r,g,b)
      implicit none
      real, intent (in)  :: xh,xs,xv
      real, intent (out) :: r,g,b
    end subroutine hsvrgb

    subroutine hsym3d(xh)
      implicit none
      real, intent (in) :: xh
    end subroutine hsym3d
 
    subroutine hsymbl(i)
      implicit none
      integer, intent (in) :: i
    end subroutine hsymbl
 
    subroutine htitle(i)
      implicit none
      integer, intent (in) :: i
    end subroutine htitle
 
    subroutine hwfont()
    end subroutine hwfont
 
    subroutine hworig(nx,ny)
      implicit none
      integer, intent (in) :: nx,ny
    end subroutine hworig
 
    subroutine hwpage(nxp,nyp)
      implicit none
      integer, intent (in) :: nxp,nyp
    end subroutine hwpage

    subroutine hwscal(x)
      implicit none
      real, intent (in) :: x
    end subroutine hwscal

    subroutine imgbox(nx,ny,nw,nh)
      implicit none
      integer, intent (in) :: nx,ny,nw,nh
    end subroutine imgbox

    subroutine imgclp(nx,ny,nw,nh)
      implicit none
      integer, intent (in) :: nx,ny,nw,nh
    end subroutine imgclp

    subroutine imgfin()
    end subroutine imgfin

    subroutine imgfmt(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine imgfmt
 
    subroutine imgini()
    end subroutine imgini

    subroutine imgmod(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine imgmod

    subroutine imgsiz(nw,nh)
      implicit none
      integer, intent (in) :: nw,nh
    end subroutine imgsiz

    subroutine inccrv(i)
      implicit none
      integer, intent (in) :: i
    end subroutine inccrv

    function incdat(id,im,iy)
      implicit none
      integer, intent (in) :: id, im, iy
      integer :: incdat
    end function incdat
 
    subroutine incfil(cstr)
      implicit none
      character (len = *), intent (in) :: cstr
    end subroutine incfil
 
    subroutine incmrk(i)
      implicit none
      integer, intent (in) :: i
    end subroutine incmrk

    function indrgb(xr,xg,xb)
      implicit none
      real, intent (in) :: xr,xg,xb
      integer :: indrgb 
    end function indrgb
 
    subroutine intax()
    end subroutine intax
 
    subroutine intcha(num,n,cnum)
      implicit none
      integer, intent (in)  :: num
      integer, intent (out) :: n
      character (len = *), intent (out) :: cnum
    end subroutine intcha
 
    subroutine intlen(nm,nlaen)
      implicit none
      integer, intent (in)  :: nm
      integer, intent (out) :: nlaen
    end subroutine intlen
 
    function intrgb(xr,xg,xb)
      implicit none
      real, intent (in) :: xr,xg,xb
      integer :: intrgb 
    end function intrgb

    subroutine intutf(iray,nray,cstr,nmax,nl)
      implicit none
      integer, intent (in) :: nray,nmax
      integer, dimension (nray), intent (in) :: iray
      character (len=*), intent (out) :: cstr
      integer, intent (out) :: nl
    end subroutine intutf

    subroutine isopts(xray,nx,yray,ny,zray,nz,wmat,wlev,  &
               &      xtri,ytri,ztri,nmax,ntri)
      implicit none
      integer, intent (in) :: nx,ny,nz,nmax
      integer, intent (out) :: ntri
      real, dimension (nx), intent (in) :: xray
      real, dimension (ny), intent (in) :: yray
      real, dimension (nz), intent (in) :: zray
      real, dimension (nx,ny,nz), intent (in) :: wmat
      real, intent (in) :: wlev
      real, dimension (nmax), intent (out) :: xtri,ytri,ztri
    end subroutine isopts

    subroutine itmcat(clis,cstr)
      implicit none
      character (len=*), intent (in out) :: clis
      character (len=*), intent (in) :: cstr
    end subroutine itmcat
 
    function itmcnt(clis)
      implicit none
      character (len=*), intent (in) :: clis
      integer :: itmcnt
    end function itmcnt
 
    subroutine itmstr(clis,nlis,cstr)
      implicit none
      character (len=*), intent (in) :: clis
      character (len=*), intent (out) :: cstr
      integer, intent (in) :: nlis
    end subroutine itmstr
 
    subroutine labclr(iclr,copt)
      implicit none
      integer, intent (in) :: iclr
      character (len = *), intent (in) :: copt
    end subroutine labclr
 
    subroutine labdig(i,cax)
      implicit none
      integer, intent (in) :: i
      character (len = *), intent (in) :: cax
    end subroutine labdig
 
    subroutine labdis(i,cax)
      implicit none
      integer, intent (in) :: i
      character (len = *), intent (in) :: cax
    end subroutine labdis
 
    subroutine labels(copt,cax)
      implicit none
      character (len = *), intent (in) :: copt,cax
    end subroutine labels
 
    subroutine labjus(copt,cax)
      implicit none
      character (len = *), intent (in) :: copt,cax
    end subroutine labjus

    subroutine labl3d(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine labl3d

    subroutine labmod(ckey,copt,cax)
      implicit none
      character (len = *), intent (in) :: ckey,copt,cax
    end subroutine labmod
 
    subroutine labpos(copt,cax)
      implicit none
      character (len = *), intent (in) :: copt,cax
    end subroutine labpos
 
    subroutine labtyp(copt,cax)
      implicit none
      character (len = *), intent (in) :: copt,cax
    end subroutine labtyp
 
    subroutine legclr()
    end subroutine legclr
 
    subroutine legend(cbf,ncor)
      implicit none
      character (len = *), intent (in) :: cbf
      integer, intent (in) :: ncor
    end subroutine legend
 
    subroutine legini(cbf,nlin,nmax)
      implicit none
      character (len = *), intent (in out) :: cbf
      integer, intent (in) :: nlin, nmax
    end subroutine legini
 
    subroutine leglin(cbf,cstr,n)
      implicit none
      character (len = *), intent (in out) :: cbf
      character (len = *), intent (in) :: cstr
      integer, intent (in) :: n
    end subroutine leglin
 
    subroutine legopt(x1,x2,x3)
      implicit none
      real, intent (in) :: x1,x2,x3
    end subroutine legopt
 
    subroutine legpat(ilin,ithk,isym,iclr,ipat,i)
      implicit none
      integer, intent (in) :: ilin,ithk,isym,iclr,ipat,i
    end subroutine legpat
 
    subroutine legpos(nx,ny)
      implicit none
      integer, intent (in) :: nx,ny
    end subroutine legpos
 
    subroutine legtit(cstr)
      implicit none
      character (len = *), intent (in) :: cstr
    end subroutine legtit
 
    subroutine legval(x,copt)
      implicit none
      real, intent (in) :: x 
      character (len = *), intent (in) :: copt
    end subroutine legval

    subroutine lfttit()
    end subroutine lfttit
 
    subroutine lincyc(i,ilin)
      implicit none
      integer, intent (in) :: i,ilin
    end subroutine lincyc
 
    subroutine line(nx,ny,nu,nv)
      implicit none
      integer, intent (in) :: nx,ny,nu,nv
    end subroutine line
 
    subroutine linesp(x)
      implicit none
      real, intent (in) :: x
    end subroutine linesp
 
    subroutine lintyp(i)
      implicit none
      integer, intent (in) :: i
    end subroutine lintyp
 
    subroutine linwid(i)
      implicit none
      integer, intent (in) :: i
    end subroutine linwid

    subroutine light(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine light

    subroutine litmod(id,copt)
      implicit none
      integer, intent (in) :: id
      character (len = *), intent (in) :: copt
    end subroutine litmod

    subroutine litopt(id,x,copt)
      implicit none
      integer, intent (in) :: id
      real, intent (in) :: x
      character (len = *), intent (in) :: copt
    end subroutine litopt

    subroutine litop3(id,xr,xg,xb,copt)
      implicit none
      integer, intent (in) :: id
      real, intent (in) :: xr,xg,xb
      character (len = *), intent (in) :: copt
    end subroutine litop3
 
    subroutine litpos(id,xp,yp,zp,copt)
      implicit none
      integer, intent (in) :: id
      real, intent (in) :: xp,yp,zp
      character (len = *), intent (in) :: copt
    end subroutine litpos

    subroutine lncap(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine lncap
 
    subroutine lnjoin(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine lnjoin
 
    subroutine lnmlt(x)
      implicit none
      real, intent (in) :: x
    end subroutine lnmlt
 
    subroutine logtic(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine logtic
 
    subroutine lsechk(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine lsechk
 
    subroutine mapbas(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine mapbas

    subroutine mapfil(cfl, copt)
      implicit none
      character (len = *), intent (in) :: cfl, copt
    end subroutine mapfil

    subroutine maplab(copt,ckey)
      implicit none
      character (len = *), intent (in) :: copt,ckey
    end subroutine maplab
 
    subroutine maplev(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine maplev
 
    subroutine mapmod(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine mapmod
 
    subroutine mappol(x,y)
      implicit none
      real, intent (in) :: x,y
    end subroutine mappol

    subroutine mapopt(copt,ckey)
      implicit none
      character (len = *), intent (in) :: copt,ckey
    end subroutine mapopt
 
    subroutine mapref(ylw,yup)
      implicit none
      real, intent (in) :: ylw,yup
    end subroutine mapref
 
    subroutine mapsph(xrad)
      implicit none
      real, intent (in) :: xrad
    end subroutine mapsph
 
    subroutine marker(i)
      implicit none
      integer, intent (in) :: i
    end subroutine marker

    subroutine matopt(x,copt)
      implicit none
      real, intent (in) :: x
      character (len = *), intent (in) :: copt
    end subroutine matopt

    subroutine matop3(xr,xg,xb,copt)
      implicit none
      real, intent (in) :: xr,xg,xb
      character (len = *), intent (in) :: copt
    end subroutine matop3
 
    subroutine mdfmat(i,j,x)
      implicit none
      integer, intent (in) :: i,j
      real, intent (in) :: x
    end subroutine mdfmat
 
    subroutine messag(cstr,nx,ny)
      implicit none
      character (len = *), intent (in) :: cstr
      integer, intent (in) :: nx,ny
    end subroutine messag
 
    subroutine metafl(ct)
      implicit none
      character (len = *), intent (in) :: ct
    end subroutine metafl
 
    subroutine mixalf()
    end subroutine mixalf
 
    subroutine mixleg()
    end subroutine mixleg
 
    subroutine moment(xray,n,copt,xv)
      implicit none
      integer, intent (in) :: n
      real, dimension (n), intent (in) :: xray
      character (len = *), intent (in) :: copt
      real, intent (out) :: xv
    end subroutine moment
 
    subroutine mpaepl(i)
      implicit none
      integer, intent (in) :: i
    end subroutine mpaepl
 
    subroutine mplang(x)
      implicit none
      real, intent (in) :: x
    end subroutine mplang
 
    subroutine mplclr(i,j)
      implicit none
      integer, intent (in) :: i,j
    end subroutine mplclr
 
    subroutine mplpos(i,j)
      implicit none
      integer, intent (in) :: i,j
    end subroutine mplpos
 
    subroutine mplsiz(i)
      implicit none
      integer, intent (in) :: i
    end subroutine mplsiz
 
    subroutine mpslogo(nx,ny,nsize,copt)
      implicit none
      integer, intent (in) :: nx,ny,nsize
      character (len = *), intent (in) :: copt
    end subroutine mpslogo

    subroutine msgbox(cstr)
      implicit none
      character (len = *), intent (in) :: cstr
    end subroutine msgbox
 
    subroutine mshclr(ic)
      implicit none
      integer, intent (in) :: ic
    end subroutine mshclr

    subroutine mshcrv(n)
      implicit none
      integer, intent (in) :: n
    end subroutine mshcrv

    subroutine mylab(cstr,i,cax)
      implicit none
      character (len = *), intent (in) :: cstr,cax
      integer, intent (in) :: i
    end subroutine mylab
 
    subroutine myline(nray,n)
      implicit none
      integer, intent (in) :: n
      integer, dimension (n), intent (in) :: nray
    end subroutine myline
 
    subroutine mypat(iang,itype,idens,icross)
      implicit none
      integer, intent (in) :: iang,itype,idens,icross
    end subroutine mypat
 
    subroutine mypie(iseg,xdat,xper,nrad,noff,ang,nvbox,idrw,iann)
      implicit none
      integer, intent (in out) :: iseg,nrad,noff,nvbox,idrw,iann
      real, intent (in out)    :: xdat,xper,ang
    end subroutine mypie

    subroutine mysymb(xray,yray,n,isym,iflag)
      implicit none
      integer, intent (in) :: n,isym,iflag
      real, dimension (n), intent (in) :: xray,yray
    end subroutine mysymb
 
    subroutine myvlt(xr,xg,xb,n)
      implicit none
      integer, intent (in) :: n
      real, dimension (n), intent (in) :: xr,xg,xb
    end subroutine myvlt
 
    subroutine namdis(i,cax)
      implicit none
      integer, intent (in) :: i
      character (len = *), intent (in) :: cax
    end subroutine namdis
 
    subroutine name(cnam,cax)
      implicit none
      character (len = *), intent (in) :: cnam,cax
    end subroutine name
 
    subroutine namjus(copt,cax)
      implicit none
      character (len = *), intent (in) :: copt,cax
    end subroutine namjus
 
    subroutine neglog(e)
      implicit none
      real, intent (in) :: e
    end subroutine neglog
 
    subroutine newmix()
    end subroutine newmix
 
    subroutine newpag()
    end subroutine newpag
 
    function nlmess(cstr)
      implicit none
      character (len = *), intent (in) :: cstr
      integer :: nlmess
    end function nlmess
 
    function nlnumb(x,ndez)
      implicit none
      real, intent (in) :: x
      integer, intent (in) :: ndez
      integer :: nlnumb
    end function nlnumb
 
    subroutine noarln()
    end subroutine noarln
 
    subroutine nobar()
    end subroutine nobar
 
    subroutine nobgd()
    end subroutine nobgd
 
    subroutine nochek()
    end subroutine nochek
 
    subroutine noclip()
    end subroutine noclip
 
    subroutine nofill()
    end subroutine nofill
 
    subroutine nograf()
    end subroutine nograf
 
    subroutine nohide()
    end subroutine nohide
 
    subroutine noline(cax)
      implicit none
      character (len = *), intent (in) :: cax
    end subroutine noline
 
    subroutine number(x,ndez,nx,ny)
      implicit none
      real, intent (in) :: x
      integer, intent (in) :: ndez,nx,ny
    end subroutine number
 
    subroutine numfmt(ct)
      implicit none
      character (len = *), intent (in) :: ct
    end subroutine numfmt
 
    subroutine numode(cdec,cgrp,cpos,cspc)
      implicit none
      character (len = *), intent (in) :: cdec,cgrp,cpos,cspc
    end subroutine numode
 
    function nwkday(id,im,iy)
      implicit none
      integer, intent (in) :: id,im,iy
      integer :: nwkday
    end function nwkday
     
    function nxlegn(cbf)
      implicit none
      character (len = *), intent (in) :: cbf
      integer :: nxlegn
    end function nxlegn

    function nxpixl(ix,iy)
      implicit none
      integer, intent (in) :: ix,iy
      integer :: nxpixl
    end function nxpixl
 
    function nxposn(x)
      implicit none
      real, intent (in) :: x
      integer :: nxposn
    end function nxposn
 
    function nylegn(cbf)
      implicit none
      character (len = *), intent (in) :: cbf
      integer :: nylegn
    end function nylegn

    function nypixl(ix,iy)
      implicit none
      integer, intent (in) :: ix,iy
      integer :: nypixl
    end function nypixl
 
    function nyposn(y)
      implicit none
      real, intent (in) :: y
      integer :: nyposn
    end function nyposn
 
    function nzposn(z)
      implicit none
      real, intent (in) :: z
      integer :: nzposn
    end function nzposn
 
    subroutine openfl(cstr,nlu,irw,istat)
      implicit none
      character (len = *), intent (in) :: cstr
      integer, intent (in) :: nlu,irw
      integer, intent (out) :: istat
    end subroutine openfl
 
    subroutine opnwin(id)
      implicit none
      integer, intent (in) :: id
    end subroutine opnwin
 
    subroutine origin(i,j)
      implicit none
      integer, intent (in) :: i,j
    end subroutine origin
 
    subroutine page(nxp,nyp)
      implicit none
      integer, intent (in) :: nxp,nyp
    end subroutine page
 
    subroutine pagera()
    end subroutine pagera
 
    subroutine pagfll(n)
      implicit none
      integer, intent (in) :: n
    end subroutine pagfll
 
    subroutine paghdr(c1,c2,iopt,idir)
      implicit none
      character (len = *), intent (in) :: c1,c2
      integer, intent (in) :: iopt,idir
    end subroutine paghdr
 
    subroutine pagmod(cmod)
      implicit none
      character (len = *), intent (in) :: cmod
    end subroutine pagmod

    subroutine pagorg(cmod)
      implicit none
      character (len = *), intent (in) :: cmod
    end subroutine pagorg

    subroutine pagwin(nxp,nyp)
      implicit none
      integer, intent (in) :: nxp,nyp
    end subroutine pagwin
 
    subroutine patcyc(i,ipat)
      implicit none
      integer, intent (in) :: i,ipat
    end subroutine patcyc

    subroutine pdfbuf(iray,nmax,nn)
      implicit none
      integer, intent (in) :: nmax
      character (len=1), intent (out), dimension (*) :: iray
      integer, intent (out) :: nn
    end subroutine pdfbuf

    subroutine pdfmod(cmod,ckey)
      implicit none
      character (len = *), intent (in) :: cmod,ckey
    end subroutine pdfmod
 
    subroutine pdfmrk(cstr,copt)
      implicit none
      character (len = *), intent (in) :: cstr,copt
    end subroutine pdfmrk

    subroutine penwid(x)
      implicit none
      real, intent (in) :: x
    end subroutine penwid
 
    subroutine pie(nx,ny,nr,a,b)
      implicit none
      integer, intent (in) :: nx,ny,nr
      real, intent (in) :: a,b
    end subroutine pie
 
    subroutine piebor(iclr)
      implicit none
      integer, intent (in) :: iclr
    end subroutine piebor

    subroutine pieclr(ic1,ic2,n)
      implicit none
      integer, intent (in) :: n
      integer, dimension (n), intent (in) :: ic1,ic2
    end subroutine pieclr
 
    subroutine pieexp()
    end subroutine pieexp
 
    subroutine piegrf(cstr,nlin,xray,n)
      implicit none
      character (len = *), intent (in) :: cstr
      integer, intent (in) :: nlin,n
      real, dimension (n), intent (in) :: xray
    end subroutine piegrf
 
    subroutine pielab(cdat,cstr)
      implicit none
      character(len = *), intent (in) :: cdat,cstr
    end subroutine pielab

    subroutine pieopt(x1,x2)
      implicit none
      real, intent (in) :: x1,x2
    end subroutine pieopt

    subroutine pietyp(ctyp)
      implicit none
      character (len = *), intent (in) :: ctyp
    end subroutine pietyp
 
    subroutine pievec(ivec,copt)
      implicit none
      integer, intent (in) :: ivec
      character (len = *), intent (in) :: copt
    end subroutine pievec

    subroutine pike3d(x1,y1,z1,x2,y2,z2,r,nsk1,nsk2)
      implicit none
      real, intent (in) :: x1,y1,z1,x2,y2,z2,r
      integer, intent (in) :: nsk1,nsk2
    end subroutine pike3d

    subroutine plat3d(x,y,z,xl,copt)
      implicit none
      real, intent (in) :: x,y,z,xl
      character (len = *), intent (in) :: copt
    end subroutine plat3d

    subroutine pngmod(cmod,ckey)
      implicit none
      character (len = *), intent (in) :: cmod,ckey
    end subroutine pngmod
 
    subroutine point(nx,ny,nb,nh,ncol)
      implicit none
      integer, intent (in) :: nx,ny,nb,nh,ncol
    end subroutine point
 
    subroutine polar(ex,orx,stepx,ory,stepy)
      implicit none
      real, intent (in) :: ex,orx,stepx,ory,stepy
    end subroutine polar

    subroutine polclp(xray,yray,n,xout,yout,nmax,nout,xv,cedge)
      implicit none
      integer, intent (in) :: n,nmax
      integer, intent (out) :: nout
      real, dimension (n), intent (in) :: xray,yray
      real, dimension (nmax), intent (out) :: xout,yout
      real, intent (in) :: xv
      character (len = *), intent (in) :: cedge
    end subroutine polclp 

    subroutine polcrv(cstr)
      implicit none
      character (len = *), intent (in) :: cstr
    end subroutine polcrv

    subroutine polmod(cpos,cdir)
      implicit none
      character (len = *), intent (in) :: cpos,cdir
    end subroutine polmod
 
    subroutine pos2pt(x,y,xp,yp)
      implicit none
      real, intent (in) :: x,y
      real, intent (out) :: xp,yp
    end subroutine pos2pt
 
    subroutine pos3pt(x,y,z,xp,yp,zp)
      implicit none
      real, intent (in) :: x,y,z
      real, intent (out) :: xp,yp,zp
    end subroutine pos3pt
 
    subroutine posifl(nlu,nbyt,istat)
      implicit none
      integer, intent (in) :: nlu,nbyt
      integer, intent (out) :: istat
    end subroutine posifl
 
    subroutine projct(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine projct
 
    subroutine psfont(cstr)
      implicit none
      character (len = *), intent (in) :: cstr
    end subroutine psfont

    subroutine psmode(cstr)
      implicit none
      character (len = *), intent (in) :: cstr
    end subroutine psmode

    subroutine pyra3d(x,y,z,xl,h1,h2,n)
      implicit none
      real, intent (in) :: x,y,z,xl,h1,h2
      integer, intent (in) :: n
    end subroutine pyra3d
 
    subroutine qplbar(x,n)
      implicit none
      integer, intent (in) :: n
      real, dimension (n), intent (in) :: x
    end subroutine qplbar

    subroutine qplclr(x,n,m)
      implicit none
      integer, intent (in) :: n,m
      real, dimension (n,m), intent (in) :: x
    end subroutine qplclr

    subroutine qplcon(x,n,m,nlev)
      implicit none
      integer, intent (in) :: n,m,nlev
      real, dimension (n,m), intent (in) :: x
    end subroutine qplcon

    subroutine qplot(x,y,n)
      implicit none
      integer, intent (in) :: n
      real, dimension (n), intent (in) :: x,y
    end subroutine qplot

    subroutine qplpie(x,n)
      implicit none
      integer, intent (in) :: n
      real, dimension (n), intent (in) :: x
    end subroutine qplpie

    subroutine qplsca(x,y,n)
      implicit none
      integer, intent (in) :: n
      real, dimension (n), intent (in) :: x,y
    end subroutine qplsca

    subroutine qplsur(x,n,m)
      implicit none
      integer, intent (in) :: n,m
      real, dimension (n,m), intent (in) :: x
    end subroutine qplsur

    subroutine quad3d(x,y,z,xl,yl,zl)
      implicit none
      real, intent (in) :: x,y,z,xl,yl,zl
    end subroutine quad3d 

    subroutine rbfpng(iray,nmax,nn)
      implicit none
      integer, intent (in) :: nmax
      character (len=1), intent (out), dimension (*) :: iray
      integer, intent (out) :: nn
    end subroutine rbfpng

    subroutine rbmp(cfl)
      implicit none
      character (len=*), intent (in) :: cfl
    end subroutine rbmp
 
    subroutine readfl(nlu,iray,nbyt,istat)
      implicit none
      integer, intent (in) :: nlu,nbyt
      character (len=1), intent (out), dimension (nbyt) :: iray
      integer, intent (out) :: istat
    end subroutine readfl
 
    subroutine reawgt()
    end subroutine reawgt

    subroutine recfll(nx,ny,nb,nh,ncol)
      implicit none
      integer, intent (in) :: nx,ny,nb,nh,ncol
    end subroutine recfll
 
    subroutine rectan(nx,ny,nb,nh)
      implicit none
      integer, intent (in) :: nx,ny,nb,nh
    end subroutine rectan
 
    subroutine rel3pt(x,y,z,xp,yp)
      implicit none
      real, intent (in) :: x,y,z
      real, intent (out) :: xp,yp
    end subroutine rel3pt
 
    subroutine resatt()
    end subroutine resatt
 
    subroutine reset(cw)
      implicit none
      character (len = *), intent (in) :: cw
    end subroutine reset
 
    subroutine revscr()
    end subroutine revscr
 
    subroutine rgbhsv(r,g,b,xh,xs,xv)
      implicit none
      real, intent (in) :: r,g,b
      real, intent (out) :: xh,xs,xv
    end subroutine rgbhsv

    subroutine rgif(cfl)
      implicit none
      character (len=*), intent (in) :: cfl
    end subroutine rgif
 
    subroutine rgtlab()
    end subroutine rgtlab
 
    subroutine rimage(cfl)
      implicit none
      character (len = *), intent (in) :: cfl
    end subroutine rimage
 
    subroutine rlarc(xm,ym,a,b,alpha,beta,theta)
      implicit none
      real, intent (in) :: xm,ym,a,b,alpha,beta,theta
    end subroutine rlarc
 
    subroutine rlarea(x,y,n)
      implicit none
      integer, intent (in) :: n
      real, dimension (n), intent (in) :: x,y
    end subroutine rlarea
 
    subroutine rlcirc(xm,ym,r)
      implicit none
      real, intent (in) :: xm,ym,r
    end subroutine rlcirc
 
    subroutine rlconn(x,y)
      implicit none
      real, intent (in) :: x,y
    end subroutine rlconn
 
    subroutine rlell(xm,ym,a,b)
      implicit none
      real, intent (in) :: xm,ym,a,b
    end subroutine rlell
 
    subroutine rline(x,y,u,v)
      implicit none
      real, intent (in) :: x,y,u,v
    end subroutine rline
 
    subroutine rlmess(cstr,x,y)
      implicit none
      character (len = *), intent (in) :: cstr
      real, intent (in) :: x,y
    end subroutine rlmess
 
    subroutine rlnumb(z,ndez,x,y)
      implicit none
      real, intent (in) :: z,x,y
      integer, intent (in) :: ndez
    end subroutine rlnumb
 
    subroutine rlpie(xm,ym,r,alpha,beta)
      implicit none
      real, intent (in) :: xm,ym,r,alpha,beta
    end subroutine rlpie
 
    subroutine rlpoin(x,y,nb,nh,ncol)
      implicit none
      real, intent (in) :: x,y
      integer, intent (in) :: nb,nh,ncol
    end subroutine rlpoin
 
    subroutine rlrec(x,y,xb,xh)
      implicit none
      real, intent (in) :: x,y,xb,xh
    end subroutine rlrec
 
    subroutine rlrnd(x,y,xb,xh,irnd)
      implicit none
      real, intent (in) :: x,y,xb,xh
      integer, intent (in) :: irnd
    end subroutine rlrnd
 
    subroutine rlsec(xm,ym,r1,r,beta,alpha,ncol)
      implicit none
      real, intent (in) :: xm,ym,r1,r,beta,alpha
      integer, intent (in) :: ncol
    end subroutine rlsec
 
    subroutine rlstrt(x,y)
      implicit none
      real, intent (in) :: x,y
    end subroutine rlstrt
 
    subroutine rlsymb(i,x,y)
      implicit none
      integer, intent (in) :: i
      real, intent (in) :: x,y
    end subroutine rlsymb
 
    subroutine rlvec(x,y,u,v,ivec)
      implicit none
      real, intent (in) :: x,y,u,v
      integer, intent (in) :: ivec
    end subroutine rlvec

    subroutine rlwind(xk,x,y,nw,a)
      implicit none
      real, intent (in) :: xk,x,y,a
      integer, intent (in) :: nw
    end subroutine rlwind
 
    subroutine rndrec(nx,ny,nb,nh,irnd)
      implicit none
      integer, intent (in) :: nx,ny,nb,nh,irnd
    end subroutine rndrec

    subroutine rot3d(xa,ya,za)
      implicit none
      real, intent (in) :: xa,ya,za
    end subroutine rot3d 
 
    subroutine rpixel(ix,iy,n)
      implicit none
      integer, intent (in)  :: ix,iy
      integer, intent (out) :: n
    end subroutine rpixel
 
    subroutine rpixls(iray,ix,iy,nw,nh)
      implicit none
      integer, intent (in) :: ix,iy,nw,nh
      character (len=1), intent (out), dimension (*) :: iray
    end subroutine rpixls

    subroutine rpng(cfl)
      implicit none
      character (len=*), intent (in) :: cfl
    end subroutine rpng

    subroutine rppm(cfl)
      implicit none
      character (len=*), intent (in) :: cfl
    end subroutine rppm
 
    subroutine rpxrow(iray,ix,iy,n)
      implicit none
      integer, intent (in) :: ix,iy,n
      character (len=1), intent (out), dimension (n) :: iray
    end subroutine rpxrow
 
    subroutine rtiff(cfl)
      implicit none
      character (len=*), intent (in) :: cfl
    end subroutine rtiff
 
    subroutine rvynam()
    end subroutine rvynam
 
    subroutine scale(copt,cax)
      implicit none
      character (len = *), intent (in) :: copt,cax
    end subroutine scale
 
    subroutine sclfac(x)
      implicit none
      real, intent (in) :: x
    end subroutine sclfac
 
    subroutine sclmod(cmod)
      implicit none
      character (len = *), intent (in) :: cmod
    end subroutine sclmod
 
    subroutine scmplx()
    end subroutine scmplx
 
    subroutine scrmod(cmod)
      implicit none
      character (len = *), intent (in) :: cmod
    end subroutine scrmod
 
    subroutine sector(nx,ny,nr1,nr2,alpha,beta,ncol)
      implicit none
      integer, intent (in) :: nx,ny,nr1,nr2,ncol
      real, intent (in) :: alpha,beta
    end subroutine sector
 
    subroutine selwin(id)
      implicit none
      integer, intent (in) :: id
    end subroutine selwin
 
    subroutine sendbf()
    end subroutine sendbf

    subroutine sendmb()
    end subroutine sendmb
 
    subroutine sendok()
    end subroutine sendok
 
    subroutine serif()
    end subroutine serif
 
    subroutine setbas(f)
      implicit none
      real, intent (in) :: f
    end subroutine setbas
 
    subroutine setcbk (callbk, copt)
      implicit none
      character (len = *), intent (in) :: copt
 
      interface
         subroutine callbk (xp, yp)
           implicit none
           real, intent (in out) :: xp,yp
         end subroutine callbk
      end interface
    end subroutine setcbk

    subroutine setclr(n)
      implicit none
      integer, intent (in) :: n
    end subroutine setclr
 
    subroutine setcsr(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine setcsr

    subroutine setexp(f)
      implicit none
      real, intent (in) :: f
    end subroutine setexp

    subroutine setfce(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine setfce
 
    subroutine setfil(ct)
      implicit none
      character (len = *), intent (in) :: ct
    end subroutine setfil
 
    subroutine setgrf(c1,c2,c3,c4)
      implicit none
      character (len = *), intent (in) :: c1,c2,c3,c4
    end subroutine setgrf
 
    subroutine setind(i,xr,xg,xb)
      implicit none
      integer, intent (in) :: i
      real, intent (in) :: xr,xg,xb
    end subroutine setind
 
    subroutine setmix(c,cstr)
      implicit none
      character (len = *), intent (in) :: c,cstr
    end subroutine setmix
 
    subroutine setpag(ct)
      implicit none
      character (len = *), intent (in) :: ct
    end subroutine setpag
 
    subroutine setres(i,j)
      implicit none
      integer, intent (in) :: i,j
    end subroutine setres
 
    subroutine setrgb(xr,xg,xb)
      implicit none
      real, intent (in) :: xr,xg,xb
    end subroutine setrgb
 
    subroutine setscl(xray,n,cstr)
      implicit none
      integer, intent (in) :: n
      real, dimension (n), intent (in) :: xray
      character (len = *), intent (in) :: cstr
    end subroutine setscl
 
    subroutine setvlt(ctab)
      implicit none
      character (len = *), intent (in) :: ctab
    end subroutine setvlt
 
    subroutine setxid(i,copt)
      implicit none
      integer, intent (in) :: i
      character (len = *), intent (in) :: copt
    end subroutine setxid

    subroutine shdafr(inat,ishd,iclr,n)
      implicit none
      integer, intent (in) :: n
      integer, dimension (n), intent (in) :: inat,ishd,iclr
    end subroutine shdafr

    subroutine shdasi(inat,ishd,iclr,n)
      implicit none
      integer, intent (in) :: n
      integer, dimension (n), intent (in) :: inat,ishd,iclr
    end subroutine shdasi

    subroutine shdaus(inat,ishd,iclr,n)
      implicit none
      integer, intent (in) :: n
      integer, dimension (n), intent (in) :: inat,ishd,iclr
    end subroutine shdaus
 
    subroutine shdcha()
    end subroutine shdcha
 
    subroutine shdcrv(x1,y1,n1,x2,y2,n2)
      implicit none
      integer, intent (in) :: n1,n2
      real, dimension (n1), intent (in) :: x1,y1
      real, dimension (n2), intent (in) :: x2,y2
    end subroutine shdcrv
 
    subroutine shdeur(inat,ishd,iclr,n)
      implicit none
      integer, intent (in) :: n
      integer, dimension (n), intent (in) :: inat,ishd,iclr
    end subroutine shdeur
 
    subroutine shdmap(cmap)
      implicit none
      character (len = *), intent (in) :: cmap
    end subroutine shdmap
 
    subroutine shdmod(copt,ctype)
      implicit none
      character (len = *), intent (in) :: copt,ctype
    end subroutine shdmod
 
    subroutine shdnor(inat,ishd,iclr,n)
      implicit none
      integer, intent (in) :: n
      integer, dimension (n), intent (in) :: inat,ishd,iclr
    end subroutine shdnor

    subroutine shdpat(i)
      implicit none
      integer, intent (in) :: i
    end subroutine shdpat

    subroutine shdsou(inat,ishd,iclr,n)
      implicit none
      integer, intent (in) :: n
      integer, dimension (n), intent (in) :: inat,ishd,iclr
    end subroutine shdsou

    subroutine shdusa(inat,ishd,iclr,n)
      implicit none
      integer, intent (in) :: n
      integer, dimension (n), intent (in) :: inat,ishd,iclr
    end subroutine shdusa
 
    subroutine shield(cblnk,cmode)
      implicit none
      character (len = *), intent (in) :: cblnk,cmode
    end subroutine shield
 
    subroutine shlcir(nx,ny,nr)
      implicit none
      integer, intent (in) :: nx,ny,nr
    end subroutine shlcir
 
    subroutine shldel(id)
      implicit none
      integer, intent (in) :: id
    end subroutine shldel
 
    subroutine shlell(nx,ny,na,nb,ang)
      implicit none
      integer, intent (in) :: nx,ny,na,nb
      real, intent (in) :: ang
    end subroutine shlell
 
    subroutine shlind(id)
      implicit none
      integer, intent (out) :: id
    end subroutine shlind
 
    subroutine shlpie(nx,ny,nr,alph,beta)
      implicit none
      integer, intent (in) :: nx,ny,nr
      real, intent (in) :: alph,beta
    end subroutine shlpie
 
    subroutine shlpol(nxray,nyray,n)
      implicit none
      integer, intent (in) :: n
      integer, dimension (n), intent (in) :: nxray,nyray
    end subroutine shlpol
 
    subroutine shlrct(nx,ny,nw,nh,ang)
      implicit none
      integer, intent (in) :: nx,ny,nw,nh
      real,    intent (in) :: ang
    end subroutine shlrct
 
    subroutine shlrec(nx,ny,nw,nh)
      implicit none
      integer, intent (in) :: nx,ny,nw,nh
    end subroutine shlrec
 
    subroutine shlres(nn)
      implicit none
      integer, intent (in) :: nn
    end subroutine shlres
 
    subroutine shlsur()
    end subroutine shlsur
 
    subroutine shlvis(id,cvis)
      implicit none
      integer, intent (in) :: id
      character (len = *), intent (in) :: cvis
    end subroutine shlvis
 
    subroutine simplx()
    end subroutine simplx
 
    subroutine skipfl(nlu,nbyt,istat)
      implicit none
      integer, intent (in) :: nlu,nbyt
      integer, intent (out) :: istat
    end subroutine skipfl
 
    subroutine smxalf(calph,ca,cb,n)
      implicit none
      character (len = *), intent (in) :: calph,ca,cb
      integer, intent (in) :: n
    end subroutine smxalf
 
    subroutine solid()
    end subroutine solid
 
    subroutine sortr1(x,n,copt)
      implicit none
      integer, intent (in) :: n
      real, dimension (n), intent (in out) :: x
      character (len = *), intent (in) :: copt
    end subroutine sortr1
 
    subroutine sortr2(x,y,n,copt)
      implicit none
      integer, intent (in) :: n
      real, dimension (n), intent (in out) :: x,y
      character (len = *), intent (in) :: copt
    end subroutine sortr2

    subroutine sphe3d(xm,ym,zm,r,nsk1,nsk2)
      implicit none
      real, intent (in) :: xm,ym,zm,r
      integer, intent (in) :: nsk1,nsk2
    end subroutine sphe3d
 
    subroutine spline(x,y,n,xray,yray,ndat)
      implicit none
      integer, intent (in) :: n
      real, dimension (n),intent (in) :: x,y
      real, dimension (n),intent (out) :: xray,yray
      integer, intent (out) :: ndat
    end subroutine spline
 
    subroutine splmod(k,n)
      implicit none
      integer, intent (in) :: k,n
    end subroutine splmod
 
    subroutine strt3d(x,y,z)
      implicit none
      real, intent (in) :: x,y,z
    end subroutine strt3d
 
    subroutine strtpt(x,y)
      implicit none
      real, intent (in) :: x,y
    end subroutine strtpt
 
    subroutine surclr(itop,ibot)
      implicit none
      integer, intent (in) :: itop,ibot
    end subroutine surclr
 
    subroutine surfce(xray,ixdim,yray,iydim,zmat)
      implicit none
      integer, intent (in) :: ixdim,iydim
      real, dimension (ixdim), intent (in) :: xray
      real, dimension (iydim), intent (in) :: yray
      real, dimension (ixdim,iydim), intent (in) :: zmat
    end subroutine surfce
 
    subroutine surfcp(zfun,a1,a2,astp,b1,b2,bstp)
      implicit none
      interface
         function zfun(x,y,iopt)
           implicit none
           real, intent (in) :: x,y
           integer, intent (in) :: iopt
           real :: zfun
         end function zfun
      end interface
      real, intent (in) :: a1,a2,astp,b1,b2,bstp
    end subroutine surfcp

    subroutine suriso(xray,nx,yray,ny,zray,nz,wmat,wlev)
      implicit none
      integer, intent (in) :: nx,ny,nz
      real, dimension (nx), intent (in) :: xray
      real, dimension (ny), intent (in) :: yray
      real, dimension (nz), intent (in) :: zray
      real, dimension (nx,ny,nz), intent (in) :: wmat
      real, intent (in) :: wlev
    end subroutine suriso
 
    subroutine surfun(zfun,ixpts,xdel,iypts,ydel)
      implicit none
      interface
         function zfun(x,y)
           implicit none
           real, intent (in) :: x,y
           real :: zfun
         end function zfun
      end interface
      integer, intent (in) :: ixpts,iypts
      real, intent (in) :: xdel,ydel
    end subroutine surfun
 
    subroutine surmat(zmat,ixdim,iydim,ixp,iyp)
      implicit none
      integer, intent (in) :: ixdim,iydim,ixp,iyp
      real, dimension (ixdim,iydim), intent (in) :: zmat
    end subroutine surmat
 
    subroutine surmsh(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine surmsh

    subroutine suropt(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine suropt
 
    subroutine surshd(xray,ixdim,yray,iydim,zmat)
      implicit none
      integer, intent (in) :: ixdim,iydim
      real, dimension (ixdim), intent (in) :: xray
      real, dimension (iydim), intent (in) :: yray
      real, dimension (ixdim,iydim), intent (in) :: zmat
    end subroutine surshd
 
    subroutine sursze(ax,ex,ay,ey)
      implicit none
      real, intent (in) :: ax,ex,ay,ey
    end subroutine sursze
 
    subroutine surtri(xray,yray,zray,n,i1ray,i2ray,i3ray,ntri)
      implicit none
      integer, intent (in) :: n,ntri
      real, dimension (n), intent (in) :: xray,yray,zray
      integer, dimension (ntri), intent (in) :: i1ray,i2ray,i3ray
    end subroutine surtri

    subroutine survis(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine survis
 
    subroutine swapi2(iray,n)
      implicit none
      integer, intent (in) :: n
      integer (kind=selected_int_kind(4)), dimension (n), &
      &       intent (in out) :: iray
    end subroutine swapi2
 
    subroutine swapi4(iray,n)
      implicit none
      integer, intent (in) :: n
      integer (kind=selected_int_kind(9)), dimension (n), &
      &       intent (in out) :: iray
    end subroutine swapi4

    subroutine swgatt (id,cval,copt)
      implicit none
      integer, intent (in) :: id
      character (len=*), intent (in) :: cval,copt
    end subroutine swgatt
 
    subroutine swgbox(id,ival)
      implicit none
      integer, intent (in) :: id,ival
    end subroutine swgbox
 
    subroutine swgbut(id,ival)
      implicit none
      integer, intent (in) :: id,ival
    end subroutine swgbut
 
    subroutine swgcb (id, callbk, iray)
      implicit none
      integer, intent (in) :: id
      integer, dimension (1), intent (in) :: iray
 
      interface
         subroutine callbk (id,iray)
           implicit none
           integer, intent (in) :: id
           integer, dimension (1), intent (in) :: iray
         end subroutine callbk
      end interface
    end subroutine swgcb

    subroutine swgcb2 (id, callbk)
      implicit none
      integer, intent (in) :: id
 
      interface
         subroutine callbk (id,irow,icol)
           implicit none
           integer, intent (in) :: id,irow,icol
         end subroutine callbk
      end interface
    end subroutine swgcb2

    subroutine swgcbk (id, callbk)
      implicit none
      integer, intent (in) :: id
 
      interface
         subroutine callbk (id)
           implicit none
           integer, intent (in) :: id
         end subroutine callbk
      end interface
    end subroutine swgcbk


    subroutine swgclr(xr,xg,xb,copt)
      implicit none
      real, intent (in) :: xr,xg,xb
      character (len=*), intent (in) :: copt
    end subroutine swgclr

    subroutine swgdrw(x)
      implicit none
      real, intent (in) :: x
    end subroutine swgdrw
 
    subroutine swgfil(id,cstr)
      implicit none
      integer, intent (in) :: id
      character (len=*), intent (in) :: cstr
    end subroutine swgfil

    subroutine swgflt(id,xval,ndig)
      implicit none
      integer, intent (in) :: id,ndig
      real, intent (in) :: xval
    end subroutine swgflt
 
    subroutine swgfnt(cstr,n)
      implicit none
      character (len=*), intent (in) :: cstr
      integer, intent (in) :: n
    end subroutine swgfnt

    subroutine swgfoc(id)
      implicit none
      integer, intent (in) :: id
    end subroutine swgfoc

    subroutine swghlp(cstr)
      implicit none
      character (len=*), intent (in) :: cstr
    end subroutine swghlp

    subroutine swgint(id,iv)
      implicit none
      integer, intent (in) :: id,iv
    end subroutine swgint
 
    subroutine swgjus (ctype,cwidg)
      implicit none
      character (len=*), intent (in) :: ctype,cwidg
    end subroutine swgjus
 
    subroutine swglis(id,ival)
      implicit none
      integer, intent (in) :: id,ival
    end subroutine swglis
 
    subroutine swgmix(c,cstr)
      implicit none
      character (len=*), intent (in) :: c,cstr
    end subroutine swgmix
 
    subroutine swgmod(cmod)
      implicit none
      character (len=*), intent (in) :: cmod
    end subroutine swgmod
 
    subroutine swgmrg(ival,cstr)
      implicit none
      integer, intent (in) :: ival
      character (len=*), intent (in) :: cstr
    end subroutine swgmrg
 
    subroutine swgoff(nx,ny)
      implicit none
      integer, intent (in) :: nx,ny
    end subroutine swgoff

    subroutine swgopt (cval,copt)
      implicit none
      character (len=*), intent (in) :: cval,copt
    end subroutine swgopt
 
    subroutine swgpop (copt)
      implicit none
      character (len=*), intent (in) :: copt
    end subroutine swgpop
 
    subroutine swgpos(nx,ny)
      implicit none
      integer, intent (in) :: nx,ny
    end subroutine swgpos

    subroutine swgray(xray,n,copt)
      implicit none
      integer, intent (in) :: n
      real, dimension (n), intent (in) :: xray 
      character (len=*), intent (in) :: copt
    end subroutine swgray
 
    subroutine swgscl(id,xval)
      implicit none
      integer, intent (in) :: id
      real, intent (in) :: xval
    end subroutine swgscl
 
    subroutine swgsiz(nx,ny)
      implicit none
      integer, intent (in) :: nx,ny
    end subroutine swgsiz

    subroutine swgspc(x, y)
      implicit none
      real, intent (in) :: x,y
    end subroutine swgspc

    subroutine swgstp(x)
      implicit none
      real, intent (in) :: x
    end subroutine swgstp

    subroutine swgtbf(id,xval,ndig,irow,icol,copt)
      implicit none
      integer, intent (in) :: id,ndig,irow,icol
      real, intent (in) :: xval 
      character (len=*), intent (in) :: copt
    end subroutine swgtbf

    subroutine swgtbi(id,ival,irow,icol,copt)
      implicit none
      integer, intent (in) :: id,ival,irow,icol
      character (len=*), intent (in) :: copt
    end subroutine swgtbi

    subroutine swgtbl(id,xray,n,ndig,idx,copt)
      implicit none
      integer, intent (in) :: id,n,ndig,idx
      real, dimension (n), intent (in) :: xray 
      character (len=*), intent (in) :: copt
    end subroutine swgtbl

    subroutine swgtbs(id,cstr,irow,icol,copt)
      implicit none
      integer, intent (in) :: id,irow,icol 
      character (len=*), intent (in) :: cstr,copt
    end subroutine swgtbs

    subroutine swgtit(cstr)
      implicit none
      character (len=*), intent (in) :: cstr
    end subroutine swgtit
 
    subroutine swgtxt(id,cstr)
      implicit none
      integer, intent (in) :: id
      character (len=*), intent (in) :: cstr
    end subroutine swgtxt
 
    subroutine swgtyp (ctype,cwidg)
      implicit none
      character (len=*), intent (in) :: ctype,cwidg
    end subroutine swgtyp

    subroutine swgval(id,xval)
      implicit none
      integer, intent (in) :: id
      real, intent (in) :: xval
    end subroutine swgval
 
    subroutine swgwin(nx,ny,nw,nh)
      implicit none
      integer, intent (in) :: nx,ny,nw,nh
    end subroutine swgwin
 
    subroutine swgwth (nwth)
      implicit none
      integer, intent (in) :: nwth
    end subroutine swgwth

    subroutine symb3d(i,x,y,z)
      implicit none
      integer, intent (in) :: i
      real, intent (in) :: x,y,z
    end subroutine symb3d
 
    subroutine symbol(i,nx,ny)
      implicit none
      integer, intent (in) :: i,nx,ny
    end subroutine symbol
 
    subroutine symfil(cdv,cst)
      implicit none
      character (len=*), intent (in) :: cdv,cst
    end subroutine symfil
 
    subroutine symrot(xrot)
      implicit none
      real, intent (in) :: xrot
    end subroutine symrot
 
    subroutine tellfl(nlu,nbyt)
      implicit none
      integer, intent (in) :: nlu
      integer, intent (out) :: nbyt
    end subroutine tellfl

    subroutine texmod(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine texmod
 
    subroutine texopt(copt,ctype)
      implicit none
      character (len = *), intent (in) :: copt,ctype
    end subroutine texopt

    subroutine texval(x,copt)
      implicit none
      real, intent (in) :: x 
      character (len = *), intent (in) :: copt
    end subroutine texval

    subroutine thkc3d(x)
      implicit none
      real, intent (in) :: x
    end subroutine thkc3d

    subroutine thkcrv(i)
      implicit none
      integer, intent (in) :: i
    end subroutine thkcrv
 
    subroutine thrfin()
      implicit none
    end subroutine thrfin

    subroutine thrini(i)
      implicit none
      integer, intent (in) :: i
    end subroutine thrini

    subroutine ticks(i,cax)
      implicit none
      integer, intent (in) :: i
      character (len = *), intent (in) :: cax
    end subroutine ticks
 
    subroutine ticlen(i1,i2)
      implicit none
      integer, intent (in) :: i1,i2
    end subroutine ticlen

    subroutine ticmod(copt,cax)
      implicit none
      character (len = *), intent (in) :: copt,cax
    end subroutine ticmod
 
    subroutine ticpos(copt,cax)
      implicit none
      character (len = *), intent (in) :: copt,cax
    end subroutine ticpos

    subroutine tifmod(n,cval,copt)
      implicit none
      character (len = *), intent (in) :: cval,copt
      integer, intent (in) :: n
    end subroutine tifmod
 
    subroutine tiforg(nx,ny)
      implicit none
      integer, intent (in) :: nx,ny
    end subroutine tiforg
 
    subroutine tifwin(nx,ny,nw,nh)
      implicit none
      integer, intent (in) :: nx,ny,nw,nh
    end subroutine tifwin
 
    subroutine timopt()
    end subroutine timopt
 
    subroutine titjus(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine titjus
 
    subroutine title()
    end subroutine title
 
    subroutine titlin(cstr,j)
      implicit none
      character (len = *), intent (in) :: cstr
      integer, intent (in) :: j
    end subroutine titlin
 
    subroutine titpos(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine titpos

    subroutine torus3d(xm,ym,zm,r1,r2,h,a1,a2,nsk1,nsk2)
      implicit none
      real, intent (in) :: xm,ym,zm,r1,r2,h,a1,a2
      integer, intent (in) :: nsk1,nsk2
    end subroutine torus3d

    subroutine tprfin()
    end subroutine tprfin

    subroutine tprini()
    end subroutine tprini

    subroutine tprmod(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine tprmod

    subroutine tprval(x)
      implicit none
      real, intent (in) :: x
    end subroutine tprval
 
    subroutine tr3res()
    end subroutine tr3res
 
    subroutine tr3rot(a,b,c)
      implicit none
      real, intent (in) :: a,b,c
    end subroutine tr3rot
 
    subroutine tr3scl(xscl,yscl,zscl)
      implicit none
      real, intent (in) :: xscl,yscl,zscl
    end subroutine tr3scl
 
    subroutine tr3shf(xshf,yshf,zshf)
      implicit none
      real, intent (in) :: xshf,yshf,zshf
    end subroutine tr3shf
 
    subroutine trfco1(xray,n,cfrom,cto)
      implicit none
      integer, intent (in) :: n
      real, dimension (n), intent (in out) :: xray
      character (len = *), intent(in) :: cfrom, cto
    end subroutine trfco1
 
    subroutine trfco2(xray,yray,n,cfrom,cto)
      implicit none
      integer, intent (in) :: n
      real, dimension (n), intent (in out) :: xray,yray
      character (len = *), intent(in) :: cfrom, cto
    end subroutine trfco2
 
    subroutine trfco3(xray,yray,zray,n,cfrom,cto)
      implicit none
      integer, intent (in) :: n
      real, dimension (n), intent (in out) :: xray,yray,zray
      character (len = *), intent(in) :: cfrom, cto
    end subroutine trfco3
 
    subroutine trfdat(ndays,id,im,iy)
      implicit none
      integer, intent (in) :: ndays
      integer, intent (out) :: id,im,iy
    end subroutine trfdat
 
    subroutine trfmat(zmat,nx,ny,zmat2,nx2,ny2)
      implicit none
      integer, intent (in) :: nx,ny,nx2,ny2
      real, dimension (nx,ny), intent (in) :: zmat
      real, dimension (nx2,ny2), intent (out) :: zmat2
    end subroutine trfmat

    subroutine trfrel(x,y,n)
      implicit none
      integer, intent (in) :: n
      real, dimension (n), intent (in out) :: x,y
    end subroutine trfrel
 
    subroutine trfres()
    end subroutine trfres
 
    subroutine trfrot(xang,nx,ny)
      implicit none
      real, intent (in) :: xang
      integer, intent (in) :: nx,ny
    end subroutine trfrot
 
    subroutine trfscl(xscl,yscl)
      implicit none
      real, intent (in) :: xscl,yscl
    end subroutine trfscl
 
    subroutine trfshf(nxshf,nyshf)
      implicit none
      integer, intent (in) :: nxshf,nyshf
    end subroutine trfshf
 
    subroutine tria3d(x,y,z)
      implicit none
      real, dimension (3), intent (in) :: x,y,z
    end subroutine tria3d

    subroutine triang(xray,yray,n,i1ray,i2ray,i3ray,nmax,ntri)
      implicit none
      integer, intent (in) :: n,nmax
      real, dimension (n), intent (in) :: xray,yray
      integer, dimension (nmax), intent (out) :: i1ray,i2ray,i3ray
      integer, intent (out) :: ntri 
    end subroutine triang

    subroutine trifll(x,y)
      implicit none
      real, dimension (3), intent (in) :: x,y
    end subroutine trifll

    subroutine triplx()
    end subroutine triplx
 
    subroutine tripts(xray,yray,zray,n,i1ray,i2ray,i3ray,ntri,zlev, &
               &      xpts, ypts, maxpts, iray, maxray, nlins)
      implicit none
      integer, intent (in) :: n,ntri,maxpts,maxray
      real, dimension (n), intent (in) :: xray,yray,zray
      real, dimension (maxpts), intent (out) :: xpts,ypts
      integer, dimension (ntri), intent (in) :: i1ray,i2ray,i3ray
      integer, dimension (maxray), intent (out) :: iray
      real, intent (in) :: zlev
      integer, intent (out) :: nlins
    end subroutine tripts

    function trmlen(cstr)
      implicit none
      character (len = *), intent (in) :: cstr
      real :: trmlen
    end function trmlen
 
    subroutine tube3d(x1,y1,z1,x2,y2,z2,r,nsk1,nsk2)
      implicit none
      real, intent (in) :: x1,y1,z1,x2,y2,z2,r
      integer, intent (in) :: nsk1,nsk2
    end subroutine tube3d

    subroutine txtjus(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine txtjus
 
    subroutine unit(i)
      implicit none
      integer, intent (in) :: i
    end subroutine unit
 
    subroutine units(cmod)
      implicit none
      character (len = *), intent (in) :: cmod
    end subroutine units

    subroutine upstr(cstr)
      implicit none
      character (len = *), intent (in out) :: cstr
    end subroutine upstr
 
    subroutine usrpie(iseg,xdat,xper,nrad,noff,ang,nvbox,idrw,iann)
      implicit none
      integer, intent (in out) :: iseg,nrad,noff,nvbox,idrw,iann
      real,    intent (in out) :: xdat,xper,ang
    end subroutine usrpie
 
    subroutine utfint(cstr,iray,n,nl)
      implicit none
      character (len=*), intent (in) :: cstr
      integer, intent (in) :: n
      integer, dimension (n), intent (out) :: iray
      integer, intent (out) :: nl
    end subroutine utfint

    subroutine vang3d(a)
      implicit none
      real, intent (in) :: a
    end subroutine vang3d
 
    subroutine vclp3d(x1,x2)
      implicit none
      real, intent (in) :: x1,x2
    end subroutine vclp3d

    subroutine vecclr(iclr)
      implicit none
      integer, intent (in) :: iclr
    end subroutine vecclr

    subroutine vecf3d(xv,yv,zv,xp,yp,zp,n,ivec)
      implicit none
      integer, intent (in) :: n,ivec
      real, dimension (n), intent (in) :: xv,yv,zv,xp,yp,zp
    end subroutine vecf3d
 
    subroutine vecfld(xv,yv,xp,yp,n,ivec)
      implicit none
      integer, intent (in) :: n,ivec
      real, dimension (n), intent (in) :: xv,yv,xp,yp
    end subroutine vecfld

    subroutine vecopt(x,copt)
      implicit none
      real, intent (in) :: x
      character (len = *), intent (in) :: copt
    end subroutine vecopt

    subroutine vector(ix1,iy1,ix2,iy2,ivec)
      implicit none
      integer, intent (in) :: ix1,iy1,ix2,iy2,ivec
    end subroutine vector
 
    subroutine vectr3(x1,y1,z1,x2,y2,z2,ivec)
      implicit none
      real, intent (in) :: x1,y1,z1,x2,y2,z2
      integer, intent (in) :: ivec
    end subroutine vectr3
 
    subroutine vfoc3d(x,y,z,cview)
      implicit none
      real, intent (in) :: x,y,z
      character (len = *), intent (in) :: cview
    end subroutine vfoc3d
 
    subroutine view3d(x,y,z,cview)
      implicit none
      real, intent (in) :: x,y,z
      character (len = *), intent (in) :: cview
    end subroutine view3d
 
    subroutine vkxbar(i)
      implicit none
      integer, intent (in) :: i
    end subroutine vkxbar
 
    subroutine vkybar(i)
      implicit none
      integer, intent (in) :: i
    end subroutine vkybar
 
    subroutine vkytit(i)
      implicit none
      integer, intent (in) :: i
    end subroutine vkytit
 
    subroutine vltfil(cfl, copt)
      implicit none
      character (len = *), intent (in) :: cfl, copt
    end subroutine vltfil

    subroutine vtx3d(xray,yray,zray,n,copt)
      implicit none
      integer, intent (in) :: n
      real, dimension (n), intent (in) :: xray,yray,zray
      character (len = *), intent (in) :: copt
    end subroutine vtx3d

    subroutine vtxc3d(xray,yray,zray,ic,n,copt)
      implicit none
      integer, intent (in) :: n
      real, dimension (n), intent (in) :: xray,yray,zray
      integer, dimension (n), intent (in) :: ic
      character (len = *), intent (in) :: copt
    end subroutine vtxc3d

    subroutine vtxn3d(xray,yray,zray,xn,yn,zn,n,copt)
      implicit none
      integer, intent (in) :: n
      real, dimension (n), intent (in) :: xray,yray,zray,xn,yn,zn
      character (len = *), intent (in) :: copt
    end subroutine vtxn3d
 
    subroutine vup3d(a)
      implicit none
      real, intent (in) :: a
    end subroutine vup3d
 
    subroutine wgapp(ip,cstr,id)
      implicit none
      character (len = *), intent (in) :: cstr
      integer, intent (in)  :: ip
      integer, intent (out) :: id
    end subroutine wgapp
 
    subroutine wgbas(ip,copt,id)
      implicit none
      character (len = *), intent (in) :: copt
      integer, intent (in)  :: ip
      integer, intent (out) :: id
    end subroutine wgbas
 
    subroutine wgbox(ip,cstr,isel,id)
      implicit none
      character (len = *), intent (in) :: cstr
      integer, intent (in)  :: ip,isel
      integer, intent (out) :: id
    end subroutine wgbox
 
    subroutine wgbut(ip,cstr,ival,id)
      implicit none
      character (len = *), intent (in) :: cstr
      integer, intent (in)  :: ip,ival
      integer, intent (out) :: id
    end subroutine wgbut
 
    subroutine wgcmd(ip,clab,cstr,id)
      implicit none
      character (len = *), intent (in) :: clab,cstr
      integer, intent (in)  :: ip
      integer, intent (out) :: id
    end subroutine wgcmd

    subroutine wgdlis(ip,cstr,isel,id)
      implicit none
      character (len = *), intent (in) :: cstr
      integer, intent (in)  :: ip,isel
      integer, intent (out) :: id
    end subroutine wgdlis

    subroutine wgdraw(ip,id)
      implicit none
      integer, intent (in)  :: ip
      integer, intent (out) :: id
    end subroutine wgdraw
 
    subroutine wgfil(ip,clab,cstr,cmask,id)
      implicit none
      character (len = *), intent (in) :: clab,cstr,cmask
      integer, intent (in)  :: ip
      integer, intent (out) :: id
    end subroutine wgfil
 
    subroutine wgfin()
    end subroutine wgfin
 
    subroutine wgini(ctype,id)
      implicit none
      character (len = *), intent (in) :: ctype
      integer, intent (out) :: id
    end subroutine wgini
 
    subroutine wglab(ip,cstr,id)
      implicit none
      character (len = *), intent (in) :: cstr
      integer, intent (in)  :: ip
      integer, intent (out) :: id
    end subroutine wglab
 
    subroutine wglis(ip,cstr,isel,id)
      implicit none
      character (len = *), intent (in) :: cstr
      integer, intent (in)  :: ip,isel
      integer, intent (out) :: id
    end subroutine wglis
 
    subroutine wgltxt(ip,clab,cstr,iper,id)
      implicit none
      character (len = *), intent (in) :: clab,cstr
      integer, intent (in)  :: ip,iper
      integer, intent (out) :: id
    end subroutine wgltxt
 
    subroutine wgok(ip,id)
      implicit none
      integer, intent (in)  :: ip
      integer, intent (out) :: id
    end subroutine wgok

    subroutine wgpbar(ip,x1,x2,xstp,id)
      implicit none
      integer, intent (in)  :: ip
      real, intent (in)     :: x1,x2,xstp
      integer, intent (out) :: id
    end subroutine wgpbar
 
    subroutine wgpbut(ip,cstr,id)
      implicit none
      character (len = *), intent (in) :: cstr
      integer, intent (in)  :: ip
      integer, intent (out) :: id
    end subroutine wgpbut
 
    subroutine wgpop(ip,cstr,id)
      implicit none
      character (len = *), intent (in) :: cstr
      integer, intent (in)  :: ip
      integer, intent (out) :: id
    end subroutine wgpop
 
    subroutine wgquit(ip,id)
      implicit none
      integer, intent (in)  :: ip
      integer, intent (out) :: id
    end subroutine wgquit
 
    subroutine wgscl(ip,cstr,x1,x2,xval,ndez,id)
      implicit none
      character (len = *), intent (in) :: cstr
      integer, intent (in)  :: ip,ndez
      real, intent (in)     :: x1,x2,xval
      integer, intent (out) :: id
    end subroutine wgscl
 
    subroutine wgstxt(ip,nsize,nmax,id)
      implicit none
      integer, intent (in)  :: ip,nsize,nmax
      integer, intent (out) :: id
    end subroutine wgstxt

    subroutine wgtbl(ip,nrows,ncols,id)
      implicit none
      integer, intent (in)  :: ip,nrows,ncols
      integer, intent (out) :: id
    end subroutine wgtbl

    subroutine wgtxt(ip,cstr,id)
      implicit none
      character (len = *), intent (in) :: cstr
      integer, intent (in)  :: ip
      integer, intent (out) :: id
    end subroutine wgtxt
 
    subroutine widbar(i)
      implicit none
      integer, intent (in) :: i
    end subroutine widbar
 
    subroutine wimage(cfl)
      implicit none
      character (len = *), intent (in) :: cfl
    end subroutine wimage
 
    subroutine winapp(copt)
      implicit none
      character (len = *), intent (in) :: copt
    end subroutine winapp

    subroutine windbr(xk,nx,ny,nw,a)
      implicit none
      real, intent (in) :: xk,a
      integer, intent (in) :: nx,ny,nw
    end subroutine windbr
 
    subroutine window(nx,ny,nw,nh)
      implicit none
      integer, intent (in) :: nx,ny,nw,nh
    end subroutine window
 
    subroutine winfnt(cfnt)
      implicit none
      character (len = *), intent (in) :: cfnt
    end subroutine winfnt
 
    subroutine winid(id)
      implicit none
      integer, intent (out) :: id
    end subroutine winid
 
   subroutine winkey(cmod)
      implicit none
      character (len = *), intent (in) :: cmod
    end subroutine winkey
  
    subroutine winmod(cmod)
      implicit none
      character (len = *), intent (in) :: cmod
    end subroutine winmod
 
    subroutine winopt(iopt,copt)
      implicit none
      character (len = *), intent (in) :: copt
      integer, intent (in) :: iopt
    end subroutine winopt

    subroutine winsiz(nw,nh)
      implicit none
      integer, intent (in) :: nw,nh
    end subroutine winsiz
 
    subroutine wintit(cstr)
      implicit none
      character (len = *), intent (in) :: cstr
    end subroutine wintit
 
    subroutine wmfmod(cmod,ckey)
      implicit none
      character (len = *), intent (in) :: cmod,ckey
    end subroutine wmfmod

    subroutine world()
    end subroutine world
 
    subroutine wpixel(ix,iy,n)
      implicit none
      integer, intent (in) :: ix,iy,n
    end subroutine wpixel
 
    subroutine wpixls(iray,ix,iy,nw,nh)
      implicit none
      integer, intent (in) :: ix,iy,nw,nh
      character (len=1), intent (in), dimension (*) :: iray
    end subroutine wpixls
 
    subroutine wpxrow(iray,ix,iy,n)
      implicit none
      integer, intent (in) :: ix,iy,n
      character (len=1), intent (in), dimension (n) :: iray
    end subroutine wpxrow
 
    subroutine writfl(nlu,iray,nbyt,istat)
      implicit none
      integer, intent (in) :: nlu,nbyt
      character (len=1), intent (in), dimension (nbyt) :: iray
      integer, intent (out) :: istat
    end subroutine writfl
 
    subroutine wtiff(cfl)
      implicit none
      character (len=*), intent (in) :: cfl
    end subroutine wtiff
 
    subroutine x11fnt(cfnt, copt)
      implicit none
      character (len = *), intent (in) :: cfnt, copt
    end subroutine x11fnt

    subroutine x11mod(cmod)
      implicit none
      character (len = *), intent (in) :: cmod
    end subroutine x11mod
 
    function x2dpos(x,y)
      implicit none
      real, intent (in) :: x,y
      real :: x2dpos
    end function x2dpos
 
    function x3dabs(x,y,z)
      implicit none
      real, intent (in) :: x,y,z
      real :: x3dabs
    end function x3dabs
 
    function x3dpos(x,y,z)
      implicit none
      real, intent (in) :: x,y,z
      real :: x3dpos
    end function x3dpos
 
    function x3drel(x,y,z)
      implicit none
      real, intent (in) :: x,y,z
      real :: x3drel
    end function x3drel
 
    subroutine xaxgit()
    end subroutine xaxgit
 
    subroutine xaxis(a,b,or,step,il,cstr,it,ix,iy)
      implicit none
      real, intent (in) :: a,b,or,step
      character (len = *), intent (in) :: cstr
      integer, intent (in) :: il,it,ix,iy
    end subroutine xaxis
 
    subroutine xaxlg(a,b,or,step,il,cstr,it,ix,iy)
      implicit none
      real, intent (in) :: a,b,or,step
      character (len = *), intent (in) :: cstr
      integer, intent (in) :: il,it,ix,iy
    end subroutine xaxlg
 
    subroutine xaxmap(a,b,or,step,cstr,it,iy)
      implicit none
      real, intent (in) :: a,b,or,step
      integer, intent (in) :: it,iy
      character (len = *), intent (in) :: cstr
    end subroutine xaxmap
 
    subroutine xcross()
    end subroutine xcross
 
    subroutine xdraw(xx,yy)
      implicit none
      real, intent (in) :: xx,yy
    end subroutine xdraw
 
    function xinvrs(i)
      implicit none
      integer, intent (in) :: i
      real :: xinvrs
    end function xinvrs
 
    subroutine xmove(x,y)
      implicit none
      real, intent (in) :: x,y
    end subroutine xmove
 
    function xposn(x)
      implicit none
      real, intent (in) :: x
      real :: xposn
    end function xposn
 
    function y2dpos(x,y)
      implicit none
      real, intent (in) :: x,y
      real :: y2dpos
    end function y2dpos
 
    function y3dabs(x,y,z)
      implicit none
      real, intent (in) :: x,y,z
      real :: y3dabs
    end function y3dabs
 
    function y3dpos(x,y,z)
      implicit none
      real, intent (in) :: x,y,z
      real :: y3dpos
    end function y3dpos
 
    function y3drel(x,y,z)
      implicit none
      real, intent (in) :: x,y,z
      real :: y3drel
    end function y3drel
 
    subroutine yaxgit()
    end subroutine yaxgit
 
    subroutine yaxis(a,b,or,step,il,cstr,it,ix,iy)
      implicit none
      real, intent (in) :: a,b,or,step
      character (len = *), intent (in) :: cstr
      integer, intent (in) :: il,it,ix,iy
    end subroutine yaxis
 
    subroutine yaxlg(a,b,or,step,il,cstr,it,ix,iy)
      implicit none
      real, intent (in) :: a,b,or,step
      character (len = *), intent (in) :: cstr
      integer, intent (in) :: il,it,ix,iy
    end subroutine yaxlg
 
    subroutine yaxmap(a,b,or,step,cstr,it,ix)
      implicit none
      real, intent (in) :: a,b,or,step
      integer, intent (in) :: it,ix
      character (len = *), intent (in) :: cstr
    end subroutine yaxmap
 
    subroutine ycross()
    end subroutine ycross
 
    function yinvrs(i)
      implicit none
      integer, intent (in) :: i
      real :: yinvrs
    end function yinvrs
 
    function yposn(y)
      implicit none
      real, intent (in) :: y
      real :: yposn
    end function yposn
 
    function z3dpos(x,y,z)
      implicit none
      real, intent (in) :: x,y,z
      real :: z3dpos
    end function z3dpos
 
    subroutine zaxis(a,b,or,step,il,cstr,it,idir,ix,iy)
      implicit none
      real, intent (in) :: a,b,or,step
      integer, intent (in) :: il,it,idir,ix,iy
      character (len = *), intent (in) :: cstr
    end subroutine zaxis
 
    subroutine zaxlg(a,b,or,step,il,cstr,it,idir,ix,iy)
      implicit none
      real, intent (in) :: a,b,or,step
      integer, intent (in) :: il,it,idir,ix,iy
      character (len = *), intent (in) :: cstr
    end subroutine zaxlg

    subroutine zbfers()
    end subroutine zbfers
 
    subroutine zbffin()
    end subroutine zbffin
 
    subroutine zbfini(iret)
      implicit none
      integer, intent (out) :: iret
    end subroutine zbfini
 
    subroutine zbfres()
    end subroutine zbfres

    subroutine zbflin(x1,y1,z1,x2,y2,z2)
      implicit none
      real, intent (in) :: x1,y1,z1,x2,y2,z2
    end subroutine zbflin
 
    subroutine zbftri(x,y,z,ic)
      implicit none
      real, dimension (3), intent (in) :: x,y,z
      integer, dimension (3), intent (in) :: ic
    end subroutine zbftri
 
    subroutine zscale(a,e)
      implicit none
      real, intent (in) :: a,e
    end subroutine zscale
  end interface
end module dislin 
#endif
