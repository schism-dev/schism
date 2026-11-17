!   Copyright 2014 College of William and Mary
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

! Program to compute nodal factors and equilibrium arguements.
! Many thanks to the author of "tide_fac" and kai.li.jx#gmail.com

! WARNING: Default is real*4

module nodalfactor
  implicit none
  private
  
  integer,public,parameter :: ncnst=37
  real,parameter :: pi180=3.14159265/180.
  character(len=8),public,save :: cname(ncnst)
  real,public,save :: nodfac(ncnst),grterm(ncnst),pcst(ncnst)
  real*8,public,save :: speed(ncnst) !for computing freq
  real, save :: ds_orb,dp_orb,dh_orb,dp1_orb,dn_orb,di_orb,dnu_orb,dxi_orb,dnup_orb,dnup2_orb,dpc_orb
  
  public :: tidalfactor
  
  contains
  
  subroutine tidalfactor(xdays,bhr,iday,imo,iyr)
    real, intent(in) :: xdays,bhr !run length (days), offset to GMT (hr)
    integer, intent(in) :: iday,imo,iyr !day,month,year
    
    !character :: tmp, arg*20
    !common /cnsnam/ cname
    real :: month,rhrs,yr,day,hrm,dayj
    integer :: nc
!    integer :: ncon(200),numcon,nc,ic
    !common /cnst/ nodfac(ncnst),grterm(ncnst),speed(ncnst),p(ncnst)
    
!    open(unit=11,file='tide_fac.out',status='unknown')
    !if (iargc().eq.0) then
    !!  write(*,*) 'Pass "-h" option for command line arguments instruction' 
    !  write(*,*) 'SCHISM trunk/src/Utility/Tide/tf_main.f90 -- calculate nodal factors and re-write bctides.in'
    !  write(*,*) 'enter length of run time (days)'
    !  read(*,*) xdays
    !  write(*,*) 'enter start time - bhr,iday,imo,iyr (iyr e.g. 1992)'
    !  read(*,*) bhr,iday,imo,iyr
    !!  write(*,*) 'generate "bctides.in.out" file for SCHISM? ( "bctides.in" is required for input )'
    !!  read(*,*)tmp
    !!  if (tmp=='y'.or.tmp=='Y') lselfe=.true.
    !else
    !  write(*,*) 'command line parsing not supported yet'
    !  stop
    !endif
    rhrs=xdays*24.
    yr=iyr
    month=imo
    day=iday
    hrm=bhr+rhrs/2.
!    write(11,10) bhr,iday,imo,iyr
!    !write(*,10) bhr,iday,imo,iyr
!    10 format(' tidal factors starting: ',' hr-',f5.2,',  day-',i3,',  month-',i3,'  year-',i5,/)
!    !write(*,11) xdays
!    write(11,11) xdays
!    11 format(' for a run lasting ',f8.2,' days',//)
    !-- determine the julian time at beginning and middle of record
    dayj=dayjul(yr,month,day)
    !-- determine node factors at middle of record
    call nfacs(yr,dayj,hrm)
    !-- determine greenwich equil. terms at beginning of record
    call gterms(yr,dayj,bhr,dayj,hrm)
    
!    !Can output all 37 
!      numcon=8
!      ncon(1)=4
!      ncon(2)=6
!      ncon(3)=30
!      ncon(4)=26
!      ncon(5)=3
!      ncon(6)=1
!      ncon(7)=2
!      ncon(8)=35
!    write(11,*) 'const   node     eq arg (ref gm)'
!    write(11,1300)
!    1300 format(' name   factor    (deg) ',//)
!    do nc=1,numcon
!      ic=ncon(nc)
!    ! equilibrium arguement is referenced to the greenwich meridian
!      if (ic==0) then
!        write(11,2001) 'z0  ',1.0,0.0
!      else
!        write(11,2001) cname(ic),nodfac(ic),grterm(ic)
!      endif
!      2001 format(1x,a4,2x,f8.5,4x,f7.2,2x,f7.4)
!    enddo !nc
!    
    !Output all
!    do nc=1,ncnst
!      write(99,*) cname(nc), nodfac(nc),grterm(nc)
!    enddo !nc
  end subroutine tidalfactor
    
  subroutine nfacs(yr,dayj,hr)
  
  !-- calculates node factors for constituent tidal signal
  
  !-- the equations used in this routine come from:
  !         "manual of harmonic analysis and prediction of tides"
  !         by paul schureman, special publication #98, us coast
  !         and geodetic survey, department of commerce (1958).
  
  !-- if daym and hrm correspond to midyear, then this routine
  !-- returns the same values as found in table 14 of schureman.
  !---------------------------------------------------------------------
  
    real, intent(in) :: yr,dayj,hr
  !  character*8   cst(37)
    real :: i,n,nu,xi,p,pc,sini,sini2,sin2i,cosi2,tani2,qainv, &
           &rainv,eq73,eq74,eq75,eq76,eq77,eq78,eq149,eq207,eq215, &
           &eq227,eq235
    
  !  common/orbitf/ds_orb,dp_orb,dh,dp1,dn,di,dnu,dxi,dnup,dnup2,dpc
  !  common/cnst/fndcst(37),eqcst(37),acst(37),pcst(37)
  !  common/cnsnam/cst
  
  !-- constituent names:
    cname=(/'M2      ','S2      ','N2      ','K1      ',&
            'M4      ','O1      ','M6      ','MK3     ',&
            'S4      ','MN4     ','nu2     ','S6      ',&
            'mu2     ','2N2     ','OO1     ','lambda2 ',&
            'S1      ','M1      ','J1      ','Mm      ',&
            'Ssa     ','Sa      ','MSf     ','Mf      ',&
            'rho1    ','Q1      ','T2      ','R2      ',&
            '2Q1     ','P1      ','2SM2    ','M3      ',&
            'L2      ','2MK3    ','K2      ','M8      ',&
            'MS4     '/)  ! to be added: 2MS6 2SM6 SK3 2MN6 eta2 ups1 MNS2 SN4
  
  !-- orbital speeds (degrees/hour):
    speed=(/28.984104252,30.0000000,28.439729568,15.041068632,&
            57.968208468,13.943035584,86.952312720,44.025172884,&
            60.0,57.423833820,28.5125831,90.0,&
            27.9682084,27.8953548,16.139101680,29.4556253,&
            15.0,14.496693984,15.5854433,0.5443747,&
            0.0821373,0.0410686,1.0158957720,1.0980331,&
            13.4715145,13.3986609,29.9589333,30.0410667,&
            12.854286252,14.9589314,31.01589576,43.476156360,&
            29.5284789,42.927139836,30.0821373,115.936416972,&
            58.984104240/)
  
  !-- number of tide cycles per day per constituent:
    pcst=(/2.,2.,2.,1.,4.,1.,6.,3.,4.,4., &
          &2.,6.,2.,2.,1.,2.,1.,1.,1.,0., &
          &0.,0.,0.,0.,1.,1.,2.,2.,1.,1., &
          &2.,3.,2.,3.,2.,8.,4./)
  
    call orbit(yr,dayj,hr)
    n=dn_orb*pi180
    i=di_orb*pi180
    nu=dnu_orb*pi180
    xi=dxi_orb*pi180
    p=dp_orb*pi180
    pc=dpc_orb*pi180
    sini=sin(i)
    sini2=sin(i/2.)
    sin2i=sin(2.*i)
    cosi2=cos(i/2.)
    tani2=tan(i/2.)
  !-- equation 197, schureman
    qainv=sqrt(2.310+1.435*cos(2.*pc))
  !-- equation 213, schureman
    rainv=sqrt(1.-12.*tani2**2*cos(2.*pc)+36.*tani2**4)
  !-- variable names refer to equation numbers in schureman
    eq73=(2./3.-sini**2)/.5021
    eq74=sini**2/.1578
    eq75=sini*cosi2**2/.37988
    eq76=sin(2*i)/.7214
    eq77=sini*sini2**2/.0164
    eq78=(cosi2**4)/.91544
    eq149=cosi2**6/.8758
    eq207=eq75*qainv
    eq215=eq78*rainv
    eq227=sqrt(.8965*sin2i**2+.6001*sin2i*cos(nu)+.1006)
    eq235=.001+sqrt(19.0444*sini**4+2.7702*sini**2*cos(2.*nu)+.0981)
  !-- node factors for 37 constituents:
    nodfac(1)=eq78
    nodfac(2)=1.0
    nodfac(3)=eq78
    nodfac(4)=eq227
    nodfac(5)=nodfac(1)**2
    nodfac(6)=eq75
    nodfac(7)=nodfac(1)**3
    nodfac(8)=nodfac(1)*nodfac(4)
    nodfac(9)=1.0
    nodfac(10)=nodfac(1)**2
    nodfac(11)=eq78
    nodfac(12)=1.0
    nodfac(13)=eq78
    nodfac(14)=eq78
    nodfac(15)=eq77
    nodfac(16)=eq78
    nodfac(17)=1.0
  !** equation 207 not producing correct answer for m1
  !**set node factor for m1 = 0 until can further research
    nodfac(18)=0.
  !  nodfac(18)=eq207
    nodfac(19)=eq76
    nodfac(20)=eq73
    nodfac(21)=1.0
    nodfac(22)=1.0
    nodfac(23)=eq78
    nodfac(24)=eq74
    nodfac(25)=eq75
    nodfac(26)=eq75
    nodfac(27)=1.0
    nodfac(28)=1.0
    nodfac(29)=eq75
    nodfac(30)=1.0
    nodfac(31)=eq78
    nodfac(32)=eq149
  !** equation 215 not producing correct answer for l2
  !** set node factor for l2 = 0 until can further research
    nodfac(33)=0.
  !  nodfac(33)=eq215
    nodfac(34)=nodfac(1)**2*nodfac(4)
    nodfac(35)=eq235
    nodfac(36)=nodfac(1)**4
    nodfac(37)=eq78
  end subroutine nfacs
  
  subroutine gterms(yr,dayj,hr,daym,hrm)
  !-- calculates equilibrium arguments v0+u for constituent tide
  
  !-- the equations used in this routine come from:
  !         "manual of harmonic analysis and prediction of tides"
  !         by paul schureman, special publication #98, us coast
  !         and geodetic survey, department of commerce (1958).
  
  !-- if daym and hrm correspond to midyear, then this routine
  !-- returns the same values as found in table 15 of schureman.
  !---------------------------------------------------------------------
    real, intent(in) :: yr,dayj,hr,daym,hrm
  
    real :: nu,nup,nup2,i,s,p,h,p1,t,xi,pc,top,bottom,q,r
    integer :: ih
  
  !  common /orbitf/ds,dp,dh,dp1,dn,di,dnu,dxi,dnup,dnup2,dpc
  !  common /cnst/ fndcst(37),eqcst(37),acst(37),pcst(37)
  !  pi180=3.14159265/180.
  !* obtaining orbital values at beginning of series for v0
    call orbit(yr,dayj,hr)
    s=ds_orb
    p=dp_orb
    h=dh_orb
    p1=dp1_orb
    t=angle(180.+hr*(360./24.))
  !** obtaining orbital values at middle of series for u
    call orbit(yr,daym,hrm)
    nu=dnu_orb
    xi=dxi_orb
    nup=dnup_orb
    nup2=dnup2_orb
  !* summing terms to obtain equilibrium arguments
    grterm(1)=2.*(t-s+h)+2.*(xi-nu)
    grterm(2)=2.*t
    grterm(3)=2.*(t+h)-3.*s+p+2.*(xi-nu)
    grterm(4)=t+h-90.-nup
    grterm(5)=4.*(t-s+h)+4.*(xi-nu)
    grterm(6)=t-2.*s+h+90.+2.*xi-nu
    grterm(7)=6.*(t-s+h)+6.*(xi-nu)
    grterm(8)=3.*(t+h)-2.*s-90.+2.*(xi-nu)-nup
    grterm(9)=4.*t
    grterm(10)=4.*(t+h)-5.*s+p+4.*(xi-nu)
    grterm(11)=2.*t-3.*s+4.*h-p+2.*(xi-nu)
    grterm(12)=6.*t
    grterm(13)=2.*(t+2.*(h-s))+2.*(xi-nu)
    grterm(14)=2.*(t-2.*s+h+p)+2.*(xi-nu)
    grterm(15)=t+2.*s+h-90.-2.*xi-nu
    grterm(16)=2.*t-s+p+180.+2.*(xi-nu)
    grterm(17)=t
    i=di_orb*pi180
    pc=dpc_orb*pi180
    top=(5.*cos(i)-1.)*sin(pc)
    bottom=(7.*cos(i)+1.)*cos(pc)
    q=arctan(top,bottom) !,1)
    grterm(18)=t-s+h-90.+xi-nu+q

    !Debug
!    write(99,*)'M1:',grterm(18),t,s,h,xi,nu,q,top,bottom

    grterm(19)=t+s+h-p-90.-nu
    grterm(20)=s-p
    grterm(21)=2.*h
    grterm(22)=h
    grterm(23)=2.*(s-h)
    grterm(24)=2.*s-2.*xi
    grterm(25)=t+3.*(h-s)-p+90.+2.*xi-nu
    grterm(26)=t-3.*s+h+p+90.+2.*xi-nu
    grterm(27)=2.*t-h+p1
    grterm(28)=2.*t+h-p1+180.
    grterm(29)=t-4.*s+h+2.*p+90.+2.*xi-nu
    grterm(30)=t-h+90.
    grterm(31)=2.*(t+s-h)+2.*(nu-xi)
    grterm(32)=3.*(t-s+h)+3.*(xi-nu)
    r=sin(2.*pc)/((1./6.)*(1./tan(.5*i))**2-cos(2.*pc))
    r=atan(r)/pi180
    grterm(33)=2.*(t+h)-s-p+180.+2.*(xi-nu)-r
    grterm(34)=3.*(t+h)-4.*s+90.+4.*(xi-nu)+nup
    grterm(35)=2.*(t+h)-2.*nup2
    grterm(36)=8.*(t-s+h)+8.*(xi-nu)
    grterm(37)=2.*(2.*t-s+h)+2.*(xi-nu)
    do ih=1,ncnst !37
      grterm(ih)=angle(grterm(ih))
    enddo !ih
  end subroutine gterms
  
  subroutine orbit(yr,dayj,hr)
  
  !-- determination of primary and secondary orbital functions
  
  !-- the equations programmed here are not represented by equations in
  !   schureman.  the coding in this routine derives from a program by
  !   the national oceanic and atmospheric administration (noaa).
  !   however, tabular values of the orbital functions can be found in
  !   table 1 of schureman.
  
  
  !---------------------------------------------------------------------
    real, intent(in) :: yr,dayj,hr
  
    real :: i,n,nu,nup,nup2,x,dyr,dday,p,xi
  !  common /orbitf/ds,dp,dh,dp1,dn,di,dnu,dxi,dnup,dnup2,dpc
  
  !  pi180=3.14159265/180.
    x=aint((yr-1901.)/4.)
    dyr=yr-1900.
    dday=dayj+x-1.
  !-- dn is the moon's node (capital n, table 1, schureman)
    dn_orb=259.1560564-19.328185764*dyr-.0529539336*dday-.0022064139*hr
    dn_orb=angle(dn_orb)
    n=dn_orb*pi180
  !-- dp is the lunar perigee (small p, table 1)
    dp_orb=334.3837214+40.66246584*dyr+.111404016*dday+.004641834*hr
    dp_orb=angle(dp_orb)
    p=dp_orb*pi180
    i=acos(.9136949-.0356926*cos(n))
    di_orb=angle(i/pi180)
    nu=asin(.0897056*sin(n)/sin(i))
    dnu_orb=nu/pi180
    xi=n-2.*atan(.64412*tan(n/2.))-nu
    dxi_orb=xi/pi180
    dpc_orb=angle(dp_orb-dxi_orb)
  !-- dh is the mean longitude of the sun (small h, table 1)
    dh_orb=280.1895014-.238724988*dyr+.9856473288*dday+.0410686387*hr
    dh_orb=angle(dh_orb)
  !-- dp1 is the solar perigee (small p1, table 1)
    dp1_orb=281.2208569+.01717836*dyr+.000047064*dday+.000001961*hr
    dp1_orb=angle(dp1_orb)
  !-- ds is the mean longitude of the moon (small s, table 1)
    ds_orb=277.0256206+129.38482032*dyr+13.176396768*dday+.549016532*hr
    ds_orb=angle(ds_orb)
    nup=atan(sin(nu)/(cos(nu)+.334766/sin(2.*i)))
    dnup_orb=nup/pi180
    nup2=atan(sin(2.*nu)/(cos(2.*nu)+.0726184/sin(i)**2))/2.
    dnup2_orb=nup2/pi180
  end subroutine orbit
  
  function angle(arg)
  !
  !*** this routine places an angle in 0-360 (+) format
  !
    real, intent(in) :: arg
    real :: angle
  
    integer :: m
  
    m=-ifix(arg/360.)
    angle=arg+float(m)*360.
    if(angle .lt. 0.) angle=angle+360.
  end function angle
  
  function arctan(top,bottom) !,key)
  !** determine arctangent and place in correct quadrant
  !   if key eq 0  no quadrant selection made
  !   if key .ne. 0 proper quadrant is selected
    real, intent(in) :: top,bottom
!    integer, intent(in) :: key
    real :: arctan
  
!    if(bottom .ne. 0.0) go to 4

    if(bottom==0.) then
!      if(top) 2,9,3
!  2   arctan=270.
!      return
!  3   arctan=90.
!      return
      if(top<0.) then
        arctan=270.
      else if(top==0.) then
        arctan=0.
      else !top>0
        arctan=90.
      endif !top
    else !bottom/=0
      arctan=atan(top/bottom)*57.2957795

!!    if(key.eq.0) return
!      if(top) 5,5,7
!  5   if(bottom) 6,9,8
!  6   arctan=arctan+180.
!      return
!  7   if(bottom) 6,3,10
!  8   arctan=arctan+360.
!      return

      if(top<=0.) then
        if(bottom<0.) then
          arctan=arctan+180.
        else if(bottom==0.) then
          arctan=0.
        else !bottom>0
          arctan=arctan+360.
        endif !bottom
      else !top>0
        if(bottom<0.) then
          arctan=arctan+180.
        else if(bottom==0.) then
          arctan=90.
        else !bottom>0
!          print*, 'WOW'
          return
        endif !bottom
      endif !top
    endif !bottom

!  9 arctan=0.
!  10 return
  end function arctan
  
  function dayjul(yr,xmonth,day)
  !
  !*** this routine computes the julian day (as a real variable)
  !
    real, intent(in) :: yr,xmonth,day
    real :: dayjul
   
    integer :: i
    real :: dayt(12),days(12),dinc,yrlp
  
    dayt=(/0.,31.,59.,90.,120.,151.,181.,212.,243.,273.,304.,334./)
    days(1:2)=(/0.,31./)
    dinc=0.
    yrlp=mod((yr-1900.),4.)
    if(yrlp .eq. 0.) dinc=1.
    do i=3,12
      days(i)=dayt(i)+dinc
    enddo !i
    dayjul=days(ifix(xmonth))+day
  end function dayjul

end module nodalfactor
