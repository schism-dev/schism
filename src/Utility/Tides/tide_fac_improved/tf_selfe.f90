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

subroutine selecttides(numcon,ncon)
  implicit none

  interface
    elemental subroutine upper_case(word)
    ! convert a word to lower case 
    character (len=*) , intent(in out) :: word
    integer :: i,ic,nlen
    end subroutine upper_case
  end interface

  integer,parameter :: ncnst=37
  real,parameter :: pi=3.1415926
  integer :: numcon,ncon,nbfr,nc,ntip
  integer :: i,istat,ncst
  real :: tip_dp
  character(48) :: start_time
  ! Boundary forcings
  logical :: lbc   !local dummy
  logical :: lerbc !flag to indicate a partitioned radiation boundary segment
  logical :: lflbc !flag to indicate existence of ifltype/=0
  logical :: ltobc !flag to indicate a patitioned temperature open boundary segment
  logical :: lsobc !flag to indicate a patitioned salinity open boundary segment
  logical,allocatable :: lelbc(:)
  integer,allocatable :: iettype(:),ifltype(:),itetype(:),isatype(:),itrtype(:),trobc(:), &
    &tobc(:),sobc(:),vobc1(:),vobc2(:)
  real,allocatable :: tamp(:),tnf(:),tfreq(:),jspc(:),tear(:)
  real,allocatable :: amig(:),ff(:),face(:)
  real,allocatable :: emo(:,:,:),efa(:,:,:),vmo(:,:),vfa(:,:)
  real,allocatable :: eth(:,:),qthcon(:),ath(:)
  real(4), allocatable :: a2th(:,:,:)
  real,allocatable :: uth(:,:),vth(:,:),uthnd(:,:,:),vthnd(:,:,:)
  real,allocatable :: carea(:),z_r2(:),eta_mean(:)
  character(48),allocatable :: tname(:),fname(:)
  character(8) cname(ncnst),ctmp
  common /cnsnam/ cname
  dimension ncon(200)
!      allocatable :: iet1lg(:),ifl1lg(:),ite1lg(:),isa1lg(:)

  ncon(:)=0
!-------------------------------------------------------------------------------
! Read in boundary condition and tidal info
!-------------------------------------------------------------------------------
  open(31,file='bctides.in',status='old')
  read(31,'(a48)') start_time
  !...  Earth tidal potential
  read(31,*) ntip,tip_dp !cut-off depth for applying tidal potential
  if(ntip>0) then
    allocate(tname(ntip),tamp(ntip),tnf(ntip),tfreq(ntip),jspc(ntip),tear(ntip),stat=istat)
    if(istat/=0) then
      write(*,*)('Allocation failure for tname tamp')
    STOP
  endif
  !'
    do i=1,ntip
      read(31,'(a48)') tname(i) !tag
      read(31,*) jspc(i),tamp(i),tfreq(i),tnf(i),tear(i)
      if(jspc(i)<0.or.jspc(i)>2) then
        write(*,*)'Illegal tidal species #',jspc(i)
        STOP
      endif
      tear(i)=tear(i)*pi/180
    enddo !i
  endif !ntip>0

  !...  Boundary forcing freqs.
  !     All b.c. arrays are global
  read(31,*) nbfr
  allocate(fname(nbfr),stat=istat)
  if(istat/=0) then
    write(*,*)('Allocation failure for bname')
    STOP
  endif
  if(nbfr>0) then
    allocate(amig(nbfr),ff(nbfr),face(nbfr),stat=istat)
    if(istat/=0) then
      write(*,*)'Allocation failure for amig etc'
      STOP
  endif
  !'
    do i=1,nbfr
      read(31,*) fname(i) !tag
      read(31,*) amig(i),ff(i),face(i) !freq., nodal factor and earth equil.
      face(i)=face(i)*pi/180
    enddo
  endif
  close(31)

  do nc=1,ntip
    call upper_case(tname(nc))
    do ncst=1,37
      ctmp=cname(ncst)
      call upper_case(ctmp)
      if(trim(adjustl(tname(nc))).eq.trim(adjustl(ctmp)))then
        ncon(nc)=ncst
      endif
    enddo
    if(ncon(nc).eq.0.and.(trim(tname(nc)).ne.'Z0'))then
      write(*,"('no match earth tidal potential constituent for: ',a5,'(',i3,'th of ',i3,' )')")trim(tname(nc)),nc,ntip
      stop
    else
      write(*,"('found match earth tidal potential constituent for: ',a5,'(',i3,'th of ',i3,' )')")trim(tname(nc)),nc,ntip
    endif
  enddo
  do nc=1,nbfr
    call upper_case(fname(nc))
    ncon(nc+ntip)=0
    do ncst=1,37
      ctmp=cname(ncst)
      call upper_case(ctmp)
      if(trim(adjustl(fname(nc))).eq.trim(adjustl(ctmp)))then
        ncon(nc+ntip)=ncst
    endif
    enddo
    if(ncon(nc+ntip).eq.0.and.(trim(fname(nc)).ne.'Z0'))then
      write(*,"('no match constituent for: ',a5,'(',i3,'th of ',i3,' )')")trim(fname(nc)),nc,nbfr
      stop
    else
      if(trim(fname(nc)).ne.'Z0')then
        write(*,"('found match constituent ',a5,' ( #',i3,' ) for: ',a5,'(',i3,'th of ',i3,' )')")&
          &trim(cname(ncon(nc+ntip))),ncon(nc+ntip),trim(fname(nc)),nc,nbfr
      else
        write(*,"('found match constituent ',a5,' ( #',i3,' ) for: ',a5,'(',i3,'th of ',i3,' )')")&
          &'Z0',ncon(nc+ntip),trim(fname(nc)),nc,nbfr
      endif
    endif
  enddo
  numcon=ntip+nbfr

end subroutine selecttides


subroutine writebctides(ncon,bhr,iday,imo,iyr)
  implicit none

  interface
    elemental subroutine upper_case(word)
    ! convert a word to lower case 
    character (len=*) , intent(in out) :: word
    integer :: i,ic,nlen
    end subroutine upper_case
  end interface

  integer,parameter :: ncnst=37
  real,parameter :: pi=3.1415926
  integer :: ncon,ntip
  integer :: i,nbfr,io,istat
  ! Boundary forcings
  logical :: lbc   !local dummy
  logical :: lerbc !flag to indicate a partitioned radiation boundary segment
  logical :: lflbc !flag to indicate existence of ifltype/=0
  logical :: ltobc !flag to indicate a patitioned temperature open boundary segment
  logical :: lsobc !flag to indicate a patitioned salinity open boundary segment
  real :: nodfac,grterm,speed,p,bhr
  real :: tip_dp
  character(48) :: start_time
  logical,allocatable :: lelbc(:)
  integer,allocatable :: iettype(:),ifltype(:),itetype(:),isatype(:),itrtype(:),trobc(:), &
    &tobc(:),sobc(:),vobc1(:),vobc2(:),jspc(:)
  real,allocatable :: tamp(:),tnf(:),tfreq(:),tear(:)
  real,allocatable :: amig(:),ff(:),face(:)
  real,allocatable :: emo(:,:,:),efa(:,:,:),vmo(:,:),vfa(:,:)
  real,allocatable :: eth(:,:),qthcon(:),ath(:)
  real(4), allocatable :: a2th(:,:,:)
  real,allocatable :: uth(:,:),vth(:,:),uthnd(:,:,:),vthnd(:,:,:)
  real,allocatable :: carea(:),z_r2(:),eta_mean(:)
  character(48),allocatable :: tname(:),fname(:)
  character(48) :: tmp
  integer :: iyr,imo,iday,ihr,imn,isc
  character(8) cname(ncnst),ctmp
  common /cnsnam/ cname
  common /cnst/ nodfac(ncnst),grterm(ncnst),speed(ncnst),p(ncnst)
  dimension ncon(200)
  !allocatable :: iet1lg(:),ifl1lg(:),ite1lg(:),isa1lg(:)

!-------------------------------------------------------------------------------
! Read in boundary condition and tidal info
!-------------------------------------------------------------------------------
  open(31,file='bctides.in',status='old')
  read(31,'(a48)') start_time
  !...  Earth tidal potential
  read(31,*) ntip,tip_dp !cut-off depth for applying tidal potential
  if(ntip>0) then
    allocate(tname(ntip),tamp(ntip),tnf(ntip),tfreq(ntip),jspc(ntip),tear(ntip),stat=istat)
    if(istat/=0) then
      write(*,*)('Allocation failure for tname tamp')
    STOP
  endif
  !'
    do i=1,ntip
      read(31,"(a48)") tname(i) !tag
      read(31,*) jspc(i),tamp(i),tfreq(i),tnf(i),tear(i)
      if(jspc(i)<0.or.jspc(i)>2) then
        write(*,*)'Illegal tidal species #',jspc(i)
        STOP
      endif
    enddo !i
  endif !ntip>0

  !...  Boundary forcing freqs.
  !     All b.c. arrays are global
  read(31,*) nbfr
  allocate(fname(nbfr),stat=istat)
  if(istat/=0) then
    write(*,*)('Allocation failure for bname')
    STOP
  endif
  if(nbfr>0) then
    allocate(amig(nbfr),ff(nbfr),face(nbfr),stat=istat)
    if(istat/=0) then
      write(*,*)'Allocation failure for amig etc'
      STOP
    endif
  !'
    do i=1,nbfr
      read(31,"(a48)") fname(i) !tag
      read(31,*) amig(i),ff(i),face(i) !freq., nodal factor and earth equil.
    enddo
  endif


  open(32,file='bctides.in.out',status='replace')
  !       write(32,'(a48)') start_time
  ihr=mod(bhr,24.0)
  imn=mod(bhr*60,60.0)
  isc=mod(bhr*3600,60.0)
  write(32,"(i2.2,'/',i2.2,'/',i4.4,1x,i2.2,':',i2.2,':',i2.2,1x,'UTC')") imo,iday,iyr,ihr,imn,isc
  !...  Earth tidal potential
  write(32,"(i2,1x,f7.3,1x,'!number of earth tidal potential, cut-off depth for applying &
    &tidal potential')") ntip,tip_dp !cut-off depth for applying tidal potential
  !'
  do i=1,ntip
    write(32,"(a48)") tname(i) !tag
    call upper_case(tname(i))
    if(trim(tname(i)).eq.'Z0')then
      write(32,"(i1,1x,f9.6,1x,e14.6,1x,f8.5,1x,f9.5)") 0, 1., 0., 1., 0.
    else
      write(32,"(i1,1x,f9.6,1x,e14.6,1x,f8.5,1x,f9.5)") jspc(i),tamp(i),tfreq(i),NODFAC(NCON(i)),GRTERM(NCON(i))
    endif
    if(jspc(i)<0.or.jspc(i)>2) then
      write(*,*)'Illegal tidal species #',jspc(i)
      STOP
    endif
  enddo !i

  !...  Boundary forcing freqs.
  !     All b.c. arrays are global
  write(32,"(i2,1x,'!number of boundary forcing freqs')") nbfr
  do i=1,nbfr
    write(32,"(a48)") fname(i) !tag
    call upper_case(fname(i))
    if(trim(fname(i)).eq.'Z0')then
      write(32,"(1x,e14.6,1x,f8.5,1x,f9.5)") amig(i),1.,0.
    else
      write(32,"(1x,e14.6,1x,f8.5,1x,f9.5)") amig(i),NODFAC(NCON(i+ntip)),GRTERM(NCON(i+ntip)) !freq., nodal factor and earth equil.
    endif
  enddo
  
  read(31,"(a48)",iostat=io)tmp
  do while (io.eq.0)
    write(32,"(a48)")tmp
    read(31,"(a48)",iostat=io)tmp
  enddo
  close(31)
  close(32)

end subroutine writebctides



elemental subroutine lower_case(word) 
! convert a word to lower case 
character (len=*) , intent(in out) :: word 
integer :: i,ic,nlen 
nlen = len(word) 
do i=1,nlen 
  ic = ichar(word(i:i)) 
  if (ic >= 65 .and. ic < 90) word(i:i) = char(ic+32) 
end do 
end subroutine lower_case 


elemental subroutine upper_case(word) 
! convert a word to lower case 
character (len=*) , intent(in out) :: word 
integer :: i,ic,nlen 
nlen = len(word) 
do i=1,nlen 
  ic = ichar(word(i:i)) 
  if (ic >= 97 .and. ic < 122) word(i:i) = char(ic-32) 
end do 
end subroutine upper_case 
