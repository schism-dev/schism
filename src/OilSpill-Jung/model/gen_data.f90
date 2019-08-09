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

!ifort -O3 -assume byterecl -Bstatic -o gen_data gen_data.f90
      program gen_data
      implicit real(kind=8)(a-h,o-z),integer(i-n)
      character(len=60), dimension(8) :: str
! ................................................
! ... spill mode
      ispl=1     ! simultaneous spill  
!      ispl=2    ! contiuous spill      

!...  open output file
      if(ispl==1) then
         open(91,file='data/spill_Wb.in',status='unknown')
        else
         open(91,file='data/spill_Ws.in',status='unknown')
      endif

!...  print basic parameters
      str(1)='Oil spill simulation in Columbia river'                                        
      str(2)='1  !nscreen'                                              
      str(3)='0  !istiff'                                          
      str(4)='2 -124.00 46.25       !ics,slam0,sfea0'                       
      str(5)='0.1 3.0 90 10 960 10  !h0,rnday,dtm,nspool,ihfskip,ndeltp'
      str(6)='1  3.0  0.2           !ihdf,hdc,horcon'
      str(7)='1  0                  !ibuoy,iwind'
      str(8)='20.0                  !pbeach'   
      write(91,'(a60)') str(1)
      write(91,'(a60)') str(2)
      write(91,'(a60)') str(3)
      write(91,'(a60)') str(4)
      write(91,'(a60)') str(5)
      write(91,'(a60)') str(6)
      write(91,'(a60)') str(7)     
      write(91,'(a60)') str(8)      

!...  print particle data
! ... Entrance
!      xp=-124.066; yp=46.2417; zp=-7.0 !sub-surface, -7.0
!      xp=-124.066; yp=46.2417; zp=0.0

! ... Tongue Point
!      xp=-123.759; yp=46.224; zp=-8.0  !sub-surface, -8.0
!      xp=-123.759; yp=46.224; zp=0.0  

! ... Wauna
      xp=-123.39; yp=46.1601; zp=-6.0 !sub-surface, -6.0
!      xp=-123.39; yp=46.1601; zp=0.0

! ... release particles 
      if(ispl.eq.1) then        ! simultaneous
        is=0                    ! time to start releasing
        nparticle=500           ! number of particles to be released
        write(91,'(i6)') nparticle
        do i=1,nparticle
!           write(91,'(i5,2x,i5,2f12.2,f8.2)') i,is,xp,yp,zp        ! ics=1
           write(91,'(i5,2x,i5,2x,2f10.4,f8.2)') i,is,xp,yp,zp    ! ics=2
        enddo
!
       else   ! continuance
        ie=24*60*60  ! spill duration(sec), 24 hours
        im=90        ! del_time 
        nt=10        ! number of particles to be released per time step
        in=0
        nparticle=ie/im*nt 
        write(91,'(i6)') nparticle
        do i=1,ie,im
           is=i
           if(i.ne.1) is=is-1
          do j=1,nt
             in=in+1
!           write(91,'(i5,2x,i5,2f12.2,f8.2)') in,is,xp,yp,zp        ! ics=1
           write(91,'(i5,2x,i5,2x,2f10.4,f8.2)') in,is,xp,yp,zp    ! ics=2   
          enddo !j
        enddo ! i
      endif
! ...
      close(91) 
      stop
      end
