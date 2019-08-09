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

	program	combine_outHA

!-----------------------------------------------------------------------
!			    ,
!	AUTHOR:		Andre Fortunato - 2009-07-07
!
!	PURPOSE:	Combine harmonic analysis output files (1 per 
!			processor) into a single file
!
!-----------------------------------------------------------------------

	include		'combine_outHA.cmn'
	real*8		ampl(MXNODL,MXFREQ),phal(MXNODL,MXFREQ)
	real*8		amp(MXNOD,MXFREQ),pha(MXNOD,MXFREQ)
	real*8		freq(MXFREQ),fft(MXFREQ),facet(MXFREQ)
	integer		i,j,k,ne_global,np_global,nvrt,nproc,np,np_local,
     #			nfreq,nodel,nodeg,ne
	character*10	namefr(MXFREQ)
	character*28	file1
	character*20	file2
	character*4	it_char
!-----------------------------------------------------------------------
! Determine number of processors
	file1='outputs/local_to_global_0000'
	open(10,file=file1,status='old')
	read(10,*)ne_global,np_global,nvrt,nproc
	close (10)

! Read files and save in global arrays
	do i = 0, nproc-1
	    write(it_char(1:4),'(i4.4)') i
	    file1='outputs/local_to_global_'//it_char(1:4)
	    open(10,file=file1,status='old')
	    read(10,*)ne_global,np_global,nvrt,nproc
	    if (np_global .gt. MXNOD) then
	    	open(9,file='core',status='unknown')
	    	write(*,*)' np_global exceeds MXNOD'
	    	write(9,*)' np_global exceeds MXNOD'
	    	close(9)
	    	stop
	    endif
	    read(10,*)

	    file2='outputs/harme.53'//it_char(1:4)
	    call r53(ampl,phal,freq,fft,facet,nfreq,np_local,file2,namefr)
! Skip element correspondences
	    write(*,*)'Rank ',i
	    read(10,*)ne
	    write(*,*)'      Elements',ne
	    do j = 1, ne
		read(10,*)
	    end do
! Read nodal correspondences
	    read(10,*)np
	    write(*,*)'      Nodes   ',np
	    if (np .ne. np_local) then
		write(*,*)'Inconsistent number of nodes',i,np,np_local
		stop
	    endif
	    do j = 1, np
		read(10,*)nodel,nodeg
		do k = 1, nfreq
		    amp(nodeg,k) = ampl(nodel,k)
		    pha(nodeg,k) = phal(nodel,k)
		end do
	    end do
	    close (10)
	end do

!  -  Write file
	call w53(amp,pha,freq,fft,facet,nfreq,np_global,namefr)

	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine	r53(amp,pha,freq,fft,facet,nfreq,nodes,file,
     #			    namefr)

c-----------------------------------------------------------------------
c			    ,
c	AUTHOR:		Andre Fortunato - 02-08-05
c
c	PURPOSE:	Read ADCIRC fort.53 files.
c
c-----------------------------------------------------------------------
	include		'combine_outHA.cmn'
	real*8		freq(MXFREQ),fft(MXFREQ),facet(MXFREQ)
	real*8		amp(MXNODL,MXFREQ),pha(MXNODL,MXFREQ)
	character*5	freqname(MXFREQ)
	integer*4	i,j,nodes,nfreq,n
	character*10	namefr(MXFREQ)
	character*20	file
c-----------------------------------------------------------------------
	open(unit=53,file=file,status='old',err=1)
	read(53,*)nfreq
	if (nfreq .gt. MXFREQ) then
	    open(9,file='core',status='unknown')
	    write(*,*)' nfreq exceeds MXFREQ',nfreq
	    write(9,*)' nfreq exceeds MXFREQ',nfreq
	    close(9)
	    stop
	endif
	do i=1,nfreq
	    read(53,3679) freq(i),fft(i),facet(i),namefr(i)
	end do
 3679	format(1x,e20.10,1x,f10.7,1x,f12.8,1x,a10)
	read(53,*)nodes
	if (MXNODL .lt. nodes) then
	    open(9,file='core',status='unknown')
	    write(*,*)' Number of nodes too large.',nodes
	    write(9,*)' Number of nodes too large.',nodes
	    close(9)
	    stop
	endif

	do i = 1, nodes
	    read(53,*)n
	    do j = 1, nfreq
		read(53,*)amp(i,j),pha(i,j)
	    end do
	end do
	close(53)
	return

1	write(*,*)' Sorry, fort.53 is not available.'
	open(9,file='core',status='unknown')
	write(9,*)' Sorry, fort.53 is not available.'
	close(9)
	stop

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine	r54(uamp,upha,vamp,vpha,freq,fft,facet,nfreq,
     #			    nodes,file)

c-----------------------------------------------------------------------
c			    ,
c	AUTHOR:		Andre Fortunato - 02-08-05
c
c	PURPOSE:	Read ADCIRC fort.54 files.
c
c-----------------------------------------------------------------------
	include		'combine_outHA.cmn'
	real*8		freq(MXFREQ),fft(MXFREQ),facet(MXFREQ)
	real*8		uamp(MXNODL,MXFREQ),upha(MXNODL,MXFREQ)
	real*8		vamp(MXNODL,MXFREQ),vpha(MXNODL,MXFREQ)
	character*5	freqname(MXFREQ)
	integer*4	i,j,nodes,nfreq,n
	character*20	file
c-----------------------------------------------------------------------
	open(unit=54,file=file,status='old',err=1)
	read(54,*)nfreq
	if (nfreq .gt. MXFREQ) then
	    open(9,file='core',status='unknown')
	    write(*,*)' nfreq exceeds MXFREQ'
	    write(9,*)' nfreq exceeds MXFREQ'
	    close(9)
	    stop
	endif
	do i = 1, nfreq
	    read(54,*)freq(i),fft(i),facet(i)
	end do
	read(54,*)nodes
	if (MXNODL .lt. nodes) then
	    open(9,file='core',status='unknown')
	    write(*,*)' Number of nodes too large.',nodes
	    write(9,*)' Number of nodes too large.',nodes
	    close(9)
	    stop
	endif

	do i = 1, nodes
	    read(54,*)n
	    do j = 1, nfreq
		read(54,*)uamp(i,j),upha(i,j),vamp(i,j),vpha(i,j)
	    end do
	end do
	close(54)
	return

1	write(*,*)' Sorry, fort.54 is not available.'
	open(9,file='core',status='unknown')
	write(9,*)' Sorry, fort.54 is not available.'
	close(9)
	stop

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine	w53(amp,pha,freq,fft,facet,nfreq,nodes,namefr)

c-----------------------------------------------------------------------
c			    ,
c	AUTHOR:		Andre Fortunato - 02-08-05
c
c	PURPOSE:	Write ADCIRC fort.53 files.
c
c-----------------------------------------------------------------------
	include		'combine_outHA.cmn'
	real*8		freq(MXFREQ),fft(MXFREQ),facet(MXFREQ)
	real*8		amp(MXNOD,MXFREQ),pha(MXNOD,MXFREQ)
	character*10	namefr(MXFREQ)
	integer*4	i,j,nodes,nfreq,n
c-----------------------------------------------------------------------
	open(unit=53,file='fort.53',status='unknown')
	write(53,*)nfreq
	do i = 1, nfreq
	    write(53,3679)freq(i),fft(i),facet(i),namefr(i)
	end do
3679	format(1x,e20.10,1x,f10.7,1x,f12.8,1x,a10)
	write(53,*)nodes

	do i = 1, nodes
	    write(53,*)i
	    do j = 1, nfreq
		write(53,6635)amp(i,j),pha(i,j)
	    end do
	end do
6635	format(2x,e16.8,1x,f11.4)
	close(53)
	return

	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine	w54(uamp,upha,vamp,vpha,freq,fft,facet,nfreq,
     #			    nodes)

c-----------------------------------------------------------------------
c			    ,
c	AUTHOR:		Andre Fortunato - 02-08-05
c
c	PURPOSE:	Write ADCIRC fort.54 files.
c
c-----------------------------------------------------------------------
	include		'combine_outHA.cmn'
	real*8		freq(MXFREQ),fft(MXFREQ),facet(MXFREQ)
	real*8		uamp(MXNOD,MXFREQ),upha(MXNOD,MXFREQ)
	real*8		vamp(MXNOD,MXFREQ),vpha(MXNOD,MXFREQ)
	character*5	freqname(MXFREQ)
	integer*4	i,j,nodes,nfreq,n
c-----------------------------------------------------------------------
	open(unit=54,file='fort.54',status='unknown')
	write(54,*)nfreq
	do i = 1, nfreq
	    write(54,3679)freq(i),fft(i),facet(i),freqname(i)
	end do
3679	format(1x,e20.10,1x,f10.7,1x,f12.8,1x,a5)
	write(54,*)nodes

	do i = 1, nodes
	    write(54,*)i
	    do j = 1, nfreq
		write(54,6636)uamp(i,j),upha(i,j),vamp(i,j),vpha(i,j)
	    end do
	end do
6636	format(2x,e16.8,1x,f11.4,2x,e16.8,1x,f11.4)
	close(54)
	return

	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
