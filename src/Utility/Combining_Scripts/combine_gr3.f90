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

!===============================================================================
! Read in *.gr3-like (rank-specific) outputs from SCHISM and combine them into one global output
! e.g. maxelev.gr3; may have multiple scalar fields
! Works with mixed hgrid

! Inputs:
!        (0) screen: filenm - name of file (e.g. maxelev);
!        (1) hgrid.gr3;
!        (2) outputs/<filenm>_0*
! Output: <filenm>.gr3
!
!  Compile on amb64xx:
!  ifort -Bstatic -O3 -assume byterecl -o combine_gr3 combine_gr3.f90
!  PGI compiler:
!  pgf90 -O2 -mcmodel=medium  -Bstatic -o combine_gr3 combine_gr3.f90
!===============================================================================

program combine_gr3
!-------------------------------------------------------------------------------

  implicit real(8)(a-h,o-z),integer(i-n)
  character(36) :: fdb,filenm 
  integer :: lfdb,lfilenm,nm(4)
  allocatable x(:),y(:),elevmax(:,:)
      
!-------------------------------------------------------------------------------
! inputs
!-------------------------------------------------------------------------------

  print*, 'Input file name (e.g.: maxelev):'
  read*, filenm
  filenm=adjustl(filenm); lfilenm=len_trim(filenm)
  print*, 'Input # of scalar fields:'
  read*, nscal
  if(nscal<=0) stop 'Wrong nscal'

  open(14,file='hgrid.gr3',status='old')
  read(14,*); read(14,*)ne,np
  allocate(x(np),y(np),elevmax(nscal,np),stat=istat)
  if(istat/=0) stop 'Allocation error: x,y'

  open(10,file='outputs/'//filenm(1:lfilenm)//'_0000',status='old')
  read(10,*)icount,nproc
  close(10)

!-------------------------------------------------------------------------------
! Combine
!-------------------------------------------------------------------------------
  fdb=filenm(1:lfilenm)//'_0000'
  lfdb=len_trim(fdb)

  do irank=0,nproc-1
    write(fdb(lfdb-3:lfdb),'(i4.4)') irank
    open(10,file='outputs/'//fdb,status='old')
    read(10,*)icount !,nproc
    do i=1,icount
      read(10,*)nd,xtmp,ytmp,elevmax(:,nd)
    enddo !i
  enddo !irank

  open(13,file=filenm(1:lfilenm)//'.gr3',status='replace')
  write(13,*); write(13,*)ne,np
  do i=1,np
    read(14,*)j,xtmp,ytmp
    write(13,'(i10,100(1x,e22.11))')i,xtmp,ytmp,elevmax(:,i)
  enddo !i
  do i=1,ne
    read(14,*)j,k,nm(1:k)
    write(13,*)j,k,nm(1:k)
  enddo !i

end program combine_gr3
