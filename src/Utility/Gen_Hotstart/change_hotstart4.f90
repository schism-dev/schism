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
! Manipulate hotstart.in: only vgrid (pure S) is
! allowed to change.

! Inputs:
!        hotstart.old (unformatted binary); screen inputs.
! Output: hotstart.in (unformatted binary). This format is different
!         between Intel and AMD!
!
!  ifort -Bstatic -O2 -assume byterecl -o change_hotstart4 change_hotstart4.f90

!===============================================================================

program combine_hotstart1
!-------------------------------------------------------------------------------
  implicit real(8)(a-h,o-z),integer(i-n)
  parameter(nbyte=4)
!  character(12) :: it_char
!  character(72) :: fgb,fgb2,fdb  ! Processor specific global output file name
!  integer :: lfgb,lfdb       ! Length of processor specific global output file name
!  allocatable ner(:),npr(:),nsr(:)
!  allocatable ielg(:,:),iplg(:,:),islg(:,:)
  allocatable idry_e(:),we(:,:),tsel(:,:,:),idry_s(:),su2(:,:),sv2(:,:)
  allocatable tsd(:,:),ssd(:,:),idry(:),eta2(:),tnd(:,:),snd(:,:)
  allocatable tem0(:,:),sal0(:,:),q2(:,:),xl(:,:),dfv(:,:),dfh(:,:)
  allocatable dfq1(:,:),dfq2(:,:),qnon(:,:),trnd(:,:,:),trel(:,:,:),trnd0(:,:,:)
  allocatable intv(:),zrat(:),swild(:,:)
!-------------------------------------------------------------------------------
      
  print*, 'Input ne, np, ns, nvrt, and ntracers (including modules):'
  read*, ne_global,np_global,ns_global,nvrt,ntracers
  if(ntracers<2) stop 'wrong ntracers'
  print*, 'Input new # of levels:'
  read*, nvrt1

  allocate(idry_e(ne_global),we(nvrt,ne_global), &
           idry_s(ns_global),su2(nvrt,ns_global),sv2(nvrt,ns_global), &
           tsd(nvrt,ns_global),ssd(nvrt,ns_global), &
           idry(np_global),eta2(np_global),trnd(ntracers,nvrt,np_global), &
           trnd0(ntracers,nvrt,np_global),trel(ntracers,nvrt,ne_global),q2(np_global,nvrt), &
           xl(np_global,nvrt),dfv(np_global,nvrt),dfh(np_global,nvrt), &
           dfq1(np_global,nvrt),dfq2(np_global,nvrt),qnon(nvrt,np_global), &
           swild(nvrt1,7+2*ntracers),intv(nvrt1),zrat(nvrt1),stat=istat)
  if(istat/=0) stop 'Allocation error (2)'

! Calculate interplation coefficients, assuming equal distance in vertical
  do k=1,nvrt1
   sig=-1+(k-1.)/(nvrt1-1)
   tmp=(sig+1)*(nvrt-1)+1
   intv(k)=int(tmp)
   if(intv(k)==nvrt) intv(k)=nvrt-1
   zrat(k)=tmp-intv(k) !from level intv(k)
  enddo !k

!-------------------------------------------------------------------------------
! Read hotstart files
!-------------------------------------------------------------------------------
  open(36,file='hotstart.old',form='unformatted',status='old')
  open(37,file='hotstart.in',form='unformatted',status='replace')
  read(36) time,it,ifile
  write(37) time,it,ifile
  do i=1,ne_global
    read(36) j,idry_e(i),(we(j,i),(trel(l,j,i),l=1,ntracers),j=1,nvrt)
    do j=1,nvrt1
      k=intv(j)
      swild(j,1)=we(k,i)*(1-zrat(k))+we(k+1,i)*zrat(k)
      swild(j,2:ntracers+1)=trel(1:ntracers,k,i)*(1-zrat(k))+trel(1:ntracers,k+1,i)*zrat(k)
    enddo !j
    write(37) i,idry_e(i),(swild(j,1:ntracers+1),j=1,nvrt1)
  enddo !i
  do i=1,ns_global
    read(36) j,idry_s(i),(su2(j,i),sv2(j,i),tsd(j,i),ssd(j,i),j=1,nvrt)
    do j=1,nvrt1
      k=intv(j)
      swild(j,1)=su2(k,i)*(1-zrat(k))+su2(k+1,i)*zrat(k)
      swild(j,2)=sv2(k,i)*(1-zrat(k))+sv2(k+1,i)*zrat(k)
      swild(j,3)=tsd(k,i)*(1-zrat(k))+tsd(k+1,i)*zrat(k)
      swild(j,4)=ssd(k,i)*(1-zrat(k))+ssd(k+1,i)*zrat(k)
    enddo !j
    write(37) i,idry_s(i),(swild(j,1:4),j=1,nvrt1)
  enddo !i
  do i=1,np_global
    read(36) j,eta2(i),idry(i),(trnd(:,j,i),trnd0(:,j,i),q2(i,j),xl(i,j), &
             dfv(i,j),dfh(i,j),dfq1(i,j),dfq2(i,j),qnon(j,i),j=1,nvrt)
    do j=1,nvrt1
      k=intv(j)
      swild(j,1:ntracers)=trnd(:,k,i)*(1-zrat(k))+trnd(:,k+1,i)*zrat(k)
      swild(j,ntracers+1:2*ntracers)=trnd0(:,k,i)*(1-zrat(k))+trnd0(:,k+1,i)*zrat(k)
      swild(j,2*ntracers+1)=q2(i,k)*(1-zrat(k))+q2(i,k+1)*zrat(k)
      swild(j,2*ntracers+2)=xl(i,k)*(1-zrat(k))+xl(i,k+1)*zrat(k)
      swild(j,2*ntracers+3)=dfv(i,k)*(1-zrat(k))+dfv(i,k+1)*zrat(k)
      swild(j,2*ntracers+4)=dfh(i,k)*(1-zrat(k))+dfh(i,k+1)*zrat(k)
      swild(j,2*ntracers+5)=dfq1(i,k)*(1-zrat(k))+dfq1(i,k+1)*zrat(k)
      swild(j,2*ntracers+6)=dfq2(i,k)*(1-zrat(k))+dfq2(i,k+1)*zrat(k)
      swild(j,2*ntracers+7)=qnon(k,i)*(1-zrat(k))+qnon(k+1,i)*zrat(k)
    enddo !j
    write(37) i,eta2(i),idry(i),(swild(j,1:2*ntracers+7),j=1,nvrt1)
  enddo !i

! Other modules (SED etc)

  close(36)
  close(37)

end program combine_hotstart1
