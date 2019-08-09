! Calc fraction for last class of sediment (to ensure that sum=1)
! May modify fractions if sum>1 (scaled down all classes proportionally)
! Inputs: screen; bed_frac_X.ic.0 where X=1,..,nsed-1
! Output: bed_frac_<1...nsed>.ic
! ifort -Bstatic -O2 -CB -mcmodel=medium -assume byterecl -o last_fraction last_fraction.f90  
  
  character(len=30) :: char2
  character(len=9) :: fname
  allocatable :: xnd(:),ynd(:),sum1(:),frac(:,:),i34(:),ielnode(:,:)

  print*, 'Input total # of classes:'
  read*, nsed

  fname='bed_frac_'
  do i=1,nsed-1
    write(char2,'(i10)')i
    char2=adjustl(char2); lchar=len_trim(char2)
    !print*, 'chars= ', char2,lchar,char2(1:lchar)
    !print*, 'input file name= ',fname//char2(1:lchar)//'.ic'
     
    open(10,file=fname//char2(1:lchar)//'.ic.0',status='old')
    read(10,*); read(10,*)ne,np
    if(i==1) then
      allocate(xnd(np),ynd(np),sum1(np),i34(ne),ielnode(4,ne),frac(np,nsed))
      sum1=0
    endif
    do m=1,np
      read(10,*)j,xnd(m),ynd(m),frac(m,i)
      if(frac(m,i)<0) stop 'frac(m,i)<0'
      sum1(m)=sum1(m)+frac(m,i)
    enddo !m
    do m=1,ne
      read(10,*)j,i34(m),ielnode(1:i34(m),m)
    enddo !m
    close(10)
  enddo !i=1,nsed-1

! Check sum
  do m=1,np
    if(sum1(m)>1) then
      do i=1,nsed-1
        frac(m,i)=frac(m,i)/sum1(m)
      enddo !i
      frac(m,nsed)=0
    else
      frac(m,nsed)=1-sum1(m)
    endif
  enddo !m

! Output
  do i=1,nsed
    write(char2,'(i10)')i
    char2=adjustl(char2); lchar=len_trim(char2)

    open(9,file=fname//char2(1:lchar)//'.ic',status='replace')
    write(9,*)'class # ',i
    write(9,*)ne,np
    do m=1,np
      write(9,*)m,xnd(m),ynd(m),frac(m,i)
    enddo !m
    do m=1,ne
      write(9,*)m,i34(m),ielnode(1:i34(m),m)
    enddo !m
    close(9)
  enddo !i=1,nsed

  stop
  end
