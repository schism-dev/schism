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


!****************************************************************************************
!			
!	Extract binary results in a part of hgrid.gr3
!       Works for mixed tri/quads outputs.
!       Inputs: binary file and hgrid.gr3; subflag.gr3 (hgrid.gr3 based)
!       Outputs: ?_sub_<file63> (new binary file); fort.13 (output grid
!       for mlab)
!
!       ifort -Bstatic -assume byterecl -O3 -o extract_subregion2 extract_subregion2.f90
!       pgf90 -O2 -mcmodel=medium -Mbounds -o extract_subregion2 extract_subregion2.f90
!		
!   History: use subflag.gr3 instead of subgrid.gr3 so much faster than extract_subregion
!
      program read_out
      parameter(nbyte=4)
      character*30 file63
      character*12 it_char
      character*48 start_time,version,variable_nm,variable_dim
      character*48 data_format
      allocatable :: sigma(:),cs(:),ztot(:),x(:),y(:),dp(:),kbp00(:)
      integer, allocatable :: i34(:),elnode(:,:),i34out(:),nmout(:,:)
      allocatable :: icum(:,:,:),eta2(:),imapelem2(:)
      allocatable :: xout(:),yout(:),dpout(:),imapnode(:),imapnode2(:),ifl(:)
      !long int for large files
      integer(kind=8) :: irec,irecout
      
      print*, 'Input file to read from (without *_):'
      read(*,'(a30)')file63
      
      print*, 'Input start and end file #:'
      read(*,*) istart,iend

!     Read hgrid.gr3
      open(14,file='hgrid.gr3',status='old')
      read(14,*)
      read(14,*)ne,np
      allocate(x(np),y(np),dp(np),kbp00(np),i34(ne),elnode(4,ne),i34out(ne),nmout(4,ne), &
     &eta2(np),xout(np),yout(np),dpout(np),imapnode(np),imapnode2(np),ifl(np),imapelem2(ne))
      do i=1,np
        read(14,*)j,x(i),y(i),dp(i)
      enddo !i
      do i=1,ne
        read(14,*)j,i34(i),elnode(1:i34(i),i)
      enddo !i
      close(14)

!     Read subflag.gr3
      open(14,file='subflag.gr3',status='old')
      read(14,*)
      read(14,*)ne2,np2
      if(np2/=np) stop 'subflag.gr3 should be based on hgrid.gr3'
!'
      do i=1,np
        read(14,*)j,xtmp,ytmp,tmp
        ifl(i)=tmp
      enddo !i
      close(14)

!...  Find nodes in the original grid; create mapping
      imapnode=0 !from new to original grid
      imapnode2=0 !from original to new grid
      imapelem2=0 !from original to new grid
      neout=0
      npout=0
      do i=1,ne
        itmp=minval(ifl(elnode(1:i34(i),i)))
        if(itmp/=0) then
          neout=neout+1
          imapelem2(i)=neout
          i34out(neout)=i34(i)
          do j=1,i34(i)
            nd=elnode(j,i)
            if(imapnode2(nd)==0) then !new node
              npout=npout+1
              imapnode2(nd)=npout
              imapnode(npout)=nd
              xout(npout)=x(nd)
              yout(npout)=y(nd)
              dpout(npout)=dp(nd)
            endif !imapnode2
          enddo !j
          nmout(1:i34out(neout),neout)=imapnode2(elnode(1:i34(i),i))
        endif !itmp
      enddo !i=1,ne

!     Output grid
      write(13,*)
      write(13,*)neout,npout
      do i=1,npout
        write(13,*)i,xout(i),yout(i),dpout(i)
      enddo !i
      do i=1,neout
        write(13,*)i,i34out(i),nmout(1:i34out(i),i)
      enddo !i
      close(13)

      print*, '# of nodes in output grid= ',npout
      print*, 'Starting time iteration...'

!...  Time iteration
!...
      it_tot=(istart-1)*nrec
      do iday=istart,iend
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
!...  Header
!...
      write(it_char,'(i12)')iday
      open(63,file=it_char//'_'//file63,status='old',access='direct',recl=nbyte)
      open(65,file=it_char//'_sub_'//file63,access='direct',recl=nbyte)

      print*, 'Doing ',it_char//'_'//file63

      irec=0
      irecout=0
      do m=1,48/nbyte
        read(63,rec=irec+m) data_format(nbyte*(m-1)+1:nbyte*m)
        write(65,rec=irecout+m) data_format(nbyte*(m-1)+1:nbyte*m)
      enddo
      if(data_format.ne.'DataFormat v5.0') then
        print*, 'This code reads only v5.0:  ',data_format
        stop
      endif
      irec=irec+48/nbyte
      irecout=irecout+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) version(nbyte*(m-1)+1:nbyte*m)
        write(65,rec=irecout+m) version(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte
      irecout=irecout+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) start_time(nbyte*(m-1)+1:nbyte*m)
        write(65,rec=irecout+m) start_time(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte
      irecout=irecout+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) variable_nm(nbyte*(m-1)+1:nbyte*m)
        write(65,rec=irecout+m) variable_nm(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte
      irecout=irecout+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) variable_dim(nbyte*(m-1)+1:nbyte*m)
        write(65,rec=irecout+m) variable_dim(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte
      irecout=irecout+48/nbyte

!      write(*,'(a48)')data_format
!      write(*,'(a48)')version
!      write(*,'(a48)')start_time
!      write(*,'(a48)')variable_nm
!      write(*,'(a48)')variable_dim

      read(63,rec=irec+1) nrec
      read(63,rec=irec+2) dtout
      read(63,rec=irec+3) nspool
      read(63,rec=irec+4) ivs
      read(63,rec=irec+5) i23d
      irec=irec+5
      write(65,rec=irecout+1) nrec
      write(65,rec=irecout+2) dtout
      write(65,rec=irecout+3) nspool
      write(65,rec=irecout+4) ivs
      write(65,rec=irecout+5) i23d
      irecout=irecout+5

!      print*, 'i23d=',i23d,' nrec= ',nrec

!     Vertical grid
      read(63,rec=irec+1) nvrt
      read(63,rec=irec+2) kz
      read(63,rec=irec+3) h0
      read(63,rec=irec+4) h_s
      read(63,rec=irec+5) h_c
      read(63,rec=irec+6) theta_b
      read(63,rec=irec+7) theta_f
      irec=irec+7
      write(65,rec=irecout+1) nvrt
      write(65,rec=irecout+2) kz
      write(65,rec=irecout+3) h0
      write(65,rec=irecout+4) h_s
      write(65,rec=irecout+5) h_c
      write(65,rec=irecout+6) theta_b
      write(65,rec=irecout+7) theta_f
      irecout=irecout+7

      if(iday==istart) then
        allocate(sigma(nvrt),cs(nvrt),ztot(nvrt),icum(np,nvrt,2))
      endif
!      print*, 'nvrt=',nvrt,' theta_f= ',theta_f

      do k=1,kz-1
        read(63,rec=irec+k) ztot(k)
        write(65,rec=irecout+k) ztot(k)
      enddo
      do k=kz,nvrt
        kin=k-kz+1
        read(63,rec=irec+k) sigma(kin)
        write(65,rec=irecout+k) sigma(kin)
!        cs(kin)=(1-theta_b)*sinh(theta_f*sigma(kin))/sinh(theta_f)+ &
!     &theta_b*(tanh(theta_f*(sigma(kin)+0.5))-tanh(theta_f*0.5))/2/tanh(theta_f*0.5)
      enddo
      irec=irec+nvrt
      irecout=irecout+nvrt

!      print*, 'Done vgrid...'

!     Horizontal grid
      read(63,rec=irec+1) np
      read(63,rec=irec+2) ne
      irec=irec+2
!      if(np.gt.mnp.or.ne.gt.mne) then
!        write(*,*)'Too many nodes/elements',np,ne
!        stop
!      endif
      do m=1,np
!        read(63,rec=irec+1)x(m)
!        read(63,rec=irec+2)y(m)
!        read(63,rec=irec+3)dp(m)
        read(63,rec=irec+4)kbp00(m)
        irec=irec+4
      enddo !m=1,np
      do m=1,ne
!        read(63,rec=irec+1)i34
        irec=irec+1
        do mm=1,i34(m)
!          read(63,rec=irec+1)elnode(mm,m)
          irec=irec+1
        enddo !mm
      enddo !m

      write(65,rec=irecout+1) npout
      write(65,rec=irecout+2) neout
      irecout=irecout+2
      do m=1,npout
        write(65,rec=irecout+1)xout(m)
        write(65,rec=irecout+2)yout(m)
        write(65,rec=irecout+3)dpout(m)
        write(65,rec=irecout+4)kbp00(imapnode(m))
        irecout=irecout+4
      enddo !m=1,npout
      do m=1,neout
        write(65,rec=irecout+1)i34out(m)
        irecout=irecout+1
        do mm=1,i34out(m)
          write(65,rec=irecout+1)nmout(mm,m)
          irecout=irecout+1
        enddo !mm
      enddo !m

!      print*, 'last element: ',nmout(1:i34(neout),neout)

!...  Compute relative record # for a node and level for 3D outputs
!...
      if(iday==istart) then
        icount=0
        do i=1,np
          do k=max0(1,kbp00(i)),nvrt
            do m=1,ivs
              icount=icount+1
              icum(i,k,m)=icount
            enddo !m
          enddo !k
        enddo !i=1,np
      endif !iday==istart

      do it1=1,nrec
        read(63,rec=irec+1) time
        read(63,rec=irec+2) it
        irec=irec+2
        it_tot=it_tot+1
        time=it_tot*dtout
        write(65,rec=irecout+1) time
        write(65,rec=irecout+2) it_tot
        irecout=irecout+2

!        print*, 'time=',time/86400

        do i=1,np
          read(63,rec=irec+i) eta2(i)
        enddo !i
        irec=irec+np
        do i=1,npout
          write(65,rec=irecout+i) eta2(imapnode(i))
        enddo !i
        irecout=irecout+npout

        if(i23d==2) then
          do i=1,npout
            do m=1,ivs
              read(63,rec=irec+(imapnode(i)-1)*ivs+m) tmp
              write(65,rec=irecout+1) tmp
              irecout=irecout+1
            enddo !m
          enddo !i
          irec=irec+np*ivs
        else !i23d=3 
          do i=1,npout
            nd=imapnode(i)
            do k=max0(1,kbp00(nd)),nvrt
              do m=1,ivs
                read(63,rec=irec+icum(nd,k,m)) tmp
                write(65,rec=irecout+1) tmp
                irecout=irecout+1
              enddo !m
            enddo !k
          enddo !i
          irec=irec+icum(np,nvrt,ivs)
        endif !i23d
      enddo !it1=1,nrec

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      enddo !iday

      stop
      end

