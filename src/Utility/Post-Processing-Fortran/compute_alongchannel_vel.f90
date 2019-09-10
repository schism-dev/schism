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

!
!****************************************************************************************
!											*
!	Read hvel. and station.xy. (x,y) read in from (build points) station.xy, and z values
!       read in from the depths in station.xy (not used for 2D variables), and then compute
!       along-channel angles (max. variance) at each point and along-channel vel.
!       Works with mixed tri/quads.
!       Inputs: station.xy; screen inputs
!       Outputs: fort.20 (along channel vel.); fort.18 (channel angle)

!       ifort -Bstatic -assume byterecl -O3 -o compute_alongchannel_vel.exe ../UtilLib/compute_zcor.f90 ../UtilLib/pt_in_poly_test.f90 compute_alongchannel_vel.f90

!											*
!****************************************************************************************
!
      program read_out
      use compute_zcor
      use pt_in_poly_test
      parameter(nbyte=4)
      character*30 file63
      character*12 it_char
      character*48 start_time,version,variable_nm,variable_dim
      character*48 data_format
      integer,allocatable:: i34(:),elnode(:,:)
      allocatable :: sigma(:),sigma_lcl(:,:),cs(:),ztot(:),x(:),y(:),dp(:),kbp00(:),kbp(:),ztmp2(:,:)
      allocatable :: out(:,:,:,:),out2(:,:,:),icum(:,:,:),eta2(:),node3(:,:),arco(:,:)
      allocatable :: ztmp(:),x00(:),y00(:),iep(:),out3(:,:,:),z00(:),theta_al(:),out4(:,:)
      dimension swild(3)
      integer :: nodel(3)
      
      pi=3.1415926
!      print*, 'Input file to read from (without *_):'
!      read(*,'(a30)')file63
      file63='hvel.64'
      
      print*, 'Input start and end file # to read:'
      read(*,*) iday1,iday2

      print*, 'Input start CORIE day (at T=0):'
      read(*,*) istart_day

      open(10,file='station.xy',status='old')
      read(10,*) 
      read(10,*) nxy
      allocate(x00(nxy),y00(nxy),z00(nxy),theta_al(nxy),stat=istat)
      if(istat/=0) stop 'Falied to allocate (1)'
      do i=1,nxy
        read(10,*)j,x00(i),y00(i),z00(i) !z00 ==0 is MSL; z00<0 is below MSL
      enddo !i
      close(10)

!      open(65,file='extract.out')
!      write(65,*)'(x,y)= ',x00,y00
      
!...  Header
!...
      open(63,file='1_'//file63,status='old',access='direct',recl=nbyte)
      irec=0
      do m=1,48/nbyte
        read(63,rec=irec+m) data_format(nbyte*(m-1)+1:nbyte*m)
      enddo
      if(data_format.ne.'DataFormat v5.0') then
        print*, 'This code reads only v5.0:  ',data_format
        stop
      endif
      irec=irec+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) version(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) start_time(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) variable_nm(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) variable_dim(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte

      write(*,'(a48)')data_format
      write(*,'(a48)')version
      write(*,'(a48)')start_time
      write(*,'(a48)')variable_nm
      write(*,'(a48)')variable_dim

      read(63,rec=irec+1) nrec
      read(63,rec=irec+2) dtout
      read(63,rec=irec+3) nspool
      read(63,rec=irec+4) ivs
      read(63,rec=irec+5) i23d
      irec=irec+5

      print*, 'i23d=',i23d,' nrec= ',nrec

!     Vertical grid (obsolete)
      read(63,rec=irec+1) nvrt
!      read(63,rec=irec+2) kz
      read(63,rec=irec+3) h0
!      read(63,rec=irec+4) h_s
!      read(63,rec=irec+5) h_c
!      read(63,rec=irec+6) theta_b
!      read(63,rec=irec+7) theta_f
!      irec=irec+7
!
!      allocate(sigma(nvrt),cs(nvrt),ztot(nvrt),ztmp(nvrt),stat=istat)
!      if(istat/=0) stop 'Falied to allocate (2)'
!      do k=1,kz-1
!        read(63,rec=irec+k) ztot(k)
!      enddo
!      do k=kz,nvrt
!        kin=k-kz+1
!        read(63,rec=irec+k) sigma(kin)
!        cs(kin)=(1-theta_b)*sinh(theta_f*sigma(kin))/sinh(theta_f)+ &
!     &theta_b*(tanh(theta_f*(sigma(kin)+0.5))-tanh(theta_f*0.5))/2/tanh(theta_f*0.5)
!      enddo
!      irec=irec+nvrt

!     Horizontal grid
      irec=irec+7+nvrt
      read(63,rec=irec+1) np
      read(63,rec=irec+2) ne
      irec=irec+2

      allocate(x(np),y(np),dp(np),kbp00(np),i34(ne),elnode(4,ne), &
     &out(nxy,3,nvrt,2),out2(nxy,nvrt,2),icum(np,nvrt,2),eta2(np),node3(nxy,3), &
     &arco(nxy,3),iep(nxy),stat=istat)
      if(istat/=0) stop 'Falied to allocate (3)'

      do m=1,np
        read(63,rec=irec+1)x(m)
        read(63,rec=irec+2)y(m)
        read(63,rec=irec+3)dp(m)
        read(63,rec=irec+4)kbp00(m)
        irec=irec+4
      enddo !m=1,np
      do m=1,ne
        read(63,rec=irec+1)i34(m)
        irec=irec+1
        do mm=1,i34(m)
          read(63,rec=irec+1)elnode(mm,m)
          irec=irec+1
        enddo !mm
      enddo !m
      irec0=irec

      print*, 'last element',elnode(i34(ne),ne)

!     Read in vgrid.in
      allocate(ztot(nvrt),sigma(nvrt),sigma_lcl(nvrt,np),kbp(np))
      call get_vgrid_single('vgrid.in',np,nvrt,ivcor,kz,h_s,h_c,theta_b,theta_f,ztot,sigma,sigma_lcl,kbp)
      allocate(ztmp(nvrt),ztmp2(nvrt,3))

!...  Find parent element for (x00,y00)
      iep=0
      do i=1,ne
        do l=1,nxy
          if(iep(l)/=0) cycle

          call pt_in_poly_single(i34(i),x(elnode(1:i34(i),i)),y(elnode(1:i34(i),i)),x00(l),y00(l),inside,arco(l,1:3),nodel)
          if(inside==1) then
            iep(l)=i
            !print*, 'Found:',l,arco(l,1:3),nodel
            node3(l,1:3)=elnode(nodel(1:3),i)
          endif !inside
        enddo !l; build pts
!          aa=0
!          ar=0 !area
!          do j=1,3
!            j1=j+1
!            j2=j+2
!            if(j1>3) j1=j1-3
!            if(j2>3) j2=j2-3
!            n0=elnode(j,i)
!            n1=elnode(j1,i)
!            n2=elnode(j2,i)
!            swild(j)=signa(x(n1),x(n2),x00(l),y(n1),y(n2),y00(l)) !temporary storage
!            aa=aa+abs(swild(j))
!            if(j==1) ar=signa(x(n1),x(n2),x(n0),y(n1),y(n2),y(n0))
!          enddo !j
!          if(ar<=0) then
!            print*, 'Negative area:',ar
!            stop
!          endif
!          ae=abs(aa-ar)/ar
!          if(ae<=1.e-5) then
!            iep(l)=i
!            node3(l,1:3)=elnode(1:3,i)
!            arco(l,1:3)=swild(1:3)/ar
!            arco(l,1)=max(0.,min(1.,arco(l,1)))
!            arco(l,2)=max(0.,min(1.,arco(l,2)))
!            if(arco(l,1)+arco(l,2)>1) then 
!              arco(l,3)=0
!              arco(l,2)=1-arco(l,1)
!            else
!              arco(l,3)=1-arco(l,1)-arco(l,2)
!            endif
!            cycle
!          endif

        ifl=0 !flag
        do l=1,nxy
          if(iep(l)==0) then
            ifl=1
            exit
          endif
        enddo !l
        if(ifl==0) exit
      enddo !i=1,ne

      do j=1,nxy
        if(iep(j)==0) then
          print*, 'Cannot find a parent for pt:',j,x00(j),y00(j)
          stop
        endif
      enddo !j

!...  Compute relative record # for a node and level for 3D outputs
!...
      icount=0
      do i=1,np
        do k=max0(1,kbp00(i)),nvrt
          do m=1,ivs
            icount=icount+1
            icum(i,k,m)=icount
          enddo !m
        enddo !k
      enddo !i=1,np

!...  Time iteration
!...
      it_tot0=(iday1-1)*nrec
      it_tot=it_tot0
      nsteps=(iday2-iday1+1)*nrec
      allocate(out3(nsteps,nxy,2),out4(nsteps,2),stat=istat)
      if(istat/=0) stop 'Failed to allocate (5)'

      out3=0
      do iday=iday1,iday2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      write(it_char,'(i12)')iday
      open(63,file=it_char//'_'//file63,status='old',access='direct',recl=nbyte)

      irec=irec0
      do it1=1,nrec
!        read(63,rec=irec+1) time
        read(63,rec=irec+2) it
        irec=irec+2
        it_tot=it_tot+1
        !if(it_tot-it_tot0>mnit) stop 'Increase mnit'
        time=it_tot*dtout

!        print*, 'time=',time/86400

        do i=1,np
          read(63,rec=irec+i) eta2(i)
        enddo !i
        irec=irec+np

        out2=0
        if(i23d.eq.2) then
          do i=1,nxy
            do j=1,3 !nodes
              nd=node3(i,j)
              do m=1,ivs
                read(63,rec=irec+(nd-1)*ivs+m) tmp
                out2(i,1,m)=out2(i,1,m)+arco(i,j)*tmp
              enddo !m
            enddo !j
          enddo !i
          irec=irec+np*ivs
          write(18,'(e16.8,12(1x,f12.3))')time,(out2(i,1,1),i=1,nxy)
        else !i23d=3 
          do i=1,nxy
            do j=1,3 !nodes
              nd=node3(i,j)
              do k=max0(1,kbp00(nd)),nvrt
                do m=1,ivs
                  read(63,rec=irec+icum(nd,k,m)) out(i,j,k,m)
                enddo !m
              enddo !k
            enddo !j
          enddo !i
          irec=irec+icum(np,nvrt,ivs)

!         Do interpolation
          do i=1,nxy
            etal=0; dep=0; idry=0
            do j=1,3
              nd=node3(i,j)
              if(eta2(nd)+dp(nd)<h0) idry=1
              etal=etal+arco(i,j)*eta2(nd)
              dep=dep+arco(i,j)*dp(nd)
      
!             Debug
!              write(11,*)i,j,nd,dp(nd),arco(i,j)

            enddo !j
            if(idry==1) then
              if(file63(1:7).eq.'hvel.64') then
                out3(it_tot-it_tot0,i,1:2)=0
              else
                out3(it_tot-it_tot0,i,1:2)=-99
              endif
!              write(65,*)'Dry'
            else !element wet
!             Compute z-coordinates
              if(ivcor==1) then !localized
                do j=1,3
                  nd=node3(i,j)
                  do k=kbp(nd)+1,nvrt-1
                    ztmp2(k,j)=(eta2(nd)+dp(nd))*sigma_lcl(k,nd)+eta2(nd)
                  enddo !k
                  ztmp2(kbp(nd),j)=-dp(nd) !to avoid underflow
                  ztmp2(nvrt,j)=eta2(nd) !to avoid underflow
                enddo !j

                ztmp=0
                kbpl=minval(kbp(node3(i,1:3)))
                do k=kbpl,nvrt
                  do j=1,3
                    nd=node3(i,j)
                    ztmp(k)=ztmp(k)+arco(i,j)*ztmp2(max(k,kbp(nd)),j)
                  enddo !j
                enddo !k

              else if(ivcor==2) then !SZ
                call zcor_SZ_single(dep,etal,h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot, &
     &sigma,ztmp(:),idry2,kbpl)
              endif

              if(1==2) then
!-------------------------------------------------------------
              do k=kz,nvrt
                kin=k-kz+1
                hmod2=min(dep,h_s)
                if(hmod2<=h_c) then
                  ztmp(k)=sigma(kin)*(hmod2+etal)+etal
                else if(etal<=-h_c-(hmod2-h_c)*theta_f/sinh(theta_f)) then
                  write(*,*)'Pls choose a larger h_c (2):',etal,h_c
                  stop
                else
                  ztmp(k)=etal*(1+sigma(kin))+h_c*sigma(kin)+(hmod2-h_c)*cs(kin)
                endif

!               Following to prevent underflow
                if(k==kz) ztmp(k)=-hmod2
                if(k==nvrt) ztmp(k)=etal
              enddo !k

              if(dep<=h_s) then
                kbpl=kz
              else !z levels
!               Find bottom index
                kbpl=0
                do k=1,kz-1
                  if(-dep>=ztot(k).and.-dep<ztot(k+1)) then
                    kbpl=k
                    exit
                  endif
                enddo !k
                if(kbpl==0) then
                  write(*,*)'Cannot find a bottom level:',dep,i
                  stop
                endif
                ztmp(kbpl)=-dep
                do k=kbpl+1,kz-1
                  ztmp(k)=ztot(k)
                enddo !k
              endif

              do k=kbpl+1,nvrt
                if(ztmp(k)-ztmp(k-1)<=0) then
                  write(*,*)'Inverted z-level:',etal,dep,ztmp(k)-ztmp(k-1)
                  stop
                endif
              enddo !k
!-------------------------------------------------------------
              endif !1==2
       
              do k=kbpl,nvrt
                do m=1,ivs
                  do j=1,3
                    nd=node3(i,j)
                    kin=max(k,kbp00(nd))
                    out2(i,k,m)=out2(i,k,m)+arco(i,j)*out(i,j,kin,m)
                  enddo !j
                enddo !m
!                write(65,*)i,k,ztmp(k),(out2(i,k,m),m=1,ivs)     
              enddo !k

!             Interplate in vertical
              k0=0
              do k=kbpl,nvrt-1
                if(z00(i)>=ztmp(k).and.z00(i)<=ztmp(k+1)) then
                  k0=k
                  rat=(z00(i)-ztmp(k))/(ztmp(k+1)-ztmp(k))
                  exit
                endif
              enddo !k
              if(k0==0) then
                write(*,*)'Warning: failed to find a vertical level:',it,i,ztmp(kbpl),ztmp(nvrt)
                stop
!                Default for out3 is 0.
              else
                do m=1,ivs
                  out3(it_tot-it_tot0,i,m)=out2(i,k0,m)*(1-rat)+out2(i,k0+1,m)*rat !time index starts from 1
                enddo !m
              endif
            endif !dry/wet
          enddo !i=1,nxy
         
        endif !i23d
      enddo !it1=1,nrec

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      enddo !iday

      it_out=it_tot-it_tot0
      if(nsteps/=it_out) stop 'Mismatch in # of records'

!     Along-channel angles
      do i=1,nxy
        av_u=sum(out3(1:it_out,i,1))/it_out
        av_v=sum(out3(1:it_out,i,2))/it_out
        out4(1:it_out,1)=out3(1:it_out,i,1)-av_u
        out4(1:it_out,2)=out3(1:it_out,i,2)-av_v
        var_u=sum(out4(1:it_out,1)*out4(1:it_out,1))/it_out
        var_v=sum(out4(1:it_out,2)*out4(1:it_out,2))/it_out
        var_uv=sum(out4(1:it_out,1)*out4(1:it_out,2))/it_out
        theta_al(i)=0.5*atan2(2*var_uv,var_u-var_v)
        write(18,*)theta_al(i)/pi*180,z00(i)

!       Debug
!        write(99,*)i,z00(i),av_u,av_v,var_u,var_v,var_uv
      enddo !i
      print*, 'done computing angles...'

      do it=1,it_out
        time=(it+it_tot0)*dtout/86400+istart_day
        write(20,'(e16.8,6000(1x,f12.3))')time, &
     &(out3(it,i,1)*cos(theta_al(i))+out3(it,i,2)*sin(theta_al(i)),i=1,nxy)
      enddo !it

      print*, 'Finished!'

      stop
      end

!      function signa(x1,x2,x3,y1,y2,y3)
!!...  Compute signed area formed by pts 1,2,3
!
!      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2
!
!      return
!      end

