
!****************************************************************************************
!
      program create_dummy_zcor
      character*80  :: hgrid    !< base file name (e.g. elev.61) of files to be read 
      character*80  :: vgrid   !< vertical input file
      real          :: eta
      
      vgrid ="vgrid.in"
      inquire(file=vgrid,exist=lexist)
      if (.not.(lexist)) then
        vgrid = '../'//vgrid
        inquire(file=vgrid,exist=lexist)
        if(.not.(lexist)) then
          write(*,*)'Unable to find vgrid.in'
          stop
        endif
      endif
      
      hgrid ="hgrid.gr3"
      inquire(file=hgrid,exist=lexist)
      if (.not.(lexist)) then
        hgrid = '../'//hgrid
        inquire(file=hgrid,exist=lexist)
        if(.not.(lexist)) then
          write(*,*)'Unable to find hgrid.in'
          stop
        endif
      endif
      
      call dummy_zcor_input(eta)
      call  create_zcor(hgrid,vgrid,eta)
      end program
      
      subroutine dummy_zcor_input(eta)
      use argparse
      implicit none
      real, intent(out) :: eta   
      integer   :: comcount
      eta    = 0.0

      cmd_name = "create_dummy_zcor"
      call cla_init(cmd_name,"create a dummy 0_zcor.63 which contains output of zcore for one time step with specified uniform free surface elevation")

      call cla_register('-e','--eta', 'sythetic uniform free surface', cla_float,'0.0')
      
      call cla_validate

      comcount = command_argument_count()
      if (comcount > 0) then
          call cla_get("--eta",eta)
      end if
     
      return
      end subroutine

      subroutine create_zcor(hgrid,vgrid,eta)
      use compute_zcor
      character(len=*),intent(in) :: hgrid
      character(len=*),intent(in) :: vgrid
      real,intent(in)  :: eta
      parameter(nbyte=4, header_len=48)
      character*12 it_char
      character*header_len start_time,version,variable_nm,variable_dim
      character*header_len data_format
      integer       :: nrec,dtout,nspool,ivs,i23d_out,ilevel
      integer       :: nvrt,kz,ivcor,itemp,idry,kbp_p,kbpl
      real          :: h0,h_s,h_c,theta_b,theta_f
      integer,allocatable :: elnode(:,:),i34(:),kbp(:)
      real,allocatable :: sigma(:),cs(:),ztot(:),dp(:),kfp(:),sigma_lcl(:,:)
      real,allocatable :: out(:,:,:,:),eta2(:)
      real,allocatable :: zcor(:,:),x(:),y(:),z2(:,:)


!...  Header

      open(65,file="0_zcor.63",status='replace',form="unformatted",access="stream")
      
      
      data_format='DataFormat v5.0'
      variable_nm= 'zcor'!not important
      variable_dim= 'zcor'
      version = 'NA'
      start_time='NA'
      
   
      write(65) data_format
      write(65) version
      write(65) start_time
      write(65) variable_nm
      write(65) variable_dim
      
      nrec=1
      dtout=120 !not important
      nspool =720 ! ditto
      write(65) nrec
      write(65) dtout
      write(65) nspool
      i23d_out=3
      ivs =1
      write(65) ivs
      write(65) i23d_out

     

      open(22,file="coor.log",status='replace')      
      open(19,file=vgrid,status='old')
      read(19,*) ivcor
      read(19,*) nvrt
      close(19)
      print*, 'nvrt=',nvrt
      if(nvrt<2) then
          write(*,*) "total level less than 2"
          stop
      endif
      
      
      
      open(14,file=hgrid,status='old')
      read(14,*);
      read(14,*) ne,np
      allocate(x(np),y(np),dp(np))
      allocate(i34(ne),elnode(4,ne))
      do m=1,np
          read(14,*) itemp,x(m),y(m),dp(m)
      enddo
       do m=1,ne
          read(14,*) itemp,i34(m),(elnode(k,m),k=1,i34(m))
      enddo
      close(14)
      allocate(ztot(nvrt),sigma(nvrt),sigma_lcl(nvrt,np),kbp(np))
      
      ! give some default value
      kz = 1
      h0 = 0.0
      h_s = 4000
      h_c= 10
      theta_f = 0.001
      tehta_b = 0
      kbp=1
      if(nvrt>2) then
        call get_vgrid_single(trim(vgrid),np,nvrt,ivcor,kz,h_s,h_c,theta_b,theta_f,ztot,sigma,sigma_lcl,kbp)
      else 
        sigma(1)=-1.0
        sigma(2)=0.0
        sigma_lcl(1,:)=-1.0
        sigma_lcl(2,:)=0.0
      endif !nvrt
      
      if (ivcor.eq.1) then ! simga value is not important in localized sigma
         sigma(:) = 0.0
      endif
      
     !Vertical grid
      write(65) nvrt
      write(65) kz
      write(65) h0
      write(65) h_s
      write(65) h_c
      write(65) theta_b
      write(65) theta_f 
      write(65) (ztot(k),k=1,kz-1)
      write(65)  (sigma(kin),kin=1,nvrt-kz+1)
      
      !Horizontal grid
      write(65) np
      write(65) ne
      
      write(65) (x(m),y(m),dp(m),kbp(m),m=1,np)
      do m=1,ne
         write(65) i34(m), (elnode(mm,m),mm=1,i34(m))
      enddo !m
      
      
      write(65) 0 ! time
      write(65) 0 ! step
    
      
      write(65) (eta,i=1,np) ! surface set to 0 for every node
      
      allocate(zcor(nvrt,np))
      
      do m=1,np
       call compute_zcoor(m,ivcor,kz,nvrt,kbp(m),h0,dp(m),eta,h_s,h_c,theta_f,theta_b,ztot,sigma,sigma_lcl(:,m),zcor(:,m),kbpl)
       if (m<8) then
          write(22,*) m,dp(m),kbp(m), "elev:", (zcor(ilevel,m),ilevel=kbp(m),nvrt)
       end if
      end do
      
      write(65) ((zcor(k,i), k=max0(1,kbp(i)),nvrt),i=1,np)
      
      close(22)
      close(65)
      deallocate(ztot,sigma,sigma_lcl,kbp,zcor)
     
      stop
      end
      
      ! copy from misc_subs.F90 zcoor
      subroutine compute_zcoor(inode,ivcor,kz,nvrt,kbp,h0,dp,eta2,h_s,h_c,theta_f,theta_b,ztot,sigma,sigma_lcl,ztmp,kbpl)

      real,intent(in) :: dp,h0,h_s,h_c,theta_f
      real,intent(inout)::eta2
      integer,intent(in) :: ivcor,kz,nvrt,kbp
      real,intent(in)    :: ztot(kz),sigma(nvrt),sigma_lcl(nvrt)
      integer, intent(out) :: kbpl
      real, intent(out) :: ztmp(nvrt)

!     Local
      integer :: k,kin,m,klev
      real    :: hmod2,z0,z_1,sp,tmp,z_pws(nvrt),z_sigma(nvrt)
      real    :: s_con1
      real, allocatable :: cs(:)
      
      !Make sure it's wet
      if(dp+eta2<=h0) then
        !write(*,*)'ZCOOR: dry location:',inode,dp,eta2
        eta2=h0-dp+0.5
      endif
     
      s_con1=sinh(theta_f)
      allocate(cs(nvrt))
      if(ivcor==2) then !SZ
        hmod2=min(dp,h_s)
        ztmp(kz)=-hmod2 !to avoid underflow
        ztmp(nvrt)=eta2
        
        do klev=kz,nvrt
           k=klev-kz+1
           cs(k)=(1-theta_b)*sinh(theta_f*sigma(k))/sinh(theta_f)+ &
           &theta_b*(tanh(theta_f*(sigma(k)+0.5))-tanh(theta_f*0.5))/2/tanh(theta_f*0.5)
        enddo !klev
        
        do k=kz+1,nvrt-1
          kin=k-kz+1
          if(hmod2<=h_c) then
            ztmp(k)=sigma(kin)*(hmod2+eta2)+eta2
        
          else if(eta2<=-h_c-(dp-h_c)*theta_f/s_con1) then
            write(*,*)'ZCOOR: Pls choose a larger h_c:',eta2,h_c
          else
            ztmp(k)=eta2*(1+sigma(kin))+h_c*sigma(kin)+(hmod2-h_c)*cs(kin)
          endif
        enddo !k

        if(dp<=h_s) then
          kbpl=kz
        else !z levels
!         Find bottom index
          kbpl=0
          do k=1,kz-1
            if(-dp>=ztot(k).and.-dp<ztot(k+1)) then
              kbpl=k
              exit
            endif
          enddo !k
        
          if(kbpl==0) then
            write(*,*)'ZCOOR: Cannot find a bottom level:',dp,inode
          endif
          ztmp(kbpl)=-dp
          do k=kbpl+1,kz-1
            ztmp(k)=ztot(k)
          enddo !k
        endif !dep<=h_s

      else if(ivcor==1) then !localized simga
        !kbpl=kbp
        do k=kbp,nvrt
          ztmp(k)=(eta2+dp)*sigma_lcl(k)+eta2
        enddo !k
      else
        write(*,*) 'ZCOOR: unknown z-coor.'
        stop
      endif !ivcor

      do k=kbp+1,nvrt
        
        if(ztmp(k)-ztmp(k-1)<=0) then
          write(*,*)'ZCOOR: Inverted z-level:'
        endif
      enddo !k
      deallocate(cs)
      end subroutine compute_zcoor




