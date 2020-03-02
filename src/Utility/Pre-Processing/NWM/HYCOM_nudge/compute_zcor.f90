!====================================================================
!  Routines to calculate various types of z-coord. in SCHISM
!  zcor_SZ
!  get_vgrid

!  ifort -Bstatic -O3 -c compute_zcor.f90
!  pgf90 -O2 -mcmodel=medium  -Bstatic -c compute_zcor.f90
!====================================================================
!====================================================================

      subroutine zcor_SZ(dp,eta,h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot,sigma,zcor,idry,kbp)
!     (Single-precision) Routine to compute z coordinates for SCHISM's SZ vertical system
!     Inputs:
!             dp: depth;
!             eta: elevation;
!             h0: min. depth;
!             h_s: transition depth between S and Z layers;
!             h_c: transition depth between S and sigma
!             theta_b, theta_f: spacing const. in S coordinate system;
!             nvrt: total # of vertical levels (S+Z);
!             kz: # of Z levels (1 if pure S);
!             ztot(1:kz):  z coordinates for Z levels; note that ztot(kz)=-h_s;
!             sigma(1:nvrt-kz+1): sigma coordinates for S (or sigma) levels; note that sigma(1)=-1, sigma(nvrt-kz+1)=0;
!     Outputs:
!             idry: wet (0) or dry (1) flag;
!             kbp: bottom index (0 if dry);
!             zcor(kbp:nvrt): z coordinates (junk if dry);    
!      implicit real*8(a-h,o-z)
      integer, intent(in) :: kz,nvrt
      real, intent(in) :: dp,eta,h0,h_s,h_c,theta_b,theta_f,ztot(nvrt),sigma(nvrt)
      integer, intent(out) :: idry,kbp
      real, intent(out) :: zcor(nvrt)

      !Local
      real :: s_con1,cs(nvrt)

!     Sanity check done before
!     Pre-compute constants
      s_con1=sinh(theta_f)
      nsig=nvrt-kz+1 !# of S levels 
      do k=1,nsig
        cs(k)=(1-theta_b)*sinh(theta_f*sigma(k))/sinh(theta_f)+ &
     &theta_b*(tanh(theta_f*(sigma(k)+0.5))-tanh(theta_f*0.5))/2/tanh(theta_f*0.5)
      enddo !k

      zcor=-99
      if(eta<=h0-h_s) then
        stop 'Deep depth dry'
      else if(eta+dp<=h0) then
        idry=1; kbp=0
      else !wet
!       S-levels
        idry=0
        hmod=min(h_s,dp)
        do k=kz,nvrt
          kin=k-kz+1
          if(hmod<=h_c) then
            zcor(k)=sigma(kin)*(hmod+eta)+eta
          else if(eta<=-h_c-(hmod-h_c)*theta_f/s_con1) then !hmod(i)>h_c>=0
            print*, 'Pls choose a larger h_c (2):',eta,h_c
            stop
          else
            zcor(k)=eta*(1+sigma(kin))+h_c*sigma(kin)+(hmod-h_c)*cs(kin)
          endif
        enddo !k=kz,nvrt

!       z-levels
        if(dp<=h_s) then
          kbp=kz
        else !bottom index 
          kbp=0 !flag
          do k=1,kz-1
            if(-dp>=ztot(k).and.-dp<ztot(k+1)) then
              kbp=k
              exit
            endif
          enddo !k
          if(kbp==0) then
            print*, 'Cannot find a bottom level for node (3):',i
            stop
          endif

          if(kbp>=kz.or.kbp<1) then
            print*, 'Impossible 92:',kbp,kz
            stop
          endif
          zcor(kbp)=-dp
          do k=kbp+1,kz-1
            zcor(k)=ztot(k)
          enddo !k
        endif

        do k=kbp+1,nvrt
          if(zcor(k)-zcor(k-1)<=0) then
            write(*,*)'Inverted z-levels at:',i,k,zcor(k)-zcor(k-1),eta,hmod
            stop
          endif
        enddo !k
      endif !wet ot dry

      end subroutine zcor_SZ

!====================================================================
      subroutine get_vgrid(vgrid,np,nvrt1,ivcor,kz,h_s,h_c,theta_b,theta_f,ztot,sigma,sigma_lcl,kbp)
!     (Single-precision) Routine to read in vgrid.in (either in the current dir or ../) 
!     Not all outputs have meaningful values, depending on ivcor.
!      implicit real*8(a-h,o-z)
      character(len=*), intent(in) :: vgrid
      integer, intent(in) :: np,nvrt1
      integer, intent(out) :: ivcor,kz,kbp(np)
      real, intent(out) :: h_s,h_c,theta_b,theta_f,ztot(nvrt1),sigma(nvrt1),sigma_lcl(nvrt1,np)

      !local
      logical :: lexist
      
      inquire(file=vgrid,exist=lexist)
      if(lexist) then
        open(19,file=vgrid,status='old')
      else
        inquire(file='../'//vgrid,exist=lexist)
        if(lexist) then
          open(19,file='../'//vgrid,status='old')
        else
          write(*,*)'Unable to find vgrid.in'
          stop
        endif
      endif

      read(19,*)ivcor
      select case(ivcor)
        case(2) !SZ
          read(19,*) nvrt,kz,h_s !kz>=1
          if(nvrt<2) stop 'nvrt<2'
          if(kz<1) then !.or.kz>nvrt-2) then
            write(*,*)'Wrong kz:',kz
            stop
          endif
          if(h_s<6) then
            write(*,*)'h_s needs to be larger:',h_s
            stop
          endif

          ! # of z-levels excluding "bottom" at h_s
          read(19,*) !for adding comment "Z levels"
          do k=1,kz-1
            read(19,*)j,ztot(k)
            if(ztot(k)>=-h_s) then
              print*, 'Illegal Z level:',k
              stop
            endif
            if(k>1) then; if(ztot(k)<=ztot(k-1)) then
              print*, 'z-level inverted:',k
              stop
            endif; endif
          enddo !k
          read(19,*) !level kz       
          ! In case kz=1, there is only 1 ztot(1)=-h_s
          ztot(kz)=-h_s

          nsig=nvrt-kz+1 !# of S levels (including "bottom" & f.s.)
          read(19,*) !for adding comment "S levels"
          read(19,*)h_c,theta_b,theta_f
          if(h_c<5.or.h_c>=h_s) then !large h_c to avoid 2nd type abnormality
            print*, 'h_c needs to be larger avoid 2nd type abnormality; &
     &do u want to continue? Enter 1 to continue, or ctrl-C to abort:'
            read*, itmp
          endif
          if(theta_b<0.or.theta_b>1) then
            write(*,*)'Wrong theta_b:',theta_b
            stop
          endif
          if(theta_f<=0) then
            write(*,*)'Wrong theta_f:',theta_f
            stop
          endif

          sigma(1)=-1 !bottom
          sigma(nsig)=0 !surface
          read(19,*) !level kz
          do k=kz+1,nvrt-1
            kin=k-kz+1
            read(19,*) j,sigma(kin)
            if(sigma(kin)<=sigma(kin-1).or.sigma(kin)>=0) then
              write(*,*)'Check sigma levels at:',k,sigma(kin),sigma(kin-1)
              stop
            endif
          enddo !k
          read(19,*) !level nvrt

        case(1) !localized sigma
          read(19,*)nvrt
          do i=1,np
            read(19,*)j,kbp(i),sigma_lcl(kbp(i):nvrt,i)
          enddo !i
        case default
          write(*,*)'Unknown ivcor:',ivcor
          stop
      end select
      close(19)      

      end subroutine get_vgrid

!====================================================================
