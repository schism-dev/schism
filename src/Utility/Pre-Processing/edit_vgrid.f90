!     Post-proc of vgrid.in after gen_vgs.f90, and do vertical interpolation if # of levels is changed
!     Inputs: nlev.gr3 (edit from output of gen_vqs.f90), vgrid.in
!     Output: vgrid.new 
!     pgf90 -O2 -mcmodel=medium  -Mbounds -Bstatic -o edit_vgrid edit_vgrid.f90
      allocatable :: sig0(:),sig(:),sigout(:)

      open(14,file='nlev.gr3',status='old')
      open(19,file='vgrid.in',status='old')
      open(20,file='vgrid.new',status='replace')
      read(14,*)
      read(14,*)ne,np
      read(19,*)
      read(19,*)nvrt
      write(20,*)1; write(20,*)nvrt
      allocate(sig0(nvrt),sig(nvrt),sigout(nvrt))
      do i=1,np
        read(14,*)j,x,y,tmp
        nlev=nint(tmp)
        read(19,*)j,kbp,sig(kbp:nvrt)

        nlev0=nvrt-kbp+1
        if(nlev0==nlev) then 
          kbp1=kbp
          sigout(kbp:nvrt)=sig(kbp:nvrt)
        else !do interpolation
          print*, 'Altering node ',i
          sig0(kbp:nvrt)=(/(-1+real(k-kbp)/(nvrt-kbp),k=kbp,nvrt)/) !even spaced base

          kbp1=nvrt-nlev+1
          sigout(kbp1)=-1; sigout(nvrt)=0
          do k=kbp1+1,nvrt-1
            sig1=-1+real(k-kbp1)/(nvrt-kbp1)
            k0=0 !flag
            do kk=kbp+1,nvrt
              if(sig1>=sig0(kk-1).and.sig1<=sig0(kk)) then
                k0=kk; exit
              endif
            enddo !kk
            if(k0==0) stop 'cannot find a level'

            rat=(sig1-sig0(k0-1))/(sig0(k0)-sig0(k0-1))
            sigout(k)=sig(k0)*rat+sig(k0-1)*(1-rat)
          enddo !k
     
          !Check
          do k=kbp1+1,nvrt
            if(sigout(k)<=sigout(k-1)) then
              print*, 'Inverted levels:',i,k
              stop
            endif
          enddo !k
        endif !nlev0

        !Output
        !write(20,'(2(1x,i10),100000(1x,f14.7))')i,kbp1,(sigout(max(k,kbp1)),k=1,nvrt)
        write(20,'(2(1x,i10),100000(1x,f14.7))')i,kbp1,sigout(kbp1:nvrt)
      enddo !i=1,np

      stop
      end
