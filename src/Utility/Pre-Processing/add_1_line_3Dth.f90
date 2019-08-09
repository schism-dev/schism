!     Add 1 more line at the start of .th to convert to new format (from R1693).
!     Input:
!            (1) th.old (old *3D.th);
!     Output: th.new (the new *3D.th)

!     pgf90 -O2 -mcmodel=medium -Mbounds -Bstatic -o add_1_line_3Dth add_1_line_3Dth.f90
      program riverforcing
      parameter(nbyte=4)
      allocatable :: th(:),th1(:),th2(:),nond(:)

      print*, 'Input # of days:'
      read*, rndays
      print*, 'Input nvrt (1 for elev.61; 2*nvrt for uv3D.th; nvrt for T,S):'
      read*, nvrt
      print*, 'Input # of open bnd nodes (lump all):'
      read*, nond0
      allocate(th(nond0*nvrt),th1(nond0*nvrt),th2(nond0*nvrt),stat=istat)
      if(istat/=0) stop 'Failed too alloc.'

      irecl=nbyte*(1+nond0*nvrt)
      open(17,file='th.old',access='direct',recl=irecl,status='old')
      open(18,file='th.new',access='direct',recl=irecl,status='replace')
      read(17,rec=1)dt0,th2(:)
      write(18,rec=1)0.,th2(:)
      print*, 'Old & new time step (sec)=',dt0
      tt0=rndays*86400
      nt0=tt0/dt0+1.e-5

      do it=1,nt0
        read(17,rec=it)time2,th2(:)
        write(18,rec=it+1)time2,th2(:)
      enddo !it=1,nt0

      stop
      end 
