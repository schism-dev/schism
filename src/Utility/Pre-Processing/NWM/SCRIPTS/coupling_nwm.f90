!This script is used to coupling the NWM streams. It returns the featureID
!of the intersected streams and the elementID of the intersected
!land boundary of the SCHISM grid.
!It assumes the land boundary is the 1st land boundary and there is only
!one land boundary.
!Inputs files: 
!      *.nc : which is converted from the shape file using the NWM streams. 
!      hgrid: which is the gr3 format SCHISM grid.
!             Note that the input files are in the same projection.
!Outputs files: 
!      source_sink.in   : contains the element ID for each
!                         intersection of the NWM streams and the land boundary.
!      msource.th       : contains the salinity and temprature
!                         of the source element along the land boundary.
!                         Salinity is set to be 0, temp = -9999.
!      vsource.th       : input of the stream flows of source elements.
!      vsind.th         : input of the stream flows of sink elements.
!
!
!ifort coupling_nwm.f90 -I$NETCDF/include -I$NETCDF_FORTRAN/include
!-L$NETCDF_FORTRAN/lib -L$NETCDF/lib -L$NETCDF/lib -lnetcdf -lnetcdff
!~/git/schism/src/Utility/UtilLib/pt_in_poly_test.f90

      program coupling_nwm
       use netcdf
       use pt_in_poly_test
       implicit real*8(a-h,o-z)
      
       character(len=*),parameter::REPODir='/sciclone/data10/feiye/&
       &DELAWARE_REPO/00_TWL_Shared/01_data/01-NWM-4-isabel-irene-sandy-&
       &13sep2018/NetCDF_HDF/Irene_output/'  
        
       character(len=*),parameter::FILE_NAME='nwm_channels_0-30m_LCC.nc'!'nwm_channels_0-30m_LCC.nc'
       integer :: ncid
       character(len =*), parameter :: lat_NAME = 'y'
       character(len =*), parameter :: lon_NAME = 'x'
       character(len =*), parameter :: featureID_NAME = 'featureID'
       character(len =*), parameter :: origID_NAME = 'ORIG_FID'
       character(len=154)::NWMfile
       character(len=1)::dd1,dd2,mm1,mm2,hh1,hh2
       character(len=4)::yy
       integer :: lat_varid, lon_varid, featureID_varid, origID_varid
       integer :: ne, np, id, id1, ele
       integer :: nID,recorddimid,latdimid,nwmdimid,n_nwm
       integer :: nwmID_varid
       integer,parameter::max_len=100000
       integer,allocatable :: nwmID(:)         
       integer,allocatable ::featureID(:),origID(:),isbnd(:),lndid(:)
       integer,allocatable::segn(:),seg_sink(:),seg_source(:),source_bnd(:),source_seg(:)
       integer,allocatable::bnd_sink(:),bnd_source(:),sink_bnd(:),sink_seg(:)
       integer,allocatable::eid(:),n1(:),n2(:),n3(:)
       integer,allocatable::uniso_ele(:),uniso_nwm(:),eleso_uni(:),nwmso_uni(:)
       integer,allocatable::unisi_ele(:),unisi_nwm(:),elesi_uni(:),nwmsi_uni(:)
       real*8,allocatable :: lats(:),lons(:),gx(:),gy(:)
       real*8,allocatable :: tstep(:)
       integer :: sf_varid,t_varid,tm1
       integer,allocatable :: tm2(:),td1(:),td2(:),th1(:),th2(:)
       real*8::x11,x12,y11,y12,x21,x22,y21,y22,x3,y3
       real*8::xp(3),yp(3),arco(3)
       integer:: inside,nodel(3)
       real*8,allocatable :: salt(:),temp(:),streamflow(:) 
       real*8,allocatable :: SF_so(:,:),SF_si(:,:)

       open(14,file='hgrid.lcc',status='old')!lambert projection, NWM has shorter precision
       read(14,*)
       read(14,*)ne,np
       write(*,*)'# of elements',ne,'# of nodes',np
       allocate(gx(np),gy(np),isbnd(np))
       allocate(eid(ne),n1(ne),n2(ne),n3(ne))       
        do i=1,np
        read(14,*)j,gx(i),gy(i),dp
        enddo     
        do i=1,ne
        read(14,*)eid(i),j,n1(i),n2(i),n3(i)
        enddo
        !     Bnd
        isbnd=0
        read(14,*) nope
        read(14,*) neta
        ntot=0
        do k=1,nope
         read(14,*) nond
         do i=1,nond
          read(14,*) nd !iond(k,i)
          isbnd(nd)=k
         enddo !i
        enddo !k
       !write(*,*)gx(1),gy(1) 
       !     Land bnds
        read(14,*) nland !total # of land segments
        read(14,*) nvel !total # of nodes for all land boundaries
        !do k=1,nland
         read(14,*)nlnd !total # of nodes on the first land boundary
         allocate(lndid(nlnd))
         do i=1,nlnd
          read(14,*)lndid(i)
         enddo !i
        !enddo !k
       close(14)
       

       write(*,*)'# of the first land boundary',nlnd
       write(*,*)'islands are excluded'
       
       call check(nf90_open(FILE_NAME, nf90_nowrite,ncid))      
       !get the size of the featureID
       call check(nf90_inq_dimid(ncid,'RecordID', latdimid))   
       call check(nf90_inquire_dimension(ncid, latdimid, len = nID))
       write(*,*)'# of NWM points',nID 
       allocate(featureID(nID))     
       allocate(origID(nID))
       allocate(lats(nID))
       allocate(lons(nID))
       call check(nf90_inq_varid(ncid,lat_NAME,lat_varid))
       call check(nf90_inq_varid(ncid,lon_NAME,lon_varid))
       call check(nf90_inq_varid(ncid,featureID_NAME,featureID_varid))
       call check(nf90_inq_varid(ncid,origID_NAME,origID_varid))

       call check(nf90_inq_varid(ncid,featureID_NAME,featureID_varid))
       call check(nf90_get_var(ncid,lat_varid,lats))
       call check(nf90_get_var(ncid,lon_varid,lons))
       call check(nf90_get_var(ncid,origID_varid,origID))
       call check(nf90_get_var(ncid,featureID_varid,featureID))
       call check(nf90_close(ncid))
       id=size(lats)
       id1=maxval(origID)
       write(*,*)'# line segments=',id1+1
       !write(*,*)lons(1),lats(1)
       
       allocate(segn(id1+1)) !total # of line segments
       do i=1,id1+1
         segn(i)=0
         do j=1,id
          if (origID(j).eq.i-1) then
          segn(i)=segn(i)+1 ! get # of points on each line segments
          
          endif
         enddo
       enddo
      
       !p1 represents point pairs on NWM streams
       !p2 represents point pairs on grid land boudnary
       ninseg=0
       n_sink=0
       n_source=0
       allocate(seg_sink(max_len)) 
       allocate(bnd_sink(max_len))
       allocate(seg_source(max_len))
       allocate(bnd_source(max_len))
      do k=1,nlnd-1 !hgrid land seg      
        x21=gx(lndid(k));   y21=gy(lndid(k));
        x22=gx(lndid(k+1)); y22=gy(lndid(k+1));
        !write(*,*)p2(1)%x,p2(2)%x,p2(1)%y,p2(2)%y
        !Find bnd elem
        do m=1,ne
         if(n1(m).eq.lndid(k).and.n2(m).eq.lndid(k+1)) then
          x3=gx(n3(m));y3=gy(n3(m));ele=m;
         elseif(n1(m).eq.lndid(k+1).and.n2(m).eq.lndid(k)) then
          x3=gx(n3(m));y3=gy(n3(m));ele=m; 
         elseif(n3(m).eq.lndid(k).and.n2(m).eq.lndid(k+1)) then
          x3=gx(n1(m));y3=gy(n1(m));ele=m;
         elseif(n3(m).eq.lndid(k+1).and.n2(m).eq.lndid(k)) then
          x3=gx(n1(m));y3=gy(n1(m));ele=m;
         elseif(n3(m).eq.lndid(k).and.n1(m).eq.lndid(k+1)) then
          x3=gx(n2(m));y3=gy(n2(m));ele=m;
         elseif(n3(m).eq.lndid(k+1).and.n1(m).eq.lndid(k)) then
          x3=gx(n2(m));y3=gy(n2(m));ele=m;
         else
          ! stop 'no element found'
         endif
        enddo !m
              
        dx2 = x21-x22
        dy2 = y21-y22        
        if(dx2==0.) then
           a2=0.
           b2=y21
         else
           a2=dy2/dx2
           b2=y21-a2*x21
         endif
       !write(*,*)'a2,b2',a2,b2
       n=1
       do j=1,id1+1 !NWM seg
         
        do i=1,segn(j)-1 !vertices
          
         x11=lons(n+i-1);   y11=lats(n+i-1);!upstream point
         x12=lons(n+i); y12=lats(n+i);!downstream point
        
         dx1 = x11-x12
         dy1 = y11-y12
         !write(*,*)p1(1)%x,p1(1)%y,dx1,dy1
         if(dx1==0.) then
           a1=0.
           b1=y11
         else
           a1=dy1/dx1 !dacter(len=*),parameterx1
           b1=y11-a1*x11
         endif
         x=(b2-b1)/(a1-a2)
         y=a1*x+b1
         !write(*,*)x,y,p1(1)%x,p2(1)%x,p1(1)%y,p2(1)%y       
         if(a1.ne.a2.and.x.ge.min(x21,x22).and.x.le.max(x21,x22).and.&
                         y.ge.min(y21,y22).and.y.le.max(y21,y22).and.&
                         x.ge.min(x11,x12).and.x.le.max(x11,x12).and.&
                         y.ge.min(y11,y12).and.y.le.max(y11,y12)) then 
            ninseg=ninseg+1
        !-----determine the source or sink by triangle area-------!
        !calculating area of the triangle in the land boundary
        !
        !    area=(x21*(y22-y3)+x22*(y3-y21)+x3*(y21-y22))/2
        !    A1= (x12*(y22-y3)+x22*(y3-y12)+x3*(y12-y22))/2
        !    A2= (x21*(y12-y3)+x12*(y3-y21)+x3*(y21-y12))/2
        !    A3= (x21*(y22-y12)+x22*(y12-y21)+x12*(y21-y22))/2
        !    Atotal=A1+A2+A3
        !    if(Atotal.eq.area) then
        !      n_source=n_source+1
        !      seg_source(n_source)=featureID(n)
        !      bnd_source(n_source)=ele
        !    else
        !      n_sink=n_sink+1
        !      seg_sink(n_sink)=featureID(n)
        !      bnd_sink(n_sink)=ele
        !    endif 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !--determine the source or sink by point in a poly-------!        
        !--select the fourth point which is calculated by adding
        ! a short distance towards the second point of NWM stream
            x4=x+(x12-x)*1.e-8; y4=y+(y12-y)*1.e-8;
            xp(1)=x21;xp(2)=x22;xp(3)=x3;
            yp(1)=y21;yp(2)=y22;yp(3)=y3;
            !write(*,*)lndid(k),lndid(k+1),ele
            call pt_in_poly_double(3,xp,yp,x12,y12,inside,arco,nodel)
            if(inside.eq.1) then
               n_source=n_source+1
               seg_source(n_source)=featureID(n)
               bnd_source(n_source)=ele
            else
               n_sink=n_sink+1
               seg_sink(n_sink)=featureID(n)
               bnd_sink(n_sink)=ele
            endif

            write(*,*)featureID(n),ele,ninseg
         else
         
         endif
        
        enddo !i
       !write(*,*)n
       n=n+segn(j)
       enddo !j
        if(n_source.gt.max_len.or.n_sink.gt.max_len)then
        stop 'n_source or n_sink is exceeds max_len'
        endif
      enddo !k
     
      
      !write(*,*)'k=',k,'j',j,'i',i
      write(*,*)'intersection#',ninseg,'sink#',n_sink,'source#',n_source
      allocate(sink_seg(max(1,n_sink))) 
      allocate(sink_bnd(max(1,n_sink)))
!      if (n_sink.ne.0) then
        do i=1,n_sink
         sink_seg(i)=seg_sink(i)
         sink_bnd(i)=bnd_sink(i)
        enddo
!      else
!         allocate(sink_seg(0))
!         allocate(sink_bnd(0))
!         write(*,*)'no sink element'
!      endif 
      allocate(source_seg(max(1,n_source)))
      allocate(source_bnd(max(1,n_source)))
!      if  (n_source.ne.0) then
      do i=1,n_source
        source_seg(i)=seg_source(i)
        source_bnd(i)=bnd_source(i)
      enddo
!      else
!         allocate(source_seg(0))
!         allocate(source_bnd(0))
!         write(*,*)'no source element'
!      endif
      
      deallocate(seg_sink,bnd_sink,seg_source,bnd_source)
!unique the element id
           
      if (n_sink.ne.0) then
       allocate(unisi_ele(max_len))
       allocate(unisi_nwm(max_len))
       j=0;nsi=0;
       do i=1,n_sink-1
          
         if (sink_bnd(i).eq.sink_bnd(i+1)) then !remove duplicate
           j=j+1
         else
           nsi=nsi+1
           unisi_ele(nsi)=sink_bnd(i) !unique hgrid elem ID
           unisi_nwm(nsi)=sink_seg(i) !NWM feature ID
          ! write(*,*)unisi_ele(nsi),unisi_nwm(nsi),nsi,j
         endif
       enddo
!       if (sink_bnd(n_sink-1).eq.sink_bnd(n_sink)) then
!          nsi=nsi-1
!       else
!          unisi_ele(nsi)=sink_bnd(n_sink)
!          unisi_nwm(nsi)=sink_seg(n_sink)
!       endif 
       allocate(elesi_uni(nsi),nwmsi_uni(nsi))
       do i=1,nsi
          elesi_uni(i)=unisi_ele(i) !element id is unique now
          nwmsi_uni(i)=unisi_nwm(i) !nwm featureID is not necessarily unique
       enddo
      else
       allocate(unisi_ele(1))
       allocate(unisi_nwm(1))
       unisi_ele(1)=0
       unisi_nwm(1)=0
      endif
     ! write(*,*)'next step'
      deallocate(unisi_ele,unisi_nwm)
     ! write(*,*)'next step'
      if (n_source.ne.0) then
       allocate(uniso_ele(max_len))
       allocate(uniso_nwm(max_len))
       j=0;nso=0;
       do i=1,n_source-1
 
         if (source_bnd(i).eq.source_bnd(i+1)) then
           j=j+1
         else
           nso=nso+1
           uniso_ele(nso)=source_bnd(i)
           uniso_nwm(nso)=source_seg(i)
       !   write(*,*)uniso_ele(nso),uniso_nwm(nso),nso,j
         endif
       enddo
!       if (source_bnd(n_source-1).eq.source_bnd(n_source)) then
!          nso=nso-1
!       else
!          uniso_ele(nso)=source_bnd(n_source)
!          uniso_nwm(nso)=source_seg(n_source)
!       endif
       allocate(eleso_uni(nso),nwmso_uni(nso))
       do i=1,nso
          eleso_uni(i)=uniso_ele(i) !element id is unique now
          nwmso_uni(i)=uniso_nwm(i) !nwm featureID is not necessarily unique
       enddo
      else
       allocate(uniso_ele(1))
       allocate(uniso_nwm(1))
       uniso_ele(1)=0
       uniso_nwm(1)=0
      endif
      deallocate(uniso_ele,uniso_nwm)
      
!write source_sink.in file
      open(13,file='source_sink.in') 
      if (n_source.ne.0) then
        write(13,*)nso
        do i=1,nso
           write(13,*)eleso_uni(i)
        enddo
      else
        write(13,*)0
      endif
      write(13,*)
      if (n_sink.ne.0) then
        write(13,*)nsi
        do i=1,nsi
           write(13,*)elesi_uni(i)
        enddo
      else
        write(13,*)0
      endif
!define time 
      allocate(tstep(80*24))
      tstep(1)=0
      do i=1,80*24
      tstep(i+1)=3600*i
      enddo
      !write(*,*)tstep(1),tstep(2)
!write msource.th file
      allocate(salt(nso),temp(nso))
      do i=1,nso
         salt(i)=0
         temp(i)=-9999
      enddo
      open(15,file='msource.th')
      do i=1,80*24
      write(15,'(10000F15.3)')tstep(i),(temp(k),k=1,nso),(salt(k),k=1,nso)
     
      enddo



 
!following is to read in the NWM output            
      allocate(tm2(80*24),&
               td1(80*24),td2(80*24),&
               th1(80*24),th2(80*24))
      tm1=0
      do i=1,80*24
         if(i<=31*24) then
          tm2(i)=7
         elseif (i>31*24.and.i<=62*24) then 
          tm2(i)=8
         else 
          tm2(i)=9
         endif
      enddo
      
       do j=1,31*24
            if(int((j-1)/24)<9) then
            td1(j)=0
            td2(j)=1+int((j-1)/24)
            else
            td=int((j-1)/24)+1
            td1(j)=int(td/10)
            td2(j)=td-10*td1(j)
            endif
       enddo
       k=31*24
       do j=k+1,k+31*24
            if(int((j-k-1)/24)<9) then
            td1(j)=0
            td2(j)=1+int((j-k-1)/24)
            else
            td=int((j-k-1)/24)+1
            td1(j)=int(td/10)
            td2(j)=td-10*td1(j)
            endif
       enddo
       k=31*24+31*24
       do j=k+1,k+18*24
            if(int((j-k-1)/24)<9) then
            td1(j)=0
            td2(j)=1+int((j-k-1)/24)
            else
            td=int((j-k-1)/24)+1
            td1(j)=int(td/10)
            td2(j)=td-10*td1(j)
            endif
       enddo
      do i=1,80
         do j=24*(i-1)+1,24*i
           if((j-24*(i-1))<=10) then
             th1(j)=0
             th2(j)=int(j-24*(i-1))-1
           else
             th=int((j-24*(i-1)))-1
             th1(j)=int(th/10)
             th2(j)=th-10*th1(j)
           endif
         enddo
       enddo
      

      allocate(SF_so(80*24,max(1,nso)))
      allocate(SF_si(80*24,max(1,nsi)))
!      if(n_source.ne.0) then
!      else
!         allocate(SF_so(1,1))
!      endif
!      if(n_sink.ne.0) then
!      else
!         allocate(SF_si(1,1))
!      endif
              
      do i=1,80*24
        write(yy,'(i4)')2011
        write(dd1,'(i1)')td1(i)
        write(dd2,'(i1)')td2(i)
        write(mm1,'(i1)')tm1
        write(mm2,'(i1)')tm2(i)
        write(hh1,'(i1)')th1(i)
        write(hh2,'(i1)')th2(i)

        NWMfile=REPODir//yy//mm1//mm2//dd1//dd2&
                     //hh1//hh2//'00.CHRTOUT_DOMAIN1.comp'
        write(*,*)NWMfile 
      
      
       call check(nf90_open(NWMfile, nf90_nowrite,ncid))
       call check(nf90_inq_dimid(ncid,'feature_id', nwmdimid))
       call check(nf90_inquire_dimension(ncid, nwmdimid, len = n_nwm))
       !write(*,*)'# of NWM output streams',n_nwm
       allocate(streamflow(n_nwm))
       allocate(nwmID(n_nwm))
       call check(nf90_inq_varid(ncid,'streamflow',sf_varid))
       call check(nf90_inq_varid(ncid,'time',t_varid))
       call check(nf90_inq_varid(ncid,'feature_id',nwmID_varid))
       call check(nf90_get_var(ncid,sf_varid,streamflow))
       call check(nf90_get_var(ncid,t_varid,time1))
       call check(nf90_get_var(ncid,nwmID_varid,nwmID))
       !write(*,*)nso,nsi,nwmso_uni(1)
       if(nso.ne.0) then
        do j=1,nso
         do k=1,n_nwm
         if(nwmID(k).eq.nwmso_uni(j)) then
            SF_so(i,j)=streamflow(k)
         endif
         enddo
        enddo
       else
         SF_so(1,1)=0
       endif
      !write(*,*)'next step'
       if(nsi.ne.0) then
        do j=1,nsi
         do k=1,n_nwm
         if(nwmID(k).eq.nwmsi_uni(j)) then
            SF_si(i,j)=-streamflow(k)
         endif
         enddo
        enddo
       else
         SF_si(1,1)=0
       endif
      ! write(*,*)'next step',i
       deallocate(streamflow,nwmID)
      enddo
      !write(*,*)'next step'      
      !write(*,*)size(SF_so),SF_so(1,1),SF_so(10,10)
!write vsource.th file
      open(16,file='vsource.th')
      if (n_source.ne.0) then
       do i=1,80*24
       write(16,'(10000F15.3)')tstep(i),(SF_so(i,k),k=1,nso)
       end do
      else
       write(16,*)0
      endif
      
!write vsink.th
      open(17,file='vsink.th')
      if (n_sink.ne.0) then
       do i=1,80*24
       write(17,'(10000F15.8)')tstep(i),(SF_si(i,k),k=1,nsi)
       end do
      else
       write(17,*)0
      endif


     
      

      call check(nf90_close(ncid))


      contains
         subroutine check(status)
          integer, intent ( in) :: status
    
          if(status /= nf90_noerr) then 
          print *, trim(nf90_strerror(status))
          stop "Stopped"
          end if
         end subroutine check  
      end


 
