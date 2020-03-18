!This script is used to coupling the NWM streams. It returns the featureID
!of the intersected streams and the elementID of the intersected
!land boundary of the SCHISM grid.
!It assumes the land boundary is the 1st land boundary and there is only
!one land boundary.
!Inputs files: 
!      NWM_shp_ll.nc : which is converted from the shape file using the NWM streams. 
!      hgrid.ll: the gr3 format SCHISM grid; please only keep land
!      bnd segments needed for coupling (islands are OK). The # of land
!      bnd nodes is only used for dimensioning so give it a large number
!      if you want.
!             Note that the input files are in the same projection.
!             Here we use lat/lon projection
!      Epsilon: a small distance used to nudge the intersection point towards the downstream vetices
!               suggested value is 1.e-3, which is about 100m when using
!               lat/lon projection. 
!      dd,mm,yyyy: are day, month, and year, respectively provided by user.
! 
!Outputs files: 
!      source_sink.in   : contains the element ID for each
!                         intersection of the NWM streams and the land boundary.
!      msource.th       : contains the salinity and temprature
!                         of the source element along the land boundary.
!                         Salinity is set to be 0, temp = -9999.
!      vsource.th       : input of the stream flows of source elements.
!      vsind.th         : input of the stream flows of sink elements.
!      fort.95 fort.96  : output the maximum and average stream flow for
!      sourch and sink nodes.
!starting dates are set by users when asking for day, month, and year
!
!Author: Wei Huang
!
! 
!
!ifort -O2 -CB -g -traceback -o coupling_nwm  ~/git/schism/src/Utility/UtilLib/julian_date.f90 ~/git/schism/src/Utility/UtilLib/schism_geometry.f90 ~/git/schism/src/Utility/UtilLib/pt_in_poly_test.f90 coupling_nwm.f90 -I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -L$NETCDF/lib -lnetcdf -lnetcdff 

       program coupling_nwm
       use netcdf
       use schism_geometry_mod
       use pt_in_poly_test
       use ifport
       implicit real*8(a-h,o-z)
      
       character(len=*),parameter::REPODir='/sciclone/home20/whuang07/&
                                             data10/NWM/CHRTOUT/'
        
       character(len=*),parameter::FILE_NAME='NWM_shp_ll.nc'
       integer :: ncid
       character(len =*), parameter :: lat_NAME = 'lat'
       character(len =*), parameter :: lon_NAME = 'lon'
       character(len =*), parameter :: featureID_NAME = 'featureID'
       character(len =*), parameter :: origID_NAME = 'ORIG_FID'
       character(len=154)::NWMfile
       character(len=2)::dd1,mm1
       character(len=1)::hh1,hh2
       character(len=4)::yy
       integer :: lat_varid, lon_varid, featureID_varid, origID_varid
       integer :: ne, np, id, id1, ele
       integer :: nID,recorddimid,latdimid,nwmdimid,n_nwm
       integer :: nwmID_varid
       integer,parameter::max_len=100000
       integer,allocatable :: nwmID(:)         
       integer,allocatable ::featureID(:),origID(:),isbnd(:),nlnd(:),lndid(:,:)
       integer,allocatable::segn(:),seg_sink(:),seg_source(:),source_bnd(:),source_seg(:)
       integer,allocatable::bnd_sink(:),bnd_source(:),sink_bnd(:),sink_seg(:)
       integer,allocatable::eid(:),i34(:),elnode(:,:),ic3(:,:),elside(:,:),isdel(:,:),isidenode(:,:), nne(:),indel(:,:)
       integer,allocatable::uniso_ele(:),uniso_nwm(:),eleso_uni(:),nwmso_uni(:)
       integer,allocatable::unisi_ele(:),unisi_nwm(:),elesi_uni(:),nwmsi_uni(:)
       real*8,allocatable :: lats(:),lons(:),gx(:),gy(:),xcj(:),ycj(:)
       real*8,allocatable :: tstep(:)
       integer :: sf_varid,t_varid,n
       integer::ntime,nday,imm,idd,iyy,stime,timearray(9)
       integer,allocatable :: tm1(:),ty1(:),td1(:),td2(:),th1(:),th2(:)
       real*8::x11,x12,y11,y12,x21,x22,y21,y22,x3,y3
       real*8::arco(3),x4,y4,dup
       integer:: inside,nodel(3)
       real*8,allocatable :: salt(:),temp(:),streamflow(:) 
       real*8,allocatable :: SF_so(:,:),SF_si(:,:),sflow(:)
       integer,parameter::nodata=-999900 !nodata value in nwm stream flows

       !Read in epsilon (nudging ratio from side center to downstream pt
       !of NWM)
       print*, 'Input nudging ratio (suggest 1.e-3):'
       read*, epsilon

       print*, 'Input rain_rate:'
       rain_rate=0.03 !m/hour
       
       print*, 'Input number of days'
       read*,nday
       ntime=nday*24

       print*, 'Enter start time - dd,mm,yyyy (e.g. 1 1 1992)'
       read(*,*) idd,imm,iyy


       open(14,file='hgrid.ll',status='old')!lambert projection, NWM has shorter precision
       read(14,*)
       read(14,*)ne,np
       write(*,*)'# of elements',ne,'# of nodes',np
       allocate(gx(np),gy(np),isbnd(np),nne(np))
       allocate(eid(ne),i34(ne),elnode(4,ne))
       do i=1,np
         read(14,*)j,gx(i),gy(i),dp
       enddo     
       do i=1,ne
         read(14,*)eid(i),i34(i),elnode(1:i34(i),i)
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
       !Land bnds
       read(14,*) nland !total # of land segments
       read(14,*) nvel !total # of nodes for all land boundaries
       if(nland<=0.or.nvel<=0) stop 'nland<=0.or.nvel<=0'
       allocate(nlnd(nland),lndid(nvel,nland))

       do k=1,nland
         read(14,*)nlnd(k) !total # of nodes on the land boundary
         if(nlnd(k)>nvel) stop 'check # of land nodes'
         do i=1,nlnd(k)
           read(14,*)lndid(i,k)
         enddo !i
       enddo !k
       close(14)
       

       write(*,*)'# of nodes on land boundary=',nlnd(:),nland

!     We also need elem ball
       nne=0
       do i=1,ne
         do j=1,i34(i)
           nd=elnode(j,i)
           nne(nd)=nne(nd)+1
         enddo
       enddo
       mnei=maxval(nne)

       allocate(indel(mnei,np),stat=istat)
       if(istat/=0) stop 'Failed to alloc. indel'
       nne=0
       do i=1,ne
         do j=1,i34(i)
           nd=elnode(j,i)
           nne(nd)=nne(nd)+1
           if(nne(nd)>mnei) then
             write(*,*)'Too many neighbors',nd
             stop
           endif
           indel(nne(nd),nd)=i
         enddo
       enddo !i

       ! Compute geometry
       call compute_nside(np,ne,i34,elnode(1:4,1:ne),ns2)
       allocate(ic3(4,ne),elside(4,ne),isdel(2,ns2),isidenode(2,ns2),xcj(ns2),ycj(ns2),stat=istat)
       if(istat/=0) stop 'Allocation error: side(0)'
       call schism_geometry_double(np,ne,ns2,gx,gy,i34,elnode(1:4,1:ne),ic3(1:4,1:ne), &
     &elside(1:4,1:ne),isdel,isidenode,xcj,ycj)

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
       id=nID !size(lats) !# of vertices
       id1=maxval(origID)+1 !# of seg's
       write(*,*)'# NWM line segments=',id1
       !write(*,*)lons(1),lats(1)
       
       allocate(segn(id1)) !list of vertices along a seg
!       do i=1,id1+1
!         segn(i)=0
!         do j=1,id
!          if (origID(j).eq.i-1) then
!            segn(i)=segn(i)+1 ! get # of points on each line segments
!          
!          endif
!         enddo
!       enddo

       segn=0
       do i=1,nID
         itmp=origID(i)+1 !seg ID
         if(itmp>id1) then
           print*, 'Overflow:',itmp,id1,i
           stop
         endif
         segn(itmp)=segn(itmp)+1 
       enddo !i
      
       !p1 represents point pairs on NWM streams
       !p2 represents point pairs on grid land boudnary
       ninseg=0
       n_sink=0
       n_source=0
       allocate(seg_sink(max_len)) 
       allocate(bnd_sink(max_len))
       allocate(seg_source(max_len))
       allocate(bnd_source(max_len))

       do k1=1,nland !loop over all land segments of hgrid
       do k2=1,nlnd(k1)-1 !land seg      
         nd1=lndid(k2,k1); nd2=lndid(k2+1,k1)
         x21=gx(nd1);  y21=gy(nd1)
         x22=gx(nd2); y22=gy(nd2)

         !Bounding box
         xmin2=min(x21,x22); xmax2=max(x21,x22)
         ymin2=min(y21,y22); ymax2=max(y21,y22)

        !write(*,*)p2(1)%x,p2(2)%x,p2(1)%y,p2(2)%y
        !Find bnd elem
         ele=0 !init as flag
         isd0=0
         loop1: do m=1,nne(nd1) !ball
                  ie=indel(m,nd1)
                  do mm=1,i34(ie)
                    isd=elside(mm,ie)
                    if((isidenode(1,isd)==nd1.and.isidenode(2,isd)==nd2).or.(isidenode(1,isd)==nd2.and.isidenode(2,isd)==nd1)) then
!                     nd3=sum(elnode(1:i34(ie),ie))-nd1-nd2
                      isd0=isd
                      ele=ie
                      exit loop1
                    endif
                  enddo !mm
                enddo loop1 !m
         if(ele==0.or.isd0==0) then
           write(*,*)'failed to find an elem'
           stop
         endif
        
         write(99,*) 'Found bnd elem:',k1,k2,nd1,nd2,ele,isd0
      
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
         do j=1,id1 !NWM seg
         
           do i=1,segn(j)-1 !vertices
          
             x11=lons(n+i-1);   y11=lats(n+i-1);!upstream point
             x12=lons(n+i); y12=lats(n+i);!downstream point
             !Bounding box
             xmin1=min(x11,x12); xmax1=max(x11,x12)
             ymin1=min(y11,y12); ymax1=max(y11,y12)

             !Skip bounding boxes that r outside each other
             if(xmin1>xmax2.or.xmin2>xmax1.or.ymin1>ymax2.or.ymin2>ymax1) cycle
        
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
             x_inter=(b2-b1)/(a1-a2)
             y_inter=a1*x_inter+b1
!            write(*,*)x_inter,y_inter       
             if(a1.ne.a2.and.x_inter.ge.min(x21,x22).and.x_inter.le.max(x21,x22).and.&
                         y_inter.ge.min(y21,y22).and.y_inter.le.max(y21,y22).and.&
                         x_inter.ge.min(x11,x12).and.x_inter.le.max(x11,x12).and.&
                         y_inter.ge.min(y11,y12).and.y_inter.le.max(y11,y12)) then 
               ninseg=ninseg+1
       
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !--determine the source or sink by point in a poly-------!        
        !--select the fourth point which is calculated by adding
        ! a short distance from side center towards the second point of NWM stream
              x4=xcj(isd0)*(1-epsilon)+epsilon*x12
              y4=ycj(isd0)*(1-epsilon)+epsilon*y12

          !    write(*,*)nd1,nd2,ele
              call pt_in_poly_double(i34(ele),gx(elnode(1:i34(ele),ele)),gy(elnode(1:i34(ele),ele)),x4,y4,inside,arco,nodel)
          
              if(inside.eq.1) then
                n_source=n_source+1
                if(n_source.gt.max_len) then
                  stop 'n_source exceeds max_len'
                endif
                seg_source(n_source)=featureID(n)
                bnd_source(n_source)=ele
              else
                n_sink=n_sink+1
                if(n_sink.gt.max_len) then
                  stop 'n_sink exceeds max_len'
                endif
                seg_sink(n_sink)=featureID(n)
                bnd_sink(n_sink)=ele
              endif

              write(98,*)featureID(n),ele,ninseg
            endif !inside bounding box
        
          enddo !i: vertices
       !write(*,*)n
          n=n+segn(j)
        enddo !j !NWM seg

      enddo !k2: each land seg 
      enddo !k1: SCHISM land bnd seg
     
      write(98,*)'intersection#',ninseg,'sink#',n_sink,'source#',n_source

      allocate(sink_seg(max(1,n_sink))) 
      allocate(sink_bnd(max(1,n_sink)))

      do i=1,n_sink
        sink_seg(i)=seg_sink(i)
        sink_bnd(i)=bnd_sink(i)
      enddo

      allocate(source_seg(max(1,n_source)))
      allocate(source_bnd(max(1,n_source)))

      do i=1,n_source
        source_seg(i)=seg_source(i)
        source_bnd(i)=bnd_source(i)
      enddo
      
      deallocate(seg_sink,bnd_sink,seg_source,bnd_source)
!count unique number of sink and source elements

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
          endif
        enddo
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
        nsi=0
      endif
      if (n_source.ne.0) then
        allocate(uniso_ele(max_len))
        allocate(uniso_nwm(max_len))
        j=0;nso=0;
        do i=1,n_source-1
          if (source_bnd(i).eq.source_bnd(i+1)) then !remove duplicate
            j=j+1
          else
            nso=nso+1
            uniso_ele(nso)=source_bnd(i)
            uniso_nwm(nso)=source_seg(i)
          endif
        enddo
        allocate(eleso_uni(nso),nwmso_uni(nso))
        do i=1,nso
          eleso_uni(i)=uniso_ele(i) !element id is unique now
          nwmso_uni(i)=uniso_nwm(i) !nwm featureID is not necessarily unique
        enddo
      else
        nso=0
        allocate(uniso_ele(1))
        allocate(uniso_nwm(1))
        uniso_ele(1)=0
        uniso_nwm(1)=0
      endif 
      deallocate(uniso_ele,uniso_nwm)
      
!define time 
      allocate(tstep(ntime))
      tstep(1)=0
      do i=2,ntime
        tstep(i)=3600*(i-1)
      enddo
      !write(*,*)tstep(1),tstep(2)
!write msource.th file
      allocate(salt(nso),temp(nso))
      do i=1,nso
         salt(i)=0
         temp(i)=-9999
      enddo
      open(15,file='msource.th')
      do i=1,ntime
      write(15,'(50000F15.3)')tstep(i),(temp(k),k=1,nso),(salt(k),k=1,nso)
     
      enddo



 
!following is to read in the NWM output            
      allocate(tm1(ntime),td1(ntime),ty1(ntime),&
               th1(ntime),th2(ntime))
!convert time from julian date to calendar date:
! for function gmtime: Numeric time data to be formatted. Number of
! seconds since 00:00:00 Greenwich mean time, January 1, 1970
! The output of gmtime is One-dimensional array with 9 elements used to
! contain numeric time data.
! tarray(1)Seconds (0-61, where 60-61 can be returned for leap seconds)
! tarray(2)Minutes (0-59)
! tarray(3)Hours (0-23)
! tarray(4)Day of month (1-31)
! tarray(5)Month (0-11)
! tarray(6) Number of years since 1900
! tarray(7)Day of week (0-6, where 0 is Sunday)
! tarray(8)Day of year (0-365)
! tarray(9)Daylight saving flag (0 if standard time, 1 if daylight
! saving time)


      do i=1,ntime
         stime=julian_date(iyy,imm,idd)*24*3600
         stime=stime-julian_date(1970,1,1)*24*3600
         stime=stime+(i-1)*3600
         call gmtime(stime, timearray)
         tm1(i)=timearray(5)+1
         td1(i)=timearray(4)
         ty1(i)=timearray(6)+1900
      enddo
    
      do i=1,nday
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
      

!      allocate(SF_so(ntime,max(1,nso)))
!      allocate(SF_si(ntime,max(1,nsi)))
      if(n_source.ne.0) then
         allocate(SF_so(ntime,max(1,nso)))
      else
         allocate(SF_so(1,1))
      endif
      if(n_sink.ne.0) then
         allocate(SF_si(ntime,max(1,nsi)))
      else
         allocate(SF_si(1,1))
      endif
              
      do i=1,ntime
        write(yy,'(i4)')ty1(i)
        write(dd1,'(i2.2)')td1(i)
        write(mm1,'(i2.2)')tm1(i)        
        write(hh1,'(i1)')th1(i)
        write(hh2,'(i1)')th2(i)
         
       ! print*, yy,dd1,mm1,hh1,hh2

        NWMfile=REPODir//yy//mm1//dd1&
                     //hh1//hh2//'00.CHRTOUT_DOMAIN1.comp'
        write(98,*)'Inside time iteration:',i,NWMfile 
      
      
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
       ilow=minval(nwmID); ihigh=maxval(nwmID)
       allocate(sflow(ilow:ihigh))
       write(98,*)'ilow,ihigh:', ilow,ihigh

       do k=1,n_nwm
         itmp=nwmID(k)
         if(itmp<ilow.or.itmp>ihigh) then
           print*, 'Overflow(2):',itmp,ilow,ihigh
           stop
         endif
         if (streamflow(k)==nodata)then
           streamflow(k)=0
         endif 
         sflow(itmp)=0.01*streamflow(k)   
         !if(nwmID(k)==2590185)then
         !  write(96,*)streamflow(k)
         !endif  
       enddo !k
       !write(96,*)2590185,860,sflow(2590185)

!unique the element id during following process
!combine flux if there are muptiple streams going through one element

       if(n_source.ne.0) then
         dup=0;nso=0;
         do j=1,n_source-1
           if(source_seg(j)<ilow.or.source_seg(j)>ihigh) then
             write(97,*)'Overflow(3):',source_seg(j),ilow,ihigh
             SF_so(i,nso)=0
           !sum streams from same source elements
           else if (source_bnd(j).eq.source_bnd(j+1)) then
             dup=dup+sflow(source_seg(j))
           else
             nso=nso+1
             SF_so(i,nso)=dup+sflow(source_seg(j))
             dup=0
           endif
         enddo !j
       else
         SF_so(1,1)=0
       endif
       !write(*,*)'next step 1'
       if(n_sink.ne.0) then
         dup=0;nsi=0;
         do j=1,n_sink-1
           if(sink_seg(j)<ilow.or.sink_seg(j)>ihigh) then
             write(97,*)'Overflow(4):',sink_seg(j),ilow,ihigh
             SF_si(i,nsi)=0
           ! sum the flux from same sink elements
           else if (sink_bnd(j).eq.sink_bnd(j+1)) then
             dup=dup+sflow(sink_seg(j))
           else
             nsi=nsi+1
             SF_si(i,nsi)=dup+sflow(sink_seg(j))
             SF_si(i,nsi)=-SF_si(i,nsi)
             dup=0
           endif
         enddo !j
       else
         SF_si(1,1)=0
       endif
       !write(*,*)'next step 2'
       deallocate(streamflow,nwmID,sflow)
      enddo !i=1,ntime
      
      !print out statistics of streamflows
      if(nso.ne.0) then
        do i=1,nso
          write(96,*)maxval(SF_so(1:ntime,i),dim=1),sum(SF_so(1:ntime,i),dim=1)/ntime  
        enddo
      endif
      if(nsi.ne.0) then
        do i=1,nsi
          write(95,*)maxval(SF_si(1:ntime,i),dim=1),sum(SF_si(1:ntime,i),dim=1)/ntime
        enddo    
      endif

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




!write vsource.th file
      if (n_source.ne.0) then
        open(16,file='vsource.th',status='replace')
        do i=1,ntime
          write(16,'(10000F15.3)')tstep(i),(SF_so(i,k),k=1,nso)
        end do
      endif
!      write(*,*)'next step 4'
!write vsink.th
      if (n_sink.ne.0) then
        open(17,file='vsink.th')
        do i=1,ntime
          write(17,'(10000F15.3)')tstep(i),(SF_si(i,k),k=1,nsi)
        end do
      endif
!      write(*,*)'next step 5'

     
      

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


 
