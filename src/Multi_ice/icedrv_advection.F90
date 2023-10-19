!=======================================================================
!
! This submodule contains the subroutines 
! for the advection of sea ice tracers
!
! Author: Lorenzo Zampieri ( lorenzo.zampieri@awi.de )
!        Qian Wang upate ICEPACK to 1.3.4
!  Modified by Qian Wang to apply to SCHISM
! Main driver: tracer_advection_icepack, which 
! calls fct_solve_icepack (sequence of lower- and higher-order FCT solver for
! each tracer)
!
! 2023.10   including upwind, central difference and TVD schemes by Qian Wang
!=======================================================================

submodule (icedrv_main) icedrv_advection

    use icedrv_kinds
    use icedrv_constants
    use icedrv_system,      only: icedrv_system_abort
    use schism_msgp,        only: parallel_abort,parallel_finalize,exchange_p2d,rtype,comm
    use icepack_intfc,      only: icepack_warnings_flush,         &
                                  icepack_warnings_aborted,       &
                                  icepack_query_tracer_indices,   &
                                  icepack_query_tracer_flags,     &
                                  icepack_query_parameters,       &
                                  icepack_query_tracer_sizes
    use schism_glbl, only: rkind,ne,nea,np,npa,elnode,nnp,indnd,time_stamp,rnday, &
   &fdb,lfdb,xnd,ynd,znd,iplg,ielg,idry,i34,isbnd,xlon,ylat,pframe,eframe,area,idry_e,nxq, &
   &elside,iself,indel,nne,distj,snx,sny,isbs,xel,yel,xctr,yctr,xcj,ycj,ics,isidenode, &
   &sframe2,mnei_p,ipgl

    implicit none

    real(kind=dbl_kind), allocatable, dimension(:)   :: &
         d_tr,        trl,                              &
         rhs_tr,      rhs_trdiv,                        &
         icepplus,    icepminus,                        &
         mass_matrix, newwork, aream,flux_ratio_casulli                

    real(kind=dbl_kind), allocatable, dimension(:,:) :: &
         icefluxes

    real(kind=dbl_kind), allocatable, dimension(:,:,:) :: &
         flux_matrix1,flux_matrix2,flux_matrix3,flux_matrix4,flux_matrix00,flux_matrix_casulli

    real(kind=dbl_kind), allocatable, dimension(:,:,:,:) :: &
         flux_matrix0
    ! Variables needed for advection

    contains

    subroutine tg_rhs_icepack(trc)
    
        !use mod_mesh
        use mice_module
        !use g_parsup
        !use o_param
        !use g_config
    
        implicit none
    
        ! Input - output
    
        !type(t_mesh),        target,        intent(in)     :: mesh
        real(kind=dbl_kind), dimension(nx), intent(inout)  :: trc  !nx=npa
    
        ! Local variables
    
        real(kind=dbl_kind)              :: diff,       um,    vm,    vol, & 
                                            entries(3), dx(3), dy(3)
        integer(kind=int_kind)           :: n,          q,     row,        &
                                            elem,       elnodes(3)
    
!#include "../associate_mesh.h"
         !scalevol=2.0e8

        ! Taylor-Galerkin (Lax-Wendroff) rhs
      
        do row = 1, nx_nh !=np (no halo)
           rhs_tr(row)=c0
        enddo
    
        ! Velocities at nodes
    
        do elem = 1, nx_elem_nh  !assembling rhs over elements (no halo)
    
           elnodes = elnode(1:3,elem)
    
           ! Derivatives
           dx  = bafux(:,elem)
           dy  = bafuy(:,elem)
           vol = voltriangle(elem)
           um=sum(uvel(elnode(1:3,elem)))
           vm=sum(vvel(elnode(1:3,elem)))
    
           ! Diffusivity
    
           !diff = ice_diff * sqrt( elem_area(elem) / scale_area )

            diff = ice_diff *sqrt(voltriangle(elem)/scalevol)
            
           do n = 1, 3
              row = elnodes(n)
               !row = elnodes(n,elem)
              do q = 1, 3
                 entries(q) = vol*dt_dyn*((dx(n)*(um+uvel(elnodes(q))) +      &
                              dy(n)*(vm+vvel(elnodes(q))))/12.0 -             &
                              diff*(dx(n)*dx(q)+ dy(n)*dy(q)) -               &
                              0.5*dt_dyn*(um*dx(n)+vm*dy(n))*(um*dx(q)+vm*dy(q))/9.0)
              enddo
              rhs_tr(row)=rhs_tr(row)+sum(entries*trc(elnode(1:3,elem)))
           enddo
        enddo
    
    end subroutine tg_rhs_icepack
    
    !=======================================================================
    
    module subroutine init_advection_icepack
    
        !use o_param
        !use o_mesh
        !use g_parsup
         use mice_module
         
        implicit none
        integer (kind=int_kind) :: ntrcr,narr

        real   (kind=dbl_kind)            :: fluxchan1,fluxchan2,suma,xctr2,yctr2,vnor1,vnor2,&
        & tmpx2,tmpx3,tmpy2,tmpy3,utmp,vtmp,xtmp,ytmp,vnorm,&
        & trctmp, trcmin, trcmax   
         integer (kind=int_kind) :: i,j,ie,id,id2,id3,isd3,isd2,     &
         m,n1,n2,jj,indx,nd
         logical (kind=log_kind) :: flag_noadv
         real(rkind) :: swild10(3,2),swild(3),swild2(3,2),trl0(nx),d_tr0(nx)

        call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
         narr   = 1 + ncat * (3 + ntrcr)
        !type(t_mesh), intent(in), target :: mesh
          
        ! Initialization of arrays necessary to implement FCT algorithm
        allocate(trl(nx))   ! low-order solutions
        allocate(d_tr(nx))  ! increments of high
                            ! order solutions
        allocate(icefluxes(nx_elem, 3))
        allocate(icepplus(nx), icepminus(nx))
        allocate(rhs_tr(nx),  rhs_trdiv(nx))
        allocate(newwork(nx))
        allocate(aream(nx))
        allocate(flux_matrix1(mnei_p,nx_nh,narr)) 
        allocate(flux_matrix2(mnei_p,nx_nh,narr))
        allocate(flux_matrix3(mnei_p,nx_nh,narr)) 
        allocate(flux_matrix4(mnei_p,nx_nh,narr)) 
        allocate(flux_matrix00(mnei_p,nx,2)) 
        allocate(flux_matrix_casulli(mnei_p,nx,2)) 
        allocate(flux_ratio_casulli(nx)) 
        allocate(flux_matrix0(3,mnei_p,nx,narr))  
        !allocate(mass_matrix(sum(nn_num(1:nx_nh))))

      
        trl(:)    = c0
        d_tr(:)   = c0
        rhs_tr(:) = c0
        rhs_trdiv(:)   = c0
        icefluxes(:,:) = c0
        icepplus(:)    = c0
        icepminus(:)   = c0
        newwork(:)     = c0 
        aream(:)     = c0  
        flux_matrix1(:,:,:) = c0
        flux_matrix2(:,:,:) = c0
        flux_matrix3(:,:,:) = c0
        flux_matrix4(:,:,:) = c0
        flux_matrix0(:,:,:,:) = c0
        !mass_matrix(:) = c0
      
        ! Fill in  the mass matrix
        !call fill_mass_matrix_icepack
        !mass_matrix=ice_matrix
        if (myrank==0) write(*,*) 'Icepack FCT is initialized'


        do i=1,npa
          aream(i) = 0
          do j=1,nne(i)
            ie=indel(j,i)
            if(ie<0) cycle
            !Wet elem
            
            id=iself(j,i)
            id3=nxq(i34(ie)-1,id,i34(ie)) !adjacent side index
            id2=nxq(i34(ie)-2,id,i34(ie)) !adjacent side index
            isd3=elside(id3,ie)
            isd2=elside(id2,ie)

            !Compute coord of side center and centroid (for ics=2)
            if(ics==1) then
              xctr2=xctr(ie); yctr2=yctr(ie)
              tmpx2=xcj(isd2); tmpy2=ycj(isd2)
              tmpx3=xcj(isd3); tmpy3=ycj(isd3)
            else !ll; use [xy]el defined in eframe
              xctr2=0.d0 !sum(xel(elnode(1:i34(ie),ie)))/i34(ie)
              yctr2=0.d0 !sum(yel(elnode(1:i34(ie),ie)))/i34(ie)
              tmpx3=(xel(id,ie)+xel(nxq(1,id,i34(ie)),ie))/2.d0
              tmpy3=(yel(id,ie)+yel(nxq(1,id,i34(ie)),ie))/2.d0
              tmpx2=(xel(id,ie)+xel(nxq(i34(ie)-1,id,i34(ie)),ie))/2.d0
              tmpy2=(yel(id,ie)+yel(nxq(i34(ie)-1,id,i34(ie)),ie))/2.d0
            endif !ics

            xtmp=yctr2-tmpy3
            ytmp=tmpx3-xctr2
            !write(12,*) i,np,iplg(i),j,nne(i),ielg(ie),id,ie,nea,xtmp,ytmp

            xtmp=tmpy2-yctr2
            ytmp=xctr2-tmpx2
            !write(12,*) i,np,iplg(i),j,nne(i),ielg(ie),id,ie,nea,xtmp,ytmp
         enddo

      enddo !i=1,np

      
        do i=1,nea
          do j=1,i34(i)
            ie=i
            !Wet elem
            
            id=j
            id3=nxq(i34(ie)-1,id,i34(ie)) !adjacent side index
            id2=nxq(i34(ie)-2,id,i34(ie)) !adjacent side index
            isd3=elside(id3,ie)
            isd2=elside(id2,ie)

            !Compute coord of side center and centroid (for ics=2)
            if(ics==1) then
              xctr2=xctr(ie); yctr2=yctr(ie)
              tmpx2=xcj(isd2); tmpy2=ycj(isd2)
              tmpx3=xcj(isd3); tmpy3=ycj(isd3)
            else !ll; use [xy]el defined in eframe
              xctr2=0.d0 !sum(xel(elnode(1:i34(ie),ie)))/i34(ie)
              yctr2=0.d0 !sum(yel(elnode(1:i34(ie),ie)))/i34(ie)
              tmpx3=(xel(j,ie)+xel(nxq(1,j,i34(ie)),ie))/2.d0
              tmpy3=(yel(j,ie)+yel(nxq(1,j,i34(ie)),ie))/2.d0
              tmpx2=(xel(j,ie)+xel(nxq(i34(ie)-1,j,i34(ie)),ie))/2.d0
              tmpy2=(yel(j,ie)+yel(nxq(i34(ie)-1,j,i34(ie)),ie))/2.d0
            endif !ics

            xtmp=yctr2-tmpy3
            ytmp=tmpx3-xctr2
            !write(12,*) elnode(j,i),npa,iplg(elnode(j,i)),i,nea,xtmp,ytmp

            xtmp=tmpy2-yctr2
            ytmp=xctr2-tmpx2
            !write(12,*) elnode(j,i),npa,iplg(elnode(j,i)),i,nea,xtmp,ytmp
         enddo

      enddo !i=1,ne
    
    end subroutine init_advection_icepack
    
    !=======================================================================
    
    subroutine fill_mass_matrix_icepack
    
        !use mod_mesh
        !use o_mesh
        !use i_param
        !use g_parsup
         use mice_module

        implicit none
      
        integer(kind=int_kind)                 :: n, n1, n2, row
        integer(kind=int_kind)                 :: elem, elnodes(3), q, offset, col, ipos
        integer(kind=int_kind), allocatable    :: col_pos(:)
        real(kind=dbl_kind)                    :: aa
        integer(kind=int_kind)                 :: flag=0 ,iflag=0
        !type(t_mesh), intent(in), target       :: mesh
      
!#include "../associate_mesh.h"
        
        !allocate(col_pos(nx))
          
        !do elem=1,nx_elem_nh
        !   elnodes=elnode(:,elem)
        !   do n = 1, 3
        !      row = elnodes(n)
        !      if ( row > nx_nh ) cycle
        !     ! Global-to-local neighbourhood correspondence
        !      do q = 1, nn_num(row)
        !         col_pos(nn_pos(q,row))=q
        !      enddo
        !      offset = ssh_stiff%rowptr(row) - ssh_stiff%rowptr(1)
        !      do q = 1, 3
        !         col  = elnodes(q)
        !         ipos = offset+col_pos(col)
        !         mass_matrix(ipos) = mass_matrix(ipos) + elem_area(elem) / 12.0_WP
        !         if ( q == n ) then
        !         mass_matrix(ipos) = mass_matrix(ipos) + elem_area(elem) / 12.0_WP
        !         end if
        !      enddo
        !   enddo
        !enddo
      
        ! TEST: area == sum of row entries in mass_matrix:
        !do q = 1, nx_nh
        !   offset = ssh_stiff%rowptr(q)   - ssh_stiff%rowptr(1) + 1
        !   n      = ssh_stiff%rowptr(q+1) - ssh_stiff%rowptr(1)
        !   aa     = sum(mass_matrix(offset:n))
        !   if ( abs(area(1,q)-aa) > p1) then
        !      iflag = q
        !      flag  = 1
        !   endif
        !enddo
      
        !if ( flag > 0 ) then
        !   offset = ssh_stiff%rowptr(iflag)   - ssh_stiff%rowptr(1)+1
        !   n      = ssh_stiff%rowptr(iflag+1) - ssh_stiff%rowptr(1)
        !   aa     = sum(mass_matrix(offset:n))
        ! 
         !  write(*,*) '#### MASS MATRIX PROBLEM', myrank, iflag, aa, area(1,iflag)
        !endif
        !deallocate(col_pos)
    
    end subroutine fill_mass_matrix_icepack
    
    !=======================================================================
    
    subroutine solve_low_order_icepack(trc)
    
        !============================
        ! Low-order solution
        !============================
        !
        ! It is assumed that m_ice, a_ice and m_snow from the previous time step
        ! are known at 1:myDim_nod2D+eDim_nod2D.
        ! We add diffusive contribution to the rhs. The diffusion operator
        ! is implemented as the difference between the consistent and lumped mass
        ! matrices acting on the field from the previous time step. The consistent
        ! mass matrix on the lhs is replaced with the lumped one.
    
        !use mod_mesh
        !use o_mesh
        !use i_param
        !use g_parsup
         use mice_module
      
        implicit none
      
        integer(kind=int_kind)           :: row,j, clo, clo2, cn, location(100)
        real   (kind=dbl_kind)           :: gamma,sum1
       !type(t_mesh),        target,        intent(in)     :: mesh
        real(kind=dbl_kind), dimension(nx), intent(inout)  :: trc
    
!#include "../associate_mesh.h"
      
        gamma = ice_gamma_fct       ! Added diffusivity parameter
                                    ! Adjust it to ensure posivity of solution
        trl(:)    = c0

        do row = 1, nx_nh
         sum1=ice_matrix(0,row)*trc(row)
            do j=1,nnp(row)
               clo=indnd(j,row)
               sum1=sum1+ice_matrix(j,row)*trc(clo)
            enddo !j
           !clo  = ssh_stiff%rowptr(row)   - ssh_stiff%rowptr(1) + 1
           !clo2 = ssh_stiff%rowptr(row+1) - ssh_stiff%rowptr(1)
           !cn   = clo2 - clo + 1
           !location(1:cn) = nn_pos(1:cn, row)
           trl(row) = (rhs_tr(row) + gamma * sum1) / lump_ice_matrix(row) +              &
                      (1.0-gamma) * trc(row)
        enddo
      
        call exchange_p2d(trl)
      
        ! Low-order solution must be known to neighbours
    
    end subroutine solve_low_order_icepack
    
    !=======================================================================
    
    subroutine solve_high_order_icepack(trc)
    
        !use mod_mesh
        !use o_mesh
        !use i_param
        !use g_parsup
         use mice_module
    
        implicit none
      
        integer(kind=int_kind)             :: n,i,j,clo,clo2,cn,location(100),row
        real   (kind=dbl_kind)          :: rhs_new,sum1
        integer(kind=int_kind), parameter  :: num_iter_solve = 3
        !type(t_mesh),        target,        intent(in)     :: mesh
        real(kind=dbl_kind), dimension(nx), intent(inout)  :: trc
      
!#include "../associate_mesh.h"
      
        ! Taylor-Galerkin solution
       
        ! First guess solution

        d_tr(:)   = c0
        trl(:)    = c0

        do row = 1, nx_nh
           d_tr(row) = rhs_tr(row) / lump_ice_matrix(row)
        end do
      
        call exchange_p2d(d_tr)
      
        ! Iterate
        do n = 1, num_iter_solve - 1
           do row = 1, nx_nh
              !clo  = ssh_stiff%rowptr(row) - ssh_stiff%rowptr(1) + 1
              !clo2 = ssh_stiff%rowptr(row+1) - ssh_stiff%rowptr(1)
              !cn   = clo2 - clo + 1
              !location(1:cn) = nn_pos(1:cn,row)
               sum1=ice_matrix(0,row)*d_tr(row)
               do j=1,nnp(row)
                  clo=indnd(j,row)
                  sum1=sum1+ice_matrix(j,row)*d_tr(clo)
               enddo !j
              rhs_new  = rhs_tr(row) - sum1
              trl(row) = d_tr(row) + rhs_new / lump_ice_matrix(row)
           enddo
           do row = 1, nx_nh
              d_tr(row) = trl(row)
           enddo
           call exchange_p2d(d_tr)
        enddo
      
    end subroutine solve_high_order_icepack
    
    !=======================================================================
    
    subroutine fem_fct_icepack(trc)
    
        !============================
        ! Flux corrected transport algorithm for tracer advection
        !============================
        !
        ! It is based on Loehner et al. (Finite-element flux-corrected
        ! transport (FEM-FCT) for the Euler and Navier-Stokes equation,
        ! Int. J. Numer. Meth. Fluids, 7 (1987), 1093--1109) as described by Kuzmin and
        ! Turek. (kuzmin@math.uni-dortmund.de)
    
        !use mod_mesh
        !use o_mesh
        !use o_param
        !use i_param
        !use g_parsup
        use mice_module
   
    
        integer(kind=int_kind)                            :: icoef(3,3), n, q, elem, elnodes(3), row
        real   (kind=dbl_kind), allocatable, dimension(:) :: tmax, tmin, trctmp
        real   (kind=dbl_kind)                            :: vol, flux, ae, gamma, tempmax
        !type(t_mesh),        target,        intent(in)     :: mesh  
        !nx=npa, nx_nh=np
        real(kind=dbl_kind), dimension(nx), intent(inout)  :: trc
    
!#include "../associate_mesh.h"
      
        gamma = ice_gamma_fct        ! It should coinside with gamma in
                                     ! ts_solve_low_order
    
        !==========================
        ! Compute elemental antidiffusive fluxes to nodes
        !==========================
        ! This is the most unpleasant part ---
        ! it takes memory and time. For every element
        ! we need its antidiffusive contribution to
        ! each of its 3 nodes
       
        allocate(tmax(nx_nh), tmin(nx_nh))
        allocate(trctmp(nx))
       
        ! Auxiliary elemental operator (mass matrix- lumped mass matrix)
       
        icoef = 1
        icefluxes(:,:) = c0
        trctmp(:) = trc(:)
        do n = 1, 3   ! three upper nodes
                      ! Cycle over rows  row=elnodes(n)
                 icoef(n,n) = -2
        enddo
       
        do elem = 1, nx_elem
           elnodes = elnode(1:3,elem)
           vol     = voltriangle(elem)
           do q = 1, 3
              icefluxes(elem,q) = -sum(icoef(:,q) * (gamma * trc(elnodes) +       &
                                  d_tr(elnodes))) * (vol / lump_ice_matrix(elnodes(q))) / 12.0
            !if(trc(elnodes(q))>10) write(12,*) elem,elnodes(q),trc(elnodes),trl(elnodes),d_tr(elnodes),vol,lump_ice_matrix(elnodes(q)),icefluxes(elem,q)
            !if(abs(icefluxes(elem,q))>0) write(12,*) elem,elnodes(q),trc(elnodes),trl(elnodes),d_tr(elnodes),icefluxes(elem,q)
           enddo
        enddo
       
        !==========================
        ! Screening the low-order solution
        !==========================
        ! TO BE ADDED IF FOUND NECESSARY
        ! Screening means comparing low-order solutions with the
        ! solution on the previous time step and using whichever
        ! is greater/smaller in computations of max/min below
       
        !==========================
        ! Cluster min/max
        !==========================
       
        do row = 1, nx_nh
           n = nnp(row)
            tmax(row) = maxval(trl(indnd(1:n,row)))
            !tempmax   = maxval(trc(indnd(1:n,row)))
            tmax(row) = max(tmax(row),trl(row))
            !tmax(row) = max(trc(row),tmax(row))
            !tmax(row) = max(tempmax,tmax(row))
            tmin(row) = minval(trl(indnd(1:n,row)))
            !tempmax   = minval(trc(indnd(1:n,row)))
            tmin(row) = min(tmin(row),trl(row))
            !tmin(row) = min(trc(row),tmin(row))
            !tmin(row) = min(tempmax,tmin(row))
           ! Admissible increments
           tmax(row) = tmax(row) - trl(row)
           tmin(row) = tmin(row) - trl(row)
        enddo
       
        !=========================
        ! Sums of positive/negative fluxes to node row
        !=========================
       
        icepplus = c0
        icepminus = c0
        do elem = 1, nx_elem
           elnodes = elnode(1:3,elem)
           do q = 1, 3
              n    = elnodes(q)
              flux = icefluxes(elem,q)
              if ( flux > 0 ) then
                 icepplus(n) = icepplus(n) + flux
              else
                 icepminus(n) = icepminus(n) + flux
              endif
           enddo
        enddo
       
        !========================
        ! The least upper bound for the correction factors
        !========================
    
        do n = 1, nx_nh
           flux = icepplus(n)
           if ( abs(flux) > 0 ) then
              icepplus(n) = min(1.0,tmax(n) / flux)
           else
              icepplus(n) = c0
           endif
       
           flux = icepminus(n)
           if ( abs(flux) > 0 ) then
              icepminus(n) = min(1.0,tmin(n) / flux)
           else
              icepminus(n) = c0
           endif
         enddo
       
         ! pminus and pplus are to be known to neighbouting PE
         call exchange_p2d(icepminus)
         call exchange_p2d(icepplus)
       
         !========================
         ! Limiting
         !========================
       
         do elem = 1, nx_elem
            elnodes = elnode(1:3,elem)
            ae = c1 !/ 10
            do q = 1, 3
               n    = elnodes(q)
               flux = icefluxes(elem,q)
               !if (sqrt(uvel(n)**2+vvel(n)**2).eq.0) ae=0
               if ( flux >=c0 ) ae = min(ae, icepplus(n))
               if ( flux < c0 ) ae = min(ae, icepminus(n))
               !if ( flux >=c0 ) then
               !   tempmax = minval(icepplus(indnd(1:nnp(n),n)))
               !   tempmax = min(icepplus(n),tempmax)
               !   ae = min(ae,tempmax)
               !endif
               !if ( flux < c0 ) then
               !   tempmax = minval(icepminus(indnd(1:nnp(n),n)))
               !   tempmax = min(icepminus(n),tempmax)
               !   ae = min(ae,tempmax)
               !endif
               !if(trc(elnodes(q))<10.and.trc(elnodes(q))>5) write(12,*) elem,n,trc(n),trl(n),rhs_tr(n)/lump_ice_matrix(n),d_tr(n),ae,icepplus(n),icepminus(n),flux
            enddo
            icefluxes(elem,:) = ae * icefluxes(elem,:)
         enddo
       
       
         !==========================
         ! Update the solution
         !==========================
         do elem = nx_elem_nh, nx_elem
            elnodes = elnode(1:3,elem)
            do q = 1, 3
               n = elnodes(q)
               !if(trc(elnodes(q))<20.and.trc(elnodes(q))>5) write(12,*) n, trc(n), trl(n), icefluxes(elem,q)
            enddo
         enddo
          do n = 1, nx_nh
             trc(n) = trl(n)
          end do
          call exchange_p2d(trc)
         
          do elem = 1, nx_elem_nh
             elnodes = elnode(1:3,elem)
             do q = 1, 3
                n = elnodes(q)
                trc(n) = trc(n) + icefluxes(elem,q)
             enddo
          enddo
          call exchange_p2d(trc)

         ! do n = 1, nx_nh
         !    if(trc(n)*trctmp(n)<0) trc(n) = trctmp(n)
         ! end do
          call exchange_p2d(trc)
       
          deallocate(tmin, tmax)
          deallocate(trctmp)
    
    end subroutine fem_fct_icepack
    
    !=======================================================================
    
    subroutine tg_rhs_div_icepack(trc)
    
        !use mod_mesh
        !use o_mesh
        !use o_param
        !use i_param
        !use g_parsup
        use mice_module
    
        implicit none
    
        real   (kind=dbl_kind)    :: diff, entries(3), um, vm, vol, dx(3), dy(3)
        integer(kind=int_kind)    :: n, q, row, elem, elnodes(3)
        real   (kind=dbl_kind)    :: c_1, c_2, c_3, c_4, c_x, entries2(3),entries3(3),utmp(3),vtmp(3)
        !type(t_mesh),        target,        intent(in)     :: mesh
        real(kind=dbl_kind), dimension(nx), intent(inout)  :: trc    

!#include "../associate_mesh.h"
    
        ! Computes the rhs in a Taylor-Galerkin way (with urrayspwind 
        ! type of correction for the advection operator).
        ! In this version I tr to split divergent term off, 
        ! so that FCT works without it.
    
        do row = 1, nx !=npa
           ! row=myList_nod2D(m)
           rhs_tr(row)    = c0
           rhs_trdiv(row) = c0
        enddo

        do elem = 1, nx_elem   !! assembling rhs over elements
                                  !! elem=myList_elem2D(m)

           elnodes = elnode(1:3,elem)
            do n = 1, 3
               !transform pframe to eframe
               if(ics==1) then
                  utmp(n) = uvel(elnodes(n))
                  vtmp(n) = vvel(elnodes(n))
               else
                  utmp(n) = uvel(elnodes(n)) * dot_product(pframe(1:3,1,elnodes(n)),eframe(1:3,1,elem)) + vvel(elnodes(n)) * dot_product(pframe(1:3,2,elnodes(n)),eframe(1:3,1,elem))
                  vtmp(n) = uvel(elnodes(n)) * dot_product(pframe(1:3,1,elnodes(n)),eframe(1:3,2,elem)) + vvel(elnodes(n)) * dot_product(pframe(1:3,2,elnodes(n)),eframe(1:3,2,elem))
               endif
              
           enddo
           ! Derivatives
           dx  = bafux(:,elem)
           dy  = bafuy(:,elem)
           vol = voltriangle(elem)
           um  = sum(utmp)
           vm  = sum(vtmp)
      
           ! This is exact computation (no assumption of u=const 
           ! on elements used in the standard version)
           c_1 = (um*um+sum(utmp(1:3)*utmp(1:3))) / 12.0_dbl_kind
           c_2 = (vm*vm+sum(vtmp(1:3)*vtmp(1:3))) / 12.0_dbl_kind
           c_3 = (um*vm+sum(vtmp(1:3)*utmp(1:3))) / 12.0_dbl_kind
           c_4 = sum(dx*utmp(1:3)+dy*vtmp(1:3))
      
           do n = 1, 3
              row = elnodes(n)
      
              do q = 1, 3
                 entries(q)  = vol*dt_dyn*((c1-p5*dt_dyn*c_4)*(dx(n)*(um+utmp(q))+ &
                               dy(n)*(vm+vtmp(q)))/12.0_dbl_kind                 - &
                               p5*dt_dyn*(c_1*dx(n)*dx(q)+c_2*dy(n)*dy(q)+c_3*(dx(n)*dy(q)+dx(q)*dy(n))))
                 entries2(q) = p5*dt_dyn*(dx(n)*(um+utmp(q))                     + &
                               dy(n)*(vm+vtmp(q))-dx(q)*(um+utmp(n))      - &
                               dy(q)*(vm+vtmp(n)))
                 entries3(q)=  vol*dt_dyn*(((dx(n)*(um+utmp(q)))+ dy(n)*(vm+vtmp(q)))/12.0_dbl_kind &
                               -p5*dt_dyn*(um*dx(n)+vm*dy(n))*(um*dx(q)+vm*dy(q))/9.0)

              enddo
              c_x = vol*dt_dyn*c_4*(sum(trc(elnodes))+trc(elnodes(n))+sum(entries2*trc(elnodes))) / 12.0_dbl_kind
              rhs_tr(row)    = rhs_tr(row) + sum(entries * trc(elnodes)) + c_x
              !rhs_tr(row)    = rhs_tr(row) + sum(entries3 * trc(elnodes)) 
              rhs_trdiv(row) = rhs_trdiv(row) - c_x
           enddo
        enddo
    
        call exchange_p2d(rhs_tr)
        call exchange_p2d(rhs_trdiv)
    end subroutine tg_rhs_div_icepack
    
    !=======================================================================
    
    subroutine update_for_div_icepack(trc)
    
        !use mod_mesh
        !use o_mesh
        !use o_param
        !use i_param
        !use g_parsup
        use mice_module
    
    
        implicit none
    
        integer(kind=int_kind)            :: n, i, clo, clo2, cn, &
                                             location(100), row,j
        real   (kind=dbl_kind)            :: rhs_new,sum1
        integer(kind=int_kind), parameter :: num_iter_solve=3
        !type(t_mesh),        target,        intent(in)     :: mesh
        real(kind=dbl_kind), dimension(nx), intent(inout)  :: trc
    
!#include "../associate_mesh.h"
    
        ! Computes Taylor-Galerkin solution
        ! first approximation

        d_tr(:)   = c0
        trl(:)    = c0

        do row = 1, nx_nh
           d_tr(row) = rhs_trdiv(row) / lump_ice_matrix(row)
        enddo
      
        call exchange_p2d(d_tr)
      
        ! Iterate
      
        do n = 1, num_iter_solve-1
           do row = 1, nx_nh
              !clo            = ssh_stiff%rowptr(row)   - ssh_stiff%rowptr(1) + 1
              !clo2           = ssh_stiff%rowptr(row+1) - ssh_stiff%rowptr(1)
              !cn             = clo2 - clo + 1
              !location(1:cn) = nn_pos(1:cn, row)
               sum1=ice_matrix(0,row)*d_tr(row)
            do j=1,nnp(row)
               clo = indnd(j,row)
               sum1=sum1+ice_matrix(j,row)*d_tr(clo)
            enddo !j
              rhs_new        = rhs_trdiv(row) - sum1
              trl(row)       = d_tr(row) + rhs_new / lump_ice_matrix(row)
           enddo
           do row = 1, nx_nh
              d_tr(row) = trl(row)
           enddo
           call exchange_p2d(d_tr)
        enddo
      
        trc = trc + d_tr

        call exchange_p2d(trc)
    
    end subroutine update_for_div_icepack 
    
    !=======================================================================
    
    subroutine fct_solve_icepack(trc)
        
        !use mod_mesh     
     
        implicit none

        real(kind=dbl_kind), dimension(nx), intent(inout) :: trc      
        !type(t_mesh),        target,        intent(in)    :: mesh
      
        ! Driving sequence
        call tg_rhs_div_icepack(trc)
        call solve_high_order_icepack(trc) ! uses arrays of low-order solutions as temp
                                                 ! storage. It should preceed the call of low
                                                 ! order solution.
        call solve_low_order_icepack(trc)
        call fem_fct_icepack(trc)
        call update_for_div_icepack(trc)

        trl(:)    = c0
        d_tr(:)   = c0
        rhs_tr(:) = c0
        rhs_trdiv(:)   = c0
        icefluxes(:,:) = c0

    end subroutine fct_solve_icepack

!=======================================================================
    
    subroutine tg_rhs_div_icepack_upwind(trc,trc_n)
    
        !use mod_mesh
        !use o_mesh
        !use o_param
        !use i_param
        !use g_parsup
        use mice_module
    
        implicit none
        integer (kind=int_kind), intent(in) :: trc_n
        real   (kind=dbl_kind)    :: diff, entries(3), um, vm, vol, dx(3), dy(3)
        integer(kind=int_kind)    :: n, q, row, elem, m, jj, indx, elnodes(3)
        real   (kind=dbl_kind)    :: c_1, c_2, c_3, c_4, c_x, entries2(3),entries3(3),utmp(3),vtmp(3),sumtmp
        !type(t_mesh),        target,        intent(in)     :: mesh
        real(kind=dbl_kind), dimension(nx), intent(inout)  :: trc    

!#include "../associate_mesh.h"
    
        ! Computes the rhs in a Taylor-Galerkin way (with urrayspwind 
        ! type of correction for the advection operator).
        ! In this version I tr to split divergent term off, 
        ! so that FCT works without it.
    
        !do row = 1, nx !=npa
           ! row=myList_nod2D(m)
        !   rhs_tr(row)    = c0
        !   rhs_trdiv(row) = c0
        !enddo

        do elem = 1, nx_elem   !! assembling rhs over elements
                                  !! elem=myList_elem2D(m)

           elnodes = elnode(1:3,elem)
            do n = 1, 3
               !transform pframe to eframe
                if(ics==1) then
                  utmp(n) = uvel(elnodes(n))
                  vtmp(n) = vvel(elnodes(n))
               else
                  utmp(n) = uvel(elnodes(n)) * dot_product(pframe(1:3,1,elnodes(n)),eframe(1:3,1,elem)) + vvel(elnodes(n)) * dot_product(pframe(1:3,2,elnodes(n)),eframe(1:3,1,elem))
                  vtmp(n) = uvel(elnodes(n)) * dot_product(pframe(1:3,1,elnodes(n)),eframe(1:3,2,elem)) + vvel(elnodes(n)) * dot_product(pframe(1:3,2,elnodes(n)),eframe(1:3,2,elem))
               endif
           enddo
           ! Derivatives
           dx  = bafux(:,elem)
           dy  = bafuy(:,elem)
           vol = voltriangle(elem)
           um  = sum(utmp)
           vm  = sum(vtmp)
      
           ! This is exact computation (no assumption of u=const 
           ! on elements used in the standard version)
           c_1 = (um*um+sum(utmp(1:3)*utmp(1:3))) / 12.0_dbl_kind
           c_2 = (vm*vm+sum(vtmp(1:3)*vtmp(1:3))) / 12.0_dbl_kind
           c_3 = (um*vm+sum(vtmp(1:3)*utmp(1:3))) / 12.0_dbl_kind
           c_4 = sum(dx*utmp(1:3)+dy*vtmp(1:3))
      
           do n = 1, 3
              row = elnodes(n)
      
              do q = 1, 3
                 entries(q)  = vol*dt_dyn*((c1-p5*dt_dyn*c_4)*(dx(n)*(um+utmp(q))+ &
                               dy(n)*(vm+vtmp(q)))/12.0_dbl_kind                 - &
                               p5*dt_dyn*(c_1*dx(n)*dx(q)+c_2*dy(n)*dy(q)+c_3*(dx(n)*dy(q)+dx(q)*dy(n))))
                 entries2(q) = p5*dt_dyn*(dx(n)*(um+utmp(q))                     + &
                               dy(n)*(vm+vtmp(q))-dx(q)*(um+utmp(n))      - &
                               dy(q)*(vm+vtmp(n)))
                 entries3(q)=  vol*dt_dyn*(((dx(n)*(um+utmp(q)))+ dy(n)*(vm+vtmp(q)))/12.0_dbl_kind &
                               -p5*dt_dyn*(um*dx(n)+vm*dy(n))*(um*dx(q)+vm*dy(q))/9.0)

              enddo
              c_x = vol*dt_dyn*c_4*(sum(trc(elnodes))+trc(elnodes(n))+sum(entries2*trc(elnodes))) / 12.0_dbl_kind
              !rhs_tr(row)    = rhs_tr(row) + sum(entries * trc(elnodes)) + c_x
              !rhs_tr(row)    = rhs_tr(row) + sum(entries3 * trc(elnodes)) 
               indx=0
                  do m=1,nne(row)
                  if(indel(m,row)==elem) then
                     indx=m; exit
                  endif
                  enddo !m
                  if(indx==0) call parallel_abort('STEP: failed to find (tg_rhs_div_icepack_upwind)')  
              flux_matrix0(1,indx,row,trc_n) = entries3(1) * trc(elnodes(1))/ lump_ice_matrix(row)
              flux_matrix0(2,indx,row,trc_n) = entries3(2) * trc(elnodes(2))/ lump_ice_matrix(row)
              flux_matrix0(3,indx,row,trc_n) = entries3(3) * trc(elnodes(3))/ lump_ice_matrix(row)
              !sumtmp = sum(flux_matrix0(:,indx,row,trc_n))
              !write(12,*) row,indx,sumtmp
              !rhs_trdiv(row) = rhs_trdiv(row) - c_x
           enddo
        enddo
    

    end subroutine tg_rhs_div_icepack_upwind
    
    !=======================================================================
    
    subroutine upwind_icepack_bynode(trc,trc_n)
        
        !use mod_mesh     
         use mice_module
        implicit none
        integer (kind=int_kind), intent(in) :: trc_n
        real(kind=dbl_kind), dimension(nx), intent(inout) :: trc  
        real   (kind=dbl_kind)            :: suma,xctr2,yctr2,vnor1,vnor2,sumtmp,areatmp,&
                                             & tmpx2,tmpx3,tmpy2,tmpy3,utmp,vtmp,xtmp,ytmp,vnorm,&
                                             & trctmp, trcmin, trcmax ,puny  
        integer (kind=int_kind) :: i,j,ie,id,id2,id3,isd3,isd2,     &
                                   m,n1,n2,jj,indx,nd,nodeid
        logical (kind=log_kind) :: flag_noadv
         real(rkind) :: swild10(3,2),swild(3),swild2(3,2),trl0(nx),d_tr0(nx),fluxchan1(nx),fluxchan2(nx)
        !type(t_mesh),        target,        intent(in)    :: mesh
         call icepack_query_parameters(puny_out=puny)

      trl0(:) = c0
      fluxchan1(:) = c0 !positive flux
      fluxchan2(:) = c0 !negative flux
      d_tr0(:) = trc(:)
      flux_matrix1(:,:,trc_n) = c0
      flux_matrix2(:,:,trc_n) = c0
      flux_matrix3(:,:,trc_n) = c0
      flux_matrix4(:,:,trc_n) = c0
      flux_matrix0(:,:,:,trc_n) = c0
        ! Driving sequence
      do i=1,nx_nh
         fluxchan2(i)=0.d0 !sum of fluxes;
         fluxchan1(i)=0.d0 !sum of fluxes;\int_\Gamma u*u_n d\Gamma
        if(idry(i)==1) cycle
        !if(isbnd(1,i)/=0) cycle
        !if(any(isbnd(1,indnd(1:nnp(i),i))/=0)) cycle
        !Wet node
          
          
          suma=0.d0 !sum of areas
          aream(i) = 0
          do j=1,nne(i)
            ie=indel(j,i)
            if(idry_e(ie)==1) cycle
            areatmp = area(ie)/6
            suma=suma+area(ie)/3 !approx for quad
            !Wet elem
            
            id=iself(j,i)
            id3=nxq(i34(ie)-1,id,i34(ie)) !adjacent side index
            id2=nxq(i34(ie)-2,id,i34(ie)) !adjacent side index
            isd3=elside(id3,ie)
            isd2=elside(id2,ie)

            do m=1,i34(ie) !side
                n1=isidenode(1,elside(m,ie))
                n2=isidenode(2,elside(m,ie))
                !swild10(m,1) = 0.5*(uvel(n1)*trc(n1)*dot_product(pframe(:,1,n1),sframe2(:,1,elside(m,ie)))+vvel(n1)*trc(n1)*dot_product(pframe(:,2,n1),sframe2(:,1,elside(m,ie)))&
                !                 &+uvel(n2)*trc(n2)*dot_product(pframe(:,1,n2),sframe2(:,1,elside(m,ie)))+vvel(n2)*trc(n2)*dot_product(pframe(:,2,n1),sframe2(:,1,elside(m,ie))))
                !swild10(m,2) = 0.5*(uvel(n1)*trc(n1)*dot_product(pframe(:,1,n1),sframe2(:,2,elside(m,ie)))+vvel(n1)*trc(n1)*dot_product(pframe(:,2,n1),sframe2(:,2,elside(m,ie)))&
                !                 &+uvel(n2)*trc(n2)*dot_product(pframe(:,1,n2),sframe2(:,2,elside(m,ie)))+vvel(n2)*trc(n2)*dot_product(pframe(:,2,n1),sframe2(:,2,elside(m,ie))))
               if(ics==1) then
                  swild10(m,1) = 0.5*(uvel(n1)+uvel(n2))
                  swild10(m,2) = 0.5*(vvel(n1)+vvel(n2))
                  swild(m) = (trc(n1)+trc(n2))*0.5
                  swild2(m,1) = uvel(elnode(m,ie))
                  swild2(m,2) = vvel(elnode(m,ie))
               else
                swild10(m,1) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),sframe2(:,1,elside(m,ie)))+vvel(n1)*dot_product(pframe(:,2,n1),sframe2(:,1,elside(m,ie)))&
                                 &+uvel(n2)*dot_product(pframe(:,1,n2),sframe2(:,1,elside(m,ie)))+vvel(n2)*dot_product(pframe(:,2,n2),sframe2(:,1,elside(m,ie))))
                swild10(m,2) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),sframe2(:,2,elside(m,ie)))+vvel(n1)*dot_product(pframe(:,2,n1),sframe2(:,2,elside(m,ie)))&
                                 &+uvel(n2)*dot_product(pframe(:,1,n2),sframe2(:,2,elside(m,ie)))+vvel(n2)*dot_product(pframe(:,2,n2),sframe2(:,2,elside(m,ie))))
                !swild10(m,1) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),eframe(:,1,ie))+vvel(n1)*dot_product(pframe(:,2,n1),eframe(:,1,ie))&
                !                 &+uvel(n2)*dot_product(pframe(:,1,n2),eframe(:,1,ie))+vvel(n2)*dot_product(pframe(:,2,n2),eframe(:,1,ie)))
                !swild10(m,2) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),eframe(:,2,ie))+vvel(n1)*dot_product(pframe(:,2,n1),eframe(:,2,ie))&
                !                 &+uvel(n2)*dot_product(pframe(:,1,n2),eframe(:,2,ie))+vvel(n2)*dot_product(pframe(:,2,n2),eframe(:,2,ie)))                      
                swild(m) = (trc(n1)+trc(n2))*0.5
                swild2(m,1) = uvel(elnode(m,ie))*dot_product(pframe(:,1,elnode(m,ie)),eframe(:,1,ie))+vvel(elnode(m,ie))*dot_product(pframe(:,2,elnode(m,ie)),eframe(:,1,ie))
                swild2(m,2) = uvel(elnode(m,ie))*dot_product(pframe(:,1,elnode(m,ie)),eframe(:,2,ie))+vvel(elnode(m,ie))*dot_product(pframe(:,2,elnode(m,ie)),eframe(:,2,ie))
               endif
              
            enddo
            !utmp = sum(swild10(1:i34(ie),1))/3 !vel @ centroid
            !vtmp = sum(swild10(1:i34(ie),2))/3 !vel @ centroid
            utmp=sum(swild2(1:i34(ie),1))/3 !vel @ centroid
            vtmp=sum(swild2(1:i34(ie),2))/3 !vel @ centroid
            trctmp = sum(swild(1:i34(ie)))/3
            trcmin = minval(trc(indnd(1:nnp(i),i)))
            trcmin = min(trcmin,trc(i))
            trcmax = maxval(trc(indnd(1:nnp(i),i)))
            trcmax = max(trcmax,trc(i))

            !Compute coord of side center and centroid (for ics=2)
            if(ics==1) then
              xctr2=xctr(ie); yctr2=yctr(ie)
              tmpx2=xcj(isd2); tmpy2=ycj(isd2)
              tmpx3=xcj(isd3); tmpy3=ycj(isd3)
            else !ll; use [xy]el defined in eframe
              xctr2=0.d0 !sum(xel(elnode(1:i34(ie),ie)))/i34(ie)
              yctr2=0.d0 !sum(yel(elnode(1:i34(ie),ie)))/i34(ie)
              tmpx3=(xel(id,ie)+xel(nxq(1,id,i34(ie)),ie))/2.d0
              tmpy3=(yel(id,ie)+yel(nxq(1,id,i34(ie)),ie))/2.d0
              tmpx2=(xel(id,ie)+xel(nxq(i34(ie)-1,id,i34(ie)),ie))/2.d0
              tmpy2=(yel(id,ie)+yel(nxq(i34(ie)-1,id,i34(ie)),ie))/2.d0
            endif !ics

            !1st segment
            !Normal dir x length
!            xtmp=yctr(ie)-ycj(isd3)
!            ytmp=xcj(isd3)-xctr(ie)
            xtmp=yctr2-tmpy3
            ytmp=tmpx3-xctr2
            vnor1=utmp*xtmp+vtmp*ytmp !normal vel x length 
            vnor2=swild10(id3,1)*xtmp+swild10(id3,2)*ytmp !normal vel@side x length 
            !fluxchan1=fluxchan1+(utmp*vnor1+swild10(id3,1)*vnor2)/2.d0
            !fluxchan2=fluxchan2+(vtmp*vnor1+swild10(id3,2)*vnor2)/2.d0
            !fluxchan1=fluxchan1+dt_dyn*(swild(id3)*vnor2)
            sumtmp = dt_dyn*(trctmp*vnor1+swild(id3)*vnor2)/2.d0
            nodeid = elnode(nxq(1,id,i34(ie)),ie)
            !fluxchan1=fluxchan1+dt_dyn*(trctmp*vnor1+swild(id3)*vnor2)/2.d0
            !fluxchan1=fluxchan1+dt_dyn*(trctmp+swild(id3))*(vnor1+vnor2)/4.d0
            !fluxchan2=fluxchan2+dt_dyn*(trctmp*vnor1+swild(id3)*vnor2)/2.d0
            flux_matrix1(j,i,trc_n) = dt_dyn*trctmp*vnor1
            flux_matrix2(j,i,trc_n) = dt_dyn*swild(id3)*vnor2
            if(sumtmp>0) then
               fluxchan2(i) = fluxchan2(i) + sumtmp
               if(sumtmp/areatmp>d_tr0(i)*0.5)then
                  fluxchan2(i) = fluxchan2(i) + areatmp*d_tr0(i)
                  flux_matrix1(j,i,trc_n) = flux_matrix1(j,i,trc_n)*areatmp*d_tr0(i)/sumtmp*0.5
                  flux_matrix2(j,i,trc_n) = flux_matrix2(j,i,trc_n)*areatmp*d_tr0(i)/sumtmp*0.5
               endif
               
            else
               fluxchan1(i) = fluxchan1(i) + sumtmp
               if(sumtmp/areatmp>d_tr0(nodeid)*0.5)then
                  fluxchan1(i) = fluxchan1(i) + areatmp*d_tr0(nodeid)
                  flux_matrix1(j,i,trc_n) = flux_matrix1(j,i,trc_n)*areatmp*d_tr0(nodeid)/sumtmp*0.5
                  flux_matrix2(j,i,trc_n) = flux_matrix2(j,i,trc_n)*areatmp*d_tr0(nodeid)/sumtmp*0.5
               endif
            endif
            !if(isbs(isd3)>0) then !open bnd
            !  !vnorm=swild10(id3,1)*sframe(1,1,isd3)+swild10(id3,2)*sframe(2,1,isd3) !outer normal vel
            !  vnorm=swild(id3)*swild10(id3,1)*snx(isd3)+swild(id3)*swild10(id3,2)*sny(isd3) !outer normal vel
            !  fluxchan1=fluxchan1+dt_dyn*vnorm*distj(isd3)/2.d0
            !  fluxchan2=fluxchan2+dt_dyn*vnorm*distj(isd3)/2.d0
            !endif !isbs>0
            
            !2nd segment
            !Normal dir x length
!            xtmp=ycj(isd2)-yctr(ie)
!            ytmp=xctr(ie)-xcj(isd2)
            xtmp=tmpy2-yctr2
            ytmp=xctr2-tmpx2
            vnor1=utmp*xtmp+vtmp*ytmp !normal vel x length
            vnor2=swild10(id2,1)*xtmp+swild10(id2,2)*ytmp !normal vel x length


            !fluxchan1=fluxchan1+dt_dyn*(swild(id2)*vnor2)
            !fluxchan1=fluxchan1+dt_dyn*(trctmp*vnor1+swild(id2)*vnor2)/2.d0
            !fluxchan1=fluxchan1+dt_dyn*(trctmp+swild(id2))*(vnor1+vnor2)/4.d0
            !fluxchan2=fluxchan2+dt_dyn*(trctmp*vnor1+swild(id2)*vnor2)/2.d0
            sumtmp = dt_dyn*(trctmp*vnor1+swild(id2)*vnor2)/2.d0
            nodeid = elnode(nxq(i34(ie)-1,id,i34(ie)),ie)
            flux_matrix3(j,i,trc_n) = dt_dyn*trctmp*vnor1
            flux_matrix4(j,i,trc_n) = dt_dyn*swild(id2)*vnor2
            if(sumtmp>0) then
               fluxchan2(i) = fluxchan2(i) + sumtmp
               if(sumtmp/areatmp>d_tr0(i)*0.5)then
                  fluxchan2(i) = fluxchan2(i) + areatmp*d_tr0(i)
                  flux_matrix3(j,i,trc_n) = flux_matrix3(j,i,trc_n)*areatmp*d_tr0(i)/sumtmp*0.5
                  flux_matrix4(j,i,trc_n) = flux_matrix4(j,i,trc_n)*areatmp*d_tr0(i)/sumtmp*0.5
               endif
               
            else
               fluxchan1(i) = fluxchan1(i) + sumtmp
               if(sumtmp/areatmp>d_tr0(nodeid)*0.5)then
                  fluxchan1(i) = fluxchan1(i) + areatmp*d_tr0(nodeid)
                  flux_matrix3(j,i,trc_n) = flux_matrix3(j,i,trc_n)*areatmp*d_tr0(nodeid)/sumtmp*0.5
                  flux_matrix4(j,i,trc_n) = flux_matrix4(j,i,trc_n)*areatmp*d_tr0(nodeid)/sumtmp*0.5
               endif
            endif
            !if(isbs(isd2)>0) then !open bnd
            ! !vnorm=swild10(id2,1)*sframe(1,1,isd2)+swild10(id2,2)*sframe(2,1,isd2) !outer normal
            !  vnorm=swild(id2)*swild10(id2,1)*snx(isd2)+swild(id2)*swild10(id2,2)*sny(isd2) !outer normal
            !  fluxchan1=fluxchan1+dt_dyn*vnorm*distj(isd2)/2.d0
            !  fluxchan2=fluxchan2+dt_dyn*vnorm*distj(isd2)/2.d0
            !endif !isbs>0

            !if(ie>ne) then
            !   flux_matrix1(j,i,trc_n) = 0
            !   flux_matrix2(j,i,trc_n) = 0
            !   flux_matrix3(j,i,trc_n) = 0
            !   flux_matrix4(j,i,trc_n) = 0
            !endif

          enddo !j=1,nne(i)
          if(suma/=0) then
            trl0(i)=(fluxchan1(i)+fluxchan2(i))/suma !m/s/s
            !trl0(i)=fluxchan1/aream(i) !m/s/s
            aream(i) = suma
            !d_tr(i)=fluxchan2/suma
          endif !suma/=0

      enddo !i=1,np
!check for negative aicen
     ! do i = 1,nx_nh
     !  trc(i) = d_tr0(i) - trl0(i) !+ d_tr(i)
     !  sumtmp = d_tr0(i) - fluxchan2(i) / aream(i)
     !  if(sumtmp<0) then
         !write(12,*)'upwind mystery1',i,trc(i),d_tr0(i),trl0(i)
         !if(trl0(i)>0) then
         !   flux_matrix1(:,i,trc_n) = flux_matrix1(:,i,trc_n)*d_tr0(i)/trl0(i)
         !   flux_matrix2(:,i,trc_n) = flux_matrix2(:,i,trc_n)*d_tr0(i)/trl0(i)
         !   flux_matrix3(:,i,trc_n) = flux_matrix3(:,i,trc_n)*d_tr0(i)/trl0(i)
         !   flux_matrix4(:,i,trc_n) = flux_matrix4(:,i,trc_n)*d_tr0(i)/trl0(i)
         !else
         !   flux_matrix1(:,i,trc_n) = 0
         !   flux_matrix2(:,i,trc_n) = 0
         !   flux_matrix3(:,i,trc_n) = 0
         !   flux_matrix4(:,i,trc_n) = 0
         !endif
     !       do j=1,nne(i)
     !         ie=indel(j,i)
     !          id=iself(j,i)
     !          !id=indnd(j,i)
     !          id3=nxq(i34(ie)-1,id,i34(ie)) !adjacent side index
     !          id2=nxq(i34(ie)-2,id,i34(ie)) !adjacent side index
     !          
     !          nd=elnode(id3,ie)
     !          indx=0
     !          do m=1,nne(nd)
     !          if(indel(m,nd)==ie) then
     !             indx=m; exit
     !          endif
     !          enddo !m
     !          if(indx==0) call parallel_abort('STEP: failed to find (9.1)')  
     !          flux_matrix3(indx,nd,trc_n) = flux_matrix1(j,i,trc_n)*-1
     !          flux_matrix4(indx,nd,trc_n) = flux_matrix2(j,i,trc_n)*-1
      
     !          nd=elnode(id2,ie)
     !          indx=0
     !          do m=1,nne(nd)
     !          if(indel(m,nd)==ie) then
     !             indx=m; exit
     !          endif
     !          enddo !m
     !          if(indx==0) call parallel_abort('STEP: failed to find (9.1)')  
     !          flux_matrix1(indx,nd,trc_n) = flux_matrix3(j,i,trc_n)*-1
     !          flux_matrix2(indx,nd,trc_n) = flux_matrix4(j,i,trc_n)*-1
     !       enddo
     !    endif
     ! enddo

     ! do i = 1,nx_elem_nh
     !    sumtmp = 0
     !    do j=1,i34(i)
     !       nd = elnode(j,i)
     !       do m=1,nne(nd)
     !          if(indel(m,nd)==i) then
     !             indx=m; exit
     !          endif
     !       enddo
     !       sumtmp = sumtmp + flux_matrix1(indx,nd,trc_n)
     !       sumtmp = sumtmp + flux_matrix2(indx,nd,trc_n)
     !       sumtmp = sumtmp + flux_matrix3(indx,nd,trc_n)
     !       sumtmp = sumtmp + flux_matrix4(indx,nd,trc_n)
     !    enddo
     !    if(sumtmp > 1) then
     !       write(12,*) i,sumtmp
     !       do j=1,i34(i)
     !          nd = elnode(j,i)
     !          do m=1,nne(nd)
     !             if(indel(m,nd)==i) then
     !                indx=m; exit
     !             endif
     !          enddo
     !          write(12,*) indx,nd,flux_matrix1(indx,nd,trc_n),flux_matrix2(indx,nd,trc_n),flux_matrix3(indx,nd,trc_n),flux_matrix4(indx,nd,trc_n)
     !       enddo
     !    endif
     ! enddo
      
      do i = 1,nx_nh
         fluxchan1(i) = sum(flux_matrix1(:,i,trc_n)+flux_matrix2(:,i,trc_n)+flux_matrix3(:,i,trc_n)+flux_matrix4(:,i,trc_n))/2
         trl0(i) = fluxchan1(i)/aream(i)
         trc(i) = d_tr0(i) - trl0(i) !+ d_tr(i)
  
         if(trc(i)<0) then
            write(12,*)'upwind mystery2',i,trc(i),d_tr0(i),trl0(i)
         endif
      enddo
      call exchange_p2d(trc)
    end subroutine upwind_icepack_bynode

    !=======================================================================
    
    subroutine upwind_icepack3(trc,trc_n)
        
      !use mod_mesh     
       use mice_module
      implicit none
      integer (kind=int_kind), intent(in) :: trc_n
      real(kind=dbl_kind), dimension(nx), intent(inout) :: trc  
      real   (kind=dbl_kind)            :: suma,xctr2,yctr2,vnor1,vnor2,sumtmp,areatmp,&
                                           & tmpx2,tmpx3,tmpy2,tmpy3,utmp,vtmp,xtmp,ytmp,vnorm,&
                                           & trctmp, trcmin, trcmax ,puny ,tmp1 ,fluxchan0, tmp2
      integer (kind=int_kind) :: i,j,ie,id,id2,id3,isd3,isd2,     &
                                 m,n1,n2,jj,indx1,indx2,nd,nodeid
      logical (kind=log_kind) :: flag_noadv
       real(kind=dbl_kind) :: swild10(3,2),swild(3),swild2(3,2),trl0(nx),d_tr0(nx),fluxchan1(nx),fluxchan2(nx)
      !type(t_mesh),        target,        intent(in)    :: mesh
       call icepack_query_parameters(puny_out=puny)

    trl0(:) = c0
    fluxchan1(:) = c0 !positive flux
    fluxchan2(:) = c0 !negative flux
    d_tr0(:) = trc(:)
    flux_matrix1(:,:,trc_n) = c0
    flux_matrix2(:,:,trc_n) = c0
    flux_matrix3(:,:,trc_n) = c0
    flux_matrix4(:,:,trc_n) = c0
    flux_matrix0(:,:,:,trc_n) = c0
      ! Driving sequence
    do ie=1,nx_elem

          if(idry_e(ie)==1) cycle
          areatmp = area(ie)/6
          !Wet elem

          do m=1,i34(ie) !side
            n1=isidenode(1,elside(m,ie))
            n2=isidenode(2,elside(m,ie))
            if(ics==1) then
               swild10(m,1) = 0.5*(uvel(n1)+uvel(n2))
               swild10(m,2) = 0.5*(vvel(n1)+vvel(n2))
               swild(m) = (trc(n1)+trc(n2))*0.5
               swild2(m,1) = uvel(elnode(m,ie))
               swild2(m,2) = vvel(elnode(m,ie))
            else
               swild10(m,1) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),eframe(:,1,ie))+vvel(n1)*dot_product(pframe(:,2,n1),eframe(:,1,ie))&
                              &+uvel(n2)*dot_product(pframe(:,1,n2),eframe(:,1,ie))+vvel(n2)*dot_product(pframe(:,2,n2),eframe(:,1,ie)))
               swild10(m,2) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),eframe(:,2,ie))+vvel(n1)*dot_product(pframe(:,2,n1),eframe(:,2,ie))&
                              &+uvel(n2)*dot_product(pframe(:,1,n2),eframe(:,2,ie))+vvel(n2)*dot_product(pframe(:,2,n2),eframe(:,2,ie)))                      
               swild(m) = (trc(n1)+trc(n2))*0.5
               swild2(m,1) = uvel(elnode(m,ie))*dot_product(pframe(:,1,elnode(m,ie)),eframe(:,1,ie))+vvel(elnode(m,ie))*dot_product(pframe(:,2,elnode(m,ie)),eframe(:,1,ie))
               swild2(m,2) = uvel(elnode(m,ie))*dot_product(pframe(:,1,elnode(m,ie)),eframe(:,2,ie))+vvel(elnode(m,ie))*dot_product(pframe(:,2,elnode(m,ie)),eframe(:,2,ie))
            endif
            !swild10(m,1) = 0.5*(uvel(n1)*trc(n1)*dot_product(pframe(:,1,n1),sframe2(:,1,elside(m,ie)))+vvel(n1)*trc(n1)*dot_product(pframe(:,2,n1),sframe2(:,1,elside(m,ie)))&
            !                 &+uvel(n2)*trc(n2)*dot_product(pframe(:,1,n2),sframe2(:,1,elside(m,ie)))+vvel(n2)*trc(n2)*dot_product(pframe(:,2,n1),sframe2(:,1,elside(m,ie))))
            !swild10(m,2) = 0.5*(uvel(n1)*trc(n1)*dot_product(pframe(:,1,n1),sframe2(:,2,elside(m,ie)))+vvel(n1)*trc(n1)*dot_product(pframe(:,2,n1),sframe2(:,2,elside(m,ie)))&
            !                 &+uvel(n2)*trc(n2)*dot_product(pframe(:,1,n2),sframe2(:,2,elside(m,ie)))+vvel(n2)*trc(n2)*dot_product(pframe(:,2,n1),sframe2(:,2,elside(m,ie))))
            !swild10(m,1) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),sframe2(:,1,elside(m,ie)))+vvel(n1)*dot_product(pframe(:,2,n1),sframe2(:,1,elside(m,ie)))&
            !                 &+uvel(n2)*dot_product(pframe(:,1,n2),sframe2(:,1,elside(m,ie)))+vvel(n2)*dot_product(pframe(:,2,n2),sframe2(:,1,elside(m,ie))))
            !swild10(m,2) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),sframe2(:,2,elside(m,ie)))+vvel(n1)*dot_product(pframe(:,2,n1),sframe2(:,2,elside(m,ie)))&
            !                 &+uvel(n2)*dot_product(pframe(:,1,n2),sframe2(:,2,elside(m,ie)))+vvel(n2)*dot_product(pframe(:,2,n2),sframe2(:,2,elside(m,ie))))
           
        enddo
        !utmp = sum(swild10(1:i34(ie),1))/3 !vel @ centroid
        !vtmp = sum(swild10(1:i34(ie),2))/3 !vel @ centroid
        utmp=sum(swild2(1:i34(ie),1))/3 !vel @ centroid
        vtmp=sum(swild2(1:i34(ie),2))/3 !vel @ centroid
        trctmp = sum(swild(1:i34(ie)))/3

          do j =1,i34(ie) !side
            n1=isidenode(1,elside(j,ie))
            n2=isidenode(2,elside(j,ie))
            
            indx1 = 0
            do m=1,nne(n1)
               if(indel(m,n1)==ie) then
                  indx1=m; exit
               endif
            enddo
            if(indx1==0) call parallel_abort('STEP: failed to find (upwind_icepack_1)') 
            indx2 = 0
            do m=1,nne(n2)
               if(indel(m,n2)==ie) then
                  indx2=m; exit
               endif
            enddo
            if(indx2==0) call parallel_abort('STEP: failed to find (upwind_icepack_2)') 

            id  = iself(indx1,n1)
            id2 = iself(indx2,n2)
            id3 = nxq(i34(ie)-2,id,i34(ie))

            if(ics==1) then
               xctr2=xctr(ie); yctr2=yctr(ie)
               tmpx3=(xnd(n1)+xnd(n2))/2.d0
               tmpy3=(ynd(n1)+ynd(n2))/2.d0
            else!ll; use [xy]el defined in eframe
               xctr2=0.d0 !sum(xel(elnode(1:i34(ie),ie)))/i34(ie)
               yctr2=0.d0 !sum(yel(elnode(1:i34(ie),ie)))/i34(ie)
               tmpx3=(xel(id,ie)+xel(id2,ie))/2.d0
               tmpy3=(yel(id,ie)+yel(id2,ie))/2.d0
            endif!ics
            !xtmp=yctr2-tmpy3
            !ytmp=tmpx3-xctr2
            !vnor1=utmp*xtmp+vtmp*ytmp !normal vel x length 
            !vnor2=swild10(j,1)*xtmp+swild10(j,2)*ytmp !normal vel@side x length 
            if(id3.eq.id2) then
               xtmp=yctr2-tmpy3
               ytmp=tmpx3-xctr2
               !vnor1=utmp*xtmp+vtmp*ytmp !normal vel x length 
               vnor2=swild10(j,1)*xtmp+swild10(j,2)*ytmp !normal vel@side x length 
               !flux_matrix1(indx1,n1,trc_n) = dt_dyn*trctmp*vnor1
               flux_matrix2(indx1,n1,trc_n) = dt_dyn*vnor2
               !flux_matrix3(indx2,n2,trc_n) = -1*dt_dyn*trctmp*vnor1
               flux_matrix4(indx2,n2,trc_n) = -1*dt_dyn*vnor2
               !sumtmp = (dt_dyn*trctmp*vnor1+dt_dyn*swild(j)*vnor2)/2
               if(vnor2>0) then
                  flux_matrix1(indx1,n1,trc_n) = n1
                  flux_matrix2(indx1,n1,trc_n) = dt_dyn*vnor2*d_tr0(n1)
                  flux_matrix3(indx2,n2,trc_n) = n1
                  flux_matrix4(indx2,n2,trc_n) = -1*dt_dyn*vnor2*d_tr0(n1)
               else
                  flux_matrix1(indx1,n1,trc_n) = n2
                  flux_matrix2(indx1,n1,trc_n) = dt_dyn*vnor2*d_tr0(n2)
                  flux_matrix3(indx2,n2,trc_n) = n2
                  flux_matrix4(indx2,n2,trc_n) = -1*dt_dyn*vnor2*d_tr0(n2)
               endif
            else
               xtmp = tmpy3 - yctr2
               ytmp = xctr2 - tmpx3
               !vnor1=utmp*xtmp+vtmp*ytmp !normal vel x length 
               vnor2=swild10(j,1)*xtmp+swild10(j,2)*ytmp !normal vel@side x length 

               !flux_matrix3(indx1,n1,trc_n) = dt_dyn*trctmp*vnor1
               flux_matrix4(indx1,n1,trc_n) = dt_dyn*vnor2
               !flux_matrix1(indx2,n2,trc_n) = -1*dt_dyn*trctmp*vnor1
               flux_matrix2(indx2,n2,trc_n) = -1*dt_dyn*vnor2
               if(vnor2>0) then
                  flux_matrix3(indx1,n1,trc_n) = n1
                  flux_matrix4(indx1,n1,trc_n) = dt_dyn*vnor2*d_tr0(n1)
                  flux_matrix1(indx2,n2,trc_n) = n1
                  flux_matrix2(indx2,n2,trc_n) = -1*dt_dyn*vnor2*d_tr0(n1)
                  
               else
                  flux_matrix3(indx1,n1,trc_n) = n2
                  flux_matrix4(indx1,n1,trc_n) = dt_dyn*vnor2*d_tr0(n2)
                  flux_matrix1(indx2,n2,trc_n) = n2
                  flux_matrix2(indx2,n2,trc_n) = -1*dt_dyn*vnor2*d_tr0(n2)
               endif
            endif
          enddo
    enddo
   do i=1,nx_nh
      suma=0.d0 !sum of areas
      fluxchan0=0.d0
      do j=1,nne(i)
         ie=indel(j,i)
         suma=suma+area(ie)/3 !approx for quad
         fluxchan0=fluxchan0 + flux_matrix2(j,i,trc_n)
         fluxchan0=fluxchan0 + flux_matrix4(j,i,trc_n)

      enddo
         if( suma/=0 ) then
            trl0(i)=fluxchan0/suma !m/s/s
         endif !suma/=0  
          !if(abs(tmp2) > puny) write(12,*)'upwind fluxchan1',i,tmp2,trl0(i),tmp1,fluxchan0/suma,fluxchan0,fluxchan1(i),suma,area_median(i)
         trc(i) = d_tr0(i) - trl0(i) !+ d_tr(i)   
         !if(trc(i)<0) then
         !   write(12,*)'upwind mystery2',i,trc(i),d_tr0(i),trl0(i)
         !endif 
   enddo  

    call exchange_p2d(trc)
  end subroutine upwind_icepack3

  !=======================================================================
subroutine upwind_icepack_other3(trc,trc_o,trc_n,trc_d,trc_d1,trc_d2)
        
        !use mod_mesh     
     
        implicit none
        integer (kind=int_kind), intent(in) :: trc_n,trc_d
        real (kind=dbl_kind), dimension(nx), intent(inout) :: trc   ,trc_o
        real (kind=dbl_kind), dimension(nx), intent(in), optional :: trc_d1,trc_d2
        real   (kind=dbl_kind)            :: fluxchan1,fluxchan2,suma,xctr2,yctr2,vnor1,vnor2,&
                                             & tmpx2,tmpx3,tmpy2,tmpy3,utmp,vtmp,xtmp,ytmp,vnorm,&
                                             & trctmp  , trcmin, trcmax   
        integer (kind=int_kind) :: i,j,ie,id,id2,id3,isd3,isd2,     &
                                   m,n1,n2,jj
        logical (kind=log_kind) :: flag_noadv
         real(rkind) :: swild10(3,2),swild(3),trl0(nx),d_tr0(nx)
        !type(t_mesh),        target,        intent(in)    :: mesh

      trl0(:) = c0
      d_tr0(:) = trc(:)
      flux_matrix1(:,:,trc_n) = c0
      flux_matrix2(:,:,trc_n) = c0
      flux_matrix3(:,:,trc_n) = c0
      flux_matrix4(:,:,trc_n) = c0
      flux_matrix0(:,:,:,trc_n) = c0
        ! Driving sequence
      do i=1,nx_nh
        if(idry(i)==1) cycle
        !if(isbnd(1,i)/=0) cycle
        !Wet node
          fluxchan1=0.d0 !sum of fluxes;\int_\Gamma u*u_n d\Gamma
          fluxchan2=0.d0 !sum of fluxes;
          suma=0.d0 !sum of areas
          do j=1,nne(i)
            ie=indel(j,i)
            if(idry_e(ie)==1) cycle
 
            !Wet elem
            suma=suma+area(ie)/3 !approx for quad
            id=iself(j,i)
            id3=nxq(i34(ie)-1,id,i34(ie)) !adjacent side index
            id2=nxq(i34(ie)-2,id,i34(ie)) !adjacent side index
            isd3=elside(id3,ie)
            isd2=elside(id2,ie)

            n1 = flux_matrix1(j,i,trc_d) 
            flux_matrix1(j,i,trc_n) = n1

            tmpx2 = trc_o(n1)
            if (present(trc_d1)) then
               tmpx2 = tmpx2 * trc_d1(n1)
            endif
            if (present(trc_d2)) then
               tmpx2 = tmpx2 * trc_d2(n1)
            endif

            flux_matrix2(j,i,trc_n) = flux_matrix2(j,i,trc_d)*tmpx2
            fluxchan1=fluxchan1+flux_matrix2(j,i,trc_n)
            
            n2 = flux_matrix3(j,i,trc_d) 
            flux_matrix3(j,i,trc_n) = n2

            tmpx2 = trc_o(n2)
            if (present(trc_d1)) then
               tmpx2 = tmpx2 * trc_d1(n2)
            endif
            if (present(trc_d2)) then
               tmpx2 = tmpx2 * trc_d2(n2)
            endif

            flux_matrix4(j,i,trc_n) = flux_matrix4(j,i,trc_d)*tmpx2
            fluxchan1=fluxchan1 + flux_matrix4(j,i,trc_n)
            !1st segment
            !flux_matrix1(j,i,trc_n) = flux_matrix1(j,i,trc_d)*trctmp
            !flux_matrix2(j,i,trc_n) = flux_matrix2(j,i,trc_d)*swild(id3)
            !fluxchan1=fluxchan1+(flux_matrix2(j,i,trc_d)*swild(id3))
            !fluxchan1=fluxchan1+(flux_matrix1(j,i,trc_d)*trctmp+flux_matrix2(j,i,trc_d)*swild(id3))/2.d0
            !fluxchan1=fluxchan1+(flux_matrix1(j,i,trc_d)+flux_matrix2(j,i,trc_d))*(trctmp+swild(id3))/4.d0
             
            
            !2nd segment
            !Normal dir x length
!            xtmp=ycj(isd2)-yctr(ie)
!            ytmp=xctr(ie)-xcj(isd2)
            !flux_matrix3(j,i,trc_n) = flux_matrix3(j,i,trc_d)*trctmp
            !flux_matrix4(j,i,trc_n) = flux_matrix4(j,i,trc_d)*swild(id2)
            !fluxchan1=fluxchan1+(flux_matrix4(j,i,trc_d)*swild(id2))
            !fluxchan1=fluxchan1+(flux_matrix3(j,i,trc_d)*trctmp+flux_matrix4(j,i,trc_d)*swild(id2))/2.d0
            !fluxchan1=fluxchan1+(flux_matrix3(j,i,trc_d)+flux_matrix4(j,i,trc_d))*(trctmp+swild(id2))/4.d0
          enddo !j=1,nne(i)
          
          if(suma/=0) then
            trl0(i)=fluxchan1/suma !m/s/s
            !trl0(i)=fluxchan1/aream(i) !m/s/s
            !d_tr(i)=fluxchan2/suma
          endif !suma/=0
      enddo !i=1,np
      call exchange_p2d(trc)
      call exchange_p2d(d_tr0)
      call exchange_p2d(trl0)
      do i = 1,nx_nh
       trc(i) = d_tr0(i) - trl0(i) !+ d_tr(i)
       !if(trc(i)*d_tr0(i)<0) write(12,*)'upwind_other mystery',i,trc(i),d_tr0(i),trl0(i)
       !if(trc(i)*d_tr(i)<0) trc(i) = d_tr(i)
      enddo

      call exchange_p2d(trc)
    end subroutine upwind_icepack_other3

    !=======================================================================
    subroutine upwind_icepack_dep3(trc,trc_o,trc_n,trc_d)
        
        !use mod_mesh     
     
        implicit none
        integer (kind=int_kind), intent(in) :: trc_n,trc_d
        real (kind=dbl_kind), dimension(nx), intent(inout) :: trc   ,trc_o
        real   (kind=dbl_kind)            :: fluxchan1,fluxchan2,suma,xctr2,yctr2,vnor1,vnor2,&
                                             & tmpx2,tmpx3,tmpy2,tmpy3,utmp,vtmp,xtmp,ytmp,vnorm,&
                                             & trctmp, trcmin, trcmax ,sumtmp,areatmp  
        integer (kind=int_kind) :: i,j,ie,id,id2,id3,isd3,isd2,     &
                                   m,n1,n2,jj,indx,nd,nodeid
        logical (kind=log_kind) :: flag_noadv
         real(rkind) :: swild10(3,2),swild(3),trl0(nx),d_tr0(nx)
        !type(t_mesh),        target,        intent(in)    :: mesh

      trl0(:) = c0
      d_tr0(:) = trc(:)
      flux_matrix1(:,:,trc_n) = c0
      flux_matrix2(:,:,trc_n) = c0
      flux_matrix3(:,:,trc_n) = c0
      flux_matrix4(:,:,trc_n) = c0
      flux_matrix0(:,:,:,trc_n) = c0
        ! Driving sequence
      do i=1,nx_nh
        if(idry(i)==1) cycle
        !Wet node
          fluxchan1=0.d0 !sum of fluxes;\int_\Gamma u*u_n d\Gamma
          fluxchan2=0.d0 !sum of fluxes;
          suma=0.d0 !sum of areas
          do j=1,nne(i)
            ie=indel(j,i)
            if(idry_e(ie)==1) cycle
 
            !Wet elem
            areatmp = area(ie)/6
            suma=suma+area(ie)/3 !approx for quad
            id=iself(j,i)
            id3=nxq(i34(ie)-1,id,i34(ie)) !adjacent side index
            id2=nxq(i34(ie)-2,id,i34(ie)) !adjacent side index
            isd3=elside(id3,ie)
            isd2=elside(id2,ie)

            do m=1,i34(ie) !side
                n1=isidenode(1,elside(m,ie))
                n2=isidenode(2,elside(m,ie))
                tmpx2 = trc_o(n1)
                tmpy2 = trc_o(n2)
                swild(m) = (tmpx2+tmpy2)*0.5
            enddo
            trctmp=sum(swild(1:i34(ie)))/3
            !1st segment
            n1 = flux_matrix1(j,i,trc_d) 
            flux_matrix1(j,i,trc_n) = n1
            flux_matrix2(j,i,trc_n) = flux_matrix2(j,i,trc_d)*trc_o(n1)
            !sumtmp = (flux_matrix1(j,i,trc_n)+flux_matrix2(j,i,trc_n))/2.d0
            !nodeid = elnode(nxq(1,id,i34(ie)),ie)
            !if(sumtmp>0) then
            !   if(sumtmp/areatmp>d_tr0(i)*0.9)then
            !      flux_matrix1(j,i,trc_n) = flux_matrix1(j,i,trc_n)*areatmp*d_tr0(i)/sumtmp*0.9
            !      flux_matrix2(j,i,trc_n) = flux_matrix2(j,i,trc_n)*areatmp*d_tr0(i)/sumtmp*0.9
            !   endif
            !   
            !else
            !   if(sumtmp/areatmp>d_tr0(nodeid)*0.9)then
            !      flux_matrix1(j,i,trc_n) = flux_matrix1(j,i,trc_n)*areatmp*d_tr0(nodeid)/sumtmp*0.9
            !      flux_matrix2(j,i,trc_n) = flux_matrix2(j,i,trc_n)*areatmp*d_tr0(nodeid)/sumtmp*0.9
            !   endif
            !endif
            !fluxchan1=fluxchan1+(flux_matrix2(j,i,trc_d)*swild(id3))
            fluxchan1=fluxchan1+flux_matrix2(j,i,trc_n)
            !fluxchan1=fluxchan1+(flux_matrix1(j,i,trc_d)+flux_matrix2(j,i,trc_d))*(trctmp+swild(id3))/4.d0
             
            
            !2nd segment
            !Normal dir x length
!            xtmp=ycj(isd2)-yctr(ie)
!            ytmp=xctr(ie)-xcj(isd2)
            n2 = flux_matrix3(j,i,trc_d) 
            flux_matrix3(j,i,trc_n) = n2
            flux_matrix4(j,i,trc_n) = flux_matrix4(j,i,trc_d)*trc_o(n2)
            !sumtmp = (flux_matrix3(j,i,trc_n) + flux_matrix4(j,i,trc_n))/2.d0
            !nodeid = elnode(nxq(i34(ie)-1,id,i34(ie)),ie)
            !if(sumtmp>0) then
            !   if(sumtmp/areatmp>d_tr0(i)*0.9)then
            !      flux_matrix3(j,i,trc_n) = flux_matrix3(j,i,trc_n)*areatmp*d_tr0(i)/sumtmp*0.9
            !      flux_matrix4(j,i,trc_n) = flux_matrix4(j,i,trc_n)*areatmp*d_tr0(i)/sumtmp*0.9
            !   endif
            !   
            !else
            !   if(sumtmp/areatmp>d_tr0(nodeid)*0.9)then
            !      flux_matrix3(j,i,trc_n) = flux_matrix3(j,i,trc_n)*areatmp*d_tr0(nodeid)/sumtmp*0.9
            !      flux_matrix4(j,i,trc_n) = flux_matrix4(j,i,trc_n)*areatmp*d_tr0(nodeid)/sumtmp*0.9
            !   endif
            !endif
            !fluxchan1=fluxchan1+(flux_matrix4(j,i,trc_d)*swild(id2))
            fluxchan1=fluxchan1 + flux_matrix4(j,i,trc_n)
            !fluxchan1=fluxchan1+(flux_matrix3(j,i,trc_d)+flux_matrix4(j,i,trc_d))*(trctmp+swild(id2))/4.d0
          enddo !j=1,nne(i)
          
          if(suma/=0) then
            trl0(i)=fluxchan1/suma !m/s/s
            !trl0(i)=fluxchan1/aream(i) !m/s/s
            !d_tr(i)=fluxchan2/suma
          endif !suma/=0
      enddo !i=1,np
      call exchange_p2d(trc)
      call exchange_p2d(d_tr0)
      call exchange_p2d(trl0)
      do i = 1,nx_nh
       trc(i) = d_tr0(i) - trl0(i) !+ d_tr(i)

      !   if(trc(i)<0) then
      !   write(12,*)'upwind mystery vice',i,trc(i),d_tr(i),trl0(i)
      !      flux_matrix1(:,i,trc_n) = 0
      !      flux_matrix2(:,i,trc_n) = 0
      !      flux_matrix3(:,i,trc_n) = 0
      !      flux_matrix4(:,i,trc_n) = 0
      !      trc(i) = d_tr0(i)
      !      do j=1,nne(i)
      !         ie=indel(j,i)
      !         id=iself(j,i)
      !         !id=indnd(j,i)
      !         do jj=1,2 !other 2 nodes
      !            nd=elnode(nxq(jj,id,i34(ie)),ie)
      !            indx=0
      !            do m=1,nne(nd)
      !            if(indel(m,nd)==ie) then
      !               indx=m; exit
      !            endif
      !            enddo !m
      !            if(indx==0) call parallel_abort('STEP: failed to find (9.1)')  
      !            flux_matrix1(indx,nd,trc_n) = 0
      !            flux_matrix2(indx,nd,trc_n) = 0
      !            flux_matrix3(indx,nd,trc_n) = 0
      !            flux_matrix4(indx,nd,trc_n) = 0
      !         enddo
      !      enddo
      !   endif
      enddo

      call exchange_p2d(trc)
    end subroutine upwind_icepack_dep3

 !=======================================================================
    
    subroutine upwind_icepack(trc,trc_n,adv_threshold0)
        
      !use mod_mesh     
       use mice_module
      implicit none
      integer (kind=int_kind), intent(in) :: trc_n
      real(kind=dbl_kind), dimension(nx), intent(inout) :: trc,adv_threshold0  
      real   (kind=dbl_kind)            :: suma,xctr2,yctr2,vnor1,vnor2,sumtmp,areatmp,&
                                           & tmpx2,tmpx3,tmpy2,tmpy3,utmp,vtmp,xtmp,ytmp,vnorm,&
                                           & trctmp, trcmin, trcmax ,puny ,tmp1 ,fluxchan0, tmp2, adv_threshold,&
                                           & tmp3
      integer (kind=int_kind) :: i,j,ie,id,id2,id3,isd3,isd2,     &
                                 m,n1,n2,jj,indx1,indx2,nd,nodeid
      logical (kind=log_kind) :: flag_noadv
       real(kind=dbl_kind) :: swild10(3,2),swild(3),swild2(3,2),trl0(nx),d_tr0(nx),fluxchan1(nx),fluxchan2(nx)
      !type(t_mesh),        target,        intent(in)    :: mesh
       call icepack_query_parameters(puny_out=puny)

    adv_threshold = 10
    trl0(:) = c0
    fluxchan1(:) = c0 !positive flux
    fluxchan2(:) = c0 !negative flux
    d_tr0(:) = trc(:)
    flux_matrix1(:,:,trc_n) = c0
    flux_matrix2(:,:,trc_n) = c0
    flux_matrix3(:,:,trc_n) = c0
    flux_matrix4(:,:,trc_n) = c0
    flux_matrix00(:,:,:) = c0
    flux_matrix0(:,:,:,trc_n) = c0
      ! Driving sequence
    do ie=1,nx_elem

          if(idry_e(ie)==1) cycle
          areatmp = area(ie)/6
          !Wet elem

          do m=1,i34(ie) !side
            n1=isidenode(1,elside(m,ie))
            n2=isidenode(2,elside(m,ie))
            if(ics==1) then 
               swild10(m,1) = 0.5*(uvel(n1)+uvel(n2))
               swild10(m,2) = 0.5*(vvel(n1)+vvel(n2))
               swild(m) = (trc(n1)+trc(n2))*0.5
               swild2(m,1) = uvel(elnode(m,ie))
               swild2(m,2) = vvel(elnode(m,ie))
            else
               swild10(m,1) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),eframe(:,1,ie))+vvel(n1)*dot_product(pframe(:,2,n1),eframe(:,1,ie))&
                              &+uvel(n2)*dot_product(pframe(:,1,n2),eframe(:,1,ie))+vvel(n2)*dot_product(pframe(:,2,n2),eframe(:,1,ie)))
               swild10(m,2) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),eframe(:,2,ie))+vvel(n1)*dot_product(pframe(:,2,n1),eframe(:,2,ie))&
                              &+uvel(n2)*dot_product(pframe(:,1,n2),eframe(:,2,ie))+vvel(n2)*dot_product(pframe(:,2,n2),eframe(:,2,ie)))                      
               swild(m) = (trc(n1)+trc(n2))*0.5
               swild2(m,1) = uvel(elnode(m,ie))*dot_product(pframe(:,1,elnode(m,ie)),eframe(:,1,ie))+vvel(elnode(m,ie))*dot_product(pframe(:,2,elnode(m,ie)),eframe(:,1,ie))
               swild2(m,2) = uvel(elnode(m,ie))*dot_product(pframe(:,1,elnode(m,ie)),eframe(:,2,ie))+vvel(elnode(m,ie))*dot_product(pframe(:,2,elnode(m,ie)),eframe(:,2,ie))

            endif
            !swild10(m,1) = 0.5*(uvel(n1)*trc(n1)*dot_product(pframe(:,1,n1),sframe2(:,1,elside(m,ie)))+vvel(n1)*trc(n1)*dot_product(pframe(:,2,n1),sframe2(:,1,elside(m,ie)))&
            !                 &+uvel(n2)*trc(n2)*dot_product(pframe(:,1,n2),sframe2(:,1,elside(m,ie)))+vvel(n2)*trc(n2)*dot_product(pframe(:,2,n1),sframe2(:,1,elside(m,ie))))
            !swild10(m,2) = 0.5*(uvel(n1)*trc(n1)*dot_product(pframe(:,1,n1),sframe2(:,2,elside(m,ie)))+vvel(n1)*trc(n1)*dot_product(pframe(:,2,n1),sframe2(:,2,elside(m,ie)))&
            !                 &+uvel(n2)*trc(n2)*dot_product(pframe(:,1,n2),sframe2(:,2,elside(m,ie)))+vvel(n2)*trc(n2)*dot_product(pframe(:,2,n1),sframe2(:,2,elside(m,ie))))
            !swild10(m,1) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),sframe2(:,1,elside(m,ie)))+vvel(n1)*dot_product(pframe(:,2,n1),sframe2(:,1,elside(m,ie)))&
            !                 &+uvel(n2)*dot_product(pframe(:,1,n2),sframe2(:,1,elside(m,ie)))+vvel(n2)*dot_product(pframe(:,2,n2),sframe2(:,1,elside(m,ie))))
            !swild10(m,2) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),sframe2(:,2,elside(m,ie)))+vvel(n1)*dot_product(pframe(:,2,n1),sframe2(:,2,elside(m,ie)))&
            !                 &+uvel(n2)*dot_product(pframe(:,1,n2),sframe2(:,2,elside(m,ie)))+vvel(n2)*dot_product(pframe(:,2,n2),sframe2(:,2,elside(m,ie))))
          enddo
        !utmp = sum(swild10(1:i34(ie),1))/3 !vel @ centroid
        !vtmp = sum(swild10(1:i34(ie),2))/3 !vel @ centroid
        utmp=sum(swild2(1:i34(ie),1))/3 !vel @ centroid
        vtmp=sum(swild2(1:i34(ie),2))/3 !vel @ centroid
        trctmp = sum(swild(1:i34(ie)))/3

          do j =1,i34(ie) !side
            n1=isidenode(1,elside(j,ie))
            n2=isidenode(2,elside(j,ie))
            
            indx1 = 0
            do m=1,nne(n1)
               if(indel(m,n1)==ie) then
                  indx1=m; exit
               endif
            enddo
            if(indx1==0) call parallel_abort('STEP: failed to find (upwind_icepack_1)') 
            indx2 = 0
            do m=1,nne(n2)
               if(indel(m,n2)==ie) then
                  indx2=m; exit
               endif
            enddo
            if(indx2==0) call parallel_abort('STEP: failed to find (upwind_icepack_2)') 

            id  = iself(indx1,n1)
            id2 = iself(indx2,n2)
            id3 = nxq(i34(ie)-2,id,i34(ie))

            if(ics==1) then
               xctr2=xctr(ie); yctr2=yctr(ie)
               tmpx3=(xnd(n1)+xnd(n2))/2.d0
               tmpy3=(ynd(n1)+ynd(n2))/2.d0
            else!ll; use [xy]el defined in eframe
               xctr2=0.d0 !sum(xel(elnode(1:i34(ie),ie)))/i34(ie)
               yctr2=0.d0 !sum(yel(elnode(1:i34(ie),ie)))/i34(ie)
               tmpx3=(xel(id,ie)+xel(id2,ie))/2.d0
               tmpy3=(yel(id,ie)+yel(id2,ie))/2.d0
            endif!ics
            !xtmp=yctr2-tmpy3
            !ytmp=tmpx3-xctr2
            !vnor1=utmp*xtmp+vtmp*ytmp !normal vel x length 
            !vnor2=swild10(j,1)*xtmp+swild10(j,2)*ytmp !normal vel@side x length 
            if(id3.eq.id2) then
               xtmp=yctr2-tmpy3
               ytmp=tmpx3-xctr2
               vnor1=utmp*xtmp+vtmp*ytmp !normal vel x length 
               vnor2=swild10(j,1)*xtmp+swild10(j,2)*ytmp !normal vel@side x length 

               sumtmp = (dt_dyn*vnor1+dt_dyn*vnor2)/2
               if(sumtmp>0) then
                  tmp3 = d_tr0(n2)-d_tr0(n1)
                  tmp2=max(c0,tmp3)
                  sumtmp = sumtmp * d_tr0(n1)
                  flux_matrix00(indx1,n1,1) = n1
                  flux_matrix00(indx2,n2,2) = n1
                  flux_matrix1(indx1,n1,trc_n) = dt_dyn * vnor1 * d_tr0(n1)
                  flux_matrix2(indx1,n1,trc_n) = dt_dyn * vnor2 * d_tr0(n1)
                  flux_matrix3(indx2,n2,trc_n) = -1 * dt_dyn * vnor1 * d_tr0(n1)
                  flux_matrix4(indx2,n2,trc_n) = -1 * dt_dyn * vnor2 * d_tr0(n1)              
                  !adv_threshold = exp(-20*tmp2)
                  !if((sumtmp/areatmp)>(d_tr0(n1)*adv_threshold))then
                  !if((sumtmp/areatmp)>(d_tr0(n1)*adv_threshold0(n1)).or.(sumtmp/areatmp)>(0.01))then
                  if((sumtmp/areatmp)>(d_tr0(n1)*adv_threshold).or.(sumtmp/areatmp)>(d_tr0(n1)*adv_threshold0(n2)))then
                     !tmp1 = areatmp*d_tr0(n1)/sumtmp*adv_threshold
                     tmp2 = min(d_tr0(n1)*adv_threshold,d_tr0(n1)*adv_threshold0(n2))
                     tmp1 = areatmp*tmp2/sumtmp
                     flux_matrix1(indx1,n1,trc_n) = flux_matrix1(indx1,n1,trc_n)*tmp1
                     flux_matrix2(indx1,n1,trc_n) = flux_matrix2(indx1,n1,trc_n)*tmp1
                     flux_matrix3(indx2,n2,trc_n) = flux_matrix3(indx2,n2,trc_n)*tmp1
                     flux_matrix4(indx2,n2,trc_n) = flux_matrix4(indx2,n2,trc_n)*tmp1
                  endif
                  
               else
                  tmp3 = d_tr0(n1)-d_tr0(n2)
                  tmp2=max(c0,tmp3)
                  sumtmp = sumtmp * d_tr0(n2)
                  flux_matrix00(indx1,n1,1) = n2
                  flux_matrix00(indx2,n2,2) = n2
                  flux_matrix1(indx1,n1,trc_n) = dt_dyn * vnor1 * d_tr0(n2)
                  flux_matrix2(indx1,n1,trc_n) = dt_dyn * vnor2 * d_tr0(n2)
                  flux_matrix3(indx2,n2,trc_n) = -1 * dt_dyn * vnor1 * d_tr0(n2)
                  flux_matrix4(indx2,n2,trc_n) = -1 * dt_dyn * vnor2 * d_tr0(n2)                    
                  !adv_threshold = exp(-20*tmp2)
                  !if((abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold))then
                  !if((abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold0(n2)).or.(abs(sumtmp)/areatmp)>(0.01))then
                  if((abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold).or.(abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold0(n1)))then
                     !tmp1 = areatmp*d_tr0(n2)/abs(sumtmp)*adv_threshold
                     tmp2 = min(d_tr0(n2)*adv_threshold,d_tr0(n2)*adv_threshold0(n1))
                     tmp1 = areatmp*tmp2/abs(sumtmp)
                     flux_matrix1(indx1,n1,trc_n) = flux_matrix1(indx1,n1,trc_n)*tmp1
                     flux_matrix2(indx1,n1,trc_n) = flux_matrix2(indx1,n1,trc_n)*tmp1
                     flux_matrix3(indx2,n2,trc_n) = flux_matrix3(indx2,n2,trc_n)*tmp1
                     flux_matrix4(indx2,n2,trc_n) = flux_matrix4(indx2,n2,trc_n)*tmp1
                  endif
               endif
            else
               xtmp = tmpy3 - yctr2
               ytmp = xctr2 - tmpx3
               vnor1=utmp*xtmp+vtmp*ytmp !normal vel x length 
               vnor2=swild10(j,1)*xtmp+swild10(j,2)*ytmp !normal vel@side x length 

               sumtmp = (dt_dyn*vnor1+dt_dyn*vnor2)/2
               if(sumtmp>0) then
                  tmp3 = d_tr0(n2)-d_tr0(n1)
                  tmp2=max(c0,tmp3)
                  sumtmp = sumtmp * d_tr0(n1)
                  flux_matrix00(indx1,n1,2) = n1
                  flux_matrix00(indx2,n2,1) = n1
                  flux_matrix3(indx1,n1,trc_n) = dt_dyn*vnor1 * d_tr0(n1)
                  flux_matrix4(indx1,n1,trc_n) = dt_dyn*vnor2 * d_tr0(n1)
                  flux_matrix1(indx2,n2,trc_n) = -1*dt_dyn*vnor1 * d_tr0(n1)
                  flux_matrix2(indx2,n2,trc_n) = -1*dt_dyn*vnor2 * d_tr0(n1)
                  !adv_threshold = exp(-20*tmp2)
                  !if((sumtmp/areatmp)>(d_tr0(n1)*adv_threshold))then
                  !if((sumtmp/areatmp)>(d_tr0(n1)*adv_threshold0(n1)).or.(sumtmp/areatmp)>(0.01))then
                  if((sumtmp/areatmp)>(d_tr0(n1)*adv_threshold).or.(sumtmp/areatmp)>(d_tr0(n1)*adv_threshold0(n2)))then
                     !tmp1 = areatmp*d_tr0(n1)/sumtmp*adv_threshold
                     tmp2 = min(d_tr0(n1)*adv_threshold,d_tr0(n1)*adv_threshold0(n2))
                     tmp1 = areatmp*tmp2/sumtmp
                     flux_matrix3(indx1,n1,trc_n) = flux_matrix3(indx1,n1,trc_n)*tmp1
                     flux_matrix4(indx1,n1,trc_n) = flux_matrix4(indx1,n1,trc_n)*tmp1
                     flux_matrix1(indx2,n2,trc_n) = flux_matrix1(indx2,n2,trc_n)*tmp1
                     flux_matrix2(indx2,n2,trc_n) = flux_matrix2(indx2,n2,trc_n)*tmp1
                  endif
                  
               else
                  tmp3 = d_tr0(n1)-d_tr0(n2)
                  tmp2=max(c0,tmp3)
                  sumtmp = sumtmp * d_tr0(n2)
                  flux_matrix00(indx1,n1,2) = n2
                  flux_matrix00(indx2,n2,1) = n2
                  flux_matrix3(indx1,n1,trc_n) = dt_dyn * vnor1 * d_tr0(n2)
                  flux_matrix4(indx1,n1,trc_n) = dt_dyn * vnor2 * d_tr0(n2)
                  flux_matrix1(indx2,n2,trc_n) = -1 * dt_dyn * vnor1 * d_tr0(n2)
                  flux_matrix2(indx2,n2,trc_n) = -1 * dt_dyn * vnor2 * d_tr0(n2)
                  !adv_threshold = exp(-20*tmp2)
                  !if((abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold))then
                  !if((abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold0(n2)).or.(abs(sumtmp)/areatmp)>0.01)then
                  if((abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold).or.(abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold0(n1)))then
                     !tmp1 = areatmp*d_tr0(n2)/abs(sumtmp)*adv_threshold
                     tmp2 = min(d_tr0(n2)*adv_threshold,d_tr0(n2)*adv_threshold0(n1))
                     tmp1 = areatmp*tmp2/abs(sumtmp)
                     flux_matrix3(indx1,n1,trc_n) = flux_matrix3(indx1,n1,trc_n)*tmp1
                     flux_matrix4(indx1,n1,trc_n) = flux_matrix4(indx1,n1,trc_n)*tmp1
                     flux_matrix1(indx2,n2,trc_n) = flux_matrix1(indx2,n2,trc_n)*tmp1
                     flux_matrix2(indx2,n2,trc_n) = flux_matrix2(indx2,n2,trc_n)*tmp1
                  endif
               endif
            endif
          enddo
    enddo
   do i=1,nx_nh
      suma=0.d0 !sum of areas
      fluxchan0=0.d0
      do j=1,nne(i)
         ie=indel(j,i)
         suma=suma+area(ie)/3 !approx for quad
         fluxchan0=fluxchan0 + (flux_matrix1(j,i,trc_n) + flux_matrix2(j,i,trc_n))/2.d0
         fluxchan0=fluxchan0 + (flux_matrix3(j,i,trc_n) + flux_matrix4(j,i,trc_n))/2.d0

      enddo
         if( suma/=0 ) then
            trl0(i)=fluxchan0/suma !m/s/s
         endif !suma/=0  
         trc(i) = d_tr0(i) - trl0(i) !+ d_tr(i)   
   enddo  
    call exchange_p2d(trc)
  end subroutine upwind_icepack

  !=======================================================================
subroutine upwind_icepack_other(trc,trc_o,trc_n,trc_d,trc_d1,trc_d2)
        
        !use mod_mesh     
     
        implicit none
        integer (kind=int_kind), intent(in) :: trc_n,trc_d
        real (kind=dbl_kind), dimension(nx), intent(inout) :: trc   ,trc_o
        real (kind=dbl_kind), dimension(nx), intent(in), optional :: trc_d1,trc_d2
        real   (kind=dbl_kind)            :: fluxchan1,fluxchan2,suma,xctr2,yctr2,vnor1,vnor2,&
                                             & tmpx2,tmpx3,tmpy2,tmpy3,utmp,vtmp,xtmp,ytmp,vnorm,&
                                             & trctmp  , trcmin, trcmax, tmpnode   
        integer (kind=int_kind) :: i,j,ie,id,id2,id3,isd3,isd2,     &
                                   m,n1,n2,jj
        logical (kind=log_kind) :: flag_noadv
         real(rkind) :: swild10(3,2),swild(3),trl0(nx),d_tr0(nx)
        !type(t_mesh),        target,        intent(in)    :: mesh

      trl0(:) = c0
      d_tr0(:) = trc(:)
      flux_matrix1(:,:,trc_n) = c0
      flux_matrix2(:,:,trc_n) = c0
      flux_matrix3(:,:,trc_n) = c0
      flux_matrix4(:,:,trc_n) = c0
      flux_matrix0(:,:,:,trc_n) = c0
        ! Driving sequence
      do i=1,nx_nh
        if(idry(i)==1) cycle
        !if(isbnd(1,i)/=0) cycle
        !Wet node
          fluxchan1=0.d0 !sum of fluxes;\int_\Gamma u*u_n d\Gamma
          fluxchan2=0.d0 !sum of fluxes;
          suma=0.d0 !sum of areas
          do j=1,nne(i)
            ie=indel(j,i)
            if(idry_e(ie)==1) cycle
 
            !Wet elem
            suma=suma+area(ie)/3 !approx for quad
            id=iself(j,i)
            id3=nxq(i34(ie)-1,id,i34(ie)) !adjacent side index
            id2=nxq(i34(ie)-2,id,i34(ie)) !adjacent side index
            isd3=elside(id3,ie)
            isd2=elside(id2,ie)

            
            do m=1,i34(ie) !side
                n1=isidenode(1,elside(m,ie))
                n2=isidenode(2,elside(m,ie))
                tmpx2 = trc_o(n1)
                tmpy2 = trc_o(n2)
                if (present(trc_d1)) then
                  tmpx2 = tmpx2 * trc_d1(n1)
                  tmpy2 = tmpy2 * trc_d1(n2)
                endif
                if (present(trc_d2)) then
                  tmpx2 = tmpx2 * trc_d2(n1)
                  tmpy2 = tmpy2 * trc_d2(n2)
                endif
                swild(m) = (tmpx2+tmpy2)*0.5
            enddo
            trctmp = sum(swild(1:i34(ie)))/3
            trcmin = minval(trc_o(indnd(1:nnp(i),i)))
            trcmin = min(trcmin,trc_o(i))
            trcmax = maxval(trc_o(indnd(1:nnp(i),i)))
            trcmax = max(trcmax,trc_o(i))
            
            tmpnode = flux_matrix00(j,i,1)
            tmpx2 = trc_o(tmpnode)
            if (present(trc_d1)) then
               tmpx2 = tmpx2 * trc_d1(tmpnode)
            endif
            if (present(trc_d2)) then
               tmpx2 = tmpx2 * trc_d2(tmpnode)
            endif
            !1st segment
            flux_matrix1(j,i,trc_n) = flux_matrix1(j,i,trc_d)*tmpx2
            flux_matrix2(j,i,trc_n) = flux_matrix2(j,i,trc_d)*tmpx2
            fluxchan1=fluxchan1+(flux_matrix1(j,i,trc_d)*tmpx2+flux_matrix2(j,i,trc_d)*tmpx2)/2.d0
             
            
            !2nd segment
            !Normal dir x length
!            xtmp=ycj(isd2)-yctr(ie)
!            ytmp=xctr(ie)-xcj(isd2)
            tmpnode = flux_matrix00(j,i,2)
            tmpx2 = trc_o(tmpnode)
            if (present(trc_d1)) then
               tmpx2 = tmpx2 * trc_d1(tmpnode)
            endif
            if (present(trc_d2)) then
               tmpx2 = tmpx2 * trc_d2(tmpnode)
            endif
            flux_matrix3(j,i,trc_n) = flux_matrix3(j,i,trc_d)*tmpx2
            flux_matrix4(j,i,trc_n) = flux_matrix4(j,i,trc_d)*tmpx2
            fluxchan1=fluxchan1+(flux_matrix3(j,i,trc_d)*tmpx2+flux_matrix4(j,i,trc_d)*tmpx2)/2.d0
          enddo !j=1,nne(i)
          
          if(suma/=0) then
            trl0(i)=fluxchan1/suma !m/s/s
          endif !suma/=0
      enddo !i=1,np
      call exchange_p2d(trc)
      call exchange_p2d(d_tr0)
      call exchange_p2d(trl0)
      do i = 1,nx_nh
       trc(i) = d_tr0(i) - trl0(i) !+ d_tr(i)
      enddo

      call exchange_p2d(trc)
    end subroutine upwind_icepack_other

    !=======================================================================
    subroutine upwind_icepack_dep(trc,trc_o,trc_n,trc_d)
        
        !use mod_mesh     
     
        implicit none
        integer (kind=int_kind), intent(in) :: trc_n,trc_d
        real (kind=dbl_kind), dimension(nx), intent(inout) :: trc   ,trc_o
        real   (kind=dbl_kind)            :: fluxchan1,fluxchan2,suma,xctr2,yctr2,vnor1,vnor2,&
                                             & tmpx2,tmpx3,tmpy2,tmpy3,utmp,vtmp,xtmp,ytmp,vnorm,&
                                             & trctmp, trcmin, trcmax ,sumtmp,areatmp, tmpnode  
        integer (kind=int_kind) :: i,j,ie,id,id2,id3,isd3,isd2,     &
                                   m,n1,n2,jj,indx,nd,nodeid
        logical (kind=log_kind) :: flag_noadv
         real(rkind) :: swild10(3,2),swild(3),trl0(nx),d_tr0(nx)
        !type(t_mesh),        target,        intent(in)    :: mesh

      trl0(:) = c0
      d_tr0(:) = trc(:)
      flux_matrix1(:,:,trc_n) = c0
      flux_matrix2(:,:,trc_n) = c0
      flux_matrix3(:,:,trc_n) = c0
      flux_matrix4(:,:,trc_n) = c0
      flux_matrix0(:,:,:,trc_n) = c0
        ! Driving sequence
      do i=1,nx_nh
        if(idry(i)==1) cycle
        !Wet node
          fluxchan1=0.d0 !sum of fluxes;\int_\Gamma u*u_n d\Gamma
          fluxchan2=0.d0 !sum of fluxes;
          suma=0.d0 !sum of areas
          do j=1,nne(i)
            ie=indel(j,i)
            if(idry_e(ie)==1) cycle
 
            !Wet elem
            areatmp = area(ie)/6
            suma=suma+area(ie)/3 !approx for quad
            id=iself(j,i)
            id3=nxq(i34(ie)-1,id,i34(ie)) !adjacent side index
            id2=nxq(i34(ie)-2,id,i34(ie)) !adjacent side index
            isd3=elside(id3,ie)
            isd2=elside(id2,ie)

            do m=1,i34(ie) !side
                n1=isidenode(1,elside(m,ie))
                n2=isidenode(2,elside(m,ie))
                tmpx2 = trc_o(n1)
                tmpy2 = trc_o(n2)
                swild(m) = (tmpx2+tmpy2)*0.5
            enddo
            trctmp=sum(swild(1:i34(ie)))/3
            !1st segment
            tmpnode = flux_matrix00(j,i,1)
            tmpx2 = trc_o(tmpnode)
            flux_matrix1(j,i,trc_n) = flux_matrix1(j,i,trc_d)*tmpx2
            flux_matrix2(j,i,trc_n) = flux_matrix2(j,i,trc_d)*tmpx2
            sumtmp = (flux_matrix1(j,i,trc_n)+flux_matrix2(j,i,trc_n))/2.d0
            nodeid = elnode(nxq(1,id,i34(ie)),ie)
            fluxchan1=fluxchan1+(flux_matrix1(j,i,trc_n)+flux_matrix2(j,i,trc_n))/2.d0
             
            
            !2nd segment
            !Normal dir x length
!            xtmp=ycj(isd2)-yctr(ie)
!            ytmp=xctr(ie)-xcj(isd2)
            tmpnode = flux_matrix00(j,i,2)
            tmpx2 = trc_o(tmpnode)
            flux_matrix3(j,i,trc_n) = flux_matrix3(j,i,trc_d)*tmpx2
            flux_matrix4(j,i,trc_n) = flux_matrix4(j,i,trc_d)*tmpx2
            sumtmp = (flux_matrix3(j,i,trc_n) + flux_matrix4(j,i,trc_n))/2.d0
            nodeid = elnode(nxq(i34(ie)-1,id,i34(ie)),ie)

            fluxchan1=fluxchan1+(flux_matrix3(j,i,trc_n) + flux_matrix4(j,i,trc_n))/2.d0
          enddo !j=1,nne(i)
          
          if(suma/=0) then
            trl0(i)=fluxchan1/suma !m/s/s
          endif !suma/=0
      enddo !i=1,np
      call exchange_p2d(trc)
      call exchange_p2d(d_tr0)
      call exchange_p2d(trl0)
      do i = 1,nx_nh
       trc(i) = d_tr0(i) - trl0(i) !+ d_tr(i)


      enddo

      call exchange_p2d(trc)
    end subroutine upwind_icepack_dep

 !=======================================================================
    
    subroutine tvd_icepack(trc,trc_n,adv_threshold0)
        
      !use mod_mesh     
       use mice_module
       use schism_glbl, only: dldxy
      implicit none
      integer (kind=int_kind), intent(in) :: trc_n
      real(kind=dbl_kind), dimension(nx), intent(inout) :: trc,adv_threshold0  
      real   (kind=dbl_kind)            :: suma,xctr2,yctr2,vnor1,vnor2,sumtmp,areatmp,&
                                           & tmpx2,tmpx3,tmpy2,tmpy3,utmp,vtmp,xtmp,ytmp,vnorm,&
                                           & trctmp, trcmin, trcmax ,puny ,tmp1 ,fluxchan0, tmp2, adv_threshold,&
                                           & tmp3, psi_r, sum1, sum2
      integer (kind=int_kind) :: i,j,ie,id,id2,id3,isd3,isd2,     &
                                 m,n1,n2,jj,indx1,indx2,nd,nodeid
      logical (kind=log_kind) :: flag_noadv
       real(kind=dbl_kind) :: swild10(3,2),swild(3),swild2(3,2),trl0(nx),d_tr0(nx),fluxchan1(nx),fluxchan2(nx), &
                              & r_fac(nx),grad_x(nx),grad_y(nx)
      !type(t_mesh),        target,        intent(in)    :: mesh
       call icepack_query_parameters(puny_out=puny)

    adv_threshold = 10
    trl0(:) = c0
    fluxchan1(:) = c0 !positive flux
    fluxchan2(:) = c0 !negative flux
    d_tr0(:) = trc(:)
    flux_matrix1(:,:,trc_n) = c0
    flux_matrix2(:,:,trc_n) = c0
    flux_matrix3(:,:,trc_n) = c0
    flux_matrix4(:,:,trc_n) = c0
    flux_matrix00(:,:,:) = c0
    flux_matrix0(:,:,:,trc_n) = c0
    r_fac(:) = c0
    grad_x(:) = c0
    grad_y(:) = c0

    do i = 1,nx_nh
    sum1=0; sum2=0
      do j=1,nne(i)
         ie=indel(j,i)
         id=iself(j,i)
         if(ics==1) then
            sum1=sum1+area(ie)/3*dot_product(trc(elnode(1:i34(ie),ie)),dldxy(1:i34(ie),1,ie))
            sum2=sum2+area(ie)/3*dot_product(trc(elnode(1:i34(ie),ie)),dldxy(1:i34(ie),2,ie))
         else
          sum1=sum1+area(ie)/3*dot_product(trc(elnode(1:i34(ie),ie)),dldxy(1:i34(ie),1,ie))*dot_product(eframe(1:3,1,ie),pframe(1:3,1,i))+&
                   &area(ie)/3*dot_product(trc(elnode(1:i34(ie),ie)),dldxy(1:i34(ie),2,ie))*dot_product(eframe(1:3,2,ie),pframe(1:3,1,i))
          sum2=sum2+area(ie)/3*dot_product(trc(elnode(1:i34(ie),ie)),dldxy(1:i34(ie),1,ie))*dot_product(eframe(1:3,1,ie),pframe(1:3,2,i))+&
                   &area(ie)/3*dot_product(trc(elnode(1:i34(ie),ie)),dldxy(1:i34(ie),2,ie))*dot_product(eframe(1:3,2,ie),pframe(1:3,2,i))
         endif
      enddo
      grad_x(i) = sum1/area_median(i)
      grad_y(i) = sum2/area_median(i)
    enddo
    call exchange_p2d(grad_x)
    call exchange_p2d(grad_y)
      ! Driving sequence
    do ie=1,nx_elem

          if(idry_e(ie)==1) cycle
          areatmp = area(ie)/6
          !Wet elem

          do m=1,i34(ie) !side
            n1=isidenode(1,elside(m,ie))
            n2=isidenode(2,elside(m,ie))
            if(ics==1) then
               swild10(m,1) = 0.5*(uvel(n1)+uvel(n2))
               swild10(m,2) = 0.5*(vvel(n1)+vvel(n2))
               swild(m) = (trc(n1)+trc(n2))*0.5
               swild2(m,1) = uvel(elnode(m,ie))
               swild2(m,2) = vvel(elnode(m,ie))
            else
               swild10(m,1) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),eframe(:,1,ie))+vvel(n1)*dot_product(pframe(:,2,n1),eframe(:,1,ie))&
                                 &+uvel(n2)*dot_product(pframe(:,1,n2),eframe(:,1,ie))+vvel(n2)*dot_product(pframe(:,2,n2),eframe(:,1,ie)))
               swild10(m,2) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),eframe(:,2,ie))+vvel(n1)*dot_product(pframe(:,2,n1),eframe(:,2,ie))&
                                 &+uvel(n2)*dot_product(pframe(:,1,n2),eframe(:,2,ie))+vvel(n2)*dot_product(pframe(:,2,n2),eframe(:,2,ie)))                      
               swild(m) = (trc(n1)+trc(n2))*0.5
               swild2(m,1) = uvel(elnode(m,ie))*dot_product(pframe(:,1,elnode(m,ie)),eframe(:,1,ie))+vvel(elnode(m,ie))*dot_product(pframe(:,2,elnode(m,ie)),eframe(:,1,ie))
               swild2(m,2) = uvel(elnode(m,ie))*dot_product(pframe(:,1,elnode(m,ie)),eframe(:,2,ie))+vvel(elnode(m,ie))*dot_product(pframe(:,2,elnode(m,ie)),eframe(:,2,ie))
            endif
            !swild10(m,1) = 0.5*(uvel(n1)*trc(n1)*dot_product(pframe(:,1,n1),sframe2(:,1,elside(m,ie)))+vvel(n1)*trc(n1)*dot_product(pframe(:,2,n1),sframe2(:,1,elside(m,ie)))&
            !                 &+uvel(n2)*trc(n2)*dot_product(pframe(:,1,n2),sframe2(:,1,elside(m,ie)))+vvel(n2)*trc(n2)*dot_product(pframe(:,2,n1),sframe2(:,1,elside(m,ie))))
            !swild10(m,2) = 0.5*(uvel(n1)*trc(n1)*dot_product(pframe(:,1,n1),sframe2(:,2,elside(m,ie)))+vvel(n1)*trc(n1)*dot_product(pframe(:,2,n1),sframe2(:,2,elside(m,ie)))&
            !                 &+uvel(n2)*trc(n2)*dot_product(pframe(:,1,n2),sframe2(:,2,elside(m,ie)))+vvel(n2)*trc(n2)*dot_product(pframe(:,2,n1),sframe2(:,2,elside(m,ie))))
            !swild10(m,1) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),sframe2(:,1,elside(m,ie)))+vvel(n1)*dot_product(pframe(:,2,n1),sframe2(:,1,elside(m,ie)))&
            !                 &+uvel(n2)*dot_product(pframe(:,1,n2),sframe2(:,1,elside(m,ie)))+vvel(n2)*dot_product(pframe(:,2,n2),sframe2(:,1,elside(m,ie))))
            !swild10(m,2) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),sframe2(:,2,elside(m,ie)))+vvel(n1)*dot_product(pframe(:,2,n1),sframe2(:,2,elside(m,ie)))&
            !                 &+uvel(n2)*dot_product(pframe(:,1,n2),sframe2(:,2,elside(m,ie)))+vvel(n2)*dot_product(pframe(:,2,n2),sframe2(:,2,elside(m,ie))))
            
         enddo
        !utmp = sum(swild10(1:i34(ie),1))/3 !vel @ centroid
        !vtmp = sum(swild10(1:i34(ie),2))/3 !vel @ centroid
        utmp=sum(swild2(1:i34(ie),1))/3 !vel @ centroid
        vtmp=sum(swild2(1:i34(ie),2))/3 !vel @ centroid
        trctmp = sum(swild(1:i34(ie)))/3

          do j =1,i34(ie) !side
            n1=isidenode(1,elside(j,ie))
            n2=isidenode(2,elside(j,ie))
            
            indx1 = 0
            do m=1,nne(n1)
               if(indel(m,n1)==ie) then
                  indx1=m; exit
               endif
            enddo
            if(indx1==0) call parallel_abort('STEP: failed to find (tvd_icepack_1)') 
            indx2 = 0
            do m=1,nne(n2)
               if(indel(m,n2)==ie) then
                  indx2=m; exit
               endif
            enddo
            if(indx2==0) call parallel_abort('STEP: failed to find (tvd_icepack_2)') 

            id  = iself(indx1,n1)
            id2 = iself(indx2,n2)
            id3 = nxq(i34(ie)-2,id,i34(ie))

            if(ics==1) then
               xctr2=xctr(ie); yctr2=yctr(ie)
               tmpx3=(xnd(n1)+xnd(n2))/2.d0
               tmpy3=(ynd(n1)+ynd(n2))/2.d0
            else!ll; use [xy]el defined in eframe
               xctr2=0.d0 !sum(xel(elnode(1:i34(ie),ie)))/i34(ie)
               yctr2=0.d0 !sum(yel(elnode(1:i34(ie),ie)))/i34(ie)
               tmpx3=(xel(id,ie)+xel(id2,ie))/2.d0
               tmpy3=(yel(id,ie)+yel(id2,ie))/2.d0
            endif!ics
            !xtmp=yctr2-tmpy3
            !ytmp=tmpx3-xctr2
            !vnor1=utmp*xtmp+vtmp*ytmp !normal vel x length 
            !vnor2=swild10(j,1)*xtmp+swild10(j,2)*ytmp !normal vel@side x length 
            if(id3.eq.id2) then
               xtmp=yctr2-tmpy3
               ytmp=tmpx3-xctr2
               vnor1=utmp*xtmp+vtmp*ytmp !normal vel x length 
               vnor2=swild10(j,1)*xtmp+swild10(j,2)*ytmp !normal vel@side x length 
               !vnor1=vnor2

               sumtmp = (dt_dyn*vnor1+dt_dyn*vnor2)/2
               if(sumtmp>0) then
                  tmp3 = d_tr0(n2)-d_tr0(n1)
                  tmp2=max(c0,tmp3)
                  !sumtmp = sumtmp * d_tr0(n1)
                  flux_matrix00(indx1,n1,1) = n1
                  flux_matrix00(indx2,n2,2) = n1
                  ! n1 center, n2 dowmstream
                  if(abs(d_tr0(n2)-d_tr0(n1))>1.e-20) then
                     !r_fac(n1) = 2*(grad_x(n1)*(xnd(n2)-xnd(n1))+grad_y(n1)*(ynd(n2)-ynd(n1)))/(d_tr0(n2)-d_tr0(n1))-1
                     r_fac(n1) = 2*(grad_x(n1)*(xel(id2,ie)-xel(id,ie))+grad_y(n1)*(yel(id2,ie)-yel(id,ie)))/(d_tr0(n2)-d_tr0(n1))-1
                     !tmp1 = 2*(grad_x(n1)*(xnd(n2)-xnd(n1))+grad_y(n1)*(ynd(n2)-ynd(n1)))
                     !write(12,*) 'case1',xnd(n1)-xnd(n2),ynd(n2)-ynd(n1),znd(n1)-znd(n2),xel(id,ie)-xel(id2,ie),yel(id,ie)-yel(id2,ie)
                     !write(12,*) 'case1',sumtmp,tmp1,d_tr0(n2),d_tr0(n1),r_fac(n1)
                  else
                     r_fac(n1) = -1
                  endif
                  psi_r = (r_fac(n1)+abs(r_fac(n1)))/(1.0d0+abs(r_fac(n1)))
                  flux_matrix1(indx1,n1,trc_n) = dt_dyn * vnor1 * (d_tr0(n1)+0.5*psi_r*(d_tr0(n2)-d_tr0(n1)))
                  flux_matrix2(indx1,n1,trc_n) = dt_dyn * vnor2 * (d_tr0(n1)+0.5*psi_r*(d_tr0(n2)-d_tr0(n1)))
                  flux_matrix3(indx2,n2,trc_n) = -1 * dt_dyn * vnor1 * (d_tr0(n1)+0.5*psi_r*(d_tr0(n2)-d_tr0(n1)))
                  flux_matrix4(indx2,n2,trc_n) = -1 * dt_dyn * vnor2 * (d_tr0(n1)+0.5*psi_r*(d_tr0(n2)-d_tr0(n1)))  
                  sumtmp = sumtmp * (d_tr0(n1)+0.5*psi_r*(d_tr0(n2)-d_tr0(n1)))             
                  !adv_threshold = exp(-20*tmp2)
                  !if((sumtmp/areatmp)>(d_tr0(n1)*adv_threshold))then
                  !if((sumtmp/areatmp)>(d_tr0(n1)*adv_threshold0(n1)).or.(sumtmp/areatmp)>(0.01))then
                  if((sumtmp/areatmp)>(d_tr0(n1)*adv_threshold).or.(sumtmp/areatmp)>(d_tr0(n1)*adv_threshold0(n2)))then
                     !tmp1 = areatmp*d_tr0(n1)/sumtmp*adv_threshold
                     tmp2 = min(d_tr0(n1)*adv_threshold,d_tr0(n1)*adv_threshold0(n2))
                     tmp1 = areatmp*tmp2/sumtmp
                     flux_matrix1(indx1,n1,trc_n) = flux_matrix1(indx1,n1,trc_n)*tmp1
                     flux_matrix2(indx1,n1,trc_n) = flux_matrix2(indx1,n1,trc_n)*tmp1
                     flux_matrix3(indx2,n2,trc_n) = flux_matrix3(indx2,n2,trc_n)*tmp1
                     flux_matrix4(indx2,n2,trc_n) = flux_matrix4(indx2,n2,trc_n)*tmp1
                  endif
                  
               else
                  tmp3 = d_tr0(n1)-d_tr0(n2)
                  tmp2=max(c0,tmp3)
                  !sumtmp = sumtmp * d_tr0(n2)
                  flux_matrix00(indx1,n1,1) = n2
                  flux_matrix00(indx2,n2,2) = n2
                  ! n2 center, n1 dowmstream
                  if(abs(d_tr0(n2)-d_tr0(n1))>1.e-20) then
                     !r_fac(n2) = 2*(grad_x(n2)*(xnd(n1)-xnd(n2))+grad_y(n2)*(ynd(n1)-ynd(n2)))/(d_tr0(n1)-d_tr0(n2))-1
                     r_fac(n2) = 2*(grad_x(n2)*(xel(id,ie)-xel(id2,ie))+grad_y(n2)*(yel(id,ie)-yel(id2,ie)))/(d_tr0(n1)-d_tr0(n2))-1
                     !tmp1 = 2*(grad_x(n2)*(xnd(n1)-xnd(n2))+grad_y(n2)*(ynd(n1)-ynd(n2)))
                     !write(12,*) 'case2',sumtmp,tmp1,d_tr0(n1),d_tr0(n2),r_fac(n2)
                  else
                     r_fac(n2) = -1
                  endif
                  psi_r = (r_fac(n2)+abs(r_fac(n2)))/(1.0d0+abs(r_fac(n2)))
                  flux_matrix1(indx1,n1,trc_n) = dt_dyn * vnor1 * (d_tr0(n2)+0.5*psi_r*(d_tr0(n1)-d_tr0(n2)))
                  flux_matrix2(indx1,n1,trc_n) = dt_dyn * vnor2 * (d_tr0(n2)+0.5*psi_r*(d_tr0(n1)-d_tr0(n2)))
                  flux_matrix3(indx2,n2,trc_n) = -1 * dt_dyn * vnor1 * (d_tr0(n2)+0.5*psi_r*(d_tr0(n1)-d_tr0(n2)))
                  flux_matrix4(indx2,n2,trc_n) = -1 * dt_dyn * vnor2 * (d_tr0(n2)+0.5*psi_r*(d_tr0(n1)-d_tr0(n2)))  
                  sumtmp = sumtmp * (d_tr0(n2)+0.5*psi_r*(d_tr0(n1)-d_tr0(n2)))               
                  !adv_threshold = exp(-20*tmp2)
                  !if((abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold))then
                  !if((abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold0(n2)).or.(abs(sumtmp)/areatmp)>(0.01))then
                  if((abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold).or.(abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold0(n1)))then
                     !tmp1 = areatmp*d_tr0(n2)/abs(sumtmp)*adv_threshold
                     tmp2 = min(d_tr0(n2)*adv_threshold,d_tr0(n2)*adv_threshold0(n1))
                     tmp1 = areatmp*tmp2/abs(sumtmp)
                     flux_matrix1(indx1,n1,trc_n) = flux_matrix1(indx1,n1,trc_n)*tmp1
                     flux_matrix2(indx1,n1,trc_n) = flux_matrix2(indx1,n1,trc_n)*tmp1
                     flux_matrix3(indx2,n2,trc_n) = flux_matrix3(indx2,n2,trc_n)*tmp1
                     flux_matrix4(indx2,n2,trc_n) = flux_matrix4(indx2,n2,trc_n)*tmp1
                  endif
               endif
            else
               xtmp = tmpy3 - yctr2
               ytmp = xctr2 - tmpx3
               vnor1=utmp*xtmp+vtmp*ytmp !normal vel x length 
               vnor2=swild10(j,1)*xtmp+swild10(j,2)*ytmp !normal vel@side x length 
               !vnor1=vnor2

               sumtmp = (dt_dyn*vnor1+dt_dyn*vnor2)/2
               if(sumtmp>0) then
                  tmp3 = d_tr0(n2)-d_tr0(n1)
                  tmp2=max(c0,tmp3)
                  !sumtmp = sumtmp * d_tr0(n1)
                  flux_matrix00(indx1,n1,2) = n1
                  flux_matrix00(indx2,n2,1) = n1
                  ! n1 center, n2 dowmstream
                  if(abs(d_tr0(n2)-d_tr0(n1))>1.e-20) then
                     !r_fac(n1) = 2*(grad_x(n1)*(xnd(n2)-xnd(n1))+grad_y(n1)*(ynd(n2)-ynd(n1)))/(d_tr0(n2)-d_tr0(n1))-1
                     r_fac(n1) = 2*(grad_x(n1)*(xel(id2,ie)-xel(id,ie))+grad_y(n1)*(yel(id2,ie)-yel(id,ie)))/(d_tr0(n2)-d_tr0(n1))-1
                     !tmp1 = 2*(grad_x(n1)*(xnd(n2)-xnd(n1))+grad_y(n1)*(ynd(n2)-ynd(n1)))
                     !write(12,*) 'case3',sumtmp,tmp1,d_tr0(n2),d_tr0(n1),r_fac(n1)
                  else
                     r_fac(n1) = -1
                  endif
                  psi_r = (r_fac(n1)+abs(r_fac(n1)))/(1.0d0+abs(r_fac(n1)))
                  flux_matrix3(indx1,n1,trc_n) = dt_dyn*vnor1 * (d_tr0(n1)+0.5*psi_r*(d_tr0(n2)-d_tr0(n1)))
                  flux_matrix4(indx1,n1,trc_n) = dt_dyn*vnor2 * (d_tr0(n1)+0.5*psi_r*(d_tr0(n2)-d_tr0(n1)))
                  flux_matrix1(indx2,n2,trc_n) = -1*dt_dyn*vnor1 * (d_tr0(n1)+0.5*psi_r*(d_tr0(n2)-d_tr0(n1)))
                  flux_matrix2(indx2,n2,trc_n) = -1*dt_dyn*vnor2 * (d_tr0(n1)+0.5*psi_r*(d_tr0(n2)-d_tr0(n1)))
                  sumtmp = sumtmp * (d_tr0(n1)+0.5*psi_r*(d_tr0(n2)-d_tr0(n1)))
                  !adv_threshold = exp(-20*tmp2)
                  !if((sumtmp/areatmp)>(d_tr0(n1)*adv_threshold))then
                  !if((sumtmp/areatmp)>(d_tr0(n1)*adv_threshold0(n1)).or.(sumtmp/areatmp)>(0.01))then
                  if((sumtmp/areatmp)>(d_tr0(n1)*adv_threshold).or.(sumtmp/areatmp)>(d_tr0(n1)*adv_threshold0(n2)))then
                     !tmp1 = areatmp*d_tr0(n1)/sumtmp*adv_threshold
                     tmp2 = min(d_tr0(n1)*adv_threshold,d_tr0(n1)*adv_threshold0(n2))
                     tmp1 = areatmp*tmp2/sumtmp
                     flux_matrix3(indx1,n1,trc_n) = flux_matrix3(indx1,n1,trc_n)*tmp1
                     flux_matrix4(indx1,n1,trc_n) = flux_matrix4(indx1,n1,trc_n)*tmp1
                     flux_matrix1(indx2,n2,trc_n) = flux_matrix1(indx2,n2,trc_n)*tmp1
                     flux_matrix2(indx2,n2,trc_n) = flux_matrix2(indx2,n2,trc_n)*tmp1
                  endif
                  
               else
                  tmp3 = d_tr0(n1)-d_tr0(n2)
                  tmp2=max(c0,tmp3)
                  !sumtmp = sumtmp * d_tr0(n2)
                  flux_matrix00(indx1,n1,2) = n2
                  flux_matrix00(indx2,n2,1) = n2
                  ! n2 center, n1 dowmstream
                  
                  if(abs(d_tr0(n1)-d_tr0(n2))>1.e-20) then
                     !r_fac(n2) = 2*(grad_x(n2)*(xnd(n1)-xnd(n2))+grad_y(n2)*(ynd(n1)-ynd(n2)))/(d_tr0(n1)-d_tr0(n2))-1
                     r_fac(n2) = 2*(grad_x(n2)*(xel(id,ie)-xel(id2,ie))+grad_y(n2)*(yel(id,ie)-yel(id2,ie)))/(d_tr0(n1)-d_tr0(n2))-1
                     !tmp1 = 2*(grad_x(n2)*(xnd(n1)-xnd(n2))+grad_y(n2)*(ynd(n1)-ynd(n2)))
                     !write(12,*) 'case4',sumtmp,tmp1,d_tr0(n1),d_tr0(n2),r_fac(n2)
                  else
                     r_fac(n2) = -1
                  endif
                  psi_r = (r_fac(n2)+abs(r_fac(n2)))/(1.0d0+abs(r_fac(n2)))
                  flux_matrix3(indx1,n1,trc_n) = dt_dyn * vnor1 * (d_tr0(n2)+0.5*psi_r*(d_tr0(n1)-d_tr0(n2)))
                  flux_matrix4(indx1,n1,trc_n) = dt_dyn * vnor2 * (d_tr0(n2)+0.5*psi_r*(d_tr0(n1)-d_tr0(n2)))
                  flux_matrix1(indx2,n2,trc_n) = -1 * dt_dyn * vnor1 * (d_tr0(n2)+0.5*psi_r*(d_tr0(n1)-d_tr0(n2)))
                  flux_matrix2(indx2,n2,trc_n) = -1 * dt_dyn * vnor2 * (d_tr0(n2)+0.5*psi_r*(d_tr0(n1)-d_tr0(n2)))
                  sumtmp = sumtmp * (d_tr0(n2)+0.5*psi_r*(d_tr0(n1)-d_tr0(n2)))
                  !adv_threshold = exp(-20*tmp2)
                  !if((abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold))then
                  !if((abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold0(n2)).or.(abs(sumtmp)/areatmp)>0.01)then
                  if((abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold).or.(abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold0(n1)))then
                     !tmp1 = areatmp*d_tr0(n2)/abs(sumtmp)*adv_threshold
                     tmp2 = min(d_tr0(n2)*adv_threshold,d_tr0(n2)*adv_threshold0(n1))
                     tmp1 = areatmp*tmp2/abs(sumtmp)
                     flux_matrix3(indx1,n1,trc_n) = flux_matrix3(indx1,n1,trc_n)*tmp1
                     flux_matrix4(indx1,n1,trc_n) = flux_matrix4(indx1,n1,trc_n)*tmp1
                     flux_matrix1(indx2,n2,trc_n) = flux_matrix1(indx2,n2,trc_n)*tmp1
                     flux_matrix2(indx2,n2,trc_n) = flux_matrix2(indx2,n2,trc_n)*tmp1
                  endif
               endif
            endif
          enddo
    enddo
   do i=1,nx_nh
      suma=0.d0 !sum of areas
      fluxchan0=0.d0
      do j=1,nne(i)
         ie=indel(j,i)
         suma=suma+area(ie)/3 !approx for quad
         fluxchan0=fluxchan0 + (flux_matrix1(j,i,trc_n) + flux_matrix2(j,i,trc_n))/2.d0
         fluxchan0=fluxchan0 + (flux_matrix3(j,i,trc_n) + flux_matrix4(j,i,trc_n))/2.d0

      enddo
         if( suma/=0 ) then
            trl0(i)=fluxchan0/suma !m/s/s
         endif !suma/=0  
         trc(i) = d_tr0(i) - trl0(i) !+ d_tr(i)   
   enddo  
    call exchange_p2d(trc)
  end subroutine tvd_icepack

  !=======================================================================
subroutine tvd_icepack_other(trc,trc_o,trc_n,trc_d,trc_d1,trc_d2)
        
        !use mod_mesh     
     
        implicit none
        integer (kind=int_kind), intent(in) :: trc_n,trc_d
        real (kind=dbl_kind), dimension(nx), intent(inout) :: trc   ,trc_o
        real (kind=dbl_kind), dimension(nx), intent(in), optional :: trc_d1,trc_d2
        real   (kind=dbl_kind)            :: fluxchan1,fluxchan2,suma,xctr2,yctr2,vnor1,vnor2,&
                                             & tmpx2,tmpx3,tmpy2,tmpy3,utmp,vtmp,xtmp,ytmp,vnorm,&
                                             & trctmp  , trcmin, trcmax, tmpnode   
        integer (kind=int_kind) :: i,j,ie,id,id2,id3,isd3,isd2,     &
                                   m,n1,n2,jj
        logical (kind=log_kind) :: flag_noadv
         real(rkind) :: swild10(3,2),swild(3),trl0(nx),d_tr0(nx)
        !type(t_mesh),        target,        intent(in)    :: mesh

      trl0(:) = c0
      d_tr0(:) = trc(:)
      flux_matrix1(:,:,trc_n) = c0
      flux_matrix2(:,:,trc_n) = c0
      flux_matrix3(:,:,trc_n) = c0
      flux_matrix4(:,:,trc_n) = c0
      flux_matrix0(:,:,:,trc_n) = c0
        ! Driving sequence
      do i=1,nx_nh
        if(idry(i)==1) cycle
        !if(isbnd(1,i)/=0) cycle
        !Wet node
          fluxchan1=0.d0 !sum of fluxes;\int_\Gamma u*u_n d\Gamma
          fluxchan2=0.d0 !sum of fluxes;
          suma=0.d0 !sum of areas
          do j=1,nne(i)
            ie=indel(j,i)
            if(idry_e(ie)==1) cycle
 
            !Wet elem
            suma=suma+area(ie)/3 !approx for quad
            id=iself(j,i)
            id3=nxq(i34(ie)-1,id,i34(ie)) !adjacent side index
            id2=nxq(i34(ie)-2,id,i34(ie)) !adjacent side index
            isd3=elside(id3,ie)
            isd2=elside(id2,ie)

            
            do m=1,i34(ie) !side
                n1=isidenode(1,elside(m,ie))
                n2=isidenode(2,elside(m,ie))
                tmpx2 = trc_o(n1)
                tmpy2 = trc_o(n2)
                if (present(trc_d1)) then
                  tmpx2 = tmpx2 * trc_d1(n1)
                  tmpy2 = tmpy2 * trc_d1(n2)
                endif
                if (present(trc_d2)) then
                  tmpx2 = tmpx2 * trc_d2(n1)
                  tmpy2 = tmpy2 * trc_d2(n2)
                endif
                swild(m) = (tmpx2+tmpy2)*0.5
            enddo
            trctmp = sum(swild(1:i34(ie)))/3
            trcmin = minval(trc_o(indnd(1:nnp(i),i)))
            trcmin = min(trcmin,trc_o(i))
            trcmax = maxval(trc_o(indnd(1:nnp(i),i)))
            trcmax = max(trcmax,trc_o(i))
            
            tmpnode = flux_matrix00(j,i,1)
            tmpx2 = trc_o(tmpnode)
            if (present(trc_d1)) then
               tmpx2 = tmpx2 * trc_d1(tmpnode)
            endif
            if (present(trc_d2)) then
               tmpx2 = tmpx2 * trc_d2(tmpnode)
            endif
            !1st segment
            flux_matrix1(j,i,trc_n) = flux_matrix1(j,i,trc_d)*tmpx2
            flux_matrix2(j,i,trc_n) = flux_matrix2(j,i,trc_d)*tmpx2
            fluxchan1=fluxchan1+(flux_matrix1(j,i,trc_d)*tmpx2+flux_matrix2(j,i,trc_d)*tmpx2)/2.d0
             
            
            !2nd segment
            !Normal dir x length
!            xtmp=ycj(isd2)-yctr(ie)
!            ytmp=xctr(ie)-xcj(isd2)
            tmpnode = flux_matrix00(j,i,2)
            tmpx2 = trc_o(tmpnode)
            if (present(trc_d1)) then
               tmpx2 = tmpx2 * trc_d1(tmpnode)
            endif
            if (present(trc_d2)) then
               tmpx2 = tmpx2 * trc_d2(tmpnode)
            endif
            flux_matrix3(j,i,trc_n) = flux_matrix3(j,i,trc_d)*tmpx2
            flux_matrix4(j,i,trc_n) = flux_matrix4(j,i,trc_d)*tmpx2
            fluxchan1=fluxchan1+(flux_matrix3(j,i,trc_d)*tmpx2+flux_matrix4(j,i,trc_d)*tmpx2)/2.d0
          enddo !j=1,nne(i)
          
          if(suma/=0) then
            trl0(i)=fluxchan1/suma !m/s/s
          endif !suma/=0
      enddo !i=1,np
      call exchange_p2d(trc)
      call exchange_p2d(d_tr0)
      call exchange_p2d(trl0)
      do i = 1,nx_nh
       trc(i) = d_tr0(i) - trl0(i) !+ d_tr(i)
      enddo

      call exchange_p2d(trc)
    end subroutine tvd_icepack_other

    !=======================================================================
    subroutine tvd_icepack_dep(trc,trc_o,trc_n,trc_d)
        
        !use mod_mesh     
     
        implicit none
        integer (kind=int_kind), intent(in) :: trc_n,trc_d
        real (kind=dbl_kind), dimension(nx), intent(inout) :: trc   ,trc_o
        real   (kind=dbl_kind)            :: fluxchan1,fluxchan2,suma,xctr2,yctr2,vnor1,vnor2,&
                                             & tmpx2,tmpx3,tmpy2,tmpy3,utmp,vtmp,xtmp,ytmp,vnorm,&
                                             & trctmp, trcmin, trcmax ,sumtmp,areatmp, tmpnode  
        integer (kind=int_kind) :: i,j,ie,id,id2,id3,isd3,isd2,     &
                                   m,n1,n2,jj,indx,nd,nodeid
        logical (kind=log_kind) :: flag_noadv
         real(rkind) :: swild10(3,2),swild(3),trl0(nx),d_tr0(nx)
        !type(t_mesh),        target,        intent(in)    :: mesh

      trl0(:) = c0
      d_tr0(:) = trc(:)
      flux_matrix1(:,:,trc_n) = c0
      flux_matrix2(:,:,trc_n) = c0
      flux_matrix3(:,:,trc_n) = c0
      flux_matrix4(:,:,trc_n) = c0
      flux_matrix0(:,:,:,trc_n) = c0
        ! Driving sequence
      do i=1,nx_nh
        if(idry(i)==1) cycle
        !Wet node
          fluxchan1=0.d0 !sum of fluxes;\int_\Gamma u*u_n d\Gamma
          fluxchan2=0.d0 !sum of fluxes;
          suma=0.d0 !sum of areas
          do j=1,nne(i)
            ie=indel(j,i)
            if(idry_e(ie)==1) cycle
 
            !Wet elem
            areatmp = area(ie)/6
            suma=suma+area(ie)/3 !approx for quad
            id=iself(j,i)
            id3=nxq(i34(ie)-1,id,i34(ie)) !adjacent side index
            id2=nxq(i34(ie)-2,id,i34(ie)) !adjacent side index
            isd3=elside(id3,ie)
            isd2=elside(id2,ie)

            do m=1,i34(ie) !side
                n1=isidenode(1,elside(m,ie))
                n2=isidenode(2,elside(m,ie))
                tmpx2 = trc_o(n1)
                tmpy2 = trc_o(n2)
                swild(m) = (tmpx2+tmpy2)*0.5
            enddo
            trctmp=sum(swild(1:i34(ie)))/3
            !1st segment
            tmpnode = flux_matrix00(j,i,1)
            tmpx2 = trc_o(tmpnode)
            flux_matrix1(j,i,trc_n) = flux_matrix1(j,i,trc_d)*tmpx2
            flux_matrix2(j,i,trc_n) = flux_matrix2(j,i,trc_d)*tmpx2
            sumtmp = (flux_matrix1(j,i,trc_n)+flux_matrix2(j,i,trc_n))/2.d0
            nodeid = elnode(nxq(1,id,i34(ie)),ie)
            fluxchan1=fluxchan1+(flux_matrix1(j,i,trc_n)+flux_matrix2(j,i,trc_n))/2.d0
             
            
            !2nd segment
            !Normal dir x length
!            xtmp=ycj(isd2)-yctr(ie)
!            ytmp=xctr(ie)-xcj(isd2)
            tmpnode = flux_matrix00(j,i,2)
            tmpx2 = trc_o(tmpnode)
            flux_matrix3(j,i,trc_n) = flux_matrix3(j,i,trc_d)*tmpx2
            flux_matrix4(j,i,trc_n) = flux_matrix4(j,i,trc_d)*tmpx2
            sumtmp = (flux_matrix3(j,i,trc_n) + flux_matrix4(j,i,trc_n))/2.d0
            nodeid = elnode(nxq(i34(ie)-1,id,i34(ie)),ie)

            fluxchan1=fluxchan1+(flux_matrix3(j,i,trc_n) + flux_matrix4(j,i,trc_n))/2.d0
          enddo !j=1,nne(i)
          
          if(suma/=0) then
            trl0(i)=fluxchan1/suma !m/s/s
          endif !suma/=0
      enddo !i=1,np
      call exchange_p2d(trc)
      call exchange_p2d(d_tr0)
      call exchange_p2d(trl0)
      do i = 1,nx_nh
       trc(i) = d_tr0(i) - trl0(i) !+ d_tr(i)


      enddo

      call exchange_p2d(trc)
    end subroutine tvd_icepack_dep

 !=======================================================================
    
    subroutine tvdu_icepack(trc,trc_n,adv_threshold0)
        
      !use mod_mesh     
       use mice_module
       use schism_glbl, only: dldxy
      implicit none
      integer (kind=int_kind), intent(in) :: trc_n
      real(kind=dbl_kind), dimension(nx), intent(inout) :: trc,adv_threshold0  
      real   (kind=dbl_kind)            :: suma,xctr2,yctr2,vnor1,vnor2,sumtmp,areatmp,&
                                           & tmpx2,tmpx3,tmpy2,tmpy3,utmp,vtmp,xtmp,ytmp,vnorm,&
                                           & trctmp, trcmin, trcmax ,puny ,tmp1 ,fluxchan0, tmp2, adv_threshold,&
                                           & tmp3, psi_r, sum1, sum2
      integer (kind=int_kind) :: i,j,ie,id,id2,id3,isd3,isd2,     &
                                 m,n1,n2,jj,indx1,indx2,nd,nodeid
      logical (kind=log_kind) :: flag_noadv
       real(kind=dbl_kind) :: swild10(3,2),swild(3),swild2(3,2),trl0(nx),d_tr0(nx),fluxchan1(nx),fluxchan2(nx), &
                              & r_fac(nx),grad_x(nx),grad_y(nx),trcmaxn(nx),trcminn(nx)  
      !type(t_mesh),        target,        intent(in)    :: mesh
       call icepack_query_parameters(puny_out=puny)

    adv_threshold = 10
    trl0(:) = c0
    fluxchan1(:) = c0 !positive flux
    fluxchan2(:) = c0 !negative flux
    d_tr0(:) = trc(:)
    flux_matrix1(:,:,trc_n) = c0
    flux_matrix2(:,:,trc_n) = c0
    flux_matrix3(:,:,trc_n) = c0
    flux_matrix4(:,:,trc_n) = c0
    flux_matrix00(:,:,:) = c0
    flux_matrix0(:,:,:,trc_n) = c0
    r_fac(:) = c0
    grad_x(:) = c0
    grad_y(:) = c0
    trcmaxn(:) = c0
    trcminn(:) = c0

    do i = 1,nx_nh
    sum1=0; sum2=0
    trcmaxn(i) = maxval(trc(indnd(1:nnp(i),i)))
    trcmaxn(i) = max(trcmaxn(i),trc(i))
    trcminn(i) = minval(trc(indnd(1:nnp(i),i)))
    trcminn(i) = min(trcminn(i),trc(i))
      do j=1,nne(i)
         ie=indel(j,i)
         id=iself(j,i)
         if(ics==1) then
            sum1=sum1+area(ie)/3*dot_product(trc(elnode(1:i34(ie),ie)),dldxy(1:i34(ie),1,ie))
            sum2=sum2+area(ie)/3*dot_product(trc(elnode(1:i34(ie),ie)),dldxy(1:i34(ie),2,ie))
         else
          sum1=sum1+area(ie)/3*dot_product(trc(elnode(1:i34(ie),ie)),dldxy(1:i34(ie),1,ie))*dot_product(eframe(1:3,1,ie),pframe(1:3,1,i))+&
                   &area(ie)/3*dot_product(trc(elnode(1:i34(ie),ie)),dldxy(1:i34(ie),2,ie))*dot_product(eframe(1:3,2,ie),pframe(1:3,1,i))
          sum2=sum2+area(ie)/3*dot_product(trc(elnode(1:i34(ie),ie)),dldxy(1:i34(ie),1,ie))*dot_product(eframe(1:3,1,ie),pframe(1:3,2,i))+&
                   &area(ie)/3*dot_product(trc(elnode(1:i34(ie),ie)),dldxy(1:i34(ie),2,ie))*dot_product(eframe(1:3,2,ie),pframe(1:3,2,i))
         endif
      enddo
      grad_x(i) = sum1/area_median(i)
      grad_y(i) = sum2/area_median(i)
    enddo
    call exchange_p2d(grad_x)
    call exchange_p2d(grad_y)
    call exchange_p2d(trcmaxn)
    call exchange_p2d(trcminn)


      ! Driving sequence
    do ie=1,nx_elem

          if(idry_e(ie)==1) cycle
          areatmp = area(ie)/6
          !Wet elem

          do m=1,i34(ie) !side
            n1=isidenode(1,elside(m,ie))
            n2=isidenode(2,elside(m,ie))
            if(ics==1) then
               swild10(m,1) = 0.5*(uvel(n1)+uvel(n2))
               swild10(m,2) = 0.5*(vvel(n1)+vvel(n2))
               swild(m) = (trc(n1)+trc(n2))*0.5
               swild2(m,1) = uvel(elnode(m,ie))
               swild2(m,2) = vvel(elnode(m,ie))
            else
               swild10(m,1) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),eframe(:,1,ie))+vvel(n1)*dot_product(pframe(:,2,n1),eframe(:,1,ie))&
                                 &+uvel(n2)*dot_product(pframe(:,1,n2),eframe(:,1,ie))+vvel(n2)*dot_product(pframe(:,2,n2),eframe(:,1,ie)))
               swild10(m,2) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),eframe(:,2,ie))+vvel(n1)*dot_product(pframe(:,2,n1),eframe(:,2,ie))&
                                 &+uvel(n2)*dot_product(pframe(:,1,n2),eframe(:,2,ie))+vvel(n2)*dot_product(pframe(:,2,n2),eframe(:,2,ie)))                      
               swild(m) = (trc(n1)+trc(n2))*0.5
               swild2(m,1) = uvel(elnode(m,ie))*dot_product(pframe(:,1,elnode(m,ie)),eframe(:,1,ie))+vvel(elnode(m,ie))*dot_product(pframe(:,2,elnode(m,ie)),eframe(:,1,ie))
               swild2(m,2) = uvel(elnode(m,ie))*dot_product(pframe(:,1,elnode(m,ie)),eframe(:,2,ie))+vvel(elnode(m,ie))*dot_product(pframe(:,2,elnode(m,ie)),eframe(:,2,ie))
            endif
            !swild10(m,1) = 0.5*(uvel(n1)*trc(n1)*dot_product(pframe(:,1,n1),sframe2(:,1,elside(m,ie)))+vvel(n1)*trc(n1)*dot_product(pframe(:,2,n1),sframe2(:,1,elside(m,ie)))&
            !                 &+uvel(n2)*trc(n2)*dot_product(pframe(:,1,n2),sframe2(:,1,elside(m,ie)))+vvel(n2)*trc(n2)*dot_product(pframe(:,2,n1),sframe2(:,1,elside(m,ie))))
            !swild10(m,2) = 0.5*(uvel(n1)*trc(n1)*dot_product(pframe(:,1,n1),sframe2(:,2,elside(m,ie)))+vvel(n1)*trc(n1)*dot_product(pframe(:,2,n1),sframe2(:,2,elside(m,ie)))&
            !                 &+uvel(n2)*trc(n2)*dot_product(pframe(:,1,n2),sframe2(:,2,elside(m,ie)))+vvel(n2)*trc(n2)*dot_product(pframe(:,2,n1),sframe2(:,2,elside(m,ie))))
            !swild10(m,1) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),sframe2(:,1,elside(m,ie)))+vvel(n1)*dot_product(pframe(:,2,n1),sframe2(:,1,elside(m,ie)))&
            !                 &+uvel(n2)*dot_product(pframe(:,1,n2),sframe2(:,1,elside(m,ie)))+vvel(n2)*dot_product(pframe(:,2,n2),sframe2(:,1,elside(m,ie))))
            !swild10(m,2) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),sframe2(:,2,elside(m,ie)))+vvel(n1)*dot_product(pframe(:,2,n1),sframe2(:,2,elside(m,ie)))&
            !                 &+uvel(n2)*dot_product(pframe(:,1,n2),sframe2(:,2,elside(m,ie)))+vvel(n2)*dot_product(pframe(:,2,n2),sframe2(:,2,elside(m,ie))))
            
         enddo
        !utmp = sum(swild10(1:i34(ie),1))/3 !vel @ centroid
        !vtmp = sum(swild10(1:i34(ie),2))/3 !vel @ centroid
        utmp=sum(swild2(1:i34(ie),1))/3 !vel @ centroid
        vtmp=sum(swild2(1:i34(ie),2))/3 !vel @ centroid
        trctmp = sum(swild(1:i34(ie)))/3

          do j =1,i34(ie) !side
            n1=isidenode(1,elside(j,ie))
            n2=isidenode(2,elside(j,ie))
            
            indx1 = 0
            do m=1,nne(n1)
               if(indel(m,n1)==ie) then
                  indx1=m; exit
               endif
            enddo
            if(indx1==0) call parallel_abort('STEP: failed to find (tvdu_icepack_1)') 
            indx2 = 0
            do m=1,nne(n2)
               if(indel(m,n2)==ie) then
                  indx2=m; exit
               endif
            enddo
            if(indx2==0) call parallel_abort('STEP: failed to find (tvdu_icepack_2)') 

            id  = iself(indx1,n1)
            id2 = iself(indx2,n2)
            id3 = nxq(i34(ie)-2,id,i34(ie))

            if(ics==1) then
               xctr2=xctr(ie); yctr2=yctr(ie)
               tmpx3=(xnd(n1)+xnd(n2))/2.d0
               tmpy3=(ynd(n1)+ynd(n2))/2.d0
            else!ll; use [xy]el defined in eframe
               xctr2=0.d0 !sum(xel(elnode(1:i34(ie),ie)))/i34(ie)
               yctr2=0.d0 !sum(yel(elnode(1:i34(ie),ie)))/i34(ie)
               tmpx3=(xel(id,ie)+xel(id2,ie))/2.d0
               tmpy3=(yel(id,ie)+yel(id2,ie))/2.d0
            endif!ics
            !xtmp=yctr2-tmpy3
            !ytmp=tmpx3-xctr2
            !vnor1=utmp*xtmp+vtmp*ytmp !normal vel x length 
            !vnor2=swild10(j,1)*xtmp+swild10(j,2)*ytmp !normal vel@side x length 
            if(id3.eq.id2) then
               xtmp=yctr2-tmpy3
               ytmp=tmpx3-xctr2
               vnor1=utmp*xtmp+vtmp*ytmp !normal vel x length 
               vnor2=swild10(j,1)*xtmp+swild10(j,2)*ytmp !normal vel@side x length 
               !vnor1=vnor2

               sumtmp = (dt_dyn*vnor1+dt_dyn*vnor2)/2
               if(sumtmp>0) then
                  tmp3 = d_tr0(n2)-d_tr0(n1)
                  tmp2=max(c0,tmp3)
                  !sumtmp = sumtmp * d_tr0(n1)
                  flux_matrix00(indx1,n1,1) = n1
                  flux_matrix00(indx2,n2,2) = n1
                  ! n1 center, n2 dowmstream
                  if(abs(d_tr0(n2)-d_tr0(n1))>1.e-20) then
                     !r_fac(n1) = 2*(grad_x(n1)*(xnd(n2)-xnd(n1))+grad_y(n1)*(ynd(n2)-ynd(n1)))/(d_tr0(n2)-d_tr0(n1))-1
                     r_fac(n1) = 2*(grad_x(n1)*(xel(id2,ie)-xel(id,ie))+grad_y(n1)*(yel(id2,ie)-yel(id,ie)))/(d_tr0(n2)-d_tr0(n1))-1
                     tmp1 = d_tr0(n2) - 2*(grad_x(n1)*(xel(id2,ie)-xel(id,ie))+grad_y(n1)*(yel(id2,ie)-yel(id,ie)))
                     !if(tmp1>trcmaxn(n1)) r_fac(n1)=(d_tr0(n1)-trcmaxn(n1))/(d_tr0(n2)-d_tr0(n1))  
                     !if(tmp1<trcminn(n1)) r_fac(n1)=(d_tr0(n1)-trcminn(n1))/(d_tr0(n2)-d_tr0(n1))   
                     if(tmp1>1.0d0) r_fac(n1)=(d_tr0(n1)-1.0d0)/(d_tr0(n2)-d_tr0(n1))  
                     if(tmp1<0.0d0) r_fac(n1)=(d_tr0(n1)-0.0d0)/(d_tr0(n2)-d_tr0(n1))             
                     !tmp1 = 2*(grad_x(n1)*(xnd(n2)-xnd(n1))+grad_y(n1)*(ynd(n2)-ynd(n1)))
                     !write(12,*) 'case1',xnd(n1)-xnd(n2),ynd(n2)-ynd(n1),znd(n1)-znd(n2),xel(id,ie)-xel(id2,ie),yel(id,ie)-yel(id2,ie)
                     !write(12,*) 'case1',sumtmp,tmp1,d_tr0(n2),d_tr0(n1),r_fac(n1)
                  else
                     r_fac(n1) = -1
                  endif
                  psi_r = (r_fac(n1)+abs(r_fac(n1)))/(1.0d0+abs(r_fac(n1)))
                  !psi_r = max(0.d0,min(1.d0,2.0d0*r_fac(n1)),min(2.d0,r_fac(n1)))
                  !psi_r = max(0.d0,min(1.d0,r_fac(n1)))
                  flux_matrix1(indx1,n1,trc_n) = dt_dyn * vnor1 * (d_tr0(n1)+0.5*psi_r*(d_tr0(n2)-d_tr0(n1)))
                  flux_matrix2(indx1,n1,trc_n) = dt_dyn * vnor2 * (d_tr0(n1)+0.5*psi_r*(d_tr0(n2)-d_tr0(n1)))
                  flux_matrix3(indx2,n2,trc_n) = -1 * dt_dyn * vnor1 * (d_tr0(n1)+0.5*psi_r*(d_tr0(n2)-d_tr0(n1)))
                  flux_matrix4(indx2,n2,trc_n) = -1 * dt_dyn * vnor2 * (d_tr0(n1)+0.5*psi_r*(d_tr0(n2)-d_tr0(n1)))  
                  sumtmp = sumtmp * (d_tr0(n1)+0.5*psi_r*(d_tr0(n2)-d_tr0(n1)))             
                  !adv_threshold = exp(-20*tmp2)
                  if((sumtmp/areatmp)>(d_tr0(n1)*adv_threshold))then
                  !if((sumtmp/areatmp)>(d_tr0(n1)*adv_threshold0(n1)).or.(sumtmp/areatmp)>(0.01))then
                  !if((sumtmp/areatmp)>(d_tr0(n1)*adv_threshold).or.(sumtmp/areatmp)>(d_tr0(n1)*adv_threshold0(n2)))then
                     tmp1 = areatmp*d_tr0(n1)/sumtmp*adv_threshold
                     !tmp2 = min(d_tr0(n1)*adv_threshold,d_tr0(n1)*adv_threshold0(n2))
                     !tmp1 = areatmp*tmp2/sumtmp
                     flux_matrix1(indx1,n1,trc_n) = flux_matrix1(indx1,n1,trc_n)*tmp1
                     flux_matrix2(indx1,n1,trc_n) = flux_matrix2(indx1,n1,trc_n)*tmp1
                     flux_matrix3(indx2,n2,trc_n) = flux_matrix3(indx2,n2,trc_n)*tmp1
                     flux_matrix4(indx2,n2,trc_n) = flux_matrix4(indx2,n2,trc_n)*tmp1
                  endif
                  
               else
                  tmp3 = d_tr0(n1)-d_tr0(n2)
                  tmp2=max(c0,tmp3)
                  !sumtmp = sumtmp * d_tr0(n2)
                  flux_matrix00(indx1,n1,1) = n2
                  flux_matrix00(indx2,n2,2) = n2
                  ! n2 center, n1 dowmstream
                  if(abs(d_tr0(n2)-d_tr0(n1))>1.e-20) then
                     !r_fac(n2) = 2*(grad_x(n2)*(xnd(n1)-xnd(n2))+grad_y(n2)*(ynd(n1)-ynd(n2)))/(d_tr0(n1)-d_tr0(n2))-1
                     r_fac(n2) = 2*(grad_x(n2)*(xel(id,ie)-xel(id2,ie))+grad_y(n2)*(yel(id,ie)-yel(id2,ie)))/(d_tr0(n1)-d_tr0(n2))-1
                     tmp1 = d_tr0(n1) - 2*(grad_x(n2)*(xel(id,ie)-xel(id2,ie))+grad_y(n2)*(yel(id,ie)-yel(id2,ie)))
                     !if(tmp1>trcmaxn(n2)) r_fac(n2)=(d_tr0(n2)-trcmaxn(n2))/(d_tr0(n1)-d_tr0(n2))
                     !if(tmp1<trcminn(n2)) r_fac(n2)=(d_tr0(n2)-trcminn(n2))/(d_tr0(n1)-d_tr0(n2))
                     if(tmp1>1.0d0) r_fac(n2)=(d_tr0(n2)-1.0d0)/(d_tr0(n1)-d_tr0(n2))
                     if(tmp1<0.0d0) r_fac(n2)=(d_tr0(n2)-0.0d0)/(d_tr0(n1)-d_tr0(n2))
                     !tmp1 = 2*(grad_x(n2)*(xnd(n1)-xnd(n2))+grad_y(n2)*(ynd(n1)-ynd(n2)))
                     !write(12,*) 'case2',sumtmp,tmp1,d_tr0(n1),d_tr0(n2),r_fac(n2)
                  else
                     r_fac(n2) = -1
                  endif
                  psi_r = (r_fac(n2)+abs(r_fac(n2)))/(1.0d0+abs(r_fac(n2)))
                  !psi_r = max(0.d0,min(1.d0,2.0d0*r_fac(n2)),min(2.d0,r_fac(n2)))
                  !psi_r = max(0.d0,min(1.d0,r_fac(n2)))
                  flux_matrix1(indx1,n1,trc_n) = dt_dyn * vnor1 * (d_tr0(n2)+0.5*psi_r*(d_tr0(n1)-d_tr0(n2)))
                  flux_matrix2(indx1,n1,trc_n) = dt_dyn * vnor2 * (d_tr0(n2)+0.5*psi_r*(d_tr0(n1)-d_tr0(n2)))
                  flux_matrix3(indx2,n2,trc_n) = -1 * dt_dyn * vnor1 * (d_tr0(n2)+0.5*psi_r*(d_tr0(n1)-d_tr0(n2)))
                  flux_matrix4(indx2,n2,trc_n) = -1 * dt_dyn * vnor2 * (d_tr0(n2)+0.5*psi_r*(d_tr0(n1)-d_tr0(n2)))  
                  sumtmp = sumtmp * (d_tr0(n2)+0.5*psi_r*(d_tr0(n1)-d_tr0(n2)))               
                  !adv_threshold = exp(-20*tmp2)
                  if((abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold))then
                  !if((abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold0(n2)).or.(abs(sumtmp)/areatmp)>(0.01))then
                  !if((abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold).or.(abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold0(n1)))then
                     tmp1 = areatmp*d_tr0(n2)/abs(sumtmp)*adv_threshold
                     !tmp2 = min(d_tr0(n2)*adv_threshold,d_tr0(n2)*adv_threshold0(n1))
                     !tmp1 = areatmp*tmp2/abs(sumtmp)
                     flux_matrix1(indx1,n1,trc_n) = flux_matrix1(indx1,n1,trc_n)*tmp1
                     flux_matrix2(indx1,n1,trc_n) = flux_matrix2(indx1,n1,trc_n)*tmp1
                     flux_matrix3(indx2,n2,trc_n) = flux_matrix3(indx2,n2,trc_n)*tmp1
                     flux_matrix4(indx2,n2,trc_n) = flux_matrix4(indx2,n2,trc_n)*tmp1
                  endif
               endif
            else
               xtmp = tmpy3 - yctr2
               ytmp = xctr2 - tmpx3
               vnor1=utmp*xtmp+vtmp*ytmp !normal vel x length 
               vnor2=swild10(j,1)*xtmp+swild10(j,2)*ytmp !normal vel@side x length 
               !vnor1=vnor2

               sumtmp = (dt_dyn*vnor1+dt_dyn*vnor2)/2
               if(sumtmp>0) then
                  tmp3 = d_tr0(n2)-d_tr0(n1)
                  tmp2=max(c0,tmp3)
                  !sumtmp = sumtmp * d_tr0(n1)
                  flux_matrix00(indx1,n1,2) = n1
                  flux_matrix00(indx2,n2,1) = n1
                  ! n1 center, n2 dowmstream
                  if(abs(d_tr0(n2)-d_tr0(n1))>1.e-20) then
                     !r_fac(n1) = 2*(grad_x(n1)*(xnd(n2)-xnd(n1))+grad_y(n1)*(ynd(n2)-ynd(n1)))/(d_tr0(n2)-d_tr0(n1))-1
                     r_fac(n1) = 2*(grad_x(n1)*(xel(id2,ie)-xel(id,ie))+grad_y(n1)*(yel(id2,ie)-yel(id,ie)))/(d_tr0(n2)-d_tr0(n1))-1
                     tmp1 = d_tr0(n2) - 2*(grad_x(n1)*(xel(id2,ie)-xel(id,ie))+grad_y(n1)*(yel(id2,ie)-yel(id,ie)))
                     !if(tmp1>trcmaxn(n1)) r_fac(n1)=(d_tr0(n1)-trcmaxn(n1))/(d_tr0(n2)-d_tr0(n1))
                     !if(tmp1<trcminn(n1)) r_fac(n1)=(d_tr0(n1)-trcminn(n1))/(d_tr0(n2)-d_tr0(n1))
                     if(tmp1>1.0d0) r_fac(n1)=(d_tr0(n1)-1.0d0)/(d_tr0(n2)-d_tr0(n1))
                     if(tmp1<0.0d0) r_fac(n1)=(d_tr0(n1)-0.0d0)/(d_tr0(n2)-d_tr0(n1))
                     !tmp1 = 2*(grad_x(n1)*(xnd(n2)-xnd(n1))+grad_y(n1)*(ynd(n2)-ynd(n1)))
                     !write(12,*) 'case3',sumtmp,tmp1,d_tr0(n2),d_tr0(n1),r_fac(n1)
                  else
                     r_fac(n1) = -1
                  endif
                  psi_r = (r_fac(n1)+abs(r_fac(n1)))/(1.0d0+abs(r_fac(n1)))
                  !psi_r = max(0.d0,min(1.d0,2.0d0*r_fac(n1)),min(2.d0,r_fac(n1)))
                  !psi_r = max(0.d0,min(1.d0,r_fac(n1)))
                  flux_matrix3(indx1,n1,trc_n) = dt_dyn*vnor1 * (d_tr0(n1)+0.5*psi_r*(d_tr0(n2)-d_tr0(n1)))
                  flux_matrix4(indx1,n1,trc_n) = dt_dyn*vnor2 * (d_tr0(n1)+0.5*psi_r*(d_tr0(n2)-d_tr0(n1)))
                  flux_matrix1(indx2,n2,trc_n) = -1*dt_dyn*vnor1 * (d_tr0(n1)+0.5*psi_r*(d_tr0(n2)-d_tr0(n1)))
                  flux_matrix2(indx2,n2,trc_n) = -1*dt_dyn*vnor2 * (d_tr0(n1)+0.5*psi_r*(d_tr0(n2)-d_tr0(n1)))
                  sumtmp = sumtmp * (d_tr0(n1)+0.5*psi_r*(d_tr0(n2)-d_tr0(n1)))
                  !adv_threshold = exp(-20*tmp2)
                  if((sumtmp/areatmp)>(d_tr0(n1)*adv_threshold))then
                  !if((sumtmp/areatmp)>(d_tr0(n1)*adv_threshold0(n1)).or.(sumtmp/areatmp)>(0.01))then
                  !if((sumtmp/areatmp)>(d_tr0(n1)*adv_threshold).or.(sumtmp/areatmp)>(d_tr0(n1)*adv_threshold0(n2)))then
                     tmp1 = areatmp*d_tr0(n1)/sumtmp*adv_threshold
                     !tmp2 = min(d_tr0(n1)*adv_threshold,d_tr0(n1)*adv_threshold0(n2))
                     !tmp1 = areatmp*tmp2/sumtmp
                     flux_matrix3(indx1,n1,trc_n) = flux_matrix3(indx1,n1,trc_n)*tmp1
                     flux_matrix4(indx1,n1,trc_n) = flux_matrix4(indx1,n1,trc_n)*tmp1
                     flux_matrix1(indx2,n2,trc_n) = flux_matrix1(indx2,n2,trc_n)*tmp1
                     flux_matrix2(indx2,n2,trc_n) = flux_matrix2(indx2,n2,trc_n)*tmp1
                  endif
                  
               else
                  tmp3 = d_tr0(n1)-d_tr0(n2)
                  tmp2=max(c0,tmp3)
                  !sumtmp = sumtmp * d_tr0(n2)
                  flux_matrix00(indx1,n1,2) = n2
                  flux_matrix00(indx2,n2,1) = n2
                  ! n2 center, n1 dowmstream
                  if(abs(d_tr0(n1)-d_tr0(n2))>1.e-20) then
                     !r_fac(n2) = 2*(grad_x(n2)*(xnd(n1)-xnd(n2))+grad_y(n2)*(ynd(n1)-ynd(n2)))/(d_tr0(n1)-d_tr0(n2))-1
                     r_fac(n2) = 2*(grad_x(n2)*(xel(id,ie)-xel(id2,ie))+grad_y(n2)*(yel(id,ie)-yel(id2,ie)))/(d_tr0(n1)-d_tr0(n2))-1
                     tmp1 = d_tr0(n1) - 2*(grad_x(n2)*(xel(id,ie)-xel(id2,ie))+grad_y(n2)*(yel(id,ie)-yel(id2,ie)))
                     !if(tmp1>trcmaxn(n2)) r_fac(n2)=(d_tr0(n2)-trcmaxn(n2))/(d_tr0(n1)-d_tr0(n2))
                     !if(tmp1<trcminn(n2)) r_fac(n2)=(d_tr0(n2)-trcminn(n2))/(d_tr0(n1)-d_tr0(n2))
                     if(tmp1>1.0d0) r_fac(n2)=(d_tr0(n2)-1.0d0)/(d_tr0(n1)-d_tr0(n2))
                     if(tmp1<0.0d0) r_fac(n2)=(d_tr0(n2)-0.0d0)/(d_tr0(n1)-d_tr0(n2))
                     !tmp1 = 2*(grad_x(n2)*(xnd(n1)-xnd(n2))+grad_y(n2)*(ynd(n1)-ynd(n2)))
                     !write(12,*) 'case4',sumtmp,tmp1,d_tr0(n1),d_tr0(n2),r_fac(n2)
                  else
                     r_fac(n2) = -1
                  endif
                  psi_r = (r_fac(n2)+abs(r_fac(n2)))/(1.0d0+abs(r_fac(n2)))
                  !psi_r = max(0.d0,min(1.d0,2.0d0*r_fac(n2)),min(2.d0,r_fac(n2)))
                  !psi_r = max(0.d0,min(1.d0,r_fac(n2)))
                  flux_matrix3(indx1,n1,trc_n) = dt_dyn * vnor1 * (d_tr0(n2)+0.5*psi_r*(d_tr0(n1)-d_tr0(n2)))
                  flux_matrix4(indx1,n1,trc_n) = dt_dyn * vnor2 * (d_tr0(n2)+0.5*psi_r*(d_tr0(n1)-d_tr0(n2)))
                  flux_matrix1(indx2,n2,trc_n) = -1 * dt_dyn * vnor1 * (d_tr0(n2)+0.5*psi_r*(d_tr0(n1)-d_tr0(n2)))
                  flux_matrix2(indx2,n2,trc_n) = -1 * dt_dyn * vnor2 * (d_tr0(n2)+0.5*psi_r*(d_tr0(n1)-d_tr0(n2)))
                  sumtmp = sumtmp * (d_tr0(n2)+0.5*psi_r*(d_tr0(n1)-d_tr0(n2)))
                  !adv_threshold = exp(-20*tmp2)
                  if((abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold))then
                  !if((abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold0(n2)).or.(abs(sumtmp)/areatmp)>0.01)then
                  !if((abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold).or.(abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold0(n1)))then
                     tmp1 = areatmp*d_tr0(n2)/abs(sumtmp)*adv_threshold
                     !tmp2 = min(d_tr0(n2)*adv_threshold,d_tr0(n2)*adv_threshold0(n1))
                     !tmp1 = areatmp*tmp2/abs(sumtmp)
                     flux_matrix3(indx1,n1,trc_n) = flux_matrix3(indx1,n1,trc_n)*tmp1
                     flux_matrix4(indx1,n1,trc_n) = flux_matrix4(indx1,n1,trc_n)*tmp1
                     flux_matrix1(indx2,n2,trc_n) = flux_matrix1(indx2,n2,trc_n)*tmp1
                     flux_matrix2(indx2,n2,trc_n) = flux_matrix2(indx2,n2,trc_n)*tmp1
                  endif
               endif
            endif
          enddo
    enddo
   do i=1,nx_nh
      suma=0.d0 !sum of areas
      fluxchan0=0.d0
      do j=1,nne(i)
         ie=indel(j,i)
         suma=suma+area(ie)/3 !approx for quad
         fluxchan0=fluxchan0 + (flux_matrix1(j,i,trc_n) + flux_matrix2(j,i,trc_n))/2.d0
         fluxchan0=fluxchan0 + (flux_matrix3(j,i,trc_n) + flux_matrix4(j,i,trc_n))/2.d0

      enddo
         if( suma/=0 ) then
            trl0(i)=fluxchan0/suma !m/s/s
         endif !suma/=0  
         trc(i) = d_tr0(i) - trl0(i) !+ d_tr(i)   
   enddo  
    call exchange_p2d(trc)
  end subroutine tvdu_icepack

  !=======================================================================
subroutine tvdu_icepack_other(trc,trc_o,trc_n,trc_d,trc_d1,trc_d2)
        
        !use mod_mesh     
     
        implicit none
        integer (kind=int_kind), intent(in) :: trc_n,trc_d
        real (kind=dbl_kind), dimension(nx), intent(inout) :: trc   ,trc_o
        real (kind=dbl_kind), dimension(nx), intent(in), optional :: trc_d1,trc_d2
        real   (kind=dbl_kind)            :: fluxchan1,fluxchan2,suma,xctr2,yctr2,vnor1,vnor2,&
                                             & tmpx2,tmpx3,tmpy2,tmpy3,utmp,vtmp,xtmp,ytmp,vnorm,&
                                             & trctmp  , trcmin, trcmax, tmpnode   
        integer (kind=int_kind) :: i,j,ie,id,id2,id3,isd3,isd2,     &
                                   m,n1,n2,jj
        logical (kind=log_kind) :: flag_noadv
         real(rkind) :: swild10(3,2),swild(3),trl0(nx),d_tr0(nx)
        !type(t_mesh),        target,        intent(in)    :: mesh

      trl0(:) = c0
      d_tr0(:) = trc(:)
      flux_matrix1(:,:,trc_n) = c0
      flux_matrix2(:,:,trc_n) = c0
      flux_matrix3(:,:,trc_n) = c0
      flux_matrix4(:,:,trc_n) = c0
      flux_matrix0(:,:,:,trc_n) = c0
        ! Driving sequence
      do i=1,nx_nh
        if(idry(i)==1) cycle
        !if(isbnd(1,i)/=0) cycle
        !Wet node
          fluxchan1=0.d0 !sum of fluxes;\int_\Gamma u*u_n d\Gamma
          fluxchan2=0.d0 !sum of fluxes;
          suma=0.d0 !sum of areas
          do j=1,nne(i)
            ie=indel(j,i)
            if(idry_e(ie)==1) cycle
 
            !Wet elem
            suma=suma+area(ie)/3 !approx for quad
            id=iself(j,i)
            id3=nxq(i34(ie)-1,id,i34(ie)) !adjacent side index
            id2=nxq(i34(ie)-2,id,i34(ie)) !adjacent side index
            isd3=elside(id3,ie)
            isd2=elside(id2,ie)

            
            do m=1,i34(ie) !side
                n1=isidenode(1,elside(m,ie))
                n2=isidenode(2,elside(m,ie))
                tmpx2 = trc_o(n1)
                tmpy2 = trc_o(n2)
                if (present(trc_d1)) then
                  tmpx2 = tmpx2 * trc_d1(n1)
                  tmpy2 = tmpy2 * trc_d1(n2)
                endif
                if (present(trc_d2)) then
                  tmpx2 = tmpx2 * trc_d2(n1)
                  tmpy2 = tmpy2 * trc_d2(n2)
                endif
                swild(m) = (tmpx2+tmpy2)*0.5
            enddo
            trctmp = sum(swild(1:i34(ie)))/3
            trcmin = minval(trc_o(indnd(1:nnp(i),i)))
            trcmin = min(trcmin,trc_o(i))
            trcmax = maxval(trc_o(indnd(1:nnp(i),i)))
            trcmax = max(trcmax,trc_o(i))
            
            tmpnode = flux_matrix00(j,i,1)
            tmpx2 = trc_o(tmpnode)
            if (present(trc_d1)) then
               tmpx2 = tmpx2 * trc_d1(tmpnode)
            endif
            if (present(trc_d2)) then
               tmpx2 = tmpx2 * trc_d2(tmpnode)
            endif
            !1st segment
            flux_matrix1(j,i,trc_n) = flux_matrix1(j,i,trc_d)*tmpx2
            flux_matrix2(j,i,trc_n) = flux_matrix2(j,i,trc_d)*tmpx2
            fluxchan1=fluxchan1+(flux_matrix1(j,i,trc_d)*tmpx2+flux_matrix2(j,i,trc_d)*tmpx2)/2.d0
             
            
            !2nd segment
            !Normal dir x length
!            xtmp=ycj(isd2)-yctr(ie)
!            ytmp=xctr(ie)-xcj(isd2)
            tmpnode = flux_matrix00(j,i,2)
            tmpx2 = trc_o(tmpnode)
            if (present(trc_d1)) then
               tmpx2 = tmpx2 * trc_d1(tmpnode)
            endif
            if (present(trc_d2)) then
               tmpx2 = tmpx2 * trc_d2(tmpnode)
            endif
            flux_matrix3(j,i,trc_n) = flux_matrix3(j,i,trc_d)*tmpx2
            flux_matrix4(j,i,trc_n) = flux_matrix4(j,i,trc_d)*tmpx2
            fluxchan1=fluxchan1+(flux_matrix3(j,i,trc_d)*tmpx2+flux_matrix4(j,i,trc_d)*tmpx2)/2.d0
          enddo !j=1,nne(i)
          
          if(suma/=0) then
            trl0(i)=fluxchan1/suma !m/s/s
          endif !suma/=0
      enddo !i=1,np
      call exchange_p2d(trc)
      call exchange_p2d(d_tr0)
      call exchange_p2d(trl0)
      do i = 1,nx_nh
       trc(i) = d_tr0(i) - trl0(i) !+ d_tr(i)
      enddo

      call exchange_p2d(trc)
    end subroutine tvdu_icepack_other

    !=======================================================================
    subroutine tvdu_icepack_dep(trc,trc_o,trc_n,trc_d)
        
        !use mod_mesh     
     
        implicit none
        integer (kind=int_kind), intent(in) :: trc_n,trc_d
        real (kind=dbl_kind), dimension(nx), intent(inout) :: trc   ,trc_o
        real   (kind=dbl_kind)            :: fluxchan1,fluxchan2,suma,xctr2,yctr2,vnor1,vnor2,&
                                             & tmpx2,tmpx3,tmpy2,tmpy3,utmp,vtmp,xtmp,ytmp,vnorm,&
                                             & trctmp, trcmin, trcmax ,sumtmp,areatmp, tmpnode  
        integer (kind=int_kind) :: i,j,ie,id,id2,id3,isd3,isd2,     &
                                   m,n1,n2,jj,indx,nd,nodeid
        logical (kind=log_kind) :: flag_noadv
         real(rkind) :: swild10(3,2),swild(3),trl0(nx),d_tr0(nx)
        !type(t_mesh),        target,        intent(in)    :: mesh

      trl0(:) = c0
      d_tr0(:) = trc(:)
      flux_matrix1(:,:,trc_n) = c0
      flux_matrix2(:,:,trc_n) = c0
      flux_matrix3(:,:,trc_n) = c0
      flux_matrix4(:,:,trc_n) = c0
      flux_matrix0(:,:,:,trc_n) = c0
        ! Driving sequence
      do i=1,nx_nh
        if(idry(i)==1) cycle
        !Wet node
          fluxchan1=0.d0 !sum of fluxes;\int_\Gamma u*u_n d\Gamma
          fluxchan2=0.d0 !sum of fluxes;
          suma=0.d0 !sum of areas
          do j=1,nne(i)
            ie=indel(j,i)
            if(idry_e(ie)==1) cycle
 
            !Wet elem
            areatmp = area(ie)/6
            suma=suma+area(ie)/3 !approx for quad
            id=iself(j,i)
            id3=nxq(i34(ie)-1,id,i34(ie)) !adjacent side index
            id2=nxq(i34(ie)-2,id,i34(ie)) !adjacent side index
            isd3=elside(id3,ie)
            isd2=elside(id2,ie)

            do m=1,i34(ie) !side
                n1=isidenode(1,elside(m,ie))
                n2=isidenode(2,elside(m,ie))
                tmpx2 = trc_o(n1)
                tmpy2 = trc_o(n2)
                swild(m) = (tmpx2+tmpy2)*0.5
            enddo
            trctmp=sum(swild(1:i34(ie)))/3
            !1st segment
            tmpnode = flux_matrix00(j,i,1)
            tmpx2 = trc_o(tmpnode)
            flux_matrix1(j,i,trc_n) = flux_matrix1(j,i,trc_d)*tmpx2
            flux_matrix2(j,i,trc_n) = flux_matrix2(j,i,trc_d)*tmpx2
            sumtmp = (flux_matrix1(j,i,trc_n)+flux_matrix2(j,i,trc_n))/2.d0
            nodeid = elnode(nxq(1,id,i34(ie)),ie)
            fluxchan1=fluxchan1+(flux_matrix1(j,i,trc_n)+flux_matrix2(j,i,trc_n))/2.d0
             
            
            !2nd segment
            !Normal dir x length
!            xtmp=ycj(isd2)-yctr(ie)
!            ytmp=xctr(ie)-xcj(isd2)
            tmpnode = flux_matrix00(j,i,2)
            tmpx2 = trc_o(tmpnode)
            flux_matrix3(j,i,trc_n) = flux_matrix3(j,i,trc_d)*tmpx2
            flux_matrix4(j,i,trc_n) = flux_matrix4(j,i,trc_d)*tmpx2
            sumtmp = (flux_matrix3(j,i,trc_n) + flux_matrix4(j,i,trc_n))/2.d0
            nodeid = elnode(nxq(i34(ie)-1,id,i34(ie)),ie)

            fluxchan1=fluxchan1+(flux_matrix3(j,i,trc_n) + flux_matrix4(j,i,trc_n))/2.d0
          enddo !j=1,nne(i)
          
          if(suma/=0) then
            trl0(i)=fluxchan1/suma !m/s/s
          endif !suma/=0
      enddo !i=1,np
      call exchange_p2d(trc)
      call exchange_p2d(d_tr0)
      call exchange_p2d(trl0)
      do i = 1,nx_nh
       trc(i) = d_tr0(i) - trl0(i) !+ d_tr(i)


      enddo

      call exchange_p2d(trc)
    end subroutine tvdu_icepack_dep

!=======================================================================
    subroutine tvd_casulli_ini
       use mice_module
       implicit none
      real   (kind=dbl_kind)            :: suma,xctr2,yctr2,vnor1,vnor2,sumtmp,areatmp,&
                                           & tmpx2,tmpx3,tmpy2,tmpy3,utmp,vtmp,xtmp,ytmp,vnorm,&
                                           & trctmp, trcmin, trcmax ,puny ,tmp1 ,fluxchan0, tmp2,&
                                           & tmp3, psi_r, sum1, sum2
      integer (kind=int_kind) :: i,j,ie,id,id2,id3,isd3,isd2,     &
                                 m,n1,n2,jj,indx1,indx2,nd,nodeid
       real(kind=dbl_kind) :: swild10(3,2),swild(3),swild2(3,2),trl0(nx),d_tr0(nx), &
                              & r_fac(nx)


      flux_matrix_casulli(:,:,:) = c0
      flux_matrix00(:,:,:) = c0


        ! Driving sequence
    do ie=1,nx_elem
          !Wet elem

          do m=1,i34(ie) !side
            n1=isidenode(1,elside(m,ie))
            n2=isidenode(2,elside(m,ie))
            if(ics==1) then
               swild10(m,1) = 0.5*(uvel(n1)+uvel(n2))
               swild10(m,2) = 0.5*(vvel(n1)+vvel(n2))
               swild2(m,1) = uvel(elnode(m,ie))
               swild2(m,2) = vvel(elnode(m,ie))
            else
               swild10(m,1) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),eframe(:,1,ie))+vvel(n1)*dot_product(pframe(:,2,n1),eframe(:,1,ie))&
                                 &+uvel(n2)*dot_product(pframe(:,1,n2),eframe(:,1,ie))+vvel(n2)*dot_product(pframe(:,2,n2),eframe(:,1,ie)))
               swild10(m,2) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),eframe(:,2,ie))+vvel(n1)*dot_product(pframe(:,2,n1),eframe(:,2,ie))&
                                 &+uvel(n2)*dot_product(pframe(:,1,n2),eframe(:,2,ie))+vvel(n2)*dot_product(pframe(:,2,n2),eframe(:,2,ie)))                      
               swild2(m,1) = uvel(elnode(m,ie))*dot_product(pframe(:,1,elnode(m,ie)),eframe(:,1,ie))+vvel(elnode(m,ie))*dot_product(pframe(:,2,elnode(m,ie)),eframe(:,1,ie))
               swild2(m,2) = uvel(elnode(m,ie))*dot_product(pframe(:,1,elnode(m,ie)),eframe(:,2,ie))+vvel(elnode(m,ie))*dot_product(pframe(:,2,elnode(m,ie)),eframe(:,2,ie))
            endif
            !swild10(m,1) = 0.5*(uvel(n1)*trc(n1)*dot_product(pframe(:,1,n1),sframe2(:,1,elside(m,ie)))+vvel(n1)*trc(n1)*dot_product(pframe(:,2,n1),sframe2(:,1,elside(m,ie)))&
            !                 &+uvel(n2)*trc(n2)*dot_product(pframe(:,1,n2),sframe2(:,1,elside(m,ie)))+vvel(n2)*trc(n2)*dot_product(pframe(:,2,n1),sframe2(:,1,elside(m,ie))))
            !swild10(m,2) = 0.5*(uvel(n1)*trc(n1)*dot_product(pframe(:,1,n1),sframe2(:,2,elside(m,ie)))+vvel(n1)*trc(n1)*dot_product(pframe(:,2,n1),sframe2(:,2,elside(m,ie)))&
            !                 &+uvel(n2)*trc(n2)*dot_product(pframe(:,1,n2),sframe2(:,2,elside(m,ie)))+vvel(n2)*trc(n2)*dot_product(pframe(:,2,n1),sframe2(:,2,elside(m,ie))))
            !swild10(m,1) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),sframe2(:,1,elside(m,ie)))+vvel(n1)*dot_product(pframe(:,2,n1),sframe2(:,1,elside(m,ie)))&
            !                 &+uvel(n2)*dot_product(pframe(:,1,n2),sframe2(:,1,elside(m,ie)))+vvel(n2)*dot_product(pframe(:,2,n2),sframe2(:,1,elside(m,ie))))
            !swild10(m,2) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),sframe2(:,2,elside(m,ie)))+vvel(n1)*dot_product(pframe(:,2,n1),sframe2(:,2,elside(m,ie)))&
            !                 &+uvel(n2)*dot_product(pframe(:,1,n2),sframe2(:,2,elside(m,ie)))+vvel(n2)*dot_product(pframe(:,2,n2),sframe2(:,2,elside(m,ie))))
            
         enddo
        !utmp = sum(swild10(1:i34(ie),1))/3 !vel @ centroid
        !vtmp = sum(swild10(1:i34(ie),2))/3 !vel @ centroid
        utmp=sum(swild2(1:i34(ie),1))/3 !vel @ centroid
        vtmp=sum(swild2(1:i34(ie),2))/3 !vel @ centroid

          do j =1,i34(ie) !side
            n1=isidenode(1,elside(j,ie))
            n2=isidenode(2,elside(j,ie))
            
            indx1 = 0
            do m=1,nne(n1)
               if(indel(m,n1)==ie) then
                  indx1=m; exit
               endif
            enddo
            if(indx1==0) call parallel_abort('STEP: failed to find (tvd_casulli_ini1)') 
            indx2 = 0
            do m=1,nne(n2)
               if(indel(m,n2)==ie) then
                  indx2=m; exit
               endif
            enddo
            if(indx2==0) call parallel_abort('STEP: failed to find (tvd_casulli_ini2)') 

            id  = iself(indx1,n1)
            id2 = iself(indx2,n2)
            id3 = nxq(i34(ie)-2,id,i34(ie))

            if(ics==1) then
               xctr2=xctr(ie); yctr2=yctr(ie)
               tmpx3=(xnd(n1)+xnd(n2))/2.d0
               tmpy3=(ynd(n1)+ynd(n2))/2.d0
            else!ll; use [xy]el defined in eframe
               xctr2=0.d0 !sum(xel(elnode(1:i34(ie),ie)))/i34(ie)
               yctr2=0.d0 !sum(yel(elnode(1:i34(ie),ie)))/i34(ie)
               tmpx3=(xel(id,ie)+xel(id2,ie))/2.d0
               tmpy3=(yel(id,ie)+yel(id2,ie))/2.d0
            endif!ics
            !xtmp=yctr2-tmpy3
            !ytmp=tmpx3-xctr2
            !vnor1=utmp*xtmp+vtmp*ytmp !normal vel x length 
            !vnor2=swild10(j,1)*xtmp+swild10(j,2)*ytmp !normal vel@side x length 
            if(id3.eq.id2) then
               xtmp=yctr2-tmpy3
               ytmp=tmpx3-xctr2
               vnor1=utmp*xtmp+vtmp*ytmp !normal vel x length 
               vnor2=swild10(j,1)*xtmp+swild10(j,2)*ytmp !normal vel@side x length 
               !vnor1=vnor2

               sumtmp = (dt_dyn*vnor1+dt_dyn*vnor2)/2
               if(sumtmp>0) then

                  !sumtmp = sumtmp * d_tr0(n1)
                  flux_matrix00(indx1,n1,1) = n1
                  flux_matrix00(indx2,n2,2) = n1
                  ! n1 center, n2 dowmstream
                  ! positive is outflow
                  flux_matrix_casulli(indx1,n1,1) = dt_dyn * (vnor1 + vnor2)/2
                  flux_matrix_casulli(indx2,n2,2) = -1 * dt_dyn * (vnor1 + vnor2)/2
           
                  
               else

                  !sumtmp = sumtmp * d_tr0(n2)
                  flux_matrix00(indx1,n1,1) = n2
                  flux_matrix00(indx2,n2,2) = n2
                  ! n2 center, n1 dowmstream
                  flux_matrix_casulli(indx1,n1,1) = dt_dyn * (vnor1 + vnor2)/2
                  flux_matrix_casulli(indx2,n2,2) = -1 * dt_dyn * (vnor1 + vnor2)/2
               endif
            else
               xtmp = tmpy3 - yctr2
               ytmp = xctr2 - tmpx3
               vnor1=utmp*xtmp+vtmp*ytmp !normal vel x length 
               vnor2=swild10(j,1)*xtmp+swild10(j,2)*ytmp !normal vel@side x length 
               !vnor1=vnor2

               sumtmp = (dt_dyn*vnor1+dt_dyn*vnor2)/2
               if(sumtmp>0) then

                  !sumtmp = sumtmp * d_tr0(n1)
                  flux_matrix00(indx1,n1,2) = n1
                  flux_matrix00(indx2,n2,1) = n1
                  ! n1 center, n2 dowmstream
                  flux_matrix_casulli(indx1,n1,2) = dt_dyn * (vnor1 + vnor2)/2
                  flux_matrix_casulli(indx2,n2,1) = -1 * dt_dyn * (vnor1 + vnor2)/2
               else

                  !sumtmp = sumtmp * d_tr0(n2)
                  flux_matrix00(indx1,n1,2) = n2
                  flux_matrix00(indx2,n2,1) = n2
                  ! n2 center, n1 dowmstream
                  flux_matrix_casulli(indx1,n1,2) = dt_dyn * (vnor1 + vnor2)/2
                  flux_matrix_casulli(indx2,n2,1) = -1 * dt_dyn * (vnor1 + vnor2)/2
               endif
            endif
          enddo
    enddo


    end subroutine tvd_casulli_ini
!=======================================================================
    subroutine tvd_casulli_icepack(trc,trc_n,adv_threshold0)
        
      !use mod_mesh     
       use mice_module
       use schism_glbl, only: dldxy
      implicit none
      integer (kind=int_kind), intent(in) :: trc_n
      real(kind=dbl_kind), dimension(nx), intent(inout) :: trc,adv_threshold0  
      real   (kind=dbl_kind)            :: suma,xctr2,yctr2,vnor1,vnor2,sumtmp,areatmp,&
                                           & tmpx2,tmpx3,tmpy2,tmpy3,utmp,vtmp,xtmp,ytmp,vnorm,&
                                           & trctmp, trcmin, trcmax ,puny ,tmp1 ,fluxchan0, tmp2, adv_threshold,&
                                           & tmp3, psi_r, sum1, sum2
      integer (kind=int_kind) :: i,j,ie,id,id2,id3,isd3,isd2,     &
                                 m,n1,n2,jj,indx1,indx2,nd,nodeid
      logical (kind=log_kind) :: flag_noadv
       real(kind=dbl_kind) :: swild10(3,2),swild(3),swild2(3,2),trl0(nx),d_tr0(nx),fluxchan1(nx),fluxchan2(nx), &
                              & r_fac(nx),grad_x(nx),grad_y(nx),trcmaxn(nx),trcminn(nx)  
      !type(t_mesh),        target,        intent(in)    :: mesh
       call icepack_query_parameters(puny_out=puny)

    adv_threshold = 10
    trl0(:) = c0
    fluxchan1(:) = c0 !positive flux
    fluxchan2(:) = c0 !negative flux
    d_tr0(:) = trc(:)
    flux_matrix1(:,:,trc_n) = c0
    flux_matrix2(:,:,trc_n) = c0
    flux_matrix3(:,:,trc_n) = c0
    flux_matrix4(:,:,trc_n) = c0
    !flux_matrix00(:,:,:) = c0
    flux_matrix0(:,:,:,trc_n) = c0
    r_fac(:) = c0

    flux_ratio_casulli(:) = c0

     do i=1,nx_nh
      tmp1=0.d0 !sum of flux
      fluxchan0=0.d0
      
      do j=1,nne(i)
         ie=indel(j,i)
         if(flux_matrix_casulli(j,i,1) < 0 ) then
            fluxchan0 = fluxchan0 + abs(flux_matrix_casulli(j,i,1)) * (trc(i)-trc(flux_matrix00(j,i,1)))
            tmp1 = tmp1 + abs(flux_matrix_casulli(j,i,1))
         endif
         
         if(flux_matrix_casulli(j,i,2) < 0 ) then
            fluxchan0 = fluxchan0 + abs(flux_matrix_casulli(j,i,2)) * (trc(i)-trc(flux_matrix00(j,i,2)))
            tmp1 = tmp1 + abs(flux_matrix_casulli(j,i,2))
         endif
      enddo
      if(tmp1 > 0) flux_ratio_casulli(i) = fluxchan0 / tmp1
      !if(trc(i) > 0) write(12,*) i,trc(i),flux_ratio_casulli(i),tmp1,flux_matrix_casulli(:,i,2)
   enddo  

   call exchange_p2d(flux_ratio_casulli)





      ! Driving sequence
    do ie=1,nx_elem

          if(idry_e(ie)==1) cycle
          areatmp = area(ie)/6
          !Wet elem

          do m=1,i34(ie) !side
            n1=isidenode(1,elside(m,ie))
            n2=isidenode(2,elside(m,ie))
            if(ics==1) then
               swild10(m,1) = 0.5*(uvel(n1)+uvel(n2))
               swild10(m,2) = 0.5*(vvel(n1)+vvel(n2))
               swild(m) = (trc(n1)+trc(n2))*0.5
               swild2(m,1) = uvel(elnode(m,ie))
               swild2(m,2) = vvel(elnode(m,ie))
            else
               swild10(m,1) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),eframe(:,1,ie))+vvel(n1)*dot_product(pframe(:,2,n1),eframe(:,1,ie))&
                                 &+uvel(n2)*dot_product(pframe(:,1,n2),eframe(:,1,ie))+vvel(n2)*dot_product(pframe(:,2,n2),eframe(:,1,ie)))
               swild10(m,2) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),eframe(:,2,ie))+vvel(n1)*dot_product(pframe(:,2,n1),eframe(:,2,ie))&
                                 &+uvel(n2)*dot_product(pframe(:,1,n2),eframe(:,2,ie))+vvel(n2)*dot_product(pframe(:,2,n2),eframe(:,2,ie)))                      
               swild(m) = (trc(n1)+trc(n2))*0.5
               swild2(m,1) = uvel(elnode(m,ie))*dot_product(pframe(:,1,elnode(m,ie)),eframe(:,1,ie))+vvel(elnode(m,ie))*dot_product(pframe(:,2,elnode(m,ie)),eframe(:,1,ie))
               swild2(m,2) = uvel(elnode(m,ie))*dot_product(pframe(:,1,elnode(m,ie)),eframe(:,2,ie))+vvel(elnode(m,ie))*dot_product(pframe(:,2,elnode(m,ie)),eframe(:,2,ie))
            endif
            !swild10(m,1) = 0.5*(uvel(n1)*trc(n1)*dot_product(pframe(:,1,n1),sframe2(:,1,elside(m,ie)))+vvel(n1)*trc(n1)*dot_product(pframe(:,2,n1),sframe2(:,1,elside(m,ie)))&
            !                 &+uvel(n2)*trc(n2)*dot_product(pframe(:,1,n2),sframe2(:,1,elside(m,ie)))+vvel(n2)*trc(n2)*dot_product(pframe(:,2,n1),sframe2(:,1,elside(m,ie))))
            !swild10(m,2) = 0.5*(uvel(n1)*trc(n1)*dot_product(pframe(:,1,n1),sframe2(:,2,elside(m,ie)))+vvel(n1)*trc(n1)*dot_product(pframe(:,2,n1),sframe2(:,2,elside(m,ie)))&
            !                 &+uvel(n2)*trc(n2)*dot_product(pframe(:,1,n2),sframe2(:,2,elside(m,ie)))+vvel(n2)*trc(n2)*dot_product(pframe(:,2,n1),sframe2(:,2,elside(m,ie))))
            !swild10(m,1) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),sframe2(:,1,elside(m,ie)))+vvel(n1)*dot_product(pframe(:,2,n1),sframe2(:,1,elside(m,ie)))&
            !                 &+uvel(n2)*dot_product(pframe(:,1,n2),sframe2(:,1,elside(m,ie)))+vvel(n2)*dot_product(pframe(:,2,n2),sframe2(:,1,elside(m,ie))))
            !swild10(m,2) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),sframe2(:,2,elside(m,ie)))+vvel(n1)*dot_product(pframe(:,2,n1),sframe2(:,2,elside(m,ie)))&
            !                 &+uvel(n2)*dot_product(pframe(:,1,n2),sframe2(:,2,elside(m,ie)))+vvel(n2)*dot_product(pframe(:,2,n2),sframe2(:,2,elside(m,ie))))
            
         enddo
        !utmp = sum(swild10(1:i34(ie),1))/3 !vel @ centroid
        !vtmp = sum(swild10(1:i34(ie),2))/3 !vel @ centroid
        utmp=sum(swild2(1:i34(ie),1))/3 !vel @ centroid
        vtmp=sum(swild2(1:i34(ie),2))/3 !vel @ centroid
        trctmp = sum(swild(1:i34(ie)))/3

          do j =1,i34(ie) !side
            n1=isidenode(1,elside(j,ie))
            n2=isidenode(2,elside(j,ie))
            
            indx1 = 0
            do m=1,nne(n1)
               if(indel(m,n1)==ie) then
                  indx1=m; exit
               endif
            enddo
            if(indx1==0) call parallel_abort('STEP: failed to find (tvd_casulli_icepack_1)') 
            indx2 = 0
            do m=1,nne(n2)
               if(indel(m,n2)==ie) then
                  indx2=m; exit
               endif
            enddo
            if(indx2==0) call parallel_abort('STEP: failed to find (tvd_casulli_icepack_2)') 

            id  = iself(indx1,n1)
            id2 = iself(indx2,n2)
            id3 = nxq(i34(ie)-2,id,i34(ie))

            if(ics==1) then
               xctr2=xctr(ie); yctr2=yctr(ie)
               tmpx3=(xnd(n1)+xnd(n2))/2.d0
               tmpy3=(ynd(n1)+ynd(n2))/2.d0
            else!ll; use [xy]el defined in eframe
               xctr2=0.d0 !sum(xel(elnode(1:i34(ie),ie)))/i34(ie)
               yctr2=0.d0 !sum(yel(elnode(1:i34(ie),ie)))/i34(ie)
               tmpx3=(xel(id,ie)+xel(id2,ie))/2.d0
               tmpy3=(yel(id,ie)+yel(id2,ie))/2.d0
            endif!ics
            !xtmp=yctr2-tmpy3
            !ytmp=tmpx3-xctr2
            !vnor1=utmp*xtmp+vtmp*ytmp !normal vel x length 
            !vnor2=swild10(j,1)*xtmp+swild10(j,2)*ytmp !normal vel@side x length 
            if(id3.eq.id2) then
               xtmp=yctr2-tmpy3
               ytmp=tmpx3-xctr2
               vnor1=utmp*xtmp+vtmp*ytmp !normal vel x length 
               vnor2=swild10(j,1)*xtmp+swild10(j,2)*ytmp !normal vel@side x length 
               !vnor1=vnor2

               sumtmp = (dt_dyn*vnor1+dt_dyn*vnor2)/2
               if(sumtmp>0) then
                  tmp3 = d_tr0(n2)-d_tr0(n1)
                  tmp2=max(c0,tmp3)
                  !sumtmp = sumtmp * d_tr0(n1)
                  !flux_matrix00(indx1,n1,1) = n1
                  !flux_matrix00(indx2,n2,2) = n1
                  ! n1 center, n2 dowmstream
                  if(abs(d_tr0(n2)-d_tr0(n1))>1.e-20) then
                     !r_fac(n1) = 2*(grad_x(n1)*(xnd(n2)-xnd(n1))+grad_y(n1)*(ynd(n2)-ynd(n1)))/(d_tr0(n2)-d_tr0(n1))-1
                     r_fac(n1) = flux_ratio_casulli(n1)/(d_tr0(n2)-d_tr0(n1))
                     !tmp1 = d_tr0(n2) - 2*(grad_x(n1)*(xel(id2,ie)-xel(id,ie))+grad_y(n1)*(yel(id2,ie)-yel(id,ie)))
                     !if(tmp1>trcmaxn(n1)) r_fac(n1)=(d_tr0(n1)-trcmaxn(n1))/(d_tr0(n2)-d_tr0(n1))  
                     !if(tmp1<trcminn(n1)) r_fac(n1)=(d_tr0(n1)-trcminn(n1))/(d_tr0(n2)-d_tr0(n1))   
                     !if(tmp1>1.0d0) r_fac(n1)=(d_tr0(n1)-1.0d0)/(d_tr0(n2)-d_tr0(n1))  
                     !if(tmp1<0.0d0) r_fac(n1)=(d_tr0(n1)-0.0d0)/(d_tr0(n2)-d_tr0(n1))             
                     !tmp1 = 2*(grad_x(n1)*(xnd(n2)-xnd(n1))+grad_y(n1)*(ynd(n2)-ynd(n1)))
                     !write(12,*) 'case1',xnd(n1)-xnd(n2),ynd(n2)-ynd(n1),znd(n1)-znd(n2),xel(id,ie)-xel(id2,ie),yel(id,ie)-yel(id2,ie)
                     !write(12,*) 'case1',sumtmp,tmp1,d_tr0(n2),d_tr0(n1),r_fac(n1)
                  else
                     r_fac(n1) = -1
                  endif
                  psi_r = (r_fac(n1)+abs(r_fac(n1)))/(1.0d0+abs(r_fac(n1)))
                  !psi_r = max(0.d0,min(1.d0,2.0d0*r_fac(n1)),min(2.d0,r_fac(n1)))
                  !psi_r = max(0.d0,min(1.d0,r_fac(n1)))
                  flux_matrix1(indx1,n1,trc_n) = dt_dyn * ( vnor1 * d_tr0(n1)+0.5*abs(vnor1)*psi_r*(d_tr0(n2)-d_tr0(n1)))
                  flux_matrix2(indx1,n1,trc_n) = dt_dyn * ( vnor2 * d_tr0(n1)+0.5*abs(vnor2)*psi_r*(d_tr0(n2)-d_tr0(n1)))
                  flux_matrix3(indx2,n2,trc_n) = -1 * dt_dyn * ( vnor1 * d_tr0(n1)+0.5*abs(vnor1)*psi_r*(d_tr0(n2)-d_tr0(n1)))
                  flux_matrix4(indx2,n2,trc_n) = -1 * dt_dyn * ( vnor2 * d_tr0(n1)+0.5*abs(vnor2)*psi_r*(d_tr0(n2)-d_tr0(n1)))  
                  sumtmp = sumtmp * (d_tr0(n1)+0.5*psi_r*(d_tr0(n2)-d_tr0(n1)))             
                  !adv_threshold = exp(-20*tmp2)
                  if((sumtmp/areatmp)>(d_tr0(n1)*adv_threshold))then
                  !if((sumtmp/areatmp)>(d_tr0(n1)*adv_threshold0(n1)).or.(sumtmp/areatmp)>(0.01))then
                  !if((sumtmp/areatmp)>(d_tr0(n1)*adv_threshold).or.(sumtmp/areatmp)>(d_tr0(n1)*adv_threshold0(n2)))then
                     tmp1 = areatmp*d_tr0(n1)/sumtmp*adv_threshold
                     !tmp2 = min(d_tr0(n1)*adv_threshold,d_tr0(n1)*adv_threshold0(n2))
                     !tmp1 = areatmp*tmp2/sumtmp
                     !flux_matrix1(indx1,n1,trc_n) = flux_matrix1(indx1,n1,trc_n)*tmp1
                     !flux_matrix2(indx1,n1,trc_n) = flux_matrix2(indx1,n1,trc_n)*tmp1
                     !flux_matrix3(indx2,n2,trc_n) = flux_matrix3(indx2,n2,trc_n)*tmp1
                     !flux_matrix4(indx2,n2,trc_n) = flux_matrix4(indx2,n2,trc_n)*tmp1
                  endif
                  
               else
                  tmp3 = d_tr0(n1)-d_tr0(n2)
                  tmp2=max(c0,tmp3)
                  !sumtmp = sumtmp * d_tr0(n2)
                  !flux_matrix00(indx1,n1,1) = n2
                  !flux_matrix00(indx2,n2,2) = n2
                  ! n2 center, n1 dowmstream
                  if(abs(d_tr0(n2)-d_tr0(n1))>1.e-20) then
                     !r_fac(n2) = 2*(grad_x(n2)*(xnd(n1)-xnd(n2))+grad_y(n2)*(ynd(n1)-ynd(n2)))/(d_tr0(n1)-d_tr0(n2))-1
                     r_fac(n2) = flux_ratio_casulli(n2)/(d_tr0(n1)-d_tr0(n2))
                     !tmp1 = d_tr0(n1) - 2*(grad_x(n2)*(xel(id,ie)-xel(id2,ie))+grad_y(n2)*(yel(id,ie)-yel(id2,ie)))
                     !if(tmp1>trcmaxn(n2)) r_fac(n2)=(d_tr0(n2)-trcmaxn(n2))/(d_tr0(n1)-d_tr0(n2))
                     !if(tmp1<trcminn(n2)) r_fac(n2)=(d_tr0(n2)-trcminn(n2))/(d_tr0(n1)-d_tr0(n2))
                     !if(tmp1>1.0d0) r_fac(n2)=(d_tr0(n2)-1.0d0)/(d_tr0(n1)-d_tr0(n2))
                     !if(tmp1<0.0d0) r_fac(n2)=(d_tr0(n2)-0.0d0)/(d_tr0(n1)-d_tr0(n2))
                     !tmp1 = 2*(grad_x(n2)*(xnd(n1)-xnd(n2))+grad_y(n2)*(ynd(n1)-ynd(n2)))
                     !write(12,*) 'case2',sumtmp,tmp1,d_tr0(n1),d_tr0(n2),r_fac(n2)
                  else
                     r_fac(n2) = -1
                  endif
                  psi_r = (r_fac(n2)+abs(r_fac(n2)))/(1.0d0+abs(r_fac(n2)))
                  !psi_r = max(0.d0,min(1.d0,2.0d0*r_fac(n2)),min(2.d0,r_fac(n2)))
                  !psi_r = max(0.d0,min(1.d0,r_fac(n2)))
                  flux_matrix1(indx1,n1,trc_n) = dt_dyn * ( abs(vnor1) * -1 * d_tr0(n2)+0.5*abs(vnor1)*psi_r*(d_tr0(n2)-d_tr0(n1)))
                  flux_matrix2(indx1,n1,trc_n) = dt_dyn * ( abs(vnor2) * -1 * d_tr0(n2)+0.5*abs(vnor2)*psi_r*(d_tr0(n2)-d_tr0(n1)))
                  flux_matrix3(indx2,n2,trc_n) = dt_dyn * ( abs(vnor1) * d_tr0(n2)+0.5*abs(vnor1)*psi_r*(d_tr0(n1)-d_tr0(n2)))
                  flux_matrix4(indx2,n2,trc_n) = dt_dyn * ( abs(vnor2) * d_tr0(n2)+0.5*abs(vnor2)*psi_r*(d_tr0(n1)-d_tr0(n2)))  
                  sumtmp = sumtmp * (d_tr0(n2)+0.5*psi_r*(d_tr0(n1)-d_tr0(n2)))               
                  !adv_threshold = exp(-20*tmp2)
                  if((abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold))then
                  !if((abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold0(n2)).or.(abs(sumtmp)/areatmp)>(0.01))then
                  !if((abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold).or.(abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold0(n1)))then
                     tmp1 = areatmp*d_tr0(n2)/abs(sumtmp)*adv_threshold
                     !tmp2 = min(d_tr0(n2)*adv_threshold,d_tr0(n2)*adv_threshold0(n1))
                     !tmp1 = areatmp*tmp2/abs(sumtmp)
                     !flux_matrix1(indx1,n1,trc_n) = flux_matrix1(indx1,n1,trc_n)*tmp1
                     !flux_matrix2(indx1,n1,trc_n) = flux_matrix2(indx1,n1,trc_n)*tmp1
                     !flux_matrix3(indx2,n2,trc_n) = flux_matrix3(indx2,n2,trc_n)*tmp1
                     !flux_matrix4(indx2,n2,trc_n) = flux_matrix4(indx2,n2,trc_n)*tmp1
                  endif
               endif
            else
               xtmp = tmpy3 - yctr2
               ytmp = xctr2 - tmpx3
               vnor1=utmp*xtmp+vtmp*ytmp !normal vel x length 
               vnor2=swild10(j,1)*xtmp+swild10(j,2)*ytmp !normal vel@side x length 
               !vnor1=vnor2

               sumtmp = (dt_dyn*vnor1+dt_dyn*vnor2)/2
               if(sumtmp>0) then
                  tmp3 = d_tr0(n2)-d_tr0(n1)
                  tmp2=max(c0,tmp3)
                  !sumtmp = sumtmp * d_tr0(n1)
                  !flux_matrix00(indx1,n1,2) = n1
                  !flux_matrix00(indx2,n2,1) = n1
                  ! n1 center, n2 dowmstream
                  if(abs(d_tr0(n2)-d_tr0(n1))>1.e-20) then
                     !r_fac(n1) = 2*(grad_x(n1)*(xnd(n2)-xnd(n1))+grad_y(n1)*(ynd(n2)-ynd(n1)))/(d_tr0(n2)-d_tr0(n1))-1
                     r_fac(n1) = flux_ratio_casulli(n1)/(d_tr0(n2)-d_tr0(n1))
                     !tmp1 = d_tr0(n2) - 2*(grad_x(n1)*(xel(id2,ie)-xel(id,ie))+grad_y(n1)*(yel(id2,ie)-yel(id,ie)))
                     !if(tmp1>trcmaxn(n1)) r_fac(n1)=(d_tr0(n1)-trcmaxn(n1))/(d_tr0(n2)-d_tr0(n1))
                     !if(tmp1<trcminn(n1)) r_fac(n1)=(d_tr0(n1)-trcminn(n1))/(d_tr0(n2)-d_tr0(n1))
                     !if(tmp1>1.0d0) r_fac(n1)=(d_tr0(n1)-1.0d0)/(d_tr0(n2)-d_tr0(n1))
                     !if(tmp1<0.0d0) r_fac(n1)=(d_tr0(n1)-0.0d0)/(d_tr0(n2)-d_tr0(n1))
                     !tmp1 = 2*(grad_x(n1)*(xnd(n2)-xnd(n1))+grad_y(n1)*(ynd(n2)-ynd(n1)))
                     !write(12,*) 'case3',sumtmp,tmp1,d_tr0(n2),d_tr0(n1),r_fac(n1)
                  else
                     r_fac(n1) = -1
                  endif
                  psi_r = (r_fac(n1)+abs(r_fac(n1)))/(1.0d0+abs(r_fac(n1)))
                  !psi_r = max(0.d0,min(1.d0,2.0d0*r_fac(n1)),min(2.d0,r_fac(n1)))
                  !psi_r = max(0.d0,min(1.d0,r_fac(n1)))
                  flux_matrix3(indx1,n1,trc_n) = dt_dyn * ( vnor1 * d_tr0(n1)+0.5*abs(vnor1)*psi_r*(d_tr0(n2)-d_tr0(n1)))
                  flux_matrix4(indx1,n1,trc_n) = dt_dyn * ( vnor2 * d_tr0(n1)+0.5*abs(vnor2)*psi_r*(d_tr0(n2)-d_tr0(n1)))
                  flux_matrix1(indx2,n2,trc_n) = -1 * dt_dyn * ( vnor1 * d_tr0(n1)+0.5*abs(vnor1)*psi_r*(d_tr0(n2)-d_tr0(n1)))
                  flux_matrix2(indx2,n2,trc_n) = -1 * dt_dyn * ( vnor2 * d_tr0(n1)+0.5*abs(vnor2)*psi_r*(d_tr0(n2)-d_tr0(n1)))
                  sumtmp = sumtmp * (d_tr0(n1)+0.5*psi_r*(d_tr0(n2)-d_tr0(n1)))
                  !adv_threshold = exp(-20*tmp2)
                  if((sumtmp/areatmp)>(d_tr0(n1)*adv_threshold))then
                  !if((sumtmp/areatmp)>(d_tr0(n1)*adv_threshold0(n1)).or.(sumtmp/areatmp)>(0.01))then
                  !if((sumtmp/areatmp)>(d_tr0(n1)*adv_threshold).or.(sumtmp/areatmp)>(d_tr0(n1)*adv_threshold0(n2)))then
                     tmp1 = areatmp*d_tr0(n1)/sumtmp*adv_threshold
                     !tmp2 = min(d_tr0(n1)*adv_threshold,d_tr0(n1)*adv_threshold0(n2))
                     !tmp1 = areatmp*tmp2/sumtmp
                     !flux_matrix3(indx1,n1,trc_n) = flux_matrix3(indx1,n1,trc_n)*tmp1
                     !flux_matrix4(indx1,n1,trc_n) = flux_matrix4(indx1,n1,trc_n)*tmp1
                     !flux_matrix1(indx2,n2,trc_n) = flux_matrix1(indx2,n2,trc_n)*tmp1
                     !flux_matrix2(indx2,n2,trc_n) = flux_matrix2(indx2,n2,trc_n)*tmp1
                  endif
                  
               else
                  tmp3 = d_tr0(n1)-d_tr0(n2)
                  tmp2=max(c0,tmp3)
                  !sumtmp = sumtmp * d_tr0(n2)
                  !flux_matrix00(indx1,n1,2) = n2
                  !flux_matrix00(indx2,n2,1) = n2
                  ! n2 center, n1 dowmstream
                  if(abs(d_tr0(n1)-d_tr0(n2))>1.e-20) then
                     !r_fac(n2) = 2*(grad_x(n2)*(xnd(n1)-xnd(n2))+grad_y(n2)*(ynd(n1)-ynd(n2)))/(d_tr0(n1)-d_tr0(n2))-1
                     r_fac(n2) = flux_ratio_casulli(n2)/(d_tr0(n1)-d_tr0(n2))
                     !tmp1 = d_tr0(n1) - 2*(grad_x(n2)*(xel(id,ie)-xel(id2,ie))+grad_y(n2)*(yel(id,ie)-yel(id2,ie)))
                     !if(tmp1>trcmaxn(n2)) r_fac(n2)=(d_tr0(n2)-trcmaxn(n2))/(d_tr0(n1)-d_tr0(n2))
                     !if(tmp1<trcminn(n2)) r_fac(n2)=(d_tr0(n2)-trcminn(n2))/(d_tr0(n1)-d_tr0(n2))
                     !if(tmp1>1.0d0) r_fac(n2)=(d_tr0(n2)-1.0d0)/(d_tr0(n1)-d_tr0(n2))
                     !if(tmp1<0.0d0) r_fac(n2)=(d_tr0(n2)-0.0d0)/(d_tr0(n1)-d_tr0(n2))
                     !tmp1 = 2*(grad_x(n2)*(xnd(n1)-xnd(n2))+grad_y(n2)*(ynd(n1)-ynd(n2)))
                     !write(12,*) 'case4',sumtmp,tmp1,d_tr0(n1),d_tr0(n2),r_fac(n2)
                  else
                     r_fac(n2) = -1
                  endif
                  psi_r = (r_fac(n2)+abs(r_fac(n2)))/(1.0d0+abs(r_fac(n2)))
                  !psi_r = max(0.d0,min(1.d0,2.0d0*r_fac(n2)),min(2.d0,r_fac(n2)))
                  !psi_r = max(0.d0,min(1.d0,r_fac(n2)))
                  flux_matrix3(indx1,n1,trc_n) = dt_dyn * ( abs(vnor1) * -1 * d_tr0(n2)+0.5*abs(vnor1)*psi_r*(d_tr0(n2)-d_tr0(n1)))
                  flux_matrix4(indx1,n1,trc_n) = dt_dyn * ( abs(vnor2) * -1 * d_tr0(n2)+0.5*abs(vnor2)*psi_r*(d_tr0(n2)-d_tr0(n1)))
                  flux_matrix1(indx2,n2,trc_n) = dt_dyn * ( abs(vnor1) * d_tr0(n2)+0.5*abs(vnor1)*psi_r*(d_tr0(n1)-d_tr0(n2)))
                  flux_matrix2(indx2,n2,trc_n) = dt_dyn * ( abs(vnor2) * d_tr0(n2)+0.5*abs(vnor2)*psi_r*(d_tr0(n1)-d_tr0(n2)))
                  sumtmp = sumtmp * (d_tr0(n2)+0.5*psi_r*(d_tr0(n1)-d_tr0(n2)))
                  !adv_threshold = exp(-20*tmp2)
                  if((abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold))then
                  !if((abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold0(n2)).or.(abs(sumtmp)/areatmp)>0.01)then
                  !if((abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold).or.(abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold0(n1)))then
                     tmp1 = areatmp*d_tr0(n2)/abs(sumtmp)*adv_threshold
                     !tmp2 = min(d_tr0(n2)*adv_threshold,d_tr0(n2)*adv_threshold0(n1))
                     !tmp1 = areatmp*tmp2/abs(sumtmp)
                     !flux_matrix3(indx1,n1,trc_n) = flux_matrix3(indx1,n1,trc_n)*tmp1
                     !flux_matrix4(indx1,n1,trc_n) = flux_matrix4(indx1,n1,trc_n)*tmp1
                     !flux_matrix1(indx2,n2,trc_n) = flux_matrix1(indx2,n2,trc_n)*tmp1
                     !flux_matrix2(indx2,n2,trc_n) = flux_matrix2(indx2,n2,trc_n)*tmp1
                  endif
               endif
            endif
          enddo
    enddo
   do i=1,nx_nh
      suma=0.d0 !sum of areas
      fluxchan0=0.d0
      do j=1,nne(i)
         ie=indel(j,i)
         suma=suma+area(ie)/3 !approx for quad
         fluxchan0=fluxchan0 + (flux_matrix1(j,i,trc_n) + flux_matrix2(j,i,trc_n))/2.d0
         fluxchan0=fluxchan0 + (flux_matrix3(j,i,trc_n) + flux_matrix4(j,i,trc_n))/2.d0

      enddo
         if( suma/=0 ) then
            trl0(i)=fluxchan0/suma !m/s/s
         endif !suma/=0  
         trc(i) = d_tr0(i) - trl0(i) !+ d_tr(i)   
   enddo  
    call exchange_p2d(trc)
  end subroutine tvd_casulli_icepack

  !=======================================================================
subroutine tvd_casulli_icepack_other(trc,trc_o,trc_n,trc_d,trc_d1,trc_d2)
        
        !use mod_mesh     
     
        implicit none
        integer (kind=int_kind), intent(in) :: trc_n,trc_d
        real (kind=dbl_kind), dimension(nx), intent(inout) :: trc   ,trc_o
        real (kind=dbl_kind), dimension(nx), intent(in), optional :: trc_d1,trc_d2
        real   (kind=dbl_kind)            :: fluxchan1,fluxchan2,suma,xctr2,yctr2,vnor1,vnor2,&
                                             & tmpx2,tmpx3,tmpy2,tmpy3,utmp,vtmp,xtmp,ytmp,vnorm,&
                                             & trctmp  , trcmin, trcmax, tmpnode   
        integer (kind=int_kind) :: i,j,ie,id,id2,id3,isd3,isd2,     &
                                   m,n1,n2,jj
        logical (kind=log_kind) :: flag_noadv
         real(rkind) :: swild10(3,2),swild(3),trl0(nx),d_tr0(nx)
        !type(t_mesh),        target,        intent(in)    :: mesh

      trl0(:) = c0
      d_tr0(:) = trc(:)
      flux_matrix1(:,:,trc_n) = c0
      flux_matrix2(:,:,trc_n) = c0
      flux_matrix3(:,:,trc_n) = c0
      flux_matrix4(:,:,trc_n) = c0
      flux_matrix0(:,:,:,trc_n) = c0
        ! Driving sequence
      do i=1,nx_nh
        if(idry(i)==1) cycle
        !if(isbnd(1,i)/=0) cycle
        !Wet node
          fluxchan1=0.d0 !sum of fluxes;\int_\Gamma u*u_n d\Gamma
          fluxchan2=0.d0 !sum of fluxes;
          suma=0.d0 !sum of areas
          do j=1,nne(i)
            ie=indel(j,i)
            if(idry_e(ie)==1) cycle
 
            !Wet elem
            suma=suma+area(ie)/3 !approx for quad
            id=iself(j,i)
            id3=nxq(i34(ie)-1,id,i34(ie)) !adjacent side index
            id2=nxq(i34(ie)-2,id,i34(ie)) !adjacent side index
            isd3=elside(id3,ie)
            isd2=elside(id2,ie)

            
            do m=1,i34(ie) !side
                n1=isidenode(1,elside(m,ie))
                n2=isidenode(2,elside(m,ie))
                tmpx2 = trc_o(n1)
                tmpy2 = trc_o(n2)
                if (present(trc_d1)) then
                  tmpx2 = tmpx2 * trc_d1(n1)
                  tmpy2 = tmpy2 * trc_d1(n2)
                endif
                if (present(trc_d2)) then
                  tmpx2 = tmpx2 * trc_d2(n1)
                  tmpy2 = tmpy2 * trc_d2(n2)
                endif
                swild(m) = (tmpx2+tmpy2)*0.5
            enddo
            trctmp = sum(swild(1:i34(ie)))/3
            trcmin = minval(trc_o(indnd(1:nnp(i),i)))
            trcmin = min(trcmin,trc_o(i))
            trcmax = maxval(trc_o(indnd(1:nnp(i),i)))
            trcmax = max(trcmax,trc_o(i))
            
            tmpnode = flux_matrix00(j,i,1)
            tmpx2 = trc_o(tmpnode)
            if (present(trc_d1)) then
               tmpx2 = tmpx2 * trc_d1(tmpnode)
            endif
            if (present(trc_d2)) then
               tmpx2 = tmpx2 * trc_d2(tmpnode)
            endif
            !1st segment
            flux_matrix1(j,i,trc_n) = flux_matrix1(j,i,trc_d)*tmpx2
            flux_matrix2(j,i,trc_n) = flux_matrix2(j,i,trc_d)*tmpx2
            fluxchan1=fluxchan1+(flux_matrix1(j,i,trc_d)*tmpx2+flux_matrix2(j,i,trc_d)*tmpx2)/2.d0
             
            
            !2nd segment
            !Normal dir x length
!            xtmp=ycj(isd2)-yctr(ie)
!            ytmp=xctr(ie)-xcj(isd2)
            tmpnode = flux_matrix00(j,i,2)
            tmpx2 = trc_o(tmpnode)
            if (present(trc_d1)) then
               tmpx2 = tmpx2 * trc_d1(tmpnode)
            endif
            if (present(trc_d2)) then
               tmpx2 = tmpx2 * trc_d2(tmpnode)
            endif
            flux_matrix3(j,i,trc_n) = flux_matrix3(j,i,trc_d)*tmpx2
            flux_matrix4(j,i,trc_n) = flux_matrix4(j,i,trc_d)*tmpx2
            fluxchan1=fluxchan1+(flux_matrix3(j,i,trc_d)*tmpx2+flux_matrix4(j,i,trc_d)*tmpx2)/2.d0
          enddo !j=1,nne(i)
          
          if(suma/=0) then
            trl0(i)=fluxchan1/suma !m/s/s
          endif !suma/=0
      enddo !i=1,np
      call exchange_p2d(trc)
      call exchange_p2d(d_tr0)
      call exchange_p2d(trl0)
      do i = 1,nx_nh
       trc(i) = d_tr0(i) - trl0(i) !+ d_tr(i)
      enddo

      call exchange_p2d(trc)
    end subroutine tvd_casulli_icepack_other

    !=======================================================================
    subroutine tvd_casulli_icepack_dep(trc,trc_o,trc_n,trc_d)
        
        !use mod_mesh     
     
        implicit none
        integer (kind=int_kind), intent(in) :: trc_n,trc_d
        real (kind=dbl_kind), dimension(nx), intent(inout) :: trc   ,trc_o
        real   (kind=dbl_kind)            :: fluxchan1,fluxchan2,suma,xctr2,yctr2,vnor1,vnor2,&
                                             & tmpx2,tmpx3,tmpy2,tmpy3,utmp,vtmp,xtmp,ytmp,vnorm,&
                                             & trctmp, trcmin, trcmax ,sumtmp,areatmp, tmpnode  
        integer (kind=int_kind) :: i,j,ie,id,id2,id3,isd3,isd2,     &
                                   m,n1,n2,jj,indx,nd,nodeid
        logical (kind=log_kind) :: flag_noadv
         real(rkind) :: swild10(3,2),swild(3),trl0(nx),d_tr0(nx)
        !type(t_mesh),        target,        intent(in)    :: mesh

      trl0(:) = c0
      d_tr0(:) = trc(:)
      flux_matrix1(:,:,trc_n) = c0
      flux_matrix2(:,:,trc_n) = c0
      flux_matrix3(:,:,trc_n) = c0
      flux_matrix4(:,:,trc_n) = c0
      flux_matrix0(:,:,:,trc_n) = c0
        ! Driving sequence
      do i=1,nx_nh
        if(idry(i)==1) cycle
        !Wet node
          fluxchan1=0.d0 !sum of fluxes;\int_\Gamma u*u_n d\Gamma
          fluxchan2=0.d0 !sum of fluxes;
          suma=0.d0 !sum of areas
          do j=1,nne(i)
            ie=indel(j,i)
            if(idry_e(ie)==1) cycle
 
            !Wet elem
            areatmp = area(ie)/6
            suma=suma+area(ie)/3 !approx for quad
            id=iself(j,i)
            id3=nxq(i34(ie)-1,id,i34(ie)) !adjacent side index
            id2=nxq(i34(ie)-2,id,i34(ie)) !adjacent side index
            isd3=elside(id3,ie)
            isd2=elside(id2,ie)

            do m=1,i34(ie) !side
                n1=isidenode(1,elside(m,ie))
                n2=isidenode(2,elside(m,ie))
                tmpx2 = trc_o(n1)
                tmpy2 = trc_o(n2)
                swild(m) = (tmpx2+tmpy2)*0.5
            enddo
            trctmp=sum(swild(1:i34(ie)))/3
            !1st segment
            tmpnode = flux_matrix00(j,i,1)
            tmpx2 = trc_o(tmpnode)
            flux_matrix1(j,i,trc_n) = flux_matrix1(j,i,trc_d)*tmpx2
            flux_matrix2(j,i,trc_n) = flux_matrix2(j,i,trc_d)*tmpx2
            sumtmp = (flux_matrix1(j,i,trc_n)+flux_matrix2(j,i,trc_n))/2.d0
            nodeid = elnode(nxq(1,id,i34(ie)),ie)
            fluxchan1=fluxchan1+(flux_matrix1(j,i,trc_n)+flux_matrix2(j,i,trc_n))/2.d0
             
            
            !2nd segment
            !Normal dir x length
!            xtmp=ycj(isd2)-yctr(ie)
!            ytmp=xctr(ie)-xcj(isd2)
            tmpnode = flux_matrix00(j,i,2)
            tmpx2 = trc_o(tmpnode)
            flux_matrix3(j,i,trc_n) = flux_matrix3(j,i,trc_d)*tmpx2
            flux_matrix4(j,i,trc_n) = flux_matrix4(j,i,trc_d)*tmpx2
            sumtmp = (flux_matrix3(j,i,trc_n) + flux_matrix4(j,i,trc_n))/2.d0
            nodeid = elnode(nxq(i34(ie)-1,id,i34(ie)),ie)

            fluxchan1=fluxchan1+(flux_matrix3(j,i,trc_n) + flux_matrix4(j,i,trc_n))/2.d0
          enddo !j=1,nne(i)
          
          if(suma/=0) then
            trl0(i)=fluxchan1/suma !m/s/s
          endif !suma/=0
      enddo !i=1,np
      call exchange_p2d(trc)
      call exchange_p2d(d_tr0)
      call exchange_p2d(trl0)
      do i = 1,nx_nh
       trc(i) = d_tr0(i) - trl0(i) !+ d_tr(i)


      enddo

      call exchange_p2d(trc)
    end subroutine tvd_casulli_icepack_dep

!=======================================================================
    
    
    subroutine upwind_icepack4(trc,trc_n,adv_threshold0)
        
      !use mod_mesh     
       use mice_module
      implicit none
      integer (kind=int_kind), intent(in) :: trc_n
      real(kind=dbl_kind), dimension(nx), intent(inout) :: trc,adv_threshold0  
      real   (kind=dbl_kind)            :: suma,xctr2,yctr2,vnor1,vnor2,sumtmp,areatmp,&
                                           & tmpx2,tmpx3,tmpy2,tmpy3,utmp,vtmp,xtmp,ytmp,vnorm,&
                                           & trctmp, trcmin, trcmax ,puny ,tmp1 ,fluxchan0, tmp2, adv_threshold,&
                                           & tmp3
      integer (kind=int_kind) :: i,j,ie,id,id2,id3,isd3,isd2,     &
                                 m,n1,n2,jj,indx1,indx2,nd,nodeid
      logical (kind=log_kind) :: flag_noadv
       real(kind=dbl_kind) :: swild10(3,2),swild(3),swild2(3,2),trl0(nx),d_tr0(nx),fluxchan1(nx),fluxchan2(nx)
      !type(t_mesh),        target,        intent(in)    :: mesh
       call icepack_query_parameters(puny_out=puny)

    adv_threshold = 10
    trl0(:) = c0
    fluxchan1(:) = c0 !positive flux
    fluxchan2(:) = c0 !negative flux
    d_tr0(:) = trc(:)
    flux_matrix1(:,:,trc_n) = c0
    flux_matrix2(:,:,trc_n) = c0
    flux_matrix3(:,:,trc_n) = c0
    flux_matrix4(:,:,trc_n) = c0
    flux_matrix0(:,:,:,trc_n) = c0
      ! Driving sequence
    do ie=1,nx_elem

          if(idry_e(ie)==1) cycle
          areatmp = area(ie)/6
          !Wet elem

          do m=1,i34(ie) !side
            n1=isidenode(1,elside(m,ie))
            n2=isidenode(2,elside(m,ie))
            if(ics==1) then
               swild10(m,1) = 0.5*(uvel(n1)+uvel(n2))
               swild10(m,2) = 0.5*(vvel(n1)+vvel(n2))
               swild(m) = (trc(n1)+trc(n2))*0.5
               swild2(m,1) = uvel(elnode(m,ie))
               swild2(m,2) = vvel(elnode(m,ie))
            else
               swild10(m,1) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),eframe(:,1,ie))+vvel(n1)*dot_product(pframe(:,2,n1),eframe(:,1,ie))&
                              &+uvel(n2)*dot_product(pframe(:,1,n2),eframe(:,1,ie))+vvel(n2)*dot_product(pframe(:,2,n2),eframe(:,1,ie)))
               swild10(m,2) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),eframe(:,2,ie))+vvel(n1)*dot_product(pframe(:,2,n1),eframe(:,2,ie))&
                              &+uvel(n2)*dot_product(pframe(:,1,n2),eframe(:,2,ie))+vvel(n2)*dot_product(pframe(:,2,n2),eframe(:,2,ie)))                      
               swild(m) = (trc(n1)+trc(n2))*0.5
               swild2(m,1) = uvel(elnode(m,ie))*dot_product(pframe(:,1,elnode(m,ie)),eframe(:,1,ie))+vvel(elnode(m,ie))*dot_product(pframe(:,2,elnode(m,ie)),eframe(:,1,ie))
               swild2(m,2) = uvel(elnode(m,ie))*dot_product(pframe(:,1,elnode(m,ie)),eframe(:,2,ie))+vvel(elnode(m,ie))*dot_product(pframe(:,2,elnode(m,ie)),eframe(:,2,ie))
            endif
            !swild10(m,1) = 0.5*(uvel(n1)*trc(n1)*dot_product(pframe(:,1,n1),sframe2(:,1,elside(m,ie)))+vvel(n1)*trc(n1)*dot_product(pframe(:,2,n1),sframe2(:,1,elside(m,ie)))&
            !                 &+uvel(n2)*trc(n2)*dot_product(pframe(:,1,n2),sframe2(:,1,elside(m,ie)))+vvel(n2)*trc(n2)*dot_product(pframe(:,2,n1),sframe2(:,1,elside(m,ie))))
            !swild10(m,2) = 0.5*(uvel(n1)*trc(n1)*dot_product(pframe(:,1,n1),sframe2(:,2,elside(m,ie)))+vvel(n1)*trc(n1)*dot_product(pframe(:,2,n1),sframe2(:,2,elside(m,ie)))&
            !                 &+uvel(n2)*trc(n2)*dot_product(pframe(:,1,n2),sframe2(:,2,elside(m,ie)))+vvel(n2)*trc(n2)*dot_product(pframe(:,2,n1),sframe2(:,2,elside(m,ie))))
            !swild10(m,1) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),sframe2(:,1,elside(m,ie)))+vvel(n1)*dot_product(pframe(:,2,n1),sframe2(:,1,elside(m,ie)))&
            !                 &+uvel(n2)*dot_product(pframe(:,1,n2),sframe2(:,1,elside(m,ie)))+vvel(n2)*dot_product(pframe(:,2,n2),sframe2(:,1,elside(m,ie))))
            !swild10(m,2) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),sframe2(:,2,elside(m,ie)))+vvel(n1)*dot_product(pframe(:,2,n1),sframe2(:,2,elside(m,ie)))&
            !                 &+uvel(n2)*dot_product(pframe(:,1,n2),sframe2(:,2,elside(m,ie)))+vvel(n2)*dot_product(pframe(:,2,n2),sframe2(:,2,elside(m,ie))))
            !swild10(m,1) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),eframe(:,1,ie))+vvel(n1)*dot_product(pframe(:,2,n1),eframe(:,1,ie))&
            !                 &+uvel(n2)*dot_product(pframe(:,1,n2),eframe(:,1,ie))+vvel(n2)*dot_product(pframe(:,2,n2),eframe(:,1,ie)))
            !swild10(m,2) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),eframe(:,2,ie))+vvel(n1)*dot_product(pframe(:,2,n1),eframe(:,2,ie))&
            !                 &+uvel(n2)*dot_product(pframe(:,1,n2),eframe(:,2,ie))+vvel(n2)*dot_product(pframe(:,2,n2),eframe(:,2,ie)))                      
            !swild(m) = (trc(n1)+trc(n2))*0.5
            !swild2(m,1) = uvel(elnode(m,ie))*dot_product(pframe(:,1,elnode(m,ie)),eframe(:,1,ie))+vvel(elnode(m,ie))*dot_product(pframe(:,2,elnode(m,ie)),eframe(:,1,ie))
            !swild2(m,2) = uvel(elnode(m,ie))*dot_product(pframe(:,1,elnode(m,ie)),eframe(:,2,ie))+vvel(elnode(m,ie))*dot_product(pframe(:,2,elnode(m,ie)),eframe(:,2,ie))
        enddo
        !utmp = sum(swild10(1:i34(ie),1))/3 !vel @ centroid
        !vtmp = sum(swild10(1:i34(ie),2))/3 !vel @ centroid
        utmp=sum(swild2(1:i34(ie),1))/3 !vel @ centroid
        vtmp=sum(swild2(1:i34(ie),2))/3 !vel @ centroid
        trctmp = sum(swild(1:i34(ie)))/3

          do j =1,i34(ie) !side
            n1=isidenode(1,elside(j,ie))
            n2=isidenode(2,elside(j,ie))
            
            indx1 = 0
            do m=1,nne(n1)
               if(indel(m,n1)==ie) then
                  indx1=m; exit
               endif
            enddo
            if(indx1==0) call parallel_abort('STEP: failed to find (upwind_icepack_1)') 
            indx2 = 0
            do m=1,nne(n2)
               if(indel(m,n2)==ie) then
                  indx2=m; exit
               endif
            enddo
            if(indx2==0) call parallel_abort('STEP: failed to find (upwind_icepack_2)') 

            id  = iself(indx1,n1)
            id2 = iself(indx2,n2)
            id3 = nxq(i34(ie)-2,id,i34(ie))

            if(ics==1) then
               xctr2=xctr(ie); yctr2=yctr(ie)
               tmpx3=(xnd(n1)+xnd(n2))/2.d0
               tmpy3=(ynd(n1)+ynd(n2))/2.d0
            else!ll; use [xy]el defined in eframe
               xctr2=0.d0 !sum(xel(elnode(1:i34(ie),ie)))/i34(ie)
               yctr2=0.d0 !sum(yel(elnode(1:i34(ie),ie)))/i34(ie)
               tmpx3=(xel(id,ie)+xel(id2,ie))/2.d0
               tmpy3=(yel(id,ie)+yel(id2,ie))/2.d0
            endif!ics
            !xtmp=yctr2-tmpy3
            !ytmp=tmpx3-xctr2
            !vnor1=utmp*xtmp+vtmp*ytmp !normal vel x length 
            !vnor2=swild10(j,1)*xtmp+swild10(j,2)*ytmp !normal vel@side x length 
            if(id3.eq.id2) then
               xtmp=yctr2-tmpy3
               ytmp=tmpx3-xctr2
               vnor1=utmp*xtmp+vtmp*ytmp !normal vel x length 
               vnor2=swild10(j,1)*xtmp+swild10(j,2)*ytmp !normal vel@side x length 
               flux_matrix1(indx1,n1,trc_n) = dt_dyn*trctmp*vnor1
               flux_matrix2(indx1,n1,trc_n) = dt_dyn*swild(j)*vnor2
               flux_matrix3(indx2,n2,trc_n) = -1*dt_dyn*trctmp*vnor1
               flux_matrix4(indx2,n2,trc_n) = -1*dt_dyn*swild(j)*vnor2
               sumtmp = (dt_dyn*trctmp*vnor1+dt_dyn*swild(j)*vnor2)/2
               if(sumtmp>0) then
                  tmp3 = d_tr0(n2)-d_tr0(n1)
                  tmp2=max(c0,tmp3)
                  !adv_threshold = exp(-20*tmp2)
                  !if((sumtmp/areatmp)>(d_tr0(n1)*adv_threshold))then
                  !if((sumtmp/areatmp)>(d_tr0(n1)*adv_threshold0(n1)).or.(sumtmp/areatmp)>(0.01))then
                  if((sumtmp/areatmp)>(d_tr0(n1)*adv_threshold).or.(sumtmp/areatmp)>(d_tr0(n1)*adv_threshold0(n2)))then
                     !tmp1 = areatmp*d_tr0(n1)/sumtmp*adv_threshold
                     tmp2 = min(d_tr0(n1)*adv_threshold,d_tr0(n1)*adv_threshold0(n2))
                     tmp1 = areatmp*tmp2/sumtmp
                     flux_matrix1(indx1,n1,trc_n) = flux_matrix1(indx1,n1,trc_n)*tmp1
                     flux_matrix2(indx1,n1,trc_n) = flux_matrix2(indx1,n1,trc_n)*tmp1
                     flux_matrix3(indx2,n2,trc_n) = flux_matrix3(indx2,n2,trc_n)*tmp1
                     flux_matrix4(indx2,n2,trc_n) = flux_matrix4(indx2,n2,trc_n)*tmp1
                  endif
                  
               else
                  tmp3 = d_tr0(n1)-d_tr0(n2)
                  tmp2=max(c0,tmp3)
                  !adv_threshold = exp(-20*tmp2)
                  !if((abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold))then
                  !if((abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold0(n2)).or.(abs(sumtmp)/areatmp)>(0.01))then
                  if((abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold).or.(abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold0(n1)))then
                     !tmp1 = areatmp*d_tr0(n2)/abs(sumtmp)*adv_threshold
                     tmp2 = min(d_tr0(n2)*adv_threshold,d_tr0(n2)*adv_threshold0(n1))
                     tmp1 = areatmp*tmp2/abs(sumtmp)
                     flux_matrix1(indx1,n1,trc_n) = flux_matrix1(indx1,n1,trc_n)*tmp1
                     flux_matrix2(indx1,n1,trc_n) = flux_matrix2(indx1,n1,trc_n)*tmp1
                     flux_matrix3(indx2,n2,trc_n) = flux_matrix3(indx2,n2,trc_n)*tmp1
                     flux_matrix4(indx2,n2,trc_n) = flux_matrix4(indx2,n2,trc_n)*tmp1
                  endif
               endif
            else
               xtmp = tmpy3 - yctr2
               ytmp = xctr2 - tmpx3
               vnor1=utmp*xtmp+vtmp*ytmp !normal vel x length 
               vnor2=swild10(j,1)*xtmp+swild10(j,2)*ytmp !normal vel@side x length 

               flux_matrix3(indx1,n1,trc_n) = dt_dyn*trctmp*vnor1
               flux_matrix4(indx1,n1,trc_n) = dt_dyn*swild(j)*vnor2
               flux_matrix1(indx2,n2,trc_n) = -1*dt_dyn*trctmp*vnor1
               flux_matrix2(indx2,n2,trc_n) = -1*dt_dyn*swild(j)*vnor2
               sumtmp = (dt_dyn*trctmp*vnor1+dt_dyn*swild(j)*vnor2)/2
               if(sumtmp>0) then
                  tmp3 = d_tr0(n2)-d_tr0(n1)
                  tmp2=max(c0,tmp3)
                  !adv_threshold = exp(-20*tmp2)
                  !if((sumtmp/areatmp)>(d_tr0(n1)*adv_threshold))then
                  !if((sumtmp/areatmp)>(d_tr0(n1)*adv_threshold0(n1)).or.(sumtmp/areatmp)>(0.01))then
                  if((sumtmp/areatmp)>(d_tr0(n1)*adv_threshold).or.(sumtmp/areatmp)>(d_tr0(n1)*adv_threshold0(n2)))then
                     !tmp1 = areatmp*d_tr0(n1)/sumtmp*adv_threshold
                     tmp2 = min(d_tr0(n1)*adv_threshold,d_tr0(n1)*adv_threshold0(n2))
                     tmp1 = areatmp*tmp2/sumtmp
                     flux_matrix3(indx1,n1,trc_n) = flux_matrix3(indx1,n1,trc_n)*tmp1
                     flux_matrix4(indx1,n1,trc_n) = flux_matrix4(indx1,n1,trc_n)*tmp1
                     flux_matrix1(indx2,n2,trc_n) = flux_matrix1(indx2,n2,trc_n)*tmp1
                     flux_matrix2(indx2,n2,trc_n) = flux_matrix2(indx2,n2,trc_n)*tmp1
                  endif
                  
               else
                  tmp3 = d_tr0(n1)-d_tr0(n2)
                  tmp2=max(c0,tmp3)
                  !adv_threshold = exp(-20*tmp2)
                  !if((abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold))then
                  !if((abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold0(n2)).or.(abs(sumtmp)/areatmp)>0.01)then
                  if((abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold).or.(abs(sumtmp)/areatmp)>(d_tr0(n2)*adv_threshold0(n1)))then
                     !tmp1 = areatmp*d_tr0(n2)/abs(sumtmp)*adv_threshold
                     tmp2 = min(d_tr0(n2)*adv_threshold,d_tr0(n2)*adv_threshold0(n1))
                     tmp1 = areatmp*tmp2/abs(sumtmp)
                     flux_matrix3(indx1,n1,trc_n) = flux_matrix3(indx1,n1,trc_n)*tmp1
                     flux_matrix4(indx1,n1,trc_n) = flux_matrix4(indx1,n1,trc_n)*tmp1
                     flux_matrix1(indx2,n2,trc_n) = flux_matrix1(indx2,n2,trc_n)*tmp1
                     flux_matrix2(indx2,n2,trc_n) = flux_matrix2(indx2,n2,trc_n)*tmp1
                  endif
               endif
            endif
          enddo
    enddo
   do i=1,nx_nh
      suma=0.d0 !sum of areas
      fluxchan0=0.d0
      do j=1,nne(i)
         ie=indel(j,i)
         suma=suma+area(ie)/3 !approx for quad
         fluxchan0=fluxchan0 + (flux_matrix1(j,i,trc_n) + flux_matrix2(j,i,trc_n))/2.d0
         fluxchan0=fluxchan0 + (flux_matrix3(j,i,trc_n) + flux_matrix4(j,i,trc_n))/2.d0

      enddo
         if( suma/=0 ) then
            trl0(i)=fluxchan0/suma !m/s/s
         endif !suma/=0  
          !fluxchan1(i) = sum(flux_matrix1(1:nne(i),i,trc_n))/2.d0 + sum(flux_matrix2(1:nne(i),i,trc_n))/2.d0 + &
          !             & sum(flux_matrix3(1:nne(i),i,trc_n))/2.d0 + sum(flux_matrix4(1:nne(i),i,trc_n))/2.d0
          !tmp1 = fluxchan1(i)/suma
          !tmp2 = trl0(i) - tmp1
          !if(abs(tmp2) > puny) write(12,*)'upwind fluxchan1',i,tmp2,trl0(i),tmp1,fluxchan0/suma,fluxchan0,fluxchan1(i),suma,area_median(i)
         trc(i) = d_tr0(i) - trl0(i) !+ d_tr(i)   
         !if(trc(i)<0) then
         !   write(12,*)'upwind mystery2',i,trc(i),d_tr0(i),trl0(i)
         !endif 
   enddo  
   !call exchange_p2d(trl0)
   ! do i = 1,nx_nh
   !    fluxchan1(i) = sum(flux_matrix1(1:nne(i),i,trc_n))/2.d0 + sum(flux_matrix2(1:nne(i),i,trc_n))/2.d0 + &
   !                 & sum(flux_matrix3(1:nne(i),i,trc_n))/2.d0 + sum(flux_matrix4(1:nne(i),i,trc_n))/2.d0
   !    trl0(i) = fluxchan1(i)/area_median(i)
       !trc(i) = d_tr0(i) - trl0(i) !+ d_tr(i)

       !if(trc(i)<0) then
       !   write(12,*)'upwind mystery2',i,trc(i),d_tr0(i),trl0(i)
       !endif
   ! enddo
    call exchange_p2d(trc)
  end subroutine upwind_icepack4

  !=======================================================================
subroutine upwind_icepack_other4(trc,trc_o,trc_n,trc_d,trc_d1,trc_d2)
        
        !use mod_mesh     
     
        implicit none
        integer (kind=int_kind), intent(in) :: trc_n,trc_d
        real (kind=dbl_kind), dimension(nx), intent(inout) :: trc   ,trc_o
        real (kind=dbl_kind), dimension(nx), intent(in), optional :: trc_d1,trc_d2
        real   (kind=dbl_kind)            :: fluxchan1,fluxchan2,suma,xctr2,yctr2,vnor1,vnor2,&
                                             & tmpx2,tmpx3,tmpy2,tmpy3,utmp,vtmp,xtmp,ytmp,vnorm,&
                                             & trctmp  , trcmin, trcmax   
        integer (kind=int_kind) :: i,j,ie,id,id2,id3,isd3,isd2,     &
                                   m,n1,n2,jj
        logical (kind=log_kind) :: flag_noadv
         real(rkind) :: swild10(3,2),swild(3),trl0(nx),d_tr0(nx)
        !type(t_mesh),        target,        intent(in)    :: mesh

      trl0(:) = c0
      d_tr0(:) = trc(:)
      flux_matrix1(:,:,trc_n) = c0
      flux_matrix2(:,:,trc_n) = c0
      flux_matrix3(:,:,trc_n) = c0
      flux_matrix4(:,:,trc_n) = c0
      flux_matrix0(:,:,:,trc_n) = c0
        ! Driving sequence
      do i=1,nx_nh
        if(idry(i)==1) cycle
        !if(isbnd(1,i)/=0) cycle
        !Wet node
          fluxchan1=0.d0 !sum of fluxes;\int_\Gamma u*u_n d\Gamma
          fluxchan2=0.d0 !sum of fluxes;
          suma=0.d0 !sum of areas
          do j=1,nne(i)
            ie=indel(j,i)
            if(idry_e(ie)==1) cycle
 
            !Wet elem
            suma=suma+area(ie)/3 !approx for quad
            id=iself(j,i)
            id3=nxq(i34(ie)-1,id,i34(ie)) !adjacent side index
            id2=nxq(i34(ie)-2,id,i34(ie)) !adjacent side index
            isd3=elside(id3,ie)
            isd2=elside(id2,ie)

            do m=1,i34(ie) !side
                n1=isidenode(1,elside(m,ie))
                n2=isidenode(2,elside(m,ie))
                tmpx2 = trc_o(n1)
                tmpy2 = trc_o(n2)
                if (present(trc_d1)) then
                  tmpx2 = tmpx2 * trc_d1(n1)
                  tmpy2 = tmpy2 * trc_d1(n2)
                endif
                if (present(trc_d2)) then
                  tmpx2 = tmpx2 * trc_d2(n1)
                  tmpy2 = tmpy2 * trc_d2(n2)
                endif
                swild(m) = (tmpx2+tmpy2)*0.5
            enddo
            trctmp = sum(swild(1:i34(ie)))/3
            trcmin = minval(trc_o(indnd(1:nnp(i),i)))
            trcmin = min(trcmin,trc_o(i))
            trcmax = maxval(trc_o(indnd(1:nnp(i),i)))
            trcmax = max(trcmax,trc_o(i))

            !1st segment
            flux_matrix1(j,i,trc_n) = flux_matrix1(j,i,trc_d)*trctmp
            flux_matrix2(j,i,trc_n) = flux_matrix2(j,i,trc_d)*swild(id3)
            !fluxchan1=fluxchan1+(flux_matrix2(j,i,trc_d)*swild(id3))
            fluxchan1=fluxchan1+(flux_matrix1(j,i,trc_d)*trctmp+flux_matrix2(j,i,trc_d)*swild(id3))/2.d0
            !fluxchan1=fluxchan1+(flux_matrix1(j,i,trc_d)+flux_matrix2(j,i,trc_d))*(trctmp+swild(id3))/4.d0
             
            
            !2nd segment
            !Normal dir x length
!            xtmp=ycj(isd2)-yctr(ie)
!            ytmp=xctr(ie)-xcj(isd2)
            flux_matrix3(j,i,trc_n) = flux_matrix3(j,i,trc_d)*trctmp
            flux_matrix4(j,i,trc_n) = flux_matrix4(j,i,trc_d)*swild(id2)
            !fluxchan1=fluxchan1+(flux_matrix4(j,i,trc_d)*swild(id2))
            fluxchan1=fluxchan1+(flux_matrix3(j,i,trc_d)*trctmp+flux_matrix4(j,i,trc_d)*swild(id2))/2.d0
            !fluxchan1=fluxchan1+(flux_matrix3(j,i,trc_d)+flux_matrix4(j,i,trc_d))*(trctmp+swild(id2))/4.d0
          enddo !j=1,nne(i)
          
          if(suma/=0) then
            trl0(i)=fluxchan1/suma !m/s/s
            !trl0(i)=fluxchan1/aream(i) !m/s/s
            !d_tr(i)=fluxchan2/suma
          endif !suma/=0
      enddo !i=1,np
      call exchange_p2d(trc)
      call exchange_p2d(d_tr0)
      call exchange_p2d(trl0)
      do i = 1,nx_nh
       trc(i) = d_tr0(i) - trl0(i) !+ d_tr(i)
       !if(trc(i)*d_tr0(i)<0) write(12,*)'upwind_other mystery',i,trc(i),d_tr0(i),trl0(i)
       !if(trc(i)*d_tr(i)<0) trc(i) = d_tr(i)
      enddo

      call exchange_p2d(trc)
    end subroutine upwind_icepack_other4

    !=======================================================================
    subroutine upwind_icepack_dep4(trc,trc_o,trc_n,trc_d)
        
        !use mod_mesh     
     
        implicit none
        integer (kind=int_kind), intent(in) :: trc_n,trc_d
        real (kind=dbl_kind), dimension(nx), intent(inout) :: trc   ,trc_o
        real   (kind=dbl_kind)            :: fluxchan1,fluxchan2,suma,xctr2,yctr2,vnor1,vnor2,&
                                             & tmpx2,tmpx3,tmpy2,tmpy3,utmp,vtmp,xtmp,ytmp,vnorm,&
                                             & trctmp, trcmin, trcmax ,sumtmp,areatmp  
        integer (kind=int_kind) :: i,j,ie,id,id2,id3,isd3,isd2,     &
                                   m,n1,n2,jj,indx,nd,nodeid
        logical (kind=log_kind) :: flag_noadv
         real(rkind) :: swild10(3,2),swild(3),trl0(nx),d_tr0(nx)
        !type(t_mesh),        target,        intent(in)    :: mesh

      trl0(:) = c0
      d_tr0(:) = trc(:)
      flux_matrix1(:,:,trc_n) = c0
      flux_matrix2(:,:,trc_n) = c0
      flux_matrix3(:,:,trc_n) = c0
      flux_matrix4(:,:,trc_n) = c0
      flux_matrix0(:,:,:,trc_n) = c0
        ! Driving sequence
      do i=1,nx_nh
        if(idry(i)==1) cycle
        !Wet node
          fluxchan1=0.d0 !sum of fluxes;\int_\Gamma u*u_n d\Gamma
          fluxchan2=0.d0 !sum of fluxes;
          suma=0.d0 !sum of areas
          do j=1,nne(i)
            ie=indel(j,i)
            if(idry_e(ie)==1) cycle
 
            !Wet elem
            areatmp = area(ie)/6
            suma=suma+area(ie)/3 !approx for quad
            id=iself(j,i)
            id3=nxq(i34(ie)-1,id,i34(ie)) !adjacent side index
            id2=nxq(i34(ie)-2,id,i34(ie)) !adjacent side index
            isd3=elside(id3,ie)
            isd2=elside(id2,ie)

            do m=1,i34(ie) !side
                n1=isidenode(1,elside(m,ie))
                n2=isidenode(2,elside(m,ie))
                tmpx2 = trc_o(n1)
                tmpy2 = trc_o(n2)
                swild(m) = (tmpx2+tmpy2)*0.5
            enddo
            trctmp=sum(swild(1:i34(ie)))/3
            !1st segment
            flux_matrix1(j,i,trc_n) = flux_matrix1(j,i,trc_d)*trctmp
            flux_matrix2(j,i,trc_n) = flux_matrix2(j,i,trc_d)*swild(id3)
            sumtmp = (flux_matrix1(j,i,trc_n)+flux_matrix2(j,i,trc_n))/2.d0
            nodeid = elnode(nxq(1,id,i34(ie)),ie)
            !if(sumtmp>0) then
            !   if(sumtmp/areatmp>d_tr0(i)*0.9)then
            !      flux_matrix1(j,i,trc_n) = flux_matrix1(j,i,trc_n)*areatmp*d_tr0(i)/sumtmp*0.9
            !      flux_matrix2(j,i,trc_n) = flux_matrix2(j,i,trc_n)*areatmp*d_tr0(i)/sumtmp*0.9
            !   endif
            !   
            !else
            !   if(sumtmp/areatmp>d_tr0(nodeid)*0.9)then
            !      flux_matrix1(j,i,trc_n) = flux_matrix1(j,i,trc_n)*areatmp*d_tr0(nodeid)/sumtmp*0.9
            !      flux_matrix2(j,i,trc_n) = flux_matrix2(j,i,trc_n)*areatmp*d_tr0(nodeid)/sumtmp*0.9
            !   endif
            !endif
            !fluxchan1=fluxchan1+(flux_matrix2(j,i,trc_d)*swild(id3))
            fluxchan1=fluxchan1+(flux_matrix1(j,i,trc_n)+flux_matrix2(j,i,trc_n))/2.d0
            !fluxchan1=fluxchan1+(flux_matrix1(j,i,trc_d)+flux_matrix2(j,i,trc_d))*(trctmp+swild(id3))/4.d0
             
            
            !2nd segment
            !Normal dir x length
!            xtmp=ycj(isd2)-yctr(ie)
!            ytmp=xctr(ie)-xcj(isd2)
            flux_matrix3(j,i,trc_n) = flux_matrix3(j,i,trc_d)*trctmp
            flux_matrix4(j,i,trc_n) = flux_matrix4(j,i,trc_d)*swild(id2)
            sumtmp = (flux_matrix3(j,i,trc_n) + flux_matrix4(j,i,trc_n))/2.d0
            nodeid = elnode(nxq(i34(ie)-1,id,i34(ie)),ie)
            !if(sumtmp>0) then
            !   if(sumtmp/areatmp>d_tr0(i)*0.9)then
            !      flux_matrix3(j,i,trc_n) = flux_matrix3(j,i,trc_n)*areatmp*d_tr0(i)/sumtmp*0.9
            !      flux_matrix4(j,i,trc_n) = flux_matrix4(j,i,trc_n)*areatmp*d_tr0(i)/sumtmp*0.9
            !   endif
            !   
            !else
            !   if(sumtmp/areatmp>d_tr0(nodeid)*0.9)then
            !      flux_matrix3(j,i,trc_n) = flux_matrix3(j,i,trc_n)*areatmp*d_tr0(nodeid)/sumtmp*0.9
            !      flux_matrix4(j,i,trc_n) = flux_matrix4(j,i,trc_n)*areatmp*d_tr0(nodeid)/sumtmp*0.9
            !   endif
            !endif
            !fluxchan1=fluxchan1+(flux_matrix4(j,i,trc_d)*swild(id2))
            fluxchan1=fluxchan1+(flux_matrix3(j,i,trc_n) + flux_matrix4(j,i,trc_n))/2.d0
            !fluxchan1=fluxchan1+(flux_matrix3(j,i,trc_d)+flux_matrix4(j,i,trc_d))*(trctmp+swild(id2))/4.d0
          enddo !j=1,nne(i)
          
          if(suma/=0) then
            trl0(i)=fluxchan1/suma !m/s/s
            !trl0(i)=fluxchan1/aream(i) !m/s/s
            !d_tr(i)=fluxchan2/suma
          endif !suma/=0
      enddo !i=1,np
      call exchange_p2d(trc)
      call exchange_p2d(d_tr0)
      call exchange_p2d(trl0)
      do i = 1,nx_nh
       trc(i) = d_tr0(i) - trl0(i) !+ d_tr(i)

      !   if(trc(i)<0) then
      !   write(12,*)'upwind mystery vice',i,trc(i),d_tr(i),trl0(i)
      !      flux_matrix1(:,i,trc_n) = 0
      !      flux_matrix2(:,i,trc_n) = 0
      !      flux_matrix3(:,i,trc_n) = 0
      !      flux_matrix4(:,i,trc_n) = 0
      !      trc(i) = d_tr0(i)
      !      do j=1,nne(i)
      !         ie=indel(j,i)
      !         id=iself(j,i)
      !         !id=indnd(j,i)
      !         do jj=1,2 !other 2 nodes
      !            nd=elnode(nxq(jj,id,i34(ie)),ie)
      !            indx=0
      !            do m=1,nne(nd)
      !            if(indel(m,nd)==ie) then
      !               indx=m; exit
      !            endif
      !            enddo !m
      !            if(indx==0) call parallel_abort('STEP: failed to find (9.1)')  
      !            flux_matrix1(indx,nd,trc_n) = 0
      !            flux_matrix2(indx,nd,trc_n) = 0
      !            flux_matrix3(indx,nd,trc_n) = 0
      !            flux_matrix4(indx,nd,trc_n) = 0
      !         enddo
      !      enddo
      !   endif
      enddo

      call exchange_p2d(trc)
    end subroutine upwind_icepack_dep4

!=======================================================================
    
    subroutine upwind_icepack2(trc,trc_n)
        
        !use mod_mesh     
     
        implicit none
        integer (kind=int_kind), intent(in) :: trc_n
        real(kind=dbl_kind), dimension(nx), intent(inout) :: trc  
        real   (kind=dbl_kind)            :: fluxchan1,fluxchan2,suma,xctr2,yctr2,vnor1,vnor2,&
                                             & tmpx2,tmpx3,tmpy2,tmpy3,utmp,vtmp,xtmp,ytmp,vnorm,&
                                             & trctmp, trcmin, trcmax ,puny  
        integer (kind=int_kind) :: i,j,ie,id,id2,id3,isd3,isd2,     &
                                   m,n1,n2,jj,indx,nd,row
        logical (kind=log_kind) :: flag_noadv
         real(rkind) :: swild10(3,2),swild(3),swild2(3,2),trl0(nx),d_tr0(nx)
        !type(t_mesh),        target,        intent(in)    :: mesh
         call icepack_query_parameters(puny_out=puny)
      trl0(:) = c0
      d_tr0(:) = trc(:)
      flux_matrix0(:,:,:,trc_n) = c0
      call tg_rhs_div_icepack_upwind(trc,trc_n)
      !do i = 1,nx_nh
         !if(isbnd(1,i)/=0) cycle
         !if(any(isbnd(1,indnd(1:nnp(i),i))/=0)) cycle
      !   fluxchan2=0.d0 !sum of fluxes;
      !   do j=1,nne(i)
      !      ie=indel(j,i)
      !      fluxchan2 = fluxchan2+sum(flux_matrix0(1:3,j,i,trc_n))
      !   enddo
      !      trc(i) = d_tr0(i) + fluxchan2 !+ d_tr(i)
      !   if(trc(i)<0) then 
      !      !write(12,*) i, trc(i), fluxchan2
      !      do j=1,nne(i)
      !         ie=indel(j,i)
      !         flux_matrix0(1:3,j,i,trc_n) = flux_matrix0(1:3,j,i,trc_n)*abs(d_tr0(i)/fluxchan2)
      !      enddo
      !      trc(i) = 0
      !   endif
      ! enddo
      ! call exchange_p2d(trc)
      ! trl0(:) = trc
      ! trc = d_tr0
        do i = 1, nx_elem
            do j = 1, i34(i)
               row = elnode(j,i)
                  indx=0
                  do m=1,nne(row)
                  if(indel(m,row)==i) then
                     indx=m; exit
                  endif
                  enddo !m
                  if(indx==0) call parallel_abort('STEP: failed to find (upwind_icepack2)')  

               flux_matrix0(:,indx,row,trc_n) = flux_matrix0(:,indx,row,trc_n)
               trc(elnode(j,i)) = trc(elnode(j,i)) + sum(flux_matrix0(1:3,indx,row,trc_n))
            enddo
   
         enddo

      call exchange_p2d(trc)
      !do i = 1, nx
      !   trctmp = trl0(i) - trc(i)
      !   if(abs(trctmp)>puny) write(12,*) 'diff in node and elem.', i,trctmp,trl0(i),trc(i)
      !enddo
    end subroutine upwind_icepack2

    !=======================================================================
subroutine upwind_icepack_other2(trc,trc_o,trc_n,trc_d,trc_d1,trc_d2)
        
        !use mod_mesh     
     
        implicit none
        integer (kind=int_kind), intent(in) :: trc_n,trc_d
        real (kind=dbl_kind), dimension(nx), intent(inout) :: trc   ,trc_o
        real (kind=dbl_kind), dimension(nx), intent(in), optional :: trc_d1,trc_d2
        real   (kind=dbl_kind)            :: fluxchan1,fluxchan2,suma,xctr2,yctr2,vnor1,vnor2,&
                                             & tmpx2,tmpx3,tmpy2,tmpy3,utmp,vtmp,xtmp,ytmp,vnorm,&
                                             & trctmp  , trcmin, trcmax   
        integer (kind=int_kind) :: i,j,ie,id,id2,id3,isd3,isd2,     &
                                   m,n1,n2,jj,indx,row
        logical (kind=log_kind) :: flag_noadv
         real(rkind) :: swild10(3,2),swild(3),trl0(nx),d_tr0(nx)
        !type(t_mesh),        target,        intent(in)    :: mesh

      trl0(:) = c0
      d_tr0(:) = trc(:)
      flux_matrix0(:,:,:,trc_n) = c0
      
        ! Driving sequence

      !do i = 1,nx_nh
         !if(isbnd(1,i)/=0) cycle
         !if(any(isbnd(1,indnd(1:nnp(i),i))/=0)) cycle
      !   fluxchan2=0.d0 !sum of fluxes;
      !   do j=1,nne(i)
      !      ie=indel(j,i)
            !fluxchan2 = fluxchan2+sum(flux_matrix0(1:3,j,i,trc_d)*trc_o(elnode(1:3,ie)))
      !      do jj=1,3
      !         fluxchan2 = fluxchan2+flux_matrix0(jj,j,i,trc_d)*trc_o(elnode(jj,ie))
      !         flux_matrix0(jj,j,i,trc_n)=flux_matrix0(jj,j,i,trc_d)*trc_o(elnode(jj,ie))
      !      enddo
      !   enddo
      !   trc(i) = d_tr0(i) + fluxchan2 !+ d_tr(i)
      ! enddo
         do i = 1,nx_elem
            do j =1 , i34(i)
               row = elnode(j,i)
                  indx=0
                  do m=1,nne(row)
                  if(indel(m,row)==i) then
                     indx=m; exit
                  endif
                  enddo !m
                  if(indx==0) call parallel_abort('STEP: failed to find (upwind_icepack_other2)')  
                  do jj=1,3
                     tmpx2 = trc_o(elnode(jj,i))
                     if (present(trc_d1)) then
                        tmpx2 = tmpx2 * trc_d1(elnode(jj,i))
                     endif
                     if (present(trc_d2)) then
                        tmpx2 = tmpx2 * trc_d2(elnode(jj,i))
                     endif
                     !flux_matrix0(jj,indx,row,trc_n) = flux_matrix0(jj,indx,row,trc_d) * trc_o(elnode(jj,i))
                     !trc(elnode(j,i)) = trc(elnode(j,i)) + flux_matrix0(jj,indx,row,trc_d) * trc_o(elnode(jj,i))
                     flux_matrix0(jj,indx,row,trc_n) = flux_matrix0(jj,indx,row,trc_d) * tmpx2
                     trc(elnode(j,i)) = trc(elnode(j,i)) + flux_matrix0(jj,indx,row,trc_d) * tmpx2
                  enddo
            enddo
   
         enddo
      call exchange_p2d(trc)
    end subroutine upwind_icepack_other2

    !=======================================================================
    subroutine upwind_icepack_dep2(trc,trc_o,trc_n,trc_d)
        
        !use mod_mesh     
     
        implicit none
        integer (kind=int_kind), intent(in) :: trc_n,trc_d
        real (kind=dbl_kind), dimension(nx), intent(inout) :: trc   ,trc_o
        real   (kind=dbl_kind)            :: fluxchan1,fluxchan2,suma,xctr2,yctr2,vnor1,vnor2,&
                                             & tmpx2,tmpx3,tmpy2,tmpy3,utmp,vtmp,xtmp,ytmp,vnorm,&
                                             & trctmp, trcmin, trcmax   
        integer (kind=int_kind) :: i,j,ie,id,id2,id3,isd3,isd2,     &
                                   m,n1,n2,jj,indx,nd,row
        logical (kind=log_kind) :: flag_noadv
         real(rkind) :: swild10(3,2),swild(3),trl0(nx),d_tr0(nx)

      trl0(:) = c0
      d_tr0(:) = trc(:)
      flux_matrix0(:,:,:,trc_n) = c0
        ! Driving sequence
      !do i = 1,nx_nh
         !if(isbnd(1,i)/=0) cycle
         !if(any(isbnd(1,indnd(1:nnp(i),i))/=0)) cycle
      !   fluxchan2=0.d0 !sum of fluxes;
      !   do j=1,nne(i)
      !      ie=indel(j,i)
      !      !fluxchan2 = fluxchan2+sum(flux_matrix0(1:3,j,i,trc_d)*trc_o(elnode(1:3,ie)))
      !      do jj=1,3
      !         fluxchan2 = fluxchan2+flux_matrix0(jj,j,i,trc_d)*trc_o(elnode(jj,ie))
      !         flux_matrix0(jj,j,i,trc_n)=flux_matrix0(jj,j,i,trc_d)*trc_o(elnode(jj,ie))
      !      enddo
      !   enddo
      !   trc(i) = d_tr0(i) + fluxchan2
      !    if(trc(i)<0) then 
      !      do j=1,nne(i)
      !         ie=indel(j,i)
      !         flux_matrix0(1:3,j,i,trc_n) = flux_matrix0(1:3,j,i,trc_n)*abs(d_tr0(i)/fluxchan2)
      !      enddo
      !      trc(i) = 0
      !   endif
      ! enddo
         do i = 1,nx_elem
            do j =1 , i34(i)
               row = elnode(j,i)
                  indx=0
                  do m=1,nne(row)
                  if(indel(m,row)==i) then
                     indx=m; exit
                  endif
                  enddo !m
                  if(indx==0) call parallel_abort('STEP: failed to find (upwind_icepack_dep2)')  
                  do jj=1,3
                     flux_matrix0(jj,indx,row,trc_n) = flux_matrix0(jj,indx,row,trc_d) * trc_o(elnode(jj,i))
                     trc(elnode(j,i)) = trc(elnode(j,i)) + flux_matrix0(jj,indx,row,trc_d) * trc_o(elnode(jj,i))
                  enddo
            enddo
   
         enddo

      call exchange_p2d(trc)
    end subroutine upwind_icepack_dep2

    !=======================================================================
    module subroutine tracer_advection_icepack

        !use mod_mesh
        use mice_module
        use icepack_intfc,        only: icepack_aggregate
        use icepack_itd,          only: cleanup_itd
        use schism_glbl,          only: dt,nstep_ice
        !use g_config,             only: dt
        use icepack_intfc,         only: icepack_init_trcr
        implicit none
            
        integer (kind=int_kind) :: ntrcr, ntrace, narr, nbtrcr, i,     &
                                   nt,   nt1,    k,   n ,ktherm,narrays,it,ierr
        integer (kind=int_kind) :: nt_Tsfc, nt_qice, nt_qsno,                          & 
                                   nt_sice, nt_fbri, nt_iage, nt_FY, nt_alvl, nt_vlvl, &
                                   nt_apnd, nt_hpnd, nt_ipnd, nt_bgc_Nit, nt_bgc_S
        logical (kind=log_kind) :: tr_pond_topo, tr_pond_lvl,                          &
                                   tr_pond,      tr_aero,     tr_FY,                   &
                                   tr_iage,      flag_old,                             &
                                   flag_old_0
        real    (kind=dbl_kind) :: puny ,Tsfc ,rhos, Lfresh, jj, kk, njj, pi, trcmin, trcmax,totalicea,totaliceh, &
        total_a,total_h
      
        ! Tracer dependencies and additional arrays
         real (kind=dbl_kind), dimension(ncat) :: &
               tmp, exc
        integer (kind=int_kind), dimension(:),    allocatable    ::    &
                tracer_type    , & ! = 1, 2, or 3 (depends on 0, 1 or 2 other tracers)
                depend             ! tracer dependencies (see below)
      
        logical (kind=log_kind), dimension (:),   allocatable   ::     &
                has_dependents    ! true if a tracer has dependent tracers
      
        real (kind=dbl_kind),    dimension (:,:), allocatable ::       &
               works,works_old, vicen_init, aicen_init, vsnon_init, aicen_tmp
        real (kind=dbl_kind),    dimension (:), allocatable ::       &
               fiso_ocn,swild,swild1,swild2, iniflag
         real (kind=dbl_kind),    dimension (:,:,:), allocatable ::       &
         trcrn_ini

         real (kind=dbl_kind), dimension(nilyr) :: &
           qin            , & ! ice enthalpy (J/m3)
           qin_max        , & ! maximum ice enthalpy (J/m3)
           zTin               ! initial ice temperature
  
        real (kind=dbl_kind), dimension(nslyr) :: &
           qsn            , & ! snow enthalpy (J/m3)
           zTsn               ! initial snow temperature

        !type(t_mesh),        target,       intent(in)    :: mesh
           character(len=*), parameter :: subname = '(tracer_advection_icepack)'

        call icepack_query_parameters(ktherm_out=ktherm,    &
                                      puny_out=puny,rhos_out=rhos,                   &
                                      Lfresh_out=Lfresh  ,pi_out=pi)
        call icepack_query_tracer_sizes(ntrcr_out=ntrcr, nbtrcr_out=nbtrcr)
        call icepack_query_tracer_flags(                                  &
               tr_iage_out=tr_iage, tr_FY_out=tr_FY,                      &
               tr_aero_out=tr_aero, tr_pond_out=tr_pond,                  &
               tr_pond_lvl_out=tr_pond_lvl, tr_pond_topo_out=tr_pond_topo)
        call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_qice_out=nt_qice, &
               nt_qsno_out=nt_qsno, nt_sice_out=nt_sice, nt_fbri_out=nt_fbri,       &
               nt_iage_out=nt_iage, nt_FY_out=nt_FY, nt_alvl_out=nt_alvl,           &
               nt_vlvl_out=nt_vlvl, nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd,       &
               nt_ipnd_out=nt_ipnd, nt_bgc_Nit_out=nt_bgc_Nit, nt_bgc_S_out=nt_bgc_S)
      
        narr   = 1 + ncat * (3 + ntrcr) ! max number of state variable arrays
      
        ! Allocate works array
        flag_old = .false.
        flag_old_0 = .false.

        if (allocated(works)) deallocate(works)
        if (allocated(works_old)) deallocate(works_old)
        allocate ( works(nx,narr) ) !(npa, max number of state variable arrays)
        allocate ( works_old(nx,narr) )
        allocate ( fiso_ocn(nx) )
        allocate ( swild(nx) )
        allocate ( swild1(nx) )
        allocate ( swild2(nx) )
        allocate ( iniflag(nx) )
        allocate (trcrn_ini(nx,ntrcr,ncat) )
        allocate (vicen_init(nx,ncat) )
        allocate (aicen_init(nx,ncat) )
        allocate (aicen_tmp(nx,ncat) )
        allocate (vsnon_init(nx,ncat) )
        works(:,:) = c0
        works_old(:,:) = c0
        fiso_ocn(:) = c0
        swild(:) = c0
        swild1(:) = c0
        swild2(:) = c0
        totalicea = 0
        totaliceh = 0
        do i=1,nx_nh
         totalicea = totalicea + aice(i)*lump_ice_matrix(i)
         !totaliceh = totaliceh + aice(i)*vice(i)*lump_ice_matrix(i)
         totaliceh = totaliceh + vice(i)*lump_ice_matrix(i)
        enddo
        !write(12,*) 'before advec. ice ', totalicea,totaliceh
        call mpi_allreduce(totalicea,total_a,1,rtype,MPI_SUM,comm,ierr)
        call mpi_allreduce(totaliceh,total_h,1,rtype,MPI_SUM,comm,ierr)
        if(myrank == 0) write(16,*) 'before advec. ice ',total_a,total_h
        trcrn_ini=trcrn
        vicen_init=vicen
        aicen_init=aicen
        vsnon_init=vsnon

        call state_to_work (ntrcr, narr, works(:,:))
         works_old=works
        ! Advect each tracer
      
        do nt = 1, narr    
              call fct_solve_icepack (works(:,nt) )
        end do
        !call fct_solve_icepack (works(:,1) )
        !narrays = 1
        ! do nt = 1, ncat
        !    call fct_solve_icepack (works(:,narrays+1) ) 
        !    call fct_solve_icepack (works(:,narrays+2) )  
        !    call fct_solve_icepack (works(:,narrays+3) )  
        !    narrays = narrays + 3 + ntrcr
        ! enddo

        !do nt = 1, narr    
        !      call upwind_icepack (works(:,nt) )
        !end do
        
        newwork(:) = c0
        call work_to_state (ntrcr, narr, works(:,:))

        iniflag(:) = c0
         do nt = 1, ncat
            do it = 1, ntrcr

               do i=1,nx_nh
                  trcmin = minval(trcrn_ini(indnd(1:nnp(i),i),it,nt))
                  trcmin = min(trcmin,trcrn_ini(i,it,nt))
                  trcmax = maxval(trcrn_ini(indnd(1:nnp(i),i),it,nt))
                  trcmax = max(trcmax,trcrn_ini(i,it,nt))
                  if(trcrn(i,it,nt) > trcmax) iniflag(i) = 1 !trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmax
                  if(trcrn(i,it,nt) < trcmin) iniflag(i) = 1 !trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmin
               enddo
            enddo
         enddo

         call exchange_p2d(iniflag)
         
         do nt = 1, ncat
            do it = 1, ntrcr

               do i=1,nx
                  if(iniflag(i) > 0) trcrn(i,it,nt) = trcrn_ini(i,it,nt)
               enddo

               !swild=trcrn(:,it,nt)
               !call exchange_p2d(swild)
               !trcrn(:,it,nt)=swild
               
            enddo
         enddo
        ! cut off icepack
      
         !call cut_off_icepack
       
        do i=1,nx
           if (ncat < 0) then ! Do we really need this?
              call cleanup_itd  (dt*nstep_ice,           ntrcr,                &
                                 nilyr,                  nslyr,                &
                                 ncat,                   hin_max(:),           &
                                 aicen(i,:),             trcrn(i,1:ntrcr,:),   &
                                 vicen(i,:),             vsnon(i,:),           &
                                 aice0(i),               aice(i),              &
                                 n_aero,                                       &
                                 nbtrcr,                 nblyr,                &
                                 tr_aero,                                      &
                                 tr_pond_topo,                                 &
                                 first_ice(i,:),                               &
                                 trcr_depend(1:ntrcr),   trcr_base(1:ntrcr,:), &
                                 n_trcr_strata(1:ntrcr), nt_strata(1:ntrcr,:), &
                                 fpond(i),               fresh(i),             &
                                 fsalt(i),               fhocn(i),             &
                                 faero_ocn(i,:),         fiso_ocn,          &
                                 flux_bio(i,1:nbtrcr))
      
              call icepack_aggregate (ncat,                    &
                                     aicen(i,:),               &
                                     trcrn(i,1:ntrcr,:),       &
                                     vicen(i,:),               &
                                     vsnon(i,:),               &
                                     aice (i),                 &
                                     trcr (i,1:ntrcr),         &
                                     vice (i),                 &
                                     vsno (i),                 &
                                     aice0(i),                 &
                                     ntrcr,                    &
                                     trcr_depend  (1:ntrcr),   &
                                     trcr_base    (1:ntrcr,:), &
                                     n_trcr_strata(1:ntrcr),   &
                                     nt_strata    (1:ntrcr,:))
           endif
        end do
        totalicea = 0
        totaliceh = 0
        do i=1,nx_nh
         totalicea = totalicea + sum(aicen(i,:))*lump_ice_matrix(i)
         !totaliceh = totaliceh + aice(i)*vice(i)*lump_ice_matrix(i)
         totaliceh = totaliceh + sum(vicen(i,:))*lump_ice_matrix(i)
        enddo
        !write(12,*) 'after advec. ice ', totalicea,totaliceh
        call mpi_allreduce(totalicea,total_a,1,rtype,MPI_SUM,comm,ierr)
        call mpi_allreduce(totaliceh,total_h,1,rtype,MPI_SUM,comm,ierr)
        if(myrank == 0) write(16,*) 'after advec. ice ',total_a,total_h
        deallocate(works)
        deallocate(works_old)
        deallocate(trcrn_ini)
        deallocate(vicen_init)
        deallocate(aicen_init)
        deallocate(vsnon_init)
        deallocate(swild)
        deallocate(swild1)
        deallocate(swild2)

    end subroutine tracer_advection_icepack

    !=======================================================================
    module subroutine tracer_advection_icepack2

        !use mod_mesh
        use mice_module
        use icepack_intfc,        only: icepack_aggregate
        use icepack_itd,          only: cleanup_itd
        use schism_glbl,          only: dt,nstep_ice
        !use g_config,             only: dt
        use icepack_intfc,         only: icepack_init_trcr
        implicit none
            
        integer (kind=int_kind) :: ntrcr, ntrace, narr, nbtrcr, i, j,     &
                                   nt,   nt1,    k,   n ,ktherm,narrays,it,ierr
        integer (kind=int_kind) :: nt_Tsfc, nt_qice, nt_qsno,                          & 
                                   nt_sice, nt_fbri, nt_iage, nt_FY, nt_alvl, nt_vlvl, &
                                   nt_apnd, nt_hpnd, nt_ipnd, nt_bgc_Nit, nt_bgc_S
        logical (kind=log_kind) :: tr_pond_topo, tr_pond_lvl,                          &
                                   tr_pond,      tr_aero,     tr_FY,                   &
                                   tr_iage,      flag_old,                             &
                                   flag_old_0
        real    (kind=dbl_kind) :: puny ,Tsfc ,rhos, Lfresh, jj, kk, njj, pi, trcmin, trcmax,totalicea,totaliceh, &
        total_a,total_h, total_e, totalenthi
      
        ! Tracer dependencies and additional arrays
         real (kind=dbl_kind), dimension(ncat) :: &
               tmp, exc
        integer (kind=int_kind), dimension(:),    allocatable    ::    &
                tracer_type    , & ! = 1, 2, or 3 (depends on 0, 1 or 2 other tracers)
                depend             ! tracer dependencies (see below)
      
        logical (kind=log_kind), dimension (:),   allocatable   ::     &
                has_dependents    ! true if a tracer has dependent tracers
      
        real (kind=dbl_kind),    dimension (:,:), allocatable ::       &
               works,works_old, vicen_init, aicen_init, vsnon_init, aicen_tmp, iniflag, hicen_init, &
               &hsnon_init
        real (kind=dbl_kind),    dimension (:), allocatable ::       &
               fiso_ocn,swild,swild1,swild2
         real (kind=dbl_kind),    dimension (:,:,:), allocatable ::       &
         trcrn_ini, temp_matrix0

         real (kind=dbl_kind), dimension(nilyr) :: &
           qin            , & ! ice enthalpy (J/m3)
           qin_max        , & ! maximum ice enthalpy (J/m3)
           zTin               ! initial ice temperature
  
        real (kind=dbl_kind), dimension(nslyr) :: &
           qsn            , & ! snow enthalpy (J/m3)
           zTsn               ! initial snow temperature

        !type(t_mesh),        target,       intent(in)    :: mesh
           character(len=*), parameter :: subname = '(tracer_advection_icepack2)'

        call icepack_query_parameters(ktherm_out=ktherm,    &
                                      puny_out=puny,rhos_out=rhos,                   &
                                      Lfresh_out=Lfresh  ,pi_out=pi)
        call icepack_query_tracer_sizes(ntrcr_out=ntrcr, nbtrcr_out=nbtrcr)
        call icepack_query_tracer_flags(                                  &
               tr_iage_out=tr_iage, tr_FY_out=tr_FY,                      &
               tr_aero_out=tr_aero, tr_pond_out=tr_pond,                  &
               tr_pond_lvl_out=tr_pond_lvl, tr_pond_topo_out=tr_pond_topo)
        call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_qice_out=nt_qice, &
               nt_qsno_out=nt_qsno, nt_sice_out=nt_sice, nt_fbri_out=nt_fbri,       &
               nt_iage_out=nt_iage, nt_FY_out=nt_FY, nt_alvl_out=nt_alvl,           &
               nt_vlvl_out=nt_vlvl, nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd,       &
               nt_ipnd_out=nt_ipnd, nt_bgc_Nit_out=nt_bgc_Nit, nt_bgc_S_out=nt_bgc_S)
      
        narr   = 1 + ncat * (3 + ntrcr) ! max number of state variable arrays
      
        ! Allocate works array
        flag_old = .false.
        flag_old_0 = .false.

        if (allocated(works)) deallocate(works)
        if (allocated(works_old)) deallocate(works_old)
        allocate ( works(nx,narr) ) !(npa, max number of state variable arrays)
        allocate ( works_old(nx,narr) )
        allocate ( fiso_ocn(nx) )
        allocate ( swild(nx) )
        allocate ( swild1(nx) )
        allocate ( swild2(nx) )
        allocate ( iniflag(nx,ncat) )
        allocate (trcrn_ini(nx,ntrcr,ncat) )
        allocate (vicen_init(nx,ncat) )
        allocate (aicen_init(nx,ncat) )
        allocate (hicen_init(nx,ncat) )
        allocate (hsnon_init(nx,ncat) )
        allocate (aicen_tmp(nx,ncat) )
        allocate (vsnon_init(nx,ncat) )
        allocate (temp_matrix0(3,mnei_p,ncat))  
        works(:,:) = c0
        works_old(:,:) = c0
        fiso_ocn(:) = c0
        swild(:) = c0
        swild1(:) = c0
        swild2(:) = c0
        temp_matrix0 = c0

        totalicea = 0
        totaliceh = 0
        totalenthi = 0
        do i=1,nx_nh
         if(associated(ipgl(iplg(i))%next)) then !interface node
            if(ipgl(iplg(i))%next%rank<myrank) cycle !already in the sum so skip
          endif
       totalicea = totalicea + sum(aicen(i,:))*lump_ice_matrix(i)
       !totaliceh = totaliceh + aice(i)*vice(i)*lump_ice_matrix(i)
       totaliceh = totaliceh + sum(vicen(i,:))*lump_ice_matrix(i)
        enddo
        call mpi_allreduce(totalicea,total_a,1,rtype,MPI_SUM,comm,ierr)
        call mpi_allreduce(totaliceh,total_h,1,rtype,MPI_SUM,comm,ierr)
        if(myrank == 0) write(16,*) 'before advec. ice ',total_a,total_h

        trcrn_ini=trcrn
        vicen_init=vicen
        aicen_init=aicen
        vsnon_init=vsnon
        hicen_init(:,:) = c0
        hsnon_init(:,:) = c0
        call state_to_work (ntrcr, narr, works(:,:))
         works_old=works
        ! Advect each tracer
      
        !do nt = 1, narr    
        !      call upwind_icepack (works(:,nt) )
        !end do
!upwind2
        call upwind_icepack2 (works(:,1),1 )
        narrays = 1
         do nt = 1, ncat
            call upwind_icepack2 (works(:,narrays+1),narrays+1 ) 
            aicen_tmp(:,nt) = works(:,narrays+1)
            narrays = narrays + 3 + ntrcr
         enddo
         
         
         !narrays = 1
         !do nt = 1, ncat
         !   totalicea = 0
         !   do i=1,nx_nh
         !      totalicea = totalicea + sum(flux_matrix0(:,:,i,narrays+1))*lump_ice_matrix(i)
         !   enddo
         !   narrays = narrays + 3 + ntrcr
         !   call mpi_allreduce(totalicea,total_a,1,rtype,MPI_SUM,comm,ierr)
         !   if(myrank == 0) write(16,*) 'flux_matrix0 advec. ice ',nt,total_a
         !enddo

         !do i = 1, nx  ! For each grid cell
         !  if (sum(aicen_tmp(i,:)) > c1) then
         !     tmp(:) = c0
         !     exc(:) = c0
         !     do n = 1, ncat
         !        if (aicen_tmp(i,n) > puny) tmp(n) = c1
         !     end do
         !     do n = 1, ncat
         !         exc(n) = max(c0,(sum(aicen_tmp(i,:)) - c1 + puny))  &
         !                  * aicen_tmp(i,n) / sum(aicen_tmp(i,:))
         !     end do
         !     do n = 1, ncat
         !        aicen_tmp(i,n) = max(c0,aicen_tmp(i,n) - exc(n))
         !     end do
         !        aice0      = max(c0,1-sum(aicen_tmp(i,:)))
         !        if (sum(aicen_tmp(i,:)) > c1) write(12,*) 'here1',sum(aicen_tmp(i,:)),aicen_tmp(i,:)
         !  end if
         !end do

         !do i = 1, nx_nh
         !   if(sum(aicen_tmp(i,:)) > c1) then
         !      temp_matrix0 = c0
         !      narrays = 1
         !      do nt = 1, ncat
         !         do j=1,nne(i)
         !            do jj=1,3 
         !               if( flux_matrix0(jj,j,i,narrays+1)>0 ) temp_matrix0(jj,j,nt) = flux_matrix0(jj,j,i,narrays+1)
         !            enddo
         !         enddo
         !         narrays = narrays + 3 + ntrcr
         !      enddo
         !
         !      narrays = 1
         !      do nt = 1, ncat
         !         do j=1,nne(i)
         !            do jj=1,3 
         !               if( flux_matrix0(jj,j,i,narrays+1)>0 ) then
         !                  works(i,narrays+1) =  works(i,narrays+1) - flux_matrix0(jj,j,i,narrays+1)
         !                  flux_matrix0(jj,j,i,narrays+1) = flux_matrix0(jj,j,i,narrays+1) * (sum(temp_matrix0) + 1 - sum(aicen_tmp(i,:))) / sum(temp_matrix0)
         !                  works(i,narrays+1) =  works(i,narrays+1) + flux_matrix0(jj,j,i,narrays+1)
         !               endif
         !            enddo
         !         enddo
         !         narrays = narrays + 3 + ntrcr
         !      enddo
         !   endif
         !enddo
         !narrays = 1
         !do nt = 1, ncat
         !   swild = works(:,narrays+1)
         !   call exchange_p2d(swild)
         !   works(:,narrays+1) = swild
         !   narrays = narrays + 3 + ntrcr
         !enddo
         narrays = 1
         do nt = 1, ncat
            do i = 1, nx

               if(aicen_init(i,nt)<puny.or.vsnon_init(i,nt)<puny) then
                  swild2(i) = 0
                  hsnon_init(i,nt) = 0
               else
                  swild2(i) = vsnon_init(i,nt)/aicen_init(i,nt)
                  hsnon_init(i,nt) = vsnon_init(i,nt)/aicen_init(i,nt)
               endif

               if(aicen_init(i,nt)<puny.or.vicen_init(i,nt)<puny) then
                  swild1(i) = 0
                  swild2(i) = 0
                  hicen_init(i,nt) = 0
                  hsnon_init(i,nt) = 0
               else
                  swild1(i) = vicen_init(i,nt)/aicen_init(i,nt)
                  hicen_init(i,nt) = vicen_init(i,nt)/aicen_init(i,nt)
               endif


            enddo
            call upwind_icepack_dep2(works(:,narrays+2),swild1,narrays+2,narrays+1)
            call upwind_icepack_dep2(works(:,narrays+3),swild2,narrays+3,narrays+1) 
            narrays = narrays + 3 + ntrcr  
         enddo
         
         
        call work_to_other2  (ntrcr, narr, works(:,:)) 
        call re_momentum2 (ntrcr, narr, works(:,:), vicen_init)

        newwork(:) = c0

        call work_to_state (ntrcr, narr, works(:,:))
        iniflag(:,:) = c0
        swild1(:) = c0
        swild2(:) = c0
         do nt = 1, ncat
            do i=1,nx_nh
               do it = 1, ntrcr
                  trcmin = minval(trcrn_ini(indnd(1:nnp(i),i),it,nt))
                  trcmin = min(trcmin,trcrn_ini(i,it,nt))
                  trcmax = maxval(trcrn_ini(indnd(1:nnp(i),i),it,nt))
                  trcmax = max(trcmax,trcrn_ini(i,it,nt))
                  if(trcrn(i,it,nt) > trcmax) iniflag(i,nt) = 1 !trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmax
                  if(trcrn(i,it,nt) < trcmin) iniflag(i,nt) = 1 !trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmin
                  !if(trcrn(i,it,nt) > trcmax) trcrn(i,it,nt) = trcmax !trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmax
                  !if(trcrn(i,it,nt) < trcmin) trcrn(i,it,nt) = trcmin !trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmin
               enddo

                  trcmax = maxval(hicen_init(indnd(1:nnp(i),i),nt))
                  trcmax = max(trcmax,hicen_init(i,nt))
                  if(aicen(i,nt)<puny.or.vicen(i,nt)<puny) then
                     swild1(i) = 0
                  else
                     swild1(i) = vicen(i,nt)/aicen(i,nt)
                  endif
                  ! reduce hicen by increase aicen
                  if(swild1(i) > trcmax .and. aicen(i,nt)<0.15 ) then
                     aicen(i,nt) = vicen(i,nt) / min(trcmax,hin_max(nt))
                     !vicen(i,nt) = aicen(i,nt) * hicen_init(i,nt)
                  endif

                  trcmax = maxval(hsnon_init(indnd(1:nnp(i),i),nt))
                  trcmax = max(trcmax,hsnon_init(i,nt))
                  if(aicen(i,nt)<puny.or.vsnon(i,nt)<puny) then
                     swild2(i) = 0
                  else
                     swild2(i) = vsnon(i,nt)/aicen(i,nt)
                  endif
                  ! reduce hsno
                  if(swild2(i) > trcmax .and. aicen(i,nt)<0.15 ) then
                     vsnon(i,nt) = trcmax * aicen(i,nt)
                     !vsnon(i,nt) = hsnon_init(i,nt) * aicen(i,nt)
                  endif


            enddo
         enddo
         do nt = 1, ncat
            swild = iniflag(:,nt)
            call exchange_p2d(swild)
            iniflag(:,nt) = swild

            swild = aicen(:,nt)
            call exchange_p2d(swild)
            aicen(:,nt) = swild

            swild = vicen(:,nt)
            call exchange_p2d(swild)
            vicen(:,nt) = swild

            swild = vsnon(:,nt)
            call exchange_p2d(swild)
            vsnon(:,nt) = swild
         enddo
         do nt = 1, ncat
               do i=1,nx
                  if(iniflag(i,nt) > 0) then
                     trcrn(i,:,nt) = trcrn_ini(i,:,nt)
                  endif
               enddo
          enddo

        ! cut off icepack
      
        !call cut_off_icepack
       
        do i=1,nx
           if (ncat < 0) then ! Do we really need this?
              call cleanup_itd  (dt*nstep_ice,           ntrcr,                &
                                 nilyr,                  nslyr,                &
                                 ncat,                   hin_max(:),           &
                                 aicen(i,:),             trcrn(i,1:ntrcr,:),   &
                                 vicen(i,:),             vsnon(i,:),           &
                                 aice0(i),               aice(i),              &
                                 n_aero,                                       &
                                 nbtrcr,                 nblyr,                &
                                 tr_aero,                                      &
                                 tr_pond_topo,                                 &
                                 first_ice(i,:),                               &
                                 trcr_depend(1:ntrcr),   trcr_base(1:ntrcr,:), &
                                 n_trcr_strata(1:ntrcr), nt_strata(1:ntrcr,:), &
                                 fpond(i),               fresh(i),             &
                                 fsalt(i),               fhocn(i),             &
                                 faero_ocn(i,:),         fiso_ocn,          &
                                 flux_bio(i,1:nbtrcr))
      
              call icepack_aggregate (ncat,                    &
                                     aicen(i,:),               &
                                     trcrn(i,1:ntrcr,:),       &
                                     vicen(i,:),               &
                                     vsnon(i,:),               &
                                     aice (i),                 &
                                     trcr (i,1:ntrcr),         &
                                     vice (i),                 &
                                     vsno (i),                 &
                                     aice0(i),                 &
                                     ntrcr,                    &
                                     trcr_depend  (1:ntrcr),   &
                                     trcr_base    (1:ntrcr,:), &
                                     n_trcr_strata(1:ntrcr),   &
                                     nt_strata    (1:ntrcr,:))
           endif
        end do
        totalicea = 0
        totaliceh = 0
        do i=1,nx_nh
         if(associated(ipgl(iplg(i))%next)) then !interface node
            if(ipgl(iplg(i))%next%rank<myrank) cycle !already in the sum so skip
          endif
         totalicea = totalicea + sum(aicen(i,:))*lump_ice_matrix(i)
         !totaliceh = totaliceh + aice(i)*vice(i)*lump_ice_matrix(i)
         totaliceh = totaliceh + sum(vicen(i,:))*lump_ice_matrix(i)
        enddo
        call mpi_allreduce(totalicea,total_a,1,rtype,MPI_SUM,comm,ierr)
        call mpi_allreduce(totaliceh,total_h,1,rtype,MPI_SUM,comm,ierr)
        if(myrank == 0) write(16,*) 'after advec. ice ',total_a,total_h
        deallocate(works)
        deallocate(works_old)
        deallocate(trcrn_ini)
        deallocate(vicen_init)
        deallocate(aicen_init)
        deallocate(vsnon_init)
        deallocate(hicen_init)
        deallocate(hsnon_init)
        deallocate(swild)
        deallocate(swild1)
        deallocate(swild2)

    end subroutine tracer_advection_icepack2

!=======================================================================
    module subroutine tracer_advection_icepack3

      !use mod_mesh
      use mice_module
      use icepack_intfc,        only: icepack_aggregate
      use icepack_itd,          only: cleanup_itd
      use schism_glbl,          only: dt,nstep_ice
      !use g_config,             only: dt
      use icepack_intfc,         only: icepack_init_trcr
      implicit none
          
      integer (kind=int_kind) :: ntrcr, ntrace, narr, nbtrcr, i,     &
                                 nt,   nt1,    k,   n ,ktherm,narrays,it,ierr
      integer (kind=int_kind) :: nt_Tsfc, nt_qice, nt_qsno,                          & 
                                 nt_sice, nt_fbri, nt_iage, nt_FY, nt_alvl, nt_vlvl, &
                                 nt_apnd, nt_hpnd, nt_ipnd, nt_bgc_Nit, nt_bgc_S
      logical (kind=log_kind) :: tr_pond_topo, tr_pond_lvl,                          &
                                 tr_pond,      tr_aero,     tr_FY,                   &
                                 tr_iage,      flag_old,                             &        
                                 flag_old_0
      real    (kind=dbl_kind) :: puny ,Tsfc ,rhos, Lfresh, jj, kk, njj, pi, trcmin, trcmax,totalicea,totaliceh, &
      total_a,total_h,momentum_u,momentum_v, total_e, totalenthi, hilyr, ridge_threshold, total_ai, totaliceai
    
      ! Tracer dependencies and additional arrays
       real (kind=dbl_kind), dimension(ncat) :: &
             tmp, exc
      integer (kind=int_kind), dimension(:),    allocatable    ::    &
              tracer_type    , & ! = 1, 2, or 3 (depends on 0, 1 or 2 other tracers)
              depend             ! tracer dependencies (see below)
    
      logical (kind=log_kind), dimension (:),   allocatable   ::     &
              has_dependents    ! true if a tracer has dependent tracers
    
      real (kind=dbl_kind),    dimension (:,:), allocatable ::       &
             works,works_old, vicen_init, aicen_init, vsnon_init, aicen_tmp, iniflag, hicen_init, &
               &hsnon_init
      real (kind=dbl_kind),    dimension (:), allocatable ::       &
             fiso_ocn,swild,swild1,swild2, adv_threshold0
       real (kind=dbl_kind),    dimension (:,:,:), allocatable ::       &
       trcrn_ini

       real (kind=dbl_kind), dimension(nilyr) :: &
         qin            , & ! ice enthalpy (J/m3)
         qin_max        , & ! maximum ice enthalpy (J/m3)
         zTin               ! initial ice temperature

      real (kind=dbl_kind), dimension(nslyr) :: &
         qsn            , & ! snow enthalpy (J/m3)
         zTsn               ! initial snow temperature

      !type(t_mesh),        target,       intent(in)    :: mesh
         character(len=*), parameter :: subname = '(tracer_advection_icepack3)'

      call icepack_query_parameters(ktherm_out=ktherm,    &
                                    puny_out=puny,rhos_out=rhos,                   &
                                    Lfresh_out=Lfresh  ,pi_out=pi)
      call icepack_query_tracer_sizes(ntrcr_out=ntrcr, nbtrcr_out=nbtrcr)
      call icepack_query_tracer_flags(                                  &
             tr_iage_out=tr_iage, tr_FY_out=tr_FY,                      &
             tr_aero_out=tr_aero, tr_pond_out=tr_pond,                  &
             tr_pond_lvl_out=tr_pond_lvl, tr_pond_topo_out=tr_pond_topo)
      call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_qice_out=nt_qice, &
             nt_qsno_out=nt_qsno, nt_sice_out=nt_sice, nt_fbri_out=nt_fbri,       &
             nt_iage_out=nt_iage, nt_FY_out=nt_FY, nt_alvl_out=nt_alvl,           &
             nt_vlvl_out=nt_vlvl, nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd,       &
             nt_ipnd_out=nt_ipnd, nt_bgc_Nit_out=nt_bgc_Nit, nt_bgc_S_out=nt_bgc_S)
    
      narr   = 1 + ncat * (3 + ntrcr) ! max number of state variable arrays
    
      ! Allocate works array
      flag_old = .false.
      flag_old_0 = .false.

      if (allocated(works)) deallocate(works)
      if (allocated(works_old)) deallocate(works_old)
      allocate ( works(nx,narr) ) !(npa, max number of state variable arrays)
      allocate ( works_old(nx,narr) )
      allocate ( fiso_ocn(nx) )
      allocate ( swild(nx) )
      allocate ( swild1(nx) )
      allocate ( swild2(nx) )
      allocate ( adv_threshold0(nx) )
      allocate ( iniflag(nx,ncat) )
      allocate (trcrn_ini(nx,ntrcr,ncat) )
      allocate (vicen_init(nx,ncat) )
      allocate (aicen_init(nx,ncat) )
      allocate (hicen_init(nx,ncat) )
      allocate (hsnon_init(nx,ncat) )
      allocate (aicen_tmp(nx,ncat) )
      allocate (vsnon_init(nx,ncat) )
      works(:,:) = c0
      works_old(:,:) = c0
      fiso_ocn(:) = c0
      swild(:) = c0
      swild1(:) = c0
      swild2(:) = c0

      totalicea = 0
      totaliceh = 0
      momentum_u = 0
      momentum_v = 0
      totalenthi = 0
      totaliceai = 0
      do i=1,nx_nh
         if(associated(ipgl(iplg(i))%next)) then !interface node
            if(ipgl(iplg(i))%next%rank<myrank) cycle !already in the sum so skip
          endif
         totalicea = totalicea + sum(aicen(i,:))*lump_ice_matrix(i)
         totaliceai = max(totaliceai,sum(aicen(i,:))+aice0(i))
         !totaliceh = totaliceh + aice(i)*vice(i)*lump_ice_matrix(i)
         totaliceh = totaliceh + sum(vicen(i,:))*lump_ice_matrix(i)
         momentum_u = momentum_u + sum(vicen(i,:))*uvel(i)*lump_ice_matrix(i)
         momentum_v = momentum_v + sum(vicen(i,:))*vvel(i)*lump_ice_matrix(i)
         do n =1,ncat
            if(aicen(i,n)>puny) then
               hilyr = vicen(i,n) / aicen(i,n) /nilyr
            else
               hilyr = 0
            endif
            do k = 1,nilyr
               totalenthi = totalenthi + aicen(i,n)*hilyr*trcrn(i,nt_qice+k-1,n)*lump_ice_matrix(i)
            enddo
         enddo 

      enddo
      call mpi_allreduce(totalicea,total_a,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(totaliceh,total_h,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(totalenthi,total_e,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(totaliceai,total_ai,1,rtype,MPI_MAX,comm,ierr)
      if(myrank == 0) write(16,*) 'before advec. ice ',total_a,total_h,total_ai
      call mpi_allreduce(momentum_u,total_a,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(momentum_v,total_h,1,rtype,MPI_SUM,comm,ierr)
      if(myrank == 0) write(16,*) 'before advec. ice momentum',total_a,total_h

      trcrn_ini=trcrn
      vicen_init=vicen
      aicen_init=aicen
      vsnon_init=vsnon
      hicen_init(:,:) = c0  
      hsnon_init(:,:) = c0

      call state_to_work (ntrcr, narr, works(:,:))
       works_old=works
      ! Advect each tracer
    
      !do nt = 1, narr    
      !      call upwind_icepack (works(:,nt) )
      !end do
!upwind2
      !do i=1,nx
      !   adv_threshold0(i)= 0.01 !max(c0,min(0.05,exp(-7*aice0(i))))
      !enddo
      do i=1,nx
         if(vice(i)>puny) then
            adv_threshold0(i)= 10 !max(c0,min(0.01,exp(-10*aice(i))/vice(i)))
         else
            adv_threshold0(i)= 10
         endif
      enddo

      call upwind_icepack (works(:,1),1,adv_threshold0 )
      
      narrays = 1
       do nt = 1, ncat
          call upwind_icepack (works(:,narrays+1),narrays+1,adv_threshold0 ) 
          aicen_tmp(:,nt) = works(:,narrays+1)
          narrays = narrays + 3 + ntrcr
       enddo

      swild = works(:,1)
       do nt = 1, ncat
         do i=1,nx
          swild(i) = swild(i) + aicen_tmp(i,nt)
         enddo
       enddo

       ridge_threshold = 1/864 + 1
       !do i=1,nx
       !   if(swild(i)>1.005) write(12,*) i,swild(i), works(i,1),aicen_tmp(i,1),aicen_tmp(i,2),aicen_tmp(i,3),&
       !   &aicen_tmp(i,4),aicen_tmp(i,5),aice0(i),aicen(i,1),aicen(i,2),aicen(i,3),aicen(i,4),aicen(i,5)
          !if(swild(i)>ridge_threshold) then
          !  works(i,1) = works(i,1) / swild(i) * ridge_threshold
            !aice0(i) = aice0(i) / swild(i) * ridge_threshold
          !  narrays = 1
          !  do nt = 1, ncat
          !     works(i,narrays+1) = works(i,narrays+1) / swild(i) * ridge_threshold
          !     !aicen(i,nt) = aicen(i,nt) / swild(i) * ridge_threshold
          !     aicen_tmp(i,nt) = works(i,narrays+1)
          !     flux_matrix3(:,i,narrays+1) = flux_matrix3(:,i,narrays+1) / swild(i) * ridge_threshold
          !     flux_matrix4(:,i,narrays+1) = flux_matrix4(:,i,narrays+1) / swild(i) * ridge_threshold
          !     flux_matrix1(:,i,narrays+1) = flux_matrix1(:,i,narrays+1) / swild(i) * ridge_threshold
          !     flux_matrix2(:,i,narrays+1) = flux_matrix2(:,i,narrays+1) / swild(i) * ridge_threshold
          !     narrays = narrays + 3 + ntrcr
          !  enddo
          !endif
       !enddo



         !narrays = 1
         !do nt = 1, ncat
         !   totalicea = 0
         !   do i=1,ne
         !      totalicea = totalicea + sum((flux_matrix1(:,elnode(1:i34(i),i),narrays+1)+flux_matrix2(:,elnode(1:i34(i),i),narrays+1))/2+ &
         !      &(flux_matrix3(:,elnode(1:i34(i),i),narrays+1)+flux_matrix4(:,elnode(1:i34(i),i),narrays+1))/2)
         !   enddo
         !   narrays = narrays + 3 + ntrcr
         !   call mpi_allreduce(totalicea,total_a,1,rtype,MPI_SUM,comm,ierr)
         !   if(myrank == 0) write(16,*) 'flux_matrix advec. ice ',nt,total_a
         !enddo
       narrays = 1
       do nt = 1, ncat
          !call upwind_icepack (works(:,narrays+2),narrays+2 )  
          !call upwind_icepack (works(:,narrays+3),narrays+3 )
          do i = 1, nx
               if(aicen_init(i,nt)<puny.or.vsnon_init(i,nt)<puny) then
                  swild2(i) = 0
                  hsnon_init(i,nt) = 0
               else
                  swild2(i) = vsnon_init(i,nt)/aicen_init(i,nt)
                  hsnon_init(i,nt) = vsnon_init(i,nt)/aicen_init(i,nt)
               endif

               if(aicen_init(i,nt)<puny.or.vicen_init(i,nt)<puny) then
                  swild1(i) = 0
                  swild2(i) = 0
                  hicen_init(i,nt) = 0
                  hsnon_init(i,nt) = 0
               else
                  swild1(i) = vicen_init(i,nt)/aicen_init(i,nt)
                  hicen_init(i,nt) = vicen_init(i,nt)/aicen_init(i,nt)
               endif
          enddo
          call upwind_icepack_dep(works(:,narrays+2),swild1,narrays+2,narrays+1)
          call upwind_icepack_dep(works(:,narrays+3),swild2,narrays+3,narrays+1) 
          narrays = narrays + 3 + ntrcr
       enddo

       

      call work_to_other  (ntrcr, narr, works(:,:)) 
      
      call re_momentum (ntrcr, narr, works(:,:), vicen_init)
      newwork(:) = c0


      call work_to_state (ntrcr, narr, works(:,:))
      iniflag(:,:) = c0
      swild1(:) = c0
      swild2(:) = c0
         do nt = 1, ncat
            do i=1,nx_nh
               do it = 1, ntrcr
                  trcmin = minval(trcrn_ini(indnd(1:nnp(i),i),it,nt))
                  trcmin = min(trcmin,trcrn_ini(i,it,nt))
                  trcmax = maxval(trcrn_ini(indnd(1:nnp(i),i),it,nt))
                  trcmax = max(trcmax,trcrn_ini(i,it,nt))
                  if(trcrn(i,it,nt) > trcmax) iniflag(i,nt) = 1 !trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmax
                  if(trcrn(i,it,nt) < trcmin) iniflag(i,nt) = 1 !trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmin
                  !if(trcrn(i,it,nt) > trcmax) trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmax
                  !if(trcrn(i,it,nt) < trcmin) trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmin
                  !if(trcrn(i,it,nt) > trcmax) trcrn(i,it,nt) = trcmax !trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmax
                  !if(trcrn(i,it,nt) < trcmin) trcrn(i,it,nt) = trcmin !trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmin
               enddo

                  trcmax = maxval(hicen_init(indnd(1:nnp(i),i),nt))
                  trcmax = max(trcmax,hicen_init(i,nt))
                  if(aicen(i,nt)<puny.or.vicen(i,nt)<puny) then
                     swild1(i) = 0
                  else
                     swild1(i) = vicen(i,nt)/aicen(i,nt)
                  endif
                  ! reduce hicen by increase aicen
                  !if(swild1(i) > trcmax .and. aicen(i,nt)<0.15 ) then
                  if(swild1(i) > trcmax) then
                     !aicen(i,nt) = vicen(i,nt) / min(trcmax,hin_max(nt))
                     vicen(i,nt) = aicen(i,nt) * min(trcmax,hin_max(nt))
                  endif

                  trcmax = maxval(hsnon_init(indnd(1:nnp(i),i),nt))
                  trcmax = max(trcmax,hsnon_init(i,nt))
                  if(aicen(i,nt)<puny.or.vsnon(i,nt)<puny) then
                     swild2(i) = 0
                  else
                     swild2(i) = vsnon(i,nt)/aicen(i,nt)
                  endif
                  ! reduce hsno
                  !if(swild2(i) > trcmax .and. aicen(i,nt)<0.15 ) then
                  if(swild2(i) > trcmax) then
                     vsnon(i,nt) = trcmax * aicen(i,nt)
                     !vsnon(i,nt) = hsnon_init(i,nt) * aicen(i,nt)
                  endif


            enddo
         enddo
         do nt = 1, ncat
            swild = iniflag(:,nt)
            call exchange_p2d(swild)
            iniflag(:,nt) = swild

            swild = aicen(:,nt)
            call exchange_p2d(swild)
            aicen(:,nt) = swild

            swild = vicen(:,nt)
            call exchange_p2d(swild)
            vicen(:,nt) = swild

            swild = vsnon(:,nt)
            call exchange_p2d(swild)
            vsnon(:,nt) = swild
         enddo

       do nt = 1, ncat
             do i=1,nx
                if(iniflag(i,nt) > 0) then
                  trcrn(i,:,nt) = trcrn_ini(i,:,nt)
      !             trcrn(i,:,:) = trcrn_ini(i,:,:)
                endif
             enddo
        enddo

      swild = c0
       
       do i=1,nx
         do nt = 1, ncat
          swild(i) = swild(i) + aicen(i,nt)
         enddo
          aice0(i) = max(c0,c1-swild(i))
       enddo


      !do i=1,nx
      !   if(isbnd(1,i)/=0) then !b.c. (including open)
      !      do nt = 1, ncat
      !         aicen(i,nt) = aicen_init(i,nt)
      !         vicen(i,nt) = vicen_init(i,nt)
      !         vsnon(i,nt) = vsnon_init(i,nt)
      !         trcrn(i,:,nt) = trcrn_ini(i,:,nt)
      !      enddo
      !   endif
      !enddo
      ! cut off icepack
    
      !call cut_off_icepack

       !write(12,*) maxval(swild(:))


      do i=1,nx
         if (ncat < 0) then ! Do we really need this?
            call cleanup_itd  (dt*nstep_ice,           ntrcr,                &
                               nilyr,                  nslyr,                &
                               ncat,                   hin_max(:),           &
                               aicen(i,:),             trcrn(i,1:ntrcr,:),   &
                               vicen(i,:),             vsnon(i,:),           &
                               aice0(i),               aice(i),              &
                               n_aero,                                       &
                               nbtrcr,                 nblyr,                &
                               tr_aero,                                      &
                               tr_pond_topo,                                 &
                               first_ice(i,:),                               &
                               trcr_depend(1:ntrcr),   trcr_base(1:ntrcr,:), &
                               n_trcr_strata(1:ntrcr), nt_strata(1:ntrcr,:), &
                               fpond(i),               fresh(i),             &
                               fsalt(i),               fhocn(i),             &
                               faero_ocn(i,:),         fiso_ocn,          &
                               flux_bio(i,1:nbtrcr))
    
            call icepack_aggregate (ncat,                    &
                                   aicen(i,:),               &
                                   trcrn(i,1:ntrcr,:),       &
                                   vicen(i,:),               &
                                   vsnon(i,:),               &
                                   aice (i),                 &
                                   trcr (i,1:ntrcr),         &
                                   vice (i),                 &
                                   vsno (i),                 &
                                   aice0(i),                 &
                                   ntrcr,                    &
                                   trcr_depend  (1:ntrcr),   &
                                   trcr_base    (1:ntrcr,:), &
                                   n_trcr_strata(1:ntrcr),   &
                                   nt_strata    (1:ntrcr,:))
         endif
      end do
      totalicea = 0
      totaliceh = 0
      momentum_u = 0
      momentum_v = 0
      totalenthi = 0
      totaliceai = 0
      do i=1,nx_nh
         if(associated(ipgl(iplg(i))%next)) then !interface node
            if(ipgl(iplg(i))%next%rank<myrank) cycle !already in the sum so skip
          endif
       totalicea = totalicea + sum(aicen(i,:))*lump_ice_matrix(i)
       totaliceai = max(totaliceai,sum(aicen(i,:))+aice0(i))
       !totaliceh = totaliceh + aice(i)*vice(i)*lump_ice_matrix(i)
       totaliceh = totaliceh + sum(vicen(i,:))*lump_ice_matrix(i)
       momentum_u = momentum_u + sum(vicen(i,:))*uvel(i)*lump_ice_matrix(i)
       momentum_v = momentum_v + sum(vicen(i,:))*vvel(i)*lump_ice_matrix(i)
        do n =1,ncat
            if(aicen(i,n)>puny) then
               hilyr = vicen(i,n) / aicen(i,n) /nilyr
            else
               hilyr = 0
            endif
            do k = 1,nilyr
               totalenthi = totalenthi + aicen(i,n)*hilyr*trcrn(i,nt_qice+k-1,n)*lump_ice_matrix(i)
            enddo
         enddo 
      enddo
      call mpi_allreduce(totalicea,total_a,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(totaliceh,total_h,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(totalenthi,total_e,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(totaliceai,total_ai,1,rtype,MPI_MAX,comm,ierr)
      if(myrank == 0) write(16,*) 'after advec. ice ',total_a,total_h,total_e,total_ai
      call mpi_allreduce(momentum_u,total_a,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(momentum_v,total_h,1,rtype,MPI_SUM,comm,ierr)
      if(myrank == 0) write(16,*) 'after advec. ice momentum',total_a,total_h
      deallocate(works)
      deallocate(works_old)
      deallocate(trcrn_ini)
      deallocate(vicen_init)
      deallocate(aicen_init)
      deallocate(vsnon_init)
      deallocate(hicen_init)
      deallocate(hsnon_init)
      deallocate(swild)
      deallocate(swild1)
      deallocate(swild2)
      deallocate(adv_threshold0)
  end subroutine tracer_advection_icepack3

!=======================================================================
    module subroutine tracer_advection_icepack_cd

      !use mod_mesh
      use mice_module
      use icepack_intfc,        only: icepack_aggregate
      use icepack_itd,          only: cleanup_itd
      use schism_glbl,          only: dt,nstep_ice
      !use g_config,             only: dt
      use icepack_intfc,         only: icepack_init_trcr
      implicit none
          
      integer (kind=int_kind) :: ntrcr, ntrace, narr, nbtrcr, i,     &
                                 nt,   nt1,    k,   n ,ktherm,narrays,it,ierr
      integer (kind=int_kind) :: nt_Tsfc, nt_qice, nt_qsno,                          & 
                                 nt_sice, nt_fbri, nt_iage, nt_FY, nt_alvl, nt_vlvl, &
                                 nt_apnd, nt_hpnd, nt_ipnd, nt_bgc_Nit, nt_bgc_S
      logical (kind=log_kind) :: tr_pond_topo, tr_pond_lvl,                          &
                                 tr_pond,      tr_aero,     tr_FY,                   &
                                 tr_iage,      flag_old,                             &
                                 flag_old_0
      real    (kind=dbl_kind) :: puny ,Tsfc ,rhos, Lfresh, jj, kk, njj, pi, trcmin, trcmax,totalicea,totaliceh, &
      total_a,total_h,momentum_u,momentum_v, total_e, totalenthi, hilyr, ridge_threshold, total_ai, totaliceai
    
      ! Tracer dependencies and additional arrays
       real (kind=dbl_kind), dimension(ncat) :: &
             tmp, exc
      integer (kind=int_kind), dimension(:),    allocatable    ::    &
              tracer_type    , & ! = 1, 2, or 3 (depends on 0, 1 or 2 other tracers)
              depend             ! tracer dependencies (see below)
    
      logical (kind=log_kind), dimension (:),   allocatable   ::     &
              has_dependents    ! true if a tracer has dependent tracers
    
      real (kind=dbl_kind),    dimension (:,:), allocatable ::       &
             works,works_old, vicen_init, aicen_init, vsnon_init, aicen_tmp, iniflag, hicen_init, &
               &hsnon_init
      real (kind=dbl_kind),    dimension (:), allocatable ::       &
             fiso_ocn,swild,swild1,swild2, adv_threshold0
       real (kind=dbl_kind),    dimension (:,:,:), allocatable ::       &
       trcrn_ini

       real (kind=dbl_kind), dimension(nilyr) :: &
         qin            , & ! ice enthalpy (J/m3)
         qin_max        , & ! maximum ice enthalpy (J/m3)
         zTin               ! initial ice temperature

      real (kind=dbl_kind), dimension(nslyr) :: &
         qsn            , & ! snow enthalpy (J/m3)
         zTsn               ! initial snow temperature

      !type(t_mesh),        target,       intent(in)    :: mesh
         character(len=*), parameter :: subname = '(tracer_advection_icepack_cd)'

      call icepack_query_parameters(ktherm_out=ktherm,    &
                                    puny_out=puny,rhos_out=rhos,                   &
                                    Lfresh_out=Lfresh  ,pi_out=pi)
      call icepack_query_tracer_sizes(ntrcr_out=ntrcr, nbtrcr_out=nbtrcr)
      call icepack_query_tracer_flags(                                  &
             tr_iage_out=tr_iage, tr_FY_out=tr_FY,                      &
             tr_aero_out=tr_aero, tr_pond_out=tr_pond,                  &
             tr_pond_lvl_out=tr_pond_lvl, tr_pond_topo_out=tr_pond_topo)
      call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_qice_out=nt_qice, &
             nt_qsno_out=nt_qsno, nt_sice_out=nt_sice, nt_fbri_out=nt_fbri,       &
             nt_iage_out=nt_iage, nt_FY_out=nt_FY, nt_alvl_out=nt_alvl,           &
             nt_vlvl_out=nt_vlvl, nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd,       &
             nt_ipnd_out=nt_ipnd, nt_bgc_Nit_out=nt_bgc_Nit, nt_bgc_S_out=nt_bgc_S)
    
      narr   = 1 + ncat * (3 + ntrcr) ! max number of state variable arrays
    
      ! Allocate works array
      flag_old = .false.
      flag_old_0 = .false.

      if (allocated(works)) deallocate(works)
      if (allocated(works_old)) deallocate(works_old)
      allocate ( works(nx,narr) ) !(npa, max number of state variable arrays)
      allocate ( works_old(nx,narr) )
      allocate ( fiso_ocn(nx) )
      allocate ( swild(nx) )
      allocate ( swild1(nx) )
      allocate ( swild2(nx) )
      allocate ( adv_threshold0(nx) )
      allocate ( iniflag(nx,ncat) )
      allocate (trcrn_ini(nx,ntrcr,ncat) )
      allocate (vicen_init(nx,ncat) )
      allocate (aicen_init(nx,ncat) )
      allocate (hicen_init(nx,ncat) )
      allocate (hsnon_init(nx,ncat) )
      allocate (aicen_tmp(nx,ncat) )
      allocate (vsnon_init(nx,ncat) )
      works(:,:) = c0
      works_old(:,:) = c0
      fiso_ocn(:) = c0
      swild(:) = c0
      swild1(:) = c0
      swild2(:) = c0

      totalicea = 0
      totaliceh = 0
      momentum_u = 0
      momentum_v = 0
      totalenthi = 0
      totaliceai = 0
      do i=1,nx_nh
         if(associated(ipgl(iplg(i))%next)) then !interface node
            if(ipgl(iplg(i))%next%rank<myrank) cycle !already in the sum so skip
          endif
         totalicea = totalicea + sum(aicen(i,:))*lump_ice_matrix(i)
         totaliceai = max(totaliceai,sum(aicen(i,:))+aice0(i))
         !totaliceh = totaliceh + aice(i)*vice(i)*lump_ice_matrix(i)
         totaliceh = totaliceh + sum(vicen(i,:))*lump_ice_matrix(i)
         momentum_u = momentum_u + sum(vicen(i,:))*uvel(i)*lump_ice_matrix(i)
         momentum_v = momentum_v + sum(vicen(i,:))*vvel(i)*lump_ice_matrix(i)
         do n =1,ncat
            if(aicen(i,n)>puny) then
               hilyr = vicen(i,n) / aicen(i,n) /nilyr
            else
               hilyr = 0
            endif
            do k = 1,nilyr
               totalenthi = totalenthi + aicen(i,n)*hilyr*trcrn(i,nt_qice+k-1,n)*lump_ice_matrix(i)
            enddo
         enddo 

      enddo
      call mpi_allreduce(totalicea,total_a,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(totaliceh,total_h,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(totalenthi,total_e,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(totaliceai,total_ai,1,rtype,MPI_MAX,comm,ierr)
      if(myrank == 0) write(16,*) 'before advec. ice ',total_a,total_h,total_ai
      call mpi_allreduce(momentum_u,total_a,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(momentum_v,total_h,1,rtype,MPI_SUM,comm,ierr)
      if(myrank == 0) write(16,*) 'before advec. ice momentum',total_a,total_h

      trcrn_ini=trcrn
      vicen_init=vicen
      aicen_init=aicen
      vsnon_init=vsnon
      hicen_init(:,:) = c0  
      hsnon_init(:,:) = c0

      call state_to_work (ntrcr, narr, works(:,:))
       works_old=works
      ! Advect each tracer
    
      !do nt = 1, narr    
      !      call upwind_icepack (works(:,nt) )
      !end do
!upwind2
      !do i=1,nx
      !   adv_threshold0(i)= 0.01 !max(c0,min(0.05,exp(-7*aice0(i))))
      !enddo
      do i=1,nx
         if(vice(i)>puny) then
            adv_threshold0(i)= 10 !max(c0,min(0.01,exp(-10*aice(i))/vice(i)))
         else
            adv_threshold0(i)= 10
         endif
      enddo

      call upwind_icepack4 (works(:,1),1,adv_threshold0 )
      
      narrays = 1
       do nt = 1, ncat
          call upwind_icepack4 (works(:,narrays+1),narrays+1,adv_threshold0 ) 
          aicen_tmp(:,nt) = works(:,narrays+1)
          narrays = narrays + 3 + ntrcr
       enddo

      swild = works(:,1)
       do nt = 1, ncat
         do i=1,nx
          swild(i) = swild(i) + aicen_tmp(i,nt)
         enddo
       enddo

       ridge_threshold = 1/864 + 1
       !do i=1,nx
       !   if(swild(i)>1.005) write(12,*) i,swild(i), works(i,1),aicen_tmp(i,1),aicen_tmp(i,2),aicen_tmp(i,3),&
       !   &aicen_tmp(i,4),aicen_tmp(i,5),aice0(i),aicen(i,1),aicen(i,2),aicen(i,3),aicen(i,4),aicen(i,5)
          !if(swild(i)>ridge_threshold) then
          !  works(i,1) = works(i,1) / swild(i) * ridge_threshold
            !aice0(i) = aice0(i) / swild(i) * ridge_threshold
          !  narrays = 1
          !  do nt = 1, ncat
          !     works(i,narrays+1) = works(i,narrays+1) / swild(i) * ridge_threshold
          !     !aicen(i,nt) = aicen(i,nt) / swild(i) * ridge_threshold
          !     aicen_tmp(i,nt) = works(i,narrays+1)
          !     flux_matrix3(:,i,narrays+1) = flux_matrix3(:,i,narrays+1) / swild(i) * ridge_threshold
          !     flux_matrix4(:,i,narrays+1) = flux_matrix4(:,i,narrays+1) / swild(i) * ridge_threshold
          !     flux_matrix1(:,i,narrays+1) = flux_matrix1(:,i,narrays+1) / swild(i) * ridge_threshold
          !     flux_matrix2(:,i,narrays+1) = flux_matrix2(:,i,narrays+1) / swild(i) * ridge_threshold
          !     narrays = narrays + 3 + ntrcr
          !  enddo
          !endif
       !enddo



         !narrays = 1
         !do nt = 1, ncat
         !   totalicea = 0
         !   do i=1,ne
         !      totalicea = totalicea + sum((flux_matrix1(:,elnode(1:i34(i),i),narrays+1)+flux_matrix2(:,elnode(1:i34(i),i),narrays+1))/2+ &
         !      &(flux_matrix3(:,elnode(1:i34(i),i),narrays+1)+flux_matrix4(:,elnode(1:i34(i),i),narrays+1))/2)
         !   enddo
         !   narrays = narrays + 3 + ntrcr
         !   call mpi_allreduce(totalicea,total_a,1,rtype,MPI_SUM,comm,ierr)
         !   if(myrank == 0) write(16,*) 'flux_matrix advec. ice ',nt,total_a
         !enddo
       narrays = 1
       do nt = 1, ncat
          !call upwind_icepack (works(:,narrays+2),narrays+2 )  
          !call upwind_icepack (works(:,narrays+3),narrays+3 )
          do i = 1, nx
               if(aicen_init(i,nt)<puny.or.vsnon_init(i,nt)<puny) then
                  swild2(i) = 0
                  hsnon_init(i,nt) = 0
               else
                  swild2(i) = vsnon_init(i,nt)/aicen_init(i,nt)
                  hsnon_init(i,nt) = vsnon_init(i,nt)/aicen_init(i,nt)
               endif

               if(aicen_init(i,nt)<puny.or.vicen_init(i,nt)<puny) then
                  swild1(i) = 0
                  swild2(i) = 0
                  hicen_init(i,nt) = 0
                  hsnon_init(i,nt) = 0
               else
                  swild1(i) = vicen_init(i,nt)/aicen_init(i,nt)
                  hicen_init(i,nt) = vicen_init(i,nt)/aicen_init(i,nt)
               endif
          enddo
          call upwind_icepack_dep4(works(:,narrays+2),swild1,narrays+2,narrays+1)
          call upwind_icepack_dep4(works(:,narrays+3),swild2,narrays+3,narrays+1) 
          narrays = narrays + 3 + ntrcr
       enddo

       

      call work_to_other_4(ntrcr, narr, works(:,:)) 
      
      call re_momentum (ntrcr, narr, works(:,:), vicen_init)
      newwork(:) = c0


      call work_to_state (ntrcr, narr, works(:,:))
      iniflag(:,:) = c0
      swild1(:) = c0
      swild2(:) = c0
         do nt = 1, ncat
            do i=1,nx_nh
               do it = 1, ntrcr
                  trcmin = minval(trcrn_ini(indnd(1:nnp(i),i),it,nt))
                  trcmin = min(trcmin,trcrn_ini(i,it,nt))
                  trcmax = maxval(trcrn_ini(indnd(1:nnp(i),i),it,nt))
                  trcmax = max(trcmax,trcrn_ini(i,it,nt))
                  if(trcrn(i,it,nt) > trcmax) iniflag(i,nt) = 1 !trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmax
                  if(trcrn(i,it,nt) < trcmin) iniflag(i,nt) = 1 !trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmin
                  !if(trcrn(i,it,nt) > trcmax) trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmax
                  !if(trcrn(i,it,nt) < trcmin) trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmin
                  !if(trcrn(i,it,nt) > trcmax) trcrn(i,it,nt) = trcmax !trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmax
                  !if(trcrn(i,it,nt) < trcmin) trcrn(i,it,nt) = trcmin !trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmin
               enddo

                  trcmax = maxval(hicen_init(indnd(1:nnp(i),i),nt))
                  trcmax = max(trcmax,hicen_init(i,nt))
                  if(aicen(i,nt)<puny.or.vicen(i,nt)<puny) then
                     swild1(i) = 0
                  else
                     swild1(i) = vicen(i,nt)/aicen(i,nt)
                  endif
                  ! reduce hicen by increase aicen
                  !if(swild1(i) > trcmax .and. aicen(i,nt)<0.15 ) then
                  if(swild1(i) > trcmax) then
                     !aicen(i,nt) = vicen(i,nt) / min(trcmax,hin_max(nt))
                     vicen(i,nt) = aicen(i,nt) * min(trcmax,hin_max(nt))
                  endif

                  trcmax = maxval(hsnon_init(indnd(1:nnp(i),i),nt))
                  trcmax = max(trcmax,hsnon_init(i,nt))
                  if(aicen(i,nt)<puny.or.vsnon(i,nt)<puny) then
                     swild2(i) = 0
                  else
                     swild2(i) = vsnon(i,nt)/aicen(i,nt)
                  endif
                  ! reduce hsno
                  !if(swild2(i) > trcmax .and. aicen(i,nt)<0.15 ) then
                  if(swild2(i) > trcmax) then
                     vsnon(i,nt) = trcmax * aicen(i,nt)
                     !vsnon(i,nt) = hsnon_init(i,nt) * aicen(i,nt)
                  endif


            enddo
         enddo
         do nt = 1, ncat
            swild = iniflag(:,nt)
            call exchange_p2d(swild)
            iniflag(:,nt) = swild

            swild = aicen(:,nt)
            call exchange_p2d(swild)
            aicen(:,nt) = swild

            swild = vicen(:,nt)
            call exchange_p2d(swild)
            vicen(:,nt) = swild

            swild = vsnon(:,nt)
            call exchange_p2d(swild)
            vsnon(:,nt) = swild
         enddo

       do nt = 1, ncat
             do i=1,nx
                if(iniflag(i,nt) > 0) then
                  trcrn(i,:,nt) = trcrn_ini(i,:,nt)
      !             trcrn(i,:,:) = trcrn_ini(i,:,:)
                endif
             enddo
        enddo

      swild = c0
       
       do i=1,nx
         do nt = 1, ncat
          swild(i) = swild(i) + aicen(i,nt)
         enddo
          aice0(i) = max(c0,c1-swild(i))
       enddo


      !do i=1,nx
      !   if(isbnd(1,i)/=0) then !b.c. (including open)
      !      do nt = 1, ncat
      !         aicen(i,nt) = aicen_init(i,nt)
      !         vicen(i,nt) = vicen_init(i,nt)
      !         vsnon(i,nt) = vsnon_init(i,nt)
      !         trcrn(i,:,nt) = trcrn_ini(i,:,nt)
      !      enddo
      !   endif
      !enddo
      ! cut off icepack
    
      !call cut_off_icepack

       !write(12,*) maxval(swild(:))


      do i=1,nx
         if (ncat < 0) then ! Do we really need this?
            call cleanup_itd  (dt*nstep_ice,           ntrcr,                &
                               nilyr,                  nslyr,                &
                               ncat,                   hin_max(:),           &
                               aicen(i,:),             trcrn(i,1:ntrcr,:),   &
                               vicen(i,:),             vsnon(i,:),           &
                               aice0(i),               aice(i),              &
                               n_aero,                                       &
                               nbtrcr,                 nblyr,                &
                               tr_aero,                                      &
                               tr_pond_topo,                                 &
                               first_ice(i,:),                               &
                               trcr_depend(1:ntrcr),   trcr_base(1:ntrcr,:), &
                               n_trcr_strata(1:ntrcr), nt_strata(1:ntrcr,:), &
                               fpond(i),               fresh(i),             &
                               fsalt(i),               fhocn(i),             &
                               faero_ocn(i,:),         fiso_ocn,          &
                               flux_bio(i,1:nbtrcr))
    
            call icepack_aggregate (ncat,                    &
                                   aicen(i,:),               &
                                   trcrn(i,1:ntrcr,:),       &
                                   vicen(i,:),               &
                                   vsnon(i,:),               &
                                   aice (i),                 &
                                   trcr (i,1:ntrcr),         &
                                   vice (i),                 &
                                   vsno (i),                 &
                                   aice0(i),                 &
                                   ntrcr,                    &
                                   trcr_depend  (1:ntrcr),   &
                                   trcr_base    (1:ntrcr,:), &
                                   n_trcr_strata(1:ntrcr),   &
                                   nt_strata    (1:ntrcr,:))
         endif
      end do
      totalicea = 0
      totaliceh = 0
      momentum_u = 0
      momentum_v = 0
      totalenthi = 0
      totaliceai = 0
      do i=1,nx_nh
         if(associated(ipgl(iplg(i))%next)) then !interface node
            if(ipgl(iplg(i))%next%rank<myrank) cycle !already in the sum so skip
          endif
       totalicea = totalicea + sum(aicen(i,:))*lump_ice_matrix(i)
       totaliceai = max(totaliceai,sum(aicen(i,:))+aice0(i))
       !totaliceh = totaliceh + aice(i)*vice(i)*lump_ice_matrix(i)
       totaliceh = totaliceh + sum(vicen(i,:))*lump_ice_matrix(i)
       momentum_u = momentum_u + sum(vicen(i,:))*uvel(i)*lump_ice_matrix(i)
       momentum_v = momentum_v + sum(vicen(i,:))*vvel(i)*lump_ice_matrix(i)
        do n =1,ncat
            if(aicen(i,n)>puny) then
               hilyr = vicen(i,n) / aicen(i,n) /nilyr
            else
               hilyr = 0
            endif
            do k = 1,nilyr
               totalenthi = totalenthi + aicen(i,n)*hilyr*trcrn(i,nt_qice+k-1,n)*lump_ice_matrix(i)
            enddo
         enddo 
      enddo
      call mpi_allreduce(totalicea,total_a,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(totaliceh,total_h,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(totalenthi,total_e,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(totaliceai,total_ai,1,rtype,MPI_MAX,comm,ierr)
      if(myrank == 0) write(16,*) 'after advec. ice ',total_a,total_h,total_e,total_ai
      call mpi_allreduce(momentum_u,total_a,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(momentum_v,total_h,1,rtype,MPI_SUM,comm,ierr)
      if(myrank == 0) write(16,*) 'after advec. ice momentum',total_a,total_h
      deallocate(works)
      deallocate(works_old)
      deallocate(trcrn_ini)
      deallocate(vicen_init)
      deallocate(aicen_init)
      deallocate(vsnon_init)
      deallocate(hicen_init)
      deallocate(hsnon_init)
      deallocate(swild)
      deallocate(swild1)
      deallocate(swild2)
      deallocate(adv_threshold0)
  end subroutine tracer_advection_icepack_cd

!=======================================================================
    module subroutine tracer_advection_icepack4

      !use mod_mesh
      use mice_module
      use icepack_intfc,        only: icepack_aggregate
      use icepack_itd,          only: cleanup_itd
      use schism_glbl,          only: dt,nstep_ice
      !use g_config,             only: dt
      use icepack_intfc,         only: icepack_init_trcr
      implicit none
          
      integer (kind=int_kind) :: ntrcr, ntrace, narr, nbtrcr, i,     &
                                 nt,   nt1,    k,   n ,ktherm,narrays,it,ierr
      integer (kind=int_kind) :: nt_Tsfc, nt_qice, nt_qsno,                          & 
                                 nt_sice, nt_fbri, nt_iage, nt_FY, nt_alvl, nt_vlvl, &
                                 nt_apnd, nt_hpnd, nt_ipnd, nt_bgc_Nit, nt_bgc_S
      logical (kind=log_kind) :: tr_pond_topo, tr_pond_lvl,                          &
                                 tr_pond,      tr_aero,     tr_FY,                   &
                                 tr_iage,      flag_old,                             &
                                 flag_old_0
      real    (kind=dbl_kind) :: puny ,Tsfc ,rhos, Lfresh, jj, kk, njj, pi, trcmin, trcmax,totalicea,totaliceh, &
      total_a,total_h,momentum_u,momentum_v, total_e, totalenthi, hilyr, ridge_threshold, total_ai, totaliceai
    
      ! Tracer dependencies and additional arrays
       real (kind=dbl_kind), dimension(ncat) :: &
             tmp, exc
      integer (kind=int_kind), dimension(:),    allocatable    ::    &
              tracer_type    , & ! = 1, 2, or 3 (depends on 0, 1 or 2 other tracers)
              depend             ! tracer dependencies (see below)
    
      logical (kind=log_kind), dimension (:),   allocatable   ::     &
              has_dependents    ! true if a tracer has dependent tracers
    
      real (kind=dbl_kind),    dimension (:,:), allocatable ::       &
             works,works_old, vicen_init, aicen_init, vsnon_init, aicen_tmp, iniflag, hicen_init, &
               &hsnon_init
      real (kind=dbl_kind),    dimension (:), allocatable ::       &
             fiso_ocn,swild,swild1,swild2, adv_threshold0
       real (kind=dbl_kind),    dimension (:,:,:), allocatable ::       &
       trcrn_ini

       real (kind=dbl_kind), dimension(nilyr) :: &
         qin            , & ! ice enthalpy (J/m3)
         qin_max        , & ! maximum ice enthalpy (J/m3)
         zTin               ! initial ice temperature

      real (kind=dbl_kind), dimension(nslyr) :: &
         qsn            , & ! snow enthalpy (J/m3)
         zTsn               ! initial snow temperature

      !type(t_mesh),        target,       intent(in)    :: mesh
         character(len=*), parameter :: subname = '(tracer_advection_icepack4)'

      call icepack_query_parameters(ktherm_out=ktherm,    &
                                    puny_out=puny,rhos_out=rhos,                   &
                                    Lfresh_out=Lfresh  ,pi_out=pi)
      call icepack_query_tracer_sizes(ntrcr_out=ntrcr, nbtrcr_out=nbtrcr)
      call icepack_query_tracer_flags(                                  &
             tr_iage_out=tr_iage, tr_FY_out=tr_FY,                      &
             tr_aero_out=tr_aero, tr_pond_out=tr_pond,                  &
             tr_pond_lvl_out=tr_pond_lvl, tr_pond_topo_out=tr_pond_topo)
      call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_qice_out=nt_qice, &
             nt_qsno_out=nt_qsno, nt_sice_out=nt_sice, nt_fbri_out=nt_fbri,       &
             nt_iage_out=nt_iage, nt_FY_out=nt_FY, nt_alvl_out=nt_alvl,           &
             nt_vlvl_out=nt_vlvl, nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd,       &
             nt_ipnd_out=nt_ipnd, nt_bgc_Nit_out=nt_bgc_Nit, nt_bgc_S_out=nt_bgc_S)
    
      narr   = 1 + ncat * (3 + ntrcr) ! max number of state variable arrays
    
      ! Allocate works array
      flag_old = .false.
      flag_old_0 = .false.

      if (allocated(works)) deallocate(works)
      if (allocated(works_old)) deallocate(works_old)
      allocate ( works(nx,narr) ) !(npa, max number of state variable arrays)
      allocate ( works_old(nx,narr) )
      allocate ( fiso_ocn(nx) )
      allocate ( swild(nx) )
      allocate ( swild1(nx) )
      allocate ( swild2(nx) )
      allocate ( adv_threshold0(nx) )
      allocate ( iniflag(nx,ncat) )
      allocate (trcrn_ini(nx,ntrcr,ncat) )
      allocate (vicen_init(nx,ncat) )
      allocate (aicen_init(nx,ncat) )
      allocate (hicen_init(nx,ncat) )
      allocate (hsnon_init(nx,ncat) )
      allocate (aicen_tmp(nx,ncat) )
      allocate (vsnon_init(nx,ncat) )
      works(:,:) = c0
      works_old(:,:) = c0
      fiso_ocn(:) = c0
      swild(:) = c0
      swild1(:) = c0
      swild2(:) = c0

      totalicea = 0
      totaliceh = 0
      momentum_u = 0
      momentum_v = 0
      totalenthi = 0
      totaliceai = 0
      do i=1,nx_nh
         if(associated(ipgl(iplg(i))%next)) then !interface node
            if(ipgl(iplg(i))%next%rank<myrank) cycle !already in the sum so skip
          endif
         totalicea = totalicea + sum(aicen(i,:))*lump_ice_matrix(i)
         totaliceai = max(totaliceai,sum(aicen(i,:))+aice0(i))
         !totaliceh = totaliceh + aice(i)*vice(i)*lump_ice_matrix(i)
         totaliceh = totaliceh + sum(vicen(i,:))*lump_ice_matrix(i)
         momentum_u = momentum_u + sum(vicen(i,:))*uvel(i)*lump_ice_matrix(i)
         momentum_v = momentum_v + sum(vicen(i,:))*vvel(i)*lump_ice_matrix(i)
         do n =1,ncat
            if(aicen(i,n)>puny) then
               hilyr = vicen(i,n) / aicen(i,n) /nilyr
            else
               hilyr = 0
            endif
            do k = 1,nilyr
               totalenthi = totalenthi + aicen(i,n)*hilyr*trcrn(i,nt_qice+k-1,n)*lump_ice_matrix(i)
            enddo
         enddo 

      enddo
      call mpi_allreduce(totalicea,total_a,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(totaliceh,total_h,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(totalenthi,total_e,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(totaliceai,total_ai,1,rtype,MPI_MAX,comm,ierr)
      if(myrank == 0) write(16,*) 'before advec. ice ',total_a,total_h,total_ai
      call mpi_allreduce(momentum_u,total_a,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(momentum_v,total_h,1,rtype,MPI_SUM,comm,ierr)
      if(myrank == 0) write(16,*) 'before advec. ice momentum',total_a,total_h

      trcrn_ini=trcrn
      vicen_init=vicen
      aicen_init=aicen
      vsnon_init=vsnon
      hicen_init(:,:) = c0  
      hsnon_init(:,:) = c0

      call state_to_work (ntrcr, narr, works(:,:))
       works_old=works
      ! Advect each tracer
    
      !do nt = 1, narr    
      !      call upwind_icepack (works(:,nt) )
      !end do
!upwind2
      !do i=1,nx
      !   adv_threshold0(i)= 0.01 !max(c0,min(0.05,exp(-7*aice0(i))))
      !enddo
      do i=1,nx
         if(vice(i)>puny) then
            adv_threshold0(i)= 10 !max(c0,min(0.01,exp(-10*aice(i))/vice(i)))
         else
            adv_threshold0(i)= 10
         endif
      enddo

      call tvd_icepack (works(:,1),1,adv_threshold0 )
      
      narrays = 1
       do nt = 1, ncat
          call tvd_icepack (works(:,narrays+1),narrays+1,adv_threshold0 ) 
          aicen_tmp(:,nt) = works(:,narrays+1)
          narrays = narrays + 3 + ntrcr
       enddo

      swild = works(:,1)
       do nt = 1, ncat
         do i=1,nx
          swild(i) = swild(i) + aicen_tmp(i,nt)
         enddo
       enddo

       ridge_threshold = 1/864 + 1
       !do i=1,nx
       !   if(swild(i)>1.005) write(12,*) i,swild(i), works(i,1),aicen_tmp(i,1),aicen_tmp(i,2),aicen_tmp(i,3),&
       !   &aicen_tmp(i,4),aicen_tmp(i,5),aice0(i),aicen(i,1),aicen(i,2),aicen(i,3),aicen(i,4),aicen(i,5)
          !if(swild(i)>ridge_threshold) then
          !  works(i,1) = works(i,1) / swild(i) * ridge_threshold
            !aice0(i) = aice0(i) / swild(i) * ridge_threshold
          !  narrays = 1
          !  do nt = 1, ncat
          !     works(i,narrays+1) = works(i,narrays+1) / swild(i) * ridge_threshold
          !     !aicen(i,nt) = aicen(i,nt) / swild(i) * ridge_threshold
          !     aicen_tmp(i,nt) = works(i,narrays+1)
          !     flux_matrix3(:,i,narrays+1) = flux_matrix3(:,i,narrays+1) / swild(i) * ridge_threshold
          !     flux_matrix4(:,i,narrays+1) = flux_matrix4(:,i,narrays+1) / swild(i) * ridge_threshold
          !     flux_matrix1(:,i,narrays+1) = flux_matrix1(:,i,narrays+1) / swild(i) * ridge_threshold
          !     flux_matrix2(:,i,narrays+1) = flux_matrix2(:,i,narrays+1) / swild(i) * ridge_threshold
          !     narrays = narrays + 3 + ntrcr
          !  enddo
          !endif
       !enddo



         !narrays = 1
         !do nt = 1, ncat
         !   totalicea = 0
         !   do i=1,ne
         !      totalicea = totalicea + sum((flux_matrix1(:,elnode(1:i34(i),i),narrays+1)+flux_matrix2(:,elnode(1:i34(i),i),narrays+1))/2+ &
         !      &(flux_matrix3(:,elnode(1:i34(i),i),narrays+1)+flux_matrix4(:,elnode(1:i34(i),i),narrays+1))/2)
         !   enddo
         !   narrays = narrays + 3 + ntrcr
         !   call mpi_allreduce(totalicea,total_a,1,rtype,MPI_SUM,comm,ierr)
         !   if(myrank == 0) write(16,*) 'flux_matrix advec. ice ',nt,total_a
         !enddo
       narrays = 1
       do nt = 1, ncat
          !call upwind_icepack (works(:,narrays+2),narrays+2 )  
          !call upwind_icepack (works(:,narrays+3),narrays+3 )
          do i = 1, nx
               if(aicen_init(i,nt)<puny.or.vsnon_init(i,nt)<puny) then
                  swild2(i) = 0
                  hsnon_init(i,nt) = 0
               else
                  swild2(i) = vsnon_init(i,nt)/aicen_init(i,nt)
                  hsnon_init(i,nt) = vsnon_init(i,nt)/aicen_init(i,nt)
               endif

               if(aicen_init(i,nt)<puny.or.vicen_init(i,nt)<puny) then
                  swild1(i) = 0
                  swild2(i) = 0
                  hicen_init(i,nt) = 0
                  hsnon_init(i,nt) = 0
               else
                  swild1(i) = vicen_init(i,nt)/aicen_init(i,nt)
                  hicen_init(i,nt) = vicen_init(i,nt)/aicen_init(i,nt)
               endif
          enddo
          call tvd_icepack_dep(works(:,narrays+2),swild1,narrays+2,narrays+1)
          call tvd_icepack_dep(works(:,narrays+3),swild2,narrays+3,narrays+1) 
          !call tvd_icepack (works(:,narrays+2),narrays+2,adv_threshold0 ) 
          !call tvd_icepack (works(:,narrays+3),narrays+3,adv_threshold0 ) 
          narrays = narrays + 3 + ntrcr
       enddo

       

      call work_to_other3  (ntrcr, narr, works(:,:)) 
      
      call re_momentum (ntrcr, narr, works(:,:), vicen_init)
      newwork(:) = c0


      call work_to_state (ntrcr, narr, works(:,:))
      iniflag(:,:) = c0
      swild1(:) = c0
      swild2(:) = c0
         do nt = 1, ncat
            do i=1,nx_nh
               do it = 1, ntrcr
                  trcmin = minval(trcrn_ini(indnd(1:nnp(i),i),it,nt))
                  trcmin = min(trcmin,trcrn_ini(i,it,nt))
                  trcmax = maxval(trcrn_ini(indnd(1:nnp(i),i),it,nt))
                  trcmax = max(trcmax,trcrn_ini(i,it,nt))
                  if(trcrn(i,it,nt) > trcmax) iniflag(i,nt) = 1 !trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmax
                  if(trcrn(i,it,nt) < trcmin) iniflag(i,nt) = 1 !trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmin
                  !if(trcrn(i,it,nt) > trcmax) trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmax
                  !if(trcrn(i,it,nt) < trcmin) trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmin
                  !if(trcrn(i,it,nt) > trcmax) trcrn(i,it,nt) = trcmax !trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmax
                  !if(trcrn(i,it,nt) < trcmin) trcrn(i,it,nt) = trcmin !trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmin
               enddo

                  trcmax = maxval(hicen_init(indnd(1:nnp(i),i),nt))
                  trcmax = max(trcmax,hicen_init(i,nt))
                  if(aicen(i,nt)<puny.or.vicen(i,nt)<puny) then
                     swild1(i) = 0
                  else
                     swild1(i) = vicen(i,nt)/aicen(i,nt)
                  endif
                  ! reduce hicen by increase aicen
                  !if(swild1(i) > trcmax .and. aicen(i,nt)<0.15 ) then
                  if(swild1(i) > trcmax) then
                     !aicen(i,nt) = vicen(i,nt) / min(trcmax,hin_max(nt))
                     vicen(i,nt) = aicen(i,nt) * min(trcmax,hin_max(nt))
                  endif

                  trcmax = maxval(hsnon_init(indnd(1:nnp(i),i),nt))
                  trcmax = max(trcmax,hsnon_init(i,nt))
                  if(aicen(i,nt)<puny.or.vsnon(i,nt)<puny) then
                     swild2(i) = 0
                  else
                     swild2(i) = vsnon(i,nt)/aicen(i,nt)
                  endif
                  ! reduce hsno
                  !if(swild2(i) > trcmax .and. aicen(i,nt)<0.15 ) then
                  if(swild2(i) > trcmax) then
                     vsnon(i,nt) = trcmax * aicen(i,nt)
                     !vsnon(i,nt) = hsnon_init(i,nt) * aicen(i,nt)
                  endif


            enddo
         enddo
         do nt = 1, ncat
            swild = iniflag(:,nt)
            call exchange_p2d(swild)
            iniflag(:,nt) = swild

            swild = aicen(:,nt)
            call exchange_p2d(swild)
            aicen(:,nt) = swild

            swild = vicen(:,nt)
            call exchange_p2d(swild)
            vicen(:,nt) = swild

            swild = vsnon(:,nt)
            call exchange_p2d(swild)
            vsnon(:,nt) = swild
         enddo

       do nt = 1, ncat
             do i=1,nx
                if(iniflag(i,nt) > 0) then
                  trcrn(i,:,nt) = trcrn_ini(i,:,nt)
      !             trcrn(i,:,:) = trcrn_ini(i,:,:)
                endif
             enddo
        enddo

      swild = c0
       
       do i=1,nx
         do nt = 1, ncat
          swild(i) = swild(i) + aicen(i,nt)
         enddo
          aice0(i) = max(c0,c1-swild(i))
       enddo


      !do i=1,nx
      !   if(isbnd(1,i)/=0) then !b.c. (including open)
      !      do nt = 1, ncat
      !         aicen(i,nt) = aicen_init(i,nt)
      !         vicen(i,nt) = vicen_init(i,nt)
      !         vsnon(i,nt) = vsnon_init(i,nt)
      !         trcrn(i,:,nt) = trcrn_ini(i,:,nt)
      !      enddo
      !   endif
      !enddo
      ! cut off icepack
    
      !call cut_off_icepack

       !write(12,*) maxval(swild(:))


      do i=1,nx
         if (ncat < 0) then ! Do we really need this?
            call cleanup_itd  (dt*nstep_ice,           ntrcr,                &
                               nilyr,                  nslyr,                &
                               ncat,                   hin_max(:),           &
                               aicen(i,:),             trcrn(i,1:ntrcr,:),   &
                               vicen(i,:),             vsnon(i,:),           &
                               aice0(i),               aice(i),              &
                               n_aero,                                       &
                               nbtrcr,                 nblyr,                &
                               tr_aero,                                      &
                               tr_pond_topo,                                 &
                               first_ice(i,:),                               &
                               trcr_depend(1:ntrcr),   trcr_base(1:ntrcr,:), &
                               n_trcr_strata(1:ntrcr), nt_strata(1:ntrcr,:), &
                               fpond(i),               fresh(i),             &
                               fsalt(i),               fhocn(i),             &
                               faero_ocn(i,:),         fiso_ocn,          &
                               flux_bio(i,1:nbtrcr))
    
            call icepack_aggregate (ncat,                    &
                                   aicen(i,:),               &
                                   trcrn(i,1:ntrcr,:),       &
                                   vicen(i,:),               &
                                   vsnon(i,:),               &
                                   aice (i),                 &
                                   trcr (i,1:ntrcr),         &
                                   vice (i),                 &
                                   vsno (i),                 &
                                   aice0(i),                 &
                                   ntrcr,                    &
                                   trcr_depend  (1:ntrcr),   &
                                   trcr_base    (1:ntrcr,:), &
                                   n_trcr_strata(1:ntrcr),   &
                                   nt_strata    (1:ntrcr,:))
         endif
      end do
      totalicea = 0
      totaliceh = 0
      momentum_u = 0
      momentum_v = 0
      totalenthi = 0
      totaliceai = 0
      do i=1,nx_nh
         if(associated(ipgl(iplg(i))%next)) then !interface node
            if(ipgl(iplg(i))%next%rank<myrank) cycle !already in the sum so skip
          endif
       totalicea = totalicea + sum(aicen(i,:))*lump_ice_matrix(i)
       totaliceai = max(totaliceai,sum(aicen(i,:))+aice0(i))
       !totaliceh = totaliceh + aice(i)*vice(i)*lump_ice_matrix(i)
       totaliceh = totaliceh + sum(vicen(i,:))*lump_ice_matrix(i)
       momentum_u = momentum_u + sum(vicen(i,:))*uvel(i)*lump_ice_matrix(i)
       momentum_v = momentum_v + sum(vicen(i,:))*vvel(i)*lump_ice_matrix(i)
        do n =1,ncat
            if(aicen(i,n)>puny) then
               hilyr = vicen(i,n) / aicen(i,n) /nilyr
               !totaliceai = max(totaliceai,vicen(i,n) / aicen(i,n))
            else
               hilyr = 0
            endif
            do k = 1,nilyr
               totalenthi = totalenthi + aicen(i,n)*hilyr*trcrn(i,nt_qice+k-1,n)*lump_ice_matrix(i)
            enddo
         enddo 
      enddo
      call mpi_allreduce(totalicea,total_a,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(totaliceh,total_h,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(totalenthi,total_e,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(totaliceai,total_ai,1,rtype,MPI_MAX,comm,ierr)
      if(myrank == 0) write(16,*) 'after advec. ice ',total_a,total_h,total_e,total_ai
      call mpi_allreduce(momentum_u,total_a,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(momentum_v,total_h,1,rtype,MPI_SUM,comm,ierr)
      if(myrank == 0) write(16,*) 'after advec. ice momentum',total_a,total_h
      deallocate(works)
      deallocate(works_old)
      deallocate(trcrn_ini)
      deallocate(vicen_init)
      deallocate(aicen_init)
      deallocate(vsnon_init)
      deallocate(hicen_init)
      deallocate(hsnon_init)
      deallocate(swild)
      deallocate(swild1)
      deallocate(swild2)
      deallocate(adv_threshold0)
  end subroutine tracer_advection_icepack4

!=======================================================================
    module subroutine tracer_advection_icepack5

      !use mod_mesh
      use mice_module
      use icepack_intfc,        only: icepack_aggregate
      use icepack_itd,          only: cleanup_itd
      use schism_glbl,          only: dt,nstep_ice
      !use g_config,             only: dt
      use icepack_intfc,         only: icepack_init_trcr
      implicit none
          
      integer (kind=int_kind) :: ntrcr, ntrace, narr, nbtrcr, i,     &
                                 nt,   nt1,    k,   n ,ktherm,narrays,it,ierr
      integer (kind=int_kind) :: nt_Tsfc, nt_qice, nt_qsno,                          & 
                                 nt_sice, nt_fbri, nt_iage, nt_FY, nt_alvl, nt_vlvl, &
                                 nt_apnd, nt_hpnd, nt_ipnd, nt_bgc_Nit, nt_bgc_S
      logical (kind=log_kind) :: tr_pond_topo, tr_pond_lvl,                           &
                                 tr_pond,      tr_aero,     tr_FY,                   &
                                 tr_iage,      flag_old,                             &
                                 flag_old_0
      real    (kind=dbl_kind) :: puny ,Tsfc ,rhos, Lfresh, jj, kk, njj, pi, trcmin, trcmax,totalicea,totaliceh, &
      total_a,total_h,momentum_u,momentum_v, total_e, totalenthi, hilyr, ridge_threshold, total_ai, totaliceai
    
      ! Tracer dependencies and additional arrays
       real (kind=dbl_kind), dimension(ncat) :: &
             tmp, exc
      integer (kind=int_kind), dimension(:),    allocatable    ::    &
              tracer_type    , & ! = 1, 2, or 3 (depends on 0, 1 or 2 other tracers)
              depend             ! tracer dependencies (see below)
    
      logical (kind=log_kind), dimension (:),   allocatable   ::     &
              has_dependents    ! true if a tracer has dependent tracers
    
      real (kind=dbl_kind),    dimension (:,:), allocatable ::       &
             works,works_old, vicen_init, aicen_init, vsnon_init, aicen_tmp, iniflag, hicen_init, &
               &hsnon_init
      real (kind=dbl_kind),    dimension (:), allocatable ::       &
             fiso_ocn,swild,swild1,swild2, adv_threshold0
       real (kind=dbl_kind),    dimension (:,:,:), allocatable ::       &
       trcrn_ini

       real (kind=dbl_kind), dimension(nilyr) :: &
         qin            , & ! ice enthalpy (J/m3)
         qin_max        , & ! maximum ice enthalpy (J/m3)
         zTin               ! initial ice temperature

      real (kind=dbl_kind), dimension(nslyr) :: &
         qsn            , & ! snow enthalpy (J/m3)
         zTsn               ! initial snow temperature

      !type(t_mesh),        target,       intent(in)    :: mesh
         character(len=*), parameter :: subname = '(tracer_advection_icepack4)'

      call icepack_query_parameters(ktherm_out=ktherm,    &
                                    puny_out=puny,rhos_out=rhos,                   &
                                    Lfresh_out=Lfresh  ,pi_out=pi)
      call icepack_query_tracer_sizes(ntrcr_out=ntrcr, nbtrcr_out=nbtrcr)
      call icepack_query_tracer_flags(                                  &
             tr_iage_out=tr_iage, tr_FY_out=tr_FY,                      &
             tr_aero_out=tr_aero, tr_pond_out=tr_pond,                  &
             tr_pond_lvl_out=tr_pond_lvl, tr_pond_topo_out=tr_pond_topo)
      call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_qice_out=nt_qice, &
             nt_qsno_out=nt_qsno, nt_sice_out=nt_sice, nt_fbri_out=nt_fbri,       &
             nt_iage_out=nt_iage, nt_FY_out=nt_FY, nt_alvl_out=nt_alvl,           &
             nt_vlvl_out=nt_vlvl, nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd,       &
             nt_ipnd_out=nt_ipnd, nt_bgc_Nit_out=nt_bgc_Nit, nt_bgc_S_out=nt_bgc_S)
    
      narr   = 1 + ncat * (3 + ntrcr) ! max number of state variable arrays
    
      ! Allocate works array
      flag_old = .false.
      flag_old_0 = .false.

      if (allocated(works)) deallocate(works)
      if (allocated(works_old)) deallocate(works_old)
      allocate ( works(nx,narr) ) !(npa, max number of state variable arrays)
      allocate ( works_old(nx,narr) )
      allocate ( fiso_ocn(nx) )
      allocate ( swild(nx) )
      allocate ( swild1(nx) )
      allocate ( swild2(nx) )
      allocate ( adv_threshold0(nx) )
      allocate ( iniflag(nx,ncat) )
      allocate (trcrn_ini(nx,ntrcr,ncat) )
      allocate (vicen_init(nx,ncat) )
      allocate (aicen_init(nx,ncat) )
      allocate (hicen_init(nx,ncat) )
      allocate (hsnon_init(nx,ncat) )
      allocate (aicen_tmp(nx,ncat) )
      allocate (vsnon_init(nx,ncat) )
      works(:,:) = c0
      works_old(:,:) = c0
      fiso_ocn(:) = c0
      swild(:) = c0
      swild1(:) = c0
      swild2(:) = c0

      totalicea = 0
      totaliceh = 0
      momentum_u = 0
      momentum_v = 0
      totalenthi = 0
      totaliceai = 0
      do i=1,nx_nh
         if(associated(ipgl(iplg(i))%next)) then !interface node
            if(ipgl(iplg(i))%next%rank<myrank) cycle !already in the sum so skip
          endif
         totalicea = totalicea + sum(aicen(i,:))*lump_ice_matrix(i)
         totaliceai = max(totaliceai,sum(aicen(i,:))+aice0(i))
         !totaliceh = totaliceh + aice(i)*vice(i)*lump_ice_matrix(i)
         totaliceh = totaliceh + sum(vicen(i,:))*lump_ice_matrix(i)
         momentum_u = momentum_u + sum(vicen(i,:))*uvel(i)*lump_ice_matrix(i)
         momentum_v = momentum_v + sum(vicen(i,:))*vvel(i)*lump_ice_matrix(i)
         do n =1,ncat
            if(aicen(i,n)>puny) then
               hilyr = vicen(i,n) / aicen(i,n) /nilyr
            else
               hilyr = 0
            endif
            do k = 1,nilyr
               totalenthi = totalenthi + aicen(i,n)*hilyr*trcrn(i,nt_qice+k-1,n)*lump_ice_matrix(i)
            enddo
         enddo 

      enddo
      call mpi_allreduce(totalicea,total_a,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(totaliceh,total_h,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(totalenthi,total_e,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(totaliceai,total_ai,1,rtype,MPI_MAX,comm,ierr)
      if(myrank == 0) write(16,*) 'before advec. ice ',total_a,total_h,total_ai
      call mpi_allreduce(momentum_u,total_a,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(momentum_v,total_h,1,rtype,MPI_SUM,comm,ierr)
      if(myrank == 0) write(16,*) 'before advec. ice momentum',total_a,total_h

      trcrn_ini=trcrn
      vicen_init=vicen
      aicen_init=aicen
      vsnon_init=vsnon
      hicen_init(:,:) = c0  
      hsnon_init(:,:) = c0

      call state_to_work (ntrcr, narr, works(:,:))
       works_old=works
      ! Advect each tracer
    
      !do nt = 1, narr    
      !      call upwind_icepack (works(:,nt) )
      !end do
!upwind2
      !do i=1,nx
      !   adv_threshold0(i)= 0.01 !max(c0,min(0.05,exp(-7*aice0(i))))
      !enddo
      do i=1,nx
         if(vice(i)>puny) then
            adv_threshold0(i)= 10 !max(c0,min(0.01,exp(-10*aice(i))/vice(i)))
         else
            adv_threshold0(i)= 10
         endif
      enddo
      
      dudicex(:,:)=c0
      dudicey(:,:)=c0
      
      call tvdu_icepack (works(:,1),1,adv_threshold0 )
      
      narrays = 1
       do nt = 1, ncat
          call tvdu_icepack (works(:,narrays+1),narrays+1,adv_threshold0) 
          aicen_tmp(:,nt) = works(:,narrays+1)
          narrays = narrays + 3 + ntrcr
       enddo

      swild = works(:,1)
       do nt = 1, ncat
         do i=1,nx
          swild(i) = swild(i) + aicen_tmp(i,nt)
         enddo
       enddo

       ridge_threshold = 1/864 + 1
       !do i=1,nx
       !   if(swild(i)>1.005) write(12,*) i,swild(i), works(i,1),aicen_tmp(i,1),aicen_tmp(i,2),aicen_tmp(i,3),&
       !   &aicen_tmp(i,4),aicen_tmp(i,5),aice0(i),aicen(i,1),aicen(i,2),aicen(i,3),aicen(i,4),aicen(i,5)
          !if(swild(i)>ridge_threshold) then
          !  works(i,1) = works(i,1) / swild(i) * ridge_threshold
            !aice0(i) = aice0(i) / swild(i) * ridge_threshold
          !  narrays = 1
          !  do nt = 1, ncat
          !     works(i,narrays+1) = works(i,narrays+1) / swild(i) * ridge_threshold
          !     !aicen(i,nt) = aicen(i,nt) / swild(i) * ridge_threshold
          !     aicen_tmp(i,nt) = works(i,narrays+1)
          !     flux_matrix3(:,i,narrays+1) = flux_matrix3(:,i,narrays+1) / swild(i) * ridge_threshold
          !     flux_matrix4(:,i,narrays+1) = flux_matrix4(:,i,narrays+1) / swild(i) * ridge_threshold
          !     flux_matrix1(:,i,narrays+1) = flux_matrix1(:,i,narrays+1) / swild(i) * ridge_threshold
          !     flux_matrix2(:,i,narrays+1) = flux_matrix2(:,i,narrays+1) / swild(i) * ridge_threshold
          !     narrays = narrays + 3 + ntrcr
          !  enddo
          !endif
       !enddo



         !narrays = 1
         !do nt = 1, ncat
         !   totalicea = 0
         !   do i=1,ne
         !      totalicea = totalicea + sum((flux_matrix1(:,elnode(1:i34(i),i),narrays+1)+flux_matrix2(:,elnode(1:i34(i),i),narrays+1))/2+ &
         !      &(flux_matrix3(:,elnode(1:i34(i),i),narrays+1)+flux_matrix4(:,elnode(1:i34(i),i),narrays+1))/2)
         !   enddo
         !   narrays = narrays + 3 + ntrcr
         !   call mpi_allreduce(totalicea,total_a,1,rtype,MPI_SUM,comm,ierr)
         !   if(myrank == 0) write(16,*) 'flux_matrix advec. ice ',nt,total_a
         !enddo
       narrays = 1
       do nt = 1, ncat
          !call upwind_icepack (works(:,narrays+2),narrays+2 )  
          !call upwind_icepack (works(:,narrays+3),narrays+3 )
          do i = 1, nx
               if(aicen_init(i,nt)<puny.or.vsnon_init(i,nt)<puny) then
                  swild2(i) = 0
                  hsnon_init(i,nt) = 0
               else
                  swild2(i) = vsnon_init(i,nt)/aicen_init(i,nt)
                  hsnon_init(i,nt) = vsnon_init(i,nt)/aicen_init(i,nt)
               endif

               if(aicen_init(i,nt)<puny.or.vicen_init(i,nt)<puny) then
                  swild1(i) = 0
                  swild2(i) = 0
                  hicen_init(i,nt) = 0
                  hsnon_init(i,nt) = 0
               else
                  swild1(i) = vicen_init(i,nt)/aicen_init(i,nt)
                  hicen_init(i,nt) = vicen_init(i,nt)/aicen_init(i,nt)
               endif
          enddo
          call tvdu_icepack_dep(works(:,narrays+2),swild1,narrays+2,narrays+1)
          call tvdu_icepack_dep(works(:,narrays+3),swild2,narrays+3,narrays+1) 
          !call tvd_icepack (works(:,narrays+2),narrays+2,adv_threshold0 ) 
          !call tvd_icepack (works(:,narrays+3),narrays+3,adv_threshold0 ) 
          narrays = narrays + 3 + ntrcr
       enddo

       

      call work_to_other_5  (ntrcr, narr, works(:,:)) 
      
      call re_momentum (ntrcr, narr, works(:,:), vicen_init)
      newwork(:) = c0


      call work_to_state (ntrcr, narr, works(:,:))
      iniflag(:,:) = c0
      swild1(:) = c0
      swild2(:) = c0
         do nt = 1, ncat
            do i=1,nx_nh
               do it = 1, ntrcr
                  trcmin = minval(trcrn_ini(indnd(1:nnp(i),i),it,nt))
                  trcmin = min(trcmin,trcrn_ini(i,it,nt))
                  trcmax = maxval(trcrn_ini(indnd(1:nnp(i),i),it,nt))
                  trcmax = max(trcmax,trcrn_ini(i,it,nt))
                  if(trcrn(i,it,nt) > trcmax) iniflag(i,nt) = 1 !trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmax
                  if(trcrn(i,it,nt) < trcmin) iniflag(i,nt) = 1 !trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmin
                  !if(trcrn(i,it,nt) > trcmax) trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmax
                  !if(trcrn(i,it,nt) < trcmin) trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmin
                  !if(trcrn(i,it,nt) > trcmax) trcrn(i,it,nt) = trcmax !trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmax
                  !if(trcrn(i,it,nt) < trcmin) trcrn(i,it,nt) = trcmin !trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmin
               enddo

                  trcmax = maxval(hicen_init(indnd(1:nnp(i),i),nt))
                  trcmax = max(trcmax,hicen_init(i,nt))
                  if(aicen(i,nt)<puny.or.vicen(i,nt)<puny) then
                     swild1(i) = 0
                  else
                     swild1(i) = vicen(i,nt)/aicen(i,nt)
                  endif
                  ! reduce hicen by increase aicen
                  !if(swild1(i) > trcmax .and. aicen(i,nt)<0.15 ) then
                  if(swild1(i) > trcmax) then
                     !aicen(i,nt) = vicen(i,nt) / min(trcmax,hin_max(nt))
                     vicen(i,nt) = aicen(i,nt) * min(trcmax,hin_max(nt))
                  endif

                  trcmax = maxval(hsnon_init(indnd(1:nnp(i),i),nt))
                  trcmax = max(trcmax,hsnon_init(i,nt))
                  if(aicen(i,nt)<puny.or.vsnon(i,nt)<puny) then
                     swild2(i) = 0
                  else
                     swild2(i) = vsnon(i,nt)/aicen(i,nt)
                  endif
                  ! reduce hsno
                  !if(swild2(i) > trcmax .and. aicen(i,nt)<0.15 ) then
                  if(swild2(i) > trcmax) then
                     vsnon(i,nt) = trcmax * aicen(i,nt)
                     !vsnon(i,nt) = hsnon_init(i,nt) * aicen(i,nt)
                  endif


            enddo
         enddo
         do nt = 1, ncat
            swild = iniflag(:,nt)
            call exchange_p2d(swild)
            iniflag(:,nt) = swild

            swild = aicen(:,nt)
            call exchange_p2d(swild)
            aicen(:,nt) = swild

            swild = vicen(:,nt)
            call exchange_p2d(swild)
            vicen(:,nt) = swild

            swild = vsnon(:,nt)
            call exchange_p2d(swild)
            vsnon(:,nt) = swild
         enddo

       do nt = 1, ncat
             do i=1,nx
                if(iniflag(i,nt) > 0) then
                  trcrn(i,:,nt) = trcrn_ini(i,:,nt)
      !             trcrn(i,:,:) = trcrn_ini(i,:,:)
                endif
             enddo
        enddo

      swild = c0
       
       do i=1,nx
         do nt = 1, ncat
          swild(i) = swild(i) + aicen(i,nt)
         enddo
          aice0(i) = max(c0,c1-swild(i))
       enddo


      !do i=1,nx
      !   if(isbnd(1,i)/=0) then !b.c. (including open)
      !      do nt = 1, ncat
      !         aicen(i,nt) = aicen_init(i,nt)
      !         vicen(i,nt) = vicen_init(i,nt)
      !         vsnon(i,nt) = vsnon_init(i,nt)
      !         trcrn(i,:,nt) = trcrn_ini(i,:,nt)
      !      enddo
      !   endif
      !enddo
      ! cut off icepack
    
      !call cut_off_icepack

       !write(12,*) maxval(swild(:))


      do i=1,nx
         if (ncat < 0) then ! Do we really need this?
            call cleanup_itd  (dt*nstep_ice,           ntrcr,                &
                               nilyr,                  nslyr,                &
                               ncat,                   hin_max(:),           &
                               aicen(i,:),             trcrn(i,1:ntrcr,:),   &
                               vicen(i,:),             vsnon(i,:),           &
                               aice0(i),               aice(i),              &
                               n_aero,                                       &
                               nbtrcr,                 nblyr,                &
                               tr_aero,                                      &
                               tr_pond_topo,                                 &
                               first_ice(i,:),                               &
                               trcr_depend(1:ntrcr),   trcr_base(1:ntrcr,:), &
                               n_trcr_strata(1:ntrcr), nt_strata(1:ntrcr,:), &
                               fpond(i),               fresh(i),             &
                               fsalt(i),               fhocn(i),             &
                               faero_ocn(i,:),         fiso_ocn,          &
                               flux_bio(i,1:nbtrcr))
    
            call icepack_aggregate (ncat,                    &
                                   aicen(i,:),               &
                                   trcrn(i,1:ntrcr,:),       &
                                   vicen(i,:),               &
                                   vsnon(i,:),               &
                                   aice (i),                 &
                                   trcr (i,1:ntrcr),         &
                                   vice (i),                 &
                                   vsno (i),                 &
                                   aice0(i),                 &
                                   ntrcr,                    &
                                   trcr_depend  (1:ntrcr),   &
                                   trcr_base    (1:ntrcr,:), &
                                   n_trcr_strata(1:ntrcr),   &
                                   nt_strata    (1:ntrcr,:))
         endif
      end do
      totalicea = 0
      totaliceh = 0
      momentum_u = 0
      momentum_v = 0
      totalenthi = 0
      totaliceai = 0
      do i=1,nx_nh
         if(associated(ipgl(iplg(i))%next)) then !interface node
            if(ipgl(iplg(i))%next%rank<myrank) cycle !already in the sum so skip
          endif
       totalicea = totalicea + sum(aicen(i,:))*lump_ice_matrix(i)
       totaliceai = max(totaliceai,sum(aicen(i,:))+aice0(i))
       !totaliceh = totaliceh + aice(i)*vice(i)*lump_ice_matrix(i)
       totaliceh = totaliceh + sum(vicen(i,:))*lump_ice_matrix(i)
       momentum_u = momentum_u + sum(vicen(i,:))*uvel(i)*lump_ice_matrix(i)
       momentum_v = momentum_v + sum(vicen(i,:))*vvel(i)*lump_ice_matrix(i)
        do n =1,ncat
            if(aicen(i,n)>puny) then
               hilyr = vicen(i,n) / aicen(i,n) /nilyr
               !totaliceai = max(totaliceai,vicen(i,n) / aicen(i,n))
            else
               hilyr = 0
            endif
            do k = 1,nilyr
               totalenthi = totalenthi + aicen(i,n)*hilyr*trcrn(i,nt_qice+k-1,n)*lump_ice_matrix(i)
            enddo
         enddo 
      enddo
      call mpi_allreduce(totalicea,total_a,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(totaliceh,total_h,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(totalenthi,total_e,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(totaliceai,total_ai,1,rtype,MPI_MAX,comm,ierr)
      if(myrank == 0) write(16,*) 'after advec. ice ',total_a,total_h,total_e,total_ai
      call mpi_allreduce(momentum_u,total_a,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(momentum_v,total_h,1,rtype,MPI_SUM,comm,ierr)
      if(myrank == 0) write(16,*) 'after advec. ice momentum',total_a,total_h
      deallocate(works)
      deallocate(works_old)
      deallocate(trcrn_ini)
      deallocate(vicen_init)
      deallocate(aicen_init)
      deallocate(vsnon_init)
      deallocate(hicen_init)
      deallocate(hsnon_init)
      deallocate(swild)
      deallocate(swild1)
      deallocate(swild2)
      deallocate(adv_threshold0)
  end subroutine tracer_advection_icepack5

!=======================================================================
    module subroutine tracer_advection_icepack6

      !use mod_mesh
      use mice_module
      use icepack_intfc,        only: icepack_aggregate
      use icepack_itd,          only: cleanup_itd
      use schism_glbl,          only: dt,nstep_ice
      !use g_config,             only: dt
      use icepack_intfc,         only: icepack_init_trcr
      implicit none
          
      integer (kind=int_kind) :: ntrcr, ntrace, narr, nbtrcr, i,     &
                                 nt,   nt1,    k,   n ,ktherm,narrays,it,ierr
      integer (kind=int_kind) :: nt_Tsfc, nt_qice, nt_qsno,                          & 
                                 nt_sice, nt_fbri, nt_iage, nt_FY, nt_alvl, nt_vlvl, &
                                 nt_apnd, nt_hpnd, nt_ipnd, nt_bgc_Nit, nt_bgc_S
      logical (kind=log_kind) :: tr_pond_topo, tr_pond_lvl,                          &
                                 tr_pond,      tr_aero,     tr_FY,                   &
                                 tr_iage,      flag_old,                             &
                                 flag_old_0
      real    (kind=dbl_kind) :: puny ,Tsfc ,rhos, Lfresh, jj, kk, njj, pi, trcmin, trcmax,totalicea,totaliceh, &
      total_a,total_h,momentum_u,momentum_v, total_e, totalenthi, hilyr, ridge_threshold, total_ai, totaliceai
    
      ! Tracer dependencies and additional arrays
       real (kind=dbl_kind), dimension(ncat) :: &
             tmp, exc
      integer (kind=int_kind), dimension(:),    allocatable    ::    &
              tracer_type    , & ! = 1, 2, or 3 (depends on 0, 1 or 2 other tracers)
              depend             ! tracer dependencies (see below)
    
      logical (kind=log_kind), dimension (:),   allocatable   ::     &
              has_dependents    ! true if a tracer has dependent tracers
    
      real (kind=dbl_kind),    dimension (:,:), allocatable ::       &
             works,works_old, vicen_init, aicen_init, vsnon_init, aicen_tmp, iniflag, hicen_init, &
               &hsnon_init
      real (kind=dbl_kind),    dimension (:), allocatable ::       &
             fiso_ocn,swild,swild1,swild2, adv_threshold0
       real (kind=dbl_kind),    dimension (:,:,:), allocatable ::       &
       trcrn_ini

       real (kind=dbl_kind), dimension(nilyr) :: &
         qin            , & ! ice enthalpy (J/m3)
         qin_max        , & ! maximum ice enthalpy (J/m3)
         zTin               ! initial ice temperature

      real (kind=dbl_kind), dimension(nslyr) :: &
         qsn            , & ! snow enthalpy (J/m3)
         zTsn               ! initial snow temperature

      !type(t_mesh),        target,       intent(in)    :: mesh
         character(len=*), parameter :: subname = '(tracer_advection_icepack6)'

      call icepack_query_parameters(ktherm_out=ktherm,    &
                                    puny_out=puny,rhos_out=rhos,                   &
                                    Lfresh_out=Lfresh  ,pi_out=pi)
      call icepack_query_tracer_sizes(ntrcr_out=ntrcr, nbtrcr_out=nbtrcr)
      call icepack_query_tracer_flags(                                  &
             tr_iage_out=tr_iage, tr_FY_out=tr_FY,                      &
             tr_aero_out=tr_aero, tr_pond_out=tr_pond,                  &
             tr_pond_lvl_out=tr_pond_lvl, tr_pond_topo_out=tr_pond_topo)
      call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_qice_out=nt_qice, &
             nt_qsno_out=nt_qsno, nt_sice_out=nt_sice, nt_fbri_out=nt_fbri,       &
             nt_iage_out=nt_iage, nt_FY_out=nt_FY, nt_alvl_out=nt_alvl,           &
             nt_vlvl_out=nt_vlvl, nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd,       &
             nt_ipnd_out=nt_ipnd, nt_bgc_Nit_out=nt_bgc_Nit, nt_bgc_S_out=nt_bgc_S)
    
      narr   = 1 + ncat * (3 + ntrcr) ! max number of state variable arrays
    
      ! Allocate works array
      flag_old = .false.
      flag_old_0 = .false.

      if (allocated(works)) deallocate(works)
      if (allocated(works_old)) deallocate(works_old)
      allocate ( works(nx,narr) ) !(npa, max number of state variable arrays)
      allocate ( works_old(nx,narr) )
      allocate ( fiso_ocn(nx) )
      allocate ( swild(nx) )
      allocate ( swild1(nx) )
      allocate ( swild2(nx) )
      allocate ( adv_threshold0(nx) )
      allocate ( iniflag(nx,ncat) )
      allocate (trcrn_ini(nx,ntrcr,ncat) )
      allocate (vicen_init(nx,ncat) )
      allocate (aicen_init(nx,ncat) )
      allocate (hicen_init(nx,ncat) )
      allocate (hsnon_init(nx,ncat) )
      allocate (aicen_tmp(nx,ncat) )
      allocate (vsnon_init(nx,ncat) )
      works(:,:) = c0
      works_old(:,:) = c0
      fiso_ocn(:) = c0
      swild(:) = c0
      swild1(:) = c0
      swild2(:) = c0

      totalicea = 0
      totaliceh = 0
      momentum_u = 0
      momentum_v = 0
      totalenthi = 0
      totaliceai = 0
      do i=1,nx_nh
         if(associated(ipgl(iplg(i))%next)) then !interface node
            if(ipgl(iplg(i))%next%rank<myrank) cycle !already in the sum so skip
          endif
         totalicea = totalicea + sum(aicen(i,:))*lump_ice_matrix(i)
         totaliceai = max(totaliceai,sum(aicen(i,:))+aice0(i))
         !totaliceh = totaliceh + aice(i)*vice(i)*lump_ice_matrix(i)
         totaliceh = totaliceh + sum(vicen(i,:))*lump_ice_matrix(i)
         momentum_u = momentum_u + sum(vicen(i,:))*uvel(i)*lump_ice_matrix(i)
         momentum_v = momentum_v + sum(vicen(i,:))*vvel(i)*lump_ice_matrix(i)
         do n =1,ncat
            if(aicen(i,n)>puny) then
               hilyr = vicen(i,n) / aicen(i,n) /nilyr
            else
               hilyr = 0
            endif
            do k = 1,nilyr
               totalenthi = totalenthi + aicen(i,n)*hilyr*trcrn(i,nt_qice+k-1,n)*lump_ice_matrix(i)
            enddo
         enddo 

      enddo
      call mpi_allreduce(totalicea,total_a,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(totaliceh,total_h,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(totalenthi,total_e,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(totaliceai,total_ai,1,rtype,MPI_MAX,comm,ierr)
      if(myrank == 0) write(16,*) 'before advec. ice ',total_a,total_h,total_ai
      call mpi_allreduce(momentum_u,total_a,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(momentum_v,total_h,1,rtype,MPI_SUM,comm,ierr)
      if(myrank == 0) write(16,*) 'before advec. ice momentum',total_a,total_h

      trcrn_ini=trcrn
      vicen_init=vicen
      aicen_init=aicen
      vsnon_init=vsnon
      hicen_init(:,:) = c0  
      hsnon_init(:,:) = c0

      call state_to_work (ntrcr, narr, works(:,:))
       works_old=works
      ! Advect each tracer
    
      !do nt = 1, narr    
      !      call upwind_icepack (works(:,nt) )
      !end do
!upwind2
      !do i=1,nx
      !   adv_threshold0(i)= 0.01 !max(c0,min(0.05,exp(-7*aice0(i))))
      !enddo
      do i=1,nx
         if(vice(i)>puny) then
            adv_threshold0(i)= 10 !max(c0,min(0.01,exp(-10*aice(i))/vice(i)))
         else
            adv_threshold0(i)= 10
         endif
      enddo
      
      dudicex(:,:)=c0
      dudicey(:,:)=c0
      call tvd_casulli_ini

      call tvd_casulli_icepack (works(:,1),1,adv_threshold0 )
      
      narrays = 1
       do nt = 1, ncat
          call tvd_casulli_icepack (works(:,narrays+1),narrays+1,adv_threshold0) 
          aicen_tmp(:,nt) = works(:,narrays+1)
          narrays = narrays + 3 + ntrcr
       enddo

      swild = works(:,1)
       do nt = 1, ncat
         do i=1,nx
          swild(i) = swild(i) + aicen_tmp(i,nt)
         enddo
       enddo

       ridge_threshold = 1/864 + 1
       !do i=1,nx
       !   if(swild(i)>1.005) write(12,*) i,swild(i), works(i,1),aicen_tmp(i,1),aicen_tmp(i,2),aicen_tmp(i,3),&
       !   &aicen_tmp(i,4),aicen_tmp(i,5),aice0(i),aicen(i,1),aicen(i,2),aicen(i,3),aicen(i,4),aicen(i,5)
          !if(swild(i)>ridge_threshold) then
          !  works(i,1) = works(i,1) / swild(i) * ridge_threshold
            !aice0(i) = aice0(i) / swild(i) * ridge_threshold
          !  narrays = 1
          !  do nt = 1, ncat
          !     works(i,narrays+1) = works(i,narrays+1) / swild(i) * ridge_threshold
          !     !aicen(i,nt) = aicen(i,nt) / swild(i) * ridge_threshold
          !     aicen_tmp(i,nt) = works(i,narrays+1)
          !     flux_matrix3(:,i,narrays+1) = flux_matrix3(:,i,narrays+1) / swild(i) * ridge_threshold
          !     flux_matrix4(:,i,narrays+1) = flux_matrix4(:,i,narrays+1) / swild(i) * ridge_threshold
          !     flux_matrix1(:,i,narrays+1) = flux_matrix1(:,i,narrays+1) / swild(i) * ridge_threshold
          !     flux_matrix2(:,i,narrays+1) = flux_matrix2(:,i,narrays+1) / swild(i) * ridge_threshold
          !     narrays = narrays + 3 + ntrcr
          !  enddo
          !endif
       !enddo



         !narrays = 1
         !do nt = 1, ncat
         !   totalicea = 0
         !   do i=1,ne
         !      totalicea = totalicea + sum((flux_matrix1(:,elnode(1:i34(i),i),narrays+1)+flux_matrix2(:,elnode(1:i34(i),i),narrays+1))/2+ &
         !      &(flux_matrix3(:,elnode(1:i34(i),i),narrays+1)+flux_matrix4(:,elnode(1:i34(i),i),narrays+1))/2)
         !   enddo
         !   narrays = narrays + 3 + ntrcr
         !   call mpi_allreduce(totalicea,total_a,1,rtype,MPI_SUM,comm,ierr)
         !   if(myrank == 0) write(16,*) 'flux_matrix advec. ice ',nt,total_a
         !enddo
       narrays = 1
       do nt = 1, ncat
          !call upwind_icepack (works(:,narrays+2),narrays+2 )  
          !call upwind_icepack (works(:,narrays+3),narrays+3 )
          do i = 1, nx
               if(aicen_init(i,nt)<puny.or.vsnon_init(i,nt)<puny) then
                  swild2(i) = 0
                  hsnon_init(i,nt) = 0
               else
                  swild2(i) = vsnon_init(i,nt)/aicen_init(i,nt)
                  hsnon_init(i,nt) = vsnon_init(i,nt)/aicen_init(i,nt)
               endif

               if(aicen_init(i,nt)<puny.or.vicen_init(i,nt)<puny) then
                  swild1(i) = 0
                  swild2(i) = 0
                  hicen_init(i,nt) = 0
                  hsnon_init(i,nt) = 0
               else
                  swild1(i) = vicen_init(i,nt)/aicen_init(i,nt)
                  hicen_init(i,nt) = vicen_init(i,nt)/aicen_init(i,nt)
               endif
          enddo
          call tvd_casulli_icepack_dep(works(:,narrays+2),swild1,narrays+2,narrays+1)
          call tvd_casulli_icepack_dep(works(:,narrays+3),swild2,narrays+3,narrays+1) 
          !call tvd_icepack (works(:,narrays+2),narrays+2,adv_threshold0 ) 
          !call tvd_icepack (works(:,narrays+3),narrays+3,adv_threshold0 ) 
          narrays = narrays + 3 + ntrcr
       enddo

       

      call work_to_other_6  (ntrcr, narr, works(:,:)) 
      
      call re_momentum (ntrcr, narr, works(:,:), vicen_init)
      newwork(:) = c0


      call work_to_state (ntrcr, narr, works(:,:))
      iniflag(:,:) = c0
      swild1(:) = c0
      swild2(:) = c0
         do nt = 1, ncat
            do i=1,nx_nh
               do it = 1, ntrcr
                  trcmin = minval(trcrn_ini(indnd(1:nnp(i),i),it,nt))
                  trcmin = min(trcmin,trcrn_ini(i,it,nt))
                  trcmax = maxval(trcrn_ini(indnd(1:nnp(i),i),it,nt))
                  trcmax = max(trcmax,trcrn_ini(i,it,nt))
                  if(trcrn(i,it,nt) > trcmax) iniflag(i,nt) = 1 !trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmax
                  if(trcrn(i,it,nt) < trcmin) iniflag(i,nt) = 1 !trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmin
                  !if(trcrn(i,it,nt) > trcmax) trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmax
                  !if(trcrn(i,it,nt) < trcmin) trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmin
                  !if(trcrn(i,it,nt) > trcmax) trcrn(i,it,nt) = trcmax !trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmax
                  !if(trcrn(i,it,nt) < trcmin) trcrn(i,it,nt) = trcmin !trcrn(i,it,nt) = trcrn_ini(i,it,nt) !trcrn_ini(i,it,nt) !trcmin
               enddo

                  trcmax = maxval(hicen_init(indnd(1:nnp(i),i),nt))
                  trcmax = max(trcmax,hicen_init(i,nt))
                  if(aicen(i,nt)<puny.or.vicen(i,nt)<puny) then
                     swild1(i) = 0
                  else
                     swild1(i) = vicen(i,nt)/aicen(i,nt)
                  endif
                  ! reduce hicen by increase aicen
                  !if(swild1(i) > trcmax .and. aicen(i,nt)<0.15 ) then
                  if(swild1(i) > trcmax) then
                     !aicen(i,nt) = vicen(i,nt) / min(trcmax,hin_max(nt))
                     vicen(i,nt) = aicen(i,nt) * min(trcmax,hin_max(nt))
                  endif

                  trcmax = maxval(hsnon_init(indnd(1:nnp(i),i),nt))
                  trcmax = max(trcmax,hsnon_init(i,nt))
                  if(aicen(i,nt)<puny.or.vsnon(i,nt)<puny) then
                     swild2(i) = 0
                  else
                     swild2(i) = vsnon(i,nt)/aicen(i,nt)
                  endif
                  ! reduce hsno
                  !if(swild2(i) > trcmax .and. aicen(i,nt)<0.15 ) then
                  if(swild2(i) > trcmax) then
                     vsnon(i,nt) = trcmax * aicen(i,nt)
                     !vsnon(i,nt) = hsnon_init(i,nt) * aicen(i,nt)
                  endif


            enddo
         enddo
         do nt = 1, ncat
            swild = iniflag(:,nt)
            call exchange_p2d(swild)
            iniflag(:,nt) = swild

            swild = aicen(:,nt)
            call exchange_p2d(swild)
            aicen(:,nt) = swild

            swild = vicen(:,nt)
            call exchange_p2d(swild)
            vicen(:,nt) = swild

            swild = vsnon(:,nt)
            call exchange_p2d(swild)
            vsnon(:,nt) = swild
         enddo

       do nt = 1, ncat
             do i=1,nx
                if(iniflag(i,nt) > 0) then
                  trcrn(i,:,nt) = trcrn_ini(i,:,nt)
      !             trcrn(i,:,:) = trcrn_ini(i,:,:)
                endif
             enddo
        enddo

      swild = c0
       
       do i=1,nx
         do nt = 1, ncat
          swild(i) = swild(i) + aicen(i,nt)
         enddo
          aice0(i) = max(c0,c1-swild(i))
       enddo


      !do i=1,nx
      !   if(isbnd(1,i)/=0) then !b.c. (including open)
      !      do nt = 1, ncat
      !         aicen(i,nt) = aicen_init(i,nt)
      !         vicen(i,nt) = vicen_init(i,nt)
      !         vsnon(i,nt) = vsnon_init(i,nt)
      !         trcrn(i,:,nt) = trcrn_ini(i,:,nt)
      !      enddo
      !   endif
      !enddo
      ! cut off icepack
    
      !call cut_off_icepack

       !write(12,*) maxval(swild(:))


      do i=1,nx
         if (ncat < 0) then ! Do we really need this?
            call cleanup_itd  (dt*nstep_ice,           ntrcr,                &
                               nilyr,                  nslyr,                &
                               ncat,                   hin_max(:),           &
                               aicen(i,:),             trcrn(i,1:ntrcr,:),   &
                               vicen(i,:),             vsnon(i,:),           &
                               aice0(i),               aice(i),              &
                               n_aero,                                       &
                               nbtrcr,                 nblyr,                &
                               tr_aero,                                      &
                               tr_pond_topo,                                 &
                               first_ice(i,:),                               &
                               trcr_depend(1:ntrcr),   trcr_base(1:ntrcr,:), &
                               n_trcr_strata(1:ntrcr), nt_strata(1:ntrcr,:), &
                               fpond(i),               fresh(i),             &
                               fsalt(i),               fhocn(i),             &
                               faero_ocn(i,:),         fiso_ocn,          &
                               flux_bio(i,1:nbtrcr))
    
            call icepack_aggregate (ncat,                    &
                                   aicen(i,:),               &
                                   trcrn(i,1:ntrcr,:),       &
                                   vicen(i,:),               &
                                   vsnon(i,:),               &
                                   aice (i),                 &
                                   trcr (i,1:ntrcr),         &
                                   vice (i),                 &
                                   vsno (i),                 &
                                   aice0(i),                 &
                                   ntrcr,                    &
                                   trcr_depend  (1:ntrcr),   &
                                   trcr_base    (1:ntrcr,:), &
                                   n_trcr_strata(1:ntrcr),   &
                                   nt_strata    (1:ntrcr,:))
         endif
      end do
      totalicea = 0
      totaliceh = 0
      momentum_u = 0
      momentum_v = 0
      totalenthi = 0
      totaliceai = 0
      do i=1,nx_nh
         if(associated(ipgl(iplg(i))%next)) then !interface node
            if(ipgl(iplg(i))%next%rank<myrank) cycle !already in the sum so skip
          endif
       totalicea = totalicea + sum(aicen(i,:))*lump_ice_matrix(i)
       totaliceai = max(totaliceai,sum(aicen(i,:))+aice0(i))
       !totaliceh = totaliceh + aice(i)*vice(i)*lump_ice_matrix(i)
       totaliceh = totaliceh + sum(vicen(i,:))*lump_ice_matrix(i)
       momentum_u = momentum_u + sum(vicen(i,:))*uvel(i)*lump_ice_matrix(i)
       momentum_v = momentum_v + sum(vicen(i,:))*vvel(i)*lump_ice_matrix(i)
        do n =1,ncat
            if(aicen(i,n)>puny) then
               hilyr = vicen(i,n) / aicen(i,n) /nilyr
               !totaliceai = max(totaliceai,vicen(i,n) / aicen(i,n))
            else
               hilyr = 0
            endif
            do k = 1,nilyr
               totalenthi = totalenthi + aicen(i,n)*hilyr*trcrn(i,nt_qice+k-1,n)*lump_ice_matrix(i)
            enddo
         enddo 
      enddo
      call mpi_allreduce(totalicea,total_a,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(totaliceh,total_h,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(totalenthi,total_e,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(totaliceai,total_ai,1,rtype,MPI_MAX,comm,ierr)
      if(myrank == 0) write(16,*) 'after advec. ice ',total_a,total_h,total_e,total_ai
      call mpi_allreduce(momentum_u,total_a,1,rtype,MPI_SUM,comm,ierr)
      call mpi_allreduce(momentum_v,total_h,1,rtype,MPI_SUM,comm,ierr)
      if(myrank == 0) write(16,*) 'after advec. ice momentum',total_a,total_h
      deallocate(works)
      deallocate(works_old)
      deallocate(trcrn_ini)
      deallocate(vicen_init)
      deallocate(aicen_init)
      deallocate(vsnon_init)
      deallocate(hicen_init)
      deallocate(hsnon_init)
      deallocate(swild)
      deallocate(swild1)
      deallocate(swild2)
      deallocate(adv_threshold0)
  end subroutine tracer_advection_icepack6

  !=======================================================================

    subroutine work_to_state (ntrcr, narr, works)

        use icepack_intfc,    only: icepack_compute_tracers,icepack_aggregate

        integer (kind=int_kind), intent(in) :: ntrcr, narr

        real (kind=dbl_kind), dimension(nx,narr), intent (inout) :: &
           works     ! work array
  
        ! local variables
  
        integer (kind=int_kind) :: &
           nt_alvl, nt_apnd, nt_fbri, nt_Tsfc, ktherm ,nt_vlvl
  
        logical (kind=log_kind) :: &
           tr_lvl,  tr_pond_lvl, tr_pond_topo, tr_brine
  
        integer (kind=int_kind) ::      &
           k, i, n, it   , & ! counting indices
           narrays       , & ! counter for number of state variable arrays
           nt_qsno       , &
           nt_qice       , &
           nt_sice
  
        real (kind=dbl_kind) :: &
           rhos       , &
           rhoi       , &
           Lfresh     , &
           Tsmelt
  
        real (kind=dbl_kind), dimension(ncat) :: &
           tmp, exc

        real (kind=dbl_kind) :: puny,small
        real (kind=dbl_kind),    dimension (:,:), allocatable ::       &
               aicen_tmp

        !real (kind=dbl_kind), parameter :: &
        !   small = 0.000001_dbl_kind
         
        character(len=*), parameter :: subname = '(state_to_work)'
  
        call icepack_query_tracer_flags(tr_pond_lvl_out=tr_pond_lvl, tr_pond_topo_out=tr_pond_topo,      &
             tr_lvl_out=tr_lvl,tr_brine_out=tr_brine)
        call icepack_query_tracer_indices(nt_alvl_out=nt_alvl, nt_apnd_out=nt_apnd, &
             nt_fbri_out=nt_fbri, nt_qsno_out=nt_qsno, nt_vlvl_out=nt_vlvl,                              &
             nt_qice_out=nt_qice, nt_sice_out=nt_sice, nt_Tsfc_out=nt_Tsfc) 
        call icepack_query_parameters(rhoi_out=rhoi,     rhos_out=rhos,                   &
                                      Lfresh_out=Lfresh, &
                                      Tsmelt_out=Tsmelt, ktherm_out=ktherm,               &
                                      puny_out=puny)
        call icepack_warnings_flush(ice_stderr)
        if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
           file=__FILE__, line=__LINE__)

        allocate (aicen_tmp(nx,ncat) )
        aicen_tmp(:,:) = c0
        trcrn(:,:,:) = c0
        aicen(:,:)   = c0
        vicen(:,:)   = c0
        vsnon(:,:)   = c0
        aice0(:)     = c0

        small       = puny
        ! Open water fraction
  
        do i = 1, nx
           !if (works(i,1) <= puny) then
           !    aice0(i) = c0
           !else if (works(i,1) >= c1) then
           !    aice0(i) = c1
           if (works(i,1) <= small) then
               aice0(i) = c0
           else
               aice0(i) = works(i,1)
           end if
        enddo
        narrays = 1
  
        ! Sea ice area and volume per unit area of ice and snow
  
        do n=1,ncat
           do i = 1, nx
              !if (works(i,narrays+1) > c1) then
              !    works(i,narrays+1) = c1
              !end if
              if (works(i,narrays+1) <= small .or. works(i,narrays+2) <= small) then
                  works(i,narrays+1) = c0
                  works(i,narrays+2) = c0
                  works(i,narrays+3) = c0
              end if
              if (works(i,narrays+3) <= small) then
                 works(i,narrays+3) = c0
              end if
              aicen(i,n) = works(i,narrays+1)
              vicen(i,n) = works(i,narrays+2)
              vsnon(i,n) = works(i,narrays+3)
           end do
           narrays = narrays + 3 + ntrcr
        end do
  
 
        narrays = 1

        do n=1, ncat
            narrays = narrays + 3

           do i = 1, nx

                  call icepack_compute_tracers(ntrcr=ntrcr,trcr_depend=trcr_depend(:),    &
                                                atrcrn = works(i,narrays+1:narrays+ntrcr), &
                                                aicen  = aicen(i,n),                       &
                                                vicen  = vicen(i,n),                       &
                                                vsnon  = vsnon(i,n),                       &
                                                trcr_base     = trcr_base(:,:),            &
                                                n_trcr_strata = n_trcr_strata(:),          &
                                                nt_strata     = nt_strata(:,:),            &
                                                trcrn  = trcrn(i,:,n))
                  trcrn(i,nt_qsno:(nt_qsno+nslyr-1),n) = trcrn(i,nt_qsno:(nt_qsno+nslyr-1),n) - rhos*Lfresh
           enddo
          
           narrays = narrays + ntrcr
        enddo
       !  aicen_tmp = aicen
       ! do i = 1, nx  ! For each grid cell
       !    if (sum(aicen(i,:)) > c1) then
       !     !write(12,*) 'sum(aicen(i,:)) > c1',sum(aicen(i,:)),aicen(i,:)
       !       tmp(:) = c0
       !       exc(:) = c0
       !       do n = 1, ncat
       !          if (aicen(i,n) > puny) tmp(n) = c1
       !       end do
       !       do n = 1, ncat
       !           exc(n) = max(c0,(sum(aicen(i,:)) - c1))  &
       !                    * aicen(i,n) / sum(aicen(i,:))
       !       end do
       !       do n = 1, ncat
       !          aicen(i,n) = max(c0,aicen(i,n) - exc(n))
       !          aice0(i)   = max(c0,1-sum(aicen(i,:)))
       !       end do
       !    end if
       ! end do

        ! do n = 1, ncat
        !    do i = 1, nx
        !       if(aicen_tmp(i,n) > puny) then
        !          vicen(i,n) = vicen(i,n) / aicen_tmp(i,n) * aicen(i,n)
        !          vsnon(i,n) = vsnon(i,n) / aicen_tmp(i,n) * aicen(i,n)
        !       else
        !          vicen(i,n) = c0
        !          vsnon(i,n) = c0
        !       endif
        !    enddo
        ! enddo

        !do i = 1, nx
        ! aice(i) = c0
        ! vice(i) = c0
        ! vsno(i) = c0
        ! do it = 1, ntrcr
        !    trcr(i,it) = c0
        ! enddo
        ! call icepack_aggregate (ncat,                    &
        !                        aicen(i,:),               &
        !                        trcrn(i,1:ntrcr,:),       &
        !                        vicen(i,:),               &
        !                        vsnon(i,:),               &
        !                        aice (i),                 &
        !                        trcr (i,1:ntrcr),         &
        !                        vice (i),                 &
        !                        vsno (i),                 &
        !                        aice0(i),                 &
        !                        ntrcr,                    &
        !                        trcr_depend  (1:ntrcr),   &
        !                        trcr_base    (1:ntrcr,:), &
        !                        n_trcr_strata(1:ntrcr),   &
        !                        nt_strata    (1:ntrcr,:))
        ! end do


!        do n=1, ncat
!  
!           narrays = narrays + 3
!  
!           do it = 1, ntrcr
!  
!              if (trcr_depend(it) == 0) then
!                 do i = 1, nx
!                    if (aicen(i,n) > c0) then
!                       if (it == nt_Tsfc) then
!                           trcrn(i,it,n) = min(c0,works(i,narrays+it)/aicen(i,n))
!                       else if (it == nt_alvl .or. it == nt_apnd) then
!                           trcrn(i,it,n) = max(c0,min(c1,works(i,narrays+it) / aicen(i,n)))
!                       endif
!                    end if
!                 enddo
!              elseif (trcr_depend(it) == 1) then
!                 do i = 1, nx
!                    if (vicen(i,n) > c0) then
!                        if (it >= nt_qice .and. it < nt_qice+nilyr) then
!                            trcrn(i,it,n) = min(c0,works(i,narrays+it)/vicen(i,n))
!                            if (.not. heat_capacity) trcrn(i,it,n) = -rhoi * Lfresh
!                        else if (it >= nt_sice .and. it < nt_sice+nilyr) then
!                            trcrn(i,it,n) = max(c0,works(i,narrays+it)/vicen(i,n))
!                        else
!                            trcrn(i,it,n) = max(c0,works(i,narrays+it)/vicen(i,n))
!                        end if
!                    end if
!                 enddo
!              elseif (trcr_depend(it) == 2) then
!                 do i = 1, nx
!                     if (vsnon(i,n) > c0) then
!                         if (it >= nt_qsno .and. it < nt_qsno+nslyr) then
!                             trcrn(i,it,n) = min(c0,works(i,narrays+it)/vsnon(i,n)) - rhos*Lfresh
!                             if (.not. heat_capacity) trcrn(i,it,n) = -rhos * Lfresh
!                         else
!                             trcrn(i,it,n) = min(c0,works(i,narrays+it)/vsnon(i,n)) - rhos*Lfresh
!                         end if
!                     end if
!                 enddo
!              elseif (trcr_depend(it) == 2+nt_alvl) then
!                 do i = 1, nx
!                    if (trcrn(i,nt_alvl,n) > small) then
!                       trcrn(i,it,n) = max( c0, works(i,narrays+it)  &
!                                           / trcrn(i,nt_alvl,n) )     
!                    endif
!                 enddo
!              elseif (trcr_depend(it) == 2+nt_apnd .and. &
!                      tr_pond_cesm .or. tr_pond_topo) then
!                 do i = 1, nx
!                    if (trcrn(i,nt_apnd,n) > small) then
!                       trcrn(i,it,n) = max( c0, works(i,narrays+it)  &
!                                           / trcrn(i,nt_apnd,n) )      
!                    endif
!                 enddo
!              elseif (trcr_depend(it) == 2+nt_apnd .and. &
!                      tr_pond_lvl) then
!                 do i = 1, nx
!                    if (trcrn(i,nt_apnd,n) > small .and. trcrn(i,nt_alvl,n) > small) then
!                       trcrn(i,it,n) = max( c0, works(i,narrays+it)  &
!                                           / trcrn(i,nt_apnd,n)      &
!                                           / trcrn(i,nt_alvl,n)      &
!                                           / aicen(i,n) )
!                    endif
!                 enddo
!              elseif (trcr_depend(it) == 2+nt_fbri) then
!                 do i = 1, nx
!                        works(i,narrays+it) = vicen(i,n) &
!                                            / trcrn(i,nt_fbri,n) &
!                                            / trcrn(i,it,n)
!                 enddo
!              endif
!           enddo
!  
!           narrays = narrays + ntrcr
!  
!        enddo                 ! number of categories 

!        if (mype == 0) write(*,*) 'Tracer salinity: ', nt_sice, ' - ', (nt_sice + nilyr - 1)
!        if (mype == 0) write(*,*) 'ktherm: ', ktherm

        if (ktherm == 1) then ! For bl99 themodynamics
                              ! always ridefine salinity
                              ! after advection
           do i = 1, nx       ! For each grid cell
              do k = 1, nilyr
                 trcrn(i,nt_sice+k-1,:) = salinz(i,k)
              end do          ! nilyr
           end do
        end if                ! ktherm==1

    end subroutine work_to_state

    !=======================================================================

    subroutine state_to_work (ntrcr, narr, works)

        integer (kind=int_kind), intent(in) :: ntrcr, narr
  
        real (kind=dbl_kind), dimension(nx,narr), intent (out) :: &
           works     ! work array
  
        ! local variables
  
        integer (kind=int_kind) :: &
           nt_alvl, nt_apnd, nt_fbri, nt_Tsfc, nt_qice
  
        logical (kind=log_kind) :: &
           tr_pond_lvl, tr_pond_topo
  
        integer (kind=int_kind) ::      &
           i, n, it   , & ! counting indices
           narrays    , & ! counter for number of state variable arrays
           nt_qsno,k,itl,ntr
  
        real (kind=dbl_kind) :: &
           rhos       , &
           Lfresh,small
  
        character(len=*), parameter :: subname = '(state_to_work)'
  
        call icepack_query_tracer_flags(tr_pond_lvl_out=tr_pond_lvl, tr_pond_topo_out=tr_pond_topo)
        call icepack_query_tracer_indices(nt_alvl_out=nt_alvl, nt_apnd_out=nt_apnd, &
             nt_fbri_out=nt_fbri, nt_qsno_out=nt_qsno, nt_Tsfc_out=nt_Tsfc,nt_qice_out=nt_qice)
        call icepack_query_parameters(rhos_out=rhos, Lfresh_out=Lfresh)
        call icepack_warnings_flush(ice_stderr)
        if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
           file=__FILE__, line=__LINE__)

         works(:,:) = c0
         small = 0.0000001
        do i = 1, nx
           works(i,1) = aice0(i)
        enddo
        narrays = 1
  
        do n=1, ncat
  
           do i = 1, nx
              works(i,narrays+1) = aicen(i,n)
              works(i,narrays+2) = vicen(i,n)
              works(i,narrays+3) = vsnon(i,n)

              !if (vicen(i,n)>10) write(12,*) 'before advec.', i, n, vicen(i,n)
           enddo                  ! i
           narrays = narrays + 3
  
           do it = 1, ntrcr
            !do i = 1, nx
            !      if (it >= nt_qsno .and. it < nt_qsno+nslyr) then
            !         works(i,narrays+it)= (trcrn(i,it,n)+ rhos*Lfresh)*(trcr_base(it,1) * aicen(i,n) &
            !                                       + trcr_base(it,2) * vicen(i,n) &
            !                                       + trcr_base(it,3) * vsnon(i,n))
            !      else
            !      works(i,narrays+it)= trcrn(i,it,n)*(trcr_base(it,1) * aicen(i,n) &
            !                                       + trcr_base(it,2) * vicen(i,n) &
            !                                       + trcr_base(it,3) * vsnon(i,n))
            !      endif
            !      if (n_trcr_strata(it) > 0) then    ! additional tracer layers
            !         do itl = 1, n_trcr_strata(it)
            !            ntr = nt_strata(it,itl)
            !            works(i,narrays+it) = works(i,narrays+it) * trcrn(i,ntr,n)
            !         enddo
            !      endif
            !   enddo

              if (trcr_depend(it) == 0) then
                 do i = 1, nx
                    works(i,narrays+it) = aicen(i,n)*trcrn(i,it,n)
                 enddo
              elseif (trcr_depend(it) == 1) then
                 do i = 1, nx
                    works(i,narrays+it) = vicen(i,n)*trcrn(i,it,n)
                 enddo
              elseif (trcr_depend(it) == 2) then
                 do i = 1, nx
                    if (it >= nt_qsno .and. it < nt_qsno+nslyr) then
                        !works(i,narrays+it) = vsnon(i,n)*trcrn(i,it,n)
                        works(i,narrays+it) = vsnon(i,n)*(trcrn(i,it,n)+ rhos*Lfresh) 
                    else
                        works(i,narrays+it) = vsnon(i,n)*trcrn(i,it,n)
                    end if
                 enddo
              elseif (trcr_depend(it) == 2+nt_alvl) then
                 do i = 1, nx
                    works(i,narrays+it) = aicen(i,n) &
                                        * trcrn(i,nt_alvl,n) &
                                        * trcrn(i,it,n)
                 enddo
              elseif (trcr_depend(it) == 2+nt_apnd .and. &
                      tr_pond_topo) then
                 do i = 1, nx
                    works(i,narrays+it) = aicen(i,n) &
                                        * trcrn(i,nt_apnd,n) &
                                        * trcrn(i,it,n)
                 enddo
              elseif (trcr_depend(it) == 2+nt_apnd .and. &
                      tr_pond_lvl) then
                 do i = 1, nx
                    works(i,narrays+it) = aicen(i,n) &
                                        * trcrn(i,nt_alvl,n) &
                                        * trcrn(i,nt_apnd,n) &
                                        * trcrn(i,it,n)
                 enddo
              elseif (trcr_depend(it) == 2+nt_fbri) then
                 do i = 1, nx
                    works(i,narrays+it) = vicen(i,n) &
                                        * trcrn(i,nt_fbri,n) &
                                        * trcrn(i,it,n)
                 enddo
              endif
           enddo
           narrays = narrays + ntrcr
  
        enddo                     ! n

        
        if (narr /= narrays .and. myrank == 0 ) write(ice_stderr,*)      &
            "Wrong number of arrays in transport bound call"

    end subroutine state_to_work

!=======================================================================

    subroutine work_to_other (ntrcr, narr, works)

        integer (kind=int_kind), intent(in) :: ntrcr, narr
  
        real (kind=dbl_kind), dimension(nx,narr), intent (inout) :: &
           works     ! work array
  
        ! local variables
  
        integer (kind=int_kind) :: &
           nt_alvl, nt_apnd, nt_fbri, nt_Tsfc, nt_qice, &
           &arr_alvl,arr_apnd,arr_fbri
  
        logical (kind=log_kind) :: &
           tr_pond_lvl, tr_pond_topo
  
        integer (kind=int_kind) ::      &
           i, n, it   , & ! counting indices
           narrays    , & ! counter for number of state variable arrays
         nt_qsno,k,j
  
        real (kind=dbl_kind) :: &
           rhos       , &
           Lfresh,small
        real (kind=dbl_kind),    dimension (:), allocatable ::       &
           swild
        character(len=*), parameter :: subname = '(work_to_other)'
  
        allocate ( swild(nx) )
        call icepack_query_tracer_flags(tr_pond_lvl_out=tr_pond_lvl, tr_pond_topo_out=tr_pond_topo)
        call icepack_query_tracer_indices(nt_alvl_out=nt_alvl, nt_apnd_out=nt_apnd, &
             nt_fbri_out=nt_fbri, nt_qsno_out=nt_qsno, nt_Tsfc_out=nt_Tsfc,nt_qice_out=nt_qice)
        call icepack_query_parameters(rhos_out=rhos, Lfresh_out=Lfresh)
        call icepack_warnings_flush(ice_stderr)
        if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
           file=__FILE__, line=__LINE__)

        small = 0.0000001
        swild(:) = c0
        narrays = 1
  
        do n=1, ncat
    
           do it = 1, ntrcr
              if (trcr_depend(it) == 0) then
                 call upwind_icepack_other(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+1)
              elseif (trcr_depend(it) == 1) then
                 call upwind_icepack_other(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+2)
              elseif (trcr_depend(it) == 2) then
               if (it >= nt_qsno .and. it < nt_qsno+nslyr) then
                     swild = trcrn(:,it,n) + rhos*Lfresh
                     call upwind_icepack_other(works(:,narrays+3+it),swild,narrays+3+it,narrays+3)
               else
                     call upwind_icepack_other(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+3)
               endif
              elseif (trcr_depend(it) == 2+nt_alvl) then
                call upwind_icepack_other(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+1,trcrn(:,nt_alvl,n))!narrays+3+nt_alvl)
              elseif (trcr_depend(it) == 2+nt_apnd .and. &
                      tr_pond_topo) then
               call upwind_icepack_other(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+1,trcrn(:,nt_apnd,n))!narrays+3+nt_apnd)
              elseif (trcr_depend(it) == 2+nt_apnd .and. &
                      tr_pond_lvl) then
                 call upwind_icepack_other(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+1,trcrn(:,nt_alvl,n),trcrn(:,nt_apnd,n))
                 !,narrays+3+nt_alvl,narrays+3+nt_apnd)
              elseif (trcr_depend(it) == 2+nt_fbri) then
                  call upwind_icepack_other(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+2,trcrn(:,nt_fbri,n))!,narrays+3+nt_fbri)
              endif
           enddo
           narrays = narrays + ntrcr +3
  
        enddo                     ! n

        
        if (narr /= narrays .and. myrank == 0 ) write(ice_stderr,*)      &
            "Wrong number of arrays in transport bound call"

    end subroutine work_to_other

!=======================================================================

    subroutine work_to_other2 (ntrcr, narr, works)

        integer (kind=int_kind), intent(in) :: ntrcr, narr
  
        real (kind=dbl_kind), dimension(nx,narr), intent (inout) :: &
           works     ! work array
  
        ! local variables
  
        integer (kind=int_kind) :: &
           nt_alvl, nt_apnd, nt_fbri, nt_Tsfc, nt_qice, &
           &arr_alvl,arr_apnd,arr_fbri
  
        logical (kind=log_kind) :: &
           tr_pond_lvl, tr_pond_topo
  
        integer (kind=int_kind) ::      &
           i, n, it   , & ! counting indices
           narrays    , & ! counter for number of state variable arrays
           nt_qsno,k,j
  
        real (kind=dbl_kind) :: &
           rhos       , &
           Lfresh,small
        real (kind=dbl_kind),    dimension (:), allocatable ::       &
           swild
        character(len=*), parameter :: subname = '(work_to_other2)'

        allocate ( swild(nx) )

        call icepack_query_tracer_flags(tr_pond_lvl_out=tr_pond_lvl, tr_pond_topo_out=tr_pond_topo)
        call icepack_query_tracer_indices(nt_alvl_out=nt_alvl, nt_apnd_out=nt_apnd, &
             nt_fbri_out=nt_fbri, nt_qsno_out=nt_qsno, nt_Tsfc_out=nt_Tsfc,nt_qice_out=nt_qice)
        call icepack_query_parameters(rhos_out=rhos, Lfresh_out=Lfresh)
        call icepack_warnings_flush(ice_stderr)
        if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
           file=__FILE__, line=__LINE__)

        small = 0.0000001
        swild(:) = c0
        narrays = 1
  
        do n=1, ncat
    
           do it = 1, ntrcr
              if (trcr_depend(it) == 0) then
                 call upwind_icepack_other2(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+1)
              elseif (trcr_depend(it) == 1) then
                 call upwind_icepack_other2(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+2)
              elseif (trcr_depend(it) == 2) then
                  if (it >= nt_qsno .and. it < nt_qsno+nslyr) then
                     swild = trcrn(:,it,n) + rhos*Lfresh
                     call upwind_icepack_other2(works(:,narrays+3+it),swild,narrays+3+it,narrays+3)
                  else
                     call upwind_icepack_other2(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+3)
                  endif 
              elseif (trcr_depend(it) == 2+nt_alvl) then
                call upwind_icepack_other2(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+1,trcrn(:,nt_alvl,n))!narrays+3+nt_alvl)
              elseif (trcr_depend(it) == 2+nt_apnd .and. &
                      tr_pond_topo) then
               call upwind_icepack_other2(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+1,trcrn(:,nt_apnd,n))!narrays+3+nt_apnd)
              elseif (trcr_depend(it) == 2+nt_apnd .and. &
                      tr_pond_lvl) then
                 call upwind_icepack_other2(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+1,trcrn(:,nt_alvl,n),trcrn(:,nt_apnd,n))
                 !,narrays+3+nt_alvl,narrays+3+nt_apnd)
              elseif (trcr_depend(it) == 2+nt_fbri) then
                  call upwind_icepack_other2(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+2,trcrn(:,nt_fbri,n))!,narrays+3+nt_fbri)
              endif
           enddo
           narrays = narrays + ntrcr +3
  
        enddo                     ! n

        
        if (narr /= narrays .and. myrank == 0 ) write(ice_stderr,*)      &
            "Wrong number of arrays in transport bound call"

    end subroutine work_to_other2

!=======================================================================

    subroutine work_to_other3 (ntrcr, narr, works)

        integer (kind=int_kind), intent(in) :: ntrcr, narr
  
        real (kind=dbl_kind), dimension(nx,narr), intent (inout) :: &
           works     ! work array
  
        ! local variables
  
        integer (kind=int_kind) :: &
           nt_alvl, nt_apnd, nt_fbri, nt_Tsfc, nt_qice, &
           &arr_alvl,arr_apnd,arr_fbri
  
        logical (kind=log_kind) :: &
           tr_pond_lvl, tr_pond_topo
  
        integer (kind=int_kind) ::      &
           i, n, it   , & ! counting indices
           narrays    , & ! counter for number of state variable arrays
         nt_qsno,k,j
  
        real (kind=dbl_kind) :: &
           rhos       , &
           Lfresh,small
        real (kind=dbl_kind),    dimension (:), allocatable ::       &
           swild
        character(len=*), parameter :: subname = '(work_to_other3)'
  
        allocate ( swild(nx) )
        call icepack_query_tracer_flags(tr_pond_lvl_out=tr_pond_lvl, tr_pond_topo_out=tr_pond_topo)
        call icepack_query_tracer_indices(nt_alvl_out=nt_alvl, nt_apnd_out=nt_apnd, &
             nt_fbri_out=nt_fbri, nt_qsno_out=nt_qsno, nt_Tsfc_out=nt_Tsfc,nt_qice_out=nt_qice)
        call icepack_query_parameters(rhos_out=rhos, Lfresh_out=Lfresh)
        call icepack_warnings_flush(ice_stderr)
        if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
           file=__FILE__, line=__LINE__)

        small = 0.0000001
        swild(:) = c0
        narrays = 1
  
        do n=1, ncat
    
           do it = 1, ntrcr
              if (trcr_depend(it) == 0) then
                 call tvd_icepack_other(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+1)
              elseif (trcr_depend(it) == 1) then
                 call tvd_icepack_other(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+2)
              elseif (trcr_depend(it) == 2) then
               if (it >= nt_qsno .and. it < nt_qsno+nslyr) then
                     swild = trcrn(:,it,n) + rhos*Lfresh
                     call tvd_icepack_other(works(:,narrays+3+it),swild,narrays+3+it,narrays+3)
               else
                     call tvd_icepack_other(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+3)
               endif
              elseif (trcr_depend(it) == 2+nt_alvl) then
                call tvd_icepack_other(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+1,trcrn(:,nt_alvl,n))!narrays+3+nt_alvl)
              elseif (trcr_depend(it) == 2+nt_apnd .and. &
                      tr_pond_topo) then
               call tvd_icepack_other(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+1,trcrn(:,nt_apnd,n))!narrays+3+nt_apnd)
              elseif (trcr_depend(it) == 2+nt_apnd .and. &
                      tr_pond_lvl) then
                 call tvd_icepack_other(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+1,trcrn(:,nt_alvl,n),trcrn(:,nt_apnd,n))
                 !,narrays+3+nt_alvl,narrays+3+nt_apnd)
              elseif (trcr_depend(it) == 2+nt_fbri) then
                  call tvd_icepack_other(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+2,trcrn(:,nt_fbri,n))!,narrays+3+nt_fbri)
              endif
           enddo
           narrays = narrays + ntrcr +3
  
        enddo                     ! n

        
        if (narr /= narrays .and. myrank == 0 ) write(ice_stderr,*)      &
            "Wrong number of arrays in transport bound call"

    end subroutine work_to_other3

!=======================================================================

    subroutine work_to_other_4 (ntrcr, narr, works)

        integer (kind=int_kind), intent(in) :: ntrcr, narr
  
        real (kind=dbl_kind), dimension(nx,narr), intent (inout) :: &
           works     ! work array
  
        ! local variables
  
        integer (kind=int_kind) :: &
           nt_alvl, nt_apnd, nt_fbri, nt_Tsfc, nt_qice, &
           &arr_alvl,arr_apnd,arr_fbri
  
        logical (kind=log_kind) :: &
           tr_pond_lvl, tr_pond_topo
  
        integer (kind=int_kind) ::      &
           i, n, it   , & ! counting indices
           narrays    , & ! counter for number of state variable arrays
         nt_qsno,k,j
  
        real (kind=dbl_kind) :: &
           rhos       , &
           Lfresh,small
        real (kind=dbl_kind),    dimension (:), allocatable ::       &
           swild
        character(len=*), parameter :: subname = '(work_to_other_4)'
  
        allocate ( swild(nx) )
        call icepack_query_tracer_flags(tr_pond_lvl_out=tr_pond_lvl, tr_pond_topo_out=tr_pond_topo)
        call icepack_query_tracer_indices(nt_alvl_out=nt_alvl, nt_apnd_out=nt_apnd, &
             nt_fbri_out=nt_fbri, nt_qsno_out=nt_qsno, nt_Tsfc_out=nt_Tsfc,nt_qice_out=nt_qice)
        call icepack_query_parameters(rhos_out=rhos, Lfresh_out=Lfresh)
        call icepack_warnings_flush(ice_stderr)
        if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
           file=__FILE__, line=__LINE__)

        small = 0.0000001
        swild(:) = c0
        narrays = 1
  
        do n=1, ncat
    
           do it = 1, ntrcr
              if (trcr_depend(it) == 0) then
                 call upwind_icepack_other4(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+1)
              elseif (trcr_depend(it) == 1) then
                 call upwind_icepack_other4(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+2)
              elseif (trcr_depend(it) == 2) then
               if (it >= nt_qsno .and. it < nt_qsno+nslyr) then
                     swild = trcrn(:,it,n) + rhos*Lfresh
                     call upwind_icepack_other4(works(:,narrays+3+it),swild,narrays+3+it,narrays+3)
               else
                     call upwind_icepack_other4(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+3)
               endif
              elseif (trcr_depend(it) == 2+nt_alvl) then
                call upwind_icepack_other4(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+1,trcrn(:,nt_alvl,n))!narrays+3+nt_alvl)
              elseif (trcr_depend(it) == 2+nt_apnd .and. &
                      tr_pond_topo) then
               call upwind_icepack_other4(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+1,trcrn(:,nt_apnd,n))!narrays+3+nt_apnd)
              elseif (trcr_depend(it) == 2+nt_apnd .and. &
                      tr_pond_lvl) then
                 call upwind_icepack_other4(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+1,trcrn(:,nt_alvl,n),trcrn(:,nt_apnd,n))
                 !,narrays+3+nt_alvl,narrays+3+nt_apnd)
              elseif (trcr_depend(it) == 2+nt_fbri) then
                  call upwind_icepack_other4(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+2,trcrn(:,nt_fbri,n))!,narrays+3+nt_fbri)
              endif
           enddo
           narrays = narrays + ntrcr +3
  
        enddo                     ! n

        
        if (narr /= narrays .and. myrank == 0 ) write(ice_stderr,*)      &
            "Wrong number of arrays in transport bound call"

    end subroutine work_to_other_4

!=======================================================================

    subroutine work_to_other_5 (ntrcr, narr, works)

        integer (kind=int_kind), intent(in) :: ntrcr, narr
  
        real (kind=dbl_kind), dimension(nx,narr), intent (inout) :: &
           works     ! work array
  
        ! local variables
  
        integer (kind=int_kind) :: &
           nt_alvl, nt_apnd, nt_fbri, nt_Tsfc, nt_qice, &
           &arr_alvl,arr_apnd,arr_fbri
  
        logical (kind=log_kind) :: &
           tr_pond_lvl, tr_pond_topo
  
        integer (kind=int_kind) ::      &
           i, n, it   , & ! counting indices
           narrays    , & ! counter for number of state variable arrays
         nt_qsno,k,j
  
        real (kind=dbl_kind) :: &
           rhos       , &
           Lfresh,small
        real (kind=dbl_kind),    dimension (:), allocatable ::       &
           swild
        character(len=*), parameter :: subname = '(work_to_other_5)'
  
        allocate ( swild(nx) )
        call icepack_query_tracer_flags(tr_pond_lvl_out=tr_pond_lvl, tr_pond_topo_out=tr_pond_topo)
        call icepack_query_tracer_indices(nt_alvl_out=nt_alvl, nt_apnd_out=nt_apnd, &
             nt_fbri_out=nt_fbri, nt_qsno_out=nt_qsno, nt_Tsfc_out=nt_Tsfc,nt_qice_out=nt_qice)
        call icepack_query_parameters(rhos_out=rhos, Lfresh_out=Lfresh)
        call icepack_warnings_flush(ice_stderr)
        if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
           file=__FILE__, line=__LINE__)

        small = 0.0000001
        swild(:) = c0
        narrays = 1
  
        do n=1, ncat
    
           do it = 1, ntrcr
              if (trcr_depend(it) == 0) then
                 call tvdu_icepack_other(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+1)
              elseif (trcr_depend(it) == 1) then
                 call tvdu_icepack_other(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+2)
              elseif (trcr_depend(it) == 2) then
               if (it >= nt_qsno .and. it < nt_qsno+nslyr) then
                     swild = trcrn(:,it,n) + rhos*Lfresh
                     call tvdu_icepack_other(works(:,narrays+3+it),swild,narrays+3+it,narrays+3)
               else
                     call tvdu_icepack_other(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+3)
               endif
              elseif (trcr_depend(it) == 2+nt_alvl) then
                call tvdu_icepack_other(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+1,trcrn(:,nt_alvl,n))!narrays+3+nt_alvl)
              elseif (trcr_depend(it) == 2+nt_apnd .and. &
                      tr_pond_topo) then
               call tvdu_icepack_other(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+1,trcrn(:,nt_apnd,n))!narrays+3+nt_apnd)
              elseif (trcr_depend(it) == 2+nt_apnd .and. &
                      tr_pond_lvl) then
                 call tvdu_icepack_other(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+1,trcrn(:,nt_alvl,n),trcrn(:,nt_apnd,n))
                 !,narrays+3+nt_alvl,narrays+3+nt_apnd)
              elseif (trcr_depend(it) == 2+nt_fbri) then
                  call tvdu_icepack_other(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+2,trcrn(:,nt_fbri,n))!,narrays+3+nt_fbri)
              endif
           enddo
           narrays = narrays + ntrcr +3
  
        enddo                     ! n

        
        if (narr /= narrays .and. myrank == 0 ) write(ice_stderr,*)      &
            "Wrong number of arrays in transport bound call"

    end subroutine work_to_other_5

!=======================================================================

    subroutine work_to_other_6 (ntrcr, narr, works)

        integer (kind=int_kind), intent(in) :: ntrcr, narr
  
        real (kind=dbl_kind), dimension(nx,narr), intent (inout) :: &
           works     ! work array
  
        ! local variables
  
        integer (kind=int_kind) :: &
           nt_alvl, nt_apnd, nt_fbri, nt_Tsfc, nt_qice, &
           &arr_alvl,arr_apnd,arr_fbri
  
        logical (kind=log_kind) :: &
           tr_pond_lvl, tr_pond_topo
  
        integer (kind=int_kind) ::      &
           i, n, it   , & ! counting indices
           narrays    , & ! counter for number of state variable arrays
         nt_qsno,k,j
  
        real (kind=dbl_kind) :: &
           rhos       , &
           Lfresh,small
        real (kind=dbl_kind),    dimension (:), allocatable ::       &
           swild
        character(len=*), parameter :: subname = '(work_to_other_5)'
  
        allocate ( swild(nx) )
        call icepack_query_tracer_flags(tr_pond_lvl_out=tr_pond_lvl, tr_pond_topo_out=tr_pond_topo)
        call icepack_query_tracer_indices(nt_alvl_out=nt_alvl, nt_apnd_out=nt_apnd, &
             nt_fbri_out=nt_fbri, nt_qsno_out=nt_qsno, nt_Tsfc_out=nt_Tsfc,nt_qice_out=nt_qice)
        call icepack_query_parameters(rhos_out=rhos, Lfresh_out=Lfresh)
        call icepack_warnings_flush(ice_stderr)
        if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
           file=__FILE__, line=__LINE__)

        small = 0.0000001
        swild(:) = c0
        narrays = 1
  
        do n=1, ncat
    
           do it = 1, ntrcr
              if (trcr_depend(it) == 0) then
                 call tvd_casulli_icepack_other(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+1)
              elseif (trcr_depend(it) == 1) then
                 call tvd_casulli_icepack_other(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+2)
              elseif (trcr_depend(it) == 2) then
               if (it >= nt_qsno .and. it < nt_qsno+nslyr) then
                     swild = trcrn(:,it,n) + rhos*Lfresh
                     call tvd_casulli_icepack_other(works(:,narrays+3+it),swild,narrays+3+it,narrays+3)
               else
                     call tvd_casulli_icepack_other(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+3)
               endif
              elseif (trcr_depend(it) == 2+nt_alvl) then
                call tvd_casulli_icepack_other(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+1,trcrn(:,nt_alvl,n))!narrays+3+nt_alvl)
              elseif (trcr_depend(it) == 2+nt_apnd .and. &
                      tr_pond_topo) then
               call tvd_casulli_icepack_other(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+1,trcrn(:,nt_apnd,n))!narrays+3+nt_apnd)
              elseif (trcr_depend(it) == 2+nt_apnd .and. &
                      tr_pond_lvl) then
                 call tvd_casulli_icepack_other(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+1,trcrn(:,nt_alvl,n),trcrn(:,nt_apnd,n))
                 !,narrays+3+nt_alvl,narrays+3+nt_apnd)
              elseif (trcr_depend(it) == 2+nt_fbri) then
                  call tvd_casulli_icepack_other(works(:,narrays+3+it),trcrn(:,it,n),narrays+3+it,narrays+2,trcrn(:,nt_fbri,n))!,narrays+3+nt_fbri)
              endif
           enddo
           narrays = narrays + ntrcr +3
  
        enddo                     ! n

        
        if (narr /= narrays .and. myrank == 0 ) write(ice_stderr,*)      &
            "Wrong number of arrays in transport bound call"

    end subroutine work_to_other_6

!=======================================================================

    subroutine re_momentum3 (ntrcr, narr, works, vicen_ini)

      use mice_module,         only:U_ice, V_ice,ice_cutoff,ice_tests
      integer (kind=int_kind), intent(in) :: ntrcr, narr

      real (kind=dbl_kind), dimension(nx,narr), intent (inout) :: &
         works     ! work array
      real (kind=dbl_kind), dimension(nx,ncat), intent (inout) :: &
         vicen_ini     ! init. vicen

      ! local variables

      integer (kind=int_kind) :: &
         nt_alvl, nt_apnd, nt_fbri, nt_Tsfc, nt_qice, &
         &arr_alvl,arr_apnd,arr_fbri

      logical (kind=log_kind) :: &
         tr_pond_lvl, tr_pond_topo

      integer (kind=int_kind) ::      &
         i, n, it   , & ! counting indices
         narrays    , & ! counter for number of state variable arrays
         nt_qsno,k,j,m, n1,n2,ie,trc_d,id,id2,id3

      real (kind=dbl_kind) :: &
         rhos       , &
         Lfresh, small, sumvice, sumvice_i, utmp, vtmp, momentum_ui, momentum_vi, &
         &suma, puny, sumaice
      
         real(rkind) :: swild10(3,2), swild(3), swild2(3,2),fluxchan1(nx),fluxchan2(nx)
      character(len=*), parameter :: subname = '(re_momentum3)'

      call icepack_query_tracer_indices(nt_alvl_out=nt_alvl, nt_apnd_out=nt_apnd, &
           nt_fbri_out=nt_fbri, nt_qsno_out=nt_qsno, nt_Tsfc_out=nt_Tsfc,nt_qice_out=nt_qice)
      call icepack_query_parameters(rhos_out=rhos, Lfresh_out=Lfresh,puny_out = puny)
      call icepack_warnings_flush(ice_stderr)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__, line=__LINE__)

      small = 0.0000001
      sumvice = 0
      fluxchan1(:) = c0
      fluxchan2(:) = c0
      do i=1,nx_nh
      fluxchan1(i) = c0
      fluxchan2(i) = c0
      suma = c0
         do j=1,nne(i)
            ie=indel(j,i)
            suma=suma+area(ie)/3 !approx for quad
            id=iself(j,i)
            id3=nxq(i34(ie)-1,id,i34(ie)) !adjacent side index
            id2=nxq(i34(ie)-2,id,i34(ie)) !adjacent side index

            do m=1,i34(ie) !side
               n1=isidenode(1,elside(m,ie))
               n2=isidenode(2,elside(m,ie))
               if(ics==1) then
                   swild10(m,1) = 0.5*(uvel(n1)+uvel(n2))
                   swild10(m,2) = 0.5*(vvel(n1)+vvel(n2))
                   swild2(m,1) = uvel(elnode(m,ie))
                   swild2(m,2) = vvel(elnode(m,ie))
               else
                  swild10(m,1) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),eframe(:,1,ie))+vvel(n1)*dot_product(pframe(:,2,n1),eframe(:,1,ie))&
                              &+uvel(n2)*dot_product(pframe(:,1,n2),eframe(:,1,ie))+vvel(n2)*dot_product(pframe(:,2,n2),eframe(:,1,ie)))
                  swild10(m,2) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),eframe(:,2,ie))+vvel(n1)*dot_product(pframe(:,2,n1),eframe(:,2,ie))&
                              &+uvel(n2)*dot_product(pframe(:,1,n2),eframe(:,2,ie))+vvel(n2)*dot_product(pframe(:,2,n2),eframe(:,2,ie)))                   
                  swild2(m,1) = uvel(elnode(m,ie))*dot_product(pframe(:,1,elnode(m,ie)),eframe(:,1,ie))+vvel(elnode(m,ie))*dot_product(pframe(:,2,elnode(m,ie)),eframe(:,1,ie))
                  swild2(m,2) = uvel(elnode(m,ie))*dot_product(pframe(:,1,elnode(m,ie)),eframe(:,2,ie))+vvel(elnode(m,ie))*dot_product(pframe(:,2,elnode(m,ie)),eframe(:,2,ie))
               endif
               !swild10(m,1) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),sframe2(:,1,elside(m,ie)))+vvel(n1)*dot_product(pframe(:,2,n1),sframe2(:,1,elside(m,ie)))&
               !               &+uvel(n2)*dot_product(pframe(:,1,n2),sframe2(:,1,elside(m,ie)))+vvel(n2)*dot_product(pframe(:,2,n2),sframe2(:,1,elside(m,ie))))
               !swild10(m,2) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),sframe2(:,2,elside(m,ie)))+vvel(n1)*dot_product(pframe(:,2,n1),sframe2(:,2,elside(m,ie)))&
               !               &+uvel(n2)*dot_product(pframe(:,1,n2),sframe2(:,2,elside(m,ie)))+vvel(n2)*dot_product(pframe(:,2,n2),sframe2(:,2,elside(m,ie)))) 
               
            enddo
         utmp=sum(swild2(1:i34(ie),1))/3 !vel @ centroid
         vtmp=sum(swild2(1:i34(ie),2))/3 !vel @ centroid
         

         trc_d = 3
         do k = 1, ncat
            n1 = flux_matrix1(j,i,trc_d) 
            fluxchan1(i) = fluxchan1(i) + flux_matrix2(j,i,trc_d)*uvel(n1)
            fluxchan2(i) = fluxchan2(i) + flux_matrix2(j,i,trc_d)*vvel(n1)
            n2 = flux_matrix3(j,i,trc_d) 
            fluxchan1(i) = fluxchan1(i) + flux_matrix4(j,i,trc_d)*uvel(n2)
            fluxchan2(i) = fluxchan2(i) + flux_matrix4(j,i,trc_d)*vvel(n2)
            trc_d = trc_d + ntrcr + 3
         enddo
         enddo
          trc_d = 3
          sumvice = 0
          sumaice = 0
         do k = 1, ncat
         sumvice = sumvice + works(i,trc_d)
         sumaice = sumaice + works(i,trc_d-1)
         trc_d = trc_d + ntrcr + 3
         enddo
         sumvice_i = sum(vicen_ini(i,:))
         momentum_ui = sumvice_i * u_ice(i)
         momentum_vi = sumvice_i * v_ice(i)
         if(sumvice>puny) then
            u_ice(i) = (momentum_ui + fluxchan1(i)/suma)/sumvice
            v_ice(i) = (momentum_vi + fluxchan2(i)/suma)/sumvice
         else
            u_ice(i) = 0
            v_ice(i) = 0
         endif
         if(isbnd(1,i)/=0) then !b.c. (including open)
            U_ice(i)=0; V_ice(i)=0
         endif
        !if(sumaice<=ice_cutoff.or.sumvice<=ice_cutoff) then !no ice
        if(sumaice<=small.or.sumvice<=small) then !no ice
            U_ice(i)=uocn(i); V_ice(i)=vocn(i)
            !U_ice(i)=0; V_ice(i)=0
        endif
        if(sqrt(U_ice(i)**2+v_ice(i)**2)>2) then
         utmp = U_ice(i)
         vtmp = V_ice(i)
         U_ice(i)=2*utmp/sqrt(utmp**2+vtmp**2)
         V_ice(i)=2*vtmp/sqrt(utmp**2+vtmp**2)
         endif  
        if(ice_tests==1) then
          U_ice(i) = 1
          V_ice(i) = 0
        endif      
      enddo
      call exchange_p2d(u_ice)
      call exchange_p2d(v_ice)
  end subroutine re_momentum3

    !=======================================================================

    subroutine re_momentum (ntrcr, narr, works, vicen_ini)

      use mice_module,         only:U_ice, V_ice,ice_cutoff,ice_tests
      integer (kind=int_kind), intent(in) :: ntrcr, narr

      real (kind=dbl_kind), dimension(nx,narr), intent (inout) :: &
         works     ! work array
      real (kind=dbl_kind), dimension(nx,ncat), intent (inout) :: &
         vicen_ini     ! init. vicen

      ! local variables

      integer (kind=int_kind) :: &
         nt_alvl, nt_apnd, nt_fbri, nt_Tsfc, nt_qice, &
         &arr_alvl,arr_apnd,arr_fbri

      logical (kind=log_kind) :: &
         tr_pond_lvl, tr_pond_topo

      integer (kind=int_kind) ::      &
         i, n, it   , & ! counting indices
         narrays    , & ! counter for number of state variable arrays
         nt_qsno,k,j,m, n1,n2,ie,trc_d,id,id2,id3

      real (kind=dbl_kind) :: &
         rhos       , &
         Lfresh, small, sumvice, sumvice_i, utmp, vtmp, momentum_ui, momentum_vi, &
         &suma, puny, sumaice
      
         real(rkind) :: swild10(3,2), swild(3), swild2(3,2),fluxchan1(nx),fluxchan2(nx)
      character(len=*), parameter :: subname = '(re_momentum)'

      call icepack_query_tracer_indices(nt_alvl_out=nt_alvl, nt_apnd_out=nt_apnd, &
           nt_fbri_out=nt_fbri, nt_qsno_out=nt_qsno, nt_Tsfc_out=nt_Tsfc,nt_qice_out=nt_qice)
      call icepack_query_parameters(rhos_out=rhos, Lfresh_out=Lfresh,puny_out = puny)
      call icepack_warnings_flush(ice_stderr)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__, line=__LINE__)

      small = 0.0000001
      sumvice = 0
      fluxchan1(:) = c0
      fluxchan2(:) = c0
      do i=1,nx_nh
      fluxchan1(i) = c0
      fluxchan2(i) = c0
      suma = c0
         do j=1,nne(i)
            ie=indel(j,i)
            suma=suma+area(ie)/3 !approx for quad
            id=iself(j,i)
            id3=nxq(i34(ie)-1,id,i34(ie)) !adjacent side index
            id2=nxq(i34(ie)-2,id,i34(ie)) !adjacent side index

            do m=1,i34(ie) !side
               n1=isidenode(1,elside(m,ie))
               n2=isidenode(2,elside(m,ie))
               if(ics==1) then
                  swild10(m,1) = 0.5*(uvel(n1)+uvel(n2))
                  swild10(m,2) = 0.5*(vvel(n1)+vvel(n2))
                  swild2(m,1) = uvel(elnode(m,ie))
                  swild2(m,2) = vvel(elnode(m,ie))
               else
                  swild10(m,1) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),eframe(:,1,ie))+vvel(n1)*dot_product(pframe(:,2,n1),eframe(:,1,ie))&
                              &+uvel(n2)*dot_product(pframe(:,1,n2),eframe(:,1,ie))+vvel(n2)*dot_product(pframe(:,2,n2),eframe(:,1,ie)))
                  swild10(m,2) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),eframe(:,2,ie))+vvel(n1)*dot_product(pframe(:,2,n1),eframe(:,2,ie))&
                              &+uvel(n2)*dot_product(pframe(:,1,n2),eframe(:,2,ie))+vvel(n2)*dot_product(pframe(:,2,n2),eframe(:,2,ie)))                   
                  swild2(m,1) = uvel(elnode(m,ie))*dot_product(pframe(:,1,elnode(m,ie)),eframe(:,1,ie))+vvel(elnode(m,ie))*dot_product(pframe(:,2,elnode(m,ie)),eframe(:,1,ie))
                  swild2(m,2) = uvel(elnode(m,ie))*dot_product(pframe(:,1,elnode(m,ie)),eframe(:,2,ie))+vvel(elnode(m,ie))*dot_product(pframe(:,2,elnode(m,ie)),eframe(:,2,ie))
               endif
               !swild10(m,1) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),sframe2(:,1,elside(m,ie)))+vvel(n1)*dot_product(pframe(:,2,n1),sframe2(:,1,elside(m,ie)))&
               !               &+uvel(n2)*dot_product(pframe(:,1,n2),sframe2(:,1,elside(m,ie)))+vvel(n2)*dot_product(pframe(:,2,n2),sframe2(:,1,elside(m,ie))))
               !swild10(m,2) = 0.5*(uvel(n1)*dot_product(pframe(:,1,n1),sframe2(:,2,elside(m,ie)))+vvel(n1)*dot_product(pframe(:,2,n1),sframe2(:,2,elside(m,ie)))&
               !               &+uvel(n2)*dot_product(pframe(:,1,n2),sframe2(:,2,elside(m,ie)))+vvel(n2)*dot_product(pframe(:,2,n2),sframe2(:,2,elside(m,ie)))) 
               
            enddo
         utmp=sum(swild2(1:i34(ie),1))/3 !vel @ centroid
         vtmp=sum(swild2(1:i34(ie),2))/3 !vel @ centroid
         

         trc_d = 3
         do k = 1, ncat
            fluxchan1(i) = fluxchan1(i) + (flux_matrix1(j,i,trc_d)*utmp+flux_matrix2(j,i,trc_d)*swild10(id3,1))/2.d0
            fluxchan2(i) = fluxchan2(i) + (flux_matrix1(j,i,trc_d)*vtmp+flux_matrix2(j,i,trc_d)*swild10(id3,2))/2.d0
            fluxchan1(i) = fluxchan1(i) + (flux_matrix3(j,i,trc_d)*utmp+flux_matrix4(j,i,trc_d)*swild10(id2,1))/2.d0
            fluxchan2(i) = fluxchan2(i) + (flux_matrix3(j,i,trc_d)*vtmp+flux_matrix4(j,i,trc_d)*swild10(id2,2))/2.d0
            trc_d = trc_d + ntrcr + 3
         enddo
         enddo
          trc_d = 3
          sumvice = 0
          sumaice = 0
         do k = 1, ncat
         sumvice = sumvice + works(i,trc_d)
         sumaice = sumaice + works(i,trc_d-1)
         trc_d = trc_d + ntrcr + 3
         enddo
         sumvice_i = sum(vicen_ini(i,:))
         momentum_ui = sumvice_i * u_ice(i)
         momentum_vi = sumvice_i * v_ice(i)
         if(sumvice>puny) then
            u_ice(i) = (momentum_ui + fluxchan1(i)/suma)/sumvice
            v_ice(i) = (momentum_vi + fluxchan2(i)/suma)/sumvice
         else
            u_ice(i) = 0
            v_ice(i) = 0
         endif
         if(isbnd(1,i)/=0) then !b.c. (including open)
            U_ice(i)=0; V_ice(i)=0
         endif
        !if(sumaice<=ice_cutoff.or.sumvice<=ice_cutoff) then !no ice
        if(sumaice<=small.or.sumvice<=small) then !no ice
            U_ice(i)=uocn(i); V_ice(i)=vocn(i)
            !U_ice(i)=0; V_ice(i)=0
        endif
        if(sqrt(U_ice(i)**2+v_ice(i)**2)>2) then
         utmp = U_ice(i)
         vtmp = V_ice(i)
         U_ice(i)=2*utmp/sqrt(utmp**2+vtmp**2)
         V_ice(i)=2*vtmp/sqrt(utmp**2+vtmp**2)
         endif       
        if(ice_tests==1) then
          U_ice(i) = 1
          V_ice(i) = 0
        endif 
      enddo
      call exchange_p2d(u_ice)
      call exchange_p2d(v_ice)
  end subroutine re_momentum

    !=======================================================================

    subroutine re_momentum2 (ntrcr, narr, works, vicen_ini)

      use mice_module,         only:U_ice, V_ice,ice_cutoff,ice_tests
      integer (kind=int_kind), intent(in) :: ntrcr, narr

      real (kind=dbl_kind), dimension(nx,narr), intent (inout) :: &
         works     ! work array
      real (kind=dbl_kind), dimension(nx,ncat), intent (inout) :: &
         vicen_ini     ! init. vicen

      ! local variables

      integer (kind=int_kind) :: &
         nt_alvl, nt_apnd, nt_fbri, nt_Tsfc, nt_qice, &
         &arr_alvl,arr_apnd,arr_fbri

      logical (kind=log_kind) :: &
         tr_pond_lvl, tr_pond_topo

      integer (kind=int_kind) ::      &
         i, n, it   , & ! counting indices
         narrays    , & ! counter for number of state variable arrays
         nt_qsno,k,j,m, n1,n2,ie,trc_d,id,id2,id3

      real (kind=dbl_kind) :: &
         rhos       , &
         Lfresh, small, sumvice, sumvice_i, utmp, vtmp, momentum_ui, momentum_vi, &
         &suma, puny, sumaice
      
         real(rkind) :: swild10(3,2), swild(3), swild2(3,2),fluxchan1(nx),fluxchan2(nx)
      character(len=*), parameter :: subname = '(re_momentum2)'

      call icepack_query_tracer_indices(nt_alvl_out=nt_alvl, nt_apnd_out=nt_apnd, &
           nt_fbri_out=nt_fbri, nt_qsno_out=nt_qsno, nt_Tsfc_out=nt_Tsfc,nt_qice_out=nt_qice)
      call icepack_query_parameters(rhos_out=rhos, Lfresh_out=Lfresh,puny_out = puny)
      call icepack_warnings_flush(ice_stderr)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         file=__FILE__, line=__LINE__)

      small = 0.0000001
      sumvice = 0
      fluxchan1(:) = c0
      fluxchan2(:) = c0
      do i=1,nx_nh
      fluxchan1(i) = c0
      fluxchan2(i) = c0
      suma = c0
         do j=1,nne(i)
            ie=indel(j,i)
            suma=suma+area(ie)/3 !approx for quad
            id=iself(j,i)
            id3=nxq(i34(ie)-1,id,i34(ie)) !adjacent side index
            id2=nxq(i34(ie)-2,id,i34(ie)) !adjacent side index


         trc_d = 3
         do k = 1, ncat
            fluxchan1(i) = fluxchan1(i) +sum(flux_matrix0(1:3,j,i,trc_d)*u_ice(elnode(1:3,ie)))
            fluxchan2(i) = fluxchan2(i) +sum(flux_matrix0(1:3,j,i,trc_d)*v_ice(elnode(1:3,ie)))
            trc_d = trc_d + ntrcr + 3
         enddo
         enddo
          trc_d = 3
          sumvice = 0
          sumaice = 0
         do k = 1, ncat
         sumvice = sumvice + works(i,trc_d)
         sumaice = sumaice + works(i,trc_d-1)
         trc_d = trc_d + ntrcr + 3
         enddo
         sumvice_i = sum(vicen_ini(i,:))
         momentum_ui = sumvice_i * u_ice(i)
         momentum_vi = sumvice_i * v_ice(i)
         if(sumvice>puny) then
            u_ice(i) = (momentum_ui + fluxchan1(i))/sumvice
            v_ice(i) = (momentum_vi + fluxchan2(i))/sumvice
         else
            u_ice(i) = 0
            v_ice(i) = 0
         endif
         if(isbnd(1,i)/=0) then !b.c. (including open)
            U_ice(i)=0; V_ice(i)=0
         endif
        if(sumaice<=small.or.sumvice<=small) then !no ice
            U_ice(i)=uocn(i); V_ice(i)=vocn(i)
            !U_ice(i)=0; V_ice(i)=0
        endif
        if(sqrt(U_ice(i)**2+v_ice(i)**2)>2) then
         utmp = U_ice(i)
         vtmp = V_ice(i)
         U_ice(i)=2*utmp/sqrt(utmp**2+vtmp**2)
         V_ice(i)=2*vtmp/sqrt(utmp**2+vtmp**2)
         endif
        if(ice_tests==1) then
          U_ice(i) = 1
          V_ice(i) = 0
        endif
      enddo
      call exchange_p2d(u_ice)
      call exchange_p2d(v_ice)
  end subroutine re_momentum2

!=======================================================================

    module subroutine cut_off_icepack

      use icepack_intfc,         only: icepack_compute_tracers
      use icepack_intfc,         only: icepack_aggregate
      use icepack_intfc,         only: icepack_init_trcr
      use icepack_intfc,         only: icepack_sea_freezing_temperature
      use icepack_therm_shared,  only: calculate_Tin_from_qin
      use icepack_mushy_physics, only: icepack_mushy_temperature_mush
      use icepack_mushy_physics, only: liquidus_temperature_mush
      use icepack_mushy_physics, only: enthalpy_mush,enthalpy_of_melting

      ! local variables

      real (kind=dbl_kind), dimension(nilyr) :: &
         qin            , & ! ice enthalpy (J/m3)
         qin_max        , & ! maximum ice enthalpy (J/m3)
         zTin           , & ! initial ice temperature
         tmltz0

      real (kind=dbl_kind), dimension(nslyr) :: &
         qsn            , & ! snow enthalpy (J/m3)
         zTsn               ! initial snow temperature
      integer (kind=int_kind) ::      &
         i, n, k, it    , & ! counting indices
         narrays        , & ! counter for number of state variable arrays
         icells         , & ! number of ocean/ice cells
         ktherm         , &
         ntrcr

       real (kind=dbl_kind), dimension(ncat) :: &
         aicecat

      real (kind=dbl_kind) ::         &
         rhos,       Lfresh,          &
         cp_ice,     cp_ocn,          &
         qrd_snow,   qrd_ice,         &
         Tsfc,       exc,             &
         depressT,   Tmin,            &
         T_air_C,    hice,            &
         puny,       Tsmelt,          &
         small,      rhoi,            &
         hpnd_max,   Tmax          
         
         

      logical (kind=log_kind) :: tr_brine, tr_lvl, flag_snow, flag_cold_ice, flag_warm_ice, &
                                 tr_pond_topo, tr_pond_lvl, tr_FY, tr_iage
      integer (kind=int_kind) :: nt_Tsfc, nt_qice, nt_qsno, nt_sice
      integer (kind=int_kind) :: nt_fbri, nt_alvl, nt_vlvl, nt_apnd, nt_hpnd, nt_ipnd, nt_FY, nt_iage

      character(len=*), parameter :: subname = '(cut_off_icepack)'

      call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
      call icepack_query_tracer_flags(tr_brine_out=tr_brine, tr_lvl_out=tr_lvl,      &
        tr_pond_topo_out=tr_pond_topo,                                               &
        tr_pond_lvl_out=tr_pond_lvl, tr_FY_out=tr_FY, tr_iage_out=tr_iage            )
      call icepack_query_tracer_indices( nt_Tsfc_out=nt_Tsfc, nt_qice_out=nt_qice,   &
        nt_qsno_out=nt_qsno, nt_sice_out=nt_sice, nt_FY_out=nt_FY,                   &
        nt_fbri_out=nt_fbri, nt_alvl_out=nt_alvl, nt_vlvl_out=nt_vlvl,               &         
        nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, nt_ipnd_out=nt_ipnd,               &
        nt_iage_out=nt_iage                                                          )
      call icepack_query_parameters(rhos_out=rhos, rhoi_out=rhoi, Lfresh_out=Lfresh, &
                                    cp_ice_out=cp_ice, cp_ocn_out=cp_ocn)
      call icepack_query_parameters(depressT_out=depressT, puny_out=puny, &
        Tsmelt_out=Tsmelt, ktherm_out=ktherm)
      call icepack_warnings_flush(ice_stderr)

      small = puny
      !small=0.005_dbl_kind        
      Tmin  = -100.0_dbl_kind

      ! only for bl99 and mushy thermodynamics

      ! Here we should implement some conditions to check the tracers
      ! when ice is present, particularly enthalpy, surface temperature
      ! and salinity.

      do n = 1, ncat                      ! For each thickness cathegory
           do i = 1, nx                   ! For each grid point

            call icepack_init_trcr(T_air(i),    Tf(i),         &
                                   !salinz(i,:), Tmltz(i,:),    &
                                  trcrn(i,nt_sice:nt_sice+nilyr-1,n), Tmltz(i,:),    &
                                   Tsfc,                       &
                                   nilyr,        nslyr,        &
                                   qin   (:),    qsn  (:))
            call icepack_warnings_flush(ice_stderr)
            if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
            file=__FILE__, line=__LINE__)  
                do k = 1, nilyr
                   if (ktherm == 2) then
                      tmltz0(k) = liquidus_temperature_mush(trcrn(i,nt_sice+k-1,n))
                      !qin(k)=enthalpy_mush(tmltz0(k), trcrn(i,nt_sice+k-1,n))
                      qin_max(k) = enthalpy_of_melting(trcrn(i,nt_sice+k-1,n)) - c1
                      !qin_max(k) = enthalpy_mush((Tmltz(i,k)-0.1), trcrn(i,nt_sice+k-1,n))
                   else
                      tmltz0(k) = -trcrn(i,nt_sice+k-1,n) * depressT
                      qin_max(k) = rhoi * cp_ocn * (tmltz0(k) - 0.1_dbl_kind) 
                      !qin_max(:) = rhoi * cp_ocn * (Tmltz(i,k) - 0.1_dbl_kind)  
                   endif
                enddo
                ! Correct qin profile for melting temperatures

                ! Maximum ice enthalpy 

                !qin_max(:) = rhoi * cp_ocn * (tmltz0(:) - 60_dbl_kind) 
                !qin_max(:) = rhoi * cp_ocn * (Tmltz(i,:) - 0.1_dbl_kind)  
                if (vicen(i,n) > small .and. aicen(i,n) > small) then

                    ! Condition on surface temperature
                    if (trcrn(i,nt_Tsfc,n) > Tsmelt .or. trcrn(i,nt_Tsfc,n) < Tmin) then
                        trcrn(i,nt_Tsfc,n) = Tsfc
                    endif

                    ! Condition on ice enthalpy

                    flag_warm_ice = .false.
                    flag_cold_ice = .false.
                    flag_snow     = .false.

                      do k = 1, nilyr  ! Check for problems

                         if (ktherm == 2) then
                            zTin(k) = icepack_mushy_temperature_mush(trcrn(i,nt_qice+k-1,n),trcrn(i,nt_sice+k-1,n))
                         else
                            zTin(k) = calculate_Tin_from_qin(trcrn(i,nt_qice+k-1,n),Tmltz(i,k))
                         endif

                         if (zTin(k) <  Tmin      ) flag_cold_ice = .true.
                         !if (zTin(k) >= Tmltz(i,k)) flag_warm_ice = .true.
                         if (zTin(k) >= tmltz0(k)) flag_warm_ice = .true.

                      enddo !nilyr

                if (flag_cold_ice) then
                    trcrn(i,nt_Tsfc,n) = Tsfc

                    do k = 1, nilyr
                         trcrn(i,nt_qice+k-1,n) = min(qin_max(k), qin(k))
                    enddo        ! nilyr

                    if (vsnon(i,n) > small) then ! Only if there is snow
                                                 ! on top of the sea ice
                        do k = 1, nslyr
                              trcrn(i,nt_qsno+k-1,n) = qsn(k)
                        enddo   
                    else                        ! No snow 
                        trcrn(i,nt_qsno:nt_qsno+nslyr-1,n) = c0
                    endif

                end if            ! flag cold ice

                if (flag_warm_ice) then         ! This sea ice should have melted already 

                    aicen(i,n)   = c0
                    vicen(i,n)   = c0
                    vsnon(i,n)   = c0
                    trcrn(i,:,n) = c0
                    trcrn(i,nt_Tsfc,n) = Tf(i)

                    !do k = 1, nilyr
                    !    trcrn(i,nt_qice+k-1,n) = qin_max(k) -c1 !'Corrected quantities'
                    !enddo        ! nilyr

                     !trcrn(i,nt_Tsfc,n) = Tsfc

                end if

                    if (vsnon(i,n) > small) then

                        flag_snow = .false.

                        do k = 1, nslyr  
                           if (trcrn(i,nt_qsno+k-1,n) >= -rhos*Lfresh) flag_snow = .true.
                           zTsn(k) = (Lfresh + trcrn(i,nt_qsno+k-1,n)/rhos)/cp_ice
                           if (zTsn(k) < Tmin) flag_snow = .true.
                        end do 
                        !flag_snow = .false.

                        if (flag_snow) then
                            trcrn(i,nt_Tsfc,n) = Tsfc
                            do k = 1, nslyr
                                 trcrn(i,nt_qsno+k-1,n) = qsn(k)
                            enddo        ! nslyr
                        endif            ! flag snow
                    endif ! vsnon(i,n) > c0

                else

                    aicen(i,n)   = c0
                    vicen(i,n)   = c0
                    vsnon(i,n)   = c0
                    trcrn(i,:,n) = c0
                    trcrn(i,nt_Tsfc,n) = Tf(i)

                endif

           enddo   ! nx
      enddo        ! ncat

      ! Melt ponds and level ice cutoff

      do n = 1, ncat
         do i = 1, nx
            if (aicen(i,n) > c0) then
               hpnd_max = 0.9_dbl_kind * vicen(i,n) / aicen(i,n)
      !         ! Sea ice age
      !         if (tr_iage) then
      !             if (trcrn(i,nt_iage,n) < 0.000001_dbl_kind) trcrn(i,nt_iage,n) = c0
      !         end if
      !         ! First year ice fraction
      !         if (tr_FY) then
      !             if (trcrn(i,nt_FY,n) < 0.000001_dbl_kind) trcrn(i,nt_FY,n) = c0
      !             if (trcrn(i,nt_FY,n) > c1) trcrn(i,nt_FY,n) = c1
      !         end if
      !         ! Level ice
      !         if (tr_lvl) then
      !             if (trcrn(i,nt_alvl,n) > aicen(i,n)) then
      !                 trcrn(i,nt_alvl,n) = aicen(i,n)
      !             elseif (trcrn(i,nt_alvl,n) < 0.000001_dbl_kind) then
      !                 trcrn(i,nt_alvl,n) = c0
      !             endif
      !             if (trcrn(i,nt_vlvl,n) < 0.000001_dbl_kind .or. trcrn(i,nt_alvl,n) < 0.000001_dbl_kind) trcrn(i,nt_vlvl,n) = c0
      !             if (trcrn(i,nt_vlvl,n) > vicen(i,n)) trcrn(i,nt_vlvl,n) = vicen(i,n)
      !         end if
               ! Topo and level melt pond parameterization
               if (tr_pond_topo .or. tr_pond_lvl) then
                   if (trcrn(i,nt_apnd,n) > c1) then
                       trcrn(i,nt_apnd,n) = c1
                   elseif (trcrn(i,nt_apnd,n) < 0.000001_dbl_kind) then
                       trcrn(i,nt_apnd,n) = c0
                   endif
                   if (trcrn(i,nt_hpnd,n) < 0.000001_dbl_kind .or. trcrn(i,nt_apnd,n) < 0.000001_dbl_kind) trcrn(i,nt_hpnd,n) = c0
                   if (trcrn(i,nt_hpnd,n) > hpnd_max) trcrn(i,nt_hpnd,n) = hpnd_max
                   if (trcrn(i,nt_ipnd,n) < 0.000001_dbl_kind .or. trcrn(i,nt_apnd,n) < 0.000001_dbl_kind) trcrn(i,nt_ipnd,n) = c0
               end if
               ! Dynamic salt
               if (tr_brine) then
                   if (trcrn(i,nt_fbri,n) < 0.000001_dbl_kind) trcrn(i,nt_fbri,n) = c0
               endif
            endif
         enddo
      enddo

      !do i = 1, nx
      !   aice(i) = c0
      !   vice(i) = c0
      !   vsno(i) = c0
      !   do it = 1, ntrcr
      !      trcr(i,it) = c0
      !   enddo
      !   call icepack_aggregate (ncat,                    &
      !                          aicen(i,:),               &
      !                          trcrn(i,1:ntrcr,:),       &
      !                          vicen(i,:),               &
      !                          vsnon(i,:),               &
      !                          aice (i),                 &
      !                          trcr (i,1:ntrcr),         &
      !                          vice (i),                 &
      !                          vsno (i),                 &
      !                          aice0(i),                 &
      !                          ntrcr,                    &
      !                          trcr_depend  (1:ntrcr),   &
      !                          trcr_base    (1:ntrcr,:), &
      !                          n_trcr_strata(1:ntrcr),   &
      !                          nt_strata    (1:ntrcr,:))
      !end do

      call icepack_warnings_flush(ice_stderr)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__, line=__LINE__)

  end subroutine cut_off_icepack


end submodule icedrv_advection


