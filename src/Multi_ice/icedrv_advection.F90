!=======================================================================
!
! This submodule contains the subroutines 
! for the advection of sea ice tracers
!
! Author: L. Zampieri ( lorenzo.zampieri@awi.de )
!  Modified by Qian Wang to apply to SCHISM
!
! Main driver: tracer_advection_icepack, which 
! calls fct_solve_icepack (sequence of lower- and higher-order FCT solver for
! each tracer)
!=======================================================================

submodule (icedrv_main) icedrv_advection

    use icedrv_kinds
    use icedrv_constants
    use icedrv_system,      only: icedrv_system_abort
    use schism_msgp,        only: parallel_abort,parallel_finalize,exchange_p2d
    use icepack_intfc,      only: icepack_warnings_flush,         &
                                  icepack_warnings_aborted,       &
                                  icepack_query_tracer_indices,   &
                                  icepack_query_tracer_flags,     &
                                  icepack_query_parameters,       &
                                  icepack_query_tracer_sizes
    use schism_glbl, only: rkind,nea,np,npa,elnode,nnp,indnd,time_stamp,rnday, &
   &fdb,lfdb,xnd,ynd,iplg,idry,i34,isbnd

    implicit none

    real(kind=dbl_kind), allocatable, dimension(:)   :: &
         d_tr,        trl,                              &
         rhs_tr,      rhs_trdiv,                        &
         icepplus,    icepminus,                        &
         mass_matrix, newwork                   

    real(kind=dbl_kind), allocatable, dimension(:,:) :: &
         icefluxes                   

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
      
        !type(t_mesh), intent(in), target :: mesh
          
        ! Initialization of arrays necessary to implement FCT algorithm
        allocate(trl(nx))   ! low-order solutions
        allocate(d_tr(nx))  ! increments of high
                            ! order solutions
        allocate(icefluxes(nx_elem, 3))
        allocate(icepplus(nx), icepminus(nx))
        allocate(rhs_tr(nx),  rhs_trdiv(nx))
        allocate(newwork(nx)) 
        !allocate(mass_matrix(sum(nn_num(1:nx_nh))))

      
        trl(:)    = c0
        d_tr(:)   = c0
        rhs_tr(:) = c0
        rhs_trdiv(:)   = c0
        icefluxes(:,:) = c0
        icepplus(:)    = c0
        icepminus(:)   = c0
        newwork(:)     = c0
        !mass_matrix(:) = c0
      
        ! Fill in  the mass matrix
        !call fill_mass_matrix_icepack
        !mass_matrix=ice_matrix
        if (myrank==0) write(*,*) 'Icepack FCT is initialized'
    
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
        real   (kind=dbl_kind), allocatable, dimension(:) :: tmax, tmin
        real   (kind=dbl_kind)                            :: vol, flux, ae, gamma
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
       
        ! Auxiliary elemental operator (mass matrix- lumped mass matrix)
       
        icoef = 1
        icefluxes(:,:) = c0

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
           tmax(row) = max(tmax(row),trl(row))
           tmin(row) = minval(trl(indnd(1:n,row)))
           tmin(row) = min(tmin(row),trl(row))
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
               if ( flux >=c0 ) ae = min(ae, icepplus(n))
               if ( flux < c0 ) ae = min(ae, icepminus(n))
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
       
          deallocate(tmin, tmax)
    
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
        real   (kind=dbl_kind)    :: c_1, c_2, c_3, c_4, c_x, entries2(3),entries3(3)
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
      
           ! Derivatives
           dx  = bafux(:,elem)
           dy  = bafuy(:,elem)
           vol = voltriangle(elem)
           um  = sum(uvel(elnodes))
           vm  = sum(vvel(elnodes))
      
           ! This is exact computation (no assumption of u=const 
           ! on elements used in the standard version)
           c_1 = (um*um+sum(uvel(elnodes)*uvel(elnodes))) / 12.0_dbl_kind
           c_2 = (vm*vm+sum(vvel(elnodes)*vvel(elnodes))) / 12.0_dbl_kind
           c_3 = (um*vm+sum(vvel(elnodes)*uvel(elnodes))) / 12.0_dbl_kind
           c_4 = sum(dx*uvel(elnodes)+dy*vvel(elnodes))
      
           do n = 1, 3
              row = elnodes(n)
      
              do q = 1, 3
                 entries(q)  = vol*dt_dyn*((c1-p5*dt_dyn*c_4)*(dx(n)*(um+uvel(elnodes(q)))+ &
                               dy(n)*(vm+vvel(elnodes(q))))/12.0_dbl_kind                 - &
                               p5*dt_dyn*(c_1*dx(n)*dx(q)+c_2*dy(n)*dy(q)+c_3*(dx(n)*dy(q)+dx(q)*dy(n))))
                 entries2(q) = p5*dt_dyn*(dx(n)*(um+uvel(elnodes(q)))                     + &
                               dy(n)*(vm+vvel(elnodes(q)))-dx(q)*(um+uvel(row))      - &
                               dy(q)*(vm+vvel(row)))
                 entries3(q)=  vol*dt_dyn*(((dx(n)*(um+uvel(elnodes(q))))+ dy(n)*(vm+vvel(elnodes(q))))/12.0_dbl_kind &
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
                                   nt,   nt1,    k,   n ,ktherm
        integer (kind=int_kind) :: nt_Tsfc, nt_qice, nt_qsno,                          & 
                                   nt_sice, nt_fbri, nt_iage, nt_FY, nt_alvl, nt_vlvl, &
                                   nt_apnd, nt_hpnd, nt_ipnd, nt_bgc_Nit, nt_bgc_S
        logical (kind=log_kind) :: tr_pond_topo, tr_pond_lvl, tr_pond_cesm,            &
                                   tr_pond,      tr_aero,     tr_FY,                   &
                                   tr_iage,      heat_capacity,  flag_old,             &
                                   flag_old_0
        real    (kind=dbl_kind) :: puny ,Tsfc ,rhos, Lfresh, jj, kk, njj
      
        ! Tracer dependencies and additional arrays
      
        integer (kind=int_kind), dimension(:),    allocatable    ::    &
                tracer_type    , & ! = 1, 2, or 3 (depends on 0, 1 or 2 other tracers)
                depend             ! tracer dependencies (see below)
      
        logical (kind=log_kind), dimension (:),   allocatable   ::     &
                has_dependents    ! true if a tracer has dependent tracers
      
        real (kind=dbl_kind),    dimension (:,:), allocatable ::       &
               works,works_old, vicen_init, aicen_init, vsnon_init
        real (kind=dbl_kind),    dimension (:), allocatable ::       &
               fiso_ocn,swild
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

        call icepack_query_parameters(heat_capacity_out=heat_capacity, ktherm_out=ktherm,    &
                                      puny_out=puny,rhos_out=rhos,                   &
                                      Lfresh_out=Lfresh)
        call icepack_query_tracer_sizes(ntrcr_out=ntrcr, nbtrcr_out=nbtrcr)
        call icepack_query_tracer_flags(                                  &
               tr_iage_out=tr_iage, tr_FY_out=tr_FY,                      &
               tr_aero_out=tr_aero, tr_pond_out=tr_pond,                  &
               tr_pond_cesm_out=tr_pond_cesm,                             &
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
        allocate (trcrn_ini(nx,ntrcr,ncat) )
        allocate (vicen_init(nx,ncat) )
        allocate (aicen_init(nx,ncat) )
        allocate (vsnon_init(nx,ncat) )
        works(:,:) = c0
        works_old(:,:) = c0
        fiso_ocn(:) = c0
        swild(:) = c0

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
        
        newwork(:) = c0

        call work_to_state (ntrcr, narr, works(:,:))

        do i=1,nx_nh
         flag_old_0 = .false.
         do n = 1, ncat   
            !flag_old = .false. 
                  do k = 1, nilyr
                     !if((trcrn_ini(i,nt_sice+k-1,n).eq.0).and.(trcrn_ini(i,nt_sice+k-1,n).ne.0)) flag_old_0= .true.
                     !if(abs(trcrn(i,nt_qice+k-1,n)/trcrn_ini(i,nt_qice+k-1,n))>1.1) flag_old_0= .true.
                     !if(abs(trcrn(i,nt_qice+k-1,n)/trcrn_ini(i,nt_qice+k-1,n))<0.9) flag_old_0= .true.
                     !if(abs(trcrn(i,nt_sice+k-1,n)/trcrn_ini(i,nt_sice+k-1,n))>1.1) flag_old_0= .true.
                     !if(abs(trcrn(i,nt_sice+k-1,n)/trcrn_ini(i,nt_sice+k-1,n))<0.9) flag_old_0= .true.
                     if(trcrn(i,nt_sice+k-1,n)>50) flag_old_0= .true.
                     if(trcrn(i,nt_sice+k-1,n)<0) flag_old_0= .true.
                  enddo
                  !do k = 1, nslyr
                  !   if(abs(trcrn(i,nt_qsno+k-1,n)/trcrn_ini(i,nt_qsno+k-1,n))<0.95) flag_old_0= .true.
                  !   if(abs(trcrn(i,nt_qsno+k-1,n)/trcrn_ini(i,nt_qsno+k-1,n))>1.05) flag_old_0= .true.
                  !enddo
                  if(flag_old_0) then
                     newwork(i) = 1
                    !write(*,*) i,nt,works(i,nt),works_old(i,nt),works(i,nt)/works_old(i,nt)
                 endif
         enddo

      end do
         do i=1,nx_nh
            if(newwork(i).eq.1)then
               do n=1,ncat
               !trcrn_ini(i,nt_sice:nt_sice+nilyr-1,n)
                  if (ktherm == 2) then
                     !vsnon(i,n) = c0
                     !trcrn(i,:,n) = trcrn_ini(i,:,n)
                     !trcrn(i,nt_Tsfc,n) = Tsfc
                     !do k = 1, nilyr
                     !  trcrn(i,nt_qice+k-1,n) = qin(k)
                        !if(any(trcrn_ini(i,nt_sice:nt_sice+nilyr-1,n).ne.0)) then
                           !trcrn(i,nt_sice:nt_sice+nilyr-1,n) = trcrn_ini(i,nt_sice:nt_sice+nilyr-1,n)
                        !else
                        !   trcrn(i,nt_sice:nt_sice+nilyr-1,n) = salinz(i,:)
                        !endif
                     !enddo
                     do k = 1, nilyr
                        njj=0
                        trcrn(i,nt_sice+k-1,n) = 0
                        do kk=1,nnp(i)
                           jj=indnd(kk,i)
                           if(trcrn_ini(jj,nt_sice+k-1,n)>0.and.trcrn_ini(jj,nt_sice+k-1,n)<50) then
                              trcrn(i,nt_sice+k-1,n) = trcrn(i,nt_sice+k-1,n) + trcrn_ini(jj,nt_sice+k-1,n)
                              njj=njj+1  
                           endif
                        enddo
                        if(njj.eq.0)then
                           trcrn(i,nt_sice+k-1,n) = 0
                        else
                           trcrn(i,nt_sice+k-1,n) = trcrn(i,nt_sice+k-1,n)/njj
                        endif
                     
                     enddo
         ! snow enthalpy
                     !do k = 1, nslyr
                     !  trcrn(i,nt_qsno+k-1,n) = -rhos * Lfresh
                     !enddo               ! nslyr
                  endif
               enddo
            endif
            do n=1,ncat
            !   if(vicen_init(i,n)<0.0001) then
                  !call icepack_init_trcr(T_air(i),    Tf(i),         &
                  !                   !salinz(i,:), Tmltz(i,:),    &
                  !                  trcrn(i,nt_sice:nt_sice+nilyr-1,n), Tmltz(i,:),    &
                  !                   Tsfc,                       &
                  !                   nilyr,        nslyr,        &
                  !                   qin   (:),    qsn  (:))

                  !call icepack_warnings_flush(ice_stderr)
                  !if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
                  !                   file=__FILE__, line=__LINE__)  

                  !trcrn(i,nt_Tsfc,n) = Tsfc
                  !trcrn(i,nt_qice:nt_qice+nilyr-1,n) = qin(:)
                  !trcrn(i,nt_qsno:nt_qsno+nslyr-1,n) = qsn(:)
                !  trcrn(i,:,n) = trcrn_ini(i,:,n)
              ! endif
               if (vicen(i,n)>10) write(12,*) 'after fct', i, n, vicen(i,n)
            enddo
        enddo
        if (ktherm == 2) then
            do n=1,ncat
               do k = 1, nilyr
                  swild=trcrn(:,nt_sice+k-1,n)
                  call exchange_p2d(swild)
                  trcrn(:,nt_sice+k-1,n)=swild
               enddo
            enddo
         endif
        ! cut off icepack
      
        call cut_off_icepack
       
        do i=1,nx
           !if (ncat < 0) then ! Do we really need this?
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
                                 heat_capacity,                                &
                                 first_ice(i,:),                               &
                                 trcr_depend(1:ntrcr),   trcr_base(1:ntrcr,:), &
                                 n_trcr_strata(1:ntrcr), nt_strata(1:ntrcr,:), &
                                 fpond(i),               fresh(i),             &
                                 fsalt(i),               fhocn(i),             &
                                 faero_ocn(i,:),         fiso_ocn,          &
                                 fzsal(i),               flux_bio(i,1:nbtrcr))
      
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
           !endif
            !flag_old_0 = .false.
            !do nt=1,ncat
            !   if( trcrn(i,nt_sice,nt)<0 .or. trcrn(i,nt_sice,nt)>50 ) flag_old_0= .true.
           ! enddo
            !if(flag_old_0) then
            !   write(*,*) i,trcrn(i,nt_sice,:),works(i,:),works_old(i,:),works(i,:)/works_old(i,:)           
            !endif
        end do

        deallocate(works)
        deallocate(works_old)
        deallocate(trcrn_ini)
        deallocate(vicen_init)
        deallocate(aicen_init)
        deallocate(vsnon_init)

    end subroutine tracer_advection_icepack

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
           tr_lvl, tr_pond_cesm, tr_pond_lvl, tr_pond_topo, heat_capacity ,tr_brine
  
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

        real (kind=dbl_kind) :: puny!,small
  
        real (kind=dbl_kind), parameter :: &
           small = 0.0000001_dbl_kind
         
        character(len=*), parameter :: subname = '(state_to_work)'
  
        call icepack_query_tracer_flags(tr_pond_cesm_out=tr_pond_cesm,              &
             tr_pond_lvl_out=tr_pond_lvl, tr_pond_topo_out=tr_pond_topo,            &
             tr_lvl_out=tr_lvl,tr_brine_out=tr_brine)
        call icepack_query_tracer_indices(nt_alvl_out=nt_alvl, nt_apnd_out=nt_apnd, &
             nt_fbri_out=nt_fbri, nt_qsno_out=nt_qsno, nt_vlvl_out=nt_vlvl,                              &
             nt_qice_out=nt_qice, nt_sice_out=nt_sice, nt_Tsfc_out=nt_Tsfc) 
        call icepack_query_parameters(rhoi_out=rhoi,     rhos_out=rhos,                   &
                                      Lfresh_out=Lfresh, heat_capacity_out=heat_capacity, &
                                      Tsmelt_out=Tsmelt, ktherm_out=ktherm,               &
                                      puny_out=puny)
        call icepack_warnings_flush(ice_stderr)
        if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
           file=__FILE__, line=__LINE__)
   
        trcrn(:,:,:) = c0
        aicen(:,:)   = c0
        vicen(:,:)   = c0
        vsnon(:,:)   = c0
        aice0(:)     = c0

        !small       = puny
        ! Open water fraction
  
        do i = 1, nx
           if (works(i,1) <= puny) then
               aice0(i) = c0
           else if (works(i,1) >= c1) then
               aice0(i) = c1
           else
               aice0(i) = works(i,1)
           end if
        enddo
        narrays = 1
  
        ! Sea ice area and volume per unit area of ice and snow
  
        do n=1,ncat
           do i = 1, nx
              if (works(i,narrays+1) > c1) then
                  works(i,narrays+1) = c1
              end if
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
  
        do i = 1, nx  ! For each grid cell
           if (sum(aicen(i,:)) > c1) then
              tmp(:) = c0
              exc(:) = c0
              do n = 1, ncat
                 if (aicen(i,n) > puny) tmp(n) = c1
              end do
              do n = 1, ncat
                  exc(n) = max(c0,(sum(aicen(i,:)) - c1))  &
                           * aicen(i,n) / sum(aicen(i,:))
              end do
              do n = 1, ncat
                 aicen(i,n) = max(c0,aicen(i,n) - exc(n))
                 aice0      = max(c0,1-sum(aicen(i,:)))
              end do
           end if
        end do
  
        narrays = 1

        do n=1, ncat
            !do i=1, nx
           !    aicen(i,n) = works(i,narrays+1)
           !    vicen(i,n) = works(i,narrays+2)
           !    vsnon(i,n) = works(i,narrays+3)
           ! enddo
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
           enddo
          
           narrays = narrays + ntrcr
        enddo

        do i = 1, nx
         aice(i) = c0
         vice(i) = c0
         vsno(i) = c0
         do it = 1, ntrcr
            trcr(i,it) = c0
         enddo
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
         end do


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
           nt_alvl, nt_apnd, nt_fbri, nt_Tsfc
  
        logical (kind=log_kind) :: &
           tr_pond_cesm, tr_pond_lvl, tr_pond_topo
  
        integer (kind=int_kind) ::      &
           i, n, it   , & ! counting indices
           narrays    , & ! counter for number of state variable arrays
           nt_qsno
  
        real (kind=dbl_kind) :: &
           rhos       , &
           Lfresh
  
        character(len=*), parameter :: subname = '(state_to_work)'
  
        call icepack_query_tracer_flags(tr_pond_cesm_out=tr_pond_cesm, &
             tr_pond_lvl_out=tr_pond_lvl, tr_pond_topo_out=tr_pond_topo)
        call icepack_query_tracer_indices(nt_alvl_out=nt_alvl, nt_apnd_out=nt_apnd, &
             nt_fbri_out=nt_fbri, nt_qsno_out=nt_qsno, nt_Tsfc_out=nt_Tsfc)
        call icepack_query_parameters(rhos_out=rhos, Lfresh_out=Lfresh)
        call icepack_warnings_flush(ice_stderr)
        if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
           file=__FILE__, line=__LINE__)
  
        do i = 1, nx
           works(i,1) = aice0(i)
        enddo
        narrays = 1
  
        do n=1, ncat
  
           do i = 1, nx
              works(i,narrays+1) = aicen(i,n)
              works(i,narrays+2) = vicen(i,n)
              works(i,narrays+3) = vsnon(i,n)

              if (vicen(i,n)>10) write(12,*) 'before fct', i, n, vicen(i,n)
           enddo                  ! i
           narrays = narrays + 3
  
           do it = 1, ntrcr
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
                        works(i,narrays+it) = vsnon(i,n)*trcrn(i,it,n) 
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
                      tr_pond_cesm .or. tr_pond_topo) then
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

    module subroutine cut_off_icepack

        use icepack_intfc,         only: icepack_compute_tracers
        use icepack_intfc,         only: icepack_aggregate
        use icepack_intfc,         only: icepack_init_trcr
        use icepack_intfc,         only: icepack_sea_freezing_temperature
        use icepack_therm_shared,  only: calculate_Tin_from_qin
        use icepack_mushy_physics, only: icepack_mushy_temperature_mush
        use icepack_mushy_physics, only: liquidus_temperature_mush
        use icepack_mushy_physics, only: enthalpy_mush
  
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
                                   tr_pond_cesm, tr_pond_topo, tr_pond_lvl, tr_FY, tr_iage,   &
                                   heat_capacity
        integer (kind=int_kind) :: nt_Tsfc, nt_qice, nt_qsno, nt_sice
        integer (kind=int_kind) :: nt_fbri, nt_alvl, nt_vlvl, nt_apnd, nt_hpnd, nt_ipnd, nt_FY, nt_iage
  
        character(len=*), parameter :: subname = '(cut_off_icepack)'
  
        call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
        call icepack_query_tracer_flags(tr_brine_out=tr_brine, tr_lvl_out=tr_lvl,      &
          tr_pond_cesm_out=tr_pond_cesm, tr_pond_topo_out=tr_pond_topo,                &
          tr_pond_lvl_out=tr_pond_lvl, tr_FY_out=tr_FY, tr_iage_out=tr_iage            )
        call icepack_query_tracer_indices( nt_Tsfc_out=nt_Tsfc, nt_qice_out=nt_qice,   &
          nt_qsno_out=nt_qsno, nt_sice_out=nt_sice, nt_FY_out=nt_FY,                   &
          nt_fbri_out=nt_fbri, nt_alvl_out=nt_alvl, nt_vlvl_out=nt_vlvl,               &         
          nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, nt_ipnd_out=nt_ipnd,               &
          nt_iage_out=nt_iage                                                          )
        call icepack_query_parameters(rhos_out=rhos, rhoi_out=rhoi, Lfresh_out=Lfresh, &
                                      cp_ice_out=cp_ice, cp_ocn_out=cp_ocn)
        call icepack_query_parameters(depressT_out=depressT, puny_out=puny, &
          Tsmelt_out=Tsmelt, ktherm_out=ktherm, heat_capacity_out=heat_capacity)
        call icepack_warnings_flush(ice_stderr)
  
        small = puny
        !small=0.005_dbl_kind        
        Tmin  = -100.0_dbl_kind
  
        if (.not. heat_capacity) then ! for 0 layer thermodynamics
           do n = 1, ncat
              do i = 1, nx
                 if (trcrn(i,nt_Tsfc,n) > Tf(i) .or. trcrn(i,nt_Tsfc,n)< Tmin) then
                     trcrn(i,nt_Tsfc,n) = min(Tf(i), (T_air(i) + 273.15_dbl_kind))
                 endif
              enddo
           enddo
        endif
  
        if (heat_capacity) then    ! only for bl99 and mushy thermodynamics
  
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
                        qin_max(k) = enthalpy_mush((tmltz0(k)-0.1), trcrn(i,nt_sice+k-1,n))
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
                     ! write(12,*) i,zTin(k),qin_max(k), qin(k), trcrn(i,nt_qice+k-1,n)
                      trcrn(i,nt_Tsfc,n) = Tsfc

                      do k = 1, nilyr
                        !write(12,*) i,zTin(k),qin_max(k), qin(k), trcrn(i,nt_qice+k-1,n),tmltz0(k), trcrn(i,nt_sice+k-1,n)
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

                  end if
  
                      if (vsnon(i,n) > small) then
  
                          flag_snow = .false.
  
                          do k = 1, nslyr  
                             if (trcrn(i,nt_qsno+k-1,n) >= -rhos*Lfresh) flag_snow = .true.
                             zTsn(k) = (Lfresh + trcrn(i,nt_qsno+k-1,n)/rhos)/cp_ice
                             if (zTsn(k) < Tmin) flag_snow = .true.
                          end do 
  
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
                 ! Sea ice age
                 if (tr_iage) then
                     if (trcrn(i,nt_iage,n) < 0.000001_dbl_kind) trcrn(i,nt_iage,n) = c0
                 end if
                 ! First year ice fraction
                 if (tr_FY) then
                     if (trcrn(i,nt_FY,n) < 0.000001_dbl_kind) trcrn(i,nt_FY,n) = c0
                     if (trcrn(i,nt_FY,n) > c1) trcrn(i,nt_FY,n) = c1
                 end if
                 ! Level ice
                 if (tr_lvl) then
                     if (trcrn(i,nt_alvl,n) > aicen(i,n)) then
                         trcrn(i,nt_alvl,n) = aicen(i,n)
                     elseif (trcrn(i,nt_alvl,n) < 0.000001_dbl_kind) then
                         trcrn(i,nt_alvl,n) = c0
                     endif
                     if (trcrn(i,nt_vlvl,n) < 0.000001_dbl_kind .or. trcrn(i,nt_alvl,n) < 0.000001_dbl_kind) trcrn(i,nt_vlvl,n) = c0
                     if (trcrn(i,nt_vlvl,n) > vicen(i,n)) trcrn(i,nt_vlvl,n) = vicen(i,n)
                 end if
                 ! CESM melt pond parameterization
                 if (tr_pond_cesm) then
                     if (trcrn(i,nt_apnd,n) > c1) then
                         trcrn(i,nt_apnd,n) = c1
                     elseif (trcrn(i,nt_apnd,n) < 0.000001_dbl_kind) then
                         trcrn(i,nt_apnd,n) = c0
                     endif
                     if (trcrn(i,nt_hpnd,n) < 0.000001_dbl_kind .or. trcrn(i,nt_apnd,n) < 0.000001_dbl_kind) trcrn(i,nt_hpnd,n) = c0
                     if (trcrn(i,nt_hpnd,n) > hpnd_max) trcrn(i,nt_hpnd,n) = hpnd_max
                 end if
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

        do i = 1, nx
           aice(i) = c0
           vice(i) = c0
           vsno(i) = c0
           do it = 1, ntrcr
              trcr(i,it) = c0
           enddo
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
        end do
  
        end if ! heat_capacity
  
        call icepack_warnings_flush(ice_stderr)
        if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
            file=__FILE__, line=__LINE__)

    end subroutine cut_off_icepack

end submodule icedrv_advection


