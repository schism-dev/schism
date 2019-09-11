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

!===============================================================================
!===============================================================================
! SCHISM HYDRAULIC SUBROUTINES

      module hydraulic_structures
      use schism_glbl, only : rkind,grav,in_dir,out_dir,len_in_dir,len_out_dir
      use schism_msgp, only : parallel_abort
      implicit none


      integer, parameter :: NAME_LEN = 32


!---- Constants for struct_type
      integer, parameter :: HYDTRANSFER = 1  !< simple prescribed transfer of flow from A to B
      integer, parameter :: ORIFICE     = 2  !< rectangular orifice
      integer, parameter :: WEIR        = 3  !< weir flow, including transition to submerged flow
      integer, parameter :: PIPE        = 4  !< pipe flow -- right now just circular orifice
      integer, parameter :: RADIAL      = 5  !< radial gate -- simplified form of HEC-RAS formula with default exponents
      integer, parameter :: RADIAL_RH   = 6  !< radial gate -- simplified form of HEC-RAS formula with default exponents
      integer, parameter :: WEIR_PIPE   = 7  !< compound structure with culvert flow through an overtoppable barrier

!---- Constants for ts_unitradia
      integer, parameter :: TS_DISABLED       = -1

!---- Type representing a structure
      type hydraulic_structure
      real(rkind) :: elev   = HUGE(1.d0)    !< invert or crest or bottom elevation wrt datum
      real(rkind) :: width  = HUGE(1.d0)    !< width or pipe radius
      real(rkind) :: height = HUGE(1.d0)    !< height or vertical length of aperture, should be > 2*radius (ie, width) for circular
      real(rkind) :: coef   = HUGE(1.d0)    !< flow coefficient (physical)
      real(rkind) :: coef_lin = HUGE(1.d0)  !< flow coefficient part that varies linearly (eg with submergence, relative gate hght)
      real(rkind) :: op_down = 1.0d0        !< operation (modulation from 0 to 1) of the gate flow in downstream direction
      real(rkind) :: op_up   = 1.d0         !< operation (modulation from 0 to 1) of the gate flow in the reverse direction
      real(rkind) :: prescribed_flow =0.0d0 !< prescribed flow (relevant for transfers)
      real(rkind) :: next_read_time         !< next time to read data for the structure              
      real(rkind),   dimension(3)  :: normal_dir = (/0.d0, 0.d0, 0.d0/)
      character(LEN=NAME_LEN)      :: struct_name = ' '     !< text name of structure
      logical     :: is_local = .FALSE.     !< structure lies on local processor
      integer     :: ts_unit = TS_DISABLED  !< unit number for time series input or TS_UNIT_NOT_INIT or TS_DISABLED     
      
      integer(rkind), allocatable  :: node_pairs(:,:)       !< corresponding (global) node pairs on either side of structure
      integer                      :: struct_type = 0       !< type of structure (pipe, weir, etc). See constants above for acceptable values

      integer     :: nduplicate = 1         !< number of identical smaller structures treated as one device. inputs are for one struct.
      integer     :: upnode                 !< upstream reference node
      integer     :: downnode               !< downstream reference node
      integer     :: npair
      logical     :: install = .true.      !< at the moment this is redundant with an array in the main code
      type(hydraulic_structure), pointer  :: substruct => null()   !< pointer to a subordinate structure. Used for compound type like WEIR_PIPE
      end type

      type (hydraulic_structure), allocatable,dimension(:), target, save :: structures
      real(rkind)   :: block_nudge
      integer       :: nhtblocks = 0            ! Number of structures


      contains  
      

!>===============================================================================
!>     Parse and process structure information
!>===============================================================================
      subroutine load_structures(filename)
      use schism_glbl, only : errmsg
      implicit none
      character(LEN=*) :: filename
      integer          :: iblock,in,npair,ipair
      integer          :: nodeup,nodedown
      integer          :: install_by_default, itimeseries
      character(LEN=NAME_LEN) :: char_struct_type
      character(LEN=NAME_LEN) :: struct_name

!...  Read in hydraulics.in,ipair
      open(31,file=in_dir(1:len_in_dir)//filename,status='old')
      
      !Specify blocks for hydraulic transfer structures (where fluxes are specified,
      !and tracers are conserved)
      read(31,*) nhtblocks !# of transfer blocks
      read(31,*) block_nudge !nudging factor
      allocate(structures(nhtblocks))
      if(block_nudge < 0.d0 .or.block_nudge>=1.d0) call parallel_abort('MAIN: wrong block_nudge')

      if(nhtblocks>0) then
        do in = 1,nhtblocks
          read(31,*)iblock,struct_name
          structures(iblock)%struct_name = trim(struct_name)
          read(31,*)npair,nodeup,nodedown
          if(npair<2) call parallel_abort('load_structures: Too few pairs of nodes in hydraulic structure/transfer')
          structures(iblock)%npair = npair
          structures(iblock)%upnode = nodeup
          structures(iblock)%downnode = nodedown
          allocate(structures(iblock)%node_pairs(2,npair))
          do ipair = 1,npair
            read(31,*)structures(iblock)%node_pairs(1,ipair),structures(iblock)%node_pairs(2,ipair)
          end do
          read(31,*)char_struct_type
          read(31,*)structures(iblock)%nduplicate 
          select case(trim(char_struct_type))
          case("transfer")
            structures(iblock)%struct_type = HYDTRANSFER
            read(31,*)structures(iblock)%prescribed_flow
          case("orifice")
            structures(iblock)%struct_type = ORIFICE
            read(31,*)structures(iblock)%elev, structures(iblock)%width, structures(iblock)%height
            read(31,*)structures(iblock)%coef, structures(iblock)%op_down, structures(iblock)%op_up
          case("radial")
            structures(iblock)%struct_type = RADIAL
            read(31,*)structures(iblock)%elev, structures(iblock)%width, structures(iblock)%height
            read(31,*)structures(iblock)%coef, structures(iblock)%op_down, structures(iblock)%op_up
          case("radial_relheight")
            structures(iblock)%struct_type = RADIAL_RH
            read(31,*)structures(iblock)%elev, structures(iblock)%width, structures(iblock)%height
            read(31,*)structures(iblock)%coef,structures(iblock)%coef_lin
            read(31,*)structures(iblock)%op_down, structures(iblock)%op_up
          case("weir")
            structures(iblock)%struct_type = WEIR
            read(31,*)structures(iblock)%elev, structures(iblock)%width
            read(31,*)structures(iblock)%coef, structures(iblock)%op_down, structures(iblock)%op_up
          case("culvert")
            structures(iblock)%struct_type = PIPE
            read(31,*)structures(iblock)%elev, structures(iblock)%width
            read(31,*)structures(iblock)%coef, structures(iblock)%op_down, structures(iblock)%op_up
          case("weir_culvert")
            structures(iblock)%struct_type = WEIR_PIPE
            allocate(structures(iblock)%substruct)
            ! read in weir params
            read(31,*)structures(iblock)%elev, structures(iblock)%width
            read(31,*)structures(iblock)%coef, structures(iblock)%op_down, structures(iblock)%op_up
            ! read in pipe params
            read(31,*)structures(iblock)%substruct%nduplicate
            read(31,*)structures(iblock)%substruct%elev, structures(iblock)%substruct%width
            read(31,*)structures(iblock)%substruct%coef, structures(iblock)%substruct%op_down, structures(iblock)%substruct%op_up
          case default
	    write(errmsg,*)"Structure type code not recognized for hydraulic structure ",structures(iblock)%struct_name
            call parallel_abort(errmsg)
          end select
          read(31,*) itimeseries
          if (itimeseries /= 0) then
            structures(iblock)%ts_unit = TS_DISABLED + 1
            structures(iblock)%install = .FALSE.
          else 
            structures(iblock)%ts_unit = TS_DISABLED
          end if
        end do
      end if
      return
      end subroutine load_structures


!>===============================================================================
!>    Initialize I/O units for hydraulic structures. The expected file name follows
!>    the pattern structure_name.th
!>===============================================================================     
      subroutine init_struct_time_series(time)
      use schism_glbl, only: rkind
      implicit none
      real(rkind) :: time
      real(rkind) :: time_next
      integer :: istruct  ! local counter
      integer :: fort_unit
      integer, parameter :: struct_base_unit = 3000

      do istruct = 1,nhtblocks
        if (structures(istruct)%ts_unit /= TS_DISABLED) then
          if (structures(istruct)%is_local ) then
            fort_unit = struct_base_unit + istruct
            structures(istruct)%ts_unit = fort_unit
            open(fort_unit,file=in_dir(1:len_in_dir)//trim(structures(istruct)%struct_name) // '.th',status='old')
            call irreg_time_history_advance(fort_unit, time, time_next)
            structures(istruct)%next_read_time = time
          end if
        end if
      end do
      return
      end subroutine init_struct_time_series



!>===============================================================================
!>    Read parameters for hydraulic structures. 
!>=============================================================================== 

      subroutine read_struct_ts(time)
      use schism_glbl, only: errmsg,rkind
      implicit none
      real(rkind), intent(in) :: time
      real(rkind)             :: time_next
      real(rkind)             :: ttt
      integer istruct
      integer fort_unit
      integer install
      do istruct = 1,nhtblocks
        fort_unit = structures(istruct)%ts_unit
        if ((fort_unit /= TS_DISABLED) .and.structures(istruct)%is_local ) then          
          if (time >= structures(istruct)%next_read_time)then
            ! forward to just at/after current time and read in new parameters
            call irreg_time_history_advance(fort_unit,time,time_next)
            structures(istruct)%next_read_time = time_next     ! pre-fetching the new read time
            select case(structures(istruct)%struct_type)
            case(HYDTRANSFER)
              read(fort_unit,*)ttt,install, structures(istruct)%prescribed_flow
            case(PIPE, WEIR)
              read(fort_unit,*)ttt,install,&
                               structures(istruct)%nduplicate,structures(istruct)%op_down, structures(istruct)%op_up, &
                               structures(istruct)%elev,structures(istruct)%width
            case(RADIAL,RADIAL_RH,ORIFICE)
              read(fort_unit,*)ttt,install,&
                             structures(istruct)%nduplicate, structures(istruct)%op_down, structures(istruct)%op_up, &
                             structures(istruct)%elev,structures(istruct)%width, structures(istruct)%height
            case(WEIR_PIPE)
              read(fort_unit,*)ttt,install,&
                               structures(istruct)%nduplicate,structures(istruct)%op_down, structures(istruct)%op_up, &
                               structures(istruct)%elev,structures(istruct)%width, &
                               structures(istruct)%substruct%nduplicate,structures(istruct)%substruct%op_down, structures(istruct)%substruct%op_up, &
                               structures(istruct)%substruct%elev,structures(istruct)%substruct%width

            case default
	      write(errmsg,*)"Structure type code not recognized for hydraulic structure ",structures(istruct)%struct_name
              call parallel_abort(errmsg)
            end select
            structures(istruct)%install = .not. (install==0)
          end if
        end if
      end do  
      return
      end subroutine



!>===============================================================================
!>     Compute thru flow at a hydraulic structure
!>===============================================================================
      subroutine calc_struc_flow(istruct,elev_up,elev_down,flow)
      use schism_glbl, only : rkind,grav,pi,errmsg
      use schism_msgp, only : parallel_abort
      implicit none
      integer, intent(in)     :: istruct !index of the struc
      real(rkind), intent(in)  :: elev_up            ! elevation at reference points '1' ('upstream')    [m]
      real(rkind), intent(in)  :: elev_down          ! elevation at reference points  '2' ('downstream') [m]
      real(rkind), intent(out) :: flow                 ! flow across the struc [m^3/s]
      
! --- locals
      real(rkind), parameter   :: sqrt2g = 4.428690551d0
      real(rkind), parameter   :: PART_SUBMERGE = 2.d0/3.d0
      real(rkind), parameter   :: FULL_SUBMERGE = 0.8d0
      real(rkind), parameter   :: expon = 1.5d0

      real(rkind) :: diff
      real(rkind) :: op
      real(rkind) :: signed_coef,depth_flow
      real(rkind) :: area,radius,angle
      real(rkind) :: subflow,subfrac,submerge_ratio,relheight
      real(rkind) :: max_elev, min_elev, dratio, attenuation
      real(rkind) :: tailw_depth,headw_energy_head
      real(rkind) :: coef_matching_factor
      logical :: main_struct_dry 
      type(hydraulic_structure), pointer :: substruct

!todo: consider whether orifice and weir need different code ... unsubmerged weir-like transition has to be handled by orifice
!todo: depth_flow being floored to struct.height is weird for weirs, and ruins the attenuation (can create NaNs)
!todo: very rudimentary pipes -- basically orifice with correct area of flow in structure for circular pipe
!todo: for pipes (really everything) might consider tabulated ratings like those produced by the Culvert Analysis Program (USGS)
!todo: upstream energy head is referenced in a few places, but in all cases this is currently just the surface, no velocity
!todo: note interpretation of the flow coefficient for radial gates, which has to provide smooth transitions at submergence 
!      coef_matching_factor" that ensures this transition is kludgy
!todo: consider energy-momentum method for radial gates


      type(hydraulic_structure), pointer :: struct
      struct => structures(istruct)

      signed_coef  = sign(struct%coef,elev_up-elev_down)
      max_elev = max(elev_up,elev_down)
      min_elev = min(elev_up,elev_down)
      

      if (signed_coef > 0.d0) then
        ! flow is up to down
         op = struct%op_down
       else
        ! down to up
          op = struct%op_up
      end if

      if(struct%struct_type .ne. HYDTRANSFER) then
        depth_flow = min(struct%height,max_elev - struct%elev)
        if (depth_flow <= 0.D0 .or. op == 0.d0)then
          flow = 0.d0
          main_struct_dry = .true.                      ! The substruct could still be wet
	  if (.not. associated(struct%substruct)) return  ! If there is no substruct bail
        else
          main_struct_dry = .false.
        end if
      end if

      select case(struct%struct_type)
      case(HYDTRANSFER)
        flow = struct%prescribed_flow
        return
      case(ORIFICE)
!       rectangular only
!       acts like orifice or weir
!       submerged weir case not handled
        depth_flow = min(struct%height, depth_flow)
        area = depth_flow*struct%width
        diff = min(max_elev - struct%elev, max_elev - min_elev)
        flow = signed_coef*area*sqrt2g*sqrt(diff)       ! m² * sqrt(m)/s * sqrt(m)  = m³/s
        flow = flow*dble(struct%nduplicate)*op
      case(WEIR)
        area = struct%width*depth_flow
        flow = signed_coef*area*sqrt2g*sqrt(depth_flow)
        if (min_elev > struct%elev) then
          dratio = (min_elev - struct%elev)/depth_flow
          attenuation = (1.d0 - (dratio)**expon)**3.85d-1
          flow = flow*attenuation
        end if 
        flow = flow*dble(struct%nduplicate)*op
      case(PIPE)
        ! todo: basically orifice with no treatment of partially submerged case with tailwater > invert
        ! todo: need an assert that guarantees height is big when not explicitly set
        radius = struct%width
        if (depth_flow < 2.d0*radius) then
!         partial flow
          angle=acos(1.-depth_flow/radius) 
          area = radius**2*angle-radius*(radius-depth_flow)*sin(angle)
        else
!         full flow
          area = radius*radius*pi
        end if
        diff = min(max_elev - struct%elev, max_elev - min_elev)
        flow = signed_coef*area*sqrt2g*sqrt(diff)        
        flow = flow*dble(struct%nduplicate)*op
      case(WEIR_PIPE)
        ! weir part of flow is the "main" structure. skip if it is dry 
        ! todo: test that this is the same as WEIR for pipe with zero radius
        if (main_struct_dry)then
            flow = 0.d0
        else
            area = struct%width*depth_flow
            flow = signed_coef*area*sqrt2g*sqrt(depth_flow)
            if (min_elev > struct%elev) then
              dratio = (min_elev - struct%elev)/depth_flow
              attenuation = (1.d0 - (dratio)**expon)**3.85d-1
              flow = flow*attenuation
            end if 
        end if
        flow = flow*dble(struct%nduplicate)*op

        ! pipe data is contained in struct%substruct
        substruct => struct%substruct
        signed_coef  = sign(substruct%coef,elev_up-elev_down)
        if (signed_coef > 0.d0) then
          ! flow is up to down
          op = substruct%op_down
        else
          ! down to up
          op = substruct%op_up
        end if
        depth_flow = min(substruct%height,max_elev - substruct%elev)
        if (depth_flow <= 0.D0 .or. op == 0.d0)then
          flow = flow + 0.d0   ! no contribution from substruct
	  return
        end if

        radius = substruct%width
        if (depth_flow < 2.d0*radius) then
!         partial flow
          angle=acos(1.-depth_flow/radius) 
          area = radius**2*angle-radius*(radius-depth_flow)*sin(angle)
        else
!         full flow
          area = radius*radius*pi
        end if
        diff = min(max_elev - substruct%elev, max_elev - min_elev)
        flow = flow + signed_coef*area*sqrt2g*sqrt(diff)*dble(substruct%nduplicate)*op

      case(RADIAL)
        headw_energy_head = max_elev - struct%elev
        tailw_depth = min_elev - struct%elev
        submerge_ratio = tailw_depth/headw_energy_head
        area = struct%width*depth_flow        
        if (submerge_ratio < PART_SUBMERGE) then
	    ! free flow
            flow = signed_coef*sqrt2g*area*sqrt(headw_energy_head)
            flow = flow*dble(struct%nduplicate)*op
        else if (submerge_ratio < FULL_SUBMERGE) then
            ! partially submerged (transitional) flow
            ! matches free case at PART_SUBMERGE and transitions linearly to submerged
            diff = max_elev - min_elev
            coef_matching_factor = sqrt(1.d0/(1.d0-PART_SUBMERGE))
            flow = signed_coef*area*sqrt2g*sqrt(diff)
            ! now weigh the two so that the flow makes a linear transition
            subfrac = (submerge_ratio-PART_SUBMERGE)/(FULL_SUBMERGE - PART_SUBMERGE)
            flow = ((1.d0 - subfrac)*coef_matching_factor + subfrac)*flow
            flow = flow*dble(struct%nduplicate)*op
        else 
            ! fully submerged flow uses orifice equational
            diff = max_elev-min_elev
            flow = signed_coef*area*sqrt2g*sqrt(diff)
            flow = flow*dble(struct%nduplicate)*op
        end if

      case(RADIAL_RH)
        headw_energy_head = max_elev - struct%elev
        tailw_depth = min_elev - struct%elev
        relheight = min(1.d0,struct%height/headw_energy_head)
        area = struct%width*depth_flow
        signed_coef = signed_coef+sign(relheight*struct%coef_lin,elev_up-elev_down)
        diff = max_elev-min_elev
        flow = signed_coef*area*sqrt2g*sqrt(diff)
        flow = flow*dble(struct%nduplicate)*op
      case default
	write(errmsg,*)"Structure type code not recognized for hydraulic structure #",struct%struct_name
        call parallel_abort(errmsg)
      end select

!-----Adjust due to gate scheduling and number of duplicate devices

      return
      end subroutine calc_struc_flow

!=============== Hardwired initialization example

      subroutine init_hydraulic_structures(nstruct)
      implicit none
      integer,intent(in) :: nstruct
      integer            :: istruct
      type(hydraulic_structure),pointer:: struct       
      nhtblocks = nstruct
      allocate(structures(nstruct))


      ! you would put in real numbers for these - JZ: I assume you'll continue on this part
      do istruct = 1,nstruct
        struct => structures(istruct)
        struct%struct_type = istruct
        struct%elev = -1.d0
        if (istruct == PIPE) then
          struct%width = 2.0d0
          struct%height = 2.d0*struct%width
        end if
	struct%width = 10.d0
        struct%height = 4.0d0
        struct%op_down = 0.5d0
        struct%op_up = 1.0d0
        struct%coef = 0.75d0
	struct%prescribed_flow = 5.d0
      end do
      
      return
      end subroutine init_hydraulic_structures

      subroutine finalize_hydraulic_structures
      implicit none
      integer :: istruct
      if (nhtblocks > 0)then

        if (allocated(structures))then
            ! first deallocate substructures
            do istruct =1,nhtblocks
                if (associated(structures(istruct)%substruct))then
                    deallocate(structures(istruct)%substruct)
                end if
            end do
          deallocate(structures)
        end if
      end if
      return
      end subroutine finalize_hydraulic_structures

!>===============================================================================
!>    Advance a time history file to the record just at (=) or after (>) the input time
!>    The routine also pre-fetches the timing of the next record. 
!>    todo: EOF issues are not currently handled
!>===============================================================================
      subroutine irreg_time_history_advance(unit, time, time_next)
      use schism_glbl, only : rkind
      implicit none

      integer, intent(in)      :: unit       !< Fortran unit to be read
      real(rkind), intent(in)  :: time       !< Time to which to advance the file
      real(rkind), intent(out) :: time_next  !< Pre-fetched time at which to expect the next value
      real(rkind) :: t  ! dummy for reading in case reading fails
      real(rkind) :: time_read 
      time_read = 0.d0

      do while(time_read <= time)
        read(unit,*,end=10)t
        time_read = t
      end do
      time_next = time_read
      backspace(unit)
      backspace(unit)
      return
 10   continue
      ! Special code for when the end of file is reached.
      time_next = huge(1.d0)
      backspace(unit)
      backspace(unit)
      return
      end subroutine


      end module hydraulic_structures
