! Bilinear interpolation of calculated NWM fluxes to SCHISM element centers.
module hydrologic_exchange_mapping
  use schism_glbl, only: rkind
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none
  private

  type, public :: soil_mapping
    integer, allocatable :: cell(:,:)
    real(rkind), allocatable :: weight(:,:)
  end type

  public :: initialize_soil_mapping, map_soil_fluxes

contains

  subroutine initialize_soil_mapping(mapping, x, y, element_x, element_y)
    ! Coordinates are in the NWM projection, prepared with nwm_io.py's
    ! lonlat_to_ldas_xy. The ascending NWM coordinate vectors index cell centers.
    type(soil_mapping), intent(out) :: mapping
    real(rkind), intent(in) :: x(:), y(:), element_x(:), element_y(:)
    integer :: i, ix, iy, nx
    real(rkind) :: rx, ry

    nx=size(x)
    allocate(mapping%cell(4,size(element_x)),mapping%weight(4,size(element_x)))
    mapping%cell=0; mapping%weight=0._rkind
    do i=1,size(element_x)
      ix=lower_cell(x,element_x(i))
      iy=lower_cell(y,element_y(i))
      if(ix==0.or.iy==0) cycle
      rx=(element_x(i)-x(ix))/(x(ix+1)-x(ix))
      ry=(element_y(i)-y(iy))/(y(iy+1)-y(iy))
      mapping%cell(:,i)=(/ix+(iy-1)*nx,ix+1+(iy-1)*nx,ix+iy*nx,ix+1+iy*nx/)
      mapping%weight(:,i)=(/(1-rx)*(1-ry),rx*(1-ry),(1-rx)*ry,rx*ry/)
    enddo
  end subroutine initialize_soil_mapping


  subroutine map_soil_fluxes(mapping, q_surf, q_perc, element_surf, element_perc)
    ! Outside the input subset and missing input data propagate NaN, as in
    ! xarray linear interpolation. Area conversion is done after interpolation.
    type(soil_mapping), intent(in) :: mapping
    real(rkind), intent(in) :: q_surf(:), q_perc(:)
    real(rkind), intent(out) :: element_surf(:), element_perc(:)
    integer :: i

    element_surf=ieee_value(0._rkind,ieee_quiet_nan)
    element_perc=element_surf
    do i=1,size(element_surf)
      if(mapping%cell(1,i)==0) cycle
      element_surf(i)=sum(mapping%weight(:,i)*q_surf(mapping%cell(:,i)))
      element_perc(i)=sum(mapping%weight(:,i)*q_perc(mapping%cell(:,i)))
    enddo
  end subroutine map_soil_fluxes


  integer function lower_cell(x, value) result(left)
    real(rkind), intent(in) :: x(:), value
    integer :: right, middle

    left=0
    if(value<x(1).or.value>x(size(x))) return
    left=1; right=size(x)
    do while(right-left>1)
      middle=(left+right)/2
      if(value>=x(middle)) then
        left=middle
      else
        right=middle
      endif
    enddo
  end function lower_cell

end module hydrologic_exchange_mapping
