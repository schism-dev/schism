module spherepack_interfaces

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: get_wavetable_size, init_wavetable
    public :: scalar_analysis, scalar_synthesis
    public :: vector_analysis, vector_synthesis

    abstract interface
        pure function get_wavetable_size(nlat, nlon) &
            result(return_value)
            import :: ip

            ! Dummy arguments
            integer(ip), intent(in) :: nlat
            integer(ip), intent(in) :: nlon
            integer(ip)             :: return_value
        end function get_wavetable_size

        subroutine init_wavetable(nlat, nlon, wavetable, error_flag)
            import :: ip, wp

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wavetable(:)
            integer(ip), intent(out) :: error_flag
        end subroutine init_wavetable

        subroutine scalar_analysis(nlat, nlon, isym, nt, g, idg, jdg, a, b, &
            mdab, ndab, wavetable, error_flag)
            import :: ip, wp

            ! Dummy arguments
            integer(ip), intent(in)   :: nlat
            integer(ip), intent(in)   :: nlon
            integer(ip), intent(in)   :: isym
            integer(ip), intent(in)   :: nt
            real(wp),    intent(in)   :: g(idg, jdg, nt)
            integer(ip), intent(in)   :: idg
            integer(ip), intent(in)   :: jdg
            real(wp),    intent(out)  :: a(mdab, ndab, nt)
            real(wp),    intent(out)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)   :: mdab
            integer(ip), intent(in)   :: ndab
            real(wp),    intent(in)   :: wavetable(:)
            integer(ip), intent(out)  :: error_flag
        end subroutine scalar_analysis

        subroutine scalar_synthesis(nlat, nlon, isym, nt, g, idg, jdg, a, b, &
            mdab, ndab, wavetable, error_flag)
            import :: ip, wp

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: g(idg, jdg, nt)
            integer(ip), intent(in)  :: idg
            integer(ip), intent(in)  :: jdg
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wavetable(:)
            integer(ip), intent(out) :: error_flag
        end subroutine scalar_synthesis

        subroutine vector_analysis(nlat, nlon, ityp, nt, v, w, idvw, jdvw, &
            br, bi, cr, ci, mdab, ndab, wavetable, error_flag)
            import :: ip, wp

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(in)  :: v(idvw, jdvw, nt)
            real(wp),    intent(in)  :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(out) :: br(mdab, ndab, nt)
            real(wp),    intent(out) :: bi(mdab, ndab, nt)
            real(wp),    intent(out) :: cr(mdab, ndab, nt)
            real(wp),    intent(out) :: ci(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wavetable(:)
            integer(ip), intent(out) :: error_flag
        end subroutine vector_analysis

        subroutine vector_synthesis(nlat, nlon, ityp, nt, v, w, idvw, jdvw, &
            br, bi, cr, ci, mdab, ndab, wavetable, error_flag)
            import :: ip, wp

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: v(idvw, jdvw, nt)
            real(wp),    intent(out) :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: br(mdab, ndab, nt)
            real(wp),    intent(in)  :: bi(mdab, ndab, nt)
            real(wp),    intent(in)  :: cr(mdab, ndab, nt)
            real(wp),    intent(in)  :: ci(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wavetable(:)
            integer(ip), intent(out) :: error_flag
        end subroutine vector_synthesis
    end interface

end module spherepack_interfaces
