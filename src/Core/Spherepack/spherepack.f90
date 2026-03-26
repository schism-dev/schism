module spherepack

    use, intrinsic :: ISO_Fortran_env, only: &
        stderr => ERROR_UNIT, &
        stdout => OUTPUT_UNIT

    use spherepack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        HALF_PI, &
        PI, &
        TWO_PI

    use divergence_routines, only: &
        divec, dives, divgc, divgs, &
        idivec, idives, idivgc, idivgs

    use gaussian_latitudes_and_weights_routines, only: &
        compute_gaussian_latitudes_and_weights

    use coordinate_transfer_routines, only: &
        geo2maths, math2geos, geo2mathv, math2geov

    use grid_transfer_routines, only: &
        trssph, sshifte, sshifti, initialize_sshifte, &
        trvsph, vshifte, vshifti, initialize_vshifte

    use gradient_routines, only: &
        gradec, grades, gradgc, gradgs, &
        igradec, igrades, igradgc, igradgs

    use module_idvtec, only: &
        idvtec

    use module_idvtes, only: &
        idvtes

    use module_idvtgc, only: &
        idvtgc

    use module_idvtgs, only: &
        idvtgs

    use icosahedral_geodesic_routines, only: &
        ihgeod

    use module_isfvpec, only: &
        isfvpec

    use module_isfvpes, only: &
        isfvpes

    use module_isfvpgc, only: &
        isfvpgc

    use module_isfvpgs, only: &
        isfvpgs

    use module_sfvpec, only: &
        sfvpec

    use module_sfvpes, only: &
        sfvpes

    use module_sfvpgc, only: &
        sfvpgc

    use module_sfvpgs, only: &
        sfvpgs

    use scalar_analysis_routines, only: &
        initialize_shaec, shaec, shaeci, &
        initialize_shaes, shaes, shaesi, &
        initialize_shagc, shagc, shagci, &
        initialize_shags, shags, shagsi

    use scalar_projection_routines, only: &
        shpe, shpei, &
        shpg, shpgi

    use scalar_synthesis_routines, only: &
        initialize_shsec, shsec, shseci, &
        initialize_shses, shses, shsesi, &
        initialize_shsgc, shsgc, shsgci, &
        initialize_shsgs, shsgs, shsgsi

    use scalar_laplacian_routines, only: &
        slapec, slapes, slapgc, slapgs, &
        islapec, islapes, islapgc, islapgs

    use vector_analysis_routines, only: &
        initialize_vhaec, vhaec, vhaeci, &
        initialize_vhaes, vhaes, vhaesi, &
        initialize_vhagc, vhagc, vhagci, &
        initialize_vhags, vhags, vhagsi

    use vector_synthesis_routines, only: &
        initialize_vhsec, vhsec, vhseci, &
        initialize_vhses, vhses, vhsesi, &
        initialize_vhsgc, vhsgc, vhsgci, &
        initialize_vhsgs, vhsgs, vhsgsi

    !    use module_visequ, only: &
    !        visequ
    !
    !    use module_visgau, only: &
    !        visgau
    !
    !    use module_visgeo, only: &
    !        visgeo

    use vector_laplacian_routines, only: &
        vlapec, vlapes, vlapgc, vlapgs, &
        ivlapec, ivlapes, ivlapgc, ivlapgs

    use vorticity_routines, only: &
        vrtec, vrtes, vrtgc, vrtgs, &
        ivrtec, ivrtes, ivrtgc, ivrtgs

    use colatitudinal_derivative_routines, only: &
        vtsec, vtseci, initialize_vtsec, &
        vtses, vtsesi, initialize_vtses, &
        vtsgc, vtsgci, initialize_vtsgc, &
        vtsgs, vtsgsi, initialize_vtsgs
    
    use type_FastFourierTransform, only: &
        FastFourierTransform

    use type_RealPeriodicFastFourierTransform, only: &
        RealPeriodicFastFourierTransform, &
        hrffti, hrfftf, hrfftb

    use type_AssociatedLegendrePolynomialGenerator, only: &
        AssociatedLegendrePolynomialGenerator, &
        alfk, lfp, lfpt, lfim, lfin

    use type_Vector3D, only: &
        Vector3D, &
        Vector3DPointer, &
        assignment(=), &
        operator(*)
    
    use type_Sphere, only: &
        Sphere

    use type_GaussianGrid, only:&
        GaussianGrid
        
    use type_GaussianSphere, only: &
        GaussianSphere
    
    use type_RegularGrid, only: &
        RegularGrid

    use type_RegularSphere, only: &
        RegularSphere

    use type_SpherepackUtility, only: &
        SpherepackUtility

    ! Explicit typing only
    implicit none
    
    ! Everything is private unless stated otherwise
    private

    ! Constants
    public :: wp, ip
    public :: HALF_PI, PI, TWO_PI

    ! Derived data types
    public :: FastFourierTransform
    public :: RealPeriodicFastFourierTransform
    public :: AssociatedLegendrePolynomialGenerator
    public :: GaussianGrid
    public :: GaussianSphere
    public :: RegularGrid
    public :: RegularSphere
    public :: Sphere
    public :: Vector3D
    public :: Vector3DPointer
    public :: assignment(=)
    public :: operator(*)

    ! Colatitude derivative
    public :: vtsec, vtseci, initialize_vtsec
    public :: vtses, vtsesi, initialize_vtses
    public :: vtsgc, vtsgci, initialize_vtsgc
    public :: vtsgs, vtsgsi, initialize_vtsgs

    ! Gradient
    public :: gradec, grades, gradgc, gradgs

    ! Inverse gradient
    public :: igradec, igrades, igradgc, igradgs

    ! Divergence
    public :: divec, dives, divgc, divgs

    ! Inverse divergence
    public :: idivec, idives, idivgc, idivgs

    ! Vorticity
    public :: vrtec, vrtes, vrtgc, vrtgs

    ! Inverse vorticity
    public :: ivrtec, ivrtes, ivrtgc, ivrtgs

    ! Gaussian wts & pts
    public :: compute_gaussian_latitudes_and_weights

    ! Geo/math coordinate transfers
    public :: geo2maths, math2geos, geo2mathv, math2geov

    ! Multiple ffts
    public :: hrffti, hrfftf, hrfftb

    public :: idvtec, idvtes, idvtgc, idvtgs

    public :: ihgeod
    public :: isfvpec, isfvpes, isfvpgc, isfvpgs
    public :: islapec, islapes, islapgc, islapgs
    public :: ivlapec, ivlapes, ivlapgc, ivlapgs

    public :: sfvpec, sfvpes, sfvpgc, sfvpgs

    ! Scalar analysis routines
    public :: shagc, shagci, initialize_shaec
    public :: shaes, shaesi, initialize_shaes
    public :: shaec, shaeci, initialize_shagc
    public :: shags, shagsi, initialize_shags

    ! Scalar synthesis routines
    public :: shsgc, shsgci, initialize_shsec
    public :: shses, shsesi, initialize_shses
    public :: shsec, shseci, initialize_shsgc
    public :: shsgs, shsgsi, initialize_shsgs

    ! Vector analysis routines
    public :: vhagc, vhagci, initialize_vhaec
    public :: vhaes, vhaesi, initialize_vhaes
    public :: vhaec, vhaeci, initialize_vhagc
    public :: vhags, vhagsi, initialize_vhags

    ! Vector synthesis routines
    public :: vhsgc, vhsgci, initialize_vhsec
    public :: vhses, vhsesi, initialize_vhses
    public :: vhsec, vhseci, initialize_vhsgc
    public :: vhsgs, vhsgsi, initialize_vhsgs

    public :: shpe, shpei
    public :: shpg, shpgi

    public :: slapec, slapes, slapgc, slapgs
    public :: vlapec, vlapes, vlapgc, vlapgs

    ! Grid transfer routines
    public :: trssph, sshifte, sshifti, initialize_sshifte
    public :: trvsph, vshifte, vshifti, initialize_vshifte

    public :: alfk, lfp, lfpt, lfim, lfin

    ! Temporary solution for testing
    public :: vecout, iout, name, vout, check_error

    interface vecout
        module procedure vecout_array
        module procedure vecout_scalar
    end interface vecout

    interface vout
        module procedure vout_array
        module procedure vout_scalar
    end interface vout

contains

    subroutine vecout_array(vec, nam, vec_size)

        ! Dummy arguments
        integer(ip),      intent(in) :: vec_size
        real(wp),         intent(in) :: vec(vec_size)
        character(len=*), intent(in) :: nam

        ! Local variables
        integer(ip) :: i

        write (stdout, 109) nam, (vec(i), i=1, vec_size)
109     format(1h a4, /(1h 8e11.4))

    end subroutine vecout_array

    subroutine vecout_scalar(vec, nam, vec_size)

        ! Dummy arguments
        real(wp),         intent(in) :: vec
        integer(ip),      intent(in) :: vec_size
        character(len=*), intent(in) ::  nam

        write (stdout, '(a, 8e11.4)') nam, vec

    end subroutine vecout_scalar

    subroutine vout_scalar(var, nam)

        ! Dummy arguments
        real(wp),         intent(in) :: var
        character(len=*), intent(in) :: nam

        write (stdout, '(a, e12.5)') nam, var

    end subroutine vout_scalar

    subroutine vout_array(var, nam)

        ! Dummy arguments
        real(wp),         intent(in) :: var(:)
        character(len=*), intent(in) :: nam

        write (stdout, *) nam, var

    end subroutine vout_array

    subroutine iout(ivar, nam)

        ! Dummy arguments
        integer(ip),      intent(in) :: ivar
        character(len=*), intent(in) :: nam

        write (stdout, '(a, i5)') nam, ivar

    end subroutine iout

    subroutine name(routine_name)

        ! Dummy arguments
        character(len=*), intent(in) :: routine_name

        write (stdout, '(a)') routine_name

    end subroutine name

    subroutine check_error(ierror)

        ! Dummy arguments
        integer(ip), intent(in) :: ierror

        if (ierror /= 0) write (stderr, '(a, i5)') '   ierror', ierror

    end subroutine check_error

end module spherepack
