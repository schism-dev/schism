# **spherepack - A Fortran library of spherical harmonic transforms**

A modernization of NCAR's SPHEREPACK3.2.

* The original work, written in fixed-from FORTRAN 77, was heavily refactored to incorporate features of free-form modern Fortran (2008+). 
* Dynamic memory allocation procedures for analysis and synthesis are now hidden from the end user.

-----------------------------------------------------------------------------

## What is spherepack?

A collection of Fortran programs for computing common spherical differential operators including divergence, vorticity, latitudinal derivatives, gradients, the Laplacian of both scalar and vector functions, and the inverses of these operators.

For example, given divergence and vorticity, the package can be used to compute velocity components, then the Laplacian inverse can be used to solve the scalar and vector Poisson equations. The package also contains routines for computing the associated Legendre functions, Gaussian points and weights, multiple fast Fourier transforms, and for converting scalar and vector fields between geophysical and mathematical spherical coordinates.

Test programs are provided for solving these equations. Each program serves two purposes: as a template to guide you in writing your own codes utilizing the spherepack routines, and as a demonstration on your computer that you can correctly produce spherepack executables.

-----------------------------------------------------------------------------

## Usage

```fortran
        	
    use spherepack, only: &
        wp, & ! working precision
        GaussianSphere

    ! Explicit typing only
    implicit none
    
    type(GaussianSphere)  :: foo
    real(wp), allocatable :: scalar_function(:,:)
    real(wp), allocatable :: laplacian(:,:)
    real(wp), allocatable :: solution(:,:)
    
    ! Initialize object
    foo = GaussianSphere(nlat=19, nlon=36)
    
    !.... generate some data
    
    ! Compute complex spectral coefficients
    call foo%perform_complex_analysis(scalar_function)
    
    ! Compute laplacian on sphere
    call foo%get_laplacian(scalar_function, laplacian)
    
    ! Invert laplacian on sphere
    call foo%invert_laplacian(laplacian, solution)
    
    ! Release memory
    call foo%destroy()

```

-----------------------------------------------------------------------------

## Requirements

* The GNU Make tool https://www.gnu.org/software/make/
* The GNU gfortran compiler https://gcc.gnu.org/wiki/GFortran

-----------------------------------------------------------------------------


## To build the project

Type the following command line arguments

```bash

	git clone https://github.com/jlokimlin/spherepack.git; cd spherepack; make all
```

-----------------------------------------------------------------------------

## Contributing

This project is still a work in progress and anyone is free to contribute under the proviso that they abstain from using the dreaded **go to**.

For bug reports or feature requests please open an issue on github.

-----------------------------------------------------------------------------


## Bibliography

[1] Swarztrauber, Paul N. "On computing the points and weights for Gauss--Legendre quadrature." *SIAM Journal on Scientific Computing* 24.3 (2003): 945-954.

[2] Swarztrauber, Paul N., and William F. Spotz. "Generalized discrete spherical harmonic transforms." *Journal of Computational Physics* 159.2 (2000): 213-230.

[3] Adams, John C., and Paul N. Swarztrauber. "SPHEREPACK 3.0: A model development facility." *Monthly Weather Review* 127.8 (1999): 1872-1878.

[4] Swarztrauber, Paul N. "Spectral transform methods for solving the shallow-water equations on the sphere." *Monthly Weather Review* 124.4 (1996): 730-744.

[5] Williamson, David L., et al. "A standard test set for numerical approximations to the shallow water equations in spherical geometry." *Journal of Computational Physics* 102.1 (1992): 211-224.



