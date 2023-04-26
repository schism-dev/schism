#Theory

The single-class ice/snow module is taken from FESOM2, courtesy of Dr. Sergey Danilov's group. 

Dynamical part solves the ice motion equation:

Here m is the ice plus snow mass per unit area, Cd the ice-ocean drag coefficient, ρo the
water density, a the ice compactness, u and uo the ice and ocean velocities, τ the wind
stress, H the sea surface elevation, g the acceleration due to gravity and Fj = ∂iσij is
the force from stresses within the ice. Here we use Cartesian coordinates for brevity, with
i, j = 1, 2 implying x and y directions, and the implementation of spherical coordinates will
be discussed further. Summation over repeating coordinate indices is implied. The stress
tensor is symmetric. The mass m is the combination of ice and snow contributions

with ρice and ρs, respectively, the densities of ice and snow and hice and hs their mean
thicknesses (volumes per unit area).
For the visco-plastic (VP) rheology (Hibler, 1979) one writes

where

is the deformation rate tensor, η and ζ are the moduli (‘viscosities’) and P the pressure.
Both the stress and deformation rate tensors are symmetric, so they are characterized by
three independent components. The standard VP rheology adopts the following scheme of
computing pressure P and moduli η and ζ:


 are default values, they are adjusted in practice. In this scheme, 
∆min serves for viscous regularization of plastic behavior in areas where ∆ is very small. We note that recent multi-category
ice implementations (such as CICE, see Hunke and Lipscomb 2008) use different parameterization for P0, which takes into account the distribution of ice over thickness categories. This
does not change the basic equations (1, 2).
In our case we deal with three tracers, the area coverage a, ice volume hice and snow
volume hs. They are advected with ice velocities and modified through thermodynamical
forcing

with Sa, Sice and Ss the sources due to the exchange with atmosphere and ocean.

#Usage
