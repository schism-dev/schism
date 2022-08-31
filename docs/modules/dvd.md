[Klingbeil et al. (2014)](#klingbeil2014) proposed a theory of estimating numerically induced mixing coefficients. At the moment we have only implemented this module for 1 tracer (salinity).

The I.C. for this module is hardwired to be 0 and the code will set it automatically. The B.C. flags should also be 0 in general, so is the nudging flag. The main output is `iof_dvd(1)` which corresponds to numerical mixing for salinity in $\text{PSU}^2/s$.

**References**

<span id="klingbeil2014">Klingbeil, K., Mohammadi-Aragh, M., Gräwe, U., Burchard, H. (2014) Quantification of spurious dissipation and mixing – Discrete variance decay in a Finite-Volume framework, Ocean Modelling, https://doi.org/10.1016/j.ocemod.2014.06.001.</span>