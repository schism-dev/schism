This module is not a tracer module and simulates the long-term migration of marshes under sea level rise. It is usually invoked together with SED (with optional morphological acceleration) and optionally with WWM; cf. [Nunez et al. (2020)](#Nunez2020).

The only parameter for this module is `slr_rate` (sea-level rise rate in mm/yr). The output flag is `iof_marsh`, which is an integer of either 0 (no marsh) or 1 (has marsh) at an element. Optionally, the user might also consider invoking the vegetation option in the code `isav` in conjunction with the marsh module to simulate the form drag and turbulence induced by the marsh vegetation.

Additional inputs are required to specify the I.C. for marsh extent (`marsh_init.prop`) as well as migration barrier info (`marsh_barrier.prop`); both use 0 or 1 to specify ‘on/off’.


**References**

<span id="Nunez2020">Nunez, K., Zhang, Y., Herman, J., Reay, W. and Hershner, C. (2020) A multi-scale approach for simulating tidal marsh evolution, Ocean Dynamics, https://doi.org/10.1007/s10236-020-01380-6</span>