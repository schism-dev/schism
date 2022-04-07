The 3D sediment model inside SCHISM (USE_SED, also referred to as ‘SED3D’ sometimes to differentiate it from the 2D sediment module ‘SED2D’) is adapted from Community Sediment Transport Model ([Warner et al. 2008](#warner2008)). The algorithm is implemented on UGs and we have also reworked several components and added a morphological module (Exner equation). Detailed are reported in [Pinto et al. (2012)](#pinto2012). Also note that for best results this module should be run in conjunction with the wave model (WWM) to account for the wave-enhanced bottom stress.

The main parameter input is sediment.in, which is a free-format file. The sample file has comments /explanations for each parameter; a few important ones are explained below:

- `sed_morph`: This parameter controls active morphology.
- `Nbed`: The # of bed layers affects the stability of the sorting, and `Nbed=1` works robustly with `sed_morph>0`.

Other inputs include `sed_class` and output flags in `param.nml`, and a few `.ic` files: `bedthick.ic`, `bed_frac_[1,2..].ic`, and I.C. for concentrations of all classes (`*_[hvar]_[1,2…].ic` or (`*_[vvar]_[1,2…].ic`). The B.C. inputs may include `SED_[1,2..].th`, `SED_3D.th.nc`. The nudging inputs may be `SED_nudge.gr3` and `SED_nu.nc`.

The outputs from SED3D can be combined and visualized just as other SCHISM outputs.

A common crash occurs when the active morphology is on. At a river inflow boundary, often the depths will decrease over time due to deposition there and eventually the boundary will become dry/blocked. A work-around is to use a combination of ‘bare-rock’ bed around the boundary (specified in `bedthick.ic`), and ‘clear-water’ inflow as B.C. (sediment concentration =0 in `bctides.in`) , and then input the incoming sediment concentration as point sources (`msource.th`) a distance away from the boundary.

**References**

<span id="pinto2012">Pinto, L., Fortunato, A.B., Zhang, Y., Oliveira, A. and Sancho, F.E.P. (2012) Development and validation of a three-dimensional morphodynamic modelling system, Ocean Mod., 57-58, 1-14.</span>
<span id="warner2008">Warner, J.C., Sherwood, C.R., Signell, R.P., Harris, C.K., Arango, H.G. 2008. Development of a three-dimensional, regional, coupled wave, current, and sediment-transport model. Comput. Geosci. 34 (10), 1284–1306. doi: 10.1016/j.cageo.2008.02.012</span>