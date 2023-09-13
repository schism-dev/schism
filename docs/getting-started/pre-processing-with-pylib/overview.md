**pylib** provides

* a matlab-style platform for python usage
* and an independent python-based SCHISM modeling workflow

please refer to `https://github.com/wzhengui/pylibs` for more information.

<br>
**some of pylib basic functions**

* database usage
* time manipulation, geometry handling (`inside_polygon, near_pts, compute_contour`)
* data analysis: `low/band-pass filters, running smooth, FFT, statistics, least-square-fit, EOF, harmonic analysis, etc`.
* handling different file formats: `SCHISM inputs/outputs, ASCII files, netcdf, shapefile, Excel, DEM files, projections, etc.`

<br>
**some of SCHISM related functions**

* geometry information: `node, element, side, etc.`
* hgrid: `plot, interplation, boundary, read/save, grid-preprocess (grd2sms, sms2grd, quads check/split, skew-element check, projection)`
* vgrid: `compute_zcor, etc.`
* point/region files
* post-process: `extract time series, profiles, slabs, fluxes, etc.`
* visualization of SCHISM inputs/outputs
