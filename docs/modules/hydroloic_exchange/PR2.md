# PR2 — Runtime soil state and forcing inputs

**Status:** Full I09y replay validated; calculation-only SCHISM lifecycle, diagnostics, and capacity hotstart implemented
**Purpose:** Make the parity-tested soil calculation available during a SCHISM run while retaining a reference mode.

## Scope

- Add an off-by-default HydrologicExchange run option using SCHISM configuration conventions.
- Define a persistent soil-active element mask independent of wet/dry status. Inundation must not erase soil identity.
- Represent available soil capacity explicitly per active element, with units and ownership documented. Support either deriving it from layer soil moisture, porosity, and thickness or supplying a precomputed capacity. In the current NWM path, derive it from `SOIL_M`, `SMCMAX`, and `DZS` unless a compatible direct capacity field is provided.
- Use the same 0.4 m control depth as the current Python driver: activate native NWM layers 1 and 2 (0.1 m + 0.3 m) for available-capacity calculation. Keep conductivity and freeze gating based on the top layer. The existing Fortran kernel accepts the active layer count; PR2 will pass two for this configuration. Retain PR1's one-layer fixture as its original regression case.
- Provide a normalized precipitation input so the soil calculation can receive rain from either NWM `QRAIN` or SCHISM `sflux`. Keep snowmelt, ET, and groundwater-related inputs explicit rather than hiding them inside a generic net-flux field.
- Preserve an NWM replay mode for accumulated `ACSNOM`, `ACCET`, and `UGDRNOFF`, with explicit timestamp handling and the previous accumulation needed to calculate restart-safe tendencies.
- Define the groundwater process boundary in terms of inputs/state and returned fluxes. Implement the NWM replay behavior only; do not add placeholder process-based groundwater physics. A later process-based implementation may own groundwater storage/head and calculate recharge or discharge through the same exchange boundary.
- Retain a calculation-only runtime mode in this stage. Do not add distributed runoff to `vsource` until PR3.

## Out of scope

- No live continuity/source application.
- No implementation of process-based groundwater transport or NextGen coupling.
- No new hydrofabric mapping framework.

## Runtime integration points

- `schism_init`: allocate and initialize local soil-active state, select input mode, and load or derive initial capacity.
- `schism_step`: update forcing/state at the configured cadence and compute element-oriented diagnostics and fluxes.
- Standard SCHISM NetCDF output and hotstart conventions: add only fields needed to inspect or restore the new runtime state.

## Source structure for review

```text
src/HydrologicExchange/
    AGENTS.md
    CMakeLists.txt
    hydrologic_exchange_soil.F90       # existing PR1 state and soil calculations
    hydrologic_exchange.F90            # runtime lifecycle, input coordination, output/restart
    hydrologic_exchange_nwm.F90        # NWM parameters and forcing records in SI units
    hydrologic_exchange_mapping.F90    # calculated NWM fluxes to local SCHISM elements
    hydrologic_exchange_run.F90        # SCHISM configuration, local mask/geometry, runtime instance
```

The runtime files are built into the existing `hydrologic_exchange` library. The PR1 calculation remains in place. `schism_init`, `schism_step`, and `schism_finalize` now initialize, advance, and release the calculation. The off-by-default `OPT` option `ihydro_exchange=1` selects prescribed NWM replay. It does not change hydrodynamic sources.

The I09y path is NWM input records -> two-layer capacity and existing soil calculation on NWM cells -> linear interpolation of calculated fluxes -> element output for comparison. The sflux entry point documents a separate element-local input path. There is no process-based groundwater stub; the current soil kernel already accepts the prescribed drainage flux through its q_perc argument.

The component API accepts owned element IDs, projected centers, areas, mask, and an NWM subset rectangle. The SCHISM wrapper now prepares these internally: it reads `soil_active.prop`, selects resident elements (excluding ghosts), averages vertex longitude/latitude consistently with the Python reference, and converts centers to the native CONUS NWM spherical Lambert coordinates. It uses SCHISM's initialized areas and a local NWM rectangle around active elements. A rank with no active soil supplies zero diagnostics without opening NWM datasets. Each active rank reads its rectangle independently; no new MPI collectives are introduced. The reader handles the actual NWM dimension order, packed scale/offset, fill values, and mm-to-SI conversions. It accepts a chronologically ordered list of LDASOUT files.

`initialize_hydrologic_exchange` accepts an enable flag, file paths, geometry, and optional precomputed NWM-cell capacity. Otherwise capacity is derived from the two prescribed moisture layers. Times default to seconds relative to the first supplied file; an explicit epoch-second origin supports SCHISM's run start and I09y's August 1 origin. `step_hydrologic_exchange` computes at native NWM intervals and linearly interpolates mapped results between records. It can reconstruct any bracket directly from that record and its predecessor.

Standalone comparison output is a separate NetCDF file with time, global element IDs, areas, `q_surf`, `q_perc`, capacity, infiltration, and diagnostic `vsource = area*q_surf`. Creation refuses to overwrite an existing file. A caller/rank must supply its own output path. The sflux helper provides element-mean rain in m/s for the separate element-local pathway; selecting that provider in a full run is deferred.

The only additional hotstart field is `hydro_exchange_capacity`, an element array in meters. It uses the existing resident-element dimension, so `combine_hotstart7` combines it without special scalar registration and permits repartitioning. SCHISM already saves time; no separate HydroExchange timestamp is written. Restart reconstructs NWM forcing brackets, accumulation predecessors, and NWM-cell capacity from the original inputs, then reads the saved element capacity. At subsequent steps, prescribed NWM snapshots continue to determine capacity. This checkpoint is not a prognostic soil update and cannot replace the cell inputs needed for reference replay. Input contents, model time origin, mask, and any explicitly supplied cell capacity must remain consistent. If an older SCHISM hotstart has no capacity field, capacity is derived from NWM. For `ihot=1`, the time-reset calculation is initialized at the new run origin. A future prognostic soil/groundwater calculation must save its actual controlling state; element capacity alone may not suffice.

## Running inside SCHISM

Set `ihydro_exchange=1` in `param.nml`'s `OPT` namelist and supply `hydro_exchange.nml` using the example in `sample_inputs/`. The `HYDRO_EXCHANGE_INPUT` namelist specifies `soil_file`, `wrf_file`, `file_list`, and `mask_file`. Paths may be absolute or relative to the SCHISM input directory. The file list starts with its record count, followed by one LDASOUT path per line in increasing time order. Include a predecessor before the run start and records covering the full simulation/restart interval. `soil_active.prop` contains one global element ID and a 0/1 flag per mesh element; choose valid NWM-covered soil elements independently of the evolving wet/dry status.

This first live mode requires `ics=2` and the native CONUS NWM projection/grid. The time origin follows SCHISM's `start_year`, `start_month`, `start_day`, and `(start_hour + utc_start)` convention. It uses NWM `QRAIN`, prescribed soil snapshots, `ACSNOM`, `ACCET`, and `UGDRNOFF`. Native 0.1 + 0.3 m layers determine capacity, with top-layer conductivity/freeze gating. Direct capacity remains an API input; selecting alternative precipitation or capacity providers through the runtime namelist is deferred. The existing `sflux` rain conversion helper remains available for a later element-local calculation path.

At ordinary `nspool` output times, the standard SCHISM element outputs include `hydroSurfaceRunoff`, `hydroDrainage`, `hydroSoilCapacity`, `hydroInfiltration`, and `hydroSoilActive`. Fluxes are in m/s, capacity in m, and the mask is dimensionless. All fields are zero on inactive elements. Missing active NWM inputs retain the reference's missing-data behavior. Scribe metadata includes units; legacy `OLDIO` uses the existing `writeout_nc` path. The separate component replay file remains available for standalone comparisons. No calculated runoff, infiltration, or drainage enters `vsource` in PR2; reach-based `qBucket` remains separate.

## Subroutine validation

The component and `pschism` were initially built in `/tmp/schism-hydrologic-exchange-femto-build-2` using Intel/OpenMPI/NetCDF. Both PR1 comparison commands passed against the unchanged Python reference. The current persistent build and full-run validation are listed under SCHISM integration evidence below.

The runtime comparison uses five NWM records on August 29 at 09, 12, 15, 18, and 21 UTC. Its default sample is eight peak-runoff sources; `--sample-size 32` adds evenly spaced source IDs. It calls the unchanged Python postprocessor rather than duplicating the equations in the test. It also checks time origin, output, hotstart reconstruction and continuation, supplied capacity, disabled mode, sflux conversion, soil masks, and mapping boundaries/missing data. The 32-element comparison passes with maximum q_surf difference `5.082197683525802e-21` m/s and maximum volume-rate difference `1.7763568394002505e-15` m3/s against both saved and reconstructed I09y.

```bash
cmake --build /tmp/schism-hydrologic-exchange-femto-build-2 --target pschism -j4
OMPI_MCA_btl=self,tcp OMPI_MCA_pml=ob1 python test/HydrologicExchange/compare_i09y_runtime.py \
  --build-dir /tmp/schism-hydrologic-exchange-femto-build-2 \
  --geometry /sciclone/schism10/feiye/CIROH/Models/STOFS-3D-Atlantic/Hydrologic_exchange_tests/I09y_python_reconstruction/mesh_geometry.npz \
  --sample-size 32
```

The earlier `1.77e-6` m3/s discrepancy was dominated by different projected areas supplied by the comparison harness. The cause is PROJ network availability: toggling only this setting reproduces both area sets. With networking disabled, the WGS84-to-ESRI:102008 transformation uses the inverse of NAD83-to-WGS84 (1), matching I09y. With networking enabled, it selects operation (23), including the Louisiana grid shift `us_noaa_lahpgn.tif`. For element 256817 these give `352031.6298560762` and `352031.8763473087` m2 respectively; for element 490704 they give `206588.4879798022` and `206588.42679961154` m2. The mesh and area formula are identical. The test now requires recorded geometry instead of letting the current projection environment choose the datum transformation. The controlled experiment and pipelines are saved as `fortran_diagnosis/diagnose_projection.py` and `projection_area_diagnosis.json` in the external test folder.

A smaller `1.4623e-9` m3/s difference remained after using the correct geometry. Python evaluates `1-IMPERV` in the original float32 parameter type before promotion by the freeze gate. Promoting IMPERV before subtraction reproduces this residual in the unchanged Python calculation. The NWM reader now supplies the pervious fraction evaluated at its stored input precision; the soil kernel accepts that optional provider-derived value and otherwise retains its rkind calculation. Infiltration-cap multiplication follows the Python operation order. No soil physics or Python reference was changed.

The first 09 UTC subset record is excluded from the saved comparison because its replay accumulation tendency is deliberately zero. `--skip-saved-check` remains available for diagnosis; the default command enforces strict baseline parity.

The corrected Fortran replay has also passed the full I09y source comparison: 1,160,344 source elements x 441 records = 511,711,704 volume rates. All IDs and native timestamps agree. Maximum absolute volume-rate difference is `4.746314452575007e-12` m3/s, RMSE is `5.642649934711366e-16` m3/s, and the maximum runoff-depth-rate difference recovered from source volumes and the common areas is `1.595309561857089e-15` m/s. No values fail `rtol=1e-8`, `atol=1e-9` m3/s. Both integrated rectangle-rule volumes are `6,652,136,072.677826` m3 at reported precision, with zero difference.

Full results and the driver are saved under `/sciclone/schism10/feiye/CIROH/Models/STOFS-3D-Atlantic/Hydrologic_exchange_tests/fortran_diagnosis/`: `full_comparison_summary.json`, `full_comparison_by_time.csv`, `run_full_replay.py`, `full_replay_provenance.json`, and `full_replay/runoff.nc`. This checks the complete saved source list using recorded Python geometry; it does not independently derive a new Fortran source-selection mask or run the hydrodynamic solver. Replay uses cached forcing brackets and compressed diagnostic output for the full-period test. The separate small SCHISM integration test exercises the live lifecycle described above.

## Native SCHISM areas

Runtime initialization now defaults to `schism_glbl::area` in the caller's local
element order, after `aquire_hgrid` has initialized geometry. An explicit area
argument remains available for reproducing the saved Python product. Source IDs
are global identifiers; they must not be used to index the local area array.

For geographic meshes, `aquire_hgrid` converts nodes to Cartesian coordinates,
averages the Cartesian vertices, constructs a local element frame, and sums
`signa` areas for triangles `(1,2,3)` and, for quadrilaterals, `(1,3,4)`. Its frame
calculation is now shared as `compute_element_frame`, with the existing equations
unchanged. The comparison driver links these actual SCHISM routines. R09y's
`param.nml` specifies `ics=2` and both Earth radii as `6378206.4` m.

Native areas have been computed for all 1,537,103 mesh elements and compared on
the complete 1,160,344-source I09y list. They are about 0.303–0.393% larger than
the Python Albers areas (median 0.339%). Reusing the previously validated Fortran
`q_surf` for all 441 records isolates this geometric effect: total rectangle-rule
volume increases from about `6.652136073e9` to `6.674806222e9` m3, or 0.340795%.
The largest individual raw `vsource` difference is about `0.047095` m3/s. This is
a different area convention, separate from the much smaller PROJ datum-operation
difference diagnosed above. Native-area volume rates should therefore be checked
against `I09y_vsource * native_area / I09y_area`, rather than requiring identical
unadjusted volumes. The full area-adjusted comparison passes at `rtol=1e-8`,
`atol=1e-9` m3/s, with maximum error about `4.8e-12` m3/s.

The geometry methods explain this systematic difference: Python uses an Albers
equal-area projection on a reference ellipsoid, while this SCHISM configuration
uses local tangent-plane areas on a sphere. Equal-area preservation on the
projection's ellipsoid does not imply agreement with SCHISM's spherical areas.
The much smaller PROJ datum-transformation difference is a separate effect.
With common areas, the Python and Fortran runoff calculations agree to roundoff.

**Accepted for PR2 validation:** retain native SCHISM areas and accept the
0.340795% integrated-volume difference from the saved I09y product. Validate
soil physics using runoff depth rates and area-adjusted volume rates. This
closes the area discrepancy; it does not validate runoff application to hydrodynamics.

An independent runtime replay also recalculates the NWM soil fluxes for 32 source
elements and five Ida records, initializing through the default SCHISM area
array. It matches unchanged Python runoff depth rates to `1.7e-21` m/s and passes
the area-adjusted saved-volume comparison, restart continuation, and existing
runtime checks. The PR1 soil regression also passes. No Python reference physics,
source selection, interpolation locations, or original I09y products are changed.
The full-period comparison reuses validated depth rates; it is not a second
full-period NWM replay or a hydrodynamic run.

Reproducible tools are `test/HydrologicExchange/native_area_driver.F90`,
`compare_native_areas.py`, and `compare_i09y_runtime.py` (use `--geometry` for the
native mesh geometry and `--reference-geometry` for the recorded Python geometry).
Results, native geometry, build logs, and per-record totals are under
`/sciclone/schism10/feiye/CIROH/Models/STOFS-3D-Atlantic/Hydrologic_exchange_tests/native_area_comparison/`.

## Verification checklist

- [x] Confirm calculation-only mode leaves hydrodynamics unchanged: enabled/default-off elevation and all existing hydrodynamic hotstart fields are identical in the MPI basin test; default-off needs no HydroExchange inputs.
- [x] Run reference mode against the unchanged Python calculation and retain the PR1 regression.
- [x] Calculate using the agreed 0.4 m configuration and compare against saved Ida forcing, preserving calculation/mapping order and documenting native-area differences.
- [x] Verify NWM units/timestamps and the separate `sflux` rain conversion helper. Runtime selection of `sflux` is deferred; NWM replay uses `QRAIN` only.
- [ ] Check missing records, accumulation resets, first-record behavior, and forcing update cadence.
- [x] Check resident-element ownership, independent soil-active mask, and a rank with no active soil.
- [x] Verify a split SCHISM run agrees with a continuous run for prescribed-state mode, including repartitioning from two compute ranks to one.
- [x] Document which values are read/derived and distinguish prescribed NWM snapshots from prognostic state.

The remaining forcing-edge checklist item requires targeted tests for missing records and accumulation resets. Existing tests cover the first-record convention, native cadence, interpolation, and missing spatial data; no missing-record recovery is implemented.

### SCHISM integration evidence

`test/HydrologicExchange/check_schism_runtime.py` creates a 32-element closed basin around a valid coastal I09y source and uses the original August 29 NWM inputs read only. The validated run uses two compute ranks plus one scribe, 108 steps of 300 s, and 36 ordinary output records. Runoff is nonzero. The enabled and default-off runs have identical elevations and identical values for every existing hydrodynamic hotstart variable. The default-off case omits the HydroExchange option and removes all HydroExchange input files.

A split run checkpoints at step 45, between native NWM records, combines capacity through the unchanged `combine_hotstart7`, and resumes on one compute rank plus one scribe. All continuation outputs, including capacity and the mask, equal the continuous outputs exactly. Another run has one compute rank with no active soil. The component replay separately verifies capacity checkpoint restoration exactly. Both ordinary scribe output and a two-rank, 45-step legacy `OLDIO` run write the five diagnostics and capacity hotstart successfully.

Against unchanged Python, maximum differences in standard single-precision SCHISM output are `4.44e-13` m/s for surface runoff, `2.23e-16` m/s for drainage, and `1.43e-14` m/s for infiltration. The double-precision component replay retains its earlier tighter agreement with the saved I09y baseline after area adjustment.

Artifacts are under `/sciclone/schism10/feiye/CIROH/Models/STOFS-3D-Atlantic/Hydrologic_exchange_tests/PR2_schism_runtime/`: `validated/summary.json`, its five run directories, `validated.log`, `component_comparison.json`, `oldio_summary.json`, and `oldio_smoke/`. The scribe build is `native_area_comparison/build`; the separate legacy build is `PR2_schism_runtime/oldio_build`.

```bash
python test/HydrologicExchange/check_schism_runtime.py \
  --build-dir /sciclone/schism10/feiye/CIROH/Models/STOFS-3D-Atlantic/Hydrologic_exchange_tests/native_area_comparison/build \
  --work-dir /path/to/new/test_directory
```

Use the Intel/OpenMPI/NetCDF environment for the build and the `stofs` Python environment for the comparison. This test runs the hydrodynamic solver but does not apply the calculated soil fluxes to its water budget.

## Saved Ida validation baseline

The Python workflow has now been reconstructed in `/sciclone/schism10/feiye/CIROH/Models/STOFS-3D-Atlantic/Hydrologic_exchange_tests/I09y_python_reconstruction/`, using unchanged snapshots of the repository's `physics.py` and `nwm_io.py`, the 0.4 m setting, and all 441 NWM records. Predecessor-overlap batching was checked for exact agreement with a single four-record calculation before regenerating the complete product.

All 1,160,344 source IDs and all 441 records match the original I09y product. The full NetCDF volume-rate comparison has maximum absolute difference `4.746203430272544e-12` m3/s and RMSE `5.640857563228592e-16` m3/s, with no failures at `rtol=1e-8`, `atol=1e-9` m3/s. Integrated rectangle-rule volumes both equal `6,652,136,072.677826` m3 at reported precision. Mass-source and dummy-sink arrays match exactly. `source_sink.in`, `vsource.th`, and `msource.th` are byte-identical, confirmed by matching SHA-256 hashes of every complete file. See `comparison_summary.json`, `comparison_by_time.csv`, and `text_file_comparison.json` in the reconstruction folder for full results and file hashes. This establishes the saved Python product as a reproduced baseline; the corrected full Fortran replay now passes as described above.

Use the existing forcing in `/sciclone/schism10/feiye/STOFS3D-v8/I09y/Source_sink/source_soil_framework/` as the PR2 acceptance baseline. Its `source.nc`, `vsource.th`, `source_sink.in`, and `msource.th` files are dated April 11, 2026. Keep these files intact; write new comparison products separately.

The current `nwm_io.py` records this output directory, the `I09y/hgrid.gr3` mesh, and NWM inputs under `NWM/CONUS/netcdf/LDASOUT/for_2021_Ida/`. That configuration is evidence of the workflow, but does not establish the exact script version used for the saved output. The agreed PR2 configuration follows the current Python driver's `depth=0.4`, selecting native layers 1 and 2 for storage; conductivity and freeze gating use layer 1. PR1's native-layer-1 (0.1 m) comparison remains a separate regression case. Compare the 0.4 m result against the saved baseline later and investigate any discrepancy.

For the Ida comparison:

- Match the original NWM records, static soil fields, active layers, mesh, spatial interpolation, source-element IDs, timestamps, and area conversion.
- Preserve the reference order of calculating runoff on the NWM grid and then interpolating runoff to element centers. Calculating from interpolated soil inputs on elements is a different discretization and requires a separate comparison.
- Compare the newly generated volume rates against saved `vsource` in m3/s by element ID and timestamp. For a depth-rate comparison, recover `q_surf` using the same element areas used by the original preprocessor.
- Report maximum absolute/relative differences, discrepancies by element/time, and integrated volume differences. Set the acceptance tolerance from the stored precision and numerical operations, and investigate any larger discrepancy.
- Describe the saved baseline's calculation mode explicitly. Matching it is required before evaluating alternative layer configurations, precipitation sources, or prognostic soil updates.

## Decisions to close in this stage

- Current replay inputs are NWM `QRAIN`, `ACSNOM`, `ACCET`, soil snapshots, and `UGDRNOFF`. Alternative provider selection remains future work.
- Capacity is recomputed from prescribed snapshots over 0.4 m; infiltration does not independently reduce a prognostic store.
- Initial capacity is derived from those snapshots; the component API also accepts prescribed cell capacity.
- Only element capacity is added to hotstart. SCHISM time and the original forcing list reconstruct bracketing records and accumulation predecessors. Keep input contents and model origin consistent across restart.
- Fluxes and capacity are computed at native NWM record times and linearly interpolated to SCHISM steps. Tendencies use native NWM intervals, not the hydrodynamic timestep.

## Exit criteria

An off-by-default runtime calculation with documented input ownership, state/restart behavior, diagnostics, parity against PR1, and matching runoff for the saved Ida forcing within a documented numerical tolerance. This stage still does not route `q_surf` through the continuity equation.
