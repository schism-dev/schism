# PR1 — Reproduce the soil postprocessor

**Status:** Complete (top-layer calculation parity; no runtime coupling)
**Purpose:** Establish a reproducible Fortran baseline for the existing soil calculation before connecting it to SCHISM hydrodynamics.

## Scope

- Add the smallest build-integrated `HydrologicExchange` soil interface that accepts forcing, soil state/parameters, and `dt`, then returns infiltration diagnostics, available capacity, `q_surf`, and the existing groundwater-related `q_perc` value.
- Keep inputs source-neutral at the calculation boundary. The parity case uses the same NWM `QRAIN`, `ACSNOM`, `ACCET`, `UGDRNOFF`, `SOIL_M`, `SOIL_W`, and static soil fields as the current postprocessor.
- Activate only the top soil layer for the reference comparison. Keep the Fortran interface able to represent multiple layers, but do not enable them in this parity case. Record the actual `DZS` thickness represented by the active layer.
- Leave the Python reference unchanged. Preserve its accumulation differencing behavior, including zero tendency at the first sample and clipping negative increments to zero.
- Preserve the postprocessor's operation order: calculate fluxes on the NWM grid, then linearly interpolate the resulting `q_surf` and `q_perc` to SCHISM element centers. A later element-native calculation is a separate spatial-discretization change because the partitioning is nonlinear.

## Reference findings

The calculation is in `src/Utility/Pre-Processing/STOFS-3D-Atl-Setup/src/stofs3d_setup/ops/Source_sink/Soil/physics.py` (`compute_fluxes`). Its NWM readers, unit conversion, layer selection, interpolation, and `vsource` preparation are in the adjacent `nwm_io.py`.

| Input | Source/reference meaning | Units at the calculation boundary |
| --- | --- | --- |
| `QRAIN` | LDASOUT liquid rain rate | m s⁻¹ |
| `ACSNOM`, `ACCET`, `UGDRNOFF` | LDASOUT accumulated melt, ET, and unconfined groundwater drainage | m; non-negative finite differences become m s⁻¹ |
| `SOIL_M`, `SOIL_W` | Total and liquid volumetric soil water by layer; only layer 1 enters the freeze gate | m³ m⁻³ |
| `SMCMAX`, `DZS` | Porosity and layer thickness used for available storage | m³ m⁻³; m |
| `DKSAT` (top layer), `IMPERV` | Conductivity rate cap and impervious fraction | m s⁻¹; fraction |
| `dt` | Interval between NWM records | s |

The calculation is stateless with respect to soil water: each timestamp's NWM `SOIL_M`/`SOIL_W` values are prescribed snapshots. The storage deficit is `S = sum_k(max(SMCMAX_k - SOIL_M_k, 0) * DZS_k)`. Accumulated-field tendencies use a zero first difference and clip negative increments (resets) to zero. Then `P_liq = max(QRAIN + dACSNOM/dt - dACCET/dt, 0)`, `q_perc = dUGDRNOFF/dt`, and `f_freeze = max(1 - max((SOIL_M_1 - SOIL_W_1)/(abs(SOIL_M_1) + 1e-6), 0), 0)^3`. The storage cap is `(S/dt + q_perc) * (1 - IMPERV) * f_freeze`; the infiltration cap is the minimum of that value and top-layer `DKSAT`; infiltration is `min(P_liq, infiltration cap)` and `q_surf = P_liq - infiltration`.

The calculation returns `q_surf` and `q_perc` in m s⁻¹ plus diagnostics for `P_liq`, infiltration, infiltration cap, freeze gate, rate cap, and effective storage cap. Available storage `S` is derived from prescribed soil snapshots; it is not an evolved state. ET contributes to `P_liq` but is not returned as an independent flux by this calculation.

`q_surf` and `q_perc` are calculated on the NWM grid in m s⁻¹, then linearly interpolated to SCHISM element centers. The current `nwm_io.main()` converts only mapped `q_surf` to `vsource` by multiplying by element area; mapped `q_perc` is diagnostic here, and is not reach-based `qBucket`. The `nwm_io` source construction omits all-zero and all-NaN elements. Smooth freezing uses the top layer's total and liquid water; `SNEQV` is not used in this mode. There is no soil-water update, groundwater transport, or soil restart state in this postprocessor.

The requested PR1 comparison activates only the first native NWM layer (`DZS[0] = 0.1 m` for the Ida inputs below). The existing script's `main()` calls `compute_from_nwm(..., depth=0.4)`, whose selector includes the first two NWM layers (`0.1 + 0.3 m`) for storage. The Python reference remains unchanged; the PR1 one-layer comparison calls its unchanged `compute_fluxes` routine with normalized one-layer inputs. It establishes top-layer calculation parity, not parity with the script's current 0.4 m default configuration.

## Minimal layout and later SCHISM call sites

- `src/HydrologicExchange/hydrologic_exchange_soil.F90` contains the soil state snapshot and calculation only. It accepts prepared flux rates and soil fields and returns diagnostics; it does not read NWM files or apply sources.
- `src/HydrologicExchange/CMakeLists.txt` builds the module; `src/Hydro/CMakeLists.txt` includes and links it when SCHISM's `Hydro` module is built. This keeps the link on the directory that owns target `hydro`, compatible with the repository's CMake 3.12 policy. PR1 has no `schism_init` or `schism_step` call.
- If PR2 adds runtime inputs/state, `schism_init` is the place for initialization and `schism_step` for the per-step calculation, following the existing routines. For PR3, `schism_step` assembles `vsource` from atmospheric precipitation/evaporation near lines 1981–1987 and later applies it to element water volume near line 6129. A soil-active-element rule must reconcile that precipitation term before adding distributed runoff, so precipitation is counted once. `qBucket` remains on its existing reach-based river-source path.

## Out of scope

- No live `vsource` contribution or change to the continuity equation.
- No prognostic soil-water update, groundwater transport, NextGen adapter, hydrofabric mapping, or two-way exchange.
- No new NWM file-reading path inside the SCHISM time step.

## Progress

- Added `src/HydrologicExchange/hydrologic_exchange_soil.F90` with an available-capacity state, capacity calculation for a caller-selected number of active layers, accumulated-field differencing, and the smooth-freeze soil flux partition.
- Added `test/HydrologicExchange/compare_soil_postprocessor.py` and a deterministic four-time, four-cell fixture. The unchanged Python routine and Fortran routine agree at `rtol=1e-12`, `atol=1e-14` for capacity and flux diagnostics. Applying the same xarray linear interpolation to three synthetic element centers also agrees; the fixture area-integrated `q_surf` volume is `0.108252885468 m^3`.
- Added `test/HydrologicExchange/soil_ida_toplayer_fixture.dat`, a small SI-normalized subset from NWM Ida files at 2021-08-29 12, 15, 18, and 21 UTC. It contains four adjacent cells at NWM indices `y=771:772`, `x=2970:2971`, including the 18 UTC rain peak; the comparison uses only native layer 1 and its 0.1 m thickness. At `rtol=1e-12`, `atol=1e-14`, the unchanged Python routine and Fortran agree at stored precision for all fields (maximum absolute difference zero), with nonzero surface runoff in 12 of 16 cell-times. Run it with `python test/HydrologicExchange/compare_soil_postprocessor.py --fixture test/HydrologicExchange/soil_ida_toplayer_fixture.dat --grid-only`.
- The Ida values come from `202108291200.LDASOUT_DOMAIN1`, `202108291500.LDASOUT_DOMAIN1`, `202108291800.LDASOUT_DOMAIN1`, and `202108292100.LDASOUT_DOMAIN1`, with static fields from `soilproperties_CONUS_FullRouting.nc` and `wrfinput_CONUS.nc`. Rows are ordered by time, then y, then x; values are converted from NWM mm s⁻¹/mm to m s⁻¹/m. No land or wet/dry mask is applied to this four-cell calculation fixture; all 16 cell-time input records are finite. The postprocessor's later source construction filters SCHISM elements that are all-zero or all-NaN in `q_surf`.
- On FEMTO with Intel 2024.0, OpenMPI 4.1.6, and NetCDF C/Fortran modules loaded, SCHISM CMake configuration succeeds. Builds of `hydrologic_exchange`, `hydro`, and `pschism` pass; `hydro` compiles `schism_init.F90` and `schism_step.F90` and links `libhydro.a` with `libhydrologic_exchange.a`, and `pschism` links successfully.
- The real-data comparison is at the soil calculation boundary. The synthetic fixture separately checks the linear interpolation operation and area conversion; no real SCHISM mesh mapping or full hydrodynamic run is part of PR1.

## Repository behavior to preserve

The input units, equations, spatial ordering, and assumptions are recorded above. `q_perc` is the `UGDRNOFF` tendency used by this calculation; it is not reach-based `qBucket` and is not itself applied as a SCHISM source by the soil calculation.

The available-capacity value is a soil-state snapshot for this mode. It may be computed from the active layer's `SMCMAX`, `SOIL_M`, and `DZS`, or supplied to the interface already computed. The current NWM path derives it from soil fields; no direct capacity field is identified in that path. PR1 does not evolve that state.

## Proposed code boundary

Use `src/HydrologicExchange/` for the calculation component, following the repository's CMake conventions. Keep the soil calculation separate from SCHISM source application. Do not add a standalone I/O module unless the parity harness demonstrates it is needed; the current Python preprocessor can prepare and map the comparison inputs.

## Verification checklist

- [x] Run `python test/HydrologicExchange/compare_soil_postprocessor.py` for the deterministic parity fixture.
- [x] Compare `q_surf`, infiltration, capacity, freeze gate, infiltration cap, and `q_perc` before spatial mapping on identical fixture inputs.
- [x] Apply the same linear grid-to-element-center interpolation and compare mapped fields on the fixture.
- [x] Compare area-integrated `q_surf` volume on the fixture.
- [x] Record real NWM input variables, units, timestamps, active-layer definition, masks, and spatial order for a data subset.
- [x] Run the Fortran and unchanged Python calculation on identical real NWM-grid inputs.
- [x] State numerical tolerance and report differences from precision, interpolation, or ordering.
- [x] Build the new component through SCHISM's configured build and build the `hydro` library that links it.
- [x] Confirm no runtime source or baseline hydrodynamic behavior changed; PR1 adds no call in `schism_init` or `schism_step`.

## Exit criteria

A documented build and comparison command, a reproducible fixture or input subset, quantitative parity results, and a clear accounting of any remaining differences. Passing PR1 establishes calculation parity only; it does not demonstrate live coupling.

## Decisions deferred

- How to reconcile the requested top-layer PR1 mode with `nwm_io.main()`'s current 0.4 m default, which selects two native layers. The Python reference remains unchanged; current parity is specifically for the top-layer mode.
- Which soil fields remain prescribed from NWM and which become SCHISM-owned state.
- Whether and how available capacity is updated by infiltration in a prognostic mode.
- Restart semantics for soil state and accumulated forcing.
- Whether rain, snowmelt, and ET are sourced together from NWM, SCHISM forcing, or a mixed configuration.
