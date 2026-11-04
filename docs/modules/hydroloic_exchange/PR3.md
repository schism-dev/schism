# PR3 — Apply distributed runoff through `vsource`

**Status:** Planned
**Purpose:** Route HydrologicExchange's element-distributed surface runoff with SCHISM hydrodynamics.

## Scope

- Feed element runoff through SCHISM's existing volume-source path. For an area-normalized rate, convert explicitly as `Q_e = area_e * q_surf_e` (m/s × m² = m³/s) and add it to `vsource` before the continuity solver consumes that array.
- Define source ownership by element. On soil-active land elements, HydrologicExchange supplies the net water entering surface routing; avoid also applying the same direct atmospheric precipitation/evaporation contribution through SCHISM's existing source calculation. Keep the existing atmospheric treatment on soil-inactive and open-water elements.
- Preserve other sources already accumulated in `vsource`, including reach-based river sources. Keep `qBucket` on its separate reach-to-river path; do not substitute local `q_perc` for it.
- Inspect and document behavior for dry elements, newly wetted elements, wet/dry transitions, and MPI ownership.
- Check freshwater tracer and heat handling for volume added through this path. Confirm ambient/source temperature and salinity behavior is appropriate instead of assuming volume-only changes are sufficient.
- Add output for the applied runoff and any source-accounting diagnostics using existing SCHISM output conventions.

## Water-accounting rule

`q_surf` includes the selected rain input, snowmelt, ET subtraction, infiltration, and capacity/freeze limitations as defined by the soil calculation. Therefore, the native net atmospheric contribution on a HydrologicExchange-active element must not be counted again. The implementation must document the selected provenance of rain, snowmelt, and ET and show that the element's water enters the surface-routing system once. Reach-based `qBucket` is kept distinct; its relationship to local NWM `UGDRNOFF` must be checked for budget overlap.

## Runtime integration points

- `schism_step`: add the area-scaled HydrologicExchange rate while assembling `vsource`, replacing the direct net atmospheric contribution only on active soil elements.
- Existing source application in the continuity calculation remains the point where the assembled `vsource` affects hydrodynamics.
- Standard output and hotstart paths carry the applied flux and any prognostic state introduced by PR2.

## Verification checklist

- [ ] Confirm the HydrologicExchange option-off path matches the pre-change baseline.
- [ ] Confirm precipitation is counted once on soil-active elements and native treatment remains on other elements.
- [ ] Verify total applied volume against `sum(q_surf * element_area * dt)`.
- [ ] Check dry/wet transitions and source application for local and shared elements across MPI partitions.
- [ ] Check salinity/temperature and other tracer behavior for freshwater additions.
- [ ] Compare a rain-on-land case and a wetting/drying case against the PR1 externally forced reference.
- [ ] Report volume-budget closure, timing, overland pathways, and timestep sensitivity.
- [ ] Keep `qBucket` independent and document any possible groundwater-budget overlap with `q_perc`.

## Exit criteria

A small SCHISM run routes distributed surface runoff through the mesh with precipitation counted once and documented volume and tracer-budget results. This does not yet establish operational regional mapping.
