# PR2 — Runtime soil state and forcing inputs

**Status:** Planned
**Purpose:** Make the parity-tested soil calculation available during a SCHISM run while retaining a reference mode.

## Scope

- Add an off-by-default HydrologicExchange run option using SCHISM configuration conventions.
- Define a persistent soil-active element mask independent of wet/dry status. Inundation must not erase soil identity.
- Represent available soil capacity explicitly per active element, with units and ownership documented. Support either deriving it from layer soil moisture, porosity, and thickness or supplying a precomputed capacity. In the current NWM path, derive it from `SOIL_M`, `SMCMAX`, and `DZS` unless a compatible direct capacity field is provided.
- Keep the calculation multi-layer capable while setting the active layer count to one for reference mode. Multi-layer behavior requires separate verification.
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

## Verification checklist

- [ ] Confirm the option-off path leaves existing behavior unchanged.
- [ ] Run reference mode against the PR1 calculation for the same input records.
- [ ] Verify selected NWM and `sflux` rain sources have consistent units and timestamps.
- [ ] Check missing records, accumulation resets, first-record behavior, and forcing update cadence.
- [ ] Check MPI local-element ownership and soil-active mask consistency.
- [ ] Verify a split run agrees with a continuous run for prescribed-state mode.
- [ ] Document which values are read, derived, or prognostic; do not describe prescribed NWM snapshots as SCHISM-evolved state.

## Decisions to close in this stage

- Which source supplies each of rain, snowmelt, ET, soil state, and groundwater drainage.
- Whether capacity is recomputed from prescribed soil snapshots or prognostically reduced by infiltration.
- How initial capacity is set when no directly supplied capacity field exists.
- Which accumulated fields and previous values must be saved for a restart.
- The cadence and interpolation used when hydrologic inputs are less frequent than SCHISM steps.

## Exit criteria

An off-by-default runtime calculation with documented input ownership, state/restart behavior, diagnostics, and parity against PR1. This stage still does not route `q_surf` through the continuity equation.
