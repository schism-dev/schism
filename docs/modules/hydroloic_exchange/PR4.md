# PR4 — Operational NWM mapping and regional comparison

**Status:** Planned
**Purpose:** Evaluate the runtime element-based calculation against existing NWM-to-SCHISM source strategies on a realistic coastal transition-zone case.

## Scope

- Keep reach-based `qBucket` mapped to resolved SCHISM river sources in all comparisons. Do not silently replace it with local `q_perc`; establish whether their budgets overlap and document the result.
- Compare three configurations using common forcing periods and documented temporal/spatial mappings:
  1. Reach-based `qSfcLatRunoff + qBucket`.
  2. Externally prepared distributed `q_surf + qBucket`.
  3. Native runtime HydrologicExchange distributed `q_surf + qBucket`.
- Begin with a manageable subdomain and advance to Ida 2021 or another regional case after the small case passes.
- Compare total input volumes, runoff hydrographs and arrival timing, wet area/ponding, overland pathways, and inland/coastal responses.
- Document the hydrofabric mapping already available and identify the catchment, flowpath, and nexus links needed to map to SCHISM land elements, channel elements, and network connections.
- Build only mapping needed for the operational comparison; do not create a general mapping framework without a concrete use case.

## Groundwater and source accounting

Keep NWM `qBucket` a distinct reach-based river source. Keep `UGDRNOFF`/`q_perc` as the soil-calculation input/output defined by the postprocessor unless a separately reviewed groundwater process changes that role. Check whether adding distributed `q_surf` alongside `qBucket` duplicates water in any configuration and describe the supported accounting explicitly.

## Verification checklist

- [ ] Use common time windows, input products, units, and documented interpolation/aggregation.
- [ ] Compare total volume and hydrographs for each source strategy.
- [ ] Compare runoff arrival at resolved channels and wetting/ponding patterns.
- [ ] Attribute differences to physics, spatial placement, input timing, or hydrodynamic routing.
- [ ] Quantify any overlap or mismatch between `qBucket` and the groundwater-related `q_perc` field.
- [ ] Record current hydrofabric coverage and unresolved coastal-transition mapping needs.
- [ ] Reproduce results from saved configurations and documented commands.

## Exit criteria

A reproducible subdomain or regional comparison with differences and budget accounting explained, plus a scoped mapping proposal for any future CTZ work.

## Later work, not part of PR4

Implement process-based groundwater physics and two-way NextGen exchange only after the one-way runoff behavior, budgets, state ownership, time synchronization, and restart semantics are stable. Treat those as separately reviewed design and implementation stages.
