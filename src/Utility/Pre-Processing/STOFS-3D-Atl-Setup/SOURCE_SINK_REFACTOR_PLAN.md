# Source/Sink Workflow Refactor Plan

## Objective

Restore a modular source/sink workflow in which each scientific operation can
be enabled, disabled, tested, and reviewed independently.

The refactor should separate responsibilities currently bundled inside
`patch_artificial_island_source_sink.py` while preserving existing numerical
behavior during the initial extraction steps.

## Working Branch

- Branch: `refactor/source-sink-workflow`
- Starting branch: `master`

## Refactoring Principles

- Work in small, reviewable steps.
- Add tests before or alongside each extraction.
- Separate mechanical code movement from behavior changes.
- Preserve function behavior, processing order, station mappings, fallback
  rules, diagnostics, and output formatting during extraction.
- Stop and investigate if an extraction changes output unexpectedly.
- Give each configuration option one clear responsibility.
- Keep artificial-island code limited to artificial-island corrections.
- Avoid simultaneous edits by multiple people to the same orchestration code.

## Current Responsibilities to Separate

The current artificial-island module includes:

- YAML parsing and validation
- Source/sink array manipulation
- Grid projection and spatial searches
- USGS downloading, caching, parsing, and interpolation
- Hudson-specific adaptive USGS downloading
- Automatic temperature correction for all eligible sources
- Delaware and Hudson flow and temperature overrides
- Region-based source zeroing
- Artificial-island source placement
- Large artificial-island sinks
- Generic source/sink exclusions
- Diagnostics and output writing
- Pre-relocation source suppression

## Proposed Ownership

```text
Source_sink/
|-- assemble_source_sink.py
|-- source_sink_components.py
|-- spatial_corrections.py
|-- correction_config.py
|-- source_sink_diagnostics.py
|
|-- Replace_with_USGS/
|   |-- download_usgs.py
|   |-- replace_with_obs.py
|   |-- usgs_series.py
|   |-- source_temperature.py
|   |-- source_overrides.py
|   `-- station_mappings.py
|
`-- Patch_artificial_island/
    |-- patch_artificial_island_source_sink.py
    `-- artificial_island_source_sink.yml
```

Expected responsibilities:

- `source_sink_components.py`: aligned source/sink element and time-series
  operations.
- `spatial_corrections.py`: nearest-element searches, radius exclusions, and
  region-based corrections.
- `correction_config.py`: loading and normalization for configured source/sink
  corrections.
- `source_sink_diagnostics.py`: generic source/sink diagnostic output.
- `usgs_series.py`: high-level USGS series retrieval, normalization, merging,
  and interpolation.
- `source_temperature.py`: automatic source-temperature correction.
- `source_overrides.py`: explicit Delaware/Hudson flow and temperature
  corrections.
- `station_mappings.py`: explicit flow, temperature, and NWM-feature station
  mappings.
- `patch_artificial_island_source_sink.py`: only artificial-island source/sink
  placement and corrections.

The exact file structure may be adjusted if an extraction reveals a cleaner
dependency boundary. Avoid creating small files that do not represent a clear
responsibility.

## Test Gates for Every Step

Run the checks relevant to each step before committing:

- [ ] Python syntax and import checks pass.
- [ ] New focused unit tests pass.
- [ ] Existing fast tests pass.
- [ ] No unexpected network access occurs in unit tests.
- [ ] Source element IDs and ordering are unchanged.
- [ ] Sink element IDs and ordering are unchanged.
- [ ] `vsource`, `msource`, and `vsink` shapes are unchanged.
- [ ] Time arrays are unchanged.
- [ ] Numerical values are unchanged unless the step explicitly changes
      behavior.
- [ ] No duplicate source or sink elements are introduced.
- [ ] Diagnostics remain equivalent where relevant.
- [ ] A cached representative workflow is compared when practical.

Use mocks or cached data for downloader tests. Do not make routine unit tests
depend on live USGS services.

## Phase 0: Establish Characterization Tests

Status: **Completed**

- [ ] Identify the smallest representative source/sink fixture.
- [ ] Record current source and sink element lists and ordering.
- [ ] Record representative `vsource`, `msource`, and `vsink` results.
- [ ] Test USGS station-ID normalization, including leading zeros.
- [ ] Test model-time to UTC-datetime conversion.
- [ ] Test exact USGS observation replacement.
- [ ] Test interpolation across permitted short gaps.
- [ ] Test fallback to original forcing across long gaps.
- [ ] Test fallback outside USGS observation coverage.
- [ ] Test cfs-to-m3/s conversion.
- [ ] Test Hudson date chunks for gaps and unintended overlaps.
- [ ] Test Hudson adaptive retries with a mocked downloader.
- [ ] Test merging and deduplication of downloaded series.

Deliverable: tests that characterize current behavior before functions move.

## Phase 1: Extract High-Level USGS Series Handling

Status: **Completed**

Create `Replace_with_USGS/usgs_series.py` and move the relevant functions
without redesigning their algorithms.

Candidate functions:

- [ ] `_normalize_usgs_station_id`
- [ ] `_as_utc_timestamp`
- [ ] `_model_datetimes`
- [ ] `_station_id_from_record`
- [ ] `_extract_usgs_values`
- [ ] `_download_usgs_series_standard`
- [ ] `_iter_hudson_chunks`
- [ ] `_merge_usgs_series`
- [ ] `_download_hudson_window`
- [ ] `_download_hudson_window_adaptive`
- [ ] `_download_hudson_usgs_series`
- [ ] `_interpolate_usgs_with_original_fallback`

Constraints:

- [ ] Preserve current signatures and defaults.
- [ ] Preserve cache naming for the first extraction.
- [ ] Preserve retry behavior.
- [ ] Preserve interpolation and fallback behavior.
- [ ] Temporarily re-export functions if needed to avoid changing callers in
      the same step.
- [ ] Run the Phase 0 tests.

Deliverable: USGS-series logic located under `Replace_with_USGS`, with no
intended behavior change.

## Phase 2: Centralize Station Mappings

Status: **Completed**

Create `Replace_with_USGS/station_mappings.py`.

- [ ] Move flow-station mappings.
- [ ] Move temperature-station mappings.
- [ ] Move manual NWM-feature-to-USGS mappings.
- [ ] Move USGS parameter IDs.
- [ ] Preserve Hudson's separate flow and temperature stations.
- [ ] Remove duplicated mapping definitions only after imports are working.
- [ ] Add exact-value tests for important mappings.

Important current mapping:

```python
"Hudson River": {
    "flow": "01358000",
    "temperature": "01359139",
}
```

Deliverable: a single authoritative location for station mappings.

## Phase 3: Extract Automatic Temperature Correction

Status: **Completed**

Create `Replace_with_USGS/source_temperature.py`.

- [ ] Move source-element-to-NWM-feature loading.
- [ ] Move NWM/USGS station association.
- [ ] Move temperature-station discovery.
- [ ] Move bulk temperature downloading.
- [ ] Move automatic temperature replacement.
- [ ] Rename relocation-specific terms only in a separate, reviewed change.
- [ ] Preserve the current behavior of scanning every source and modifying
      only sources with usable observations.

Tests:

- [ ] Eligible source temperature is modified.
- [ ] Unassociated source is unchanged.
- [ ] Missing observations retain the original temperature.
- [ ] Long gaps retain the original temperature.
- [ ] Source ordering and flow are unchanged.

Deliverable: general source-temperature correction independent of the
artificial-island module.

## Phase 4: Extract Delaware and Hudson Overrides

Status: **Completed**

Create `Replace_with_USGS/source_overrides.py`.

- [ ] Move the explicit existing-source replacement logic.
- [ ] Preserve coordinate matching and search radii initially.
- [ ] Preserve direct replacement and fallback behavior.
- [ ] Preserve negative-flow handling.
- [ ] Keep flow and temperature station selection parameter-specific.

Tests:

- [ ] Correct source is selected within the configured radius.
- [ ] Missing source within the radius produces the current error behavior.
- [ ] Flow-only replacement works.
- [ ] Temperature-only replacement works.
- [ ] Existing forcing is retained outside observation coverage.
- [ ] Negative-flow handling is unchanged.

Do not move these stations into `manual_nwm2usgs` during the mechanical
extraction. The existing post-patch uses direct USGS replacement, while
`source_nwm2usgs()` applies an `original + USGS - NWM` correction. Choosing
between those algorithms is a later scientific decision.

Deliverable: Delaware/Hudson logic under `Replace_with_USGS`, with the current
algorithm unchanged.

## Phase 5: Extract Generic Source/Sink Components

Status: **Completed**

Create `Source_sink/source_sink_components.py` if the extracted modules need a
shared representation for aligned arrays.

Candidate operations:

- [ ] Copy source components.
- [ ] Copy sink components.
- [ ] Build `TimeHistory` objects.
- [ ] Validate matching time arrays.
- [ ] Add and remove source columns.
- [ ] Add and remove sink columns.
- [ ] Split negative source flow into a sink.
- [ ] Rebuild a `source_sink` object.

Tests:

- [ ] Element and data-column alignment is preserved.
- [ ] All tracer columns remain aligned with sources.
- [ ] Removing the first, middle, and last columns works.
- [ ] Empty source or sink collections work.
- [ ] Duplicate element IDs are rejected where appropriate.

Deliverable: reusable, tested source/sink transformations without scientific
station or location policy.

## Phase 6: Extract Generic Spatial Corrections

Status: **Completed**

Create `Source_sink/spatial_corrections.py`.

- [ ] Move grid-center calculation.
- [ ] Move coordinate transformation.
- [ ] Move nearest-element search.
- [ ] Move radius search.
- [ ] Move region-based source zeroing.
- [ ] Move generic source/sink exclusions.

Tests:

- [ ] Nearest-element selection is correct.
- [ ] Radius boundary behavior is defined and tested.
- [ ] Sources inside a region are selected correctly.
- [ ] Sources outside a region remain unchanged.
- [ ] Region zeroing preserves source IDs and tracer columns.
- [ ] Exclusions preserve remaining array alignment.

Deliverable: spatial corrections independent of artificial-island and USGS
policy.

## Phase 7: Reduce the Artificial-Island Module

Status: **Completed**

After the preceding extractions, retain only genuinely artificial-island
behavior:

- [x] Forced Wando source placement.
- [x] Forced Turkey source placement.
- [x] Forced Buffalo Bluff source placement.
- [x] Forced Dunns Creek source placement.
- [x] Buffalo Bluff special sink.
- [x] Artificial-island-specific negative-flow handling.
- [x] Pre-relocation suppression temporarily, until relocation is decoupled.
- [x] Artificial-island-specific diagnostics.

Remove from the artificial-island module:

- [x] General USGS series downloading.
- [x] Automatic all-source temperature replacement.
- [x] Delaware/Hudson replacement.
- [x] Generic region zeroing.
- [x] Generic exclusions.
- [x] Duplicated station mappings.

Deliverable: a small artificial-island module with a defensible name and
scope.

## Phase 8: Restore Independent Workflow Controls

Status: **Completed**

This phase changes orchestration and configuration, so it must remain separate
from mechanical extraction commits.

Target structure:

```python
source_sink = generate_nwm_source_sink(...)

if config.replace_nwm_flow_with_usgs:
    source_sink = replace_flow_with_usgs(source_sink, ...)

if config.relocate_source:
    source_sink = relocate_source_sink(source_sink, ...)

if config.replace_source_temperature_with_usgs:
    source_sink = replace_temperature_with_usgs(source_sink, ...)

if config.source_overrides:
    source_sink = apply_source_overrides(source_sink, ...)

if config.zero_source_regions:
    source_sink = zero_sources_in_regions(source_sink, ...)

if config.artificial_island_corrections:
    source_sink = apply_artificial_island_corrections(source_sink, ...)
```

- [x] Each option enables one clear operation.
- [x] Temperature processing does not depend on artificial-island settings.
- [x] Delaware/Hudson processing does not depend on artificial-island
      settings.
- [x] Generic spatial corrections do not depend on artificial-island settings.
- [x] Disabled operations do not require their input files.
- [x] Defaults are documented and tested.

Current incremental state:

- [x] `assemble_source_sink.py` calls public temperature, explicit override,
      region-zeroing, and artificial-island stages in visible order.
- [x] Temperature processing selects `original_source_sink/sources.json` when
      relocation is disabled.
- [x] The artificial-island stage no longer repeats the extracted explicit
      overrides or region zeroing when called by the main workflow.
- [x] Added separate configuration switches and configuration paths for
      selected-source overrides, region zeroing, and artificial-island work.

Deliverable: independently configurable workflow stages.

## Phase 8A: Integrate Extracted Post-Relocation Stages

Status: **Completed**

Scope this phase to the extracted post-relocation stages beginning at
`apply post-generation corrections` in `assemble_source_sink.py`. Preserve the
earlier NWM-flow replacement, pre-relocation artificial-island suppression,
and relocation behavior until Phase 9. The caller order is the integration
order; do not combine the following steps into one commit.

Baseline gate for every step:

- Run focused unit tests for the stage being integrated.
- Run the fast non-MPI repository tests.
- Run the cached I102a post-relocation reproduction, which currently takes
  about two minutes.
- Compare `source_sink.in`, `vsource.th`, `msource.th`, and `vsink.th`
  byte-for-byte with HJ's reference.
- Compare every `source.nc` dimension and variable, using exact equality where
  possible and an absolute `vsource` tolerance of `1e-10` m3/s for the
  already-observed serialization roundoff.
- Confirm that no constant-sink directory is created by the focused test.

### Step 8A.1: Integrate Selected Delaware/Hudson Flow Overrides

This is the first post-relocation caller in `assemble_source_sink.py`.

- [x] Add a dedicated `apply_source_flow_overrides()` stage under
      `Replace_with_USGS`.
- [x] Split named USGS flow retrieval from temperature retrieval so the flow
      stage does not download or modify temperature.
- [x] Call selected-source flow replacement before temperature processing in
      `assemble_source_sink.py`.
- [x] Preserve direct replacement with original relocated forcing in data
      gaps; do not substitute the different pre-relocation
      USGS-minus-NWM-difference algorithm from `source_nwm2usgs()`.
- [x] Verify that only relocated elements 2 (Hudson) and 948 (Delaware) change
      and both temperature tracers remain exact.
- [x] Pass the cached I102a reproduction with byte-identical text forcing and
      the established `source.nc` tolerance.

Deliverable: Delaware/Hudson flow correction is an independent
`Replace_with_USGS` stage.

### Step 8A.2: Integrate Automatic Source Temperature Replacement

This is the second post-relocation caller in `assemble_source_sink.py`.

- [x] Make `source_temperature.py` the sole owner of automatic all-source
      temperature replacement.
- [x] Keep the direct call from `assemble_source_sink.py` and its independent
      configuration switch.
- [x] Make the temperature stage create and own its diagnostics directory;
      prevent `invalid_usgs_station_coordinates.csv` or other diagnostics
      from leaking into the process working directory.
- [x] Remove `replace_all_source_temperature`, the automatic-temperature
      branch, and the corresponding private import from
      `patch_artificial_island_source_sink()`.
- [x] Verify that this operation changes only temperature tracer 1 and retains
      source flow, sinks, element IDs, column order, and time arrays.
- [x] Reproduce the current count of 25 temperature-adjusted I102a sources.

Deliverable: temperature replacement has one production implementation and
one caller-visible workflow stage.

### Step 8A.3: Integrate Selected Delaware/Hudson Temperature Overrides

This is the third post-relocation caller. Its explicit station values override
the preceding automatic temperature result for Delaware and Hudson.

- [x] Make `source_overrides.py` own loading and normalizing its selected-source
      configuration, while retaining the array-level helper for focused tests.
- [x] Keep the direct, independently controlled call from
      `assemble_source_sink.py`.
- [x] Remove compatibility imports and parameters that are no longer used by
      the artificial-island module.
- [x] Remove `replace_only_source_locations` handling and its generic loop from
      `patch_artificial_island_source_sink()` after both flow and temperature
      have independent production callers.
- [x] Verify that only temperature tracer 1 changes for Delaware and Hudson;
      source flow, sinks, source IDs, and ordering remain unchanged.
- [x] Verify the cached temperature fallback behavior for the known Hudson
      data gap.

Deliverable: Delaware/Hudson temperature correction is independent of its
flow correction and has no artificial-island fallback path.

### Step 8A.4: Integrate Configured Region Zeroing

This is the fourth post-relocation caller.

- [x] Make the spatial-correction stage own loading and resolving its region
      configuration, while retaining the polygon/array helper for focused
      tests.
- [x] Keep the direct, independently controlled call from
      `assemble_source_sink.py`.
- [x] Remove `zero_source_regions` parsing and execution from
      `patch_artificial_island_source_sink()`.
- [x] Verify that only `vsource` values change; retain source elements,
      tracers, sinks, ordering, and time arrays.
- [x] Reproduce the current I102a result of three zeroed Savannah-region
      sources.

Deliverable: generic region zeroing is owned completely by the spatial
correction module.

### Step 8A.5: Integrate the Reduced Artificial-Island Stage

This is the fifth and final post-relocation correction caller.

- [x] Rename or expose the production entry point as
      `apply_artificial_island_corrections()`; retain a temporary compatibility
      alias only if an external caller requires it.
- [x] Reduce its interface after Steps 8A.1-8A.3: remove the relocated mapping,
      NWM shapefile, state list, and generic-stage mode flags that are no longer
      needed.
- [x] Retain only forced Wando, Turkey, Buffalo Bluff, and Dunns Creek source
      handling, Buffalo Bluff's large sink, island-specific negative-flow
      handling, forcing-template access, and island diagnostics.
- [x] Characterize any configured island exclusion before deciding whether to
      delete it or move it to the generic spatial stage; do not silently retain
      a generic responsibility in this module.
- [x] Reproduce the current I102a final count of 1,343 sources and 3 sinks.

Deliverable: the artificial-island stage has a narrow name, interface, and
responsibility.

### Step 8A.6: Clean Up the Caller

- [x] Remove obsolete compatibility variables, imports, comments, and YAML
      parsing from `assemble_source_sink.py`.
- [x] Keep the five stage calls visibly ordered as selected-source flow,
      automatic temperature, selected-source temperature, region zeroing, and
      artificial-island corrections.
- [x] Resolve the USGS cache path once and only when a USGS-backed stage is
      enabled.
- [x] Verify independently that each stage can be disabled without requiring
      its YAML, cache, mapping, or diagnostics path.
- [x] Run the complete focused reproduction once more and record timings and
      comparisons in the progress log.

Deliverable: `assemble_source_sink.py` is a readable workflow caller, while
scientific behavior and configuration ownership reside in the stage modules.

## Phase 9: Decouple Artificial-Island Corrections from Relocation

Status: **Not started**

Clarify intended scientific behavior before implementation:

- [ ] For each forced source, define whether the existing source is moved,
      replaced, retained, or zeroed.
- [ ] Define the forcing fallback when USGS data are missing.
- [ ] Define whether zero-flow placeholders should remain.
- [ ] Define behavior when the target element already contains a source.
- [ ] Define behavior with and without relocation.
- [ ] Confirm which grid supplies source element IDs.

Implementation goals:

- [ ] Operate on `base_ss` regardless of relocation state.
- [ ] Do not require `relocated_source_sink/sources.json` when relocation is
      disabled.
- [ ] Preserve an unmodified forcing template before replacing a source.
- [ ] Avoid duplicate source elements.
- [ ] Remove the pre/post split when relocation is disabled.

Deliverable: artificial-island corrections that do not silently require
relocation.

## Phase 10: Evaluate Generalizing the Hudson Downloader

Status: **Completed**

Only begin after the current Hudson behavior is covered by tests and moved
under `Replace_with_USGS`.

- [x] Compare Hudson retry behavior with `download_stations()`.
- [x] Determine whether `days_per_chunk` and existing fallback are sufficient.
- [x] If necessary, add a configurable retry schedule to the generic layer.
- [x] Verify identical merged series for cached representative data.
- [x] Remove Hudson-specific downloader functions only after equivalence is
      established.

Possible future interface:

```python
download_usgs_series(
    station_id,
    parameter_id,
    start_time,
    end_time,
    policy=UsgsDownloadPolicy(
        primary_chunk_days=100,
        retry_chunk_days=(20, 10, 5),
    ),
)
```

Deliverable: Hudson-specific behavior represented as downloader configuration
if equivalence can be demonstrated.

## Phase 11: Make Source/Sink Stage Outputs Explicit

Status: **Deferred until all behavior-preserving reproduction tests pass**

The name `original_source_sink` is currently misleading. It represents the
current pre-relocation source/sink state, rather than an immutable copy of the
original NWM result. The workflow recreates the directory at startup, replaces
its `vsource.th` during NWM-to-USGS adjustment, and may subsequently modify the
linked adjusted flow during pre-relocation zeroing.

Treat this as the final cleanup step so that changing file ownership and data
lineage does not complicate the reproduction test or the earlier scientific
decoupling work.

- [ ] Preserve the initial NWM result as an immutable, clearly named stage
      output.
- [ ] Give USGS flow adjustment, pre-relocation corrections, relocation, and
      post-generation corrections separate explicit outputs.
- [ ] Stop replacing `original_source_sink/vsource.th` with a symlink to a
      mutable downstream result.
- [ ] Make each stage consume the preceding stage explicitly rather than
      relying on the evolving meaning of `original_source_sink`.
- [ ] Define safe rerun behavior so that an existing stage output is not
      silently removed.
- [ ] Update diagnostics, documentation, configuration, and tests to reflect
      the new stage names and data lineage.
- [ ] Verify final source/sink products against the behavior-preserving
      reference run before removing compatibility paths.

Deliverable: explicit, independently inspectable stage directories with an
immutable original NWM source/sink result.

## Collaboration and Review

- Structural reorganization and workflow interfaces are owned on this branch.
- Scientific intent should be confirmed for the artificial-island locations.
- Avoid open-ended simultaneous edits to
  `patch_artificial_island_source_sink.py`.
- Review one responsibility per commit or small pull request.
- Record intentional behavior changes in the decision log below.

## Scientific Behavior Questions

Complete before Phase 9:

| Location | Existing source action | Flow source | Temperature source | Sink handling |
|---|---|---|---|---|
| Wando | TBD | Current NWM template | Current/default | None currently |
| Turkey | TBD | USGS with current fallback | TBD | None currently |
| Buffalo Bluff | TBD | USGS with current fallback | USGS | Negative flow and constant sink |
| Dunns Creek | TBD | USGS with current fallback | USGS | Negative flow |
| Clyo | Zero or remove: TBD | N/A | Unchanged currently | N/A |
| Moultrie | Zero or remove: TBD | N/A | Unchanged currently | N/A |

## Decision Log

| Date | Decision | Reason |
|---|---|---|
| 2026-09-30 | Use `refactor/source-sink-workflow` for the work. | Isolate incremental refactoring from `master`. |
| 2026-09-30 | Separate mechanical extraction from behavior changes. | Make regressions easier to identify and review. |
| 2026-09-30 | Test after every extraction. | Preserve scientific behavior and array alignment. |
| 2026-09-30 | Do not immediately merge Hudson into the generic downloader. | Current adaptive retry behavior is not yet proven equivalent. |
| 2026-09-30 | Do not immediately move Delaware/Hudson into `manual_nwm2usgs`. | The current and proposed flow algorithms differ. |
| 2026-09-30 | Give temperature, selected-source overrides, region zeroing, and artificial-island work independent configuration controls. | Removing the artificial-island gate restores the original workflow's stage-level control. |
| 2026-09-30 | Split the shared correction YAML by stage ownership. | Each operation should load only the configuration and supporting files it owns. |
| 2026-09-30 | Integrate extracted post-relocation stages in their caller order and remove the corresponding compatibility branch after each stage passes the cached I102a gate. | Small commits isolate numerical regressions and progressively make each module the sole owner of its responsibility. |
| 2026-09-30 | Integrate Delaware/Hudson flow before temperature, but keep it post-relocation for the reference reproduction. | Direct relocated-source replacement and NWM gap fallback reproduce the reference; the pre-relocation `source_nwm2usgs()` algorithm is not byte-equivalent. |
| 2026-09-30 | Use the soft pre-relocation NWM-to-USGS associations for Delaware and Hudson in the eventual simplified workflow. | The I102a comparison showed nearly identical Delaware flow and a defensible Hudson difference because the soft method retains the other relocated FeatureID contributions. |
| 2026-09-30 | Do not enable the Delaware/Hudson soft links while reproducing HJ's I102a; retain the post-relocation overrides for now. | A post-relocation override replaces soft flow where USGS is supported but retains the soft result during long gaps, so enabling both would not be byte-identical to HJ. Add the two manual mappings and disable their post-relocation flow overrides together in a later transition. |
| 2026-10-01 | Keep `first_usable` as the v7.4 replication default and add `discharge_weighted` as an opt-in temperature mode. | This preserves HJ's output while allowing multiple upstream temperature stations to be combined using discharge from each station's associated NWM feature. |
| 2026-10-01 | Apply no spatial discharge-coverage threshold, but require a complete finite time series before writing a weighted temperature column. | A lone small-creek temperature may still represent nearby water temperature; atomic replacement prevents SCHISM from interpolating unpredictably between numeric values and the `-9999` ambient-temperature sentinel. |
| 2026-10-01 | Keep generic region zeroing immediately before artificial-island corrections. | Zeroing remains an optional manual enforcement step, while the following island stage can restore required island sources so they are not erased accidentally. |
| 2026-10-01 | Remove post-generation island exclusion handling. | The v7.4 island exclusion list was empty, so retaining the generic removal machinery added responsibility without affecting output. Non-empty generic spatial sections now fail explicitly in the island stage. |
| 2026-10-01 | Represent adaptive USGS retrieval as source configuration rather than a Hudson code branch. | A generic primary/retry-window policy preserves the cached Hudson behavior and can be assigned to any configured source. |
| 2026-10-01 | Place correction parsing and diagnostics with their owning stages, but retain shared SCHISM component operations at the source/sink level. | This removes mixed-responsibility modules without creating a generic `utils` bucket or changing numerical operations. |
| 2026-10-01 | Store manual NWM FeatureID-to-USGS links in stage-local YAML with separate flow-adjustment and temperature-search scopes. | The temperature search has three additional links; preserving that distinction keeps the I102a flow workflow unchanged. |
| 2026-10-01 | Disable selected Delaware/Hudson post-relocation overrides in the v7.4 preferred profile and apply their FeatureID links during NWM-to-USGS flow correction instead. | This retains the preferred soft association workflow after the HJ byte-identity milestone while keeping the HJ replay inputs available for reproduction. |

## Progress Log

| Date | Phase | Work completed | Tests/results | Commit |
|---|---|---|---|---|
| 2026-09-30 | Setup | Created refactor branch and tracking plan. | Pending | Pending |
| 2026-09-30 | 0-7 | Added characterization tests and extracted USGS series, station mappings, automatic temperature processing, Delaware/Hudson overrides, aligned source/sink components, spatial helpers, configuration, diagnostics, and background-sink overlap handling. | 22 focused tests and 87 repository tests passed; compile checks passed; 47 extracted function bodies match `master` by AST comparison. | `08252989` |
| 2026-09-30 | 8 | Added public stage APIs and made the main assembler call temperature replacement, Delaware/Hudson overrides, region zeroing, and artificial-island corrections explicitly. | 26 focused tests and 91 non-MPI repository tests passed; compile and diff checks passed. | `d545bd60` |
| 2026-09-30 | 8 | Moved the extracted stages outside the artificial-island conditional and added independent configuration switches with v7.4 compatibility settings. | 28 focused tests and 93 non-MPI repository tests passed; compile and diff checks passed. | `8503e7b3` |
| 2026-09-30 | 8 | Split Delaware/Hudson overrides, Savannah region zeroing, and artificial-island settings into separate YAML files and moved the region asset beside its configuration. | YAML values and region bytes match the previous combined configuration; 28 focused tests and 93 non-MPI repository tests passed. | `8c47653a` |
| 2026-09-30 | 8 | Reproduced the post-relocation stages from HJ's cached I102a relocated forcing, stopping before constant-sink generation. | Completed in 122.5 seconds. `source_sink.in`, `vsource.th`, `msource.th`, and `vsink.th` are byte-identical to the reference. NetCDF dimensions and variables match; only `vsource` has floating-point roundoff, with maximum absolute difference `1.46e-11` m3/s. | Pending |
| 2026-09-30 | 8A.1 | Split selected Delaware/Hudson flow from temperature retrieval and made flow the first post-relocation correction stage. | 29 focused tests and 94 non-MPI repository tests passed. The cached I102a chain completed in 119.5 seconds with byte-identical text forcing and maximum NetCDF `vsource` difference `1.46e-11` m3/s. | Pending |
| 2026-09-30 | 8A investigation | Compared manual pre-relocation NWM-to-USGS associations for Hudson (`6186156 -> 01358000`) and Delaware (`2590217 -> 01463500`) with HJ's direct post-relocation replacements using pinned cached observations and the same relocation mapping. | Both manual corrections ran successfully across 9,505 records. During records where HJ's direct USGS replacement was active, the soft method differed by a mean/MAE of `15.54/15.54` m3/s for Hudson and `1.25/1.26` m3/s for Delaware. The retained contributions from the other relocated feature IDs explain the predominantly positive differences. Full-period discrepancies are larger around direct-replacement fallback gaps. Results are under `manual_nwm2usgs_delaware_hudson/` in the I102a test area. | Pending |
| 2026-10-01 | 8A.2-8A.3 | Made automatic temperature and selected Delaware/Hudson temperature overrides independent production stages, removed both compatibility paths from the artificial-island module, and added opt-in discharge-weighted pooling with atomic all-numeric/all-ambient columns. | 100 non-MPI tests passed. Cached I102a legacy reproduction completed in 118.7 seconds with 25 automatic replacements, 2 selected temperature overrides, 3 region-zeroed sources, and 1,343 final sources/3 sinks. All text forcing files are byte-identical; NetCDF dimensions/variables match and the maximum `vsource` difference is `1.46e-11` m3/s. Output: `I102a_temperature_integration/`. | Pending |
| 2026-10-01 | 8A.4-8A.6 | Gave the spatial stage ownership of region configuration, reduced the island entry point to island responsibilities, removed the empty exclusion path, and cleaned the five-stage caller. | 100 non-MPI tests passed. Final cached I102a gate completed in 118.1 seconds with three region-zeroed sources and 1,343 final sources/3 sinks. All text forcing files are byte-identical; NetCDF dimensions and variables match with maximum `vsource` difference `1.46e-11` m3/s. No constant-sink directory was created. Output: `I102a_final_integration/`. | Pending |
| 2026-10-01 | 10 | Replaced the Hudson-only downloader branch and helpers with a generic configurable USGS download policy. The Hudson override selects `100 -> 20 -> 10 -> 5` day windows in YAML while retaining its existing cache prefix. | Cached I102a gate completed in 118.5 seconds with unchanged operation counts. All text forcing files and the prior integration `source.nc` are byte-identical. Output: `I102a_generic_usgs_downloader/`. | Pending |
| 2026-10-01 | Milestone | Reproduced I102a from the completed initial NWM/USGS stage through the island post-patch, before constant sinks. | Pre-zeroed flow, relocation mapping/forcing, and post-patch text are byte-identical to HJ's reference; every NetCDF variable matches exactly. Output: `I102a_full_chain_milestone/`. | `7cd9b818` |
| 2026-10-01 | Final module sweep | Split stage-specific configuration, shared spatial selection, and region zeroing; moved island diagnostics and feeder patch to their owning folders; moved manual FeatureID-to-USGS links to YAML. Kept `source_sink_components.py` shared. | 39 source/sink tests pass. Full cached I102a reproduction completed in 338.5 seconds with byte-identical text forcing and exact NetCDF variable values; no constant-sink directory. Scratch output: `I102a_structure_cleanup_milestone/`. | `d6dfec35`, `9efeaba3`, `579f0209`, `f1daa581`, `2f0c82a7` |
