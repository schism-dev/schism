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

Status: **In progress**

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

Status: **In progress**

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

Status: **In progress**

After the preceding extractions, retain only genuinely artificial-island
behavior:

- [ ] Forced Wando source placement.
- [ ] Forced Turkey source placement.
- [ ] Forced Buffalo Bluff source placement.
- [ ] Forced Dunns Creek source placement.
- [ ] Buffalo Bluff special sink.
- [ ] Artificial-island-specific negative-flow handling.
- [ ] Pre-relocation suppression temporarily, until relocation is decoupled.
- [ ] Artificial-island-specific diagnostics.

Remove from the artificial-island module:

- [ ] General USGS series downloading.
- [ ] Automatic all-source temperature replacement.
- [ ] Delaware/Hudson replacement.
- [ ] Generic region zeroing.
- [ ] Generic exclusions.
- [ ] Duplicated station mappings.

Deliverable: a small artificial-island module with a defensible name and
scope.

## Phase 8: Restore Independent Workflow Controls

Status: **Not started**

This phase changes orchestration and configuration, so it must remain separate
from mechanical extraction commits.

Target structure:

```python
source_sink = generate_nwm_source_sink(...)

if config.replace_nwm_flow_with_usgs:
    source_sink = replace_flow_with_usgs(source_sink, ...)

if config.replace_source_temperature_with_usgs:
    source_sink = replace_temperature_with_usgs(source_sink, ...)

if config.source_overrides:
    source_sink = apply_source_overrides(source_sink, ...)

if config.relocate_source:
    source_sink = relocate_source_sink(source_sink, ...)

if config.zero_source_regions:
    source_sink = zero_sources_in_regions(source_sink, ...)

if config.artificial_island_corrections:
    source_sink = apply_artificial_island_corrections(source_sink, ...)
```

- [ ] Each option enables one clear operation.
- [ ] Temperature processing does not depend on artificial-island settings.
- [ ] Delaware/Hudson processing does not depend on artificial-island
      settings.
- [ ] Generic spatial corrections do not depend on artificial-island settings.
- [ ] Disabled operations do not require their input files.
- [ ] Defaults are documented and tested.

Deliverable: independently configurable workflow stages.

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

Status: **Not started**

Only begin after the current Hudson behavior is covered by tests and moved
under `Replace_with_USGS`.

- [ ] Compare Hudson retry behavior with `download_stations()`.
- [ ] Determine whether `days_per_chunk` and existing fallback are sufficient.
- [ ] If necessary, add a configurable retry schedule to the generic layer.
- [ ] Verify identical merged series for cached representative data.
- [ ] Remove Hudson-specific downloader functions only after equivalence is
      established.

Possible future interface:

```python
download_station_series(
    station_id,
    parameter_id,
    start_time,
    end_time,
    primary_chunk_days=100,
    retry_chunk_days=(20, 10, 5),
)
```

Deliverable: Hudson-specific behavior represented as downloader configuration
if equivalence can be demonstrated.

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

## Progress Log

| Date | Phase | Work completed | Tests/results | Commit |
|---|---|---|---|---|
| 2026-09-30 | Setup | Created refactor branch and tracking plan. | Pending | Pending |
| 2026-09-30 | 0-7 | Added characterization tests and extracted USGS series, station mappings, automatic temperature processing, Delaware/Hudson overrides, aligned source/sink components, spatial helpers, configuration, diagnostics, and background-sink overlap handling. | 22 focused tests and 87 repository tests passed; compile checks passed; 47 extracted function bodies match `master` by AST comparison. | `08252989` |
