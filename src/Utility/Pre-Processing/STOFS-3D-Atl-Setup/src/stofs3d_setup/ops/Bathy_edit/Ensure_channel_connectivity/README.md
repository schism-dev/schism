# Hydrofabric matching prototype

`hydrofabric_match.py` tests partial-overlap matching between RiverMapper river
centerlines and hydrofabric reaches. It does not change the production
channel-connectivity bathymetry workflow.

The default invocation uses the Savannah-area test bounds and input files:

```bash
python src/stofs3d_setup/ops/Bathy_edit/Ensure_channel_connectivity/hydrofabric_match.py
```

By default, generated products are written outside the Git worktree to:

```text
/sciclone/schism10/feiye/CIROH/Channel_Geometry/Connectivity_test/savannah
```

Use one subdirectory under `Connectivity_test` for each additional test region
so its caches and diagnostics do not overwrite another region.

The cache directory contains:

- `hydrofabric.parquet`: regional hydrofabric reaches in EPSG:5070.
- `river_arcs.parquet`: regional RiverMapper arcs.
- `river_centerlines.parquet`: centerlines assembled from those arcs.
- `hgrid_nodes.npz`: regional mesh-node IDs, coordinates, and depths.
- `hydrofabric_river_matches.gpkg`: GIS diagnostic layers for candidate
  overlaps, selected COMID intervals, station assignments, intervals flagged
  for review, unmatched stations, and separate unmatched-cause layers.
- `match_overview.png`: overview grouped by interval review flag and unmatched
  station reason.
- `summary.json`: counts and regional match coverage.

The cache metadata includes input fingerprints, bounds, and buffer distance.
Use `--force-cache` to rebuild all regional inputs explicitly. Run `--help` for
the distance, direction, sampling, and overlap controls.

The current matching prototype resolves overlapping candidates from local
distance and direction. Network-topology continuity and full per-river
channel-polygon expansion across ordinary reaches are intentionally not
implemented. The later dredging stage does use a channel polygon to screen
reverse-search candidates specifically at multi-river junctions.

The acceptance limits are 75 degrees, 25 m contiguous overlap, or 10%
coverage of the shorter line. Review flags retain the stricter 30-degree,
100 m, and 250 m thresholds so relaxed matches remain easy to inspect.

After line matching, a continuity fallback fills a bounded unmatched station
run only when its nearest matched stations on both sides have the same COMID.
It also bridges gaps no longer than 25 m when the selected intervals directly
on both sides have the same COMID. These additions are labeled
`continuity_fill`, record their `continuity_basis`, and are written to dedicated
point and interval layers in the diagnostic GeoPackage.

## Full-domain MPI run

`hydrofabric_match_mpi.py` applies the same matching and continuity functions
to the complete RiverMapper map. A 20-rank full-domain launch is:

```bash
srun -n 20 python \
  src/stofs3d_setup/ops/Bathy_edit/Ensure_channel_connectivity/hydrofabric_match_mpi.py \
  --overwrite
```

The default full-domain output directory is:

```text
/sciclone/schism10/feiye/CIROH/Channel_Geometry/Connectivity_test/full_domain
```

Rank 0 builds two reusable inputs under `full_domain/input_cache`:

- `river_arcs.parquet` and `river_centerlines.parquet`, parsed once from the
  complete RiverMapper SMS map.
- `hydrofabric_indexed.gpkg`, a valid-`bnk_depth`, EPSG:5070 spatial cache made
  with a streamed GDAL conversion. Each rank then reads only the reaches in a
  500 m-expanded bounding box around its current work tile.

Complete `river_idx` groups are assigned to 100 km tiles. Tiles are greedily
balanced among ranks by RiverMapper station count, so a river and all its
stations always remain on one rank and the continuity fallback cannot be split
at an MPI boundary. Each rank writes its own temporary GeoPackage; only rank 0
merges the final `hydrofabric_river_matches.gpkg` and writes `summary.json`.
Temporary rank files are removed after a successful merge unless
`--keep-parts` is specified.

The indexed hydrofabric remains beside the final diagnostic GeoPackage rather
than being duplicated into it. Use `--include-hydrofabric-layer` if a single,
larger standalone GeoPackage is required. Use `--force-input-cache` only when
the source data changed or a cache rebuild is explicitly needed.

For a bounded MPI regression, add `--bbox MIN_LON MIN_LAT MAX_LON MAX_LAT` and
use a separate `--output-dir`. The hgrid is intentionally not read by this
matching stage: it produces the vertex-to-COMID/`bnk_depth` mapping consumed by
the serial and MPI mesh-assimilation drivers documented below.

## Hydrofabric-depth KD-tree dredging

`hydrofabric_dredge.py` is the first mesh-assimilation implementation. It keeps
the existing nearest-node approach while using `selected_intervals` as the
authoritative piecewise COMID and bankfull-depth mapping. The established
scalar-depth production function remains unchanged while this path is tested.

For each complete RiverMapper river, the driver:

1. Materializes the selected intervals at centerline stations using half-open
   interval ownership. Primary line matches take precedence; continuity-fill
   intervals assign only stations that remain unmatched.
2. Propagates each station's COMID and `bnk_depth` across all parallel arcs.
3. Reconstructs the production connectivity mask from `watershed.shp` minus
   the configured Maine exclusion polygon. The optional test bounding box is
   intersected with this mask; it does not replace it.
4. Builds one EPSG:5070 `cKDTree` from hgrid nodes inside the effective
   watershed and maps every arc vertex to its nearest eligible mesh node. This
   prevents an interior RiverMapper vertex near the mask boundary from
   requesting an outside-watershed mesh node.
5. Computes the existing high-bank or lower-bank reference at each station.
6. Uses `max(min_channel_depth, bnk_depth)` for matched stations and, by
   default, the scalar minimum depth for unmatched stations.
7. Requests changes only from inner arc vertices inside the effective
   watershed and skips mappings farther than 500 m.
8. Reduces duplicate vertex requests with `numpy.maximum.at`, ensuring the
   deepest request wins and no mesh node becomes shallower.
9. Runs reverse screening to recover mesh nodes around intersections. The 200 m search
   first finds mesh-node candidates near at least two distinct `river_idx`
   values. For each contributing river, the local centerline distance is then
   compared with half the maximum RiverMapper bank-to-bank width across the
   nearest station and its immediate neighbors. A final target requires at
   least two rivers to pass `distance <= half width + 10 m`. Each river must
   also lie farther than 5% of its local width from both continuous bank arcs.
   The junction is retained when at least one half-width-qualified river passes
   this bank-line test, allowing a tributary to cross another river's bank.
   Only bank-safe rivers contribute the proposed target depth. Bank-vertex
   distance is retained as an audit field. Bank-line distance is signed using
   the complete polygon between the two bank arcs: positive inside the channel
   and negative outside. The broad search is watershed constrained, but
   centerline, bank-line, and channel-polygon calculations use complete river
   geometries so they are not truncated at the watershed boundary. Strict
   half-width and pre-bank-screen results remain available for comparison.
   Accepted intersection targets are merged with the forward vertex requests;
   the deepest request wins when both paths select the same mesh node.

The 500 m bank-mapping limit is also an intentional scope guard. At some
watershed boundaries, a geometrically valid intersection target belongs to a
channel so wide that only one of its two RiverMapper banks is within 500 m of
an eligible in-watershed mesh node. Such a station fails the required two-bank
reference check. The point remains in `intersection_target_mesh_nodes` for
diagnosis with `target_dp_available = false`, but it contributes no dredging
request. These are no longer small watershed rivers, so leaving them untouched
is the desired behavior rather than relaxing the distance limit.

The driver is diagnostic-only unless `--write-hgrid` is given. A Savannah test
run is:

```bash
python \
  src/stofs3d_setup/ops/Bathy_edit/Ensure_channel_connectivity/hydrofabric_dredge.py \
  --bbox -81.37073 31.77 -81.0104 32.10649286 \
  --output-dir /path/to/dredge_tests/savannah \
  --write-gpkg
```

Each test directory contains:

- `river_vertex_mesh_mapping.parquet`: full vertex, interval, KD-tree, bank,
  target, and final-node diagnostics.
- `dredge_requested_mesh_nodes.parquet`: unique mesh nodes requested by forward
  RiverMapper vertices, accepted intersection targets, or both, before checking
  whether the requested depth is deeper than the existing mesh. The raw
  `requested_dp` and `requested_dredging_delta_m` remain visible even when a
  request is limited by the maximum-delta safety cap.
- `dredged_mesh_nodes.parquet`: one row per mesh node that deepens, including
  requests whose applied change is capped.
- `intersection_candidates_200m.parquet`: all broad multi-river candidates,
  including per-river centerline distances and half-width results.
- `intersection_target_mesh_nodes.parquet`: candidates for which at least two
  rivers pass the tolerant half-width screen and at least one of those rivers
  passes the signed channel-polygon/bank-proximity screen. It also records
  whether each node is new relative to the forward nearest-neighbor requests
  and the depth proposed by the bank-safe contributing rivers.
- `hydrofabric_dredge_diagnostics.gpkg`: request vertices, combined
  `dredge_requested_mesh_nodes`, post-depth-screen `changed_mesh_nodes`,
  intersection targets, true in-watershed distance flags, outside-watershed
  vertices, the effective watershed, and the test boundary.
- `summary.json`: mapping-distance, bankfull-depth, collision, forward-only
  request, combined forward/intersection request, and bathymetry-change
  statistics. `unique_requested_mesh_nodes` is the legacy forward-only count;
  `dredge_requested_mesh_nodes` is the combined unique-node count.

By default, planned deepening is limited to 6 m per mesh node:

```text
applied_dredging_delta_m = min(requested_dredging_delta_m, 6 m)
final_dp = original_dp + applied_dredging_delta_m
```

The 6 m value is the rounded, slightly conservative form of the Test 1 P99
planned-deepening value (6.156 m). Unlike the former 10 m rejection rule, the
cap does not leave an unchanged node between dredged neighbors. Requests above
the cap remain in `dredged_mesh_nodes` with
`passes_max_dredging_delta = false`, `dredging_delta_capped = true`, and
`depth_screen_status = capped_at_max_dredging_delta`. Their uncapped value
remains in `requested_dredging_delta_m`, while `dredging_delta_m` records the
6 m applied change. When a GPKG is requested, capped nodes are written to
`capped_dredging_delta`; the raw greater-than-10 m subset is also written to
`capped_raw_delta_gt10m` for severity review.

This choice intentionally allows major rivers to receive up to 6 m of
deepening at a small number of nodes. That is accepted as preferable to an
isolated undredged node that could interrupt connectivity. The threshold is
configurable with `--max-dredging-delta-m` and is fixed at 6 m for
reproducibility rather than recomputed on every run.

### Full-domain review layers

The compact review package is
`dredge_full_v32e/full_domain_review.gpkg`. Its layers describe independent
safety checks; membership in a review layer does not by itself mean that a
depth will be written:

- `invalid_distance_vertices`: RiverMapper vertices inside the effective
  dredging region whose nearest *eligible in-watershed* mesh node is farther
  than the 500 m mapping limit. They cannot make a forward dredging request.
  Near the watershed boundary, a much closer unrestricted mesh node can lie
  only slightly outside the polygon; the in-watershed KD-tree deliberately
  excludes it. These boundary-sensitive cases remain rejected rather than
  weakening the watershed constraint or distance guard.
- `invalid_bank_stations`: unique river stations with fewer than two valid
  bank mappings. A valid bank must map to an eligible mesh node within 500 m.
  Without both bank elevations, the station has no reliable high-bank or
  low-bank reference, so its inner RiverMapper vertices do not request
  dredging. This commonly identifies very wide channels and rivers crossing
  the effective-watershed boundary.
- `large_delta_gt10m`: requested mesh nodes whose proposed deepening exceeds
  10 m. These usually occur on major rivers with unusually large hydrofabric
  bankfull depths. This layer records the uncapped request; under the current
  policy its applied change would be limited to 6 m rather than rejected.
- `extreme_delta_gt20m`: the greater-than-20 m subset of
  `large_delta_gt10m`. It is a severity-focused view, not a separate decision
  rule; every member receives the same 6 m maximum applied change.
- `intersection_target_no_depth`: mesh nodes that pass the multi-river,
  half-width, and bank-line intersection geometry screens but have no finite
  depth from a bank-safe contributing river. In the reviewed full-domain run,
  this occurs because those contributing stations fail the two-valid-bank
  requirement, generally where a wide or boundary channel has one bank more
  than 500 m from an eligible mesh node. The target remains visible for
  geometry review but contributes no intersection dredging request. If the
  same node independently has a valid forward request, that request is still
  evaluated normally.

The reviewed snapshot contains 3 invalid-distance vertices, 2,560 invalid
bank stations, 158 proposed changes greater than 10 m (including 3 greater
than 20 m), and 125 intersection targets without a finite station depth.

`--unmatched-policy skip` disables the scalar fallback. Use
`--measured-from-lower-bank` to select the more aggressive legacy bank rule.
Use `--channel-depth-source constant` to ignore hydrofabric `bnk_depth` and
apply `--min-channel-depth-m` at every valid station. This provides a direct
scalar-depth comparison through the modern mapping and safety pipeline. Add
`--disable-intersection-recovery` for forward RiverMapper-vertex requests
only, which is the closest safeguarded analogue of the original
`ensure_channel_connectivity` targeting method.
The 500 m nearest-distance limit is configurable, but should not be relaxed
without reviewing `flagged_nearest_distance`; a distant nearest node can belong
to an unrelated part of the mesh. Relaxing it would also bring very wide
boundary channels into a workflow intended for smaller watershed rivers.

### v32e comparison base

The default hgrid for both the serial and MPI dredging drivers is the
bathymetry-loaded/edited, pre-connectivity grid from the verified v32e channel
comparison set:

```text
/sciclone/schism10/Hgrid_projects/STOFS3D-v7.4/v32e/Bathy_edit_channel_variations/Channel_variations_verified_pre_channel/pre_channel/hgrid_pre_channel.gr3
```

Using this grid isolates the hydrofabric-depth method from the later channel
connectivity edits and makes its result directly comparable with the existing
v32e channel variants.

### Full-domain MPI dredging

`hydrofabric_dredge_mpi.py` runs the same forward KD-tree mapping and reverse
intersection screening over the full domain. A 32-rank diagnostic launch on an
allocated compute node is:

```bash
srun -n 32 python \
  src/stofs3d_setup/ops/Bathy_edit/Ensure_channel_connectivity/hydrofabric_dredge_mpi.py \
  --output-dir /sciclone/schism10/feiye/CIROH/Channel_Geometry/Connectivity_test/dredge_full_v32e \
  --tile-size-m 50000 \
  --overwrite
```

The first run builds reusable memory-mapped hgrid arrays and a Parquet copy of
`selected_intervals` under `OUTPUT_DIR/input_cache`. To retain the same input
cache while testing another output directory, pass
`--input-cache-dir /path/to/existing/input_cache`. Use
`--force-input-cache` only after an input file changes.

Mesh tiles have unique spatial owners for intersection candidates. Every
owner receives complete river geometry intersecting its tile plus a 500 m
halo, so bank lines and channel polygons are never clipped at rank boundaries.
Rank 0 broadcasts the expected centerline, arc, and vertex counts for every
touched `river_idx`; each rank validates its loaded geometry against that
inventory and aborts rather than screening with a clipped river.
Complete rivers also receive one separate balanced owner for forward requests
and vertex diagnostics. Rank 0 deterministically reduces duplicate mesh-node
requests, retains the maximum requested depth, and writes the same five
Parquet products as the serial driver.

The default run computes complete forward and intersection node changes but
does not write a modified hgrid. Review `summary.json` and the Parquet products
first; add `--write-gpkg` for compact QGIS layers and add `--write-hgrid` only
when the reviewed dredged hgrid should be produced.
Per-rank temporary Parquet files are removed after a successful merge unless
`--keep-parts` is supplied. An interrupted job can be rerun with `--overwrite`;
the validated input cache is reused.

### v32e lower-bank comparison runs

The full-domain lower-bank comparison is under
`Connectivity_test/dredge_full_runs_v32e`. It uses a 2 m minimum/fallback and
the 6 m applied-deepening cap:

- `01_ml_lowerbank`: hydrofabric depths, with 2 m used as both the matched
  minimum and unmatched fallback; forward and intersection recovery enabled.
- `02_constant2_lowerbank`: constant 2 m depths at every valid station;
  forward and intersection recovery enabled.
- `03_constant2_lowerbank_forward`: constant 2 m depths with intersection
  recovery disabled. This is compared with the original
  `Channel_variations_verified_pre_channel/lower_bank_2m` product.

The regenerated 32-rank results are:

| Test | Changed nodes | Capped at 6 m | Runtime |
| --- | ---: | ---: | ---: |
| Test 1: ML lower bank | 485,123 | 5,247 | 1:13 |
| Test 2: constant 2 m with intersections | 476,471 | 2,556 | 1:11 |
| Test 3: constant 2 m, forward only | 439,863 | 1,813 | 0:34 |

Among the 476,471 nodes changed by Test 2, Test 1 applies additional depth at
196,988 nodes (41.3%) and more than 0.5 m of additional depth at 93,184 nodes
(19.6%). Test 1 also changes 8,652 nodes that Test 2 leaves unchanged.
Intersection recovery accounts for 36,608 Test 2 changes absent from Test 3.
Test 3 overlaps 438,057 of the 439,308 legacy changed nodes (99.7%), providing
a regression check on the modern forward mapping and lower-bank depth rule.

`hydrofabric_dredge_compare.py` reproducibly builds `comparison_review.gpkg`
and `comparison_summary.json`. The GPKG contains Test 1 versus Test 2 depth
differences, Test 1-only changes, Test 2 intersection additions, capped nodes
and raw greater-than-10 m requests for each test, plus Test 3 versus legacy
forward comparisons. An empty comparison category is recorded with a zero in
the JSON but omitted from the GPKG; for example, the capped runs currently
have no Test 2-only changed nodes. These diagnostic runs do not contain a
modified hgrid unless it is materialized explicitly.

The current `comparison_review.gpkg` layers are:

| Layer | Meaning | Main cause of large differences |
| --- | --- | --- |
| `t1_t2_diff_gt0p5m` | Test 1 and Test 2 final depths differ by more than 0.5 m. Positive `test1_minus_test2_final_m` means Test 1 is deeper. | Test 1 uses `max(bnk_depth, 2 m)`, while Test 2 always uses 2 m. Large positive differences therefore identify hydrofabric bankfull depths well above 2 m; request collisions can select the deepest ML-informed river at junctions. The 6 m cap bounds the applied difference. |
| `t1_only_changed` | Changed by Test 1 but not Test 2. | The existing node is already deep enough for the constant 2 m target, but not for the larger ML target. |
| `t2_intersection_added` | Changed by Test 2 but not forward-only Test 3, isolating intersection recovery. | Junction mesh nodes may be displaced from RiverMapper arc vertices and receive no forward request. Reverse mapping recovers them when multiple rivers pass the radius, half-width, and bank-line screens. |
| `t1_capped_at6m` | Test 1 planned change exceeded 6 m and was applied at the 6 m cap. | Large ML bankfull depths, deep lower-bank references, or the deepest of several colliding requests can produce the upper tail. Major-river nodes are common, but the applied change remains 6 m. |
| `t2_capped_at6m` | Test 2 planned change capped at 6 m. | Even a constant 2 m bank-relative target can require more than 6 m of change when the selected lower bank is much deeper than the existing inner node; boundary mappings and request collisions can amplify this. |
| `t3_capped_at6m` | Test 3 planned change capped at 6 m. | Same lower-bank-versus-inner-node contrast as Test 2, but only for forward requests; no intersection request contributes. |
| `t1_raw_delta_gt10m` | Test 1 uncapped planned change exceeded 10 m; applied change remains 6 m. | Extreme ML depth estimates and major-river lower-bank references dominate; this is a severity view of `t1_capped_at6m`, not an additional rule. |
| `t2_raw_delta_gt10m` | Test 2 uncapped planned change exceeded 10 m. | Usually an anomalously deep lower-bank reference relative to the mapped inner node, sometimes selected through a collision or boundary mapping. |
| `t3_raw_delta_gt10m` | Test 3 uncapped planned change exceeded 10 m. | Same forward lower-bank mapping cause as Test 2, without intersection recovery. |
| `legacy_only_changed` | Changed by the original lower-bank 2 m product but not Test 3. | Legacy used an unrestricted longitude/latitude KD-tree. Test 3 uses projected metric distances, only eligible in-watershed nodes, a 500 m limit, two-bank validation, and the 6 m cap. These differences are concentrated near watershed boundaries or ambiguous nearest-node mappings. |
| `t3_only_changed` | Changed by Test 3 but not the original product. | Modern projected/in-watershed mapping can assign an inner vertex to a different node, and its deterministic deepest-request collision reduction can select nodes not changed by legacy. |
| `legacy_t3_depth_diff` | Changed by both legacy and Test 3, but assigned different final depths. Positive `test3_minus_legacy_m` means Test 3 is deeper. | Small values arise from different distance metrics and nearest-node choices. Large values usually occur near watershed boundaries, where bank vertices map to different eligible nodes, or at collisions where the two methods associate different station requests with the same mesh node. Test 3 capping at 6 m can add another bounded difference from uncapped legacy output. |

The omitted `t2_only_changed` category has zero nodes: after replacing
rejection with capping, every Test 2 changed node is also changed by Test 1.

The reviewed 6 m capped changes have been materialized in both SCHISM GR3 and
SMS 2DM formats:

- Test 1: `01_ml_lowerbank/hgrid_test1_ml_lowerbank.gr3` and
  `hgrid_test1_ml_lowerbank.2dm`.
- Test 2:
  `02_constant2_lowerbank/hgrid_test2_constant2_lowerbank_intersections.gr3`
  and `hgrid_test2_constant2_lowerbank_intersections.2dm`.
- Test 3:
  `03_constant2_lowerbank_forward/hgrid_test3_constant2_lowerbank_forward.gr3`
  and `hgrid_test3_constant2_lowerbank_forward.2dm`.

All six products retain the complete pre-channel topology and 3,080,551-node
coordinate set. Their depths were checked exhaustively against the respective
`dredged_mesh_nodes.parquet` products, including unchanged nodes, and the GR3
and 2DM representations agree.

#### Presentation summary

Using the same v32e pre-connectivity mesh and lower-bank reference, the
ML-informed workflow deepened 485,123 mesh nodes compared with 476,471 for a
constant 2 m baseline. ML estimates supplied additional depth at 41% of the
baseline-changed nodes and activated 8,652 additional nodes, while reverse
intersection recovery added 36,608 changes around multi-river junctions. A
P99-derived 6 m cap limits the upper tail without leaving isolated undredged
nodes, and the modern forward-only result reproduces 99.7% of the legacy
changed-node set, supporting implementation correctness.
