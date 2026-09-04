# Atlas skeleton-to-mesh mapping

Open `atlas_skeleton_vertex_mesh_matches.gpkg` in QGIS. The main layers are:

- `river_centerlines`: coarsened v19.1 river-skeleton LineStrings.
- `skeleton_vertex_mapping`: one point for every centerline coordinate and its
  selected Atlas mesh node.
- `mapping_links`: lines from skeleton vertices to their selected mesh nodes.
- `in_river_mesh_nodes`: mesh nodes identified from RiverMapper inner arcs.
- `outside_hgrid_vertices`: skeleton vertices outside the hgrid footprint.
- `long_links_gt_500m`: mappings longer than 500 m that deserve review.
- `skeleton_rivermapper_correspondence`: the RiverMapper river assigned to
  each skeleton LineString.

## Associate a point with its centerline

Use these fields:

```text
river_centerlines.river_idx
    = skeleton_vertex_mapping.skeleton_river_idx
```

`vertex_idx` is the zero-based coordinate position along the LineString. The
unique identifier for a skeleton vertex is therefore:

```text
(skeleton_river_idx, vertex_idx)
```

In QGIS, create a one-to-many project relation from
`river_centerlines.river_idx` to
`skeleton_vertex_mapping.skeleton_river_idx`. Selecting a centerline then
provides access to all its mapped vertex points.

## Interpret a mapped vertex

The most useful fields in `skeleton_vertex_mapping` are:

- `mapping_method = strict_inner_river_node`: mapped through the corresponding
  RiverMapper river and one of its inner arcs.
- `mapping_method = nearest_mesh_node`: no viable corresponding inner node was
  found, so the nearest unrestricted mesh node was used.
- `mapping_method = outside_hgrid`: outside the mesh footprint; no mesh node
  was assigned.
- `mesh_node_id`: one-based SCHISM node ID; `-1` means unmapped.
- `mesh_position`: zero-based array position used by the Python workflow.
- `distance_m`: distance from the skeleton vertex to the selected mesh node.
- `source_arc_river_idx`, `source_local_arc_idx`, and `source_station_idx`:
  corresponding RiverMapper inner-arc location.

Bank arcs are not used to seed the strict-inner candidate set. Shared nodes at
river intersections remain valid so tributaries can connect through receiving
rivers. No watershed or effective-watershed mask is applied to this mapping.

For a quick review, symbolize `skeleton_vertex_mapping` by `mapping_method`,
overlay `mapping_links`, and inspect `long_links_gt_500m` and
`outside_hgrid_vertices` separately.
