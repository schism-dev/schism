"""Apply ordered GeoPackage regional tweaks to SCHISM hgrid depths."""

from __future__ import annotations

import copy
from pathlib import Path
import warnings

import geopandas as gpd
import numpy as np
from shapely import intersects_xy

from stofs3d_setup.utils.reg2gpkg import find_polygon_overlaps


MODULE_DIR = Path(__file__).resolve().parent
DEFAULT_REGIONAL_TWEAKS_GPKG = MODULE_DIR / "default_regional_tweaks.gpkg"
DEFAULT_REGIONAL_TWEAKS_V7P2_GPKG = (
    MODULE_DIR / "default_regional_tweaks_v7p2.gpkg"
)

REQUIRED_COLUMNS = {
    "region",
    "apply_order",
    "operation",
    "value_m",
    "min_input_depth_m",
    "max_input_depth_m",
}

DEFAULT_REGIONAL_TWEAKS_v7p2 = {  # incoorporating SECOFS updates, for SMS v27 and after
    'min_5m_ll_noPR': 5,
    'SabinePass': 7,
    'BergenPoint': 5,
     # 'Washington_3': 15,  # deleted
    'Elk_river': 2,
    'Hudson_river': 16,
    'James_river': 2,  # changed
    'NorthEast_river': 5,
    'Rappahannock_river': 6,
    'Susquehanna_river': 10,
    'York_river': 4,  # changed
    'Androscoggin_Kennebec_rivers': 3,
    'Merrimack_river': 3,
    'Patuxent_river': 5,
    'Penobscot_river': 5,
    'Saco_river': 3,
    'StCroix_river': 5,
    'Oyster_landing': 1,
    'st_lawrence1': 10,
    'st_lawrence2': 10,
    'st_lawrence3': 10,
}

DEFAULT_REGIONAL_TWEAKS_shp_v7p4 = {
    'Savannah_upstream': ['set', 0.6],
    'Cooper_upstream': ['deepen', 5.0],
    'Wando_upstream': ['deepen', 2.0],
    'Wachapreague': ['deepen', 2.0],
    'Savannah_marsh_block': ['shallow', 2.0, ('between', -2.0, 0.0)],
    'Savannah_connectivity': ['set', 0.2, ('<', 0.0)],
    'Turkey': ['set', 0.6, ('<', 0.0)],
    'Baltimore': ['set', 2.0, ('<', 2.0)],
    'Fall_River': ['set', 1.5, ('<', 1.5)],
    'Providence': ['set', 1.5, ('<', 1.5)],
}

DEFAULT_REGION_DIR = '/sciclone/schism10/Hgrid_projects/DEMs/regions/'

def tweak_shp_hgrid_depth(hgrid, regional_tweaks=None, regions_dir=DEFAULT_REGION_DIR,):
    """Apply regional depth tweaks using shapefiles."""

    import geopandas as gpd

    if regional_tweaks is None:
        regional_tweaks = DEFAULT_REGIONAL_TWEAKS_shp_v7p4

    hgrid_gdf = gpd.GeoDataFrame(
        geometry=gpd.points_from_xy(hgrid.x, hgrid.y),
        crs='EPSG:4326'
    )

    for region, tweak in regional_tweaks.items():
        mode = tweak[0]
        value = tweak[1]
        condition = tweak[2] if len(tweak) > 2 else None

        # Read region polygon
        gdf = gpd.read_file(f'{regions_dir}/{region}.shp')

        if gdf.crs is None:
            gdf = gdf.set_crs('EPSG:4326')
        else:
            gdf = gdf.to_crs('EPSG:4326')

        # Find hgrid nodes inside the polygon
        idx = gpd.sjoin(
            hgrid_gdf,
            gdf[['geometry']],
            how='inner',
            predicate='within'
        ).index.unique().to_numpy(dtype=int)

        if len(idx) == 0:
            print(f'Warning: no hgrid nodes found in {region}.')
            continue

        # Apply condition
        if condition is not None:
            dp = hgrid.dp[idx]

            if condition[0] == '<':
                idx = idx[dp < condition[1]]

            elif condition[0] == 'between':
                idx = idx[
                    (dp > condition[1]) &
                    (dp < condition[2])
                ]

            else:
                raise ValueError(
                    f'Unknown condition {condition[0]} in {region}.'
                )

        if len(idx) == 0:
            print(f'Warning: no nodes satisfy condition in {region}.')
            continue

        # Apply depth tweak
        if mode == 'set':
            hgrid.dp[idx] = value

        elif mode == 'deepen':
            hgrid.dp[idx] += value

        elif mode == 'shallow':
            hgrid.dp[idx] -= value

        else:
            raise ValueError(
                f'Unknown mode {mode} in {region}.'
            )

        print(
            f'{region}: {mode} by {value} m '
            f'for {len(idx)} nodes.'
        )

    return hgrid

class RegionalTweakWarning(UserWarning):
    """Warn about empty or overlapping regional tweak selections."""


def _numeric_column(regions, column, *, allow_null=False):
    """Return a float array and give a clear error for invalid numeric fields.

    Parameters
    ----------
    regions : geopandas.GeoDataFrame
        Regional tweak features read from the GeoPackage.
    column : str
        Name of the numeric field to convert.
    allow_null : bool, optional
        Whether ``NULL``/``NaN`` means that the field is unbounded.
    """

    try:
        values = regions[column].astype(float).to_numpy()
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{column} must contain numeric values") from exc

    null = np.isnan(values)
    if not allow_null and null.any():
        raise ValueError(f"{column} cannot contain NULL values")
    if not np.isfinite(values[~null]).all():
        raise ValueError(f"{column} must contain finite values")
    return values


def _load_regions(gpkg_file, layer, hgrid_crs):
    """Read, validate, order, and reproject regional tweak polygons.

    The modern schema is required in full so configuration errors are found
    before the hgrid is copied or changed. Extra provenance columns, such as
    ``source_file``, are allowed and ignored by the tweak logic.
    """

    read_options = {"layer": layer} if layer is not None else {}
    regions = gpd.read_file(Path(gpkg_file), **read_options)
    if regions.empty:
        raise ValueError("The regional-tweak GeoPackage contains no features")
    if regions.crs is None:
        raise ValueError("The regional-tweak GeoPackage has no CRS")

    missing = sorted(REQUIRED_COLUMNS - set(regions.columns))
    if missing:
        raise ValueError(
            "The regional-tweak GeoPackage is missing required columns: "
            + ", ".join(missing)
        )

    if regions["region"].isna().any():
        raise ValueError("region cannot contain NULL values")
    regions["region"] = regions["region"].astype(str).str.strip()
    if (regions["region"] == "").any():
        raise ValueError("region cannot contain empty names")
    if regions["region"].duplicated().any():
        duplicates = sorted(
            regions.loc[regions["region"].duplicated(False), "region"].unique()
        )
        raise ValueError("region names must be unique: " + ", ".join(duplicates))

    apply_order = _numeric_column(regions, "apply_order")
    if not np.equal(apply_order, np.floor(apply_order)).all():
        raise ValueError("apply_order must contain integers")
    if len(np.unique(apply_order)) != len(apply_order):
        raise ValueError("apply_order values must be unique")
    regions["apply_order"] = apply_order.astype(np.int64)

    if regions["operation"].isna().any():
        raise ValueError("operation cannot contain NULL values")
    regions["operation"] = regions["operation"].astype(str).str.strip().str.lower()
    invalid_operations = sorted(set(regions["operation"]) - SUPPORTED_OPERATIONS)
    if invalid_operations:
        raise ValueError(
            "Unsupported operations: "
            + ", ".join(invalid_operations)
            + f"; expected {sorted(SUPPORTED_OPERATIONS)}"
        )

    regions["value_m"] = _numeric_column(regions, "value_m")
    for column in ("min_input_depth_m", "max_input_depth_m"):
        regions[column] = _numeric_column(regions, column, allow_null=True)
    invalid_bounds = (
        regions["min_input_depth_m"].notna()
        & regions["max_input_depth_m"].notna()
        & (regions["min_input_depth_m"] > regions["max_input_depth_m"])
    )
    if invalid_bounds.any():
        names = ", ".join(regions.loc[invalid_bounds, "region"])
        raise ValueError(f"Minimum input depth exceeds maximum for: {names}")

    for index, geometry in regions.geometry.items():
        name = regions.at[index, "region"]
        if geometry is None or geometry.is_empty:
            raise ValueError(f"Region {name!r} has empty geometry")
        if geometry.geom_type not in POLYGON_TYPES:
            raise ValueError(
                f"Region {name!r} has {geometry.geom_type} geometry; "
                "expected Polygon or MultiPolygon"
            )
        if not geometry.is_valid:
            raise ValueError(f"Region {name!r} has invalid geometry")

    regions = regions.sort_values("apply_order").reset_index(drop=True)
    return regions.to_crs(hgrid_crs)


def _select_nodes(regions, x, y, input_depth):
    """Select spatial and depth-filtered hgrid nodes for every region.

    Polygon boundaries are included. Depth bounds are inclusive and always
    use the original input depths, so earlier operations do not change which
    nodes a later feature selects. The separate spatial selections preserve
    the original diagnostic-mask behavior.
    """

    selections = []
    spatial_selections = []
    for row in regions.itertuples():
        spatial_indices = np.flatnonzero(intersects_xy(row.geometry, x, y))
        if spatial_indices.size == 0:
            warnings.warn(
                f"Region {row.region!r} contains no hgrid nodes",
                RegionalTweakWarning,
                stacklevel=3,
            )
        indices = spatial_indices
        if not np.isnan(row.min_input_depth_m):
            indices = indices[input_depth[indices] >= row.min_input_depth_m]
        if not np.isnan(row.max_input_depth_m):
            indices = indices[input_depth[indices] <= row.max_input_depth_m]
        spatial_selections.append(spatial_indices)
        selections.append(indices)
    return selections, spatial_selections


def _warn_about_overlaps(regions, selections, node_count):
    """Warn when polygon or selected-node overlaps make order significant.

    Positive-area polygon overlaps and shared selected nodes are both checked.
    The warning prints the actual ascending execution order for easy review.
    """

    polygon_overlaps = find_polygon_overlaps(regions)
    selection_count = np.zeros(node_count, dtype=np.uint16)
    for indices in selections:
        selection_count[indices] += 1
    repeated_nodes = selection_count > 1
    if not polygon_overlaps and not repeated_nodes.any():
        return

    order_by_region = dict(zip(regions["region"], regions["apply_order"]))
    lines = [
        "Overlapping regional tweaks are order-dependent. Verify apply_order "
        "before using this output:"
    ]
    for overlap in polygon_overlaps:
        pair = sorted(
            (overlap.first_region, overlap.second_region),
            key=order_by_region.get,
        )
        lines.append(
            f"  {pair[0]} [{order_by_region[pair[0]]}] -> "
            f"{pair[1]} [{order_by_region[pair[1]]}]"
        )
    if repeated_nodes.any():
        affected = [
            f"{row.region} [{row.apply_order}]"
            for row, indices in zip(regions.itertuples(), selections)
            if np.any(repeated_nodes[indices])
        ]
        lines.append(
            f"  {int(repeated_nodes.sum())} hgrid nodes are selected more than "
            "once by: " + ", ".join(affected)
        )
    warnings.warn("\n".join(lines), RegionalTweakWarning, stacklevel=3)


def shape_tweak(
    hgrid,
    gpkg_file,
    *,
    layer=None,
    hgrid_crs="EPSG:4326",
):
    """Apply ordered GeoPackage depth operations and return two hgrid copies.

    Required GeoPackage fields
    --------------------------
    region : str
        Unique diagnostic name for the polygon.
    apply_order : int
        Ascending execution order. Larger values apply later.
    operation : {"set", "add"}
        ``set`` assigns ``value_m``; ``add`` adds it to the current depth.
    value_m : float
        Assignment or increment in meters. For ``set`` only, values greater
        than 999999 set every selected node to the selected region's current
        maximum depth; values less than -999999 use its current minimum.
    min_input_depth_m, max_input_depth_m : float or NULL
        Optional inclusive bounds evaluated against the original hgrid depth.

    Parameters
    ----------
    hgrid : object
        SCHISM grid-like object with one-dimensional ``x``, ``y``, and ``dp``
        arrays. The input object is not modified.
    gpkg_file : path-like
        GeoPackage containing the regional tweak polygons and fields above.
    layer : str, optional
        GeoPackage layer name. The first layer is used when omitted.
    hgrid_crs : str, optional
        CRS of the hgrid coordinates. Defaults to longitude/latitude WGS 84.

    Returns
    -------
    tweaked_hgrid, touched_hgrid : tuple
        A depth-modified copy and a copy whose depth is 1 at every selected
        polygon-covered node and 0 elsewhere. The diagnostic mask includes
        nodes excluded by input-depth bounds, matching the original function.
    """

    x = np.asarray(hgrid.x, dtype=float)
    y = np.asarray(hgrid.y, dtype=float)
    input_depth = np.asarray(hgrid.dp).copy()
    if x.ndim != 1 or y.ndim != 1 or input_depth.ndim != 1:
        raise ValueError("hgrid x, y, and dp must be one-dimensional arrays")
    if not (x.size == y.size == input_depth.size):
        raise ValueError("hgrid x, y, and dp must have equal lengths")
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("hgrid coordinates must be finite")
    if not np.isfinite(input_depth).all():
        raise ValueError("hgrid depths must be finite")

    regions = _load_regions(gpkg_file, layer, hgrid_crs)
    selections, spatial_selections = _select_nodes(regions, x, y, input_depth)
    _warn_about_overlaps(regions, selections, input_depth.size)

    tweaked = copy.deepcopy(hgrid)
    touched = copy.deepcopy(hgrid)
    touched.dp[:] = 0
    for row, indices, spatial_indices in zip(
        regions.itertuples(), selections, spatial_selections
    ):
        touched.dp[spatial_indices] = 1
        if indices.size == 0:
            print(f"Skipped {row.region} [{row.apply_order}]: no selected nodes")
            continue
        if row.operation == "set" and row.value_m > 999999:
            regional_max = float(np.max(tweaked.dp[indices]))
            tweaked.dp[indices] = regional_max
            description = f"set depth to regional maximum {regional_max:g} m"
        elif row.operation == "set" and row.value_m < -999999:
            regional_min = float(np.min(tweaked.dp[indices]))
            tweaked.dp[indices] = regional_min
            description = f"set depth to regional minimum {regional_min:g} m"
        elif row.operation == "set":
            tweaked.dp[indices] = row.value_m
            description = f"set depth to {row.value_m:g} m"
        else:
            tweaked.dp[indices] += row.value_m
            description = f"added {row.value_m:g} m"
        print(
            f"Applied {row.region} [{row.apply_order}]: {description} "
            f"at {indices.size} nodes"
        )

    print(f"Total points tweaked: {int(np.count_nonzero(touched.dp))}")
    return tweaked, touched
