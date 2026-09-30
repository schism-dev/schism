"""Generic spatial selection and region corrections for source/sink data."""

from pathlib import Path

import numpy as np
from pyproj import Transformer
from scipy.spatial import cKDTree


def _compute_grid_centers(hgrid) -> tuple[np.ndarray, np.ndarray]:
    """
    Return element-center coordinates without leaving new xctr/yctr
    attributes on the caller's hgrid object.
    """
    had_xctr = hasattr(hgrid, "xctr")
    had_yctr = hasattr(hgrid, "yctr")

    old_xctr = (
        np.asarray(hgrid.xctr).copy()
        if had_xctr
        else None
    )
    old_yctr = (
        np.asarray(hgrid.yctr).copy()
        if had_yctr
        else None
    )

    try:
        hgrid.compute_ctr()

        xctr = np.asarray(hgrid.xctr, dtype=float).copy()
        yctr = np.asarray(hgrid.yctr, dtype=float).copy()

    finally:
        if had_xctr:
            hgrid.xctr = old_xctr
        elif hasattr(hgrid, "xctr"):
            delattr(hgrid, "xctr")

        if had_yctr:
            hgrid.yctr = old_yctr
        elif hasattr(hgrid, "yctr"):
            delattr(hgrid, "yctr")

    return xctr, yctr


def _zero_sources_inside_regions(
    regions: list[dict],
    source_eles: list[int],
    source_time_and_data: tuple[np.ndarray, np.ndarray] | None,
    xctr: np.ndarray,
    yctr: np.ndarray,
) -> tuple[tuple[np.ndarray, np.ndarray] | None, int]:
    """
    Set all vsource records to zero for source centers inside SCHISM regions.

    Source element IDs are retained. This function changes neither msource
    nor sink forcing.
    """
    if not regions:
        return source_time_and_data, 0

    if source_time_and_data is None or not source_eles:
        print(
            "[REGION SOURCE ZERO] no source forcing exists; "
            "nothing was changed."
        )
        return source_time_and_data, 0

    from pylib import inside_polygon, read_schism_reg

    source_time, source_data = source_time_and_data
    source_data = np.asarray(source_data, dtype=float).copy()
    source_eles_array = np.asarray(source_eles, dtype=int)
    source_xy = np.c_[
        xctr[source_eles_array - 1],
        yctr[source_eles_array - 1],
    ]
    zeroed_elements: set[int] = set()

    for entry in regions:
        region_file = Path(entry["region_file"])
        region = read_schism_reg(str(region_file))
        inside = np.asarray(
            inside_polygon(source_xy, region.x, region.y)
        ).reshape(-1) == 1
        affected_indices = np.flatnonzero(inside)

        if affected_indices.size == 0:
            print(
                f"[REGION SOURCE ZERO] {entry['name']}: no source "
                f"center found inside {region_file.name}."
            )
            continue

        means_before = np.mean(
            source_data[:, affected_indices],
            axis=0,
        )
        source_data[:, affected_indices] = 0.0
        affected_elements = source_eles_array[affected_indices]
        zeroed_elements.update(int(ele) for ele in affected_elements)

        print(
            f"[REGION SOURCE ZERO] {entry['name']}: set vsource=0 "
            f"for {affected_indices.size} source(s) inside "
            f"{region_file.name}."
        )
        for ele, mean_before in zip(affected_elements, means_before):
            print(
                f"[REGION SOURCE ZERO] element {int(ele)}: "
                f"mean vsource {float(mean_before):.6f} -> 0.0 m3/s"
            )

    return (source_time, source_data), len(zeroed_elements)

def _make_transformer() -> Transformer:
    """
    Build a lon/lat-to-meter transformer.

    EPSG:3857 is used only for local nearest-neighbor screening. For the
    O(1 km) search radii used here, it is adequate and consistent with the
    earlier user-defined patch utilities.
    """
    return Transformer.from_crs(
        "EPSG:4326",
        "EPSG:3857",
        always_xy=True,
    )


def _element_center_tree(
    element_ids: list[int] | np.ndarray,
    xctr: np.ndarray,
    yctr: np.ndarray,
    transformer: Transformer,
) -> tuple[np.ndarray, cKDTree | None]:
    """Build a KDTree for selected SCHISM element centers."""
    element_ids = np.asarray(element_ids, dtype=int).reshape(-1)

    if len(element_ids) == 0:
        return element_ids, None

    xp, yp = transformer.transform(
        xctr[element_ids - 1],
        yctr[element_ids - 1],
    )

    return element_ids, cKDTree(np.c_[xp, yp])


def _nearest_element(
    x: float,
    y: float,
    candidate_element_ids: list[int] | np.ndarray,
    xctr: np.ndarray,
    yctr: np.ndarray,
    transformer: Transformer,
) -> tuple[int | None, float]:
    """Find the nearest candidate element center and return distance in meters."""
    element_ids, tree = _element_center_tree(
        candidate_element_ids,
        xctr,
        yctr,
        transformer,
    )

    if tree is None:
        return None, np.inf

    xq, yq = transformer.transform(x, y)
    distance_m, idx = tree.query([xq, yq])

    return int(element_ids[int(idx)]), float(distance_m)


def _elements_within_radius(
    x: float,
    y: float,
    radius_m: float,
    candidate_element_ids: list[int] | np.ndarray,
    xctr: np.ndarray,
    yctr: np.ndarray,
    transformer: Transformer,
) -> list[int]:
    """Return candidate element IDs whose centers are within radius_m."""
    element_ids, tree = _element_center_tree(
        candidate_element_ids,
        xctr,
        yctr,
        transformer,
    )

    if tree is None:
        return []

    xq, yq = transformer.transform(x, y)
    idxs = tree.query_ball_point([xq, yq], r=float(radius_m))

    return [int(element_ids[int(idx)]) for idx in idxs]
