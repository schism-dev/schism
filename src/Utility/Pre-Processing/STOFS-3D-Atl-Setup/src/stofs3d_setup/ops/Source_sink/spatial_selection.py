"""Shared SCHISM element-center selection for source/sink operations."""

import numpy as np
from pyproj import Transformer
from scipy.spatial import cKDTree


def _compute_grid_centers(hgrid) -> tuple[np.ndarray, np.ndarray]:
    """Return element centers without changing the caller's grid attributes."""
    had_xctr = hasattr(hgrid, "xctr")
    had_yctr = hasattr(hgrid, "yctr")

    old_xctr = np.asarray(hgrid.xctr).copy() if had_xctr else None
    old_yctr = np.asarray(hgrid.yctr).copy() if had_yctr else None

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


def _make_transformer() -> Transformer:
    """Build the lon/lat-to-meter transformer used for local screening."""
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
    """Find the nearest candidate center and return distance in meters."""
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
    """Return candidate element IDs with centers within radius_m."""
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
