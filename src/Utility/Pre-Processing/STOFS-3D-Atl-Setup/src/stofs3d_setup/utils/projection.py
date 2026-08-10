"""Reliable pyproj-backed projection with a persistent local grid cache."""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import tempfile
import time
from urllib.parse import urlparse
from urllib.request import urlopen

import numpy as np
from pyproj import CRS, Transformer, datadir, network
from pyproj.aoi import AreaOfInterest
from pyproj.transformer import TransformerGroup


PROJ_DATA_DIR = Path(__file__).resolve().parent / "proj_data"


class ProjectionError(RuntimeError):
    """Raised when a coordinate transformation fails or yields invalid values."""


def validate_finite_xy(x, y, context="coordinate transformation"):
    """Raise ``ProjectionError`` when transformed x/y coordinates are invalid."""
    try:
        x_array, y_array = np.broadcast_arrays(
            np.asarray(x, dtype=float), np.asarray(y, dtype=float)
        )
    except ValueError as exc:
        raise ProjectionError(
            f"{context} returned incompatible x/y shapes: "
            f"{np.shape(x)} and {np.shape(y)}."
        ) from exc

    finite = np.isfinite(x_array) & np.isfinite(y_array)
    if np.all(finite):
        return

    bad_coordinates = np.flatnonzero(~finite)
    preview = "; ".join(
        f"index {np.unravel_index(index, finite.shape)}: "
        f"x={x_array.flat[index]!r}, y={y_array.flat[index]!r}"
        for index in bad_coordinates[:5]
    )
    raise ProjectionError(
        f"{context} produced non-finite coordinates at "
        f"{bad_coordinates.size} of {finite.size} positions. "
        f"First failures: {preview}."
    )


def validate_finite_coordinates(gd, context="coordinate transformation"):
    """Raise a descriptive error if a SCHISM grid has invalid x/y coordinates."""
    validate_finite_xy(gd.x, gd.y, context)


def _geometry_coordinate_pairs(geometry):
    """Yield x/y pairs from a Shapely geometry without requiring Shapely v2."""
    if geometry is None or geometry.is_empty:
        return

    mapping = geometry.__geo_interface__
    if mapping.get("type") == "GeometryCollection":
        for child in geometry.geoms:
            yield from _geometry_coordinate_pairs(child)
        return
    coordinates = mapping.get("coordinates", ())

    def walk(value):
        if isinstance(value, np.ndarray):
            value = value.tolist()
        if not value:
            return
        first = value[0]
        if np.isscalar(first):
            if len(value) >= 2:
                yield value[0], value[1]
            return
        for child in value:
            yield from walk(child)

    yield from walk(coordinates)


def validate_geometries_finite(geometries, context="geometry transformation"):
    """Validate transformed GeoPandas/Shapely geometries.

    Empty and ``None`` geometries are retained because they can be legitimate
    input. Any non-finite x/y coordinate in a non-empty geometry raises.
    """
    items = geometries.items() if hasattr(geometries, "items") else enumerate(geometries)
    failures = []
    checked = 0
    for label, geometry in items:
        for x, y in _geometry_coordinate_pairs(geometry):
            checked += 1
            if not np.isfinite(x) or not np.isfinite(y):
                failures.append((label, x, y))
                if len(failures) == 5:
                    break
        if len(failures) == 5:
            break

    if failures:
        preview = "; ".join(
            f"geometry {label!r}: x={x!r}, y={y!r}"
            for label, x, y in failures
        )
        raise ProjectionError(
            f"{context} produced non-finite geometry coordinates after "
            f"checking {checked} coordinate pairs. First failures: {preview}."
        )


def project_geodataframe(data, target_crs, context="GeoPandas transformation"):
    """Project GeoPandas data and reject any non-finite output coordinates.

    This intentionally provides validation only. Unlike ``project_grid()``, it
    does not download or retry datum grids because GeoPandas owns the transform.
    """
    try:
        projected = data.to_crs(target_crs)
    except Exception as exc:
        raise ProjectionError(
            f"{context} failed while transforming to {target_crs!r}: {exc}"
        ) from exc
    geometries = projected.geometry if hasattr(projected, "geometry") else projected
    validate_geometries_finite(geometries, context)
    return projected


def transform_points(transformer, x, y, context="point transformation"):
    """Transform raw point arrays and convert failures to ``ProjectionError``.

    This is validation only; it intentionally does not download or retry PROJ
    grids. ``transformer`` can be any object exposing ``transform(x, y)``.
    """
    try:
        transformed_x, transformed_y = transformer.transform(x, y)
    except Exception as exc:
        raise ProjectionError(f"{context} failed: {exc}") from exc
    validate_finite_xy(transformed_x, transformed_y, context)
    return transformed_x, transformed_y


def _register_proj_data(proj_data_dir):
    """Add downloaded grids to the current pyproj search path once."""
    proj_data_dir = str(Path(proj_data_dir).resolve())
    search_paths = datadir.get_data_dir().split(os.pathsep)
    if proj_data_dir not in search_paths:
        datadir.append_data_dir(proj_data_dir)


def _source_geographic_bbox(x, y, source_crs):
    """Compute an AOI without performing the failing cross-datum operation."""
    source_crs = CRS.from_user_input(source_crs)
    transformer = Transformer.from_crs(
        source_crs, source_crs.geodetic_crs, always_xy=True
    )
    lon, lat = transformer.transform(x, y)
    finite = np.isfinite(lon) & np.isfinite(lat)
    if not np.all(finite):
        raise ProjectionError(
            "Cannot determine the projection area of interest because the "
            "source CRS-to-geodetic conversion produced non-finite values."
        )
    return AreaOfInterest(
        west_lon_degree=float(np.min(lon)),
        south_lat_degree=float(np.min(lat)),
        east_lon_degree=float(np.max(lon)),
        north_lat_degree=float(np.max(lat)),
    )


def _find_missing_grids(source_crs, target_crs, area_of_interest):
    """Return directly downloadable grids for unavailable candidate operations."""
    network_was_enabled = network.is_network_enabled()
    network.set_network_enabled(False)
    try:
        group = TransformerGroup(
            CRS.from_user_input(source_crs),
            CRS.from_user_input(target_crs),
            always_xy=True,
            area_of_interest=area_of_interest,
        )
    finally:
        network.set_network_enabled(network_was_enabled)

    grids = {}
    for operation in group.unavailable_operations:
        for grid in operation.grids:
            if grid.available:
                continue
            if not grid.direct_download or not grid.url:
                continue
            grids[grid.short_name] = grid.url
    return sorted(grids.items())


def _download_grid(name, url, proj_data_dir, attempts=3):
    """Download one official PROJ grid atomically."""
    parsed = urlparse(url)
    if parsed.scheme != "https" or parsed.netloc != "cdn.proj.org":
        raise ProjectionError(f"Refusing unexpected PROJ grid URL: {url}")
    if Path(parsed.path).name != name or Path(name).name != name:
        raise ProjectionError(f"Unsafe PROJ grid filename: {name}")

    destination = proj_data_dir / name
    if destination.is_file():
        return False

    last_error = None
    for attempt in range(1, attempts + 1):
        temporary_name = None
        try:
            with tempfile.NamedTemporaryFile(
                prefix=f".{name}.", suffix=".part",
                dir=proj_data_dir, delete=False,
            ) as temporary:
                temporary_name = Path(temporary.name)
                with urlopen(url, timeout=120) as response:
                    shutil.copyfileobj(response, temporary)
            temporary_name.replace(destination)
            return True
        except Exception as exc:
            last_error = exc
            if temporary_name is not None:
                temporary_name.unlink(missing_ok=True)
            if attempt < attempts:
                time.sleep(attempt)

    raise ProjectionError(
        f"Failed to download required PROJ grid {name} from {url} "
        f"after {attempts} attempts: {last_error}"
    ) from last_error


def _download_missing_grids(grids, proj_data_dir):
    proj_data_dir.mkdir(parents=True, exist_ok=True)
    print(
        f"Projection needs {len(grids)} PROJ transformation grid(s). "
        f"Downloading missing files to {proj_data_dir}"
    )
    print(
        "NOTE: files under this proj_data directory are a local runtime cache "
        "and are intentionally not tracked by Git."
    )
    downloaded = 0
    for index, (name, url) in enumerate(grids, start=1):
        if _download_grid(name, url, proj_data_dir):
            downloaded += 1
            print(f"Downloaded PROJ grid {index}/{len(grids)}: {name}")
    print(f"PROJ grid cache ready; {downloaded} new file(s) downloaded.")


def project_grid(gd, source_crs, target_crs, wdir=None):
    """Project a SCHISM grid, downloading missing PROJ grids and retrying once.

    The first projection uses the existing pipeline unchanged. If it raises or
    produces non-finite coordinates, the original coordinates are restored,
    required grids are downloaded into ``utils/proj_data``, and the projection
    is retried. ``wdir`` is accepted for caller compatibility and diagnostics;
    the shared package cache is intentionally used for all working directories.
    """
    del wdir  # downloads are shared under this module's utils directory
    source_x = np.asarray(gd.x, dtype=float).copy()
    source_y = np.asarray(gd.y, dtype=float).copy()
    context = f"projection from {source_crs} to {target_crs}"

    if PROJ_DATA_DIR.is_dir():
        _register_proj_data(PROJ_DATA_DIR)

    try:
        gd.proj(prj0=source_crs, prj1=target_crs)
        validate_finite_coordinates(gd, context)
        return gd
    except Exception as initial_error:
        gd.x, gd.y = source_x.copy(), source_y.copy()
        print(f"{context} failed: {initial_error}")

    try:
        area_of_interest = _source_geographic_bbox(
            source_x, source_y, source_crs
        )
        grids = _find_missing_grids(
            source_crs, target_crs, area_of_interest
        )
        if not grids:
            raise ProjectionError(
                "PROJ reported no directly downloadable missing grids for "
                f"{context}."
            )
        _download_missing_grids(grids, PROJ_DATA_DIR)
        _register_proj_data(PROJ_DATA_DIR)

        gd.x, gd.y = source_x.copy(), source_y.copy()
        gd.proj(prj0=source_crs, prj1=target_crs)
        validate_finite_coordinates(gd, f"{context} retry")
        print(f"{context} succeeded after installing PROJ grids.")
        return gd
    except Exception as retry_error:
        gd.x, gd.y = source_x, source_y
        if isinstance(retry_error, ProjectionError):
            raise
        raise ProjectionError(
            f"{context} failed after installing available PROJ grids: "
            f"{retry_error}"
        ) from retry_error
