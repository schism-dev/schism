"""Zero existing SCHISM source flow inside configured regions."""

from pathlib import Path

import numpy as np

from stofs3d_setup.ops.Source_sink.correction_config import (
    load_source_sink_corrections,
)
from stofs3d_setup.ops.Source_sink.source_sink_components import (
    _build_source_sink,
    _copy_sink_components,
    _copy_source_components,
)
from stofs3d_setup.ops.Source_sink.spatial_selection import _compute_grid_centers
from stofs3d_setup.ops.Source_sink.Spatial_corrections.region_config import (
    zero_source_regions as configured_zero_source_regions,
)


def _zero_sources_inside_regions(
    regions: list[dict],
    source_eles: list[int],
    source_time_and_data: tuple[np.ndarray, np.ndarray] | None,
    xctr: np.ndarray,
    yctr: np.ndarray,
) -> tuple[tuple[np.ndarray, np.ndarray] | None, int]:
    """Zero source flows inside regions without changing tracers or sinks."""
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


def zero_sources_in_regions(base_ss, hgrid, regions: list[dict]):
    """Return a copy with source flow zeroed inside configured regions."""
    source_eles, source_values, msource_values = _copy_source_components(
        base_ss
    )
    sink_eles, sink_values = _copy_sink_components(base_ss)
    xctr, yctr = _compute_grid_centers(hgrid)

    source_values, zeroed_count = _zero_sources_inside_regions(
        regions=regions,
        source_eles=source_eles,
        source_time_and_data=source_values,
        xctr=xctr,
        yctr=yctr,
    )
    corrected_ss = _build_source_sink(
        source_eles=source_eles,
        source_time_and_data=source_values,
        msource_data_list=msource_values,
        sink_eles=sink_eles,
        sink_time_and_data=sink_values,
    )
    return corrected_ss, zeroed_count


def zero_configured_source_regions(base_ss, hgrid, correction_info):
    """Load one region-zeroing configuration and apply it to source flow."""
    corrections = load_source_sink_corrections(correction_info)
    config_dir = (
        Path.cwd()
        if isinstance(correction_info, dict)
        else Path(correction_info).expanduser().resolve().parent
    )
    regions = configured_zero_source_regions(
        corrections,
        config_dir=config_dir,
    )
    return zero_sources_in_regions(
        base_ss=base_ss,
        hgrid=hgrid,
        regions=regions,
    )
