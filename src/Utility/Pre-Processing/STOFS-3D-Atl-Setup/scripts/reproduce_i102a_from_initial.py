#!/usr/bin/env python3
"""Reproduce I102a source/sink stages after cached NWM and USGS generation.

Run manually on a compute node. The completed initial results and HJ's
reference are read-only inputs; all generated files go beneath --output.
The run stops after artificial-island corrections, before constant sinks.
"""

from __future__ import annotations

import argparse
import filecmp
import json
import shutil
import time
from datetime import datetime
from pathlib import Path

import netCDF4
import numpy as np
from pylib import schism_grid
from pylib_experimental.schism_file import source_sink

from stofs3d_setup.config.stofs3d_atl_config import ConfigStofs3dAtlantic
from stofs3d_setup.ops.Source_sink.Patch_artificial_island.patch_artificial_island_source_sink import (
    apply_artificial_island_corrections,
    zero_artificial_island_sources_after_replace_USGS_before_relocation,
)
from stofs3d_setup.ops.Source_sink.Relocate.relocate_source_feeder import (
    relocate_sources2,
)
from stofs3d_setup.ops.Source_sink.Replace_with_USGS.source_overrides import (
    apply_source_flow_overrides,
    apply_source_temperature_overrides,
    load_source_override_points,
)
from stofs3d_setup.ops.Source_sink.Replace_with_USGS.source_temperature import (
    replace_source_temperatures_with_usgs,
)
from stofs3d_setup.ops.Source_sink.assemble_source_sink import (
    gen_relocated_source,
)
from stofs3d_setup.ops.Source_sink.Spatial_corrections.region_zeroing import (
    zero_configured_source_regions,
)
from stofs3d_setup.utils.utils import STOFS3D_ATL_STATES


NWM_SHAPEFILE = Path("/sciclone/schism10/Hgrid_projects/NWM/ecgc/ecgc.shp")
TEXT_OUTPUTS = ("source_sink.in", "vsource.th", "msource.th", "vsink.th")
STATIC_INITIAL_FILES = (
    "source_sink.in",
    "sources.json",
    "sinks.json",
    "msource.th",
    "vsink.th",
)
VSOURCE_ATOL_M3_PER_S = 1e-10


def _require_files(paths: list[Path]) -> None:
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Required input file(s) missing: {missing}")


def _compare_files(
    actual_dir: Path,
    reference_dir: Path,
    names: tuple[str, ...],
) -> dict[str, bool]:
    return {
        name: (
            (actual_dir / name).is_file()
            and (reference_dir / name).is_file()
            and filecmp.cmp(
                actual_dir / name,
                reference_dir / name,
                shallow=False,
            )
        )
        for name in names
    }


def _compare_source_nc(actual: Path, reference: Path) -> dict:
    """Compare every dimension and variable, tolerating known flow roundoff."""
    result = {"dimensions_match": False, "variables": {}}
    with netCDF4.Dataset(actual) as current, netCDF4.Dataset(reference) as expected:
        current_dims = {name: len(dim) for name, dim in current.dimensions.items()}
        expected_dims = {name: len(dim) for name, dim in expected.dimensions.items()}
        result["dimensions_match"] = current_dims == expected_dims
        result["actual_dimensions"] = current_dims
        result["reference_dimensions"] = expected_dims
        for name in sorted(set(current.variables) | set(expected.variables)):
            if name not in current.variables or name not in expected.variables:
                result["variables"][name] = {"match": False, "reason": "missing"}
                continue
            actual_values = np.asarray(current.variables[name][:])
            reference_values = np.asarray(expected.variables[name][:])
            same_shape = actual_values.shape == reference_values.shape
            if name == "vsource" and same_shape:
                max_abs = float(np.max(np.abs(actual_values - reference_values)))
                match = bool(
                    np.allclose(
                        actual_values,
                        reference_values,
                        rtol=0,
                        atol=VSOURCE_ATOL_M3_PER_S,
                        equal_nan=True,
                    )
                )
                result["variables"][name] = {
                    "match": match,
                    "max_abs_difference_m3_per_s": max_abs,
                    "absolute_tolerance_m3_per_s": VSOURCE_ATOL_M3_PER_S,
                }
            else:
                result["variables"][name] = {
                    "match": bool(
                        same_shape and np.array_equal(actual_values, reference_values)
                    )
                }
    result["match"] = result["dimensions_match"] and all(
        item["match"] for item in result["variables"].values()
    )
    return result


def _timed(label: str, operation, timings: dict[str, float]):
    started = time.monotonic()
    result = operation()
    timings[label] = time.monotonic() - started
    print(f"[TIMING] {label}: {timings[label]:.1f} s", flush=True)
    return result


def _stage_initial_results(
    initial_results: Path,
    original_dir: Path,
    hgrid: Path,
) -> None:
    """Keep static files shared, but isolate the flow that zeroing rewrites."""
    original_dir.mkdir()
    (original_dir / "hgrid.gr3").symlink_to(hgrid)
    initial_original = initial_results / "original_source_sink"
    for name in STATIC_INITIAL_FILES:
        (original_dir / name).symlink_to((initial_original / name).resolve())
    shutil.copyfile(
        initial_results / "USGS_adjusted_sources" / "adjusted_vsource.th",
        original_dir / "vsource.th",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--initial-results", type=Path, required=True)
    parser.add_argument("--reference-input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--preflight", action="store_true")
    args = parser.parse_args()

    initial_results = args.initial_results.expanduser().resolve(strict=True)
    reference_input = args.reference_input.expanduser().resolve(strict=True)
    output = args.output.expanduser().absolute()
    if output.exists() or output.is_symlink():
        parser.error(f"Output already exists; refusing to overwrite: {output}")

    initial_manifest_file = initial_results / "run_result.json"
    _require_files([initial_manifest_file])
    initial_manifest = json.loads(initial_manifest_file.read_text(encoding="utf-8"))
    if initial_manifest.get("completed_through") != "NWM-to-USGS replacement":
        raise ValueError("Initial results did not complete NWM-to-USGS replacement")
    run_config = initial_manifest["config"]
    if run_config.get("profile") != "v7.4_2017":
        raise ValueError(f"Expected v7.4_2017 input, got {run_config.get('profile')}")

    hgrid_file = reference_input / "hgrid.gr3"
    usgs_cache = Path(initial_manifest["resolved_inputs"]["usgs_cache"])
    nwm_cache = Path(initial_manifest["resolved_inputs"]["nwm_cache"])
    reference_ss = reference_input / "Source_sink"
    initial_original = initial_results / "original_source_sink"
    initial_pristine = initial_results / "original_source_sink_before_USGS_adjustment"
    adjusted_flow = initial_results / "USGS_adjusted_sources" / "adjusted_vsource.th"
    required = [
        hgrid_file,
        NWM_SHAPEFILE,
        adjusted_flow,
        initial_pristine / "vsource.th",
        *(initial_original / name for name in STATIC_INITIAL_FILES),
        *(
            reference_ss / "original_source_sink" / name
            for name in STATIC_INITIAL_FILES
        ),
        reference_ss / "original_source_sink_before_USGS_adjustment" / "vsource.th",
        reference_ss
        / "zeroed_sources_after_replace_USGS_before_relocation"
        / "vsource_zeroed.th",
        reference_ss / "relocated_source_sink" / "sources.json",
        *(
            reference_ss / "patch_artificial_island_source_sink" / name
            for name in (*TEXT_OUTPUTS, "source.nc")
        ),
    ]
    _require_files(required)
    if not usgs_cache.is_dir() or not nwm_cache.is_dir():
        raise FileNotFoundError("NWM or USGS cache directory is unavailable")
    initial_hgrid = Path(initial_manifest["resolved_inputs"]["hgrid"])
    if hgrid_file.resolve() != initial_hgrid.resolve():
        raise ValueError("Initial run and reference use different hgrid files")

    start_time = datetime.fromisoformat(run_config["startdate"])
    config = ConfigStofs3dAtlantic.v7p4().model_copy(
        update={
            "startdate": start_time,
            "rnday": int(run_config["rnday"]),
            "nwm_cache_folder": nwm_cache,
            "usgs_cache_folder": usgs_cache,
        }
    )
    if not config.relocate_source or config.reuse_source_json:
        raise ValueError(
            "This reproduction requires v7.4 relocation without source JSON reuse"
        )
    _require_files(
        [
            Path(config.feeder_info_file),
            Path(config.selected_source_override_info),
            Path(config.zero_source_region_info),
            Path(config.artificial_island_source_sink_info),
        ]
    )
    if args.preflight:
        print("Preflight passed; no output was created.")
        return 0

    initial_comparison = _compare_files(
        initial_original,
        reference_ss / "original_source_sink",
        STATIC_INITIAL_FILES,
    )
    initial_comparison["pristine_vsource.th"] = filecmp.cmp(
        initial_pristine / "vsource.th",
        reference_ss / "original_source_sink_before_USGS_adjustment" / "vsource.th",
        shallow=False,
    )

    output.mkdir(parents=True, exist_ok=False)
    (output / "hgrid.gr3").symlink_to(hgrid_file)
    source_sink_dir = output / "Source_sink"
    source_sink_dir.mkdir()
    original_dir = source_sink_dir / "original_source_sink"
    _stage_initial_results(initial_results, original_dir, hgrid_file)

    timings: dict[str, float] = {}
    counts: dict[str, int] = {}
    comparisons: dict[str, dict[str, bool]] = {}
    started = time.monotonic()
    prezero_dir = (
        source_sink_dir / "zeroed_sources_after_replace_USGS_before_relocation"
    )
    _timed(
        "pre_relocation_zeroing",
        lambda: zero_artificial_island_sources_after_replace_USGS_before_relocation(
            source_sink_dir=original_dir,
            hgrid_file=hgrid_file,
            patch_info_file=config.artificial_island_source_sink_info,
            output_dir=prezero_dir,
        ),
        timings,
    )
    comparisons["pre_relocation_zeroing"] = _compare_files(
        prezero_dir,
        reference_ss / prezero_dir.name,
        ("vsource_zeroed.th",),
    )

    relocated_dir = source_sink_dir / "relocated_source_sink"
    relocated_dir.mkdir()
    (relocated_dir / "hgrid.gr3").symlink_to(hgrid_file)
    _timed(
        "relocation_mapping",
        lambda: relocate_sources2(
            old_ss_dir=str(original_dir),
            outdir=str(relocated_dir),
            main_hgrid_has_feeder=False,
            feeder_info_file=str(config.feeder_info_file),
            hgrid_fname=str(hgrid_file),
            allow_neglection=False,
            max_search_radius=2100,
            mandatory_sources_coor=np.asarray(config.mandatory_sources_coor).copy(),
        ),
        timings,
    )
    comparisons["relocation_mapping"] = _compare_files(
        relocated_dir,
        reference_ss / relocated_dir.name,
        ("sources.json",),
    )

    def generate_relocated_forcing():
        vsource, msource = gen_relocated_source(
            original_source_sink_dir=str(original_dir),
            relocated_source_sink_dir=str(relocated_dir),
        )
        relocated_ss = source_sink(vsource=vsource, vsink=None, msource=msource)
        relocated_ss.writer(str(relocated_dir))
        return relocated_ss

    base_ss = _timed("relocated_forcing", generate_relocated_forcing, timings)
    comparisons["relocated_forcing"] = _compare_files(
        relocated_dir,
        reference_ss / relocated_dir.name,
        ("source_sink.in", "vsource.th", "msource.th"),
    )

    hgrid = schism_grid(str(hgrid_file))
    override_points = load_source_override_points(config.selected_source_override_info)
    base_ss, counts["selected_source_flow_overrides"] = _timed(
        "selected_source_flow_overrides",
        lambda: apply_source_flow_overrides(
            base_ss=base_ss,
            hgrid=hgrid,
            points=override_points,
            start_time=config.startdate,
            usgs_cache_folder=usgs_cache,
        ),
        timings,
    )
    base_ss, counts["automatic_temperature_replacements"] = _timed(
        "automatic_temperature_replacements",
        lambda: replace_source_temperatures_with_usgs(
            base_ss=base_ss,
            hgrid=hgrid,
            source_mapping_dir=relocated_dir,
            start_time=config.startdate,
            usgs_cache_folder=usgs_cache,
            nwm_shapefile=NWM_SHAPEFILE,
            states=STOFS3D_ATL_STATES,
            diagnostics_dir=source_sink_dir / "source_temperature",
            pooling=config.source_temperature_pooling,
            nwm_data_dir=nwm_cache,
        ),
        timings,
    )
    base_ss, counts["selected_source_temperature_overrides"] = _timed(
        "selected_source_temperature_overrides",
        lambda: apply_source_temperature_overrides(
            base_ss=base_ss,
            hgrid=hgrid,
            points=override_points,
            start_time=config.startdate,
            usgs_cache_folder=usgs_cache,
            temperature_pooling=config.source_temperature_pooling,
        ),
        timings,
    )
    base_ss, counts["region_zeroed_sources"] = _timed(
        "region_zeroing",
        lambda: zero_configured_source_regions(
            base_ss=base_ss,
            hgrid=hgrid,
            correction_info=config.zero_source_region_info,
        ),
        timings,
    )
    patch_dir = source_sink_dir / "patch_artificial_island_source_sink"
    patch_dir.mkdir()
    (patch_dir / "hgrid.gr3").symlink_to(hgrid_file)
    patched_ss = _timed(
        "artificial_island_post_patch",
        lambda: apply_artificial_island_corrections(
            base_ss=base_ss,
            hgrid=hgrid,
            original_source_sink_dir=original_dir,
            patch_info_file=config.artificial_island_source_sink_info,
            start_time=config.startdate,
            usgs_cache_folder=usgs_cache,
            output_dir=patch_dir,
        ),
        timings,
    )
    comparisons["post_patch"] = _compare_files(
        patch_dir,
        reference_ss / patch_dir.name,
        TEXT_OUTPUTS,
    )
    source_nc = _compare_source_nc(
        patch_dir / "source.nc",
        reference_ss / patch_dir.name / "source.nc",
    )
    timings["total"] = time.monotonic() - started
    constant_sink_created = (source_sink_dir / "constant_sink").exists()
    passed = (
        all(initial_comparison.values())
        and all(all(stage.values()) for stage in comparisons.values())
        and source_nc["match"]
        and not constant_sink_created
    )
    report = {
        "initial_results": str(initial_results),
        "reference_input": str(reference_input),
        "output": str(output),
        "timings_seconds": timings,
        "initial_stage_byte_identical_to_reference": initial_comparison,
        "operation_counts": counts,
        "final_source_count": len(patched_ss.source_eles),
        "final_sink_count": len(patched_ss.sink_eles),
        "byte_identical_to_reference": comparisons,
        "source_nc_comparison": source_nc,
        "constant_sink_created": constant_sink_created,
        "passed": passed,
    }
    report_file = output / "reproduction_result.json"
    report_file.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2), flush=True)
    return 0 if passed else 2


if __name__ == "__main__":
    raise SystemExit(main())
