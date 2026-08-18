"""Convert SCHISM ``.reg`` polygon files to a GeoPackage.

The converter intentionally does not repair overlapping regions.  Use
``find_polygon_overlaps`` or ``assert_non_overlapping`` to validate a completed
collection before using it for order-independent regional operations.
"""

from __future__ import annotations

import argparse
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
import tempfile

import geopandas as gpd
import numpy as np
from shapely.geometry import MultiPolygon, Polygon
from shapely.validation import explain_validity


DEFAULT_CRS = "EPSG:4326"
DEFAULT_LAYER = "regional_tweaks"
APPLY_ORDER_COLUMN = "apply_order"


class RegFileError(ValueError):
    """Raised when a SCHISM region file cannot be converted to a polygon."""


class OverlappingRegionsError(ValueError):
    """Raised when two or more region polygons overlap in area."""


@dataclass(frozen=True)
class PolygonOverlap:
    """A positive-area intersection between two region features."""

    first_index: object
    first_region: str
    second_index: object
    second_region: str
    area: float


def _nonempty_lines(path: Path) -> list[tuple[int, str]]:
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as exc:
        raise RegFileError(f"Could not read region file {path}: {exc}") from exc
    return [
        (line_number, line.strip())
        for line_number, line in enumerate(text.splitlines(), start=1)
        if line.strip()
    ]


def _parse_integer(value: str, path: Path, line_number: int, label: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise RegFileError(
            f"{path}:{line_number}: invalid {label} {value!r}"
        ) from exc
    if parsed < 1:
        raise RegFileError(
            f"{path}:{line_number}: {label} must be positive, got {parsed}"
        )
    return parsed


def read_reg_polygon(reg_file: str | Path) -> Polygon | MultiPolygon:
    """Read polygon geometry from an ACE/SMS SCHISM ``.reg`` file.

    Files containing more than one polygon block are returned as a
    :class:`~shapely.geometry.MultiPolygon`.  Invalid, empty, and zero-area
    polygons are rejected rather than silently repaired.
    """

    path = Path(reg_file)
    lines = _nonempty_lines(path)
    if len(lines) < 4:
        raise RegFileError(f"{path}: incomplete region file")

    # The first non-empty line is a free-form title.
    region_count_line, region_count_text = lines[1]
    region_count = _parse_integer(
        region_count_text.split()[0], path, region_count_line, "polygon count"
    )
    cursor = 2
    polygons: list[Polygon] = []

    for polygon_number in range(1, region_count + 1):
        if cursor >= len(lines):
            raise RegFileError(
                f"{path}: missing header for polygon {polygon_number}"
            )
        header_line, header = lines[cursor]
        cursor += 1
        fields = header.split()
        if not fields:
            raise RegFileError(
                f"{path}:{header_line}: missing polygon point count"
            )
        point_count = _parse_integer(
            fields[0], path, header_line, "polygon point count"
        )
        if point_count < 3:
            raise RegFileError(
                f"{path}:{header_line}: polygon {polygon_number} has only "
                f"{point_count} points"
            )
        if cursor + point_count > len(lines):
            raise RegFileError(
                f"{path}:{header_line}: polygon {polygon_number} declares "
                f"{point_count} points but the file ends early"
            )

        coordinates: list[tuple[float, float]] = []
        for line_number, coordinate_text in lines[cursor : cursor + point_count]:
            coordinate_fields = coordinate_text.split()
            if len(coordinate_fields) < 2:
                raise RegFileError(
                    f"{path}:{line_number}: expected x and y coordinates"
                )
            try:
                x, y = map(float, coordinate_fields[:2])
            except ValueError as exc:
                raise RegFileError(
                    f"{path}:{line_number}: invalid coordinates "
                    f"{coordinate_fields[:2]!r}"
                ) from exc
            if not np.isfinite(x) or not np.isfinite(y):
                raise RegFileError(
                    f"{path}:{line_number}: coordinates must be finite"
                )
            coordinates.append((x, y))
        cursor += point_count

        if len(set(coordinates)) < 3:
            raise RegFileError(
                f"{path}:{header_line}: polygon {polygon_number} has fewer "
                "than three distinct points"
            )
        polygon = Polygon(coordinates)
        if not polygon.is_valid:
            raise RegFileError(
                f"{path}:{header_line}: polygon {polygon_number} is invalid: "
                f"{explain_validity(polygon)}"
            )
        if polygon.is_empty or polygon.area == 0:
            raise RegFileError(
                f"{path}:{header_line}: polygon {polygon_number} has zero area"
            )
        polygons.append(polygon)

    if cursor != len(lines):
        line_number, _ = lines[cursor]
        raise RegFileError(
            f"{path}:{line_number}: unexpected content after the declared polygons"
        )

    geometry: Polygon | MultiPolygon
    if len(polygons) == 1:
        geometry = polygons[0]
    else:
        geometry = MultiPolygon(polygons)
        if not geometry.is_valid:
            raise RegFileError(
                f"{path}: polygon collection is invalid: "
                f"{explain_validity(geometry)}"
            )
    return geometry


def _resolve_reg_files(
    reg_files: str | Path | Iterable[str | Path],
) -> list[Path]:
    if isinstance(reg_files, (str, Path)):
        path = Path(reg_files)
        files = sorted(path.glob("*.reg")) if path.is_dir() else [path]
    else:
        files = [Path(path) for path in reg_files]

    if not files:
        raise ValueError("No .reg files were supplied")
    missing = [path for path in files if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "Region files do not exist: " + ", ".join(map(str, missing))
        )
    wrong_suffix = [path for path in files if path.suffix.lower() != ".reg"]
    if wrong_suffix:
        raise ValueError(
            "Expected .reg files, got: " + ", ".join(map(str, wrong_suffix))
        )

    names = [path.stem for path in files]
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise ValueError(
            "Region names must be unique; duplicate file stems: "
            + ", ".join(duplicates)
        )
    return files


def reg_files_to_geodataframe(
    reg_files: str | Path | Iterable[str | Path],
    *,
    crs: str = DEFAULT_CRS,
    attributes: Mapping[str, Mapping[str, object]] | None = None,
) -> gpd.GeoDataFrame:
    """Convert one or more ``.reg`` files into an in-memory GeoDataFrame.

    ``attributes`` is keyed by region name (the input filename stem).  When it
    is supplied, it must contain exactly one attribute mapping per input file.
    ``apply_order`` is assigned from one in the supplied file order; directory
    inputs are sorted by filename first.
    """

    files = _resolve_reg_files(reg_files)
    regions = [path.stem for path in files]
    if attributes is not None:
        expected = set(regions)
        supplied = set(attributes)
        missing = sorted(expected - supplied)
        extra = sorted(supplied - expected)
        if missing or extra:
            details = []
            if missing:
                details.append("missing: " + ", ".join(missing))
            if extra:
                details.append("unknown: " + ", ".join(extra))
            raise ValueError(
                "Attribute regions do not match inputs ("
                + "; ".join(details)
                + ")"
            )

    records = []
    for apply_order, (path, region) in enumerate(zip(files, regions), start=1):
        record: dict[str, object] = {
            "region": region,
            APPLY_ORDER_COLUMN: apply_order,
            "source_file": path.name,
            "geometry": read_reg_polygon(path),
        }
        if attributes is not None:
            reserved = {
                "region",
                APPLY_ORDER_COLUMN,
                "source_file",
                "geometry",
            } & set(attributes[region])
            if reserved:
                raise ValueError(
                    f"Attributes for {region!r} use reserved columns: "
                    + ", ".join(sorted(reserved))
                )
            record.update(attributes[region])
        records.append(record)
    return gpd.GeoDataFrame(records, geometry="geometry", crs=crs)


def find_polygon_overlaps(
    regions: gpd.GeoDataFrame,
    *,
    region_column: str = "region",
    area_tolerance: float = 0.0,
) -> list[PolygonOverlap]:
    """Return every pair of polygons whose intersection exceeds *area_tolerance*.

    Areas are expressed in the GeoDataFrame CRS units squared.  Shared edges
    and vertices have zero area and are not considered overlaps.
    """

    if area_tolerance < 0:
        raise ValueError("area_tolerance must be non-negative")
    if region_column not in regions:
        raise ValueError(f"GeoDataFrame has no {region_column!r} column")

    overlaps: list[PolygonOverlap] = []
    geometries = regions.geometry
    for position, (first_index, first_geometry) in enumerate(geometries.items()):
        if first_geometry is None or first_geometry.is_empty:
            continue
        candidate_positions = regions.sindex.query(
            first_geometry, predicate="intersects"
        )
        for candidate_position in candidate_positions:
            if int(candidate_position) <= position:
                continue
            second_index = regions.index[int(candidate_position)]
            second_geometry = geometries.iloc[int(candidate_position)]
            area = float(first_geometry.intersection(second_geometry).area)
            if area <= area_tolerance:
                continue
            overlaps.append(
                PolygonOverlap(
                    first_index=first_index,
                    first_region=str(regions.at[first_index, region_column]),
                    second_index=second_index,
                    second_region=str(regions.at[second_index, region_column]),
                    area=area,
                )
            )
    return overlaps


def assert_non_overlapping(
    regions: gpd.GeoDataFrame,
    *,
    region_column: str = "region",
    area_tolerance: float = 0.0,
) -> None:
    """Raise :class:`OverlappingRegionsError` for positive-area overlaps."""

    overlaps = find_polygon_overlaps(
        regions,
        region_column=region_column,
        area_tolerance=area_tolerance,
    )
    if not overlaps:
        return
    conflicts = "; ".join(
        f"{overlap.first_region!r} / {overlap.second_region!r} "
        f"(area={overlap.area:g})"
        for overlap in overlaps
    )
    raise OverlappingRegionsError(f"Overlapping region polygons: {conflicts}")


def reg2gpkg(
    reg_files: str | Path | Iterable[str | Path],
    output_file: str | Path,
    *,
    crs: str = DEFAULT_CRS,
    layer: str = DEFAULT_LAYER,
    attributes: Mapping[str, Mapping[str, object]] | None = None,
    overwrite: bool = False,
) -> gpd.GeoDataFrame:
    """Convert ``.reg`` files to a single GeoPackage layer.

    The output is staged beside its destination and moved into place only after
    GeoPandas closes the completed GeoPackage.
    """

    output_path = Path(output_file)
    if output_path.suffix.lower() != ".gpkg":
        raise ValueError(f"GeoPackage output must end in .gpkg: {output_path}")
    if output_path.exists() and not overwrite:
        raise FileExistsError(
            f"Output already exists: {output_path}; pass overwrite=True to replace it"
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)

    regions = reg_files_to_geodataframe(
        reg_files,
        crs=crs,
        attributes=attributes,
    )
    with tempfile.TemporaryDirectory(
        prefix=f".{output_path.stem}.", dir=output_path.parent
    ) as staging_dir:
        staged_path = Path(staging_dir) / output_path.name
        regions.to_file(staged_path, layer=layer, driver="GPKG", index=False)
        staged_path.replace(output_path)
    return regions


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "inputs",
        nargs="+",
        type=Path,
        help="One or more .reg files, or one directory containing .reg files",
    )
    parser.add_argument("--output", "-o", required=True, type=Path)
    parser.add_argument("--crs", default=DEFAULT_CRS)
    parser.add_argument("--layer", default=DEFAULT_LAYER)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--check-overlap",
        action="store_true",
        help="Fail without writing when any polygons overlap in area",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    inputs: Path | list[Path]
    inputs = args.inputs[0] if len(args.inputs) == 1 else args.inputs
    regions = reg_files_to_geodataframe(inputs, crs=args.crs)
    if args.check_overlap:
        assert_non_overlapping(regions)
    reg2gpkg(
        inputs,
        args.output,
        crs=args.crs,
        layer=args.layer,
        overwrite=args.overwrite,
    )
    print(f"Wrote {len(regions)} regions to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
