#!/usr/bin/env python3
"""Set the standard STOFS-3D open and land boundaries on a SCHISM grid.

The mesh body is copied byte-for-byte.  Any boundary section already present
after the last element is discarded, then the newly computed ``grd.bnd`` is
appended.  A standalone copy of the same boundary section is always written.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import tempfile
from typing import Mapping, Sequence

import numpy as np

from .constants import STOFS_BOUNDARIES


def convert_boundary_dict(
    boundary_dict: Mapping[str, Sequence[Sequence[float]]],
) -> list[list[float]]:
    """Convert named endpoint pairs to pylib's ``[x1, x2, y1, y2]`` form."""
    boundary_list = []
    for name, points in boundary_dict.items():
        boundary_points = np.asarray(points, dtype=float)
        if boundary_points.shape != (2, 2):
            raise ValueError(
                f"boundary {name!r} must contain two [x, y] endpoint pairs"
            )
        boundary_list.append(
            [
                boundary_points[0, 0],
                boundary_points[1, 0],
                boundary_points[0, 1],
                boundary_points[1, 1],
            ]
        )
    if not boundary_list:
        raise ValueError("at least one open boundary must be specified")
    return boundary_list


def _read_hgrid(path: Path, reader: str = "auto"):
    """Read a grid, preferring pylib's accelerated reader where available."""
    if reader not in {"auto", "cpp", "python"}:
        raise ValueError(f"unknown reader: {reader}")

    if reader in {"auto", "cpp"}:
        try:
            from pylib_experimental.schism_file import cread_schism_hgrid

            print("Reading hgrid with pylib's accelerated C++ reader")
            return cread_schism_hgrid(str(path))
        except (ImportError, OSError) as exc:
            if reader == "cpp":
                raise RuntimeError("the accelerated hgrid reader is unavailable") from exc
            print(f"Accelerated reader unavailable ({exc}); using Python reader")

    from pylib import read_schism_hgrid

    print("Reading hgrid with pylib's Python reader")
    return read_schism_hgrid(str(path))


def compute_boundary_file(
    input_grid: Path,
    boundary_output: Path,
    boundary_dict: Mapping[str, Sequence[Sequence[float]]] = STOFS_BOUNDARIES,
    reader: str = "auto",
) -> None:
    """Compute boundaries and atomically write the standalone boundary file."""
    boundary_output.parent.mkdir(parents=True, exist_ok=True)
    grid = _read_hgrid(input_grid, reader=reader)
    grid.compute_bnd(convert_boundary_dict(boundary_dict))

    fd, temporary_name = tempfile.mkstemp(
        prefix=f".{boundary_output.name}.", dir=boundary_output.parent
    )
    os.close(fd)
    temporary_path = Path(temporary_name)
    try:
        grid.write_bnd(str(temporary_path))
        os.replace(temporary_path, boundary_output)
    finally:
        temporary_path.unlink(missing_ok=True)


def mesh_body_line_count(input_grid: Path) -> tuple[int, int, int]:
    """Return ``(body_lines, ne, np)`` from a SCHISM grid header."""
    with input_grid.open("rb") as stream:
        if not stream.readline():
            raise ValueError(f"empty input grid: {input_grid}")
        counts_line = stream.readline()
    try:
        fields = counts_line.split()
        ne, np_ = int(fields[0]), int(fields[1])
    except (IndexError, ValueError) as exc:
        raise ValueError(
            f"invalid SCHISM element/node count line in {input_grid}: "
            f"{counts_line!r}"
        ) from exc
    if ne < 1 or np_ < 1:
        raise ValueError(f"invalid SCHISM counts in {input_grid}: ne={ne}, np={np_}")
    return 2 + np_ + ne, ne, np_


def replace_boundary_section(
    input_grid: Path, boundary_file: Path, output_grid: Path
) -> tuple[int, int]:
    """Copy the mesh body and replace its optional boundary tail.

    Returns the number of discarded lines and bytes from the old tail.
    """
    body_lines, _, _ = mesh_body_line_count(input_grid)
    output_grid.parent.mkdir(parents=True, exist_ok=True)

    fd, temporary_name = tempfile.mkstemp(
        prefix=f".{output_grid.name}.", dir=output_grid.parent
    )
    os.close(fd)
    temporary_path = Path(temporary_name)
    discarded_lines = 0
    discarded_bytes = 0
    try:
        with input_grid.open("rb") as source, temporary_path.open("wb") as target:
            for line_number in range(1, body_lines + 1):
                line = source.readline()
                if not line:
                    raise ValueError(
                        f"{input_grid} ended at line {line_number - 1}; "
                        f"header requires {body_lines} mesh-body lines"
                    )
                target.write(line)

            for line in source:
                discarded_lines += 1
                discarded_bytes += len(line)

            with boundary_file.open("rb") as boundary_stream:
                first = boundary_stream.read(1)
                if not first:
                    raise ValueError(f"empty boundary file: {boundary_file}")
                target.write(first)
                while True:
                    chunk = boundary_stream.read(1024 * 1024)
                    if not chunk:
                        break
                    target.write(chunk)

        os.replace(temporary_path, output_grid)
    finally:
        temporary_path.unlink(missing_ok=True)
    return discarded_lines, discarded_bytes


def make_stofsv8_boundary(
    hgrid_obj,
    output_dir: str | Path = "./",
    write_hgrid: bool = False,
) -> None:
    """Compatibility wrapper for callers of the original script."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    hgrid_obj.compute_bnd(convert_boundary_dict(STOFS_BOUNDARIES))
    hgrid_obj.write_bnd(str(output_dir / "grd.bnd"))
    if write_hgrid:
        hgrid_obj.save(str(output_dir / "hgrid_with_bnd.gr3"), fmt=1)


def command_line_interface(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Set the standard STOFS-3D open/land boundaries, replacing any "
            "existing hgrid boundary section without rewriting the mesh body."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("grid", type=Path, help="input SCHISM hgrid (.gr3 or .ll)")
    parser.add_argument(
        "-o", "--output-grid", type=Path,
        help="complete output grid (default: hgrid_with_bnd.gr3 beside GRID)",
    )
    parser.add_argument(
        "-b", "--boundary-output", type=Path,
        help="standalone boundary file (default: grd.bnd beside output grid)",
    )
    parser.add_argument(
        "--in-place", action="store_true",
        help="replace GRID atomically; cannot be used when GRID is a symlink",
    )
    parser.add_argument(
        "--reader", choices=("auto", "cpp", "python"), default="auto",
        help="pylib grid reader",
    )
    args = parser.parse_args(argv)

    args.grid = args.grid.expanduser().absolute()
    if not args.grid.is_file():
        parser.error(f"input grid does not exist or is not a file: {args.grid}")
    if args.grid.suffix.lower() not in {".gr3", ".ll"}:
        parser.error("GRID must have a .gr3 or .ll suffix")
    if args.in_place and args.output_grid is not None:
        parser.error("--in-place and --output-grid are mutually exclusive")
    if args.in_place and args.grid.is_symlink():
        parser.error("--in-place is refused for a symlink; name its target explicitly")

    if args.in_place:
        args.output_grid = args.grid
    elif args.output_grid is None:
        args.output_grid = args.grid.parent / "hgrid_with_bnd.gr3"
    else:
        args.output_grid = args.output_grid.expanduser().absolute()

    if args.boundary_output is None:
        args.boundary_output = args.output_grid.parent / "grd.bnd"
    else:
        args.boundary_output = args.boundary_output.expanduser().absolute()

    def aliases(first: Path, second: Path) -> bool:
        if first == second:
            return True
        try:
            return first.exists() and second.exists() and os.path.samefile(first, second)
        except OSError:
            return False

    if not args.in_place and aliases(args.output_grid, args.grid):
        parser.error("--output-grid aliases GRID; use --in-place explicitly")
    if aliases(args.boundary_output, args.grid):
        parser.error("--boundary-output must not alias GRID")
    if aliases(args.boundary_output, args.output_grid):
        parser.error("the standalone boundary file and complete grid must differ")
    return args


def main(argv: Sequence[str] | None = None) -> int:
    args = command_line_interface(argv)
    _, ne, np_ = mesh_body_line_count(args.grid)
    print(f"Input: {args.grid} ({np_} nodes, {ne} elements)")
    compute_boundary_file(
        args.grid, args.boundary_output, reader=args.reader
    )
    discarded_lines, discarded_bytes = replace_boundary_section(
        args.grid, args.boundary_output, args.output_grid
    )
    print(f"Boundary file: {args.boundary_output}")
    print(f"Complete grid: {args.output_grid}")
    if discarded_lines:
        print(
            f"Replaced existing boundary tail: {discarded_lines} lines, "
            f"{discarded_bytes} bytes"
        )
    else:
        print("Input had no existing boundary section; appended the new boundary")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
