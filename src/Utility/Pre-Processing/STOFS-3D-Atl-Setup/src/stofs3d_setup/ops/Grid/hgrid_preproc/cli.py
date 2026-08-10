#!/usr/bin/env python3
"""Run STOFS-3D hgrid preprocessing and prepare offline partition checking.

Offline partitioning itself remains a separate, explicitly launched step.
"""

from __future__ import annotations

import argparse
import gc
from pathlib import Path
import sys
from typing import Sequence

from .partition_check import prepare_partition_check as create_partition_check


STAGES = ("improve", "split", "boundary", "partition")
DEFAULT_PARTITION_SCRIPT = Path("~/bin/partition_offline.pl")


def _refuse_existing(paths: Sequence[Path], force: bool) -> None:
    existing = [path for path in paths if path.exists() or path.is_symlink()]
    if existing and not force:
        joined = "\n  ".join(str(path) for path in existing)
        raise FileExistsError(
            "stage output already exists; use --force to replace it:\n  " + joined
        )


def run_improvement(args: argparse.Namespace, output_dir: Path) -> Path:
    print("\n=== Stage 1/4: improve grid ===", flush=True)
    _refuse_existing(
        [output_dir / f"hgrid{suffix}" for suffix in (".2dm", ".gr3", ".ll")],
        args.force,
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    from .improve import main as improve_main

    stage_args = [
        str(args.grid),
        "--input-crs", str(args.input_crs),
        "--skewness-threshold", str(args.skewness_threshold),
        "--area-threshold", str(args.area_threshold),
        "--max-rounds", str(args.max_rounds),
        "--intersection-fix-rounds", str(args.intersection_fix_rounds),
        "--output-dir", str(output_dir),
        "--output-stem", "hgrid",
    ]
    if args.diagnostic_iterations:
        stage_args.append("--diagnostic-iterations")
    improve_main(stage_args)

    output = output_dir / "hgrid.ll"
    if not output.is_file():
        raise RuntimeError(f"improvement did not create expected output: {output}")
    gc.collect()
    return output


def run_quad_split(args: argparse.Namespace, input_grid: Path, output_dir: Path) -> Path:
    print("\n=== Stage 2/4: split bad quads ===", flush=True)
    _refuse_existing(
        [output_dir / "hgrid.gr3.new", output_dir / "split_loc.bp"], args.force
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    from .bad_quads import fix_bad_quads
    from .improve import read_schism_hgrid

    grid = read_schism_hgrid(str(input_grid))
    fix_bad_quads(
        grid,
        angle_ratio=args.quad_angle_ratio,
        aspect_ratio=args.quad_aspect_ratio,
        output_dir=output_dir,
        legacy_header=True,
    )
    del grid
    gc.collect()

    output = output_dir / "hgrid.gr3.new"
    if not output.is_file():
        raise RuntimeError(f"quad splitting did not create expected output: {output}")
    return output


def run_boundaries(args: argparse.Namespace, input_grid: Path, output_dir: Path) -> Path:
    print("\n=== Stage 3/4: set open and land boundaries ===", flush=True)
    boundary_file = output_dir / "grd.bnd"
    output_grid = output_dir / "hgrid_with_bnd.gr3"
    _refuse_existing([boundary_file, output_grid], args.force)
    output_dir.mkdir(parents=True, exist_ok=True)

    from .boundaries import compute_boundary_file, replace_boundary_section

    compute_boundary_file(
        input_grid, boundary_file, reader=args.boundary_reader
    )
    discarded_lines, discarded_bytes = replace_boundary_section(
        input_grid, boundary_file, output_grid
    )
    if discarded_lines:
        print(
            f"Replaced an existing boundary tail: {discarded_lines} lines, "
            f"{discarded_bytes} bytes"
        )
    if not output_grid.is_file():
        raise RuntimeError(f"boundary stage did not create expected output: {output_grid}")
    gc.collect()
    return output_grid


def prepare_partition_check(
    args: argparse.Namespace, input_grid: Path, output_dir: Path
) -> None:
    print("\n=== Stage 4/4: prepare standalone partition check ===", flush=True)
    create_partition_check(
        input_grid=input_grid,
        output_dir=output_dir,
        partition_script=args.partition_script,
        partition_count=args.partition_count,
        overwrite=args.force,
    )


def command_line_interface(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Improve a STOFS-3D grid, split bad quads, set boundaries, and "
            "prepare—but do not run—the offline partition check."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("grid", type=Path, help="input grid (.2dm, .gr3, or .ll)")
    parser.add_argument(
        "--output-root", type=Path, required=True,
        help="parent of Improve, Split_quads, Bnd, and Partition_check",
    )
    parser.add_argument("--input-crs", default="esri:102008")
    parser.add_argument("--skewness-threshold", type=float, default=35.0)
    parser.add_argument("--area-threshold", type=float, default=1.0)
    parser.add_argument("--max-rounds", type=int, default=4)
    parser.add_argument("--intersection-fix-rounds", type=int, default=0)
    parser.add_argument("--diagnostic-iterations", action="store_true")
    parser.add_argument("--quad-angle-ratio", type=float, default=0.5)
    parser.add_argument("--quad-aspect-ratio", type=float, default=5.0)
    parser.add_argument(
        "--boundary-reader", choices=("auto", "cpp", "python"), default="auto"
    )
    parser.add_argument("--partition-count", type=int, default=3994)
    parser.add_argument(
        "--partition-script", type=Path, default=DEFAULT_PARTITION_SCRIPT,
        help="standalone partition driver printed in the final instructions",
    )
    parser.add_argument(
        "--start-at", choices=STAGES, default="improve",
        help="resume at this stage using products in the standard directories",
    )
    parser.add_argument(
        "--stop-after", choices=STAGES, default="partition",
        help="stop after this stage",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="replace expected stage products; unrelated files are preserved",
    )
    args = parser.parse_args(argv)

    args.grid = args.grid.expanduser().absolute()
    args.output_root = args.output_root.expanduser().absolute()
    if not args.grid.is_file():
        parser.error(f"input grid does not exist or is not a file: {args.grid}")
    if args.grid.suffix.lower() not in {".2dm", ".gr3", ".ll"}:
        parser.error("GRID must have a .2dm, .gr3, or .ll suffix")
    if args.skewness_threshold <= 0 or args.area_threshold <= 0:
        parser.error("quality thresholds must be positive")
    if args.max_rounds < 0 or args.intersection_fix_rounds < 0:
        parser.error("repair-round counts cannot be negative")
    if args.quad_angle_ratio <= 0 or args.quad_aspect_ratio <= 0:
        parser.error("quad thresholds must be positive")
    if args.partition_count < 1:
        parser.error("--partition-count must be positive")
    if STAGES.index(args.start_at) > STAGES.index(args.stop_after):
        parser.error("--start-at cannot follow --stop-after")
    return args


def main(argv: Sequence[str] | None = None) -> int:
    args = command_line_interface(argv)
    directories = {
        "improve": args.output_root / "Improve",
        "split": args.output_root / "Split_quads",
        "boundary": args.output_root / "Bnd",
        "partition": args.output_root / "Partition_check",
    }
    products = {
        "improve": directories["improve"] / "hgrid.ll",
        "split": directories["split"] / "hgrid.gr3.new",
        "boundary": directories["boundary"] / "hgrid_with_bnd.gr3",
    }

    first = STAGES.index(args.start_at)
    last = STAGES.index(args.stop_after)
    selected = STAGES[first:last + 1]
    print(f"Input grid: {args.grid}")
    print(f"Output root: {args.output_root}")
    print(f"Stages: {', '.join(selected)}")

    if "improve" in selected:
        products["improve"] = run_improvement(args, directories["improve"])

    if "split" in selected:
        if not products["improve"].is_file():
            raise FileNotFoundError(f"missing prior-stage product: {products['improve']}")
        products["split"] = run_quad_split(
            args, products["improve"], directories["split"]
        )

    if "boundary" in selected:
        if not products["split"].is_file():
            raise FileNotFoundError(f"missing prior-stage product: {products['split']}")
        products["boundary"] = run_boundaries(
            args, products["split"], directories["boundary"]
        )

    if "partition" in selected:
        if not products["boundary"].is_file():
            raise FileNotFoundError(f"missing prior-stage product: {products['boundary']}")
        prepare_partition_check(
            args, products["boundary"], directories["partition"]
        )

    print("\nRequested preprocessing stages completed successfully.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileExistsError, FileNotFoundError, RuntimeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
