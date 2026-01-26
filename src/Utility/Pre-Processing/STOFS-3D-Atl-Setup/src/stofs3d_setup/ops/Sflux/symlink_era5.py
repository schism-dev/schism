#!/usr/bin/env python3
from __future__ import annotations
import argparse
from datetime import datetime, timedelta
from pathlib import Path
import os
import sys

TYPES = ("air", "prc", "rad")


def daterange_inclusive(d0: datetime, d1: datetime):
    cur = d0
    while cur <= d1:
        yield cur
        cur = cur + timedelta(days=1)


def make_link(src: Path, dst: Path, force: bool, dry_run: bool, absolute: bool):
    if dst.exists() or dst.is_symlink():
        if force:
            if not dry_run:
                try:
                    dst.unlink()
                except FileNotFoundError:
                    pass
        else:
            print(f"[skip] target exists: {dst}", file=sys.stderr)
            return False

    if not src.exists():
        print(f"[warn] missing source: {src}", file=sys.stderr)
        return False

    if dry_run:
        link_target = str(src.resolve() if absolute else os.path.relpath(src, start=dst.parent))
        print(f"[dry-run] ln -s {link_target} -> {dst}")
    else:
        link_target = src.resolve() if absolute else os.path.relpath(src, start=dst.parent)
        dst.symlink_to(link_target)
        print(f"[link] {dst} -> {link_target}")
    return True


def main():
    """
    Sample usage:
    under the run dir:
    python symlink_era5.py 2020-12-01 2022-01-01 --base-dir {ERA5_dir} --sflux-dir sflux --force --dry-run
    2020-12-01 to 2022-01-01 (inclusive)

    ERA5_dir on sciclone: /sciclone/schism10/hyu05/NOAA_NWM/EnOI/sflux/ERA5/cnv/
    python /sciclone/data10/feiye/stofs3d-setup/src/stofs3d_setup/ops/Sflux/symlink_era5.py \
            2022-12-01 2024-01-01 --base-dir /sciclone/schism10/hyu05/NOAA_NWM/EnOI/sflux/ERA5/cnv/ \
                --sflux-dir sflux --force --absolute --dry-run

    ERA5_dir on Hercules: /work2/noaa/nos-surge/feiye/ERA5/
    python ./symlink_era5.py 2021-12-01 2023-01-01 --base-dir ./ --sflux-dir sflux --force --dry-run
    """

    ap = argparse.ArgumentParser(
        description="Create SCHISM sflux symlinks from monthly ERA5 files."
    )
    ap.add_argument("start_date", help="Start date (YYYY-MM-DD)")
    ap.add_argument("end_date", help="End date (YYYY-MM-DD)")
    ap.add_argument("--base-dir", default=".", help="Root dir containing YYYY_MM folders (default: .)")
    ap.add_argument("--sflux-dir", default="sflux", help="Output sflux directory (default: sflux)")
    ap.add_argument("--force", action="store_true", help="Overwrite existing links")
    ap.add_argument("--dry-run", action="store_true", help="Show what would be done without changing files")
    ap.add_argument("--absolute", action="store_true", help="Use absolute path for symlink targets instead of relative")

    args = ap.parse_args()

    try:
        start = datetime.strptime(args.start_date, "%Y-%m-%d").date()
        end = datetime.strptime(args.end_date, "%Y-%m-%d").date()
    except ValueError:
        print("Dates must be in YYYY-MM-DD format.", file=sys.stderr)
        sys.exit(2)

    if end < start:
        print("end_date must be >= start_date.", file=sys.stderr)
        sys.exit(2)

    # Pad one day on each side
    start_p = start - timedelta(days=1)
    end_p = end + timedelta(days=1)

    base = Path(args.base_dir).resolve()
    sflux = Path(args.sflux_dir).resolve()

    if not args.dry_run:
        sflux.mkdir(parents=True, exist_ok=True)

    # Build the full inclusive date list (sorted)
    dates = list(daterange_inclusive(datetime.combine(start_p, datetime.min.time()),
                                     datetime.combine(end_p, datetime.min.time())))

    print(f"[info] Date window (padded): {start_p} .. {end_p} ({len(dates)} days)")
    print(f"[info] Base dir:  {base}")
    print(f"[info] sflux dir: {sflux}")
    print(f"[info] Types: {', '.join(TYPES)}")
    print(f"[info] Mode: {'dry-run' if args.dry_run else 'apply'}; force={args.force}")

    # Link order: air, then prc, then rad; numbering restarts at 1 for each type
    totals = {t: 0 for t in TYPES}
    for t in TYPES:
        seq = 1
        for d in dates:
            y_m = d.strftime("%Y_%m")
            y_m_d = d.strftime("%Y_%m_%d")
            src = base / y_m / f"era5_{t}.{y_m_d}.nc"
            dst = sflux / f"sflux_{t}_1.{seq}.nc"
            if make_link(src, dst, force=args.force, dry_run=args.dry_run, absolute=args.absolute):
                totals[t] += 1
                seq += 1
        print(f"[summary] {t}: linked {totals[t]} file(s) into {sflux}")
    
    # write a dummy sflux_inputs.txt
    if not args.dry_run:
        with open(sflux / "sflux_inputs.txt", "w") as f:
            f.write("&sflux_inputs\n")
            f.write("/\n")

    # Non-zero exit if nothing linked at all (and not just dry-run)
    if not args.dry_run and sum(totals.values()) == 0:
        print("[error] No links were created. Check paths and filenames.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
