#!/usr/bin/env python3

import os
from pathlib import Path
import numpy as np
import xarray as xr

###############################################################################
# EDIT THESE THREE LINES BEFORE RUNNING
###############################################################################

ref_file = Path("/sciclone/schism10/feiye/STOFS3D-v7.3/I21g/NWM_20221201/nwm.t00z.medium_range.channel_rt_1.2023092000.conus.nc")
original_nwm_folder = Path("/sciclone/schism10/feiye/STOFS3D-v7.3/I21g/NWM_20221201/")
output_folder = Path("/sciclone/schism10/feiye/STOFS3D-v7.3/I21g/NWM_20221201_realigned/")

###############################################################################


def load_feature_ids(nc_path):
    ds = xr.open_dataset(nc_path)
    try:
        fids = ds["feature_id"].values
    finally:
        ds.close()
    return fids


def reindex_and_write(src_path, dst_path, target_feature_ids):
    ds = xr.open_dataset(src_path)
    try:
        if "feature_id" not in ds.coords:
            ds = ds.set_coords("feature_id")
        ds_aligned = ds.reindex(feature_id=target_feature_ids)
        ds_aligned.to_netcdf(dst_path)
    finally:
        ds.close()


if __name__ == "__main__":

    # make sure output folder exists
    output_folder.mkdir(parents=True, exist_ok=True)

    # read reference feature_id array
    print(f"Loading reference feature IDs from: {ref_file}")
    ref_ids = load_feature_ids(ref_file)

    # loop over input .nc files
    for src_path in sorted(original_nwm_folder.glob("*.nc")):
        dst_path = output_folder / src_path.name

        # skip if already exists
        if dst_path.exists():
            print(f"[SKIP] {dst_path.name} already exists.")
            continue

        print(f"Processing: {src_path.name}")
        cur_ids = load_feature_ids(src_path)

        same_shape = cur_ids.shape == ref_ids.shape
        same_order = same_shape and np.array_equal(cur_ids, ref_ids)

        if same_order:
            # symlink to save time & disk
            try:
                # Compute path from output file to source file
                rel_target = os.path.relpath(src_path, start=dst_path.parent)
                os.symlink(rel_target, dst_path)
                print(f"  -> identical IDs, relative symlink: {rel_target}")
            except Exception as e:
                print(f"  -> symlink failed ({e}), copying instead.")
                reindex_and_write(src_path, dst_path, ref_ids)
        else:
            print(
                f"  -> feature_id differs "
                f"(src={len(cur_ids)}, ref={len(ref_ids)}), reindexing..."
            )
            reindex_and_write(src_path, dst_path, ref_ids)
            print("  -> new aligned file written.")

    print("\nDone.")
