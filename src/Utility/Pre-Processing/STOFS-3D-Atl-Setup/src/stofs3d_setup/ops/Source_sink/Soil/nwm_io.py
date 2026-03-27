"""I/O helpers for NWM LDASOUT + soil parameter + WRF input files (Noah-MP style).

This module reads NWM LDASOUT, soil parameter, and wrfinput files, harmonizes units
(m, s, m s⁻¹), selects layers within a target control depth, and outputs
inputs for the `physics` infiltration partitioning module.

Functions:
  - open_ldasout(path_pattern): read precipitation, soil, and snow variables.
  - open_soilparams(path): read static soil parameters.
  - open_wrfinput(path): read DZS layer thickness.
  - prepare_control_depth_inputs(ds_ldas, ds_soil, ds_wrf, depth): prepare arrays.
  - compute_from_nwm(ldas_pattern, soil_param_path, wrfinput_path, depth): high-level wrapper.
"""
from __future__ import annotations

import numpy as np
import xarray as xr
from typing import Dict
import functools
from glob import glob

from physics import SoilParams, compute_fluxes

# -----------------------------------------------------------------------------
# Variable alias mapping
# -----------------------------------------------------------------------------
ALIAS = {
    "QRAIN": ["QRAIN", "RAINRATE"],
    "ACSNOM": ["ACSNOM", "FSNO_ACCUM"],
    "UGDRNOFF": ["UGDRNOFF"],
    "SOIL_M": ["SOIL_M", "SMOIS"],
    "SOIL_W": ["SOIL_W"],
    "SNEQV": ["SNEQV"],
    "DZS": ["DZS", "DZSOIL"],
    "smcmax": ["smcmax", "SMCMAX"],
    "smcref": ["smcref", "SMCREF"],
    "smcwlt": ["smcwlt", "SMCWLT"],
    "bexp": ["bexp", "BEXP"],
    "dksat": ["dksat", "DKSAT"],
    "imperv": ["imperv", "IMPERV"],
}


# -----------------------------------------------------------------------------
# Utility functions
# -----------------------------------------------------------------------------

def _find_var(ds: xr.Dataset, key: str) -> str:
    for name in ALIAS.get(key, [key]):
        if name in ds.variables:
            return name
    raise KeyError(f"Variable {key} not found in dataset.")


def _to_m(da: xr.DataArray) -> xr.DataArray:
    units = (da.attrs.get("units", "") or "").lower()
    if "kg m" in units:
        return da / 1000.0
    return da


def _to_rate(da: xr.DataArray) -> xr.DataArray:
    units = (da.attrs.get("units", "") or "").lower()
    if "kg m" in units:
        return da / 1000.0
    return da


# -----------------------------------------------------------------------------
# File readers
# -----------------------------------------------------------------------------

PROJ4_LDAS = ("+proj=lcc +units=m +a=6370000.0 +b=6370000.0 "
              "+lat_1=30.0 +lat_2=60.0 +lat_0=40.0 +lon_0=-97.0 "
              "+x_0=0 +y_0=0 +k_0=1.0 +nadgrids=@null +wktext  +no_defs")


def lonlat_to_ldas_xy(lon, lat, proj4_str=None):
    """
    lon, lat: arrays in degrees (EPSG:4326)
    proj4_str: e.g., ds.proj4 from LDASOUT
    returns x, y in meters (LDAS projection)
    """
    from pyproj import CRS, Transformer
    if proj4_str is None:
        proj4_str = PROJ4_LDAS

    crs_wgs84 = CRS.from_epsg(4326)
    crs_ldas = CRS.from_proj4(proj4_str)
    tf = Transformer.from_crs(crs_wgs84, crs_ldas, always_xy=True)
    x, y = tf.transform(lon, lat)   # returns in meters
    return np.asarray(x), np.asarray(y)


def _find_xy_names(ds: xr.Dataset):
    xname = next((d for d in ("x", "west_east") if d in ds.dims), None)
    yname = next((d for d in ("y", "south_north") if d in ds.dims), None)
    if xname is None or yname is None:
        raise ValueError(f"Could not find x/y dims in {list(ds.dims)}")
    return xname, yname


def _slice_val(coord_vals, v0, v1):
    """Value-based slice that respects ascending or descending coords."""
    lo, hi = (min(v0, v1), max(v0, v1))
    asc = bool(coord_vals[0] <= coord_vals[-1])
    return slice(lo, hi) if asc else slice(hi, lo)


def _subset_xy(ds, x_rng=None, y_rng=None):
    if x_rng is None or y_rng is None:
        return ds
    # pick dim names that exist in this file
    xname = "x" if "x" in ds.dims else ("west_east" if "west_east" in ds.dims else None)
    yname = "y" if "y" in ds.dims else ("south_north" if "south_north" in ds.dims else None)
    if xname and yname:
        # Note: many LDAS grids have ascending x, ascending y; if yours has descending y, swap the slice bounds
        x0, x1 = x_rng
        y0, y1 = y_rng
        ds = ds.sel({xname: slice(x0, x1), yname: slice(y0, y1)})
    if ds.proj4 != PROJ4_LDAS:
        raise ValueError(f"Unexpected projection: {ds.proj4}")
    return ds


def _find_ldas_xy_names(ds: xr.Dataset):
    xdim = next((d for d in ("x", "west_east") if d in ds.dims), None)
    ydim = next((d for d in ("y", "south_north") if d in ds.dims), None)
    if xdim is None or ydim is None:
        raise ValueError(f"LDAS: could not find x/y dims in {list(ds.dims)}")
    # coord variables with same name as dims are expected in LDAS
    xcoord = ds.coords.get(xdim, None)
    ycoord = ds.coords.get(ydim, None)
    if xcoord is None or ycoord is None:
        raise ValueError("LDAS: expected 1D coordinate variables for x/y.")
    return xdim, ydim, xcoord, ycoord


def _index_slice_for_range(coord: xr.DataArray, v0: float, v1: float) -> slice:
    vals = np.asarray(coord.values)
    lo, hi = (min(v0, v1), max(v0, v1))
    # clamp to data range (in case bbox extends past edges)
    lo = max(lo, float(vals.min()))
    hi = min(hi, float(vals.max()))
    if hi <= lo:
        raise ValueError(f"No overlap: requested [{v0}, {v1}] vs coord [{vals.min()}, {vals.max()}]")
    if vals[0] <= vals[-1]:  # ascending
        i0 = np.searchsorted(vals, lo, side="left")
        i1 = np.searchsorted(vals, hi, side="right")
        return slice(i0, i1)
    else:  # descending
        vrev = vals[::-1]
        i0d = np.searchsorted(vrev, lo, side="left")
        i1d = np.searchsorted(vrev, hi, side="right")
        return slice(len(vals) - i1d, len(vals) - i0d)


def compute_ldas_xy_slices(sample_path: str, x_rng, y_rng):
    """
    Open a single LDASOUT file, compute index slices for (x_rng, y_rng) in projected meters.
    Returns dict: {'xdim', 'ydim', 'x_slice', 'y_slice'}
    """
    with xr.open_dataset(sample_path, decode_times=False, engine="netcdf4") as ds0:
        xdim, ydim, xcoord, ycoord = _find_ldas_xy_names(ds0)
        xs = _index_slice_for_range(xcoord, x_rng[0], x_rng[1])
        ys = _index_slice_for_range(ycoord, y_rng[0], y_rng[1])
    return {"xdim": xdim, "ydim": ydim, "x_slice": xs, "y_slice": ys}


def open_ldasout(path_pattern: str, x_rng=None, y_rng=None, chunks=None) -> xr.Dataset:
    preprocess = functools.partial(_subset_xy, x_rng=x_rng, y_rng=y_rng)
    ds = xr.open_mfdataset(
        path_pattern,
        combine="by_coords",
        decode_times=True,
        preprocess=preprocess,   # <- subsetting happens before dask graph builds
        chunks=chunks,           # e.g., {'time': 8, 'y': 512, 'x': 512}
        engine="netcdf4",
    )
    keep = []
    for k in ["QRAIN", "ACSNOM", "UGDRNOFF", "SOIL_M", "SOIL_W", "SNEQV"]:
        try:
            keep.append(_find_var(ds, k))
        except KeyError:
            pass
    ds = ds[keep]

    # ALSO compute index slices once (from any one file in the pattern)
    # so we can isel() static datasets that lack coords:
    sample_file = glob(path_pattern)[0]
    if x_rng is not None and y_rng is not None:
        # take the first file in the pattern (xarray stores it internally)
        # safest: pass an explicit sample path you already have; otherwise:
        if sample_file:
            xy_slices = compute_ldas_xy_slices(sample_file, x_rng, y_rng)
        else:
            raise FileNotFoundError("Could not determine sample file from dataset for computing xy slices.")
    else:
        xy_slices = None

    return ds, xy_slices


def _apply_xy_isel(ds: xr.Dataset, xy_slices: dict) -> xr.Dataset:
    if xy_slices is None:
        return ds
    xdim = next((d for d in (xy_slices['xdim'], "x", "west_east") if d in ds.dims), None)
    ydim = next((d for d in (xy_slices['ydim'], "y", "south_north") if d in ds.dims), None)
    if xdim is None or ydim is None:
        # dataset truly lacks spatial dims—just return as is
        return ds
    return ds.isel({xdim: xy_slices["x_slice"], ydim: xy_slices["y_slice"]})


def open_soilparams(path: str, xy_slices: dict | None = None, chunks=None) -> xr.Dataset:
    ds = xr.open_dataset(path, chunks=chunks or {'south_north': 512, 'west_east': 512}, engine='netcdf4')
    ds = _apply_xy_isel(ds, xy_slices)
    keep = []
    for k in ["smcmax", "smcref", "smcwlt", "bexp", "dksat", "imperv"]:
        try:
            keep.append(_find_var(ds, k))
        except KeyError:
            pass
    return ds[keep]


def open_wrfinput(path: str, xy_slices: dict | None = None, chunks=None) -> xr.Dataset:
    ds = xr.open_dataset(path, chunks=chunks or {'south_north': 512, 'west_east': 512}, engine='netcdf4')
    ds = _apply_xy_isel(ds, xy_slices)
    if not any(v in ds.variables for v in ALIAS["DZS"]):
        raise KeyError("No DZS variable found in WRF input file.")
    return ds[[_find_var(ds, "DZS")]]


# -----------------------------------------------------------------------------
# Layer selection
# -----------------------------------------------------------------------------

def select_layers_by_depth(DZS: xr.DataArray, depth: float, layer_dim='soil_layers_stag') -> slice:
    if layer_dim not in DZS.dims:
        raise ValueError(f"Expected layer dimension '{layer_dim}' in DZS data array.")

    dz_mean = DZS.mean(dim=[d for d in DZS.dims if d != layer_dim])
    csum = np.array(dz_mean.cumsum(dim=layer_dim))
    nlayer = int((csum < depth).sum())
    nlayer = max(1, min(nlayer + 1, DZS.sizes[layer_dim]))
    return slice(0, nlayer)


# -----------------------------------------------------------------------------
# Data preparation
# -----------------------------------------------------------------------------

def prepare_control_depth_inputs(ds_ldas: xr.Dataset, ds_soil: xr.Dataset, ds_wrf: xr.Dataset, depth: float = 0.4):
    """Prepare LDASOUT + soil param + WRF input data at specified control depth."""

    v = {k: _find_var(ds_ldas, k) for k in ["QRAIN", "ACSNOM", "UGDRNOFF", "SOIL_M", "SOIL_W", "SNEQV"]}
    t = ds_ldas.indexes["time"]
    dt = np.diff(t.values).astype("timedelta64[s]").astype(float)
    dt = np.r_[dt[0], dt]

    QRAIN = _to_rate(ds_ldas[v["QRAIN"]])
    ACSNOM = _to_m(ds_ldas[v["ACSNOM"]])
    UGDRNOFF = _to_m(ds_ldas[v["UGDRNOFF"]])
    SNEQV = _to_m(ds_ldas[v["SNEQV"]])
    SOIL_M = ds_ldas[v["SOIL_M"]]
    SOIL_W = ds_ldas[v["SOIL_W"]]

    # DZS now comes from wrfinput
    DZS = ds_wrf[_find_var(ds_wrf, "DZS")]

    layer_dim = [d for d in DZS.dims if d not in QRAIN.dims and d.lower() != "time"][0]
    layer_slice = select_layers_by_depth(DZS, depth)

    SM_ctl = SOIL_M.isel({layer_dim: layer_slice})
    SW_ctl = SOIL_W.isel({layer_dim: layer_slice})
    DZS_ctl = DZS.isel({layer_dim: layer_slice})

    soil_vars = {k: _find_var(ds_soil, k) for k in ["smcmax", "dksat", "imperv"]}
    smcmax = ds_soil[soil_vars["smcmax"]].isel({layer_dim: layer_slice})
    dksat_top = ds_soil[soil_vars["dksat"]].isel({layer_dim: 0})
    imperv = ds_soil[soil_vars["imperv"]]

    soil_params = SoilParams(
        smcmax=smcmax.values,
        DZS=DZS_ctl.values,
        dksat_top=dksat_top.values,
        imperv=imperv.values,
    )

    return dt, QRAIN, ACSNOM, UGDRNOFF, SM_ctl, SW_ctl, SNEQV, soil_params


# -----------------------------------------------------------------------------
# High-level driver
# -----------------------------------------------------------------------------

def compute_from_nwm(
    x_rng, y_rng,
    ldas_pattern: str, soil_param_path: str, wrfinput_path: str,
    depth: float = 0.4, use_smooth_freeze: bool = True
):
    ds_ldas, xy_slices = open_ldasout(ldas_pattern, x_rng=x_rng, y_rng=y_rng)
    ds_soil = open_soilparams(soil_param_path, xy_slices=xy_slices)
    ds_wrf = open_wrfinput(wrfinput_path, xy_slices=xy_slices)
    dt, QRAIN, ACS, UG, SM, SW, SNEQV, soil = prepare_control_depth_inputs(ds_ldas, ds_soil, ds_wrf, depth)

    # import matplotlib.pyplot as plt

    # sm2d = SM.isel(reference_time=1, time=1, soil_layers_stag=0)
    # preview = sm2d.coarsen(y=4, x=4, boundary="trim").mean()
    # plt.figure(figsize=(9,7))
    # # xarray handles dask lazy computation here
    # h = preview.plot.imshow(
    #     x="x", y="y", origin="lower", robust=True,
    #     cmap="viridis"  # or leave default
    # )
    # plt.gca().set_aspect("equal")
    # plt.title("Soil Moisture (time=0, layer=0)")
    # plt.xlabel("x (m)")
    # plt.ylabel("y (m)")
    # plt.show()

    time = QRAIN["time"]
    space_dims = [d for d in QRAIN.dims if d.lower() != "time"]
    shape = tuple(QRAIN.sizes[d] for d in space_dims)

    def _transpose_safe(da: xr.DataArray, desired: list[str]) -> xr.DataArray:
        # keep only dims that exist, preserve order, no duplicates
        ordered = []
        seen = set()
        for d in desired:
            if d in da.dims and d not in seen:
                ordered.append(d)
                seen.add(d)
        # append any remaining dims at the end
        tail = [d for d in da.dims if d not in seen]
        return da.transpose(*(ordered + tail))

    def _flat(da: xr.DataArray):
        desired = ["time", *space_dims]
        da_t = _transpose_safe(da, desired)
        T = da_t.sizes.get("time", 1)
        rest = int(da_t.size // max(T, 1))
        return da_t.values.reshape((T, rest))

    def _flat_layers(da: xr.DataArray):
        # detect layer dimension robustly
        layer_alias = {"soil_layers_stag", "soil_layers", "nsoil", "soil_layers_stag_stag"}
        dims = list(da.dims)
        ld = next((d for d in dims if d in layer_alias), None)
        if ld is None:
            ld = next(d for d in dims if d not in ("time", "reference_time", *space_dims))
        desired = [ld, "time", *space_dims]
        da_t = _transpose_safe(da, desired)
        L = da_t.sizes[ld]
        T = da_t.sizes.get("time", 1)
        rest = int(da_t.size // max(L * T, 1))
        return da_t.values.reshape((L, T, rest))

    QRAIN_f, ACS_f, UG_f = _flat(QRAIN), _flat(ACS), _flat(UG)
    SM_f, SW_f, SNEQV_f = _flat_layers(SM), _flat_layers(SW), _flat(SNEQV)

    soil_flat = SoilParams(
        smcmax=soil.smcmax.reshape((soil.smcmax.shape[0], -1)),
        DZS=soil.DZS.reshape((soil.DZS.shape[0], -1)),
        dksat_top=soil.dksat_top.reshape((-1,)),
        imperv=soil.imperv.reshape((-1,)),
    )

    q_surf, q_perc, diag = compute_fluxes(dt, QRAIN_f, ACS_f, UG_f, SM_f, SW_f, SNEQV_f, soil_flat, use_smooth_freeze)

    shp = (QRAIN.sizes["time"],) + shape
    q_surf = xr.DataArray(q_surf.reshape(shp), dims=("time", *space_dims),
                          coords={"time": time, **{d: QRAIN[d] for d in space_dims}}, name="q_surface")
    q_perc = xr.DataArray(q_perc.reshape(shp), dims=("time", *space_dims),
                          coords=q_surf.coords, name="q_perc")

    return q_surf, q_perc, diag


if __name__ == "__main__":
    from pylib import read
    from matplotlib import pyplot as plt
    schism_gd = read('/sciclone/schism10/feiye/STOFS3D-v8/I29f/hgrid.gr3')
    x_lower_left, y_lower_left = lonlat_to_ldas_xy(schism_gd.x.min(), schism_gd.y.min())
    x_upper_right, y_upper_right = lonlat_to_ldas_xy(schism_gd.x.max(), schism_gd.y.max())

    q_surf, q_perc, diag = compute_from_nwm(
        (x_lower_left, x_upper_right),
        (y_lower_left, y_upper_right),
        ldas_pattern="/sciclone/schism10/feiye/STOFS3D-v8/NWM/CONUS/netcdf/LDASOUT/2018/*.LDASOUT_DOMAIN1",
        soil_param_path=("/sciclone/schism10/feiye/STOFS3D-v8/NWM/"
                         "Parameters/v3.0_full_parameters/soilproperties_CONUS_FullRouting.nc"),
        wrfinput_path="/sciclone/schism10/feiye/STOFS3D-v8/NWM/Parameters/v3.0_full_parameters/wrfinput_CONUS.nc",
        depth=0.4,                # control depth in meters
        use_smooth_freeze=True,   # recommended
    )

    q_surf.isel(time=0).plot()
    q_perc.isel(time=0).plot()

    schism_gd.compute_ctr()
    schism_xctr_ldas, schism_yctr_ldas = lonlat_to_ldas_xy(schism_gd.xctr, schism_gd.yctr)
    # interpolate to SCHISM grid
    q_surf_schism = q_surf.interp(
        x=xr.DataArray(schism_xctr_ldas, dims="nface"),
        y=xr.DataArray(schism_yctr_ldas, dims="nface"),
        method="linear"
    )
    schism_gd.plot(fmt=1, value=q_surf_schism.isel(time=0).values)
    plt.show()
    print('Done.')
