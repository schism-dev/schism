'''
Sample functions for post-processing SCHISM output
'''

from copy import deepcopy
import numpy as np
import xarray as xr


def split_quads(elements=None):
    '''
    Split quad elements to triangles and
    append additional elements to the end of the element table

    elements is read from SCHISM output file as:
        elements = ds['SCHISM_hgrid_face_nodes'][:]
    '''

    if elements is None:
        raise ValueError('No elements provided')

    if elements.shape[1] == 3:  # already triangles
        return elements
    elif elements.shape[1] != 4:
        raise ValueError('elements should be a numpy array of (n,3) or (n,4)')

    triangles = deepcopy(elements)
    quad_idx = ~elements[:, -1].mask

    # split quads into two triangles
    quads = elements[quad_idx]
    upper_triangle = np.c_[  # last node is -1 (not applicable)
        quads[:, 0], quads[:, 1], quads[:, 3], -np.ones((quads.shape[0], 1))]
    lower_triangle = np.c_[  # last node is -1 (not applicable)
        quads[:, 1], quads[:, 2], quads[:, 3], -np.ones((quads.shape[0], 1))]

    # replace quads with upper triangle
    triangles[quad_idx, :] = upper_triangle
    # append lower triangle a the end
    triangles = np.ma.concatenate([triangles, lower_triangle], axis=0)
    # mask the last node, because all quads have been changed to triangles
    triangles.mask[:, -1] = True

    return triangles[:, :3]  # only return the first 3 nodes of each element


def cal_maxelev(depth: np.ndarray, elev: np.ndarray, fill_value: float = -99999.0):
    """
    Calculate the maximum elevation and the time of the maximum elevation

    - inputs:
        elev: np.ndarray of shape (ntimes, npoints)
        depth: np.ndarray of shape (npoints,)
        fill_value: the value to fill the masked values in elev
    - outputs:
        maxelev: np.ndarray of shape (npoints,)
        time_idx_maxelev: np.ndarray of shape (npoints,)
    """
    if fill_value > -1000:
        raise ValueError("fill_value should be a large negative number")

    # mask large values
    elev[np.where(elev > 100000)] = fill_value
    # mask dry nodes
    elev[elev + depth <= 1e-6] = fill_value  # native schism outputs
    elev[np.isnan(elev)] = fill_value  # deprecated, adcirc format

    maxelev = np.max(elev, axis=0)
    time_idx_maxelev = np.argmax(elev, axis=0)

    return maxelev, time_idx_maxelev


def cal_disturbance(
    depth: np.ndarray, elev: np.ndarray,
    city_node_idx_file: str = None,
    fillvalue: float = -99999.0
):
    """
    Calculate the maximum disturbance

    - inputs:
        depth: np.ndarray of shape (npoints,)
        elev: np.ndarray of shape (ntimes, npoints)
    - outputs:
        maxdist: np.ndarray of shape (npoints,)
    """
    if fillvalue > -1000:
        raise ValueError("fillvalue should be a large negative number")

    elev[np.isnan(elev)] = fillvalue  # deprecated, adcirc format

    disturb = elev.copy()  # same as elev in the ocean, so initialize with elev

    # read handle city node indices
    if city_node_idx_file is not None:
        city_node_idx = np.loadtxt(city_node_idx_file, encoding='utf-8').astype(bool)
    else:
        city_node_idx = np.zeros_like(depth, dtype=bool)
    # define land nodes, including city nodes
    land_node_idx = (depth < 0) | city_node_idx

    # Define disturbance:
    # On land, disturbance is the sum of elevation and depth, i.e., the water depth.
    # Also, disturbance is zero if the water depth is negative.
    disturb[:, land_node_idx] = np.maximum(0, elev[:, land_node_idx] + depth[land_node_idx])

    max_disturb = np.max(disturb, axis=0)
    time_idx_max_disturb = np.argmax(disturb, axis=0)

    # mask small max disturbance (< 0.3 m) on land (including cities)
    small_dist_on_land = (max_disturb < 0.3) * land_node_idx  # True if both conditions are met
    max_disturb[small_dist_on_land] = fillvalue

    return max_disturb, time_idx_max_disturb


def sample_usage():
    """
    Sample usage of cal_maxelev and cal_disturbance
    """

    # read in multiple SCHISM output files at a time,
    # useful for operational forecast where each day contains two files
    ds = xr.open_mfdataset([
        './outputs/schout_adcirc_20240926.nc',
        './outputs/schout_adcirc_20240927.nc',
    ], concat_dim='time', combine='nested')

    # read in the depth and elevation
    elev = ds['zeta'].values
    depth = ds['depth'].values
    # accomodate for a bug which set depth's dimension to (time, x, y)
    depth = depth[0, :] if depth.ndim == 2 else depth

    # calculate the maximum elevation and the time of the maximum elevation
    max_elev, time_idx_max_elev = cal_maxelev(depth, elev)
    # calculate the maximum disturbance and the time of the maximum disturbance
    max_disturb, time_idx_max_disturb = cal_disturbance(depth, elev, city_node_idx_file=None)

    # Extra detail: mask small disturbance in cities,
    # which need city node indices from a file (different between operation and shadow forecasts)
    # max_disturb, time_idx_max_disturb = cal_disturbance(
    #     depth, elev, city_node_idx_file='./inputs/city_poly.node_id.oper.txt')

    # print some diagnostics
    max_elev = np.ma.masked_values(max_elev, -99999.0)
    max_disturb = np.ma.masked_values(max_disturb, -99999.0)
    print(f'time_idx_max_elev: min={np.min(time_idx_max_elev)}, max={np.max(time_idx_max_elev)}')
    print(f'max_disturb: min={np.nanmin(max_disturb)}, max={np.nanmax(max_disturb)}')


if __name__ == '__main__':
    sample_usage()
