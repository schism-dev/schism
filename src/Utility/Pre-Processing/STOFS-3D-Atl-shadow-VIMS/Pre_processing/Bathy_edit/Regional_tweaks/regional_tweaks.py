"""
Apply temporary regional tweaks to hgrid depth,
only in regions where the channel representation is not accurate.
"""

import numpy as np
from pylib import inside_polygon, read_schism_bpfile
try:
    # Try to import from pylib_experimental to speed up the grid reading
    from pylib_experimental.schism_file import cread_schism_hgrid as read_hgrid
    print("Using 'cread_schism_hgrid' from 'pylib_experimental.schism_file'")
except ImportError:
    # If that fails, fall back to pylib
    from pylib import read as read_hgrid
    print("Using 'read' from 'pylib'")


REGIONAL_TWEAKS = {
    'min_5m_ll_noPR': 5,
    'SabinePass': 7,
    'BergenPoint': 5,
    'Washington_3': 15,
    'Elk_river': 2,
    'Hudson_river': 16,
    'James_river': 14,
    'NorthEast_river': 5,
    'Rappahannock_river': 6,
    'Susquehanna_river': 10,
    'York_river': 10,
    'Androscoggin_Kennebec_rivers': 3,
    'Merrimack_river': 3,
    'Patuxent_river': 5,
    'Penobscot_river': 5,
    'Saco_river': 3,
    'StCroix_river': 5,
    'Oyster_landing': 1,
    'st_lawrence1': 10,
    'st_lawrence2': 10,
    'st_lawrence3': 10,
}


def tweak_hgrid_depth(hgrid, regional_tweaks, regions_dir):
    '''
    Set the minimum depth in the regions specified in regional_tweaks.
    '''
    for region, depth in regional_tweaks.items():
        region_file = f'{regions_dir}/{region}.reg'
        bp = read_schism_bpfile(region_file, fmt=1)
        in_region = inside_polygon(np.c_[hgrid.x, hgrid.y], bp.x, bp.y).astype(bool)
        tweak_idx = in_region * (hgrid.dp < depth)
        hgrid.dp[tweak_idx] = depth
        print(f'Applied min depth {depth} in region {region}')
    return hgrid


def sample():
    '''
    Sample usage of the tweak_hgrid_depth function.
    '''
    regions_dir = '/sciclone/schism10/Hgrid_projects/DEMs/regions/'
    wdir = '/sciclone/schism10/feiye/STOFS3D-v7/Inputs/v7_test/Bathy_edit/DEM_loading/'

    hgrid = read_hgrid(f'{wdir}/hgrid.ll.dem_loaded.mpi.gr3')
    hgrid = tweak_hgrid_depth(hgrid, REGIONAL_TWEAKS, regions_dir)
    hgrid.write(f'{wdir}/hgrid.tweaked.gr3')


if __name__ == '__main__':
    sample()
    print('Done')
