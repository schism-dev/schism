from schism_py_pre_post.Geometry.inpoly import find_node_in_shpfiles
from pylib import schism_grid, proj
import numpy as np
import os


def set_additional_dp_v11_91(wdir='./', levee_info_dir='./Levee_info/Additional_Polygons/'):
    # v11.91, set upstream mississippi river levee height to 25 m, to prevent overtopping near the source
    #         set Bonnet Carre Spill Way gate to 8.5 m
    #         revert levee height near Ostrica, LA to DEM values

    gd_ll = schism_grid(f'{wdir}/hgrid.ll')
    gd_utm = schism_grid(f'{wdir}/hgrid.utm.gr3')
    gd_dem = schism_grid(f'{wdir}/hgrid.DEM_loaded.ll')

    for set_depth, shapefile in zip([-9, -20], [f'{levee_info_dir}/BonnetCarre_utm26918.shp', f'{levee_info_dir}/la_levee_center_line_upstream_missi_13m_buffer.shp']):
        i_inpoly = find_node_in_shpfiles(shapefile_names=[shapefile], gd=gd_utm)
        gd_ll.dp[i_inpoly] = np.minimum(gd_ll.dp[i_inpoly], set_depth)

    i_inpoly = find_node_in_shpfiles(shapefile_names=[f'{levee_info_dir}/la_levee_center_line_Ostrica_13m_buffer.shp'], gd=gd_utm)
    gd_ll.dp[i_inpoly] = gd_dem.dp[i_inpoly]

    gd_ll.save(f'{wdir}/hgrid.additional_dp.ll')

    os.system(f"cp {wdir}/hgrid.additional_dp.ll {wdir}/hgrid.ll")
    proj(
        f'{wdir}/hgrid.ll', 0, 'epsg:4326',
        f'{wdir}/hgrid.utm.gr3', 0, 'epsg:26918',
    )
    os.system(f"cp {wdir}/hgrid.utm.gr3 {wdir}/hgrid.utm.26918")

if __name__ == "__main__":
    # needs hgrid.ll and hgrid.utm.gr3 (epsg:26918) in wdir
    # needs Levee_info/, copy it into the wdir
    set_additional_dp_v11_91(
        wdir='./',
        levee_info_dir='./Levee_info/Additional_Polygons/'
    )