import numpy as np
from scipy.spatial import cKDTree
from pylib_experimental.schism_file import cread_schism_hgrid

hg_tranfer = cread_schism_hgrid(
    '/sciclone/schism10/Hgrid_projects/STOFS3D-v8/v28/v27_to_v28_splice.gr3')
hg_tranfer.proj(prj0='esri:102008', prj1='epsg:4326')
hg1 = cread_schism_hgrid(
    '/sciclone/schism10/feiye/STOFS3D-v8/I13y_v7/hgrid.gr3')
hg2 = cread_schism_hgrid(
    '/sciclone/schism10/feiye/STOFS3D-v8/I13y/Bathy_edit/RiverArc_Dredge/hgrid_dredged.gr3')

transfer_dist, transfer_idx1 = cKDTree(np.c_[hg1.x, hg1.y]).query(np.c_[hg_tranfer.x, hg_tranfer.y], k=1)
if np.any(transfer_dist > 1e-5):
    print(f'Warning: {sum(transfer_dist > 1e-5)} points are not matched between hg_tranfer and hg1')

transfer_dist, transfer_idx2 = cKDTree(np.c_[hg2.x, hg2.y]).query(np.c_[hg_tranfer.x, hg_tranfer.y], k=1)
if np.any(transfer_dist > 1e-5):
    print(f'Warning: {sum(transfer_dist > 1e-5)} points are not matched between hg_tranfer and hg1')

hg2.dp[transfer_idx2] = hg1.dp[transfer_idx1]
hg2.save(
    '/sciclone/schism10/feiye/STOFS3D-v8/I13y/Bathy_edit/hgrid_dredged_transferred_R27_in_splice_region.gr3', fmt=1)

print('Done')