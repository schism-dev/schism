"""
Check if any feeders are isolated from the main hgrid domain
"""


from schism_py_pre_post.Utilities.import_util import get_hgrid_reader

read_hgrid = get_hgrid_reader()


hg = read_hgrid('/sciclone/schism10/feiye/STOFS3D-v8/I15d4_v7/hgrid.gr3')
hg.compute_bnd()

n_true_land_boundary = hg.nob if hg.nob > 0 else 1

n_isolated = sum(hg.island == 0) - n_true_land_boundary
# which is true, because there should be same number of land and open boundaries
print(f'Number of isolated feeders: {n_isolated}')

with open('/sciclone/schism10/feiye/STOFS3D-v8/I15c_v7/isolated_feeders.txt', 'w') as f:
    f.write('lon lat\n')
    for i, is_island_boundary in enumerate(hg.island):
        if i >= n_true_land_boundary and is_island_boundary == 0:  # first nob land boundary is true
            nodes = hg.ilbn[i]
            for lon, lat in zip(hg.x[nodes], hg.y[nodes]):
                f.write(f'{lon} {lat} {i-n_true_land_boundary+1}\n')

print('Done')
