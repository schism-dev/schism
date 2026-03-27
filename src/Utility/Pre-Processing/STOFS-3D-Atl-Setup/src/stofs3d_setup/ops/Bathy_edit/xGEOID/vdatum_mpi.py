import numpy as np
import mpi4py.MPI as MPI

from pylib_essentials.schism_file import cread_schism_hgrid
from convert2xgeoid import point_conversion


comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


if rank == 0:
    hgrid_obj = cread_schism_hgrid('/sciclone/schism10/feiye/STOFS3D-v7/Runs/R00a/hgrid.gr3')
    dp_original = hgrid_obj.dp.copy()

    if hgrid_obj.np % size != 0:
        # append dummy data to the xyz, so that the number of points can be divided by size
        sub_size = hgrid_obj.np // size + 1
        n_dummy = sub_size * size - hgrid_obj.np
        # fill the dummy data with the last point
        x = np.r_[hgrid_obj.x, np.tile(hgrid_obj.x[-1], n_dummy)]
        y = np.r_[hgrid_obj.y, np.tile(hgrid_obj.y[-1], n_dummy)]
        z = np.r_[hgrid_obj.dp, np.tile(hgrid_obj.dp[-1], n_dummy)]
    else:
        sub_size = hgrid_obj.np // size
        x = hgrid_obj.x
        y = hgrid_obj.y
        z = hgrid_obj.dp
else:
    sub_size = None
    x = None; y = None; z = None
    

sub_size = comm.bcast(sub_size, root=0)

x_local = np.empty((sub_size,  ), dtype=float)
comm.Scatter(x, x_local, root=0)
x_local = x_local.flatten()

y_local = np.empty((sub_size,  ), dtype=float)
comm.Scatter(y, y_local, root=0)
y_local = y_local.flatten()

z_local = np.empty((sub_size,  ), dtype=float)
comm.Scatter(z, z_local, root=0)
z_local = z_local.flatten()

print('rank', rank, 'z_local.shape', z_local.shape, 'xyz_local:', np.c_[x_local, y_local, z_local][:3], "\n\n", np.c_[x_local, y_local, z_local][-3:])

z_local_updated = point_conversion(x_local, y_local, z_local, print_info=f'rank {rank}: ')

comm.Barrier()

# gather the updated z_local to the root rank
comm.Gather(z_local_updated, z, root=0)

if rank == 0:
    # remove the dummy data
    x = x[:-n_dummy]; y = y[:-n_dummy]; z = z[:-n_dummy]
    # find all changed points
    changed = np.any(np.abs(z - dp_original) > 1e-4)
    i_changed = np.where(changed)[0]
    np.savetxt('z_local_updated.txt', np.c_[i_changed+1, x[changed], y[changed], z[changed], dp_original[changed]])

pass


