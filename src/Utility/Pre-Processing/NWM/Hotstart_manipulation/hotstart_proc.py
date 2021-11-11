from pylib import schism_grid, schism_vgrid, WriteNC, inside_polygon  # from ZG's pylib: git@github.com:wzhengui/pylibs.git
from scipy import spatial
from scipy.interpolate import griddata
from scipy.sparse import csc_matrix
import numpy as np
from time import time
import os
import xarray as xr
# import pymp  # pip install pymp-pypi; setenv OMP_PROC_BIND true; srun --cpus-per-task=20 --time=01:00:00 python ./interp_hotstart4.py
import matplotlib.pyplot as plt
import gc
import shapefile


def nearest_neighbour(points_a, points_b):
    tree = spatial.cKDTree(points_b)
    return tree.query(points_a)[1]


class zdata():
    def __init__(self):
        pass  # dummy class to adhere to the style in pylib


class Hotstart():
    __slots__ = ('time', 'iths', 'ifile', 'idry_e', 'idry_s', 'idry', 'eta2', 'we', 'tr_el',
                 'tr_nd', 'tr_nd0', 'su2', 'sv2', 'q2', 'xl', 'dfv', 'dfh', 'dfq1', 'dfq2',
                 'nsteps_from_cold', 'cumsum_eta', 'file_format', 'dimname', 'vars', 'dims',
                 'vars_0d', 'vars_1d_node_based', 'vars_2d_node_based', 'vars_3d_node_based', 'var_wild')

    def __init__(self, ne, np, ns, nvrt, ntracers=2):
        self.file_format = 'NETCDF4'
        self.dimname = ['node', 'elem', 'side', 'nVert', 'ntracers', 'one']
        self.dims = [np, ne, ns, nvrt, ntracers, 1]
        self.vars = ['time', 'iths', 'ifile', 'idry_e', 'idry_s', 'idry', 'eta2', 'we', 'tr_el',
                     'tr_nd', 'tr_nd0', 'su2', 'sv2', 'q2', 'xl', 'dfv', 'dfh', 'dfq1', 'dfq2',
                     'nsteps_from_cold', 'cumsum_eta']
        self.vars_0d = ['time', 'iths', 'ifile', 'nsteps_from_cold']
        self.vars_1d_node_based = ['cumsum_eta']
        self.vars_2d_node_based = ['q2', 'xl', 'dfv', 'dfh', 'dfq1', 'dfq2']
        self.vars_3d_node_based = ['tr_nd', 'tr_nd0']
        self.var_wild = None
        pass

    def set_var(self, var_str="", val=None):
        vi = zdata()
        if var_str in ['time', 'iths', 'ifile', 'nsteps_from_cold']:
            vi.dimname = ('one',)
        elif var_str == 'idry_e':
            vi.dimname = ('elem',)
        elif var_str == 'idry_s':
            vi.dimname = ('side',)
        elif var_str in ['idry', 'eta2', 'cumsum_eta']:
            vi.dimname = ('node',)
        elif var_str == 'we':
            vi.dimname = ('elem', 'nVert')
        elif var_str in ['su2', 'sv2']:
            vi.dimname = ('side', 'nVert')
        elif var_str in ['q2', 'xl', 'dfv', 'dfh', 'dfq1', 'dfq2']:
            vi.dimname = ('node', 'nVert')
        elif var_str == 'tr_el':
            vi.dimname = ('elem', 'nVert', 'ntracers')
        elif var_str in ['tr_nd', 'tr_nd0']:
            vi.dimname = ('node', 'nVert', 'ntracers')
        else:
            raise Exception(f'{var_str} is not a variable in hotstart.nc')

        if var_str == 'time':
            vi.val = np.array(val+0.0)
        elif var_str in ['iths', 'ifile', 'nsteps_from_cold']:
            vi.val = np.array(val).astype('int')
        elif var_str in ['idry_e', 'idry_s', 'idry']:
            vi.val = val.astype('int32')
        else:
            vi.val = val

        exec(f"self.{var_str} = vi")

    def writer(self, fname):
        WriteNC(fname, self)


def gather_grid_info(grid_in, eta, dir):
    grid_in.compute_ctr()
    grid_in.isidenode = grid_in.compute_side(fmt=2)
    grid_in.side_x = (grid_in.x[grid_in.isidenode[:, 0]] + grid_in.x[grid_in.isidenode[:, 1]]) / 2.0
    grid_in.side_y = (grid_in.y[grid_in.isidenode[:, 0]] + grid_in.y[grid_in.isidenode[:, 1]]) / 2.0

    # temporary measure to reorder sides
    # tmp = np.loadtxt(f'{dir}/sidecenters.bp')
    # side_x = tmp[:, 1]
    # side_y = tmp[:, 2]
    # neighbor = nearest_neighbour(np.c_[side_x, side_y], np.c_[grid_in.side_x, grid_in.side_y])
    # grid_in.isidenode = grid_in.isidenode[neighbor, :]
    # grid_in.side_x = (grid_in.x[grid_in.isidenode[:, 0]] + grid_in.x[grid_in.isidenode[:, 1]]) / 2.0
    # grid_in.side_y = (grid_in.y[grid_in.isidenode[:, 0]] + grid_in.y[grid_in.isidenode[:, 1]]) / 2.0

    if not hasattr(grid_in, 'vgrid'):
        grid_in.vgrid = schism_vgrid()
        grid_in.vgrid.read_vgrid(f'{dir}/vgrid.in')

        grid_in.vgrid.zcor = grid_in.vgrid.compute_zcor(grid_in.dp, eta=eta)

        grid_in.vgrid.zcor_s = (grid_in.vgrid.zcor[grid_in.isidenode[:, 0]] + grid_in.vgrid.zcor[grid_in.isidenode[:, 1]]) / 2.0
        grid_in.vgrid.kbp_s = np.min(grid_in.vgrid.kbp[grid_in.isidenode[:, :]], axis=1)

        grid_in.vgrid.zcor_e = np.zeros((grid_in.ne, grid_in.vgrid.nvrt))
        grid_in.vgrid.kbp_e = np.zeros((grid_in.ne)).astype(int)
        for i in [3, 4]:
            II = (grid_in.i34 == i)
            for j in range(i):
                grid_in.vgrid.zcor_e[II, :] += grid_in.vgrid.zcor[grid_in.elnode[II, j], :] / i
            grid_in.vgrid.kbp_e[II] = np.min(grid_in.vgrid.kbp[grid_in.elnode[II, :i]], axis=1)

    return grid_in


def GetVerticalWeight(zcor_in, kbp_in, zcor_out, neighbors):
    # at grid_out's node, save vertical interp weights
    z_weight_lower = 0.0 * zcor_out
    z_idx_lower = (0.0 * zcor_out).astype(int)
    z_idx_upper = (0.0 * zcor_out).astype(int)
    zcor_tmp = zcor_in[neighbors]

    n_points = zcor_out.shape[0]
    nvrt = zcor_out.shape[1]
    dz = zcor_tmp[:, 1:] - zcor_tmp[:, :-1]

    # openmp-like parallelism, which can save some time
    # but the operations inside the loop are probably efficient enough
    # because they are already vectorized with respect to the vertical dimension
    # >>>
    # ncpus = int(os.getenv('SLURM_CPUS_PER_TASK', 1))
    # with pymp.Parallel(ncpus) as p:
    #    for i in p.range(n_points):
    for i in range(n_points):
            l_idx = np.searchsorted(zcor_tmp[i], zcor_out[i]) - 1
            below = l_idx == -1
            interior = (l_idx >= 0) & (l_idx < nvrt-1)
            above = (l_idx == (nvrt-1))

            z_weight_lower[i, interior] = (zcor_tmp[i, l_idx[interior]+1]-zcor_out[i, interior]) / dz[i, l_idx[interior]]
            z_idx_lower[i, interior] = l_idx[interior]
            z_idx_upper[i, interior] = l_idx[interior] + 1

            z_weight_lower[i, below] = 0.0
            z_idx_lower[i, below] = kbp_in[neighbors[i]]
            z_idx_upper[i, below] = kbp_in[neighbors[i]]

            z_weight_lower[i, above] = 1.0
            z_idx_lower[i, above] = l_idx[above]
            z_idx_upper[i, above] = l_idx[above]

    return [z_weight_lower, z_idx_lower, z_idx_upper]


def interp_hot(in_dir, out_dir, iplot=False, i_vert_interp=True):
    '''
    in_dir:
        must contain: hgrid.gr3, vgrid.in, and hotstart.nc
    out_dir:
        must contain: hgrid.gr3 and vgrid.in
    i_plot:
        True: save diagnostic plots of the interplated hotstart.nc
    i_vert_interp:
        True: linear interpolation in the vertical dimension
        False: assuming the vertical discretization is unchanged
    '''

    print('Using %d processors' % int(os.getenv('SLURM_CPUS_PER_TASK', 1)), flush=True)
    print('Using %d threads' % int(os.getenv('OMP_NUM_THREADS', 1)), flush=True)
    print('Using %d tasks' % int(os.getenv('SLURM_NTASKS', 1)), flush=True)

    h0 = 1e-5  # minimum water depth

    t = time()
    t0 = t

    # ----------------------- read base hotstart ----------------------------
    hot_in = xr.open_dataset(f'{in_dir}/hotstart.nc')
    eta2_in = np.array(hot_in.eta2)

    # ----------------------- read hgrids and vgrids ----------------------------
    grid_in = schism_grid(f'{in_dir}/hgrid.pkl')
    if not hasattr(grid_in, 'vgrid'):
        grid_in = gather_grid_info(grid_in, eta=eta2_in, dir=in_dir)
        grid_in.save(f'{in_dir}/hgrid.pkl')

    grid_out = schism_grid(f'{out_dir}/hgrid.pkl')
    eta2 = griddata(np.c_[grid_in.x, grid_in.y], hot_in.eta2, (grid_out.x, grid_out.y), method='linear')
    eta2_tmp = griddata(np.c_[grid_in.x, grid_in.y], hot_in.eta2, (grid_out.x, grid_out.y), method='nearest')
    eta2[np.isnan(eta2)] = eta2_tmp[np.isnan(eta2)]
    if not hasattr(grid_out, 'vgrid'):
        grid_out = gather_grid_info(grid_out, eta=eta2, dir=out_dir)
        grid_out.save(f'{out_dir}/hgrid.pkl')

    print(f'Reading inputs and calculating geometry took {time()-t} seconds', flush=True)
    t = time()

    # ---------------find the nearest neighbors in the horizontal dimension ----------------------------
    neighbor = nearest_neighbour(np.c_[grid_out.x, grid_out.y], np.c_[grid_in.x, grid_in.y])
    neighbor_e = nearest_neighbour(np.c_[grid_out.xctr, grid_out.yctr], np.c_[grid_in.xctr, grid_in.yctr])
    neighbor_s = nearest_neighbour(np.c_[grid_out.side_x, grid_out.side_y], np.c_[grid_in.side_x, grid_in.side_y])
    print(f'Finding nearest neighbors for nodes/elements/sides took {time()-t} seconds', flush=True)
    t = time()

    if iplot:
        plt.scatter(grid_in.side_x, grid_in.side_y, s=20, c=hot_in.su2.data[:, -1], cmap='jet', vmin=0, vmax=2)
        plt.savefig(f'{out_dir}/su2_in2.png', dpi=700)
        plt.scatter(grid_in.x, grid_in.y, s=20, c=hot_in.tr_nd.data[:, -1, 0], cmap='jet', vmin=0, vmax=33)
        plt.savefig(f'{out_dir}/trnd_in2.png', dpi=700)
        plt.scatter(grid_in.xctr, grid_in.yctr, s=20, c=hot_in.tr_el.data[:, -1, 0], cmap='jet', vmin=0, vmax=33)
        plt.savefig(f'{out_dir}/trel_in2.png', dpi=700)

    # ----------------------- calculate vertical interpolation weights ----------------------------
    if i_vert_interp:
        if not hasattr(grid_out.vgrid, "z_weight_lower"):
            [grid_out.vgrid.z_weight_lower, grid_out.vgrid.z_idx_lower, grid_out.vgrid.z_idx_upper] = GetVerticalWeight(
                grid_in.vgrid.zcor, grid_in.vgrid.kbp, grid_out.vgrid.zcor, neighbors=neighbor)
        print(f'Calculating node-based vertical interpolation weights took {time()-t} seconds', flush=True)
        t = time()
        if not hasattr(grid_out.vgrid, "ze_weight_lower"):
            [grid_out.vgrid.ze_weight_lower, grid_out.vgrid.ze_idx_lower, grid_out.vgrid.ze_idx_upper] = GetVerticalWeight(
                grid_in.vgrid.zcor_e, grid_in.vgrid.kbp_e, grid_out.vgrid.zcor_e, neighbors=neighbor_e)
        print(f'Calculating element-based vertical interpolation weights took {time()-t} seconds', flush=True)
        t = time()
        if not hasattr(grid_out.vgrid, "zs_weight_lower"):
            [grid_out.vgrid.zs_weight_lower, grid_out.vgrid.zs_idx_lower, grid_out.vgrid.zs_idx_upper] = GetVerticalWeight(
                grid_in.vgrid.zcor_s, grid_in.vgrid.kbp_s, grid_out.vgrid.zcor_s, neighbors=neighbor_s)
        print(f'Calculating side-based vertical interpolation weights took {time()-t} seconds', flush=True)
        t = time()
        grid_out.save(f'{out_dir}/hgrid.pkl')

    # set idry, idry_e and idry_s, based on eta2 rather than from interpolation
    idry = (eta2 < -grid_out.dp+h0)
    # An element is wet if and only if depths at all nodes >h0
    idry_e = np.ones(grid_out.ne)
    for i34 in [3, 4]:
        idx = (grid_out.i34 == i34)
        idry_e[idx] = np.amax(idry[grid_out.elnode[idx, 0: i34]], axis=1).astype(int)
    # np.savetxt(f'{out_dir}/idry_e.prop', np.c_[list(range(1, grid_out.ne+1)), idry_e], fmt='%i')

    # slightly different from SCHISM:
    # SCHISM: A node is wet if and only if at least one surrounding element is wet. This script: skipped
    # SCHISM: A side is wet if and only if at least one surrounding element is wet. This script: changed to both nodes are wet
    grid_out.isidenode = grid_out.compute_side(fmt=2)
    idry_s = np.amax(idry[grid_out.isidenode], axis=1).astype(int)

    print(f'Setting dry indicators took {time()-t} seconds', flush=True)
    t = time()

    # ------------ trivial variables --------------------------------
    ntracers = hot_in.tr_nd.data.shape[-1]
    nvrt_out = grid_out.vgrid.nvrt
    hot_out = Hotstart(grid_out.ne, grid_out.np, grid_out.ns, nvrt_out, ntracers)
    for var_str in hot_out.vars:
        if var_str in ['eta2', 'idry', 'idry_e', 'idry_s']:
            exec(f"hot_out.var_wild = {var_str}")  # already set
        elif var_str in hot_out.vars_0d:  # single number
            exec(f'hot_out.var_wild = hot_in.{var_str}.data')
        elif var_str in hot_out.vars_1d_node_based:
            exec(f'hot_out.var_wild = hot_in.{var_str}.data[neighbor]')
        elif var_str in hot_out.vars_2d_node_based + hot_out.vars_3d_node_based:
            if i_vert_interp:
                continue  # to be dealt with later
            else:
                exec(f'hot_out.var_wild = hot_in.{var_str}.data[neighbor]')
        elif var_str in ['we', 'tr_el']:
            if i_vert_interp:
                continue  # to be dealt with later
            else:
                exec(f'hot_out.var_wild = hot_in.{var_str}.data[neighbor_e]')
        elif var_str in ['su2', 'sv2']:
            if i_vert_interp:
                continue  # to be dealt with later
            else:
                exec(f'hot_out.var_wild = hot_in.{var_str}.data[neighbor_s]')
        else:
            raise Exception(f'{var_str} not in list')
        hot_out.set_var(var_str=var_str, val=hot_out.var_wild)
        print(f'Processing {var_str} took {time()-t} seconds', flush=True)
        t = time()

    if i_vert_interp:
        # ----------------------- we and tr_el ----------------------------
        we_out = np.zeros((hot_out.dims[1], hot_out.dims[3]))
        we_tmp = hot_in.we.data[neighbor_e]
        trel_out = np.zeros((hot_out.dims[1], hot_out.dims[3], hot_out.dims[4]))
        trel_tmp = hot_in.tr_el.data[neighbor_e]
        print(f'Reading we and tr_el took {time()-t} seconds', flush=True)
        t = time()
        row = np.r_[np.array(range(hot_out.dims[1])), np.array(range(hot_out.dims[1]))]
        for k in range(hot_out.dims[3]):
            col = np.r_[grid_out.vgrid.ze_idx_lower[:, k], grid_out.vgrid.ze_idx_upper[:, k]]
            data = np.r_[grid_out.vgrid.ze_weight_lower[:, k], 1.0 - grid_out.vgrid.ze_weight_lower[:, k]]
            weights = csc_matrix((data, (row, col)), shape=(hot_out.dims[1], hot_out.dims[3])).toarray()
            we_out[:, k] = np.sum(we_tmp * weights, axis=1)
            for j in range(hot_out.dims[4]):
                trel_out[:, k, j] = np.sum(trel_tmp[:, :, j] * weights, axis=1)
            print(f'Processing Layer {k+1} of {hot_out.dims[3]} for we and tr_el took {time()-t} seconds', flush=True)
            t = time()
        hot_out.set_var(var_str='we', val=we_out)
        hot_out.set_var(var_str='tr_el', val=trel_out)
        if iplot:
            plt.scatter(grid_out.xctr, grid_out.yctr, s=20, c=trel_out[:, -1, 0], cmap='jet', vmin=0, vmax=33)
            # plt.show()
            plt.savefig(f'{out_dir}/trel_out2.png', dpi=700)

        for delete_var in ['we_out', 'we_tmp', 'trel_out', 'trel_tmp']:
            exec(f'del {delete_var}')
        gc.collect()

        # ----------------------- su2, sv2 ----------------------------
        row = np.r_[np.array(range(hot_out.dims[2])), np.array(range(hot_out.dims[2]))]
        for var_str in ['su2', 'sv2']:
            exec(f'{var_str}_out = np.zeros((hot_out.dims[2], hot_out.dims[3]))')
            exec(f'{var_str}_tmp = hot_in.{var_str}.data[neighbor_s]')
        print(f'Reading su2 and sv2 took {time()-t} seconds', flush=True)
        t = time()
        for k in range(hot_out.dims[3]):
            col = np.r_[grid_out.vgrid.zs_idx_lower[:, k], grid_out.vgrid.zs_idx_upper[:, k]]
            data = np.r_[grid_out.vgrid.zs_weight_lower[:, k], 1.0 - grid_out.vgrid.zs_weight_lower[:, k]]
            weights = csc_matrix((data, (row, col)), shape=(hot_out.dims[2], hot_out.dims[3])).toarray()
            su2_out[:, k] = np.sum(su2_tmp * weights, axis=1)
            sv2_out[:, k] = np.sum(sv2_tmp * weights, axis=1)
            print(f'Processing Layer {k+1} of {hot_out.dims[3]} for su2 and sv2 took {time()-t} seconds', flush=True)
            t = time()
        hot_out.set_var(var_str='su2', val=su2_out)
        hot_out.set_var(var_str='sv2', val=sv2_out)
        if iplot:
            plt.scatter(grid_out.side_x, grid_out.side_y, s=20, c=su2_out[:, -1], cmap='jet', vmin=0, vmax=2)
            plt.savefig(f'{out_dir}/su2_out2.png', dpi=700)
        for delete_var in ['su2_out', 'su2_tmp', 'sv2_out', 'sv2_tmp']:
            exec(f'del {delete_var}')
        gc.collect()

        # ----------------------- vars(np, nvrt) and vars(np, nvrt, ntracers) ----------------------------
        for var_str in hot_out.vars_2d_node_based:
            exec(f'{var_str}_out = np.zeros((hot_out.dims[0], hot_out.dims[3]))')
            exec(f'{var_str}_tmp = hot_in.{var_str}.data[neighbor]')
            print(f'Reading {var_str} took {time()-t} seconds', flush=True)
            t = time()
        trnd_out = np.zeros((hot_out.dims[0], hot_out.dims[3], hot_out.dims[4]))
        trnd_tmp = hot_in.tr_nd.data[neighbor]
        print(f'Reading tr_nd took {time()-t} seconds', flush=True)
        t = time()

        row = np.r_[np.array(range(hot_out.dims[0])), np.array(range(hot_out.dims[0]))]
        for k in range(hot_out.dims[3]):
            col = np.r_[grid_out.vgrid.z_idx_lower[:, k], grid_out.vgrid.z_idx_upper[:, k]]
            data = np.r_[grid_out.vgrid.z_weight_lower[:, k], 1.0 - grid_out.vgrid.z_weight_lower[:, k]]
            weights = csc_matrix((data, (row, col)), shape=(hot_out.dims[0], hot_out.dims[3])).toarray()
            for var_str in hot_out.vars_2d_node_based:
                exec(f'{var_str}_out[:, k] = np.sum({var_str}_tmp * weights, axis=1)')

            for j in range(hot_out.dims[4]):  # loop ntracers
                trnd_out[:, k, j] = np.sum(trnd_tmp[:, :, j] * weights, axis=1)
                pass

            print(f'Processing Layer {k+1} of {hot_out.dims[3]} for all N-dimensional node-based variables took {time()-t} seconds', flush=True)
            t = time()

        for var_str in hot_out.vars_2d_node_based:
            exec(f'hot_out.set_var(var_str=var_str, val={var_str}_out)')
            exec(f'del {var_str}_out')
            exec(f'del {var_str}_tmp')
        gc.collect()

        hot_out.set_var(var_str='tr_nd', val=trnd_out)
        hot_out.set_var(var_str='tr_nd0', val=trnd_out)
        if iplot:
            plt.scatter(grid_out.x, grid_out.y, s=20, c=trnd_out[:, -1, 0], cmap='jet', vmin=0, vmax=33)
            plt.savefig(f'{out_dir}/trnd_out2.png', dpi=700)
        del trnd_out
        del trnd_tmp
        gc.collect()

    # ----------------------- write new hotstart.nc ----------------------------
    hot_out.writer(f'{out_dir}/new_hotstart.nc')
    print(f'Total time: {time()-t0} seconds', flush=True)


def find_ele_node_in_shpfile(shapefile_name, grid):
    '''
    Find element/node index within polygons of a list of shapefiles
    '''
    grid.compute_ctr()

    sf = shapefile.Reader(shapefile_name)
    shapes = sf.shapes()

    ele_ind_list = []
    for shp in shapes:
        poly_xy = np.array(shp.points).T
        ind = inside_polygon(np.c_[grid.xctr, grid.yctr], poly_xy[0], poly_xy[1])  # 1: true; 0: false
        ind = ind.astype('bool')
        ele_ind_list.append(ind)

    node_ind_list = []
    for shp in shapes:
        poly_xy = np.array(shp.points).T
        ind = inside_polygon(np.c_[grid.x, grid.y], poly_xy[0], poly_xy[1])  # 1: true; 0: false
        ind = ind.astype('bool')
        node_ind_list.append(ind)

    return [ele_ind_list, node_ind_list]


def replace_hot_vars(infile, outfile, grid, vars_list=[], shapefile_name=None, n_points=0):
    hot_in = xr.open_dataset(infile)
    hot_out = xr.open_dataset(outfile)

    for dim_name in ['elem', 'side', 'node', 'nVert', 'ntracers']:
        if hot_in.dims[dim_name] != hot_out.dims[dim_name]:
            raise Exception(f'{dim_name} not equal')
    if grid.ne != hot_in.dims['elem']:
        raise Exception('grid not consistent with hotstart.nc')

    if shapefile_name is not None:
        [ele_idx_list, node_idx_list] = find_ele_node_in_shpfile(
            # order matters, latter prevails
            shapefile_name=shapefile_name,
            grid=grid,
        )

    for var in vars_list:
        exec(f'grid.n_points = hot_out.{var}.data.shape[0]')
        if grid.n_points == hot_in.dims['elem']:
            for ind in ele_idx_list:
                exec(f'hot_out.{var}.data[ind] = hot_in.{var}.data[ind]')
        elif grid.n_points == hot_in.dims['node']:
            for ind in node_idx_list:
                exec(f'hot_out.{var}.data[ind] = hot_in.{var}.data[ind]')
        else:
            raise Exception(f'unknown dimension {grid.n_points}')

    hot_out.to_netcdf(outfile + '.new')


if __name__ == "__main__":
    # Sample usage
    replace_hot_vars(
        infile='/sciclone/schism10/feiye/From_Nabi/RUN02/Test_Hot/hotstart_it=31104.nc',
        outfile='/sciclone/schism10/feiye/From_Nabi/RUN02/Test_Hot/hotstart.nc.0',
        grid=schism_grid('/sciclone/schism10/feiye/From_Nabi/RUN02/hgrid.ll'),
        vars_list=['tr_nd', 'tr_nd0', 'tr_el'],
        shapefile_name='/sciclone/schism10/feiye/From_Nabi/RUN02/Test_Hot/ocean.shp'
    )

    interp_hot(
        in_dir='/sciclone/schism10/feiye/ICOGS/Ida01/',
        out_dir='/sciclone/schism10/feiye/ICOGS/RUN20f/',
        iplot=False,
        i_vert_interp=False
    )