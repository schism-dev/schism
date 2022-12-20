from pylib import schism_grid, schism_vgrid, WriteNC, inside_polygon  # from ZG's pylib: pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple pylibs4schism==0.1.10
from scipy import spatial
from scipy.interpolate import griddata
from scipy.sparse import csc_matrix
import numpy as np
from time import time
import os
import xarray as xr
# import pymp  # (The speed gain is probably not worth the trouble) pip install pymp-pypi; setenv OMP_PROC_BIND true; srun --cpus-per-task=20 --time=01:00:00 python ./interp_hotstart4.py
import matplotlib.pyplot as plt
import shapefile
import copy


def nearest_neighbour(points_a, points_b):
    tree = spatial.cKDTree(points_b)
    return tree.query(points_a)[1]


class zdata():
    def __init__(self):
        pass  # dummy class to adhere to pylib's style


class Hotstart():
    '''
    A class for SCHISM's hotstart.nc.

    Different ways to instantiate a hotstart instance:

        (1) Hotstart(grid_info=some_dir), where some_dir contains hgrid.gr3 and vgrid.in; the values of the hotstart variables are initialized to 0 mostly.

        (2) Hotstart(grid_info=some_dir, hot_file=some_file), same as (1), but with the hotstart variable values from some_file.

        (3) Hotstart(grid_info={'ne': 2, 'np': 4, 'ns': 5, 'nvrt': 4}); the values of the hotstart variables are initialized to 0 mostly.
    '''

    __slots__ = ('time', 'iths', 'ifile', 'idry_e', 'idry_s', 'idry', 'eta2', 'we', 'tr_el',
                 'tr_nd', 'tr_nd0', 'su2', 'sv2', 'q2', 'xl', 'dfv', 'dfh', 'dfq1', 'dfq2',
                 'nsteps_from_cold', 'cumsum_eta', 'file_format', 'dimname', 'vars', 'dims',
                 'vars_0d', 'vars_1d_node_based', 'vars_2d_node_based', 'vars_3d_node_based', 'var_wild',
                 'grid', 'source_dir', 'var_dict')

    def __init__(self, grid_info, ntracers=2, hot_file=None):
        self.file_format = 'NETCDF4'
        self.dimname = ['node', 'elem', 'side', 'nVert', 'ntracers', 'one']
        if isinstance(grid_info, dict):
            self.dims = [grid_info['np'], grid_info['ne'], grid_info['ns'], grid_info['nvrt'], ntracers, 1]
        elif isinstance(grid_info, str):
            self.source_dir = grid_info
            self.grid = zdata()
            self.grid.hgrid = schism_grid(f'{self.source_dir}/hgrid.gr3')
            self.grid.hgrid.compute_side(fmt=1)
            self.grid.vgrid = schism_vgrid()
            self.grid.vgrid.read_vgrid(f'{self.source_dir}/vgrid.in')
            self.dims = [self.grid.hgrid.np, self.grid.hgrid.ne, self.grid.hgrid.ns, self.grid.vgrid.nvrt] + [ntracers, 1]
            if hot_file is not None:
                my_hot = xr.open_dataset(hot_file)
                if my_hot.dims['elem'] != self.dims[1]:
                    raise Exception(f'Inconsistent geometry: {self.source_dir}')
                if ntracers != my_hot.dims['ntracers']:
                    ntracers = my_hot.dims['ntracers']
                    print(f'Warning: inconsistent ntracers, setting ntracers={ntracers} based on the value in the hotstart file')

        self.vars = ['time', 'iths', 'ifile', 'idry_e', 'idry_s', 'idry', 'eta2', 'we', 'tr_el',
                     'tr_nd', 'tr_nd0', 'su2', 'sv2', 'q2', 'xl', 'dfv', 'dfh', 'dfq1', 'dfq2',
                     'nsteps_from_cold', 'cumsum_eta']

        self.vars_0d = ['time', 'iths', 'ifile', 'nsteps_from_cold']
        self.vars_1d_node_based = ['cumsum_eta']
        self.vars_2d_node_based = ['q2', 'xl', 'dfv', 'dfh', 'dfq1', 'dfq2']
        self.vars_3d_node_based = ['tr_nd', 'tr_nd0']
        self.var_wild = None

        if hot_file is None:  # initialize with 0
            self.set_var('time', 0.0)
            self.set_var('iths', 0)
            self.set_var('ifile', 1)

            self.set_var('idry_e', np.zeros(self.dims[1]))  # all wet
            self.set_var('idry_s', np.zeros(self.dims[2]))  # all wet
            self.set_var('idry', np.zeros(self.dims[0]))  # all wet

            self.set_var('eta2', np.zeros(self.dims[0]))
            self.set_var('we', np.zeros((self.dims[1], self.dims[3])))
            self.set_var('tr_el', np.zeros((self.dims[1], self.dims[3], ntracers)))
            self.set_var('tr_nd', np.zeros((self.dims[0], self.dims[3], ntracers)))
            self.set_var('tr_nd0', np.zeros((self.dims[0], self.dims[3], ntracers)))
            self.set_var('su2', np.zeros((self.dims[2], self.dims[3])))
            self.set_var('sv2', np.zeros((self.dims[2], self.dims[3])))

            self.set_var('q2', np.zeros((self.dims[0], self.dims[3])))
            self.set_var('xl', np.zeros((self.dims[0], self.dims[3])))
            self.set_var('dfv', np.zeros((self.dims[0], self.dims[3])))
            self.set_var('dfh', np.zeros((self.dims[0], self.dims[3])))
            self.set_var('dfq1', np.zeros((self.dims[0], self.dims[3])))
            self.set_var('dfq2', np.zeros((self.dims[0], self.dims[3])))

            self.set_var('nsteps_from_cold', 0)
            self.set_var('cumsum_eta', np.zeros(self.dims[0]))
        else:  # read from existing hotstart.nc
            for var_str in self.vars:
                try:
                    self.set_var(var_str, my_hot[var_str].data)
                except KeyError:
                    if var_str in ['nsteps_from_cold', 'cumsum_eta']:
                        self.set_var(var_str, 0)

        self.var_dict = {}
        for var_str in self.vars:
            exec(f'self.var_dict[var_str] = self.{var_str}')

    def set_var(self, var_str="", val=None):
        vi = zdata()
        if var_str in self.vars_0d:
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
        elif var_str in self.vars_2d_node_based:
            vi.dimname = ('node', 'nVert')
        elif var_str == 'tr_el':
            vi.dimname = ('elem', 'nVert', 'ntracers')
        elif var_str in self.vars_3d_node_based:
            vi.dimname = ('node', 'nVert', 'ntracers')
        else:
            raise Exception(f'{var_str} is not a variable in hotstart.nc')

        if var_str == 'time':
            vi.val = np.array(val).astype('float64')
        elif var_str in ['iths', 'ifile', 'nsteps_from_cold']:
            vi.val = np.array(val).astype('int32')
        elif var_str in ['idry_e', 'idry_s', 'idry']:
            vi.val = np.array(val).astype('int32')
        else:
            vi.val = np.array(val).astype('float64')

        exec(f"self.{var_str} = vi")

    def interp_from_existing_hotstart(self, hot_in, iplot=False, i_vert_interp=True):
        '''
        The interpolation method in this routine sacrifices some accuracy for efficiency,
        by using the nearest neighbor method in the horizontal dimension.
        hot_in:
            Background hotstart.nc as an instance of the <Hostart> class
        i_plot:
            True: save diagnostic plots of the interplated hotstart.nc; can take a few minutes for larger grid
            False (default): don't plot
        i_vert_interp:
            True: linear interpolation in the vertical dimension
            False: assuming the vertical discretization is unchanged (faster)
        '''

        print('Using %d processors' % int(os.getenv('SLURM_CPUS_PER_TASK', 1)), flush=True)
        print('Using %d threads' % int(os.getenv('OMP_NUM_THREADS', 1)), flush=True)
        print('Using %d tasks' % int(os.getenv('SLURM_NTASKS', 1)), flush=True)

        h0 = 1e-5  # minimum water depth
        plot_layer = -1  # surface
        t = time()
        t0 = t

        in_dir = hot_in.source_dir
        out_dir = self.source_dir

        # ----------------------- read base hotstart ----------------------------
        eta2_in = np.array(hot_in.eta2.val)

        # ----------------------- read hgrids and vgrids ----------------------------
        hot_in.grid = gather_grid_info(hot_in.grid, eta=eta2_in, dir=in_dir)
        grid_in = hot_in.grid

        grid_out = self.grid
        eta2 = griddata(np.c_[grid_in.hgrid.x, grid_in.hgrid.y], hot_in.eta2.val, (grid_out.hgrid.x, grid_out.hgrid.y), method='linear')
        eta2_tmp = griddata(np.c_[grid_in.hgrid.x, grid_in.hgrid.y], hot_in.eta2.val, (self.grid.hgrid.x, self.grid.hgrid.y), method='nearest')
        eta2[np.isnan(eta2)] = eta2_tmp[np.isnan(eta2)]
        grid_out = gather_grid_info(grid_out, eta=eta2, dir=out_dir)

        print(f'Reading inputs and calculating geometry took {time()-t} seconds', flush=True)
        t = time()

        # ---------------find the nearest neighbors in the horizontal dimension ----------------------------
        neighbor = nearest_neighbour(np.c_[grid_out.hgrid.x, grid_out.hgrid.y], np.c_[grid_in.hgrid.x, grid_in.hgrid.y])
        neighbor_e = nearest_neighbour(np.c_[grid_out.hgrid.xctr, grid_out.hgrid.yctr], np.c_[grid_in.hgrid.xctr, grid_in.hgrid.yctr])
        neighbor_s = nearest_neighbour(np.c_[grid_out.hgrid.side_x, grid_out.hgrid.side_y], np.c_[grid_in.hgrid.side_x, grid_in.hgrid.side_y])
        print(f'Finding nearest neighbors for nodes/elements/sides took {time()-t} seconds', flush=True)
        t = time()

        if iplot:
            plt.scatter(grid_in.hgrid.side_x, grid_in.hgrid.side_y, s=5, c=hot_in.su2.val[:, plot_layer], cmap='jet', vmin=0, vmax=2)
            plt.savefig(f'{out_dir}/su2_in2.png', dpi=700)
            plt.scatter(grid_in.hgrid.x, grid_in.hgrid.y, s=5, c=hot_in.tr_nd.val[:, plot_layer, 0], cmap='jet', vmin=0, vmax=33)
            plt.savefig(f'{out_dir}/trnd_in2.png', dpi=700)
            plt.scatter(grid_in.hgrid.xctr, grid_in.hgrid.yctr, s=5, c=hot_in.tr_el.val[:, plot_layer, 0], cmap='jet', vmin=0, vmax=33)
            plt.savefig(f'{out_dir}/trel_in2.png', dpi=700)
            print(f'Generating diagnostic outputs for hot_in took {time()-t} seconds', flush=True)
            t = time()

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

        # set idry, idry_e and idry_s, based on eta2 rather than from interpolation
        self.idry.val = (eta2 < -grid_out.hgrid.dp + h0).astype('int32')
        # An element is wet if and only if depths at all nodes >h0
        self.idry_e.val = np.ones(grid_out.hgrid.ne).astype('int32')
        for i34 in [3, 4]:
            idx = (grid_out.hgrid.i34 == i34)
            self.idry_e.val[idx] = np.amax(self.idry.val[grid_out.hgrid.elnode[idx, 0: i34]], axis=1).astype(int)
        # np.savetxt(f'{out_dir}/idry_e.prop', np.c_[list(range(1, grid_out.hgrid.ne+1)), idry_e], fmt='%i')

        # slightly different from SCHISM:
        # SCHISM: A node is wet if and only if at least one surrounding element is wet. This script: skipped
        # SCHISM: A side is wet if and only if at least one surrounding element is wet. This script: changed to both nodes are wet
        self.idry_s.val = np.amax(self.idry.val[grid_out.hgrid.isidenode], axis=1).astype('int32')
        self.eta2.val = eta2

        print(f'Setting dry indicators took {time()-t} seconds', flush=True)
        t = time()

        # ------------ trivial variables --------------------------------
        for var_str in self.vars:
            if var_str in ['eta2', 'idry', 'idry_e', 'idry_s']:
                continue  # already set
            elif var_str in self.vars_0d:  # single number
                self.var_dict[var_str].val = hot_in.var_dict[var_str].val
            elif var_str in self.vars_1d_node_based:
                self.var_dict[var_str].val = hot_in.var_dict[var_str].val[neighbor]
            elif var_str in self.vars_2d_node_based + self.vars_3d_node_based:
                if i_vert_interp:
                    continue  # to be dealt with later
                else:
                    self.var_dict[var_str].val = hot_in.var_dict[var_str].val[neighbor]
            elif var_str in ['we', 'tr_el']:
                if i_vert_interp:
                    continue  # to be dealt with later
                else:
                    self.var_dict[var_str].val = hot_in.var_dict[var_str].val[neighbor_e]
            elif var_str in ['su2', 'sv2']:
                if i_vert_interp:
                    continue  # to be dealt with later
                else:
                    self.var_dict[var_str].val = hot_in.var_dict[var_str].val[neighbor_s]
            else:
                raise Exception(f'{var_str} not in list')
            print(f'Processing {var_str} took {time()-t} seconds', flush=True)
            t = time()

        if i_vert_interp:
            # ----------------------- we and tr_el ----------------------------
            we_tmp = hot_in.we.val[neighbor_e]
            trel_tmp = hot_in.tr_el.val[neighbor_e]
            print(f'Reading we and tr_el took {time()-t} seconds', flush=True)
            t = time()
            row = np.r_[np.array(range(self.dims[1])), np.array(range(self.dims[1]))]
            for k in range(self.dims[3]):
                col = np.r_[grid_out.vgrid.ze_idx_lower[:, k], grid_out.vgrid.ze_idx_upper[:, k]]
                data = np.r_[grid_out.vgrid.ze_weight_lower[:, k], 1.0 - grid_out.vgrid.ze_weight_lower[:, k]]
                weights = csc_matrix((data, (row, col)), shape=(self.dims[1], hot_in.dims[3])).toarray()
                self.we.val[:, k] = np.sum(we_tmp * weights, axis=1)
                for j in range(self.dims[4]):
                    self.tr_el.val[:, k, j] = np.sum(trel_tmp[:, :, j] * weights, axis=1)
                print(f'Processing Layer {k+1} of {self.dims[3]} for we and tr_el took {time()-t} seconds', flush=True)
                t = time()
            if iplot:
                plt.scatter(grid_out.hgrid.xctr, grid_out.hgrid.yctr, s=5, c=self.tr_el.val[:, plot_layer, 0], cmap='jet', vmin=0, vmax=33)
                plt.savefig(f'{out_dir}/trel_out2.png', dpi=700)
                print(f'Generating diagnostic outputs for tr_el took {time()-t} seconds', flush=True)
                t = time()

            # ----------------------- su2, sv2 ----------------------------
            row = np.r_[np.array(range(self.dims[2])), np.array(range(self.dims[2]))]
            for k in range(self.dims[3]):
                col = np.r_[grid_out.vgrid.zs_idx_lower[:, k], grid_out.vgrid.zs_idx_upper[:, k]]
                data = np.r_[grid_out.vgrid.zs_weight_lower[:, k], 1.0 - grid_out.vgrid.zs_weight_lower[:, k]]
                weights = csc_matrix((data, (row, col)), shape=(self.dims[2], hot_in.dims[3])).toarray()
                self.su2.val[:, k] = np.sum(hot_in.var_dict['su2'].val[neighbor_s] * weights, axis=1)
                self.sv2.val[:, k] = np.sum(hot_in.var_dict['sv2'].val[neighbor_s] * weights, axis=1)
                print(f'Processing Layer {k+1} of {self.dims[3]} for su2 and sv2 took {time()-t} seconds', flush=True)
                t = time()
            if iplot:
                plt.scatter(grid_out.hgrid.side_x, grid_out.hgrid.side_y, s=5, c=self.su2.val[:, plot_layer], cmap='jet', vmin=0, vmax=2)
                plt.savefig(f'{out_dir}/su2_out2.png', dpi=700)
                print(f'Generating diagnostic outputs for su2 took {time()-t} seconds', flush=True)
                t = time()

            # ----------------------- vars(np, nvrt) and vars(np, nvrt, ntracers) ----------------------------
            trnd_tmp = hot_in.tr_nd.val[neighbor]
            row = np.r_[np.array(range(self.dims[0])), np.array(range(self.dims[0]))]
            for k in range(self.dims[3]):
                col = np.r_[grid_out.vgrid.z_idx_lower[:, k], grid_out.vgrid.z_idx_upper[:, k]]
                data = np.r_[grid_out.vgrid.z_weight_lower[:, k], 1.0 - grid_out.vgrid.z_weight_lower[:, k]]
                weights = csc_matrix((data, (row, col)), shape=(self.dims[0], hot_in.dims[3])).toarray()
                for var_str in self.vars_2d_node_based:
                    self.var_dict[var_str].val[:, k] = np.sum(hot_in.var_dict[var_str].val[neighbor] * weights, axis=1)

                for j in range(self.dims[4]):  # loop ntracers
                    self.tr_nd.val[:, k, j] = np.sum(trnd_tmp[:, :, j] * weights, axis=1)
                    pass

                print(f'Processing Layer {k+1} of {self.dims[3]} for all N-dimensional node-based variables took {time()-t} seconds', flush=True)
                t = time()
            self.tr_nd0.val = self.tr_nd.val[:]

            if iplot:
                plt.scatter(grid_out.hgrid.x, grid_out.hgrid.y, s=5, c=self.tr_nd.val[:, plot_layer, 0], cmap='jet', vmin=0, vmax=33)
                plt.savefig(f'{out_dir}/trnd_out2.png', dpi=700)
                print(f'Generating diagnostic outputs for trnd took {time()-t} seconds', flush=True)
                t = time()
            print(f'Total time for interpolation: {time()-t0} seconds', flush=True)

    def trnd_propogate(self):
        '''
        Propogate trnd to tr_el and tr_nd0
        '''
        tmp_grid = copy.deepcopy(self.grid.hgrid)
        for i, _ in enumerate(['tem', 'sal']):
            for j in range(0, self.grid.vgrid.nvrt):
                tmp_ele_vals = tmp_grid.interp_node_to_elem(value=self.tr_nd.val[:, j, i])
                self.tr_el.val[:, j, i] = copy.deepcopy(tmp_ele_vals)

        self.tr_nd0.val[:] = self.tr_nd.val[:]

    def writer(self, fname):
        WriteNC(fname, self)

    def replace_vars(self, var_dict=[], shapefile_name=None):
        '''
        var_dict:
                keys: variable names whose values are to be replaced;
                values: replacement values
        shapefile_name:
                contains one or more polygons inside which the value replacement will occur
        '''
        if not hasattr(self.grid, 'hgrid'):
            raise Exception('missing hgrid, initialize Hotstart instance with hgrid and try again.')
        if shapefile_name is not None:
            [ele_idx_list, node_idx_list] = find_ele_node_in_shpfile(
                # order matters, latter prevails
                shapefile_name=shapefile_name,
                grid=self.grid.hgrid,
            )

        for var in var_dict.keys():
            if self.var_dict[var].val.shape != var_dict[var].shape:
                raise Exception('inconsistent dimensions')
            if self.var_dict[var].dimname[0] == 'elem':
                for ind in ele_idx_list:
                    self.var_dict[var].val[ind] = var_dict[var][ind]
            elif self.var_dict[var].dimname[0] == 'node':
                for ind in node_idx_list:
                    self.var_dict[var].val[ind] = var_dict[var][ind]
            else:
                raise Exception(f'operation not implemented for dimension {self.var_dict[var].dimname}')


def gather_grid_info(grid_in, eta, dir):
    if not hasattr(grid_in, 'hgrid'):
        grid_in.hgrid = schism_grid(f'{dir}/hgrid.gr3')

    grid_in.hgrid.compute_ctr()
    grid_in.hgrid.compute_side(fmt=1)
    grid_in.hgrid.side_x = (grid_in.hgrid.x[grid_in.hgrid.isidenode[:, 0]] + grid_in.hgrid.x[grid_in.hgrid.isidenode[:, 1]]) / 2.0
    grid_in.hgrid.side_y = (grid_in.hgrid.y[grid_in.hgrid.isidenode[:, 0]] + grid_in.hgrid.y[grid_in.hgrid.isidenode[:, 1]]) / 2.0

    grid_in.vgrid = schism_vgrid()
    grid_in.vgrid.read_vgrid(f'{dir}/vgrid.in')

    grid_in.vgrid.zcor = grid_in.vgrid.compute_zcor(grid_in.hgrid.dp, eta=eta)

    grid_in.vgrid.zcor_s = (grid_in.vgrid.zcor[grid_in.hgrid.isidenode[:, 0]] + grid_in.vgrid.zcor[grid_in.hgrid.isidenode[:, 1]]) / 2.0
    grid_in.vgrid.kbp_s = np.min(grid_in.vgrid.kbp[grid_in.hgrid.isidenode[:, :]], axis=1)

    grid_in.vgrid.zcor_e = np.zeros((grid_in.hgrid.ne, grid_in.vgrid.nvrt))
    grid_in.vgrid.kbp_e = np.zeros((grid_in.hgrid.ne)).astype(int)
    for i in [3, 4]:
        II = (grid_in.hgrid.i34 == i)
        for j in range(i):
            grid_in.vgrid.zcor_e[II, :] += grid_in.vgrid.zcor[grid_in.hgrid.elnode[II, j], :] / i
        grid_in.vgrid.kbp_e[II] = np.min(grid_in.vgrid.kbp[grid_in.hgrid.elnode[II, :i]], axis=1)

    return grid_in


def GetVerticalWeight(zcor_in, kbp_in, zcor_out, neighbors):
    # at grid_out's node, save vertical interp weights
    z_weight_lower = 0.0 * zcor_out
    z_idx_lower = (0.0 * zcor_out).astype(int)
    z_idx_upper = (0.0 * zcor_out).astype(int)
    zcor_tmp = zcor_in[neighbors]

    n_points = zcor_out.shape[0]
    nvrt_in = zcor_in.shape[1]
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
        interior = (l_idx >= 0) & (l_idx < nvrt_in - 1)
        above = (l_idx == (nvrt_in - 1))

        z_weight_lower[i, interior] = (zcor_tmp[i, l_idx[interior] + 1] - zcor_out[i, interior]) / dz[i, l_idx[interior]]
        z_idx_lower[i, interior] = l_idx[interior]
        z_idx_upper[i, interior] = l_idx[interior] + 1

        z_weight_lower[i, below] = 0.0
        z_idx_lower[i, below] = kbp_in[neighbors[i]]
        z_idx_upper[i, below] = kbp_in[neighbors[i]]

        z_weight_lower[i, above] = 1.0
        z_idx_lower[i, above] = l_idx[above]
        z_idx_upper[i, above] = l_idx[above]

    return [z_weight_lower, z_idx_lower, z_idx_upper]


def find_ele_node_in_shpfile(shapefile_name, grid):
    '''
    Find element/node index within one or more polygons defined in a shapefile
        shapefile_name: file contains polygon(s)
        grid: schism_grid instance
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


def replace_hot_vars(infile, outfile, grid, vars_list=[], shapefile_name=None):
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
        grid.n_points = hot_out[var].data.shape[0]
        if grid.n_points == hot_in.dims['elem']:
            for ind in ele_idx_list:
                hot_out[var].data[ind] = hot_in[var].data[ind]
        elif grid.n_points == hot_in.dims['node']:
            for ind in node_idx_list:
                hot_out[var].data[ind] = hot_in[var].data[ind]
        else:
            raise Exception(f'unknown dimension {grid.n_points}')

    hot_out.to_netcdf(outfile + '.new')


if __name__ == "__main__":

    '''
    Sample usage:
    '''

    # # Sample 1: automatically fill missing vars with 0 for hotstart.nc of newer versions
    # my_hot = Hotstart(
    #     grid_info='/sciclone/data10/feiye/vims20/work/ChesBay/RUN110y/',
    #     hot_file='/sciclone/data10/feiye/vims20/work/ChesBay/RUN110y/hotstart.nc'
    # )
    # my_hot.writer(fname='/sciclone/data10/feiye/vims20/work/ChesBay/RUN110y/New_fmt_convert/new_hotstart.nc')

    # Sample 2: replacing variable values within a region
    # Sample 2.1: overwrite an existing hotstart.nc from another one on the same hgrid/vgrid
    replace_hot_vars(
        infile='/sciclone/schism10/feiye/From_Nabi/RUN02/Test_Hot/hotstart_it=31104.nc',
        outfile='/sciclone/schism10/feiye/From_Nabi/RUN02/Test_Hot/hotstart.nc.0',
        grid=schism_grid('/sciclone/schism10/feiye/From_Nabi/RUN02/hgrid.ll'),
        vars_list=['tr_nd', 'tr_nd0', 'tr_el'],
        shapefile_name='/sciclone/schism10/feiye/From_Nabi/RUN02/Test_Hot/ocean.shp'
    )

    # Sample 2.2: similar as 2.1, but using a function provided by the Hotstart class
    bg_hot = Hotstart(
        grid_info='/sciclone/schism10/feiye/ICOGS/v3_hot_20211115/',
        hot_file='/sciclone/schism10/feiye/ICOGS/v3_hot_20211115/hotstart.nc.hycom'
    )
    fg_hot = Hotstart(
        grid_info='/sciclone/schism10/feiye/ICOGS/v3_hot_20211115/',
        hot_file='/sciclone/schism10/feiye/ICOGS/v3_hot_20211115/hotstart.nc.combined_from_fcst'
    )
    fg_hot.replace_vars(
        var_dict={'tr_nd': bg_hot.tr_nd.val, 'tr_nd0': bg_hot.tr_nd0.val, 'tr_el': bg_hot.tr_el.val},
        shapefile_name='./ocean.shp',
    )  # sample *.shp is provided under the same folder as this script in SCHISM GIT

    # Sample 2.3 tweaking a single variable directly
    [_, node_idx_list] = find_ele_node_in_shpfile(
        shapefile_name="/sciclone/schism10/feiye/ICOGS/Ida04b/Elev_IC/city_polys_from_v10_lonlat.shp",
        grid=fg_hot.grid.hgrid
    )
    for ind in node_idx_list:
        fg_hot.eta2.val[ind] = -fg_hot.grid.hgrid.dp[ind] - 0.1  # set water surface to 0.1 m below ground
    # Note: do fg_hot.trnd_propogate() before writing if trnd is changed; it propogates trnd values to trnd0 and tr_el
    fg_hot.writer(f'{fg_hot.source_dir}/hotstart.nc')

    # Sample 3: interpolating one hotstart.nc to another
    hot_background = Hotstart(
        grid_info='/sciclone/schism10/feiye/Coastal_Act/RUN11c/',  # contains hgrid.gr3 and vgrid.in
        hot_file='/sciclone/schism10/feiye/Coastal_Act/RUN11c/hotstart.nc'
    )  # create a Hotstart instance with existing values
    my_hot = Hotstart(
        grid_info='/sciclone/schism10/feiye/Coastal_Act/Interp_Hot/11d/',
        ntracers=hot_background.dims[4]  # dims: [np, ne, ns, nvrt, ntracers, one]
    )  # create a Hotstart instance with empty values
    my_hot.interp_from_existing_hotstart(hot_in=hot_background, iplot=True, i_vert_interp=True)
    my_hot.writer(f'{my_hot.source_dir}/interp_hotstart.nc')
