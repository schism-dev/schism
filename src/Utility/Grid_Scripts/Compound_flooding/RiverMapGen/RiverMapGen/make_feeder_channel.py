import numpy as np
from pylib import schism_grid
from schism_py_pre_post.Grid.SMS import get_all_points_from_shp, \
    SMS_ARC, SMS_MAP, curvature, cpp2lonlat, lonlat2cpp, dl_cpp2lonlat, \
    get_perpendicular_angle, merge_maps
import pickle
from time import time

class Feeder():
    def __init__(self, points_x:np.ndarray, points_y:np.ndarray, base_id) -> None:
        self.points_x, self.points_y = points_x, points_y
        self.head = np.c_[np.mean(points_x[:, 0]), np.mean(points_y[:, 0])] 

        base_id = min(self.points_x.shape[1], base_id)
        self.base = np.c_[np.mean(points_x[:, base_id]), np.mean(points_y[:, base_id])] 

def find_inlets_outlets(line_map:SMS_MAP, gd:schism_grid):
    '''
    Find the index of the last inlet point (index is outside the grid, index-1 is inside) of all rivers.
    The projection of the river and the hgrid must be consistent.
    '''
    timer = time()

    if not hasattr(line_map, 'l2g') or not hasattr(line_map, 'xyz'):
        line_map.get_xyz()
    # ------------------------- find in-grid river points ---------------------------
    inside = gd.inside_grid(line_map.xyz[:, :2])
    # with open('inside.pkl', 'wb') as file:
    #     pickle.dump(inside, file)
    # with open('inside.pkl', 'rb') as file:
    #     inside = pickle.load(file)
    print(f'finding inside points took {time()-timer} seconds'); timer = time()

    inlets = -np.ones(len(line_map.l2g), dtype=int)
    outlets = -np.ones(len(line_map.l2g), dtype=int)
    for i, ids in enumerate(line_map.l2g):
        arc_points_inside = inside[ids]
        if np.any(arc_points_inside) and not np.all(arc_points_inside):  # intersecting with mesh boundary
            inlet = np.argwhere(arc_points_inside==0)[0][0]  # last inlet, assuming the river points are arranged from down to upstream
            if inlet > 0:
                inlets[i] = inlet
            outlet = np.argwhere(arc_points_inside)[0][0]  # last inlet, assuming the river points are arranged from down to upstream
            if outlet > 0:
                outlets[i] = outlet

    print(f'{len(line_map.l2g)} rivers, {sum(inlets>0)} inlets, {sum(outlets>0)} outlets')
    return inlets, outlets

if __name__ == "__main__":

    output_dir = '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/Parallel/SMS_proj/feeder/'

    timer = time()

    # river_centerlines_fname = f'{output_dir}/total_centerlines.map'
    rivermap_fname = f'{output_dir}/total_inner_arcs.map'
    grid_fname = '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/Parallel/SMS_proj/feeder/hgrid.ll'

    '''

    # river_map = SMS_MAP(rivermap_fname)
    # xyz, l2g = river_map.get_xyz()
    # with open('river_map.pkl', 'wb') as file:
    #     pickle.dump([river_map, xyz, l2g], file)

    with open('river_map.pkl', 'rb') as file:
        river_map, xyz, l2g = pickle.load(file)
    
    n_rivers = 0
    narcs_rivers = -np.ones((0, 1), dtype=int)  # number of river arcs for each river
    n = 0
    centerlines = []
    while (n < len(river_map.arcs)):
        # Number of river arcs (number of lines in the along-river direction)
        # is recorded in the integer part of the z-field of the river_map.
        # Get narcs from the first point of the first arc of a river
        narcs = int(river_map.arcs[n].points[0, -1].reshape(1, 1))

        # calculate centerline point coordinates by averaging river arcs
        centerline_points = river_map.arcs[n].points * 0.0
        for m in range(n, n+narcs):
            centerline_points += river_map.arcs[m].points
        centerline_points /= narcs
        centerlines.append(SMS_ARC(points=centerline_points, src_prj='epsg:4326'))

        n += narcs
        narcs_rivers = np.append(narcs_rivers, narcs)
        n_rivers += 1

    centerline_map = SMS_MAP(arcs=np.array(centerlines).reshape((-1, 1)))
    xyz_c, l2g_c = centerline_map.get_xyz()
        
    print(f'reading river map took {time()-timer} seconds'); timer = time()
    
    gd = schism_grid(grid_fname)
    # with open(grid_fname, 'rb') as file:
    #     gd = pickle.load(file)
    # gd.x, gd.y = gd.lon, gd.lat
    print(f'reading hgrid took {time()-timer} seconds'); timer = time()

    inlets, outlets = find_inlets_outlets(line_map=centerline_map, gd=gd)
    print(f'finding inlets/outlets took {time()-timer} seconds'); timer = time()

    with open('tmp.pkl', 'wb') as file:
        pickle.dump([river_map, centerline_map, n_rivers, narcs_rivers, gd, inlets, outlets], file)
    '''

    with open('tmp.pkl', 'rb') as file:
       [river_map, centerline_map, n_rivers, narcs_rivers, gd, inlets, outlets] = pickle.load(file)

    feeder_channel_extension = np.array([-10.0, -5.0, 0])
    i_inlet_options = [2, 1, 0]

    # get total number of river-arcs from all inlets
    n_inlets = sum(inlets > -1)
    narcs_rivers_inlets = sum(narcs_rivers[inlets>-1])

    # the number of feeder arcs = along-river arcs + cross-river arcs
    max_n_follow = 8
    feeder_arcs = np.empty((narcs_rivers_inlets + n_inlets * (max_n_follow + len(feeder_channel_extension)), 1), dtype=object)
    # the number of outlet arcs (only considering the base cross-river arc, since no pseudo channel is needed)
    # is set as the same of the number of inlets, but most inlets don't have an outlet
    outlet_arcs = np.empty((n_inlets, 1), dtype=object) # left bank and right bank for each thalweg

    i = 0; i_river = 0
    n_feeder = 0; n_outlet = 0
    feeders = []
    while (i < len(river_map.arcs)):
        this_inlets = inlets[i_river]  # inlets[range(i, i+narcs_rivers[i_river])] 
        if this_inlets > 0:  # any(this_inlets>0)
            for i_inlet_option in i_inlet_options:
                feeder_base_pts = np.zeros((0,3), dtype=float)
                feeder_follow_pts = np.zeros((0,3), dtype=float)

                inlet = min(this_inlets+i_inlet_option, len(river_map.arcs[i].points)-1) # max(0, np.min(this_inlets[this_inlets!=-1]) - 1)
                i_follow = range(inlet-1, max(inlet-max_n_follow, 0), -1)
                n_follow = len(i_follow)
                for j in range(i, i+narcs_rivers[i_river]):  # inner arcs of a river
                    feeder_base_pts = np.r_[feeder_base_pts, river_map.arcs[j].points[inlet, :].reshape(1,3)]
                    feeder_follow_pts = np.r_[feeder_follow_pts, river_map.arcs[j].points[i_follow, :].reshape(-1,3)]
                perp = np.mean(get_perpendicular_angle(line=feeder_base_pts[[1, -1], :2]))
                width = ((feeder_base_pts[0, 0] - feeder_base_pts[1, 0])**2 + (feeder_base_pts[0, 1] - feeder_base_pts[1, 1])**2) ** 0.5 
                feeder_channel_length = feeder_channel_extension * width

                xt = np.zeros((narcs_rivers[i_river], len(feeder_channel_extension)+n_follow), dtype=float)
                yt = np.zeros((narcs_rivers[i_river], len(feeder_channel_extension)+n_follow), dtype=float)
                for k in range(len(feeder_channel_extension)):
                    xt[:, k] = feeder_base_pts[:, 0] + feeder_channel_length[k] * np.cos(perp)
                    yt[:, k] = feeder_base_pts[:, 1] + feeder_channel_length[k] * np.sin(perp)
                
                ingrid_feeders = gd.inside_grid(np.c_[xt[:, i_inlet_option+1:len(feeder_channel_extension)].reshape(-1, 1),
                                                      yt[:, i_inlet_option+1:len(feeder_channel_extension)].reshape(-1, 1)])
                if sum(ingrid_feeders) == 0:
                    break  # found valid feeder; the worse case is non of the i_inlet_option works, in which case the last one is kept
                print(f'unclean connection at arc {i+1}')
                
            if n_follow > 0:
                xt[:, -n_follow:] = feeder_follow_pts[:, 0].reshape(-1, n_follow)
                yt[:, -n_follow:] = feeder_follow_pts[:, 1].reshape(-1, n_follow)
            
            for k in range(xt.shape[0]):
                feeder_arcs[n_feeder] = SMS_ARC(points=np.c_[xt[k, :], yt[k, :]], src_prj='epsg:4326')
                n_feeder += 1

            for k in range(xt.shape[1]):
                feeder_arcs[n_feeder] = SMS_ARC(points=np.c_[xt[:, k], yt[:, k]], src_prj='epsg:4326')
                n_feeder += 1

            feeders.append(Feeder(points_x=xt, points_y=yt, base_id=len(feeder_channel_extension)+1))

        this_outlets = outlets[i_river]  # outlets[range(i, i+narcs_rivers[i_river])] 
        if this_outlets>0: # any(this_outlets>0):
            outlet_base_pts = np.zeros((0,3), dtype=float)
            outlet = this_outlets  # max(0, np.min(this_outlets[this_outlets!=-1]))
            for j in range(i, i+narcs_rivers[i_river]):  # inner arcs of a river
                outlet_base_pts = np.r_[outlet_base_pts, river_map.arcs[j].points[outlet, :].reshape(1,3)]
                outlet_arcs[n_outlet] = SMS_ARC(points=outlet_base_pts, src_prj='epsg:4326')
                n_outlet += 1

        i += narcs_rivers[i_river]
        i_river += 1
    
    if len(feeders) != n_inlets:
        raise Exception("Inconsistent number of inlets and feeder channels")

    SMS_MAP(arcs=feeder_arcs.reshape((-1, 1))).writer(filename=f'{output_dir}/feeders.map')
    SMS_MAP(arcs=outlet_arcs.reshape((-1, 1))).writer(filename=f'{output_dir}/outlets.map')

    # save feeders info in a *.pkl
    feeder_heads = np.zeros((len(feeders), 3), dtype=float)
    feeder_bases = np.zeros((len(feeders), 3), dtype=float)
    feeder_l2g = [None] * len(feeders); npts = 0
    feeder_points = np.empty((0, 2), dtype=float)
    feeder_arrays_x = [None] * len(feeders)
    feeder_arrays_y = [None] * len(feeders)
    for i, feeder in enumerate(feeders):
        feeder_heads[i, :2] = feeder.head[:]
        feeder_bases[i, :2] = feeder.base[:]
        feeder_l2g[i] = np.array(np.arange(npts, npts+feeder.points_x.size))
        feeder_points = np.r_[feeder_points, np.c_[feeder.points_x.reshape(-1,1), feeder.points_y.reshape(-1,1)]]
        feeder_arrays_x[i] = feeder.points_x
        feeder_arrays_y[i] = feeder.points_y
        npts += feeder.points_x.size

    with open(f'{output_dir}/feeder.pkl', 'wb') as file:
        pickle.dump([feeder_l2g, feeder_points, feeder_heads, feeder_bases], file)

    with open(f'{output_dir}/feeder_arrays.pkl', 'wb') as file:
        pickle.dump([feeder_arrays_x, feeder_arrays_y], file)
    
    SMS_MAP(detached_nodes=np.r_[feeder_heads, feeder_bases]).writer(filename=f'{output_dir}/feeder_hb.map')
    
    pass

