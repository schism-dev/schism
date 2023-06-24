from schism_py_pre_post.Grid.SourceSinkIn import source_sink, SourceSinkIn
from schism_py_pre_post.Grid.SMS import lonlat2cpp
from schism_py_pre_post.Timeseries.TimeHistory import TimeHistory
from schism_py_pre_post.Grid.Hgrid_extended import read_schism_hgrid_cached
import numpy as np
from scipy import spatial
import pickle


def nearest_neighbour(points_a, points_b):
    tree = spatial.cKDTree(points_b)
    return np.array(tree.query(points_a)[1]).reshape(-1,), np.array(tree.query(points_a)[0]).reshape(-1,)

def dist(points_group_A, points_group_B):
    points_A = np.squeeze(points_group_A.view(np.complex128))
    points_B = np.squeeze(points_group_B.view(np.complex128))
    return np.absolute(points_A-points_B)


def relocate_sources(
    old_ss_dir=None,
    feeder_info_file=None,
    hgrid_fname=None,
    outdir=None,
    max_search_radius = 1000.0,
    mandatory_sources_coor = np.empty((0, 4)),
    relocate_map = None,
    #--------------------------- end inputs -------------------------
):

    # Read original source/sink based on NWM
    old_source_sink = source_sink(source_dir=old_ss_dir)

    if relocate_map is None:
        # get old hgrid (without feeders)
        old_gd = read_schism_hgrid_cached(f'{old_ss_dir}/hgrid.gr3', overwrite_cache=False)

        # get old source coordinates
        old_gd.compute_ctr()

        # Note: source_eles starts from 1
        old_sources_coor = np.c_[old_gd.xctr[old_source_sink.source_eles-1], old_gd.yctr[old_source_sink.source_eles-1],
                                old_source_sink.vsource.get_time_average([])]
        np.savetxt(f'{old_ss_dir}/original_sources.xy', old_sources_coor)

        # Read feeder channel info
        with open(feeder_info_file, 'rb') as file:
            [feeder_l2g, feeder_points, feeder_heads, feeder_bases] = pickle.load(file)
        np.savetxt(f'{outdir}/feeder_heads.xy', feeder_heads)
        np.savetxt(f'{outdir}/feeder_bases.xy', feeder_bases)

        # get new hgrid (with feeders)
        new_gd = read_schism_hgrid_cached(hgrid_fname, overwrite_cache=False)
        new_gd.compute_ctr()


        # process nan values in mandatory_sources_coor
        for i, row in enumerate(mandatory_sources_coor):
            if np.isnan(row[2]):
                mandatory_sources_coor[i, 2] = mandatory_sources_coor[i, 0]
            if np.isnan(row[3]):
                mandatory_sources_coor[i, 3] = mandatory_sources_coor[i, 1]

        # find matching source point at mandatory_sources_coor and feeders' base points
        # these are the desired new source locations
        new_sources_coor = np.r_[mandatory_sources_coor[:, :2], feeder_heads[:, :2]]
        # these are the locations used to search for the closest old source
        # , in other words, to link the new source to the old source
        new_sources_search_point_coor = np.r_[mandatory_sources_coor[:, 2:4], feeder_bases[:, :2]]

        # get new source location (ele id, index starting from 0)
        new_sources_eleids, _ = nearest_neighbour(new_sources_coor, np.c_[new_gd.xctr, new_gd.yctr])

        # projection from lon/lat to meters
        old_sources_x, old_sources_y = lonlat2cpp(old_sources_coor[:, 0], lat=old_sources_coor[:, 1])
        new_sources_x, new_sources_y = lonlat2cpp(new_sources_search_point_coor[:, 0], lat=new_sources_search_point_coor[:, 1])

        # link each new source to the closest old source
        # mandatory sources are always relocated despite the distance
        new2old_sources = - np.ones(len(new_sources_x), dtype=int)
        relocation_distance = np.zeros(len(new_sources_x))
        for i, [x, y] in enumerate(zip(new_sources_x, new_sources_y)):
            # find the distance of each new source to all old sources
            dist = np.sqrt((old_sources_x - x)**2 + (old_sources_y - y)**2)
            valid_old_sources = np.argwhere(dist < max_search_radius).reshape(-1,)
            if len(valid_old_sources) > 0:
                # take the old source with the largest vsource within the search radius
                valid_old_vsources = old_source_sink.vsource.get_time_average(valid_old_sources)
                target_old_source = valid_old_sources[np.argmax(valid_old_vsources)]
                new2old_sources[i] = target_old_source
                relocation_distance[i] = dist[target_old_source]
            else:
                if i < len(mandatory_sources_coor):  # mandatory sources must be relocated
                    raise ValueError(f'mandatory new source {i}: {x, y} cannot be mapped to an old source')

        for old_source in np.unique(new2old_sources):
            if old_source == -1: continue  # skip invalid new sources

            ids = np.argwhere(new2old_sources == old_source)  # find all new sources mapped to the same old source
            if len(ids) == 0:
                print(f'old source {old_source} cannot be mapped to a new source')
            elif len(ids) == 1:  # exact match
                pass
            else:  # multiple new sources mapped to the same old source, pick the closest new source
                min_dist_id = np.argmin(relocation_distance[ids])
                new2old_sources[ids] = -1
                new2old_sources[ids[min_dist_id]] = old_source

        valid_relocation = new2old_sources >= 0
        valid_new_sources_eleids = new_sources_eleids[valid_relocation] + 1
        valid_new2old_sources = new2old_sources[valid_relocation]
        np.savetxt(f'{outdir}/relocate_map.txt', np.c_[valid_new_sources_eleids, valid_new2old_sources], fmt='%d %d')
    else:
        valid_new_sources_eleids = relocate_map[:, 0]; valid_new2old_sources = relocate_map[:, 1]

    # Assemble new source/sink files
    nsources = len(valid_new_sources_eleids)

    source_sink_in = SourceSinkIn(filename=None, number_of_groups=2, ele_groups=[valid_new_sources_eleids.tolist(),[]])
    vsource = TimeHistory(
        file_name=None,
        data_array=np.c_[old_source_sink.vsource.time, old_source_sink.vsource.data[:, valid_new2old_sources]],
        columns=['datetime'] + valid_new_sources_eleids.astype('str').tolist()
    )
    msource = TimeHistory(
        file_name=None,
        data_array=np.c_[
            old_source_sink.msource.time,
            -9999*np.ones([old_source_sink.msource.n_time, nsources]),
            np.zeros([old_source_sink.msource.n_time, nsources])
        ],
        columns=['datetime'] + valid_new_sources_eleids.astype('str').tolist() + valid_new_sources_eleids.astype('str').tolist()
    )
    # Note: source_eles starts from 1
    new_sources_coor = np.c_[
        new_gd.xctr[source_sink_in.ip_group[0]-1],
        new_gd.yctr[source_sink_in.ip_group[0]-1],
        vsource.get_time_average([])
    ]
    np.savetxt(f'{outdir}/relocated_sources.xyz', new_sources_coor)

    source_sink_in.writer(f'{outdir}/source_sink.in')
    vsource.writer(f'{outdir}/vsource.th')
    msource.writer(f'{outdir}/msource.th')

if __name__ == "__main__":
    #--------------------------- inputs -------------------------
    old_ss_dir = '../original_source_sink/'
    feeder_info_file = f'/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/Parallel/SMS_proj/feeder/feeder.pkl'
    hgrid_fname = './hgrid.gr3'
    outdir = './'
    max_search_radius = 2000.0

    # Some major rivers may not have a feeder channel, or the feeder channel doesn't match the main river inputs.
    # In such cases, manually relocate sources with the mandatory_sources_coor array.
    # (e.g., by visualizing the map generated by see viz_source.py).
    # The first two columns are longitude and latitude of the mandatory (relocated) source locations.
    # The third and fourth columns are longitude and latitude of the corresponding original source locations.
    # If the target original source locations are the same as the mandatory source locations,
    # the third and fourth columns can be set as np.nan
    mandatory_sources_coor = np.array([
        [-69.77256, 44.31494, np.nan, np.nan],  # Kennebec River, ME
        [-73.90869, 42.13509, np.nan, np.nan],  # Hudson River, NY
        [-76.12905866666667, 39.60955966666666, np.nan, np.nan],  # Susquehanna River, VA
        [-74.94442, 40.34478, np.nan, np.nan],  # Delaware River, NJ
        [-78.425288, 34.508177, np.nan, np.nan],  # Cape Fear River, NC
        [-91.72306, 31.04462, np.nan, np.nan],  # Red River (upstream of Atchafalaya River), LA
        [-80.10808, 33.50005, np.nan, np.nan],  # Santee River, SC
        [-79.81703, 33.59694, np.nan, np.nan],  # Black River, SC
        [-79.57210, 33.71223, np.nan, np.nan],  # Black Mingo Creek, SC
        [-79.49997, 33.84686, np.nan, np.nan],  # Lynches River, SC
        [-79.48467, 33.93939, np.nan, np.nan],  # Pee Dee River, SC
        [-79.33247, 33.98196, np.nan, np.nan],  # Little Pee Dee River, SC
        [-77.917829, 34.749979, np.nan, np.nan],  # Northeast Cape Fear River, NC
        [-87.9523, 30.8472, np.nan, np.nan],  # Mobile River, AL
        [-96.695401, 28.968284, -96.69652166667, 28.990345],  # Lavaca River, TX
        [-96.548436, 28.999706, -96.554498, 29.024612666667],  # Lake Texana, TX
        [-93.83342666667, 30.355123333333, -93.83342666667, 30.355123333333],  # Cypress Creek, TX
        [-89.764476, 30.551926, -89.76781133333, 30.538070666667],  # Lotts Creek, LA
        [-87.219805, 30.567296, -87.24471466667, 30.601442333333],  # Escambia River, FL
        [-83.987035, 30.331327, np.nan, np.nan],  # Horsehead Creek and Little River, FL
        [-83.928038, 30.30404, np.nan, np.nan],  # Bailey Mill Creek, FL
        [-82.950913, 29.958097, -82.99605566667, 30.007415],  # Suwannee River, FL
        [-81.02370433333333, 27.315079666666666, np.nan, np.nan],  # Kissimmee River, FL
        [-81.997572, 30.786870, -82.040457, 30.74494233333333],  # St Marys River, FL
        [-79.43425, 33.84487, -79.50974266666667, 33.85385866666667],  # Lyches River, SC
        [-74.74868, 39.47915, -74.75470666666668, 39.485390333333335],  # Great Egg Harbor River, NJ
        [-73.94009733333333, 42.06972966666667, np.nan, np.nan],  # Saugeties Creek, NY
        [-73.971293, 41.920595999999996, np.nan, np.nan],  # Hudson River branch, NY
        [-73.92918633333333, 41.592421333333334, np.nan, np.nan],  # Hudson River branch, NY
        [-73.07229533333333, 41.303546000000004, np.nan, np.nan],  # Housatonic River, CT
        [-72.625735, 41.656137666666666, np.nan, np.nan],  # Connecticut River, CT
        [-72.64970633333333, 41.572111666666665, np.nan, np.nan],  # Mattabesset River, CT
        [-72.470818, 41.47020933333334, np.nan, np.nan],  # Salmon River, CT
        [-72.11158266666666, 41.455657333333335, np.nan, np.nan],  # Stony Brook, CT
        [-72.090553, 41.535118000000004, np.nan, np.nan],  # Yantic River, CT
        [-72.06195833333334, 41.525600000000004, np.nan, np.nan],  # Quinebaug River, CT
    ]).reshape(-1, 4)

    relocate_map = None  # if not None, use the prepared relocate_map to relocate sources

    relocate_sources(
        old_ss_dir=old_ss_dir,
        feeder_info_file=feeder_info_file,
        hgrid_fname=hgrid_fname,
        outdir=outdir,
        max_search_radius=max_search_radius,
        mandatory_sources_coor=mandatory_sources_coor,
        relocate_map=relocate_map
    )
