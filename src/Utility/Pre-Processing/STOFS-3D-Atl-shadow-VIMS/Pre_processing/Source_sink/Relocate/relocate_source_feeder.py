#!/usr/bin/env python3

'''
Relocate sources/sinks based on the feeder channel information.
See __main__ for sample usage.
'''

import os
from pathlib import Path
import pickle

import json
import numpy as np
from scipy import spatial

from pylib import schism_grid
from pylib_experimental.schism_file import source_sink, TimeHistory
from pylib_essentials.schism_file import schism_bpfile
from pylib_essentials.utility_functions import inside_polygon

# Global Var
# [About the mandatory_sources_coor array]
# Some major rivers may not have a auto-arc/feeder channel,
# or the auto-arc/feeder channel doesn't match the location of the main stem.
# To find them, after generating the original sources/sinks,
# overlay sources (e.g., see how "original_sources.xy" is generated)
# with resolved channels (total_river_arcs_extra.map).
#
# If a major river source does not flow into any resolved channel,
# record them in the mandatory_sources_coor array.
# The first two columns are the lon/lat to be relocated to.
# The third and fourth columns are lon/lat of original source location.
#
# Some rivers are based on manual polygons so there is no nearby auto arcs.
# They also need to be recorded in the mandatory_sources_coor array,
# otherwise the sources will be relocated to the nearest auto arc.
# In this case, the original source is already in the correct location,
# ,so the third and fourth columns can be set as np.nan, which means they will take
# the values of the first and second columns.
v16_mandatory_sources_coor = np.array([
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

v19p2_mandatory_sources_coor = np.array([
    [-73.90869, 42.13509, np.nan, np.nan],  # Hudson River, NY
    [-74.94442, 40.34478, np.nan, np.nan],  # Delaware River, NJ
    [-76.132095, 39.603897, np.nan, np.nan],  # Susquehanna River, VA
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
    [-91.56184, 31.05043, np.nan, np.nan],  # Mississippi River
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

v19p2_for_sms_v27_mandatory_sources_coor = np.array([
    [-73.90869, 42.13509, np.nan, np.nan],  # Hudson River, NY
    [-74.94442, 40.34478, np.nan, np.nan],  # Delaware River, NJ
    [-76.1079760188,39.5871749053, np.nan, np.nan],  # Susquehanna River, VA
    [-78.425288, 34.508177, np.nan, np.nan],  # Cape Fear River, NC
    [-91.72306, 31.04462, np.nan, np.nan],  # Red River (upstream of Atchafalaya River), LA
    [-80.10808, 33.50005, np.nan, np.nan],  # Santee River, SC
    [-79.83866896,33.65649571, np.nan, np.nan],  # Black River, SC; note: channel not resolved
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
    [-91.56184, 31.05043, np.nan, np.nan],  # Mississippi River
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

v23p3_mandatory_sources_coor = np.array([
    [-73.90869, 42.13509, np.nan, np.nan],  # Hudson River, NY
    [-74.94442, 40.34478, np.nan, np.nan],  # Delaware River, NJ
    [-76.1011, 39.5802, np.nan, np.nan],  # Susquehanna River, VA
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
    [-91.56184, 31.05043, np.nan, np.nan],  # Mississippi River
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

v24p4_mandatory_sources_coor = np.array([
    [-73.90869, 42.13509, np.nan, np.nan],  # Hudson River, NY
    [-74.94442, 40.34478, np.nan, np.nan],  # Delaware River, NJ
    [-76.1011, 39.5802, np.nan, np.nan],  # Susquehanna River, VA
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
    [-91.56184, 31.05043, np.nan, np.nan],  # Mississippi River
    [-79.43425, 33.84487, -79.50974266666667, 33.85385866666667],  # Lyches River, SC
    [-74.74868, 39.47915, -74.75470666666668, 39.485390333333335],  # Great Egg Harbor River, NJ
    [-73.94009733333333, 42.06972966666667, np.nan, np.nan],  # Saugeties Creek, NY
    [-73.971293, 41.920595999999996, np.nan, np.nan],  # Hudson River branch, NY
    [-73.07229533333333, 41.303546000000004, np.nan, np.nan],  # Housatonic River, CT
    [-72.625735, 41.656137666666666, np.nan, np.nan],  # Connecticut River, CT
    [-72.64970633333333, 41.572111666666665, np.nan, np.nan],  # Mattabesset River, CT
    [-72.470818, 41.47020933333334, np.nan, np.nan],  # Salmon River, CT
    [-72.11158266666666, 41.455657333333335, np.nan, np.nan],  # Stony Brook, CT
    [-72.090553, 41.535118000000004, np.nan, np.nan],  # Yantic River, CT
    [-72.06195833333334, 41.525600000000004, np.nan, np.nan],  # Quinebaug River, CT
]).reshape(-1, 4)


def nearest_neighbour(points_a, points_b):
    '''A wrapper for scipy.spatial.cKDTree to find the nearest neighbour of points_a in points_b'''
    tree = spatial.cKDTree(points_b)
    return np.array(tree.query(points_a)[1]).reshape(-1,), np.array(tree.query(points_a)[0]).reshape(-1,)


def cal_dist(points_group_a, points_group_b):
    """Calculate the distance between two groups of points"""
    points_a = np.squeeze(points_group_a.view(np.complex128))
    points_b = np.squeeze(points_group_b.view(np.complex128))
    return np.absolute(points_a-points_b)


def lonlat2cpp(lon, lat, lon0=0, lat0=0):
    """Convert lon/lat to Cartesian coordinates in meters"""
    R_EARTH = 6378206.4

    lon_radian, lat_radian = lon/180*np.pi, lat/180*np.pi
    lon0_radian, lat0_radian = lon0/180*np.pi, lat0/180*np.pi

    xout = R_EARTH * (lon_radian - lon0_radian) * np.cos(lat0_radian)
    yout = R_EARTH * lat_radian

    return [xout, yout]


def cpp2lonlat(x, y, lon0=0, lat0=0):
    """
    Convert cpp in meters to lon/lat
    """
    R_EARTH = 6378206.4
    lon0_radian, lat0_radian = lon0/180*np.pi, lat0/180*np.pi
    lon_radian = x / R_EARTH / np.cos(lat0_radian) + lon0_radian
    lat0_radian = y / R_EARTH
    return [lon_radian/np.pi*180, lat0_radian/np.pi*180]


def relocate_sources2(
    old_ss_dir=None, outdir=None,
    no_feeder=False,
    feeder_info_file=None, hgrid_fname=None, allow_neglection=True,
    max_search_radius=1000.0, mandatory_sources_coor=np.empty((0, 4)),
    region_list=None,
):
    """
    Relocate sources to feeder channels of the resolved rivers

    :old_ss_dir: the directory of the original source_sink files,
        which should contain the sources.json file generated by pyschism
    :no_feeder: if true, feeder_bases (which should still be present in a mesh without feeders)
        will be used; otherwise, feeder_heads will be used.
    :outdir: the output directory
    :feeder_info_file: the feeder channel information generated by make_feeder_channel.py
    :hgrid_fname: the hgrid file name of the new grid with feeders
    :allow_neglection: Allow neglecting channels that cannot be matched to a new source location
        within a specified search radius.
    :max_search_radius: the maximum search radius to find the closest old source to a new source
    :mandatory_sources_coor: a 2D array (n, 4) of source locations
        that must be included in the relocation process.
    :region_list: a list of polygons, each specified by a 2D array (npts, 2) of coordinates.
        This is used to restrict the relocation to specific regions.
    """

    if (not allow_neglection) and region_list is not None:
        raise ValueError('must allow neglection if region_list is specified')

    # make the output directory if it doesn't exist
    os.makedirs(outdir, exist_ok=True)

    # -------------------------------------info from old source_sink-------------------------------------
    # Read original source/sink generated by pyschism based on NWM
    old_source_sink = source_sink.from_files(source_dir=old_ss_dir)
    old_gd = schism_grid(f'{old_ss_dir}/hgrid.gr3')
    old_gd.compute_ctr()  # computer element center coordinates

    # read source to fid mapping from json
    # assuming ordered dictionary keys, i.e., python >= 3.7
    with open(f'{old_ss_dir}/sources.json', 'r', encoding='utf-8') as f:
        old_sources2fids = json.load(f)

    # this may not be in the same order as in source_sink.in
    old_sources_eles = np.array(list(old_sources2fids.keys()), dtype=int)

    # get mean flow of old sources, using the same order as the keys from sources.json
    old_vsource_mean = old_source_sink.vsource.df[old_sources_eles.astype(str)].mean().values
    # get old sources' coordinates, in the order of sources.json's keys
    old_sources_coor = np.c_[
        old_gd.xctr[old_sources_eles-1], old_gd.yctr[old_sources_eles-1], old_vsource_mean
    ]
    # diagnostic output
    np.savetxt(f'{old_ss_dir}/original_sources.xy', old_sources_coor)

    # -------------------------------------info on feeders-------------------------------------
    # Read feeder channel info, which is generated by make_feeder_channel.py
    # ".pkl" and "*.xy" are both generated; the latter is preferred because it is human readable
    if Path(feeder_info_file).suffix == '.pkl':
        with open(feeder_info_file, 'rb') as file:
            [_, _, feeder_heads, feeder_bases] = pickle.load(file)  # neglected returns: feeder_l2g and feeder_points
        np.savetxt(f'{outdir}/feeder_heads_bases.xy', np.c_[feeder_heads[:, :2], feeder_bases[:, :2]])
    elif Path(feeder_info_file).suffix == '.xy':
        tmp = np.loadtxt(feeder_info_file)
        feeder_heads = tmp[:, :2]
        feeder_bases = tmp[:, 2:]
    else:
        raise ValueError('feeder info file must be .pkl or .xy')

    # -------------------------------------new hgrid-------------------------------------
    # get new hgrid (with feeders)
    new_gd = schism_grid(hgrid_fname)
    new_gd.compute_ctr()

    # -------------------- process nan values in mandatory_sources_coor ------------------
    # nan values means the mandatory relocation is the same as the old source location
    for i, row in enumerate(mandatory_sources_coor):
        if np.isnan(row[2]):
            mandatory_sources_coor[i, 2] = mandatory_sources_coor[i, 0]
        if np.isnan(row[3]):
            mandatory_sources_coor[i, 3] = mandatory_sources_coor[i, 1]

    # -------------------- restrict the relocation to regions if specified --------------------
    if region_list is not None:
        idx = np.zeros(len(feeder_heads), dtype=bool)
        for region in region_list:
            idx += inside_polygon(feeder_heads[:, :2], region[:, 0], region[:, 1]).astype(bool)
        if any(idx):  # subset if any feeder is in regions
            feeder_heads = feeder_heads[idx]
            feeder_bases = feeder_bases[idx]

        for region in region_list:
            idx = inside_polygon(mandatory_sources_coor[:, :2], region[:, 0], region[:, 1]).astype(bool)
        if any(idx):
            mandatory_sources_coor = mandatory_sources_coor[idx]

    # -------------------------------------relocation-------------------------------------
    # Find matching source point at mandatory_sources_coor and feeders
    # These are the desired new source locations
    if no_feeder:
        # If feeders are not implemented in the mesh, use "feeder_bases" instead of "feeder_heads",
        # because "feeder_bases" are the locations where the feeders would be placed.
        # The relocation still works in the sense that only resolved rivers receive new sources.
        new_sources_coor = np.r_[mandatory_sources_coor[:, :2], feeder_bases[:, :2]]
    else:
        new_sources_coor = np.r_[mandatory_sources_coor[:, :2], feeder_heads[:, :2]]
    # These are the locations used to search for the closest old source
    # ; in other words, to link the new source to the old source.
    # For the mandatory sources, these may be different from the new source locations
    # to allow for more reasonable matching between the old source location and the new source location
    # than the direct distance between the new sources and the old sources.
    # For example, it makes sense to use the feeder base to search for the closest old source,
    # but the new source should be placed at the feeder head.
    new_sources_search_point_coor = np.r_[mandatory_sources_coor[:, 2:4], feeder_bases[:, :2]]

    # get new source elements (ele_idx, 0-based indexing)
    new_sources_ele_idx, _ = nearest_neighbour(new_sources_coor, np.c_[new_gd.xctr, new_gd.yctr])
    new_sources_eles = new_sources_ele_idx + 1  # 1-based indexing

    # projection from lon/lat to meters
    old_sources_x, old_sources_y = lonlat2cpp(old_sources_coor[:, 0], lat=old_sources_coor[:, 1])
    # Note: new_sources_x and y are based on new_sources_search_point_coor
    new_sources_x, new_sources_y = lonlat2cpp(new_sources_search_point_coor[:, 0], new_sources_search_point_coor[:, 1])

    # link each new source to the largest old source within the search radius
    # mandatory sources are the first len(mandatory_sources_coor) new_sources
    new2old_sources = - np.ones(len(new_sources_x), dtype=int)
    relocation_distance = np.zeros(len(new_sources_x), dtype=float)
    # for each new source, find at most 1 old source within the search radius
    for i, [x, y] in enumerate(zip(new_sources_x, new_sources_y)):
        distance = cal_dist(np.c_[x, y], np.c_[old_sources_x, old_sources_y])
        valid_old_sources_idx = np.argwhere(distance < max_search_radius).reshape(-1,)
        if len(valid_old_sources_idx) > 0:
            # within valid_old_sources_idx, choose the one with the largest mean flow,
            # it is likely the correct old source because the new sources will be in larger channels
            target_old_source = valid_old_sources_idx[np.argmax(old_vsource_mean[valid_old_sources_idx])]
            # this holds the index of the old source in sources.json, not element idx or id
            new2old_sources[i] = target_old_source
            relocation_distance[i] = distance[target_old_source]
        else:
            if i < len(mandatory_sources_coor):  # still processing mandatory sources
                # mandatory sources are the first len(manadatory_sources_coor) new_sources,
                # which must be relocated
                raise ValueError(f'mandatory new source {i}: {cpp2lonlat(x, y)} cannot be mapped to an old source')

    # -------------------------------------clean up odd cases-------------------------------------
    # check if multiple new sources are linked to the same old source, i.e., double-counting an old source
    new2old_sources_copy = new2old_sources.copy()
    for old_source in np.unique(new2old_sources_copy):
        if old_source == -1:
            continue  # new source not linked to any old source, skip

        ids = np.argwhere(new2old_sources_copy == old_source).flatten()  # find all new sources mapped to the same old source
        if len(ids) == 0:  # impossible
            raise ValueError(f'Impossible: old source {old_source} in new2old mapping not mapped to any new source')
        if len(ids) > 1:  # multiple new sources mapped to the same old source, only retain the closest new source
            min_dist_id = np.argmin(relocation_distance[ids])
            new2old_sources[ids] = -1  # remove all corresponding new sources first
            new2old_sources[ids[min_dist_id]] = old_source  # reset the closest new source as the only one

    # add 1 to get element ids (1-based indexing)
    valid_new_sources_eleids = new_sources_eles[new2old_sources >= 0]  # 1-based indexing
    # 0-based index of the old source in sources.json, not element id
    valid_new2old_sources = new2old_sources[new2old_sources >= 0]

    # sanity check
    if np.unique(valid_new_sources_eleids).shape[0] != valid_new_sources_eleids.shape[0]:
        duplicates = np.unique(valid_new_sources_eleids)[np.unique(valid_new_sources_eleids, return_counts=True)[1] > 1]
        for duplicate in duplicates:
            print(f'duplicate new source: {duplicate} at {new_gd.xctr[int(duplicate)-1], new_gd.yctr[int(duplicate)-1]}')
        raise ValueError('Duplicated new sources can occur if a mandatory source has already been resolved by auto arcs')
    if np.unique(valid_new2old_sources).shape[0] != valid_new2old_sources.shape[0]:
        raise ValueError('Duplicated old sources')

    # At this point, all new sources are linked to at most one old source
    # , while valid_new_sources are linked to exactly one old source.
    # Build a dict with valid new source element ids as keys,
    # other potential new sources will be added to the dict if any
    # remaining old sources can be linked
    new2fid = {str(k): set() for k in valid_new_sources_eleids}
    recorded_fid = set()
    for new_source_eleid, old_source_idx in zip(valid_new_sources_eleids, valid_new2old_sources):
        # add the fids of the old source to the new source
        fids = old_sources2fids[str(old_sources_eles[old_source_idx])]
        # remove already assigned fids
        # duplicates can occur because multiple old sources can be linked to the same fid
        fids = [fid for fid in fids if fid not in recorded_fid]
        # some new sources have all fids already assigned to other new sources, i.e., fids ==[]
        if fids != []:
            new2fid[str(new_source_eleid)].update(fids)
            recorded_fid.update(fids)
        else:
            # remove the new source from new2fid if it has no fids
            new2fid.pop(str(new_source_eleid))

    # check for duplicates in new2fid
    tmp = [fid for fids in new2fid.values() for fid in fids]
    if len(tmp) != len(set(tmp)):
        raise ValueError('Duplicated fids in new2fid')

    # ------------------ find the remaining unassigned old sources ------------------
    if not allow_neglection:
        # indices (0-based) in the old source_sink.in
        unassigned_old_sources_idx = [i for i in range(len(old_sources_coor)) if i not in valid_new2old_sources]
        # find the closest new source to each unassigned old source;
        # all new source locations are considered, including unassigned ones, i.e., new2old_sources == -1
        unassigned2new_idx, _ = nearest_neighbour(
            np.c_[old_sources_x[unassigned_old_sources_idx], old_sources_y[unassigned_old_sources_idx]],
            np.c_[new_sources_x, new_sources_y])

        # test for duplicates after assigning the unassigned old sources
        for unassigned_old_source_idx, candidate_new_idx in zip(unassigned_old_sources_idx, unassigned2new_idx):
            old_fids = old_sources2fids[str(old_sources_eles[unassigned_old_source_idx])]
            old_fids = [fid for fid in old_fids if fid not in recorded_fid]  # remove already assigned fids

            if old_fids != []:
                target_newele_id = new_sources_eles[candidate_new_idx]  # 1-based indexing
                if str(target_newele_id) not in new2fid.keys():  # some new sources have not been linked to any old source
                    new2fid[str(target_newele_id)] = set(old_fids)
                else:
                    new2fid[str(target_newele_id)].update(old_fids)  # this prevents duplicates in each new source's fids
                recorded_fid.update(old_fids)  # update the recorded fids

    # check duplicated fids in new2fid
    tmp = [fid for fids in new2fid.values() for fid in fids]
    if len(tmp) != len(set(tmp)):
        raise ValueError('Duplicated fids in new2fid')

    # -------------------------------------write -------------------------------------
    # convert set to list for json serialization
    new2fid = {k: list(v) for k, v in new2fid.items()}
    json.dump(new2fid, open(f'{outdir}/sources.json', 'w', encoding='utf-8'), indent=4)
    # other diagnostic outputs
    with open(f'{outdir}/sources.txt', 'w', encoding='utf-8') as f:
        f.write('lon, lat, fids\n')
        for k, v in new2fid.items():
            f.write(f"{new_gd.xctr[int(k)-1]}, {new_gd.yctr[int(k)-1]}, {'-'.join(map(str, v))}\n")
    with open(f'{outdir}/old_sources.txt', 'w', encoding='utf-8') as f:
        f.write('lon, lat, fids\n')
        for k, v in old_sources2fids.items():
            f.write(f"{old_gd.xctr[int(k)-1]}, {old_gd.yctr[int(k)-1]}, {'-'.join(map(str, v))}\n")

    return new2fid  # same format as sources.json


def relocate_sources(
    old_ss_dir=None, outdir=None, relocate_map=None,
    no_feeder=False,
    allow_neglection=True, feeder_info_file=None, hgrid_fname=None,
    max_search_radius=1000.0, mandatory_sources_coor=np.empty((0, 4)),
    region_list=None,
):
    """
    Relocate sources to feeder channels of the resolved rivers

    :old_ss_dir: the directory of the original source_sink files
    :outdir: the output directory
    :relocate_map: a 2D array (n, 2) of new source element ids (1-based)
        and old source indices in source_sink.in (0 based)
        If relocate_map is provided, the relocation will be based on the map,
        which is useful for the operational forecast to save time,
        in this case no need to provide the arguments below;
        otherwise, the relocation map will be generated.
    :no_feeder: if true, feeder_bases (which should still be present in a mesh without feeders)
        will be used; otherwise, feeder_heads will be used.
    :allow_neglection: Allow neglecting channels that cannot be matched to a new source location.
        This is useful when the relocation is restricted to a specific region.
        If False, any unassigned old sources will be linked to the closest new source.
    :feeder_info_file:
        Feeder channel information generated by make_feeder_channel.py after RiverMapper
        This can be None if relocation_map is provided.
    :hgrid_fname: the hgrid file name of the new grid with feeders
    :mandatory_sources_coor:
        Major rivers that must be included in the relocation process.
        See more details in the global variables "*_mandatory_sources_coor" at the beginning of this script.
    :region_list: a list of polygons, each specified by a 2D array (npts, 2) of coordinates.
        This is used to restrict the relocation to specific regions.
    """

    if region_list is None:
        region_list = []

    # make the output directory if it doesn't exist
    os.makedirs(outdir, exist_ok=True)

    # Read original source/sink based on NWM
    old_source_sink = source_sink.from_files(source_dir=old_ss_dir)

    if relocate_map is None:  # if the relocation map is not provided, it needs to be generated here

        # -------------------------------------info from old source_sink-------------------------------------
        # get old hgrid (without feeders); hgrid with feeders may also be used
        old_gd = schism_grid(f'{old_ss_dir}/hgrid.gr3')

        # get old source coordinates
        old_gd.compute_ctr()
        # Note: source_eles starts from 1; mean vsource is added as the z field for diagnostics
        old_sources_coor = np.c_[old_gd.xctr[old_source_sink.source_eles-1], old_gd.yctr[old_source_sink.source_eles-1],
                                 old_source_sink.vsource.df.mean().values]
        np.savetxt(f'{old_ss_dir}/original_sources.xy', old_sources_coor)

        # -------------------------------------info on feeders-------------------------------------
        # Read feeder channel info, which is generated by make_feeder_channel.py
        # ".pkl" and "*.xy" are both generated; the latter is preferred because it is human readable
        if Path(feeder_info_file).suffix == '.pkl':
            with open(feeder_info_file, 'rb') as file:
                [_, _, feeder_heads, feeder_bases] = pickle.load(file)
            np.savetxt(f'{outdir}/feeder_heads_bases.xy', np.c_[feeder_heads[:, :2], feeder_bases[:, :2]])
        elif Path(feeder_info_file).suffix == '.xy':
            tmp = np.loadtxt(feeder_info_file)
            feeder_heads = tmp[:, :2]
            feeder_bases = tmp[:, 2:]
        else:
            raise ValueError('feeder info file must be .pkl or .xy')

        # -------------------------------------new hgrid-------------------------------------
        # get new hgrid (with feeders)
        new_gd = schism_grid(hgrid_fname)
        new_gd.compute_ctr()

        # -------------------------------------manual intervetion as pre-requisits -------------------------------------
        # process nan values in mandatory_sources_coor
        for i, row in enumerate(mandatory_sources_coor):
            if np.isnan(row[2]):
                mandatory_sources_coor[i, 2] = mandatory_sources_coor[i, 0]
            if np.isnan(row[3]):
                mandatory_sources_coor[i, 3] = mandatory_sources_coor[i, 1]

        # restrict the relocation to regions if specified
        if region_list != []:
            idx = np.zeros(len(feeder_heads), dtype=bool)
            for region in region_list:
                idx += inside_polygon(feeder_heads[:, :2], region[:, 0], region[:, 1]).astype(bool)
            if any(idx):  # subset if any feeder is in regions
                feeder_heads = feeder_heads[idx]
                feeder_bases = feeder_bases[idx]

            for region in region_list:
                idx = inside_polygon(mandatory_sources_coor[:, :2], region[:, 0], region[:, 1]).astype(bool)
            if any(idx):
                mandatory_sources_coor = mandatory_sources_coor[idx]

        # -------------------------------------relocation-------------------------------------
        # Find matching source point at mandatory_sources_coor and feeders
        # These are the desired new source locations
        if no_feeder:
            # If feeders are not implemented in the mesh, use "feeder_bases" instead of "feeder_heads",
            # because "feeder_bases" are the locations where the feeders are to be placed.
            # The relocation still works in the sense that only resolved rivers receive new sources.
            new_sources_coor = np.r_[mandatory_sources_coor[:, :2], feeder_bases[:, :2]]
        else:
            new_sources_coor = np.r_[mandatory_sources_coor[:, :2], feeder_heads[:, :2]]
        # These are the locations used to search for the closest old source
        # ; in other words, to link the new source to the old source.
        # For the mandatory sources, these may be different from the new source locations
        # to allow for more reasonable matching between the old source location and the new source location
        # than the direct distance between the new sources and the old sources.
        # For example, it makes sense to use the feeder base to search for the closest old source,
        # but the new source should be placed at the feeder head.
        new_sources_search_point_coor = np.r_[mandatory_sources_coor[:, 2:4], feeder_bases[:, :2]]

        # get new source elements (ele_idx, 0-based indexing)
        new_sources_ele_idx, _ = nearest_neighbour(new_sources_coor, np.c_[new_gd.xctr, new_gd.yctr])

        # projection from lon/lat to meters
        old_sources_x, old_sources_y = lonlat2cpp(old_sources_coor[:, 0], lat=old_sources_coor[:, 1])
        # Note: new_sources_x and y are based on new_sources_search_point_coor
        new_sources_x, new_sources_y = lonlat2cpp(
            new_sources_search_point_coor[:, 0], new_sources_search_point_coor[:, 1])

        # link each new source to the closest old source
        # mandatory sources are always relocated despite the distance
        new2old_sources = - np.ones(len(new_sources_x), dtype=int)
        relocation_distance = np.zeros(len(new_sources_x))
        old_vsource_mean = old_source_sink.vsource.df.mean().values
        for i, [x, y] in enumerate(zip(new_sources_x, new_sources_y)):
            # for each new source, find the old sources within the search radius
            distance = np.sqrt((old_sources_x - x)**2 + (old_sources_y - y)**2)
            valid_old_sources = np.argwhere(distance < max_search_radius).reshape(-1,)
            if len(valid_old_sources) > 0:
                # within valid_old_sources, choose the one with the largest mean flow,
                # it is likely the correct old source because the new sources are on larger channels
                target_old_source = valid_old_sources[np.argmax(old_vsource_mean[valid_old_sources])]
                # this holds the index of the old source in source_sink.in, not element idx or id
                new2old_sources[i] = target_old_source
                relocation_distance[i] = distance[target_old_source]
            else:
                if i < len(mandatory_sources_coor):
                    raise ValueError(f'mandatory new source {i}: {x, y} cannot be mapped to an old source')

        # -------------------------------------clean up odd cases-------------------------------------
        # check if multiple new sources are linked to the same old source, i.e., double-counting an old source
        for old_source in np.unique(new2old_sources):
            if old_source == -1:
                continue  # skip invalid new sources

            ids = np.argwhere(new2old_sources == old_source)  # find all new sources mapped to the same old source
            if len(ids) == 0:
                print(f'old source {old_source} cannot be mapped to a new source')
            elif len(ids) == 1:  # exact match
                pass
            else:  # multiple new sources mapped to the same old source, only retain the closest new source
                min_dist_id = np.argmin(relocation_distance[ids])
                new2old_sources[ids] = -1  # remove all corresponding new sources first
                new2old_sources[ids[min_dist_id]] = old_source  # keep the closest new source
        # add 1 to get element ids (1-based indexing)
        valid_new_sources_eleids = new_sources_ele_idx[new2old_sources >= 0] + 1
        # 0-based index of the old source in source_sink.in, not element id
        valid_new2old_sources = new2old_sources[new2old_sources >= 0]

        # sanity check
        if np.unique(valid_new_sources_eleids).shape[0] != valid_new_sources_eleids.shape[0]:
            raise ValueError('Multiple new sources')
        if np.unique(valid_new2old_sources).shape[0] != valid_new2old_sources.shape[0]:
            raise ValueError('Multiple old sources')

        # find the remaining unassigned old sources
        if not allow_neglection:
            # indices (0-based) in the old source_sink.in
            unassigned_old_sources_idx = [i for i in range(len(old_sources_coor)) if i not in valid_new2old_sources]
            # find the closest new source to each unassigned old source;
            # all new source locations are considered, including unassigned ones, i.e., new2old_sources == -1
            unassigned2new_idx, _ = nearest_neighbour(
                np.c_[old_sources_x[unassigned_old_sources_idx], old_sources_y[unassigned_old_sources_idx]],
                np.c_[new_sources_x, new_sources_y])

            # test for duplicates
            for unassigned_old_source_idx, candidate_new_idx in zip(unassigned_old_sources_idx, unassigned2new_idx):
                target_ele_id = new_sources_ele_idx[candidate_new_idx] + 1

                is_dup = False
                if target_ele_id in valid_new_sources_eleids:
                    # new source location already assigned, need to check for duplicates
                    dup_indices = np.argwhere(valid_new_sources_eleids == target_ele_id).reshape(-1)
                    for existing_idx in dup_indices:
                        vs_existing = old_source_sink.vsource.data[:, valid_new2old_sources[existing_idx]]
                        vs_test = old_source_sink.vsource.data[:, unassigned_old_source_idx]
                        if np.allclose(vs_test, vs_existing, atol=1e-6, rtol=1e-6):
                            print(f'unassigned old sources {unassigned_old_sources_idx}'
                                  f' is a duplicate of the new source {valid_new_sources_eleids[existing_idx]}')
                            is_dup = True
                            break

                if not is_dup:
                    # add the unassigned old sources to the valid relocation
                    # new_sources_ele_idx starts from 0, add 1 to get element ids (1-based indexing)
                    valid_new_sources_eleids = np.r_[valid_new_sources_eleids, target_ele_id]
                    valid_new2old_sources = np.r_[valid_new2old_sources, unassigned_old_source_idx]

        # -------------------------------------write relocation map-------------------------------------
        # this can also be used as an input for other runs with the same grid setup to bypass the relocation process
        np.savetxt(f'{outdir}/relocate_map.txt', np.c_[valid_new_sources_eleids, valid_new2old_sources], fmt='%d %d')

    else:
        valid_new_sources_eleids = relocate_map[:, 0]
        valid_new2old_sources = relocate_map[:, 1]

    # -------------------------------------write all outputs-------------------------------------
    # assemble source_sink object
    msource_list = []
    for old_msource in old_source_sink.msource:
        msource_list.append(
            TimeHistory(
                data_array=np.c_[old_msource.time, old_msource.data[:, valid_new2old_sources]],
                columns=['datetime'] + valid_new_sources_eleids.astype('str').tolist()
            )
        )

    vsource = TimeHistory(
        data_array=np.c_[old_source_sink.vsource.time, old_source_sink.vsource.data[:, valid_new2old_sources]],
        columns=['datetime'] + valid_new_sources_eleids.astype('str').tolist()
    )

    new_source_sink = source_sink(vsource=vsource, msource=msource_list, vsink=None)
    new_source_sink.writer(outdir)

    # Other diagnostic outputs. Note: source_eles starts from 1
    new_sources_coor = np.c_[
        new_gd.xctr[new_source_sink.source_sink_in.ip_group[0]-1],
        new_gd.yctr[new_source_sink.source_sink_in.ip_group[0]-1],
        vsource.df.mean().values
    ]
    np.savetxt(f'{outdir}/relocated_sources.xyz', new_sources_coor)

    return new_source_sink


def samples():
    '''
    Sample usage: relocate sources/sinks in a region of the Florence mesh.
    The mesh is same as the mesh of STOFS-3D-Atl v2.1 in the Florence region and coarser elsewhere.
    '''
    schism_git_dir = '/sciclone/data10/feiye/SCHISM_REPOSITORY/schism/'
    # inputs
    wdir = '/sciclone/schism10/feiye/Test_Project/Runs/R99a/Source_sink/'
    # this file is prepared during mesh generation
    # '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/Parallel/SMS_proj/feeder/feeder.pkl'
    feeder_info_file = (f'{schism_git_dir}/src/Utility/Pre-Processing/STOFS-3D-Atl-shadow-VIMS/'
                        'Pre_processing/Source_sink/feeder_heads_bases_v2.1.xy')
    region_file = (f'{schism_git_dir}/src/Utility/Pre-Processing/STOFS-3D-Atl-shadow-VIMS/'
                   'Pre_processing/Source_sink/relocate_florence.reg')
    old_ss_dir = f'{wdir}/original_source_sink/'
    hgrid_fname = f'{old_ss_dir}/hgrid.gr3'
    # end inputs

    # copy data files from git to the working directory for record if needed

    # read hgrid
    hgrid = schism_grid(hgrid_fname)

    # read region
    region = schism_bpfile()
    region.read_reg(fname=region_file)

    # read original source/sink
    original_ss = source_sink.from_files(source_dir=old_ss_dir)
    # original_ss.diag_writer(hgrid, old_ss_dir)

    # split source/sink into inside and outside region
    _, outside_ss = original_ss.clip_by_polygons(hgrid=hgrid, polygons_xy=[np.c_[region.x, region.y]])

    # relocate sources
    relocated_ss = relocate_sources(
        # strictly speaking, this should be the hgrid without feeders,
        # but an hgrid with feeders is also applicable with minor differences in the results
        old_ss_dir=old_ss_dir,
        feeder_info_file=feeder_info_file,
        hgrid_fname=hgrid_fname,  # strictly speaking, this is the with feeder hgrid
        outdir=f'{wdir}/relocated_source_sink/',
        max_search_radius=2000,  # search radius (in meters) for relocating sources
        mandatory_sources_coor=v16_mandatory_sources_coor,
        relocate_map=None,
        region_list=[np.c_[region.x, region.y]]
    )

    # combine outside and relocated sources
    combined_ss = outside_ss + relocated_ss
    combined_ss.writer(f'{wdir}/combined_source_sink/')
    # combined_ss.diag_writer(hgrid, f'{wdir}/combined_source_sink/')


def test():
    '''temporary tests'''
    wdir = '/sciclone/schism10/feiye/STOFS3D-v7/Inputs/I12z/Source_sink/relocated_source_sink/'
    relocate_sources2(
        old_ss_dir='/sciclone/schism10/feiye/STOFS3D-v7/Inputs/I12z/Source_sink/original_source_sink/',
        feeder_info_file='/sciclone/schism10/Hgrid_projects/STOFS3D-v7/v20.0/Feeder/feeder_heads_bases.xy',
        hgrid_fname=f'{wdir}/hgrid.gr3', outdir=wdir,
        max_search_radius=2100, mandatory_sources_coor=v23p3_mandatory_sources_coor,
        allow_neglection=False
    )


if __name__ == "__main__":
    test()
    print('Done!')
