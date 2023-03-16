"""
This is the main worker script for generating river maps

Usage: 
Import and call the function "make_river_map",
either directly, see sample_serial.py in the installation directory;
or via the parralel driver, see sample_parallel.py in the installation directory
"""


import os
import pandas as pd
from copy import deepcopy
import numpy as np
from scipy import interpolate
import math
from shapely.geometry import LineString, Point
from shapely.ops import polygonize, unary_union, split, snap
import geopandas as gpd
from RiverMapper.SMS import get_all_points_from_shp, curvature, get_perpendicular_angle
from RiverMapper.SMS import SMS_ARC, SMS_MAP, cpp2lonlat, lonlat2cpp, dl_cpp2lonlat, dl_lonlat2cpp
from RiverMapper.river_map_tif_preproc import Tif2XYZ, get_elev_from_tiles


np.seterr(all='raise')


class Bombs():
    def __init__(self, xyz: np.ndarray, crs='epsg:4326'):
        if xyz.size == 0:
            raise ValueError('Input array empty')
        elif xyz.shape[1] != 3:
            raise ValueError('Input array must have 3 columns, corresponding to x, y, z (bomb range)')

        self.i_cleaned = False  # if the bomb points are cleaned (removing points too close to each other)
        self.points = np.squeeze(deepcopy(xyz[:, :2]).view(np.complex128))
        self.res = xyz[:, 2]
        self.crs = crs

    def __add__(self, other):
        if other is not None:
            self.points = np.r_[self.points, other.points]
            self.res = np.r_[self.res, other.res]
            self.i_cleaned *= other.i_cleaned
        return self

    def clean(self, bomb_radius_coef=0.5):
        i_valid = np.ones((len(self.points),), dtype=bool)
        for i, _ in enumerate(self.points):
            if i_valid[i]:
                dist = abs(self.points[i] - self.points)
                nearby_idx = np.where(dist<bomb_radius_coef*self.res[i])[0]
                nearby_idx = nearby_idx[nearby_idx!=i]
                i_valid[nearby_idx] = False

        if sum(i_valid) == 0:
            raise ValueError("impossible: all points bombed.")

        self.points = self.points[i_valid]
        self.res = self.res[i_valid]
        self.i_cleaned = True

        return np.c_[self.points.real, self.points.imag, self.res]
    
    def get_convex_hull(self):
        ch = gpd.GeoSeries(map(Point, np.c_[self.points.real, self.points.imag]), crs=self.crs).unary_union.convex_hull
        ch_buffered = ch.buffer(distance=0.3*np.mean(self.res))
        return ch_buffered

class Geoms_XY():
    def __init__(self, geom_list, crs='epsg:4326'):
        try:
            if not isinstance(geom_list[-1], LineString):
                raise TypeError()
        except:
            raise ValueError('Input must be an iterable of geometries')
        
        self.geom_list = [geom for geom in geom_list]
        self.xy, self.xy_idx = self.geomlist2xy()
        self.crs = crs

    def geomlist2xy(self):
        geoms_xy_idx = np.ones((len(self.geom_list), 2), dtype=int) * -9999

        idx = 0
        for i, geom in enumerate(self.geom_list):
            geoms_xy_idx[i, 0] = idx
            geoms_xy_idx[i, 1] = idx + len(geom.xy[0])
            idx += len(geom.xy[0])

        geoms_xy = np.empty((geoms_xy_idx[-1, 1], 2), dtype=float)
        for i, geom in enumerate(self.geom_list):
            geoms_xy[geoms_xy_idx[i, 0]:geoms_xy_idx[i, 1], :] = np.array(geom.xy).T

        return geoms_xy, geoms_xy_idx

    def xy2geomlist(self):
        for i, _ in enumerate(self.geom_list):
            self.geom_list[i] = LineString(self.xy[self.xy_idx[i, 0]:self.xy_idx[i, 1], :])

    def snap_to_points(self, snap_points, target_poly_gdf=None):
        geoms_xy_gdf = points2GeoDataFrame(self.xy, crs=self.crs)
        
        if target_poly_gdf is None:
            i_target_points = np.ones((self.xy.shape[0], ), dtype=bool)
        else:
            pointInPolys = gpd.tools.sjoin(geoms_xy_gdf, target_poly_gdf, predicate="within", how='left')
            _, idx = np.unique(pointInPolys.index, return_index=True)  # some points belong to multiple polygons
            i_target_points = np.array(~np.isnan(pointInPolys.index_right))[idx]

        _, idx = nearest_neighbour(self.xy[i_target_points, :], np.c_[snap_points[:, 0], snap_points[:, 1]])
        self.xy[i_target_points, :] = snap_points[idx, :2]

        self.xy2geomlist()
    
def moving_average(a, n=10, self_weights=0):
    if len(a) <= n:
        ret2 = a * 0.0 + np.mean(a)
        return ret2
    else:
        ret = np.cumsum(a, axis=0, dtype=float)
        ret[n:] = ret[n:] - ret[:-n]
        ret[n-1:] = ret[n-1:] / n

        # re-align time series
        ret1 = ret * 0.0
        m = int(np.floor(n/2))
        ret1[m:-m] = ret[2*m:]

        # fill the first and last few records
        ret1[:m] = ret1[m]
        ret1[-m:] = ret1[-m-1]

        # put more weights on self
        ret2 = (ret1 + self_weights * a) / (1 + self_weights)

        return ret2

def getAngle(a, b, c):
    ang = math.degrees(math.atan2(c[1]-b[1], c[0]-b[0]) - math.atan2(a[1]-b[1], a[0]-b[0]))
    return ang + 360 if ang < 0 else ang

def get_angle_diffs(xs, ys):
    line = np.c_[xs, ys]
    line_cplx = np.squeeze(line.view(np.complex128))
    angles = np.angle(np.diff(line_cplx))
    angle_diff0 = np.diff(angles)
    angle_diff = np.diff(angles)
    angle_diff[angle_diff0 > np.pi] -= 2 * np.pi
    angle_diff[angle_diff0 < -np.pi] += 2 * np.pi

    return angle_diff

def geos2SmsArcList(geoms=None):
    sms_arc_list = []
    for i, line in enumerate(geoms):
        sms_arc_list.append(SMS_ARC(points=np.c_[line.xy[0], line.xy[1]], src_prj='epsg:4326'))
    
    return sms_arc_list

def nearest_neighbour(points_a, points_b):
    from scipy import spatial
    tree = spatial.cKDTree(points_b)
    return tree.query(points_a)[0], tree.query(points_a)[1]

def split_line_by_point(line, point, tolerance: float=1.0e-12):
    # return split(snap(line, point, tolerance), point)
    return split(line, point)

def snap_vertices(line, thalweg_resolution):
    dist_along_thalweg = get_dist_increment(line)

    idx = 0
    original_seg_length = dist_along_thalweg[0]
    while idx < len(dist_along_thalweg)-1:
        if original_seg_length < thalweg_resolution[idx]:
            line[idx+1, :] = line[idx, :]  # snap current point to the previous one
            original_seg_length += dist_along_thalweg[idx+1]
        else:
            original_seg_length = dist_along_thalweg[idx+1]
        idx += 1  # move original arc forward

    return line

def get_dist_increment(line):
    line_copy = deepcopy(line)
    line_cplx = np.squeeze(line_copy.view(np.complex128))
    dist = np.absolute(line_cplx[1:] - line_cplx[:-1])

    # return np.r_[0.0, np.cumsum(dist)]
    return dist

def ccw(A,B,C):
    # A is a point with the coordinate x=A[0], y=A[1]
    return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0])

# Return true if line segments AB and CD intersect
def intersect(A,B,C,D):
    return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)

def smooth_thalweg(line, ang_diff_shres=np.pi/2.4, nmax=100, smooth_coef=0.2):
    xs = line[:, 0]
    ys = line[:, 1]

    n = 0
    while n < nmax:
        angle_diffs = np.r_[0.0, get_angle_diffs(xs, ys), 0.0]
        sharp_turns = np.argwhere(abs(angle_diffs) > ang_diff_shres)[:, 0]
        if not sharp_turns.size:
            break
        else:
            # Step 1: plan to move the sharp turn point to a new point, based on
            # the average coordinates of the two adjacent points and the original location
            line[sharp_turns, 0] = np.array((xs[sharp_turns-1] + xs[sharp_turns+1]) / 2) * smooth_coef + xs[sharp_turns] * (1-smooth_coef)
            line[sharp_turns, 1] = np.array((ys[sharp_turns-1] + ys[sharp_turns+1]) / 2) * smooth_coef + ys[sharp_turns] * (1-smooth_coef)

            n += 1

    if n == nmax:
        print(f'warning: smooth_thalweg did not converge in {n} steps\n')
    else:
        print(f'smooth_thalweg converged in {n} steps\n')

    perp = get_perpendicular_angle(line)

    return line, perp

def river_quality(xs, ys, idx, ang_diff_shres=np.pi/2.4):
    # identify channels that are too ambiguous
    if sum(idx) < 2:
        return False

    for [x, y] in zip(xs, ys):
        angle_diffs = abs(np.r_[0.0, get_angle_diffs(x[idx], y[idx]), 0.0])
        if np.sum(angle_diffs)/len(x[idx]) > 0.8:
            print(f'discarding arc based on the number of sharp turns\n')
            return False

    return True

def smooth_bank(line, xs, ys, xs_other_side, ys_other_side, ang_diff_shres=np.pi/2.4, nmax=100, smooth_coef=0.2):
    n = 0
    while n < nmax:
        angle_diffs = np.r_[0.0, get_angle_diffs(xs, ys), 0.0]
        sharp_turns = np.argwhere(abs(angle_diffs) > ang_diff_shres)[:, 0]

        if len(sharp_turns)==0:
            break
        else:
            # Step 1: plan to move the sharp turn point to a new point, based on
            # the average coordinates of the two adjacent points and the original location
            xs_moved = np.array((xs[sharp_turns-1] + xs[sharp_turns+1]) / 2) * smooth_coef + xs[sharp_turns] * (1-smooth_coef)
            ys_moved = np.array((ys[sharp_turns-1] + ys[sharp_turns+1]) / 2) * smooth_coef + ys[sharp_turns] * (1-smooth_coef)

            # Step 2: decide if the planned move is too far (i.e., across the thalweg)
            insert_turns_idx = []
            for i, [x_moved, y_moved, sharp_turn] in enumerate(zip(xs_moved, ys_moved, sharp_turns)):
                invalid_move = intersect([xs[sharp_turn], ys[sharp_turn]], [x_moved, y_moved], line[sharp_turn, :], line[sharp_turn-1, :]) + \
                               intersect([xs[sharp_turn], ys[sharp_turn]], [x_moved, y_moved], line[sharp_turn, :], line[sharp_turn+1, :])
                if invalid_move:
                    # prepare to insert 2 more points around the sharp turn
                    insert_turns_idx.append(i)
                else:
                    xs[sharp_turns] = xs_moved
                    ys[sharp_turns] = ys_moved

            if len(insert_turns_idx) > 0:
                # replace the point at sharp bend with two more points adjacent to it
                idx = sharp_turns[insert_turns_idx]

                tmp = np.c_[line, xs, ys, xs_other_side, ys_other_side]
                tmp_original = deepcopy(tmp)
                tmp[idx, :] = (tmp_original[idx, :] + tmp_original[idx+1, :])/2
                tmp = np.insert(tmp, idx, (tmp_original[idx, :] + tmp_original[idx-1, :])/2, axis=0)

                line = tmp[:, :2]
                xs = tmp[:, 2]
                ys = tmp[:, 3]
                xs_other_side = tmp[:, 4]
                ys_other_side = tmp[:, 5]

            n += 1

    if n == nmax:
        print(f'warning: smooth_bank did not converge in {n} steps\n')
    else:
        print(f'smooth_bank converged in {n} steps\n')

    perp = get_perpendicular_angle(line)

    return line, xs, ys, xs_other_side, ys_other_side, perp

def nudge_bank(line, perp, xs, ys, dist=None):
    if dist is None:
        dist=np.array([35, 500])

    ds = ((line[:, 0] - xs)**2 + (line[:, 1] - ys)**2)**0.5

    idx = ds < dist[0]
    xs[idx] = line[idx, 0] + dist[0] * np.cos(perp[idx])
    ys[idx] = line[idx, 1] + dist[0] * np.sin(perp[idx])

    idx = ds > dist[1]
    xs[idx] = line[idx, 0] + dist[1] * np.cos(perp[idx])
    ys[idx] = line[idx, 1] + dist[1] * np.sin(perp[idx])

    return xs, ys

def set_eta(x, y):
    eta = np.zeros(x.shape, dtype=float)

    # thalweg_eta = np.maximum(0.0, (y - 3313760.0))/(3367300.0 - 3313760.0) * 1.2
    # thalweg_eta = np.ones(y.shape) * 0.5

    # y0 = [0, 23388517, 3461404, 9e9]
    # eta0 = [0, 0, 0, 0]
    # eta = np.interp(y, y0, eta0)
    return eta

def set_eta_thalweg(x, y, z):
    # approximation based on bathymetry
    eta = np.zeros(x.shape, dtype=float)

    # smooth bathymetry along thalweg because the elevation is smoother than bathymetry
    mean_dl = np.mean(get_dist_increment(np.c_[x, y]))
    z_smooth = moving_average(z, n=int(max(100.0/(mean_dl+1e-6), 2)), self_weights=2)

    const_depth = 1.0
    coastal_z = [0.0, 3.0]

    # coastal (deep): assume zero
    idx = z_smooth <= coastal_z[0]
    # do nothing, use intial value 0

    # upland (high): assume constant depth
    idx = z_smooth >= coastal_z[1]
    eta[idx] = z_smooth[idx] + const_depth

    # transitional zone: assume linear transition
    idx = (z_smooth>coastal_z[0])*(z_smooth<coastal_z[1])
    eta[idx] = (z_smooth[idx]+const_depth) * (z_smooth[idx]-coastal_z[0])/(coastal_z[1]-coastal_z[0])

    return eta


def get_two_banks(S_list, thalweg, thalweg_eta, search_length, search_steps, min_width):
    range_arcs = []
    # find perpendicular direction along thalweg at each point
    perp = get_perpendicular_angle(thalweg)

    # find search area for a thalweg, consisting of two lines on each side
    xt_right = thalweg[:, 0] + search_length * np.cos(perp)
    yt_right = thalweg[:, 1] + search_length * np.sin(perp)
    xt_left = thalweg[:, 0] + search_length * np.cos(perp + np.pi)
    yt_left = thalweg[:, 1] + search_length * np.sin(perp + np.pi)

    # Diagnostic: save search area as SMS arcs
    range_arcs += [SMS_ARC(points=np.c_[xt_left, yt_left], src_prj='cpp'), SMS_ARC(points=np.c_[xt_right, yt_right], src_prj='cpp')]

    # find two banks
    x_banks_left, y_banks_left = \
        get_bank(S_list, thalweg[:, 0], thalweg[:, 1], thalweg_eta, xt_left, yt_left, search_steps)
    x_banks_right, y_banks_right = \
        get_bank(S_list, thalweg[:, 0], thalweg[:, 1], thalweg_eta, xt_right, yt_right, search_steps)

    # get attributes of the initial banks
    # average width, for deciding nudging distance
    if x_banks_left is None or x_banks_right is None:
        print('warning: failed to find banks ... ')
        return None, None, None, None, None, None

    bank2bank_width = ( (x_banks_left - x_banks_right)**2 + (y_banks_left - y_banks_right)**2 ) **0.5

    # deal with very small widths
    ismall = bank2bank_width < min_width
    x_banks_right[ismall] = thalweg[ismall, 0] + min_width/2 * np.cos(perp[ismall])
    y_banks_right[ismall] = thalweg[ismall, 1] + min_width/2 * np.sin(perp[ismall])
    x_banks_left[ismall] = thalweg[ismall, 0] + min_width/2 * np.cos(perp[ismall] + np.pi)
    y_banks_left[ismall] = thalweg[ismall, 1] + min_width/2 * np.sin(perp[ismall] + np.pi)
    bank2bank_width = ( (x_banks_left - x_banks_right)**2 + (y_banks_left - y_banks_right)**2 ) **0.5

    # SMS_MAP(arcs=range_arcs).writer(filename=f'{output_dir}/bank_range.map')

    return x_banks_left, y_banks_left, x_banks_right, y_banks_right, perp, bank2bank_width

def redistribute_arc(line, line_smooth, channel_width, this_nrow_arcs, smooth_option=1, R_coef=0.4, length_width_ratio=6.0, reso_thres=(5, 300), endpoints_scale=1.0, idryrun=False):

    cross_channel_length_scale = channel_width/this_nrow_arcs

    retained_points = np.ones((line.shape[0]), dtype=bool)

    if (line.shape[0] < 2):
        print(f'warning: line only has one point, no need for redistributing')
        return line, line_smooth, np.zeros(line.shape[0], dtype=float), retained_points

    # along-river distance, for redistribution
    dist_along_thalweg = get_dist_increment(line)

    # # smooth curvature to get rid of small-scale zig-zags
    if smooth_option == -1:  # existing line_smooth
        curv = curvature(line_smooth)
    elif smooth_option == 0:  # no smoothing
        line_smooth = line
        curv = curvature(line)
    elif smooth_option == 1:  # Option 1: moving average
        line_smooth = moving_average(line, n=30, self_weights=2)
        curv = curvature(line_smooth)
    elif smooth_option == 2:  # Option 2: spline (slow and doesn't really work because original points are preserved)
        smooth_factor = 4
        #create spline function
        f, u = interpolate.splprep([line[:, 0], line[:, 1]], s=10, per=0)
        #create interpolated lists of points
        uint = np.interp(np.arange(0, len(u)-1+1/smooth_factor, 1/smooth_factor), np.arange(0, len(u)), u)
        xint, yint = interpolate.splev(uint, f)
        line_smooth = np.c_[xint, yint]
        curv_sp = curvature(line_smooth)
        curv = curv_sp[::smooth_factor]
    '''
    plt.scatter(line_smooth[:, 0], line_smooth[:, 1])
    plt.scatter(line[:, 0], line[:, 1], s=0.3)
    plt.show()
    '''

    # use different interval along the line to calculate curvature
    if smooth_option != 2:
        for i in [1, 2]:  # more combinations -> more conservative, i.e., larger curvature
            for j in range(i):
                curv[j::i] = np.maximum(curv[j::i], curvature(line_smooth[j::i]))

    R = 1.0/(curv+1e-10)
    # resolution at points
    river_resolution = np.minimum(R_coef * R, length_width_ratio * cross_channel_length_scale)
    river_resolution = np.minimum(np.maximum(reso_thres[0], river_resolution), reso_thres[1])

    if idryrun:
        river_resolution /= 1.0  # testing

    # deprecated: increase resolution near endpoints for better intersections
    if endpoints_scale != 1.0:
        for k in [0.5, 1, 2]:
            starting_points = np.r_[0.0, np.cumsum(dist_along_thalweg)] < k * np.mean(channel_width)  # tweak within the length of k river widths
            river_resolution[starting_points] /= endpoints_scale
            ending_points = np.flip(np.r_[0.0, np.cumsum(np.flip(dist_along_thalweg))]) <  k * np.mean(channel_width)  # tweak within the length of 3 river widths 
            river_resolution[ending_points] /= endpoints_scale

    # resolution between two points
    river_resolution_seg = (river_resolution[:-1]+river_resolution[1:])/2  # resolution between two points

    idx = 0
    this_seg_length = dist_along_thalweg[0]  # dist between pt0 and pt1
    while idx < len(dist_along_thalweg)-1:
        if this_seg_length < river_resolution_seg[idx]:  # resolution of the seg between pt0 and pt1
            retained_points[idx+1] = False  # remove point
            this_seg_length += dist_along_thalweg[idx+1]
        else:
            this_seg_length = dist_along_thalweg[idx+1]
        idx += 1  # move original arc forward
    # last point should be retained
    retained_points[-1] = True

    return line[retained_points, :], line_smooth, river_resolution, retained_points

def get_thalweg_neighbors(thalwegs, thalweg_endpoints):
    # rounding_digits = 12
    # thalweg_endpoints_unique = np.unique(np.array(thalweg_endpoints).round(decimals=rounding_digits), axis=0)
    # thalweg_bombs = np.zeros((len(thalwegs), 2), dtype=bool)

    thalweg_neighbors = [None] * len(thalwegs) * 2
    for i, thalweg in enumerate(thalwegs):
        dist = thalweg_endpoints[:, :] - thalweg[0, :]
        same_points = dist[:, 0]**2 + dist[:, 1]**2 < 1e-10**2
        if sum(same_points) > 1:  # intersection point
            thalweg_neighbors[2*i] = np.argwhere(same_points)[0]

    return thalweg_neighbors

def bomb_line(line, blast_radius):
    valid_idx = np.ones(len(line[:, 0]), dtype=bool)
    valid_idx_headtail = np.ones((len(line[:, 0]), 2), dtype=bool)

    for k in [0, -1]:
        valid_idx_headtail[:, -k] = abs((line[:, 0] - line[k, 0]) + 1j* (line[:, 1] - line[k, 1])) >= blast_radius[k]
        valid_idx *= valid_idx_headtail[:, -k]

    return valid_idx, valid_idx_headtail

def width2narcs(width):
    return int(3 + np.floor(0.35*width**0.25))

def points2GeoDataFrame(points: np.ndarray, crs='epsg:4326'):
    '''
    Takes a (n, 2) array of (lon, lat) coordinates and convert to GeoPanda's GeoDataFrame
    '''
    df = pd.DataFrame({'lon':points[:, 0], 'lat':points[:, 1]})
    df['coords'] = list(zip(df['lon'],df['lat']))
    df['coords'] = df['coords'].apply(Point)
    return gpd.GeoDataFrame(df, geometry='coords', crs=crs)

def list2gdf(obj_list, crs='epsg:4326'):
    if isinstance(obj_list, list) or isinstance(obj_list, np.ndarray):
        return gpd.GeoDataFrame(index=range(len(obj_list)), crs=crs, geometry=obj_list)
    else:
        raise TypeError('Input objects must be in a list or np.array')

def clean_intersections(arcs=None, target_polygons=None, snap_points: np.ndarray=None, buffer_coef=0.2):
    '''
    Clean arcs (LineStringList, a list of Shapely's LineString objects)
    by first intersecting them (by unary_union),
    then snapping target points (within 'target_polygons') to 'snap_points',
    '''
    if isinstance(arcs, list):
        arcs_gdf = gpd.GeoDataFrame({'index': range(len(arcs)),'geometry': arcs})
    elif isinstance(arcs, gpd.GeoDataFrame):
        arcs_gdf = arcs
    else:
        raise TypeError()
        
    if isinstance(target_polygons, list):
        target_poly_gdf = list2gdf(target_polygons)
    elif isinstance(target_polygons, gpd.GeoDataFrame):
        target_poly_gdf = target_polygons
    else:
        raise TypeError()

    # ------------------snapping -------------------------------------------
    # mapping between arcs to arc points
    arcs = arcs_gdf.geometry.unary_union.geoms
    arc_points = Geoms_XY(geom_list=arcs, crs='epsg:4326')
    arc_points.snap_to_points(snap_points[:, :2], target_poly_gdf)

    snapped_arcs = np.array([arc for arc in gpd.GeoDataFrame(geometry=arc_points.geom_list).unary_union.geoms])

    # # ------------------ further cleaning intersection arcs -------------------------------------------
    cleaned_arcs_gdf = gpd.GeoDataFrame(crs='epsg:4326', index=range(len(snapped_arcs)), geometry=snapped_arcs)
    lineInPolys = gpd.tools.sjoin(cleaned_arcs_gdf, target_poly_gdf, predicate="within", how='left')
    _, idx = np.unique(lineInPolys.index, return_index=True)  # some points belong to multiple polygons
    i_cleaned_noninter_arcs =  np.array(np.isnan(lineInPolys.index_right))[idx]
    i_cleaned_intersection_arcs = np.array(~np.isnan(lineInPolys.index_right))[idx]

    # keep non intersection arcs as is
    cleaned_noninter_arcs = snapped_arcs[i_cleaned_noninter_arcs]
    cleaned_inter_arcs = snapped_arcs[i_cleaned_intersection_arcs]

    # clean intersection arcs that are too close to each other
    if np.any(i_cleaned_intersection_arcs):
        snapped_intersection_arcs_gdf = gpd.GeoDataFrame(crs='epsg:4326', index=range(len(cleaned_inter_arcs)), geometry=cleaned_inter_arcs)
        tmp_gdf = snapped_intersection_arcs_gdf.to_crs("esri:102008")
        reso_gs = gpd.GeoSeries(map(Point, snap_points[:, :2]), crs='epsg:4326').to_crs('esri:102008')
        _, idx = nearest_neighbour(np.c_[tmp_gdf.centroid.x, tmp_gdf.centroid.y], np.c_[reso_gs.x, reso_gs.y])

        arc_buffer = dl_lonlat2cpp(snap_points[idx, 2] * buffer_coef, snap_points[idx, 1])
        line_buf_gdf = snapped_intersection_arcs_gdf.to_crs('esri:102008').buffer(distance=arc_buffer)

        arc_points = Geoms_XY(geom_list=cleaned_inter_arcs, crs='epsg:4326')

        arc_pointsInPolys = gpd.tools.sjoin(
            points2GeoDataFrame(arc_points.xy, crs='epsg:4326').to_crs('esri:102008'),
            gpd.GeoDataFrame(geometry=line_buf_gdf), predicate="within", how='left'
        )
        arcs_buffer = [None] * len(cleaned_inter_arcs)
        for i in range(len(arcs_buffer)):
            arcs_buffer[i] = []
        for index, row in arc_pointsInPolys.iterrows():
            arcs_buffer[row.index_right].append(index)
        
        invalid_point_idx = np.empty((0, ), dtype=int)
        for i, pt in enumerate(arcs_buffer):
            arc_points_inbuffer = arc_points.xy[pt, :]
            dist, idx = nearest_neighbour(arc_points_inbuffer, np.array(arc_points.geom_list[i].xy).T)
            invalid = dist > 1e-10
            if np.any(invalid):
                invalid_point_idx = np.r_[invalid_point_idx, np.array(pt)[invalid]]
        invalid_point_idx = np.unique(invalid_point_idx)
        i_invalid = np.zeros((len(arc_points.xy), ), dtype=bool)
        i_invalid[invalid_point_idx] = True
        arc_points.snap_to_points(snap_points=arc_points.xy[~i_invalid, :])
        cleaned_inter_arcs = [arc for arc in gpd.GeoDataFrame(geometry=arc_points.geom_list).unary_union.geoms]

        # gdf = gpd.GeoDataFrame(geometry=cleaned_inter_arcs, crs='epsg:4326')
        # arcs_in_poly = gpd.tools.sjoin(gdf.to_crs('esri:102008'), target_poly_gdf.to_crs('esri:102008'), predicate="within", how='left')
        # arcs_in_poly = arcs_in_poly.to_crs('epsg:4326')
        # polys = np.unique(arcs_in_poly.index_right)
        # hull_list = [arcs_in_poly[arcs_in_poly.index_right==poly].unary_union.convex_hull.boundary for poly in polys]
        # cleaned_inter_arcs = np.array(cleaned_inter_arcs + hull_list)

    # assemble
    cleaned2_arcs = np.r_[cleaned_noninter_arcs, cleaned_inter_arcs]
    return [arc for arc in gpd.GeoDataFrame(geometry=cleaned2_arcs).unary_union.geoms]

def clean_river_arcs(river_arcs=None, total_arcs=None):

    total_arcs_geomsxy = Geoms_XY(total_arcs, crs="epsg:4326")

    river_arcs_cleaned = deepcopy(river_arcs)
    for i, river in enumerate(river_arcs_cleaned):
        for j, arc in enumerate(river):
            if arc is not None:
                _, idx = nearest_neighbour(arc.points[:, :2], total_arcs_geomsxy.xy)
                river_arcs_cleaned[i, j].points[:, :2] = total_arcs_geomsxy.xy[idx, :]

    return river_arcs_cleaned

def generate_river_outline_polys(river_arcs=None):
    river_outline_polygons = []
    river_polygons = [None] * river_arcs.shape[0]
    for i, river in enumerate(river_arcs):
        # save river polygon (enclosed by two out-most arcs and two cross-river transects at both ends)
        if sum(river[:] != None) >= 2:  # at least two rows of arcs to make a polygon
            river_polygons[i] = []
            idx = np.argwhere(river_arcs[i, :] != None).squeeze()
            valid_river_arcs = river_arcs[i, idx]
            for j in range(1):  # range(len(valid_river_arcs)-1):
                mls_uu = unary_union(LineString(np.r_[valid_river_arcs[0].points[:, :2], np.flipud(valid_river_arcs[-1].points[:, :2]), valid_river_arcs[0].points[0, :2].reshape(-1,2)]))
                for polygon in polygonize(mls_uu):
                    river_polygons[i].append(polygon)

    # convert to 1D list for later conversion to shapefile
    for river_polygon in river_polygons:
        if river_polygon is not None:
            river_outline_polygons += river_polygon
    
    return river_outline_polygons

def clean_arcs(LineStringList, blast_center, blast_radius, minimum_arc_length=10):
    uu = gpd.GeoDataFrame({'index': range(len(LineStringList)),'geometry': LineStringList}).geometry.unary_union

    # # sample
    # uu = GeoSeries(uu).set_crs('epsg:4326')
    # uu_meter = uu.to_crs('esri:102008')

    # ls1 = LineString([(0,0), (1,1)])
    # ls2 = LineString([(0,1), (1,0)])
    # inter = ls1.intersection(ls2)
    # gc = split_line_by_point(ls1, MultiPoint([inter, Point(0.25, 0.25)]))
    
    uu_xy = np.empty((0, 2), dtype=float)
    uu_xy_idx = np.empty((0, 2), dtype=int)
    idx = 0
    for i, line in enumerate(uu.geoms):
        uu_xy = np.r_[uu_xy, np.array(line.xy).T]
        uu_xy_idx = np.r_[uu_xy_idx, np.c_[idx, idx+len(line.xy[0])]]
        idx += len(line.xy[0])
    
    i_snap = np.zeros((len(uu_xy), ), dtype=bool)
    for center, radius in zip(blast_center, blast_radius):
        idx = abs(uu_xy[:, 0]+1j*uu_xy[:, 1] - (center[0]+1j*center[1])) < radius
        i_snap[idx] = True
    
    _, idx = nearest_neighbour(uu_xy[i_snap], np.c_[blast_center[:, 0], blast_center[:, 1]])
    uu_xy[i_snap, :] = blast_center[idx, :]

    LineStringList_Cleaned = []
    for i, line in enumerate(uu.geoms):
        line_xy = uu_xy[uu_xy_idx[i, 0]:uu_xy_idx[i, 1], :]
        LineStringList_Cleaned.append(LineString(np.c_[line_xy[:, 0], line_xy[:, 1]]))

    river_arcs_gdf = gpd.GeoDataFrame({'index':range(len(LineStringList_Cleaned)),'geometry':LineStringList_Cleaned})
    uu = river_arcs_gdf.geometry.unary_union

    cleaned_arc_list = []
    for i, line in enumerate(uu.geoms):
        line = LineString(np.array(line.xy).T)
        if line.length > dl_cpp2lonlat(minimum_arc_length, line.xy[1][0]):
            cleaned_arc_list.append(line)

    return cleaned_arc_list

def improve_thalwegs(S_list, dl, line, search_length, perp, mpi_print_prefix):
    x = line[:, 0]
    y = line[:, 1]

    xt_right = line[:, 0] + search_length * np.cos(perp)
    yt_right = line[:, 1] + search_length * np.sin(perp)
    xt_left = line[:, 0] + search_length * np.cos(perp + np.pi)
    yt_left = line[:, 1] + search_length * np.sin(perp + np.pi)

    __search_steps = int(np.max(search_length/dl))
    __search_steps = max(5, __search_steps)  # give it at least something to search for, i.e., the length of 5 grid points

    x_new = np.empty((len(x), 2), dtype=float); x_new.fill(np.nan)
    y_new = np.empty((len(x), 2), dtype=float); y_new.fill(np.nan)
    elev_new = np.empty((len(x), 2), dtype=float); y_new.fill(np.nan)
    thalweg_idx = np.ones(xt_right.shape) * 9999
    for k, [xt, yt] in enumerate([[xt_left, yt_left], [xt_right, yt_right]]):
        xts = np.linspace(x, xt, __search_steps, axis=1)
        yts = np.linspace(y, yt, __search_steps, axis=1)

        ''' One tile
        [jj, ii], valid = Sidx(S, xts[:], yts[:])
        if valid.all():
            elevs = S.elev[ii, jj]
            low = np.argpartition(elevs, min(10, elevs.shape[1]-1), axis=1)
            thalweg_idx = np.median(low[:, :10], axis=1).astype(int)

            x_new[:, k] = xts[range(len(x)), thalweg_idx]
            y_new[:, k] = yts[range(len(x)), thalweg_idx]
            elev_new[:, k] = elevs[range(len(x)), thalweg_idx]
        else:
            return np.c_[x, y], False  # return orignial thalweg
        '''
        # multiple tiles in a tile list
        elevs = get_elev_from_tiles(xts, yts, S_list)
        if np.isnan(elevs).any():
            print(f'{mpi_print_prefix} Warning: nan found in elevs\n' + \
                  f'when trying to improve Thalweg: {np.c_[x, y]}')
            return np.c_[x, y], False  # return orignial thalweg

        if elevs is not None:
            if elevs.shape[1] < 2:
                print(f'{mpi_print_prefix} Warning: elevs shape[1] < 2')
                return np.c_[x, y], False  # return orignial thalweg

            n_low = min(10, elevs.shape[1]-1)
            low = np.argpartition(elevs, n_low, axis=1)
            thalweg_idx = np.median(low[:, :n_low], axis=1).astype(int)

            if any(thalweg_idx<0) or any(thalweg_idx>=len(x)):
                return np.c_[x, y], False  # return orignial thalweg

            x_new[:, k] = xts[range(len(x)), thalweg_idx]
            y_new[:, k] = yts[range(len(x)), thalweg_idx]
            elev_new[:, k] = elevs[range(len(x)), thalweg_idx]
        else:
            return np.c_[x, y], False  # return orignial thalweg

    left_or_right = elev_new[:, 0] > elev_new[:, 1]
    x_real = x_new[range(len(x)), left_or_right.astype(int)]  # if left higher than right, then use right
    y_real = y_new[range(len(x)), left_or_right.astype(int)]

    return np.c_[x_real, y_real], True

def get_bank(S_list, x, y, thalweg_eta, xt, yt, search_steps=100, search_tolerance=5):
    '''
    Get a bank on one side of the thalweg (x, y)
    Inputs:
        x, y, eta along a thalweg
        parameter deciding the search area: search_stps
    '''

    # search_steps_tile = np.repeat(np.arange(search_steps).reshape(1, -1), len(x), axis=0)  # expanded to the search area

    # form a search area between thalweg and search limit
    xts = np.linspace(x, xt, search_steps, axis=1)
    yts = np.linspace(y, yt, search_steps, axis=1)

    eta_stream = np.tile(thalweg_eta.reshape(-1, 1), (1, search_steps))  # expanded to the search area

    elevs = get_elev_from_tiles(xts, yts, S_list)
    if elevs is None:
        return None, None

    R = (elevs - eta_stream)  # closeness to target depth
    bank_idx = np.argmax(R>0, axis=1)

    invalid = bank_idx == 0
    bank_idx[invalid] = np.argmin(abs(R[invalid, :]), axis=1)

    # R_sort_idx = np.argsort(R)
    # bank_idx = np.min(R_sort_idx[:, :min(search_steps, search_tolerance)], axis=1)

    # x_banks = S.lon[jj[range(0, len(x)), bank_idx]]
    # y_banks = S.lat[ii[range(0, len(x)), bank_idx]]
    x_banks = xts[range(len(x)), bank_idx]
    y_banks = yts[range(len(x)), bank_idx]

    return x_banks, y_banks


def make_river_map(
    tif_fnames: list, thalweg_shp_fname = '',
    selected_thalweg = None, output_dir = './', output_prefix = '', mpi_print_prefix = '',
    MapUnit2METER = 1, river_threshold = (5, 400), 
    outer_arcs_positions = (), length_width_ratio = 12.0,
    i_close_poly = True, i_blast_intersection = False,
    blast_radius_scale = 0.7, bomb_radius_coef = 0.2,
    i_DEM_cache = True
):
    '''
    [Core routine for making river maps]

    tif_fnames: a list of TIF file names. These TIFs should cover the area of interest and be arranged by priority (higher priority ones in front)

    thalweg_shp_fname: name of a polyline shapefile containing the thalwegs

    selected_thalweg: indices of selected thalwegs for which the river arcs will be sought.

    output_dir: must specify one.

    output_prefix: a prefix of the output files, mainly used by the caller of this script; can be empty

    mpi_print_prefix: a prefix string to identify the calling mpi processe in the output messages; can be empty

    MapUnit2METER = 1:  to be replaced by projection code, e.g., epsg: 4326, esri: 120008, etc.

    river_threshold:  minimum and maximum river widths (in meters) to be resolved

    outer_arc_positions: relative position of outer arcs, e.g., (0.1, 0.2) will add 2 outer arcs on each side of the river (4 in total),
        at 0.1*riverwidth and 0.2*riverwidth from the banks.

    length_width_ratio: a ratio of element length in the along-channel direction to river width;
                        when a river is narrower than the lower limit, the bank will be nudged (see next parameter) to widen the river
    
    i_close_poly: whether to add cross-channel arcs to enclose river arcs into a polygon

    i_blast_intersection: whether to replace intersecting arcs (often noisy) at river intersections with scatter points (cleaner)
    
    blast_radius_scale:  coef controlling the blast radius at intersections, a larger number leads to more intersection features being deleted

    bomb_radius_coef:  coef controlling the spacing among intersection joints, a larger number leads to sparser intersection joints

    i_DEM_cache : Whether or not to read DEM info from cache.
                  Reading from original *.tif files can be slow, so the default option is True
    '''

    # ------------------------- other input parameters not exposed to user ---------------------------
    iAdditionalOutputs = False
    nudge_ratio = np.array((0.3, 2.0))  # ratio between nudging distance to mean half-channel-width
    thalweg_smooth_shp_fname = None  # deprecated: name of a polyline shapefile containing the smoothed thalwegs (e.g., pre-processed by GIS tools or SMS)
    # ------------------------- end other inputs ---------------------------

    # ----------------------   pre-process some inputs -------------------------
    river_threshold = np.array(river_threshold) / MapUnit2METER

    outer_arcs_positions = np.array(outer_arcs_positions).reshape(-1, )  # example: np.array([-0.1, 0.1])
    if np.any(outer_arcs_positions <= 0.0):
        raise ValueError('outer arcs position must > 0, a pair of arcs (one on each side of the river) will be placed for each position value')

    # maximum number of arcs to resolve a channel (including bank arcs, inner arcs and outer arcs)
    max_nrow_arcs = width2narcs(river_threshold[-1]) + 2 * outer_arcs_positions.size
    # refine toward endpoints to better fit intersections
    endpoints_scale = 1.1
    # ---------------------- end pre-processing some inputs -------------------------

    # ------------------------- read DEM ---------------------------
    main_dem_id = 1
    S_list = []

    nvalid_tile = 0
    for i, tif_fname in enumerate(tif_fnames):
        if tif_fname is None:
            continue
        else:
            nvalid_tile += 1

        S, is_new_cache = Tif2XYZ(tif_fname=tif_fname, cache=i_DEM_cache)

        print(f'{mpi_print_prefix} [{os.path.basename(tif_fname)}] DEM box: {min(S.lon)}, {min(S.lat)}, {max(S.lon)}, {max(S.lat)}')
        S_list.append(S)

        if nvalid_tile == main_dem_id:
            S_x, S_y = lonlat2cpp(S.lon[:2], S.lat[:2])
            dx = S_x[1] - S_x[0]
            dy = S_y[1] - S_y[0]
            dl = (abs(dx) + abs(dy)) / 2
            search_length = river_threshold[-1] * 1.1
            search_steps = int(river_threshold[-1] / dl)

    if nvalid_tile == 0:
        raise ValueError('Fatal Error: no valid DEM tiles')

    # ------------------------- read thalweg ---------------------------
    xyz, l2g, curv, _ = get_all_points_from_shp(thalweg_shp_fname)
    xyz[:, 0], xyz[:, 1] = lonlat2cpp(xyz[:, 0], xyz[:, 1])

    # Optional (deprecated): provide a smoothed thalweg (on the 2D plane, not smoothed in z) to guide vertices distribution.
    # The default option is to let the script do the smoothing
    if thalweg_smooth_shp_fname is not None:
        xyz_s, l2g_s, curv_s, _ = get_all_points_from_shp(thalweg_smooth_shp_fname)
        xyz_s[:, 0], xyz_s[:, 1] = lonlat2cpp(xyz_s[:, 0], xyz_s[:, 1])


    thalwegs = []
    thalwegs_smooth = []
    thalwegs_curv = []
    thalweg_endpoints = np.empty((0, 2), dtype=float)
    if selected_thalweg is None:
        selected_thalweg = np.arange(len(l2g))
    for i, idx in enumerate(l2g):
        if i in selected_thalweg:
            # print(f'Arc {i+1} of {len(thalwegs)}')
            thalwegs.append(xyz[idx, :])
            thalwegs_curv.append(curv[idx])
            thalweg_endpoints = np.r_[thalweg_endpoints, np.reshape(xyz[idx][0, :], (1, 2))]
            thalweg_endpoints = np.r_[thalweg_endpoints, np.reshape(xyz[idx][-1, :], (1, 2))]
            if thalweg_smooth_shp_fname is not None:
                thalwegs_smooth.append(xyz_s[idx, :])
            else:
                thalwegs_smooth.append(None)

    # ------------------------- Dry run (finding approximate locations of two banks) ---------------------------
    print(f'{mpi_print_prefix} Dry run')

    thalweg_endpoints_width = np.empty((len(thalwegs)*2, 1), dtype=float); thalweg_endpoints_width.fill(np.nan)
    thalwegs_neighbors = np.empty((len(thalwegs), 2), dtype=object)  # [, 0] is head, [, 1] is tail
    thalweg_widths = [None] * len(thalwegs)
    valid_thalwegs = [True] * len(thalwegs)
    original_banks = [None] * len(thalwegs) * 2

    for i, [thalweg, curv] in enumerate(zip(thalwegs, thalwegs_curv)):
        # print(f'Dry run: Arc {i+1} of {len(thalwegs)}')

        elevs = get_elev_from_tiles(thalweg[:, 0], thalweg[:, 1], S_list)
        if elevs is None:
            print(f"{mpi_print_prefix} warning: some elevs not found on thalweg {i+1}, the thalweg will be neglected ...")
            valid_thalwegs[i] = False
            continue

        # set water level at each point along the thalweg, based on observation, simulation, estimation, etc.
        thalweg_eta = set_eta_thalweg(thalweg[:, 0], thalweg[:, 1], elevs)

        x_banks_left, y_banks_left, x_banks_right, y_banks_right, _, width = \
            get_two_banks(S_list, thalweg, thalweg_eta, search_length, search_steps, min_width=river_threshold[0])
        thalweg_widths[i] = width
        if width is None:
            thalweg_endpoints_width[i*2] = 0.0
            thalweg_endpoints_width[i*2+1] = 0.0
        else:
            thalweg_endpoints_width[i*2] = width[0]
            thalweg_endpoints_width[i*2+1] = width[-1]

        if len(thalweg[:, 0]) < 2:
            print(f"{mpi_print_prefix} warning: thalweg {i+1} only has one point, neglecting ...")
            valid_thalwegs[i] = False
            continue

        if x_banks_left is None or x_banks_right is None:
            print(f"{mpi_print_prefix} warning: thalweg {i+1} out of DEM coverage, neglecting ...")
            valid_thalwegs[i] = False
            continue

        original_banks[2*i] = SMS_ARC(points=np.c_[x_banks_left, y_banks_left], src_prj='cpp')
        original_banks[2*i+1] = SMS_ARC(points=np.c_[x_banks_right, y_banks_right], src_prj='cpp')

    # End Dry run: found valid river segments; record approximate channel width

    # ------------------------- Wet run ---------------------------
    # initialize some lists and arrays
    bank_arcs = np.empty((len(thalwegs), 2), dtype=object) # left bank and right bank for each thalweg
    bank_arcs_raw = np.empty((len(thalwegs), 2), dtype=object)
    bank_arcs_final = np.empty((len(thalwegs), 2), dtype=object)
    cc_arcs = np.empty((len(thalwegs), 2), dtype=object)  # [, 0] is head, [, 1] is tail
    river_arcs = np.empty((len(thalwegs), max_nrow_arcs), dtype=object)
    blast_radius = -np.ones((len(thalwegs), 2), dtype=float)
    blast_center = np.zeros((len(thalwegs), 2), dtype=complex)
    smoothed_thalwegs = [None] * len(thalwegs)
    redistributed_thalwegs = [None] * len(thalwegs)
    corrected_thalwegs = [None] * len(thalwegs)
    centerlines = [None] * len(thalwegs)
    thalwegs_cc_reso = [None] * len(thalwegs)
    final_thalwegs = [None] * len(thalwegs)
    real_bank_width = np.zeros((len(thalwegs), 2), dtype=float)  # [, 0] is head, [, 1] is tail
    bombed_points = np.empty((0, 3), dtype=float)  # left bank and right bank for each thalweg
    bombs = [None]* len(thalwegs) *2 # left bank and right bank for each thalweg
    bomb_polygons = []

    # enumerate each thalweg
    for i, [thalweg, curv, width, valid_thalweg, thalweg_smooth] in enumerate(zip(thalwegs, thalwegs_curv, thalweg_widths, valid_thalwegs, thalwegs_smooth)):
        # print(f'{mpi_print_prefix} Wet run: Arc {i+1} of {len(thalwegs)}')

        if not valid_thalweg:
            print(f"{mpi_print_prefix} Thalweg {i} marked as invalid in dry run, skipping ...")
            continue

        # Redistribute thalwegs vertices
        this_nrow_arcs = min(max_nrow_arcs, width2narcs(np.mean(width)))
        thalweg, thalweg_smooth, reso, retained_idx = redistribute_arc(thalweg, thalweg_smooth, width, this_nrow_arcs, length_width_ratio=length_width_ratio, smooth_option=1, endpoints_scale=endpoints_scale, idryrun=True)
        smoothed_thalwegs[i] = SMS_ARC(points=np.c_[thalweg_smooth[:, 0], thalweg_smooth[:, 1]], src_prj='cpp')
        redistributed_thalwegs[i] = SMS_ARC(points=np.c_[thalweg[:, 0], thalweg[:, 1]], src_prj='cpp')

        if len(thalweg[:, 0]) < 2:
            print(f"{mpi_print_prefix} warning: thalweg {i+1} only has one point after redistribution, neglecting ...")
            continue

        # update thalweg info
        elevs = get_elev_from_tiles(thalweg[:, 0],thalweg[:, 1], S_list)
        thalweg_eta = set_eta_thalweg(thalweg[:, 0], thalweg[:, 1], elevs)

        # re-make banks based on redistributed thalweg
        x_banks_left, y_banks_left, x_banks_right, y_banks_right, perp, width = \
            get_two_banks(S_list, thalweg, thalweg_eta, search_length, search_steps, min_width=river_threshold[0])

        # correct thalwegs
        width_moving_avg = moving_average(width, n=10)
        thalweg, is_corrected= improve_thalwegs(S_list, dl, thalweg, width_moving_avg*0.5, perp, mpi_print_prefix)
        if not is_corrected:
            print(f"{mpi_print_prefix} warning: thalweg {i+1} (head: {thalweg[0]}; tail: {thalweg[-1]}) failed to correct, using original thalweg ...")
        corrected_thalwegs[i] = SMS_ARC(points=np.c_[thalweg[:, 0], thalweg[:, 1]], src_prj='cpp')

        # Redistribute thalwegs vertices
        this_nrow_arcs = min(max_nrow_arcs, width2narcs(np.mean(width)))
        thalweg, thalweg_smooth, reso, retained_idx = redistribute_arc(thalweg, thalweg_smooth[retained_idx], width, this_nrow_arcs, length_width_ratio=length_width_ratio, smooth_option=1, endpoints_scale=endpoints_scale)
        smoothed_thalwegs[i] = SMS_ARC(points=np.c_[thalweg_smooth[:, 0], thalweg_smooth[:, 1]], src_prj='cpp')
        redistributed_thalwegs[i] = SMS_ARC(points=np.c_[thalweg[:, 0], thalweg[:, 1]], src_prj='cpp')

        # Smooth thalweg
        thalweg, perp = smooth_thalweg(thalweg, ang_diff_shres=np.pi/2.4)
        final_thalwegs[i] = SMS_ARC(points=np.c_[thalweg[:, 0], thalweg[:, 1]], src_prj='cpp')

        # update thalweg info
        elevs = get_elev_from_tiles(thalweg[:, 0],thalweg[:, 1], S_list)
        thalweg_eta = set_eta_thalweg(thalweg[:, 0], thalweg[:, 1], elevs)

        # re-make banks based on corrected thalweg
        x_banks_left, y_banks_left, x_banks_right, y_banks_right, perp, width = \
            get_two_banks(S_list, thalweg, thalweg_eta, search_length, search_steps, min_width=river_threshold[0])

        # touch-ups on the two banks
        if x_banks_left is None or x_banks_right is None:
            print(f'{mpi_print_prefix} warning: cannot find banks for thalweg {i+1} after redistribution, neglecting ... ')
            continue

        # nudge banks
        x_banks_left, y_banks_left = nudge_bank(thalweg, perp+np.pi, x_banks_left, y_banks_left, dist=nudge_ratio*0.5*np.mean(width))
        x_banks_right, y_banks_right = nudge_bank(thalweg, perp, x_banks_right, y_banks_right, dist=nudge_ratio*0.5*np.mean(width))

        # smooth banks
        thalweg, x_banks_left, y_banks_left, x_banks_right, y_banks_right, perp = smooth_bank(thalweg, x_banks_left, y_banks_left, x_banks_right, y_banks_right)
        if thalweg is None:
            continue
        thalweg, x_banks_right, y_banks_right, x_banks_left, y_banks_left, perp = smooth_bank(thalweg, x_banks_right, y_banks_right, x_banks_left, y_banks_left)
        if thalweg is None:
            continue

        # update width
        width = ((x_banks_left-x_banks_right)**2 + (y_banks_left-y_banks_right)**2)**0.5

        # get actual resolution along redistributed/smoothed thalweg
        # thalweg_resolutions[i] = np.c_[(thalweg[:-1, :]+thalweg[1:, :])/2, get_dist_increment(thalweg)]

        # make inner arcs between two banks
        this_nrow_arcs = min(max_nrow_arcs, width2narcs(np.mean(width)))
        thalwegs_cc_reso[i] = np.c_[thalweg, width / this_nrow_arcs]  # record cross-channel (cc) resolution at each thalweg point
        arc_position = np.linspace(0.0, 1.0, this_nrow_arcs)
        arc_position = np.r_[-outer_arcs_positions, arc_position, 1.0+outer_arcs_positions].reshape(-1, 1)
        x_river_arcs = x_banks_left.reshape(1, -1) + np.matmul(arc_position, (x_banks_right-x_banks_left).reshape(1, -1))
        y_river_arcs = y_banks_left.reshape(1, -1) + np.matmul(arc_position, (y_banks_right-y_banks_left).reshape(1, -1))

        z_centerline = this_nrow_arcs + width/1e4

        # determine blast radius based on mean channel width at an intersection
        valid_points = np.ones(x_banks_left.shape, dtype=bool)
        for k in [0, -1]:  # head and tail
            dist = thalweg_endpoints[:, :] - thalweg[k, :]
            # keep the index of self to facilitate the calculation of the blast_radius, which is based on the widths of all intersecting rivers
            thalwegs_neighbors[i, k] = np.argwhere(dist[:, 0]**2 + dist[:, 1]**2 < 200**2)
            if len(thalwegs_neighbors[i, k]) > 1:  # at least two rivers, i.e., at least one neighbor (except self)
                blast_radius[i, k] = blast_radius_scale * np.mean(thalweg_endpoints_width[thalwegs_neighbors[i, k]])
                blast_center[i, k] = (x_banks_left[k]+x_banks_right[k])/2 + 1j*(y_banks_left[k]+y_banks_right[k])/2

        # bomb intersections
        valid_l, valid_l_headtail = bomb_line(np.c_[x_banks_left, y_banks_left], blast_radius[i, :])
        valid_r, valid_r_headtail = bomb_line(np.c_[x_banks_right, y_banks_right], blast_radius[i, :])
        valid_points = valid_l * valid_r
        valid_points_headtail = valid_l_headtail * valid_r_headtail

        bombed_idx = ~valid_points
        if not i_blast_intersection:
            valid_points[:] = True

        # assemble banks
        for k, line in enumerate([np.c_[x_banks_left, y_banks_left], np.c_[x_banks_right, y_banks_right]]):
            bank_arcs_raw[i, k] = SMS_ARC(points=np.c_[line[:, 0], line[:, 1]], src_prj='cpp')
            if sum(valid_points) > 0:
                bank_arcs[i, k] = SMS_ARC(points=np.c_[line[valid_points, 0], line[valid_points, 1]], src_prj='cpp')

        if sum(valid_points) > 0:
            for j in [0, -1]:
                real_bank_width[i, j] = ((x_banks_left[valid_points][j]-x_banks_right[valid_points][j])**2 + (y_banks_left[valid_points][j]-y_banks_right[valid_points][j])**2)**0.5

        # quality check river arcs
        if river_quality(x_river_arcs, y_river_arcs, valid_points):
            # save centerlines
            # assemble inner arcs
            bombs_xyz = [np.empty((0,3), dtype=float)] * 2
            for k, [x_river_arc, y_river_arc] in enumerate(zip(x_river_arcs, y_river_arcs)):
                line = np.c_[x_river_arc, y_river_arc]
                if sum(valid_points) > 0:
                    # snap vertices too close to each other
                    line[valid_points, :] = snap_vertices(line[valid_points, :], width[valid_points] * 0.3)  # optional: thalweg_resolution*0.75

                    # ----------Save-------
                    # save final bank arcs
                    if k == 0:  # left bank
                        bank_arcs_final[i, 0] = SMS_ARC(points=np.c_[line[valid_points, 0], line[valid_points, 1], z_centerline[valid_points]], src_prj='cpp')
                    elif k == len(x_river_arcs)-1:  # right bank
                        bank_arcs_final[i, 1] = SMS_ARC(points=np.c_[line[valid_points, 0], line[valid_points, 1], z_centerline[valid_points]], src_prj='cpp')
                    # save inner arcs
                    river_arcs[i, k] = SMS_ARC(points=np.c_[line[valid_points, 0], line[valid_points, 1], z_centerline[valid_points]], src_prj='cpp')
                    # save centerline
                    if k == int(len(x_river_arcs)/2):
                        centerlines[i] = SMS_ARC(points=np.c_[line[valid_points, 0], line[valid_points, 1], z_centerline[valid_points]], src_prj='cpp')
                    # save bombed points
                    bombed_points = np.r_[bombed_points, np.c_[line[bombed_idx, 0], line[bombed_idx, 1], width[bombed_idx]/this_nrow_arcs]]

                # test bombs
                for l in [0, 1]:
                    if sum(~valid_points_headtail[:, l]) > 0:
                        bombs_xyz[l] = np.r_[bombs_xyz[l], np.c_[line[~valid_points_headtail[:, l]][:, :2], width[~valid_points_headtail[:, l]]/this_nrow_arcs]]

            for l in [0, 1]:
                if len(bombs_xyz[l]) > 0:
                    lon, lat = cpp2lonlat(bombs_xyz[l][:, 0], bombs_xyz[l][:, 1])
                    dl_lonlat = dl_cpp2lonlat(bombs_xyz[l][:, 2], lat0=lat)
                    # dl_cpp = dl_lonlat2cpp(dl_lonlat, lat0=lat)
                    bombs[2*i+l] = Bombs(xyz=np.c_[lon, lat, dl_lonlat], crs='epsg:4326')
            
            # save bombing polygons
            if np.all(bombed_idx):
                poly_xy = np.r_[
                    np.c_[x_river_arcs[0], y_river_arcs[0]],
                    np.flipud(np.c_[x_river_arcs[-1], y_river_arcs[-1]]),
                    np.c_[x_river_arcs[0][0], y_river_arcs[0][0]],
                ]
                poly_lon, poly_lat = cpp2lonlat(poly_xy[:, 0], poly_xy[:, 1])
                mls_uu = unary_union(LineString(np.c_[poly_lon, poly_lat]))
                for polygon in polygonize(mls_uu):
                    bomb_polygons.append(polygon)

            elif np.any(bombed_idx):
                x_river_midpoints = 0.5 * (x_river_arcs[:, :-1] + x_river_arcs[:, 1:])
                y_river_midpoints = 0.5 * (y_river_arcs[:, :-1] + y_river_arcs[:, 1:])
                for l in range(2):
                    poly_xy = None
                    if len(thalwegs_neighbors[i, l]) <= 1:  # one river, i.e., no intersections
                        continue
                    if l == 0:
                        bomb_start = 0
                        try:
                            bomb_end = np.where(~bombed_idx)[0][0]
                        except IndexError:
                            bomb_end = len(bombed_idx)
                        if bomb_end >= bomb_start:
                            poly_xy = np.r_[
                                np.c_[x_river_arcs[0][bomb_start:bomb_end], y_river_arcs[0][bomb_start:bomb_end]],
                                np.c_[x_river_midpoints[0][bomb_end-1], y_river_midpoints[0][bomb_end-1]],
                                np.c_[x_river_midpoints[-1][bomb_end-1], y_river_midpoints[-1][bomb_end-1]],
                                np.flipud(np.c_[x_river_arcs[-1][bomb_start:bomb_end], y_river_arcs[-1][bomb_start:bomb_end]]),
                                np.c_[x_river_arcs[0][bomb_start], y_river_arcs[0][bomb_start]],
                            ]
                    else:
                        try:
                            bomb_start = np.where(~bombed_idx)[0][-1]
                        except IndexError:
                            bomb_start = len(bombed_idx)
                        bomb_end = len(bombed_idx) - 1
                        if bomb_end >= bomb_start:
                            poly_xy = np.r_[
                                np.c_[x_river_midpoints[0][bomb_start], y_river_midpoints[0][bomb_start]],
                                np.c_[x_river_arcs[0][bomb_start+1:bomb_end], y_river_arcs[0][bomb_start+1:bomb_end]],
                                np.flipud(np.c_[x_river_arcs[-1][bomb_start+1:bomb_end], y_river_arcs[-1][bomb_start+1:bomb_end]]),
                                np.c_[x_river_midpoints[-1][bomb_start], y_river_midpoints[-1][bomb_start]],
                                np.c_[x_river_midpoints[0][bomb_start], y_river_midpoints[0][bomb_start]],
                            ]
                    
                    poly_lon, poly_lat = cpp2lonlat(poly_xy[:, 0], poly_xy[:, 1])
                    mls_uu = unary_union(LineString(np.c_[poly_lon, poly_lat]))
                    for polygon in polygonize(mls_uu):
                        bomb_polygons.append(polygon)

            if sum(valid_points) > 0:
                # assemble cross-channel arcs
                for j in [0, -1]:
                    cc_arcs[i, j] = SMS_ARC(points=np.c_[x_river_arcs[:, valid_points][:, j], y_river_arcs[:, valid_points][:, j]], src_prj='cpp')

    # end enumerating each thalweg

    # -----------------------------------diagnostic outputs ----------------------------
    if any(bank_arcs.flatten()):  # not all arcs are None
        SMS_MAP(arcs=bank_arcs.reshape((-1, 1))).writer(filename=f'{output_dir}/{output_prefix}bank.map')
        SMS_MAP(arcs=bank_arcs_raw.reshape((-1, 1))).writer(filename=f'{output_dir}/{output_prefix}bank_raw.map')
        SMS_MAP(arcs=cc_arcs.reshape((-1, 1))).writer(filename=f'{output_dir}/{output_prefix}cc_arcs.map')
        SMS_MAP(arcs=river_arcs.reshape((-1, 1))).writer(filename=f'{output_dir}/{output_prefix}river_arcs.map')
        SMS_MAP(detached_nodes=bombed_points).writer(filename=f'{output_dir}/{output_prefix}relax_points.map')
        SMS_MAP(arcs=smoothed_thalwegs).writer(filename=f'{output_dir}/{output_prefix}smoothed_thalweg.map')
        SMS_MAP(arcs=redistributed_thalwegs).writer(filename=f'{output_dir}/{output_prefix}redist_thalweg.map')
        SMS_MAP(arcs=corrected_thalwegs).writer(filename=f'{output_dir}/{output_prefix}corrected_thalweg.map')
    else:
        print(f'{mpi_print_prefix} No arcs found, aborted writing to *.map')
    
    del smoothed_thalwegs[:], redistributed_thalwegs[:], corrected_thalwegs[:]
    del bank_arcs, bank_arcs_raw, bombed_points, smoothed_thalwegs, redistributed_thalwegs, corrected_thalwegs

    # ------------------------- Clean up and finalize ---------------------------
    for i, thalweg_neighbors in enumerate(thalwegs_neighbors):
        for j, neibs in enumerate(thalweg_neighbors):  # head and tail
            if bombs[2*i+j] is not None:
                for nei in neibs:
                    nei = int(nei)
                    if (2*i+j) != nei:
                        bombs[2*i+j] += bombs[nei]
                        bombs[nei] = None

    # gather bomb info
    bombed_xyz = np.empty((0,3), dtype=float)
    bomb_polygons = []
    for i, bomb in enumerate(bombs):
        if bomb is not None:
            bomb_polygons.append(bomb.get_convex_hull())  # make convex hull before cleaning to include more space
            bomb.clean(bomb_radius_coef=bomb_radius_coef)
            xyz = np.c_[bomb.points.real, bomb.points.imag, bomb.res]
            bombed_xyz = np.r_[bombed_xyz, xyz]

    # Clean river intersections
    total_arcs = []
    if i_close_poly:
        arc_groups = [river_arcs, cc_arcs]
    else:
        arc_groups = [river_arcs]
    total_arcs = [LineString(line.points[:, :2]) for arcs in arc_groups for arc in arcs for line in arc if line is not None]

    if len(total_arcs) > 0:
        total_arcs_cleaned = clean_intersections(arcs=total_arcs, target_polygons=bomb_polygons, snap_points=bombed_xyz)
    else:
        total_arcs_cleaned = []
    del total_arcs[:]; del total_arcs, cc_arcs

    # ------------------------- main outputs ---------------------------
    if len(total_arcs_cleaned) > 0:
        SMS_MAP(arcs=centerlines).writer(filename=f'{output_dir}/{output_prefix}centerlines.map')
        del centerlines[:]; del centerlines
        SMS_MAP(arcs=final_thalwegs).writer(filename=f'{output_dir}/{output_prefix}final_thalweg.map')
        del final_thalwegs[:]; del final_thalwegs
        SMS_MAP(arcs=bank_arcs_final.reshape((-1, 1))).writer(filename=f'{output_dir}/{output_prefix}bank_final.map')
        del bank_arcs_final

        if len(bomb_polygons) > 0:
            gpd.GeoDataFrame(index=range(len(bomb_polygons)), crs='epsg:4326', geometry=bomb_polygons).to_file(filename=f'{output_dir}/{output_prefix}bomb_polygons.shp', driver="ESRI Shapefile")
        else:
            print(f'{mpi_print_prefix} Warning: bomb_polygons empty')
        del bomb_polygons[:]; del bomb_polygons

        if len(total_arcs_cleaned) > 0:
            SMS_MAP(arcs=geos2SmsArcList(geoms=total_arcs_cleaned)).writer(filename=f'{output_dir}/{output_prefix}total_arcs.map')
        else:
            print(f'{mpi_print_prefix} Warning: total_sms_arcs_cleaned empty')

        SMS_MAP(detached_nodes=bombed_xyz).writer(f'{output_dir}/{output_prefix}total_intersection_joints.map')

        if iAdditionalOutputs:
            total_arcs_cleaned_polys = [poly for poly in polygonize(gpd.GeoSeries(total_arcs_cleaned))]
            if len(total_arcs_cleaned_polys) > 0:
                gpd.GeoDataFrame(index=range(len(total_arcs_cleaned_polys)), crs='epsg:4326', geometry=total_arcs_cleaned_polys).to_file(filename=f'{output_dir}/{output_prefix}river_arc_polygons.shp', driver="ESRI Shapefile")
            else:
                print(f'{mpi_print_prefix} Warning: total_arcs_cleaned_polys empty')
            del total_arcs_cleaned_polys[:]; del total_arcs_cleaned_polys

            cleaned_river_arcs = clean_river_arcs(river_arcs=river_arcs, total_arcs=total_arcs_cleaned)
            del river_arcs
            total_river_outline_polys = generate_river_outline_polys(river_arcs=cleaned_river_arcs)
            del cleaned_river_arcs
            if len(total_river_outline_polys) > 0:
                gpd.GeoDataFrame(index=range(len(total_river_outline_polys)), crs='epsg:4326', geometry=total_river_outline_polys).to_file(filename=f'{output_dir}/{output_prefix}river_outline.shp', driver="ESRI Shapefile")
            else:
                print(f'{mpi_print_prefix} Warning: total_river_outline_polys empty')
            del total_river_outline_polys

        del total_arcs_cleaned[:]; del total_arcs_cleaned
