import os
from copy import deepcopy
import numpy as np
from scipy import interpolate
import math
from osgeo import gdal
from dataclasses import dataclass
import pathlib
import pickle
from RiverMapGen.SMS import get_all_points_from_shp, SMS_ARC, SMS_MAP, \
    curvature, cpp2lonlat, lonlat2cpp, dl_cpp2lonlat, \
    get_perpendicular_angle


np.seterr(all='raise')

@dataclass
class dem_data():
    x: np.ndarray
    y: np.ndarray
    lon: np.ndarray
    lat: np.ndarray
    elev: np.ndarray
    dx: float
    dy: float

@dataclass
class Bombs():
    def __init__(self, center_x=None, center_y=None,
                 x=np.empty((0,1), dtype=float),
                 y=np.empty((0,1), dtype=float),
                 res=np.empty((0,1), dtype=float)):
        self.center = None
        self.points = None
        self.res = None

        if center_x is not None and center_y is not None:
            self.center = complex(center_x, center_y)
        if bool(x.shape) and x.shape==y.shape:
            self.points = np.squeeze(np.c_[x, y].view(np.complex128))
        if len(res)>0:
            self.res = res

    def __add__(self, other):
        if other is not None:
            self.points = np.r_[self.points, other.points]
            self.res = np.r_[self.res, other.res]
        return self

    def clean(self, bomb_radius_coef=0.5):
        # lon, lat = cpp2lonlat(self.points.real, self.points.imag)
        # if (abs(lon+93.72869978753)<0.00001).any():
        #     print()
        i_valid = np.ones((len(self.points),), dtype=bool)
        for i, _ in enumerate(self.points):
            if i_valid[i]:
                dist = abs(self.points[i] - self.points)
                nearby_idx = np.where(dist<bomb_radius_coef*self.res[i])[0]
                nearby_idx = nearby_idx[nearby_idx!=i]
                i_valid[nearby_idx] = False

        if sum(i_valid) == 0:
            raise Exception("impossible: all points bombed.")

        self.points = self.points[i_valid]

        return np.c_[self.points.real, self.points.imag]

def getAngle(a, b, c):
    ang = math.degrees(math.atan2(c[1]-b[1], c[0]-b[0]) - math.atan2(a[1]-b[1], a[0]-b[0]))
    return ang + 360 if ang < 0 else ang


def get_elev_from_tiles(x_cpp, y_cpp, tile_list):
    '''
    x: vector of x coordinates, assuming cpp;
    y: vector of x coordinates, assuming cpp;
    tile_list: list of tiles
    '''

    lon, lat = cpp2lonlat(x_cpp, y_cpp)

    elevs = np.empty(lon.shape, dtype=float); elevs.fill(np.nan)
    for S in tile_list:
        [j, i], in_box = Sidx(S, lon, lat)
        idx = (np.isnan(elevs) * in_box).astype(bool)  # only update valid entries that are not already set and in DEM box
        elevs[idx] = S.elev[i[idx], j[idx]]

    if np.isnan(elevs).any():
        pass
        # raise Exception('failed to find elevation')
        return None
    else:
        return elevs

def Sidx(S, lon, lat):
    '''
    return nearest index (i, j) in DEM mesh for point (x, y),
    assuming lon/lat, not projected coordinates
    '''
    dSx = S.lon[1] - S.lon[0]
    dSy = S.lat[1] - S.lat[0]
    i = (np.round((lon - S.lon[0]) / dSx)).astype(int)
    j = (np.round((lat - S.lat[0]) / dSy)).astype(int)

    valid = (i < S.lon.shape) * (j < S.lat.shape) * (i >= 0) * (j >= 0)
    return [i, j], valid

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

# %%
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

    ''' one tile
    [j, i], valid = Sidx(S, x, y)
    if all(valid):
        elevs = S.elev[i, j]
    else:
        return None, None
    '''

    eta_stream = np.tile(thalweg_eta.reshape(-1, 1), (1, search_steps))  # expanded to the search area

    ''' one tile
    [jj, ii], valid = Sidx(S, xts[:], yts[:])
    if valid.all:
        elevs = S.elev[ii, jj]
    else:
        return None, None
    '''
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

def get_dist_increment(line):
    line_copy = deepcopy(line)
    line_cplx = np.squeeze(line_copy.view(np.complex128))
    dist = np.absolute(line_cplx[1:] - line_cplx[:-1])

    # return np.r_[0.0, np.cumsum(dist)]
    return dist

def get_angle_diffs(xs, ys):
    line = np.c_[xs, ys]
    line_cplx = np.squeeze(line.view(np.complex128))
    angles = np.angle(np.diff(line_cplx))
    angle_diff0 = np.diff(angles)
    angle_diff = np.diff(angles)
    angle_diff[angle_diff0 > np.pi] -= 2 * np.pi
    angle_diff[angle_diff0 < -np.pi] += 2 * np.pi

    return angle_diff

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
        if np.sum(angle_diffs)/len(x[idx]) > 0.75:
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

def nudge_bank(line, perp, xs, ys, dist=np.array([35, 500])):
    ds = ((line[:, 0] - xs)**2 + (line[:, 1] - ys)**2)**0.5

    idx = ds < dist[0]
    xs[idx] = line[idx, 0] + dist[0] * np.cos(perp[idx])
    ys[idx] = line[idx, 1] + dist[0] * np.sin(perp[idx])

    idx = ds > dist[1]
    xs[idx] = line[idx, 0] + dist[1] * np.cos(perp[idx])
    ys[idx] = line[idx, 1] + dist[1] * np.sin(perp[idx])

    return xs, ys

def Tif2XYZ(tif_fname=None, cache=True):
    cache_name = tif_fname + '.pkl'
    if not cache:
        if os.path.exists(cache_name):
            os.remove(cache_name)

    if os.path.exists(cache_name):
        with open(cache_name, 'rb') as f:
            S = pickle.load(f)
        pass
    else:
        ds = gdal.Open(tif_fname, gdal.GA_ReadOnly)
        band = ds.GetRasterBand(1)

        width = ds.RasterXSize
        height = ds.RasterYSize

        gt = ds.GetGeoTransform()
        TL_x, TL_y = gt[0], gt[3]

        #showing a 2D image of the topo
        # plt.imshow(elevation, cmap='gist_earth',extent=[minX, maxX, minY, maxY])
        # plt.show()

        z = band.ReadAsArray()

        dx = gt[1]
        dy = gt[5]
        if gt[2] != 0 or gt[4] != 0:
            raise Exception()

        x_idx = np.array(range(width))
        y_idx = np.array(range(height))
        xp = dx * x_idx + TL_x + dx/2
        yp = dy * y_idx + TL_y + dy/2

        S = dem_data(xp, yp, xp, yp, z, dx, dy)
        with open(cache_name, 'wb') as f:
            pickle.dump(S, f, protocol=pickle.HIGHEST_PROTOCOL)

        ds = None  # close dataset

    return S

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

class river_arc():
    def __init__(self, line=None):
        self.thalweg = line
        self.perp = get_perpendicular_angle(line)
        self.endpoints = np.array([line[0, :], line[-1, :]])

class river_seg():
    def __init__(self, thalweg: river_arc):
        self.thalweg = thalweg

    def make_inner_arcs(self, narcs=2):
        inner_arcs = np.linspace(self.left_bank, self.right_bank, narcs)

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

def redistribute_arc(line, line_smooth, channel_width, smooth_option=1, R_coef=0.4, length_width_ratio=6.0, reso_thres=[3, 300]):
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
    river_resolution = np.minimum(R_coef * R, length_width_ratio * channel_width)
    river_resolution = np.minimum(np.maximum(reso_thres[0], river_resolution), reso_thres[1])
    # increase resolution near endpoints for better intersections
    endpoints_scale = 3.0
    starting_points = np.r_[0.0, np.cumsum(dist_along_thalweg)] < sum(river_resolution[:2])
    river_resolution[starting_points] /= endpoints_scale
    ending_points = np.flip(np.r_[0.0, np.cumsum(np.flip(dist_along_thalweg))]) < sum(river_resolution[-2:])
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

def nearest_neighbour(points_a, points_b):
    from scipy import spatial
    tree = spatial.cKDTree(points_b)
    return tree.query(points_a)[0], tree.query(points_a)[1]

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

def bomb_line(line, blast_radius, thalweg_id, i_check_valid=False):
    valid_idx = np.ones(len(line[:, 0]), dtype=bool)
    valid_idx_headtail = np.ones((len(line[:, 0]), 2), dtype=bool)
    if i_check_valid:
        for k in [0, -1]:
            valid_idx_headtail[:, -k] = (line[:, 0] - line[k, 0])**2 + (line[:, 1] - line[k, 1])**2 >= blast_radius[k]**2
            valid_idx *= valid_idx_headtail[:, -k]
        if sum(valid_idx) < 1:
            print(f'warning: thalweg {thalweg_id+1} has less than 1 points after bombing, neglecting ...')
            return valid_idx, valid_idx_headtail  # no changes

    return valid_idx, valid_idx_headtail

def width2narcs(width):
    width /= 4.0
    return int(3 + np.floor(0.5*width**0.25))

def make_river_map(
    tif_fnames = [],
    thalweg_shp_fname = '',
    thalweg_smooth_shp_fname = None,
    selected_thalweg = None,
    cache_folder=None,
    output_dir = './',
    output_prefix = '',
    mpi_print_prefix = ''
):
    '''
    [Core routine for making river maps]

    tif_fnames: a list of TIF file names. These TIFs should cover the area of interest and be arranged by priority (higher priority ones in front)

    thalweg_shp_fname: name of a polyline shapefile containing the thalwegs

    thalweg_smooth_shp_fname: name of a polyline shapefile containing the smoothed thalwegs (e.g., pre-processed by GIS tools or SMS)

    selected_thalweg: indices of selected thalwegs for which the river arcs will be sought.

    cache_folder: must specify one. Some intermediate information based on the input files will be saved for a faster execution in the future.

    output_dir: must specify one.

    output_prefix: a prefix of the output files, mainly used by the caller of this script; can be empty

    mpi_print_prefix: a prefix string to identify the calling mpi processe in the output messages; can be empty
    '''

    # ------------------------- other input parameters not exposed to user ---------------------------
    MapUnit2METER = 1  #00000

    river_threshold = np.array([16, 400]) / MapUnit2METER
    length_width_ratio = 8.0

    nudge_ratio = np.array((0.3, 2.0))  # ratio between nudging distance to mean half-channel-width

    i_close_poly = True

    i_blast_intersection = True
    blast_radius_scale = 0.6  # coef controlling the blast radius at intersections
    intersect_res_scale  = 0.4  # coef controlling the resolution of the paved mesh at intersections
    bomb_radius_coef = 1.2  # coef controlling the spacing among intersection joints

    starndard_watershed_resolution = 400.0  # meters
    max_nrow_arcs = width2narcs(river_threshold[-1])  # maximum number of arcs to resolve a channel (including bank arcs)

    i_thalweg_cache = False
    i_dem_cache = False
    # ------------------------- end basic inputs ---------------------------

    if i_thalweg_cache:
        if cache_folder is None or not os.path.exists(cache_folder):
            raise Exception("Cache folder not found.")

    # ------------------------- read DEM ---------------------------
    main_dem_id = 1
    S_list = []

    nvalid_tile = 0
    for i, tif_fname in enumerate(tif_fnames):
        if tif_fname is None:
            continue
        else:
            nvalid_tile += 1

        if pathlib.Path(tif_fname).suffix == ".tif" :
            S = Tif2XYZ(tif_fname=tif_fname, cache=i_dem_cache)
        else:
            raise Exception("Unknown DEM format.")
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
        raise Exception('Fatal Error: no valid DEM tiles')

    # ------------------------- read thalweg ---------------------------
    xyz, l2g, curv, _ = get_all_points_from_shp(thalweg_shp_fname, iCache=i_thalweg_cache, cache_folder=cache_folder)
    xyz[:, 0], xyz[:, 1] = lonlat2cpp(xyz[:, 0], xyz[:, 1])

    # Optional: provide a smoothed thalweg (on the 2D plane, not smoothed in z) to guide vertices distribution.
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
    thalweg_widths = [None] * len(thalwegs)
    valid_thalwegs = [True] * len(thalwegs)
    original_banks = [None] * len(thalwegs) * 2

    for i, [line, curv] in enumerate(zip(thalwegs, thalwegs_curv)):
        # print(f'Dry run: Arc {i+1} of {len(thalwegs)}')

        elevs = get_elev_from_tiles(line[:, 0], line[:, 1], S_list)
        if elevs is None:
            print(f"{mpi_print_prefix} warning: some elevs not found on thalweg {i+1}, the thalweg will be neglected ...")
            valid_thalwegs[i] = False
            continue

        # set water level at each point along the thalweg, based on observation, simulation, estimation, etc.
        thalweg_eta = set_eta_thalweg(line[:, 0], line[:, 1], elevs)

        x_banks_left, y_banks_left, x_banks_right, y_banks_right, _, width = \
            get_two_banks(S_list, line, thalweg_eta, search_length, search_steps, min_width=river_threshold[0])
        thalweg_widths[i] = width
        if width is None:
            thalweg_endpoints_width[i*2] = 0.0
            thalweg_endpoints_width[i*2+1] = 0.0
        else:
            thalweg_endpoints_width[i*2] = width[0]
            thalweg_endpoints_width[i*2+1] = width[-1]

        if len(line[:, 0]) < 2:
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
    # initialize some lists and array to hold the arc information
    bank_arcs = np.empty((len(thalwegs), 2), dtype=object) # left bank and right bank for each thalweg
    bank_arcs_raw = deepcopy(bank_arcs)
    bank_arcs_final = deepcopy(bank_arcs)
    inner_arcs = np.empty((len(thalwegs), max_nrow_arcs), dtype=object)
    cc_arcs = deepcopy(bank_arcs)  # [, 0] is head, [, 1] is tail
    smoothed_thalwegs = [None] * len(thalwegs)
    redistributed_thalwegs = [None] * len(thalwegs)
    corrected_thalwegs = [None] * len(thalwegs)
    centerlines = [None] * len(thalwegs)
    final_thalwegs = [None] * len(thalwegs)
    intersection_res_scatters = []
    thalwegs_neighbors = deepcopy(bank_arcs)  # [, 0] is head, [, 1] is tail
    real_bank_width = np.zeros((len(thalwegs), 2), dtype=float)  # [, 0] is head, [, 1] is tail
    bombed_points = np.empty((0, 3), dtype=float)  # left bank and right bank for each thalweg
    bombs = [None]* len(thalwegs) *2 # left bank and right bank for each thalweg

    # enumerate each thalweg
    for i, [thalweg, curv, width, valid_thalweg, thalweg_smooth] in enumerate(zip(thalwegs, thalwegs_curv, thalweg_widths, valid_thalwegs, thalwegs_smooth)):
        # print(f'{mpi_print_prefix} Wet run: Arc {i+1} of {len(thalwegs)}')

        if not valid_thalweg:
            print(f"{mpi_print_prefix} Thalweg {i} marked as invalid in dry run, skipping ...")
            continue

        # Redistribute thalwegs vertices
        thalweg, thalweg_smooth, reso, retained_idx = redistribute_arc(thalweg, thalweg_smooth, width, length_width_ratio=length_width_ratio, smooth_option=1)
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
        thalweg, thalweg_smooth, reso, retained_idx = redistribute_arc(thalweg, thalweg_smooth[retained_idx], width, length_width_ratio=length_width_ratio, smooth_option=1)
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
        thalweg_resolution = get_dist_increment(thalweg)

        # make inner arcs between two banks
        this_nrow_arcs = min(max_nrow_arcs, width2narcs(np.mean(width)))
        x_inner_arcs = np.linspace(x_banks_left, x_banks_right, this_nrow_arcs)
        y_inner_arcs = np.linspace(y_banks_left, y_banks_right, this_nrow_arcs)

        z_centerline = this_nrow_arcs + width/1e4

        # determine blast radius based on mean channel width at an intersection
        blast_radius = np.array([0.0, 0.0])
        valid_points = np.ones(x_banks_left.shape, dtype=bool)
        for k in [0, -1]:  # head and tail
            dist = thalweg_endpoints[:, :] - thalweg[k, :]
            neighbor_thalwegs_endpoints = np.argwhere(dist[:, 0]**2 + dist[:, 1]**2 < 200**2)
            thalwegs_neighbors[i, k] = neighbor_thalwegs_endpoints
            if len(neighbor_thalwegs_endpoints) > 1:
                blast_radius[k] = blast_radius_scale * np.mean(thalweg_endpoints_width[neighbor_thalwegs_endpoints])
            else:  # no intersections
                blast_radius[k] = 0.0

        # bomb intersections
        valid_l, valid_l_headtail = bomb_line(np.c_[x_banks_left, y_banks_left], blast_radius, i, i_check_valid=True)
        valid_r, valid_r_headtail = bomb_line(np.c_[x_banks_right, y_banks_right], blast_radius, i, i_check_valid=True)
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
        if river_quality(x_inner_arcs, y_inner_arcs, valid_points):
            # save centerlines
            # assemble inner arcs
            bombs_xyz = [np.empty((0,3), dtype=float)] * 2
            for k, [x_inner_arc, y_inner_arc] in enumerate(zip(x_inner_arcs, y_inner_arcs)):
                line = np.c_[x_inner_arc, y_inner_arc]
                if sum(valid_points) > 0:
                    # snap vertices too close to each other
                    line[valid_points, :] = snap_vertices(line[valid_points, :], width[valid_points] * 0.3)  # optional: thalweg_resolution*0.75

                    # ----------Save-------
                    # save final bank arcs
                    if k == 0:  # left bank
                        bank_arcs_final[i, 0] = SMS_ARC(points=np.c_[line[valid_points, 0], line[valid_points, 1], z_centerline[valid_points]], src_prj='cpp')
                    elif k == len(x_inner_arcs)-1:  # right bank
                        bank_arcs_final[i, 1] = SMS_ARC(points=np.c_[line[valid_points, 0], line[valid_points, 1], z_centerline[valid_points]], src_prj='cpp')
                    # save inner arcs
                    inner_arcs[i, k] = SMS_ARC(points=np.c_[line[valid_points, 0], line[valid_points, 1], z_centerline[valid_points]], src_prj='cpp')
                    # save centerline
                    if k == int(len(x_inner_arcs)/2):
                        centerlines[i] = SMS_ARC(points=np.c_[line[valid_points, 0], line[valid_points, 1], z_centerline[valid_points]], src_prj='cpp')
                    # save bombed points
                    bombed_points = np.r_[bombed_points, np.c_[line[bombed_idx, 0], line[bombed_idx, 1], width[bombed_idx]/this_nrow_arcs]]

                # test bombs
                for l in [0, 1]:
                    if sum(~valid_points_headtail[:, l]) > 0:
                        bombs_xyz[l] = np.r_[bombs_xyz[l], np.c_[line[~valid_points_headtail[:, l]][:, :2], width[~valid_points_headtail[:, l]]/this_nrow_arcs]]

            for l in [0, 1]:
                if len(bombs_xyz[l]) > 0:
                    bombs[2*i+l] = Bombs(x=bombs_xyz[l][:, 0], y=bombs_xyz[l][:, 1], res=bombs_xyz[l][:, 2])

            # assemble cross-channel arcs
            if sum(valid_points) > 0:
                for j in [0, -1]:
                    cc_arcs[i, j] = SMS_ARC(points=np.c_[x_inner_arcs[:, valid_points][:, j], y_inner_arcs[:, valid_points][:, j]], src_prj='cpp')
                    pass

    # assemble intersection resolution scatters
    for i, thalweg_neighbors in enumerate(thalwegs_neighbors):
        for j, neibs in enumerate(thalweg_neighbors):  # head and tail
            if neibs is not None and len(neibs) > 1:
                # intersect_ring is the effective refining region, consisting of bank endpoints at the intersection of multiple thalwegs
                intersect_ring = np.zeros((0, 3), dtype=float)
                for nei in neibs:
                    id = int(nei/2)
                    i_tail_head = int(nei%2)
                    if bank_arcs[id, 0] is not None:
                        intersect_ring = np.r_[intersect_ring, np.c_[bank_arcs[id, 0].nodes[i_tail_head, :2].reshape(1, -1), intersect_res_scale * real_bank_width[id, i_tail_head]]]
                    if bank_arcs[id, 1] is not None:
                        intersect_ring = np.r_[intersect_ring, np.c_[bank_arcs[id, 1].nodes[i_tail_head, :2].reshape(1, -1), intersect_res_scale * real_bank_width[id, i_tail_head]]]

                # make a scatter set consisting of a center point (where further refining is optional) and an outer ring (where resolution reverts to normal value)
                if len(intersect_ring) > 0:
                    center_point = np.mean(intersect_ring, axis=0).reshape(1, -1)
                    if bombs[2*i+j] is not None:
                        bombs[2*i+j].center = complex(center_point[0, 0], center_point[0, 1])

                    diff = intersect_ring - center_point
                    radius = np.maximum((diff[:, 0]**2 + diff[:, 1]**2)**0.5, 3.0)  # set minimum to 3.0 m to avoid division by 0
                    stretch = (np.maximum(1.2, radius/np.mean(radius))/(radius/np.mean(radius))).reshape(-1, 1)

                    # outter ring is nudged to 120% based on the intersection ring
                    outter_ring = (intersect_ring - center_point) * 1.2 + center_point
                    # nudge the outter ring points to make the outter ring more or less a circle, avoiding outter_ring points inside intersection ring
                    diff = outter_ring - center_point
                    outter_ring = np.tile(stretch, (1,3)) * (outter_ring - center_point) + center_point

                    outter_ring[:, -1] = starndard_watershed_resolution

                    # assemble all rings
                    all = np.r_[intersect_ring, center_point, outter_ring]
                    intersection_res_scatters.append(all)

    intersect_res = None
    if intersection_res_scatters:  # len > 0
        intersect_res = np.r_[np.concatenate(intersection_res_scatters, axis=0), bombed_points]
        valid = intersect_res[:, -1] > 0
        intersect_res = intersect_res[valid, :]
        # conver to lon/lat for outputs
        intersect_res_lonlat = np.c_[intersect_res[:, :2], dl_cpp2lonlat(intersect_res[:, 2], intersect_res[:, 1])]

    # assemble river arcs
    # river_arcs = np.r_[bank_arcs.reshape((-1, 1)), inner_arcs.reshape((-1, 1))]
    # End wet run

    # bombed_lon, bomed_lat = cpp2lonlat(bombed_points[:, 0], bombed_points[:, 1])

    for i, thalweg_neighbors in enumerate(thalwegs_neighbors):
        for j, neibs in enumerate(thalweg_neighbors):  # head and tail
            if bombs[2*i+j] is not None:
                for nei in neibs:
                    nei = int(nei)
                    if (2*i+j) != nei:
                        bombs[2*i+j] += bombs[nei]
                        bombs[nei] = None


    bombed_xy = np.empty((0,2), dtype=float)
    for i, bomb in enumerate(bombs):
        if bomb is not None:
            bomb.clean(bomb_radius_coef=bomb_radius_coef)
            bombed_xy = np.r_[bombed_xy, np.c_[bomb.points.real, bomb.points.imag]]
    bombed_lon, bomed_lat = cpp2lonlat(bombed_xy[:, 0], bombed_xy[:, 1])

    # ------------------------- write SMS maps ---------------------------
    if any(bank_arcs.flatten()):  # not all arcs are None
        SMS_MAP(arcs=bank_arcs.reshape((-1, 1))).writer(filename=f'{output_dir}/{output_prefix}bank.map')
        SMS_MAP(arcs=cc_arcs.reshape((-1, 1))).writer(filename=f'{output_dir}/{output_prefix}cc_arcs.map')
        SMS_MAP(arcs=bank_arcs_raw.reshape((-1, 1))).writer(filename=f'{output_dir}/{output_prefix}bank_raw.map')
        SMS_MAP(arcs=bank_arcs_final.reshape((-1, 1))).writer(filename=f'{output_dir}/{output_prefix}bank_final.map')
        SMS_MAP(arcs=inner_arcs.reshape((-1, 1))).writer(filename=f'{output_dir}/{output_prefix}inner_arcs.map')
        SMS_MAP(detached_nodes=bombed_points).writer(filename=f'{output_dir}/{output_prefix}relax_points.map')
        SMS_MAP(arcs=smoothed_thalwegs).writer(filename=f'{output_dir}/{output_prefix}smoothed_thalweg.map')
        SMS_MAP(arcs=redistributed_thalwegs).writer(filename=f'{output_dir}/{output_prefix}redist_thalweg.map')
        SMS_MAP(arcs=corrected_thalwegs).writer(filename=f'{output_dir}/{output_prefix}corrected_thalweg.map')
        SMS_MAP(arcs=centerlines).writer(filename=f'{output_dir}/{output_prefix}centerlines.map')
        SMS_MAP(arcs=final_thalwegs).writer(filename=f'{output_dir}/{output_prefix}final_thalweg.map')
        if intersect_res is not None:
            np.savetxt(f'{output_dir}/{output_prefix}intersection_res.xyz', intersect_res_lonlat)

        if i_close_poly:
            total_arcs = np.r_[inner_arcs.reshape((-1, 1)), cc_arcs.reshape((-1, 1))]
        else:
            total_arcs = np.r_[inner_arcs.reshape((-1, 1))]
        SMS_MAP(arcs=total_arcs).writer(filename=f'{output_dir}/{output_prefix}total_arcs.map')
        SMS_MAP(detached_nodes=np.c_[bombed_lon, bomed_lat]).writer(f'{output_dir}/{output_prefix}total_intersection_joints.map')
    else:
        print(f'{mpi_print_prefix} No arcs found, aborting writing to *.map')
