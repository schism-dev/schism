from logging import raiseExceptions
import pickle
import os
import numpy as np
from shapely.geometry import LineString
# import matplotlib.pyplot as plt
import glob
import numpy as np
import re
import shapefile
import geopandas as gpd
from pathlib import Path


def lonlat2cpp(lon, lat, lon0=0, lat0=0):
    R = 6378206.4

    lon_radian, lat_radian = lon/180*np.pi, lat/180*np.pi
    lon0_radian, lat0_radian = lon0/180*np.pi, lat0/180*np.pi

    xout = R * (lon_radian - lon0_radian) * np.cos(lat0_radian)
    yout = R * lat_radian

    return [xout, yout]

def dl_cpp2lonlat(dl, lat0=0):
    R = 6378206.4

    lat0_radian = lat0/180*np.pi

    dlon_radian = dl/R/np.cos(lat0_radian)
    
    dlon = dlon_radian*180/np.pi

    return dlon

def cpp2lonlat(x, y, lon0=0, lat0=0):
    R = 6378206.4

    lon0_radian, lat0_radian = lon0/180*np.pi, lat0/180*np.pi

    lon_radian = lon0_radian + x/R/np.cos(lat0_radian)
    lat_radian = y / R
    
    lon, lat = lon_radian*180/np.pi, lat_radian*180/np.pi 

    return [lon, lat]

def normalizeVec(x, y):
    distance = np.sqrt(x*x+y*y)
    return x/distance, y/distance

def makeOffsetPoly(oldX, oldY, offset, outer_ccw = 1):
    num_points = len(oldX)
    newX = []
    newY = []

    for curr in range(num_points):
        prev = (curr + num_points - 1) % num_points
        next = (curr + 1) % num_points

        vnX =  oldX[next] - oldX[curr]
        vnY =  oldY[next] - oldY[curr]
        vnnX, vnnY = normalizeVec(vnX,vnY)
        nnnX = vnnY
        nnnY = -vnnX

        vpX =  oldX[curr] - oldX[prev]
        vpY =  oldY[curr] - oldY[prev]
        vpnX, vpnY = normalizeVec(vpX,vpY)
        npnX = vpnY * outer_ccw
        npnY = -vpnX * outer_ccw

        bisX = (nnnX + npnX) * outer_ccw
        bisY = (nnnY + npnY) * outer_ccw

        bisnX, bisnY = normalizeVec(bisX,  bisY)
        bislen = offset /  np.sqrt(1 + nnnX*npnX + nnnY*npnY)

        newX.append(oldX[curr] + bislen * bisnX)
        newY.append(oldY[curr] + bislen * bisnY)

    return newX, newY

def redistribute(x, y, length=None, num_points=None, iplot=False):
    line = LineString(np.c_[x, y])

    if length is None and num_points is None:
      raise Exception("Needs to specify either length or num_points")

    if length is not None:
        num_points = max(2, int(line.length / length))
    
    new_points = [line.interpolate(i/float(num_points - 1), normalized=True) for i in range(num_points)]
    x_subsampled = [p.x for p in new_points]
    y_subsampled = [p.y for p in new_points]

    # if iplot:
    #     plt.plot(x, y, '+')
    #     plt.plot(x_subsampled, y_subsampled, 'o')
    #     plt.axis('equal')
    #     plt.show()

    return x_subsampled, y_subsampled, new_points

'''
<Sample map only containing arcs>
MAP VERSION 8
BEGCOV
COVFLDR "Area Property"
COVNAME "Area Property"
COVELEV 0.000000
COVID 26200
COVGUID 57a1fdc1-d908-44d3-befe-8785288e69e7
COVATTS VISIBLE 1
COVATTS ACTIVECOVERAGE Area Property
COV_WKT GEOGCS["GCS_WGS_1984",DATUM["WGS84",SPHEROID["WGS84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]END_COV_WKT 
COV_VERT_DATUM 0
COV_VERT_UNITS 0
COVATTS PROPERTIES MESH_GENERATION
NODE
XY -28.76 73.25 0.0
ID 1
END
NODE
XY -49.76 43.09 0.0
ID 2
END
NODE
XY -24.91 56.11 0.0
ID 3
END
NODE
XY -31.34 44.28 0.0
ID 4
END
ARC
ID 1
ARCELEVATION 0.000000
NODES        1        2
ARCVERTICES 6
-4.020000000000001 75.980000000000004 0.000000000000000
27.030000000000001 68.600000000000009 0.000000000000000
45.730000000000004 45.450000000000003 0.000000000000000
40.829999999999998 26.420000000000002 0.000000000000000
8.830000000000000 18.789999999999999 0.000000000000000
-37.609999999999999 30.300000000000001 0.000000000000000
DISTNODE 0
NODETYPE 0
ARCBIAS 1
MERGED 0 0 0 0 
END
ARC
ID 2
ARCELEVATION 0.000000
NODES        3        4
ARCVERTICES 3
-4.300000000000000 57.250000000000000 0.000000000000000
-0.160000000000000 45.130000000000003 0.000000000000000
-20.059999999999999 39.300000000000004 0.000000000000000
DISTNODE 0
NODETYPE 0
ARCBIAS 1
MERGED 0 0 0 0 
END
ENDCOV
BEGTS
LEND
'''

def merge_maps(mapfile_glob_str, merged_fname):
    map_file_list = glob.glob(mapfile_glob_str) 
    if len(map_file_list) > 0:
        map_objects = [SMS_MAP(filename=map_file) for map_file in map_file_list]

        total_map = map_objects[0]
        for map_object in map_objects[1:]:
            total_map += map_object
        total_map.writer(merged_fname)
    else:
        print(f'failed to combine {mapfile_glob_str}, no files found')


class SMS_ARC():
    '''class for manipulating arcs in SMS maps''' 
    def __init__(self, points=None, node_idx=[0, -1], src_prj=None, dst_prj='epsg:4326'):
        # self.isDummy = (len(points) == 0)

        if src_prj is None:
            raise Exception('source projection not specified when initializing SMS_ARC')

        if src_prj == 'cpp' and dst_prj == 'epsg:4326':
            points[:, 0], points[: ,1] = cpp2lonlat(points[:, 0], points[: ,1])

        npoints, ncol = points.shape
        self.points = np.zeros((npoints, 3), dtype=float)
        self.points[:, :min(3, ncol)] = points[:, :min(3, ncol)]

        self.nodes = self.points[node_idx, :]
        self.arcvertices = np.delete(self.points, node_idx, axis=0)
        self.arcnode_glb_ids = np.empty(self.nodes[:, 0].shape, dtype=int)

        self.arc_hats = np.zeros((4, 3), dtype=float)
        self.arc_hat_length = -1

    def make_hats(self, arc_hat_length=-1):
        if arc_hat_length <= 0:
            raise Exception('Arc hat length <= 0')
        else:
            self.arc_hat_length = arc_hat_length
        
        # make hats (a perpendicular line at each of the arc ends)
        for i, [x0, y0, xx, yy] in enumerate([
            [self.points[0, 0], self.points[0, 1], self.points[1, 0], self.points[1, 1]],
            [self.points[-1, 0], self.points[-1, 1], self.points[-2, 0], self.points[-2, 1]],
        ]):
            xt = xx - x0
            yt = yy - y0
            st = (xt**2 + yt**2)**0.5
            xt = xt/st*arc_hat_length/2
            yt = yt/st*arc_hat_length/2

            self.arc_hats[2*i, 0] = x0 - yt
            self.arc_hats[2*i, 1] = y0 + xt
            self.arc_hats[2*i+1, 0] = x0 + yt
            self.arc_hats[2*i+1, 1] = y0 - xt
        
        # import matplotlib.pyplot as plt
        # plt.scatter(self.points[:, 0], self.points[:, 1])
        # plt.scatter(self.arc_hats[:, 0], self.arc_hats[:, 1])
        # plt.gca().set_aspect('equal', adjustable='box')
        # plt.show()

        return [SMS_ARC(points=self.arc_hats[:2, :]), SMS_ARC(points=self.arc_hats[2:, :])]
   
class SMS_MAP():
    '''class for manipulating SMS maps''' 
    def __init__(self, filename=None, arcs=[], detached_nodes=[], epsg=4326):
        self.epsg = None
        self.arcs = []
        self.nodes = np.zeros((0, 3))
        self.detached_nodes = np.zeros((0, 3))
        self.valid = True

        self.epsg = epsg

        if filename is not None:
            self.reader(filename=filename)
        else:
            if type(arcs) == np.ndarray:
                arcs = np.squeeze(arcs).tolist()
            arcs = list(filter(lambda item: item is not None, arcs))
            if arcs == [] and detached_nodes==[]:
                self.valid = False
            
            self.arcs = arcs
            self.detached_nodes = detached_nodes
    
    def __add__(self, other):
        self.arcs = self.arcs + other.arcs
        self.detached_nodes = np.r_[self.detached_nodes, other.detached_nodes]
        return SMS_MAP(arcs=self.arcs, detached_nodes=self.detached_nodes, epsg=self.epsg)
    
    def get_xyz(self):
        self.n_xyz = 0
        self.l2g = []

        for arc in self.arcs:
            self.l2g.append(np.array(np.arange(self.n_xyz, self.n_xyz+len(arc.points))))
            self.n_xyz += len(arc.points)
        
        self.xyz = np.zeros((self.n_xyz, 3), dtype=float)
        for ids, arc in zip(self.l2g, self.arcs):
            self.xyz[ids, :] = arc.points

        return self.xyz, self.l2g
    
    def reader(self, filename='test.map'):
        self.n_glb_nodes = 0
        self.n_arcs = 0
        self.n_detached_nodes = 0

        arc_nodes = []
        with open(filename) as f:
            while True:
                line = f.readline()
                if not line:
                    break

                strs = re.split(' +', line.strip())
                if strs[0] == 'COV_WKT':
                    if "GCS_WGS_1984" in line:
                        self.epsg = 4326
                    else:
                        raiseExceptions('unkown epsg')
                elif strs[0] == 'NODE':
                    line = f.readline()
                    strs = re.split(' +', line.strip())
                    self.n_glb_nodes += 1
                    self.nodes = np.append(self.nodes, np.reshape([float(strs[1]), float(strs[2]), float(strs[3])], (1,3)), axis=0)
                elif line.strip() == 'POINT':
                    line = f.readline()
                    strs = re.split(' +', line.strip())
                    self.n_detached_nodes += 1
                    self.detached_nodes = np.append(self.detached_nodes, np.reshape([float(strs[1]), float(strs[2]), float(strs[3])], (1,3)), axis=0)
                elif line.strip() == 'ARC':
                    self.n_arcs += 1
                elif strs[0] == 'NODES':
                    this_arc_node_idx = np.array([int(strs[1]), int(strs[2])])-1
                elif strs[0] == 'ARCVERTICES':
                    this_arc_nvert = int(strs[1])
                    this_arc_verts = np.zeros((this_arc_nvert, 3), dtype=float)
                    for i in range(this_arc_nvert):
                        strs = f.readline().strip().split(' ')
                        this_arc_verts[i, :] = np.array([strs[0], strs[1], strs[2]])
                    node_1 = np.reshape(self.nodes[this_arc_node_idx[0], :], (1, 3))
                    node_2 = np.reshape(self.nodes[this_arc_node_idx[1], :], (1, 3))
                    this_arc = SMS_ARC(points=np.r_[node_1, this_arc_verts, node_2], src_prj=f'epsg: {self.epsg}')
                    self.arcs.append(this_arc)
                    arc_nodes.append(this_arc_node_idx[0])
                    arc_nodes.append(this_arc_node_idx[1])
        pass
    
    def writer(self, filename='test.map'):
        import os

        if not self.valid:
            print(f'No arcs found in map, aborting writing to *.map')
            return

        fpath = os.path.dirname(filename)
        if not os.path.exists(fpath):
            os.makedirs(fpath, exist_ok=True)

        with open(filename, 'w') as f:
            # write header
            f.write('MAP VERSION 8\n')
            f.write('BEGCOV\n')
            # f.write('COVFLDR "Area Property"\n')
            # f.write('COVNAME "Area Property"\n')
            # f.write('COVELEV 0.000000\n')
            f.write('COVID 26200\n')
            f.write('COVGUID 57a1fdc1-d908-44d3-befe-8785288e69e7\n')
            f.write('COVATTS VISIBLE 1\n')
            f.write('COVATTS ACTIVECOVERAGE Area Property\n')
            if self.epsg == 4326:
                f.write('COV_WKT GEOGCS["GCS_WGS_1984",DATUM["WGS84",SPHEROID["WGS84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]END_COV_WKT \n')
            elif self.epsg == 26918:
                f.write('COV_WKT PROJCS["NAD83 / UTM zone 18N",GEOGCS["NAD83",DATUM["North_American_Datum_1983",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6269"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4269"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-75],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","26918"]]END_COV_WKT')
            else:
                raiseExceptions("Projection not supported.")
            f.write('COV_VERT_DATUM 0\n')
            f.write('COV_VERT_UNITS 0\n')
            f.write('COVATTS PROPERTIES MESH_GENERATION\n')
            
            node_counter = 0
            for i, arc in enumerate(self.arcs):
                for j, node in enumerate(arc.nodes):
                    node_counter += 1
                    self.arcs[i].arcnode_glb_ids[j] = node_counter
                    f.write('NODE\n')
                    f.write(f'XY {node[0]} {node[1]} {node[2]}\n')
                    f.write(f'ID {node_counter}\n')
                    f.write('END\n')

            for i, node in enumerate(self.detached_nodes):
                node_counter += 1
                f.write('POINT\n')
                f.write(f'XY {node[0]} {node[1]} 0.0\n')
                f.write(f'ID {node_counter}\n')
                f.write('END\n')

            for i, arc in enumerate(self.arcs):
                f.write('ARC\n')
                f.write(f'ID {i+1}\n')
                f.write('ARCELEVATION 0.00\n')
                f.write(f'NODES {" ".join(arc.arcnode_glb_ids.astype(str))}\n')
                f.write(f'ARCVERTICES {len(arc.arcvertices)}\n')
                for vertex in arc.arcvertices:
                    f.write(f'{vertex[0]} {vertex[1]} {vertex[2]}\n')
                f.write('END\n')
                    
            f.write('ENDCOV\n')
            f.write('BEGTS\n')
            f.write('LEND\n')
            pass

class Levee_SMS_MAP(SMS_MAP):
    def __init__(self, arcs=[], epsg=4326):
        super().__init__(arcs=arcs, epsg=epsg)
        self.centerline_list = arcs
        self.subsampled_centerline_list = []
        self.offsetline_list = []
    
    def make_levee_maps(self, offset_list=[-5, 5, -15, 15], subsample=[300, 10]):
        for arc in self.centerline_list:
            x_sub, y_sub, _ = redistribute(x=arc.points[:, 0], y=arc.points[:, 1], length=subsample[0])
            self.subsampled_centerline_list.append(SMS_ARC(points=np.c_[x_sub, y_sub]))
            
            for offset in offset_list:
                x_off, y_off = makeOffsetPoly(x_sub, y_sub, offset)
                self.offsetline_list.append(SMS_ARC(points=np.c_[x_off, y_off]))
        return SMS_MAP(arcs=self.subsampled_centerline_list), SMS_MAP(arcs=self.offsetline_list)

def get_perpendicular_angle(line):
    line_cplx = np.squeeze(line.copy().view(np.complex128))
    angles = np.angle(np.diff(line_cplx))
    angle_diff0 = np.diff(angles)
    angle_diff = np.diff(angles)
    angle_diff[angle_diff0 > np.pi] -= 2 * np.pi
    angle_diff[angle_diff0 < -np.pi] += 2 * np.pi
    perp = angles[:-1] + angle_diff / 2 - np.pi / 2
    perp = np.r_[angles[0] - np.pi / 2, perp, angles[-1] - np.pi / 2]

    return perp

def curvature(pts):
    if len(pts[:, 0]) < 3:
        cur = np.zeros((len(pts[:, 0])))
    else:
        dx = np.gradient(pts[:,0]) # first derivatives
        dy = np.gradient(pts[:,1])

        d2x = np.gradient(dx) #second derivatives
        d2y = np.gradient(dy)

        cur = np.abs(dx * d2y - d2x * dy) / ((dx * dx + dy * dy)**1.5 + 1e-20)

    return cur

def get_all_points_from_shp(fname, iNoPrint=True, iCache=False, cache_folder=None):
    if not iNoPrint: print(f'reading shapefile: {fname}')

    if cache_folder is None:
        cache_folder = ''

    cache_name = cache_folder + Path(fname).stem + '.pkl'

    if iCache == False:
        if os.path.exists(cache_name):
            os.remove(cache_name)

    if os.path.exists(cache_name):
        with open(cache_name, 'rb') as file:
            tmp_dict = pickle.load(file)
            xyz = tmp_dict['xyz']
            shape_pts_l2g = tmp_dict['shape_pts_l2g']
            curv = tmp_dict['curv']
            perp = tmp_dict['perp']
        if not iNoPrint: print(f'Cache of the shapefile loaded.')
    else:
        '''using pyshp
        sf = shapefile.Reader(fname)
        shapes = sf.shapes()

        shape_pts_l2g = []
        xyz = np.empty((0, 2), dtype=float)
        curv = np.empty((0, ), dtype=float)
        n = 0
        for i, shp in enumerate(shapes):
            pts = np.array(shp.points)
            curv = np.r_[curv, curvature(pts)]
            # pts_cplx = np.array(pts).view(np.complex128)
            # dl = abs(pts_cplx[2:-1] - pts_cplx[1:-2])

            # if not iNoPrint: print(f'shp {i+1} of {len(shapes)}, {len(pts)} points')

            xyz = np.append(xyz, shp.points, axis=0)
            shape_pts_l2g.append(np.array(np.arange(n, n+len(shp.points))))
            n += len(shp.points)
        '''
        
        # using geopandas, which seems more efficient than pyshp
        shapefile = gpd.read_file(fname)
        npts = 0
        nvalid_shps = 0
        for i in range(shapefile.shape[0]):
            try:
                shp_points = np.array(shapefile.iloc[i, :]['geometry'].coords.xy).shape[1]
            except:
                print(f"warning: shape {i+1} of {shapefile.shape[0]} is invalid")
                continue
            npts += shp_points
            nvalid_shps += 1

        xyz = np.zeros((npts, 2), dtype=float)
        shape_pts_l2g =[None] * nvalid_shps
        ptr = 0; ptr_shp = 0
        for i in range(shapefile.shape[0]):
            try:
                shp_points = np.array(shapefile.iloc[i, :]['geometry'].coords.xy).shape[1]
            except:
                print(f"warning: shape {i+1} of {shapefile.shape[0]} is invalid")
                continue

            xyz[ptr:ptr+shp_points] = np.array(shapefile.iloc[i, :]['geometry'].coords.xy).T
            shape_pts_l2g[ptr_shp] = np.array(np.arange(ptr, ptr+shp_points))
            ptr += shp_points; ptr_shp += 1;
        if ptr != npts or ptr_shp != nvalid_shps:
            raise Exception("number of shapes/points does not match")

        curv = np.empty((npts, ), dtype=float)
        perp = np.empty((npts, ), dtype=float)
        for i, _ in enumerate(shape_pts_l2g):
            line = xyz[shape_pts_l2g[i], :]
            curv[shape_pts_l2g[i]] = curvature(line)
            perp[shape_pts_l2g[i]] = get_perpendicular_angle(line)
        
        # for i in range(shapefile.shape[0]):
        #     try:
        #         shp_points = np.array(shapefile.iloc[i, :]['geometry'].coords.xy).shape[1]
        #     except NotImplementedError:
        #         print(f"Warning: multi-part geometries, neglecting ...")
        #         continue
        #     except:
        #         raiseExceptions(f'Undefined error reading shapefile {fname}')

        #     shape_pts_l2g.append(np.array(np.arange(npts, npts+shp_points)))
        #     npts += shp_points
        # xyz = np.empty((npts, 2), dtype=float)

        # curv = np.empty((npts, ), dtype=float)
        # perp = np.empty((npts, ), dtype=float)
        # for i in range(shapefile.shape[0]):
        #     xyz[shape_pts_l2g[i], :] = np.array(shapefile.iloc[i, :]['geometry'].coords.xy).T
        #     curv[shape_pts_l2g[i]] = curvature(xyz[shape_pts_l2g[i], :])
        #     perp[shape_pts_l2g[i]] = get_perpendicular_angle(xyz[shape_pts_l2g[i], :2])

        # if not iNoPrint: print(f'Number of shapes read: {len(shapes)}')

        with open(cache_name, 'wb') as file:
            tmp_dict = {'xyz': xyz, 'shape_pts_l2g': shape_pts_l2g, 'curv': curv, 'perp': perp}
            pickle.dump(tmp_dict, file)

    return xyz, shape_pts_l2g, curv, perp

def replace_shp_pts(inshp_fname, pts, l2g, outshp_fname):
    sf = shapefile.Reader(inshp_fname)
    shapes = sf.shapes()

    with shapefile.Writer(outshp_fname) as w:
        w.fields = sf.fields[1:] # skip first deletion field
        for i, feature in enumerate(sf.iterShapeRecords()): # iteration on both record and shape for a feature
            if len(l2g[i]) > 0:
                w.record(*feature.record) # * for unpacking tuple
                feature.shape.points = pts[l2g[i]]
                w.shape(feature.shape)

def extract_quad_polygons(input_fname='test.map', output_fname=None):
    if output_fname is None:
        output_fname = os.path.splitext(input_fname)[0] + '.quad.map'
    
    with open(output_fname, 'w') as fout:
        lines_buffer = []
        iPatch = False

        n_read = 0
        with open(input_fname) as f:
            while True:
                line = f.readline()
                if not line: break
                strs = re.split(' +', line.strip())

                if strs[0] != 'POLYGON':
                    fout.write(line)
                else:
                    lines_buffer.append(line)
                    while True:
                        line = f.readline()
                        strs = re.split(' +', line.strip())
                        lines_buffer.append(line)
                        if strs[0] == 'PATCH':
                            iPatch = True
                        elif strs[0] == 'END':
                            if iPatch:
                                for line_buffer in lines_buffer:
                                    fout.write(line_buffer)
                                iPatch = False
                                fout.flush()
                            lines_buffer = []
                            break


if __name__ == '__main__':
    # my_map = SMS_MAP(filename='test_z.map')
    # my_map.get_xyz()
    # my_map.writer('./test.map')

    merge_maps(f'/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/Parallel/Outputs/CUDEM_merged_thalwegs_1e6_single_fix_simple_sms_cleaned_32cores/*corrected_thalweg*.map', merged_fname=f'/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/Parallel/Outputs/CUDEM_merged_thalwegs_1e6_single_fix_simple_sms_cleaned_32cores/total_corrected_thalwegs.map')

    extract_quad_polygons(input_fname='/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/Parallel/SMS_proj/feeder/14.33.map')
