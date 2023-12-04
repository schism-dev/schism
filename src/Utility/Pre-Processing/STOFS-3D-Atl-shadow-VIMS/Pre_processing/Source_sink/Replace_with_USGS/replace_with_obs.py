#!/usr/bin/env python

import os
import time
import copy

import numpy as np
import shapefile
import pandas as pd
from datetime import datetime  #, timedelta
import matplotlib.pyplot as plt
import xarray as xr
import geopandas as gpd

import json
from schism_py_pre_post.Grid.Prop import Prop
from schism_py_pre_post.Timeseries.TimeHistory import TimeHistory
from schism_py_pre_post.Download.Data import ObsData
from schism_py_pre_post.Download.download_usgs_with_api import download_stations, usgs_var_dict, convert_to_ObsData
from pylib import schism_grid
import pickle
  

def BinA(A=None, B=None):
    import numpy as np
    if A is None and B is None:
        A = np.array([3,5,7,1,9,8,6,6])
        B = np.array([3,1,5,8,6])
    else:
        A = np.array(A)
        B = np.array(B)

    index = np.argsort(A)
    sorted_A = A[index]
    sorted_index = np.searchsorted(sorted_A, B)

    Bindex = np.take(index, sorted_index, mode="raise")
    mask = A[Bindex] != B

    result = np.ma.array(Bindex, mask=mask)
    return result

def read_nwm_file(filename):
    n = 0
    while n < 10:   # try 10 times
        try:
            d = xr.open_dataset(filename)
            fID_all, stream_flow_all  = copy.deepcopy(d.feature_id.data), copy.deepcopy(d.streamflow.data)
            d.close()
            return [fID_all, stream_flow_all]
        except RuntimeError:
            print(f"warning: failed to read netcdf file {filename}, Attempt {n}")
        except:
            raise Exception(f'failed to read netcdf file {filename}')
        n += 1
    
    raise Exception(f'Failed to read {filename} after 10 trys')



def read_nwm_data(data_dir=None, fIDs=[], date_range=pd.date_range("2015-06-01 00:00:00", "2015-08-01 00:00:00", freq="60min")):
    # assemble data file list
    datafiles = []
    for this_date in date_range:
        datafiles.append(f'{data_dir}/{datetime.strftime(this_date, "%Y%m%d%H%M")}.CHRTOUT_DOMAIN1.comp')

    # check file availability
    for this_file in datafiles:
        if not os.path.exists(this_file):
            raise Exception(f'{this_file} not found\n')

    # read data
    stream_flow = np.empty((len(datafiles), len(fIDs)))
    stream_flow[:] = np.nan
    for i, this_file in enumerate(datafiles):
        print(f'reading {this_file}\n')

        fID_all, stream_flow_all  = read_nwm_file(this_file)

        fIDs_in_All = BinA(A=fID_all, B=fIDs) 
        if np.ma.is_masked(fIDs_in_All):
            raise Exception('failed to find some segments in NWM data')

        stream_flow[i, :] = stream_flow_all[fIDs_in_All]

    return stream_flow


class usgs_st():
    def __init__(self, st_id=None, nearby_nwm_fid=None):
        self.st_id = st_id
        self.nearby_nwm_fid = nearby_nwm_fid


class vsource:
    def __init__(self, xyz=[0, 0, 0], hgrid_ie=0, nwm_fid=0, df=pd.DataFrame()):
        self.xyz = xyz
        self.usgs_st = []
        self.nwm_fid = nwm_fid

        self.upstream_path = search_path()
        self.downstream_path = search_path()
        self.hgrid_ie = hgrid_ie
        self.df = df  # time series from vsource.th

    @property
    def df(self):
        return self._df

    @df.setter
    def df(self, dataframe):
        if isinstance(dataframe, pd.DataFrame):
            if 'datetime' in dataframe.keys():
                dataframe = dataframe.set_index(pd.DatetimeIndex(pd.to_datetime(dataframe["datetime"], utc=True)))
                self._df = dataframe
        else:
            raise Exception("needs an instance of the 'TimeHistory' class.")


class search_path:
    '''
    This class is to search along the upstream or downstream path of a vsource.
    '''
    def __init__(self):
        self._segID = []
        self._order = []

    def forward_path(self, seg_id, seg_order):
        self._segID.append(seg_id)
        self._order.append(seg_order)

    def reset_path(self):
        self._segID = []
        self._order = []

    def writer(self, fname, shp_gpd=None):
        with open(fname, 'w') as fout:
            fout.write(f'lon lat \n')
            for fid in self._segID:
                fout.write(str(shp_gpd.loc[fid].lon) + " " + str(shp_gpd.loc[fid].lat) + "\n")


def find_usgs_along_nwm(iupstream=True, starting_seg_id=None, vsource=None, total_seg_length=0.0, order=0, nwm_shp=None):
    """This is a recursive function to find the up/down-stream USGS stations,
    The search stops when the total length of the path exceeds 30 km.
    The order of the stream is also recorded, but not used as a stopping criterion at this moment."""

    if starting_seg_id is None:
        return # some starting segs may not be in nwm_shp, e.g., near the boundary of the subset domain

    try:
        starting_seg_rec = nwm_shp.loc[starting_seg_id]
    except KeyError:
        return  # some starting segs may not be in nwm_shp, e.g., near the boundary of the subset domain

    if iupstream:
        vsource.upstream_path.forward_path(seg_id=starting_seg_id, seg_order=order)
    else:
        vsource.upstream_path.forward_path(seg_id=starting_seg_id, seg_order=order)

    total_seg_length += nwm_shp.loc[starting_seg_id].Length
    if (total_seg_length > 30e3):  # 30 km
        return

    if nwm_shp.loc[starting_seg_id].gages is not None and nwm_shp.loc[starting_seg_id].gages != '':  # gage found
        st = usgs_st(st_id=nwm_shp.loc[starting_seg_id].gages, nearby_nwm_fid=starting_seg_id)
        vsource.usgs_st.append(st)
        return
    else:  # search upstream or downstream
        if iupstream:
            uids = nwm_shp.loc[starting_seg_id]['from']
        else:
            uids = [nwm_shp.loc[starting_seg_id]['to']]  # 'to' is a single value, make it a list for the for loop

        for uid in uids:
            find_usgs_along_nwm(iupstream=iupstream, starting_seg_id=uid,
                                vsource=vsource, total_seg_length=total_seg_length, order=order+1, nwm_shp=nwm_shp)


def source_nwm2usgs(
    start_time_str="2017-12-01 00:00:00",
    f_shapefile="/sciclone/schism10/feiye/schism20/REPO/NWM/Shapefiles/ecgc/ecgc.shp",
    run_dir='/sciclone/schism10/feiye/Requests/RUN02a_JZ/src/NWM/',
    nwm_data_dir='/sciclone/schism10/whuang07/schism20/NWM_v2.1/',
    output_dir='/sciclone/schism10/feiye/Requests/RUN02a_JZ/src/NWM/USGS_adjusted_sources/',
):
    '''
    This script is to adjust the vsource.th file based on USGS data.
    The adjustment is done by comparing the NWM data near the vsource and the USGS data near the vsource.
    The USGS data is downloaded from USGS website using the API.
    The NWM data is pre-downloaded during the original NWM to source/sink conversion by gen_source_sink.py.
    See the Atlas paper draft for details of how USGS stations are selected for each vsource.
    In short, the USGS stations are selected along the upstream and downstream branches of the vsource.
    The NWM shapefile provides the information of the upstream and downstream branches,
    as well as the location of the USGS stations.
    However, only NWM v1's shapefile has the USGS station information.
    This is acceptable because the segment IDs are the same between NWM v1 and later versions.
    '''

    usgs_data_fname = f'{output_dir}/usgs_data.pkl'

    ''' test
    my_obs = ObsData(usgs_data_fname, from_raw_data=False)  # testing pickle
    ids = [station.id for station in my_obs.stations]
    with open(f'{run_dir}/ele_fid_dict') as json_file:
        mysrc_nwm_fid = json.load(json_file)
    
    with open(f'/sciclone/schism10/lcui01/schism20/ICOGS/ICOGS3D/Scripts/RUN24/NWM/sources.json') as json_file:
        mysrc_nwm_fid = json.load(json_file)    
    n = 0
    for key, value in mysrc_nwm_fid.items():
        if len(value) > 1:
            n += 1

        mysrc_nwm_fid[key] = int(value[0])
    '''
    

    # -----------------------------------------------------------------------------------
    # ------Read and process NWM shapefile------
    # -----------------------------------------------------------------------------------
    t1 = time.time()

    if os.path.exists(f'{f_shapefile}.pkl'):  # cached
        with open(f'{f_shapefile}.pkl', 'rb') as file:
            nwm_shp = pickle.load(file)
    else:
        nwm_shp = gpd.read_file(f_shapefile)

        # append "from" property to each seg
        nwm_shp['from'] = None

        index_g2l = {}  # global to local index
        for i, rec in nwm_shp.iterrows():
            index_g2l[rec.featureID] = i

        from_list = [[] for _ in range(len(nwm_shp))]  # directly use nwm_shp.loc is slow
        for i, rec in nwm_shp.iterrows():
            try:
                from_list[index_g2l[rec.to]].append(rec.featureID)  # index is featureID
            except KeyError:
                pass  # index_g2l[rec.to] can be None, then 0 is pre-assigned
        nwm_shp['from'] = from_list

        # set featureID as index
        nwm_shp = nwm_shp.set_index('featureID')
        # duplicate featureID as a column
        nwm_shp['featureID'] = nwm_shp.index

        with open(f'{f_shapefile}.pkl', 'wb') as file:
            pickle.dump(nwm_shp, file)

    print(f'---------------loading and processing shapefile took: {time.time()-t1} s ---------------\n')

    # -----------------------------------------------------------------------------------
    # ------load vsource------
    # -----------------------------------------------------------------------------------

    # read hgrid to determine vsource locations
    t1 = time.time()
    hg = schism_grid(f'{run_dir}/hgrid.gr3')  # hg.save()
    hg.compute_ctr()

    # read featureID at vsource elements
    if os.path.exists(f'{run_dir}/ele_fid_dict.json'):
        with open(f'{run_dir}/ele_fid_dict.json') as json_file:
            mysrc_nwm_fid = json.load(json_file)
    elif os.path.exists(f'{run_dir}/sources.json'):
        with open(f'{run_dir}/sources.json') as json_file:
            mysrc_nwm_fid = json.load(json_file)    
        for key, value in mysrc_nwm_fid.items():
            mysrc_nwm_fid[key] = int(value[0])
    else:  
        raise Exception('no ele_fid_dict or sources.json found')

    source_eles = Prop(f'{run_dir}/source_sink.in', 2).ip_group[0]
    my_th = TimeHistory(f'{run_dir}/vsource.th', start_time_str=start_time_str, mask_val=None)
    end_time_str = my_th.df.iloc[-1, 0].strftime("%m/%d/%Y, %H:%M:%S")

    my_sources = np.empty((len(source_eles)), dtype=vsource)
    for i, source_ele in enumerate(source_eles):
        my_sources[i] = vsource(
            xyz=[hg.xctr[source_ele-1], hg.yctr[source_ele-1], hg.dpe[source_ele-1]],
            hgrid_ie=source_ele,
            nwm_fid=mysrc_nwm_fid[str(source_ele)],
            df=pd.DataFrame({'datetime': my_th.datetime, 'Data': my_th.data[:, i]}),
        )
        pass
    print(f'---------------loading vsource took: {time.time()-t1} s ---------------\n')

    # -----------------------------------------------------------------------------------
    # ------search USGS stations for each vsource------
    # -----------------------------------------------------------------------------------
    t1 = time.time()
    for i, my_source in enumerate(my_sources):
        # recursively go through upstream and downstream branches to find available USGS stations
        find_usgs_along_nwm(iupstream=True, starting_seg_id=my_source.nwm_fid, vsource=my_source, total_seg_length=0.0, order=0, nwm_shp=nwm_shp)
        my_source.upstream_path.writer(f'{output_dir}/Ele{str(my_source.hgrid_ie)}_upstream', shp_gpd=nwm_shp)

        find_usgs_along_nwm(iupstream=False, starting_seg_id=my_source.nwm_fid, vsource=my_source, total_seg_length=0.0, order=0, nwm_shp=nwm_shp)
        my_source.downstream_path.writer(f'{output_dir}/Ele{str(my_source.hgrid_ie)}_downstream', shp_gpd=nwm_shp)

    usgs_stations = [st.st_id for my_source in my_sources for st in my_source.usgs_st]
    usgs_stations = list(set(usgs_stations))  # remove duplicates
    print(f'---------------searching usgs stations took: {time.time()-t1} s ---------------\n')

    # -----------------------------------------------------------------------------------
    # Download USGS data
    # -----------------------------------------------------------------------------------
    t1 = time.time()
    # pad one day before and after the start and end time
    padded_start_time = pd.to_datetime(start_time_str) - pd.Timedelta('1 day')
    padded_end_time = pd.to_datetime(end_time_str) + pd.Timedelta('1 day')
    if os.path.exists(usgs_data_fname):
        my_obs = ObsData(usgs_data_fname, from_raw_data=False)
    else:
        usgs_data = download_stations(
            param_id=usgs_var_dict['streamflow']['id'],
            station_ids=usgs_stations,
            datelist=pd.date_range(start=padded_start_time, end=padded_end_time),
        )
        my_obs = convert_to_ObsData(usgs_data, cache_fname=usgs_data_fname)

    print(f'---------------collecting usgs data took: {time.time()-t1} s ---------------\n')

    # -----------------------------------------------------------------------------------
    #  reading NWM data
    # -----------------------------------------------------------------------------------
    t1 = time.time()
    # assemble subset of fIDs
    fID_subset = []
    for my_source in my_sources:
        for st in my_source.usgs_st:
            fID_subset.append(st.nearby_nwm_fid)

    # read NWM data
    my_date_range = pd.date_range(start_time_str, end_time_str, freq="60min")
    stream_flow = read_nwm_data(
        data_dir=nwm_data_dir,
        fIDs=fID_subset,
        date_range=my_date_range
    )
    # build df
    nwm_df = pd.DataFrame({'datetime': my_date_range}).join(pd.DataFrame(stream_flow, columns=fID_subset))
    nwm_df = nwm_df.loc[:, ~nwm_df.columns.duplicated()]  # drop duplicate columns
    nwm_df = nwm_df.set_index('datetime').resample("H").mean()

    print(f'---------------reading nwm data took: {time.time()-t1} s ---------------\n')

    # -----------------------------------------------------------------------------------
    # plot
    # -----------------------------------------------------------------------------------
    linestyle1 = ['-b', '-r', '-g', '-c', '-m', '-y', '-k']
    linestyle2 = ['-.b', '-.r', '-.g', '-.c', '-.m', '-.y', '-.k']
    for i, my_source in enumerate(my_sources):

        if len(my_source.usgs_st) == 0:  # no usgs to correct nwm
            print(f'warning: source {i}_{my_source.hgrid_ie} does not have nearby usgs obs')
            continue

        # close figure if there is one and make a new one
        plt.close('all')
        plt.figure()

        # original vsource
        df = my_sources[i].df[start_time_str:end_time_str]  # .resample("D").mean()
        nrec = len(df.index)

        # NWM, USGS, adjustment
        # list usgs stations and nwm segs associated with my_source
        mysrc_usgs_id = []
        mysrc_usgs_idx = []
        mysrc_nwm_fid = []
        for st in my_source.usgs_st:
            if st.st_id in my_obs.fID2id.keys():
                idx = my_obs.fID2id[st.st_id]
                tmp_df = my_obs.stations[idx].df[start_time_str:end_time_str]
                if len(tmp_df.index) == 0:
                    print(f'warning: USGS {st.st_id} do not have any data in the selected period {start_time_str} {end_time_str}')
                else:
                    nrec_valid_usgs = len(tmp_df.resample("H").mean())
                    if nrec_valid_usgs == nrec:  # proceed only if there is enough usgs data
                        if st.st_id not in mysrc_usgs_id:  # skip duplicates
                            mysrc_usgs_id.append(st.st_id)
                            mysrc_usgs_idx.append(idx)
                            mysrc_nwm_fid.append(st.nearby_nwm_fid)
                    else:
                        print(f'warning: USGS {st.st_id} do not have enough data in the selected period {start_time_str} {end_time_str}')
            else:
                print(f'warning: USGS Station {st.st_id} is not in USGS database')
        if len(mysrc_usgs_id) == 0:  # no usgs to correct nwm
            print(f'warning: source {i}_{my_source.hgrid_ie} does not have nearby usgs obs')
            continue

        # get data from pairs of usgs and nwm
        this_vs_usgs_array = np.array(
            [
                my_obs.stations[idx].df[start_time_str:end_time_str].iloc[:, 1]
                .resample("H").mean().to_numpy()*0.028316847
                for idx in mysrc_usgs_idx
            ]
        ).T
        this_vsource_nwm_array = nwm_df[mysrc_nwm_fid].to_numpy()

        plt.plot(df.index, df['Data'], '.k', label="original vsource")
        for k, id in enumerate(mysrc_usgs_id):
            plt.plot(df.index, this_vs_usgs_array[:, k], linestyle1[k], label=f'usgs {id}')
            plt.plot(df.index, this_vsource_nwm_array[:, k], linestyle2[k], label=f'NWM near usgs {id}')

        # take the diff (including all relavant usgs-nwm pairs found)
        diff_usgs_nwm = this_vs_usgs_array - this_vsource_nwm_array
        # use maximum to prevent negative sources
        adjusted_source = np.maximum(df['Data'].to_numpy() + np.sum(diff_usgs_nwm, axis=1), 0.0)
        plt.plot(df.index, adjusted_source, '+k', label="vsource adjusted by USGS-NWM-diff")

        plt.legend()
        # plt.show()
        plt.savefig(f'{output_dir}/{i}_{my_source.hgrid_ie}.png', dpi=400)
        plt.clf()

        # replace vsource at this source
        my_th.df.iloc[:, i] = adjusted_source

        print(f'---------------successfully adjusted vsource {i}_{my_source.hgrid_ie} ---------------\n')

    my_th.df_propagate()
    my_th.writer(f'{output_dir}/adjusted_vsource.th')


if __name__ == "__main__":
    # sample usage (the final product is "adjusted_vsource.th" in the output_dir):
    source_nwm2usgs(
        start_time_str="2017-12-01 00:00:00",
        f_shapefile="/sciclone/schism10/feiye/schism20/REPO/NWM/Shapefiles/ecgc/ecgc.shp",
        run_dir='/sciclone/schism10/feiye/Requests/RUN02a_JZ/src/NWM/',
        nwm_data_dir='/sciclone/schism10/whuang07/schism20/NWM_v2.1/',
        output_dir='/sciclone/schism10/feiye/Requests/RUN02a_JZ/src/NWM/USGS_adjusted_sources/',
    )
