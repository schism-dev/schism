from datetime import datetime
import glob
import argparse

import numpy as np
import pandas as pd
from netCDF4 import Dataset

def read_featureID_file(filename):

    with open(filename) as f:
        lines = f.readlines()
        feature_ids = []
        for line in lines:
            feature_ids.append(line.split('\n')[0])
    return feature_ids

def add_pump_to_sink(sinks, pump):

    sinks_all = []
    for row in sinks:
        #[sinks_all.extend([-i]) for i in row]
        [row.extend([i]) for i in pump.tolist()]
        sinks_all.append(row)

    return sinks_all

def write_th_file(dataset, timeinterval, fname, issource=True):

    data = []
    for values, interval in zip(dataset, timeinterval):
        if issource:
            data.append(" ".join([f"{interval:G}", *[f'{x: .4f}' for x in values], '\n']))
        else:
            data.append(" ".join([f"{interval:G}", *[f'{-x: .4f}' for x in values], '\n']))

    with open(fname, 'w+') as fid:
        fid.writelines(data)


def get_aggregated_features(nc_feature_id, features):

    in_file=[]
    for feature in features:
        idx=np.where(nc_feature_id == int(feature))[0]
        in_file.append(idx.item())

    return in_file

def streamflow_lookup(file, indexes, threshold=-1e-5):
    nc = Dataset(file)
    streamflow = nc["streamflow"][:]
    streamflow[np.where(streamflow < threshold)] = 0.0
    #change masked value to zero
    streamflow[np.where(streamflow.mask)] = 0.0
    data = []
    for indxs in indexes:
        # Note: Dataset already consideres scale factor and offset.
        data.append(np.sum(streamflow[indxs]))
    nc.close()
    return data

if __name__ == '__main__':
    '''
    Usage: python extract2asci.py yyyy-mm-dd
    Run this script in oper_3D/NWM/ directory. Inputs are in the same directory:
        1. featureID_source.ix
        2. featureID_sink.idx
        3. ./combine/*.nc
        4. pump_sinks.txt

    '''

    #input paramters 
    argparser = argparse.ArgumentParser()
    argparser.add_argument('date', type=datetime.fromisoformat, help='input file date')
    args=argparser.parse_args()
    startdate=args.date
    #startdate = datetime(2022, 3, 29, 0)

    #1. and 2.
    fname_source = 'featureID_source.idx'
    fname_sink   = 'featureID_sink.idx'

    #3.
    # zy: files = glob.glob('./combine/*.nc')
    files = glob.glob('nwm*.nc')
    files.sort()

    #4.
    fname_pump = 'pump_sinks.txt'

    #5.
    fname_scale = 'source_scale.txt'

    sources_fid = read_featureID_file(fname_source)
    sinks_fid   = read_featureID_file(fname_sink)

    sources = []
    sinks = []
    times = []
    nc_fid0 = Dataset(files[0])["feature_id"][:]
    src_idxs = get_aggregated_features(nc_fid0, sources_fid)
    snk_idxs = get_aggregated_features(nc_fid0, sinks_fid)

    start = datetime.now()
    for fname in files:
        ds = Dataset(fname)
        ncfeatureid=ds['feature_id'][:]
        if not np.all(ncfeatureid == nc_fid0):
            print(f'Indexes of feature_id are changed in  {fname}')
            src_idxs=get_aggregated_features(ncfeatureid, sources_fid)
            snk_idxs=get_aggregated_features(ncfeatureid, sinks_fid)
            nc_fid0 = ncfeatureid

        sources.append(streamflow_lookup(fname, src_idxs))
        sinks.append(streamflow_lookup(fname, snk_idxs))

        model_time = datetime.strptime(ds.model_output_valid_time, "%Y-%m-%d_%H:%M:%S")
        times.append((model_time - startdate).total_seconds())
        ds.close()

    sources = np.array(sources)
    #source scaling
    with open(fname_scale) as f:
       total =  f.readline().split(' ')[0]
       for line in f.read().splitlines():
           scale_idx = int(line.split(',')[-2])
           scale_value = float(line.split(',')[-1])

           print(f'Pre-scaling is {sources[:, scale_idx-1]}')
           sources[:, scale_idx-1] = sources[:, scale_idx-1]*scale_value
           print(f'Post-scaling is {sources[:, scale_idx-1]}')

    #add pump to sinks
    data = np.loadtxt(fname_pump)
    #ids = data[:, 0]
    #streamflows = data[:, 1]
    sinks = add_pump_to_sink(sinks, -data[:, 1])

    #write to file
    write_th_file(sources, times, 'vsource.th', issource=True)
    write_th_file(sinks, times, 'vsink.th', issource=False)
