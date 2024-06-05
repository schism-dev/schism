import os
from datetime import datetime, timedelta
import glob
import argparse
import json

import numpy as np
import pandas as pd
from netCDF4 import Dataset

from relocate_source_feeder_lean import SourceSinkIn, relocate_sources

def write_th_file(dataset, timeinterval, fname, issource=True):

    data = []
    for values, interval in zip(dataset, timeinterval):
        if issource:
            data.append(" ".join([f"{interval:G}", *[f'{x: .4f}' for x in values], '\n']))
        else:
            data.append(" ".join([f"{interval:G}", *[f'{-x: .4f}' for x in values], '\n']))

    with open(fname, 'w+') as fid:
        fid.writelines(data)

def write_mth_file(temp, salinity, fname):

    data = []
    #first record
    dt = 0.
    line = [dt]
    [line.append(x) for x in temp]
    [line.append(x) for x in salinity]
    data.append(line)

    #last record
    dt = 1555200.
    line = [dt]
    [line.append(x) for x in temp]
    [line.append(x) for x in salinity]
    data.append(line)
    newset = np.array(data)
    #data.append(" ".join([f"{dt}", *[f'{x: .1f}' for x in salinity], '\n']))

    #with open(fname, 'w+') as fid:
    #    fid.write(str(line))
    np.savetxt(fname, newset, fmt='%.1f')

def get_aggregated_features(nc_feature_id, features):

    # get unique featureIDs from features
    aggregated_features = []
    for source_feats in features:
        aggregated_features.extend(list(source_feats))

    # get the indexes of aggregated_features in nc_feature_id
    in_file=[]
    for feature in aggregated_features:
        idx=np.where(nc_feature_id == int(feature))[0]
        in_file.append(idx.item())

    in_file_2 = []
    sidx = 0
    for source_feats in features:
        eidx = sidx + len(source_feats)
        in_file_2.append(in_file[sidx:eidx])
        sidx = eidx
    return in_file_2

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
    Usage: python gen_sourcesink.py "yyyy-mm-dd" or "yyyy-mm-dd hh:mm:ss"
    '''

    #cmd line input paramters 
    argparser = argparse.ArgumentParser()
    argparser.add_argument('date', 
        type=datetime.fromisoformat, 
        help="The date format 'YYYY-MM-DD HH:MM:SS'", 
    )
    args=argparser.parse_args()
    startdate=args.date
    #startdate = datetime(2022, 3, 29, 0)

    # ------------------------- hardwired inputs for operation--------------------------
    # folder to save outputs (e.g., outdir='./')
    outdir = './'

    # folder holding outputs from the original source/sink files (e.g., outputs from hindcast's NWM/gen_source_sink.py):
    # source_sink.in, sources.json
    original_source_sink_folder = './original_source_sink'

    # mapping from original to final source locations (e.g., from relocate_map.txt, which is generated from relocate_source_feeder.py)
    relocate_map = np.loadtxt(f'./relocate_map.txt', dtype=int)

    # source scaling file, containing element id and scale value, for temporary solution to correct bias
    fname_scale = 'source_scale.txt'

    # folder where nwm*.nc saved (
    fdir = 'combine'

    layer = 'conus'
    # ------------------------- end hardwired inputs--------------------------

    eid_sources = []
    times = []
    dates = []

    # read source2fid_dict, note that the order of elements may not be the same as in source_sink.in
    source2fid_dict = json.load(open(f'./{original_source_sink_folder}/sources.json'))

    #add to the final list, use element order from source_sink.in
    old_source_sink_in = SourceSinkIn(filename = f'{original_source_sink_folder}/source_sink.in')
    eid_sources.extend(list(old_source_sink_in.ip_group[0]))

    #read nc files
    files = glob.glob(f'./{fdir}/nwm*.{layer}.nc')
    files.sort()
    print(f'file 0 is {files[0]}')
    nc_fid0 = Dataset(files[0])["feature_id"][:]

    # get the index of featureID in the netcdf file in the order of source_sink.in's source element order
    source2fid = [source2fid_dict[str(x)] for x in eid_sources]
    src_ncidxs = [[np.where(nc_fid0 == int(fid))[0].item() for fid in fids] for fids in source2fid]
    # src_idxs = get_aggregated_features(nc_fid0, source2fid_dict.values())

    sources = []
    for fname in files:
        ds = Dataset(fname)
        ncfeatureid=ds['feature_id'][:]
        if not np.all(ncfeatureid == nc_fid0):  # accommodate for product switch (should be rare)
            print(f'Indexes of feature_id are changed in  {fname}')
            src_ncidxs = [[np.where(nc_fid0 == int(fid))[0].item() for fid in fids] for fids in source2fid]
            # src_idxs=get_aggregated_features(ncfeatureid, source2fid_dict.values())
            nc_fid0 = ncfeatureid

        sources.append(streamflow_lookup(fname, src_ncidxs))

        model_time = datetime.strptime(ds.model_output_valid_time, "%Y-%m-%d_%H:%M:%S")
        dates.append(str(model_time))
        times.append((model_time - startdate).total_seconds())
        ds.close()

    #relocate sources; output: vsource.th and msources.th
    df_vsources = relocate_sources(
        old_source_sink_in=old_source_sink_in,
        old_vsource=sources,
        times=np.array(times),
        outdir=outdir,
        relocate_map=relocate_map
    )
    
    #source scaling (temporary solution to correct bias)
    with open(fname_scale) as f:
       total =  f.readline().split(' ')[0]
       print(f'Total sources need to be scaled: {total}!')
       for line in f.read().splitlines():
           scale_idx = line.split(',')[-2].strip()
           scale_value = float(line.split(',')[-1])
           print(scale_idx)
           print(scale_value)
    
           print(f'Pre-scaling is {df_vsources[scale_idx]}')
           df_vsources[scale_idx] = df_vsources[scale_idx]*scale_value
           print(f'Post-scaling is {df_vsources[scale_idx]}')
    #overwrite vsource.th with scaled values
    df_vsources.to_csv(f'{outdir}/vsource.th', index=False, header=False, sep=' ', float_format='%.2f')
