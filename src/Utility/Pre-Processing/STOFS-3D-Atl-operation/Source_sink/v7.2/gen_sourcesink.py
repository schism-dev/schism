"""
Generate source and sink files for the NWM model for operational forecast.
No dependency on pyschism.
"""

from datetime import datetime
import glob
import argparse
import json

import numpy as np
from netCDF4 import Dataset


def streamflow_lookup(nwm_file, indexes, threshold=-1e-5):
    '''look up streamflow in NWM's netcdf files, sum flow under each index'''
    nc = Dataset(nwm_file)
    streamflow = nc["streamflow"][:]
    streamflow[np.where(streamflow < threshold)] = 0.0
    # change masked value to zero
    streamflow[np.where(streamflow.mask)] = 0.0
    data = []
    for indxs in indexes:
        # Note: Dataset already consideres scale factor and offset.
        data.append(np.sum(streamflow[indxs]))
    nc.close()
    return data


def main():
    '''
    Usage: python gen_sourcesink.py yyyy-mm-dd or "yyyy-mm-dd hh:mm:ss"
    '''

    argparser = argparse.ArgumentParser()
    argparser.add_argument('date', type=datetime.fromisoformat,
                           help="The date format 'YYYY-MM-DD HH:MM:SS'")
    args = argparser.parse_args()
    startdate = args.date

    # startdate = datetime(2024, 3, 5)

    # ------------------------- hardwired inputs for operation--------------------------
    working_dir = './'
    nwm_folder = './combine/'
    layer = 'conus'
    # ------------------------- end hardwired inputs--------------------------

    # read source2fid_dict, note that the order of elements must be the same as in source_sink.in
    source2fid_dict = json.load(open(f'{working_dir}/sources.json', encoding='utf-8'))

    # read nc files
    files = sorted(glob.glob(f'{nwm_folder}/nwm*.{layer}.nc'))
    print(f'file 0 is {files[0]}')
    nc_fid0 = np.array(Dataset(files[0])["feature_id"])

    # Get the index of featureID in the netcdf file in the order of sources.json's keys
    # "src_ncidxs" is similar to source2fid_dict's values,
    # but with the fids replaced by their indexes in the netcdf file
    src_ncidxs = []
    for _, fids in source2fid_dict.items():
        src_ncidxs.append([np.where(nc_fid0 == int(fid))[0].item() for fid in fids])

    times = np.zeros((len(files), 1), dtype=float)  # ntimes x 1
    sources = np.zeros((len(files), len(src_ncidxs)), dtype=float)  # ntimes x nsrc
    for it, fname in enumerate(files):
        ds = Dataset(fname)

        # accommodate for different nc_fid0 due to product switch (should be rare)
        this_nc_fid0 = np.array(ds['feature_id'])
        if not np.array_equal(nc_fid0, this_nc_fid0):
            print(f'Indexes of feature_id are changed in {fname}, regenerating src_ncidxs ...')
            nc_fid0 = this_nc_fid0  # update nc_fid0 for the rest of the files
            src_ncidxs = []  # update src_ncidxs for the rest of the files
            for _, fids in source2fid_dict.items():
                src_ncidxs.append([np.where(nc_fid0 == int(fid))[0].item() for fid in fids])

        sources[it, :] = streamflow_lookup(fname, src_ncidxs)
        model_time = datetime.strptime(ds.model_output_valid_time, "%Y-%m-%d_%H:%M:%S")
        times[it] = (model_time - startdate).total_seconds()
        ds.close()

    # idx = np.argsort(np.array([int(key) for key in source2fid_dict.keys()]))
    # write source time history file, write sources[:, idx] if sorted idx is needed
    np.savetxt(f'{working_dir}/vsource.th',
               np.c_[times, sources[:, :]], fmt='%10.4f', delimiter=' ')


if __name__ == '__main__':
    main()
