'''
Subset inputs by time
'''

import xarray as xr
# import numpy as np
import pandas as pd
from pathlib import Path
from glob import glob
from pylib_experimental.schism_file import TimeHistory, SourceSink


def subset_nc(wdir, base_time: pd.Timestamp, start_time: pd.Timestamp, end_time: pd.Timestamp):
    '''
    Subset nudging files
    The nudging files to be subsetted are:
        - *nu.nc
        - *D.th.nc
    Assuming time is in days since base_time for *nu.nc
    Assuming time is in seconds since base_time for *D.th.nc
    '''
    file_dict = {}
    file_names = glob(f'{wdir}/*D.th.nc')
    for file in file_names:
        file_dict[Path(file).stem] = {
            'file_path': Path(file),
            'time_unit': 's',
        }
    file_names = glob(f'{wdir}/*nu.nc')
    for file in file_names:
        file_dict[Path(file).stem] = {
            'file_path': Path(file),
            'time_unit': 'D',
        }

    for file, file_info in file_dict.items():
        if Path(file_info['file_path']).exists() is False:
            print(f'File {file} does not exist. Skipping.')
            continue

        print('subsetting file:', file)

        ds = xr.open_dataset(file_info['file_path'])
        time_array = base_time + pd.to_timedelta(ds['time'].values, unit=file_info['time_unit'])

        subset_mask = (time_array >= start_time) & (time_array <= end_time)

        subset = ds.sel(time=subset_mask)

        # realign the time coordinate
        new_time = subset['time'].values - subset['time'].values[0]
        subset = subset.assign_coords(time=new_time)

        output_file = (
            wdir + Path(file_info['file_path']).stem +
            f'_{start_time.strftime("%Y%m%d")}_{end_time.strftime("%Y%m%d")}.nc'
        )
        subset.to_netcdf(output_file)


def subset_th(wdir, base_time: pd.Timestamp, start_time: pd.Timestamp, end_time: pd.Timestamp):
    '''
    Subset TH files
    The TH files to be subsetted are:
        - *.th
    Assuming time is in days since base_time
    '''
    file_names = glob(f'{wdir}/*.th')
    for file in file_names:
        if Path(file).exists() is False:
            print(f'File {file} does not exist. Skipping.')
            continue

        print('subsetting file:', file)
        th = TimeHistory.from_file(file, start_time_str=base_time.strftime('%Y%m%d %H%M%S'))

        time_array = base_time + pd.to_timedelta(th.time, unit='s')
        subset_mask = (time_array >= start_time) & (time_array <= end_time)
        th = th[subset_mask]

        time_slice = slice(start_time, end_time)
        th = th[time_slice]

        output_file = (
            wdir + Path(file).stem +
            f'_{start_time.strftime("%Y%m%d")}_{end_time.strftime("%Y%m%d")}.th'
        )
        th.writer(output_file)


def subset_source_nc(wdir, base_time: pd.Timestamp, start_time: pd.Timestamp, end_time: pd.Timestamp):
    '''
    Subset source/sink files
    The source/sink files to be subsetted are:
        - source.nc
    If source/sink is in ascii format, they can be handled by subset_th
    '''
    file_name = f'{wdir}/source.nc'

    if Path(file_name).exists() is False:
        print(f'File {file_name} does not exist. Skipping subsetting source.nc.')
        return

    print('subsetting file:', file_name)
    my_ss = SourceSink.from_ncfile(file_name, start_time_str=base_time.strftime('%Y-%m-%d %H:%M:%S'))

    subset_ss = my_ss.subset_by_time(start_time, end_time)
    subset_ss.nc_writer(
        f'{wdir}/source_{start_time.strftime("%Y%m%d")}_{end_time.strftime("%Y%m%d")}.nc'
    )


def subset_inputs_by_time():
    '''
    Subset inputs by time.
    The time-dependent inputs to be subsetted are:
        - *D.th.nc
        - *nu.nc
        - *.th
        - source/sink
    '''


if __name__ == '__main__':
    # Example usage
    WDIR = '/sciclone/schism10/feiye/STOFS3D-v7.3/I21i/Time_subset_inputs/'
    base_date = pd.Timestamp('2019-12-01')
    start_date = pd.Timestamp('2019-12-31')
    end_date = pd.Timestamp('2021-01-02')

    subset_nc(WDIR, base_date, start_date, end_date)
    subset_source_nc(WDIR, base_date, start_date, end_date)
    subset_th(WDIR, base_date, start_date, end_date)
