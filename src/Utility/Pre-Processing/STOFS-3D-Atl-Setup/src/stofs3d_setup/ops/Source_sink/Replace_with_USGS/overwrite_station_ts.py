"""
Replace the time series of a station with a new time series.
The time stamps of the new time series should be the same as the old one.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from schism_py_pre_post.Grid.SourceSinkIn import source_sink
from schism_py_pre_post.Download.download_usgs_with_api import \
    download_stations, usgs_var_dict


def get_usgs_ts(usgs_stations, start_time_str, end_time_str):
    """
    Get the USGS time series of the stations.
    """
    padded_start_time = pd.to_datetime(start_time_str) - pd.Timedelta('1 day')
    padded_end_time = pd.to_datetime(end_time_str) + pd.Timedelta('1 day')
    usgs_data = download_stations(
        param_id=usgs_var_dict['streamflow']['id'],
        station_ids=usgs_stations,
        datelist=pd.date_range(start=padded_start_time, end=padded_end_time),
    )
    return usgs_data


def replace_station_ts_with_ts():
    """
    Replace the time series of a station with a new time series.
    """

    source_sink_dataset1 = source_sink(
        '/sciclone/schism10/feiye/STOFS3D-v8/I13x_v7/Source_sink/original_source_sink/'
    )
    source_sink_dataset2 = source_sink(
        '/sciclone/schism10/feiye/STOFS3D-v8/I13x_v7/Source_sink/'
    )

    element_id_mapping = {
        '5772778': '5773144'
    }

    if not all(
        [ele_id in source_sink_dataset1.vsource.df.columns for ele_id in element_id_mapping.keys()]
    ):
        raise ValueError("The element ids should be in the base source_sink dataset.")
    if not all(
        [ele_id in source_sink_dataset2.vsource.df.columns for ele_id in element_id_mapping.values()]
    ):
        raise ValueError("The element ids should be in the target source_sink dataset.")
    if source_sink_dataset1.vsource.data.shape[1] == source_sink_dataset2.vsource.data.shape[1]:
        raise ValueError("The number of time steps in the two datasets should be the same.")

    for base_ele_id, target_ele_id in element_id_mapping.items():
        source_sink_dataset2.vsource.df.loc[target_ele_id] = source_sink_dataset1.vsource.df.loc[base_ele_id]


def replace_station_ts_with_usgs(start_time: pd.Timestamp, source_usgs_mapping: dict):
    """
    Replace the time series of a station with the USGS time series.

    Example inputs:
        start_time = pd.to_datetime('2017-12-01 00:00:00')

        source_usgs_mapping = {
            '5773144': '07381490',
            '1234': '0123456'
        }
    """

    source_sink_dataset = source_sink(
        '/sciclone/schism10/feiye/STOFS3D-v8/I13x_v7/Source_sink/',
        start_time_str=start_time.strftime('%Y-%m-%d %H:%M:%S')
    )
    end_time = source_sink_dataset.vsource.df['datetime'].iloc[-1]

    for ele_id, usgs_station in source_usgs_mapping.items():
        usgs_data = get_usgs_ts([usgs_station], start_time, end_time)
        df = usgs_data[0].df
        df['date'] = pd.to_datetime(df['date'], utc=True)
        # convert time to seconds since start_time
        df_time = (df['date'] - pd.Timestamp(start_time, tz='UTC')).dt.total_seconds()
        plt.plot(
            df['date'], df['value'] * 0.0283168,
            'r', linewidth=3.5, label=f'USGS {usgs_station}')  # convert to m^3/s

        plt.plot(
            source_sink_dataset.vsource.df['datetime'],
            source_sink_dataset.vsource.df[ele_id], label=f'Source/Sink {ele_id}')
        plt.grid()

        interp_ts = np.interp(source_sink_dataset.vsource.time, df_time.values, df['value'].values * 0.0283168)
        plt.plot(
            source_sink_dataset.vsource.df['datetime'], interp_ts,
            '-.k', linewidth=1, label=f'Interpolated {ele_id}')

        plt.legend()
        plt.show()

        source_sink_dataset.vsource.df[ele_id] = interp_ts

    source_sink_dataset.nc_writer(
        '/sciclone/schism10/feiye/STOFS3D-v8/I13x_v7/Source_sink/Adjusted_source_sink/'
    )
    source_sink_dataset.writer(
        '/sciclone/schism10/feiye/STOFS3D-v8/I13x_v7/Source_sink/Adjusted_source_sink/'
    )
    print("The new source_sink file has been written.")


if __name__ == '__main__':
    replace_station_ts_with_usgs(pd.to_datetime('2017-12-01 00:00:00'), {'5773144': '07381490'})
    replace_station_ts_with_ts()
