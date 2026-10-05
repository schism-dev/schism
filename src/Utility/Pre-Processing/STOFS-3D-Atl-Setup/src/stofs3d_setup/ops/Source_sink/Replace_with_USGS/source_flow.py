#!/usr/bin/env python

"""This script is to adjust the vsource.th file based on USGS data."""

import os

os.environ["MPLBACKEND"] = "Agg"   # or set in your shell: export MPLBACKEND=Agg
os.environ.pop("DISPLAY", None)    # optional: avoid accidental GUI use

from pathlib import Path
import time
import json

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from schism_py_pre_post.Grid.Prop import Prop
from schism_py_pre_post.Timeseries.TimeHistory import TimeHistory
from stofs3d_setup.ops.USGS.ogc.source_sink import \
    download_stations, usgs_var_dict, convert_to_ObsData, \
    detect_data_gap
from stofs3d_setup.utils.utils import STOFS3D_ATL_STATES
from stofs3d_setup.ops.Source_sink.Replace_with_USGS.station_mappings import (
    NWM_TO_USGS_FLOW_ADJUSTMENT,
    OPTIONAL_FLOW_STATIONS,
    REQUIRED_FLOW_STATIONS,
)
from stofs3d_setup.ops.USGS.ogc.download_policy import (
    USGS_MAX_MISSING_GAP_DAYS,
)
from stofs3d_setup.ops.Source_sink.workflow_constants import (
    NWM_ECGC_SHAPEFILE,
)
from stofs3d_setup.ops.Source_sink.Replace_with_USGS.nwm_usgs_station_search import (
    Vsource,
    find_usgs_along_nwm,
    prepare_nwm_usgs_station_search,
    read_nwm_data,
)
from pylib import schism_grid


















def source_nwm2usgs(
    start_time_str="2017-12-01 00:00:00",
    states=['LA', 'MS', 'TX'],
    f_shapefile=NWM_ECGC_SHAPEFILE,
    original_ss_dir='/sciclone/schism10/feiye/Requests/RUN02a_JZ/src/NWM/',
    nwm_data_dir='/sciclone/schism10/whuang07/schism20/NWM_v2.1/',
    output_dir='/sciclone/schism10/feiye/Requests/RUN02a_JZ/src/NWM/USGS_adjusted_sources/',
    usgs_cache_folder=None,
    warning_log_path=None,
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
    NWM v1's shapefile is used to associate USGS stations with NWM segments. Later versions of NWM has this information in RouteLink_DA_pts_CONUS of the hydrofabric.
    This is acceptable because the segment IDs are the same between NWM v1 and later versions.
    And for the STOFS3D-ATL domain, the NWM segments seem to be the same through different NWM versions as of v3
    The output is adjusted vsource.th written to output_dir.
    It will also link the adjusted vsource.th to original_ss_dir/vsource.th and save the original_ss_dir to
    original_source_sink_before_USGS_adjustment
    '''

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(f'{output_dir}/Diag/', exist_ok=True)
    os.makedirs(f'{output_dir}/Pngs/', exist_ok=True)
    usgs_cache_folder = Path(
        usgs_cache_folder or Path(output_dir).resolve().parents[1] / 'USGS_cache'
    )
    os.makedirs(usgs_cache_folder, exist_ok=True)

    nwm_shp = prepare_nwm_usgs_station_search(
        states=states,
        nwm_shapefile=f_shapefile,
        usgs_cache_folder=usgs_cache_folder,
        manual_feature_station_links=NWM_TO_USGS_FLOW_ADJUSTMENT,
        station_diagnostic_path=f'{output_dir}/usgs_stations.csv',
        association_diagnostic_path=f'{output_dir}/auto_nearby_gages.txt',
    )

    idx = nwm_shp['gages'].notnull()
    print(f'{nwm_shp.loc[idx]}')

    # -----------------------------------------------------------------------------------
    # ------load vsource------
    # -----------------------------------------------------------------------------------

    # read hgrid to determine vsource locations
    t1 = time.time()
    hg = schism_grid(f'{original_ss_dir}/hgrid.gr3')  # hg.save()
    hg.compute_ctr()

    # read featureID at vsource elements
    if os.path.exists(f'{original_ss_dir}/ele_fid_dict.json'):
        with open(f'{original_ss_dir}/ele_fid_dict.json', encoding='utf-8') as json_file:
            mysrc_nwm_fid = json.load(json_file)
    elif os.path.exists(f'{original_ss_dir}/sources.json'):
        with open(f'{original_ss_dir}/sources.json', encoding='utf-8') as json_file:
            mysrc_nwm_fid = json.load(json_file)
        for key, value in mysrc_nwm_fid.items():
            mysrc_nwm_fid[key] = int(value[0])
    else:
        raise FileNotFoundError('no ele_fid_dict or sources.json found')

    source_eles = Prop(f'{original_ss_dir}/source_sink.in', 2).ip_group[0]
    my_th = TimeHistory(f'{original_ss_dir}/vsource.th', start_time_str=start_time_str, mask_val=None)
    end_time_str = my_th.df.iloc[-1, 0].strftime("%m/%d/%Y, %H:%M:%S")

    my_sources = np.empty((len(source_eles)), dtype=Vsource)
    for i, source_ele in enumerate(source_eles):
        my_sources[i] = Vsource(
            xyz=[hg.xctr[source_ele-1], hg.yctr[source_ele-1], hg.dpe[source_ele-1]],
            hgrid_ie=source_ele,
            nwm_fid=mysrc_nwm_fid[str(source_ele)],
            df=pd.DataFrame({'datetime': my_th.datetime, 'Data': my_th.data[:, i]}),
        )
    print(f'---------------loading vsource took: {time.time()-t1} s ---------------\n')

    # -----------------------------------------------------------------------------------
    # ------search USGS stations for each vsource------
    # -----------------------------------------------------------------------------------
    t1 = time.time()
    for i, my_source in enumerate(my_sources):
        # recursively go through upstream branches to find available USGS stations
        find_usgs_along_nwm(
            iupstream=True, starting_seg_id=my_source.nwm_fid,
            vsource=my_source, total_seg_length=0.0, order=0, nwm_shp=nwm_shp
        )
        my_source.upstream_path.writer(
            f'{output_dir}/Diag/Ele{str(my_source.hgrid_ie)}_upstream', shp_gpd=nwm_shp)

        # disable downstream search for now, because it may double count a USGS flow
        # when multiple upstream branches merge into one downstream branch
        # recursively go through downstream branches to find available USGS stations
        # find_usgs_along_nwm(
        #     iupstream=False, starting_seg_id=my_source.nwm_fid,
        #     vsource=my_source, total_seg_length=0.0, order=0, nwm_shp=nwm_shp
        # )
        # my_source.downstream_path.writer(f'{output_dir}/Ele{str(my_source.hgrid_ie)}_downstream', shp_gpd=nwm_shp)

    usgs_stations = [st.st_id for my_source in my_sources for st in my_source.usgs_st]
    usgs_stations = sorted(list(set(usgs_stations)))  # remove duplicates
    print(f'---------------searching usgs stations took: {time.time()-t1} s ---------------\n')

    # -----------------------------------------------------------------------------------
    # Download USGS data
    # -----------------------------------------------------------------------------------
    t1 = time.time()

    # pad one day before and after the start and end time
    padded_start_time = pd.to_datetime(start_time_str) - pd.Timedelta('1 day')
    padded_end_time = pd.to_datetime(end_time_str) + pd.Timedelta('1 day')
    param_dict = {
        "param_id": usgs_var_dict['streamflow']['id'],
        "station_ids": sorted(usgs_stations),
        "cache_fname": (
            f'{usgs_cache_folder}/usgs_stofs3d_atl_{usgs_var_dict["streamflow"]["id"]}_'
            f'{padded_start_time.strftime("%Y%m%d")}_{padded_end_time.strftime("%Y%m%d")}.pq'
        ),
        "padded_start_time": padded_start_time.strftime("%Y-%m-%d %H:%M:%S"),
        "padded_end_time": padded_end_time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    param_save = f'{output_dir}/Diag/usgs_download_params1.json'
    with open(param_save, 'w', encoding='utf-8') as f:
        json.dump(param_dict, f, indent=4)
    usgs_data = download_stations(
        param_id=param_dict["param_id"],
        station_ids=param_dict["station_ids"],
        cache_fname=param_dict["cache_fname"],
        datelist=pd.date_range(start=param_dict["padded_start_time"], end=param_dict["padded_end_time"]),
        required_station_ids=(
            REQUIRED_FLOW_STATIONS
            & set(param_dict["station_ids"])
        ),
        optional_station_ids=(
            OPTIONAL_FLOW_STATIONS & set(param_dict["station_ids"])
        ),
        validation_start=start_time_str,
        validation_end=end_time_str,
        warning_log_path=warning_log_path,
        warning_input_files=(
            "Source_sink/vsource.th",
            "Source_sink/source.nc",
        ),
        warning_forcing_name="streamflow",
    )
    # check data gap
    for i, data in enumerate(usgs_data.copy()):
        gap = detect_data_gap(
            df=data.df.set_index(pd.DatetimeIndex(pd.to_datetime(data.df["date"], utc=True))),
            start=pd.Timestamp(start_time_str, tz='UTC'),
            end=pd.Timestamp(end_time_str, tz='UTC'),
            freq=pd.Timedelta('1 hour')
        )
        if gap > pd.Timedelta(
            days=USGS_MAX_MISSING_GAP_DAYS
        ):
            print(f"warning: USGS station {i} ({data.station_info['id']}) has a data gap of {gap}, "
                  f'removing it from the list of available stations.\n')
            usgs_data.remove(data)

    my_obs = convert_to_ObsData(usgs_data)

    print(f'---------------collecting usgs data took: {time.time()-t1} s ---------------\n')

    # -----------------------------------------------------------------------------------
    #  reading NWM data
    # -----------------------------------------------------------------------------------
    t1 = time.time()
    test_pq = False
    if test_pq:
        nwm_df = pd.read_parquet("nwm.pq", engine="pyarrow")
    else:
        # assemble subset of fids
        fid_subset = []
        for my_source in my_sources:
            for st in my_source.usgs_st:
                fid_subset.append(st.nearby_nwm_fid)
        my_date_range = pd.date_range(start_time_str, end_time_str, freq="3H")
        stream_flow = read_nwm_data(
            data_dir=nwm_data_dir,
            fids=fid_subset,
            date_range=my_date_range
        )
        # build df
        nwm_df = pd.DataFrame({'datetime': my_date_range}).join(pd.DataFrame(stream_flow, columns=fid_subset))
        # assign utc time zone
        nwm_df['datetime'] = nwm_df['datetime'].dt.tz_localize('UTC')
        nwm_df = nwm_df.loc[:, ~nwm_df.columns.duplicated()]  # drop duplicate columns
        # nwm_df = nwm_df.set_index('datetime').resample("H").mean()
        # nwm_df = nwm_df.interpolate(method='time')  # fill in missing data

        # set datetime as index
        nwm_df = nwm_df.set_index('datetime')

    print(f'---------------reading nwm data took: {time.time()-t1} s ---------------\n')

    # -----------------------------------------------------------------------------------
    # Run adjustment for each vsource
    # -----------------------------------------------------------------------------------
    linestyle1 = ['-b', '-g', '-c', '-m', '-y', '-k']
    linestyle2 = ['-.b', '-.g', '-.c', '-.m', '-.y', '-.k']
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
                tmp_df = my_obs.stations[idx].df
                if len(tmp_df.index) == 0:
                    print(f'warning: USGS {st.st_id} do not have any data'
                          'in the selected period {start_time_str} {end_time_str}')
                else:
                    # nrec_valid_usgs = len(tmp_df.resample("H").mean())
                    # if nrec_valid_usgs == nrec:  # proceed only if there is enough usgs data
                    if tmp_df.index[0] <= pd.to_datetime(start_time_str, utc=True) and \
                       tmp_df.index[-1] >= pd.to_datetime(end_time_str, utc=True):
                        if st.st_id not in mysrc_usgs_id:  # skip duplicates already processed
                            mysrc_usgs_id.append(st.st_id)
                            mysrc_usgs_idx.append(idx)
                            mysrc_nwm_fid.append(st.nearby_nwm_fid)
                    else:
                        print(f'warning: USGS {st.st_id} do not have enough data'
                              'in the selected period {start_time_str} {end_time_str}')
            else:
                print(f'warning: USGS Station {st.st_id} is not in USGS database')
        if len(mysrc_usgs_id) == 0:  # no usgs to correct nwm
            print(f'warning: source {i}_{my_source.hgrid_ie} does not have nearby usgs obs')
            continue

        # interpolate usgs to vsource time
        this_vs_usgs_array = np.zeros((nrec, len(mysrc_usgs_id)), dtype=float)
        for k, idx in enumerate(mysrc_usgs_idx):
            this_vs_usgs_array[:, k] = np.interp(
                (df.index - df.index[0]).total_seconds().to_numpy(),
                (my_obs.stations[idx].df.index - df.index[0]).total_seconds().to_numpy(),
                my_obs.stations[idx].df.iloc[:, 1].to_numpy()*0.028316847  # convert cfs to m3/s
            )

        # this_vsource_nwm_array = nwm_df[mysrc_nwm_fid].to_numpy()
        this_vsource_nwm_array = np.zeros((nrec, len(mysrc_nwm_fid)), dtype=float)
        # interpolate nwm to vsource time
        for k, idx in enumerate(mysrc_nwm_fid):
            this_vsource_nwm_array[:, k] = np.interp(
                (df.index - df.index[0]).total_seconds().to_numpy(),
                (nwm_df.index - df.index[0]).total_seconds().to_numpy(),
                nwm_df[idx].to_numpy().squeeze()
            )

        plt.figure(figsize=(10, 6))
        plt.plot(df.index, df['Data'], '.k', linewidth=0.5, label="original vsource")
        for k, usgs_id in enumerate(mysrc_usgs_id):
            plt.plot(
                df.index, this_vs_usgs_array[:, k],
                linestyle1[k], linewidth=0.5, label=f'usgs {usgs_id}'
            )
            plt.plot(
                df.index, this_vsource_nwm_array[:, k],
                linestyle2[k], linewidth=0.5, label=f'NWM near usgs {usgs_id}'
            )

        # take the diff (including all relavant usgs-nwm pairs found)
        diff_usgs_nwm = this_vs_usgs_array - this_vsource_nwm_array
        # use maximum to prevent negative sources
        adjusted_source = np.maximum(df['Data'].to_numpy() + np.sum(diff_usgs_nwm, axis=1), 0.0)
        plt.plot(df.index, adjusted_source, '--r', label="vsource adjusted by USGS-NWM-diff")

        plt.legend()
        # plt.show()
        plt.savefig(f'{output_dir}/Pngs/{i}_{my_source.hgrid_ie}.png', dpi=400)
        plt.clf()

        # replace vsource at this source
        my_th.df.iloc[:, i+1] = adjusted_source

        print(f'-----------successfully adjusted vsource {i}_{my_source.hgrid_ie}-------------\n')

    my_th.df_propagate()
    my_th.writer(f'{output_dir}/adjusted_vsource.th')

    # link adjusted vsource.th to original source_sink dir and backup original vsource.th
    pwd = os.getcwd()
    os.system(f'cp -a {original_ss_dir} {Path(output_dir).parent}/original_source_sink_before_USGS_adjustment/')
    os.chdir(original_ss_dir)
    os.system('rm vsource.th')
    os.system(f'ln -s ../{Path(output_dir).name}/adjusted_vsource.th vsource.th')
    os.chdir(pwd)




def sample_LA_domain():
    """Sample LA domain for testing purposes."""

    # sample usage (the final product is "adjusted_vsource.th" in the output_dir):

    source_nwm2usgs(
        start_time_str="2024-03-05 00:00:00",
        states=['LA', 'MS', 'TX'],
        f_shapefile="/sciclone/schism10/Hgrid_projects/STOFS3D-v8/v20p2s2_RiverMapper/shapefiles/LA_nwm_v1p2.shp",
        original_ss_dir='/sciclone/schism10/feiye/STOFS3D-v8/I19/Source_sink/original_source_sink/',
        nwm_data_dir='/sciclone/schism10/feiye/STOFS3D-v8/I09/Source_sink/original_source_sink/20240305/',
        output_dir='/sciclone/schism10/feiye/STOFS3D-v8/I19/Source_sink/USGS_adjusted_sources/',
    )

    # source_nwm2usgs(
    #     start_time_str="2021-08-01 00:00:00",
    #     f_shapefile="/sciclone/schism10/Hgrid_projects/STOFS3D-v8/v20p2s2_RiverMapper/shapefiles/LA_nwm_v1p2.shp",
    #     original_ss_dir='/sciclone/schism10/feiye/STOFS3D-v8/I09g/Source_sink/original_source_sink/',
    #     nwm_data_dir='/sciclone/schism10/feiye/STOFS3D-v8/I09g/Source_sink/original_source_sink/20210801/',
    #     output_dir='/sciclone/schism10/feiye/STOFS3D-v8/I09j/Source_sink/USGS_adjusted_sources/',
    # )

    # source_nwm2usgs(
    #     start_time_str="2017-12-01 00:00:00",
    #     f_shapefile="/sciclone/schism10/Hgrid_projects/STOFS3D-v8/v20p2s2_RiverMapper/shapefiles/LA_nwm_v1p2.shp",
    #     original_ss_dir='/sciclone/schism10/feiye/STOFS3D-v8/I14/Source_sink/USGS_adjusted_sources/original/',
    #     nwm_data_dir='/sciclone/schism10/feiye/STOFS3D-v8/I13/Source_sink/original_source_sink/20171201/',
    #     output_dir='/sciclone/schism10/feiye/STOFS3D-v8/I14/Source_sink/USGS_adjusted_sources/',
    # )

    # source_nwm2usgs(
    #     start_time_str="2017-12-01 00:00:00",
    #     f_shapefile=NWM_ECGC_SHAPEFILE,
    #     original_ss_dir='/sciclone/schism10/feiye/STOFS3D-v8/I15_v7/Source_sink/original_source_sink/',
    #     nwm_data_dir='/sciclone/schism10/feiye/STOFS3D-v8/I13/Source_sink/original_source_sink/20171201/',
    #     output_dir='/sciclone/schism10/feiye/STOFS3D-v8/I15_v7/Source_sink/USGS_adjusted_sources/',
    # )


def sample_stofs_v7v8():
    """
    Sample STOFS v7/v8 domain for testing purposes.
    """
    source_nwm2usgs(
        start_time_str="2017-12-01 00:00:00",
        states=STOFS3D_ATL_STATES,
        f_shapefile=NWM_ECGC_SHAPEFILE,
        original_ss_dir='/sciclone/schism10/feiye/STOFS3D-v8/I15n_v7/original_source_sink0/',
        nwm_data_dir='/sciclone/schism10/feiye/STOFS3D-v8/NWM/CONUS/netcdf/CHRTOUT/for_2018_hindcast/',
        output_dir='/sciclone/schism10/feiye/STOFS3D-v8/I15n_v7/Source_sink/USGS_adjusted_sources/',
    )


if __name__ == "__main__":
    sample_LA_domain()
    print('done')
