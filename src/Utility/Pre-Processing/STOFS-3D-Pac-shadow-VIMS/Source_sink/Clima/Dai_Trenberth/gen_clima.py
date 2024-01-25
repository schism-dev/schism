import numpy as np
import pandas as pd
import xarray as xr

from pylib import read_schism_hgrid

def read_station_file(station_file_name):
    with open(station_file_name) as f:
        f.readline()
        f.readline()
        id =[]
        lon=[]
        lat=[]
        rank = []
        for line in f.read().splitlines():
            id.append(int(line.split(' ')[0]))
            lon.append(float(line.split(' ')[1]))
            lat.append(float(line.split(' ')[2]))
            rank.append(int(float(line.split(' ')[3])))
    return np.array(id), np.array(lon), np.array(lat), np.array(rank)

if __name__ == '__main__':

    #year = 2024 (value is independant with year. date will be changed based on NWM's date, gen_sourcesink.py)
    date_forecast = pd.date_range('2024-01-01', '2025-02-01', freq='MS')
    river_filename = 'out_river'
    hgrid = 'hgrid.gr3'

    #From WorldRivers/Dai_Trenberth
    ds = xr.open_dataset('coastal-stns-Vol-monthly.updated-May2019.nc', decode_times=False)
    times = ds.time
    stations = ds.station.values.astype(int)
    lon = ds.lon
    lat = ds.lat
    riv_name = ds.riv_name
    stn_name = ds.stn_name
    flow = ds.FLOW.values
    
    #Fraser river is id=30, index=29
    #fl = flow[:, 29]
    dates = pd.date_range(start='1900-01-01', end = '2018-12-31', freq='m')
    month = [pd.to_datetime(date).month for date in dates]
    df = pd.DataFrame(data=flow, columns=stations, index=dates)
    df['Month'] = month
    #calculate monthly average
    month_mean = df.groupby(df['Month']).mean() 
    #print(month_mean.head())

    for i in np.arange(len(stations)):
        #check nan for the first time and fill with month_mean
        mask = month_mean.loc[:, i].isnull()
        if mask.sum() > 0:
            month_mean.loc[mask == 1, i] = 0

    #add January and Feb to the end
    month_mean.loc[13] = month_mean.loc[1]
    month_mean.loc[14] = month_mean.loc[2]

    #climatology monthly
    df2024 = pd.DataFrame(data=month_mean.values, columns=stations, index=date_forecast)
    #df_2024_daily = df2024.resample('D').asfreq()
    #df_2024_daily.interpolate('linear', inplace=True)
    df_2024_hourly = df2024.resample('H').asfreq()
    df_2024_hourly.interpolate('linear', inplace=True)

    #to daily
    # tmp = month_mean.loc[df.index.month]
    # tmp['Date'] = dates
    # tmp.set_index('Date', inplace=True)
    # out_flow = df.loc[:, 0:924].copy()
    # for i in np.arange(len(stations)):
    #     #check nan for the first time and fill with month_mean
    #     mask = out_flow.loc[:, i].isnull()
    #     if mask.sum() > 0:
    #         out_flow.loc[mask == 1, i] = tmp.loc[mask == 1, i]
    #     #check nan for the second time and fill with zero
    #     mask = out_flow.loc[:, i].isnull()
    #     if mask.sum() > 0:
    #         out_flow.loc[mask == 1, i] = 0


    # print(out_flow.head())

    #read out_river
    id, lon, lat, rank = read_station_file(river_filename)
    #print(rank)

    #select rivers
    df_clima_hourly = df_2024_hourly.loc[:, rank-1]

    #Find nearest elem
    gd = read_schism_hgrid(hgrid)
    gd.compute_ctr()
    xctr = gd.xctr
    yctr = gd.yctr
    ie_river = abs((lon + 1j*lat)[:, None] - (xctr + 1j*yctr)[None, :]).argmin(axis=1)
    df_clima_hourly_transposed = df_clima_hourly.T
    #set index start with 1
    df_clima_hourly_transposed['ie_river'] = ie_river + 1
    out_flow = df_clima_hourly_transposed.groupby('ie_river').sum()
    #print(out_flow.head())
    out_flow_transposed = out_flow.T

    #check
    # values, counts = np.unique(ie_river, return_counts=True)
    # ie_r = values[15]
    # idxs = np.where(ie_river == ie_r)[0]
    # cols = rank[idxs] - 1
    # tmp = df_2024_hourly.loc[:, cols[0]]+df_2024_hourly.loc[:, cols[1]]+df_2024_hourly.loc[:, cols[2]]
    # tmp2 = out_flow_transposed.loc[:, ie_r]
    
    
    out_flow_transposed.to_csv('Dai_Trenberth_climatology_1990-2018_hourly.csv', index_label='date')
    #out_flow_transposed.to_csv('Dai_Trenberth_climatology_1990-2018_hourly.csv', index_label='date', mode='a', header=False)
