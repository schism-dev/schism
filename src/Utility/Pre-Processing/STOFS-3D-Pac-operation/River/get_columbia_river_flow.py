from datetime import datetime, timedelta
import pandas as pd

unit_conv = 0.028316846

def get_columbia_river_obs(fname, datevectors):
    df = pd.read_csv(fname, sep=',', na_values=' ')
    df= df.rename(columns={df.columns[1]: 'flow'})
    #convert PST to UTC (be careful with DST)
    ts = pd.to_datetime(pd.Series(df['timestamp'].values))
    #ts2 = ts.dt.tz_localize('US/Pacific', ambiguous='infer', nonexistent='shift_forward')
    ts3 = ts.dt.tz_convert(None)
    df.insert(1, 'date_utc', ts3)

    df.set_index('date_utc', inplace=True)
    df.drop(columns=['timestamp'], inplace=True)
    data = []
    for dt in datevectors:
        data.append(round(float(df.loc[dt]['flow'])*unit_conv, 3))
    data.append(data[-1])
    data.append(data[-1])
    return data

if __name__ == '__main__':
    #date = datetime.now() - timedelta(days=1)
    date = datetime(2023, 11, 16)
    startDT = datetime(date.year, date.month, date.day)
    endDT = startDT + timedelta(days=1)
    datevectors = pd.date_range(start=startDT.strftime('%Y-%m-%d'), end=endDT.strftime('%Y-%m-%d'))
    endDT2 = startDT + timedelta(days=3)
    datevectors2 = pd.date_range(start=startDT.strftime('%Y-%m-%d'), end=endDT2.strftime('%Y-%m-%d'))

    flow = {}
    fname_columbia = 'BON.Flow-Out.Ave.1Hour.1Hour.CBT-REV'
    #get columbia river
    flow['Columbia'] = get_columbia_river_obs(fname_columbia, datevectors)
    print(flow['Columbia'])
