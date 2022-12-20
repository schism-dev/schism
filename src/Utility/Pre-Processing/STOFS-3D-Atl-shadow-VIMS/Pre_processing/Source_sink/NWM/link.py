from datetime import datetime, timedelta
import os

import numpy as np

if __name__ == "__main__":
    '''
    This script is used to link VIMS archived NWM
    '''
    startdate = datetime(2022, 1, 1)
    rnday = timedelta(days=171)
    timevector = np.arange(startdate, rnday, timedelta(hours=3)).astype(datetime)
    print(timevector)

    basepath = '/sciclone/scr10/lcui01/ICOGS3D/RUN_20220101/NWM'

    for i, date in enumerate(timevector):
        if date.hour == 0:
            date2 = date - timedelta(days=1)
            src = f'{basepath}/{date2.strftime("%Y%m%d")}/nwm.t00z.medium_range.channel_rt_1.f024.conus.nc'
        else:
            src = f'{basepath}/{date.strftime("%Y%m%d")}/nwm.t00z.medium_range.channel_rt_1.f{int(date.hour):03d}.conus.nc'
        dst = f'nwm.t00z.medium_range.channel_rt_1.{date.strftime("%Y%m%d%H")}.conus.nc'
        print(src)
        print(dst)
        os.symlink(src, dst)
