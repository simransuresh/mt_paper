import os
from datetime import datetime
import pandas as pd

dt_range = pd.date_range(start=datetime(2014, 11, 1), end=datetime(2015, 5, 31))

for dt in dt_range:
    year = dt.year
    month = dt.month
    day = dt.day

    if month in [11, 12, 1, 2, 3, 4, 5] or day in [1, 2, 3, 4, 5, 6, 7, 8, 9]:
        if len(str(month)) == 1:
            month  = '0' + str(month)
        if len(str(day)) == 1:
            day = '0' + str(day)
        
        os.system(f' wget https://www.globsnow.info/swe/nrt/{year}/data/GlobSnow_SWE_L3A_{year}{month}{day}_v.1.0.nc.gz -P ./swe/ ' )
        os.system(f' wget https://www.globsnow.info/swe/nrt/{year}/data/GlobSnow_SWE_L3A_{year}{month}{day}_v.2.0.nc.gz -P ./swe/ ' )

