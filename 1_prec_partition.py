#####################
# Using T, P of winter season, finds SFF and RFF for whole grid
#####################

import xarray as xr
import numpy as np
import pandas as pd
import geopandas as gpd
import os

# precipitation_partition from VIC model default src:https://www.engr.scu.edu/~emaurer/chile/vic_taller/01_vic_training_overview_processes.pdf
t_rain_min = 0
t_snow_max = 2

def prec_part(T, P):
        
    sff = 0
    rff = 0
    sf = 0
    rf = 0
    
    if P > 0:
        if T <= t_rain_min:
            sff = 1   
        elif T >= t_snow_max:
            sff = 0   
        else:
            # linear interpolation y=-0.5t+1, slope=-0.5, intercept=1 
            sff = -0.5*T + 1
            
        # min and max cap
        if sff > 1: sff = 1
        if sff < 0: sff = 0
        
        rff = 1-sff

        sf = sff*P
        rf = rff*P
        
    # return float(sff), float(rff) 
    return float(sf), float(rf), float(sff), float(rff) 

pp_vecfunc = np.vectorize(prec_part)

# SFF
# WINTER
# nc file with DJF (winter) daily temperature data
tds = xr.open_dataset('tmean_winter_2000_2024.nc')
pds = xr.open_dataset('pmean_winter_2000_2024.nc')

temp = xr.DataArray(tds['tg'])
prec = xr.DataArray(pds['rr'])
print('Temp shape', temp.shape)
print('Prec shape', prec.shape)

sf, rf, sff, rff = pp_vecfunc(temp, prec)
# sf = snowfall_approximation(prec, temp)
print('Snow fall fraction shape', sff.shape, sf.shape)
print('Rain fall fraction shape', rff.shape, rf.shape)

# get snow percent for each day on all grid cells
sff_da = xr.DataArray(
    sff,
    dims=("time", "latitude", "longitude"),
    coords={"time": tds['time'], "latitude": tds['latitude'], "longitude": tds['longitude']},
)

# # Save the DataArray to a NetCDF file
sff_da.to_netcdf('sff_winter_2000_2024.nc')

# get snow percent for each day on all grid cells
rff_da = xr.DataArray(
    rff,
    dims=("time", "latitude", "longitude"),
    coords={"time": tds['time'], "latitude": tds['latitude'], "longitude": tds['longitude']},
)

# Save the DataArray to a NetCDF file
rff_da.to_netcdf('rff_winter_2000_2024.nc')

# get snow percent for each day on all grid cells
sf_da = xr.DataArray(
    sf,
    dims=("time", "latitude", "longitude"),
    coords={"time": tds['time'], "latitude": tds['latitude'], "longitude": tds['longitude']},
)

# Save the DataArray to a NetCDF file
sf_da.to_netcdf('sf_winter_2000_2024.nc')

# get snow percent for each day on all grid cells
rf_da = xr.DataArray(
    rf,
    dims=("time", "latitude", "longitude"),
    coords={"time": tds['time'], "latitude": tds['latitude'], "longitude": tds['longitude']},
)

# Save the DataArray to a NetCDF file
rf_da.to_netcdf('rf_winter_2000_2024.nc')
