import geopandas as gpd
import numpy as np
import netCDF4 as nc
from datetime import date, datetime, timedelta
import pandas as pd
import matplotlib.pyplot as plt
import regionmask
import os
import csv
import xarray as xr

def convert_time(time):
    # Format as 'YYYY-MM' from days since 1700-01-01
    base_date = datetime(1700, 1, 1)
    target_date = base_date + timedelta(days=int(time))
    return target_date.strftime('%Y-%m')

###### actual runoff info file
# file_path = 'runoff/GRDC-Monthly.nc'
# ncfile = nc.Dataset(file_path, 'r')
# # print(ncfile.variables)
# wids = ncfile.variables['id'][:]
# lat = ncfile.variables['geo_y'][:]
# lon = ncfile.variables['geo_x'][:]
# D = ncfile.variables['geo_z'][:]
# A = ncfile.variables['area'][:]
# t = np.array([convert_time(idt) for idt in ncfile.variables['time'][:]])[-296:] # since 2000
# runoff = ncfile.variables['runoff_mean'][-296:, :]

# print(t, runoff.shape)

# # plot for one station
# # plt.plot(t, runoff[:,12])
# # plt.show()


# ###### watershed station info
poly = 'runoff/stationbasins.geojson'
shpd = gpd.read_file(poly)
# print(shpd.head)  

# minx, miny, maxx, maxy = shpd.total_bounds  # Get bounding box of the shapefile
# print(minx, miny, maxx, maxy)


##### creating shapefile for each wids
file_path = 'p/p_spring_2000_2024.nc'
ds = nc.Dataset(file_path, 'r')
lat = ds.variables['latitude'][:]   # for masking
lon = ds.variables['longitude'][:]

# Step 2: Iterate through each row in the GeoDataFrame
for index, row in shpd.iterrows():
    # Step 3: Extract the geometry and other relevant data (e.g., ID)
    geometry = row['geometry']
    region_id = row['grdc_no']  # Example: Assuming the 'grdc_no' is your region identifier

    # Step 4: Create a new GeoDataFrame for the individual geometry
    single_region = gpd.GeoDataFrame(
        {'grdc_no': [region_id], 'geometry': [geometry]},
        crs=shpd.crs  # Retain the CRS of the original GeoDataFrame
    )

    # Step 5: Define the output shapefile path
    output_shapefile = f'masks/region_{region_id}.shp'
    
    # Step 6: Save the individual geometry as a shapefile
    single_region.to_file(output_shapefile)

    print(f"Saved shapefile for region {region_id} at {output_shapefile}")
    
#     ###### creatinf mask from shp for each wid
#     mask = gpd.read_file(output_shapefile)
#     masked_shp = regionmask.mask_geopandas(mask, lon, lat)
#     masked_shp.to_netcdf(f'masks/{region_id}.nc')


#     os.system(f'rm masks/region_{region_id}.*')
    
#     print('done', region_id)

    
##### mask out t, p, sf, rf for each wid

##### winter, spring mean for each wid
# tf = open('sf/sff_winter_19s.csv', 'w', newline='')
# tfp = csv.writer(tf, delimiter=',')
# tfp.writerow(['WID']+list(str(year) for year in range(1950, 2025)))

# for index, row in shpd.iterrows():
#     # Step 3: Extract the geometry and other relevant data (e.g., ID)
#     wid = row['grdc_no'] 
#     print(row)
    
#     # mask out yearly avg of each poly and find its spatial mean
#     os.system(f' cdo ifnotthen masks/{wid}.nc sf/sff_winter_1950_2024.nc ymax_{wid}.nc ')
#     os.system(f' cdo fldmean ymax_{wid}.nc ymax_st_{wid}.nc ')

#     ds = xr.open_dataset(f'ymax_st_{wid}.nc')

#     # for each year get mean
#     means = []
#     for i in range(75):
#         temp = ds['__xarray_dataarray_variable__'].values[i][0][0]
#         means.append(float(temp))

#     # write into file
#     tfp.writerow([wid]+ means)

#     os.system(f' rm ymax_{wid}.nc ymax_st_{wid}.nc')

# tf.close()


