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
# file_path = 'runoff/GRDC-Monthly.nc'  # wid to be mapped from station basins -> mrb_basins
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
# poly = 'runoff/stationbasins.geojson'   # station basins=1060, has wid to be mapped with Runoff.nc
poly = '../GRDC_Riverbasins/mrb_basins.json' # rivers=125, river basins=77 , no wid, map using river names with stationbasins
shpd = gpd.read_file(poly)
shpd = shpd.to_crs("EPSG:6933")
print(shpd.head, shpd.crs)  
# river_counts = shpd['MRBID'].nunique()
# print(river_counts)

# gdf = shpd.to_crs(epsg=3857)


# minx, miny, maxx, maxy = shpd.total_bounds  # Get bounding box of the shapefile
# print(minx, miny, maxx, maxy)


##### creating shapefile for each wids
# file_path = 'clim/p_spring_2000_2024.nc'
# ds = nc.Dataset(file_path, 'r')
# lat = ds.variables['latitude'][:]   # for masking
# lon = ds.variables['longitude'][:]

# # Step 2: Iterate through each row in the GeoDataFrame
# for index, row in shpd.iterrows():
#     # Step 3: Extract the geometry and other relevant data (e.g., ID)
#     geometry = row['geometry']
#     # region_id = row['MRBID']  # Example: Assuming the 'grdc_no' is your region identifier
#     region_id = row['grdc_no']  # Example: Assuming the 'grdc_no' is your region identifier

#     # Step 4: Create a new GeoDataFrame for the individual geometry
#     single_region = gpd.GeoDataFrame(
#         {'MRBID': [region_id], 'geometry': [geometry]},
#         crs=shpd.crs  # Retain the CRS of the original GeoDataFrame
#     )

#     # Step 5: Define the output shapefile path
#     output_shapefile = f'masks/stationbasins/region_{region_id}.shp'
    
#     # Step 6: Save the individual geometry as a shapefile
#     single_region.to_file(output_shapefile)

#     print(f"Saved shapefile for region {region_id} at {output_shapefile}")
    
#     ###### creatinf mask from shp for each wid
#     mask = gpd.read_file(output_shapefile)
#     masked_shp = regionmask.mask_geopandas(mask, lon, lat)
#     masked_shp.to_netcdf(f'masks/stationbasins/{region_id}.nc')


#     os.system(f'rm masks/stationbasins/region_{region_id}.*')
    
#     print('done', region_id)

    
##### mask out t, p, sf, rf for each wid

##### winter, spring mean for each wid

# tf = open('snow_days_basin_wise_2000_2014.csv', 'w', newline='')
# tfp = csv.writer(tf, delimiter=',')
# tfp.writerow(['MRBID', 'RIVERBASIN', 'AREA']+list(str(year) for year in range(2000, 2015)))

# for index, row in shpd.iterrows():
#     # Step 3: Extract the geometry and other relevant data (e.g., ID)
#     wid = row['MRBID'] 
#     river = row['RIVERBASIN'] 
#     area = row['geometry'].area/1e6
    
    
#     # mask out yearly avg of each poly and find its spatial mean
#     os.system(f' cdo ifnotthen ../masks/basins/{wid}.nc ../final/snow_days_2000_2014.nc ymax_{wid}.nc ')
#     os.system(f' cdo fldmean ymax_{wid}.nc ymax_st_{wid}.nc ')

#     ds = xr.open_dataset(f'ymax_st_{wid}.nc')
    
#     # for each year get mean
#     means = []
#     for i in range(14):
#         temp = ds['__xarray_dataarray_variable__'].values[i][0][0]
#         means.append(float(temp))

#     # write into file
#     tfp.writerow([wid, river, area]+ means)
#     print(river, area, means)

#     os.system(f' rm ymax_{wid}.nc ymax_st_{wid}.nc')

# tf.close()



tf = open('../final/rivers/snow_days_anom.csv', 'w', newline='')
tfp = csv.writer(tf, delimiter=',')
tfp.writerow(['MRBID', 'RIVERBASIN', 'AREA']+list(str(year) for year in range(2015, 2025)))

for index, row in shpd.iterrows():
    # Step 3: Extract the geometry and other relevant data (e.g., ID)
    wid = row['MRBID'] 
    river = row['RIVERBASIN'] 
    area = row['geometry'].area/1e6
    
    # mask out yearly avg of each poly and find its spatial mean
    os.system(f' cdo ifnotthen ../masks/basins/{wid}.nc ../final/snow_days_anom.nc ymax_{wid}.nc ')
    os.system(f' cdo fldmean ymax_{wid}.nc ymax_st_{wid}.nc ')

    ds = xr.open_dataset(f'ymax_st_{wid}.nc')
    
    # for each year get mean
    means = []
    for i in range(10):
        temp = ds['__xarray_dataarray_variable__'].values[i][0][0]
        means.append(float(temp))

    # write into file
    tfp.writerow([wid, river, area]+ means)
    print(river, area, means)

    os.system(f' rm ymax_{wid}.nc ymax_st_{wid}.nc')

tf.close()


