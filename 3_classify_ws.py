import pandas as pd
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
import cartopy.crs as ccrs
from matplotlib.colors import Normalize

##### read sff file formed after 2_filter_ws.py - sff_winter_2k.csv
# sff = pd.read_csv('sf/sff_winter_2k.csv')
# sff = sff.drop(columns=['Mean', 'Class'])

# ####### find average over 2000-2024
# sff["Mean"] = sff.iloc[:, 1:].mean(axis=1)
# # sff.to_csv('sff_winter_2k.csv', index=False)

# ####### classificaiton based on https://www.nature.com/articles/nclimate2246
# # if sff > 0.15, SD 
# # if sff > 0 and <0.15, SRD, remaining RD
# sff = sff[sff['Mean']>=0]
# sff["Class"] = sff["Mean"].apply(lambda x: "SD" if x >= 0.15 else "SRD")
# sff.loc[sff["Mean"] == 0, "Class"] = "RD"
# sff = sff.dropna(subset=['Mean', 'Class'])
# sff.to_csv('sff_winter_2k.csv', index=False)

####### tabulate 
#     1950-2024   2000-2024 # 0.5
# SD    679           606
# SRD   311           302
# RD    12            94

#     1950-2024   2000-2024 # 0.15
# SD    679           606
# SRD   311           302
# RD    12            94
# class_counts = sff["Class"].value_counts()
# print(class_counts)


####### plot classification with sff in colorbar, SD in dots, SRD in dot-dash
sff1 = pd.read_csv('sf/sff_winter_19s.csv')
print(sff1["Class"].value_counts())
sff2 = pd.read_csv('sf/sff_winter_2k.csv')
print(sff2["Class"].value_counts())

min_lat = 35
max_lat = 71
min_lon = -10
max_lon = 48
clon = (min_lon + max_lon) / 2
clat = (min_lat + max_lat) / 2
world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
europe_bbox = [min_lon, max_lon, min_lat, max_lat]

fig, axs = plt.subplots(1,2, figsize=(8,8), subplot_kw={'projection': ccrs.Orthographic(central_latitude=clat, central_longitude=clon)})
cax = fig.add_axes([0.3, 0.06, 0.4, 0.03])  # [left, bottom, width, height] 

# ax1 = axs
ax1 = axs[0]
ax2 = axs[1]
ax1.set_extent(europe_bbox, crs=ccrs.PlateCarree())
ax2.set_extent(europe_bbox, crs=ccrs.PlateCarree())
ax1.set_title('1950-2024', fontsize=12)
ax2.set_title('2000-2024', fontsize=12)

world.boundary.plot(ax=ax1, linewidth=0.5, edgecolor='black', transform=ccrs.PlateCarree())
world.boundary.plot(ax=ax2, linewidth=0.5, edgecolor='black', transform=ccrs.PlateCarree())

cmap = plt.cm.jet
       
for idx, row in sff1.iterrows():
    wid = row['WID']
    sff = row['Mean']
    classi = row['Class']
    hatch='..' if classi == 'SRD' else ''
    if classi != 'RD':
        gdf = gpd.read_file(f'masks/region_{wid}.shp')
        gdf.plot(ax=ax1, color=cmap(float(sff)), edgecolor='black', linewidth=0, hatch=hatch, transform=ccrs.PlateCarree())
  
for idx, row in sff2.iterrows():
    wid = row['WID']
    sff = row['Mean']
    classi = row['Class']
    hatch='..' if classi == 'SRD' else ''
    if classi != 'RD':
        gdf = gpd.read_file(f'masks/region_{wid}.shp')
        gdf.plot(ax=ax2, color=cmap(float(sff)), edgecolor='black', linewidth=0, hatch=hatch, transform=ccrs.PlateCarree())
     
# Set the colormap and normalization
sm = ScalarMappable(cmap=cmap, norm=Normalize(vmin=0.0, vmax=1.0))  # Updated vmax to 0.7
sm.set_array([])

# Create the horizontal colorbar
cbar = fig.colorbar(sm, cax=cax, orientation='horizontal')

# Add a label to the colorbar
cbar.ax.tick_params(labelsize=10)  # Adjust tick size

# Add optional custom text (if required)
cbar.ax.text(0.5, 1.5, 'SFF', fontsize=12, ha='center', va='center', transform=cbar.ax.transAxes)

# Adjust the figure spacing
fig.subplots_adjust(wspace=0.01, hspace=0)

# # axs.set_aspect('equal')
# # axs.gridlines()
    
ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax1.gridlines()
ax2.gridlines()
plt.show()
# plt.savefig('FIGURES/eu_ws.png', dpi=500, bbox_tight=True)