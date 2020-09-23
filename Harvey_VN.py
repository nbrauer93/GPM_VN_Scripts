#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 13:19:05 2020

@author: noahbrauer
"""


import netCDF4
import gzip
import numpy as np
import matplotlib.pyplot as plt

import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import shiftgrid
import pyart



import os
import shutil
import tempfile



def open_netcdf(fname):
    if fname.endswith(".gz"):
        infile = gzip.open(fname, 'rb')
        tmp = tempfile.NamedTemporaryFile(delete=False)
        shutil.copyfileobj(infile, tmp)
        infile.close()
        tmp.close()
        data = netCDF4.Dataset(tmp.name)
        os.unlink(tmp.name)
    else:
        data = netCDF4.Dataset(fname)
    return data

file = 'GRtoDPR.KHGX.170827.19870.V06A.KU.NS.1_21.nc.gz'


nc = open_netcdf(file)


latitude = nc.variables['latitude'][:].T
longitude = nc.variables['longitude'][:].T
z = nc.variables['ZFactorCorrected'][:]
nw = nc.variables['Nw'][:]/10
dm = nc.variables['Dm'][:]

#Remove bad values

nw[nw<=0] = np.nan
dm[dm<=0] = np.nan







elev_angle = 0
z_lowest = z[elev_angle,:]
nw_lowest = nw[elev_angle,:]
dm_lowest = dm[elev_angle,:]


#Determine beam height and distance from the radar; Equations from Doviak and Zrnic

beam_height = np.zeros()





#%%

#from matplotlib.colors import ListedColormap

colormap = ['white','dodgerblue', 'deepskyblue', 'lawngreen', 'lightgreen', 'green', 'gold', 'darkorange', 'red', 'firebrick']

z_color = np.empty(667,dtype = 'str')
z_color = []

#for i in range(len(z_color)):
for i in range(len(z_lowest)):

    if z_lowest[i]<=5:
        #z_color[i] = colormap[0]
        z_color.append(colormap[0])
    
    elif z_lowest[i]>5 and z_lowest[i]<=10:
        #z_color[i] = colormap[1]
        z_color.append(colormap[1])
    
    elif z_lowest[i]>10 and z_lowest[i]<=15:
        #z_color[i] = colormap[2]
        z_color.append(colormap[2])
    
    elif z_lowest[i]>15 and z_lowest[i]<=20:
        #z_color[i] = colormap[3]
        z_color.append(colormap[3])
    
    elif z_lowest[i]>20 and z_lowest[i]<=25:
        #z_color[i] = colormap[4]
        z_color.append(colormap[4])
    
    elif z_lowest[i]>25 and z_lowest[i]<=30:
        #z_color[i] = colormap[5]
        z_color.append(colormap[5])
        
    elif z_lowest[i]>30 and z_lowest[i]<=35:
        #z_color[i] = colormap[6]
        z_color.append(colormap[6])
    
    elif z_lowest[i]>35 and z_lowest[i]<=40:
        #z_color[i] = colormap[7]
        z_color.append(colormap[7])
        
    elif z_lowest[i]>40 and z_lowest[i]<=45:
        #z_color[i] = colormap[8]
        z_color.append(colormap[8])
        
    elif z_lowest[i] == -100:
        #z_color[i] = colormap[0]
        z_color.append(colormap[0])
        
        
from matplotlib.colors import ListedColormap
cmap_z = ListedColormap(colormap)        
    

print(np.nanmax(nw))
print(np.nanmin(nw))
#Now do the same thing for log10(Nw) 


nw_color = []

for i in range(len(nw)):

    if nw_lowest[i]<=1:
        
        nw_color.append(colormap[0])
    
    elif nw_lowest[i]>1 and nw_lowest[i]<=1.5:
        
        nw_color.append(colormap[1])
    
    elif nw_lowest[i]>1.5 and nw_lowest[i]<=2:
        
        nw_color.append(colormap[2])
    
    elif nw_lowest[i]>2 and nw_lowest[i]<=2.5:
        
        nw_color.append(colormap[3])
    
    elif nw_lowest[i]>2.5 and nw_lowest[i]<=3:
        
        nw_color.append(colormap[4])
    
    elif nw_lowest[i]>3 and nw_lowest[i]<=3.5:
        
        nw_color.append(colormap[5])
        
    elif nw_lowest[i]>3.5 and nw_lowest[i]<=4:
        
        nw_color.append(colormap[6])
    
    elif nw_lowest[i]>4 and nw_lowest[i]<=4.5:
       
        nw_color.append(colormap[7])
        
    elif nw_lowest[i]>4.5 and nw_lowest[i]<=5:
        
        nw_color.append(colormap[8])
        
    
#Now for dm
        
        
print(np.nanmin(dm))
print(np.nanmax(dm))  

dm_colormap = ['white','dodgerblue', 'deepskyblue', 'lawngreen', 'lightgreen', 'green', 'yellow', 'gold', 'darkorange', 'red', 'firebrick', 'blueviolet', 'magenta']      
        
dm_color = []    

for i in range(len(dm)):

    if dm_lowest[i]<=0.5:
        
        dm_color.append(dm_colormap[0])
    
    elif dm_lowest[i]>0.5 and dm_lowest[i]<=0.75:
        
        dm_color.append(dm_colormap[1])
    
    elif dm_lowest[i]>0.75 and dm_lowest[i]<=1.0:
        
        dm_color.append(dm_colormap[2])
    
    elif dm_lowest[i]>1.0 and dm_lowest[i]<=1.25:
        
        dm_color.append(dm_colormap[3])
    
    elif dm_lowest[i]>1.25 and dm_lowest[i]<=1.5:
        
        dm_color.append(dm_colormap[4])
    
    elif dm_lowest[i]>1.5 and dm_lowest[i]<=1.75:
        
        dm_color.append(dm_colormap[5])
        
    elif dm_lowest[i]>1.75 and dm_lowest[i]<=2.0:
        
        dm_color.append(dm_colormap[6])
    
    elif dm_lowest[i]>2.0 and dm_lowest[i]<=2.25:
       
        dm_color.append(dm_colormap[7])
        
    elif dm_lowest[i]>2.25 and dm_lowest[i]<=2.5:
        
        dm_color.append(dm_colormap[8])
        
    elif dm_lowest[i]>2.5 and dm_lowest[i]<=2.75:
        
        dm_color.append(dm_colormap[9])
        
    elif dm_lowest[i]>2.75 and dm_lowest[i]<=3.0:
        
        dm_color.append(dm_colormap[10])
        
    elif dm_lowest[i]>3.0 and dm_lowest[i]<=3.25:
        
        dm_color.append(dm_colormap[11])
        
    elif dm_lowest[i]>3.25 and dm_lowest[i]<=3.5:
        
        dm_color.append(dm_colormap[12])
            
            
cmap_dm = ListedColormap(dm_colormap)           






#%%
import matplotlib

elev_angle = 0
label_size = 24
title_size = 28
tick_label_size = 20

#Setup plotting 
cmin = 0; cmax = 50; cint = 5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name=cmap_z,lut=nlevs)

colour_norm_object = matplotlib.colors.Normalize(vmin=cmin, vmax=cmax, clip=False)
scalar_mappable_object = plt.cm.ScalarMappable(cmap=cmap, norm=colour_norm_object)
scalar_mappable_object.set_array(z_color)

figure_object, axes_object = plt.subplots(1, 1, figsize=(10, 10))

c = plt.scatter(longitude[:, elev_angle], latitude[:, elev_angle], c = z_color, vmin = 0, vmax = 50, cmap = cmap, edgecolors = 'none')


for i in range(len(z_color)):
        if z_color[i] ==0: 
            c2 = plt.scatter(longitude[i,elev_angle], latitude[i,elev_angle], edgecolors = 'k', facecolors = 'none')
    
    
plt.xlabel('Longitude', size = label_size)
plt.ylabel('Latitude', size = label_size)
plt.xticks(size = tick_label_size)
plt.yticks(size = tick_label_size)
plt.title(r'KHGX $0.5^{o}$ $Z_{H}$ 8/27 2023 UTC ',name='Calibri',size=26)

color_bar_object = plt.colorbar(ax=axes_object, mappable=scalar_mappable_object, orientation='vertical')

color_bar_object.set_label('dBZ', size = title_size)
color_bar_object.ax.tick_params(labelsize = label_size)

plt.show()
plt.close(figure_object)


#Now plot log10(nw)

cmin = 1; cmax = 5; cint = 0.5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name=cmap_z,lut=nlevs)

colour_norm_object = matplotlib.colors.Normalize(vmin=cmin, vmax=cmax, clip=False)
scalar_mappable_object = plt.cm.ScalarMappable(cmap=cmap, norm=colour_norm_object)
scalar_mappable_object.set_array(z_color)

figure_object, axes_object = plt.subplots(1, 1, figsize=(10, 10))

c = plt.scatter(longitude[:, elev_angle], latitude[:, elev_angle], c = nw_color, vmin = 1, vmax = 5, cmap = cmap, edgecolors = 'none')


for i in range(len(nw_color)):
        if nw_color[i] ==0: 
            c2 = plt.scatter(longitude[i,elev_angle], latitude[i,elev_angle], edgecolors = 'k', facecolors = 'none')
    
    
plt.xlabel('Longitude', size = label_size)
plt.ylabel('Latitude', size = label_size)
plt.xticks(size = tick_label_size)
plt.yticks(size = tick_label_size)
plt.title(r'KHGX $0.5^{o}$ $log_{10}(N_{w})$ 8/27 2023 UTC ',name='Calibri',size=26)

color_bar_object = plt.colorbar(ax=axes_object, mappable=scalar_mappable_object, orientation='vertical')

color_bar_object.set_label(r'[mm $m^{-3}$]', size = title_size)
color_bar_object.ax.tick_params(labelsize = label_size)

plt.show()
plt.close(figure_object)

###Now for Dm

cmin = 0.5; cmax = 3.5; cint = 0.25; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name=cmap_dm,lut=nlevs)

colour_norm_object = matplotlib.colors.Normalize(vmin=cmin, vmax=cmax, clip=False)
scalar_mappable_object = plt.cm.ScalarMappable(cmap=cmap, norm=colour_norm_object)
scalar_mappable_object.set_array(z_color)

figure_object, axes_object = plt.subplots(1, 1, figsize=(10, 10))

c = plt.scatter(longitude[:, elev_angle], latitude[:, elev_angle], c = dm_color, vmin = 0.5, vmax = 3.5, cmap = cmap, edgecolors = 'none')


for i in range(len(dm_color)):
        if dm_color[i] ==0: 
            c2 = plt.scatter(longitude[i,elev_angle], latitude[i,elev_angle], edgecolors = 'k', facecolors = 'none')
    
    
plt.xlabel('Longitude', size = label_size)
plt.ylabel('Latitude', size = label_size)
plt.xticks(size = tick_label_size)
plt.yticks(size = tick_label_size)
plt.title(r'KHGX $0.5^{o}$ $D_{m}$ 8/27 2023 UTC ',name='Calibri',size=26)

color_bar_object = plt.colorbar(ax=axes_object, mappable=scalar_mappable_object, orientation='vertical')

color_bar_object.set_label(r'[mm]', size = title_size)
color_bar_object.ax.tick_params(labelsize = label_size)

plt.show()
plt.close(figure_object)
