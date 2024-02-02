#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 10:29:34 2023

@author: tjor
"""

import sys
import os
import numpy as np
import glob
#from monda.sorad import access, plots, qc
sys.path.append('/users/rsg-new/tjor/TaraOcean/monda/MONDA/src/monda/sorad') 

import datetime
import logging
import pandas as pd

import matplotlib 
import matplotlib.pyplot as plt
import matplotlib.cm as cm


import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm

import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

import io
import scipy.io


def plot_coveragemap_total(lat, lon , q, file_id, target, map_resolution=7):
    """ coverage map showing quality control filtered data: color scheme matches
    `results' plot"""
    if np.sum(q) > 0:
        colors= cm.cool(np.linspace(0, 1, int(sum(q))))
        
 #   breakpoint()     
    

    plt.figure(figsize=(15,10))
    extent = [-12,27,35,68]
    request = cimgt.GoogleTiles(style='satellite')
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.add_image(request, 6)
    ax.set_extent(extent, ccrs.PlateCarree())
    gl = ax.gridlines(draw_labels=True)
    gl.top_labels_top = gl.right_labels = False
    gl.xformatter =  LONGITUDE_FORMATTER
    gl.yformatter =  LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 12,  'rotation':45}
    gl.ylabel_style = {'size': 12,  'rotation': 0}
        
     
      
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.tick_params(labelsize=10)
    # plt.scatter(lon,lat,s=8,color='gray',transform=ccrs.PlateCarree(), label='Failed QC')
    lon = lon[q==1]
    lat = lat[q==1]
        
    ax.scatter(lon[0:4438],lat[0:4438],s=30,color='magenta',transform=ccrs.PlateCarree(),label='So-Rad')
    ax.scatter(lon[4439:],lat[4439:],s=30,color='yellow',transform=ccrs.PlateCarree(),label='So-Rad and HSP-1')
    plt.legend()
         # for i in range(int(sum(q))):
              #  print(i)
      #   plt.scatter(lon[q==1][i],lat[q==1][i],s=30,color='magenta',transform=ccrs.PlateCarree())
    
    plt.rc('font', size=20)
    #plt.title(str(file_id))
    # plt.legend()
    
    plt.savefig(os.path.join(target, file_id + 'total_coverage-map.png'), format='png', dpi=300)

    return



def plot_coveragemap_total_SoRad(lat, lon, q3, q4, q5, map_resolution=7):
    """ coverage map showing quality control filtered data: color scheme matches
    `results' plot"""

        
    #  breakpoint()     
    
    plt.figure(figsize=(15,10))
    extent = [-12,27,35,68]
    request = cimgt.GoogleTiles(style='satellite')
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.add_image(request, 6)
    ax.set_extent(extent, ccrs.PlateCarree())
    gl = ax.gridlines(draw_labels=True)
    gl.top_labels_top = gl.right_labels = False
    gl.xformatter =  LONGITUDE_FORMATTER
    gl.yformatter =  LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 12,  'rotation':45}
    gl.ylabel_style = {'size': 12,  'rotation': 0}
         
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.tick_params(labelsize=10)
    # plt.scatter(lon,lat,s=8,color='gray',transform=ccrs.PlateCarree(), label='Failed QC')
    
        
    ax.scatter(lon,lat,s=10,color='gray',transform=ccrs.PlateCarree(),label='So-Rad: all data')
    ax.scatter(lon[q3==1], lat[q3==1], s=10, color='green', transform=ccrs.PlateCarree(),label='So-Rad: q_3 (MONDA default)')
    ax.scatter(lon[q4==1], lat[q4==1], s=10, color='blue', transform=ccrs.PlateCarree(),label='So-Rad: q_4 (Inc. wind, tilt, rel. azimuth filter)')
    ax.scatter(lon[q5==1], lat[q5==1], s=10, color='white', transform=ccrs.PlateCarree(),label='So-Rad: q_5 (As q_4, but with atmospheric conditions filter')

    #  for i in range(int(sum(q))):
    #  print(i)
    #  plt.scatter(lon[q==1][i],lat[q==1][i],s=30,color='magenta',transform=ccrs.PlateCarree())
    
    #  plt.rc('font', size=20)
    # plt.title(str(file_id))
    # plt.legend()
    
    plt.legend(fontsize =16)
   # plt.savefig(os.path.join(target, file_id + 'total_coverage-map.png'), format='png', dpi=300)

    return



def plot_coveragemap_total_HSP(lat, lon, IDR, map_resolution=7):
    """ coverage map showing quality control filtered data: color scheme matches
    `results' plot"""

        
    #  breakpoint()
    q_HSP = ~np.isnan(IDR)      
    
    plt.figure(figsize=(15,10))
    extent = [-12,27,35,68]
    request = cimgt.GoogleTiles(style='satellite')
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.add_image(request, 6)
    ax.set_extent(extent, ccrs.PlateCarree())
    gl = ax.gridlines(draw_labels=True)
    gl.top_labels_top = gl.right_labels = False
    gl.xformatter =  LONGITUDE_FORMATTER
    gl.yformatter =  LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 12,  'rotation':45}
    gl.ylabel_style = {'size': 12,  'rotation': 0}
         
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.tick_params(labelsize=10)
    # plt.scatter(lon,lat,s=8,color='gray',transform=ccrs.PlateCarree(), label='Failed QC')
    
        
    ax.scatter(lon[q_HSP==1],lat[q_HSP==1],s=10,color='magenta',transform=ccrs.PlateCarree(),label='HSP-1: all data')

    #  for i in range(int(sum(q))):
    #  print(i)
    #  plt.scatter(lon[q==1][i],lat[q==1][i],s=30,color='magenta',transform=ccrs.PlateCarree())
    
    #  plt.rc('font', size=20)
    # plt.title(str(file_id))
    # plt.legend()
    
    plt.legend(fontsize =16)
   # plt.savefig(os.path.join(target, file_id + 'total_coverage-map.png'), format='png', dpi=300)

    return



def plot_coveragemap_total_ACS(lat, lon, map_resolution=7):
    """ coverage map showing quality control filtered data: color scheme matches
    `results' plot"""

        
    #  breakpoint()

    plt.figure(figsize=(15,10))
    extent = [-12,27,35,68]
    request = cimgt.GoogleTiles(style='satellite')
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.add_image(request, 6)
    ax.set_extent(extent, ccrs.PlateCarree())
    gl = ax.gridlines(draw_labels=True)
    gl.top_labels_top = gl.right_labels = False
    gl.xformatter =  LONGITUDE_FORMATTER
    gl.yformatter =  LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 12,  'rotation':45}
    gl.ylabel_style = {'size': 12,  'rotation': 0}
         
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.tick_params(labelsize=10)
    # plt.scatter(lon,lat,s=8,color='gray',transform=ccrs.PlateCarree(), label='Failed QC')
    
        
    ax.scatter(lon,lat,s=10,color='magenta',transform=ccrs.PlateCarree(),label='Inline underway GPS: all data')

    #  for i in range(int(sum(q))):
    #  print(i)
    #  plt.scatter(lon[q==1][i],lat[q==1][i],s=30,color='magenta',transform=ccrs.PlateCarree())
    
    #  plt.rc('font', size=20)
    # plt.title(str(file_id))
    # plt.legend()
    
    plt.legend(fontsize =16)
   # plt.savefig(os.path.join(target, file_id + 'total_coverage-map.png'), format='png', dpi=300)

    return


if __name__ == '__main__':
    
    ###########################################################################
    # '/users/rsg/tjor/TaraOcean/monda/MONDA/src/monda/tests/Lt_experiments' ##
    dir_main = '/users/rsg/tjor/TaraOcean/monda/MONDA/src/monda/tests/Tara_output_newcals/So-Rad_meta_with_TaraUnderway/So-Rad_meta_V2/'
    meta_files = sorted(glob.glob(os.path.join(dir_main, '*meta*'))) # 
        
    #
    meta_sr_list = []
    for i in range(len(meta_files)):  
            meta_sr = pd.read_csv(os.path.join(meta_files[i])) #
            meta_sr_list.append(meta_sr)
           
    #
    meta_total =  meta_sr_list[0]
    for d in meta_sr_list[1:]:
        meta_total = pd.concat([meta_total,d])
    
    lat = np.array(meta_total['lat'])
    lon = np.array(meta_total['lon'])
   
    q3 = np.array(meta_total['q_3'])
    q4 = np.array(meta_total['q_4'])
    q5 = np.array(meta_total['q_5'])
    
    IDR = np.array(meta_total['IDR_hsp'])
    
    plot_coveragemap_total_SoRad(lat, lon, q3, q4, q5, map_resolution=7)
    plot_coveragemap_total_HSP(lat, lon, IDR, map_resolution=7)
    
    file_ACS = '/users/rsg/tjor/TaraOcean/monda/MONDA/src/monda/tests/Underway_gps_lat_lon_final.mat'

    data_ACS = scipy.io.loadmat(file_ACS)
    lat_ACS = (data_ACS['data1']['lat']).item().squeeze()
    lon_ACS = (data_ACS['data1']['lon']).item().squeeze()

    plot_coveragemap_total_ACS(lat_ACS, lon_ACS, map_resolution=7)
    
    #ax = plt.axes(projection=ccrs.PlateCarree())
    #ax.coastlines()
    
    #plt.figure()
    #plt.scatter(lon,lat)
    
    #  for i in range(len(meta_total)):
    #   if np.array(meta_total['timestamp'])[q==1][i][0:7] =='2023-06':
    #       print('hi')     
    #       print(i) 
    #       break
        
     # plots.plot_coveragemap_total(lat, lon, q_4, '_', dir_main)
  #  plot_coveragemap_total_V2(lat, lon, q3, q4, q5, file_id, target, map_resolution=7)
 #   extend = [np.floor(np.min(lon*10))/10, np.ceil(np.max(lon*10))/10, np.floor(np.min(lat*10))/10, np.ceil(np.max(lat*10))/10]