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
sys.path.append('/users/rsg-new/tjor/TaraOcean/monda/MONDA/src/monda/sorad') # standard (non-diffuse) implementation of 3C using latest 2021 version of git repo
import access, plots, qc
import monda.sorad
import datetime
import logging
import pandas as pd
import argparse
from math import ceil


import matplotlib 
import matplotlib.pyplot as plt
import matplotlib.cm as cm


import os
import sys
import logging
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm

import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

log = logging.getLogger('sorad-plotter')
myFormat = '%(asctime)s | %(name)s | %(levelname)s | %(message)s'
formatter = logging.Formatter(myFormat)
logging.basicConfig(level = 'INFO', format = myFormat, stream = sys.stdout)

if __name__ == '__main__':
    
    ###########################################################################
    # '/users/rsg/tjor/TaraOcean/monda/MONDA/src/monda/tests/Lt_experiments' ##
    dir_main = '/users/rsg/tjor/TaraOcean/monda/MONDA/src/monda/tests/Tara_output_newcals/So-Rad_meta/'
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
    q = np.array(meta_total['q_1'])



    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    
    plt.figure()
    plt.scatter(lon,lat)

    



#  for i in range(len(meta_total)):
     #   if np.array(meta_total['timestamp'])[q==1][i][0:7] =='2023-06':
     #       print('hi')     
     #       print(i) 
     #       break
        
    plots.plot_coveragemap_total(np.array(meta_total['lat']), np.array(meta_total['lon']), np.array(meta_total['q_1']), '_', dir_main)
   
 #   extend = [np.floor(np.min(lon*10))/10, np.ceil(np.max(lon*10))/10, np.floor(np.min(lat*10))/10, np.ceil(np.max(lat*10))/10]