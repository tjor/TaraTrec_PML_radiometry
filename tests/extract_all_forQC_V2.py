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

from statistics import mode

import os
import sys
import logging
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm

import scipy
from scipy import stats

import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

log = logging.getLogger('sorad-plotter')
myFormat = '%(asctime)s | %(name)s | %(levelname)s | %(message)s'
formatter = logging.Formatter(myFormat)
logging.basicConfig(level = 'INFO', format = myFormat, stream = sys.stdout)


def plot_rrs_qc_tilthresholds(rrs, time, wl, q_3, t_1, t_2, t_3, file_id, target):
    """ rrs plot function showing sequential quality control filters"""

    if np.sum(q_3) > 0:
        plt.figure()
        plt.figure(figsize=(12, 8))
        plt.suptitle(str(file_id))

        ymax = np.ceil(np.nanmax(rrs.T[:, q_3==1]*1000))/1000;

        plt.subplot(2, 2, 1)
        plt.title(f"Default MONDA QC (n = {int(np.sum(q_3))})")
        plt.plot(wl, rrs.T[:,t_1==1], linewidth=0.4, alpha=0.6)
        plt.xlim(350, 800)
        plt.ylim(-0.0005, ymax) # force axis limits
        plt.grid()
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('$R_{rs}$ [sr$^{-1}$]')

        plt.subplot(2, 2, 2)
        plt.title(f" Tilt > 5 degs removed  (n = {int(np.sum(t_1))})")
        plt.plot(wl, rrs.T[:,t_1==1], linewidth=0.4, alpha=0.8)
        plt.xlim(350, 800)
        plt.ylim(-0.0005, ymax) # force axis limits
        plt.grid()
        plt.xlabel('Wavelength [nm]')


        plt.subplot(2, 2, 3)
        plt.title(f"Tilt > 3 degs removed  (n = {int(np.sum(t_2))})")
        plt.plot(wl, rrs.T[:,t_2==1], linewidth=0.4, alpha=0.8)
        plt.xlim(380,800 )
        plt.ylim(-0.0005, ymax) # force axis limits
        plt.grid()
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('$R_{rs}$ [sr$^{-1}$]')
        
        plt.subplot(2, 2, 4)
        plt.title(f"Tilt > 2 degs removed (n = {int(np.sum(t_3))})")
        plt.plot(wl, rrs.T[:,t_3==1], linewidth=0.4, alpha=0.8)
        plt.xlim(380,800)
        plt.ylim(-0.0005, ymax) # force axis limits
        plt.grid()
        plt.xlabel('Wavelength [nm]')

        plt.subplots_adjust(hspace=0.5)

        plt.savefig(os.path.join(target, file_id + '_tilt_QC.png'), format='png', dpi=200)
 
        return

def plot_rrs_qc_tilstdthresholds(rrs, time, wl, q_3, t_1, t_2, t_3, file_id, target):
    """ rrs plot function showing sequential quality control filters"""

    if np.sum(q_3) > 0:
        plt.figure()
        plt.figure(figsize=(12, 8))
        plt.suptitle(str(file_id))

        ymax = np.ceil(np.nanmax(rrs.T[:, q_3==1]*1000))/1000;

        plt.subplot(2, 2, 1)
        plt.title(f"Default MONDA QC (n = {int(np.sum(q_3))})")
        plt.plot(wl, rrs.T[:,t_1==1], linewidth=0.4, alpha=0.6)
        plt.xlim(350, 800)
        plt.ylim(-0.0005, ymax) # force axis limits
        plt.grid()
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('$R_{rs}$ [sr$^{-1}$]')

        plt.subplot(2, 2, 2)
        plt.title(f" Tilt STD > 3 degs removed  (n = {int(np.sum(t_1))})")
        plt.plot(wl, rrs.T[:,t_1==1], linewidth=0.4, alpha=0.8)
        plt.xlim(350, 800)
        plt.ylim(-0.0005, ymax) # force axis limits
        plt.grid()
        plt.xlabel('Wavelength [nm]')


        plt.subplot(2, 2, 3)
        plt.title(f"Tilt STD > 2 degs removed  (n = {int(np.sum(t_2))})")
        plt.plot(wl, rrs.T[:,t_2==1], linewidth=0.4, alpha=0.8)
        plt.xlim(380,800 )
        plt.ylim(-0.0005, ymax) # force axis limits
        plt.grid()
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('$R_{rs}$ [sr$^{-1}$]')
        
        plt.subplot(2, 2, 4)
        plt.title(f"Tilt STD > 1 degs removed (n = {int(np.sum(t_3))})")
        plt.plot(wl, rrs.T[:,t_3==1], linewidth=0.4, alpha=0.8)
        plt.xlim(380,800)
        plt.ylim(-0.0005, ymax) # force axis limits
        plt.grid()
        plt.xlabel('Wavelength [nm]')

        plt.subplots_adjust(hspace=0.5)

        plt.savefig(os.path.join(target, file_id + '_tiltSTD_QC.png'), format='png', dpi=200)
        
        plt.tight_layout()
 
        return



def plot_rrs_qc_windthresholds(rrs, time, wl, q_3, t_1, t_2, t_3, file_id, target):
    """ rrs plot function showing sequential quality control filters"""

    if np.sum(q_3) > 0:
        plt.figure()
        plt.figure(figsize=(12, 8))
        plt.suptitle(str(file_id))

        ymax = np.ceil(np.nanmax(rrs.T[:, q_3==1]*1000))/1000;

        plt.subplot(2, 2, 1)
        plt.title(f"Default MONDA QC (n = {int(np.sum(q_3))})")
        plt.plot(wl, rrs.T[:,t_1==1], linewidth=0.4, alpha=0.6)
        plt.xlim(350, 800)
        plt.ylim(-0.0005, ymax) # force axis limits
        plt.grid()
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('$R_{rs}$ [sr$^{-1}$]')

        plt.subplot(2, 2, 2)
        plt.title(f" Windspeed > 18 m/s removed  (n = {int(np.sum(t_1))})")
        plt.plot(wl, rrs.T[:,t_1==1], linewidth=0.4, alpha=0.8)
        plt.xlim(350, 800)
        plt.ylim(-0.0005, ymax) # force axis limits
        plt.grid()
        plt.xlabel('Wavelength [nm]')


        plt.subplot(2, 2, 3)
        plt.title(f"Windspeed > 14 m/s removed  (n = {int(np.sum(t_2))})")
        plt.plot(wl, rrs.T[:,t_2==1], linewidth=0.4, alpha=0.8)
        plt.xlim(380,800 )
        plt.ylim(-0.0005, ymax) # force axis limits
        plt.grid()
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('$R_{rs}$ [sr$^{-1}$]')
        
        plt.subplot(2, 2, 4)
        plt.title(f"Windspeed > 10 m/s removed (n = {int(np.sum(t_3))})")
        plt.plot(wl, rrs.T[:,t_3==1], linewidth=0.4, alpha=0.8)
        plt.xlim(380,800)
        plt.ylim(-0.0005, ymax) # force axis limits
        plt.grid()
        plt.xlabel('Wavelength [nm]')


        plt.subplots_adjust(hspace=0.5)
        plt.savefig(os.path.join(target, file_id + '_QC.png'), format='png', dpi=200)
 
        return
    
    
    
def plot_rrs_qc_IDRthresholds(rrs, time, wl, q_3, t_1, t_2, t_3, file_id, target):
    """ rrs plot function showing sequential quality control filters"""

    if np.sum(q_3) > 0:
        plt.figure()
        plt.figure(figsize=(12, 8))
        plt.suptitle(str(file_id))

        ymax = np.ceil(np.nanmax(rrs.T[:, q_3==1]*1000))/1000;

        plt.subplot(2, 2, 1)
        plt.title(f"Default MONDA QC (n = {int(np.sum(q_3))})")
        plt.plot(wl, rrs.T[:,t_1==1], linewidth=0.4, alpha=0.6)
        plt.xlim(350, 800)
        plt.ylim(-0.0005, ymax) # force axis limits
        plt.grid()
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('$R_{rs}$ [sr$^{-1}$]')

        plt.subplot(2, 2, 2)
        plt.title(f" IDR > 0.3 removed (n = {int(np.sum(t_1))})")
        plt.plot(wl, rrs.T[:,t_1==1], linewidth=0.4, alpha=0.8)
        plt.xlim(350, 800)
        plt.ylim(-0.0005, ymax) # force axis limits
        plt.grid()
        plt.xlabel('Wavelength [nm]')


        plt.subplot(2, 2, 3)
        plt.title(f"IDR > 0.25 removed (n = {int(np.sum(t_2))})")
        plt.plot(wl, rrs.T[:,t_2==1], linewidth=0.4, alpha=0.8)
        plt.xlim(380,800 )
        plt.ylim(-0.0005, ymax) # force axis limits
        plt.grid()
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('$R_{rs}$ [sr$^{-1}$]')
        
        plt.subplot(2, 2, 4)
        plt.title(f"IDR > 0.2 removed (n = {int(np.sum(t_3))})")
        plt.plot(wl, rrs.T[:,t_3==1], linewidth=0.4, alpha=0.8)
        plt.xlim(380,800)
        plt.ylim(-0.0005, ymax) # force axis limits
        plt.grid()
        plt.xlabel('Wavelength [nm]')

        plt.subplots_adjust(hspace=0.5)
        plt.savefig(os.path.join(target, file_id + '_QC.png'), format='png', dpi=200)
 
        return

def plot_rrs_qc_combined_thresholds(rrs, time, wl, q_3, q_az, q_tiltavg, q_tiltstd, q_windbatos, q_windgir, q_IDR, q_4, q_5, file_id, target):
  
  
        plt.figure()
        plt.figure(figsize=(25, 10))
        plt.suptitle(str(file_id))
      
        if np.sum(q_3)>1:
            
            ymax = np.ceil(np.nanmax(rrs.T[:, q_3==1]*1000))/1000;
            if np.sum(~np.isnan(q_az))>1:
                plt.subplot(2, 4, 1)
                plt.title(f"Default azimuth QC (n = {int(np.sum(q_az))})")
                plt.plot(wl, rrs.T[:,q_az==1], linewidth=0.4, alpha=0.6)
                plt.xlim(350, 800)
                plt.ylim(-0.0005, ymax) # force axis limits
                plt.grid()
                plt.xlabel('Wavelength [nm]')
                plt.ylabel('$R_{rs}$ [sr$^{-1}$]')
    
            if np.sum(~np.isnan(q_tiltavg))>1:
                plt.subplot(2, 4, 2)
                plt.title(f"Tilt avg < 3 deg (n = {int(np.sum(q_tiltavg))})")
                plt.plot(wl, rrs.T[:,q_tiltavg==1], linewidth=0.4, alpha=0.8)
                plt.xlim(350, 800)
                plt.ylim(-0.0005, ymax) # force axis limits
                plt.grid()
                plt.xlabel('Wavelength [nm]')
            
            if np.sum(~np.isnan(q_tiltstd))>1:
                plt.subplot(2, 4, 3)
                plt.title(f"Tilt std < 2 deg (n = {int(np.sum(q_tiltavg))})")
                plt.plot(wl, rrs.T[:,q_tiltstd==1], linewidth=0.4, alpha=0.8)
                plt.xlim(350, 800)
                plt.ylim(-0.0005, ymax) # force axis limits
                plt.grid()
                plt.xlabel('Wavelength [nm]')
            
            q_wind = np.logical_and(q_windbatos,q_windgir)
            q_wind = np.ravel(q_wind)
            if np.sum(~np.isnan(q_windbatos))>1:
                plt.subplot(2, 4, 4)
                plt.title(f"Wind < 10 m/s (n = {int(np.sum(q_windbatos))})")
                plt.plot(wl, rrs.T[:,q_wind==1], linewidth=0.4, alpha=0.8)
                plt.xlim(380,800)
                plt.ylim(-0.0005, ymax) # force axis limits
                plt.grid()
                plt.xlabel('Wavelength [nm]')
            
            if np.sum(~np.isnan(q_IDR))>1:
                plt.subplot(2, 4, 5)
                plt.title(f"IDR < 0.25 (n = {int(np.sum(q_IDR))})")
                plt.plot(wl, rrs.T[:,q_IDR==1], linewidth=0.4, alpha=0.8)
                plt.xlim(380,800)
                plt.ylim(-0.0005, ymax) # force axis limits
                plt.grid()
                plt.xlabel('Wavelength [nm]')
    
    
            if np.sum(~np.isnan(q_4))>1:
                   plt.subplot(2, 4, 6)
                   plt.title(f"q_4 (all fields except IDR) (n = {int(np.sum(q_4))})")
                   plt.plot(wl, rrs.T[:,q_4==1], linewidth=0.4, alpha=0.8)
                   plt.xlim(380,800)
                   plt.ylim(-0.0005, ymax) # force axis limits
                   plt.grid()
                   plt.xlabel('Wavelength [nm]')
            
    
            if np.sum(~np.isnan(q_5))>1:
                 plt.subplot(2, 4, 7)
                 plt.title(f"q_5 (all fields ) (n = {int(np.sum(q_5))})")
                 plt.plot(wl, rrs.T[:,q_5==1], linewidth=0.4, alpha=0.8)
                 plt.xlim(380,800)
                 plt.ylim(-0.0005, ymax) # force axis limits
                 plt.grid()
                 plt.xlabel('Wavelength [nm]')

    
            plt.subplots_adjust(hspace=0.5)
            plt.savefig(os.path.join(target, file_id + '_QCall.png'), format='png', dpi=200)
  
        return  
  
def plot_rrs_qc_azimuth_thresholds(rrs, time, wl, q_3, t_3,file_id, target):
    """ rrs plot function showing sequential quality control filters"""

    if np.sum(q_3) > 0:
        plt.figure()
        plt.figure(figsize=(12, 8))
        plt.suptitle(str(file_id))

        ymax = np.ceil(np.nanmax(rrs.T[:, q_3==1]*1000))/1000;

        plt.subplot(2,1,1)
        plt.title(f"Default MONDA QC (n = {int(np.sum(q_3))})")
        plt.plot(wl, rrs.T[:,q_3==1], linewidth=0.4, alpha=0.6)
        plt.xlim(350, 800)
        plt.ylim(-0.0005, ymax) # force axis limits
        plt.grid()
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('$R_{rs}$ [sr$^{-1}$]')

        plt.subplot(2, 1, 2)
        plt.title(f" With azi. filt [reject > 170, < 100 degs] (n = {int(np.sum(t_3))})")
        plt.plot(wl, rrs.T[:,t_3==1], linewidth=0.4, alpha=0.8)
        plt.xlim(350, 800)
        plt.ylim(-0.0005, ymax) # force axis limits
        plt.grid()
        plt.xlabel('Wavelength [nm]')


        plt.subplots_adjust(hspace=0.5)
        plt.savefig(os.path.join(target, file_id + '_QC.png'), format='png', dpi=200)
 
        return


def tilt_avfilters():
    
    wl =  np.arange(340,900,1)
    target = os.path.join('.', 'So-Rad_test-output')
    for i in range(0,len(rrs_files)):
        print(i)
        print(str(meta_files[i]))
        print(str(rrs_files[i]))
        meta_i = meta_files[i]
        rrs_i = rrs_files[i]
        
        rrs_i = np.loadtxt(open(rrs_files[i],"rb"), delimiter=",") #
        meta_i = pd.read_csv(os.path.join(meta_files[i])) #
        
        time_i = [datetime.datetime.strptime(meta_i['timestamp'][j][0:19],'%Y-%m-%d %H:%M:%S') for j in range(len(meta_i))] # convert
        time_i = np.array(time_i)  # convert to np array for masking
        
        
        tilt_i = np.array(meta_i['tilt_avg'])
        q3_i = np.array(meta_i['q_3'])
        
        t1_i = ~np.isnan(np.where(tilt_i < 5, tilt_i, np.nan))
        t1_i = np.logical_and(t1_i, q3_i)
             
        t2_i = ~np.isnan(np.where(tilt_i < 3, tilt_i, np.nan))
        t2_i = np.logical_and(t2_i, q3_i)
            
        t3_i = ~np.isnan(np.where(tilt_i < 2, tilt_i, np.nan))
        t3_i = np.logical_and(t3_i, q3_i)
            
    
        date_i = str(time_i[0])[0:10]
        if np.std(tilt_i) > 0.05:
            plot_rrs_qc_tilthresholds(rrs_i, time_i, wl, q3_i, t1_i, t2_i, t3_i, date_i, target)

    return


def tilt_stdfilters():
    
    wl =  np.arange(340,900,1)
    target = os.path.join('.', 'So-Rad_test-output')
    for i in range(0,len(rrs_files)):
        print(i)
        print(str(meta_files[i]))
        print(str(rrs_files[i]))
        meta_i = meta_files[i]
        rrs_i = rrs_files[i]
        
        rrs_i = np.loadtxt(open(rrs_files[i],"rb"), delimiter=",") #
        meta_i = pd.read_csv(os.path.join(meta_files[i])) #
        
        time_i = [datetime.datetime.strptime(meta_i['timestamp'][j][0:19],'%Y-%m-%d %H:%M:%S') for j in range(len(meta_i))] # convert
        time_i = np.array(time_i)  # convert to np array for masking
        
        
        tilt_i = np.array(meta_i['tilt_std'])
        q3_i = np.array(meta_i['q_3'])
        
        t1_i = ~np.isnan(np.where(tilt_i < 3, tilt_i, np.nan))
        t1_i = np.logical_and(t1_i, q3_i)
             
        t2_i = ~np.isnan(np.where(tilt_i < 2, tilt_i, np.nan))
        t2_i = np.logical_and(t2_i, q3_i)
            
        t3_i = ~np.isnan(np.where(tilt_i < 1, tilt_i, np.nan))
        t3_i = np.logical_and(t3_i, q3_i)
            
    
        date_i = str(time_i[0])[0:10]
        if np.std(tilt_i) > 0.05:
            plot_rrs_qc_tilstdthresholds(rrs_i, time_i, wl, q3_i, t1_i, t2_i, t3_i, date_i, target)

    return


def wind_filters():
   
       wl =  np.arange(340,900,1)
       target = os.path.join('.', 'So-Rad_test-output')
       for i in range(0,len(rrs_files)):
           print(i)
           print(str(meta_files[i]))
           print(str(rrs_files[i]))
           meta_i = meta_files[i]
           rrs_i = rrs_files[i]
         
           rrs_i = np.loadtxt(open(rrs_files[i],"rb"), delimiter=",") #
           meta_i = pd.read_csv(os.path.join(meta_files[i])) #
                     
           if 'wind_speed_batos'  in  meta_i.keys(): 
               time_i = [datetime.datetime.strptime(meta_i['timestamp'][j][0:19],'%Y-%m-%d %H:%M:%S') for j in range(len(meta_i))] # convert
               time_i = np.array(time_i)  # convert to np array for masking
               
               
               wind_i = np.array(meta_i['wind_speed_batos'])
               q3_i = np.array(meta_i['q_3'])
               
               t1_i = ~np.isnan(np.where(wind_i < 18, wind_i, np.nan))
               t1_i = np.logical_and(t1_i, q3_i)
                    
               t2_i = ~np.isnan(np.where(wind_i < 14, wind_i, np.nan))
               t2_i = np.logical_and(t2_i, q3_i)
                   
               t3_i = ~np.isnan(np.where(wind_i < 10, wind_i, np.nan))
               t3_i = np.logical_and(t3_i, q3_i)
                   
           
               date_i = str(time_i[0])[0:10]
     
               plot_rrs_qc_windthresholds(rrs_i, time_i, wl, q3_i, t1_i, t2_i, t3_i, date_i, target)

       return
   
    
    
def IDR_filters():
   
       wl =  np.arange(340,900,1)
       target = os.path.join('.', 'So-Rad_test-output')
       

     
       for i in range(0,len(rrs_files)):
               print(i)
               print(str(meta_files[i]))
               print(str(rrs_files[i]))
               
               print(str(meta_files[i]))
               print(str(rrs_files[i]))
               meta_i = meta_files[i]
               rrs_i = rrs_files[i]
             
               rrs_i = np.loadtxt(open(rrs_files[i],"rb"), delimiter=",") #
               meta_i = pd.read_csv(os.path.join(meta_files[i])) #
                         
               
               if 'IDR_hsp'  in  meta_i.keys():    
                   meta_i = meta_files[i]
                   rrs_i = rrs_files[i]
                 
                   rrs_i = np.loadtxt(open(rrs_files[i],"rb"), delimiter=",") #
                   meta_i = pd.read_csv(os.path.join(meta_files[i])) #
                             
            
                   time_i = [datetime.datetime.strptime(meta_i['timestamp'][j][0:19],'%Y-%m-%d %H:%M:%S') for j in range(len(meta_i))] # convert
                   time_i = np.array(time_i)  # convert to np array for masking
                    
                    
                   IDR_i = np.array(meta_i['IDR_hsp'])
                   q3_i = np.array(meta_i['q_3'])
                    
                   t1_i = ~np.isnan(np.where(IDR_i < 0.3, IDR_i, np.nan))
                   t1_i = np.logical_and(t1_i, q3_i)
                         
                   t2_i = ~np.isnan(np.where(IDR_i < 0.25, IDR_i, np.nan))
                   t2_i = np.logical_and(t2_i, q3_i)
                        
                   t3_i = ~np.isnan(np.where(IDR_i < 0.2, IDR_i, np.nan))
                   t3_i = np.logical_and(t3_i, q3_i)
                        
            
                   date_i = str(time_i[0])[0:10]
                  
                   plot_rrs_qc_IDRthresholds(rrs_i, time_i, wl, q3_i, t1_i, t2_i, t3_i, date_i, target)

       return


    
    
def azimuth_filters():
   
       wl =  np.arange(340,900,1)
       target = os.path.join('.', 'So-Rad_test-output')
       

     
       for i in range(0,len(rrs_files)):
               print(i)
               print(str(meta_files[i]))
               print(str(rrs_files[i]))
               
               print(str(meta_files[i]))
               print(str(rrs_files[i]))
               meta_i = meta_files[i]
               rrs_i = rrs_files[i]
             
               rrs_i = np.loadtxt(open(rrs_files[i],"rb"), delimiter=",") #
               meta_i = pd.read_csv(os.path.join(meta_files[i])) #
                         
               
               if 'rel_view_az'  in  meta_i.keys():    
                   meta_i = meta_files[i]
                   rrs_i = rrs_files[i]
                 
                   rrs_i = np.loadtxt(open(rrs_files[i],"rb"), delimiter=",") #
                   meta_i = pd.read_csv(os.path.join(meta_files[i])) #
                             
            
                   time_i = [datetime.datetime.strptime(meta_i['timestamp'][j][0:19],'%Y-%m-%d %H:%M:%S') for j in range(len(meta_i))] # convert
                   time_i = np.array(time_i)  # convert to np array for masking
                    
                    
                   theta_i = np.array(meta_i['rel_view_az'])
                   q3_i = np.array(meta_i['q_3'])
                    
                   t1_i = ~np.isnan(np.where(np.abs(theta_i) < 170, theta_i, np.nan))
                   t2_i = ~np.isnan(np.where(np.abs(theta_i) > 100, theta_i, np.nan))
                  
                   t3_i = np.logical_and(t1_i, t2_i)  
                   
                   t3_i = np.logical_and(t3_i, q3_i)  
                  
                   date_i = str(time_i[0])[0:10]
                  
                   plot_rrs_qc_azimuth_thresholds(rrs_i, time_i, wl, q3_i, t3_i, date_i, target)

       return


def season_summaries(tilt_av, tilt_std, wind_speed_batos, rel_view_az, gps_speed, IDR, just_underway=False, just_stationary=False):
    
    plt.figure(figsize=(15,10))
    
    plt.suptitle('All data ' +  'N =' +  str(np.sum(~np.isnan(tilt_av))))
    plt.rcParams.update({'font.size': 16})
    
    if just_underway == True:
        
     
        tilt_av = tilt_av[gps_speed >3]
        tilt_std = tilt_std[gps_speed >3]
    
        wind_speed_batos = wind_speed_batos[gps_speed >3]
        rel_view_az = rel_view_az[gps_speed >3]

        IDR = IDR[gps_speed >3]
        gps_speed = gps_speed[gps_speed >3]
        
        plt.suptitle('GPS_speed > 3 m/s (approx. underway): ' +  'N =' +  str(np.sum(~np.isnan(tilt_av))))
        
    if just_stationary == True:
        
 
        tilt_av = tilt_av[gps_speed <1]
        tilt_std = tilt_std[gps_speed <1]
    
        wind_speed_batos = wind_speed_batos[gps_speed <1]
        rel_view_az = rel_view_az[gps_speed < 1]
        
        IDR = IDR[gps_speed < 1]

        gps_speed = gps_speed[gps_speed < 1]        
           
        plt.suptitle('GPS_speed < 3 m/s (approx. st breakpoint()ationary)' +  'N =' +  str(np.sum(~np.isnan(tilt_av))))
        
        
    plt.subplot(2,3,1)
    plt.title('Average tilt')
    plt.hist(tilt_av,bins=15,edgecolor='black')
    plt.xlim(0,10)
    plt.ylabel('Frequency')
    plt.xlabel('Deg.')
    
    plt.subplot(2,3,2)
    plt.title('Tilt STD')
    plt.hist(tilt_std,bins=15,edgecolor='black', color='red')
    plt.xlim(0,10)
    plt.ylabel('Frequency')
    plt.xlabel('Deg.')
    
    plt.subplot(2,3,3)
    plt.title('So-Rad GPS speed ')
    plt.hist(gps_speed,bins=15,edgecolor='black',color='orange')
    plt.xlim(0,16)
    plt.ylabel('Frequency')
    plt.xlabel('m/s')
        
    plt.subplot(2,3,4)
    plt.title('Relative_viewing azimuth')
    plt.hist(rel_view_az, bins=30,edgecolor='black',color='green')
    plt.xlim(0,180)
    plt.ylabel('Frequency')
    plt.xlabel('Deg.')
           
    plt.subplot(2,3,5)
    plt.title('Windspeed (BATOS)') # N =' +   str(np.sum(~np.isnan(wind_speed_batos))))
    plt.hist(wind_speed_batos, bins=20,edgecolor='black',color='yellow')
    plt.xlim(0,20)
    plt.ylabel('Frequency')
    plt.xlabel('Deg.')
    
    plt.subplot(2,3,6)
    plt.title('IDR (HSP)')
    plt.hist(IDR, bins=10,edgecolor='black',color='gray')
    plt.xlim(0,1)
    plt.ylabel('Frequency')
    plt.xlabel('[-]')
    
    plt.tight_layout()

    return

    

if __name__ == '__main__':
    
        ###########################################################################
    dir_main = '/users/rsg/tjor/TaraOcean/monda/MONDA/src/monda/tests/Tara_output_newcals/So-Rad_meta_with_TaraUnderway/Files/'
    #dir_main = '/users/rsg/tjor/TaraOcean/monda/MONDA/src/monda/tests/Tara_output_newcals/So-Rad_meta/'
   
    meta_files = sorted(glob.glob(os.path.join(dir_main, '*meta*'))) # 
        
    
    dir_rrs = '/users/rsg/tjor/TaraOcean/monda/MONDA/src/monda/tests/Tara_output_newcals/So-Rad_L2_MONDA/'
    rrs_files = sorted(glob.glob(os.path.join(dir_rrs, '*csv*'))) # 

    #
    meta_sr_list = []
    for i in range(len(meta_files)):  
            meta_sr = pd.read_csv(os.path.join(meta_files[i])) #
            meta_sr_list.append(meta_sr)
           
    # deployment wide-plots 
    meta_total =  meta_sr_list[0]
    for d in meta_sr_list[1:]:
        meta_total = pd.concat([meta_total,d])

        
    ###### global properties ####
    meta_total.keys()
    
    lat = np.array(meta_total['lat'])
    lon = np.array(meta_total['lon'])
    q = np.array(meta_total['q_3'])
    gps_speed = np.array(meta_total['gps_speed'])
    tilt_av = np.array(meta_total['tilt_avg'])
    tilt_std = np.array(meta_total['tilt_std'])
    rel_view_az = np.array(meta_total['rel_view_az'])
    wind_speed_batos = np.array(meta_total['wind_speed_batos'])
    IDR =  np.array(meta_total['IDR_hsp'])
    
    season_summaries(tilt_av, tilt_std, wind_speed_batos, rel_view_az, gps_speed, IDR, just_underway = False, just_stationary = False)
    season_summaries(tilt_av, tilt_std, wind_speed_batos, rel_view_az, gps_speed, IDR, just_underway = True, just_stationary = False)
    season_summaries(tilt_av, tilt_std, wind_speed_batos, rel_view_az, gps_speed, IDR, just_underway = False, just_stationary = True)


     
def all_filters_savemask():
   
       wl =  np.arange(340,900,1)
       target = os.path.join('.', 'So-Rad_test-output')
       for i in range(0,len(rrs_files)):
            print(i)
            print(str(meta_files[i]))
            print(str(rrs_files[i]))
    
            meta_i = meta_files[i]
            rrs_i = rrs_files[i]
              
            rrs_i = np.loadtxt(open(rrs_files[i],"rb"), delimiter=",") #
            meta_i = pd.read_csv(os.path.join(meta_files[i])) #
               
            time_i = [datetime.datetime.strptime(meta_i['timestamp'][j][0:19],'%Y-%m-%d %H:%M:%S') for j in range(len(meta_i))] # convert
            time_i = np.array(time_i)  # convert to np array for masking
                               
            q3_i = np.array(meta_i['q_3'])
                         
    
            # tilt_av 
            tiltavg_i = np.array(meta_i['tilt_avg'])
            q_tiltavg_i = ~np.isnan(np.where(tiltavg_i  < 3, tiltavg_i, np.nan))
            
            mode_tiltavg = mode(tiltavg_i)   # set constant values to nan 
            len_mode = len(tiltavg_i[mode_tiltavg==tiltavg_i])
            if len_mode > 5 :
                q_tiltavg_i[mode_tiltavg==tiltavg_i] = np.nan
                tiltavg_i[mode_tiltavg==tiltavg_i] = np.nan
                
            q_tiltavg_i = q_tiltavg_i.astype(int) 
            q_tiltavg_i = np.logical_and (q3_i,  q_tiltavg_i)
            q_tiltavg_i = np.ravel(q_tiltavg_i)
            
            
            # tilt_std
            tiltstd_i = np.array(meta_i['tilt_std'])
            q_tiltstd_i = ~np.isnan(np.where(tiltstd_i  < np.sqrt(3), tiltstd_i, np.nan))
            
            mode_tiltstd = mode(tiltstd_i)   # set constant values to nan 
            len_mode = len(tiltavg_i[mode_tiltstd==tiltstd_i])
            if len_mode > 5 :
                q_tiltstd_i[mode_tiltstd==tiltstd_i] = np.nan
                tiltstd_i[mode_tiltstd==tiltstd_i] = np.nan
                
            q_tiltstd_i = q_tiltstd_i.astype(int) 
            q_tiltstd_i = np.logical_and (q3_i,  q_tiltstd_i)
            q_tiltstd_i = np.ravel(q_tiltstd_i)
            
            # Wind batos
            q_windbatos_i  = np.nan*np.ones([1,len(q_tiltavg_i)])
            if 'wind_speed_batos'  in  meta_i.keys():    
                wind_batos_i = np.array(meta_i['wind_speed_batos'])
                q_windbatos_i = ~np.isnan(np.where(wind_batos_i < 10, wind_batos_i, np.nan))
                q_windbatos_i =  q_windbatos_i.astype(int) 
                q_windbatos_i = np.logical_and (q3_i, q_windbatos_i)
            q_windbatos_i = np.ravel(q_windbatos_i)
                
                
            # Wind gir
            q_windgir_i  = np.nan*np.ones([1,len(q_tiltavg_i)])
            if 'wind_speed_gir'  in  meta_i.keys():    
                wind_gir_i = np.array(meta_i['wind_speed_gir'])
                q_windgir_i = ~np.isnan(np.where(wind_gir_i < 10, wind_gir_i, np.nan))
                q_windgir_i =  q_windgir_i.astype(int) 
                q_windgir_i = np.logical_and (q3_i, q_windgir_i)
            q_windgir_i = np.ravel(q_windgir_i)
                
            # IDR_HSP
            q_IDR_i  = np.nan*np.ones([1,len(q_tiltavg_i)])
            if 'IDR_hsp'  in  meta_i.keys():    
                IDR_i = np.array(meta_i['IDR_hsp'])
                q_IDR_i = ~np.isnan(np.where(IDR_i < 0.25, IDR_i, np.nan))
                q_IDR_i =  q_IDR_i.astype(int) 
                q_IDR_i = np.logical_and (q3_i, q_IDR_i)
            q_IDR_i = np.ravel(q_IDR_i)
                    
           # azimuth        
            q_az_i  = np.nan*np.ones([1,len(q_tiltavg_i)])
            if  i > 21: # accounts for days with no reliable azimuth data   
               theta_i = np.array(meta_i['rel_view_az'])
               t1_i = ~np.isnan(np.where(np.abs(theta_i) < 170, theta_i, np.nan))
               t2_i = ~np.isnan(np.where(np.abs(theta_i) > 100, theta_i, np.nan))
               q_az_i = np.logical_and(t1_i, t2_i)
               q_az_i =  q_az_i.astype(int) 
               q_az_i = np.logical_and (q3_i, q_az_i)
            q_az_i = np.ravel(q_az_i)
                
            q4_i = np.logical_and(q_az_i, np.logical_and(q_tiltavg_i,  np.logical_and(q_tiltstd_i, np.logical_and(q_windbatos_i, q_windgir_i))))
            q5_i = np.logical_and(q_az_i, np.logical_and(q_tiltavg_i,  np.logical_and(q_tiltstd_i, np.logical_and(q_windbatos_i, np.logical_and(q_windgir_i, q_IDR_i)))))
                  
            q4_i =  q4_i.astype(int) 
            q5_i =  q5_i.astype(int) 
            
            q4_i = np.ravel(q4_i)
            q5_i = np.ravel(q5_i)
        
        
            # Add extra columnns to meta_i
            meta_i['q_azimuth'] = q_az_i
            meta_i['q_tiltavg'] = q_tiltavg_i
            meta_i['q_tiltstd'] = q_tiltstd_i
            meta_i['q_windspd_batos'] = q_windbatos_i
            meta_i['q_windspd_gir'] = q_windgir_i
            meta_i['q_IDR'] = q_IDR_i
            meta_i['q_4'] = q4_i
            meta_i['q_5'] = q5_i
            
      
            save_dir = '/users/rsg/tjor/TaraOcean/monda/MONDA/src/monda/tests/Tara_output_newcals/So-Rad_meta_with_TaraUnderway/'
            date_i= str(time_i[0])[0:10]  
            meta_i.to_csv(save_dir + date_i + '_3C_metadata_V2.csv')
            plot_rrs_qc_combined_thresholds(rrs_i, time_i, wl, q3_i, q_az_i, q_tiltavg_i, q_tiltstd_i, q_windbatos_i, q_windgir_i, q_IDR_i, q4_i, q5_i, date_i, target)

       return

    


# plt.scatter(wind_speed_batos,tilt_std)

# tilt_avg =tilt_avg[np.isnan(wind_speed_batos)==0]
# tilt_std = tilt_std[np.isnan(wind_speed_batos)==0]
# wind_speed_batos = wind_speed_batos[np.isnan(wind_speed_batos)==0]

# scipy.stats.pearsonr(wind_speed_batos,  tilt_avg)
# scipy.stats.pearsonr(np.log10(wind_speed_batos),  np.log10(tilt_avg))

# plt.scatter(np.log10(wind_speed_batos),  np.log10(tilt_avg))