#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 11:35:43 2023

@author: tjor
Tara Ocean - meta data

"""

import os
import sys
import csv
import netCDF4 as nc
import pandas as pd

import numpy as np
import glob   
import pandas as pd
import ephem
import xarray as xr

import matplotlib 
import matplotlib.pyplot as plt
import matplotlib.cm as cm


import math
import matplotlib.dates as mdates
import datetime as dt
from scipy.interpolate import interp1d
import scipy.io as io
import datetime


import pandas as pd
import glob

import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from xml.dom.minidom import *

import scipy
from scipy import stats
from scipy import signal as sg

import os

import warnings

import pyproj

from True_wind import truew


def plot_meta(df_sr):
    
     plt.figure(figsize=(15,12),dpi=300)
     plt.suptitle(str(df_sr.index[0])[0:10])
     
     # plt.rc('font', size=14)
     # plt.suptitle(str(df_sr.index.date[0]))
     # plt.subplot(3,3,1)  
     # plt.title('QC mask')
     # plt.plot_date(df_sr.index, df_sr['q_3'], label= 'MONDA default', color='blue') 
     # plt.xlabel('Time')
     # plt.xticks(rotation=45)
     
     plt.subplot(2,3,1)  
     plt.title('Tilt')
     plt.plot_date(df_sr.index, df_sr['tilt_avg'], color='blue') 
     plt.xlabel('Time')
     plt.ylabel('Deg')
     plt.xticks(rotation=45)
     ax = plt.gca()
     ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%HH') )
     
     plt.subplot(2,3,2)  
     plt.title('Relative viewing aziumuth')
     plt.plot_date(df_sr.index, df_sr['rel_view_az'], color='blue') 
     plt.xlabel('Time')
     plt.ylabel('Deg')
     plt.xticks(rotation=45)
     plt.xticks(rotation=45)
     ax = plt.gca()
     ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%HH'))
     
     
     plt.subplot(2,3,3)  
     plt.title('Wind speed')
     if 'wind_speed_batos'  in df_sr.keys():
         plt.plot_date(df_sr.index, df_sr['wind_speed_batos'],label='BATOS', color='red') 
     if 'wind_speed_gir'  in df_sr.keys():
        plt.plot_date(df_sr.index, df_sr['wind_speed_gir'],label='gir',color='green') 
     plt.xlabel('Time')
     plt.ylabel('m/s')
     plt.xticks(rotation=45)
     plt.legend()
     plt.xticks(rotation=45)
     ax = plt.gca()
     ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%HH') )
     
     
     plt.subplot(2,3,4)  
     plt.title('Heading & wind angle')
     if 'heading_gps32'  in df_sr.keys():
         plt.plot_date(df_sr.index, df_sr['heading_gps32'], label='heading: GPS32', color='black') 
     if 'heading_sr'  in df_sr.keys():
             plt.plot_date(df_sr.index, df_sr['heading_sr'], label='heading: So_rad', color='blue') 
     if 'wind_angle_batos'  in df_sr.keys():
         plt.plot_date(df_sr.index, df_sr['wind_angle_batos'], label='wind: BATOS', color='red') 
     if 'wind_angle_gir'  in df_sr.keys():
         plt.plot_date(df_sr.index, df_sr['wind_angle_gir'], label='wind: gir',color='green') 
     plt.xlabel('Time')
     plt.ylabel('Deg')
     plt.xticks(rotation=45)
     plt.legend()
     plt.xticks(rotation=45)
     ax = plt.gca()
     ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%HH') )
     

     plt.subplot(2,3,5)  
     plt.title('Platform speed')
     plt.plot_date(df_sr.index, df_sr['gps_speed'], label='SoRad', color='blue') 
     if 'spd_over_grnd_batos'  in df_sr.keys():
         plt.plot_date(df_sr.index, df_sr['spd_over_grnd_batos'], label='BATOS', color='red') 
     if 'spd_over_grnd_gps32'  in df_sr.keys():
         plt.plot_date(df_sr.index, df_sr['spd_over_grnd_gps32'], label='GPS32', color='black')
     plt.xlabel('Time')
     plt.ylabel('m/s')
     plt.xticks(rotation=45)     
     plt.legend()
     plt.tight_layout()
     plt.xticks(rotation=45)
     ax = plt.gca()
     ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%HH') )
     
     plt.subplot(2,3,6)  
     plt.title('Diffuse light fraction')
     if 'IDR_hsp' in df_sr.keys():
             plt.plot_date(df_sr.index, df_sr['IDR_hsp'], color='black') 
             plt.xlabel('Time')
             plt.ylabel('IDR')
             plt.xticks(rotation=45)     
             plt.legend()
             plt.xticks(rotation=45)
             plt.ylim(0,1)
     
     ax = plt.gca()
     ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%HH') )
             
     plt.tight_layout()
     
     filename =     dir_combined_meta + files_sr[i].split("/")[-1][:-4] + '_'  'metaplot.png'
     plt.savefig(filename, dpi=200)
     plt.close() 
     
     return
 
    
def append_batos(files_batos, df_sr, date_str):
    
    # loads batos data stream, and re-samples to so-rad timestamps    
    
    df_batos = pd.read_csv(files_batos[0], skiprows=[1]) # open zeroth file
    df_batos['time_batos'] = [datetime.datetime.strptime(str(df_batos['time'][i][0:19]),'%Y/%m/%d %H:%M:%S') for i in range(len(df_batos))]
    df_batos = df_batos.set_index('time_batos')
    df_batos = df_batos.groupby(df_batos.index, group_keys=True, dropna=False).agg(np.nanmedian) # group repeatad time in
    
    for j in range(1,len(files_batos)): # loop over and concatanate remaining batos files in day       
        df_batos_j = pd.read_csv(files_batos[j], skiprows=[1])
        df_batos_j['time_batos'] = [datetime.datetime.strptime(str(df_batos_j['time'][i][0:19]),'%Y/%m/%d %H:%M:%S') for i in range(len(df_batos_j))]
        df_batos_j = df_batos_j.set_index('time_batos')
        df_batos_j = df_batos_j.groupby(df_batos_j.index, group_keys=True, dropna=False).agg(np.nanmedian) # group repeatad time indices    
        df_batos = pd.concat([df_batos, df_batos_j]) 
    
    # re-sample batos data to so-rad timestamps
    df_batos = df_batos[~df_batos.index.duplicated(keep='first')]
    df_batos = df_batos.reindex(df_sr.index,method='nearest',tolerance='15s') 
    
    # append relevant batos fields to so-rad dataframe
    if ' truecourse'  in df_batos.keys():
        df_sr['true_course_batos'] = df_batos[' true_course'] 
    else:
        df_sr['true_course_batos'] = np.nan*np.ones(len(df_sr))
    
    if ' spd_over_grnd'  in df_batos.keys():
        df_sr['spd_over_grnd_batos'] = df_batos[' spd_over_grnd'] 
    else:
        df_sr['spd_over_grnd_batos'] = np.nan*np.ones(len(df_sr))
    
    if ' wind_angle'  in df_batos.keys():
        df_sr['wind_angle_batos'] =   df_batos[' wind_angle'] 
    else:
        df_sr['wind_angle_batos'] = np.nan*np.ones(len(df_sr))

    if ' wind_speed'  in df_batos.keys():
        df_sr['wind_speed_batos'] =  df_batos[' wind_speed']**0.5144444444 
    else:
        df_sr['wind_speed_batos'] = np.nan*np.ones(len(df_sr))

    if ' temperature'  in df_batos.keys():
        df_sr['temperature_batos'] =  df_batos[' temperature'] 
    else:
        df_sr['temperature_batos'] = np.nan*np.ones(len(df_sr))

    return df_sr


def append_gps32(files_gps, df_sr, date_str):
    
     # loads gps 32 datastream (file format 1), and re-samples to so-rad timestamps    
     df_gps32 = pd.read_csv(files_gps[0], skiprows=[1]) # open zeroth file   
     df_gps32['time_gps32'] = [datetime.datetime.strptime(str(df_gps32['time'][i][0:19]),'%Y/%m/%d %H:%M:%S') for i in range(len(df_gps32))]
     df_gps32 = df_gps32.set_index('time_gps32')
     df_gps32 = df_gps32.groupby(df_gps32.index, group_keys=True, dropna=False).agg(np.nanmedian)
     for j in range(1,len(files_gps)): # loop over and concatanate remaining batos files in day       
         df_gps32_j = pd.read_csv(files_gps[j], skiprows=[1]) # open zeroth file   
         df_gps32_j['time_gps32'] = [datetime.datetime.strptime(str(df_gps32_j['time'][i][0:19]),'%Y/%m/%d %H:%M:%S') for i in range(len(df_gps32_j))]
         df_gps32_j = df_gps32_j.set_index('time_gps32')
         df_gps32_j = df_gps32_j.groupby(df_gps32_j.index,group_keys=True, dropna=False).agg(np.nanmedian)  
         df_gps32 = pd.concat([df_gps32, df_gps32_j]) 
     
     df_gps32 = df_gps32[~df_gps32.index.duplicated(keep='first')]
     df_gps32 = df_gps32.reindex(df_sr.index,method='nearest', tolerance='15s') # resample so-rad timestam      
    
     if ' latitude'  in df_gps32.keys():
        df_sr['lat_gps32'] = df_gps32[' latitude'] 
     else:
        df_sr['lat_gps32'] = np.nan*np.ones(len(df_sr))
    
     if ' longitude'  in df_gps32.keys():
        df_sr['lon_gps32'] = df_gps32[' longitude'] 
     else:
        df_sr['lon_gps32'] = np.nan*np.ones(len(df_sr))
     
     if ' gps_qual'  in df_gps32.keys():
        df_sr['gps_qual_gps32'] = df_gps32[' gps_qual'] 
     else:
        df_sr['lon_gps32'] = np.nan*np.ones(len(df_sr))
     
     if ' heading'  in df_gps32.keys():
        df_sr['heading_gps32'] = df_gps32[' heading'] 
     else:
        df_sr['heading_gps32'] = np.nan*np.ones(len(df_sr))
     
     if ' true_course'  in df_gps32.keys():
        df_sr['true_course_gps32'] = df_gps32[' true_course'] 
     else:
        df_sr['true_course_gps32'] = np.nan*np.ones(len(df_sr))      
     
     if ' spd_over_grnd'  in  df_gps32.keys():
        df_sr['spd_over_grnd_gps32'] = df_gps32[' spd_over_grnd'] 
     else:
        df_sr['spd_over_grnd_gps32'] = np.nan*np.ones(len(df_sr))

     return df_sr


def append_gps32_viaship(files_gps_viaship, df_sr, date_str):
    
    # loads gps 32 datastream (file format 1), and re-samples to so-rad timestamps
    # for now, just considers heading (other datafields are amiguous)    

    df_gps32_index = pd.read_csv(files_gps_viaship[0], sep =',',skiprows=2, header=None)
    df_gps32 = pd.DataFrame()          
    
    
    #df_gps32_index.iloc[:,0]  - Lat: neglect for now
    #df_gps32_index.iloc[:,1]  - Lon: neglect for now
    time = np.nan*np.ones(len(df_gps32_index), dtype=object)
    for i in range(len(df_gps32_index)):
        year_i =df_gps32_index.iloc[:,2][i]
        month_i = df_gps32_index.iloc[:,3][i]
        day_i = df_gps32_index.iloc[:,4][i]
        hour_i = df_gps32_index.iloc[:,5][i]
        min_i = df_gps32_index.iloc[:,6][i]
        sec_i = df_gps32_index.iloc[:,7][i]
        time[i] = dt.datetime(year_i, month_i, day_i, hour_i, min_i, sec_i)
    
    df_gps32['time_gps32'] = time
    df_gps32['heading_gps32'] =  df_gps32_index.iloc[:,10]
    
    df_gps32 = df_gps32.set_index('time_gps32')
    df_gps32 = df_gps32.groupby(df_gps32.index, group_keys=True, dropna=False).agg(np.nanmedian)
        
    df_gps32 = df_gps32[~df_gps32.index.duplicated(keep='first')]
    df_gps32 = df_gps32.reindex(df_sr.index,method='nearest', tolerance='15s') # resample so-rad timestam      
          
    if 'heading_gps32'  in df_gps32.keys():
       df_sr['heading_gps32'] = df_gps32['heading_gps32'] 
    else:
       df_sr['heading_gps32'] = np.nan*np.ones(len(df_sr))
       
    return df_sr


def append_gir(files_gir, df_sr, date_str):
    
    df_gir = pd.read_csv(files_gir[0], skiprows=[1]) # open zeroth file   
    df_gir['time_gir'] = [datetime.datetime.strptime(str(df_gir['time'][i][0:19]),'%Y/%m/%d %H:%M:%S') for i in range(len(df_gir))]
    df_gir = df_gir.set_index('time_gir')
    df_gir = df_gir.groupby(df_gir.index, group_keys=True, dropna=False).agg(np.nanmedian)
                
    for j in range(1,len(files_gir)): # loop over and concatanate remaining batos files in day       
        df_gir_j = pd.read_csv(files_gir[j], skiprows=[1]) # open zeroth file   
        df_gir_j['time_gir'] = [datetime.datetime.strptime(str(df_gir_j['time'][i][0:19]),'%Y/%m/%d %H:%M:%S') for i in range(len(df_gir_j))]
        df_gir_j = df_gir_j.set_index('time_gir')
        df_gir_j = df_gir_j.groupby(df_gir_j.index,group_keys=True, dropna=False).agg(np.nanmedian)  
        df_gir = pd.concat([df_gir, df_gir_j]) 
  
    df_gir = df_gir[~df_gir.index.duplicated(keep='first')]
    df_gir = df_gir.reindex(df_sr.index,method='nearest', tolerance='15s') 
    # breakpoint()  
      
    df_sr['wind_speed_gir'] =   df_gir[' wind_speed']*0.5144444444  # converts to m/s
    df_sr['wind_angle_gir'] =   df_gir[' wind_angle']

    return df_sr


def append_heading_fromlatlon(df_sr):
    
    if(len(df_sr))>1:
        heading = np.nan*np.ones(len(df_sr))
        geodesic = pyproj.Geod(ellps='WGS84')
        if len(df_sr) > 3:
            for i in range(1,len(df_sr) -1):
                heading[i], back_azimuth, distance = geodesic.inv(df_sr['lon'].iloc[i-1], df_sr['lat'].iloc[i-1], df_sr['lon'].iloc[i+1], df_sr['lat'].iloc[i+1])
                if heading [i] < 0:
                    heading[i] = heading[i] +360
        heading[0] = heading[1]
        heading[-1] = heading[-2]
        
        df_sr['heading_sr'] = heading          
   
    return df_sr




def angular_diff(df_sr, field1, field2):
    
   # calculates offsets between headings and wind directions 

   #  field1 = 'wind_angle_batos'
   #  field2 = 'heading_sr'    
   
    newfield = field1 + '_minus_' + field2
    
    theta_r = np.array(df_sr[field1].values -df_sr[field2].values)
    for i in range(len(theta_r)):
        if theta_r[i] < 0:
            theta_r[i] = theta_r[i] + 360
    
    df_sr[newfield] = theta_r

   # plt.figure()
   # plt.plot(np.array(df_sr[field1].values))
   # plt.plot(np.array(df_sr[field2].values))

    return df_sr


def calc_true_wind_sr_batos(df_sr):
     
    # calculates true_wind between headings and wind directions 
    # approximate course of ship by heading 
    
    true_wind_angle = np.nan*np.ones(len(df_sr))
    true_wind_speed = np.nan*np.ones(len(df_sr))
    
    for i in range(len(df_sr)):
            crse_i = df_sr['heading_sr'].iloc[i] # assume heading == course
            cspd_i = df_sr['gps_speed'].iloc[i]
            wdir_i = df_sr['wind_angle_batos'].iloc[i]
            zlr_i = float(78)
            hd_i = df_sr['heading_sr'].iloc[i]
            wspd_i = df_sr['wind_speed_batos'].iloc[i]
            
            true_wind_vec = truew( crse_i,   cspd_i ,  wdir_i, zlr_i, hd_i , wspd_i)
            true_wind_angle[i] = true_wind_vec[0]
            true_wind_speed[i] = true_wind_vec[1]
    
    df_sr['true_wind_angle'] = true_wind_angle
    df_sr['true_wind_speed'] = true_wind_speed

    return df_sr


def import_HSP_radiance(radiancefile,year):  # function to import radiance files (diffuse or total)
    date_time_list=[]        # list which stores date time objects
    spectra_list=[] # list which stores radiance spectra
 
    filename=open(radiancefile,"r",encoding='utf-8-sig')
    print(filename)
    for i, row in enumerate(filename):
            line = filename.readline()
            if line[:4] == year: # tests for line which contains radiance entries (first 4 elements are year)
                date_time_string=line[0:19]
                if not date_time_string:  # break conditon for empty string
                    break
                spectra=line[20:-1]
                try:
                     all(isinstance(x, float) for x in spectra) # checks for encoding errors in data (issues with 2019)
                     spectra_list.append(np.array(spectra.split("\t"),dtype=float))  # radiance
                     date_time_list.append(datetime.datetime.strptime(date_time_string, '%Y-%m-%d %H:%M:%S'))
                except ValueError:
                    break  
                
    return  spectra_list, date_time_list 

       
def append_IDR_hsp(folder_hsp ,df_sr):
           
        file_hsp_eds =  folder_hsp + '/Diffuse.txt' # total set of so-rad files
        file_hsp_ed = folder_hsp +'/Total.txt' # total set of so-rad files
        
        eds_hsp, time_hsp = import_HSP_radiance(file_hsp_eds , '2023')           
        ed_hsp, time_hsp = import_HSP_radiance(file_hsp_ed , '2023')
        
        eds_hsp = np.array(eds_hsp)
        ed_hsp = np.array(ed_hsp)

        IDR = np.sum(eds_hsp[:,86:601],axis=1)/np.sum(ed_hsp[:,86:601],axis=1)
           
        df_hsp = pd.DataFrame()
        df_hsp['time_hsp'] = time_hsp
        df_hsp['IDR'] = IDR
        df_hsp = df_hsp.set_index('time_hsp')
        df_hsp = df_hsp.reindex(df_sr.index,method='nearest',tolerance='120s')
        
        df_sr['IDR_hsp'] = df_hsp['IDR'] 
        
        return df_sr


if __name__ == '__main__':    
    
    warnings.filterwarnings("ignore") # turns off empty slice error
    
    dir_Tara_gps32 = '/users/rsg/tjor/TaraOcean/TO_Uway_metadata/2023/GPS32Tara/'
    dir_Tara_gps32_viaship = '/users/rsg/tjor/TaraOcean/TO_Uway_metadata/2023/GPS32Tara/GPS32 DATA_via_ship_software/'
    dir_Tara_batos = '/users/rsg/tjor/TaraOcean/TO_Uway_metadata/2023/BATOSTara/'
    dir_Tara_gir = '/users/rsg/tjor/TaraOcean/TO_Uway_metadata/2023/GirouetteTara/'
    dir_hsp ='/users/rsg/tjor/TaraOcean/monda/MONDA/src/monda/tests/HSP-1'
        
    dir_SoRad_meta = '/users/rsg/tjor/TaraOcean/monda/MONDA/src/monda/tests/Tara_output_newcals/So-Rad_meta/'
    dir_combined_meta = '/users/rsg/tjor/TaraOcean/monda/MONDA/src/monda/tests/Tara_output_newcals/So-Rad_meta_with_TaraUnderway/'
    
    files_sr = sorted(glob.glob(dir_SoRad_meta +'/*csv')) # total set of so-rad files for 2023
    
    # OFFSET_gir = [] #
    # OFFSET_batos = [] #
    for i in range(1,len(files_sr)):
        
        #### So-RAD ####
        # load sr data on ith day
        df_sr = pd.read_csv(files_sr[i])
        df_sr['time_sr'] = [datetime.datetime.strptime(str(df_sr['timestamp'][i][0:19]),'%Y-%m-%d %H:%M:%S') for i in range(len(df_sr))]
        df_sr = df_sr.set_index('time_sr')
        df_sr = df_sr.drop('Unnamed: 0', axis=1)  
        df_sr = append_heading_fromlatlon(df_sr)
        
        # caulculate bearing        
        date_str = str(df_sr.index[0])[0:4] + str(df_sr.index[0])[5:7] + str(df_sr.index[0])[8:10]  
        print(date_str)
        
        #### BATOS #######
        # load batos data on ith day - 24 files per day 
        files_batos = sorted(glob.glob(dir_Tara_batos +  '/*' +  date_str +'*csv')) 
        if len(files_batos) > 1:
            print('Appending BATOS')
            df_sr = append_batos(files_batos, df_sr, date_str)
            
        #### GPS #######
        files_gps = sorted(glob.glob(dir_Tara_gps32 +  '/*' +  date_str +'*csv')) # total set of so-rad files
        files_gps_viaship = sorted(glob.glob(dir_Tara_gps32_viaship +  '/*' +  date_str +'*csv')) 
        if len(files_gps) > 1: 
            print('Appending GPS')
            # print(len(files_gps))
            df_sr = append_gps32(files_gps, df_sr, date_str)
        elif len(files_gps_viaship) == 1: # use ship GPS as back-up #
            print('Appending GPS via ship')
            # print(files_gps_viaship) #
            df_sr = append_gps32_viaship(files_gps_viaship, df_sr, date_str)
            
        #### GIR ####### 
        date_str_gir = str(df_sr.index[0])[0:4]  + str(df_sr.index[0])[5:7] + str(df_sr.index[0])[8:10]  
        files_gir = sorted(glob.glob(dir_Tara_gir +  '/*' +  date_str +'*csv')) # 
        if len(files_gir) > 1:
           print('Appending Girouette')
           df_sr = append_gir(files_gir, df_sr, date_str)
           
        ### HSP1 ####
        date_str_hsp = date_str[0:4] + '-' + date_str[4:6] + '-' + date_str[6:8]
        folder_hsp =  dir_hsp +  '/' +  date_str_hsp 
        if os.path.isdir(folder_hsp) ==True:
            print('Appending HSP')
            append_IDR_hsp(folder_hsp, df_sr)
        

        # output
        plot_meta(df_sr)
        
        new_file_name = files_sr[i].split('/')[12][:-4] + '_withTO.csv'
        df_sr.to_csv('/users/rsg/tjor/TaraOcean/monda/MONDA/src/monda/tests/Tara_output_newcals/So-Rad_meta_with_TaraUnderway/' + new_file_name)
        
 
        ######## True Wind #######   
        #  df_sr = calc_true_wind_sr_batos(df_sr)
        # plot_meta(df_sr)
        # if 'wind_angle_batos'  in df_sr.keys() and 'heading_gps32' in df_sr.keys():
           # df_sr = angular_diff(df_sr, 'wind_angle_batos','heading_gps32')
           # offset_batos = np.nanmedian(np.array(df_sr['wind_angle_batos_minus_heading_gps32'].values)[np.array(df_sr['gps_speed'])<1])
           # print('wind_angle_batos_minus_heading_gps32 median  = ' + str(offset_batos))
           # OFFSET_batos.append(offset_batos)
            # if 'wind_angle_gir'  in df_sr.keys() and 'heading_sr' in df_sr.keys():
            #   df_sr = angular_diff(df_sr, 'wind_angle_gir','heading_sr')
            #  C = np.nanmedian(np.array(df_sr['wind_angle_gir_minus_heading_sr'].values)[np.array(df_sr['gps_speed'])<1])
            # print('wind_angle_gir_minus_heading_sr  median  = ' + str(C))
             
        # if 'wind_angle_gir'  in df_sr.keys() and 'heading_gps32' in df_sr.keys():
         #   df_sr = angular_diff(df_sr, 'wind_angle_gir','heading_gps32')
          #  offset_gir= np.nanmedian(np.array(df_sr['wind_angle_gir_minus_heading_gps32'].values)[np.array(df_sr['gps_speed'])<1])
           # print('wind_angle_gir_minus_heading_gps32  median  = ' + str(offset_gir)) # use 161?
            #OFFSET_gir.append(offset_gir)
        

