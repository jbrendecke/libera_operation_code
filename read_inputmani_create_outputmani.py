#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 29 13:45:35 2025

@author: jbrendecke
"""
import numpy as np
import xarray as xr
import datetime
from libera_utils import Manifest, smart_open, DataProductConfig
import glob
import os


#def read_manifest(input_manifest_file):
    #//////////////////////////////////////
    # Reads in json input manifest file
    # and returns xarray dataset sets for
    # SSF properties, ERA5 levels, and ERA5 surface
    #/////////////////////////////////////
    
    
# # # # # Read Input Manifest Files # # # # #
input_manifest_path=glob.glob("/home/jbrendecke/L2/data/*.json")[0]

input_manifest = Manifest.from_file(input_manifest_path)

#put filename and datasets in manifest files into dictionary
all_data={}
for i, file_info in enumerate(input_manifest.files):
    try: 
        with smart_open(file_info.filename) as file_handle:
            dataset = xr.open_dataset(file_handle)
            dataset.load()
            all_data[file_info.filename] = dataset
    except NameError:
        print("Error")

#save individual dataset by product 
for mfile, mds in all_data.items():
    if 'CERES_SSF' in mfile:
        ds_ssf = mds
    if 'ERA5_level' in mfile:
        ds_era5lev = mds
    if 'ERA5_surface' in mfile:
        ds_era5sfc = mds
   
#Make sure all files are there
try:
    ds_ssf
except NameError:
    print('CERES SSF file not found')

try:
    ds_era5lev
except NameError:
    print('ERA5 levels file not found')
    
try:
    ds_era5sfc
except NameError:
    print('ERA5 surface file not found')

#return ds_ssf, ds_era5lev, ds_era5sfc


#%%
#make/do the science output variables as an example
time_stamp = np.full((100,), np.datetime64('2020-07-11T10:00:00'))
for i in range(99):
    time_stamp[i+1] = time_stamp[i] + np.timedelta64(1,'h') 
start_time = time_stamp[0].astype(datetime.datetime)
end_time = time_stamp[-1].astype(datetime.datetime)
latitude = np.linspace(-90,90,100)
longitude = np.linspace(0,360,100)
SWD_TOA_VIS = np.full(100, 816)
SWD_TOA_NIR = np.full(100, 545)

#%%
# # # # # Write the NetCDF file(s) # # # # #
#path to .yml file
yml_filepath = '/home/jbrendecke/L2/L2_calculations.yml'
# use libera-utils to config object using .yml file
product_config = DataProductConfig.from_data_config_file(yml_filepath)

#add science variables to config object, "name variable" must be same as in .yml
product_config.add_data_to_variable("time_stamp", time_stamp)
product_config.add_data_to_variable("latitude", latitude)
product_config.add_data_to_variable("longitude", longitude)
product_config.add_data_to_variable("SWD_TOA_VIS", SWD_TOA_VIS)
product_config.add_data_to_variable("SWD_TOA_NIR", SWD_TOA_NIR)

#in terminal set environment vairale: PROCESSING_PATH:
# # export PROCESSING_PATH="/home/jbrendecke/rad_global/ex_libera_op/ex_output_netcdf/"
#
output_folder = os.getenv("PROCESSING_PATH")
if not output_folder:
    raise ValueError("PROCESSING_PATH environment variable is not set")
    
output_path = product_config.write( 
                    folder_location=output_folder,
                    start_time = start_time,
                    end_time = end_time,
                    allow_incomplete=True) 

#%%
# # # # # Write the Manifest output file(s) # # # # #
# 1. Create a new, empty OUTPUT manifest linked to our INPUT manifest
output_manifest = Manifest.output_manifest_from_input_manifest(input_manifest)

# 2. Add the new files you created to it (can be local or S3 paths)
# (Here, we'll use the filename we generated earlier)
new_data_file_path = glob.glob("/home/jbrendecke/L2/data/*.nc")[0]
output_manifest.add_files(new_data_file_path)

# 3. Write the final manifest to a specified location
output_destination = os.getenv("PROCESSING_PATH")
output_manifest.write(output_destination)



#////////////////////////////////////
# convert input dictionary to  netcdf file and output manifest file
# uses .yml file will list of expected outputs definitions/units for netcdf file
# dictionary variable names must mathc .yml variable names
# 
#
# Inputs:
#       dictionary of CCCma outputs for entire day
#       includes: SWD/SWU TOA, SWD/SWU SFC, Direct/Diffuse All-sky & Clear Sky
#
# Output:
#   xarray dataset for masked ERA5 surface dataset
#       indcludes: pressure, height(km), temp, rh, specific humidity, ozone
#       dims: time, lev, lat, lon
#/////////////////////////////////////



