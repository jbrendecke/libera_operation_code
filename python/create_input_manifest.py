#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 29 13:01:39 2025
Example of creating Manifest File for files needed to run operational code.
@author: jbrendecke
"""
import glob
from libera_utils import Manifest, ManifestType

#make intial manifest file
input_manifest = Manifest(manifest_type= ManifestType.INPUT)

# load Netcdf files
yymmdd='20180322'
filessf = glob.glob('/home/jbrendecke/L2/input_data/CERES/*'+yymmdd+'*.nc')[0]
file_lev = glob.glob('/home/jbrendecke/L2/input_data/ERA5/levels/*'+yymmdd+'*.nc')[0]
file_sfc = glob.glob('/home/jbrendecke/L2/input_data/ERA5/surface/*'+yymmdd+'*.nc')[0]

#filessf = glob.glob('/home/CCCma/L2/input_data/CERES/*'+yymmdd+'*.nc')[0]
#file_lev = glob.glob('/home/CCCma/L2/input_data/ERA5/levels/*'+yymmdd+'*.nc')[0]
#file_sfc = glob.glob('/home/CCCma/L2/input_data/ERA5/surface/*'+yymmdd+'*.nc')[0]

# add Netcdf files to manifest file
input_manifest.add_files(filessf)
input_manifest.add_files(file_lev)
input_manifest.add_files(file_sfc)

# write input manifest files to path
input_manifest.write(out_path='/home/jbrendecke/L2/data/')
#input_manifest.write(out_path='/home/CCCma/L2/data/')






