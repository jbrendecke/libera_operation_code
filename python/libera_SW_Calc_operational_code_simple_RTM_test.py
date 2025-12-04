#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sept 10 12:45:19 2024
#load input for CCCma with the inclusion of clouds from CERES SSF running the 
RTM for each footprint
First iteration for Libera operation, no CCCma to more easily set things up with docker/AWS
@author: jbrendecke
"""

import numpy as np
import subprocess
import os
import xarray as xr
from metpy.calc import specific_humidity_from_dewpoint, relative_humidity_from_dewpoint
from metpy.units import units
from libera_utils import Manifest, smart_open, DataProductConfig
import glob
import logging
import warnings
warnings.filterwarnings('ignore')


def set_jday(iyr, imo, iday):
#    ; Determine Julian Day from Calendar Day
#    ; INPUT:
#    ;; ------------------------------
#    ;; iyr = 4-digit year (INTEGER)
#    ;; imo = 2-digit month (INTEGER)
#    ;; iday = 2-digit day (INTEGER)
#    ;; 
#    ;; OUTPUT:
#    ;; ------------------------------
#    ;; iddd = Julian Day (0-365) (INTEGER)

    idymon = [[0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334],  # non-leap year
              [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]]  # leap year

    leap = iyr - (iyr // 4) * 4

    if leap == 0:  # Leap Year
        irow = 1
    else:  # Non-leap Year
        irow = 0

    iddd = idymon[irow][imo-1] + iday
    return iddd
#//////////////////////////////////////////////////////////////////////////////

def read_manifest(input_manifest_file):
    #//////////////////////////////////////
    # Reads in json input manifest file
    # and returns xarray dataset sets for
    # SSF properties, ERA5 levels, and ERA5 surface
    #/////////////////////////////////////
    
    input_manifest = Manifest.from_file(input_manifest_file)
    
    #put filename and datasets in manifest files into dictionary
    all_data={}
    for i, file_info in enumerate(input_manifest.files):
        try: 
            with smart_open(file_info.filename) as file_handle:
                dataset = xr.open_dataset(file_handle)
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
        
    return ds_ssf, ds_era5lev, ds_era5sfc
#//////////////////////////////////////////////////////////////////////////////

def mask_era5_sfc(ds_sfc, ds_cer):
    #////////////////////////////////////
    # identifies era5 grids that are closest SSF location/times
    # and sets any era5 grid that are not within .2 degrees and 30 minutes
    # to Nan values
    # Should reduce processing time of era5 profiles
    #
    # Inputs:
    #   xarray dataset for ERA5 surface
    #       includes: surface pressure, surface goepotential, 2m temperature, 2m dew point
    #   xarray dataset of footprint
    #       must include time, latitude, and longitude
    #
    # Output:
    #   xarray dataset for masked ERA5 surface dataset
    #       indcludes: pressure, height(km), temp, rh, specific humidity, ozone
    #       dims: time, lev, lat, lon
    #/////////////////////////////////////
    

    MAX_TIME_DIFF = np.timedelta64(31, 'm')  # e.g., 31 minutes
    MAX_LAT_DIFF = 0.2  # e.g., 0.2 degrees
    MAX_LON_DIFF = 0.2  # e.g., 0.2 degrees
    
    # Initialize a mask with all False (no grid points are close yet)
    mask_shape = (len(ds_sfc['valid_time']), len(ds_sfc['latitude']), len(ds_sfc['longitude']))
    proximity_mask = np.full(mask_shape, False)
    
    # Get the indices of the grid that correspond to the satellite locations
    time_indices = ds_sfc.indexes['valid_time'].get_indexer(ds_cer['time'].values, method='nearest')
    lat_indices = ds_sfc.indexes['latitude'].get_indexer(ds_cer['lat'].values, method='nearest')
    lon_indices = ds_sfc.indexes['longitude'].get_indexer(ds_cer['lon'].values, method='nearest')
    
    # Iterate over all satellite samples
    for i in range(len(ds_cer['time'])):
        t_idx = time_indices[i]
        lat_idx = lat_indices[i]
        lon_idx = lon_indices[i]
            
        # Calculate the actual difference between the satellite point and the nearest grid point
        time_diff = np.abs(ds_cer['time'].values[i] - ds_sfc['valid_time'].values[t_idx])
        lat_diff = np.abs(ds_cer['lat'].values[i] - ds_sfc['latitude'].values[lat_idx])
        lon_diff = np.abs(ds_cer['lon'].values[i] - ds_sfc['longitude'].values[lon_idx])
    
        # If all differences are within the defined tolerance, mark the grid point as True
        if (time_diff < MAX_TIME_DIFF and 
            lat_diff < MAX_LAT_DIFF and 
            lon_diff < MAX_LON_DIFF):
            
            proximity_mask[t_idx, lat_idx, lon_idx] = True
    
    
    # Convert the NumPy mask into an xarray DataArray
    proximity_da = xr.DataArray(
        proximity_mask,
        coords={'valid_time': ds_sfc['valid_time'], 'latitude': ds_sfc['latitude'], 'longitude': ds_sfc['longitude']},
        dims=['valid_time', 'latitude', 'longitude']
    )
    
    # Apply the inverse mask to set unwanted points to NaN
    ds_sfc_masked = ds_sfc.where(proximity_da, drop=False)
    return ds_sfc_masked
#//////////////////////////////////////////////////////////////////////////////
    
def combine_era5_reanalysis(ds_lev, ds_sfc, NLEV):
    #////////////////////////////////////
    # combines era5 surface product with era5 levels product to produce
    # profiles from surface to TOA for pressure, hgt, temp, qv/rh, o3
    # profiles are interpolated to NLEV 
    #
    # Inputs:
    #   ds_lev: 
    #       xarray data set for ERA5 levels
    #       includes: Levels, geopotential heights, temp, specific & relative humidity, ozone
    #   ds_sfc:
    #       xarray data set for ERA5 surface
    #       includes: surface pressure, surface goepotential, 2m temperature, 2m dew point
    #   NLEV: 
    #       interger for number of levels to interpolate to 
    #
    # Output:
    #   xarray data set for combined product
    #       indcludes: pressure, height(km), temp, rh, specific humidity, ozone
    #       dims: time, lev, lat, lon
    #/////////////////////////////////////
    
    # Get altitude from geopotential
    earth_radius = 6378137 
    ds_lev['z'] = ds_lev['z']/9.8
    ds_sfc['z'] = ds_sfc['z']/9.8
    ds_lev['z'] = (earth_radius * ds_lev['z'] / (earth_radius - ds_lev['z']))/1000
    ds_sfc['z'] = (earth_radius * ds_sfc['z'] / (earth_radius - ds_sfc['z']))/1000
    ds_sfc['z'] = ds_sfc['z'].round(decimals=5)
    ds_lev['z'] = ds_lev['z'].round(decimals=5)

    #determine specifc humidity from dewpoint/pressure and it to sfc dataset
    q_sfc = specific_humidity_from_dewpoint(ds_sfc['sp']/100 * units.hPa, ds_sfc['d2m'] * units.kelvin)
    ds_sfc['q'] = q_sfc
    rh_sfc = relative_humidity_from_dewpoint(ds_sfc['t2m'] * units.kelvin, ds_sfc['d2m'] * units.kelvin)*100
    ds_sfc['rh'] = rh_sfc
   
    #make sure variables are not negative/ within range
    ds_lev['r'] = ds_lev['r'].clip(min=0, max=100)
    ds_lev['q'] = ds_lev['q'].where(ds_lev['q'] >= 0, 0)
    ds_lev['o3'] = ds_lev['o3'].where(ds_lev['o3'] >= 0, 0)
    ds_lev['t'] = ds_lev['t'].where(ds_lev['t'] >= 0, 300)
    
    #set variables equal to NaN where pressure levels in levels dataset are higher than surface pressure(i.e. below ground)
    sp = ds_sfc['sp']/100
    expanded_sp = sp.expand_dims({'pressure_level': ds_lev.pressure_level}, axis=1)
    mask = ds_lev.pressure_level > expanded_sp
    ds_lev = ds_lev.where(~mask)
        
    #extract variables from era5 
    time_era5 = ds_lev.valid_time.data
    lat_era5 = ds_lev.latitude.data 
    lon_era5 = ds_lev.longitude.data  
    p_lev = ds_lev['pressure_level'].data
    z_lev = ds_lev['z'].data
    t_lev = ds_lev['t'].data
    q_lev = ds_lev['q'].data
    rh_lev = ds_lev['r'].data
    o3_lev = ds_lev['o3'].data
    p_sfc = ds_sfc['sp'].data/100
    z_sfc = ds_sfc['z'].data
    t_sfc = ds_sfc['t2m'].data
    rh_sfc = ds_sfc['rh'].data
    q_sfc = ds_sfc['q'].data
    ntime = len(time_era5)
    nlat = len(lat_era5)
    nlon = len(lon_era5)
    
    presx = np.zeros((ntime,NLEV,nlat,nlon))
    hgtx = np.zeros((ntime,NLEV,nlat,nlon))
    tempx = np.zeros((ntime,NLEV,nlat,nlon))
    qvx = np.zeros((ntime,NLEV,nlat,nlon))
    rhx = np.zeros((ntime,NLEV,nlat,nlon))
    o3x = np.zeros((ntime,NLEV,nlat,nlon))
    for itime in range(ntime):
        for ilat in range(nlat):
            for ilon in range(nlon):
                # only process where temp surface values is true
                #(False when not within time/distance range of footprint)
                tempnan = ~np.isnan(t_sfc[itime, ilat, ilon])               
                if tempnan:
                
                    indz = np.where((~np.isnan(z_lev[itime,:, ilat, ilon])) & 
                                    (z_lev[itime, :, ilat, ilon] >= 0))
        
                    #pull out variables for time-step
                    p_l = p_lev[indz]
                    z_l = np.squeeze(z_lev[itime,indz,ilat,ilon])
                    t_l = np.squeeze(t_lev[itime,indz,ilat,ilon])
                    q_l = np.squeeze(q_lev[itime,indz,ilat,ilon])
                    rh_l = np.squeeze(rh_lev[itime,indz,ilat,ilon])
                    o3_l = np.squeeze(o3_lev[itime,indz,ilat,ilon])
                    p_s = p_sfc[itime,ilat,ilon]
                    z_s = z_sfc[itime,ilat,ilon]
                    t_s = t_sfc[itime,ilat,ilon]
                    q_s = q_sfc[itime,ilat,ilon]
                    rh_s = rh_sfc[itime,ilat,ilon]
                    
                    #check height/pressure increases/decreases with height
                    if z_s >= z_l[0]: #sometimes surface altitude is higher than level-1 height
                        z_s = z_l[0] - ((z_s -z_l[0])+.03) # make surface 3 meters lower than lowest level
                    if p_s <= p_l[0]:
                        p_s = p_l[0] + ((p_l[0] - p_s)+2) # add 2 hps
                        
                    if p_s <= p_l[0]:
                        print('p:', itime,ilat,ilon)
                    if z_s >= z_l[0]: 
                        print('z lower:', itime,ilat,ilon, z_s - z_l[0])
                    
                    
                    #combine surface level and pressure level
                    press = np.concatenate(([p_s], p_l))
                    hts = np.concatenate(([z_s], z_l))
                    temps = np.concatenate(([t_s], t_l))
                    qvs = np.concatenate(([q_s], q_l))
                    rhs = np.concatenate(([rh_s], rh_l))
                    o3s = np.concatenate(([o3_l[0]], o3_l))
                    
                    #interpolate to 41 layers                    
                    x = np.linspace(0, 1, len(press))  # Example x-coordinates, adjust if needed.
                    new_x = np.linspace(0, 1, NLEV)  # New x-coordinates for interpolation
                    
                    presx[itime, :, ilat, ilon] = np.interp(new_x, x, press)  
                    hgtx[itime, :, ilat, ilon] = np.interp(new_x, x, hts)
                    rhx[itime, :, ilat, ilon] = np.interp(new_x, x, rhs)
                    tempx[itime, :, ilat, ilon] = np.interp(new_x, x, temps)
                    qvx[itime, :, ilat, ilon] = np.interp(new_x, x, qvs)
                    o3x[itime, :, ilat, ilon] = np.interp(new_x, x, o3s)
                    
    lev_arbt = np.linspace(1, NLEV, NLEV, dtype=int)
    dims = ('time', 'lev', 'lat', 'lon')
    
    ds_new = xr.Dataset(
        data_vars= {
            'p' : (dims, presx, {'units' : 'mb'}),
            'z' : (dims, hgtx, {'units' : 'km'}),
            'rh' : (dims, rhx, {'units' : '%'}),
            't' : (dims, tempx, {'units' : 'kelvin'}),
            'q' : (dims, qvx, {'units' : 'kg/kg'}),
            'o3' : (dims, o3x, {'units' : 'g/kg'})
            },
            coords={
                'time' : time_era5,
                'lev' : lev_arbt,
                'lat' : lat_era5,
                'lon' : lon_era5
                }
            )
    
    return ds_new   
#//////////////////////////////////////////////////////////////////////////////

def cloud_profile_constant(hgt, pres, re, lwp, CB_pres, CT_pres):
    #Coputes constant RC and LWC with height
    #Inputs:
        #hgt: profile of heights from surface to TOA [km]
        #pres: profile of pressure form surface to TOA [hPa]
        #re: single value, either liquid or ice particle size [um]
        #LWP: total LWP or IWP [g m-2]
        #CB_pres: Cloud Base pressure level [hPa]
        #CT_pres: Cloud Top pressure level [hPa]
    #Outputs:
        #LWC_profile: Profile of LWC or IWC 
        #Rc_profile: Profile of Rc or Ri
        
    Re_profile_lev = np.full(len(hgt), np.nan)  
    WC_profile_lev = np.full(len(hgt), np.nan)
    Re_profile_lay = np.full(len(hgt)-1, np.nan)  
    WC_profile_lay = np.full(len(hgt)-1, np.nan)
    
    cldlevb = np.argmin(abs(CB_pres-pres))
    cldlevt = np.argmin(abs(CT_pres-pres))
    
    #lwp = (tau * 5 * re*1e-6 *100**3) / 9    
    
    if (cldlevb == cldlevt):
        cldthick = (hgt[cldlevt+1]-hgt[cldlevb])*1000 
        Re_profile_lev[cldlevb:cldlevt+2] = re
        WC_profile_lev[cldlevb:cldlevt+2] = lwp / cldthick
    else:
        cldthick = (hgt[cldlevt] - hgt[cldlevb])*1000
        Re_profile_lev[cldlevb:cldlevt+1] = re
        WC_profile_lev[cldlevb:cldlevt+1] = lwp / cldthick
    
    #average into layer average
    if cldlevb==cldlevt:
        cldlevt +=1
    for i in range(cldlevb,cldlevt): 
        WC_profile_lay[i] = (WC_profile_lev[i+1]+WC_profile_lev[i])/2
        Re_profile_lay[i] = (Re_profile_lev[i+1]+Re_profile_lev[i])/2
        
    #check that LWC profile equals LWP
    check=0
    for i in range(len(hgt)-1):
        check += (hgt[i+1]-hgt[i])*WC_profile_lay[i]*1000
    if abs(check -lwp) > 0.1:
        logging.error("CLOUD PROFILE NOT CORRECT, LWP OFF: {CHECK} NOT EQUAL TO {LWP)")
                
    return WC_profile_lay, Re_profile_lay
#//////////////////////////////////////////////////////////////////////////////

def cloud_fraction_profile(pres, CB_pres, CT_pres, CF):
    #creates a profile for cloud fraction base on both liquid and ice properties
    CF_profile_lev = np.full(len(pres), np.nan) 
    CF_profile_lay = np.full(len(pres)-1, np.nan) 

    cldlevb = np.argmin(abs(CB_pres-pres))
    cldlevt = np.argmin(abs(CT_pres-pres))
    if CF >= 0.95:
        CF = 1.0
        
    if (cldlevb == cldlevt):
         CF_profile_lev[cldlevb:cldlevb+2] = CF
    else:
        CF_profile_lev[cldlevb:cldlevt+1] = CF
        
    if cldlevb==cldlevt:
        cldlevt +=1
    for i in range(cldlevb,cldlevt): 
        CF_profile_lay[i] = (CF_profile_lev[i+1]+CF_profile_lev[i])/2
        
    return CF_profile_lay
#//////////////////////////////////////////////////////////////////////////////

def get_cloud_profile(dataset_cloud, hgt, pres, lat_fp, lon_fp, time_fp, lyr):
    # determine cloud profiles bases on CERES SSF Input
    # Inputs: 
    # CERES Nc file dataset
    # Reanalysis heights and pressure levels
    # Footprint latitude, longitdue, and time
    # number of atmospheric levels
    #//////////////////////////////////////////////////////
    nlay = lyr -1
    
    dataset_cld = dataset_cloud.sel(time = time_fp, method='nearest' )
    cf_overlap = (dataset_cld['Clear_layer_overlap_percent_coverages'].data)#Columns: Clear | Bottom CF | Top CF | Overlap % (always=0)
    ct_ssf = (dataset_cld['Mean_cloud_top_pressure_for_cloud_layer'].data)
    cb_ssf = (dataset_cld['Mean_cloud_base_pressure_for_cloud_layer'].data)
    rc_21_ssf = (dataset_cld['Mean_water_particle_radius_for_cloud_layer__2_1_'].data)
    rc_37_ssf = (dataset_cld['Mean_water_particle_radius_for_cloud_layer__3_7_'].data)
    ric_21_ssf = (dataset_cld['Mean_ice_particle_effective_radius_for_cloud_layer__2_1_'].data)
    ric_37_ssf = (dataset_cld['Mean_ice_particle_effective_radius_for_cloud_layer__3_7_'].data)
    cod_ssf = (dataset_cld['Mean_visible_optical_depth_for_cloud_layer'].data)
          
    #determine cloud properties using MERRA2 Hieghts and CERES SYN Properties
    #cf_ssf = 100 - cf_ssf
    cf_btm = cf_overlap[1]/100
    cf_top = cf_overlap[2]/100
    cf_tot = cf_btm + cf_top
    if np.isnan(cf_tot):
        cf_tot = 1- (cf_overlap[0]/100)
    if np.isnan(cf_tot):
        cf_tot = 0
    
    #Rc_2.1um not report for small COD values, can subsitute for 3.7um channel
    if ~np.isnan(rc_21_ssf[0]):
        rcb = rc_21_ssf[0]
    elif ~np.isnan(rc_37_ssf[0]):
        rcb = rc_37_ssf[0]
    else:
        rcb = np.nan
    if ~np.isnan(rc_21_ssf[1]):
        rct = rc_21_ssf[1]
    elif ~np.isnan(rc_37_ssf[1]):
        rct = rc_37_ssf[1]
    else:
        rct = np.nan
        
    if ~np.isnan(ric_21_ssf[0]):
        ricb = ric_21_ssf[0]
    elif ~np.isnan(rc_37_ssf[0]):
        ricb = ric_37_ssf[0]
    else:
        ricb = np.nan
    if ~np.isnan(ric_21_ssf[1]):
        rict = ric_21_ssf[1]
    elif ~np.isnan(rc_37_ssf[1]):
        rict = ric_37_ssf[1]
    else:
        rict = np.nan
    
    if (cf_tot > 0.01):
        #Upper Layer
        if (cf_top > 0.05) & (~np.isnan(rct)) | (~np.isnan(rict)): 
            cf_hg = cloud_fraction_profile(pres, cb_ssf[1], ct_ssf[1], cf_top)
            if ~np.isnan(rct):  
                lwp_top = (cod_ssf[1] * 5 * rct*1e-6 *100**3) / 9  
                lwc_hg, rc_hg = cloud_profile_constant(hgt, pres, rct, lwp_top, 
                                                         cb_ssf[1], ct_ssf[1])
            else: 
                lwc_hg = np.full(nlay, np.nan)
                rc_hg = np.full(nlay, np.nan)
                
            if ~np.isnan(rict):
                iwp_top = (cod_ssf[1] * 2 * rict*1e-6 * 96**3) / 3  
                iwc_hg, ri_hg = cloud_profile_constant(hgt, pres, rict, iwp_top, 
                                                         cb_ssf[1], ct_ssf[1])
            else: 
                iwc_hg = np.full(nlay, np.nan)
                ri_hg = np.full(nlay, np.nan)
        else:
            lwc_hg = np.full(nlay, np.nan); iwc_hg = np.full(nlay, np.nan)
            rc_hg = np.full(nlay, np.nan); ri_hg = np.full(nlay, np.nan)
            cf_hg = np.full(nlay, np.nan)


        #Low Layer
        if (cf_btm > 0.05) & (~np.isnan(rcb)) | (~np.isnan(ricb)): 
            cf_lw = cloud_fraction_profile(pres, cb_ssf[0], ct_ssf[0], cf_btm)
            if ~np.isnan(rcb):  
                lwp_btm = (cod_ssf[0] * 5 * rcb*1e-6 *100**3) / 9  
                lwc_lw, rc_lw = cloud_profile_constant(hgt, pres, rcb, lwp_btm, 
                                                         cb_ssf[0], ct_ssf[0])
            else: 
                lwc_lw = np.full(nlay, np.nan)
                rc_lw = np.full(nlay, np.nan)
                
            if ~np.isnan(ricb):
                iwp_btm = (cod_ssf[0] * 2 * ricb*1e-6 * 96**3) / 3  
                iwc_lw, ri_lw = cloud_profile_constant(hgt, pres, ricb, iwp_btm, 
                                                         cb_ssf[0], ct_ssf[0])
            else: 
                iwc_lw = np.full(nlay, np.nan)
                ri_lw = np.full(nlay, np.nan)
        else:
            lwc_lw = np.full(nlay, np.nan); iwc_lw = np.full(nlay, np.nan)
            rc_lw = np.full(nlay, np.nan); ri_lw = np.full(nlay, np.nan)
            cf_lw = np.full(nlay, np.nan)
        
        
        #####
        #still some weird cases where bottom liq/ice are record for either
        #top/bottom layer or where top/bottom layer thickness overlaps with
        #bottom/top layer thickness then mean cf is weird and you have 
        #liq/ice in the same layer
        #####
                                   
        #combine profiles into one
        rcm = np.nanmean((rc_lw, rc_hg), axis =0)
        lwcm = np.nanmean((lwc_lw, lwc_hg), axis=0)
        rim = np.nanmean((ri_lw, ri_hg), axis=0)
        iwcm = np.nanmean((iwc_lw, iwc_hg), axis=0)
        cfm = np.nanmean((cf_lw, cf_hg), axis =0)
        #replace NaNs
        rcm[np.isnan(rcm)] = 0.0
        lwcm[np.isnan(lwcm)] = 0.0
        rim[np.isnan(rim)] = 0.0
        iwcm[np.isnan(iwcm)] = 0.0
        cfm[np.isnan(cfm)] = 0.0
        
    else:
        rcm = np.zeros(nlay)
        lwcm = np.zeros(nlay)
        rim = np.zeros(nlay)
        iwcm = np.zeros(nlay)
        cfm = np.zeros(nlay)
        
    return rcm, lwcm, rim, iwcm, cfm, cf_tot
#//////////////////////////////////////////////////////////////////////////////

def simple_RTM(sza, cf, COD, AOD):
    # very simple RTM to replace CCCma for now & keep everything in python
    nlen = len(sza)
    S0 = 1361 * np.cos(sza * (np.pi/180))
    w0 = .85
    asym =.6
    
    swd_dir = np.zeros(nlen); swd_dif = np.zeros(nlen)
    swd_dir_clr = np.zeros(nlen); swd_dif_clr = np.zeros(nlen)
    toa = np.zeros(nlen); toa_clr = np.zeros(nlen)
    for i in range(nlen):
        if cf[i] > .5:
            tau = max(COD[i]) + 5
        else:
            tau = AOD[i]
            
        swd_dir[i] = S0[i] * np.exp(-1 * tau * (1/np.cos(sza[i]*np.pi/180))) 
        swd_dif[i] = S0[i] * w0 * (1 - np.exp(-1 * tau * (1/np.cos(sza[i]*np.pi/180)))) * asym
        
        swd = swd_dir[i] +swd_dif[i]
        
        if cf[i] > .5:
            albedo = .85
            toa = S0[i] * albedo
        else:
            albedo = .5
            toa = (swd * albedo) * np.exp(-1 * tau) 
        
        if cf[i] > .5:
            tau = AOD[i]
            swd_dir_clr[i] = S0[i] * np.exp(-1 * tau * (1/np.cos(sza[i]*np.pi/180))) 
            swd_dif_clr[i] = S0[i] * w0 * (1 - np.exp(-1 * tau * (1/np.cos(sza[i]*np.pi/180)))) * asym
            
            swd = swd_dir[i] +swd_dif[i]
            albedo = .5
            toa_clr[i] = (swd * albedo) * np.exp(-1 * tau)
        else:
            swd_dir_clr[i] = swd_dir[i]
            swd_dif_clr[i] = swd_dif[i]
            toa_clr[i] = toa
    
    swd_toa_nir = 1361 *np.cos(sza *(np.pi/180)) * .55
    swd_toa_vis = 1361 *np.cos(sza *(np.pi/180)) * .45
    swu_toa_all_nir = toa *.55
    swu_toa_all_vis = toa * .45
    swu_toa_clr_nir = toa_clr *.55
    swu_toa_clr_vis = toa_clr * .45
    swd_sfc_all_nir = (swd_dir + swd_dif) * .55
    swd_sfc_all_vis = (swd_dir + swd_dif) * .45
    swd_sfc_clr_nir = (swd_dir_clr + swd_dif_clr) * .55
    swd_sfc_clr_vis = (swd_dir_clr + swd_dif_clr) * .45
    swu_sfc_all_nir = swd_sfc_all_nir * .5
    swu_sfc_all_vis = swd_sfc_all_vis * .5
    swu_sfc_clr_nir = swd_sfc_all_nir * .5
    swu_sfc_clr_vis = swd_sfc_all_vis * .5
    direct_sfc_all = swd_dir
    diffuse_sfc_all = swd_dif
    direct_sfc_clr = swd_dir_clr
    diffuse_sfc_clr = swd_dif_clr
   
    swu_toa_all_nir = (swu_toa_all_nir * cf) + (swu_toa_clr_nir * (1-cf))
    swu_toa_all_vis = (swu_toa_all_vis * cf) + (swu_toa_clr_vis * (1-cf))
    swd_sfc_all_nir = (swd_sfc_all_nir * cf) + (swd_sfc_clr_nir * (1-cf))
    swd_sfc_all_vis = (swd_sfc_all_vis * cf) + (swd_sfc_clr_vis * (1-cf))
    swu_sfc_all_nir = (swu_sfc_all_nir * cf) + (swu_sfc_clr_nir * (1-cf))
    swu_sfc_all_vis = (swu_sfc_all_vis * cf) + (swu_sfc_clr_vis * (1-cf))
    direct_sfc_all = (direct_sfc_all * cf) + (direct_sfc_clr * (1-cf))     
    diffuse_sfc_all = (diffuse_sfc_all * cf) + (diffuse_sfc_clr * (1-cf))
    
    rtm_dict = {
        'swd_toa_nir': swd_toa_nir, 'swd_toa_vis': swd_toa_vis,
        'swu_toa_all_nir': swu_toa_all_nir, 'swu_toa_all_vis': swu_toa_all_vis,
	    'swu_toa_clr_nir': swu_toa_clr_nir, 'swu_toa_clr_vis': swu_toa_clr_vis,
	    'swd_sfc_all_nir': swd_sfc_all_nir, 'swd_sfc_all_vis': swd_sfc_all_vis,
	    'swd_sfc_clr_nir': swd_sfc_clr_nir, 'swd_sfc_clr_vis': swd_sfc_clr_vis,
	    'swu_sfc_all_nir': swu_sfc_all_nir, 'swu_sfc_all_vis': swu_sfc_all_vis,
	    'swu_sfc_clr_nir': swu_sfc_clr_nir, 'swu_sfc_clr_vis': swu_sfc_clr_vis,
	    'direct_sfc_all': direct_sfc_all, 'diffuse_sfc_all': diffuse_sfc_all,
	    'direct_sfc_clr': direct_sfc_clr, 'diffuse_sfc_clr': diffuse_sfc_clr
	}
    return rtm_dict    
#//////////////////////////////////////////////////////////////////////////////


#def L2_calc_main(iyr,imon, iday):
xx=0
if xx==0: 
    iyr=2018
    imon=3
    iday=22
    # Main function 
    #//////////////////////////////////////////////////////////////////////////////
    NLEV = 41
    NLAY = NLEV - 1
    
    
    iyr = int(iyr)
    imon = int(imon)
    iday = int(iday)
    
    smonth0 = f"{imon:02}"
    sday0 = f"{iday:02}"
    yymmdd = f"{iyr}{smonth0}{sday0}"
    
    jday_int = set_jday(int(yymmdd[0:4]), int(yymmdd[4:6]), int(yymmdd[6:]))   
    if (int(yymmdd[4:6]) >=3) & (int(yymmdd[4:6]) <= 8):
        season = 1
    else:
        season = 0
        
        
        
        
    #//////////////////LOAD SSF Non-cloud properties///////////////////////////////
    path_ceres ='/home/CCCma/L2/input_data/CERES/'
    filessf = glob.glob(path_ceres+'*'+yymmdd+'*.nc')
    ds_ssf = xr.open_dataset(filessf[0])
    
    #only testing 1 hour for now
    ds_ssf = ds_ssf.isel(time=slice(0,99400)) #1hr
    
    time_ssf=np.array(ds_ssf.time).copy()
    lon_ssf = np.array(ds_ssf.lon)
    lat_ssf = np.array(ds_ssf.lat)
    cod_ssf = (ds_ssf['Mean_visible_optical_depth_for_cloud_layer'].data)
    sza_ssf = (ds_ssf['CERES_solar_zenith_at_surface'].data)
    aod_land = (ds_ssf['PSF_wtd_MOD04_deep_blue_aerosol_optical_depth_land__0_550_'].data)
    aod_ocean = (ds_ssf['PSF_wtd_MOD04_effective_optical_depth_average_ocean__0_550_'].data)
    igbp = (ds_ssf['Surface_type_index'].data)
    ntime=len(time_ssf)
    
    #time in CERES file is number of days from 1970 and because of that, datatime is wrong, this fixes it
    correct_date = np.datetime64(f'{iyr}-{smonth0}-{sday0}T00:00:00')
    time_ssf = correct_date + (time_ssf - time_ssf[0])
    time_check = str(time_ssf[0])
    if (time_check[:4] != str(iyr)) & (time_check[5:7] != smonth0) & (time_check[8:10] != sday0):
        raise ValueError('SSF date does not match input date')
    
    #make longitude from 0-360
    new_long = np.mod(ds_ssf['lon'], 360) 
    ds_ssf = ds_ssf.assign_coords(lon=new_long)
    
    #right now assume most dominant scene type is used/includes snow
    igbp = igbp[:,0]
    #adjust to CCCma ID
    igbp = np.int_(igbp + 40)
    
    
    
    
    #////////////////////////ERA5 Atmospheric Profiles////////////////////////////
    file_lev = glob.glob('/home/CCCma/L2/input_data/ERA5/levels/*level*'+yymmdd+'*.nc')[0]
    file_sfc = glob.glob('/home/CCCma/L2/input_data/ERA5/surface/*surface*'+yymmdd+'*.nc')[0]
    
    if not os.path.isfile(file_lev) | os.path.isfile(file_sfc):
        logging.error('ERA5 File Not Found')
        raise ValueError('ERA5 File Not Found')
    
    #open dataset and combine with function
    ds_lev = xr.open_dataset(file_lev)
    ds_sfc = xr.open_dataset(file_sfc)
    
    
    #for first hour     
    ds_lev = ds_lev.isel(valid_time=[0])
    ds_sfc = ds_sfc.isel(valid_time=[0])
    
    #process era5 profiles
    ds_sfc = mask_era5_sfc(ds_sfc, ds_ssf)
    ds_era5 = combine_era5_reanalysis(ds_lev, ds_sfc, NLEV)
    
    
    #///////////////////////Output Variables to be saved///////////////////////////
    swd_toa_vis = np.zeros(ntime); swu_toa_all_vis = np.zeros(ntime); swu_toa_clr_vis = np.zeros(ntime);
    swd_toa_nir = np.zeros(ntime); swu_toa_all_nir = np.zeros(ntime); swu_toa_clr_nir = np.zeros(ntime);
    swd_sfc_all_vis = np.zeros(ntime); swd_sfc_clr_vis = np.zeros(ntime);
    swu_sfc_all_vis = np.zeros(ntime); swu_sfc_clr_vis = np.zeros(ntime);
    swd_sfc_all_nir = np.zeros(ntime); swd_sfc_clr_nir = np.zeros(ntime);
    swu_sfc_all_nir = np.zeros(ntime); swu_sfc_clr_nir = np.zeros(ntime);
    direct_sfc_all = np.zeros(ntime); direct_sfc_clr = np.zeros(ntime);
    diffuse_sfc_all = np.zeros(ntime); diffuse_sfc_clr = np.zeros(ntime);
    
    
    #//////////////////RUN CKD CODE FOR EACH SSF////////////
    #select out daylight SSFs
    ind_day = np.where(sza_ssf <= 89.5)[0]
    max0 = 20000 #max number of ilg loops that is ran in fortran
    nloop = int(np.floor(len(ind_day) / max0))
    remainder = len(ind_day) % max0
  
    # # # # # Run Loop for each CCCma ilg loop # # # # #
    for ilg in range(nloop+1):
       
        #determine footprints with daylight
        if ilg == nloop:
            #for last remaining footprints
            indi = ind_day[len(ind_day)-remainder:]
        else:
            # most footprints; length equals max0
            indi = ind_day[max0*ilg:max0*(ilg+1)]
        nfoots = len(indi)
        
        
        #save properties for RTM Input
        cod = cod_ssf[indi]
        sza = sza_ssf[indi]
        igbptyp = igbp[indi]
        
        aod = np.zeros(nfoots)
        aodtype = np.zeros(nfoots, dtype=int)
        ind20km = np.zeros(nfoots, dtype=int)
        presx = np.zeros((nfoots, NLEV))
        htx = np.zeros((nfoots, NLEV))
        rhx = np.zeros((nfoots, NLEV))
        tempx = np.zeros((nfoots, NLAY))
        qvx = np.zeros((nfoots, NLAY))
        o3x = np.zeros((nfoots, NLAY))
        relx = np.zeros((nfoots, NLAY))
        reix = np.zeros((nfoots, NLAY))
        lwcx = np.zeros((nfoots, NLAY))
        iwcx = np.zeros((nfoots, NLAY))
        cf_px = np.zeros((nfoots, NLAY))
        cft = np.zeros(nfoots)
    
        # # # # # RUN LOOP FOR EACH FOOTPRINT # # # # #
        kk=0
        for ii in indi:
            # Determine Aerosol Properties
            #AOD
            if igbptyp[kk] ==57:
                aod[kk] = aod_ocean[ii]
            else:
                aod[kk] = aod_land[ii]
            if np.isnan(aod[kk]):
                aod[kk] = 0.10
            #AOD TYPE
            if igbptyp[kk] == 57:
                aodtype[kk] = 1
            elif igbptyp[kk] == 56:
                aodtype[kk] = 4
            elif abs(lat_ssf[ii]) > 70:
                aodtype[kk] = 5
            else:
                aodtype[kk] = 2
            
            #/////////////Save ERA5 Profile/////////////
            #Pressure, Hieght, RH should have 41 levels
            #Temperature, QV, O3 should have 40 levels
            ds_era5ii = ds_era5.sel(time = time_ssf[ii], lat=lat_ssf[ii], 
                                    lon = lon_ssf[ii], method='nearest')
            
            pres = ds_era5ii['p'].data
            ht = ds_era5ii['z'].data
            temp = ds_era5ii['t'].data
            qv = ds_era5ii['q'].data
            rh = ds_era5ii['rh'].data
            o3 = ds_era5ii['o3'].data
            
            if ((np.any(np.isnan(pres))) | (np.any(np.isnan(ht))) | (np.any(np.isnan(temp))) | 
                (np.any(np.isnan(qv))) | (np.any(np.isnan(rh))) | (np.any(np.isnan(o3)))):
                raise ValueError('Nan values present in ERA5 profile extraction')
            
            #these variables need layer average
            tempm = np.zeros(NLAY)
            qvm = np.zeros(NLAY)
            o3m = np.zeros(NLAY)
            for zz in range(1, NLEV):
                tempm[zz-1] = (temp[zz-1]+temp[zz])/2
                qvm[zz-1] = (qv[zz-1]+qv[zz])/2
                o3m[zz-1] = (o3[zz-1]+o3[zz])/2
            
            #identify index closest to 20km altitude
            ind20km[kk] = np.argmin(abs(20 - ht))
            
            #/////////////Save Cloud Profile///////////
            #40 levels of Re liq, Re ice, LWC & IWC, and CF
            rc, lwc, ri, iwc, cf_p, cf_tot = get_cloud_profile(ds_ssf, ht, pres, 
                                                               lat_ssf[ii], lon_ssf[ii],
                                                               time_ssf[ii], NLEV)
            
            if ((np.any(np.isnan(rc))) | (np.any(np.isnan(ri))) | (np.any(np.isnan(lwc))) | 
                (np.any(np.isnan(iwc))) | (np.any(np.isnan(cf_p))) | (np.isnan(cf_tot))):
                raise ValueError('Nan values present in cloud profile')
                
            relx[kk, :] = np.flip(rc)
            reix[kk, :] = np.flip(ri)
            lwcx[kk, :] = np.flip(lwc)
            iwcx[kk, :] = np.flip(iwc)
            cf_px[kk, :] = np.flip(cf_p)
            cft[kk] = cf_tot
            
            presx[kk, :] = np.flip(pres)
            htx[kk, :] = np.flip(ht)
            tempx[kk, :] = np.flip(tempm)
            qvx[kk, :] = np.flip(qvm)
            rhx[kk, :] = np.flip(rh)
            o3x[kk, :] = np.flip(o3m)
            #/////////////
            kk += 1
    
        # # # # # Run Simple RTM # # # # #
        rtm_dic= simple_RTM(sza, cft, cod, aod)
        
        swd_toa_nir[indi] = rtm_dic['swd_toa_nir']
        swd_toa_vis[indi] = rtm_dic['swd_toa_vis']
        swu_toa_all_nir[indi] = rtm_dic['swu_toa_all_nir']
        swu_toa_all_vis[indi] = rtm_dic['swu_toa_all_vis']
        swu_toa_clr_nir[indi] = rtm_dic['swu_toa_clr_nir']
        swu_toa_clr_vis[indi] = rtm_dic['swu_toa_clr_vis']
        swd_sfc_all_nir[indi] = rtm_dic['swd_sfc_all_nir']
        swd_sfc_all_vis[indi] = rtm_dic['swd_sfc_all_vis']
        swd_sfc_clr_nir[indi] = rtm_dic['swd_sfc_clr_nir']
        swd_sfc_clr_vis[indi] = rtm_dic['swd_sfc_clr_vis']
        swu_sfc_all_nir[indi] = rtm_dic['swu_sfc_all_nir']
        swu_sfc_all_vis[indi] = rtm_dic['swu_sfc_all_vis']
        swu_sfc_clr_nir[indi] = rtm_dic['swu_sfc_clr_nir']
        swu_sfc_clr_vis[indi] = rtm_dic['swu_sfc_clr_vis']
        direct_sfc_all[indi] = rtm_dic['direct_sfc_all']
        diffuse_sfc_all[indi] = rtm_dic['diffuse_sfc_all']
        direct_sfc_clr[indi] = rtm_dic['direct_sfc_clr']
        diffuse_sfc_clr[indi] = rtm_dic['diffuse_sfc_clr']
        
        
        # # RUN CCCma
        # ctx=nfoots
        # os.chdir('/home/CCCma/L2/CCCma_code/ver10/')
        # pathi ='/home/CCCma/L2/CCCma_code/ver10/input/'
        # patho ='/home/CCCma/L2/CCCma_code/ver10/output/'
        
        # #WRITE INPUT FILE AND ADJUST UA_MAIN.F90
        # with open('UA_main.F90', 'r') as file:
        #     dmain = file.readlines()
    
        # dmain[4] = f'   parameter (lay = {NLAY}, lev = lay +1, ilg={ctx}, il1 = 1, il2 = {ctx}, modlay=34, nxloc=100)\n'
    
        # with open('UA_main.F90', 'w') as file:
        #     file.writelines(dmain)
        
        # with open('input.dat', 'w') as file:
        #     for k in range(nfoots):
        #         file.write(f"{sza[k]:11.4f} {aod[k]:11.4f} {igbptyp[k]:4d} {season:4d} {aodtype[k]:4d} {jday_int:5d} {cft[k]:7.4f}\n")
        #         for ilyr in range(NLAY):
        #             file.write(f"{htx[k,ilyr]:10.4f} {rhx[k,ilyr]:9.4f} {presx[k,ilyr]:9.4f} {tempx[k,ilyr]:9.4f}" 
        #                        f"{qvx[k,ilyr]:11.4e} {o3x[k,ilyr]:11.4e} {lwcx[k,ilyr]:9.4f} {iwcx[k,ilyr]:9.4f}"
        #                        f"{relx[k,ilyr]:9.4f} {reix[k,ilyr]:9.4f} {cf_px[k,ilyr]:9.4f}\n")
                    
        #         file.write(f"{htx[k,ilyr+1]:10.4f} {rhx[k,ilyr+1]:9.4f} {presx[k,ilyr+1]:9.4f}\n")
            
        # #RUN THE CCCMA
        # subprocess.run(['gfortran *.F90'], shell=True)
        # subprocess.run(['./a.out'], shell=True)
        # subprocess.run(['cp input.dat '+pathi+'input.dat'], shell=True)
        # subprocess.run(['cp output.dat '+patho+'output.dat'], shell=True)
    
        # #read CCCma output file
        # output_file = glob.glob(patho+'output.dat')
        # nlines=49
        # swd_toa_vis_ = np.zeros(nfoots)
        # swd_toa_nir_ = np.zeros(nfoots)
        # swu_toa_all_vis_ = np.zeros(nfoots)
        # swu_toa_all_nir_ = np.zeros(nfoots)
        # swu_toa_clr_vis_ = np.zeros(nfoots)
        # swu_toa_clr_nir_ = np.zeros(nfoots)
        # swd_sfc_all_vis_ = np.zeros(nfoots)
        # swd_sfc_all_nir_ = np.zeros(nfoots)
        # swd_sfc_clr_vis_ = np.zeros(nfoots)
        # swd_sfc_clr_nir_ = np.zeros(nfoots)
        # swu_sfc_all_vis_ = np.zeros(nfoots)
        # swu_sfc_all_nir_ = np.zeros(nfoots)
        # swu_sfc_clr_vis_ = np.zeros(nfoots)
        # swu_sfc_clr_nir_ = np.zeros(nfoots)
        # direct_sfc_all_ = np.zeros(nfoots)
        # diffuse_sfc_all_ = np.zeros(nfoots)
        # direct_sfc_clr_ = np.zeros(nfoots)
        # diffuse_sfc_clr_ = np.zeros(nfoots)
        
        # with open(output_file[0], 'r') as ofile:
        #     lines = ofile.readlines()
        #     for i in range(nfoots):
        #         swd_toa_vis_[i] =(float(lines[1+(i*nlines)].split()[3]))
        #         swd_toa_nir_[i] =(float(lines[1+(i*nlines)].split()[4]))
        #         swu_toa_all_vis_[i] =(float(lines[ind20km[i]+(i*nlines)+1].split()[1]))
        #         swu_toa_all_nir_[i] =(float(lines[ind20km[i]+(i*nlines)+1].split()[2]))
        #         swu_toa_clr_vis_[i] =(float(lines[ind20km[i]+(i*nlines)+1].split()[5]))
        #         swu_toa_clr_nir_[i] =(float(lines[ind20km[i]+(i*nlines)+1].split()[6]))
        #         swd_sfc_all_vis_[i] =(float(lines[41+(i*nlines)].split()[3]))
        #         swd_sfc_all_nir_[i] =(float(lines[41+(i*nlines)].split()[4]))
        #         swd_sfc_clr_vis_[i] =(float(lines[41+(i*nlines)].split()[7]))
        #         swd_sfc_clr_nir_[i] =(float(lines[41+(i*nlines)].split()[8]))
        #         swu_sfc_all_vis_[i] =(float(lines[41+(i*nlines)].split()[1]))
        #         swu_sfc_all_nir_[i] =(float(lines[41+(i*nlines)].split()[2]))
        #         swu_sfc_clr_vis_[i] =(float(lines[41+(i*nlines)].split()[5]))
        #         swu_sfc_clr_nir_[i] =(float(lines[41+(i*nlines)].split()[6]))
        #         direct_sfc_all_[i] =(float(lines[(nlines-5)+(i*nlines)].split()[0]))
        #         diffuse_sfc_all_[i] =(float(lines[(nlines-5)+(i*nlines)].split()[1]))
        #         direct_sfc_clr_[i] =(float(lines[(nlines-2)+(i*nlines)].split()[0]))
        #         diffuse_sfc_clr_[i] =(float(lines[(nlines-2)+(i*nlines)].split()[1]))
                
                
        # swd_toa_vis[indi] = swd_toa_vis_
        # swd_toa_nir[indi] = swd_toa_nir_       
        # swu_toa_all_nir[indi] = (swu_toa_all_nir_ * cft) + (swu_toa_clr_nir_ * (1-cft))
        # swu_toa_all_vis[indi] = (swu_toa_all_vis_ * cft) + (swu_toa_clr_vis_ * (1-cft))
        # swd_sfc_all_nir[indi] = (swd_sfc_all_nir_ * cft) + (swd_sfc_clr_nir_ * (1-cft))
        # swd_sfc_all_vis[indi] = (swd_sfc_all_vis_ * cft) + (swd_sfc_clr_vis_ * (1-cft))
        # swu_sfc_all_nir[indi] = (swu_sfc_all_nir_ * cft) + (swu_sfc_clr_nir_ * (1-cft))
        # swu_sfc_all_vis[indi] = (swu_sfc_all_vis_ * cft) + (swu_sfc_clr_vis_ * (1-cft))
        # direct_sfc_all[indi] = (direct_sfc_all_ * cft) + (direct_sfc_clr_ * (1-cft))     
        # diffuse_sfc_all[indi] = (diffuse_sfc_all_ * cft) + (diffuse_sfc_clr_ * (1-cft))
        # swu_toa_clr_nir[indi] = swu_toa_clr_nir_
        # swu_toa_clr_vis[indi] = swu_toa_clr_vis_
        # swd_sfc_clr_nir[indi] = swd_sfc_clr_nir_
        # swd_sfc_clr_vis[indi] = swd_sfc_clr_vis_
        # swu_sfc_clr_nir[indi] = swu_sfc_clr_nir_
        # swu_sfc_clr_vis[indi] = swu_sfc_clr_vis_
        # direct_sfc_clr[indi] = direct_sfc_clr_
        # diffuse_sfc_clr[indi] = diffuse_sfc_all_
          
        
    #SAVE AS NETCDF
    output_dict ={
        'SWD_TOA_VIS': (('time'), swd_toa_vis),
        'SWD_TOA_NIR': (('time'), swd_toa_nir),
        'SWU_TOA_VIS_ALL': (('time'), swu_toa_all_vis),
        'SWU_TOA_NIR_ALL': (('time'), swu_toa_all_nir),
        'SWU_TOA_VIS_CLR': (('time'), swu_toa_clr_vis),
        'SWU_TOA_NIR_CLR': (('time'), swu_toa_clr_nir),
        'SWD_SFC_VIS_ALL': (('time'), swd_sfc_all_vis),
        'SWD_SFC_NIR_ALL': (('time'), swd_sfc_all_nir),
        'SWD_SFC_VIS_CLR': (('time'), swd_sfc_clr_vis),
        'SWD_SFC_NIR_CLR': (('time'), swd_sfc_clr_nir),
        'SWU_SFC_VIS_ALL': (('time'), swu_sfc_all_vis),
        'SWU_SFC_NIR_ALL': (('time'), swu_sfc_all_nir),
        'SWU_SFC_VIS_CLR': (('time'), swu_sfc_clr_vis),
        'SWU_SFC_NIR_CLR': (('time'), swu_sfc_clr_nir),
        'DIRECT_SFC_ALL': (('time'), direct_sfc_all),
        'DIFFUSE_SFC_ALL': (('time'), diffuse_sfc_all),
        'DIRECT_SFC_CLR': (('time'), direct_sfc_clr),
        'DIFFUSE_SFC_CLR': (('time'), diffuse_sfc_clr)
        }
    
    coords_dict = {
        'time': time_ssf
        }
    
    ds_output = xr.Dataset(output_dict, coords=coords_dict)
    #ds_output.to_netcdf(f'/home/CCCma/L2/output/Allsky_SimpleRTM_era5_SSF_{yymmdd}.nc')
    
  #  return ds_output



# import time 
# start_time = time.time()

# year = 2018
# month = 3 
# day = 22

# ds_out = L2_calc_main(year,month, day)


# end_time = time.time()
# print('Run Time: ',end_time - start_time)
