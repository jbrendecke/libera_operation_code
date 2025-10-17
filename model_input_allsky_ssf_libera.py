#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sept 10 12:45:19 2024
#load imput for CCCma with the inclusion of clouds from CERES SSF running the 
RTM for each footprint
First iteration for Libera operation, no CCCma to more easily set things up with docker/AWS
@author: jbrendecke
"""

import numpy as np
import subprocess
import os
from netCDF4 import Dataset
import xarray as xr
from metpy.calc import specific_humidity_from_dewpoint, relative_humidity_from_dewpoint
from metpy.units import units
import glob
import sys
import logging
import warnings
warnings.filterwarnings("ignore")

# get the doy
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

# get profiles from reanalysis
def get_atmo_profiles(dataset_lev, dataset_sfc, lat_fp, lon_fp, time_fp, lyr):
    #///////////////////////////////////////////////////////////////////
    #extract out reanalysis profiles for nearest footprint/time of global reanalysis
    # Inputs: reanalysis Nc file dataset for levels and surface
    # Footprint latitude, longitdue, and time
    # number of atmospheric levels to be used in CCCma
    #//////////////////////////////////////////////////////
    #///////////////////////////////////////////////////////////////////
    dl = dataset_lev.sel(latitude = lat_fp, 
                         longitude = lon_fp, 
                         valid_time = time_fp,
                         method='nearest')
    ds = dataset_sfc.sel(latitude = lat_fp, 
                         longitude = lon_fp, 
                         valid_time = time_fp,
                         method='nearest')
    

    #convert geopotential to altitude
    earth_radius = 6378137 
    dl['z'] = dl['z']/9.8
    ds['z'] = ds['z']/9.8
    dl['z'] = (earth_radius * dl['z'] / (earth_radius - dl['z']))/1000
    ds['z'] = (earth_radius * ds['z'] / (earth_radius - ds['z']))/1000

    #determine specifc humidity from dewpoint/pressure
    q_sfc = specific_humidity_from_dewpoint(ds['sp']/100 * units.hPa, ds['d2m'] * units.kelvin)
    rh_sfc = relative_humidity_from_dewpoint(ds['t2m'] * units.kelvin, ds['d2m'] * units.kelvin)*100

    #make sure variables are not negative
    dl['r'] = xr.where(dl['r'] < 0, 0, dl['r'])
    dl['r'] = xr.where(dl['r'] > 100, 100, dl['r'])
    dl['q'] = xr.where(dl['q'] < 0, 0, dl['q'])
    dl['o3'] = xr.where(dl['o3'] < 0, 0, dl['o3'])
    dl['t'] = xr.where(dl['t'] < 0, 300, dl['t'])

    #set variables equal to NaN where pressure levels in levels dataset are higher than surface pressure(i.e. below ground)
    sp = ds['sp']/100
    mask = dl.pressure_level > sp
    dl = dl.where(~mask)

    #extract variables from era5    
    p_lev = dl['pressure_level'].data
    z_lev = dl['z'].data
    t_lev = dl['t'].data
    q_lev = dl['q'].data
    rh_lev = dl['r'].data
    o3_lev = dl['o3'].data
    p_sfc = ds['sp'].data/100
    z_sfc = ds['z'].data
    t_sfc = ds['t2m'].data
    
    #remove NaN values
    indz = np.where((~np.isnan(z_lev)) & (z_lev >= 0))
    p_lev = p_lev[indz]
    z_lev = z_lev[indz]
    t_lev = t_lev[indz]
    q_lev = q_lev[indz]
    rh_lev = rh_lev[indz]
    o3_lev = o3_lev[indz]
    
    #check height/pressure increases/decreases with height
    if z_sfc < 0.0: #CCCma doesn't allow negative height/below sea-level 
        z_sfc = 0.
        
    if z_sfc >= z_lev[0]: #sometimes surface altitude is higher than lowest level height
        z_sfc = z_lev[0] - ((z_sfc -z_lev[0])+.03) # make surface 3 meters lower\
            
    if p_sfc <= p_lev[0]: #ensure sfc pressure is higher than lowest level pressure
        p_sfc = p_lev[0] + ((p_lev[0] - p_sfc)+2) # add 2 hps
      #  print('p:', itime,ilat,ilon)
        
    #if (z_s >= z_l[0]): 
    #    print('z lower:', itime,ilat,ilon, z_s - z_l[0])
        
    #if (z_s <0.):
    #    print('z less than 0', itime,ilat,ilon, z_s)
    
    #combine surface level and pressure level
    press = np.concatenate(([p_sfc], p_lev))
    hts = np.concatenate(([z_sfc], z_lev))
    temps = np.concatenate(([t_sfc], t_lev))
    qvs = np.concatenate(([q_sfc], q_lev))
    rhs = np.concatenate(([rh_sfc], rh_lev))
    o3s = np.concatenate(([o3_lev[0]], o3_lev))
    
    #interpolate profile to 41 layers
    x = np.linspace(0, 1, len(press))  # Example x-coordinates, adjust if needed.
    new_x = np.linspace(0, 1, lyr)  # New x-coordinates for interpolation
    
    presx = np.interp(new_x, x, press)  
    htx = np.interp(new_x, x, hts)
    rhx = np.interp(new_x, x, rhs)
    tempx = np.interp(new_x, x, temps)
    qvx = np.interp(new_x, x, qvs)
    o3x = np.interp(new_x, x, o3s)
   
    #these variables use layer average
    tempa = np.zeros(lyr-1)
    qva = np.zeros(lyr-1)
    o3a = np.zeros(lyr-1)
    for zz in range(1,nlev):
        tempa[zz-1] = (tempx[zz-1]+tempx[zz])/2
        qva[zz-1] = (qvx[zz-1]+qvx[zz])/2
        o3x[zz-1] = (o3x[zz-1]+o3x[zz])/2
    
    return presx, htx, tempa, qva, rhx, o3a
#//////////////////////////////////////////////////////////////////////////////

# determine cloud profiles of Rc/Ri and LWC/IWC
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

# determine cloud profiles of CF
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

# get cloud profiles 
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

# very simple RTM to replace CCCma for now & keep everything in python
def simple_RTM(sza, cf, COD, AOD):
    nlen = len(sza)
    S0 = 1361 * np.cos(sza * (np.pi/180))
    w0 = .9
    asym =.9
    
    swd_dir = np.zeros(nlen); swd_dif = np.zeros(nlen)
    swd_dir_clr = np.zeros(nlen); swd_dif_clr = np.zeros(nlen)
    toa = np.zeros(nlen); toa_clr = np.zeros(nlen)
    for i in range(nlen):
        if cf[i] > .5:
            tau = max(COD[i]) + 5
        else:
            tau = AOD[i] + 5
            
        swd_dir[i] = S0[i] * np.exp(-1 * tau * (1/np.cos(sza[i]*np.pi/180))) 
        swd_dif[i] = S0[i] * w0 * (1 - np.exp(-1 * tau * (1/np.cos(sza[i]*np.pi/180)))) * asym
        
        swd = swd_dir[i] +swd_dif[i]
        
        if cf[i] > .5:
            albedo = .85
            toa = S0 * albedo
        else:
            albedo = .5
            toa = (swd * albedo) * np.exp(-1 * tau) 
        
        if cf[i] > .5:
            tau = AOD[i] + 5
            swd_dir_clr[i] = S0[i] * np.exp(-1 * tau * (1/np.cos(sza[i]*np.pi/180))) 
            swd_dif_clr[i] = S0[i] * w0 * (1 - np.exp(-1 * tau * (1/np.cos(sza[i]*np.pi/180)))) * asym
            
            swd = swd_dir[i] +swd_dif[i]
            albedo = .5
            toa_clr[i] = (swd * albedo) * np.exp(-1 * tau)
        else:
            swd_dir_clr[i] = swd_dir[i]
            swd_dif_clr[i] = swd_dif[i]
            toa_clr[i] = toa[i]


    return swd_dir, swd_dif, toa, swd_dir_clr, swd_dif_clr, toa_clr    
#//////////////////////////////////////////////////////////////////////////////
# Main function 
#//////////////////////////////////////////////////////////////////////////////

yymmdd='20180322'
jday_int = set_jday(int(yymmdd[0:4]), int(yymmdd[4:6]), int(yymmdd[6:]))

    
if (int(yymmdd[4:6]) >=3) & (int(yymmdd[4:6]) <= 8):
    season = 1
else:
    season = 0
    
#//////////////////LOAD SSF Non-cloud properties///////////////////////////////
path_ceres ='/home/jbrendecke/rad_global/SSF_cloud/data/Terra/'
filessf = glob.glob(path_ceres+'*'+yymmdd+'*.nc')
dc = xr.open_dataset(filessf[0])
time_ssf=np.array(dc.time).copy()
lon_ssf = np.array(dc.lon)
lat_ssf = np.array(dc.lat)
aod_ssf = (dc['PSF_wtd_MOD04_deep_blue_aerosol_optical_depth_land__0_550_'].data)
cod_ssf = (dc['Mean_visible_optical_depth_for_cloud_layer'].data)
sza_ssf = (dc['CERES_solar_zenith_at_surface'].data)
aod_land = (dc['PSF_wtd_MOD04_deep_blue_aerosol_optical_depth_land__0_550_'].data)
aod_ocean = (dc['PSF_wtd_MOD04_effective_optical_depth_average_ocean__0_550_'].data)
igbp = (dc['Surface_type_index'].data)

#time is number of days from 1970 and because of that, datatime is wrong on the data, this fixes it
correct_date = np.datetime64('2018-03-22T00:00:00')
time_ssf = correct_date + (time_ssf - time_ssf[0])

#right now assume most dominant scene type is used
igbp = igbp[:,0]
#adjust to CCCma ID
igbp = np.int_(igbp + 40)

ntime=len(time_ssf)
nlev = 41
nlay = nlev-1


#////////////////////////ERA5 Atmospheric Profiles////////////////////////////
file_lev = glob.glob('/home/jbrendecke/rad_global/ERA5/*level*'+yymmdd+'*.nc')[0]
file_sfc = glob.glob('/home/jbrendecke/rad_global/ERA5/*surface*'+yymmdd+'*.nc')[0]

if not os.path.isfile(file_lev) | os.path.isfile(file_sfc):
    logging.error('ERA5 File Not Found')
    
#open dataset and select 1 deg positions
dl = xr.open_dataset(file_lev)
ds = xr.open_dataset(file_sfc)


#///////////////////////Output Variables to be saved///////////////////////////
swd_toa_vis = np.zeros(ntime); swu_toa_all_vis = np.zeros(ntime); swu_toa_clr_vis = np.zeros(ntime);
swd_toa_nir = np.zeros(ntime); swu_toa_all_nir = np.zeros(ntime); swu_toa_clr_nir = np.zeros(ntime);
swu_20km_all_vis = np.zeros(ntime); swu_20km_clr_vis = np.zeros(ntime);
swu_20km_all_nir = np.zeros(ntime); swu_20km_clr_nir = np.zeros(ntime);
swd_sfc_all_vis = np.zeros(ntime); swd_sfc_clr_vis = np.zeros(ntime);
swu_sfc_all_vis = np.zeros(ntime); swu_sfc_clr_vis = np.zeros(ntime);
swd_sfc_all_nir = np.zeros(ntime); swd_sfc_clr_nir = np.zeros(ntime);
swu_sfc_all_nir = np.zeros(ntime); swu_sfc_clr_nir = np.zeros(ntime);
direct_sfc_all = np.zeros(ntime); direct_sfc_clr = np.zeros(ntime);
diffuse_sfc_all = np.zeros(ntime); diffuse_sfc_clr = np.zeros(ntime);

#%%
#//////////////////RUN CKD CODE FOR EACH SSF////////////
#select out daylight SSFs
nfoots = len(time_ssf)
ind_day = np.where(sza_ssf <= 89.5)[0]
max0 = 100 #30000 #max number of ilg loops that is ran in fortran
nloop = np.floor(len(ind_day) / max0)
remainder = ind_day % max0

# # # # # Run Loop for each CCCma ilg loop # # # # #
for ilg in range(5):#range(nloop+1):
    
    #determine footprints with daylight
    if ilg == nloop:
        #for last remaining footprints
        indi = ind_day[len(ind_day)-remainder:]
    else:
        # most footprints; length equals max0
        indi = ind_day[max0*ilg:max0*(ilg+1)]
    nfoots = len(indi)

    #save properties for RTM Input
    aod = aod_ssf[indi]
    cod = cod_ssf[indi]
    aodtype = np.zeros(nfoots, dtype=int)
    sza = sza_ssf[indi]
    igbptyp = igbp[indi]
    
    tigbp = np.zeros(nfoots)
    tmpaerotyp = np.zeros(nfoots, dtype=int) 
    ind20km = np.zeros(nfoots)
    presx = np.zeros((nfoots, nlev))
    htx = np.zeros((nfoots, nlev))
    rhx = np.zeros((nfoots, nlev))
    tempx = np.zeros((nfoots, nlay))
    qvx = np.zeros((nfoots, nlay))
    o3x = np.zeros((nfoots, nlay))
    cf_totalx = np.zeros(nfoots)
    relx = np.zeros((nfoots, nlay))
    reix = np.zeros((nfoots, nlay))
    lwcx = np.zeros((nfoots, nlay))
    iwcx = np.zeros((nfoots, nlay))
    cf_px = np.zeros((nfoots,nlay))
    cft = np.zeros(nfoots)
    kk=0
    
    
    dl2 = dl.sel(latitude = lat_ssf[indi[0]], 
                         longitude = lon_ssf[indi[0]], 
                         valid_time = time_ssf[indi[0]],
                         method='nearest')
    
    t_lev1 = dl['t'].data
    t_lev2 = dl2['t'].data

    
    sys.exit()
    # # # # # RUN LOOP FOR EACH FOOTPRINT # # # # #
    for ii in indi:
        # Determine Aerosol Properties
        if np.isnan(aod[kk]):
            aod[kk] = 0.10
        if igbptyp[kk] == 57:
            aodtype[kk] = 1
        elif igbptyp[kk] == 56:
            aodtype[kk] = 4
        elif abs(lat_ssf[indi[kk]]) > 70:
            aodtype[kk] = 5
        else:
            aodtype[kk] = 2
        
        #/////////////Save ERA5 Profile/////////////
        #Pressure, Hieght, RH should have 41 levels
        #Temperature, QV, O3 should have 40 levels
        pres, ht, temp, qv, rh, o3 = get_atmo_profiles(dl, ds, 
                                                       lat_ssf[ii], lon_ssf[ii],
                                                       time_ssf[ii], nlev)
        presx[kk, :] = np.flip(pres)
        htx[kk, :] = np.flip(ht)
        tempx[kk, :] = np.flip(temp)
        qvx[kk, :] = np.flip(qv)
        rhx[kk, :] = np.flip(rh)
        o3x[kk, :] = np.flip(o3)
        
        #identify index closest to 20km altitude
        ind20km[kk] = np.argmin(abs(20 - ht))
        
        
        #/////////////Save Cloud Profile///////////
        #40 levels of Re liq, Re ice, LWC & IWC, and CF
        rc, lwc, ri, iwc, cf_p, cf_tot = get_cloud_profile(dc, ht, pres, 
                                                           lat_ssf[ii], lon_ssf[ii],
                                                           time_ssf[ii], nlev)
        relx[kk, :] = np.flip(rc)
        reix[kk, :] = np.flip(ri)
        lwcx[kk, :] = np.flip(lwc)
        iwcx[kk, :] = np.flip(iwc)
        cf_px[kk, :] = np.flip(cf_p)
        cft[kk] = cf_tot
        
        #/////////////
        kk += 1

    # # # # # Run RTM # # # # #
    swd_dir, swd_dif, toa, swd_dir_clr, swd_dif_clr, toa_clr = simple_RTM(sza, cft, cod, aod)
    
    swd = swd_dir + swd_dif
    swd_clr = swd_dir_clr + swd_dif_clr    
    
    swd_toa_nir[indi] = 1365 *np.cos(sza *(np.pi/180)) * .55
    swd_toa_vis[indi] = 1365 *np.cos(sza *(np.pi/180)) * .45
    swu_toa_all_nir[indi] = toa *.55
    swu_toa_all_vis[indi] = toa * .45
    swu_toa_clr_nir[indi] = toa_clr *.55
    swu_toa_clr_vis[indi] = toa_clr * .45
    swu_20km_all_nir[indi] = toa *.55
    swu_20km_all_vis[indi] = toa * .45
    swu_20km_clr_nir[indi] = toa_clr *.55
    swu_20km_clr_vis[indi] = toa_clr * .45
    swd_sfc_all_nir[indi] = (swd_dir + swd_dif) * .55
    swd_sfc_all_vis[indi] = (swd_dir + swd_dif) * .45
    swd_sfc_clr_nir[indi] = (swd_dir_clr + swd_dif_clr) * .55
    swd_sfc_clr_vis[indi] = (swd_dir_clr + swd_dif_clr) * .45
    swu_sfc_all_nir[indi] = swd_sfc_all_nir[indi] * .5
    swu_sfc_all_vis[indi] = swd_sfc_all_nir[indi] * .5
    swu_sfc_clr_nir[indi] = swd_sfc_all_nir[indi] * .5
    swu_sfc_clr_vis[indi] = swd_sfc_all_nir[indi] * .5
    direct_sfc_all[indi] = swd_dir
    diffuse_sfc_all[indi] = swd_dif
    direct_sfc_clr[indi] = swd_dir_clr
    diffuse_sfc_clr[indi] = swd_dif_clr
    
    
    # # RUN CCCma
    # ctx=nfoots
    # os.chdir('/home/jbrendecke/CKD/ver10_global/')
    # pathi ='/home/jbrendecke/CKD/ver10_global/input/'
    # patho ='/home/jbrendecke/CKD/ver10_global/output/'
    
    # #WRITE INPUT FILE AND ADJUST UA_MAIN.F90
    # with open('UA_main.F90', 'r') as file:
    #     dmain = file.readlines()

    # dmain[4] = f'   parameter (lay = {nlay}, lev = lay +1, ilg={ctx}, il1 = 1, il2 = {ctx}, modlay=34, nxloc=100)\n'

    # with open('UA_main.F90', 'w') as file:
    #     file.writelines(dmain)
    
    # with open('input.dat', 'w') as file:
    #     for k in range(nfoots):
    #         file.write(f"{sza[k]:11.4f} {aod[k]:11.4f} {igbptyp[k]:4d} {season:4d} {aodtype[k]:4d} {jday_int:5d} {cft[k]:7.4f}\n")
    #         for ilyr in range(nlay):
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
    # swd_toa_vis_ = np.zeros(max0)
    # swd_toa_nir_ = np.zeros(max0)
    # swu_toa_all_vis_ = np.zeros(max0)
    # swu_toa_all_nir_ = np.zeros(max0)
    # swu_toa_clr_vis_ = np.zeros(max0)
    # swu_toa_clr_nir_ = np.zeros(max0)
    # swu_20km_all_vis_ = np.zeros(max0)
    # swu_20km_all_nir_ = np.zeros(max0)
    # swu_20km_clr_vis_ = np.zeros(max0)
    # swu_20km_clr_nir_ = np.zeros(max0)
    # swd_sfc_all_vis_ = np.zeros(max0)
    # swd_sfc_all_nir_ = np.zeros(max0)
    # swd_sfc_clr_vis_ = np.zeros(max0)
    # swd_sfc_clr_nir_ = np.zeros(max0)
    # swu_sfc_all_vis_ = np.zeros(max0)
    # swu_sfc_all_nir_ = np.zeros(max0)
    # swu_sfc_clr_vis_ = np.zeros(max0)
    # swu_sfc_clr_nir_ = np.zeros(max0)
    # direct_sfc_all_ = np.zeros(max0)
    # diffuse_sfc_all_ = np.zeros(max0)
    # direct_sfc_clr_ = np.zeros(max0)
    # diffuse_sfc_clr_ = np.zeros(max0)
    
    # with open(output_file[0], 'r') as ofile:
    #     lines = ofile.readlines()
    #     for i in range(nfoots):
    #         swd_toa_vis_w[i] =(float(lines[1+(i*nlines)].split()[3]))
    #         swd_toa_nir_w[i] =(float(lines[1+(i*nlines)].split()[4]))
    #         swu_toa_all_vis_w[i] =(float(lines[1+(i*nlines)].split()[1]))
    #         swu_toa_all_nir_w[i] =(float(lines[1+(i*nlines)].split()[2]))
    #         swu_toa_clr_vis_w[i] =(float(lines[1+(i*nlines)].split()[5]))
    #         swu_toa_clr_nir_w[i] =(float(lines[1+(i*nlines)].split()[6]))
    #         swu_20km_all_vis_w[i] =(float(lines[ind20km[i]+(i*nlines)+1].split()[1]))
    #         swu_20km_all_nir_w[i] =(float(lines[ind20km[i]+(i*nlines)+1].split()[2]))
    #         swu_20km_clr_vis_w[i] =(float(lines[ind20km[i]+(i*nlines)+1].split()[5]))
    #         swu_20km_clr_nir_w[i] =(float(lines[ind20km[i]+(i*nlines)+1].split()[6]))
    #         swd_sfc_all_vis_w[i] =(float(lines[41+(i*nlines)].split()[3]))
    #         swd_sfc_all_nir_w[i] =(float(lines[41+(i*nlines)].split()[4]))
    #         swd_sfc_clr_vis_w[i] =(float(lines[41+(i*nlines)].split()[7]))
    #         swd_sfc_clr_nir_w[i] =(float(lines[41+(i*nlines)].split()[8]))
    #         swu_sfc_all_vis_w[i] =(float(lines[41+(i*nlines)].split()[1]))
    #         swu_sfc_all_nir_w[i] =(float(lines[41+(i*nlines)].split()[2]))
    #         swu_sfc_clr_vis_w[i] =(float(lines[41+(i*nlines)].split()[5]))
    #         swu_sfc_clr_nir_w[i] =(float(lines[41+(i*nlines)].split()[6]))
    #         direct_sfc_all_w[i] =(float(lines[(nlines-5)+(i*nlines)].split()[0]))
    #         diffuse_sfc_all_w[i] =(float(lines[(nlines-5)+(i*nlines)].split()[1]))
    #         direct_sfc_clr_w[i] =(float(lines[(nlines-2)+(i*nlines)].split()[0]))
    #         diffuse_sfc_clr_w[i] =(float(lines[(nlines-2)+(i*nlines)].split()[1]))
            
   
    #SAVE AS NETCDF
    output_dict ={
        'SWD_TOA_VIS': (('time'), swd_toa_vis),
        'SWD_TOA_NIR': (('time'), swd_toa_nir),
        'SWU_TOA_VIS_ALL': (('time'), swu_toa_all_vis),
        'SWU_TOA_NIR_ALL': (('time'), swu_toa_all_nir),
        'SWU_TOA_VIS_CLR': (('time'), swu_toa_clr_vis),
        'SWU_TOA_NIR_CLR': (('time'), swu_toa_clr_nir),
        'SWU_20km_VIS_ALL': (('time'), swu_20km_all_vis),
        'SWU_20km_NIR_ALL': (('time'), swu_20km_all_nir),
        'SWU_20km_VIS_CLR': (('time'), swu_20km_clr_vis),
        'SWU_20km_NIR_CLR': (('time'), swu_20km_clr_nir),
        'SWD_SFC_VIS_ALL': (('time'), swd_sfc_all_vis),
        'SWD_SFC_NIR_ALL': (('time'), swd_sfc_all_nir),
        'SWD_SFC_VIS_CLR': (('time'), swd_sfc_clr_vis),
        'SWD_SFC_NIR_CLR': (('time'), swd_sfc_clr_nir),
        'SWU_SFC_VIS_ALL': (('time'), swu_sfc_all_vis),
        'SWU_SFC_NIR_ALL': (('time'), swu_sfc_all_nir),
        'SWU_SFC_VIS_CLR': (('time'), swu_sfc_clr_vis),
        'SWU_SFC_NIR_CLR': (('time'), swu_sfc_clr_nir),
        'DIRECT_SFC_ALL': (('time'), direct_sfc_clr),
        'DIFFUSE_SFC_ALL': (('time'), diffuse_sfc_clr),
        'DIRECT_SFC_CLR': (('time'), direct_sfc_clr),
        'DIFFUSE_SFC_CLR': (('time'), diffuse_sfc_clr)
        }
    
    coords_dict = {
        'time': time_ssf
        }
    
    #ds_output = xr.Dataset(output_dict, coords=coords_dict)
    #ds_output.to_netcdf(f'/home/jbrendecke/rad_global/SSF_cloud/SimpleRTM_output/Allsky_SimpleRTM_era5_SSF_{yymmdd}____.nc')
    
    
