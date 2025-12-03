# Python and Fortran Script to run CCCma RTM

## Fortran Code:
- Reads inputs sza, IGBP type, AOD, day of year, aerosol type, total cloud fraction,
  ERA5 profiles (pressure, height, RH, temp, specific humidity, O3), and 
  cloud profile properties (Re liquid, Re ice, LWC, IWC, and CF)
- Hiearchy of code UA_main.F90 -> radiation.F90 -> raddriv11.F90 -> etc.
- Output profile fluxes for upward/downard all-sky and clear-sky conditions.
  Also outputs surface all/clear-sky direct and diffuse flux.

IGBP Interger Number:
41- Evergreen Needle Forest ! (= Mosart 14, pine forest)
42- Evergreen Broadleaf Forest ! (= Mosart 18, broadleaf-pine forest)
43- Deciduous Needle Forest ! (= Mosart 18, broadleaf-pine forest)
44- Deciduous Broadleaf Forest ! (= Mosart 6, broadleaf forest)
45- Mixed Forest ! (= Mosart 25, broadleaf 70-pine 30)
46- Closed Shrubs ! (= Mosart 22, pine-brush)
47- Open/Shrubs ! (= Mosart 40, broadleaf-brush)
48- Woody Savanna ! (= Mosart 20, soil-grass-scrub)
49- Savanna ! (= Mosart 19, grass-scrub)
50- Grassland ! (= Mosart 13, meadow grass)
51- Wetland ! (= Mosart 51, wetland)
52- Cropland ! (= Mosart 45, crop)
53- Urban ! (= Mosart 21, urban commercial)
54- Crop Mosaic ! (= Mosart 46, mixed-veg)
55- Antarctic Snow ! (= Mosart 9, old snow 1000 micron radius)
56- Barren/Desert ! (= Mosart 28, mixture of material (rock & silt-sand))
57- Ocean Water ! (= Mosart 1, water)
58- Tundra ! (= Mosart 16, tundra)
59- Fresh Snow ! (= Mosart 43, fresh snow (50 micron radius))
60- Sea Ice ! (= Mosart 10 sea ice, 3 meters thick)

Aerosol Integer Type:
0- No aerosols optical properties(Extinction/Absorption =0.0)
1- Oceanic Aerosols
2- Rural Aerosols
3- Urban Aerosols
4- Desert Aerosols
5- Background Tropospheric

## Python Code:
- input argurments:
        manifest: absolute path of input manifest file
        -v or --verbose: detail debug logging
- reads in .json input file to determine input file locations
- Loads CERES SSF properties from .nc file
- Loads ERA5 surface and levels .nc file
- Processes these files to get to column profiles
- loads profiles in txt file.
- MAXI is used to load different number profiles per txt file
  (CCCma code will read these as a loop)
- tells the computer to run fortran code with new input text file
- reads CCCma output text file and saves fluxes
- outputs both .nc file and .json file for output location



