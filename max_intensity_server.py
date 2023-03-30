#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 17:02:02 2020

@author: margaret

a bit of code for identifying the point of max intensity along a track,
saving the point to a netCDF file, and then making a density plot of those points
"""

# # # # # # # # # #
# import the things
# # # # # # # # # #
import glob
import sys
import matplotlib.pyplot as plt
#import pandas as pd
import numpy as np
import cartopy.crs as ccrs
import seaborn as sns

# netCDF stuff
import netCDF4
from netCDF4 import Dataset

level = int(sys.argv[1])
plot   = 0 #needs to be manually changed

#these capture all of the files for when we do a plot after all the loop stuff
if plot == 1:
    binlon = float(sys.argv[2])
    binlat = float(sys.argv[3])
    max_intensities_all = []
    max_lats_all = []
    max_lons_all = []


# # # # # # #
#do the thing
# # # # # # #
trackfiles = sorted(glob.glob(f'/data3/TRACK/TRACKOutput/MERRA2/FilteredNoTCUpdated/{level}hPa/post_filter.*.NoTC.nc'))

for trackfile in trackfiles:
    tracks = Dataset(trackfile, mode='r+')

    TrackID = tracks['TRACK_ID'][:]
    FirstPt = tracks['FIRST_PT'][:]
    NumPts  = tracks['NUM_PTS'][:]

    #the tracks themselves
    lats = tracks['latitude'][:]
    lons = tracks['longitude'][:]

    intensities = tracks['curvature_vorticity'][:]
    
    max_intensities = []
    max_lats = []
    max_lons = []

    # # # # # # # # # # # # # # # #
    #find location of max intensity
    # # # # # # # # # # # # # # # #

#probs going to need to make this more complex once I get the whole zero length track handled

    for i in range(len(TrackID)):
        firstpt = FirstPt[i]
        tracklons = lons[firstpt:(firstpt+NumPts[i])]
        tracklats = lats[firstpt:(firstpt+NumPts[i])]
        trackintensities = intensities[firstpt:(firstpt+NumPts[i])]
        max_inten = 0
        max_lat = 0
        max_lon = 0
        for k in range(len(trackintensities)-1):
            if trackintensities[k] > max_inten:
                max_inten = trackintensities[k]
                max_lat = tracklats[k]
                max_lon = tracklons[k]
            else:
                continue
        max_intensities.append(max_inten)
        max_lats.append(max_lat)
        max_lons.append(max_lon)
        
    # # # # # # # #
    # save the data
    # # # # # # # #
    
    # going to be fancy here and save the max intensities as another variable
    # rather than new files, so I don't have to pick where to put the files
    # and because this could be useful for seasonal breakdowns and other stuff
    
    #create a new variable in which to store these max intensities
    max_intensity = tracks.createVariable('max_int', 'f8', ('tracks'))
    max_int_lat   = tracks.createVariable('max_lat', 'f8', ('tracks'))
    max_int_lon   = tracks.createVariable('max_lon', 'f8', ('tracks'))
    
    #give the new variables some attributes
    max_intensity.standard_name = 'max intensity'
    max_intensity.long_name = 'Max Intensity'
    max_intensity.units = 's**-1'
    max_intensity.scale_factor = 1.e-05
    
    max_int_lat.standard_name = 'latitude'
    max_int_lat.long_name = 'Latitude'
    max_int_lat.units = 'degrees_north'
    
    max_int_lon.standard_name = 'longitude'
    max_int_lon.long_name = 'Longitude'
    max_int_lon.units = 'degrees_east'
    
    #now actually write the things
    max_intensity[:] = max_intensities
    max_int_lat[:] = max_lats
    max_int_lon[:] = max_lons
    
    #I think the following block was to fix something that didn't run right
    #where I just needed to rewrite the values and not remake the variable
# =============================================================================
#     tracks['max_int'][:] = max_intensities
#     tracks['max_lat'][:] = max_lats
#     tracks['max_lon'][:] = max_lons
#idk why it ran all of these but I think it did so we're going to hope for the best
# =============================================================================
    
    
    
    #also save to variable for plotting all the things
    if plot == 1:
        max_intensities_all.extend(max_intensities)
        max_lats_all.extend(max_lats)
        max_lons_all.extend(max_lons)
    
    #and close the file before we restart the loop
    tracks.close()
        
#fancy idea: put 850 shaded and 700 contoured
#and if I wanted to get fancy with the intensity boxplots, I'd do 850 and 700 next to each other
#I suppose that would involve saving the data somehow
#or just doing it all in here


# =============================================================================
# # # # # # # # # #
# # make some plots
# # # # # # # # # #
# if plot == 1:
#     fig = plt.figure(dpi=300)
#     ax = plt.axes(projection=ccrs.PlateCarree())
#     ax.coastlines(resolution='110m')
#     
#     #sns.kdeplot(max_lons, max_lats, shade=True, cbar=True, ax=ax,
#     #            levels=np.arange(start=0, stop=2.5e-4, step=2e-5), shade_lowest=False)
#     ax.set_title(f'Maximum Intensities {level}hPa')
#     #plt.show()
#     sns.histplot(x=max_lons, y=max_lats, stat='count', ax=ax, binwidth=(binlon, binlat), cbar=True)
#     plt.savefig(f'max_intensity_map.{level}.v2.png', bbox_inches = "tight", dpi = 300)
#     
#     fig2 = plt.figure(figsize=(8,5),dpi=300)
#     ax2 = plt.axes()
#     sns.boxplot(max_intensities)
#     ax2.set_xlim(0,0.0003)
#     ax2.set_title(f'TEW Intensities {level}hPa ')
#     ax2.set_xlabel('CV Intensity (s^-1)')
#     #plt.show()
#     plt.savefig(f'max_intensity_boxplot.{level}.png', bbox_inches='tight', dpi=300)
# =============================================================================
