#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Weds Dec 30 2020 0226 UTC

@author: margaret

remove_TCs_server.py
Code for removing TC (and recurving post-TC) portions of tracks
Uses IBTrACS.since1980.v04r00.nc for the TC tracks to match against
Matches TCs by date/time and distance from track
To match, wave and track must match time and be within 1000km

Features on server version:
    year takes argument from user
    add a level argument to be taken from user
    hemisphere stuff
"""

# packages to import
import sys
import numpy as np
import numpy.ma as ma
import math
from math import radians, degrees, sin, cos, asin, acos, sqrt

# netCDF stuff
import netCDF4
from netCDF4 import Dataset, num2date

# # # # # # # # # # # # # #
# User-supplied information
# # # # # # # # # # # # # #

year = int(sys.argv[1])
level = int(sys.argv[2])
hemi = int(sys.argv[3])
currentYear = 2022

# # # # # # # # # # # # # # # # # # # #
# opening files and accessing variables
# # # # # # # # # # # # # # # # # # # #

if hemi == 1:
    hemitext = 'northern'
    hemishort = 'NH'
elif hemi == 2:
    hemitext = 'southern'
    hemishort = 'SH'
else:
    print('Hemisphere is not valid.')
    sys.exit()

# First the TEW tracks
trackfile = f'/data3/TRACK/TRACKOutput/MERRA2/FilteredTracksUpdated/{level}hPa/post_filter.{year}.{hemishort}.nc'
alltracks = Dataset(trackfile, mode='r')

# the variables we need, and time converted to python datetime
TrackID = alltracks['TRACK_ID'][:]
NumPts  = alltracks['NUM_PTS'][:]
TrackTimeUnits = alltracks.variables['time'].units
numrecs = alltracks.dimensions['record'].size

print(f'There are {len(TrackID)} tracks in this file')

# Then IBTrACS
IBTrACSfile = '/data3/Data_Processed/IBTRACS/IBTrACS.since1980.v04r00.nc'
IBTrACS     = Dataset(IBTrACSfile, mode='r')

# let's not put all of IBTrACS into memory and only read in the year we need
IBTrACStime = IBTrACS['time'][:]
IBTrACStimeUnits = IBTrACS.variables['time'].units
IBTrACSDateTime  = num2date(IBTrACStime, units=IBTrACStimeUnits)

# let's not have to compare to all of the years that aren't in the year of interest
#find index of level of interest
first = 0
last = 0
for storm in range(len(IBTrACSDateTime[:,0])):
    if first == 0:
        if IBTrACSDateTime[storm,0].year == year:
            first = storm
    elif first != 0:
        if IBTrACSDateTime[storm,0].year == currentYear:
            last = storm
        elif IBTrACSDateTime[storm,0].year == (year+2):
            last = storm
IBTrACSDateTime = IBTrACSDateTime[first:last, :]
IBTrACSlat  = IBTrACS['lat'][first:last, :]
IBTrACSlon  = IBTrACS['lon'][first:last, :]

# # # # # # # # # # # #
# Function for Distance
# # # # # # # # # # # #
# found at https://medium.com/@petehouston/calculate-distance-of-two-locations-on-earth-using-python-1501b1944d97

def great_circle(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    return 6371 * (acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2)))


# # # # # # # # # # # # # # # # # # # # #
# set up a file for saving TC-less tracks
# # # # # # # # # # # # # # # # # # # # #

maskfile = f'/data3/TRACK/TRACKOutput/MERRA2/FilteredNoTCUpdated/{level}hPa/post_filter.{year}.{hemishort}.NoTC.nc'
masktracks = Dataset(maskfile, mode='w')

# Dimensions
out_tracks = masktracks.createDimension('tracks', alltracks.dimensions['tracks'].size)
out_record = masktracks.createDimension('record', None)

# Create variables
out_TrackID = masktracks.createVariable('TRACK_ID', 'i8', ('tracks'))
out_firstpt = masktracks.createVariable('FIRST_PT', 'i8', ('tracks'))
out_numpts  = masktracks.createVariable('NUM_PTS', 'i8', ('tracks'))
out_index   = masktracks.createVariable('index', 'i8', ('record'))
out_time    = masktracks.createVariable('time', 'f8', ('record'))
out_lons    = masktracks.createVariable('longitude', 'f8', ('record'))
out_lats    = masktracks.createVariable('latitude', 'f8', ('record'))
out_cv      = masktracks.createVariable('curvature_vorticity', 'f8', ('record'))

# Variable Attributes, Units, etc.
out_TrackID.add_fld_num = 0
out_TrackID.tot_add_fld_num = 0
out_TrackID.loc_flags = ""
out_TrackID.cf_role = "trajectory_id"

#out_firstpt does not have any attributes

out_numpts.long_name = 'number of obs for this trajectory'
out_numpts.sample_dimension = 'record'

#out_index does not have any attributes

out_time.standard_name = 'time'
out_time.long_name = 'Time'
out_time.units = TrackTimeUnits
out_time.time_calendar = 'gregorian'
out_time.start = alltracks.variables['time'].start
out_time.step = '3'

out_lons.standard_name = 'longitude'
out_lons.long_name = 'Longitude'
out_lons.units = 'degrees_east'

out_lats.standard_name = 'latitude'
out_lats.long_name = 'Latitude'
out_lats.units = 'degrees_north'

out_cv.standard_name = 'curvature_vorticity'
out_cv.long_name = 'Relative Vorticity'
out_cv.units = 's**-1'
out_cv.scale_factor = 1.e-05
out_cv.coordinates = 'time latitude longitude TRACK_ID'

# Global attributes
masktracks.realm = 'atmos'
masktracks.history = 'testing'
masktracks.Conventions = 'CF-1.7'
masktracks.featureType = 'trajectory'
masktracks.comment = f'Filtered tracks for {year}, with TCs removed. Filtered with version 5.1a, TCs removed with v2'

# # # # # # # # # # # # # # # # # # # # # # # # # # #
# Go through filtered tracks and try to match IBTrACS
# # # # # # # # # # # # # # # # # # # # # # # # # # #

# for each TrackID, cycle through the times. Attempt to match each time to one
# in IBTrACS. If times match, then check the lat/lon coordinates. If the two
# points are within 1000km of each other, call it a match. At least five
# points must match for the TC to count as matched tracks. Then, when saving,
# all points after the first match are considered TC, and will be removed.
# Repeat for all tracks.

#so the way to break out of this many loops is a function
def trackmatch():       
    for storm in range(len(IBTrACSDateTime)):
        for t in range(len(IBTrACSDateTime[storm,:])):
            if type(IBTrACSDateTime[storm,t]) == np.ma.core.MaskedConstant:
                break #break is if it reaches masked IBTrACS data, so that we proceed to the next storm
            elif track_y == (IBTrACSDateTime[storm,t].year):
                if track_m == (IBTrACSDateTime[storm,t].month):                        
                    if track_d == (IBTrACSDateTime[storm,t].day):
                        if track_h == (IBTrACSDateTime[storm,t].hour):
                            dist = great_circle(track_lon, track_lat, IBTrACSlon[storm,t], IBTrACSlat[storm,t])
                            if dist <= 1000:
                                #finally, we have a space *and* time match
                                return True
    #when we make it down here, it means that all points along the track have been compared to IBTrACS
    return False
                                    

out_TrackID[:] = alltracks['TRACK_ID'][:]
start = 0
for i in range(len(TrackID)):
    firstpt    = alltracks['FIRST_PT'][i]
    tracklons  = alltracks['longitude'][firstpt:(firstpt+NumPts[i])]
    tracklats  = alltracks['latitude'][firstpt:(firstpt+NumPts[i])]
    tracktimes = alltracks['time'][firstpt:(firstpt+NumPts[i])]
    trackdates = num2date(tracktimes, units= TrackTimeUnits)
    trackcv    = alltracks['curvature_vorticity'][firstpt:(firstpt+NumPts[i])]
    print(f'Now attempting to match the {i}th track of {len(TrackID)} total tracks')
    
    #reset firstmatch index and count of matches
    firstmatch = -1
    matches = 0
        
    for k in range(len(tracklons)-1):
        track_y   = trackdates[k].year
        track_m   = trackdates[k].month
        track_d   = trackdates[k].day
        track_h   = trackdates[k].hour
        track_lon = tracklons[k]
        track_lat = tracklats[k]
    
        match = trackmatch()
        if match == True and firstmatch == -1:
            firstmatch = k
        matches += match
        
    if matches >= 5:
        out_numpts[i] = firstmatch
        #make sure the *very* first point is right
        if i == 0:
            out_firstpt[i] = 0
        else:
            out_firstpt[i] = (out_firstpt[i-1] + out_numpts[i-1])
        
        #write the rest of the data
        end = start + out_numpts[i]
        
        all_start = alltracks['FIRST_PT'][i]
        all_end   = all_start + firstmatch
        
        out_index[start:end] = alltracks['index'][all_start:all_end]
        out_time[start:end]  = alltracks['time'][all_start:all_end]
        out_lons[start:end]  = alltracks['longitude'][all_start:all_end]
        out_lats[start:end]  = alltracks['latitude'][all_start:all_end]
        out_cv[start:end]    = alltracks['curvature_vorticity'][all_start:all_end]
    
        start = end
    else:
        out_numpts[i] = alltracks['NUM_PTS'][i]
        #make sure the *very* first point is right
        if i == 0:
            out_firstpt[i] = 0
        else:
            out_firstpt[i] = (out_firstpt[i-1] + out_numpts[i-1])
        
        #write the rest of the data
        end = start + out_numpts[i]
        
        all_start = alltracks['FIRST_PT'][i]
        all_end   = all_start + out_numpts[i]
        
        out_index[start:end] = alltracks['index'][all_start:all_end]
        out_time[start:end]  = alltracks['time'][all_start:all_end]
        out_lons[start:end]  = alltracks['longitude'][all_start:all_end]
        out_lats[start:end]  = alltracks['latitude'][all_start:all_end]
        out_cv[start:end]    = alltracks['curvature_vorticity'][all_start:all_end]
    
        start = end

#close the files
IBTrACS.close()
alltracks.close()
masktracks.close()