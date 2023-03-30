#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sept 02 Z2020

@author: margaret

Version 6 of filtering TRACK output for *tropical* easterly waves and saving
a netCDF file of just TEWs.
New this version:
    Implemented restriction that easterly distance must be within tropics
    
5.1b:
    - dropping easterly/westerly and just using the trop/subtrop eastdist requirement
    - updated line 113 to abs(tracklats[k]) < 40 to try to eliminate SH clutter

5.1b.s:
    Server version changes:
        - changed year, level, hemisphere to take input from sys.argv
        - changed track plotting to be optional (this must be manually changed)
        - changed filepaths for those on Aardvark/Narwhal/Platypus
        
6:
    - added simple elif lines to remove features on the wrong side of the equator (lat < 0 for NH files
    and lat > 0 for SH files)   
"""

# # # # # # # # #
# import packages
# # # # # # # # #

# packages imported in my initial attempt to get things to plot
import numpy as np

# netCDF stuff
import netCDF4
from netCDF4 import Dataset

#packages for easy of running from terminal on the server
import sys

# # # # # # # # # # # 
# initial setup stuff
# # # # # # # # # # #

# variables from arguments when running code
year = int(sys.argv[1])
level = int(sys.argv[2]) #currently only 850 or 700hPa
hemi = int(sys.argv[3]) #use 1 for NH, 2 for SH

plot = 0 #change to 1 to create plots of TEW tracks with cfplot

# set human readable information about the sys.argv input for hemisphere, level
if hemi == 1:
    hemitext = 'northern'
    hemishort = 'NH'
elif hemi == 2:
    hemitext = 'southern'
    hemishort = 'SH'
else:
    print('Hemisphere is not valid.')
    sys.exit()

print(f'The provided year is {year}. The provided level is {level}hPa. The provided hemisphere is the {hemitext} hemisphere.')

# where is the track file
# will abort if a valid level or hemisphere is not given and it didn't earlier
if hemi == 1:
    trackfile = f'/data3/TRACK/TRACKOutput/MERRA2/{level}hPa/{year}/MERRA_CVOR{level}_T63FILT/ff_trs_pos.new.nc'
elif hemi == 2:
    trackfile = f'/data3/TRACK/TRACKOutput/MERRA2/{level}hPa/{year}/MERRA_CVOR{level}_T63FILT_SH/ff_trs_pos.new.nc'
else:
    print(f'Level provided was {level}hPa but no valid hemisphere was provided.')
    sys.exit()

alltracks  = Dataset(trackfile, mode='r')

TrackID = alltracks['TRACK_ID'][:]
NumPts = alltracks['NUM_PTS'][:]

# # # # # # # # # #
# filter track file
# # # # # # # # # #

# set up list in which to put the TrackIDs of easterly waves
eastID = []

# for each TrackID, set up a counter for whether the next lon is east or west
# of the current lon. Then iterate through the lons along the track, updating
# the counter. Once the track has been processed, if the track has mostly
# easterly motion, append the ID of the current track to the list of easterly
# tracks. Repeat for all tracks.
for i in TrackID:
    
    #eliminate things that originate poleward of 40
    #because even if it's easterly, if it originated that far north
    #it's not what we're looking for
    if abs(alltracks['latitude'][alltracks['FIRST_PT'][i]]) > 40:
        #print(f'Now skipping TRACKID {i} because it is poleward of 40')
        continue
    elif hemi == 1 and alltracks['latitude'][alltracks['FIRST_PT'][i]] < 0:
        continue
    elif hemi == 2 and alltracks['latitude'][alltracks['FIRST_PT'][i]] > 0:
        continue
    else:
        #now for trop/subtrop things, the basic filtering is just
        #are there more points where it's easterly or westerly
        #and record the TRACK ID index for the easterly things
        easterly = 0
        westerly = 0
        eastdist = 0
        firstpt = alltracks['FIRST_PT'][i]
        tracklons = alltracks['longitude'][firstpt:(firstpt+NumPts[i])]
        tracklats = alltracks['latitude'][firstpt:(firstpt+NumPts[i])]
        for k in range(len(tracklons)-1):
            #if (tracklons[k] < tracklons[k+1]): #I want better filtering for things that cross 0/360 # so this commented bit is leftover from a previous version of this
            #    westerly += 1
            #else:
            #    easterly +=1
            #    print('now for distance')
            if abs(tracklats[k]) < 40:
                if hemi == 1 and tracklats[k] < 0:
                    continue
                elif hemi == 2 and tracklats[k] > 0:
                    continue
                elif -15 < tracklons[k+1] - tracklons[k] < 0:
                    eastdist += abs(tracklons[k+1] - tracklons[k])
                    #print('normal easterly distance')
                elif tracklons[k+1] - tracklons[k] > 300:
                    eastdist += abs(tracklons[k+1] - (tracklons[k]+360))
                    #print('cross prime meridian easterly distance')
                else:           
                    #print('did not append easterly distance')
                    continue
        #if easterly >= westerly and eastdist >= 15:
        if eastdist >= 15: #may play with this distance
            eastID.append (i)
            #print(f'TRACKID {i} found to be tropical easterly, appended to eastID')
        else:
            #print(f'TRACKID {i} is not a TEW, continuing')
            continue

NumEast =  len(eastID)

# # # # # # # # # # # # # # #
# set up, save to output file
# # # # # # # # # # # # # # #

# Things that need to be the same from the old file to the new one
timeunits = alltracks.variables['time'].units
firsttime = alltracks.variables['time'].start

# Create file and open it
outfile = f'/data3/TRACK/TRACKOutput/MERRA2/FilteredTracksUpdated/{level}hPa/post_filter.{year}.{hemishort}.nc'
east_out = Dataset(outfile, 'w', format='NETCDF4')

# Dimensions
out_tracks = east_out.createDimension('tracks', NumEast)
out_record = east_out.createDimension('record', None)

# Create variables
out_TrackID = east_out.createVariable('TRACK_ID', 'i8', ('tracks'))
out_firstpt = east_out.createVariable('FIRST_PT', 'i8', ('tracks'))
out_numpts  = east_out.createVariable('NUM_PTS', 'i8', ('tracks'))
out_index   = east_out.createVariable('index', 'i8', ('record'))
out_time    = east_out.createVariable('time', 'f8', ('record'))
out_lons    = east_out.createVariable('longitude', 'f8', ('record'))
out_lats    = east_out.createVariable('latitude', 'f8', ('record'))
out_cv      = east_out.createVariable('curvature_vorticity', 'f8', ('record'))

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
out_time.units = timeunits
out_time.time_calendar = 'gregorian'
out_time.start = firsttime
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
east_out.realm = 'atmos'
east_out.history = 'testing'
east_out.Conventions = 'CF-1.7'
east_out.featureType = 'trajectory'
east_out.comment = f'Filtered tracks for {year}. Created with filtering version 6'

# Write the data
out_TrackID[:] = eastID
start = 0
for j in range(len(eastID)):
    out_numpts[j]  = alltracks['NUM_PTS'][eastID[j]]
    if j == 0:
        out_firstpt[j] = 0
    else:
        out_firstpt[j] = (out_firstpt[j-1] + out_numpts[j-1])
          
    end = start + out_numpts[j]
    
    all_start = alltracks['FIRST_PT'][eastID[j]]
    all_end   = all_start + out_numpts[j]
    
    out_index[start:end] = alltracks['index'][all_start:all_end]
    out_time[start:end]  = alltracks['time'][all_start:all_end]
    out_lons[start:end]  = alltracks['longitude'][all_start:all_end]
    out_lats[start:end]  = alltracks['latitude'][all_start:all_end]
    out_cv[start:end]    = alltracks['curvature_vorticity'][all_start:all_end]
    
    start = end

#close the files
east_out.close()
alltracks.close()
print('Tracks have been filtered and saved.')
print(f'The filtered tracks are located at {outfile}')

# # # # # # # # # # # # # #
# Do a plotting with cfplot
# # # # # # # # # # # # # #

if plot == 1:
    # packages for traj to work
    import cf
    import cfplot as cfp
    
    f = cf.read(outfile)[1]
    cfp.traj(f, title=f'Easterly Tracks {year}', marker=None)
    f.close()