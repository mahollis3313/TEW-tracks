#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 16:22:31 2020

@author: margaret

code for calculating things about wave propagation:
    the total distance a wave tracks and the distance between points, both in kilometers
    the speed of the wave, both between points and an overall average
        NOTE: Timestep must be manually changed. Originally set to 3 hours, in seconds, for MERRA-2
        
9 March 2019: fixed the units finally on the wave speeds, realized that the distances are actually fine
"""
# import packages
import numpy as np
from math import radians, degrees, sin, cos, asin, acos, sqrt
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import sys
import glob
import seaborn as sns

# netCDF stuff
import netCDF4
from netCDF4 import Dataset

#things we will need
level = int(sys.argv[1])
timestep = 3600*3 

# that handy distance function
# found at https://medium.com/@petehouston/calculate-distance-of-two-locations-on-earth-using-python-1501b1944d97

def great_circle(lon1, lat1, lon2, lat2):
    '''

    Parameters
    ----------
    lon1 : float
        The longitude of the first point, in degrees
    lat1 : float
        The latitude of the first point, in degrees
    lon2 : float
        The longitude of the second point, in degrees
    lat2 : float
        The latitude of the second point, in degrees

    Returns
    -------
    float
        The distance between the two points, in kilometers

    '''
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    return 6371 * (acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2)))

# make list of data file paths

trackfiles = sorted(glob.glob(f'/data3/TRACK/TRACKOutput/MERRA2/FilteredNoTCUpdated/{level}hPa/post_filter.*.nc'))

for trackfile in trackfiles:
    tracks = Dataset(trackfile, mode='r+')

    TrackID = tracks['TRACK_ID'][:]
    FirstPt = tracks['FIRST_PT'][:]
    NumPts  = tracks['NUM_PTS'][:]

    #the tracks themselves
    lats = tracks['latitude'][:]
    lons = tracks['longitude'][:]
 
    # lists for the things to go into
    net_travel_2d = []
    all_travel_2d = []
    #net_travel_1d = []
    #all_travel_1d = []
    all_speed = []
    avg_speed = []
    
    # # # # # # # # # # # # # # # #
    #find location of max intensity
    # # # # # # # # # # # # # # # #

    for i in range(len(TrackID)):
        firstpt = FirstPt[i]
        tracklons = lons[firstpt:(firstpt+NumPts[i])]
        tracklats = lats[firstpt:(firstpt+NumPts[i])]
        distances2d = []
        #distances1d = []
                
        for k in range(len(tracklons)-1):
            if k == 0:
                net2d = 0
                distances2d.append(0)
                #distances1d.append(0)
                #startlat = tracklats[0]
            travel2d = great_circle(tracklons[k], tracklats[k], tracklons[k+1], tracklats[k+1]) #THIS IS IN KILOMETERS
            #travel1d = great_circle(tracklons[k], startlat, tracklons[k+1], startlat)
            distances2d.append(travel2d)
            net2d += travel2d
            
        net_travel_2d.append(net2d)
        all_travel_2d.extend(distances2d)
        
        distances2d = np.asarray(distances2d)
        speeds = distances2d*1000/timestep #units are m/s, thanks to multiplying by 1000
        all_speed.extend(speeds)
        avg_speed.append(np.average(speeds))
        
    # # # # # # # #
    # save the data
    # # # # # # # #
    
    # going to be fancy here and save the things as additional variables
    # rather than new files, so I don't have to pick where to put the files
    # and because this could be useful for seasonal breakdowns and other stuff
    
    #create new variables in which to store these distances and speeds
    out_net_travel_2d = tracks.createVariable('net_dist_2d', 'f8', ('tracks'))
    out_all_travel_2d = tracks.createVariable('all_dist_2d', 'f8', ('record'))
    out_speed         = tracks.createVariable('speed', 'f8', ('record'))
    out_avg_speed     = tracks.createVariable('avg_spd', 'f8', ('tracks'))
    
    #give the new variables some attributes
    out_net_travel_2d.standard_name = '2d net travel'
    out_net_travel_2d.long_name = 'Net Travel Distance'
    out_net_travel_2d.units = 'km'
    
    out_all_travel_2d.standard_name = '2d travel'
    out_all_travel_2d.long_name = '2d travel between points'
    out_all_travel_2d.units = 'km'
    
    out_speed.standard_name = 'speed'
    out_speed.long_name = 'wave propagation speed  over preceding timestep'
    out_speed.units = 'm s**-1'
    
    out_avg_speed.standard_name = 'avg speed'
    out_avg_speed.long_name = 'average wave propagation speed over preceding timestep'
    out_avg_speed.units = 'm s**-1'
    
    #now actually write the things
    tracks['net_dist_2d'][:] = net_travel_2d
    tracks['all_dist_2d'][:] = all_travel_2d
    tracks['speed'][:]     = all_speed
    tracks['avg_spd'][:] = avg_speed
    
    #and close the file before we restart the loop
    tracks.close()
