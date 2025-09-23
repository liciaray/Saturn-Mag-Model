#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 27 January 2025

@author: liciaray
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from astropy import units as u
from astropy.coordinates.representation import CartesianRepresentation, SphericalRepresentation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.colors import Normalize


df = pd.read_csv('../output/aargh_part2.txt', sep=r'\s+')

cassini_location = np.sqrt(df.x_Cassini**2 + df.y_Cassini**2 + df.z_Cassini**2)

car = CartesianRepresentation(df.x_footprint, df.y_footprint, df.z_footprint)
spherical = car.represent_as(SphericalRepresentation)



def plot_latlon_polar(latitudes, longitudes, data_lat, data_lon, cass_loc):
    # Create the figure and polar subplot
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    
    # Define the number of latitude rings and longitude radial lines
    num_latitude_rings = len(latitudes)
    num_longitude_lines = len(longitudes)
    
    # Plot latitude rings (as circles)
    for lat in latitudes:
        # Convert latitude to radial distance (0° corresponds to the edge of the circle, 90° at center)
        radius = 90-lat  # Latitude is 0° at equator, 90° at the pole
        ax.plot(np.linspace(0, 2*np.pi, 100), [radius] * 100, color='blue', lw=0.5)  # Circle
    
    # Plot longitude lines (as radial lines)
    for lon in longitudes:
        ax.plot([np.radians(lon), np.radians(lon)], [0, 90], color='red', lw=0.5)  # Radial lines

        
    data_radius = 90 - np.degrees(data_lat)
    sc = ax.scatter(data_lon, data_radius, c=cass_loc, cmap='viridis')
   


    # Set plot limits and labels
    ax.set_ylim(0, 90)  # Limit radius to the northern hemisphere
    ax.set_rticks([0, 30, 60, 90])  # Latitudes as ticks (in degrees)
    ax.set_yticklabels([90, 60, 30, 0])  # Labels for latitude rings
    
    ax.set_xticks(np.radians(longitudes))  # Set angular ticks for longitude
    ax.set_xticklabels(longitudes)  # Labels for longitude lines
    cbar = plt.colorbar(sc,ax=ax)

    
    # Add labels and title
    ax.set_title("Cassini Footprint Locations with Latitude and Longitude", va='bottom')
    plt.show()


latitudes = np.linspace(0,90,10)
longitudes = np.linspace(0,330,12)

north = spherical.lat.value>0
south = spherical.lat.value<0

plot_latlon_polar(latitudes,longitudes, spherical.lat.value[north], spherical.lon.value[north], cassini_location.values.astype(int)[north])

plot_latlon_polar(latitudes,longitudes, np.abs(spherical.lat.value[south]), spherical.lon.value[south], cassini_location.values.astype(int)[south] )
