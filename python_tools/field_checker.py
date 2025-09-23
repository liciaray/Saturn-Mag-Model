#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 26 January 2025

@author: liciaray
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.colors import Normalize


field_line = np.loadtxt('../output/event_2006-07-12T07h23m.txt', skiprows=1)

x_field = field_line[:,0]
y_field = field_line[:,1]
z_field = field_line[:,2]
r_field = field_line[:,3]
s_field = field_line[:,4]
b_mag   = field_line[:,5]


#make a oblate sphere to represent saturn
u, v = np.mgrid[0:2*np.pi:180j,0:np.pi:90j]
r_sphere = 1.01 - 1./11.1*np.abs(np.cos(v))**2
x_sphere = r_sphere*np.cos(u)*np.sin(v)
y_sphere = r_sphere*np.sin(u)*np.sin(v)
z_sphere = r_sphere*np.cos(v)

#make some rings!
ring_radius = 3  # Adjust this based on your sphere's size
ring_radius, ring_theta = np.mgrid[1.2:2.2:20j,0:2*np.pi:180j]
ring_x = ring_radius * np.cos(ring_theta)
ring_y = ring_radius * np.sin(ring_theta)
ring_z = np.zeros_like(ring_x)


# Create a 3D scatter plot
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Normalize the magnetic field data to use it for color scaling
norm = Normalize(vmin=np.min(np.log(b_mag)), vmax=np.max(np.log(b_mag)))
cmap = cm.viridis  # You can choose any colormap like 'plasma', 'inferno', etc.


# Create the scatter plot
sc = ax.scatter(x_field, y_field, z_field, c=np.log(b_mag), cmap=cmap, norm=norm)
sc1 = ax.plot_surface(x_sphere, y_sphere, z_sphere, alpha = 0.8, cmap='YlOrBr')

# Add color bar
cbar = plt.colorbar(sc)
cbar.set_label('Log(B$_{Mag}$)')

# Plot the ring
sc2=ax.plot(ring_x, ring_y, ring_z, color='#B0B0B0', lw=1, alpha=0.5)

# Set labels
ax.set_xlabel('X$_{KSM}$')
ax.set_ylabel('Y$_{KSM}$')
ax.set_zlabel('Z$_{KSM}$')
ax.set_title('Magnetic Field Trace from Cassini to Footprint at Saturn')

#Find the range of the field line to set a range for the scatter
#plot that will make an evenly spaced cube
max_range = np.array([x_field.max()-x_field.min(), y_field.max()-y_field.min(), z_field.max()-z_field.min()]).max()

mid_x = (x_field.max() + x_field.min()) * 0.5
mid_y = (y_field.max() + y_field.min()) * 0.5
mid_z = (z_field.max() + z_field.min()) * 0.5


ax.set_xlim(mid_x - max_range/2, mid_x + max_range/2)
ax.set_ylim(mid_y - max_range/2, mid_y + max_range/2)
ax.set_zlim(mid_z - max_range/2, mid_z + max_range/2)


plt.savefig("OddPoint.png")

# Show the plot
plt.show(block=True)
plt.pause(1)
