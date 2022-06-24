#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 11:46:13 2022

@author: photon-zb
"""
import meep as mp
import math
import cmath

import matplotlib.pyplot as plt
import matplotlib.patches as patches

# from meep.materials import Si

### ----------------------------------------------------------------
# Resolution in simulation
resolution = 200     # pixels/um

# Input parameters
# Source parameters
lmin = 0.4        # source min wavelength
lmax = 0.8        # source max wavelength
fmin = 1/lmax       # source min frequency
fmax = 1/lmin       # source max frequency
fcen = 0.5*(fmin+fmax) # center frequency
df = fmax-fmin     # Frequency span
theta_in = 0    # Plane wave incident angle
theta = math.radians(theta_in)  # convert degree to radian

# Structure parameters
ax = 0.1        # lattice periodicity
#ay = 3**0.5*ax
ay = 0.1
Bw = 0.09
Tw = 0.05

h = 0.5         # nanocylinder height

#tsub = 1.0          # substrate thickness
#tabs = 1.0          # PML thickness
#tair = 1.0          # air thickness
tsub = .1          # substrate thickness
tabs = .10          # PML thickness
tair = .10          # air thickness
sz = tabs+tair+h+tsub+tabs # Total length along the z-axis

# Defined Cell for Calculation
cell_size = mp.Vector3(ay,ax,sz)  # Unit cell


# material parameter
Cyl =  mp.Medium(index=1.49)
Sub = mp.Medium(index=3.5)

# PML Layer
pml_layers = [mp.PML(thickness=tabs,direction=mp.Z,side=mp.High),
              mp.Absorber(thickness=tabs,direction=mp.Z,side=mp.Low)]
  
geometry = [ mp.Sphere(material=Cyl, radius=Tw/2, 
                         center=mp.Vector3(0,0,0.5*sz-tabs-tair)),
            mp.Cone(material=Cyl, radius=Bw/2, radius2=Tw/2, height=h, 
                     center=mp.Vector3(0,0,0.5*sz-tabs-tair-0.5*h)),
                mp.Block(material=Sub, size=mp.Vector3(mp.inf,mp.inf,tsub+tabs),
                              center=mp.Vector3(0,0,0.5*sz-tabs-tair-h-0.5*(tsub+tabs))) ]   
                                   
## Geometry Plot
sim = mp.Simulation(resolution=resolution,cell_size=cell_size,geometry=geometry)
sim.init_sim()
eps_data = sim.get_epsilon()

import numpy as np
II = np.asarray(eps_data)
xo,yo,zo = np.shape(II)
print("Geometry shape:",xo,yo,zo)

#fig, (ax0,ax1, ax2) = plt.subplots(1, 3)
#ax0.imshow(II[:,:,int(zo-(tabs+tair+0.5*h)*resolution)])
#ax0.set_title("XY plane in pixel unit")
#ax0.set_xlabel("X"); ax0.set_ylabel("Y")
#
#ax1.imshow(np.flipud(np.transpose(II[int(xo/2),:,:])),aspect ='auto')
#ax1.set_title("YZ plane at X=0 in pixel unit")
#ax1.set_xlabel("Y"); ax1.set_ylabel("Z")
#
#ax2.imshow(np.flipud(np.transpose(II[:,int(yo/2),:])),aspect ='auto')
#ax2.set_title("XZ plane at Y=0 in pixel unit")
#ax2.set_xlabel("X"); ax2.set_ylabel("Z")
#plt.show()

from mpl_toolkits.axes_grid1 import make_axes_locatable
x = np.linspace(-ay/2,ay/2,xo)
y = np.ones(xo)
fig, (ax0,ax1, ax2) = plt.subplots(1, 3)
im0=ax0.imshow(II[:,:,int(zo-(tabs+tair+0.5*h)*resolution)],extent=[-ax/2,ax/2,-ay/2,ay/2])
ax0.set_title("XY plane in pixel unit")
ax0.set_xlabel("X (um)"); ax0.set_ylabel("Y (um)")
#divider = make_axes_locatable(ax0)
#cax = divider.append_axes('right', size='5%', pad=0.05)
#fig.colorbar(im0, cax=cax, orientation='vertical')
plt.colorbar(im0, ax=ax0)

ax1.imshow(np.flipud(np.transpose(II[int(xo/2),:,:])),aspect ='auto',
           extent=[-ax/2,ax/2,-sz/2,sz/2],)
ax1.set_title("YZ plane at X=0")
ax1.set_xlabel("Y (um)"); ax1.set_ylabel("Z (um)")

im2=ax2.imshow(np.flipud(np.transpose(II[:,int(yo/2),:])),aspect ='auto',
           extent=[-ay/2,ay/2,-sz/2,sz/2])
ax2.set_title("XZ plane at Y=0")
ax2.set_xlabel("X (um)"); ax2.set_ylabel("Z (um)")
plt.colorbar(im2, ax=ax2)

ax2.plot(x,y*(0.5*sz-tabs-0.2*tair),'r')
ax2.plot(x,y*(0.5*sz-tabs-0.4*tair),'g')#,linewidth='3')
#ax2.plot(x,y*(0.5*sz-tabs-0.6*tair),'g:')#,linewidth='3')
ax2.plot(x,y*(-0.5*sz+tabs+0.2*tsub),'b')

# Create PML as a Rectangle patch
rect_top = patches.Rectangle((-ay/2, 0.5*sz-tabs), ay, tabs, edgecolor='none', facecolor='r')
ax2.add_patch(rect_top)
rect_bot = patches.Rectangle((-ay/2, -0.5*sz), ay, tabs, edgecolor='none', facecolor='r')
ax2.add_patch(rect_bot)

plt.show()

