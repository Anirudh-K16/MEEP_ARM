#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 11:46:13 2022

@author: photon-zb
"""
import meep as mp
import math
import cmath
import numpy as np

import matplotlib.pyplot as plt
from meep.materials import Al

### ----------------------------------------------------------------
# Resolution in simulation
resolution = 30     # pixels/um
autoshut = 1E-4

# Input parameters
# Source parameters
lmin = 0.4      # source min wavelength
lmax = 0.8      # source max wavelength
fmin = 1/lmax       # source min frequency
fmax = 1/lmin       # source max frequency
fcen = 0.5*(fmin+fmax) # center frequency
df = fmax-fmin     # Frequency span
theta_in = 0    # Plane wave incident angle
theta = math.radians(theta_in)  # convert degree to radian

# Structure parameters
a = 0.45         # lattice periodicity
r = 0.1       # nanocylinder radius
h = 0.15            # nanocylinder height
tsub = 4.0          # substrate thickness
tabs = 5.0          # PML thickness
tair = 4.0          # air thickness
sz = tabs+tair+h+tsub+tabs # Total length along the z-axis

# Defined Cell for Calculation
cell_size = mp.Vector3(a,a,sz)  # Unit cell


# material parameter
#Sub =  mp.Medium(index=1,D_conductivity=10)
Cyl =  mp.Medium(index=2.5)
Sub = mp.Medium(index=1.5)

# PML Layer
pml_layers = [mp.PML(thickness=tabs,direction=mp.Z,side=mp.High),
              mp.Absorber(thickness=tabs,direction=mp.Z,side=mp.Low)]

# Structure Geometry 
geometry = [ mp.Cylinder(material=Cyl, radius=r, height=h, 
                         center=mp.Vector3(0,0,0.5*sz-tabs-tair-0.5*h)),
            mp.Block(material=Sub, size=mp.Vector3(mp.inf,mp.inf,tsub+tabs),
                center=mp.Vector3(0,0,0.5*sz-tabs-tair-h-0.5*(tsub+tabs))) ]
            


 # k with correct length (plane of incidence: XZ) 

k = mp.Vector3(math.sin(theta),0,math.cos(theta))*2*math.pi/(lmax+lmin)/2


src_pos = 0.5*sz-tabs-0.4*tair
sources = [mp.Source(mp.GaussianSource(fcen,fwidth=df),component=mp.Ey,
                     center=mp.Vector3(0,0,src_pos),size=mp.Vector3(a,a,0))]    

sim = mp.Simulation(cell_size=cell_size,
#     geometry=geometry,
     sources=sources,
     boundary_layers=pml_layers,
     k_point = k,
     resolution=resolution)

nfreq = 4001
 # reflected flux
refl = sim.add_flux(fcen, df, nfreq, mp.FluxRegion(center=mp.Vector3(0,0,0.5*sz-tabs-0.2*tair),size=mp.Vector3(a,a,0)))
tran = sim.add_flux(fcen, df, nfreq, mp.FluxRegion(center=mp.Vector3(0,0,-0.5*sz+tabs+0.2*tsub),size=mp.Vector3(a,a,0)))


sim.run(until_after_sources=mp.stop_when_fields_decayed(25, mp.Ey, mp.Vector3(0,0,0.5*sz-tabs-0.6*tair), autoshut))



# for normalization run, save flux fields data for reflection plane
# initial / without structure
Rin_flux = sim.get_flux_data(refl)
Rin = mp.get_fluxes(refl)      # Reflection data
Tin = mp.get_fluxes(tran)     # Transmission


#
sim.reset_meep()

#-------------------- Main Calculation
sim = mp.Simulation(cell_size=cell_size,
     geometry=geometry,
     sources=sources,
     boundary_layers=pml_layers,
     k_point = k,
     resolution=resolution)


 # reflected flux
refl = sim.add_flux(fcen, df, nfreq, mp.FluxRegion(center=mp.Vector3(0,0,0.5*sz-tabs-0.2*tair),size=mp.Vector3(a,a,0)))
tran = sim.add_flux(fcen, df, nfreq, mp.FluxRegion(center=mp.Vector3(0,0,-0.5*sz+tabs+0.2*tsub),size=mp.Vector3(a,a,0)))

# for normal run, load negated fields to subtract incident from refl. fields
sim.load_minus_flux_data(refl,Rin_flux)

sim.run(until_after_sources=mp.stop_when_fields_decayed(25, mp.Ey, mp.Vector3(0,0,0.5*sz-tabs-0.6*tair), autoshut))


fs = mp.get_flux_freqs(refl)
Rf = mp.get_fluxes(refl)      # Reflection data
Tr = mp.get_fluxes(tran)     # Transmission

R  = np.asarray(Rf)
T  = np.asarray(Tr)
f  = np.asarray(fs)

Ri = np.asarray(Rin)
Ti = np.asarray(Tin)
np.savetxt('test1.txt', (f,Ri,Ti,R,T)) 



R = -R/Ti
T = T/Ti

Lm = 1000/f
plt.figure(10);
plt.plot(Lm,R,'r',label="Reflection")
plt.plot(Lm,T,'b',label="Transmission")
plt.plot(Lm,R+T,'g:',label="Total")
plt.xlabel("Wavelength (nm)")
plt.legend()
plt.show()


