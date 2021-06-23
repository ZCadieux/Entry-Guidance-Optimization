from openmdao.api import Problem
from omtools.api import Group
import omtools.api as ot
import numpy as np

from entry_functions import EntryFunctions
from entry_system import EntrySystem

# Test Case 1
if 1:
    alt = 423819.5625 # geodetic altitude (ft)
    lon = -175.80 # longitude (deg)
    lat = -22.4 # geodetic latitude (deg)
    vel = 4694.1 # relative velocity (m/s)
    fpa = -10.3228 # relative flight path angle (deg)
    hea = -0.0522 # relative heading angle (deg)

 

    alt = alt * 0.3048 # geodetic altitude (m)
    if lon < 0.0:
        lon = lon + 360.0

 

# Test Case 2
if 0:
    alt = 125000.0  # geodetic altitude (m)
    lon = -176.40167 # longitude (deg)
    lat = -21.3 # geodetic latitude (deg)
    vel = 4705.6 # relative velocity (m/s)
    fpa = -9.9706 # relative flight path angle (deg)
    hea = -2.8758 # relative heading angle (deg)

 

    if lon < 0.0:
        lon = lon + 360.0
# =====================================


#case for testing entry_functions on its own
num_nodes = 1
R0 = 3.39619e6
flatten = 0.005886007555526
area = 159.94 #from vehicle data

prob1 = Problem(EntryFunctions())
prob1.setup()
prob1.set_val('r', alt)
prob1.set_val('V', vel)
prob1.run_model()

print('')
print('C_L: ', prob1['C_L'])
print('C_D: ', prob1['C_D'])
print('lift: ', prob1['lift'])
print('drag: ', prob1['drag'])
print('density: ', prob1['density'])
print('Vsound: ', prob1['Vsound'])

D = prob1['drag']
L = prob1['lift']

#Dynamics
prob = Problem(EntrySystem())
prob.setup()
prob.set_val('r', alt)
#prob.set_val('theta', lon) dynamics does not rely on longitude
prob.set_val('phi', lat)
prob.set_val('V', vel)
prob.set_val('gamma', fpa)
prob.set_val('psi', hea)

#control variable, needs to be set with control algorithm (set to 0 for test case)
prob.set_val('sigma', 0)
prob.run_model()

print('')
print('dr_dt: ', prob['dr_dt'])
print('dtheta_dt: ', prob['dtheta_dt'])
print('dphi_dt: ', prob['dphi_dt'])
print('dV_dt: ', prob['dV_dt'])
print('dgamma_dt: ', prob['dgamma_dt'])
print('dpsi_dt: ', prob['dpsi_dt'])