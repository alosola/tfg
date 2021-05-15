# -*- coding: utf-8 -*-
"""
Created on Fri May 14 11:11:08 2021

@author: alondra sola

Project: turbine design tool (TFG)
Program: 1D analysis (imitating the "automobile gas generator" design tool)
File: main.py

"""

import numpy as np                   # library for math and calculation functions
# import matplotlib.pyplot as plt      # library for plots


###################### GIVEN PARAMETERS AND ASSUMTIOSN ########################

# GIVEN PARAMETERS (TURBINE)
# P01 = 101325            # inlet total pressure [Pa]
# T01 = 288               # inlet total temperature [K]
# pi_compressor = 4       # compressor pressure ratio [-]
# T03 = 1400              # exit combustor temperature [K]
# NPw = 150000            # net power to wheels [W]
T03 = 1400                # inlet total temperature = T4t [K]
P03 = 397194              # inlet total pressure [Pa]
P05 = 101325              # outlet total pressure [Pa]
T05 = 1083.529            # outlet total temperature [K]



# FLUID PROPERTIES (TURBINE)
gamma = 1.3              # [-]
cp = 1240                # [J/kg/K]
R = 286.1538462          # [J/kg/K]


# ASSUMPTIONS
# note: 3 is rotor exit, 2 in my report
Mach3 = 0.420              # rotor/turbine exit Mach number [-]
GR = 0.32                  # reaction degree [-]
phi = 1.75                  # loading factor [-]
eta_stator = 0.888          # stator efficiency [-]
eta_rotor = 0.808           # rotor efficiency [-]
alpha2  = np.radians(75.3)  # stator angle [deg->rad]
beta3 = np.radians(-65)     # rotor angle [deg->rad]
RHT = 0.9                   # ratio hub/tip radius




# using mach, find static pressure
P3 = P05/(1+(gamma-1)/2*Mach3**2)**(gamma/(gamma-1))
T3s = T05*(P3/P05)**((gamma-1)/gamma)

# find intermediate static pressure
P2 = GR*(P03-P3)+P3 # Siverding rp definition
P2 = P03*((GR-1)*(1-(P3/P03)**((gamma-1)/gamma) ) +1)**(gamma/(gamma-1))  # GP rp definition
### where does that definition come from? doublecheck GR definition

# isentropic evolution in stator: absolute total temperature is same
T02 = T03
T2s = T03*(P2/P03)**((gamma-1)/gamma)
V2s = np.sqrt(2*cp*(T02 - T2s))

# assume rotor efficiency
T2 = T02 - eta_stator*(T02-T2s)
V2 = np.sqrt(2*cp*(T02 - T2))

# stator kinetic loss coefficient
xi_stator = (V2s**2 - V2**2)/V2s**2
### equation is different in report and excel sheet

# check loss inpose ??? ###
# ck = V2**2/V2s**2

# density (from static quantities)
rho2 = P2/T2/R

# speed of sound and Mach
a2 = np.sqrt(gamma*R*T2)
Mach2 = V2/a2

# check if pressure calculated meets the restriction
T02b = T2+V2**2/2/cp
### include this in the check loop

P02 = P2*(T02/T2)**(gamma/(gamma-1))

# Total pressure loss
tpl = (P03 - P02)/(P03 - P2)

# stator pressure loss coefficient, Total pressure loss down 02
omega_stator = (P03 - P02)/(P02 - P2)
