# -*- coding: utf-8 -*-
"""
Created on Fri May 14 11:11:08 2021

@author: alondra sola

Project: turbine design tool (TFG)
Program: 1D analysis (imitating the "automobile gas generator" design tool)
File: main.py

"""

import numpy as np                   # library for math and calculation functions
import matplotlib.pyplot as plt      # library for plots
import pandas as pd                  # library for tables and data visualisation

from aux_functions import rho1_converge, blade_geometry, turbine

###################### GIVEN PARAMETERS AND ASSUMTIOSN ########################

# GIVEN PARAMETERS (TURBINE)
# P01 = 101325            # inlet total pressure [Pa]
# T01 = 288               # inlet total temperature [K]
# pi_compressor = 4       # compressor pressure ratio [-]
# NPw = 150000            # net power to wheels [W]
T01 = 1400                     # inlet total temperature (exit combustor) = T4t [K]
P01 = 397194                   # inlet total pressure (exit combustor) [Pa]
P03 = 101325                   # outlet total pressure [Pa]
T03 = 1083.529                 # outlet total temperature [K]
DeltaH_prod = 392423.1109      # enthalpy produced by turbine [J/Kg]



# FLUID PROPERTIES (TURBINE)
gamma = 1.3              # [-]
cp = 1240                # [J/kg/K]
R = 286.1538462          # [J/kg/K]


# ASSUMPTIONS
# note: 3 is rotor exit, 2 in my report
Mach3 = 0.420               # rotor/turbine exit Mach number [-]
GR = 0.32                   # reaction degree [-]
phi = 1.75                  # loading factor [-]
eta_stator = 0.888          # stator efficiency [-]
eta_rotor = 0.808           # rotor efficiency [-]
alpha2  = np.radians(75.3)  # stator angle [deg->rad]
beta3 = np.radians(-65)     # rotor angle [deg->rad]
RHT = 0.9                   # ratio hub/tip radius
mdot = 0.777482308          # total mass flow [kg/s]
rho1 = 0.98317              # inlet density [kg/m^3]
R3hub = 0.106               # radius of hub at rotor exit [m]
alpha1 = 0                  # stator inlet angle (0 if asusmed axial) [rad]







# using mach, find static pressure
P3 = P03/(1+(gamma-1)/2*Mach3**2)**(gamma/(gamma-1))
T3s = T03*(P3/P03)**((gamma-1)/gamma)

# find intermediate static pressure
P2 = GR*(P01-P3)+P3 # Siverding rp definition
P2 = P01*((GR-1)*(1-(P3/P01)**((gamma-1)/gamma) ) +1)**(gamma/(gamma-1))  # GP rp definition
### where does that definition come from? doublecheck GR definition

# isentropic evolution in stator: absolute total temperature is same
T02 = T01
T2s = T01*(P2/P01)**((gamma-1)/gamma)
V2s = np.sqrt(2*cp*(T02 - T2s))

# assume rotor efficiency
T2 = T02 - eta_stator*(T02-T2s)
V2 = np.sqrt(2*cp*(T02 - T2))

# stator kinetic loss coefficient
xi_stator = (V2s**2 - V2**2)/V2s**2
### equation is different in report and excel sheet (corregir)

# check loss inpose ??? ###
loss_stator = V2**2/V2s**2

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
tpl_stator = (P01 - P02)/(P01 - P2)

# stator pressure loss coefficient, Total pressure loss down 02
omega_stator = (P01 - P02)/(P02 - P2)

# assume an alpha2 and project velocities
V2x = V2*np.cos(alpha2)
V2u = V2*np.sin(alpha2)

# assume a loading factor and calculate peripheral speed
U2 = np.sqrt(DeltaH_prod/phi)

# velocity triangle (pithagoras)
W2u = V2u - U2
W2x = V2x
W2 = np.sqrt(W2u**2 + W2x**2)

# relative inlet angle and relative Mach number
beta2 = np.arctan(W2u/V2x)
Mach2r = W2/a2

# relative total quantities at rotor inlet
# (add relative speed to static conditions)
T02r = T2*(1 + Mach2r**2*(gamma-1)/2) # total temperature is conserved
P02r = P2*(1 + (gamma-1)/2*Mach2r**2)**(gamma/(gamma-1))


############ ROTOR #############
# in the rotor, the relative total temperature is constant (se conservan rotalpías)
T03r = T02r

# ideal outlet temperature, real outlet temperature
# assuming the expanse to P3, calculated with the assumed M3
T3s = T03r*(P3/P02r)**((gamma-1)/gamma)
T3 = T03r - eta_rotor*(T03r - T3s)

# velocities
W3 = np.sqrt(2*cp*(T03r - T3))
W3s = np.sqrt(2*cp*(T03r - T3s))

# rotor kinetic loss
xi_rotor = (W3s**2 - W3**2)/W3s**2
### denominator no coincide entre excel and JS report

# check loss inpose ??? ###
loss_rotor = W3**2/W3s**2

# density
rho3 = P3/T3/R

# speed of sound and relative mach number
a3 = np.sqrt(gamma*R*T3)
Mach3r = W3/a3

# total absolute pressure
P03r = P3*(T03r/T3)**(gamma/(gamma-1))

# Total pressure loss
tpl_rotor = (P02r - P03r)/(P02r - P3)
# equation where does it come from?

# stator pressure loss coefficient, Total pressure loss down 02
omega_rotor = (P02r - P03r)/(P03r - P3)

# assume a value for beta3

# velocity projections, axial and tangential
W3x = W3*np.cos(beta3)
W3u = W3*np.sin(beta3)

# constant radius, constant peripheral speed
U3 = U2

# tangential outlet speed
V3u = W3u + U3
V3u_euler = -DeltaH_prod/U2 + V2u

# assume axial speed constant
V3x = W3x
V3 = np.sqrt(V3x**2 + V3u**2)
alpha3 = np.arctan(V3u/V3x)
Mach3c = V3/a3
DeltaH_T1 = U3*(V2u-V3u)    # euler formulation

# compute total condtions
T03 = T3*(1+(gamma-1)/2*Mach3**2)
P03 = P3*(1+(gamma-1)/2*Mach3**2)**(gamma/(gamma-1))
DeltaH_T2 = cp*(T01-T03)

# blade geometry at 2
R2tip, R2hub, h2, R2mean, D2mean, A2 = blade_geometry(mdot, rho2, V2x, RHT)

# rotational velocity
Omega2 = U2/R2mean
RPM = 60*Omega2/2/np.pi

# massflow and characteristics at inlet
rho01 = P01/T01/R
rho1 = rho1_converge(0.9999*rho01, rho01, T01, mdot, A2, gamma, cp) # iteraciones hasta converger
V1 = np.sqrt(2*cp*T01*(1-(rho1/rho01)**(gamma-1)))
T1 = T01 - V1**2/2/cp
P1 = P01*(T1/T01)**(gamma/(gamma-1)) # isentropic
a1 = np.sqrt(T1*gamma*R)
Mach1 = V1/a1
A1 = mdot/rho1/V1

# blade geometry at 1
R1tip, R1hub, h1, R1mean, D1mean, A1 = blade_geometry(mdot, rho1, V1, RHT)


# blade geometry at 3
########## NO ENTIENDO NADAAAAA
A3 = mdot/V3x/rho3
h3 = A3/np.pi/2/R2mean
R3tip = R2mean+h3/2
D3mean = R3tip*2 - h3
R3mean = D3mean/2


Omega3 = U3/R3mean








## print data to tables
data1 = {'Variable':  ['T01', 'P01', 'V1', 'T1', 'P1', 'rho1', 'A1', 'Mach1', 'alpha1', 'mdot', 'R1tip', 'R1hub', 'h1', 'mdot'],
        'Value':     [T01, P01, V1, T1, P1, rho1, A1, Mach1, np.degrees(alpha1), mdot, R1tip, R1hub, h1, mdot],
        'Unit':      ['K', 'Pa', 'm/s', 'K', 'Pa', 'kg/m^3', 'm^2', '-', 'deg', 'kg/s', 'm', 'm', 'm', 'kg/s']}
df1 = pd.DataFrame (data1, columns = ['Variable','Value','Unit'])
print (df1)

data2 = {'Variable': ['T02', 'P02', 'V2', 'T2','P2','rho2','T2is','M2','Alpha2','V2x','V2u','W2','T02r','P02r','Beta2','W2x','W2u','Mw2','v2is','U2','Total pressure loss (stator)','A2','rh2','rt2','rm2','h2','omega'],
         'Value':    [T02, P02, V2, T2, P2,rho2,T2s,Mach2,np.degrees(alpha2),V2x,V2u,W2,T02r,P02r,np.degrees(beta2),W2x,W2u,Mach2r,V2s,U2,tpl_stator,A2,R2hub,R2tip,R2mean,h2,Omega2]}
# data2 = {'Variable': ['T02', 'P02', 'V2', 'T2','P2','rho2','T2is','M2','Alpha2','V2x','V2u','W2','T02r','P02r','Beta2','W2x','W2u','Mw2','v2is','U2','Total pressure loss (rotor?)','A2','rh2','rt2','rm2','h2','omega','H01','epsp','kinetic loss stator','optimim pitch/chord sorderberg','chord','h/c','Re dh','dh','Re c','nblades'
# ],
#          'Value':    [T02, P02, V2, T2, P2,rho2,T2s,Mach2,np.degrees(alpha2),V2x,V2u,W2,T02r,P02r,np.degrees(beta2),W2x,W2u,Mach2r,V2s,U2,tpl_rotor,A2,R2hub,R2tip,R2mean,h2,Omega2,'--','---',kinetic loss stator,optimim pitch/chord sorderberg,chord,h/c,Re dh,dh,Re c,nblades]}
df2 = pd.DataFrame (data2, columns = ['Variable','Value'])
print (df2)

data3 = {'Variable': ['M3 impose','P3 impose','T03','P03','V3','T3','P3','rho3','T3is','T3iss','M3 achieve','P3 achieve','Alpha3','V3x','V3u','W3','T03r','P03r','Beta3','W3x','W3u','Mw3','U3','Total pressure loss (rotor)','A3','rh3','rt3','rm3'],
         'Value':     [Mach3,P3,T03,P03,V3,T3,P3,rho3,T3s,'---','---','---',np.degrees(alpha3),V3x,V3u,W3,T03r,P03r,np.degrees(beta3),W3x,W3u,Mach3r,U3,tpl_rotor,A3,R2hub,R3tip,R3mean]}
df3 = pd.DataFrame (data3, columns = ['Variable','Value'])
print (df3)

turb = turbine()
turb.P.Rtip = 20
turb.rho1_init = 1.1
turb.rho01 = 1.2
turb.T01 = 1000
turb.mdot = 7
turb.A2 = 0.05
turb.gamma = 1.4
turb.cp = 1000

print(turb.P.Rtip)

# plt.plot([1, 2, 3], [P01, P02, P3], 'k--')
# plt.title("P0")


# plt.plot([1, 2, 3], [T01, T02b, P3], 'k--')
# plt.title("T0")