# -*- coding: utf-8 -*-
"""
Created on Fri May 14 11:11:08 2021

@author: alondra sola

Project: turbine design tool (TFG)
Program: 1D analysis (imitating the "automobile gas generator" design tool)
File: main.py

"""

## IMPORT NECESSARY LIBRARIES
import numpy as np                   # library for math and calculation functions
import matplotlib.pyplot as plt      # library for plots
import pandas as pd                  # library for tables and data visualisation

## IMPORT FUNCTIONS AND DEFINITIONS FROM OTHER PYTHON FILES
import aux_functions as f
import plot_functions as graphs
from definitions import component









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
P3 = f.static_pressure(P03, gamma, Mach3)    # = P03/(1+(gamma-1)/2*Mach3**2)**(gamma/(gamma-1))
T3s = f.P2T(T03, P3, P03, gamma)             # = T03*(P3/P03)**((gamma-1)/gamma)

# find intermediate static pressure
P2 = GR*(P01-P3)+P3 # Siverding rp definition
P2 = P01*((GR-1)*(1-(P3/P01)**((gamma-1)/gamma) ) +1)**(gamma/(gamma-1))  # GP rp definition
### where does that definition come from? doublecheck GR definition

# isentropic evolution in stator: absolute total temperature is same
T02 = T01
T2s = f.P2T(T01, P2, P01, gamma)             # =  T01*(P2/P01)**((gamma-1)/gamma)
V2s = f.isen_velocity(cp, T02, T2s)          # = np.sqrt(2*cp*(T02 - T2s))

# assume rotor efficiency
T2 = f.stage_efficiency(eta_stator, T02, T2s)    # = T02 - eta_stator*(T02-T2s)
V2 = f.isen_velocity(cp, T02, T2)                # = np.sqrt(2*cp*(T02 - T2))


# density (from static quantities)
rho2 = f.static_density(P2,T2,R)                 # = P2/T2/R

# speed of sound and Mach
a2, Mach2 = f.sonic(gamma, R, T2, V2)            # = np.sqrt(gamma*R*T2),    = V2/a2

# check if pressure calculated meets the restriction
T02b = T2+V2**2/2/cp
### include this in the check loop

P02 = f.T2P(P2, T02, T2, gamma)                  # = P2*(T02/T2)**(gamma/(gamma-1))

# stator kinetic loss coefficient   # check loss inpose ??? ###   # Total pressure loss   # stator pressure loss coefficient, Total pressure loss down 02
xi_stator, loss_stator, tpl_stator, omega_stator = f.losses(V2, V2s, P01, P02, P2)
# = (V2s**2 - V2**2)/V2s**2,   = V2**2/V2s**2,    = (P01 - P02)/(P01 - P2),     = (P01 - P02)/(P02 - P2)
### equation for xi is different in report and excel sheet (corregir)

# assume an alpha2 and project velocities
V2x, V2u = f.velocity_projections(V2, alpha2)

# assume a loading factor and calculate peripheral speed
U2 = np.sqrt(DeltaH_prod/phi)

# velocity triangle (pithagoras)
W2u = V2u - U2
W2x = V2x
W2 = f.mag(W2u, W2x)

# relative inlet angle and relative Mach number
beta2 = np.arctan(W2u/V2x)
Mach2r = W2/a2

# relative total quantities at rotor inlet
# (add relative speed to static conditions)
# T02r = T2*(1 + Mach2r**2*(gamma-1)/2)              # total temperature is conserved
# P02r = P2*(1 + (gamma-1)/2*Mach2r**2)**(gamma/(gamma-1))
T02r, P02r = f.relative_temperature_pressure(T2, P2, Mach2r, gamma)

############ ROTOR #############
# in the rotor, the relative total temperature is constant (se conservan rotalpías)
T03r = T02r

# ideal outlet temperature, real outlet temperature
# assuming the expanse to P3, calculated with the assumed M3
T3s = f.P2T(T03r, P3, P02r, gamma)                  # = T03r*(P3/P02r)**((gamma-1)/gamma)
T3 = f.stage_efficiency(eta_rotor, T03r, T3s)       # = T03r - eta_rotor*(T03r - T3s)

# velocities
W3 = f.isen_velocity(cp, T03r, T3)                  # = np.sqrt(2*cp*(T03r - T3))
W3s = f.isen_velocity(cp, T03r, T3s)                # = np.sqrt(2*cp*(T03r - T3s))

# density
rho3 = f.static_density(P3,T3,R)                    # = P3/T3/R

# speed of sound and relative mach number
a3, Mach3r = f.sonic(gamma, R, T3, W3)                      # = np.sqrt(gamma*R*T3)

# total absolute pressure
P03r = f.T2P(P3, T03r, T3, gamma)                           # = P3*(T03r/T3)**(gamma/(gamma-1))

# rotor kinetic loss coefficient   # check loss inpose ??? ###   # Total pressure loss   # rotor pressure loss coefficient, Total pressure loss down 02
xi_rotor, loss_rotor, tpl_rotor, omega_rotor = f.losses(W3, W3s, P02r, P03r, P3)
# = (V2s**2 - V2**2)/V2s**2,   = V2**2/V2s**2,    = (P01 - P02)/(P01 - P2),     = (P01 - P02)/(P02 - P2)
### equation for xi is different in report and excel sheet (corregir)
### denominator no coincide entre excel and JS report


# assume a value for beta3

# velocity projections, axial and tangential
W3x, W3u = f.velocity_projections(W3, beta3)

# constant radius, constant peripheral speed
U3 = U2

# tangential outlet speed
V3u = W3u + U3
V3u_euler = -DeltaH_prod/U2 + V2u

# assume axial speed constant
V3x = W3x
V3 = f.mag(V3x, V3u)                           # = np.sqrt(V3x**2 + V3u**2)
alpha3 = np.arctan(V3u/V3x)
Mach3c = V3/a3
DeltaH_T1 = U3*(V2u-V3u)    # euler formulation

# compute total condtions
T03, P03 = f.relative_temperature_pressure(T3, P3, Mach3, gamma)
# T03 = T3*(1+(gamma-1)/2*Mach3**2)
# P03 = P3*(1+(gamma-1)/2*Mach3**2)**(gamma/(gamma-1))
DeltaH_T2 = cp*(T01-T03)






# blade geometry at 2
R2tip, R2hub, h2, R2mean, D2mean, A2 = f.blade_geometry(mdot, rho2, V2x, RHT)

# rotational velocity
Omega2 = U2/R2mean
RPM = f.RPM(Omega2)

# massflow and characteristics at inlet
rho01 = f.static_density(P01, T01, R)                                  # = P01/T01/R
rho1 = f.rhox_converge(0.9999*rho01, rho01, T01, mdot, A2, gamma, cp)  # iteraciones hasta converger
V1 = np.sqrt(2*cp*T01*(1-(rho1/rho01)**(gamma-1)))
T1 = T01 - V1**2/2/cp
P1 = f.T2P(P01, T1, T01, gamma)                     # = P01*(T1/T01)**(gamma/(gamma-1)) # isentropic
a1, Mach1 = f.sonic(gamma, R, T1, V1)
A1 = mdot/rho1/V1

# blade geometry at 1
R1tip, R1hub, h1, R1mean, D1mean, A1 = f.blade_geometry(mdot, rho1, V1, RHT)

# blade geometry at 3
########## NO ENTIENDO NADAAAAA
A3 = mdot/V3x/rho3
h3 = A3/np.pi/2/R2mean
R3tip = R2mean+h3/2
D3mean = R3tip*2 - h3
R3mean = D3mean/2

Omega3 = U3/R3mean

GR_enthalpy = (T2 - T3)/(T02 - T03) # Assume v1 = v3, not right
Deltabeta = abs(beta2 - beta3) # por qué no negativo?
DeltaW = W3 - W2
heightratio = h3/h2




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
          'Value':     [Mach3,P3,T03,P03,V3,T3,P3,rho3,T3s,'---',Mach3c,'---',np.degrees(alpha3),V3x,V3u,W3,T03r,P03r,np.degrees(beta3),W3x,W3u,Mach3r,U3,tpl_rotor,A3,R2hub,R3tip,R3mean]}
df3 = pd.DataFrame (data3, columns = ['Variable','Value'])
print (df3)


graphs.velocity_triangle(V2u, V2x, V3u, V3x, W2u, W2x, W3u, W3x)


# plt.plot([1, 2, 3], [P01, P02, P3], 'k--')
# plt.title("P0")


# plt.plot([1, 2, 3], [T01, T02b, P3], 'k--')
# plt.title("T0")