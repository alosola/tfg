# -*- coding: utf-8 -*-
"""
Created on Fri May 14 11:11:08 2021

@autwor: alondra sola

Project: turbine design tool (TFG)
Program: 1D analysis (imitating the "automobile gas generator" design tool)
File: main.py

"""

## IMPORT NECESSARY LIBRARIES
import numpy as np                   # library for math and calculation functions
import pandas as pd                  # library for tables and data visualisation

## IMPORT FUNCTIONS AND DEFINITIONS FROM OTHER PYTHON FILES
import aux_functions as f
import plot_functions as graphs
from definitions import plane, component
from solve_functions import solve_rhox
from converge_mach3 import converge_mach3



one = plane()
two = plane()
thr = plane()

rotor = component()
stator = component()


###################### GIVEN PARAMETERS AND ASSUMTIOSN ########################

# GIVEN PARAMETERS (TURBINE)
one.T0 = 1400                     # inlet total temperature (exit combustor) = T4t [K]
one.P0 = 397194                   # inlet total pressure (exit combustor) [Pa]
thr.P0 = 101325                   # outlet total pressure [Pa]
thr.T0 = 1083.529                 # outlet total temperature [K]
DeltaH_prod = 392423.1109         # enthalpy produced by turbine [J/Kg]


# ASSUMPTIONS
GR = 0.32                       # reaction degree [-]


# FLUID PROPERTIES (TURBINE)
gamma = 1.3              # [-]
cp = 1240                # [J/kg/K]
R = 286.1538462          # [J/kg/K]


# ASSUMPTIONS
Mach3_init = 0.50               # rotor/turbine exit Mach number [-]
phi = 1.75                      # loading factor [-]
stator.eta = 0.888              # stator efficiency [-]
rotor.eta = 0.808               # rotor efficiency [-]
two.alpha  = np.radians(75.3)   # stator angle [deg->rad]
thr.beta = np.radians(-65)      # rotor angle [deg->rad]

RHT = 0.9                       # ratio hub/tip radius
mdot = 0.777482308              # total mass flow [kg/s]
one.rho = 0.98317               # inlet density [kg/m^3]
thr.geo.Rh = 0.106              # radius of hub at rotor exit [m]
one.alpha = 0                   # stator inlet angle (0 if asusmed axial) [rad]



one, two, thr, stator, rotor = converge_mach3(Mach3_init, one, two, thr, stator, rotor, gamma, cp, R, GR, phi, DeltaH_prod)



# stator kinetic loss coefficient   # check loss inpose ??? ###   # Total pressure loss   # stator pressure loss coefficient, Total pressure loss down 02
xi_stator, loss_stator, tpl_stator, omega_stator = f.losses(two.vel.V, two.vel.Vs, one.P0, two.P0, two.P)
# = (two.vel.Vs**2 - two.vel.V**2)/two.vel.Vs**2,   = two.vel.V**2/two.vel.Vs**2,    = (one.P0 - two.P0)/(one.P0 - two.P),     = (one.P0 - two.P0)/(two.P0 - two.P)
### equation for xi is different in report and excel sheet (corregir)


# rotor kinetic loss coefficient   # check loss inpose ??? ###   # Total pressure loss   # rotor pressure loss coefficient, Total pressure loss down 02
xi_rotor, loss_rotor, tpl_rotor, omega_rotor = f.losses(thr.vel.W, thr.vel.Ws, two.P0r, thr.P0r, thr.P)
# = (two.vel.Vs**2 - two.vel.V**2)/two.vel.Vs**2,   = two.vel.V**2/two.vel.Vs**2,    = (one.P0 - two.P0)/(one.P0 - two.P),     = (one.P0 - two.P0)/(two.P0 - two.P)
### equation for xi is different in report and excel sheet (corregir)
### denominator no coincide entre excel and JS report




DeltaH_T1 = thr.vel.U*(two.vel.Vu-thr.vel.Vu)    # euler formulation

# compute total condtions
# thr.T0, thr.P0 = f.relative_temperature_pressure(thr.T, thr.P, thr.vel.M, gamma)
# thr.T0 = thr.T*(1+(gamma-1)/2*thr.vel.M**2)
# thr.P0 = thr.P*(1+(gamma-1)/2*thr.vel.M**2)**(gamma/(gamma-1))

thr.T0, thr.P0 = f.total_conditions(thr.T, thr.vel.V, thr.P, cp, gamma)

DeltaH_T2 = cp*(one.T0-thr.T0)






# blade geometry at 2
two.geo.Rt, two.geo.Rh, two.geo.h, two.geo.Rm, two.geo.Dm, A2 = f.blade_geometry(mdot, two.rho, two.vel.Vx, RHT)

# rotational velocity
two.vel.Omega = two.vel.U/two.geo.Rm
two.vel.RPM = f.RPM(two.vel.Omega)

# massflow and characteristics at inlet
one.rho0 = f.static_density(one.P0, one.T0, R)                                  # = one.P0/one.T0/R
one.rho = solve_rhox(0.9999*one.rho0, one.rho0, one.T0, mdot, A2, gamma, cp)  # iteraciones hasta converger
one.vel.V = np.sqrt(2*cp*one.T0*(1-(one.rho/one.rho0)**(gamma-1)))
one.T = one.T0 - one.vel.V**2/2/cp
one.P = f.T2P(one.P0, one.T, one.T0, gamma)                     # = one.P0*(one.T/one.T0)**(gamma/(gamma-1)) # isentropic
one.vel.a, one.vel.M = f.sonic(gamma, R, one.T, one.vel.V)
one.geo.A = mdot/one.rho/one.vel.V

# blade geometry at 1
one.geo.Rt, one.geo.Rh, one.geo.h, one.geo.Rm, one.geo.Dm, one.geo.A = f.blade_geometry(mdot, one.rho, one.vel.V, RHT)

# blade geometry at 3
########## NO ENTIENDO NADAAAAA
thr.geo.A = mdot/thr.vel.Vx/thr.rho
thr.geo.h = thr.geo.A/np.pi/2/two.geo.Rm
thr.geo.Rt = two.geo.Rm+thr.geo.h/2
thr.geo.Dm = thr.geo.Rt*2 - thr.geo.h
thr.geo.Rm = thr.geo.Dm/2

thr.vel.Omega = thr.vel.U/thr.geo.Rm

GR_enthalpy = (two.T - thr.T)/(two.T0 - thr.T0) # Assume v1 = v3, not right
Deltabeta = two.beta - thr.beta # por qué no negativo?
Deltaalpha = two.alpha - thr.alpha
GR_final = (two.T - thr.T)/(one.T - thr.T)
DeltaW = thr.vel.W - two.vel.W
h3h2 = thr.geo.h/two.geo.h


# SOODERBERG CORRELATION LOSSES
# rotor
thr.H = thr.T*cp
thr.Hs = thr.Ts*cp
epsr = 2*(thr.H - thr.Hs)/thr.vel.W**2
#stator
H2 = two.T*cp
two.Hs = two.Ts*cp
epse = 2*(H2 - two.Hs)/two.vel.V**2

# efficiencies
one.H0 = cp*one.T0
thr.H0 = cp*thr.T0
thr.Ts = one.T0*(thr.P/one.P0)**((gamma-1)/gamma)
thr.Hs = cp*thr.Ts
etaTT = (one.H0 - thr.H0)/(one.H0 - (thr.Hs+thr.vel.V**2/2))
etaTS = (one.H0 - thr.H0)/(one.H0 - thr.Hs)

epsp = (two.vel.Vs**2 - two.vel.V**2)/two.vel.V**2
epspp = (thr.vel.Ws**2 - thr.vel.W**2)/thr.vel.W**2

knls_stator = epsp*(1+gamma*two.vel.M**2/2)
knls_rotor = epspp*(1+gamma*thr.vel.Mr**2/2)



stg = np.radians(60)
opt_pitch_chord = 0.6*np.cos(stg)/2/np.cos(two.alpha)**2/(one.vel.V/two.vel.Vx*np.tan(one.alpha) + np.tan(two.alpha))*(one.geo.h/two.geo.h)*one.P0/two.P0




# ## print data to tables
# datone.vel.a = {'Variable':  ['one.T0', 'one.P0', 'one.vel.V', 'one.T', 'one.P', 'one.rho', 'one.geo.A', 'one.vel.M', 'one.alpha', 'mdot', 'one.geo.Rt', 'one.geo.Rh', 'one.geo.h', 'mdot'],
#         'Value':     [one.T0, one.P0, one.vel.V, one.T, one.P, one.rho, one.geo.A, one.vel.M, np.degrees(one.alpha), mdot, one.geo.Rt, one.geo.Rh, one.geo.h, mdot],
#         'Unit':      ['K', 'Pa', 'm/s', 'K', 'Pa', 'kg/m^3', 'm^2', '-', 'deg', 'kg/s', 'm', 'm', 'm', 'kg/s']}
# df1 = pd.DataFrame (datone.vel.a, columns = ['Variable','Value','Unit'])
# print (df1)

# dattwo.vel.a = {'Variable': ['two.T0', 'two.P0', 'two.vel.V', 'two.T','two.P','two.rho','two.Tis','M2','two.alpha','two.vel.Vx','two.vel.Vu','two.vel.W','two.T0r','two.P0r','Bettwo.vel.a','two.vel.Wx','two.vel.Wu','Mw2','v2is','two.vel.U','Total pressure loss (stator)','A2','rtwo.geo.h','rt2','rm2','two.geo.h','omega'],
#           'Value':    [two.T0, two.P0, two.vel.V, two.T, two.P,two.rho,two.Ts,two.vel.M,np.degrees(two.alpha),two.vel.Vx,two.vel.Vu,two.vel.W,two.T0r,two.P0r,np.degrees(two.beta),two.vel.Wx,two.vel.Wu,two.vel.Mr,two.vel.Vs,two.vel.U,tpl_stator,A2,two.geo.Rh,two.geo.Rt,two.geo.Rm,two.geo.h,two.vel.Omega]}
# # dattwo.vel.a = {'Variable': ['two.T0', 'two.P0', 'two.vel.V', 'two.T','two.P','two.rho','two.Tis','M2','two.alpha','two.vel.Vx','two.vel.Vu','two.vel.W','two.T0r','two.P0r','Bettwo.vel.a','two.vel.Wx','two.vel.Wu','Mw2','v2is','two.vel.U','Total pressure loss (rotor?)','A2','rtwo.geo.h','rt2','rm2','two.geo.h','omega','one.H0','epsp','kinetic loss stator','optimim pitch/chord sorderberg','chord','h/c','Re dh','dh','Re c','nblades'
# # ],
# #          'Value':    [two.T0, two.P0, two.vel.V, two.T, two.P,two.rho,two.Ts,two.vel.M,np.degrees(two.alpha),two.vel.Vx,two.vel.Vu,two.vel.W,two.T0r,two.P0r,np.degrees(two.beta),two.vel.Wx,two.vel.Wu,two.vel.Mr,two.vel.Vs,two.vel.U,tpl_rotor,A2,two.geo.Rh,two.geo.Rt,two.geo.Rm,two.geo.h,two.vel.Omega,'--','---',kinetic loss stator,optimim pitch/chord sorderberg,chord,h/c,Re dh,dh,Re c,nblades]}
# df2 = pd.DataFrame (dattwo.vel.a, columns = ['Variable','Value'])
# print (df2)

# datthr.vel.a = {'Variable': ['M3 impose','thr.P impose','thr.T0','thr.P0','thr.vel.V','thr.T','thr.P','thr.rho','thr.Tis','thr.Tiss','M3 achieve','thr.P achieve','Alphthr.vel.a','thr.vel.Vx','thr.vel.Vu','thr.vel.W','thr.T0r','thr.P0r','thr.beta','thr.vel.Wx','thr.vel.Wu','Mw3','thr.vel.U','Total pressure loss (rotor)','thr.geo.A','rthr.geo.h','rt3','rm3'],
#           'Value':     [thr.vel.M,thr.P,thr.T0,thr.P0,thr.vel.V,thr.T,thr.P,thr.rho,thr.Ts,'---',thr.vel.Mc,'---',np.degrees(thr.alpha),thr.vel.Vx,thr.vel.Vu,thr.vel.W,thr.T0r,thr.P0r,np.degrees(thr.beta),thr.vel.Wx,thr.vel.Wu,thr.vel.Mr,thr.vel.U,tpl_rotor,thr.geo.A,two.geo.Rh,thr.geo.Rt,thr.geo.Rm]}
# df3 = pd.DataFrame (datthr.vel.a, columns = ['Variable','Value'])
# print (df3)

# ORIGINAL GRAPH FROM JS STUDY
graphs.velocity_triangle(753.18, 197.59, -71.35, 254.09 ,279.64, 197.6, -544.89, 254.09)


graphs.velocity_triangle(two.vel.Vu, two.vel.Vx, thr.vel.Vu, thr.vel.Vx, two.vel.Wu, two.vel.Wx, thr.vel.Wu, thr.vel.Wx)


# plt.plot([1, 2, 3], [one.P0, two.P0, thr.P], 'k--')
# plt.title("P0")


# plt.plot([1, 2, 3], [one.T0, two.T0b, thr.P], 'k--')
# plt.title("T0")
