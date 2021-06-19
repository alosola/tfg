# -*- coding: utf-8 -*-
"""
Created on Fri May 14 11:11:08 2021

@author: alondra sola

Project: turbine design tool (TFG)
Program: inclusion of losses for 1D turbine
File: main.py

"""

## IMPORT NECESSARY LIBRARIES
import numpy as np                   # library for math and calculation functions

## IMPORT FUNCTIONS AND DEFINITIONS FROM OTHER PYTHON FILES
import aux_functions as f
import plot_functions as graphs
import tables_results as tab
from definitions import plane, component
from solve_functions import solve_rhox
from converge_mach3 import converge_mach3


# initialize class variables
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
DeltaH_prod = 390000              # enthalpy produced by turbine [J/Kg]


# ASSUMPTIONS
GR = 0.4                            # reaction degree [-]
psi = 1.75                           # loading factor [-]
eta_stator_init = 0.888              # stator efficiency [-]
eta_rotor_init = 0.808               # rotor efficiency [-]

# upper and lower bounds for alpha2 and beta3
# in this format: [alpha2_min, beta3_min], [alpha2_max, beta3_max]
bounds_angles = ([np.radians(70),np.radians(-70)], [np.radians(75), np.radians(-60)])


# FLUID PROPERTIES (TURBINE)
gamma = 1.3              # [-]
cp = 1240                # [J/kg/K]
R = 286.1538462          # [J/kg/K]

# INITIAL GUESSES
Mach3_init = 0.5                  # rotor/turbine exit Mach number [-]
alpha2_init  = np.radians(75)     # stator angle [deg->rad]
beta3_init = np.radians(-65)      # rotor angle [deg->rad]


RHT = 0.9                       # ratio hub/tip radius
mdot = 0.777482308              # total mass flow [kg/s]
one.rho = 0.98317               # inlet density [kg/m^3]
thr.geo.Rh = 0.106              # radius of hub at rotor exit [m]
one.alpha = 0                   # stator inlet angle (0 if asusmed axial) [rad]



## ok ok here we go



# use initial values to begin study
thr.vel.M = Mach3_init
two.alpha = alpha2_init
thr.beta = beta3_init

stator.eta = eta_stator_init
rotor.eta = eta_rotor_init


# mach 3 converge until calculated value meets the initial guess
# inside of this function, alpha2 and beta 3 are also set so that the DeltaH is met
one, two, thr, stator, rotor = converge_mach3(one, two, thr, stator, rotor, gamma, cp, R, GR, psi, DeltaH_prod, bounds_angles)


# loss coefficients (stator)
stator = f.losses(two.vel.V, two.vel.Vs, one.P0, two.P0, two.P, stator)

# loss coefficients (rotor)
rotor = f.losses(thr.vel.W, thr.vel.Ws, two.P0r, thr.P0r, thr.P, rotor)

# DeltaH_calculations
DeltaH_calc = thr.vel.U*(two.vel.Vu-thr.vel.Vu)    # euler formulation
DeltaH_T = cp*(one.T0-thr.T0)                      # total temperatur formulation


# compute total condtions
thr.T0, thr.P0 = f.total_conditions(thr.T, thr.vel.V, thr.P, cp, gamma)

# blade geometry at 2
two.geo.Rt, two.geo.Rh, two.geo.h, two.geo.Rm, two.geo.Dm, two.geo.A = f.blade_geometry(mdot, two.rho, two.vel.Vx, RHT)

# rotational velocity
two.vel.Omega = two.vel.U/two.geo.Rm
two.vel.RPM = f.RPM(two.vel.Omega)

# massflow and characteristics at inlet
one.rho0 = f.static_density(one.P0, one.T0, R)                                  # = one.P0/one.T0/R
one.rho = solve_rhox(0.9999*one.rho0, one.rho0, one.T0, mdot,  two.geo.A, gamma, cp)  # iteraciones hasta converger
one.vel.V = np.sqrt(2*cp*one.T0*(1-(one.rho/one.rho0)**(gamma-1)))
one.vel.Vx = one.vel.V
one.T = one.T0 - one.vel.V**2/2/cp
one.P = f.T2P(one.P0, one.T, one.T0, gamma)                     # = one.P0*(one.T/one.T0)**(gamma/(gamma-1)) # isentropic
one.vel.a, one.vel.M = f.sonic(gamma, R, one.T, one.vel.V)
one.geo.A = mdot/one.rho/one.vel.V

# blade geometry at 1
one.geo.Rt, one.geo.Rh, one.geo.h, one.geo.Rm, one.geo.Dm, one.geo.A = f.blade_geometry(mdot, one.rho, one.vel.V, RHT)

# blade geometry at 3
# what
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
thr.Tss = one.T0*(thr.P/one.P0)**((gamma-1)/gamma)
thr.Hs = cp*thr.Ts
stator.eta = (one.H0 - thr.H0)/(one.H0 - (thr.Hs+thr.vel.V**2/2))
rotor.eta = (one.H0 - thr.H0)/(one.H0 - thr.Hs)

epsp = (two.vel.Vs**2 - two.vel.V**2)/two.vel.V**2
epspp = (thr.vel.Ws**2 - thr.vel.W**2)/thr.vel.W**2

knls_stator = epsp*(1+gamma*two.vel.M**2/2)
knls_rotor = epspp*(1+gamma*thr.vel.Mr**2/2)

stg = np.radians(60)
opt_pitch_chord = 0.6*np.cos(stg)/2/np.cos(two.alpha)**2/(one.vel.V/two.vel.Vx*np.tan(one.alpha) + np.tan(two.alpha))*(one.geo.h/two.geo.h)*one.P0/two.P0






# PLOT TABLE OF RESULTS
#tab.convergence_parameters(one, two, thr, Mach3_init, alpha2_init, beta3_init, DeltaH_prod, DeltaH_calc)
tab.print_all_tables(one, two, thr, Mach3_init, alpha2_init, beta3_init, DeltaH_prod, DeltaH_calc, DeltaH_T, mdot, Deltabeta, psi, GR, h3h2, stator, rotor, eta_stator_init, eta_rotor_init)

# ORIGINAL GRAPH FROM JS GAS GENERATOR STUDY
# graphs.velocity_triangle(753.18, 197.59, -71.35, 254.09 ,279.64, 197.6, -544.89, 254.09)

# GRAPH VELOCITY TRIANGLE FROM RESULTS
graphs.velocity_triangle(two.vel.Vu, two.vel.Vx, thr.vel.Vu, thr.vel.Vx, two.vel.Wu, two.vel.Wx, thr.vel.Wu, thr.vel.Wx)