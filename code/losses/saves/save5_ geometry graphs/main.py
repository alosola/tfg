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
import matplotlib.pyplot as plt

## IMPORT FUNCTIONS AND DEFINITIONS FROM OTHER PYTHON FILES
import aux_functions as f
import plot_functions as graphs
import tables_results as tab
from definitions import plane, component
from converge_mach3 import converge_mach3
from solve_functions import solve_geometry
from check_limits import check_limits
from converge_efficiencies import converge_efficiencies


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


mdot = 0.777482308             # total mass flow [kg/s]
one.rho = 0.98317                 # inlet density [kg/m^3]


# FLUID PROPERTIES (TURBINE)
gamma = 1.3              # [-]
cp = 1240                # [J/kg/K]
R = 286.1538462          # [J/kg/K]


# INITIAL GUESSES
Mach3_init = 0.5                  # rotor/turbine exit Mach number [-]
alpha2_init  = np.radians(75)     # stator angle [deg->rad]
beta3_init = np.radians(-65)      # rotor angle [deg->rad]
eta_stator_init = 0.919              # stator efficiency [-]
eta_rotor_init = 0.827               # rotor efficiency [-]
# upper and lower bounds for alpha2 and beta3
# in this format: [alpha2_min, beta3_min], [alpha2_max, beta3_max]
bounds_angles = ([np.radians(70),np.radians(-70)], [np.radians(75), np.radians(-60)])



one.alpha = 0                   # stator inlet angle (0 because flow is axial) [rad]



# ASSUMPTIONS
GR = 0.4                          # reaction degree [-]
psi = 1.75                        # loading factor [-]
RHT = 0.7                         # ratio hub/tip radius
h_c_stator = 0.7                  # height/chord ratio stator [-]
h_c_rotor = 1.4                   # height/chord ratio rotor [-]
t_o = 0.25                        # trailing-egde thickness to throat opening ratio [-]

## ok ok here we go

# use initial values to begin study
thr.vel.M = Mach3_init
two.alpha = alpha2_init
thr.beta = beta3_init
two.geo.hc = h_c_stator
thr.geo.hc = h_c_rotor
two.geo.to = t_o
thr.geo.to = t_o

# using the initial efficiencies, calculate the cycle
etas = np.array([eta_stator_init, eta_rotor_init])
converge_efficiencies(etas, stator, rotor, one, two, thr, gamma, cp, R, GR, psi, DeltaH_prod, bounds_angles, RHT, mdot)



# 1.3
# DeltaH_calculations
DeltaH_calc = thr.vel.U*(two.vel.Vu-thr.vel.Vu)    # euler formulation
DeltaH_T = cp*(one.T0-thr.T0)                      # total temperature formulation





epsp = (two.vel.Vs**2 - two.vel.V**2)/two.vel.V**2
epspp = (thr.vel.Ws**2 - thr.vel.W**2)/thr.vel.W**2

knls_stator = epsp*(1+gamma*two.vel.M**2/2)
knls_rotor = epspp*(1+gamma*thr.vel.Mr**2/2)

stg = np.radians(60)
opt_pitch_chord = 0.6*np.cos(stg)/2/np.cos(two.alpha)**2/(one.vel.V/two.vel.Vx*np.tan(one.alpha) + np.tan(two.alpha))*(one.geo.h/two.geo.h)*one.P0/two.P0



# calculate pressure loss coefficient YP = tpl using Kacker-Okapuu
# rotor.tplKC = kackerokapuu(component, S, alpha2, alpha3, t, c, bx, h, two.vel.M, thr.vel.M, two.P, thr.P, gamma, RHT, Rec, t_o)



# PLOT TABLE OF RESULTS
#tab.convergence_parameters(one, two, thr, Mach3_init, alpha2_init, beta3_init, DeltaH_prod, DeltaH_calc)

# extra calculations fro print results (rearrange this later)
GR_enthalpy = (two.T - thr.T)/(two.T0 - thr.T0) # Assume v1 = v3, not right
Deltabeta = two.beta - thr.beta # por qué no negativo?
Deltaalpha = two.alpha - thr.alpha
GR_final = (two.T - thr.T)/(one.T - thr.T)
DeltaW = thr.vel.W - two.vel.W
h3h2 = thr.geo.h/two.geo.h

# %% Tables and results
tab.print_all_tables(one, two, thr, Mach3_init, alpha2_init, beta3_init, DeltaH_prod, DeltaH_calc, DeltaH_T, mdot, Deltabeta, psi, GR, h3h2, stator, rotor, eta_stator_init, eta_rotor_init)

print("Stator pressure loss: ", round(stator.omegaKC,2))
print("Rotor pressure loss: ", round(rotor.omegaKC,2))

print("Stator efficiency: ", round(stator.eta,2))
print("Rotor efficiency: ", round(rotor.eta,2))

# %% GRAPHS

# ORIGINAL GRAPH FROM JS GAS GENERATOR STUDY
# graphs.velocity_triangle(753.18, 197.59, -71.35, 254.09 ,279.64, 197.6, -544.89, 254.09)

# GRAPH VELOCITY TRIANGLE FROM RESULTS
graphs.velocity_triangle(two.vel.Vu, two.vel.Vx, thr.vel.Vu, thr.vel.Vx, two.vel.Wu, two.vel.Wx, thr.vel.Wu, thr.vel.Wx)


# %% GRAPH ROTOR GEOMETRY
graphs.geometry(one, two, thr)

# fig2 = plt.figure()
# ax2 = fig2.add_subplot(111)

# Xstator = 0
# Wstator = two.geo.c
# Xrotor = Xstator + Wstator + Wstator/5
# Wrotor = thr.geo.c


# stator = [[Xstator, one.geo.Rh], [Xstator, one.geo.Rt], [Xstator+Wstator, two.geo.Rt], [Xstator+Wstator, two.geo.Rh], [Xstator,one.geo.Rh]] #the points to trace the edges.
# polygon_stator = plt.Polygon(stator,  fill=True, edgecolor=None, facecolor=(151/222, 178/222, 1))
# ax2.add_patch(polygon_stator)


# rotor = [[Xrotor, two.geo.Rh], [Xrotor, two.geo.Rt], [Xrotor+Wrotor, thr.geo.Rt], [Xrotor+Wrotor, thr.geo.Rh], [Xrotor,two.geo.Rh]] #the points to trace the edges.
# polygon_rotor = plt.Polygon(rotor,  fill=True, edgecolor=None, facecolor=(173/217, 1, 165/217))
# ax2.add_patch(polygon_rotor)

# antistator = [[Xstator, -one.geo.Rt], [Xstator, -one.geo.Rh], [Xstator+Wstator, -two.geo.Rh], [Xstator+Wstator, -two.geo.Rt], [Xstator,-one.geo.Rt]] #the points to trace the edges.
# polygon_antistator = plt.Polygon(antistator,  fill=True, edgecolor=None, facecolor=(151/222, 178/222, 1))
# ax2.add_patch(polygon_antistator)

# antirotor = [[Xrotor, -two.geo.Rt], [Xrotor, -two.geo.Rh], [Xrotor+Wrotor, -thr.geo.Rh], [Xrotor+Wrotor, -thr.geo.Rt], [Xrotor,-two.geo.Rt]] #the points to trace the edges.
# polygon_antirotor = plt.Polygon(antirotor,  fill=True, edgecolor=None, facecolor=(173/217, 1, 165/217))
# ax2.add_patch(polygon_antirotor)


# plt.axis([-0.3*Xrotor,(Xrotor + Wrotor + 0.3*Xrotor),1.3*thr.geo.Rt,-1.3*thr.geo.Rt])


# plt.show()




