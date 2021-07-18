# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 00:51:22 2021

@author: alondra sola

Project: turbine design tool (TFG)
Program: parametric analysis of efficiency
File: runloss.py

"""

## IMPORT NECESSARY LIBRARIES
import numpy as np                   # library for math and calculation functions

## IMPORT FUNCTIONS AND DEFINITIONS FROM OTHER PYTHON FILES
import plot_functions as graphs
import tables_results as tab
from definitions import plane, component
from converge_efficiencies import converge_efficiencies
from no_losses import no_losses



def runloss(psi):



    # initialize class variables
    one = plane()
    two = plane()
    thr = plane()
    rotor = component()
    stator = component()




    ###################### GIVEN PARAMETERS AND ASSUMTIOSN ########################

    # turbine model from gas generator for automobile study
    # jorge saavedra

    # GIVEN PARAMETERS (TURBINE)
    one.T0 = 1400                     # inlet total temperature (exit combustor) = T4t [K]
    one.P0 = 397194                  # inlet total pressure (exit combustor) [Pa]
    thr.P0 = 101372                    # outlet static pressure [Pa]
    thr.T0 = 1085.12                  #  # outlet total temperature [K]
    DeltaH_prod = 390000           #  # enthalpy produced by turbine [J/Kg]

    # # DESIGN VARIABLES
    mdot = 0.777                        # total mass flow [kg/s]
    one.alpha = 0                     # stator inlet angle (0 because flow is axial) [rad]
    GR = 0.32                         # reaction degree [-]
    # psi = 1.75                         # loading factor [-]
    RHT = 0.9                      #  # ratio hub/tip radius


    # # FLUID PROPERTIES (TURBINE)
    gamma = 1.3              # [-]
    cp = 1240                # [J/kg/K]
    R = 286.1538462          # [J/kg/K]



    # INITIAL GUESSES
    Mach3_init = 0.5                 #  # rotor/turbine exit Mach number [-]
    alpha2_init  = np.radians(72.5)     # stator angle [deg->rad]
    beta3_init = np.radians(-69.5)   #  # rotor angle [deg->rad]
    eta_stator_init = 0.8888           #  # stator efficiency [-]
    eta_rotor_init = 0.808            #  # rotor efficiency [-]
    # upper and lower bounds for alpha2 and beta3
    # in this format: [alpha2_min, beta3_min], [alpha2_max, beta3_max]
    bounds_angles = ([np.radians(70),np.radians(-70)], [np.radians(75), np.radians(-65)])


    # # ASSUMPTIONS
    h_c_stator = 0.7                  # height/chord ratio stator [-]
    h_c_rotor = 1.4                     # height/chord ratio rotor [-]
    t_o = 0.22                        # trailing-egde thickness to throat opening ratio [-]


    ## ok ok here we goet


    # use initial values to begin study
    thr.vel.M = Mach3_init
    two.alpha = alpha2_init
    thr.beta = beta3_init
    two.geo.hc = h_c_stator
    thr.geo.hc = h_c_rotor
    two.geo.to = t_o
    thr.geo.to = t_o




    etas = np.array([eta_stator_init, eta_rotor_init])
    converge_efficiencies(etas, stator, rotor, one, two, thr, gamma, cp, R, GR, psi, DeltaH_prod, bounds_angles, RHT, mdot)


    # print("Stator pressure loss: ", round(stator.omegaKO,2))
    # print("Rotor pressure loss: ", round(rotor.omegaKO,2))
    # print("Stator efficiency: ", round(stator.eta,2))
    # print("Rotor efficiency: ", round(rotor.eta,2))
    # print("Total stage efficiency: ", round(stator.eta*rotor.eta,4),'\n')

    results = [stator.eta, rotor.eta, stator.omegaKO, rotor.omegaKO]

    return results