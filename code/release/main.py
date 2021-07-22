# -*- coding: utf-8 -*-
"""
Created on Fri May 14 11:11:08 2021

@author: alondra sola

Project: turbine design tool (TFG)
Program: release of turbine optimisation tool with loss mechanisms
File: main.py

"""

## IMPORT NECESSARY LIBRARIES
import numpy as np                   # library for math and calculation functions

## IMPORT FUNCTIONS AND DEFINITIONS FROM OTHER PYTHON FILES
import plot_functions as graphs
import tables_results as tab
from definitions import plane, component
from converge_efficiencies import converge_efficiencies

def main():

    # initialize class variables
    one = plane()
    two = plane()
    thr = plane()
    rotor = component()
    stator = component()


    ###################### GIVEN PARAMETERS AND ASSUMTIOSN ########################

    # turbine model from Design Study for Single Stage High Pressure Turbine of Gas
    # Turbine Engines
    # Ajoko, Tolumoye John
    # LYULKA AL-2LF-3


    # GIVEN PARAMETERS (TURBINE)
    one.T0 = 562                     # inlet total temperature (exit combustor) = T4t [K]
    one.P0 = 346325                  # inlet total pressure (exit combustor) [Pa]
    thr.P0 = 105305                    # outlet static pressure [Pa]
    thr.T0 = 562-118                  #  # outlet total temperature [K]
    DeltaH_prod = 180000 #3.84e+06/25.5            #  # enthalpy produced by turbine [J/Kg]

    # # DESIGN VARIABLES
    mdot =  25.5                        # total mass flow [kg/s]
    one.alpha = 0                     # stator inlet angle (0 because flow is axial) [rad]
    GR = 0.4                      # reaction degree [-]
    psi = 2.065                         # loading factor [-]
    RHT = 0.98                      #  # ratio hub/tip radius

    # FLUID PROPERTIES (TURBINE)
    gamma = 1.29              # [-]
    cp = 1277                # [J/kg/K]
    R = 286.1538462          # [J/kg/K]

    # INITIAL GUESSES
    Mach3_init = 0.46                 #  # rotor/turbine exit Mach number [-]
    alpha2_init  = np.radians(63.5)     # stator angle [deg->rad]
    beta3_init = np.radians(-50)   #  # rotor angle [deg->rad]
    tau = 1710/1757
    etat = (1-tau)/(1-tau**(1/0.9))
    eta_stator_init = np.sqrt(etat)           #  # stator efficiency [-]
    eta_rotor_init = np.sqrt(etat)            #  # rotor efficiency [-]
    # upper and lower bounds for alpha2 and beta3
    # in this format: [alpha2_min, beta3_min], [alpha2_max, beta3_max]
    bounds_angles = ([np.radians(60),np.radians(-70)], [np.radians(75), np.radians(-50)])


    # # ASSUMPTIONS
    h_c_stator = 0.75                  # height/chord ratio stator [-]
    h_c_rotor = 0.75                     # height/chord ratio rotor [-]
    t_o = 0.22                        # trailing-egde thickness to throat opening ratio [-]

    # expected results
    expected = [etat, 0.503, 1.18, 0.25, 0.04]



    ## ok ok here we go

    # use initial values to begin study
    thr.vel.M = Mach3_init
    two.alpha = alpha2_init
    thr.beta = beta3_init
    two.geo.hc = h_c_stator
    thr.geo.hc = h_c_rotor
    two.geo.to = t_o
    thr.geo.to = t_o


    # find efficiency
    etas = np.array([eta_stator_init, eta_rotor_init])
    converge_efficiencies(etas, stator, rotor, one, two, thr, gamma, cp, R, GR, psi, DeltaH_prod, bounds_angles, RHT, mdot)

    inputs = [DeltaH_prod, mdot, GR, psi, RHT]
    tab.print_testcase_terminal(one, two, thr, inputs, expected, round(stator.eta*rotor.eta,4), stator, rotor)


    # %% GRAPHS

    # ORIGINAL GRAPH FROM JS GAS GENERATOR STUDY
    # graphs.velocity_triangle(753.18, 197.59, -71.35, 254.09 ,279.64, 197.6, -544.89, 254.09)

    # GRAPH VELOCITY TRIANGLE FROM RESULTS
    graphs.velocity_triangle(two.vel.Vu, two.vel.Vx, thr.vel.Vu, thr.vel.Vx, two.vel.Wu, two.vel.Wx, thr.vel.Wu, thr.vel.Wx)





if __name__ == '__main__':
    main()
