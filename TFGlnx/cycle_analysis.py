#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 19:23:36 2021

@author: alondra sola

Engine cycle analysis following the design analysis example (Saavedra, 2014)


"""

def main():

    # IMPORT NECESSARY LIBRARIES
    # import numpy as np



    # DEFINE DESIGN CONSTRAINTS

    pi_compressor = 4       # compressor pressure ratio [-]
    P01 = 101325            # inlet total pressure [Pa]
    T01 = 288               # inlet total temperature [K]
    NPw = 150000            # net power to wheels [W]
    T03 = 1400              # exit combustor temperature [K]


    # FIX FLUID PROPERTIES

    gamma_c = 1.4  # [-]
    cp_c = 1000    # [J/kg/K]
    gamma_t = 1.3  # [-]
    cp_t = 1240    # [J/kg/K]

    # COMPUTE REMAINING NECESSARY FLOW QUANTITIES
    cycle_analysis(P01, T01, T03, pi_compressor, gamma_c, cp_c, gamma_t, cp_t, NPw, True)

def cycle_analysis(P01, T01, T03, pi_compressor, gamma_c, cp_c, gamma_t, cp_t, NPw, print_results):

    # ASSUMPTIONS:
    eta_compressor = 0.82  # compressor efficiency [-]
    eta_shaft = 0.94       # shaft efficiency [-]
    eta_turbine = 0.836     # turbine efficiency [-]
    eta_burner = 0.88      # burner efficiency [-]
    eta_combustor_p = 0.98 # compressor pressure efficiency [-]
    eta_trans = 0.89       # transmission efficiency [-]

    # DEFINE VALUES
    L = 44400000             # gasoline energy density [J/kg]


    ####################
    # COMPRESSOR STAGE #

    # the compression ratio provides the pressure after the compressor stage, ahead of the combustion chamber
    P02 = P01*pi_compressor

    # assuming an isentropic evolution, the temperature is easily found
    T02s = T01*pi_compressor**((gamma_c-1)/gamma_c)

    # the calculation of the total temperature takes into account the compressor efficiency
    T02 = T01 + (T02s - T01)/eta_compressor

    # the work required is founfrom the temperature increase
    DeltaH_compressor = cp_c*(T02 - T01)


    ###################
    # COMBUSTOR STAGE #

    # the pressure at the outlet can be assumed with the efficiency
    P03 = P02*eta_combustor_p

    # the exit temperature is a design constraint, and the temperature increase is calculated
#    DeltaT_combustor = T03 - T02

    # using the temperatures and energy density, the fuel-to-mass-flow ratio can be estimated
    f = (cp_c*T03 - cp_c*T02)/(eta_burner*L - cp_c*T03)


    ############################
    # TURBINE EXPLANSION STAGE #

    # the flow is expanded to atmospheric pressure (total)
    P04 = P01

    # the expansion ratio can be calculated immediately
    pi_turbine = P03/P04

    # assuming no losses. the expansion is isentropic
    T04s = T03*(1/pi_turbine)**((gamma_t - 1)/gamma_t)

    # incorporating turbine efficiency
    T04 = T03 - eta_turbine*(T03-T04s)


    ##################
    # WORK AND POWER #

    # the available work is
    DeltaH_turbine = cp_t*(T03-T04)

    # the work for transmission is the work not needed for compressor
    DeltaH_av_trans = DeltaH_turbine - DeltaH_compressor/eta_shaft/(1+f)

    # since the net power needed by the wheels is a design constraint, the massflow is calculated
    ma = NPw/(DeltaH_av_trans*eta_trans*(1+f))

    if print_results==True:
        print('Massflow: ', round(ma,4) ,' kg/s')
        print('T04: ', round(T04,2) , ' K')
        print('P03: ', round(P03,0) , ' Pa')
        print('Expansion ratio: ', round(pi_turbine,2) )
        print('Delta H compressor: ', round(DeltaH_compressor/1000, 2), ' kJ/kg')
        print('Delta H available for transmission: ', round(DeltaH_av_trans/1000, 2), ' kJ/kg')


    return ma, pi_turbine, DeltaH_turbine, P03, P04, T03, T04;





if __name__ == '__main__':
    main()
