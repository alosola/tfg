# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 19:39:36 2021

@author: alondra sola

Project: turbine design tool (TFG)
Program: inclusion of losses for 1D turbine
File: KC_profile_losses.py
Descpription: functions for calculation of profile losses within Kacker-Okapuu

"""

## IMPORT NECESSARY LIBRARIES
import numpy as np                   # library for math and calculation functions
import scipy.interpolate


def profile_losses_beta1_0(SC, alpha2):
    # for figure 2: when beta1 = 0

    alpha2 = np.degrees(alpha2) # NOTE: this entire function is defined for alpha2 in DEGREES, not RADIANS

    x = np.genfromtxt('fig1_x.csv', delimiter=',')               # import SC vector
    y = np.flip(np.genfromtxt('fig1_y.csv', delimiter=','))      # import alpha2 vector
    Z = np.flip(np.genfromtxt('fig1_Zm.csv', delimiter=','), 1)  # import YP mesh

    YP_spline =  scipy.interpolate.RectBivariateSpline(x, y, Z)  # create spline for evaluation
    YP = YP_spline.ev(SC, alpha2)                                # evaluate at desired point

    # print warking messages if values are outside recommended area
    if ( SC < np.amin(x) ) or ( SC > np.amax(x) ) :
        print("Warning: S/C value is outside the data range for profile pressure loss calculations. Program will contine with data extrapolation ")

    if ( alpha2 < np.amin(y) ) or ( alpha2 > np.amax(y) ) :
        print("Warning: alpha2 value is outside the data range for profile pressure loss calculations. Program will contine with data extrapolation ")

    return YP




def profile_losses_beta1_alpha2(SC, alpha2):
    # for figure 2: when beta1 = alpha2

    alpha2 = np.degrees(alpha2) # NOTE: this entire function is defined for alpha2 in DEGREES, not RADIANS

    x = np.genfromtxt('fig2_x.csv', delimiter=',')               # import SC vector
    y = np.genfromtxt('fig2_y.csv', delimiter=',')      # import alpha2 vector
    Z = np.flip(np.genfromtxt('fig2_Zm.csv', delimiter=','), 1)  # import YP mesh

    YP_spline =  scipy.interpolate.RectBivariateSpline(x, y, Z)  # create spline for evaluation
    YP = YP_spline.ev(SC, alpha2)                                # evaluate at desired point

    # print warking messages if values are outside recommended area
    if ( SC < np.amin(x) ) or ( SC > np.amax(x) ) :
        print("Warning: S/C value is outside the data range for profile pressure loss calculations. Program will contine with data extrapolation ")

    if ( alpha2 < np.amin(y) ) or ( alpha2 > np.amax(y) ) :
        print("Warning: alpha2 value is outside the data range for profile pressure loss calculations. Program will contine with data extrapolation ")

    return YP


def profile_losses_AMDC(S_C, alpha2, t_c, beta1):

    # AMDC profile loss correlation

    YP_beta1_0 = profile_losses_beta1_0(S_C, alpha2)
    YP_beta1_alpha2 = profile_losses_beta1_alpha2(S_C, alpha2)

    BA = beta1/alpha2;

    YP_AMDC = (YP_beta1_0 + np.absolute(BA)*BA*(YP_beta1_alpha2 - YP_beta1_0))*(t_c/0.2)**BA

    return YP_AMDC


def profile_losses_YShock(M2, M3, K, RHT, gamma, P2, P3):

    M2H = M2*(1+K*np.absolute(RHT-1)**2.2)

    YShock = 0.75*(M2H-0.4)**1.75*RHT*(P2/P3)*(1-(1+(gamma-1)*M2**2/2)**(gamma/(gamma-1)))/((1-(1+(gamma-1)*M3**2/2)**(gamma/(gamma-1))))

    return YShock


def losses_KP(M2, M3):
    # M2 is Mach in
    # M3 is Mach out
    # (does not match nomenclature in Kacker-Okapuu)

    KP = 1 - 1.25(M3 - 0.2)*(M2- M3)**2

    return KP


def profile_losses(S_C, alpha2, t_c, beta1, M2, M3, K, P2, P3, gamma, RHT):

    YP_AMDC = profile_losses_AMDC(S_C, alpha2, t_c, beta1)
    KP = losses_KP(M2, M3)
    YShock = profile_losses_YShock(M2, M3, K, RHT, gamma, P2, P3)
    YP = 0.914*(2/3*YP_AMDC*KP + YShock)

    return YP
