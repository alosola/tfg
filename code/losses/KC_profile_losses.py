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


def profile_losses_alpha2_0(S_C, alpha3):
    """
    for figure 1: when alpha2 = 0
    NOTE: in Kacker-Okapuu paper, the entry angle is beta1 instead of alpha2
    """

    alpha3 = np.degrees(alpha3) # NOTE: this entire function is defined for alpha3 in DEGREES, not RADIANS

    x = np.genfromtxt('fig1_x.csv', delimiter=',')               # import SC vector
    y = np.flip(np.genfromtxt('fig1_y.csv', delimiter=','))      # import alpha2 vector
    Z = np.flip(np.genfromtxt('fig1_Zm.csv', delimiter=','), 1)  # import YP mesh

    YP_spline =  scipy.interpolate.RectBivariateSpline(x, y, Z)  # create spline for evaluation
    YP = YP_spline.ev(S_C, alpha3)                                # evaluate at desired point

    # print warking messages if values are outside recommended area
    if ( S_C < np.amin(x) ) or ( S_C > np.amax(x) ) :
        print("Warning: S/C value is outside the data range for profile pressure loss calculations. Program will contine with data extrapolation ")

    if ( alpha3 < np.amin(y) ) or ( alpha3 > np.amax(y) ) :
        print("Warning: alpha2 value is outside the data range for profile pressure loss calculations. Program will contine with data extrapolation ")

    return YP




def profile_losses_alpha2_alpha3(S_C, alpha3):
    """
    for figure 2: when alpha2 = alpha3
    NOTE: in Kacker-Okapuu paper, the entry angle is beta1 instead of alpha2,
    and the outlet angle is alpha2 instead of alpha3
    """

    alpha3 = np.degrees(alpha3) # NOTE: this entire function is defined for alpha2 in DEGREES, not RADIANS

    x = np.genfromtxt('fig2_x.csv', delimiter=',')               # import SC vector
    y = np.genfromtxt('fig2_y.csv', delimiter=',')      # import alpha2 vector
    Z = np.flip(np.genfromtxt('fig2_Zm.csv', delimiter=','), 1)  # import YP mesh

    YP_spline =  scipy.interpolate.RectBivariateSpline(x, y, Z)  # create spline for evaluation
    YP = YP_spline.ev(S_C, alpha3)                                # evaluate at desired point

    # print warking messages if values are outside recommended area
    if ( S_C < np.amin(x) ) or ( S_C > np.amax(x) ) :
        print("Warning: S/C value is outside the data range for profile pressure loss calculations. Program will contine with data extrapolation ")

    if ( alpha3 < np.amin(y) ) or ( alpha3 > np.amax(y) ) :
        print("Warning: alpha2 value is outside the data range for profile pressure loss calculations. Program will contine with data extrapolation ")

    return YP





def thickness_chord_ratio(deltabeta):

    deltabeta = abs(np.degrees(deltabeta))

    if deltabeta > 120:
        t_c = 0.25
    if deltabeta < 40:
        t_c = 0.15
    else:
        t_c = 0.15 + 1.25e-03*(deltabeta - 40)

    return t_c



def profile_losses_AMDC(S_C, alpha3, alpha2):

    t_c = thickness_chord_ratio(alpha2 - alpha3)

    # AMDC profile loss correlation

    YP_alpha2_0 = profile_losses_alpha2_0(S_C, alpha3)
    YP_alpha2_alpha3 = profile_losses_alpha2_alpha3(S_C, alpha3)

    B_A = alpha2/alpha3; # in reference, beta1/alpha2

    YP_AMDC = (YP_alpha2_0 + np.absolute(B_A)*B_A*(YP_alpha2_alpha3 - YP_alpha2_0))*(t_c/0.2)**B_A

    return YP_AMDC



def profile_losses_YShock(M2, M3, K, RHT, gamma, P2, P3):

    M2H = M2*(1+K*np.absolute(RHT-1)**2.2)
    dp_hub = 0.5 #0.75*(M2H - 0.4)**1.75
    dp_shock = dp_hub*RHT

    def f(M, gamma):
        return 1- (1+(gamma-1)/2*M**2)**(gamma/(gamma-1))

    YShock = dp_shock*(P2/P3)*f(M2,gamma)/f(M3,gamma)

    return YShock/1000


def losses_KP(M2, M3):
    # M2 is Mach in
    # M3 is Mach out
    # (does not match nomenclature in Kacker-Okapuu)

    KP = 1 - 1.25*(M3 - 0.2)*(M2- M3)**2

    return KP


def profile_losses(S_C, alpha3, alpha2, M2, M3, K, P2, P3, gamma, RHT):

    YP_AMDC = profile_losses_AMDC(S_C, alpha3, alpha2)
    KP = losses_KP(M2, M3)
    YShock = profile_losses_YShock(M2, M3, K, RHT, gamma, P2, P3)
    YP = 0.914*(2/3*YP_AMDC*KP + YShock)

    return YP
