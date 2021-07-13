# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 19:39:36 2021

@author: alondra sola

Project: turbine design tool (TFG)
Program: inclusion of losses for 1D turbine
File: KO_profile_losses.py
Descpription: functions for calculation of profile losses within Kacker-Okapuu

"""

## IMPORT NECESSARY LIBRARIES
import numpy as np                   # library for math and calculation functions
import scipy.interpolate


def profile_losses_alpha_in_0(S_C, alpha_out):
    """
    for figure 1: when alpha_in = 0
    NOTE: in Kacker-Okapuu paper, the entry angle is beta1 instead of alpha_in
    """

    alpha_out = np.degrees(alpha_out) # NOTE: this entire function is defined for alpha_out in DEGREES, not RADIANS

    x = np.genfromtxt('fig1_x.csv', delimiter=',')               # import SC vector
    y = np.flip(np.genfromtxt('fig1_y.csv', delimiter=','))      # import alpha_in vector
    Z = np.flip(np.genfromtxt('fig1_Zm.csv', delimiter=','), 1)  # import YP mesh

    YP_spline =  scipy.interpolate.RectBivariateSpline(x, y, Z)  # create spline for evaluation
    YP = YP_spline.ev(S_C, alpha_out)                               # evaluate at desired point

    return YP




def profile_losses_alpha_in_alpha_out(S_C, alpha_out):
    """
    for figure 2: when alpha_in = alpha_out
    NOTE: in Kacker-Okapuu paper, the entry angle is beta1 instead of alpha_in,
    and the outlet angle is alpha_in instead of alpha_out
    """

    alpha_out = np.degrees(alpha_out) # NOTE: this entire function is defined for alpha_in in DEGREES, not RADIANS

    x = np.genfromtxt('fig2_x.csv', delimiter=',')               # import SC vector
    y = np.genfromtxt('fig2_y.csv', delimiter=',')      # import alpha_in vector
    Z = np.flip(np.genfromtxt('fig2_Zm.csv', delimiter=','), 1)  # import YP mesh

    YP_spline =  scipy.interpolate.RectBivariateSpline(x, y, Z)  # create spline for evaluation
    YP = YP_spline.ev(S_C, alpha_out)                               # evaluate at desired point

    return YP





def thickness_chord_ratio(sumtheta):

    sumtheta = abs(np.degrees(sumtheta))

    if sumtheta > 120:
        t_c = 0.25
    elif sumtheta < 40:
        t_c = 0.15
    else:
        t_c = 0.15 + 1.25e-03*(sumtheta - 40)

    return t_c



def profile_losses_AMDC(S_C, alpha_out, alpha_in):

    t_c = thickness_chord_ratio(alpha_in - alpha_out)

    # AMDC profile loss correlation

    YP_alpha_in_0 = profile_losses_alpha_in_0(S_C, alpha_out)
    YP_alpha_in_alpha_out = profile_losses_alpha_in_alpha_out(S_C, alpha_out)

    B_A = alpha_in/alpha_out; # in reference, beta1/alpha_in

    YP_AMDC = (YP_alpha_in_0 + np.absolute(B_A)*B_A*(YP_alpha_in_alpha_out - YP_alpha_in_0))*(t_c/0.2)**B_A

    return YP_AMDC



def profile_losses_YShock(M_in, M_out, K, RHT, gamma, P_in, P_out):

    M_inH = M_in*(1+K*np.absolute(RHT-1)**2.2)

    if M_inH < 0.4:
        dp_hub = 0
    else:
        dp_hub = 0.75*(M_inH - 0.4)**1.75
    dp_shock = dp_hub*RHT

    def f(M, gamma):
        return 1- (1+(gamma-1)/2*M**2)**(gamma/(gamma-1))

    YShock = dp_shock*(P_in/P_out)*f(M_in,gamma)/f(M_out,gamma)

    return YShock/1000


def losses_KP(MA, MB):

    K2 = (MA/MB)**2

    if MB < 0.2:
        K1 = 1
    else:
        K1 = 1 - 1.25*(MB-0.2)

    KP = 1-K2*(1-K1)

    return KP


def profile_losses(S_C, alpha_out, alpha_in, M_in, M_out, K, P_in, P_out, gamma, RHT):

    YP_AMDC = profile_losses_AMDC(S_C, alpha_out, alpha_in)
    KP = losses_KP(M_in, M_out)
    YShock = profile_losses_YShock(M_in, M_out, K, RHT, gamma, P_in, P_out)
    YP = 0.914*(2/3*YP_AMDC*KP + YShock)

    return YP
