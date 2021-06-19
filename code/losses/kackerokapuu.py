# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 12:22:02 2021

@author: alondra sola

Project: turbine design tool (TFG)
Program: inclusion of losses for 1D turbine
File: kackerokapuu.py
Descpription: functions for calculation of pressure loss coefficient by means
of the Ainley/Mathieson/Dunham/Crane/Kacker/Okapuu method
Source: A Mean Line Prediction Method for Axial Flow Turbine Efficiency

"""

## IMPORT NECESSARY LIBRARIES
import numpy as np                   # library for math and calculation functions
import matplotlib.pyplot as plt      # library for plots
import scipy.interpolate



# def kackerokapuu():

#     Y = pressure_losses()

#     Y = total_pressure_loss_coefficient()

#     return Y.tot


def profile_losses_beta1_0(SC, alpha2):

    # alpha2 = np.degrees(alpha2)

    ## NOTE: this entire function is defined for alpha2 in DEGREES, not RADIANS

    x = np.genfromtxt('fig1_x.csv', delimiter=',')               # import SC vector
    y = np.flip(np.genfromtxt('fig1_y.csv', delimiter=','))      # import alpha2 vector
    Z = np.flip(np.genfromtxt('fig1_Zm.csv', delimiter=','), 1)  # import YP mesh

    YP_spline =  scipy.interpolate.RectBivariateSpline(x, y, Z)  # create spline for evaluation
    YP = YP_spline.ev(SC, alpha2)                                # evaluate at desired point

    # print warking messages if values are outside recommended area
    if ( SC < np.minimum(x) ) or ( SC > np.maximum(x) ) :
        print("Warning: S/C value is outside the data range for profile pressure loss calculations. Program will contine with data extrapolation ")

    if ( alpha2 < np.minimum(y) ) or ( alpha2 > np.maximum(y) ) :
        print("Warning: alpha2 value is outside the data range for profile pressure loss calculations. Program will contine with data extrapolation ")

    return YP



profile_losses_beta1_0(0.5,50)