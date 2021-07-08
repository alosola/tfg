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

## import the secondary fucntions defined in other files
from KC_profile_losses import profile_losses
from KC_secondary_losses import secondary_losses
import scipy.interpolate


def trailing_edge(M3, gamma, t_o):

    # for axial entry blades (not impulse blades)
    dataset = np.genfromtxt('dataset_trailing_edge.csv', delimiter=',')
    x = dataset[:,0]
    y = dataset[:,1]

    spline = scipy.interpolate.InterpolatedUnivariateSpline(x,y)

    dphi2 = np.float64(spline(t_o))

    pw = -gamma/(gamma-1)
    YTET = ( (1-(gamma-1)/2*M3**2*(1/(1-dphi2)-1))**pw -1)/(1-(1+(gamma-1)/2*M3**2)**pw)

    return YTET




def tip_clearance():
    losses = 0
    return losses





def reynolds_correction(Rec):
    # Rec = Reynolds number based on true chord and exit gas conditions

    mean = 2*10**5

    if Rec <= mean:
        fRe = (Rec/mean)**-0.4
    elif Rec < 10**6:
        fRe = 1
    else:
        fRe = (Rec/10**6)**-0.2

    return fRe





def kackerokapuu(component, s, alpha2, alpha3, c, bx, h, M2, M3, P2, P3, gamma, RHT, Re, t_o):

    """
    Function for calculation of total losses using Kacker-Okapuu formulation
    Y.P = profile losses
    Y.Re = Reynolds correction parameter
    Y.S = secondary losses
    Y.TET = trailing edge losses
    Y.TC = tip leakage losses

    Total losses = Ytot = Y.P * Y.Re + Y.S + Y.TET + Y.TC
    """

    if component == 'stator' :
        K = 1.8
    elif component == 'rotor' :
        K = 5.2
    else :
        print('Error: a turbine component must be specified for K estimation in Kacker-Okapuu losses.')
        K = 5.2


    YP = profile_losses(s/c, alpha3, alpha2, M2, M3, K, P2, P3, gamma, RHT)
    fRe = reynolds_correction(Re)
    YS = secondary_losses(alpha2, alpha3, h/c, bx, h, M2, M3)
    YTC = tip_clearance()
    YTET = trailing_edge(M3, gamma, t_o)

    Ytot = YP*fRe + YS + YTET + YTC

    return Ytot



# S = 0.008549314    # pitch [m]
# alpha2 = np.radians(75.3)
# alpha3 = np.radians(65)
# c = 0.017007665 # chord
# h = 0.011906101
# bx = 0.6*h
# M2 = 1.11
# M3 = 0.46
# P2 = 153781.890
# P3 = 90483.172
# gamma = 1.3
# RHT  = 0.9
# Re = 42838
# t_o = 0.15

# loss_rotor = kackerokapuu('rotor', S, alpha2, alpha3, c, bx, h, M2, M3, P2, P3, gamma, RHT, Re, t_o)





