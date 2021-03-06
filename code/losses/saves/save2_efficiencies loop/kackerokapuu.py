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
from scipy import interpolate


def trailing_edge(M3, gamma, t_o):

    # for axial entry blades (not impulse blades)
    dataset = np.genfromtxt('dataset_trailing_edge.csv', delimiter=',')
    x = dataset[:,0]
    y = dataset[:,1]

    spline = interpolate.InterpolatedUnivariateSpline(x,y)

    dphi2 = np.float64(spline(t_o))

    pw = -gamma/(gamma-1)
    YTET = ( (1-(gamma-1)/2*M3**2*(1/(1-dphi2)-1))**pw -1)/(1-(1+(gamma-1)/2*M3**2)**pw)

    return YTET




def tip_clearance():
    return 0





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




def kackerokapuu(component, S, alpha2, alpha3, t, c, bx, h, M2, M3, P2, P3, gamma, RHT, Rec, t_o):

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


    YP = profile_losses(S/c, alpha3, t/c, alpha2, M2, M3, K, P2, P3, gamma, RHT)
    fRe = reynolds_correction(Rec)
    YS = secondary_losses(alpha2, alpha3, h/c, bx, h, M2, M3)
    YTC = tip_clearance()
    YTET = trailing_edge(M3, gamma, t_o)

    Ytot = YP*fRe + YS + YTET + YTC

    return Ytot













