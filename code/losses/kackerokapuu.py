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



def kackerokapuu(component, S, alpha2, alpha3, t, c, bx, h, beta1, M2, M3, K, P2, P3, gamma, RHT, Rec):

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


    # Y = pressure_losses()

    YP = profile_losses(S/c, alpha2, t/c, beta1, M2, M3, K, P2, P3, gamma, RHT)
    fRe = reynolds_correction(Rec)
    YS = secondary_losses(alpha2, alpha3, h/c, bx, h, M2, M3)

    Ytot = YP*fRe + YS + YTET + YTC

    return Ytot
