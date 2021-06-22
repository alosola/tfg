# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 19:41:47 2021

@author: alondra sola

Project: turbine design tool (TFG)
Program: inclusion of losses for 1D turbine
File: KC_secondary_losses.py
Descpription: functions for calculation of secondary losses within Kacker-Okapuu

"""

## IMPORT NECESSARY LIBRARIES
import numpy as np                   # library for math and calculation functions

from KC_profile_losses import losses_KP

def secondary_losses_fAR(h_c):

    if h_c > 2:
        fAR = 1/h_c
    else:
        fAR = (1-0.25*np.sqrt(2-h_c))/h_c

    return fAR


def secondary_losses_AMDC(alpha2, alpha3, h_c):

    fAR = secondary_losses_fAR(h_c)

    alpham = np.arctan(0.5*(np.tan(alpha2)-np.tan(alpha3)))
    C1 = 2*(np.tan(alpha2)+np.tan(alpha3))*np.cos(alpham)

    YS_AMDC = 0.0334*fAR*np.cos(alpha3)/np.cos(alpha2)*C1**2*np.cos(alpha3)**2/np.cos(alpham)**3

    return YS_AMDC


def secondary_losses_Ks(bx, h, M2, M3):

    K3 = (bx/h)**2
    KP = losses_KP(M2, M3)

    Ks = 1 - K3*(1-KP)

    return Ks


def secondary_losses(alpha2, alpha3, h_c, bx, h, M2, M3):


    YS_AMDC = secondary_losses_AMDC(alpha2, alpha3, h_c)
    Ks = secondary_losses_Ks(bx, h, M2, M3)


    YS = 1.2*YS_AMDC*Ks

    return YS
