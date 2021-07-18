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

from KO_profile_losses import losses_KP

def secondary_losses_fAR(h_c):

    if h_c > 2:
        fAR = 1/h_c
    else:
        fAR = (1-0.25*np.sqrt(2-h_c))/h_c

    return fAR


def secondary_losses_alpham(alpha_in, alpha_out):
    alpham = np.arctan(0.5*(np.tan(alpha_in)-np.tan(alpha_out)))
    return alpham


def secondary_losses_AMDC(alpha_in, alpha_out, h_c):

    fAR = secondary_losses_fAR(h_c)

    alpham = secondary_losses_alpham(alpha_in, alpha_out)
    C1 = 2*(np.tan(alpha_in)+np.tan(alpha_out))*np.cos(alpham)

    YS_AMDC = 0.0334*fAR*np.cos(alpha_out)/np.cos(alpha_in)*C1**2*np.cos(alpha_out)**2/np.cos(alpham)**3

    return YS_AMDC


def secondary_losses_Ks(bx, h, M_in, M_out):

    K3 = (bx/h)**2
    KP = losses_KP(M_in, M_out)

    Ks = 1 - K3*(1-KP)

    return Ks


def secondary_losses(alpha_in, alpha_out, h_c, bx, h, M_in, M_out):


    YS_AMDC = secondary_losses_AMDC(alpha_in, alpha_out, h_c)
    Ks = secondary_losses_Ks(bx, h, M_in, M_out)


    YS = 1.2*YS_AMDC*Ks

    return YS