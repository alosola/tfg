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
from KO_profile_losses import profile_losses
from KO_secondary_losses import secondary_losses
import scipy.interpolate


def trailing_edge(M_out, gamma, t_o, alpha_in, alpha_out):

    dataset = np.genfromtxt('dataset_trailing_edge.csv', delimiter=',')
    x = dataset[:,0]
    y = dataset[:,1]

    spline = scipy.interpolate.InterpolatedUnivariateSpline(x,y)

    dphi2_0 = np.float64(spline(t_o))
    dphi2_2 = 0.5*dphi2_0

    dphi = dphi2_0 + abs(alpha_in/alpha_out)*(alpha_in/-alpha_out)*(dphi2_2 - dphi2_0)

    pw = -gamma/(gamma-1)
    YTET = ( (1-(gamma-1)/2*M_out**2*(1/(1-dphi)-1))**pw -1)/(1-(1+(gamma-1)/2*M_out**2)**pw)

    return YTET



def reynolds_correction(Re):
    # Rec = Reynolds number based on true chord and exit gas conditions

    mean = 2*10**5

    if Re <= mean:
        fRe = (Re/mean)**-0.4
    elif Re < 10**6:
        fRe = 1
    else:
        fRe = (Re/10**6)**-0.2

    return fRe


def calculateKO(K, s, alpha_in, alpha_out, c, bx, h, M_in, M_out, P_in, P_out, gamma, RHT, Re, t_o):
    """
    Function for calculation of total losses using Kacker-Okapuu formulation
    Y.P = profile losses
    Y.Re = Reynolds correction parameter
    Y.S = secondary losses
    Y.TET = trailing edge losses

    Total losses = Ytot = Y.P * Y.Re + Y.S + Y.TET
    """


    YP = profile_losses(s/c, alpha_out, alpha_in, M_in, M_out, K, P_in, P_out, gamma, RHT)
    fRe = reynolds_correction(Re)
    YS = secondary_losses(alpha_in, alpha_out, h/c, bx, h, M_in, M_out)
    YTET = trailing_edge(M_out, gamma, t_o, alpha_in, alpha_out)

    Ytot = YP*fRe + YS + YTET
    return Ytot


def kackerokapuu(one, two, thr, stator, rotor, gamma, RHT):
    """
    Function to run the Kacker-Okapuu losses on the stator and rotor
    """

    # run loss correlaations on stator
    K = 1.8
    stator.omegaKO = calculateKO(K, two.geo.s, abs(one.alpha), abs(two.alpha), two.geo.c, two.geo.bx, two.geo.h, one.vel.M, two.vel.M, one.P, two.P, gamma, RHT, two.Re, two.geo.to)
    K = 5.2
    rotor.omegaKO  = calculateKO(K, thr.geo.s, abs(two.beta), abs(thr.beta), thr.geo.c, thr.geo.bx, thr.geo.h, two.vel.Mr, thr.vel.Mr, two.P, thr.P, gamma, RHT, thr.Re, thr.geo.to)
