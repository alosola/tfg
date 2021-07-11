# -*- coding: utf-8 -*-
"""
Created on Sat Jul 10 20:27:14 2021

@author: alondra sola

Project: turbine design tool (TFG)
Program: inclusion of losses for 1D turbine
File: denton.py
Descpription: functions for calculation of tip clearance losses by means
of the Denton method (contour plots)
Source: Loss Mechanisms in Turbomachines

"""

import numpy as np                   # library for math and calculation functions
import matplotlib.pyplot as plt      # library for plots
import scipy.interpolate

def unshrouded(alphain, alphaout):
    x = np.genfromtxt('unshrouded.x.csv', delimiter=',')               # import alphaout vector
    y = np.flip(np.genfromtxt('unshrouded.y.csv', delimiter=','))      # import alphain vector
    z = np.flip(np.genfromtxt('unshrouded.Zm.csv', delimiter=','), 1)  # import loss data

    Y_spline =  scipy.interpolate.RectBivariateSpline(x, y, z)  # create spline for evaluation
    Y = Y_spline.ev(alphaout, alphain)

    return Y

Y_unshrouded = unshrouded(-40,10)

print(Y_unshrouded)