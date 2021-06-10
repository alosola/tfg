# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 15:05:31 2021

@author: alondra sola

Project: turbine design tool (TFG)
Program: 1D analysis (imitating the "automobile gas generator" design tool)
File: converge_rhox.py
Description: contains the implementation of all the calculations within the
rhox convergence loop as well as the function with the non-linear solver

"""

import numpy as np                   # library for math and calculation functions
from scipy.optimize import fsolve


def solve_rhox(rhox_init, rho0x, T0x, mdot, Ay, gamma, cp):
    """
    Function to find the density rho for a certain point x, from the
    temperature, mass flux, and area data.
    Using an initial guess rhox_init, which is always less than rho0x
    """

    def f_rhox(rhox, rho0x, T0x, mdot, Ay, gamma, cp):
        Vx = np.sqrt(2*cp*T0x*(1-(rhox/rho0x)**(gamma-1)))
        Ax = mdot/rhox/Vx
        return Ay-Ax

    rhox = fsolve(f_rhox,rhox_init,args=(rho0x, T0x, mdot, Ay, gamma, cp))[0]

    return rhox;
