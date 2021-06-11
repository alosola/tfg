# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 15:05:31 2021

@author: alondra sola

Project: turbine design tool (TFG)
Program: 1D analysis (imitating the "automobile gas generator" design tool)
File: solve_functions.py
Description: contains the implementation of all the functions needed for
convergence with non-linear solvers

"""

import numpy as np                   # library for math and calculation functions
import matplotlib as plt
from scipy.optimize import fsolve, least_squares
import aux_functions as f


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



def solve_angles(angles_init, one, two, thr, stator, rotor, gamma, cp, R, GR, psi, DeltaH_prod, bounds_angles):

    def f_angles(angles, one, two, thr, stator, rotor, gamma, cp, R, GR, psi, DeltaH_prod):

        alpha = angles[0]
        beta = angles[1]

        # assume an two.alpha and project velocities

        two.vel.Vx, two.vel.Vu = f.velocity_projections(two.vel.V, alpha)

        # assume a loading factor and calculate peripheral speed
        two.vel.U = np.sqrt(DeltaH_prod/psi)

        # velocity triangle (pithagoras)
        two.vel.Wu = two.vel.Vu - two.vel.U
        two.vel.Wx = two.vel.Vx
        two.vel.W = f.mag(two.vel.Wu, two.vel.Wx)

        # relative inlet angle and relative Mach number
        two.vel.Mr = two.vel.W/two.vel.a

        # relative total quantities at rotor inlet
        two.T0r, two.P0r = f.relative_temperature_pressure(two.T, two.P, two.vel.Mr, gamma)

        ############ ROTOR #############
        # in the rotor, the relative total temperature is constant (se conservan rotalpías)
        thr.T0r = two.T0r

        # ideal outlet temperature, real outlet temperature
        # assuming the expanse to thr.P, calculated with the assumed M3
        thr.Ts = f.P2T(thr.T0r, thr.P, two.P0r, gamma)                  # = thr.T0r*(thr.P/two.P0r)**((gamma-1)/gamma)
        thr.T = f.stage_efficiency(rotor.eta, thr.T0r, thr.Ts)       # = thr.T0r - rotor.eta*(thr.T0r - thr.Ts)

        # velocities
        thr.vel.W = f.isen_velocity(cp, thr.T0r, thr.T)                  # = np.sqrt(2*cp*(thr.T0r - thr.T))

        # velocity projections, axial and tangential
        thr.vel.Wx, thr.vel.Wu = f.velocity_projections(thr.vel.W, beta)

        # constant radius, constant peripheral speed
        thr.vel.U = two.vel.U

        # tangential outlet speed
        thr.vel.Vu = thr.vel.Wu + thr.vel.U
        # V3u_euler = -DeltaH_prod/two.vel.U + two.vel.Vu


        # find work, check for angle iteration
        DeltaH_calculated = two.vel.U*(two.vel.Vu - thr.vel.Vu)

        diff = DeltaH_calculated - DeltaH_prod


        return np.array([diff[0], 0])

    # angles  = fsolve(f_angles,angles_init,args=(one, two, thr, stator, rotor, gamma, cp, R, GR, psi, DeltaH_prod))
    results  = least_squares(f_angles,angles_init,args=(one, two, thr, stator, rotor, gamma, cp, R, GR, psi, DeltaH_prod),bounds=bounds_angles, method='trf')
    # (the least_squares method allows for the inclusion of bounds)

    # if results.success != False:
    #     print('ERROR: The angle conversion was not successful (function solve_angles in solve_functions.py)')

    return results.x

