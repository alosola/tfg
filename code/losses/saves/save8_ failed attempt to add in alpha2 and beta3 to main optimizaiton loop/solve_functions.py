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
from scipy.optimize import fsolve, least_squares
import aux_functions as f


def solve_rhox(rho_init, rho0, T0, mdot, Ay, gamma, cp):
    """
    Function to find the density rho for a certain point x, from the
    temperature, mass flux, and area data.
    Using an initial guess rhox_init, which is always less than rho0x
    """

    def f_rho(rho, rho0, T0, mdot, Ay, gamma, cp):
        V = np.sqrt(2*cp*T0*(1-(rho/rho0)**(gamma-1)))
        A = mdot/rho/V
        return Ay-A

    rho = fsolve(f_rho,rho_init,args=(rho0, T0, mdot, Ay, gamma, cp))[0]
    V = np.sqrt(2*cp*T0*(1-(rho/rho0)**(gamma-1)))
    A = mdot/rho/V

    return rho, V, A;



def solve_angles(one, two, thr, stator, rotor, gamma, cp, R, GR, psi, DeltaH_prod):

    two.vel.Vx, two.vel.Vu = f.velocity_projections(two.vel.V, two.alpha)

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
    thr.vel.Wx, thr.vel.Wu = f.velocity_projections(thr.vel.W, thr.beta)

    # constant radius, constant peripheral speed
    thr.vel.U = two.vel.U

    # tangential outlet speed
    thr.vel.Vu = thr.vel.Wu + thr.vel.U
    # V3u_euler = -DeltaH_prod/two.vel.U + two.vel.Vu




