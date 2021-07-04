# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 18:49:13 2021

@author: alondra sola

Project: turbine design tool (TFG)
Program: inclusion of losses for 1D turbine
File: main.py

"""

import numpy as np
import aux_functions as f
from converge_mach3 import converge_mach3
from solve_functions import solve_geometry
from check_limits import check_limits
from scipy.optimize import fsolve
from kackerokapuu import kackerokapuu

def converge_efficiencies(efficiencies_init, stator, rotor, one, two, thr, gamma, cp, R, GR, psi, DeltaH_prod, bounds_angles, RHT, mdot):

    def f_etas(etas, stator, rotor, one, two, thr, gamma, cp, R, GR, psi, DeltaH_prod, bounds_angles, RHT, mdot):
        stator.eta = etas[0]
        rotor.eta = etas[1]


        # 1.
        # mach 3 converge until calculated value meets the initial guess
        # inside of this function, alpha2 and beta 3 are also set so that the DeltaH is met
        one, two, thr, stator, rotor = converge_mach3(one, two, thr, stator, rotor, gamma, cp, R, GR, psi, DeltaH_prod, bounds_angles)

        # 1.2
        # calculate loss coefficients
        stator = f.losses(two.vel.V, two.vel.Vs, one.P0, two.P0, two.P, stator)
        rotor = f.losses(thr.vel.W, thr.vel.Ws, two.P0r, thr.P0r, thr.P, rotor)

        # 1.4
        # compute total condtions
        thr.T0, thr.P0 = f.total_conditions(thr.T, thr.vel.V, thr.P, cp, gamma)


        # 2.
        # establish geometry after assuming RHT
        solve_geometry(RHT, one, two, thr, mdot, R, gamma, cp)

        # efficiency calculations
        stator.eta = (one.T0 - thr.T0)/(one.T0 - (thr.Ts+thr.vel.V**2/2/cp))
        rotor.eta = (one.T0 - thr.T0)/(one.T0 - thr.Ts)

        # calculate reynolds numbers
        f.reynolds(two, thr)

        # use kacker-okapuu to calculate losses in stator and rotor
        stator.omegaKC = kackerokapuu('stator', two.geo.s, abs(one.alpha), abs(two.beta), two.geo.c, two.geo.bx, two.geo.h, one.vel.M, two.vel.M, one.P, two.P, gamma, RHT, two.Re, two.geo.to)
        rotor.omegaKC  = kackerokapuu('rotor',  thr.geo.s, abs(two.alpha), abs(thr.beta), thr.geo.c, thr.geo.bx, thr.geo.h, two.vel.Mr, thr.vel.Mr, two.P, thr.P, gamma, RHT, thr.Re, thr.geo.to)


        diff1 = abs(stator.omega - stator.omegaKC)
        diff2 = abs(rotor.omega  - rotor.omegaKC)

        return np.array([diff1, diff2])


    etas = fsolve(f_etas,efficiencies_init,args=(stator, rotor, one, two, thr, gamma, cp, R, GR, psi, DeltaH_prod, bounds_angles, RHT, mdot))

    stator.eta = etas[0]
    rotor.eta = etas[1]


    # 1.
    # mach 3 converge until calculated value meets the initial guess
    # inside of this function, alpha2 and beta 3 are also set so that the DeltaH is met
    one, two, thr, stator, rotor = converge_mach3(one, two, thr, stator, rotor, gamma, cp, R, GR, psi, DeltaH_prod, bounds_angles)

    # 1.2
    # calculate loss coefficients
    stator = f.losses(two.vel.V, two.vel.Vs, one.P0, two.P0, two.P, stator)
    rotor = f.losses(thr.vel.W, thr.vel.Ws, two.P0r, thr.P0r, thr.P, rotor)

    # 1.4
    # compute total condtions
    thr.T0, thr.P0 = f.total_conditions(thr.T, thr.vel.V, thr.P, cp, gamma)


    # 2.
    # establish geometry after assuming RHT
    solve_geometry(RHT, one, two, thr, mdot, R, gamma, cp)

    # 2.2
    # check limits
    pass_limits = check_limits(two, thr)

    stator.eta = (one.T0 - thr.T0)/(one.T0 - (thr.Ts+thr.vel.V**2/2/cp))
    rotor.eta = (one.T0 - thr.T0)/(one.T0 - thr.Ts)


