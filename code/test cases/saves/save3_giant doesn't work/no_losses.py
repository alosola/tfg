# -*- coding: utf-8 -*-
"""
Created on Sun Jul  4 15:08:15 2021

@author: alond
@author: alondra sola

Project: turbine design tool (TFG)
Program: inclusion of losses for 1D turbine
File: no_losses.py

"""

import numpy as np
import aux_functions as f
from converge_mach3 import converge_mach3
from check_limits import check_limits


def no_losses(stator, rotor, one, two, thr, gamma, cp, R, GR, psi, DeltaH_prod, bounds_angles, RHT, mdot):

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
    f.geometry(RHT, one, two, thr, mdot, R, gamma, cp)

    # 2.2
    # check limits
    pass_limits = check_limits(two, thr)

    stator.eta = (one.T0 - thr.T0)/(one.T0 - (thr.Ts+thr.vel.V**2/2/cp))
    rotor.eta = (one.T0 - thr.T0)/(one.T0 - thr.Ts)


