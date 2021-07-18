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
from check_limits import check_limits
from scipy.optimize import fsolve, minimize, NonlinearConstraint, Bounds
from kackerokapuu import kackerokapuu

# from trust_constr import minimize, NonlinearConstraint, LinearConstraint, Bounds, check_grad


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
        f.geometry(RHT, one, two, thr, mdot, R, gamma, cp)

        # calculate reynolds numbers
        f.reynolds(two, thr)

        # use kacker-okapuu to calculate losses in stator and rotor
        kackerokapuu(one, two, thr, stator, rotor, gamma, RHT)

        diff1 = abs(stator.omega - stator.omegaKO)
        diff2 = abs(rotor.omega  - rotor.omegaKO)

        return np.array([diff1, diff2])


    stator.eta, rotor.eta = fsolve(f_etas,efficiencies_init,args=(stator, rotor, one, two, thr, gamma, cp, R, GR, psi, DeltaH_prod, bounds_angles, RHT, mdot))


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


