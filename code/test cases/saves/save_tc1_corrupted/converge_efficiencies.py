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

    # 2.2
    # check limits
    lim_checks, lim_b, lim_ub, lim_values = check_limits(two, thr)

    # stator.eta = (one.T0 - thr.T0)/(one.T0 - (thr.Ts+thr.vel.V**2/2/cp))
    # rotor.eta = (one.T0 - thr.T0)/(one.T0 - thr.Ts)



def converge_efficiencies_limits(efficiencies_init, stator, rotor, one, two, thr, gamma, cp, R, GR, psi, DeltaH_prod, bounds_angles, RHT, mdot):

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

        return diff1 + diff2



    # def nlc_values(etas, stator, rotor, one, two, thr, gamma, cp, R, GR, psi, DeltaH_prod, bounds_angles, RHT, mdot):
    def nlc_values(etas):

        values = np.ones((1,8))

        values[0,0] = np.degrees(abs(two.alpha - thr.alpha))
        values[0,1] = np.degrees(abs(two.beta - thr.beta))
        values[0,2] = thr.geo.A/two.geo.A
        values[0,3] = two.vel.M
        values[0,4] = two.vel.Mr
        values[0,5] = thr.vel.M
        values[0,6] = thr.vel.Mr
        values[0,7] = np.degrees(two.beta)

        arr = values.tolist()

        # arr = [np.degrees(abs(two.alpha - thr.alpha)), np.degrees(abs(two.beta - thr.beta)), thr.geo.A/two.geo.A]

        return arr[0]



    def nlc_limits():

        limits = np.ones((2,8))

        limits[0:2,0] = [0,120]        # turning alpha (alpha 2 - alpha 3)  [deg]
        limits[0:2,1] = [0,120]        # turning beta  (beta 2 - beta 3)    [deg]
        limits[0:2,2] = [1,1.36]        # height ratio  (h3 / h2)            [-]
        limits[0:2,3] = [0.85,1.2]     # absolute mach at 2                 [-]
        limits[0:2,4] = [0,0.50]        # relative mach at 2                 [-]
        limits[0:2,5] = [0,0.60]       # absolute mach at 3                 [-]
        limits[0:2,6] = [0.85,1.2]     # relative mach at 3                 [-]
        limits[0:2,7] = [40,50]        # relative rotor inlet angle beta2   [-]

        lb = limits[0,:].tolist()
        ub = limits[1,:].tolist()

        # lb = [0,0,1]
        # ub = [120, 120,1.34]

        return lb, ub


    f_etas(efficiencies_init, stator, rotor, one, two, thr, gamma, cp, R, GR, psi, DeltaH_prod, bounds_angles, RHT, mdot)

    lb, ub = nlc_limits()

    # con = lambda x: [np.degrees(abs(two.alpha - thr.alpha)), np.degrees(abs(two.beta - thr.beta)), thr.geo.A/two.geo.A, two.vel.M, two.vel.Mr, thr.vel.M, thr.vel.Mr, np.degrees(two.beta) ]


    constraints = NonlinearConstraint(nlc_values, lb, ub)
    bounds = Bounds([0,0], [1, 1])

    res = minimize(f_etas, efficiencies_init, bounds=bounds, constraints=constraints, args=(stator, rotor, one, two, thr, gamma, cp, R, GR, psi, DeltaH_prod, bounds_angles, RHT, mdot), options={'maxiter': 200, 'disp': True})
    etas = res.x

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

    # 2.2
    # check limits
    lim_checks, lim_lb, lim_ub, lim_values = check_limits(two, thr)

    # stator.eta = (one.T0 - thr.T0)/(one.T0 - (thr.Ts+thr.vel.V**2/2/cp))
    # rotor.eta = (one.T0 - thr.T0)/(one.T0 - thr.Ts)

    return res
