# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 20:12:17 2021

@author: alondra sola

Project: turbine design tool (TFG)
Program: inclusion of losses for 1D turbine
File: check_limits.py

"""
import numpy as np
import tables_results as tab


def check_limits(two, thr):

    n = 8 # number of limit checks
    N = np.arange(n)

    limits = np.ones((2,n))

    limits[0:2,0] = [0,120]        # turning alpha (alpha 2 - alpha 3)  [deg]
    limits[0:2,1] = [0,120]        # turning beta  (beta 2 - beta 3)    [deg]
    limits[0:2,2] = [1,1.2]        # height ratio  (h3 / h2)            [-]
    limits[0:2,3] = [0.85,1.2]     # absolute mach at 2                 [-]
    limits[0:2,4] = [0,0.5]        # relative mach at 2                 [-]
    limits[0:2,5] = [0,0.45]       # absolute mach at 3                 [-]
    limits[0:2,6] = [0.85,1.2]     # relative mach at 3                 [-]
    limits[0:2,7] = [45,50]        # relative rotor inlet angle beta2   [-]

    values = np.ones((1,n))

    values[0,0] = np.degrees(abs(two.alpha - thr.alpha))
    values[0,1] = np.degrees(abs(two.beta - thr.beta))
    values[0,2] = thr.geo.A/two.geo.A
    values[0,3] = two.vel.M
    values[0,4] = two.vel.Mr
    values[0,5] = thr.vel.M
    values[0,6] = thr.vel.Mr
    values[0,7] = np.degrees(two.beta)

    checks = np.ones((1,n), dtype=object)

    for i in N:
        if (values[0,i] < limits[0,i]) or (values[0,i] > limits[1,i]):
            checks[0,i] = 'No'
        else:
            checks[0,i] = 'Yes'
        values[0,i] = round(values[0,i],2)


    tab.print_limits(limits, values, checks)


    return checks