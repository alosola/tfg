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


# def turning(two, thr):
#     """
#     use this function to define the limits of what is acceptable for
#     turning values (will be a design constraint in main program)

#     NOTE: need to find references for these values.

#     """

#     # MAX LIMIT for ALPHA TURNING:
#     limit_alpha = np.radians(120)

#     # MAX LIMIT for BETA TURNING:
#     limit_beta = np.radians(100)


#     turning_alpha = abs(thr.alpha - two.alpha)
#     turning_beta = abs(thr.beta - two.beta)

#     result = np.ones((1,2))
#     if turning_alpha > limit_alpha:
#         print('Alpha turning angle is outside acceptable bounds (', round(np.degrees(turning_alpha),2),'is above', round(np.degrees(limit_alpha),2), ')' )
#         result[0,0] = False
#     if  turning_beta > limit_beta:
#         print('Beta turning angle is outside acceptable bounds (', round(np.degrees(turning_beta),2),'is above', round(np.degrees(limit_beta),2), ')' )
#         result[0,1] = False

#     return result


# def area_ratio(two, thr):
#     """
#     use this function to define the limits of what is acceptable for
#     area ratio (will be a design constraint in main program)

#     NOTE: need to find references for this value.

#     """

#     # MAX LIMIT for AREA RATIO:
#     limit = 1.2

#     ratio = thr.geo.A/two.geo.A

#     result = True
#     if ratio > limit:
#         print('Area ratio is outside acceptable bounds (', round(ratio,2),'is above', round(limit,2), ')' )
#         result = False

#     return result


# def machs(two, thr):
#     """
#     use this function to define the limits of what is acceptable for
#     mach numbers in turbine (will be a design constraint in main program)

#     NOTE: need to find references for these values.

#     """

#     # MIN/MAX LIMIT for ABSOLUTE MACH at point 2
#     limit_M2 = np.array([0.85,1.2])

#     # MAX LIMIT for RELATIVE MACH at point 2
#     limit_M2r = 0.5

#     # MIN/MAX LIMIT for RELATIVE MACH at point 3
#     limit_M3r = np.array([0.85,1.2])


#     result = np.ones((1,3))
#     if two.vel.M < limit_M2[0] or two.vel.M > limit_M2[1]:
#         print('Rotor entry absolute Mach is outside acceptable bounds')
#         result[0,0] = False
#     if two.vel.Mr > limit_M2r:
#         print('Rotor entry relative Mach is above acceptable bounds')
#         result[0,1] = False
#     if thr.vel.Mr < limit_M3r[0] or thr.vel.Mr > limit_M3r[1]:
#         print('Rotor outlet relative Mach is outside acceptable bounds')
#         result[0,2] = False

#     return result


# def check_limits(two, thr):

#     results = np.zeros([1,6], dtype=bool)

#     results[0,0:2] = turning(two, thr)
#     results[0,2] = area_ratio(two, thr)
#     results[0,3:6] = machs(two, thr)

#     pass_all = np.amin(results)
#     if pass_all == 0:
#         print('Since one or more limits were not met, the one-dimensional design must be reesablished.')

#     return results



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