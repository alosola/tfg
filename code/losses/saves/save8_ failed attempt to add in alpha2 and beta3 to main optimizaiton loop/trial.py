# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 23:48:18 2021

@author: alond
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate



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


lb = limits[0,:].tolist()
ub = limits[1,:].tolist()

limits = limits.tolist()


# RHT = np.linspace(0.5, 1, 100)
# MH = np.zeros((1,100))
# K = 1.8

# MH = (1+K*np.absolute(RHT-1)**2.2)

# plt.plot(RHT,MH)

# K = 5.2

# MH = (1+K*np.absolute(RHT-1)**2.2)

# plt.plot(RHT,MH)







# dataset = np.genfromtxt('dataset_trailing_edge.csv', delimiter=',')
# t_o = dataset[:,0]
# dphi = dataset[:,1]

# s = interpolate.InterpolatedUnivariateSpline(t_o, dphi)
# xnew = np.linspace(np.amin(t_o), np.amax(t_o), 500)
# ynew = s(xnew)

# Dphi = np.float64(s(0.2105))


# plt.plot(t_o,dphi)
# plt.plot(xnew, ynew)