# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 23:48:18 2021

@author: alond
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate




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