# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 23:48:18 2021

@author: alond
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate




dataset = np.genfromtxt('dataset_trailing_edge.csv', delimiter=',')
t_o = dataset[:,0]
dphi = dataset[:,1]

s = interpolate.InterpolatedUnivariateSpline(t_o, dphi)
xnew = np.linspace(np.amin(t_o), np.amax(t_o), 500)
ynew = s(xnew)

Dphi = np.float64(s(0.2105))


plt.plot(t_o,dphi)
plt.plot(xnew, ynew)