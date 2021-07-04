# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 21:48:08 2021

@author: alond
"""

import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt      # library for plots

dataset = np.genfromtxt('dataset_machhub.csv', delimiter=',')


M2H = dataset[:,0]
ysol = dataset[:,1]



# def fdp(V, M2H, ysol):
#     dp = V[0]*(M2H-V[1])**V[2]
#     mm = np.amax(np.absolute(dp-ysol))
#     return np.array([mm,mm,mm])

# V = fsolve(fdp,np.array([7.5, 0.4, 1.75]),args=(M2H, ysol))

V = np.array([7.5, 0.4, 1.75])
dp = V[0]*(M2H-V[1])**V[2]


plt.plot(M2H,dp)
plt.plot(M2H,ysol)