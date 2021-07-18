# -*- coding: utf-8 -*-
"""
Created on Fri May 14 11:11:08 2021

@author: alondra sola

Project: turbine design tool (TFG)
Program: inclusion of losses for 1D turbine
File: main.py

"""

## IMPORT NECESSARY LIBRARIES
import numpy as np                   # library for math and calculation functions

## IMPORT FUNCTIONS AND DEFINITIONS FROM OTHER PYTHON FILES
from runloss import runloss
import matplotlib.pyplot as plt

# %%
# number of steps
N = 10

GR = np.linspace(0.3,0.8,N)
results = np.zeros((N,5))

for i in range(N):
    results[i,1:5] = runloss(GR[i])
    print(i,results[i,1:5])

results[:,0] = results[:,1]*results[:,2]

# %%

ax1 = plt.subplot(131)
ax1.plot(GR, results[:,0],color='k')
ax1.set_ylabel(r'Total efficiency $\eta$')

ax2 = plt.subplot(132)
ax2.plot(GR, results[:,1:3])
ax2.set_ylabel(r'Rotor/stator efficiency $\eta_x$')

ax3 = plt.subplot(133)
ax3.plot(GR, results[:,3:6])
ax3.set_ylabel(r'Rotor/stator pressure loss $Y_x$')

ax3.legend([r'$Y_{stator}$',r'$Y_{rotor}$'],loc='center left', bbox_to_anchor=(1, 0.5))
# ax2.set_ylabel('')

ax2.set_xlabel(r'Degree of reaction $GR$')
plt.show()
