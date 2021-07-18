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

GR = np.linspace(0.1,0.8,N)
results = np.zeros((N,5))

for i in range(N):
    results[i,1:5] = runloss(GR[i])
    print(i,results[i,1:5])

results[:,0] = results[:,1]*results[:,2]

# %%
# res = np.zeros((10,1))
# gr = np.zeros((10,1))

# for i in range(10):
#     res[i] = results[i*5,0]
#     gr[i] = GR[i*5]

# %%
fig, (ax1, ax2) = plt.subplots(2)
ax1.plot(GR, results[:,0:3],label = [])
ax1.legend([r'$\eta_{total}$',r'$\eta_{stator}$',r'$\eta_{rotor}$'],loc='center left', bbox_to_anchor=(1, 0.5))

ax2.plot(GR, results[:,3:6])
ax2.legend([r'$Y_{stator}$',r'$Y_{rotor}$'],loc='center left', bbox_to_anchor=(1, 0.5))
# ax2.set_ylabel('')

ax2.set_xlabel(r'Degree of reaction $GR$')
plt.show()

# %%
######################### stage loading

# number of steps
N = 10

psi = np.linspace(0.8, 3,N)
results = np.zeros((N,5))

for i in range(N):
    results[i,1:5] = runloss(psi[i])
    print(i,results[i,1:5])

results[:,0] = results[:,1]*results[:,2]


fig, (ax1, ax2) = plt.subplots(2)
ax1.plot(psi, results[:,0:3],label = [])
ax1.legend([r'$\eta_{total}$',r'$\eta_{stator}$',r'$\eta_{rotor}$'],loc='center left', bbox_to_anchor=(1, 0.5))

ax2.plot(psi, results[:,3:6])
ax2.legend([r'$Y_{stator}$',r'$Y_{rotor}$'],loc='center left', bbox_to_anchor=(1, 0.5))
# ax2.set_ylabel('')

ax2.set_xlabel(r'Stage loading factor $\psi$')
plt.show()