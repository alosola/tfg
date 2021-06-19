# -*- coding: utf-8 -*-
"""
Created on Sat Jun 19 16:22:47 2021

@author: alondra sola

Project: turbine design tool (TFG)
Program: inclusion of losses for 1D turbine
File: kackerokapuu_plot.py
Descpription: simple depiction of surface plots fro kacker-okapuu correlations
(used to confirm that data has been imported well)

"""
## IMPORT NECESSARY LIBRARIES
import numpy as np                   # library for math and calculation functions
import matplotlib.pyplot as plt      # library for plots

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

X = np.genfromtxt('fig1_Xm.csv', delimiter=',')
Y = np.genfromtxt('fig1_Ym.csv', delimiter=',')
Z = np.flip(np.genfromtxt('fig1_Zm.csv', delimiter=','), 1)

ax.plot_surface(X, Y, Z,cmap='viridis', edgecolor='none')
ax.set_title('Surface plot')
plt.show()
