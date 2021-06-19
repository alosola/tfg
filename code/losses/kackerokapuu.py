# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 12:22:02 2021

@author: alondra sola

Project: turbine design tool (TFG)
Program: inclusion of losses for 1D turbine
File: kackerokapuu.py
Descpription: functions for calculation of pressure loss coefficient by means
of the Ainley/Mathieson/Dunham/Crane/Kacker/Okapuu method
Souce: A Mean Line Prediction Method for Axial Flow Turbine Efficiency

"""

## IMPORT NECESSARY LIBRARIES
import numpy as np                   # library for math and calculation functions
import pandas as pd                  # library for tables and data import
import matplotlib.pyplot as plt      # library for plots


# change default font for plots
plt.rcParams["font.family"] = "sans"

## IMPORT FUNCTIONS AND DEFINITIONS FROM OTHER PYTHON FILES
from definitions import pressure_losses


P = np.genfromtxt('fig1.csv', delimiter=',')
plt.scatter(P[:,0],P[:,1],s=2)

# def kackerokapuu():

#     Y = pressure_losses()

#     Y = total_pressure_loss_coefficient()

#     return Y.tot


# def profile_losses():

#     P = pd.read_csv("fig1.csv")