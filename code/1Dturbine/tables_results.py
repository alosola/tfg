# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 23:16:50 2021

@author: alondra sola

Project: turbine design tool (TFG)
Program: 1D analysis (imitating the "automobile gas generator" design tool)
File: tables_results.py
Description: contains functions for plotting results to cmd window in tables

"""

import numpy as np
import pandas as pd


def convergence_parameters(one, two, thr, Mach3_init, alpha2_init, beta3_init, DeltaH_prod, DeltaH_calculated):
    data1 = {'Parameter':  ['Mach3 initial guess', 'Mach3 final value', 'Alpha2 initial guess', 'Alpha2 final value', 'Beta3 initial guess', 'Beta3 final value', 'DeltaH necessary', 'DeltaH calculated'],
        'Value':     [Mach3_init, thr.vel.M, np.degrees(alpha2_init), np.degrees(two.alpha), np.degrees(beta3_init), np.degrees(thr.beta), DeltaH_prod/1000, DeltaH_calculated/1000],
        'Unit':      ['-', '-',  'deg', 'deg',  'deg', 'deg', 'kJ/kg', 'kJ/kg']}
    df1 = pd.DataFrame (data1, columns = ['Parameter','Value','Unit'])
    print (df1)
