# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 23:48:18 2021

@author: alond
"""

import numpy as np
import matplotlib.pyplot as plt



def secondary_losses_fAR(h_c):

    if h_c > 2:
        fAR = 1/h_c
    else:
        fAR = (1-0.25*np.sqrt(2-h_c))/h_c

    return fAR


A = np.linspace(0.9,2.5,num=100)

B = np.array([])

for i in A:
    B = np.append(B,secondary_losses_fAR(1/i))


plt.plot(A,B)