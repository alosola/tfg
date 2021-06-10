# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 15:15:05 2021

@author: alondra sola

Project: turbine design tool (TFG)
Program: 1D analysis (imitating the "automobile gas generator" design tool)
File: plot_functions.py
Description: contains definitions of functions used to plot results

"""

## IMPORT NECESSARY LIBRARIES
import numpy as np                   # library for math and calculation functions
import matplotlib.pyplot as plt      # library for plots


# change default font for plots
plt.rcParams["font.family"] = "cursive"


def velocity_triangle(V2u, V2x, V3u, V3x, W2u, W2x, W3u, W3x):
    y0 = W3x-V2x
    arrayV2 = np.array([[0,V2u],[V2x + y0, y0]])
    arrayW2 = np.array([[0,W2u],[W2x + y0, y0]])
    arrayU2 = np.array([[W2u,V2u],[y0,y0]])
    arrayV3 = np.array([[0,V3u],[V3x,0]])
    arrayW3 = np.array([[0,W3u],[W3x,0]])
    arrayU3 = np.array([[V3u,W3u],[0,0]])

    arrayN = np.array(['V2', 'W2', 'U2', 'V3', 'W3', 'U3'])
    arrayT = [arrayV2, arrayW2, arrayU2, arrayV3, arrayW3, arrayU3]


    x = np.linspace(-1000,1000,100)
    y = np.linspace(-10,300,100)

    i = 0
    for n in arrayT:
        plt.plot(n[0],n[1],label=arrayN[i])
        plt.axis([max(x),min(x),min(y),max(y)])
        plt.legend()
        i+=1
    plt.savefig('velocity_triangle.png', dpi=1000)
    plt.show()