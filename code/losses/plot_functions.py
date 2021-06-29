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
import matplotlib.patches as patches

# change default font for plots
plt.rcParams["font.family"] = "sans"


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



def turbine_geometry(one, two, thr):
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, aspect='equal')

    width = 2
    Xstator = 0
    Xmiddle = Xstator + width
    Xrotor = Xstator + 2*width

    # points = [[Xstator, one.geo.Rh], [Xstator, one.geo.Rt], [Xmiddle, one.geo.Rt], [Xmiddle, one.geo.Rh], [Xstator, one.geo.Rh]] #the points to trace the edges.
    # stator= plt.Polygon(points,  fill=True, edgecolor='r', facecolor=(151/222, 178/222, 1))
    # ax2.add_patch(stator)
    # points = [[Xmiddle, two.geo.Rh], [Xmiddle, two.geo.Rt], [Xrotor, thr.geo.Rt], [Xrotor, thr.geo.Rh], [Xmiddle, two.geo.Rh]] #the points to trace the edges.
    # rotor= plt.Polygon(points,  fill=True, edgecolor=None, facecolor=(173/217, 1, 165/217))
    # ax2.add_patch(rotor)

    points = [[0.2, 0.4], [0.4, 0.8], [0.8, 0.8], [0.6, 0.4], [0.2,0.4]] #the points to trace the edges.
    polygon= plt.Polygon(points,  fill=None, edgecolor='r')
    ax2.add_patch(polygon)
    fig2.savefig('graph-geometry.png', dpi=90, bbox_inches='tight')
    plt.show()


def geometry(two, thr):

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, aspect='equal')

    Xrotor = 0.1
    width = 0.5

    rotor = [[Xrotor, two.geo.Rh], [Xrotor, two.geo.Rt], [Xrotor+width, thr.geo.Rt], [Xrotor+width, thr.geo.Rh], [Xrotor,two.geo.Rh]] #the points to trace the edges.
    # polygon= plt.Polygon(rotor,  fill=True, edgecolor=None, facecolor=(173/217, 1, 165/217))
    # ax2.add_patch(polygon)

    antirotor = rotor
    for i in np.array([0,1,2,3,4]):
        antirotor[i][1] = -rotor[i][1]

    antipolygon= plt.Polygon(antirotor,  fill=True, edgecolor=None, facecolor=(173/217, 1, 165/217))
    ax2.add_patch(antipolygon)


    fig2.savefig('graph_geometry.png', dpi=90, bbox_inches='tight')
    plt.show()


