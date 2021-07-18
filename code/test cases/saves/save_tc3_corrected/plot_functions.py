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
from mpl_toolkits.axes_grid.axislines import Subplot

# change default font for plots
plt.rcParams["font.family"] = "sans"


def velocity_triangle(V2u, V2x, V3u, V3x, W2u, W2x, W3u, W3x):
    if W3x>V2x:
        y0 = W3x-V2x
        z0 = 0
    else:
        y0 = 0
        z0 = V2x-W3x
    arrayV2 = np.array([[0,V2u],[V2x + y0, y0]])
    arrayW2 = np.array([[0,W2u],[W2x + y0, y0]])
    arrayU2 = np.array([[W2u,V2u],[y0,y0]])
    arrayV3 = np.array([[0,V3u],[V3x + z0,z0]])
    arrayW3 = np.array([[0,W3u],[W3x + z0,z0]])
    arrayU3 = np.array([[V3u,W3u],[z0,z0]])

    arrayN = np.array([r'$V_2$', r'$W_2$', r'$U_2$', r'$V_3$', r'$W_3$', r'$U_3$'])
    arrayT = [arrayV2, arrayW2, arrayU2, arrayV3, arrayW3, arrayU3]

    X = 1.2*max([V2u, W2u, V3u, W3u])
    Y = 1.2*max([W3x,V2x])
    x = np.linspace(-X,X,100)
    y = np.linspace(-0.1*Y,Y,100)

    i = 0
    for n in arrayT:
        plt.plot(n[0],n[1],label=arrayN[i])
        plt.axis([max(x),min(x),1.3*min(y),max(y)])
        plt.legend()
        i+=1
    plt.savefig('velocity_triangle.png', dpi=1000)
    plt.xlabel(r'Tangential velocity $v_u$ [m/s]')
    plt.ylabel(r'Axial velocity $v_x$ [m/s]')
    plt.show()



def geometry(one, two, thr):

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.set_ylabel('Radial distance [m]')
    ax2.set_xlabel('Turbine axis [m]')

    Xstator = 0
    Wstator = np.asscalar(two.geo.c)
    Xrotor = Xstator + Wstator + Wstator/5
    Wrotor = thr.geo.c

    High = 1.1*thr.geo.Rt

    background = [[Xstator, -High], [Xstator, High], [Xrotor+Wrotor, High], [Xrotor+Wrotor, -High], [Xstator,-High]] #the points to trace the edges.
    polygon_background = plt.Polygon(background,  fill=True, edgecolor=None, facecolor='#ededed')
    ax2.add_patch(polygon_background)

    h = one.geo.Rh
    t = two.geo.Rt

    stator = [[Xstator, h], [Xstator, t], [Xstator+Wstator, t], [Xstator+Wstator, h], [Xstator,h]] #the points to trace the edges.
    polygon_stator = plt.Polygon(stator,  fill=True, edgecolor=None, facecolor='#7896bf')
    ax2.add_patch(polygon_stator)

    h3 = thr.geo.Rh
    t3 = thr.geo.Rt

    rotor = [[Xrotor, h], [Xrotor, t], [Xrotor+Wrotor, t3], [Xrotor+Wrotor, h3], [Xrotor,h]] #the points to trace the edges.
    polygon_rotor = plt.Polygon(rotor,  fill=True, edgecolor=None, facecolor='#7cbf78')
    ax2.add_patch(polygon_rotor)

    antistator = [[Xstator, -t], [Xstator, -h], [Xstator+Wstator, -h], [Xstator+Wstator, -t], [Xstator,-t]]
    polygon_antistator = plt.Polygon(antistator,  fill=True, edgecolor=None, facecolor='#7896bf')
    ax2.add_patch(polygon_antistator)

    antirotor = [[Xrotor, -t], [Xrotor, -h], [Xrotor+Wrotor, -h3], [Xrotor+Wrotor, -t3], [Xrotor,-t]] #the points to trace the edges.
    polygon_antirotor = plt.Polygon(antirotor,  fill=True, edgecolor=None, facecolor='#7cbf78')
    ax2.add_patch(polygon_antirotor)

    Rm_top = np.ones(100)*two.geo.Rm
    Rm_bot = -Rm_top
    Xrm = np.linspace(0,Xstator+Wstator,100)
    plt.plot(Xrm, Rm_top, Xrm, Rm_bot, color='black', linewidth = 0.5)

    Xrm = np.linspace(Xrotor,Xrotor+Wrotor,100)
    plt.plot(Xrm, Rm_top, Xrm, Rm_bot, color='black', linewidth = 0.5)

    Zero = np.zeros(100)
    Xz = np.linspace(0,Xrotor+Wrotor,100)
    plt.plot(Xz, Zero, '#6b6b6b', linewidth = 7)

    plt.gca().set_aspect('equal', adjustable='box')
    plt.axis('scaled')
    fig2.patch.set_visible(False)
    plt.axis('off')
    plt.show()

    fig2.savefig('graph_geometry.png', bbox_inches='tight')



