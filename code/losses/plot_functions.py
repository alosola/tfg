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

    X = 1.2*np.amax(arrayT)
    Y = 1.2*V3x
    x = np.linspace(-X,X,100)
    y = np.linspace(-0.1*Y,Y,100)

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


def geometry(one, two, thr):

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, aspect=0.9)
    ax2.set_ylabel('Radial distance [m]')
    ax2.set_xlabel('Turbine axis [m]')

    Xstator = 0
    Wstator = two.geo.c
    Xrotor = Xstator + Wstator + Wstator/5
    Wrotor = thr.geo.c

    background = [[Xstator, -1.1*thr.geo.Rt], [Xstator, 1.1*thr.geo.Rt], [Xrotor+Wrotor, 1.1*thr.geo.Rt], [Xrotor+Wrotor, -1.1*thr.geo.Rt], [Xstator,-1.1*thr.geo.Rt]] #the points to trace the edges.
    polygon_background = plt.Polygon(background,  fill=True, edgecolor=None, facecolor='#ededed')
    ax2.add_patch(polygon_background)


    stator = [[Xstator, one.geo.Rh], [Xstator, one.geo.Rt], [Xstator+Wstator, two.geo.Rt], [Xstator+Wstator, two.geo.Rh], [Xstator,one.geo.Rh]] #the points to trace the edges.
    polygon_stator = plt.Polygon(stator,  fill=True, edgecolor=None, facecolor='#7896bf')
    ax2.add_patch(polygon_stator)


    rotor = [[Xrotor, two.geo.Rh], [Xrotor, two.geo.Rt], [Xrotor+Wrotor, thr.geo.Rt], [Xrotor+Wrotor, thr.geo.Rh], [Xrotor,two.geo.Rh]] #the points to trace the edges.
    polygon_rotor = plt.Polygon(rotor,  fill=True, edgecolor=None, facecolor='#7cbf78')
    ax2.add_patch(polygon_rotor)

    antistator = [[Xstator, -one.geo.Rt], [Xstator, -one.geo.Rh], [Xstator+Wstator, -two.geo.Rh], [Xstator+Wstator, -two.geo.Rt], [Xstator,-one.geo.Rt]] #the points to trace the edges.
    polygon_antistator = plt.Polygon(antistator,  fill=True, edgecolor=None, facecolor='#7896bf')
    ax2.add_patch(polygon_antistator)

    antirotor = [[Xrotor, -two.geo.Rt], [Xrotor, -two.geo.Rh], [Xrotor+Wrotor, -thr.geo.Rh], [Xrotor+Wrotor, -thr.geo.Rt], [Xrotor,-two.geo.Rt]] #the points to trace the edges.
    polygon_antirotor = plt.Polygon(antirotor,  fill=True, edgecolor=None, facecolor='#7cbf78')
    ax2.add_patch(polygon_antirotor)

    # plt.axis([-0.3*Xrotor,(Xrotor + Wrotor + 0.3*Xrotor),-1.3*thr.geo.Rt,1.3*thr.geo.Rt])
    plt.axis([-0.3*Xrotor,(2.6*thr.geo.Rt - 0.3*Xrotor)/2,-1.3*thr.geo.Rt,1.3*thr.geo.Rt])

    Rm_top = np.ones(100)*two.geo.Rm
    Rm_bot = -Rm_top
    Xrm = np.linspace(0,Xstator+Wstator,100)
    plt.plot(Xrm, Rm_top, Xrm, Rm_bot, color='black', linewidth = 0.5)

    Xrm = np.linspace(Xrotor,Xrotor+Wrotor,100)
    plt.plot(Xrm, Rm_top, Xrm, Rm_bot, color='black', linewidth = 0.5)

    Zero = np.zeros(100)
    Xz = np.linspace(0,Xrotor+Wrotor,100)
    plt.plot(Xz, Zero, '#6b6b6b', linewidth = 7)

    fig2.savefig('graph_geometry.png', dpi=90, bbox_inches='tight')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.axis('scaled')
    plt.show()
