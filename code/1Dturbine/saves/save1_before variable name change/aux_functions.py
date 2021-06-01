# -*- coding: utf-8 -*-
"""
Created on Thu May 27 19:15:23 2021

@author: alond
"""

import numpy as np                   # library for math and calculation functions
from scipy.optimize import newton, fsolve


def rhox_converge(rhox_init, rho0x, T0x, mdot, Ay, gamma, cp):
    """
    Function to find the density rho for a certain point x, from the
    temperature, mass flux, and area data.
    Using an initial guess rhox_init, which is always less than rho0x
    """

    def f(rhox, rho0x, T0x, mdot, Ay, gamma, cp):
        Vx = np.sqrt(2*cp*T0x*(1-(rhox/rho0x)**(gamma-1)))
        Ax = mdot/rhox/Vx
        return Ay-Ax

    rhox = fsolve(f,rhox_init,args=(rho0x, T0x, mdot, Ay, gamma, cp))[0]

    return rhox;


def blade_geometry(mdot, rho, Vx, RHT):

    Rtip = np.sqrt(mdot/(np.pi*rho*Vx*(1-RHT**2)))
    Rhub = Rtip*RHT
    h = Rtip-Rhub
    Rmean = Rtip - h/2
    Dmean = Rmean*2
    A = np.pi*(Rtip**2-Rhub**2)

    return Rtip, Rhub, h, Rmean, Dmean, A




# def Newton_Raphson(func, M, interations = 10,*args):

#         #first guess
#         x_n = M

#         epsilon = 0.0001
#         for n in range(interations):

#             fn   = func(x_n,*args)
#             fnn  = (func((x_n + epsilon),*args) - func(x_n,*args))/epsilon


#             delta_n = -fn/fnn
#             x_n    += delta_n

#         print(x_n)
#         return x_n

















# def rho1_converge(rho1_init, rho01, P01, T01, gamma, cp, R):

#     def f(rho1, rho01, P01, T01, gamma, cp, R):
#         V1 = np.sqrt(2*cp*T01*(1-(rho1/rho01)**(gamma-1)))
#         T1 = T01 - V1**2/2/cp
#         P1 = P01*(T1/T01)**(gamma/(gamma-1)) # isentropic
#         rho1_new = P1/T1/R
#         print(rho1_new)
#         return (rho1 - rho1_new)

#     rho1 = fsolve(f,rho1_init,args=(rho01, P01, T01, gamma, cp, R))

#     V1 = np.sqrt(2*cp*T01*(1-(rho1/rho01)**(gamma-1)))
#     T1 = T01 - V1**2/2/cp
#     P1 = P01*(T1/T01)**(gamma/(gamma-1)) # isentropic

#     return rho1, V1, T1, P1;
