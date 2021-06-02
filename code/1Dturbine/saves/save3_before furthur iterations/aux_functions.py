# -*- coding: utf-8 -*-
"""
Created on Thu May 27 19:15:23 2021

@author: alondra sola

Project: turbine design tool (TFG)
Program: 1D analysis (imitating the "automobile gas generator" design tool)
File: aux_functions.py
Description: contains definitions of thermodynamic functions used in the
             turbine model

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


def static_pressure(P0, gamma, Mach):
    P = P0/(1+(gamma-1)/2*Mach**2)**(gamma/(gamma-1))
    return P

def P2T(T0, P, P0, gamma):
    Ts = T0*(P/P0)**((gamma-1)/gamma)
    return Ts

def T2P(P, T0, T, gamma):
    P0 = P*(T0/T)**(gamma/(gamma-1))
    return P0

def stage_efficiency(eta, T0, Ts):
    T = T0 - eta*(T0-Ts)
    return T

def isen_velocity(cp, T0, T):
    V = np.sqrt(2*cp*(T0 - T))
    return V

def losses(V, Vs, P0x, P0y, Py):
    xi = (Vs**2 - V**2)/Vs**2
    lrt = V**2/Vs**2
    tpl = (P0x - P0y)/(P0x - Py)
    omega = (P0x - P0y)/(P0y - Py)
    return xi, lrt, tpl, omega

def static_density(P, T, R):
    rho = P/T/R
    return rho

def sonic(gamma, R, T, V):
    a = np.sqrt(gamma*R*T)
    M = V/a
    return a, M

def mag(x, y):
    return np.sqrt(x**2 + y**2)

def velocity_projections(mag, angle):
    x = mag*np.cos(angle)
    u = mag*np.sin(angle)
    return x ,u

def relative_temperature_pressure(T, P, Mach, gamma):
    T0 = T*(1+(gamma-1)/2*Mach**2)
    P0 = P*(1+(gamma-1)/2*Mach**2)**(gamma/(gamma-1))
    return T0, P0

def RPM(Omega):
    RPM = 60*Omega/2/np.pi
    return RPM






# def isen_evolution(x, M, gamma):
#     x0 = x*(1 + M**2*(gamma-1)/2)
#     return x0

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
