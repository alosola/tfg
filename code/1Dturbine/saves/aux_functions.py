# -*- coding: utf-8 -*-
"""
Created on Thu May 27 19:15:23 2021

@author: alond
"""

import numpy as np                   # library for math and calculation functions
from scipy.optimize import newton, fsolve


def rho1_converge(rho1_init, rho01, T01, mdot, A2, gamma, cp):

    def f(rho1, rho01, T01, mdot, A2, gamma, cp):
        V1 = np.sqrt(2*cp*T01*(1-(rho1/rho01)**(gamma-1)))
        A1 = mdot/rho1/V1
        return A2-A1

    rho1 = fsolve(f,rho1_init,args=(rho01, T01, mdot, A2, gamma, cp))[0]
    # rho1 = Newton_Raphson(f, rho1_init, rho01, T01, mdot, A2, gamma, cp)

    return rho1;



def blade_geometry(mdot, rho, Vx, RHT):
    Rtip = np.sqrt(mdot/(np.pi*rho*Vx*(1-RHT**2)))
    Rhub = Rtip*RHT
    h = Rtip-Rhub
    Rmean = Rtip - h/2
    Dmean = Rmean*2
    A = np.pi*(Rtip**2-Rhub**2)

    return Rtip, Rhub, h, Rmean, Dmean, A


class turbine():

    def __init__(self):
        self.P = geometry()
        self.rho1_init = None
        self.rho01 = 0
        self.T01 = 0
        self.mdot = 0
        self.A2 = 0
        self.gamma = 0
        self.cp = 0

    # def f(self):
    #     V1 = np.sqrt(2*self.cp*self.T01*(1-(self.rho1/self.rho01)**(self.gamma-1)))
    #     A1 = self.mdot/self.rho1/V1
    #     return self.A2-A1

    # def rho1_converge(self):
    #     args=(self.rho01, self.T01, self.mdot, self.A2, self.gamma, self.cp)
    #     self.rho1 = fsolve(self.f)[0]


class geometry():
    def __init__(self):
        self.Rheight = 0
        self.Rtip = 20

turb = turbine()
turb.P.Rtip = 20
turb.rho1_init = 1.1
turb.rho01 = 1.2
turb.T01 = 1000
turb.mdot = 7
turb.A2 = 0.05
turb.gamma = 1.4
turb.cp = 1000

print(turb.P.Rtip)


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
