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
import scipy.interpolate
from solve_functions import solve_rhox

def blade_geometry(mdot, rho, Vx, RHT, stage):

    Rtip = np.sqrt(mdot/(np.pi*rho*Vx*(1-RHT**2)))
    Rhub = Rtip*RHT
    h = Rtip-Rhub
    Rmean = Rtip - h/2
    Dmean = Rmean*2
    A = np.pi*(Rtip**2-Rhub**2)

    stage.geo.Rt = Rtip
    stage.geo.Rh = Rhub
    stage.geo.h = h
    stage.geo.Rm = Rmean
    stage.geo.Dm = Dmean
    stage.geo.A = A

    return


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

def losses(V, Vs, P0x, P0y, Py, component):

    if (P0y == Py):
        P0y = 1.00001*P0y
    if (P0x == Py):
        P0x = 1.00001*P0x
    component.xi = (Vs**2 - V**2)/Vs**2
    component.lrt = V**2/Vs**2
    component.tpl = (P0x - P0y)/(P0x - Py)
    component.omega = (P0x - P0y)/(P0y - Py) # pressure loss coefficient, ainly + mathiesen
    return component

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


def total_conditions(T3, V3, P3, cp, gamma):
    T03 = T3 + V3**2/2/cp
    P03 = P3*(T03/T3)**(gamma/(gamma-1))
    return T03, P03

def efficiency(I, F, Fs):
    eta = (I - F)/(I - Fs)
    return eta



def yes_or_no(question):
    while "the answer is invalid":
        reply = str(input(question+' (y/n): ')).lower().strip()
        if reply[:1] == 'y':
            return True
        if reply[:1] == 'n':
            return False



def stagger_angle(alpha2, alpha3):
    x = np.genfromtxt('dataset_stg_x.csv', delimiter=',')               # import SC vector
    y = np.genfromtxt('dataset_stg_y.csv', delimiter=',')      # import alpha2 vector
    z = np.flip(np.genfromtxt('dataset_stg_Zm.csv', delimiter=','),0)  # import YP mesh

    phi_spline = scipy.interpolate.RectBivariateSpline(x, y, z)  # create spline for evaluation
    phi = phi_spline.ev(abs(np.degrees(alpha2)), np.degrees(abs(alpha3)))

    return np.radians(phi)


def axial_chord(c, alpha2, alpha3):

    phi = stagger_angle(alpha2, alpha3)

    cx = c*np.cos(phi)

    return cx, phi


def zweiffel(beta2, alpha2, alpha3):

    alpha2 = np.degrees(alpha2)
    alpha3 = np.degrees(alpha3)

    x = np.genfromtxt('stg_x.csv', delimiter=',')
    y = np.genfromtxt('stg_y.csv', delimiter=',')
    z = np.flip(np.genfromtxt('stg_Zm.csv', delimiter=','),0)

    phi_spline = scipy.interpolate.RectBivariateSpline(x, y, z)  # create spline for evaluation
    phi = phi_spline.ev(alpha2, alpha3)


    sc0 = 0.427 + (90-alpha2)/58 - ((90-alpha2)/93)**2
    sc1 = 0.224 + (1.575 - (90-alpha2)/90)*(90-alpha2)/90

    xi = beta2/alpha2

    sc = sc0 + (sc1 + sc0)*xi**2*np.sign(xi)



def pitch(Z, angle_in, angle_out, V_in, V_out, P0_in, P0_out, h_in, h_out, c):

    stagger = stagger_angle(angle_in, angle_out)

    s_c = Z*np.cos(stagger)/2/np.cos(angle_out)**2/(V_in/V_out*np.tan(angle_in)+np.tan(abs(angle_out)))*h_in/h_out*P0_in/P0_out
    pitch = c*s_c

    return pitch



def reynolds(two, thr):

    dataset = np.genfromtxt('dataset_viscosity.csv', delimiter=',')
    x = dataset[:,0]
    y = dataset[:,1]

    spline = scipy.interpolate.InterpolatedUnivariateSpline(x,y)
    nu = np.float64(spline(np.array([two.T,thr.T])+273.15)) # data referred to celcius
    Re = np.array([two.vel.V,thr.vel.W])*np.array([two.geo.c,thr.geo.c])/nu
    two.Re = Re[0]
    thr.Re = Re[1]



def geometry(RHT, one, two, thr, mdot, R, gamma, cp):

    # def f_geometry(geometry, one, two, thr, mdot, R, gamma, cp):


    # blade geometry at 2
    blade_geometry(mdot, two.rho, two.vel.Vx, RHT, two)
    two.geo.c = two.geo.h/two.geo.hc
    two.geo.bx, two.geo.phi = axial_chord(two.geo.c, one.alpha, two.alpha)

    # rotational velocity
    two.vel.Omega = two.vel.U/two.geo.Rm
    two.vel.RPM = RPM(two.vel.Omega)

    # massflow and characteristics at inlet
    one.rho0 = static_density(one.P0, one.T0, R)                                  # = one.P0/one.T0/R
    one.rho, one.vel.V, one.geo.A = solve_rhox(0.9999*one.rho0, one.rho0, one.T0, mdot,  two.geo.A, gamma, cp)  # iteraciones hasta converger
    one.vel.Vx = one.vel.V
    one.T = one.T0 - one.vel.V**2/2/cp
    one.P = T2P(one.P0, one.T, one.T0, gamma)                     # = one.P0*(one.T/one.T0)**(gamma/(gamma-1)) # isentropic
    one.vel.a, one.vel.M = sonic(gamma, R, one.T, one.vel.V)

    # blade geometry at 1
    blade_geometry(mdot, one.rho, one.vel.V, RHT, one)

    # blade geometry at 3
    # what
    thr.geo.A = mdot/thr.vel.Vx/thr.rho
    thr.geo.h = thr.geo.A/np.pi/2/two.geo.Rm
    thr.geo.Rt = two.geo.Rm+thr.geo.h/2
    thr.geo.Dm = thr.geo.Rt*2 - thr.geo.h
    thr.geo.Rm = thr.geo.Dm/2
    thr.geo.c = thr.geo.h/thr.geo.hc
    thr.geo.bx, thr.geo.phi = axial_chord(thr.geo.c, thr.alpha, thr.alpha)
    thr.geo.Rh = thr.geo.Rt - thr.geo.h

    thr.vel.Omega = thr.vel.U/thr.geo.Rm
    thr.vel.RPM = RPM(thr.vel.Omega)

    # compute pitch with soderberg/zweiffel
    two.geo.s = pitch(0.6, one.alpha, two.alpha, one.vel.Vx, two.vel.Vx, one.P0, two.P0, one.geo.h, two.geo.h, two.geo.c)
    thr.geo.s = pitch(0.8, two.beta, thr.beta, two.vel.Wx, thr.vel.Wx, two.P0r, thr.P0r, two.geo.h, thr.geo.h, thr.geo.c)




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
