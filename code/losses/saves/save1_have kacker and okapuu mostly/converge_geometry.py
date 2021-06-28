# -*- coding: utf-8 -*-
"""
Created on Sun Jun 27 13:37:26 2021

@author: alondra sola

Project: turbine design tool (TFG)
Program: inclusion of losses for 1D turbine
File: converge_mach3.py
Description: contains functions for convergence of geometry with loss estimations

"""
import numpy as np
import aux_functions as f
from solve_functions import solve_rhox



def converge_geometry(RHT, one, two, thr, mdot, R, gamma, cp):

    # def f_geometry(geometry, one, two, thr, mdot, R, gamma, cp):


    # blade geometry at 2
    f.blade_geometry(mdot, two.rho, two.vel.Vx, RHT, two)

    # rotational velocity
    two.vel.Omega = two.vel.U/two.geo.Rm
    two.vel.RPM = f.RPM(two.vel.Omega)

    # massflow and characteristics at inlet
    one.rho0 = f.static_density(one.P0, one.T0, R)                                  # = one.P0/one.T0/R
    one.rho = solve_rhox(0.9999*one.rho0, one.rho0, one.T0, mdot,  two.geo.A, gamma, cp)  # iteraciones hasta converger
    one.vel.V = np.sqrt(2*cp*one.T0*(1-(one.rho/one.rho0)**(gamma-1)))
    one.vel.Vx = one.vel.V
    one.T = one.T0 - one.vel.V**2/2/cp
    one.P = f.T2P(one.P0, one.T, one.T0, gamma)                     # = one.P0*(one.T/one.T0)**(gamma/(gamma-1)) # isentropic
    one.vel.a, one.vel.M = f.sonic(gamma, R, one.T, one.vel.V)
    one.geo.A = mdot/one.rho/one.vel.V

    # blade geometry at 1
    f.blade_geometry(mdot, one.rho, one.vel.V, RHT, one)

    # blade geometry at 3
    # what
    thr.geo.A = mdot/thr.vel.Vx/thr.rho
    thr.geo.h = thr.geo.A/np.pi/2/two.geo.Rm
    thr.geo.Rt = two.geo.Rm+thr.geo.h/2
    thr.geo.Dm = thr.geo.Rt*2 - thr.geo.h
    thr.geo.Rm = thr.geo.Dm/2

    thr.vel.Omega = thr.vel.U/thr.geo.Rm


