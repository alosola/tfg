# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 15:03:56 2021

@author: alondra sola

Project: turbine design tool (TFG)
Program: 1D analysis (imitating the "automobile gas generator" design tool)
File: converge_mach3.py
Description: contains the implementation of all the calculations within the
Mach3 convergence loop as well as the function with the non-linear solver

"""

import numpy as np                   # library for math and calculation functions
from scipy.optimize import fsolve
import aux_functions as f


def converge_mach3(Mach3_init, one, two, thr, stator, rotor, gamma, cp, R, GR, phi, DeltaH_prod):
    """
    Function to find the density rho for a certain point x, from the
    temperature, mass flux, and area data.
    Using an initial guess rhox_init, which is always less than rho0x
    """

    def f_mach3(Mach3, one, two, thr, stator, rotor, gamma, cp, R, GR, phi, DeltaH_prod):

        thr.vel.M = Mach3

        # using mach, find static pressure
        thr.P = f.static_pressure(thr.P0, gamma, thr.vel.M)    # = thr.P0/(1+(gamma-1)/2*thr.vel.M**2)**(gamma/(gamma-1))
        thr.Ts = f.P2T(thr.T0, thr.P, thr.P0, gamma)              # = thr.T0*(thr.P/thr.P0)**((gamma-1)/gamma)

        # find intermediate static pressure
        two.P = GR*(one.P0-thr.P)+thr.P # Siverding rp definition
        two.P = one.P0*((GR-1)*(1-(thr.P/one.P0)**((gamma-1)/gamma) ) +1)**(gamma/(gamma-1))  # GP rp definition
        ### where does that definition come from? doublecheck GR definition

        # isentropic evolution in stator: absolute total temperature is same
        two.T0 = one.T0
        two.Ts = f.P2T(one.T0, two.P, one.P0, gamma)             # =  one.T0*(two.P/one.P0)**((gamma-1)/gamma)
        two.vel.Vs = f.isen_velocity(cp, two.T0, two.Ts)                # = np.sqrt(2*cp*(two.T0 - two.Ts))

        # assume rotor efficiency
        two.T = f.stage_efficiency(stator.eta, two.T0, two.Ts)    # = two.T0 - stator.eta*(two.T0-two.Ts)
        two.vel.V = f.isen_velocity(cp, two.T0, two.T)                # = np.sqrt(2*cp*(two.T0 - two.T))


        # density (from static quantities)
        two.rho = f.static_density(two.P,two.T,R)                 # = two.P/two.T/R

        # speed of sound and Mach
        two.vel.a, two.vel.M = f.sonic(gamma, R, two.T, two.vel.V)            # = np.sqrt(gamma*R*two.T),    = two.vel.V/two.vel.a

        # check if pressure calculated meets the restriction
        two.T0b = two.T+two.vel.V**2/2/cp
        ### include this in the check loop

        two.P0 = f.T2P(two.P, two.T0, two.T, gamma)                  # = two.P*(two.T0/two.T)**(gamma/(gamma-1))

        # assume an two.alpha and project velocities
        two.vel.Vx, two.vel.Vu = f.velocity_projections(two.vel.V, two.alpha)

        # assume a loading factor and calculate peripheral speed
        two.vel.U = np.sqrt(DeltaH_prod/phi)

        # velocity triangle (pithagoras)
        two.vel.Wu = two.vel.Vu - two.vel.U
        two.vel.Wx = two.vel.Vx
        two.vel.W = f.mag(two.vel.Wu, two.vel.Wx)

        # relative inlet angle and relative Mach number
        two.beta = np.arctan(two.vel.Wu/two.vel.Vx)
        two.vel.Mr = two.vel.W/two.vel.a

        # relative total quantities at rotor inlet
        # (add relative speed to static conditions)
        # two.T0r = two.T*(1 + two.vel.Mr**2*(gamma-1)/2)              # total temperature is conserved
        # two.P0r = two.P*(1 + (gamma-1)/2*two.vel.Mr**2)**(gamma/(gamma-1))
        two.T0r, two.P0r = f.relative_temperature_pressure(two.T, two.P, two.vel.Mr, gamma)

        ############ ROTOR #############
        # in the rotor, the relative total temperature is constant (se conservan rotalpías)
        thr.T0r = two.T0r

        # ideal outlet temperature, real outlet temperature
        # assuming the expanse to thr.P, calculated with the assumed M3
        thr.Ts = f.P2T(thr.T0r, thr.P, two.P0r, gamma)                  # = thr.T0r*(thr.P/two.P0r)**((gamma-1)/gamma)
        thr.T = f.stage_efficiency(rotor.eta, thr.T0r, thr.Ts)       # = thr.T0r - rotor.eta*(thr.T0r - thr.Ts)

        # velocities
        thr.vel.W = f.isen_velocity(cp, thr.T0r, thr.T)                  # = np.sqrt(2*cp*(thr.T0r - thr.T))
        thr.vel.Ws = f.isen_velocity(cp, thr.T0r, thr.Ts)                # = np.sqrt(2*cp*(thr.T0r - thr.Ts))

        # density
        thr.rho = f.static_density(thr.P,thr.T,R)                    # = thr.P/thr.T/R

        # speed of sound and relative mach number
        thr.vel.a, thr.vel.Mr = f.sonic(gamma, R, thr.T, thr.vel.W)                      # = np.sqrt(gamma*R*thr.T)

        # total absolute pressure
        thr.P0r = f.T2P(thr.P, thr.T0r, thr.T, gamma)                           # = thr.P*(thr.T0r/thr.T)**(gamma/(gamma-1))

        # assume a value for thr.beta

        # velocity projections, axial and tangential
        thr.vel.Wx, thr.vel.Wu = f.velocity_projections(thr.vel.W, thr.beta)

        # constant radius, constant peripheral speed
        thr.vel.U = two.vel.U

        # tangential outlet speed
        thr.vel.Vu = thr.vel.Wu + thr.vel.U
        # V3u_euler = -DeltaH_prod/two.vel.U + two.vel.Vu

        # assume axial speed constant
        thr.vel.Vx = thr.vel.Wx
        thr.vel.V = f.mag(thr.vel.Vx, thr.vel.Vu)                           # = np.sqrt(thr.vel.Vx**2 + thr.vel.Vu**2)
        thr.alpha = np.arctan(thr.vel.Vu/thr.vel.Vx)
        Mach3_calculated = thr.vel.V/thr.vel.a

        return Mach3_calculated - Mach3

    thr.vel.M  = fsolve(f_mach3,Mach3_init,args=(one, two, thr, stator, rotor, gamma, cp, R, GR, phi, DeltaH_prod))[0]

    # using mach, find static pressure
    thr.P = f.static_pressure(thr.P0, gamma, thr.vel.M)    # = thr.P0/(1+(gamma-1)/2*thr.vel.M**2)**(gamma/(gamma-1))
    thr.Ts = f.P2T(thr.T0, thr.P, thr.P0, gamma)              # = thr.T0*(thr.P/thr.P0)**((gamma-1)/gamma)

    # find intermediate static pressure
    two.P = GR*(one.P0-thr.P)+thr.P # Siverding rp definition
    two.P = one.P0*((GR-1)*(1-(thr.P/one.P0)**((gamma-1)/gamma) ) +1)**(gamma/(gamma-1))  # GP rp definition
    ### where does that definition come from? doublecheck GR definition

    # isentropic evolution in stator: absolute total temperature is same
    two.T0 = one.T0
    two.Ts = f.P2T(one.T0, two.P, one.P0, gamma)             # =  one.T0*(two.P/one.P0)**((gamma-1)/gamma)
    two.vel.Vs = f.isen_velocity(cp, two.T0, two.Ts)                # = np.sqrt(2*cp*(two.T0 - two.Ts))

    # assume rotor efficiency
    two.T = f.stage_efficiency(stator.eta, two.T0, two.Ts)    # = two.T0 - stator.eta*(two.T0-two.Ts)
    two.vel.V = f.isen_velocity(cp, two.T0, two.T)                # = np.sqrt(2*cp*(two.T0 - two.T))


    # density (from static quantities)
    two.rho = f.static_density(two.P,two.T,R)                 # = two.P/two.T/R

    # speed of sound and Mach
    two.vel.a, two.vel.M = f.sonic(gamma, R, two.T, two.vel.V)            # = np.sqrt(gamma*R*two.T),    = two.vel.V/two.vel.a

    # check if pressure calculated meets the restriction
    two.T0b = two.T+two.vel.V**2/2/cp
    ### include this in the check loop

    two.P0 = f.T2P(two.P, two.T0, two.T, gamma)                  # = two.P*(two.T0/two.T)**(gamma/(gamma-1))

    # assume an two.alpha and project velocities
    two.vel.Vx, two.vel.Vu = f.velocity_projections(two.vel.V, two.alpha)

    # assume a loading factor and calculate peripheral speed
    two.vel.U = np.sqrt(DeltaH_prod/phi)

    # velocity triangle (pithagoras)
    two.vel.Wu = two.vel.Vu - two.vel.U
    two.vel.Wx = two.vel.Vx
    two.vel.W = f.mag(two.vel.Wu, two.vel.Wx)

    # relative inlet angle and relative Mach number
    two.beta = np.arctan(two.vel.Wu/two.vel.Vx)
    two.vel.Mr = two.vel.W/two.vel.a

    # relative total quantities at rotor inlet
    # (add relative speed to static conditions)
    # two.T0r = two.T*(1 + two.vel.Mr**2*(gamma-1)/2)              # total temperature is conserved
    # two.P0r = two.P*(1 + (gamma-1)/2*two.vel.Mr**2)**(gamma/(gamma-1))
    two.T0r, two.P0r = f.relative_temperature_pressure(two.T, two.P, two.vel.Mr, gamma)

    ############ ROTOR #############
    # in the rotor, the relative total temperature is constant (se conservan rotalpías)
    thr.T0r = two.T0r

    # ideal outlet temperature, real outlet temperature
    # assuming the expanse to thr.P, calculated with the assumed M3
    thr.Ts = f.P2T(thr.T0r, thr.P, two.P0r, gamma)                  # = thr.T0r*(thr.P/two.P0r)**((gamma-1)/gamma)
    thr.T = f.stage_efficiency(rotor.eta, thr.T0r, thr.Ts)       # = thr.T0r - rotor.eta*(thr.T0r - thr.Ts)

    # velocities
    thr.vel.W = f.isen_velocity(cp, thr.T0r, thr.T)                  # = np.sqrt(2*cp*(thr.T0r - thr.T))
    thr.vel.Ws = f.isen_velocity(cp, thr.T0r, thr.Ts)                # = np.sqrt(2*cp*(thr.T0r - thr.Ts))

    # density
    thr.rho = f.static_density(thr.P,thr.T,R)                    # = thr.P/thr.T/R

    # speed of sound and relative mach number
    thr.vel.a, thr.vel.Mr = f.sonic(gamma, R, thr.T, thr.vel.W)                      # = np.sqrt(gamma*R*thr.T)

    # total absolute pressure
    thr.P0r = f.T2P(thr.P, thr.T0r, thr.T, gamma)                           # = thr.P*(thr.T0r/thr.T)**(gamma/(gamma-1))

    # assume a value for thr.beta

    # velocity projections, axial and tangential
    thr.vel.Wx, thr.vel.Wu = f.velocity_projections(thr.vel.W, thr.beta)

    # constant radius, constant peripheral speed
    thr.vel.U = two.vel.U

    # tangential outlet speed
    thr.vel.Vu = thr.vel.Wu + thr.vel.U
    # V3u_euler = -DeltaH_prod/two.vel.U + two.vel.Vu

    # assume axial speed constant
    thr.vel.Vx = thr.vel.Wx
    thr.vel.V = f.mag(thr.vel.Vx, thr.vel.Vu)                           # = np.sqrt(thr.vel.Vx**2 + thr.vel.Vu**2)
    thr.alpha = np.arctan(thr.vel.Vu/thr.vel.Vx)

    return one, two, thr, stator, rotor
