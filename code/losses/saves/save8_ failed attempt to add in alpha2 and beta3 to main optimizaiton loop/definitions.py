# -*- coding: utf-8 -*-
"""
Created on Sun May 30 13:59:43 2021

@author: alondra sola

Project: turbine design tool (TFG)
Program: inclusion of losses for 1D turbine
File: definitinos.py
Description: containts definitions for additional clases and variable types

"""

class thermo():
    def __init__(self):
        self.gamma = None
        self.cp = None
        self.R = None

class plane():
    def __init__(self):
        self.P = None
        self.P0 = None
        self.P0r = None
        self.T = None
        self.T0 = None
        self.T0r = None
        self.Ts = None
        self.Tss = None
        self.H = None
        self.H0 = None
        self.Hs = None
        self.rho = None
        self.rho0 = None
        self.geo = geometry()
        self.vel = velocity()
        self.alpha = None
        self.beta = None
        self.Re = None

class component():
    def __init__(self):
        self.eta = None    # stage efficiency
        self.xi = None     # kinetic loss coefficient
        self.lrt = None    # loss ??
        self.tpl = None    # total pressure loss coefficient (calculated with pressures)
        self.omegaKC = None  # total pressure loss coefficient (calculated with Kacker-Okapuu)
        self.omega = None


class geometry():
    def __init__(self):
        self.Rh = None     # radius hub [m]
        self.Rt = None     # radius tip [m]
        self.h = None      # blade height [m]
        self.Rm = None     # mean radius [m]
        self.Dm = None     # mean diameter [m]
        self.A = None      # cross-sectional area [m^2]
        self.c = None      # chord [m]
        self.s = None      # pitch [m]
        self.t = None      # trailing edge thickness [m]
        self.hc = None     # height/chord ratio [-]
        self.phi = None    # stagger angle [º]
        self.cx = None     # axial chord [m]
        self.to = None     # trailing-edge-thickness to throat-opening ratio [-]

    # def mean_radius(self):
    #     self.mean = (self.Rh + self.Rt)/2

class velocity():
    def __init__(self):
        self.V = None
        self.Vs = None
        self.Vx = None
        self.Vu = None
        self.a = None
        self.M = None
        self.Mr = None
        self.U = None
        self.W = None
        self.Wx = None
        self.Wu = None
        self.Ws = None
        self.Omega = None
        self.RPM = None

class pressure_losses():
    def __init__(self):
        self.P = None        # profile losses
        self.S = None        # secondary losses
        self.TET = None      # trailing edge losses
        self.TC = None       # tip leakage losses
        self.Re = None       # Reynolds effect correction (for profile losses)

    def total_pressure_loss_coefficient(self):
        self.tot = self.P*self.Re + self.S + self.TET + self.TC




