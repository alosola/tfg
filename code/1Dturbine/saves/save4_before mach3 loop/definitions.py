# -*- coding: utf-8 -*-
"""
Created on Sun May 30 13:59:43 2021

@author: alondra sola

Project: turbine design tool (TFG)
Program: 1D analysis (imitating the "automobile gas generator" design tool)
File: definitinos.py
Description: containts definitions for additional clases and variable types

"""


class stage():
    def __init__(self):
        self.P = None
        self.P0 = None
        self.P0r = None
        self.T = None
        self.T0 = None
        self.T0r = None
        self.Ts = None
        self.rho = None
        self.rho0 = None
        self.A = None
        self.geo = geometry()
        self.vel = velocity()
        self.alpha = None
        self.beta = None

class component():
    def __init__(self):
        self.eta = None    # stage efficiency
        self.xi = None     # kinetic loss coefficient
        self.loss = None   # loss ??
        self.tpl = None    # total pressure loss
        self.omega = None  # rotational speed
        self.RPM = None    # revolutions per minute

class geometry():
    def __init__(self):
        self.Rh = None     # radius hub [m]
        self.Rt = None     # radius tip [m]
        self.h = None      # blade height [m]
        self.Rm = None     # mean radius [m]
        self.Dm = None     # mean diameter [m]
        self.A = None      # cross-sectional area [m^2]

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




