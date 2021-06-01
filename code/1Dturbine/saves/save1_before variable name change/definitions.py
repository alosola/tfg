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
        self.P = geometry()
        self.rho1_init = None
        self.rho01 = 0
        self.T01 = 0
        self.mdot = 0
        self.A2 = 0
        self.gamma = 0
        self.cp = 0


class geometry():
    def __init__(self):
        self.Rheight = 0
        self.Rtip = 20