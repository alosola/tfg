# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 21:53:50 2021

@author: alond
"""

def soderberg_zweiffel():

    # SOODERBERG CORRELATION LOSSES
    # rotor
    thr.H = thr.T*cp
    thr.Hs = thr.Ts*cp
    epsr = 2*(thr.H - thr.Hs)/thr.vel.W**2
    #stator
    H2 = two.T*cp
    two.Hs = two.Ts*cp
    epse = 2*(H2 - two.Hs)/two.vel.V**2

    # efficiencies
    one.H0 = cp*one.T0
    thr.H0 = cp*thr.T0
    thr.Tss = one.T0*(thr.P/one.P0)**((gamma-1)/gamma)
    thr.Hs = cp*thr.Ts
    stator.eta = (one.H0 - thr.H0)/(one.H0 - (thr.Hs+thr.vel.V**2/2))
    rotor.eta = (one.H0 - thr.H0)/(one.H0 - thr.Hs)

    epsp = (two.vel.Vs**2 - two.vel.V**2)/two.vel.V**2
    epspp = (thr.vel.Ws**2 - thr.vel.W**2)/thr.vel.W**2

    knls_stator = epsp*(1+gamma*two.vel.M**2/2)
    knls_rotor = epspp*(1+gamma*thr.vel.Mr**2/2)