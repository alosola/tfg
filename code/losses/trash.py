# -*- coding: utf-8 -*-
"""
Created on Thu May 27 20:44:42 2021

@author: alond
"""

rhos = np.linspace(0,1,100)
As = []

def f(rho1, rho01, T01, mdot, A2, gamma, cp):
    V1 = np.sqrt(2*cp*T01*(1-(rho1/rho01)**(gamma-1)))
    A1 = mdot/rho1/V1

    return A1

i = 0
for rho1 in rhos:
    As.append(f(rho1, rho01, T01, mdot, A2, gamma, cp))

plt.plot(rhos,As,rhos,[A2]*len(rhos))