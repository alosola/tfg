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








# for kacker-okapuu tests


S = 0.008549314    # pitch [m]
alpha2 = np.radians(75.3)
alpha3 = np.radians(65)
c = 0.017007665 # chord
t = c/10 # max thickness
h = 0.011906101
bx = 0.6*h
Mw2 = 0.522
Mw3 = 0.959
P2 = 153781.890
P3 = 90483.172
gamma = 1.3
RHT  = 0.9
mu = 0.00001748 + 0.0000000431*1155.515021
D = h*4
v = 342.4038851
rho = 0.465082515

Rec = 6e+05 #rho*v*D/mu

t_o = 0.2

alpha2_vec = np.linspace(np.radians(70), np.radians(75.3), num=100)
loss = np.ones(100)
i = 0

for alpha2 in alpha2_vec:

    loss[i]  = kackerokapuu('rotor', S, alpha2, alpha3, t, c, bx, h, Mw2, Mw3, P2, P3, gamma, RHT, Rec, t_o)
    i +=1


import matplotlib.pyplot as plt

plt.plot(alpha2_vec,loss)



