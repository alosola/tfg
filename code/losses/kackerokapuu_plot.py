# -*- coding: utf-8 -*-
"""
Created on Sat Jun 19 16:22:47 2021

@author: alondra sola

Project: turbine design tool (TFG)
Program: inclusion of losses for 1D turbine
File: kackerokapuu_plot.py
Descpription: simple depiction of surface plots fro kacker-okapuu correlations
(used to confirm that data has been imported well)

"""
## IMPORT NECESSARY LIBRARIES
import numpy as np                   # library for math and calculation functions
import matplotlib.pyplot as plt      # library for plots
import scipy.interpolate
   # greek letters

############ PLOT FOR FIGURE YP ALPHA IN = 0
# x = np.genfromtxt('fig1_x.csv', delimiter=',')               # import SC vector
# y = np.flip(np.genfromtxt('fig1_y.csv', delimiter=','))      # import alpha2 vector
# z = np.flip(np.genfromtxt('fig1_Zm.csv', delimiter=','), 1)  # import YP mesh

# YP_spline =  scipy.interpolate.RectBivariateSpline(x, y, z)  # create spline for evaluation


# N = 500
# X = np.linspace(np.amin(x),np.amax(x),N)
# Y = np.linspace(np.amin(y),np.amax(y),N)
# Z = np.zeros((N,N))

# for i in range(N):
#     for j in range(N):
#         Z[i,j] = YP_spline(X[i],Y[j])

# Xgrid = np.zeros((N,N))
# Ygrid = np.zeros((N,N))
# for i in range(N):
#     Xgrid[:,i] = X
#     Ygrid[i,:] = Y

# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.set_ylabel('Pitch-to-chord ratio s/c [-]')
# ax.set_xlabel(r'Inlet angle $\alpha$ [deg]')
# ax.set_zlabel(r'Pressure loss coefficient $Y_P$ [-]')
# # ax.plot_wireframe(Xgrid, Ygrid, Z, color='black')
# ax.plot_surface(Xgrid, Ygrid, Z, cstride=1, rstride=1, cmap='viridis', edgecolor='none')
# ax.set_title('Surface plot')
# plt.show()






############ PLOT FOR FIGURE YP ALPHA IN = ALPHA OUT
x = np.genfromtxt('fig2_x.csv', delimiter=',')               # import SC vector
y = np.genfromtxt('fig2_y.csv', delimiter=',')      # import alpha2 vector
z = np.flip(np.genfromtxt('fig2_Zm.csv', delimiter=','), 1)  # import YP mesh

YP_spline =  scipy.interpolate.RectBivariateSpline(x, y, z)  # create spline for evaluation


N = 200
X = np.linspace(np.amin(x),np.amax(x),N)
Y = np.linspace(np.amin(y),np.amax(y),N)
Z = np.zeros((N,N))

for i in range(N):
    for j in range(N):
        Z[i,j] = YP_spline(X[i],Y[j])

Xgrid = np.zeros((N,N))
Ygrid = np.zeros((N,N))
for i in range(N):
    Xgrid[:,i] = X
    Ygrid[i,:] = Y

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set_ylabel('Pitch-to-chord ratio s/c [-]')
ax.set_xlabel(r'Inlet angle $\alpha$ [deg]')
ax.set_zlabel(r'Pressure loss coefficient $Y_P$ [-]')
# ax.plot_wireframe(Xgrid, Ygrid, Z, color='black')
ax.plot_surface(Xgrid, Ygrid, Z, cstride=1, rstride=1, cmap='viridis', edgecolor='none')
ax.set_title('Surface plot')
plt.show()












######## STAGGER ANGLE

# x = np.genfromtxt('stg_x.csv', delimiter=',')               # import SC vector
# y = np.genfromtxt('stg_y.csv', delimiter=',')      # import alpha2 vector
# z = np.flip(np.genfromtxt('stg_Zm.csv', delimiter=','),0)  # import YP mesh

# YP_spline =  scipy.interpolate.RectBivariateSpline(x, y, z)  # create spline for evaluation

# N = 60
# X = np.linspace(np.amin(x),np.amax(x),N)
# Y = np.linspace(np.amin(y),np.amax(y),N)
# Z = np.zeros((N,N))

# for i in range(N):
#     for j in range(N):
#         Z[i,j] = YP_spline(X[i],Y[j])

# Xgrid = np.zeros((N,N))
# Ygrid = np.zeros((N,N))
# for i in range(N):
#     Xgrid[:,i] = X
#     Ygrid[i,:] = Y

# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.plot_surface(Xgrid, Ygrid, Z, cstride=1, rstride=1, cmap='viridis', edgecolor='none')
# ax.set_title('Surface plot')
# plt.show()