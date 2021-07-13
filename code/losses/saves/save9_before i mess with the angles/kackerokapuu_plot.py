# -*- coding: utf-8 -*-
"""
Created on Sat Jun 19 16:22:47 2021

@author: alondra sola

Project: turbine design tool (TFG)
Program: inclusion of losses for 1D turbine
File: kackerokapuu_plot.py
Descpription: simple depiction of surface plots for kacker-okapuu correlations
(used to confirm that data has been imported well)

"""
## IMPORT NECESSARY LIBRARIES
import numpy as np                   # library for math and calculation functions
import matplotlib.pyplot as plt      # library for plots
import scipy.interpolate


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
# ax.set_xlabel('Pitch-to-chord ratio s/c [-]')
# ax.set_ylabel(r'Inlet angle $\alpha$ [deg]')
# ax.set_zlabel(r'Pressure loss coefficient $Y_P$ [-]')
# # ax.plot_wireframe(Xgrid, Ygrid, Z, color='black')
# ax.plot_surface(Xgrid, Ygrid, Z, cstride=1, rstride=1, cmap='viridis', edgecolor='none')
# ax.set_title('Surface plot')
# plt.show()






############ PLOT FOR FIGURE YP ALPHA IN = ALPHA OUT
# x = np.genfromtxt('fig2_x.csv', delimiter=',')               # import SC vector
# y = np.genfromtxt('fig2_y.csv', delimiter=',')      # import alpha2 vector
# z = np.flip(np.genfromtxt('fig2_Zm.csv', delimiter=','), 1)  # import YP mesh

# YP_spline =  scipy.interpolate.RectBivariateSpline(x, y, z)  # create spline for evaluation


# N = 200
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
# ax.set_xlabel('Pitch-to-chord ratio s/c [-]')
# ax.set_ylabel(r'Inlet angle $\alpha$ [deg]')
# ax.set_zlabel(r'Pressure loss coefficient $Y_P$ [-]')
# # ax.plot_wireframe(Xgrid, Ygrid, Z, color='black')
# ax.plot_surface(Xgrid, Ygrid, Z, cstride=1, rstride=1, cmap='viridis', edgecolor='none')
# ax.set_title('Surface plot')
# plt.show()




# ######### THICKNESS TO CHORD RATIO
# def thickness_chord_ratio(sumtheta):

#     if sumtheta > 120:
#         t_c = 0.25
#     elif sumtheta < 40:
#         t_c = 0.15
#     else:
#         t_c = 0.15 + 1.25e-03*(sumtheta - 40)

#     return t_c


# N = 200
# X = np.linspace(5,150,N)
# Y = 0.15+1.25e-03*(X-40)

# for i in range(N):
#     Y[i] = thickness_chord_ratio(X[i])

# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.set_ylabel('Max thickness to chord ratio $t_{max}/c$ [-]')
# ax.set_xlabel(r'Sum of angles $\alpha_{in} + \alpha_{out}$ [deg]')
# # ax.plot_wireframe(Xgrid, Ygrid, Z, color='black')
# plt.plot(X, Y)
# ax.set_title('Surface plot')
# plt.show()







######### K1 factor
def losses_KP(X):
    if X < 0.2:
        K1 = 1
    else:
        K1 = 1 - 1.25*(X-0.2)

    return K1


N = 200
X = np.linspace(0,1,N)
K1 = 1 - 1.25*(X-0.2)
K2 = X**2

for i in range(N):
    K1[i] = losses_KP(X[i])

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel(r'Oulet Mach $M_{out}$ [-]')
ax.set_ylabel(r'Outlet Mach correction factor $K_1$ [-]')
# ax.plot_wireframe(Xgrid, Ygrid, Z, color='black')
plt.plot(X, K1)
plt.show()


fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.set_xlabel(r'Mach ratio $M_{in}/M_{out}$ [-]')
ax2.set_ylabel(r'Channel acceleration correction factor $K_2$ [-]')
# ax.plot_wireframe(Xgrid, Ygrid, Z, color='black')
plt.plot(X, K2)
plt.show()

# KP = 1-K2*(1-K1)
# fig3 = plt.figure()
# ax3 = fig3.add_subplot(111)
# ax3.set_xlabel(r'Mach ratio $M_{in}/M_{out}$ [-]')
# ax3.set_ylabel(r'Channel acceleration correction factor $K_2$ [-]')
# # ax.plot_wireframe(Xgrid, Ygrid, Z, color='black')
# plt.plot(X, KP)
# plt.show()










######### Mach 1 hub tip to hub ratio
# N = 200
# X = np.linspace(0.5,1,N)
# M_rotor = (1+5.2*np.absolute(X-1)**2.2)
# M_stator = (1+1.8*np.absolute(X-1)**2.2)

# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.set_xlabel(r'Tip-to-hub radius ratio $R_{th}$ [-]')
# ax.set_ylabel(r'$M_{1,hub}$ /$M_1$  [-]')
# plt.plot(X, M_rotor, label='Rotors')
# plt.plot(X, M_stator, label='Stators')
# plt.legend()
# plt.show()






# ############# REYNOLDS NUMBER CORRECTION
# def reynolds_correction(Rec):
#     # Rec = Reynolds number based on true chord and exit gas conditions

#     mean = 2*10**5

#     if Rec <= mean:
#         fRe = (Rec/mean)**-0.4
#     elif Rec < 10**6:
#         fRe = 1
#     else:
#         fRe = (Rec/10**6)**-0.2

#     return fRe


# N = 200
# X = np.linspace(4e+04,1e+08,N)
# fRev = np.linspace(1e+03,1e+07,N)

# for i in range(N):
#     fRev[i] = reynolds_correction(X[i])

# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.set_xlabel(r'Reynolds number $Re$ [-]')
# ax.set_ylabel(r'Correction factor $f_{Re}$  [-]')
# ax.set_xscale('log')
# plt.plot(X, fRev)
# plt.show()





# ############# ASPECT RATIO CORRECTION
def secondary_losses_fAR(h_c):

    if h_c > 2:
        fAR = 1/h_c
    else:
        fAR = (1-0.25*np.sqrt(2-h_c))/h_c

    return fAR

N = 200
X = np.linspace(0.2,4,N)
fAR = np.linspace(0.2,4,N)

for i in range(N):
    fAR[i] = secondary_losses_fAR(X[i])

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel(r'Aspect ratio $h/c$ [-]')
ax.set_ylabel(r'Aspect ratio correctionn $f_{AR}$  [-]')
plt.plot(X, fAR)
plt.show()





# ############# SUBSONIC MACH CORRECTION K3
# N = 200
# X = np.linspace(0,1.4,N)
# K3 = X**2


# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.set_xlabel(r'Inverse axial aspect ratio $c_x/h$  [-]')
# ax.set_ylabel(r'Correction factor $K_3$ [-]')
# plt.plot(X, K3)
# plt.show()







###################### IMPULSE BLADE ENERGY LOSS COEFFICIENT
# N = 200
# X = np.linspace(0,0.4,N)
# PHITET = np.linspace(0,0.4,N)


# # for axial entry blades (not impulse blades)
# dataset = np.genfromtxt('dataset_trailing_edge.csv', delimiter=',')
# x = dataset[:,0]
# y = dataset[:,1]
# z = 0.5*y

# spline = scipy.interpolate.InterpolatedUnivariateSpline(x,y)

# dphi = np.float64(spline(PHITET))


# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.set_xlabel(r'$t_{et}/o_t$  [-]')
# ax.set_ylabel(r'$\Delta\phi^2$  for $\alpha_{in} = 0$ [-]')
# plt.plot(x,y,label='axial entry nozzle')
# plt.plot(x,z,label='impulse blading')
# plt.legend()
# plt.show()





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