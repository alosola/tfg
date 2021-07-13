# -*- coding: utf-8 -*-
"""
Created on Sat Jul 10 21:15:19 2021

@author: alondra sola

Project: turbine design tool (TFG)
Program: inclusion of losses for 1D turbine
File: denton_plot.py
Descpription: simple depiction of surface plots for denton correlations
(used to confirm that data has been imported well)
"""

########### PLOT FOR CONTOUR OR SURFACE FIGURE OF UNSHROUDED BLADE DATA
x = np.genfromtxt('unshrouded.x.csv', delimiter=',')               # import SC vector
y = np.flip(np.genfromtxt('unshrouded.y.csv', delimiter=','))      # import alpha2 vector
z = np.flip(np.genfromtxt('unshrouded.Zm.csv', delimiter=','), 0)  # import YP mesh

# YP_spline =  scipy.interpolate.RectBivariateSpline(x, y, z)  # create spline for evaluation


# N = 20
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

fig = plt.figure()
ax = fig.add_subplot(1, 2, 1)
ax.set_xlabel('Pitch-to-chord ratio s/c [-]')
ax.set_ylabel(r'Inlet angle $\alpha$ [deg]')
# ax.set_zlabel(r'Pressure loss coefficient $Y_P$ [-]')
# ax.plot_wireframe(Xgrid, Ygrid, Z, color='black')
# ax.plot_surface(Xgrid, Ygrid, Z, cstride=1, rstride=1, cmap='viridis', edgecolor='none')
plt.contour(x,y,z,12)
ax.set_title('Unshrouded blade')
plt.show()