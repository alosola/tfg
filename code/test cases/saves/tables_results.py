# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 23:16:50 2021

@author: alondra sola

Project: turbine design tool (TFG)
Program: 1D analysis (imitating the "automobile gas generator" design tool)
File: tables_results.py
Description: contains functions for plotting results to cmd window in tables

"""

import numpy as np
import pandas as pd

def print_all_tables(one, two, thr, Mach3_init, alpha2_init, beta3_init, DeltaH_prod, DeltaH_calc, DeltaH_T, mdot, Deltabeta, psi, GR, h3h2, stator, rotor, eta_stator_init, eta_rotor_init):
    results_1(one, mdot)
    results_2(two, mdot, stator.tpl)
    results_3(thr, Mach3_init, rotor.tpl)
    global_parameters(one, two, thr, DeltaH_calc, DeltaH_T, Deltabeta, psi, GR, h3h2, eta_stator_init, eta_rotor_init, stator, rotor)
    convergence_parameters(one, two, thr, Mach3_init, alpha2_init, beta3_init, DeltaH_prod, DeltaH_calc)

def convergence_parameters(one, two, thr, Mach3_init, alpha2_init, beta3_init, DeltaH_prod, DeltaH_calc):
    " print table of convergence parameters to command window "

    print('CONVERGENCE PARAMETERS AND INITIAL GUESSES')

    data1 = {'Parameter':  ['Mach3 initial guess', 'Mach3 final value', 'Alpha2 initial guess', 'Alpha2 final value', 'Beta3 initial guess', 'Beta3 final value', 'DeltaH necessary', 'DeltaH calc'],
        'Value':     [Mach3_init, thr.vel.M, np.degrees(alpha2_init), np.degrees(two.alpha), np.degrees(beta3_init), np.degrees(thr.beta), DeltaH_prod/1000, DeltaH_calc/1000],
        'Unit':      ['-', '-',  'deg', 'deg',  'deg', 'deg', 'kJ/kg', 'kJ/kg']}
    df1 = pd.DataFrame (data1, columns = ['Parameter','Value','Unit'])
    print (df1 ,'\n')


def global_parameters(one, two, thr, DeltaH_calc, DeltaH_T, Deltabeta, psi, GR, h3h2, eta_stator_init, eta_rotor_init, stator, rotor):
    " print table of global parameters to command window "

    print('GLOBAL PARAMETERS')

    data1 = {'Parameter':  ['Delta H (euler)', 'Delta H (T0)','GR', 'Delta Beta', 'h3/h2','A3/A2','v2x/v1x','v3x/v2x','eta stator init','eta TT','eta rotor init','eta TS','w3/w2','loading factor','flow coefficient'],
             'Value':     [DeltaH_calc/1000,DeltaH_T/1000,GR,np.degrees(Deltabeta), h3h2, thr.geo.A/two.geo.A ,two.vel.Vx/one.vel.Vx, thr.vel.Vx/two.vel.Vx, eta_stator_init, stator.eta, eta_rotor_init, rotor.eta, thr.vel.W/two.vel.W, psi, thr.vel.Vx/two.vel.U],
             'Unit':      ['kJ/kg','kJ/kg','-','-','-','-','-','-','-','-','-','-','-','-','-']}
    df1 = pd.DataFrame (data1, columns = ['Parameter','Value','Unit'])
    print (df1 ,'\n')

def results_1(one, mdot):
    " print results from plane 1 to command window "

    print('FINAL VALUES PLANE 1')

    data1 = {'Variable':  ['T01', 'P01', 'V1', 'T1', 'P1', 'rho1', 'A1', 'Mach1', 'alpha1', 'mdot', 'R1tip', 'R1hub', 'h1', 'mdot'],
             'Value':     [one.T0, one.P0, one.vel.V, one.T, one.P, one.rho, one.geo.A, one.vel.M, np.degrees(one.alpha), mdot, one.geo.Rt, one.geo.Rh, one.geo.h, mdot],
             'Unit':      ['K', 'Pa', 'm/s', 'K', 'Pa', 'kg/m^3', 'm^2', '-', 'deg', 'kg/s', 'm', 'm', 'm', 'kg/s']}
    df1 = pd.DataFrame (data1, columns = ['Variable','Value','Unit'])
    print (df1 ,'\n')



def results_2(two, mdot, tpl_stator):
    " print results from plane 2 to command window "

    print('FINAL VALUES PLANE 2')

    data2 = {'Variable': ['T02', 'P02', 'V2', 'T2','P2','rho2','T2is','M2','Alpha2','V2x','V2u','W2','T02r','P02r','Beta2','W2x','W2u','Mw2','v2is','U2','TPL (stator)','A2','rh2','rt2','rm2','h2','omega'],
             'Value':    [two.T0, two.P0, two.vel.V, two.T, two.P,two.rho,two.Ts,two.vel.M,np.degrees(two.alpha),two.vel.Vx,two.vel.Vu,two.vel.W,two.T0r,two.P0r,np.degrees(two.beta),two.vel.Wx,two.vel.Wu,two.vel.Mr,two.vel.Vs,two.vel.U,tpl_stator,two.geo.A,two.geo.Rh,two.geo.Rt,two.geo.Rm,two.geo.h,two.vel.Omega],
             'Unit':     ['K', 'Pa', 'm/s', 'K', 'Pa', 'kg/m^3', 'K','-','deg','m/s^2', 'm/s^2','m/s^2','K', 'Pa', 'deg' ,'m/s^2','m/s^2','-','m/s^2','m/s^2','-','m^2','m','m','m','m','rad/s^2']}
    df2 = pd.DataFrame (data2, columns = ['Variable','Value','Unit'])
    print (df2 ,'\n')



def results_3(thr, Mach3_init,tpl_rotor):
    " print results from plane 3 to command window "

    print('FINAL VALUES PLANE 3')

    data3 = {'Variable': ['M3 impose','P3 impose','T03','P03','V3','T3','P3','rho3','T3is','T3iss','M3 achieve','P3 achieve','Alpha3','V3x','V3u','W3','T03r','P03r','Beta3','W3x','W3u','Mw3','U3','TPL (rotor)','A3','rh3','rt3','rm3'],
             'Value':    [Mach3_init,thr.P,thr.T0,thr.P0,thr.vel.V,thr.T,thr.P,thr.rho,thr.Ts,thr.Tss,thr.vel.M,thr.P,np.degrees(thr.alpha),thr.vel.Vx,thr.vel.Vu,thr.vel.W,thr.T0r,thr.P0r,np.degrees(thr.beta),thr.vel.Wx,thr.vel.Wu,thr.vel.Mr,thr.vel.U,tpl_rotor,thr.geo.A,thr.geo.Rh,thr.geo.Rt,thr.geo.Rm],
             'Unit':     ['-','Pa','K','Pa','m/s^2', 'K','Pa','kg/m^3','K','K','-','Pa','deg','m/s^2','m/s^2','m/s^2','K','Pa','deg','m/s^2','m/s^2','-','m/s^2','-','m^2','m','m','m']}
    df3 = pd.DataFrame (data3, columns = ['Variable','Value','Unit'])
    print (df3 ,'\n')



def print_limits(limits, values, checks):
    data = {'Variable': ['Turning (alpha)','Turning (beta)','Height ratio','Mach 2', 'Mach 2r','Mach 3','Mach3r','Beta 2'],
             'Min limit': limits[0,:],
             'Max limit': limits[1,:],
             'Value':     values[0,:],
             'Unit':     ['deg','deg','-','-','-','-','-','deg'],
             'Within limits?': checks[0,:]
             }
    df = pd.DataFrame (data, columns = ['Variable','Min limit','Max limit','Value','Unit','Within limits?'])
    print('LIMIT CHECKS')
    print (df ,'\n')


def print_mattignly(one,two,thr,inputs,printtype):
    datamat = {'Stage':    ['$T_0$','$T$','$P_0$','$P$','$M$','$v$','$w$','$alpha$','$beta$'],
            'One':      [round(one.T0,2), round(one.T,2), round(one.P0/1000,2), round(one.P/1000,2), round(one.vel.M,2), round(one.vel.V,2), '-', round(np.degrees(one.alpha),2), '-'],
            'Two':      [round(two.T0,2), round(two.T,2), round(two.P0/1000,2), round(two.P/1000,2), round(two.vel.M,2), round(two.vel.V,2), '-', round(np.degrees(two.alpha),2), '-'],
            'Two r':    [round(two.T0r,2), round(two.T,2), round(two.P0r/1000,2), round(two.P/1000,2), round(two.vel.Mr,2), '-', round(two.vel.W,2), '-', round(np.degrees(two.beta),2)],
            'Thr r':    [round(thr.T0r,2), round(thr.T,2), round(thr.P0r/1000,2), round(thr.P/1000,2), round(thr.vel.Mr,2), '-', round(thr.vel.W,2), '-', round(np.degrees(thr.beta),2)],
            'Thr':      [round(thr.T0,2), round(thr.T,2), round(thr.P0/1000,2), round(thr.P/1000,2), round(thr.vel.M,2), round(thr.vel.V,2), '-', round(thr.alpha,2), '-'],
            'Unit':     ['K','K','Pa','Pa','-','$\si{\meter\per\second}$','$\si{\meter\per\second}$','deg','deg'],
             }
    dfmat = pd.DataFrame(datamat)
    print('LIMIT CHECKS')
    if printtype=='terminal':
        print(dfmat.to_string(index=False))
    if printtype=='latex':
        print(dfmat.to_latex(index=False, escape=False))
    print ('\n')

def print_mattignly_terminal(one,two,thr,inputs):
    datamat = {'Stage':    ['$T_0$','$T$','$P_0$','$P$','$M$','$v$','$w$','$alpha$','$beta$'],
            'One':      [round(one.T0,2), round(one.T,2), round(one.P0/1000,2), round(one.P/1000,2), round(one.vel.M,2), round(one.vel.V,2), '-', round(np.degrees(one.alpha),2), '-'],
            'Two':      [round(two.T0,2), round(two.T,2), round(two.P0/1000,2), round(two.P/1000,2), round(two.vel.M,2), round(two.vel.V,2), '-', round(np.degrees(two.alpha),2), '-'],
            'Two r':    [round(two.T0r,2), round(two.T,2), round(two.P0r/1000,2), round(two.P/1000,2), round(two.vel.Mr,2), '-', round(two.vel.W,2), '-', round(np.degrees(two.beta),2)],
            'Thr r':    [round(thr.T0r,2), round(thr.T,2), round(thr.P0r/1000,2), round(thr.P/1000,2), round(thr.vel.Mr,2), '-', round(thr.vel.W,2), '-', round(np.degrees(thr.beta),2)],
            'Thr':      [round(thr.T0,2), round(thr.T,2), round(thr.P0/1000,2), round(thr.P/1000,2), round(thr.vel.M,2), round(thr.vel.V,2), '-', round(thr.alpha,2), '-'],
            'Unit':     ['K','K','Pa','Pa','-','$\si{\meter\per\second}$','$\si{\meter\per\second}$','deg','deg'],
             }
    dfmat = pd.DataFrame(datamat)
    print('LIMIT CHECKS')
    print(dfmat.to_string(index=False))
    print ('\n')


# def print_inputs():
#         datamat = {'Parameter':    ['$T_{01}$','$T_{03}$','$P_{01}$','$P_{03}$','$delta H$','$dot{m}$','$GR$','$psi$','$R_{ht}$'],
#             'One':      [round(one.T0,2), round(one.T,2), round(one.P0/1000,2), round(one.P/1000,2), round(one.vel.M,2), round(one.vel.V,2), '-', round(np.degrees(one.alpha),2), '-'],
#             'Two':      [round(two.T0,2), round(two.T,2), round(two.P0/1000,2), round(two.P/1000,2), round(two.vel.M,2), round(two.vel.V,2), '-', round(np.degrees(two.alpha),2), '-'],
#             'Two r':    [round(two.T0r,2), round(two.T,2), round(two.P0r/1000,2), round(two.P/1000,2), round(two.vel.Mr,2), '-', round(two.vel.W,2), '-', round(np.degrees(two.beta),2)],
#             'Thr r':    [round(thr.T0r,2), round(thr.T,2), round(thr.P0r/1000,2), round(thr.P/1000,2), round(thr.vel.Mr,2), '-', round(thr.vel.W,2), '-', round(np.degrees(thr.beta),2)],
#             'Thr':      [round(thr.T0,2), round(thr.T,2), round(thr.P0/1000,2), round(thr.P/1000,2), round(thr.vel.M,2), round(thr.vel.V,2), '-', round(thr.alpha,2), '-'],
#             'Unit':     ['K','K','Pa','Pa','-','$\si{\meter\per\second}$','$\si{\meter\per\second}$','deg','deg'],
#              }
#     dfmat = pd.DataFrame(datamat)
#     print('LIMIT CHECKS')
#     print(dfmat.to_string(index=False))
#     print(dfmat.to_latex(index=False, escape=False))
#     print ('\n')