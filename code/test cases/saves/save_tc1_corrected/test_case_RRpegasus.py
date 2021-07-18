# -*- coding: utf-8 -*-
"""
Created on Thu Jul 15 18:29:40 2021

@author: alondra sola

Project: turbine design tool (TFG)
Program: test case: Rolls Royce Pegasis 11 Mk.101
File: test_Case_RRpegasus.py
Description:

    Available data for a two-shaft gas turbine engine. The complete engine
    cycle is calculated to obtain the necessary parameters for input into the
    turbine design tool.

"""


eta_diffusor = 0.91             # efficiency of inlet diffusor
eta_fan_p    = 0.89             # efficiency of fan, primary flow
eta_fan_s    = 0.88             # efficiency of fan, secondary flow
eta_HPC      = 0.88             # efficiency of high pressure compressor
eta_cc       = 0.95             # efficiency in combustion chamber
eta_HPT      = 0.92             # efficiency of high pressure turbine
eta_LPT      = 0.90             # efficiency of low pressure turbine
eta_HPS      = 0.985            # mechanical efficiency of high pressure shaft
eta_LPS      = 0.99             # mechanical efficiency of low pressure shaft
DPcc         = 0.02             # pressure loss in combustion chamber

pi_fan_p     = 2                # compression ratio in fan, primary flow
pi_fan_s     = 2                # compression ratio in fan, secondary flow
pi_HPC       = 7.9              # compression ratio in high pressure compressor
bypass_ratio = 1.3              # ratio of secondary flow to primary flow

T04          = 1450             # total temperature at cc exit [K]
P04          = 1.43e+06         # total pressure at cc exit [Pa]

mdot_a       = 190              # total massflow [kg/s]
mdot_HTP     = 82.6             # massflow through HPT stage 1


deltaH_HPT1  = 16.8e+06         # enthalpy change from HPT stage 1

[21/06 10:59] Jorge Saavedra García
    Temperatura de salida de cámara de combustión 1450 K
​[21/06 11:00] Jorge Saavedra García
    Rendimiento aproximado de turbina de alta 92%, y de la de baja 90%
​[21/06 11:01] Jorge Saavedra García
    Rendimiento mecánico del eje de alta 98.5%
​[21/06 11:01] Jorge Saavedra García
    y del de baja 99%
​[21/06 11:01] Jorge Saavedra García
    En el eje de alta, la turbina de alta presión tiene 2 etapas, pero podemos centrarnos en la primera que es la que aporta el 60% de potencia al eje
​[21/06 11:03] Jorge Saavedra García

Por lo que esta turbina tendría un gasto de 82.6 kg/s, con una presión y temperatura total a la entrada de 1.43 MPa y 1450K
​[21/06 11:03] Jorge Saavedra García
    y esta turbina tiene que generar aproximadamente 16.8 MW
Edited