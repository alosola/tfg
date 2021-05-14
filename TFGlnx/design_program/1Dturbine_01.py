

import numpy as np
import matplotlib.pyplot as plt
from cycle_analysis import cycle_analysis

def main():


    ############################  DESIGN ANALYSIS ##############################
    # DEFINE DESIGN CONSTRAINTS

    pi_compressor = 4       # compressor pressure ratio [-]
    P01 = 101325            # inlet total pressure [Pa]
    T01 = 288               # inlet total temperature [K]
    NPw = 150000            # net power to wheels [W]
    T03 = 1400              # exit combustor temperature [K]


    # FIX FLUID PROPERTIES

    gamma_c = 1.4  # [-]
    cp_c = 1000    # [J/kg/K]
    gamma_t = 1.3  # [-]
    cp_t = 1240    # [J/kg/K]
    R = 287        # [J/kg/K]

    # COMPUTE REMAINING NECESSARY FLOW QUANTITIES
    mdot, pi_turbine, DeltaH_turbine, P03, P04, T03, T04 = cycle_analysis(P01, T01, T03, pi_compressor, gamma_c, cp_c, gamma_t, cp_t, NPw, False)

    ########################## END OF DESIGN ANALYSIS ##########################


    # VARIABLE REDEFINITION
    # 1 is inlet to stator (inlet to turbine)
    # 2 is between stator and rotor
    # 3 is outlet to rotor (outlet to turbine)

    # engine variables will be overwritten by turbine variables
    P01 = P03
    T01 = T03
    P03 = P04
    T03 = T04


    # CONSTRAINTS (from cycle analysis):
    # mass flow (mdot), expansion ratio (pi_turbine), required specific work (DeltaH_turbine)

    # SELECT INITIAL VALUES
    # note: 3 is rotor exit, 2 in my report
    Mach30 = 0.42             # initial guess for exit Mach number (3 is rotor exit, 2 in notes) [-]
    GR0 = 0.32                # initial guess for degreee of reaction
    phi = 1.75                # initial guess for loading factor
    eta_stator = 0.98
    eta_rotor = 0.98

def oneD_stage(mdot, pi_turbine, Wreq, P01, T01, P03, T03, gamma, cp, R, Mach3, GR, phi, eta_stator, eta_rotor, DeltaH):

    # using mach, find static pressure
    P3 = P03/(1+(gamma-1)/2*Mach3**2)**(gamma/(gamma-1))

    # find intermediate static pressure
    P2 = P01*( 1 - (1 - GR)*(1 - (P3/P01)**((gamma-1)/gamma) ) )**(gamma/(gamma-1))

    # isentropic evolution in stator: absolute total temperature is same
    T02 = T01
    T2s = T02*(P2/P01)**((gamma-1)/gamma)
    T2 = T02 - eta_stator*(T02-T2s)
    P02 = P2*(T02/T2)**(gamma/(gamma-1))

    # velocities
    V2 = np.sqrt(2*cp*(T02 - T2))
    V2s = np.sqrt(2*cp*(T02 - T2s))

    # density (from static quantities)
    rho2 = P2/T2/R

    # speed of sound and Mach
    a2 = np.sqrt(gamma*R*T2)
    Mach2 = V2/a2

    # stator kinetic loss coefficient
    xi_stator = (V2s**2 - V2**2)/V2**2

    # stator pressure loss coefficient
    omega_stator = (P01 - P02)/(P02 - P2)

    # rotor convergence
    alpha2 = 0.5

    oneD_rotor(alpha2, V2, P3, a2, phi, DeltaH, eta_rotor, cp, gamma, R)





def oneD_rotor(alpha2, beta3, V2, P3, a2, phi, DeltaH, eta_rotor, cp, gamma, R):

    # decompose velocity into axial and tangential
    V2x = V2*np.cos(alpha2)
    V2u = V2*np.sin(alpha2)

    # blade speed or peripheral speed
    U2 = np.sqrt(DeltaH/phi)

    # velocity triangle
    W2u = V2u - U2
    W2 = np.sqrt(W2u**2 + V2x**2)

    # relative inlet angle and relative Mach number
    beta2 = np.arctan(W2u/V2x)
    Mach2r = W2/a2

    # relative total quantities at rotor inlet
    T02r = 1 + (gamma-1)/2*Mach2r**2
    P02r = (1 + (gamma-1)/2*Mach2r**2)**(gamma/(gamma-1))

    # in the rotor, the relative total temperature is constant (se conservan rotalpías)
    T03r = T02r

    # ideal outlet temperature, real outlet temperature
    T3s = T03r*(P3/P02r)**((gamma-1)/gamma)
    T3 = T03r - eta_rotor*(T03r - T3s)

    # velocities
    W3 = np.sqrt(2*cp*(T03r - T3))
    W3s = np.sqrt(2*cp*(T03r - T3s))

    # density
    rho3 = P3/T3/R

    # speed of sound and relative mach number
    a3 = np.sqrt(gamma*R*T3)
    Mach3r = W3/a3

    # total absolute pressure
    P03r = P3*(T03r/T3)**(gamma/(gamma-1))

    # rotor kinetic loss
    xi_rotor = (W3s**2 - W3**2)/W3**2

    # rotor pressure loss coefficient
    omega_rotor = (P02r - P03r)/(P03r - P3)



    beta3_range = range(np.radians(60),np.radians(65),100)
    i = 0

    DeltaH_calc = np.zeros(100)

    for beta3 in beta3_range:

        DeltaH_calc[i] = iter_beta3(beta3, W3, U2, V2u)

        i = i+1

    plt.plot(np.degrees(beta3_range), DeltaH_calc)
    plt.plot(np.degrees(beta3_range), DeltaH*np.ones_like(np.degrees(beta3_range)), 'k--')
    plt.ylabel('DeltaH')
    plt.xlabel('Beta3 [degrees]')
    plt.title('Work variation with beta3')


def iter_beta3(beta3, W3, U2, V2u):

    # velocity projections, axial and tangential
    W3x = W3*np.cos(beta3)
    W3u = W3*np.sin(beta3)

    # constant radius, constant peripheral speed
    U3 = U2

    # tangential outlet speed
    V3u = W3u - U3

    # euler formulation to find work
#    DeltaH_calc = U2*(V2u - V3u)

    return U2*(V2u - V3u)












if __name__ == '__main__':
    main()
