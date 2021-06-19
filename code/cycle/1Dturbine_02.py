import numpy as np
import matplotlib.pyplot as plt

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
    Mach30 = 0.421            # initial guess for exit Mach number (3 is rotor exit, 2 in notes) [-]
    GR0 = 0.32                # initial guess for degreee of reaction
    phi = 1.75                # initial guess for loading factor
    eta_stator = 0.8877       # initial guess for stator efficiency
    eta_rotor = 0.88          # initial guess for rotor efficiency


    oneD_stage(mdot, pi_turbine, P01, T01, P03, T03, gamma_t, cp_t, R, Mach30, GR0, phi, eta_stator, eta_rotor, DeltaH_turbine)


def oneD_stage(mdot, pi_turbine, P01, T01, P03, T03, gamma, cp, R, Mach3, GR, phi, eta_stator, eta_rotor, DeltaH):


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
    # alpha2 = np.radians(75)
    # V2x, U2 = oneD_rotor(alpha2, V2, T2, P2, P3, a2, phi, DeltaH, eta_rotor, cp, gamma, R)

    # cycling alpha 2 to find variation of DeltaH
    alpha2_range = np.linspace(np.radians(65),np.radians(75),num=5)
    for alpha2 in alpha2_range:
        V2x, U2 = oneD_rotor(alpha2, T01, V2, T2, P2, P3, a2, phi, DeltaH, eta_rotor, cp, gamma, R)





def oneD_rotor(alpha2, T01, V2, T2, P2, P3, a2, phi, DeltaH, eta_rotor, cp, gamma, R):

    # decompose velocity into axial and tangential
    V2x = V2*np.cos(alpha2)
    V2u = V2*np.sin(alpha2)

    # blade speed or peripheral speed
    U2 = np.sqrt(DeltaH/phi)
    U2 = np.sqrt(DeltaH/phi)

    # velocity triangle
    W2u = V2u - U2
    W2 = np.sqrt(W2u**2 + V2x**2)

    # relative inlet angle and relative Mach number
    beta2 = np.arctan(W2u/V2x)
    Mach2r = W2/a2

    # relative total quantities at rotor inlet
    T02r = T2*(1 + Mach2r**2*(gamma-1)/2)
    P02r = P2*(1 + (gamma-1)/2*Mach2r**2)**(gamma/(gamma-1))

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


    # cycle through beta iterations to find which meets DeltaH requirement
    beta3_range = -np.linspace(np.radians(60),np.radians(65),num=100)
    i = 0
    DeltaH_calc = np.zeros(len(beta3_range))

    for beta3 in beta3_range:
        DeltaH_calc[i], V3u, W3x = iter_beta3(beta3, W3, U2, V2u)
        i = i+1

    plt.plot(np.degrees(beta3_range), DeltaH_calc)
    plt.plot(np.degrees(beta3_range), DeltaH*np.ones_like(np.degrees(beta3_range)), 'k--')
    plt.plot(np.degrees(beta3_range), 1.01*DeltaH*np.ones_like(np.degrees(beta3_range)), 'b--')
    plt.plot(np.degrees(beta3_range), 0.99*DeltaH*np.ones_like(np.degrees(beta3_range)), 'b--')
    plt.ylabel('DeltaH [J/kg]')
    plt.xlabel('Beta3 [deg]')
    plt.title("Varying alpha2")

    # use a fixed falue of beta3
    # beta3 = np.radians(-65)
    # DeltaH_calc, V3u, W3x = iter_beta3(beta3, W3, U2, V2u)

    V3x = W3x
    V3 = np.sqrt(V3x**2 + V3u**2)

    # complete 1D analysis
    alpha3 = np.sin(V3u/V3)
    Mach3 = V3/a3

    # alternative thermodynamic analysis
    T03 = T3*(1+(gamma-1)/2*Mach3**2)
    P03 = P3*(1+(gamma-1)/2*Mach3**2)**(gamma/(gamma-1))
    DeltaH_thermo = cp*(T01-T03)

    return V2x, U2



def iter_beta3(beta3, W3, U2, V2u):

    # velocity projections, axial and tangential
    W3x = W3*np.cos(beta3)
    W3u = W3*np.sin(beta3)

    # constant radius, constant peripheral speed
    U3 = U2

    # tangential outlet speed
    V3u = W3u + U3

    # euler formulation to find work
    DeltaH = U2*(V2u - V3u)

    return DeltaH, V3u, W3x


def geometrical(mdot, rho2, V2x, U2, rho01, cp, gamma, T01):

    # ASSUMPTIONS
    RHT = 0.9            # ratio hub/tip radius
    rho10 = 0.98317      # initial guess for rho1

    # tip radius
    Rtip = np.sqrt(mdot/(np.pi*rho2*V2x*(1-RHT)))

    # hub radius
    Rhub = Rtip*RHT

    # blade height
    h2 = Rtip - Rhub

    # mean radius
    Rmean = Rtip - h2/2

    # angular speed rotor (< 500 m/s)
    Omega = U2/Rmean

    # area stage 2
    A2 = np.pi*(Rtip**2 - Rhub**2)

    # global flow speed
    V1 = np.sqrt(2*cp*T01*(1-(rho10/rho01)**(gamma-1)))

    # area stage 1
    A1 = mdot/(rho10*V1)

    # iterate untill A1 = A2 (cylindrical chamber)



def cycle_analysis(P01, T01, T03, pi_compressor, gamma_c, cp_c, gamma_t, cp_t, NPw, print_results):

    # ASSUMPTIONS:
    eta_compressor = 0.82  # compressor efficiency [-]
    eta_shaft = 0.94       # shaft efficiency [-]
    eta_turbine = 0.836     # turbine efficiency [-]
    eta_burner = 0.88      # burner efficiency [-]
    eta_combustor_p = 0.98 # compressor pressure efficiency [-]
    eta_trans = 0.89       # transmission efficiency [-]

    # DEFINE VALUES
    L = 44400000             # gasoline energy density [J/kg]


    ####################
    # COMPRESSOR STAGE #

    # the compression ratio provides the pressure after the compressor stage, ahead of the combustion chamber
    P02 = P01*pi_compressor

    # assuming an isentropic evolution, the temperature is easily found
    T02s = T01*pi_compressor**((gamma_c-1)/gamma_c)

    # the calculation of the total temperature takes into account the compressor efficiency
    T02 = T01 + (T02s - T01)/eta_compressor

    # the work required is founfrom the temperature increase
    DeltaH_compressor = cp_c*(T02 - T01)


    ###################
    # COMBUSTOR STAGE #

    # the pressure at the outlet can be assumed with the efficiency
    P03 = P02*eta_combustor_p

    # the exit temperature is a design constraint, and the temperature increase is calculated
#    DeltaT_combustor = T03 - T02

    # using the temperatures and energy density, the fuel-to-mass-flow ratio can be estimated
    f = (cp_c*T03 - cp_c*T02)/(eta_burner*L - cp_c*T03)


    ############################
    # TURBINE EXPLANSION STAGE #

    # the flow is expanded to atmospheric pressure (total)
    P04 = P01

    # the expansion ratio can be calculated immediately
    pi_turbine = P03/P04

    # assuming no losses. the expansion is isentropic
    T04s = T03*(1/pi_turbine)**((gamma_t - 1)/gamma_t)

    # incorporating turbine efficiency
    T04 = T03 - eta_turbine*(T03-T04s)


    ##################
    # WORK AND POWER #

    # the available work is
    DeltaH_turbine = cp_t*(T03-T04)

    # the work for transmission is the work not needed for compressor
    DeltaH_av_trans = DeltaH_turbine - DeltaH_compressor/eta_shaft/(1+f)

    # since the net power needed by the wheels is a design constraint, the massflow is calculated
    ma = NPw/(DeltaH_av_trans*eta_trans*(1+f))

    if print_results==True:
        print('Massflow: ', round(ma,4) ,' kg/s')
        print('T04: ', round(T04,2) , ' K')
        print('P03: ', round(P03,0) , ' Pa')
        print('Expansion ratio: ', round(pi_turbine,2) )
        print('Delta H compressor: ', round(DeltaH_compressor/1000, 2), ' kJ/kg')
        print('Delta H available for transmission: ', round(DeltaH_av_trans/1000, 2), ' kJ/kg')


    return ma, pi_turbine, DeltaH_turbine, P03, P04, T03, T04;





if __name__ == '__main__':
    main()
