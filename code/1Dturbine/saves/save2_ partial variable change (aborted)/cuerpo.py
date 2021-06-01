
# using mach, find static pressure
thr.P = thr.P0/(1+(gamma-1)/2*thr.vel.M**2)**(gamma/(gamma-1))
thr.Ts = T03*(thr.P/thr.P0)**((gamma-1)/gamma)

# find intermediate static pressure
two.P = GR*(one.P0-thr.P)+thr.P # Siverding rp definition
two.P = one.P0*((GR-1)*(1-(thr.P/one.P0)**((gamma-1)/gamma) ) +1)**(gamma/(gamma-1))  # GP rp definition
### where does that definition come from? doublecheck GR definition

# isentropic evolution in stator: absolute total temperature is same
two.T0 = one.T0
two.Ts = one.T0*(two.P/one.P0)**((gamma-1)/gamma)
two.vel.Vs = np.sqrt(2*cp*(two.T0 - two.Ts))

# assume rotor efficiency
two.T = two.T0 - eta_stator*(two.T0-two.Ts)
two.vel.V = np.sqrt(2*cp*(two.T0 - two.T))

# stator kinetic loss coefficient
xi_stator = (two.vel.Vs**2 - two.vel.V**2)/two.vel.Vs**2
### equation is different in report and excel sheet (corregir)

# check loss inpose ??? ###
loss_stator = two.vel.V**2/two.vel.Vs**2

# density (from static quantities)
rho2 = two.P/two.T/R

# speed of sound and Mach
a2 = np.sqrt(gamma*R*two.T)
Mach2 = two.vel.V/a2

# check if pressure calculated meets the restriction
two.T0b = two.T+two.vel.V**2/2/cp
### include this in the check loop

P02 = two.P*(two.T0/two.T)**(gamma/(gamma-1))

# Total pressure loss
tpl_stator = (one.P0 - P02)/(one.P0 - two.P)

# stator pressure loss coefficient, Total pressure loss down 02
omega_stator = (one.P0 - P02)/(P02 - two.P)

# assume an alpha2 and project velocities
two.vel.Vx = two.vel.V*np.cos(alpha2)
two.vel.Vu = two.vel.V*np.sin(alpha2)

# assume a loading factor and calculate peripheral speed
U2 = np.sqrt(DeltaH_prod/phi)

# velocity triangle (pithagoras)
W2u = two.vel.Vu - U2
W2x = two.vel.Vx
W2 = np.sqrt(W2u**2 + W2x**2)

# relative inlet angle and relative Mach number
beta2 = np.arctan(W2u/two.vel.Vx)
Mach2r = W2/a2

# relative total quantities at rotor inlet
# (add relative speed to static conditions)
two.T0r = two.T*(1 + Mach2r**2*(gamma-1)/2) # total temperature is conserved
P02r = two.P*(1 + (gamma-1)/2*Mach2r**2)**(gamma/(gamma-1))


############ ROTOR #############
# in the rotor, the relative total temperature is constant (se conservan rotalpías)
T03r = two.T0r

# ideal outlet temperature, real outlet temperature
# assuming the expanse to thr.P, calculated with the assumed M3
thr.Ts = T03r*(thr.P/P02r)**((gamma-1)/gamma)
thr.T = T03r - eta_rotor*(T03r - thr.Ts)

# velocities
W3 = np.sqrt(2*cp*(T03r - thr.T))
W3s = np.sqrt(2*cp*(T03r - thr.Ts))

# rotor kinetic loss
xi_rotor = (W3s**2 - W3**2)/W3s**2
### denominator no coincide entre excel and JS report

# check loss inpose ??? ###
loss_rotor = W3**2/W3s**2

# density
rho3 = thr.P/thr.T/R

# speed of sound and relative mach number
a3 = np.sqrt(gamma*R*thr.T)
thr.vel.Mr = W3/a3

# total absolute pressure
thr.P0r = thr.P*(T03r/thr.T)**(gamma/(gamma-1))

# Total pressure loss
tpl_rotor = (P02r - thr.P0r)/(P02r - thr.P)
# equation where does it come from?

# stator pressure loss coefficient, Total pressure loss down 02
omega_rotor = (P02r - thr.P0r)/(thr.P0r - thr.P)

# assume a value for beta3

# velocity projections, axial and tangential
W3x = W3*np.cos(beta3)
W3u = W3*np.sin(beta3)

# constant radius, constant peripheral speed
U3 = U2

# tangential outlet speed
V3u = W3u + U3
V3u_euler = -DeltaH_prod/U2 + two.vel.Vu

# assume axial speed constant
V3x = W3x
V3 = np.sqrt(V3x**2 + V3u**2)
alpha3 = np.arctan(V3u/V3x)
thr.vel.Mc = V3/a3
DeltaH_T1 = U3*(two.vel.Vu-V3u)    # euler formulation

# compute total condtions
T03 = thr.T*(1+(gamma-1)/2*thr.vel.M**2)
thr.P0 = thr.P*(1+(gamma-1)/2*thr.vel.M**2)**(gamma/(gamma-1))
DeltaH_two.T = cp*(one.T0-T03)






# blade geometry at 2
R2tip, R2hub, h2, R2mean, D2mean, A2 = blade_geometry(mdot, rho2, two.vel.Vx, RHT)

# rotational velocity
Omega2 = U2/R2mean
RPM = 60*Omega2/2/np.pi

# massflow and characteristics at inlet
rho01 = one.P0/one.T0/R
rho1 = rhox_converge(0.9999*rho01, rho01, one.T0, mdot, A2, gamma, cp) # iteraciones hasta converger
V1 = np.sqrt(2*cp*one.T0*(1-(rho1/rho01)**(gamma-1)))
T1 = one.T0 - V1**2/2/cp
P1 = one.P0*(T1/one.T0)**(gamma/(gamma-1)) # isentropic
a1 = np.sqrt(T1*gamma*R)
Mach1 = V1/a1
A1 = mdot/rho1/V1

# blade geometry at 1
R1tip, R1hub, h1, R1mean, D1mean, A1 = blade_geometry(mdot, rho1, V1, RHT)

# blade geometry at 3
########## NO ENTIENDO NADAAAAA
A3 = mdot/V3x/rho3
h3 = A3/np.pi/2/R2mean
R3tip = R2mean+h3/2
D3mean = R3tip*2 - h3
R3mean = D3mean/2

Omega3 = U3/R3mean

GR_enthalpy = (two.T - thr.T)/(two.T0 - T03) # Assume v1 = v3, not right
Deltabeta = abs(beta2 - beta3) # por qué no negativo?
DeltaW = W3 - W2
heightratio = h3/h2
