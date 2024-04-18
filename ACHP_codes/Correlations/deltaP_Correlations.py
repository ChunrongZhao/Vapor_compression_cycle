from __future__         import division, print_function, absolute_import
from math               import pi,log,sqrt,exp,cos,sin,tan,log10
from scipy.integrate    import quad,quadrature,trapz,simps,fixed_quad
from scipy.optimize     import brentq,fsolve
import numpy            as np
import CoolProp         as CP

# Machine precision
machine_eps = np.finfo(float).eps
print(machine_eps, '; num=2')


# ------------------------------------------------------------------
# accelerational pressure drop
def AccelPressureDrop(x_min, x_max, AS, G, T_bubble, T_dew, D=None, rho_sat_L=None, rho_sat_V=None, slipModel='Zivi'):
    """
    Accelerational pressure drop

    From -dpdz|A=G^2*d[x^2v_g/alpha+(1-x)^2*v_f/(1-alpha)]/dz

    Integrating over z from 0 to L where x=x_1 at z=0 and x=x_2 at z=L

    Maxima code:
        alpha:1/(1+S*rho_g/rho_f*(1-x)/x)$
        num1:x^2/rho_g$
        num2:(1-x)^2/rho_f$
        subst(num1/alpha+num2/(1-alpha),x,1);
        subst(num1/alpha+num2/(1-alpha),x,0);
    """
    if rho_sat_L == None or rho_sat_V == None:
        AS.update(CP.QT_INPUTS, 1.0, T_dew)
        rho_sat_V               = AS.rhomass()      # [kg/m^3]
        AS.update(CP.QT_INPUTS, 0.0, T_bubble)
        rho_sat_L               = AS.rhomass()      # [kg/m^3]

    def f(x, AS, G, D, T_bubble, T_dew, rho_L, rho_V):
        if abs(x) < 1e-12:
            return 1 / rho_sat_L
        elif abs(1 - x) < 1e-12:
            return 1 / rho_sat_V
        else:
            if slipModel == 'Premoli':
                S               = Premoli(x, AS, G, D, T_bubble, T_dew, rho_L, rho_V)
            elif slipModel == 'Zivi':
                S               = pow(rho_L / rho_V, 1/3)
            elif slipModel == 'Homogeneous':
                S               = 1
            else:
                raise ValueError("slipModel must be either 'Premoli', 'Zivi' or 'Homogeneous'")

            alpha               = 1 / (1 + S * rho_V / rho_L * (1 - x) / x)

            return x**2 / rho_V / alpha + (1 - x)**2 / rho_L / (1 - alpha)

    return G**2 * (f(x_min, AS, G, D, T_bubble, T_dew, rho_sat_L, rho_sat_V) - f(x_max, AS, G, D, T_bubble, T_dew, rho_sat_L, rho_sat_V))


def Premoli(x, AS, G, D, T_bubble, T_dew, rho_L=None, rho_V=None):
    """
    return Premoli (1970) slip flow factor: Void fraction
    function copied from ACMODEL source code
    same correlations can be found in the Appendix A2 of Petterson (2000)
    """
    if rho_L == None or rho_V == None:
        AS.update(CP.QT_INPUTS, 1.0, T_dew)
        rho_V                   = AS.rhomass()      # [kg/m^3]
        AS.update(CP.QT_INPUTS, 0.0, T_bubble)
        rho_L                   = AS.rhomass()      # [kg/m^3]

    AS.update(CP.QT_INPUTS, 0.0, T_bubble)
    mu_L                        = AS.viscosity()    # [Pa-s]
    p_sat                       = AS.p()            # [Pa]
    AS.update(CP.PQ_INPUTS, p_sat, x)
    sigma                       = AS.surface_tension() # [N/m]

    PI1                         = rho_V / rho_L
    We                          = pow(G, 2) * D / (sigma * rho_L)
    Re_L                        = G * D / mu_L
    F_1                         = 1.578 * pow(Re_L, -0.19) * pow(PI1, -0.22)
    F_2                         = 0.0273 * We * pow(Re_L, -0.51) * pow(PI1, 0.08)
    Y                           = (x / (1 - x)) * 1 / PI1
    S                           = 1 + F_1 * pow((Y / (1 + F_2 * Y) - F_2 * Y), 0.5)

    return S


# averaged pressure drop over [x_min, x_max] or [x=0, x=L]
def LMPressureGradientAvg(x_min, x_max, AS, G, D, T_bubble, T_dew, C=None, satTransport=None):
    """
    Returns the average pressure gradient (or average pressure drop) between qualities of x_min and x_max.

    To obtain the pressure gradient for a given value of x, pass it in as x_min and x_max

    Required parameters:
    * x_min : The minimum quality for the range [-]
    * x_max : The maximum quality for the range [-]
    * AS : AbstractState with the refrigerant name and backend
    * G : Mass flux [kg/m^2/s] rho * u
    * D : Diameter of tube [m]
    * T_bubble : Bubble point temperature of refrigerant [K]
    * T_dew : Dewpoint temperature of refrigerant [K]

    Optional parameters:
    * C : The coefficient in the pressure drop
    * satTransport : A dictionary with the keys 'mu_f','mu_g,'v_f','v_g' for the saturation properties.
                     So they can be calculated once and passed in for a slight improvement in efficiency
    """
    def LMFunc(x):
        dp_dz, alpha     = LockhartMartinelli(AS, G, D, x, T_bubble, T_dew, C, satTransport)
        return dp_dz

    # Use Simpson's Rule to calculate the average pressure gradient
    # Can't use adaptive quadrature since function is not sufficiently smooth
    # Not clear why not sufficiently smooth at x > 0.9

    if x_min == x_max:
        return LMFunc(x_min)
    else:
        # Calculate the transport properties once
        satTransport                    = {}
        AS.update(CP.QT_INPUTS, 0.0, T_bubble)
        satTransport['v_f']             = 1 / AS.rhomass()  # [m^3/kg]
        satTransport['mu_f']            = AS.viscosity()    # [Pa-s]
        AS.update(CP.QT_INPUTS, 1.0, T_dew)
        satTransport['v_g']             = 1 / AS.rhomass()  # [m^3/kg]
        satTransport['mu_g']            = AS.viscosity()    # [Pa-s]

        xx                              = np.linspace(x_min, x_max, 30)
        DP                              = np.zeros_like(xx)
        for i in range(len(xx)):
            DP[i]                       = LMFunc(xx[i])
        return -simps(DP, xx) / (x_max - x_min)


def LockhartMartinelli(AS, G, D, x, T_bubble, T_dew, C=None, satTransport=None):
    # Following the method laid out in ME506 notes on
    # Separated Flow pressure drop calculations

    # Convert the quality, which might come in as a single numpy float value, to a float
    # With the conversion, >20x speedup in the LockhartMartinelli function, not clear why
    x = float(x)

    if satTransport == None:
        # Calculate Necessary saturation properties
        AS.update(CP.QT_INPUTS, 0.0, T_bubble)
        v_f                 = 1 / AS.rhomass()  # [m^3/kg]
        mu_f                = AS.viscosity()    # [Pa-s]
        AS.update(CP.QT_INPUTS, 1.0, T_dew)
        v_g                 = 1 / AS.rhomass()  # [m^3/kg]
        mu_g                = AS.viscosity()    # [Pa-s]
    else:
        # Pull out of the dictionary
        v_f                 = satTransport['v_f']
        v_g                 = satTransport['v_g']
        mu_f                = satTransport['mu_f']
        mu_g                = satTransport['mu_g']

    # 1. Find the Reynolds Number for each phase based on the actual flow rate of the individual phase
    Re_g                    = G * x * D / mu_g          # vapor
    Re_f                    = G * (1 - x) * D / mu_f    # liquid

    # 2. Friction factor for each phase
    if x == 1:                                  # No liquid
        f_f                 = 0                 # Just to be ok until next step
    elif Re_f < 1000:                           # Laminar todo fanning friction factor = 1 / 4 * Darcy fraction factor
        f_f                 = 16 / Re_f         # https://engineerexcel.com/darcy-vs-fanning-friction-factor/
    elif Re_f > 2000:                           # Fully-Turbulent
        f_f                 = 0.046 / (Re_f**0.2)
    else:                                       # Mixed
        # Weighting factor
        w                   = (Re_f - 1000) / (2000 - 1000)
        # Linear interpolation between laminar and turbulent
        f_f                 = (1 - w) * 16.0 / Re_f + w * 0.046 / (Re_f**0.2)

    if x == 0:                                  # No gas
        f_g                 = 0                 # Just to be ok until next step
    elif Re_g < 1000:                           # Laminar
        f_g                 = 16.0 / Re_g
    elif Re_g > 2000:                           # Fully-Turbulent
        f_g                 = 0.046 / (Re_g**0.2)
    else:                                       # Mixed
        # Weighting factor
        w                   = (Re_g - 1000) / (2000 - 1000)
        # Linear interpolation between laminar and turbulent
        f_g                 = (1 - w) * 16.0 / Re_g + w * 0.046 / (Re_g**0.2)

    # 3. Frictional pressure drop based on actual flow rate of each phase
    dp_dz_f                 = 2 * f_f * G**2 * (1-x)**2 * v_f / D
    dp_dz_g                 = 2 * f_g * G**2 * x**2 * v_g / D

    if x <= 0:                                  # Entirely liquid
        alpha               = 0.0
        dp_dz               = dp_dz_f
        return dp_dz, alpha
    if x >= 1:                                  # Entirely vapor
        alpha               = 1.0
        dp_dz               = dp_dz_g
        return dp_dz, alpha

    # 4. Lockhart-Martinelli parameter
    X                       = np.sqrt(dp_dz_f / dp_dz_g)

    # 5. Find the Constant based on the flow Re of each phase (using 1500 as the transitional Re to ensure continuity)
    # Calculate C if not passed in: https://www.scribd.com/document/423501015/Thesis-Po-Ya-Chuang-pdf, Page (95) 120/271
    if C == None:
        if Re_f > 1500 and Re_g > 1500:
            C               = 20.0
        elif Re_f < 1500 and Re_g > 1500:
            C               = 12.0
        elif Re_f > 1500 and Re_g < 1500:
            C               = 10.0
        else:
            C               = 5.0

    # 6. Two-phase multipliers for each phase
    phi_g2                  = 1 + C * X + X**2
    phi_f2                  = 1 + C / X + 1 / X**2

    # 7. Find gradient
    if dp_dz_g * phi_g2 > dp_dz_f * phi_f2:
        dp_dz               = dp_dz_g * phi_g2
    else:
        dp_dz               = dp_dz_f * phi_f2

    # 8. Void Fraction
    alpha                   = 1 - X / sqrt(X * X + 20 * X + 1)

    return dp_dz, alpha


def PHE_1phase_hdP(Inputs, JustGeo=False):
    """
    Based on the single-phase pressure drop and heat transfer correlations
    in VDI Heat Atlas Chapter N6: Pressure Drop and Heat Transfer in Plate Heat
    Exchangers by Holger Martin DOI: 10.1007/978-3-540-77877-6_66 Springer Verlag
    Outputs: for JustGeo=True >> Ap, Vchannel, Aflow, Dh, PHI
             for JustGeo=False >> Dh, h, Ap, DELTAP, Re_g, w_g, k_g, cp_g, Vchannel, Aflow
    ::
        =============================
        ||   __               __    ||
        ||  /  \             /  \   ||
        || |    |           |    |  ||  ===
        ||  \__/             \__/   ||   |
        ||                          ||   |
        ||             | <-  B   -> ||   |
        ||                          ||   |
        ||                          ||   |
        ||                          ||
        ||                          ||
        ||             |\           ||
        ||             | \          ||   Lp
        ||             |  \         ||
        ||             |   \        ||
        ||             |phi \       ||
        ||             |     \      ||   |
        ||                          ||   |
        ||   __               __    ||   |
        ||  /  \             /  \   ||   |
        || |    |           |    |  ||  ===
        ||  \__/             \__/   ||
        ||                          ||
        =============================
         | -----      Bp  --------- |

         phi is the inclination angle
    """

    # Plate parameters
    PlateAmplitude                          = Inputs['PlateAmplitude']  # mean flow channel gap, m
    PlateWavelength                         = Inputs['PlateWavelength'] # corrugation pitch, m
    InclinationAngle                        = Inputs['InclinationAngle']
    Bp                                      = Inputs['Bp']
    Lp                                      = Inputs['Lp']
    if JustGeo == False:
        AS                                  = Inputs['AS']
        T                                   = Inputs['T']
        p                                   = Inputs['p']
        m_dot_gap                           = Inputs['m_dot_gap']                   # mass flow rate per channel

    X                                       = 2 * pi * PlateAmplitude / PlateWavelength
    PHI                                     = 1 / 6 * (1 + sqrt(1 + X**2) + 4 * sqrt(1 + X**2 / 2))

    # The plane surface between the ports
    A0                                      = Bp * Lp

    # The plane surface of one plate
    Ap                                      = PHI * A0

    # The volume of one channel
    Vchannel                                = Bp * Lp * 2 * PlateAmplitude

    # Hydraulic diameter
    dh                                      = 4 * PlateAmplitude / PHI

    if JustGeo == True:
        return {'Ap': Ap, 'Vchannel': Vchannel, 'Aflow': 2 * PlateAmplitude * Bp, 'Dh': dh, 'PHI': PHI}
    else:
        # Also calculate the thermodynamics and pressure drop
        # Single phase Fluid properties
        AS.update(CP.PT_INPUTS, p, T)
        rho_g                               = AS.rhomass()                          # [kg/m^3]
        eta_g                               = AS.viscosity()                        # Viscosity[Pa-s]
        cp_g                                = AS.cpmass()                           # [J/kg-K]
        k_g                                 = AS.conductivity()                     # Thermal conductivity[W/m/K]

        Pr_g                                = cp_g * eta_g / k_g

        eta_g_w                             = eta_g                                 # TODO: allow for temperature dependence?
        w_g                                 = m_dot_gap / rho_g /(2 * PlateAmplitude * Bp)
        Re_g                                = rho_g*w_g*dh/eta_g

        # Calculate the friction factor zeta
        phi                                 = InclinationAngle

        if Re_g < 2000:
            zeta0                           = 64 / Re_g
            zeta1_0                         = 597 / Re_g + 3.85
        else:
            zeta0                           = (1.8 * log(Re_g) - 1.5)**(-2)
            zeta1_0                         = 39 / Re_g**0.289

        a                                   = 3.8
        b                                   = 0.18
        c                                   = 0.36

        zeta1                               = a * zeta1_0
        # RHS from Equation 18
        RHS                                 = cos(phi) / sqrt(b * tan(phi) + c * sin(phi) + zeta0 / cos(phi)) + (1 - cos(phi)) / sqrt(zeta1)
        zeta                                = 1 / RHS**2
        # Hagen number
        Hg                                  = zeta * Re_g**2 / 2

        # Constants for Nu correlation
        c_q                                 = 0.122
        q                                   = 0.374                                 # q=0.39
        # Nusselt number [-]
        Nu                                  = c_q * Pr_g**(1/3) * (eta_g / eta_g_w)**(1/6) * (2 * Hg * sin(2 * phi))**(q)

        # Heat transfer coefficient [W/m^2-K]
        h                                   = Nu * k_g / dh

        # Pressure drop
        DELTAP                              = Hg * eta_g**2 * Lp / (rho_g * dh**3)

        # There are quite a lot of things that might be useful to have access to
        # in outer functions, so pack up parameters into a dictionary
        Outputs = {
             'Dh':                          dh,                         # Hydraulic diameter [m]
             'h':                           h,                          # Heat transfer coefficient [W/m^2-K]
             'Ap':                          Ap,                         # Area of one plate [m^2]
             'DELTAP':                      DELTAP,                     # Pressure drop [Pa]
             'Re':                          Re_g,                       # Reynold number
             'U':                           w_g,                        # Velocity of fluid in channel [m/s]
             'k':                           k_g,                        # Thermal conductivity of fluid [W/m-K]
             'cp':                          cp_g,                       # Specific heat of fluid [J/kg-K]
             'Vchannel':                    Vchannel,                   # Volume of one channel [m^3]
             'Aflow':                       2 * PlateAmplitude * Bp     # Area of flow [m^2]
        }
        return Outputs


# todo check
def Kim_Mudawar_condensing_DPDZ_h(AS, G, Dh, x, Tbubble, Tdew, p, beta, C=None, satTransport=None):
    """
    This function return the pressure gradient and heat transfer coefficient for
    two phase fluid inside Micro-channel tube while CONDENSATION
    Correlations Based on:
    Kim and Mudawar (2012) "Universal approach to predicting two-phase
    frictional pressure drop and condensing mini/micro-channel flows", Int. J Heat Mass, 55, 3246-3261
    and
    Kim and Mudawar (2013) "Universal approach to predicting heat transfer coefficient
    for condensing min/micro-channel flow", Int. J Heat Mass, 56, 238-250
    """

    # Convert the quality, which might come in as a single numpy float value, to a float
    # With the conversion, >20x speedup in the LockhartMartinelli function, not clear why
    x       = float(x)

    if satTransport==None:
        # Calculate Necessary saturation properties
        AS.update(CP.QT_INPUTS,0.0,Tbubble)
        rho_f=AS.rhomass() #[kg/m^3]
        mu_f=AS.viscosity() #[Pa-s OR kg/m-s]
        cp_f=AS.cpmass() #[J/kg-K]
        k_f=AS.conductivity() #[W/m-K]
        AS.update(CP.QT_INPUTS,1.0,Tdew)
        rho_g=AS.rhomass() #[kg/m^3]
        mu_g=AS.viscosity() #[Pa-s OR kg/m-s]
    else:
        # Pull out of the dictionary
        rho_f=satTransport['rho_f']
        rho_g=satTransport['rho_g']
        mu_f=satTransport['mu_f']
        mu_g=satTransport['mu_g']
        cp_f=satTransport['cp_f']
        k_f=satTransport['k_f']

    AS.update(CP.PQ_INPUTS,p,x)
    sigma=AS.surface_tension() #surface tesnion [N/m]

    Pr_f = cp_f * mu_f / k_f #[-]

    Re_f = G*(1-x)*Dh/mu_f
    Re_g = G*x*Dh/mu_g


    if x==1: #No liquid
        f_f = 0 #Just to be ok until next step
    elif (Re_f<2000): #Laminar
        f_f = 16.0/Re_f
        if (beta<1):
            f_f = 24*(1-1.3553*beta+1.9467*beta*beta-1.7012*pow(beta,3)+0.9564*pow(beta,4)-0.2537*pow(beta,5))/Re_f
    elif (Re_f>=20000): #Fully-Turbulent
        f_f = 0.046*pow(Re_f,-0.2)
    else: #Transient
        f_f = 0.079*pow(Re_f,-0.25)

    if x==0: #No gas
        f_g = 0 #Just to be ok until next step
    elif (Re_g<2000): #Laminar
        f_g=16.0/Re_g
        if (beta<1):
            f_g = 24*(1-1.3553*beta+1.9467*beta*beta-1.7012*pow(beta,3)+0.9564*pow(beta,4)-0.2537*pow(beta,5))/Re_g
    elif (Re_g>=20000): #Fully-Turbulent
        f_g = 0.046*pow(Re_g,-0.2)
    else: #Transient
        f_g = 0.079*pow(Re_g,-0.25)

    Re_fo = G*Dh/mu_f
    Su_go = rho_g*sigma*Dh/pow(mu_g,2)

    dpdz_f = 2*f_f/rho_f*pow(G*(1-x),2)/Dh
    dpdz_g = 2*f_g/rho_g*pow(G*x,2)/Dh

    if x<=0:
        # Entirely liquid
        dpdz = dpdz_f
        AS.update(CP.QT_INPUTS,0.0,Tbubble)
        psat = AS.p() #pressure [Pa]
        h = f_h_1phase_MicroTube(G, Dh, Tbubble, psat, AS, Phase='SatLiq')[1]
        return dpdz, h
    if x>=1:
        #Entirely vapor
        dpdz = dpdz_g
        AS.update(CP.QT_INPUTS,1.0,Tdew)
        psat = AS.p() #pressure [Pa]
        h = f_h_1phase_MicroTube(G, Dh, Tdew, psat, AS, Phase='SatVap')[1]
        return dpdz, h

    X = sqrt(dpdz_f/dpdz_g)

    # Find the C coefficient (Calculate C if not passed, otherwise use the set value of C)
    if C==None:
        if (Re_f<2000 and Re_g<2000):
            C = 3.5e-5*pow(Re_fo,0.44)*pow(Su_go,0.50)*pow(rho_f/rho_g,0.48)
        elif (Re_f<2000 and Re_g>=2000):
            C = 0.0015*pow(Re_fo,0.59)*pow(Su_go,0.19)*pow(rho_f/rho_g,0.36)
        elif (Re_f>=2000 and Re_g<2000):
            C = 8.7e-4*pow(Re_fo,0.17)*pow(Su_go,0.50)*pow(rho_f/rho_g,0.14)
        elif (Re_f>=2000 and Re_g>=2000):
            C = 0.39*pow(Re_fo,0.03)*pow(Su_go,0.10)*pow(rho_f/rho_g,0.35)
    else:
        pass

    # Two-phase multiplier
    phi_f_square = 1.0 + C/X + 1.0/X**2
    phi_g_square = 1.0 + C*X + X**2

    # Find Condensing pressure drop griendient
    if dpdz_g*phi_g_square > dpdz_f*phi_f_square:
        dpdz=dpdz_g*phi_g_square
    else:
        dpdz=dpdz_f*phi_f_square

    #Use calculated Lockhart-Martinelli parameter
    Xtt = X
    # Simplified Lockhart-Martinelli paramter from Kim & Mudawar (2013) "Universal approach to predict HTC for condensing mini/micro-channel flow"
    #Xtt = pow(mu_f/mu_g,0.1) * pow((1-x)/x,0.9) * pow(rho_g/rho_f,0.5)

    # Modified Weber number
    if (Re_f <= 1250):
        We_star = 2.45 * pow(Re_g,0.64) / (pow(Su_go,0.3) * pow(1 + 1.09*pow(Xtt,0.039),0.4))
    else:
        We_star = 0.85 * pow(Re_g,0.79) * pow(Xtt,0.157) / (pow(Su_go,0.3) * pow(1 + 1.09*pow(Xtt,0.039),0.4)) * pow(pow(mu_g/mu_f,2) * (rho_f/rho_g),0.084)

    # Condensation Heat transfer coefficient
    if (We_star > 7*Xtt**0.2): ##for annual flow (smooth-annular, wavy-annular, transition)
        h = k_f/Dh * 0.048 * pow(Re_f,0.69) * pow(Pr_f,0.34) * sqrt(phi_g_square) / Xtt
    else: ##for slug and bubbly flow
        h = k_f/Dh *pow((0.048 * pow(Re_f,0.69) * pow(Pr_f,0.34) * sqrt(phi_g_square) / Xtt)**2 + (3.2e-7 * pow(Re_f,-0.38) * pow(Su_go,1.39))**2 ,0.5)

    return dpdz, h


def f_h_1phase_MicroTube(G, Dh, T, p, AS, Phase='Single'):
    """
    This function return the friction factor, heat transfer coefficient,
    and Reynold's number for single phase fluid inside flat plate tube
    Micro-channel HX
    """
    if Phase=="SatVap":
        AS.update(CP.QT_INPUTS,1.0,T)
        mu = AS.viscosity() #[Pa-s OR kg/m-s]
        cp = AS.cpmass() #[J/kg-K]
        k = AS.conductivity() #[W/m-K]
        rho = AS.rhomass() #[kg/m^3]
    elif Phase =="SatLiq":
        AS.update(CP.QT_INPUTS,0.0,T)
        mu = AS.viscosity() #[Pa-s OR kg/m-s]
        cp = AS.cpmass() #[J/kg-K]
        k = AS.conductivity() #[W/m-K]
        rho = AS.rhomass() #[kg/m^3]
    else:
        AS.update(CP.PT_INPUTS,p,T)
        mu = AS.viscosity() #[Pa-s OR kg/m-s]
        cp = AS.cpmass() #[J/kg-K]
        k = AS.conductivity() #[W/m-K]
        rho = AS.rhomass() #[kg/m^3]

    Pr = cp * mu / k #[-]

    Re=G*Dh/mu

    # Friction factor of Churchill (Darcy Friction factor where f_laminar=64/Re)
    e_D = 0.0
    A = ((-2.457 * log( (7.0 / Re)**(0.9) + 0.27 * e_D)))**16
    B = (37530.0 / Re)**16
    f = 8.0 * ((8.0/Re)**12.0 + 1.0 / (A + B)**(1.5))**(1/12)

    # Heat Transfer coefficient of Gnielinski
    Nu = (f/8.0)*(Re-1000.0)*Pr/(1.0+12.7*sqrt(f/8.0)*(Pr**(2/3)-1)) #[-]
    h = k*Nu/Dh #W/m^2-K
    return (f, h, Re)


# ------------------------------------------------------------------
if __name__ == '__main__':
    # test simps
    x = np.array([1., 2.])
    y = np.array([1., 4.])
    res = simps(y, x)
    print(x, y, res)
