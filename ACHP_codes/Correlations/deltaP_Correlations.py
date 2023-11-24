from __future__         import division, print_function, absolute_import
from math               import pi,log,sqrt,exp,cos,sin,tan,log10
from scipy.integrate    import quad,quadrature,trapz,simps,fixed_quad
from scipy.optimize     import brentq,fsolve
import numpy            as np
import CoolProp         as CP

# Machine precision
machine_eps = np.finfo(float).eps
print(machine_eps)


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


# ------------------------------------------------------------------
if __name__ == '__main__':
    # test simps
    x = np.array([1., 2.])
    y = np.array([1., 4.])
    res = simps(y, x)
    print(x, y, res)
