from __future__         import division, print_function, absolute_import
from math               import pi,log,sqrt,exp,cos,sin,tan,log10
from scipy.integrate    import quad,quadrature,trapz,simps,fixed_quad
from scipy.optimize     import brentq,fsolve
import numpy            as np
import CoolProp         as CP

# Machine precision
machine_eps = np.finfo(float).eps
print(machine_eps, '; num=1')


# ------------------------------------------------------------------
def Phase_ph(AS, p, h, T_bubble, T_dew, rho_sat_L, rho_sat_V):
    """
    Convenience function to return just the Phase without temperature or density
    AS : AbstractState with the refrigerant name and backend
    T_dew: dew point of the air temperature, cooler than which the airborne water will condense,from water vapor to liquid
    the dew point is affected by the air's humidity, the more moisture the air contains, the higher its dew point
    """
    (T, rho, Phase) = TrhoPhase_ph(AS, p, h, T_bubble, T_dew, rho_sat_L, rho_sat_V)
    return Phase


def TrhoPhase_ph(AS, p, h, T_bubble, T_dew, rho_sat_L=None, rho_sat_V=None):
    """
    Convenience function to find temperature, density, and phase of fluid as a function of pressure and enthalpy
    AS : AbstractState with the refrigerant name and backend
    """

    if 'IncompressibleBackend' in AS.backend_name():
        # It is subcooled;
        # This backend provides the thermophysical properties for incompressible pure fluids, incompressible mixtures, and brines
        AS.update(CP.HmassP_INPUTS, h, p)
        T           = AS.T()        # [K]
        rho         = AS.rhomass()  # [kg/m^3]
        return T, rho, 'Subcooled'

    # Check if it is supercritical
    if p >= AS.p_critical():
        AS.update(CP.HmassP_INPUTS, h, p)
        T           = AS.T()        # [K]
        rho         = AS.rhomass()  # [kg/m^3]
        if T >= AS.T_critical():
            return T, rho, 'Supercritical'
        else:
            return T, rho, 'Supercrit_liq'
    else:                           # It is subcritical
        if rho_sat_L == None:
            AS.update(CP.QT_INPUTS, 0.0, T_bubble)
            rho_sat_L   = AS.rhomass()  # [kg/m^3]
            AS.update(CP.QT_INPUTS, 1.0, T_dew)
            rho_sat_V   = AS.rhomass()  # [kg/m^3]
        v_sat_L         = 1 / rho_sat_L
        v_sat_V         = 1 / rho_sat_V
        AS.update(CP.DmassT_INPUTS, rho_sat_L, T_bubble)
        h_sat_L         = AS.hmass()    # [J/kg]
        AS.update(CP.DmassT_INPUTS, rho_sat_V, T_dew)
        h_sat_V         = AS.hmass()    # [J/kg]

        if h > h_sat_V:                 # It's superheated
            AS.update(CP.HmassP_INPUTS, h, p)
            cp          = AS.cpmass()   # [J/kg-]
            T           = AS.T()        # [K]
            rho         = AS.rhomass()  # [kg/m^3]
            return T, rho, 'Superheated'
        elif h < h_sat_L:               # It's subcooled
            AS.update(CP.HmassP_INPUTS, h, p)
            cp          = AS.cpmass()   # [J/kg-]
            T           = AS.T()        # [K]
            rho         = AS.rhomass()  # [kg/m^3]
            return T, rho, 'Subcooled'
        else:                           # It's two-phase
            x           = (h - h_sat_L) / (h_sat_V - h_sat_L) # [-] quality
            v           = x * v_sat_V + (1 - x) * v_sat_L   # [m^3/kg] specific volume
            T           = x * T_dew + (1 - x) * T_bubble    # [K] absolute temperature
            rho         = 1 / v         # [kg/m^3]
            return T, rho, 'TwoPhase'


def TwoPhaseDensity(AS, x_min, x_max, T_dew, T_bubble, slipModel='Zivi'):
    """
    function to obtain the average density in the two-phase region
    AS : AbstractState with the refrigerant name and backend
    bubble point:
    the bubble point is the temperature (at a given pressure)
    where the first bubble of vapor is formed when heating a liquid consisting of two or more components;
    For a single component the bubble point and the dew point are the same and are referred to as the boiling point.
    """
    AS.update(CP.QT_INPUTS, 1.0, T_dew)
    rho_g               = AS.rhomass()          # [kg/m^3]
    AS.update(CP.QT_INPUTS, 0.0, T_bubble)
    rho_f               = AS.rhomass()          # [kg/m^3]

    if slipModel == 'Zivi':
        S               = pow(rho_f / rho_g, 0.3333)
    elif slipModel == 'Homogeneous':
        S               = 1
    else:
        raise ValueError("slipModel must be either 'Zivi' or 'Homogeneous'")
    C                   = S * rho_g / rho_f

    if x_min + 5 * machine_eps < 0 or x_max - 10 * machine_eps > 1.0:
        raise ValueError('Quality must be between 0 and 1, ' + 'x_min: ' + str(x_min + 5 * machine_eps) + ', x_max: '+str(x_max - 10 * machine_eps))
    # Avoid the zero and one qualities (undefined integral)
    if x_min == x_max:
        alpha_average   = 1 / (1 + C * (1 - x_min) / x_min)
    else:
        if x_min >= 1.0:
            alpha_average = 1.0
        elif x_max <= 0.0:
            alpha_average = 0.0
        else:
            alpha_average = -(C * (np.log(((x_max - 1.0) * C - x_max) / ((x_min - 1.0) * C - x_min)) + x_max - x_min) - x_max + x_min) \
                            / (C**2 - 2 * C + 1) / (x_max - x_min) # todo find the source?

    return alpha_average * rho_g + (1 - alpha_average) * rho_f


# ------------------------------------------------------------------
if __name__ == '__main__':
    # test simps
    x = np.array([1., 2.])
    y = np.array([1., 4.])
    res = simps(y, x)
    print(x, y, res)
