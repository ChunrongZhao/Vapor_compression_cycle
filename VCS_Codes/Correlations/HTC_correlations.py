from math               import pi,log,sqrt,exp,cos,sin,tan,log10
from scipy.integrate    import quad,quadrature,trapz,simps,fixed_quad
from scipy.optimize     import brentq,fsolve
import numpy            as np
import CoolProp         as CP

# Machine precision
machine_eps = np.finfo(float).eps
print(machine_eps, '; num=3')


# ------------------------------------------------------------------------------------------
def f_h_1phase_Channel(m_dot, d_H, U, T, p, AS, Phase):

    if Phase == "SatVap":
        AS.update(CP.QT_INPUTS, 1.0, T)
        mu                          = AS.viscosity()    # [Pa-s OR kg/m-s]
        cp                          = AS.cpmass()       # [J/kg-K]
        k                           = AS.conductivity() # [W/m-K]
        rho                         = AS.rhomass()      # [kg/m^3]
    elif Phase == "SatLiq":
        AS.update(CP.QT_INPUTS, 0.0, T)
        mu                          = AS.viscosity()    # [Pa-s OR kg/m-s]
        cp                          = AS.cpmass()       # [J/kg-K]
        k                           = AS.conductivity() # [W/m-K]
        rho                         = AS.rhomass()      # [kg/m^3]
    else:
        AS.update(CP.PT_INPUTS, p, T)
        mu                          = AS.viscosity()    # [Pa-s OR kg/m-s]
        cp                          = AS.cpmass()       # [J/kg-K]
        k                           = AS.conductivity() # [W/m-K]
        rho                         = AS.rhomass()      # [kg/m^3]

    Pr                              = cp * mu / k   # [-]
    Re                              = rho * U * d_H / mu

    # Friction factor of Churchill (Darcy Friction factor where f_laminar=64/Re)
    e_D                             = 0
    A                               = (-2.457 * log((7.0 / Re) ** 0.9 + 0.27 * e_D)) ** 16
    B                               = (37530.0 / Re)**16
    f                               = 8 * ((8/Re) ** 12.0 + 1 / (A + B) ** 1.5) ** (1 / 12)

    # Heat Transfer coefficient of Gnielinski for high Re,
    if Re > 1000:
        Nu = (f / 8) * (Re - 1000.) * Pr / (1 + 12.7 * sqrt(f / 8) * (Pr ** 0.66666 - 1)) # [-]
    else:
        Nu = 3.66
    h                               = k * Nu / d_H #W/m^2-K
    return (f, h, Re)


# ------------------------------------------------------------------
# evaporation
def ShahEvaporation_Average(x_min, x_max, AS, G, D, p, q_flux, T_bubble, T_dew):
    """
    Returns the average heat transfer coefficient between qualities of x_min and x_max.

    Required parameters:
    * x_min : The minimum quality for the range [-]
    * x_max : The maximum quality for the range [-]
    * AS : AbstractState with the refrigerant name and backend
    * G : Mass flux [kg/m^2/s]
    * D : Diameter of tube [m]
    * p : Pressure [Pa]
    * q_flux : Heat transfer flux [W/m^2]
    * T_bubble : Bubble point temperature of refrigerant [K]
    * T_dew : Dew point temperature of refrigerant [K]
    """
    # ********************************
    #        Necessary Properties
    # ********************************
    AS.update(CP.QT_INPUTS, 0.0, T_bubble)      # x = 0, pure liquid
    rho_f                   = AS.rhomass()      # [kg/m^3]
    mu_f                    = AS.viscosity()    # [Pa-s OR kg/m-s]
    cp_f                    = AS.cpmass()       # [J/kg-K]
    k_f                     = AS.conductivity() # [W/m-K]
    h_l                     = AS.hmass()        # [J/kg]

    AS.update(CP.QT_INPUTS, 1.0, T_dew)
    rho_g                   = AS.rhomass()      # [kg/m^3]
    mu_g                    = AS.viscosity()    # [Pa-s OR kg/m-s]
    cp_g                    = AS.cpmass()       # [J/kg-K]
    k_g                     = AS.conductivity() # [W/m-K]
    h_v                     = AS.hmass()        # [J/kg]

    h_fg                    = h_v - h_l         # [J/kg]
    Pr_f                    = cp_f * mu_f / k_f # [-]
    Pr_g                    = cp_g * mu_g / k_g # [-]

    g_grav                  = 9.81              # [m/s^2]

    # Shah evaporation correlation
    Fr_L                    = G**2 / (rho_f * rho_f * g_grav * D)   # [-]
    Bo                      = q_flux / (G * h_fg)                   # [-]

    if Bo < 0:
        raise ValueError('Heat flux for Shah Evaporation must be positive')

    if Bo > 0.0011:
        F                   = 14.7
    else:
        F                   = 15.43

    # Pure vapor single-phase heat transfer coefficient
    h_g                     = 0.023 * (G * D / mu_g)**0.8 * Pr_g**0.4 * k_g / D     # [W/m^2-K]

    def ShahEvaporation(x):
        """Shah, M.M. (1979). A general correlation for heat transfer during film condensation inside
            pipes. International Journal of Heat and Mass Transfer, Vol. 22, No. 4, 547-556. (3)"""

        if abs(1-x) < 5 * machine_eps:
            return h_g

        # If the quality is above 0.999, linearly interpolate to avoid division by zero
        if x > 0.999:
            h_1             = ShahEvaporation(1.0)                                  # Fully fry
            h_999           = ShahEvaporation(0.999)                                # At a quality of 0.999
            return (h_1 - h_999) / (1.0 - 0.999) * (x - 0.999) + h_999              # Linear interpolation

        if abs(x) < 5 * machine_eps:
            # https://www.scribd.com/document/423501015/Thesis-Po-Ya-Chuang-pdf, Page (91) 116/271
            h_L             = 0.023 * (G * (1 - x) * D / mu_f)**0.8 * Pr_f**0.4 * k_f / D   # [W/m^2-K] fully liquid
            return h_L
        else:
            h_L             = 0.023 * (G * (1 - x) * D / mu_f)**0.8 * Pr_f**0.4 * k_f / D   # [W/m^2-K]

        Co                  = (1 / x - 1)**(0.8) * (rho_g / rho_f)**(0.5)                   # [-]

        if Fr_L >= 0.04:
            N               = Co
        else:
            N               = 0.38 * Fr_L**(-0.3) * Co

        # the ratio of h_TP / h_l
        psi_cb              = 1.8 / N**(0.8)

        if 0.1 < N and N <= 1.0:
            psi_bs          = F * Bo**0.5 * exp(2.74 * N**(-0.1))
            psi             = max([psi_bs, psi_cb])
        else:
            if N > 1.0:
                if Bo > 0.00003:
                    psi_nb  = 230 * Bo**0.5
                else:
                    psi_nb  = 1.0 + 46.0 * Bo**0.5
                psi         = max([psi_nb, psi_cb])
            else:
                psi_bs      = F * (Bo)**(0.5) * exp(2.47 * N**(-0.15))
                psi         = max([psi_bs, psi_cb])

        return psi * h_L    # HTC_2phase [W/m^2-K]

    # Calculate h over the range of x
    x                       = np.linspace(x_min, x_max, 100)
    h                       = np.zeros_like(x)

    for i in range(len(x)):
        h[i]                = ShahEvaporation(x[i])

    # if x_min == x_max, or they are really close to being the same
    if abs(x_max - x_min) < 5 * machine_eps:
        # return just one of the edge values
        return h[0]
    else:
        # Use Simpson's rule to carry out numerical integration to get average
        return simps(h, x) / (x_max - x_min)
