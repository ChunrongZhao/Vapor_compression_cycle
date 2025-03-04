from __future__         import division, print_function, absolute_import
from math               import pi,log,sqrt,exp,cos,sin,tan,log10
from scipy.integrate    import quad,quadrature,trapz,simps,fixed_quad
from scipy.optimize     import brentq,fsolve
import numpy            as np
import CoolProp         as CP

# Machine precision
machine_eps = np.finfo(float).eps
print(machine_eps, '; num=3')


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


def KandlikarEvaporation_average(x_min, x_max, AS, G, D, p, q_flux, T_bubble, T_dew):
    """
    Kandlikar (1990) recommended by Petterson et al. (2000) for CO2, Heat transfer and pressure drop for flow supercritical
    and subcritical CO2 in microchannel tubes
    All details for this correlation are available in Ding Li Thesis (Appendix C):
    "INVESTIGATION OF AN EJECTOR-EXPANSION DEVICE IN A TRANSCRITICAL CARBON DIOXIDE CYCLE FOR MILITARY ECU APPLICATIONS"

    Returns the average heat transfer coefficient between qualities of x_min and x_max.

    Required parameters:
    * x_min : The minimum quality for the range [-]
    * x_max : The maximum quality for the range [-]
    * AS : AbstractState with the refrigerant name and backend
    * G : Mass flux [kg/m^2/s]
    * D : Diameter of tube [m]
    * p : Pressure [Pa]
    * q_flux : Heat transfer flux [W/m^2]
    * T_bubble : Bubblepoint temperature of refrigerant [K]
    * T_dew : Dewpoint temperature of refrigerant [K]
    """
    # ********************************
    #        Necessary Properties
    # ********************************
    AS.update(CP.QT_INPUTS, 0.0, T_bubble)
    rho_f                       = AS.rhomass()          # [kg/m^3]
    mu_f                        = AS.viscosity()        # [Pa-s OR kg/m-s]
    cp_f                        = AS.cpmass()           # [J/kg-K]
    k_f                         = AS.conductivity()     # [W/m-K]
    h_l                         = AS.hmass()            # [J/kg]

    AS.update(CP.QT_INPUTS, 1.0, T_dew)
    rho_g                       = AS.rhomass()          # [kg/m^3]
    mu_g                        = AS.viscosity()        # [Pa-s OR kg/m-s]
    cp_g                        = AS.cpmass()           # [J/kg-K]
    k_g                         = AS.conductivity()     # [W/m-K]
    h_v                         = AS.hmass()            # [J/kg]

    h_fg                        = h_v - h_l             # [J/kg]
    Pr_f                        = cp_f * mu_f / k_f     # [-]
    Pr_g                        = cp_g * mu_g / k_g     # [-]

    g_grav                      = 9.81                  # [m/s^2]

    # Petterson evaporation correlation
    Fr_L                        = G**2 / (rho_f * rho_f * g_grav * D)   # [-]
    Bo                          = q_flux / (G * h_fg)                   # [-]

    if Bo < 0:
        raise ValueError('Heat flux for Petterson Evaporation must be positive')

    # Forster and Zuber multiplier depend on fluid type. CO2 is not available, therefore F_fl=1 (for water) is selected.
    """ the parameter FFl is a fluid-surface specific parameter that is recommended to be set at 1 
    for surface-fluid combinations that have not yet been researched or in general for preliminary analysis."""
    F_fl                        = 1
    # Kandlikar correlation constants for CO2
    c_c_1                       = 1.1360
    c_c_2                       = -0.9
    c_c_3                       = 667.2
    c_c_4                       = 0.7
    c_n_1                       = 0.6683
    c_n_2                       = -0.2
    c_n_3                       = 1058.0
    c_n_4                       = 0.7

    if Fr_L > 0.4:
        c_c_5                   = 0.0
        c_n_5                   = 0.0
    else:
        c_c_5                   = 0.3
        c_n_5                   = 0.3

    # Pure vapor single-phase heat transfer coefficient
    h_g                         = 0.023 * (G * D / mu_g)**(0.8) * Pr_g**(0.4) * k_g / D     # [W/m^2-K]

    def KandlikarEvaporation(x):

        if abs(1-x) < 5 * machine_eps:
            return h_g
        # If the quality is above 0.999, linearly interpolate to avoid division by zero
        if x > 0.999:
            h_1                 = KandlikarEvaporation(1.0)                                 # Fully fry
            h_999               = KandlikarEvaporation(0.999)                               # At a quality of 0.999
            return (h_1 - h_999) / (1.0 - 0.999) * (x - 0.999) + h_999                      # Linear interpolation
        if abs(x) < 5 * machine_eps:
            h_L                 = 0.023 * (G * (1 - x) * D / mu_f)**(0.8) * Pr_f**(0.4) * k_f / D   # [W/m^2-K]
            return h_L
        else:
            h_L                 = 0.023 * (G * (1 - x) * D / mu_f)**(0.8) * Pr_f**(0.4) * k_f / D   # [W/m^2-K]

        Co                      = (1 / x - 1)**(0.8) * (rho_g / rho_f)**(0.5)                       # [-]

        """https://repository.tudelft.nl/islandora/object/uuid:9a16e0f8-722d-4a17-9f79-96cea2de6906 Page (40) 62/132"""
        # HTC due to convective boiling
        h_c                     = h_L * (c_c_1 * pow(Co, c_c_2) * pow((25.0 * Fr_L), c_c_5) + c_c_3 * pow(Bo, c_c_4) * F_fl)
        # HTC due to nucleate boiling
        h_n                     = h_L * (c_n_1 * pow(Co, c_n_2) * pow((25.0 * Fr_L), c_n_5) + c_n_3 * pow(Bo, c_n_4) * F_fl)

        # This was found in ACCO2 model, however Petterson (2000) recommends to take the max of h_n and h_c
        h                       = max(h_c,h_n)
        return h        # HTC_2phase

    # Calculate h over the range of x
    x                           = np.linspace(x_min, x_max, 100)
    h                           = np.zeros_like(x)
    for i in range(len(x)):
        h[i]                    = KandlikarEvaporation(x[i])

    # if x_min == x_max, or they are really really close to being the same
    if abs(x_max - x_min) < 5 * machine_eps:
        # return just one of the edge values
        return h[0]
    else:
        # Use Simpson's rule to carry out numerical integration to get average
        return simps(h, x) / (x_max - x_min)


# condensation
def LongoCondensation(x_avg, G, dh, AS, T_sat_L, T_sat_V):

    AS.update(CP.QT_INPUTS, 1.0, T_sat_V)
    rho_V                       = AS.rhomass()      # [kg/m^3]
    AS.update(CP.QT_INPUTS, 0.0, T_sat_L)
    rho_L                       = AS.rhomass()      # [kg/m^3]
    mu_L                        = AS.viscosity()    # [Pa-s OR kg/m-s]
    cp_L                        = AS.cpmass()       # [J/kg-K]
    k_L                         = AS.conductivity() # [W/m-K]
    Pr_L                        = cp_L * mu_L / k_L # [-]

    Re_eq                       = G * ((1 - x_avg) + x_avg * sqrt(rho_L / rho_V)) * dh / mu_L

    if Re_eq < 1750:
        Nu                      = 60 * Pr_L**(1/3)
    else:
        Nu                      = ((75 - 60) / (3000 - 1750) * (Re_eq - 1750) + 60) * Pr_L**(1/3)
    h                           = Nu * k_L / dh
    return h


def ShahCondensation_Average(x_min, x_max, AS, G, D, p, T_sat_L, T_sat_V):
    # ********************************
    #    Necessary Properties
    #    Calculated outside the quadrature integration for speed
    # ********************************
    AS.update(CP.QT_INPUTS, 0.0, T_sat_L)
    mu_f                        = AS.viscosity()            # [Pa-s OR kg/m-s]
    cp_f                        = AS.cpmass()               # [J/kg-K]
    k_f                         = AS.conductivity()         # [W/m-K]
    Pr_f                        = cp_f * mu_f / k_f         # [-]
    pcrit                       = AS.p_critical()           # [Pa]
    Pstar                       = p / pcrit
    h_L                         = 0.023 * (G * D / mu_f)**(0.8) * Pr_f**(0.4) * k_f / D     # [W/m^2-K]

    def ShahCondensation(x):
        return h_L * ((1 - x)**(0.8) + (3.8 * x**(0.76) * (1 - x)**(0.04)) / (Pstar**(0.38)) )

    if not x_min == x_max:
        # A proper range is given
        return quad(ShahCondensation, x_min, x_max)[0] / (x_max - x_min)
    else:
        # A single value is given
        return ShahCondensation(x_min)


# supercritical
def Petterson_supercritical_average(Tout, Tin, T_w, AS, G, OD, ID, D_l, mdot, p, q_flux_w):
    """
    Petterson et al. (2000), Heat transfer and pressure drop for flow supercritical and subcritical CO2 in microchannel tubes
    All details for this correlation are available in Ding Li Thesis (Appendix B):
    "INVESTIGATION OF AN EJECTOR-EXPANSION DEVICE IN A TRANSCRITICAL CARBON DIOXIDE CYCLE FOR MILITARY ECU APPLICATIONS"
    """

    def SuperCriticalCondensation_h(T,   T_w, AS, G, OD, ID, D_l, mdot, p, q_flux_w):
        """return h value"""
        return Petterson_supercritical(T, T_w, AS, G, OD, ID, D_l, mdot, p, q_flux_w)[0]

    def SuperCriticalCondensation_f(T,   T_w, AS, G, OD, ID, D_l, mdot, p, q_flux_w):
        """return f value"""
        return Petterson_supercritical(T, T_w, AS, G, OD, ID, D_l, mdot, p, q_flux_w)[1]

    def SuperCriticalCondensation_cp(T,  T_w, AS, G, OD, ID, D_l, mdot, p, q_flux_w):
        """return cp value"""
        return Petterson_supercritical(T, T_w, AS, G, OD, ID, D_l, mdot, p, q_flux_w)[2]

    def SuperCriticalCondensation_rho(T, T_w, AS, G, OD, ID, D_l, mdot, p, q_flux_w):
        """return rho value"""
        return Petterson_supercritical(T, T_w, AS, G, OD, ID, D_l, mdot, p, q_flux_w)[3]

    if not Tout == Tin:
        # A proper range is given
        h   = quad(SuperCriticalCondensation_h,   Tin, Tout, args=(T_w, AS, G, OD, ID, D_l, mdot, p, q_flux_w))[0]/(Tout-Tin)
        f   = quad(SuperCriticalCondensation_f,   Tin, Tout, args=(T_w, AS, G, OD, ID, D_l, mdot, p, q_flux_w))[0]/(Tout-Tin)
        cp  = quad(SuperCriticalCondensation_cp,  Tin, Tout, args=(T_w, AS, G, OD, ID, D_l, mdot, p, q_flux_w))[0]/(Tout-Tin)
        rho = quad(SuperCriticalCondensation_rho, Tin, Tout, args=(T_w, AS, G, OD, ID, D_l, mdot, p, q_flux_w))[0]/(Tout-Tin)
        return (h, f, cp, rho)
    else:
        # A single value is given
        return Petterson_supercritical(Tout, T_w, AS, G, OD, ID, D_l, mdot, p, q_flux_w)


def Petterson_supercritical(T, T_w, AS, G, OD, ID, D_l, mdot, p, q_flux_w):
    AS.update(CP.PT_INPUTS, p, T_w)
    h_w                             = AS.hmass()                    # [J/kg]
    mu_w                            = AS.viscosity()                # [Pa-s OR kg/m-s]
    cp_w                            = AS.cpmass()                   # [J/kg-K]
    k_w                             = AS.conductivity()             # [W/m-K]
    rho_w                           = AS.rhomass()                  # [kg/m^3]
    Pr_w                            = cp_w * mu_w / k_w             # [-]

    AS.update(CP.PT_INPUTS, p, T)
    h                               = AS.hmass()                    # [J/kg]
    mu                              = AS.viscosity()                # [Pa-s OR kg/m-s]
    cp                              = AS.cpmass()                   # [J/kg-K]
    k                               = AS.conductivity()             # [W/m-K]
    rho                             = AS.rhomass()                  # [kg/m^3]
    Pr                              = cp * mu / k                   # [-]

    if mdot == 0:                                                   # For the case of Micro-channel
        Dh                          = OD
        Re                          = G * Dh / mu
        Re_w                        = G * Dh / mu_w
    else:                                                           # for the case of fin-and-tube
        Dh                          = OD - ID
        Area                        = pi * (OD**2 - ID**2) / 4.0
        u                           = mdot / (Area * rho)
        Re                          = rho * u * Dh / mu
        Re_w                        = Re                            # rho_w*u*Dh/mu_w

    if G > 350:
        e_D                         = 0                             # smooth pipe
        f                           = (-1.8 * log10(6.9 / Re + (1 / 3.7 * e_D)**1.11))**(-2) / 4
        Nu_m                        = (f/8) * (Re-1000) * Pr / (1 + 12.7 * sqrt(f/8) * (Pr**(2/3) - 1)) * (1 + D_l**(2/3))
        Nu                          = Nu_m * (Pr/Pr_w)**0.11

    else:                                                           # G < 350

        M                           = 0.001                         # [kg/J]
        K                           = 0.00041                       # [kg/J]

        cp_avg                      = (h - h_w) / (T - T_w)

        if cp_avg / cp_w <= 1:
            n                       = 0.66 - K * (q_flux_w / G)
        else:                                                       # cp_avg/cp_w <1
            n                       = 0.9 - K * (q_flux_w / G)

        f0                          = (0.79 * log(Re) - 1.64)**(-2)

        g                           = 9.81
        # coefficient of thermal expansion
        beta                        = AS.isobaric_expansion_coefficient()   # [1/K]
        # Grashoff number
        Gr                          = g * beta * (T_w - T) * Dh**3 / (mu/rho)**2
        if Gr / Re**2 < 5e-4:
            f                       = f0 * (mu_w/mu)**0.22
        elif Gr/Re**2 >= 5e-4 and G/Re**2 < 0.3:
            f                       = 2.15 * f0 * (mu_w/mu)**0.22 * (Gr/Re)**0.1
        else:                                                       # use f0 for friction factor
            f                       = f0

        Nu_w_ppk                    = (f0/8) * Re_w * Pr_w / (1.07 + 12.7 * sqrt(f/8) * (Pr_w**(2/3) - 1))

        Nu                          = Nu_w_ppk * (1 - M * q_flux_w / G) * (cp_avg / cp_w)**n

    h                               = k * Nu / Dh                   # [W/m^2-K]

    return (h, f, cp, rho)


# single phase flow
def f_h_1phase_Tube(m_dot, ID, T, p, AS):
    """
    Convenience function to run annular model for tube.  Tube is a degenerate case of annulus with inner diameter of 0

    """
    return f_h_1phase_Annulus(m_dot, ID, 0.0, T, p, AS, Phase='Single')


def f_h_1phase_Annulus(m_dot, OD, ID, T, p, AS, Phase='Single'):
    """
    This function return the friction factor, heat transfer coefficient,
    and Reynolds number for single phase fluid inside annular pipe
    """
    if Phase == "SatVap":
        AS.update(CP.QT_INPUTS, 1.0, T)
        mu                      = AS.viscosity()            # [Pa-s OR kg/m-s]
        cp                      = AS.cpmass()               # [J/kg-K]
        k                       = AS.conductivity()         # [W/m-K]
        rho                     = AS.rhomass()              # [kg/m^3]
    elif Phase == "SatLiq":
        AS.update(CP.QT_INPUTS, 0.0, T)
        mu                      = AS.viscosity()            # [Pa-s OR kg/m-s]
        cp                      = AS.cpmass()               # [J/kg-K]
        k                       = AS.conductivity()         # [W/m-K]
        rho                     = AS.rhomass()              # [kg/m^3]
    else:
        AS.update(CP.PT_INPUTS, p, T)
        mu                      = AS.viscosity()            # [Pa-s OR kg/m-s]
        cp                      = AS.cpmass()               # [J/kg-K]
        k                       = AS.conductivity()         # [W/m-K]
        rho                     = AS.rhomass()              # [kg/m^3]

    Pr                          = cp * mu / k               # [-]

    Dh                          = OD - ID
    Area                        = pi * (OD**2 - ID**2) / 4.0
    u                           = m_dot / (Area * rho)
    Re                          = rho * u * Dh / mu

    # Friction factor of Churchill (Darcy Friction factor where f_laminar=64/Re)
    e_D                         = 0
    A                           = (-2.457 * log((7.0 / Re)**0.9 + 0.27 * e_D)) ** 16
    B                           = (37530.0 / Re)**16
    f                           = 8 * ((8 / Re) ** 12.0 + 1 / (A + B)**1.5) ** (1 / 12)

    # Heat Transfer coefficient of Gnielinski
    Nu                          = (f / 8) * (Re - 1000) * Pr / (1 + 12.7 * sqrt(f / 8) * (Pr**(2/3) - 1)) # [-]
    h                           = k * Nu / Dh                       # W/m^2-K
    return (f, h, Re)


def f_h_1phase_Channel(mdot, W, H, T, p, AS, Phase):

    if Phase == "SatVap":
        AS.update(CP.QT_INPUTS,1.0,T)
        mu = AS.viscosity() #[Pa-s OR kg/m-s]
        cp = AS.cpmass() #[J/kg-K]
        k = AS.conductivity() #[W/m-K]
        rho = AS.rhomass() #[kg/m^3]
    elif Phase == "SatLiq":
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

    Pr      = cp * mu / k #[-]

    Dh      = 2 * H * W / (H + W)
    Area    = W * H
    u       = mdot / (Area * rho)
    Re      = rho * u * (Dh) / mu

    # Friction factor of Churchill (Darcy Friction factor where f_laminar=64/Re)
    e_D = 0
    A = ((-2.457 * log( (7.0 / Re)**(0.9) + 0.27 * e_D)))**16
    B = (37530.0 / Re)**16
    f = 8 * ((8/Re)**12.0 + 1 / (A + B)**(1.5))**(1/12)

    # Heat Transfer coefficient of Gnielinski for high Re,
    if Re > 1000:
        Nu = (f / 8) * (Re - 1000.) * Pr / (1 + 12.7 * sqrt(f / 8) * (Pr**(0.66666) - 1)) #[-]
    else:
        Nu = 3.66
    h = k * Nu / Dh #W/m^2-K

    print('Re = ', Re, "h_vap = ", h)
    return (f, h, Re)


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


# todo check
def Bertsch_MC(x, AS, G, Dh, q_flux, L, T_bubble, T_dew):
    """
    This function return the heat transfer coefficient for two phase fluid
    inside Micro-channel tube
    Correlatation is based on Bertsch (2009)
    """
    # Define properties
    AS.update(CP.QT_INPUTS, 0.0, T_bubble)
    k_L                                     = AS.conductivity()             # [W/m/K]
    cp_L                                    = AS.cpmass()                   # [J/kg-K]
    mu_L                                    = AS.viscosity()                # [Pa-s OR kg/m-s]
    rho_L                                   = AS.rhomass()                  # [kg/m^3]

    AS.update(CP.QT_INPUTS, 1.0, T_dew)
    k_G                                     = AS.conductivity()             # [W/m/K]
    cp_G                                    = AS.cpmass()                   # [J/kg-K]
    mu_G                                    = AS.viscosity()                # [Pa-s OR kg/m-s]
    rho_G                                   = AS.rhomass()                  # [kg/m^3]

    p                                       = AS.p()                        # saturation pressure [Pa] @ Tdew
    pc                                      = AS.p_critical()               # critical pressure [Pa]
    pr                                      = p / pc
    M                                       = AS.molar_mass()               # molar mass [kg/mol]

    AS.update(CP.PQ_INPUTS, p, x)
    sig                                     = AS.surface_tension()          # surface tension [N/m]
    g                                       = 9.81

    # if Ref=='R290':
    #    sig=55.28*(1-T_dew/369.818)**(1.258)/1000.
    # elif Ref=='R410A':
    #  From Okada 1999 "Surface Tension of HFC Refrigerant Mixtures"
    #    sig=62.38*(1-Tdew/344.56)**(1.246)/1000.

    Re_L                                    = G * Dh / mu_L
    Re_G                                    = G * Dh / mu_G
    Pr_L                                    = cp_L * mu_G / k_L
    Pr_G                                    = cp_G * mu_G / k_G

    h_nb        = 55 * pr** 0.12 * (-log10(pr))**(-0.55) * M**(-0.5) * q_flux**0.67
    h_conv_l    = (3.66 + (0.0668 * Dh / L * Re_L * Pr_L) / (1 + 0.04 * (Dh / L * Re_L * Pr_L)**(2.0/3.0))) * k_L / Dh
    h_conv_g    = (3.66 + (0.0668 * Dh / L * Re_G * Pr_G) / (1 + 0.04 * (Dh / L * Re_G * Pr_G)**(2.0/3.0))) * k_G / Dh
    h_conv_tp   = h_conv_l * (1 - x) + h_conv_g * x
    Co          = sqrt(sig / (g * (rho_L - rho_G) * Dh**2))
    h_TP        = h_nb * (1 - x) + h_conv_tp * (1.0 + 80.0 * (x**2 - x**6) * exp(-0.6 * Co))

    return h_TP


def Bertsch_MC_Average(x_min, x_max, AS, G, Dh, q_flux, L, T_sat_L, T_sat_V):
    """
    Returns the average heat transfer coefficient
    between qualities of x_min and x_max.
    for Bertsch two-phase evaporation in mico-channel HX
    """
    if not x_min == x_max:
        # A proper range is given
        return quad(Bertsch_MC, x_min, x_max, args=(AS, G, Dh, q_flux, L, T_sat_L, T_sat_V))[0] / (x_max - x_min)
    else:
        # A single value is given
        return Bertsch_MC(x_min, AS, G, Dh, q_flux, L, T_sat_L, T_sat_V)


def KM_Cond_Average(x_min, x_max, AS, G, Dh, T_bubble, T_dew, p, beta, C=None, satTransport=None):
    """
    Returns the average pressure gradient and average heat transfer coefficient
    between qualities of x_min and x_max.
    for Kim&Mudawar two-phase condensation in mico-channel HX

    To obtain the pressure gradient for a given value of x, pass it in as x_min and x_max

    Required parameters:
    * x_min : The minimum quality for the range [-]
    * x_max : The maximum quality for the range [-]
    * AS : AbstractState with the refrigerant name and backend
    * G : Mass flux [kg/m^2/s]
    * Dh : Hydraulic diameter of tube [m]
    * Tbubble : Bubblepoint temperature of refrigerant [K]
    * Tdew : Dewpoint temperature of refrigerant [K]
    * beta: channel aspect ratio (=width/height)

    Optional parameters:
    * satTransport : A dictionary with the keys 'mu_f','mu_g,'rho_f','rho_g', 'sigma' for the saturation properties.  So they can be calculated once and passed in for a slight improvement in efficiency
    """

    def KMFunc(x):
        dpdz, h     = Kim_Mudawar_condensing_DPDZ_h(AS, G, Dh, x, T_bubble, T_dew, p, beta, C, satTransport)
        return dpdz , h

    # Use Simpson's Rule to calculate the average pressure gradient
    # Can't use adapative quadrature since function is not sufficiently smooth
    # Not clear why not sufficiently smooth at x>0.9
    if x_min == x_max:
        return KMFunc(x_min)
    else:
        # Calculate the transport properties once
        satTransport    = {}
        AS.update(CP.QT_INPUTS, 0.0, T_bubble)
        satTransport['rho_f']                       = AS.rhomass()          # [kg/m^3]
        satTransport['mu_f']                        = AS.viscosity()        # [Pa-s OR kg/m-s]
        satTransport['cp_f']                        = AS.cpmass()           # [J/kg-K]
        satTransport['k_f']                         = AS.conductivity()     # [W/m-K]
        AS.update(CP.QT_INPUTS, 1.0, T_dew)
        satTransport['rho_g']                       = AS.rhomass()          # [kg/m^3]
        satTransport['mu_g']                        = AS.viscosity()        # [Pa-s OR kg/m-s]

        # Calculate Dp and h over the range of xx
        xx                                          = np.linspace(x_min, x_max, 100)
        DP                                          = np.zeros_like(xx)
        h                                           = np.zeros_like(xx)
        for i in range(len(xx)):
            DP[i]                                   = KMFunc(xx[i])[0]
            h[i]                                    = KMFunc(xx[i])[1]

        # Use Simpson's rule to carry out numerical integration to get average DP and average h
        if abs(x_max - x_min) < 5 * machine_eps:
            # return just one of the edge values
            return -DP[0], h[0]
        else:
            # Use Simpson's rule to carry out numerical integration to get average DP and average h
            return -simps(DP, xx) / (x_max - x_min), simps(h, xx) / (x_max - x_min)


def Kim_Mudawar_condensing_DPDZ_h(AS, G, Dh, x, T_bubble, T_dew, p, beta, C=None, satTransport=None):
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
    x                                       = float(x)

    if satTransport == None:
        # Calculate Necessary saturation properties
        AS.update(CP.QT_INPUTS, 0.0, T_bubble)
        rho_f                               = AS.rhomass()                  # [kg/m^3]
        mu_f                                = AS.viscosity()                # [Pa-s OR kg/m-s]
        cp_f                                = AS.cpmass()                   # [J/kg-K]
        k_f                                 = AS.conductivity()             # [W/m-K]
        AS.update(CP.QT_INPUTS, 1.0, T_dew)
        rho_g                               = AS.rhomass()                  # [kg/m^3]
        mu_g                                = AS.viscosity()                # [Pa-s OR kg/m-s]
    else:
        # Pull out of the dictionary
        rho_f                               = satTransport['rho_f']
        rho_g                               = satTransport['rho_g']
        mu_f                                = satTransport['mu_f']
        mu_g                                = satTransport['mu_g']
        cp_f                                = satTransport['cp_f']
        k_f                                 = satTransport['k_f']

    AS.update(CP.PQ_INPUTS, p, x)
    sigma                                   = AS.surface_tension()          # surface tesnion [N/m]

    Pr_f                                    = cp_f * mu_f / k_f             # [-]

    Re_f                                    = G * (1 - x) * Dh / mu_f
    Re_g                                    = G * x * Dh / mu_g


    if x == 1:                                      # No liquid
        f_f = 0                                     # Just to be ok until next step
    elif Re_f < 2000:                               # Laminar
        f_f                                 = 16.0 / Re_f
        if beta < 1:
            f_f = 24*(1-1.3553*beta + 1.9467*beta*beta - 1.7012*pow(beta, 3) + 0.9564*pow(beta, 4) - 0.2537*pow(beta, 5))/Re_f
    elif Re_f >= 20000:                             # Fully-Turbulent
        f_f                                 = 0.046 * pow(Re_f, -0.2)
    else:                                           # Transient
        f_f                                 = 0.079 * pow(Re_f, -0.25)

    if x == 0:                                      # No gas
        f_g                                 = 0     # Just to be ok until next step
    elif Re_g < 2000:                               # Laminar
        f_g                                 = 16.0 / Re_g
        if beta < 1:
            f_g = 24 * (1 - 1.3553 * beta + 1.9467 * beta*beta - 1.7012 * pow(beta, 3) + 0.9564*pow(beta, 4) - 0.2537 * pow(beta, 5)) / Re_g
    elif Re_g >= 20000:                             # Fully-Turbulent
        f_g                                 = 0.046 * pow(Re_g, -0.2)
    else:                                           # Transient
        f_g                                 = 0.079 * pow(Re_g, -0.25)

    Re_fo                                   = G * Dh / mu_f
    Su_go                                   = rho_g * sigma * Dh / pow(mu_g, 2)

    dpdz_f                                  = 2 * f_f / rho_f * pow(G * (1 - x), 2) / Dh
    dpdz_g                                  = 2 * f_g / rho_g * pow(G * x, 2) / Dh

    if x <= 0:
        # Entirely liquid
        dpdz                                = dpdz_f
        AS.update(CP.QT_INPUTS, 0.0, T_bubble)
        p_sat                               = AS.p()                # pressure [Pa]
        h = f_h_1phase_MicroTube(G, Dh, T_bubble, p_sat, AS, Phase='SatLiq')[1]
        return dpdz, h
    if x >= 1:
        # Entirely vapor
        dpdz                                = dpdz_g
        AS.update(CP.QT_INPUTS, 1.0, T_dew)
        p_sat                               = AS.p()                # pressure [Pa]
        h = f_h_1phase_MicroTube(G, Dh, T_dew, p_sat, AS, Phase='SatVap')[1]
        return dpdz, h

    X                                       = sqrt(dpdz_f/dpdz_g)

    # Find the C coefficient (Calculate C if not passed, otherwise use the set value of C)
    if C == None:
        if Re_f < 2000 and Re_g < 2000:
            C = 3.5e-5 * pow(Re_fo, 0.44) * pow(Su_go, 0.50) * pow(rho_f / rho_g, 0.48)
        elif Re_f < 2000 and Re_g >= 2000:
            C = 0.0015 * pow(Re_fo, 0.59) * pow(Su_go, 0.19) * pow(rho_f / rho_g, 0.36)
        elif Re_f >= 2000 and Re_g < 2000:
            C = 8.7e-4 * pow(Re_fo, 0.17) * pow(Su_go, 0.50) * pow(rho_f / rho_g, 0.14)
        elif Re_f >= 2000 and Re_g >= 2000:
            C = 0.39 * pow(Re_fo, 0.03) * pow(Su_go, 0.10) * pow(rho_f / rho_g, 0.35)
    else:
        pass

    # Two-phase multiplier
    phi_f_square                            = 1.0 + C / X + 1.0 / X**2
    phi_g_square                            = 1.0 + C * X + X**2

    # Find Condensing pressure drop griendient
    if dpdz_g * phi_g_square > dpdz_f * phi_f_square:
        dpdz                                = dpdz_g * phi_g_square
    else:
        dpdz                                = dpdz_f * phi_f_square

    # Use calculated Lockhart-Martinelli parameter
    Xtt                                     = X
    # Simplified Lockhart-Martinelli parameter from Kim & Mudawar (2013) "Universal approach to predict HTC for condensing mini/micro-channel flow"
    # Xtt = pow(mu_f/mu_g,0.1) * pow((1-x)/x,0.9) * pow(rho_g/rho_f,0.5)

    # Modified Weber number
    if Re_f <= 1250:
        We_star     = 2.45 * pow(Re_g, 0.64) / (pow(Su_go, 0.3) * pow(1 + 1.09 * pow(Xtt, 0.039),0.4))
    else:
        We_star     = 0.85 * pow(Re_g, 0.79) * pow(Xtt, 0.157) / (pow(Su_go, 0.3) * pow(1 + 1.09*pow(Xtt, 0.039), 0.4)) \
                      * pow(pow(mu_g / mu_f, 2) * (rho_f/rho_g), 0.084)

    # Condensation Heat transfer coefficient
    if We_star > 7 * Xtt**0.2:                                  # for annual flow (smooth-annular, wavy-annular, transition)
        h = k_f / Dh * 0.048 * pow(Re_f, 0.69) * pow(Pr_f, 0.34) * sqrt(phi_g_square) / Xtt
    else:                                                       # for slug and bubbly flow
        h = k_f / Dh * pow((0.048 * pow(Re_f, 0.69) * pow(Pr_f, 0.34) * sqrt(phi_g_square) / Xtt)**2
                           + (3.2e-7 * pow(Re_f, -0.38) * pow(Su_go, 1.39))**2, 0.5)

    return dpdz, h


def KM_Evap_Average(x_min, x_max, AS, G, Dh, T_bubble, T_dew, p, beta, q_fluxH, PH_PF=1, C=None, satTransport=None):
    """
    Returns the average pressure gradient and average heat transfer coefficient
    between qualities of x_min and x_max.
    for Kim&Mudawar two-phase evaporation in mico-channel HX

    To obtain the pressure gradient for a given value of x, pass it in as x_min and x_max

    Required parameters:
    * x_min : The minimum quality for the range [-]
    * x_max : The maximum quality for the range [-]
    * AS : AbstractState with the refrigerant name and backend
    * G : Mass flux [kg/m^2/s]
    * Dh : Hydraulic diameter of tube [m]
    * Tbubble : Bubblepoint temperature of refrigerant [K]
    * Tdew : Dewpoint temperature of refrigerant [K]
    * p : pressure [Pa]
    * beta: channel aspect ratio (=width/height)
    * q_fluxH: heat flux [W/m^2]
    * PH_PF: ratio of PH over PF where PH: heated perimeter of channel, PF: wetted perimeter of channel

    Optional parameters:
    * satTransport : A dictionary with the keys 'mu_f','mu_g,'rho_f','rho_g', 'sigma' for the saturation properties.  So they can be calculated once and passed in for a slight improvement in efficiency
    """

    def KMFunc(x):
        dpdz, h     = Kim_Mudawar_boiling_DPDZ_h(AS, G, Dh, x, T_bubble, T_dew, p, beta, q_fluxH, PH_PF, C, satTransport)
        return dpdz, h

    # Use Simpson's Rule to calculate the average pressure gradient
    # Can't use adapative quadrature since function is not sufficiently smooth
    # Not clear why not sufficiently smooth at x>0.9
    if x_min == x_max:
        return KMFunc(x_min)
    else:
        # Calculate the tranport properties once
        satTransport    = {}
        AS.update(CP.QT_INPUTS, 0.0, T_bubble)
        satTransport['rho_f']                           = AS.rhomass()              # [kg/m^3]
        satTransport['mu_f']                            = AS.viscosity()            # [Pa-s OR kg/m-s]
        h_f                                             = AS.hmass()                # [J/kg]
        satTransport['cp_f']                            = AS.cpmass()               # [J/kg-K]
        satTransport['k_f']                             = AS.conductivity()         # [W/m-K]
        AS.update(CP.QT_INPUTS, 1.0, T_dew)
        satTransport['rho_g']                           = AS.rhomass()              # [kg/m^3]
        satTransport['mu_g']                            = AS.viscosity()            # [Pa-s OR kg/m-s]
        h_g                                             = AS.hmass()                # [J/kg]
        satTransport['h_fg']                            = h_g - h_f                 # [J/kg]

        # Calculate Dp and h over the range of xx
        xx                                              = np.linspace(x_min, x_max, 100)
        DP                                              = np.zeros_like(xx)
        h                                               = np.zeros_like(xx)
        for i in range(len(xx)):
            DP[i]                                       = KMFunc(xx[i])[0]
            h[i]                                        = KMFunc(xx[i])[1]

        # Use Simpson's rule to carry out numerical integration to get average DP and average h
        if abs(x_max - x_min) < 5 * machine_eps:
            # return just one of the edge values
            return -DP[0], h[0]
        else:
            # Use Simpson's rule to carry out numerical integration to get average DP and average h
            return -simps(DP, xx) / (x_max - x_min), simps(h, xx) / (x_max - x_min)


def Kim_Mudawar_boiling_DPDZ_h(AS, G, Dh, x, T_bubble, T_dew, p, beta, q_fluxH, PH_PF=1, C=None, satTransport=None):
    """
    This function return the pressure gradient and heat transfer coefficient for
    two phase fluid inside Micro-channel tube while BOILING (EVAPORATION)

    Correlations of DPDZ Based on: Kim and Mudawar (2013) "Universal approach to predicting
    two-phase frictional pressure drop for mini/micro-channel saturated flow boiling"

    Correlations of HTC Based on: Kim and Mudawar (2013) "Universal approach to predicting
    saturated flow boiling heat transfer in mini/micro-channels - Part II. Two-heat heat transfer coefficient"
    """
    # Convert the quality, which might come in as a single numpy float value, to a float
    # With the conversion, >20x speedup in the LockhartMartinelli function, not clear why
    x                                           = float(x)

    if satTransport == None:
        # Calculate Necessary saturation properties
        AS.update(CP.QT_INPUTS, 0.0, T_bubble)
        rho_f                                   = AS.rhomass()              # [kg/m^3]
        mu_f                                    = AS.viscosity()            # [Pa-s OR kg/m-s]
        h_f                                     = AS.hmass()                # [J/kg]
        cp_f                                    = AS.cpmass()               # [J/kg-K]
        k_f                                     = AS.conductivity()         # [W/m-K]
        AS.update(CP.QT_INPUTS, 1.0, T_dew)
        rho_g                                   = AS.rhomass()              # [kg/m^3]
        mu_g                                    = AS.viscosity()            # [Pa-s OR kg/m-s]
        h_g                                     = AS.hmass()                # [J/kg]
        h_fg                                    = h_g - h_f                 # [J/kg]
    else:
        # Pull out of the dictionary
        rho_f                                   = satTransport['rho_f']
        rho_g                                   = satTransport['rho_g']
        mu_f                                    = satTransport['mu_f']
        mu_g                                    = satTransport['mu_g']
        h_fg                                    = satTransport['h_fg']
        cp_f                                    = satTransport['cp_f']
        k_f                                     = satTransport['k_f']

    pc                                          = AS.p_critical()           # critical pressure [Pa]
    pr                                          = p / pc                    # reduced pressure [-]
    AS.update(CP.PQ_INPUTS, p, x)
    sigma                                       = AS.surface_tension()      # surface tension [N/m]

    Re_f                                        = G * (1 - x) * Dh / mu_f
    Re_g                                        = G * x * Dh / mu_g
    Pr_f                                        = cp_f * mu_f / k_f

    if x == 1:                                                              # No liquid
        f_f                                     = 0                         # Just to be ok until next step
    elif Re_f < 2000:                                                       # Laminar
        f_f                                     = 16.0 / Re_f
        if beta < 1:
            f_f     = 24 * (1 - 1.3553 * beta + 1.9467 * beta*beta - 1.7012 * pow(beta, 3) + 0.9564 * pow(beta, 4) - 0.2537 * pow(beta, 5)) / Re_f
    elif Re_f >= 20000:                                                     # Fully-Turbulent
        f_f                                     = 0.046 * pow(Re_f, -0.2)
    else:                                                                   # Transient
        f_f                                     = 0.079 * pow(Re_f, -0.25)

    if x == 0:                                                              # No gas
        f_g                                     = 0                         # Just to be ok until next step
    elif Re_g < 2000:                                                       # Laminar
        f_g                                     = 16.0 / Re_g
        if beta < 1:
            f_g = 24 * (1 - 1.3553 * beta + 1.9467 * beta*beta - 1.7012 * pow(beta, 3) + 0.9564 * pow(beta, 4) - 0.2537 * pow(beta, 5)) / Re_g
    elif Re_g >= 20000:                                                     # Fully-Turbulent
        f_g                                     = 0.046 * pow(Re_g, -0.2)
    else:                                                                   # Transient
        f_g                                     = 0.079 * pow(Re_g, -0.25)

    Re_fo                                       = G * Dh / mu_f
    Su_go                                       = rho_g * sigma * Dh / pow(mu_g, 2)

    dpdz_f                                      = 2 * f_f / rho_f * pow(G * (1 - x), 2) / Dh
    dpdz_g                                      = 2 * f_g / rho_g * pow(G * x, 2) / Dh

    if x <= 0:
        # Entirely liquid
        dpdz                                    = dpdz_f
        AS.update(CP.QT_INPUTS, 0.0, T_bubble)
        p_sat                                   = AS.p()                # pressure [Pa]
        h   = f_h_1phase_MicroTube(G, Dh, T_bubble, p_sat, AS, Phase='SatLiq')[1]
        return dpdz, h
    if x >= 1:
        # Entirely vapor
        dpdz                                    = dpdz_g
        AS.update(CP.QT_INPUTS, 1.0, T_dew)
        p_sat                                   = AS.p()                # pressure [Pa]
        h = f_h_1phase_MicroTube(G, Dh, T_dew, p_sat, AS, Phase='SatVap')[1]
        return dpdz, h

    X                                           = sqrt(dpdz_f / dpdz_g)

    We_fo                                       = G * G * Dh / rho_f / sigma
    Bo                                          = q_fluxH / (G * h_fg)

    # Find the C coefficient (Calculate C if not passed, otherwise use the set value of C)
    if C == None:
        # Calculate C (non boiling)
        if Re_f < 2000 and Re_g < 2000:
            Cnon_boiling = 3.5e-5 * pow(Re_fo, 0.44) * pow(Su_go, 0.50) * pow(rho_f / rho_g, 0.48)
        elif Re_f < 2000 and Re_g >= 2000:
            Cnon_boiling = 0.0015 * pow(Re_fo, 0.59) * pow(Su_go, 0.19) * pow(rho_f / rho_g, 0.36)
        elif Re_f >= 2000 and Re_g < 2000:
            Cnon_boiling = 8.7e-4 * pow(Re_fo, 0.17) * pow(Su_go, 0.50) * pow(rho_f / rho_g, 0.14)
        elif Re_f >= 2000 and Re_g >= 2000:
            Cnon_boiling = 0.39 * pow(Re_fo, 0.03) * pow(Su_go, 0.10) * pow(rho_f / rho_g, 0.35)
        # Calculate actual C
        if Re_f >= 2000:
            C       = Cnon_boiling * (1 + 60 * pow(We_fo, 0.32) * pow(Bo * PH_PF, 0.78))
        else:
            C       = Cnon_boiling * (1 + 530 * pow(We_fo, 0.52) * pow(Bo * PH_PF, 1.09))
    else:
        pass

    # Two-phase multiplier
    phi_f_square                                = 1 + C / X + 1 / X**2
    phi_g_square                                = 1 + C * X + X**2

    # Find Boiling pressure drop gradient
    if dpdz_g * phi_g_square > dpdz_f * phi_f_square:
        dpdz                                    = dpdz_g * phi_g_square
    else:
        dpdz                                    = dpdz_f * phi_f_square

    # Use calculated Lockhart-Martinelli parameter
    Xtt = X
    # Simplified X_tt from Kim and Mudawar (2013) "Universal approach to predicting .... Part II. Two-phase heat transfer coefficient"
    # Xtt = pow(mu_f/mu_g,0.1)*pow((1-x)/x,0.9)*pow(rho_g/rho_f,0.5)

    # Pre-dryout saturated flow boiling Heat transfer coefficient
    h_nb    = (2345 * pow(Bo * PH_PF, 0.7) * pow(pr, 0.38) * pow(1 - x, -0.51)) * (0.023 * pow(Re_f, 0.8) * pow(Pr_f, 0.4) * k_f / Dh)
    h_cb    = (5.2 * pow(Bo * PH_PF, 0.08) * pow(We_fo, -0.54) + 3.5*pow(1/Xtt, 0.94) * pow(rho_g / rho_f, 0.25)) * (0.023 * pow(Re_f, 0.8) * pow(Pr_f, 0.4) * k_f / Dh)
    h       = pow(h_nb**2 +h_cb**2, 0.5)

    return dpdz, h


def Cooper_PoolBoiling(p_star, Rp, q, M):
    """
    Cooper M.G., 1984, "Heat flow rates in saturated nucleate boiling - A wide-ranging
    examination using reduced properties. Advances in Heat Transfer. Vol. 16,
    Eds. J.P. Harnett and T.F. Irvine Jr., Academic Press, Orlando, Florida. pp 157-239"

    Rp : surface roughness in microns
    """
    return 55 * p_star**(0.12 - 0.2 * log10(Rp)) * (-log10(p_star))**(-0.55) * q**0.67 * M**(-0.5)


def KandlikarPHE(AS, x_mean, G, D, q, T_bubble, T_dew):
    """
    From http://www.rit.edu/kgcoe/mechanical/taleme/Papers/Conference%20Papers/C041.pdf

    Not recommended for fluids other than R134a
    """
    AS.update(CP.QT_INPUTS, 0.0, T_bubble)
    rhoL                                        = AS.rhomass()              # [kg/m^3]
    mu_f                                        = AS.viscosity()            # [Pa-s OR kg/m-s]
    cp_f                                        = AS.cpmass()               # [J/kg-K]
    k_f                                         = AS.conductivity()         # [W/m/K]
    h_L                                         = AS.hmass()                # [J/kg]

    AS.update(CP.QT_INPUTS, 1.0, T_dew)
    rhoG                                        = AS.rhomass()              # [kg/m^3]
    h_G                                         = AS.hmass()                # [J/kg]

    Pr_f                                        = cp_f * mu_f / k_f         # [-]

    h_LG                                        = h_G-h_L                   # [J/kg]
    alpha_L     = 0.023 * (G * D / mu_f)**0.8 * Pr_f**0.4 * k_f / D         # [W/m^2-K]
    Co          = (rhoG / rhoL)**0.5 * ((1 - x_mean) / x_mean)**0.8
    Bo          = q / (G * h_LG)
#    E_CB=0.512
#    E_NB=0.338
#    F_fl=1.0
    # alpha_r=(2.312*Co**(-0.3)*E_CB+667.3*Bo**(2.8)*F_fl*E_NB)*(1-xmean)**(0.003)*alpha_L
    alpha_r     = 1.055 * (1.056 * Co**(-0.4) + 1.02 * Bo**0.9) * x_mean**(-0.12) * alpha_L**0.98
    return alpha_r


# natural convection for the receiver
def NaturalConv_HTC(AS, HTCat, T_wall, T_inf, P_film, L,D_pipe=None, PlateNum=None):

    """
    Nusselt number for different heat transfer categories;
    find heat transfer coefficient based on Nusselt number;
    characteristic length, A/P;
    Based on: Incropera et al. "Fundamentals of Heat and Mass Transfer"

    Return heat transfer rate by natural convection

    Parameters
    ----------
    AS :            AbstractState with the refrigerant name and backend
    HTCat :         'horizontal_pipe' or 'vertical_pipe' or 'vertical_plate' or 'horizontal_plate'
    T_wall [K]:     surface temperature
    P_film [Pa]:    pressure at the film
    Tinf [K]:       surrounding temperature
    L [m]:          characteristic length
    D_pipe [m]:     pipe diameter
    PlateNum :      'upper_surface' or 'lower_surface'
    ----------
    Return
    h [W/m^2 K]:    Heat transfer coefficient
    """
    # Gravity acceleration
    g                       = 9.81                      # [m/s^2]

    # Film temperature, used to calculate thermal properties
    T_film                  = (T_wall + T_inf) / 2      # [K]

    # thermal expansion coefficient, assumed ideal gas
    beta                    = 1 / T_film                # [1/K]

    # Transport properties calculation film
    AS.update(CP.PT_INPUTS, P_film, T_film)   # use the film temperature to find the outer convective coefficient

    rho_film                = AS.rhomass()              # [kg/m3]
    k_film                  = AS.conductivity()         # [W/m-K]
    mu_film                 = AS.viscosity()            # [Pa-s OR kg/m-s]
    nu_film                 = mu_film/rho_film          # [m^2/s]
    cp_film                 = AS.cpmass()               # [J/kg/K]
    Pr_film                 = cp_film*mu_film/k_film    # [-]

    if T_wall < T_inf:
        # Grashof number, beta-thermal expansion coefficient
        Gr = g * beta * (T_inf - T_wall) * L**3 / nu_film**2    # [-]
    else:
        # Grashof number, beta-thermal expansion coefficient
        Gr = g * beta * (T_wall - T_inf) * L**3 / nu_film**2    # [-]

    # Rayleigh number
    RaL                     = Gr * Pr_film              # [-]

    if HTCat == 'horizontal_pipe':
        RaD                 = RaL * D_pipe**3 / L**3    # [-]

    if RaL < 1e-3:
        Nu                  = 0.0                       # [-]
    else:
        if HTCat == 'vertical_plate':
            if RaL > 1e9:
                Nu = (0.825 + 0.387 * RaL**(1/6) / (1 + (0.492/Pr_film)**(9/16))**(8/27))**2    # [-]
            else:
                Nu = 0.68 + 0.670 * RaL**(1/4) / (1 + (0.492/Pr_film)**(9/16))**(4/9)           # [-]

        elif HTCat == 'horizontal_plate':
            if PlateNum == 'upper_surface':
                if T_wall > T_inf:          # hot plate
                    if RaL >= 1e4 and RaL <= 1e7:
                        Nu = 0.54 * RaL**(1/4)          # [-]
                    elif RaL >= 1e7 and RaL <= 1e11:
                        Nu = 0.15 * RaL**(1/3)          # [-]
                    else:
                        Nu = 0.71 * RaL**(1/4)          # [-]

                else: # cold plate
                    if RaL >= 1e5 and RaL <= 1e10:
                        Nu = 0.27*RaL**(1/4)
                    else:
                        Nu = 0.71*RaL**(1/4) #[-]

            elif PlateNum == 'lower_surface':
                if T_wall > T_inf:              # hot plate
                    if RaL >= 1e5 and RaL <= 1e10:
                        Nu = 0.27 * RaL**(1/4)          # [-]
                    else:
                        Nu = 0.25 * RaL**(1/4)          # [-]

                else: # cold plate
                    if RaL >= 1e4 and RaL <= 1e7:
                        Nu = 0.54 * RaL**(1/4)          # [-]
                    elif RaL >= 1e7 and RaL <= 1e11:
                        Nu = 0.15 * RaL**(1/3)          # [-]
                    else:
                        Nu = 0.25 * RaL**(1/4)          # [-]
            else:
                raise Exception('PlateNum must be either upper_surface or lower_surface')

        elif HTCat == 'horizontal_pipe':
            if RaD <= 1e12:
                # Churchill and Chu, 1975, RaL->RaD
                Nu = (0.60 + 0.387 * RaD**(1/6) / (1 + (0.559/Pr_film)**(9/16))**(8/27))**2   # [-]

            else:                       # Kuehn and Goldstein, 1976.
                temp    = ((0.518 * (RaD**0.25) * (1 + (0.559/Pr_film)**0.6)**(-5/12))**15 + (0.1 * RaD**(1/3))**15)**(1/15)
                Nu      = 2 / (log(1 + 2 / temp))        # [-]

        elif HTCat == 'vertical_pipe':
            if (D_pipe/L) < 35/Gr**(1/4):
                F = 1/3 * ((L/D_pipe) / (1/Gr))**(1/4) + 1      # [-]
            else:
                F = 1.0                                         # [-]
            Nu = F * (0.825 + 0.387 * RaL**(1/6) / (1 + (0.492/Pr_film)**(9/16))**(8/27))**2    # [-]

        else:
            raise Exception('not implemented')

    # Convective heat transfer coefficient
    if HTCat == 'horizontal_pipe':
        h                   = Nu * k_film / D_pipe      # [W/m^2 K]
    else:
        h                   = Nu * k_film / L           # [W/m^2 K]

    return h


# ------------------------------------------------------------------
if __name__ == '__main__':
    # test simps
    x = np.array([1., 2.])
    y = np.array([1., 4.])
    res = simps(y, x)
    print(x, y, res)
