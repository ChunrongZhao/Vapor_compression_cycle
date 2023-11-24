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
def Petterson_supercritical_average(Tout,Tin,T_w,AS,G,OD,ID,D_l,mdot,p,q_flux_w):
    '''
    Petterson et al. (2000), Heat transfer and pressure drop for flow supercritical and subcritical CO2 in microchannel tubes
    All details for this correlation are available in Ding Li Thesis (Appendix B):
    "INVESTIGATION OF AN EJECTOR-EXPANSION DEVICE IN A TRANSCRITICAL CARBON DIOXIDE CYCLE FOR MILITARY ECU APPLICATIONS"
    '''

    def SuperCriticalCondensation_h(T,T_w,AS,G,OD,ID,D_l,mdot,p,q_flux_w):
        '''return h value'''
        return Petterson_supercritical(T,T_w,AS,G,OD,ID,D_l,mdot,p,q_flux_w)[0]
    def SuperCriticalCondensation_f(T,T_w,AS,G,OD,ID,D_l,mdot,p,q_flux_w):
        '''return f value'''
        return Petterson_supercritical(T,T_w,AS,G,OD,ID,D_l,mdot,p,q_flux_w)[1]
    def SuperCriticalCondensation_cp(T,T_w,AS,G,OD,ID,D_l,mdot,p,q_flux_w):
        '''return cp value'''
        return Petterson_supercritical(T,T_w,AS,G,OD,ID,D_l,mdot,p,q_flux_w)[2]
    def SuperCriticalCondensation_rho(T,T_w,AS,G,OD,ID,D_l,mdot,p,q_flux_w):
        '''return rho value'''
        return Petterson_supercritical(T,T_w,AS,G,OD,ID,D_l,mdot,p,q_flux_w)[3]

    if not Tout==Tin:
        #A proper range is given
        h = quad(SuperCriticalCondensation_h,Tin,Tout,args=(T_w,AS,G,OD,ID,D_l,mdot,p,q_flux_w))[0]/(Tout-Tin)
        f = quad(SuperCriticalCondensation_f,Tin,Tout,args=(T_w,AS,G,OD,ID,D_l,mdot,p,q_flux_w))[0]/(Tout-Tin)
        cp = quad(SuperCriticalCondensation_cp,Tin,Tout,args=(T_w,AS,G,OD,ID,D_l,mdot,p,q_flux_w))[0]/(Tout-Tin)
        rho = quad(SuperCriticalCondensation_rho,Tin,Tout,args=(T_w,AS,G,OD,ID,D_l,mdot,p,q_flux_w))[0]/(Tout-Tin)
        return (h,f,cp,rho)
    else:
        #A single value is given
        return Petterson_supercritical(Tout,T_w,AS,G,OD,ID,D_l,mdot,p,q_flux_w)


def Petterson_supercritical(T,T_w,AS,G,OD,ID,D_l,mdot,p,q_flux_w):
    AS.update(CP.PT_INPUTS,p,T_w)
    h_w = AS.hmass() #[J/kg]
    mu_w = AS.viscosity() #[Pa-s OR kg/m-s]
    cp_w = AS.cpmass() #[J/kg-K]
    k_w = AS.conductivity() #[W/m-K]
    rho_w = AS.rhomass() #[kg/m^3]
    Pr_w = cp_w * mu_w / k_w #[-]

    AS.update(CP.PT_INPUTS,p,T)
    h = AS.hmass() #[J/kg]
    mu = AS.viscosity() #[Pa-s OR kg/m-s]
    cp = AS.cpmass() #[J/kg-K]
    k = AS.conductivity() #[W/m-K]
    rho = AS.rhomass() #[kg/m^3]
    Pr = cp * mu / k #[-]

    if mdot == 0: #For the case of Michro-channel
        Dh = OD
        Re=G*Dh/mu
        Re_w=G*Dh/mu_w
    else: #for the case of fin-and-tube
        Dh = OD - ID
        Area=pi*(OD**2-ID**2)/4.0
        u=mdot/(Area*rho)
        Re=rho*u*Dh/mu
        Re_w=Re#rho_w*u*Dh/mu_w

    if G > 350:
        e_D = 0 #smooth pipe
        f = (-1.8*log10(6.9/Re + (1/3.7*e_D)**1.11))**(-2)/4
        Nu_m = (f/8)*(Re-1000)*Pr/(1+12.7*sqrt(f/8)*(Pr**(2/3)-1)) *(1+(D_l)**(2/3))
        Nu = Nu_m * (Pr/Pr_w)**0.11

    else: # G<350

        M = 0.001 #[kg/J]
        K = 0.00041 #[kg/J]

        cp_avg = (h-h_w)/(T-T_w)

        if cp_avg/cp_w <= 1:
            n = 0.66 - K*(q_flux_w/G)
        else: #cp_avg/cp_w <1
            n = 0.9 - K*(q_flux_w/G)

        f0 = (0.79*log(Re)-1.64)**(-2)

        g =9.81
        #coefficient of thermal expansion
        beta = AS.isobaric_expansion_coefficient() #[1/K]
        #Grashoff number
        Gr = g*beta*(T_w-T)*Dh**3/(mu/rho)**2
        if Gr/Re**2 < 5e-4:
            f = f0 * (mu_w/mu)**0.22
        elif  Gr/Re**2 >= 5e-4 and G/Re**2 < 0.3:
            f = 2.15 * f0 * (mu_w/mu)**0.22 * (Gr/Re)**0.1
        else: #use f0 for friction factor
            f = f0

        Nu_w_ppk = (f0/8)*Re_w*Pr_w/(1.07+12.7*sqrt(f/8)*(Pr_w**(2/3)-1))

        Nu = Nu_w_ppk * (1-M*q_flux_w/G) * (cp_avg/cp_w)**n

    h = k*Nu/Dh #[W/m^2-K]

    return (h,f,cp,rho)


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


def f_h_1phase_Channel(mdot,W,H,T,p,AS,Phase):

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

    Dh = 2*H*W/(H+W)
    Area=W*H
    u=mdot/(Area*rho)
    Re=rho*u*(Dh)/mu

    # Friction factor of Churchill (Darcy Friction factor where f_laminar=64/Re)
    e_D = 0
    A = ((-2.457 * log( (7.0 / Re)**(0.9) + 0.27 * e_D)))**16
    B = (37530.0 / Re)**16
    f = 8 * ((8/Re)**12.0 + 1 / (A + B)**(1.5))**(1/12)

    # Heat Transfer coefficient of Gnielinski for high Re,
    if Re>1000:
        Nu = (f / 8) * (Re - 1000.) * Pr / (1 + 12.7 * sqrt(f / 8) * (Pr**(0.66666) - 1)) #[-]
    else:
        Nu = 3.66
    h = k * Nu / Dh #W/m^2-K
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


# ------------------------------------------------------------------
if __name__ == '__main__':
    # test simps
    x = np.array([1., 2.])
    y = np.array([1., 4.])
    res = simps(y, x)
    print(x, y, res)
