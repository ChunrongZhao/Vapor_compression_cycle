from __future__                                             import division, print_function, absolute_import
from math                                                   import log,exp
from CoolProp.CoolProp                                      import HAPropsSI, cair_sat
from ACHP_codes.Correlations.FinStructure_Correlations      import WavyLouveredFins, HerringboneFins, PlainFins
from ACHP_codes.Correlations.MicroFinCorrelations           import MultiLouveredMicroFins
import numpy                                                as np
# ----------------------------------------------------------------------
class DWSVals():
    """
    Empty Class for passing data with DryWetSegment
    """
    def __init__(self):
        # Don't do anything
        pass


def DryWetSegment(DWS):
    """
    Generic solver function for dry-wet mixed surface conditions for a given element.
    Can handle superheated, subcooled and two-phase regions.
    Does not handle the pressure drops, only HT required to get the dry/wet interface
    """

    # List of required parameters
    RequiredParameters = ['T_in_a', 'h_a', 'cp_da', 'eta_a', 'A_a', 'p_in_a', 'RH_in_a', 'T_in_r', 'p_in_r',
                          'h_r', 'cp_r', 'A_r', 'R_w', 'm_dot_r', 'Fins', 'FinsType']

    # Check that all the parameters are included, raise exception otherwise
    for param in RequiredParameters:
        if not hasattr(DWS, param):
            raise AttributeError("Parameter " + param + " is required for DWS class in DryWetSegment")

    # Retrieve values from structures defined above
    T_in_a                      = DWS.T_in_a
    if DWS.h_a < 0.000000001:
        print("Warning: Dws.h_a was constrained to 0.001, original value: ", DWS.h_a)
        h_a                     = 0.000000001
    else:
        h_a                     = DWS.h_a

    cp_da                       = DWS.cp_da
    eta_a                       = DWS.eta_a                 # from fin correlations, overall airside surface effectiveness
    A_a                         = DWS.A_a
    p_in_a                      = DWS.p_in_a
    RH_in_a                     = DWS.RH_in_a
    m_dot_da                    = DWS.m_dot_da

    T_in_r                      = DWS.T_in_r
    p_in_r                      = DWS.p_in_r

    if DWS.h_r < 0.000000001:
        print("Warning: Dws.h_r was constrained to 0.001, original value: ", DWS.h_r)
        h_r                     = 0.000000001
    else:
        h_r                     = DWS.h_r

    cp_r                        = DWS.cp_r
    A_r                         = DWS.A_r
    m_dot_r                     = DWS.m_dot_r
    R_w                         = DWS.R_w

    # Calculate the dewpoint (amongst others); http://www.coolprop.org/fluid_properties/HumidAir.html
    # kg water/kg dry air; Humidity Ratio
    omega_in                    = HAPropsSI('W', 'T', T_in_a, 'P', p_in_a, 'R', RH_in_a)    # Relative humidity in [0, 1]
    T_dp                        = HAPropsSI('D', 'T', T_in_a, 'P', p_in_a, 'W', omega_in)   # dew-point temperature
    h_in_a                      = HAPropsSI('H', 'T', T_in_a, 'P', p_in_a, 'W', omega_in)   # Mixture enthalpy per dry air [J/kg_da]

    # Internal UA between fluid flow and outside surface (neglecting tube conduction)
    UA_i                        = h_r * A_r                                     # [W/K], from Shah or f_h_1phase_Tube-fct -> Correlations.py
    # External UA between wall and free stream
    UA_o                        = eta_a * h_a * A_a                             # [W/K], from fin correlations
    # wall UA
    UA_w                        = 1 / R_w
    # Internal Ntu
    Ntu_i                       = UA_i / (m_dot_r * cp_r)                       # [-]
    # External Ntu (multiplied by eta_a since surface is finned and has lower effectiveness)
    Ntu_o                       = eta_a * h_a * A_a / (m_dot_da * cp_da)        # [-]

    # --------------------------------------------------------------------------------------
    # (Two-Phase analysis)
    if DWS.IsTwoPhase:
        UA                      = 1 / (1 / (h_a * A_a * eta_a) + 1 / (h_r * A_r) + 1 / UA_w)    # overall heat transfer coefficient
        Ntu_dry                 = UA / (m_dot_da * cp_da)                                       # Number of transfer units
        epsilon_dry             = 1 - exp(-Ntu_dry)             # since Cr=0, e.g. see Incropera - Fundamentals of Heat and Mass Transfer, 2007, p. 690
        Q_dry                   = epsilon_dry * m_dot_da * cp_da * (T_in_a - T_in_r)
        T_out_a                 = T_in_a - Q_dry / (m_dot_da * cp_da)                           # outlet temperature, dry fin

        T_so_a                  = (UA_o * T_in_a + UA_i * T_in_r) / (UA_o + UA_i)               # inlet surface temperature (neglect wall thermal conductance)
        T_so_b                  = (UA_o * T_out_a+ UA_i * T_in_r) / (UA_o + UA_i)               # outlet surface temperature (neglect wall thermal conductance)

        if T_so_b > T_dp:
            # All dry, since surface at outlet dry
            f_dry               = 1.0
            Q                   = Q_dry                                                         # [W]
            Q_sensible          = Q                                                             # [W]
            h_out_a             = h_in_a - Q / m_dot_da                                         # [J/kg_da]
            # Air outlet humidity ratio
            DWS.omega_out       = omega_in                                                      # [kg/kg]
        else:
            if T_so_a < T_dp:
                # All wet, since surface at inlet wet
                f_dry           = 0.0
                Q_dry           = 0.0
                T_ac            = T_in_a                                                        # temp at onset of wetted wall
                h_ac            = h_in_a                                                        # enthalpy at onset of wetted surface
            else:
                # Partially wet and dry (i.e T_so_b < Tdp < T_so_a)
                # Air temperature at the interface between wet and dry surface
                # Based on equating heat fluxes at the wall which is at dew point UA_i*(Tw-Ti)=UA_o*(To-Tw)
                T_ac            = T_dp + UA_i / UA_o * (T_dp - T_in_r)
                # Dry effectiveness (minimum capacitance on the air side by definition)
                epsilon_dry     = (T_in_a - T_ac) / (T_in_a - T_in_r)
                # Dry fraction found by solving epsilon=1-exp(-f_dry*Ntu) for known epsilon from above equation
                f_dry           = -1.0 / Ntu_dry * log(1.0 - epsilon_dry)
                # Enthalpy, using air humidity at the interface between wet and dry surfaces, which is same humidity ratio as inlet
                h_ac            = HAPropsSI('H', 'T', T_ac, 'P', p_in_a, 'W', omega_in)         # [J/kg_da]
                # Dry heat transfer
                Q_dry           = m_dot_da * cp_da * (T_in_a - T_ac)

            # Saturation specific heat at mean water temp (c_s : partial derivative dh_sat/dT @ T_sat_r)
            c_s                 = cair_sat(T_in_r) * 1000                                       # [J/kg-K]
            # Find new, effective fin efficiency since cs/cp is changed from wetting
            # Ratio of specific heats [-]
            DWS.Fins.Air.cs_cp  = c_s / cp_da
            DWS.Fins.WetDry     = 'Wet'

            # Compute the fin efficiency based on the user choice of FinsType
            if DWS.FinsType == 'WavyLouveredFins':
                WavyLouveredFins(DWS.Fins)
            elif DWS.FinsType == 'HerringboneFins':
                HerringboneFins(DWS.Fins)
            elif DWS.FinsType == 'PlainFins':
                PlainFins(DWS.Fins)
            elif DWS.FinsType == 'MultiLouveredMicroFins':
                MultiLouveredMicroFins(DWS.Fins)

            eta_a_wet           = DWS.Fins.eta_a_wet
            UA_o                = eta_a_wet * h_a * A_a
            Ntu_o               = eta_a_wet * h_a * A_a / (m_dot_da * cp_da)

            # Wet analysis overall Ntu for two-phase refrigerant
            # Minimum capacitance rate is by definition on the air side
            # Ntu_wet is the NTU if the entire two-phase region were to be wetted
            UA_wet              = 1 / (c_s / UA_i + cp_da / UA_o + c_s / UA_w)
            Ntu_wet             = UA_wet / m_dot_da
            # Wet effectiveness [-]
            epsilon_wet         = 1 - exp(-(1 - f_dry) * Ntu_wet)
            # Air saturated at refrigerant saturation temp [J/kg]
            h_s_s_o             = HAPropsSI('H', 'T', T_in_r, 'P', p_in_a, 'R', 1.0)    # [kJ/kg_da]

            # Wet heat transfer [W]
            Q_wet               = epsilon_wet * m_dot_da * (h_ac - h_s_s_o)
            # Total heat transfer [W]
            Q                   = Q_wet + Q_dry
            # Air exit enthalpy [J/kg]
            h_out_a             = h_ac - Q_wet / m_dot_da
            # Saturated air temp at effective surface temp [J/kg_da]
            h_s_s_e             = h_ac - (h_ac - h_out_a) / (1 - exp(-(1 - f_dry) * Ntu_o))
            # Effective surface temperature [K]
            T_s_e               = HAPropsSI('T', 'H', h_s_s_e, 'P', p_in_a, 'R', 1.0)
            # Outlet dry-bulb temp [K]
            T_out_a             = T_s_e + (T_ac - T_s_e) * exp(-(1 - f_dry) * Ntu_o)
            # Sensible heat transfer rate [kW]
            Q_sensible          = m_dot_da * cp_da * (T_in_a - T_out_a)
        # Outlet is saturated vapor
        T_out_r                  = DWS.T_dew_r

    # (Single-Phase analysis)
    else:
        # Overall UA
        UA                      = 1 / (1 / UA_i + 1 / UA_o + 1 / UA_w)
        # Min and max capacitance rates [W/K]
        C_min                   = min([cp_r * m_dot_r, cp_da * m_dot_da])
        C_max                   = max([cp_r * m_dot_r, cp_da * m_dot_da])
        # Capacitance rate ratio [-]
        C_star                  = C_min / C_max
        # Ntu overall [-]
        Ntu_dry                 = UA / C_min

        if Ntu_dry < 0.0000001:
            print("warning:  NTU_dry in dry wet segment was negative. forced it to positive value of 0.001!")
            Ntu_dry             = 0.0000001

        # Counterflow effectiveness [-]
        # epsilon_dry = ((1 - exp(-Ntu_dry * (1 - C_star))) /
        #   (1 - C_star * exp(-Ntu_dry * (1 - C_star))))

        # Crossflow effectiveness (e.g. see Incropera - Fundamentals of Heat and Mass Transfer, 2007, p. 662)
        if (cp_r * m_dot_r) < (cp_da * m_dot_da):
            epsilon_dry         = 1 - exp(-C_star**(-1) * (1 - exp(-C_star * Ntu_dry)))
            # Cross flow, single phase, cmax is airside, which is unmixed
        else:
            epsilon_dry         = (1 / C_star) * (1 - exp(-C_star * (1 - exp(-Ntu_dry))))
            # Cross flow, single phase, cmax is refrigerant side, which is mixed

        # Dry heat transfer [W]
        Q_dry                   = epsilon_dry * C_min * (T_in_a - T_in_r)
        # Dry-analysis air outlet temp [K]
        T_out_a_dry             = T_in_a - Q_dry / (m_dot_da * cp_da)
        # Dry-analysis outlet temp [K]
        T_out_r                 = T_in_r + Q_dry / (m_dot_r * cp_r)
        # Dry-analysis air outlet enthalpy from energy balance [J/kg]
        h_out_a                 = h_in_a - Q_dry / m_dot_da
        # Dry-analysis surface outlet temp [K] (neglect wall thermal conductance)
        T_out_s                 = (UA_o * T_out_a_dry + UA_i * T_in_r) / (UA_o + UA_i)
        # Dry-analysis surface inlet temp [K] (neglect wall thermal conductance)
        T_in_s                  = (UA_o * T_in_a + UA_i * T_out_r) / (UA_o + UA_i)
        # Dry-analysis outlet refrigerant temp [K]
        T_out_r_dry             = T_out_r
        # Dry fraction [-]
        f_dry                   = 1.0
        # Air outlet humidity ratio [-]
        DWS.omega_out           = omega_in

        # If inlet surface temp below dewpoint, whole surface is wetted
        if T_in_s < T_dp:
            isFullyWet          = True
        else:
            isFullyWet          = False

        if T_out_s < T_dp or isFullyWet:
            # There is some wetting, either the coil is fully wetted or partially wetted
            # Loop to get the correct c_s
            # Start with the inlet temp as the outlet temp
            x1                  = T_in_r + 1                            # Lowest possible outlet temperature
            x2                  = T_in_a - 1                            # Highest possible outlet temperature
            eps                 = 1e-8
            iter = 1
            change = 999
            while (iter <= 3 or change > eps) and iter < 100:
                if iter == 1:
                    T_out_r     = x1
                if iter > 1:
                    T_out_r     = x2

                T_out_r_start   = T_out_r
                # Saturated air enthalpy at the inlet water temperature [J/kg]
                h_s_w_i         = HAPropsSI('H', 'T', T_in_r, 'P', p_in_a, 'R', 1.0)            # [J/kg_da]
                # Saturation specific heat at mean water temp [J/kg]
                c_s             = cair_sat((T_in_r + T_out_r) / 2.0) * 1000
                # Ratio of specific heats [-]
                DWS.Fins.Air.cs_cp  = c_s / cp_da
                # Find new, effective fin efficiency since cs/cp is changed from wetting
                # Based on the user choice of FinsType
                if DWS.FinsType == 'WavyLouveredFins':
                    WavyLouveredFins(DWS.Fins)
                elif DWS.FinsType == 'HerringboneFins':
                    HerringboneFins(DWS.Fins)
                elif DWS.FinsType == 'PlainFins':
                    PlainFins(DWS.Fins)
                elif DWS.FinsType == 'MultiLouveredMicroFins':
                    MultiLouveredMicroFins(DWS.Fins)

                # Effective humid air mass flow ratio
                m_star          = m_dot_da / (m_dot_r * (cp_r / c_s))

                # compute the new Ntu_owet
                Ntu_owet        = eta_a * h_a * A_a / (m_dot_da * cp_da)
                m_star          = min([cp_r * m_dot_r / c_s, m_dot_da]) / max([cp_r * m_dot_r / c_s, m_dot_da])
                m_dot_min       = min([cp_r * m_dot_r / c_s, m_dot_da])
                # Wet-analysis overall Ntu [-] (neglect wall thermal conductance)
                Ntu_wet         = Ntu_o / (1 + m_star * (Ntu_owet / Ntu_i))
                if cp_r * m_dot_r > c_s * m_dot_da:
                    Ntu_wet     = Ntu_o / (1 + m_star * (Ntu_owet / Ntu_i))
                else:
                    Ntu_wet     = Ntu_i / (1 + m_star * (Ntu_i / Ntu_owet))

                # Counterflow effectiveness for wet analysis
                epsilon_wet     = ((1 - exp(-Ntu_wet * (1 - m_star))) / (1 - m_star * exp(-Ntu_wet * (1 - m_star))))
                # Wet-analysis heat transfer rate
                Q_wet           = epsilon_wet * m_dot_min * (h_in_a - h_s_w_i)
                # Air outlet enthalpy [J/kg_da]
                h_out_a         = h_in_a - Q_wet / m_dot_da
                # Water outlet temp [K]
                T_out_r         = T_in_r + m_dot_da / (m_dot_r * cp_r) * (h_in_a - h_out_a)
                # Water outlet saturated surface enthalpy [J/kg_da]
                h_s_w_o         = HAPropsSI('H', 'T', T_out_r, 'P', p_in_a, 'R', 1.0)               # [J/kg_da]
                # Local UA* and c_s
                UA_star         = 1 / (cp_da / (eta_a * h_a * A_a) + cair_sat((T_in_a + T_out_r) / 2.0) * 1000 * (1 / (h_r * A_r) + 1 / UA_w))
                # Wet-analysis surface temperature [K]
                T_in_s          = T_out_r + UA_star / (h_r * A_r) * (h_in_a - h_s_w_o)
                # Wet-analysis saturation enthalpy [J/kg_da]
                h_s_s_e         = h_in_a + (h_out_a - h_in_a) / (1 - exp(-Ntu_owet))
                # Surface effective temperature [K]
                T_s_e           = HAPropsSI('T', 'H', h_s_s_e, 'P', p_in_a, 'R', 1.0)
                # Air outlet temp based on effective temp [K]
                T_out_a         = T_s_e + (T_in_a - T_s_e) * exp(-Ntu_o)
                # Sensible heat transfer rate [W]
                Q_sensible      = m_dot_da * cp_da * (T_in_a - T_out_a)
                # Error between guess and recalculated value [K]
                error_T_out_r   = T_out_r - T_out_r_start

                if iter > 500:
                    print("Superheated region wet analysis T_out_r convergence failed")
                    DWS.Q       = Q_dry
                    return
                if iter == 1:
                    y1          = error_T_out_r
                if iter > 1:
                    y2          = error_T_out_r
                    x3          = x2 - y2 / (y2 - y1) * (x2 - x1)
                    change      = abs(y2 / (y2 - y1) * (x2 - x1))
                    y1          = y2
                    x1          = x2
                    x2          = x3
                if hasattr(DWS, 'Verbosity') and DWS.Verbosity > 7:
                    print("Full wet iter %d T_out_r %0.5f dT %g" % (iter, T_out_r, error_T_out_r))
                # Update loop counter
                iter            += 1

            # Fully wetted outlet temperature [K]
            T_out_r_wet         = T_out_r
            # Dry fraction
            f_dry               = 0.0

            if T_in_s > T_dp and not isFullyWet:

                # Partially wet and partially dry with single-phase on refrigerant side
                """
                -----------------------------------------------------------
                                            |
                * T_out_a   <----            * T_a,x                 <---- * T_in_a
                                            |
                 ____________Wet____________|              Dry 
                ----------------------------o T_dp ------------------------
                                            |
                * T_in_r    ---->            * T_w,x                 ----> * T_out_r               
                                            |
                -----------------------------------------------------------
                """

                iter            = 1

                # Now do an iterative solver to find the fraction of the coil that is wetted
                x1              = 0.0001
                x2              = 0.9999
                eps             = 1e-8
                while (iter <= 3 or error > eps) and iter < 100:
                    if iter == 1:
                        f_dry   = x1
                    if iter > 1:
                        f_dry   = x2

                    K                       = Ntu_dry * (1.0 - C_star)
                    expk                    = exp(-K * f_dry)
                    if cp_da * m_dot_da < cp_r * m_dot_r:
                        T_out_r_guess       = (T_dp + C_star * (T_in_a - T_dp) - expk * (1 - K / Ntu_o) * T_in_a) / (1 - expk * (1 - K / Ntu_o))
                    else:
                        T_out_r_guess       = (expk * (T_in_a + (C_star - 1) * T_dp) - C_star * (1 + K / Ntu_o) * T_in_a) / (expk * C_star - C_star * (1 + K / Ntu_o))

                    # Wet and dry effective effectiveness
                    epsilon_dry             = ((1 - exp(-f_dry * Ntu_dry * (1 - C_star))) / (1 - C_star * exp(-f_dry * Ntu_dry * (1 - C_star))))
                    epsilon_wet             = ((1 - exp(-(1 - f_dry) * Ntu_wet * (1 - m_star))) / (1 - m_star * exp(-(1 - f_dry) * Ntu_wet * (1 - m_star))))

                    # Temperature of "water" where condensation begins
                    T_w_x       = (T_in_r + m_dot_min / (cp_r * m_dot_r) * epsilon_wet * (h_in_a - h_s_w_i - epsilon_dry * C_min / m_dot_da * T_in_a)) /\
                                  (1 - (C_min * m_dot_min) / (cp_r * m_dot_r * m_dot_da) * epsilon_wet * epsilon_dry)
                    # Temperature of air where condensation begins [K]
                    # Obtained from energy balance on air side
                    T_a_x       = T_in_a - epsilon_dry * C_min * (T_in_a - T_w_x) / (m_dot_da * cp_da)
                    # Enthalpy of air where condensation begins
                    h_a_x       = h_in_a - cp_da * (T_in_a - T_a_x)
                    # New "water" temperature (stored temporarily to be able to build change
                    T_out_r     = C_min / (cp_r * m_dot_r) * epsilon_dry * T_in_a + (1 - C_min / (cp_r * m_dot_r) * epsilon_dry) * T_w_x
                    # Difference between initial guess and outlet
                    error       = T_out_r - T_out_r_guess

                    if iter > 500:
                        print("Superheated region wet analysis f_dry convergence failed")
                        DWS.Q               = Q_dry
                        return
                    if iter == 1:
                        y1                  = error
                    if iter > 1:
                        y2                  = error
                        x3                  = x2 - y2 / (y2 - y1) * (x2 - x1)
                        change              = abs(y2 / (y2 - y1) * (x2 - x1))
                        y1                  = y2
                        x1                  = x2
                        x2                  = x3
                    if hasattr(DWS, 'Verbosity') and DWS.Verbosity > 7:
                        print("Part-wet iter %d T_out_r_guess %0.5f diff %g f_dry: %g" %(iter, T_out_r_guess, error, f_dry))
                    # Update loop counter
                    iter                    += 1

                # Wet-analysis saturation enthalpy [J/kg]
                h_s_s_e                     = h_a_x + (h_out_a - h_a_x) / (1 - exp(-(1 - f_dry) * Ntu_owet))
                # Surface effective temperature [K]
                T_s_e                       = HAPropsSI('T', 'H', h_s_s_e, 'P', p_in_a, 'R', 1.0)
                # Air outlet temp based on effective surface temp [K]
                T_out_a                     = T_s_e + (T_a_x - T_s_e) * exp(-(1 - f_dry) * Ntu_o)
                # Heat transferred [W]
                Q                           = m_dot_r * cp_r * (T_out_r - T_in_r)
                # Dry-analysis air outlet enthalpy from energy balance [J/kg]
                h_out_a                     = h_in_a - Q / m_dot_da
                # Sensible heat transfer rate [kW]
                Q_sensible                  = m_dot_da * cp_da * (T_in_a - T_out_a)
            else:
                Q                           = Q_wet
        else:
            # Coil is fully dry
            T_out_a                         = T_out_a_dry
            Q                               = Q_dry
            Q_sensible                      = Q_dry

    DWS.f_dry                               = f_dry     # f_fry = 0, fully wet
    DWS.omega_out                           = HAPropsSI('W', 'T', T_out_a, 'P', 101325, 'H', h_out_a)
    DWS.RH_out_a                            = HAPropsSI('R', 'T', T_out_a, 'P', 101325, 'W', DWS.omega_out)
    DWS.T_out_a                             = T_out_a
    DWS.Q                                   = Q
    DWS.Q_sensible                          = Q_sensible

    DWS.h_out_a                             = h_out_a
    DWS.h_in_a                              = h_in_a
    DWS.T_out_r                             = T_out_r
    DWS.T_wall_s                            = T_out_r - Q / UA_i            # inner wall temperature for gas cooler model
