from __future__                         import division, print_function, absolute_import
from CoolProp.CoolProp                  import HAPropsSI
from scipy.optimize                     import fsolve, minimize
from ACHP_codes                         import Correlations
from math                               import pi
from ACHP_codes.ACHP_Tools.Solvers      import MultiDimNewtRaph
import numpy                            as np
import CoolProp                         as CP


# ------------------------------------------------------------------------------
def DXPreconditioner(Cycle, epsilon=0.96):
    # Assume the heat exchangers are highly effective
    # Condensing heat transfer rate from enthalpies
    rho_air                             = 1.1
    Cp_air                              = 1005                      # [J/kg/K]

    # AbstractState
    AS                                  = Cycle.AS

    def OBJECTIVE(x):
        T_evap                          = x[0]
        T_cond                          = x[1]
        DT_sh                           = abs(x[2])

        # Use fixed effectiveness to get a guess for the condenser capacity
        Q_cond      = epsilon * Cycle.Condenser.Fins.Air.V_dot_ha * rho_air * (Cycle.Condenser.Fins.Air.T_db - T_cond) * Cp_air

        AS.update(CP.QT_INPUTS, 1.0, T_evap)
        p_evap                          = AS.p()                    # [pa]
        AS.update(CP.QT_INPUTS, 1.0, T_cond)
        p_cond                          = AS.p()                    # [pa]
        Cycle.Compressor.p_in_r         = p_evap
        Cycle.Compressor.p_out_r        = p_cond
        Cycle.Compressor.T_in_r         = T_evap + DT_sh
        Cycle.Compressor.Ref            = Cycle.Ref
        Cycle.Compressor.Oil            = Cycle.Oil
        Cycle.Compressor.AS             = Cycle.AS
        Cycle.Compressor.Calculate()
        W                               = Cycle.Compressor.W

        # Evaporator fully-dry analysis
        Q_evap_dry      = epsilon * Cycle.Evaporator.Fins.Air.V_dot_ha * rho_air * (Cycle.Evaporator.Fins.Air.T_db - T_evap) * Cp_air

        # Air-side heat transfer UA
        Evap                            = Cycle.Evaporator
        Evap.m_dot_r                    = Cycle.Compressor.m_dot_r
        Evap.p_sat_r                    = p_evap
        Evap.AS                         = Cycle.AS
        Evap.Initialize()

        UA_a                            = Evap.Fins.h_a * Evap.Fins.A_a * Evap.Fins.eta_a
        T_in_a                          = Evap.Fins.Air.T_db
        T_out_a                         = T_in_a + Q_evap_dry / (Evap.Fins.m_dot_da * Evap.Fins.cp_da)
        # Refrigerant-side heat transfer UA
        UA_r    = Evap.A_r_wetted * Correlations.ShahEvaporation_Average(0.5, 0.5, AS, Evap.G_r, Evap.ID, Evap.p_sat_r,
                                            Q_evap_dry/Evap.A_r_wetted, Evap.T_bubble_r, Evap.T_dew_r)
        # Get wall temperatures at inlet and outlet from energy balance
        T_so_a                          = (UA_a * Evap.T_in_a + UA_r * T_evap) / (UA_a + UA_r)
        T_so_b                          = (UA_a * T_out_a + UA_r * T_evap) / (UA_a + UA_r)

        T_dewpoint  = HAPropsSI('D', 'T', Cycle.Evaporator.Fins.Air.T_db, 'P', 101325, 'R', Evap.Fins.Air.RH)

        # Now calculate the fully-wet analysis
        # Evaporator is bounded by saturated air at the refrigerant temperature.
        h_ai        = HAPropsSI('H', 'T', Cycle.Evaporator.Fins.Air.T_db, 'P', 101325, 'R', Cycle.Evaporator.Fins.Air.RH)   # [J/kg_da]
        h_s_w_o     = HAPropsSI('H', 'T', T_evap, 'P', 101325, 'R', 1.0)        # [J/kg_da]
        Q_evap_wet  = epsilon * Cycle.Evaporator.Fins.Air.V_dot_ha * rho_air * (h_ai - h_s_w_o)

        # Coil is either fully-wet, fully-dry or partially wet, partially dry
        if T_so_a > T_dewpoint and T_so_b > T_dewpoint:
            # Fully dry, use dry Q
            f_dry                       = 1.0
        elif T_so_a < T_dewpoint and T_so_b < T_dewpoint:
            # Fully wet, use wet Q
            f_dry                       = 0.0
        else:
            f_dry                       = 1 - (T_dewpoint - T_so_a) / (T_so_b - T_so_a)
        Q_evap      = f_dry * Q_evap_dry + (1 - f_dry) * Q_evap_wet

        if Cycle.ImposedVariable == 'Subcooling':               # if Subcooling impose
            AS.update(CP.PT_INPUTS, p_cond, T_cond - Cycle.DT_sc_target)
            h_target                    = AS.hmass()            # [J/kg]
            Q_cond_enthalpy             = Cycle.Compressor.m_dot_r * (Cycle.Compressor.h_out_r - h_target)
        else:                                                   # otherwise, if Charge impose
            AS.update(CP.PT_INPUTS, p_cond, T_cond - 5)
            h_target                    = AS.hmass()            # [J/kg]
            Q_cond_enthalpy             = Cycle.Compressor.m_dot_r * (Cycle.Compressor.h_out_r - h_target)

        Q_evap_enthalpy                 = Cycle.Compressor.m_dot_r * (Cycle.Compressor.h_in_r - h_target)

        resids      = [Q_evap + W + Q_cond, Q_cond + Q_cond_enthalpy, Q_evap - Q_evap_enthalpy]     # [f_dry]
        return resids

    T_evap_init                         = Cycle.Evaporator.Fins.Air.T_db - 15
    T_cond_init                         = Cycle.Condenser.Fins.Air.T_db + 8
    DT_sh_init                          = Cycle.Evaporator.DT_sh
    x               = fsolve(OBJECTIVE, [T_evap_init, T_cond_init, DT_sh_init])
    DT_evap                             = Cycle.Evaporator.Fins.Air.T_db - x[0]
    DT_cond                             = x[1] - Cycle.Condenser.Fins.Air.T_db
    DT_sh                               = abs(x[2])

    return DT_evap-1, DT_cond+1, DT_sh


def SecondaryLoopPreconditioner(Cycle, epsilon=0.9):
    rho_air                             = 1.1                   # [kg/m^2]
    Cp_air                              = 1005                  # [J/kg-K]

    # AbstractState
    AS                                  = Cycle.AS
    # AbstractState for SecLoopFluid
    AS_SLF                              = Cycle.AS_SLF

    def OBJECTIVE(x):
        T_evap                          = x[0]
        T_cond                          = x[1]
        T_in_CC                         = x[2]
        if Cycle.Mode == 'AC':
            # Condenser heat transfer rate
            Q_cond      = epsilon * Cycle.Condenser.Fins.Air.V_dot_ha * rho_air * (Cycle.Condenser.Fins.Air.T_db - T_cond) * Cp_air

            # Compressor power
            AS.update(CP.QT_INPUTS, 1.0, T_evap)
            p_evap                      = AS.p()                # [Pa]
            AS.update(CP.QT_INPUTS, 1.0, T_cond)
            p_cond                      = AS.p()                # [Pa]
            Cycle.Compressor.p_in_r     = p_evap
            Cycle.Compressor.p_out_r    = p_cond
            Cycle.Compressor.T_in_r     = T_evap + Cycle.Compressor.DT_sh
            Cycle.Compressor.Ref        = Cycle.Ref
            Cycle.Compressor.AS         = Cycle.AS
            Cycle.Compressor.Calculate()
            W                           = Cycle.Compressor.W

            Q_coolingcoil_dry   = epsilon * Cycle.CoolingCoil.Fins.Air.V_dot_ha * rho_air * (Cycle.CoolingCoil.Fins.Air.T_db - T_in_CC) * Cp_air

            # Air-side heat transfer UA
            CC                          = Cycle.CoolingCoil
            # CC.AS_g=AS_SLF
            CC.Initialize()
            UA_a                        = CC.Fins.h_a * CC.Fins.A_a * CC.Fins.eta_a
            # Air outlet temp from dry analysis
            T_out_a                     = CC.T_in_a - Q_coolingcoil_dry / (CC.Fins.m_dot_a * CC.Fins.cp_a)

            # Refrigerant side UA
            f, h, Re    = Correlations.f_h_1phase_Tube(Cycle.Pump.m_dot_g/CC.N_circuits, CC.ID, T_in_CC, CC.p_in_g, CC.AS_g)
            UA_r                        = CC.A_g_wetted * h
            # Glycol specific heat
            AS_SLF.update(CP.PT_INPUTS, Cycle.Pump.p_in_g, T_in_CC)
            cp_g                        = AS_SLF.cpmass()                   # [J/kg-K]
            # Refrigerant outlet temp
            T_out_CC                    = T_in_CC + Q_coolingcoil_dry / (Cycle.Pump.m_dot_g * cp_g)

            # Get wall temperatures at inlet and outlet from energy balance
            T_so_a                      = (UA_a * CC.T_in_a + UA_r * T_out_CC) / (UA_a + UA_r)
            T_so_b                      = (UA_a * T_out_a + UA_r * T_in_CC) / (UA_a + UA_r)

            T_dewpoint=HAPropsSI('D','T',CC.Fins.Air.T_db,'P',101325,'R',CC.Fins.Air.RH)
            # Now calculate the fully-wet analysis
            # Evaporator is bounded by saturated air at the refrigerant temperature.
            h_ai                        = HAPropsSI('H', 'T', CC.Fins.Air.T_db, 'P', 101325, 'R', CC.Fins.Air.RH)
            h_s_w_o                     = HAPropsSI('H', 'T', T_in_CC, 'P', 101325, 'R', 1.0)
            Q_coolingcoil_wet           = epsilon * CC.Fins.Air.V_dot_ha * rho_air * (h_ai - h_s_w_o)

            # Coil is either fully-wet, fully-dry or partially wet, partially dry
            if T_so_a > T_dewpoint and T_so_b > T_dewpoint:
                # Fully dry, use dry Q
                f_dry                   = 1.0
            elif T_so_a < T_dewpoint and T_so_b < T_dewpoint:
                # Fully wet, use wet Q
                f_dry                   = 0.0
            else:
                f_dry                   = 1 - (T_dewpoint - T_so_a) / (T_so_b - T_so_a)
            Q_coolingcoil               = f_dry * Q_coolingcoil_dry + (1 - f_dry) * Q_coolingcoil_wet

            T_in_IHX                    = T_in_CC + Q_coolingcoil / (Cycle.Pump.m_dot_g * cp_g)
            Q_IHX                       = epsilon * Cycle.Pump.m_dot_g * cp_g * (T_in_IHX - T_evap)

            AS.update(CP.PT_INPUTS, p_cond, T_cond - Cycle.DT_sc_target)
            h_target                    = AS.hmass()                # [J/kg]
            Q_cond_enthalpy             = Cycle.Compressor.m_dot_r * (Cycle.Compressor.h_out_r - h_target)
            resids                      = [Q_IHX + W + Q_cond, Q_cond + Q_cond_enthalpy, Q_coolingcoil - Q_IHX]
            return resids

        elif Cycle.Mode == 'HP':
            # Evaporator heat transfer rate
            Q_evap      = epsilon * Cycle.Evaporator.Fins.Air.V_dot_ha * rho_air * (Cycle.Evaporator.Fins.Air.T_db - T_evap) * Cp_air

            # Compressor power
            AS.update(CP.QT_INPUTS, 1.0, T_evap)
            p_evap                      = AS.p()                    # [pa]
            AS.update(CP.QT_INPUTS, 1.0, T_cond)
            p_cond                      = AS.p()                    # [pa]
            Cycle.Compressor.p_in_r     = p_evap
            Cycle.Compressor.p_out_r    = p_cond
            Cycle.Compressor.T_in_r     = T_evap + Cycle.Evaporator.DT_sh
            Cycle.Compressor.Ref        = Cycle.Ref
            Cycle.Compressor.AS         = Cycle.AS
            Cycle.Compressor.Calculate()
            W                           = Cycle.Compressor.W

            # Evaporator will be dry
            Q_coolingcoil   = epsilon * Cycle.CoolingCoil.Fins.Air.V_dot_ha * rho_air * (T_in_CC - Cycle.CoolingCoil.Fins.Air.T_db) * Cp_air

            # Glycol specifi heat
            AS_SLF.update(CP.PT_INPUTS, Cycle.Pump.p_in_g, T_in_CC)
            cp_g                        = AS_SLF.cpmass()           # [J/kg/K]

            T_in_IHX                    = T_in_CC - Q_coolingcoil / (Cycle.Pump.m_dot_g * cp_g)
            Q_IHX                       = epsilon * Cycle.Pump.m_dot_g * cp_g * (T_in_IHX - T_cond)

            AS.update(CP.PT_INPUTS, p_cond, T_cond - Cycle.DT_sc_target)
            h_target                    = AS.hmass()                # [J/kg]
            Q_IHX_enthalpy              = Cycle.Compressor.m_dot_r * (Cycle.Compressor.h_out_r - h_target)

            resids                      = [Q_IHX + W + Q_evap, Q_IHX + Q_IHX_enthalpy, Q_coolingcoil + Q_IHX]
            return resids

    solverFunc                          = fsolve
    if Cycle.Mode == 'AC':
        T_evap_init                     = Cycle.CoolingCoil.Fins.Air.T_db - 15
        T_cond_init                     = Cycle.Condenser.Fins.Air.T_db + 8
        T_in_CC                         = T_evap_init + 1
        # First try using the fsolve algorithm
        try:
            x                           = fsolve(OBJECTIVE, [T_evap_init, T_cond_init, T_in_CC])
        except:
            # If that doesnt work, try the Mult-Dimensional Newton-raphson solver
            try:
                x                       = MultiDimNewtRaph(OBJECTIVE, [T_evap_init, T_cond_init, T_in_CC])
            except:
                x                       = [T_evap_init, T_cond_init, 284]
        DT_evap                         = Cycle.CoolingCoil.Fins.Air.T_db - x[0]
        DT_cond                         = x[1] - Cycle.Condenser.Fins.Air.T_db
        T_in_CC                         = x[2]
    elif Cycle.Mode == 'HP':
        T_evap_init                     = Cycle.Evaporator.Fins.Air.T_db - 8
        T_cond_init                     = Cycle.CoolingCoil.Fins.Air.T_db + 15
        T_in_CC                         = T_cond_init - 1
        x                               = solverFunc(OBJECTIVE, [T_evap_init, T_cond_init, T_in_CC])
        DT_evap                         = Cycle.Evaporator.Fins.Air.T_db - x[0]
        T_in_CC                         = x[2]
        DT_cond                         = x[1] - T_in_CC
    else:
        raise ValueError()

    return DT_evap - 2, DT_cond + 2, T_in_CC
