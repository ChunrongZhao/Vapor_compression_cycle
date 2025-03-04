"""the main components of the vapor compression cycle are
compressor, evaporator, expansion valve and condenser"""

# VCS_components.py
#
# Created: Oct 2023, C.R. Zhao

# ---------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import numpy as np
import CoolProp.CoolProp as CoolProp
from matplotlib import pyplot as plt
from scipy.optimize import fsolve
import pandas as pd
from openpyxl.workbook import Workbook
import pickle
from pathlib import Path
from more_itertools import chunked
from VCS_Codes.Correlations.HTC_correlations import f_h_1phase_Channel
import CoolProp as CP
from ACHP_codes.Correlations.FinStructure_Correlations import HerringboneFins_condenser
from collections import defaultdict
from ACHP_codes.Correlations.HTC_Correlations import f_h_1phase_Tube, ShahCondensation_Average, f_h_1phase_Channel
from ACHP_codes.Correlations.HTC_Correlations import ShahEvaporation_Average


# ----------------------------------------------------------------------
class AttributeDict(defaultdict):
    def __init__(self):
        super(AttributeDict, self).__init__(AttributeDict)

    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError(key)

    def __setattr__(self, key, value):
        self[key] = value


# ---------------------------------------------------------------------
#   Condenser
# ---------------------------------------------------------------------
# finned tube
class VCS_condenser_Finned_Tube_3Zones:
    """
    Class of condenser.
    Inputs:
        refName:    Refrigerant fluid
        T_air:      K,      Outdoor (ambient) air temperature
        m_dot:      kg/s,   refrigerant flowrate
        m_dot_air:  kg/s,   air side flowrate
        P_in:       Pq,     fluid pressure in the inlet
        T_in:       K,      fluid temperature in the inlet
        subcooling  K,      Degree of subcooling exiting the condenser
    """

    def __init__(self, refrigerant, HTF, T_air_inlet, m_dot_ref, m_dot_air, P_i, T_i, SC_degree, HTC_air):
        # working fluid types
        self.refrigerant_name               = refrigerant
        self.heat_transfer_fluid            = HTF           # air
        # air side parameters
        self.T_air_inlet                    = T_air_inlet         # K, inlet air temperature              273
        self.m_dot_air                      = m_dot_air     # kg/s, air mass flow rate              0.4
        self.HTC_air                        = HTC_air       # W/m2K, heat transfer coefficient: refrigerant liquid phase
        # refrigerant side parameters
        self.subcooling_degree              = SC_degree    # -8.4 K, level of refrigerant subcooling 8.4
        self.m_dot_refrigerant              = m_dot_ref     # kg/s, refrigerant mass flow rate      0.02
        # inlet refrigerant states
        self.P_inlet                        = P_i           # Pa, inlet pressure
        self.T_inlet                        = T_i           # K, inlet temperature
        self.h_inlet                        = CoolProp.PropsSI('H', 'P', self.P_inlet, 'T', T_i, self.refrigerant_name)
        self.vapor_quality_inlet            = 1.0
        self.rho_inlet                      = CoolProp.PropsSI('D', 'P', self.P_inlet, 'T', T_i, self.refrigerant_name)
        # outlet refrigerant states
        self.P_outlet                       = None
        self.T_outlet                       = None
        self.h_outlet                       = None
        self.vapor_quality_outlet           = None
        self.Q_condenser                    = None
        # 3 zone lengths
        self.L_1                            = None
        self.L_2                            = None
        self.L_3                            = None
        # outlet HTF(air) states
        self.T_air_out                      = None
        self.UA_tot                         = None
        self.effectiveness                  = None

    def thermodynamics_calculation(self):
        # Geometric parameters, finned-tube structure
        L_tube                              = 10            # m, tube length
        Number_of_tubes                     = 5             # number of tubes in parallel
        ID_tube                             = 0.01          # m, tube inner diameter    0.01
        fin_efficiency                      = 0.8           # fin efficiency
        fin_ratio                           = 7.2           # fin ratio

        HTC_liq                             = 600           # W/m2K, heat transfer coefficient: refrigerant liquid phase
        HTC_2ph                             = 3000          # W/m2K, heat transfer coefficient: refrigerant two phases
        HTC_vap                             = 600           # W/m2K, heat transfer coefficient: refrigerant vapor phase

        #   ------------------------------------------------------------------------------
        # mass flow rate of refrigerant in single tube
        mfr_channel                         = self.m_dot_refrigerant / Number_of_tubes

        # obtaining overall heat transfer coefficient, U, for each zone
        U_vap                               = 1 / (1 / HTC_vap + 1 / (self.HTC_air * fin_efficiency * fin_ratio))
        U_2ph                               = 1 / (1 / HTC_2ph + 1 / (self.HTC_air * fin_efficiency * fin_ratio))
        U_liq                               = 1 / (1 / HTC_liq + 1 / (self.HTC_air * fin_efficiency * fin_ratio))

        # Pressure drop within the condenser estimation
        f_refrigerant                       = 0.0015   # fractional loss
        u_refrigerant                       = 4 * mfr_channel / (np.pi * ID_tube**2 * self.rho_inlet)    # refrigerant flow velocity
        delta_p_refrigerant                 = f_refrigerant * (ID_tube / L_tube) * (u_refrigerant**2 / 2)  # pressure drop in the tube
        # outlet pressure is set equal to inlet pressure with a refrigerant pump
        P_outlet                            = self.P_inlet

        # 3 zones energy balance
        U_1                                 = U_vap
        U_2                                 = U_2ph
        U_3                                 = U_liq

        # estimate the specific heat based on the inlet pressure and enthalpy
        cp_1                                = CoolProp.PropsSI('C', 'P', self.P_inlet, 'H', self.h_inlet, self.refrigerant_name)
        # estimate enthalpy of saturated vapor, x = 1
        h_sat_vap                           = CoolProp.PropsSI('H', 'P', self.P_inlet, 'Q', 1, self.refrigerant_name)
        T_sat_vap                           = CoolProp.PropsSI('T', 'P', self.P_inlet, 'Q', 1, self.refrigerant_name)
        # estimate enthalpy of saturated liquid, x = 0
        h_sat_liq                           = CoolProp.PropsSI('H', 'P', P_outlet, 'Q', 0, self.refrigerant_name)
        T_sat_liq                           = CoolProp.PropsSI('T', 'P', P_outlet, 'Q', 0, self.refrigerant_name)
        # estimate the specific heat based on the outlet pressure and enthalpy
        T_outlet_design                     = T_sat_liq - self.subcooling_degree
        cp_3                                = CoolProp.PropsSI('C', 'P', P_outlet, 'T', T_outlet_design, self.refrigerant_name)

        """algorithm to get 3-zone info;  
        5 thermophysical outputs are: [pressure, temperature, enthalpy, vapor quality, transferred heat]; and
        3 zone lengths are: [L_1, L_2, L_3] """

        # if T_sat_vap < T_air, the superheated vapor will be cooled down but will not condense
        if T_sat_vap - self.T_air_inlet < 0:
            self.L_1                        = L_tube
            self.L_2                        = 0.
            self.L_3                        = 0.
            # output parameters; scenario - 1
            self.P_outlet                   = P_outlet
            self.T_outlet                   = self.getT(self.L_1, U_1, self.T_inlet, cp_1, self.T_air_inlet, mfr_channel, ID_tube)
            self.h_outlet                   = CoolProp.PropsSI('H', 'P', P_outlet, 'T', self.T_outlet, self.refrigerant_name)
            self.vapor_quality_outlet       = 1.0
            self.Q_condenser                = self.m_dot_refrigerant * (self.h_outlet - self.h_inlet)

        else:
            # if T_sat_vap > T_air, then the vapor will condense;
            # first to determine the characteristic temperature at the superheated (SH) zone
            T_SH                            = (self.T_inlet + T_sat_vap) / 2
            delta_h_SH                      = abs(self.h_inlet - h_sat_vap)
            self.L_1                        = self.getL(T_SH, delta_h_SH, U_1, self.T_air_inlet, mfr_channel, ID_tube)
            self.L_2                        = L_tube - self.L_1

            if self.L_2 < 0:
                # only 1 zone exists;
                self.L_1                    = L_tube
                self.L_2                    = 0.
                self.L_3                    = 0.
                # output parameters; scenario - 2
                self.P_outlet               = P_outlet
                self.T_outlet               = self.getT(self.L_1, U_1, self.T_inlet, cp_1, self.T_air_inlet, mfr_channel, ID_tube)
                self.h_outlet               = CoolProp.PropsSI('H', 'P', P_outlet, 'T', self.T_outlet, self.refrigerant_name)
                self.vapor_quality_outlet   = 1.0
                self.Q_condenser            = self.m_dot_refrigerant * (self.h_outlet - self.h_inlet)

            else:
                # check if 2 zones exist; if P_inlet = P_outlet, then T_sat_vap = T_sat_liq
                T_2phase                    = (T_sat_vap + T_sat_liq) / 2
                # enthalpy change in the two-phase zone
                delta_h_2phase              = abs(h_sat_vap - h_sat_liq)
                self.L_2                    = self.getL(T_2phase, delta_h_2phase, U_2, self.T_air_inlet, mfr_channel, ID_tube)
                self.L_3                    = L_tube - self.L_1 - self.L_2

                if self.L_3 < 0:
                    # only 2 zones exist
                    self.L_2                = L_tube - self.L_1
                    self.L_3                = 0
                    h_2                     = self.getHout(T_2phase, self.L_2, U_2, self.T_air_inlet, mfr_channel, ID_tube, h_sat_vap)
                    T_2                     = CoolProp.PropsSI('T', 'P', P_outlet, 'H', h_2, self.refrigerant_name)
                    # output parameters; scenario - 3
                    self.T_outlet           = T_2
                    self.P_outlet           = P_outlet
                    self.h_outlet           = h_2
                    self.vapor_quality_outlet = CoolProp.PropsSI('Q', 'P', P_outlet, 'H', self.h_outlet, self.refrigerant_name)
                    self.Q_condenser        = self.m_dot_refrigerant * (h_2 - self.h_inlet)

                else:
                    T_3                     = self.getT(self.L_3, U_3, T_sat_liq, cp_3, self.T_air_inlet, mfr_channel, ID_tube)
                    h_3                     = CoolProp.PropsSI('H', 'P', P_outlet, 'T', T_3, self.refrigerant_name)
                    # output parameters; scenario - 4
                    self.T_outlet           = T_3
                    self.P_outlet           = P_outlet
                    self.h_outlet           = h_3
                    self.vapor_quality_outlet = 0.
                    self.Q_condenser        = self.m_dot_refrigerant * (h_3 - self.h_inlet)

            # calculate air outlet temperature
            cp_air                          = CoolProp.PropsSI('C', 'P', 101325, 'T', self.T_air_inlet, self.heat_transfer_fluid)
            C_air                           = cp_air * self.m_dot_air
            self.T_air_out                  = self.T_air_inlet - self.Q_condenser / C_air
            # condenser overall conductance and effectiveness
            self.UA_tot                     = self.L_3 / L_tube * U_3 + self.L_2 / L_tube * U_2 + self.L_1 / L_tube * U_1
            self.effectiveness              = self.Q_condenser / (C_air * (self.T_air_inlet - self.T_inlet))

        condenser_results          = results_Data_VCS_condenser_Finned_Tube_3Zones()
        condenser_results.T_outlet                       = self.T_outlet
        condenser_results.P_outlet                       = self.P_outlet
        condenser_results.h_outlet                       = self.h_outlet
        condenser_results.x_outlet                       = self.vapor_quality_outlet
        condenser_results.Q_condenser                    = self.Q_condenser
        condenser_results.T_air_out                      = self.T_air_out
        condenser_results.UA_tot                         = self.UA_tot
        condenser_results.effectiveness                  = self.effectiveness

        condenser_results.L_tot                          = L_tube
        condenser_results.L_1                            = self.L_1
        condenser_results.L_2                            = self.L_2
        condenser_results.L_3                            = self.L_3

        return condenser_results

    def getL(self, T_char, dh, U, T_air, mfr_channel, ID_tube):
        dT                          = abs(T_char - T_air)
        L                           = (mfr_channel * dh) / (dT * U * np.pi * ID_tube)
        return L

    def getHout(self, Tchar, length, U, T_air, mfr_channel, ID, h_in):
        UA                          = U * (np.pi * ID * length)
        dT                          = (Tchar - T_air)
        h                           = h_in - UA * dT / mfr_channel
        return h

    def getT(self, L_zone, U, T_in, cp, T_air, mfr_channel, ID_tube):
        # get outlet temperature of a zone
        UA                          = U * (np.pi * ID_tube * L_zone)
        T_outlet                    = (UA * (T_air - 0.5 * T_in) + mfr_channel * cp * T_in) / (mfr_channel * cp + 0.5 * UA)
        # cannot be lower than the air temperature
        if (T_outlet - T_air) < 0:
            T_outlet                = T_air

        return T_outlet


# with herringbone fin details
class VCS_condenser_Finned_Tube_3Zones_herringbone:
    """
    Class of condenser.
    Inputs:
        refName:    Refrigerant fluid
        T_air:      K,      Outdoor (ambient) air temperature
        m_dot:      kg/s,   refrigerant flowrate
        m_dot_air:  kg/s,   air side flowrate
        P_in:       Pq,     fluid pressure in the inlet
        T_in:       K,      fluid temperature in the inlet
        subcooling  K,      Degree of subcooling exiting the condenser
    """

    def __init__(self, refrigerant, HTF, T_air_inlet, m_dot_ref, m_dot_air, P_i, T_i, SC_degree, HTC_air, eta_fin, fin_ratio):
        # working fluid types
        self.refrigerant_name               = refrigerant
        self.heat_transfer_fluid            = HTF           # air
        # air side parameters
        self.T_air_inlet                    = T_air_inlet         # K, inlet air temperature              273
        self.m_dot_air                      = m_dot_air     # kg/s, air mass flow rate              0.4
        self.HTC_air                        = HTC_air       # W/m2K, heat transfer coefficient: refrigerant liquid phase
        self.eta_fin                        = eta_fin
        self.fin_ratio                      = fin_ratio
        # refrigerant side parameters
        self.subcooling_degree              = SC_degree    # -8.4 K, level of refrigerant subcooling 8.4
        self.m_dot_refrigerant              = m_dot_ref     # kg/s, refrigerant mass flow rate      0.02
        # inlet refrigerant states
        self.P_inlet                        = P_i           # Pa, inlet pressure
        self.T_inlet                        = T_i           # K, inlet temperature
        self.h_inlet                        = CoolProp.PropsSI('H', 'P', self.P_inlet, 'T', T_i, self.refrigerant_name)
        self.vapor_quality_inlet            = 1.0
        self.rho_inlet                      = CoolProp.PropsSI('D', 'P', self.P_inlet, 'T', T_i, self.refrigerant_name)
        # outlet refrigerant states
        self.P_outlet                       = None
        self.T_outlet                       = None
        self.h_outlet                       = None
        self.vapor_quality_outlet           = None
        self.Q_condenser                    = None
        # 3 zone lengths
        self.L_1                            = None
        self.L_2                            = None
        self.L_3                            = None
        # outlet HTF(air) states
        self.T_air_out                      = None
        self.UA_tot                         = None
        self.effectiveness                  = None

    def thermodynamics_calculation(self):
        # Geometric parameters, finned-tube structure
        L_tube                              = 1            # m, tube length
        Number_of_tubes                     = 50             # number of tubes in parallel
        ID_tube                             = 0.01          # m, tube inner diameter    0.01
        fin_efficiency                      = self.eta_fin           # fin efficiency
        fin_ratio                           = self.fin_ratio           # fin ratio

        HTC_liq                             = 1414          # W/m2K, heat transfer coefficient: refrigerant liquid phase
        HTC_2ph                             = 4274        # W/m2K, heat transfer coefficient: refrigerant two phases
        HTC_vap                             = 1291        # W/m2K, heat transfer coefficient: refrigerant vapor phase

        # --------------------------------------------------------------------------------
        rho_chan                            = 2719
        k_chan                              = 202.4

        #   ------------------------------------------------------------------------------
        # mass flow rate of refrigerant in single tube
        mfr_channel                         = self.m_dot_refrigerant / Number_of_tubes

        # obtaining overall heat transfer coefficient, U, for each zone
        U_vap                               = 1 / (1 / HTC_vap + 1 / (self.HTC_air * fin_efficiency * fin_ratio))
        U_2ph                               = 1 / (1 / HTC_2ph + 1 / (self.HTC_air * fin_efficiency * fin_ratio))
        U_liq                               = 1 / (1 / HTC_liq + 1 / (self.HTC_air * fin_efficiency * fin_ratio))

        # Pressure drop within the condenser estimation
        f_refrigerant                       = 0.0015   # fractional loss
        u_refrigerant                       = 4 * mfr_channel / (np.pi * ID_tube**2 * self.rho_inlet)    # refrigerant flow velocity
        delta_p_refrigerant                 = f_refrigerant * (ID_tube / L_tube) * (u_refrigerant**2 / 2)  # pressure drop in the tube
        # outlet pressure is set equal to inlet pressure with a refrigerant pump
        P_outlet                            = self.P_inlet

        # 3 zones energy balance
        U_1                                 = U_vap
        U_2                                 = U_2ph
        U_3                                 = U_liq

        # estimate the specific heat based on the inlet pressure and enthalpy
        cp_1                                = CoolProp.PropsSI('C', 'P', self.P_inlet, 'H', self.h_inlet, self.refrigerant_name)
        # estimate enthalpy of saturated vapor, x = 1
        h_sat_vap                           = CoolProp.PropsSI('H', 'P', self.P_inlet, 'Q', 1, self.refrigerant_name)
        T_sat_vap                           = CoolProp.PropsSI('T', 'P', self.P_inlet, 'Q', 1, self.refrigerant_name)
        # estimate enthalpy of saturated liquid, x = 0
        h_sat_liq                           = CoolProp.PropsSI('H', 'P', P_outlet, 'Q', 0, self.refrigerant_name)
        T_sat_liq                           = CoolProp.PropsSI('T', 'P', P_outlet, 'Q', 0, self.refrigerant_name)
        # estimate the specific heat based on the outlet pressure and enthalpy
        T_outlet_design                     = T_sat_liq - self.subcooling_degree
        cp_3                                = CoolProp.PropsSI('C', 'P', P_outlet, 'T', T_outlet_design, self.refrigerant_name)

        """algorithm to get 3-zone info;  
        5 thermophysical outputs are: [pressure, temperature, enthalpy, vapor quality, transferred heat]; and
        3 zone lengths are: [L_1, L_2, L_3] """

        # if T_sat_vap < T_air, the superheated vapor will be cooled down but will not condense
        if T_sat_vap - self.T_air_inlet < 0:
            self.L_1                        = L_tube
            self.L_2                        = 0.
            self.L_3                        = 0.
            # output parameters; scenario - 1
            self.P_outlet                   = P_outlet
            self.T_outlet                   = self.getT(self.L_1, U_1, self.T_inlet, cp_1, self.T_air_inlet, mfr_channel, ID_tube)
            self.h_outlet                   = CoolProp.PropsSI('H', 'P', P_outlet, 'T', self.T_outlet, self.refrigerant_name)
            self.vapor_quality_outlet       = 1.0
            self.Q_condenser                = self.m_dot_refrigerant * (self.h_outlet - self.h_inlet)

        else:
            # if T_sat_vap > T_air, then the vapor will condense;
            # first to determine the characteristic temperature at the superheated (SH) zone
            T_SH                            = (self.T_inlet + T_sat_vap) / 2
            delta_h_SH                      = abs(self.h_inlet - h_sat_vap)
            self.L_1                        = self.getL(T_SH, delta_h_SH, U_1, self.T_air_inlet, mfr_channel, ID_tube)
            self.L_2                        = L_tube - self.L_1

            if self.L_2 < 0:
                # only 1 zone exists;
                self.L_1                    = L_tube
                self.L_2                    = 0.
                self.L_3                    = 0.
                # output parameters; scenario - 2
                self.P_outlet               = P_outlet
                self.T_outlet               = self.getT(self.L_1, U_1, self.T_inlet, cp_1, self.T_air_inlet, mfr_channel, ID_tube)
                self.h_outlet               = CoolProp.PropsSI('H', 'P', P_outlet, 'T', self.T_outlet, self.refrigerant_name)
                self.vapor_quality_outlet   = 1.0
                self.Q_condenser            = self.m_dot_refrigerant * (self.h_outlet - self.h_inlet)

            else:
                # check if 2 zones exist; if P_inlet = P_outlet, then T_sat_vap = T_sat_liq
                T_2phase                    = (T_sat_vap + T_sat_liq) / 2
                # enthalpy change in the two-phase zone
                delta_h_2phase              = abs(h_sat_vap - h_sat_liq)
                self.L_2                    = self.getL(T_2phase, delta_h_2phase, U_2, self.T_air_inlet, mfr_channel, ID_tube)
                self.L_3                    = L_tube - self.L_1 - self.L_2

                if self.L_3 < 0:
                    # only 2 zones exist
                    self.L_2                = L_tube - self.L_1
                    self.L_3                = 0
                    h_2                     = self.getHout(T_2phase, self.L_2, U_2, self.T_air_inlet, mfr_channel, ID_tube, h_sat_vap)
                    T_2                     = CoolProp.PropsSI('T', 'P', P_outlet, 'H', h_2, self.refrigerant_name)
                    # output parameters; scenario - 3
                    self.T_outlet           = T_2
                    self.P_outlet           = P_outlet
                    self.h_outlet           = h_2
                    self.vapor_quality_outlet = CoolProp.PropsSI('Q', 'P', P_outlet, 'H', self.h_outlet, self.refrigerant_name)
                    self.Q_condenser        = self.m_dot_refrigerant * (h_2 - self.h_inlet)

                else:
                    T_3                     = self.getT(self.L_3, U_3, T_sat_liq, cp_3, self.T_air_inlet, mfr_channel, ID_tube)
                    h_3                     = CoolProp.PropsSI('H', 'P', P_outlet, 'T', T_3, self.refrigerant_name)
                    # output parameters; scenario - 4
                    self.T_outlet           = T_3
                    self.P_outlet           = P_outlet
                    self.h_outlet           = h_3
                    self.vapor_quality_outlet = 0.
                    self.Q_condenser        = self.m_dot_refrigerant * (h_3 - self.h_inlet)

            # calculate air outlet temperature
            cp_air                          = CoolProp.PropsSI('C', 'P', 101325, 'T', self.T_air_inlet, self.heat_transfer_fluid)
            C_air                           = cp_air * self.m_dot_air
            self.T_air_out                  = self.T_air_inlet - self.Q_condenser / C_air
            # condenser overall conductance and effectiveness
            self.UA_tot                     = self.L_3 / L_tube * U_3 + self.L_2 / L_tube * U_2 + self.L_1 / L_tube * U_1
            self.effectiveness              = self.Q_condenser / (C_air * (self.T_air_inlet - self.T_inlet))

        condenser_results          = results_Data_VCS_condenser_Finned_Tube_3Zones()
        condenser_results.T_outlet                       = self.T_outlet
        condenser_results.P_outlet                       = self.P_outlet
        condenser_results.h_outlet                       = self.h_outlet
        condenser_results.x_outlet                       = self.vapor_quality_outlet
        condenser_results.Q_condenser                    = self.Q_condenser
        condenser_results.T_air_out                      = self.T_air_out
        condenser_results.UA_tot                         = self.UA_tot
        condenser_results.effectiveness                  = self.effectiveness

        condenser_results.L_tot                          = L_tube
        condenser_results.L_1                            = self.L_1
        condenser_results.L_2                            = self.L_2
        condenser_results.L_3                            = self.L_3

        return condenser_results

    def getL(self, T_char, dh, U, T_air, mfr_channel, ID_tube):
        dT                          = abs(T_char - T_air)
        L                           = (mfr_channel * dh) / (dT * U * np.pi * ID_tube)
        return L

    def getHout(self, Tchar, length, U, T_air, mfr_channel, ID, h_in):
        UA                          = U * (np.pi * ID * length)
        dT                          = (Tchar - T_air)
        h                           = h_in - UA * dT / mfr_channel
        return h

    def getT(self, L_zone, U, T_in, cp, T_air, mfr_channel, ID_tube):
        # get outlet temperature of a zone
        UA                          = U * (np.pi * ID_tube * L_zone)
        T_outlet                    = (UA * (T_air - 0.5 * T_in) + mfr_channel * cp * T_in) / (mfr_channel * cp + 0.5 * UA)
        # cannot be lower than the air temperature
        if (T_outlet - T_air) < 0:
            T_outlet                = T_air

        return T_outlet


# with herringbone fin details
class VCS_condenser_Finned_Tube_3Zones_herringbone001:
    """
    Class of condenser.
    Inputs:
        refName:    Refrigerant fluid
        T_air:      K,      Outdoor (ambient) air temperature
        m_dot:      kg/s,   refrigerant flowrate
        m_dot_air:  kg/s,   air side flowrate
        P_in:       Pq,     fluid pressure in the inlet
        T_in:       K,      fluid temperature in the inlet
        subcooling  K,      Degree of subcooling exiting the condenser
    """

    def __init__(self, refrigerant, HTF, m_dot_ref, P_i, T_i, SC_degree, T_air_inlet, V_dot_air):
        # working fluid types
        self.refrigerant_name               = refrigerant
        self.heat_transfer_fluid            = HTF           # air
        # air side parameters
        self.T_air_inlet                    = T_air_inlet   # K, inlet air temperature              273
        self.V_dot_ha                       = V_dot_air     # m3/s, rated volumetric flowrate
        # refrigerant side parameters
        self.subcooling_degree              = SC_degree    # -8.4 K, level of refrigerant subcooling 8.4
        self.m_dot_refrigerant              = m_dot_ref     # kg/s, refrigerant mass flow rate      0.02
        # inlet refrigerant states
        self.P_inlet                        = P_i           # Pa, inlet pressure
        self.T_inlet                        = T_i           # K, inlet temperature
        self.h_inlet                        = CoolProp.PropsSI('H', 'P', self.P_inlet, 'T', T_i, self.refrigerant_name)
        self.vapor_quality_inlet            = 1.0
        self.rho_inlet                      = CoolProp.PropsSI('D', 'P', self.P_inlet, 'T', T_i, self.refrigerant_name)
        # outlet refrigerant states
        self.P_outlet                       = None
        self.T_outlet                       = None
        self.h_outlet                       = None
        self.vapor_quality_outlet           = None
        self.Q_condenser                    = None
        # 3 zone lengths
        self.L_1                            = None
        self.L_2                            = None
        self.L_3                            = None
        # outlet HTF(air) states
        self.T_air_out                      = None
        self.UA_tot                         = None
        self.effectiveness                  = None

    def HerringboneFin(self):
        Fins        = AttributeDict()
        Fins.Tubes.N_Tubes_per_bank     = 5                    # number of tubes per bank=row
        Fins.Tubes.N_bank               = 1                     # number of banks/rows
        Fins.Tubes.L_tube               = 10
        Fins.Tubes.OD                   = 0.011
        Fins.Tubes.ID                   = 0.01
        Fins.Tubes.Pl                   = 0.0191                # distance between center of tubes in flow direction
        Fins.Tubes.Pt                   = 0.0254                # distance between center of tubes orthogonal to flow direction
        Fins.Tubes.k_w                  = 237                   # wall thermal conductivity (i.e. pipe material)

        Fins.Fins.FPI                   = 25                    # Number of fins per inch
        Fins.Fins.Pd                    = 0.001                 # 2* amplitude of wavy fin
        Fins.Fins.xf                    = 0.001                 # 1/2 period of fin
        Fins.Fins.t                     = 0.00011               # Thickness of fin material
        Fins.Fins.k_fin                 = 237                   # Thermal conductivity of fin material

        Fins.Air.V_dot_ha               = self.V_dot_ha         # rated volumetric flowrate
        Fins.Air.T_mean                 = self.T_air_inlet
        Fins.Air.T_db                   = self.T_air_inlet
        Fins.Air.p                      = 101325                # Condenser Air pressures in Pa
        Fins.Air.RH                     = 0.51
        Fins.Air.RH_mean                = 0.51
        Fins.Air.FanPower               = 260

        A, eta_o, h_a, DeltaP_air, m_dot_da, m_dot_ha = HerringboneFins_condenser(Fins)

        self.m_dot_air_cal              = m_dot_ha
        self.A_fin                      = A
        self.eta_fin                    = eta_o
        self.HTC_air_cal                = h_a

        # print(m_dot_ha, m_dot_da)

    def thermodynamics_calculation(self):
        # Geometric parameters, finned-tube structure
        L_tube                              = 10            # m, tube length
        Number_of_tubes                     = 5             # number of tubes in parallel
        ID_tube                             = 0.01          # m, tube inner diameter    0.01
        # overall fin efficiency
        self.HerringboneFin()
        fin_efficiency_cal                  = self.eta_fin
        HTC_air_cal                         = self.HTC_air_cal
        m_dot_air_cal                       = self.m_dot_air_cal

        # ----------------------------------------------
        fin_efficiency                      = 0.8 * 0.72
        m_dot_air                           = 0.47

        HTC_liq                             = 600           # W/m2K, heat transfer coefficient: refrigerant liquid phase
        HTC_2ph                             = 3000          # W/m2K, heat transfer coefficient: refrigerant two phases
        HTC_vap                             = 600           # W/m2K, heat transfer coefficient: refrigerant vapor phase

        HTC_air                             = 100           # W/m2K, heat transfer coefficient: air
        #   ------------------------------------------------------------------------------
        # mass flow rate of refrigerant in single tube
        mfr_channel                         = self.m_dot_refrigerant / Number_of_tubes

        # obtaining overall heat transfer coefficient, U, for each zone
        U_vap                               = 1 / (1 / HTC_vap + 1 / (HTC_air * fin_efficiency))
        U_2ph                               = 1 / (1 / HTC_2ph + 1 / (HTC_air * fin_efficiency))
        U_liq                               = 1 / (1 / HTC_liq + 1 / (HTC_air * fin_efficiency))

        # Pressure drop within the condenser estimation
        f_refrigerant                       = 0.0015   # fractional loss
        u_refrigerant                       = 4 * mfr_channel / (np.pi * ID_tube**2 * self.rho_inlet)    # refrigerant flow velocity
        delta_p_refrigerant                 = f_refrigerant * (ID_tube / L_tube) * (u_refrigerant**2 / 2)  # pressure drop in the tube
        # outlet pressure is set equal to inlet pressure with a refrigerant pump
        P_outlet                            = self.P_inlet

        # 3 zones energy balance
        U_1                                 = U_vap
        U_2                                 = U_2ph
        U_3                                 = U_liq

        # estimate the specific heat based on the inlet pressure and enthalpy
        cp_1                                = CoolProp.PropsSI('C', 'P', self.P_inlet, 'H', self.h_inlet, self.refrigerant_name)
        # estimate enthalpy of saturated vapor, x = 1
        h_sat_vap                           = CoolProp.PropsSI('H', 'P', self.P_inlet, 'Q', 1, self.refrigerant_name)
        T_sat_vap                           = CoolProp.PropsSI('T', 'P', self.P_inlet, 'Q', 1, self.refrigerant_name)
        # estimate enthalpy of saturated liquid, x = 0
        h_sat_liq                           = CoolProp.PropsSI('H', 'P', P_outlet, 'Q', 0, self.refrigerant_name)
        T_sat_liq                           = CoolProp.PropsSI('T', 'P', P_outlet, 'Q', 0, self.refrigerant_name)
        # estimate the specific heat based on the outlet pressure and enthalpy
        T_outlet_design                     = T_sat_liq - self.subcooling_degree
        cp_3                                = CoolProp.PropsSI('C', 'P', P_outlet, 'T', T_outlet_design, self.refrigerant_name)

        """algorithm to get 3-zone info;  
        5 thermophysical outputs are: [pressure, temperature, enthalpy, vapor quality, transferred heat]; and
        3 zone lengths are: [L_1, L_2, L_3] """

        # if T_sat_vap < T_air, the superheated vapor will be cooled down but will not condense
        if T_sat_vap - self.T_air_inlet < 0:
            self.L_1                        = L_tube
            self.L_2                        = 0.
            self.L_3                        = 0.
            # output parameters; scenario - 1
            self.P_outlet                   = P_outlet
            self.T_outlet                   = self.getT(self.L_1, U_1, self.T_inlet, cp_1, self.T_air_inlet, mfr_channel, ID_tube)
            self.h_outlet                   = CoolProp.PropsSI('H', 'P', P_outlet, 'T', self.T_outlet, self.refrigerant_name)
            self.vapor_quality_outlet       = 1.0
            self.Q_condenser                = self.m_dot_refrigerant * (self.h_outlet - self.h_inlet)

        else:
            # if T_sat_vap > T_air, then the vapor will condense;
            # first to determine the characteristic temperature at the superheated (SH) zone
            T_SH                            = (self.T_inlet + T_sat_vap) / 2
            delta_h_SH                      = abs(self.h_inlet - h_sat_vap)
            self.L_1                        = self.getL(T_SH, delta_h_SH, U_1, self.T_air_inlet, mfr_channel, ID_tube)
            self.L_2                        = L_tube - self.L_1

            if self.L_2 < 0:
                # only 1 zone exists;
                self.L_1                    = L_tube
                self.L_2                    = 0.
                self.L_3                    = 0.
                # output parameters; scenario - 2
                self.P_outlet               = P_outlet
                self.T_outlet               = self.getT(self.L_1, U_1, self.T_inlet, cp_1, self.T_air_inlet, mfr_channel, ID_tube)
                self.h_outlet               = CoolProp.PropsSI('H', 'P', P_outlet, 'T', self.T_outlet, self.refrigerant_name)
                self.vapor_quality_outlet   = 1.0
                self.Q_condenser            = self.m_dot_refrigerant * (self.h_outlet - self.h_inlet)

            else:
                # check if 2 zones exist; if P_inlet = P_outlet, then T_sat_vap = T_sat_liq
                T_2phase                    = (T_sat_vap + T_sat_liq) / 2
                # enthalpy change in the two-phase zone
                delta_h_2phase              = abs(h_sat_vap - h_sat_liq)
                self.L_2                    = self.getL(T_2phase, delta_h_2phase, U_2, self.T_air_inlet, mfr_channel, ID_tube)
                self.L_3                    = L_tube - self.L_1 - self.L_2

                if self.L_3 < 0:
                    # only 2 zones exist
                    self.L_2                = L_tube - self.L_1
                    self.L_3                = 0
                    h_2                     = self.getHout(T_2phase, self.L_2, U_2, self.T_air_inlet, mfr_channel, ID_tube, h_sat_vap)
                    T_2                     = CoolProp.PropsSI('T', 'P', P_outlet, 'H', h_2, self.refrigerant_name)
                    # output parameters; scenario - 3
                    self.T_outlet           = T_2
                    self.P_outlet           = P_outlet
                    self.h_outlet           = h_2
                    self.vapor_quality_outlet = CoolProp.PropsSI('Q', 'P', P_outlet, 'H', self.h_outlet, self.refrigerant_name)
                    self.Q_condenser        = self.m_dot_refrigerant * (h_2 - self.h_inlet)

                else:
                    T_3                     = self.getT(self.L_3, U_3, T_sat_liq, cp_3, self.T_air_inlet, mfr_channel, ID_tube)
                    h_3                     = CoolProp.PropsSI('H', 'P', P_outlet, 'T', T_3, self.refrigerant_name)
                    # output parameters; scenario - 4
                    self.T_outlet           = T_3
                    self.P_outlet           = P_outlet
                    self.h_outlet           = h_3
                    self.vapor_quality_outlet = 0.
                    self.Q_condenser        = self.m_dot_refrigerant * (h_3 - self.h_inlet)

            # calculate air outlet temperature
            cp_air                          = CoolProp.PropsSI('C', 'P', 101325, 'T', self.T_air_inlet, self.heat_transfer_fluid)
            C_air                           = cp_air * m_dot_air
            self.T_air_out                  = self.T_air_inlet - self.Q_condenser / C_air
            # condenser overall conductance and effectiveness
            self.UA_tot                     = self.L_3 / L_tube * U_3 + self.L_2 / L_tube * U_2 + self.L_1 / L_tube * U_1
            self.effectiveness              = self.Q_condenser / (C_air * (self.T_air_inlet - self.T_inlet))

        condenser_results          = results_Data_VCS_condenser_Finned_Tube_3Zones()
        condenser_results.T_outlet                       = self.T_outlet
        condenser_results.P_outlet                       = self.P_outlet
        condenser_results.h_outlet                       = self.h_outlet
        condenser_results.x_outlet                       = self.vapor_quality_outlet
        condenser_results.Q_condenser                    = self.Q_condenser
        condenser_results.T_air_out                      = self.T_air_out
        condenser_results.UA_tot                         = self.UA_tot
        condenser_results.effectiveness                  = self.effectiveness

        condenser_results.L_tot                          = L_tube
        condenser_results.L_1                            = self.L_1
        condenser_results.L_2                            = self.L_2
        condenser_results.L_3                            = self.L_3

        return condenser_results

    def getL(self, T_char, dh, U, T_air, mfr_channel, ID_tube):
        dT                          = abs(T_char - T_air)
        L                           = (mfr_channel * dh) / (dT * U * np.pi * ID_tube)
        return L

    def getHout(self, Tchar, length, U, T_air, mfr_channel, ID, h_in):
        UA                          = U * (np.pi * ID * length)
        dT                          = (Tchar - T_air)
        h                           = h_in - UA * dT / mfr_channel
        return h

    def getT(self, L_zone, U, T_in, cp, T_air, mfr_channel, ID_tube):
        # get outlet temperature of a zone
        UA                          = U * (np.pi * ID_tube * L_zone)
        T_outlet                    = (UA * (T_air - 0.5 * T_in) + mfr_channel * cp * T_in) / (mfr_channel * cp + 0.5 * UA)
        # cannot be lower than the air temperature
        if (T_outlet - T_air) < 0:
            T_outlet                = T_air

        return T_outlet


# data storage
class results_Data_VCS_condenser_Finned_Tube_3Zones:
    def __init__(self):
        self.T_outlet                       = []
        self.P_outlet                       = []
        self.h_outlet                       = []
        self.x_outlet                       = []
        self.Q_condenser                    = []
        self.T_air_out                      = []
        self.UA_tot                         = []
        self.effectiveness                  = []

        self.L_tot                          = []
        self.L_1                            = []
        self.L_2                            = []
        self.L_3                            = []

    def save_data(self, filename):
        with open(filename, 'wb') as f:
            pickle.dump(self, f)


# ---------------------------------------------------------------------
#   Evaporator
# ---------------------------------------------------------------------

# finned tube
class VCS_evaporator_Finned_Tube_2Zones:
    """
    Class of evaporator
    .
    Inputs:
        refName:        Refrigerant fluid
        T_air:          K,      indoor (room) air temperature
        m_dot:          kg/s,   refrigerant flow rate
        m_dot_liq:      kg/s,   liq coolant side flow rate
        P_in:           Pq,     fluid pressure in the inlet
        T_in:           K,      fluid temperature in the inlet
        superheating    K,      Degree of superheating exiting the condenser
    """

    def __init__(self, refrigerant, T_coolant_inlet, m_dot_ref, m_dot_coolant, P_i, h_i, SH_degree, HTC_coolant, eta_fin_c, fin_ratio_c):
        # working fluid types
        self.refrigerant_name               = refrigerant

        # air side parameters
        self.T_coolant_inlet                = T_coolant_inlet           # K, inlet coolant temperature              273
        self.m_dot_coolant                  = m_dot_coolant             # kg/s, coolant mass flow rate              0.4
        self.HTC_coolant                    = HTC_coolant               # W/m2K, heat transfer coefficient: liquid coolant
        # refrigerant side parameters
        self.superheating_degree            = SH_degree                 # 11.2 K, level of refrigerant superheating
        self.m_dot_refrigerant              = m_dot_ref                 # kg/s, refrigerant mass flow rate      0.02
        # inlet refrigerant states
        self.P_inlet                        = P_i                       # Pa, inlet pressure
        self.h_inlet                        = h_i                       # K,  inlet temperature
        self.T_inlet                        = CoolProp.PropsSI('T', 'P', self.P_inlet, 'H', self.h_inlet, self.refrigerant_name)
        self.vapor_quality_inlet            = CoolProp.PropsSI('Q', 'P', self.P_inlet, 'H', self.h_inlet, self.refrigerant_name)
        self.rho_inlet                      = CoolProp.PropsSI('D', 'P', self.P_inlet, 'H', self.h_inlet, self.refrigerant_name)

        # outlet refrigerant states
        self.P_outlet                       = None
        self.T_outlet                       = None
        self.h_outlet                       = None
        self.vapor_quality_outlet           = None
        self.Q_evaporator                   = None
        # 3 zone lengths
        self.L_1                            = None
        self.L_2                            = None
        self.L_3                            = None
        # outlet HTF(air) states
        self.T_coolant_out                  = None
        self.UA_tot                         = None
        self.effectiveness                  = None

        # ---------------added at 20/11/2024---------------------------------
        self.eta_fin_c                      = eta_fin_c
        self.fin_ratio_c                    = fin_ratio_c

    def thermodynamics_calculation(self):
        # Geometric parameters, finned-tube structure
        L_tube                              = 10            # m, tube length
        Number_of_tubes                     = 5             # number of tubes in parallel
        ID_tube                             = 0.01          # m, tube inner diameter    0.01
        fin_efficiency                      = self.eta_fin_c           # fin efficiency
        fin_ratio                           = self.fin_ratio_c           # fin ratio

        # Coolant Properties
        rho_c               = 1075
        k_c                 = 0.387
        cp_c                = 3300
        mu_c                = 0.0019

        # HTC_liq                             = 800           # W/m2K, heat transfer coefficient: refrigerant liquid phase
        HTC_2ph                             = 3706          # W/m2K, heat transfer coefficient: refrigerant two phases
        HTC_vap                             = 1076           # W/m2K, heat transfer coefficient: refrigerant vapor phase

        #   ------------------------------------------------------------------------------
        # mass flow rate of refrigerant in single tube
        mfr_channel                         = self.m_dot_refrigerant / Number_of_tubes

        # obtaining overall heat transfer coefficient, U, for each zone
        U_2ph                               = 1 / (1 / HTC_2ph + 1 / (self.HTC_coolant * fin_efficiency * fin_ratio))
        U_vap                               = 1 / (1 / HTC_vap + 1 / (self.HTC_coolant * fin_efficiency * fin_ratio))

        # Pressure drop within the condenser estimation
        f_refrigerant                       = 0.0015   # fractional loss
        u_refrigerant                       = 4 * mfr_channel / (np.pi * ID_tube**2 * self.rho_inlet)    # refrigerant flow velocity
        delta_p_refrigerant                 = f_refrigerant * (ID_tube / L_tube) * (u_refrigerant**2 / 2)  # pressure drop in the tube
        # outlet pressure is set equal to inlet pressure with a refrigerant pump
        P_outlet                            = self.P_inlet - delta_p_refrigerant

        # 3 zones energy balance
        # [1, 2] = [2phase, superheated vapor]
        U_1                                 = U_2ph
        U_2                                 = U_vap

        # estimate the specific heat based on the inlet pressure and enthalpy
        cp_1                                = CoolProp.PropsSI('C', 'P', self.P_inlet, 'H', self.h_inlet, self.refrigerant_name)
        # estimate enthalpy of saturated vapor, x = 1
        h_sat_vap                           = CoolProp.PropsSI('H', 'P', self.P_inlet, 'Q', 1, self.refrigerant_name)
        T_sat_vap                           = CoolProp.PropsSI('T', 'P', self.P_inlet, 'Q', 1, self.refrigerant_name)

        # estimate the specific heat based on the outlet pressure and enthalpy
        T_outlet_design                     = T_sat_vap + self.superheating_degree
        cp_3                                = CoolProp.PropsSI('C', 'P', P_outlet, 'T', T_outlet_design, self.refrigerant_name)

        """algorithm to get 3-zone info;  
        5 thermophysical outputs are: [pressure, temperature, enthalpy, vapor quality, transferred heat]; and
        3 zone lengths are: [L_1, L_2, L_3] """

        # if T_sat_vap > T_air, then the vapor will condense;
        # first to determine the characteristic temperature at the two-phase zone
        T_2phase                        = (self.T_inlet + T_sat_vap) / 2
        delta_h_2phase                  = abs(self.h_inlet - h_sat_vap)
        # calculate 2phase length
        self.L_1                        = self.getL(T_2phase, delta_h_2phase, U_1, self.T_coolant_inlet, mfr_channel, ID_tube)
        # calculate superheated vapor length
        self.L_2                        = L_tube - self.L_1

        if self.L_2 < 0:
            # only 1 zones exist
            self.L_1                = L_tube
            self.L_2                = 0.
            h_1                     = self.getHout(T_2phase, self.L_1, U_1, self.T_coolant_inlet, mfr_channel, ID_tube, self.h_inlet)
            T_1                     = CoolProp.PropsSI('T', 'P', P_outlet, 'H', h_1, self.refrigerant_name)
            # output parameters; scenario - 1
            self.T_outlet           = T_1
            self.P_outlet           = P_outlet
            self.h_outlet           = h_1
            self.vapor_quality_outlet = CoolProp.PropsSI('Q', 'P', P_outlet, 'H', self.h_outlet, self.refrigerant_name)
            self.Q_evaporator       = self.m_dot_refrigerant * (h_1 - self.h_inlet)

        else:
            T_1                     = self.getT(self.L_2, U_2, T_sat_vap, cp_3, self.T_coolant_inlet, mfr_channel, ID_tube)
            h_1                     = CoolProp.PropsSI('H', 'P', P_outlet, 'T', T_1, self.refrigerant_name)
            # output parameters; scenario - 2
            self.T_outlet           = T_1
            self.P_outlet           = P_outlet
            self.h_outlet           = h_1
            self.vapor_quality_outlet = 1.0
            self.Q_evaporator        = self.m_dot_refrigerant * (h_1 - self.h_inlet)

        # calculate air outlet temperature
        # cp_coolant                      = CoolProp.PropsSI('C', 'P', 101325, 'T', self.T_coolant_inlet, self.heat_transfer_fluid)
        C_coolant                       = cp_c * self.m_dot_coolant
        self.T_coolant_out              = self.T_coolant_inlet - self.Q_evaporator / C_coolant
        # condenser overall conductance and effectiveness
        self.UA_tot                     = self.L_2 / L_tube * U_2 + self.L_1 / L_tube * U_1
        self.effectiveness              = self.Q_evaporator / (C_coolant * (self.T_coolant_inlet - self.T_inlet))

        evaporator_results          = results_Data_VCS_evaporator_Finned_Tube_2Zones()
        evaporator_results.T_outlet                       = self.T_outlet
        evaporator_results.P_outlet                       = self.P_outlet
        evaporator_results.h_outlet                       = self.h_outlet
        evaporator_results.x_inlet                        = self.vapor_quality_inlet
        evaporator_results.x_outlet                       = self.vapor_quality_outlet
        evaporator_results.Q_evaporator                   = self.Q_evaporator
        evaporator_results.T_coolant_out                  = self.T_coolant_out
        evaporator_results.UA_tot                         = self.UA_tot
        evaporator_results.effectiveness                  = self.effectiveness

        evaporator_results.L_tot                          = L_tube
        evaporator_results.L_1                            = self.L_1     # 2phase zone
        evaporator_results.L_2                            = self.L_2     # superheated zone

        return evaporator_results

    def getL(self, T_char, dh, U, T_air, mfr_channel, ID_tube):
        dT                          = abs(T_char - T_air)
        L                           = (mfr_channel * dh) / (dT * U * np.pi * ID_tube)
        return L

    def getHout(self, Tchar, length, U, T_air, mfr_channel, ID, h_in):
        UA                          = U * (np.pi * ID * length)
        dT                          = (Tchar - T_air)
        h                           = h_in - UA * dT / mfr_channel
        return h

    def getT(self, L_zone, U, T_in, cp, T_coolant, mfr_channel, ID_tube):
        # get outlet temperature of a zone
        UA                          = U * (np.pi * ID_tube * L_zone)
        T_outlet                    = (UA * (T_coolant - 0.5 * T_in) + mfr_channel * cp * T_in) / (mfr_channel * cp + 0.5 * UA)
        # cannot be lower than the air temperature
        if (T_outlet - T_coolant) > 0:
            T_outlet                = T_coolant

        return T_outlet


# data storage
class results_Data_VCS_evaporator_Finned_Tube_2Zones:
    def __init__(self):
        self.T_outlet                       = []
        self.P_outlet                       = []
        self.h_outlet                       = []
        self.x_inlet                        = []
        self.x_outlet                       = []
        self.Q_evaporator                   = []
        self.T_coolant_out                  = []
        self.UA_tot                         = []
        self.effectiveness                  = []

        self.L_tot                          = []
        self.L_1                            = []
        self.L_2                            = []
        self.L_3                            = []

    def save_data(self, filename):
        with open(filename, 'wb') as f:
            pickle.dump(self, f)


class results_Data_Battery_SinglePhaseCooler:
    def __init__(self):
        self.T_bat_updated  = []
        self.T_o            = []
        self.Q_convec       = []
        self.eff_WavyChan   = []
        self.P_pump_CHAN_coolant = []

    def save_data(self, filename):
        with open(filename, 'wb') as f:
            pickle.dump(self, f)


# wavychannel evaporator
class Battery_Wavychannel_Evaporator_2Zones:
    """
    https://arc.aiaa.org/doi/pdf/10.2514/1.C037404

    A battery pack consists of 150 s x 130 p; the pack has 10 'sheet'; each 'sheet' has 15 wavy-channel flow path
    so that m_dot_chan = m_dot_h / (15 * 10)
    .
    Inputs:
        refName:        Refrigerant fluid
        T_air:          K,      indoor (room) air temperature
        m_dot:          kg/s,   refrigerant flow rate
        m_dot_liq:      kg/s,   liq coolant side flow rate
        P_in:           Pq,     fluid pressure in the inlet
        T_in:           K,      fluid temperature in the inlet
        superheating    K,      Degree of superheating exiting the condenser
    """

    def __init__(self, refrigerant, dt, Q_bat, T_bat, m_dot_ref, P_i, h_i, SH_degree):
        # working fluid types
        self.refrigerant_name               = refrigerant
        # battery side parameters
        self.T_bat                          = T_bat                     # K, inlet coolant temperature              273
        self.Q_bat                          = Q_bat
        self.dt                             = dt
        # refrigerant side parameters
        self.superheating_degree            = SH_degree                 # 11.2 K, level of refrigerant superheating
        self.m_dot_refrigerant              = m_dot_ref                 # kg/s, refrigerant mass flow rate      0.02
        # inlet refrigerant states
        self.P_inlet                        = P_i                       # Pa, inlet pressure
        self.h_inlet                        = h_i                       # K,  inlet temperature
        self.T_inlet                        = CoolProp.PropsSI('T', 'P', self.P_inlet, 'H', self.h_inlet, self.refrigerant_name)
        self.vapor_quality_inlet            = CoolProp.PropsSI('Q', 'P', self.P_inlet, 'H', self.h_inlet, self.refrigerant_name)
        self.rho_inlet                      = CoolProp.PropsSI('D', 'P', self.P_inlet, 'H', self.h_inlet, self.refrigerant_name)

        # outlet refrigerant states
        self.P_outlet                       = None
        self.T_outlet                       = None
        self.h_outlet                       = None
        self.vapor_quality_outlet           = None
        self.Q_evaporator                   = None
        # 3 zone lengths
        self.L_1                            = None
        self.L_2                            = None
        self.L_3                            = None
        # outlet HTF(air) states
        self.T_coolant_out                  = None
        self.UA_tot                         = None
        self.effectiveness                  = None

        # battery
        self.tag                            = "LIB_NMC_18650"
        self.diameter                       = 0.018  # [m]
        self.height                         = 0.065  # [m]
        self.mass                           = 0.048  # [kg]
        self.surface_area                   = (np.pi * self.height * self.diameter) + (0.5 * np.pi * self.diameter ** 2)  # [m^2]
        self.volume                         = np.pi * (0.5 * self.diameter) ** 2 * self.height
        self.density                        = self.mass / self.volume  # [kg/m^3]
        self.electrode_area                 = 0.0342  # [m^2]

        self.max_voltage                    = 4.2  # [V]
        self.nominal_capacity               = 3.55  # [Amp-Hrs]
        self.nominal_voltage                = 3.6  # [V]
        self.charging_voltage               = self.nominal_voltage  # [V]

        self.watt_hour_rating               = self.nominal_capacity * self.nominal_voltage  # [Watt-hours]
        self.specific_energy                = self.watt_hour_rating * 3600 / self.mass  # [J/kg]
        self.specific_power                 = self.specific_energy / self.nominal_capacity  # [W/kg]  need to be required
        self.resistance                     = 0.025  # [Ohms]

        self.specific_heat_capacity         = 1108  # [J/kgK]
        self.radial_thermal_conductivity    = 0.4  # [J/kgK]
        self.axial_thermal_conductivity     = 32.2  # [J/kgK] # estimated

    def thermodynamics_calculation(self):

        """ source: Zhao C, Sousa A C M, Jiang F. Minimization of thermal non-uniformity in lithium-ion battery pack
        cooled by channeled liquid flow[J]. International journal of heat and mass transfer, 2019, 129: 660-670."""
        Theta_inlet                         = 44.5  # contact angle at inlet
        Theta_outlet                        = 50.5
        # wavy channel parameters
        a                                   = 1e-3
        b                                   = 5e-4
        c                                   = 6.3e-2
        d                                   = 2e-3  # todo modified to microchannel
        # hydraulic diameter
        d_H                                 = 2 * (c * d) / (c + d)
        # channel aspect ratio
        gamma_chan                          = d / c
        # channel cross-section area
        A_c_chan                            = c * d
        # channel perimeter
        peri_chan                           = 2 * (c + d)

        rho_chan                            = 2719
        k_chan                              = 202.4
        """
        https://arc.aiaa.org/doi/pdf/10.2514/1.C037404
        A battery pack consists of 150 s x 130 p; the pack has 10 'sheet'; each 'sheet' has 15 wavy-channel flow path
        so that m_dot_chan = m_dot_h / (15 * 10); the length of each channel is 13 * 5 = 65
        """
        # a single pipe with two rows of batteries: 13 * 5 * 2 = 130 cells
        A_surf = (Theta_inlet + Theta_outlet) / 2 / 360 * np.pi * self.diameter * self.height * 65 * 2
        L_surf = (Theta_inlet + Theta_outlet) / 2 / 360 * np.pi * (self.diameter + b + d / 2) * 65 * 2

        L_tube                              = L_surf / 2    # m, tube length
        Number_of_tubes                     = 15 * 10       # number of tubes in parallel
        ID_tube                             = d_H           # m, tube inner diameter

        # mass flow rate of refrigerant in single tube
        mfr_channel                         = self.m_dot_refrigerant / Number_of_tubes
        # u_refrigerant                       = 4 * mfr_channel / (A_c_chan * self.rho_inlet)    # refrigerant flow velocity, why *4?
        u_refrigerant                       = mfr_channel / (A_c_chan * self.rho_inlet)    # refrigerant flow velocity
        # print("channel flow rate", mfr_channel, " kg/s")
        # ------------------------------------------------------------------------------
        # Calculate heat transfer coefficient
        # ------------------------------------------------------------------------------
        # HTC_liq                             = 600           # W/m2K, heat transfer coefficient: refrigerant liquid phase
        HTC_vap                             = 200         # W/m2K, heat transfer coefficient: refrigerant vapor phase
        # f, h, Re                             = f_h_1phase_Channel(mfr_channel, d_H, u_refrigerant, T_sat_vap, self.P_inlet, Phase='Single')
        HTC_2ph                             = 3000        # W/m2K, heat transfer coefficient: refrigerant two phases
        #   ------------------------------------------------------------------------------

        # obtaining overall heat transfer coefficient, U, for each zone
        U_2ph                               = 1 / (1 / HTC_2ph + b / k_chan)
        U_vap                               = 1 / (1 / HTC_vap + b / k_chan)

        # Pressure drop within the condenser estimation
        f_refrigerant                       = 0.0015   # fractional loss
        delta_p_refrigerant                 = f_refrigerant * (ID_tube / L_tube) * (u_refrigerant**2 / 2)  # pressure drop in the tube
        # outlet pressure is set equal to inlet pressure with a refrigerant pump
        P_outlet                            = self.P_inlet - delta_p_refrigerant

        # 3 zones energy balance
        # [1, 2] = [2phase, superheated vapor]
        U_1                                 = U_2ph
        U_2                                 = U_vap

        # estimate the specific heat based on the inlet pressure and enthalpy
        cp_1                                = CoolProp.PropsSI('C', 'P', self.P_inlet, 'H', self.h_inlet, self.refrigerant_name)
        # estimate enthalpy of saturated vapor, x = 1
        h_sat_vap                           = CoolProp.PropsSI('H', 'P', self.P_inlet, 'Q', 1, self.refrigerant_name)
        T_sat_vap                           = CoolProp.PropsSI('T', 'P', self.P_inlet, 'Q', 1, self.refrigerant_name)

        # estimate the specific heat based on the outlet pressure and enthalpy
        T_outlet_design                     = T_sat_vap + self.superheating_degree
        cp_3                                = CoolProp.PropsSI('C', 'P', P_outlet, 'T', T_outlet_design, self.refrigerant_name)

        """algorithm to get 3-zone info;  
        5 thermophysical outputs are: [pressure, temperature, enthalpy, vapor quality, transferred heat]; and
        3 zone lengths are: [L_1, L_2, L_3] """

        # if T_sat_vap > T_air, then the vapor will condense;
        # first to determine the characteristic temperature at the two-phase zone
        T_2phase                        = (self.T_inlet + T_sat_vap) / 2
        delta_h_2phase                  = abs(self.h_inlet - h_sat_vap)
        # calculate 2phase length
        self.L_1                        = self.getL(T_2phase, delta_h_2phase, U_1, self.T_bat, mfr_channel, peri_chan)
        # calculate superheated vapor length
        self.L_2                        = L_tube - self.L_1

        if self.L_2 < 0:
            # only 1 zones exist
            self.L_1                = L_tube
            self.L_2                = 0.
            h_1                     = self.getHout(T_2phase, self.L_1, U_1, self.T_bat, mfr_channel, peri_chan, self.h_inlet)
            T_1                     = CoolProp.PropsSI('T', 'P', P_outlet, 'H', h_1, self.refrigerant_name)
            # output parameters; scenario - 1
            self.T_outlet           = T_1
            self.P_outlet           = P_outlet
            self.h_outlet           = h_1
            self.vapor_quality_outlet = CoolProp.PropsSI('Q', 'P', P_outlet, 'H', self.h_outlet, self.refrigerant_name)
            self.Q_evaporator       = self.m_dot_refrigerant * (h_1 - self.h_inlet)

        else:
            T_1                     = self.getT(self.L_2, U_2, T_sat_vap, cp_3, self.T_bat, mfr_channel, ID_tube)
            h_1                     = CoolProp.PropsSI('H', 'P', P_outlet, 'T', T_1, self.refrigerant_name)
            # output parameters; scenario - 2
            self.T_outlet           = T_1
            self.P_outlet           = P_outlet
            self.h_outlet           = h_1
            self.vapor_quality_outlet = 1.0
            self.Q_evaporator        = self.m_dot_refrigerant * (h_1 - self.h_inlet)

        # calculate updated battery temperature
        C_bat                       = self.mass * self.specific_heat_capacity * 150 * 130
        self.T_bat_updated              = self.T_bat + (self.Q_bat - self.Q_evaporator) * self.dt * 60 / C_bat
        # condenser overall conductance and effectiveness
        self.UA_tot                     = self.L_2 / L_tube * U_2 + self.L_1 / L_tube * U_1
        # self.effectiveness              = self.Q_evaporator / (C_coolant * (self.T_coolant_inlet - self.T_inlet))

        evaporator_results          = results_Data_VCS_evaporator_Finned_Tube_2Zones()
        evaporator_results.T_outlet                       = self.T_outlet
        evaporator_results.P_outlet                       = self.P_outlet
        evaporator_results.h_outlet                       = self.h_outlet
        evaporator_results.x_inlet                        = self.vapor_quality_inlet
        evaporator_results.x_outlet                       = self.vapor_quality_outlet
        evaporator_results.Q_evaporator                   = self.Q_evaporator
        evaporator_results.T_coolant_out                  = self.T_coolant_out
        evaporator_results.UA_tot                         = self.UA_tot
        # evaporator_results.effectiveness                  = self.effectiveness
        evaporator_results.T_bat_updated                  = self.T_bat_updated
        evaporator_results.L_tot                          = L_tube
        evaporator_results.L_1                            = self.L_1     # 2phase zone
        evaporator_results.L_2                            = self.L_2     # superheated zone

        return evaporator_results

    def getL(self, T_char, dh, U, T_air, mfr_channel, peri_chan):
        dT                          = abs(T_char - T_air)
        L                           = (mfr_channel * dh) / (dT * U * peri_chan)
        return L

    def getHout(self, Tchar, length, U, T_air, mfr_channel, peri_chan, h_in):
        UA                          = U * (peri_chan * length)
        dT                          = (Tchar - T_air)
        h                           = h_in - UA * dT / mfr_channel
        return h

    def getT(self, L_zone, U, T_in, cp, T_coolant, mfr_channel, peri_chan):
        # get outlet temperature of a zone
        UA                          = U * (peri_chan * L_zone)
        T_outlet                    = (UA * (T_coolant - 0.5 * T_in) + mfr_channel * cp * T_in) / (mfr_channel * cp + 0.5 * UA)
        # cannot be lower than the air temperature
        if (T_outlet - T_coolant) > 0:
            T_outlet                = T_coolant

        return T_outlet

    def weight_cal(self):
        # wavy channel parameters
        a                                   = 1e-3
        b                                   = 5e-4
        c                                   = 6.3e-2
        d                                   = 2e-3  # todo modified to microchannel
        Theta_inlet                         = 44.5  # contact angle at inlet
        Theta_outlet                        = 50.5

        A_c                                 = 2 * (a * d + b * (c + 2 * a))     # m2
        rho_chan                            = 2719                              # kg/m3
        L_surf = (Theta_inlet + Theta_outlet) / 2 / 360 * np.pi * (self.diameter + b + d / 2) * 65 * 2 # m

        L_tube                              = L_surf / 2        # m, tube length
        Number_of_tubes                     = 15 * 10       # number of tubes in parallel

        # weight (dry, without refrigerant)
        m_evap_dry                          = rho_chan * Number_of_tubes * L_tube * A_c

        return m_evap_dry


# wavychannel single phase cooler
class Battery_Wavychannel_SinglePhaseCooler:
    """
    https://arc.aiaa.org/doi/pdf/10.2514/1.C037404

    A battery pack consists of 150 s x 130 p; the pack has 10 'sheet'; each 'sheet' has 15 wavy-channel flow path
    so that m_dot_chan = m_dot_h / (15 * 10)
    .
    Inputs:
        refName:        Refrigerant fluid
        T_air:          K,      indoor (room) air temperature
        m_dot:          kg/s,   refrigerant flow rate
        m_dot_liq:      kg/s,   liq coolant side flow rate
        P_in:           Pq,     fluid pressure in the inlet
        T_in:           K,      fluid temperature in the inlet
        superheating    K,      Degree of superheating exiting the condenser
    """

    def __init__(self, dt, Q_bat, T_bat, m_dot_c, T_i):
        # battery side parameters
        self.T_bat                          = T_bat                     # K, inlet coolant temperature              273
        self.T_i                            = T_i                       # coolant inlet temperature for wavychannel
        self.Q_bat                          = Q_bat
        self.dt                             = dt
        # refrigerant side parameters
        self.m_dot_c                        = m_dot_c                   # kg/s, coolant mass flow rate      0.02

        # outlet refrigerant states
        self.P_outlet                       = None
        self.T_outlet                       = None
        self.h_outlet                       = None
        self.vapor_quality_outlet           = None
        self.Q_evaporator                   = None
        # 3 zone lengths
        self.L_1                            = None
        self.L_2                            = None
        self.L_3                            = None
        # outlet HTF(air) states
        self.T_coolant_out                  = None
        self.UA_tot                         = None
        self.effectiveness                  = None

        # battery
        self.tag                            = "LIB_NMC_18650"
        self.diameter                       = 0.018  # [m]
        self.height                         = 0.065  # [m]
        self.mass                           = 0.048  # [kg]
        self.surface_area                   = (np.pi * self.height * self.diameter) + (0.5 * np.pi * self.diameter ** 2)  # [m^2]
        self.volume                         = np.pi * (0.5 * self.diameter) ** 2 * self.height
        self.density                        = self.mass / self.volume  # [kg/m^3]
        self.electrode_area                 = 0.0342  # [m^2]

        self.max_voltage                    = 4.2  # [V]
        self.nominal_capacity               = 3.55  # [Amp-Hrs]
        self.nominal_voltage                = 3.6  # [V]
        self.charging_voltage               = self.nominal_voltage  # [V]

        self.watt_hour_rating               = self.nominal_capacity * self.nominal_voltage  # [Watt-hours]
        self.specific_energy                = self.watt_hour_rating * 3600 / self.mass  # [J/kg]
        self.specific_power                 = self.specific_energy / self.nominal_capacity  # [W/kg]  need to be required
        self.resistance                     = 0.025  # [Ohms]

        self.specific_heat_capacity         = 1108  # [J/kgK]
        self.radial_thermal_conductivity    = 0.4  # [J/kgK]
        self.axial_thermal_conductivity     = 32.2  # [J/kgK] # estimated

    def thermodynamics_calculation(self):

        """ source: Zhao C, Sousa A C M, Jiang F. Minimization of thermal non-uniformity in lithium-ion battery pack
        cooled by channeled liquid flow[J]. International journal of heat and mass transfer, 2019, 129: 660-670."""
        Theta_inlet                         = 44.5  # contact angle at inlet
        Theta_outlet                        = 50.5
        # wavy channel parameters
        a                                   = 1e-3
        b                                   = 5e-4
        c                                   = 6.3e-2
        d                                   = 2e-3  # todo modified to microchannel
        # hydraulic diameter
        d_H                                 = 2 * (c * d) / (c + d)
        # channel aspect ratio
        gamma_chan                          = d / c
        # channel cross-section area
        A_c_chan                            = c * d
        # channel perimeter
        peri_chan                           = 2 * (c + d)

        rho_chan                            = 2719
        k_chan                              = 202.4
        """
        https://arc.aiaa.org/doi/pdf/10.2514/1.C037404
        A battery pack consists of 150 s x 130 p; the pack has 10 'sheet'; each 'sheet' has 15 wavy-channel flow path
        so that m_dot_chan = m_dot_h / (15 * 10); the length of each channel is 13 * 5 = 65
        """
        # a single pipe with two rows of batteries: 13 * 5 * 2 = 130 cells
        A_surf = (Theta_inlet + Theta_outlet) / 2 / 360 * np.pi * self.diameter * self.height * 65 * 2
        L_surf = (Theta_inlet + Theta_outlet) / 2 / 360 * np.pi * (self.diameter + b + d / 2) * 65 * 2

        L_tube                              = L_surf / 2    # m, tube length
        Number_of_tubes                     = 15 * 10       # number of tubes in parallel
        ID_tube                             = d_H           # m, tube inner diameter

        # Coolant Properties
        rho_c               = 1075
        k_c                 = 0.387
        cp_c                = 3300
        mu_c                = 0.0019
        # mass flow rate of refrigerant in single tube
        mfr_channel                         = self.m_dot_c / Number_of_tubes
        u_coolant                           = mfr_channel / (A_c_chan * rho_c)    # refrigerant flow velocity
        # print("channel flow rate", mfr_channel, " kg/s")

        # Pressure drop within the condenser estimation
        f_coolant                           = 0.0015   # fractional loss
        delta_p_coolant                     = f_coolant * (ID_tube / L_tube) * (u_coolant**2 / 2)  # pressure drop in the tube
        eff_pump                            = 0.7
        P_pump_CHAN_coolant                 = mfr_channel * delta_p_coolant / rho_c / eff_pump

        Re_coolant                          = rho_c * u_coolant * d_H / mu_c
        Pr_coolant                          = cp_c * mu_c / k_c

        if Re_coolant <= 2300:
          # Nusselt number
            Nu = 8.235 * (1 - 2.0421 * gamma_chan + 3.0853 * np.power(gamma_chan, 2) - 2.4765 * np.power(gamma_chan, 3)
                          + 1.0578 * np.power(gamma_chan, 4) - 0.1861 * np.power(gamma_chan, 5))
            # Colburn factor
            j_coolant = Nu / Re_coolant * np.power(Pr_coolant, -1 / 3)
        else:
            # Colburn factor
            f_coolant = 1 / (4 * (1.8 * np.log10(Re_coolant / 7.7)) ** 2)
            # use Gnielinski equation to calculate Nu for 0.5 < Pr < 2000, 3000 < Re < 5e6
            Nu = (f_coolant / 2) * (Re_coolant - 1000) * Pr_coolant / (
                    1 + 12.7 * np.power((f_coolant / 2), 0.5) * (np.power(Pr_coolant, 2 / 3) - 1))
            # Colburn factor
            j_coolant = Nu / Re_coolant * np.power(Pr_coolant, -1 / 3)

        # heat transfer coefficient of channeled coolant fluid
        h_coolant = k_c * Nu / d_H
        # total heat transfer coefficient
        h_tot = 1 / (1 / h_coolant + b / k_chan)

        # number of transfer units and effectiveness
        NTU                     = h_tot * A_surf / (mfr_channel * cp_c)
        self.eff_WavyChan            = 1 - np.exp(-NTU)

        # log mean temperature
        Tw_Ti                   = self.T_bat - self.T_i
        Tw_To                   = Tw_Ti * np.exp(-NTU)
        dT_lm                   = (Tw_Ti - Tw_To) / np.log(Tw_Ti / Tw_To)

        # convective heat from the battery
        self.Q_convec                = h_tot * A_surf * dT_lm * Number_of_tubes

        # check the heat generated
        self.T_o                     = self.T_bat - Tw_To
        Q_conv_check            = mfr_channel * cp_c * (self.T_o - self.T_i) * Number_of_tubes
        delta_Q_conv = np.abs(self.Q_convec - Q_conv_check)

        # check the wavy channel effectiveness
        eff_WavyChan_check      = (self.T_o - self.T_i) / (self.T_bat - self.T_i)
        delta_eff_WavyChan      = np.abs(self.eff_WavyChan - eff_WavyChan_check)

        # net heat stored in the battery
        P_net                   = self.Q_bat - self.Q_convec

        # temperature rise [update cell temperature]
        dT_dt                   = P_net / (self.mass * self.specific_energy * 150*130)
        self.T_bat_updated      = self.T_bat + dT_dt * self.dt

        battery_results          = results_Data_Battery_SinglePhaseCooler()
        battery_results.T_bat_updated                   = self.T_bat_updated
        battery_results.T_o                             = self.T_o
        battery_results.Q_convec                        = self.Q_convec
        battery_results.eff_WavyChan                    = self.eff_WavyChan
        battery_results.P_pump_CHAN_coolant             = P_pump_CHAN_coolant

        return battery_results

    def weight_cal(self):
        # wavy channel parameters
        a                                   = 1e-3
        b                                   = 5e-4
        c                                   = 6.3e-2
        d                                   = 2e-3  # todo modified to microchannel
        Theta_inlet                         = 44.5  # contact angle at inlet
        Theta_outlet                        = 50.5

        A_c                                 = 2 * (a * d + b * (c + 2 * a))     # m2
        rho_chan                            = 2719                              # kg/m3
        L_surf = (Theta_inlet + Theta_outlet) / 2 / 360 * np.pi * (self.diameter + b + d / 2) * 65 * 2 # m

        L_tube                              = L_surf / 2        # m, tube length
        Number_of_tubes                     = 15 * 10       # number of tubes in parallel

        # weight (dry, without refrigerant)
        m_evap_dry                          = rho_chan * Number_of_tubes * L_tube * A_c

        return m_evap_dry


# -----------------------------------------------------
if __name__ == '__main__':
    condenser   = VCS_condenser_Finned_Tube_3Zones_herringbone('R134a', 'water', 0.5, 5e5, 323, 8.5,308)
    condenser.HerringboneFin()
