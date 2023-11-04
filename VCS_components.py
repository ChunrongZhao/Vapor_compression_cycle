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


# ---------------------------------------------------------------------
#   Condenser
# ---------------------------------------------------------------------
class VCS_condenser_PHE_design:
    "James Bull, et al. Heat Exchanger Sizing for Organic Rankine Cycle"
    def __init__(self):
        "nominal plate geometrical parameters"
        # parameters                            # value             # unit          # symbol
        self.effective_length                   = 0.4               # m             # L
        self.effective_width                    = 0.2               # m             # W
        self.plate_thickness                    = 0.001             # m             # t
        self.port_diameter                      = 0.018             # m             # D_p
        self.vertical_port_distance             = 0.01              # m             # L_v
        self.horizontal_port_distance           = 0.02              # m             # L_c
        self.corrugation_wavelength             = 0.009             # m             # lamda
        self.corrugation_angle                  = 60                # o             # beta
        self.compressed_plate_pack_length       = 0.3               # m             # L_c
        self.surface_enlargement_factor         = 1.25              # -             # phi
        self.thermal_conductivity               = 16                # W/mK          # k (stainless steel)


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

    def __init__(self, refrigerant, coolant, T_coolant_inlet, m_dot_ref, m_dot_coolant, P_i, h_i, SH_degree, HTC_coolant):
        # working fluid types
        self.refrigerant_name               = refrigerant
        self.heat_transfer_fluid            = coolant                   # liquid coolant
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
        cp_coolant                      = CoolProp.PropsSI('C', 'P', 101325, 'T', self.T_coolant_inlet, self.heat_transfer_fluid)
        C_coolant                       = cp_coolant * self.m_dot_coolant
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
