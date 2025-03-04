# VCS_Cycle_Performance.py
#
# Created: Oct 2023, C.R. Zhao

# ---------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import numpy                        as np
import CoolProp.CoolProp            as CoolProp
from matplotlib                     import pyplot as plt
from scipy.optimize                 import fsolve
import pandas                       as pd
from openpyxl.workbook              import Workbook
import pickle
from pathlib                        import Path
from more_itertools                 import chunked
from scipy                          import optimize

from VCS_Codes.Components.VCS_components         import VCS_condenser_Finned_Tube_3Zones, Battery_Wavychannel_Evaporator_2Zones, \
    VCS_condenser_Finned_Tube_3Zones_herringbone
from ACHP_codes.Components.WavyChan import WavyChan_NMC
from matplotlib import rcParams, rc, font_manager
import matplotlib as mpl
import matplotlib.pyplot as plt
import os

import CoolProp as CP
from ACHP_codes.Correlations.HTC_Correlations import f_h_1phase_Tube, ShahCondensation_Average, f_h_1phase_Channel
from ACHP_codes.Correlations.HTC_Correlations import ShahEvaporation_Average
from collections import defaultdict

from ACHP_codes.Correlations.FinStructure_Correlations import HerringboneFins_condenser


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
#   Imports
# ---------------------------------------------------------------------
class Cycle_Performance_BatteryAir:
    """"Heat pump class:
        This Heat pump considers a 3 zone model for the condenser and the evaporator.
        It is suitable for air-air heat pumps.

        Inputs:
            rpm:                rpm, compressor rotation
            eta_C_iso:          isoentropic efficiency
            Disp:               m3,displacement
            eta_C_vol:          volumetric efficiency
            m_dot_w:            kg/s, indoor fluid mass flow rate
            T_8:                K, indoor fluid temperature
            m_dot_air:          kg/s, outdoor fluid mass flow rate
            T_air_in:           K, outdoor fluid temperature
            Condfluid='water':  indoor fluid (water or air)
            fluid='R134a':      refrigerant fluid

        Methods:
            state_:             [Pressure, Temperature, Enthalpy, Quality]
            Q_cond:             W, heat transfered in the condenser
            T_5:                K, fluid temperature exiting condenser
            Q_evap:             W, heat transfered in the evaporator
            T_air_out:          K, fluid temperature exiting evaporator
            P_1:                Pa, pressure in the compressor inlet
            P_2:                Pa, pressure in the compressor outlet
            W_c:                W, compressor work
            R_c:                compression rate
            COP:                real Coefficient of Performance
            E_balance:          residual from calculation
            Capacity:           Btu/h, refrigeration capacity
            """

    def __init__(self, SysParams, dt, Q_bat, T_bat, m_dot_air, T_air_inlet, air, refrigerant):
        self.rpm                        = SysParams['Comp: rpm']
        self.eta_C_iso                  = SysParams['Comp: eta_C_iso']
        self.Displacement               = SysParams['Comp: Disp']
        self.eta_C_vol                  = SysParams['Comp: eta_C_vol']
        self.eta_C_mec                  = SysParams['Comp: eta_C_mec']

        self.superheating               = SysParams['Sys: Superheat']
        self.subcooling                 = SysParams['Sys: Subcool']

        self.W_Fan_COND                 = SysParams['Cond: Fan']

        self.T_bat                      = T_bat
        self.Q_bat                      = Q_bat
        self.dt                         = dt
        self.m_dot_air                  = m_dot_air
        self.T_air_inlet                = T_air_inlet
        self.m_dot_refrigerant          = None
        self.refrigerant                = refrigerant
        self.heat_transfer_fluid_cond   = air   # condenser side; air

        self.dec                        = 10    # decimal for some PropsSI functions
        self.CD_HTC_air                 = 100   # W/m2K, heat transfer coefficient of Condenser: air side
        self.EV_HTC_coolant             = 600   # W/m2K, heat transfer coefficient of Evaporator: liquid coolant side
        self.C_min_evap                 = None
        self.C_min_cond                 = None
        # results
        self.state_1                    = None
        self.state_2                    = None
        self.state_3                    = None
        self.state_4                    = None
        self.P_1                        = None
        self.P_2                        = None

        self.flag1                      = False
        self.Q_evaporator               = None
        self.Q_condenser                = None

    def March_Compressor_Inlet_and_Outlet_Pressures(self, P_guess):
        # Flag True if the pressure guesses are negative
        if P_guess[0] < 0 or P_guess[1] < 0:
            self.flag1                  = True
            return (0, 0)

        else:
            # self.C_min_evap             = CoolProp.PropsSI('C', 'P', 101325, 'T', self.T_coolant_inlet, self.refrigerant) * self.m_dot_coolant
            self.C_min_cond             = CoolProp.PropsSI('C', 'P', 101325, 'T', self.T_air_inlet, self.heat_transfer_fluid_cond) * self.m_dot_air

            # ---------------------------------------------------------
            # Point 1
            # ---------------------------------------------------------
            self.T_1                    = round(CoolProp.PropsSI('T', 'P', P_guess[0], 'Q', 1, self.refrigerant) + self.superheating, self.dec)
            self.h_1                    = CoolProp.PropsSI('H', 'P', P_guess[0], 'T', self.T_1, self.refrigerant)
            self.x_1                    = abs(CoolProp.PropsSI('Q', 'P', P_guess[0], 'T', self.T_1, self.refrigerant))  # Q is the Quality
            self.s_1                    = CoolProp.PropsSI('S', 'P', P_guess[0], 'T', self.T_1, self.refrigerant)
            rho_1                       = CoolProp.PropsSI('D', 'P', P_guess[0], 'T', self.T_1, self.refrigerant)

            self.m_dot_refrigerant      = rho_1 * self.Displacement * (self.rpm / 60) * self.eta_C_vol
            self.state_1                = [P_guess[0], self.h_1, self.T_1, self.s_1, self.x_1]

            # ---------------------------------------------------------
            # Point 2
            # ---------------------------------------------------------
            self.s_2                    = self.s_1
            self.h_2s                   = CoolProp.PropsSI('H', 'P', P_guess[1], 'S', self.s_2, self.refrigerant)
            self.h_2                    = (self.h_2s - self.h_1) / self.eta_C_iso + self.h_1
            self.T_2                    = round(CoolProp.PropsSI('T', 'P', P_guess[1], 'H', self.h_2, self.refrigerant), self.dec)
            self.x_2                    = abs(CoolProp.PropsSI('Q', 'P', P_guess[1], 'T', self.T_2, self.refrigerant))
            self.state_2                = [P_guess[1], self.h_2, self.T_2, self.s_2, self.x_2]
            # ---------------------------------------------------------
            # Point 3
            # ---------------------------------------------------------
            self.T_3                    = round(CoolProp.PropsSI('T', 'P', P_guess[1], 'Q', 0, self.refrigerant) - self.subcooling, self.dec)  # 8.4
            self.h_3                    = CoolProp.PropsSI('H', 'P', P_guess[1], 'T', self.T_3, self.refrigerant)
            self.s_3                    = CoolProp.PropsSI('S', 'P', P_guess[1], 'T', self.T_3, self.refrigerant)
            # Get outlet results of condenser through 3-zone model
            condenser                   = VCS_condenser_Finned_Tube_3Zones(self.refrigerant, self.heat_transfer_fluid_cond, self.T_air_inlet,
                                          self.m_dot_refrigerant, self.m_dot_air, P_guess[1], self.T_2, self.subcooling, self.CD_HTC_air)
            condenser_results           = condenser.thermodynamics_calculation()
            # ---------------------------------------------------------------
            Q_cond                      = condenser_results.Q_condenser
            T_3_cal                     = condenser_results.T_outlet
            P_outlet                    = condenser_results.P_outlet
            h_outlet                    = condenser_results.h_outlet
            x_outlet                    = condenser_results.x_outlet
            T_air_out                   = condenser_results.T_air_out
            UA_tot                      = condenser_results.UA_tot
            effectiveness               = condenser_results.effectiveness
            self.cond_L_tot                       = condenser_results.L_tot
            self.cond_L_1                         = condenser_results.L_1
            self.cond_L_2                         = condenser_results.L_2
            self.cond_L_3                         = condenser_results.L_3
            # ---------------------------------------------------------------
            h_3_cal                     = Q_cond / self.m_dot_refrigerant + self.h_2
            T_3_cal                     = T_3_cal

            self.T_air_outlet           = T_air_out
            self.x_3                    = x_outlet
            self.state_3                = [P_guess[1], self.h_3, self.T_3, self.s_3, self.x_3]
            # ---------------------------------------------------------------
            # Point 4
            self.h_4                    = self.h_3
            self.T_4                    = round(CoolProp.PropsSI('T', 'P', P_guess[0], 'H', self.h_4, self.refrigerant), self.dec)
            self.x_4                    = CoolProp.PropsSI('Q', 'P', P_guess[0], 'H', self.h_4, self.refrigerant)
            # Get outlet results of condenser through 2-zone model, assuming that the inlet refrigerant is in 2phase state
            evaporator                   = Battery_Wavychannel_Evaporator_2Zones(self.refrigerant, self.dt, self.Q_bat, self.T_bat, self.m_dot_refrigerant, round(P_guess[0], 2), self.h_4, self.superheating)
            evaporator_results           = evaporator.thermodynamics_calculation()
            # ---------------------------------------------------------------
            Q_evap                      = evaporator_results.Q_evaporator
            T_1_cal                     = evaporator_results.T_outlet
            P_outlet                    = evaporator_results.P_outlet
            h_outlet                    = evaporator_results.h_outlet
            x_inlet                     = evaporator_results.x_inlet
            x_outlet                    = evaporator_results.x_outlet
            UA_tot                      = evaporator_results.UA_tot
            effectiveness               = evaporator_results.effectiveness
            self.evap_L_tot                       = evaporator_results.L_tot
            self.evap_L_1                         = evaporator_results.L_1
            self.evap_L_2                         = evaporator_results.L_2
            # ---------------------------------------------------------------
            h_1_cal                     = Q_evap / self.m_dot_refrigerant + self.h_4
            T_1_cal                     = T_1_cal
            self.T_bat_updated          = evaporator_results.T_bat_updated
            self.s_4                    = CoolProp.PropsSI('S', 'P', P_guess[0], 'Q', abs(self.x_4), self.refrigerant)
            self.state_4                = [P_guess[0], self.h_4, self.T_4, self.s_4, self.x_4]
            # ---------------------------------------------------------
            # store parameters
            self.Q_evaporator           = Q_evap
            self.Q_condenser            = Q_cond
            # ---------------------------------------------------------
            errP1                       = (self.h_1 - h_1_cal)
            errP2                       = (self.h_3 - h_3_cal)
            return (errP1, errP2)

    def Cycle_Solver(self):
        # ---------------------------------------------------------
        # Initial guess considering the saturation temperatures with superheating and subcooling
        P_1_initial                     = int(CoolProp.PropsSI('P', 'Q', 1, 'T', self.T_bat - self.subcooling, self.refrigerant))
        P_2_initial                     = int(CoolProp.PropsSI('P', 'Q', 1, 'T', self.T_air_inlet + self.superheating, self.refrigerant))

        P_guess                         = np.array([P_1_initial, P_2_initial])

        self.flag1                      = False  # error flag when guess is P<0
        result                          = optimize.root(self.March_Compressor_Inlet_and_Outlet_Pressures, P_guess, method='hybr')  # Minimization function
        self.P_1                        = result.x[0]
        self.P_2                        = result.x[1]

    def PostProcessing(self):
        # Post Processing
        # ---------------------------------------------------------
        self.W_c                         = self.m_dot_refrigerant * (self.h_2 - self.h_1)
        self.R_c                         = self.P_2 / self.P_1
        # refrigerating COP
        self.COP                         = self.Q_evaporator / ((self.W_c / self.eta_C_mec) + self.W_Fan_COND)

        self.E_balance                   = self.Q_condenser + self.Q_evaporator + self.W_c
        # 1W = 3.41 Btu/h
        self.Capacity                    = - self.Q_evaporator * 3.41

        print('Phase: ', CoolProp.PhaseSI('P', self.P_1, 'Q', self.x_4, self.refrigerant))
        self.s_4                            = CoolProp.PropsSI('S', 'P', self.P_1, 'Q', self.x_4, self.refrigerant)

        self.s_4_liq                        = CoolProp.PropsSI('S', 'P', self.P_1, 'Q', 0, self.refrigerant)
        self.s_4_vap                        = CoolProp.PropsSI('S', 'P', self.P_1, 'Q', 1, self.refrigerant)
        print(self.x_4)
        print(self.s_4_vap, self.s_4_liq, (self.s_4_vap * self.x_4 + self.s_4_liq * (1 - self.x_4)))

        return {'W_c': self.W_c, 'R_c': self.R_c, 'COP': self.COP, 'E_balance': self.E_balance, 'Capacity': self.Capacity}

    def Interpolation_Singularity_Points(self):
        # Linear interpolation to predict values in singularity points
        # ---------------------------------------------------------
        if self.flag1 == False:
            res1 = self.PostProcessing()
            self.W_c                    = res1['W_c']
            self.R_c                    = res1['R_c']
            self.COP                    = res1['COP']
            self.E_balance              = res1['E_balance']
            self.Capacity               = res1['Capacity']
        else:
            dif                         = 1  # value away from the singularity to perform the interpolation
            # ---------------------Simulation with T_amb - 1 K--------------------------------
            P_1_initial                 = int(CoolProp.PropsSI('P', 'Q', 1, 'T', self.T_bat - self.subcooling, self.refrigerant))
            P_2_initial                 = int(CoolProp.PropsSI('P', 'Q', 1, 'T', self.T_air_inlet - dif + self.superheating, self.refrigerant))

            P_guess                     = np.array([P_1_initial, P_2_initial])

            result2                     = optimize.root(self.March_Compressor_Inlet_and_Outlet_Pressures, P_guess, method='hybr')  # Minimization function
            self.P_1                    = result2.x[0]
            self.P_2                    = result2.x[1]
            res2                        = self.PostProcessing()

            # --------------------Simulation with T_amb + 1 K-------------------------------------
            P_1_initial                 = int(CoolProp.PropsSI('P', 'Q', 1, 'T', self.T_bat - self.subcooling, self.refrigerant))
            P_2_initial                 = int(CoolProp.PropsSI('P', 'Q', 1, 'T', self.T_air_inlet + dif + self.superheating, self.refrigerant))

            P_guess                     = np.array([P_1_initial, P_2_initial])

            result3                     = optimize.root(self.March_Compressor_Inlet_and_Outlet_Pressures, P_guess, method='hybr')  # Minimization function
            self.P_1                    = result3.x[0]
            self.P_2                    = result3.x[1]
            res3                        = self.PostProcessing()

            self.W_c                    = (res2['W_c'] + res3['W_c']) / 2
            self.R_c                    = (res2['R_c'] + res3['R_c']) / 2
            self.COP                    = (res2['COP'] + res3['COP']) / 2
            self.E_balance              = (res2['E_balance'] + res3['E_balance']) / 2
            self.Capacity               = (res2['Capacity'] + res3['Capacity']) / 2


class Cycle_Performance_BatteryAir_herringbone:
    """"Heat pump class:
        This Heat pump considers a 3 zone model for the condenser and the evaporator.
        It is suitable for air-air heat pumps.

        Inputs:
            rpm:                rpm, compressor rotation
            eta_C_iso:          isoentropic efficiency
            Disp:               m3,displacement
            eta_C_vol:          volumetric efficiency
            m_dot_w:            kg/s, indoor fluid mass flow rate
            T_8:                K, indoor fluid temperature
            m_dot_air:          kg/s, outdoor fluid mass flow rate
            T_air_in:           K, outdoor fluid temperature
            Condfluid='water':  indoor fluid (water or air)
            fluid='R134a':      refrigerant fluid

        Methods:
            state_:             [Pressure, Temperature, Enthalpy, Quality]
            Q_cond:             W, heat transfered in the condenser
            T_5:                K, fluid temperature exiting condenser
            Q_evap:             W, heat transfered in the evaporator
            T_air_out:          K, fluid temperature exiting evaporator
            P_1:                Pa, pressure in the compressor inlet
            P_2:                Pa, pressure in the compressor outlet
            W_c:                W, compressor work
            R_c:                compression rate
            COP:                real Coefficient of Performance
            E_balance:          residual from calculation
            Capacity:           Btu/h, refrigeration capacity
            """

    def __init__(self, SysParams, dt, Q_bat, T_bat, m_dot_air, T_air_inlet, eta_fin, fin_ratio, air, refrigerant):
        self.rpm                        = SysParams['Comp: rpm']
        self.eta_C_iso                  = SysParams['Comp: eta_C_iso']
        self.Displacement               = SysParams['Comp: Disp']
        self.eta_C_vol                  = SysParams['Comp: eta_C_vol']
        self.eta_C_mec                  = SysParams['Comp: eta_C_mec']

        self.superheating               = SysParams['Sys: Superheat']
        self.subcooling                 = SysParams['Sys: Subcool']

        self.W_Fan_COND                 = SysParams['Cond: Fan']

        self.T_bat                      = T_bat
        self.Q_bat                      = Q_bat
        self.dt                         = dt
        self.m_dot_air                  = m_dot_air
        self.T_air_inlet                = T_air_inlet
        self.m_dot_refrigerant          = None
        self.refrigerant                = refrigerant
        self.heat_transfer_fluid_cond   = air   # condenser side; air

        self.eta_fin                    = eta_fin
        self.fin_ratio                  = fin_ratio

        self.dec                        = 10    # decimal for some PropsSI functions
        self.CD_HTC_air                 = 100   # W/m2K, heat transfer coefficient of Condenser: air side
        self.EV_HTC_coolant             = 600   # W/m2K, heat transfer coefficient of Evaporator: liquid coolant side
        self.C_min_evap                 = None
        self.C_min_cond                 = None
        # results
        self.state_1                    = None
        self.state_2                    = None
        self.state_3                    = None
        self.state_4                    = None
        self.P_1                        = None
        self.P_2                        = None

        self.flag1                      = False
        self.Q_evaporator               = None
        self.Q_condenser                = None

    def March_Compressor_Inlet_and_Outlet_Pressures(self, P_guess):
        # Flag True if the pressure guesses are negative
        if P_guess[0] < 0 or P_guess[1] < 0:
            self.flag1                  = True
            return (0, 0)

        else:
            # self.C_min_evap             = CoolProp.PropsSI('C', 'P', 101325, 'T', self.T_coolant_inlet, self.refrigerant) * self.m_dot_coolant
            self.C_min_cond             = CoolProp.PropsSI('C', 'P', 101325, 'T', self.T_air_inlet, self.heat_transfer_fluid_cond) * self.m_dot_air

            # ---------------------------------------------------------
            # Point 1
            # ---------------------------------------------------------
            self.T_1                    = round(CoolProp.PropsSI('T', 'P', P_guess[0], 'Q', 1, self.refrigerant) + self.superheating, self.dec)
            self.h_1                    = CoolProp.PropsSI('H', 'P', P_guess[0], 'T', self.T_1, self.refrigerant)
            self.x_1                    = abs(CoolProp.PropsSI('Q', 'P', P_guess[0], 'T', self.T_1, self.refrigerant))  # Q is the Quality
            self.s_1                    = CoolProp.PropsSI('S', 'P', P_guess[0], 'T', self.T_1, self.refrigerant)
            rho_1                       = CoolProp.PropsSI('D', 'P', P_guess[0], 'T', self.T_1, self.refrigerant)

            self.m_dot_refrigerant      = rho_1 * self.Displacement * (self.rpm / 60) * self.eta_C_vol
            self.state_1                = [P_guess[0], self.h_1, self.T_1, self.s_1, self.x_1]

            # ---------------------------------------------------------
            # Point 2
            # ---------------------------------------------------------
            self.s_2                    = self.s_1
            self.h_2s                   = CoolProp.PropsSI('H', 'P', P_guess[1], 'S', self.s_2, self.refrigerant)
            self.h_2                    = (self.h_2s - self.h_1) / self.eta_C_iso + self.h_1
            self.T_2                    = round(CoolProp.PropsSI('T', 'P', P_guess[1], 'H', self.h_2, self.refrigerant), self.dec)
            self.x_2                    = abs(CoolProp.PropsSI('Q', 'P', P_guess[1], 'T', self.T_2, self.refrigerant))
            self.state_2                = [P_guess[1], self.h_2, self.T_2, self.s_2, self.x_2]
            # ---------------------------------------------------------
            # Point 3
            # ---------------------------------------------------------
            self.T_3                    = round(CoolProp.PropsSI('T', 'P', P_guess[1], 'Q', 0, self.refrigerant) - self.subcooling, self.dec)  # 8.4
            self.h_3                    = CoolProp.PropsSI('H', 'P', P_guess[1], 'T', self.T_3, self.refrigerant)
            self.s_3                    = CoolProp.PropsSI('S', 'P', P_guess[1], 'T', self.T_3, self.refrigerant)
            # Get outlet results of condenser through 3-zone model
            condenser                   = VCS_condenser_Finned_Tube_3Zones_herringbone(self.refrigerant, self.heat_transfer_fluid_cond, self.T_air_inlet,
                                          self.m_dot_refrigerant, self.m_dot_air, P_guess[1], self.T_2, self.subcooling, self.CD_HTC_air, self.eta_fin, self.fin_ratio)
            condenser_results           = condenser.thermodynamics_calculation()
            # ---------------------------------------------------------------
            Q_cond                      = condenser_results.Q_condenser
            T_3_cal                     = condenser_results.T_outlet
            P_outlet                    = condenser_results.P_outlet
            h_outlet                    = condenser_results.h_outlet
            x_outlet                    = condenser_results.x_outlet
            T_air_out                   = condenser_results.T_air_out
            UA_tot                      = condenser_results.UA_tot
            effectiveness               = condenser_results.effectiveness
            self.cond_L_tot                       = condenser_results.L_tot
            self.cond_L_1                         = condenser_results.L_1
            self.cond_L_2                         = condenser_results.L_2
            self.cond_L_3                         = condenser_results.L_3
            # ---------------------------------------------------------------
            h_3_cal                     = Q_cond / self.m_dot_refrigerant + self.h_2
            T_3_cal                     = T_3_cal

            self.T_air_outlet           = T_air_out
            self.x_3                    = x_outlet
            self.state_3                = [P_guess[1], self.h_3, self.T_3, self.s_3, self.x_3]
            # ---------------------------------------------------------------
            # Point 4
            self.h_4                    = self.h_3
            self.T_4                    = round(CoolProp.PropsSI('T', 'P', P_guess[0], 'H', self.h_4, self.refrigerant), self.dec)
            self.x_4                    = CoolProp.PropsSI('Q', 'P', P_guess[0], 'H', self.h_4, self.refrigerant)
            # Get outlet results of condenser through 2-zone model, assuming that the inlet refrigerant is in 2phase state
            evaporator                   = Battery_Wavychannel_Evaporator_2Zones(self.refrigerant, self.dt, self.Q_bat, self.T_bat, self.m_dot_refrigerant, round(P_guess[0], 2), self.h_4, self.superheating)
            evaporator_results           = evaporator.thermodynamics_calculation()
            # ---------------------------------------------------------------
            Q_evap                      = evaporator_results.Q_evaporator
            T_1_cal                     = evaporator_results.T_outlet
            P_outlet                    = evaporator_results.P_outlet
            h_outlet                    = evaporator_results.h_outlet
            x_inlet                     = evaporator_results.x_inlet
            x_outlet                    = evaporator_results.x_outlet
            UA_tot                      = evaporator_results.UA_tot
            effectiveness               = evaporator_results.effectiveness
            self.evap_L_tot                       = evaporator_results.L_tot
            self.evap_L_1                         = evaporator_results.L_1
            self.evap_L_2                         = evaporator_results.L_2
            # ---------------------------------------------------------------
            h_1_cal                     = Q_evap / self.m_dot_refrigerant + self.h_4
            T_1_cal                     = T_1_cal
            self.T_bat_updated          = evaporator_results.T_bat_updated
            self.s_4                    = CoolProp.PropsSI('S', 'P', P_guess[0], 'Q', abs(self.x_4), self.refrigerant)
            self.state_4                = [P_guess[0], self.h_4, self.T_4, self.s_4, self.x_4]
            # ---------------------------------------------------------
            # store parameters
            self.Q_evaporator           = Q_evap
            self.Q_condenser            = Q_cond
            # ---------------------------------------------------------
            errP1                       = (self.h_1 - h_1_cal)
            errP2                       = (self.h_3 - h_3_cal)
            return (errP1, errP2)

    def Cycle_Solver(self):
        # ---------------------------------------------------------
        # Initial guess considering the saturation temperatures with superheating and subcooling
        P_1_initial                     = int(CoolProp.PropsSI('P', 'Q', 1, 'T', self.T_bat - self.subcooling, self.refrigerant))
        P_2_initial                     = int(CoolProp.PropsSI('P', 'Q', 1, 'T', self.T_air_inlet + self.superheating, self.refrigerant))

        P_guess                         = np.array([P_1_initial, P_2_initial])
        print(P_guess)
        self.flag1                      = False  # error flag when guess is P<0
        result                          = optimize.root(self.March_Compressor_Inlet_and_Outlet_Pressures, P_guess, method='hybr')  # Minimization function
        self.P_1                        = result.x[0]
        self.P_2                        = result.x[1]
        print(self.P_1, "Pa")

    def PostProcessing(self):
        # Post Processing
        # ---------------------------------------------------------
        self.W_c                         = self.m_dot_refrigerant * (self.h_2 - self.h_1)
        self.R_c                         = self.P_2 / self.P_1
        # refrigerating COP
        self.COP                         = self.Q_evaporator / ((self.W_c / self.eta_C_mec) + self.W_Fan_COND)

        self.E_balance                   = self.Q_condenser + self.Q_evaporator + self.W_c
        # 1W = 3.41 Btu/h
        self.Capacity                    = - self.Q_evaporator * 3.41

        # print('Phase: ', CoolProp.PhaseSI('P', self.P_1, 'Q', self.x_4, self.refrigerant))
        self.s_4                            = CoolProp.PropsSI('S', 'P', self.P_1, 'Q', self.x_4, self.refrigerant)

        self.s_4_liq                        = CoolProp.PropsSI('S', 'P', self.P_1, 'Q', 0, self.refrigerant)
        self.s_4_vap                        = CoolProp.PropsSI('S', 'P', self.P_1, 'Q', 1, self.refrigerant)
        # print(self.x_4)
        # print(self.s_4_vap, self.s_4_liq, (self.s_4_vap * self.x_4 + self.s_4_liq * (1 - self.x_4)))

        return {'W_c': self.W_c, 'R_c': self.R_c, 'COP': self.COP, 'E_balance': self.E_balance, 'Capacity': self.Capacity}

    def Interpolation_Singularity_Points(self):
        # Linear interpolation to predict values in singularity points
        # ---------------------------------------------------------
        if self.flag1 == False:
            res1 = self.PostProcessing()
            self.W_c                    = res1['W_c']
            self.R_c                    = res1['R_c']
            self.COP                    = res1['COP']
            self.E_balance              = res1['E_balance']
            self.Capacity               = res1['Capacity']
        else:
            dif                         = 1  # value away from the singularity to perform the interpolation
            # ---------------------Simulation with T_amb - 1 K--------------------------------
            P_1_initial                 = int(CoolProp.PropsSI('P', 'Q', 1, 'T', self.T_bat - self.subcooling, self.refrigerant))
            P_2_initial                 = int(CoolProp.PropsSI('P', 'Q', 1, 'T', self.T_air_inlet - dif + self.superheating, self.refrigerant))

            P_guess                     = np.array([P_1_initial, P_2_initial])

            result2                     = optimize.root(self.March_Compressor_Inlet_and_Outlet_Pressures, P_guess, method='hybr')  # Minimization function
            self.P_1                    = result2.x[0]
            self.P_2                    = result2.x[1]
            res2                        = self.PostProcessing()

            # --------------------Simulation with T_amb + 1 K-------------------------------------
            P_1_initial                 = int(CoolProp.PropsSI('P', 'Q', 1, 'T', self.T_bat - self.subcooling, self.refrigerant))
            P_2_initial                 = int(CoolProp.PropsSI('P', 'Q', 1, 'T', self.T_air_inlet + dif + self.superheating, self.refrigerant))

            P_guess                     = np.array([P_1_initial, P_2_initial])

            result3                     = optimize.root(self.March_Compressor_Inlet_and_Outlet_Pressures, P_guess, method='hybr')  # Minimization function
            self.P_1                    = result3.x[0]
            self.P_2                    = result3.x[1]
            res3                        = self.PostProcessing()

            self.W_c                    = (res2['W_c'] + res3['W_c']) / 2
            self.R_c                    = (res2['R_c'] + res3['R_c']) / 2
            self.COP                    = (res2['COP'] + res3['COP']) / 2
            self.E_balance              = (res2['E_balance'] + res3['E_balance']) / 2
            self.Capacity               = (res2['Capacity'] + res3['Capacity']) / 2


def HerringboneFinCal(V_dot_ha, T_air_inlet, p=101325, RH=0.51):
        Fins        = AttributeDict()
        Fins.Tubes.N_Tubes_per_bank     = 5                    # number of tubes per bank=row
        Fins.Tubes.N_bank               = 1                     # number of banks/rows
        Fins.Tubes.L_tube               = 10
        Fins.Tubes.OD                   = 0.011
        Fins.Tubes.ID                   = 0.01
        Fins.Tubes.Pl                   = 0.0191                # distance between center of tubes in flow direction
        Fins.Tubes.Pt                   = 0.0254                # distance between center of tubes orthogonal to flow direction
        Fins.Tubes.k_w                  = 237                   # wall thermal conductivity (i.e. pipe material)

        Fins.Fins.FPI                   = 8                    # Number of fins per inch
        Fins.Fins.Pd                    = 0.001                 # 2* amplitude of wavy fin
        Fins.Fins.xf                    = 0.001                 # 1/2 period of fin
        Fins.Fins.t                     = 0.00011               # Thickness of fin material
        Fins.Fins.k_fin                 = 237                   # Thermal conductivity of fin material

        Fins.Air.V_dot_ha               = V_dot_ha         # rated volumetric flowrate
        Fins.Air.T_mean                 = T_air_inlet
        Fins.Air.T_db                   = T_air_inlet
        Fins.Air.p                      = p                # Condenser Air pressures in Pa
        Fins.Air.RH                     = RH
        Fins.Air.RH_mean                = RH
        Fins.Air.FanPower               = 260

        A, eta_o, h_a, DeltaP_air, m_dot_da, m_dot_ha, fin_ratio = HerringboneFins_condenser(Fins)

        m_dot_air_cal              = m_dot_ha
        A_fin                      = A
        eta_fin                    = eta_o
        HTC_air_cal                = h_a

        # print(m_dot_air_cal, A_fin, eta_fin, HTC_air_cal, fin_ratio)

        return m_dot_air_cal, A_fin, eta_fin, HTC_air_cal, fin_ratio


# ---------------------------------------------------------------------
#   Functions
# ---------------------------------------------------------------------
def BatteryData(i=80, DT_amb=35):
    address                         = os.getcwd()
    with open(address + '\BatInputs_for_VCS.pkl', 'rb') as f:
        BatInputs                   = pickle.load(f)

    # print(BatInputs['Q_bat'])

    WavyChan                        = WavyChan_NMC()
    # results_chan                    = WavyChan_data()
    T_air                           = BatInputs['T_amb'][i] + DT_amb
    # T_bat                           = BatInputs['T_b_air'][i] + DT_amb
    Q_bat                           = BatInputs['Q_bat'][i] * 10
    dt                              = BatInputs['time'][i] - BatInputs['time'][i-1]

    # print(Q_bat, 'W', T_bat, 'K', T_air, 'K')

    return T_air, Q_bat, dt


def VCC_weight_estimation(W_comp, m_dot_ref):
    # evaporator, wavy channel
    m_evap_dry          = evap_weight()

    # compressor
    m_comp              = 0.3448 * W_comp + 6.4655

    # condenser
    m_cond_dry          = cond_weight()

    # expansion valve
    # density at 20 oC; 1225.5 kg/m3 (liq) 27.773 kg/m3 (vap)
    m_throttle          = 0.026 * m_dot_ref / 1225.5 * 1e3 / 60 + 0.73996

    m_VCC               = m_evap_dry + m_comp + m_cond_dry + m_throttle

    return m_VCC


def evap_weight():
    # wavy channel parameters
    a                                   = 1e-3
    b                                   = 5e-4
    c                                   = 6.3e-2
    d                                   = 2e-3
    Theta_inlet                         = 44.5  # contact angle at inlet
    Theta_outlet                        = 50.5

    A_c                                 = 2 * (a * d + b * (c + 2 * a))     # m2
    rho_chan                            = 2719                              # kg/m3
    L_surf = (Theta_inlet + Theta_outlet) / 2 / 360 * np.pi * (0.018 + b + d / 2) * 65 * 2 # m

    L_tube                              = L_surf        # m, tube length
    Number_of_tubes                     = 15 * 10       # number of tubes in parallel

    # weight (dry, without refrigerant)
    m_evap_dry                          = rho_chan * Number_of_tubes * L_tube * A_c

    return m_evap_dry


def cond_weight():
    # Geometric parameters, finned-tube structure
    L_tube                              = 10            # m, tube length
    Number_of_tubes                     = 5             # number of tubes in parallel
    ID_tube                             = 0.01          # m, tube inner diameter    0.01
    OD_tube                             = 0.011

    # condenser dry
    m_cond_dry                          = np.pi * (OD_tube**2 - ID_tube**2) * L_tube * Number_of_tubes

    return m_cond_dry


# -----HTC ---------------
# 2->3
def HTC_cond_2ph_comparison(Tsat, m_dot_r, Number_of_tubes=5, D=0.01, HTC=3000):
    x                           = np.linspace(0, 1, 1000)
    h                           = np.zeros_like(x)
    AS                          = CP.AbstractState('HEOS', 'R134a')

    # --------------------------------------------------------
    G = m_dot_r / Number_of_tubes / (np.pi * D**2 / 4)

    AS.update(CP.QT_INPUTS, 0.0, Tsat)     # saturated liq
    p_sat_r                     = AS.p()  # [Pa]

    AS.update(CP.PQ_INPUTS, p_sat_r, 0.0)
    h_l                         = AS.hmass()                # [J/kg]
    AS.update(CP.PQ_INPUTS, p_sat_r, 1.0)
    h_v                         = AS.hmass()                # [J/kg]

    h_fg                        = h_v - h_l                 # [J/kg]

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

    for i in range(len(x)):
        h[i] = ShahCondensation_Average(x[i], x[i], AS, G, D, p_sat_r, Tsat, Tsat)

    h_2ph_ave = np.trapz(h, x=x)     # to integrate along the axis
    err_h_2ph_cond  = (h_2ph_ave - HTC) / h_2ph_ave * 100
    # print('HTC_2phase_cond error: ', err_h_2ph_cond, '%')

    return h_2ph_ave, err_h_2ph_cond


def HTC_cond_vap_comparison(Tsat, T_in_r, m_dot_r, Number_of_tubes=5, D=0.01, HTC=600):
    AS                          = CP.AbstractState('HEOS', 'R134a')
    # --------------------------------------------------------
    m_dot               = m_dot_r / Number_of_tubes

    AS.update(CP.QT_INPUTS, 0.0, Tsat)     # saturated liq
    p_sat_r                     = AS.p()  # [Pa]

    T_vap_ave           = (Tsat + T_in_r) / 2

    f_vap, h_vap_ave, Re_vap           = f_h_1phase_Tube(m_dot,D, T_vap_ave, p_sat_r, AS)
    err_h_vap_cond      = (h_vap_ave - HTC) / h_vap_ave * 100
    # print('HTC_superheat_vap_cond error: ', err_h_vap_cond, '%')

    return h_vap_ave, err_h_vap_cond


def HTC_cond_liq_comparison(Tsat, T_out_r, m_dot_r, Number_of_tubes=5, D=0.01, HTC=600):
    AS                          = CP.AbstractState('HEOS', 'R134a')
    # --------------------------------------------------------
    m_dot               = m_dot_r / Number_of_tubes

    AS.update(CP.QT_INPUTS, 1.0, Tsat)     # saturated liq
    p_sat_r                     = AS.p()  # [Pa]

    T_liq_ave           = (Tsat + T_out_r) / 2

    f_liq, h_liq_ave, Re_liq           = f_h_1phase_Tube(m_dot, D, T_liq_ave, p_sat_r, AS)
    err_h_liq_cond      = (h_liq_ave - HTC) / h_liq_ave * 100
    # print('HTC_subcooled_liq_cond error: ', err_h_liq_cond, '%')

    return h_liq_ave, err_h_liq_cond


# 4->1
def HTC_evap_2ph_comparison(Tsat, m_dot_r, x_in, L_2ph, Number_of_tubes=150, HTC=3000):
    x                           = np.linspace(x_in, 1, 1000)
    h                           = np.zeros_like(x)
    AS                          = CP.AbstractState('HEOS', 'R134a')

    # --------------------------------------------------------
    # wavy channel parameters
    a                                   = 1e-3
    b                                   = 5e-4
    c                                   = 6.3e-2
    d                                   = 2e-3  # todo modified to microchannel
    A_c                         =  c * d
    G                           = m_dot_r / Number_of_tubes / A_c
    m_dot                       = m_dot_r / Number_of_tubes

    AS.update(CP.QT_INPUTS, 0.0, Tsat)     # saturated liq
    p_sat_r                     = AS.p()  # [Pa]

    AS.update(CP.PQ_INPUTS, p_sat_r, 0.0)
    h_l                         = AS.hmass()                # [J/kg]
    AS.update(CP.PQ_INPUTS, p_sat_r, 1.0)
    h_v                         = AS.hmass()                # [J/kg]

    h_fg                        = h_v - h_l                 # [J/kg]
    Q_targe                     = m_dot * (1.0 - x_in) * h_fg   # [W]
    Ar                          = 2 * (c + d) * L_2ph
    q_flux                      = Q_targe / Ar              # [W/m^2]

    # value of W/cm^2
    Q_flux_cm = q_flux /1e4
    # hydraulic diameter
    d_H                                 = 2 * (c * d) / (c + d)
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

    for i in range(len(x)):
        h[i] = ShahEvaporation_Average(x[i], x[i], AS, G, d_H, p_sat_r, q_flux, Tsat, Tsat)

    h_2ph_ave = np.trapz(h, x=x)     # to integrate along the axis
    err_h_2ph_evap  = (h_2ph_ave - HTC) / h_2ph_ave * 100
    # print('HTC_2phase_evap error: ', err_h_2ph_evap, '%')

    return h_2ph_ave, err_h_2ph_evap


def HTC_evap_vap_comparison(Tsat, T_out_r, m_dot_r, Number_of_tubes=150, HTC=200):
    AS                          = CP.AbstractState('HEOS', 'R134a')
    # --------------------------------------------------------
    m_dot                               = m_dot_r / Number_of_tubes

    AS.update(CP.QT_INPUTS, 0.0, Tsat)     # saturated liq
    p_sat_r                             = AS.p()  # [Pa]
    T_vap_ave                           = (Tsat + T_out_r) / 2

    # wavy channel parameters
    a                                   = 1e-3
    b                                   = 5e-4
    c                                   = 6.3e-2
    d                                   = 2e-3  # todo modified to microchannel

    f_vap, h_vap_ave, Re_vap            = f_h_1phase_Channel(m_dot, d, c, T_vap_ave, p_sat_r, AS, Phase='Single')
    err_h_vap_evap                      = (h_vap_ave - HTC) / h_vap_ave * 100
    # print('HTC_superheat_vap_evap error: ', err_h_vap_evap, '%')

    return h_vap_ave, err_h_vap_evap


# air flow HTC
def HTC_cond_airside_finned_tube(HTC_air_cal, HTC=100):
    err_h_air_cond      = (HTC_air_cal - HTC) / HTC_air_cal * 100
    # print('HTC_superheat_vap_evap error: ', err_h_air_cond, '%')

    return HTC_air_cal, err_h_air_cond


# ---------------------------------------------------------------------
#   Data storage
# ---------------------------------------------------------------------
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


def save_file(self, filename):
    with open(filename, 'wb') as f:
        pickle.dump(self, f)


def plot_BTMS_results(BTMS_results, width=10, height=9):
    # ---------------------------------------------------------------------
    #   set the font style
    # ---------------------------------------------------------------------
    # Let's set the font to the one we want
    font_path = ['/Users/Chunrong/PycharmProjects/Library/Fonts' ]
    font_files = font_manager.findSystemFonts(fontpaths=font_path)

    for font_file in font_files:
        font_manager.fontManager.addfont(font_file)
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    mpl.rcParams['font.family'] = 'Gulliver-Regular'
    # ---------------------------------------------------------------------
    #   todo flags need to be played with
    # ---------------------------------------------------------------------
    set_ticks = False
    set_axis_limits = False
    save_figure = False
    # ---------------------------------------------------------------------
    #   main plot settings
    # ---------------------------------------------------------------------
    line_width = 2.25
    line_width_soft = 2
    marker_size_exp = 7
    marker_size_exp = str(marker_size_exp)
    marker_size = 8
    marker_size = str(marker_size)
    fontsize = 14
    fontsize_legend = 12
    axes_linewidth = 1.75
    alpha = 1.0
    # journals usually want 600 dpi or more but that means large file sizes. I recommend at least 300
    image_resolution = 300  # in dots per inch (dpi).

    mpl.rc('axes', linewidth=axes_linewidth, labelsize=fontsize, titlesize=fontsize)
    mpl.rc('xtick', labelsize=fontsize)
    mpl.rc('ytick', labelsize=fontsize)
    mpl.rc('legend', fontsize=fontsize_legend)
    mpl.rc('savefig', dpi=image_resolution, format='png', bbox='tight')
    # ---------------------------------------------------------------------
    #   set the plot colors
    # ---------------------------------------------------------------------
    # defining the colors
    myblack = [0, 0, 0]
    myblue = '#0F95D7'
    myred = '#e41a1c'
    myyellow = [255 / 255, 194 / 255, 10 / 255]
    mygreen = '#4daf4a'
    mybrown = '#a65628'
    mydarkblue = '#377eb8'
    mypurple = '#984ea3'
    myorange = '#ff7f00'
    mygray = [89 / 255, 89 / 255, 89 / 255]

    # setting up color and marker sequence
    # note there are only 10 but you should not have more than 10 lines on a plot
    colors = [myblack, myblue, myred, myyellow, mydarkblue, myorange, mybrown, mygreen, mypurple, mygray]
    markers = ['o', '^', 's', 'p', 'v', '*', 'x']

    # ---------------------------------------------------------------------
    #   load results
    # ---------------------------------------------------------------------
    blob_size                           = 100 /1                 # for solutions in plot
    x_label                             = 'Time (min)'
    y1_label                            = '$T_{bat}$ ($^\circ$C)'
    y2_label                            = '$L_{2ph,evap} / L_{tot}$ (%)'
    y3_label                            = 'COP (-)'
    y4_label                            = '$W_{comp}$ (kW)'
    y5_label                            = '$Q_{bat}$ (kW)'
    y6_label                            = '$Q_{conv}$ (kW)'
    # ---------------------------------------------------------------------
    #   plot results
    # ---------------------------------------------------------------------
    fig1 = plt.figure(None)
    fig1.set_size_inches(width, height)
    fig1.suptitle(None)

    plt.subplot(3, 2, 1)
    plt.scatter(BTMS_results['mission_time'][:-10], BTMS_results['T_bat_liq'][:-10], marker=markers[0], alpha=alpha, s=blob_size, c='none', edgecolors=colors[0], linewidths=line_width, label='$T_{bat,vcs}$')
    plt.scatter(BTMS_results['mission_time'][:-10], BTMS_results['T_bat_air'][:-10], marker=markers[0], alpha=alpha, s=blob_size, c='none', edgecolors=colors[1], linewidths=line_width, label='$T_{bat,air}$')
    plt.grid(visible=True, which='major', linestyle='-', linewidth=0.75, color=colors[0])
    plt.minorticks_on()
    plt.grid(visible=True, which='minor', linestyle='--', linewidth=0.15, color=colors[0])
    plt.ylabel(y1_label)
    plt.legend()

    plt.subplot(3, 2, 2)
    plt.scatter(BTMS_results['mission_time'][:-10], BTMS_results['L_2ph'][:-10], marker=markers[1], alpha=alpha, s=blob_size, c='none', edgecolors=colors[2], linewidths=line_width)
    plt.grid(visible=True, which='major', linestyle='-', linewidth=0.75, color=colors[0])
    plt.minorticks_on()
    plt.grid(visible=True, which='minor', linestyle='--', linewidth=0.15, color=colors[0])
    plt.ylabel(y2_label)

    plt.subplot(3, 2, 3)
    plt.scatter(BTMS_results['mission_time'][:-10], BTMS_results['COP'][:-10], marker=markers[3], alpha=alpha, s=blob_size, c='none', edgecolors=colors[7], linewidths=line_width)
    plt.grid(visible=True, which='major', linestyle='-', linewidth=0.75, color=colors[0])
    plt.minorticks_on()
    plt.grid(visible=True, which='minor', linestyle='--', linewidth=0.15, color=colors[0])
    plt.xlabel(x_label)
    plt.ylabel(y3_label)

    plt.subplot(3, 2, 4)
    plt.scatter(BTMS_results['mission_time'][:-10], BTMS_results['Power'][:-10], marker=markers[3], alpha=alpha, s=blob_size, c='none', edgecolors=colors[8], linewidths=line_width)
    plt.grid(visible=True, which='major', linestyle='-', linewidth=0.75, color=colors[0])
    plt.minorticks_on()
    plt.grid(visible=True, which='minor', linestyle='--', linewidth=0.15, color=colors[0])
    plt.xlabel(x_label)
    plt.ylabel(y4_label)

    plt.subplot(3, 2, 5)
    plt.scatter(BTMS_results['mission_time'][:-10], BTMS_results['Q_bat'][:-10], marker=markers[3], alpha=alpha, s=blob_size, c='none', edgecolors=colors[4], linewidths=line_width)
    plt.grid(visible=True, which='major', linestyle='-', linewidth=0.75, color=colors[0])
    plt.minorticks_on()
    plt.grid(visible=True, which='minor', linestyle='--', linewidth=0.15, color=colors[0])
    plt.xlabel(x_label)
    plt.ylabel(y5_label)

    plt.subplot(3, 2, 6)
    plt.scatter(BTMS_results['mission_time'][:-10], BTMS_results['Q_conv'][:-10], marker=markers[3], alpha=alpha, s=blob_size, c='none', edgecolors=colors[5], linewidths=line_width)
    plt.grid(visible=True, which='major', linestyle='-', linewidth=0.75, color=colors[0])
    plt.minorticks_on()
    plt.grid(visible=True, which='minor', linestyle='--', linewidth=0.15, color=colors[0])
    plt.xlabel(x_label)
    plt.ylabel(y6_label)

    plt.tight_layout()
    plt.show()
    # plt.savefig('fig10.png', dpi=image_resolution, format='png')


def plot_BTMS_results2(BTMS_results, width=10, height=9):
    # ---------------------------------------------------------------------
    #   set the font style
    # ---------------------------------------------------------------------
    # Let's set the font to the one we want
    font_path = ['/Users/Chunrong/PycharmProjects/Library/Fonts' ]
    font_files = font_manager.findSystemFonts(fontpaths=font_path)

    for font_file in font_files:
        font_manager.fontManager.addfont(font_file)
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    mpl.rcParams['font.family'] = 'Gulliver-Regular'
    # ---------------------------------------------------------------------
    #   todo flags need to be played with
    # ---------------------------------------------------------------------
    set_ticks = False
    set_axis_limits = False
    save_figure = False
    # ---------------------------------------------------------------------
    #   main plot settings
    # ---------------------------------------------------------------------
    line_width = 2.25
    line_width_soft = 2
    marker_size_exp = 7
    marker_size_exp = str(marker_size_exp)
    marker_size = 8
    marker_size = str(marker_size)
    fontsize = 14
    fontsize_legend = 12
    axes_linewidth = 1.75
    alpha = 1.0
    # journals usually want 600 dpi or more but that means large file sizes. I recommend at least 300
    image_resolution = 300  # in dots per inch (dpi).

    mpl.rc('axes', linewidth=axes_linewidth, labelsize=fontsize, titlesize=fontsize)
    mpl.rc('xtick', labelsize=fontsize)
    mpl.rc('ytick', labelsize=fontsize)
    mpl.rc('legend', fontsize=fontsize_legend)
    mpl.rc('savefig', dpi=image_resolution, format='png', bbox='tight')
    # ---------------------------------------------------------------------
    #   set the plot colors
    # ---------------------------------------------------------------------
    # defining the colors
    myblack = [0, 0, 0]
    myblue = '#0F95D7'
    myred = '#e41a1c'
    myyellow = [255 / 255, 194 / 255, 10 / 255]
    mygreen = '#4daf4a'
    mybrown = '#a65628'
    mydarkblue = '#377eb8'
    mypurple = '#984ea3'
    myorange = '#ff7f00'
    mygray = [89 / 255, 89 / 255, 89 / 255]

    # setting up color and marker sequence
    # note there are only 10 but you should not have more than 10 lines on a plot
    colors = [myblack, myblue, myred, myyellow, mydarkblue, myorange, mybrown, mygreen, mypurple, mygray]
    markers = ['o', '^', 's', 'p', 'v', '*', 'x']

    # ---------------------------------------------------------------------
    #   load results
    # ---------------------------------------------------------------------
    blob_size                           = 100 /1                 # for solutions in plot
    x_label                             = 'Time (min)'
    y1_label                            = '$T_{bat}$ ($^\circ$C)'
    y2_label                            = '$L_{2ph,evap} / L_{tot}$ (%)'
    y3_label                            = 'COP (-)'
    y4_label                            = '$W_{comp}$ (kW)'
    y5_label                            = '$Q_{bat}$ (kW)'
    y6_label                            = '$Q_{conv}$ (kW)'
    y7_label                            = '$\dot{m}_{ref}$ (kg/s)'
    y8_label                            = '$m_{VCC}$ (kg)'
    y9_label                            = '$P_1$ (kPa)'
    y10_label                            = '$T_{sat,1}$ (K)'
    y11_label                            = '$P_2$ (kPa)'
    y12_label                            = '$T_{sat,2}$ (K)'
    # ---------------------------------------------------------------------
    #   plot results
    # ---------------------------------------------------------------------
    fig1 = plt.figure(None)
    fig1.set_size_inches(width, height)
    fig1.suptitle(None)

    plt.subplot(3, 2, 1)
    plt.scatter(BTMS_results['mission_time'][:-10], BTMS_results['m_dot_ref'][:-10], marker=markers[3], alpha=alpha, s=blob_size, c='none', edgecolors=colors[0], linewidths=line_width)
    plt.grid(visible=True, which='major', linestyle='-', linewidth=0.75, color=colors[0])
    plt.minorticks_on()
    plt.grid(visible=True, which='minor', linestyle='--', linewidth=0.15, color=colors[0])
    plt.ylabel(y7_label)

    plt.subplot(3, 2, 2)
    plt.scatter(BTMS_results['mission_time'][:-10], BTMS_results['m_VCC'][:-10], marker=markers[3], alpha=alpha, s=blob_size, c='none', edgecolors=colors[1], linewidths=line_width)
    plt.grid(visible=True, which='major', linestyle='-', linewidth=0.75, color=colors[0])
    plt.minorticks_on()
    plt.grid(visible=True, which='minor', linestyle='--', linewidth=0.15, color=colors[0])
    plt.ylabel(y8_label)

    plt.subplot(3, 2, 3)
    plt.scatter(BTMS_results['mission_time'][:-10], BTMS_results['P1'][:-10], marker=markers[3], alpha=alpha, s=blob_size, c='none', edgecolors=colors[2], linewidths=line_width)
    plt.grid(visible=True, which='major', linestyle='-', linewidth=0.75, color=colors[0])
    plt.minorticks_on()
    plt.grid(visible=True, which='minor', linestyle='--', linewidth=0.15, color=colors[0])
    plt.ylabel(y9_label)

    plt.subplot(3, 2, 4)
    plt.scatter(BTMS_results['mission_time'][:-10], BTMS_results['T1_sat'][:-10], marker=markers[3], alpha=alpha, s=blob_size, c='none', edgecolors=colors[3], linewidths=line_width)
    plt.grid(visible=True, which='major', linestyle='-', linewidth=0.75, color=colors[0])
    plt.minorticks_on()
    plt.grid(visible=True, which='minor', linestyle='--', linewidth=0.15, color=colors[0])
    plt.ylabel(y10_label)

    plt.subplot(3, 2, 5)
    plt.scatter(BTMS_results['mission_time'][:-10], BTMS_results['P2'][:-10], marker=markers[3], alpha=alpha, s=blob_size, c='none', edgecolors=colors[4], linewidths=line_width)
    plt.grid(visible=True, which='major', linestyle='-', linewidth=0.75, color=colors[0])
    plt.minorticks_on()
    plt.grid(visible=True, which='minor', linestyle='--', linewidth=0.15, color=colors[0])
    plt.xlabel(x_label)
    plt.ylabel(y11_label)

    plt.subplot(3, 2, 6)
    plt.scatter(BTMS_results['mission_time'][:-10], BTMS_results['T2_sat'][:-10], marker=markers[3], alpha=alpha, s=blob_size, c='none', edgecolors=colors[5], linewidths=line_width)
    plt.grid(visible=True, which='major', linestyle='-', linewidth=0.75, color=colors[0])
    plt.minorticks_on()
    plt.grid(visible=True, which='minor', linestyle='--', linewidth=0.15, color=colors[0])
    plt.xlabel(x_label)
    plt.ylabel(y12_label)

    plt.tight_layout()
    plt.show()
    # plt.savefig('fig10.png', dpi=image_resolution, format='png')


def plot_BTMS_results3(BTMS_results, width=10, height=9):
    # ---------------------------------------------------------------------
    #   set the font style
    # ---------------------------------------------------------------------
    # Let's set the font to the one we want
    font_path = ['/Users/Chunrong/PycharmProjects/Library/Fonts' ]
    font_files = font_manager.findSystemFonts(fontpaths=font_path)

    for font_file in font_files:
        font_manager.fontManager.addfont(font_file)
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    mpl.rcParams['font.family'] = 'Gulliver-Regular'
    # ---------------------------------------------------------------------
    #   todo flags need to be played with
    # ---------------------------------------------------------------------
    set_ticks = False
    set_axis_limits = False
    save_figure = False
    # ---------------------------------------------------------------------
    #   main plot settings
    # ---------------------------------------------------------------------
    line_width = 2.25
    line_width_soft = 2
    marker_size_exp = 7
    marker_size_exp = str(marker_size_exp)
    marker_size = 8
    marker_size = str(marker_size)
    fontsize = 14
    fontsize_legend = 12
    axes_linewidth = 1.75
    alpha = 1.0
    # journals usually want 600 dpi or more but that means large file sizes. I recommend at least 300
    image_resolution = 300  # in dots per inch (dpi).

    mpl.rc('axes', linewidth=axes_linewidth, labelsize=fontsize, titlesize=fontsize)
    mpl.rc('xtick', labelsize=fontsize)
    mpl.rc('ytick', labelsize=fontsize)
    mpl.rc('legend', fontsize=fontsize_legend)
    mpl.rc('savefig', dpi=image_resolution, format='png', bbox='tight')
    # ---------------------------------------------------------------------
    #   set the plot colors
    # ---------------------------------------------------------------------
    # defining the colors
    myblack = [0, 0, 0]
    myblue = '#0F95D7'
    myred = '#e41a1c'
    myyellow = [255 / 255, 194 / 255, 10 / 255]
    mygreen = '#4daf4a'
    mybrown = '#a65628'
    mydarkblue = '#377eb8'
    mypurple = '#984ea3'
    myorange = '#ff7f00'
    mygray = [89 / 255, 89 / 255, 89 / 255]

    # setting up color and marker sequence
    # note there are only 10 but you should not have more than 10 lines on a plot
    colors = [myblack, myblue, myred, myyellow, mydarkblue, myorange, mybrown, mygreen, mypurple, mygray]
    markers = ['o', '^', 's', 'p', 'v', '*', 'x']

    # ---------------------------------------------------------------------
    #   load results
    # ---------------------------------------------------------------------
    blob_size                           = 100 /1                 # for solutions in plot
    x_label                             = 'Time (min)'
    y1_label                            = '$h_{air,cond}$ (W/m$^2$K)'
    y2_label                            = '$h_{liq,cond}$ (W/m$^2$K)'
    y3_label                            = '$h_{2ph,cond}$ (W/m$^2$K)'
    y4_label                            = '$h_{vap,cond}$ (W/m$^2$K)'
    y5_label                            = '$h_{2ph,evap}$ (W/m$^2$K)'
    y6_label                            = '$h_{vap,evap}$ (W/m$^2$K)'
    y7_label                            = '$\sigma_{air,cond}$ (%)'
    y8_label                            = '$\sigma_{liq,cond}$ (%)'
    y9_label                            = '$\sigma_{2ph,cond}$ (%)'
    y10_label                            = '$\sigma_{vap,cond}$ (%)'
    y11_label                            = '$\sigma_{2ph,evap}$ (%)'
    y12_label                            = '$\sigma_{vap,evap}$ (%)'
    # ---------------------------------------------------------------------
    #   plot results
    # ---------------------------------------------------------------------
    fig1 = plt.figure(None)
    fig1.set_size_inches(width, height)
    fig1.suptitle(None)

    plt.subplot(3, 2, 1)
    plt.scatter(BTMS_results['mission_time'][:-10], BTMS_results['err_HTC_air'][:-10], marker=markers[3], alpha=alpha, s=blob_size, c='none', edgecolors=colors[0], linewidths=line_width)
    plt.grid(visible=True, which='major', linestyle='-', linewidth=0.75, color=colors[0])
    plt.minorticks_on()
    plt.grid(visible=True, which='minor', linestyle='--', linewidth=0.15, color=colors[0])
    plt.ylabel(y7_label)

    plt.subplot(3, 2, 2)
    plt.scatter(BTMS_results['mission_time'][:-10], BTMS_results['err_HTC_liq_cond'][:-10], marker=markers[3], alpha=alpha, s=blob_size, c='none', edgecolors=colors[1], linewidths=line_width)
    plt.grid(visible=True, which='major', linestyle='-', linewidth=0.75, color=colors[0])
    plt.minorticks_on()
    plt.grid(visible=True, which='minor', linestyle='--', linewidth=0.15, color=colors[0])
    plt.ylabel(y8_label)

    plt.subplot(3, 2, 3)
    plt.scatter(BTMS_results['mission_time'][:-10], BTMS_results['err_HTC_2ph_cond'][:-10], marker=markers[3], alpha=alpha, s=blob_size, c='none', edgecolors=colors[2], linewidths=line_width)
    plt.grid(visible=True, which='major', linestyle='-', linewidth=0.75, color=colors[0])
    plt.minorticks_on()
    plt.grid(visible=True, which='minor', linestyle='--', linewidth=0.15, color=colors[0])
    plt.ylabel(y9_label)

    plt.subplot(3, 2, 4)
    plt.scatter(BTMS_results['mission_time'][:-10], BTMS_results['err_HTC_vap_cond'][:-10], marker=markers[3], alpha=alpha, s=blob_size, c='none', edgecolors=colors[3], linewidths=line_width)
    plt.grid(visible=True, which='major', linestyle='-', linewidth=0.75, color=colors[0])
    plt.minorticks_on()
    plt.grid(visible=True, which='minor', linestyle='--', linewidth=0.15, color=colors[0])
    plt.ylabel(y10_label)

    plt.subplot(3, 2, 5)
    plt.scatter(BTMS_results['mission_time'][:-10], BTMS_results['err_HTC_2ph_evap'][:-10], marker=markers[3], alpha=alpha, s=blob_size, c='none', edgecolors=colors[4], linewidths=line_width)
    plt.grid(visible=True, which='major', linestyle='-', linewidth=0.75, color=colors[0])
    plt.minorticks_on()
    plt.grid(visible=True, which='minor', linestyle='--', linewidth=0.15, color=colors[0])
    plt.xlabel(x_label)
    plt.ylabel(y11_label)

    plt.subplot(3, 2, 6)
    plt.scatter(BTMS_results['mission_time'][:-10], BTMS_results['err_HTC_vap_evap'][:-10], marker=markers[3], alpha=alpha, s=blob_size, c='none', edgecolors=colors[5], linewidths=line_width)
    plt.grid(visible=True, which='major', linestyle='-', linewidth=0.75, color=colors[0])
    plt.minorticks_on()
    plt.grid(visible=True, which='minor', linestyle='--', linewidth=0.15, color=colors[0])
    plt.xlabel(x_label)
    plt.ylabel(y12_label)

    plt.tight_layout()
    plt.show()
    # plt.savefig('fig10.png', dpi=image_resolution, format='png')


# ---------------------------------------------------------------------
#   Main
# ---------------------------------------------------------------------
def main_1():
    # Parameters of a generic HP to be used in the VC model
    # ------------------------------------------------------------------------------
    Comp        = {'Comp: rpm': 3.6e3*1,      'Comp: eta_C_iso': 0.63,    'Comp: Disp': 4.25e-5, 'Comp: eta_C_vol': 0.95, 'Comp: eta_C_mec': 0.88}
    # Evap        = {'Evap: TubeL': 10,       'Evap: TubeN': 5,           "Evap: TubeID": 0.011, 'Evap: FinEff': 0.75,    'Evap: FinRatio': 7.2, 'Evap: Pump': 220}
    Cond        = {'Cond: TubeL': 10,       'Cond: TubeN': 5,           "Cond: TubeID": 0.011, 'Cond: FinEff': 0.75,    'Cond: FinRatio': 7.2, 'Cond: Fan': 365}
    System      = {'Sys: Superheat': 8.5,   'Sys: Subcool': 6.5}
    SysParams   = dict(Comp, **Cond, **System)   # add the dicts together

    m_dot_air                       = 0.47 * 1
    # -----------------------------------------------------
    address                         = os.getcwd()
    with open(address + '\BatInputs_for_VCS.pkl', 'rb') as f:
        BatInputs                   = pickle.load(f)
    DT_amb                          = 35
    T_bat                           = BatInputs['T_b_air'][0] + DT_amb

    VCS_results = {'mission_time': [], 'T_bat_liq': [], 'T_bat_air': [], 'COP': [], 'Power': [], 'L_2ph': [], 'Q_bat': [],
                   'Q_conv': [], 'm_dot_ref': [], 'm_VCC': [], 'P1': [], 'T1_sat': [], 'P2': [], 'T2_sat': [],
                   'HTC_2ph_cond': [], 'HTC_liq_cond': [], 'HTC_vap_cond': [], 'HTC_2ph_evap': [], 'HTC_vap_evap': [], 'HTC_air': []}
    # ------------------------------------------------------
    for i in range(1, len(BatInputs['T_amb'])):
        T_air, Q_bat, dt            = BatteryData(i=i, DT_amb=DT_amb)
        T_bat_air                   = BatInputs['T_b_air'][i] + DT_amb
        HP        = Cycle_Performance_BatteryAir(SysParams, dt, Q_bat, T_bat, m_dot_air, T_air, 'air', 'R134a')
        HP.Cycle_Solver()
        HP.PostProcessing()
        # HP.Interpolation_Singularity_Points()
        # show the results
        print('********************i = {:.2f}**********************'.format(i))
        print('COP:                 {:.2f} [-]'.format(HP.COP))
        print('Q_cond:              {:.2f} [W]'.format(HP.Q_condenser))
        print('Q_evap:              {:.2f} [W]'.format(HP.Q_evaporator))
        print('W_comp:              {:.2f} [W]'.format(HP.W_c))
        print('T_bat_updated:       {:.2f} [K]'.format(HP.T_bat_updated))
        print('Energy Balance:      {:.2f} [-]'.format(HP.E_balance))

        print('-----------------------------------------------------------------------------------')
        print('         pressure [Pa],   enthalpy [J/kg],    temperature [K],    entropy [J/K],  mixture[-]')
        print('Point-1: {:.2f}       {:.2f}           {:.2f}              {:.2f}         {:.2f}'.format(HP.state_1[0], HP.state_1[1], HP.state_1[2], HP.state_1[3], HP.state_1[4]))
        print('Point-2: {:.2f}       {:.2f}           {:.2f}              {:.2f}         {:.2f}'.format(HP.state_2[0], HP.state_2[1], HP.state_2[2], HP.state_2[3], HP.state_2[4]))
        print('Point-3: {:.2f}       {:.2f}           {:.2f}              {:.2f}         {:.2f}'.format(HP.state_3[0], HP.state_3[1], HP.state_3[2], HP.state_3[3], HP.state_3[4]))
        print('Point-4: {:.2f}       {:.2f}           {:.2f}              {:.2f}         {:.2f}'.format(HP.state_4[0], HP.state_4[1], HP.state_4[2], HP.state_4[3], HP.state_4[4]))
        # ----------
        # 2-> 3
        h_vap_cond, err_h_vap_cond = HTC_cond_vap_comparison(HP.state_3[2] + HP.subcooling, HP.state_2[2], HP.m_dot_refrigerant)
        h_2ph_cond, err_h_2ph_cond = HTC_cond_2ph_comparison(HP.state_3[2] + HP.subcooling, HP.m_dot_refrigerant)
        h_liq_cond, err_h_liq_cond = HTC_cond_liq_comparison(HP.state_3[2] + HP.subcooling, HP.state_3[2], HP.m_dot_refrigerant)
        # 4->1
        h_2ph_evap, err_h_2ph_evap = HTC_evap_2ph_comparison(HP.state_1[2] - HP.superheating, HP.m_dot_refrigerant, HP.state_4[4], HP.evap_L_1)
        h_vap_evap, err_h_vap_evap = HTC_evap_vap_comparison(HP.state_1[2] - HP.superheating, HP.state_1[2], HP.m_dot_refrigerant)
        # ----------
        T_bat   = HP.T_bat_updated
        VCS_results['T_bat_liq'].append(T_bat - 273.15)
        VCS_results['T_bat_air'].append(T_bat_air - 273.15)
        VCS_results['L_2ph'].append(HP.evap_L_1 / HP.evap_L_tot * 100)
        time    = BatInputs['time'][i]
        VCS_results['mission_time'].append(time)
        VCS_results['COP'].append(HP.COP)
        VCS_results['Power'].append(HP.W_c/1e3)
        VCS_results['Q_bat'].append(HP.Q_bat/1e3)
        VCS_results['Q_conv'].append(HP.Q_evaporator/1e3)
        VCS_results['m_dot_ref'].append(HP.m_dot_refrigerant)
        m_VCC       = VCC_weight_estimation(HP.W_c/1e3, HP.m_dot_refrigerant)
        VCS_results['m_VCC'].append(m_VCC)

        VCS_results['P1'].append(HP.state_1[0]/1e3)
        VCS_results['T1_sat'].append(HP.state_1[2] - HP.superheating)
        VCS_results['P2'].append(HP.state_3[0]/1e3)
        VCS_results['T2_sat'].append(HP.state_3[2] + HP.subcooling)

    # plot_BTMS_results(VCS_results)
    plot_BTMS_results2(VCS_results)

    # # ----------------------------------------------------------
    # T_bat, T_air, Q_bat, dt         = BatteryData(i=97)
    #
    # HP        = Cycle_Performance_BatteryAir(SysParams, dt, Q_bat, T_bat, m_dot_air, T_air, 'air', 'R134a')
    # HP.Cycle_Solver()
    # HP.PostProcessing()
    # # HP.Interpolation_Singularity_Points()
    # # show the results
    # print('i = :{:.2f}'.format(97))
    # print('COP:         {:.2f}'.format(HP.COP))
    # print('Q_cond:      {:.2f} [W]'.format(HP.Q_condenser))
    # print('Q_evap:      {:.2f} [W]'.format(HP.Q_evaporator))
    # print('W_comp:      {:.2f} [W]'.format(HP.W_c))
    # print('T_bat_updated: {:.2f} [K]'.format(HP.T_bat_updated))
    # print('Energy Balance: {:.2f} [-]'.format(HP.E_balance))
    #
    # print('-----------------------------------------------------------------------------------')
    # print('         pressure [Pa],   enthalpy [J/kg],    temperature [K],    entropy [J/K],  mixture[-]')
    # print('Point-1: {:.2f}       {:.2f}           {:.2f}              {:.2f}         {:.2f}'.format(HP.state_1[0], HP.state_1[1], HP.state_1[2], HP.state_1[3], HP.state_1[4]))
    # print('Point-2: {:.2f}       {:.2f}           {:.2f}              {:.2f}         {:.2f}'.format(HP.state_2[0], HP.state_2[1], HP.state_2[2], HP.state_2[3], HP.state_2[4]))
    # print('Point-3: {:.2f}       {:.2f}           {:.2f}              {:.2f}         {:.2f}'.format(HP.state_3[0], HP.state_3[1], HP.state_3[2], HP.state_3[3], HP.state_3[4]))
    # print('Point-4: {:.2f}       {:.2f}           {:.2f}              {:.2f}         {:.2f}'.format(HP.state_4[0], HP.state_4[1], HP.state_4[2], HP.state_4[3], HP.state_4[4])


def main_2():
    # Parameters of a generic HP to be used in the VC model
    # ------------------------------------------------------------------------------
    Comp        = {'Comp: rpm': 3.6e3*1,      'Comp: eta_C_iso': 0.63,    'Comp: Disp': 4.25e-5, 'Comp: eta_C_vol': 0.95, 'Comp: eta_C_mec': 0.88}
    # Evap        = {'Evap: TubeL': 10,       'Evap: TubeN': 5,           "Evap: TubeID": 0.011, 'Evap: FinEff': 0.75,    'Evap: FinRatio': 7.2, 'Evap: Pump': 220}
    Cond        = {'Cond: TubeL': 10,       'Cond: TubeN': 5,           "Cond: TubeID": 0.011, 'Cond: FinEff': 0.75,    'Cond: FinRatio': 7.2, 'Cond: Fan': 365}
    System      = {'Sys: Superheat': 8.5,   'Sys: Subcool': 6.5}
    SysParams   = dict(Comp, **Cond, **System)   # add the dicts together

    m_dot_air                       = 0.23 * 1
    # -----------------------------------------------------
    address                         = os.getcwd()
    with open(address + '\BatInputs_for_VCS.pkl', 'rb') as f:
        BatInputs                   = pickle.load(f)
    DT_amb                          = 35
    T_bat                           = BatInputs['T_b_air'][0] + DT_amb

    VCS_results = {'mission_time': [], 'T_bat_liq': [], 'T_bat_air': [], 'COP': [], 'Power': [], 'L_2ph': [], 'Q_bat': [],
                   'Q_conv': [], 'm_dot_ref': [], 'm_VCC': [], 'P1': [], 'T1_sat': [], 'P2': [], 'T2_sat': [],
                   'HTC_2ph_cond': [], 'HTC_liq_cond': [], 'HTC_vap_cond': [], 'HTC_2ph_evap': [], 'HTC_vap_evap': [], 'HTC_air': []}
    # ------------------------------------------------------
    for i in range(1, len(BatInputs['T_amb'])):
        T_air, Q_bat, dt            = BatteryData(i=i, DT_amb=DT_amb)
        T_bat_air                   = BatInputs['T_b_air'][i] + DT_amb
        HP        = Cycle_Performance_BatteryAir(SysParams, dt, Q_bat, T_bat, m_dot_air, T_air, 'air', 'R134a')
        HP.Cycle_Solver()
        HP.PostProcessing()
        # HP.Interpolation_Singularity_Points()
        # show the results
        print('********************i = {:.2f}**********************'.format(i))
        print('COP:                 {:.2f} [-]'.format(HP.COP))
        print('Q_cond:              {:.2f} [W]'.format(HP.Q_condenser))
        print('Q_evap:              {:.2f} [W]'.format(HP.Q_evaporator))
        print('W_comp:              {:.2f} [W]'.format(HP.W_c))
        print('T_bat_updated:       {:.2f} [K]'.format(HP.T_bat_updated))
        print('Energy Balance:      {:.2f} [-]'.format(HP.E_balance))

        print('-----------------------------------------------------------------------------------')
        print('         pressure [Pa],   enthalpy [J/kg],    temperature [K],    entropy [J/K],  mixture[-]')
        print('Point-1: {:.2f}       {:.2f}           {:.2f}              {:.2f}         {:.2f}'.format(HP.state_1[0], HP.state_1[1], HP.state_1[2], HP.state_1[3], HP.state_1[4]))
        print('Point-2: {:.2f}       {:.2f}           {:.2f}              {:.2f}         {:.2f}'.format(HP.state_2[0], HP.state_2[1], HP.state_2[2], HP.state_2[3], HP.state_2[4]))
        print('Point-3: {:.2f}       {:.2f}           {:.2f}              {:.2f}         {:.2f}'.format(HP.state_3[0], HP.state_3[1], HP.state_3[2], HP.state_3[3], HP.state_3[4]))
        print('Point-4: {:.2f}       {:.2f}           {:.2f}              {:.2f}         {:.2f}'.format(HP.state_4[0], HP.state_4[1], HP.state_4[2], HP.state_4[3], HP.state_4[4]))
        # ----------
        # 2-> 3
        h_vap_cond, err_h_vap_cond = HTC_cond_vap_comparison(HP.state_3[2] + HP.subcooling, HP.state_2[2], HP.m_dot_refrigerant)
        h_2ph_cond, err_h_2ph_cond = HTC_cond_2ph_comparison(HP.state_3[2] + HP.subcooling, HP.m_dot_refrigerant)
        h_liq_cond, err_h_liq_cond = HTC_cond_liq_comparison(HP.state_3[2] + HP.subcooling, HP.state_3[2], HP.m_dot_refrigerant)
        # 4->1
        h_2ph_evap, err_h_2ph_evap = HTC_evap_2ph_comparison(HP.state_1[2] - HP.superheating, HP.m_dot_refrigerant, HP.state_4[4], HP.evap_L_1)
        h_vap_evap, err_h_vap_evap = HTC_evap_vap_comparison(HP.state_1[2] - HP.superheating, HP.state_1[2], HP.m_dot_refrigerant)
        # ----------
        T_bat   = HP.T_bat_updated
        VCS_results['T_bat_liq'].append(T_bat - 273.15)
        VCS_results['T_bat_air'].append(T_bat_air - 273.15)
        VCS_results['L_2ph'].append(HP.evap_L_1 / HP.evap_L_tot * 100)
        time    = BatInputs['time'][i]
        VCS_results['mission_time'].append(time)
        VCS_results['COP'].append(HP.COP)
        VCS_results['Power'].append(HP.W_c/1e3)
        VCS_results['Q_bat'].append(HP.Q_bat/1e3)
        VCS_results['Q_conv'].append(HP.Q_evaporator/1e3)
        VCS_results['m_dot_ref'].append(HP.m_dot_refrigerant)
        m_VCC       = VCC_weight_estimation(HP.W_c/1e3, HP.m_dot_refrigerant)
        VCS_results['m_VCC'].append(m_VCC)

        VCS_results['P1'].append(HP.state_1[0]/1e3)
        VCS_results['T1_sat'].append(HP.state_1[2] - HP.superheating)
        VCS_results['P2'].append(HP.state_3[0]/1e3)
        VCS_results['T2_sat'].append(HP.state_3[2] + HP.subcooling)

    # plot_BTMS_results(VCS_results)
    plot_BTMS_results2(VCS_results)


def main_herringbone001():
    # Parameters of a generic HP to be used in the VC model
    # ------------------------------------------------------------------------------
    Comp        = {'Comp: rpm': 3.6e3*1,      'Comp: eta_C_iso': 0.63,    'Comp: Disp': 4.25e-5, 'Comp: eta_C_vol': 0.95, 'Comp: eta_C_mec': 0.88}
    # Evap        = {'Evap: TubeL': 10,       'Evap: TubeN': 5,           "Evap: TubeID": 0.011, 'Evap: FinEff': 0.75,    'Evap: FinRatio': 7.2, 'Evap: Pump': 220}
    Cond        = {'Cond: TubeL': 10,       'Cond: TubeN': 5,           "Cond: TubeID": 0.011, 'Cond: FinEff': 0.75,    'Cond: FinRatio': 7.2, 'Cond: Fan': 365}
    System      = {'Sys: Superheat': 8.5,   'Sys: Subcool': 6.5}
    SysParams   = dict(Comp, **Cond, **System)   # add the dicts together

    V_dot_air                       = 1.7934 * 0.3
    # -----------------------------------------------------
    address                         = os.getcwd()
    with open(address + '\BatInputs_for_VCS.pkl', 'rb') as f:
        BatInputs                   = pickle.load(f)
    DT_amb                          = 35
    T_bat                           = BatInputs['T_b_air'][0] + DT_amb

    VCS_results = {'mission_time': [], 'T_bat_liq': [], 'T_bat_air': [], 'COP': [], 'Power': [], 'L_2ph': [], 'Q_bat': [],
                   'Q_conv': [], 'm_dot_ref': [], 'm_VCC': [], 'P1': [], 'T1_sat': [], 'P2': [], 'T2_sat': [],
                   'HTC_2ph_cond': [], 'HTC_liq_cond': [], 'HTC_vap_cond': [], 'HTC_2ph_evap': [], 'HTC_vap_evap': [], 'HTC_air': []}
    # ------------------------------------------------------
    for i in range(1, len(BatInputs['T_amb'])):
        T_air, Q_bat, dt            = BatteryData(i=i, DT_amb=DT_amb)
        T_bat_air                   = BatInputs['T_b_air'][i] + DT_amb
        # ------------------------------------------------------
        m_dot_air_cal, A_fin, eta_fin, HTC_air_cal = HerringboneFinCal(V_dot_air, T_air)
        print(m_dot_air_cal, A_fin, eta_fin, HTC_air_cal)
        m_dot_air   = 0.47    # todo check
        eta_fin     = 0.6
        # ------------------------------------------------------
        HP        = Cycle_Performance_BatteryAir_herringbone(SysParams, dt, Q_bat, T_bat, m_dot_air, T_air, eta_fin, 'air', 'R134a')
        HP.Cycle_Solver()
        HP.PostProcessing()
        # HP.Interpolation_Singularity_Points()
        # show the results
        print('********************i = {:.2f}**********************'.format(i))
        print('COP:                 {:.2f} [-]'.format(HP.COP))
        print('Q_cond:              {:.2f} [W]'.format(HP.Q_condenser))
        print('Q_evap:              {:.2f} [W]'.format(HP.Q_evaporator))
        print('W_comp:              {:.2f} [W]'.format(HP.W_c))
        print('T_bat_updated:       {:.2f} [K]'.format(HP.T_bat_updated))
        print('Energy Balance:      {:.2f} [-]'.format(HP.E_balance))

        print('-----------------------------------------------------------------------------------')
        print('         pressure [Pa],   enthalpy [J/kg],    temperature [K],    entropy [J/K],  mixture[-]')
        print('Point-1: {:.2f}       {:.2f}           {:.2f}              {:.2f}         {:.2f}'.format(HP.state_1[0], HP.state_1[1], HP.state_1[2], HP.state_1[3], HP.state_1[4]))
        print('Point-2: {:.2f}       {:.2f}           {:.2f}              {:.2f}         {:.2f}'.format(HP.state_2[0], HP.state_2[1], HP.state_2[2], HP.state_2[3], HP.state_2[4]))
        print('Point-3: {:.2f}       {:.2f}           {:.2f}              {:.2f}         {:.2f}'.format(HP.state_3[0], HP.state_3[1], HP.state_3[2], HP.state_3[3], HP.state_3[4]))
        print('Point-4: {:.2f}       {:.2f}           {:.2f}              {:.2f}         {:.2f}'.format(HP.state_4[0], HP.state_4[1], HP.state_4[2], HP.state_4[3], HP.state_4[4]))
        # ----------
        # 2-> 3
        h_vap_cond, err_h_vap_cond = HTC_cond_vap_comparison(HP.state_3[2] + HP.subcooling, HP.state_2[2], HP.m_dot_refrigerant)
        h_2ph_cond, err_h_2ph_cond = HTC_cond_2ph_comparison(HP.state_3[2] + HP.subcooling, HP.m_dot_refrigerant)
        h_liq_cond, err_h_liq_cond = HTC_cond_liq_comparison(HP.state_3[2] + HP.subcooling, HP.state_3[2], HP.m_dot_refrigerant)
        # 4->1
        h_2ph_evap, err_h_2ph_evap = HTC_evap_2ph_comparison(HP.state_1[2] - HP.superheating, HP.m_dot_refrigerant, HP.state_4[4], HP.evap_L_1)
        h_vap_evap, err_h_vap_evap = HTC_evap_vap_comparison(HP.state_1[2] - HP.superheating, HP.state_1[2], HP.m_dot_refrigerant)
        # ----------
        T_bat   = HP.T_bat_updated
        VCS_results['T_bat_liq'].append(T_bat - 273.15)
        VCS_results['T_bat_air'].append(T_bat_air - 273.15)
        VCS_results['L_2ph'].append(HP.evap_L_1 / HP.evap_L_tot * 100)
        time    = BatInputs['time'][i]
        VCS_results['mission_time'].append(time)
        VCS_results['COP'].append(HP.COP)
        VCS_results['Power'].append(HP.W_c/1e3)
        VCS_results['Q_bat'].append(HP.Q_bat/1e3)
        VCS_results['Q_conv'].append(HP.Q_evaporator/1e3)
        VCS_results['m_dot_ref'].append(HP.m_dot_refrigerant)
        m_VCC       = VCC_weight_estimation(HP.W_c/1e3, HP.m_dot_refrigerant)
        VCS_results['m_VCC'].append(m_VCC)

        VCS_results['P1'].append(HP.state_1[0]/1e3)
        VCS_results['T1_sat'].append(HP.state_1[2] - HP.superheating)
        VCS_results['P2'].append(HP.state_3[0]/1e3)
        VCS_results['T2_sat'].append(HP.state_3[2] + HP.subcooling)

    # plot_BTMS_results(VCS_results)
    plot_BTMS_results2(VCS_results)


def main_herringbone():
    # Parameters of a generic HP to be used in the VC model
    # ------------------------------------------------------------------------------
    Comp        = {'Comp: rpm': 3.6e3*2,      'Comp: eta_C_iso': 0.63,    'Comp: Disp': 4.25e-5*1, 'Comp: eta_C_vol': 0.95, 'Comp: eta_C_mec': 0.88}
    # Evap        = {'Evap: TubeL': 10,       'Evap: TubeN': 5,           "Evap: TubeID": 0.011, 'Evap: FinEff': 0.75,    'Evap: FinRatio': 7.2, 'Evap: Pump': 220}
    Cond        = {'Cond: TubeL': 10,       'Cond: TubeN': 5,           "Cond: TubeID": 0.011, 'Cond: FinEff': 0.75,    'Cond: FinRatio': 7.2, 'Cond: Fan': 365}
    System      = {'Sys: Superheat': 8.5,   'Sys: Subcool': 6.5}
    SysParams   = dict(Comp, **Cond, **System)   # add the dicts together

    # fin structure inputs
    V_dot_ha                        = 1.7934 * 2

    # -----------------------------------------------------
    address                         = os.getcwd()
    with open(address + '\BatInputs_for_VCS.pkl', 'rb') as f:
        BatInputs                   = pickle.load(f)
    DT_amb                          = 35
    T_bat                           = BatInputs['T_b_air'][0] + DT_amb

    VCS_results = {'mission_time': [], 'T_bat_liq': [], 'T_bat_air': [], 'COP': [], 'Power': [], 'L_2ph': [], 'Q_bat': [],
                   'Q_conv': [], 'm_dot_ref': [], 'm_VCC': [], 'P1': [], 'T1_sat': [], 'P2': [], 'T2_sat': [],
                   'HTC_2ph_cond': [], 'HTC_liq_cond': [], 'HTC_vap_cond': [], 'HTC_2ph_evap': [], 'HTC_vap_evap': [], 'HTC_air': [],
                   'err_HTC_2ph_cond': [], 'err_HTC_liq_cond': [], 'err_HTC_vap_cond': [], 'err_HTC_2ph_evap': [], 'err_HTC_vap_evap': [], 'err_HTC_air': []}
    # ------------------------------------------------------
    for i in range(1, len(BatInputs['T_amb'])):
        T_air, Q_bat, dt            = BatteryData(i=i, DT_amb=DT_amb)
        T_bat_air                   = BatInputs['T_b_air'][i] + DT_amb

        m_dot_air_cal, A_fin, eta_fin, HTC_air_cal, fin_ratio = HerringboneFinCal(V_dot_ha, T_air, p=101325, RH=0.51)
        m_dot_air   = m_dot_air_cal
        eta_fin     = eta_fin   # overall efficiency
        fin_ratio   = fin_ratio # outside finned area (total wet: A_ubfinned + A_finned) / inside tube wet area (A_i)
        HP        = Cycle_Performance_BatteryAir_herringbone(SysParams, dt, Q_bat, T_bat, m_dot_air, T_air, eta_fin, fin_ratio, 'air', 'R134a')
        HP.Cycle_Solver()
        # HP.PostProcessing()
        HP.Interpolation_Singularity_Points()
        # HP.Interpolation_Singularity_Points()
        # show the results
        print('********************i = {:.2f}**********************'.format(i))
        # print('COP:                 {:.2f} [-]'.format(HP.COP))
        # print('Q_cond:              {:.2f} [W]'.format(HP.Q_condenser))
        # print('Q_evap:              {:.2f} [W]'.format(HP.Q_evaporator))
        # print('W_comp:              {:.2f} [W]'.format(HP.W_c))
        # print('T_bat_updated:       {:.2f} [K]'.format(HP.T_bat_updated))
        print('Energy Balance:      {:.2f} [-]'.format(HP.E_balance))

        print('-----------------------------------------------------------------------------------')
        print('         pressure [Pa],   enthalpy [J/kg],    temperature [K],    entropy [J/K],  mixture[-]')
        print('Point-1: {:.2f}       {:.2f}           {:.2f}              {:.2f}         {:.2f}'.format(HP.state_1[0], HP.state_1[1], HP.state_1[2], HP.state_1[3], HP.state_1[4]))
        print('Point-2: {:.2f}       {:.2f}           {:.2f}              {:.2f}         {:.2f}'.format(HP.state_2[0], HP.state_2[1], HP.state_2[2], HP.state_2[3], HP.state_2[4]))
        print('Point-3: {:.2f}       {:.2f}           {:.2f}              {:.2f}         {:.2f}'.format(HP.state_3[0], HP.state_3[1], HP.state_3[2], HP.state_3[3], HP.state_3[4]))
        print('Point-4: {:.2f}       {:.2f}           {:.2f}              {:.2f}         {:.2f}'.format(HP.state_4[0], HP.state_4[1], HP.state_4[2], HP.state_4[3], HP.state_4[4]))
        # ----------
        # 2-> 3
        h_vap_cond, err_h_vap_cond = HTC_cond_vap_comparison(HP.state_3[2] + HP.subcooling, HP.state_2[2], HP.m_dot_refrigerant)
        h_2ph_cond, err_h_2ph_cond = HTC_cond_2ph_comparison(HP.state_3[2] + HP.subcooling, HP.m_dot_refrigerant)
        h_liq_cond, err_h_liq_cond = HTC_cond_liq_comparison(HP.state_3[2] + HP.subcooling, HP.state_3[2], HP.m_dot_refrigerant)
        # 4->1
        h_2ph_evap, err_h_2ph_evap = HTC_evap_2ph_comparison(HP.state_1[2] - HP.superheating, HP.m_dot_refrigerant, HP.state_4[4], HP.evap_L_1)
        h_vap_evap, err_h_vap_evap = HTC_evap_vap_comparison(HP.state_1[2] - HP.superheating, HP.state_1[2], HP.m_dot_refrigerant)
        # condenser, airside
        h_air_cond, err_h_air_cond = HTC_cond_airside_finned_tube(HTC_air_cal)

        print('h_2ph_evap: ', h_2ph_evap, 'h_vap_evap: ', h_vap_evap)
        # ----------
        T_bat   = HP.T_bat_updated
        VCS_results['T_bat_liq'].append(T_bat - 273.15)
        VCS_results['T_bat_air'].append(T_bat_air - 273.15)
        VCS_results['L_2ph'].append(HP.evap_L_1 / HP.evap_L_tot * 100)
        time    = BatInputs['time'][i]
        VCS_results['mission_time'].append(time)
        VCS_results['COP'].append(HP.COP)
        VCS_results['Power'].append(HP.W_c/1e3)
        VCS_results['Q_bat'].append(HP.Q_bat/1e3)
        VCS_results['Q_conv'].append(HP.Q_evaporator/1e3)
        VCS_results['m_dot_ref'].append(HP.m_dot_refrigerant)
        m_VCC       = VCC_weight_estimation(HP.W_c/1e3, HP.m_dot_refrigerant)
        VCS_results['m_VCC'].append(m_VCC)

        VCS_results['P1'].append(HP.state_1[0]/1e3)
        VCS_results['T1_sat'].append(HP.state_1[2] - HP.superheating)
        VCS_results['P2'].append(HP.state_3[0]/1e3)
        VCS_results['T2_sat'].append(HP.state_3[2] + HP.subcooling)

        # -----29/07/2024 ----------------------------------------------------
        VCS_results['HTC_air'].append(h_air_cond)
        VCS_results['HTC_liq_cond'].append(h_liq_cond)
        VCS_results['HTC_2ph_cond'].append(h_2ph_cond)
        VCS_results['HTC_vap_cond'].append(h_vap_cond)
        VCS_results['HTC_2ph_evap'].append(h_2ph_evap)
        VCS_results['HTC_vap_evap'].append(h_vap_evap)

        VCS_results['err_HTC_air'].append(err_h_air_cond)
        VCS_results['err_HTC_liq_cond'].append(err_h_liq_cond)
        VCS_results['err_HTC_2ph_cond'].append(err_h_2ph_cond)
        VCS_results['err_HTC_vap_cond'].append(err_h_vap_cond)
        VCS_results['err_HTC_2ph_evap'].append(err_h_2ph_evap)
        VCS_results['err_HTC_vap_evap'].append(err_h_vap_evap)
    plot_BTMS_results(VCS_results)
    # plot_BTMS_results2(VCS_results)
    plot_BTMS_results3(VCS_results)


# ---------------------------------------------------------------------
if __name__ == '__main__':

    # V_dot_air                       = 1.7934 * 0.25
    # T_air                           = 15 + 273.15 + 35
    # HerringboneFinCal(V_dot_air, T_air)

    main_herringbone()
    # main_2()
    # main_1()
