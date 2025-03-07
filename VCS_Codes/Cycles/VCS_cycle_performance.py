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

from VCS_Codes.Components.VCS_components         import VCS_condenser_Finned_Tube_3Zones, VCS_evaporator_Finned_Tube_2Zones


# ---------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
class Cycle_Performance_Air:
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

    def __init__(self, SysParams, m_dot_coolant, T_coolant_inlet, m_dot_air, T_air_inlet, coolant, air, refrigerant):
        self.rpm                        = SysParams['Comp: rpm']
        self.eta_C_iso                  = SysParams['Comp: eta_C_iso']
        self.Displacement               = SysParams['Comp: Disp']
        self.eta_C_vol                  = SysParams['Comp: eta_C_vol']
        self.eta_C_mec                  = SysParams['Comp: eta_C_mec']

        self.superheating               = SysParams['Sys: Superheat']
        self.subcooling                 = SysParams['Sys: Subcool']

        self.W_Fan_COND                 = SysParams['Cond: Fan']
        self.W_Pump_EVAP                = SysParams['Evap: Pump']

        self.m_dot_coolant              = m_dot_coolant
        self.T_coolant_inlet            = T_coolant_inlet
        self.m_dot_air                  = m_dot_air
        self.T_air_inlet                = T_air_inlet
        self.m_dot_refrigerant          = None
        self.refrigerant                = refrigerant
        self.heat_transfer_fluid_evap   = coolant
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
            self.C_min_evap             = CoolProp.PropsSI('C', 'P', 101325, 'T', self.T_coolant_inlet, self.heat_transfer_fluid_evap) * self.m_dot_coolant
            self.C_min_cond             = CoolProp.PropsSI('C', 'P', 101325, 'T', self.T_air_inlet, self.heat_transfer_fluid_cond) * self.m_dot_air

            # ---------------------------------------------------------
            # Point 1
            # ---------------------------------------------------------
            self.T_1                    = round(CoolProp.PropsSI('T', 'P', P_guess[0], 'Q', 1, self.refrigerant) + self.superheating, self.dec) # old
            # self.T_1                    = round(CoolProp.PropsSI('T', 'P', P_guess[0], 'Q', 1, self.refrigerant), self.dec) # todo 22/07/2024
            self.h_1                    = CoolProp.PropsSI('H', 'P', P_guess[0], 'T', self.T_1, self.refrigerant)
            self.x_1                    = CoolProp.PropsSI('Q', 'P', P_guess[0], 'T', self.T_1, self.refrigerant)  # Q is the Quality
            self.s_1                    = CoolProp.PropsSI('S', 'P', P_guess[0], 'T', self.T_1, self.refrigerant)
            rho_1                       = CoolProp.PropsSI('D', 'P', P_guess[0], 'T', self.T_1, self.refrigerant)

            self.m_dot_refrigerant      = rho_1 * self.Displacement * (self.rpm / 60) * self.eta_C_vol
            self.state_1                = [P_guess[0], self.T_1, self.h_1, self.x_1]

            # ---------------------------------------------------------
            # Point 2
            # ---------------------------------------------------------
            self.s_2                    = self.s_1
            self.h_2s                   = CoolProp.PropsSI('H', 'P', P_guess[1], 'S', self.s_2, self.refrigerant)
            self.h_2                    = (self.h_2s - self.h_1) / self.eta_C_iso + self.h_1
            self.T_2                    = round(CoolProp.PropsSI('T', 'P', P_guess[1], 'H', self.h_2, self.refrigerant), self.dec)
            self.x_2                    = CoolProp.PropsSI('Q', 'P', P_guess[1], 'T', self.T_2, self.refrigerant)
            self.state_2                = [P_guess[1], self.T_2, self.h_2, self.x_2]
            # ---------------------------------------------------------
            # Point 3
            # ---------------------------------------------------------
            self.T_3                    = round(CoolProp.PropsSI('T', 'P', P_guess[1], 'Q', 0, self.refrigerant) - self.subcooling, self.dec)  # 8.4
            self.h_3                    = CoolProp.PropsSI('H', 'P', P_guess[1], 'T', self.T_3, self.refrigerant)

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
            L_tot                       = condenser_results.L_tot
            L_1                         = condenser_results.L_1
            L_2                         = condenser_results.L_2
            L_3                         = condenser_results.L_3
            # ---------------------------------------------------------------
            h_3_cal                     = Q_cond / self.m_dot_refrigerant + self.h_2
            T_3_cal                     = T_3_cal

            self.T_air_outlet           = T_air_out
            self.x_3                    = x_outlet
            self.state_3                = [P_guess[1], self.T_3, self.h_3, self.x_3]
            # ---------------------------------------------------------
            # Point 4
            self.h_4                    = self.h_3
            self.T_4                    = round(CoolProp.PropsSI('T', 'P', P_guess[0], 'H', self.h_4, self.refrigerant), self.dec)
            self.x_4                    = CoolProp.PropsSI('Q', 'P', P_guess[0], 'H', self.h_4, self.refrigerant)

            # Get outlet results of condenser through 2-zone model, assuming that the inlet refrigerant is in 2phase state
            evaporator                   = VCS_evaporator_Finned_Tube_2Zones(self.refrigerant, self.heat_transfer_fluid_evap, self.T_coolant_inlet,
                                          self.m_dot_refrigerant, self.m_dot_coolant, round(P_guess[0], 2), self.h_4, self.superheating, self.EV_HTC_coolant)
            evaporator_results           = evaporator.thermodynamics_calculation()
            # ---------------------------------------------------------------
            Q_evap                      = evaporator_results.Q_evaporator
            T_1_cal                     = evaporator_results.T_outlet
            P_outlet                    = evaporator_results.P_outlet
            h_outlet                    = evaporator_results.h_outlet
            x_inlet                     = evaporator_results.x_inlet
            x_outlet                    = evaporator_results.x_outlet
            T_coolant_out               = evaporator_results.T_coolant_out
            UA_tot                      = evaporator_results.UA_tot
            effectiveness               = evaporator_results.effectiveness
            L_tot                       = evaporator_results.L_tot
            L_1                         = evaporator_results.L_1
            L_2                         = evaporator_results.L_2
            # ---------------------------------------------------------------
            h_1_cal                     = Q_evap / self.m_dot_refrigerant + self.h_4
            T_1_cal                     = T_1_cal

            self.T_coolant_outlet       = T_coolant_out
            self.state_4                = [P_guess[0], self.T_4, self.h_4, self.x_4]
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
        P_1_initial                     = int(CoolProp.PropsSI('P', 'Q', 1, 'T', self.T_coolant_inlet - self.subcooling, self.refrigerant))
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
        self.COP                         = self.Q_evaporator / ((self.W_c / self.eta_C_mec) + self.W_Fan_COND + self.W_Pump_EVAP)

        self.E_balance                   = self.Q_condenser + self.Q_evaporator + self.W_c
        # 1W = 3.41 Btu/h
        self.Capacity                    = - self.Q_evaporator * 3.41

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
            P_1_initial                 = int(CoolProp.PropsSI('P', 'Q', 1, 'T', self.T_coolant_inlet - self.subcooling, self.refrigerant))
            P_2_initial                 = int(CoolProp.PropsSI('P', 'Q', 1, 'T', self.T_air_inlet - dif + self.superheating, self.refrigerant))

            P_guess                     = np.array([P_1_initial, P_2_initial])

            result2                     = optimize.root(self.March_Compressor_Inlet_and_Outlet_Pressures, P_guess, method='hybr')  # Minimization function
            self.P_1                    = result2.x[0]
            self.P_2                    = result2.x[1]
            res2                        = self.PostProcessing()

            # --------------------Simulation with T_amb + 1 K-------------------------------------
            P_1_initial                 = int(CoolProp.PropsSI('P', 'Q', 1, 'T', self.T_coolant_inlet - self.subcooling, self.refrigerant))
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


# ---------------------------------------------------------------------
#   Main
# ---------------------------------------------------------------------

if __name__ == '__main__':
    # Parameters of a generic HP to be used in the VC model
    # ------------------------------------------------------------------------------
    Comp        = {'Comp: rpm': 3.6e3,      'Comp: eta_C_iso': 0.63,    'Comp: Disp': 4.25e-5, 'Comp: eta_C_vol': 0.95, 'Comp: eta_C_mec': 0.88}
    Evap        = {'Evap: TubeL': 10,       'Evap: TubeN': 5,           "Evap: TubeID": 0.011, 'Evap: FinEff': 0.75,    'Evap: FinRatio': 7.2, 'Evap: Pump': 220}
    Cond        = {'Cond: TubeL': 10,       'Cond: TubeN': 5,           "Cond: TubeOD": 0.011, 'Cond: FinEff': 0.75,    'Cond: FinRatio': 7.2, 'Cond: Fan': 365}
    System      = {'Sys: Superheat': 8.5,   'Sys: Subcool': 6.5}
    SysParams   = dict(Comp, **Evap, **Cond, **System)   # add the dicts together

    m_dot_coolant                   = 0.1
    m_dot_air                       = 0.47
    T_coolant_inlet                 = 273 + 30
    T_air_inlet                     = 273 + 50

    HP        = Cycle_Performance_Air(SysParams, m_dot_coolant, T_coolant_inlet, m_dot_air, T_air_inlet, 'water', 'air', 'R134a')
    HP.Cycle_Solver()
    HP.PostProcessing()
    # HP.Interpolation_Singularity_Points()
    # show the results
    print('COP:         {:.2f}'.format(HP.COP))
    print('Q_cond:      {:.2f} [W]'.format(HP.Q_condenser))
    print('Q_evap:      {:.2f} [W]'.format(HP.Q_evaporator))
    print('W_comp:      {:.2f} [W]'.format(HP.W_c))
    print('T_coolant_outlet: {:.2f} [K]'.format(HP.T_coolant_outlet))
    print('Energy Balance: {:.2f} [-]'.format(HP.E_balance))
    print('P1:          {:.2f} [kPa]'.format(HP.P_1/1e3))
    print('T1_SAT:          {:.2f} [K]'.format(HP.T_1 - HP.superheating))
    print('P2:          {:.2f} [kPa]'.format(HP.P_2/1e3))
    print('T2_SAT:          {:.2f} [K]'.format(HP.T_3 + HP.subcooling))
