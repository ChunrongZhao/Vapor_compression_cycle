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

from VCS_Codes.Components.VCS_components         import VCS_condenser_Finned_Tube_3Zones, Battery_Wavychannel_Evaporator_2Zones
from ACHP_codes.Components.WavyChan import WavyChan_NMC
from matplotlib import rcParams, rc, font_manager
import matplotlib as mpl
import matplotlib.pyplot as plt
import os


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
            L_tot                       = evaporator_results.L_tot
            L_1                         = evaporator_results.L_1
            L_2                         = evaporator_results.L_2
            # ---------------------------------------------------------------
            h_1_cal                     = Q_evap / self.m_dot_refrigerant + self.h_4
            T_1_cal                     = T_1_cal
            self.T_bat_updated          = evaporator_results.T_bat_updated
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
#   Functions
# ---------------------------------------------------------------------
def BatteryData(i=1):
    address                         = os.getcwd()
    with open(address + '\BatInputs_for_VCS.pkl', 'rb') as f:
        BatInputs                   = pickle.load(f)

    # print(BatInputs['Q_bat'])

    WavyChan                        = WavyChan_NMC()
    # results_chan                    = WavyChan_data()
    T_air                           = BatInputs['T_amb'][i]+35
    T_bat                           = BatInputs['T_b_air'][i]+35
    Q_bat                           = BatInputs['Q_bat'][i] * 10
    dt                              = BatInputs['time'][i]-BatInputs['time'][i-1]

    print(Q_bat, 'W', T_bat, 'K', T_air, 'K')

    return T_bat, T_air, Q_bat, dt


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


def plot_BTMS_results(width=10, height=8):
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
    y2_label                            = '$T_{chan}$ ($^\circ$C)'
    y3_label                            = 'COP (-)'
    y4_label                            = 'Power (kW)'

    address     = os.getcwd()
    with open(address + '\BTMS_Case2_results.pkl', 'rb') as f:
        BTMS_results = pickle.load(f)

    # ---------------------------------------------------------------------
    #   plot results
    # ---------------------------------------------------------------------
    fig1 = plt.figure(None)
    fig1.set_size_inches(width, height)
    fig1.suptitle(None)

    plt.subplot(2, 2, 1)
    plt.scatter(BTMS_results['mission_time'], BTMS_results['T_bat'], marker=markers[0], alpha=alpha, s=blob_size, c='none', edgecolors=colors[0], linewidths=line_width)
    plt.grid(visible=True, which='major', linestyle='-', linewidth=0.75, color=colors[0])
    plt.minorticks_on()
    plt.grid(visible=True, which='minor', linestyle='--', linewidth=0.15, color=colors[0])
    plt.ylabel(y1_label)

    plt.subplot(2, 2, 2)
    plt.scatter(BTMS_results['mission_time'], BTMS_results['T_o'] - 273.15, marker=markers[1], alpha=alpha, s=blob_size, c='none', edgecolors=colors[1], linewidths=line_width, label='$T_{o}$')
    plt.scatter(BTMS_results['mission_time'], BTMS_results['T_i'] - 273.15, marker=markers[2], alpha=alpha, s=blob_size, c='none', edgecolors=colors[2], linewidths=line_width, label='$T_{i}$')
    plt.grid(visible=True, which='major', linestyle='-', linewidth=0.75, color=colors[0])
    plt.minorticks_on()
    plt.grid(visible=True, which='minor', linestyle='--', linewidth=0.15, color=colors[0])
    plt.ylabel(y2_label)

    plt.subplot(2, 2, 3)
    plt.scatter(BTMS_results['mission_time'], BTMS_results['COP'], marker=markers[3], alpha=alpha, s=blob_size, c='none', edgecolors=colors[7], linewidths=line_width)
    plt.grid(visible=True, which='major', linestyle='-', linewidth=0.75, color=colors[0])
    plt.minorticks_on()
    plt.grid(visible=True, which='minor', linestyle='--', linewidth=0.15, color=colors[0])
    plt.xlabel(x_label)
    plt.ylabel(y3_label)

    plt.subplot(2, 2, 4)
    plt.scatter(BTMS_results['mission_time'], BTMS_results['Power'], marker=markers[3], alpha=alpha, s=blob_size, c='none', edgecolors=colors[8], linewidths=line_width)
    plt.grid(visible=True, which='major', linestyle='-', linewidth=0.75, color=colors[0])
    plt.minorticks_on()
    plt.grid(visible=True, which='minor', linestyle='--', linewidth=0.15, color=colors[0])
    plt.xlabel(x_label)
    plt.ylabel(y4_label)

    plt.tight_layout()
    plt.show()
    # plt.savefig('fig10.png', dpi=image_resolution, format='png')


# ---------------------------------------------------------------------
#   Main
# ---------------------------------------------------------------------

if __name__ == '__main__':
    # Parameters of a generic HP to be used in the VC model
    # ------------------------------------------------------------------------------
    Comp        = {'Comp: rpm': 3.6e3*2,      'Comp: eta_C_iso': 0.63,    'Comp: Disp': 4.25e-5, 'Comp: eta_C_vol': 0.95, 'Comp: eta_C_mec': 0.88}
    # Evap        = {'Evap: TubeL': 10,       'Evap: TubeN': 5,           "Evap: TubeID": 0.011, 'Evap: FinEff': 0.75,    'Evap: FinRatio': 7.2, 'Evap: Pump': 220}
    Cond        = {'Cond: TubeL': 10,       'Cond: TubeN': 5,           "Cond: TubeID": 0.011, 'Cond: FinEff': 0.75,    'Cond: FinRatio': 7.2, 'Cond: Fan': 365}
    System      = {'Sys: Superheat': 8.5,   'Sys: Subcool': 6.5}
    SysParams   = dict(Comp, **Cond, **System)   # add the dicts together

    m_dot_air                       = 0.47 * 1
    # battery pack layout
    N_series                        = 150
    N_parallel                      = 130

    T_bat, T_air, Q_bat, dt         = BatteryData()

    HP        = Cycle_Performance_BatteryAir(SysParams, dt, Q_bat, T_bat, m_dot_air, T_air, 'air', 'R134a')
    HP.Cycle_Solver()
    HP.PostProcessing()
    # HP.Interpolation_Singularity_Points()
    # show the results
    print('COP:         {:.2f}'.format(HP.COP))
    print('Q_cond:      {:.2f} [W]'.format(HP.Q_condenser))
    print('Q_evap:      {:.2f} [W]'.format(HP.Q_evaporator))
    print('W_comp:      {:.2f} [W]'.format(HP.W_c))
    print('T_bat_updated: {:.2f} [K]'.format(HP.T_bat_updated))
    print('Energy Balance: {:.2f} [-]'.format(HP.E_balance))

