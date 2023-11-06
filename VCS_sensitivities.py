# VCS_sensitivities.py
#
# Created: Nov 2023, C.R. Zhao

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

# Matplotlib Formatting
import matplotlib as mpl
from matplotlib import rcParams, rc, font_manager
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle
from VCS_cycle_performance import Cycle_Performance_Air


# ---------------------------------------------------------------------
#   Component sensitivity: VCS
# ---------------------------------------------------------------------
def get_sensitivity_variables():

    sensitivity_variables = ['m_dot_coolant', 'm_dot_air', 'T_coolant_inlet', 'T_air_inlet']

    return sensitivity_variables


def VCS_sensitivity_parameters(variable_name=None, numbers_of_points=None):

    input_parameters                = input_parameter_Data()
    # parameters                     values                 units
    # ----------------------------------------------------------------------------
    # nominal values
    # ----------------------------------------------------------------------------
    # TEM geometric parameters
    if not variable_name == 'm_dot_coolant':
        m_dot_coolant                       = 1.                  # m^2
    else:
        m_dot_coolant                       = np.linspace(0.5, 1.5, numbers_of_points)

    if not variable_name == 'm_dot_air':
        m_dot_air                           = 5.              # m
    else:
        m_dot_air                           = np.linspace(3., 7., numbers_of_points)

    if not variable_name == 'T_coolant_inlet':
        T_coolant_inlet                     = 20 + 273.15                # -
    else:
        T_coolant_inlet                     = np.linspace(5 + 273.15, 35 + 273.15, numbers_of_points)

    # TEM operating parameters
    if not variable_name == 'T_air_inlet':
        T_air_inlet                         = 40 + 273.15                    # K
    else:
        T_air_inlet                         = np.linspace(30 + 273.15, 50 + 273.15, numbers_of_points)    # 47 - 62 - 77

    # save VCS data
    input_parameters.m_dot_coolant          = m_dot_coolant
    input_parameters.m_dot_air              = m_dot_air
    input_parameters.T_coolant_inlet        = T_coolant_inlet
    input_parameters.T_air_inlet            = T_air_inlet

    return input_parameters


def get_parameter_combinations_for_input(j,input_parameters=None):

    # VCS cooler
    m_dot_coolant                           = input_parameters.m_dot_coolant
    m_dot_air                               = input_parameters.m_dot_air
    T_coolant_inlet                         = input_parameters.T_coolant_inlet
    T_air_inlet                             = input_parameters.T_air_inlet

    parameter_combination               = np.array([m_dot_coolant, m_dot_air, T_coolant_inlet, T_air_inlet])

    for xx in range(len(parameter_combination)):
        if isinstance(parameter_combination[xx], float):
            parameter_combination[xx] = parameter_combination[xx]
        else:
            parameter_combination[xx] = parameter_combination[xx][j]

    return parameter_combination


def run_VCS_sensitivity_analysis(parameter_combination=None):
    """input parameters"""
    # VCS cooler
    m_dot_coolant                           = parameter_combination[0]
    m_dot_air                               = parameter_combination[1]
    T_coolant_inlet                         = parameter_combination[2]
    T_air_inlet                             = parameter_combination[3]

    """constant parameters"""
    Comp        = {'Comp: rpm': 3.6e3,      'Comp: eta_C_iso': 0.63,    'Comp: Disp': 4.25e-5, 'Comp: eta_C_vol': 0.95, 'Comp: eta_C_mec': 0.88}
    Evap        = {'Evap: TubeL': 10,       'Evap: TubeN': 5,           "Evap: TubeID": 0.011, 'Evap: FinEff': 0.75,    'Evap: FinRatio': 7.2, 'Evap: Pump': 220}
    Cond        = {'Cond: TubeL': 10,       'Cond: TubeN': 5,           "Cond: TubeID": 0.011, 'Cond: FinEff': 0.75,    'Cond: FinRatio': 7.2, 'Cond: Fan': 365}
    System      = {'Sys: Superheat': 8.5,   'Sys: Subcool': 6.5}
    SysParams   = dict(Comp, **Evap, **Cond, **System)   # add the dicts together

    # -----------------------------------------------------------------------------------------
    # Calculation
    # -----------------------------------------------------------------------------------------
    HP          = Cycle_Performance_Air(SysParams, m_dot_coolant, T_coolant_inlet, m_dot_air, T_air_inlet, 'water', 'air', 'R134a')
    HP.Cycle_Solver()
    HP.PostProcessing()

    COP                                     = HP.COP
    Q_evaporator                            = HP.Q_evaporator
    W_comp                                  = HP.W_c
    T_coolant_outlet                        = HP.T_coolant_outlet

    # save the outputs
    output_parameters               = results_Data()
    output_parameters.COP                            = COP
    output_parameters.Q_evaporator                   = Q_evaporator
    output_parameters.W_comp                         = W_comp
    output_parameters.T_coolant_outlet               = T_coolant_outlet

    return output_parameters


def obtain_resutls_of_VCS_sensitivities(numbers_of_points=5):

    sensitivity_variables = get_sensitivity_variables()
    print(sensitivity_variables[3])

    for i in range(len(sensitivity_variables)):
        print('variable index is {index}\n'.format(index = i))
        variable_name           = sensitivity_variables[i]
        input_parameters        = VCS_sensitivity_parameters(variable_name=variable_name, numbers_of_points=numbers_of_points)

        save_output_data        = results_Data()
        COP                     = save_output_data.COP
        Q_evaporator            = save_output_data.Q_evaporator
        W_comp                  = save_output_data.W_comp
        T_coolant_outlet        = save_output_data.T_coolant_outlet

        for j in range(numbers_of_points):
            parameter_combination   = get_parameter_combinations_for_input(j, input_parameters=input_parameters)
            output_parameters       = run_VCS_sensitivity_analysis(parameter_combination=parameter_combination)

            COP                     = np.append(COP, output_parameters.COP)
            Q_evaporator            = np.append(Q_evaporator, output_parameters.Q_evaporator)
            W_comp                  = np.append(W_comp, output_parameters.W_comp)
            T_coolant_outlet        = np.append(T_coolant_outlet, output_parameters.T_coolant_outlet)

        # save the data ---------------------------------------------
        save_output_data.COP                           = COP
        save_output_data.Q_evaporator                  = Q_evaporator
        save_output_data.W_comp                        = W_comp
        save_output_data.T_coolant_outlet              = T_coolant_outlet

        # save the results into a pickle file
        pickle_save_name_TEM_HEX = 'VCS sensitives of ' + variable_name + "_3points.pkl"
        save_output_data.save_data(pickle_save_name_TEM_HEX)


def plot_VCS_sensitivities_results(results=None, width=9, height=4):
    # ---------------------------------------------------------------------
    #   Plot Styles
    # ---------------------------------------------------------------------
    # plt.style.use('seaborn-talk')
    mpl.font_manager._get_fontconfig_fonts.cache_clear()

    font_path                           = ['/Users/Chunrong/PycharmProjects/Library/Fonts' ]
    font_files                          = mpl.font_manager.findSystemFonts(fontpaths=font_path)
    for font in font_files:
        mpl.font_manager.fontManager.addfont(font)
    mpl.rc('font', family='Gulliver-Regular')
    mpl_font                            = {'fontname': 'Gulliver-Regular'}

    mpl.rcParams['pdf.fonttype']        = 42
    mpl.rcParams['ps.fonttype']         = 42
    mpl.rcParams['font.family']         = 'Gulliver-Regular'

    fontsize                            = 14
    fontsize_legend                     = 12
    mpl.rc('axes', linewidth=1.75, labelsize=fontsize, titlesize=fontsize)
    mpl.rc('xtick', labelsize=fontsize)
    mpl.rc('ytick', labelsize=fontsize)
    mpl.rc('legend', fontsize=fontsize_legend)
    mpl.rcParams.update({'font.size': fontsize})

    # ---------------------------------------------------------------------
    #   main plot settings
    # ---------------------------------------------------------------------
    line_width                  = 2.25 + 0.5
    line_width_soft             = 2
    marker_size_exp             = 7
    marker_size_exp             = str(marker_size_exp)
    marker_size                 = 8
    marker_size                 = str(marker_size)
    fontsize                    = 14
    fontsize_legend             = 12
    axes_linewidth              = 1.75
    alpha                       = 1.0
    blob_size                   = 1e2 * 1.2
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
    colors = [myblack, myblue, myred, myyellow,  mydarkblue, myorange, mybrown, mygreen, mypurple, mygray]
    markers = ['o', '^', 's', 'p', 'v', '*', 'x']

    # -------------------------------------------------------------------------
    # input results data
    # --------------------------------------------------------------------------
    with open('VCS sensitives of m_dot_coolant_25points.pkl', 'rb') as f:
        results_m_dot_coolant = pickle.load(f)
    with open('VCS sensitives of m_dot_air_25points.pkl', 'rb') as f:
        results_m_dot_air = pickle.load(f)
    with open('VCS sensitives of T_coolant_inlet_25points.pkl', 'rb') as f:
        results_T_coolant_inlet = pickle.load(f)
    with open('VCS sensitives of T_air_inlet_25points.pkl', 'rb') as f:
        results_T_air_inlet = pickle.load(f)

    with open('VCS sensitives of m_dot_coolant_3points.pkl', 'rb') as f:
        results_3points_m_dot_coolant = pickle.load(f)
    with open('VCS sensitives of m_dot_air_3points.pkl', 'rb') as f:
        results_3points_m_dot_air = pickle.load(f)
    with open('VCS sensitives of T_coolant_inlet_3points.pkl', 'rb') as f:
        results_3points_T_coolant_inlet = pickle.load(f)
    with open('VCS sensitives of T_air_inlet_3points.pkl', 'rb') as f:
        results_3points_T_air_inlet= pickle.load(f)

    # -------------------------------------------------------------------------
    # plot the figure
    # --------------------------------------------------------------------------

    # All fronts
    fig        = plt.figure()
    fig.set_size_inches(width, height)

    # legend labels
    label1 = '$\dot{m}_{coolant}$'
    label2 = '$\dot{m}_{air}$'
    label3 = '$T_{coolant,in}$'
    label4 = '$T_{air,in}$'

    # plot tick labels
    x_label_1       = '$Q_{EVAP}$ [kW]'
    y_label_1       = '$W_{COMP}$ [kW]'
    x_label_2       = '$COP$ [-]'
    y_label_2       = '$T_{coolant,out}$ [$^\circ$C]'

    # ---------------------------------------------------------------------------------------------

    plt.subplot(1,2,1)
    plt.scatter(results_3points_m_dot_coolant.Q_evaporator / 1e3,   results_3points_m_dot_coolant.W_comp / 1e3, marker=markers[0], c='none', edgecolors=colors[0], linewidth=line_width, alpha=1, s=blob_size, label=label1)
    plt.scatter(results_3points_m_dot_air.Q_evaporator / 1e3,   results_3points_m_dot_air.W_comp / 1e3, marker=markers[1], c='none', edgecolors=colors[1], linewidth=line_width, alpha=1, s=blob_size, label=label2)
    plt.scatter(results_3points_T_coolant_inlet.Q_evaporator / 1e3,     results_3points_T_coolant_inlet.W_comp / 1e3, marker=markers[2], c='none', edgecolors=colors[2], linewidth=line_width, alpha=1, s=blob_size, label=label3)
    plt.scatter(results_3points_T_air_inlet.Q_evaporator / 1e3,     results_3points_T_air_inlet.W_comp / 1e3, marker=markers[3], c='none', edgecolors=colors[3], linewidth=line_width, alpha=1, s=blob_size, label=label4)
    plt.plot(results_m_dot_coolant.Q_evaporator / 1e3,   results_m_dot_coolant.W_comp / 1e3,        linestyle='solid',       color=colors[0], linewidth=line_width)
    plt.plot(results_m_dot_air.Q_evaporator / 1e3,   results_m_dot_air.W_comp / 1e3,    linestyle='dashdot',     color=colors[1], linewidth=line_width)
    plt.plot(results_T_coolant_inlet.Q_evaporator / 1e3,     results_T_coolant_inlet.W_comp / 1e3,      linestyle='dashed',      color=colors[2], linewidth=line_width)
    plt.plot(results_T_air_inlet.Q_evaporator / 1e3,     results_T_air_inlet.W_comp / 1e3,      linestyle='dotted',      color=colors[3], linewidth=line_width)
    plt.grid(visible=True, which='major', linestyle='-', linewidth=0.75, color=colors[0])
    plt.minorticks_on()
    plt.grid(visible=True, which='minor', linestyle='--', linewidth=0.15, color=colors[0])
    plt.xlabel(x_label_1, fontsize=fontsize+4)
    plt.ylabel(y_label_1, fontsize=fontsize+4)
    plt.legend(ncol=1)


    plt.subplot(1,2,2)
    plt.scatter(results_3points_m_dot_coolant.COP,   results_3points_m_dot_coolant.T_coolant_outlet - 273.15, marker=markers[0], c='none', edgecolors=colors[0], linewidth=line_width, alpha=1, s=blob_size, label=label1)
    plt.scatter(results_3points_m_dot_air.COP,   results_3points_m_dot_air.T_coolant_outlet - 273.15, marker=markers[1], c='none', edgecolors=colors[1], linewidth=line_width, alpha=1, s=blob_size, label=label2)
    plt.scatter(results_3points_T_coolant_inlet.COP,     results_3points_T_coolant_inlet.T_coolant_outlet - 273.15,  marker=markers[2], c='none', edgecolors=colors[2], linewidth=line_width, alpha=1, s=blob_size, label=label3)
    plt.scatter(results_3points_T_air_inlet.COP,     results_3points_T_air_inlet.T_coolant_outlet - 273.15,  marker=markers[3], c='none', edgecolors=colors[3], linewidth=line_width, alpha=1, s=blob_size, label=label4)
    plt.plot(results_m_dot_coolant.COP,   results_m_dot_coolant.T_coolant_outlet - 273.15,        linestyle='solid',       color=colors[0], linewidth=line_width)
    plt.plot(results_m_dot_air.COP,   results_m_dot_air.T_coolant_outlet - 273.15,    linestyle='dashdot',     color=colors[1], linewidth=line_width)
    plt.plot(results_T_coolant_inlet.COP,     results_T_coolant_inlet.T_coolant_outlet - 273.15,      linestyle='dashed',      color=colors[2], linewidth=line_width)
    plt.plot(results_T_air_inlet.COP,     results_T_air_inlet.T_coolant_outlet - 273.15,      linestyle='dotted',      color=colors[3], linewidth=line_width)
    plt.grid(visible=True, which='major', linestyle='-', linewidth=0.75, color=colors[0])
    plt.minorticks_on()
    plt.grid(visible=True, which='minor', linestyle='--', linewidth=0.15, color=colors[0])
    plt.xlabel(x_label_2, fontsize=fontsize+4)
    plt.ylabel(y_label_2, fontsize=fontsize+4)

    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------------
#   Data storage
# ---------------------------------------------------------------------
class input_parameter_Data:
    def __init__(self):
        self.m_dot_coolant                              = []
        self.m_dot_air                                  = []
        self.T_coolant_inlet                            = []
        self.T_air_inlet                                = []

    def save_data(self, filename):
        with open(filename, 'wb') as f:
            pickle.dump(self, f)


class results_Data:
    def __init__(self):
        self.COP                              = []
        self.Q_evaporator                     = []
        self.W_comp                           = []
        self.T_coolant_outlet                 = []

    def save_data(self, filename):
        with open(filename, 'wb') as f:
            pickle.dump(self, f)


# ---------------------------------------------------------------------
#   Main
# ---------------------------------------------------------------------
if __name__ == '__main__':

    numbers_of_points           = 3
    # obtain_resutls_of_VCS_sensitivities(numbers_of_points=numbers_of_points)
    plot_VCS_sensitivities_results()
