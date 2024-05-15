# Paper_4_VCS_design.py
#
# Created: Feb. 2024, C.R. Zhao

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import os
import pickle
from ACHP_codes.Components.WavyChan import WavyChan_NMC, WavyChan_data
from ACHP_codes.Components.ThermalStorage import LatentHeatThermalEnergyStorage
from ACHP_codes.CycleTests.SL_BatteryThermalManagementSystem import BTMS_VCS_SL
from matplotlib import rcParams, rc, font_manager
import matplotlib as mpl
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------
#   Functions
# ---------------------------------------------------------------------
def Main_Design(m_dot_h, N_series, N_parallel):
    address                         = os.getcwd()
    with open(address + '\BatInputs_for_VCS.pkl', 'rb') as f:
        BatInputs                   = pickle.load(f)

    # print(BatInputs['Q_bat'])

    WavyChan                        = WavyChan_NMC()
    results_chan                    = WavyChan_data()
    T_i                             = BatInputs['T_amb'][0]+35
    T_amb                           = BatInputs['T_amb']+35
    T_current                       = BatInputs['T_b_air'][0]+35
    # Q_heat_gen                      = BatInputs['Q_bat'] * 10
    VCS_results                     = {'T_bat': [], 'T_i': [], 'T_o': [], 'COP': [], 'Power': [], 'mission_time': []}

    LHTES                           = LatentHeatThermalEnergyStorage()
    r_melt_front                    = 0.
    Heat_accumulated                = 0.
    for i in range(1, len(BatInputs['Q_bat'])):
        # ---------------------------------------------------------------------
        #   Heat Exchange between the Battery and Wavy Channel
        # ---------------------------------------------------------------------
        """perform wavy channel calculation"""
        T_current, Q_conv, T_o, delta_P_chan                = WavyChan.compute_net_generated_battery_heat_module(Q_heat_gen=BatInputs['Q_bat'][i],
                                                                   m_dot_coolant=m_dot_h, T_cell=T_current,
                                                                   T_i=T_i, dt=BatInputs['time'][i]-BatInputs['time'][i-1],
                                                                   N_series=N_series, N_parallel=N_parallel)

        P_chan, m_chan, eff_WavyChan                        = WavyChan.compute_power_and_mass(m_dot_coolant=m_dot_h,
                                                                              N_series=N_series, N_parallel=N_parallel)

        Q_wavychannel                                       = Q_conv * 10
        T_bat                                               = T_current

        T_o_coolant, Heat_accumulated, m_LHTES \
            = LHTES.SimplifiedLHTES(m_dot_coolant=m_dot_h, Q_wavychannel=Q_wavychannel, T_i_coolant=T_o,
                                    dt=BatInputs['time'][i]-BatInputs['time'][i-1], Heat_accumulated=Heat_accumulated, Q_evap_design=1e4)

        COP, Power, T_chan_o         = BTMS_VCS_SL(Q_wavychannel=Q_wavychannel, m_dot_g=m_dot_h, T_bat=T_bat, T_amb=T_amb[i])
        # -----------------------------------------------------------------
        VCS_results['mission_time'].append(BatInputs['time'][i]-BatInputs['time'][0])
        VCS_results['T_bat'].append(T_current)
        VCS_results['T_i'].append(T_i)
        VCS_results['T_o'].append(T_o)
        VCS_results['COP'].append(COP)
        VCS_results['Power'].append(Power)
        # -----------------------------------------------------------------
        T_i                                                 = T_chan_o

    address     = os.getcwd()
    save_file(VCS_results, address + '\BTMS_Case2_results.pkl')


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


# ------------------------------------------------------------------------------------
if __name__ == '__main__':
    # battery pack layout
    N_series                        = 150
    N_parallel                      = 130
    # coolant mass flow rate
    m_dot_h                         = 1.

    Main_Design(m_dot_h, N_series, N_parallel)

    plot_BTMS_results()
