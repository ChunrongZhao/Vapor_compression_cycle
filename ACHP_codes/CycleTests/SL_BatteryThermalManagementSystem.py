"""This code is for secondary loop cycle in AC Mode"""
# from __future__                                 import division, absolute_import, print_function
from ACHP_codes.Cycles.Cycle                    import  BTMS_SecondaryCycleClass
from ACHP_codes.ACHP_Tools.Plots                import PlotsClass, SL_PlotsClass, BTMS_SL_PlotsClass
from time                                       import time
import pickle
import CoolProp.CoolProp                        as CP
from pathlib                        import Path
import os
from ACHP_codes.Components.WavyChan import WavyChan_NMC, WavyChan_data
from ACHP_codes.Components.ThermalStorage import LatentHeatThermalEnergyStorage
# ------------------------------------------------------------------------------------
def CoolProp_AbstractState_test1():
    Backend_SLF             = 'INCOMP'
    SecLoopFluid            = 'MEG'
    AS_SLF                  = CP.AbstractState(Backend_SLF, SecLoopFluid)

    T_max                   = AS_SLF.Tmax()
    T_min                   = AS_SLF.Tmin()

    print('T_min = ', T_min, ' [K] or ', T_min-273.15, '[oC]')
    print('T_max = ', T_max, ' [K] or ', T_max-273.15, '[oC]')

    AS_SLF.update(CP.PT_INPUTS, 101e3, 284.150285)
    h                       = AS_SLF.hmass()
    print('*********************************************')
    print('h = ', h, ' [J/kg]')


def BTMS_VCS_SL(Q_wavychannel=5000, m_dot_g=0.38, T_bat=289.15, T_amb=308.15):
    # Instantiate the class
    Cycle                                   = BTMS_SecondaryCycleClass()

    # ----------------------------------------------------------------------------------
    # Cycle parameters
    # ----------------------------------------------------------------------------------

    Cycle.Verbosity                         = 0                     # the idea here is to have different levels of debug output
    Cycle.ImposedVariable                   = 'Subcooling'          # or this could be 'Charge' for imposed charge
    Cycle.DT_sc_target                      = 7.0
    Cycle.CycleType                         = 'Secondary'
    Cycle.Ref                               = 'R410A'               # todo variable
    Cycle.Backend                           = 'HEOS'        # Backend for refrigerant properties calculation: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    Cycle.Oil                               = 'POE32'
    Cycle.SecLoopFluid                      = 'MEG'                 # todo variable
    Cycle.MassFrac_SLF                      = 0.21                  # Mass fraction of incompressible SecLoopFluid [i.e MEG-20%]
    Cycle.Backend_SLF                       = 'INCOMP'              # backend of SecLoopFluid
    Cycle.IHXType                           = 'PHE'                 # or could be 'Coaxial'
    Cycle.shell_pressure                    = 'low-pressure'
    Cycle.Mode                              = 'AC'
    Cycle.DT_sh                             = 5                     # superheat
    Cycle.T_bat                             = T_bat                 # todo check
    Cycle.EvapSolver                        = 'Moving-Boundary'     # choose the type of Evaporator solver scheme ('Moving-Boundary' or 'Finite-Element')
    Cycle.EvapType                          = 'Fin-tube'            # if EvapSolver = 'Moving-Boundary', choose the type of evaporator ('Fin-tube' or 'Micro-channel')
    Cycle.CondSolver                        = 'Moving-Boundary'     # choose the type of Condenser solver scheme ('Moving-Boundary' or 'Finite-Element')
    Cycle.CondType                          = 'Fin-tube'            # if CondSolver = 'Moving-Boundary', choose the type of condenser ('Fin-tube' or 'Micro-channel')
    Cycle.Update()

    # ----------------------------------------------------------------------------------
    #     Charge correction parameters (activate by setting Cycle.ImposedVariable to 'Charge' and Cycle.ChargeMethod to either 'One-point' or 'Two-point')
    # ----------------------------------------------------------------------------------
    Cycle.C                                 = 0                     # [kg]
    Cycle.K                                 = 0                     # [kg]
    Cycle.w_ref                             = 0                     # [-]

    # ----------------------------------------------------------------------------------
    #       Compressor parameters
    # ----------------------------------------------------------------------------------

    # A 3 ton cooling capacity compressor map
    if Cycle.Ref == 'R410A':
        M   = [217.3163128, 5.094492028, -0.593170311, 4.38E-02, -2.14E-02, 1.04E-02, 7.90E-05, -5.73E-05, 1.79E-04, -8.08E-05]
        P   = [-561.3615705, -15.62601841, 46.92506685, -0.217949552, 0.435062616, -0.442400826, 2.25E-04, 2.37E-03, -3.32E-03, 2.50E-03]
    elif Cycle.Ref == 'R134a':
        M   = [217.3163128, 5.094492028, -0.593170311, 4.38E-02, -2.14E-02, 1.04E-02, 7.90E-05, -5.73E-05, 1.79E-04, -8.08E-05],
        P   = [-561.3615705, -15.62601841, 46.92506685, -0.217949552, 0.435062616, -0.442400826, 2.25E-04, 2.37E-03, -3.32E-03, 2.50E-03],

    params  = {
            'M':                            M,
            'P':                            P,
            'Ref':                          Cycle.Ref,              # refrigerant
            'Oil':                          Cycle.Oil,              # Compressor lubricant oil
            'shell_pressure':               Cycle.shell_pressure,   # Compressor shell pressure
            'fp':                           0.15,                   # Fraction of electrical power lost as heat to ambient
            'V_dot_ratio':                  1.0,                    # Displacement Scale factor to up- or downsize compressor (1=original)
            'V_oil_sump':                   0,                      # Volume of oil in the sump
            'Verbosity':                    0,                      # How verbose should the debugging statements be [0 to 10]
            }
    Cycle.Compressor.Update(**params)

    # ----------------------------------------------------------------------------------
    # Condenser parameters
    # ----------------------------------------------------------------------------------
    Cycle.Condenser.Fins.Tubes.N_Tubes_per_bank     = 24            # number of tubes per bank=row
    Cycle.Condenser.Fins.Tubes.N_bank               = 1             # number of banks/rows
    Cycle.Condenser.Fins.Tubes.N_circuits           = 3
    Cycle.Condenser.Fins.Tubes.L_tube               = 2.252
    Cycle.Condenser.Fins.Tubes.OD                   = 0.00913
    Cycle.Condenser.Fins.Tubes.ID                   = 0.00849
    Cycle.Condenser.Fins.Tubes.Pl                   = 0.0191        # distance between center of tubes in flow direction
    Cycle.Condenser.Fins.Tubes.Pt                   = 0.0254        # distance between center of tubes orthogonal to flow direction
    Cycle.Condenser.Fins.Tubes.k_w                  = 237           # wall thermal conductivity (i.e pipe material)

    Cycle.Condenser.Fins.Fins.FPI                   = 25            # Number of fins per inch
    Cycle.Condenser.Fins.Fins.Pd                    = 0.001         # 2* amplitude of wavy fin
    Cycle.Condenser.Fins.Fins.xf                    = 0.001         # 1/2 period of fin
    Cycle.Condenser.Fins.Fins.t                     = 0.00011       # Thickness of fin material
    Cycle.Condenser.Fins.Fins.k_fin                 = 237           # Thermal conductivity of fin material

    Cycle.Condenser.Fins.Air.V_dot_ha               = 1.7934*1      # rated volumetric flowrate Cycle.Condenser.Fins.Air.Tmean=308.15
    Cycle.Condenser.Fins.Air.T_db                   = T_amb         # Dry Bulb Temperature
    Cycle.Condenser.Fins.Air.p                      = 101325        # Condenser Air pressures in Pa
    Cycle.Condenser.Fins.Air.RH                     = 0.51          # Relative Humidity
    Cycle.Condenser.Fins.Air.RH_mean                = 0.51
    Cycle.Condenser.Fins.Air.FanPower               = 260

    Cycle.Condenser.FinsType                        = 'WavyLouveredFins'   # Choose fin Type: 'WavyLouveredFins' or 'HerringboneFins' or 'PlainFins'
    Cycle.Condenser.Verbosity                       = 0

    # -----------------------------------------------------------------------------------------
    #       Expansion device Parameters
    # -----------------------------------------------------------------------------------------
    params  = {
            'ExpType':                              'Linear-TXV',           # expansion device type
            'T_sh_static':                          4,                      # static superheat
            'T_sh_max':                             6,                      # maximum superheat
            'D':                                    0.0006604,              # inside diameter [m]
            'C':                                    1.2656e-6,              # constant from manufacturer [m^2/K]
            'Adj':                                  0.7630,                 # Adjust the diameter (tuning factor)
        }
    Cycle.ExpDev.Update(**params)

    # ----------------------------------------------------------------------------------
    #       IHX Plate Heat Exchanger
    # ----------------------------------------------------------------------------------
    params  = {
            'p_in_h':                               300000,
            # Geometric parameters
            'HXType':                               'Plate-HX',
            'Bp':                                   0.117,
            'Lp':                                   0.300,          # Center-to-center distance between ports
            'Nplates':                              46,             # todo check
            'PlateAmplitude':                       0.001,          # [m]
            'PlateThickness':                       0.0003,         # [m]
            'PlateConductivity':                    15.0,           # [W/m-K]
            'Rp':                                   1.0,            # [microns] Surface roughness
            'PlateWavelength':                      0.00628,        # [m]
            'InclinationAngle':                     3.14159/3,      # [rad]
            'MoreChannels':                         'Hot',          # Which stream gets the extra channel, 'Hot' or 'Cold'
            'Verbosity':                            0,
            }
    Cycle.PHEIHX.Update(**params)

    # ----------------------------------------------------------------------------------
    # Battery Pack and WavyChannel parameters
    # ----------------------------------------------------------------------------------
    params  = {
            'Q_wavychannel':                        Q_wavychannel,  # todo variable
            'p_in_g':                               101325,
            'DP_g':                                 2000,
            'eta':                                  0.95,
            'Verbosity':                            0,
            'm_dot_g':                              m_dot_g,        # todo variable
            }
    Cycle.WavyChannel.Update(**params)

    # ----------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------
    #       Secondary Loop Pump
    # ----------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------
    params  = {
            'eta':                                  0.5,            # Pump+motor efficiency
            'm_dot_g':                              m_dot_g,           # Flow Rate kg/s
            'p_in_g':                               300000,
            'Verbosity':                            0,
            }
    Cycle.Pump.Update(**params)

    # ----------------------------------------------------------------------------------
    #       Line Set Supply and Return Parameters
    # ----------------------------------------------------------------------------------
    params  = {
            'L':                                    5,
            'k_tube':                               0.19,
            't_insul':                              0.02,
            'k_insul':                              0.036,
            'T_air':                                297,
            'h_air':                                0.0000000001,
            'LineSetOption':                        'On'
            }
    Cycle.LineSetSupply.Update(**params)
    Cycle.LineSetReturn.Update(**params)

    Cycle.LineSetSupply.OD                          = 0.01905
    Cycle.LineSetSupply.ID                          = 0.017526
    Cycle.LineSetReturn.OD                          = 0.01905
    Cycle.LineSetReturn.ID                          = 0.017526

    # ----------------------------------------------------------------------------------
    #       Line Set Parameters Suction and Liquid
    # ----------------------------------------------------------------------------------
    params  = {
            'L':                                    7.6,
            'k_tube':                               0.19,
            't_insul':                              0.02,
            'k_insul':                              0.036,
            'T_air':                                297,
            'h_air':                                0.0000000001,
            'LineSetOption':                        'Off'
            }

    Cycle.LineSetLiquid.Update(**params)
    Cycle.LineSetSuction.Update(**params)

    Cycle.LineSetLiquid.OD                          = 0.009525
    Cycle.LineSetLiquid.ID                          = 0.007986
    Cycle.LineSetSuction.OD                         = 0.01905
    Cycle.LineSetSuction.ID                         = 0.017526

    # ----------------------------------------------------------------------------------
    #       Line Set Discharge Parameters
    # ----------------------------------------------------------------------------------
    params  = {
            'L':                                    0.3,                    # tube length in m
            'k_tube':                               0.19,
            't_insul':                              0,                      # no insulation
            'k_insul':                              0.036,
            'T_air':                                297,
            'h_air':                                0.0000000001,
            'LineSetOption':                        'Off'
            }

    Cycle.LineSetDischarge.Update(**params)

    Cycle.LineSetDischarge.OD                       = 0.009525
    Cycle.LineSetDischarge.ID                       = 0.007986

    # ----------------------------------------------------------------------------------
    # Now solve
    t1                                              = time()
    Cycle.PreconditionedSolve()

    # Outputs
    print('Took ' + str(time() - t1) + ' seconds to run Cycle model')
    print('Cycle coefficient of system performance is ' + str(Cycle.COSP))
    print('Cycle refrigerant charge is ' + str(Cycle.Charge) + ' kg')

    # Now do cycle plotting
    Plot_result     = 0
    if Plot_result:
        plot                                            = BTMS_SL_PlotsClass()
        plot.TSOverlay(Cycle)
        plot.PHOverlay(Cycle)

    return Cycle.COSP, Cycle.Power, Cycle.PHEIHX.T_out_h, Cycle.PHEIHX.Q


def Main_Design(m_dot_h=1., N_series=150, N_parallel=130):
    address                         = os.getcwd()
    parent2          = os.path.dirname(os.path.dirname(address))
    with open(parent2 + '\BatInputs_for_VCS.pkl', 'rb') as f:
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

    # for i in range(1, len(BatInputs['Q_bat'])):
    i = 1
    # ---------------------------------------------------------------------
    #   Heat Exchange between the Battery and Wavy Channel
    # ---------------------------------------------------------------------
    """perform wavy channel calculation"""
    T_current, Q_conv, T_o, delta_P_chan                = WavyChan.compute_net_generated_battery_heat_module(Q_heat_gen=BatInputs['Q_bat'][i],
                                                               m_dot_coolant=m_dot_h, T_cell=T_current,
                                                               T_i=T_i, dt=BatInputs['time'][i]-BatInputs['time'][i-1],
                                                               N_series=N_series, N_parallel=N_parallel)

    Q_wavychannel                                       = Q_conv * 10
    T_bat                                               = T_current

    T_o_coolant, Heat_accumulated, m_LHTES \
        = LHTES.SimplifiedLHTES(m_dot_coolant=m_dot_h, Q_wavychannel=Q_wavychannel, T_i_coolant=T_o,
                                dt=BatInputs['time'][i]-BatInputs['time'][i-1], Heat_accumulated=Heat_accumulated, Q_evap_design=1e4)

    COP, Power, T_chan_o, Q         = BTMS_VCS_SL(Q_wavychannel=Q_wavychannel, m_dot_g=m_dot_h, T_bat=T_bat, T_amb=T_amb[i])

    return T_o_coolant, T_chan_o, m_LHTES, COP, Q


# ------------------------------------------------------------------------------------
if __name__ == '__main__':
    # print('why')
    # CoolProp_AbstractState_test1()

    # T_o_coolant, T_chan_o, m_LHTES, COP, Q = Main_Design()
    COP, Power, T_o_coolant, Q = BTMS_VCS_SL()
    print(COP, Power, T_o_coolant, Q)

    # AS_SLF         = CP.AbstractState('INCOMP', 'MEG')
    # AS_SLF.update(CP.PT_INPUTS, 101325, 293.15)
    # h              = AS_SLF.hmass()
    # print(h)
