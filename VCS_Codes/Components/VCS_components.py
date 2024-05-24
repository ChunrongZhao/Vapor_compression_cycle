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

# ---------------------------------------------------------------------
#   Condenser
# ---------------------------------------------------------------------
class CondenserClass():
    def __init__(self, **kwargs):
        # Load the parameters passed in using the dictionary
        self.__dict__.update(kwargs)

    def OutputList(self):
        """
            Return a list of parameters for this component for further output

            It is a list of tuples, and each tuple is formed of items:
                [0] Description of value
                [1] Units of value
                [2] The value itself
        """
        return [
            ('Volumetric flow rate of humid air',   'm^3/s',        self.Fins.Air.V_dot_ha),
            ('Inlet Dry bulb temp',                 'K',            self.T_in_a),
            ('Inlet Air pressure',                  'Pa',           self.Fins.Air.p),
            ('Inlet Air Relative Humidity',         '-',            self.Fins.Air.RH),
            ('Tubes per bank',                      '-',            self.Fins.Tubes.N_Tubes_per_bank),
            ('Number of banks',                     '-',            self.Fins.Tubes.N_bank),
            ('Number circuits',                     '-',            self.Fins.Tubes.N_circuits),
            ('Length of tube',                      'm',            self.Fins.Tubes.L_tube),
            ('Tube OD',                             'm',            self.OD),
            ('Tube ID',                             'm',            self.ID),
            ('Tube Long. Pitch',                    'm',            self.Fins.Tubes.Pl),
            ('Tube Transverse Pitch',               'm',            self.Fins.Tubes.Pt),
            ('Tube Conductivity',                   'W/m-K',        self.Fins.Tubes.k_w),
            ('Fins per inch',                       '1/in',         self.Fins.Fins.FPI),
            ('Fin waviness pd',                     'm',            self.Fins.Fins.Pd),
            ('Fin waviness xf',                     'm',            self.Fins.Fins.xf),
            ('Fin thickness',                       'm',            self.Fins.Fins.t),
            ('Fin Conductivity',                    'W/m-K',        self.Fins.Fins.k_fin),
            ('Fins Type',                           '-',            self.FinsType),
            ('Q Total',                             'W',            self.Q),
            ('Q Superheat',                         'W',            self.Q_superheat),
            ('Q Two-Phase',                         'W',            self.Q_2phase),
            ('Q Subcool',                           'W',            self.Q_subcool),
            ('Inlet Temp',                          'K',            self.T_in_r),
            ('Outlet Temp',                         'K',            self.T_out_r),
            ('Pressure Drop Total',                 'Pa',           self.DP_r),
            ('Pressure Drop Superheat',             'Pa',           self.DP_r_superheat),
            ('Pressure Drop Two-Phase',             'Pa',           self.DP_r_2phase),
            ('Pressure Drop Subcool',               'Pa',           self.DP_r_subcool),
            ('Charge Total',                        'kg',           self.Charge),
            ('Charge Superheat',                    'kg',           self.Charge_superheat),
            ('Charge Two-Phase',                    'kg',           self.Charge_2phase),
            ('Charge Subcool',                      'kg',           self.Charge_subcool),
            ('Mean HTC Superheat',                  'W/m^2-K',      self.h_r_superheat),
            ('Mean HTC Two-phase',                  'W/m^2-K',      self.h_r_2phase),
            ('Mean HTC Subcool',                    'W/m^2-K',      self.h_r_subcool),
            ('Wetted Area Fraction Superheat',      '-',            self.w_superheat),
            ('Wetted Area Fraction Two-phase',      '-',            self.w_2phase),
            ('Wetted Area Fraction Subcool',        '-',            self.w_subcool),
            ('Mean Air HTC',                        'W/m^2-K',      self.Fins.h_a),
            ('Surface Effectiveness',               '-',            self.Fins.eta_a),
            ('Air-side area (fin+tubes)',           'm^2',          self.Fins.A_a*self.h_a_tuning),
            ('Mass Flow rate of Dry Air',           'kg/s',         self.Fins.m_dot_da),
            ('Mass Flow rate of Humid Air',         'kg/s',         self.Fins.m_dot_ha),
            ('Pressure Drop Air-side',              'Pa',           self.Fins.dP_a),
            ('Subcooling',                          'K',            self.DT_sc)
        ]

    def Update(self, **kwargs):
        # Update the parameters passed in using the dictionary
        self.__dict__.update(kwargs)

    def Calculate(self):
        # Set tuning factors to 1 in case not given by user
        if not hasattr(self, 'h_a_tuning'):
            self.h_a_tuning                     = 1
        if not hasattr(self, 'h_tp_tuning'):
            self.h_tp_tuning                    = 1
        if not hasattr(self, 'DP_tuning'):
            self.DP_tuning                      = 1

        # AbstractState
        AS                                      = self.AS

        # Retrieve some parameters from nested structures
        # for code compactness
        self.ID                                 = self.Fins.Tubes.ID
        self.OD                                 = self.Fins.Tubes.OD
        self.L_tube                             = self.Fins.Tubes.L_tube
        self.N_Tubes_per_bank                   = self.Fins.Tubes.N_Tubes_per_bank
        self.N_bank                             = self.Fins.Tubes.N_bank
        self.N_circuits                         = self.Fins.Tubes.N_circuits
        self.P_l                                = self.Fins.Tubes.Pl
        self.T_in_a                             = self.Fins.Air.T_db
        self.k_w                                = self.Fins.Tubes.k_w           # thermal conductivity of tube wall

        # Bubble and dew temperatures (same for fluids without glide)
        AS.update(CP.PQ_INPUTS, self.p_sat_r, 0.0)
        self.T_bubble                           = AS.T()                    # [K]
        self.h_l                                = AS.hmass()                # [J/kg]
        self.cp_sat_L                           = AS.cpmass()               # [J/kg-K]

        AS.update(CP.PQ_INPUTS, self.p_sat_r, 1.0)
        self.T_dew                              = AS.T()                    # [K]
        self.h_v                                = AS.hmass()                # [J/kg]

        # Calculate an effective length of circuit if circuits are not all the same length
        TotalLength                             = self.L_tube * self.N_Tubes_per_bank * self.N_bank
        self.L_circuit                          = TotalLength / self.N_circuits

        self.V_r                                = pi * self.ID**2 / 4.0 * self.L_circuit * self.N_circuits  # total refrigerant volume
        self.A_r_wetted                         = pi * self.ID * self.N_circuits * self.L_circuit
        self.G_r                                = self.m_dot_r / (self.N_circuits * pi * self.ID**2 / 4.0)

        # Thermal resistance at the wall
        self.R_w                                = log(self.OD / self.ID) / (2 * pi * self.k_w * self.L_circuit * self.N_circuits)

        # Define known parameters
        AS.update(CP.PT_INPUTS, self.p_sat_r, self.T_in_r)
        self.h_in_r                             = AS.hmass()                # [J/kg]
        self.s_in_r                             = AS.smass()                # [J/kg-K]

        # Definitely have a superheated portion
        # First try to run with a full superheated section from Tin_r to Tdew
        self._Superheat_Forward(T_sh_out_r=self.T_dew)
        # Maybe have a full two-phase section
        # First try to run with a full two-phase section from quality of 1 to quality of 0
        self._TwoPhase_Forward()
        # If we have already used too much of the HX (max possible sum of w is 1.0)
        if self.w_superheat >= 1:
            # There is partial superheated region
            brentq(self._Superheat_Forward, self.T_dew, self.T_in_r)
            # Zero out all the two-phase parameters
            self.Q_2phase                       = 0.0
            self.DP_r_2phase                    = 0.0
            self.Charge_2phase                  = 0.0
            self.w_2phase                       = 0.0
            self.h_r_2phase                     = 0.0
            # Zero out all the subcooled parameters
            self.Q_subcool                      = 0.0
            self.DP_r_subcool                   = 0.0
            self.Charge_subcool=0.0
            self.w_subcool                      = 0.0
            self.h_r_subcool                    = 0.0
            self.existsSubcooled                = False
        elif self.w_2phase + self.w_superheat > 1:
            # There is no subcooled portion, solve for outlet quality
            brentq(self._TwoPhase_Forward, 0.0000001, 0.9999999)
            # Zero out all the subcooled parameters
            self.Q_subcool                      = 0.0
            self.DP_r_subcool                   = 0.0
            self.Charge_subcool                 = 0.0
            self.w_subcool                      = 0.0
            self.h_r_subcool                    = 0.0
            self.existsSubcooled                = False
        else:
            # By definition then we have a subcooled portion, solve for it
            self.existsSubcooled                = True
            self._Subcool_Forward()

        # Overall calculations
        self.Q                                  = self.Q_superheat + self.Q_2phase + self.Q_subcool
        self.DP_r                               = self.DP_r_superheat + self.DP_r_2phase + self.DP_r_subcool
        self.DP_r                               = self.DP_r * self.DP_tuning        # correcting the pressure drop
        self.Charge                             = self.Charge_2phase + self.Charge_subcool + self.Charge_superheat

        if self.existsSubcooled == True:
            # Outlet condition is in the subcooled zone
            AS.update(CP.PT_INPUTS, self.p_sat_r, self.T_out_r)
            self.h_out_r                        = AS.hmass()            # [J/kg]
            self.s_out_r                        = AS.smass()            # [J/kg-K]
        elif self.w_superheat >= 1:
            # Outlet condition is in the superheated zone
            AS.update(CP.PT_INPUTS, self.p_sat_r, self.T_out_r)
            self.h_out_r                        = AS.hmass()            # [J/kg]
            self.s_out_r                        = AS.smass()            # [J/kg-K]
            self.DT_sc                          = 0.0                   # [K]
        else:
            # Outlet conditions is in the two-phase zone
            self.T_out_r                        = self.x_out_2phase * self.T_dew + (1 - self.x_out_2phase) * self.T_bubble
            AS.update(CP.QT_INPUTS, 0.0, self.T_out_r)
            h_l                                 = AS.hmass()            # [J/kg]
            s_l                                 = AS.smass()            # [J/kg-K]
            AS.update(CP.QT_INPUTS, 1.0, self.T_out_r)
            h_v                                 = AS.hmass()            # [J/kg]
            s_v                                 = AS.smass()            # [J/kg-K]
            self.h_out_r                        = h_l + self.x_out_2phase * (h_v - h_l)
            self.s_out_r                        = s_l + self.x_out_2phase * (s_v - s_l)
            # Use the effective subcooling
            self.DT_sc                          = self.DT_sc_2phase

        # Calculate the mean outlet air temperature [K]
        self.T_out_a        = self.T_in_a - self.Q / (self.Fins.cp_da * self.Fins.m_dot_da)
        self.h_mean_r       = self.w_2phase * self.h_r_2phase + self.w_superheat * self.h_r_superheat + self.w_subcool * self.h_r_subcool
        self.UA_r           = self.h_mean_r * self.A_r_wetted
        self.UA_a           = (self.Fins.h_a * self.h_a_tuning) * self.Fins.A_a * self.Fins.eta_a
        self.UA_w           = 1 / self.R_w

        # added by Chunrong at 24/11/2023
        # rho_pipe            = 2710          # aluminium
        # rho_fin             = 8960          # copper
        # self.L_COND, self.H_COND, self.W_COND, self.m_COND, self.V_COND = HerringboneFins_Dimensions(self.height_cond, self.L_tube, self.N_bank,
        #                 self.P_l, self.secTheta, self.t_f, self.N_fin, rho_fin, self.ID, self.OD, self.N_Tubes_per_bank, rho_pipe, self.Charge)

        return self.T_out_a

    def _Superheat_Forward(self, T_sh_out_r):
        # **********************************************************************
        #                      SUPERHEATED PART
        # **********************************************************************
        # AbstractState
        AS                  = self.AS
        # Dew temperature for constant pressure cooling to saturation
        T_dew               = T_sh_out_r    # self.T_dew
        T_bubble            = self.T_bubble

        # Average fluid temps are used for the calculation of properties
        # Average temp of refrigerant is average of sat. temp and outlet temp
        # Secondary fluid is air over the fins
        """the superheated region (SH)"""
        self.f_r_superheat, self.h_r_superheat, self.Re_r_superheat = f_h_1phase_Tube(self.m_dot_r / self.N_circuits,
            self.ID, (T_dew + self.T_in_r) / 2.0, self.p_sat_r, self.AS)

        AS.update(CP.PT_INPUTS, self.p_sat_r, (T_dew + self.T_in_r) / 2)
        cp_r                    = AS.cpmass()       # [J/kg-K]
        rho_superheat           = AS.rhomass()      # [kg/m^3]

        # Compute Fins Efficiency based on FinsType
        if self.FinsType == 'WavyLouveredFins':
            WavyLouveredFins(self.Fins)
        elif self.FinsType == 'HerringboneFins':
            self.height_cond, self.secTheta, self.N_fin, self.t_f  = HerringboneFins(self.Fins)
        elif self.FinsType == 'PlainFins':
            PlainFins(self.Fins)

        self.m_dot_da           = self.Fins.m_dot_da

        # todo start here tomorrow
        # Cross-flow in the superheated region.
        # Using effectiveness-Ntu relationships for cross flow with non-zero Cr.
        UA_overall              = 1 / (1 / (self.Fins.eta_a * self.Fins.h_a * self.Fins.A_a * self.h_a_tuning)
                                       + 1 / (self.h_r_superheat * self.A_r_wetted) + self.R_w)

        epsilon_superheat       = (T_dew - self.T_in_r) / (self.T_in_a - self.T_in_r)

        Ntu                     = UA_overall / (self.m_dot_da * self.Fins.cp_da)
        if epsilon_superheat > 1.0:
            epsilon_superheat   = 1.0-1e-12
        self.w_superheat        = -log(1 - epsilon_superheat) * self.m_dot_r * cp_r / ((1 - exp(-Ntu)) * self.m_dot_da * self.Fins.cp_da)

        # Positive Q is heat input to the refrigerant, negative Q is heat utpout from refrigerant.
        # Heat is removed here from the refrigerant since it is being cooled
        self.Q_superheat        = self.m_dot_r * cp_r * (T_dew - self.T_in_r)
        self.DT_sh              = -self.Q_superheat / (self.m_dot_r * cp_r)     # temperature difference between the inlet and outlet of the SH region
        self.T_out_r            = self.T_in_r - self.DT_sh                      # the outlet temperature of refrigerant of the SH region, should be T_dew

        # Pressure drop calculations for superheated refrigerant
        v_r                     = 1. / rho_superheat
        # Pressure gradient using Darcy friction factor
        dpdz_r                  = -self.f_r_superheat * v_r * self.G_r**2 / (2 * self.ID)   # Pressure gradient
        self.DP_r_superheat     = dpdz_r * (self.L_circuit * self.w_superheat)
        self.Charge_superheat   = self.w_superheat * (self.V_r * rho_superheat)             # total mass of refrigerant in the SH region

        # Latent heat needed for pseudo-quality calc
        AS.update(CP.QT_INPUTS, 0.0, T_bubble)
        h_l                     = AS.hmass()                                                # [J/kg]
        AS.update(CP.QT_INPUTS, 1.0, T_dew)
        h_v                     = AS.hmass()                                                # [J/kg]
        h_fg                    = h_v - h_l                                                 # [J/kg]
        self.xin_r              = 1.0 + cp_r * (self.T_in_r - T_dew) / h_fg

        return self.w_superheat - 1

    def _TwoPhase_Forward(self, x_out_r_2phase=0.0):
        """
            x_out_r_2phase: quality of refrigerant at end of two-phase portion
                default value is 0.0 (full two phase region)
        """
        # Bubble and dew temperatures (same for fluids without glide)
        T_bubble                    = self.T_bubble
        T_dew                       = self.T_dew

        # Mean temperature for use in HT relationships
        T_sat_r                     = (T_bubble + T_dew) / 2

        h_l                         = self.h_l                  # [J/kg]
        h_v                         = self.h_v                  # [J/kg]
        h_fg                        = h_v - h_l                 # [J/kg]

        # This block calculates the average refrigerant heat transfer coefficient by
        # integrating the local heat transfer coefficient between
        # a quality of 1.0 and the outlet quality
        h_r_2phase                  = ShahCondensation_Average(x_out_r_2phase, 1.0, self.AS, self.G_r, self.ID, self.p_sat_r, T_bubble, T_dew)
        self.h_r_2phase             = h_r_2phase * self.h_tp_tuning
        # todo check as UA considered all areas
        UA_overall          = 1 / (1 / (self.Fins.eta_a * self.Fins.h_a * self.Fins.A_a * self.h_a_tuning) + 1 / (self.h_r_2phase * self.A_r_wetted) + self.R_w)
        self.epsilon_2phase = 1 - exp(-UA_overall / (self.m_dot_da * self.Fins.cp_da))

        self.w_2phase       = -self.m_dot_r * h_fg * (1.0 - x_out_r_2phase) / (self.m_dot_da * self.Fins.cp_da * (self.T_in_a - T_sat_r) * self.epsilon_2phase)

        # Positive Q is heat input to the refrigerant, negative Q is heat output from refrigerant.
        # Heat is removed here from the refrigerant since it is condensing
        # todo check the air temperature does not change for each bank?, or circuit
        self.Q_2phase               = self.epsilon_2phase * self.Fins.cp_da * self.m_dot_da * self.w_2phase * (self.T_in_a - T_sat_r)
        self.x_out_2phase           = x_out_r_2phase

        # Frictional pressure drop component
        DP_frict            = LMPressureGradientAvg(self.x_out_2phase, 1.0, self.AS, self.G_r, self.ID, T_bubble, T_dew) * self.L_circuit * self.w_2phase
        # Accelerational pressure drop component
        DP_accel            = -AccelPressureDrop(self.x_out_2phase, 1.0, self.AS, self.G_r, T_bubble, T_dew, slipModel='Zivi') * self.L_circuit * self.w_2phase
        # Total pressure drop is the sum of accelerational and frictional components (neglecting gravitational effects)
        self.DP_r_2phase    = DP_frict + DP_accel

        rho_average         = TwoPhaseDensity(self.AS, self.x_out_2phase, 1.0, self.T_dew, self.T_bubble, slipModel='Zivi')
        self.Charge_2phase  = rho_average * self.w_2phase * self.V_r

        if self.Verbosity > 7:
            print('2phase cond resid', self.w_2phase - (1 - self.w_superheat))
            print('h_r_2phase', self.h_r_2phase)

        # Calculate an effective pseudo-subcooling based on the equality
        #     cp*DT_sc=-dx*h_fg
        cp_sat_L                    = self.cp_sat_L                     # [J/kg-K]
        self.DT_sc_2phase           = -self.x_out_2phase * h_fg / cp_sat_L  # 0 as x_out_2phase = 0

        # If the quality is being solved for, the length of the two-phase and subcooled
        # sections should add to the length of the HX.  Return the residual

        return self.w_2phase-(1-self.w_superheat)

    def _Subcool_Forward(self):
        self.w_subcool              = 1 - self.w_2phase - self.w_superheat

        if self.w_subcool < 0:
            raise ValueError('w_subcool in Condenser cannot be less than zero')

        # AbstractState
        AS                          = self.AS
        # Bubble and dew temperatures (same for fluids without glide)
        T_bubble                    = self.T_bubble

        # Based on the construction of the cycle model there is guaranteed to be a
        # two-phase portion of the heat exchanger
        A_a_subcool                 = self.Fins.A_a * self.w_subcool
        m_dot_da_subcool            = self.m_dot_da * self.w_subcool
        A_r_subcool                 = self.A_r_wetted * self.w_subcool

        # Friction factor and HTC in the refrigerant portions.
        # Average fluid temps are used for the calculation of properties
        # Average temp of refrigerant is average of sat. temp and outlet temp
        # Secondary fluid is air over the fins

        self.f_r_subcool, self.h_r_subcool, self.Re_r_subcool = f_h_1phase_Tube(
          self.m_dot_r / self.N_circuits, self.ID, T_bubble-1.0, self.p_sat_r, self.AS)

        AS.update(CP.PT_INPUTS, self.p_sat_r, T_bubble-1)
        cp_r                        = AS.cpmass()               # [J/kg-K]

        # Cross-flow in the subcooled region.
        R_a                         = 1. / (self.Fins.eta_a * self.Fins.h_a * self.Fins.A_a * self.h_a_tuning)
        R_r                         = 1. / (self.h_r_subcool * self.A_r_wetted)
        UA_subcool                  = self.w_subcool / (R_a + R_r + self.R_w)
        C_min                       = min([self.m_dot_da*self.Fins.cp_da*self.w_subcool, self.m_dot_r*cp_r])
        C_max                       = max([self.m_dot_da*self.Fins.cp_da*self.w_subcool, self.m_dot_r*cp_r])
        Cr                          = C_min / C_max
        NTU                         = UA_subcool / C_min

        if(self.m_dot_da*self.Fins.cp_da*self.w_subcool > self.m_dot_r*cp_r):
            # Minimum capacitance rate on refrigerant side
            epsilon_subcool = 1. - exp(-1. / Cr * (1. - exp(-Cr * NTU)))
        else:
            # Minimum capacitance rate on air side
            epsilon_subcool = 1 / Cr * (1 - exp(-Cr * (1 - exp(-NTU))))

        # Positive Q is heat input to the refrigerant, negative Q is heat output from refrigerant.
        # Heat is removed here from the refrigerant since it is condensing
        self.Q_subcool              = -epsilon_subcool * C_min * (T_bubble - self.T_in_a)
        self.DT_sc                  = -self.Q_subcool / (self.m_dot_r*cp_r)
        self.T_out_r                = T_bubble - self.DT_sc

        AS.update(CP.PT_INPUTS, self.p_sat_r, (T_bubble + self.T_out_r) / 2)
        rho_subcool                 = AS.rhomass()                  # [kg/m^3]
        self.Charge_subcool         = self.w_subcool * self.V_r * rho_subcool

        # Pressure drop calculations for subcooled refrigerant
        v_r                         = 1 / rho_subcool
        # Pressure gradient using Darcy friction factor
        dpdz_r                      = -self.f_r_subcool * v_r * self.G_r**2 / (2 * self.ID)  # Pressure gradient
        self.DP_r_subcool           = dpdz_r * self.L_circuit * self.w_subcool


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
        d                                   = 2e-3
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

        L_tube                              = L_surf        # m, tube length
        Number_of_tubes                     = 15 * 10       # number of tubes in parallel
        ID_tube                             = d_H           # m, tube inner diameter

        # mass flow rate of refrigerant in single tube
        mfr_channel                         = self.m_dot_refrigerant / Number_of_tubes
        u_refrigerant                       = 4 * mfr_channel / (A_c_chan * self.rho_inlet)    # refrigerant flow velocity

        # ------------------------------------------------------------------------------
        # Calculate heat transfer coefficient
        # ------------------------------------------------------------------------------
        # HTC_liq                             = 600           # W/m2K, heat transfer coefficient: refrigerant liquid phase
        HTC_vap                             = 600           # W/m2K, heat transfer coefficient: refrigerant vapor phase
        # f, h, Re                             = f_h_1phase_Channel(mfr_channel, d_H, u_refrigerant, T_sat_vap, self.P_inlet, Phase='Single')
        HTC_2ph                             = 3000          # W/m2K, heat transfer coefficient: refrigerant two phases
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
        self.T_bat_updated              = self.T_bat + (self.Q_bat - self.Q_evaporator) * self.dt / C_bat
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
