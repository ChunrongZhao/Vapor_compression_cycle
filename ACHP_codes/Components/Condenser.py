from __future__                                         import division, print_function, absolute_import
from math                                               import pi,log,exp
from scipy.optimize                                     import brentq
import CoolProp                                         as CP

from ACHP_codes.Correlations.HTC_Correlations           import f_h_1phase_Tube, ShahCondensation_Average
from ACHP_codes.Correlations.deltaP_Correlations        import LMPressureGradientAvg, AccelPressureDrop
from ACHP_codes.Correlations.PhaseState_Correlations    import TwoPhaseDensity

from ACHP_codes.Correlations.FinStructure_Correlations  import WavyLouveredFins, FinInputs, HerringboneFins, PlainFins, IsFinsClass, \
                                                    HerringboneFins_Dimensions, FinInputs_Iter, IsFinsClass_Iter, HerringboneFins_Iter, \
                                                    HerringboneFins_Dimensions_Iter
from ACHP_codes.ACHP_Tools.ACHP_Tools                   import ValidateFields


# ------------------------------------------------------------------------
# one bank layer
# ------------------------------------------------------------------------
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
        # Only validate the first time
        if not hasattr(self, 'IsValidated'):
            self.Fins.Validate()
            reqFields = [
               ('Fins',         IsFinsClass,    None,   None),
               ('FinsType',     str,            None,   None),
               ('m_dot_r',      float,          0.00001,20),
               ('T_in_r',       float,          200,    500),
               ('p_sat_r',      float,          0.01,   20000000)]
            optFields               = ['Verbosity', 'AS', 'h_a_tuning', 'h_tp_tuning', 'DP_tuning']
            ValidateFields(self.__dict__, reqFields, optFields)
            self.IsValidated        = True

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
        rho_pipe            = 2710          # aluminium
        rho_fin             = 8960          # copper
        self.L_COND, self.H_COND, self.W_COND, self.m_COND, self.V_COND = HerringboneFins_Dimensions(self.height_cond, self.L_tube, self.N_bank,
                        self.P_l, self.secTheta, self.t_f, self.N_fin, rho_fin, self.ID, self.OD, self.N_Tubes_per_bank, rho_pipe, self.Charge)

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


def SampleCondenser(AS, T=41.37):
    Fins               = FinInputs()
    Fins.Tubes.N_Tubes_per_bank     = 45                    # number of tubes per bank or row
    Fins.Tubes.N_bank               = 1                     # number of banks or rows
    Fins.Tubes.N_circuits           = 2                     # number of circuits
    Fins.Tubes.L_tube               = 1.7                   # one tube length
    Fins.Tubes.OD                   = 0.007
    Fins.Tubes.ID                   = 0.0063904
    Fins.Tubes.Pl                   = 0.0191                # distance between center of tubes in flow direction
    Fins.Tubes.Pt                   = 0.0222                # distance between center of tubes orthogonal to flow direction
    Fins.Tubes.k_w                  = 237                   # Wall thermal conductivity

    Fins.Fins.FPI                   = 25                    # Number of fins per inch
    Fins.Fins.Pd                    = 0.001                 # 2* amplitude of wavy fin
    Fins.Fins.xf                    = 0.001                 # 1/2 period of fin
    Fins.Fins.t                     = 0.00011               # Thickness of fin material
    Fins.Fins.k_fin                 = 237                   # Thermal conductivity of fin material

    Fins.Air.V_dot_ha               = 1.7934                # rated volumetric flow rate
    Fins.Air.T_mean                 = 308.15
    # Fins.Air.T_db                   = 308.15                # Dry Bulb Temperature
    Fins.Air.T_db                   = 34.638755258533024 + 273.15
    Fins.Air.p                      = 101325                # Air pressure in Pa
    Fins.Air.RH                     = 0.51                  # Relative Humidity
    Fins.Air.RH_mean                = 0.51
    Fins.Air.FanPower               = 160

    params = {
        'AS':                       AS,                     # Abstract State
        'm_dot_r':                  0.0708,
        'T_in_r':                   T + 20 + 273.15,
        'p_sat_r':                  2500076.19,             # PropsSI('P','T',T+273.15,'Q',1.0,'R410A')
        'Fins':                     Fins,
        'FinsType':                 'HerringboneFins',      # Choose fin Type: 'WavyLouveredFins' or 'HerringboneFins' or 'PlainFins'
        'Verbosity':                0,
        'h_a_tuning':               1,
        'h_tp_tuning':              1,
        'DP_tuning':                1
    }

    Cond                            = CondenserClass(**params)

    Cond.Calculate()

    return Cond


def SampleIntercooler(AS,T=41.37):
    # Only superheated region
    Fins=FinInputs()
    Fins.Tubes.NTubes_per_bank=41       #number of tubes per bank or row
    Fins.Tubes.Nbank=1                  #number of banks or rows
    Fins.Tubes.Ncircuits=5              #number of circuits
    Fins.Tubes.Ltube=3                  #one tube length
    Fins.Tubes.OD=0.007
    Fins.Tubes.ID=0.0063904
    Fins.Tubes.Pl=0.0191                #distance between center of tubes in flow direction
    Fins.Tubes.Pt=0.0222                #distance between center of tubes orthogonal to flow direction
    Fins.Tubes.kw=237                   #Wall thermal conductivity

    Fins.Fins.FPI=25                    #Number of fins per inch
    Fins.Fins.Pd=0.001                  #2* amplitude of wavy fin
    Fins.Fins.xf=0.001                  #1/2 period of fin
    Fins.Fins.t=0.00011                 #Thickness of fin material
    Fins.Fins.k_fin=237                 #Thermal conductivity of fin material

    Fins.Air.Vdot_ha=8.998             #rated volumetric flowrate
    Fins.Air.Tmean=48.9+273.15
    Fins.Air.Tdb=48.9+273.15                 #Dry Bulb Temperature
    Fins.Air.p=101325                   #Air pressure in Pa
    Fins.Air.RH=0.51                    #Relative Humidity
    Fins.Air.RHmean=0.51
    Fins.Air.FanPower=160

    params={
        'AS': AS, #Abstract State
        'mdot_r': 1.209,
        'Tin_r': T+273.15,
        'psat_r': 716700,
        'Fins': Fins,
        'FinsType': 'HerringboneFins',  #Choose fin Type: 'WavyLouveredFins' or 'HerringboneFins'or 'PlainFins'
        'Verbosity':0,
        'h_a_tuning':1,
        'h_tp_tuning':1,
        'DP_tuning':1,
    }
    Cond=CondenserClass(**params)
    Cond.Calculate()
    return Cond


class FinVals():
    def __init__(self):
        pass


# ------------------------------------------------------------------------
# New Testing with multi-layer bank
# ------------------------------------------------------------------------
class CondenserClass_Iterate():
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
            # ('Inlet Dry bulb temp',                 'K',            self.T_in_a),
            ('Inlet Air pressure',                  'Pa',           self.Fins.Air.p),
            ('Inlet Air Relative Humidity',         '-',            self.Fins.Air.RH),
            ('Tubes per bank',                      '-',            self.Fins.Tubes.N_Tubes_per_bank),
            # ('Number of banks',                     '-',            self.Fins.Tubes.N_bank),
            # ('Number circuits',                     '-',            self.Fins.Tubes.N_circuits),
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

    def Calculate(self, T_in_a):
        # Only validate the first time
        if not hasattr(self, 'IsValidated'):
            self.Fins.Validate()
            reqFields = [
               ('Fins',         IsFinsClass_Iter,    None,   None),
               ('FinsType',     str,            None,   None),
               ('m_dot_r',      float,          0.00001,20),
               ('T_in_r',       float,          200,    500),
               ('p_sat_r',      float,          0.01,   20000000)]
            optFields               = ['Verbosity', 'AS', 'h_a_tuning', 'h_tp_tuning', 'DP_tuning']
            ValidateFields(self.__dict__, reqFields, optFields)
            self.IsValidated        = True

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
        self.N_circuits                         = self.Fins.Tubes.N_circuits_per_bank
        self.P_l                                = self.Fins.Tubes.Pl
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
        TotalLength                             = self.L_tube * self.N_Tubes_per_bank
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
        self._Superheat_Forward(T_sh_out_r=self.T_dew, T_in_a=T_in_a)
        # Maybe have a full two-phase section
        # First try to run with a full two-phase section from quality of 1 to quality of 0
        self._TwoPhase_Forward(T_in_a=T_in_a)
        # If we have already used too much of the HX (max possible sum of w is 1.0)
        if self.w_superheat >= 1:
            # There is partial superheated region
            self.T_in_a                         = T_in_a
            brentq(self._Superheat_Forward_Iter, self.T_dew, self.T_in_r)
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
            brentq(self._TwoPhase_Forward(T_in_a=T_in_a), 0.0000001, 0.9999999)
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
            self._Subcool_Forward(T_in_a=T_in_a)

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
        T_out_a             = T_in_a - self.Q / (self.Fins.cp_da * self.Fins.m_dot_da)
        self.h_mean_r       = self.w_2phase * self.h_r_2phase + self.w_superheat * self.h_r_superheat + self.w_subcool * self.h_r_subcool
        self.UA_r           = self.h_mean_r * self.A_r_wetted
        self.UA_a           = (self.Fins.h_a * self.h_a_tuning) * self.Fins.A_a * self.Fins.eta_a
        self.UA_w           = 1 / self.R_w

        # added by Chunrong at 24/11/2023
        rho_pipe            = 2710          # aluminium
        rho_fin             = 8960          # copper
        self.L_COND, self.H_COND, self.W_COND, self.m_COND, self.V_COND = HerringboneFins_Dimensions_Iter(self.height_cond, self.L_tube,
                        self.P_l, self.secTheta, self.t_f, self.N_fin, rho_fin, self.ID, self.OD, self.N_Tubes_per_bank, rho_pipe, self.Charge)

        return T_out_a

    def _Superheat_Forward(self, T_sh_out_r, T_in_a=None):
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
            self.height_cond, self.secTheta, self.N_fin, self.t_f  = HerringboneFins_Iter(self.Fins, T_in_a)
        elif self.FinsType == 'PlainFins':
            PlainFins(self.Fins)

        self.m_dot_da           = self.Fins.m_dot_da

        # todo start here tomorrow
        # Cross-flow in the superheated region.
        # Using effectiveness-Ntu relationships for cross flow with non-zero Cr.
        UA_overall              = 1 / (1 / (self.Fins.eta_a * self.Fins.h_a * self.Fins.A_a * self.h_a_tuning)
                                       + 1 / (self.h_r_superheat * self.A_r_wetted) + self.R_w)

        epsilon_superheat       = (T_dew - self.T_in_r) / (T_in_a - self.T_in_r)

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

    def _Superheat_Forward_Iter(self, T_sh_out_r):
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
            self.height_cond, self.secTheta, self.N_fin, self.t_f  = HerringboneFins_Iter(self.Fins, self.T_in_a)
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

    def _TwoPhase_Forward(self, x_out_r_2phase=0.0, T_in_a=None):
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

        self.w_2phase       = -self.m_dot_r * h_fg * (1.0 - x_out_r_2phase) / (self.m_dot_da * self.Fins.cp_da * (T_in_a - T_sat_r) * self.epsilon_2phase)

        # Positive Q is heat input to the refrigerant, negative Q is heat output from refrigerant.
        # Heat is removed here from the refrigerant since it is condensing
        # todo check the air temperature does not change for each bank?, or circuit
        self.Q_2phase               = self.epsilon_2phase * self.Fins.cp_da * self.m_dot_da * self.w_2phase * (T_in_a - T_sat_r)
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

    def _Subcool_Forward(self, T_in_a=None):
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
        self.Q_subcool              = -epsilon_subcool * C_min * (T_bubble - T_in_a)
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


def SampleCondenser_Iter(AS, T=41.37):
    Fins               = FinInputs_Iter()
    Fins.Tubes.N_Tubes_per_bank     = 45                    # number of tubes per bank or row
    # Fins.Tubes.N_bank               = 1                     # number of banks or rows
    Fins.Tubes.N_circuits_per_bank  = 2                     # number of circuits
    Fins.Tubes.L_tube               = 1.7                   # one tube length
    Fins.Tubes.OD                   = 0.007
    Fins.Tubes.ID                   = 0.0063904
    Fins.Tubes.Pl                   = 0.0191                # distance between center of tubes in flow direction
    Fins.Tubes.Pt                   = 0.0222                # distance between center of tubes orthogonal to flow direction
    Fins.Tubes.k_w                  = 237                   # Wall thermal conductivity

    Fins.Fins.FPI                   = 25                    # Number of fins per inch
    Fins.Fins.Pd                    = 0.001                 # 2* amplitude of wavy fin
    Fins.Fins.xf                    = 0.001                 # 1/2 period of fin
    Fins.Fins.t                     = 0.00011               # Thickness of fin material
    Fins.Fins.k_fin                 = 237                   # Thermal conductivity of fin material

    Fins.Air.V_dot_ha               = 1.7934                # rated volumetric flow rate
    # Fins.Air.T_mean                 = 308.15
    # Fins.Air.T_db                   = 308.15                # Dry Bulb Temperature
    Fins.Air.p                      = 101325                # Air pressure in Pa
    Fins.Air.RH                     = 0.51                  # Relative Humidity
    Fins.Air.RH_mean                = 0.51
    Fins.Air.FanPower               = 160

    params = {
        'AS':                       AS,                     # Abstract State
        'm_dot_r':                  0.0708,
        'T_in_r':                   T + 20 + 273.15,
        'p_sat_r':                  2500076.19,             # PropsSI('P','T',T+273.15,'Q',1.0,'R410A')
        'Fins':                     Fins,
        'FinsType':                 'HerringboneFins',      # Choose fin Type: 'WavyLouveredFins' or 'HerringboneFins' or 'PlainFins'
        'Verbosity':                0,
        'h_a_tuning':               1,
        'h_tp_tuning':              1,
        'DP_tuning':                1
    }

    Cond                            = CondenserClass_Iterate(**params)

    # Number of banks
    N_bank                          = 3
    T_in_a                          = 20. + 273.15
    for i in range(1, N_bank+1):
        T_out_a                     = Cond.Calculate(T_in_a)
        print("Bank number = ", i, "; inlet and outlet air temperatures are {} oC".format([T_in_a-273.15, T_out_a-273.15]))
        T_in_a                      = T_out_a
        print('Heat transfer rate in condenser is',                     Cond.Q,             'W')
        print('Heat transfer rate in condenser (superheat section) is', Cond.Q_superheat,   'W')
        print('Heat transfer rate in condenser (two-phase section) is', Cond.Q_2phase,      'W')
        print('Heat transfer rate in condenser (subcooled section) is', Cond.Q_subcool,     'W')
        print('Fraction of circuit length in superheated section is',   Cond.w_superheat)
        print('Fraction of circuit length in two-phase section is',     Cond.w_2phase)
        print('Fraction of circuit length in subcooled section is',     Cond.w_subcool)
        print('***************************************************************')
        print('Condenser [Height, Length, and Width] are {} m, respectively'.format([Cond.H_COND, Cond.L_COND, Cond.W_COND]))
        print('Condenser weight is {} kg'.format(Cond.m_COND))

    m_COND_total        = N_bank * Cond.m_COND
    W_COND_total        = N_bank * Cond.W_COND

    print('total weight and width of the condenser are {}'.format([m_COND_total, W_COND_total]))

    return Cond


# -----------------------------------------------------------------------------------------
if __name__=='__main__':
    # This runs if you run this file directly
    Ref         = 'R410A'
    Backend     = 'HEOS'  # choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    AS          = CP.AbstractState(Backend, Ref)    # Abstract State
    Cond        = SampleCondenser_Iter(AS, 43.3)
    # Cond        = SampleCondenser(AS, 43.3)
    print(Cond.OutputList())


