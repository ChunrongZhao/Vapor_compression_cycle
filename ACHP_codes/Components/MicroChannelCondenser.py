from __future__                                         import division, print_function, absolute_import
from math                                               import pi, log, exp
from CoolProp.CoolProp                                  import HAPropsSI
from ACHP_codes.Correlations.HTC_Correlations           import f_h_1phase_MicroTube, KM_Cond_Average
from ACHP_codes.Correlations.PhaseState_Correlations    import TwoPhaseDensity
from ACHP_codes.Correlations.deltaP_Correlations        import AccelPressureDrop
from ACHP_codes.Correlations.MicroFinCorrelations       import MultiLouveredMicroFins, MicroFinInputs, IsFinsClass
from scipy.optimize                                     import brentq, fsolve
from ACHP_codes.ACHP_Tools.ACHP_Tools                   import ValidateFields
import CoolProp                                         as CP


# -----------------------------------------------------------------------------------------------------------
class FinVals():
    def __init__(self):
        pass


class MicroCondenserClass():
    def __init__(self, **kwargs):
        # Load the parameters passed in
        # using the dictionary
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
            ('Volumetric flow rate of humid air',           'm^3/s',            self.Fins.Air.V_dot_ha),
            ('Inlet Dry bulb temp',                         'K',                self.T_in_a),
            ('Inlet Air pressure',                          'Pa',               self.Fins.Air.p),
            ('Inlet Air Relative Humidity',                 '-',                self.Fins.Air.RH),
            ('Tubes per bank',                              '-',                self.Fins.Tubes.N_Tubes),
            ('Number of banks',                             '-',                self.Fins.Tubes.N_bank),
            ('Number of passes',                            '-',                self.Fins.Tubes.N_pass),
            ('Number of ports',                             '-',                self.Fins.Tubes.N_ports),
            ('Length of tube',                              'm',                self.Fins.Tubes.L_tube),
            ('Tube width',                                  'm',                self.Td),
            ('Tube height',                                 'm',                self.Ht),
            ('Tube spacing',                                'm',                self.Fins.Tubes.b),
            ('Tube thickness',                              'm',                self.Fins.Tubes.tw),
            ('Wall port thickness',                         'm',                self.Fins.Tubes.twp),
            ('Channel aspect ratio',                        '-',                self.Fins.Tubes.beta),
            ('Tube Conductivity',                           'W/m-K',            self.Fins.Tubes.k_w),
            ('Fins per inch',                               '1/in',             self.Fins.Fins.FPI),
            ('Fin length',                                  'm',                self.Fins.Fins.Lf),
            ('Fin thickness',                               'm',                self.Fins.Fins.t),
            ('Fin Conductivity',                            'W/m-K',            self.Fins.Fins.k_fin),
            ('Louver angle',                                'degree',           self.Fins.Louvers.Lalpha),
            ('Louver pitch',                                'm',                self.Fins.Louvers.lp),
            ('Louver cut length',                           'm',                self.Fins.Llouv),
            ('Fins Type',                                   '-',                self.FinsType),
            ('Q Total',                                     'W',                self.Q),
            ('Q Superheat',                                 'W',                self.Q_superheat),
            ('Q Two-Phase',                                 'W',                self.Q_2phase),
            ('Q Subcool',                                   'W',                self.Q_subcool),
            ('Inlet Temp',                                  'K',                self.T_in_r),
            ('Outlet Temp',                                 'K',                self.T_out_r),
            ('Pressure Drop Total',                         'Pa',               self.DP_r),
            ('Pressure Drop Superheat',                     'Pa',               self.DP_r_superheat),
            ('Pressure Drop Two-Phase',                     'Pa',               self.DP_r_2phase),
            ('Pressure Drop Subcool',                       'Pa',               self.DP_r_subcool),
            ('Charge Total',                                'kg',               self.Charge),
            ('Charge Superheat',                            'kg',               self.Charge_superheat),
            ('Charge Two-Phase',                            'kg',               self.Charge_2phase),
            ('Charge Subcool',                              'kg',               self.Charge_subcool),
            ('Mean HTC Superheat',                          'W/m^2-K',          self.h_r_superheat),
            ('Mean HTC Two-phase',                          'W/m^2-K',          self.h_r_2phase),
            ('Mean HTC Subcool',                            'W/m^2-K',          self.h_r_subcool),
            ('Wetted Area Fraction Superheat',              '-',                self.w_superheat),
            ('Wetted Area Fraction Two-phase',              '-',                self.w_2phase),
            ('Wetted Area Fraction Subcool',                '-',                self.w_subcool),
            ('Mean Air HTC',                                'W/m^2-K',          self.Fins.h_a * self.h_a_tuning),
            ('Surface Effectiveness',                       '-',                self.Fins.eta_a),
            ('Air-side area (fin+tubes)',                   'm^2',              self.Fins.A_a),
            ('Mass Flow rate of Dry Air',                   'kg/s',             self.Fins.m_dot_da),
            ('Mass Flow rate of Humid Air',                 'kg/s',             self.Fins.m_dot_ha),
            ('Pressure Drop Air-side (core only)',          'Pa',               self.Fins.dP_a),
            ('Pressure Drop Air-side (total)',              'Pa',               self.dP_a),
            ('Subcooling',                                  'K',                self.DT_sc),
            ('Number of Circuits',                          '-',                self.N_circuits)
        ]

    def Update(self, **kwargs):
        # Update the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)

    def Calculate(self):
        # Only validate the first time
        if not hasattr(self, 'IsValidated'):
            self.Fins.Validate()
            reqFields       = [
               ('Fins',             IsFinsClass,        None,           None),
               ('FinsType',         str,                None,           None),
               ('m_dot_r',          float,              0.00001,        20),
               ('T_in_r',           float,              200,            500),
               ('p_sat_r',          float,              0.01,           20000000)
               ]
            optFields       = ['Verbosity', 'AS', 'h_a_tuning', 'h_tp_tuning', 'DP_tuning']
            ValidateFields(self.__dict__, reqFields, optFields)
            self.IsValidated                            = True

        # set tuning factors to 1 in case not given by user
        if not hasattr(self, 'h_a_tuning'):
            self.h_a_tuning                             = 1
        if not hasattr(self, 'h_tp_tuning'):
            self.h_tp_tuning                            = 1
        if not hasattr(self, 'DP_tuning'):
            self.DP_tuning                              = 1

        # AbstractState
        AS                                              = self.AS

        # Retrieve some parameters from nested structures
        # for code compactness
        self.L_tube                                     = self.Fins.Tubes.L_tube            # tube length
        self.N_Tubes                                    = self.Fins.Tubes.N_Tubes           # number of tube (per bank)
        self.N_bank                                     = self.Fins.Tubes.N_bank            # number of banks
        self.T_in_a                                     = self.Fins.Air.T_db                # inlet air temperature
        self.P_in_a                                     = self.Fins.Air.p                   # inlet air pressure
        self.RH_in_a                                    = self.Fins.Air.RH                  # inlet air relative humidity
        self.Td                                         = self.Fins.Tubes.Td                # tube outside width (depth)
        self.Ht                                         = self.Fins.Tubes.Ht                # tube outside height (major diameter)
        self.b                                          = self.Fins.Tubes.b                 # tube spacing
        self.tw                                         = self.Fins.Tubes.tw                # tube thickness
        self.N_pass                                     = self.Fins.Tubes.N_pass            # Number of passes on ref-side (per bank)
        self.k_w                                        = self.Fins.Tubes.k_w               # thermal conductivity of tube wall
        self.N_ports                                    = self.Fins.Tubes.N_ports           # Number of rectangular ports
        self.twp                                        = self.Fins.Tubes.twp               # Port wall thickness
        self.beta                                       = self.Fins.Tubes.beta              # channel (port) aspect ratio (=width/height)

        # Bubble and dew temperatures (same for fluids without glide)
        AS.update(CP.PQ_INPUTS, self.p_sat_r, 0.0)
        self.T_bubble                                   = AS.T()                            # [K]
        self.h_l                                        = AS.hmass()                        # [J/kg]
        self.cp_sat_L                                   = AS.cpmass()                       # [J/kg-K]

        AS.update(CP.PQ_INPUTS, self.p_sat_r, 1.0)
        self.T_dew                                      = AS.T()                            # [K]
        self.h_v                                        = AS.hmass()                        # [J/kg]

        # Define Number of circuits (=number of tubes per pass)
        self.N_circuits                                 = self.N_Tubes / self.N_pass
        # Calculate an effective length of circuit if circuits are
        # not all the same length
        TotalLength                                     = self.L_tube * self.N_Tubes * self.N_bank
        self.L_circuit                                  = TotalLength / self.N_circuits

        # Volume of refrigerant = rectangle of tube + circular part at the ends - thickness between ports
        self.V_r            = ((self.Td - self.Ht) * (self.Ht - 2.0 * self.tw) + (pi / 4.0) * (self.Ht - 2.0 * self.tw)**2
                               - (self.Ht - 2.0 * self.tw) * self.twp * (self.N_ports - 1)) * self.L_circuit * self.N_circuits
        # Tube wetted area = tube straight length + circular shape at the ends - horizontal port thickness  + vertical thickness between ports
        self.A_r_wetted     = (2.0 * (self.Td - self.Ht) + pi * (self.Ht - 2.0 * self.tw) - 2.0 * self.twp * (self.N_ports - 1)
                               + 2.0 * (self.Ht - 2.0 * self.tw) * (self.N_ports - 1)) * self.L_circuit * self.N_circuits
        # Free-flow area on refrigerant-side = area of rectangle tube + circular parts at end - area thickness between ports
        self.A_c            = ((self.Td - self.Ht) * (self.Ht - 2.0 * self.tw) + (pi / 4.0) * (self.Ht - 2.0 * self.tw)**2
                               - self.twp * (self.Ht - 2.0 * self.tw) * (self.N_ports - 1)) * self.N_circuits
        # Hydraulic diameter on ref-side
        self.Dh             = 4 * self.A_c * self.L_circuit / self.A_r_wetted
        # Mass flux ref-side
        self.G_r            = self.m_dot_r / self.A_c

        # Total conduction area (exclude port's thickness)
        self.Aw             = 2 * (self.Td - self.twp * (self.N_ports - 1)) * self.L_circuit * self.N_circuits
        # Thermal resistance at the wall
        self.Rw             = self.tw / (self.k_w * self.Aw)

        # Definitely have a superheated portion
        self._Superheat_Forward()
        # Maybe have a full two-phase section
        # First try to run with a full two-phase section from quality of 1 to quality of 0
        self._TwoPhase_Forward()
        # If we have already used too much of the HX (max possible sum of w is 1.0)
        if self.w_2phase + self.w_superheat > 1:
            # There is no subcooled portion, solve for outlet quality
            brentq(self._TwoPhase_Forward, 0.0000001, 0.9999999)
            # Zero out all the subcooled parameters
            self.Q_subcool                              = 0.0
            self.DP_r_subcool                           = 0.0
            self.Charge_subcool                         = 0.0
            self.w_subcool                              = 0.0
            self.h_r_subcool                            = 0.0
            self.existsSubcooled                        = False
        else:
            # By definition then we have a subcooled portion, solve for it
            self.existsSubcooled                        = True
            self._Subcool_Forward()

        # Overall calculations
        self.Q                                          = self.Q_superheat + self.Q_2phase + self.Q_subcool
        self.DP_r                                       = self.DP_r_superheat + self.DP_r_2phase + self.DP_r_subcool
        self.DP_r                                       = self.DP_r * self.DP_tuning            # correcting the pressure drop
        self.Charge                                     = self.Charge_2phase + self.Charge_subcool + self.Charge_superheat

        # define known parameters
        AS.update(CP.PT_INPUTS, self.p_sat_r, self.T_in_r)
        self.h_in_r                                     = AS.hmass()                            # [J/kg]
        self.s_in_r                                     = AS.smass()                            # [J/kg-K]

        if self.existsSubcooled == True:
            AS.update(CP.PT_INPUTS, self.p_sat_r, self.T_out_r)
            self.h_out_r                                = AS.hmass()                            # [J/kg]
            self.s_out_r                                = AS.smass()                            # [J/kg-K]
        else:
            self.T_out_r                                = self.x_out_2phase * self.T_dew + (1 - self.x_out_2phase) * self.T_bubble
            AS.update(CP.QT_INPUTS, 0.0, self.T_bubble)
            h_l                                         = AS.hmass()                            # [J/kg]
            s_l                                         = AS.smass()                            # [J/kg-K]
            AS.update(CP.QT_INPUTS, 1.0, self.T_dew)
            h_v                                         = AS.hmass()                            # [J/kg]
            s_v                                         = AS.smass()                            # [J/kg-K]
            self.h_out_r                                = h_l + self.x_out_2phase * (h_v - h_l)
            self.s_out_r                                = s_l + self.x_out_2phase * (s_v - s_l)
            # Use the effective subcooling
            self.DT_sc                                  = self.DT_sc_2phase

        # Calculate the mean outlet air temperature [K]
        self.T_out_a                                    = self.T_in_a - self.Q / (self.Fins.cp_da * self.Fins.m_dot_da)
        self.h_mean_r                                   = self.w_2phase * self.h_r_2phase + self.w_superheat * self.h_r_superheat \
                                                          + self.w_subcool * self.h_r_subcool
        self.UA_r                                       = self.h_mean_r * self.A_r_wetted
        self.UA_a                                       = (self.Fins.h_a * self.h_a_tuning) * self.Fins.A_a * self.Fins.eta_a
        self.UA_w                                       = 1 / self.Rw

        # Update air-side pressure drop based on the outlet air temperature
        # the air-side pressure drop here include momentum, contraction and expansion effects
        # Objective function
        def OBJECTIVE(x):
            Pair_o                                      = x[0]
            W                                           = x[1]
            if W < 0:                                   # to ensure that the humidity ratio is
                print('Microchannel Condenser -- Humidity ratio for air pressure drop is less than zero. Humidity ratio is set to 0.0')
                W                                       = 0.0
            v_da                                        = HAPropsSI('V', 'T', self.T_out_a, 'P', Pair_o, 'W', W)
            W_new                                       = HAPropsSI('W', 'T', self.T_out_a, 'P', Pair_o, 'V', v_da)

            # outlet air density
            rho_o                                       = 1 / v_da * (1 + W_new)                # [m^3/kg_ha]
            # mean air density
            rho_m                                       = pow(0.5 * (1 / self.Fins.rho_i_air + 1/rho_o), -1)
            # air-side pressure drop including momentum, expansion and contraction effects
            DeltaP_air          = self.Fins.G_air**2 / 2 / self.Fins.rho_i_air * ((1 - self.Fins.sigma**2 + self.Fins.Kc_tri)
                                    + 2 * (self.Fins.rho_i_air / rho_o - 1) + self.Fins.f_a * self.Fins.A_a / self.Fins.A_a_c
                                    * (self.Fins.rho_i_air / rho_m) - (1 - self.Fins.sigma**2 - self.Fins.Ke_tri) * (self.Fins.rho_i_air / rho_o))

            resids                                      = [(self.P_in_a - Pair_o) - DeltaP_air, W - W_new]
            return resids

        # Initial guesses
        P_init                                          = self.P_in_a
        w_init                                          = HAPropsSI('W', 'T', self.T_in_a, 'P', self.P_in_a, 'R', self.RH_in_a)
        # solve for outlet air pressure and outlet humidity ratio
        x                                               = fsolve(OBJECTIVE, [P_init, w_init])
        # update the air-side pressure drop
        self.dP_a                                       = self.P_in_a - x[0]

    def _Superheat_Forward(self):
        # **********************************************************************
        #                      SUPERHEATED PART
        # **********************************************************************
        # AbstractState
        AS                                              = self.AS
        # Dew temperature for constant pressure cooling to saturation
        T_dew                                           = self.T_dew
        T_bubble                                        = self.T_bubble

        # Average fluid temps are used for the calculation of properties
        # Average temp of refrigerant is average of sat. temp and outlet temp
        # Secondary fluid is air over the fins
        self.f_r_superheat, self.h_r_superheat, self.Re_r_superheat = \
            f_h_1phase_MicroTube(self.G_r, self.Dh, (T_dew + self.T_in_r) / 2.0, self.p_sat_r, self.AS, "Single")

        AS.update(CP.PT_INPUTS, self.p_sat_r, (T_dew + self.T_in_r) / 2)
        cp_r                                            = AS.cpmass()                       # [J/kg-K]

        # Compute Fins Efficiency based on FinsType
        if self.FinsType == 'MultiLouveredMicroFins':
            MultiLouveredMicroFins(self.Fins)

        self.m_dot_da                                   = self.Fins.m_dot_da

        # Cross-flow in the superheated region.
        # Using effectiveness-Ntu relationships for cross flow with non-zero Cr.
        UA_overall              = 1. / (1. / (self.Fins.eta_a * self.Fins.h_a * self.Fins.A_a * self.h_a_tuning)
                                        + 1. / (self.h_r_superheat * self.A_r_wetted) + self.Rw)
        epsilon_superheat                               = (T_dew - self.T_in_r) / (self.T_in_a - self.T_in_r)
        Ntu                                             = UA_overall / (self.m_dot_da * self.Fins.cp_da)
        if epsilon_superheat > 1.0:
            epsilon_superheat                           = 1.0 - 1e-12
        self.w_superheat        = -log(1 - epsilon_superheat) * self.m_dot_r * cp_r / ((1 - exp(-Ntu)) * self.m_dot_da * self.Fins.cp_da)

        # Positive Q is heat input to the refrigerant, negative Q is heat output from refrigerant.
        # Heat is removed here from the refrigerant since it is being cooled
        self.Q_superheat                                = self.m_dot_r * cp_r * (T_dew - self.T_in_r)

        AS.update(CP.PT_INPUTS, self.p_sat_r, (self.T_in_r + T_dew) / 2.0)
        rho_superheat                                   = AS.rhomass()                      # [kg/m^3]
        # Pressure drop calculations for superheated refrigerant
        v_r                                             = 1. / rho_superheat
        # Pressure gradient using Darcy friction factor
        dpdz_r                                          = -self.f_r_superheat * v_r * self.G_r**2 / (2. * self.Dh)  # Pressure gradient
        self.DP_r_superheat                             = dpdz_r * self.L_circuit * self.w_superheat
        self.Charge_superheat                           = self.w_superheat * self.V_r * rho_superheat

        # Latent heat needed for pseudo-quality calc
        AS.update(CP.QT_INPUTS, 0.0, T_bubble)
        h_l                                             = AS.hmass()                        # [J/kg]
        AS.update(CP.QT_INPUTS, 1.0, T_dew)
        h_v                                             = AS.hmass()                        # [J/kg]
        h_fg                                            = h_v - h_l                         # [J/kg]
        self.x_in_r                                     = 1.0 + cp_r * (self.T_in_r - T_dew) / h_fg

    def _TwoPhase_Forward(self, x_out_r_2phase=0.0):
        # **********************************************************************
        #                      TWO-PHASE PART
        # **********************************************************************
        """
            x_out_r_2phase: quality of refrigerant at end of two-phase portion
                default value is 0.0 (full two phase region)
        """
        # Bubble and dew temperatures (same for fluids without glide)
        T_bubble                                        = self.T_bubble
        T_dew                                           = self.T_dew
        # Mean temperature for use in HT relationships
        T_sat_r                                         = (T_bubble + T_dew) / 2

        h_l                                             = self.h_l                          # [J/kg]
        h_v                                             = self.h_v                          # [J/kg]
        h_fg                                            = h_v - h_l                         # [J/kg]

        # This block calculates the average frictional pressure drop gradient
        # and average refrigerant heat transfer coefficient by
        # integrating the local heat transfer coefficient between
        # a quality of 1.0 and the outlet quality
        DPDZ_frict_2phase, h_r_2phase   \
            = KM_Cond_Average(x_out_r_2phase, 1.0, self.AS, self.G_r, self.Dh, T_bubble, T_dew, self.p_sat_r, self.beta)

        self.h_r_2phase                                 = h_r_2phase * self.h_tp_tuning

        UA_overall                  = 1 / (1 / (self.Fins.eta_a * self.Fins.h_a * self.Fins.A_a * self.h_a_tuning)
                                           + 1 / (self.h_r_2phase * self.A_r_wetted) + self.Rw)
        self.epsilon_2phase         = 1 - exp(-UA_overall / (self.m_dot_da * self.Fins.cp_da))
        self.w_2phase               = -self.m_dot_r * h_fg * (1.0 - x_out_r_2phase) / \
                                      (self.m_dot_da * self.Fins.cp_da * (self.T_in_a - T_sat_r) * self.epsilon_2phase)

        # Positive Q is heat input to the refrigerant, negative Q is heat output from refrigerant.
        # Heat is removed here from the refrigerant since it is condensing
        self.Q_2phase               = self.epsilon_2phase * self.Fins.cp_da * self.m_dot_da * self.w_2phase * (self.T_in_a - T_sat_r)

        self.x_out_2phase                               = x_out_r_2phase

        # Frictional pressure drop component
        DP_frict                    = DPDZ_frict_2phase * self.L_circuit * self.w_2phase
        # Accelerational pressure drop component
        DP_accel                    = -AccelPressureDrop(self.x_out_2phase, 1.0, self.AS, self.G_r, T_bubble, T_dew, slipModel='Zivi') \
                                      * self.L_circuit * self.w_2phase
        # Total pressure drop is the sum of accelerational and frictional components (neglecting gravitational effects)
        self.DP_r_2phase            = DP_frict + DP_accel

        rho_average                 = TwoPhaseDensity(self.AS, self.x_out_2phase, 1.0, self.T_dew, self.T_bubble, slipModel='Zivi')
        self.Charge_2phase          = rho_average * self.w_2phase * self.V_r

        if self.Verbosity > 7:
            print('2phase cond resid', self.w_2phase - (1 - self.w_superheat))
            print('h_r_2phase', self.h_r_2phase)

        # Calculate an effective pseudo-subcooling based on the equality
        cp_sat_L                    = self.cp_sat_L
        self.DT_sc_2phase           = -self.x_out_2phase * h_fg / cp_sat_L

        # If the quality is being solved for, the length of the two-phase and subcooled
        # sections should add to the length of the HX.  Return the residual
        return self.w_2phase - (1 - self.w_superheat)

    def _Subcool_Forward(self):
        # **********************************************************************
        #                      SUBCOOLED PART
        # **********************************************************************
        self.w_subcool              = 1 - self.w_2phase - self.w_superheat

        if self.w_subcool < 0:
            raise ValueError('w_subcool in Condenser cannot be less than zero')
        # AbstractState
        AS                          = self.AS
        # Bubble temperature
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
        self.f_r_subcool, self.h_r_subcool, self.Re_r_subcool = f_h_1phase_MicroTube(
          self.G_r, self.Dh, T_bubble - 1.0, self.p_sat_r, self.AS, "Single")

        AS.update(CP.PT_INPUTS, self.p_sat_r, T_bubble - 1)
        cp_r                        = AS.cpmass()                                           # [J/kg-K]

        # Cross-flow in the subcooled region.
        R_a                         = 1. / (self.Fins.eta_a * self.Fins.h_a * self.Fins.A_a * self.h_a_tuning)
        R_r                         = 1. / (self.h_r_subcool * self.A_r_wetted)
        UA_subcool                  = self.w_subcool / (R_a + R_r + self.Rw)
        C_min                       = min([self.m_dot_da * self.Fins.cp_da * self.w_subcool, self.m_dot_r * cp_r])
        C_max                       = max([self.m_dot_da * self.Fins.cp_da * self.w_subcool, self.m_dot_r * cp_r])
        Cr                          = C_min / C_max
        NTU                         = UA_subcool / C_min

        if self.m_dot_da * self.Fins.cp_da * self.w_subcool > self.m_dot_r * cp_r:
            # Minimum capacitance rate on refrigerant side
            epsilon_subcool         = 1. - exp(-1. / Cr * (1. - exp(-Cr * NTU)))
        else:
            # Minimum capacitance rate on air side
            epsilon_subcool         = 1 / Cr * (1 - exp(-Cr * (1 - exp(-NTU))))

        # Effectiveness for both fluids unmixed:
        # epsilon_subcool = 1 - exp(1/Cr * NTU**0.22 * (exp(-Cr*NTU**0.78)-1))

        # Positive Q is heat input to the refrigerant, negative Q is heat output from refrigerant.
        # Heat is removed here from the refrigerant since it is condensing
        self.Q_subcool              = -epsilon_subcool * C_min * (T_bubble - self.T_in_a)
        self.DT_sc                  = -self.Q_subcool / (self.m_dot_r * cp_r)
        self.T_out_r                = T_bubble - self.DT_sc

        AS.update(CP.PT_INPUTS, self.p_sat_r, (T_bubble + self.T_out_r) / 2)
        rho_subcool                 = AS.rhomass()                                          # [kg/m^3]
        self.Charge_subcool         = self.w_subcool * self.V_r * rho_subcool

        # Pressure drop calculations for subcooled refrigerant
        v_r                         = 1 / rho_subcool
        # Pressure gradient using Darcy friction factor
        dpdz_r                      = -self.f_r_subcool * v_r * self.G_r**2 / (2 * self.Dh)  # Pressure gradient
        self.DP_r_subcool           = dpdz_r * self.L_circuit * self.w_subcool


# -----------------------------------------------------------------------------------------
def SampleMicroCondenser(AS, T=95):
    Fins                                                = MicroFinInputs()
    Fins.Tubes.N_Tubes                                  = 61.354                                # Number of tubes (per bank for now!)
    Fins.Tubes.N_bank                                   = 1                                     # Number of banks (set to 1 for now!)
    Fins.Tubes.N_pass                                   = 3                                     # Number of passes (per bank-averaged)
    Fins.Tubes.N_ports                                  = 1                                     # Number of rectangular ports
    Fins.Tubes.L_tube                                   = 0.30213                               # length of a single tube
    Fins.Tubes.Td                                       = 0.0333                                # Tube outside width (depth)
    Fins.Tubes.Ht                                       = 0.002                                 # Tube outside height (major diameter)
    Fins.Tubes.b                                        = 0.00635                               # Tube spacing
    Fins.Tubes.tw                                       = 0.0003                                # Tube wall thickness
    Fins.Tubes.twp                                      = 0.0003                                # Port (channel) wall thickness
    Fins.Tubes.beta                                     = 1                                     # Port (channel) aspect ratio (=width/height)
    Fins.Tubes.k_w                                      = 117                                   # wall thermal conductivity

    Fins.Fins.FPI                                       = 11.0998                               # Fin per inch
    Fins.Fins.Lf                                        = 0.0333                                # Fin length
    Fins.Fins.t                                         = 0.000152                              # Fin thickness
    Fins.Fins.k_fin                                     = 117                                   # Fin thermal conductivity

    Fins.Air.V_dot_ha                                   = 1.05                                  # Air volume flow rate in m^3/s
    Fins.Air.T_mean                                     = 298
    Fins.Air.T_db                                       = 298                                   # Air inlet temperature, K
    Fins.Air.p                                          = 100000                                # Air pressure in Pa
    Fins.Air.RH_mean                                    = 0.5
    Fins.Air.RH                                         = 0.5                                   # Air inlet relative humidity
    Fins.Air.FanPower                                   = 854.9                                 # Fan power, Watts

    Fins.Louvers.Lalpha                                 = 20                                    # Louver angle, in degree
    Fins.Louvers.lp                                     = 0.001                                 # Louver pitch
    Fins.Louvers.Llouv                                  = 0.005737                              # Louver cut length

    params      = {
        'AS':                                           AS,
        'm_dot_r':                                      0.0683,
        'T_in_r':                                       T+273.15,
        'p_sat_r':                                      3500000,
        'Fins':                                         Fins,
        'FinsType':                                     'MultiLouveredMicroFins',
        'Verbosity':                                    0,
        'h_a_tuning':                                   1,
        'h_tp_tuning':                                  1,
        'DP_tuning':                                    1
    }
    MicroCond                                           = MicroCondenserClass(**params)
    MicroCond.Calculate()
    return MicroCond


# -----------------------------------------------------------------------------------------
if __name__ == '__main__':
    # This runs if you run this file directly
    Ref                                                 = 'R410A'
    Backend                                             = 'HEOS'        # choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    AS                                                  = CP.AbstractState(Backend, Ref)        # Abstract State
    MicroCond                                           = SampleMicroCondenser(AS, 95)
    # print (MicroCond.OutputList())
    print('Heat transfer rate in condenser is',                     MicroCond.Q,            'W')
    print('Heat transfer rate in condenser (superheat section) is', MicroCond.Q_superheat,  'W')
    print('Heat transfer rate in condenser (twophase section) is',  MicroCond.Q_2phase,     'W')
    print('Heat transfer rate in condenser (subcooled section) is', MicroCond.Q_subcool,    'W')
    print('Fraction of circuit length in superheated section is',   MicroCond.w_superheat)
    print('Fraction of circuit length in twophase section is',      MicroCond.w_2phase)
    print('Fraction of circuit length in subcooled section is',     MicroCond.w_subcool)
