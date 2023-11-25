from __future__                                             import division, print_function, absolute_import
from math                                                   import pi,log,exp
from scipy.optimize                                         import brentq # solver to find roots (zero points) of functions
from scipy.interpolate                                      import interp1d
import numpy                                                as np
import CoolProp                                             as CP

from ACHP_codes.Correlations.HTC_Correlations               import f_h_1phase_Tube, ShahEvaporation_Average, KandlikarEvaporation_average
from ACHP_codes.Correlations.deltaP_Correlations            import LMPressureGradientAvg, AccelPressureDrop
from ACHP_codes.Correlations.PhaseState_Correlations        import TwoPhaseDensity
from ACHP_codes.Correlations.FinStructure_Correlations      import WavyLouveredFins, FinInputs, IsFinsClass, HerringboneFins, PlainFins
from ACHP_codes.Correlations.DryWetSegment                  import DWSVals, DryWetSegment
from ACHP_codes.ACHP_Tools.ACHP_Tools                       import ValidateFields


# ---------------------------------------------------------------------------------------
class EvaporatorClass():
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def Update(self, **kwargs):
        self.__dict__.update(kwargs)

    def OutputList(self):
        """
            Return a list of parameters for this component for further output

            It is a list of tuples, and each tuple is formed of items:
                [0] Description of value
                [1] Units of value
                [2] The value itself
        """
        Output_List             = []
        # append optional parameters, if applicable
        if hasattr(self, 'TestName'):
            Output_List.append(('Name', 'N/A', self.TestName))
        if hasattr(self, 'TestDescription'):
            Output_List.append(('Description', 'N/A', self.TestDescription))
        if hasattr(self, 'TestDetails'):
            Output_List.append(('Details', 'N/A', self.TestDetails))

        Output_List_default = [                      # default output list
            ('Volumetric flow rate',                'm^3/s',                self.Fins.Air.V_dot_ha),
            ('Inlet Dry bulb temp',                 'K',                    self.T_in_a),
            ('Inlet Air pressure',                  'Pa',                   self.Fins.Air.p),
            ('Inlet Air Relative Humidity',         '-',                    self.Fins.Air.RH),
            ('Tubes per bank',                      '-',                    self.Fins.Tubes.N_Tubes_per_bank),
            ('Number of banks',                     '-',                    self.Fins.Tubes.N_bank),
            ('Number circuits',                     '-',                    self.Fins.Tubes.N_circuits),
            ('Length of tube',                      'm',                    self.Fins.Tubes.L_tube),
            ('Tube OD',                             'm',                    self.Fins.Tubes.OD),
            ('Tube ID',                             'm',                    self.Fins.Tubes.ID),
            ('Tube Long. Pitch',                    'm',                    self.Fins.Tubes.Pl),
            ('Tube Transverse Pitch',               'm',                    self.Fins.Tubes.Pt),
            ('Tube Conductivity',                   'W/m-K',                self.Fins.Tubes.k_w),
            ('Outlet superheat',                    'K',                    self.T_out_r - self.T_dew_r),
            ('Fins per inch',                       '1/in',                 self.Fins.Fins.FPI),
            ('Fin waviness pd',                     'm',                    self.Fins.Fins.Pd),
            ('Fin waviness xf',                     'm',                    self.Fins.Fins.xf),
            ('Fin thickness',                       'm',                    self.Fins.Fins.t),
            ('Fin Conductivity',                    'W/m-K',                self.Fins.Fins.k_fin),
            ('Fins Type',                           '-',                    self.FinsType),
            ('Q Total',                             'W',                    self.Q),
            ('Q Superheat',                         'W',                    self.Q_superheat),
            ('Q Two-Phase',                         'W',                    self.Q_2phase),
            ('Inlet ref. temp',                     'K',                    self.T_in_r),
            ('Outlet ref. temp',                    'K',                    self.T_out_r),
            ('Outlet air temp',                     'K',                    self.T_out_a),
            ('Evaporator P_sat in',                 'Pa',                   self.p_sat_r),
            ('Evaporator inlet quality',            '-',                    self.x_in_r),
            ('Evaporator ref. flowrate',            'kg/s',                 self.m_dot_r),
            ('Pressure Drop Total',                 'Pa',                   self.DP_r),
            ('Pressure Drop Superheat',             'Pa',                   self.DP_r_superheat),
            ('Pressure Drop Two-Phase',             'Pa',                   self.DP_r_2phase),
            ('Charge Total',                        'kg',                   self.Charge),
            ('Charge Superheat',                    'kg',                   self.Charge_superheat),
            ('Charge Two-Phase',                    'kg',                   self.Charge_2phase),
            ('Mean HTC Superheat',                  'W/m^2-K',              self.h_r_superheat),
            ('Mean HTC Two-phase',                  'W/m^2-K',              self.h_r_2phase),
            ('Wetted Area Fraction Superheat',      '-',                    self.w_superheat),
            ('Wetted Area Fraction Two-phase',      '-',                    self.w_2phase),
            ('Mean Air HTC',                        'W/m^2-K',              self.Fins.h_a * self.h_a_tuning),
            ('Surface Effectiveness',               '-',                    self.Fins.eta_a),
            ('Air-side area (fin+tubes)',           'm^2',                  self.Fins.A_a),
            ('Mass Flow rate dry Air',              'kg/s',                 self.Fins.m_dot_da),
            ('Mass Flow rate humid Air',            'kg/s',                 self.Fins.m_dot_ha),
            ('Pressure Drop Air-side',              'Pa',                   self.Fins.dP_a),
            ('Sensible Heat Ratio',                 '-',                    self.SHR)
            # ('Bend Temperature profile','K',self.T_bends)
        ]

        for i in range(0, len(Output_List_default)):         # append default parameters to output list
            Output_List.append(Output_List_default[i])
        return Output_List

    def AirSideCalcs(self):
        # Update with user FinType
        if self.FinsType == 'WavyLouveredFins':
            WavyLouveredFins(self.Fins)
        elif self.FinsType == 'HerringboneFins':
            HerringboneFins(self.Fins)
        elif self.FinsType == 'PlainFins':
            PlainFins(self.Fins)

    def Initialize(self):
        # Input validation the first call of Initialize
        # if False:
        if not hasattr(self, 'IsValidated'):
            self.Fins.Validate()
            reqFields = [
                       ('p_sat_r',      float,          0.001,          100000000),
                       ('Fins',         IsFinsClass,    None,           None),
                       ('FinsType',     str,            None,           None),
                       ('h_in_r',       float,          -100000,        10000000),
                       ('m_dot_r',      float,          0.000001,       10),
                       ]
            optFields           = ['Verbosity', 'AS', 'h_a_tuning', 'h_tp_tuning', 'DP_tuning']
            d                   = self.__dict__                 # Current fields in model
            ValidateFields(d, reqFields, optFields)
            self.IsValidated    = True

        # set tuning factors to 1 in case not given by user
        if not hasattr(self, 'h_a_tuning'):
            self.h_a_tuning     = 1
        if not hasattr(self, 'h_tp_tuning'):
            self.h_tp_tuning    = 1
        if not hasattr(self, 'DP_tuning'):
            self.DP_tuning      = 1

        # Make sure AbstractState is passed
        assert hasattr(self, 'AS'), 'Please specify the Abstract State'

        # Retrieve some parameters from nested structures for code compactness
        self.ID                 = self.Fins.Tubes.ID
        self.OD                 = self.Fins.Tubes.OD
        self.L_tube             = self.Fins.Tubes.L_tube
        self.N_Tubes_per_bank   = self.Fins.Tubes.N_Tubes_per_bank
        self.N_bank             = self.Fins.Tubes.N_bank
        self.N_circuits         = self.Fins.Tubes.N_circuits
        self.T_in_a             = self.Fins.Air.T_db
        self.k_w                = self.Fins.Tubes.k_w           # thermal conductivity of tube wall

        # Calculate an effective length of circuit if circuits are not all the same length
        TotalLength             = self.L_tube * self.N_Tubes_per_bank * self.N_bank
        self.L_circuit          = TotalLength / self.N_circuits
        # Wetted area on the refrigerant side
        self.A_r_wetted         = self.N_circuits * pi * self.ID * self.L_circuit
        self.V_r                = self.N_circuits * self.L_circuit * pi * self.ID**2 / 4.0
        # Average mass flux of refrigerant in circuit
        self.G_r                = self.m_dot_r / (self.N_circuits * pi * self.ID**2 / 4.0)  # [kg/m^2-s]

        # Thermal resistance at the wall
        self.R_w                = log(self.OD / self.ID) / (2 * pi * self.k_w * self.L_circuit * self.N_circuits)

        # Bubble and dew temperatures (same for fluids without glide)
        self.AS.update(CP.PQ_INPUTS, self.p_sat_r, 0.0)
        self.T_bubble_r         = self.AS.T()                   # [K]
        h_l                     = self.AS.hmass()               # [J/kg]
        self.AS.update(CP.PQ_INPUTS, self.p_sat_r, 1.0)
        self.T_dew_r            = self.AS.T()                   # [K]
        h_v                     = self.AS.hmass()               # [J/kg]
        # Mean temperature for use in HT relationships
        self.T_sat_r            = (self.T_bubble_r + self.T_dew_r) / 2
        # Latent heat
        self.h_fg               = h_v - h_l                     # [J/kg]

        self.Fins.Air.RH_mean   = self.Fins.Air.RH

        # Update with user FinType
        if self.FinsType == 'WavyLouveredFins':
            WavyLouveredFins(self.Fins)
        elif self.FinsType == 'HerringboneFins':
            HerringboneFins(self.Fins)
        elif self.FinsType == 'PlainFins':
            PlainFins(self.Fins)

        self.m_dot_ha           = self.Fins.m_dot_ha            # [kg_ha/s]
        self.m_dot_da           = self.Fins.m_dot_da            # [kg_da/s]

    def Calculate(self):
        # Initialize
        self.Initialize()
        AS                          = self.AS

        # Input and output thermodynamic properties
        AS.update(CP.QT_INPUTS, 0.0, self.T_bubble_r)
        s_sat_L                     = AS.smass()                    # [J/kg-K]
        h_sat_L                     = AS.hmass()                    # [J/kg]
        AS.update(CP.QT_INPUTS, 1.0, self.T_dew_r)
        s_sat_V                     = AS.smass()                    # [J/kg-K]
        h_sat_V                     = AS.hmass()                    # [J/kg]

        # if we give enthalpy and pressure as inputs
        if hasattr(self, 'h_in_r'):
            self.x_in_r             = (self.h_in_r - h_sat_L) / (h_sat_V - h_sat_L)
            self.s_in_r             = self.x_in_r * s_sat_V + (1 - self.x_in_r) * s_sat_L
            self.T_in_r             = self.x_in_r * self.T_dew_r + (1 - self.x_in_r) * self.T_bubble_r
        # if given quality and pressure as inputs
        elif hasattr(self, 'x_in_r'):
            self.h_in_r             = self.x_in_r * h_sat_V + (1 - self.x_in_r) * h_sat_L
            self.s_in_r             = self.x_in_r * s_sat_V + (1 - self.x_in_r) * s_sat_L
            self.T_in_r             = self.x_in_r * self.T_dew_r + (1 - self.x_in_r) * self.T_bubble_r

        if self.x_in_r > 1.0:
            raise ValueError
        # Begin by assuming that you go all the way to saturated vapor at least
        self.x_out_2phase           = 1.0

        if self._TwoPhase_Forward(1.0) < 0:
            # Evaporator outlet is in two-phase region, use all the area and find the outlet quality
            existsSuperheat         = False
            self.w_2phase           = 1.0
            def OBJECTIVE(x_out):
                self.x_out_2phase   = x_out
                Q_target            = self.m_dot_r * (x_out - self.x_in_r) * (h_sat_V - h_sat_L)
                self._TwoPhase_Forward(self.w_2phase)
                return self.Q_2phase - Q_target
            # Use a solver to find outlet quality
            brentq(OBJECTIVE, self.x_in_r, 1.0)
            self.w_superheat                = 0.0
            self.Q_superheat                = 0.0
            self.h_r_superheat              = 0.0
            self.Re_r_superheat             = 0.0
            self.Charge_superheat           = 0.0
            self.Q_sensible_superheat       = 0.0
            self.T_out_a_superheat          = 0.0
            self.DP_r_superheat             = 0.0
            self.fdry_superheat             = 0.0
        else:
            existsSuperheat         = True
            # Evaporator outlet is in superheated region, everything is ok
            self.w_2phase                   = brentq(self._TwoPhase_Forward, 0.00000000001, 0.9999999999)
            self._Superheat_Forward(1 - self.w_2phase)

        self.Q                              = self.Q_superheat + self.Q_2phase
        self.Charge                         = self.Charge_superheat + self.Charge_2phase
        if self.Verbosity > 4:
            print(self.Q, "Evaporator.Q")
        self.Capacity                       = self.Q - self.Fins.Air.FanPower

        # Sensible heat ratio [-]
        self.SHR                            = (self.Q_sensible_2phase + self.Q_sensible_superheat) / self.Q
        # Average air outlet temperature (area fraction weighted average) [K]
        self.T_out_a                        = self.w_superheat * self.T_out_a_superheat + self.w_2phase * self.T_out_a_2phase
        self.DP_r                           = self.DP_r_superheat + self.DP_r_2phase
        self.DP_r                           = self.DP_r * self.DP_tuning            # correcting the total pressure drop

        # Outlet enthalpy obtained from energy balance
        self.h_out_r                        = self.h_in_r + self.Q / self.m_dot_r

        # Outlet entropy
        if existsSuperheat == True:
            AS.update(CP.PT_INPUTS, self.p_sat_r, self.T_out_r)
            self.s_out_r                    = AS.smass()                            # [J/kg-K]
        else:
            x_out_r                         = (self.h_out_r - h_sat_L) / (h_sat_V - h_sat_L)
            self.s_out_r                    = s_sat_V * x_out_r + (1 - x_out_r) * s_sat_L

        # Outlet superheat an temperature (in case of two phase)
        if existsSuperheat:
            self.DT_sh_calc                 = self.T_out_r - self.T_dew_r
        else:
            AS.update(CP.QT_INPUTS, 1.0, self.T_dew_r)
            cp_sh                           = AS.cpmass()                           # [J/kg-K]
            self.DT_sh_calc                 = (self.h_out_r - h_sat_V) / cp_sh      # Effective superheat
            AS.update(CP.PQ_INPUTS, self.p_sat_r + self.DP_r, x_out_r)
            self.T_out_r                    = AS.T()                                # saturated temperature at outlet quality [K]
        self.h_mean_r                       = self.w_2phase * self.h_r_2phase + self.w_superheat * self.h_r_superheat
        self.UA_r                           = self.h_mean_r * self.A_r_wetted
        self.UA_a                           = (self.Fins.h_a * self.h_a_tuning) * self.Fins.A_a * self.Fins.eta_a
        self.UA_w                           = 1 / self.R_w

        # Build a vector of temperatures at each point where there is a phase transition along the averaged circuit
        if existsSuperheat:
            # Insert the shoulder point
            Tv                              = [self.T_in_r, self.T_dew_r, self.T_out_r]
            x                               = [0, self.w_2phase, 1]
        else:
            Tv                              = [self.T_in_r, x_out_r * self.T_dew_r + (1 - x_out_r) * self.T_bubble_r]
            x                               = [0, 1]

        # Determine each bend temperature by interpolation
        # ------------------------------------------------
        # Number of bends (including inlet and outlet of coil)
        N_bends                             = 1 + self.L_circuit / self.L_tube
        # x-position of each point
        xv                                  = np.linspace(0, 1, int(N_bends))

        self.T_bends                        = interp1d(x, Tv)(xv)

    def _TwoPhase_Forward(self, w_2phase):

        DWS                                 = DWSVals()                             # DryWetSegment structure

        # Store temporary values to be passed to DryWetSegment
        DWS.Fins                            = self.Fins
        DWS.FinsType                        = self.FinsType
        DWS.A_a                             = self.Fins.A_a * w_2phase
        DWS.cp_da                           = self.Fins.cp_da
        DWS.eta_a                           = self.Fins.eta_a
        DWS.h_a                             = self.Fins.h_a * self.h_a_tuning       # Heat transfer coefficient, not enthalpy
        DWS.m_dot_da                        = self.m_dot_da * w_2phase
        DWS.p_in_a                          = self.Fins.Air.p
        DWS.T_dew_r                         = self.T_dew_r
        DWS.T_bubble_r                      = self.T_bubble_r

        DWS.T_in_a                          = self.T_in_a
        DWS.RH_in_a                         = self.Fins.Air.RH

        DWS.T_in_r                          = self.T_sat_r
        DWS.A_r                             = self.A_r_wetted * w_2phase
        DWS.R_w                             = self.R_w / w_2phase
        DWS.cp_r                            = 1.0e15        # In the two-phase region the cp is infinite, use 1e15 as a big number;
        DWS.p_in_r                          = self.p_sat_r
        DWS.m_dot_r                         = self.m_dot_r
        DWS.IsTwoPhase                      = True

        # Target heat transfer to go from inlet quality to saturated vapor
        Q_target                            = self.m_dot_r * (self.x_out_2phase - self.x_in_r) * self.h_fg

        if Q_target < 0:
            raise ValueError('Q_target in Evaporator must be positive')

        # Average Refrigerant heat transfer coefficient
        try:
            if self.AS.name() in 'CarbonDioxide':
                h_r = KandlikarEvaporation_average(self.x_in_r, self.x_out_2phase, self.AS, self.G_r, self.ID, self.p_sat_r, Q_target/DWS.A_r, self.T_bubble_r, self.T_dew_r)
            else:
                h_r = ShahEvaporation_Average(self.x_in_r, self.x_out_2phase, self.AS, self.G_r, self.ID, self.p_sat_r, Q_target/DWS.A_r, self.T_bubble_r, self.T_dew_r)
        except:
                h_r = ShahEvaporation_Average(self.x_in_r, self.x_out_2phase, self.AS, self.G_r, self.ID, self.p_sat_r, Q_target/DWS.A_r, self.T_bubble_r, self.T_dew_r)

        DWS.h_r                             = h_r * self.h_tp_tuning        # correct refrigerant side convection heat transfer

        # Run the DryWetSegment to carry out the heat and mass transfer analysis
        DryWetSegment(DWS)

        self.Q_2phase                       = DWS.Q
        self.Q_sensible_2phase              = DWS.Q_sensible
        self.h_r_2phase                     = DWS.h_r
        self.fdry_2phase                    = DWS.f_dry
        self.T_out_a_2phase                 = DWS.T_out_a

        rho_average                         = TwoPhaseDensity(self.AS, self.x_in_r, self.x_out_2phase, self.T_dew_r, self.T_bubble_r, slipModel='Zivi')
        self.Charge_2phase                  = rho_average * w_2phase * self.V_r

        # Frictional pressure drop component
        DP_frict    = LMPressureGradientAvg(self.x_in_r, self.x_out_2phase, self.AS, self.G_r, self.ID, self.T_bubble_r, self.T_dew_r) * self.L_circuit * w_2phase
        # Accelerational pressure drop component
        try:
            if self.AS.name() in 'CarbonDioxide':
                DP_accel    = AccelPressureDrop(self.x_in_r, self.x_out_2phase, self.AS, self.G_r, self.T_bubble_r, self.T_dew_r, D=self.ID, slipModel='Premoli') * self.L_circuit * w_2phase
            else:
                DP_accel    = AccelPressureDrop(self.x_in_r, self.x_out_2phase, self.AS, self.G_r, self.T_bubble_r, self.T_dew_r, slipModel='Zivi') * self.L_circuit * w_2phase
        except:
                DP_accel    = AccelPressureDrop(self.x_in_r, self.x_out_2phase, self.AS, self.G_r, self.T_bubble_r, self.T_dew_r, slipModel='Zivi') * self.L_circuit * w_2phase

        self.DP_r_2phase                    = DP_frict + DP_accel

        if self.Verbosity > 7:
            print(w_2phase, DWS.Q, Q_target, self.x_in_r, "w_2phase, DWS.Q, Q_target, self.x_in_r")

        return DWS.Q - Q_target

    def _Superheat_Forward(self, w_superheat):
        self.w_superheat                    = w_superheat
        DWS                                 = DWSVals()                         # DryWetSegment structure
        AS                                  = self.AS                           # AbstractState

        # Store temporary values to be passed to DryWetSegment
        DWS.A_a                             = self.Fins.A_a * w_superheat
        DWS.cp_da                           = self.Fins.cp_da
        DWS.eta_a                           = self.Fins.eta_a
        DWS.h_a                             = self.Fins.h_a * self.h_a_tuning   # Heat transfer coefficient
        DWS.m_dot_da                        = self.m_dot_da * w_superheat
        DWS.p_in_a                          = self.Fins.Air.p
        DWS.Fins                            = self.Fins
        DWS.FinsType                        = self.FinsType

        # Inputs on the air side to two phase region are inlet air again
        DWS.T_in_a                          = self.T_in_a
        DWS.RH_in_a                         = self.Fins.Air.RH

        DWS.T_in_r                          = self.T_dew_r
        DWS.A_r                             = self.A_r_wetted * w_superheat
        DWS.R_w                             = self.R_w / w_superheat

        AS.update(CP.PT_INPUTS, self.p_sat_r, self.T_dew_r+2.5)
        DWS.cp_r                            = AS.cpmass()       # Use a guess value of 6K superheat to calculate cp [J/kg-K]
        DWS.p_in_r                          = self.p_sat_r
        DWS.m_dot_r                         = self.m_dot_r
        DWS.IsTwoPhase                      = False

        # Use a guess value of 6K superheat to calculate the properties
        self.f_r_superheat, self.h_r_superheat, self.Re_r_superheat = f_h_1phase_Tube(self.m_dot_r / self.N_circuits, self.ID,
            self.T_dew_r+3, self.p_sat_r, self.AS)

        # Average Refrigerant heat transfer coefficient
        DWS.h_r                             = self.h_r_superheat

        # Run DryWetSegment
        DryWetSegment(DWS)

        AS.update(CP.PT_INPUTS, self.p_sat_r, (DWS.T_out_r + self.T_dew_r) / 2.0)
        rho_superheat                       = AS.rhomass()                          # [kg/m^3]
        self.Charge_superheat               = w_superheat * self.V_r * rho_superheat

        # Pressure drop calculations for superheated refrigerant
        v_r                                 = 1 / rho_superheat
        # Pressure gradient using Darcy friction factor
        dpdz_r                              = -self.f_r_superheat * v_r * self.G_r**2 / (2 * self.ID)   # Pressure gradient
        self.DP_r_superheat                 = dpdz_r * self.L_circuit * self.w_superheat

        # Set values
        self.Q_superheat                    = DWS.Q
        self.Q_sensible_superheat           = DWS.Q_sensible
        self.fdry_superheat                 = DWS.f_dry
        self.T_out_a_superheat              = DWS.T_out_a
        self.T_out_r                        = DWS.T_out_r


# -----------------------------------------------------------------------------------------
if __name__ == '__main__':
    # Example usage for a parametric study
    from CoolProp.CoolProp import PropsSI
    import pylab

    num_points                  = 101
    T_dews                      = np.linspace(273, 299.7, num_points)
    TT                          = np.empty(num_points)
    Q_2p                        = np.empty(num_points)
    w_2p                        = np.empty(num_points)
    w_sh                        = np.empty(num_points)
    Q_tot                       = np.empty(num_points)
    h_2p                        = np.empty(num_points)
    h_sh                        = np.empty(num_points)

    FinsTubes                   = FinInputs()

    FinsTubes.Tubes.N_Tubes_per_bank            = 32
    FinsTubes.Tubes.N_circuits                  = 5
    FinsTubes.Tubes.N_bank                      = 3
    FinsTubes.Tubes.L_tube                      = 0.452
    FinsTubes.Tubes.OD                          = 0.009525
    FinsTubes.Tubes.ID                          = 0.0089154
    FinsTubes.Tubes.Pl                          = 0.0254
    FinsTubes.Tubes.Pt                          = 0.0219964
    FinsTubes.Tubes.k_w                         = 237                   # wall thermal conductivity (i.e pipe material)

    FinsTubes.Fins.FPI                          = 14.5
    FinsTubes.Fins.Pd                           = 0.001
    FinsTubes.Fins.xf                           = 0.001
    FinsTubes.Fins.t                            = 0.00011
    FinsTubes.Fins.k_fin                        = 237

    FinsTubes.Air.V_dot_ha                      = 0.5663
    FinsTubes.Air.T_mean                        = 299.9
    FinsTubes.Air.T_db                          = 299.9
    FinsTubes.Air.p                             = 101325
    FinsTubes.Air.RH                            = 0.51
    FinsTubes.Air.RH_mean                       = 0.51
    FinsTubes.Air.FanPower                      = 438

    # Abstract State
    Ref             = 'R410A'
    Backend         = 'HEOS'    # choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    AS              = CP.AbstractState(Backend, Ref)

    kwargs = {'AS': AS,
            'm_dot_r':          0.0708,
            'p_sat_r':          PropsSI('P', 'T', T_dews[0], 'Q', 1.0, Ref),
            'Fins':             FinsTubes,
            'FinsType':         'WavyLouveredFins',  # Choose fin Type: 'WavyLouveredFins' or 'HerringboneFins'or 'PlainFins'
            'h_in_r':           PropsSI('H', 'P', PropsSI('P', 'T', 282, 'Q', 1.0, Ref), 'Q', 0.15, Ref),
            'Verbosity':        0,
            'h_a_tuning':       1,
            'h_tp_tuning':      1,
            'DP_tuning':        1}

    Evap                                        = EvaporatorClass(**kwargs)     # generate new evaporator instance and update kwargs

    for i in range(0, len(T_dews)):
        kwargs                                  = {'p_sat_r':  PropsSI('P', 'T', T_dews[i], 'Q', 1.0, Ref)}
        Evap.Update(**kwargs)
        Evap.Calculate()

        Q_tot[i]                                = Evap.Q
        Q_2p[i]                                 = Evap.Q_2phase
        w_2p[i]                                 = Evap.w_2phase
        w_sh[i]                                 = Evap.w_superheat
        h_2p[i]                                 = Evap.h_r_2phase
        h_sh[i]                                 = Evap.h_r_superheat

    print("Demonstrate output list")
    # print (Evap.OutputList())
    for id, unit, value in Evap.OutputList():
        print(str(id) + ' = ' + str(value) + ' ' + str(unit))

    pylab.plot(T_dews, Q_2p, T_dews, Q_tot)
    pylab.title('Parametric Study With Fixed flow-rates - Capacity')
    pylab.legend(['two-phase', 'total'], loc='best')
    pylab.title('Parametric Study With Fixed flow-rates - Capacity')
    pylab.xlabel('Evaporation Dew Temperature in Kelvin')
    pylab.ylabel('Capacity in Watt')
    # pylab.savefig('Evaporator_py_capacity.pdf')
    pylab.show()
    pylab.plot(T_dews, h_2p, T_dews, h_sh)
    pylab.title('Parametric Study with fixed flow-rates - Heat Transfer Coefficients')
    pylab.legend(['two-phase', 'superheat'], loc='best')
    pylab.xlabel('Evaporation Dew Temperature in Kelvin')
    pylab.ylabel('Heat Transfer Coefficient in W/m2-K')
    # pylab.savefig('Evaporator_py_HTC.pdf')
    pylab.show()
    pylab.plot(T_dews, w_2p, T_dews, w_sh)
    pylab.title('Parametric Study with fixed flow-rates - Area Fraction')
    pylab.legend(['two-phase', 'superheat'], loc='best')
    pylab.xlabel('Evaporation Dew Temperature in Kelvin')
    pylab.ylabel('Two-phase Wetted Area Fraction')
    pylab.ylim(-0.01, 1.01)
    # pylab.savefig('Evaporator_py_wetted_area.pdf')
    pylab.show()
