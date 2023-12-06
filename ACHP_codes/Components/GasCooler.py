from __future__                                         import division, print_function, absolute_import
from math                                               import pi, log, exp

from ACHP_codes.Correlations.HTC_Correlations           import f_h_1phase_Tube, Petterson_supercritical_average
from ACHP_codes.Correlations.FinStructure_Correlations  import WavyLouveredFins, FinInputs, IsFinsClass, HerringboneFins, PlainFins
from ACHP_codes.Correlations.DryWetSegment              import DWSVals, DryWetSegment
from ACHP_codes.ACHP_Tools.ACHP_Tools                   import ValidateFields
from time                                               import time
from scipy.optimize                                     import brentq
import CoolProp                                         as CP


# ----------------------------------------------------------------------------
class FinVals():
    def __init__(self):
        pass


class GasCoolerClass():
    def __init__(self, **kwargs):
        # Load the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)

    def Update(self, **kwargs):
        # Update the parameters passed in
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
            ('Volumetric flow rate of humid air',               'm^3/s',                self.Fins.Air.Vdot_ha),
            ('Inlet Dry bulb temp',                             'K',                    self.Tin_a),
            ('Inlet Air pressure',                              'Pa',                   self.Fins.Air.p),
            ('Inlet Air Relative Humidity',                     '-',                    self.Fins.Air.RH),
            ('Tubes per bank',                                  '-',                    self.Fins.Tubes.N_Tubes_per_bank),
            ('Number of banks',                                 '-',                    self.Fins.Tubes.N_bank),
            ('Number circuits',                                 '-',                    self.Fins.Tubes.N_circuits),
            ('Length of tube',                                  'm',                    self.Fins.Tubes.L_tube),
            ('Tube OD',                                         'm',                    self.OD),
            ('Tube ID',                                         'm',                    self.ID),
            ('Tube Long. Pitch',                                'm',                    self.Fins.Tubes.Pl),
            ('Tube Transverse Pitch',                           'm',                    self.Fins.Tubes.Pt),
            ('Tube Conductivity',                               'W/m-K',                self.Fins.Tubes.k_w),
            ('Fins per inch',                                   '1/in',                 self.Fins.Fins.FPI),
            ('Fin waviness pd',                                 'm',                    self.Fins.Fins.Pd),
            ('Fin waviness xf',                                 'm',                    self.Fins.Fins.xf),
            ('Fin thickness',                                   'm',                    self.Fins.Fins.t),
            ('Fin Conductivity',                                'W/m-K',                self.Fins.Fins.k_fin),
            ('Fins Type',                                       '-',                    self.FinsType),
            ('Q Total',                                         'W',                    self.Q),
            ('Q Supercritical',                                 'W',                    self.Q_supercritical),
            ('Q Supercritical_liquid',                          'W',                    self.Q_supercrit_liq),
            ('Inlet Temp',                                      'K',                    self.T_in_r),
            ('Outlet Temp',                                     'K',                    self.T_out_r),
            ('Pressure Drop Total',                             'Pa',                   self.DP_r),
            ('Pressure Drop Supercritical',                     'Pa',                   self.DP_r_supercritical),
            ('Pressure Drop Supercritical_liquid',              'Pa',                   self.DP_r_supercrit_liq),
            ('Charge Total',                                    'kg',                   self.Charge),
            ('Charge Supercritical',                            'kg',                   self.Charge_supercritical),
            ('Charge Supercritical_liquid',                     'kg',                   self.Charge_supercrit_liq),
            ('Mean HTC Superheat',                              'W/m^2-K',              self.h_r_supercritical),
            ('Mean HTC Supercritical_liquid',                   'W/m^2-K',              self.h_r_supercrit_liq),
            ('Wetted Area Fraction Supercritical',              '-',                    self.w_supercritical),
            ('Wetted Area Fraction Supercritical_liquid',       '-',                    self.w_supercrit_liq),
            ('Mean Air HTC',                                    'W/m^2-K',              self.Fins.h_a * self.h_a_tuning),
            ('Surface Effectiveness',                           '-',                    self.Fins.eta_a),
            ('Air-side area (fin+tubes)',                       'm^2',                  self.Fins.A_a),
            ('Mass Flow rate of Dry Air',                       'kg/s',                 self.Fins.m_dot_da),
            ('Mass Flow rate of Humid Air',                     'kg/s',                 self.Fins.m_dot_ha),
            ('Pressure Drop Air-side',                          'Pa',                   self.Fins.dP_a),
            ('Approach temperature degree',                     'K',                    self.DT_app),
        ]

    def Initialize(self):
        # Only validate the first time
        if not hasattr(self, 'IsValidated'):
            self.Fins.Validate()
            reqFields   = [
               ('Fins',         IsFinsClass,    None,       None),
               ('FinsType',     str,            None,       None),
               ('m_dot_r',      float,          0.00001,    20),
               ('T_in_r',       float,          200,        500),
               ('p_sat_r',      float,          0.01,       20000000)]

            optFields           = ['Verbosity', 'AS', 'h_a_tuning', 'h_r_tuning', 'DP_tuning']
            ValidateFields(self.__dict__, reqFields, optFields)
            self.IsValidated    = True

        # set tuning factors to 1 in case not given by user
        if not hasattr(self, 'h_a_tuning'):
            self.h_a_tuning                     = 1
        if not hasattr(self, 'h_r_tuning'):
            self.h_r_tuning                     = 1
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
        self.T_in_a                             = self.Fins.Air.T_db
        self.k_w                                = self.Fins.Tubes.k_w           # thermal conductivity of tube wall

        # Calculate an effective length of circuit if circuits are
        # not all the same length
        TotalLength                             = self.L_tube * self.N_Tubes_per_bank * self.N_bank
        self.L_circuit                          = TotalLength / self.N_circuits
        self.V_r                                = pi * self.ID**2 / 4.0 * self.L_circuit * self.N_circuits
        self.A_r_wetted                         = pi * self.ID * self.N_circuits * self.L_circuit
        self.G_r                                = self.m_dot_r / (self.N_circuits * pi * self.ID**2 / 4.0)

        # Thermal resistance at the wall
        self.R_w                                = log(self.OD / self.ID) / (2 * pi * self.k_w * self.L_circuit * self.N_circuits)

        # Define known parameters
        AS.update(CP.PT_INPUTS, self.p_sat_r, self.T_in_r)
        self.h_in_r                             = AS.hmass()                    # [J/kg]
        self.s_in_r                             = AS.smass()                    # [J/kg-K]

        # Define critical pressure and temperature
        self.P_cr                               = AS.p_critical()               # [Pa]
        self.T_cr                               = AS.T_critical()               # [K]

        # critical enthalpy at defined pressure
        AS.update(CP.PT_INPUTS, self.p_sat_r, self.T_cr)
        self.h_cr                               = AS.hmass()                    # [J/kg]

        # triple temperature
        self.T_triple                           = AS.Ttriple()

        self.Fins.Air.RH_mean                   = self.Fins.Air.RH

        # Update with user FinType
        if self.FinsType == 'WavyLouveredFins':
            WavyLouveredFins(self.Fins)
        elif self.FinsType == 'HerringboneFins':
            HerringboneFins(self.Fins)
        elif self.FinsType == 'PlainFins':
            PlainFins(self.Fins)

        self.m_dot_ha                            = self.Fins.m_dot_ha           # [kg_ha/s]
        self.m_dot_da                            = self.Fins.m_dot_da           # [kg_da/s]

    def Calculate(self):
        # Initialize
        self.Initialize()
        AS                                      = self.AS

        # assume we have all supercritical region
        self.T_out_r_cr                         = self.T_cr
        # give an initial guess for the inner wall temperature
        self.T_w                                = (self.T_out_r_cr + self.T_in_a) / 2
        # If we have already used too much of the HX (max possible sum of w is 1.0)
        if self._Supercritical_Forward(1.0) < 0:
            self.existsSubcooled                = False
            self.w_supercritical                = 1.0

            def OBJECTIVE(T_out_r_cr):
                self.T_out_r_cr                 = T_out_r_cr
                AS.update(CP.PT_INPUTS, self.p_sat_r, self.T_out_r_cr)
                h_out                           = AS.hmass()
                Q_target                        = self.m_dot_r * (h_out - self.h_in_r)
                self._Supercritical_Forward(self.w_supercritical)
                return self.Q_supercritical - Q_target

            brentq(OBJECTIVE, self.T_in_r, self.T_cr)
            # Zero out all the supercritical_liquid parameters
            self.Q_supercrit_liq                = 0.0
            self.DP_r_supercrit_liq             = 0.0
            self.Charge_supercrit_liq           = 0.0
            self.w_supercrit_liq                = 0.0
            self.h_r_supercrit_liq              = 0.0
            self.Re_r_supercrit_liq             = 0.0
            self.T_out_a_supercrit_liq          = 0.0
            self.fdry_supercrit_liq             = 0.0
        else:
            # By definition then we have a supercritical_liquid portion, solve for it
            self.existsSubcooled                = True
            self.w_supercritical                = brentq(self._Supercritical_Forward, 0.00000000001, 0.9999999999)

            def OBJECTIVE(T_out_r_sc):
                self.T_out_r_sc                 = T_out_r_sc
                AS.update(CP.PT_INPUTS, self.p_sat_r, self.T_out_r_sc)
                h_out                           = AS.hmass()
                Q_target                        = self.m_dot_r * (h_out - self.h_cr)
                self._Supercrit_liq_Forward(1 - self.w_supercritical)
                return self.Q_supercrit_liq - Q_target

            brentq(OBJECTIVE, self.T_cr, self.T_triple)

        # Overall calculations
        self.Q                                  = self.Q_supercritical + self.Q_supercrit_liq
        self.DP_r                               = (self.DP_r_supercritical + self.DP_r_supercrit_liq) * self.DP_tuning
        self.Charge                             = self.Charge_supercritical + self.Charge_supercrit_liq

        # Average air outlet temperature (area fraction weighted average) [K]
        self.T_out_a        = self.w_supercritical * self.T_out_a_supercritical + self.w_supercrit_liq * self.T_out_a_supercrit_liq

        # Outlet enthalpy obtained from energy balance
        self.h_out_r                            = self.h_in_r + self.Q / self.m_dot_r

        AS.update(CP.HmassP_INPUTS, self.h_out_r, self.p_sat_r)
        self.T_out_r                            = AS.T()                            # [K]
        self.s_out_r                            = AS.smass()                        # [J/kg-K]

        # Approach temperature
        self.DT_app                             = self.T_out_r - self.T_in_a

        self.h_mean_r       = self.w_supercritical * self.h_r_supercritical + self.w_supercrit_liq * self.h_r_supercrit_liq
        self.UA_r                               = self.h_mean_r * self.A_r_wetted
        self.UA_a                               = (self.Fins.h_a * self.h_a_tuning) * self.Fins.A_a * self.Fins.eta_a
        self.UA_w                               = 1 / self.R_w

    def _Supercritical_Forward(self, w_supercritical):
        # **********************************************************************
        #                      SUPERCRITICAL PART
        # **********************************************************************
        # AbstractState
        AS                          = self.AS

        DWS                         = DWSVals()                         # DryWetSegment structure (only dry-analysis, single phase is used)

        # Store temporary values to be passed to DryWetSegment
        DWS.Fins                    = self.Fins
        DWS.FinsType                = self.FinsType
        DWS.A_a                     = self.Fins.A_a*w_supercritical
        DWS.cp_da                   = self.Fins.cp_da
        DWS.eta_a                   = self.Fins.eta_a
        DWS.h_a                     = self.Fins.h_a * self.h_a_tuning   # Heat transfer coefficient, not enthalpy
        DWS.m_dot_da                = self.m_dot_da * w_supercritical
        DWS.p_in_a                  = self.Fins.Air.p

        DWS.T_in_a                  = self.T_in_a
        DWS.RH_in_a                 = self.Fins.Air.RH  # relative humidity in [0, 1]

        DWS.T_in_r                  = self.T_in_r
        DWS.A_r                     = self.A_r_wetted * w_supercritical
        DWS.R_w                     = self.R_w / w_supercritical
        DWS.p_in_r                  = self.p_sat_r
        DWS.m_dot_r                 = self.m_dot_r
        DWS.IsTwoPhase              = False

        # ----------------------------------------------------------------------------
        AS.update(CP.PT_INPUTS, self.p_sat_r, self.T_out_r_cr)
        h_out                       = AS.hmass()                        # [J/kg]
        # Target heat transfer to go from inlet temperature to iterative outlet temperature
        Q_target                    = self.m_dot_r * (h_out - self.h_in_r)

        if Q_target > 0:
            raise ValueError('Q_target in Gas cooler must be negative')

        # This block calculates the average refrigerant heat transfer coefficient,
        # average friction factor, average specific heat, and average density
        h_r, f_r_supercritical, DWS.cp_r, rho_supercritical = Petterson_supercritical_average(self.T_out_r_cr, self.T_in_r,
            self.T_w, self.AS, self.G_r, self.ID, 0, self.ID/self.L_circuit, self.m_dot_r / self.N_circuits, self.p_sat_r, -Q_target/DWS.A_r)

        DWS.h_r                 = h_r * self.h_r_tuning             # Correct refrigerant HTC with tuning factor

        # Compute Fins Efficiency based on FinsType
        DryWetSegment(DWS)

        self.T_w                    = DWS.T_wall_s                  # inner surface wall temperature (refrigerant)
        self.Q_supercritical        = DWS.Q
        self.h_r_supercritical      = DWS.h_r
        self.fdry_supercritical     = DWS.f_dry
        self.T_out_a_supercritical  = DWS.T_out_a

        # Pressure drop calculations for supercritical refrigerant
        v_r                         = 1. / rho_supercritical
        # Pressure gradient using Darcy friction factor
        dpdz_r                      = -f_r_supercritical * v_r * self.G_r**2 / (2 * self.ID)    # Pressure gradient
        self.DP_r_supercritical     = dpdz_r * self.L_circuit * w_supercritical
        # charge for the supercritical portion
        self.Charge_supercritical   = w_supercritical * self.V_r * rho_supercritical

        if self.Verbosity > 7:
            print(w_supercritical, DWS.Q, Q_target, "w_supercritical, DWS.Q, Q_target")

        return Q_target - DWS.Q

    def _Supercrit_liq_Forward(self, w_supercrit_liq):
        # **********************************************************************
        #                      SUPERCRITICAL_LIQUID PART
        # **********************************************************************
        self.w_supercrit_liq                    = w_supercrit_liq

        if self.w_supercrit_liq < 0:
            raise ValueError('w_supercrit_liq in Gas cooler cannot be less than zero')

        # AbstractState
        AS                          = self.AS

        DWS                         = DWSVals()                         # DryWetSegment structure

        # Store temporary values to be passed to DryWetSegment
        DWS.A_a                     = self.Fins.A_a * w_supercrit_liq
        DWS.cp_da                   = self.Fins.cp_da
        DWS.eta_a                   = self.Fins.eta_a
        DWS.h_a                     = self.Fins.h_a * self.h_a_tuning   # Heat transfer coefficient
        DWS.m_dot_da                = self.m_dot_da * w_supercrit_liq
        DWS.p_in_a                  = self.Fins.Air.p
        DWS.Fins                    = self.Fins
        DWS.FinsType                = self.FinsType

        # Inputs on the air side to two phase region are inlet air again
        DWS.T_in_a                  = self.T_in_a
        DWS.RH_in_a                 = self.Fins.Air.RH

        DWS.T_in_r                  = self.T_cr
        DWS.A_r                     = self.A_r_wetted * w_supercrit_liq
        DWS.R_w                     = self.R_w / w_supercrit_liq

        AS.update(CP.PT_INPUTS, self.p_sat_r, self.T_cr-1)
        DWS.cp_r                    = AS.cpmass()                       # [J/kg-K]

        DWS.p_in_r                  = self.p_sat_r
        DWS.m_dot_r                 = self.m_dot_r
        DWS.IsTwoPhase              = False

        AS.update(CP.PT_INPUTS, self.p_sat_r, self.T_out_r_sc)
        h_out                       = AS.hmass()                        # [J/kg]
        # Target heat transfer to go from inlet temperature (critical) to iterative outlet temperature
        Q_target                    = self.m_dot_r * (h_out - self.h_cr)

        if Q_target > 0:
            raise ValueError('Q_target in Gas cooler must be negative')

        # Friction factor and HTC in the refrigerant portions.
        # Average fluid temps are used for the calculation of properties
        # Average temp of refrigerant is average of sat. temp and outlet temp
        # Secondary fluid is air over the fins
        h_r, self.f_r_supercrit_liq, DWS.cp_r, rho_supercrit_liq = Petterson_supercritical_average(self.T_out_r_sc, self.T_cr,
            self.T_w, self.AS, self.G_r, self.ID, 0, self.ID/self.L_circuit, self.m_dot_r / self.N_circuits, self.p_sat_r, -Q_target/DWS.A_r)
        DWS.h_r                     = h_r * self.h_r_tuning             # Correct refrigerant HTC with tuning factor

        # Run DryWetSegment
        DryWetSegment(DWS)

        # Positive Q is heat input to the refrigerant, negative Q is heat output from refrigerant.
        # Heat is removed here from the refrigerant since it is condensing
        self.T_w                    = DWS.T_wall_s                      # inner surface wall temperature (refrigerant)
        self.Q_supercrit_liq        = DWS.Q
        self.h_r_supercrit_liq      = DWS.h_r
        self.fdry_supercrit_liq     = DWS.f_dry
        self.T_out_a_supercrit_liq  = DWS.T_out_a
        self.T_out_r                = DWS.T_out_r

        self.Charge_supercrit_liq   = self.w_supercrit_liq * self.V_r * rho_supercrit_liq

        # Pressure drop calculations for supercrit_liqed refrigerant
        v_r                         = 1 / rho_supercrit_liq
        # Pressure gradient using Darcy friction factor
        dpdz_r                      = -self.f_r_supercrit_liq * v_r * self.G_r**2 / (2 * self.ID)   # Pressure gradient
        self.DP_r_supercrit_liq     = dpdz_r * self.L_circuit * self.w_supercrit_liq

        return Q_target - DWS.Q


def SampleGasCooler(AS):
    Fins                                = FinInputs()
    Fins.Tubes.N_Tubes_per_bank         = 18                    # number of tubes per bank or row
    Fins.Tubes.N_bank                   = 3                     # number of banks or rows
    Fins.Tubes.N_circuits               = 1                     # number of circuits
    Fins.Tubes.L_tube                   = 0.61                  # one tube length
    Fins.Tubes.OD                       = 7.9 / 1000
    Fins.Tubes.ID                       = 7.5 / 1000
    Fins.Tubes.Pl                       = 19 / 1000             # distance between center of tubes in flow direction
    Fins.Tubes.Pt                       = 25/1000               # distance between center of tubes orthogonal to flow direction
    Fins.Tubes.k_w                      = 237                   # wall thermal conductivity (i.e. pipe material)

    Fins.Fins.FPI                       = 1 / (1.5 / 1000 / 0.0254)     # Number of fins per inch
    Fins.Fins.Pd                        = 0.001                         # 2* amplitude of wavy fin
    Fins.Fins.xf                        = 0.001                         # 1/2 period of fin
    Fins.Fins.t                         = 0.13 / 1000                   # Thickness of fin material
    Fins.Fins.k_fin                     = 237                           # Thermal conductivity of fin material

    Fins.Air.V_dot_ha                   = 0.281                         # rated volumetric flowrate (m^3/s)
    Fins.Air.T_mean                     = 29.4+273.15
    Fins.Air.T_db                       = 29.4+273.15                   # Dry Bulb Temperature
    Fins.Air.p                          = 101325                        # Air pressure in Pa
    Fins.Air.RH                         = 0.5                           # Relative Humidity
    Fins.Air.RH_mean                    = 0.5
    Fins.Air.FanPower                   = 160

    params  = {
        'AS':                           AS,
        'm_dot_r':                      0.076,
        'T_in_r':                       110.6+273.15,
        'p_sat_r':                      11000000,
        'Fins':                         Fins,
        'FinsType':                     'WavyLouveredFins',     # Choose fin Type: 'WavyLouveredFins' or 'HerringboneFins' or 'PlainFins'
        'Verbosity':                    0,
        'h_a_tuning':                   1,
        'h_r_tuning':                   1,
        'DP_tuning':                    1}

    Cond                                = GasCoolerClass(**params)
    Cond.Calculate()
    return Cond


# ----------------------------------------------------------------------------------------
if __name__ == '__main__':
    # This runs if you run this file directly

    t1                                  = time()
    Ref                                 = 'R744'
    Backend                             = 'HEOS'    # choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    AS                                  = CP.AbstractState(Backend, Ref)            # Abstract State
    Cond                                = SampleGasCooler(AS)

    # print(Cond.OutputList())

    print('Heat transfer rate in gas cooler is',                                Cond.Q,                 'W')
    print('Heat transfer rate in gas cooler (supercritical section) is',        Cond.Q_supercritical,   'W')
    print('Heat transfer rate in gas cooler (supercritical_liquid section) is', Cond.Q_supercrit_liq,   'W')
    print('Fraction of circuit length in supercritical section is',             Cond.w_supercritical)
    print('Fraction of circuit length in supercritical_liquid section is',      Cond.w_supercrit_liq)
    print('Refrigerant outlet temperature is',                                  Cond.T_out_r-273.15,     'C')
    print('Air outlet temperature is',                                          Cond.T_out_a-273.15,     'C')
    print('Took ' + str(time() - t1) + ' seconds to run Gas Cooler model')
