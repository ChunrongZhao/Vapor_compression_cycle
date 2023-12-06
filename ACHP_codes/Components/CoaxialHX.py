from __future__                                     import division, print_function, absolute_import
from CoolProp.CoolProp                              import PropsSI
from ACHP_codes.Correlations.HTC_Correlations       import f_h_1phase_Annulus, f_h_1phase_Tube, ShahEvaporation_Average
from ACHP_codes.Correlations.deltaP_Correlations    import LMPressureGradientAvg, AccelPressureDrop
from ACHP_codes.Correlations.PhaseState_Correlations import TwoPhaseDensity
from math                                           import pi, exp, log
from scipy.optimize                                 import brentq
import numpy                                        as np
import CoolProp                                     as CP


# ----------------------------------------------------------------------------------
class CoaxialHXClass():
    def __init__(self, **kwargs):
        # Load the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)

    def Update(self, **kwargs):
        # Update the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)

        # Wetted area on the refrigerant side
        self.A_r_wetted=pi*self.ID_i*self.L
        # Wetted area of the glycol side (not including outer tube)
        self.A_g_wetted=pi*self.OD_i*self.L

        self.V_r=self.L*pi*self.ID_i**2/4.0
        self.V_g=self.L*pi*(self.ID_o**2-self.OD_i**2)/4.0

    def OutputList(self):
        """
            Return a list of parameters for this component for further output

            It is a list of tuples, and each tuple is formed of items:
                [0] Description of value
                [1] Units of value
                [2] The value itself
        """
        return [
            ('Length of tube',                      'm',                    self.L),
            ('Annulus wetted OD',                   'm',                    self.ID_o),
            ('Tube wetted OD/Annulus wetted ID',    'm',                    self.OD_i),
            ('Tube wetted ID',                      'm',                    self.ID_i),
            ('Outlet Superheat',                    'K',                    self.T_in_r - self.T_dew_r),
            ('Q Total',                             'W',                    self.Q),
            ('Q Superheat',                         'W',                    self.Q_superheat),
            ('Q Two-Phase',                         'W',                    self.Q_2phase),
            ('Q Subcooled',                         'W',                    self.Q_subcool),
            ('Inlet glycol temp',                   'K',                    self.T_in_g),
            ('Outlet glycol temp',                  'K',                    self.T_out_g),
            ('Inlet ref. temp',                     'K',                    self.T_in_r),
            ('Outlet ref. temp',                    'K',                    self.T_out_r),
            ('Charge Total',                        'kg',                   self.Charge_r),
            ('Charge Superheat',                    'kg',                   self.Charge_r_superheat),
            ('Charge Two-Phase',                    'kg',                   self.Charge_r_2phase),
            ('Charge Subcool',                      'kg',                   self.Charge_r_subcool),
            ('Mean HTC Ref. Superheat',             'W/m^2-K',              self.h_r_superheat),
            ('Mean HTC Ref. Two-Phase',             'W/m^2-K',              self.h_r_2phase),
            ('Mean HTC Ref. Subcool',               'W/m^2-K',              self.h_r_subcool),
            ('Mean HTC Gly. Superheat',             'W/m^2-K',              self.h_g),
            ('Mean Reynolds # Gly. Superheat',      '-',                    self.Re_g),
            ('Pressure Drop Gly.',                  'Pa',                   self.DP_g),
            ('Pressure Drop Ref.',                  'Pa',                   self.DP_r),
            ('Pressure Drop Ref. Superheat',        'Pa',                   self.DP_r_superheat),
            ('Pressure Drop Ref. Two-Phase',        'Pa',                   self.DP_r_2phase),
            ('Pressure Drop Ref. Subcool',          'Pa',                   self.DP_r_subcool),
            ('Area fraction Superheat',             '-',                    self.w_superheat),
            ('Area fraction Two-Phase',             '-',                    self.w_2phase),
            ('Area fraction Subcooled',             '-',                    self.w_subcool)
         ]

    def Calculate(self):
        # AbstractState (Ref)
        AS_r                                = self.AS_r
        # AbstractState (Glycol) 乙二醇
        AS_g                                = self.AS_g
        if hasattr(self, 'MassFrac_g'):
            AS_g.set_mass_fractions([self.MassFrac_g])
        elif hasattr(self, 'VoluFrac_g'):
            AS_g.set_volu_fractions([self.VoluFrac_g])

        # set tuning factors to 1 in case not given by user
        if not hasattr(self,'h_g_tuning'):
            self.h_g_tuning                 = 1
        if not hasattr(self, 'h_tp_tuning'):
            self.h_tp_tuning                = 1
        if not hasattr(self, 'DP_g_tuning'):
            self.DP_g_tuning                = 1
        if not hasattr(self, 'DP_r_tuning'):
            self.DP_r_tuning                = 1

        # Update the parameters
        self.Update()

        # Average mass flux of refrigerant [kg/m^2-s]
        self.G_r                            = self.m_dot_r / (pi * self.ID_i**2 / 4.0)
        # Average mass flux of glycol [kg/m^2-s]
        self.G_g                            = self.m_dot_g / (pi * (self.ID_o**2 - self.OD_i**2) / 4.0)
        # Hydraulic diameter
        self.Dh_g                           = self.ID_o - self.OD_i
        # Evaporation hydraulic diameter [m]
        self.Dh_r                           = self.ID_i
        # Thermal conductivity of the intermediate wall pipe
        self.k                              = self.Conductivity

        # Thermal Conduction Resistance of the intermediate wall
        self.R_w                            = log(self.OD_i / self.ID_i) / (2 * pi * self.k * self.L)

        AS_r.update(CP.PQ_INPUTS, self.p_in_r, 0.0)
        self.T_bubble_r                     = AS_r.T()                          # [K]
        h_sat_L                             = AS_r.hmass()                      # [J/kg]
        s_sat_L                             = AS_r.smass()                      # [J/kg-K]
        AS_r.update(CP.PQ_INPUTS, self.p_in_r, 1.0)
        self.T_dew_r                        = AS_r.T()                          # [K]
        h_sat_V                             = AS_r.hmass()                      # [J/kg]
        s_sat_V                             = AS_r.smass()                      # [J/kg-K]

        # Saturation temperature
        self.T_sat_r                        = (self.T_bubble_r + self.T_dew_r) / 2.0

        # Inlet quality
        self.x_in_r                         = (self.h_in_r - h_sat_L) / (h_sat_V - h_sat_L)

        # Change in enthalpy through two-phase region [J/kg]
        self.h_fg                           = h_sat_V - h_sat_L
        self.T_in_r                         = self.x_in_r * self.T_dew_r + (1 - self.x_in_r) * self.T_bubble_r
        # Inlet entropy
        self.s_in_r                         = self.x_in_r * s_sat_V + (1 - self.x_in_r) * s_sat_L

        # Mean values for the glycol side based on average of inlet temperatures
        T_avg_g                             = (self.T_sat_r + self.T_in_g) / 2.0
        self.f_g, self.h_g, self.Re_g       = f_h_1phase_Annulus(self.m_dot_g, self.ID_o, self.OD_i, T_avg_g, self.p_in_g, self.AS_g)
        self.h_g                            = self.h_g * self.h_g_tuning        # correct h_g with tuning factor

        AS_g.update(CP.PT_INPUTS, self.p_in_g, T_avg_g)
        self.cp_g                           = AS_g.cpmass()                     # [J/kg-K]
        v_g                                 = 1 / AS_g.rhomass()                # [m^3/kg]

        # Glycol pressure drop
        dpdz_g                              = -self.f_g * v_g * self.G_g**2 / (2. * self.Dh_g)  # Pressure gradient
        self.DP_g                           = dpdz_g * self.L * self.DP_g_tuning

        def OBJECTIVE(w_superheat):
            """Nested function for driving the Brent's method solver"""
            # Run the superheated portion
            self._Superheat_Forward(w_superheat)
            # Run the two-phase portion and return residual
            return self._TwoPhase_Forward(1 - w_superheat)

        def OBJECTIVE_2phase(x_out_2phase):
            """Nested function for finding outlet quality"""
            # Need to pass in outlet quality but still maintain the full w_2phase
            return self._TwoPhase_Forward(1.0, x_out_2phase)

        # First see if you have a superheated portion.  Try to use the entire HX
        # for the two-phase portion
        # --------------------------------------------------------------------------
        # Intermediate glycol temp between superheated and 2phase sections [K]
        self.T_g_x                          = self.T_in_g
        # Call two-phase forward method
        error                               = self._TwoPhase_Forward(1.0)
        # If HT greater than required
        if error > 0:
            # Too much HT if all is 2phase, there is a superheated section
            existsSuperheat                 = True
            # Solve for the break between 2phase and superheated parts
            w_superheat                     = brentq(OBJECTIVE, 0.00001, 0.99999)
            self.w_2phase                   = 1 - w_superheat
            self.w_superheat                = w_superheat
        else:
            existsSuperheat                 = False
            # Solve for outlet quality in 2phase section, the lowest possible outlet
            # quality is the inlet quality
            self.x_out_2phase               = brentq(OBJECTIVE_2phase, self.x_in_r, 0.99999)
            # Dummy variables for the superheated section which doesn't exist
            self.Q_superheat                = 0.0
            self.Charge_r_superheat         = 0.0
            self.h_r_superheat              = 0.0
            self.DP_r_superheat             = 0.0
            self.w_superheat                = 0.0
            self.w_2phase                   = 1.0

        self.Charge_r                       = self.Charge_r_2phase + self.Charge_r_superheat
        self.Q                              = self.Q_2phase + self.Q_superheat
        self.T_out_g                        = self.T_in_g - self.Q / (self.cp_g * self.m_dot_g)

        self.DP_r                           = (self.DP_r_2phase + self.DP_r_superheat) * self.DP_r_tuning

        if existsSuperheat == True:
            AS_r.update(CP.PT_INPUTS, self.p_in_r, self.T_out_r)
            self.h_out_r                    = AS_r.hmass()                      # [J/kg]
            self.s_out_r                    = AS_r.smass()                      # [J/kg-K]
        else:
            self.T_out_r                    = self.x_out_2phase * self.T_dew_r + (1 - self.x_out_2phase) * self.T_bubble_r
            AS_r.update(CP.QT_INPUTS, self.x_out_2phase, self.T_out_r)
            self.h_out_r                    = AS_r.hmass()                      # [J/kg]
            self.s_out_r                    = AS_r.smass()                      # [J/kg-K]

        # #Dummy variables for the subcooled section which doesn't exist
        # self.Q_subcool=0.0
        # self.DP_r_subcool=0.0
        # self.h_r_subcool=0.0
        # self.Re_r_subcool=0.0
        # self.Charge_r_subcool=0.0
        # self.w_subcool=0.0

    def _Superheat_Forward(self, w_superheat):
        # Superheated portion
        # Mean temperature for superheated part can be taken to be average
        # of dew and glycol inlet temps
        T_avg_sh_r                          = (self.T_dew_r + self.T_in_g) / 2.0
        self.f_r_superheat, self.h_r_superheat, self.Re_r_superheat \
            = f_h_1phase_Tube(self.m_dot_r, self.ID_i, T_avg_sh_r, self.p_in_r, self.AS_r)
        # Refrigerant specific heat
        self.AS_r.update(CP.PT_INPUTS, self.p_in_r, T_avg_sh_r)
        cp_r_superheat                      = self.AS_r.cpmass()                # [J/kg-K]

        # Overall conductance of heat transfer surface in superheated portion
        UA_superheat    = w_superheat / (1 / (self.h_g * self.A_g_wetted) + 1 / (self.h_r_superheat * self.A_r_wetted) + self.R_w)
        # List of capacitance rates [W/K]
        C                                   = [cp_r_superheat * self.m_dot_r, self.cp_g * self.m_dot_g]
        C_min                               = min(C)
        Cr                                  = C_min / max(C)
        Ntu_superheat                       = UA_superheat / C_min

        # for Cr<1 and pure counter flow (Incropera. Table 11.3)
        epsilon_superheat       = ((1 - exp(-Ntu_superheat * (1 - Cr))) / (1 - Cr * exp(-Ntu_superheat * (1 - Cr))))

        self.Q_superheat        = epsilon_superheat * C_min * (self.T_in_g - self.T_dew_r)

        self.T_out_r            = self.T_dew_r + self.Q_superheat / (self.m_dot_r * cp_r_superheat)
        # Refrigerant density (superheated)
        self.AS_r.update(CP.PT_INPUTS, self.p_in_r, (self.T_in_g + self.T_dew_r) / 2.0)
        rho_superheat           = self.AS_r.rhomass()                           # [kg/m^3]
        # Refrigerant charge (superheated)
        self.Charge_r_superheat = w_superheat * self.V_r * rho_superheat

        # Pressure drop calculations for superheated refrigerant
        v_r                     = 1. / rho_superheat
        # Pressure gradient using Darcy friction factor
        dpdz_r                  = -self.f_r_superheat * v_r * self.G_r**2 / (2. * self.Dh_r)    # Pressure gradient
        self.DP_r_superheat     = dpdz_r * self.L * w_superheat

        # Temperature of "glycol" at the point where the refrigerant is at a quality of 1.0 [K]
        self.T_g_x              = self.T_in_g - self.Q_superheat / (self.cp_g * self.m_dot_g)

    def _TwoPhase_Forward(self, w_2phase=1.0, x_out_2phase=1.0):
        # Update outlet quality field [-]
        self.x_out_2phase                   = x_out_2phase
        # Heat transfer rate based on inlet quality [W]
        self.Q_2phase                       = self.m_dot_r * (self.x_out_2phase - self.x_in_r) * self.h_fg
        # Heat flux in 2phase section (for Shah correlation) [W/m^2]
        q_flux                              = self.Q_2phase / (w_2phase * self.A_r_wetted)

        h_r_2phase              = ShahEvaporation_Average(self.x_in_r, 1.0, self.AS_r, self.G_r, self.Dh_r, self.p_in_r,
                                                          q_flux, self.T_bubble_r, self.T_dew_r)
        self.h_r_2phase         = h_r_2phase * self.h_tp_tuning
        UA_2phase               = w_2phase / (1 / (self.h_g * self.A_g_wetted) + 1 / (self.h_r_2phase * self.A_r_wetted) + self.R_w)
        C_g                     = self.cp_g * self.m_dot_g
        Ntu_2phase              = UA_2phase / C_g

        # for Cr=0 and counter flow in two-phase:
        epsilon_2phase          = 1 - exp(-Ntu_2phase)
        Q_2phase_eNTU           = epsilon_2phase * C_g * (self.T_g_x - self.T_sat_r)

        rho_average             = TwoPhaseDensity(self.AS_r, self.x_in_r, x_out_2phase, self.T_dew_r, self.T_bubble_r, slipModel='Zivi')
        self.Charge_r_2phase    = rho_average * w_2phase * self.V_r

        # Frictional pressure drop component
        DP_frict                = LMPressureGradientAvg(self.x_in_r, x_out_2phase, self.AS_r, self.G_r, self.Dh_r, self.T_bubble_r, self.T_dew_r)\
                                  * w_2phase * self.L
        # Accelerational prssure drop component
        DP_accel                = AccelPressureDrop(self.x_in_r, x_out_2phase, self.AS_r, self.G_r, self.T_bubble_r, self.T_dew_r)\
                                  * w_2phase * self.L
        self.DP_r_2phase        = DP_frict + DP_accel

        if self.Verbosity > 4:
            print(Q_2phase_eNTU - self.Q_2phase)
        return Q_2phase_eNTU - self.Q_2phase


# -----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    TT                          = []
    QQ                          = []
    Q1                          = []
    # refrigerant Abstract State
    Ref_r                       = 'R290'
    Backend_r                   = 'HEOS'    # choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    AS_r                        = CP.AbstractState(Backend_r, Ref_r)
    # glycol Abstract State
    Ref_g                       = 'Water'
    Backend_g                   = 'HEOS'    # choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    AS_g                        = CP.AbstractState(Backend_g, Ref_g)
    for T_dew_evap in np.linspace(270, 290.4):
        T_dew_cond              = 317.73
        # T_dew_evap             = 285.42
        p_dew_cond              = PropsSI('P', 'T', T_dew_cond, 'Q', 1.0, Ref_r)
        h                       = PropsSI('H', 'T', T_dew_cond-7, 'P', p_dew_cond, Ref_r)
        params  = {
                'ID_i':         0.0278,         # inner tube, Internal Diameter (ID)
                'OD_i':         0.03415,        # inner tube, Outer Diameter (OD)
                'ID_o':         0.045,          # outer tube (annulus), Internal Diameter (ID)
                'L':            50,
                'm_dot_r':      0.040,
                'm_dot_g':      0.38,
                'h_in_r':       h,
                'p_in_r':       PropsSI('P', 'T', T_dew_evap, 'Q', 1.0, Ref_r),
                'p_in_g':       300000,         # p_in_g in Pa
                'T_in_g':       290.52,
                'AS_r':         AS_r,           # Abstract state of refrigerant
                'AS_g':         AS_g,           # Abstract state of glycol
                'Verbosity':    0,
                'Conductivity': 237,            # [W/m-K]
                'h_g_tuning':   1,
                'h_tp_tuning':  1,
                'DP_g_tuning':  1,
                'DP_r_tuning':  1
                }
        IHX                     = CoaxialHXClass(**params)
        IHX.Calculate()

        TT.append(T_dew_evap)
        QQ.append(IHX.h_r_2phase)               # IHX.Q
        Q1.append(IHX.h_r_superheat)

        print(IHX.Q)
