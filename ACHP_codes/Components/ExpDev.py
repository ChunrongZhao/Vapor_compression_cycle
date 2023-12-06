from __future__                             import division, absolute_import, print_function
from CoolProp.CoolProp                      import PropsSI
import CoolProp                             as CP

from scipy.optimize                         import brentq
from math                                   import pi, exp, log, sqrt, tan, cos, sin, pow, atan
from ACHP_codes.ACHP_Tools.convert_units    import cms2gpm, psi2kPa, C2K, in2m


# ------------------------------------------------------------------------------
class ExpDevClass():
    """
    Expansion devices models
    """
    def __init__(self, **kwargs):
        # Load the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)

    def Update(self, **kwargs):
        # Update the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)

    def OutputList(self):               # TODO: fix this list of outputs
        """
            Return a list of parameters for this component for further output

            It is a list of tuples, and each tuple is formed of items:
                [0] Description of value
                [1] Units of value
                [2] The value itself
        """

        return [
            ('Expansion Device Type',       '-',        self.ExpType),
            ('Upstream Pressure',           'Pa',       self.p_in_r),
            ('Upstream Enthalpy',           'j/kg',     self.h_in_r),
            ('Downstream Pressure',         'Pa',       self.p_out_r),
            ('Downstream Quality',          '-',        self.x_out_r),
            ('Mass flow rate',              'kg/s',     self.m_dot_r)]

    def Initialize(self):

        # AbstractState
        assert hasattr(self, 'AS'), 'Please specify the Abstract State'

        # If the user doesn't include the ExpType, fail
        assert hasattr(self, 'ExpType'), 'Please specify the type of the expansion device'

    def Calculate(self):

        # Initialize
        self.Initialize()
        # AbstractState
        AS                                      = self.AS

        if self.ExpType == 'Ideal':
            # ===================================================================
            # No information about expansion device is given
            # ===================================================================
            # inlet state
            if self.p_in_r > AS.p_critical():                           # Supercritical
                AS.update(CP.HmassP_INPUTS, self.h_in_r, self.p_in_r)
                self.s_in_r                     = AS.smass()            # [J/kg-K]
                self.T_in_r                     = AS.T()                # [K]
            else:                                                       # other refrigerants
                AS.update(CP.PQ_INPUTS, self.p_in_r, 0.0)
                T_bubble_in                     = AS.T()                # [K]
                h_l_in                          = AS.hmass()            # [J/kg]
                s_l_in                          = AS.smass()            # [J/kg-K]
                AS.update(CP.PQ_INPUTS, self.p_in_r, 1.0)
                T_dew_in                        = AS.T()                # [K]
                h_v_in                          = AS.hmass()            # [J/kg]
                s_v_in                          = AS.smass()            # [J/kg-K]

                self.x_in_r         = (self.h_in_r - h_l_in) / (h_v_in - h_l_in)

                if self.x_in_r > 0.999:
                    print("ExpDev :: Upstream state in the expansion device is superheated")
                    raise
                if self.x_in_r > 0.0 and self.x_in_r < 1.0:             # two-phase state at the inlet
                    self.s_in_r     = self.x_in_r * s_v_in + (1 - self.x_in_r) * s_l_in         # [J/kg-K]
                    self.T_in_r     = self.x_in_r * T_dew_in + (1 - self.x_in_r) * T_bubble_in  # [K]
                else:                                                   # liquid state at the inlet
                    AS.update(CP.HmassP_INPUTS, self.h_in_r, self.p_in_r)
                    self.s_in_r                 = AS.smass()            # [J/kg-K]
                    self.T_in_r                 = AS.T()                # [K]

            # outlet state (assume h = constant)
            self.h_out_r                        = self.h_in_r           # [J/kg]

            AS.update(CP.PQ_INPUTS, self.p_out_r, 0.0)
            T_bubble_out                        = AS.T()                # [K]
            h_l_out                             = AS.hmass()            # [J/kg]
            s_l_out                             = AS.smass()            # [J/kg-K]
            AS.update(CP.PQ_INPUTS, self.p_out_r, 1.0)
            T_dew_out                           = AS.T()                # [K]
            h_v_out                             = AS.hmass()            # [J/kg]
            s_v_out                             = AS.smass()            # [J/kg-K]

            # outlet state (two-phase)
            self.x_out_r = (self.h_out_r - h_l_out) / (h_v_out - h_l_out)                   # [-]
            self.T_out_r = self.x_out_r * T_dew_out + (1 - self.x_out_r) * T_bubble_out     # [K]
            self.s_out_r = self.x_out_r * s_v_out + (1 - self.x_out_r) * s_l_out            # [J/kg-K]

            # mass flow rate
            self.m_dot_r                        = 'N/A'

            # heat losses
            self.Q_amb                          = 0.0                   # [W]

        if self.ExpType == 'Linear-TXV':
            # ===================================================================
            # Global Linear TxV model from Haorong Li paper (2004)
            # paper title: Modeling Adjustable throat-Area Expansion Valves
            # ===================================================================
            D                                   = self.D                # inside diameter [m]
            T_sh_static                         = self.T_sh_static      # [K]
            T_sh_max                            = self.T_sh_max         # [K]
            Adj                                 = self.Adj              # [-]
            C                                   = self.C                # [m^2/K]

            T_sup                               = self.T_sup            # superheat value (user defined)

            P_up                                = self.p_in_r
            P_down                              = self.p_out_r

            A                                   = (T_sup - T_sh_static)
            if A > T_sh_max:
                A                               = T_sh_max

            # thermodynamic properties
            AS.update(CP.PQ_INPUTS, self.p_in_r, 0.0)
            T_bubble_in                         = AS.T()                # [K]
            h_l_in                              = AS.hmass()            # [J/kg]
            s_l_in                              = AS.smass()            # [J/kg-K]
            rho_l_in                            = AS.rhomass()          # [kg/m^3]
            AS.update(CP.PQ_INPUTS, self.p_in_r, 1.0)
            T_dew_in                            = AS.T()                # [K]
            h_v_in                              = AS.hmass()            # [J/kg]
            s_v_in                              = AS.smass()            # [J/kg-K]
            rho_v_in                            = AS.rhomass()          # [kg/m^3]

            # inlet state
            self.x_in_r                         = (self.h_in_r - h_l_in) / (h_v_in - h_l_in)
            if self.x_in_r > 0.999:
                print("ExpDev :: Upstream state in the expansion device is superheated")
                raise
            if self.x_in_r > 0.0 and self.x_in_r < 1.0:
                # 2phase upstream state
                print("ExpDev :: Upstream state in the expansion device is 2-phase")
                self.s_in_r = self.x_in_r * s_v_in + (1 - self.x_in_r) * s_l_in             # [J/kg-K]
                self.T_in_r = self.x_in_r * T_dew_in + (1 - self.x_in_r) * T_bubble_in      # [K]
            else:                                               # liquid state at the inlet
                AS.update(CP.HmassP_INPUTS, self.h_in_r, self.p_in_r)
                self.s_in_r = AS.smass()                                # [J/kg-K]
                self.T_in_r = AS.T()                                    # [K]

            # upstream saturated liquid density
            rho_up                              = rho_l_in

            # calculate mass flow rate
            m_dot_r                             = C * A * pow(rho_up * (P_up-P_down), 0.5)

            # adjust the mass flow rate via adjustment factor related with geometry (tuning factor)
            self.m_dot_r                        = m_dot_r * Adj

            # outlet state (assume h = constant)
            self.h_out_r                        = self.h_in_r           # [J/kg]

            AS.update(CP.PQ_INPUTS, self.p_out_r, 0.0)
            T_bubble_out                        = AS.T()                # [K]
            h_l_out                             = AS.hmass()            # [J/kg]
            s_l_out                             = AS.smass()            # [J/kg-K]
            AS.update(CP.PQ_INPUTS, self.p_out_r, 1.0)
            T_dew_out                           = AS.T()                # [K]
            h_v_out                             = AS.hmass()            # [J/kg]
            s_v_out                             = AS.smass()            # [J/kg-K]

            # outlet state (two-phase)
            self.x_out_r            = (self.h_out_r - h_l_out) / (h_v_out - h_l_out)                # [-]
            self.T_out_r            = self.x_out_r * T_dew_out + (1 - self.x_out_r) * T_bubble_out  # [K]
            self.s_out_r            = self.x_out_r * s_v_out + (1 - self.x_out_r) * s_l_out         # [J/kg-K]

            # heat losses
            self.Q_amb                          = 0.0                   # [W]

        if self.ExpType == 'Nonlinear-TXV':
            # ===================================================================
            # Nonlinear TxV model from Haorong Li paper (2004)
            # paper title: Modeling Adjustable throat-Area Expansion Valves
            # ===================================================================
            D                                   = self.D                    # inside diameter [m]
            T_sh_static                         = self.T_sh_static          # [K]
            T_sh_max                            = self.T_sh_max             # [K]
            Adj                                 = self.Adj                  # [-]
            C                                   = self.C                    # [m^2/K]

            T_sup                               = self.T_sup                # superheat value (user defined)

            P_up                                = self.p_in_r
            P_down                              = self.p_out_r

            A                                   = (T_sup - T_sh_static) / T_sh_max
            if A > 1:
                A                               = 1

            # thermodynamic properties
            AS.update(CP.PQ_INPUTS, self.p_in_r, 0.0)
            T_bubble_in                         = AS.T()                    # [K]
            h_l_in                              = AS.hmass()                # [J/kg]
            s_l_in                              = AS.smass()                # [J/kg-K]
            rho_l_in                            = AS.rhomass()              # [kg/m^3]
            AS.update(CP.PQ_INPUTS, self.p_in_r, 1.0)
            T_dew_in                            = AS.T()                    # [K]
            h_v_in                              = AS.hmass()                # [J/kg]
            s_v_in                              = AS.smass()                # [J/kg-K]
            rho_v_in                            = AS.rhomass()              # [kg/m^3]

            # inlet state
            self.x_in_r                         = (self.h_in_r - h_l_in) / (h_v_in - h_l_in)

            if self.x_in_r > 0.999:
                print("ExpDev :: Upstream state in the expansion device is superheated")
                raise
            if self.x_in_r > 0.0 and self.x_in_r < 1.0:
                # 2phase upstream state
                print("ExpDev :: Upstream state in the expansion device is 2-phase")
                self.s_in_r     = self.x_in_r * s_v_in + (1 - self.x_in_r) * s_l_in         # [J/kg-K]
                self.T_in_r     = self.x_in_r * T_dew_in + (1 - self.x_in_r) * T_bubble_in  # [K]
            else:                                                           # liquid state at the inlet
                AS.update(CP.HmassP_INPUTS, self.h_in_r, self.p_in_r)
                self.s_in_r                     = AS.smass()                # [J/kg-K]
                self.T_in_r                     = AS.T()                    # [K]

            # upstream saturated liquid density
            rho_up                              = rho_l_in

            # calculate mass flow rate
            m_dot_r                             = C * (2 * A - A * A) * pow(rho_up * (P_up - P_down), 0.5)

            # adjust the mass flow rate via adjustment factor related with geometry (tuning factor)
            self.m_dot_r                        = m_dot_r * Adj

            # outlet state (assume h = constant)
            self.h_out_r                        = self.h_in_r               # [J/kg]

            AS.update(CP.PQ_INPUTS, self.p_out_r, 0.0)
            T_bubble_out                        = AS.T()                    # [K]
            h_l_out                             = AS.hmass()                # [J/kg]
            s_l_out                             = AS.smass()                # [J/kg-K]
            AS.update(CP.PQ_INPUTS, self.p_out_r, 1.0)
            T_dew_out                           = AS.T()                    # [K]
            h_v_out                             = AS.hmass()                # [J/kg]
            s_v_out                             = AS.smass()                # [J/kg-K]

            # outlet state (two-phase)
            self.x_out_r        = (self.h_out_r - h_l_out) / (h_v_out - h_l_out)                # [-]
            self.T_out_r        = self.x_out_r * T_dew_out + (1 - self.x_out_r) * T_bubble_out  # [K]
            self.s_out_r        = self.x_out_r * s_v_out + (1 - self.x_out_r) * s_l_out         # [J/kg-K]

            # heat losses
            self.Q_amb                          = 0.0                       # [W]

        if self.ExpType == 'Short-tube':
            # ===================================================================
            # Short tube expansion from Payne and O'Neal (2004)
            # paper title: A Mass Flowrate Correlation for Refrigerants and Refrigerant Mixtures, Journal of HVAC
            # based on empirical dimensionless PI correlation, recommended for R-12, R-134a, R-502, R-22, R-407C, and R-410A
            # ===================================================================
            D                                   = self.D                        # inside diameter of the short-tube[m]
            L                                   = self.L                        # length of the short-tube [m]
            Adj                                 = self.Adj                      # adjusting the inside diameter [-];
            L_c                                 = self.L_c                      # chamfered length [m]
            Ang_c                               = self.Ang_c                    # chamfered angle [degree]
            BranNum                             = int(self.BranNum)             # Number of Paralelled expansion devices

            A_s                                 = pi / 4 * D * D

            P_up                                = self.p_in_r
            P_down                              = self.p_out_r

            # critical point of refrigerant
            P_c                                 = AS.p_critical()               # [Pa]
            T_c                                 = AS.T_critical()               # [K]

            # orifice adjustment parameter
            C_c                                 = Adj

            # thermodynamic properties
            AS.update(CP.PQ_INPUTS, self.p_in_r, 0.0)
            T_bubble_in                         = AS.T()                        # [K]
            h_l_in                              = AS.hmass()                    # [J/kg]
            s_l_in                              = AS.smass()                    # [J/kg-K]
            rho_l_in                            = AS.rhomass()                  # [kg/m^3]
            AS.update(CP.PQ_INPUTS, self.p_in_r, 1.0)
            T_dew_in                            = AS.T()                        # [K]
            h_v_in                              = AS.hmass()                    # [J/kg]
            s_v_in                              = AS.smass()                    # [J/kg-K]
            rho_v_in                            = AS.rhomass()                  # [kg/m^3]

            # inlet state
            self.x_in_r                         = (self.h_in_r - h_l_in) / (h_v_in - h_l_in)
            if self.x_in_r > 0.999:
                print("ExpDev :: Upstream state in the expansion device is superheated")
                raise
            if self.x_in_r > 0.0 and self.x_in_r < 1.0:                         # two-phase state at the inlet
                self.s_in_r         = self.x_in_r * s_v_in + (1 - self.x_in_r) * s_l_in         # [J/kg-K]
                self.T_in_r         = self.x_in_r * T_dew_in + (1 - self.x_in_r) * T_bubble_in  # [K]
            else:                                                               # liquid state at the inlet
                AS.update(CP.HmassP_INPUTS, self.h_in_r, self.p_in_r)
                self.s_in_r                     = AS.smass()                    # [J/kg-K]
                self.T_in_r                     = AS.T()                        # [K]

            AS.update(CP.QT_INPUTS, 0, self.T_in_r)
            P_sat                               = AS.p()    # P_sat corresponding to upstream temperature (liquid saturation pressure) [Pa]
            T_sat                               = T_bubble_in
            T_sub                               = T_sat - self.T_in_r

            # upstream saturated liquid density
            AS.update(CP.PQ_INPUTS, P_sat, 0.0)
            rho_f                               = AS.rhomass()                  # [kg/m^3]
            # upstream saturated gas density
            AS.update(CP.PQ_INPUTS, P_sat, 1.0)
            rho_g                               = AS.rhomass()                  # [kg/m^3]

            # non-dimensional groups
            pi_3                                = (P_up - P_sat) / P_c
            pi_6                                = rho_g / rho_f
            pi_9                                = T_sub / T_c
            pi_10                               = L / D

            # coefficients
            a1                                  = 3.8811e-1
            a2                                  = 1.1427e1
            a3                                  = -1.4194e1
            a4                                  = 1.0703e0
            a5                                  = -9.1928e-2
            a6                                  = 2.1425e1
            a7                                  = -5.8195e2

            # single-phase flow rate
            pi_1            = (a1 + a2 * pi_3 + a3 * pi_9 + a4 * pi_6 + a5 * log(pi_10)) / (1 + a6 * pi_3 + a7 * pi_9 * pi_9)
            G               = pi_1 * pow((rho_f * P_c), 0.5)

            # mass flow rate
            m_dot_r                             = G * A_s

            if self.x_in_r < 0.000001:                                          # subcooled upstream state
                m_dot_r                         = m_dot_r

            else:                                                               # two-phase upstream state
                x_up                            = self.x_in_r
                rho_mup                         = 1 / ((1 - x_up) / rho_f + x_up / rho_g)

                # non-dimensional groups
                tp6                             = rho_mup / rho_f
                tp35                            = (P_c - P_sat) / P_c
                tp32                            = (P_c - P_up) / P_c
                tp27                            = L / D
                tp34                            = x_up / (1 - x_up) * pow((rho_f / rho_g), 0.5)
                tp28                            = P_up / P_c

                # coefficients
                b1                              = 1.1831e0
                b2                              = -1.468e0
                b3                              = -1.5285e-1
                b4                              = -1.4639e1
                b5                              = 9.8401e0
                b6                              = -1.9798e-2
                b7                              = -1.5348e0
                b8                              = -2.0533e0
                b9                              = -1.7195e1

                numer       = (b1 * tp6 + b2 * pow(tp6, 2.0) + b3 * pow(log(tp6), 2.0) + b4 * pow(log(tp35), 2.0)
                               + b5 * pow(log(tp32), 2.0) + b6 * pow(log(tp27), 2.0))
                deno        = (1 + b7 * tp6 + b8 * tp34 + b9 * pow(tp28, 3.0))
                C_tp        = numer / deno                                      # two-phase flow rate adjustment

                if C_tp > 1:
                    C_tp                        = 1                             # since C_tp>1 is not right

                # correct the mass flow rate by two-phase entrance
                m_dot_r                         = m_dot_r * C_tp

            # adjust the mass flow rate via adjustment factor related with geometry (tuning factor)
            self.m_dot_r                        = m_dot_r * C_c

            # multiply mass flow rate by the number of parallel branches
            if BranNum == 0:
                self.m_dot_r                    = self.m_dot_r
            else:
                self.m_dot_r                    = self.m_dot_r * BranNum

            # outlet state (assume h = constant)
            self.h_out_r                        = self.h_in_r                   # [J/kg]

            AS.update(CP.PQ_INPUTS, self.p_out_r, 0.0)
            T_bubble_out                        = AS.T()                        # [K]
            h_l_out                             = AS.hmass()                    # [J/kg]
            s_l_out                             = AS.smass()                    # [J/kg-K]
            AS.update(CP.PQ_INPUTS, self.p_out_r, 1.0)
            T_dew_out                           = AS.T()                        # [K]
            h_v_out                             = AS.hmass()                    # [J/kg]
            s_v_out                             = AS.smass()                    # [J/kg-K]

            # outlet state (two-phase)
            self.x_out_r    = (self.h_out_r - h_l_out) / (h_v_out - h_l_out)                # [-]
            self.T_out_r    = self.x_out_r * T_dew_out + (1 - self.x_out_r) * T_bubble_out  # [K]
            self.s_out_r    = self.x_out_r * s_v_out + (1 - self.x_out_r) * s_l_out         # [J/kg-K]

        if self.ExpType == 'Expander':
            # ===================================================================
            # General expander with given isentropic efficiency
            # ===================================================================
            # inlet state
            if self.p_in_r > AS.p_critical():                           # Supercritical
                AS.update(CP.HmassP_INPUTS, self.h_in_r, self.p_in_r)
                self.s_in_r                     = AS.smass()            # [J/kg-K]
                self.T_in_r                     = AS.T()                # [K]
            else:                                                       # other refrigerants
                AS.update(CP.PQ_INPUTS, self.p_in_r, 0.0)
                T_bubble_in                     = AS.T()                # [K]
                h_l_in                          = AS.hmass()            # [J/kg]
                s_l_in                          = AS.smass()            # [J/kg-K]
                AS.update(CP.PQ_INPUTS, self.p_in_r, 1.0)
                T_dew_in                        = AS.T()                # [K]
                h_v_in                          = AS.hmass()            # [J/kg]
                s_v_in                          = AS.smass()            # [J/kg-K]

                self.x_in_r                     = (self.h_in_r - h_l_in) / (h_v_in - h_l_in)
                if self.x_in_r > 0.999:
                    print("ExpDev :: Upstream state in the expansion device is superheated")
                    raise
                if self.x_in_r > 0.0 and self.x_in_r < 1.0:             # two-phase state at the inlet
                    self.s_in_r     = self.x_in_r * s_v_in + (1 - self.x_in_r) * s_l_in             # [J/kg-K]
                    self.T_in_r     = self.x_in_r * T_dew_in + (1 - self.x_in_r) * T_bubble_in      # [K]
                else:                                                   # liquid state at the inlet
                    AS.update(CP.HmassP_INPUTS, self.h_in_r, self.p_in_r)
                    self.s_in_r                 = AS.smass()            # [J/kg-K]
                    self.T_in_r                 = AS.T()                # [K]

            # isentropic outlet state
            AS.update(CP.PSmass_INPUTS, self.p_out_r, self.s_in_r)
            self.h_out_s_r                      = AS.hmass()            # [J/kg]

            # outlet state (assume eta_is = given)
            self.h_out_r            = self.h_in_r - self.eta_is * (self.h_in_r - self.h_out_s_r)    # [J/kg]

            AS.update(CP.PQ_INPUTS, self.p_out_r, 0.0)
            T_bubble_out                        = AS.T()                # [K]
            h_l_out                             = AS.hmass()            # [J/kg]
            s_l_out                             = AS.smass()            # [J/kg-K]
            AS.update(CP.PQ_INPUTS, self.p_out_r, 1.0)
            T_dew_out                           = AS.T()                # [K]
            h_v_out                             = AS.hmass()            # [J/kg]
            s_v_out                             = AS.smass()            # [J/kg-K]

            # outlet state (two-phase)
            self.x_out_r        = (self.h_out_r - h_l_out) / (h_v_out - h_l_out)                    # [-]
            self.T_out_r        = self.x_out_r * T_dew_out + (1 - self.x_out_r) * T_bubble_out      # [K]
            self.s_out_r        = self.x_out_r * s_v_out + (1 - self.x_out_r) * s_l_out             # [J/kg-K]

            # adjust the mass flow rate via adjustment factor related with geometry (tuning factor)
            # TODO: need to add a mass flow model
            m_dot_r                             = self.m_dot
            self.m_dot_r                        = m_dot_r * self.C_exp

            # heat losses
            self.Q_amb                          = 0.0                   # [W]


# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # Abstract State
    Ref                         = 'R410A'
    Backend                     = 'HEOS'        # choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    AS                          = CP.AbstractState(Backend, Ref)

    print('Example for Ideal expansion device')
    params  = {
        'AS':                       AS,
        'ExpType':                  'Ideal',                                                                    # expansion device type
        'p_in_r':                   PropsSI('P', 'T', 60+273.15, 'Q', 0, Ref),                                  # upstream pressure
        'h_in_r':                   PropsSI('H', 'P', PropsSI('P', 'T', 60+273.15, 'Q', 0, Ref), 'Q', 0, Ref),  # upstream enthalpy
        'p_out_r':                  PropsSI('P', 'T', 10+273.15, 'Q', 0, Ref),                                  # downstream pressure
    }
    Exp                         = ExpDevClass(**params)
    Exp.Calculate()
    print(Exp.OutputList())
    print()

    # ------------------------------------------------------------------------------
    print('Example for Linear expansion device')
    params  = {
        'AS':                       AS,
        'ExpType':                  'Linear-TXV',                           # expansion device type
        'T_sh_static':              4,                                      # static superheat
        'T_sh_max':                 6,                                      # maximum superheat
        'D':                        0.0006604,                              # inside diameter [m]
        'C':                        1.2656e-6,                              # constant from manufacturer [m^2/K]
        'Adj':                      0.7630,                                 # Adjust the diameter (tuning factor)
        'T_sup':                    5,                                      # superheat value (user defined)
        'p_in_r':                   PropsSI('P', 'T', 60+273.15, 'Q', 0, Ref),                                  # upstream pressure
        'h_in_r':                   PropsSI('H', 'P', PropsSI('P', 'T', 60+273.15, 'Q', 0, Ref), 'Q', 0, Ref),  # upstream enthalpy
        'p_out_r':                  PropsSI('P', 'T', 10+273.15, 'Q', 0, Ref),                                  # downstream pressure
    }
    Exp                         = ExpDevClass(**params)
    Exp.Calculate()
    print('T_out =',                Exp.T_out_r,        '[K]')
    print('h_out =',                Exp.h_out_r,        '[J/kg]')
    print('x_out =',                Exp.x_out_r,        '[-]')
    print('m_dot_r =',              Exp.m_dot_r,        '[kg/s]')
    print()

    # ---------------------------------------------------------------------------------------------
    print('Example for short-tube expansion device')
    params  = {
        'AS':                       AS,
        'ExpType':                  'Short-tube',                           # expansion device type
        'D':                        0.0006604,                              # inside diameter [m]
        'L':                        0.0052324,                              # length of short-tube[m]
        'L_c':                      0.0001524,                              # chamfered length [m] (P.S. not included in the solver yet)
        'Ang_c':                    45,                                     # chamfered angle [degree] (P.S. not included in the solver yet)
        'BranNum':                  12,                                     # Number of Paralelled short-tubes (0 -- default for 1 short-tube only)
        'Adj':                      1.094,                                  # Adjust the diameter (tuning factor)
        'p_in_r':                   PropsSI('P', 'T', 60+273.15, 'Q', 0, Ref),                                  # upstream pressure
        'h_in_r':                   PropsSI('H', 'P', PropsSI('P', 'T', 60+273.15, 'Q', 0, Ref), 'Q', 0, Ref),  # upstream enthalpy
        'p_out_r':                  PropsSI('P', 'T', 10+273.15, 'Q', 0, Ref),                                  # downstream pressure
    }
    Exp.Update(**params)
    Exp.Calculate()

    print('T_out =',                Exp.T_out_r,        '[K]')
    print('h_out =',                Exp.h_out_r,        '[J/kg]')
    print('x_out =',                Exp.x_out_r,        '[-]')
    print('m_dot_r =',              Exp.m_dot_r,        '[kg/s]')
    print()

    # ---------------------------------------------------------------------------------------------
    print('Example for expander device')
    params  = {
        'AS':                       AS,
        'ExpType':                  'Expander',                             # expansion device type
        'eta_is':                   0.8,                                    # isentropic efficiency [-]
        'C_exp':                    1,                                      # flow factor [-]
        'm_dot':                    0.01,                                   # mass flow rate [kg/s]
        'p_in_r':                   PropsSI('P', 'T', 60+273.15, 'Q', 0, Ref),                                  # upstream pressure
        'h_in_r':                   PropsSI('H', 'P', PropsSI('P', 'T', 60+273.15, 'Q', 0, Ref), 'Q', 0, Ref),  # upstream enthalpy
        'p_out_r':                  PropsSI('P', 'T', 10+273.15, 'Q', 0, Ref),                                  # downstream pressure
    }
    Exp.Update(**params)
    Exp.Calculate()

    print('T_out =',                Exp.T_out_r,        '[K]')
    print('h_out =',                Exp.h_out_r,        '[J/kg]')
    print('x_out =',                Exp.x_out_r,        '[-]')
    print('h_out_s =',              Exp.h_out_s_r,      '[J/kg]')
    print('m_dot_r =',              Exp.m_dot_r,        '[kg/s]')
    print()
