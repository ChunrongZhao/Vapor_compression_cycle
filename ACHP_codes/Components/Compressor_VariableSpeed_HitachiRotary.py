from __future__                                 import division, print_function, absolute_import
from math                                       import fabs
from CoolProp.CoolProp                          import PropsSI
from ACHP_codes.Correlations.OilPropLib         import *
import CoolProp                                 as CP


# -------------------------------------------------------------------------------------------------
class HitachiVariableSpeedRotaryCompressorClass():
    """
    Hitachi variable-speed rolling piston compressor with R401A.

    Required Parameters:

    ===========    ==========  ========================================================================
    Variable       Units       Description
    ===========    ==========  ========================================================================
    a_etav          --         A numpy-like list of compressor map coefficients for volumetric eff.
    a_etais         --         A numpy-like list of compressor map coefficients for isentropic eff.
    a_etaoi          --        A numpy-like list of compressor map coefficients for overall eff.
    Ref            N/A         A string representing the refrigerant
    Oil            N/A         A string representing the lubricant oil
    T_in_r         K           Refrigerant inlet temperature
    p_in_r         Pa          Refrigerant suction pressure (absolute)
    p_out_r        Pa          Refrigerant discharge pressure (absolute)
    V_disp         m^3         Displacement volume
    V_dot_ratio    --          Displacement Scale factor
    V_oil_sump     m^3         Total volume of oil sump inside the compressor shell
    shell_pressure N/A         A string defining the shell pressure of the compressor
    N              rpm         Compressor rotational speed
    ===========   ==========  ========================================================================
    """

    def __init__(self, **kwargs):
        # Load up the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)

    def Update(self, **kwargs):
        # Update the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)

    def OutputList(self):
        """
            Return a list of parameters for this component for further output

            It is a list of tuples, and each tuple is formed of items with indices:
                [0] Description of value
                [1] Units of value
                [2] The value itself
        """
        return [
            ('Heat Loss Fraction',              '-',            self.fp),
            ('Displacement scale factor',       '-',            self.V_dot_ratio),
            ('Power',                           'W',            self.W),
            ('Mass flow rate',                  'kg/s',         self.m_dot_r),
            ('Inlet Temperature',               'K',            self.T_in_r),
            ('Outlet Temperature',              'K',            self.T_out_r),
            ('Inlet Enthalpy',                  'J/kg',         self.h_in_r),
            ('Outlet Enthalpy',                 'J/kg',         self.h_out_r),
            ('Overall isentropic efficiency',   '-',            self.eta_oi),
            ('Pumped flow rate',                'm^3/s',        self.V_dot_pumped),
            ('Ambient heat loss',               'W',            self.Q_amb),
            ('Refrigerant change in oil sump',  'kg',           self.Charge)
         ]

    def Calculate(self):
        # AbstractState
        AS                              = self.AS

        # Reference molar mass of R410A
        MM_R410A                        = 72.5854                   # [kg/kmol]

        # Molar mass of working fluid
        MM                              = AS.molar_mass()           # [kg/mol]
        MM                              *= 1000                     # [kg/kmol]

        # Ambient conditions
        T_amb                           = self.T_amb

        # Compressor suction conditions
        p_in                            = self.p_in_r
        T_in                            = self.T_in_r

        AS.update(CP.PT_INPUTS, p_in, T_in)
        rho_in                          = AS.rhomass()              # [kg/m^3]
        h_in                            = AS.hmass()                # [J/kg]
        s_in                            = AS.smass()                # [J/kg-K]

        # Compressor discharge conditions
        p_out                           = self.p_out_r

        AS.update(CP.PSmass_INPUTS, p_out, s_in)
        h_out_is                        = AS.hmass()                # [J/kg]
        T_out_is                        = AS.T()                    # [K]

        # Reference rotational speed
        N_r                             = 60                        # [rps]

        # Rotational speed
        N                               = self.N / 60               # [rps]

        # Compressor displacement with scale factor
        V_disp                          = self.V_disp * self.V_dot_ratio

        # Isentropic enthalpy difference
        delta_h_is                      = h_out_is - h_in

        # Local copies of coefficients
        a_etav                          = self.a_etav
        a_etais                         = self.a_etais
        a_etaoi                         = self.a_etaoi

        # pi-groups isentropic efficiency
        pi2_etais                       = p_out / p_in
        pi3_etais                       = N_r / N
        pi4_etais                       = (N**3 * V_disp) / delta_h_is
        pi5_etais                       = delta_h_is * rho_in / p_in
        # The math.fabs() method returns the absolute value of a number, as a float.
        pi6_etais                       = 1 / fabs((T_in + T_out_is) / 2 - T_amb)
        pi7_etais                       = MM_R410A / MM

        # pi-groups volumetric efficiency
        pi2_etav                        = p_out / p_in
        pi3_etav                        = (rho_in / p_in)**1.5 * N**3 * V_disp
        pi4_etav                        = N_r / N
        pi5_etav                        = MM_R410A / MM

        # Expressions for volumetric efficiency, isentropic efficiency and overall isentropic efficiency
        eta_v   = a_etav[0] * pi2_etav**a_etav[1]   * pi4_etav**a_etav[2]   * pi5_etav**a_etav[3]
        eta_is  = a_etais[0]* pi2_etais**a_etais[1] * pi3_etais**a_etais[2] * pi4_etais**a_etais[3] * pi6_etais**a_etais[4]
        eta_oi  = a_etaoi[0]* pi2_etais**a_etaoi[1] * pi3_etais**a_etaoi[2] * pi6_etais**a_etaoi[3] * pi7_etais**a_etaoi[4]

        # Refrigerant mass flow rate
        m_dot                           = eta_v * rho_in * V_disp * N       # [kg/s]

        # Discharge enthalpy
        h_out                           = h_in + (h_out_is - h_in) / eta_is

        # Discharge temperature and specific entropy
        AS.update(CP.HmassP_INPUTS, h_out, p_out)
        T_out                           = AS.T()                            # [K]
        s_out                           = AS.smass()                        # [J/kg-K]

        # Power input
        W_dot_el                        = m_dot * (h_out_is - h_in) / eta_oi

        # Heat loss to ambient
        self.Q_amb                      = m_dot * (h_out - h_in) - W_dot_el

        # Outputs
        self.h_in_r                     = h_in
        self.s_in_r                     = s_in
        self.h_out_r                    = h_out
        self.T_out_r                    = T_out
        self.s_out_r                    = s_out
        self.m_dot_r                    = m_dot
        self.W                          = W_dot_el
        self.eta_v                      = eta_v
        self.eta_is                     = eta_is
        self.eta_oi                     = eta_oi
        self.fp                         = -self.Q_amb / W_dot_el
        self.CycleEnergyIn              = W_dot_el * (1 - self.fp)
        self.V_dot_pumped               = m_dot / rho_in

        # Estimate refrigerant dissolved in the oil sump
        T_ave                           = (T_in + self.T_out_r) / 2
        if self.shell_pressure == 'high-pressure':
            p_shell                     = p_out
        elif self.shell_pressure == 'low-pressure':
            p_shell                     = p_in

        # Solubility fraction
        self.x_Ref, error               = Solubility_Ref_in_Liq(self.Ref, self.Oil, T_ave, p_shell/1000)

        AS.update(CP.PT_INPUTS, p_shell, T_ave)
        rho_shell                       = AS.rhomass()                      # [kg/m^3]

        rho_mass_oil                    = rho_oil(self.Oil, T_ave-273.15)
        self.m_oil                      = self.V_oil_sump * rho_mass_oil

        # Amount of refrigerant dissolved in the oil sump
        self.Charge                     = self.m_oil * self.x_Ref / (1 - self.x_Ref)

        if self.Verbosity > 0:
            print('Compressor pin [Pa]:',               p_in)
            print('Compressor pout [Pa]:',              p_out)
            print('Compressor Tout [C]:',               T_out-273.15)
            print('Compressor eta_oi [-]:',             self.eta_oi)
            print('Compressor eta_is [-]:',             self.eta_is)
            print('Compressor eta_v [-]:',              self.eta_v)
            print('Compressor mass flow rate [g/s]:',   self.m_dot_r*1000)
            print('Compressor Power [W]:',              self.W)


# -------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # Abstract State
    Ref                                 = 'R410a'
    Backend                             = 'HEOS'    # choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    AS                                  = CP.AbstractState(Backend, Ref)

    # Compressor inlet temperature
    T_in_r                              = 19.9                          # oC

    # Compressor inlet and outlet pressures
    p_in_r                              = (157.6*6.89476)               # psia to kPa
    p_out_r                             = (398.6*6.89476)               # psia to kPa

    for i in range(1):
        kwds    = {
                'AS':                   AS,                             # Abstract state
                'Ref':                  Ref,
                'T_amb':                35 + 273.15,
                'T_in_r':               T_in_r + 273.15,                # K
                'p_in_r':               p_in_r*1000,                    # Pa
                'p_out_r':              p_out_r*1000,                   # Pa
                'V_disp':               47e-6,                          # Displacement volume, m^3
                'V_dot_ratio':          1.0,                            # Displacement Scale factor
                'N':                    3600,                           # Rotational speed, rpm
                'a_etav':               [1., -0.32094114, 0.00170576, -0.03695206],             # Volumetric eff. coeff.
                'a_etais':              [1., 0.24148598,  0.37491376, 0.07996134, -0.03366503], # Isentropic eff. coeff.
                'a_etaoi':              [1., -0.01360404, -0.07435644, 0.18584579, 0.63589724], # Overall eff. coeff.
                'shell_pressure':       'low-pressure',
                'Oil':                  'POE32',
                'V_oil_sump':           0.0,
                'Verbosity':            0.0
                }
        Comp                            = HitachiVariableSpeedRotaryCompressorClass(**kwds)
        Comp.Calculate()
        print('Power:',                                 Comp.W,         'W')
        print('Mass flow rate:',                        Comp.m_dot_r,   'kg/s')
        print('Volumetric efficiency:',                 Comp.eta_v,     '-')
        print('Isentropic efficiency:',                 Comp.eta_is,    '-')
        print('Overall efficiency:',                    Comp.eta_oi,    '-')
        print('Refrigerant dissolved in oil sump:',     Comp.Charge,    'kg')
