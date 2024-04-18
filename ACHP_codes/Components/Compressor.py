from __future__                         import division, print_function, absolute_import
from CoolProp.CoolProp                  import PropsSI
from ACHP_codes.Correlations.OilPropLib import *
import CoolProp                         as CP


# -------------------------------------------------------------------------------------------
class CompressorClass():
    """
    Compressor Model based on 10-coefficient Model from `ANSI/AHRI standard 540
    <http://www.ahrinet.org/App_Content/ahri/files/standards%20pdfs/ANSI%20standards%20pdfs/ANSI-ARI-540-2004%20latest.pdf>`_

    Required Parameters:

    ===========    ==========  ========================================================================
    Variable       Units       Description
    ===========    ==========  ========================================================================
    M              Ibm/hr      A numpy-like list of compressor map coefficients for mass flow
    P              Watts       A numpy-like list of compressor map coefficients for electrical power
    Ref            N/A         A string representing the refrigerant
    Oil            N/A         A string representing the lubricant oil
    T_in_r         K           Refrigerant inlet temperature
    p_in_r         Pa          Refrigerant suction pressure (absolute)
    p_out_r        Pa          Refrigerant discharge pressure (absolute)
    fp             --          Fraction of electrical power lost as heat to ambient
    V_dot_ratio    --          Displacement Scale factor
    V_oil_sump     m^3         Total volume of oil sump inside the compressor shell
    shell_pressure N/A         A string defining the shell pressure of the compressor
    ===========   ==========  ========================================================================

    All variables are of double-type unless otherwise specified

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
            ('M1',                          '-',        self.M[0]),
            ('M2',                          '-',        self.M[1]),
            ('M3',                          '-',        self.M[2]),
            ('M4',                          '-',        self.M[3]),
            ('M5',                          '-',        self.M[4]),
            ('M6',                          '-',        self.M[5]),
            ('M7',                          '-',        self.M[6]),
            ('M8',                          '-',        self.M[7]),
            ('M9',                          '-',        self.M[8]),
            ('M10',                         '-',        self.M[9]),
            ('P1',                          '-',        self.P[0]),
            ('P2',                          '-',        self.P[1]),
            ('P3',                          '-',        self.P[2]),
            ('P4',                          '-',        self.P[3]),
            ('P5',                          '-',        self.P[4]),
            ('P6',                          '-',        self.P[5]),
            ('P7',                          '-',        self.P[6]),
            ('P8',                          '-',        self.P[7]),
            ('P9',                          '-',        self.P[8]),
            ('P10',                         '-',        self.P[9]),
            ('Heat Loss Fraction',          '-',        self.fp),
            ('Displacement scale factor',   '-',        self.V_dot_ratio),
            ('Power',                       'W',        self.W),
            ('Mass flow rate',              'kg/s',     self.m_dot_r),
            ('Inlet Temperature',           'K',        self.T_in_r),
            ('Outlet Temperature',          'K',        self.T_out_r),
            ('Inlet Enthalpy',              'J/kg',     self.h_in_r),
            ('Outlet Enthalpy',             'J/kg',     self.h_out_r),
            ('Overall isentropic efficiency','-',       self.eta_oi),
            ('Pumped flow rate',            'm^3/s',    self.V_dot_pumped),
            ('Ambient heat loss',           'W',        self.Q_amb),
            ('Refrigerant change in oil sump','kg',     self.Charge)
         ]

    def Calculate(self):
        # AbstractState
        AS                                  = self.AS

        # Local copies of coefficients
        P                                   = self.P
        M                                   = self.M

        # Calculate suction superheat and dew temperatures
        AS.update(CP.PQ_INPUTS, self.p_in_r, 1.0)
        self.T_sat_s_K                      = AS.T()                    # [K]

        AS.update(CP.PQ_INPUTS, self.p_out_r, 1.0)
        self.T_sat_d_K                      = AS.T()                    # [K]
        self.DT_sh_K                        = self.T_in_r-self.T_sat_s_K

        # Convert saturation temperatures in K to F
        T_sat_s                             = self.T_sat_s_K * 9/5 - 459.67
        T_sat_d                             = self.T_sat_d_K * 9/5 - 459.67

        # Apply the 10 coefficient ARI map to saturation temps in F
        # p=P[0][9]
        power_map       = P[0] + P[1] * T_sat_s + P[2] * T_sat_d + P[3] * T_sat_s**2 + P[4] * T_sat_s * T_sat_d + P[5] * T_sat_d**2 \
                          + P[6] * T_sat_s**3 + P[7] * T_sat_d * T_sat_s**2 + P[8] * T_sat_d**2*T_sat_s + P[9] * T_sat_d**3
        m_dot_map       = M[0] + M[1] * T_sat_s + M[2] * T_sat_d + M[3] * T_sat_s**2 + M[4] * T_sat_s * T_sat_d + M[5] * T_sat_d**2 \
                          + M[6] * T_sat_s**3 + M[7] * T_sat_d * T_sat_s**2 + M[8] * T_sat_d**2*T_sat_s + M[9] * T_sat_d**3

        # Convert mass flow rate to kg/s from lbm/h
        m_dot_map                           *= 0.000125998

        # Add more mass flow rate to scale
        m_dot_map                           *= self.V_dot_ratio
        power_map                           *= self.V_dot_ratio

        P1                                  = self.p_in_r
        P2                                  = self.p_out_r
        T1_actual                           = self.T_sat_s_K + self.DT_sh_K     # T_sat_d_K
        T1_map                              = self.T_sat_s_K + 20 * 5 / 9

        AS.update(CP.PT_INPUTS, P1, T1_map)
        v_map                               = 1 / AS.rhomass()                      # [m^3/kg]
        s1_map                              = AS.smass()                            # [J/kg-K]
        h1_map                              = AS.hmass()                            # [J/kg]

        AS.update(CP.PT_INPUTS, P1, T1_actual)
        s1_actual                           = AS.smass()                            # [J/kg-K]
        h1_actual                           = AS.hmass()                            # [J/kg]
        v_actual                            = 1 / AS.rhomass()                      # [m^3/kg]
        F                                   = 0.75
        m_dot                               = (1 + F * (v_map / v_actual - 1)) * m_dot_map

        AS.update(CP.PSmass_INPUTS, P2, s1_map)
        h2s_map                             = AS.hmass()                            # [J/kg]

        AS.update(CP.PSmass_INPUTS, P2, s1_actual)
        h2s_actual                          = AS.hmass()                            # [J/kg]

        # Shaft power based on 20F superheat calculation from fit overall isentropic efficiency
        power                               = power_map * (m_dot / m_dot_map) * (h2s_actual - h1_actual) / (h2s_map - h1_map)

        h2                                  = power * (1 - self.fp) / m_dot + h1_actual
        self.eta_oi                         = m_dot * (h2s_actual - h1_actual) / power

        AS.update(CP.HmassP_INPUTS, h2, P2)
        self.T_out_r                        = AS.T()                                # [K]
        self.s_out_r                        = AS.smass()                            # [J/kg-K]

        self.s_in_r                         = s1_actual
        self.h_out_r                        = h2
        self.h_in_r                         = h1_actual
        self.m_dot_r                        = m_dot
        self.W                              = power
        self.CycleEnergyIn                  = power * (1 - self.fp)
        self.V_dot_pumped                   = m_dot * v_actual
        self.Q_amb                          = -self.fp * power

        # Estimate refrigerant dissolved in the oil sump
        T_ave                               = (T1_actual + self.T_out_r) / 2
        if self.shell_pressure == 'high-pressure':
            p_shell                         = P2        # outlet pressure
        elif self.shell_pressure == 'low-pressure':
            p_shell                         = P1        # inlet pressure

        self.x_Ref, error                   = Solubility_Ref_in_Liq(self.Ref, self.Oil, T_ave, p_shell/1000)

        AS.update(CP.PT_INPUTS, p_shell, T_ave)
        rho_shell                           = AS.rhomass()                          # [kg/m^3]

        rhomass_oil                         = rho_oil(self.Oil, T_ave - 273.15)
        self.m_oil                          = self.V_oil_sump * rhomass_oil

        # Amount of refrigerant dissolved in the oil sump
        self.Charge                         = self.m_oil * self.x_Ref / (1 - self.x_Ref)


# --------------------------------------------------------------------------------------
if __name__ == '__main__':
    # Abstract State
    Ref                                     = 'R134a'
    Backend                                 = 'HEOS'                                # choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    AS = CP.AbstractState(Backend, Ref)
    for i in range(1):
        kwds={
              'M':                          [217.3163128, 5.094492028, -0.593170311, 4.38E-02, -2.14E-02, 1.04E-02, 7.90E-05, -5.73E-05, 1.79E-04, -8.08E-05],
              'P':                          [-561.3615705, -15.62601841, 46.92506685, -0.217949552, 0.435062616, -0.442400826, 2.25E-04, 2.37E-03, -3.32E-03, 2.50E-03],
              'AS':                         AS,                                     # Abstract state
              'Ref':                        Ref,
              'T_in_r':                     280,
              'p_in_r':                     PropsSI('P', 'T', 279, 'Q', 1, Ref),
              'p_out_r':                    PropsSI('P', 'T', 315, 'Q', 1, Ref),
              'fp':                         0.15,                                   # Fraction of electrical power lost as heat to ambient
              'V_dot_ratio':                1.0,                                    # Displacement Scale factor
              'shell_pressure':             'low-pressure',
              'Oil':                        'POE32',
              'V_oil_sump':                 0.0,
              }
        Comp                                = CompressorClass(**kwds)
        Comp.Calculate()
        print('Power:',                     Comp.W,                 'W')
        print('Flow rate:',                 Comp.V_dot_pumped,      'm^3/s')
        print('Heat loss rate:',            Comp.Q_amb,             'W')
        print('Refrigerant dissolved in oil sump:', Comp.Charge,    'kg')
