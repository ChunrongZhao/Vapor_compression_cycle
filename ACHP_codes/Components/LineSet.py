from __future__                                         import division, print_function, absolute_import
from CoolProp.CoolProp                                  import PropsSI
from ACHP_codes.Correlations.HTC_Correlations           import f_h_1phase_Tube
from ACHP_codes.Correlations.PhaseState_Correlations    import TrhoPhase_ph
from math                                               import log, pi, exp
import CoolProp                                         as CP


# ------------------------------------------------------------------------------------------
class LineSetClass():
    def __init__(self, **kwargs):
        # Load the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)

    def Update(self, **kwargs):
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
            ('Length of tube',              'm',            self.L),
            ('Supply line OD',              'm',            self.OD),
            ('Supply line ID',              'm',            self.ID),
            ('Tube Conductivity',           'W/m-K',        self.k_tube),
            ('Insulation thickness',        'm',            self.t_insul),
            ('Insulation conductivity',     'W/m-K',        self.k_insul),
            ('Air overall HTC',             'W/m^2-K',      self.h_air),
            ('Air Temperature',             'K',            self.T_air),
            ('Q Total',                     'W',            self.Q),
            ('Pressure drop ',              'Pa',           self.DP),
            ('Reynolds # Fluid',            '-',            self.Re_fluid),
            ('Mean HTC Fluid',              'W/m^2-K',      self.h_fluid),
            ('Charge',                      'kg',           self.Charge),
            ('Inlet Temperature',           'K',            self.T_in),
            ('Outlet Temperature',          'K',            self.T_out)
         ]

    def Calculate(self):
        # AbstractState
        if hasattr(self, 'Backend'):                        # check if backend is given
            AS                      = CP.AbstractState(self.Backend, self.Ref)
            if hasattr(self, 'MassFrac'):
                AS.set_mass_fractions([self.MassFrac])
        else:                                               # otherwise, use the defualt backend
            AS = CP.AbstractState('HEOS', self.Ref)

        self.AS                     = AS

        if not 'IncompressibleBackend' in AS.backend_name():
            # Figure out the inlet state
            AS.update(CP.PQ_INPUTS, self.p_in, 0.0)
            self.T_bubble           = AS.T()                # [K]
            AS.update(CP.PQ_INPUTS, self.p_in, 1.0)
            self.T_dew              = AS.T()                # [K]
        else:
            # It is a brine
            self.T_bubble           = None
            self.T_dew              = None

        self.T_in,self.rho_in,self.Phase_in = TrhoPhase_ph(self.AS, self.p_in, self.h_in, self.T_bubble, self.T_dew)
        # Solver shows TwoPhase in the first iteration,
        # the following if statement just to avoid ValueError with CoolProp for pseudo-pure refrigerants
        if self.Phase_in == 'TwoPhase':
            self.f_fluid, self.h_fluid, self.Re_fluid   = f_h_1phase_Tube(self.m_dot, self.ID, self.T_in-1, self.p_in, self.AS)
            AS.update(CP.PT_INPUTS, self.p_in, self.T_in-1)
            # Specific heat capacity [J/kg-K]
            cp                      = AS.cpmass()
            # Density [kg/m^3]
            rho                     = AS.rhomass()
        else: #Single phase
            self.f_fluid, self.h_fluid, self.Re_fluid   = f_h_1phase_Tube(self.m_dot, self.ID, self.T_in, self.p_in, self.AS)
            AS.update(CP.PT_INPUTS, self.p_in, self.T_in)
            # Specific heat capacity [J/kg-K]
            cp                      = AS.cpmass()
            # Density [kg/m^3]
            rho                     = AS.rhomass()

        # Thermal resistance of tube
        R_tube                      = log(self.OD / self.ID) / (2 * pi * self.L * self.k_tube)
        # Thermal resistance of insulation
        R_insul                     = log((self.OD + 2.0 * self.t_insul) / self.OD) / (2 * pi * self.L * self.k_insul)
        # Convective UA for inside the tube
        UA_i                        = pi * self.ID * self.L * self.h_fluid
        # Convective UA for the air-side
        UA_o                        = pi * (self.OD + 2 * self.t_insul) * self.L * self.h_air

        # Avoid the possibility of division by zero if h_air is zero
        if UA_o < 1e-12:
            UA_o                    = 1e-12

        # Overall UA value
        UA                          = 1 / (1 / UA_i + R_tube + R_insul + 1 / UA_o)

        # Outlet fluid temperature [K]
        self.T_out                  = self.T_air - exp(-UA / (self.m_dot * cp)) * (self.T_air - self.T_in)
        # Overall heat transfer rate [W]
        self.Q                      = self.m_dot * cp * (self.T_out - self.T_in)
        self.h_out                  = self.h_in + self.Q / self.m_dot

        # Pressure drop calculations for single phase refrigerant
        v                           = 1. / rho
        G                           = self.m_dot / (pi * self.ID**2 / 4.0)
        # Pressure gradient using Darcy friction factor
        dpdz                        = -self.f_fluid * v * G**2 / (2. * self.ID)     # Pressure gradient
        self.DP                     = dpdz * self.L

        # Charge in Line set [kg]
        self.Charge                 = pi * self.ID**2 / 4.0 * self.L * rho


class LineSetOptionClass():
    def __init__(self, **kwargs):
        # Load the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)

    def Update(self, **kwargs):
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
            ('Length of tube',          'm',            self.L),
            ('Supply line OD',          'm',            self.OD),
            ('Supply line ID',          'm',            self.ID),
            ('Tube Conductivity',       'W/m-K',        self.k_tube),
            ('Insulation thickness',    'm',            self.t_insul),
            ('Insulation conductivity', 'W/m-K',        self.k_insul),
            ('Air overall HTC',         'W/m^2-K',      self.h_air),
            ('Air Temperature',         'K',            self.T_air),
            ('Q Total',                 'W',            self.Q),
            ('Pressure drop ',          'Pa',           self.DP),
            ('Reynolds # Fluid',        '-',            self.Re_fluid),
            ('Mean HTC Fluid',          'W/m^2-K',      self.h_fluid),
            ('Charge',                  'kg',           self.Charge),
            ('Inlet Temperature',       'K',            self.T_in),
            ('Outlet Temperature',      'K',            self.T_out)]

    def Calculate(self):

        # AbstractState
        AS                                  = self.AS
        if hasattr(self, 'MassFrac'):
            AS.set_mass_fractions([self.MassFrac])
        elif hasattr(self, 'VoluFrac'):
            AS.set_volu_fractions([self.VoluFrac])

        if not 'IncompressibleBackend' in AS.backend_name():
            # Figure out the inlet state
            AS.update(CP.PQ_INPUTS, self.p_in, 0.0)
            self.T_bubble                   = AS.T()                    # [K]
            AS.update(CP.PQ_INPUTS, self.p_in, 1.0)
            self.T_dew                      = AS.T()                    # [K]

            # Check for supercritical state
            if self.p_in > AS.p_critical():
                # Supercritical
                self.T_bubble               = None
                self.T_dew                  = None

        else:
            # It is a brine or incompressible
            self.T_bubble                   = None
            self.T_dew                      = None

        if hasattr(self, 'LineSetOption'):                              # check if LineSetOption is given

            if self.LineSetOption == 'On':

                self.T_in, self.rho_in, self.Phase_in       = TrhoPhase_ph(self.AS, self.p_in, self.h_in, self.T_bubble, self.T_dew)
                # Solver shows TwoPhase in the first iteration,
                # the following if statement just to avoid ValueError with CoolProp for pseudo-pure refrigerants
                if self.Phase_in == 'Supercritical':
                    # TO DO: Need to be UPDATED with Petterson et al. (2000) correlation for transcritical CO2
                    self.f_fluid, self.h_fluid, self.Re_fluid   = f_h_1phase_Tube(self.m_dot, self.ID, self.T_in, self.p_in, self.AS)
                    AS.update(CP.PT_INPUTS, self.p_in, self.T_in)
                    # Specific heat capacity [J/kg-K]
                    cp                      = AS.cpmass()
                    # Density [kg/m^3]
                    rho                     = AS.rhomass()
                    # Specific entropy [J/kg-K]
                    s_in                    = AS.smass()
                elif self.Phase_in == 'TwoPhase':
                    self.f_fluid, self.h_fluid, self.Re_fluid   = f_h_1phase_Tube(self.m_dot, self.ID, self.T_in-1, self.p_in, self.AS)
                    AS.update(CP.PT_INPUTS, self.p_in, self.T_in-1)
                    # Specific heat capacity [J/kg-K]
                    cp                      = AS.cpmass()
                    # Density [kg/m^3]
                    rho                     = AS.rhomass()
                    # Specific entropy [J/kg-K]
                    s_in                    = AS.smass()
                else:                                                       # Single phase
                    self.f_fluid, self.h_fluid, self.Re_fluid   = f_h_1phase_Tube(self.m_dot, self.ID, self.T_in, self.p_in, self.AS)
                    AS.update(CP.PT_INPUTS, self.p_in, self.T_in)
                    # Specific heat capacity [J/kg-K]
                    cp                      = AS.cpmass()
                    # Density [kg/m^3]
                    rho                     = AS.rhomass()
                    # Specific entropy [J/kg-K]
                    s_in                    = AS.smass()

                # Inlet specific entropy
                self.s_in                   = s_in
                # Thermal resistance of tube
                R_tube                      = log(self.OD / self.ID) / (2 * pi * self.L * self.k_tube)
                # Thermal resistance of insulation
                R_insul                     = log((self.OD + 2.0 * self.t_insul) / self.OD) / (2 * pi * self.L * self.k_insul)
                # Convective UA for inside the tube
                UA_i                        = pi * self.ID * self.L * self.h_fluid
                # Convective UA for the air-side
                UA_o                        = pi * (self.OD + 2 * self.t_insul) * self.L * self.h_air

                # Avoid the possibility of division by zero if h_air is zero
                if UA_o < 1e-12:
                    UA_o                    = 1e-12

                # Overall UA value
                UA                          = 1 / (1 / UA_i + R_tube + R_insul + 1 / UA_o)

                # Outlet fluid temperature [K]
                self.T_out                  = self.T_air - exp(-UA / (self.m_dot * cp)) * (self.T_air - self.T_in)
                # Overall heat transfer rate [W]
                self.Q                      = self.m_dot * cp * (self.T_out - self.T_in)
                self.h_out                  = self.h_in + self.Q / self.m_dot

                # Pressure drop calculations for single phase refrigerant
                v                           = 1. / rho
                G                           = self.m_dot / (pi * self.ID**2 / 4.0)
                # Pressure gradient using Darcy friction factor
                dpdz                        = -self.f_fluid * v * G**2 / (2. * self.ID)     # Pressure gradient
                self.DP                     = dpdz * self.L

                # Outlet specific entropy
                AS.update(CP.HmassP_INPUTS, self.h_out, self.p_in + self.DP)
                self.s_out                  = AS.smass()

                # Charge in Line set [kg]
                self.Charge                 = pi * self.ID**2 / 4.0 * self.L * rho

            else:                                           # LineSetOption == 'Off'
                print('lineSetOption is off for ' + str(self.Name))
                self.T_out                  = self.T_in

                # Inlet specific entropy
                AS.update(CP.HmassP_INPUTS, self.h_in, self.p_in)
                self.s_in                   = AS.smass()

                # Overall heat transfer rate [W]
                self.Q                      = 0.0
                self.h_out                  = self.h_in

                self.Re_fluid               = 0.0
                self.h_fluid                = 0.0

                # Pressure gradient using Darcy friction factor
                self.DP                     = 0.0

                # Outlet specific entropy
                AS.update(CP.HmassP_INPUTS, self.h_out, self.p_in + self.DP)
                self.s_out                  = AS.smass()

                # Charge in Line set [kg]
                self.Charge                 = 0.0

        else:
            print(str(self.Name) + ' is neglected')
            self.T_out                      = self.T_in

            # Inlet specific entropy
            AS.update(CP.HmassP_INPUTS, self.h_in, self.p_in)
            self.s_in                       = AS.smass()

            # Overall heat transfer rate [W]
            self.Q                          = 0.0
            self.h_out                      = self.h_in

            self.Re_fluid                   = 0.0
            self.h_fluid                    = 0.0

            # Pressure gradient using Darcy friction factor
            self.DP                         = 0.0

            # Outlet specific entropy
            AS.update(CP.HmassP_INPUTS, self.h_out, self.p_in + self.DP)
            self.s_out                      = AS.smass()

            # Charge in Line set [kg]
            self.Charge                     = 0.0


# -----------------------------------------------------------------------------------------
if __name__ == '__main__':
    # Abstract State
    Ref                                     = 'R410A'
    Backend                                 = 'HEOS'    # choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    AS                                      = CP.AbstractState(Backend, Ref)

    kwargs  = {
            'L':                            7.6,
            'k_tube':                       0.19,
            't_insul':                      0.02,
            'k_insul':                      0.036,
            'T_air':                        297,
            'AS':                           AS,
            'h_air':                        0.0000000001,
            'T_in':                         PropsSI('T', 'P', 500000, 'Q', 0, Ref) - 10,
            'p_in':                         500000,                 # Pressure of the fluid at the inlet i.e 500kPa
            'h_in':                         PropsSI('H', 'P', 500000, 'T', PropsSI('T', 'P', 500000, 'Q', 0, Ref) - 10, Ref),
                                            # Enthalpy of the fluid at the inlet i.e subcooled 10K below bubble
            'm_dot':                        0.03,                   # fluid mass flow rate
            'LineSetOption':                'On',
            'Name':                         'LineSet',
            }

    LineSetSupply                           = LineSetOptionClass(**kwargs)
    LineSetSupply.OD                        = 0.009525
    LineSetSupply.ID                        = 0.007986
    LineSetSupply.Calculate()

    LineSetReturn                           = LineSetOptionClass(**kwargs)
    LineSetReturn.OD                        = 0.01905
    LineSetReturn.ID                        = 0.017526
    LineSetReturn.Calculate()

    print(LineSetSupply.OutputList())
    print(LineSetReturn.OutputList())
