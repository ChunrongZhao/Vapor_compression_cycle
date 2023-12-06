from __future__                                         import division, print_function, absolute_import
from CoolProp.CoolProp                                  import PropsSI
from ACHP_codes.Correlations.HTC_Correlations           import f_h_1phase_Tube
from ACHP_codes.Correlations.PhaseState_Correlations    import TrhoPhase_ph
from math                                               import log, pi, exp
from ACHP_codes.ACHP_Tools.convert_units                import *
import CoolProp                                         as CP


# ------------------------------------------------------------------------------------------
class SightGlassFilterDrierMicroMotionClass():
    def __init__(self, **kwargs):
        # Load the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)

    def Update(self,**kwargs):
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
            ('height of Sight Glass',       'm',        self.h),
            ('Sight glass diameter',        'm',        self.D),
            ('Charge',                      'kg',       self.Charge),
            ('Volume filter drier',         'm^3',      self.V),
            ('Micromotion tube diameter',   'm',        self.D_Micro),
            ('Micormotion tube length',     'm',        self.L_Micro),
            ('Micormotion number of tubes', '-',        self.n_Micro),

         ]

    def Calculate(self):
        # AbstractState
        AS                              = self.AS
        if hasattr(self,'MassFrac'):
            AS.set_mass_fractions([self.MassFrac])
        elif hasattr(self, 'VoluFrac'):
            AS.set_volu_fractions([self.VoluFrac])

        if not 'IncompressibleBackend' in AS.backend_name():
            # Figure out the inlet state
            AS.update(CP.PQ_INPUTS, self.p_in, 0.0)
            self.T_bubble               = AS.T()                    # [K]
            AS.update(CP.PQ_INPUTS, self.p_in, 1.0)
            self.T_dew                  = AS.T()                    # [K]
        else:
            # It is a brine
            self.T_bubble               = None
            self.T_dew                  = None

        self.T_in, self.rho_in, self.Phasein = TrhoPhase_ph(self.AS, self.p_in, self.h_in, self.T_bubble, self.T_dew)

        # TODO: Need to be UPDATED with Petterson et al. (2000) correlation for transcritical CO2
        if self.Phasein == 'Supercritical':
            print("Cauation::Transcritical phase at the inlet of SightGlass during iteration")
            self.f_fluid, self.h_fluid, self.Re_fluid = f_h_1phase_Tube(self.m_dot, self.ID, self.T_in, self.p_in, self.AS)

            AS.update(CP.PT_INPUTS, self.p_in, self.T_in)
            # Specific heat capacity [J/kg-K]
            cp                          = AS.cpmass()
            # Density [kg/m^3]
            rho                         = AS.rhomass()
        elif self.Phasein == 'TwoPhase':
            print("Caution::two phase at the inlet of SightGlass during iteration")
            self.f_fluid, self.h_fluid, self.Re_fluid = f_h_1phase_Tube(self.m_dot, self.ID, self.T_in-1, self.p_in, self.AS)
            AS.update(CP.PT_INPUTS, self.p_in, self.T_in-1)
            # Specific heat capacity [J/kg-K]
            cp                          = AS.cpmass()
            # Density [kg/m^3]
            rho                         = AS.rhomass()
        else:                           # Single phase
            self.f_fluid, self.h_fluid, self.Re_fluid = f_h_1phase_Tube(self.m_dot, self.ID, self.T_in, self.p_in, self.AS)
            AS.update(CP.PT_INPUTS, self.p_in, self.T_in)
            # Specific heat capacity [J/kg-K]
            cp                          = AS.cpmass()
            # Density [kg/m^3]
            rho                         = AS.rhomass()

        # Pressure drop calculations for single phase refrigerant
        v                               = 1. / rho

        # Charge in Sight Glass [kg]
        self.SightGlassCharge           = pi * self.D**2 / 4.0 * self.h * rho
        # Charge in FilterDrier [kg]
        self.FilterDrierCharge          = self.V * rho
        # Charge in MicroMotion [kg]
        self.MicroMotionCharge          = self.n_Micro * pi * self.D_Micro**2 / 4.0 * self.L_Micro * rho
        # Total Charge [kg]
        self.Charge                     = self.n_sight * self.SightGlassCharge + self.FilterDrierCharge + self.MicroMotionCharge


class SightGlassFilterDrierMicroMotionOptionClass():
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
            ('height of Sight Glass',           'm',        self.h),
            ('Sight glass diameter',            'm',        self.D),
            ('Charge',                          'kg',       self.Charge),
            ('Volume filter drier',             'm^3',      self.V),
            ('Micromotion tube diameter',       'm',        self.D_Micro),
            ('Micormotion tube length',         'm',        self.L_Micro),
            ('Micormotion number of tubes',     '-',        self.n_Micro)]

    def Calculate(self):
        # AbstractState
        AS                                      = self.AS
        if hasattr(self, 'MassFrac'):
            AS.set_mass_fractions([self.MassFrac])
        elif hasattr(self, 'VoluFrac'):
            AS.set_volu_fractions([self.VoluFrac])

        if not 'IncompressibleBackend' in AS.backend_name():
            # Figure out the inlet state
            AS.update(CP.PQ_INPUTS, self.p_in, 0.0)
            self.T_bubble                       = AS.T()            # [K]
            AS.update(CP.PQ_INPUTS, self.p_in, 1.0)
            self.T_dew                          = AS.T()            # [K]
        else:
            # It is a brine
            self.T_bubble                       = None
            self.T_dew                          = None

        if hasattr(self, 'SystemVolumeOption'):                     # check if LineSetOption is given
            if self.SystemVolumeOption == 'On':
                self.T_in, self.rho_in, self.Phasein = TrhoPhase_ph(self.AS, self.p_in, self.h_in, self.T_bubble, self.T_dew)

                # TODO: Need to be UPDATED with Petterson et al. (2000) correlation for transcritical CO2
                if self.Phasein == 'Supercritical':
                    print("Caution::Transcritical phase at the inlet of SightGlass during iteration")
                    self.f_fluid, self.h_fluid, self.Re_fluid   = f_h_1phase_Tube(self.m_dot, self.ID, self.T_in, self.p_in, self.AS)
                    AS.update(CP.PT_INPUTS, self.p_in, self.T_in)
                    # Specific heat capacity [J/kg-K]
                    cp                          = AS.cpmass()
                    # Density [kg/m^3]
                    rho                         = AS.rhomass()
                elif self.Phasein == 'TwoPhase':
                    print("Caution::two phase at the inlet of SightGlass during iteration")
                    self.f_fluid, self.h_fluid, self.Re_fluid   = f_h_1phase_Tube(self.m_dot, self.ID, self.T_in-1, self.p_in, self.AS)
                    AS.update(CP.PT_INPUTS, self.p_in, self.T_in-1)
                    # Specific heat capacity [J/kg-K]
                    cp                          = AS.cpmass()
                    # Density [kg/m^3]
                    rho                         = AS.rhomass()
                else:                           # Single phase
                    self.f_fluid, self.h_fluid, self.Re_fluid   = f_h_1phase_Tube(self.m_dot, self.ID, self.T_in, self.p_in, self.AS)
                    AS.update(CP.PT_INPUTS, self.p_in, self.T_in)
                    # Specific heat capacity [J/kg-K]
                    cp                          = AS.cpmass()
                    # Density [kg/m^3]
                    rho                         = AS.rhomass()

                # Pressure drop calculations for single phase refrigerant
                v                               = 1. / rho

                # Charge in Sight Glass [kg]
                self.SightGlassCharge           = pi * self.D**2 / 4.0 * self.h * rho
                # Charge in FilterDrier [kg]
                self.FilterDrierCharge          = self.V * rho
                # Charge in MicroMotion [kg]
                self.MicroMotionCharge          = self.n_Micro * pi * self.D_Micro**2 / 4.0 * self.L_Micro * rho
                # Total Charge [kg]
                self.Charge                     = self.n_sight * self.SightGlassCharge + self.FilterDrierCharge + self.MicroMotionCharge

            else:                               # SystemVolumeOption == 'Off'
                print('SystemVolumeOption is off for ' + str(self.Name))
                # Charge in Sight Glass [kg]
                self.SightGlassCharge           = 0
                # Charge in FilterDrier [kg]
                self.FilterDrierCharge          = 0
                # Charge in MicroMotion [kg]
                self.MicroMotionCharge          = 0
                # Total Charge [kg]
                self.Charge                     = 0
        else:
            print(str(self.Name) + ' is neglected')
            # Charge in Sight Glass [kg]
            self.SightGlassCharge               = 0
            # Charge in FilterDrier [kg]
            self.FilterDrierCharge              = 0
            # Charge in MicroMotion [kg]
            self.MicroMotionCharge              = 0
            # Total Charge [kg]
            self.Charge                         = 0


# -------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # Abstract State
    Ref                             = 'R410A'
    Backend                         = 'HEOS'    # choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    AS                              = CP.AbstractState(Backend, Ref)

    kwargs  = {
            'AS':                   AS,
            'p_in':                 500000,
            'h_in':                 PropsSI('H', 'P', 500000, 'T', PropsSI('T', 'P', 500000, 'Q', 0, Ref) - 10, Ref),
            'm_dot':                0.03,
            'ID':                   in2m(3.0/8.0)-mm2m(2),      # pipe diameter
            'h':                    in2m(1.370),                # height of sight glass in m
            'D':                    in2m(1.110),                # diameter of sight glass in m
            'n_sight':              2,                          # number of sight glasses
            'V':                    cubin2cubm(16),             # volume of filter drier (website = 13.74in^3 calculated) (manual = 16in^3)
            'D_Micro':              in2m(0.21),                 # micromotion tube diameter
            'L_Micro':              in2m(14.6),                 # micormotion tube length
            'n_Micro':              2,                          # micormotion number of tubes
            }

    SightGlassFilterDrierMicroMotionClass   = SightGlassFilterDrierMicroMotionClass(**kwargs)
    SightGlassFilterDrierMicroMotionClass.Calculate()
    print(SightGlassFilterDrierMicroMotionClass.OutputList())
