from __future__                 import division, print_function, absolute_import
import CoolProp                 as CP


# ----------------------------------------------------------------------
class WavyChannelCoolerClass():
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
            ('Inlet pressure',          'kPa',          self.p_in_g),
            ('Overall Efficiency',      '-',            self.eta),
            ('Q Wavy Channel',          'W',            self.Q_wavychannel),
            ('Inlet coolant temp',      'K',            self.T_in_g),
            ('Outlet coolant temp',     'K',            self.T_out_g),
            ('Pressure drop',           'Pa',           self.DP_g),
            ('Mass flow rate',          'kg/s',         self.m_dot_g)
         ]

    def Calculate(self):
        # AbstractState
        AS_g                            = self.AS_g
        if hasattr(self, 'MassFrac_g'):
            AS_g.set_mass_fractions([self.MassFrac_g])
        elif hasattr(self, 'VoluFrac_SLF'):
            AS_g.set_volu_fractions([self.VoluFrac_g])

        AS_g.update(CP.PT_INPUTS, self.p_in_g, self.T_in_g)
        cp_g                            = AS_g.cpmass()            # [J/kg K]

        self.T_out_g = self.T_in_g + self.Q_wavychannel / cp_g

        print('T_in_wc = ', self.T_in_g, 'T_out_wc = ', self.T_out_g)




