# WavyChan.py
#
# Created: Oct. 2022, C.R. Zhao

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import numpy as np
import pickle


# ----------------------------------------------------------------------
#  heat generation rate of a single battery cell
# ----------------------------------------------------------------------
class WavyChan_NMC():

    def __init__(self):
        self.tag = "LIB_NMC_18650"
        self.diameter = 0.018  # [m]
        self.height = 0.065  # [m]
        self.mass = 0.048  # [kg]
        self.surface_area = (np.pi * self.height * self.diameter) + (0.5 * np.pi * self.diameter ** 2)  # [m^2]
        self.volume = np.pi * (0.5 * self.diameter) ** 2 * self.height
        self.density = self.mass / self.volume  # [kg/m^3]
        self.electrode_area = 0.0342  # [m^2]

        self.max_voltage = 4.2  # [V]
        self.nominal_capacity = 3.55  # [Amp-Hrs]
        self.nominal_voltage = 3.6  # [V]
        self.charging_voltage = self.nominal_voltage  # [V]

        self.watt_hour_rating = self.nominal_capacity * self.nominal_voltage  # [Watt-hours]
        self.specific_energy = self.watt_hour_rating * 3600 / self.mass  # [J/kg]
        self.specific_power = self.specific_energy / self.nominal_capacity  # [W/kg]  need to be required
        self.resistance = 0.025  # [Ohms]

        self.specific_heat_capacity = 1108  # [J/kgK]
        self.radial_thermal_conductivity = 0.4  # [J/kgK]
        self.axial_thermal_conductivity = 32.2  # [J/kgK] # estimated

    def compute_flow_resistance_and_parameters(self, m_dot_coolant=None, dry_mass=False, Tiltwing=True, a=1e-3, b=5e-4, c=6.3e-2, d=2e-3,
                                               single_side_contact=True, N_sheet=10, N_series=None, N_parallel=None):
        """ source: Zhao C, Sousa A C M, Jiang F. Minimization of thermal non-uniformity in lithium-ion battery pack
        cooled by channeled liquid flow[J]. International journal of heat and mass transfer, 2019, 129: 660-670."""
        Theta_inlet = 44.5  # contact angle at inlet
        Theta_outlet = 50.5

        # thermopysical properties
        rho_coolant = 1075
        k_coolant = 0.387
        cp_coolant = 3300
        mu_coolant = 0.0019
        Pr_coolant = cp_coolant * mu_coolant / k_coolant

        rho_chan = 2719
        k_chan = 202.4
        cp_chan = 871

        # -------------------------------------------------------------
        N_s = N_series / N_sheet
        N_p = N_parallel      # todo check before run
        # -------------------------------------------------------------

        N_tot = N_s * N_p
        # assumptions
        L_u_turn = 2 * self.diameter  # U-turn's length for single side contact scenario
        L_res = self.diameter  # the rest length to connect the manifold of each inlet/outlet

        # contact surface of a single battery cell
        if single_side_contact:
            A_surf = (Theta_inlet + Theta_outlet) / 2 / 360 * np.pi * self.diameter * self.height
            L_surf = (Theta_inlet + Theta_outlet) / 2 / 360 * np.pi * (self.diameter + b + d / 2)
        else:
            A_surf = 2 * (Theta_inlet + Theta_outlet) / 2 / 360 * np.pi * self.diameter * self.height
            L_surf = 2 * (Theta_inlet + Theta_outlet) / 2 / 360 * np.pi * (self.diameter + b + d / 2)

        # contact surface of a module
        if single_side_contact:
            A_chan_module = N_tot * A_surf  # 10 'modules', means 1 'sheet' here
            L_chan_module = N_tot * L_surf + L_u_turn + L_res
        else:
            A_chan_module = 2 * N_tot * A_surf
            L_chan_module = 2 * N_tot * L_surf + 2 * L_res

        # hydraulic diameter
        d_H = 2 * (c * d) / (c + d)
        # channel aspect ratio
        gamma_chan = d / c
        # channel cross-section area
        A_c_chan = c * d
        "********************************************************************************"
        # flow velocity, 65 cells a line, each channel contact 2 lines (each sheet has 15 )
        m_dot_unit  = m_dot_coolant / (N_s * N_sheet)    # 10 sheets
        u_coolant   = m_dot_unit / (rho_coolant * A_c_chan)
        "********************************************************************************"
        # Reynolds number
        Re_coolant = u_coolant * d_H / (mu_coolant / rho_coolant)

        if Re_coolant <= 2300:
            # Fanning friction factor; c_geom constant: f * Re
            f_coolant = 24 * (1 - 1.3553 * gamma_chan + 1.9467 * np.power(gamma_chan, 2) - 1.7012 * np.power(gamma_chan, 3)
                           + 0.9564 * np.power(gamma_chan, 4) - 0.2537 * np.power(gamma_chan, 5)) / Re_coolant
            # Nusselt number
            Nu = 8.235 * (1 - 2.0421 * gamma_chan + 3.0853 * np.power(gamma_chan, 2) - 2.4765 * np.power(gamma_chan, 3)
                          + 1.0578 * np.power(gamma_chan, 4) - 0.1861 * np.power(gamma_chan, 5))
            # Colburn factor
            j_coolant = Nu / Re_coolant * np.power(Pr_coolant, -1 / 3)
        else:
            # Fanning friction factor
            f_coolant = 1 / (4 * (1.8 * np.log10(Re_coolant / 7.7)) ** 2)
            # use Gnielinski equation to calculate Nu for 0.5 < Pr < 2000, 3000 < Re < 5e6
            Nu = (f_coolant / 2) * (Re_coolant - 1000) * Pr_coolant / (
                    1 + 12.7 * np.power((f_coolant / 2), 0.5) * (np.power(Pr_coolant, 2 / 3) - 1))
            # Colburn factor
            j_coolant = Nu / Re_coolant * np.power(Pr_coolant, -1 / 3)

        # heat transfer coefficient of channeled coolant fluid
        h_coolant = k_coolant * Nu / d_H
        # total heat transfer coefficient
        h_tot = 1 / (1 / h_coolant + b / k_chan)

        # the pressure drop due to friction
        """hydraulic radius: r_h = A_o / P = d_H / 4"""
        delta_P_chan = 2 * f_coolant * (rho_coolant * u_coolant * u_coolant) * (L_chan_module / d_H)

        # line density of a cross-section: dry; with liquid coolant: + rho_coolant * A_c_chan
        if dry_mass:
            rho_line = rho_chan * (2 * a * (2 * b + d) + 2 * b * c)
        else:
            rho_line = rho_chan * (2 * a * (2 * b + d) + 2 * b * c) + rho_coolant * A_c_chan

        # NTU value
        NTU = h_tot * A_chan_module / ((m_dot_unit * N_s) * cp_coolant)
        # liquid wavy mini-channel efficiency
        eff_WavyChan = 1 - np.exp(-NTU)
        # eff_WavyChan_check = (T_o - T_i) / (T_current - T_i)

        return A_chan_module, L_chan_module, h_tot, rho_line, eff_WavyChan, delta_P_chan

    def compute_net_generated_battery_heat_module(self, Q_heat_gen=None, m_dot_coolant=None, T_cell=None, T_i=None,
                                                Tiltwing=True, dry_mass=False, heat_transfer_efficiency=1, dt=None,
                                                  N_sheet=10, N_series=None, N_parallel=None):
        """Computes the net heat generated in a battery module during cycling.
        Assumptions:
        1) Battery pack cell heat transfer can be modelled as a cooling columns in a cross-flow
        2) Isothermal battery cell - the temperature at the center of the cell is the same at
        the surface of the cell

        Inputs:
            battery.
                  h                         (heat transfer coefficient)  [W/(m^2*K)]
                  As_cell                   (battery cell surface area)  [meters^2]
                  H_cell                    (battery cell height)        [meters]
                  T_ambient                 (ambient temperature)        [Kelvin]
                  T_current                 (pack temperature)           [Kelvin]
                  T_cell                    (battery cell temperature)   [Kelvin]
                  heat_transfer_efficiency                               [unitless]

          Outputs:
            battery.
                 net_power                                               [Watts]
        """
        # thermopysical properties
        rho_coolant = 1075
        k_coolant = 0.387
        cp_coolant = 3300
        mu_coolant = 0.0019
        Pr_coolant = cp_coolant * mu_coolant / k_coolant

        rho_chan = 2719
        k_chan = 202.4
        cp_chan = 871

        # -------------------------------------------------------------
        N_s = N_series / N_sheet
        N_p = N_parallel      # todo check before run

        # -------------------------------------------------------------
        N_tot = N_s * N_p

        "*************************************************************************************"
        A_chan_module, L_chan_module, h_tot, rho_line, eff_WavyChan, delta_P_chan = \
            self.compute_flow_resistance_and_parameters(m_dot_coolant=m_dot_coolant, dry_mass=dry_mass, Tiltwing=Tiltwing,N_series=N_series, N_parallel=N_parallel)
        "*************************************************************************************"

        T_current = T_cell
        Tw_Ti = T_current - T_i
        "*************************************************************************************"
        # a unit has 65 * 2 cells, thus total mass flow rate is (150*130/(65*2)) times higher than a distributed mass flow rate within a unit
        if Tiltwing:    # check if it is tiltwing aircraft
            Number_of_flow_channels = 150*130/(65*2)
        else:           # otherwise it is a stopped rotor aircraft
            Number_of_flow_channels = 140*100/(50*2)

        m_dot_unit  = m_dot_coolant / Number_of_flow_channels
        "*************************************************************************************"
        # a module has (130*15/(65*2)) units
        if Tiltwing:    # check if it is tiltwing aircraft
            Channels_in_a_module = 130*15/(65*2)
        else:           # otherwise it is a stopped rotor aircraft
            Channels_in_a_module = 100*14/(50*2)

        NTU = h_tot * A_chan_module / (m_dot_unit * Channels_in_a_module * cp_coolant)
        Tw_To = Tw_Ti * np.exp(-NTU)

        # outlet temperature of the liquid coolant
        T_o = T_current - Tw_To

        # log mean temperature difference
        dT_lm = (Tw_Ti - Tw_To) / np.log(Tw_Ti / Tw_To)

        # convective heat transfer for a module
        Q_conv = h_tot * A_chan_module * dT_lm
        # Q_conv[Tw_Ti == 0.] = 0.  # avoid being divided by zero
        Q_conv_right = m_dot_unit * Channels_in_a_module * cp_coolant * (T_o - T_i)

        delta_Q_conv = np.abs(Q_conv - Q_conv_right)

        eff_WavyChan_check = (T_o - T_i) / (T_current - T_i)

        # temperature rise
        P_net = Q_heat_gen - Q_conv
        dT_dt = P_net / (N_tot * self.mass * self.specific_heat_capacity)

        # update T_current (cell temperature)
        T_current = T_current + dT_dt * dt * 60

        return T_current, Q_conv, T_o, delta_P_chan

    def compute_power_and_mass(self, m_dot_coolant=None, rho_coolant=1075, eff_PUMP=0.7, dry_mass=False, Tiltwing=True, N_sheet=10, N_series=None, N_parallel=None):

        A_chan_module, L_chan_module, h_tot, rho_line, eff_WavyChan, delta_P_chan = \
            self.compute_flow_resistance_and_parameters(m_dot_coolant=m_dot_coolant, dry_mass=dry_mass, Tiltwing=Tiltwing, N_series=N_series, N_parallel=N_parallel)

        # pump power consumption for the battery pack side channel
        P_chan = N_sheet * m_dot_coolant * delta_P_chan / (1e3 * rho_coolant * eff_PUMP)
        # mass of channelled liquid of the pack
        m_chan = N_sheet * rho_line * L_chan_module

        return P_chan, m_chan, eff_WavyChan


class WavyChan_data:
    def __init__(self):
        self.T_current      = []
        self.Q_conv         = []
        self.T_o            = []
        # new values
        self.eff_WavyChan   = []
        self.P_chan         = []
        self.m_chan         = []

    def save_data(self, filename):
        with open(filename, 'wb') as f:
            pickle.dump(self, f)
