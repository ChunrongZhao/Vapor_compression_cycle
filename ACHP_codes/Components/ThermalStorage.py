# ThermalStorage.py
#
# Created: May 2024, C.R. Zhao

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import numpy as np
import pickle


# ----------------------------------------------------------------------
#  Main
# ----------------------------------------------------------------------
def PCM_composite():
    # PCM properties; RT21
    # https://www.rubitherm.eu/media/products/datasheets/Techdata_-RT21_EN_09042024.PDF
    PCM_density_s               = 840           # kg/m3
    PCM_density_l               = 770           # kg/m3
    PCM_thermal_conductivity    = 0.2           # W/m.K
    PCM_heat_capacity           = 2000          # J/kg.K
    PCM_latent_heat             = 165e3         # J/kg
    PCM_melting_temperature     = 21 + 273.15   # K

    # graphite periodic structure
    GPS_density                 = 641           # kg/m3
    GPS_thermal_conductivity    = 2000          # W/m.K
    GPS_heat_capacity           = 710           # J/kg.K
    GPS_porosity                = 0.95

    # composite
    LHTES_composite_density                 = (1 - GPS_porosity) * GPS_density + GPS_porosity * (PCM_density_s + PCM_density_l) / 2
    LHTES_composite_thermal_conductivity_p  = (1 - GPS_porosity) * GPS_thermal_conductivity + GPS_porosity * PCM_thermal_conductivity
    LHTES_composite_thermal_conductivity_s  = 1 / ((1 - GPS_porosity) / GPS_thermal_conductivity + GPS_porosity / PCM_thermal_conductivity)
    LHTES_composite_thermal_conductivity    = 0.41 * LHTES_composite_thermal_conductivity_p + 0.59 * LHTES_composite_thermal_conductivity_s
    LHTES_composite_heat_capacity           = (1 - GPS_porosity) * GPS_heat_capacity + GPS_porosity * PCM_heat_capacity
    LHTES_composite_latent_heat             = GPS_porosity * PCM_latent_heat

    return LHTES_composite_density, LHTES_composite_thermal_conductivity, LHTES_composite_heat_capacity, LHTES_composite_latent_heat, PCM_melting_temperature


class LatentHeatThermalEnergyStorage:

    def __init__(self):
        self.tag                = "PCM_GF_composite"
        self._subchannel_width  = 2e-3  # d
        self._subchannel_height = 0.2   # c
        self._subchannel_tw     = 1e-4  # b
        self._subchannel_rho    = 2719  # kg/m3
        self._subchannel_k      = 202.4 # W/m.K
        self._subchannel_number = 40    # todo check
        self._subchannel_length = 0.4   # todo check
        self._LHTES_thickness   = 1e-2  # todo check
        self._coolant_rho       = 1075  # kg/m3
        self._coolant_k         = 0.387 # W/m.K
        self._coolant_cp        = 3300
        self._coolant_mu        = 0.0019
        self._coolant_Pr        = self._coolant_cp * self._coolant_mu / self._coolant_k

    def LHTES_composite_properties(self):
        rho_eff, k_eff, cp_eff, h_sf_eff, T_melt = PCM_composite()

        self.rho_eff            = rho_eff
        self.cp_eff             = cp_eff
        self.k_eff              = k_eff
        self.h_sf_eff           = h_sf_eff
        self.T_melt             = T_melt

    def LHTES_fluid_cal(self, m_dot_coolant=None, T_i_coolant=None, r_melt_front=None, dt=None, eff_PUMP=0.7):
        # todo check: r_melt_front < self._LHTES_thickness
        # --------------Reynold number ---------------------------------
        # hydraulic diameter
        d_H             = 2 * (self._subchannel_height * self._subchannel_width) / (self._subchannel_height + self._subchannel_width)
        # channel aspect ratio
        gamma_chan      = self._subchannel_width / self._subchannel_height
        # channel cross-section area
        A_c_chan        = self._subchannel_height * self._subchannel_width

        m_dot_unit      = m_dot_coolant / self._subchannel_number
        u_coolant       = m_dot_unit / (self._coolant_rho * A_c_chan)

        # Reynolds number
        Re_coolant = u_coolant * d_H / (self._coolant_mu / self._coolant_rho)

        # ------------Nusselt number & friction factor ---------------
        if Re_coolant <= 2300:
            # Fanning friction factor; c_geom constant: f * Re
            f_coolant = 24 * (1 - 1.3553 * gamma_chan + 1.9467 * np.power(gamma_chan, 2) - 1.7012 * np.power(gamma_chan, 3)
                           + 0.9564 * np.power(gamma_chan, 4) - 0.2537 * np.power(gamma_chan, 5)) / Re_coolant
            # Nusselt number
            Nu = 8.235 * (1 - 2.0421 * gamma_chan + 3.0853 * np.power(gamma_chan, 2) - 2.4765 * np.power(gamma_chan, 3)
                          + 1.0578 * np.power(gamma_chan, 4) - 0.1861 * np.power(gamma_chan, 5))
            # Colburn factor
            j_coolant = Nu / Re_coolant * np.power(self._coolant_Pr, -1 / 3)
        else:
            # Fanning friction factor
            f_coolant = 1 / (4 * (1.8 * np.log10(Re_coolant / 7.7)) ** 2)
            # use Gnielinski equation to calculate Nu for 0.5 < Pr < 2000, 3000 < Re < 5e6
            Nu = (f_coolant / 2) * (Re_coolant - 1000) * self._coolant_Pr / (
                    1 + 12.7 * np.power((f_coolant / 2), 0.5) * (np.power(self._coolant_Pr, 2 / 3) - 1))
            # Colburn factor
            j_coolant = Nu / Re_coolant * np.power(self._coolant_Pr, -1 / 3)

        # ------------overall heat transfer coefficient -------------------------
        # heat transfer coefficient of channeled coolant fluid
        h_coolant = self._coolant_k * Nu / d_H
        # total heat transfer coefficient
        self.LHTES_composite_properties()
        h_tot = 1 / (1 / h_coolant + self._subchannel_tw / self._subchannel_k + r_melt_front / self.k_eff)

        # -------------- heat transfer area --------------------------------------
        A_LHTES         = 2* self._subchannel_height * self._subchannel_length * self._subchannel_number
        # NTU value
        NTU             = h_tot * A_LHTES / (m_dot_coolant * self._coolant_cp)
        # efficiency
        eff_LHTES       = 1 - np.exp(-NTU)

        Tm_Ti = self.T_melt - T_i_coolant
        Tm_To = Tm_Ti * np.exp(-NTU)
        # log mean temperature difference
        dT_lm = (Tm_Ti - Tm_To) / np.log(Tm_Ti / Tm_To)

        # convective heat transfer
        Q_conv = h_tot * A_LHTES * dT_lm

        # outlet temperature of the liquid coolant
        T_o_coolant = self.T_melt - Tm_To

        # todo check
        eff_LHTES_check = (T_o_coolant - T_i_coolant) / (self.T_melt - T_i_coolant) - eff_LHTES
        Q_conv_check    = m_dot_coolant * self._coolant_cp * (T_o_coolant - T_i_coolant) - Q_conv

        # --------------- PCM melting distance ----------------------------------
        # 0 < r_melt_front < self._LHTES_thickness
        # melt rate: kg/s
        m_dot_melt      = - Q_conv / self.h_sf_eff
        dr_melt_front   = m_dot_melt / (self.rho_eff * A_LHTES)
        # linear rate: m/s
        dr_by_dt        = dr_melt_front / dt

        # melt location update
        r_melt_front = r_melt_front + dr_melt_front

        # ----------- pressure drop & LHTES Power --------------------------------
        # the pressure drop due to friction
        """hydraulic radius: r_h = A_o / P = d_H / 4"""
        delta_P_chan    = 2 * f_coolant * (self._coolant_rho * u_coolant * u_coolant) * (self._subchannel_length / d_H) * self._subchannel_number
        # pump Power
        P_LHTES         = m_dot_coolant * delta_P_chan / (1e3 * self._coolant_rho * eff_PUMP)

        # ----------- LHTES Mass -------------------------------------------------
        # line density of a cross-section
        rho_line = (self._subchannel_rho * (self._subchannel_height * self._subchannel_tw * 2) + self._coolant_rho * A_c_chan
                    + self.rho_eff * (self._subchannel_height * self._LHTES_thickness)) * self._subchannel_number
        # total mass
        m_LHTES     = rho_line * self._subchannel_length

        return A_LHTES, P_LHTES, m_LHTES, r_melt_front, T_o_coolant, Q_conv, eff_LHTES_check, Q_conv_check

    def SimplifiedLHTES(self, m_dot_coolant=None, Q_wavychannel=None, T_i_coolant=None, dt=None, Heat_accumulated=None, Q_evap_design=1e4):
        T_o_coolant         = T_i_coolant + (Q_evap_design - Q_wavychannel) / (m_dot_coolant * self._coolant_cp)
        Heat_accumulated    = Heat_accumulated + (Q_wavychannel - Q_evap_design) * dt
        self.LHTES_composite_properties()
        m_LHTES             = Heat_accumulated / self.h_sf_eff

        return T_o_coolant, Heat_accumulated, m_LHTES


if __name__ == '__main__':
    LHTES = LatentHeatThermalEnergyStorage()
    T_i_coolant = 35 + 273.15
    r_melt_front = 0.001
    m_dot_coolant = 1.
    dt = 5
    A_LHTES, P_LHTES, m_LHTES, r_melt_front, T_o_coolant, Q_conv, eff_LHTES_check, Q_conv_check \
        = LHTES.LHTES_fluid_cal(m_dot_coolant=m_dot_coolant, r_melt_front=r_melt_front, T_i_coolant=T_i_coolant, dt=dt)
    print(eff_LHTES_check, Q_conv_check)
    print(A_LHTES, P_LHTES, m_LHTES)
    print(r_melt_front, T_o_coolant, Q_conv)
