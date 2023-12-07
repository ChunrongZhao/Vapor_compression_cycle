from __future__                     import division, print_function, absolute_import
import sys, os

from matplotlib.figure              import Figure
from matplotlib                     import pyplot as plt
from matplotlib.ticker              import ScalarFormatter
import numpy                        as np
from math                           import tan, atan, pi
from CoolProp.CoolProp              import PropsSI


# ----------------------------------------------------------------------------
class PlotsClass():

    def __init__(self,  *args,  **kwds):
        pass

    # T-s diagram
    def TSOverlay(self, Cycle, **kwargs):

        Ref                         = Cycle.Ref

        T_min                       = PropsSI('Tmin', Ref)
        T_max                       = PropsSI('Tmax', Ref)
        p_min                       = PropsSI('pmin', Ref)
        p_max                       = PropsSI('pmax', Ref)

        fig                         = plt.figure(figsize=(6, 4))
        self.axes                   = fig.add_subplot(111)

        T_sat                       = np.linspace(T_min, PropsSI(Ref, "Tcrit") - 0.0000000001, 1000)
        (s_sat_L, s_sat_V)          = (0.0 * T_sat, 0.0 * T_sat)
        h_crit                      = np.linspace(0, 0, 2)
        p_crit                      = np.linspace(0, 0, 2)

        for i in np.arange(len(T_sat)):
            s_sat_L[i]              = PropsSI('S', 'T', T_sat[i], 'Q', 0, Ref) / 1000
            s_sat_V[i]              = PropsSI('S', 'T', T_sat[i], 'Q', 1, Ref) / 1000

        # Create critical iso-therm
        p_critical                  = Cycle.AS.p_critical()
        T_critical                  = Cycle.AS.T_critical()

        T_lim                       = np.linspace(T_min+10, T_max, 10000)
        s_lim                       = 0.0 * T_lim

        for i in np.arange(len(T_lim)):
            s_lim[i]                = PropsSI('S', 'T', T_lim[i], 'P', p_critical * 0.999, Ref) / 1000

        # -------------------------------------------------------------------------------------
        # Check if the cycle is transcritical
        # -------------------------------------------------------------------------------------
        if Cycle.Compressor.T_out_r > T_critical and Cycle.Compressor.p_out_r > p_critical:
            # Average saturated temperatures and pressures
            T_sat_evap              = Cycle.T_dew_evap
            p_sat_evap              = PropsSI('P', 'T', T_sat_evap, 'Q', 1, Ref)

            # Saturated conditions at the evaporator
            s_dew_evap              = PropsSI('S', 'P', p_sat_evap, 'Q', 1, Ref)
            T_dew_evap              = PropsSI('T', 'P', p_sat_evap, 'Q', 1, Ref)
            s_bubble_evap           = PropsSI('S', 'P', p_sat_evap, 'Q', 0, Ref)
            T_bubble_evap           = PropsSI('T', 'P', p_sat_evap, 'Q', 0, Ref)

            T_evap_range1           = np.linspace(T_min, T_bubble_evap, 1000)
            T_evap_range2           = np.linspace(T_dew_evap, T_max, 1000)
            s_evap_1                = 0.0 * T_evap_range1
            s_evap_2                = 0.0 * T_evap_range2

            for i in range(len(T_evap_range1)):
                try:
                    s_evap_1[i]     = PropsSI('S', 'P', p_sat_evap, 'T', T_evap_range1[i], Ref) / 1000
                except:
                    s_evap_1[i]     = PropsSI('S', 'P', p_sat_evap, 'Q', 0, Ref) / 1000

                try:
                    s_evap_2[i]     = PropsSI('S', 'P', p_sat_evap, 'T', T_evap_range2[i], Ref) / 1000
                except:
                    s_evap_2[i]     = PropsSI('S', 'P', p_sat_evap, 'Q', 1, Ref) / 1000

        else:
            # Average saturated temperatures and pressures
            T_sat_evap              = Cycle.T_dew_evap
            T_sat_cond              = Cycle.T_dew_cond
            p_sat_evap              = PropsSI('P', 'T', T_sat_evap, 'Q', 1, Ref)
            p_sat_cond              = PropsSI('P', 'T', T_sat_cond, 'Q', 1, Ref)

            # Saturated conditions at the evaporator
            s_dew_evap              = PropsSI('S', 'P', p_sat_evap, 'Q', 1, Ref)
            T_dew_evap              = PropsSI('T', 'P', p_sat_evap, 'Q', 1, Ref)
            s_bubble_evap           = PropsSI('S', 'P', p_sat_evap, 'Q', 0, Ref)
            T_bubble_evap           = PropsSI('T', 'P', p_sat_evap, 'Q', 0, Ref)

            # Saturated conditions at the condenser
            s_dew_cond              = PropsSI('S', 'P', p_sat_cond, 'Q', 1, Ref)
            T_dew_cond              = PropsSI('T', 'P', p_sat_cond, 'Q', 1, Ref)
            s_bubble_cond           = PropsSI('S', 'P', p_sat_cond, 'Q', 0, Ref)
            T_bubble_cond           = PropsSI('T', 'P', p_sat_cond, 'Q', 0, Ref)

            T_cond_range1           = np.linspace(T_min, T_bubble_cond, 1000)
            T_cond_range2           = np.linspace(T_dew_cond, T_max, 1000)
            s_cond_1                = 0.0 * T_cond_range1
            s_cond_2                = 0.0 * T_cond_range2
            T_evap_range1           = np.linspace(T_min, T_bubble_evap, 1000)
            T_evap_range2           = np.linspace(T_dew_evap, T_max, 1000)
            s_evap_1                = 0.0 * T_evap_range1
            s_evap_2                = 0.0 * T_evap_range2

            # ---------------------------------------------------------------------------------
            for i in range(len(T_cond_range1)):
                try:
                    s_cond_1[i]     = PropsSI('S', 'P', p_sat_cond, 'T', T_cond_range1[i], Ref) / 1000
                except:
                    s_cond_1[i]     = PropsSI('S', 'P', p_sat_cond, 'Q', 0, Ref) / 1000

                try:
                    s_cond_2[i]     = PropsSI('S', 'P', p_sat_cond, 'T', T_cond_range2[i], Ref) / 1000
                except:
                    s_cond_2[i]     = PropsSI('S', 'P', p_sat_cond, 'Q', 1, Ref) / 1000

                try:
                    s_evap_1[i]     = PropsSI('S', 'P', p_sat_evap, 'T', T_evap_range1[i], Ref) / 1000
                except:
                    s_evap_1[i]     = PropsSI('S', 'P', p_sat_evap, 'Q', 0, Ref) / 1000

                try:
                    s_evap_2[i]     = PropsSI('S', 'P', p_sat_evap, 'T', T_evap_range2[i], Ref) / 1000
                except:
                    s_evap_2[i]     = PropsSI('S', 'P', p_sat_evap, 'Q', 1, Ref) / 1000

            # Isolines condensing temperature
            self.axes.plot(s_cond_1, T_cond_range1, 'k-.', lw=1, alpha=0.3)
            self.axes.plot(s_cond_2, T_cond_range2, 'k-.', lw=1, alpha=0.3)
            self.axes.plot([s_dew_cond / 1000, s_bubble_cond / 1000], [T_dew_cond, T_bubble_cond], 'k-.', lw=1, alpha=0.3)

        # Critical isothermal
        self.axes.plot(s_lim, T_lim, 'k-.', lw=1, alpha=0.3)

        # Ts saturated lines
        self.axes.plot(s_sat_L, T_sat, 'k', lw=1.5, alpha=0.5)
        self.axes.plot(s_sat_V, T_sat, 'k', lw=1.5, alpha=0.5)

        # Isolines evaporating temperature
        self.axes.plot(s_evap_1, T_evap_range1, 'k-.', lw=1, alpha=0.3)
        self.axes.plot(s_evap_2, T_evap_range2, 'k-.', lw=1, alpha=0.3)
        self.axes.plot([s_dew_evap / 1000, s_bubble_evap / 1000], [T_dew_evap, T_bubble_evap], 'k-.', lw=1, alpha=0.3)

        # Axes labels
        self.axes.set_xlabel('Entropy [kJ/(kg-K)]')
        self.axes.set_ylabel('Temperature [K]')

        # -------------------------------------------------------------------------------------
        # Check if the cycle is transcritical
        # -------------------------------------------------------------------------------------
        if Cycle.Compressor.T_out_r > T_critical  and Cycle.Compressor.p_out_r > p_critical:
            # Axes limits
            self.axes.set_xlim(0.9 * Cycle.GasCooler.s_out_r / 1000, 1.1 * Cycle.Compressor.s_out_r / 1000)
            self.axes.set_ylim(0.95 * Cycle.Evaporator.T_out_r, 1.1 * Cycle.Compressor.T_out_r)

        else:
            # Axes limits
            self.axes.set_xlim(0.7 * Cycle.Condenser.s_out_r / 1000, 1.2 * Cycle.Compressor.s_out_r / 1000)
            self.axes.set_ylim(0.9 * Cycle.Evaporator.T_out_r, 1.1 * Cycle.Compressor.T_out_r)

        # VI cycles compressor
        if Cycle.CycleType == '1INJ_Econ' or Cycle.CycleType == '1INJ_FlashTank':
            s_comp = np.r_[Cycle.Compressor.s_in_r / 1000., Cycle.Compressor.s_31 / 1000, Cycle.Compressor.s_inj_r / 1000., Cycle.Compressor.s_out_r / 1000.]
            T_comp = np.r_[Cycle.Compressor.T_in_r, Cycle.Compressor.T_31, Cycle.Compressor.T_inj_r, Cycle.Compressor.T_out_r]
            self.axes.plot(s_comp, T_comp, 'k', lw=2)
            self.axes.plot(s_comp[0], T_comp[0], 'ko', mfc='w')
            self.axes.plot(s_comp[1], T_comp[1], 'ko', mfc='w')
            self.axes.plot(s_comp[2], T_comp[2], 'ko', mfc='w')
            self.axes.plot(s_comp[3], T_comp[3], 'ko', mfc='w')
            self.axes.text(s_comp[0], T_comp[0], ' 1', ha='left', va='top')
            # self.axes.text(h_comp[1],p_comp[1],' 11',ha='left',va='bottom')
            self.axes.text(s_comp[2], T_comp[2], ' 1inj', ha='right', va='top')
            self.axes.text(s_comp[3], T_comp[3], ' 2', ha='left', va='bottom')

            s_s_comp = np.r_[Cycle.Compressor.s_in_r / 1000., Cycle.Compressor.s_31s / 1000, Cycle.Compressor.s_41s / 1000]
            T_s_comp = np.r_[Cycle.Compressor.T_in_r, Cycle.Compressor.T_inj_r, Cycle.Compressor.T_out_r]
            # self.axes.plot(s_s_comp,T_s_comp,'r--')

        else:
            s_comp = np.r_[Cycle.Compressor.s_in_r / 1000., Cycle.Compressor.s_out_r / 1000.]
            T_comp = np.r_[Cycle.Compressor.T_in_r, Cycle.Compressor.T_out_r]
            self.axes.plot(s_comp, T_comp, 'k', lw=2)
            self.axes.plot(s_comp[0], T_comp[0], 'ko', mfc='w')
            self.axes.plot(s_comp[1], T_comp[1], 'ko', mfc='w')
            self.axes.text(s_comp[0], T_comp[0], ' 1', ha='left', va='top')
            self.axes.text(s_comp[1], T_comp[1], ' 2', ha='left', va='bottom')

        # -------------------------------------------------------------------------------------
        if (Cycle.CycleType == 'Secondary' and Cycle.Mode == 'AC') or Cycle.CycleType == 'DX':
            # DX systems and secondary loop in cooling mode have condenser
            sL = PropsSI('S', 'T', Cycle.T_dew_cond, 'Q', 0, Ref) / 1000
            sV = PropsSI('S', 'T', Cycle.T_dew_cond, 'Q', 1, Ref) / 1000
            s_cond = np.r_[Cycle.Condenser.s_in_r / 1000., s_dew_cond / 1000, s_bubble_cond / 1000, Cycle.Condenser.s_out_r / 1000.]
            T_cond = np.r_[Cycle.Condenser.T_in_r, T_dew_cond, T_bubble_cond, Cycle.Condenser.T_out_r]
            self.axes.plot(s_cond, T_cond, 'k', lw=2)
            self.axes.plot(s_cond[3], T_cond[3], 'ko', mfc='w')
            self.axes.text(s_cond[3], T_cond[3], '3$\quad\quad$', ha='right', va='bottom')

            self.axes.plot([s_cond[3], s_cond[0]], [Cycle.Condenser.T_in_a, Cycle.Condenser.T_out_a], 'r-')

            self.axes.text(0.5 * s_cond[0] + 0.5 * s_cond[3], Cycle.Condenser.T_in_a, 'Outdoor Air', backgroundcolor='w', ha='center', va='center')

        elif Cycle.CycleType == 'CO2-DX' or Cycle.CycleType == 'CO2-LSHX-DX' or Cycle.CycleType == 'Transcrit-DX':
            # Transcritical cycles
            s_gascool = np.linspace(Cycle.GasCooler.s_in_r, Cycle.GasCooler.s_out_r, 10000, endpoint=True)
            p_gascool = Cycle.GasCooler.p_sat_r

            T_gascool = 0.0 * s_gascool

            for i in np.arange(len(s_gascool)):
                T_gascool[i] = PropsSI('T', 'S', s_gascool[i], 'P', p_gascool, Ref)

            self.axes.plot(s_gascool / 1000, T_gascool, 'k', lw=2)

        elif Cycle.CycleType == '1INJ_Econ' or Cycle.CycleType == '1INJ_FlashTank':
            # VI cycles condenser
            sL      = PropsSI('S', 'T', Cycle.T_dew_cond, 'Q', 0, Ref) / 1000
            sV      = PropsSI('S', 'T', Cycle.T_dew_cond, 'Q', 1, Ref) / 1000
            s_cond  = np.r_[Cycle.Condenser.s_in_r/1000., s_dew_cond / 1000, s_bubble_cond / 1000, Cycle.Condenser.s_out_r/1000.]
            T_cond  = np.r_[Cycle.Condenser.T_in_r, T_dew_cond, T_bubble_cond, Cycle.Condenser.T_out_r]
            self.axes.plot(s_cond, T_cond, 'k', lw=2)
            self.axes.plot(s_cond[3], T_cond[3], 'ko', mfc='w')
            # self.axes.text(s_cond[3],T_cond[3],'3$\quad\quad$',ha='right',va='bottom')

            self.axes.plot([s_cond[3], s_cond[0]], [Cycle.Condenser.T_in_a, Cycle.Condenser.T_out_a], 'r-')

            self.axes.text(0.5 * s_cond[0] + 0.5 * s_cond[3], Cycle.Condenser.T_in_a, 'Outdoor Air', backgroundcolor='w', a='center', va='center')

        elif Cycle.CycleType == 'Secondary' and Cycle.Mode == 'HP':
            sV      = PropsSI('S', 'T', Cycle.T_dew_evap, 'Q', 1, Ref) / 1000
            s_evap  = np.r_[Cycle.Evaporator.s_in_r / 1000., sV, Cycle.Evaporator.s_out_r / 1000.]
            T_evap  = np.r_[Cycle.Evaporator.T_in_r, Cycle.T_dew_evap, Cycle.Evaporator.T_out_r]
            self.axes.plot(s_evap, T_evap, 'k', lw=2)
            self.axes.plot(s_evap[0], T_evap[0], 'ko', mfc='w')
            self.axes.text(s_evap[0], T_evap[0], '4$\quad\quad$', ha='right', va='bottom')

            self.axes.plot([s_evap[2], s_evap[0]], [Cycle.Evaporator.T_in_a, Cycle.Evaporator.T_out_a], 'r-')

            self.axes.text(0.5 * s_evap[0] + 0.5 * s_evap[2], Cycle.Evaporator.T_in_a, 'Outdoor Air', backgroundcolor='w', ha='center', va='center')

        # -------------------------------------------------------------------------------------
        if Cycle.CycleType == "Secondary":
            if Cycle.IHXType == 'Coaxial':
                IHX                 = Cycle.CoaxialIHX
                sV                  = PropsSI('S', 'T', Cycle.T_dew_evap, 'Q', 1, Ref) / 1000
                s_IHX               = np.r_[IHX.s_in_r / 1000., sV, IHX.s_out_r / 1000.]
                T_IHX               = np.r_[IHX.T_in_r, Cycle.T_dew_evap, IHX.T_out_r]
                s_XV                = np.r_[IHX.s_in_r / 1000., Cycle.Condenser.s_out_r / 1000.]
                T_XV                = np.r_[IHX.T_in_r, Cycle.Condenser.T_out_r]
                T_IHXg              = np.r_[IHX.T_out_g, IHX.T_in_g]
                s_IHXg              = np.r_[s_IHX[0], s_IHX[len(s_IHX)-1]]

            elif Cycle.IHXType == 'PHE':
                if Cycle.Mode == 'AC':
                    IHX             = Cycle.PHEIHX
                    sV              = PropsSI('S', 'T', Cycle.T_dew_evap, 'Q', 1, Ref) / 1000
                    s_IHX           = np.r_[IHX.s_in_c / 1000., sV, IHX.s_out_c / 1000.]
                    T_IHX           = np.r_[IHX.T_in_c, Cycle.T_dew_evap, IHX.T_out_c]
                    s_XV            = np.r_[IHX.s_in_c / 1000., Cycle.Condenser.s_out_r / 1000.]
                    T_XV            = np.r_[IHX.T_in_c, Cycle.Condenser.T_out_r]
                    T_IHXg          = np.r_[IHX.T_out_h, IHX.T_in_h]
                    s_IHXg          = np.r_[s_IHX[0], s_IHX[len(s_IHX)-1]]

                    self.axes.text(s_IHX[0], T_IHX[0], '4$\quad\quad$', ha='right', va='top')
                    self.axes.plot(s_IHX[0], T_IHX[0], 'ko', mfc='w')
                else:
                    IHX = Cycle.PHEIHX
                    sL              = PropsSI('S', 'T', Cycle.T_bubble_cond, 'Q', 0, Ref) / 1000
                    sV              = PropsSI('S', 'T', Cycle.T_dew_cond, 'Q', 1, Ref) / 1000
                    s_IHX           = np.r_[IHX.s_in_h / 1000., sV, sL, IHX.s_out_h / 1000.]
                    T_IHX           = np.r_[IHX.T_in_h, Cycle.T_dew_cond, Cycle.T_bubble_cond, IHX.T_out_h]
                    s_XV            = np.r_[Cycle.Evaporator.s_in_r / 1000., IHX.s_out_h / 1000.]
                    T_XV            = np.r_[Cycle.Evaporator.T_in_r, IHX.T_out_h]
                    T_IHXg          = np.r_[IHX.T_out_c, IHX.T_in_c]
                    s_IHXg          = np.r_[s_IHX[0], s_IHX[len(s_IHX)-1]]

                    self.axes.text(s_IHX[3], T_IHX[3], '3$\quad\quad$', ha='right' ,va='top')
                    self.axes.plot(s_IHX[3], T_IHX[3], 'ko', mfc='w')
            else:
                raise ValueError('Secondary loop system must have a coaxial or PHE heat exchanger')

            self.axes.plot(s_IHX, T_IHX, 'k', lw=2)

            self.axes.plot(s_XV, T_XV, 'k', lw=2)
            self.axes.plot([s_IHX[2], s_IHX[0]], [Cycle.CoolingCoil.T_in_a, Cycle.CoolingCoil.T_out_a], 'b-')
            self.axes.text(0.5 * s_IHX[2] + 0.5 * s_IHX[0], Cycle.CoolingCoil.T_in_a, 'Indoor Air', backgroundcolor='w', ha='center', va='center')
            self.axes.text(0.5 * s_IHXg[0] + 0.5 * s_IHXg[1], np.mean(T_IHXg), 'Glycol IHX', backgroundcolor='w', ha='center', va='center')
            self.axes.plot(s_IHXg, T_IHXg, 'g-.')

        elif Cycle.CycleType == "DX":
            sL                  = PropsSI('S', 'T', Cycle.T_dew_evap, 'Q', 0, Ref) / 1000
            sV                  = PropsSI('S', 'T', Cycle.T_dew_evap, 'Q', 1, Ref) / 1000
            s_evap              = np.r_[Cycle.Evaporator.s_in_r / 1000., sV, Cycle.Evaporator.s_out_r / 1000.]
            T_evap              = np.r_[Cycle.Evaporator.T_sat_r, Cycle.T_dew_evap, Cycle.Evaporator.T_out_r]
            self.axes.plot(s_evap, T_evap, 'k', lw=2)
            self.axes.plot(s_evap[0], T_evap[0], 'ko', mfc='w')
            self.axes.text(s_evap[0], T_evap[0], '4$\quad\quad$', ha='right', va='top')

            s_XV                = np.r_[Cycle.Evaporator.s_in_r / 1000., Cycle.Condenser.s_out_r / 1000.]
            T_XV                = np.r_[Cycle.Evaporator.T_sat_r, Cycle.Condenser.T_out_r]
            self.axes.plot(s_XV, T_XV, 'k', lw=2)

            self.axes.plot([s_evap[2], s_evap[0]], [Cycle.Evaporator.T_in_a, Cycle.Evaporator.T_out_a], 'b-')
            self.axes.text(0.5 * s_evap[2] + 0.5 * s_evap[0], Cycle.Evaporator.T_in_a, 'Indoor Air', backgroundcolor='w', ha='center', va='center')

        elif Cycle.CycleType == 'CO2-DX' or Cycle.CycleType == 'CO2-LSHX-DX' or Cycle.CycleType == 'Transcrit-DX':
            sL                      = PropsSI('S', 'T', Cycle.T_dew_evap, 'Q', 0, Ref) / 1000
            sV                      = PropsSI('S', 'T', Cycle.T_dew_evap, 'Q', 1, Ref) / 1000
            s_evap                  = np.r_[Cycle.Evaporator.s_in_r / 1000., sV, Cycle.Evaporator.s_out_r / 1000.]
            T_evap                  = np.r_[Cycle.Evaporator.T_sat_r, Cycle.Tdew_evap, Cycle.Evaporator.T_out_r]
            self.axes.plot(s_evap, T_evap, 'k', lw=2)
            self.axes.plot(s_evap[0], T_evap[0], 'ko', mfc='w')
            self.axes.text(s_evap[0], T_evap[0], '4$\quad\quad$', ha='right', va='top')

            s_XV                    = np.r_[Cycle.Evaporator.s_in_r / 1000., Cycle.GasCooler.s_out_r / 1000.]
            T_XV                    = np.r_[Cycle.Evaporator.T_sat_r, Cycle.GasCooler.T_out_r]
            self.axes.plot(s_XV, T_XV, 'k', lw=2)

            self.axes.plot([s_evap[2], s_evap[0]], [Cycle.Evaporator.T_in_a, Cycle.Evaporator.T_out_a], 'b-')
            self.axes.text(0.5 * s_evap[2] + 0.5 * s_evap[0], Cycle.Evaporator.T_in_a, 'Indoor Air', backgroundcolor='w', ha='center', va='center')

        elif Cycle.CycleType == '1INJ_Econ':

            s_evap                  = np.r_[Cycle.Evaporator.s_in_r / 1000., s_dew_evap / 1000, Cycle.Evaporator.s_out_r / 1000.]
            T_evap                  = np.r_[Cycle.Evaporator.T_in_r, Cycle.T_dew_evap, Cycle.Evaporator.T_out_r]
            self.axes.plot(s_evap, T_evap, 'k', lw=2)
            self.axes.plot(s_evap[0], T_evap[0], 'ko', mfc='w')
            self.axes.text(s_evap[0], T_evap[0], '7$\quad\quad$', ha='right', va='top')

            s_PHEHX                 = np.r_[Cycle.PHEHX.s_out_h / 1000., Cycle.Condenser.s_out_r / 1000.]
            T_PHEHX                 = np.r_[Cycle.PHEHX.T_out_h, Cycle.Condenser.T_out_r]
            self.axes.plot(s_PHEHX, T_PHEHX, 'k', lw=2)
            self.axes.plot(s_PHEHX[0], T_PHEHX[0], 'ko', mfc='w')
            self.axes.text(s_PHEHX[0], T_PHEHX[0], '4$\quad\quad$', ha='right', va='bottom')
            s_XV                    = np.r_[Cycle.Evaporator.s_in_r / 1000., Cycle.PHEHX.s_out_h / 1000.]
            T_XV                    = np.r_[Cycle.Evaporator.T_in_r, Cycle.PHEHX.T_out_h]
            self.axes.plot(s_XV, T_XV, 'k', lw=2)

            # h_XV2=np.r_[Cycle.Condenser.hout_r/1000.,Cycle.Condenser.hout_r/1000.]
            # p_XV2=np.r_[Cycle.Condenser.psat_r/1000.,Cycle.PHEHX.pin_c/1000.]
            # self.axes.plot(h_XV[1],p_XV[1],'bo')
            # self.axes.plot(h_XV2,p_XV2,'b')

            s_dew_inj   = PropsSI('S', 'P', Cycle.Compressor.p_inj_r, 'Q', 1, Cycle.Ref)
            T_dew_inj   = PropsSI('T', 'P', Cycle.Compressor.p_inj_r, 'Q', 1, Cycle.Ref)

            s_1INJ      = np.r_[Cycle.PHEHX.s_out_h / 1000, s_dew_inj / 1000, Cycle.PHEHX.s_out_c / 1000.]
            T_1INJ      = np.r_[Cycle.PHEHX.T_in_c, T_dew_inj, Cycle.PHEHX.T_out_c]
            self.axes.plot(s_1INJ, T_1INJ, 'k', lw=2)
            self.axes.plot(s_1INJ[0], T_1INJ[0], 'ko', mfc='w')
            self.axes.text(s_1INJ[0], T_1INJ[0], '5$\quad\quad$', ha='right', va='bottom')

            self.axes.plot([s_evap[2], s_evap[0]], [Cycle.Evaporator.T_in_a, Cycle.Evaporator.T_out_a], 'b-')
            self.axes.text(0.5 * s_evap[2] + 0.5 * s_evap[0], Cycle.Evaporator.T_in_a, 'Indoor Air', backgroundcolor='w', ha='center', va='center')

        elif Cycle.CycleType == '1INJ_FlashTank':

            s_evap          = np.r_[Cycle.Evaporator.s_in_r / 1000., s_dew_evap / 1000, Cycle.Evaporator.s_out_r / 1000.]
            T_evap          = np.r_[Cycle.Evaporator.T_in_r, T_dew_evap, Cycle.Evaporator.T_out_r]
            self.axes.plot(s_evap, T_evap, 'k', lw=2)
            self.axes.plot(s_evap[0], T_evap[0], 'ko', mfc='w')
            self.axes.text(s_evap[0], T_evap[0], '9$\quad\quad$', ha='right', va='top')

            s_XV            = np.r_[Cycle.FlashTank.s_in / 1000., Cycle.Condenser.s_out_r / 1000.]
            T_XV            = np.r_[Cycle.FlashTank.T_in, Cycle.Condenser.T_out_r]
            self.axes.plot(s_XV, T_XV, 'k', lw=2)

            s_FlashTank_liquid = np.r_[Cycle.FlashTank.s_in / 1000., Cycle.FlashTank.s_out / 1000.]
            T_FlashTank_liquid = np.r_[Cycle.FlashTank.T_in, Cycle.FlashTank.T_out]
            self.axes.plot(s_FlashTank_liquid, T_FlashTank_liquid, 'k', lw=2)
            self.axes.plot(s_FlashTank_liquid[0], T_FlashTank_liquid[0], 'ko', mfc='w')
            self.axes.text(s_FlashTank_liquid[0], T_FlashTank_liquid[0], '6$\quad\quad$', ha='right', va='top')

            s_1INJ              = np.r_[Cycle.FlashTank.s_in / 1000, Cycle.Compressor.s_inj_r / 1000.]
            T_1INJ              = np.r_[Cycle.FlashTank.T_in, Cycle.Compressor.T_inj_r]
            self.axes.plot(s_1INJ, T_1INJ, 'k', lw=2)
            self.axes.plot(s_1INJ[0], T_1INJ[0], 'ro', mfc='w')
            self.axes.text(s_1INJ[0], T_1INJ[0], '8$\quad\quad$', ha='right', va='bottom')

            s_XV2               = np.r_[Cycle.FlashTank.s_out / 1000., Cycle.Evaporator.s_in_r / 1000]
            T_XV2               = np.r_[Cycle.FlashTank.T_in, Cycle.Evaporator.T_sat_r]
            self.axes.plot(s_XV[1], T_XV[1], 'ko', mfc='w')
            self.axes.plot(s_XV2, T_XV2, 'k', lw=2)
            self.axes.plot(s_XV2[0], T_XV2[0], 'ko', mfc='w')
            self.axes.text(s_XV2[0], T_XV2[0], '8$\quad\quad$', ha='right', va='top')

            self.axes.plot([s_evap[2], s_evap[0]], [Cycle.Evaporator.T_in_a, Cycle.Evaporator.T_out_a], 'b-')
            self.axes.text(0.5 * s_evap[2] + 0.5 * s_evap[0], Cycle.Evaporator.T_in_a, 'Indoor Air', backgroundcolor='w', ha='center', va='center')

        plt.show()

    # p-h diagram
    def PHOverlay(self, Cycle, **kwargs):

        Ref                             = Cycle.Ref

        T_min                           = PropsSI('Tmin', Ref)
        T_max                           = PropsSI('Tmax', Ref)
        p_min                           = PropsSI('pmin', Ref)
        p_max                           = PropsSI('pmax', Ref)

        fig                             = plt.figure(figsize=(6,4))
        self.axes                       = fig.add_subplot(111)

        # Create saturated lines
        T_sat                           = np.linspace(T_min, PropsSI(Ref, "Tcrit") - 0.0000000001, 1000)
        (h_sat_L, h_sat_V)              = (0.0 * T_sat, 0.0 * T_sat)
        (p_sat_L, p_sat_V)              = (0.0 * T_sat, 0.0 * T_sat)
        for i in np.arange(len(T_sat)):
            p_sat_L[i]                  = PropsSI('P', 'T', T_sat[i], 'Q', 0, Ref) / 1000
            p_sat_V[i]                  = PropsSI('P', 'T', T_sat[i], 'Q', 1, Ref) / 1000
            h_sat_L[i]                  = PropsSI('H', 'T', T_sat[i], 'Q', 0, Ref) / 1000
            h_sat_V[i]                  = PropsSI('H', 'T', T_sat[i], 'Q', 1, Ref) / 1000

        # Create critical iso-therm
        p_critical                      = Cycle.AS.p_critical()
        T_critical                      = Cycle.AS.T_critical()

        p_lim                           = np.linspace(p_min, p_max, 10000)
        h_lim                           = 0.0 * p_lim

        for i in np.arange(len(p_lim)):
            h_lim[i]                    = PropsSI('H', 'P', p_lim[i], 'T', T_critical - 0.0000000001, Ref) / 1000

        # --------------------------------------------------------------------------------
        # Check if the cycle is transcritical
        # --------------------------------------------------------------------------------
        if Cycle.Compressor.T_out_r > T_critical  and Cycle.Compressor.p_out_r > p_critical:
            # Average saturated temperatures and pressures
            T_dew_evap                  = Cycle.T_dew_evap
            p_sat_evap                  = PropsSI('P', 'T', T_dew_evap, 'Q', 1, Ref)

            # Saturated conditions at the evaporator
            h_dew_evap                  = PropsSI('H', 'P', p_sat_evap, 'Q', 1, Ref)
            T_dew_evap                  = PropsSI('T', 'P', p_sat_evap, 'Q', 1, Ref)
            h_bubble_evap               = PropsSI('H', 'P', p_sat_evap, 'Q', 0, Ref)
            T_bubble_evap               = PropsSI('T', 'P', p_sat_evap, 'Q', 0, Ref)

            p_evap_range1               = np.linspace(p_min, p_sat_evap, 1000)
            p_evap_range2               = np.linspace(p_sat_evap, p_max, 1000)
            h_evap_1                    = 0.0 * p_evap_range1
            h_evap_2                    = 0.0 * p_evap_range2

            # --------------------------------------------------------------------------------
            for i in range(len(p_evap_range1)):
                try:
                    h_evap_1[i]         = PropsSI('H', 'T', T_dew_evap, 'P', p_evap_range1[i], Ref) / 1000
                except:
                    h_evap_1[i]         = PropsSI('H', 'P', p_sat_evap, 'Q', 1, Ref) / 1000

                try:
                    h_evap_2[i]         = PropsSI('H', 'T', T_bubble_evap, 'P', p_evap_range2[i], Ref) / 1000
                except:
                    h_evap_2[i]         = PropsSI('H', 'P', p_sat_evap, 'Q', 0, Ref) / 1000

        else:
            # Average saturated temperatures and pressures
            T_dew_evap                  = Cycle.T_dew_evap
            T_dew_cond                  = Cycle.T_dew_cond
            p_sat_evap                  = PropsSI('P', 'T', T_dew_evap, 'Q', 1, Ref)
            p_sat_cond                  = PropsSI('P', 'T', T_dew_cond, 'Q', 1, Ref)

            # Saturated conditions at the evaporator
            h_dew_evap                  = PropsSI('H', 'P', p_sat_evap, 'Q', 1, Ref)
            T_dew_evap                  = PropsSI('T', 'P', p_sat_evap, 'Q', 1, Ref)
            h_bubble_evap               = PropsSI('H', 'P', p_sat_evap, 'Q', 0, Ref)
            T_bubble_evap               = PropsSI('T', 'P', p_sat_evap, 'Q', 0, Ref)

            # Saturated conditions at the condenser
            h_dew_cond                  = PropsSI('H', 'P', p_sat_cond, 'Q', 1, Ref)
            T_dew_cond                  = PropsSI('T', 'P', p_sat_cond, 'Q', 1, Ref)
            h_bubble_cond               = PropsSI('H', 'P', p_sat_cond, 'Q', 0, Ref)
            T_bubble_cond               = PropsSI('T', 'P', p_sat_cond, 'Q', 0, Ref)

            p_cond_range1               = np.linspace(p_min, p_sat_cond, 1000)
            p_cond_range2               = np.linspace(p_sat_cond, p_max, 1000)
            h_cond_1                    = 0.0 * p_cond_range1
            h_cond_2                    = 0.0 * p_cond_range2
            p_evap_range1               = np.linspace(p_min, p_sat_evap, 1000)
            p_evap_range2               = np.linspace(p_sat_evap, p_max, 1000)
            h_evap_1                    = 0.0 * p_evap_range1
            h_evap_2                    = 0.0 * p_evap_range2

            for i in range(len(p_cond_range1)):
                try:
                    h_cond_1[i]         = PropsSI('H', 'T', T_dew_cond, 'P', p_cond_range1[i], Ref) / 1000
                except:
                    h_cond_1[i]         = PropsSI('H', 'P', p_sat_cond, 'Q', 1, Ref) / 1000

                try:
                    h_cond_2[i]         = PropsSI('H', 'T', T_bubble_cond, 'P', p_cond_range2[i], Ref) / 1000
                except:
                    h_cond_2[i]         = PropsSI('H', 'P', p_sat_cond, 'Q', 0, Ref) / 1000

                try:
                    h_evap_1[i]         = PropsSI('H', 'T', T_dew_evap, 'P', p_evap_range1[i], Ref) / 1000
                except:
                    h_evap_1[i]         = PropsSI('H', 'P', p_sat_evap, 'Q', 1, Ref) / 1000

                try:
                    h_evap_2[i]         = PropsSI('H', 'T', T_bubble_evap, 'P', p_evap_range2[i], Ref) / 1000
                except:
                    h_evap_2[i]         = PropsSI('H', 'P', p_sat_evap, 'Q', 0, Ref) / 1000

            # Isolines condensing temperature
            self.axes.semilogy(h_cond_1, p_cond_range1 / 1000, 'k-.', lw=1, alpha=0.3)
            self.axes.semilogy(h_cond_2, p_cond_range2 / 1000, 'k-.', lw=1, alpha=0.3)
            self.axes.semilogy([h_dew_cond / 1000, h_bubble_cond / 1000], [p_sat_cond / 1000, p_sat_cond / 1000], 'k-.', lw=1, alpha=0.3)

        # Isolines evaporating temperature
        self.axes.semilogy(h_evap_1, p_evap_range1 / 1000, 'k-.', lw=1, alpha=0.3)
        self.axes.semilogy(h_evap_2, p_evap_range2 / 1000, 'k-.', lw=1, alpha=0.3)
        self.axes.semilogy([h_dew_evap / 1000, h_bubble_evap / 1000], [p_sat_evap / 1000, p_sat_evap / 1000], 'k-.', lw=1, alpha=0.3)

        # Critical isothermal
        self.axes.semilogy(h_lim, p_lim / 1000, 'k-.', lw=1, alpha=0.3)

        # ph saturated line plot
        self.axes.semilogy(h_sat_L, p_sat_L, 'k', lw=1.5, alpha=0.5)
        self.axes.semilogy(h_sat_V, p_sat_V, 'k', lw=1.5, alpha=0.5)

        # Axes labels
        self.axes.set_xlabel('Enthalpy [kJ/kg]')
        self.axes.set_ylabel('Pressure [kPa]')

        # --------------------------------------------------------------------------------
        # Check if the cycle is transcritical
        # --------------------------------------------------------------------------------
        if Cycle.Compressor.T_out_r > T_critical  and Cycle.Compressor.p_out_r > p_critical:
            # Axes limits
            self.axes.set_xlim(0.7 * Cycle.GasCooler.h_out_r / 1000, 1.1 * Cycle.Compressor.h_out_r / 1000)
            self.axes.set_ylim(10 * p_min / 1000, 2 * Cycle.GasCooler.p_sat_r / 1000)

        else:
            # Axes limits
            self.axes.set_xlim(0.7 * Cycle.Condenser.h_out_r / 1000, 1.2 * Cycle.Compressor.h_out_r / 1000)
            self.axes.set_ylim(10 * p_min / 1000, 0.5 * p_max / 1000)

        # Check cycle type to plot the cycle
        if Cycle.CycleType == '1INJ_Econ' or Cycle.CycleType == '1INJ_FlashTank':
            h_comp  = np.r_[Cycle.Compressor.h_in_r / 1000., Cycle.Compressor.h_31 / 1000, Cycle.Compressor.h_inj_r / 1000., Cycle.Compressor.h_out_r / 1000.]
            p_comp  = np.r_[Cycle.Compressor.p_in_r / 1000, Cycle.Compressor.p_inj_r / 1000, Cycle.Compressor.p_inj_r / 1000, Cycle.Compressor.p_out_r / 1000]
            self.axes.semilogy(h_comp, p_comp, 'k', lw=2)
            self.axes.semilogy(h_comp[0], p_comp[0], 'ko', mfc='w')
            self.axes.semilogy(h_comp[1], p_comp[1], 'ko', mfc='w')
            self.axes.semilogy(h_comp[2], p_comp[2], 'ko', mfc='w')
            self.axes.semilogy(h_comp[3], p_comp[3], 'ko', mfc='w')
            self.axes.text(h_comp[0], p_comp[0], ' 1', ha='left', va='top')
            # self.axes.text(h_comp[1],p_comp[1],' 11',ha='left',va='bottom')
            self.axes.text(h_comp[2], p_comp[2], ' 1inj', ha='left', va='bottom')
            self.axes.text(h_comp[3], p_comp[3], ' 2', ha='left', va='bottom')

            h_s_comp    = np.r_[Cycle.Compressor.h_in_r / 1000., Cycle.Compressor.h_31s / 1000, Cycle.Compressor.h_41s / 1000]
            p_s_comp    = np.r_[Cycle.Compressor.p_in_r / 1000., Cycle.Compressor.p_inj_r / 1000, Cycle.Compressor.p_out_r / 1000.]

            # self.axes.plot(h_s_comp,p_s_comp,'r--')
            # self.axes.plot(Cycle.Compressor.h_31/1000,Cycle.Compressor.pinj_r/1000,'bo')
            # self.axes.plot(Cycle.Compressor.h_31s/1000,Cycle.Compressor.pinj_r/1000,'bo')
            # self.axes.plot(Cycle.Compressor.h_22s/1000,Cycle.Compressor.pinj_r/1000,'bo')
            # self.axes.plot(Cycle.Compressor.h_22/1000,Cycle.Compressor.pinj_r/1000,'bo')
            # self.axes.plot(Cycle.Compressor.h_32s/1000,Cycle.Compressor.pout_r/1000,'bo')
            # self.axes.plot(Cycle.Compressor.h_41s/1000,Cycle.Compressor.pout_r/1000,'ro')

        else:
            h_comp      = np.r_[Cycle.Compressor.h_in_r / 1000., Cycle.Compressor.h_out_r / 1000.]
            p_comp      = np.r_[Cycle.Compressor.p_in_r / 1000, Cycle.Compressor.p_out_r / 1000]
            self.axes.semilogy(h_comp, p_comp, 'k', lw=2)
            self.axes.semilogy(h_comp[0], p_comp[0], 'ko', mfc='w')
            self.axes.semilogy(h_comp[1], p_comp[1], 'ko', mfc='w')
            self.axes.text(h_comp[0], p_comp[0], ' 1', ha='left', va='top')
            self.axes.text(h_comp[1], p_comp[1], ' 2', ha='left', va='bottom')

        # --------------------------------------------------------------------------------
        if Cycle.CycleType == 'Secondary' and Cycle.Mode == 'HP':
            h_evap      = np.r_[Cycle.Evaporator.h_in_r / 1000., Cycle.Evaporator.h_out_r / 1000]
            p_evap      = np.r_[Cycle.Evaporator.p_sat_r / 1000, Cycle.Evaporator.p_sat_r / 1000]
            self.axes.semilogy(h_evap, p_evap, 'k', lw=2)
            self.axes.semilogy(h_evap[0], p_evap[0], 'ko', mfc='w')
            self.axes.text(h_evap[0], p_evap[0], '4$\quad\quad$', ha='right', va='bottom')

        elif Cycle.CycleType == '1INJ_Econ' or Cycle.CycleType == '1INJ_FlashTank':
            h_cond      = np.r_[Cycle.Condenser.h_in_r / 1000., Cycle.Condenser.h_out_r / 1000.]
            p_cond      = np.r_[Cycle.Compressor.p_out_r / 1000, Cycle.Compressor.p_out_r / 1000]
            self.axes.semilogy(h_cond, p_cond, 'k', lw=2)
            self.axes.semilogy(h_cond[1], p_cond[1], 'ko', mfc='w')
            # self.axes.semilogy(h_cond[1],p_cond[1],'3$\quad\quad$',ha='right',va='bottom')

        elif Cycle.CycleType == 'CO2-DX' or Cycle.CycleType == 'CO2-LSHX-DX' or Cycle.CycleType == 'Transcrit-DX':
            h_gascool       = np.r_[Cycle.Compressor.h_out_r / 1000., Cycle.GasCooler.h_out_r / 1000.]
            p_gascool       = np.r_[Cycle.GasCooler.p_sat_r / 1000, (Cycle.GasCooler.p_sat_r - Cycle.GasCooler.DP_r) / 1000]
            self.axes.semilogy(h_gascool, p_gascool, 'k', lw=2)
            self.axes.semilogy(h_gascool[1], p_gascool[1], 'ko', mfc='w')
            self.axes.text(h_gascool[1], p_gascool[1], '3$\quad\quad$', ha='right', va='bottom')

        else:
            h_cond      = np.r_[Cycle.Condenser.h_in_r / 1000., Cycle.Condenser.h_out_r / 1000.]
            p_cond      = np.r_[Cycle.Condenser.p_sat_r / 1000, (Cycle.Condenser.p_sat_r - Cycle.Condenser.DP_r) / 1000]
            self.axes.semilogy(h_cond, p_cond, 'k', lw=2)
            self.axes.semilogy(h_cond[1], p_cond[1], 'ko', mfc='w')
            self.axes.text(h_cond[1], p_cond[1], '3$\quad\quad$', ha='right', va='bottom')

        # --------------------------------------------------------------------------------
        if Cycle.CycleType == "Secondary":
            if Cycle.IHXType == 'Coaxial':
                IHX                                 = Cycle.CoaxialIHX
                h_IHX                               = np.r_[IHX.h_in_r / 1000., IHX.h_out_r / 1000.]
                p_IHX                               = np.r_[IHX.p_in_r / 1000, IHX.p_in_r]
                h_XV                                = np.r_[IHX.h_in_r / 1000., Cycle.Condenser.h_out_r / 1000.]
                p_XV                                = np.r_[IHX.p_in_r / 1000, Cycle.Condenser.p_sat_r / 1000]
            elif Cycle.IHXType == 'PHE':
                IHX                                 = Cycle.PHEIHX

                if Cycle.Mode == 'AC':
                    h_IHX                           = np.r_[IHX.h_in_c / 1000., IHX.h_out_c / 1000.]
                    p_IHX                           = np.r_[IHX.p_in_c / 1000, IHX.p_in_c / 1000]
                    h_XV                            = np.r_[IHX.h_in_c / 1000., Cycle.Condenser.h_out_r / 1000.]
                    p_XV                            = np.r_[IHX.p_in_c / 1000, Cycle.Condenser.p_sat_r / 1000]
                    self.axes.semilogy(h_IHX[0], p_IHX[0], 'ko', mfc='w')
                    self.axes.text(h_IHX[0], p_IHX[0], '4$\quad\quad$', ha='right', va='top')
                else:
                    h_IHX                           = np.r_[IHX.h_in_h / 1000., IHX.h_out_h / 1000.]
                    p_IHX                           = np.r_[IHX.p_in_h / 1000, IHX.p_in_h / 1000]
                    h_XV                            = np.r_[Cycle.Evaporator.h_in_r / 1000., IHX.h_out_h / 1000.]
                    p_XV                            = np.r_[Cycle.Evaporator.p_sat_r / 1000, IHX.p_in_h / 1000]
                    self.axes.semilogy(h_IHX[1], p_IHX[1], 'ko', mfc='w')
                    self.axes.text(h_IHX[1], p_IHX[1], '3$\quad\quad$', ha='right', va='top')
            else:
                raise ValueError('Secondary loop system must have a coaxial or PHE heat exchanger')

            self.axes.semilogy(h_IHX, p_IHX, 'k', lw=2)
            self.axes.semilogy(h_XV, p_XV, 'b', lw=2)

        elif Cycle.CycleType == "DX":
            h_evap                      = np.r_[Cycle.Evaporator.h_in_r / 1000., Cycle.Evaporator.h_out_r / 1000.]
            p_evap                      = np.r_[Cycle.Evaporator.p_sat_r / 1000, Cycle.Evaporator.p_sat_r / 1000]
            self.axes.semilogy(h_evap, p_evap, 'k', lw=2)
            self.axes.semilogy(h_evap[0], p_evap[0], 'ko', mfc='w')
            self.axes.text(h_evap[0], p_evap[0], '4$\quad\quad$', ha='right', va='top')

            h_XV                        = np.r_[Cycle.Evaporator.h_in_r / 1000., Cycle.Condenser.h_out_r / 1000.]
            p_XV                        = np.r_[Cycle.Evaporator.p_sat_r / 1000, Cycle.Condenser.p_sat_r / 1000]
            self.axes.semilogy(h_XV, p_XV, 'k', lw=2)

        elif Cycle.CycleType == 'CO2-DX' or Cycle.CycleType == 'CO2-LSHX-DX' or Cycle.CycleType == 'Transcrit-DX':
            h_evap                      = np.r_[Cycle.Evaporator.h_in_r / 1000., Cycle.Evaporator.h_out_r / 1000.]
            p_evap                      = np.r_[Cycle.Evaporator.p_sat_r / 1000, Cycle.Evaporator.p_sat_r / 1000]
            self.axes.semilogy(h_evap, p_evap, 'k', lw=2)
            self.axes.semilogy(h_evap[0], p_evap[0], 'ko', mfc='w')
            self.axes.text(h_evap[0], p_evap[0], '4$\quad\quad$', ha='right', va='top')

            h_XV                        = np.r_[Cycle.Evaporator.h_in_r / 1000., Cycle.GasCooler.h_out_r / 1000.]
            p_XV                        = np.r_[Cycle.Evaporator.p_sat_r / 1000, Cycle.GasCooler.p_sat_r / 1000]
            self.axes.semilogy(h_XV, p_XV, 'k', lw=2)

        elif Cycle.CycleType == '1INJ_Econ':
            h_evap                      = np.r_[Cycle.Evaporator.h_in_r / 1000., Cycle.Evaporator.h_out_r / 1000.]
            p_evap                      = np.r_[Cycle.Evaporator.p_sat_r / 1000, Cycle.Evaporator.p_sat_r / 1000]
            self.axes.semilogy(h_evap, p_evap, 'k', lw=2)
            self.axes.semilogy(h_evap[0], p_evap[0], 'ko', mfc='w')
            self.axes.text(h_evap[0], p_evap[0], '7$\quad\quad$', ha='right', va='top')

            h_PHEHX = np.r_[Cycle.PHEHX.hout_h/1000.,Cycle.Condenser.hout_r/1000.]
            p_PHEHX = np.r_[Cycle.Compressor.pout_r/1000,Cycle.Compressor.pout_r/1000]
            self.axes.semilogy(h_PHEHX,p_PHEHX,'k', lw = 2)
            self.axes.semilogy(h_PHEHX[0],p_PHEHX[0],'ko', mfc = 'w')
            self.axes.text(h_PHEHX[0],p_PHEHX[0],'4$\quad\quad$',ha='right',va='top')
            h_XV=np.r_[Cycle.Evaporator.hin_r/1000.,Cycle.PHEHX.hout_h/1000.]
            p_XV=np.r_[Cycle.Evaporator.psat_r/1000,(Cycle.PHEHX.pin_h-Cycle.PHEHX.DP_h)/1000.]
            self.axes.plot(h_XV,p_XV,'k', lw = 2)

            h_XV2=np.r_[Cycle.Condenser.hout_r/1000.,Cycle.Condenser.hout_r/1000.]
            p_XV2=np.r_[Cycle.Condenser.psat_r/1000.,Cycle.PHEHX.pin_c/1000.]
            #self.axes.plot(h_XV[1],p_XV[1],'bo')
            #self.axes.plot(h_XV2,p_XV2,'b')

            h_1INJ=np.r_[Cycle.PHEHX.hout_h/1000,Cycle.PHEHX.hout_c/1000.]
            p_1INJ=np.r_[Cycle.PHEHX.pin_c/1000,(Cycle.PHEHX.pin_c-Cycle.PHEHX.DP_c)/1000.]
            self.axes.semilogy(h_1INJ,p_1INJ,'k', lw = 2)
            self.axes.semilogy(h_1INJ[0],p_1INJ[0],'ko', mfc = 'w')
            self.axes.text(h_1INJ[0],p_1INJ[0],'5$\quad\quad$',ha='right',va='top')

        elif Cycle.CycleType == '1INJ_FlashTank':

            h_evap                      = np.r_[Cycle.Evaporator.h_in_r / 1000., Cycle.Evaporator.h_out_r / 1000.]
            p_evap                      = np.r_[Cycle.Evaporator.p_sat_r / 1000, Cycle.Evaporator.p_sat_r / 1000]
            self.axes.semilogy(h_evap, p_evap, 'k', lw=2)
            self.axes.semilogy(h_evap[0], p_evap[0], 'ko', mfc='w')
            self.axes.text(h_evap[0], p_evap[0], '9$\quad\quad$', ha='right', va='top')

            h_XV                        = np.r_[Cycle.FlashTank.h_in / 1000., Cycle.Condenser.h_out_r / 1000.]
            p_XV                        = np.r_[Cycle.FlashTank.p_in / 1000, Cycle.Compressor.p_out_r / 1000]
            self.axes.plot(h_XV, p_XV, 'k', lw=2)

            h_FlashTank_liquid          = np.r_[Cycle.FlashTank.h_in / 1000., Cycle.FlashTank.h_out / 1000.]
            p_FlashTank_liquid          = np.r_[Cycle.FlashTank.p_in / 1000, Cycle.FlashTank.p_in / 1000]
            self.axes.semilogy(h_FlashTank_liquid, p_FlashTank_liquid, 'k', lw=2)
            self.axes.semilogy(h_FlashTank_liquid[0], p_FlashTank_liquid[0], 'ko', mfc='w')
            self.axes.text(h_FlashTank_liquid[0], p_FlashTank_liquid[0], '6$\quad\quad$', ha='right', va='top')

            h_1INJ                      = np.r_[Cycle.FlashTank.h_in / 1000, Cycle.Compressor.h_inj_r / 1000.]
            p_1INJ                      = np.r_[Cycle.FlashTank.p_in / 1000, Cycle.Compressor.p_inj_r / 1000.]
            self.axes.semilogy(h_1INJ, p_1INJ, 'k', lw=2)
            self.axes.semilogy(h_1INJ[0], p_1INJ[0], 'ro', mfc='w')
            self.axes.text(h_1INJ[0], p_1INJ[0], '8$\quad\quad$', ha='right', va='top')

            h_XV2                       = np.r_[Cycle.FlashTank.h_out / 1000., Cycle.Evaporator.h_in_r / 1000.]
            p_XV2                       = np.r_[Cycle.FlashTank.p_in / 1000., Cycle.Evaporator.p_sat_r / 1000.]
            self.axes.semilogy(h_XV[1], p_XV[1], 'ko', mfc='w')
            self.axes.semilogy(h_XV2, p_XV2, 'k', lw=2)
            self.axes.semilogy(h_XV2[0], p_XV2[0], 'ko', mfc='w')
            self.axes.text(h_XV2[0], p_XV2[0], '8$\quad\quad$', ha='right', va='top')

        plt.show()

    def PHEjector(self, Ejector, **kwargs):

        Ref                             = Ejector.Ref

        T_min                           = PropsSI('Tmin', Ref)
        T_max                           = PropsSI('Tmax', Ref)
        p_min                           = PropsSI('pmin', Ref)
        p_max                           = PropsSI('pmax', Ref)

        fig                             = plt.figure(figsize=(6, 4))
        self.axes                       = fig.add_subplot(111)

        T_sat                           = np.linspace(T_min, PropsSI(Ref, "Tcrit") - 0.0000000001, 1000)
        (h_sat_L, h_sat_V)              = (0.0 * T_sat, 0.0 * T_sat)
        (p_sat_L, p_sat_V)              = (0.0 * T_sat, 0.0 * T_sat)

        for i in np.arange(len(T_sat)):
            p_sat_L[i]                  = PropsSI('P', 'T', T_sat[i], 'Q', 0, Ref) / 1000
            p_sat_V[i]                  = PropsSI('P', 'T', T_sat[i], 'Q', 1, Ref) / 1000
            h_sat_L[i]                  = PropsSI('H', 'T', T_sat[i], 'Q', 0, Ref) / 1000
            h_sat_V[i]                  = PropsSI('H', 'T', T_sat[i], 'Q', 1, Ref) / 1000

        h_xd                            = np.zeros_like((T_sat))
        p_sat                           = np.linspace(p_min, 0.999 * PropsSI(Ref, "pcrit"), 1000)
        for i in np.arange(len(p_sat)):
            h_xd[i]                     = PropsSI('H', 'P', p_sat[i], 'Q', Ejector.x_d, Ref) / 1000

        # ph saturated line plot
        self.axes.semilogy(h_sat_L, p_sat_L,    'k', lw=1.5, alpha=0.5)
        self.axes.semilogy(h_sat_V, p_sat_V,    'k', lw=1.5, alpha=0.5)
        self.axes.semilogy(h_xd, p_sat / 1000,  'k', lw=1.5, alpha=0.3)

        self.axes.text(Ejector.h_d / 1000, 0.8 * Ejector.p_sat_ev / 1000, '%0.2f$\quad\quad$'
                       % Ejector.x_d, ha='right', va='top', backgroundcolor='w', color='k', alpha=0.3)

        # Axes labels
        self.axes.set_xlabel('h [kJ/kg]')
        self.axes.set_ylabel('p [kPa]')

        # Axes limits
        self.axes.set_xlim(0.7 * Ejector.h_out_gc / 1000, 1.1 * Ejector.h_out_ev / 1000)
        self.axes.set_ylim(0.7 * Ejector.p_sat_ev / 1000, 1.2 * Ejector.p_out_gc / 1000)

        # Gas cooler outlet - ejector inlet
        h_gc_mb                     = np.r_[Ejector.h_out_gc / 1000, Ejector.h_mb / 1000]
        p_gc_mb                     = np.r_[Ejector.p_out_gc / 1000, Ejector.p_mb / 1000]
        self.axes.semilogy(h_gc_mb, p_gc_mb, 'ko-', mfc='w')
        self.axes.text(h_gc_mb[0], p_gc_mb[0], 'gc$\quad\quad$', ha='left', va='bottom')
        self.axes.text(h_gc_mb[1], p_gc_mb[1], 'mb$\quad\quad$', ha='right', va='top')

        # Evaporator outlet
        h_out_ev                    = Ejector.h_out_ev / 1000
        p_sat_ev                    = Ejector.p_sat_ev / 1000
        self.axes.semilogy(h_out_ev, p_sat_ev, 'ko', mfc='w')
        self.axes.text(h_out_ev, p_sat_ev, 'ev$\quad\quad$', ha='left', va='bottom')

        # Ejector suction stream to mix
        h_sb_mix                    = np.r_[Ejector.h_sb / 1000, Ejector.h_mix / 1000]
        p_sb_mix                    = np.r_[Ejector.p_sb / 1000, Ejector.p_mix / 1000]
        self.axes.semilogy(h_sb_mix, p_sb_mix, 'ko-', mfc='w')
        self.axes.text(h_sb_mix[0], p_sb_mix[0], 'sb$\quad\quad$', ha='right', va='top')
        self.axes.text(h_sb_mix[1], p_sb_mix[1], 'mix$\quad\quad$', ha='right', va='top')

        # Ejector motion stream to mix
        h_mb_mix                    = np.r_[Ejector.h_mb / 1000, Ejector.h_mix / 1000]
        p_mb_mix                    = np.r_[Ejector.p_mb / 1000, Ejector.p_mix / 1000]
        self.axes.semilogy(h_mb_mix, p_mb_mix, 'ko-', mfc='w')

        # Ejector discharge
        h_d                         = Ejector.h_d / 1000
        p_d                         = Ejector.p_d / 1000
        self.axes.semilogy(h_d, p_d, 'or', mfc='w')
        self.axes.text(h_d, p_d, 'd$\quad\quad$', ha='left', va='bottom')

        self.axes.yaxis.set_major_formatter(ScalarFormatter())
        self.axes.yaxis.set_minor_formatter(ScalarFormatter())

        fig.savefig('Ejector.png', dpi=600)
        plt.show()
