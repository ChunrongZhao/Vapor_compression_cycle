from __future__         import division, print_function, absolute_import
import sys
from scipy.optimize import brentq, fsolve,newton
import numpy            as np                  # NumPy is fundamental scientific package
import CoolProp         as CP

from ACHP_codes.ACHP_Tools.convert_units                        import W2BTUh
from ACHP_codes.Components.Compressor                           import CompressorClass                      # Compressor
from ACHP_codes.Components.Condenser                            import CondenserClass                       # Condenser
from ACHP_codes.Components.MicroChannelCondenser                import MicroCondenserClass                  # MicroChannelCondenser
from ACHP_codes.Components.Evaporator                           import EvaporatorClass                      # Evaporator
from ACHP_codes.Components.CoolingCoil                          import CoolingCoilClass                     # Cooling Coil
from ACHP_codes.Components.MultiCircuitEvaporator               import MultiCircuitEvaporatorClass
from ACHP_codes.Components.CoaxialHX                            import CoaxialHXClass                       # Coaxial internal heat exchanger
from ACHP_codes.Components.PHE_Evaporator                       import PHEHXClass                           # Plate-Heat-Exchanger
from ACHP_codes.Components.LineSet                              import LineSetClass, LineSetOptionClass     # Line set class
from ACHP_codes.Components.ExpDev                               import ExpDevClass
from ACHP_codes.Components.Pump                                 import PumpClass                            # Secondary loop pump class
from ACHP_codes.Components.WavyChannelCooler                    import WavyChannelCoolerClass
from ACHP_codes.Components.GasCooler                            import GasCoolerClass                       # GasCooler
from ACHP_codes.Components.MicroChannelGasCooler                import MicroChannelGasCoolerClass           # Micro-channel GasCooler
from ACHP_codes.Components.MicroChannelEvaporator               import MicroChannelEvaporatorClass          # Micro-channel Evaporator
from ACHP_codes.Components.System_Volumes                       import SightGlassFilterDrierMicroMotionClass, SightGlassFilterDrierMicroMotionOptionClass
from ACHP_codes.Components.Compressor_VariableSpeed_HitachiRotary import HitachiVariableSpeedRotaryCompressorClass
from ACHP_codes.Components.Receiver                             import SuctionAccumulatorClass

from ACHP_codes.Correlations.PhaseState_Correlations            import TrhoPhase_ph
from ACHP_codes.ACHP_Tools.Solvers                              import MultiDimNewtRaph, Broyden, DE
from ACHP_codes.Cycles.Preconditioners                          import DXPreconditioner, SecondaryLoopPreconditioner, BTMS_SecondaryLoopPreconditioner
from ACHP_codes.Correlations.FinStructure_Correlations          import FinInputs                            # fin correlations
from ACHP_codes.Correlations.MicroFinCorrelations               import MicroFinInputs                       # Micro-fin correlations


# ----------------------------------------------------------------------
class DXCycleClass():
    def __init__(self):
        """
        Load up the necessary sub-structures to be filled with
        the code that follows
        """
        self.Compressor                 = CompressorClass()
        self.LineSetSuction             = LineSetOptionClass()
        self.LineSetDischarge           = LineSetOptionClass()
        self.LineSetLiquid              = LineSetOptionClass()
        self.ExpDev                     = ExpDevClass()

    def Update(self):
        """
        Update cycle class with selected HX type
        Update cycle class with Abstract State
        """
        if self.EvapSolver == 'Moving-Boundary':
            if self.EvapType == 'Fin-tube':
                self.Evaporator         = EvaporatorClass()
                self.Evaporator.Fins    = FinInputs()
            elif self.EvapType == 'Micro-channel':
                self.Evaporator         = MicroChannelEvaporatorClass()
                self.Evaporator.Fins    = MicroFinInputs()
            else:
                raise
        elif self.EvapSolver == 'Finite-Element':
            raise Exception('Discretized HX model not implemented')
        else:
            raise

        if self.CondSolver == 'Moving-Boundary':
            if self.CondType == 'Fin-tube':
                self.Condenser          = CondenserClass()
                self.Condenser.Fins     = FinInputs()
            elif self.CondType == 'Micro-channel':
                self.Condenser          = MicroCondenserClass()
                self.Condenser.Fins     = MicroFinInputs()
            else:
                raise
        elif self.CondSolver == 'Finite-Element':
            raise Exception('Discretized HX model not implemented')
        else:
            raise

        # Abstract State
        self.AS                         = CP.AbstractState(self.Backend, self.Ref)
        if hasattr(self, 'MassFrac'):
            self.AS.set_mass_fractions([self.MassFrac, 1 - self.MassFrac])
            self.AS.build_phase_envelope("dummy")
            self.AS.get_phase_envelope_data()
        elif hasattr(self, 'VoluFrac'):
            self.AS.set_volu_fractions([self.VoluFrac])
            self.AS.set_mass_fractions([self.VoluFrac, 1 - self.VoluFrac])
            self.AS.build_phase_envelope("dummy")
            self.AS.get_phase_envelope_data()

    def OutputList(self):
        """
            Return a list of parameters for this component for further output

            It is a list of tuples, and each tuple is formed of items:
                [0] Description of value
                [1] Units of value
                [2] The value itself
        """
        Output_List                     = []
        # append optional parameters, if applicable
        if hasattr(self, 'TestName'):
            Output_List.append(('Name',                     'N/A',  self.TestName))
        if hasattr(self, 'TestDescription'):
            Output_List.append(('Description',              'N/A',  self.TestDescription))
        if hasattr(self, 'TestDetails'):
            Output_List.append(('Details',                  'N/A',  self.TestDetails))
        if self.ImposedVariable == 'Charge':
            Output_List.append(('Charge w/no Correction',   'kg',   self.Charge_noCorrect))
            Output_List.append(('Charge w/ Correction',     'kg',   self.Charge_correct))
        Output_List_default = [                                     # default output list
            ('Charge',                      'kg',                   self.Charge),
            ('Condensation temp (dew)',     'K',                    self.T_dew_cond),
            ('Evaporation temp (dew)',      'K',                    self.T_dew_evap),
            ('Condenser Subcooling',        'K',                    self.DT_sc),
            ('Evaporator Superheat',        'K',                    self.DT_sh),
            ('Primary Ref.',                '-',                    self.Ref),
            ('COP',                         '-',                    self.COP),
            ('COSP',                        '-',                    self.COSP),
            ('Net Capacity',                'W',                    self.Capacity),
            ('Net Power',                   'W',                    self.Power),
            ('EER',                         'BTU/(W-hr)',           self.EER),
            ('SHR',                         '-',                    self.SHR),
            ('Imposed Variable',            '-',                    self.ImposedVariable),
         ]
        for i in range(0, len(Output_List_default)):                # append default parameters to output list
            Output_List.append(Output_List_default[i])
        return Output_List

    def Calculate(self, DT_evap, DT_cond, DT_sh):
        """
        Inputs are differences in temperature [K] between HX air inlet temperature
        and the dew temperature for the heat exchanger.

        Required Inputs:
            DT_evap:
                Difference in temperature [K] between evaporator air inlet temperature and refrigerant dew temperature
            DT_cond:
                Difference in temperature [K] between condenser air inlet temperature and refrigerant dew temperature
            DT_sh:
                Superheat value
        """
        if self.Verbosity > 1:
            print('DT_evap %7.4f DT_cond %7.4f DT_sh %7.4f,' %(DT_evap, DT_cond, DT_sh))

        # AbstractState
        AS                                  = self.AS

        # Surrounding
        if hasattr(self, 'T_0') and hasattr(self, 'p_0'):
            T_0                             = self.T_0
            p_0                             = self.p_0
        else:
            T_0                             = 25+273.15                 # K
            p_0                             = 101325                    # Pa
            self.T_0                        = T_0
            self.p_0                        = p_0

        AS.update(CP.PT_INPUTS, p_0, T_0)
        h_0                                 = AS.hmass()
        s_0                                 = AS.smass()

        self.h_0                            = h_0
        self.s_0                            = s_0

        # Condenser and evaporator dew temperature (guess)
        T_dew_cond                          = self.Condenser.Fins.Air.T_db + DT_cond    # the values (T_in_a,..) come from line 128ff
        T_dew_evap                          = self.Evaporator.Fins.Air.T_db - DT_evap
        # Condenser and evaporator saturation pressures
        AS.update(CP.QT_INPUTS, 1.0, T_dew_cond)
        p_sat_cond                          = AS.p()                    # [Pa]
        AS.update(CP.QT_INPUTS, 1.0, T_dew_evap)
        p_sat_evap                          = AS.p()                    # [Pa]
        # evaporator bubble temperature
        AS.update(CP.PQ_INPUTS, p_sat_evap, 0.0)
        T_bubble_evap                       = AS.T()                    # [T]

        self.T_dew_cond                     = T_dew_cond
        self.T_dew_evap                     = T_dew_evap

        # If the user doesn't include the Mode, fail
        assert hasattr(self, 'Mode')

        # --------------------------------------------------------------------------------
        # Cycle Solver in 'AC' model
        # --------------------------------------------------------------------------------
        if self.Mode == 'AC':
            # -------------------------------------------------------------------------------------
            if not hasattr(self.Compressor, 'm_dot_r') or self.Compressor.m_dot_r < 0.00001:
                # The first run of model, run the compressor just so you can get a preliminary value
                # for the mass flow rate for the line set
                params  = {                     # dictionary -> key:value, e.g. 'key':2345,
                    'p_in_r':               p_sat_evap,
                    'p_out_r':              p_sat_cond,
                    'T_in_r':               T_dew_evap + DT_sh,
                    'AS':                   AS,
                }
                self.Compressor.Update(**params)
                self.Compressor.Calculate()

            # -------------------------------------------------------------------------------------
            # Calculate inlet enthalpy
            AS.update(CP.PT_INPUTS, p_sat_evap, T_dew_evap + DT_sh)
            h_in                            = AS.hmass()                # [J/kg]
            params  = {
                'T_in':                     T_dew_evap + DT_sh,
                'p_in':                     p_sat_evap,
                'h_in':                     h_in,
                'm_dot':                    self.Compressor.m_dot_r,
                'AS':                       AS,
                'Name': '                   Suction Line'
            }
            self.LineSetSuction.Update(**params)
            self.LineSetSuction.Calculate()

            # -------------------------------------------------------------------------------------
            params  = {                     # dictionary -> key:value, e.g. 'key':2345,
                'p_in_r':                   p_sat_evap - self.DP_low,
                'p_out_r':                  p_sat_cond + self.DP_high,
                'T_in_r':                   TrhoPhase_ph(self.AS, p_sat_evap, self.LineSetSuction.h_out, T_bubble_evap, T_dew_evap)[0],
                'AS':                       AS,
            }
            self.Compressor.Update(**params)
            self.Compressor.Calculate()

            if self.Verbosity > 1:
                print('Comp DP L H', self.DP_low, self.DP_high)

            # -------------------------------------------------------------------------------------
            params  = {
                'T_in':                     self.Compressor.T_out_r,
                'p_in':                     p_sat_cond,
                'h_in':                     self.Compressor.h_out_r,
                'm_dot':                    self.Compressor.m_dot_r,
                'AS':                       AS,
                'Name':                     'Discharge Line'
            }
            self.LineSetDischarge.Update(**params)
            self.LineSetDischarge.Calculate()

            # -------------------------------------------------------------------------------------
            params  = {
                'm_dot_r':                  self.Compressor.m_dot_r,
                'T_in_r':                   self.LineSetDischarge.T_out,
                'p_sat_r':                  p_sat_cond,
                'AS':                       AS,
            }
            self.Condenser.Update(**params)
            self.Condenser.Calculate()

            # -------------------------------------------------------------------------------------
            params  = {
                'T_in':                     self.Condenser.T_out_r,
                'p_in':                     p_sat_cond,
                'h_in':                     self.Condenser.h_out_r,
                'm_dot':                    self.Compressor.m_dot_r,
                'AS':                       AS,
                'Name':                     'Liquid Line'
            }
            self.LineSetLiquid.Update(**params)
            self.LineSetLiquid.Calculate()

            # -------------------------------------------------------------------------------------
            params  = {
                'p_in_r':                   p_sat_cond,
                'h_in_r':                   self.LineSetLiquid.h_out,
                'T_sup':                    DT_sh,
                'p_out_r':                  p_sat_evap,
                'AS':                       AS,
            }
            self.ExpDev.Update(**params)
            self.ExpDev.Calculate()

            # -------------------------------------------------------------------------------------
            params  = {
                'm_dot_r':                  self.Compressor.m_dot_r,
                'p_sat_r':                  p_sat_evap,
                'h_in_r':                   self.ExpDev.h_out_r,
                'AS':                       AS,
            }
            self.Evaporator.Update(**params)
            self.Evaporator.Calculate()

            # -------------------------------------------------------------------------------------
            # Energy and Mass conservation
            self.EnergyBalance = self.Compressor.CycleEnergyIn + self.Condenser.Q + self.Evaporator.Q + self.LineSetSuction.Q \
                                 + self.LineSetDischarge.Q + self.LineSetLiquid.Q

            Charge              = self.Compressor.Charge + self.Condenser.Charge + self.Evaporator.Charge + self.LineSetSuction.Charge \
                                  + self.LineSetDischarge.Charge + self.LineSetLiquid.Charge

            # Charge with no correction
            self.Charge_noCorrect           = Charge

            # Correct the charge with Bo Shen's two point regression method
            if self.ImposedVariable == 'Charge':
                C                           = self.C                                # [kg]
                if self.ChargeMethod == 'One-point':
                    self.Charge_correct = Charge + C                                # one-point charge correction
                elif self.ChargeMethod == 'Two-point':
                    K                       = self.K                                # [kg]
                    w_ref                   = self.w_ref                            # [-]
                    w_liq                   = self.Condenser.w_subcool
                    self.delta_charge       = C + K * (w_liq - w_ref)
                    self.Charge_correct     = Charge + self.delta_charge            # two-point charge correction
                self.Charge                 = self.Charge_correct
            else:
                self.Charge                 = Charge

            # -------------------------------------------------------------------
            if self.ExpDev.ExpType == 'Ideal':
                resid                       = np.zeros((3))
                resid[0]                    = self.Compressor.m_dot_r * (self.Evaporator.h_out_r - self.LineSetSuction.h_in)
                resid[2]                    = DT_sh - self.Evaporator.DT_sh
            else:
                resid                       = np.zeros((3))
                resid[0]                    = self.Evaporator.h_out_r - self.LineSetSuction.h_in
                resid[2]                    = self.Compressor.m_dot_r - self.ExpDev.m_dot_r

            # -------------------------------------------------------------------
            if self.ImposedVariable == 'Subcooling':
                resid[1]                    = self.Condenser.DT_sc - self.DT_sc_target
            elif self.ImposedVariable == 'Charge':
                resid[1]                    = self.Charge - self.Charge_target

            if self.Verbosity > 1:
                print(resid)

            self.Capacity                   = self.Evaporator.Capacity
            self.Power                      = self.Compressor.W + self.Evaporator.Fins.Air.FanPower + self.Condenser.Fins.Air.FanPower
            self.COP                        = self.Evaporator.Q / self.Compressor.W
            self.COSP                       = self.Evaporator.Capacity / self.Power
            self.EER                        = W2BTUh(self.Capacity) / self.Power        # [BTU/(W-hr)]
            self.SHR                        = self.Evaporator.SHR
            self.DT_sc                      = self.Condenser.DT_sc
            self.DP_HighPressure            = self.Condenser.DP_r + self.LineSetLiquid.DP + self.LineSetDischarge.DP
            self.DP_LowPressure             = self.Evaporator.DP_r + self.LineSetSuction.DP

        # --------------------------------------------------------------------------------
        # Cycle Solver in 'HP' model
        # --------------------------------------------------------------------------------
        elif self.Mode == 'HP':
            if not hasattr(self.Compressor, 'm_dot_r') or self.Compressor.m_dot_r < 0.00001:
                # The first run of model, run the compressor just so you can get a preliminary value
                # for the mass flow rate for the line set
                # ------------------------------------------------------------------------
                params  = {                     # dictionary -> key:value, e.g. 'key':2345,
                    'p_in_r':               p_sat_evap,
                    'p_out_r':              p_sat_cond,
                    'T_in_r':               T_dew_evap + DT_sh,
                    'AS':                   AS,
                }
                self.Compressor.Update(**params)
                self.Compressor.Calculate()

            # -----------------------------------------------------------------------------
            # Calculate inlet enthalpy
            AS.update(CP.PT_INPUTS, p_sat_evap, T_dew_evap + DT_sh)
            h_in                            = AS.hmass()                    # [J/kg]
            params  = {
                'T_in':                     T_dew_evap + DT_sh,
                'p_in':                     p_sat_evap,
                'h_in':                     h_in,
                'm_dot':                    self.Compressor.m_dot_r,
                'AS':                       AS,
                'Name':                     'Suction Line'
            }
            self.LineSetSuction.Update(**params)
            self.LineSetSuction.Calculate()

            # -----------------------------------------------------------------------------
            params  = {                     # dictionary -> key:value, e.g. 'key':2345,
                'p_in_r':                   p_sat_evap - self.DP_low,
                'p_out_r':                  p_sat_cond + self.DP_high,
                'T_in_r':                   TrhoPhase_ph(self.AS, p_sat_evap, self.LineSetSuction.h_out, T_bubble_evap, T_dew_evap)[0],
                'Ref':                      self.Ref,
                'AS':                       AS
            }
            self.Compressor.Update(**params)
            self.Compressor.Calculate()

            # -----------------------------------------------------------------------------
            params  = {
                'T_in':                     self.Compressor.T_out_r,
                'p_in':                     p_sat_cond,
                'h_in':                     self.Compressor.h_out_r,
                'm_dot':                    self.Compressor.m_dot_r,
                'AS':                       AS,
                'Name':                     'Discharge Line'
            }
            self.LineSetDischarge.Update(**params)
            self.LineSetDischarge.Calculate()

            # -----------------------------------------------------------------------------
            params  = {
                'm_dot_r':                  self.Compressor.m_dot_r,
                'T_in_r':                   self.LineSetDischarge.T_out,
                'p_sat_r':                  p_sat_cond,
                'AS':                       AS,
            }
            self.Condenser.Update(**params)
            self.Condenser.Calculate()

            # -----------------------------------------------------------------------------
            params  = {
                'T_in':                     self.Condenser.T_out_r,
                'p_in':                     p_sat_cond,
                'h_in':                     self.Condenser.h_out_r,
                'm_dot':                    self.Compressor.m_dot_r,
                'AS':                       AS,
                'Name':                     'Liquid Line'
            }
            self.LineSetLiquid.Update(**params)
            self.LineSetLiquid.Calculate()

            # -----------------------------------------------------------------------------
            params  = {
                'p_in_r':                   p_sat_cond,
                'h_in_r':                   self.LineSetLiquid.h_out,
                'T_sup':                    DT_sh,
                'p_out_r':                  p_sat_evap,
                'AS':                       AS,
            }
            self.ExpDev.Update(**params)
            self.ExpDev.Calculate()

            # -----------------------------------------------------------------------------
            params  = {
                'm_dot_r':                  self.Compressor.m_dot_r,
                'p_sat_r':                  p_sat_evap,
                'h_in_r':                   self.ExpDev.h_out_r,
                'AS':                       AS,
            }
            self.Evaporator.Update(**params)
            self.Evaporator.Calculate()

            # -----------------------------------------------------------------------------
            # Energy and Mass conservation

            self.EnergyBalance  = self.Compressor.CycleEnergyIn + self.Condenser.Q + self.Evaporator.Q + self.LineSetSuction.Q \
                                  + self.LineSetDischarge.Q + self.LineSetLiquid.Q

            Charge              = self.Compressor.Charge + self.Condenser.Charge + self.Evaporator.Charge + self.LineSetSuction.Charge \
                                  + self.LineSetDischarge.Charge + self.LineSetLiquid.Charge

            self.Charge_noCorrect           = Charge                    # Charge with no correction
            # Correct the charge with Bo Shen's two point regression method
            if self.ImposedVariable == 'Charge':
                C                           = self.C                    # [kg]
                if self.ChargeMethod == 'One-point':
                    self.Charge_correct     = Charge + C                # one-point charge correction
                elif self.ChargeMethod == 'Two-point':
                    K                       = self.K                    # [kg]
                    w_ref                   = self.w_ref                # [-]
                    w_liq                   = self.Condenser.w_subcool
                    self.delta_charge       = C + K * (w_liq - w_ref)
                    self.Charge_correct     = Charge + self.delta_charge # two-point charge correction
                self.Charge                 = self.Charge_correct
            else:
                self.Charge                 = Charge

            # -----------------------------------------------------------------------------
            if self.ExpDev.ExpType == 'Ideal':
                resid                       = np.zeros((3))
                resid[0]                    = self.Compressor.m_dot_r * (self.Evaporator.h_out_r - self.LineSetSuction.h_in)
                resid[2]                    = DT_sh - self.Evaporator.DT_sh
            else:
                resid                       = np.zeros((3))
                resid[0]                    = self.Evaporator.h_out_r - self.LineSetSuction.h_in
                resid[2]                    = self.Compressor.m_dot_r - self.ExpDev.m_dot_r

            # -----------------------------------------------------------------------------
            if self.ImposedVariable == 'Subcooling':
                resid[1]                    = self.Condenser.DT_sc - self.DT_sc_target
            elif self.ImposedVariable == 'Charge':
                resid[1]                    = self.Charge - self.Charge_target

            # -----------------------------------------------------------------------------
            self.Capacity                   = -self.Condenser.Q + self.Condenser.Fins.Air.FanPower
            self.DT_sc                      = self.Condenser.DT_sc
            self.Power                      = self.Compressor.W + self.Evaporator.Fins.Air.FanPower + self.Condenser.Fins.Air.FanPower
            self.COP                        = -self.Condenser.Q / self.Compressor.W
            self.COSP                       = self.Capacity / self.Power
            self.EER                        = W2BTUh(self.Capacity) / self.Power            # [BTU/(W-hr)]
            self.SHR                        = self.Evaporator.SHR
            self.DP_HighPressure            = self.Condenser.DP_r + self.LineSetDischarge.DP + self.LineSetLiquid.DP
            self.DP_LowPressure             = self.Evaporator.DP_r + self.LineSetSuction.DP

        # --------------------------------------------------------------------------------
        # Other options
        # --------------------------------------------------------------------------------
        else:
            ValueError("DX Cycle mode must be 'AC', or 'HP'")

        if self.Verbosity > 1:
            print('DT_evap {:4.4f} DT_cond {:4.4f} DT_sh {:4.4f} resid[0] {:12.6e} resid[1]: {:12.6e} resid[3]: {:12.6e} Charge {:4.4f} SC: {:4.4f}'
                   .format(DT_evap, DT_cond, DT_sh, resid[0], resid[1], resid[2], self.Charge, self.Condenser.DT_sc))
        self.DT_evap                        = DT_evap
        self.DT_cond                        = DT_cond
        self.DT_sh                          = DT_sh

        return resid

    def PreconditionedSolve(self):
        """
        Solver that will precondition by trying a range of DeltaT until the model
        can solve, then will kick into 2-D Newton Raphson solve

        The two input variables for the system solver are the differences in
        temperature between the inlet air temperature of the heat exchanger and the
        dew temperature of the refrigerant.  This is important for refrigerant blends
        with temperature glide during constant-pressure evaporation or condensation.
        Good examples of common working fluid with glide would be R404A or R410A.
        """
        def OBJECTIVE_DXCycle(x):
            """
            A wrapper function to convert input vector for fsolve to the proper form for the solver
            """
            try:
                resids  = self.Calculate(DT_evap=float(x[0]), DT_cond=float(x[1]), DT_sh=float(x[2])) #,DP_low=float(x[2]),DP_high=float(x[3]))
            except ValueError:
                raise Exception('Something not correct!')
            return resids

        # Use the preconditioner to determine a reasonably good starting guess
        print ('-------------------------------------')
        print ('            Running Preconditioner   ')
        print ('-------------------------------------')
        DT_evap_init, DT_cond_init, DT_sh_init = DXPreconditioner(self)

        print ('-------------------------------------')
        print ('           Preconditioner Completed  ')
        print ('             Starting Main Cycle     ')
        print ('-------------------------------------')


        GoodRun                             = False
        while GoodRun == False:
            try:
                self.DP_low                 = 0
                self.DP_high                = 0
                DP_converged                = False
                while DP_converged == False:
                    # Actually run the Newton-Raphson solver to get the solution
                    x                       = Broyden(OBJECTIVE_DXCycle, [DT_evap_init, DT_cond_init, DT_sh_init])
                    delta_low               = abs(self.DP_low - abs(self.DP_LowPressure))
                    delta_high              = abs(self.DP_high - abs(self.DP_HighPressure))
                    self.DP_low             = abs(self.DP_LowPressure)
                    self.DP_high            = abs(self.DP_HighPressure)
                    # Update the guess values based on last converged values
                    DT_evap_init            = self.DT_evap
                    DT_cond_init            = self.DT_cond
                    DT_sh_init              = self.DT_sh
                    if delta_low < 1 and delta_high < 1:
                        DP_converged        = True
                    if self.Verbosity > 4:
                        print(self.DP_HighPressure, self.DP_LowPressure, 'DPHP')
                    GoodRun                 = True
            except AttributeError:
                # This will be a fatal error !! Should never have attribute error
                raise
            except:
                print ("--------------  Exception Caught ---------------- ")
                print ("Error of type", sys.exc_info()[0], " is: " + sys.exc_info()[1].message)
                raise

        if self.Verbosity > 0:

            print ("Debugging Results")
            print ()
            print ('Capacity: ', self.Capacity)
            print ('COP: ',self.COP)
            print ('COP (w/ both fans): ',self.COSP)
            print ('SHR: ',self.SHR)
            print ('UA_r_evap',self.Evaporator.UA_r)
            print ('UA_a_evap',self.Evaporator.UA_a)
            print ('UA_r_cond',self.Condenser.UA_r)
            print ('UA_a_cond',self.Condenser.UA_a)

        print ('-------------------------------------')
        print ('     Simulation Completed            ')
        print ('-------------------------------------')


class SecondaryCycleClass():
    def __init__(self):
        """
        Load up the necessary sub-structures to be filled with
        the code that follows
        """
        self.Compressor                 = CompressorClass()
        # Outdoor coil is a Condenser in cooling mode and evaporator in heating mode
        self.CoolingCoil                = CoolingCoilClass()
        self.CoolingCoil.Fins           = FinInputs()
        self.Pump                       = PumpClass()
        # Add both types of internal heat exchangers
        self.CoaxialIHX                 = CoaxialHXClass()
        self.PHEIHX                     = PHEHXClass()
        self.LineSetSupply              = LineSetOptionClass()
        self.LineSetReturn              = LineSetOptionClass()
        self.LineSetSuction             = LineSetOptionClass()
        self.LineSetDischarge           = LineSetOptionClass()
        self.LineSetLiquid              = LineSetOptionClass()

        # Make IHX an empty class for holding parameters common to PHE and Coaxial IHX
        class struct:
            pass
        self.IHX                        = struct()

    def Update(self):
        """
        Update cycle class with selected HX type
        Update cycle class with Abstract State
        """
        # -----------------------------------------------------------------------------
        if self.EvapSolver == 'Moving-Boundary':
            if self.EvapType == 'Fin-tube':
                self.Evaporator                 = EvaporatorClass()
                self.Evaporator.Fins            = FinInputs()
            elif self.EvapType == 'Micro-channel':
                self.Evaporator                 = MicroChannelEvaporatorClass()
                self.Evaporator.Fins            = MicroFinInputs()
            else:
                raise
        elif self.EvapSolver == 'Finite-Element':
            raise Exception('Discretized HX model not implemented')
        else:
            raise

        # -----------------------------------------------------------------------------
        if self.CondSolver == 'Moving-Boundary':
            if self.CondType == 'Fin-tube':
                self.Condenser                  = CondenserClass()
                self.Condenser.Fins             = FinInputs()
            elif self.CondType == 'Micro-channel':
                self.Condenser                  = MicroCondenserClass()
                self.Condenser.Fins             = MicroFinInputs()
            else:
                raise
        elif self.CondSolver == 'Finite-Element':
            raise Exception('Discretized HX model not implemented')
        else:
            raise

        # -----------------------------------------------------------------------------
        # Abstract State
        self.AS                     = CP.AbstractState(self.Backend, self.Ref)
        if hasattr(self, 'MassFrac'):
            self.AS.set_mass_fractions([self.MassFrac, 1 - self.MassFrac])
            self.AS.build_phase_envelope("dummy")
            self.AS.get_phase_envelope_data()
        elif hasattr(self, 'VoluFrac'):
            self.AS.set_volu_fractions([self.VoluFrac, 1 - self.VoluFrac])
            self.AS.build_phase_envelope("dummy")
            self.AS.get_phase_envelope_data()

        # Abstract State for SecLoopFluid
        self.AS_SLF                 = CP.AbstractState(self.Backend_SLF, self.SecLoopFluid)
        if hasattr(self, 'MassFrac_SLF'):
            self.AS_SLF.set_mass_fractions([self.MassFrac_SLF])
        elif hasattr(self, 'VoluFrac_SLF'):
            self.AS_SLF.set_volu_fractions([self.VoluFrac_SLF])

    def OutputList(self):
        """
            Return a list of parameters for this component for further output

            It is a list of tuples, and each tuple is formed of items:
                [0] Description of value
                [1] Units of value
                [2] The value itself
        """
        return [
            ('Charge',                  'kg',           self.Charge),
            ('Condenser Subcooling',    'K',            self.DT_sc),
            ('Primary Ref.',            '-',            self.Ref),
            ('Secondary Ref.',          '-',            self.SecLoopFluid),
            ('Imposed Variable',        '-',            self.ImposedVariable),
            ('IHX Type',                '-',            self.IHXType),
            ('COP',                     '-',            self.COP),
            ('COSP',                    '-',            self.COSP),
            ('Net Capacity',            'W',            self.CoolingCoil.Capacity),
            ('Net Power',               'W',            self.Power),
            ('EER',                     'BTU/(W-hr)',   self.EER),
            ('SHR',                     '-',            self.SHR),
            ('Condensation temp (dew)', 'K',            self.T_dew_cond),
            ('Evaporation temp (dew)',  'K',            self.T_dew_evap)]

    def Calculate(self, DT_evap, DT_cond, T_in_CC):
        """
        Inputs are differences in temperature [K] between HX air inlet temperature
        and the dew temperature for the heat exchanger.

        Required Inputs:
            DT_evap:
                Difference in temperature [K] between cooling coil air inlet temperature and refrigerant dew temperature
            DT_cond:
                Difference in temperature [K] between condenser air inlet temperature and refrigerant dew temperature
            Tin_CC:
                Inlet "glycol" temperature to line set feeding cooling coil
        """
        if self.Verbosity > 1:
            print('Inputs: DT_evap %7.4f DT_cond %7.4f fT_IHX %7.4f' %(DT_evap, DT_cond, T_in_CC))

        # AbstractState
        AS                          = self.AS
        # AbstractState for SecLoopFluid
        AS_SLF                      = self.AS_SLF

        # Surrounding
        if hasattr(self, 'T_0') and hasattr(self, 'p_0'):
            T_0                     = self.T_0
            p_0                     = self.p_0
        else:
            T_0                     = 25+273.15                     # K
            p_0                     = 101325                        # Pa

        AS.update(CP.PT_INPUTS, p_0, T_0)
        h_0                         = AS.hmass()
        s_0                         = AS.smass()

        # Store the values to save on computation for later
        self.DT_evap                = DT_evap
        self.DT_cond                = DT_cond
        self.T_in_CC                = T_in_CC

        # If the user doesn't include the Mode, set it to Air Conditioning
        if not hasattr(self, 'Mode'):
            self.Mode               = 'AC'

        if self.Mode == 'AC':
            self.T_dew_cond         = self.Condenser.Fins.Air.T_db + DT_cond
            self.T_dew_evap         = self.CoolingCoil.Fins.Air.T_db - DT_evap
        elif self.Mode == 'HP':
            self.T_dew_cond         = T_in_CC + DT_cond
            self.T_dew_evap         = self.Evaporator.Fins.Air.T_db - DT_evap
        else:
            raise ValueError('Mode must be AC or HP')

        # -----------------------------------------------------------------------------
        # Evaporator and condenser saturation pressures
        AS.update(CP.QT_INPUTS, 1.0, self.T_dew_cond)
        p_sat_cond                  = AS.p()                        # [Pa]
        AS.update(CP.QT_INPUTS, 1.0, self.T_dew_evap)
        p_sat_evap                  = AS.p()                        # [Pa]

        # Evaporator and condenser bubble temperatures
        AS.update(CP.PQ_INPUTS, p_sat_evap, 0.0)
        self.T_bubble_evap          = AS.T()                        # [K]
        AS.update(CP.PQ_INPUTS, p_sat_cond, 0.0)
        self.T_bubble_cond          = AS.T()                        # [K]

        # ----------------------------------------------------------------------------------
        # Cycle solver for 'AC' mode
        # ----------------------------------------------------------------------------------
        if self.Mode == 'AC':
            if not hasattr(self.Compressor, 'm_dot_r') or self.Compressor.m_dot_r < 0.00001:
                # The first run of model, run the compressor just so you can get a preliminary value
                # for the mass flow rate for the line set
                # -----------------------------------------------------------------------------
                params  = {                 # dictionary -> key:value, e.g. 'key':2345,
                    'p_in_r':               p_sat_evap,
                    'p_out_r':              p_sat_cond,
                    'T_in_r':               self.T_dew_evap + self.DT_sh,
                    'AS':                   AS,
                }
                self.Compressor.Update(**params)
                self.Compressor.Calculate()

            # -----------------------------------------------------------------------------
            # Calculate inlet enthalpy
            AS.update(CP.PT_INPUTS, p_sat_evap, self.T_dew_evap + self.DT_sh)
            h_in                            = AS.hmass()                # [J/kg]
            params  = {
                'T_in':                     self.T_dew_evap + self.DT_sh,
                'p_in':                     p_sat_evap,
                'h_in':                     h_in,
                'm_dot':                    self.Compressor.m_dot_r,
                'AS':                       AS,
                'Name':                     'Suction Line'
            }
            self.LineSetSuction.Update(**params)
            self.LineSetSuction.Calculate()

            # -----------------------------------------------------------------------------
            params  = {                     # dictionary -> key:value, e.g. 'key':2345,
                'p_in_r':                   p_sat_evap + self.DP_low,
                'p_out_r':                  p_sat_cond - self.DP_high,
                'T_in_r':                   TrhoPhase_ph(self.AS, p_sat_evap, self.LineSetSuction.h_out, self.T_bubble_evap, self.T_dew_evap)[0],
                'AS':                       AS
            }
            self.Compressor.Update(**params)
            self.Compressor.Calculate()

            # -----------------------------------------------------------------------------
            params  = {
                'T_in':                     self.Compressor.T_out_r,
                'p_in':                     p_sat_cond,
                'h_in':                     self.Compressor.h_out_r,
                'm_dot':                    self.Compressor.m_dot_r,
                'AS':                       AS,
                'Name':                     'Discharge Line'
            }
            self.LineSetDischarge.Update(**params)
            self.LineSetDischarge.Calculate()

            # -----------------------------------------------------------------------------
            params  = {
                'm_dot_r':                  self.Compressor.m_dot_r,
                'T_in_r':                   self.Compressor.T_out_r,
                'p_sat_r':                  p_sat_cond,
                'AS':                       AS,
                # 'Oil': self.Oil,
                # 'shell_pressure':self.shell_pressure,
            }
            self.Condenser.Update(**params)
            self.Condenser.Calculate()

            # -----------------------------------------------------------------------------
            params  = {
                'T_in':                     self.Condenser.T_out_r,
                'p_in':                     p_sat_cond,
                'h_in':                     self.Condenser.h_out_r,
                'm_dot':                    self.Compressor.m_dot_r,
                'AS':                       AS,
                'Name':                     'Liquid Line'
            }
            self.LineSetLiquid.Update(**params)
            self.LineSetLiquid.Calculate()

            # -----------------------------------------------------------------------------
            # Inlet enthalpy to LineSetSupply
            AS_SLF.update(CP.PT_INPUTS, self.Pump.p_in_g, T_in_CC)
            h_in_LineSetSupply              = AS_SLF.hmass()            # [J/kg]
            params  = {
                'T_in':                     T_in_CC,
                'm_dot':                    self.Pump.m_dot_g,
                'h_in':                     h_in_LineSetSupply,
                'p_in':                     self.Pump.p_in_g,
                'AS':                       AS_SLF,
                'Name':                     'Supply Line'
            }
            self.LineSetSupply.Update(**params)
            self.LineSetSupply.Calculate()

            # -----------------------------------------------------------------------------
            # Now run CoolingCoil to predict inlet glycol temperature to IHX
            params  = {
                'm_dot_g':                  self.Pump.m_dot_g,
                'T_in_g':                   self.LineSetSupply.T_out,
                'p_in_g':                   self.CoolingCoil.p_in_g,
                'AS_g':                     AS_SLF
            }
            self.CoolingCoil.Update(**params)
            self.CoolingCoil.Calculate()

            # -----------------------------------------------------------------------------
            # Inlet enthalpy to LineSetReturn
            AS_SLF.update(CP.PT_INPUTS, self.Pump.p_in_g, self.CoolingCoil.T_out_g)
            h_in_LineSetReturn              = AS_SLF.hmass()            # [J/kg]
            params  = {
                'T_in':                     self.CoolingCoil.T_out_g,
                'm_dot':                    self.Pump.m_dot_g,
                'h_in':                     h_in_LineSetReturn,
                'p_in':                     self.CoolingCoil.p_in_g,
                'AS':                       AS_SLF,
                'Name':                     'Return Line'
            }
            self.LineSetReturn.Update(**params)
            self.LineSetReturn.Calculate()

            # -----------------------------------------------------------------------------
            if self.IHXType == 'Coaxial':
                params  = {
                    'm_dot_g':              self.Pump.m_dot_g,
                    'T_in_g':               self.CoolingCoil.T_out_g,
                    'p_in_g':               self.Pump.p_in_g,
                    'AS_g':                 AS_SLF,
                    'p_in_r':               p_sat_evap,
                    'h_in_r':               self.LineSetLiquid.h_out,
                    'AS_r':                 AS,
                    'm_dot_r':              self.Compressor.m_dot_r,
                }
                self.CoaxialIHX.Update(**params)
                self.CoaxialIHX.Calculate()

                self.IHX.Charge_r           = self.CoaxialIHX.Charge_r
                self.IHX.Q                  = self.CoaxialIHX.Q
                self.IHX.T_out_g            = self.CoaxialIHX.T_out_g
                self.IHX.DP_g               = self.CoaxialIHX.DP_g
                self.IHX.h_out_r            = self.CoaxialIHX.h_out_r
                self.IHX.T_out_r            = self.CoaxialIHX.T_out_r
                self.IHX.DP_r               = self.CoaxialIHX.DP_r
                if hasattr(self, 'PHEIHX'):
                    del self.PHEIHX

            # -----------------------------------------------------------------------------
            elif self.IHXType == 'PHE':
                # Inlet enthalpy to PHEIHX
                AS_SLF.update(CP.PT_INPUTS, self.PHEIHX.p_in_h, self.CoolingCoil.T_out_g)
                h_in_PHEIHX                 = AS_SLF.hmass()            # [J/kg]
                params  = {
                    'm_dot_h':              self.Pump.m_dot_g,
                    'h_in_h':               h_in_PHEIHX,
                    'p_in_h':               self.PHEIHX.p_in_h,
                    'AS_h':                 AS_SLF,
                    'm_dot_c':              self.Compressor.m_dot_r,
                    'p_in_c':               p_sat_evap,
                    'h_in_c':               self.LineSetLiquid.h_out,
                    'AS_c':                 AS
                }
                self.PHEIHX.Update(**params)
                self.PHEIHX.Calculate()

                self.IHX.Charge_r           = self.PHEIHX.Charge_c
                self.IHX.Q                  = self.PHEIHX.Q
                self.IHX.T_out_g            = self.PHEIHX.T_out_h
                self.IHX.DP_g               = self.PHEIHX.DP_h
                self.IHX.DP_r               = self.PHEIHX.DP_c
                self.IHX.h_out_r            = self.PHEIHX.h_out_c
                self.IHX.T_out_r            = self.PHEIHX.T_out_c
                if hasattr(self, 'CoaxialIHX'):
                    del self.CoaxialIHX

            # -----------------------------------------------------------------------------
            params  = {
                'DP_g':     self.IHX.DP_g + self.CoolingCoil.DP_g + self.LineSetSupply.DP + self.LineSetReturn.DP,
                'T_in_g':                   self.CoolingCoil.T_out_g,
                'p_in_g':                   self.Pump.p_in_g,
                'AS_g':                     AS_SLF,
            }
            self.Pump.Update(**params)
            self.Pump.Calculate()

            # -----------------------------------------------------------------------------
            # Energy and Mass conservation
            self.EnergyBalance      = self.Compressor.CycleEnergyIn + self.Condenser.Q + self.IHX.Q + self.LineSetSuction.Q \
                                      + self.LineSetDischarge.Q + self.LineSetLiquid.Q

            Charge                  = self.Compressor.Charge + self.Condenser.Charge + self.IHX.Charge_r + self.LineSetSuction.Charge \
                                      + self.LineSetDischarge.Charge + self.LineSetLiquid.Charge

            self.Charge_noCorrect           = Charge                            # Charge with no correction
            # Correct the charge with Bo Shen's two point regression method
            if self.ImposedVariable == 'Charge':
                C                           = self.C                            # [kg]
                if self.ChargeMethod == 'One-point':
                    self.Charge_correct     = Charge + C                        # one-point charge correction
                elif self.ChargeMethod == 'Two-point':
                    K                       = self.K                            # [kg]
                    w_ref                   = self.w_ref                        # [-]
                    w_liq                   = self.Condenser.w_subcool
                    self.delta_charge       = C + K * (w_liq - w_ref)
                    self.Charge_correct     = Charge + self.delta_charge        # two-point charge correction
                self.Charge                 = self.Charge_correct
            else:
                self.Charge                 = Charge

            # -----------------------------------------------------------------------------
            # Calculate properties
            AS.update(CP.QT_INPUTS, 0.0, self.T_bubble_cond)
            h_L                             = AS.hmass()                        # [J/kg]
            cp_L                            = AS.cpmass()                       # [J/kg-K]
            AS.update(CP.PT_INPUTS, p_sat_cond, self.T_bubble_cond - self.DT_sc_target)
            h_target                        = AS.hmass()                        # [J/kg]
            self.DT_sc                      = (h_L - self.Condenser.h_out_r) / cp_L
            delta_H_sc                      = self.Compressor.m_dot_r * (h_L - h_target)

            # -----------------------------------------------------------------------------
            resid                           = np.zeros((3))
            resid[0]                        = self.Compressor.m_dot_r * (self.IHX.h_out_r - self.LineSetSuction.h_in)

            if self.ImposedVariable == 'Subcooling':
                resid[1]                    = self.Condenser.DT_sc - self.DT_sc_target
            elif self.ImposedVariable == 'Charge':
                resid[1]                    = self.Charge - self.Charge_target
            # resid[2]=self.IHX.Q-self.CoolingCoil.Q+self.Pump.W

            self.residSL    = self.IHX.Q - self.CoolingCoil.Q + self.Pump.W + self.LineSetSupply.Q + self.LineSetReturn.Q
            resid[2]                        = self.residSL

            if self.Verbosity > 7:
                print('W_comp % 12.6e Q_cond: % 12.6e Q_PHE %10.4f ' %(self.Compressor.W, self.Condenser.Q, self.IHX.Q))
            if self.Verbosity > 1:
                print('Q_res % 12.6e Resid2: % 12.6e ResSL %10.4f  Charge %10.4f SC: %8.4f'
                      % (resid[0], resid[1], self.residSL, self.Charge, self.Condenser.DT_sc))

            # -----------------------------------------------------------------------------
            self.Capacity                   = self.CoolingCoil.Capacity
            self.COP                        = self.CoolingCoil.Q / self.Compressor.W
            self.COSP       = self.CoolingCoil.Capacity / (self.Compressor.W + self.Pump.W + self.CoolingCoil.Fins.Air.FanPower
                                                           + self.Condenser.Fins.Air.FanPower)
            self.SHR                        = self.CoolingCoil.SHR
            self.Power      = self.Compressor.W + self.Pump.W + self.CoolingCoil.Fins.Air.FanPower + self.Condenser.Fins.Air.FanPower
            self.EER                        = W2BTUh(self.Capacity) / self.Power                                        # [BTU/(W-hr)]
            self.DP_high_Model              = self.Condenser.DP_r + self.LineSetDischarge.DP + self.LineSetLiquid.DP    # [Pa]
            self.DP_low_Model               = self.IHX.DP_r + self.LineSetSuction.DP                                    # [Pa]

        # ----------------------------------------------------------------------------------
        # Cycle solver for 'HP' mode
        # ----------------------------------------------------------------------------------
        elif self.Mode == 'HP':
            if p_sat_evap + self.DP_low < 0:
                raise ValueError('Compressor inlet pressure less than zero ['+str(p_sat_evap + self.DP_low) +
                                 ' Pa] - is low side pressure drop too high?')

            if not hasattr(self.Compressor, 'm_dot_r') or self.Compressor.m_dot_r < 0.00001:
                # The first run of model, run the compressor just so you can get a preliminary value
                # for the mass flow rate for the line set
                # ----------------------------------------------------------------------------------
                params  = {                     # dictionary -> key:value, e.g. 'key':2345,
                    'p_in_r':                   p_sat_evap,
                    'p_out_r':                  p_sat_cond,
                    'T_in_r':                   self.T_dew_evap+self.DT_sh,
                    'AS':                       AS,
                }
                self.Compressor.Update(**params)
                self.Compressor.Calculate()

            # ----------------------------------------------------------------------------------
            # Calculate inlet enthalpy
            AS.update(CP.PT_INPUTS, p_sat_evap, self.T_dew_evap + self.DT_sh)
            h_in                                = AS.hmass()                        # [J/kg]
            params  = {
                'T_in':                         self.T_dew_evap+self.DT_sh,
                'p_in':                         p_sat_evap,
                'h_in':                         h_in,
                'm_dot':                        self.Compressor.m_dot_r,
                'AS':                           AS,
                'Name':                         'Suction Line'
            }
            self.LineSetSuction.Update(**params)
            self.LineSetSuction.Calculate()

            # ----------------------------------------------------------------------------------
            params  = {                         # dictionary -> key:value, e.g. 'key':2345,
                'p_in_r':                       p_sat_evap + self.DP_low,
                'p_out_r':                      p_sat_cond - self.DP_high,
                'T_in_r': TrhoPhase_ph(self.AS, p_sat_evap, self.LineSetSuction.h_out, self.T_bubble_evap, self.T_dew_evap)[0],
                'AS':                           AS
            }
            self.Compressor.Update(**params)
            self.Compressor.Calculate()

            # ----------------------------------------------------------------------------------
            params  = {
                'T_in':                         self.Compressor.T_out_r,
                'p_in':                         p_sat_cond,
                'h_in':                         self.Compressor.h_out_r,
                'm_dot':                        self.Compressor.m_dot_r,
                'AS':                           AS,
                'Name':                         'Discharge Line'
            }
            self.LineSetDischarge.Update(**params)
            self.LineSetDischarge.Calculate()

            # ----------------------------------------------------------------------------------
            # Inlet enthalpy to LineSetSupply
            AS_SLF.update(CP.PT_INPUTS, self.Pump.p_in_g, T_in_CC)
            h_in_LineSetSupply                  = AS_SLF.hmass()                # [J/kg]
            params  = {
                'T_in':                         T_in_CC,
                'm_dot':                        self.Pump.m_dot_g,
                'h_in':                         h_in_LineSetSupply,
                'AS':                           AS_SLF,
                'Name':                         'Supply Line'
            }
            self.LineSetSupply.Update(**params)
            self.LineSetSupply.Calculate()

            # ----------------------------------------------------------------------------------
            # Now run CoolingCoil to predict inlet glycol temperature to IHX
            params  = {
                'm_dot_g':                      self.Pump.m_dot_g,
                'T_in_g':                       self.LineSetSupply.T_out,
                'p_in_g':                       self.PHEIHX.p_in_c,
                'AS_g':                         AS_SLF
            }
            self.CoolingCoil.Update(**params)
            self.CoolingCoil.Calculate()

            # ----------------------------------------------------------------------------------
            # Inlet enthalpy to LineSetReturn
            AS_SLF.update(CP.PT_INPUTS, self.Pump.p_in_g, self.CoolingCoil.T_out_g)
            h_in_LineSetReturn                  = AS_SLF.hmass()                # [J/kg]
            params  = {
                'T_in':                         self.CoolingCoil.T_out_g,
                'm_dot':                        self.Pump.m_dot_g,
                'h_in':                         h_in_LineSetReturn,
                'AS':                           AS_SLF,
                'Name':                         'Return Line'
            }
            self.LineSetReturn.Update(**params)
            self.LineSetReturn.Calculate()

            # ----------------------------------------------------------------------------------
            # Inlet enthalpy to PHEIHX
            AS_SLF.update(CP.PT_INPUTS, self.PHEIHX.p_in_c, self.LineSetReturn.T_out)
            h_in_PHEIHX                         = AS_SLF.hmass()                # [J/kg]
            params  = {
                'm_dot_h':                      self.Compressor.m_dot_r,
                'h_in_h':                       self.LineSetDischarge.h_out,
                'p_in_h':                       p_sat_cond,
                'AS_h':                         AS,
                'm_dot_c':                      self.Pump.m_dot_g,
                'h_in_c':                       h_in_PHEIHX,
                'p_in_c':                       self.Pump.p_in_g,
                'AS_c':                         AS_SLF,

            }
            self.PHEIHX.Update(**params)
            self.PHEIHX.Calculate()

            # ----------------------------------------------------------------------------------
            params  = {
                'T_in':                         self.PHEIHX.T_out_h,
                'p_in':                         p_sat_cond,
                'h_in':                         self.PHEIHX.h_out_h,
                'm_dot':                        self.Compressor.m_dot_r,
                'AS':                           AS,
                'Name':                         'Liquid Line'
            }
            self.LineSetLiquid.Update(**params)
            self.LineSetLiquid.Calculate()

            # ----------------------------------------------------------------------------------
            params  = {
                'm_dot_r':                      self.Compressor.m_dot_r,
                'p_sat_r':                      p_sat_evap,
                'h_in_r':                       self.LineSetLiquid.h_out,
                'AS':                           AS,
            }
            self.Evaporator.Update(**params)
            self.Evaporator.Calculate()

            # ----------------------------------------------------------------------------------
            params  = {
                'DP_g':      self.PHEIHX.DP_c + self.CoolingCoil.DP_g + self.LineSetSupply.DP + self.LineSetReturn.DP,
                'T_in_g':                       self.CoolingCoil.T_out_g,
                'p_in_g':                       self.Pump.p_in_g,
                'AS_g':                         AS_SLF
            }
            self.Pump.Update(**params)
            self.Pump.Calculate()

            # ----------------------------------------------------------------------------------
            # No Energy conservation???????
            Charge          = self.Compressor.Charge + self.Evaporator.Charge + self.PHEIHX.Charge_h + self.LineSetSuction.Charge \
                              + self.LineSetDischarge.Charge + self.LineSetLiquid.Charge

            self.Charge_noCorrect               = Charge                            # Charge with no correction

            # Correct the charge with Bo Shen's two point regression method
            if self.ImposedVariable == 'Charge':
                C                               = self.C                            # [kg]
                if self.ChargeMethod == 'One-point':
                    self.Charge_correct         = Charge + C                        # one-point charge correction
                elif self.ChargeMethod == 'Two-point':
                    K                           = self.K                            # [kg]
                    w_ref                       = self.w_ref                        # [-]
                    w_liq                       = self.Condenser.w_subcool
                    self.delta_charge           = C + K * (w_liq - w_ref)
                    self.Charge_correct         = Charge + self.delta_charge        # two-point charge correction
                self.Charge                     = self.Charge_correct
            else:
                self.Charge                     = Charge

            # ----------------------------------------------------------------------------------
            # Calculate properties
            AS.update(CP.QT_INPUTS, 0.0, self.T_bubble_cond)
            h_L                                 = AS.hmass()                        # [J/kg]
            cp_L                                = AS.cpmass()                       # [J/kg-K]
            AS.update(CP.PT_INPUTS, p_sat_cond, self.T_bubble_cond - self.DT_sc_target)
            h_target                            = AS.hmass()                        # [J/kg]
            # Calculate an effective subcooling amount by delta_h/cp_satL
            # Can be positive or negative (negative is quality at outlet
            # (PropsSI('H','T',self.Tbubble_cond,'Q',0,self.Ref)-self.PHEIHX.hout_h)/
            # (PropsSI('C','T',self.Tbubble_cond,'Q',0,self.Ref)) #*1000 #*1000
            self.DT_sc                          = self.PHEIHX.DT_sc_h
            delta_H_sc                          = self.Compressor.m_dot_r * (h_L - h_target)

            # ----------------------------------------------------------------------------------
            resid                   = np.zeros((3))
            resid[0]                = self.Compressor.m_dot_r * (self.Evaporator.h_out_r - self.LineSetSuction.h_in)

            if self.ImposedVariable == 'Subcooling':
                resid[1]            = self.DT_sc - self.DT_sc_target
            elif self.ImposedVariable == 'Charge':
                resid[1]            = self.Charge - self.Charge_target
            self.residSL            = self.PHEIHX.Q + self.CoolingCoil.Q + self.Pump.W + self.LineSetSupply.Q + self.LineSetReturn.Q
            resid[2]                = self.residSL

            if self.Verbosity > 1:
                print('Q_res % 12.6e Resid2: % 12.6e ResSL %10.4f Charge %10.4f SC: %8.4f'
                      %(resid[0], resid[1], self.residSL, self.Charge, self.DT_sc))

            # ----------------------------------------------------------------------------------
            self.Capacity           = -self.CoolingCoil.Q + self.CoolingCoil.Fins.Air.FanPower
            self.COP                = -self.CoolingCoil.Q / self.Compressor.W
            self.COSP   = self.Capacity / (self.Compressor.W + self.Pump.W + self.CoolingCoil.Fins.Air.FanPower + self.Evaporator.Fins.Air.FanPower)
            self.Power  = self.Compressor.W + self.Pump.W + self.CoolingCoil.Fins.Air.FanPower + self.Evaporator.Fins.Air.FanPower
            self.EER                = W2BTUh(self.Capacity) / self.Power                                    # [BTU/(W-hr)]
            self.SHR                = -self.CoolingCoil.SHR
            self.DP_high_Model      = self.PHEIHX.DP_h + self.LineSetDischarge.DP + self.LineSetLiquid.DP   # [Pa]
            self.DP_low_Model       = self.Evaporator.DP_r + self.LineSetSuction.DP                         # [Pa]

        self.DT_evap                = DT_evap
        self.DT_cond                = DT_cond

        return resid

    def PreconditionedSolve(self, PrecondValues=None):
        """
        PrecondValues = dictionary of values DT_evap, DT_cond and Tin_CC
        """

        def OBJECTIVE(x):
            """
            Takes the place of a lambda function since lambda functions do not bubble error properly
            """
            return self.Calculate(x[0], x[1], x[2])

        def OBJECTIVE2(x, T_in):
            """
            Takes the place of a lambda function since lambda functions do not bubble error properly
            """
            return self.Calculate(x[0], x[1], T_in)

        def OBJECTIVE_SL(T_in_CC):
            """
            Objective function for the inner loop of the vapor compression system

            Using the MultiDimNewtRaph function will re-evaluate the Jacobian at
            every step.  Slower, but more robust since the solution surfaces aren't
            smooth enough

            Note: This function is not currently used!
            """
            x   = MultiDimNewtRaph(OBJECTIVE2, [self.DT_evap, self.DT_cond], args=(T_in_CC,))

            # Update the guess values for Delta Ts starting
            # at the third step (after at least one update away
            # from the boundaries)
            if self.OBJ_SL_counter >= 0:
                self.DT_evap                    = x[0]
                self.DT_cond                    = x[1]
                pass
            self.OBJ_SL_counter                 += 1
            return self.residSL

        def PrintDPs():
            print('DP_LP :: Input:', self.DP_low,  'Pa / Model calc:', self.DP_low_Model,  'Pa')
            print('DP_HP :: Input:', self.DP_high, 'Pa / Model calc:', self.DP_high_Model, 'Pa')

        # Some variables need to be initialized
        self.DP_low                             = 0     # The actual low-side pressure drop to be used in Pa
        self.DP_high                            = 0     # The actual low-side pressure drop to be used in Pa
        self.OBJ_SL_counter                     = 0

        # Run the preconditioner to get guess values for the temperatures
        if PrecondValues is None:
            self.DT_evap, self.DT_cond, T_in_CC = SecondaryLoopPreconditioner(self)
        else:
            self.DT_evap                        = PrecondValues['DT_evap']
            self.DT_cond                        = PrecondValues['DT_cond']
            T_in_CC                             = PrecondValues['T_in_CC']

        # Remove the other, non-used IHX class if found
        if self.IHXType == 'PHE':
            if hasattr(self, 'CoaxialIHX'):
                del self.CoaxialIHX
        else:
            if hasattr(self, 'PHEIHX'):
                del self.PHEIHX

        # Remove the condenser if in heating mode and condenser found
        if self.Mode == 'HP':
            if hasattr(self, 'Condenser'):
                del self.Condenser

        # ----------------------------------------------------------------------------------
        iter                    = 1
        max_error_DP            = 999
        # Outer loop with a more relaxed convergence criterion
        while max_error_DP > 0.5:
            iter_inner          = 1
            # Inner loop to determine pressure drop for high and low sides
            while max_error_DP > 0.05 and iter_inner < 10:

                # Run to calculate the pressure drop as starting point
                OBJECTIVE([self.DT_evap, self.DT_cond, T_in_CC])

                # Calculate the max error
                max_error_DP    = max([abs(self.DP_low_Model - self.DP_low), abs(self.DP_high_Model - self.DP_high)])

                if self.Verbosity > 0:
                    PrintDPs()
                    print('Max pressure drop error [inner loop] is', max_error_DP, 'Pa')

                # Update the pressure drop terms
                self.DP_low     = self.DP_low_Model         # /1000
                self.DP_high    = self.DP_high_Model        # /1000

                iter_inner      += 1

            if self.Verbosity > 0:
                print("Done with the inner loop on pressure drop")

            # Use Newton-Raphson solver
            (self.DT_evap, self.DT_cond, T_in_CC)   = MultiDimNewtRaph(OBJECTIVE, [self.DT_evap, self.DT_cond, T_in_CC], dx=0.1)

            # Calculate the error
            max_error_DP    = max([abs(self.DP_low_Model - self.DP_low), abs(self.DP_high_Model - self.DP_high)])

            if self.Verbosity > 0:
                PrintDPs()
                print('Max pressure drop error [outer loop] is', max_error_DP, 'Pa')

        if self.Verbosity > 1:
            print('Capacity: ',             self.Capacity)
            print('COP: ',                  self.COP)
            print('COP (w/ both fans): ',   self.COSP)
            print('SHR: ',                  self.SHR)

        print('-------------------------------------')
        print('     Simulation Completed            ')
        print('-------------------------------------')

        return


class VariableSpeedHPClass():
    def __init__(self):
        """
        Load up the necessary sub-structures to be filled with
        the code that follows
        """
        self.LineSetEvapAccumulator             = LineSetOptionClass()
        self.LineSetSuction                     = LineSetOptionClass()
        self.LineSetDischarge                   = LineSetOptionClass()
        self.LineSetLiquid                      = LineSetOptionClass()
        self.ExpDev                             = ExpDevClass()
        self.SuctionAccumulator                 = SuctionAccumulatorClass()

    def Update(self):
        """
        Update cycle class with selected compressor type
        Update cycle class with selected HX type
        Update cycle class with Abstract State
        """
        # ----------------------------------------------------------------------------------
        if self.CompModel == 'Miranda-Mendoza':
            if self.CompType == 'Hitachi-RollingPiston':
                self.Compressor                 = HitachiVariableSpeedRotaryCompressorClass()
            else:
                raise
        else:
            raise

        # ----------------------------------------------------------------------------------
        if self.EvapSolver == 'Moving-Boundary':
            if self.EvapType == 'Fin-tube':
                self.Evaporator                 = EvaporatorClass()
                self.Evaporator.Fins            = FinInputs()
            elif self.EvapType == 'Micro-channel':
                self.Evaporator                 = MicroChannelEvaporatorClass()
                self.Evaporator.Fins            = MicroFinInputs()
            else:
                raise
        elif self.EvapSolver == 'Finite-Element':
            raise Exception('Discretized HX model not implemented')
        else:
            raise

        # ----------------------------------------------------------------------------------
        if self.CondSolver == 'Moving-Boundary':
            if self.CondType == 'Fin-tube':
                self.Condenser                  = CondenserClass()
                self.Condenser.Fins             = FinInputs()
            elif self.CondType == 'Micro-channel':
                self.Condenser                  = MicroCondenserClass()
                self.Condenser.Fins             = MicroFinInputs()
            else:
                raise
        elif self.CondSolver == 'Finite-Element':
            raise Exception('Discretized HX model not implemented')
        else:
            raise

        # ----------------------------------------------------------------------------------
        # Abstract State
        self.AS                                 = CP.AbstractState(self.Backend, self.Ref)
        if hasattr(self, 'MassFrac'):
            self.AS.set_mass_fractions([self.MassFrac,1-self.MassFrac])
            self.AS.build_phase_envelope("dummy")
            self.AS.get_phase_envelope_data()
        elif hasattr(self, 'VoluFrac'):
            self.AS.set_volu_fractions([self.VoluFrac])
            self.AS.set_mass_fractions([self.VoluFrac,1-self.VoluFrac])
            self.AS.build_phase_envelope("dummy")
            self.AS.get_phase_envelope_data()

    def OutputList(self):
        """
            Return a list of parameters for this component for further output

            It is a list of tuples, and each tuple is formed of items:
                [0] Description of value
                [1] Units of value
                [2] The value itself
        """
        Output_List                             = []
        # append optional parameters, if applicable
        if hasattr(self, 'TestName'):
            Output_List.append(('Name',                     'N/A',  self.TestName))
        if hasattr(self, 'TestDescription'):
            Output_List.append(('Description',              'N/A',  self.TestDescription))
        if hasattr(self, 'TestDetails'):
            Output_List.append(('Details',                  'N/A',  self.TestDetails))
        if self.ImposedVariable == 'Charge':
            Output_List.append(('Charge w/no Correction',   'kg',   self.Charge_noCorrect))
            Output_List.append(('Charge w/ Correction',     'kg',   self.Charge_correct))

        Output_List_default = [                             # default output list
            ('Charge',                      'kg',           self.Charge),
            ('Condensation temp (dew)',     'K',            self.T_dew_cond),
            ('Evaporation temp (dew)',      'K',            self.T_dew_evap),
            ('Condenser Subcooling',        'K',            self.DT_sc),
            ('Evaporator Superheat',        'K',            self.DT_sh),
            ('Primary Ref.',                '-',            self.Ref),
            ('COP',                         '-',            self.COP),
            ('COSP',                        '-',            self.COSP),
            ('Net Capacity',                'W',            self.Capacity),
            ('Net Power',                   'W',            self.Power),
            ('EER',                         'BTU/(W-hr)',   self.EER),
            ('SHR',                         '-',            self.SHR),
            ('Imposed Variable',            '-',            self.ImposedVariable),
         ]
        for i in range(0,len(Output_List_default)):         # append default parameters to output list
            Output_List.append(Output_List_default[i])
        return Output_List

    def Calculate(self, DT_evap, DT_cond, DT_sh):
        """
        Inputs are differences in temperature [K] between HX air inlet temperature
        and the dew temperature for the heat exchanger.

        Required Inputs:
            DT_evap:
                Difference in temperature [K] between evaporator air inlet temperature and refrigerant dew temperature
            DT_cond:
                Difference in temperature [K] between condenser air inlet temperature and refrigeant dew temperature
            DT_sh:
                Superheat value
        """
        if self.Verbosity > 1:
            print('DT_evap %7.4f DT_cond %7.4f DTsh %7.4f,' % (DT_evap, DT_cond, DT_sh))

        # AbstractState
        AS                          = self.AS

        # Surrounding
        if hasattr(self, 'T_0') and hasattr(self, 'p_0'):
            T_0                     = self.T_0
            p_0                     = self.p_0
        else:
            T_0                     = 25+273.15                 # K
            p_0                     = 101325                    # Pa
            self.T_0                = T_0
            self.p_0                = p_0

        AS.update(CP.PT_INPUTS, p_0, T_0)
        h_0                         = AS.hmass()
        s_0                         = AS.smass()

        self.h_0                    = h_0
        self.s_0                    = s_0

        # Condenser and evaporator dew temperature (guess)
        T_dew_cond                  = self.Condenser.Fins.Air.T_db + DT_cond    # the values (T_in_a,..) come from line 128ff
        T_dew_evap                  = self.Evaporator.Fins.Air.T_db - DT_evap
        # Condenser and evaporator saturation pressures
        AS.update(CP.QT_INPUTS, 1.0, T_dew_cond)
        p_sat_cond                  = AS.p()                    # [Pa]
        AS.update(CP.QT_INPUTS, 1.0, T_dew_evap)
        p_sat_evap                  = AS.p()                    # [Pa]
        # evaporator bubble temperature
        AS.update(CP.PQ_INPUTS, p_sat_evap, 0.0)
        T_bubble_evap               = AS.T()                    # [T]

        self.T_dew_cond             = T_dew_cond
        self.T_dew_evap             = T_dew_evap

        # If the user doesn't include the Mode, fail
        assert hasattr(self, 'Mode')

        # -------------------------------------------------------------------------------
        # Cycle Solver in 'AC' model
        # -------------------------------------------------------------------------------
        if self.Mode == 'AC':
            if not hasattr(self.Compressor, 'm_dot_r') or self.Compressor.m_dot_r < 0.00001:
                # The first run of model, run the compressor just so you can get a preliminary value
                # for the mass flow rate for the line set
                # -------------------------------------------------------------------------------
                params  = {                 # dictionary -> key:value, e.g. 'key':2345,
                    'p_in_r':               p_sat_evap,
                    'p_out_r':              p_sat_cond,
                    'T_in_r':               T_dew_evap + DT_sh,
                    'AS':                   AS,
                }
                self.Compressor.Update(**params)
                self.Compressor.Calculate()

            # -------------------------------------------------------------------------------
            # Calculate inlet enthalpy
            AS.update(CP.PT_INPUTS, p_sat_evap, T_dew_evap + DT_sh)
            h_in                            = AS.hmass()                # [J/kg]
            params  = {
                'T_in':                     T_dew_evap+DT_sh,
                'p_in':                     p_sat_evap,
                'h_in':                     h_in,
                'm_dot':                    self.Compressor.m_dot_r,
                'AS':                       AS,
                'Name':                     'Evap Accumulator Line'
            }
            self.LineSetEvapAccumulator.Update(**params)
            self.LineSetEvapAccumulator.Calculate()

            # -------------------------------------------------------------------------------
            params  = {
                'T_in_r':                   self.LineSetEvapAccumulator.T_out,
                'p_in_r':                   p_sat_evap,
                'AS':                       AS,
            }
            self.SuctionAccumulator.Update(**params)
            self.SuctionAccumulator.Calculate()

            # -------------------------------------------------------------------------------
            params  = {
                'T_in':                     self.SuctionAccumulator.T_out_r,
                'p_in':                     p_sat_evap,
                'h_in':                     self.SuctionAccumulator.h_out_r,
                'm_dot':                    self.Compressor.m_dot_r,
                'AS':                       AS,
                'Name':                     'Suction Line'
            }
            self.LineSetSuction.Update(**params)
            self.LineSetSuction.Calculate()

            # -------------------------------------------------------------------------------
            params  = {                     # dictionary -> key:value, e.g. 'key':2345,
                'p_in_r':                   p_sat_evap - self.DP_low,
                'p_out_r':                  p_sat_cond + self.DP_high,
                'T_in_r':   TrhoPhase_ph(self.AS, p_sat_evap, self.LineSetSuction.h_out, T_bubble_evap, T_dew_evap)[0],
                'AS':                       AS,
            }
            self.Compressor.Update(**params)
            self.Compressor.Calculate()

            if self.Verbosity > 1:
                print('Comp DP L H', self.DP_low, self.DP_high)

            # -------------------------------------------------------------------------------
            params  = {
                'T_in':                     self.Compressor.T_out_r,
                'p_in':                     p_sat_cond,
                'h_in':                     self.Compressor.h_out_r,
                'm_dot':                    self.Compressor.m_dot_r,
                'AS':                       AS,
                'Name':                     'Discharge Line'
            }
            self.LineSetDischarge.Update(**params)
            self.LineSetDischarge.Calculate()

            # -------------------------------------------------------------------------------
            params  = {
                'm_dot_r':                  self.Compressor.m_dot_r,
                'T_in_r':                   self.LineSetDischarge.T_out,
                'p_sat_r':                  p_sat_cond,
                'AS':                       AS,
            }
            self.Condenser.Update(**params)
            self.Condenser.Calculate()

            # -------------------------------------------------------------------------------
            params  = {
                'T_in':                     self.Condenser.T_out_r,
                'p_in':                     p_sat_cond,
                'h_in':                     self.Condenser.h_out_r,
                'm_dot':                    self.Compressor.m_dot_r,
                'AS':                       AS,
                'Name':                     'Liquid Line'
            }
            self.LineSetLiquid.Update(**params)
            self.LineSetLiquid.Calculate()

            # -------------------------------------------------------------------------------
            params  = {
                'p_in_r':                   p_sat_cond,
                'h_in_r':                   self.LineSetLiquid.h_out,
                'T_sup':                    DT_sh,
                'p_out_r':                  p_sat_evap,
                'AS':                       AS,
            }
            self.ExpDev.Update(**params)
            self.ExpDev.Calculate()

            # -------------------------------------------------------------------------------
            params  = {
                'm_dot_r':                  self.Compressor.m_dot_r,
                'p_sat_r':                  p_sat_evap,
                'h_in_r':                   self.ExpDev.h_out_r,
                'AS':                       AS,
            }
            self.Evaporator.Update(**params)
            self.Evaporator.Calculate()

            # -------------------------------------------------------------------------------
            self.EnergyBalance  = self.Compressor.CycleEnergyIn + self.Condenser.Q + self.Evaporator.Q + self.LineSetEvapAccumulator.Q \
                                  + self.LineSetSuction.Q + self.LineSetDischarge.Q + self.LineSetLiquid.Q + self.SuctionAccumulator.Q

            Charge              = self.Compressor.Charge + self.Condenser.Charge + self.Evaporator.Charge + self.LineSetEvapAccumulator.Charge \
                                  + self.LineSetSuction.Charge + self.LineSetDischarge.Charge + self.LineSetLiquid.Charge + self.SuctionAccumulator.Charge

            self.Charge_noCorrect           = Charge                    # Charge with no correction
            # Correct the charge with Bo Shen's two point regression method
            if self.ImposedVariable == 'Charge':
                C                           = self.C                    # [kg]
                if self.ChargeMethod == 'One-point':
                    self.Charge_correct     = Charge + C                # one-point charge correction
                elif self.ChargeMethod == 'Two-point':
                    K                       = self.K                    # [kg]
                    w_ref                   = self.w_ref                # [-]
                    w_liq                   = self.Condenser.w_subcool
                    self.delta_charge       = C + K * (w_liq - w_ref)
                    self.Charge_correct     = Charge + self.delta_charge # two-point charge correction
                self.Charge                 = self.Charge_correct
            else:
                self.Charge                 = Charge

            # -------------------------------------------------------------------------------
            if self.ExpDev.ExpType == 'Ideal':
                resid                       = np.zeros((3))
                resid[0]                    = self.Compressor.m_dot_r * (self.Evaporator.h_out_r - self.LineSetSuction.h_in)
                resid[2]                    = DT_sh - self.Evaporator.DT_sh
            else:
                resid                       = np.zeros((3))
                resid[0]                    = self.Evaporator.h_out_r - self.LineSetEvapAccumulator.h_in
                resid[2]                    = self.Compressor.m_dot_r - self.ExpDev.m_dot_r

            if self.ImposedVariable == 'Subcooling':
                resid[1]                    = self.Condenser.DT_sc - self.DT_sc_target
            elif self.ImposedVariable == 'Charge':
                resid[1]                    = self.Charge - self.Charge_target

            if self.Verbosity > 1:
                print(resid)

            # -------------------------------------------------------------------------------
            self.Capacity                   = self.Evaporator.Capacity
            self.Power                      = self.Compressor.W + self.Evaporator.Fins.Air.FanPower + self.Condenser.Fins.Air.FanPower
            self.COP                        = self.Evaporator.Q / self.Compressor.W
            self.COSP                       = self.Evaporator.Capacity / self.Power
            self.EER                        = W2BTUh(self.Capacity) / self.Power                # [BTU/(W-hr)]
            self.SHR                        = self.Evaporator.SHR
            self.DT_sc                      = self.Condenser.DT_sc
            self.DP_HighPressure            = self.Condenser.DP_r + self.LineSetLiquid.DP + self.LineSetDischarge.DP
            self.DP_LowPressure             = self.Evaporator.DP_r + self.LineSetEvapAccumulator.DP + self.LineSetSuction.DP

        # -------------------------------------------------------------------------------
        # Cycle Solver in 'HP' model
        # -------------------------------------------------------------------------------
        elif self.Mode == 'HP':
            if not hasattr(self.Compressor, 'm_dot_r') or self.Compressor.m_dot_r < 0.00001:
                # The first run of model, run the compressor just so you can get a preliminary value
                # for the mass flow rate for the line set
                params  = {                 # dictionary -> key:value, e.g. 'key':2345,
                    'p_in_r':               p_sat_evap,
                    'p_out_r':              p_sat_cond,
                    'T_in_r':               T_dew_evap + DT_sh,
                    'AS':                   AS,
                }
                self.Compressor.Update(**params)
                self.Compressor.Calculate()

            # -------------------------------------------------------------------------------
            # Calculate inlet enthalpy
            AS.update(CP.PT_INPUTS, p_sat_evap, T_dew_evap + DT_sh)
            h_in                            = AS.hmass()                # [J/kg]

            params  = {
                'T_in':                     T_dew_evap+DT_sh,
                'p_in':                     p_sat_evap,
                'h_in':                     h_in,
                'm_dot':                    self.Compressor.m_dot_r,
                'AS':                       AS,
                'Name':                     'Evap Accumulator Line'
            }
            self.LineSetEvapAccumulator.Update(**params)
            self.LineSetEvapAccumulator.Calculate()

            # -------------------------------------------------------------------------------
            params  = {
                'T_in':                     self.LineSetEvapAccumulator.T_out,
                'p_in_r':                   p_sat_evap,
                'AS':                       AS,
            }
            self.SuctionAccumulator.Update(**params)
            self.SuctionAccumulator.Calculate()

            # -------------------------------------------------------------------------------
            params  = {
                'T_in':                     self.SuctionAccumulator.T_out_r,
                'p_in':                     p_sat_evap,
                'h_in':                     self.SuctionAccumulator.h_out_r,
                'm_dot':                    self.Compressor.m_dot_r,
                'AS':                       AS,
                'Name':                     'Suction Line'
            }
            self.LineSetSuction.Update(**params)
            self.LineSetSuction.Calculate()

            # -------------------------------------------------------------------------------
            params  = {                     # dictionary -> key:value, e.g. 'key':2345,
                'p_in_r':                   p_sat_evap - self.DP_low,
                'p_out_r':                  p_sat_cond + self.DP_high,
                'T_in_r':   TrhoPhase_ph(self.AS, p_sat_evap, self.LineSetSuction.h_out, T_bubble_evap, T_dew_evap)[0],
                'Ref':                      self.Ref,
                'AS':                       AS
            }
            self.Compressor.Update(**params)
            self.Compressor.Calculate()

            # -------------------------------------------------------------------------------
            params  = {
                'T_in':                     self.Compressor.T_out_r,
                'p_in':                     p_sat_cond,
                'h_in':                     self.Compressor.h_out_r,
                'm_dot':                    self.Compressor.m_dot_r,
                'AS':                       AS,
                'Name':                     'Discharge Line'
            }
            self.LineSetDischarge.Update(**params)
            self.LineSetDischarge.Calculate()

            # -------------------------------------------------------------------------------
            params  = {
                'm_dot_r':                  self.Compressor.m_dot_r,
                'T_in_r':                   self.LineSetDischarge.T_out,
                'p_sat_r':                  p_sat_cond,
                'AS':                       AS,
            }
            self.Condenser.Update(**params)
            self.Condenser.Calculate()

            # -------------------------------------------------------------------------------
            params  = {
                'T_in':                     self.Condenser.T_out_r,
                'p_in':                     p_sat_cond,
                'h_in':                     self.Condenser.h_out_r,
                'm_dot':                    self.Compressor.m_dot_r,
                'AS':                       AS,
                'Name':                     'Liquid Line'
            }
            self.LineSetLiquid.Update(**params)
            self.LineSetLiquid.Calculate()

            # -------------------------------------------------------------------------------
            params  = {
                'p_in_r':                   p_sat_cond,
                'h_in_r':                   self.LineSetLiquid.h_out,
                'T_sup':                    DT_sh,
                'p_out_r':                  p_sat_evap,
                'AS':                       AS,
            }
            self.ExpDev.Update(**params)
            self.ExpDev.Calculate()

            # -------------------------------------------------------------------------------
            params  = {
                'm_dot_r':                  self.Compressor.m_dot_r,
                'p_sat_r':                  p_sat_evap,
                'h_in_r':                   self.ExpDev.h_out_r,
                'AS':                       AS,
            }
            self.Evaporator.Update(**params)
            self.Evaporator.Calculate()

            # -------------------------------------------------------------------------------
            # Energy and Mass conservation
            self.EnergyBalance  = self.Compressor.CycleEnergyIn + self.Condenser.Q + self.Evaporator.Q + self.LineSetEvapAccumulator.Q \
                                  + self.LineSetSuction.Q + self.LineSetDischarge.Q + self.LineSetLiquid.Q + self.SuctionAccumulator.Q

            Charge              = self.Compressor.Charge + self.Condenser.Charge + self.Evaporator.Charge + self.LineSetEvapAccumulator.Charge \
                                  + self.LineSetSuction.Charge + self.LineSetDischarge.Charge + self.LineSetLiquid.Charge + self.SuctionAccumulator.Charge

            self.Charge_noCorrect           = Charge                    # Charge with no correction
            # Correct the charge with Bo Shen's two point regression method
            if self.ImposedVariable == 'Charge':
                C                           = self.C                    # [kg]
                if self.ChargeMethod == 'One-point':
                    self.Charge_correct     = Charge + C                # one-point charge correction
                elif self.ChargeMethod=='Two-point':
                    K                       = self.K                    # [kg]
                    w_ref                   = self.w_ref                # [-]
                    w_liq                   = self.Condenser.w_subcool
                    self.delta_charge       = C + K * (w_liq - w_ref)
                    self.Charge_correct     = Charge + self.delta_charge # two-point charge correction
                self.Charge                 = self.Charge_correct
            else:
                self.Charge                 = Charge

            # -------------------------------------------------------------------------------
            if self.ExpDev.ExpType == 'Ideal':
                resid                       = np.zeros((3))
                resid[0]                    = self.Compressor.m_dot_r * (self.Evaporator.h_out_r - self.LineSetSuction.h_in)
                resid[2]                    = DT_sh - self.Evaporator.DT_sh
            else:
                resid                       = np.zeros((3))
                resid[0]                    = self.Evaporator.h_out_r - self.LineSetEvapAccumulator.h_in
                resid[2]                    = self.Compressor.m_dot_r - self.ExpDev.m_dot_r

            if self.ImposedVariable == 'Subcooling':
                resid[1]                    = self.Condenser.DT_sc - self.DT_sc_target
            elif self.ImposedVariable == 'Charge':
                resid[1]                    = self.Charge - self.Charge_target

            # -------------------------------------------------------------------------------
            self.Capacity                   = -self.Condenser.Q + self.Condenser.Fins.Air.FanPower
            self.DT_sc                      = self.Condenser.DT_sc
            self.Power  = self.Compressor.W + self.Evaporator.Fins.Air.FanPower + self.Condenser.Fins.Air.FanPower
            self.COP                        = -self.Condenser.Q / self.Compressor.W
            self.COSP                       = self.Capacity / self.Power
            self.EER                        = W2BTUh(self.Capacity) / self.Power        # [BTU/(W-hr)]
            self.SHR                        = self.Evaporator.SHR
            self.DP_HighPressure            = self.Condenser.DP_r + self.LineSetLiquid.DP + self.LineSetDischarge.DP
            self.DP_LowPressure             = self.Evaporator.DP_r + self.LineSetEvapAccumulator.DP + self.LineSetSuction.DP

        # -------------------------------------------------------------------------------
        # Other options
        # -------------------------------------------------------------------------------
        else:
            ValueError("DX Cycle mode must be 'AC', or 'HP'")

        # -------------------------------------------------------------------------------
        if self.Verbosity > 1:
            print('DT_evap {:4.4f} DT_cond {:4.4f} DTsh {:4.4f} resid[0] {:12.6e} resid[1]: {:12.6e} resid[3]: {:12.6e} Charge {:4.4f} SC: {:4.4f}'
                  .format(DT_evap, DT_cond, DT_sh, resid[0], resid[1], resid[2], self.Charge, self.Condenser.DT_sc))
        self.DT_evap            = DT_evap
        self.DT_cond            = DT_cond
        self.DT_sh              = DT_sh

        return resid

    def PreconditionedSolve(self):
        """
        Solver that will precondition by trying a range of DeltaT until the model
        can solve, then will kick into 2-D Newton Raphson solve

        The two input variables for the system solver are the differences in
        temperature between the inlet air temperature of the heat exchanger and the
        dew temperature of the refrigerant.  This is important for refrigerant blends
        with temperature glide during constant-pressure evaporation or condensation.
        Good examples of common working fluid with glide would be R404A or R410A.
        """
        def OBJECTIVE_DXCycle(x):
            """
            A wrapper function to convert input vector for fsolve to the proper form for the solver
            """
            try:
                resids=self.Calculate(DT_evap=float(x[0]),DT_cond=float(x[1]),DT_sh=float(x[2]))#,DP_low=float(x[2]),DP_high=float(x[3]))
            except ValueError:
                raise
            return resids

        # Use the preconditioner to determine a reasonably good starting guess
        print ('-------------------------------------')
        print ('            Running Preconditioner   ')
        print ('-------------------------------------')
        DT_evap_init,DT_cond_init,DT_sh_init=DXPreconditioner(self)

        print ('-------------------------------------')
        print ('           Preconditioner Completed  ')
        print ('             Starting Main Cycle     ')
        print ('-------------------------------------')


        GoodRun=False
        while GoodRun==False:
            try:
                self.DP_low=0
                self.DP_high=0
                DP_converged=False
                while DP_converged==False:
                    # Actually run the Newton-Raphson solver to get the solution
                    x=Broyden(OBJECTIVE_DXCycle,[DT_evap_init,DT_cond_init,DT_sh_init])
                    delta_low=abs(self.DP_low-abs(self.DP_LowPressure))
                    delta_high=abs(self.DP_high-abs(self.DP_HighPressure))
                    self.DP_low=abs(self.DP_LowPressure)
                    self.DP_high=abs(self.DP_HighPressure)
                    # Update the guess values based on last converged values
                    DT_evap_init=self.DT_evap
                    DT_cond_init=self.DT_cond
                    DT_sh_init=self.DT_sh
                    if delta_low<1 and delta_high<1:
                        DP_converged=True
                    if self.Verbosity>4:
                        print (self.DP_HighPressure,self.DP_LowPressure,'DPHP')
                    GoodRun=True
            except AttributeError:
                # This will be a fatal error !! Should never have attribute error
                raise
            except:
                print ("--------------  Exception Caught ---------------- ")
                print ("Error of type",sys.exc_info()[0]," is: " + sys.exc_info()[1].message)
                raise

        if self.Verbosity>0:

            print ("Debugging Results")
            print ()
            print ('Q Evaporator:', self.Evaporator.Q)
            print ('Capacity: ', self.Capacity)
            print ('Compressor Power:',self.Compressor.W)
            print ('Total Power', self.Power)
            print ('COP: ',self.COP)
            print ('COP (w/ both fans): ',self.COSP)
            print ('SHR: ',self.SHR)

        print ('-------------------------------------')
        print ('     Simulation Completed            ')
        print ('-------------------------------------')


# ------------------------------------------------------------------------
class BTMS_DXCycleClass():
    def __init__(self):
        """
        Load up the necessary sub-structures to be filled with
        the code that follows
        """
        self.Compressor                 = CompressorClass()
        self.LineSetSuction             = LineSetOptionClass()
        self.LineSetDischarge           = LineSetOptionClass()
        self.LineSetLiquid              = LineSetOptionClass()
        self.ExpDev                     = ExpDevClass()

    def Update(self):
        """
        Update cycle class with selected HX type
        Update cycle class with Abstract State
        """
        if self.EvapSolver == 'Moving-Boundary':
            if self.EvapType == 'Fin-tube':
                self.Evaporator         = EvaporatorClass()
                self.Evaporator.Fins    = FinInputs()
            elif self.EvapType == 'Micro-channel':
                self.Evaporator         = MicroChannelEvaporatorClass()
                self.Evaporator.Fins    = MicroFinInputs()
            else:
                raise
        elif self.EvapSolver == 'Finite-Element':
            raise Exception('Discretized HX model not implemented')
        else:
            raise

        if self.CondSolver == 'Moving-Boundary':
            if self.CondType == 'Fin-tube':
                self.Condenser          = CondenserClass()
                self.Condenser.Fins     = FinInputs()
            elif self.CondType == 'Micro-channel':
                self.Condenser          = MicroCondenserClass()
                self.Condenser.Fins     = MicroFinInputs()
            else:
                raise
        elif self.CondSolver == 'Finite-Element':
            raise Exception('Discretized HX model not implemented')
        else:
            raise

        # Abstract State
        self.AS                         = CP.AbstractState(self.Backend, self.Ref)
        if hasattr(self, 'MassFrac'):
            self.AS.set_mass_fractions([self.MassFrac, 1 - self.MassFrac])
            self.AS.build_phase_envelope("dummy")
            self.AS.get_phase_envelope_data()
        elif hasattr(self, 'VoluFrac'):
            self.AS.set_volu_fractions([self.VoluFrac])
            self.AS.set_mass_fractions([self.VoluFrac, 1 - self.VoluFrac])
            self.AS.build_phase_envelope("dummy")
            self.AS.get_phase_envelope_data()

    def OutputList(self):
        """
            Return a list of parameters for this component for further output

            It is a list of tuples, and each tuple is formed of items:
                [0] Description of value
                [1] Units of value
                [2] The value itself
        """
        Output_List                     = []
        # append optional parameters, if applicable
        if hasattr(self, 'TestName'):
            Output_List.append(('Name',                     'N/A',  self.TestName))
        if hasattr(self, 'TestDescription'):
            Output_List.append(('Description',              'N/A',  self.TestDescription))
        if hasattr(self, 'TestDetails'):
            Output_List.append(('Details',                  'N/A',  self.TestDetails))
        if self.ImposedVariable == 'Charge':
            Output_List.append(('Charge w/no Correction',   'kg',   self.Charge_noCorrect))
            Output_List.append(('Charge w/ Correction',     'kg',   self.Charge_correct))
        Output_List_default = [                                     # default output list
            ('Charge',                      'kg',                   self.Charge),
            ('Condensation temp (dew)',     'K',                    self.T_dew_cond),
            ('Evaporation temp (dew)',      'K',                    self.T_dew_evap),
            ('Condenser Subcooling',        'K',                    self.DT_sc),
            ('Evaporator Superheat',        'K',                    self.DT_sh),
            ('Primary Ref.',                '-',                    self.Ref),
            ('COP',                         '-',                    self.COP),
            ('COSP',                        '-',                    self.COSP),
            ('Net Capacity',                'W',                    self.Capacity),
            ('Net Power',                   'W',                    self.Power),
            ('EER',                         'BTU/(W-hr)',           self.EER),
            ('SHR',                         '-',                    self.SHR),
            ('Imposed Variable',            '-',                    self.ImposedVariable),
         ]
        for i in range(0, len(Output_List_default)):                # append default parameters to output list
            Output_List.append(Output_List_default[i])
        return Output_List

    def Calculate(self, DT_evap, DT_cond, DT_sh):
        """
        Inputs are differences in temperature [K] between HX air inlet temperature
        and the dew temperature for the heat exchanger.

        Required Inputs:
            DT_evap:
                Difference in temperature [K] between evaporator air inlet temperature and refrigerant dew temperature
            DT_cond:
                Difference in temperature [K] between condenser air inlet temperature and refrigerant dew temperature
            DT_sh:
                Superheat value
        """
        if self.Verbosity > 1:
            print('DT_evap %7.4f DT_cond %7.4f DT_sh %7.4f,' %(DT_evap, DT_cond, DT_sh))

        # AbstractState
        AS                                  = self.AS

        # Surrounding
        if hasattr(self, 'T_0') and hasattr(self, 'p_0'):
            T_0                             = self.T_0
            p_0                             = self.p_0
        else:
            T_0                             = 25+273.15                 # K
            p_0                             = 101325                    # Pa
            self.T_0                        = T_0
            self.p_0                        = p_0

        AS.update(CP.PT_INPUTS, p_0, T_0)
        h_0                                 = AS.hmass()
        s_0                                 = AS.smass()

        self.h_0                            = h_0
        self.s_0                            = s_0

        # Condenser and evaporator dew temperature (guess)
        T_dew_cond                          = self.Condenser.Fins.Air.T_db + DT_cond    # the values (T_in_a,..) come from line 128ff
        T_dew_evap                          = self.Evaporator.Fins.Air.T_db - DT_evap
        # Condenser and evaporator saturation pressures
        AS.update(CP.QT_INPUTS, 1.0, T_dew_cond)
        p_sat_cond                          = AS.p()                    # [Pa]
        AS.update(CP.QT_INPUTS, 1.0, T_dew_evap)
        p_sat_evap                          = AS.p()                    # [Pa]
        # evaporator bubble temperature
        AS.update(CP.PQ_INPUTS, p_sat_evap, 0.0)
        T_bubble_evap                       = AS.T()                    # [T]

        self.T_dew_cond                     = T_dew_cond
        self.T_dew_evap                     = T_dew_evap

        # If the user doesn't include the Mode, fail
        assert hasattr(self, 'Mode')

        # --------------------------------------------------------------------------------
        # Cycle Solver in 'AC' model
        # --------------------------------------------------------------------------------
        if self.Mode == 'AC':
            # -------------------------------------------------------------------------------------
            if not hasattr(self.Compressor, 'm_dot_r') or self.Compressor.m_dot_r < 0.00001:
                # The first run of model, run the compressor just so you can get a preliminary value
                # for the mass flow rate for the line set
                params  = {                     # dictionary -> key:value, e.g. 'key':2345,
                    'p_in_r':               p_sat_evap,
                    'p_out_r':              p_sat_cond,
                    'T_in_r':               T_dew_evap + DT_sh,
                    'AS':                   AS,
                }
                self.Compressor.Update(**params)
                self.Compressor.Calculate()

            # -------------------------------------------------------------------------------------
            # Calculate inlet enthalpy
            AS.update(CP.PT_INPUTS, p_sat_evap, T_dew_evap + DT_sh)
            h_in                            = AS.hmass()                # [J/kg]
            params  = {
                'T_in':                     T_dew_evap + DT_sh,
                'p_in':                     p_sat_evap,
                'h_in':                     h_in,
                'm_dot':                    self.Compressor.m_dot_r,
                'AS':                       AS,
                'Name': '                   Suction Line'
            }
            self.LineSetSuction.Update(**params)
            self.LineSetSuction.Calculate()

            # -------------------------------------------------------------------------------------
            params  = {                     # dictionary -> key:value, e.g. 'key':2345,
                'p_in_r':                   p_sat_evap - self.DP_low,
                'p_out_r':                  p_sat_cond + self.DP_high,
                'T_in_r':                   TrhoPhase_ph(self.AS, p_sat_evap, self.LineSetSuction.h_out, T_bubble_evap, T_dew_evap)[0],
                'AS':                       AS,
            }
            self.Compressor.Update(**params)
            self.Compressor.Calculate()

            if self.Verbosity > 1:
                print('Comp DP L H', self.DP_low, self.DP_high)

            # -------------------------------------------------------------------------------------
            params  = {
                'T_in':                     self.Compressor.T_out_r,
                'p_in':                     p_sat_cond,
                'h_in':                     self.Compressor.h_out_r,
                'm_dot':                    self.Compressor.m_dot_r,
                'AS':                       AS,
                'Name':                     'Discharge Line'
            }
            self.LineSetDischarge.Update(**params)
            self.LineSetDischarge.Calculate()

            # -------------------------------------------------------------------------------------
            params  = {
                'm_dot_r':                  self.Compressor.m_dot_r,
                'T_in_r':                   self.LineSetDischarge.T_out,
                'p_sat_r':                  p_sat_cond,
                'AS':                       AS,
            }
            self.Condenser.Update(**params)
            self.Condenser.Calculate()

            # -------------------------------------------------------------------------------------
            params  = {
                'T_in':                     self.Condenser.T_out_r,
                'p_in':                     p_sat_cond,
                'h_in':                     self.Condenser.h_out_r,
                'm_dot':                    self.Compressor.m_dot_r,
                'AS':                       AS,
                'Name':                     'Liquid Line'
            }
            self.LineSetLiquid.Update(**params)
            self.LineSetLiquid.Calculate()

            # -------------------------------------------------------------------------------------
            params  = {
                'p_in_r':                   p_sat_cond,
                'h_in_r':                   self.LineSetLiquid.h_out,
                'T_sup':                    DT_sh,
                'p_out_r':                  p_sat_evap,
                'AS':                       AS,
            }
            self.ExpDev.Update(**params)
            self.ExpDev.Calculate()

            # ------------------------------change to wavychannel batteries---------------------------------------------------
            params  = {
                'm_dot_r':                  self.Compressor.m_dot_r,
                'p_sat_r':                  p_sat_evap,
                'h_in_r':                   self.ExpDev.h_out_r,
                'AS':                       AS,
            }
            self.Evaporator.Update(**params)
            self.Evaporator.Calculate()

            # -------------------------------------------------------------------------------------
            # Energy and Mass conservation
            self.EnergyBalance = self.Compressor.CycleEnergyIn + self.Condenser.Q + self.Evaporator.Q + self.LineSetSuction.Q \
                                 + self.LineSetDischarge.Q + self.LineSetLiquid.Q

            Charge              = self.Compressor.Charge + self.Condenser.Charge + self.Evaporator.Charge + self.LineSetSuction.Charge \
                                  + self.LineSetDischarge.Charge + self.LineSetLiquid.Charge

            # Charge with no correction
            self.Charge_noCorrect           = Charge

            # Correct the charge with Bo Shen's two point regression method
            if self.ImposedVariable == 'Charge':
                C                           = self.C                                # [kg]
                if self.ChargeMethod == 'One-point':
                    self.Charge_correct = Charge + C                                # one-point charge correction
                elif self.ChargeMethod == 'Two-point':
                    K                       = self.K                                # [kg]
                    w_ref                   = self.w_ref                            # [-]
                    w_liq                   = self.Condenser.w_subcool
                    self.delta_charge       = C + K * (w_liq - w_ref)
                    self.Charge_correct     = Charge + self.delta_charge            # two-point charge correction
                self.Charge                 = self.Charge_correct
            else:
                self.Charge                 = Charge

            # -------------------------------------------------------------------
            if self.ExpDev.ExpType == 'Ideal':
                resid                       = np.zeros((3))
                resid[0]                    = self.Compressor.m_dot_r * (self.Evaporator.h_out_r - self.LineSetSuction.h_in)
                resid[2]                    = DT_sh - self.Evaporator.DT_sh
            else:
                resid                       = np.zeros((3))
                resid[0]                    = self.Evaporator.h_out_r - self.LineSetSuction.h_in
                resid[2]                    = self.Compressor.m_dot_r - self.ExpDev.m_dot_r

            # -------------------------------------------------------------------
            if self.ImposedVariable == 'Subcooling':
                resid[1]                    = self.Condenser.DT_sc - self.DT_sc_target
            elif self.ImposedVariable == 'Charge':
                resid[1]                    = self.Charge - self.Charge_target

            if self.Verbosity > 1:
                print(resid)

            self.Capacity                   = self.Evaporator.Capacity
            self.Power                      = self.Compressor.W + self.Evaporator.Fins.Air.FanPower + self.Condenser.Fins.Air.FanPower
            self.COP                        = self.Evaporator.Q / self.Compressor.W
            self.COSP                       = self.Evaporator.Capacity / self.Power
            self.EER                        = W2BTUh(self.Capacity) / self.Power        # [BTU/(W-hr)]
            self.SHR                        = self.Evaporator.SHR
            self.DT_sc                      = self.Condenser.DT_sc
            self.DP_HighPressure            = self.Condenser.DP_r + self.LineSetLiquid.DP + self.LineSetDischarge.DP
            self.DP_LowPressure             = self.Evaporator.DP_r + self.LineSetSuction.DP

        # --------------------------------------------------------------------------------
        # Cycle Solver in 'HP' model
        # --------------------------------------------------------------------------------
        elif self.Mode == 'HP':
            if not hasattr(self.Compressor, 'm_dot_r') or self.Compressor.m_dot_r < 0.00001:
                # The first run of model, run the compressor just so you can get a preliminary value
                # for the mass flow rate for the line set
                # ------------------------------------------------------------------------
                params  = {                     # dictionary -> key:value, e.g. 'key':2345,
                    'p_in_r':               p_sat_evap,
                    'p_out_r':              p_sat_cond,
                    'T_in_r':               T_dew_evap + DT_sh,
                    'AS':                   AS,
                }
                self.Compressor.Update(**params)
                self.Compressor.Calculate()

            # -----------------------------------------------------------------------------
            # Calculate inlet enthalpy
            AS.update(CP.PT_INPUTS, p_sat_evap, T_dew_evap + DT_sh)
            h_in                            = AS.hmass()                    # [J/kg]
            params  = {
                'T_in':                     T_dew_evap + DT_sh,
                'p_in':                     p_sat_evap,
                'h_in':                     h_in,
                'm_dot':                    self.Compressor.m_dot_r,
                'AS':                       AS,
                'Name':                     'Suction Line'
            }
            self.LineSetSuction.Update(**params)
            self.LineSetSuction.Calculate()

            # -----------------------------------------------------------------------------
            params  = {                     # dictionary -> key:value, e.g. 'key':2345,
                'p_in_r':                   p_sat_evap - self.DP_low,
                'p_out_r':                  p_sat_cond + self.DP_high,
                'T_in_r':                   TrhoPhase_ph(self.AS, p_sat_evap, self.LineSetSuction.h_out, T_bubble_evap, T_dew_evap)[0],
                'Ref':                      self.Ref,
                'AS':                       AS
            }
            self.Compressor.Update(**params)
            self.Compressor.Calculate()

            # -----------------------------------------------------------------------------
            params  = {
                'T_in':                     self.Compressor.T_out_r,
                'p_in':                     p_sat_cond,
                'h_in':                     self.Compressor.h_out_r,
                'm_dot':                    self.Compressor.m_dot_r,
                'AS':                       AS,
                'Name':                     'Discharge Line'
            }
            self.LineSetDischarge.Update(**params)
            self.LineSetDischarge.Calculate()

            # -----------------------------------------------------------------------------
            params  = {
                'm_dot_r':                  self.Compressor.m_dot_r,
                'T_in_r':                   self.LineSetDischarge.T_out,
                'p_sat_r':                  p_sat_cond,
                'AS':                       AS,
            }
            self.Condenser.Update(**params)
            self.Condenser.Calculate()

            # -----------------------------------------------------------------------------
            params  = {
                'T_in':                     self.Condenser.T_out_r,
                'p_in':                     p_sat_cond,
                'h_in':                     self.Condenser.h_out_r,
                'm_dot':                    self.Compressor.m_dot_r,
                'AS':                       AS,
                'Name':                     'Liquid Line'
            }
            self.LineSetLiquid.Update(**params)
            self.LineSetLiquid.Calculate()

            # -----------------------------------------------------------------------------
            params  = {
                'p_in_r':                   p_sat_cond,
                'h_in_r':                   self.LineSetLiquid.h_out,
                'T_sup':                    DT_sh,
                'p_out_r':                  p_sat_evap,
                'AS':                       AS,
            }
            self.ExpDev.Update(**params)
            self.ExpDev.Calculate()

            # -----------------------------------------------------------------------------
            params  = {
                'm_dot_r':                  self.Compressor.m_dot_r,
                'p_sat_r':                  p_sat_evap,
                'h_in_r':                   self.ExpDev.h_out_r,
                'AS':                       AS,
            }
            self.Evaporator.Update(**params)
            self.Evaporator.Calculate()

            # -----------------------------------------------------------------------------
            # Energy and Mass conservation

            self.EnergyBalance  = self.Compressor.CycleEnergyIn + self.Condenser.Q + self.Evaporator.Q + self.LineSetSuction.Q \
                                  + self.LineSetDischarge.Q + self.LineSetLiquid.Q

            Charge              = self.Compressor.Charge + self.Condenser.Charge + self.Evaporator.Charge + self.LineSetSuction.Charge \
                                  + self.LineSetDischarge.Charge + self.LineSetLiquid.Charge

            self.Charge_noCorrect           = Charge                    # Charge with no correction
            # Correct the charge with Bo Shen's two point regression method
            if self.ImposedVariable == 'Charge':
                C                           = self.C                    # [kg]
                if self.ChargeMethod == 'One-point':
                    self.Charge_correct     = Charge + C                # one-point charge correction
                elif self.ChargeMethod == 'Two-point':
                    K                       = self.K                    # [kg]
                    w_ref                   = self.w_ref                # [-]
                    w_liq                   = self.Condenser.w_subcool
                    self.delta_charge       = C + K * (w_liq - w_ref)
                    self.Charge_correct     = Charge + self.delta_charge # two-point charge correction
                self.Charge                 = self.Charge_correct
            else:
                self.Charge                 = Charge

            # -----------------------------------------------------------------------------
            if self.ExpDev.ExpType == 'Ideal':
                resid                       = np.zeros((3))
                resid[0]                    = self.Compressor.m_dot_r * (self.Evaporator.h_out_r - self.LineSetSuction.h_in)
                resid[2]                    = DT_sh - self.Evaporator.DT_sh
            else:
                resid                       = np.zeros((3))
                resid[0]                    = self.Evaporator.h_out_r - self.LineSetSuction.h_in
                resid[2]                    = self.Compressor.m_dot_r - self.ExpDev.m_dot_r

            # -----------------------------------------------------------------------------
            if self.ImposedVariable == 'Subcooling':
                resid[1]                    = self.Condenser.DT_sc - self.DT_sc_target
            elif self.ImposedVariable == 'Charge':
                resid[1]                    = self.Charge - self.Charge_target

            # -----------------------------------------------------------------------------
            self.Capacity                   = -self.Condenser.Q + self.Condenser.Fins.Air.FanPower
            self.DT_sc                      = self.Condenser.DT_sc
            self.Power                      = self.Compressor.W + self.Evaporator.Fins.Air.FanPower + self.Condenser.Fins.Air.FanPower
            self.COP                        = -self.Condenser.Q / self.Compressor.W
            self.COSP                       = self.Capacity / self.Power
            self.EER                        = W2BTUh(self.Capacity) / self.Power            # [BTU/(W-hr)]
            self.SHR                        = self.Evaporator.SHR
            self.DP_HighPressure            = self.Condenser.DP_r + self.LineSetDischarge.DP + self.LineSetLiquid.DP
            self.DP_LowPressure             = self.Evaporator.DP_r + self.LineSetSuction.DP

        # --------------------------------------------------------------------------------
        # Other options
        # --------------------------------------------------------------------------------
        else:
            ValueError("DX Cycle mode must be 'AC', or 'HP'")

        if self.Verbosity > 1:
            print('DT_evap {:4.4f} DT_cond {:4.4f} DT_sh {:4.4f} resid[0] {:12.6e} resid[1]: {:12.6e} resid[3]: {:12.6e} Charge {:4.4f} SC: {:4.4f}'
                   .format(DT_evap, DT_cond, DT_sh, resid[0], resid[1], resid[2], self.Charge, self.Condenser.DT_sc))
        self.DT_evap                        = DT_evap
        self.DT_cond                        = DT_cond
        self.DT_sh                          = DT_sh

        return resid

    def PreconditionedSolve(self):
        """
        Solver that will precondition by trying a range of DeltaT until the model
        can solve, then will kick into 2-D Newton Raphson solve

        The two input variables for the system solver are the differences in
        temperature between the inlet air temperature of the heat exchanger and the
        dew temperature of the refrigerant.  This is important for refrigerant blends
        with temperature glide during constant-pressure evaporation or condensation.
        Good examples of common working fluid with glide would be R404A or R410A.
        """
        def OBJECTIVE_DXCycle(x):
            """
            A wrapper function to convert input vector for fsolve to the proper form for the solver
            """
            try:
                resids  = self.Calculate(DT_evap=float(x[0]), DT_cond=float(x[1]), DT_sh=float(x[2])) #,DP_low=float(x[2]),DP_high=float(x[3]))
            except ValueError:
                raise Exception('Something not correct!')
            return resids

        # Use the preconditioner to determine a reasonably good starting guess
        print ('-------------------------------------')
        print ('            Running Preconditioner   ')
        print ('-------------------------------------')
        DT_evap_init, DT_cond_init, DT_sh_init = DXPreconditioner(self)

        print ('-------------------------------------')
        print ('           Preconditioner Completed  ')
        print ('             Starting Main Cycle     ')
        print ('-------------------------------------')


        GoodRun                             = False
        while GoodRun == False:
            try:
                self.DP_low                 = 0
                self.DP_high                = 0
                DP_converged                = False
                while DP_converged == False:
                    # Actually run the Newton-Raphson solver to get the solution
                    x                       = Broyden(OBJECTIVE_DXCycle, [DT_evap_init, DT_cond_init, DT_sh_init])
                    delta_low               = abs(self.DP_low - abs(self.DP_LowPressure))
                    delta_high              = abs(self.DP_high - abs(self.DP_HighPressure))
                    self.DP_low             = abs(self.DP_LowPressure)
                    self.DP_high            = abs(self.DP_HighPressure)
                    # Update the guess values based on last converged values
                    DT_evap_init            = self.DT_evap
                    DT_cond_init            = self.DT_cond
                    DT_sh_init              = self.DT_sh
                    if delta_low < 1 and delta_high < 1:
                        DP_converged        = True
                    if self.Verbosity > 4:
                        print(self.DP_HighPressure, self.DP_LowPressure, 'DPHP')
                    GoodRun                 = True
            except AttributeError:
                # This will be a fatal error !! Should never have attribute error
                raise
            except:
                print ("--------------  Exception Caught ---------------- ")
                print ("Error of type", sys.exc_info()[0], " is: " + sys.exc_info()[1].message)
                raise

        if self.Verbosity > 0:

            print ("Debugging Results")
            print ()
            print ('Capacity: ', self.Capacity)
            print ('COP: ',self.COP)
            print ('COP (w/ both fans): ',self.COSP)
            print ('SHR: ',self.SHR)
            print ('UA_r_evap',self.Evaporator.UA_r)
            print ('UA_a_evap',self.Evaporator.UA_a)
            print ('UA_r_cond',self.Condenser.UA_r)
            print ('UA_a_cond',self.Condenser.UA_a)

        print ('-------------------------------------')
        print ('     Simulation Completed            ')
        print ('-------------------------------------')


class BTMS_SecondaryCycleClass():
    def __init__(self):
        """
        Load up the necessary sub-structures to be filled with
        the code that follows
        """
        self.Compressor                 = CompressorClass()
        self.Pump                       = PumpClass()
        self.WavyChannel                = WavyChannelCoolerClass()
        self.ExpDev                     = ExpDevClass()
        # Add both types of internal heat exchangers
        self.CoaxialIHX                 = CoaxialHXClass()
        self.PHEIHX                     = PHEHXClass()
        self.LineSetSupply              = LineSetOptionClass()
        self.LineSetReturn              = LineSetOptionClass()
        self.LineSetSuction             = LineSetOptionClass()
        self.LineSetDischarge           = LineSetOptionClass()
        self.LineSetLiquid              = LineSetOptionClass()

        # Make IHX an empty class for holding parameters common to PHE and Coaxial IHX
        class struct:
            pass
        self.IHX                        = struct()

    def Update(self):
        """
        Update cycle class with selected HX type
        Update cycle class with Abstract State
        """
        # -----------------------------------------------------------------------------
        if self.EvapSolver == 'Moving-Boundary':
            if self.EvapType == 'Fin-tube':
                self.Evaporator                 = EvaporatorClass()
                self.Evaporator.Fins            = FinInputs()
            elif self.EvapType == 'Micro-channel':
                self.Evaporator                 = MicroChannelEvaporatorClass()
                self.Evaporator.Fins            = MicroFinInputs()
            else:
                raise
        elif self.EvapSolver == 'Finite-Element':
            raise Exception('Discretized HX model not implemented')
        else:
            raise

        # -----------------------------------------------------------------------------
        if self.CondSolver == 'Moving-Boundary':
            if self.CondType == 'Fin-tube':
                self.Condenser                  = CondenserClass()
                self.Condenser.Fins             = FinInputs()
            elif self.CondType == 'Micro-channel':
                self.Condenser                  = MicroCondenserClass()
                self.Condenser.Fins             = MicroFinInputs()
            else:
                raise
        elif self.CondSolver == 'Finite-Element':
            raise Exception('Discretized HX model not implemented')
        else:
            raise

        # -----------------------------------------------------------------------------
        # Abstract State
        self.AS                     = CP.AbstractState(self.Backend, self.Ref)
        if hasattr(self, 'MassFrac'):
            self.AS.set_mass_fractions([self.MassFrac, 1 - self.MassFrac])
            self.AS.build_phase_envelope("dummy")
            self.AS.get_phase_envelope_data()
        elif hasattr(self, 'VoluFrac'):
            self.AS.set_volu_fractions([self.VoluFrac, 1 - self.VoluFrac])
            self.AS.build_phase_envelope("dummy")
            self.AS.get_phase_envelope_data()

        # Abstract State for SecLoopFluid
        self.AS_SLF                 = CP.AbstractState(self.Backend_SLF, self.SecLoopFluid)
        if hasattr(self, 'MassFrac_SLF'):
            self.AS_SLF.set_mass_fractions([self.MassFrac_SLF])
        elif hasattr(self, 'VoluFrac_SLF'):
            self.AS_SLF.set_volu_fractions([self.VoluFrac_SLF])

    def OutputList(self):
        """
            Return a list of parameters for this component for further output

            It is a list of tuples, and each tuple is formed of items:
                [0] Description of value
                [1] Units of value
                [2] The value itself
        """
        return [
            ('Charge',                  'kg',           self.Charge),
            ('Condenser Subcooling',    'K',            self.DT_sc),
            ('Primary Ref.',            '-',            self.Ref),
            ('Secondary Ref.',          '-',            self.SecLoopFluid),
            ('Imposed Variable',        '-',            self.ImposedVariable),
            ('IHX Type',                '-',            self.IHXType),
            ('COP',                     '-',            self.COP),
            ('COSP',                    '-',            self.COSP),
            ('Net Capacity',            'W',            self.CoolingCoil.Capacity),
            ('Net Power',               'W',            self.Power),
            ('EER',                     'BTU/(W-hr)',   self.EER),
            ('SHR',                     '-',            self.SHR),
            ('Condensation temp (dew)', 'K',            self.T_dew_cond),
            ('Evaporation temp (dew)',  'K',            self.T_dew_evap),
            ('Wavychannel coolant temp (dew)', 'K',     self.WavyChannel.T_out_g),
        ]

    def Calculate(self, DT_evap, DT_cond, T_in_WC):
        """
        Inputs are differences in temperature [K] between HX air inlet temperature
        and the dew temperature for the heat exchanger.

        Required Inputs:
            DT_evap:
                Difference in temperature [K] between cooling coil air inlet temperature and refrigerant dew temperature
            DT_cond:
                Difference in temperature [K] between condenser air inlet temperature and refrigerant dew temperature
            Tin_CC:
                Inlet "glycol" temperature to line set feeding cooling coil
        """
        if self.Verbosity > 1:
            print('Inputs: DT_evap %7.4f DT_cond %7.4f fT_IHX %7.4f' %(DT_evap, DT_cond, T_in_WC))

        # AbstractState
        AS                          = self.AS
        # AbstractState for SecLoopFluid
        AS_SLF                      = self.AS_SLF

        # Surrounding
        if hasattr(self, 'T_0') and hasattr(self, 'p_0'):
            T_0                     = self.T_0
            p_0                     = self.p_0
        else:
            T_0                     = 25+273.15                     # K
            p_0                     = 101325                        # Pa

        AS.update(CP.PT_INPUTS, p_0, T_0)
        h_0                         = AS.hmass()
        s_0                         = AS.smass()

        # Store the values to save on computation for later
        self.DT_evap                = DT_evap
        self.DT_cond                = DT_cond
        self.T_in_WC                = T_in_WC

        # If the user doesn't include the Mode, set it to Air Conditioning
        if not hasattr(self, 'Mode'):
            self.Mode               = 'AC'

        if self.Mode == 'AC':
            self.T_dew_cond         = self.Condenser.Fins.Air.T_db + DT_cond
            self.T_dew_evap         = self.T_bat - DT_evap
        elif self.Mode == 'HP':
            self.T_dew_cond         = T_in_WC + DT_cond
            self.T_dew_evap         = self.Evaporator.Fins.Air.T_db - DT_evap
        else:
            raise ValueError('Mode must be AC or HP')

        # -----------------------------------------------------------------------------
        # Evaporator and condenser saturation pressures
        AS.update(CP.QT_INPUTS, 1.0, self.T_dew_cond)
        p_sat_cond                  = AS.p()                        # [Pa]
        AS.update(CP.QT_INPUTS, 1.0, self.T_dew_evap)
        p_sat_evap                  = AS.p()                        # [Pa]

        print('p_evap = ', p_sat_evap/1e3, ' [kPa]; p_cond = ', p_sat_cond/1e3, ' [kPa]')

        # Evaporator and condenser bubble temperatures
        AS.update(CP.PQ_INPUTS, p_sat_evap, 0.0)
        self.T_bubble_evap          = AS.T()                        # [K]
        AS.update(CP.PQ_INPUTS, p_sat_cond, 0.0)
        self.T_bubble_cond          = AS.T()                        # [K]

        # ----------------------------------------------------------------------------------
        # Cycle solver for 'AC' mode; transfer heat to outdoors
        # ----------------------------------------------------------------------------------
        if self.Mode == 'AC':
            if not hasattr(self.Compressor, 'm_dot_r') or self.Compressor.m_dot_r < 0.00001:
                # The first run of model, run the compressor just so you can get a preliminary value
                # for the mass flow rate for the line set
                # -----------------------------------------------------------------------------
                params  = {                 # dictionary -> key:value, e.g. 'key':2345,
                    'p_in_r':               p_sat_evap,
                    'p_out_r':              p_sat_cond,
                    'T_in_r':               self.T_dew_evap + self.DT_sh,
                    'AS':                   AS,
                }
                self.Compressor.Update(**params)
                self.Compressor.Calculate()

            # ----- (1) ---------------------------------------------------------------------
            # Calculate inlet enthalpy; AS denotes the refrigerant
            # right after the IHX, the evaporator; point 1
            AS.update(CP.PT_INPUTS, p_sat_evap, self.T_dew_evap + self.DT_sh)
            h_in                            = AS.hmass()                # [J/kg]
            params  = {
                'T_in':                     self.T_dew_evap + self.DT_sh,
                'p_in':                     p_sat_evap,
                'h_in':                     h_in,
                'm_dot':                    self.Compressor.m_dot_r,
                'AS':                       AS,
                'Name':                     'Suction Line'  # vapor line
            }
            self.LineSetSuction.Update(**params)
            self.LineSetSuction.Calculate()

            # ------(1) -> (2)--------------------------------------------------------------------
            # point 1 - 2, the compressor
            params  = {                     # dictionary -> key:value, e.g. 'key':2345,
                'p_in_r':                   p_sat_evap + self.DP_low,
                'p_out_r':                  p_sat_cond - self.DP_high,
                'T_in_r':                   TrhoPhase_ph(self.AS, p_sat_evap, self.LineSetSuction.h_out, self.T_bubble_evap, self.T_dew_evap)[0],
                'AS':                       AS
            }
            self.Compressor.Update(**params)
            self.Compressor.Calculate()

            # ------- (2) -------------------------------------------------------------------
            params  = {
                'T_in':                     self.Compressor.T_out_r,
                'p_in':                     p_sat_cond,
                'h_in':                     self.Compressor.h_out_r,
                'm_dot':                    self.Compressor.m_dot_r,
                'AS':                       AS,
                'Name':                     'Discharge Line'
            }
            self.LineSetDischarge.Update(**params)
            self.LineSetDischarge.Calculate()

            # -------(2) -> (3)-------------------------------------------------------------------
            params  = {
                'm_dot_r':                  self.Compressor.m_dot_r,
                'T_in_r':                   self.Compressor.T_out_r,
                'p_sat_r':                  p_sat_cond,
                'AS':                       AS,
                # 'Oil': self.Oil,
                # 'shell_pressure':self.shell_pressure,
            }
            self.Condenser.Update(**params)
            self.Condenser.Calculate()

            # ------- (3) ------------------------------------------------------------------
            params  = {
                'T_in':                     self.Condenser.T_out_r,
                'p_in':                     p_sat_cond,
                'h_in':                     self.Condenser.h_out_r,
                'm_dot':                    self.Compressor.m_dot_r,
                'AS':                       AS,
                'Name':                     'Liquid Line'
            }
            self.LineSetLiquid.Update(**params)
            self.LineSetLiquid.Calculate()

            # --------(3) -> (4)-------------------------------------------------------------------------
            params  = {
                'p_in_r':                   p_sat_cond,
                'h_in_r':                   self.LineSetLiquid.h_out,
                'T_sup':                    self.DT_sh,
                'p_out_r':                  p_sat_evap,
                'AS':                       AS,
            }
            self.ExpDev.Update(**params)
            self.ExpDev.Calculate()

            # -----------------------------------------------------------------------------
            # -----------------------------------------------------------------------------
            # --------(8) -> (5)-----------------------------------------------------------
            # Inlet enthalpy to LineSetSupply
            AS_SLF.update(CP.PT_INPUTS, self.Pump.p_in_g, T_in_WC)
            h_in_LineSetSupply              = AS_SLF.hmass()            # [J/kg]
            params  = {
                'T_in':                     T_in_WC,
                'm_dot':                    self.Pump.m_dot_g,
                'h_in':                     h_in_LineSetSupply,
                'p_in':                     self.Pump.p_in_g,
                'AS':                       AS_SLF,
                'Name':                     'Supply Line'
            }
            self.LineSetSupply.Update(**params)
            self.LineSetSupply.Calculate()

            # -----------------------------------------------------------------------------
            # Now run WavyChannel to predict inlet glycol temperature to IHX
            params  = {
                'm_dot_g':                  self.Pump.m_dot_g,
                'T_in_g':                   self.LineSetSupply.T_out,
                'p_in_g':                   self.WavyChannel.p_in_g,
                'AS_g':                     AS_SLF
            }
            self.WavyChannel.Update(**params)
            self.WavyChannel.Calculate()

            # -----------------------------------------------------------------------------
            # Inlet enthalpy to LineSetReturn
            AS_SLF.update(CP.PT_INPUTS, self.Pump.p_in_g, self.WavyChannel.T_out_g)
            h_in_LineSetReturn              = AS_SLF.hmass()            # [J/kg]
            params  = {
                'T_in':                     self.WavyChannel.T_out_g,
                'm_dot':                    self.Pump.m_dot_g,
                'h_in':                     h_in_LineSetReturn,
                'p_in':                     self.WavyChannel.p_in_g - self.WavyChannel.DP_g,
                'AS':                       AS_SLF,
                'Name':                     'Return Line'
            }
            self.LineSetReturn.Update(**params)
            self.LineSetReturn.Calculate()

            # -----------------------------------------------------------------------------
            if self.IHXType == 'Coaxial':
                params  = {
                    'm_dot_g':              self.Pump.m_dot_g,
                    'T_in_g':               self.LineSetReturn.T_out_g,
                    'p_in_g':               self.Pump.p_in_g,
                    'AS_g':                 AS_SLF,
                    'p_in_r':               p_sat_evap,
                    'h_in_r':               self.LineSetLiquid.h_out,
                    'AS_r':                 AS,
                    'm_dot_r':              self.Compressor.m_dot_r,
                }
                self.CoaxialIHX.Update(**params)
                self.CoaxialIHX.Calculate()

                self.IHX.Charge_r           = self.CoaxialIHX.Charge_r
                self.IHX.Q                  = self.CoaxialIHX.Q
                self.IHX.T_out_g            = self.CoaxialIHX.T_out_g
                self.IHX.DP_g               = self.CoaxialIHX.DP_g
                self.IHX.h_out_r            = self.CoaxialIHX.h_out_r
                self.IHX.T_out_r            = self.CoaxialIHX.T_out_r
                self.IHX.DP_r               = self.CoaxialIHX.DP_r
                if hasattr(self, 'PHEIHX'):
                    del self.PHEIHX

            # -----------------------------------------------------------------------------
            elif self.IHXType == 'PHE':
                # Inlet enthalpy to PHEIHX
                AS_SLF.update(CP.PT_INPUTS, self.PHEIHX.p_in_h, self.LineSetReturn.T_out)
                h_in_PHEIHX                 = AS_SLF.hmass()            # [J/kg]
                params  = {
                    'm_dot_h':              self.Pump.m_dot_g,
                    'h_in_h':               h_in_PHEIHX,
                    'p_in_h':               self.PHEIHX.p_in_h,
                    'AS_h':                 AS_SLF,
                    'm_dot_c':              self.Compressor.m_dot_r,
                    'p_in_c':               p_sat_evap,
                    'h_in_c':               self.LineSetLiquid.h_out,
                    'AS_c':                 AS
                }
                self.PHEIHX.Update(**params)
                self.PHEIHX.Calculate()

                self.IHX.Charge_r           = self.PHEIHX.Charge_c
                self.IHX.Q                  = self.PHEIHX.Q
                self.IHX.T_out_g            = self.PHEIHX.T_out_h
                self.IHX.DP_g               = self.PHEIHX.DP_h
                self.IHX.DP_r               = self.PHEIHX.DP_c
                self.IHX.h_out_r            = self.PHEIHX.h_out_c
                self.IHX.T_out_r            = self.PHEIHX.T_out_c
                if hasattr(self, 'CoaxialIHX'):
                    del self.CoaxialIHX

            # -----------------------------------------------------------------------------
            params  = {
                'DP_g':     self.IHX.DP_g + self.WavyChannel.DP_g + self.LineSetSupply.DP + self.LineSetReturn.DP,
                'T_in_g':                   self.LineSetReturn.T_out,
                'p_in_g':                   self.Pump.p_in_g,
                'AS_g':                     AS_SLF,
            }
            self.Pump.Update(**params)
            self.Pump.Calculate()

            # -----------------------------------------------------------------------------
            # Energy and Mass conservation
            self.EnergyBalance      = self.Compressor.CycleEnergyIn + self.Condenser.Q + self.IHX.Q + self.LineSetSuction.Q \
                                      + self.LineSetDischarge.Q + self.LineSetLiquid.Q

            Charge                  = self.Compressor.Charge + self.Condenser.Charge + self.IHX.Charge_r + self.LineSetSuction.Charge \
                                      + self.LineSetDischarge.Charge + self.LineSetLiquid.Charge

            self.Charge_noCorrect           = Charge                            # Charge with no correction
            # Correct the charge with Bo Shen's two point regression method
            if self.ImposedVariable == 'Charge':
                C                           = self.C                            # [kg]
                if self.ChargeMethod == 'One-point':
                    self.Charge_correct     = Charge + C                        # one-point charge correction
                elif self.ChargeMethod == 'Two-point':
                    K                       = self.K                            # [kg]
                    w_ref                   = self.w_ref                        # [-]
                    w_liq                   = self.Condenser.w_subcool
                    self.delta_charge       = C + K * (w_liq - w_ref)
                    self.Charge_correct     = Charge + self.delta_charge        # two-point charge correction
                self.Charge                 = self.Charge_correct
            else:
                self.Charge                 = Charge

            # -----------------------------------------------------------------------------
            # Calculate properties
            AS.update(CP.QT_INPUTS, 0.0, self.T_bubble_cond)
            h_L                             = AS.hmass()                        # [J/kg]
            cp_L                            = AS.cpmass()                       # [J/kg-K]
            AS.update(CP.PT_INPUTS, p_sat_cond, self.T_bubble_cond - self.DT_sc_target)
            h_target                        = AS.hmass()                        # [J/kg]
            self.DT_sc                      = (h_L - self.Condenser.h_out_r) / cp_L
            delta_H_sc                      = self.Compressor.m_dot_r * (h_L - h_target)

            # -----------------------------------------------------------------------------
            resid                           = np.zeros((3))
            resid[0]                        = self.Compressor.m_dot_r * (self.IHX.h_out_r - self.LineSetSuction.h_in)

            if self.ImposedVariable == 'Subcooling':
                resid[1]                    = self.Condenser.DT_sc - self.DT_sc_target
            elif self.ImposedVariable == 'Charge':
                resid[1]                    = self.Charge - self.Charge_target
            # resid[2]=self.IHX.Q-self.CoolingCoil.Q+self.Pump.W

            self.resid_SL    = self.IHX.Q - self.WavyChannel.Q_wavychannel + self.Pump.W + self.LineSetSupply.Q + self.LineSetReturn.Q
            resid[2]                        = self.resid_SL

            if self.Verbosity > 7:
                print('W_comp % 12.6e Q_cond: % 12.6e Q_PHE %10.4f ' %(self.Compressor.W, self.Condenser.Q, self.IHX.Q))
            if self.Verbosity > 1:
                print('Q_res % 12.6e Resid2: % 12.6e ResSL %10.4f  Charge %10.4f SC: %8.4f'
                      % (resid[0], resid[1], self.resid_SL, self.Charge, self.Condenser.DT_sc))

            # -----------------------------------------------------------------------------
            # self.Capacity                   = self.WavyChannel.Q_wavychannel
            self.Capacity                   = self.PHEIHX.Q
            self.COP                        = self.WavyChannel.Q_wavychannel / self.Compressor.W
            self.COSP       = self.Capacity / (self.Compressor.W + self.Pump.W + self.Condenser.Fins.Air.FanPower)
            self.Power      = self.Compressor.W + self.Pump.W + self.Condenser.Fins.Air.FanPower
            self.EER                        = W2BTUh(self.Capacity) / self.Power                                        # [BTU/(W-hr)]
            self.DP_high_Model              = self.Condenser.DP_r + self.LineSetDischarge.DP + self.LineSetLiquid.DP    # [Pa]
            self.DP_low_Model               = self.IHX.DP_r + self.LineSetSuction.DP                                    # [Pa]

        # ----------------------------------------------------------------------------------
        # Cycle solver for 'HP' mode; transfer heat to indoors
        # ----------------------------------------------------------------------------------
        elif self.Mode == 'HP':
            if p_sat_evap + self.DP_low < 0:
                raise ValueError('Compressor inlet pressure less than zero ['+str(p_sat_evap + self.DP_low) +
                                 ' Pa] - is low side pressure drop too high?')

            if not hasattr(self.Compressor, 'm_dot_r') or self.Compressor.m_dot_r < 0.00001:
                # The first run of model, run the compressor just so you can get a preliminary value
                # for the mass flow rate for the line set
                # ----------------------------------------------------------------------------------
                params  = {                     # dictionary -> key:value, e.g. 'key':2345,
                    'p_in_r':                   p_sat_evap,
                    'p_out_r':                  p_sat_cond,
                    'T_in_r':                   self.T_dew_evap+self.DT_sh,
                    'AS':                       AS,
                }
                self.Compressor.Update(**params)
                self.Compressor.Calculate()

            # ----------------------------------------------------------------------------------
            # Calculate inlet enthalpy
            AS.update(CP.PT_INPUTS, p_sat_evap, self.T_dew_evap + self.DT_sh)
            h_in                                = AS.hmass()                        # [J/kg]
            params  = {
                'T_in':                         self.T_dew_evap+self.DT_sh,
                'p_in':                         p_sat_evap,
                'h_in':                         h_in,
                'm_dot':                        self.Compressor.m_dot_r,
                'AS':                           AS,
                'Name':                         'Suction Line'
            }
            self.LineSetSuction.Update(**params)
            self.LineSetSuction.Calculate()

            # ----------------------------------------------------------------------------------
            params  = {                         # dictionary -> key:value, e.g. 'key':2345,
                'p_in_r':                       p_sat_evap + self.DP_low,
                'p_out_r':                      p_sat_cond - self.DP_high,
                'T_in_r': TrhoPhase_ph(self.AS, p_sat_evap, self.LineSetSuction.h_out, self.T_bubble_evap, self.T_dew_evap)[0],
                'AS':                           AS
            }
            self.Compressor.Update(**params)
            self.Compressor.Calculate()

            # ----------------------------------------------------------------------------------
            params  = {
                'T_in':                         self.Compressor.T_out_r,
                'p_in':                         p_sat_cond,
                'h_in':                         self.Compressor.h_out_r,
                'm_dot':                        self.Compressor.m_dot_r,
                'AS':                           AS,
                'Name':                         'Discharge Line'
            }
            self.LineSetDischarge.Update(**params)
            self.LineSetDischarge.Calculate()

            # ----------------------------------------------------------------------------------
            # Inlet enthalpy to LineSetSupply
            AS_SLF.update(CP.PT_INPUTS, self.Pump.p_in_g, T_in_WC)
            h_in_LineSetSupply                  = AS_SLF.hmass()                # [J/kg]
            params  = {
                'T_in':                         T_in_WC,
                'm_dot':                        self.Pump.m_dot_g,
                'h_in':                         h_in_LineSetSupply,
                'AS':                           AS_SLF,
                'Name':                         'Supply Line'
            }
            self.LineSetSupply.Update(**params)
            self.LineSetSupply.Calculate()

            # ----------------------------------------------------------------------------------
            # Now run WavyChannel to predict inlet glycol temperature to IHX
            params  = {
                'm_dot_g':                      self.Pump.m_dot_g,
                'T_in_g':                       self.LineSetSupply.T_out,
                'p_in_g':                       self.PHEIHX.p_in_c,
                'AS_g':                         AS_SLF
            }
            self.WavyChannel.Update(**params)
            self.WavyChannel.Calculate()

            # ----------------------------------------------------------------------------------
            # Inlet enthalpy to LineSetReturn
            AS_SLF.update(CP.PT_INPUTS, self.Pump.p_in_g, self.WavyChannel.T_out_g)
            h_in_LineSetReturn                  = AS_SLF.hmass()                # [J/kg]
            params  = {
                'T_in':                         self.WavyChannel.T_out_g,
                'm_dot':                        self.Pump.m_dot_g,
                'h_in':                         h_in_LineSetReturn,
                'AS':                           AS_SLF,
                'Name':                         'Return Line'
            }
            self.LineSetReturn.Update(**params)
            self.LineSetReturn.Calculate()

            # ----------------------------------------------------------------------------------
            # Inlet enthalpy to PHEIHX
            AS_SLF.update(CP.PT_INPUTS, self.PHEIHX.p_in_c, self.LineSetReturn.T_out)
            h_in_PHEIHX                         = AS_SLF.hmass()                # [J/kg]
            params  = {
                'm_dot_h':                      self.Compressor.m_dot_r,
                'h_in_h':                       self.LineSetDischarge.h_out,
                'p_in_h':                       p_sat_cond,
                'AS_h':                         AS,
                'm_dot_c':                      self.Pump.m_dot_g,
                'h_in_c':                       h_in_PHEIHX,
                'p_in_c':                       self.Pump.p_in_g,
                'AS_c':                         AS_SLF,

            }
            self.PHEIHX.Update(**params)
            self.PHEIHX.Calculate()

            # ----------------------------------------------------------------------------------
            params  = {
                'T_in':                         self.PHEIHX.T_out_h,
                'p_in':                         p_sat_cond,
                'h_in':                         self.PHEIHX.h_out_h,
                'm_dot':                        self.Compressor.m_dot_r,
                'AS':                           AS,
                'Name':                         'Liquid Line'
            }
            self.LineSetLiquid.Update(**params)
            self.LineSetLiquid.Calculate()

            # ----------------------------------------------------------------------------------
            params  = {
                'm_dot_r':                      self.Compressor.m_dot_r,
                'p_sat_r':                      p_sat_evap,
                'h_in_r':                       self.LineSetLiquid.h_out,
                'AS':                           AS,
            }
            self.Evaporator.Update(**params)
            self.Evaporator.Calculate()

            # ----------------------------------------------------------------------------------
            params  = {
                'DP_g':      self.PHEIHX.DP_c + self.WavyChannel.DP_g + self.LineSetSupply.DP + self.LineSetReturn.DP,
                'T_in_g':                       self.WavyChannel.T_out_g,
                'p_in_g':                       self.Pump.p_in_g,
                'AS_g':                         AS_SLF
            }
            self.Pump.Update(**params)
            self.Pump.Calculate()

            # ----------------------------------------------------------------------------------
            # No Energy conservation???????
            Charge          = self.Compressor.Charge + self.Evaporator.Charge + self.PHEIHX.Charge_h + self.LineSetSuction.Charge \
                              + self.LineSetDischarge.Charge + self.LineSetLiquid.Charge

            self.Charge_noCorrect               = Charge                            # Charge with no correction

            # Correct the charge with Bo Shen's two point regression method
            if self.ImposedVariable == 'Charge':
                C                               = self.C                            # [kg]
                if self.ChargeMethod == 'One-point':
                    self.Charge_correct         = Charge + C                        # one-point charge correction
                elif self.ChargeMethod == 'Two-point':
                    K                           = self.K                            # [kg]
                    w_ref                       = self.w_ref                        # [-]
                    w_liq                       = self.Condenser.w_subcool
                    self.delta_charge           = C + K * (w_liq - w_ref)
                    self.Charge_correct         = Charge + self.delta_charge        # two-point charge correction
                self.Charge                     = self.Charge_correct
            else:
                self.Charge                     = Charge

            # ----------------------------------------------------------------------------------
            # Calculate properties
            AS.update(CP.QT_INPUTS, 0.0, self.T_bubble_cond)
            h_L                                 = AS.hmass()                        # [J/kg]
            cp_L                                = AS.cpmass()                       # [J/kg-K]
            AS.update(CP.PT_INPUTS, p_sat_cond, self.T_bubble_cond - self.DT_sc_target)
            h_target                            = AS.hmass()                        # [J/kg]
            # Calculate an effective subcooling amount by delta_h/cp_satL
            # Can be positive or negative (negative is quality at outlet
            # (PropsSI('H','T',self.Tbubble_cond,'Q',0,self.Ref)-self.PHEIHX.hout_h)/
            # (PropsSI('C','T',self.Tbubble_cond,'Q',0,self.Ref)) #*1000 #*1000
            self.DT_sc                          = self.PHEIHX.DT_sc_h
            delta_H_sc                          = self.Compressor.m_dot_r * (h_L - h_target)

            # ----------------------------------------------------------------------------------
            resid                   = np.zeros((3))
            resid[0]                = self.Compressor.m_dot_r * (self.Evaporator.h_out_r - self.LineSetSuction.h_in)

            if self.ImposedVariable == 'Subcooling':
                resid[1]            = self.DT_sc - self.DT_sc_target
            elif self.ImposedVariable == 'Charge':
                resid[1]            = self.Charge - self.Charge_target
            self.residSL            = self.PHEIHX.Q + self.WavyChannel.Q_wavychannel + self.Pump.W + self.LineSetSupply.Q + self.LineSetReturn.Q
            resid[2]                = self.residSL

            if self.Verbosity > 1:
                print('Q_res % 12.6e Resid2: % 12.6e ResSL %10.4f Charge %10.4f SC: %8.4f'
                      %(resid[0], resid[1], self.residSL, self.Charge, self.DT_sc))

            # ----------------------------------------------------------------------------------
            self.Capacity           = -self.WavyChannel.Q_wavychannel + self.Condenser.Fins.Air.FanPower
            self.COP                = -self.WavyChannel.Q_wavychannel / self.Compressor.W
            self.COSP   = self.Capacity / (self.Compressor.W + self.Pump.W + self.Evaporator.Fins.Air.FanPower)
            self.Power  = self.Compressor.W + self.Pump.W + self.Evaporator.Fins.Air.FanPower
            self.EER                = W2BTUh(self.Capacity) / self.Power                                    # [BTU/(W-hr)]
            self.DP_high_Model      = self.PHEIHX.DP_h + self.LineSetDischarge.DP + self.LineSetLiquid.DP   # [Pa]
            self.DP_low_Model       = self.Evaporator.DP_r + self.LineSetSuction.DP                         # [Pa]

        self.DT_evap                = DT_evap
        self.DT_cond                = DT_cond

        return resid

    def PreconditionedSolve(self, PrecondValues=None):
        """
        PrecondValues = dictionary of values DT_evap, DT_cond and Tin_CC
        """

        def OBJECTIVE(x):
            """
            Takes the place of a lambda function since lambda functions do not bubble error properly
            """
            return self.Calculate(x[0], x[1], x[2])

        def OBJECTIVE2(x, T_in):
            """
            Takes the place of a lambda function since lambda functions do not bubble error properly
            """
            return self.Calculate(x[0], x[1], T_in)

        def OBJECTIVE_SL(T_in_WC):
            """
            Objective function for the inner loop of the vapor compression system

            Using the MultiDimNewtRaph function will re-evaluate the Jacobian at
            every step.  Slower, but more robust since the solution surfaces aren't
            smooth enough

            Note: This function is not currently used!
            """
            x   = MultiDimNewtRaph(OBJECTIVE2, [self.DT_evap, self.DT_cond], args=(T_in_WC,))

            # Update the guess values for Delta Ts starting
            # at the third step (after at least one update away
            # from the boundaries)
            if self.OBJ_SL_counter >= 0:
                self.DT_evap                    = x[0]
                self.DT_cond                    = x[1]
                pass
            self.OBJ_SL_counter                 += 1
            return self.residSL

        def PrintDPs():
            print('DP_LP :: Input:', self.DP_low,  'Pa / Model calc:', self.DP_low_Model,  'Pa')
            print('DP_HP :: Input:', self.DP_high, 'Pa / Model calc:', self.DP_high_Model, 'Pa')

        # Some variables need to be initialized
        self.DP_low                             = 0     # The actual low-side pressure drop to be used in Pa
        self.DP_high                            = 0     # The actual low-side pressure drop to be used in Pa
        self.OBJ_SL_counter                     = 0

        # Run the preconditioner to get guess values for the temperatures
        if PrecondValues is None:
            self.DT_evap, self.DT_cond, T_in_WC = BTMS_SecondaryLoopPreconditioner(self)    # The preconditioner has been changed
        else:
            self.DT_evap                        = PrecondValues['DT_evap']
            self.DT_cond                        = PrecondValues['DT_cond']
            T_in_WC                             = PrecondValues['T_in_WC']

        # Remove the other, non-used IHX class if found
        if self.IHXType == 'PHE':
            if hasattr(self, 'CoaxialIHX'):
                del self.CoaxialIHX
        else:
            if hasattr(self, 'PHEIHX'):
                del self.PHEIHX

        # Remove the condenser if in heating mode and condenser found
        if self.Mode == 'HP':
            if hasattr(self, 'Condenser'):
                del self.Condenser

        # ----------------------------------------------------------------------------------
        iter                    = 1
        max_error_DP            = 999
        # Outer loop with a more relaxed convergence criterion
        while max_error_DP > 0.5:
            iter_inner          = 1
            # Inner loop to determine pressure drop for high and low sides
            while max_error_DP > 0.05 and iter_inner < 10:

                # Run to calculate the pressure drop as starting point
                OBJECTIVE([self.DT_evap, self.DT_cond, T_in_WC])

                # Calculate the max error
                max_error_DP    = max([abs(self.DP_low_Model - self.DP_low), abs(self.DP_high_Model - self.DP_high)])

                if self.Verbosity > 0:
                    PrintDPs()
                    print('Max pressure drop error [inner loop] is', max_error_DP, 'Pa')

                # Update the pressure drop terms
                self.DP_low     = self.DP_low_Model         # /1000
                self.DP_high    = self.DP_high_Model        # /1000

                iter_inner      += 1

            if self.Verbosity > 0:
                print("Done with the inner loop on pressure drop")

            # Use Newton-Raphson solver
            (self.DT_evap, self.DT_cond, T_in_WC)   = MultiDimNewtRaph(OBJECTIVE, [self.DT_evap, self.DT_cond, T_in_WC], dx=0.1)

            # Calculate the error
            max_error_DP    = max([abs(self.DP_low_Model - self.DP_low), abs(self.DP_high_Model - self.DP_high)])

            if self.Verbosity > 0:
                PrintDPs()
                print('Max pressure drop error [outer loop] is', max_error_DP, 'Pa')

        if self.Verbosity > 1:
            print('Capacity: ',             self.Capacity)
            print('COP: ',                  self.COP)
            print('COP (w/ both fans): ',   self.COSP)

        print('-------------------------------------')
        print('     Simulation Completed            ')
        print('-------------------------------------')

        return
