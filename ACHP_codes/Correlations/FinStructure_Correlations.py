from __future__                                 import division, print_function, absolute_import

import numpy as np
from CoolProp.CoolProp                          import HAPropsSI
from math                                       import sqrt,pi,tanh,exp,cos,log
from ACHP_codes.ACHP_Tools.ACHP_Tools           import ValidateFields


class FinInputs:
    """
    Empty Class for fin data
    """
    def __init__(self):
        # Don't do anything
        self.Tubes      = TubesVals()
        self.Air        = AirVals()
        self.Fins       = FinsVals()

    def __repr__(self):
        string          = "Tubes::\n"
        for field in self.Tubes.__dict__.keys():
            string      += field + ": " + repr(self.Tubes.__dict__[field]) + "\n"

        string          += "Fins::\n"
        for field in self.Fins.__dict__.keys():
            string      += field + ": " + repr(self.Fins.__dict__[field]) + "\n"

        string          += "Air::\n"
        for field in self.Air.__dict__.keys():
            string      += field + ": " + repr(self.Air.__dict__[field]) + "\n"

        return string

    def Validate(self):
        """
        Check that all required fields are included, and no extra fields not listed here are added
        This is quite strict, but for the best to avoid typos

        =========   =============================================================
        Variable    Type & Description
        =========   =============================================================
        d           dict of values that are part of structure
        reqFields   list of tuples in the form (fieldname, typepointer, min, max)
        optFields   list of other fieldnames
        =========   =============================================================
        """
        # for air check
        reqFields = [
            ('RH',          float,  0,          1),
            ('T_db',        float,  -80+273.15, 200+273.15),
            ('FanPower',    float,  0,          4000),
            ('p',           float,  10,         10000000),
            ('V_dot_ha',    float,  0.001,      10)]
        optFields       = ['RH_mean', 'T_mean']
        d               = dict(self.Air.__dict__)   # The current values
        ValidateFields(d, reqFields, optFields)

        # for fin check
        reqFields = [
            ('FPI',         float,  0.1,        100),
            ('Pd',          float,  0.000001,   0.1),
            ('xf',          float,  0.000001,   0.1),
            ('t',           float,  0.00001,    0.01),
            ('k_fin',       float,  0.01,       10000)]
        optFields       = None
        d               = dict(self.Fins.__dict__)  # The current values
        ValidateFields(d, reqFields, optFields)

        # for fin check
        reqFields = [
            ('N_Tubes_per_bank',    int,    0.1,    100),
            ('N_circuits',          int,    1,      50),
            ('N_bank',              float,  1,      50),
            ('L_tube',              float,  0.001,  10),
            ('ID',                  float,  0.0001, 1),
            ('OD',                  float,  0.0001, 1),
            ('Pl',                  float,  0.0001, 1),
            ('Pt',                  float,  0.0001, 1),
            ('k_w',                 float,  0.01,   10000)]
        optFields       = None
        d               = dict(self.Tubes.__dict__)  # The current values
        ValidateFields(d, reqFields, optFields)


class FinInputs_Iter:
    """
    Empty Class for fin data
    """
    def __init__(self):
        # Don't do anything
        self.Tubes      = TubesVals()
        self.Air        = AirVals()
        self.Fins       = FinsVals()

    def __repr__(self):
        string          = "Tubes::\n"
        for field in self.Tubes.__dict__.keys():
            string      += field + ": " + repr(self.Tubes.__dict__[field]) + "\n"

        string          += "Fins::\n"
        for field in self.Fins.__dict__.keys():
            string      += field + ": " + repr(self.Fins.__dict__[field]) + "\n"

        string          += "Air::\n"
        for field in self.Air.__dict__.keys():
            string      += field + ": " + repr(self.Air.__dict__[field]) + "\n"

        return string

    def Validate(self):
        """
        Check that all required fields are included, and no extra fields not listed here are added
        This is quite strict, but for the best to avoid typos

        =========   =============================================================
        Variable    Type & Description
        =========   =============================================================
        d           dict of values that are part of structure
        reqFields   list of tuples in the form (fieldname, typepointer, min, max)
        optFields   list of other fieldnames
        =========   =============================================================
        """
        # for air check
        reqFields = [
            ('RH',          float,  0,          1),
            ('FanPower',    float,  0,          4000),
            ('p',           float,  10,         10000000),
            ('V_dot_ha',    float,  0.001,      10)]
        optFields       = ['RH_mean', 'T_mean']
        d               = dict(self.Air.__dict__)   # The current values
        ValidateFields(d, reqFields, optFields)

        # for fin check
        reqFields = [
            ('FPI',         float,  0.1,        100),
            ('Pd',          float,  0.000001,   0.1),
            ('xf',          float,  0.000001,   0.1),
            ('t',           float,  0.00001,    0.01),
            ('k_fin',       float,  0.01,       10000)]
        optFields       = None
        d               = dict(self.Fins.__dict__)  # The current values
        ValidateFields(d, reqFields, optFields)

        # for fin check
        reqFields = [
            ('N_Tubes_per_bank',    int,    0.1,    100),
            ('N_circuits_per_bank', int,    1,      50),
            ('L_tube',              float,  0.001,  10),
            ('ID',                  float,  0.0001, 1),
            ('OD',                  float,  0.0001, 1),
            ('Pl',                  float,  0.0001, 1),
            ('Pt',                  float,  0.0001, 1),
            ('k_w',                 float,  0.01,   10000)]
        optFields       = None
        d               = dict(self.Tubes.__dict__)  # The current values
        ValidateFields(d, reqFields, optFields)


# --------------------------------------------------------------------------
def WavyLouveredFins(Inputs):
    """
    # Correlations from:
    # Chi-Chuan Wang and Yu-Min Tsai and Ding-Chong Lu, 1998, "Comprehensive
    # Study of Convex-Louver and Wavy Fin-and-Tube Heat Exchangers", Journal
    # of Thermophysics and Heat Transfer

    || -    xf    ->
    ^              ==                          ==
    |           ==  |  ==                  ==
    Pd       ==     |     ==            ==
    |     ==        |        ==     ==
    =  ==           s             ==
                    |
                    |
                    |
                   ==                        ==
                ==     ==                 ==
             ==           ==           ==
          ==                 ==     ==
       ==                        ==

     t: thickness of fin plate
     Pf: fin pitch (centerline-to-centerline distance between fins)
     Pd: indentation for waviness (not including fin thickness)
     s: fin spacing (free space between fins) = Pf-t



                 |--       Pl      -|
                ___                 |
              /     \               |
       =     |       |              |
       |      \ ___ /               |
       |                            |
       |                           ___
       |                         /     \      |
      Pt                        |       |     D
       |                         \ ___ /      |
       |
       |        ___
       |      /     \
       =     |       |
              \ ___ /
    """

    N_tubes_bank        = Inputs.Tubes.N_Tubes_per_bank         # tubes per bank
    N_bank              = Inputs.Tubes.N_bank                   # Number of banks
    L_tube              = Inputs.Tubes.L_tube                   # length of a single tube
    D                   = Inputs.Tubes.OD                       # Outer diameter of tube
    Pl                  = Inputs.Tubes.Pl                       # Horizontal spacing between banks (center to center)
    Pt                  = Inputs.Tubes.Pt                       # Vertical spacing between tubes in a bank (center to center)

    FPI                 = Inputs.Fins.FPI
    pd                  = Inputs.Fins.Pd
    xf                  = Inputs.Fins.xf
    t                   = Inputs.Fins.t
    k_fin               = Inputs.Fins.k_fin

    V_dot_ha            = Inputs.Air.V_dot_ha
    p                   = Inputs.Air.p

    if hasattr(Inputs, 'T_mean'):
        Temp            = Inputs.Air.T_mean
    else:
        Temp            = Inputs.Air.T_db

    if hasattr(Inputs, 'RH_mean'):
        RH_in           = Inputs.Air.RH_mean
    else:
        RH_in           = Inputs.Air.RH

    # Check that cs_cp is defined, if so, set it to the value passed in
    if (hasattr(Inputs, 'cs_cp') and Inputs.cs_cp > 0) or (hasattr(Inputs, 'WetDry') and Inputs.WetDry == 'Wet'):
        isWet           = True
        cs_cp           = Inputs.Air.cs_cp
    else:
        isWet           = False
        cs_cp           = 1.0

    # Film temperature [K]
    T_film              = Temp
    # Fins per meter [1/m]
    FPM                 = FPI / 0.0254
    # Fin pitch (distance between centerlines of fins)
    pf = 1 / FPM
    # Spacing between fins
    s                   = 1 / FPM - t

    # Height of heat exchanger [m]
    Height              = Pt * (N_tubes_bank + 1)       # assuming that fin extends 1/2 pt above/below last tube in bundle
    # A_duct is the face area [m^2] equivalent to the duct cross-section
    A_duct              = Height * L_tube               # neglecting the additional height of the fins above/below the last tubes in the bundle
    # Number of fins in the tube sheet [-]
    N_fin               = L_tube * FPM
    # Secant of theta is the area enhancement factor [-]
    # It captures the increase in area due to the waviness of the fins
    sec_theta           = sqrt(xf*xf + pd*pd) / xf
    # Duct cross-sectional area that is not fin or tube [m^2]
    Ac                  = A_duct - t * N_fin * (Height - D * N_tubes_bank) - N_tubes_bank * D * L_tube
    # Total outer area of the tubes [m^2]
    A_tube              = N_tubes_bank * N_bank * pi * D * L_tube
    # Wetted Area of a single fin [m^2]
    A_1fin              = 2.0 * (Height * Pl * (N_bank + 1) * sec_theta - N_tubes_bank * N_bank * pi * D*D / 4)
                                                        # assuming that fin extends 1/2 pt in front/after last tube in bundle
    # Total wetted area of the fins [m^2]
    Af                  = N_fin * A_1fin
    # Total area including tube and fins [m^2]
    A = Af + N_tubes_bank * N_bank * pi * D * (L_tube - N_fin * t)

    # Evaluate the mass flow rate based on inlet conditions
    # To convert a parameter from per kg_{dry air} to per kg_{humid air}, divide by (1+W)
    W                   = HAPropsSI('W', 'T', Inputs.Air.T_db, 'P', p, 'R', Inputs.Air.RH)
    v_da                = HAPropsSI('V', 'T', Inputs.Air.T_db, 'P', p, 'W', W)
    h_da                = HAPropsSI('H', 'T', Inputs.Air.T_db, 'P', p, 'W', W)          # [J/kg]
    rho_ha              = 1 / v_da * (1 + W)                                            # [kg_ha/m^3]
    rho_da              = 1 / v_da                                                      # [kg_da/m^3]
    m_dot_ha            = V_dot_ha * rho_ha                                             # [kg_ha/s]
    m_dot_da            = V_dot_ha * rho_da                                             # [kg_da/s]

    umax                = m_dot_ha / (rho_ha * Ac)                                      # [m/s]

    # Use a forward difference to calculate cp from cp=dh/dT
    dT                  = 0.0001                                                        # [K]
    cp_da               = (HAPropsSI('H', 'T', Inputs.Air.T_db+dT, 'P', p, 'W', W) - h_da) / dT     # [J/kg_da/K]
    cp_ha               = cp_da / (1 + W)                                               # [J/kg_ha/K]

    # Transport properties of humid air from CoolProp
    mu_ha               = HAPropsSI('M', 'T', Inputs.Air.T_db, 'P', p, 'W', W)
    k_ha                = HAPropsSI('K', 'T', Inputs.Air.T_db, 'P', p, 'W', W)

    # Dimensionless Groups
    Pr                  = cp_ha * mu_ha / k_ha
    Re_D                = rho_ha * umax * D / mu_ha

    # Heat transfer
    j   = 16.06 * pow(Re_D, -1.02 * (pf / D) - 0.256) * pow(A / A_tube, -0.601) * pow(N_bank, -0.069) * pow(pf / D, 0.84) # Colburn j-Factor
    h_a = j * rho_ha * umax * cp_ha / pow(Pr, 2.0/3.0)                                  # air side mean heat transfer coefficient

    # air side pressure drop correlations
    if Re_D < 1e3:
        fa_total        = 0.264 * (0.105 + 0.708 * exp(-Re_D / 225.0)) * pow(Re_D, -0.637) * pow(A / A_tube, 0.263) * pow(pf / D, -0.317)
    else:
        fa_total        = 0.768 * (0.0494 + 0.142 * exp(-Re_D / 1180.0)) * pow(A / A_tube, 0.0195) * pow(pf / D, -0.121)

    # calcs needed for specific fin types
    r                   = D / 2
    X_D                 = sqrt(Pl*Pl + Pt*Pt / 4) / 2
    X_T                 = Pt / 2
    rf_r                = 1.27 * X_T / r * sqrt(X_D / X_T - 0.3)
    m                   = sqrt(2 * h_a * cs_cp / (k_fin * t))                           # cs_cp is the correction for heat/mass transfer

    # Using the circular fin correlation of Schmidt
    phi                 = (rf_r - 1) * (1 + 0.35 * log(rf_r))
    eta_f               = tanh(m * r * phi) / (m * r * phi)

    # Fin efficiency based on analysis in
    # "FIN EFFICIENCY CALCULATION IN ENHANCED FIN-AND-TUBE HEAT EXCHANGERS IN DRY CONDITIONS"
    # by Thomas PERROTIN, Denis CLODIC, International Congress of Refrigeration 2006
    # In the paper, there is no 0.1 in the cosine term, but if the cosine term is used without
    # the correction, the results are garbage for wet analysis
    # Using the offset fins correlation
    phi                 = (rf_r - 1) * (1 + (0.3 + pow(m * (rf_r * r - r) / 2.5, 1.5 - rf_r / 12.0) * (0.26 * pow(rf_r, 0.3) - 0.3)) * log(rf_r))

    # finned surface efficiency
    eta_f               = tanh(m * r * phi) / (m * r * phi) * cos(0.1 * m * r * phi)

    # overall surface efficiency
    eta_o               = 1 - Af / A * (1 - eta_f)

    G_c                 = m_dot_ha / Ac                                                 # air mass flux
    DeltaP_air          = A / Ac / rho_ha * G_c**2 / 2.0 * fa_total                     # airside pressure drop

    # write necessary values back into the given structure
    Inputs.A_a          = A
    Inputs.cp_da        = cp_da
    Inputs.cp_ha        = cp_ha

    if isWet == True:
        Inputs.eta_a_wet    = eta_o
    else:
        Inputs.eta_a        = eta_o

    Inputs.h_a              = h_a
    Inputs.m_dot_ha         = m_dot_ha
    Inputs.m_dot_da         = m_dot_da
    Inputs.f_a              = fa_total
    Inputs.dP_a             = DeltaP_air
    Inputs.Re               = Re_D



# --------------------------------------------------------------------------
def HerringboneFins(Inputs):
    # Source:
    # Empirical correlations for heat transfer and flow friction characteristics of herringbone wavy fin-and-tube heat exchangers
    # Chi-Chuan Wang, Young-Ming Hwang, Yur-Tsai Lin
    # International Journal of Refrigeration, 25, 2002, 637-680

    # --------------------------------------------------------
    # Air Properties
    p                   = Inputs.Air.p # air pressure
    # http://www.coolprop.org/fluid_properties/HumidAir.html; da means dry air
    W                   = HAPropsSI('W', 'T', Inputs.Air.T_db, 'P', p, 'R', Inputs.Air.RH)  # Humidity ratio, kg water/kg dry air

    # Transport properties of humid air from CoolProp; ha means humid air
    mu_ha               = HAPropsSI('M', 'T', Inputs.Air.T_db, 'P', p, 'W', W)
    k_ha                = HAPropsSI('K', 'T', Inputs.Air.T_db, 'P', p, 'W', W)

    # Evaluate the mass flow rate based on inlet conditions
    V_dot_ha            = Inputs.Air.V_dot_ha

    # To convert a parameter from per kg_{dry air} to per kg_{humid air}, divide by (1+W)
    W                   = HAPropsSI('W', 'T', Inputs.Air.T_db, 'P', p, 'R', Inputs.Air.RH)
    v_da                = HAPropsSI('V', 'T', Inputs.Air.T_db, 'P', p, 'W', W)
    h_da                = HAPropsSI('H', 'T', Inputs.Air.T_db, 'P', p, 'W', W)  # [J/kg]
    rho_ha              = 1 / v_da * (1 + W)    # [kg_ha/m^3]
    rho_da              = 1 / v_da              # [kg_da/m^3]
    m_dot_ha            = V_dot_ha * rho_ha     # [kg_ha/s]
    m_dot_da            = V_dot_ha * rho_da     # [kg_da/s]

    # Use a forward difference to calculate cp from cp=dh/dT
    dT                  = 0.0001                # [K]
    cp_da               = (HAPropsSI('H', 'T', Inputs.Air.T_db+dT, 'P', p, 'W', W) - h_da) / dT  # [J/kg_da/K]
    cp_ha               = cp_da / (1 + W)       # [J/kg_ha/K]

    # Check that cs_cp is defined, if so, set it to the value passed in
    if (hasattr(Inputs, 'cs_cp') and Inputs.cs_cp > 0) or (hasattr(Inputs, 'WetDry') and Inputs.WetDry == 'Wet'):
        isWet           = True
        cs_cp           = Inputs.Air.cs_cp      # c_s / c_p
    else:
        isWet           = False
        cs_cp           = 1.0

    # --------------------------------------------------------
    # Fin structure parameters
    # Dimensions and values used for both Reynolds number ranges:
    delta_f             = Inputs.Fins.t         # fin thickness(m)
    FPI                 = Inputs.Fins.FPI       # fins per inch
    FPM                 = FPI / 0.0254          # Fins per meter [1/m]
    pf                  = 1 / FPM               # Fin pitch (distance between center-lines of fins)
    F_s                 = 1 / FPM - delta_f     # Fin spacing(m)

    P_t                 = Inputs.Tubes.Pt       # transverse tube pitch (m), along the airflow direction
    P_L                 = Inputs.Tubes.Pl       # longitudinal tube pitch (m)

    D_o                 = Inputs.Tubes.OD       # Outer diameter of tube (m)
    D_c                 = D_o + 2 * delta_f     # fin collar outside diameter (m) (tube + 2*fin thickness)

    L_tube              = Inputs.Tubes.L_tube   # length of a single tube [m]
    N_fin               = L_tube * FPM          # Number of fins in the tube sheet [-]
    N_tubes_bank        = Inputs.Tubes.N_Tubes_per_bank     # tubes per bank
    Height              = P_t * (N_tubes_bank + 1)          # Height of heat exchanger [m] # assuming that fin extends 1/2 pt above/below last tube in bundle
    A_duct              = Height * L_tube                   # A_duct is the face area [m^2] - equivalent to the duct cross-section
                                                            # neglecting the additional height of the fins above/below the last tubes in the bundle
    Ac                  = A_duct - delta_f * N_fin * (Height - D_c * N_tubes_bank) - N_tubes_bank * D_c * L_tube
                                                            # Minimum duct cross-sectional area that is not fin or tube(-collar) [m^2]
    u_max               = m_dot_ha / (rho_ha * Ac)          # maximum airside velocity [m/s]
    Re_Dc               = rho_ha * u_max * D_c / mu_ha      # from reference #4 of Wang et al. Slightly different notation used in [4]

    P_d                 = Inputs.Fins.Pd                    # wave height
    X_f                 = Inputs.Fins.xf                    # projected fin length (m)
    tanTheta            = P_d / X_f                         # tangens of corrugation angle

    secTheta            = sqrt(X_f*X_f + P_d*P_d) / X_f     # !!used wavy louvered fins definition - there seems to be a bug in the paper !!
    # secTheta            = sqrt(4*X_f*X_f + P_d*P_d) / (2*X_f)  # todo check the correlations; sec = 1 / cos
    beta                = (pi * D_c**2) / (4.0 * P_t * P_L)
    oneMbeta            = 1.0 - beta                        # save one operation
    D_h                 = 2.0 * F_s * oneMbeta / (oneMbeta * secTheta + 2 * F_s * beta / D_c)   # hydraulic diameter (m)

    N                   = Inputs.Tubes.N_bank               # number of longitudinal tube rows (#Number of banks in ACHP-notation)

    if Re_Dc < 1000.0:
        # Heat transfer
        J1 = 0.0045-0.491*pow(Re_Dc,-0.0316-0.0171*log(N*tanTheta))*pow(P_L/P_t,-0.109*log(N*tanTheta))*pow(D_c/D_h,0.542+0.0471*N)*pow(F_s/D_c,0.984)*pow(F_s/P_t,-0.349)
        J2 = -2.72+6.84*tanTheta
        J3 = 2.66*tanTheta
        j  = 0.882*pow(Re_Dc,J1)*pow(D_c/D_h,J2)*pow(F_s/P_t,J3)*pow(F_s/D_c,-1.58)*pow(tanTheta,-0.2)

        # Friction
        F1 = -0.574-0.137*pow(log(Re_Dc)-5.26,0.245)*pow(P_t/D_c,-0.765)*pow(D_c/D_h,-0.243)*pow(F_s/D_h,-0.474)*pow(tanTheta,-0.217)*pow(N,0.035)
        F2 = -3.05*tanTheta
        F3 = -0.192*N
        F4 = -0.646*tanTheta
        f  = 4.37*pow(Re_Dc,F1)*pow(F_s/D_h,F2)*pow(P_L/P_t,F3)*pow(D_c/D_h,0.2054)*pow(N,F4)
    else:
        # Heat transfer
        j1 = -0.0545-0.0538*tanTheta-0.302*pow(N,-0.24)*pow(F_s/P_L,-1.3)*pow(P_L/P_t,0.379)*pow(P_L/D_h,-1.35)*pow(tanTheta,-0.256)
        j2 = -1.29*pow(P_L/P_t,1.77-9.43*tanTheta)*pow(D_c/D_h,0.229-1.43*tanTheta)*pow(N,-0.166-1.08*tanTheta)*pow(F_s/P_t,-0.174*log(0.5*N))
        j  = 0.0646*pow(Re_Dc,j1)*pow(D_c/D_h,j2)*pow(F_s/P_t,-1.03)*pow(P_L/D_c,0.432)*pow(tanTheta,-0.692)*pow(N,-0.737)

        # Friction
        f1 = -0.141*pow(F_s/P_L,0.0512)*pow(tanTheta,-0.472)*pow(P_L/P_t,0.35)*pow(P_t/D_h,0.449*tanTheta)*pow(N,-0.049+0.237*tanTheta)
        f2 = -0.562*pow(log(Re_Dc),-0.0923)*pow(N,0.013)
        f3 = 0.302*pow(Re_Dc,0.03)*pow(P_t/D_c,0.026)
        f4 = -0.306+3.63*tanTheta
        f  = 0.228*pow(Re_Dc,f1)*pow(tanTheta,f2)*pow(F_s/P_L,f3)*pow(P_L/D_c,f4)*pow(D_c/D_h,0.383)*pow(P_L/P_t,-0.247)

    Pr                  = cp_ha * mu_ha / k_ha
    h_a                 = j * rho_ha * u_max * cp_ha / pow(Pr,2.0/3.0)  # air side mean heat transfer coefficient using Colborn j-factor

    # calcs needed for specific fin types
    # additional parameters needed
    k_fin               = Inputs.Fins.k_fin
    N_bank              = Inputs.Tubes.N_bank              # Number of banks
    # secTheta=sqrt(X_f*X_f + P_d*P_d) / X_f  #secTheta : already calculated, no need to re-calculate it
    # Wetted Area of a single fin [m^2]
    A_1fin              = 2.0 * (Height * P_L * (N_bank+1) * secTheta  - N_tubes_bank * N_bank * pi * D_o * D_o / 4)
                                                            # assuming that fin extends 1/2 pt in front/after last tube in bundle
    # Total wetted area of the fins [m^2]
    Af                  = N_fin * A_1fin        # N_fin is the number of fin sheet in the transverse direction
    # Total wetted area including tube and fins [m^2]
    A                   = Af + N_tubes_bank * N_bank * pi * D_o * (L_tube - N_fin * delta_f)

    r                   = D_o / 2
    X_D                 = sqrt(P_L*P_L + P_t*P_t / 4) / 2
    X_T                 = P_t / 2
    rf_r                = 1.27 * X_T / r * sqrt(X_D / X_T - 0.3)
    m                   = sqrt(2 * h_a * cs_cp / (k_fin * delta_f))     # cs_cp is the correction for heat/mass transfer for a wetted surface

    # Using the circular fin correlation of Schmidt
    phi                 = (rf_r - 1) * (1 + 0.35 * log(rf_r))
    eta_f               = tanh(m * r * phi) / (m * r * phi)

    # Fin efficiency based on analysis in
    # "FIN EFFICIENCY CALCULATION IN ENHANCED FIN-AND-TUBE HEAT EXCHANGERS IN DRY CONDITIONS"
    # by Thomas PERROTIN, Denis CLODIC, International Congress of Refrigeration 2006
    # In the paper, there is no 0.1 in the cosine term, but if the cosine term is used without
    # the correction, the results are garbage for wet analysis
    # Using the offset fins correlation
    phi                 = (rf_r - 1) * (1 + (0.3+pow(m*(rf_r*r-r)/2.5, 1.5-rf_r/12.0)*(0.26*pow(rf_r,0.3)-0.3)) * log(rf_r))

    # finned surface efficiency
    eta_f               = tanh(m * r * phi) / (m * r * phi) * cos(0.1 * m * r * phi)

    # overall surface efficiency
    eta_o               = 1 - Af / A * (1 - eta_f)

    G_c                 = m_dot_ha / Ac                                 # air mass flux
    A_tube              = N_tubes_bank * N_bank * pi * D_o * L_tube     # Total outer area of the tubes [m^2]
    DeltaP_air          = (A / Ac) / rho_ha * G_c**2 / 2.0 * f          # airside pressure drop,
                            # todo: Is A/Ac = 4* L/d_H, where A is the total area; Ac is the minimum free flow area
                            # Yes, D_h = 4 * Ac * L / A, leading to A / Ac = 4 * L / D_h; thus f is the Fanning friction factor, not the Darcy one

    # write necessary values back into the given structure
    Inputs.A_a          = A
    Inputs.cp_da        = cp_da
    Inputs.cp_ha        = cp_ha

    if isWet == True:
        Inputs.eta_a_wet= eta_o
    else:
        Inputs.eta_a    = eta_o

    Inputs.h_a          = h_a
    Inputs.m_dot_ha     = m_dot_ha
    Inputs.m_dot_da     = m_dot_da
    Inputs.f_a          = f
    Inputs.dP_a         = DeltaP_air
    Inputs.Re           = Re_Dc

    return Height, secTheta, N_fin, delta_f


def HerringboneFins_Iter(Inputs, T_db, N_bank=1):
    # Source:
    # Empirical correlations for heat transfer and flow friction characteristics of herringbone wavy fin-and-tube heat exchangers
    # Chi-Chuan Wang, Young-Ming Hwang, Yur-Tsai Lin
    # International Journal of Refrigeration, 25, 2002, 637-680

    # --------------------------------------------------------
    # Air Properties
    p                   = Inputs.Air.p # air pressure
    # http://www.coolprop.org/fluid_properties/HumidAir.html; da means dry air
    W                   = HAPropsSI('W', 'T', T_db, 'P', p, 'R', Inputs.Air.RH)  # Humidity ratio, kg water/kg dry air

    # Transport properties of humid air from CoolProp; ha means humid air
    mu_ha               = HAPropsSI('M', 'T', T_db, 'P', p, 'W', W)
    k_ha                = HAPropsSI('K', 'T', T_db, 'P', p, 'W', W)

    # Evaluate the mass flow rate based on inlet conditions
    V_dot_ha            = Inputs.Air.V_dot_ha

    # To convert a parameter from per kg_{dry air} to per kg_{humid air}, divide by (1+W)
    W                   = HAPropsSI('W', 'T', T_db, 'P', p, 'R', Inputs.Air.RH)
    v_da                = HAPropsSI('V', 'T', T_db, 'P', p, 'W', W)
    h_da                = HAPropsSI('H', 'T', T_db, 'P', p, 'W', W)  # [J/kg]
    rho_ha              = 1 / v_da * (1 + W)    # [kg_ha/m^3]
    rho_da              = 1 / v_da              # [kg_da/m^3]
    m_dot_ha            = V_dot_ha * rho_ha     # [kg_ha/s]
    m_dot_da            = V_dot_ha * rho_da     # [kg_da/s]

    # Use a forward difference to calculate cp from cp=dh/dT
    dT                  = 0.0001                # [K]
    cp_da               = (HAPropsSI('H', 'T', T_db+dT, 'P', p, 'W', W) - h_da) / dT  # [J/kg_da/K]
    cp_ha               = cp_da / (1 + W)       # [J/kg_ha/K]

    # Check that cs_cp is defined, if so, set it to the value passed in
    if (hasattr(Inputs, 'cs_cp') and Inputs.cs_cp > 0) or (hasattr(Inputs, 'WetDry') and Inputs.WetDry == 'Wet'):
        isWet           = True
        cs_cp           = Inputs.Air.cs_cp      # c_s / c_p
    else:
        isWet           = False
        cs_cp           = 1.0

    # --------------------------------------------------------
    # Fin structure parameters
    # Dimensions and values used for both Reynolds number ranges:
    delta_f             = Inputs.Fins.t         # fin thickness(m)
    FPI                 = Inputs.Fins.FPI       # fins per inch
    FPM                 = FPI / 0.0254          # Fins per meter [1/m]
    pf                  = 1 / FPM               # Fin pitch (distance between center-lines of fins)
    F_s                 = 1 / FPM - delta_f     # Fin spacing(m)

    P_t                 = Inputs.Tubes.Pt       # transverse tube pitch (m), along the airflow direction
    P_L                 = Inputs.Tubes.Pl       # longitudinal tube pitch (m)

    D_o                 = Inputs.Tubes.OD       # Outer diameter of tube (m)
    D_c                 = D_o + 2 * delta_f     # fin collar outside diameter (m) (tube + 2*fin thickness)

    L_tube              = Inputs.Tubes.L_tube   # length of a single tube [m]
    N_fin               = L_tube * FPM          # Number of fins in the tube sheet [-]
    N_tubes_bank        = Inputs.Tubes.N_Tubes_per_bank     # tubes per bank
    Height              = P_t * (N_tubes_bank + 1)          # Height of heat exchanger [m] # assuming that fin extends 1/2 pt above/below last tube in bundle
    A_duct              = Height * L_tube                   # A_duct is the face area [m^2] - equivalent to the duct cross-section
                                                            # neglecting the additional height of the fins above/below the last tubes in the bundle
    Ac                  = A_duct - delta_f * N_fin * (Height - D_c * N_tubes_bank) - N_tubes_bank * D_c * L_tube
                                                            # Minimum duct cross-sectional area that is not fin or tube(-collar) [m^2]
    u_max               = m_dot_ha / (rho_ha * Ac)          # maximum airside velocity [m/s]
    Re_Dc               = rho_ha * u_max * D_c / mu_ha      # from reference #4 of Wang et al. Slightly different notation used in [4]

    P_d                 = Inputs.Fins.Pd                    # wave height
    X_f                 = Inputs.Fins.xf                    # projected fin length (m)
    tanTheta            = P_d / X_f                         # tangens of corrugation angle

    secTheta            = sqrt(X_f*X_f + P_d*P_d) / X_f     # !!used wavy louvered fins definition - there seems to be a bug in the paper !!
    # secTheta            = sqrt(4*X_f*X_f + P_d*P_d) / (2*X_f)  # todo check the correlations; sec = 1 / cos
    beta                = (pi * D_c**2) / (4.0 * P_t * P_L)
    oneMbeta            = 1.0 - beta                        # save one operation
    D_h                 = 2.0 * F_s * oneMbeta / (oneMbeta * secTheta + 2 * F_s * beta / D_c)   # hydraulic diameter (m)

    N                   = N_bank              # number of longitudinal tube rows (#Number of banks in ACHP-notation)

    if Re_Dc < 1000.0:
        # Heat transfer
        J1 = 0.0045-0.491*pow(Re_Dc,-0.0316-0.0171*log(N*tanTheta))*pow(P_L/P_t,-0.109*log(N*tanTheta))*pow(D_c/D_h,0.542+0.0471*N)*pow(F_s/D_c,0.984)*pow(F_s/P_t,-0.349)
        J2 = -2.72+6.84*tanTheta
        J3 = 2.66*tanTheta
        j  = 0.882*pow(Re_Dc,J1)*pow(D_c/D_h,J2)*pow(F_s/P_t,J3)*pow(F_s/D_c,-1.58)*pow(tanTheta,-0.2)

        # Friction
        F1 = -0.574-0.137*pow(log(Re_Dc)-5.26,0.245)*pow(P_t/D_c,-0.765)*pow(D_c/D_h,-0.243)*pow(F_s/D_h,-0.474)*pow(tanTheta,-0.217)*pow(N,0.035)
        F2 = -3.05*tanTheta
        F3 = -0.192*N
        F4 = -0.646*tanTheta
        f  = 4.37*pow(Re_Dc,F1)*pow(F_s/D_h,F2)*pow(P_L/P_t,F3)*pow(D_c/D_h,0.2054)*pow(N,F4)
    else:
        # Heat transfer
        j1 = -0.0545-0.0538*tanTheta-0.302*pow(N,-0.24)*pow(F_s/P_L,-1.3)*pow(P_L/P_t,0.379)*pow(P_L/D_h,-1.35)*pow(tanTheta,-0.256)
        j2 = -1.29*pow(P_L/P_t,1.77-9.43*tanTheta)*pow(D_c/D_h,0.229-1.43*tanTheta)*pow(N,-0.166-1.08*tanTheta)*pow(F_s/P_t,-0.174*log(0.5*N))
        j  = 0.0646*pow(Re_Dc,j1)*pow(D_c/D_h,j2)*pow(F_s/P_t,-1.03)*pow(P_L/D_c,0.432)*pow(tanTheta,-0.692)*pow(N,-0.737)

        # Friction
        f1 = -0.141*pow(F_s/P_L,0.0512)*pow(tanTheta,-0.472)*pow(P_L/P_t,0.35)*pow(P_t/D_h,0.449*tanTheta)*pow(N,-0.049+0.237*tanTheta)
        f2 = -0.562*pow(log(Re_Dc),-0.0923)*pow(N,0.013)
        f3 = 0.302*pow(Re_Dc,0.03)*pow(P_t/D_c,0.026)
        f4 = -0.306+3.63*tanTheta
        f  = 0.228*pow(Re_Dc,f1)*pow(tanTheta,f2)*pow(F_s/P_L,f3)*pow(P_L/D_c,f4)*pow(D_c/D_h,0.383)*pow(P_L/P_t,-0.247)

    Pr                  = cp_ha * mu_ha / k_ha
    h_a                 = j * rho_ha * u_max * cp_ha / pow(Pr,2.0/3.0)  # air side mean heat transfer coefficient using Colborn j-factor

    # calcs needed for specific fin types
    # additional parameters needed
    k_fin               = Inputs.Fins.k_fin
    N_bank              = N_bank              # Number of banks, cal one bank per iteration
    # secTheta=sqrt(X_f*X_f + P_d*P_d) / X_f  #secTheta : already calculated, no need to re-calculate it
    # Wetted Area of a single fin [m^2]
    A_1fin              = 2.0 * (Height * P_L * (N_bank+1) * secTheta  - N_tubes_bank * N_bank * pi * D_o * D_o / 4)
                                                            # assuming that fin extends 1/2 pt in front/after last tube in bundle
    # Total wetted area of the fins [m^2]
    Af                  = N_fin * A_1fin        # N_fin is the number of fin sheet in the transverse direction
    # Total wetted area including tube and fins [m^2]
    A                   = Af + N_tubes_bank * N_bank * pi * D_o * (L_tube - N_fin * delta_f)

    r                   = D_o / 2
    X_D                 = sqrt(P_L*P_L + P_t*P_t / 4) / 2
    X_T                 = P_t / 2
    rf_r                = 1.27 * X_T / r * sqrt(X_D / X_T - 0.3)
    m                   = sqrt(2 * h_a * cs_cp / (k_fin * delta_f))     # cs_cp is the correction for heat/mass transfer for a wetted surface

    # Using the circular fin correlation of Schmidt
    phi                 = (rf_r - 1) * (1 + 0.35 * log(rf_r))
    eta_f               = tanh(m * r * phi) / (m * r * phi)

    # Fin efficiency based on analysis in
    # "FIN EFFICIENCY CALCULATION IN ENHANCED FIN-AND-TUBE HEAT EXCHANGERS IN DRY CONDITIONS"
    # by Thomas PERROTIN, Denis CLODIC, International Congress of Refrigeration 2006
    # In the paper, there is no 0.1 in the cosine term, but if the cosine term is used without
    # the correction, the results are garbage for wet analysis
    # Using the offset fins correlation
    phi                 = (rf_r - 1) * (1 + (0.3+pow(m*(rf_r*r-r)/2.5, 1.5-rf_r/12.0)*(0.26*pow(rf_r,0.3)-0.3)) * log(rf_r))

    # finned surface efficiency
    eta_f               = tanh(m * r * phi) / (m * r * phi) * cos(0.1 * m * r * phi)

    # overall surface efficiency
    eta_o               = 1 - Af / A * (1 - eta_f)

    G_c                 = m_dot_ha / Ac                                 # air mass flux
    A_tube              = N_tubes_bank * N_bank * pi * D_o * L_tube     # Total outer area of the tubes [m^2]
    DeltaP_air          = (A / Ac) / rho_ha * G_c**2 / 2.0 * f          # airside pressure drop,
                            # todo: Is A/Ac = 4* L/d_H, where A is the total area; Ac is the minimum free flow area
                            # Yes, D_h = 4 * Ac * L / A, leading to A / Ac = 4 * L / D_h; thus f is the Fanning friction factor, not the Darcy one

    # write necessary values back into the given structure
    Inputs.A_a          = A
    Inputs.cp_da        = cp_da
    Inputs.cp_ha        = cp_ha

    if isWet == True:
        Inputs.eta_a_wet= eta_o
    else:
        Inputs.eta_a    = eta_o

    Inputs.h_a          = h_a
    Inputs.m_dot_ha     = m_dot_ha
    Inputs.m_dot_da     = m_dot_da
    Inputs.f_a          = f
    Inputs.dP_a         = DeltaP_air
    Inputs.Re           = Re_Dc

    return Height, secTheta, N_fin, delta_f


def HerringboneFins_Dimensions(Height, L_tube, N_bank, P_L, secTheta, t_fin, N_fin, rho_fin, D_i, D_o, N_tubes_per_bank,
                               rho_pipe, Charge_ref):
    H_total                 = Height
    L_total                 = L_tube
    W_total                 = N_bank * P_L
    Volume                  = H_total * L_total * W_total

    # Fin volume / mass calculation
    W_fin                   = W_total * secTheta
    H_fin                   = H_total
    V_fin                   = W_fin * H_fin * t_fin * N_fin
    m_fin                   = rho_fin * V_fin

    # pipe volume / mass with and without refrigerant
    # pipe
    V_1pipe                 = np.pi * (D_o**2 - D_i**2) / 4 * L_tube
    V_pipe_total            = V_1pipe * N_tubes_per_bank * N_bank
    m_pipe                  = rho_pipe * V_pipe_total
    # refrigerant, charge total
    V_1pipe_ref             = np.pi * D_i**2 / 4 * L_tube
    V_ref_total             = V_1pipe_ref * N_tubes_per_bank * N_bank
    m_ref                   = Charge_ref

    # output
    m_COND                  = m_ref + m_pipe + m_fin
    V_COND                  = Volume
    L_COND                  = L_total
    W_COND                  = W_total
    H_COND                  = H_total

    return L_COND, H_COND, W_COND, m_COND, V_COND


def HerringboneFins_Dimensions_Iter(Height, L_tube, P_L, secTheta, t_fin, N_fin, rho_fin, D_i, D_o, N_tubes_per_bank,
                               rho_pipe, Charge_ref, N_bank=1):
    H_total                 = Height
    L_total                 = L_tube
    W_total                 = N_bank * P_L
    Volume                  = H_total * L_total * W_total

    # Fin volume / mass calculation
    W_fin                   = W_total * secTheta
    H_fin                   = H_total
    V_fin                   = W_fin * H_fin * t_fin * N_fin
    m_fin                   = rho_fin * V_fin

    # pipe volume / mass with and without refrigerant
    # pipe
    V_1pipe                 = np.pi * (D_o**2 - D_i**2) / 4 * L_tube
    V_pipe_total            = V_1pipe * N_tubes_per_bank * N_bank
    m_pipe                  = rho_pipe * V_pipe_total
    # refrigerant, charge total
    V_1pipe_ref             = np.pi * D_i**2 / 4 * L_tube
    V_ref_total             = V_1pipe_ref * N_tubes_per_bank * N_bank
    m_ref                   = Charge_ref

    # output
    m_COND                  = m_ref + m_pipe + m_fin
    V_COND                  = Volume
    L_COND                  = L_total
    W_COND                  = W_total
    H_COND                  = H_total

    return L_COND, H_COND, W_COND, m_COND, V_COND


# --------------------------------------------------------------------------
def PlainFins(Inputs):
    #Source:
    #Heat transfer and friction characteristics of plain fin-and-tube heat exchangers, part II: Correlation
    #Chi-Chuan Wang, Kuan-Yu Chi, Chun-Jung Chang

        #Properties:
    p =  Inputs.Air.p #air pressure in kPa
    W=HAPropsSI('W','T',Inputs.Air.Tdb,'P',p,'R',Inputs.Air.RH)
    #Transport properties of humid air from CoolProp
    mu_ha=HAPropsSI('M','T',Inputs.Air.Tdb,'P',p,'W',W)
    k_ha=HAPropsSI('K','T',Inputs.Air.Tdb,'P',p,'W',W)
    #Evaluate the mass flow rate based on inlet conditions
    Vdot_ha =     Inputs.Air.Vdot_ha
    # To convert a parameter from per kg_{dry air} to per kg_{humid air}, divide by (1+W)
    W=HAPropsSI('W','T',Inputs.Air.Tdb,'P',p,'R',Inputs.Air.RH)
    v_da=HAPropsSI('V','T',Inputs.Air.Tdb,'P',p,'W',W)
    h_da=HAPropsSI('H','T',Inputs.Air.Tdb,'P',p,'W',W)
    rho_ha = 1 / v_da*(1+W) #[kg_ha/m^3]
    rho_da = 1 / v_da #[kg_da/m^3]
    mdot_ha = Vdot_ha * rho_ha #[kg_ha/s]
    mdot_da = Vdot_ha * rho_da #[kg_da/s]
    #Use a forward difference to calculate cp from cp=dh/dT
    dT=0.0001 #[K]
    cp_da=(HAPropsSI('H','T',Inputs.Air.Tdb+dT,'P', p, 'W',W)-h_da)/dT #[J/kg_da/K]
    cp_ha=cp_da/(1+W) #[J/kg_ha/K]

    # Check that cs_cp is defined, if so, set it to the value passed in
    if (hasattr(Inputs,'cs_cp') and Inputs.cs_cp>0) or (hasattr(Inputs,'WetDry') and Inputs.WetDry=='Wet'):
        isWet=True
        cs_cp=Inputs.Air.cs_cp
    else:
        isWet=False
        cs_cp=1.0

    #Dimensions and values used for both Reynolds number ranges:
    delta_f = Inputs.Fins.t     #fin thickness(m)
    FPI = Inputs.Fins.FPI       #fins per inch
    FPM = FPI / 0.0254         #Fins per meter [1/m]
    F_p = 1 / FPM                  #Fin pitch (distance between centerlines of fins)
    F_s = 1 / FPM - delta_f       #Fin spacing(m)

    P_t=Inputs.Tubes.Pt   #transverse tube pitch (m)
    P_L=Inputs.Tubes.Pl  #longitudinal tube pitch (m)

    D_o = Inputs.Tubes.OD      #Outer diameter of tube (m)
    D_c=D_o+2*delta_f      #fin collar outside diameter (m) (tube + 2*fin thickness)

    Ltube = Inputs.Tubes.Ltube       #length of a single tube [m]
    Nfin = Ltube * FPM #Number of fins in the tube sheet [-]
    Ntubes_bank = Inputs.Tubes.NTubes_per_bank #tubes per bank
    Height = P_t * (Ntubes_bank+1)  #Height of heat exchanger [m] # assuming that fin extends 1/2 pt above/below last tube in bundle
    A_duct = Height * Ltube  #A_duct is the face area [m^2] - equivalent to the duct cross-section #neglecting the additional height of the fins above/below the last tubes in the bundle
    A_c = A_duct - delta_f * Nfin * (Height-D_c*Ntubes_bank) - Ntubes_bank * D_c * Ltube #Minimum duct cross-sectional area that is not fin or tube(-collar) [m^2]
    u_max = mdot_ha / (rho_ha * A_c) #maximum airside velocity [m/s]
    Re_Dc=rho_ha*u_max*D_c/mu_ha #from reference #4 of Wang et all. Slightly different notation used in [4]

    N=Inputs.Tubes.Nbank  #number of lomgitudunal tube rows (#Number of banks in ACHP-notation)
    L=(N+1)*P_L #depth of heat exchanger in airflow direction, assuming fins extends 1/2 bank oer the first and last tube

    Nbank=N #to be able to reuse equations
    #Wetted Area of a single fin [m^2]
    A_1fin = 2.0 * (Height * P_L * (Nbank+1)   - Ntubes_bank*Nbank * pi*D_c*D_c/4) #assuming that fin extends 1/2 pt in front/after last tube in bundle
    # Total wetted area of the fins [m^2]
    Af = Nfin * A_1fin
    #Total area including tube and fins [m^2]
    A_o = Af + Ntubes_bank * Nbank * pi * D_c * (Ltube-Nfin*delta_f)

    D_h=4*A_c*L/A_o #Hydraulic diameter

    #heat transfer
    if N==1:
        P1=1.9-0.23*log(Re_Dc)
        P2=-0.236+0.126*log(Re_Dc)
        j=0.108*pow(Re_Dc,-0.29)*pow(P_t/P_L,P1)*pow(F_p/D_c,-1.084)*pow(F_p/D_h,-0.786)*pow(F_p/P_t,P2)
    if N>=2:
        P3=-0.361-0.042*N/log(Re_Dc)+0.158*log(N*(F_p/D_c)**0.41)
        P4=-1.224-0.076*pow(P_L/D_h,1.42)/log(Re_Dc)
        P5=-0.083+0.058*N/log(Re_Dc)
        P6=-5.735+1.21*log(Re_Dc/N)
        j=0.086*pow(Re_Dc,P3)*pow(N,P4)*pow(F_p/D_c,P5)*pow(F_p/D_h,P6)*pow(F_p/P_t,-0.93)
    #pressure drop
    F1=-0.764+0.739*P_t/P_L+0.177*F_p/D_c-0.00758/N
    F2=-15.689+64.021/log(Re_Dc)
    F3=1.696-15.695/log(Re_Dc)
    f=0.0267*pow(Re_Dc,F1)*pow(P_t/P_L,F2)*pow(F_p/D_c,F3)


    Pr = cp_ha * mu_ha / k_ha
    h_a = j * rho_ha * u_max * cp_ha / pow(Pr,2.0/3.0) #air side mean heat transfer coefficient using colborn j-factor

    #calcs needed for specific fin types
    #additional parameters needed
    k_fin =  Inputs.Fins.k_fin

    #Total area including tube and fins [m^2]
    A = A_o

    r = D_o / 2
    X_D = sqrt(P_L*P_L + P_t*P_t / 4) / 2
    X_T = P_t / 2
    rf_r = 1.27 * X_T / r * sqrt(X_D / X_T - 0.3)
    m = sqrt(2 * h_a * cs_cp / (k_fin * delta_f)) #cs_cp is the correction for heat/mass transfer

    #Using the circular fin correlation of Schmidt
    phi = (rf_r - 1) * (1 + 0.35 * log(rf_r))
    eta_f = tanh(m * r * phi) / (m * r * phi)

    # Fin efficiency based on analysis in
    # "FIN EFFICIENCY CALCULATION IN ENHANCED FIN-AND-TUBE HEAT EXCHANGERS IN DRY CONDITIONS"
    # by Thomas PERROTIN, Denis CLODIC, International Congress of Refrigeration 2006
    # In the paper, there is no 0.1 in the cosine term, but if the cosine term is used without
    # the correction, the results are garbage for wet analysis
    # Using the offset fins correlation
    phi = (rf_r - 1) * (1 + (0.3+pow(m*(rf_r*r-r)/2.5,1.5-rf_r/12.0)*(0.26*pow(rf_r,0.3)-0.3)) * log(rf_r))

    #finned surface efficiency
    eta_f = tanh(m * r * phi) / (m * r * phi)*cos(0.1 * m * r * phi)

    #overall surface efficiency
    eta_o = 1 - Af / A * (1 - eta_f)

    G_c=mdot_ha/A_c #air mass flux
    Atube = Ntubes_bank * Nbank * pi * D_o * Ltube# Total outer area of the tubes [m^2]
    DeltaP_air=A/A_c/rho_ha*G_c**2/2.0*f #airside pressure drop

    #write necessary values back into the given structure
    Inputs.A_a=A;
    Inputs.cp_da=cp_da
    Inputs.cp_ha=cp_ha
    if isWet==True:
        Inputs.eta_a_wet=eta_o
    else:
        Inputs.eta_a=eta_o
    Inputs.h_a=h_a
    Inputs.mdot_ha=mdot_ha
    Inputs.mdot_da=mdot_da
    Inputs.f_a=f
    Inputs.dP_a=DeltaP_air
    Inputs.Re=Re_Dc


# ----------------------------------------------------
# empty class
# ------------------------------------------------------
class FinsVals:
    pass


class TubesVals:
    pass


class AirVals:
    pass


def IsFinsClass(Fins):
    """
    Returns the Fins class if Fins is an instance of the FinInputs class, False otherwise
    Convenience function for the Validator
    """
    if isinstance(Fins, FinInputs):
        return Fins
    else:
        return False

def IsFinsClass_Iter(Fins):
    """
    Returns the Fins class if Fins is an instance of the FinInputs class, False otherwise
    Convenience function for the Validator
    """
    if isinstance(Fins, FinInputs_Iter):
        return Fins
    else:
        return False
