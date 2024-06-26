from __future__ import division, print_function  # Make integer 3/2 give 1.5 in python 2.x

import CoolProp as CP
import numpy as np
import pylab as pylab
from math import pi
# from PyACHP.Correlations import AccelPressureDrop, LockhartMartinelli
from ACHP_codes.Correlations.deltaP_Correlations import AccelPressureDrop, LockhartMartinelli, \
    Kim_Mudawar_condensing_DPDZ_h

# flow, dimensions and heat flux related variables
G_r = 300  # Mass velocity [kg/m-s]
D = 0.01  # tube diameter[m]
Q = 200  # W/m^s
L = 1  # length of tube
m_dot = (G_r * pi * D ** 2) / 4.0
Ref = "R134a"
Backend = 'HEOS'  # choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
psat_r = 702800  # Pa
AS = CP.AbstractState(Backend, Ref)
AS.update(CP.PQ_INPUTS, psat_r, 0.0)
Tbubble_r = AS.T()  # [K]
rho_l = AS.rhomass()  # [kg/m^3]
h_l = AS.hmass()  # [J/kg]
AS.update(CP.PQ_INPUTS, psat_r, 1.0)
Tdew_r = AS.T()  # [K]
rho_v = AS.rhomass()  # [kg/m^3]
h_v = AS.hmass()  # [J/kg]
print(AS.p_critical())  # [Pa]
print(Tdew_r, Tbubble_r)
# print(PropsSI('H','P',702800,'Q',1.0,Ref))
print(rho_l)
h_fg = h_v - h_l
print("h_fg", h_fg)
v_f = 1 / rho_l  # liquid density [m^3/kg]
v_g = 1 / rho_v  # gas density [m^3/kg]
v_fg = v_g - v_f  # difference between the two of them [m^3/kg]
DELTA_x = Q * D * pi * L / (h_fg * m_dot)  # change of quality in section

# generate data
x = np.linspace(0, 1, 1000)
dp_acc_zivi = np.zeros_like(x)  # accelerational pressure drop, Zivi slip flow model
dp_acc_hom = np.zeros_like(x)  # accelerational pressure drop, homogeneous equilibrium model
dp_lm = np.zeros_like(x)  # frictional pressure drop (Lockhardt-Martinelli)
dp_km = np.zeros_like(x)  # frictional condensation pressure drop (Kim-Mudawar)
for i in range(len(x)):
    # determine inlet/outlet quality
    x_in = x[i]
    x_out = x_in + DELTA_x
    if x_out > 1.0:  # last point in graph needs to be moved to the left
        x_in = 1.0 - DELTA_x
        x_out = 1.0
    # accelerational pressure drop using ACHP
    dp_acc_zivi[i] = -AccelPressureDrop(x_in, x_out, AS, G_r, Tbubble_r, Tdew_r, None, None, 'Zivi') * L
    # accelerational pressure drop using homogeneous equilibrium model:
    dp_acc_hom[i] = -AccelPressureDrop(x_in, x_out, AS, G_r, Tbubble_r, Tdew_r, None, None, 'Homogeneous') * L
    # frictional pressure drop using Lockhart-Martinelli
    C = None
    satTransport = None
    dp_lm[i] = LockhartMartinelli(AS, G_r, D, x[i], Tbubble_r, Tdew_r, C, satTransport)[0] * L
    beta = 1  # channel aspect ratio (=width/height)
    dp_km[i] = Kim_Mudawar_condensing_DPDZ_h(AS, G_r, D, x[i], Tbubble_r, Tdew_r, psat_r, beta, satTransport)[0] * L
    # havg=np.trapz(h,x=x)


def plot(KM=None, LM=None, Zivi=None, Hom=None, scl='log', txtx=0.1, txty=5, pos='best'):
    pylab.figure(figsize=(7, 5))
    pylab.text(txtx, txty,
               '%s\nG=%0.1f kg/m$^2$\nD=%0.3f, L=%0.3f m\nq\"=%0.1f W/m$^2$\np=%0.1f Pa' % (Ref, G_r, D, L, Q, psat_r))
    if LM: pylab.plot(x, dp_lm, label="Frictional: Lockhart-Martinelli")
    if Zivi: pylab.plot(x, dp_acc_zivi, label="Accelerational: Zivi")
    if Hom: pylab.plot(x, dp_acc_hom, '--', label="Accelerational: HEM")
    if KM: pylab.plot(x, dp_km, '--y', label="Frictional Condensation: Kim-Mudawar")
    leg = pylab.legend(loc=pos, fancybox=True)
    leg.get_frame().set_alpha(0.5)
    pylab.title('Pressure drop components as a function of quality')
    pylab.gca().set_xlabel('Quality in [-]')
    pylab.gca().set_ylabel(r'Pressure drop in [Pa]')
    pylab.gca().set_yscale(scl)
    pylab.gca().set_ylim(0, None)
    # pylab.subplots_adjust(left=0.18,bottom=0.1, right=0.95, top=0.9


plot('KM', 'LM', 'Zivi', 'Hom', 'linear', 0.5, 2000)
plot('KM', 'LM', 'Zivi', 'Hom', 'log', 0.08, 5.5)
plot(None, None, 'Zivi', 'Hom', 'linear', 0.6, 1.5)
pylab.show()
