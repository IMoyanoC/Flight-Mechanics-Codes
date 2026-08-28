
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import g
import math as mt
import sympy as sp
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
from scipy.integrate import quad

#==========================================================================
#------------------------ AIRCRAFT PARAMETERS -----------------------------
#==========================================================================

#------Parametros De Aeronave-----
from Duke import *
#------Parametros Generales-----

w = M * g # Peso de la aeronave
def rho(h):
    rho = 1.225 * (1-2.2558e-5 * h)**4.2559
    return rho

#CURVA CL(al_pha) del avión completo
K_I = (1 + 2.15 * (b_f / b)) * (S_net / S) + (mt.pi/2) * (b_f ** 2) / (CLa_w * S)
K_II = (1 + 0.7 * b_f/b) * S_net/S
CLa_wf = K_I * CLa_w # pendiente sust ala fuselaje
print('Pendiente de sustentación combinación ala-fuselaje: CLa_wf = ', CLa_wf)
CLa = CLa_wf * (1 + (x_cg- x_ac) * CAM / l_h) #Pendiente de sustentacion del avion completo
CL0 = (w) / (0.5 * rho(hc) * V_crucero**2 * S) # Coeficiente de sustentación en condiciones de crucero
al_pha_0 = -CL0 / CLa
def CL(al_pha):
  CL = CL0 + CLa * al_pha
  return CL
al_pha_values = np.linspace(np.radians(-5), np.radians(18), 100)
al_pha_values_degr = np.linspace(-5, 18, 100)
cl_al_pha_values = [CL(al_pha) for al_pha in al_pha_values]
plt.figure(figsize=(10, 5))
plt.plot(al_pha_values_degr, cl_al_pha_values, label=r'$C_L(\alpha)$', color='blue')
plt.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
plt.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
plt.xlabel(r'$\alpha$ (grados)')
plt.ylabel('$C_L$')
plt.title('Curva de sustentación para vuelo trimado')
plt.legend()
plt.grid(True)
plt.savefig("CL_al_pha.pdf", format='pdf', dpi=300)
plt.show()
#CURVA CM(CL) interferencia ala-fuselaje
D_f_CMac =-1.8*(1- (2.5 * b_f)/l_f) * (np.pi * b_f * h_f * l_f)/(4 * S * CAM) * (CL0/(CLa_wf))
D_f_CMac = D_f_CMac*(np.pi/4)*b_f*h_f
Cmac_wf = Cm_ac_w + D_f_CMac
D_f1 =-(1.8/CLa_wf)*(b_f*h_f*l_fn/(S*CAM))
D_f2 = (0.273/(1+taper)) * ((b_f * x_cg *(b- b_f)) / (CAM**2 *(b + 2.15* b_f))) * mt.tan(mt.radians(sweep1_4))
x_ac_wf = x_ac + D_f1 + D_f2
def Cm_wf(CL_wf):
  Cm_wf = Cmac_wf + CL_wf * (x_cg- x_ac_wf)
  return Cm_wf
CL_wf_values = np.linspace(-2, 15)
Cm_wf_values = [Cm_wf(CL_wf) for CL_wf in CL_wf_values]
#Curva CM(CL) de la configuración completa
D_n =-4*((b_n**2*l_n)/(S*CAM*CLa_wf)) #Correccion del centro aerodinamico por presencia de barquillas
D_p=-0.05*((N_p*dp**2*l_p)/(S*CAM*CLa_wf)) #Cambio del centro aerodinamico por presencia de la planta propulsora
x_ac_completa = x_ac_wf + D_n + D_p
CLa2 = CLa_wf + CLa_h * (1- de_da) * (S_h/S) * n_h
xn_c = x_ac_completa + (CLa_h/CLa2)*(1-de_da)*((S_h * l_h)/(S*CAM))*n_h
i_hf = ((Cm_ac_w + CL0*(x_cg-x_ac_completa))/((CLa_h*S_h * l_h * n_h)/(S*CAM))) + (de_da / CLa_w) * CL0
i_h = i_hf- CL0/CLa
Cm0 = Cm_ac_w- CLa_h*(i_h)* (S_h * l_h) * n_h / (S * CAM)
dCm_dCL =- (xn_c- x_cg)
def CM(CL):
  C_M = Cm0 + dCm_dCL * CL
  return C_M
CL_values = np.linspace(-2, 15)
CM_values = [CM(CL) for CL in CL_values]
#Graficos de CM(Cl)
#Ala:
def CM_w(CL):
  CM_w = Cm_ac_w + CL * (x_cg- x_ac)
  return CM_w
#Ala-fuselaje
def CM_wf(CL):
  Cm_wf = Cmac_wf + CL * (x_cg- x_ac_wf)
  return Cm_wf
#Avión sin EH
def CM_seh(CL):
  CM_A_h = Cmac_wf + CL* (x_cg- x_ac_completa)
  return CM_A_h
#Configuración completa
def CM(CL):
  C_M = Cm0 + dCm_dCL * CL
  return C_M
#graficos
CL_values = np.linspace(0, 1.6)
CM_w_values = [CM_w(CL) for CL in CL_values]
CM_wf_values = [CM_wf(CL) for CL in CL_values]
CM_seh_values = [CM_seh(CL) for CL in CL_values]
CM_values = [CM(CL) for CL in CL_values]
plt.figure(figsize=(12,6))
plt.plot(CL_values,CM_w_values, label='Ala')
plt.plot(CL_values, CM_wf_values, label = 'Ala-fuselaje' )
plt.plot(CL_values, CM_seh_values, label='Avión sin Empenaje Horizonal')
plt.plot(CL_values, CM_values, label='Configuración Completa')
plt.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
plt.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
plt.xlabel('CL')
plt.ylabel('CM')
plt.title('Curvas de momento para distintas configuraciones')
plt.legend()
plt.grid(True)
plt.savefig("CM(CL).pdf", format='pdf', dpi=300)
plt.show()

########################################################################################################

########################################################################################################

########################################################################################################

########################################################################################################

########################################################################################################

#Polares


#VARIABLES AGREGADAS POR NOSOTROS
#Consideracion de efectos de compresibilidad
Mach=0.39
sweep1_2 = sweep1_4-(1/A)*((1-taper)/(1+taper))
tg_sweep_compresible= sweep1_2/(np.sqrt(1-Mach**2))
Wto=4500 #[kg]

visc = 1.62e-5
beta_graficos =np.sqrt(1-Mach**2)
Re_w = (rho(hc) * V_crucero * CAM) / visc
Re_f = (rho(hc) * V_crucero * l_f) / visc
Re_h = (rho(hc) * V_crucero * CAMh) / visc
Re_b = (rho(hc) * V_crucero * l_n) / visc
Re_v = (rho(hc) * V_crucero * CAMv) / visc
Ref_local = (rho(hc) * V_crucero * l_fn) / visc

#Parametros de flap
Cla_flap = 6.45
cd_po=0.00587 #RESISTENCIA PERFIL FLAP
bf= 3.4 #envergadura flap
Sf= 4.32 #superficie flap ES EL SWF DE LA GRAFICA
bfi= 0.9 #estacion interna flap
bfo= 4.3 #estacion externa flap
k2_3 = 1.15
cf=Sf/bf
cf_c= cf/CAM
Dc_cf= 0.250
cprim_c= 1 + Dc_cf * cf_c
cf_cprim= cf_c / cprim_c
kd= 0.12 #DE TABLA PAG3
deltaf = np.radians(30)
kl= 1
eta_delta = 0.57
Kc= 1.05
Kb= 0.32
#DE TABLA
#DE TABLA
#DE TABLA
#RESISTENCIA POR SUSTENTACION DE LOS COMPONENTES
#RESISTENCIA INDUCIDA POR VORTICES- ALA SIN ALABEO (METODO B)
Cte1=0.405
Cte2=0.278
Cte3=0.335
F=(2*np.pi * A) /(Cla * np.cos(sweep1_4))
integ1= (1+2*taper)/(3*(1+taper))
integ2= 4/(3*np.pi)+0.001*mt.atan(tg_sweep_compresible)
eta_cp= Cte1*integ1 + Cte2*4/(3*np.pi) +Cte3*integ2
delta=46.264 * (eta_cp-4/(3*np.pi))**2
C_vorticoso_sin_alabeo= (1+delta)/(np.pi*A)
#RESISTENCIA POR ALABEO A SUSTENTACION NULA
e_tip =-0.0084
#Alabeo AERO EN GRADOS
C_0_1= 2.6e-5
C_1_1= 3.7e-4
A_vorticoso_con_alabeo=e_tip**2*C_0_1
B_vorticoso_con_alabeo=e_tip*C_1_1
#RESISTENCIA POR SUSTENTACION DL FUSELAJE
CL= sp.symbols('CL')
al_pha_f= (CL-CLOw)/CLaw
Cd_S_f = 0.15 * al_pha_f**2 * Vf**(2/3) / S
Cd_S_f_desarrollado = sp.expand(Cd_S_f)
polinomio = sp.Poly(Cd_S_f_desarrollado, CL)
coeficientes = polinomio.coeffs()
C_vorticos_fuselaje = coeficientes[0] if len(coeficientes) > 0 else 0
B_vorticos_fuselaje = coeficientes[1] if len(coeficientes) > 1 else 0
A_vorticos_fuselaje = coeficientes[2] if len(coeficientes) > 2 else 0
#RESISTENCIA POR SUSTENTACION DEL EMPENAJE
x_ac= 0.16751622805445915
Cmac=-0.25479659918747544
Cl_h= (Cmac+CL*(x_cg-x_ac))/(Sh*l_h/(S*CAM))
CD_S_h= (1.02*Cl_h**2 * Sh/(np.pi * Ah))/S
CD_S_h_desarrollado = sp.expand(CD_S_h)
polinomio = sp.Poly(CD_S_h_desarrollado, CL)
coeficientes = polinomio.coeffs()
C_vorticos_empenaje = coeficientes[0] if len(coeficientes) > 0 else 0
B_vorticos_empenaje = coeficientes[1] if len(coeficientes) > 1 else 0
A_vorticos_empenaje = coeficientes[2] if len(coeficientes) > 2 else 0
#RESISTENCIA DE PERFIL
#RESISTENCIA DE LOS PERFILES - ALA
CF_lam_w= 1.33/np.sqrt(Re_w)
CF_turb_w= 0.455/np.log10(Re_w)**2.58
x_trans = 0.175*CAM
phi_w = 2.7*t_c + 100*t_c**4
Cdp_min_lam= 2*CF_lam_w*(1+phi_w*np.cos(sweep1_2)**2)
Cdp_min_turb= 2*CF_turb_w*(1+phi_w*np.cos(sweep1_2)**2)
Cdp_min= 0.175*Cdp_min_lam + (1-0.175)*Cdp_min_turb
Cl = sp.symbols('Cl')
DlCdp_ref = (67*Clmax)/(np.log10(Re_w))**4.5- 0.0046*(1 + 2.75*t_c)
DlCdp = 0.75 * DlCdp_ref * ((CL- Cli)/(Clmax- Cli))**2
DlCdp_desarrollado = sp.expand(DlCdp)
polinomio = sp.Poly(DlCdp_desarrollado, CL)
coeficientes = polinomio.coeffs()
C_perfiles_ala = coeficientes[0] if len(coeficientes) > 0 else 0
B_perfiles_ala = coeficientes[1] if len(coeficientes) > 1 else 0
A_perfiles_ala = Cdp_min*(Swet_w/S) + coeficientes[2] if len(coeficientes) > 2 else 0
#RESISTENCIA DE LOS PERFILES- FUSELAJE
Ac = 2.5
l_n = 2.5
la = 4.5
Deff = np.sqrt(4/np.pi * Ac)
taper_eff = min(l_f/Deff , (l_n + la)/Deff +2)
phi_f = 2.2/taper_eff**1.5 + 3.8/taper_eff**3
CF_turb_f= 0.455/np.log10(Re_f)**2.58
CD_S_b = CF_turb_f * Swet_f * (1 + phi_f)
A_perfiles_fuselaje = 1/S * CD_S_b
#CORRECCION POR ANGULO DE INCIDENCIA Y DE COLA
A_I = 16
A_II = 6.5
beta = np.radians(13)
#al_pha_f_values = np.linspace(np.rad(-5), np.rad(15))
Dab_CD_S = A_I * sp.Abs(sp.sin(al_pha_f)** 3) + A_II / sp.cos(beta) * (sp.Abs(sp.sin(al_pha_f- beta)** 3))
D1 = 0.0005
D2 = 1.7124
al_pha_fmin = (beta*np.sqrt(A_I/A_II)-1)/(A_I/A_II-1)
Dab_CD_S_ajust = (D1 + D2 * (al_pha_f- al_pha_fmin)**2)
Dab_CD_S_ajust_exp = sp.expand(Dab_CD_S_ajust) /S
#sp.plot(Dab_CD_S_ajust_exp, (al_pha_f- al_pha_fmin) ** 2, xlabel='', ylabel='')
A_perfilesf_upsweep= 0 #0.001
B_perfilesf_upsweep= 0 #-0.004902
C_perfilesf_upsweep= 0 #0.006
#RESISTENCIA DE LAS BARQUILLAS
CF_turb_b= 0.455/np.log10(Re_b)**2.58
taper_effb = l_n / b_n
#NO LO PUDE HACER FUNCIONAR
CD_S_n = CF_turb_b *(1 + 2.2 / taper_effb**1.5 + 3.8 / taper_effb**3.8 ) * Swet_barqui
A_perfiles_barquilla = 2 * CD_S_n / S
#RESISTENCIA DE PERFIL DEL EMPENAJE HORIZONTAL
CF_turb_h= 0.455/np.log10(Re_h)**2.58
CD_S_hbas = 2*CF_turb_h*(1+2.75*t_c_h*np.cos(sweep1_2h)**2)*Sh
Dl_CD_S_h = 0.33 * Cl_h**2 / (np.pi*Ah*(np.cos(sweep1_2h))) * Sh
Dl_CD_S_h_desarrollado = sp.expand(Dl_CD_S_h)
CD_S_h_perf = (CD_S_hbas + Dl_CD_S_h_desarrollado) / S
polinomio = sp.Poly(CD_S_h_perf, CL)
coeficientes = polinomio.coeffs()
C_perfiles_empenaje = coeficientes[0] if len(coeficientes) > 0 else 0
B_perfiles_empenaje = coeficientes[1] if len(coeficientes) > 1 else 0
A_perfiles_empenaje = coeficientes[2] if len(coeficientes) > 2 else 0
#RESISTENCIA DE PERFIL DEL EMPENAJE VERTICAL
CF_turb_v = 0.455/np.log10(Re_v)**2.58
CD_S_vbas = 2*CF_turb_v* (1+ 2.75*t_c_v*np.cos(sweep1_2v)**2) * Sv
A_perfiles_EV = CD_S_vbas / S
#INTERFERENCIA Y CORRECCIONES
#INTERFERENCIA ALA- FUSELAJE
#Inducida por vortices
Di_CDv= (0.55* (Deff /b)) / (1+taper) * (2- np.pi * (Deff /b)) * CLOw**2 / (np.pi * A)
A_interferencia_af_vort=Di_CDv
#Inducida por viscosidad
C_ci=4.5* Cr
CF_turb_f_local= 0.455/np.log10(Ref_local)**2.58
Di_CD_Sp= 1.5*CF_turb_f_local*tr*C_ci*np.cos(sweep1_2)**2
A_interferencia_af_visc = Di_CD_Sp/S
B_interferencia_af_visc = CF_turb_f_local*Cr*Deff / S
#Fuselaje por presencia del ala
Di_CD_Sp_fa=2*beta*np.cos(sweep1_2)/A * D2*CL/CLaw
B_interferencia_fa= Di_CD_Sp_fa / S
B_interferencia_fa_desarrollado = sp.expand(B_interferencia_fa)
polinomio = sp.Poly(B_interferencia_fa_desarrollado, CL)
coeficientes = polinomio.coeffs()
B_interferencia_fa = coeficientes[0] if len(coeficientes) > 0 else 0
#Interferencia avion-Barquilla
#Avion a helice
Di_CD_Sn= 0.004*Sfront_barqui * 2
A_interferencia_avb= Di_CD_Sn / S
#INTERFERENCIA AVION-EMPENAJE
de_dCL= de_da*(1/CLaw)
Di_CD_Sh= Cl_h*CL * (de_dCL-2/(np.pi*A)) * Sh / S
Di_CD_Sh_desarrollado = sp.expand(Di_CD_Sh)
polinomio = sp.Poly(Di_CD_Sh, CL)
coeficientes = polinomio.coeffs()
C_avion_empenaje = coeficientes[0] if len(coeficientes) > 0 else 0
B_avion_empenaje = coeficientes[1] if len(coeficientes) > 1 else 0
#SUMAS
#RESSITENCIA INDUCIDA POR VORTICES
#Resistencia ala sin alabeo
C_total_sustentacion=C_vorticoso_sin_alabeo+C_vorticos_fuselaje+C_vorticos_empenaje
B_total_sustentacion= B_vorticos_fuselaje + B_vorticos_empenaje + B_vorticoso_con_alabeo
A_total_sustentacion= A_vorticos_fuselaje+A_vorticos_empenaje+A_vorticoso_con_alabeo
#Resitencia de perfiles
C_total_perfiles= C_perfiles_ala + C_perfilesf_upsweep + C_perfiles_empenaje
B_total_perfiles= B_perfiles_ala+ B_perfilesf_upsweep+ B_perfiles_empenaje
A_total_perfiles= A_perfiles_ala + A_perfiles_fuselaje + A_perfilesf_upsweep + A_perfiles_barquilla +A_perfiles_empenaje+A_perfiles_EV
#Interferencia
C_total_interferencia= C_avion_empenaje
B_total_interferencia= B_interferencia_af_visc + B_interferencia_fa + B_avion_empenaje
A_total_interferencia= A_interferencia_af_vort+ A_interferencia_af_visc + A_interferencia_avb
#PROTUBERANCIAS
C_imperf_Ala= 0.06 * (C_perfiles_ala)
B_imperf_Ala=0.06 * (B_perfiles_ala)
A_imperf_Ala=0.06 * (A_perfiles_ala)
C_imperf_FH= 0.07 * (C_perfilesf_upsweep)
B_imperf_FH=0.07 * (B_perfilesf_upsweep)
A_imperf_FH=0.07 * (A_perfiles_fuselaje+A_perfilesf_upsweep)
A_inst_motor=0.15 * (A_perfiles_barquilla)
C_Instalaciones= 0.06 * (C_perfilesf_upsweep)
B_Instalaciones=0.06 * (B_perfilesf_upsweep)
A_Instalaciones=0.06 * (A_perfiles_fuselaje+A_perfilesf_upsweep)
#SUMA
C_total_protu= C_imperf_Ala+C_imperf_FH+C_Instalaciones
B_total_protu=B_imperf_Ala+B_imperf_FH+B_Instalaciones
A_total_protu=A_imperf_Ala+A_imperf_FH+A_inst_motor+A_Instalaciones
#SUMA TOTAL
C_total= C_total_protu+ C_total_interferencia+C_total_perfiles+C_total_sustentacion
B_total=B_total_protu+ B_total_interferencia+B_total_perfiles+B_total_sustentacion
A_total=A_total_protu+ A_total_interferencia+A_total_perfiles+A_total_sustentacion
Polar_crucero= C_total*(CL)**2 + B_total*CL + A_total
#sp.plot(CL, Polar_crucero, xlabel='', ylabel='', xlim=(-1.5,1.5),ylim=(0,0.18) )
Polar_crucero= C_total * CL**2 + B_total * CL + A_total
Polar_crucero_func = sp.lambdify(CL, Polar_crucero, 'numpy')
#Gráfica Polar Crucero
CL_vals = np.linspace(-1.5, 1.5, 400)
Polar_crucero_vals = Polar_crucero_func(CL_vals)
plt.figure(figsize=(10, 6))
plt.plot(CL_vals, Polar_crucero_vals, label='Curva Polar Crucero', color='red')
plt.xlabel('Coeficiente de Sustentación (CL)', fontsize=14)
plt.ylabel('Coeficiente de Drag (CD)', fontsize=14)
plt.title('Curva Polar de vuelo crucero', fontsize=16, fontweight='bold')
plt.xlim(-1.5, 1.5)
plt.ylim(0, 0.18)
plt.axhline(0, color='black', linewidth=0.5)
plt.axvline(0, color='black', linewidth=0.5)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend(loc='lower left', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
plt.savefig("PolarCrucero.pdf", format='pdf', dpi=300)
plt.show()
#---------------------------POLAR BAJA VELOCIDAD------------------------------
#Incremento de resistencia de perfil del ala por deflexion de flap
titaf= np.arccos(2*cf_cprim-1)
al_pha_delta_prima= 1-(titaf-np.sin(titaf))/np.pi
Df_cd_po= kd*Cla_flap* al_pha_delta_prima* cf_c * deltaf * np.sin (deltaf) + cd_po*(cprim_c- 1)
Df_Clo_prim= Cla_flap*eta_delta* al_pha_delta_prima*deltaf
Df_Clo = Df_Clo_prim * cprim_c + clor * (cprim_c-1)
Df_CLO= Df_Clo * CLaw / Clar * Kc * Kb
Df_CD_perfil= k2_3 * Sf/Swet_w * Df_cd_po * np.cos(sweep1_4)- cd_po * kl*Df_CLO*(CL-(CLOw+0.25*Df_CLO))
Df_CD_perfil_desarrollado = sp.expand(Df_CD_perfil)
polinomio = sp.Poly(Df_CD_perfil, CL)
coeficientes = polinomio.coeffs()
B_perfil_flap = coeficientes[0] if len(coeficientes) > 0 else 0
A_perfil_flap = coeficientes[1] if len(coeficientes) > 1 else 0
#RESISTENCIA INDUCIDA POR VORTICES
Kff=1/3
w= 0.007 #DE TABLA
z= 0.07/(1+taper)*(1-Kff)**2 * bfi/b
v= 0.0095 #DE TABLA
Df_Cd_v=(w+z)*Df_Clo**2 + v*CL * Df_Clo
Df_Cd_v_desarrollado = sp.expand(Df_Cd_v)
polinomio = sp.Poly(Df_Cd_v, CL)
coeficientes = polinomio.coeffs()
B_vorticoso_flap = coeficientes[0] if len(coeficientes) > 0 else 0
A_vorticoso_flap = coeficientes[1] if len(coeficientes) > 1 else 0
#RESISTENCIA DE TRIMADO
mhu1=0.225 #DE TABLA
mhu2= 0.41 #DE TABLA
mhu3= 0.04 #DE TABLA
cl=CL + Df_Clo * (1-Sf/S)
Df_cm1_4=-mhu1*Df_Clo * cprim_c- cl/8 * cprim_c* (cprim_c-1)
Df_Cm1_4= mhu2 * Df_cm1_4 + 0.7*A/(1+2/A)*mhu3*Df_Clo * np.tan(sweep1_4)
eh= 1-(0.25/np.cos(sweep1_4h)**2)
Dtrim_CD= Df_Cm1_4*(Df_Cm1_4 + 2*CMO)/(np.pi*Ah*eh*(l_h/CAMh)**2 * (Sh/S))
Dtrim_CD_desarrollado = sp.expand(Dtrim_CD)
polinomio = sp.Poly(Dtrim_CD, CL)
coeficientes = polinomio.coeffs()
C_trimado = coeficientes[0] if len(coeficientes) > 0 else 0
B_trimado = coeficientes[1] if len(coeficientes) > 1 else 0
A_trimado = coeficientes[2] if len(coeficientes) > 2 else 0
#4. Resistencia debido al tren de aterrizaje
CD_uc=7e-4 * Wto**(0.7859)/S
A_tren= CD_uc
#SUMA
C_total_flap= C_trimado
B_total_flap= B_perfil_flap+B_vorticoso_flap +B_trimado
A_total_flap= A_perfil_flap+ A_vorticoso_flap+ A_trimado+A_tren
#SUMA FINAL PARA AVION A BAJA VELOCIDAD
C_total_LS= C_total_flap + C_total
B_total_LS= B_total_flap + B_total
A_total_LS= A_total_flap + A_total
Polar_LS= C_total_LS*(CL)**2 + B_total_LS*CL + A_total_LS
#sp.plot(CL, Polar_LS, xlabel='', ylabel='', xlim=(-1.5,1.5),ylim=(0,0.18) )
#Gráfico Polar Baja Velocidad (delta=30grados)
polar_LS_func = sp.lambdify(CL, Polar_LS, 'numpy')
Polar_LS_vals = polar_LS_func(CL_vals)
plt.figure(figsize=(10, 6))
plt.plot(CL_vals, Polar_LS_vals, label='Curva Polar Baja Velocidad', color='green')
plt.xlabel('Coeficiente de Sustentación (CL)', fontsize=14)
plt.ylabel('Coeficiente de Drag (CD)', fontsize=14)
plt.title('Curva Polar de vuelo crucero', fontsize=16, fontweight='bold')
plt.xlim(-1.5, 1.5)
plt.ylim(0, 0.18)
plt.axhline(0, color='black', linewidth=0.5)
plt.axvline(0, color='black', linewidth=0.5)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend(loc='lower left', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
plt.savefig("PolarLS.pdf", format='pdf', dpi=300)
plt.show()  