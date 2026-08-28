import numpy as np
#------Parametros Geometricos-----
#------Parametros Planta Alar-----
S = 19.72 # superficie alar
S_net = 23.368 # Superficie expuesta
Swet_w= 23.368 # Superficie mojada ala
CAM = 1.92 # Cuerda aerodinamica media
CGM=1.781 # Cuerda media geometrica
b = 15.33 # Envergadura
taper = 0.374 # Ahusamiento
A = 8.61 # Alargamiento
sweep1_4= np.radians(0) # Angulo sweep cuarto cuerda
clt_alf_2 = 6.5855 # Pendiente sustentacion (puntera)
cl_alf_1 = 6.3810 # Pendiente sustentacion (raiz)
Cm_ac_w =-0.016 # Coeficiente momento libre ala
alabeo = np.radians(-3.95) # Alabeo aerodinamico total
delta_z_CL = 0.1 #Porque es ala baja

#-----------Parámetros condición vuelo crucero---------
M = 4000 # Masa aeronave
V_crucero = 113.18#100.55 # velocidad crucero en m/s
hc = 8850 #7925 # altura crucero
x_cg = 0.3
x_ac = 0.248 # Posición centro de gravedad adimensionada con la CAM


#------------Parámetros del fuselaje------------
l_f = 10.8  # longitud del fuselaje
b_f = 1.65  # ancho del fuselaje
h_f = 1.71  # alto fuselaje

#------------Parámetros del empenaje----------------
CLa_h = 4.05  # Pendiente sustentacion empenaje
S_h = 5.448   # Superficie empenaje
n_h = 0.95    # rendimiento empenaje

#--------Parámetros de la planta motora--------------
#PARAMETROS PLANTA MOTORA
l_n=2.74     # Parametro geometrico l_n
l_fn=2.83    # Parametro geometrico l_fn
Sfront_barqui = 0.55
Swet_barqui = 8 
b_n = 0.755  # parametro geometrico Bn
l_p = 2.91   # Parametro geometrico l_p
l_h = 5.14   # parametro geometrico l_h
mb2 = 1.36   # Parametro geometrico mb2
N_p = 4 # número de palas
Nn = 2 # número de barquillas
dp = 2.286  # Diametro helice

#DEFINICION DE PARAMETROS
#PARAMETROS PLANTA ALAR (EQUIVALENTE)
Clat=6.3810
Clar=6.5855
Clmax = (1.6 + 1.8)/2
Cli = 0.3
Cla=(Clat + Clar)/2
CMO=-0.01573702756129474
CLaw= 5.172
CLOw= 0.35697719620
t_c= 0.18
Cr= 2.67 #cuerda raiz
tr= t_c*Cr #espesor cuerda raiz
clor= 0.1451

#PARAMETROS FUSELAJE
zf=1.71
Swet_f = 59.22
Vf = 16.21
#PARAMETROS EMPENAJE HORIZONTAL
Clah = 4.05
Sh = 5.44
rendimientoh = 0.85
Ah=5.25
sweep1_4h = 0
taperh = 0.6
t_c_h = 0.12
CAMh=1.27
de_da=0.38
sweep1_2h = sweep1_4h-(1/Ah)*((1-taperh)/(1+taperh))
#PARAMETRO DEL EMPENAJE VERTICAL
Sv=2.5
Av=1.52
sweep1_4v=np.radians(50)
taperv=0.45
t_c_v=0.12
CAMv=1.45
sweep1_2v = sweep1_4v-(1/Av)*((1-taperv)/(1+taperv))



#-------------Parametros extraidos de otros TPS---------
#Parámetros del ala:
CLa_w = 5.1726 # Pendiente de sustentación del ala ' Sale de Multhopp
de_da = 0.38 #0.54
Cm_ac_w =-0.01547 #sale de Multhopp
#CURVA CL(al_pha) del avión completo
K_I = (1 + 2.15 * (b_f / b)) * (S_net / S) + (mt.pi/2) * (b_f ** 2) / (CLa_w * S)
K_II = (1 + 0.7 * b_f/b) * S_net/S
CLa_wf = K_I * CLa_w # pendiente sust ala fuselaje
print('Pendiente de sustentación combinación ala-fuselaje: CLa_wf = ', CLa_wf)
CLa = CLa_wf * (1 + (x_cg- x_ac) * CAM / l_h) #Pendiente de sustentacion del avion completo
CL0 = (w) / (0.5 * rho_c * V_crucero**2 * S) # Coeficiente de sustentación en condiciones de crucero
# CL0 = 0.356977
al_pha_0 =-CL0 / CLa