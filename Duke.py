import numpy as np
#------Parametros Geometricos-----
#------Parametros Planta Alar-----
S = 19.780                  # superficie alar [m2]
S_net = 15.023              # Superficie expuesta [m2]
Swet_w= 23.368              # Superficie mojada ala [m2]
CAM = 1.850                 # Cuerda aerodinamica media [m]
CGM= 2.731                  # Cuerda media geometrica [m]
b = 11.970                  # Envergadura [m]
taper = 0.321               # Ahusamiento
A = 7.243                   # Alargamiento
sweep1_4= np.radians(0)     # Angulo sweep cuarto cuerda [grados]
clt_alf_2 = 6.5855          # Pendiente sustentacion (puntera)
cl_alf_1 = 6.3810           # Pendiente sustentacion (raiz)
Cm_ac_w = -0.016            # Coeficiente momento libre ala
alabeo = np.radians(-3.95)  # Alabeo aerodinamico total [grados]
delta_z_CL = 0.1            # Porque es ala baja

#-----------Parámetros condición vuelo crucero---------
M = 3073            # Masa aeronave [kg]
V_crucero = 113.18  # Velocidad crucero en m/s
hc = 8850           # Altura crucero 
x_cg = 0.3          # Posición centro de gravedad adimensionada con la CAM
x_ac = 0.248        # Posición centro aerodinámico adimensionada con la CAM


#------------Parámetros del fuselaje------------
l_f = 10.8    # longitud del fuselaje
b_f = 1.65    # ancho del fuselaje
h_f = 1.71    # alto fuselaje

#------------Parámetros del empenaje----------------
CLa_h = 4.05  # Pendiente sustentacion empenaje
S_h = 5.448   # Superficie empenaje
n_h = 0.95    # rendimiento empenaje

#--------Parámetros de la planta motora--------------
#PARAMETROS PLANTA MOTORA
l_n=2.74      # Parametro geometrico l_n
l_fn=2.83     # Parametro geometrico l_fn
Sfront_barqui = 0.55
Swet_barqui = 8 
b_n = 0.755   # parametro geometrico Bn
l_p = 2.91    # Parametro geometrico l_p
l_h = 5.14    # parametro geometrico l_h
mb2 = 1.36    # Parametro geometrico mb2
N_p = 4       # número de palas
Nn = 2        # número de barquillas
dp = 2.286    # Diametro helice

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
Cr= 2.67       # cuerda raiz
tr= t_c*Cr     # espesor cuerda raiz
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
de_da = 0.38 #
Cm_ac_w =-0.01547 #sale de Multhopp
