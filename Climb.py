import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import g
IMPORTAR EXCEL
# ============================================================
# 1. CONSTANTES Y PARÁMETROS DE LA AERONAVE
# ============================================================


# --- Beechcraft 60 Duke ---
from Duke import *

W = M * g                    # [N]


# ============================================================
# 2. ATMÓSFERA ISA
# ============================================================

def atmosfera_isa(h):
    """
    Atmósfera ISA para troposfera (h < 11000 m).

    Entrada:
        h : altitud [m]

    Salida:
        rho : densidad [kg/m³]
        a   : velocidad del sonido [m/s]
    """

    T0 = 288.15
    p0 = 101325.0
    L = 0.0065
    R = 287.05
    gamma_air = 1.4

    T = T0 - L * h

    p = p0 * (T / T0) ** (g / (R * L))

    rho = p / (R * T)

    a = np.sqrt(gamma_air * R * T)

    return rho, a


# ============================================================
# 3. POTENCIA DE LOS MOTORES
# ============================================================
#
# IMPORTANTE:
# P_eje debe ser la potencia TOTAL de ambos motores.
#
# Estos números son SOLAMENTE DE EJEMPLO.
# Aquí debes colocar los valores que ya calculaste/obtuviste
# para el Duke.
# ============================================================

altitud_motor = np.arange(0, 10001, 1000)  # [m]
print(altitud_motor)


potencia_total_hp = []
for h in altitud_motor:

    potencia_total_hp.append(2 * Pot_eje_alt[h])
print(potencia_total_hp)

HP_TO_W = 745.7

potencia_total_W = np.array(potencia_total_hp) * HP_TO_W
    
def potencia_eje(h):
    """
    Interpolación de Potencia TOTAL disponible al eje de ambos motores [W]
    en función de la altitud.
    """

    return np.interp(
        h,
        altitud_motor,
        potencia_total_W
    )


# ============================================================
# 4. RENDIMIENTO DE HÉLICE
# ============================================================
#
# Puedes reemplazar esta función directamente por tus
# curvas de Solies η(V,h).
#
# Por ahora dejo una curva analítica de ejemplo.
# ============================================================

def eta_helice(h, V):
    """
    Rendimiento de hélice.

    h : [m]
    V : TAS [m/s]

    Devuelve eta_p [-]
    """

    # EJEMPLO PROVISORIO
    # Sustituir por tus curvas reales η(V,h)

    eta = (
        0.84
        - 0.00010 * (V - 70.0)**2
        - 0.000005 * h
    )

    return np.clip(eta, 0.0, 0.90)


# ============================================================
# 5. POLAR AERODINÁMICA
# ============================================================

def coeficiente_CD(CL):

    CD = CD0 + K * CL**2

    return CD


# ============================================================
# 6. CÁLCULO ITERATIVO PARA UNA h Y UNA V
# ============================================================

def calcular_trepada(h, V,
                     tolerancia=1e-8,
                     max_iter=200):
    """
    Implementa el algoritmo iterativo de la clase para
    una altitud h y una velocidad TAS V.

    Devuelve:
        ROC      [m/s]
        gamma    [rad]
        CL
        CD
        D        [N]
        P_req    [W]
        P_disp   [W]

    Si el punto de vuelo no es físicamente válido,
    devuelve NaN.
    """

    rho, a = atmosfera_isa(h)

    q = 0.5 * rho * V**2

    M = V / a

    # -----------------------------------
    # Paso 2 de la clase:
    # gamma inicial
    # -----------------------------------

    gamma = 0.0

    for _ in range(max_iter):

        # -----------------------------------
        # Paso 3:
        #
        #            W cos(gamma)
        # CL = -------------------------
        #       1/2 rho V² S
        # -----------------------------------

        CL = W * np.cos(gamma) / (q * S)

        # El avión no puede mantener esta
        # condición si supera CLmax
        if CL > CL_max:
            return (np.nan,) * 7

        # -----------------------------------
        # Paso 4: polar
        # -----------------------------------

        CD = coeficiente_CD(CL, M)

        # -----------------------------------
        # Paso 5: resistencia y potencia
        # requerida
        # -----------------------------------

        D = q * S * CD

        P_req = D * V

        # -----------------------------------
        # Potencia propulsiva disponible
        #
        # P_A = eta_p * P_eje
        # -----------------------------------

        P_eje = potencia_eje(h)

        eta_p = eta_helice(h, V)

        P_disp = eta_p * P_eje

        # -----------------------------------
        # Paso 6:
        #
        #            P_A - D V
        # sin(gamma) = ----------
        #                W V
        # -----------------------------------

        sin_gamma = (P_disp - P_req) / (W * V)

        # Protección numérica
        if abs(sin_gamma) > 1.0:
            return (np.nan,) * 7

        gamma_nuevo = np.arcsin(sin_gamma)

        # -----------------------------------
        # Convergencia
        # -----------------------------------

        if abs(gamma_nuevo - gamma) < tolerancia:
            gamma = gamma_nuevo
            break

        gamma = gamma_nuevo

    else:
        # No convergió
        return (np.nan,) * 7

    # -----------------------------------
    # Paso 7:
    #
    # ROC = V sin(gamma)
    # -----------------------------------

    ROC = V * np.sin(gamma)

    return ROC, gamma, CL, CD, D, P_req, P_disp


# ============================================================
# 7. BUSCAR VELOCIDAD DE ROC MÁXIMA A UNA ALTITUD
# ============================================================

def buscar_Vy(h,
              V_max=130.0,
              dV=0.1):
    """
    Busca por barrido la velocidad TAS que produce el
    máximo ROC a una determinada altitud.

    h     : [m]
    V_max : velocidad máxima a analizar [m/s]
    dV    : paso de velocidades [m/s]

    Devuelve:
        Vy
        ROC_max
        gamma
    """

    rho, _ = atmosfera_isa(h)

    # Velocidad de pérdida aproximada
    # para vuelo horizontal
    Vs = np.sqrt(
        2.0 * W /
        (rho * S * CL_max)
    )

    # Evitamos trabajar exactamente en pérdida
    V_min = 1.02 * Vs

    velocidades = np.arange(
        V_min,
        V_max + dV,
        dV
    )

    ROC = np.full_like(velocidades, np.nan)

    gamma = np.full_like(velocidades, np.nan)

    for i, V in enumerate(velocidades):

        resultado = calcular_trepada(h, V)

        ROC[i] = resultado[0]
        gamma[i] = resultado[1]

    # Si no hay ningún punto válido
    if np.all(np.isnan(ROC)):
        return np.nan, np.nan, np.nan

    i_max = np.nanargmax(ROC)

    Vy = velocidades[i_max]

    ROC_max = ROC[i_max]

    gamma_Vy = gamma[i_max]

    return Vy, ROC_max, gamma_Vy


# ============================================================
# 8. ENVOLVENTE CON ALTITUD
# ============================================================

h_inicial = 0.0
h_final = 9000.0
dh = 250.0

altitudes = np.arange(
    h_inicial,
    h_final + dh,
    dh
)

Vy = np.full_like(altitudes, np.nan)
ROC_max = np.full_like(altitudes, np.nan)
gamma_Vy = np.full_like(altitudes, np.nan)


for i, h in enumerate(altitudes):

    Vy[i], ROC_max[i], gamma_Vy[i] = buscar_Vy(h)

    print(
        f"h = {h:6.0f} m | "
        f"Vy = {Vy[i]:6.2f} m/s | "
        f"Vy = {Vy[i] * 1.94384:6.1f} kt | "
        f"ROCmax = {ROC_max[i]:6.2f} m/s | "
        f"{ROC_max[i] * 196.8504:7.0f} ft/min"
    )


# ============================================================
# 9. GRÁFICA Vy vs ALTITUD
# ============================================================

plt.figure(figsize=(8, 5))

plt.plot(
    Vy * 1.94384,
    altitudes,
    marker="o"
)

plt.xlabel("Velocidad de máxima tasa de ascenso Vy [kt TAS]")
plt.ylabel("Altitud [m]")
plt.title("Beechcraft 60 Duke - Vy en función de la altitud")

plt.grid(True)

plt.tight_layout()
plt.show()


# ============================================================
# 10. GRÁFICA ROC MÁXIMA vs ALTITUD
# ============================================================

plt.figure(figsize=(8, 5))

plt.plot(
    ROC_max * 196.8504,
    altitudes,
    marker="o"
)

plt.xlabel("Tasa de ascenso máxima [ft/min]")
plt.ylabel("Altitud [m]")
plt.title("Beechcraft 60 Duke - Tasa de ascenso máxima")

plt.grid(True)

plt.tight_layout()
plt.show()