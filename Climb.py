import numpy as np
import matplotlib.pyplot as plt

from scipy.constants import g

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
from motor_tio541_2900rpm import potencia_total_duke

HP_TO_W = 745.7
    
def potencia_eje(h):
    """
    Cálculo de la potencia [W]TOTAL disponible (2 motores) 
    en el eje para una altitud h [m], con MAP 40 y 
    """
    return potencia_total_duke(
        h,
        40,
    ) * HP_TO_W


# ============================================================
# 4. RENDIMIENTO DE HÉLICE
# ============================================================
from eficiencia_helice import eficiencia_helice

def eta_helice(h, V):

    eta = eficiencia_helice(V, h)

    return eta


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
    Implementa el algoritmo iterativo para
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

        CD = coeficiente_CD(CL)

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
    V_min = 1.001 * Vs

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
h_final = 10000.0
dh = 1000.0

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

# ============================================================
# CURVAS DE PERFORMANCE DE ASCENSO
# ROC(V) y gamma(V) para diferentes altitudes
# ============================================================

rho0, _ = atmosfera_isa(0.0)
altitudes_grafico = [
    0,
    3000,
    6000,
    9000
]

def TAS_a_EAS(V_tas, h):
    """
    Convierte True Airspeed a Equivalent Airspeed.

    V_E = V_TAS * sqrt(rho/rho0)

    Entradas:
        V_tas : [m/s]
        h     : [m]

    Salida:
        V_eas : [m/s]
    """

    rho, _ = atmosfera_isa(h)

    return V_tas * np.sqrt(rho / rho0)


def EAS_a_TAS(V_eas, h):
    """
    Convierte Equivalent Airspeed a True Airspeed.

    V_TAS = V_E / sqrt(rho/rho0)
    """

    rho, _ = atmosfera_isa(h)

    return V_eas / np.sqrt(rho / rho0)


def calcular_curva_ascenso(
    h,
    Ve_max=150.0,
    dVe=1,
    solo_ascenso=True
):
    """
    Calcula ROC y gamma en función de la velocidad equivalente
    para una altitud determinada.

    Parámetros
    ----------
    h : float
        Altitud [m]

    Ve_max : float
        Máxima velocidad equivalente a analizar [m/s]

    dVe : float
        Paso de velocidades equivalentes [m/s]

    solo_ascenso : bool
        Si True, elimina los puntos con ROC < 0.
        Esto hace que las curvas terminen cuando ROC = 0,
        como en los gráficos de clase.

    Devuelve
    --------
    Ve : ndarray
        Velocidad equivalente [m/s]

    Vtas : ndarray
        Velocidad verdadera [m/s]

    ROC : ndarray
        Tasa de ascenso [m/s]

    gamma : ndarray
        Ángulo de ascenso [rad]
    """

    rho, _ = atmosfera_isa(h)

    # --------------------------------------------------------
    # Velocidad de pérdida equivalente
    #
    # Al trabajar en EAS:
    #
    # Vs_E = sqrt(2 W / (rho0 S CLmax))
    #
    # Es prácticamente independiente de la altitud.
    # --------------------------------------------------------

    Vs_eas = np.sqrt(
        2.0 * W /
        (rho0 * S * CL_max)
    )

    # Empezamos ligeramente por encima de pérdida
    Ve_min = 1.02 * Vs_eas

    Ve = np.arange(
        Ve_min,
        Ve_max + dVe,
        dVe
    )

    ROC = np.full_like(Ve, np.nan, dtype=float)
    gamma = np.full_like(Ve, np.nan, dtype=float)
    Vtas = np.full_like(Ve, np.nan, dtype=float)

    for i, Ve_i in enumerate(Ve):

        # La aerodinámica debe calcularse con TAS
        Vtas_i = EAS_a_TAS(Ve_i, h)

        Vtas[i] = Vtas_i

        resultado = calcular_trepada(
            h,
            Vtas_i
        )

        ROC_i = resultado[0]
        gamma_i = resultado[1]

        # Si el punto no es válido, se deja NaN
        if np.isnan(ROC_i):
            continue

        # Podemos eliminar la región de descenso
        # para reproducir el estilo de la clase.
        if solo_ascenso and ROC_i < 0:
            continue

        ROC[i] = ROC_i
        gamma[i] = gamma_i

    return Ve, Vtas, ROC, gamma



# ============================================================
# GRÁFICAS ROC(Ve) Y gamma(Ve)
# ============================================================

def graficar_performance_ascenso(
    altitudes_grafico,
    Ve_max=150.0,
    dVe=1,
    roc_unidad="m/min",
    mostrar_maximos=True
):

    """
    Genera:

        1) ROC vs velocidad equivalente
        2) gamma vs velocidad equivalente

    para varias altitudes.

    Parámetros
    ----------
    altitudes_grafico : iterable
        Altitudes a representar [m]

        Ejemplo:
            [0, 1500, 3000, 4500]

    Ve_max : float
        Máxima EAS del gráfico [m/s]

    dVe : float
        Paso de velocidad [m/s]

    roc_unidad : str
        "m/s", "m/min" o "ft/min"

    mostrar_maximos : bool
        Marca sobre cada curva el punto de ROC máxima.
    """

    # Velocidad de pérdida equivalente
    Vs_eas = np.sqrt(
        2.0 * W /
        (rho0 * S * CL_max)
    )

    # ========================================================
    # FIGURA 1 - ROC vs Ve
    # ========================================================

    plt.figure(figsize=(8, 6))

    for h in altitudes_grafico:

        Ve, Vtas, ROC, gamma = calcular_curva_ascenso(
            h,
            Ve_max=Ve_max,
            dVe=dVe,
            solo_ascenso=True
        )

        # --------------------------------------------
        # Conversión de unidades para ROC
        # --------------------------------------------

        if roc_unidad == "m/s":

            ROC_plot = ROC
            ylabel = r"$ROC$ [m/s]"

        elif roc_unidad == "m/min":

            ROC_plot = ROC * 60.0
            ylabel = r"$ROC$ [m/min]"

        elif roc_unidad == "ft/min":

            ROC_plot = ROC * 196.8504
            ylabel = r"$ROC$ [ft/min]"

        else:
            raise ValueError(
                "roc_unidad debe ser "
                "'m/s', 'm/min' o 'ft/min'"
            )

        plt.plot(
            Ve,
            ROC_plot,
            label=f"{h:.0f} m"
        )

        # --------------------------------------------
        # Punto de máxima tasa de ascenso
        # --------------------------------------------

        if mostrar_maximos and not np.all(np.isnan(ROC)):

            i_max = np.nanargmax(ROC)

            plt.plot(
                Ve[i_max],
                ROC_plot[i_max],
                "o"
            )

    # Línea de pérdida equivalente
    plt.axvline(
        Vs_eas,
        color="black",
        linestyle="--",
        linewidth=1.0,
        label=r"$V_{S,E}$"
    )

    plt.xlabel(
        r"Velocidad equivalente $V_e$ [m/s]"
    )

    plt.ylabel(ylabel)

    plt.title(
        "Performance de ascenso - Tasa de ascenso"
    )

    plt.grid(
        True,
        alpha=0.3
    )

    plt.legend()

    plt.tight_layout()

    plt.show()


    # ========================================================
    # FIGURA 2 - gamma vs Ve
    # ========================================================

    plt.figure(figsize=(8, 6))

    for h in altitudes_grafico:

        Ve, Vtas, ROC, gamma = calcular_curva_ascenso(
            h,
            Ve_max=Ve_max,
            dVe=dVe,
            solo_ascenso=True
        )

        gamma_deg = np.degrees(gamma)

        plt.plot(
            Ve,
            gamma_deg,
            label=f"{h:.0f} m"
        )

        # --------------------------------------------
        # Marcar gamma correspondiente a ROC máxima
        # --------------------------------------------

        if mostrar_maximos and not np.all(np.isnan(ROC)):

            i_max = np.nanargmax(ROC)

            plt.plot(
                Ve[i_max],
                gamma_deg[i_max],
                "o"
            )

    plt.axvline(
        Vs_eas,
        color="black",
        linestyle="--",
        linewidth=1.0,
        label=r"$V_{S,E}$"
    )

    plt.xlabel(
        r"Velocidad equivalente $V_e$ [m/s]"
    )

    plt.ylabel(
        r"Ángulo de ascenso $\gamma$ [°]"
    )

    plt.title(
        "Performance de ascenso - Ángulo de ascenso"
    )

    plt.grid(
        True,
        alpha=0.3
    )

    plt.legend()

    plt.tight_layout()

    plt.show()

graficar_performance_ascenso(
    altitudes_grafico,
    Ve_max=150.0,
    dVe=1,
    roc_unidad="m/min",
    mostrar_maximos=True
)