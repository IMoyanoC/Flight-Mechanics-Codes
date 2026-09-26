"""Modelo digitalizado del Lycoming TIO-541-E a 2900 RPM.

Fuente: gráfico provisto por el usuario, "Lycoming Aircraft Engine Performance
Data", TIO-541-E Series, altitude performance, zero ram, mezcla 0.690 lb/BHP/h.

El gráfico da potencia AL FRENO/AL EJE de UN motor en función de la altitud de
presión y de la presión absoluta seca de admisión (MAP).  Los datos de este
módulo son una lectura aproximada del gráfico original, no datos tabulados del
fabricante. La tolerancia indicada por el propio gráfico es +/-2 %.

Uso rápido desde otro archivo::

    from motor_tio541_2900rpm import potencia_motor, potencia_total_duke

    p_un_motor = potencia_motor(15000, 40)       # hp, h en m
    p_dos_motores = potencia_total_duke(15000, 40)

Para trabajar en pies::

    p_un_motor = potencia_motor(4500, 40, unidad_altitud="ft")

La interpolación PCHIP reproduce toda la curva dibujada, incluido el tramo de
caída de potencia próximo a la altitud crítica. La función
``potencia_cuadratica`` es una simplificación analítica y sólo se considera
válida entre 0 y 22000 ft.

Para extender exclusivamente la curva de 42 inHg por encima de 27 000 ft debe
activarse explícitamente ``extrapolar_42=True``. La extrapolación conserva la
pendiente del último tramo digitalizado (26--27 kft), por lo que constituye una
hipótesis de cálculo y no información respaldada por la carta del fabricante.
"""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import numpy as np

try:
    from scipy.interpolate import PchipInterpolator
except ImportError:  # El modelo sigue funcionando mediante interpolación lineal.
    PchipInterpolator = None


RPM = 2900
# Filas: MAP = 34, 36, 38, 40 y 42 inHg.
MAPS_INHG = np.array([34.0, 36.0, 38.0, 40.0, 42.0]) 
# Columnas: altitud de presión = 0, 1, ..., 27 miles de ft.
ALTITUD_KFT = np.arange(0.0, 28.0, 1.0)


# Potencia aproximada en hp para UN motor.
POTENCIA_HP = np.array(
    [
        [
            322, 324, 326, 328, 329, 330, 331, 331, 331, 330, 329, 328,
            326, 324, 322, 320, 317, 314, 311, 308, 304, 300, 296, 291,
            287, 282, 279, 275,
        ],
        [
            337, 340, 342, 344, 346, 347, 348, 348, 348, 347, 345, 343,
            341, 339, 336, 334, 331, 328, 325, 322, 318, 314, 309, 304,
            298, 292, 284, 275,
        ],
        [
            353, 356, 359, 361, 363, 364, 365, 365, 364, 363, 362, 360,
            358, 356, 353, 350, 347, 344, 340, 336, 332, 327, 322, 316,
            309, 301, 289, 275,
        ],
        [
            370, 373, 375, 377, 379, 380, 380, 380, 380, 379, 378, 376,
            374, 372, 370, 367, 364, 361, 357, 353, 349, 344, 339, 333,
            325, 313, 292, 275,
        ],
        [
            389, 391, 392, 393, 394, 394, 394, 394, 393, 392, 391, 390,
            388, 386, 384, 382, 379, 376, 372, 368, 363, 358, 353, 346,
            335, 318, 294, 275,
        ],
    ],
    dtype=float,
)

# Hipótesis para h > 27 kft: continuación de la recta que une los dos últimos
# puntos de la curva de 42 inHg: (26 kft, 294 hp) y (27 kft, 275 hp).
PENDIENTE_EXTRAPOLACION_42_HP_POR_KFT = (
    POTENCIA_HP[-1, -1] - POTENCIA_HP[-1, -2]
) / (ALTITUD_KFT[-1] - ALTITUD_KFT[-2])
ALTITUD_POTENCIA_NULA_42_KFT = ALTITUD_KFT[-1] - (
    POTENCIA_HP[-1, -1] / PENDIENTE_EXTRAPOLACION_42_HP_POR_KFT
)

# P(h) = a*h^2 + b*h + c, con h en miles de ft y P en hp.
# Ajustes por mínimos cuadrados de los puntos entre 0 y 22 kft.
COEFICIENTES_CUADRATICOS = {
    34.0: np.array([-0.1565217391, 2.1964426877, 322.5782608696]),
    36.0: np.array([-0.1667419537, 2.2671372106, 338.7478260870]),
    38.0: np.array([-0.1869565217, 2.5754940711, 354.6913043478]),
    40.0: np.array([-0.1789102202, 2.4488706945, 371.1478260870]),
    42.0: np.array([-0.1548277809, 1.8113495200, 388.9260869565]),
}

AJUSTE_CUADRATICO_RMSE_HP = {
    34.0: 0.524,
    36.0: 1.135,
    38.0: 0.905,
    40.0: 0.580,
    42.0: 0.497,
}


def _altitud_a_kft(
    altitud: float | np.ndarray,
    unidad_altitud: Literal["ft", "kft", "m"],
) -> np.ndarray:
    """Convierte una altitud escalar o vectorial a miles de pies."""
    h = np.asarray(altitud, dtype=float)
    if unidad_altitud == "ft":
        return h / 1000.0
    if unidad_altitud == "kft":
        return h
    if unidad_altitud == "m":
        return h * 3.280839895 / 1000.0
    raise ValueError("unidad_altitud debe ser 'ft', 'kft' o 'm'.")


def _validar_dominio(
    h_kft: np.ndarray,
    map_inhg: np.ndarray,
    extrapolar_42: bool,
) -> None:
    if np.any(h_kft < ALTITUD_KFT[0]):
        raise ValueError("La altitud no puede ser negativa.")
    if np.any((map_inhg < MAPS_INHG[0]) | (map_inhg > MAPS_INHG[-1])):
        raise ValueError("La digitalización sólo cubre MAP entre 34 y 42 inHg.")

    fuera_de_carta = h_kft > ALTITUD_KFT[-1]
    if np.any(fuera_de_carta) and not extrapolar_42:
        raise ValueError(
            "La carta termina en 27 000 ft. Para continuar linealmente la curva "
            "de 42 inHg usa extrapolar_42=True."
        )
    if np.any(fuera_de_carta & ~np.isclose(map_inhg, 42.0)):
        raise ValueError(
            "Por encima de 27 000 ft sólo se implementó la extrapolación de MAP=42 inHg."
        )
    if np.any(h_kft > ALTITUD_POTENCIA_NULA_42_KFT):
        raise ValueError(
            "La extrapolación lineal alcanzaría potencia negativa por encima de "
            f"{ALTITUD_POTENCIA_NULA_42_KFT:.2f} kft."
        )


def _interpolar_entre_maps(
    valores_por_map: np.ndarray,
    map_objetivo: np.ndarray,
) -> np.ndarray:
    """Interpola linealmente en MAP para cada punto solicitado."""
    map_plano = map_objetivo.ravel()
    j = np.searchsorted(MAPS_INHG, map_plano, side="right") - 1
    j = np.clip(j, 0, len(MAPS_INHG) - 2)

    map_inf = MAPS_INHG[j]
    map_sup = MAPS_INHG[j + 1]
    fraccion = (map_plano - map_inf) / (map_sup - map_inf)

    columnas = np.arange(map_plano.size)
    p_inf = valores_por_map[j, columnas]
    p_sup = valores_por_map[j + 1, columnas]
    return (p_inf + fraccion * (p_sup - p_inf)).reshape(map_objetivo.shape)


def potencia_motor(
    altitud: float | np.ndarray,
    map_inhg: float | np.ndarray,
    *,
    unidad_altitud: Literal["ft", "kft", "m"] = "m",
    metodo: Literal["pchip", "lineal"] = "pchip",
    extrapolar_42: bool = True,
) -> float | np.ndarray:
    """Devuelve la potencia aproximada de UN motor, en hp.

    Los argumentos ``altitud`` y ``map_inhg`` pueden ser escalares o arrays
    compatibles mediante broadcasting. Primero se interpola cada curva en
    altitud y luego linealmente entre las curvas de MAP.

    Parameters
    ----------
    altitud:
        Altitud de presión; la unidad se selecciona con ``unidad_altitud``.
    map_inhg:
        Presión absoluta seca de admisión, entre 34 y 42 inHg.
    unidad_altitud:
        ``"ft"`` (por defecto), ``"kft"`` o ``"m"``.
    metodo:
        ``"pchip"`` para curva suave y sin oscilaciones; ``"lineal"`` para
        interpolar directamente entre los puntos digitalizados.
    extrapolar_42:
        Si es ``True``, permite continuar la curva de 42 inHg por encima de
        27 000 ft con una recta de pendiente -19 hp/kft. Esta opción sólo se
        admite para MAP=42 inHg y mientras la potencia calculada sea positiva.
    """
    h_kft, map_array = np.broadcast_arrays(
        _altitud_a_kft(altitud, unidad_altitud),
        np.asarray(map_inhg, dtype=float),
    )
    _validar_dominio(h_kft, map_array, extrapolar_42)
    h_plano = h_kft.ravel()
    h_dentro_carta = np.minimum(h_plano, ALTITUD_KFT[-1])

    if metodo == "pchip":
        if PchipInterpolator is None:
            raise ImportError(
                "El método 'pchip' requiere scipy. Usa metodo='lineal' o instala scipy."
            )
        valores_por_map = np.vstack(
            [
                PchipInterpolator(ALTITUD_KFT, fila)(h_dentro_carta)
                for fila in POTENCIA_HP
            ]
        )
    elif metodo == "lineal":
        valores_por_map = np.vstack(
            [np.interp(h_dentro_carta, ALTITUD_KFT, fila) for fila in POTENCIA_HP]
        )
    else:
        raise ValueError("metodo debe ser 'pchip' o 'lineal'.")

    resultado = _interpolar_entre_maps(valores_por_map, map_array)
    fuera_de_carta = h_kft > ALTITUD_KFT[-1]
    if np.any(fuera_de_carta):
        resultado = np.asarray(resultado)
        resultado[fuera_de_carta] = POTENCIA_HP[-1, -1] + (
            PENDIENTE_EXTRAPOLACION_42_HP_POR_KFT
            * (h_kft[fuera_de_carta] - ALTITUD_KFT[-1])
        )
    return float(resultado) if resultado.ndim == 0 else resultado


def potencia_total_duke(
    altitud: float | np.ndarray,
    map_inhg: float | np.ndarray,
    *,
    unidad_altitud: Literal["ft", "kft", "m"] = "m",
    metodo: Literal["pchip", "lineal"] = "pchip",
    extrapolar_42: bool = True,
) -> float | np.ndarray:
    """Potencia total en el eje de los DOS motores del Beechcraft Duke, en hp."""
    return 2.0 * potencia_motor(
        altitud,
        map_inhg,
        unidad_altitud=unidad_altitud,
        metodo=metodo,
        extrapolar_42=extrapolar_42,
    )


def potencia_cuadratica(
    altitud: float | np.ndarray,
    map_inhg: float | np.ndarray,
    *,
    unidad_altitud: Literal["ft", "kft", "m"] = "ft",
) -> float | np.ndarray:
    """Aproximación cuadrática, válida únicamente entre 0 y 22 000 ft.

    Para MAP intermedios se interpolan linealmente los coeficientes de las
    parábolas vecinas. No usar esta función para reproducir la caída pronunciada
    del extremo derecho del gráfico; para eso se debe usar ``potencia_motor``.
    """
    h_kft, map_array = np.broadcast_arrays(
        _altitud_a_kft(altitud, unidad_altitud),
        np.asarray(map_inhg, dtype=float),
    )
    if np.any((h_kft < 0.0) | (h_kft > 22.0)):
        raise ValueError("La aproximación cuadrática sólo es válida entre 0 y 22 000 ft.")
    if np.any((map_array < 34.0) | (map_array > 42.0)):
        raise ValueError("La aproximación sólo cubre MAP entre 34 y 42 inHg.")

    coef = np.vstack([COEFICIENTES_CUADRATICOS[m] for m in MAPS_INHG])
    h_plano = h_kft.ravel()
    valores_por_map = np.vstack(
        [a * h_plano**2 + b * h_plano + c for a, b, c in coef]
    )
    resultado = _interpolar_entre_maps(valores_por_map, map_array)
    return float(resultado) if resultado.ndim == 0 else resultado


def graficar_ajustes(
    archivo: str | Path | None = None,
    *,
    mostrar_cuadraticas: bool = True,
    mostrar: bool = True,
):
    """Grafica puntos digitalizados, PCHIP y, opcionalmente, las cuadráticas.

    Si ``archivo`` termina en ``.png`` o ``.pdf``, la figura se exporta en ese
    formato. La función devuelve ``(fig, ax)`` para poder modificar el gráfico.
    """
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    fig, ax = plt.subplots(figsize=(9.0, 5.8))
    h_fina = np.linspace(0.0, 27.0, 600)
    h_cuad = np.linspace(0.0, 22.0, 400)
    colores = plt.get_cmap("viridis")(np.linspace(0.08, 0.92, len(MAPS_INHG)))

    for i, (map_i, color) in enumerate(zip(MAPS_INHG, colores)):
        ax.scatter(
            ALTITUD_KFT,
            POTENCIA_HP[i],
            s=15,
            color=color,
            zorder=3,
        )
        ax.plot(
            h_fina,
            potencia_motor(h_fina, map_i, unidad_altitud="kft"),
            color=color,
            linewidth=1.8,
        )
        if mostrar_cuadraticas:
            ax.plot(
                h_cuad,
                potencia_cuadratica(h_cuad, map_i, unidad_altitud="kft"),
                color=color,
                linestyle="--",
                linewidth=1.0,
                alpha=0.9,
            )

    ax.set(
        xlabel="Altitud de presión, h [10³ ft]",
        ylabel="Potencia al eje de un motor [hp]",
        title="Lycoming TIO-541-E — familia de curvas a 2900 RPM",
        xlim=(0, 27),
        ylim=(270, 405),
    )
    ax.grid(True, alpha=0.3)

    leyenda_map = [
        Line2D([0], [0], color=color, linewidth=2, label=f"MAP={map_i:.0f} inHg")
        for map_i, color in zip(MAPS_INHG, colores)
    ]
    leyenda_estilo = [
        Line2D([0], [0], marker="o", linestyle="none", color="black", label="Puntos digitalizados"),
        Line2D([0], [0], linestyle="-", color="black", label="Interpolación PCHIP"),
    ]
    if mostrar_cuadraticas:
        leyenda_estilo.append(
            Line2D(
                [0], [0], linestyle="--", color="black",
                label="Ajuste cuadrático (0–22 kft)",
            )
        )
    primera_leyenda = ax.legend(handles=leyenda_map, fontsize=8, ncol=3, loc="upper right")
    ax.add_artist(primera_leyenda)
    ax.legend(handles=leyenda_estilo, fontsize=8, loc="lower left")
    fig.tight_layout()

    if archivo is not None:
        fig.savefig(archivo, dpi=300, bbox_inches="tight")
    if mostrar:
        plt.show()
    return fig, ax


if __name__ == "__main__":
    print(f"Un motor:  h=15 000 ft, MAP=40 inHg -> {potencia_motor(5000, 42):.1f} hp")
    print(f"Dos motores: h=15 000 ft, MAP=40 inHg -> {potencia_total_duke(5000, 42):.1f} hp")
    print(
        "Un motor extrapolado: h=30 000 ft, MAP=42 inHg -> "
        f"{potencia_motor(3000, 42, extrapolar_42=True):.1f} hp"
    )
    graficar_ajustes("curvas_tio541_2900rpm.png", mostrar=False)
    #graficar_ajustes("curvas_tio541_2900rpm.pdf", mostrar=False)

#Test
print (potencia_motor(7000, 42,unidad_altitud="ft"))