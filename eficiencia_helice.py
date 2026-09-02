"""Eficiencia de hélice en función de velocidad y altitud.

Esta función es la forma continua (y algebraicamente equivalente) del modelo Solies (1994). Acepta escalares o arrays
de NumPy, con broadcasting entre velocidad y altitud.

Unidades:
    V: velocidad verdadera [m/s]
    h: altitud ISA [m]
    resultado: eficiencia adimensional [-]

Como usar:
    from eficiencia_helice import eficiencia_helice

    eta = eficiencia_helice(120.0, 3000.0)
    print(eta)  # 0.856770278...

o para arrays y broadcasting:
    eta = eficiencia_helice(
        V=[50.0, 100.0],
        h=2000.0,
    )
"""

from __future__ import annotations

import numpy as np


# Parámetros del modelo original
D = 1.8796                    # diámetro de hélice [m]
P = 380.0 * 745.699872        # potencia al eje [W]
ETA_0 = 0.88                  # eficiencia de perfil empírica [-]
RHO_0 = 1.225                 # densidad ISA al nivel del mar [kg/m³]
A = np.pi * D**2 / 4.0        # área del disco [m²]


def densidad_isa(h: float | np.ndarray) -> float | np.ndarray:
    """Densidad ISA según la misma ley empleada para generar la tabla."""
    h_arr = np.asarray(h, dtype=float)
    if np.any((h_arr < 0.0) | (h_arr > 10_000.0)):
        raise ValueError("h debe estar entre 0 y 10 000 m (dominio de la tabla).")

    rho = RHO_0 * (1.0 - 2.2558e-5 * h_arr) ** 4.2559
    return float(rho) if rho.ndim == 0 else rho


def eficiencia_helice(
    V: float | np.ndarray,
    h: float | np.ndarray,
) -> float | np.ndarray:
    """Devuelve la eficiencia eta(V, h) para V > 0 y 0 <= h <= 10 000.

    La solución positiva ``x = vi/V`` satisface

        x (1 + x)^2 = ETA_0 P / (2 rho A V^3),

    y la eficiencia resulta ``eta = ETA_0 / (1 + x)``. La raíz se evalúa
    directamente mediante Cardano, sin iteraciones ni dependencia de SciPy.

    Ejemplos
    --------
    >>> eficiencia_helice(120.0, 3000.0)
    0.8567702780653651
    >>> eficiencia_helice([50.0, 100.0], 2000.0)
    array([0.71289871, 0.84519247])
    """
    V_arr = np.asarray(V, dtype=float)
    h_arr = np.asarray(h, dtype=float)
    V_b, h_b = np.broadcast_arrays(V_arr, h_arr)

    if np.any(V_b <= 0.0):
        raise ValueError("V debe ser estrictamente positiva.")

    rho = np.asarray(densidad_isa(h_b))
    q = ETA_0 * P / (2.0 * rho * A * V_b**3)

    # Raíz real positiva de x^3 + 2*x^2 + x - q = 0.
    delta = q**2 / 4.0 + q / 27.0
    centro = q / 2.0 + 1.0 / 27.0
    x = (
        np.cbrt(centro + np.sqrt(delta))
        + np.cbrt(centro - np.sqrt(delta))
        - 2.0 / 3.0
    )
    eta = ETA_0 / (1.0 + x)

    return float(eta) if eta.ndim == 0 else eta


if __name__ == "__main__":
    print(f"eta(120 m/s, 3000 m) = {eficiencia_helice(120.0, 3000.0):.8f}")
