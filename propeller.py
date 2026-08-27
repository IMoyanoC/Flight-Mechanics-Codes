import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import pandas as pd

# ============================================================================
# INPUTS
# ============================================================================

# Propeller geometry
D = 1.8796              # Propeller diameter [m]

# Engine / shaft data
RPM = 2900             # Propeller rotational speed [rev/min]
P_hp = 380             # Shaft power [hp]
Pc_hp = 585.0          # Shaft power in cruise condition [HP]

# Atmospheric data
rho0 = 1.225           # Air density [kg/m^3]              
rho_cruise = 0.905216    # Cruise air density (10.000 ft) [kg/m^3]

#Cruise conditions
V_cruise = 120       # Cruise speed [m/s]

# Empirical profile efficiency
eta0 = 0.88           #0.88 for consertive estimate, 0.92 for high-performance 

# Speed range
V_min = 5            # Minimum speed [m/s] (not zero)
V_max = 150.0          # Maximum speed [m/s]
n_delta = 5         # Number of points in speed vector


# ============================================================================
# UNIT CONVERSIONS
# ============================================================================

HP_TO_W = 745.699872
P = P_hp * HP_TO_W             # Power [W]
N = Pc_hp * HP_TO_W
n = RPM / 60.0                 # Rev/s
A = np.pi * D**2 / 4.0         # Propeller disk area [m²]

# Speed vector
V_array = np.arange(V_min, V_max, n_delta)


# ============================================================================
# SOLIES ITERATIVE SOLVER
# ============================================================================

def propeller_efficiency(V, P, rho, A, eta0, eta_guess=0.9, tol=1e-8,max_iter=100):

    eta = eta_guess

    for _ in range(max_iter):
        # Eq. (2)
        T = eta * P / V

        # Eq. (3)
        dV_half = -V / 2.0 + np.sqrt(V**2 / 4.0 + T / (2.0 * rho* A))

        # Eq. (4)
        eta_i = V / (V + dV_half)

        # Eq. (1)
        eta_new = eta0 * eta_i

        if abs(eta_new - eta) < tol:
            eta = eta_new
            break

        eta = eta_new

    # Final thrust
    T = eta * P / V

    return eta, T, eta_i


# ============================================================================
# CALCULATIONS
# ============================================================================

eta_list = []
T_list = []
J_list = []
J_listc = []
eta_listc = []
T_listc = []
for V in V_array:
    eta, T, eta_i = propeller_efficiency(V, P, rho0, A, eta0)

    # Advance ratio
    J = V / (n * D)

    eta_list.append(eta)
    T_list.append(T)
    J_list.append(J)
for V in V_array:
    eta, T, eta_i = propeller_efficiency(V, N, rho_cruise, A, eta0)

    # Advance ratio
    J = V / (n * D)

    eta_listc.append(eta)
    T_listc.append(T)
    J_listc.append(J)
eta_array = np.array(eta_list)
T_array = np.array(T_list)
eta_arrayc = np.array(eta_listc)
T_arrayc = np.array(T_listc)
J_arrayc = np.array(J_listc)
J_array = np.array(J_list)

eta_max_index = np.argmax(eta_array)
eta_max = eta_array[eta_max_index]

# ============================================================================
# LINEAL MOMENTUM THEORY
# ============================================================================
# Para la condición estática, el empuje se puede calcular usando la teoría de momento lineal:
T_staticCdM = (2 * rho0 * A)**(1/3) * (P)**(2/3)

# ============================================================
# RESOLUCIÓN POR TEORÍA DE CANTIDAD DE MOVIMIENTO
# ============================================================
eta_ideal = []
J_vals = []
for V in V_array:

    # Ecuación:
    # P = 2*rho*A*vi*(V + vi)^2

    def ecuacion(vi):
        return 2 * rho0 * A * vi * (V + vi)**2 - P

    # Valor inicial razonable
    vi_sol = fsolve(ecuacion, 5)[0]

    # Eficiencia ideal
    eta = V / (V + vi_sol)

    # Relación de avance
    J = V / (n * D)

    eta_ideal.append(eta)
    J_vals.append(J)

# ============================================================
# CONVERTIR A ARRAYS
# ============================================================

eta_ideal = np.array(eta_ideal)
J_vals = np.array(J_vals)


# ============================================================================
# PLOTS
# ============================================================================

plt.style.use("default")
fig, axes = plt.subplots(3, 1, figsize=(8, 18))


# --------------------------------------------------------------------------
# 1) η vs V
# --------------------------------------------------------------------------
axes[0].plot(V_array, eta_array, linewidth=2)
axes[0].set_xlabel("Velocidad V [m/s]")
axes[0].set_ylabel("Eficiencia η")
axes[0].set_title("Eficiencia de la hélice vs Velocidad")
axes[0].grid(True) 

# --------------------------------------------------------------------------
# 2) η vs J
# --------------------------------------------------------------------------
axes[1].plot(J_array, eta_array, linewidth=2)
axes[1].plot(J_vals, eta_ideal, linewidth=2, linestyle='--')
axes[1].set_xlabel("Relación de avance J")
axes[1].set_ylabel("Eficiencia η")
axes[1].set_title("Eficiencia de la hélice vs Relación de avance")
axes[1].legend(["Método Solies", "Teoría CdM"])
axes[1].grid(True)

# --------------------------------------------------------------------------
# 3) T vs V
# --------------------------------------------------------------------------
axes[2].plot(V_array, T_array, linewidth=2)
axes[2].plot(V_array, T_arrayc, linewidth=2)
axes[2].set_xlabel("Velocidad V [m/s]")
axes[2].set_ylabel("Empuje T [N]")
axes[2].set_title("Empuje vs Velocidad")
axes[2].legend(["Condición estática", "Condición de crucero"])
axes[2].grid(True)

plt.tight_layout()
plt.show()


print("===== PROPELLER PARAMETERS =====")
print(f"eta0                = {eta0:.2f}")
print(f"n                   = {n} rev/s")
print(f"eta_max             = {eta_max:.4f} at V = {V_array[eta_max_index]:.2f} m/s (J = {J_array[eta_max_index]:.3f})")
# ============================================================================
# STATIC THRUST
# ============================================================================
print("===== STATIC CONDITION =====")
T_static = (2.0 * rho0 * A)**(1/3) * (eta0 * P)**(2/3)  # Eq. (5)
print(f"Static thrust       = {T_static:.1f} N ({T_static*0.224809:.1f} lbf)")
print(f"Total static thrust = {T_static * 2:.1f} N ({T_static * 2 * 0.224809:.1f} lbf)")


# ============================================================================
# CRUISE THRUST
# ============================================================================

J_cruise = V_cruise / (n * D)                   # Advance ratio at cruise
Cp_cruise = N / (rho_cruise * n**3 * D**5)      # Power coefficient
eta_cruise, T_cruise, eta_i_cruise = propeller_efficiency(    V_cruise, N, rho_cruise, A, eta0)
Ct_cruise = eta_cruise * Cp_cruise / J_cruise   # Thrust coefficient

print("===== CRUISE CONDITION =====")
print(f"V_cruise            = {V_cruise:.2f} m/s ({V_cruise*1.94384:.1f} kt)")
print(f"Cruise J            = {J_cruise:.3f}")
print(f"Cruise Cp           = {Cp_cruise:.6f}")
print(f"Cruise Ct           = {Ct_cruise:.6f}") 
print(f"eta_cruise           = {eta_cruise:.4f}")
print(f"T_cruise (/motor)   = {T_cruise:.1f} N ({T_cruise*0.224809:.1f} lbf)")
print(f"T_Cruise (Total)    = {T_cruise * 2:.1f} N ({T_cruise * 2 * 0.224809:.1f} lbf)")

# ============================================================================
# EFFICIENCY VS ALTITUDE
# ============================================================================

#Altitude vector
altitude = np.arange(0, 10001, 500)  # Altitude [m]

table_data = {}

for h in altitude:

    # ISA density at altitude h
    rho_h = rho0 * (1 - 2.2558e-5 * h)**4.2559

    eta_h = []

    for V in V_array:
        eta, T, eta_i = propeller_efficiency(V, P, rho_h, A, eta0)
        eta_h.append(eta)

    table_data[(f"{h:.0f} m", "V [m/s]")] = V_array
    table_data[(f"{h:.0f} m", "eta")] = eta_h

eta_altitude_table = pd.DataFrame(table_data)

print("\n===== EFFICIENCY VS VELOCITY AND ALTITUDE =====")
print(eta_altitude_table.to_string(index=False))
tabla_excel = eta_altitude_table.copy()

tabla_excel.columns = [
    f"{altitud} - {variable}"
    for altitud, variable in tabla_excel.columns
]

tabla_excel.to_excel(
    "tabla_eficiencia_altitud.xlsx",
    index=False
)