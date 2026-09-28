import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
df = pd.read_excel("tabla_eficiencia_altitud.xlsx")
from Duke import Pot_eje_alt,CL_max,K,CD0
# ============================================================================
# AIRCRAFT DATA
# ============================================================================

rho0 = 1.225                 # Sea-level density [kg/m^3]

W = 30146                     # Aircraft weight [N]
S = 19.78                      # Wing area [m^2]

CLmax = CL_max                  # Maximum lift coefficient

w = W / S                    # Wing loading [N/m^2]

HP_TO_W = 745.699872


# ============================================================================
# ALTITUDES
# ============================================================================

# These MUST coincide with the altitudes used in tabla_eficiencia_altitud

altitude = np.arange(0, 10001, 500)
print("Altitudes:", altitude, "m")
print("Number of altitudes:", len(altitude))
# ============================================================================
# ENGINE SHAFT POWER
# ============================================================================

# Shaft power developed by ONE engine at each altitude.
# Constant with flight speed at each altitude.
#
# Example:
# P_engine_hp = np.array([
#     620,
#     610,
#     595,
#     ...
# ])

P_engine_hp = Pot_eje_alt

print(len(P_engine_hp), "engine power values provided for altitudes:", altitude, "m")
P_engine = P_engine_hp * HP_TO_W 


if len(P_engine) != len(altitude):
    raise ValueError(
        "P_engine_hp must contain one value for each altitude."
    )


# ============================================================================
# ATMOSPHERE
# ============================================================================

def density_ISA(h):

    return rho0 * (1 - 2.2558e-5 * h)**4.2559


# ============================================================================
# MINIMUM REQUIRED POWER - EQ. (28)
# ============================================================================

def minimum_required_power(rho):

    sigma = rho / rho0

    Pmin_dimensionless = (
        np.sqrt(2 / sigma)
        * (4 / 3)
        * (3 * CD0 * K**3)**0.25
    )

    Pmin = (
        W
        * np.sqrt(w / rho0)
        * Pmin_dimensionless
    )

    return Pmin


# ============================================================================
# ITERATIVE SOLUTION - EQ. (32)
# ============================================================================

def solve_velocity(
    v0,
    V_MD,
    P_engine_h,
    P_required_min,
    V_eta,
    eta_values,
    tol=1e-8,
    max_iter=100
):

    v = v0

    for _ in range(max_iter):

        # ----------------------------------------------------
        # Dimensionless velocity -> true airspeed
        # ----------------------------------------------------

        V = v * V_MD


        # ----------------------------------------------------
        # Propeller efficiency from Solies table
        # ----------------------------------------------------

        eta_p = np.interp(
            V,
            V_eta,
            eta_values
        )


        # ----------------------------------------------------
        # Available propulsive power
        # ----------------------------------------------------

        P_available = eta_p * P_engine_h


        # ----------------------------------------------------
        # epsilon_P
        # ----------------------------------------------------

        epsilon_P = (
            P_available
            / P_required_min
        )


        # ----------------------------------------------------
        # Equation (32)
        # ----------------------------------------------------

        numerator = (
            v**4
            - (4 / 3**0.75) * epsilon_P * v
            + 1
        )

        denominator = (
            4 * v**3
            - (4 / 3**0.75) * epsilon_P
        )


        # Avoid division by zero
        if abs(denominator) < 1e-12:
            return np.nan


        v_new = (
            v
            - numerator / denominator
        )


        # Physically impossible solution
        if v_new <= 0:
            return np.nan


        # ----------------------------------------------------
        # Convergence
        # ----------------------------------------------------

        if abs(v_new - v) < tol:

            return v_new


        v = v_new


    return np.nan


# ============================================================================
# FLIGHT ENVELOPE
# ============================================================================

V_min_power = []
V_max = []
V_stall = []
V_min_real = []

valid_altitude = []


for h, Pe in zip(altitude, P_engine):


    # ========================================================================
    # ATMOSPHERIC CONDITIONS
    # ========================================================================

    rho = density_ISA(h)


    # ========================================================================
    # MINIMUM REQUIRED POWER
    # ========================================================================

    P_required_min = minimum_required_power(rho)


    # ========================================================================
    # REFERENCE VELOCITY
    # ========================================================================

    CL_MD = np.sqrt(CD0 / K)

    V_MD = np.sqrt(
        2 * W
        /
        (rho * S * CL_MD)
    )


    # Speed corresponding to minimum required power
    V_MP = V_MD / 3**0.25


    # ========================================================================
    # GET SOLIES ETA(V) FOR THIS ALTITUDE
    # ========================================================================

    column_name = f"{h:.0f} m"

    V_eta = df[
        (df["Altitude"] == column_name) & (df["Variable"] == "V [m/s]")
    ]["Value"].to_numpy()

    eta_values = df[
        (df["Altitude"] == column_name) & (df["Variable"] == "eta")
    ]["Value"].to_numpy()


    # ========================================================================
    # INITIAL epsilon_P
    # ========================================================================

    # Evaluate propeller efficiency near minimum-power speed

    eta_initial = np.interp(
        V_MP,
        V_eta,
        eta_values
    )

    P_available_initial = eta_initial * Pe

    epsilon_initial = (
        P_available_initial
        / P_required_min
    )


    # ========================================================================
    # CHECK AVAILABLE POWER
    # ========================================================================

    if epsilon_initial <= 0:

        continue


    # ========================================================================
    # INITIAL VALUES - EQ. (33)
    # ========================================================================

    v0_max = (
        (4 * epsilon_initial)**(1/3)
        /
        3**0.25
    )

    v0_min = (
        3**0.75
        /
        (4 * epsilon_initial)
    )


    # ========================================================================
    # MAXIMUM SPEED
    # ========================================================================

    v_max = solve_velocity(
        v0_max,
        V_MD,
        Pe,
        P_required_min,
        V_eta,
        eta_values
    )


    # ========================================================================
    # MINIMUM SPEED DUE TO POWER
    # ========================================================================

    v_min = solve_velocity(
        v0_min,
        V_MD,
        Pe,
        P_required_min,
        V_eta,
        eta_values
    )


    # ========================================================================
    # CONVERT TO TRUE AIRSPEED
    # ========================================================================

    if np.isnan(v_max) or np.isnan(v_min):
        continue


    Vmax_h = v_max * V_MD
    Vmin_power_h = v_min * V_MD


    # ========================================================================
    # STALL SPEED
    # ========================================================================

    Vs_h = np.sqrt(
        2 * W
        /
        (rho * S * CLmax)
    )


    # ========================================================================
    # ACTUAL MINIMUM SPEED
    # ========================================================================

    Vmin_real_h = max(
        Vmin_power_h,
        Vs_h
    )


    # No level-flight envelope if stall exceeds maximum speed
    if Vmin_real_h >= Vmax_h:
        continue


    # ========================================================================
    # STORE RESULTS
    # ========================================================================

    valid_altitude.append(h)

    V_min_power.append(
        Vmin_power_h
    )

    V_max.append(
        Vmax_h
    )

    V_stall.append(
        Vs_h
    )

    V_min_real.append(
        Vmin_real_h
    )


# ============================================================================
# ARRAYS
# ============================================================================

valid_altitude = np.array(valid_altitude)

V_min_power = np.array(V_min_power)
V_max = np.array(V_max)

V_stall = np.array(V_stall)
V_min_real = np.array(V_min_real)


# ============================================================================
# RESULTS
# ============================================================================

print("\n===== FLIGHT ENVELOPE =====")

for h, vs, vpmin, vmin, vmax in zip(
    valid_altitude,
    V_stall,
    V_min_power,
    V_min_real,
    V_max
):

    print(
        f"h = {h:5.0f} m | "
        f"Vs = {vs:6.2f} m/s | "
        f"Vmin(P) = {vpmin:6.2f} m/s | "
        f"Vmin(real) = {vmin:6.2f} m/s | "
        f"Vmax = {vmax:6.2f} m/s"
    )


# ============================================================================
# FLIGHT ENVELOPE PLOT
# ============================================================================

plt.figure(figsize=(9, 8))


# Minimum speed due to power
plt.plot(
    V_min_power,
    valid_altitude,
    "--",
    linewidth=1.5,
    label="V mínima por potencia"
)


# Stall speed
plt.plot(
    V_stall,
    valid_altitude,
    linewidth=2,
    label="Velocidad de pérdida $V_s$"
)


# Actual minimum envelope
plt.plot(
    V_min_real,
    valid_altitude,
    linewidth=2.5,
    label="Límite inferior real"
)


# Maximum speed
plt.plot(
    V_max,
    valid_altitude,
    linewidth=2.5,
    label="Velocidad máxima"
)


# Flight region
plt.fill_betweenx(
    valid_altitude,
    V_min_real,
    V_max,
    alpha=0.15
)


plt.xlabel("Velocidad verdadera [m/s]")
plt.ylabel("Altitud [m]")

plt.title(
    "Envolvente de vuelo"
)

plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()