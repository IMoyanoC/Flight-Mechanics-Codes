"""
Distancia de despegue: clase + Roskam, capitulo 10.
Las interfaces locales de potencia y helice reproducen la 
forma minima de los modelos existentes para que este proyecto
pueda ejecutarse sin NumPy/SciPy.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import acos, asin, atan, cos, isfinite, pi, sin, sqrt, tan


G = 9.80665  # m/s2
R_AIR = 287.05287  # J/(kg K)
GAMMA_AIR = 1.4
T0_ISA = 288.15  # K
P0_ISA = 101_325.0  # Pa
LAPSE = 0.0065  # K/m
HP_TO_W = 745.699872
FT_TO_M = 0.3048


@dataclass(frozen=True)
class AirportConditions:
    altitude_m: float
    temperature_c: float
    headwind_mps: float
    runway_slope_percent: float
    rolling_friction: float
    runway_surface: str
    runway_available_m: float
    pressure_pa: float | None = None


@dataclass(frozen=True)
class Aircraft:
    mass_kg: float
    wing_area_m2: float
    cl_max_takeoff: float
    cd0: float
    induced_drag_k: float
    cl_ground: float
    cd_ground: float
    engine_map_inhg: float
    engine_count: int
    propeller_diameter_m: float
    propeller_profile_efficiency: float
    v_mc_mps: float | None = None


@dataclass(frozen=True)
class ModelOptions:
    ground_method: str = "general"
    intervals: int = 100
    rotation_time_s: float = 1.0
    vr_over_vs: float = 1.10
    vlof_over_vs: float = 1.15
    v50_over_vs: float = 1.20
    obstacle_height_m: float = 50.0 * FT_TO_M


@dataclass(frozen=True)
class Atmosphere:
    temperature_k: float
    pressure_pa: float
    pressure_altitude_m: float
    density_kgm3: float
    speed_of_sound_mps: float


def atmosphere_at_airport(
    altitude_m: float,
    temperature_c: float,
    pressure_pa: float | None = None,
) -> Atmosphere:
    """
    Troposfera minima: presion ISA a elevacion y densidad con OAT.
    Correspondencia futura: Flight-Mechanics-Codes/Climb.py -> atmosfera_isa.
    La version local admite temperatura real y presion indicada opcional.
    """
    if not -500.0 <= altitude_m < 11_000.0:
        raise ValueError("La altitud debe estar entre -500 y 11 000 m.")
    temperature_k = temperature_c + 273.15
    if temperature_k <= 0.0:
        raise ValueError("La temperatura absoluta debe ser positiva.")
    if pressure_pa is None:
        t_isa = T0_ISA - LAPSE * altitude_m
        pressure_pa = P0_ISA * (t_isa / T0_ISA) ** (G / (R_AIR * LAPSE))
    if pressure_pa <= 0.0:
        raise ValueError("La presion debe ser positiva.")
    rho = pressure_pa / (R_AIR * temperature_k)
    pressure_altitude_m = T0_ISA / LAPSE * (
        1.0 - (pressure_pa / P0_ISA) ** (R_AIR * LAPSE / G)
    )
    return Atmosphere(
        temperature_k=temperature_k,
        pressure_pa=pressure_pa,
        pressure_altitude_m=pressure_altitude_m,
        density_kgm3=rho,
        speed_of_sound_mps=sqrt(GAMMA_AIR * R_AIR * temperature_k),
    )



# Flight-Mechanics-Codes/motor_tio541_2900rpm.py -> potencia_total_duke.
# Coeficientes P(h)=a*h_kft^2+b*h_kft+c,(valida solamente entre 0 y 22 kft).
_ENGINE_POWER_COEFFICIENTS = {
    34.0: (-0.1565217391, 2.1964426877, 322.5782608696),
    36.0: (-0.1667419537, 2.2671372106, 338.7478260870),
    38.0: (-0.1869565217, 2.5754940711, 354.6913043478),
    40.0: (-0.1789102202, 2.4488706945, 371.1478260870),
    42.0: (-0.1548277809, 1.8113495200, 388.9260869565),
}


def shaft_power_available(
    altitude_m: float,
    map_inhg: float,
    engine_count: int,
) -> tuple[float, float]:
    """Potencia total al eje [W] e individual [hp], interpolando en MAP."""
    if engine_count <= 0:
        raise ValueError("El numero de motores debe ser positivo.")
    h_kft = altitude_m / FT_TO_M / 1000.0
    if not 0.0 <= h_kft <= 22.0:
        raise ValueError("La curva cuadratica local es valida entre 0 y 22 kft.")
    maps = sorted(_ENGINE_POWER_COEFFICIENTS)
    if not maps[0] <= map_inhg <= maps[-1]:
        raise ValueError("MAP debe estar entre 34 y 42 inHg.")
    upper = next(value for value in maps if value >= map_inhg)
    lower = next(value for value in reversed(maps) if value <= map_inhg)

    def evaluate(map_value: float) -> float:
        a, b, c = _ENGINE_POWER_COEFFICIENTS[map_value]
        return a * h_kft**2 + b * h_kft + c

    if upper == lower:
        hp_per_engine = evaluate(lower)
    else:
        fraction = (map_inhg - lower) / (upper - lower)
        hp_per_engine = evaluate(lower) + fraction * (
            evaluate(upper) - evaluate(lower)
        )
    return hp_per_engine * engine_count * HP_TO_W, hp_per_engine


def _real_cuberoot(value: float) -> float:
    return value ** (1.0 / 3.0) if value >= 0.0 else -(-value) ** (1.0 / 3.0)


def propeller_efficiency_and_thrust(
    airspeed_mps: float,
    density_kgm3: float,
    total_shaft_power_w: float,
    engine_count: int,
    diameter_m: float,
    profile_efficiency: float,
) -> tuple[float, float]:
    """Modelo Solies usado en Flight-Mechanics-Codes/eficiencia_helice.py.

    En V=0 se usa su limite de empuje estatico por cantidad de movimiento.
    Para viento de cola al inicio, la formulacion 1-D de clase admite V con
    signo; el modelo propulsivo se evalua con su magnitud.
    """
    if density_kgm3 <= 0.0 or total_shaft_power_w <= 0.0:
        raise ValueError("Densidad y potencia deben ser positivas.")
    if engine_count <= 0 or diameter_m <= 0.0:
        raise ValueError("Numero de motores y diametro deben ser positivos.")
    if not 0.0 < profile_efficiency <= 1.0:
        raise ValueError("La eficiencia de perfil debe estar en (0, 1].")
    speed = abs(airspeed_mps)
    disk_area = pi * diameter_m**2 / 4.0
    power_per_engine = total_shaft_power_w / engine_count
    if speed < 1.0e-9:
        thrust_per_engine = (2.0 * density_kgm3 * disk_area) ** (1.0 / 3.0) * (
            profile_efficiency * power_per_engine
        ) ** (2.0 / 3.0)
        return 0.0, engine_count * thrust_per_engine

    q = profile_efficiency * power_per_engine / (
        2.0 * density_kgm3 * disk_area * speed**3
    )
    delta = q**2 / 4.0 + q / 27.0
    center = q / 2.0 + 1.0 / 27.0
    x = (
        _real_cuberoot(center + sqrt(delta))
        + _real_cuberoot(center - sqrt(delta))
        - 2.0 / 3.0
    )
    eta = profile_efficiency / (1.0 + x)
    thrust = eta * total_shaft_power_w / speed
    return eta, thrust


def drag_coefficient(cl: float, aircraft: Aircraft) -> float:
    """Polar preliminar de Flight-Mechanics-Codes/Duke.py."""
    return aircraft.cd0 + aircraft.induced_drag_k * cl**2


def validate_inputs(
    airport: AirportConditions,
    aircraft: Aircraft,
    options: ModelOptions,
) -> None:
    if aircraft.mass_kg <= 0.0 or aircraft.wing_area_m2 <= 0.0:
        raise ValueError("Masa y superficie alar deben ser positivas.")
    if aircraft.cl_max_takeoff <= 0.0:
        raise ValueError("CL_max_TO debe ser positivo.")
    if not 0.0 <= airport.rolling_friction <= 0.30:
        raise ValueError("mu_g debe estar entre 0 y 0.30.")
    if abs(airport.headwind_mps) > 25.0:
        raise ValueError("El viento debe estar entre -25 y +25 m/s.")
    if airport.runway_available_m <= 0.0:
        raise ValueError("La longitud de pista disponible debe ser positiva.")
    if options.obstacle_height_m <= 0.0:
        raise ValueError("La altura de obstaculo debe ser positiva.")
    if options.rotation_time_s <= 0.0:
        raise ValueError("El tiempo de rotacion debe ser positivo.")
    if options.intervals <= 0:
        raise ValueError("El numero de intervalos debe ser positivo.")
    if options.ground_method not in {"general", "aproximado"}:
        raise ValueError("ground_method debe ser 'general' o 'aproximado'.")


def calculate_takeoff(
    airport: AirportConditions,
    aircraft: Aircraft,
    options: ModelOptions,
) -> dict[str, object]:
    """Calcula s_TO=s_NGR+s_R+s_TR+s_CL segun las diapositivas."""
    validate_inputs(airport, aircraft, options)
    atmosphere = atmosphere_at_airport(
        airport.altitude_m, airport.temperature_c, airport.pressure_pa
    )
    rho = atmosphere.density_kgm3
    weight_n = aircraft.mass_kg * G
    wing_loading = weight_n / aircraft.wing_area_m2
    vs = sqrt(
        2.0 * weight_n
        / (rho * aircraft.wing_area_m2 * aircraft.cl_max_takeoff)
    )
    vr = options.vr_over_vs * vs
    vlof = options.vlof_over_vs * vs
    v50 = options.v50_over_vs * vs
    if not (vr > vs and vlof > vr and v50 >= vlof):
        raise ValueError("Se requiere V_R > V_S, V_LOF > V_R y V_50 >= V_LOF.")
    if airport.headwind_mps >= vr:
        raise ValueError("El viento de frente no puede alcanzar V_R en este modelo.")

    total_power_w, hp_per_engine = shaft_power_available(
        atmosphere.pressure_altitude_m,
        aircraft.engine_map_inhg,
        aircraft.engine_count,
    )
    slope_rad = atan(airport.runway_slope_percent / 100.0)

    def ground_state(v_air: float) -> dict[str, float]:
        q = 0.5 * rho * v_air**2
        lift = q * aircraft.wing_area_m2 * aircraft.cl_ground
        drag = q * aircraft.wing_area_m2 * aircraft.cd_ground
        eta, thrust = propeller_efficiency_and_thrust(
            v_air,
            rho,
            total_power_w,
            aircraft.engine_count,
            aircraft.propeller_diameter_m,
            aircraft.propeller_profile_efficiency,
        )
        # Clase Performance Despegue/Aterrizaje, Ec. (4): pendiente positiva
        # cuesta arriba y, por lo tanto, desfavorable.
        net_force = (
            thrust
            - drag
            - airport.rolling_friction * (weight_n - lift)
            - weight_n * slope_rad
        )
        acceleration = G * net_force / weight_n
        return {
            "V_air_mps": v_air,
            "V_ground_mps": v_air - airport.headwind_mps,
            "lift_n": lift,
            "drag_n": drag,
            "eta_prop": eta,
            "thrust_n": thrust,
            "net_force_n": net_force,
            "acceleration_mps2": acceleration,
        }

    # Clase, Ec. (8): V es velocidad respecto del aire y V_suelo=V-V_frente.
    v_start = airport.headwind_mps
    profile: list[dict[str, float]] = []
    if options.ground_method == "general":
        dv = (vr - v_start) / options.intervals
        s_ngr = 0.0
        for index in range(options.intervals):
            v_mid = v_start + (index + 0.5) * dv
            state = ground_state(v_mid)
            acceleration = state["acceleration_mps2"]
            if acceleration <= 0.0 or not isfinite(acceleration):
                raise ValueError(
                    f"Aceleracion no positiva durante s_NGR en V={v_mid:.3f} m/s."
                )
            ds = state["V_ground_mps"] * dv / acceleration
            if ds < 0.0:
                raise ValueError("La integracion produjo una distancia negativa.")
            s_ngr += ds
            state["cumulative_distance_m"] = s_ngr
            profile.append(state)
    else:
        # Clase, Ecs. (9)-(10): para helice evaluar a_g en 0.74 V_R;
        # esta aproximacion fue presentada expresamente para viento nulo.
        if abs(airport.headwind_mps) > 1.0e-12:
            raise ValueError("El metodo aproximado de clase requiere viento nulo.")
        state = ground_state(0.74 * vr)
        if state["acceleration_mps2"] <= 0.0:
            raise ValueError("Aceleracion no positiva en 0.74 V_R.")
        s_ngr = vr**2 / (2.0 * state["acceleration_mps2"])
        state["cumulative_distance_m"] = s_ngr
        profile.append(state)

    # Clase, Ec. (11), conservando literalmente su termino de viento dentro
    # del factor 1/2: s_R=0.5[(V_R+V_LOF)-V_frente]t_R.
    s_rotation = 0.5 * (
        vr + vlof - airport.headwind_mps
    ) * options.rotation_time_s
    if s_rotation <= 0.0:
        raise ValueError("La distancia de rotacion no es positiva.")

    # Clase, ejemplo 2: velocidad media entre lift-off y obstaculo.
    v_transition = 0.5 * (vlof + v50)
    ratio = v_transition / vs
    delta_cl = 0.5 * (ratio**2 - 1.0) * (
        aircraft.cl_max_takeoff * ((1.0 / ratio) ** 2 - 0.53) + 0.38
    )
    if delta_cl <= 0.0:
        raise ValueError("Delta CL de transicion no positivo.")

    # Clase, Ec. (12). Se conserva V_LOF en esta magnitud de chequeo.
    cl_transition = (
        2.0 * weight_n
        / (rho * vlof**2 * aircraft.wing_area_m2)
        + delta_cl
    )
    # Para T-D en V_TR, el ejemplo 2 usa CL de equilibrio a V_TR.
    cl_at_v_transition = (
        2.0 * weight_n
        / (rho * v_transition**2 * aircraft.wing_area_m2)
    )
    cd_transition = drag_coefficient(cl_at_v_transition, aircraft)
    drag_transition = (
        0.5
        * rho
        * v_transition**2
        * aircraft.wing_area_m2
        * cd_transition
    )
    eta_transition, thrust_transition = propeller_efficiency_and_thrust(
        v_transition,
        rho,
        total_power_w,
        aircraft.engine_count,
        aircraft.propeller_diameter_m,
        aircraft.propeller_profile_efficiency,
    )
    sin_gamma = (thrust_transition - drag_transition) / weight_n
    if not 0.0 < sin_gamma < 1.0:
        raise ValueError("No existe un angulo de trepada positivo en V_TR.")
    gamma_c = asin(sin_gamma)

    # Clase, Ecs. (15)-(23): transicion circular y trepada rectilinea.
    radius_transition = (
        2.0 * wing_loading / (rho * G * delta_cl)
    )
    h_transition_full = radius_transition * (1.0 - cos(gamma_c))
    if h_transition_full < options.obstacle_height_m:
        s_transition = radius_transition * sin(gamma_c)
        h_transition = h_transition_full
        s_climb = (options.obstacle_height_m - h_transition) / tan(gamma_c)
        obstacle_phase = "trepada"
    else:
        gamma_at_obstacle = acos(
            1.0 - options.obstacle_height_m / radius_transition
        )
        s_transition = radius_transition * sin(gamma_at_obstacle)
        h_transition = options.obstacle_height_m
        s_climb = 0.0
        obstacle_phase = "transicion"

    s_ground = s_ngr + s_rotation
    s_air = s_transition + s_climb
    s_takeoff = s_ground + s_air
    far23_factor = 1.0
    s_far23 = far23_factor * s_takeoff

    v_mc = aircraft.v_mc_mps
    if v_mc is None:
        far23_speed_compliant: bool | None = None
    else:
        far23_speed_compliant = vr >= v_mc and v50 >= 1.10 * v_mc

    state_static = ground_state(0.0)
    state_074vr = ground_state(0.74 * vr)
    state_vr = ground_state(vr)
    accelerations = [state["acceleration_mps2"] for state in profile]
    return {
        "atmosphere": atmosphere,
        "weight_n": weight_n,
        "wing_loading_npm2": wing_loading,
        "power_total_w": total_power_w,
        "power_per_engine_hp": hp_per_engine,
        "V_S_mps": vs,
        "V_R_mps": vr,
        "V_LOF_mps": vlof,
        "V_50_mps": v50,
        "V_TR_mps": v_transition,
        "delta_CL": delta_cl,
        "CL_TR": cl_transition,
        "CL_at_V_TR": cl_at_v_transition,
        "CD_TR": cd_transition,
        "eta_prop_TR": eta_transition,
        "thrust_TR_n": thrust_transition,
        "drag_TR_n": drag_transition,
        "gamma_c_rad": gamma_c,
        "R_TR_m": radius_transition,
        "h_TR_m": h_transition,
        "obstacle_phase": obstacle_phase,
        "s_NGR_m": s_ngr,
        "s_R_m": s_rotation,
        "s_TR_m": s_transition,
        "s_CL_m": s_climb,
        "s_G_m": s_ground,
        "s_A_m": s_air,
        "s_TO_m": s_takeoff,
        "far23_factor": far23_factor,
        "s_FAR23_m": s_far23,
        "far23_speed_compliant": far23_speed_compliant,
        "runway_margin_m": airport.runway_available_m - s_far23,
        "ground_profile": profile,
        "state_static": state_static,
        "state_074VR": state_074vr,
        "state_VR": state_vr,
        "min_ground_acceleration_mps2": min(accelerations),
        "max_ground_acceleration_mps2": max(accelerations),
    }
