#!/usr/bin/env python3
"""Calculo editable de distancia de despegue del Beechcraft Duke 60."""

from math import degrees

from takeoff_model import (
    Aircraft,
    AirportConditions,
    FT_TO_M,
    ModelOptions,
    calculate_takeoff,
)


# ==========================================================
# 1. CONDICIONES DEL AEROPUERTO
# ==========================================================

# CASO DE ESTUDIO: nivel del mar, ISA, pista seca de concreto, sin viento.
altitud_aeropuerto_m = 489.0
temperatura_c = 20.0
presion_aeropuerto_pa = None  # None: presion ISA a la altitud indicada
viento_frente_mps = 0.0       # positivo: frente; negativo: cola
pendiente_pista_porcentaje = 0.0  # positivo: subida en sentido de despegue
tipo_superficie = "concreto seco"
mu_rodadura = 0.025           # clase/Roskam: concreto 0.02-0.03
longitud_pista_disponible_m = 3000


# ==========================================================
# 2. CONDICION DEL AVION
# ==========================================================

# Correspondencia futura: Flight-Mechanics-Codes/Duke.py
masa_despegue_kg = 3073.0  # M [kg], condicion del proyecto existente
S_m2 = 19.780               # S [m2]
CL_max_TO = 1.27            # CL_max preliminar del avion
CD0 = 0.026
K = 0.0596

CL_g = 0.40
CD_g = CD0 + K * CL_g**2

MAP_despegue_inhg = 40.0
numero_motores = 2
diametro_helice_m = 1.8796
eta_perfil_helice = 0.88

# V_MC certificada del Duke en esta configuracion [m/s].
# Sin ella se calcula la distancia fisica al obstaculo, pero el chequeo completo
# de velocidades FAR 23 queda marcado como pendiente.
V_MC_mps = None


# ==========================================================
# 3. OPCIONES DEL MODELO
# ==========================================================

metodo_rodaje = "general"  # "general" (punto medio) o "aproximado"
numero_intervalos = 100
tiempo_rotacion_s = 1.0     # clase/Roskam: avion liviano
factor_VR_VS = 1.10         # Roskam 10.3.3.1, estimacion preliminar
factor_VLOF_VS = 1.15       # Roskam 10.3.3.1, estimacion preliminar
factor_V50_VS = 1.20        # FAR 23 multiengine, Roskam tabla 10.1
altura_obstaculo_m = 50.0 * FT_TO_M  # clase: 50 ft para FAR 23


def _kt(value_mps: float) -> float:
    return value_mps / 0.514444


def main() -> None:
    airport = AirportConditions(
        altitude_m=altitud_aeropuerto_m,
        temperature_c=temperatura_c,
        pressure_pa=presion_aeropuerto_pa,
        headwind_mps=viento_frente_mps,
        runway_slope_percent=pendiente_pista_porcentaje,
        rolling_friction=mu_rodadura,
        runway_surface=tipo_superficie,
        runway_available_m=longitud_pista_disponible_m,
    )
    aircraft = Aircraft(
        mass_kg=masa_despegue_kg,
        wing_area_m2=S_m2,
        cl_max_takeoff=CL_max_TO,
        cd0=CD0,
        induced_drag_k=K,
        cl_ground=CL_g,
        cd_ground=CD_g,
        engine_map_inhg=MAP_despegue_inhg,
        engine_count=numero_motores,
        propeller_diameter_m=diametro_helice_m,
        propeller_profile_efficiency=eta_perfil_helice,
        v_mc_mps=V_MC_mps,
    )
    options = ModelOptions(
        ground_method=metodo_rodaje,
        intervals=numero_intervalos,
        rotation_time_s=tiempo_rotacion_s,
        vr_over_vs=factor_VR_VS,
        vlof_over_vs=factor_VLOF_VS,
        v50_over_vs=factor_V50_VS,
        obstacle_height_m=altura_obstaculo_m,
    )
    result = calculate_takeoff(airport, aircraft, options)
    atmosphere = result["atmosphere"]

    print("CASO DE VERIFICACION")
    print("====================")
    print("\nCONDICIONES")
    print("-----------")
    print(f"Altitud aeropuerto:       {altitud_aeropuerto_m:.1f} m")
    print(f"Temperatura:              {temperatura_c:.2f} degC")
    print(f"Presion:                  {atmosphere.pressure_pa:.2f} Pa")
    print(f"Altitud de presion:       {atmosphere.pressure_altitude_m:.2f} m")
    print(f"Densidad rho:             {atmosphere.density_kgm3:.4f} kg/m3")
    print(f"Velocidad del sonido:     {atmosphere.speed_of_sound_mps:.2f} m/s")
    print(f"Viento de frente:         {viento_frente_mps:.2f} m/s")
    print(f"Pendiente (subida +):     {pendiente_pista_porcentaje:.3f} %")
    print(f"Superficie:               {tipo_superficie}")
    print(f"mu_g:                     {mu_rodadura:.4f}")
    print(f"Masa de despegue:         {masa_despegue_kg:.2f} kg")
    print(f"Peso de despegue W:       {result['weight_n']:.2f} N")

    print("\nVELOCIDADES")
    print("-----------")
    for label, key in (
        ("V_S", "V_S_mps"),
        ("V_R", "V_R_mps"),
        ("V_LOF", "V_LOF_mps"),
        ("V_50", "V_50_mps"),
        ("V_TR", "V_TR_mps"),
    ):
        value = result[key]
        print(f"{label:8s}:                 {value:7.3f} m/s ({_kt(value):7.3f} kt)")

    print("\nRESULTADOS INTERMEDIOS")
    print("----------------------")
    print(f"CL_g / CD_g:              {CL_g:.4f} / {CD_g:.5f}")
    print(f"Potencia por motor:       {result['power_per_engine_hp']:.3f} hp")
    print(f"Potencia total al eje:    {result['power_total_w'] / 1000:.3f} kW")
    for label, key in (
        ("V=0", "state_static"),
        ("V=0.74 V_R", "state_074VR"),
        ("V=V_R", "state_VR"),
    ):
        state = result[key]
        print(
            f"{label:10s}: T={state['thrust_n']:8.2f} N, "
            f"L={state['lift_n']:8.2f} N, D={state['drag_n']:7.2f} N, "
            f"a_g={state['acceleration_mps2']:.4f} m/s2"
        )
    print(
        "Rango a_g integracion:    "
        f"{result['min_ground_acceleration_mps2']:.4f} a "
        f"{result['max_ground_acceleration_mps2']:.4f} m/s2"
    )
    print(f"Delta CL transicion:      {result['delta_CL']:.5f}")
    print(f"CL_TR (Ec. 12):           {result['CL_TR']:.5f}")
    print(f"CL / CD en V_TR:          {result['CL_at_V_TR']:.5f} / {result['CD_TR']:.5f}")
    print(f"eta helice en V_TR:       {result['eta_prop_TR']:.5f}")
    print(
        f"T / D en V_TR:           {result['thrust_TR_n']:.2f} / "
        f"{result['drag_TR_n']:.2f} N"
    )
    print(
        f"Angulo trepada gamma_c:   {degrees(result['gamma_c_rad']):.3f} deg "
        f"({result['gamma_c_rad']:.5f} rad)"
    )
    print(f"Radio transicion R_TR:    {result['R_TR_m']:.3f} m")
    print(f"Altura al fin de s_TR:    {result['h_TR_m']:.3f} m")
    print(f"Obstaculo alcanzado en:   {result['obstacle_phase']}")

    print("\nDISTANCIAS")
    print("----------")
    print(f"Nariz en pista s_NGR:     {result['s_NGR_m']:.3f} m")
    print(f"Rotacion s_R:             {result['s_R_m']:.3f} m")
    print(f"Transicion s_TR:          {result['s_TR_m']:.3f} m")
    print(f"Trepada s_CL:             {result['s_CL_m']:.3f} m")
    print(f"Distancia en suelo s_G:   {result['s_G_m']:.3f} m")
    print(f"Distancia en aire s_A:    {result['s_A_m']:.3f} m")

    print("\nDISTANCIA TOTAL DE DESPEGUE")
    print("---------------------------")
    print(f"s_TO =                    {result['s_TO_m']:.3f} m")

    print("\nDISTANCIA REQUERIDA FAR 23")
    print("--------------------------")
    print(f"Altura de obstaculo:      {altura_obstaculo_m:.3f} m (50 ft)")
    print(f"Factor reglamentario:     {result['far23_factor']:.3f}")
    print(f"s_FAR23 =                 {result['s_FAR23_m']:.3f} m")
    if result["far23_speed_compliant"] is None:
        print("Chequeo V_MC:              PENDIENTE (TODO_USER_INPUT)")
    else:
        status = "CUMPLE" if result["far23_speed_compliant"] else "NO CUMPLE"
        print(f"Chequeo V_MC:              {status}")
    print(f"Pista disponible:         {longitud_pista_disponible_m:.3f} m")
    print(f"Margen de pista:          {result['runway_margin_m']:.3f} m")


if __name__ == "__main__":
    main()
