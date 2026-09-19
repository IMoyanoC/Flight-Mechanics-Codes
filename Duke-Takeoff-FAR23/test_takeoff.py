"""Comprobaciones numericas del CASO DE VERIFICACION, no datos operativos."""

from dataclasses import replace
import unittest

from takeoff_model import Aircraft, AirportConditions, ModelOptions, calculate_takeoff


AIRPORT_VERIFICATION = AirportConditions(
    altitude_m=0.0,
    temperature_c=15.0,
    pressure_pa=None,
    headwind_mps=0.0,
    runway_slope_percent=0.0,
    rolling_friction=0.025,
    runway_surface="concreto seco",
    runway_available_m=1_000.0,
)
AIRCRAFT_VERIFICATION = Aircraft(
    mass_kg=2_795.0,
    wing_area_m2=19.780,
    cl_max_takeoff=1.27,
    cd0=0.026,
    induced_drag_k=0.0596,
    cl_ground=0.40,
    cd_ground=0.026 + 0.0596 * 0.40**2,
    engine_map_inhg=40.0,
    engine_count=2,
    propeller_diameter_m=1.8796,
    propeller_profile_efficiency=0.88,
)


class TakeoffVerificationTests(unittest.TestCase):
    def test_phase_sums(self) -> None:
        result = calculate_takeoff(
            AIRPORT_VERIFICATION, AIRCRAFT_VERIFICATION, ModelOptions(intervals=100)
        )
        self.assertAlmostEqual(
            result["s_TO_m"],
            result["s_NGR_m"]
            + result["s_R_m"]
            + result["s_TR_m"]
            + result["s_CL_m"],
            places=10,
        )

    def test_midpoint_convergence(self) -> None:
        result_100 = calculate_takeoff(
            AIRPORT_VERIFICATION, AIRCRAFT_VERIFICATION, ModelOptions(intervals=100)
        )
        result_1000 = calculate_takeoff(
            AIRPORT_VERIFICATION, AIRCRAFT_VERIFICATION, ModelOptions(intervals=1000)
        )
        self.assertLess(abs(result_100["s_TO_m"] - result_1000["s_TO_m"]), 0.01)

    def test_wind_sign(self) -> None:
        calm = calculate_takeoff(
            AIRPORT_VERIFICATION, AIRCRAFT_VERIFICATION, ModelOptions()
        )["s_TO_m"]
        headwind = calculate_takeoff(
            replace(AIRPORT_VERIFICATION, headwind_mps=5.0),
            AIRCRAFT_VERIFICATION,
            ModelOptions(),
        )["s_TO_m"]
        tailwind = calculate_takeoff(
            replace(AIRPORT_VERIFICATION, headwind_mps=-5.0),
            AIRCRAFT_VERIFICATION,
            ModelOptions(),
        )["s_TO_m"]
        self.assertLess(headwind, calm)
        self.assertGreater(tailwind, calm)

    def test_slope_sign(self) -> None:
        downhill = calculate_takeoff(
            replace(AIRPORT_VERIFICATION, runway_slope_percent=-2.0),
            AIRCRAFT_VERIFICATION,
            ModelOptions(),
        )["s_TO_m"]
        uphill = calculate_takeoff(
            replace(AIRPORT_VERIFICATION, runway_slope_percent=2.0),
            AIRCRAFT_VERIFICATION,
            ModelOptions(),
        )["s_TO_m"]
        self.assertLess(downhill, uphill)


if __name__ == "__main__":
    unittest.main()
