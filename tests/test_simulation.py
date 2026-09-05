"""
Test suite to verify Python FD-PBD simulation physics, speed, and numerical accuracy.
"""

import sys
import os
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import time
import numpy as np

from fdpbd import (
    Material,
    AIR,
    WATER_VAPOR,
    PMMA,
    ALUMINUM,
    FUSED_SILICA,
    Sample,
    Layer,
    LaserBeamConfig,
    DetectorConfig,
    MassTransportConfig,
    SimulationGrid,
    FDPBDSimulator,
    solve_thermal_field,
    solve_fresnel_reflection,
    calculate_steady_state_heating,
)


def test_thermal_field():
    sample = Sample.create_standard_fdpbd_sample(
        polymer_thickness=100e-9,
        al_thickness=80e-9,
    )
    frequencies = np.array([10.0, 100.0, 1000.0, 10000.0])
    W0 = 11.3e-6
    k = np.arange(0.01 / W0, 5.0 / W0 + 0.005 / W0, 0.01 / W0)

    t_start = time.perf_counter()
    res = solve_thermal_field(
        sample=sample,
        frequencies=frequencies,
        k_vectors=k,
        heating_layer_idx=sample.heating_layer_idx,
    )
    t_elapsed = time.perf_counter() - t_start

    assert res.G_tot.shape == (len(k), len(frequencies))
    assert not np.any(np.isnan(res.G_tot))
    assert not np.any(np.isinf(res.G_tot))
    assert t_elapsed < 0.1
    print(f"[PASS] Thermal field evaluation time: {t_elapsed * 1000:.2f} ms")


def test_fresnel_reflection():
    sample = Sample.create_standard_fdpbd_sample(
        polymer_thickness=100e-9,
        al_thickness=80e-9,
    )
    R, T, phi = solve_fresnel_reflection(sample)
    print(f"[PASS] Reflectance R = {R:.4f}, Phase phi = {phi:.4f} rad")
    assert 0.70 < R < 0.95
    assert -np.pi <= phi <= np.pi


def test_steady_state_heating():
    sample = Sample.create_standard_fdpbd_sample(
        polymer_thickness=100e-9,
        al_thickness=80e-9,
    )
    W0 = 11.3e-6
    k = np.arange(0.01 / W0, 5.0 / W0 + 0.005 / W0, 0.01 / W0)
    R, _, _ = solve_fresnel_reflection(sample)
    beam = LaserBeamConfig(pump_power=1.5e-3, pump_radius=W0)
    A_pump = beam.absorbed_power_amplitude(R)
    A_SS = np.pi * A_pump

    dT_ss = calculate_steady_state_heating(
        sample=sample,
        k_vectors=k,
        pump_power_absorbed=A_SS,
        pump_radius=W0,
    )
    print(f"[PASS] Steady state dT_SS (at P=1.5mW) = {dT_ss:.2f} K")
    assert 5.0 < dT_ss < 20.0


def test_full_fdpbd_simulation():
    sample = Sample.create_standard_fdpbd_sample(
        polymer_thickness=114e-9,
        al_thickness=80e-9,
    )
    frequencies = np.logspace(np.log10(1), np.log10(100e3), 100)
    beam = LaserBeamConfig(
        pump_power=3.0e-3,
        wavelength=783e-9,
        pump_radius=11.3e-6,
        beam_offset=12.0e-6,
        focal_length=40.0e-3,
    )
    detector = DetectorConfig(
        v0_multimeter=830e-3,
        gain2=1.0,
        free_space_radius=1.55e-3 / 2.0,
    )
    mass = MassTransportConfig(
        Dm=1.1e-12,
        dndT_mass=-1.0e-4,
        aCET_mass=-2.0e-5,
    )

    t_start = time.perf_counter()
    sim_result = FDPBDSimulator.simulate(
        sample=sample,
        frequencies=frequencies,
        beam_config=beam,
        detector_config=detector,
        mass_config=mass,
        decompose_components=True,
    )
    t_elapsed = time.perf_counter() - t_start

    print(f"[PASS] Full 100-frequency simulation time: {t_elapsed * 1000:.2f} ms")
    assert len(sim_result.signal_total) == 100
    assert not np.any(np.isnan(sim_result.signal_total))
    assert sim_result.signal_heat is not None
    assert sim_result.signal_mass is not None

    assert np.all(sim_result.Vin > 0)
    print(f"[PASS] Max Vin: {np.max(sim_result.Vin):.2f} uV, Min Vin: {np.min(sim_result.Vin):.2f} uV")
    print(f"[PASS] Min Vout: {np.min(sim_result.Vout):.2f} uV, Max Vout: {np.max(sim_result.Vout):.2f} uV")


if __name__ == "__main__":
    test_thermal_field()
    test_fresnel_reflection()
    test_steady_state_heating()
    test_full_fdpbd_simulation()
    print("\nAll FD-PBD core simulation tests passed successfully!")
