"""
Detailed benchmark and profiler for FD-PBD.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import time
import numpy as np
from fdpbd import (
    Sample,
    PMMA,
    ALUMINUM,
    FUSED_SILICA,
    WATER_VAPOR,
    LaserBeamConfig,
    DetectorConfig,
    MassTransportConfig,
    SimulationGrid,
    solve_thermal_field,
    solve_fresnel_reflection,
    solve_multilayer_optics_perturbed,
    solve_displacements,
    calculate_beam_deflection,
    FDPBDSimulator,
)

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
grid = SimulationGrid()
k = grid.get_k_grid(beam.pump_radius)

print(f"Grid size: Nk = {len(k)}, Nf = {len(frequencies)}")

# 1. Profile Fresnel reflection
t0 = time.perf_counter()
R, T, phi = solve_fresnel_reflection(sample)
t1 = time.perf_counter()
print(f"Fresnel reflection: {(t1-t0)*1000:.3f} ms, R = {R:.4f}")

# 2. Profile Thermal field
t0 = time.perf_counter()
thermal_res = solve_thermal_field(sample, frequencies, k, heating_layer_idx=sample.heating_layer_idx)
t1 = time.perf_counter()
print(f"Thermal field (heating): {(t1-t0)*1000:.3f} ms")

# 3. Profile Perturbed Optics
t0 = time.perf_counter()
optics_res = solve_multilayer_optics_perturbed(
    sample, frequencies, k, thermal_res.G_tot, mass, beam.wavelength
)
t1 = time.perf_counter()
print(f"Perturbed optics: {(t1-t0)*1000:.3f} ms")

# 4. Profile Displacements
t0 = time.perf_counter()
disp = solve_displacements(sample, frequencies, k, mass, beam.wavelength, optics_res)
t1 = time.perf_counter()
print(f"Displacements (Z0..Z4): {(t1-t0)*1000:.3f} ms")

# 5. Profile Beam Deflection
t0 = time.perf_counter()
theta, sig = calculate_beam_deflection(disp.Z_tot_m, k, beam, detector, R)
t1 = time.perf_counter()
print(f"Deflection quadrature: {(t1-t0)*1000:.3f} ms")

print("\nSignal values at f = [1, 10, 100, 1k, 10k, 100k] Hz:")
test_freqs = [1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0]
for tf in test_freqs:
    idx = np.argmin(np.abs(frequencies - tf))
    print(f"f = {frequencies[idx]:8.1f} Hz: Vin = {np.real(sig[idx]):8.2f} uV, Vout = {np.imag(sig[idx]):8.2f} uV, Phase = {np.angle(sig[idx], deg=True):6.1f} deg")

print("\nBreakdown of components at f = 100 Hz (k_idx=0):")
f_100_idx = np.argmin(np.abs(frequencies - 100.0))
print(f"Z0_m: {disp.Z0_m[0, f_100_idx]*1e12:.3e} pm")
print(f"Z1_m: {disp.Z1_m[0, f_100_idx]*1e12:.3e} pm")
print(f"Z2_m: {disp.Z2_m[0, f_100_idx]*1e12:.3e} pm")
print(f"Z3_m: {disp.Z3_m[0, f_100_idx]*1e12:.3e} pm")
print(f"Z4_m: {disp.Z4_m[0, f_100_idx]*1e12:.3e} pm")
print(f"Z_tot_m: {disp.Z_tot_m[0, f_100_idx]*1e12:.3e} pm")
