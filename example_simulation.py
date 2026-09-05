"""
Example: Running an FD-PBD simulation on a PMMA/Al/SiO2 sample.

This script demonstrates how to:
1. Configure materials and multilayer sample stack.
2. Set up optical beam and detector settings.
3. Simulate frequency-domain probe beam deflection signals.
4. Decompose signals into heat and mass transport contributions.
5. Plot amplitude, phase, in-phase (Vin), and out-of-phase (Vout) response.
"""

import numpy as np
import matplotlib.pyplot as plt

from fdpbd import (
    Material,
    PMMA,
    ALUMINUM,
    FUSED_SILICA,
    WATER_VAPOR,
    Sample,
    LaserBeamConfig,
    DetectorConfig,
    MassTransportConfig,
    SimulationGrid,
    FDPBDSimulator,
    plot_all,
    plot_displacement_breakdown,
)


def main():
    print("==================================================")
    print("FD-PBD Heat & Mass Transport Simulation")
    print("==================================================")

    # 1. Define Sample Stack
    # Standard 4-layer structure: Water Vapor / PMMA (114 nm) / Al (80 nm) / Fused Silica
    sample = Sample.create_standard_fdpbd_sample(
        polymer_thickness=114.0e-9,      # 114 nm PMMA
        al_thickness=80.0e-9,           # 80 nm Al transducer
        polymer_material=PMMA,
        ambient_material=WATER_VAPOR,
        transducer_material=ALUMINUM,
        substrate_material=FUSED_SILICA,
    )

    # 2. Configure Laser Optics & Detector
    # 5X objective lens configuration (W0 = 11.3 um, r0 = 12.0 um, F = 40 mm)
    beam = LaserBeamConfig(
        pump_power=3.0e-3,              # 3.0 mW pump power
        wavelength=783.0e-9,            # 783 nm Ti:sapphire laser
        pump_radius=11.3e-6,            # 11.3 um 1/e^2 beam radius
        beam_offset=12.0e-6,            # 12.0 um pump-probe offset
        focal_length=40.0e-3,           # 40 mm focal length
        trans_objective=0.90,           # Objective transmission
        trans_window=0.92,              # Chamber window transmission
    )

    detector = DetectorConfig(
        v0_multimeter=830.0e-3,         # 830 mV multimeter DC sum voltage
        gain2=1.0,                      # PSD second-stage gain
        free_space_radius=1.55e-3 / 2.0,# Free-space beam radius at PSD
    )

    # 3. Configure Mass Transport Properties
    mass = MassTransportConfig(
        Dm=1.1e-12,                     # Water diffusivity in PMMA (m^2/s)
        dndT_mass=-1.0e-4,              # Mass thermo-optic equivalent (1/K)
        aCET_mass=-2.0e-5,              # Mass expansion equivalent (1/K)
    )

    # 4. Frequency Grid (10 Hz to 100 kHz)
    frequencies = np.logspace(np.log10(1.0), np.log10(100.0e3), 100)

    # 5. Run Simulation
    print("Running simulation...")
    result = FDPBDSimulator.simulate(
        sample=sample,
        frequencies=frequencies,
        beam_config=beam,
        detector_config=detector,
        mass_config=mass,
        decompose_components=True,
    )

    print(f"-> Sample Reflectance: {result.reflectance:.4f}")
    print(f"-> Absorbed Pump Amplitude: {result.absorbed_power_A * 1e3:.2f} mW")
    print(f"-> Steady-State Temperature Rise (dT_SS): {result.steady_state_dT:.2f} K")
    print(f"-> Max In-Phase Signal (Vin): {np.max(result.Vin):.2f} uV")
    print(f"-> Min Out-of-Phase Signal (Vout): {np.min(result.Vout):.2f} uV")

    # 6. Plot Comprehensive 4-Panel Results
    fig, axes = plot_all(
        sim_result=result,
        title="FD-PBD Simulation: PMMA (114 nm) / Al (80 nm) / Fused Silica in Water Vapor",
    )
    plt.savefig("fdpbd_simulation_results.png", dpi=200)
    print("Saved simulation plot: fdpbd_simulation_results.png")

    # 7. Plot Displacement Breakdown (Z0 to Z4)
    fig_disp, axes_disp = plot_displacement_breakdown(result, k_index=0)
    plt.savefig("fdpbd_displacement_breakdown.png", dpi=200)
    print("Saved displacement breakdown plot: fdpbd_displacement_breakdown.png")
    plt.close("all")


if __name__ == "__main__":
    main()
