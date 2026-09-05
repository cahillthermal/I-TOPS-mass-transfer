"""
Test suite for FD-PBD fitting routines.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import time
import numpy as np
from fdpbd import (
    Sample,
    PMMA,
    LaserBeamConfig,
    DetectorConfig,
    MassTransportConfig,
    FDPBDSimulator,
    FDPBDFitter,
    ExperimentData,
)


def test_fitter_recovery():
    # 1. Generate synthetic experimental data with known Dm
    true_Dm = 1.2e-12
    sample = Sample.create_standard_fdpbd_sample(
        polymer_thickness=114e-9,
        al_thickness=80e-9,
    )
    frequencies = np.logspace(np.log10(10), np.log10(10e3), 30)
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
    true_mass = MassTransportConfig(
        Dm=true_Dm,
        dndT_mass=-1.0e-4,
        aCET_mass=-2.0e-5,
    )

    sim_synth = FDPBDSimulator.simulate(
        sample=sample,
        frequencies=frequencies,
        beam_config=beam,
        detector_config=detector,
        mass_config=true_mass,
        decompose_components=False,
    )

    synth_exp = ExperimentData(
        frequencies=frequencies,
        Vin=sim_synth.Vin,
        Vout=sim_synth.Vout,
    )

    # 2. Setup Fitter with offset initial guess
    fitter = FDPBDFitter(
        sample=sample,
        exp_data=synth_exp,
        beam_config=beam,
        detector_config=detector,
        mass_config=MassTransportConfig(Dm=0.5e-12),
    )
    fitter.add_parameter("Dm", value=0.5e-12, bounds=(1e-14, 1e-10))

    t0 = time.perf_counter()
    fit_res = fitter.fit(method="least_squares", loss="phase", max_nfev=50)
    t_elapsed = time.perf_counter() - t0

    fitted_Dm = fit_res.params["Dm"]
    print(f"Fitting elapsed time: {t_elapsed:.2f} s")
    print(f"True Dm: {true_Dm:.3e}, Fitted Dm: {fitted_Dm:.3e}")
    print(f"Optimizer success: {fit_res.success}, cost: {fit_res.cost:.4e}")

    assert abs(fitted_Dm - true_Dm) / true_Dm < 0.05
    print("[PASS] Parameter recovery test passed!")


if __name__ == "__main__":
    test_fitter_recovery()
