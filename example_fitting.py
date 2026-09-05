"""
Example: Loading experimental lock-in data, calibrating, and fitting parameters with FD-PBD.

This script demonstrates how to:
1. Load raw lock-in data and detector calibration files.
2. Apply frequency-dependent phase and amplitude calibration.
3. Configure sample parameters and select free parameters (e.g. Dm, dndT_mass, aCET_mass).
4. Run least-squares optimization.
5. Print fitted values, standard errors, and plot the fit against experiment.
"""

import os
import numpy as np
import matplotlib.pyplot as plt

from fdpbd import (
    Sample,
    PMMA,
    ALUMINUM,
    FUSED_SILICA,
    WATER_VAPOR,
    LaserBeamConfig,
    DetectorConfig,
    MassTransportConfig,
    ExperimentData,
    CalibrationData,
    FDPBDFitter,
    plot_fit,
)


def main():
    print("==================================================")
    print("FD-PBD Experimental Data Calibration & Fitting")
    print("==================================================")

    # File names from MATLAB project
    exp_file = "091917_2_R23_Scan_Vacuum(0.04psi)_PMMA(114nm)_Al_SiO2_5X_pump3mW_probe6mW_TC3s_12dB_Sens200uV_PSD920mV_Wait15s_dewell15s_12um.txt"
    cali_file = "091917_6_Cali_R23_Before.txt"

    # Check if files exist in current directory
    if os.path.exists(exp_file) and os.path.exists(cali_file):
        print(f"Loading experiment data: {exp_file}")
        raw_exp = ExperimentData.load_from_file(exp_file, fix_legacy_labview=True)
        print(f"Loading calibration data: {cali_file}")
        cali = CalibrationData.load_from_file(cali_file)

        print("Applying phase and amplitude calibration...")
        calibrated_exp = raw_exp.apply_calibration(cali, correct_phase=True, correct_amplitude=True)
        # Filter frequency range to 10 Hz - 20 kHz
        exp_data = calibrated_exp.filter_frequency_range(f_min=10.0, f_max=20000.0)
    else:
        print("Experiment files not found, generating synthetic experimental data with noise...")
        sample_true = Sample.create_standard_fdpbd_sample(
            polymer_thickness=114e-9,
            al_thickness=78e-9,
        )
        f_exp = np.logspace(np.log10(10), np.log10(20e3), 40)
        beam_synth = LaserBeamConfig(pump_power=3.0e-3, pump_radius=11.3e-6, beam_offset=12.0e-6)
        detector_synth = DetectorConfig(v0_multimeter=830e-3)
        mass_synth = MassTransportConfig(Dm=1.2e-12, dndT_mass=-1.1e-4, aCET_mass=-2.2e-5)

        from fdpbd import FDPBDSimulator
        sim_true = FDPBDSimulator.simulate(
            sample=sample_true,
            frequencies=f_exp,
            beam_config=beam_synth,
            detector_config=detector_synth,
            mass_config=mass_synth,
            decompose_components=False,
        )
        # Add slight 1% Gaussian noise
        noise_in = np.random.normal(0, 0.01 * np.max(sim_true.Vin), len(f_exp))
        noise_out = np.random.normal(0, 0.01 * np.abs(np.min(sim_true.Vout)), len(f_exp))
        exp_data = ExperimentData(
            frequencies=f_exp,
            Vin=sim_true.Vin + noise_in,
            Vout=sim_true.Vout + noise_out,
        )

    # 1. Define sample for fitting
    sample = Sample.create_standard_fdpbd_sample(
        polymer_thickness=114.0e-9,
        al_thickness=78.0e-9,
    )

    beam = LaserBeamConfig(
        pump_power=3.0e-3,
        wavelength=783.0e-9,
        pump_radius=11.3e-6,
        beam_offset=12.0e-6,
        focal_length=40.0e-3,
    )
    detector = DetectorConfig(
        v0_multimeter=830.0e-3,
        gain2=1.0,
        free_space_radius=1.55e-3 / 2.0,
    )

    # 2. Setup Fitter and parameters
    fitter = FDPBDFitter(
        sample=sample,
        exp_data=exp_data,
        beam_config=beam,
        detector_config=detector,
    )

    # Configure free parameters to fit with initial guesses and physical bounds
    fitter.add_parameter("Dm", value=0.8e-12, bounds=(1.0e-14, 1.0e-10))
    fitter.add_parameter("dndT_mass", value=-0.8e-4, bounds=(-5.0e-4, 0.0))

    print("\nInitial Parameter Guesses:")
    for name, p in fitter.params.items():
        print(f"  {name}: {p.value:.3e} (bounds: {p.bounds})")

    # 3. Run Non-linear Least Squares Fit
    print("\nExecuting least-squares optimization...")
    fit_result = fitter.fit(method="least_squares", loss="phase", max_nfev=60)

    print("\n================ Fit Results ================")
    print(f"Convergence Success: {fit_result.success}")
    print(f"Optimizer Message:   {fit_result.message}")
    print(f"Final Residual Cost: {fit_result.cost:.4e}")
    print("\nFitted Parameters:")
    for name, val in fit_result.params.items():
        err = fit_result.errors.get(name)
        err_str = f" +/- {err:.3e}" if err is not None else ""
        print(f"  {name:12s} = {val:.4e}{err_str}")

    # 4. Plot Fit Comparison
    fig, axes = plot_fit(fit_result, exp_data)
    plt.savefig("fdpbd_fit_comparison.png", dpi=200)
    print("\nSaved fit comparison plot: fdpbd_fit_comparison.png")
    plt.close("all")


if __name__ == "__main__":
    main()
