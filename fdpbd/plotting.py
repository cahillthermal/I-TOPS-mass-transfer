"""
Visualization and plotting routines for FD-PBD simulations, displacements, and fits.
"""

from __future__ import annotations
from typing import Optional, Union, Tuple, Sequence
import numpy as np
import matplotlib.pyplot as plt

from .simulator import SimulationResult
from .calibration import ExperimentData
from .fitting import FitResult


def plot_bode(
    sim_result: SimulationResult,
    exp_data: Optional[ExperimentData] = None,
    title: Optional[str] = "FD-PBD Frequency Response (Bode Plot)",
    figsize: Tuple[float, float] = (10.0, 4.5),
) -> Tuple[plt.Figure, Tuple[plt.Axes, plt.Axes]]:
    """
    Plot signal amplitude and phase as a function of modulation frequency.
    """
    fig, (ax_amp, ax_phase) = plt.subplots(1, 2, figsize=figsize)

    f = sim_result.frequencies
    ax_amp.semilogx(f, sim_result.amplitude, "-k", label="Model Total", lw=2)
    if sim_result.signal_heat is not None:
        ax_amp.semilogx(f, np.abs(sim_result.signal_heat), "--b", label="Heat Only", lw=1.5)
    if sim_result.signal_mass is not None:
        ax_amp.semilogx(f, np.abs(sim_result.signal_mass), ":r", label="Mass Only", lw=1.5)

    if exp_data is not None:
        ax_amp.semilogx(exp_data.frequencies, exp_data.amplitude, "ok", label="Experiment", ms=4)

    ax_amp.set_xlabel("Modulation Frequency (Hz)")
    ax_amp.set_ylabel(r"Amplitude ($\mu$V)")
    ax_amp.grid(True, which="both", ls="--", alpha=0.5)
    ax_amp.legend()

    ax_phase.semilogx(f, sim_result.phase, "-g", label="Model Phase", lw=2)
    if sim_result.signal_heat is not None:
        ax_phase.semilogx(f, np.angle(sim_result.signal_heat, deg=True), "--b", label="Heat Phase", lw=1.5)
    if sim_result.signal_mass is not None:
        ax_phase.semilogx(f, np.angle(sim_result.signal_mass, deg=True), ":r", label="Mass Phase", lw=1.5)

    if exp_data is not None:
        ax_phase.semilogx(exp_data.frequencies, exp_data.phase_deg, "og", label="Experiment", ms=4)

    ax_phase.set_xlabel("Modulation Frequency (Hz)")
    ax_phase.set_ylabel("Phase (degrees)")
    ax_phase.grid(True, which="both", ls="--", alpha=0.5)
    ax_phase.legend()

    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig, (ax_amp, ax_phase)


def plot_inphase_outphase(
    sim_result: SimulationResult,
    exp_data: Optional[ExperimentData] = None,
    title: Optional[str] = "FD-PBD Lock-in Quadratures",
    figsize: Tuple[float, float] = (10.0, 4.5),
) -> Tuple[plt.Figure, Tuple[plt.Axes, plt.Axes]]:
    """
    Plot in-phase (Vin) and out-of-phase (Vout) signals as a function of modulation frequency.
    """
    fig, (ax_in, ax_out) = plt.subplots(1, 2, figsize=figsize)

    f = sim_result.frequencies
    ax_in.semilogx(f, sim_result.Vin, "-b", label="Model Total", lw=2)
    if sim_result.signal_heat is not None:
        ax_in.semilogx(f, sim_result.Vin_heat, "--b", label="Heat", lw=1.5, alpha=0.7)
    if sim_result.signal_mass is not None:
        ax_in.semilogx(f, sim_result.Vin_mass, ":b", label="Mass", lw=1.5, alpha=0.7)

    if exp_data is not None:
        ax_in.semilogx(exp_data.frequencies, exp_data.Vin, "ob", label="Experiment", ms=4)

    ax_in.set_xlabel("Modulation Frequency (Hz)")
    ax_in.set_ylabel(r"$V_{\mathrm{in}}$ ($\mu$V)")
    ax_in.grid(True, which="both", ls="--", alpha=0.5)
    ax_in.legend()

    ax_out.semilogx(f, sim_result.Vout, "-r", label="Model Total", lw=2)
    if sim_result.signal_heat is not None:
        ax_out.semilogx(f, sim_result.Vout_heat, "--r", label="Heat", lw=1.5, alpha=0.7)
    if sim_result.signal_mass is not None:
        ax_out.semilogx(f, sim_result.Vout_mass, ":r", label="Mass", lw=1.5, alpha=0.7)

    if exp_data is not None:
        ax_out.semilogx(exp_data.frequencies, exp_data.Vout, "or", label="Experiment", ms=4)

    ax_out.set_xlabel("Modulation Frequency (Hz)")
    ax_out.set_ylabel(r"$V_{\mathrm{out}}$ ($\mu$V)")
    ax_out.grid(True, which="both", ls="--", alpha=0.5)
    ax_out.legend()

    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig, (ax_in, ax_out)


def plot_all(
    sim_result: SimulationResult,
    exp_data: Optional[ExperimentData] = None,
    title: str = "FD-PBD Measurement & Simulation",
    figsize: Tuple[float, float] = (11.0, 8.5),
) -> Tuple[plt.Figure, Sequence[plt.Axes]]:
    """
    Generate 4-panel comprehensive figure:
    (a) Amplitude, (b) Phase, (c) In-Phase (Vin), (d) Out-of-Phase (Vout).
    """
    fig, ((ax_amp, ax_phase), (ax_in, ax_out)) = plt.subplots(2, 2, figsize=figsize)
    f = sim_result.frequencies

    # (a) Amplitude
    ax_amp.semilogx(f, sim_result.amplitude, "-k", label="Sim Total", lw=2)
    if sim_result.signal_heat is not None:
        ax_amp.semilogx(f, np.abs(sim_result.signal_heat), "--k", label="Heat", lw=1.5)
    if sim_result.signal_mass is not None:
        ax_amp.semilogx(f, np.abs(sim_result.signal_mass), ":k", label="Mass", lw=1.5)
    if exp_data is not None:
        ax_amp.semilogx(exp_data.frequencies, exp_data.amplitude, "ok", label="Exp Data", ms=4)
    ax_amp.set_xlabel("Frequency (Hz)")
    ax_amp.set_ylabel(r"Amplitude ($\mu$V)")
    ax_amp.set_title("(a) Amplitude")
    ax_amp.grid(True, which="both", ls="--", alpha=0.5)
    ax_amp.legend()

    # (b) Phase
    ax_phase.semilogx(f, sim_result.phase, "-g", label="Sim Total", lw=2)
    if sim_result.signal_heat is not None:
        ax_phase.semilogx(f, np.angle(sim_result.signal_heat, deg=True), "--g", label="Heat", lw=1.5)
    if sim_result.signal_mass is not None:
        ax_phase.semilogx(f, np.angle(sim_result.signal_mass, deg=True), ":g", label="Mass", lw=1.5)
    if exp_data is not None:
        ax_phase.semilogx(exp_data.frequencies, exp_data.phase_deg, "og", label="Exp Data", ms=4)
    ax_phase.set_xlabel("Frequency (Hz)")
    ax_phase.set_ylabel("Phase (degrees)")
    ax_phase.set_title("(b) Phase")
    ax_phase.grid(True, which="both", ls="--", alpha=0.5)
    ax_phase.legend()

    # (c) In-Phase (Vin)
    ax_in.semilogx(f, sim_result.Vin, "-b", label="Sim Total", lw=2)
    if sim_result.signal_heat is not None:
        ax_in.semilogx(f, sim_result.Vin_heat, "--b", label="Heat", lw=1.5)
    if sim_result.signal_mass is not None:
        ax_in.semilogx(f, sim_result.Vin_mass, ":b", label="Mass", lw=1.5)
    if exp_data is not None:
        ax_in.semilogx(exp_data.frequencies, exp_data.Vin, "ob", label="Exp Data", ms=4)
    ax_in.set_xlabel("Frequency (Hz)")
    ax_in.set_ylabel(r"$V_{\mathrm{in}}$ ($\mu$V)")
    ax_in.set_title(r"(c) In-Phase $V_{\mathrm{in}}$")
    ax_in.grid(True, which="both", ls="--", alpha=0.5)
    ax_in.legend()

    # (d) Out-of-Phase (Vout)
    ax_out.semilogx(f, sim_result.Vout, "-r", label="Sim Total", lw=2)
    if sim_result.signal_heat is not None:
        ax_out.semilogx(f, sim_result.Vout_heat, "--r", label="Heat", lw=1.5)
    if sim_result.signal_mass is not None:
        ax_out.semilogx(f, sim_result.Vout_mass, ":r", label="Mass", lw=1.5)
    if exp_data is not None:
        ax_out.semilogx(exp_data.frequencies, exp_data.Vout, "or", label="Exp Data", ms=4)
    ax_out.set_xlabel("Frequency (Hz)")
    ax_out.set_ylabel(r"$V_{\mathrm{out}}$ ($\mu$V)")
    ax_out.set_title(r"(d) Out-of-Phase $V_{\mathrm{out}}$")
    ax_out.grid(True, which="both", ls="--", alpha=0.5)
    ax_out.legend()

    fig.suptitle(title, fontsize=13)
    fig.tight_layout()
    return fig, (ax_amp, ax_phase, ax_in, ax_out)


def plot_displacement_breakdown(
    sim_result: SimulationResult,
    k_index: int = 0,
    figsize: Tuple[float, float] = (12.0, 5.0),
) -> Tuple[plt.Figure, Tuple[plt.Axes, plt.Axes]]:
    """
    Plot the frequency dependence of the 5 displacement mechanisms (Z0 to Z4)
    reproducing Figure 2 from Xie et al., RSI 2018.
    """
    fig, (ax_in, ax_out) = plt.subplots(1, 2, figsize=figsize)
    f = sim_result.frequencies
    disp = sim_result.displacements

    z0 = disp.Z0_m[k_index, :]
    z1 = disp.Z1_m[k_index, :]
    z2 = disp.Z2_m[k_index, :]
    z3 = disp.Z3_m[k_index, :]
    z4 = disp.Z4_m[k_index, :]

    # In-Phase
    ax_in.semilogx(f, np.real(z0) * 1e12, label=r"$Z_0$ (Fresnel / Optics)", lw=1.8)
    ax_in.semilogx(f, np.real(z1) * 1e12, label=r"$Z_1$ (Al Expansion)", lw=1.8)
    ax_in.semilogx(f, np.real(z2) * 1e12, label=r"$Z_2$ (Substrate Stress)", lw=1.8)
    ax_in.semilogx(f, np.real(z3) * 1e12, label=r"$Z_3$ (Substrate Expansion)", lw=1.8)
    ax_in.semilogx(f, np.real(z4) * 1e12, label=r"$Z_4$ (Ambient Thermo-optic)", lw=1.8)
    ax_in.set_xlabel("Frequency (Hz)")
    ax_in.set_ylabel(r"In-Phase Displacement (pm/$\Delta$T)")
    ax_in.set_title("In-Phase Displacements")
    ax_in.grid(True, which="both", ls="--", alpha=0.5)
    ax_in.legend()

    # Out-of-Phase
    ax_out.semilogx(f, np.imag(z0) * 1e12, label=r"$Z_0$ (Fresnel / Optics)", lw=1.8)
    ax_out.semilogx(f, np.imag(z1) * 1e12, label=r"$Z_1$ (Al Expansion)", lw=1.8)
    ax_out.semilogx(f, np.imag(z2) * 1e12, label=r"$Z_2$ (Substrate Stress)", lw=1.8)
    ax_out.semilogx(f, np.imag(z3) * 1e12, label=r"$Z_3$ (Substrate Expansion)", lw=1.8)
    ax_out.semilogx(f, np.imag(z4) * 1e12, label=r"$Z_4$ (Ambient Thermo-optic)", lw=1.8)
    ax_out.set_xlabel("Frequency (Hz)")
    ax_out.set_ylabel(r"Out-of-Phase Displacement (pm/$\Delta$T)")
    ax_out.set_title("Out-of-Phase Displacements")
    ax_out.grid(True, which="both", ls="--", alpha=0.5)
    ax_out.legend()

    fig.tight_layout()
    return fig, (ax_in, ax_out)


def plot_fit(
    fit_result: FitResult,
    exp_data: ExperimentData,
    figsize: Tuple[float, float] = (11.0, 8.5),
) -> Tuple[plt.Figure, Sequence[plt.Axes]]:
    """
    Plot fitted simulation against experimental data with parameter legend.
    """
    param_str = ", ".join([f"{k} = {v:.3e}" for k, v in fit_result.params.items()])
    title = f"FD-PBD Fit Result: {param_str}"
    return plot_all(
        sim_result=fit_result.simulation,
        exp_data=exp_data,
        title=title,
        figsize=figsize,
    )
