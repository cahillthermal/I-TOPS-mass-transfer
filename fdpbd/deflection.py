"""
Beam deflection angle integration and lock-in voltage conversion.
"""

from __future__ import annotations
from typing import Union, Tuple, Optional
import numpy as np
from scipy.special import j1
from .config import LaserBeamConfig, DetectorConfig


def calculate_beam_deflection(
    displacement_m: np.ndarray,
    k_vectors: np.ndarray,
    beam_config: LaserBeamConfig,
    detector_config: DetectorConfig,
    reflectance: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Integrate equivalent measurable surface displacement into average probe beam deflection angle
    and convert to lock-in voltage signal (microvolts uV).

    Equation (23) and (25) from Xie et al., RSI 2018:
    theta(omega) = 8*pi^2 * integral[ P(k) * Z_m(k, omega) * S(k) * J1(2*pi*k*r0) * k^2 dk ]
    V_lockin = theta * prefactor

    Parameters
    ----------
    displacement_m : np.ndarray
        Measurable equivalent surface displacement matrix of shape (Nk, Nf).
    k_vectors : np.ndarray
        Spatial frequencies in 1/m (1D array of length Nk).
    beam_config : LaserBeamConfig
        Laser beam parameters (pump_radius, probe_radius, offset, power, focal length).
    detector_config : DetectorConfig
        Photodetector and lock-in parameters (V0, gain, beam radius, etc.).
    reflectance : float
        Nominal power reflectance of the sample stack (used to compute absorbed power amplitude A).

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        (theta_ave, signal_uV)
        - theta_ave : complex beam deflection angle in radians, shape (Nf,)
        - signal_uV : complex lock-in voltage in microvolts (uV), shape (Nf,)
    """
    k = np.asarray(k_vectors, dtype=float).reshape(-1, 1)  # (Nk, 1)
    W0 = beam_config.pump_radius
    W1 = beam_config.probe_radius if beam_config.probe_radius is not None else W0
    r0 = beam_config.beam_offset
    A = beam_config.absorbed_power_amplitude(reflectance)

    # Pump spatial distribution (Hankel transform)
    Pk = A * np.exp(-np.pi**2 * (k**2) * (W0**2) / 2.0)  # (Nk, 1)

    # Probe spatial distribution (Hankel transform)
    Sk = np.exp(-np.pi**2 * (k**2) * (W1**2) / 2.0)  # (Nk, 1)

    # Bessel function of the first kind J1(2*pi*k*r0)
    J1_term = j1(2.0 * np.pi * k * r0)  # (Nk, 1)

    # Integrand over k: 8*pi^2 * Pk * Z_m * Sk * J1 * k^2
    integrand = 8.0 * (np.pi**2) * (Pk * displacement_m) * Sk * J1_term * (k**2)  # (Nk, Nf)

    # Vectorized numerical integration over spatial frequency k
    k_1d = k.ravel()
    theta_ave = np.trapezoid(integrand, k_1d, axis=0)  # (Nf,)

    # Lock-in voltage conversion
    prefactor = detector_config.lockin_voltage_prefactor(beam_config.focal_length)
    signal_uV = theta_ave * prefactor  # (Nf,)

    return theta_ave, signal_uV
