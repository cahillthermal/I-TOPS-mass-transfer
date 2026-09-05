"""
Configuration dataclasses for FD-PBD simulations.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional, Union, Sequence
import numpy as np


@dataclass
class LaserBeamConfig:
    """
    Optical pump and probe beam parameters.

    Attributes
    ----------
    pump_power : float
        Incident pump laser power in Watts (W) measured at back-focal plane.
    wavelength : float
        Laser wavelength in meters (m). Default is 783 nm (783e-9 m).
    pump_radius : float
        1/e^2 intensity radius of focused pump beam (W0) in meters.
        E.g., 11.3 um for 5X objective, 2.7 um for 20X objective.
    probe_radius : Optional[float]
        1/e^2 intensity radius of focused probe beam (W1) in meters.
        Defaults to pump_radius if None.
    beam_offset : float
        Lateral offset / separation (r0) between pump and probe beam centers in meters.
        Optimal is typically ~ 1.1-1.2 * pump_radius (e.g. 12 um for 5X, 4 um for 20X).
    focal_length : float
        Focal length of objective lens (F) in meters.
        (40 mm for 5X lens, 10 mm for 20X lens).
    trans_objective : float
        Optical transmission of objective lens (0.90 for 5X, 0.70 for 20X).
    trans_window : float
        Optical transmission of vacuum/vapor chamber optical window (e.g. 0.92).
    trans_factor : float
        Absorption correction factor for transducer (typically 1.08 for Al).
    """

    pump_power: float = 3.0e-3
    wavelength: float = 783.0e-9
    pump_radius: float = 11.3e-6
    probe_radius: Optional[float] = None
    beam_offset: float = 12.0e-6
    focal_length: float = 40.0e-3
    trans_objective: float = 0.90
    trans_window: float = 0.92
    trans_factor: float = 1.08

    def __post_init__(self):
        if self.probe_radius is None:
            self.probe_radius = self.pump_radius

    @property
    def total_transmission(self) -> float:
        """Overall optical transmission to sample."""
        return self.trans_objective * self.trans_window

    def absorbed_power_amplitude(self, reflectance: float) -> float:
        """
        Calculate absorbed pump power amplitude parameter A (in Watts).
        A = (4/pi * P) * (1 - R) * trans_factor * total_transmission
        """
        return (4.0 / np.pi * self.pump_power) * (1.0 - reflectance) * self.trans_factor * self.total_transmission


@dataclass
class DetectorConfig:
    """
    Position-Sensitive Detector (PSD) and Lock-in detection configuration.

    Attributes
    ----------
    v0_multimeter : float
        PSD multimeter sum voltage (V0) in Volts (e.g. 0.830 V = 830 mV).
    gain2 : float
        Second-stage electronic gain of PSD (e.g. 1, 3, or 10).
    free_space_radius : float
        Free space laser beam radius in front of PSD in meters.
        (e.g. 1.55e-3 / 2 = 0.775 mm, or computed as wavelength * F / (pi * W0)).
    calib_factor : float
        Calibration geometric prefactor for PSD response (typically 0.65).
    """

    v0_multimeter: float = 830.0e-3
    gain2: float = 1.0
    free_space_radius: float = 1.55e-3 / 2.0
    calib_factor: float = 0.65

    def lockin_voltage_prefactor(self, focal_length: float) -> float:
        """
        Conversion prefactor from beam deflection angle (radians) to lock-in signal (microvolts uV).
        Prefactor = F * V0 * G2 / (0.65 * radii * sqrt(2)) * 1e6
        """
        return (
            focal_length
            * self.v0_multimeter
            * self.gain2
            / (self.calib_factor * self.free_space_radius * np.sqrt(2.0))
            * 1.0e6
        )


@dataclass
class MassTransportConfig:
    """
    Mass diffusion and sorption parameters in polymer film.

    Attributes
    ----------
    Dm : float
        Mass diffusivity in m^2/s (e.g., ~ 1e-12 m^2/s for PMMA).
    dndT_mass : float
        Refractive index change coefficient (dn/dC * dC/dT) in 1/K.
    aCET_mass : float
        Expansion coefficient due to water concentration ((dL/L/dC) * (dC/dT)) in 1/K.
    """

    Dm: float = 1.1e-12
    dndT_mass: float = -1.0e-4
    aCET_mass: float = -2.0e-5


@dataclass
class SimulationGrid:
    """
    Spatial frequency (k) and temporal modulation frequency (f) grid settings.

    Attributes
    ----------
    frequencies : np.ndarray
        Heating/modulation frequencies in Hz (1D array).
    k_vectors : Optional[np.ndarray]
        Spatial frequencies in 1/m (1D array).
        If None, automatically constructed based on pump_radius.
    k_min_ratio : float
        k_min = k_min_ratio / W0 (default 0.01).
    k_max_ratio : float
        k_max = k_max_ratio / W0 (default 5.0).
    k_step_ratio : float
        k_step = k_step_ratio / W0 (default 0.01).
    """

    frequencies: np.ndarray = field(
        default_factory=lambda: np.logspace(np.log10(1.0), np.log10(100.0e3), 100)
    )
    k_vectors: Optional[np.ndarray] = None
    k_min_ratio: float = 0.01
    k_max_ratio: float = 5.0
    k_step_ratio: float = 0.01

    def get_k_grid(self, pump_radius: float) -> np.ndarray:
        """Get or generate spatial frequency grid k (1/m)."""
        if self.k_vectors is not None:
            return np.asarray(self.k_vectors, dtype=float)
        k_min = self.k_min_ratio / pump_radius
        k_max = self.k_max_ratio / pump_radius
        k_step = self.k_step_ratio / pump_radius
        return np.arange(k_min, k_max + k_step * 0.5, k_step)


@dataclass
class ExperimentConfig:
    """
    Complete configuration bundling optical, detection, mass transport, and grid setups.
    """

    beam: LaserBeamConfig = field(default_factory=LaserBeamConfig)
    detector: DetectorConfig = field(default_factory=DetectorConfig)
    mass: MassTransportConfig = field(default_factory=MassTransportConfig)
    grid: SimulationGrid = field(default_factory=SimulationGrid)
