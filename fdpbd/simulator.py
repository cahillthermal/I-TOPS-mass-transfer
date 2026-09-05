"""
High-level simulation engine for Frequency-Domain Probe Beam Deflection (FD-PBD).
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional, Union, Sequence, Tuple
import numpy as np

from .sample import Sample, Layer
from .config import (
    LaserBeamConfig,
    DetectorConfig,
    MassTransportConfig,
    SimulationGrid,
    ExperimentConfig,
)
from .thermal import solve_thermal_field, calculate_steady_state_heating
from .optics import solve_fresnel_reflection, solve_multilayer_optics_perturbed, OpticalReflectionResult
from .mechanics import solve_displacements, DisplacementBreakdown
from .deflection import calculate_beam_deflection


@dataclass
class SimulationResult:
    """
    Comprehensive output from an FD-PBD simulation.

    Attributes
    ----------
    frequencies : np.ndarray
        Heating/modulation frequencies in Hz, shape (Nf,).
    signal_total : np.ndarray
        Total complex lock-in deflection signal in microvolts (uV), shape (Nf,).
    signal_heat : Optional[np.ndarray]
        Thermal-only deflection signal in uV, shape (Nf,).
    signal_mass : Optional[np.ndarray]
        Mass-transport-only deflection signal in uV, shape (Nf,).
    displacements : DisplacementBreakdown
        Full 2D breakdown of displacement components Z0..Z4 over (k, f).
    theta_total : np.ndarray
        Total complex beam deflection angle in radians, shape (Nf,).
    reflectance : float
        Nominal sample stack power reflectance.
    absorbed_power_A : float
        Absorbed pump power amplitude in Watts.
    steady_state_dT : float
        Maximum steady-state temperature rise in Kelvin.
    """

    frequencies: np.ndarray
    signal_total: np.ndarray
    signal_heat: Optional[np.ndarray]
    signal_mass: Optional[np.ndarray]
    displacements: DisplacementBreakdown
    theta_total: np.ndarray
    reflectance: float
    absorbed_power_A: float
    steady_state_dT: float

    @property
    def Vin(self) -> np.ndarray:
        """In-phase lock-in signal (uV)."""
        return np.real(self.signal_total)

    @property
    def Vout(self) -> np.ndarray:
        """Out-of-phase lock-in signal (uV)."""
        return np.imag(self.signal_total)

    @property
    def amplitude(self) -> np.ndarray:
        """Signal amplitude (uV)."""
        return np.abs(self.signal_total)

    @property
    def phase(self) -> np.ndarray:
        """Signal phase in degrees."""
        return np.angle(self.signal_total, deg=True)

    @property
    def Vin_heat(self) -> Optional[np.ndarray]:
        return np.real(self.signal_heat) if self.signal_heat is not None else None

    @property
    def Vout_heat(self) -> Optional[np.ndarray]:
        return np.imag(self.signal_heat) if self.signal_heat is not None else None

    @property
    def Vin_mass(self) -> Optional[np.ndarray]:
        return np.real(self.signal_mass) if self.signal_mass is not None else None

    @property
    def Vout_mass(self) -> Optional[np.ndarray]:
        return np.imag(self.signal_mass) if self.signal_mass is not None else None


class FDPBDSimulator:
    """
    Simulation orchestrator for FD-PBD experiments.
    """

    @staticmethod
    def simulate(
        sample: Sample,
        frequencies: Optional[np.ndarray] = None,
        k_vectors: Optional[np.ndarray] = None,
        beam_config: Optional[LaserBeamConfig] = None,
        detector_config: Optional[DetectorConfig] = None,
        mass_config: Optional[MassTransportConfig] = None,
        grid_config: Optional[SimulationGrid] = None,
        experiment_config: Optional[ExperimentConfig] = None,
        decompose_components: bool = True,
    ) -> SimulationResult:
        """
        Execute full FD-PBD simulation.

        Parameters
        ----------
        sample : Sample
            Multilayer stack definition.
        frequencies : Optional[np.ndarray]
            Array of modulation frequencies (Hz).
        k_vectors : Optional[np.ndarray]
            Spatial frequency array (1/m).
        beam_config : Optional[LaserBeamConfig]
            Optical beam parameters.
        detector_config : Optional[DetectorConfig]
            PSD / lock-in parameters.
        mass_config : Optional[MassTransportConfig]
            Mass transfer parameters.
        grid_config : Optional[SimulationGrid]
            Simulation grid options.
        experiment_config : Optional[ExperimentConfig]
            Pre-bundled configuration object (overrides individual config objects if passed).
        decompose_components : bool
            Whether to additionally compute isolated thermal-only and mass-only signals (default True).

        Returns
        -------
        SimulationResult
        """
        # Resolve configurations
        if experiment_config is not None:
            if beam_config is None:
                beam_config = experiment_config.beam
            if detector_config is None:
                detector_config = experiment_config.detector
            if mass_config is None:
                mass_config = experiment_config.mass
            if grid_config is None:
                grid_config = experiment_config.grid

        if beam_config is None:
            beam_config = LaserBeamConfig()
        if detector_config is None:
            detector_config = DetectorConfig()
        if grid_config is None:
            grid_config = SimulationGrid()

        if frequencies is None:
            frequencies = np.asarray(grid_config.frequencies, dtype=float)
        else:
            frequencies = np.asarray(frequencies, dtype=float)

        if k_vectors is None:
            k_vectors = grid_config.get_k_grid(beam_config.pump_radius)
        else:
            k_vectors = np.asarray(k_vectors, dtype=float)

        # 1. Base multilayer optics and nominal reflectance
        R_nom, _, _ = solve_fresnel_reflection(
            sample=sample,
            theta_deg=0.0,
            wavelength=beam_config.wavelength,
        )

        A_pump = beam_config.absorbed_power_amplitude(R_nom)

        # 2. Steady-state temperature rise calculation
        # Steady state heating: A_SS = (2*A + 2*A)*pi/4 = pi*A
        A_SS = np.pi * A_pump
        dT_ss = calculate_steady_state_heating(
            sample=sample,
            k_vectors=k_vectors,
            pump_power_absorbed=A_SS,
            pump_radius=beam_config.pump_radius,
            heating_layer_idx=sample.heating_layer_idx,
            temp_layer_idx=sample.heating_layer_idx,
            z_offset=0.0,
        )

        # 3. Calculate full surface displacements (Z0 to Z4)
        disp_breakdown = solve_displacements(
            sample=sample,
            frequencies=frequencies,
            k_vectors=k_vectors,
            mass_config=mass_config,
            wavelength=beam_config.wavelength,
            heating_layer_idx=sample.heating_layer_idx,
        )

        # 4. Integrate total deflection angle and convert to lock-in signal
        theta_total, signal_total = calculate_beam_deflection(
            displacement_m=disp_breakdown.Z_tot_m,
            k_vectors=k_vectors,
            beam_config=beam_config,
            detector_config=detector_config,
            reflectance=R_nom,
        )

        signal_heat = None
        signal_mass = None

        # 5. Optional decomposition into heat-only and mass-only signals
        if decompose_components and sample.heating_layer_idx > 1:
            poly_idx = sample.heating_layer_idx - 1
            # Thermal only: zero out mass transport
            zero_mass = MassTransportConfig(
                Dm=mass_config.Dm if mass_config else sample.layers[poly_idx].Dm,
                dndT_mass=0.0,
                aCET_mass=0.0,
            )
            disp_heat = solve_displacements(
                sample=sample,
                frequencies=frequencies,
                k_vectors=k_vectors,
                mass_config=zero_mass,
                wavelength=beam_config.wavelength,
                heating_layer_idx=sample.heating_layer_idx,
            )
            _, signal_heat = calculate_beam_deflection(
                displacement_m=disp_heat.Z_tot_m,
                k_vectors=k_vectors,
                beam_config=beam_config,
                detector_config=detector_config,
                reflectance=R_nom,
            )

            # Mass only: total minus thermal
            signal_mass = signal_total - signal_heat

        return SimulationResult(
            frequencies=frequencies,
            signal_total=signal_total,
            signal_heat=signal_heat,
            signal_mass=signal_mass,
            displacements=disp_breakdown,
            theta_total=theta_total,
            reflectance=R_nom,
            absorbed_power_A=A_pump,
            steady_state_dT=dT_ss,
        )
