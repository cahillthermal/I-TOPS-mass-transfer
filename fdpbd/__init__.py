"""
FD-PBD (Frequency-Domain Probe Beam Deflection) Heat and Mass Transfer Simulation & Fitting.

Based on the physical model described in:
Xie, Dennison, Shin, Diao, Cahill, Rev. Sci. Instrum. 89, 104904 (2018).
"""

from .materials import (
    Material,
    AIR,
    WATER_VAPOR,
    PMMA,
    ALUMINUM,
    FUSED_SILICA,
    POLYIMIDE,
    PVP,
    CELLULOSE_ACETATE,
    PIP_TMC,
    create_aluminum_thermo_optic,
)
from .sample import Layer, Sample
from .config import (
    LaserBeamConfig,
    DetectorConfig,
    MassTransportConfig,
    SimulationGrid,
    ExperimentConfig,
)
from .simulator import FDPBDSimulator, SimulationResult
from .thermal import solve_thermal_field, calculate_steady_state_heating, calculate_temperature_rise
from .optics import solve_fresnel_reflection, solve_multilayer_optics_perturbed
from .mechanics import solve_displacements, DisplacementBreakdown
from .deflection import calculate_beam_deflection
from .calibration import ExperimentData, CalibrationData
from .fitting import FDPBDFitter, FitParameter, FitResult
from .plotting import (
    plot_bode,
    plot_inphase_outphase,
    plot_displacement_breakdown,
    plot_fit,
    plot_all,
)

__version__ = "1.0.0"

__all__ = [
    "Material",
    "AIR",
    "WATER_VAPOR",
    "PMMA",
    "ALUMINUM",
    "FUSED_SILICA",
    "POLYIMIDE",
    "PVP",
    "CELLULOSE_ACETATE",
    "PIP_TMC",
    "create_aluminum_thermo_optic",
    "Layer",
    "Sample",
    "LaserBeamConfig",
    "DetectorConfig",
    "MassTransportConfig",
    "SimulationGrid",
    "ExperimentConfig",
    "FDPBDSimulator",
    "SimulationResult",
    "solve_thermal_field",
    "calculate_steady_state_heating",
    "calculate_temperature_rise",
    "solve_fresnel_reflection",
    "solve_multilayer_optics_perturbed",
    "solve_displacements",
    "DisplacementBreakdown",
    "calculate_beam_deflection",
    "ExperimentData",
    "CalibrationData",
    "FDPBDFitter",
    "FitParameter",
    "FitResult",
    "plot_bode",
    "plot_inphase_outphase",
    "plot_displacement_breakdown",
    "plot_fit",
    "plot_all",
]
