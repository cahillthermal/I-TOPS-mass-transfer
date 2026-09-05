"""
Material definitions and database for FD-PBD simulations.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional, Union
import numpy as np


@dataclass
class Material:
    """
    Physical properties of a material layer in the FD-PBD model.

    Attributes
    ----------
    name : str
        Descriptive name of the material.
    kth : float
        Thermal conductivity in W/(m*K).
    Cv : float
        Volumetric heat capacity in J/(m^3*K).
    n : complex or float
        Complex refractive index at probe/pump wavelength.
    dndT : complex or float
        Thermo-optic coefficient (dn/dT) in 1/K.
    aCET : float
        Linear coefficient of thermal expansion (CTE) in 1/K.
    nu : float
        Poisson's ratio (dimensionless).
    Y : float
        Young's modulus in Pa.
    Vl : float
        Longitudinal acoustic sound velocity in m/s (used for substrate).
    Vt : float
        Transverse acoustic sound velocity in m/s (used for substrate).
    eta : float
        Thermal conductivity anisotropy ratio (kx/kz), default 1.0 (isotropic).
    dndT_mass : float
        Mass-induced thermo-optic equivalent coefficient (dn/dC * dC/dT) in 1/K.
    aCET_mass : float
        Mass-induced thermal expansion equivalent coefficient (dL/L/dC * dC/dT) in 1/K.
    Dm : float
        Mass diffusivity in m^2/s (for permeable layers, e.g. polymers).
    """

    name: str
    kth: float = 0.0
    Cv: float = 0.0
    n: Union[complex, float] = 1.0
    dndT: Union[complex, float] = 0.0
    aCET: float = 0.0
    nu: float = 0.0
    Y: float = 0.0
    Vl: float = 0.0
    Vt: float = 0.0
    eta: float = 1.0
    dndT_mass: float = 0.0
    aCET_mass: float = 0.0
    Dm: float = 0.0

    @property
    def thermal_diffusivity(self) -> float:
        """Thermal diffusivity D = kth / Cv in m^2/s."""
        if self.Cv <= 0:
            return 0.0
        return self.kth / self.Cv

    def copy_with(self, **kwargs) -> Material:
        """Create a modified copy of the material."""
        d = {
            "name": self.name,
            "kth": self.kth,
            "Cv": self.Cv,
            "n": self.n,
            "dndT": self.dndT,
            "aCET": self.aCET,
            "nu": self.nu,
            "Y": self.Y,
            "Vl": self.Vl,
            "Vt": self.Vt,
            "eta": self.eta,
            "dndT_mass": self.dndT_mass,
            "aCET_mass": self.aCET_mass,
            "Dm": self.Dm,
        }
        d.update(kwargs)
        return Material(**d)


def create_aluminum_thermo_optic(
    n_al: complex = 1.96 + 7.10j,
    dRdT: float = 1.0e-4,
    dphi_dT: float = 3.0e-5,
    R_ref: float = 0.87,
) -> complex:
    """
    Compute the complex dn/dT of Aluminum transducer from TDTR/FDTR calibration.

    Parameters
    ----------
    n_al : complex
        Complex refractive index of aluminum.
    dRdT : float
        Temperature derivative of reflectivity dR/dT in 1/K.
    dphi_dT : float
        Temperature derivative of Fresnel reflection phase in rad/K.
    R_ref : float
        Nominal reflectance.

    Returns
    -------
    dndT_Al : complex
    """
    dMdT = 0.5 * dRdT / R_ref
    dndT_Al = 0.5 * (n_al**2 - 1.0) * (dMdT + 1j * dphi_dT)
    return dndT_Al


# Standard preset materials from Table I and manuscript (Xie et al., RSI 2018)

AIR = Material(
    name="Air",
    kth=0.026,
    Cv=0.0012 * 1e6,  # 1.2e3 J/m^3/K
    n=1.0,
    dndT=-9.1e-7,
    aCET=0.0,
    nu=0.0,
    Y=0.0,
)

WATER_VAPOR = Material(
    name="Water Vapor",
    kth=0.019,
    Cv=40.0e-6 * 1e6,  # 40 J/m^3/K (at 21 Torr, 23 C)
    n=1.0,
    dndT=-2.1e-8,
    aCET=0.0,
    nu=0.0,
    Y=0.0,
)

PMMA = Material(
    name="PMMA",
    kth=0.19,
    Cv=1.60 * 1e6,  # 1.6e6 J/m^3/K
    n=1.49,
    dndT=-0.78e-4,
    aCET=33.0e-6,
    nu=0.37,
    Y=2.5e9,
    Dm=1.1e-12,
    dndT_mass=-1.0e-4,
    aCET_mass=-2.0e-5,
)

ALUMINUM = Material(
    name="Aluminum",
    kth=160.0,
    Cv=2.42 * 1e6,  # 2.42e6 J/m^3/K
    n=1.96 + 7.10j,
    dndT=create_aluminum_thermo_optic(),
    aCET=23.0e-6,
    nu=0.33,
    Y=69.0e9,
)

FUSED_SILICA = Material(
    name="Fused Silica",
    kth=1.32,
    Cv=1.63 * 1e6,  # 1.63e6 J/m^3/K
    n=1.4536,
    dndT=1.0e-5,
    aCET=0.54e-6,
    nu=0.17,
    Y=73.0e9,
    Vl=5970.0,
    Vt=3770.0,
)

POLYIMIDE = Material(
    name="Polyimide",
    kth=0.20,
    Cv=1.50 * 1e6,
    n=1.70,
    dndT=-1.0e-4,
    aCET=35.0e-6,
    nu=0.34,
    Y=3.0e9,
    Dm=7.0e-13,
    dndT_mass=-1.0e-4,
    aCET_mass=-2.0e-5,
)

PVP = Material(
    name="PVP",
    kth=0.20,
    Cv=1.50 * 1e6,
    n=1.52,
    dndT=-1.0e-4,
    aCET=35.0e-6,
    nu=0.35,
    Y=2.0e9,
    Dm=1.7e-12,
    dndT_mass=-1.0e-4,
    aCET_mass=-2.0e-5,
)

CELLULOSE_ACETATE = Material(
    name="Cellulose Acetate",
    kth=0.20,
    Cv=1.50 * 1e6,
    n=1.47,
    dndT=-1.0e-4,
    aCET=35.0e-6,
    nu=0.35,
    Y=2.0e9,
    Dm=2.6e-11,
    dndT_mass=-1.0e-4,
    aCET_mass=-2.0e-5,
)

PIP_TMC = Material(
    name="PIP/TMC Polyamide",
    kth=0.20,
    Cv=1.50 * 1e6,
    n=1.55,
    dndT=-1.0e-4,
    aCET=35.0e-6,
    nu=0.35,
    Y=2.5e9,
    Dm=9.0e-11,
    dndT_mass=-1.0e-4,
    aCET_mass=-2.0e-5,
)
