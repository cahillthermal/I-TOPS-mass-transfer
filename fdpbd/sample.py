"""
Layer and Sample multilayer stack representations for FD-PBD.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Optional, Union, Sequence
import numpy as np
from .materials import Material, WATER_VAPOR, AIR, PMMA, ALUMINUM, FUSED_SILICA


@dataclass
class Layer:
    """
    A single layer in the sample stack.

    Parameters
    ----------
    material : Material
        The base material providing physical and optical properties.
    thickness : float
        Layer thickness in meters (m). For semi-infinite layers (e.g. ambient or substrate),
        use a large number like 1.0 (1 m) or np.inf.
    name : Optional[str]
        Optional custom name override.
    overrides : dict
        Optional property overrides for this specific layer instance.
    """

    material: Material
    thickness: float
    name: Optional[str] = None
    overrides: dict = field(default_factory=dict)

    def __post_init__(self):
        if self.name is None:
            self.name = self.material.name

    def get_property(self, prop_name: str):
        """Get property from overrides if present, otherwise from material."""
        if prop_name in self.overrides:
            return self.overrides[prop_name]
        return getattr(self.material, prop_name)

    @property
    def kth(self) -> float:
        return self.get_property("kth")

    @property
    def Cv(self) -> float:
        return self.get_property("Cv")

    @property
    def n(self) -> Union[complex, float]:
        return self.get_property("n")

    @property
    def dndT(self) -> Union[complex, float]:
        return self.get_property("dndT")

    @property
    def aCET(self) -> float:
        return self.get_property("aCET")

    @property
    def nu(self) -> float:
        return self.get_property("nu")

    @property
    def Y(self) -> float:
        return self.get_property("Y")

    @property
    def Vl(self) -> float:
        return self.get_property("Vl")

    @property
    def Vt(self) -> float:
        return self.get_property("Vt")

    @property
    def eta(self) -> float:
        return self.get_property("eta")

    @property
    def dndT_mass(self) -> float:
        return self.get_property("dndT_mass")

    @property
    def aCET_mass(self) -> float:
        return self.get_property("aCET_mass")

    @property
    def Dm(self) -> float:
        return self.get_property("Dm")

    @property
    def D(self) -> float:
        """Thermal diffusivity in m^2/s."""
        Cv_val = self.Cv
        return self.kth / Cv_val if Cv_val > 0 else 0.0


@dataclass
class Sample:
    """
    Multilayer sample stack for FD-PBD thermal, mass, and optical modeling.

    Layers are indexed from 0 (top ambient/vapor) to N-1 (bottom substrate).
    The heating interface is defined at the top surface of `heating_layer_idx`.

    Parameters
    ----------
    layers : List[Layer]
        Sequence of layers from top (ambient) to bottom (substrate).
    heating_layer_idx : int
        0-based index of the layer whose top surface absorbs pump heat (typically the transducer layer, e.g. Al).
        Default is 2 for standard 4-layer stack [Vapor, Polymer, Al, SiO2] (corresponding to MATLAB nheat = 3).
    """

    layers: List[Layer] = field(default_factory=list)
    heating_layer_idx: int = 2

    def add_layer(
        self,
        material: Material,
        thickness: float,
        name: Optional[str] = None,
        **overrides,
    ) -> Sample:
        """Add a layer to the bottom of the current stack."""
        self.layers.append(Layer(material=material, thickness=thickness, name=name, overrides=overrides))
        return self

    def __len__(self) -> int:
        return len(self.layers)

    def __getitem__(self, idx: int) -> Layer:
        return self.layers[idx]

    @property
    def num_layers(self) -> int:
        return len(self.layers)

    @property
    def kth(self) -> np.ndarray:
        return np.array([layer.kth for layer in self.layers], dtype=float)

    @property
    def Cv(self) -> np.ndarray:
        return np.array([layer.Cv for layer in self.layers], dtype=float)

    @property
    def thickness(self) -> np.ndarray:
        return np.array([layer.thickness for layer in self.layers], dtype=float)

    @property
    def n_refrac(self) -> np.ndarray:
        return np.array([layer.n for layer in self.layers], dtype=complex)

    @property
    def dndT(self) -> np.ndarray:
        return np.array([layer.dndT for layer in self.layers], dtype=complex)

    @property
    def aCET(self) -> np.ndarray:
        return np.array([layer.aCET for layer in self.layers], dtype=float)

    @property
    def nu(self) -> np.ndarray:
        return np.array([layer.nu for layer in self.layers], dtype=float)

    @property
    def Y(self) -> np.ndarray:
        return np.array([layer.Y for layer in self.layers], dtype=float)

    @property
    def Vl(self) -> np.ndarray:
        return np.array([layer.Vl for layer in self.layers], dtype=float)

    @property
    def Vt(self) -> np.ndarray:
        return np.array([layer.Vt for layer in self.layers], dtype=float)

    @property
    def eta(self) -> np.ndarray:
        return np.array([layer.eta for layer in self.layers], dtype=float)

    @property
    def D(self) -> np.ndarray:
        """Thermal diffusivities for all layers in m^2/s."""
        return self.kth / self.Cv

    # 1-based index helper for MATLAB compatibility
    @property
    def nheat_matlab(self) -> int:
        """1-based MATLAB nheat index."""
        return self.heating_layer_idx + 1

    @classmethod
    def create_standard_fdpbd_sample(
        cls,
        polymer_thickness: float = 100e-9,
        al_thickness: float = 80e-9,
        polymer_material: Material = PMMA,
        ambient_material: Material = WATER_VAPOR,
        transducer_material: Material = ALUMINUM,
        substrate_material: Material = FUSED_SILICA,
    ) -> Sample:
        """
        Factory to create standard 4-layer FD-PBD sample:
        Layer 0: Ambient (Water Vapor / Air) [semi-infinite]
        Layer 1: Polymer Thin Film (e.g. PMMA, PI, PVP, CA)
        Layer 2: Metal Transducer (e.g. Al, heating layer)
        Layer 3: Substrate (e.g. Fused Silica) [semi-infinite]
        """
        sample = cls(heating_layer_idx=2)
        sample.add_layer(ambient_material, thickness=1.0, name="Ambient")
        sample.add_layer(polymer_material, thickness=polymer_thickness, name="Polymer")
        sample.add_layer(transducer_material, thickness=al_thickness, name="Transducer")
        sample.add_layer(substrate_material, thickness=1.0, name="Substrate")
        return sample

    @classmethod
    def create_control_sample(
        cls,
        al_thickness: float = 80e-9,
        ambient_material: Material = WATER_VAPOR,
        transducer_material: Material = ALUMINUM,
        substrate_material: Material = FUSED_SILICA,
    ) -> Sample:
        """
        Factory for 3-layer control sample without polymer:
        Layer 0: Ambient (Water Vapor / Air / Vacuum)
        Layer 1: Transducer (Al, heating layer)
        Layer 2: Substrate (Fused Silica)
        """
        sample = cls(heating_layer_idx=1)
        sample.add_layer(ambient_material, thickness=1.0, name="Ambient")
        sample.add_layer(transducer_material, thickness=al_thickness, name="Transducer")
        sample.add_layer(substrate_material, thickness=1.0, name="Substrate")
        return sample
