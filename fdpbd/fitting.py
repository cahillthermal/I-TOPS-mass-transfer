"""
Nonlinear parameter estimation and curve fitting for FD-PBD.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Callable, Union
import numpy as np
from scipy.optimize import minimize, least_squares

from .sample import Sample
from .config import LaserBeamConfig, DetectorConfig, MassTransportConfig, SimulationGrid, ExperimentConfig
from .simulator import FDPBDSimulator, SimulationResult
from .calibration import ExperimentData


@dataclass
class FitParameter:
    """
    Definition of a parameter to fit in the FD-PBD model.
    """

    name: str
    value: float
    bounds: Tuple[Optional[float], Optional[float]] = (None, None)
    scale: float = 1.0
    vary: bool = True

    @property
    def normalized_value(self) -> float:
        return self.value / self.scale

    def unnormalize(self, x_norm: float) -> float:
        return x_norm * self.scale


@dataclass
class FitResult:
    """
    Results of an FD-PBD model fit to experimental data.
    """

    params: Dict[str, float]
    errors: Dict[str, Optional[float]]
    cost: float
    success: bool
    message: str
    simulation: SimulationResult


class FDPBDFitter:
    """
    Flexible curve fitting engine for FD-PBD data.
    """

    def __init__(
        self,
        sample: Sample,
        exp_data: ExperimentData,
        beam_config: Optional[LaserBeamConfig] = None,
        detector_config: Optional[DetectorConfig] = None,
        mass_config: Optional[MassTransportConfig] = None,
        k_vectors: Optional[np.ndarray] = None,
    ):
        self.sample = sample
        self.exp_data = exp_data
        self.beam_config = beam_config or LaserBeamConfig()
        self.detector_config = detector_config or DetectorConfig()
        self.mass_config = mass_config or MassTransportConfig()
        self.k_vectors = k_vectors
        self.params: Dict[str, FitParameter] = {}

    def add_parameter(
        self,
        name: str,
        value: Optional[float] = None,
        bounds: Tuple[Optional[float], Optional[float]] = (None, None),
        scale: Optional[float] = None,
        vary: bool = True,
    ) -> FDPBDFitter:
        """Add or configure a parameter to be fitted."""
        if value is None:
            if name == "Dm":
                value = self.mass_config.Dm
            elif name == "dndT_mass":
                value = self.mass_config.dndT_mass
            elif name == "aCET_mass":
                value = self.mass_config.aCET_mass
            elif name in ("dndT_poly", "dndT_PMMA"):
                poly_idx = max(0, self.sample.heating_layer_idx - 1)
                value = float(np.real(self.sample.layers[poly_idx].dndT))
            elif name in ("kth_poly", "kth_PMMA"):
                poly_idx = max(0, self.sample.heating_layer_idx - 1)
                value = self.sample.layers[poly_idx].kth
            else:
                value = 1.0

        if scale is None:
            if name == "Dm":
                scale = 1.0e-12
            elif name in ("dndT_mass", "aCET_mass", "dndT_poly", "dndT_PMMA"):
                scale = 1.0e-4
            elif name in ("thickness_poly", "thickness_PMMA"):
                scale = 1.0e-9
            else:
                scale = abs(value) if abs(value) > 1e-15 else 1.0

        self.params[name] = FitParameter(
            name=name,
            value=value,
            bounds=bounds,
            scale=scale,
            vary=vary,
        )
        return self

    def _apply_params_to_model(self, param_dict: Dict[str, float]) -> Tuple[Sample, MassTransportConfig]:
        poly_idx = max(0, self.sample.heating_layer_idx - 1)
        mass_copy = MassTransportConfig(
            Dm=param_dict.get("Dm", self.mass_config.Dm),
            dndT_mass=param_dict.get("dndT_mass", self.mass_config.dndT_mass),
            aCET_mass=param_dict.get("aCET_mass", self.mass_config.aCET_mass),
        )

        sample_copy = Sample(
            heating_layer_idx=self.sample.heating_layer_idx,
        )
        for i, layer in enumerate(self.sample.layers):
            overrides = layer.overrides.copy()
            if i == poly_idx:
                if "dndT_poly" in param_dict or "dndT_PMMA" in param_dict:
                    val = param_dict.get("dndT_poly", param_dict.get("dndT_PMMA"))
                    overrides["dndT"] = val
                if "kth_poly" in param_dict or "kth_PMMA" in param_dict:
                    val = param_dict.get("kth_poly", param_dict.get("kth_PMMA"))
                    overrides["kth"] = val
                if "aCET_poly" in param_dict:
                    overrides["aCET"] = param_dict["aCET_poly"]
                t_val = param_dict.get("thickness_poly", layer.thickness)
            elif i == self.sample.heating_layer_idx:
                if "dndT_al" in param_dict:
                    overrides["dndT"] = param_dict["dndT_al"]
                t_val = layer.thickness
            else:
                t_val = layer.thickness

            sample_copy.add_layer(
                material=layer.material,
                thickness=t_val,
                name=layer.name,
                **overrides,
            )

        return sample_copy, mass_copy

    def fit(
        self,
        method: str = "least_squares",
        loss: str = "phase",
        max_nfev: int = 200,
        tol: float = 1e-4,
        diff_step: float = 1e-3,
    ) -> FitResult:
        """
        Perform model fit to experimental data.
        """
        active_params = [p for p in self.params.values() if p.vary]
        if not active_params:
            raise ValueError("No active parameters marked with vary=True to fit.")

        x0 = np.array([p.normalized_value for p in active_params], dtype=float)

        lower_bounds = []
        upper_bounds = []
        for p in active_params:
            lb = p.bounds[0] / p.scale if p.bounds[0] is not None else -np.inf
            ub = p.bounds[1] / p.scale if p.bounds[1] is not None else np.inf
            lower_bounds.append(lb)
            upper_bounds.append(ub)

        f_exp = self.exp_data.frequencies
        Vin_exp = self.exp_data.Vin
        Vout_exp = self.exp_data.Vout
        phase_exp = self.exp_data.phase_deg

        def _residual_vector(x_norm: np.ndarray) -> np.ndarray:
            p_dict = {}
            for i, p in enumerate(active_params):
                p_dict[p.name] = p.unnormalize(x_norm[i])

            sample_eval, mass_eval = self._apply_params_to_model(p_dict)
            sim = FDPBDSimulator.simulate(
                sample=sample_eval,
                frequencies=f_exp,
                k_vectors=self.k_vectors,
                beam_config=self.beam_config,
                detector_config=self.detector_config,
                mass_config=mass_eval,
                decompose_components=False,
            )

            if loss == "phase":
                return sim.phase - phase_exp
            elif loss == "complex":
                res_in = sim.Vin - Vin_exp
                res_out = sim.Vout - Vout_exp
                return np.concatenate([res_in, res_out])
            elif loss == "vout":
                return sim.Vout - Vout_exp
            elif loss == "all":
                res_in = sim.Vin - Vin_exp
                res_out = sim.Vout - Vout_exp
                res_phase = (sim.phase - phase_exp) * np.std(Vin_exp) / 10.0
                return np.concatenate([res_in, res_out, res_phase])
            else:
                raise ValueError(f"Unknown loss '{loss}'. Use 'phase', 'complex', 'vout', or 'all'.")

        def _scalar_objective(x_norm: np.ndarray) -> float:
            res = _residual_vector(x_norm)
            return float(np.sum(res**2))

        opt_bounds = [(lb, ub) for lb, ub in zip(lower_bounds, upper_bounds)]

        if method == "least_squares":
            res_opt = least_squares(
                _residual_vector,
                x0=x0,
                bounds=(lower_bounds, upper_bounds),
                diff_step=diff_step,
                max_nfev=max_nfev,
                ftol=tol,
                xtol=tol,
                gtol=tol,
            )
            x_sol = res_opt.x
            cost = res_opt.cost
            success = res_opt.success
            msg = res_opt.message

            errors = {}
            try:
                J = res_opt.jac
                cov = np.linalg.inv(J.T @ J) * (2 * cost / max(1, len(f_exp) - len(x0)))
                for i, p in enumerate(active_params):
                    err_norm = np.sqrt(np.diag(cov))[i]
                    errors[p.name] = float(err_norm * p.scale)
            except Exception:
                for p in active_params:
                    errors[p.name] = None
        else:
            res_opt = minimize(
                _scalar_objective,
                x0=x0,
                method=method,
                bounds=opt_bounds if method in ("L-BFGS-B", "Powell", "Nelder-Mead") else None,
                options={"maxiter": max_nfev, "xatol": tol, "fatol": tol}
                if method == "Nelder-Mead"
                else {"maxiter": max_nfev},
            )
            x_sol = res_opt.x
            cost = res_opt.fun
            success = res_opt.success
            msg = res_opt.message
            errors = {p.name: None for p in active_params}

        final_params = {}
        for p in self.params.values():
            if p.vary:
                idx = active_params.index(p)
                final_params[p.name] = float(p.unnormalize(x_sol[idx]))
            else:
                final_params[p.name] = float(p.value)

        final_sample, final_mass = self._apply_params_to_model(final_params)
        final_sim = FDPBDSimulator.simulate(
            sample=final_sample,
            frequencies=f_exp,
            k_vectors=self.k_vectors,
            beam_config=self.beam_config,
            detector_config=self.detector_config,
            mass_config=final_mass,
            decompose_components=True,
        )

        return FitResult(
            params=final_params,
            errors=errors,
            cost=cost,
            success=success,
            message=msg,
            simulation=final_sim,
        )
