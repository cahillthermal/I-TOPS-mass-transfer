"""
Multilayer bi-directional heat transport solver in spatial and temporal frequency domains.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Tuple, Optional, Union, Dict, List
import numpy as np
from .sample import Sample


@dataclass
class ThermalFieldResult:
    """
    Result of thermal wave transfer matrix calculation.

    Attributes
    ----------
    G_tot : np.ndarray
        Temperature propagator G(k, omega) at the requested evaluation position. Shape (Nk, Nf).
    TempBplus : np.ndarray
        Forward thermal wave amplitude T+ at top surface of target layer. Shape (Nk, Nf).
    TempBminus : np.ndarray
        Backward thermal wave amplitude T- at top surface of target layer. Shape (Nk, Nf).
    TempAplus : np.ndarray
        Forward thermal wave amplitude T+ at bottom surface of layer above target. Shape (Nk, Nf).
    TempAminus : np.ndarray
        Backward thermal wave amplitude T- at bottom surface of layer above target. Shape (Nk, Nf).
    """

    G_tot: np.ndarray
    TempBplus: np.ndarray
    TempBminus: np.ndarray
    TempAplus: np.ndarray
    TempAminus: np.ndarray


@dataclass
class AllLayersThermalResult:
    """
    Complete thermal wave amplitudes for all layers in the sample stack.
    """

    G_heat: np.ndarray
    Gup: np.ndarray
    Gdown: np.ndarray
    A1plus: np.ndarray
    BNminus: np.ndarray
    TempBplus: List[np.ndarray]   # Top surface of each layer i: 0..N-1
    TempBminus: List[np.ndarray]
    TempAplus: List[np.ndarray]   # Bottom surface of each layer i: 0..N-1
    TempAminus: List[np.ndarray]


def solve_all_layer_thermal_fields(
    sample: Sample,
    frequencies: np.ndarray,
    k_vectors: np.ndarray,
    heating_layer_idx: Optional[int] = None,
) -> AllLayersThermalResult:
    """
    Compute thermal wave amplitudes for ALL layers simultaneously in a single optimized pass.
    """
    if heating_layer_idx is None:
        heating_layer_idx = sample.heating_layer_idx

    n_layers = sample.num_layers
    k = np.asarray(k_vectors, dtype=float).reshape(-1, 1)  # (Nk, 1)
    f = np.asarray(frequencies, dtype=float).reshape(1, -1)  # (1, Nf)
    Nk = k.shape[0]
    Nf = f.shape[1]

    kth = sample.kth
    Cv = sample.Cv
    t = sample.thickness
    eta = sample.eta
    D = sample.D

    omega = 2.0 * np.pi * f  # (1, Nf)
    k_term2 = 4.0 * np.pi**2 * (k**2)  # (Nk, 1)
    ii = 1j

    # Compute un, gamman for all layers
    u_all = []
    gamma_all = []
    for i in range(n_layers):
        q2_i = (ii * omega) / D[i]
        u_i = np.sqrt(eta[i] * k_term2 + q2_i)
        u_all.append(u_i)
        gamma_all.append(kth[i] * u_i)

    # 1. Upward iteration from layer 0 up to heating_layer_idx - 1
    # Store Aplus_history, Aminus_history at bottom of each layer
    Aplus_list = [np.ones((Nk, Nf), dtype=complex)]
    Aminus_list = [np.zeros((Nk, Nf), dtype=complex)]

    Aplus = Aplus_list[0]
    Aminus = Aminus_list[0]
    gamman = gamma_all[0]

    for n in range(1, heating_layer_idx):
        un = u_all[n]
        gamman_next = gamma_all[n]
        CC = gamman_next + gamman
        DD = gamman_next - gamman
        temp1 = CC * Aplus + DD * Aminus
        temp2 = DD * Aplus + CC * Aminus
        expterm = np.exp(un * t[n])
        inv_gam = 0.5 / gamman_next
        Aplus = (inv_gam * expterm) * temp1
        Aminus = (inv_gam / expterm) * temp2

        pen_mask = (t[n] * np.abs(un)) > 100.0
        Aplus[pen_mask] = 1.0
        Aminus[pen_mask] = 0.0

        Aplus_list.append(Aplus)
        Aminus_list.append(Aminus)
        gamman = gamman_next

    gamma_up = gamman
    Gup = (Aplus + Aminus) / ((Aminus - Aplus) * gamma_up)

    # 2. Downward iteration from layer N-1 down to heating_layer_idx
    Bplus_list = [None] * n_layers
    Bminus_list = [None] * n_layers

    Bplus = np.zeros((Nk, Nf), dtype=complex)
    Bminus = np.ones((Nk, Nf), dtype=complex)
    Bplus_list[n_layers - 1] = Bplus
    Bminus_list[n_layers - 1] = Bminus
    gamman = gamma_all[n_layers - 1]

    for n in range(n_layers - 2, heating_layer_idx - 1, -1):
        un = u_all[n]
        gamman_prev = gamma_all[n]
        AA = gamman_prev + gamman
        BB = gamman_prev - gamman
        temp1 = AA * Bplus + BB * Bminus
        temp2 = BB * Bplus + AA * Bminus
        expterm = np.exp(un * t[n])
        inv_gam = 0.5 / gamman_prev
        Bplus = (inv_gam / expterm) * temp1
        Bminus = (inv_gam * expterm) * temp2

        pen_mask = (t[n] * np.abs(un)) > 100.0
        Bplus[pen_mask] = 0.0
        Bminus[pen_mask] = 1.0

        Bplus_list[n] = Bplus
        Bminus_list[n] = Bminus
        gamman = gamman_prev

    gamma_down = gamman
    Gdown = (Bplus + Bminus) / ((Bminus - Bplus) * gamma_down)

    G_heat = 1.0 / (-1.0 / Gup + 1.0 / Gdown)

    # 3. Boundary scaling amplitudes
    A1plus = (Gup * Gdown) / ((Gup - Gdown) * (Aplus_list[-1] + Aminus_list[-1]))
    BNminus = (Gup * Gdown) / ((Gup - Gdown) * (Bplus_list[heating_layer_idx] + Bminus_list[heating_layer_idx]))

    # 4. Scale all wave amplitudes
    TempAplus_all = [None] * n_layers
    TempAminus_all = [None] * n_layers
    TempBplus_all = [None] * n_layers
    TempBminus_all = [None] * n_layers

    # For layers above heating
    for idx in range(heating_layer_idx):
        TempAplus_all[idx] = Aplus_list[idx] * A1plus
        TempAminus_all[idx] = Aminus_list[idx] * A1plus

    # For layers at and below heating
    for idx in range(heating_layer_idx, n_layers):
        TempBplus_all[idx] = Bplus_list[idx] * BNminus
        TempBminus_all[idx] = Bminus_list[idx] * BNminus

    # Cross-interface calculations
    # Heating layer top surface
    TempBplus_all[heating_layer_idx] = Bplus_list[heating_layer_idx] * BNminus
    TempBminus_all[heating_layer_idx] = Bminus_list[heating_layer_idx] * BNminus

    # For intermediate layers compute other surface
    for idx in range(1, heating_layer_idx):
        # Convert bottom of idx-1 (TempA) to top of idx (TempB)
        g_prev = gamma_all[idx - 1]
        g_curr = gamma_all[idx]
        AA = g_curr + g_prev
        BB = g_curr - g_prev
        t1 = AA * TempAplus_all[idx - 1] + BB * TempAminus_all[idx - 1]
        t2 = BB * TempAplus_all[idx - 1] + AA * TempAminus_all[idx - 1]
        TempBplus_all[idx] = (0.5 / g_curr) * t1
        TempBminus_all[idx] = (0.5 / g_curr) * t2

    for idx in range(heating_layer_idx, n_layers - 1):
        # Convert top of idx+1 (TempB) to bottom of idx (TempA)
        g_curr = gamma_all[idx]
        g_next = gamma_all[idx + 1]
        AA = g_curr + g_next
        BB = g_curr - g_next
        t1 = AA * TempBplus_all[idx + 1] + BB * TempBminus_all[idx + 1]
        t2 = BB * TempBplus_all[idx + 1] + AA * TempBminus_all[idx + 1]
        TempAplus_all[idx] = (0.5 / g_curr) * t1
        TempAminus_all[idx] = (0.5 / g_curr) * t2

    # Layer 0 has no top surface interface from above; layer N-1 has no bottom interface to below
    TempBplus_all[0] = np.ones((Nk, Nf), dtype=complex)
    TempBminus_all[0] = np.zeros((Nk, Nf), dtype=complex)
    TempAplus_all[n_layers - 1] = np.zeros((Nk, Nf), dtype=complex)
    TempAminus_all[n_layers - 1] = BNminus

    return AllLayersThermalResult(
        G_heat=G_heat,
        Gup=Gup,
        Gdown=Gdown,
        A1plus=A1plus,
        BNminus=BNminus,
        TempBplus=TempBplus_all,
        TempBminus=TempBminus_all,
        TempAplus=TempAplus_all,
        TempAminus=TempAminus_all,
    )


def solve_thermal_field(
    sample: Sample,
    frequencies: np.ndarray,
    k_vectors: np.ndarray,
    heating_layer_idx: Optional[int] = None,
    temp_layer_idx: Optional[int] = None,
    z_offset: float = 0.0,
) -> ThermalFieldResult:
    """
    Compute temperature field at target layer and depth offset.
    """
    if heating_layer_idx is None:
        heating_layer_idx = sample.heating_layer_idx
    if temp_layer_idx is None:
        temp_layer_idx = heating_layer_idx

    all_res = solve_all_layer_thermal_fields(
        sample=sample,
        frequencies=frequencies,
        k_vectors=k_vectors,
        heating_layer_idx=heating_layer_idx,
    )

    k = np.asarray(k_vectors, dtype=float).reshape(-1, 1)
    f = np.asarray(frequencies, dtype=float).reshape(1, -1)
    k_term2 = 4.0 * np.pi**2 * (k**2)
    omega = 2.0 * np.pi * f
    ii = 1j

    TempBplus = all_res.TempBplus[temp_layer_idx]
    TempBminus = all_res.TempBminus[temp_layer_idx]
    TempAplus = all_res.TempAplus[temp_layer_idx - 1] if temp_layer_idx >= 1 else all_res.TempAplus[0]
    TempAminus = all_res.TempAminus[temp_layer_idx - 1] if temp_layer_idx >= 1 else all_res.TempAminus[0]

    if z_offset <= 0:
        prev_idx = max(0, temp_layer_idx - 1)
        q2_prev = (ii * omega) / sample.D[prev_idx]
        u_prev = np.sqrt(sample.eta[prev_idx] * k_term2 + q2_prev)
        G_tot = TempAplus * np.exp(u_prev * z_offset) + TempAminus * np.exp(-u_prev * z_offset)
    else:
        q2_curr = (ii * omega) / sample.D[temp_layer_idx]
        u_curr = np.sqrt(sample.eta[temp_layer_idx] * k_term2 + q2_curr)
        G_tot = TempBplus * np.exp(u_curr * z_offset) + TempBminus * np.exp(-u_curr * z_offset)

    return ThermalFieldResult(
        G_tot=G_tot,
        TempBplus=TempBplus,
        TempBminus=TempBminus,
        TempAplus=TempAplus,
        TempAminus=TempAminus,
    )


def calculate_steady_state_heating(
    sample: Sample,
    k_vectors: np.ndarray,
    pump_power_absorbed: float,
    pump_radius: float,
    heating_layer_idx: Optional[int] = None,
    temp_layer_idx: Optional[int] = None,
    z_offset: float = 0.0,
    f_ss: float = 0.01,
) -> float:
    res = solve_thermal_field(
        sample=sample,
        frequencies=np.array([f_ss]),
        k_vectors=k_vectors,
        heating_layer_idx=heating_layer_idx,
        temp_layer_idx=temp_layer_idx,
        z_offset=z_offset,
    )
    k = np.asarray(k_vectors, dtype=float)
    G = res.G_tot[:, 0]
    kernel = 2.0 * np.pi * pump_power_absorbed * G * np.exp(-np.pi**2 * (k**2) * (pump_radius**2)) * k
    dT_ss = np.real(np.trapezoid(kernel, k))
    return float(dT_ss)


def calculate_temperature_rise(
    sample: Sample,
    frequencies: np.ndarray,
    k_vectors: np.ndarray,
    pump_power_absorbed: float,
    pump_radius: float,
    heating_layer_idx: Optional[int] = None,
    temp_layer_idx: Optional[int] = None,
    z_offset: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray]:
    res = solve_thermal_field(
        sample=sample,
        frequencies=frequencies,
        k_vectors=k_vectors,
        heating_layer_idx=heating_layer_idx,
        temp_layer_idx=temp_layer_idx,
        z_offset=z_offset,
    )
    k = np.asarray(k_vectors, dtype=float)
    k_col = k.reshape(-1, 1)
    kernel = (
        2.0
        * np.pi
        * pump_power_absorbed
        * res.G_tot
        * np.exp(-np.pi**2 * (k_col**2) * (pump_radius**2))
        * k_col
    )
    T_signal = np.trapezoid(kernel, k, axis=0)
    return np.real(T_signal), np.imag(T_signal)
