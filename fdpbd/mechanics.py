"""
Surface displacement calculations (Z0 to Z4) for FD-PBD.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional
import numpy as np
from .sample import Sample
from .config import MassTransportConfig
from .thermal import solve_all_layer_thermal_fields
from .optics import OpticalReflectionResult, solve_multilayer_optics_perturbed


@dataclass
class DisplacementBreakdown:
    """
    Breakdown of the 5 equivalent surface displacement mechanisms.
    """

    Z0: np.ndarray
    Z1: np.ndarray
    Z2: np.ndarray
    Z3: np.ndarray
    Z4: np.ndarray
    Z_tot: np.ndarray
    Z_tot_m: np.ndarray
    Z0_m: np.ndarray
    Z1_m: np.ndarray
    Z2_m: np.ndarray
    Z3_m: np.ndarray
    Z4_m: np.ndarray
    phi0: np.ndarray
    phi1: np.ndarray
    phi2: np.ndarray
    phi3: np.ndarray
    phi4: np.ndarray
    phi_tot: np.ndarray


def solve_displacements(
    sample: Sample,
    frequencies: np.ndarray,
    k_vectors: np.ndarray,
    mass_config: Optional[MassTransportConfig] = None,
    wavelength: float = 783.0e-9,
    optics_res: Optional[OpticalReflectionResult] = None,
    heating_layer_idx: Optional[int] = None,
) -> DisplacementBreakdown:
    """
    Calculate equivalent surface displacements Z0..Z4 across spatial (k) and temporal (f) frequencies.
    """
    if heating_layer_idx is None:
        heating_layer_idx = sample.heating_layer_idx

    k = np.asarray(k_vectors, dtype=float).reshape(-1, 1)  # (Nk, 1)
    f = np.asarray(frequencies, dtype=float).reshape(1, -1)  # (1, Nf)
    Nk = k.shape[0]
    Nf = f.shape[1]

    n_layers = sample.num_layers
    kth = sample.kth
    Cv = sample.Cv
    t = sample.thickness
    eta = sample.eta
    D = sample.D
    n_refrac = sample.n_refrac
    dndT = sample.dndT
    aCET = sample.aCET
    nu = sample.nu
    Y = sample.Y
    Vl = sample.Vl
    Vt = sample.Vt

    n1_real = float(np.real(n_refrac[0])) if np.real(n_refrac[0]) > 0 else 1.0
    omega = 2.0 * np.pi * f
    k_term2 = 4.0 * np.pi**2 * (k**2)
    ii = 1j

    poly_idx = heating_layer_idx - 1
    if mass_config is not None:
        Dm = mass_config.Dm
        aCET_mass = mass_config.aCET_mass
    elif poly_idx >= 0:
        Dm = sample.layers[poly_idx].Dm
        aCET_mass = sample.layers[poly_idx].aCET_mass
    else:
        Dm = 1.0e-12
        aCET_mass = 0.0

    # 1. Single-pass fast thermal wave solver for all layers
    all_thermal = solve_all_layer_thermal_fields(
        sample=sample,
        frequencies=frequencies,
        k_vectors=k_vectors,
        heating_layer_idx=heating_layer_idx,
    )
    G_nheat = all_thermal.G_heat

    # ------------------ (0) Fresnel Reflection (Z0) ------------------
    if optics_res is None:
        optics_res = solve_multilayer_optics_perturbed(
            sample=sample,
            frequencies=frequencies,
            k_vectors=k_vectors,
            G_nheat=G_nheat,
            mass_config=mass_config,
            wavelength=wavelength,
            heating_layer_idx=heating_layer_idx,
        )

    dphi0 = -optics_res.dphi
    deta_Z0 = -optics_res.dphi * wavelength / (4.0 * np.pi * n1_real)
    deta_Z0m = deta_Z0 * n1_real

    # ------------------ (1) Thermal expansion of thin films (Z1) ------------------
    deta_Z1 = np.zeros((Nk, Nf), dtype=complex)
    deta_phi1 = np.zeros((Nk, Nf), dtype=complex)

    if heating_layer_idx < (n_layers - 1):
        for i in range(heating_layer_idx, n_layers - 1):
            q2_i = (ii * omega) / D[i]
            un_i = np.sqrt(eta[i] * k_term2 + q2_i)
            TempBplus_i = all_thermal.TempBplus[i]
            TempBminus_i = all_thermal.TempBminus[i]

            poisson_factor = (1.0 + nu[i]) / (1.0 - nu[i]) if abs(1.0 - nu[i]) > 1e-6 else 1.0
            integral_term = (
                (np.exp(un_i * t[i]) - 1.0) * TempBplus_i
                - (np.exp(-un_i * t[i]) - 1.0) * TempBminus_i
            )
            dZ1_i = poisson_factor * aCET[i] / un_i * integral_term
            deta_Z1 += dZ1_i
            deta_phi1 += (4.0 * np.pi * n1_real / wavelength) * dZ1_i

    deta_Z1m = deta_Z1 * n1_real

    # ------------------ (2) Substrate deformation due to stresses (Z2) ------------------
    gT = np.zeros((Nk, Nf), dtype=complex)
    for i in range(1, n_layers - 1):
        q2_i = (ii * omega) / D[i]
        un_i = np.sqrt(eta[i] * k_term2 + q2_i)
        TempBplus_i = all_thermal.TempBplus[i]
        TempBminus_i = all_thermal.TempBminus[i]

        poisson_factor = Y[i] / (1.0 - nu[i]) if abs(1.0 - nu[i]) > 1e-6 else 0.0
        integral_term = (
            (np.exp(un_i * t[i]) - 1.0) * TempBplus_i
            - (np.exp(-un_i * t[i]) - 1.0) * TempBminus_i
        )
        gT += poisson_factor * aCET[i] / un_i * integral_term

    sub_idx = n_layers - 1
    if poly_idx >= 0 and Dm > 0 and abs(aCET_mass) > 1e-15:
        q2mass = (ii * omega) / Dm
        unmass = np.sqrt(eta[0] * k_term2 + q2mass)
        l_poly = t[poly_idx]
        poisson_poly = Y[poly_idx] / (1.0 - nu[poly_idx]) if abs(1.0 - nu[poly_idx]) > 1e-6 else 0.0
        gTM = gT + poisson_poly * aCET_mass * G_nheat / unmass * np.sinh(unmass * l_poly) / np.cosh(unmass * l_poly)
    else:
        gTM = gT

    W2 = omega**2  # (1, Nf)
    Vl_sub = Vl[sub_idx] if Vl[sub_idx] > 0 else 5970.0
    Vt_sub = Vt[sub_idx] if Vt[sub_idx] > 0 else 3770.0

    apha = np.sqrt(eta[sub_idx] * k_term2 - (W2 / (Vl_sub**2)) + 0j)
    beta = np.sqrt(eta[sub_idx] * k_term2 - (W2 / (Vt_sub**2)) + 0j)

    with np.errstate(divide="ignore", invalid="ignore"):
        num_z2 = k_term2 * (beta**2 + k_term2 - 2.0 * apha * beta)
        den_z2 = (beta**2 + k_term2) ** 2 - 4.0 * apha * beta * k_term2
        deta_Z2_part = np.where(np.abs(den_z2) > 1e-30, num_z2 / den_z2, 0.0)

    # Low frequency approximation
    low_freq_mask = (W2 * 10000.0) < (Vt_sub**2 * k_term2)
    deta_Z2_part[low_freq_mask] = 0.5 * (Vt_sub**2) / (Vt_sub**2 - Vl_sub**2)

    Y_sub = Y[sub_idx] if Y[sub_idx] > 0 else 73.0e9
    nu_sub = nu[sub_idx] if nu[sub_idx] > 0 else 0.17
    deta_Z2 = -deta_Z2_part * 2.0 * (1.0 + nu_sub) / Y_sub * gTM
    deta_Z2m = deta_Z2 * n1_real
    deta_phi2 = deta_Z2 * (4.0 * np.pi * n1_real) / wavelength

    # ------------------ (3) Thermal expansion of the substrate (Z3) ------------------
    q2_sub = (ii * omega) / D[sub_idx]
    un_sub = np.sqrt(eta[sub_idx] * k_term2 + q2_sub)
    TempBminusN = all_thermal.TempBminus[sub_idx]

    xi = np.sqrt(apha**2 + q2_sub)
    with np.errstate(divide="ignore", invalid="ignore"):
        num_z3 = beta**4 - 16.0 * (np.pi**4) * (k**4)
        den_z3 = (beta**2 + 4.0 * np.pi**2 * (k**2)) ** 2 - 16.0 * apha * beta * (np.pi**2) * (k**2)
        deta_Z3_part = np.where(np.abs(den_z3) > 1e-30, num_z3 / den_z3, 0.0)

    deta_Z3_part[low_freq_mask] = (Vl_sub**2) / (Vl_sub**2 - Vt_sub**2)

    poisson_sub_factor = (1.0 + nu_sub) / (1.0 - nu_sub)
    deta_Z3 = (
        deta_Z3_part
        * poisson_sub_factor
        * aCET[sub_idx]
        * (un_sub / q2_sub)
        * (1.0 - apha / xi)
        * TempBminusN
    )
    deta_Z3m = deta_Z3 * n1_real
    deta_phi3 = deta_Z3 * (4.0 * np.pi * n1_real) / wavelength

    # ------------------ (4) Thermo-optic effect in ambient/vapor (Z4) ------------------
    q2_amb = (ii * omega) / D[0]
    un_amb = np.sqrt(eta[0] * k_term2 + q2_amb)
    TempAplus2 = all_thermal.TempAplus[0]

    dndT_amb = dndT[0] if not np.isnan(dndT[0]) else 0.0
    deta_phi4 = -4.0 * np.pi / wavelength * dndT_amb * TempAplus2 / un_amb
    deta_Z4 = deta_phi4 * wavelength / (4.0 * np.pi * n1_real)
    deta_Z4m = deta_Z4 * n1_real

    # ------------------ Total Values ------------------
    deta_Ztot = deta_Z0 + deta_Z1 + deta_Z2 + deta_Z3 + deta_Z4
    deta_Ztotm = deta_Z0m + deta_Z1m + deta_Z2m + deta_Z3m + deta_Z4m
    deta_phitot = dphi0 + deta_phi1 + deta_phi2 + deta_phi3 + deta_phi4

    return DisplacementBreakdown(
        Z0=deta_Z0,
        Z1=deta_Z1,
        Z2=deta_Z2,
        Z3=deta_Z3,
        Z4=deta_Z4,
        Z_tot=deta_Ztot,
        Z_tot_m=deta_Ztotm,
        Z0_m=deta_Z0m,
        Z1_m=deta_Z1m,
        Z2_m=deta_Z2m,
        Z3_m=deta_Z3m,
        Z4_m=deta_Z4m,
        phi0=dphi0,
        phi1=deta_phi1,
        phi2=deta_phi2,
        phi3=deta_phi3,
        phi4=deta_phi4,
        phi_tot=deta_phitot,
    )
