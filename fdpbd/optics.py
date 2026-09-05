"""
High-performance multilayer optical reflection and perturbed transfer matrix solver.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Tuple, Optional
import numpy as np
from .sample import Sample
from .config import MassTransportConfig


@dataclass
class OpticalReflectionResult:
    """
    Results of multilayer optical reflection calculation.

    Attributes
    ----------
    reflectance : float
        Nominal unperturbed power reflectance (0 to 1).
    phase : float
        Nominal unperturbed Fresnel reflection phase (radians).
    dphi : np.ndarray
        Complex differential phase dphi = dphi_real + 1j * dphi_imag of shape (Nk, Nf).
    dR : np.ndarray
        Complex differential reflectance dR = dR_real + 1j * dR_imag of shape (Nk, Nf).
    """

    reflectance: float
    phase: float
    dphi: np.ndarray
    dR: np.ndarray


def solve_fresnel_reflection(
    sample: Sample,
    theta_deg: float = 0.0,
    wavelength: float = 783.0e-9,
) -> Tuple[float, float, float]:
    """
    Calculate optical reflectance and transmittance for unperturbed multilayer stack.
    """
    c0 = 3.0e8
    omega_L = 2.0 * np.pi * c0 / wavelength
    n_refrac = np.asarray(sample.n_refrac, dtype=complex)
    thickness = np.asarray(sample.thickness, dtype=float)
    n_layers = len(sample)

    theta = np.zeros(n_layers, dtype=complex)
    theta[0] = np.radians(theta_deg)
    for i in range(1, n_layers):
        theta[i] = np.arcsin(np.sin(theta[0]) * n_refrac[0] / n_refrac[i])

    n_sub_plus1 = 1.0
    theta_sub_plus1 = np.arcsin(np.sin(theta[0]) * n_refrac[0] / n_sub_plus1)
    gammam = n_sub_plus1 / np.cos(theta_sub_plus1)

    Fplus = 1.0 + 0.0j
    Fminus = 0.0 + 0.0j

    for m in range(n_layers, 1, -1):
        idx = m - 1
        vmminus = 1j * omega_L * n_refrac[idx] * np.cos(theta[idx]) / c0
        gammamminus = n_refrac[idx] / np.cos(theta[idx])

        AA = gammamminus + gammam
        BB = gammamminus - gammam
        temp1 = AA * Fplus + BB * Fminus
        temp2 = BB * Fplus + AA * Fminus
        expterm = np.exp(vmminus * thickness[idx])

        Fplus = (0.5 / (gammamminus * expterm)) * temp1
        Fminus = (0.5 / gammamminus) * expterm * temp2

        im_n = np.imag(n_refrac[idx])
        if abs(im_n) > 1e-12:
            L_pen = wavelength / (2.0 * np.pi * im_n * np.real(np.cos(theta[idx])))
            if abs(L_pen) > 0 and (thickness[idx] / abs(L_pen)) > 10.0:
                Fplus = 1.0 + 0.0j
                Fminus = 0.0 + 0.0j

        gammam = gammamminus

    vmminus = 1j * omega_L * n_refrac[0] * np.cos(theta[0]) / c0
    gammamminus = n_refrac[0] / np.cos(theta[0])
    AA = gammamminus + gammam
    BB = gammamminus - gammam
    temp1 = AA * Fplus + BB * Fminus
    temp2 = BB * Fplus + AA * Fminus
    Cplus = (0.5 / gammamminus) * temp1
    Cminus = (0.5 / gammamminus) * temp2

    r = Cminus / Cplus
    R = float(np.abs(r) ** 2)
    T = float(np.abs(1.0 / Cplus) ** 2)
    phi = float(np.angle(r))

    return R, T, phi


def solve_multilayer_optics_perturbed(
    sample: Sample,
    frequencies: np.ndarray,
    k_vectors: np.ndarray,
    G_nheat: np.ndarray,
    mass_config: Optional[MassTransportConfig] = None,
    wavelength: float = 783.0e-9,
    heating_layer_idx: Optional[int] = None,
    max_sublayers: int = 50,
) -> OpticalReflectionResult:
    """
    Optimized vectorized computation of differential Fresnel reflection dphi and dR.
    """
    if heating_layer_idx is None:
        heating_layer_idx = sample.heating_layer_idx

    k = np.asarray(k_vectors, dtype=float).reshape(-1, 1)  # (Nk, 1)
    f = np.asarray(frequencies, dtype=float).reshape(1, -1)  # (1, Nf)
    Nk = k.shape[0]
    Nf = f.shape[1]

    c0 = 3.0e8
    omega_L = 2.0 * np.pi * c0 / wavelength
    omega_L_c0 = omega_L / c0
    omega_H = 2.0 * np.pi * f  # (1, Nf)
    k_term2 = 4.0 * np.pi**2 * (k**2)  # (Nk, 1)

    t_orig = sample.thickness
    n_orig = sample.n_refrac
    v_orig = sample.nu
    aCET_orig = sample.aCET
    dndT_orig = sample.dndT

    poly_idx = heating_layer_idx - 1
    if mass_config is not None:
        Dm = mass_config.Dm
        dndT_mass = mass_config.dndT_mass
        aCET_mass = mass_config.aCET_mass
    elif poly_idx >= 0:
        Dm = sample.layers[poly_idx].Dm
        dndT_mass = sample.layers[poly_idx].dndT_mass
        aCET_mass = sample.layers[poly_idx].aCET_mass
    else:
        Dm = 1.0e-12
        dndT_mass = 0.0
        aCET_mass = 0.0

    # Polymer sublayer discretization
    if heating_layer_idx == 1 or poly_idx < 0 or Dm <= 0:
        N_sub = 1
        t_new = [t_orig[i] for i in range(len(t_orig))]
        n_new = [n_orig[i] for i in range(len(n_orig))]
        v_new = [v_orig[i] for i in range(len(v_orig))]
        aCET_new = [aCET_orig[i] for i in range(len(aCET_orig))]
        dndT_new = [dndT_orig[i] for i in range(len(dndT_orig))]
        poly_sub_range = range(0, 0)
    else:
        l_poly = t_orig[poly_idx]
        omega_max = 2.0 * np.pi * np.max(frequencies)
        l_diff = np.sqrt(Dm / omega_max) if omega_max > 0 else l_poly
        raw_N_sub = int(np.round(l_poly / (l_diff / 5.0))) + 1
        N_sub = max(1, min(raw_N_sub, max_sublayers))

        N_total = len(t_orig) - 1 + N_sub
        t_new = [0.0] * N_total
        n_new = [0.0 + 0.0j] * N_total
        v_new = [0.0] * N_total
        aCET_new = [0.0] * N_total
        dndT_new = [0.0 + 0.0j] * N_total

        poly_start = poly_idx
        poly_end = poly_idx + N_sub
        poly_sub_range = range(poly_start, poly_end)

        for i in range(poly_start):
            t_new[i] = t_orig[i]
            n_new[i] = n_orig[i]
            v_new[i] = v_orig[i]
            aCET_new[i] = aCET_orig[i]
            dndT_new[i] = dndT_orig[i]

        for i in poly_sub_range:
            t_new[i] = l_poly / N_sub
            n_new[i] = n_orig[poly_idx]
            v_new[i] = v_orig[poly_idx]
            aCET_new[i] = aCET_orig[poly_idx]
            dndT_new[i] = dndT_orig[poly_idx]

        for i in range(poly_end, N_total):
            orig_i = i - N_sub + 1
            t_new[i] = t_orig[orig_i]
            n_new[i] = n_orig[orig_i]
            v_new[i] = v_orig[orig_i]
            aCET_new[i] = aCET_orig[orig_i]
            dndT_new[i] = dndT_orig[orig_i]

    # ------------------ Base Unperturbed Reflection ------------------
    al_idx = (heating_layer_idx - 1 + N_sub) if heating_layer_idx > 1 else heating_layer_idx
    q2_al = (n_new[al_idx] * omega_L_c0) ** 2
    vm_al = np.sqrt(k_term2 - q2_al + 0j)
    gamma_curr = vm_al

    Fplus = np.ones((Nk, Nf), dtype=complex)
    Fminus = np.zeros((Nk, Nf), dtype=complex)

    for idx in range(al_idx - 1, 0, -1):
        q2_m = (n_new[idx] * omega_L_c0) ** 2
        vm_minus = np.sqrt(k_term2 - q2_m + 0j)
        gamma_minus = vm_minus

        AA = gamma_minus + gamma_curr
        BB = gamma_minus - gamma_curr
        temp1 = AA * Fplus + BB * Fminus
        temp2 = BB * Fplus + AA * Fminus
        expterm = np.exp(vm_minus * t_new[idx])

        inv_gam = 0.5 / gamma_minus
        Fplus = (inv_gam / expterm) * temp1
        Fminus = inv_gam * expterm * temp2
        gamma_curr = gamma_minus

    q2_0 = (n_new[0] * omega_L_c0) ** 2
    vm_minus = np.sqrt(k_term2 - q2_0 + 0j)
    gamma_minus = vm_minus
    AA = gamma_minus + gamma_curr
    BB = gamma_minus - gamma_curr
    temp1 = AA * Fplus + BB * Fminus
    temp2 = BB * Fplus + AA * Fminus
    inv_gam0 = 0.5 / gamma_minus
    Cplus = inv_gam0 * temp1
    Cminus = inv_gam0 * temp2

    r_base = Cminus / Cplus
    R_base = np.abs(r_base) ** 2
    phi_base = np.angle(r_base)

    # ------------------ Perturbation Precomputations ------------------
    qmass2 = (1j * omega_H) / Dm if Dm > 0 else 0j
    lambda_mass = np.sqrt(k_term2 + qmass2 + 0j)
    l_poly = t_orig[poly_idx] if poly_idx >= 0 else 0.0
    cosh_l = np.cosh(lambda_mass * l_poly) if l_poly > 0 else 1.0
    inv_cosh_l = 1.0 / cosh_l if l_poly > 0 else 1.0

    cosh_factors = []
    if len(poly_sub_range) > 0:
        for i in poly_sub_range:
            sub_i = i - poly_start
            Z = (N_sub - sub_i - 0.5) * t_new[i]
            cosh_factors.append(np.cosh(lambda_mass * Z) * inv_cosh_l)

    def _eval_perturbed_optics(deltaT: np.ndarray, is_real: bool) -> Tuple[np.ndarray, np.ndarray]:
        dndT_al = dndT_new[al_idx] if not np.isnan(dndT_new[al_idx]) else 0.0
        n_pert_al = n_new[al_idx] + dndT_al * deltaT if abs(dndT_al) > 0 else n_new[al_idx]
        q2_al = (n_pert_al * omega_L_c0) ** 2
        vm_al = np.sqrt(k_term2 - q2_al + 0j)
        gamma_curr = vm_al

        Fp = np.ones((Nk, Nf), dtype=complex)
        Fm = np.zeros((Nk, Nf), dtype=complex)

        for idx in range(al_idx - 1, 0, -1):
            dndT_val = dndT_new[idx] if not np.isnan(dndT_new[idx]) else 0.0

            if idx in poly_sub_range:
                sub_i = idx - poly_start
                cf = cosh_factors[sub_i]
                deltaC = np.real(G_nheat * cf) if is_real else np.imag(G_nheat * cf)
                nu_val = v_new[idx]
                poisson_fac = (1.0 + nu_val) / (1.0 - nu_val) if abs(1.0 - nu_val) > 1e-6 else 1.0
                t_pert = t_new[idx] * (1.0 + aCET_new[idx] * poisson_fac * deltaT + aCET_mass * poisson_fac * deltaC)
                n_pert = n_new[idx] + dndT_val * deltaT + dndT_mass * deltaC
            else:
                t_pert = t_new[idx]
                n_pert = n_new[idx] + dndT_val * deltaT if abs(dndT_val) > 0 else n_new[idx]

            q2_m = (n_pert * omega_L_c0) ** 2
            vm_minus = np.sqrt(k_term2 - q2_m + 0j)
            gamma_minus = vm_minus

            AA = gamma_minus + gamma_curr
            BB = gamma_minus - gamma_curr
            temp1 = AA * Fp + BB * Fm
            temp2 = BB * Fp + AA * Fm
            expterm = np.exp(vm_minus * t_pert)

            inv_gam = 0.5 / gamma_minus
            Fp = (inv_gam / expterm) * temp1
            Fm = inv_gam * expterm * temp2
            gamma_curr = gamma_minus

        # Ambient layer 0
        dndT_0 = dndT_new[0] if not np.isnan(dndT_new[0]) else 0.0
        n_pert_0 = n_new[0] + dndT_0 * deltaT if abs(dndT_0) > 0 else n_new[0]
        q2_0 = (n_pert_0 * omega_L_c0) ** 2
        vm_minus = np.sqrt(k_term2 - q2_0 + 0j)
        gamma_minus = vm_minus
        AA = gamma_minus + gamma_curr
        BB = gamma_minus - gamma_curr
        temp1 = AA * Fp + BB * Fm
        temp2 = BB * Fp + AA * Fm
        inv_gam0 = 0.5 / gamma_minus
        Cp = inv_gam0 * temp1
        Cm = inv_gam0 * temp2

        r_pert = Cm / Cp
        R_pert = np.abs(r_pert) ** 2
        phi_pert = np.angle(r_pert)

        return (phi_pert - phi_base), (R_pert - R_base)

    dphi_real, dR_real = _eval_perturbed_optics(np.real(G_nheat), is_real=True)
    dphi_imag, dR_imag = _eval_perturbed_optics(np.imag(G_nheat), is_real=False)

    dphi_complex = dphi_real + 1j * dphi_imag
    dR_complex = dR_real + 1j * dR_imag

    return OpticalReflectionResult(
        reflectance=float(np.mean(R_base)),
        phase=float(np.mean(phi_base)),
        dphi=dphi_complex,
        dR=dR_complex,
    )
