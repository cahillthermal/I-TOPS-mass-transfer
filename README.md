# Frequency-Domain Probe Beam Deflection (FD-PBD)

A high-performance Python package for simulating and analyzing Frequency-Domain Probe Beam Deflection (FD-PBD) measurements of mass diffusion, thermal transport, and thermomechanical deformations in nanoscale thin films and multilayer stacks.

The physical formulation and experimental methodology are described in:
> **X. Xie, J. M. Dennison, J. Shin, Z. Diao, and D. G. Cahill**, *Measurement of water vapor diffusion in nanoscale polymer films by frequency-domain probe beam deflection*, **Rev. Sci. Instrum. 89, 104904 (2018)**. [doi:10.1063/1.5039731](https://doi.org/10.1063/1.5039731)

---

## Overview and Physical Formalism

Frequency-domain probe beam deflection (FD-PBD) measures the mutual-diffusion of vapor in nanoscale polymer films with microsecond temporal resolution. In FD-PBD, an intensity-modulated pump laser beam heats an absorbing transducer layer. The resulting periodic temperature excursion alters the equilibrium concentration of vapor at the film surface, driving molecular diffusion through the polymer layer.

A displaced probe laser beam detects the lateral gradients in the optical phase and surface displacements. Beam deflection originates from five distinct physical mechanisms:
1. **$Z_0$ (Multilayer Fresnel Optics)**: Complex optical phase change $\Delta\psi$ created by thermo-optic variations ($\mathrm{d}n/\mathrm{d}T$), thermal expansion ($\alpha_T$), and water concentration perturbations ($\mathrm{d}n/\mathrm{d}C$, $\alpha_m$).
2. **$Z_1$ (Transducer Thermal Expansion)**: Effective vertical thermal expansion of the metal transducer and underlying thin films.
3. **$Z_2$ (Substrate Stress-Induced Deformation)**: Elastic deformation of the substrate produced by lateral shear stresses in the heated and swelled thin films above.
4. **$Z_3$ (Substrate Thermal Expansion)**: Vertical thermal expansion of the bulk substrate integrated across the thermal penetration depth.
5. **$Z_4$ (Ambient Thermo-Optic Effect)**: Refractive index gradient formed in the ambient vapor or gas layer above the sample.

```
                              [ Incident Probe Beam ]
                                      \     ^
                                       \   /  [ Reflected Deflected Beam ]
                                        v /
=================================================================================
 Layer 0: Ambient Gas / Vapor (Water Vapor / Air)                   [Z4: dn/dT]
---------------------------------------------------------------------------------
 Layer 1: Polymer Thin Film (PMMA, PI, PVP, CA, PIP/TMC)            [Z0, Z2: Dm, dn/dC, am]
---------------------------------------------------------------------------------
 Layer 2: Metal Transducer (Al, 80 nm, Heat Absorption Plane)       [Z1, Z2: aT, Y, nu]
---------------------------------------------------------------------------------
 Layer 3: Bulk Substrate (Fused Silica)                             [Z2, Z3: aT, Y, nu, Vl, Vt]
=================================================================================
```

---

## Mathematical Formulation

### 1. Thermal Transport
The sample and ambient vapor are modeled as an $N$-layer structure. Heat diffusion in cylindrical coordinates is transformed to algebraic equations via Hankel ($k$) and Fourier ($\omega = 2\pi f$) transforms:
$$\frac{\partial^2 T_j}{\partial z^2} - u_j^2 T_j = 0, \quad u_j = \sqrt{4\pi^2 \eta_j k^2 + \frac{i\omega}{D_j}}$$
where $D_j = \Lambda_j / C_{v,j}$ is the thermal diffusivity, $\Lambda_j$ is the thermal conductivity, and $C_{v,j}$ is the volumetric heat capacity.

Bi-directional thermal transfer matrices yield the upward ($G_{\mathrm{up}}$) and downward ($G_{\mathrm{down}}$) thermal propagators. The combined thermal response at the heating interface is:
$$G_{\mathrm{tot}}(k, \omega) = \left( \frac{1}{G_{\mathrm{down}}} - \frac{1}{G_{\mathrm{up}}} \right)^{-1}$$

### 2. Mass Transport
Equilibrium concentration boundary conditions are imposed at the polymer surface ($z = 0$). The mass diffusion equation in the polymer film of thickness $L_p$ yields:
$$C(k, \omega, z) = \frac{\mathrm{d}C_{\mathrm{eq}}}{\mathrm{d}T} T_p \frac{\cosh(\lambda_m (L_p - z))}{\cosh(\lambda_m L_p)}, \quad \lambda_m = \sqrt{4\pi^2 k^2 + \frac{i\omega}{D_m}}$$
where $D_m$ is the mutual mass diffusivity and $T_p$ is the temperature excursion.

### 3. Surface Displacements and Beam Deflection
The equivalent surface displacement is obtained by summing all five contributions:
$$Z_{\mathrm{tot}, m}(k, \omega) = \sum_{i=0}^{4} Z_{i, m}(k, \omega)$$

The average deflection angle $\theta(\omega)$ is calculated by Hankel convolution with the pump ($P$) and probe ($S$) spatial distributions:
$$\theta(\omega) = 8\pi^2 \int_0^\infty P(k) Z_{\mathrm{tot}, m}(k, \omega) S(k) J_1(2\pi k r_0) k^2 \, \mathrm{d}k$$
where $r_0$ is the pump-probe beam separation and $J_1$ is the first-order Bessel function.

The deflection angle is converted to root-mean-square lock-in voltage:
$$V_{\mathrm{lock-in}} = \frac{F \theta}{\sqrt{2} w_0} \frac{G_2}{0.65} V_{\mathrm{sum}}$$
where $F$ is the objective focal length, $w_0$ is the beam radius at the detector, $G_2$ is the electronic gain, and $V_{\mathrm{sum}}$ is the detector sum voltage.

---

## Package Architecture

The implementation is partitioned into modular components:

```
fdpbd/
├── __init__.py         # Package interface and top-level exports
├── materials.py        # Material definitions and thermophysical property database
├── sample.py           # Multilayer stack and Layer data structures
├── config.py           # Optical beam, detector, and simulation configurations
├── thermal.py          # Vectorized bi-directional thermal wave matrix solver
├── optics.py           # Multilayer optical reflection and perturbed transfer matrix
├── mechanics.py        # Surface displacement solvers (Z0 through Z4)
├── deflection.py       # Hankel quadrature and lock-in voltage conversion
├── simulator.py        # High-level simulation engine (FDPBDSimulator)
├── calibration.py      # Lock-in data loading, phase rotation, and amplitude calibration
├── fitting.py          # Non-linear parameter estimation (FDPBDFitter)
└── plotting.py         # Quadrature and displacement breakdown plotting routines
```

### Performance Enhancements
- **Vectorized 2D $(k, \omega)$ Grid Operations**: Inner numerical loops over spatial frequency $k$ and temporal frequency $\omega$ are vectorized across NumPy arrays.
- **Single-Pass Thermal Matrix Solver**: Forward and backward thermal wave amplitudes for all layers are evaluated in a single sweep, avoiding redundant matrix inversions.
- **Optimized Optical Layer Propagation**: Matrix multiplication begins directly at the optically opaque transducer layer, eliminating unneeded calculations in the substrate.
- **Execution Speed**: A 100-frequency simulation executes in $\approx 400\text{ ms}$, representing a $>50\times$ acceleration compared to unvectorized scripts.

---

## Usage Examples

### 1. Forward Simulation

```python
import numpy as np
from fdpbd import (
    Sample, PMMA, ALUMINUM, FUSED_SILICA, WATER_VAPOR,
    LaserBeamConfig, DetectorConfig, MassTransportConfig,
    FDPBDSimulator, plot_all
)

# 1. Define sample stack
sample = Sample.create_standard_fdpbd_sample(
    polymer_thickness=114e-9,
    al_thickness=80e-9,
    polymer_material=PMMA,
    ambient_material=WATER_VAPOR,
    transducer_material=ALUMINUM,
    substrate_material=FUSED_SILICA,
)

# 2. Configure optical parameters for 5X objective lens
beam = LaserBeamConfig(
    pump_power=3.0e-3,       # 3.0 mW pump power
    pump_radius=11.3e-6,     # 11.3 um 1/e^2 radius
    beam_offset=12.0e-6,     # 12.0 um separation
    wavelength=783e-9,       # 783 nm probe wavelength
    focal_length=40.0e-3,    # 40 mm focal length
)

# 3. Set mass transport properties
mass = MassTransportConfig(
    Dm=1.1e-12,              # Diffusivity (m^2/s)
    dndT_mass=-1.0e-4,       # (dn/dC)*(dC/dT) (1/K)
    aCET_mass=-2.0e-5,       # (dL/L/dC)*(dC/dT) (1/K)
)

# 4. Execute simulation across 10 Hz - 100 kHz
frequencies = np.logspace(np.log10(10.0), np.log10(100.0e3), 100)
result = FDPBDSimulator.simulate(
    sample=sample,
    frequencies=frequencies,
    beam_config=beam,
    mass_config=mass,
    decompose_components=True,
)

print(f"Reflectance: {result.reflectance:.4f}")
print(f"Steady-State Temperature Rise: {result.steady_state_dT:.2f} K")
print(f"In-phase Signal (10 Hz): {result.Vin[0]:.2f} uV")

# 5. Plot 4-panel response
fig, axes = plot_all(result)
```

### 2. Experimental Data Calibration & Parameter Fitting

```python
from fdpbd import (
    Sample, PMMA, LaserBeamConfig, DetectorConfig,
    ExperimentData, CalibrationData, FDPBDFitter, plot_fit
)

# Load raw experimental data and detector reference calibration
raw_exp = ExperimentData.load_from_file("exp_data.txt", fix_legacy_labview=True)
cali = CalibrationData.load_from_file("cali_data.txt")

# Apply phase rotation and amplitude correction
exp_data = raw_exp.apply_calibration(cali).filter_frequency_range(10.0, 20000.0)

# Configure sample and optimization bounds
sample = Sample.create_standard_fdpbd_sample(polymer_thickness=114e-9, al_thickness=78e-9)
fitter = FDPBDFitter(sample=sample, exp_data=exp_data)
fitter.add_parameter("Dm", value=0.5e-12, bounds=(1e-14, 1e-10))
fitter.add_parameter("dndT_mass", value=-0.8e-4, bounds=(-5e-4, 0.0))

# Execute least-squares fit
fit_result = fitter.fit(method="least_squares", loss="phase")

print(f"Fitted Dm: {fit_result.params['Dm']:.3e} m^2/s")
fig, axes = plot_fit(fit_result, exp_data)
```

---

## Test Verification

Unit tests and benchmark profiles are provided in the `tests/` directory:

```bash
# Run physical simulation test suite
python tests/test_simulation.py

# Run parameter fitting recovery test
python tests/test_fitting.py

# Run computational benchmark profile
python tests/benchmark_profile.py
```

---

## License

This software is released under the **GNU General Public License v3.0 (GPL-3.0)**. See the [`LICENSE`](LICENSE) file for details.
