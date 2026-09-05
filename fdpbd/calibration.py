"""
Experimental data ingestion, phase/amplitude calibration, and frequency filtering.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional, Union, Tuple
import numpy as np
from scipy.interpolate import interp1d


@dataclass
class CalibrationData:
    """
    Reference system calibration data (photodetector & optical path frequency response).

    Attributes
    ----------
    frequencies : np.ndarray
        Calibration modulation frequencies in Hz.
    Vin : np.ndarray
        In-phase calibration voltage (uV).
    Vout : np.ndarray
        Out-of-phase calibration voltage (uV).
    """

    frequencies: np.ndarray
    Vin: np.ndarray
    Vout: np.ndarray

    @property
    def complex_signal(self) -> np.ndarray:
        return self.Vin + 1j * self.Vout

    @property
    def amplitude(self) -> np.ndarray:
        return np.abs(self.complex_signal)

    @property
    def phase_deg(self) -> np.ndarray:
        return np.angle(self.complex_signal, deg=True)

    @classmethod
    def load_from_file(cls, filepath: str) -> CalibrationData:
        """
        Load calibration data from 7-column whitespace/tab delimited text file.
        Columns: [f, t, Vin, Vout, ..., V0, ...]
        """
        data = np.loadtxt(filepath)
        if data.ndim == 1:
            data = data.reshape(1, -1)
        # If transposed (7 rows, N cols), adjust to (N cols, 7 rows)
        if data.shape[0] == 7 and data.shape[1] > 7:
            data = data.T

        frequencies = data[:, 0]
        Vin = data[:, 2]
        Vout = data[:, 3]
        return cls(frequencies=frequencies, Vin=Vin, Vout=Vout)


@dataclass
class ExperimentData:
    """
    Experimental lock-in measurement data for FD-PBD.

    Attributes
    ----------
    frequencies : np.ndarray
        Modulation frequencies in Hz, shape (N,).
    Vin : np.ndarray
        In-phase lock-in signal in uV, shape (N,).
    Vout : np.ndarray
        Out-of-phase lock-in signal in uV, shape (N,).
    V0 : Optional[np.ndarray]
        Measured DC / sum voltage (V), optional.
    """

    frequencies: np.ndarray
    Vin: np.ndarray
    Vout: np.ndarray
    V0: Optional[np.ndarray] = None

    @property
    def complex_signal(self) -> np.ndarray:
        return self.Vin + 1j * self.Vout

    @property
    def amplitude(self) -> np.ndarray:
        return np.abs(self.complex_signal)

    @property
    def phase_deg(self) -> np.ndarray:
        return np.angle(self.complex_signal, deg=True)

    @classmethod
    def load_from_file(
        cls,
        filepath: str,
        fix_legacy_labview: bool = False,
    ) -> ExperimentData:
        """
        Load experimental data from text file.

        Parameters
        ----------
        filepath : str
            Path to data file.
        fix_legacy_labview : bool
            Whether to apply 1-point frequency shift for legacy LabVIEW files (pre Oct 2017).
        """
        data = np.loadtxt(filepath)
        if data.ndim == 1:
            data = data.reshape(1, -1)
        if data.shape[0] == 7 and data.shape[1] > 7:
            data = data.T

        f = data[:, 0]
        Vin = data[:, 2]
        Vout = data[:, 3]
        V0 = data[:, 5] if data.shape[1] > 5 else None

        if fix_legacy_labview and len(f) > 1:
            f = f[1:]
            Vin = Vin[:-1]
            Vout = Vout[:-1]
            if V0 is not None:
                V0 = V0[:-1]

        return cls(frequencies=f, Vin=Vin, Vout=Vout, V0=V0)

    def apply_calibration(
        self,
        calibration: CalibrationData,
        correct_phase: bool = True,
        correct_amplitude: bool = True,
    ) -> ExperimentData:
        """
        Calibrate experimental phase shift and amplitude transfer function against calibration data.

        Parameters
        ----------
        calibration : CalibrationData
            Calibration measurement.
        correct_phase : bool
            Apply phase rotation theta(f) (default True).
        correct_amplitude : bool
            Apply normalized amplitude scaling (default True).

        Returns
        -------
        ExperimentData
            Calibrated dataset.
        """
        # Interpolate calibration phase
        f_cal = calibration.frequencies
        phase_cal = calibration.phase_deg
        interp_phase = interp1d(
            f_cal, phase_cal, kind="linear", bounds_error=False, fill_value="extrapolate"
        )
        theta_deg = interp_phase(self.frequencies)
        theta_rad = np.radians(theta_deg)

        if correct_phase:
            # Rotate complex phasor: Vin' = Vin*cos(theta) + Vout*sin(theta), Vout' = -Vin*sin(theta) + Vout*cos(theta)
            Vin_corr = self.Vin * np.cos(theta_rad) + self.Vout * np.sin(theta_rad)
            Vout_corr = -self.Vin * np.sin(theta_rad) + self.Vout * np.cos(theta_rad)
        else:
            Vin_corr = self.Vin.copy()
            Vout_corr = self.Vout.copy()

        if correct_amplitude:
            amp_cal = calibration.amplitude
            amp_cal_norm = amp_cal / np.mean(amp_cal)
            interp_amp = interp1d(
                f_cal, amp_cal_norm, kind="linear", bounds_error=False, fill_value="extrapolate"
            )
            amp_norm = interp_amp(self.frequencies)
            # Avoid divide by zero
            amp_norm = np.where(amp_norm <= 0, 1.0, amp_norm)
            Vin_corr = Vin_corr / amp_norm
            Vout_corr = Vout_corr / amp_norm

        return ExperimentData(
            frequencies=self.frequencies.copy(),
            Vin=Vin_corr,
            Vout=Vout_corr,
            V0=self.V0.copy() if self.V0 is not None else None,
        )

    def filter_frequency_range(self, f_min: float, f_max: float) -> ExperimentData:
        """
        Filter data to frequency interval [f_min, f_max].
        Ensures ascending frequency ordering.
        """
        mask = (self.frequencies >= f_min) & (self.frequencies <= f_max)
        f_sub = self.frequencies[mask]
        Vin_sub = self.Vin[mask]
        Vout_sub = self.Vout[mask]
        V0_sub = self.V0[mask] if self.V0 is not None else None

        # Ensure ascending order
        sort_idx = np.argsort(f_sub)
        return ExperimentData(
            frequencies=f_sub[sort_idx],
            Vin=Vin_sub[sort_idx],
            Vout=Vout_sub[sort_idx],
            V0=V0_sub[sort_idx] if V0_sub is not None else None,
        )
