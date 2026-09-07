"""Rubidium spectroscopy methods recovered from the surviving laboratory notes.

The CSV import, interpolation and Gaussian/Lorentzian model structure follow
gpt.pdf. Validation, bounded peak construction and output helpers complete the
missing sections. Historical measurements are shown in the report figures.
"""
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from lmfit.models import GaussianModel, LorentzianModel

SPEED_OF_LIGHT = 299_792_458.0
REPORT_CAVITY_LENGTH_M = 0.0897
# Eight fitted centres transcribed from gpt.pdf, page 15 (milliseconds to seconds).
RECOVERED_CAVITY_PEAKS_S = np.array([
    8.04217287, 13.9016315, 19.5455640, 25.0937040,
    30.6160871, 36.0764452, 41.6331799, 47.5381153,
]) * 1e-3


@dataclass(frozen=True)
class Peak:
    """Peak guess in the x-axis units; amplitude is integrated area.

    sigma is the Gaussian standard deviation or Lorentzian half-width.
    Optional centre bounds keep neighbouring components from exchanging labels.
    """
    center: float
    amplitude: float
    sigma: float
    center_bounds: tuple[float, float] | None = None


def _trace(x, y):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if x.ndim != 1 or x.size < 2 or y.shape != x.shape:
        raise ValueError("A trace needs matching one-dimensional x and y arrays.")
    if not (np.isfinite(x).all() and np.isfinite(y).all()):
        raise ValueError("Trace values must be finite.")
    if np.any(np.diff(x) <= 0):
        raise ValueError("The x axis must be strictly increasing.")
    return x, y


def load_scope(path, *, time_column="Time", signal_column="CH2",
               sample_interval_s=None, skiprows=0, time_offset_s=0.0,
               signal_scale=1.0, signal_offset=0.0):
    """Read an oscilloscope CSV with explicit column names and scale corrections.

    Set time_column=None and provide sample_interval_s when only channels were
    saved. skiprows applies before the CSV header; numeric metadata rows are
    never guessed away. Sample i is at i*dt, not at i*n*dt/(n-1).
    """
    frame = pd.read_csv(Path(path), skiprows=skiprows)
    y = pd.to_numeric(frame[signal_column], errors="raise").to_numpy(float)
    if time_column is None:
        if sample_interval_s is None or not np.isfinite(sample_interval_s) or sample_interval_s <= 0:
            raise ValueError("A positive sample interval is required without a time column.")
        x = np.arange(len(y), dtype=float) * sample_interval_s
    else:
        x = pd.to_numeric(frame[time_column], errors="raise").to_numpy(float)
    return _trace(x + time_offset_s, y * signal_scale + signal_offset)


def align_traces(*traces, sample_interval_s=None):
    """Interpolate aligned traces onto their common time interval only.

    Recorded time offsets must be applied before calling this function. By
    default use the coarsest median sample spacing, avoiding false precision.
    Interpolation does not create independent observations.
    """
    if len(traces) < 2:
        raise ValueError("Supply at least two traces.")
    traces = [_trace(*trace) for trace in traces]
    lower = max(x[0] for x, _ in traces)
    upper = min(x[-1] for x, _ in traces)
    if sample_interval_s is None:
        sample_interval_s = max(np.median(np.diff(x)) for x, _ in traces)
    if not np.isfinite(sample_interval_s) or sample_interval_s <= 0 or upper <= lower:
        raise ValueError("Traces must overlap and the sample interval must be positive.")
    count = int(np.floor((upper - lower) / sample_interval_s)) + 1
    if count < 2:
        raise ValueError("The common interval must contain at least two samples.")
    time = lower + np.arange(count) * sample_interval_s
    return time, [np.interp(time, x, y) for x, y in traces]


def separate_components(sweep, doppler, saturated, *, sample_interval_s=None):
    """Return time, Doppler-broadened signal and Doppler-free signal.

    Each input is (time_seconds, detector_signal). Subtraction signs follow
    the surviving code and Figure 3: doppler - sweep; saturated - doppler.
    """
    time, (baseline, broad, narrow) = align_traces(
        sweep, doppler, saturated, sample_interval_s=sample_interval_s
    )
    return time, broad - baseline, narrow - broad


def time_to_frequency(time_s, cavity_peaks_s, *, cavity_length_m=REPORT_CAVITY_LENGTH_M,
                      refractive_index=1.0, zero_time_s=None):
    """Piecewise-linear cavity calibration, returning relative frequency in GHz.

    Each pair of adjacent markers must be consecutive cavity resonances.
    Values outside the calibrated interval are rejected, rather than clamped.
    """
    peaks = np.asarray(cavity_peaks_s, dtype=float)
    _trace(peaks, np.zeros_like(peaks))
    t = np.asarray(time_s, dtype=float)
    if not np.isfinite(t).all() or np.any(t < peaks[0]) or np.any(t > peaks[-1]):
        raise ValueError("All times must lie between the first and last cavity markers.")
    if not np.isfinite([cavity_length_m, refractive_index]).all() or min(cavity_length_m, refractive_index) <= 0:
        raise ValueError("Cavity length and refractive index must be positive.")
    origin = peaks[0] if zero_time_s is None else zero_time_s
    if not np.isfinite(origin) or not peaks[0] <= origin <= peaks[-1]:
        raise ValueError("The frequency origin must lie inside the calibrated interval.")
    order = np.arange(peaks.size)
    fsr_ghz = SPEED_OF_LIGHT / (2 * refractive_index * cavity_length_m * 1e9)
    return (np.interp(t, peaks, order) - np.interp(origin, peaks, order)) * fsr_ghz


def fit_peaks(x, signal, peaks, *, kind="lorentzian", background_degree=2, uncertainty=None):
    """Fit a sum of line profiles and an optional polynomial background.

    Use Gaussian profiles for cavity resonances and broad troughs; use
    Lorentzians for resolved hyperfine regions. Negative amplitudes describe
    absorption troughs. Polynomial fitting uses a centred, scaled coordinate
    for numerical stability while peak parameters retain the original units.
    """
    x, signal = _trace(x, signal)
    if kind not in {"gaussian", "lorentzian"} or not peaks:
        raise ValueError("Choose gaussian or lorentzian and supply peak guesses.")
    if background_degree is not None and (not isinstance(background_degree, int) or not 0 <= background_degree <= 3):
        raise ValueError("Use no background or a polynomial degree from 0 to 3.")
    profile = GaussianModel if kind == "gaussian" else LorentzianModel
    model, parameters = None, None
    for i, peak in enumerate(peaks, 1):
        if not np.isfinite([peak.center, peak.amplitude, peak.sigma]).all() or peak.sigma <= 0 or peak.amplitude == 0:
            raise ValueError("Peak guesses need finite values, positive width and nonzero area.")
        prefix = f"peak{i}_"
        component = profile(prefix=prefix)
        guess = component.make_params(center=peak.center, amplitude=peak.amplitude, sigma=peak.sigma)
        lo, hi = peak.center_bounds or (x[0], x[-1])
        if not x[0] <= lo < hi <= x[-1] or not lo <= peak.center <= hi:
            raise ValueError("Peak centre bounds must contain the guess and lie inside the data.")
        guess[prefix + "center"].set(min=lo, max=hi)
        guess[prefix + "sigma"].set(min=np.finfo(float).eps, max=np.ptp(x))
        guess[prefix + "amplitude"].set(min=0 if peak.amplitude > 0 else -np.inf,
                                         max=np.inf if peak.amplitude > 0 else 0)
        model = component if model is None else model + component
        if parameters is None:
            parameters = guess
        else:
            parameters.update(guess)
    if background_degree is not None:
        from lmfit import Model
        origin, scale = float(np.mean(x)), float(np.ptp(x))
        def background(x, c0=0.0, c1=0.0, c2=0.0, c3=0.0):
            z = (x - origin) / scale
            return c0 + c1*z + c2*z*z + c3*z*z*z
        floor = Model(background, prefix="bg_")
        guess = floor.make_params()
        for degree in range(4):
            guess[f"bg_c{degree}"].set(vary=degree <= background_degree)
        model = model + floor
        parameters.update(guess)
    weights = None
    if uncertainty is not None:
        uncertainty = np.asarray(uncertainty, dtype=float)
        if uncertainty.shape != signal.shape or not np.isfinite(uncertainty).all() or np.any(uncertainty <= 0):
            raise ValueError("Provide a positive uncertainty for every observation.")
        weights = 1 / uncertainty
    if x.size <= sum(p.vary for p in parameters.values()):
        raise ValueError("The fit needs more observations than free parameters.")
    return model.fit(signal, parameters, x=x, weights=weights,
                     scale_covar=uncertainty is None, nan_policy="raise")


def subtract_fitted_background(x, signal, result):
    """Return baseline-corrected observations and each fitted line component."""
    x, signal = _trace(x, signal)
    components = result.eval_components(x=x)
    floor = components.pop("bg_", np.zeros_like(signal))
    return signal - floor, components


def ground_state_constant(splitting_mhz, nuclear_spin):
    """For a J=1/2 ground state, A = splitting / (I + 1/2), in MHz."""
    if not np.isfinite([splitting_mhz, nuclear_spin]).all() or min(splitting_mhz, nuclear_spin) <= 0:
        raise ValueError("Splitting and nuclear spin must be positive.")
    return splitting_mhz / (nuclear_spin + 0.5)
