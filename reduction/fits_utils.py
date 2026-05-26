"""Shared helpers for the Vienna 0.8 m ETC reduction pipeline.

Loads FITS frames, builds master calibration frames, and provides robust
statistics on central regions (to avoid vignetting at the field edges).
"""
from __future__ import annotations

import glob
import os

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(HERE)
DATA_DIR = os.path.join(PROJECT_ROOT, "data")

FILTERS = ["B", "V", "R", "I"]
CLUSTERS = ["M50", "M67", "NGC2281"]

# Detector / optics constants (from FITS headers)
SATURATION_ADU = 65535
DARK_EXPTIME = 30.0  # seconds, all dark frames


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------
def load(path: str) -> np.ndarray:
    """Return image data as float64 in physical ADU (BZERO/BSCALE applied)."""
    with fits.open(path) as hdul:
        data = hdul[0].data.astype(np.float64)
    return data


def header(path: str) -> fits.Header:
    with fits.open(path) as hdul:
        return hdul[0].header.copy()


def find(pattern: str) -> list[str]:
    """Glob inside the data directory, sorted."""
    return sorted(glob.glob(os.path.join(DATA_DIR, pattern)))


# Frame collections -------------------------------------------------------
def bias_files(read_mode: str) -> list[str]:
    tag = "8MHz" if "8" in read_mode else "1MHz"
    return find(f"HD103095-*_bias_{tag}.fits")


def dark_files() -> list[str]:
    return find("HD103095-*_dark.fits")


def flat_files(filt: str) -> list[str]:
    return find(f"M50-*_flat{filt}.fits")


def light_files(cluster: str, filt: str) -> list[str]:
    # Note: one batch of R-long files for HD103095 lacks the underscore, but
    # cluster science frames are uniformly named "<cluster>-NNNN_light<F>.fits".
    return find(f"{cluster}-*_light{filt}.fits")


# ---------------------------------------------------------------------------
# Stacking
# ---------------------------------------------------------------------------
def median_stack(paths: list[str]) -> np.ndarray:
    """Median combine a list of frames into one master frame."""
    cube = np.stack([load(p) for p in paths], axis=0)
    return np.median(cube, axis=0)


def central(data: np.ndarray, half: int = 400) -> np.ndarray:
    """Central (2*half) square cutout, to avoid vignetted/edge pixels."""
    cy, cx = data.shape[0] // 2, data.shape[1] // 2
    return data[cy - half:cy + half, cx - half:cx + half]


def robust_stats(data: np.ndarray, sigma: float = 3.0):
    """Sigma-clipped (mean, median, std)."""
    return sigma_clipped_stats(data, sigma=sigma, maxiters=5)


# ---------------------------------------------------------------------------
# Master calibration frames (cached on disk to speed up repeated runs)
# ---------------------------------------------------------------------------
_CACHE: dict[str, np.ndarray] = {}


def master_bias(read_mode: str = "8 MHz") -> np.ndarray:
    key = f"bias_{read_mode}"
    if key not in _CACHE:
        _CACHE[key] = median_stack(bias_files(read_mode))
    return _CACHE[key]


def master_dark() -> np.ndarray:
    """Bias-subtracted master dark (ADU over DARK_EXPTIME seconds)."""
    if "dark" not in _CACHE:
        _CACHE["dark"] = median_stack(dark_files()) - master_bias("8 MHz")
    return _CACHE["dark"]


def master_flat(filt: str) -> np.ndarray:
    """Bias+dark subtracted, normalized (median == 1) master flat."""
    key = f"flat_{filt}"
    if key not in _CACHE:
        raw = median_stack(flat_files(filt))
        # scale the (bias-subtracted) dark to the flat's exposure time
        hdr = header(flat_files(filt)[0])
        flat_exp = float(hdr["EXPTIME"])
        dark_scaled = master_dark() * (flat_exp / DARK_EXPTIME)
        cal = raw - master_bias("8 MHz") - dark_scaled
        norm = np.median(central(cal))
        flat = cal / norm
        flat[flat < 0.1] = np.nan  # mask dead/edge pixels (avoid div-by-zero)
        _CACHE[key] = flat
    return _CACHE[key]


def calibrate_light(path: str) -> tuple[np.ndarray, fits.Header]:
    """Bias+dark+flat calibrate a science frame; return (data_ADU, header)."""
    data = load(path)
    hdr = header(path)
    filt = hdr["FILTER"].strip()
    exp = float(hdr["EXPTIME"])
    dark_scaled = master_dark() * (exp / DARK_EXPTIME)
    cal = (data - master_bias("8 MHz") - dark_scaled) / master_flat(filt)
    return cal, hdr
