"""Derive detector constants: read noise, gain (PTC), dark current.

All three come straight from the calibration frames; no internet needed.

Method
------
* Read noise: for disjoint bias pairs, RN = std(b1 - b2) / sqrt(2)  [ADU],
  converted to electrons with the measured gain.
* Gain (photon transfer): for flat pairs, with the level rescaled to remove
  twilight drift,  gain = S / (var(f1 - f2)/2 - RN^2)  [e-/ADU],
  where S is the bias-subtracted mean signal.
* Dark current: sigma-clipped mean of the master dark (ADU) / exptime * gain.
"""
from __future__ import annotations

import numpy as np

import fits_utils as fu


def read_noise_adu(read_mode: str) -> float:
    """Read noise in ADU from disjoint bias-frame pairs."""
    files = fu.bias_files(read_mode)
    rns = []
    for i in range(0, len(files) - 1, 2):
        b1 = fu.central(fu.load(files[i]))
        b2 = fu.central(fu.load(files[i + 1]))
        _, _, std = fu.robust_stats(b1 - b2)
        rns.append(std / np.sqrt(2.0))
    return float(np.median(rns))


def gain_e_per_adu(rn_adu_8mhz: float) -> tuple[float, dict]:
    """Photon-transfer gain, median over all flat pairs and filters."""
    bias_mean, _, _ = fu.robust_stats(fu.central(fu.master_bias("8 MHz")))
    rn_var = rn_adu_8mhz ** 2
    per_filter = {}
    all_gains = []
    for filt in fu.FILTERS:
        files = fu.flat_files(filt)
        gains = []
        for i in range(0, len(files) - 1, 2):
            f1 = fu.central(fu.load(files[i]))
            f2 = fu.central(fu.load(files[i + 1]))
            m1 = np.median(f1)
            m2 = np.median(f2)
            # rescale f2 to f1's level to remove twilight drift (keeps shot noise)
            f2s = f2 * (m1 / m2)
            signal = 0.5 * (m1 + m2) - bias_mean
            _, _, std = fu.robust_stats(f1 - f2s)
            var_shot = std ** 2 / 2.0 - rn_var
            if var_shot > 0:
                gains.append(signal / var_shot)
        if gains:
            per_filter[filt] = float(np.median(gains))
            all_gains.extend(gains)
    return float(np.median(all_gains)), per_filter


def dark_current_e_per_s(gain: float) -> float:
    mean, _, _ = fu.robust_stats(fu.central(fu.master_dark()))
    return float(mean / fu.DARK_EXPTIME * gain)


def derive() -> dict:
    rn8 = read_noise_adu("8 MHz")
    rn1 = read_noise_adu("1 MHz")
    gain, gain_per_filter = gain_e_per_adu(rn8)
    dark = dark_current_e_per_s(gain)
    result = {
        "gain_e_per_adu": round(gain, 4),
        "gain_per_filter": {k: round(v, 4) for k, v in gain_per_filter.items()},
        "read_noise_e": {
            "8 MHz": round(rn8 * gain, 3),
            "1 MHz": round(rn1 * gain, 3),
        },
        "read_noise_adu": {"8 MHz": round(rn8, 3), "1 MHz": round(rn1, 3)},
        "dark_current_e_per_s_at_-5C": round(dark, 4),
        "saturation_adu": fu.SATURATION_ADU,
    }
    return result


if __name__ == "__main__":
    import json
    print(json.dumps(derive(), indent=2))
