"""Derive the sky background rate per filter from the cluster science frames.

For each calibrated cluster frame we take the sigma-clipped background median
(ADU), convert to electrons per second per binned pixel, and average over all
frames/clusters in that filter. This is representative of the conditions on the
night of 2026-02-18 (it depends on moon/sky brightness).
"""
from __future__ import annotations

import numpy as np

import fits_utils as fu


def derive(gain: float) -> dict:
    out = {}
    for filt in fu.FILTERS:
        rates = []
        for cluster in fu.CLUSTERS:
            for path in fu.light_files(cluster, filt):
                cal, hdr = fu.calibrate_light(path)
                exp = float(hdr["EXPTIME"])
                # median of the sigma-clipped image approximates the sky level
                _, median, _ = fu.robust_stats(fu.central(cal))
                rates.append(median / exp * gain)
        out[filt] = round(float(np.median(rates)), 4)
    return out


if __name__ == "__main__":
    import json
    import derive_detector
    g = derive_detector.derive()["gain_e_per_adu"]
    print(json.dumps(derive(g), indent=2))
