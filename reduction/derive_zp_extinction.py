"""Bouguer fit: derive zero point and extinction coefficient per filter.

For each filter, pool matched stars from all three clusters and fit
    y = m_inst - m_cat = C + k * X
with iterative sigma-clipping. Then:
    extinction_k        = k            (slope; should be positive)
    zeropoint_mag_1es   = -C           (magnitude giving 1 e-/s, for the ETC)

A diagnostic Bouguer plot is saved per filter under reduction/plots/.
"""
from __future__ import annotations

import os

import numpy as np

import fits_utils as fu

PLOT_DIR = os.path.join(fu.HERE, "plots")


def robust_linfit(x, y, nsig=2.5, iters=4):
    keep = np.ones(len(x), bool)
    slope = intercept = 0.0
    for _ in range(iters):
        slope, intercept = np.polyfit(x[keep], y[keep], 1)
        resid = y - (slope * x + intercept)
        std = np.std(resid[keep])
        if std == 0:
            break
        keep = np.abs(resid) < nsig * std
        if keep.sum() < 3:
            break
    resid = y[keep] - (slope * x[keep] + intercept)
    rms = float(np.sqrt(np.mean(resid ** 2))) if keep.sum() else float("nan")
    return float(slope), float(intercept), rms, keep


def _plot(filt, x, y, keep, slope, intercept, clusters):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(PLOT_DIR, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 5))
    colors = {"M50": "#4f46e5", "M67": "#f59e0b", "NGC2281": "#10b981"}
    for cl in fu.CLUSTERS:
        m = np.array([c == cl for c in clusters]) & keep
        ax.scatter(x[m], y[m], s=18, color=colors[cl], label=cl, alpha=0.7)
    rej = ~keep
    ax.scatter(x[rej], y[rej], s=14, facecolors="none", edgecolors="grey",
               label="rejected")
    xs = np.linspace(min(x), max(x), 50)
    ax.plot(xs, slope * xs + intercept, "k-",
            label=f"k={slope:.3f}, ZP={-intercept:.3f}")
    ax.set_xlabel("airmass X")
    ax.set_ylabel(r"$m_{inst} - m_{cat}$")
    ax.set_title(f"Bouguer line - filter {filt}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, f"bouguer_{filt}.png"), dpi=110)
    plt.close(fig)


def derive(records: list[dict], exclude_clusters=()) -> dict:
    if exclude_clusters:
        records = [r for r in records if r["cluster"] not in exclude_clusters]
    out = {}
    for filt in fu.FILTERS:
        recs = [r for r in records if r["filter"] == filt]
        if len(recs) < 3:
            out[filt] = {"zeropoint_mag_1es": None, "extinction_k": None,
                         "n_stars": len(recs), "fit_rms": None}
            continue
        x = np.array([r["airmass"] for r in recs])
        y = np.array([r["m_inst"] - r["m_cat"] for r in recs])
        clusters = [r["cluster"] for r in recs]
        slope, intercept, rms, keep = robust_linfit(x, y)
        _plot(filt, x, y, keep, slope, intercept, clusters)
        out[filt] = {
            "zeropoint_mag_1es": round(-intercept, 4),
            "extinction_k": round(slope, 4),
            "n_stars": int(keep.sum()),
            "fit_rms": round(rms, 4),
        }
    return out


def derive_typical(records, typical_k: dict, anchor_cluster: str) -> dict:
    """Interim option: adopt typical (literature) extinction coefficients and
    re-anchor the zero point on a well-measured reference field so the ETC still
    reproduces the actual measured throughput at that field's airmass.

    From  m_inst - m_cat = k*X - ZP  =>  ZP = k*X_anchor - offset_anchor.
    """
    out = {}
    for filt in fu.FILTERS:
        recs = [r for r in records
                if r["filter"] == filt and r["cluster"] == anchor_cluster]
        if not recs:
            out[filt] = {"zeropoint_mag_1es": None, "extinction_k": None,
                         "n_stars": 0}
            continue
        X = recs[0]["airmass"]
        offs = np.array([r["m_inst"] - r["m_cat"] for r in recs])
        med = np.median(offs)
        keep = np.abs(offs - med) < 3 * np.std(offs)
        off = float(np.median(offs[keep]))
        k = typical_k[filt]
        out[filt] = {
            "zeropoint_mag_1es": round(k * X - off, 4),
            "extinction_k": k,
            "n_stars": int(keep.sum()),
            "anchor_cluster": anchor_cluster,
            "anchor_airmass": round(float(X), 3),
        }
    return out


if __name__ == "__main__":
    import json
    with open(os.path.join(fu.HERE, "photometry_records.json")) as f:
        recs = json.load(f)
    print(json.dumps(derive(recs), indent=2))
