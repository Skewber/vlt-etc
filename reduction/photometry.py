"""Aperture photometry + astrometric calibration for the cluster fields.

For one frame per (cluster, filter):
  1. calibrate (bias/dark/flat)            -> fits_utils.calibrate_light
  2. detect stars (DAOStarFinder)
  3. plate-solve via astrometry.net        -> WCS  (needs ASTROMETRY_API_KEY)
  4. aperture photometry (local sky)       -> instrumental mag (e-/s)
  5. query APASS (Vizier) and cross-match  -> catalog B, V, R_c, I_c
  6. compute airmass from DATE-OBS + site + cluster coordinates

Returns a flat list of matched-star records consumed by derive_zp_extinction.
"""
from __future__ import annotations

import os
import warnings

import numpy as np
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.stats import sigma_clipped_stats
from astropy.time import Time
from astropy.wcs import WCS

import fits_utils as fu

warnings.simplefilter("ignore")

# Vienna observatory (from FITS SITELAT/SITELONG headers)
SITE = EarthLocation(
    lat=(48 + 13 / 60 + 55 / 3600) * u.deg,
    lon=(16 + 20 / 60 + 1 / 3600) * u.deg,
    height=240 * u.m,
)

# Cluster centres (J2000), used as plate-solve hints and for airmass
CLUSTER_COORDS = {
    "M50": SkyCoord("07h02m42s", "-08d23m00s", frame="icrs"),
    "M67": SkyCoord("08h51m18s", "+11d48m00s", frame="icrs"),
    "NGC2281": SkyCoord("06h48m18s", "+41d05m00s", frame="icrs"),
}

PIXSCALE_BINNED = 0.5591  # arcsec / binned pixel
FWHM_PX = 3.5            # ~2" seeing at the binned scale
APER_R = 7.0            # aperture radius [px]
ANN_IN, ANN_OUT = 10.0, 16.0


def compute_airmass(hdr, cluster: str) -> float:
    t = Time(hdr["DATE-OBS"], format="isot", scale="utc")
    altaz = CLUSTER_COORDS[cluster].transform_to(AltAz(obstime=t, location=SITE))
    return float(altaz.secz)


def detect_sources(data: np.ndarray):
    from photutils.detection import DAOStarFinder

    mean, median, std = sigma_clipped_stats(data, sigma=3.0)
    finder = DAOStarFinder(fwhm=FWHM_PX, threshold=8.0 * std)
    tbl = finder(data - median)
    if tbl is None:
        return None
    tbl.sort("flux")
    tbl.reverse()  # brightest first (helps the solver)
    return tbl


def solve_wcs(tbl, shape, cluster: str):
    from astroquery.astrometry_net import AstrometryNet

    ast = AstrometryNet()
    ast.api_key = os.environ["ASTROMETRY_API_KEY"]
    c = CLUSTER_COORDS[cluster]
    try:
        hdr = ast.solve_from_source_list(
            tbl["xcentroid"], tbl["ycentroid"], shape[1], shape[0],
            scale_units="arcsecperpix", scale_type="ev",
            scale_est=PIXSCALE_BINNED, scale_err=15,
            center_ra=float(c.ra.deg), center_dec=float(c.dec.deg),
            radius=1.0, parity=2, solve_timeout=300,
        )
    except Exception as exc:  # noqa: BLE001
        print(f"    solve failed: {exc}")
        return None
    if not hdr:
        return None
    return WCS(hdr)


def photometer(data, tbl, exp, gain):
    """Return arrays (x, y, m_inst) for non-saturated, positive-flux sources."""
    from photutils.aperture import (ApertureStats, CircularAnnulus,
                                     CircularAperture, aperture_photometry)

    pos = np.transpose((tbl["xcentroid"], tbl["ycentroid"]))
    aper = CircularAperture(pos, r=APER_R)
    ann = CircularAnnulus(pos, r_in=ANN_IN, r_out=ANN_OUT)
    sky = ApertureStats(data, ann).median
    peak = ApertureStats(data, aper).max
    phot = aperture_photometry(data, aper)
    net_adu = phot["aperture_sum"] - sky * aper.area_overlap(data)
    rate_e = net_adu * gain / exp  # electrons per second
    good = (net_adu > 0) & (peak < 0.9 * fu.SATURATION_ADU) & np.isfinite(rate_e)
    x = np.asarray(tbl["xcentroid"])[good]
    y = np.asarray(tbl["ycentroid"])[good]
    m_inst = -2.5 * np.log10(np.asarray(rate_e)[good])
    return x, y, m_inst


def query_apass(center: SkyCoord, radius_deg=0.35):
    """APASS DR9 over the field, with Cousins R/I via Lupton (2005)."""
    from astroquery.vizier import Vizier

    v = Vizier(columns=["RAJ2000", "DEJ2000", "Bmag", "Vmag",
                        "g'mag", "r'mag", "i'mag"], row_limit=-1)
    res = v.query_region(center, radius=radius_deg * u.deg, catalog="II/336/apass9")
    if not res:
        return None
    t = res[0]
    g, r, i = t["g'mag"], t["r'mag"], t["i'mag"]
    cat = {
        "coord": SkyCoord(t["RAJ2000"], t["DEJ2000"], unit="deg"),
        "B": np.asarray(t["Bmag"], float),
        "V": np.asarray(t["Vmag"], float),
        # Lupton 2005 SDSS -> Cousins transforms
        "R": np.asarray(r - 0.1837 * (g - r) - 0.0971, float),
        "I": np.asarray(r - 1.2444 * (r - i) - 0.3820, float),
    }
    return cat


def process_field(cluster: str, filt: str):
    """Return list of records for one (cluster, filter), or [] on failure."""
    from astropy.coordinates import match_coordinates_sky

    files = fu.light_files(cluster, filt)
    gain = process_field.gain
    for path in files[:3]:  # try a few frames until one solves
        cal, hdr = fu.calibrate_light(path)
        tbl = detect_sources(cal)
        if tbl is None or len(tbl) < 10:
            continue
        wcs = solve_wcs(tbl, cal.shape, cluster)
        if wcs is None:
            continue
        exp = float(hdr["EXPTIME"])
        x, y, m_inst = photometer(cal, tbl, exp, gain)
        if len(x) == 0:
            continue
        src = SkyCoord.from_pixel(x, y, wcs)
        cat = query_apass(CLUSTER_COORDS[cluster])
        if cat is None:
            continue
        idx, sep, _ = match_coordinates_sky(src, cat["coord"])
        airmass = compute_airmass(hdr, cluster)
        records = []
        for k in range(len(src)):
            if sep[k].arcsec > 2.0:
                continue
            m_cat = cat[filt][idx[k]]
            if not np.isfinite(m_cat):
                continue
            records.append({
                "cluster": cluster, "filter": filt, "airmass": airmass,
                "m_inst": float(m_inst[k]), "m_cat": float(m_cat),
            })
        print(f"  {cluster:8s} {filt}: {len(records):3d} matched stars "
              f"(X={airmass:.3f}, frame={os.path.basename(path)})")
        return records
    print(f"  {cluster:8s} {filt}: FAILED (no solve/match)")
    return []


def run(gain: float) -> list[dict]:
    process_field.gain = gain
    records = []
    for filt in fu.FILTERS:
        for cluster in fu.CLUSTERS:
            records.extend(process_field(cluster, filt))
    return records


if __name__ == "__main__":
    import json
    import derive_detector
    g = derive_detector.derive()["gain_e_per_adu"]
    recs = run(g)
    print(f"\nTotal matched stars: {len(recs)}")
    with open(os.path.join(fu.HERE, "photometry_records.json"), "w") as f:
        json.dump(recs, f, indent=2)
