"use strict";

// ===========================================================================
// sources.js — source spectra (SED color helper) and 2D extended-source models.
// Pure, self-contained helpers used by app.js. No dependency on constants (C);
// all instrument values are passed in as arguments.
// ===========================================================================

// ---------------------------------------------------------------------------
// SED color helper (APPROXIMATE — no instrument throughput curves available).
// Converts a magnitude given in one band to the observing band via colors.
// ---------------------------------------------------------------------------
const SED = (() => {
  const EFF = { B: 445, V: 551, R: 658, I: 806 }; // effective wavelength [nm]
  const T_VEGA = 9700; // K  (defines the zero-color point)
  const HCK_NM = 1.438777e7; // hc/k in nm*K

  // Planck B_lambda up to constants that cancel in a ratio.
  function planck(lam_nm, T) {
    const x = HCK_NM / (lam_nm * T);
    return 1.0 / (Math.pow(lam_nm, 5) * (Math.exp(x) - 1.0));
  }

  // Blackbody color (m_x - m_y), Vega-zeroed so an ~A0V star has zero colors.
  function colorBB(x, y, T) {
    const raw = -2.5 * Math.log10(planck(EFF[x], T) / planck(EFF[y], T));
    const veg = -2.5 * Math.log10(planck(EFF[x], T_VEGA) / planck(EFF[y], T_VEGA));
    return raw - veg;
  }

  // Approximate broadband colors (B-V, V-Rc, V-Ic) for each template. Sources:
  //   * stellar main sequence  : Pecaut & Mamajek (2013)
  //   * SN Ia (near max)       : Hsiao+ 2007 templates, Filippenko 1997 review
  //   * SN II (near peak)      : Filippenko 1997
  //   * white dwarf (DA, ~10 kK): synthetic-photometry estimates
  //   * AGN / quasar composite : Vanden Berk+ 2001
  //   * galaxies               : Fukugita+ 1995
  //   * reflection nebula      : scattered starlight, illuminator ~ B-type
  //   * planetary nebula       : line-dominated, broadband-averaged (very rough)
  //   * HII region (Hα)        : strong Hα emission boosts R band
  // Values for emission-line objects are filter-dependent and approximate.
  const SPEC = {
    "O5": [-0.33, -0.15, -0.32], "B0": [-0.30, -0.13, -0.29],
    "B5": [-0.16, -0.06, -0.16], "A0": [ 0.00,  0.00,  0.00],
    "A5": [ 0.15,  0.10,  0.21], "F0": [ 0.30,  0.18,  0.38],
    "F5": [ 0.44,  0.26,  0.53], "G0": [ 0.59,  0.33,  0.66],
    "G2 (Sun)": [ 0.63,  0.35,  0.70], "K0": [ 0.82,  0.43,  0.86],
    "K5": [ 1.15,  0.63,  1.28], "M0": [ 1.40,  0.78,  1.62],
    "M2": [ 1.50,  0.95,  1.95], "M5": [ 1.64,  1.20,  2.85],
    "SN Ia (max)":      [ 0.05, -0.05, -0.05],
    "SN II (peak)":     [ 0.10,  0.05,  0.10],
    "White dwarf (DA)": [ 0.05,  0.00,  0.05],
    "AGN / quasar":     [ 0.30,  0.40,  0.70],
    // extended-object templates
    "Elliptical galaxy": [ 0.95,  0.55,  1.10],
    "Spiral galaxy":     [ 0.65,  0.45,  0.90],
    "Irregular galaxy":  [ 0.40,  0.35,  0.65],
    "Reflection nebula": [ 0.10,  0.05,  0.10],
    "Planetary nebula":  [ 0.10, -0.40, -0.30],
    "HII region (Hα)":   [ 0.00, -0.70, -0.50],
  };
  const POINT_TYPES = [
    "O5","B0","B5","A0","A5","F0","F5","G0","G2 (Sun)","K0","K5","M0","M2","M5",
    "SN Ia (max)","SN II (peak)","White dwarf (DA)","AGN / quasar",
  ];
  const EXT_TYPES = [
    "Elliptical galaxy","Spiral galaxy","Irregular galaxy",
    "Reflection nebula","Planetary nebula","HII region (Hα)",
  ];

  function specColors(type) {
    const row = SPEC[type] || SPEC["A0"];
    return { BV: row[0], VR: row[1], VI: row[2] };
  }

  // Colors (relative to V) for the chosen SED.
  function colorsFor(sed) {
    if (sed.mode === "bb") {
      return { BV: colorBB("B", "V", sed.T), VR: colorBB("V", "R", sed.T),
               VI: colorBB("V", "I", sed.T) };
    }
    if (sed.mode === "spec") return specColors(sed.type);
    return { BV: 0, VR: 0, VI: 0 }; // flat
  }

  // Magnitude in obsBand given a magnitude in refBand and an SED. V is the hub.
  function bandMag(mRef, refBand, obsBand, sed) {
    const c = colorsFor(sed);
    const toV = { B: -c.BV, V: 0, R: c.VR, I: c.VI };   // mV = mRef + toV[refBand]
    const fromV = { B: c.BV, V: 0, R: -c.VR, I: -c.VI }; // mObs = mV + fromV[obsBand]
    const mV = mRef + toV[refBand];
    const mObs = mV + fromV[obsBand];
    return { mObs, mV, colors: c };
  }

  return { EFF, colorBB, specColors, colorsFor, bandMag,
           pointTypes: POINT_TYPES, extTypes: EXT_TYPES };
})();

// ---------------------------------------------------------------------------
// Extended source profiles — radially symmetric analytic light distributions.
// SRC.encFracBuilder(profile, sizeArcsec, n) returns a function r -> fraction
// of total flux enclosed within radius r. `size` has a profile-dependent
// meaning (see SRC.SIZE_HINT).
// ---------------------------------------------------------------------------
const SRC = (() => {
  const SIZE_HINT = {
    uniform: "diameter",
    gaussian: "FWHM",
    sersic: "half-light diameter (2 Re)",
  };

  // Uniform disk of diameter `size`.
  function buildUniform(size) {
    const R = size / 2;
    return (r) => (r <= 0 ? 0 : (r >= R ? 1 : (r / R) * (r / R)));
  }

  // Gaussian with FWHM = size.
  function buildGaussian(size) {
    const sigma = size / 2.3548;
    const s2 = 2 * sigma * sigma;
    return (r) => 1 - Math.exp(-r * r / s2);
  }

  // Sérsic n: I(r) = exp(-bn * ((r/Re)^(1/n) - 1)), half-light radius Re = size/2.
  // Enclosed flux is precomputed on a grid out to 10 Re; lookup is linear.
  function buildSersic(size, n) {
    const Re = Math.max(size / 2, 1e-6);
    const bn = 2 * n - 1 / 3 + 4 / (405 * n) + 46 / (25515 * n * n);
    const Rmax = Re * 12;
    const N = 600;
    const dr = Rmax / N;
    const cum = new Float64Array(N + 1);
    let prev = 0;
    for (let i = 1; i <= N; i++) {
      const rho = i * dr;
      const I = Math.exp(-bn * (Math.pow(rho / Re, 1 / n) - 1));
      const integ = 2 * Math.PI * rho * I;
      cum[i] = cum[i - 1] + 0.5 * (prev + integ) * dr;
      prev = integ;
    }
    const total = cum[N];
    return (r) => {
      if (r <= 0) return 0;
      if (r >= Rmax) return 1;
      const idx = r / dr;
      const i = Math.floor(idx), f = idx - i;
      const v = cum[i] * (1 - f) + cum[Math.min(N, i + 1)] * f;
      return Math.min(1, v / total);
    };
  }

  function encFracBuilder(profile, sizeArcsec, sersicN) {
    if (profile === "gaussian") return buildGaussian(sizeArcsec);
    if (profile === "sersic") return buildSersic(sizeArcsec, sersicN || 4);
    return buildUniform(sizeArcsec); // default uniform
  }

  return { encFracBuilder, SIZE_HINT };
})();
