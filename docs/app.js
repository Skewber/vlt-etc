"use strict";

// Loaded from constants.json (derived from the 2026-02-18 observing run).
let C = null;

// ---------------------------------------------------------------------------
// Core ETC math  (mirrors the corrected notebook)
// ---------------------------------------------------------------------------

// Number of pixels in the photometry aperture.
// `apDiam` is the aperture DIAMETER in units of FWHM.
function nPix(fwhm, binning, apDiam) {
  const sBinned = C.pixel_scale_arcsec_unbinned * binning;
  const r = (apDiam / 2) * fwhm / sBinned; // radius in binned px
  return Math.PI * r * r;
}

// Per-binned-pixel rates, scaled from the bin-2 measurement.
function skyRate(filter, binning, skyMult) {
  return C.filters[filter].sky_rate_e_per_s_per_pix *
    Math.pow(binning / 2, 2) * skyMult;
}
function darkRate(binning) {
  return C.camera["dark_current_e_per_s_at_-5C"] * Math.pow(binning / 2, 2);
}

// Star electron rate [e-/s] from magnitude, zero point and extinction.
function starRate(mag, filter, airmass) {
  const f = C.filters[filter];
  return Math.pow(10, 0.4 * (f.zeropoint_mag_1es - (mag + f.extinction_k * airmass)));
}

function snrForTime(Rstar, t, np, Rsky, Rdark, Nread) {
  const signal = Rstar * t;
  const noiseVar = signal + np * (t * (Rsky + Rdark) + Nread * Nread);
  return signal / Math.sqrt(noiseVar);
}

// Solve quadratic  R*^2 t^2 - SNR^2 (R* + npix(Rsky+Rdark)) t - SNR^2 npix Nread^2 = 0
function timeForSnr(Rstar, snr, np, Rsky, Rdark, Nread) {
  const a = Rstar * Rstar;
  const b = -(snr * snr) * (Rstar + np * (Rsky + Rdark));
  const c = -(snr * snr) * np * Nread * Nread;
  const disc = b * b - 4 * a * c;
  if (a <= 0 || disc < 0) return NaN;
  return (-b + Math.sqrt(disc)) / (2 * a);
}

// Peak pixel value (e-) at exposure t, assuming a Gaussian PSF.
function peakElectrons(Rstar, t, fwhm, binning, Rsky, Rdark) {
  const sBinned = C.pixel_scale_arcsec_unbinned * binning;
  const sigmaPx = (fwhm / sBinned) / 2.3548;
  const peakStarRate = Rstar / (2 * Math.PI * sigmaPx * sigmaPx);
  return (peakStarRate + Rsky + Rdark) * t;
}

// Gather everything for a given set of inputs.
function compute(inp) {
  const np = nPix(inp.seeing, inp.binning, inp.aperture);
  const Rsky = skyRate(inp.filter, inp.binning, inp.skymult);
  const Rdark = darkRate(inp.binning);
  const Nread = C.camera.read_noise_e[inp.readmode];
  const Rstar = starRate(inp.mag, inp.filter, inp.airmass);

  let t, snr;
  if (inp.mode === "time") {
    snr = inp.snr;
    t = timeForSnr(Rstar, snr, np, Rsky, Rdark, Nread);
  } else {
    t = inp.exptime;
    snr = snrForTime(Rstar, t, np, Rsky, Rdark, Nread);
  }

  const gain = C.camera.gain_e_per_adu;
  const peakE = peakElectrons(Rstar, t, inp.seeing, inp.binning, Rsky, Rdark);
  const peakADU = peakE / gain;
  const sBinned = C.pixel_scale_arcsec_unbinned * inp.binning;

  return {
    np, Rsky, Rdark, Nread, Rstar, t, snr, gain, peakE, peakADU,
    satPct: peakADU / C.camera.saturation_adu * 100,
    fwhmSampling: inp.seeing / sBinned,
    signalE: Rstar * t,
    skyPerPixE: Rsky * t,
    noise: {
      star: Rstar * t,
      sky: np * Rsky * t,
      dark: np * Rdark * t,
      read: np * Nread * Nread,
    },
  };
}

// ---------------------------------------------------------------------------
// UI plumbing
// ---------------------------------------------------------------------------
function readInputs() {
  const mode = document.querySelector('input[name="mode"]:checked').value;
  return {
    filter: val("filter"),
    mag: num("mag"),
    mode,
    snr: num("snr"),
    exptime: num("exptime"),
    airmass: num("airmass"),
    seeing: num("seeing"),
    skymult: num("skymult"),
    binning: parseInt(val("binning"), 10),
    readmode: val("readmode"),
    aperture: num("aperture"),
  };
}
const val = (id) => document.getElementById(id).value;
const num = (id) => parseFloat(document.getElementById(id).value);

function fmtTime(t) {
  if (!isFinite(t) || t <= 0) return "—";
  if (t < 1) return t.toFixed(3) + " s";
  if (t < 120) return t.toFixed(1) + " s";
  if (t < 3600) return (t / 60).toFixed(1) + " min";
  return (t / 3600).toFixed(2) + " h";
}

function renderResults(inp, r) {
  const label = document.getElementById("headline-label");
  const value = document.getElementById("headline-value");
  if (inp.mode === "time") {
    label.textContent = `Exposure time for SNR ${inp.snr}`;
    value.textContent = fmtTime(r.t);
  } else {
    label.textContent = `SNR in ${fmtTime(r.t)}`;
    value.textContent = isFinite(r.snr) ? r.snr.toFixed(1) : "—";
  }

  const rows = [
    ["Filter", inp.filter],
    ["Star rate R&#42;", r.Rstar.toExponential(3) + " e&#8315;/s"],
    ["Exposure time", fmtTime(r.t)],
    ["SNR", isFinite(r.snr) ? r.snr.toFixed(2) : "—"],
    ["Peak pixel", `${r.peakE.toFixed(0)} e&#8315; (${r.peakADU.toFixed(0)} ADU)`],
    ["Saturation", r.satPct.toFixed(1) + " %"],
    ["FWHM sampling", r.fwhmSampling.toFixed(2) + " px"],
    ["Aperture area", r.np.toFixed(1) + " px"],
    ["Signal in aperture", `${r.signalE.toFixed(0)} e&#8315; (${(r.signalE / r.gain).toFixed(0)} ADU)`],
    ["Sky / pixel", `${r.skyPerPixE.toFixed(1)} e&#8315; (${(r.skyPerPixE / r.gain).toFixed(1)} ADU)`],
  ];
  document.getElementById("result-table").innerHTML =
    rows.map(([k, v]) => `<tr><td>${k}</td><td>${v}</td></tr>`).join("");

  const warn = document.getElementById("sat-warning");
  warn.classList.remove("hidden", "ok");
  if (r.satPct >= 100) {
    warn.textContent = `⚠ Saturated: peak ≈ ${r.peakADU.toFixed(0)} ADU exceeds the ` +
      `${C.camera.saturation_adu} ADU full well. Shorten the exposure or use a brighter-star setting.`;
  } else if (r.satPct >= 70) {
    warn.textContent = `⚠ Near saturation (${r.satPct.toFixed(0)} % of full well). Stay below ~70 % for linear photometry.`;
  } else {
    warn.classList.add("ok");
    warn.textContent = `✓ Peak at ${r.satPct.toFixed(0)} % of full well — safely in the linear regime.`;
  }
}

// ---------------------------------------------------------------------------
// Charts
// ---------------------------------------------------------------------------
const charts = {};
const COL = { line: "#6366f1", pt: "#22d3ee", grid: "#2b3040", txt: "#9aa6bf" };

function baseOpts(xlabel, ylabel) {
  return {
    responsive: true, maintainAspectRatio: false,
    plugins: { legend: { display: false },
      title: { display: true, color: COL.txt, text: `${ylabel} vs ${xlabel}` } },
    scales: {
      x: { title: { display: true, text: xlabel, color: COL.txt },
           ticks: { color: COL.txt }, grid: { color: COL.grid } },
      y: { title: { display: true, text: ylabel, color: COL.txt },
           ticks: { color: COL.txt }, grid: { color: COL.grid } },
    },
  };
}

function lineChart(id, xs, ys, xlabel, ylabel, marker) {
  const opts = baseOpts(xlabel, ylabel);
  opts.parsing = false;
  opts.scales.x.type = "linear";
  const datasets = [
    { data: xs.map((x, i) => ({ x, y: ys[i] })), borderColor: COL.line,
      borderWidth: 2, pointRadius: 0, tension: .2 },
  ];
  if (marker && isFinite(marker.x) && isFinite(marker.y)) {
    datasets.push({ data: [{ x: marker.x, y: marker.y }], showLine: false,
      pointBackgroundColor: COL.pt, pointRadius: 5 });
  }
  if (charts[id]) charts[id].destroy();
  charts[id] = new Chart(document.getElementById(id),
    { type: "line", data: { datasets }, options: opts });
}

function updateCharts(inp, r) {
  const np = r.np, Rsky = r.Rsky, Rdark = r.Rdark, Nread = r.Nread;

  // SNR vs exposure time
  const tMax = Math.max(r.t * 2.5, 5);
  const ts = lin(0.1, tMax, 80);
  lineChart("chart-time",
    ts, ts.map((t) => snrForTime(r.Rstar, t, np, Rsky, Rdark, Nread)),
    "exposure time [s]", "SNR", { x: r.t, y: r.snr });

  // SNR vs magnitude (fixed time r.t)
  const mags = lin(inp.mag - 4, inp.mag + 4, 80);
  lineChart("chart-mag",
    mags, mags.map((m) => snrForTime(starRate(m, inp.filter, inp.airmass), r.t, np, Rsky, Rdark, Nread)),
    "magnitude", "SNR", { x: inp.mag, y: r.snr });

  // SNR vs airmass (fixed time r.t)
  const ams = lin(1.0, 2.5, 60);
  lineChart("chart-airmass",
    ams, ams.map((X) => snrForTime(starRate(inp.mag, inp.filter, X), r.t, np, Rsky, Rdark, Nread)),
    "airmass", "SNR", { x: inp.airmass, y: r.snr });

  // Peak ADU vs time (saturation)
  const adus = ts.map((t) => peakElectrons(r.Rstar, t, inp.seeing, inp.binning, Rsky, Rdark) / r.gain);
  satChart("chart-sat", ts, adus, r.t, r.peakADU);

  // Noise breakdown doughnut
  noiseChart("chart-noise", r.noise);
}

function satChart(id, xs, ys, tMark, aduMark) {
  const opts = baseOpts("exposure time [s]", "peak [ADU]");
  opts.parsing = false;
  opts.scales.x.type = "linear";
  const sat = C.camera.saturation_adu;
  const data = {
    datasets: [
      { data: xs.map((x, i) => ({ x, y: ys[i] })), borderColor: COL.line,
        borderWidth: 2, pointRadius: 0, tension: .2 },
      { data: [{ x: tMark, y: aduMark }], showLine: false,
        pointBackgroundColor: COL.pt, pointRadius: 5 },
      { data: [{ x: xs[0], y: sat }, { x: xs[xs.length - 1], y: sat }],
        borderColor: "#ef4444", borderWidth: 1.5, borderDash: [6, 4], pointRadius: 0 },
    ],
  };
  if (charts[id]) charts[id].destroy();
  charts[id] = new Chart(document.getElementById(id), { type: "line", data, options: opts });
}

function noiseChart(id, noise) {
  const labels = ["Star (Poisson)", "Sky", "Dark current", "Read noise"];
  const vals = [noise.star, noise.sky, noise.dark, noise.read];
  const opts = {
    responsive: true, maintainAspectRatio: false,
    plugins: {
      legend: { position: "right", labels: { color: COL.txt, boxWidth: 12 } },
      title: { display: true, color: COL.txt, text: "Noise variance breakdown" },
    },
  };
  if (charts[id]) charts[id].destroy();
  charts[id] = new Chart(document.getElementById(id), {
    type: "doughnut",
    data: { labels, datasets: [{ data: vals,
      backgroundColor: ["#6366f1", "#9aa6bf", "#f59e0b", "#ef4444"], borderWidth: 0 }] },
    options: opts,
  });
}

function lin(a, b, n) {
  const out = [];
  for (let i = 0; i < n; i++) out.push(a + (b - a) * i / (n - 1));
  return out;
}

// ---------------------------------------------------------------------------
// Wiring
// ---------------------------------------------------------------------------
function recalc() {
  if (!C) return;
  const inp = readInputs();
  const r = compute(inp);
  renderResults(inp, r);
  updateCharts(inp, r);
}

function setupMode() {
  document.querySelectorAll('input[name="mode"]').forEach((el) =>
    el.addEventListener("change", () => {
      const mode = document.querySelector('input[name="mode"]:checked').value;
      document.getElementById("snr-input").classList.toggle("hidden", mode !== "time");
      document.getElementById("time-input").classList.toggle("hidden", mode !== "snr");
      recalc();
    }));
}

function setupAltitudeLink() {
  const alt = document.getElementById("altitude");
  const air = document.getElementById("airmass");
  alt.addEventListener("input", () => {
    const a = parseFloat(alt.value);
    if (isFinite(a) && a > 0 && a <= 90) {
      air.value = (1 / Math.cos((90 - a) * Math.PI / 180)).toFixed(3);
      recalc();
    }
  });
  air.addEventListener("input", () => { alt.value = ""; });
}

async function init() {
  C = await (await fetch("constants.json")).json();
  const sel = document.getElementById("filter");
  sel.innerHTML = Object.keys(C.filters)
    .filter((f) => C.filters[f].zeropoint_mag_1es != null)
    .map((f) => `<option value="${f}">${f}</option>`).join("");
  sel.value = C.filters.V ? "V" : sel.options[0].value;

  setupMode();
  setupAltitudeLink();
  document.querySelectorAll("input, select").forEach((el) =>
    el.addEventListener("input", recalc));
  recalc();
}

init();
