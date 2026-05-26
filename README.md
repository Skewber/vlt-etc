# Exposure Time Calculator — Vienna 0.8 m Telescope

An exposure time calculator (ETC) for the University of Vienna 0.8 m telescope
(f/8.3, 6640 mm focal length, FLI CCD, B/V/R/I filters). The instrumental
constants are **derived from real data** taken on 2026-02-18 (three open clusters
at different airmasses plus bias/dark/flat calibration frames), not guessed.

The calculator is a static web app in [`docs/`](docs/) that can be hosted for free
on GitHub Pages.

## Repository layout

```
docs/            static web ETC (deploy this to GitHub Pages)
  index.html
  app.js         ETC math + Chart.js plots
  style.css
  constants.json derived instrumental constants (produced by the pipeline)
  chart.min.js   vendored Chart.js (offline)
reduction/       Python pipeline that derives constants.json from data/
etc.ipynb        the original reference notebook (now reads constants.json)
data/            raw FITS frames (not committed — large)
```

## How the constants are derived

`reduction/run_all.py` produces `reduction/constants.json`:

| Quantity | Method |
|---|---|
| **Gain** | photon-transfer (mean–variance) on flat pairs |
| **Read noise** (8 MHz, 1 MHz) | difference of bias pairs |
| **Dark current** | sigma-clipped master dark / exptime |
| **Sky rate** per filter | sigma-clipped background of the cluster frames |
| **Zero point + extinction** per filter | aperture photometry on plate-solved cluster frames, cross-matched to APASS DR9, then a Bouguer fit (instrumental − catalog mag vs airmass) over the three clusters |

R, I catalog magnitudes come from APASS g,r,i via the Lupton (2005) SDSS→Cousins
transforms. The narrowband filters (OIII, Hα, SII) are **not** characterised — no
calibration data exists for them yet.

### Re-running the pipeline

```bash
python -m venv .venv
.venv/Scripts/python -m pip install -r reduction/requirements.txt   # Windows
# source .venv/bin/activate; pip install -r reduction/requirements.txt  # *nix

# astrometry.net key is needed only for the photometry (ZP + extinction) step
set ASTROMETRY_API_KEY=your_key_here            # Windows
# export ASTROMETRY_API_KEY=your_key_here       # *nix

cd reduction
python run_all.py
cp constants.json ../docs/constants.json
```

The detector and sky steps need no internet; only the photometry step contacts
astrometry.net and Vizier.

## Hosting on GitHub Pages

1. Create a GitHub repo and push this project (the `.gitignore` keeps `data/` and
   `.venv/` out).
2. In the repo: **Settings → Pages → Build and deployment**, set
   **Source = Deploy from a branch**, **Branch = `main`**, **folder = `/docs`**.
3. After a minute the ETC is live at
   `https://<your-username>.github.io/<repo-name>/`.

To preview locally:

```bash
cd docs
python -m http.server 8000
# open http://localhost:8000
```

## Using the calculator

Pick a filter and enter the target brightness. Choose whether to solve for the
**exposure time** needed to reach a target SNR, or the **SNR** reached in a given
exposure time. Set airmass (or altitude), seeing, binning, readout mode and a sky
brightness multiplier (1 = the derived night; raise it for moonlight). The app
reports the result, peak pixel / saturation, and diagnostic plots.

### Target options

- **Point source** (default): enter the magnitude in the observing filter; the
  aperture is set by the seeing.
- **Extended source** (galaxies / nebulae): choose a *light profile* — **uniform
  disk**, **Gaussian**, or **Sérsic** (with adjustable index n; n=1 ≈ exponential
  disk, n=4 ≈ de Vaucouleurs elliptical) — plus an angular size (arcsec). The
  size's meaning depends on the profile (diameter for uniform, FWHM for Gaussian,
  half-light diameter for Sérsic; shown inline). Brightness is given as a **total
  magnitude** or **mean surface brightness** (mag/arcsec²) over a circular
  **extraction aperture** (diameter in arcsec, or 0 for a single peak pixel).
  Profiles are evaluated analytically; the peak surface brightness is
  seeing-limited (mean flux within a seeing-FWHM disc).
- **Source spectrum** (optional): pick a *blackbody* temperature (point sources
  only) or a *template* and the band your magnitude is given in; the app converts
  it to the observing filter using approximate colors. Templates available:
  - **Point**: stellar main-sequence O5…M5, plus *SN Ia (max)*, *SN II (peak)*,
    *White dwarf (DA)*, *AGN / quasar*.
  - **Extended**: *Elliptical / Spiral / Irregular galaxy*, *Reflection nebula*,
    *Planetary nebula*, *HII region (Hα)*.

  This is **approximate** — a cross-band color helper, not a true
  synthetic-photometry integration (we have no per-filter throughput curves), and
  emission-line nebula colors in broadband filters are particularly rough. For an
  in-filter magnitude with no conversion, leave it on *Flat (none)*.

- **Limiting magnitude (S/N = 5)**: shown for point sources — the faintest
  in-band magnitude detectable at S/N = 5 in the same exposure / sub structure /
  conditions you've set. Solved from `S² = (Rt)² / (Rt + n_pix(t(R_S+R_D) +
  N·N_R²))`. Useful for "what's the deepest I can go tonight?".

### Sub-exposures, cooling, and the noise breakdown

- **Exposure structure**: choose how the integration is built —
  - *Single exposure*.
  - *N subs of fixed length*: give a sub length and N (total = N × length).
  - *Split total into N subs*: give a total and N (sub length = total ÷ N).

  The complementary quantity (total or sub length) is shown inline under the field.
  Read noise is paid once per sub (`N·N_read²`) and **saturation is checked per
  sub** — the relevant limit for stacking. When solving for *exposure time* from a
  target SNR, the total is solved and the structure determines the rest: with a
  fixed sub length the number of subs is derived (rounded up to whole subs); with a
  fixed N the sub length is derived.
- **Cooler temperature**: scales the dark current with
  `D(T) = D_ref · 2^((T − T_ref)/ΔT_double)`, with `D_ref` the measured dark at the
  −5 °C setpoint and `ΔT_double` a **provisional 6.3 K** doubling
  (`camera.dark_doubling_kelvin` in `constants.json`). Replace it by fitting darks
  taken at several temperatures.
- **Noise breakdown chart**: the four slices are the source-Poisson, sky, dark,
  and read-noise variance contributions, shown as **percentages** of the total
  noise variance (the term under √ in the SNR formula). Whichever dominates tells
  you what helps: source → photon-limited; sky → sky-limited; read → use longer
  subs / 1 MHz / binning; dark → cool the camera.
- **Moon contribution**: enter the moon illumination (0=new … 1=full) and the
  moon–target separation (deg). The app applies a Krisciunas & Schaefer (1991)
  V-band model to scale up the sky rate; the result is shown as `Sky × (moon)` in
  the results table. (Approximate — assumes moon at the same altitude as the
  target and uses the V-band factor for all bands.)
- **Wall-clock total**: the results also report the wall-clock time including
  per-frame readout + download (`camera.readout_time_s` in `constants.json`,
  estimated 1.0 s at 8 MHz and 8.0 s at 1 MHz — refine by timing real frames).
  Useful when planning many short subs: e.g. 100 × 5 s science integration with
  1 s readout each = 600 s on the clock, 50 % overhead.
- **Export .txt**: click *Export .txt* in the Results header to download a plain-
  text file with all inputs, the headline results, and the noise breakdown in %.
