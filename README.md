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

Pick a filter and enter the target magnitude. Choose whether to solve for the
**exposure time** needed to reach a target SNR, or the **SNR** reached in a given
exposure time. Set airmass (or altitude), seeing, binning, readout mode and a sky
brightness multiplier (1 = the derived night; raise it for moonlight). The app
reports the result, peak pixel / saturation, and diagnostic plots.
