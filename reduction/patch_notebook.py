"""Patch etc.ipynb to load derived constants from reduction/constants.json
instead of the placeholder values, and document the corrections."""
import json
import os

NB = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "etc.ipynb")


def lines(s):
    return s.splitlines(keepends=True)


def code(s):
    return {"cell_type": "code", "metadata": {}, "execution_count": None,
            "outputs": [], "source": lines(s)}


def md(s):
    return {"cell_type": "markdown", "metadata": {}, "source": lines(s)}


nb = json.load(open(NB, encoding="utf-8"))
cells = nb["cells"]


def find(substr):
    for i, c in enumerate(cells):
        if substr in "".join(c["source"]):
            return i
    raise KeyError(substr)


# 1) replace placeholder constant cells in place ----------------------------
cells[find('read_mode_vals = {"8 MHz"')] = code(
    "# read noise per readout mode (electrons) -- derived from bias frames\n"
    "read_mode_vals = CONST['camera']['read_noise_e']\n"
    "N_readout = read_mode_vals[read_mode]")

cells[find("R_dark = 21.085")] = code(
    "# dark current (e-/s @ -5C) -- derived from dark frames\n"
    "R_dark = CONST['camera']['dark_current_e_per_s_at_-5C']")

cells[find("R_sky_vals = {'B':80")] = code(
    "# sky background rate (e-/s/binned pixel) -- derived from the cluster frames\n"
    "R_sky_vals = {f: CONST['filters'][f]['sky_rate_e_per_s_per_pix']\n"
    "              for f in CONST['filters']}\n"
    "R_sky = R_sky_vals[used_filter]")

cells[find("ZP_vals = {'B':18")] = code(
    "# zero points and extinction coefficients -- derived via Bouguer fit\n"
    "ZP_vals = {f: CONST['filters'][f]['zeropoint_mag_1es'] for f in CONST['filters']}\n"
    "k_vals  = {f: CONST['filters'][f]['extinction_k']      for f in CONST['filters']}\n\n"
    "def get_R_star(mag, filt, airmass):\n"
    "    ZP = ZP_vals[filt]\n"
    "    k = k_vals[filt]\n"
    "    return 10**(0.4*(ZP - (mag + k*airmass)))\n\n"
    "R_star = get_R_star(mag, used_filter, airmass)")

# 2) point the two "gain = 1.56" cells at the derived gain -------------------
for i, c in enumerate(cells):
    src = "".join(c["source"])
    if "gain = 1.56" in src:
        src = src.replace("gain = 1.56      # e-/ADU (from your camera data)",
                          "gain = CONST['camera']['gain_e_per_adu']")
        src = src.replace("gain = 1.56  # e-/ADU",
                          "gain = CONST['camera']['gain_e_per_adu']")
        cells[i]["source"] = lines(src)

# 3) insert a constants-loading cell right after the imports -----------------
imp = find("import numpy")
load = code(
    "# --- Derived instrumental constants (reduction/run_all.py) ---\n"
    "import json, os\n"
    "with open(os.path.join('reduction', 'constants.json')) as _f:\n"
    "    CONST = json.load(_f)\n"
    "gain = CONST['camera']['gain_e_per_adu']\n"
    "SAT_ADU = CONST['camera']['saturation_adu']\n"
    "CONST['camera']")
cells.insert(imp + 1, load)

# 4) documentation markdown at the very top ---------------------------------
cells.insert(0, md(
    "## Constants now derived from data\n"
    "The placeholder constants were replaced with values measured from the "
    "2026-02-18 run (see `reduction/`). Corrections found in the original "
    "notebook:\n\n"
    "| Constant | old (placeholder) | derived |\n"
    "|---|---|---|\n"
    "| dark current | 21.085 e-/s | ~0.38 e-/s |\n"
    "| read noise (8 MHz) | 14 e- | ~19.5 e- |\n"
    "| gain | 1.56 e-/ADU | ~1.50 e-/ADU |\n"
    "| ZP / k / R_sky | round numbers | from photometry + Bouguer fit |\n\n"
    "Run `python reduction/run_all.py` to regenerate `constants.json`. "
    "The `air_quality` input is no longer used (extinction k is now the "
    "measured per-filter value for the night)."))

json.dump(nb, open(NB, "w", encoding="utf-8"), indent=1)
print(f"Patched {NB}: now {len(cells)} cells")
