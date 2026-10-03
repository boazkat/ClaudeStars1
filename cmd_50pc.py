#!/usr/bin/env python3
"""Colour-magnitude diagram of stars within 50 pc.

Colour      : Tycho-2 BT - VT
Magnitude   : M_Ks = Ks(2MASS) + 5 log10(parallax[mas]) - 10
Distance    : Gaia DR3 parallax > 20 mas (d < 50 pc)

Cross-matches use the Gaia DR3 pre-computed best-neighbour tables
(tycho2tdsc_merge_best_neighbour and tmass_psc_xsc_best_neighbour), queried
through the Gaia archive TAP service. Results are cached to a CSV so the
plot can be regenerated offline.

Usage:
    python cmd_50pc.py                 # query (if no cache) and plot
    python cmd_50pc.py --refresh       # force a new archive query
"""
import argparse
import os

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

TAP = "https://gea.esac.esa.int/tap-server/tap"
CACHE = "data/stars_50pc_tycho2_2mass_gaia.csv"

QUERY = """
SELECT g.source_id, g.parallax, g.parallax_error, g.ruwe, g.phot_g_mean_mag,
       t.id AS tycho_id, t.hip, t.bt_mag, t.vt_mag, t.e_bt_mag, t.e_vt_mag,
       tm.ks_m, tm.ks_msigcom, tm.ph_qual
FROM gaiadr3.gaia_source AS g
JOIN gaiadr3.tycho2tdsc_merge_best_neighbour AS tx
     ON tx.source_id = g.source_id
JOIN gaiadr3.tycho2tdsc_merge AS t
     ON t.id = tx.original_ext_source_id
JOIN gaiadr3.tmass_psc_xsc_best_neighbour AS xm
     ON xm.source_id = g.source_id
JOIN gaiadr3.tmass_psc_xsc_join AS xj
     ON xm.original_ext_source_id = xj.original_psc_source_id
JOIN gaiadr1.tmass_original_valid AS tm
     ON xj.original_psc_source_id = tm.designation
WHERE g.parallax > 20
"""


def _col(tab, *names):
    """Return the first column present (Tycho column names vary by table)."""
    for n in names:
        if n in tab.colnames:
            return np.asarray(tab[n].filled(np.nan) if hasattr(tab[n], "filled")
                              else tab[n], dtype=float)
    raise KeyError(f"none of {names} in table columns {tab.colnames}")


def fetch(refresh=False):
    if os.path.exists(CACHE) and not refresh:
        return Table.read(CACHE, format="ascii.csv")
    import time
    import requests
    # Asynchronous TAP job: the sync endpoint's statement timeout is too short.
    r = requests.post(TAP + "/async",
                      data=dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv",
                                PHASE="RUN", QUERY=QUERY),
                      allow_redirects=False, timeout=60)
    r.raise_for_status()
    job = r.headers["Location"]
    while True:
        try:
            phase = requests.get(job + "/phase", timeout=120).text.strip()
        except requests.exceptions.RequestException as exc:
            print(f"polling failed ({exc}); retrying")
            time.sleep(15)
            continue
        print(f"Gaia archive job {job.rsplit('/', 1)[-1]}: {phase}")
        if phase == "COMPLETED":
            break
        if phase in ("ERROR", "ABORTED"):
            raise RuntimeError(requests.get(job + "/error", timeout=60).text)
        time.sleep(15)
    res = requests.get(job + "/results/result", timeout=600)
    res.raise_for_status()
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    with open(CACHE, "w") as f:
        f.write(res.text)
    return Table.read(CACHE, format="ascii.csv")


def select(tab, max_plx_frac_err=0.1, max_ruwe=1.4, max_color_err=0.1):
    plx = _col(tab, "parallax")
    plx_e = _col(tab, "parallax_error")
    ruwe = _col(tab, "ruwe")
    bt = _col(tab, "bt_mag", "BTmag")
    vt = _col(tab, "vt_mag", "VTmag")
    e_bt = _col(tab, "bt_mag_error", "e_bt_mag", "e_BTmag")
    e_vt = _col(tab, "vt_mag_error", "e_vt_mag", "e_VTmag")
    ks = _col(tab, "ks_m")
    qual = np.array([str(q)[2] if len(str(q)) == 3 else "X" for q in tab["ph_qual"]])

    color = bt - vt
    color_err = np.hypot(e_bt, e_vt)
    mk = ks + 5 * np.log10(plx) - 10

    good = (np.isfinite(color) & np.isfinite(mk)
            & (plx_e / plx < max_plx_frac_err)
            & (np.nan_to_num(ruwe, nan=0) < max_ruwe)
            & (color_err < max_color_err)
            & (qual == "A"))           # 2MASS Ks photometric quality A
    return color, mk, color_err, good


def plot(color, mk, good, out="cmd_50pc.png"):
    fig, ax = plt.subplots(figsize=(7, 8), dpi=150)
    ax.scatter(color[~good], mk[~good], s=3, c="0.75", lw=0,
               label=f"lower quality ({(~good).sum()})", rasterized=True)
    ax.scatter(color[good], mk[good], s=4, c="#2a5d9f", lw=0,
               label=f"good photometry & astrometry ({good.sum()})",
               rasterized=True)
    ax.set_xlabel(r"$B_T - V_T$ (Tycho-2)")
    ax.set_ylabel(r"$M_{K_s}$ (2MASS + Gaia DR3 parallax)")
    ax.set_title("Stars within 50 pc")
    ax.set_xlim(-0.4, 2.3)
    ax.set_ylim(9, -6)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right", markerscale=3, fontsize=8)
    fig.tight_layout()
    fig.savefig(out)
    print(f"wrote {out}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--refresh", action="store_true", help="re-query the archive")
    p.add_argument("--out", default="cmd_50pc.png")
    args = p.parse_args()

    tab = fetch(args.refresh)
    color, mk, _, good = select(tab)
    print(f"{len(tab)} Gaia+Tycho-2+2MASS stars with parallax > 20 mas; "
          f"{good.sum()} pass quality cuts")
    plot(color, mk, good, args.out)


if __name__ == "__main__":
    main()
