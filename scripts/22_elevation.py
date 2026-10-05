"""Step 5b: flood area and flood fraction by elevation band.

Implements rules E1-E8 of scripts/step5_rules.md with the existing DEM
(Copernicus GLO-30 on the Landsat grid), the step 3 flood classes and the
step 3 analysis area. Nothing is resampled; all three rasters share one grid.

Descriptive only: elevation is not presented as a cause of flooding.
Copernicus GLO-30 is a surface model (heights include buildings and
vegetation), its vertical datum is not stated in the files and is not
established here, and the analysis stops at 1,800 m because the analysis
area was defined from this DEM.

Outputs: outputs/elevation_bands.csv, outputs/elevation.json,
outputs/figures/elevation.png

Usage:  .venv/bin/python scripts/22_elevation.py
"""

import csv
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from rasterio.warp import transform as warp_transform
from scipy import ndimage

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "outputs"
PIXEL_KM2 = 0.0009
DEM_NODATA = -9999.0
LOWER_EDGE = 1570.0                                          # E3
UPPER_EDGE = 1800.0
EDGES = [LOWER_EDGE + 10 * k for k in range(0, 24)]          # 1570, 1580, ..., 1800
SMALL_BAND_KM2 = 10.0                                        # E4
PCTS = [1, 5, 10, 25, 75, 90, 95, 99]                        # E5
FLOODED, PERSISTENT, NA = 1, 2, 255
FLAT_MIN_KM2 = 1.0          # post-hoc diagnostic (not in step5_rules.md): flat patch size


def stats(v):
    v = v.astype(np.float64)
    out = {"pixels": int(v.size), "min_m": round(float(v.min()), 2), "max_m": round(float(v.max()), 2),
           "mean_m": round(float(v.mean()), 2), "median_m": round(float(np.median(v)), 2),
           "std_m": round(float(v.std()), 2)}
    out.update({f"p{p}_m": round(float(x), 2) for p, x in zip(PCTS, np.percentile(v, PCTS))})
    return out


def main():
    (OUT / "figures").mkdir(parents=True, exist_ok=True)
    fr = json.loads((ROOT / "data" / "flood_report.json").read_text())
    srcs = {}
    for name, rel in (("dem", "data/clipped/dem_cop30_bilinear.tif"),
                      ("cls", "data/flood/flood_classes.tif"),
                      ("area", "data/flood/analysis_area.tif")):
        with rasterio.open(ROOT / rel) as s:
            srcs[name] = (s.read(1), s.crs.to_string(), tuple(s.transform)[:6], s.shape)
    grids = {(c, t, sh) for _, c, t, sh in srcs.values()}
    if len(grids) != 1:
        sys.exit(f"STOP: DEM, flood classes and analysis area are not on one grid: {grids}")
    dem, cls = srcs["dem"][0], srcs["cls"][0]
    area = srcs["area"][0].astype(bool)
    if round(int(area.sum()) * PIXEL_KM2, 4) != fr["analysis_area_km2"]:
        sys.exit("STOP: analysis area differs from step 3")
    nodata = area & (dem == DEM_NODATA)
    if nodata.any():
        sys.exit(f"STOP: {int(nodata.sum())} analysis-area pixels have no DEM value")
    if dem[area].max() > UPPER_EDGE:
        sys.exit("STOP: analysis-area elevation above 1,800 m; the analysis area is not as recorded")

    assessed = area & np.isin(cls, (0, 1, 2, 3))
    flooded = area & (cls == FLOODED)
    persistent = area & (cls == PERSISTENT)
    not_assessed = area & (cls == NA)
    if round(int(flooded.sum()) * PIXEL_KM2, 4) != fr["primary"]["flooded_km2"]:
        sys.exit("STOP: flooded pixels do not reproduce the step 3 result")

    lowest = float(dem[area].min())
    bands = [(lowest, LOWER_EDGE, f"below {LOWER_EDGE:.0f} m ({lowest:.1f}-{LOWER_EDGE:.0f})")]
    bands += [(a, b, f"{a:.0f}-{b:.0f} m") for a, b in zip(EDGES[:-1], EDGES[1:])]
    total_flooded = int(flooded.sum())
    rows = []
    for i, (lo, hi, label) in enumerate(bands):
        last = i == len(bands) - 1
        inb = area & (dem >= lo) & ((dem <= hi) if last else (dem < hi))
        n_area, n_ass = int(inb.sum()), int((inb & assessed).sum())
        n_fl = int((inb & flooded).sum())
        rows.append({
            "band": label, "lower_m": round(lo, 2), "upper_m": hi,
            "upper_inclusive": last,
            "analysis_area_pixels": n_area, "assessed_pixels": n_ass,
            "not_assessed_pixels": int((inb & not_assessed).sum()),
            "flooded_pixels": n_fl, "flooded_km2": round(n_fl * PIXEL_KM2, 4),
            "assessed_km2": round(n_ass * PIXEL_KM2, 4),
            "flood_fraction_pct": round(100 * n_fl / n_ass, 3) if n_ass else None,
            "persistent_water_pixels": int((inb & persistent).sum()),
            "share_of_all_flooded_pct": round(100 * n_fl / total_flooded, 3),
            "small_band": n_ass * PIXEL_KM2 < SMALL_BAND_KM2})
    if sum(r["analysis_area_pixels"] for r in rows) != int(area.sum()):
        sys.exit("STOP: bands do not cover the analysis area exactly once")
    if sum(r["flooded_pixels"] for r in rows) != total_flooded:
        sys.exit("STOP: bands do not account for every flooded pixel")

    # Post-hoc diagnostic, added after the band results were seen and NOT part
    # of step5_rules.md: contiguous patches (8-connected) of exactly constant
    # DEM value covering at least FLAT_MIN_KM2. Large flat patches inside
    # mapped lakes are edited water surfaces, not ground heights.
    osm = rasterio.open(ROOT / "data" / "osm" / "water_wetland_mask.tif").read(1).astype(bool)
    vals, counts = np.unique(dem[area], return_counts=True)
    flat = np.zeros(area.shape, bool)
    patches = []
    with rasterio.open(ROOT / "data" / "clipped" / "dem_cop30_bilinear.tif") as s0:
        T, crs = s0.transform, s0.crs
    for val in vals[counts * PIXEL_KM2 >= FLAT_MIN_KM2]:
        lab, n = ndimage.label(area & (dem == val), np.ones((3, 3)))
        sizes = np.bincount(lab.ravel())
        for k in np.nonzero(sizes * PIXEL_KM2 >= FLAT_MIN_KM2)[0]:
            if k == 0:
                continue
            m = lab == k
            flat |= m
            rr, cc = np.nonzero(m)
            lon, lat = warp_transform(crs, "EPSG:4326", [T.c + (cc.mean() + .5) * T.a],
                                      [T.f + (rr.mean() + .5) * T.e])
            patches.append({"elevation_m": float(val), "km2": round(int(m.sum()) * PIXEL_KM2, 4),
                            "centroid_lat_lon": [round(lat[0], 3), round(lon[0], 3)],
                            "flooded_km2": round(int((m & flooded).sum()) * PIXEL_KM2, 4),
                            "persistent_km2": round(int((m & persistent).sum()) * PIXEL_KM2, 4),
                            "inside_osm_water_wetland_pct": round(100 * float((m & osm).sum()) / m.sum(), 1)})
    patches.sort(key=lambda d: -d["km2"])
    high = flooded & (dem >= 1600)
    hl, hn = ndimage.label(high, np.ones((3, 3)))
    hs = np.bincount(hl.ravel())[1:]
    diagnostics = {
        "note": "Post-hoc diagnostics added after the band results were seen; not part of step5_rules.md. "
                "They do not change the band results above.",
        "flat_surfaces": {"rule": f"8-connected patches of identical DEM value, at least {FLAT_MIN_KM2:g} km2",
                          "patches": patches,
                          "total_km2": round(int(flat.sum()) * PIXEL_KM2, 4),
                          "flooded_on_flat_km2": round(int((flat & flooded).sum()) * PIXEL_KM2, 4),
                          "flooded_on_flat_pct_of_footprint": round(100 * float((flat & flooded).sum()) / total_flooded, 2),
                          "elevation_flooded_excluding_flat": stats(dem[flooded & ~flat]),
                          "elevation_assessed_excluding_flat": stats(dem[assessed & ~flat])},
        "flooded_at_or_above_1600m": {"km2": round(int(high.sum()) * PIXEL_KM2, 4),
                                      "pct_of_footprint": round(100 * float(high.sum()) / total_flooded, 2),
                                      "patches": int(hn), "median_patch_pixels": int(np.median(hs)) if hn else 0,
                                      "inside_osm_water_wetland_km2": round(int((high & osm).sum()) * PIXEL_KM2, 4)},
    }

    result = {
        "script": "scripts/22_elevation.py", "rules": "scripts/step5_rules.md (B)",
        "inputs": {"dem": "data/clipped/dem_cop30_bilinear.tif (Copernicus GLO-30, bilinear onto the "
                          "Landsat grid in step 2; see METHODS.md section 1)",
                   "flood_classes": "data/flood/flood_classes.tif (step 3)",
                   "analysis_area": "data/flood/analysis_area.tif (step 3)",
                   "crs": srcs["dem"][1], "pixel_m": 30, "shape": list(srcs["dem"][3])},
        "caveats": ["Copernicus GLO-30 is a digital surface model: heights include buildings and vegetation.",
                    "The vertical datum is not stated in the DEM files and is not independently established here.",
                    "The analysis area was defined as DEM <= 1,800 m connected to Srinagar, so the analysis "
                    "is truncated at 1,800 m by construction.",
                    "Descriptive only: elevation and distance to the river are spatially related; no causal "
                    "claim is made."],
        "elevation_flooded": stats(dem[flooded]),
        "elevation_assessed": stats(dem[assessed]),
        "bands": rows,
        "diagnostics": diagnostics,
        "totals": {"analysis_area_pixels": int(area.sum()), "assessed_pixels": int(assessed.sum()),
                   "flooded_pixels": total_flooded, "flooded_km2": round(total_flooded * PIXEL_KM2, 4)},
    }
    (OUT / "elevation.json").write_text(json.dumps(result, indent=1))
    with open(OUT / "elevation_bands.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    # Figure: flood fraction and flooded area by band; assessed area for comparison.
    labels = ["<1570"] + [f"{a:.0f}" for a in EDGES[:-1]]
    x = np.arange(len(rows))
    frac = [r["flood_fraction_pct"] for r in rows]
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(11, 7.5), sharex=True,
                                 gridspec_kw={"height_ratios": [1.1, 1]})
    a1.bar(x, frac, color=["0.75" if r["small_band"] else "#1F5E59" for r in rows], width=0.8)
    a1.set_ylabel("flooded share of assessed area (%)")
    a1.set_title("Detected flood footprint, 10 September 2014, by DEM elevation band (descriptive)")
    a1.grid(axis="y", alpha=0.3)
    a2.bar(x - 0.2, [r["assessed_km2"] for r in rows], width=0.4, color="0.7", label="assessed area")
    a2.bar(x + 0.2, [r["flooded_km2"] for r in rows], width=0.4, color="#E69F00", label="flooded area")
    a2.set_ylabel("area (km2)")
    a2.set_xticks(x)
    a2.set_xticklabels(labels, rotation=60, ha="right", fontsize=8.5)
    a2.set_xlabel("band lower edge, m (Copernicus GLO-30 surface model; 10 m bands; last band 1790-1800)")
    a2.legend(frameon=False)
    a2.grid(axis="y", alpha=0.3)
    fig.text(0.01, 0.005, "Grey bars in the upper panel: bands with less than 10 km2 assessed (unstable "
             "fractions). Analysis area stops at 1,800 m by definition.", fontsize=8.5, color="0.3")
    fig.savefig(OUT / "figures" / "elevation.png", dpi=150, bbox_inches="tight")

    print(json.dumps({"flooded": result["elevation_flooded"], "assessed": result["elevation_assessed"],
                      "diagnostics": {k: v for k, v in diagnostics["flat_surfaces"].items() if k != "patches"},
                      "flat_patches": diagnostics["flat_surfaces"]["patches"][:6],
                      "high": diagnostics["flooded_at_or_above_1600m"]}, indent=1))
    for r in rows:
        print(f"{r['band']:28s} assessed {r['assessed_km2']:9.2f} km2  flooded {r['flooded_km2']:8.2f} km2  "
              f"{r['flood_fraction_pct']:6.2f} %  share {r['share_of_all_flooded_pct']:6.2f} %"
              f"{'  small' if r['small_band'] else ''}")


if __name__ == "__main__":
    main()
