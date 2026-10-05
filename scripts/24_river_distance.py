"""Step 6b: distance to the Jhelum centreline and the flood footprint by band.

Implements rules J5 (as amended, amendment 1) and D1-D8 of scripts/step6_rules.md. Distance is the
Euclidean distance in metres (EPSG:32643 plane) from each pixel centre to the
centreline selected by 23_jhelum.py, from an exact Euclidean distance
transform on the 30 m grid padded by 20 km, after rasterising the centreline
(all touched pixels); the error is at most half a pixel diagonal (21.2 m).
It is not hydrological distance. Pixels inside Wular get their distance to
the mapped centreline, not to the lake shore.

Descriptive only: no regression and no causal claim. Distance and elevation
are spatially related.

Outputs: data/river/jhelum_distance_m.tif (intermediate, not committed),
outputs/river_distance_bands.csv, outputs/river_distance.json,
outputs/figures/river_distance.png

Usage:  .venv/bin/python scripts/24_river_distance.py
"""

import csv
import json
import math
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
import shapely
from rasterio.features import rasterize
from rasterio.transform import Affine
from rasterio.warp import transform_geom
from scipy import ndimage
from shapely.geometry import shape

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "outputs"
RIVER_DIR = ROOT / "data" / "river"
CRS = "EPSG:32643"
PIXEL_KM2 = 0.0009
PAD_PX = 667                                                    # D1: 20 km padding
NODATA = -9999.0                                                # D5
EDGES = [0, 250, 500, 1000, 2000, 3000, 5000, 7500, 10000, 15000]   # D2, metres
PCTS = [1, 5, 10, 25, 75, 90, 95, 99]                            # D4
FLOODED, PERSISTENT, NA = 1, 2, 255


def stop(msg):
    sys.exit(f"STOP: {msg}")


def stats(v):
    v = v.astype(np.float64)
    out = {"pixels": int(v.size), "min_m": round(float(v.min()), 1), "max_m": round(float(v.max()), 1),
           "mean_m": round(float(v.mean()), 1), "median_m": round(float(np.median(v)), 1),
           "std_m": round(float(v.std()), 1)}
    out.update({f"p{p}_m": round(float(x), 1) for p, x in zip(PCTS, np.percentile(v, PCTS))})
    return out


def bands_from_rules():
    """The D2 band edges as written in the committed rules file."""
    text = (Path(__file__).parent / "step6_rules.md").read_text()
    m = re.search(r"- D2\. Bands \(lower bound inclusive, upper exclusive\):(.*?)Reason:", text, re.S)
    if not m:
        stop("D2 not found in step6_rules.md")
    nums = [int(n.replace(",", "")) for n in re.findall(r"\d[\d,]*", m.group(1))]
    edges = sorted(set(nums))
    return edges


def edt_to(geom, shape_, transform):
    """Distance (m) from every pixel centre to the pixels touched by geom."""
    burn = rasterize([(geom, 1)], out_shape=shape_, transform=transform, all_touched=True, dtype="uint8")
    if not burn.any():
        stop("geometry does not touch the padded grid")
    return ndimage.distance_transform_edt(burn == 0) * transform.a


def main():
    (OUT / "figures").mkdir(parents=True, exist_ok=True)
    RIVER_DIR.mkdir(parents=True, exist_ok=True)
    fr = json.loads((ROOT / "data" / "flood_report.json").read_text())
    jh = json.loads((OUT / "jhelum.json").read_text())
    gj_path = OUT / "jhelum_centreline.geojson"
    if not gj_path.exists():
        stop("outputs/jhelum_centreline.geojson is missing; run 23_jhelum.py")
    gj = json.loads(gj_path.read_text())

    with rasterio.open(ROOT / "data" / "flood" / "flood_classes.tif") as src:
        cls, T, H, W, crs = src.read(1), src.transform, src.height, src.width, src.crs.to_string()
        profile = src.profile
    with rasterio.open(ROOT / "data" / "flood" / "analysis_area.tif") as src:
        if (src.transform, src.shape, src.crs.to_string()) != (T, (H, W), crs):
            stop("analysis area and flood classes are not on one grid")
        area = src.read(1).astype(bool)
    if crs != CRS:
        stop(f"flood grid CRS {crs} is not {CRS}")
    if round(int(area.sum()) * PIXEL_KM2, 4) != fr["analysis_area_km2"]:
        stop("analysis area differs from step 3")

    lines = [shape(transform_geom("EPSG:4326", CRS, f["geometry"])) for f in gj["features"]]
    river = shapely.union_all(lines)

    # D1: padded grid, exact EDT to the rasterised centreline.
    Tp = Affine(T.a, 0, T.c - PAD_PX * T.a, 0, T.e, T.f - PAD_PX * T.e)
    shp_p = (H + 2 * PAD_PX, W + 2 * PAD_PX)
    d_river = edt_to(river, shp_p, Tp)[PAD_PX:PAD_PX + H, PAD_PX:PAD_PX + W]

    if bands_from_rules() != EDGES:
        stop(f"band edges in the script {EDGES} differ from step6_rules.md D2 {bands_from_rules()}")

    # J5 (amended): pixels closer to the extract boundary than to the mapped
    # line are boundary-affected; their true distance lies in [d_bound,
    # d_river] and they may proceed only if that range stays in one band.
    index = json.loads((ROOT / "data" / "osm" / "geofabrik" / "index-v1.json").read_text())
    ext_id = json.loads((ROOT / "data" / "osm" / "geofabrik" / "extract.json").read_text())["extract_id"]
    feat = [f for f in index["features"] if f["properties"]["id"] == ext_id]
    if len(feat) != 1:
        stop(f"extract {ext_id} not found once in the saved Geofabrik index")
    boundary = shape(transform_geom("EPSG:4326", CRS, feat[0]["geometry"])).boundary
    d_bound = edt_to(boundary, shp_p, Tp)[PAD_PX:PAD_PX + H, PAD_PX:PAD_PX + W]
    affected = area & (d_river > d_bound)
    lo_band = np.searchsorted(EDGES, d_bound[affected], side="right")
    hi_band = np.searchsorted(EDGES, d_river[affected], side="right")
    crossing = int(np.count_nonzero(lo_band != hi_band))
    j5 = {"rule": "J5 as amended (amendment 1)",
          "boundary_affected_pixels": int(affected.sum()),
          "boundary_affected_km2": round(int(affected.sum()) * PIXEL_KM2, 4),
          "pixels_whose_range_crosses_a_band_edge": crossing}
    if affected.any():
        j5.update({"lower_bound_m_range": [round(float(d_bound[affected].min()), 1),
                                           round(float(d_bound[affected].max()), 1)],
                   "upper_bound_m_range": [round(float(d_river[affected].min()), 1),
                                           round(float(d_river[affected].max()), 1)],
                   "bands_of_affected_pixels": sorted({EDGES[i - 1] for i in hi_band.tolist()})})
    if crossing:
        stop(f"J5 (amended): {crossing} boundary-affected pixels have a distance range that crosses a band edge")

    # D5 / D6: raster and checks.
    dist = np.where(area, d_river, NODATA).astype(np.float32)
    diag = math.hypot(*shp_p) * T.a
    if (dist[area] < 0).any() or dist[area].max() > diag:
        stop("distance values outside the physically possible range")
    out_tif = RIVER_DIR / "jhelum_distance_m.tif"
    prof = {**profile, "dtype": "float32", "nodata": NODATA, "compress": "deflate", "tiled": True,
            "blockxsize": 256, "blockysize": 256}
    with rasterio.open(out_tif, "w", **prof) as dst:
        dst.write(dist, 1)
        dst.update_tags(UNITS="metres", CRS=CRS,
                        DEFINITION="Euclidean distance in the EPSG:32643 plane from the pixel centre to the "
                                   "Jhelum centreline rasterised on this grid (all touched), exact EDT, "
                                   "max error 21.2 m; not hydrological distance",
                        SOURCE="outputs/jhelum_centreline.geojson: OSM relation 8306691 main_stream "
                               f"members, {jh['source']['file']}",
                        OSM_TIMESTAMP=jh["source"]["osm_timestamp"],
                        NODATA="-9999 outside the analysis area")
    flag_tif = RIVER_DIR / "jhelum_distance_flag.tif"
    with rasterio.open(flag_tif, "w", **{**prof, "dtype": "uint8", "nodata": None}) as dst:
        dst.write(affected.astype(np.uint8), 1)
        dst.update_tags(VALUES="1 = boundary-affected (J5 amended): the distance raster value is an upper "
                               "bound; the true distance is at least the distance to the extract boundary")
    for path in (out_tif, flag_tif):
        with rasterio.open(path) as chk:
            if (chk.transform, chk.shape, chk.crs.to_string()) != (T, (H, W), CRS):
                stop(f"written {path.name} does not match the flood grid")

    # D3 / D4: bands and statistics.
    assessed = area & np.isin(cls, (0, 1, 2, 3))
    flooded = area & (cls == FLOODED)
    persistent = area & (cls == PERSISTENT)
    not_assessed = area & (cls == NA)
    total_flooded = int(flooded.sum())
    if round(total_flooded * PIXEL_KM2, 4) != fr["primary"]["flooded_km2"]:
        stop("flooded pixels do not reproduce the step 3 result")
    rows = []
    bounds = EDGES + [None]
    for lo, hi in zip(bounds[:-1], bounds[1:]):
        inb = area & (d_river >= lo) & ((d_river < hi) if hi is not None else True)
        n_ass, n_fl = int((inb & assessed).sum()), int((inb & flooded).sum())
        rows.append({"band": f"{lo:,}-{hi:,} m" if hi is not None else f">= {lo:,} m",
                     "lower_m": lo, "upper_m": hi,
                     "analysis_area_pixels": int(inb.sum()), "assessed_pixels": n_ass,
                     "assessed_km2": round(n_ass * PIXEL_KM2, 4),
                     "not_assessed_km2": round(int((inb & not_assessed).sum()) * PIXEL_KM2, 4),
                     "flooded_pixels": n_fl, "flooded_km2": round(n_fl * PIXEL_KM2, 4),
                     "flood_fraction_pct": round(100 * n_fl / n_ass, 3) if n_ass else None,
                     "share_of_all_flooded_pct": round(100 * n_fl / total_flooded, 3),
                     "persistent_water_km2": round(int((inb & persistent).sum()) * PIXEL_KM2, 4),
                     "boundary_affected_km2": round(int((inb & affected).sum()) * PIXEL_KM2, 4)})
    if sum(r["analysis_area_pixels"] for r in rows) != int(area.sum()):
        stop("bands do not cover the analysis area exactly once")
    if sum(r["flooded_pixels"] for r in rows) != total_flooded:
        stop("bands do not account for every flooded pixel")

    result = {
        "script": "scripts/24_river_distance.py", "rules": "scripts/step6_rules.md (J5, D)",
        "band_edges_m": EDGES,
        "distance_raster": {"file": "data/river/jhelum_distance_m.tif",
                            "flag_raster": "data/river/jhelum_distance_flag.tif", "crs": CRS, "units": "metres",
                            "pixel_m": T.a, "shape": [H, W], "nodata": NODATA,
                            "padding_m": PAD_PX * T.a, "max_error_m": round(T.a * math.sqrt(2) / 2, 1),
                            "source": "outputs/jhelum_centreline.geojson",
                            "osm_timestamp": jh["source"]["osm_timestamp"]},
        "j5_missing_members_check": j5,
        "boundary_affected": {
            "note": "Exact river distances of these pixels are uncertain: the extract does not cover all "
                    "downstream river geometry. They are assigned to a band only because their whole possible "
                    "range lies in it, and they are excluded from the exact-distance statistics.",
            "assessed_km2": round(int((affected & assessed).sum()) * PIXEL_KM2, 4),
            "flooded_km2": round(int((affected & flooded).sum()) * PIXEL_KM2, 4),
            "not_assessed_km2": round(int((affected & not_assessed).sum()) * PIXEL_KM2, 4),
            "distance_lower_bound_m": {"min": j5.get("lower_bound_m_range", [None])[0],
                                       "max": j5.get("lower_bound_m_range", [None, None])[1]},
            "distance_upper_bound_m": {"min": j5.get("upper_bound_m_range", [None])[0],
                                       "max": j5.get("upper_bound_m_range", [None, None])[1]},
        },
        "caveats": ["Euclidean distance to the mapped centreline, not hydrological distance.",
                    "Pixels inside Wular Lake get their distance to the centreline through the lake, "
                    "not to the lake shore.",
                    "The centreline is today's OpenStreetMap (snapshot above), not the 2014 channel.",
                    "Descriptive only: distance and elevation are spatially related; no causal claim "
                    "and no regression."],
        "distance_flooded_exact_only": stats(d_river[flooded & ~affected]),
        "distance_assessed_exact_only": stats(d_river[assessed & ~affected]),
        "exact_only_note": "Statistics exclude boundary-affected pixels, whose exact distance is not known.",
        "bands": rows,
        "totals": {"analysis_area_pixels": int(area.sum()), "assessed_pixels": int(assessed.sum()),
                   "flooded_pixels": total_flooded, "flooded_km2": round(total_flooded * PIXEL_KM2, 4)},
    }
    (OUT / "river_distance.json").write_text(json.dumps(result, indent=1))
    with open(OUT / "river_distance_bands.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    labels = [r["band"].replace(" m", "") for r in rows]
    x = np.arange(len(rows))
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(10.5, 7.5), sharex=True,
                                 gridspec_kw={"height_ratios": [1.1, 1]})
    a1.bar(x, [r["flood_fraction_pct"] for r in rows], color="#1F5E59", width=0.75)
    a1.set_ylabel("flooded share of assessed area (%)")
    a1.set_title("Detected flood footprint, 10 September 2014, by distance to the mapped Jhelum "
                 "centreline (descriptive)")
    a1.grid(axis="y", alpha=0.3)
    a2.bar(x - 0.2, [r["assessed_km2"] for r in rows], width=0.4, color="0.7", label="assessed area")
    a2.bar(x + 0.2, [r["flooded_km2"] for r in rows], width=0.4, color="#E69F00", label="flooded area")
    a2.set_ylabel("area (km2)")
    a2.set_xticks(x)
    a2.set_xticklabels(labels, rotation=35, ha="right")
    a2.set_xlabel("Euclidean distance to the centreline (m); bands of unequal width; not hydrological distance")
    a2.legend(frameon=False)
    a2.grid(axis="y", alpha=0.3)
    fig.savefig(OUT / "figures" / "river_distance.png", dpi=150, bbox_inches="tight")

    print(json.dumps({"j5": j5, "boundary_affected": result["boundary_affected"],
                      "flooded_exact_only": result["distance_flooded_exact_only"],
                      "assessed_exact_only": result["distance_assessed_exact_only"]}, indent=1))
    for r in rows:
        print(f"{r['band']:18s} assessed {r['assessed_km2']:8.2f}  flooded {r['flooded_km2']:7.2f} km2  "
              f"{r['flood_fraction_pct']:6.2f} %  share {r['share_of_all_flooded_pct']:6.2f} %  "
              f"persistent {r['persistent_water_km2']:6.2f}  n/a {r['not_assessed_km2']:.2f}")


if __name__ == "__main__":
    main()
