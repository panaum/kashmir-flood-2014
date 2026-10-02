"""Step 4b-c: Microsoft footprints on the Srinagar flood map.

Base: the step 3 flood map (2014-09-10 false colour, flooded pixels red,
persistent water blue), as in 14_exposure_preview.py. Left: the Srinagar box
with footprint centroids (yellow = on flooded ground, white = not). Right: a
6 km x 6 km zoom on the city centre with the footprint outlines themselves
and the flooded motor-road pieces from step 4 (black).

Usage:  .venv/bin/python scripts/17_msbuildings_preview.py
"""

import csv
import importlib.util
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyarrow.parquet as pq
import rasterio
import shapely
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D
from rasterio.warp import transform as warp_transform
from rasterio.warp import transform_bounds
from rasterio.windows import from_bounds

ROOT = Path(__file__).resolve().parent.parent
SRINAGAR_BBOX = (74.65, 33.95, 74.95, 34.20)
ZOOM_CENTRE, ZOOM_HALF_M = (34.07, 74.81), 3000
FLOOD_RGB, PERSIST_RGB = (1.0, 0.15, 0.1), (0.1, 0.35, 1.0)

_spec = importlib.util.spec_from_file_location("ex", Path(__file__).parent / "13_exposure.py")
ex = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ex)


def main():
    lims = list(json.loads((ROOT / "previews" / "qa_summary.json").read_text())
                ["stretch_2_98_percentile_raw_DN"].values())
    rep = json.loads((ROOT / "data" / "msbuildings_report.json").read_text())
    with rasterio.open(ROOT / "data" / "flood" / "flood_classes.tif") as src:
        cls = src.read(1)
        T, crs = src.transform, src.crs
    area = rasterio.open(ROOT / "data" / "flood" / "analysis_area.tif").read(1).astype(bool)
    rgb = np.stack([rasterio.open(ROOT / "data" / "clipped" / f"flood_2014-09-10_{b}.tif")
                    .read(1).astype(np.float32) for b in ("swir1", "nir", "red")], axis=-1)
    for c, (lo, hi) in enumerate(lims):
        rgb[..., c] = np.clip((rgb[..., c] - lo) / (hi - lo), 0, 1)
    rgb[~area] *= 0.4
    rgb[cls == 1] = 0.55 * np.array(FLOOD_RGB) + 0.45 * rgb[cls == 1]
    rgb[cls == 2] = PERSIST_RGB
    sb = transform_bounds("EPSG:4326", crs, *SRINAGAR_BBOX, densify_pts=21)
    win = from_bounds(*sb, transform=T).round_offsets(op="floor").round_lengths(op="ceil")
    img = rgb[win.toslices()]
    x0, y1 = T * (int(win.col_off), int(win.row_off))
    extent = (x0, x0 + img.shape[1] * 30, y1 - img.shape[0] * 30, y1)

    pts = {"flooded": [], "dry": []}
    with open(ROOT / "data" / "msbuildings" / "footprints.csv") as f:
        for row in csv.DictReader(f):
            if row["in_srinagar_box"] == "1" and row["state_centroid"] in pts:
                pts[row["state_centroid"]].append((float(row["centroid_x"]), float(row["centroid_y"])))

    # Footprint outlines for the zoom window.
    zx, zy = warp_transform("EPSG:4326", crs, [ZOOM_CENTRE[1]], [ZOOM_CENTRE[0]])
    zoom = (zx[0] - ZOOM_HALF_M, zx[0] + ZOOM_HALF_M, zy[0] - ZOOM_HALF_M, zy[0] + ZOOM_HALF_M)
    zx_ll, zy_ll = warp_transform(crs, "EPSG:4326", [zoom[0], zoom[1]], [zoom[2], zoom[3]])
    record = json.loads((ROOT / "data" / "msbuildings" / "record.json").read_text())
    geoms = []
    for it in record["items"]:
        b = it["bbox"]
        if b[0] > zx_ll[1] or b[2] < zx_ll[0] or b[1] > zy_ll[1] or b[3] < zy_ll[0]:
            continue
        for f in it["files"]:
            g = shapely.from_wkb(pq.read_table(ROOT / f["local"], columns=["geometry"])
                                 .column("geometry").to_pylist())
            bb = shapely.bounds(g)
            keep = (bb[:, 2] >= zx_ll[0]) & (bb[:, 0] <= zx_ll[1]) & (bb[:, 3] >= zy_ll[0]) & (bb[:, 1] <= zy_ll[1])
            geoms.extend(g[keep])
    utm = ex.to_utm(geoms, crs)
    state, code, *_ = ex.state_rasters()
    _, _, _, _, b_state, _ = ex.classify_buildings(utm, code, state, T, state.shape)
    polys = {"flooded": [], "dry": []}
    for g, s in zip(utm, b_state):
        key = "flooded" if s == ex.FLOODED else "dry" if s == ex.DRY else None
        if key:
            for part in getattr(g, "geoms", [g]):
                polys[key].append(np.asarray(part.exterior.coords))

    pieces = json.loads((ROOT / "data" / "exposure" / "flooded_road_pieces.json").read_text())
    motor = []
    for gj in pieces["motor"]:
        g = shapely.from_geojson(gj)
        for part in getattr(g, "geoms", [g]):
            if part.geom_type == "LineString":
                motor.append(np.asarray(part.coords))

    c = rep["counts"]["centroid"]["srinagar_box"]
    fig, axes = plt.subplots(1, 2, figsize=(16, 8.6))
    axes[0].imshow(img, extent=extent, interpolation="nearest")
    for key, colour, s, a in (("dry", "white", 0.15, 0.5), ("flooded", "yellow", 0.3, 0.9)):
        xy = np.asarray(pts[key])
        if len(xy):
            axes[0].scatter(xy[:, 0], xy[:, 1], s=s, c=colour, alpha=a, linewidths=0, zorder=3)
    axes[0].set_title("Srinagar box: footprint centroids")
    axes[1].imshow(img, extent=extent, interpolation="nearest")
    axes[1].add_collection(PolyCollection(polys["dry"], facecolors="0.85", edgecolors="none", zorder=3))
    axes[1].add_collection(PolyCollection(polys["flooded"], facecolors="yellow", edgecolors="black",
                                          linewidths=0.15, zorder=4))
    for seg in motor:
        axes[1].plot(seg[:, 0], seg[:, 1], color="black", lw=0.8, zorder=5)
    axes[1].set_xlim(zoom[0], zoom[1])
    axes[1].set_ylim(zoom[2], zoom[3])
    axes[1].set_title("City centre, 6 km x 6 km: footprints and flooded motor roads")
    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle(f"Microsoft building footprints (imagery 2017-2021) standing on ground flooded on "
                 f"2014-09-10\nSrinagar box: {c['flooded']} of {c['total_in_analysis_area']} footprints "
                 "(centroid rule). Not a count of buildings flooded in 2014.")
    fig.legend(handles=[Line2D([], [], marker="s", ls="", mfc="yellow", mec="black", label="footprint on flooded ground"),
                        Line2D([], [], marker="s", ls="", mfc="0.85", mec="none", label="footprint, not flooded"),
                        Line2D([], [], color="black", label="flooded motor road (OSM)"),
                        Line2D([], [], marker="s", ls="", mfc=FLOOD_RGB, mec="none", label="flooded (step 3)"),
                        Line2D([], [], marker="s", ls="", mfc=PERSIST_RGB, mec="none", label="persistent water")],
               loc="lower center", ncol=5, frameon=False)
    fig.savefig(ROOT / "previews" / "msbuildings_srinagar.png", dpi=150, bbox_inches="tight")
    print(f"wrote previews/msbuildings_srinagar.png ({len(polys['flooded'])} flooded, "
          f"{len(polys['dry'])} dry outlines in zoom)")


if __name__ == "__main__":
    main()
