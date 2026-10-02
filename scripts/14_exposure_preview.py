"""Step 4c: flooded buildings and roads over the Srinagar flood map.

Base: the step 3 flood map (2014-09-10 false colour with flooded pixels in
red and persistent water in blue). On top: OSM building centroids (yellow =
flooded, white = assessed and not flooded) and the flooded parts of roads
(black = motor roads, grey = other highway types). Left: the Srinagar box;
right: a 6 km x 6 km zoom on the city centre (around 34.07 N, 74.81 E).

Usage:  .venv/bin/python scripts/14_exposure_preview.py
"""

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
import shapely
from matplotlib.lines import Line2D
from rasterio.warp import transform as warp_transform
from rasterio.warp import transform_bounds
from rasterio.windows import from_bounds

ROOT = Path(__file__).resolve().parent.parent
SRINAGAR_BBOX = (74.65, 33.95, 74.95, 34.20)
ZOOM_CENTRE, ZOOM_HALF_M = (34.07, 74.81), 3000
FLOOD_RGB, PERSIST_RGB = (1.0, 0.15, 0.1), (0.1, 0.35, 1.0)


def main():
    lims = list(json.loads((ROOT / "previews" / "qa_summary.json").read_text())
                ["stretch_2_98_percentile_raw_DN"].values())
    rep = json.loads((ROOT / "data" / "exposure_report.json").read_text())
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
    r0, c0 = int(win.row_off), int(win.col_off)
    img = rgb[win.toslices()]
    x0, y1 = T * (c0, r0)
    extent = (x0, x0 + img.shape[1] * 30, y1 - img.shape[0] * 30, y1)

    pts = {"flooded": [], "dry": []}
    with open(ROOT / "data" / "exposure" / "buildings.csv") as f:
        for row in csv.DictReader(f):
            if row["in_srinagar_box"] == "1" and row["state_centroid"] in pts:
                pts[row["state_centroid"]].append((float(row["centroid_x"]), float(row["centroid_y"])))
    pieces = json.loads((ROOT / "data" / "exposure" / "flooded_road_pieces.json").read_text())

    def lines(key):
        out = []
        for gj in pieces[key]:
            g = shapely.from_geojson(gj)
            for part in getattr(g, "geoms", [g]):
                if part.geom_type == "LineString":
                    out.append(np.asarray(part.coords))
        return out

    motor, other = lines("motor"), lines("other")
    zx, zy = warp_transform("EPSG:4326", crs, [ZOOM_CENTRE[1]], [ZOOM_CENTRE[0]])
    zoom = (zx[0] - ZOOM_HALF_M, zx[0] + ZOOM_HALF_M, zy[0] - ZOOM_HALF_M, zy[0] + ZOOM_HALF_M)

    b = rep["buildings"]["centroid"]["srinagar_box"]
    m = rep["roads"]["motor"]["srinagar_box"]
    fig, axes = plt.subplots(1, 2, figsize=(16, 8.6))
    for ax, ext, size in ((axes[0], None, (4, 1.5)), (axes[1], zoom, (16, 6))):
        ax.imshow(img, extent=extent, interpolation="nearest")
        for seg in other:
            ax.plot(seg[:, 0], seg[:, 1], color="0.55", lw=0.5)
        for seg in motor:
            ax.plot(seg[:, 0], seg[:, 1], color="black", lw=0.6 if ext is None else 1.2)
        for state, colour, s in (("dry", "white", size[1]), ("flooded", "yellow", size[0])):
            if pts[state]:
                xy = np.asarray(pts[state])
                ax.scatter(xy[:, 0], xy[:, 1], s=s, c=colour, edgecolors="black" if state == "flooded" else "none",
                           linewidths=0.3, zorder=3)
        if ext:
            ax.set_xlim(ext[0], ext[1])
            ax.set_ylim(ext[2], ext[3])
        ax.set_xticks([])
        ax.set_yticks([])
    axes[0].set_title("Srinagar box")
    axes[1].set_title("City centre, 6 km x 6 km")
    fig.suptitle(f"Srinagar box: {b['flooded']} of {b['total_in_analysis_area']} mapped OSM buildings "
                 f"flooded (centroid rule); {m['flooded']:.1f} of {m['total_in_analysis_area']:.1f} km "
                 "motor road flooded.\nOSM map of today, not 2014; most buildings are not mapped in OSM.")
    fig.legend(handles=[Line2D([], [], marker="o", ls="", mfc="yellow", mec="black", label="flooded building"),
                        Line2D([], [], marker="o", ls="", mfc="white", mec="0.5", label="building, not flooded"),
                        Line2D([], [], color="black", label="flooded motor road"),
                        Line2D([], [], color="0.55", label="flooded other way"),
                        Line2D([], [], marker="s", ls="", mfc=FLOOD_RGB, mec="none", label="flooded (step 3)"),
                        Line2D([], [], marker="s", ls="", mfc=PERSIST_RGB, mec="none", label="persistent water")],
               loc="lower center", ncol=6, frameon=False)
    fig.savefig(ROOT / "previews" / "exposure_srinagar.png", dpi=150, bbox_inches="tight")
    print("wrote previews/exposure_srinagar.png")


if __name__ == "__main__":
    main()
