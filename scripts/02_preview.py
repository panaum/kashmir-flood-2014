"""Step 1c: cloud check and false-colour previews of the mosaicked data.

Reads the mosaics made by 01b_mosaic.py and:
  - counts cloud, cloud-shadow, cirrus and dilated-cloud pixels from QA_PIXEL,
    over the Srinagar box and over the whole valid area, for both dates;
  - renders SWIR1 / NIR / Red as R / G / B: Srinagar for each date, a
    Srinagar side-by-side, and a whole-valley side-by-side.

All previews use one shared colour stretch, so a colour means the same thing
in every image. The stretch is for viewing only; no data files are modified.
Numbers go to previews/qa_summary.json.

Usage:  .venv/bin/python scripts/02_preview.py
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from rasterio.warp import transform_bounds
from rasterio.windows import from_bounds

ROOT = Path(__file__).resolve().parent.parent
CLIP_DIR = ROOT / "data" / "clipped"
OUT_DIR = ROOT / "previews"
SCENES = {"pre": "2014-08-25", "flood": "2014-09-10"}
RGB = ("swir1", "nir", "red")

# Srinagar city and surroundings, WGS84 lon/lat (west, south, east, north).
# Centre of Srinagar is roughly 74.80 E, 34.08 N.
SRINAGAR_BBOX = (74.65, 33.95, 74.95, 34.20)

# Landsat Collection 2 QA_PIXEL bit positions, as listed in the qa_pixel
# asset's classification:bitfields in the STAC item (see DATA.md).
QA_BITS = {"fill": 0, "dilated_cloud": 1, "cirrus": 2, "cloud": 3,
           "cloud_shadow": 4, "snow": 5, "clear": 6}


def read(label, band):
    """Whole mosaic for a band, plus the Srinagar window into it."""
    path = CLIP_DIR / f"{label}_{SCENES[label]}_{band}.tif"
    with rasterio.open(path) as src:
        b = transform_bounds("EPSG:4326", src.crs, *SRINAGAR_BBOX, densify_pts=21)
        win = from_bounds(*b, transform=src.transform)
        win = win.round_offsets(op="floor").round_lengths(op="ceil")
        return src.read(1), win.toslices(), src.nodata


def qa_stats(qa):
    valid = ((qa >> QA_BITS["fill"]) & 1) == 0
    q = qa[valid]
    out = {"valid_pixels": int(q.size)}
    for name, bit in QA_BITS.items():
        if name != "fill":
            out[f"{name}_fraction"] = round(float(np.mean((q >> bit) & 1)), 4)
    out["cloud_or_shadow_fraction"] = round(
        float(np.mean(((q >> 3) | (q >> 4)) & 1)), 4)
    out["cloud_shadow_cirrus_or_dilated_fraction"] = round(
        float(np.mean(((q >> 1) | (q >> 2) | (q >> 3) | (q >> 4)) & 1)), 4)
    return out


def to_display(rgb, valid, lims):
    out = np.zeros(rgb.shape, dtype=np.float32)
    for c, (lo, hi) in enumerate(lims):
        out[..., c] = np.clip((rgb[..., c] - lo) / (hi - lo), 0, 1)
    out[~valid] = 1.0  # no data shown as white
    return out


def main():
    OUT_DIR.mkdir(exist_ok=True)
    data = {}
    for label in SCENES:
        bands = {}
        for band in RGB + ("qa_pixel",):
            arr, srinagar, nodata = read(label, band)
            bands[band] = arr
        valid = ((bands["qa_pixel"] & 1) == 0)
        rgb = np.stack([bands[b].astype(np.float32) for b in RGB], axis=-1)
        data[label] = {"rgb": rgb, "valid": valid, "qa": bands["qa_pixel"],
                       "srinagar": srinagar}

    # One stretch for every preview: 2nd-98th percentile per channel of all
    # valid valley pixels, both dates pooled.
    lims = []
    for c in range(3):
        vals = np.concatenate([d["rgb"][..., c][d["valid"]] for d in data.values()])
        lims.append(np.percentile(vals, [2, 98]))

    summary = {"srinagar_bbox_wgs84": SRINAGAR_BBOX,
               "stretch_2_98_percentile_raw_DN": {
                   b: [float(x) for x in l] for b, l in zip(RGB, lims)},
               "qa_pixel": {}}
    valley, srinagar = {}, {}
    for label, date in SCENES.items():
        d = data[label]
        sl = d["srinagar"]
        summary["qa_pixel"][date] = {
            "srinagar": qa_stats(d["qa"][sl]),
            "srinagar_pixels_total": int(d["qa"][sl].size),
            "whole_valid_area": qa_stats(d["qa"]),
        }
        valley[label] = to_display(d["rgb"], d["valid"], lims)
        srinagar[label] = valley[label][sl]

        fig, ax = plt.subplots(figsize=(7, 7))
        ax.imshow(srinagar[label], interpolation="nearest")
        ax.set_title(f"Srinagar, {date}\nfalse colour R=SWIR1, G=NIR, B=Red")
        ax.axis("off")
        fig.savefig(OUT_DIR / f"srinagar_{label}_{date}.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

    for name, imgs, size in (("srinagar", srinagar, (14, 7.5)),
                             ("valley", valley, (14, 8.5))):
        fig, axes = plt.subplots(1, 2, figsize=size)
        for ax, (label, date) in zip(axes, SCENES.items()):
            ax.imshow(imgs[label], interpolation="antialiased")
            ax.set_title(f"{date} ({'before flood' if label == 'pre' else 'during flood'})")
            ax.axis("off")
        where = "Srinagar" if name == "srinagar" else "Kashmir valley (rows 036+037 mosaic)"
        fig.suptitle(f"{where}, Landsat 8 false colour (R=SWIR1, G=NIR, B=Red), "
                     "same stretch; white = no data")
        fig.savefig(OUT_DIR / f"{name}_side_by_side.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

    (OUT_DIR / "qa_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
