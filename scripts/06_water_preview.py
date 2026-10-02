"""Step 2c: water masks drawn over the false-colour images.

For both dates, draws the water mask from 05_water.py over the SWIR1 / NIR /
Red false-colour image (same stretch as 02_preview.py), for Srinagar and for
the whole valley, as two side-by-side figures. Also draws the three
sanity-check boxes on both dates so the site numbers can be checked by eye.

Colours: cyan = water; magenta = masked (cloud, shadow, dilated cloud, snow,
reflectance <= 0, no data); darkened = above the 1800 m valley floor (outside
the analysis area; not counted).

Usage:  .venv/bin/python scripts/06_water_preview.py
"""

import importlib.util
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from matplotlib.patches import Patch, Rectangle
from rasterio.warp import transform_bounds
from rasterio.windows import from_bounds

ROOT = Path(__file__).resolve().parent.parent
CLIP_DIR = ROOT / "data" / "clipped"
WATER_DIR = ROOT / "data" / "water"
OUT_DIR = ROOT / "previews"
DATES = {"pre": "2014-08-25", "flood": "2014-09-10"}
RGB = ("swir1", "nir", "red")
WATER_RGB = (0.0, 0.9, 1.0)
MASKED_RGB = (1.0, 0.0, 1.0)
SRINAGAR_BBOX = (74.65, 33.95, 74.95, 34.20)  # as in 02_preview.py

_spec = importlib.util.spec_from_file_location("wa", Path(__file__).parent / "05_water.py")
wa = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(wa)


def composite(label, date, lims, floor):
    rgb = np.stack([rasterio.open(CLIP_DIR / f"{label}_{date}_{b}.tif").read(1).astype(np.float32)
                    for b in RGB], axis=-1)
    for c, (lo, hi) in enumerate(lims):
        rgb[..., c] = np.clip((rgb[..., c] - lo) / (hi - lo), 0, 1)
    water = rasterio.open(WATER_DIR / f"{label}_{date}_water_mask.tif").read(1)
    rgb[~floor] *= 0.45
    on = floor & (water == 1)
    rgb[on] = WATER_RGB
    masked = floor & (water == 255)
    rgb[masked] = 0.4 * rgb[masked] + 0.6 * np.array(MASKED_RGB)
    return rgb


def legend(fig):
    fig.legend(handles=[Patch(color=WATER_RGB, label="water"),
                        Patch(color=MASKED_RGB, label="masked (cloud, shadow, snow, no data)"),
                        Patch(color="0.25", label="above 1800 m (not analysed)")],
               loc="lower center", ncol=3, frameon=False)


def main():
    lims = list(json.loads((OUT_DIR / "qa_summary.json").read_text())
                ["stretch_2_98_percentile_raw_DN"].values())
    report = json.loads((ROOT / "data" / "water_report.json").read_text())
    with rasterio.open(WATER_DIR / "valley_floor.tif") as src:
        floor = src.read(1).astype(bool)
        transform, crs = src.transform, src.crs
    imgs = {label: composite(label, date, lims, floor) for label, date in DATES.items()}

    sb = transform_bounds("EPSG:4326", crs, *SRINAGAR_BBOX, densify_pts=21)
    sr = from_bounds(*sb, transform=transform).round_offsets(op="floor") \
        .round_lengths(op="ceil").toslices()  # same window code as 02_preview.py
    for name, cut, size in (("srinagar", sr, (14, 7.8)), ("valley", (slice(None), slice(None)), (14, 8.8))):
        fig, axes = plt.subplots(1, 2, figsize=size)
        for ax, (label, date) in zip(axes, DATES.items()):
            t = report["dates"][date]["thresholds"]["otsu"]
            km = report["dates"][date]["water_km2_floor"]["otsu"]
            ax.imshow(imgs[label][cut], interpolation="nearest" if name == "srinagar" else "antialiased")
            ax.set_title(f"{date}: MNDWI > {t:.3f} (Otsu)\nvalley-floor water {km:.1f} km2")
            ax.axis("off")
        where = "Srinagar" if name == "srinagar" else "Kashmir valley"
        fig.suptitle(f"{where}: water mask over false colour (R=SWIR1, G=NIR, B=Red)")
        legend(fig)
        fig.savefig(OUT_DIR / f"water_{name}_side_by_side.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

    # Sanity-check sites: each box drawn on a 0.1 degree crop, both dates.
    fig, axes = plt.subplots(2, len(wa.SITES), figsize=(5 * len(wa.SITES), 10.5))
    for j, (label, date) in enumerate(DATES.items()):
        for i, (site, (lat, lon)) in enumerate(wa.SITES.items()):
            rs, cs = wa.lonlat_window(transform, crs, lat, lon, wa.SITE_HALF_DEG)
            big = wa.lonlat_window(transform, crs, lat, lon, 0.05)
            ax = axes[j, i]
            ax.imshow(imgs[label][big], interpolation="nearest")
            ax.add_patch(Rectangle((cs.start - big[1].start - 0.5, rs.start - big[0].start - 0.5),
                                   cs.stop - cs.start, rs.stop - rs.start,
                                   fill=False, ec="yellow", lw=1.8))
            s = report["dates"][date]["sites"][site]
            ax.set_title(f"{date} {site}\nwater {100 * s['water_fraction_of_valid']:.1f} % "
                         f"of valid, masked {100 * s['masked_fraction']:.1f} %")
            ax.axis("off")
    fig.suptitle(f"Sanity-check boxes (yellow, +/- {wa.SITE_HALF_DEG} deg) on 0.1 deg crops")
    legend(fig)
    fig.savefig(OUT_DIR / "water_sanity_sites.png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    print("wrote previews/water_srinagar_side_by_side.png, water_valley_side_by_side.png, "
          "water_sanity_sites.png")


if __name__ == "__main__":
    main()
