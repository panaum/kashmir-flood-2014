"""Step 3d: flood extent previews for Srinagar and the valley.

Left: 2014-09-10 false colour (R=SWIR1, G=NIR, B=Red, step 1 stretch).
Right: the same image with the flood classes from 10_flood.py drawn on top.
Colours: red = flooded, blue = persistent water, yellow = water lost,
magenta = not assessed, darkened = outside the analysis area.

Usage:  .venv/bin/python scripts/11_flood_preview.py
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from matplotlib.patches import Patch
from rasterio.warp import transform_bounds
from rasterio.windows import from_bounds

ROOT = Path(__file__).resolve().parent.parent
CLIP_DIR = ROOT / "data" / "clipped"
OUT_DIR = ROOT / "previews"
SRINAGAR_BBOX = (74.65, 33.95, 74.95, 34.20)  # as in 02_preview.py
COLOURS = {1: ((1.0, 0.15, 0.1), "flooded"), 2: ((0.1, 0.35, 1.0), "persistent water"),
           3: ((1.0, 0.9, 0.0), "water lost")}
NA_RGB = (1.0, 0.0, 1.0)


def main():
    lims = list(json.loads((OUT_DIR / "qa_summary.json").read_text())
                ["stretch_2_98_percentile_raw_DN"].values())
    rep = json.loads((ROOT / "data" / "flood_report.json").read_text())
    with rasterio.open(ROOT / "data" / "flood" / "flood_classes.tif") as src:
        cls = src.read(1)
        transform, crs = src.transform, src.crs
    with rasterio.open(ROOT / "data" / "flood" / "analysis_area.tif") as src:
        area = src.read(1).astype(bool)
    base = np.stack([rasterio.open(CLIP_DIR / f"flood_2014-09-10_{b}.tif").read(1).astype(np.float32)
                     for b in ("swir1", "nir", "red")], axis=-1)
    for c, (lo, hi) in enumerate(lims):
        base[..., c] = np.clip((base[..., c] - lo) / (hi - lo), 0, 1)
    over = base.copy()
    over[~area] *= 0.4
    for v, (rgb, _) in COLOURS.items():
        over[cls == v] = rgb
    na = area & (cls == 255)
    over[na] = NA_RGB

    sb = transform_bounds("EPSG:4326", crs, *SRINAGAR_BBOX, densify_pts=21)
    sr = from_bounds(*sb, transform=transform).round_offsets(op="floor") \
        .round_lengths(op="ceil").toslices()
    handles = [Patch(color=rgb, label=lab) for rgb, lab in COLOURS.values()] + \
              [Patch(color=NA_RGB, label="not assessed"),
               Patch(color="0.25", label="outside analysis area")]
    p = rep["primary"]
    for name, cut, size, title in (
            ("srinagar", sr, (14, 7.8),
             f"Srinagar box: flooded {rep['srinagar']['flooded_km2']:.1f} km2"),
            ("valley", (slice(None), slice(None)), (14, 8.8),
             f"Analysis area: flooded {p['flooded_km2']:.1f} km2, "
             f"persistent water {p['persistent_km2']:.1f} km2")):
        fig, axes = plt.subplots(1, 2, figsize=size)
        interp = "nearest" if name == "srinagar" else "antialiased"
        axes[0].imshow(base[cut], interpolation=interp)
        axes[0].set_title("2014-09-10 false colour (R=SWIR1, G=NIR, B=Red)")
        axes[1].imshow(over[cut], interpolation=interp)
        axes[1].set_title(f"Flood extent, MNDWI > {p['t_post']:.3f} on both dates")
        for ax in axes:
            ax.axis("off")
        fig.suptitle(title)
        fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False)
        fig.savefig(OUT_DIR / f"flood_{name}.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
    print("wrote previews/flood_srinagar.png, previews/flood_valley.png")


if __name__ == "__main__":
    main()
