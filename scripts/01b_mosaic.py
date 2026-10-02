"""Step 1b: check the scenes line up, then mosaic rows 036 + 037 per date.

Reads data/metadata.json (from 01_download.py) and the raw band files, and
writes one mosaicked GeoTIFF per band and date to data/clipped/, covering the
Kashmir valley clip box. Pixel values are copied unchanged: no scaling, no
resampling. A report with the checks, the overlap comparison and the remaining
empty fraction is written to data/mosaic_report.json.

Checks, each of which stops the script instead of resampling or guessing:
  - each row 037 scene is from the same satellite pass as its row 036
    partner (STAC acquisition times within 60 s)
  - all scenes share one CRS, a 30 m north-up pixel size, and the same pixel
    grid (origins differ by a whole number of pixels)

Mosaic rule: a pixel is taken from row 036 if row 036 has valid data there,
otherwise from row 037. "Valid" is decided once per pixel for all five bands
together: QA_PIXEL fill bit (bit 0) is 0 and none of the four reflectance
bands equals its fill value 0. Using one decision for all bands means every
band of a given pixel comes from the same scene.

Usage:  .venv/bin/python scripts/01b_mosaic.py
"""

import json
import math
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import rasterio
from rasterio.features import rasterize
from rasterio.transform import Affine
from rasterio.warp import transform_bounds, transform_geom

ROOT = Path(__file__).resolve().parent.parent
CLIP_DIR = ROOT / "data" / "clipped"
DATES = {"pre": "2014-08-25", "flood": "2014-09-10"}
PRIMARY_ROW, SECONDARY_ROW = "036", "037"
REFLECTANCE = ["green", "red", "nir", "swir1"]
BANDS = REFLECTANCE + ["qa_pixel"]
PIXEL = 30.0
MAX_PASS_SECONDS = 60

# Kashmir valley, WGS84 lon/lat (west, south, east, north)
BBOX = (74.2, 33.5, 75.4, 34.6)


def stop(msg):
    sys.exit(f"STOP: {msg}")


def parse_time(s):
    return datetime.fromisoformat(s.replace("Z", "+00:00"))


def check_grids(scenes):
    """Return (crs, phase_x, phase_y) shared by every band of every scene."""
    crs = phase = None
    for name, sc in scenes.items():
        for band, b in sc["bands"].items():
            g = b["grid"]
            a, _b, c, _d, e, f = g["transform"]
            if a != PIXEL or e != -PIXEL or _b != 0 or _d != 0:
                stop(f"{name} {band}: pixel grid is not 30 m north-up: {g['transform']}")
            p = (c % PIXEL, f % PIXEL)
            if crs is None:
                crs, phase = g["crs"], p
            if g["crs"] != crs:
                stop(f"{name} {band}: CRS {g['crs']} differs from {crs}")
            if not all(math.isclose(u, v, abs_tol=1e-6) for u, v in zip(p, phase)):
                stop(f"{name} {band}: grid origin offset {p} differs from {phase}; "
                     "grids are not aligned, would need resampling")
    return crs, phase


def target_grid(crs, phase):
    """Smallest grid-aligned rectangle (in the scene CRS) enclosing BBOX."""
    minx, miny, maxx, maxy = transform_bounds("EPSG:4326", crs, *BBOX, densify_pts=101)
    px, py = phase
    x0 = px + math.floor((minx - px) / PIXEL) * PIXEL
    x1 = px + math.ceil((maxx - px) / PIXEL) * PIXEL
    y0 = py + math.floor((miny - py) / PIXEL) * PIXEL
    y1 = py + math.ceil((maxy - py) / PIXEL) * PIXEL
    width, height = round((x1 - x0) / PIXEL), round((y1 - y0) / PIXEL)
    return Affine(PIXEL, 0, x0, 0, -PIXEL, y1), width, height


def read_on_grid(path, transform, width, height, fill):
    """Copy the part of path that falls on the target grid; rest = fill.

    Offsets are whole pixels (checked), so this is a straight array copy.
    """
    out = np.full((height, width), fill, dtype=np.uint16)
    with rasterio.open(path) as src:
        col = (transform.c - src.transform.c) / PIXEL
        row = (src.transform.f - transform.f) / PIXEL
        if col != int(col) or row != int(row):
            stop(f"{path.name}: non-integer pixel offset {col}, {row}")
        col, row = int(col), int(row)
        # Intersection of target window with the source raster, in source pixels.
        c0, r0 = max(col, 0), max(row, 0)
        c1, r1 = min(col + width, src.width), min(row + height, src.height)
        if c1 > c0 and r1 > r0:
            data = src.read(1, window=((r0, r1), (c0, c1)))
            out[r0 - row:r1 - row, c0 - col:c1 - col] = data
    return out


def overlap_stats(a, b, scale):
    a = a.astype(np.float64)
    b = b.astype(np.float64)
    d = a - b
    return {
        "pixels": int(a.size),
        "identical_fraction": round(float(np.mean(d == 0)), 6),
        "mean_diff_dn_036_minus_037": round(float(d.mean()), 3),
        "median_diff_dn": float(np.median(d)),
        "mean_abs_diff_dn": round(float(np.abs(d).mean()), 3),
        "p1_p99_diff_dn": [float(x) for x in np.percentile(d, [1, 99])],
        "max_abs_diff_dn": float(np.abs(d).max()),
        "mean_abs_diff_reflectance": round(float(np.abs(d).mean()) * scale, 6),
        "mean_value_dn_036": round(float(a.mean()), 1),
        "pearson_r": round(float(np.corrcoef(a, b)[0, 1]), 6),
    }


def write(path, data, transform, crs, fill, tags):
    profile = {
        "driver": "GTiff", "width": data.shape[1], "height": data.shape[0],
        "count": 1, "dtype": data.dtype, "crs": crs, "transform": transform,
        "nodata": fill, "compress": "deflate", "tiled": True,
        "blockxsize": 256, "blockysize": 256,
    }
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(data, 1)
        dst.update_tags(**tags)


def main():
    meta = json.loads((ROOT / "data" / "metadata.json").read_text())
    scenes = meta["scenes"]
    CLIP_DIR.mkdir(parents=True, exist_ok=True)

    report = {"bbox_wgs84": BBOX, "mosaic_rule": (
        "row 036 where row 036 is valid, otherwise row 037; valid = QA_PIXEL "
        "bit 0 (fill) is 0 and all four reflectance bands != 0; one decision "
        "per pixel shared by all bands"), "pass_check": {}, "dates": {}}

    # 1. Same pass.
    for label, date in DATES.items():
        p = scenes[f"{label}_r{PRIMARY_ROW}"]
        s = scenes[f"{label}_r{SECONDARY_ROW}"]
        dt = (parse_time(s["datetime_utc_stac"]) - parse_time(p["datetime_utc_stac"])).total_seconds()
        report["pass_check"][date] = {
            "row036": p["datetime_utc_stac"], "row037": s["datetime_utc_stac"],
            "seconds_apart": round(dt, 3)}
        print(f"{date}: row 037 is {dt:.1f} s after row 036")
        if abs(dt) > MAX_PASS_SECONDS:
            stop(f"{date}: rows 036 and 037 are {dt:.1f} s apart, not the same pass")

    # 2. Same CRS and grid.
    crs, phase = check_grids(scenes)
    transform, width, height = target_grid(crs, phase)
    report["grid"] = {"crs": crs, "pixel_size_m": PIXEL,
                      "grid_origin_offset_mod_30m": list(phase),
                      "transform": list(transform)[:6],
                      "width": width, "height": height}
    print(f"grid: {crs}, {width} x {height} px, origin "
          f"({transform.c}, {transform.f}), all scenes aligned")

    # Pixels whose centre lies inside the lon/lat clip box (the box is not a
    # rectangle in UTM, so the grid rectangle extends slightly beyond it).
    w, s_, e, n = BBOX
    steps = 200
    ring = ([(w + (e - w) * i / steps, s_) for i in range(steps)]
            + [(e, s_ + (n - s_) * i / steps) for i in range(steps)]
            + [(e - (e - w) * i / steps, n) for i in range(steps)]
            + [(w, n - (n - s_) * i / steps) for i in range(steps)] + [(w, s_)])
    box_utm = transform_geom("EPSG:4326", crs, {"type": "Polygon", "coordinates": [ring]})
    in_box = rasterize([(box_utm, 1)], out_shape=(height, width),
                       transform=transform, fill=0, dtype="uint8").astype(bool)

    for label, date in DATES.items():
        rows = {r: scenes[f"{label}_r{r}"] for r in (PRIMARY_ROW, SECONDARY_ROW)}
        arrays, valid, band_disagree = {}, {}, {}
        for r, sc in rows.items():
            arrays[r] = {}
            for band in BANDS:
                b = sc["bands"][band]
                arrays[r][band] = read_on_grid(ROOT / b["raw_file"], transform,
                                               width, height, b["fill_value"])
            qa_ok = (arrays[r]["qa_pixel"] & 1) == 0
            refl_ok = [arrays[r][bd] != sc["bands"][bd]["fill_value"] for bd in REFLECTANCE]
            valid[r] = qa_ok & np.logical_and.reduce(refl_ok)
            # How often do the bands disagree about whether a pixel has data?
            any_ok = qa_ok | np.logical_or.reduce(refl_ok)
            band_disagree[r] = int(np.sum(any_ok & ~valid[r]))

        source = np.zeros((height, width), dtype=np.uint8)
        source[valid[SECONDARY_ROW]] = 2
        source[valid[PRIMARY_ROW]] = 1
        tags = {"MOSAIC_RULE": report["mosaic_rule"],
                "SOURCE_ROW036": rows[PRIMARY_ROW]["scene_id"],
                "SOURCE_ROW037": rows[SECONDARY_ROW]["scene_id"]}
        for band in BANDS:
            fill = rows[PRIMARY_ROW]["bands"][band]["fill_value"]
            if rows[SECONDARY_ROW]["bands"][band]["fill_value"] != fill:
                stop(f"{date} {band}: rows use different fill values")
            out = np.full((height, width), fill, dtype=np.uint16)
            out[source == 1] = arrays[PRIMARY_ROW][band][source == 1]
            out[source == 2] = arrays[SECONDARY_ROW][band][source == 2]
            path = CLIP_DIR / f"{label}_{date}_{band}.tif"
            write(path, out, transform, crs, fill, {**tags, "BAND": band})
            print(f"  wrote {path.relative_to(ROOT)}")
        write(CLIP_DIR / f"{label}_{date}_source.tif", source, transform, crs, 0,
              {**tags, "VALUES": "0 = no valid data, 1 = row 036, 2 = row 037"})

        # Overlap comparison: pixels valid in both scenes.
        both = valid[PRIMARY_ROW] & valid[SECONDARY_ROW]
        qa36, qa37 = arrays[PRIMARY_ROW]["qa_pixel"], arrays[SECONDARY_ROW]["qa_pixel"]
        clear_both = both & (((qa36 >> 6) & 1) == 1) & (((qa37 >> 6) & 1) == 1)
        cloud36, cloud37 = (qa36 >> 3) & 1, (qa37 >> 3) & 1
        overlap = {
            "pixels_valid_in_both": int(both.sum()),
            "pixels_clear_in_both": int(clear_both.sum()),
            "qa_pixel_identical_fraction": round(float(np.mean(qa36[both] == qa37[both])), 6),
            "cloud_bit_agreement_fraction": round(float(np.mean(cloud36[both] == cloud37[both])), 6),
            "bands_all_valid_pixels": {}, "bands_clear_pixels": {},
        }
        for band in REFLECTANCE:
            scale = rows[PRIMARY_ROW]["bands"][band]["scale_stac"]
            a, b = arrays[PRIMARY_ROW][band], arrays[SECONDARY_ROW][band]
            overlap["bands_all_valid_pixels"][band] = overlap_stats(a[both], b[both], scale)
            overlap["bands_clear_pixels"][band] = overlap_stats(a[clear_both], b[clear_both], scale)

        empty = source == 0
        report["dates"][date] = {
            "row036_scene": rows[PRIMARY_ROW]["scene_id"],
            "row037_scene": rows[SECONDARY_ROW]["scene_id"],
            "pixels_from_row036": int((source == 1).sum()),
            "pixels_from_row037": int((source == 2).sum()),
            "pixels_partial_band_data_treated_as_invalid": band_disagree,
            "empty_pixels_grid_rectangle": int(empty.sum()),
            "empty_fraction_grid_rectangle": round(float(empty.mean()), 8),
            "pixels_in_lonlat_box": int(in_box.sum()),
            "empty_pixels_lonlat_box": int(empty[in_box].sum()),
            "empty_fraction_lonlat_box": round(float(empty[in_box].mean()), 8),
            "empty_fraction_row036_only_lonlat_box":
                round(float((~valid[PRIMARY_ROW])[in_box].mean()), 6),
            "overlap": overlap,
        }
        print(f"  {date}: empty after mosaic {empty[in_box].mean():.4%} of lon/lat box "
              f"(row 036 alone: {(~valid[PRIMARY_ROW])[in_box].mean():.4%})")
        print(f"  overlap {both.sum()} px; red mean |diff| "
              f"{overlap['bands_all_valid_pixels']['red']['mean_abs_diff_dn']} DN")

    out = ROOT / "data" / "mosaic_report.json"
    out.write_text(json.dumps(report, indent=2))
    print(f"wrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
