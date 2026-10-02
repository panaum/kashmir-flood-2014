"""Final step: collect headline numbers and display layers for the website.

Reads only recorded outputs (data/*_report.json, data/metadata.json, the
rasters in data/flood/ and data/clipped/, scripts/published_figures.json) and
writes, for the static site in docs/:

  docs/data/results.json   every number the page shows, with its source
  docs/data/layers.json    how each display layer was made, and its bounds
  docs/img/flood_display.png          flood classes, Web Mercator (display only)
  docs/img/base_2014-08-25.jpg        Landsat false-colour base map (display only)
  docs/img/srinagar_2014-08-25.jpg    before/after slider images
  docs/img/srinagar_2014-09-10.jpg
  docs/data/analysis_area.geojson     simplified outline of the analysis area

Nothing here changes an analysis result. Display layers are reprojected for
the web map only; all counts and areas come from the analysis files.

Usage:  .venv/bin/python scripts/18_site_data.py
"""

import hashlib
import json
from pathlib import Path

import numpy as np
import rasterio
import shapely
from PIL import Image
from rasterio.features import shapes
from rasterio.warp import (Resampling, calculate_default_transform, reproject,
                           transform_bounds, transform_geom)
from rasterio.windows import from_bounds
from shapely.geometry import mapping, shape

ROOT = Path(__file__).resolve().parent.parent
DOCS = ROOT / "docs"
REPO = "https://github.com/panaum/kashmir-flood-2014"
SRINAGAR_BBOX = (74.65, 33.95, 74.95, 34.20)
BASE_RES_M = 80.0          # base-map pixel size in Web Mercator metres (display only)
OUTLINE_TOL_M = 30.0       # outline simplification tolerance (display only)

# Display palette: Okabe-Ito colours (colour-blind safe), index -> RGB.
CLASSES = {1: ("flooded", (230, 159, 0)),          # orange
           2: ("persistent_water", (0, 114, 178)),   # blue
           3: ("water_lost", (240, 228, 66)),        # yellow
           4: ("not_assessed", (204, 121, 167))}     # reddish purple

WORDS = {10000: "ten thousand", 20000: "twenty thousand", 30000: "thirty thousand",
         40000: "forty thousand", 50000: "fifty thousand", 60000: "sixty thousand",
         70000: "seventy thousand", 80000: "eighty thousand", 90000: "ninety thousand"}


def load(name):
    p = ROOT / "data" / name
    return json.loads(p.read_text()), hashlib.sha256(p.read_bytes()).hexdigest()


def about(n):
    """Round to the nearest ten thousand, in words (for plain-language text)."""
    r = int(round(n / 10000.0) * 10000)
    if r not in WORDS:
        raise SystemExit(f"STOP: no wording for {n} (rounds to {r})")
    return r, WORDS[r]


def hhmm(iso):
    """Acquisition time rounded to the nearest minute (UTC), as HH:MM."""
    from datetime import datetime, timedelta
    d = datetime.fromisoformat(iso.replace("Z", "+00:00")) + timedelta(seconds=30)
    return d.strftime("%H:%M")


def stretch(paths, lims):
    rgb = np.stack([rasterio.open(p).read(1).astype(np.float32) for p in paths], axis=-1)
    for c, (lo, hi) in enumerate(lims):
        rgb[..., c] = np.clip((rgb[..., c] - lo) / (hi - lo), 0, 1)
    return (rgb * 255).round().astype(np.uint8)


def to_mercator(arr, src_transform, src_crs, resampling, res=None, fill=0):
    h, w = arr.shape[:2]
    left, top = src_transform * (0, 0)
    right, bottom = src_transform * (w, h)
    kw = {"resolution": res} if res else {}
    dst_t, dw, dh = calculate_default_transform(src_crs, "EPSG:3857", w, h,
                                                left, bottom, right, top, **kw)
    bands = arr[..., None] if arr.ndim == 2 else arr
    out = np.full((dh, dw, bands.shape[2]), fill, dtype=arr.dtype)
    for b in range(bands.shape[2]):
        reproject(bands[..., b], out[..., b], src_transform=src_transform, src_crs=src_crs,
                  dst_transform=dst_t, dst_crs="EPSG:3857", resampling=resampling,
                  src_nodata=None, dst_nodata=fill)
    west, south, east, north = transform_bounds("EPSG:3857", "EPSG:4326",
                                                dst_t.c, dst_t.f + dh * dst_t.e,
                                                dst_t.c + dw * dst_t.a, dst_t.f)
    return (out[..., 0] if arr.ndim == 2 else out), dst_t, [[south, west], [north, east]]


def main():
    (DOCS / "data").mkdir(parents=True, exist_ok=True)
    (DOCS / "img").mkdir(parents=True, exist_ok=True)
    meta, h_meta = load("metadata.json")
    wr, h_wr = load("water_report.json")
    fr, h_fr = load("flood_report.json")
    osm, h_osm = load("osm_report.json")
    ex, h_ex = load("exposure_report.json")
    ms, h_ms = load("msbuildings_report.json")
    dem, h_dem = load("dem_report.json")
    pub = json.loads((ROOT / "scripts" / "published_figures.json").read_text())
    lims = list(json.loads((ROOT / "previews" / "qa_summary.json").read_text())
                ["stretch_2_98_percentile_raw_DN"].values())

    p, v = fr["primary"], fr["variants"]
    variants = {"primary": p, **v}
    lo = min(variants.values(), key=lambda x: x["flooded_km2"])
    hi = max(variants.values(), key=lambda x: x["flooded_km2"])
    msc = ms["counts"]["centroid"]
    about_n, about_words = about(msc["analysis_area"]["flooded"])
    sc = {k: s for k, s in meta["scenes"].items()}
    results = {
        "generated_by": "scripts/18_site_data.py",
        "sources": {"data/metadata.json": h_meta, "data/water_report.json": h_wr,
                    "data/flood_report.json": h_fr, "data/osm_report.json": h_osm,
                    "data/exposure_report.json": h_ex, "data/msbuildings_report.json": h_ms,
                    "data/dem_report.json": h_dem},
        "dates": {"before": sc["pre_r036"]["date"], "flood": sc["flood_r036"]["date"],
                  "flood_time_utc": hhmm(sc["flood_r036"]["datetime_utc_stac"]),
                  "before_time_utc": hhmm(sc["pre_r036"]["datetime_utc_stac"])},
        "scenes": [s["scene_id"] for s in meta["scenes"].values()],
        "method": {"pixel_m": 30, "valley_floor_max_m": wr["valley_floor_max_m"],
                   "threshold": round(p["t_post"], 3),
                   "analysis_area_km2": fr["analysis_area_km2"]},
        "flood": {"flooded_km2": p["flooded_km2"],
                  "range_km2": [lo["flooded_km2"], hi["flooded_km2"]],
                  "range_thresholds": [[lo["t_pre"], lo["t_post"]], [hi["t_pre"], hi["t_post"]]],
                  "variants": {k: {"t_pre": round(x["t_pre"], 3), "t_post": round(x["t_post"], 3),
                                   "flooded_km2": x["flooded_km2"]} for k, x in variants.items()},
                  "persistent_km2": p["persistent_km2"], "lost_km2": p["lost_km2"],
                  "not_assessed_km2": fr["not_assessed"]["km2"],
                  "not_assessed_pct": round(100 * fr["not_assessed"]["fraction_of_area"], 3),
                  "inside_lakes_wetlands_km2": fr["osm"]["flooded_inside_km2"],
                  "outside_lakes_wetlands_km2": fr["osm"]["flooded_outside_km2"],
                  "srinagar_box_km2": fr["srinagar"]["flooded_km2"],
                  "unflagged_cloud_u2_km2": fr["unflagged_clouds"]["u2_bright_white_km2"],
                  "usgs_water_bit_kappa_flood_date": fr["qa_water_bit"][sc["flood_r036"]["date"]]["cohen_kappa"]},
        "roads": {"motor_flooded_km": ex["roads"]["motor"]["analysis_area"]["flooded"],
                  "motor_total_km": ex["roads"]["motor"]["analysis_area"]["total_in_analysis_area"],
                  "motor_srinagar_flooded_km": ex["roads"]["motor"]["srinagar_box"]["flooded"],
                  "motor_srinagar_total_km": ex["roads"]["motor"]["srinagar_box"]["total_in_analysis_area"],
                  "osm_timestamp": ex["osm_source"]["extract"]["osm_timestamp_header"]},
        "buildings": {"footprints_flooded": msc["analysis_area"]["flooded"],
                      "footprints_flooded_overlap": ms["counts"]["overlap"]["analysis_area"]["flooded"],
                      "footprints_total": msc["analysis_area"]["total_in_analysis_area"],
                      "footprints_srinagar_flooded": msc["srinagar_box"]["flooded"],
                      "footprints_srinagar_total": msc["srinagar_box"]["total_in_analysis_area"],
                      "about": about_n, "about_words": about_words,
                      "imagery_start": min(i["start_datetime"] for i in ms["items"])[:4],
                      "imagery_end": max(i["end_datetime"] for i in ms["items"])[:4],
                      "ms_version": ms["source"]["dataset_version_processing_date"],
                      "osm_mapped_flooded": ex["buildings"]["centroid"]["analysis_area"]["flooded"],
                      "osm_cells_roads_no_buildings": ex["completeness"]["c1_cells_roads_no_buildings"],
                      "ms_cells_roads_no_buildings": ms["completeness"]["c1_cells_roads_no_buildings"],
                      "cells_with_roads": ms["completeness"]["cells_with_ge_2km_motor_road"]},
        "published": pub,
        "links": {"repository": REPO, "methods": f"{REPO}/blob/main/METHODS.md",
                  "data": f"{REPO}/blob/main/DATA.md", "readme": f"{REPO}#readme"},
    }

    layers = {"note": "Display layers for the web map only. Counts and areas come from the "
                      "analysis files, not from these images.", "layers": {}}

    # Flood classes -> paletted PNG in Web Mercator, nearest neighbour.
    with rasterio.open(ROOT / "data" / "flood" / "flood_classes.tif") as src:
        cls = src.read(1)
        T, crs = src.transform, src.crs
    area = rasterio.open(ROOT / "data" / "flood" / "analysis_area.tif").read(1).astype(bool)
    # Background gets its own value (5) during reprojection: with no source
    # nodata given, the warp treats source zeros as missing and fills them.
    disp = np.full(cls.shape, 5, np.uint8)
    disp[cls == 1] = 1
    disp[cls == 2] = 2
    disp[cls == 3] = 3
    disp[area & (cls == 255)] = 4
    merc, mt, bounds = to_mercator(disp, T, crs, Resampling.nearest)
    merc[merc == 5] = 0
    disp[disp == 5] = 0
    img = Image.fromarray(merc, mode="P")
    pal = [0, 0, 0]
    for k in range(1, 5):
        pal += list(CLASSES[k][1])
    img.putpalette(pal + [0, 0, 0] * (256 - 5))
    img.save(DOCS / "img" / "flood_display.png", optimize=True, transparency=0)
    counts_src = {CLASSES[k][0]: int((disp == k).sum()) for k in CLASSES}
    counts_disp = {CLASSES[k][0]: int((merc == k).sum()) for k in CLASSES}
    layers["layers"]["flood"] = {
        "file": "img/flood_display.png", "bounds_latlon": bounds,
        "source": "data/flood/flood_classes.tif + data/flood/analysis_area.tif (EPSG:32643, 30 m)",
        "method": "classes 1 flooded, 2 persistent water, 3 water lost, 4 not assessed; reprojected "
                  "to EPSG:3857 with nearest-neighbour resampling (rasterio calculate_default_transform), "
                  "written as a paletted PNG with index 0 transparent",
        "size_px": [int(merc.shape[1]), int(merc.shape[0])],
        "pixel_m_mercator": round(mt.a, 2),
        "colours": {CLASSES[k][0]: "#%02X%02X%02X" % CLASSES[k][1] for k in CLASSES},
        "pixels_source": counts_src, "pixels_display": counts_disp}

    # Base map: 2014-08-25 false colour, step 1 stretch, bilinear, 60 m.
    base = stretch([ROOT / "data" / "clipped" / f"pre_2014-08-25_{b}.tif" for b in ("swir1", "nir", "red")], lims)
    bm, bt, bb = to_mercator(base, T, crs, Resampling.bilinear, res=BASE_RES_M, fill=246)
    Image.fromarray(bm, mode="RGB").save(DOCS / "img" / "base_2014-08-25.jpg", quality=80, optimize=True, progressive=True)
    layers["layers"]["base"] = {
        "file": "img/base_2014-08-25.jpg", "bounds_latlon": bb,
        "source": "data/clipped/pre_2014-08-25_{swir1,nir,red}.tif",
        "method": f"R=SWIR1 G=NIR B=Red, linear stretch as DATA.md (2nd-98th percentile), 8-bit; "
                  f"reprojected to EPSG:3857 at {BASE_RES_M:g} m with bilinear resampling; progressive JPEG quality 80",
        "size_px": [int(bm.shape[1]), int(bm.shape[0])]}

    # Before/after slider: Srinagar box, native 30 m grid, no reprojection.
    sb = transform_bounds("EPSG:4326", crs, *SRINAGAR_BBOX, densify_pts=21)
    win = from_bounds(*sb, transform=T).round_offsets(op="floor").round_lengths(op="ceil").toslices()
    for label, date in (("pre", "2014-08-25"), ("flood", "2014-09-10")):
        rgb = stretch([ROOT / "data" / "clipped" / f"{label}_{date}_{b}.tif" for b in ("swir1", "nir", "red")], lims)
        Image.fromarray(rgb[win]).save(DOCS / "img" / f"srinagar_{date}.jpg", quality=88, optimize=True)
    layers["layers"]["slider"] = {
        "files": ["img/srinagar_2014-08-25.jpg", "img/srinagar_2014-09-10.jpg"],
        "bbox_wgs84": SRINAGAR_BBOX,
        "method": "R=SWIR1 G=NIR B=Red, same stretch as DATA.md, native 30 m UTM grid, no "
                  "reprojection, JPEG quality 88",
        "size_px": [int(rgb[win].shape[1]), int(rgb[win].shape[0])]}

    # Analysis area outline, simplified for display.
    polys = [shape(g) for g, val in shapes(area.astype(np.uint8), mask=area, transform=T) if val == 1]
    outline = shapely.simplify(shapely.union_all(polys), OUTLINE_TOL_M)
    gj = transform_geom(crs, "EPSG:4326", mapping(outline))
    (DOCS / "data" / "analysis_area.geojson").write_text(json.dumps(
        {"type": "Feature", "properties": {"name": "analysis area (display outline)"}, "geometry": gj},
        separators=(",", ":")))
    layers["layers"]["outline"] = {
        "file": "data/analysis_area.geojson",
        "method": f"pixel outlines of data/flood/analysis_area.tif, merged, simplified with a "
                  f"{OUTLINE_TOL_M:g} m tolerance (shapely), converted to WGS84"}
    srin = transform_bounds("EPSG:4326", "EPSG:4326", *SRINAGAR_BBOX)
    layers["views"] = {"srinagar": [[srin[1], srin[0]], [srin[3], srin[2]]],
                       "valley": bounds}

    (DOCS / "data" / "results.json").write_text(json.dumps(results, indent=1))
    (DOCS / "data" / "layers.json").write_text(json.dumps(layers, indent=1))
    for f in sorted((DOCS / "img").iterdir()):
        print(f"{f.relative_to(ROOT)}  {f.stat().st_size / 1024:.0f} KB")
    print(json.dumps({k: results[k] for k in ("flood", "roads", "buildings")}, indent=1))
    print(json.dumps(layers["layers"]["flood"]["pixels_source"]), json.dumps(layers["layers"]["flood"]["pixels_display"]))


if __name__ == "__main__":
    main()
