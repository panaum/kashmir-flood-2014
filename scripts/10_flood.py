"""Step 3c: flood extent under the rules fixed in scripts/step3_rules.md.

Implements rules R1-R9, the unflagged-cloud estimate U1-U5, the OSM split
O1-O3 and the QA water-bit check I1, exactly as written in step3_rules.md
(committed before this script was run). Inputs: the Landsat mosaics, the
resampled DEM, QA_RADSAT mosaics (08_radsat.py), the OSM mask (09_osm.py)
and the step 2 record data/water_report.json.

Outputs: data/flood/flood_classes.tif, data/flood/flood_mask.tif and the
record data/flood_report.json (all numbers in METHODS.md come from it).

Usage:  .venv/bin/python scripts/10_flood.py
"""

import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
import rasterio
from rasterio.warp import transform as warp_transform
from rasterio.warp import transform_bounds
from rasterio.windows import from_bounds
from scipy import ndimage

ROOT = Path(__file__).resolve().parent.parent
CLIP_DIR = ROOT / "data" / "clipped"
FLOOD_DIR = ROOT / "data" / "flood"
DATES = {"pre": "2014-08-25", "flood": "2014-09-10"}
PRE, POST = DATES["pre"], DATES["flood"]
PIXEL_KM2 = 0.0009                                   # R8
SRINAGAR_POINT = (34.08, 74.80)                       # R2
SRINAGAR_BBOX = (74.65, 33.95, 74.95, 34.20)          # R9
FILL, DILATED, CIRRUS, CLOUD, SHADOW, SNOW, CLEAR, WATER_BIT = 0, 1, 2, 3, 4, 5, 6, 7
SAT_GREEN_BIT, SAT_SWIR1_BIT, OCCLUSION_BIT = 2, 5, 11  # R3
# U2 / U3 (fixed in step3_rules.md)
BRIGHT_MIN, WHITENESS_MAX = 0.20, 0.7
HEIGHTS_M = np.arange(500, 5000 + 1, 250)
DARK_NIR_MAX, DARK_SWIR1_MAX = 0.15, 0.10

_spec = importlib.util.spec_from_file_location("wa", Path(__file__).parent / "05_water.py")
wa = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(wa)


def km2(n):
    return round(float(np.count_nonzero(n) if isinstance(n, np.ndarray) else n) * PIXEL_KM2, 4)


def bit(a, b):
    return ((a >> b) & 1).astype(bool)


def read(path):
    with rasterio.open(path) as src:
        return src.read(1)


def load_date(meta, label, date, area):
    bands = {b: read(CLIP_DIR / f"{label}_{date}_{b}.tif")
             for b in ("green", "red", "nir", "swir1", "qa_pixel", "qa_radsat")}
    refl = {}
    for b in ("green", "red", "nir", "swir1"):
        scale, offset = wa.scale_offset(meta, label, b)
        refl[b] = bands[b].astype(np.float64) * scale + offset
    mndwi, g, s, _ = wa.mndwi_from_dn(bands["green"], bands["swir1"],
                                      wa.scale_offset(meta, label, "green"),
                                      wa.scale_offset(meta, label, "swir1"))
    qa, rs = bands["qa_pixel"], bands["qa_radsat"]
    nodata = bit(qa, FILL) | (bands["green"] == 0) | (bands["swir1"] == 0)
    flags = {"no_data": nodata, "dilated_cloud": bit(qa, DILATED), "cloud": bit(qa, CLOUD),
             "cloud_shadow": bit(qa, SHADOW), "snow": bit(qa, SNOW),
             "saturated_green_or_swir1": bit(rs, SAT_GREEN_BIT) | bit(rs, SAT_SWIR1_BIT),
             "green_le0_swir1_le0": (g <= 0) & (s <= 0),
             "green_le0_swir1_gt0": (g <= 0) & (s > 0)}
    invalid = np.zeros(area.shape, bool)
    for f in flags.values():
        invalid |= f
    valid = area & ~invalid                                         # R5
    swir_le0_water = valid & (s <= 0) & (g > 0)                     # rule 4
    return {"mndwi": mndwi, "refl": refl, "qa": qa, "radsat": rs, "valid": valid,
            "flags": flags, "swir_le0_water": swir_le0_water}


def water(d, t):
    with np.errstate(invalid="ignore"):
        return d["valid"] & ((d["mndwi"] > t) | d["swir_le0_water"])   # R5


def outcome(pre, post, t_pre, t_post, both):
    wp, wq = water(pre, t_pre) & both, water(post, t_post) & both
    return {"t_pre": t_pre, "t_post": t_post,
            "flooded_km2": km2(wq & ~wp), "persistent_km2": km2(wp & wq),
            "lost_km2": km2(wp & ~wq), "water_pre_km2": km2(wp),
            "water_post_km2": km2(wq)}, wp, wq


def kappa(a, b):
    n = a.size
    po = np.mean(a == b)
    pe = (a.mean() * b.mean()) + ((1 - a.mean()) * (1 - b.mean()))
    return round(float((po - pe) / (1 - pe)), 4) if pe < 1 else None


def main():
    FLOOD_DIR.mkdir(parents=True, exist_ok=True)
    meta = json.loads((ROOT / "data" / "metadata.json").read_text())
    wr = json.loads((ROOT / "data" / "water_report.json").read_text())
    rules = (Path(__file__).parent / "step3_rules.md").read_text()

    with rasterio.open(CLIP_DIR / "dem_cop30_bilinear.tif") as src:
        dem = src.read(1)
        profile, transform, crs = src.profile, src.transform, src.crs
    floor = (dem != -9999.0) & (dem <= wr["valley_floor_max_m"])
    labels, _ = ndimage.label(floor, structure=np.ones((3, 3)))
    xs, ys = warp_transform("EPSG:4326", crs, [SRINAGAR_POINT[1]], [SRINAGAR_POINT[0]])
    c, r = ~transform * (xs[0], ys[0])
    area = labels == labels[int(r), int(c)]                        # R2
    if km2(area) != wr["floor_connected_to_srinagar_km2"]:
        sys.exit("STOP: analysis area differs from the step 2 connected patch")

    d = {PRE: load_date(meta, "pre", PRE, area), POST: load_date(meta, "flood", POST, area)}
    both = d[PRE]["valid"] & d[POST]["valid"]

    t_sep = wr["dates"][POST]["thresholds"]["otsu"]                # R1
    t_aug = wr["dates"][PRE]["thresholds"]["otsu"]
    if round(t_sep, 3) != 0.127:
        sys.exit(f"STOP: recorded 2014-09-10 Otsu value is {t_sep}, rules say 0.127")
    primary, wp, wq = outcome(d[PRE], d[POST], t_sep, t_sep, both)
    flooded = wq & ~wp
    variants = {"a_own_otsu": outcome(d[PRE], d[POST], t_aug, t_sep, both)[0],
                "b_zero": outcome(d[PRE], d[POST], 0.0, 0.0, both)[0],
                "c_minus_0.05": outcome(d[PRE], d[POST], round(t_sep - 0.05, 6),
                                        round(t_sep - 0.05, 6), both)[0],
                "c_plus_0.05": outcome(d[PRE], d[POST], round(t_sep + 0.05, 6),
                                       round(t_sep + 0.05, 6), both)[0]}

    # Not assessed (R6), with reasons per date.
    not_assessed = area & ~both
    reasons = {date: {k: km2(f & area) for k, f in d[date]["flags"].items()} for date in d}
    na = {"km2": km2(not_assessed),
          "fraction_of_area": round(float(not_assessed.sum() / area.sum()), 6),
          "invalid_pre_km2": km2(area & ~d[PRE]["valid"]),
          "invalid_post_km2": km2(area & ~d[POST]["valid"]),
          "reasons_km2_by_date": reasons,
          "terrain_occluded_km2": {date: km2(bit(d[date]["radsat"], OCCLUSION_BIT) & area)
                                   for date in d}}

    # Rule 4 effect: water only because SWIR1 <= 0 and green > 0.
    rule4 = {date: km2(d[date]["swir_le0_water"]) for date in d}
    rule4_flooded = km2(flooded & d[POST]["swir_le0_water"])

    # OSM split (O3).
    osm = read(ROOT / "data" / "osm" / "water_wetland_mask.tif").astype(bool)
    osm_split = {"flooded_inside_km2": km2(flooded & osm),
                 "flooded_outside_km2": km2(flooded & ~osm),
                 "persistent_inside_km2": km2(wp & wq & osm),
                 "persistent_outside_km2": km2(wp & wq & ~osm),
                 "osm_polygons_in_analysis_area_km2": km2(osm & area),
                 "osm_polygons_in_area_not_assessed_km2": km2(osm & not_assessed)}

    # Srinagar box (R9).
    sb = transform_bounds("EPSG:4326", crs, *SRINAGAR_BBOX, densify_pts=21)
    win = from_bounds(*sb, transform=transform).round_offsets(op="floor") \
        .round_lengths(op="ceil").toslices()
    srinagar = {"pixels_in_box": int(area[win].size),
                "analysis_area_km2": km2(area[win]),
                "flooded_km2": km2(flooded[win]),
                "persistent_km2": km2((wp & wq)[win]),
                "not_assessed_km2": km2(not_assessed[win])}

    # QA_PIXEL water bit (I1), primary threshold.
    qa_check = {}
    for date in d:
        v = d[date]["valid"]
        ours = water(d[date], t_sep)[v]
        usgs = bit(d[date]["qa"], WATER_BIT)[v]
        qa_check[date] = {
            "pixels": int(v.sum()), "agreement": round(float(np.mean(ours == usgs)), 6),
            "ours_water_usgs_not_km2": km2(ours & ~usgs),
            "usgs_water_ours_not_km2": km2(~ours & usgs),
            "both_water_km2": km2(ours & usgs),
            "ours_water_km2": km2(ours), "usgs_water_km2": km2(usgs),
            "cohen_kappa": kappa(ours, usgs)}
        # Spectral make-up of agreeing and disagreeing pixels.
        R = d[date]["refl"]
        with np.errstate(invalid="ignore", divide="ignore"):
            ndvi = ((R["nir"] - R["red"]) / (R["nir"] + R["red"]))[v]
        groups = {"both": ours & usgs, "ours_only": ours & ~usgs, "usgs_only": ~ours & usgs}
        qa_check[date]["spectra"] = {
            k: {"pixels": int(m.sum()),
                **{f"median_{b}": round(float(np.median(R[b][v][m])), 4)
                   for b in ("green", "red", "nir", "swir1")},
                "median_ndvi": round(float(np.median(ndvi[m])), 4),
                "median_mndwi": round(float(np.nanmedian(d[date]["mndwi"][v][m])), 4),
                "fraction_nir_lt_0.11": round(float(np.mean(R["nir"][v][m] < 0.11)), 4)}
            for k, m in groups.items() if m.any()}
    fl_v = flooded[both]
    usgs_flood = (bit(d[POST]["qa"], WATER_BIT) & ~bit(d[PRE]["qa"], WATER_BIT))[both]
    qa_check["flood_from_usgs_water_bit_km2"] = km2(usgs_flood)
    qa_check["flood_agreement"] = round(float(np.mean(fl_v == usgs_flood)), 6)
    qa_check["flood_cohen_kappa"] = kappa(fl_v, usgs_flood)

    # Unflagged clouds on 2014-09-10 (U1-U5).
    post = d[POST]
    vpost = post["valid"]
    conf = (post["qa"] >> 8) & 3
    u1 = vpost & (conf == 2) & ~bit(post["qa"], CLOUD)

    def bright_white(dd, base):
        R = dd["refl"]
        bright = (R["green"] > BRIGHT_MIN) & (R["red"] > BRIGHT_MIN) & \
                 (R["nir"] > BRIGHT_MIN) & (R["swir1"] > BRIGHT_MIN)
        vis = np.stack([R["green"], R["red"], R["nir"]])
        with np.errstate(invalid="ignore", divide="ignore"):
            white = (vis.max(0) - vis.min(0)) / vis.mean(0) < WHITENESS_MAX
        return base & bright & white

    u2 = bright_white(post, vpost)
    mtl = json.loads((ROOT / meta["scenes"]["flood_r036"]["mtl_file"]).read_text())
    attrs = mtl["LANDSAT_METADATA_FILE"]["IMAGE_ATTRIBUTES"]
    az, el = float(attrs["SUN_AZIMUTH"]), float(attrs["SUN_ELEVATION"])
    zone = np.zeros(area.shape, bool)
    rr, cc = np.nonzero(u2)
    for h in HEIGHTS_M:
        dist = h / math.tan(math.radians(el))
        dcol = -dist * math.sin(math.radians(az)) / 30.0   # away from the sun
        drow = dist * math.cos(math.radians(az)) / 30.0
        r2 = np.round(rr + drow).astype(int)
        c2 = np.round(cc + dcol).astype(int)
        ok = (r2 >= 0) & (r2 < area.shape[0]) & (c2 >= 0) & (c2 < area.shape[1])
        zone[r2[ok], c2[ok]] = True
    zone &= vpost & ~u2
    dark = zone & (post["refl"]["nir"] < DARK_NIR_MAX) & (post["refl"]["swir1"] < DARK_SWIR1_MAX)
    pre_clear = d[PRE]["valid"] & bit(d[PRE]["qa"], CLEAR)
    u4 = bright_white(d[PRE], pre_clear)
    unflagged = {
        "sun_azimuth_deg": az, "sun_elevation_deg": el,
        "valid_post_km2": km2(vpost),
        "u1_medium_conf_cloud_km2": km2(u1), "u1_flooded_km2": km2(u1 & flooded),
        "u2_bright_white_km2": km2(u2), "u2_fraction_of_valid": round(float(u2.sum() / vpost.sum()), 6),
        "u2_flooded_km2": km2(u2 & flooded),
        "u3_shadow_zone_km2": km2(zone), "u3_shadow_zone_flooded_km2": km2(zone & flooded),
        "u3_dark_km2": km2(dark), "u3_dark_flooded_km2": km2(dark & flooded),
        "u4_control_pre_clear_km2": km2(pre_clear),
        "u4_control_bright_white_km2": km2(u4),
        "u4_control_fraction": round(float(u4.sum() / pre_clear.sum()), 6),
    }

    # Output rasters.
    classes = np.full(area.shape, 255, dtype=np.uint8)
    classes[both] = 0
    classes[flooded] = 1
    classes[wp & wq] = 2
    classes[wp & ~wq] = 3
    base = {**profile, "dtype": "uint8", "nodata": 255, "compress": "deflate",
            "tiled": True, "blockxsize": 256, "blockysize": 256}
    with rasterio.open(FLOOD_DIR / "flood_classes.tif", "w", **base) as dst:
        dst.write(classes, 1)
        dst.update_tags(VALUES="0 = not water on either date, 1 = flooded, 2 = persistent water, "
                               "3 = water lost, 255 = outside analysis area or not assessed",
                        THRESHOLD=f"{t_sep} both dates", RULES="scripts/step3_rules.md")
    with rasterio.open(FLOOD_DIR / "analysis_area.tif", "w", **{**base, "nodata": None}) as dst:
        dst.write(area.astype(np.uint8), 1)
        dst.update_tags(VALUES="1 = valley floor (DEM <= 1800 m) connected to Srinagar (rule R2)")
    fm = np.full(area.shape, 255, dtype=np.uint8)
    fm[both] = 0
    fm[flooded] = 1
    with rasterio.open(FLOOD_DIR / "flood_mask.tif", "w", **base) as dst:
        dst.write(fm, 1)
        dst.update_tags(VALUES="1 = flooded, 0 = assessed and not flooded, "
                               "255 = outside analysis area or not assessed",
                        THRESHOLD=f"{t_sep} both dates", RULES="scripts/step3_rules.md")

    report = {"rules_file": "scripts/step3_rules.md", "rules_sha256":
              __import__("hashlib").sha256(rules.encode()).hexdigest(),
              "analysis_area_km2": km2(area), "valid_both_km2": km2(both),
              "primary": primary, "variants": variants, "not_assessed": na,
              "rule4_swir_le0_water_km2": rule4, "rule4_flooded_km2": rule4_flooded,
              "saturated_in_area_km2": {date: km2(d[date]["flags"]["saturated_green_or_swir1"] & area)
                                        for date in d},
              "osm": osm_split, "srinagar": srinagar, "qa_water_bit": qa_check,
              "unflagged_clouds": unflagged}
    (ROOT / "data" / "flood_report.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
