"""Step 2b: MNDWI, masks, valley floor, Otsu threshold and water masks.

Reads the Landsat mosaics (01b_mosaic.py) and the resampled DEM (04_dem.py)
and, for each date:
  - converts green and SWIR1 to surface reflectance with the recorded scale
    and offset, and computes MNDWI = (green - swir1) / (green + swir1);
  - builds the exclusion mask from QA_PIXEL (fill, cloud, cloud shadow,
    dilated cloud, snow) plus pixels where green or SWIR1 reflectance <= 0,
    for which MNDWI is not meaningful. Cirrus is reported, not masked;
  - restricts the analysis to the valley floor (DEM <= 1800 m);
  - picks the water threshold with Otsu's method on valid valley-floor
    pixels and reports the water area at several thresholds;
  - writes the water mask (0 = not water, 1 = water, 255 = masked/no data);
  - runs the lake / airport sanity checks and the row 036 / 037 seam check.

Outputs: data/water/*.tif, previews/mndwi_histogram_*.png, and the record
data/water_report.json (all numbers in METHODS.md come from it).

Usage:  .venv/bin/python scripts/05_water.py
"""

import importlib.util
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from rasterio.warp import transform as warp_transform
from scipy import ndimage

ROOT = Path(__file__).resolve().parent.parent
CLIP_DIR = ROOT / "data" / "clipped"
WATER_DIR = ROOT / "data" / "water"
PREVIEW_DIR = ROOT / "previews"
DATES = {"pre": "2014-08-25", "flood": "2014-09-10"}

VALLEY_FLOOR_MAX_M = 1800.0          # analysis choice
OTSU_BINS = 2000                     # histogram bins over [-1, 1], width 0.001
SENSITIVITY = 0.05                   # Otsu +/- this
PIXEL_AREA_KM2 = 30 * 30 / 1e6       # UTM map area of one pixel

# QA_PIXEL bits (STAC classification:bitfields, see DATA.md)
FILL, DILATED, CIRRUS, CLOUD, SHADOW, SNOW = 0, 1, 2, 3, 4, 5
EXCLUDED_FLAGS = {"cloud": CLOUD, "cloud_shadow": SHADOW,
                  "dilated_cloud": DILATED, "snow": SNOW}

# Sanity-check sites, WGS84 (lat, lon); box = +/- SITE_HALF_DEG in both axes.
SITES = {"dal_lake": (34.12, 74.86), "wular_lake": (34.36, 74.56),
         "srinagar_airport": (33.99, 74.77)}
SITE_HALF_DEG = 0.01
SRINAGAR_POINT = (34.08, 74.80)

# Reuse the exact grid-copy code of the mosaic step for the seam check.
_spec = importlib.util.spec_from_file_location("mo", Path(__file__).parent / "01b_mosaic.py")
mo = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mo)


def bit(qa, b):
    return ((qa >> b) & 1).astype(bool)


def km2(n):
    return round(float(n) * PIXEL_AREA_KM2, 4)


def read(path):
    with rasterio.open(path) as src:
        return src.read(1), src.profile


def scale_offset(meta, label, band):
    vals = {(sc["bands"][band]["scale_stac"], sc["bands"][band]["offset_stac"])
            for sc in meta["scenes"].values() if sc["label"] == label}
    if len(vals) != 1:
        raise SystemExit(f"STOP: {label} {band}: scenes disagree on scale/offset {vals}")
    return vals.pop()


def mndwi_from_dn(g_dn, s_dn, g_so, s_so):
    """Reflectance, MNDWI, and a mask of pixels where MNDWI is not meaningful."""
    g = g_dn.astype(np.float64) * g_so[0] + g_so[1]
    s = s_dn.astype(np.float64) * s_so[0] + s_so[1]
    bad = (g <= 0) | (s <= 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        m = (g - s) / (g + s)
    m[bad] = np.nan
    return m, g, s, bad


def otsu(values):
    """Otsu threshold on a fixed histogram over [-1, 1] (OTSU_BINS bins).

    Returns the bin edge that maximises the between-class variance; values
    above it are the 'water' class.
    """
    counts, edges = np.histogram(values, bins=OTSU_BINS, range=(-1.0, 1.0))
    centres = (edges[:-1] + edges[1:]) / 2
    w0 = np.cumsum(counts)
    w1 = w0[-1] - w0
    s0 = np.cumsum(counts * centres)
    m0 = np.divide(s0, w0, out=np.zeros_like(s0), where=w0 > 0)
    m1 = np.divide(s0[-1] - s0, w1, out=np.zeros_like(s0), where=w1 > 0)
    between = w0 * w1 * (m0 - m1) ** 2
    k = int(np.argmax(between[:-1]))
    return float(edges[k + 1]), counts, edges


def stats(d):
    d = d.astype(np.float64)
    if d.size == 0:
        return {"pixels": 0}
    return {"pixels": int(d.size), "mean": round(float(d.mean()), 5),
            "median": round(float(np.median(d)), 5),
            "mean_abs": round(float(np.abs(d).mean()), 5),
            "p1_p99": [round(float(x), 5) for x in np.percentile(d, [1, 99])]}


def write(path, data, profile, dtype, nodata, tags):
    p = {**profile, "dtype": dtype, "nodata": nodata, "count": 1,
         "compress": "deflate", "tiled": True, "blockxsize": 256, "blockysize": 256}
    with rasterio.open(path, "w", **p) as dst:
        dst.write(data.astype(dtype), 1)
        dst.update_tags(**tags)


def lonlat_window(transform, crs, lat, lon, half):
    """Row/col slices of grid pixels whose centres fall inside a lon/lat box."""
    xs, ys = warp_transform("EPSG:4326", crs, [lon - half, lon + half, lon - half, lon + half],
                            [lat - half, lat - half, lat + half, lat + half])
    inv = ~transform
    cols, rows = zip(*[inv * (x, y) for x, y in zip(xs, ys)])
    return (slice(int(np.ceil(min(rows) - 0.5)), int(np.floor(max(rows) - 0.5)) + 1),
            slice(int(np.ceil(min(cols) - 0.5)), int(np.floor(max(cols) - 0.5)) + 1))


def main():
    WATER_DIR.mkdir(parents=True, exist_ok=True)
    meta = json.loads((ROOT / "data" / "metadata.json").read_text())
    dem, profile = read(CLIP_DIR / "dem_cop30_bilinear.tif")
    profile.pop("nodata", None)
    dem_ok = dem != -9999.0
    floor = dem_ok & (dem <= VALLEY_FLOOR_MAX_M)
    transform, crs = profile["transform"], profile["crs"]
    write(WATER_DIR / "valley_floor.tif", floor, profile, "uint8", None,
          {"VALUES": f"1 = DEM <= {VALLEY_FLOOR_MAX_M:g} m, 0 = higher or no DEM"})

    # Diagnostic: which parts of the <= 1800 m area connect to Srinagar?
    labels, n_comp = ndimage.label(floor, structure=np.ones((3, 3)))
    xs, ys = warp_transform("EPSG:4326", crs, [SRINAGAR_POINT[1]], [SRINAGAR_POINT[0]])
    c_, r_ = ~transform * (xs[0], ys[0])
    sr_label = labels[int(r_), int(c_)]
    if sr_label == 0:
        raise SystemExit(f"STOP: Srinagar point {SRINAGAR_POINT} is not on the valley floor")
    connected = labels == sr_label
    sizes = np.bincount(labels.ravel())
    sizes[0] = 0
    sizes[sr_label] = 0
    others = []
    for lab_id in np.argsort(sizes)[::-1][:3]:
        if sizes[lab_id] == 0:
            break
        rr, cc = np.nonzero(labels == lab_id)
        x, y = transform * (cc.mean() + 0.5, rr.mean() + 0.5)
        lon, lat = warp_transform(crs, "EPSG:4326", [x], [y])
        others.append({"km2": km2(sizes[lab_id]),
                       "centroid_lat_lon": [round(lat[0], 3), round(lon[0], 3)]})
    report = {
        "valley_floor_max_m": VALLEY_FLOOR_MAX_M,
        "pixel_area_km2": PIXEL_AREA_KM2,
        "grid_pixels": int(floor.size),
        "valley_floor_pixels": int(floor.sum()),
        "valley_floor_km2": km2(floor.sum()),
        "floor_components_8conn": int(n_comp),
        "floor_connected_to_srinagar_km2": km2(connected.sum()),
        "floor_not_connected_to_srinagar_km2": km2((floor & ~connected).sum()),
        "largest_other_floor_patches": others,
        "otsu_bins": OTSU_BINS, "sensitivity": SENSITIVITY,
        "site_half_deg": SITE_HALF_DEG, "sites": {k: list(v) for k, v in SITES.items()},
        "dates": {},
    }
    print(f"valley floor {report['valley_floor_km2']} km2, of which "
          f"{report['floor_connected_to_srinagar_km2']} km2 connected to Srinagar")

    excluded_any = np.zeros(floor.shape, bool)
    for label, date in DATES.items():
        g_dn, _ = read(CLIP_DIR / f"{label}_{date}_green.tif")
        s_dn, _ = read(CLIP_DIR / f"{label}_{date}_swir1.tif")
        r_dn, _ = read(CLIP_DIR / f"{label}_{date}_red.tif")
        n_dn, _ = read(CLIP_DIR / f"{label}_{date}_nir.tif")
        r_so, n_so = scale_offset(meta, label, "red"), scale_offset(meta, label, "nir")
        red = r_dn.astype(np.float64) * r_so[0] + r_so[1]
        nir = n_dn.astype(np.float64) * n_so[0] + n_so[1]
        with np.errstate(divide="ignore", invalid="ignore"):
            ndvi = (nir - red) / (nir + red)
        qa, _ = read(CLIP_DIR / f"{label}_{date}_qa_pixel.tif")
        g_so, s_so = scale_offset(meta, label, "green"), scale_offset(meta, label, "swir1")
        nodata = bit(qa, FILL) | (g_dn == 0) | (s_dn == 0) | ~dem_ok
        mndwi, g, s, bad = mndwi_from_dn(g_dn, s_dn, g_so, s_so)
        bad &= ~nodata
        mndwi[nodata] = np.nan

        flags = {k: bit(qa, b) & ~nodata for k, b in EXCLUDED_FLAGS.items()}
        flags["reflectance_le_0"] = bad
        excluded = nodata.copy()
        for f in flags.values():
            excluded |= f
        cirrus = bit(qa, CIRRUS) & ~nodata
        valid = ~excluded
        excluded_any |= excluded & floor

        write(WATER_DIR / f"{label}_{date}_mndwi.tif", mndwi, profile, "float32", np.nan,
              {"FORMULA": "(green - swir1) / (green + swir1), surface reflectance",
               "SCALE_OFFSET_GREEN": str(g_so), "SCALE_OFFSET_SWIR1": str(s_so),
               "NAN": "no data, or green/swir1 reflectance <= 0"})

        # How the masks split between valley floor and higher ground.
        split = {}
        for name, f in {**flags, "cirrus_not_masked": cirrus,
                        "cirrus_not_already_excluded": cirrus & ~excluded,
                        "all_excluded_incl_nodata": excluded,
                        "no_data": nodata}.items():
            split[name] = {"floor_km2": km2((f & floor).sum()),
                           "floor_fraction": round(float((f & floor).sum() / floor.sum()), 6),
                           "connected_floor_km2": km2((f & connected).sum()),
                           "connected_floor_fraction": round(float((f & connected).sum()
                                                                   / connected.sum()), 6),
                           "higher_km2": km2((f & ~floor & dem_ok).sum()),
                           "higher_fraction": round(float((f & ~floor & dem_ok).sum()
                                                          / (~floor & dem_ok).sum()), 6)}

        sample = mndwi[valid & floor]
        t_otsu, counts, edges = otsu(sample)
        thresholds = {"otsu": t_otsu, "otsu_minus": round(t_otsu - SENSITIVITY, 6),
                      "otsu_plus": round(t_otsu + SENSITIVITY, 6), "zero": 0.0}
        areas = {k: km2((sample > t).sum()) for k, t in thresholds.items()}
        # If cirrus were also masked, at the Otsu threshold:
        valid_nc = valid & floor & ~cirrus
        area_nc = km2((mndwi[valid_nc] > t_otsu).sum())
        # Pixels with SWIR1 <= 0 but green > 0 on the floor (MNDWI would be > 1).
        swir_le0_floor = (s <= 0) & (g > 0) & floor & ~nodata
        for f in EXCLUDED_FLAGS.values():
            swir_le0_floor &= ~bit(qa, f)

        water = np.full(floor.shape, 255, dtype=np.uint8)
        water[valid] = (mndwi[valid] > t_otsu).astype(np.uint8)
        write(WATER_DIR / f"{label}_{date}_water_mask.tif", water, profile, "uint8", 255,
              {"VALUES": "0 = not water, 1 = water (MNDWI > threshold), 255 = masked or no data",
               "THRESHOLD": f"{t_otsu} (Otsu on valid valley-floor pixels)",
               "MASKED": "fill, cloud, cloud shadow, dilated cloud, snow, "
                         "green or swir1 reflectance <= 0; cirrus NOT masked"})

        # Histogram figure.
        fig, ax = plt.subplots(figsize=(9, 5))
        ax.stairs(counts, edges, fill=True, color="0.55")
        for name, t, ls in (("Otsu", t_otsu, "-"), ("Otsu - 0.05", t_otsu - SENSITIVITY, "--"),
                            ("Otsu + 0.05", t_otsu + SENSITIVITY, "--"), ("0", 0.0, ":")):
            ax.axvline(t, color="k", ls=ls, lw=1.2, label=f"{name} = {t:.3f}")
        ax.set_yscale("log")
        ax.set_xlim(-1, 1)
        ax.set_xlabel("MNDWI")
        ax.set_ylabel(f"pixels per bin (bin width {2 / OTSU_BINS:g}, log scale)")
        ax.set_title(f"{date}: MNDWI of valid valley-floor pixels (DEM <= "
                     f"{VALLEY_FLOOR_MAX_M:g} m), n = {sample.size}")
        ax.legend(loc="upper left", frameon=False)
        fig.savefig(PREVIEW_DIR / f"mndwi_histogram_{date}.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

        # Sanity sites.
        sites = {}
        for name, (lat, lon) in SITES.items():
            rs, cs = lonlat_window(transform, crs, lat, lon, SITE_HALF_DEG)
            w = water[rs, cs]
            n_valid = int((w != 255).sum())
            ok = w != 255
            spectra = {"median_reflectance": {
                           k: round(float(np.median(v[rs, cs][ok])), 4)
                           for k, v in (("green", g), ("red", red), ("nir", nir), ("swir1", s))},
                       "median_ndvi": round(float(np.median(ndvi[rs, cs][ok])), 4),
                       "mndwi_p10_p50_p90": [round(float(x), 4) for x in
                                             np.percentile(mndwi[rs, cs][ok], [10, 50, 90])]} \
                if n_valid else {}
            sites[name] = {"pixels": int(w.size), "valid_pixels": n_valid,
                           "masked_fraction": round(float((w == 255).mean()), 4),
                           "water_fraction_of_valid":
                               round(float((w == 1).sum() / n_valid), 4) if n_valid else None,
                           "on_valley_floor_fraction": round(float(floor[rs, cs].mean()), 4),
                           **spectra}

        report["dates"][date] = {
            "scale_offset_green": g_so, "scale_offset_swir1": s_so,
            "valid_floor_pixels": int(sample.size),
            "valid_floor_km2": km2(sample.size),
            "mndwi_valid_floor": {"min": round(float(sample.min()), 4),
                                  "max": round(float(sample.max()), 4),
                                  "median": round(float(np.median(sample)), 4)},
            "mask_split": split,
            "thresholds": thresholds,
            "water_km2_floor": areas,
            "water_km2_floor_otsu_if_cirrus_masked": area_nc,
            "swir1_le0_green_gt0_clear_floor_pixels": int(swir_le0_floor.sum()),
            "swir1_le0_green_gt0_clear_floor_km2": km2(swir_le0_floor.sum()),
            "water_km2_floor_otsu_connected_to_srinagar":
                km2(((water == 1) & floor & connected).sum()),
            "water_km2_floor_otsu_not_connected":
                km2(((water == 1) & floor & ~connected).sum()),
            "sites": sites,
        }
        print(f"{date}: Otsu {t_otsu:.3f}; water {areas} km2; sites "
              f"{ {k: v['water_fraction_of_valid'] for k, v in sites.items()} }")

        report["dates"][date]["seam"] = seam_check(meta, label, date, transform,
                                                   floor, t_otsu, g_so, s_so)

    report["connected_floor_masked_on_either_date_km2"] = km2((excluded_any & connected).sum())
    report["connected_floor_masked_on_either_date_fraction"] = round(
        float((excluded_any & connected).sum() / connected.sum()), 6)
    report["floor_masked_on_either_date_km2"] = km2(excluded_any.sum())
    report["floor_masked_on_either_date_fraction"] = round(
        float(excluded_any.sum() / floor.sum()), 6)
    (ROOT / "data" / "water_report.json").write_text(json.dumps(report, indent=2))
    print("wrote data/water_report.json")


def seam_check(meta, label, date, transform, floor, t_otsu, g_so, s_so):
    """MNDWI of row 036 vs row 037 where both are clear, on the mosaic grid."""
    h, w = floor.shape
    out = {}
    per_row = {}
    for r in ("036", "037"):
        sc = meta["scenes"][f"{label}_r{r}"]
        b = sc["bands"]
        g = mo.read_on_grid(ROOT / b["green"]["raw_file"], transform, w, h, 0)
        s = mo.read_on_grid(ROOT / b["swir1"]["raw_file"], transform, w, h, 0)
        qa = mo.read_on_grid(ROOT / b["qa_pixel"]["raw_file"], transform, w, h, 1)
        m, _, _, bad = mndwi_from_dn(g, s, g_so, s_so)
        ok = ~bit(qa, FILL) & (g != 0) & (s != 0) & ~bad
        for f in EXCLUDED_FLAGS.values():
            ok &= ~bit(qa, f)
        per_row[r] = (m, ok)
    m36, ok36 = per_row["036"]
    m37, ok37 = per_row["037"]
    both = ok36 & ok37
    for name, area in (("all_overlap", both), ("valley_floor_overlap", both & floor)):
        a, b_ = m36[area], m37[area]
        cls_diff = (a > t_otsu) != (b_ > t_otsu)
        out[name] = {"mndwi_diff_036_minus_037": stats(a - b_),
                     "pearson_r": round(float(np.corrcoef(a, b_)[0, 1]), 6) if a.size > 1 else None,
                     "class_differs_fraction": round(float(cls_diff.mean()), 6) if a.size else None,
                     "class_differs_pixels": int(cls_diff.sum()),
                     "water_fraction_row036": round(float((a > t_otsu).mean()), 6) if a.size else None,
                     "water_fraction_row037": round(float((b_ > t_otsu).mean()), 6) if a.size else None}
    print(f"  seam {date}: {out['all_overlap']['mndwi_diff_036_minus_037']} "
          f"class differs {out['all_overlap']['class_differs_fraction']}")
    return out


if __name__ == "__main__":
    main()
