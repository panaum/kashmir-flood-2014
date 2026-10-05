"""Step 7a: draw the blinded validation sample and make the labelling material.

Implements rules S, P and B of scripts/step7_rules.md: seven strata over the
assessed analysis area, 56 samples per stratum (392), seed 20140910, at least
300 m apart, shuffled and numbered V001-V392. The answer key (stratum, map
class) goes to data/validation/key.csv, which is not committed and not shown
to the interpreter; only its sha256 is committed.

The labelling material shows only sample ids, coordinates and image chips
made from the Landsat reflectance bands. A check scans it for any
classification information and stops if it finds some. No validation metric
is computed here.

Outputs:
  data/validation/key.csv                      hidden answer key (not committed)
  data/validation/labelling/                   folder the interpreter opens
      index.html, samples.js, labels_template.csv, README.txt, chips/*.png
  outputs/validation/samples_blind.csv         id, lat, lon (committed)
  outputs/validation/labels_template.csv       empty labels (committed)
  outputs/validation/key_sha256.txt            checksum of the hidden key (committed)

Usage:  .venv/bin/python scripts/25_validation_sample.py
"""

import csv
import hashlib
import json
import re
import shutil
import sys
from pathlib import Path

import numpy as np
import rasterio
from PIL import Image
from rasterio.warp import transform as warp_transform
from scipy import ndimage

ROOT = Path(__file__).resolve().parent.parent
CLIP = ROOT / "data" / "clipped"
VAL = ROOT / "data" / "validation"
LAB = VAL / "labelling"
OUT = ROOT / "outputs" / "validation"
SEED = 20140910                       # P2
PER_STRATUM = 56                      # P1
MIN_SPACING_M = 300.0                 # P2
EDGE_M = 500.0                        # S6
BUILT_WINDOW, BUILT_MIN = 5, 10       # built-up proxy
DETAIL, CONTEXT = 51, 201             # B2, pixels
STRATA = ["S1", "S2", "S3", "S4", "S5", "S6", "S7"]
DATES = {"pre": "2014-08-25", "flood": "2014-09-10"}
FORBIDDEN = ["stratum", "strata", "mndwi", "flood", "class", "osm", "openstreetmap",
             "building", "footprint", "threshold", "persistent", "water mask", "map label"]
SURFACES = ["open water", "water with vegetation", "built-up", "crops or grass", "trees",
            "bare soil", "wetland or marsh", "cloud or shadow"]


def stop(msg):
    sys.exit(f"STOP: {msg}")


def read(path):
    with rasterio.open(path) as s:
        return s.read(1), s.transform, s.crs


def strata_map(cls, area, osm, built, near):
    """Stratum index 1-7 per pixel (first match wins), 0 = not in the population."""
    assessed = area & np.isin(cls, (0, 1, 2, 3))
    st = np.zeros(cls.shape, np.uint8)
    rules = [(1, (cls == 1) & built), (2, cls == 1), (3, np.isin(cls, (2, 3))),
             (4, (cls == 0) & osm), (5, (cls == 0) & built), (6, (cls == 0) & near), (7, cls == 0)]
    for k, m in rules:
        st[assessed & (st == 0) & m] = k
    return st, assessed


def stretch_limits(band, label_dates):
    """2nd-98th percentile of valid pixels, both dates pooled (step 1 rule)."""
    vals = []
    for label, date in label_dates.items():
        b = rasterio.open(CLIP / f"{label}_{date}_{band}.tif").read(1)
        qa = rasterio.open(CLIP / f"{label}_{date}_qa_pixel.tif").read(1)
        vals.append(b[(qa & 1) == 0].astype(np.float64))
    return [float(x) for x in np.percentile(np.concatenate(vals), [2, 98])]


def chip(bands, lims, r, c, size):
    h = size // 2
    H, W = bands[0].shape
    out = np.full((size, size, 3), 128, np.uint8)
    r0, r1, c0, c1 = max(r - h, 0), min(r + h + 1, H), max(c - h, 0), min(c + h + 1, W)
    for k, (b, (lo, hi)) in enumerate(zip(bands, lims)):
        v = np.clip((b[r0:r1, c0:c1].astype(np.float32) - lo) / (hi - lo), 0, 1)
        out[r0 - (r - h):r1 - (r - h), c0 - (c - h):c1 - (c - h), k] = (v * 255).round().astype(np.uint8)
    return Image.fromarray(out)


def visible_text(html):
    html = re.sub(r"<script.*?</script>|<style.*?</style>", " ", html, flags=re.S)
    return re.sub(r"<[^>]+>", " ", html)


def main():
    rules = (Path(__file__).parent / "step7_rules.md").read_text()
    if "20140910" not in rules or "56 samples per stratum" not in rules:
        stop("step7_rules.md does not match this script")
    fr = json.loads((ROOT / "data" / "flood_report.json").read_text())
    cls, T, crs = read(ROOT / "data" / "flood" / "flood_classes.tif")
    area = read(ROOT / "data" / "flood" / "analysis_area.tif")[0].astype(bool)
    osm = read(ROOT / "data" / "osm" / "water_wetland_mask.tif")[0].astype(bool)
    if round(int(area.sum()) * 0.0009, 4) != fr["analysis_area_km2"]:
        stop("analysis area differs from step 3")
    H, W = cls.shape

    # Built-up proxy: footprint centroids per 5 x 5 window (stratification only).
    cnt = np.zeros((H, W), np.int32)
    with open(ROOT / "data" / "msbuildings" / "footprints.csv") as f:
        for row in csv.DictReader(f):
            cc, rr = ~T * (float(row["centroid_x"]), float(row["centroid_y"]))
            cc, rr = int(np.floor(cc)), int(np.floor(rr))
            if 0 <= rr < H and 0 <= cc < W:
                cnt[rr, cc] += 1
    win = ndimage.convolve(cnt, np.ones((BUILT_WINDOW, BUILT_WINDOW), np.int32), mode="constant")
    built = win >= BUILT_MIN
    near = ndimage.distance_transform_edt(cls != 1) * T.a <= EDGE_M
    st, assessed = strata_map(cls, area, osm, built, near)

    # Sampling (P2): strata in order, random candidate order, 300 m spacing.
    rng = np.random.default_rng(SEED)
    cell = MIN_SPACING_M
    grid = {}
    accepted = []

    def far_enough(x, y):
        gx, gy = int(x // cell), int(y // cell)
        for i in (-1, 0, 1):
            for j in (-1, 0, 1):
                for (px, py) in grid.get((gx + i, gy + j), ()):
                    if (px - x) ** 2 + (py - y) ** 2 < MIN_SPACING_M ** 2:
                        return False
        return True

    stratum_km2 = {}
    for k, name in enumerate(STRATA, start=1):
        rr, cc = np.nonzero(st == k)
        stratum_km2[name] = round(rr.size * 0.0009, 4)
        order = rng.permutation(rr.size)
        got = 0
        for i in order:
            x, y = T * (cc[i] + 0.5, rr[i] + 0.5)
            if far_enough(x, y):
                grid.setdefault((int(x // cell), int(y // cell)), []).append((x, y))
                accepted.append({"stratum": name, "row": int(rr[i]), "col": int(cc[i]),
                                 "x": x, "y": y, "map_class": int(cls[rr[i], cc[i]])})
                got += 1
                if got == PER_STRATUM:
                    break
        if got < PER_STRATUM:
            stop(f"{name}: only {got} samples at {MIN_SPACING_M:g} m spacing")
    if round(sum(stratum_km2.values()), 4) != fr["valid_both_km2"]:
        stop("strata do not cover the assessed area exactly")

    # P3: shuffle and number.
    order = rng.permutation(len(accepted))
    samples = [accepted[i] for i in order]
    lon, lat = warp_transform(crs, "EPSG:4326", [s["x"] for s in samples], [s["y"] for s in samples])
    for n, (s, lo, la) in enumerate(zip(samples, lon, lat), start=1):
        s["id"] = f"V{n:03d}"
        s["lat"], s["lon"] = round(la, 6), round(lo, 6)

    # Hidden key (P4).
    VAL.mkdir(parents=True, exist_ok=True)
    key_path = VAL / "key.csv"
    with open(key_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sample_id", "stratum", "map_class", "row", "col", "x_utm43n", "y_utm43n"])
        for s in sorted(samples, key=lambda s: s["id"]):
            w.writerow([s["id"], s["stratum"], s["map_class"], s["row"], s["col"],
                        round(s["x"], 2), round(s["y"], 2)])
    key_sha = hashlib.sha256(key_path.read_bytes()).hexdigest()

    # Chips (B2) from reflectance bands only.
    if LAB.exists():
        shutil.rmtree(LAB)
    (LAB / "chips").mkdir(parents=True)
    lims = {b: stretch_limits(b, DATES) for b in ("swir1", "nir", "red", "green")}
    step1 = json.loads((ROOT / "previews" / "qa_summary.json").read_text())["stretch_2_98_percentile_raw_DN"]
    for b in ("swir1", "nir", "red"):
        if [round(v, 6) for v in lims[b]] != [round(v, 6) for v in step1[b]]:
            stop(f"{b} stretch differs from step 1")
    bands = {(lab, b): rasterio.open(CLIP / f"{lab}_{d}_{b}.tif").read(1)
             for lab, d in DATES.items() for b in ("swir1", "nir", "red", "green")}
    tag = {"pre": "aug", "flood": "sep"}
    for s in samples:
        for lab in DATES:
            swir = [bands[(lab, b)] for b in ("swir1", "nir", "red")]
            cir = [bands[(lab, b)] for b in ("nir", "red", "green")]
            l_swir = [lims[b] for b in ("swir1", "nir", "red")]
            l_cir = [lims[b] for b in ("nir", "red", "green")]
            chip(swir, l_swir, s["row"], s["col"], DETAIL).save(LAB / "chips" / f"{s['id']}_{tag[lab]}_a.png")
            chip(cir, l_cir, s["row"], s["col"], DETAIL).save(LAB / "chips" / f"{s['id']}_{tag[lab]}_b.png")
            chip(swir, l_swir, s["row"], s["col"], CONTEXT).save(LAB / "chips" / f"{s['id']}_{tag[lab]}_c.png")

    # Labelling files: ids, coordinates and chips only.
    blind = [{"id": s["id"], "lat": s["lat"], "lon": s["lon"]} for s in sorted(samples, key=lambda s: s["id"])]
    (LAB / "samples.js").write_text("window.SAMPLES = " + json.dumps(blind) + ";\n")
    cols = ["sample_id", "lat", "lon", "aug_water", "aug_surface", "sep_water", "sep_surface",
            "confidence", "other_source", "notes"]
    OUT.mkdir(parents=True, exist_ok=True)
    for path in (LAB / "labels_template.csv", OUT / "labels_template.csv"):
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(cols)
            for b in blind:
                w.writerow([b["id"], b["lat"], b["lon"]] + [""] * (len(cols) - 3))
    with open(OUT / "samples_blind.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sample_id", "lat", "lon"])
        for b in blind:
            w.writerow([b["id"], b["lat"], b["lon"]])
    (OUT / "key_sha256.txt").write_text(f"{key_sha}  data/validation/key.csv\n")
    shutil.copy(Path(__file__).parent / "validation_page.html", LAB / "index.html")
    (LAB / "README.txt").write_text(
        "Water labelling, Landsat 8, 25 August 2014 and 10 September 2014\n\n"
        "Open index.html in a web browser. It works offline from this folder.\n"
        "For each sample, look at the centre pixel (outlined) on each date and record whether it is\n"
        "water, what surface you see, and how confident you are. Use 'cannot tell' whenever the\n"
        "image does not let you decide. Your answers are kept in the browser as you go; use\n"
        "'Download CSV' to save them, and 'Load CSV' to continue later.\n\n"
        "Chips: a = shortwave infrared / near infrared / red; b = near infrared / red / green;\n"
        "c = wider view (6 km) in shortwave infrared / near infrared / red. 30 m pixels.\n"
        "Do not consult any project map or result while labelling.\n")

    # B3: blinding check on everything in the labelling folder.
    found = []
    for path in [LAB / "samples.js", LAB / "labels_template.csv", LAB / "README.txt",
                 OUT / "samples_blind.csv", OUT / "labels_template.csv"]:
        text = path.read_text().lower()
        found += [(path.name, t) for t in FORBIDDEN if t in text]
    vis = visible_text((LAB / "index.html").read_text()).lower()
    found += [("index.html (visible text)", t) for t in FORBIDDEN if t in vis]
    js = re.findall(r"<script>(.*?)</script>", (LAB / "index.html").read_text(), re.S)
    if not js:
        stop("no inline script found in the labelling page to check")
    strings = " ".join(a or b for a, b in re.findall(r"\"([^\"]*)\"|'([^']*)'", " ".join(js)))
    found += [("index.html (script strings)", t) for t in FORBIDDEN if t in strings.lower()]
    names = " ".join(p.name for p in (LAB / "chips").iterdir()).lower()
    found += [("chip file names", t) for t in FORBIDDEN if t in names]
    if found:
        stop(f"classification information found in the labelling material: {found}")

    # Checks on the sample (item 4 of the request; rules P).
    key = list(csv.DictReader(open(key_path)))
    ids = [k["sample_id"] for k in key]
    if len(set(ids)) != 392 or set(ids) != {b["id"] for b in blind}:
        stop("sample ids do not join one-to-one with the key")
    for k in key:
        r, c = int(k["row"]), int(k["col"])
        if not assessed[r, c] or int(cls[r, c]) not in (0, 1, 2, 3) or int(k["map_class"]) != int(cls[r, c]):
            stop(f"{k['sample_id']} is not on an assessed pixel")
    xy = np.array([[float(k["x_utm43n"]), float(k["y_utm43n"])] for k in key])
    d = np.sqrt(((xy[:, None, :] - xy[None, :, :]) ** 2).sum(-1)) + np.eye(len(xy)) * 1e9
    counts = {}
    for k in key:
        counts[k["stratum"]] = counts.get(k["stratum"], 0) + 1
    if any(counts.get(s) != PER_STRATUM for s in STRATA):
        stop(f"stratum counts {counts} are not {PER_STRATUM} each")
    n_chips = len(list((LAB / "chips").iterdir()))
    print(json.dumps({"samples": len(key), "per_stratum": counts, "stratum_km2": stratum_km2,
                      "min_spacing_m": round(float(d.min()), 1), "all_on_assessed_pixels": True,
                      "ids_join_key": True, "chips": n_chips, "blinding_check": "passed",
                      "key_sha256": key_sha, "stretch_limits": lims}, indent=1))


if __name__ == "__main__":
    main()
