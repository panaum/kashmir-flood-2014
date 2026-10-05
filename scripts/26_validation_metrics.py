"""Step 7b: validation metrics from the frozen reference labels (rule M).

Implements rules C and M of scripts/step7_rules.md with amendments 2 and 3,
using only the frozen labels outputs/validation/labels_v1_final.csv (checked
against its committed sha256) and the hidden v1 key. The reference set is
the completed rows (August water, September water and confidence present);
V005 and every other incomplete row are excluded. The 140-point v2 sample is
not read. Reference labels come only from the human file; map labels come
only from the key.

Estimators (Stehman 2014 stratified estimator, strata differ from map
classes): W_h = stratum area / total assessed area; Y = sum W_h mean_h(y),
V(Y) = sum W_h^2 s_h^2(y) / n_h; ratios R = Y/X with the linearised variance;
F1 = Y/X with y = 2 [map and ref], x = [map] + [ref]; sample variances with
n_h - 1; 95 % interval = estimate +/- 1.96 SE, not clipped, out-of-range
bounds flagged.

The reference is image-based interpretation of the same Landsat scenes, not
field truth. The reference-based area is secondary to the 376.6077 km2
result and never replaces it.

Outputs: outputs/validation/metrics_v1.json, outputs/validation/metrics_v1_summary.csv

Usage:  .venv/bin/python scripts/26_validation_metrics.py
"""

import csv
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parent.parent
VOUT = ROOT / "outputs" / "validation"
LABELS = VOUT / "labels_v1_final.csv"
KEY = ROOT / "data" / "validation" / "key.csv"
STRATA = ["S1", "S2", "S3", "S4", "S5", "S6", "S7"]
EXPECTED_N = {"S1": 18, "S2": 9, "S3": 21, "S4": 16, "S5": 19, "S6": 17, "S7": 16}
Z = 1.96
PIXEL_KM2 = 0.0009
PRIMARY_KM2 = 376.6077
SMALL_N = 10

_spec = importlib.util.spec_from_file_location("vs", Path(__file__).parent / "25_validation_sample.py")
vs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(vs)


def stop(msg):
    sys.exit(f"STOP: {msg}")


def stratum_areas():
    """Rule S strata on the step 3 rasters, exactly as the sampling script built them."""
    import csv as _csv
    from scipy import ndimage
    cls, T, _ = vs.read(ROOT / "data" / "flood" / "flood_classes.tif")
    area = vs.read(ROOT / "data" / "flood" / "analysis_area.tif")[0].astype(bool)
    osm = vs.read(ROOT / "data" / "osm" / "water_wetland_mask.tif")[0].astype(bool)
    H, W = cls.shape
    cnt = np.zeros((H, W), np.int32)
    with open(ROOT / "data" / "msbuildings" / "footprints.csv") as f:
        for row in _csv.DictReader(f):
            cc, rr = ~T * (float(row["centroid_x"]), float(row["centroid_y"]))
            cc, rr = int(np.floor(cc)), int(np.floor(rr))
            if 0 <= rr < H and 0 <= cc < W:
                cnt[rr, cc] += 1
    win = ndimage.convolve(cnt, np.ones((vs.BUILT_WINDOW, vs.BUILT_WINDOW), np.int32), mode="constant")
    near = ndimage.distance_transform_edt(cls != 1) * T.a <= vs.EDGE_M
    st, assessed = vs.strata_map(cls, area, osm, win >= vs.BUILT_MIN, near)
    return st, assessed, cls


def strat_mean(groups, W, key):
    """Y and V(Y) for indicator `key` over strata (M1)."""
    Y = V = 0.0
    for h in STRATA:
        y = np.array([s[key] for s in groups[h]], float)
        n = y.size
        if n < 2:
            stop(f"{h}: fewer than 2 usable samples for {key}")
        Y += W[h] * y.mean()
        V += W[h] ** 2 * y.var(ddof=1) / n
    return Y, V


def strat_ratio(groups, W, ykey, xkey):
    """R = Y/X and its linearised variance (M2, amendment 3)."""
    Y, _ = strat_mean(groups, W, ykey)
    X, _ = strat_mean(groups, W, xkey)
    if X == 0:
        return None, None
    R = Y / X
    V = 0.0
    for h in STRATA:
        y = np.array([s[ykey] for s in groups[h]], float)
        x = np.array([s[xkey] for s in groups[h]], float)
        n = y.size
        cov = np.cov(x, y, ddof=1)[0, 1]
        V += W[h] ** 2 * (y.var(ddof=1) + R ** 2 * x.var(ddof=1) - 2 * R * cov) / n
    return R, V / X ** 2


def interval(est, var):
    se = float(np.sqrt(max(var, 0.0)))
    lo, hi = est - Z * se, est + Z * se
    return {"estimate": round(float(est), 6), "se": round(se, 6),
            "ci95": [round(float(lo), 6), round(float(hi), 6)],
            "bound_outside_0_1": bool(lo < 0 or hi > 1)}


def evaluate(samples, W, target):
    """Rule M3 for one target: 'aug', 'sep' or 'flood'."""
    groups = {h: [] for h in STRATA}
    for s in samples:
        m, r = s[f"map_{target}"], s[f"ref_{target}"]
        if r is None:                      # cannot tell (rule C3)
            continue
        groups[s["stratum"]].append({
            "agree": int(m == r), "map": int(m), "ref": int(r), "tp": int(m and r),
            "c11": int(m and r), "c10": int(m and not r), "c01": int((not m) and r), "c00": int((not m) and (not r)),
            "f1y": 2 * int(m and r), "f1x": int(m) + int(r)})
    n_used = {h: len(groups[h]) for h in STRATA}
    out = {"n_used_by_stratum": n_used, "n_used": sum(n_used.values()),
           "excluded_cannot_tell": len(samples) - sum(n_used.values())}
    out["error_matrix_proportions"] = {
        name: interval(*strat_mean(groups, W, k))
        for name, k in (("map_positive_ref_positive", "c11"), ("map_positive_ref_negative", "c10"),
                        ("map_negative_ref_positive", "c01"), ("map_negative_ref_negative", "c00"))}
    out["overall_accuracy"] = interval(*strat_mean(groups, W, "agree"))
    p, vp = strat_ratio(groups, W, "tp", "map")
    r, vr = strat_ratio(groups, W, "tp", "ref")
    f, vf = strat_ratio(groups, W, "f1y", "f1x")
    out["precision_users_accuracy"] = interval(p, vp) if p is not None else None
    out["recall_producers_accuracy"] = interval(r, vr) if r is not None else None
    out["f1"] = interval(f, vf) if f is not None else None
    if p is not None and r is not None and f is not None and abs(f - 2 * p * r / (p + r)) > 1e-9:
        stop(f"{target}: F1 ratio form does not equal 2PR/(P+R)")
    out["map_positive_proportion"] = interval(*strat_mean(groups, W, "map"))
    out["ref_positive_proportion"] = interval(*strat_mean(groups, W, "ref"))
    out["raw_counts_by_stratum"] = {
        h: {k: int(sum(s[k] for s in groups[h])) for k in ("c11", "c10", "c01", "c00")} for h in STRATA}
    out["zero_variance_strata_for_accuracy"] = [h for h in STRATA
                                                if np.var([s["agree"] for s in groups[h]]) == 0]
    return out, groups


def main():
    # --- Inputs and safeguards -------------------------------------------------
    sha = hashlib.sha256(LABELS.read_bytes()).hexdigest()
    if sha != (VOUT / "labels_v1_final_sha256.txt").read_text().split()[0]:
        stop("frozen labels do not match their committed sha256")
    if hashlib.sha256(KEY.read_bytes()).hexdigest() != (VOUT / "key_sha256.txt").read_text().split()[0]:
        stop("v1 key does not match its committed sha256")
    rows = list(csv.DictReader(open(LABELS, newline="")))
    if [r["sample_id"] for r in rows] != [f"V{i:03d}" for i in range(1, 393)]:
        stop("labels file ids are not V001-V392")
    key = {k["sample_id"]: k for k in csv.DictReader(open(KEY))}
    complete = [r for r in rows if all(r[f].strip() for f in ("aug_water", "sep_water", "confidence"))]
    if any(r["sample_id"] == "V005" for r in complete):
        stop("V005 is in the reference set")
    counts = {h: sum(1 for r in complete if key[r["sample_id"]]["stratum"] == h) for h in STRATA}
    if counts != EXPECTED_N:
        stop(f"reference counts {counts} differ from the frozen record {EXPECTED_N}")

    st, assessed, cls = stratum_areas()
    area_px = {h: int((st == i).sum()) for i, h in enumerate(STRATA, start=1)}
    total_px = int(assessed.sum())
    if sum(area_px.values()) != total_px:
        stop("strata do not cover the assessed area")
    W = {h: area_px[h] / total_px for h in STRATA}
    total_km2 = round(total_px * PIXEL_KM2, 4)

    WATER = {"yes": True, "no": False, "cannot tell": None}
    samples = []
    for r in complete:
        k = key[r["sample_id"]]
        rr, cc, c = int(k["row"]), int(k["col"]), int(k["map_class"])
        if STRATA[st[rr, cc] - 1] != k["stratum"] or int(cls[rr, cc]) != c:
            stop(f"{r['sample_id']}: key does not match the rasters")
        a, s_ = WATER[r["aug_water"]], WATER[r["sep_water"]]
        samples.append({
            "id": r["sample_id"], "stratum": k["stratum"], "confidence": r["confidence"],
            "map_aug": c in (2, 3), "map_sep": c in (1, 2), "map_flood": c == 1,
            "ref_aug": a, "ref_sep": s_,
            "ref_flood": None if (a is None or s_ is None) else (s_ and not a)})

    # --- Rule C: cannot tell ---------------------------------------------------
    ct = {"aug": sum(s["ref_aug"] is None for s in samples),
          "sep": sum(s["ref_sep"] is None for s in samples)}
    ct_by = {h: {d: sum(1 for s in samples if s["stratum"] == h and s[f"ref_{d}"] is None) for d in ("aug", "sep")}
             for h in STRATA}
    c4 = any(ct[d] > 0.05 * len(samples) for d in ct) or any(
        ct_by[h][d] > 0.10 * counts[h] for h in STRATA for d in ("aug", "sep"))
    if c4:
        stop("rule C4 sensitivity is triggered; not implemented in this run")

    # --- Rule M3, M4, M5 ---------------------------------------------------------
    results = {}
    for target in ("aug", "sep", "flood"):
        results[target], groups = evaluate(samples, W, target)
    Yf = results["flood"]["ref_positive_proportion"]
    area = {k: round(v * total_km2, 4) for k, v in (("estimate", Yf["estimate"]),
                                                    ("se", Yf["se"]))}
    area["ci95"] = [round(Yf["ci95"][0] * total_km2, 4), round(Yf["ci95"][1] * total_km2, 4)]
    area["bound_below_0"] = area["ci95"][0] < 0
    # Consistency check on the unrounded estimate: S1 and S2 are exactly the
    # flooded pixels, so the sample must reproduce the primary area.
    y_map, _ = strat_mean(groups, W, "map")
    map_area = round(y_map * total_px * PIXEL_KM2, 4)
    if abs(map_area - PRIMARY_KM2) > 0.001:
        stop(f"map flood proportion from the sample ({map_area} km2) does not reproduce the primary result")

    # --- Rule M6: without low-confidence labels ---------------------------------
    low = [s for s in samples if s["confidence"] == "low"]
    m6 = {"low_confidence_rows": len(low)}
    if low:
        stop("low-confidence labels present; M6 sensitivity run not implemented in this run")
    m6["note"] = "No low-confidence labels; the sensitivity without them is identical to the primary metrics."

    report = {
        "script": "scripts/26_validation_metrics.py",
        "rules": "scripts/step7_rules.md, rules C and M with amendments 2 and 3",
        "nature_of_reference": "Image-based reference interpretation of the same Landsat 8 scenes by a human "
                               "interpreter who knows the project; not field truth or independent ground truth.",
        "inputs": {"labels": "outputs/validation/labels_v1_final.csv", "labels_sha256": sha,
                   "key": "data/validation/key.csv (not committed; sha256 in outputs/validation/key_sha256.txt)"},
        "reference_rows": len(samples), "reference_rows_by_stratum": counts,
        "excluded": {"incomplete_rows": len(rows) - len(complete),
                     "V005": "stray confidence value, no water answers; not in the reference set (freeze record)",
                     "v2_sample": "not read"},
        "stratum_area_km2": {h: round(area_px[h] * PIXEL_KM2, 4) for h in STRATA},
        "stratum_weights": {h: round(W[h], 6) for h in STRATA},
        "total_assessed_km2": total_km2,
        "cannot_tell": {"by_date": ct, "by_stratum": ct_by, "c4_flag": c4},
        "small_strata_n_lt_10": [h for h in STRATA if counts[h] < SMALL_N],
        "water_2014_08_25": results["aug"], "water_2014_09_10": results["sep"], "flood": results["flood"],
        "reference_based_flooded_area_km2": {**area, "primary_result_km2": PRIMARY_KM2,
                                             "note": "secondary; does not replace the primary result"},
        "map_flood_area_from_sample_km2": map_area,
        "m6_low_confidence": m6,
    }
    (VOUT / "metrics_v1.json").write_text(json.dumps(report, indent=1))
    with open(VOUT / "metrics_v1_summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["target", "metric", "estimate", "se", "ci95_low", "ci95_high", "bound_outside_0_1", "n_used"])
        for t, name in (("aug", "water 2014-08-25"), ("sep", "water 2014-09-10"), ("flood", "flood")):
            for m in ("overall_accuracy", "precision_users_accuracy", "recall_producers_accuracy", "f1",
                      "map_positive_proportion", "ref_positive_proportion"):
                v = results[t][m]
                if v:
                    w.writerow([name, m, v["estimate"], v["se"], v["ci95"][0], v["ci95"][1],
                                v["bound_outside_0_1"], results[t]["n_used"]])
        w.writerow(["flood", "reference_based_flooded_area_km2", area["estimate"], area["se"],
                    area["ci95"][0], area["ci95"][1], area["bound_below_0"], results["flood"]["n_used"]])
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
