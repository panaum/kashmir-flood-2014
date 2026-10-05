"""Step 5a: continuous threshold sensitivity of the flood footprint.

Implements rules T1-T7 of scripts/step5_rules.md. The step 3 functions in
10_flood.py (load_date, water, outcome) are used unchanged: only the MNDWI
threshold varies, applied to both dates. Everything else (inputs, masks,
analysis area, flood definition) is exactly step 3.

The step 3 result (376.6077 km2 at 0.127) stays the primary result; this
sweep only describes how the footprint changes with the threshold.

Outputs: outputs/threshold_sweep.csv, outputs/threshold_sweep.json,
outputs/figures/threshold_sweep.png

Usage:  .venv/bin/python scripts/21_threshold_sweep.py
"""

import csv
import importlib.util
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "outputs"
GRID = [round(k / 100, 2) for k in range(0, 31)]          # T2
WINDOW = (0.08, 0.17)                                       # T6: grid points within +/- 0.05 of 0.127
SHADE = (0.077, 0.177)                                      # T7: step 3 +/- 0.05 variants
BASELINE_KM2 = 376.6077                                     # T3

_spec = importlib.util.spec_from_file_location("fl", Path(__file__).parent / "10_flood.py")
fl = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fl)


def main():
    (OUT / "figures").mkdir(parents=True, exist_ok=True)
    meta = json.loads((ROOT / "data" / "metadata.json").read_text())
    wr = json.loads((ROOT / "data" / "water_report.json").read_text())
    fr = json.loads((ROOT / "data" / "flood_report.json").read_text())

    area = rasterio.open(ROOT / "data" / "flood" / "analysis_area.tif").read(1).astype(bool)
    if fl.km2(area) != fr["analysis_area_km2"]:
        sys.exit("STOP: analysis_area.tif differs from the step 3 analysis area")
    pre = fl.load_date(meta, "pre", fl.PRE, area)
    post = fl.load_date(meta, "flood", fl.POST, area)
    both = pre["valid"] & post["valid"]
    if fl.km2(both) != fr["valid_both_km2"]:
        sys.exit("STOP: valid-on-both-dates area differs from step 3")

    def run(t):
        res, wp, wq = fl.outcome(pre, post, t, t, both)
        n = int(np.count_nonzero(wq & ~wp))
        return {"threshold": t, "flooded_pixels": n, "flooded_km2": res["flooded_km2"],
                "flooded_pct_of_analysis_area": round(100 * n / int(area.sum()), 4),
                "persistent_km2": res["persistent_km2"], "lost_km2": res["lost_km2"],
                "water_aug_km2": res["water_pre_km2"], "water_sep_km2": res["water_post_km2"]}

    t0 = wr["dates"][fl.POST]["thresholds"]["otsu"]
    base = run(t0)
    if base["flooded_km2"] != BASELINE_KM2:
        sys.exit(f"STOP: threshold {t0} gives {base['flooded_km2']} km2, not the step 3 {BASELINE_KM2}")
    rows = [run(t) for t in GRID]
    print(f"baseline {t0}: {base['flooded_km2']} km2 (reproduced)")

    areas = np.array([r["flooded_km2"] for r in rows])
    win = [r for r in rows if WINDOW[0] <= r["threshold"] <= WINDOW[1]]
    wa = np.array([r["flooded_km2"] for r in win])
    steps = [{"from": a["threshold"], "to": b["threshold"],
              "change_km2": round(b["flooded_km2"] - a["flooded_km2"], 4)}
             for a, b in zip(rows, rows[1:])]
    summary = {
        "grid": {"start": GRID[0], "stop": GRID[-1], "step": 0.01, "n": len(GRID)},
        "baseline": base,
        "grid_min": rows[int(areas.argmin())], "grid_max": rows[int(areas.argmax())],
        "window": {"range": list(WINDOW), "n": len(win),
                   "min_km2": float(wa.min()), "max_km2": float(wa.max()),
                   "min_pct_of_baseline": round(100 * wa.min() / base["flooded_km2"], 2),
                   "max_pct_of_baseline": round(100 * wa.max() / base["flooded_km2"], 2)},
        "largest_step_change": max(steps, key=lambda s: abs(s["change_km2"])),
        "step_changes": steps,
        "monotonic_decreasing": bool(np.all(np.diff(areas) <= 0)),
    }

    with open(OUT / "threshold_sweep.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    (OUT / "threshold_sweep.json").write_text(json.dumps(
        {"script": "scripts/21_threshold_sweep.py", "rules": "scripts/step5_rules.md (A)",
         "flood_definition": "water on 2014-09-10 and not water on 2014-08-25 (MNDWI > t on both "
                             "dates), pixels valid on both dates, step 3 analysis area",
         "summary": summary, "rows": rows}, indent=1))

    fig, ax = plt.subplots(figsize=(9, 5.2))
    ax.axvspan(*SHADE, facecolor="0.88", edgecolor="0.55", lw=0.8, zorder=0, label="step 3 variants, 0.077 to 0.177")
    ax.plot([r["threshold"] for r in rows], areas, "-o", color="#1F5E59", ms=3.5, lw=1.4,
            label="flooded area, same threshold on both dates")
    ax.plot([t0], [base["flooded_km2"]], "D", color="#8A3A1B", ms=8, zorder=5,
            label=f"primary result: {t0:g} -> {base['flooded_km2']:.1f} km2")
    ax.axvline(t0, color="#8A3A1B", lw=0.8, ls="--")
    ax.set_xlabel("MNDWI water threshold")
    ax.set_ylabel("flooded area (km2)")
    ax.set_xlim(-0.005, 0.305)
    ax.set_title("Detected flood footprint, 10 September 2014, against the water threshold")
    ax.grid(alpha=0.3)
    ax.legend(frameon=False, loc="upper right")
    fig.savefig(OUT / "figures" / "threshold_sweep.png", dpi=150, bbox_inches="tight")
    print(json.dumps({k: summary[k] for k in ("grid_min", "grid_max", "window", "largest_step_change",
                                             "monotonic_decreasing")}, indent=1))
    for r in rows:
        print(f"{r['threshold']:.2f}  {r['flooded_km2']:9.4f} km2  {r['flooded_pct_of_analysis_area']:.3f} %")


if __name__ == "__main__":
    main()
