"""Step 4b-b: count Microsoft footprints on flooded ground (rules M3-M8).

Reads the parquet files recorded by 15_msbuildings_download.py and applies the
step 4 counting rules unchanged (B3 centroid, B4 any-overlap, S1 Srinagar
window, W1 step 3 water/wetland mask), using the same functions as
13_exposure.py. Adds the Srinagar 1 km grid (M5), the OSM cross-check (M6)
and the completeness test (M7).

Every figure is "footprints in the Microsoft data (imagery 2017-2021 per the
catalogue) standing on ground flooded on 2014-09-10" (M8).

Outputs: data/msbuildings_report.json, data/msbuildings/footprints.csv.

Usage:  .venv/bin/python scripts/16_msbuildings_count.py
"""

import csv
import importlib.util
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
import rasterio
import rasterio.windows
import shapely
from rasterio.features import shapes
from shapely.geometry import box

ROOT = Path(__file__).resolve().parent.parent
MS_DIR = ROOT / "data" / "msbuildings"
CLIP_BOX = (74.2, 33.5, 75.4, 34.6)
GRID_M = 1000.0
SHARE_MIN_ROAD_KM = 1.0          # M5
MATCH_SHARE = 0.5                # M6

_spec = importlib.util.spec_from_file_location("ex", Path(__file__).parent / "13_exposure.py")
ex = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ex)
FLOODED, NOT_ASSESSED, OUTSIDE = ex.FLOODED, ex.NOT_ASSESSED, ex.OUTSIDE


def load_footprints(record):
    """All rows -> unique WKB -> shapely geometries in lon/lat (M3)."""
    seen, total, crs_meta = {}, 0, set()
    for it in record["items"]:
        for f in it["files"]:
            pf = pq.ParquetFile(ROOT / f["local"])
            geo = (pf.schema_arrow.metadata or {}).get(b"geo")
            if geo:
                g = json.loads(geo)
                col = g.get("primary_column", "geometry")
                crs_meta.add(json.dumps(g["columns"][col].get("crs", "absent (GeoParquet default OGC:CRS84)")))
            for wkb in pf.read(columns=["geometry"]).column("geometry").to_pylist():
                total += 1
                seen.setdefault(wkb, None)
    wkbs = list(seen)
    geoms = shapely.from_wkb(wkbs, on_invalid="ignore")
    bad = np.array([g is None or g.is_empty for g in geoms])
    return np.asarray(geoms, dtype=object)[~bad], {"rows": total, "unique_geometries": len(wkbs),
                                                   "duplicates_removed": total - len(wkbs),
                                                   "unparseable_or_empty": int(bad.sum()),
                                                   "geoparquet_crs": sorted(crs_meta)}


def main():
    record = json.loads((MS_DIR / "record.json").read_text())
    geoms, load_stats = load_footprints(record)
    crs_set = load_stats["geoparquet_crs"]
    if any("4326" not in c and "CRS84" not in c and "absent" not in c for c in crs_set):
        sys.exit(f"STOP: unexpected footprint CRS {crs_set}")
    # Keep footprints whose bounds touch the clip box before projecting.
    b = shapely.bounds(geoms)
    near = (b[:, 2] >= CLIP_BOX[0]) & (b[:, 0] <= CLIP_BOX[2]) & (b[:, 3] >= CLIP_BOX[1]) & (b[:, 1] <= CLIP_BOX[3])
    geoms = geoms[near]
    print(f"{load_stats}; {len(geoms)} footprints near the clip box")

    state, code, transform, crs, shp, sb, win = ex.state_rasters()
    utm = ex.to_utm(list(geoms), crs)
    cent, cx, cy, b_code, b_state, ov_state = ex.classify_buildings(utm, code, state, transform, shp)
    in_osm, in_srin = (b_code & 4) > 0, (b_code & 8) > 0
    every = np.ones(len(utm), bool)
    counts = {
        "centroid": {"analysis_area": ex.bcount(every, b_state),
                     "srinagar_box": ex.bcount(in_srin, b_state),
                     "inside_water_wetland": ex.bcount(in_osm, b_state),
                     "outside_water_wetland": ex.bcount(~in_osm, b_state),
                     "srinagar_inside_water_wetland": ex.bcount(in_srin & in_osm, b_state),
                     "srinagar_outside_water_wetland": ex.bcount(in_srin & ~in_osm, b_state)},
        "overlap": {"analysis_area": ex.bcount(every, ov_state),
                    "srinagar_box": ex.bcount(in_srin, ov_state)}}
    area_m2 = shapely.area(utm)
    sel = b_state != OUTSIDE
    footprint_area = {"median_m2": round(float(np.median(area_m2[sel])), 1),
                      "p10_p90_m2": [round(float(x), 1) for x in np.percentile(area_m2[sel], [10, 90])]}

    # OSM roads and buildings (step 4 extract), projected the same way.
    _, elements = ex.load_elements()
    o_els, o_geoms, _, r_els, r_geoms, _ = ex.osm_geometries(elements)
    r_utm = ex.to_utm(r_geoms, crs)
    hw = np.array([e["tags"]["highway"] for e in r_els])
    motor = np.isin(hw, list(ex.MOTOR))
    o_utm = ex.to_utm(o_geoms, crs)
    o_cent, _, _, _, o_state, _ = ex.classify_buildings(o_utm, code, state, transform, shp)

    grid = srinagar_grid(state, transform, crs, win, cent, b_state, r_utm[motor], o_cent, o_state)
    cross = cross_check(grid["osm_top_cells"], o_els, o_utm, o_cent, o_state, utm)
    completeness = ex.completeness_check(transform, shp, state, cent, b_state, r_utm, motor, in_srin, sb)

    report = {"rules_file": "scripts/step4b_rules.md",
              "rules_sha256": __import__("hashlib").sha256(
                  (Path(__file__).parent / "step4b_rules.md").read_bytes()).hexdigest(),
              "source": {k: record[k] for k in ("stac_api", "collection", "collection_title", "license",
                                                "dataset_version_processing_date", "delta_table",
                                                "delta_log", "started_utc", "finished_utc", "tls_ca_bundle",
                                                "rows_total")},
              "items": [{k: i[k] for k in ("id", "region", "quadkey", "start_datetime", "end_datetime",
                                           "rows", "files_in_directory_not_in_log",
                                           "files_in_log_not_in_directory")} |
                        {"files": len(i["files"]), "bytes": sum(f["bytes"] for f in i["files"]),
                         "md5_verified_files": sum(f["md5_verified"] for f in i["files"])}
                        for i in record["items"]],
              "load": load_stats, "footprints_near_clip_box": int(len(geoms)),
              "footprint_area_in_analysis_area": footprint_area,
              "counts": counts, "grid": {k: v for k, v in grid.items() if k != "osm_top_cells"},
              "cross_check": cross, "completeness": completeness}
    (ROOT / "data" / "msbuildings_report.json").write_text(json.dumps(report, indent=2))
    with open(MS_DIR / "footprints.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["centroid_x", "centroid_y", "area_m2", "state_centroid", "state_overlap",
                    "in_water_wetland", "in_srinagar_box"])
        for x, y, a, s, o, wo, sr in zip(cx, cy, area_m2, b_state, ov_state, in_osm, in_srin):
            if s != OUTSIDE or o != OUTSIDE:
                w.writerow([round(x, 2), round(y, 2), round(a, 1), ex.STATE_NAMES[s],
                            ex.STATE_NAMES[o], int(wo), int(sr)])
    print(json.dumps({k: report[k] for k in ("counts", "grid", "cross_check", "completeness")}, indent=1)[:9000])


def place_points(crs):
    data = json.loads((ex.OSM_DIR / "srinagar_areas.json").read_text())
    places = []
    for el in data["elements"]:
        tags = el.get("tags", {})
        name = tags.get("name:en") or tags.get("name")
        if tags.get("place") in ("suburb", "neighbourhood", "quarter") and name:
            if el["type"] == "node":
                pt = ex.to_utm([shapely.Point(el["lon"], el["lat"])], crs)[0]
            else:
                poly = ex.osm3.element_polygon(el)
                if poly is None:
                    continue
                pt = shapely.centroid(ex.to_utm([poly], crs)[0])
            places.append((name, tags["place"], pt))
    if not places:
        sys.exit("STOP: no named places for grid labels")
    return places


def srinagar_grid(state, transform, crs, win, cent, b_state, motor_utm, o_cent, o_state):
    """1 km cells over the Srinagar window (M5)."""
    r0, c0 = int(win.row_off), int(win.col_off)
    wx0, wy1 = transform * (c0, r0)
    wx1, wy0 = transform * (c0 + int(win.width), r0 + int(win.height))
    window_box = box(wx0, wy0, wx1, wy1)
    places = place_points(crs)
    px = np.array([shapely.get_x(p) for *_, p in places])
    py = np.array([shapely.get_y(p) for *_, p in places])

    # Motor-road pieces inside the window, by pixel state (exact pixel outlines).
    sub = state[win.toslices()]
    sub_t = rasterio.windows.transform(win, transform)
    polys, vals = [], []
    for g, v in shapes(sub, transform=sub_t, connectivity=4):
        if int(v) != OUTSIDE:
            polys.append(shapely.geometry.shape(g))
            vals.append(int(v))
    polys, vals = np.asarray(polys, dtype=object), np.asarray(vals)
    tree = shapely.STRtree(polys)
    li, pj = tree.query(motor_utm, predicate="intersects")
    pieces = shapely.intersection(motor_utm[li], polys[pj])
    pstate = vals[pj]

    gx0 = math.floor(wx0 / GRID_M) * GRID_M
    gy0 = math.floor(wy0 / GRID_M) * GRID_M
    nx = int(math.ceil((wx1 - gx0) / GRID_M))
    ny = int(math.ceil((wy1 - gy0) / GRID_M))
    cells = {}
    for i in range(nx):
        for j in range(ny):
            cell = box(gx0 + i * GRID_M, gy0 + j * GRID_M, gx0 + (i + 1) * GRID_M,
                       gy0 + (j + 1) * GRID_M).intersection(window_box)
            if cell.is_empty:
                continue
            c = shapely.centroid(cell)
            d = np.hypot(px - shapely.get_x(c), py - shapely.get_y(c))
            n = int(np.argmin(d))
            cells[(i, j)] = {"cell_lower_left_utm": [gx0 + i * GRID_M, gy0 + j * GRID_M],
                             "cell_area_km2": round(cell.area / 1e6, 4),
                             "nearest_place": places[n][0], "place_type": places[n][1],
                             "distance_to_place_m": round(float(d[n]), 1),
                             "ms_flooded": 0, "ms_total": 0, "ms_not_assessed": 0,
                             "motor_km_flooded": 0.0, "motor_km_total": 0.0,
                             "osm_buildings": 0, "_geom": cell}

    def cell_of(x, y):
        return (int(math.floor((x - gx0) / GRID_M)), int(math.floor((y - gy0) / GRID_M)))

    def inside_window(x, y):
        return wx0 <= x < wx1 and wy0 < y <= wy1

    for pt, s in zip(cent, b_state):
        if s == OUTSIDE:
            continue
        x, y = shapely.get_x(pt), shapely.get_y(pt)
        if not inside_window(x, y):
            continue
        c = cells.get(cell_of(x, y))
        if c is None:
            continue
        c["ms_total"] += 1
        c["ms_flooded"] += int(s == FLOODED)
        c["ms_not_assessed"] += int(s == NOT_ASSESSED)
    for pt, s in zip(o_cent, o_state):
        if s == OUTSIDE:
            continue
        x, y = shapely.get_x(pt), shapely.get_y(pt)
        if inside_window(x, y) and cell_of(x, y) in cells:
            cells[cell_of(x, y)]["osm_buildings"] += 1
    cgeoms = np.asarray([c["_geom"] for c in cells.values()], dtype=object)
    ckeys = list(cells)
    ctree = shapely.STRtree(cgeoms)
    pi, ci = ctree.query(pieces, predicate="intersects")
    lens = shapely.length(shapely.intersection(pieces[pi], cgeoms[ci]))
    for k, c_idx, L in zip(pi, ci, lens):
        c = cells[ckeys[c_idx]]
        c["motor_km_total"] += L / 1000.0
        if pstate[k] == FLOODED:
            c["motor_km_flooded"] += L / 1000.0
    rows = []
    for c in cells.values():
        c.pop("_geom")
        c["motor_km_flooded"] = round(c["motor_km_flooded"], 4)
        c["motor_km_total"] = round(c["motor_km_total"], 4)
        c["motor_flooded_share"] = round(c["motor_km_flooded"] / c["motor_km_total"], 4) \
            if c["motor_km_total"] > 0 else None
        rows.append(c)
    top_b = sorted(rows, key=lambda r: (-r["ms_flooded"], -r["ms_total"]))[:10]
    top_r = sorted(rows, key=lambda r: (-r["motor_km_flooded"], -r["motor_km_total"]))[:10]
    top_s = sorted([r for r in rows if r["motor_km_total"] >= SHARE_MIN_ROAD_KM],
                   key=lambda r: (-r["motor_flooded_share"], -r["motor_km_flooded"]))[:10]
    osm_top = sorted(rows, key=lambda r: -r["osm_buildings"])[:10]
    return {"cells": len(rows), "cell_m": GRID_M, "share_min_road_km": SHARE_MIN_ROAD_KM,
            "ms_flooded_in_grid": sum(r["ms_flooded"] for r in rows),
            "ms_total_in_grid": sum(r["ms_total"] for r in rows),
            "top10_flooded_buildings": top_b, "top10_flooded_motor_km": top_r,
            "top10_flooded_motor_share": top_s, "osm_top_cells": osm_top}


def cross_check(osm_top, o_els, o_utm, o_cent, o_state, ms_utm):
    """M6: OSM building polygons in the ten OSM-densest cells vs Microsoft footprints."""
    tree = shapely.STRtree(ms_utm)
    out_cells, tot = [], defaultdict(int)
    for c in osm_top:
        x0, y0 = c["cell_lower_left_utm"]
        cell = box(x0, y0, x0 + GRID_M, y0 + GRID_M)
        idx = [k for k, (pt, s) in enumerate(zip(o_cent, o_state))
               if s != OUTSIDE and cell.contains(pt)]
        polys = [k for k in idx if o_utm[k].geom_type in ("Polygon", "MultiPolygon")]
        nodes = len(idx) - len(polys)
        matched = touched = 0
        for k in polys:
            g = o_utm[k]
            cand = tree.query(g, predicate="intersects")
            if len(cand) == 0:
                continue
            touched += 1
            inter = shapely.area(shapely.intersection(g, ms_utm[cand]))
            if np.any((inter >= MATCH_SHARE * g.area) | (inter >= MATCH_SHARE * shapely.area(ms_utm[cand]))):
                matched += 1
        out_cells.append({"cell_lower_left_utm": c["cell_lower_left_utm"], "nearest_place": c["nearest_place"],
                          "osm_buildings": len(idx), "osm_polygons": len(polys), "osm_nodes_skipped": nodes,
                          "matched": matched, "touched_by_any_footprint": touched,
                          "ms_footprints_in_cell": c["ms_total"]})
        for k, v in (("osm_buildings", len(idx)), ("osm_polygons", len(polys)),
                     ("osm_nodes_skipped", nodes), ("matched", matched), ("touched", touched),
                     ("ms_footprints", c["ms_total"])):
            tot[k] += v
    tot = dict(tot)
    tot["matched_share"] = round(tot["matched"] / tot["osm_polygons"], 4) if tot.get("osm_polygons") else None
    tot["touched_share"] = round(tot["touched"] / tot["osm_polygons"], 4) if tot.get("osm_polygons") else None
    return {"match_rule": f"intersection >= {MATCH_SHARE:.0%} of the OSM polygon or of the footprint",
            "cells": out_cells, "totals": tot}


if __name__ == "__main__":
    main()
