"""Step 4b: count flooded buildings and road length under scripts/step4_rules.md.

Reads the OSM elements from 12b_read_extract.py (Geofabrik extract, amendment
1) and the step 3 rasters,
and implements rules B2-B4, R1-R3, S1, W1, A2-A4 and C1-C2 as written in
step4_rules.md (committed before this script was run).

Pixel state, used for every count: flooded (flood class 1), dry (classes 0,
2, 3), not assessed (analysis-area pixel with class 255), outside (not in the
analysis area). Roads are cut along the exact pixel outlines of a combined
state raster (state x water/wetland mask x Srinagar window), so every metre of
road is assigned to exactly one combination.

Outputs: data/exposure_report.json, data/exposure/buildings.csv,
data/exposure/roads.csv, data/exposure/flooded_road_pieces.json (for previews).

Usage:  .venv/bin/python scripts/13_exposure.py
"""

import csv
import importlib.util
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import rasterio
import shapely
from rasterio.features import rasterize, shapes
from rasterio.warp import transform as warp_transform
from rasterio.warp import transform_bounds
from rasterio.windows import from_bounds
from shapely.geometry import Point, box, shape

ROOT = Path(__file__).resolve().parent.parent
EXPO_DIR = ROOT / "data" / "exposure"
OSM_DIR = ROOT / "data" / "osm" / "exposure"
SRINAGAR_BBOX = (74.65, 33.95, 74.95, 34.20)
MOTOR = {"motorway", "motorway_link", "trunk", "trunk_link", "primary", "primary_link",
         "secondary", "secondary_link", "tertiary", "tertiary_link", "residential",
         "unclassified", "service", "living_street"}                       # R2
OUTSIDE, DRY, FLOODED, NOT_ASSESSED = 0, 1, 2, 3
STATE_NAMES = {OUTSIDE: "outside", DRY: "dry", FLOODED: "flooded", NOT_ASSESSED: "not_assessed"}
COVERAGE_MIN, MIN_POLYGONS = 0.90, 5                                       # A2, A3
GRID_KM = 1000.0                                                            # A3
CELL_C = 2000.0                                                             # C1
ROAD_MIN_KM = 2.0                                                           # C1

_spec = importlib.util.spec_from_file_location("osm3", Path(__file__).parent / "09_osm.py")
osm3 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(osm3)


def km(x):
    return round(float(x) / 1000.0, 4)


def load_elements():
    manifest = json.loads((OSM_DIR / "manifest.json").read_text())
    seen = {}
    for t in manifest["tiles"]:
        for el in json.loads((ROOT / t["file"]).read_text())["elements"]:
            seen[(el["type"], el["id"])] = el
    return manifest, list(seen.values())


def to_utm(geoms, crs):
    """Project WGS84 shapely geometries to the mosaic CRS in one call."""
    def f(coords):
        xs, ys = warp_transform("EPSG:4326", crs, coords[:, 0].tolist(), coords[:, 1].tolist())
        return np.column_stack([xs, ys])
    return shapely.transform(np.asarray(geoms, dtype=object), f)


def pixel_index(transform, xs, ys, shape_):
    inv = ~transform
    cols, rows = inv * (np.asarray(xs), np.asarray(ys))
    rows, cols = np.floor(rows).astype(int), np.floor(cols).astype(int)
    inside = (rows >= 0) & (rows < shape_[0]) & (cols >= 0) & (cols < shape_[1])
    return rows, cols, inside


def main():
    EXPO_DIR.mkdir(parents=True, exist_ok=True)
    manifest, elements = load_elements()
    with rasterio.open(ROOT / "data" / "flood" / "flood_classes.tif") as src:
        cls = src.read(1)
        transform, crs, shp = src.transform, src.crs, (src.height, src.width)
    area = rasterio.open(ROOT / "data" / "flood" / "analysis_area.tif").read(1).astype(bool)
    osm_wet = rasterio.open(ROOT / "data" / "osm" / "water_wetland_mask.tif").read(1).astype(bool)

    state = np.full(shp, OUTSIDE, dtype=np.uint8)
    state[area & np.isin(cls, (0, 2, 3))] = DRY
    state[area & (cls == 1)] = FLOODED
    state[area & (cls == 255)] = NOT_ASSESSED
    sb = transform_bounds("EPSG:4326", crs, *SRINAGAR_BBOX, densify_pts=21)
    win = from_bounds(*sb, transform=transform).round_offsets(op="floor").round_lengths(op="ceil")
    srin = np.zeros(shp, bool)
    srin[win.toslices()] = True
    code = state + 4 * osm_wet.astype(np.uint8) + 8 * srin.astype(np.uint8)

    # ---- Build geometries (B2, R1)
    b_els, b_geoms, unbuilt_b = [], [], Counter()
    r_els, r_geoms, area_ways = [], [], Counter()
    for el in elements:
        tags = el.get("tags", {})
        if "building" in tags:
            if el["type"] == "node":
                g = Point(el["lon"], el["lat"])
            else:
                g = osm3.element_polygon(el)
            if g is None or g.is_empty:
                unbuilt_b[el["type"]] += 1
            else:
                b_els.append(el)
                b_geoms.append(g)
        if el["type"] == "way" and "highway" in tags:
            if tags.get("area") == "yes":
                area_ways[tags["highway"]] += 1
                continue
            c = osm3.way_coords(el)
            if len(c) >= 2:
                r_els.append(el)
                r_geoms.append(shapely.LineString(c))
    print(f"{len(b_geoms)} buildings, {len(r_geoms)} road ways")
    b_utm = to_utm(b_geoms, crs)
    r_utm = to_utm(r_geoms, crs)

    # ---- Buildings: centroid rule (B3)
    cent = shapely.centroid(b_utm)
    cx, cy = shapely.get_x(cent), shapely.get_y(cent)
    rr, cc, ok = pixel_index(transform, cx, cy, shp)
    b_code = np.full(len(b_utm), 0, dtype=np.uint8)
    b_code[ok] = code[rr[ok], cc[ok]]
    b_state = b_code & 3

    # ---- Buildings: overlap rule (B4)
    ov_state = np.full(len(b_utm), OUTSIDE, dtype=np.uint8)
    for i, g in enumerate(b_utm):
        if g.geom_type == "Point":
            ov_state[i] = b_state[i]
            continue
        minx, miny, maxx, maxy = g.bounds
        c0, r1_ = ~transform * (minx, miny)
        c1, r0_ = ~transform * (maxx, maxy)
        r0, r1 = max(int(math.floor(r0_)), 0), min(int(math.floor(r1_)) + 1, shp[0])
        c0, c1 = max(int(math.floor(c0)), 0), min(int(math.floor(c1)) + 1, shp[1])
        if r1 <= r0 or c1 <= c0:
            continue
        wt = rasterio.windows.transform(rasterio.windows.Window(c0, r0, c1 - c0, r1 - r0), transform)
        touched = rasterize([(g, 1)], out_shape=(r1 - r0, c1 - c0), transform=wt,
                            fill=0, all_touched=True, dtype="uint8").astype(bool)
        st = state[r0:r1, c0:c1][touched]
        if (st == FLOODED).any():
            ov_state[i] = FLOODED
        elif (st == NOT_ASSESSED).any():
            ov_state[i] = NOT_ASSESSED
        elif (st == DRY).any():
            ov_state[i] = DRY

    def bcount(mask_codes, states):
        out = {STATE_NAMES[s]: int(np.sum(mask_codes & (states == s))) for s in (DRY, FLOODED, NOT_ASSESSED)}
        out["total_in_analysis_area"] = sum(out.values())
        return out

    in_osm = (b_code & 4) > 0
    in_srin = (b_code & 8) > 0
    every = np.ones(len(b_utm), bool)
    buildings = {
        "centroid": {"analysis_area": bcount(every, b_state), "srinagar_box": bcount(in_srin, b_state),
                     "inside_water_wetland": bcount(in_osm, b_state),
                     "outside_water_wetland": bcount(~in_osm, b_state),
                     "srinagar_inside_water_wetland": bcount(in_srin & in_osm, b_state),
                     "srinagar_outside_water_wetland": bcount(in_srin & ~in_osm, b_state)},
        "overlap": {"analysis_area": bcount(every, ov_state), "srinagar_box": bcount(in_srin, ov_state)},
        "elements_by_type": dict(Counter(e["type"] for e in b_els)),
        "not_turned_into_geometry": dict(unbuilt_b),
        "building_values_top": Counter(e["tags"]["building"] for e, s in zip(b_els, b_state)
                                       if s != OUTSIDE).most_common(15),
        "building_no_in_area": int(sum(1 for e, s in zip(b_els, b_state)
                                       if s != OUTSIDE and e["tags"]["building"] == "no")),
    }

    # ---- Roads (R3, S1, W1): cut along exact pixel outlines of the code raster.
    near = shapely.box(*shapely.total_bounds(r_utm))
    polys, pcodes = [], []
    for geom, val in shapes(code, mask=None, transform=transform, connectivity=4):
        if int(val) & 3 == OUTSIDE:
            continue          # outside the analysis area: not needed for lengths
        polys.append(shape(geom))
        pcodes.append(int(val))
    polys = np.asarray(polys, dtype=object)
    pcodes = np.asarray(pcodes)
    tree = shapely.STRtree(polys)
    li, pj = tree.query(r_utm, predicate="intersects")
    pieces = shapely.intersection(r_utm[li], polys[pj])
    plen = shapely.length(pieces)
    hw = np.array([e["tags"]["highway"] for e in r_els])
    motor = np.isin(hw, list(MOTOR))
    piece_state = pcodes[pj] & 3
    piece_osm = (pcodes[pj] & 4) > 0
    piece_srin = (pcodes[pj] & 8) > 0

    def rlen(sel):
        out = {STATE_NAMES[s]: km(plen[sel & (piece_state == s)].sum())
               for s in (DRY, FLOODED, NOT_ASSESSED)}
        out["total_in_analysis_area"] = round(sum(out.values()), 4)
        return out

    pm = motor[li]
    allp = np.ones(len(li), bool)
    roads = {
        "motor": {"analysis_area": rlen(pm), "srinagar_box": rlen(pm & piece_srin),
                  "inside_water_wetland": rlen(pm & piece_osm),
                  "outside_water_wetland": rlen(pm & ~piece_osm),
                  "srinagar_inside_water_wetland": rlen(pm & piece_srin & piece_osm),
                  "srinagar_outside_water_wetland": rlen(pm & piece_srin & ~piece_osm)},
        "other": {"analysis_area": rlen(~pm), "srinagar_box": rlen(~pm & piece_srin),
                  "inside_water_wetland": rlen(~pm & piece_osm),
                  "outside_water_wetland": rlen(~pm & ~piece_osm)},
        "other_by_value": {v: rlen((hw[li] == v)) for v in sorted(set(hw[~motor]))},
        "motor_by_value": {v: rlen((hw[li] == v)) for v in sorted(set(hw[motor]))},
        "area_yes_ways_excluded": dict(area_ways),
        "ways": int(len(r_utm)),
    }
    roads["other_by_value"] = {k: v for k, v in roads["other_by_value"].items()
                               if v["total_in_analysis_area"] > 0}

    # ---- Areas within Srinagar (A2-A4)
    areas = srinagar_areas(crs, cent, b_state, in_srin, li, pieces, piece_state, pm, sb)

    # ---- Completeness (C1, C2)
    completeness = completeness_check(transform, shp, state, cent, b_state, r_utm, motor, in_srin, sb)

    report = {"rules_file": "scripts/step4_rules.md",
              "rules_sha256": __import__("hashlib").sha256(
                  (Path(__file__).parent / "step4_rules.md").read_bytes()).hexdigest(),
              "osm_source": {"extract": manifest["extract"],
                             "elements_file": manifest["tiles"][0],
                             "areas_file": manifest["areas"],
                             "read_stats": manifest["read_stats"],
                             "elements_unique": len(elements),
                             "water_wetland_snapshot_step3": json.loads(
                                 (ROOT / "data" / "osm_report.json").read_text())["osm_base_timestamp"]},
              "buildings": buildings, "roads": roads, "areas": areas,
              "completeness": completeness}
    (ROOT / "data" / "exposure_report.json").write_text(json.dumps(report, indent=2))

    with open(EXPO_DIR / "buildings.csv", "w", newline="") as f:
        wcsv = csv.writer(f)
        wcsv.writerow(["osm_type", "osm_id", "building", "centroid_x", "centroid_y",
                       "state_centroid", "state_overlap", "in_water_wetland", "in_srinagar_box"])
        for e, x, y, s, o, wo, sr in zip(b_els, cx, cy, b_state, ov_state, in_osm, in_srin):
            if s != OUTSIDE or o != OUTSIDE:
                wcsv.writerow([e["type"], e["id"], e["tags"]["building"], round(x, 2), round(y, 2),
                               STATE_NAMES[s], STATE_NAMES[o], int(wo), int(sr)])
    per_road = defaultdict(lambda: defaultdict(float))
    for i, s, L in zip(li, piece_state, plen):
        per_road[i][STATE_NAMES[s]] += L
    with open(EXPO_DIR / "roads.csv", "w", newline="") as f:
        wcsv = csv.writer(f)
        wcsv.writerow(["osm_id", "highway", "motor", "m_dry", "m_flooded", "m_not_assessed"])
        for i, d in sorted(per_road.items()):
            wcsv.writerow([r_els[i]["id"], hw[i], int(motor[i]), round(d["dry"], 2),
                           round(d["flooded"], 2), round(d["not_assessed"], 2)])
    fl = (piece_state == FLOODED) & piece_srin
    (EXPO_DIR / "flooded_road_pieces.json").write_text(json.dumps(
        {"crs": str(crs), "motor": [shapely.to_geojson(g) for g in pieces[fl & pm]],
         "other": [shapely.to_geojson(g) for g in pieces[fl & ~pm]]}))
    print(json.dumps({k: report[k] for k in ("buildings", "roads")}, indent=1)[:6000])
    print(json.dumps(report["areas"], indent=1)[:3000])
    print(json.dumps(report["completeness"], indent=1))


def srinagar_areas(crs, cent, b_state, in_srin, li, pieces, piece_state, pm, sb):
    data = json.loads((OSM_DIR / "srinagar_areas.json").read_text())
    box_utm = box(*sb)
    cand = defaultdict(list)       # set name -> [(name, polygon)]
    places = []                    # (name, kind, point) for grid labels
    for el in data["elements"]:
        tags = el.get("tags", {})
        name = tags.get("name:en") or tags.get("name")
        # Boundary relations are tagged type=boundary; they are assembled from
        # their outer/inner ways exactly like multipolygons.
        as_mp = {**el, "tags": {k: v for k, v in tags.items() if k != "type"}} \
            if el["type"] == "relation" and tags.get("type") == "boundary" else el
        poly = None if el["type"] == "node" else osm3.element_polygon(as_mp)
        if poly is not None and not poly.is_empty:
            poly = to_utm([poly], crs)[0]
        if tags.get("boundary") == "administrative" and poly is not None:
            cand[f"admin_level_{tags.get('admin_level', 'none')}"].append((name, poly))
        if tags.get("place") in ("suburb", "neighbourhood", "quarter"):
            if poly is not None:
                cand["place_polygons"].append((name, poly))
            pt = to_utm([Point(el["lon"], el["lat"])], crs)[0] if el["type"] == "node" \
                else (shapely.centroid(poly) if poly is not None else None)
            if name and pt is not None:
                places.append((name, tags["place"], pt))

    sel = in_srin & (b_state != OUTSIDE)
    pts = cent[sel]
    evaluation = {}
    for set_name, items in cand.items():
        in_box = [(n, p) for n, p in items if p.intersects(box_utm)]
        cover = 0.0
        if in_box and len(pts):
            u = shapely.union_all([p for _, p in in_box])
            cover = float(np.mean(shapely.contains(u, pts)))
        evaluation[set_name] = {"polygons_in_box": len(in_box), "coverage_of_box_buildings": round(cover, 4)}
    qualifying = [k for k, v in evaluation.items()
                  if v["polygons_in_box"] >= MIN_POLYGONS and v["coverage_of_box_buildings"] >= COVERAGE_MIN]

    def rank(k):
        return (1, 0) if k == "place_polygons" else (0, int(k.split("_")[-1]) if k.split("_")[-1].isdigit() else -1)

    states_sel = b_state[sel]
    fl_pieces = pieces[(piece_state == FLOODED) & pm]
    result = {"evaluation": evaluation, "places_with_names": len(places)}
    if qualifying:
        chosen = max(qualifying, key=rank)
        items = [(n, p) for n, p in cand[chosen] if p.intersects(box_utm)]
        geoms = np.asarray([p for _, p in items], dtype=object)
        tree = shapely.STRtree(geoms)
        pi, gi = tree.query(pts, predicate="within")
        assign = {}
        for a, b in zip(pi, gi):     # smallest polygon wins if they overlap
            if a not in assign or geoms[b].area < geoms[assign[a]].area:
                assign[a] = b
        rows = []
        for gidx, (n, p) in enumerate(items):
            members = [a for a, b in assign.items() if b == gidx]
            st = states_sel[members] if members else np.array([], int)
            rows.append({"name": n, "flooded_buildings": int(np.sum(st == FLOODED)),
                         "buildings": int(len(members)),
                         "not_assessed_buildings": int(np.sum(st == NOT_ASSESSED)),
                         "flooded_motor_road_km": km(shapely.length(shapely.intersection(fl_pieces, p)).sum())
                         if len(fl_pieces) else 0.0})
        result.update({"method": "polygons", "chosen_set": chosen,
                       "buildings_unassigned": int(len(pts) - len(assign))})
    else:
        rows = grid_rows(places, pts, states_sel, fl_pieces, sb)
        result.update({"method": "1 km grid", "chosen_set": None})
    rows.sort(key=lambda r: (-r["flooded_buildings"], -r["buildings"]))
    result["top10"] = rows[:10]
    result["units_with_buildings"] = sum(1 for r in rows if r["buildings"])
    result["flooded_buildings_in_top10"] = sum(r["flooded_buildings"] for r in rows[:10])
    if result["method"] == "polygons" and places:
        # Supplementary, added after the rule's result was seen (owner's
        # decision): the A3 fallback grid, reported alongside, not instead.
        g = grid_rows(places, pts, states_sel, fl_pieces, sb)
        g.sort(key=lambda r: (-r["flooded_buildings"], -r["buildings"]))
        result["supplementary_grid"] = {
            "note": "added after the committed rule selected districts; not the rule result",
            "top10": g[:10], "cells_with_buildings": len(g),
            "flooded_buildings_in_top10": sum(r["flooded_buildings"] for r in g[:10])}
    return result


def grid_rows(places, pts, states_sel, fl_pieces, sb):
    """1 km UTM cells over the Srinagar box, labelled with the nearest named place (A3)."""
    if not places:
        sys.exit("STOP: no named suburb/neighbourhood/quarter places for the 1 km grid")
    minx, miny, maxx, maxy = sb
    x0, y0 = math.floor(minx / GRID_KM) * GRID_KM, math.floor(miny / GRID_KM) * GRID_KM
    px = np.array([shapely.get_x(p) for *_, p in places])
    py = np.array([shapely.get_y(p) for *_, p in places])
    bx, by = shapely.get_x(pts), shapely.get_y(pts)
    ci = np.floor((bx - x0) / GRID_KM).astype(int)
    cj = np.floor((by - y0) / GRID_KM).astype(int)
    cells = defaultdict(list)
    for k, key in enumerate(zip(ci, cj)):
        cells[key].append(k)
    rows = []
    for (i, j), members in cells.items():
        cxm, cym = x0 + (i + 0.5) * GRID_KM, y0 + (j + 0.5) * GRID_KM
        dist = np.hypot(px - cxm, py - cym)
        n = int(np.argmin(dist))
        st = states_sel[members]
        cell = box(x0 + i * GRID_KM, y0 + j * GRID_KM, x0 + (i + 1) * GRID_KM, y0 + (j + 1) * GRID_KM)
        rows.append({"cell_lower_left_utm": [x0 + i * GRID_KM, y0 + j * GRID_KM],
                     "nearest_place": places[n][0], "place_type": places[n][1],
                     "distance_to_place_m": round(float(dist[n]), 1),
                     "flooded_buildings": int(np.sum(st == FLOODED)), "buildings": len(members),
                     "not_assessed_buildings": int(np.sum(st == NOT_ASSESSED)),
                     "flooded_motor_road_km": km(shapely.length(shapely.intersection(fl_pieces, cell)).sum())
                     if len(fl_pieces) else 0.0})
    return rows


def completeness_check(transform, shp, state, cent, b_state, r_utm, motor, in_srin, sb):
    in_area = state != OUTSIDE
    rows, cols = np.nonzero(in_area)
    xs, ys = transform * (cols + 0.5, rows + 0.5)
    x0, y0 = math.floor(xs.min() / CELL_C) * CELL_C, math.floor(ys.min() / CELL_C) * CELL_C
    ci = np.floor((xs - x0) / CELL_C).astype(int)
    cj = np.floor((ys - y0) / CELL_C).astype(int)
    nx, ny = ci.max() + 1, cj.max() + 1
    area_px = np.zeros((nx, ny), int)
    np.add.at(area_px, (ci, cj), 1)
    flooded_px = np.zeros((nx, ny), int)
    fl = state[rows, cols] == FLOODED
    np.add.at(flooded_px, (ci[fl], cj[fl]), 1)
    # buildings per cell (centroid in analysis area)
    sel = b_state != OUTSIDE
    bi = np.floor((shapely.get_x(cent[sel]) - x0) / CELL_C).astype(int)
    bj = np.floor((shapely.get_y(cent[sel]) - y0) / CELL_C).astype(int)
    okb = (bi >= 0) & (bi < nx) & (bj >= 0) & (bj < ny)
    nb = np.zeros((nx, ny), int)
    np.add.at(nb, (bi[okb], bj[okb]), 1)
    # motor road length per cell, clipped to cells that contain analysis area
    cell_ids = np.argwhere(area_px > 0)
    cell_geoms = np.asarray([box(x0 + i * CELL_C, y0 + j * CELL_C, x0 + (i + 1) * CELL_C,
                                 y0 + (j + 1) * CELL_C) for i, j in cell_ids], dtype=object)
    tree = shapely.STRtree(cell_geoms)
    mlines = r_utm[motor]
    a, b = tree.query(mlines, predicate="intersects")
    lens = shapely.length(shapely.intersection(mlines[a], cell_geoms[b]))
    road = np.zeros((nx, ny))
    np.add.at(road, (cell_ids[b, 0], cell_ids[b, 1]), lens)
    road_km = road / 1000.0
    roads_ok = (area_px > 0) & (road_km >= ROAD_MIN_KM)
    flagged = roads_ok & (nb == 0)
    density = nb[roads_ok] / road_km[roads_ok]
    # Srinagar box cells vs rest (cell centre inside the box)
    cx = x0 + (np.arange(nx)[:, None] + 0.5) * CELL_C + 0 * np.arange(ny)[None, :]
    cy = y0 + (np.arange(ny)[None, :] + 0.5) * CELL_C + 0 * np.arange(nx)[:, None]
    in_box = (cx >= sb[0]) & (cx <= sb[2]) & (cy >= sb[1]) & (cy <= sb[3])
    def pcts(v):
        return [round(float(x), 2) for x in np.percentile(v, [10, 25, 50, 75, 90])] if v.size else None
    return {
        "cell_m": CELL_C, "road_min_km": ROAD_MIN_KM,
        "cells_with_analysis_area": int((area_px > 0).sum()),
        "cells_with_ge_2km_motor_road": int(roads_ok.sum()),
        "c1_cells_roads_no_buildings": int(flagged.sum()),
        "c1_analysis_area_km2_in_flagged": round(float(area_px[flagged].sum()) * 0.0009, 4),
        "c1_flooded_km2_in_flagged": round(float(flooded_px[flagged].sum()) * 0.0009, 4),
        "c1_motor_road_km_in_flagged": round(float(road_km[flagged].sum()), 3),
        "c2_buildings_per_motor_km_p10_p25_p50_p75_p90": pcts(density),
        "c2_srinagar_box_cells": pcts(nb[roads_ok & in_box] / road_km[roads_ok & in_box]),
        "c2_rest_cells": pcts(nb[roads_ok & ~in_box] / road_km[roads_ok & ~in_box]),
        "c2_cells_lt_1_building_per_km": int(np.sum(density < 1)),
        "grid_origin_utm": [x0, y0],
    }


if __name__ == "__main__":
    main()
