"""Step 6a: select and check the Jhelum centreline (rules J1-J4, J6-J9).

From the existing Geofabrik extract (step 4), takes the ways with role
main_stream in OSM relation 8306691 (the Jhelum River, wikidata Q214874) at
full length, and checks that inside the clip box they form one connected
chain. Lists named Jhelum ways outside the relation, members missing from the
extract, and every river way that touches the chain without being selected.
Checks that the centreline follows the 2014 persistent-water class.

Rule J5 (no missing member could be closer than the mapped centreline) needs
the distance grid and is checked in 24_river_distance.py before any result.

The mapped centreline through Wular Lake is kept; the lake is not treated as
the river (rule J8).

Outputs: outputs/jhelum_centreline.geojson, outputs/jhelum.json,
outputs/figures/jhelum_centreline.png

Usage:  .venv/bin/python scripts/23_jhelum.py
"""

import json
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import osmium
import rasterio
import rasterio.features
import shapely
from rasterio.warp import transform_geom
from scipy import ndimage
from shapely.geometry import LineString, box, mapping, shape

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "outputs"
PBF = ROOT / "data" / "osm" / "geofabrik" / "northern-zone-latest.osm.pbf"
RELATION_ID = 8306691                     # J2
EXPECTED = {"type": "waterway", "waterway": "river", "name:en": "Jhelum River", "wikidata": "Q214874"}
CLIP_BOX = (74.2, 33.5, 75.4, 34.6)
CRS = "EPSG:32643"
NEAR_M = (30, 60, 90)                     # J7
MIN_SHARE_90M = 50.0                      # J7 stop rule
STEP_M = 10.0                             # centreline sampling step for J7


def stop(msg):
    sys.exit(f"STOP: {msg}")


def to_utm(coords):
    return shape(transform_geom("EPSG:4326", CRS, {"type": "LineString", "coordinates": coords}))


def jname(tags):
    return tags.get("name:en") or tags.get("name")


def components(way_ids, refs):
    """Connected groups of ways that share at least one OSM node."""
    by_node = defaultdict(set)
    for w in way_ids:
        for n in refs[w]:
            by_node[n].add(w)
    seen, comps = set(), []
    for w in way_ids:
        if w in seen:
            continue
        stack, comp = [w], set()
        while stack:
            x = stack.pop()
            if x in comp:
                continue
            comp.add(x)
            for n in refs[x]:
                stack.extend(by_node[n] - comp)
        seen |= comp
        comps.append(sorted(comp))
    return comps


def main():
    (OUT / "figures").mkdir(parents=True, exist_ok=True)
    extract = json.loads((ROOT / "data" / "osm" / "geofabrik" / "extract.json").read_text())
    if not PBF.exists():
        stop(f"{PBF} is missing; run 12_geofabrik_download.py")

    # The relation (J2).
    rel = None
    for r in osmium.FileProcessor(str(PBF), osmium.osm.RELATION).with_filter(osmium.filter.IdFilter([RELATION_ID])):
        rel = {"tags": dict(r.tags), "members": [(m.type, m.ref, m.role) for m in r.members]}
    if rel is None:
        stop(f"relation {RELATION_ID} not in the extract")
    wrong = {k: rel["tags"].get(k) for k, v in EXPECTED.items() if rel["tags"].get(k) != v}
    if wrong:
        stop(f"relation {RELATION_ID} tags differ from the rules: {wrong}")
    main_ids = [ref for t, ref, role in rel["members"] if t == "w" and role == "main_stream"]
    side_ids = [ref for t, ref, role in rel["members"] if t == "w" and role == "side_stream"]

    # All river ways (for the main stem, the tributary check and named-Jhelum ways).
    rivers = {}
    for o in osmium.FileProcessor(str(PBF), osmium.osm.NODE | osmium.osm.WAY).with_locations() \
            .with_filter(osmium.filter.TagFilter(("waterway", "river"))) \
            .with_filter(osmium.filter.EntityFilter(osmium.osm.WAY)):
        rivers[o.id] = {"tags": dict(o.tags), "refs": [n.ref for n in o.nodes],
                        "coords": [(n.lon, n.lat) for n in o.nodes if n.location.valid()],
                        "missing_locations": sum(1 for n in o.nodes if not n.location.valid())}
    # Main-stream members are selected by relation membership (J2), whatever
    # their own waterway tag, so load them by id.
    members = {}
    for o in osmium.FileProcessor(str(PBF), osmium.osm.NODE | osmium.osm.WAY).with_locations() \
            .with_filter(osmium.filter.IdFilter(main_ids)) \
            .with_filter(osmium.filter.EntityFilter(osmium.osm.WAY)):
        members[o.id] = {"tags": dict(o.tags), "refs": [n.ref for n in o.nodes],
                         "coords": [(n.lon, n.lat) for n in o.nodes if n.location.valid()],
                         "missing_locations": sum(1 for n in o.nodes if not n.location.valid())}
    present = [w for w in main_ids if w in members]
    missing = [w for w in main_ids if w not in members]
    for w in present:
        rivers.setdefault(w, members[w])
    if not present:
        stop("no main_stream member of the relation is in the extract")
    if any(rivers[w]["missing_locations"] for w in present):
        stop("a selected way has nodes without location in the extract")

    refs = {w: rivers[w]["refs"] for w in present}
    utm = {w: to_utm(rivers[w]["coords"]) for w in present}
    clip_utm = shape(transform_geom("EPSG:4326", CRS, mapping(box(*CLIP_BOX))))
    in_box = [w for w in present if utm[w].intersects(clip_utm)]

    # Continuity (J4): components before joining = ways; after = shared-node groups.
    comps_box = components(in_box, refs)
    comps_all = components(present, refs)
    if len(comps_box) != 1:
        stop(f"the selected ways inside the clip box form {len(comps_box)} components, not one")
    union = shapely.line_merge(shapely.union_all([utm[w] for w in present]))
    in_box_geom = union.intersection(clip_utm)

    # Analysis area (pixel outlines) for lengths and the persistent-water check.
    with rasterio.open(ROOT / "data" / "flood" / "flood_classes.tif") as src:
        cls = src.read(1)
        T, shp = src.transform, src.shape
        if src.crs.to_string() != CRS:
            stop(f"flood raster CRS {src.crs} is not {CRS}")
    area = rasterio.open(ROOT / "data" / "flood" / "analysis_area.tif").read(1).astype(bool)
    area_poly = shapely.union_all([shape(g) for g, v in rasterio.features.shapes(
        area.astype(np.uint8), mask=area, transform=T) if v == 1])
    in_area_geom = union.intersection(area_poly)

    # J7: distance from points along the centreline (inside the analysis area)
    # to the nearest persistent-water pixel, and to water on either date.
    def dist_to(mask):
        return ndimage.distance_transform_edt(~mask) * T.a

    d_pers = dist_to(area & (cls == 2))
    d_water = dist_to(area & np.isin(cls, (1, 2, 3)))
    parts = getattr(in_area_geom, "geoms", [in_area_geom])
    pts = []
    for g in parts:
        if g.length == 0:
            continue
        n = max(2, int(g.length // STEP_M) + 1)
        pts.extend(shapely.line_interpolate_point(g, np.linspace(0, g.length, n)))
    xs, ys = shapely.get_x(np.asarray(pts)), shapely.get_y(np.asarray(pts))
    cols, rows = ~T * (xs, ys)
    rows, cols = np.floor(rows).astype(int), np.floor(cols).astype(int)
    ok = (rows >= 0) & (rows < shp[0]) & (cols >= 0) & (cols < shp[1])
    dp, dw = d_pers[rows[ok], cols[ok]], d_water[rows[ok], cols[ok]]
    check = {f"within_{m}m_of_persistent_water_pct": round(100 * float(np.mean(dp <= m)), 2) for m in NEAR_M}
    check["within_90m_of_water_either_date_pct"] = round(100 * float(np.mean(dw <= 90)), 2)
    check["sample_points"] = int(ok.sum())
    check["sample_step_m"] = STEP_M
    check["median_distance_to_persistent_water_m"] = round(float(np.median(dp)), 1)
    if check["within_90m_of_persistent_water_pct"] < MIN_SHARE_90M:
        stop(f"only {check['within_90m_of_persistent_water_pct']} % of the centreline is within 90 m of "
             "persistent water; the geometry looks misplaced")

    # J6: river ways touching the chain but not selected.
    chain_nodes = {n for w in present for n in refs[w]}
    touching = sorted((w for w, r in rivers.items()
                       if w not in present and chain_nodes & set(r["refs"])),
                      key=lambda w: w)
    touching_rec = [{"way_id": w, "name": jname(rivers[w]["tags"]),
                     "length_km": round(to_utm(rivers[w]["coords"]).length / 1000, 2),
                     "relation_role": "side_stream" if w in side_ids else None}
                    for w in touching]
    named_outside = [{"way_id": w, "name": jname(r["tags"]),
                      "length_km": round(to_utm(r["coords"]).length / 1000, 2),
                      "shares_nodes_with_chain": bool(chain_nodes & set(r["refs"])),
                      "min_distance_to_chain_m": round(to_utm(r["coords"]).distance(union), 1)}
                     for w, r in rivers.items()
                     if w not in present and any(s in (jname(r["tags"]) or "").lower()
                                                 for s in ("jhelum", "jehlum", "veth"))
                     and to_utm(r["coords"]).intersects(clip_utm)]

    # Outputs.
    features = []
    for w in present:
        features.append({"type": "Feature",
                         "properties": {"osm_way_id": w, "name": jname(rivers[w]["tags"]),
                                        "relation": RELATION_ID, "role": "main_stream",
                                        "length_m": round(utm[w].length, 1),
                                        "in_clip_box": w in in_box},
                         "geometry": {"type": "LineString", "coordinates": rivers[w]["coords"]}})
    (OUT / "jhelum_centreline.geojson").write_text(json.dumps(
        {"type": "FeatureCollection",
         "properties": {"source": extract["file"], "osm_timestamp": extract["osm_timestamp_header"],
                        "relation": RELATION_ID, "selection": "main_stream members present in the extract"},
         "features": features}, separators=(",", ":")))
    record = {
        "script": "scripts/23_jhelum.py", "rules": "scripts/step6_rules.md (J)",
        "source": {"file": extract["file"], "sha256": extract["sha256"],
                   "osm_timestamp": extract["osm_timestamp_header"]},
        "relation": {"id": RELATION_ID, "tags": {k: rel["tags"].get(k) for k in
                                                 ("type", "waterway", "name", "name:en", "wikidata")},
                     "main_stream_members": len(main_ids), "side_stream_members": len(side_ids)},
        "selected_ways": [{"way_id": w, "name": jname(rivers[w]["tags"]),
                           "waterway_tag": rivers[w]["tags"].get("waterway"), "in_clip_box": w in in_box,
                           "length_km": round(utm[w].length / 1000, 3)} for w in present],
        "main_stream_members_missing_from_extract": missing,
        "lengths_km": {"selected_total": round(sum(utm[w].length for w in present) / 1000, 3),
                       "inside_clip_box": round(in_box_geom.length / 1000, 3),
                       "inside_analysis_area": round(in_area_geom.length / 1000, 3)},
        "continuity": {"ways_inside_clip_box": len(in_box),
                       "components_before_joining_inside_box": len(in_box),
                       "components_after_joining_inside_box": len(comps_box),
                       "components_after_joining_all_selected": len(comps_all),
                       "merged_geometry_type": union.geom_type,
                       "merged_parts": len(getattr(union, "geoms", [union]))},
        "jhelum_named_ways_not_selected": named_outside,
        "river_ways_touching_chain_not_selected": touching_rec,
        "persistent_water_check": check,
        "wular": "The mapped centreline through Wular Lake is kept; the lake is not treated as the "
                 "river and distances are not measured to its shore.",
    }
    (OUT / "jhelum.json").write_text(json.dumps(record, indent=1))

    # Figure: centreline over the analysis area and flood classes.
    fig, ax = plt.subplots(figsize=(9, 9.5))
    rgb = np.full(shp + (3,), 0.96)
    rgb[area] = (0.86, 0.84, 0.80)
    rgb[area & (cls == 1)] = (0.90, 0.62, 0.0)
    rgb[area & (cls == 2)] = (0.0, 0.45, 0.70)
    ext = (T.c, T.c + shp[1] * T.a, T.f + shp[0] * T.e, T.f)
    ax.imshow(rgb, extent=ext, interpolation="nearest")
    for g in getattr(union, "geoms", [union]):
        x, y = g.xy
        ax.plot(x, y, color="black", lw=1.1)
    ax.set_xlim(ext[0], ext[1])
    ax.set_ylim(ext[2], ext[3])
    ax.set_title("Selected Jhelum centreline (OSM relation 8306691, main stream) over the analysis area\n"
                 "grey: analysis area; orange: flooded 10 Sep 2014; blue: water on both dates")
    ax.set_xlabel("UTM 43N easting (m)")
    ax.set_ylabel("UTM 43N northing (m)")
    fig.savefig(OUT / "figures" / "jhelum_centreline.png", dpi=130, bbox_inches="tight")
    print(json.dumps({k: record[k] for k in ("lengths_km", "continuity", "persistent_water_check",
                                            "main_stream_members_missing_from_extract",
                                            "jhelum_named_ways_not_selected",
                                            "river_ways_touching_chain_not_selected")}, indent=1))


if __name__ == "__main__":
    main()
