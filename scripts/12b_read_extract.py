"""Step 4a (cont.): read buildings, roads and Srinagar areas from the extract.

Implements amendment rule G3: from the verified Geofabrik file
(12_geofabrik_download.py), keeps
  - nodes, ways and relations with a `building` key, and ways with a
    `highway` key, if any of their nodes lies in the clip box;
  - boundary=administrative relations and place=suburb/neighbourhood/quarter
    nodes, ways and relations, if any of their nodes lies in the Srinagar box.
Ways get their node coordinates; relations get their member ways'
coordinates, so the geometry is built exactly as in step 3 / rule B2.

The elements are written in the same JSON layout the Overpass API uses, so
13_exposure.py counts them without change:
data/osm/exposure/extract_elements.json, data/osm/exposure/srinagar_areas.json
and data/osm/exposure/manifest.json (with the extract record).

Usage:  .venv/bin/python scripts/12b_read_extract.py
"""

import hashlib
import json
from collections import Counter
from pathlib import Path

import osmium

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "osm" / "exposure"
CLIP_BOX = (74.2, 33.5, 75.4, 34.6)
SRINAGAR_BBOX = (74.65, 33.95, 74.95, 34.20)
PLACES = {"suburb", "neighbourhood", "quarter"}


def inside(lon, lat, b):
    return b[0] <= lon <= b[2] and b[1] <= lat <= b[3]


def is_area_tags(tags):
    return tags.get("boundary") == "administrative" or tags.get("place") in PLACES


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rec = json.loads((ROOT / "data" / "osm" / "geofabrik" / "extract.json").read_text())
    pbf = str(ROOT / rec["file"])
    stats = Counter()

    # Pass 1: relations we may need, and the ways they reference.
    rels, member_ways = {}, set()
    for o in osmium.FileProcessor(pbf, osmium.osm.RELATION):
        tags = dict(o.tags)
        if "building" in tags or is_area_tags(tags):
            rels[o.id] = {"type": "relation", "id": o.id, "tags": tags,
                          "members": [{"type": {"w": "way", "n": "node", "r": "relation"}[m.type],
                                       "ref": m.ref, "role": m.role} for m in o.members]}
            member_ways.update(m.ref for m in o.members if m.type == "w")
    print(f"pass 1: {len(rels)} candidate relations, {len(member_ways)} member ways")

    # Pass 2: ways with node locations.
    ways_out, area_ways_out, member_geom = [], [], {}
    for o in osmium.FileProcessor(pbf, osmium.osm.NODE | osmium.osm.WAY) \
            .with_locations().with_filter(osmium.filter.EntityFilter(osmium.osm.WAY)):
        tags = dict(o.tags)
        want = "building" in tags or "highway" in tags
        want_area = is_area_tags(tags)
        if not (want or want_area or o.id in member_ways):
            continue
        geom, missing = [], 0
        for n in o.nodes:
            if n.location.valid():
                geom.append({"lat": n.location.lat, "lon": n.location.lon})
            else:
                missing += 1
        if missing:
            stats["ways_with_missing_node_locations"] += 1
        if o.id in member_ways:
            member_geom[o.id] = geom
        el = {"type": "way", "id": o.id, "tags": tags, "geometry": geom}
        if want and any(inside(p["lon"], p["lat"], CLIP_BOX) for p in geom):
            ways_out.append(el)
        if want_area and any(inside(p["lon"], p["lat"], SRINAGAR_BBOX) for p in geom):
            area_ways_out.append(el)
    print(f"pass 2: {len(ways_out)} building/highway ways in clip box")

    # Pass 3: tagged nodes.
    nodes_out, area_nodes_out = [], []
    for o in osmium.FileProcessor(pbf, osmium.osm.NODE) \
            .with_filter(osmium.filter.KeyFilter("building", "place")):
        tags = dict(o.tags)
        if not o.location.valid():
            continue
        el = {"type": "node", "id": o.id, "lat": o.location.lat, "lon": o.location.lon, "tags": tags}
        if "building" in tags and inside(el["lon"], el["lat"], CLIP_BOX):
            nodes_out.append(el)
        if tags.get("place") in PLACES and inside(el["lon"], el["lat"], SRINAGAR_BBOX):
            area_nodes_out.append(el)

    # Relations: attach member way geometry, keep by box rule.
    rels_out, area_rels_out = [], []
    for r in rels.values():
        pts = []
        for m in r["members"]:
            if m["type"] == "way":
                g = member_geom.get(m["ref"])
                if g is None:
                    stats["relation_members_missing_from_file"] += 1
                    continue
                m["geometry"] = g
                pts.extend(g)
        if "building" in r["tags"] and any(inside(p["lon"], p["lat"], CLIP_BOX) for p in pts):
            rels_out.append(r)
        if is_area_tags(r["tags"]) and any(inside(p["lon"], p["lat"], SRINAGAR_BBOX) for p in pts):
            area_rels_out.append(r)

    def dump(name, els):
        raw = json.dumps({"source": rec["file"], "osm_timestamp": rec["osm_timestamp_header"],
                          "elements": els}).encode()
        (OUT / name).write_bytes(raw)
        return {"file": f"data/osm/exposure/{name}", "elements": len(els),
                "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest(),
                "osm_base_timestamp": rec["osm_timestamp_header"], "tile_deg": None,
                "downloaded_utc": rec["downloaded_utc"]}

    main_rec = dump("extract_elements.json", nodes_out + ways_out + rels_out)
    area_rec = dump("srinagar_areas.json", area_nodes_out + area_ways_out + area_rels_out)
    manifest = {"source": "geofabrik", "extract": rec, "clip_box": CLIP_BOX,
                "started_utc": rec["downloaded_utc"], "finished_utc": rec["downloaded_utc"],
                "tiles": [main_rec], "areas": area_rec, "read_stats": dict(stats),
                "counts": {"building_nodes": sum("building" in e["tags"] for e in nodes_out),
                           "ways": len(ways_out), "building_relations": len(rels_out),
                           "area_elements": len(area_nodes_out) + len(area_ways_out) + len(area_rels_out)}}
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(json.dumps({k: manifest[k] for k in ("read_stats", "counts")}, indent=2))


if __name__ == "__main__":
    main()
