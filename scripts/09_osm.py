"""Step 3b: OpenStreetMap water bodies and wetlands for the clip box.

Queries the Overpass API for ways and relations tagged natural=water,
natural=wetland, landuse=reservoir or waterway=riverbank in the clip box
(rule O1), saves the raw response with its download time, builds polygons
(closed ways; multipolygon relations from their outer/inner members), and
rasterises them onto the mosaic grid: a pixel is inside if its centre is
inside a polygon (rule O2). Stops if the download fails (rule O3).

Note: OSM data are current (download date), not from 2014.

Outputs: data/osm/overpass_response.json, data/osm/water_wetland_mask.tif
(1 = inside, 0 = outside), and the record data/osm_report.json.

Usage:  .venv/bin/python scripts/09_osm.py
"""

import hashlib
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import rasterio
import requests
from rasterio.features import rasterize
from rasterio.warp import transform_geom
from shapely.geometry import LineString, Polygon, mapping
from shapely.ops import linemerge, polygonize, unary_union

ROOT = Path(__file__).resolve().parent.parent
OSM_DIR = ROOT / "data" / "osm"
GRID_REF = ROOT / "data" / "water" / "valley_floor.tif"
ENDPOINT = "https://overpass-api.de/api/interpreter"
TAGS = [("natural", "water"), ("natural", "wetland"),
        ("landuse", "reservoir"), ("waterway", "riverbank")]
ATTEMPTS = 3
# Overpass rejects the generic python-requests client name (HTTP 406).
HEADERS = {"User-Agent": "kashmir-flood-2014-landsat-study/1.0 (research script; python-requests)",
           "Accept": "application/json"}


def query(bbox):
    w, s, e, n = bbox
    b = f"({s},{w},{n},{e})"
    parts = "".join(f'way["{k}"="{v}"]{b};relation["{k}"="{v}"]{b};' for k, v in TAGS)
    return f"[out:json][timeout:300];({parts});out geom;"


def fetch(q):
    last = None
    for attempt in range(1, ATTEMPTS + 1):
        try:
            r = requests.post(ENDPOINT, data={"data": q}, headers=HEADERS, timeout=400)
            if r.status_code == 200:
                return r
            last = f"HTTP {r.status_code}: {r.text[:300]}"
        except requests.RequestException as e:
            last = repr(e)
        print(f"  attempt {attempt} failed: {last}")
        time.sleep(30 * attempt)
    sys.exit(f"STOP: Overpass download failed after {ATTEMPTS} attempts: {last}")


def way_coords(el):
    return [(p["lon"], p["lat"]) for p in el.get("geometry", []) if p]


def element_polygon(el):
    """Polygon for a closed way or a multipolygon relation, else None."""
    if el["type"] == "way":
        c = way_coords(el)
        if len(c) >= 4 and c[0] == c[-1]:
            return Polygon(c).buffer(0)
        return None
    if el.get("tags", {}).get("type") not in ("multipolygon", None):
        return None
    rings = {"outer": [], "inner": []}
    for m in el.get("members", []):
        if m.get("type") != "way" or "geometry" not in m:
            continue
        c = [(p["lon"], p["lat"]) for p in m["geometry"] if p]
        if len(c) >= 2:
            rings["inner" if m.get("role") == "inner" else "outer"].append(LineString(c))
    if not rings["outer"]:
        return None
    outer = unary_union(list(polygonize(linemerge(rings["outer"]))))
    if rings["inner"]:
        inner = unary_union(list(polygonize(linemerge(rings["inner"]))))
        outer = outer.difference(inner)
    return outer.buffer(0) if not outer.is_empty else None


def main():
    OSM_DIR.mkdir(parents=True, exist_ok=True)
    bbox = json.loads((ROOT / "data" / "mosaic_report.json").read_text())["bbox_wgs84"]
    q = query(bbox)
    started = datetime.now(timezone.utc).isoformat(timespec="seconds")
    r = fetch(q)
    raw = r.content
    data = r.json()
    (OSM_DIR / "overpass_response.json").write_bytes(raw)
    els = data.get("elements", [])
    if not els:
        sys.exit("STOP: Overpass returned no elements for the clip box")

    polys, counts, unbuilt = [], {}, {"way": 0, "relation": 0}
    for el in els:
        tags = el.get("tags", {})
        key = next((f"{k}={v}" for k, v in TAGS if tags.get(k) == v), "other")
        counts.setdefault(key, {"way": 0, "relation": 0, "polygons": 0})
        counts[key][el["type"]] += 1
        p = element_polygon(el)
        if p is None or p.is_empty:
            unbuilt[el["type"]] += 1
            continue
        counts[key]["polygons"] += 1
        polys.append(p)
    union = unary_union(polys)

    with rasterio.open(GRID_REF) as ref:
        transform, crs, shape, profile = ref.transform, ref.crs, (ref.height, ref.width), ref.profile
    geom = transform_geom("EPSG:4326", crs, mapping(union))
    mask = rasterize([(geom, 1)], out_shape=shape, transform=transform, fill=0, dtype="uint8")
    with rasterio.open(OSM_DIR / "water_wetland_mask.tif", "w",
                       **{**profile, "dtype": "uint8", "nodata": None}) as dst:
        dst.write(mask, 1)
        dst.update_tags(SOURCE=f"OpenStreetMap via {ENDPOINT}, downloaded {started}",
                        TAGS=", ".join(f"{k}={v}" for k, v in TAGS),
                        VALUES="1 = pixel centre inside a polygon")

    report = {
        "endpoint": ENDPOINT, "user_agent": HEADERS["User-Agent"], "query": q, "download_started_utc": started,
        "osm_base_timestamp": data.get("osm3s", {}).get("timestamp_osm_base"),
        "response_bytes": len(raw), "response_sha256": hashlib.sha256(raw).hexdigest(),
        "elements": len(els), "by_tag": counts,
        "elements_not_turned_into_polygons": unbuilt,
        "polygon_area_deg2": round(union.area, 8),
        "pixels_inside": int(mask.sum()), "km2_inside": round(int(mask.sum()) * 0.0009, 4),
        "licence": "OpenStreetMap data, Open Database License (ODbL); "
                   "attribution: (c) OpenStreetMap contributors",
    }
    (ROOT / "data" / "osm_report.json").write_text(json.dumps(report, indent=2))
    print(json.dumps({k: v for k, v in report.items() if k != "query"}, indent=2))


if __name__ == "__main__":
    main()
