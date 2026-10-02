"""Step 4a: download the Geofabrik OpenStreetMap extract (amendment 1, G1-G2).

Reads Geofabrik's extract index, picks the smallest extract (by area of its
published boundary) that contains the whole analysis area, downloads its
.osm.pbf with the same size checks as the Landsat files, verifies it against
Geofabrik's published MD5, and records file name, URL, size, MD5, sha256 and
the OpenStreetMap data timestamp in the file header.

Stops if no extract contains the analysis area, or if the MD5 does not match.
Record: data/osm/geofabrik/extract.json.

Usage:  .venv/bin/python scripts/12_geofabrik_download.py
"""

import hashlib
import importlib.util
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import osmium
import rasterio
import requests
import shapely
from rasterio.features import shapes
from rasterio.warp import transform_geom
from shapely.geometry import shape

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "osm" / "geofabrik"
INDEX_URL = "https://download.geofabrik.de/index-v1.json"
HEADERS = {"User-Agent": "kashmir-flood-2014-landsat-study/1.0 (research script; python-requests)"}

_spec = importlib.util.spec_from_file_location("dl", Path(__file__).parent / "01_download.py")
dl = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dl)


def analysis_area_lonlat():
    with rasterio.open(ROOT / "data" / "flood" / "analysis_area.tif") as src:
        a = src.read(1)
        polys = [shape(transform_geom(src.crs, "EPSG:4326", g))
                 for g, v in shapes(a, mask=a == 1, transform=src.transform) if v == 1]
    return shapely.union_all(polys)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fetched = datetime.now(timezone.utc).isoformat(timespec="seconds")
    r = requests.get(INDEX_URL, headers=HEADERS, timeout=120)
    r.raise_for_status()
    (OUT / "index-v1.json").write_bytes(r.content)
    index = r.json()

    area = analysis_area_lonlat()
    candidates = []
    for f in index["features"]:
        p = f["properties"]
        pbf = p.get("urls", {}).get("pbf")
        if not pbf or f.get("geometry") is None:
            continue
        g = shape(f["geometry"])
        if g.covers(area):
            candidates.append((g.area, p["id"], p.get("name"), pbf, g))
    if not candidates:
        sys.exit("STOP: no Geofabrik extract contains the whole analysis area")
    candidates.sort(key=lambda c: c[0])
    deg2, ext_id, name, url, geom = candidates[0]
    print(f"chosen: {ext_id} ({name}) {url}; {len(candidates)} extracts contain the area: "
          f"{[c[1] for c in candidates[:6]]}")
    # Double check with a small inward tolerance-free test: every vertex inside.
    outside = shapely.difference(area, geom)
    if not outside.is_empty and outside.area > 0:
        sys.exit(f"STOP: part of the analysis area ({outside.area} deg2) is outside {ext_id}")

    md5_txt = requests.get(url + ".md5", headers=HEADERS, timeout=60)
    md5_txt.raise_for_status()
    published_md5 = md5_txt.text.split()[0].strip().lower()
    dest = OUT / Path(url).name
    check = dl.download(url, dest)
    h = hashlib.md5()
    with open(dest, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    if h.hexdigest() != published_md5:
        sys.exit(f"STOP: {dest.name} MD5 {h.hexdigest()} != published {published_md5}")
    header = osmium.io.Reader(str(dest), osmium.osm.osm_entity_bits.NOTHING).header()
    record = {
        "index_url": INDEX_URL, "index_fetched_utc": fetched,
        "index_sha256": hashlib.sha256(r.content).hexdigest(),
        "extract_id": ext_id, "extract_name": name, "url": url,
        "extracts_containing_area": [c[1] for c in candidates],
        "boundary_area_deg2": round(deg2, 4),
        "file": f"data/osm/geofabrik/{dest.name}", **check,
        "md5_published": published_md5, "md5_computed": h.hexdigest(),
        "osm_timestamp_header": header.get("osmosis_replication_timestamp"),
        "osm_replication_sequence": header.get("osmosis_replication_sequence_number"),
        "header_generator": header.get("generator"),
        "downloaded_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "analysis_area_outside_boundary_deg2": 0.0,
    }
    (OUT / "extract.json").write_text(json.dumps(record, indent=2))
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
