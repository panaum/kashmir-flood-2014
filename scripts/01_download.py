"""Step 1a: find and download the Landsat 8 scenes.

Searches the Microsoft Planetary Computer STAC API for Landsat 8 Collection 2
Level-2 scenes over WRS path 149, rows 036 and 037, on the two target dates,
and downloads the full original band files to data/raw/ unchanged. Everything
later steps and DATA.md need (IDs, times, cloud cover, grid, scale/offset,
fill values, checksums) is written to data/metadata.json.

Checks, each of which stops the script instead of guessing:
  - exactly one Landsat 8 item per row and date
  - each band found by its meaning (eo:bands common name / asset title),
    not by an assumed asset key, and matching exactly one asset
  - downloads must match the server's Content-Length (and Content-MD5 when
    the server provides one); files already on disk are re-checked against
    the server the same way before being reused
  - STAC scale/offset must equal the scene's own MTL metadata
  - a fill value in the file header must equal the STAC fill value

Usage:  .venv/bin/python scripts/01_download.py
"""

import base64
import hashlib
import json
import sys
from pathlib import Path

import planetary_computer
import pystac_client
import rasterio
import requests

STAC_URL = "https://planetarycomputer.microsoft.com/api/stac/v1"
COLLECTION = "landsat-c2-l2"
PLATFORM = "landsat-8"
WRS_PATH = "149"
WRS_ROWS = ["036", "037"]
DATES = {"pre": "2014-08-25", "flood": "2014-09-10"}

# What we want, described by meaning rather than by asset key. The key is
# looked up from each item's own asset metadata (see find_asset_key).
# Note: the asset keyed "qa" is the surface *temperature* QA band, not pixel QA.
WANTED = {
    "green": {"eo_common_name": "green"},   # OLI band 3
    "red": {"eo_common_name": "red"},       # OLI band 4
    "nir": {"eo_common_name": "nir08"},     # OLI band 5
    "swir1": {"eo_common_name": "swir16"},  # OLI band 6
    "qa_pixel": {"title": "Pixel Quality Assessment Band"},
}

ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = ROOT / "data" / "raw"


def find_scene(catalog, row, date):
    search = catalog.search(
        collections=[COLLECTION],
        datetime=date,
        query={
            "platform": {"eq": PLATFORM},
            "landsat:wrs_path": {"eq": WRS_PATH},
            "landsat:wrs_row": {"eq": row},
        },
    )
    items = list(search.items())
    if len(items) != 1:
        sys.exit(f"STOP: expected exactly 1 Landsat 8 scene for path {WRS_PATH} "
                 f"row {row} on {date}, found {len(items)}: {[i.id for i in items]}")
    return items[0]


def find_asset_key(item, spec):
    """Return the single asset key matching spec, or stop if 0 or >1 match."""
    matches = []
    for key, asset in item.assets.items():
        if "eo_common_name" in spec:
            names = [b.get("common_name") for b in asset.extra_fields.get("eo:bands", [])]
            if spec["eo_common_name"] in names:
                matches.append(key)
        elif asset.title == spec["title"]:
            matches.append(key)
    if len(matches) != 1:
        sys.exit(f"STOP: {item.id}: asset lookup {spec} matched {matches}")
    return matches[0]


def file_md5_b64(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return base64.b64encode(h.digest()).decode()


def file_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify(path, size, md5):
    """Stop unless path matches the server's size (and MD5 if given)."""
    if size is not None and path.stat().st_size != size:
        sys.exit(f"STOP: {path.name}: {path.stat().st_size} bytes on disk, "
                 f"server says {size}")
    if md5 is not None and file_md5_b64(path) != md5:
        sys.exit(f"STOP: {path.name}: MD5 does not match server Content-MD5")


def download(url, dest):
    """Download url to dest, or re-verify an existing dest against the server.

    Returns a dict with sha256, size and how the file was verified.
    """
    if dest.exists():
        r = requests.head(url, timeout=60)
        r.raise_for_status()
        size = int(r.headers["Content-Length"]) if "Content-Length" in r.headers else None
        md5 = r.headers.get("Content-MD5")
        verify(dest, size, md5)
        print(f"    exists, re-verified {dest.relative_to(ROOT)}")
    else:
        tmp = dest.with_suffix(dest.suffix + ".part")
        with requests.get(url, stream=True, timeout=120) as r:
            r.raise_for_status()
            size = int(r.headers["Content-Length"]) if "Content-Length" in r.headers else None
            md5 = r.headers.get("Content-MD5")
            with open(tmp, "wb") as f:
                for chunk in r.iter_content(1 << 20):
                    f.write(chunk)
        verify(tmp, size, md5)
        tmp.rename(dest)
        print(f"    downloaded {dest.relative_to(ROOT)}")
    if size is None:
        sys.exit(f"STOP: server gave no Content-Length for {dest.name}; "
                 f"cannot confirm the download is complete")
    return {"bytes": dest.stat().st_size, "sha256": file_sha256(dest),
            "size_verified": True, "md5_verified": md5 is not None}


def mtl_scale_offset(mtl, eo_band_names):
    """Surface reflectance MULT/ADD for an OLI band, from the scene's MTL."""
    if len(eo_band_names) != 1 or not eo_band_names[0].startswith("OLI_B"):
        return None, None
    n = eo_band_names[0].removeprefix("OLI_B")
    p = mtl["LEVEL2_SURFACE_REFLECTANCE_PARAMETERS"]
    return float(p[f"REFLECTANCE_MULT_BAND_{n}"]), float(p[f"REFLECTANCE_ADD_BAND_{n}"])


def main():
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    catalog = pystac_client.Client.open(STAC_URL)
    metadata = {"stac_api": STAC_URL, "collection": COLLECTION, "scenes": {}}

    for label, date in DATES.items():
        for row in WRS_ROWS:
            item = find_scene(catalog, row, date)
            key_name = f"{label}_r{row}"
            print(f"[{key_name}] {item.id}  {item.datetime.isoformat()}")
            print(f"    all asset keys: {sorted(item.assets)}")

            # Save the unsigned STAC item for the record.
            (RAW_DIR / f"{item.id}_stac_item.json").write_text(
                json.dumps(item.to_dict(), indent=2))
            signed = planetary_computer.sign(item)

            mtl_path = RAW_DIR / f"{item.id}_MTL.json"
            mtl_check = download(signed.assets["mtl.json"].href, mtl_path)
            mtl = json.loads(mtl_path.read_text())["LANDSAT_METADATA_FILE"]
            attrs = mtl["IMAGE_ATTRIBUTES"]

            p = item.properties
            scene = {
                "label": label,
                "date": date,
                "wrs_path": WRS_PATH,
                "wrs_row": row,
                "scene_id": item.id,
                "landsat_product_id": mtl["PRODUCT_CONTENTS"]["LANDSAT_PRODUCT_ID"],
                "landsat_scene_id": p.get("landsat:scene_id"),
                "platform": p.get("platform"),
                "datetime_utc_stac": p["datetime"],
                "date_acquired_mtl": attrs["DATE_ACQUIRED"],
                "scene_center_time_mtl": attrs["SCENE_CENTER_TIME"],
                "eo_cloud_cover_stac": p.get("eo:cloud_cover"),
                "cloud_cover_land_stac": p.get("landsat:cloud_cover_land"),
                "cloud_cover_mtl": float(attrs["CLOUD_COVER"]),
                "cloud_cover_land_mtl": float(attrs["CLOUD_COVER_LAND"]),
                "crs_stac": p.get("proj:code"),
                "collection_category": p.get("landsat:collection_category"),
                "mtl_file": f"data/raw/{mtl_path.name}",
                "mtl_sha256": mtl_check["sha256"],
                "bands": {},
            }

            for band, spec in WANTED.items():
                key = find_asset_key(item, spec)
                asset = signed.assets[key]
                raw_name = Path(asset.href.split("?")[0]).name
                raw_path = RAW_DIR / raw_name
                print(f"  {band:9s} -> asset '{key}' ({asset.title})")
                check = download(asset.href, raw_path)

                raster = (asset.extra_fields.get("raster:bands") or [{}])[0]
                eo_names = [b.get("name") for b in asset.extra_fields.get("eo:bands", [])]
                mtl_scale, mtl_offset = mtl_scale_offset(mtl, eo_names)
                if mtl_scale is not None and (mtl_scale != raster.get("scale")
                                              or mtl_offset != raster.get("offset")):
                    sys.exit(f"STOP: {raw_name}: STAC scale/offset "
                             f"{raster.get('scale')}/{raster.get('offset')} != MTL "
                             f"{mtl_scale}/{mtl_offset}")

                with rasterio.open(raw_path) as src:
                    file_nodata = src.nodata
                    grid = {
                        "crs": src.crs.to_string(),
                        "transform": list(src.transform)[:6],
                        "width": src.width,
                        "height": src.height,
                        "dtype": src.dtypes[0],
                    }
                stac_nodata = raster.get("nodata")
                if file_nodata is not None and stac_nodata is not None \
                        and file_nodata != stac_nodata:
                    sys.exit(f"STOP: {raw_name}: file nodata {file_nodata} != "
                             f"STAC nodata {stac_nodata}")
                fill = file_nodata if file_nodata is not None else stac_nodata
                if fill is None:
                    sys.exit(f"STOP: no fill value known for {raw_name}")

                scene["bands"][band] = {
                    "asset_key": key,
                    "asset_title": asset.title,
                    "eo_band_name": eo_names,
                    "raw_file": f"data/raw/{raw_name}",
                    **check,
                    "grid": grid,
                    "scale_stac": raster.get("scale"),
                    "offset_stac": raster.get("offset"),
                    "scale_mtl": mtl_scale,
                    "offset_mtl": mtl_offset,
                    "fill_value": fill,
                    "fill_in_file_header": file_nodata,
                    "fill_in_stac": stac_nodata,
                }
            metadata["scenes"][key_name] = scene

    out = ROOT / "data" / "metadata.json"
    out.write_text(json.dumps(metadata, indent=2))
    print(f"wrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
