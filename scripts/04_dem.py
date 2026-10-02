"""Step 2a: download the Copernicus DEM GLO-30 and put it on the Landsat grid.

Searches the same Planetary Computer STAC API for collection cop-dem-glo-30
tiles intersecting the Kashmir valley clip box, downloads each tile unchanged
to data/raw/dem/ with the same checks as the Landsat files (01_download.py),
joins the tiles on their shared 1-arcsecond grid (a straight array copy), and
resamples the result onto the Landsat mosaic grid with bilinear interpolation.

This is the only resampled layer in the project. The output is
data/clipped/dem_cop30_bilinear.tif; the record is data/dem_report.json.

Usage:  .venv/bin/python scripts/04_dem.py
"""

import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
import planetary_computer
import pystac_client
import rasterio
from rasterio.transform import Affine
from rasterio.warp import Resampling, reproject

ROOT = Path(__file__).resolve().parent.parent
STAC_URL = "https://planetarycomputer.microsoft.com/api/stac/v1"
COLLECTION = "cop-dem-glo-30"
DEM_DIR = ROOT / "data" / "raw" / "dem"
OUT = ROOT / "data" / "clipped" / "dem_cop30_bilinear.tif"
GRID_REF = ROOT / "data" / "clipped" / "pre_2014-08-25_green.tif"
DST_NODATA = -9999.0

# Same download/verification code as the Landsat files.
_spec = importlib.util.spec_from_file_location("dl", Path(__file__).parent / "01_download.py")
dl = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dl)


def stop(msg):
    sys.exit(f"STOP: {msg}")


def data_asset_key(item):
    keys = [k for k, a in item.assets.items()
            if "data" in (a.roles or []) and (a.media_type or "").startswith("image/tiff")]
    if len(keys) != 1:
        stop(f"{item.id}: expected one GeoTIFF data asset, found {keys}")
    return keys[0]


def main():
    DEM_DIR.mkdir(parents=True, exist_ok=True)
    mosaic = json.loads((ROOT / "data" / "mosaic_report.json").read_text())
    bbox = mosaic["bbox_wgs84"]

    catalog = pystac_client.Client.open(STAC_URL)
    try:
        coll = catalog.get_collection(COLLECTION)
    except Exception as e:  # noqa: BLE001 - any failure means "not there"
        stop(f"collection {COLLECTION} not available: {e}")
    items = sorted(catalog.search(collections=[COLLECTION], bbox=bbox).items(),
                   key=lambda i: i.id)
    if not items:
        stop(f"no {COLLECTION} tiles intersect {bbox}")
    print(f"{len(items)} tiles: {[i.id for i in items]}")

    report = {"stac_api": STAC_URL, "collection": COLLECTION,
              "collection_title": coll.title,
              "collection_description_first_sentence": coll.description.split(". ")[0] + ".",
              "license": coll.license,
              "license_links": [l.href for l in coll.links if l.rel == "license"],
              "providers": [p.name for p in coll.providers or []],
              "bbox_wgs84": bbox, "tiles": []}

    tiles = []
    for item in items:
        key = data_asset_key(item)
        signed = planetary_computer.sign(item)
        name = Path(signed.assets[key].href.split("?")[0]).name
        path = DEM_DIR / name
        print(f"  {item.id} -> asset '{key}'")
        check = dl.download(signed.assets[key].href, path)
        (DEM_DIR / f"{item.id}_stac_item.json").write_text(json.dumps(item.to_dict(), indent=2))
        with rasterio.open(path) as src:
            t = src.transform
            info = {"crs": src.crs.to_string(), "transform": list(t)[:6],
                    "width": src.width, "height": src.height,
                    "dtype": src.dtypes[0], "nodata_in_file_header": src.nodata,
                    "crs_wkt": src.crs.to_wkt(), "tags": src.tags()}
            arr = src.read(1)
        tiles.append((path, t, arr, src.nodata))
        report["tiles"].append({"item_id": item.id, "asset_key": key,
                                "item_datetime": item.properties.get("datetime"),
                                "file": f"data/raw/dem/{name}", **check, **info,
                                "min": float(arr.min()), "max": float(arr.max())})

    # Join tiles on their shared grid; stop if they do not share one.
    t0 = tiles[0][1]
    res_x, res_y = t0.a, t0.e
    crs0 = report["tiles"][0]["crs"]
    nodatas = {nd for *_, nd in tiles}
    if len(nodatas) != 1:
        stop(f"DEM tiles declare different nodata values: {nodatas}")
    src_nodata = nodatas.pop()
    for (path, t, arr, _), rec in zip(tiles, report["tiles"]):
        if rec["crs"] != crs0 or t.a != res_x or t.e != res_y or t.b or t.d:
            stop(f"{path.name}: grid differs from the first tile")
    west = min(t.c for _, t, _, _ in tiles)
    north = max(t.f for _, t, _, _ in tiles)
    east = max(t.c + a.shape[1] * res_x for _, t, a, _ in tiles)
    south = min(t.f + a.shape[0] * res_y for _, t, a, _ in tiles)
    width = round((east - west) / res_x)
    height = round((south - north) / res_y)
    fill = src_nodata if src_nodata is not None else np.nan
    joined = np.full((height, width), fill, dtype=np.float32)
    for path, t, arr, _ in tiles:
        col, row = (t.c - west) / res_x, (t.f - north) / res_y
        if not (math.isclose(col, round(col), abs_tol=1e-6) and
                math.isclose(row, round(row), abs_tol=1e-6)):
            stop(f"{path.name}: not on a whole-pixel offset from the joined grid")
        col, row = round(col), round(row)
        joined[row:row + arr.shape[0], col:col + arr.shape[1]] = arr
    if np.isnan(joined).any():
        stop("gap between DEM tiles inside their joined extent")
    joined_transform = Affine(res_x, 0, west, 0, res_y, north)

    # Resample (bilinear) onto the Landsat mosaic grid.
    with rasterio.open(GRID_REF) as ref:
        dst_transform, dst_crs = ref.transform, ref.crs
        dst_shape = (ref.height, ref.width)
    dem = np.full(dst_shape, DST_NODATA, dtype=np.float32)
    reproject(source=joined, destination=dem,
              src_transform=joined_transform, src_crs=crs0,
              src_nodata=src_nodata, dst_transform=dst_transform, dst_crs=dst_crs,
              dst_nodata=DST_NODATA, resampling=Resampling.bilinear)
    empty = dem == DST_NODATA
    profile = {"driver": "GTiff", "width": dst_shape[1], "height": dst_shape[0],
               "count": 1, "dtype": "float32", "crs": dst_crs,
               "transform": dst_transform, "nodata": DST_NODATA,
               "compress": "deflate", "tiled": True, "blockxsize": 256, "blockysize": 256}
    with rasterio.open(OUT, "w", **profile) as dst:
        dst.write(dem, 1)
        dst.update_tags(SOURCE=", ".join(i.id for i in items),
                        RESAMPLING="bilinear from EPSG:4326 1 arcsec to Landsat 30 m UTM grid",
                        UNITS="metres")
    v = dem[~empty]
    report["joined_grid"] = {"crs": crs0, "transform": list(joined_transform)[:6],
                             "width": width, "height": height}
    report["resampled"] = {
        "file": f"data/clipped/{OUT.name}", "resampling": "bilinear",
        "dst_crs": str(dst_crs), "dst_transform": list(dst_transform)[:6],
        "dst_width": dst_shape[1], "dst_height": dst_shape[0],
        "dst_nodata": DST_NODATA, "nodata_pixels": int(empty.sum()),
        "min_m": round(float(v.min()), 2), "max_m": round(float(v.max()), 2),
        "median_m": round(float(np.median(v)), 2),
        "src_nodata": src_nodata}
    print(json.dumps(report["resampled"], indent=2))
    (ROOT / "data" / "dem_report.json").write_text(json.dumps(report, indent=2))
    print("wrote data/dem_report.json")


if __name__ == "__main__":
    main()
