"""Step 3a: download QA_RADSAT for the four scenes and put it on the mosaic grid.

QA_RADSAT marks pixels where a band was radiometrically saturated (and terrain
occlusion). It is found by its STAC title, downloaded with the same checks as
the other bands (01_download.py), and placed on the mosaic grid by taking each
pixel from the same row the reflectance mosaic used (data/clipped/*_source.tif),
so the saturation flags belong to the pixel values actually used.

QA_RADSAT declares no fill value (neither in the file nor in STAC). Where the
mosaic has no data (source = 0) the output holds 0; those pixels are already
fill in QA_PIXEL. The record is data/radsat_report.json.

Usage:  .venv/bin/python scripts/08_radsat.py
"""

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import planetary_computer
import pystac_client
import rasterio

ROOT = Path(__file__).resolve().parent.parent
CLIP_DIR = ROOT / "data" / "clipped"
TITLE = "Radiometric Saturation and Terrain Occlusion Quality Assessment Band"
DATES = {"pre": "2014-08-25", "flood": "2014-09-10"}


def load(name, file):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).parent / file)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


dl = load("dl", "01_download.py")
mo = load("mo", "01b_mosaic.py")


def main():
    meta = json.loads((ROOT / "data" / "metadata.json").read_text())
    catalog = pystac_client.Client.open(dl.STAC_URL)
    report = {"title": TITLE, "scenes": {}, "mosaics": {}}

    for name, sc in meta["scenes"].items():
        item = dl.find_scene(catalog, sc["wrs_row"], sc["date"])
        if item.id != sc["scene_id"]:
            sys.exit(f"STOP: {name}: search returned {item.id}, step 1 used {sc['scene_id']}")
        key = dl.find_asset_key(item, {"title": TITLE})
        asset = planetary_computer.sign(item).assets[key]
        raw_name = Path(asset.href.split("?")[0]).name
        raw_path = dl.RAW_DIR / raw_name
        print(f"[{name}] {item.id} -> asset '{key}'")
        check = dl.download(asset.href, raw_path)
        with rasterio.open(raw_path) as src:
            grid = {"crs": src.crs.to_string(), "transform": list(src.transform)[:6],
                    "width": src.width, "height": src.height, "dtype": src.dtypes[0]}
            nodata = src.nodata
        ref = sc["bands"]["green"]["grid"]
        if grid["crs"] != ref["crs"] or grid["transform"] != ref["transform"] \
                or (grid["width"], grid["height"]) != (ref["width"], ref["height"]):
            sys.exit(f"STOP: {raw_name}: grid differs from the scene's green band")
        bitfields = item.assets[key].extra_fields.get("classification:bitfields", [])
        report["scenes"][name] = {
            "scene_id": item.id, "asset_key": key, "raw_file": f"data/raw/{raw_name}",
            **check, "grid": grid, "fill_in_file_header": nodata,
            "fill_in_stac": (item.assets[key].extra_fields.get("raster:bands") or [{}])[0].get("nodata"),
            "bitfields": [{"offset": b["offset"], "name": b["name"],
                           "description": b.get("description")} for b in bitfields]}

    for label, date in DATES.items():
        with rasterio.open(CLIP_DIR / f"{label}_{date}_source.tif") as src:
            source = src.read(1)
            profile = src.profile
        h, w = source.shape
        out = np.zeros((h, w), dtype=np.uint16)
        for code, row in ((1, "036"), (2, "037")):
            sc = report["scenes"][f"{label}_r{row}"]
            arr = mo.read_on_grid(ROOT / sc["raw_file"], profile["transform"], w, h, 0)
            out[source == code] = arr[source == code]
        path = CLIP_DIR / f"{label}_{date}_qa_radsat.tif"
        p = {**profile, "dtype": "uint16", "nodata": None}
        with rasterio.open(path, "w", **p) as dst:
            dst.write(out, 1)
            dst.update_tags(VALUES="QA_RADSAT bit field, copied from the row in *_source.tif; "
                                   "0 where source = 0 (no data)")
        report["mosaics"][date] = {"file": f"data/clipped/{path.name}",
                                   "pixels_green_saturated": int(((out >> 2) & 1).sum()),
                                   "pixels_swir1_saturated": int(((out >> 5) & 1).sum()),
                                   "pixels_terrain_occluded": int(((out >> 11) & 1).sum())}
        print(f"  {date}: {report['mosaics'][date]}")

    (ROOT / "data" / "radsat_report.json").write_text(json.dumps(report, indent=2))
    print("wrote data/radsat_report.json")


if __name__ == "__main__":
    main()
