"""Step 4b-a: download Microsoft building footprints (rules M1-M2).

Searches the Planetary Computer for `ms-buildings` items intersecting the clip
box, keeps only the newest `msbuildings:processing-date` (the dataset
version), checks that the chosen items together contain the whole analysis
area, and downloads every active parquet file of their partitions. Active
files are taken from the Delta table's transaction log (files added and not
later removed), not from a directory listing, and the two are compared.

Each file is checked against the server's size (and Content-MD5 where given)
and a sha256 is computed. Every file is opened with pyarrow and its row count
recorded, so a file that cannot be read completely stops the work.

Record: data/msbuildings/record.json. Files: data/msbuildings/<quadkey>/.

TLS: this Python build has no default CA bundle, so SSL_CERT_FILE is pointed
at certifi's bundle (the one `requests` already uses). Verification stays on.

Usage:  .venv/bin/python scripts/15_msbuildings_download.py
"""

import os

import certifi

os.environ.setdefault("SSL_CERT_FILE", certifi.where())

import base64
import hashlib
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import adlfs
import planetary_computer
import pyarrow.parquet as pq
import pystac_client
import shapely
from shapely.geometry import shape

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "msbuildings"
STAC_URL = "https://planetarycomputer.microsoft.com/api/stac/v1"
COLLECTION = "ms-buildings"
CLIP_BOX = [74.2, 33.5, 75.4, 34.6]


def stop(msg):
    sys.exit(f"STOP: {msg}")


def retry(fn, what, attempts=4):
    for k in range(1, attempts + 1):
        try:
            return fn()
        except Exception as e:  # noqa: BLE001
            print(f"  {what}: attempt {k} failed: {type(e).__name__}: {str(e)[:150]}")
            if k == attempts:
                stop(f"{what} failed after {attempts} attempts")
            time.sleep(15 * k)


def analysis_area_lonlat():
    import importlib.util
    spec = importlib.util.spec_from_file_location("gf", Path(__file__).parent / "12_geofabrik_download.py")
    gf = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gf)
    return gf.analysis_area_lonlat()


def delta_active_files(fs, table_root):
    """Replay the Delta log (checkpoint + later JSON commits) -> active paths."""
    log = f"{table_root}/_delta_log"
    names = sorted(p.split("/")[-1] for p in fs.ls(log))
    active, start = {}, -1
    if "_last_checkpoint" in names:
        last = json.loads(fs.cat_file(f"{log}/_last_checkpoint"))
        start = int(last["version"])
        parts = sorted(n for n in names if n.startswith(f"{start:020d}.checkpoint"))
        for n in parts:
            t = pq.read_table(fs.open(f"{log}/{n}", "rb")).to_pylist()
            for row in t:
                if row.get("add"):
                    active[row["add"]["path"]] = row["add"]
    commits = sorted(n for n in names if n.endswith(".json") and int(n[:20]) > start)
    for n in commits:
        for line in fs.cat_file(f"{log}/{n}").decode().splitlines():
            if not line.strip():
                continue
            act = json.loads(line)
            if "add" in act:
                active[act["add"]["path"]] = act["add"]
            elif "remove" in act:
                active.pop(act["remove"]["path"], None)
    return active, {"checkpoint_version": start, "json_commits_replayed": len(commits),
                    "log_files": len(names)}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    started = datetime.now(timezone.utc).isoformat(timespec="seconds")
    cat = retry(lambda: pystac_client.Client.open(STAC_URL), "open STAC API")
    coll = retry(lambda: cat.get_collection(COLLECTION), f"collection {COLLECTION}")
    items = retry(lambda: list(cat.search(collections=[COLLECTION], bbox=CLIP_BOX).items()), "search")
    dated = [i for i in items if i.properties.get("msbuildings:processing-date")]
    if not dated:
        stop("no partitioned ms-buildings items intersect the clip box")
    version = max(i.properties["msbuildings:processing-date"] for i in dated)
    chosen = sorted((i for i in dated if i.properties["msbuildings:processing-date"] == version),
                    key=lambda i: i.id)
    print(f"version {version}: {[i.id for i in chosen]}")

    area = analysis_area_lonlat()
    union = shapely.union_all([shape(i.geometry) for i in chosen])
    missing = shapely.difference(area, union)
    if not missing.is_empty and missing.area > 0:
        stop(f"chosen items do not contain the analysis area ({missing.area:.6f} deg2 outside)")

    # Read token for the storage container named in the collection metadata.
    account = coll.extra_fields["msft:storage_account"]
    container = coll.extra_fields["msft:container"]
    token = retry(lambda: planetary_computer.sas.get_token(account, container).token, "SAS token")
    table_root = coll.assets["delta"].href.replace("az://", "").rstrip("/")
    fs = adlfs.AzureBlobFileSystem(account_name=account, credential=token)
    active, log_info = retry(lambda: delta_active_files(fs, table_root), "read Delta log")
    print(f"Delta log: {log_info}; {len(active)} active files in the table")

    record = {"stac_api": STAC_URL, "collection": COLLECTION, "collection_title": coll.title,
              "license": coll.license, "dataset_version_processing_date": version,
              "delta_table": table_root, "delta_log": log_info, "started_utc": started,
              "tls_ca_bundle": "certifi " + __import__("importlib.metadata").metadata.version("certifi"),
              "items": []}
    for it in chosen:
        p = it.properties
        q = str(p["msbuildings:quadkey"])
        region = p["msbuildings:region"]
        prefix = f"RegionName={region}/quadkey={q}/"
        files_log = sorted(k for k in active if k.startswith(prefix))
        listing = retry(lambda: fs.find(f"{table_root}/{prefix}", detail=True), f"list {q}")
        files_dir = sorted(k.split(f"{table_root}/", 1)[1] for k in listing)
        if not files_log:
            stop(f"quadkey {q}: no active files in the Delta log")
        qdir = OUT / q
        qdir.mkdir(exist_ok=True)
        frec = []
        for rel in files_log:
            remote = f"{table_root}/{rel}"
            info = retry(lambda: fs.info(remote), f"info {rel}")
            cs = info.get("content_settings")
            md5 = getattr(cs, "content_md5", None) if cs is not None else None
            md5_b64 = base64.b64encode(bytes(md5)).decode() if md5 else None
            dest = qdir / rel.split("/")[-1]
            if not (dest.exists() and dest.stat().st_size == info["size"]):
                tmp = dest.with_suffix(".part")
                retry(lambda: fs.get(remote, str(tmp)), f"download {rel}")
                if tmp.stat().st_size != info["size"]:
                    stop(f"{rel}: {tmp.stat().st_size} bytes, server says {info['size']}")
                tmp.rename(dest)
            raw = dest.read_bytes()
            if md5_b64 and base64.b64encode(hashlib.md5(raw).digest()).decode() != md5_b64:
                stop(f"{rel}: MD5 does not match server Content-MD5")
            try:
                meta = pq.read_metadata(dest)
                table = pq.read_table(dest)
            except Exception as e:  # noqa: BLE001
                stop(f"{rel}: cannot be read completely: {e}")
            if table.num_rows != meta.num_rows:
                stop(f"{rel}: read {table.num_rows} rows, metadata says {meta.num_rows}")
            frec.append({"path": rel, "local": str(dest.relative_to(ROOT)), "bytes": info["size"],
                         "content_md5": md5_b64, "md5_verified": md5_b64 is not None,
                         "sha256": hashlib.sha256(raw).hexdigest(), "rows": table.num_rows,
                         "columns": table.column_names,
                         "log_num_records": (active[rel].get("stats") and
                                             json.loads(active[rel]["stats"]).get("numRecords"))})
            print(f"  {rel}: {table.num_rows} rows")
        record["items"].append({
            "id": it.id, "region": region, "quadkey": q,
            "start_datetime": p.get("start_datetime"), "end_datetime": p.get("end_datetime"),
            "table_row_count_property": p.get("table:row_count"),
            "columns": [c["name"] for c in p.get("table:columns", [])],
            "bbox": it.bbox, "files": frec,
            "files_in_directory_not_in_log": sorted(set(files_dir) - set(files_log)),
            "files_in_log_not_in_directory": sorted(set(files_log) - set(files_dir)),
            "rows": sum(f["rows"] for f in frec)})
        for f in frec:
            if f["log_num_records"] is not None and f["log_num_records"] != f["rows"]:
                stop(f"{f['path']}: {f['rows']} rows, Delta log says {f['log_num_records']}")
    record["finished_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    record["rows_total"] = sum(i["rows"] for i in record["items"])
    (OUT / "record.json").write_text(json.dumps(record, indent=2))
    print(f"{record['rows_total']} footprint rows in {sum(len(i['files']) for i in record['items'])} files")
    for i in record["items"]:
        if i["files_in_directory_not_in_log"] or i["files_in_log_not_in_directory"]:
            print(f"  note {i['quadkey']}: dir-only {len(i['files_in_directory_not_in_log'])}, "
                  f"log-only {len(i['files_in_log_not_in_directory'])}")


if __name__ == "__main__":
    main()
