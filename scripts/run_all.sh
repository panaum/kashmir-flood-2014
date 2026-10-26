#!/usr/bin/env bash
# Run the full analysis in order and stop on the first error.
# Script 07 (write_methods_md) runs near the end because it reads all
# analysis outputs; it cannot run until steps 08-17 have completed.
set -euo pipefail
P=.venv/bin/python
$P scripts/01_download.py             # Landsat 8 scenes
$P scripts/01b_mosaic.py              # rows 036 + 037 mosaic
$P scripts/02_preview.py              # cloud check, previews
$P scripts/03_write_data_md.py        # DATA.md
$P scripts/04_dem.py                  # Copernicus DEM
$P scripts/05_water.py                # MNDWI, masks, Otsu, water masks
$P scripts/06_water_preview.py        # water previews
$P scripts/08_radsat.py               # saturation band
$P scripts/09_osm.py                  # OSM lakes and wetlands
$P scripts/10_flood.py                # flood extent
$P scripts/11_flood_preview.py        # flood previews
$P scripts/12_geofabrik_download.py   # OSM extract
$P scripts/12b_read_extract.py        # OSM buildings, roads, places
$P scripts/13_exposure.py             # OSM counts
$P scripts/14_exposure_preview.py     # OSM preview
$P scripts/15_msbuildings_download.py # Microsoft footprints
$P scripts/16_msbuildings_count.py    # footprint counts
$P scripts/17_msbuildings_preview.py  # footprint preview
$P scripts/07_write_methods_md.py     # METHODS.md (reads all outputs above)
$P scripts/18_site_data.py            # website data and layers
$P scripts/19_write_readme.py         # README.md
