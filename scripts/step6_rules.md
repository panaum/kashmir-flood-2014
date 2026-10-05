## Step 6 rules: Jhelum centreline and distance analysis (fixed before any distance result)

These rules were written and committed before the Jhelum geometry script or
the distance analysis was run (see the git history of `scripts/step6_rules.md`).
The step 3 result (376.6077 km2 at MNDWI 0.127) remains the primary result.
The threshold sweep and elevation results (step 5) are not changed.

This is a descriptive spatial analysis: how the detected 2014 flood footprint
is distributed with respect to the mapped Jhelum centreline. It does not ask
whether distance caused flooding, whether the river connected every flooded
pixel, how depth varies with distance, or how water travelled.

**Exploration done before writing these rules** (scratch files only, not
committed, no flood information used): the river ways in the existing OSM
extract were listed, the OSM waterway relation for the Jhelum was found, and
the distance distribution of the analysis area (not of the flood) was
measured to choose the bands in D2.

### J. The Jhelum centreline (`scripts/23_jhelum.py`)

- J1. Source: the existing Geofabrik extract
  `data/osm/geofabrik/northern-zone-latest.osm.pbf` (OSM timestamp and
  checksums in `data/osm/geofabrik/extract.json`). No river data are
  downloaded.
- J2. Selection: the ways that are members, with role `main_stream`, of OSM
  relation 8306691 (`type=waterway`, `waterway=river`, `name:en=Jhelum River`,
  `wikidata=Q214874`). Members with role `side_stream` are not used. Ways
  named Jhelum that are not members of the relation are not used and are
  listed (when these rules were written: 56241496, 1.6 km, touching the
  main stem at one node near Khanabal; 102436463 and 102436474, "Jhelum
  (Veth)", 14.3 km, at least 888 m from the main stem). No unnamed way is
  added.
- J3. Geometry: all `main_stream` members present in the extract, at full
  length (not clipped to the clip box), so that pixels near the box edge get
  their true nearest distance. Members absent from the extract are listed.
- J4. Continuity: inside the clip box (lon 74.2 to 75.4, lat 33.5 to 34.6)
  the selected ways must form one connected chain through shared OSM nodes.
  Components before and after joining are reported (joining = union of the
  ways; no gap is bridged). If more than one component lies inside the clip
  box, the work stops; gaps are not filled by guessing.
- J5. Missing members: for every analysis-area pixel, the distance to the
  available centreline must not exceed the pixel's distance to the boundary
  of the Geofabrik extract (from the Geofabrik index saved in step 4). If it
  does, a missing member could be closer and the work stops.
- J6. Tributary check: every river way in the extract that touches the
  selected chain but is not selected is listed with its name, to show that
  no tributary is included.
- J7. Persistent-water check: the share of the centreline length inside the
  analysis area that lies within 30, 60 and 90 m of a 2014 persistent-water
  pixel (step 3 flood class 2), and within 90 m of a pixel that is water on
  either date. If less than 50 % lies within 90 m of persistent water, the
  geometry is taken to be misplaced and the work stops.
- J8. Wular Lake: the mapped centreline through Wular is kept as it is. The
  lake surface is not treated as the river, and distances are not measured
  to the lake shore. Pixels inside Wular get their distance to the
  centreline although they lie in a broad water body. This is a stated
  limitation.
- J9. Output for audit: `outputs/jhelum_centreline.geojson` (selected
  geometry, WGS84, with way ids and names) and `outputs/jhelum.json` (ids,
  names, lengths, components, continuity, missing members, tributary list,
  persistent-water check).

### D. Distance to the centreline (`scripts/24_river_distance.py`)

- D1. Distance definition: Euclidean distance in metres in the EPSG:32643
  plane from each pixel centre of the analysis grid to the nearest point of
  the selected centreline, computed with an exact Euclidean distance
  transform (scipy) on the 30 m grid, padded by 20 km on every side, after
  rasterising the centreline onto that grid (all touched pixels). The
  rasterisation limits the error to at most half a pixel diagonal (21.2 m).
  This is not hydrological distance.
- D2. Bands (lower bound inclusive, upper exclusive): 0-250, 250-500,
  500-1,000, 1,000-2,000, 2,000-3,000, 3,000-5,000, 5,000-7,500,
  7,500-10,000, 10,000-15,000 m and 15,000 m or more. Reason: every band
  holds at least 100 km2 of the analysis area (from the distance
  distribution of the area, not of the flood); narrow bands near the river,
  where the floodplain is, and wider ones further out; the narrowest band is
  more than 8 pixels wide.
- D3. Pixels as in step 5 (E2): analysis-area pixels; flooded = class 1;
  assessed = classes 0 to 3; not assessed = class 255. Per band: assessed
  km2, flooded km2, flood fraction (% of assessed), share of the whole
  footprint (%), not-assessed km2, and persistent water (class 2) km2,
  reported separately as context. Persistent water is not removed from or
  added to the flood result.
- D4. Statistics of distance for flooded pixels and for all assessed pixels:
  count, minimum, maximum, mean, median, standard deviation, and the 1st,
  5th, 10th, 25th, 75th, 90th, 95th and 99th percentiles.
- D5. Output raster `data/river/jhelum_distance_m.tif`: float32, same CRS,
  transform and shape as `data/flood/flood_classes.tif`; distance in metres
  for analysis-area pixels; -9999 (nodata) outside the analysis area. Tags
  record CRS, units, source geometry, OSM snapshot and the distance
  definition. The work stops if the grid differs from the flood raster.
- D6. Checks before results are accepted: the centreline lies in the
  intended area (J4); dimensions, CRS and transform match the flood raster;
  outside-analysis pixels are nodata; distances are non-negative and no
  larger than the diagonal of the padded grid.
- D7. Figures: `outputs/figures/jhelum_centreline.png` (selected centreline
  over the analysis area and flood classes) and
  `outputs/figures/river_distance.png` (flood fraction and areas by band).
- D8. Descriptive only. No regression is fitted. No causal language is used.
  Distance and elevation are spatially related.

**Amendment 1 to J5 (committed before any distance-band result was produced)**

- Reason: the first run of `24_river_distance.py` stopped at J5, before any
  band result was computed. 90,526 analysis-area pixels (81.5 km2, all in
  the north-west corner, lon 74.19 to 74.35, lat 34.43 to 34.60) are closer
  to the edge of the Geofabrik extract (19.1 to 30.2 km) than to the mapped
  Jhelum centreline (26.2 to 37.6 km). Nine downstream main_stream members
  are outside the extract (the mapped line leaves it at 34.108 N, 73.900 E),
  so for these pixels a missing member could in principle be closer. The
  diagnosis used distances only, no flood information. No data were
  downloaded.
- Amended J5: a pixel whose distance to the extract boundary is smaller than
  its distance to the mapped centreline is "boundary-affected". Its true
  distance to the full river lies between those two values (lower bound =
  distance to the extract boundary, upper bound = distance to the mapped
  centreline). Such a pixel may proceed only if that whole range lies inside
  one D2 band; it is then assigned to that band for band-level area
  calculations. If the range of any boundary-affected pixel crosses a D2
  band edge, the work stops.
- For this run all 81.5 km2 have both bounds at or above 15,000 m, so they
  stay in the ">= 15,000 m" band. The D2 bands are unchanged.
- Exact distances of boundary-affected pixels are not treated as known. The
  exact-distance statistics (D4) are computed only from pixels that are not
  boundary-affected and are labelled so; the boundary-affected pixels are
  reported separately with their lower and upper bounds. A flag raster
  `data/river/jhelum_distance_flag.tif` (1 = value in the distance raster is
  an upper bound) accompanies the distance raster.
