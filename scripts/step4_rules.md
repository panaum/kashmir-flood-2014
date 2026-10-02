## Step 4 rules: buildings and roads (fixed before any count was computed)

These rules were written and committed before any building or road count was
computed (see the git history of `scripts/step4_rules.md`). They are copied
into METHODS.md unchanged.

**Rules set by the project owner**

1. Source: OpenStreetMap buildings (any building tag) and roads (highway
   tag) for the analysis area, fetched now. Record the download time, save
   the raw responses with checksums, and split the request into tiles if
   needed.
2. A building counts as flooded if its centroid falls in a flooded pixel of
   the step 3 flood mask. Also report the count using any overlap with a
   flooded pixel, as a sensitivity.
3. Road length is the length of road inside flooded pixels, in kilometres.
   Report motor roads (motorway down to residential, plus service,
   unclassified and living_street) as the main figure, and all other highway
   types separately.
4. Buildings and roads in pixels that were not assessed are reported as not
   assessed, not as dry.
5. The map is today's OpenStreetMap, not 2014's. This limitation is stated
   and not corrected for.

**How the rules are implemented (fixed at the same time)**

- B1. Download: Overpass API (same endpoint and client name as step 3),
  over the clip box (lon 74.2-75.4, lat 33.5-34.6) in tiles of 0.25 x 0.25
  degrees. A tile that fails is split into four, down to 0.0625 degrees; if
  a piece still fails after 3 attempts, the work stops. Each raw response is
  saved with its sha256 and download time. Elements appearing in more than
  one tile are kept once (by OSM type and id).
- B2. Buildings = nodes, ways and relations with a `building` key of any
  value (counts by value are reported, including `building=no` if present).
  Ways must be closed to form a polygon; relations are built from their
  outer/inner members as in step 3. Nodes are points.
- B3. Geometry is projected to the mosaic CRS (EPSG:32643). Centroid =
  polygon centroid (point itself for nodes). The centroid's pixel decides:
  flooded (flood class 1), assessed not flooded (0, 2, 3), not assessed
  (analysis-area pixel with class 255), or outside the analysis area.
  Totals are buildings whose centroid is in the analysis area.
- B4. Overlap sensitivity: the pixels a building touches (all-touched
  rasterisation on the 30 m grid; the centroid pixel for nodes). Flooded if
  any touched pixel is flooded; otherwise not assessed if any touched pixel
  is not assessed; otherwise dry. Total = buildings touching the analysis
  area.
- R1. Roads = ways with a `highway` key of any value. Ways tagged area=yes
  are areas, not lines: they are excluded from lengths and counted
  separately.
- R2. Motor roads (main figure): highway = motorway, motorway_link, trunk,
  trunk_link, primary, primary_link, secondary, secondary_link, tertiary,
  tertiary_link, residential, unclassified, service, living_street. Every
  other highway value is reported separately, by value.
- R3. Lengths: lines are clipped to the exact pixel outlines of the flooded,
  not-assessed and analysis-area pixels, in EPSG:32643 map units, reported
  in km. No correction for projection scale (under 0.05 % here) or slope.
- S1. Srinagar box: the step 3 pixel window (lon 74.65-74.95, lat
  33.95-34.20); buildings by centroid pixel, roads clipped to the window
  outline.
- W1. Inside / outside water and wetland: the step 3 OSM pixel mask
  (`data/osm/water_wetland_mask.tif`); buildings by centroid pixel, roads
  clipped to the mask's pixel outlines.

**Areas within Srinagar**

- A1. Query the Srinagar box for boundary=administrative relations (all
  admin_level values) and for place=suburb, neighbourhood or quarter (nodes,
  ways, relations).
- A2. Candidate polygon sets: each admin_level separately, and the set of
  place polygons. Coverage of a set = share of the analysis-area buildings
  (by centroid) in the Srinagar box that fall inside one of its polygons.
  Only sets with at least 5 polygons in the box count as areas within the
  city.
- A3. If a candidate set has coverage >= 90 %, the finest one (highest
  admin_level; place polygons preferred over admin if both qualify) is used.
  Otherwise a 1 km grid is used: UTM cells aligned to multiples of 1,000 m,
  each labelled with the nearest named suburb, neighbourhood or quarter
  (node position, or polygon centroid), with the distance recorded.
- A4. The ten areas or cells with the most flooded buildings (centroid rule)
  are reported, with flooded buildings, total buildings and flooded motor
  road length.

**Completeness warning**

- C1. On a 2 km UTM grid over the analysis area, cells are flagged "roads
  but no buildings" if they hold >= 2 km of motor road (R2) and no building
  centroid. Reported: number of such cells, their area, and the flooded area
  inside them.
- C2. Buildings per km of motor road, per 2 km cell, summarised as
  percentiles over cells with >= 2 km of motor road, and the same for the
  Srinagar box against the rest of the analysis area.
- C3. Nothing here corrects counts; it only warns where OSM may be missing
  buildings.
