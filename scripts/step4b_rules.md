## Step 4b rules: building count from Microsoft footprints (fixed before any count)

These rules were written and committed before any footprint was downloaded
or counted (see the git history of `scripts/step4b_rules.md`). They are
copied into METHODS.md unchanged. The step 4 OpenStreetMap building count and
its completeness warning stay in METHODS.md as they are.

**Rules set by the project owner**

1. Source: Microsoft Global ML Building Footprints from the Planetary
   Computer catalogue (collection `ms-buildings`), for the analysis area. If
   the collection is not there or does not cover the area, the work stops.
   No other source is used.
2. The dataset version and any imagery dates it provides are recorded.
3. Counting rules are the same as step 4: centroid rule as the main figure,
   any-overlap as a sensitivity, not-assessed pixels reported separately.
4. The result is described as buildings in today's footprint data standing
   on ground that was flooded on 2014-09-10, never as buildings flooded in
   2014.
5. The OpenStreetMap building count stays in METHODS.md as it is.

**How the rules are implemented (fixed at the same time)**

- M1. Items: `ms-buildings` items intersecting the clip box with the newest
  `msbuildings:processing-date` (the dataset version). Only that version is
  used; older, unpartitioned items are not mixed in. The work stops if the
  union of the chosen items' geometries does not contain the whole analysis
  area, or if any partition cannot be read completely.
- M2. Recorded: processing date, item ids, region, quadkeys,
  `start_datetime` / `end_datetime` (the imagery date range the catalogue
  gives), row counts, and for every parquet file its name, size, the
  server's Content-MD5 where given, and a locally computed sha256.
- M3. Each row's WKB geometry is one footprint. Footprints with byte-identical
  geometry in more than one file are counted once; the number is reported.
  Footprints that cannot be parsed or are empty are reported, not counted.
- M4. Counting: rules B3 (centroid pixel decides flooded / not flooded / not
  assessed / outside), B4 (any-overlap sensitivity, all-touched
  rasterisation), S1 (Srinagar box window) and W1 (step 3 water/wetland
  mask) of step 4, unchanged.
- M5. 1 km grid over the Srinagar box: UTM cells aligned to multiples of
  1,000 m, clipped to the Srinagar box window. Per cell: footprints by
  centroid (flooded, total), and motor-road length (step 4 rule R2, OSM
  extract) flooded and total within the analysis area; share = flooded /
  total motor-road length. Labels: nearest named OSM suburb, neighbourhood
  or quarter (the step 4 place points), with the distance. Three top-ten
  lists: by flooded footprints, by flooded motor-road km, and by flooded
  share among cells with at least 1 km of motor road.
- M6. Cross-check: the ten 1 km cells with the most OSM buildings (step 4,
  centroid in the analysis area). For each OSM building polygon there (nodes
  are skipped and counted), a match is a Microsoft footprint whose
  intersection with it is at least 50 % of the OSM polygon's area or at least
  50 % of the footprint's area. Also reported: OSM buildings touched by any
  footprint.
- M7. Completeness: step 4 rules C1 and C2 repeated with the Microsoft
  footprints (2 km cells, at least 2 km of motor road, footprints per km of
  motor road).
- M8. Wording: every reported figure is "footprints in the Microsoft data
  (imagery 2017-2021 per the catalogue) standing on ground flooded on
  2014-09-10".
