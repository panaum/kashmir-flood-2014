## Step 5 rules: threshold sweep and elevation (fixed before any result was computed)

These rules were written and committed before the threshold sweep or the
elevation analysis was run (see the git history of `scripts/step5_rules.md`).
They cover only those two analyses. The step 3 result (376.6077 km2 at MNDWI
0.127) remains the primary result and is not replaced.

**Baseline confirmed before writing these rules** (from `data/flood_report.json`):
flooded 376.6077 km2 at threshold 0.127 on both dates; analysis area
3696.7941 km2; valid on both dates 3696.2415 km2.

### A. Continuous threshold sensitivity (`scripts/21_threshold_sweep.py`)

- T1. Implementation: the step 3 functions in `scripts/10_flood.py`
  (`load_date`, `water`, `outcome`) are imported and used unchanged. No second
  MNDWI or flood implementation is written. Inputs, mosaic, CRS, pixel grid,
  QA and saturation masking, analysis area (`data/flood/analysis_area.tif`,
  checked against the step 3 area), the SWIR1 <= 0 rule, the "valid on both
  dates" rule and the flood definition (water on 2014-09-10 and not water on
  2014-08-25) are exactly those of step 3.
- T2. Grid: thresholds t = k / 100 for k = 0, 1, ..., 30 (0.00 to 0.30
  inclusive, 31 values), the same threshold applied to both dates. Water =
  MNDWI > t, as in step 3.
- T3. The primary threshold is evaluated separately with the exact value
  recorded in `data/water_report.json` (0.127). If it does not reproduce
  376.6077 km2, the work stops.
- T4. Negative thresholds are not part of the grid. The earlier variant with
  -0.058 for August used different thresholds on the two dates and stays
  reported in METHODS.md as it is.
- T5. Reported for every threshold: flooded pixel count, flooded area (km2),
  flooded share of the analysis area (%), and for context persistent water,
  water lost, August water and September water (km2), all on pixels valid on
  both dates.
- T6. Stability is described, not optimised: the minimum and maximum flooded
  area over the grid; the same over the grid points within +/- 0.05 of 0.127
  (0.08 to 0.17); and the change in area per 0.01 step. No threshold is chosen
  or recommended from the sweep.
- T7. Figure: threshold vs flooded area, with 0.127 marked and the range
  0.077 to 0.177 (the step 3 +/- 0.05 variants) shaded.

### B. Flood area and fraction by elevation (`scripts/22_elevation.py`)

- E1. Inputs: `data/clipped/dem_cop30_bilinear.tif` (Copernicus GLO-30, the
  project's one resampled layer, bilinear onto the Landsat grid),
  `data/flood/flood_classes.tif` and `data/flood/analysis_area.tif`. All three
  share EPSG:32643, 30 m pixels, 3716 x 4077. Nothing is resampled.
- E2. Pixels: analysis-area pixels. Flooded = flood class 1. Assessed = flood
  classes 0 to 3 (valid on both dates). Not assessed = analysis-area pixels
  with class 255. DEM nodata (-9999) pixels are excluded and counted; the
  work stops if any analysis-area pixel has DEM nodata.
- E3. Bands (fixed before computing, from the elevation distribution of the
  analysis area only): one lower band "below 1,570 m" (from the lowest
  analysis-area elevation, 1495.2 m, up to 1570 m), then 10 m bands
  [1570, 1580), [1580, 1590), ..., [1780, 1790), and [1790, 1800] (the last
  one closed at 1800 m, the analysis-area limit). 24 bands in all. Reason:
  the valley floor is concentrated between 1570 and 1600 m and needs 10 m
  detail; below 1570 m there are only 4.9 km2 in total, too little for
  stable per-band fractions.
- E4. Per band: total analysis-area pixels, assessed pixels, not-assessed
  pixels, flooded pixels, flooded area (km2), flood fraction (% of assessed
  pixels), persistent-water pixels, and the band's share of all flooded area.
  Bands with less than 10 km2 assessed are flagged "small band".
- E5. Statistics of DEM elevation for flooded pixels, and for all assessed
  analysis-area pixels as the baseline: count, minimum, maximum, mean,
  median, standard deviation, and the 1st, 5th, 10th, 25th, 75th, 90th, 95th
  and 99th percentiles.
- E6. Figure: flood fraction (%) by band, with flooded area (km2) by band,
  and the assessed-area distribution for comparison.
- E7. Interpretation is descriptive only. Elevation is not presented as a
  cause of flooding: elevation and distance to the river are spatially
  related, and the DEM is a surface model.
- E8. Stated with the results: Copernicus GLO-30 is a digital surface model,
  so heights include buildings and vegetation; its vertical datum is not
  stated in the files and is not independently established here; the
  analysis is truncated at 1,800 m by the existing analysis-area definition,
  which was itself made from this DEM.

### Outputs

New results go to `outputs/` (CSV and JSON) and `outputs/figures/` (PNG),
separate from the step 1 to 4b outputs, which are not modified.
