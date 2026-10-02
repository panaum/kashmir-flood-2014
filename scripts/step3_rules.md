## Step 3 rules: flood extent (fixed before any flood area was computed)

These rules were written and committed before any flood area was computed
(see the git history of `scripts/step3_rules.md`). They are copied into
METHODS.md unchanged. The result is not to be adjusted toward any published
figure.

**Rules set by the project owner**

1. Water threshold: the 2014-09-10 Otsu value from step 2 (0.127), applied
   to both dates.
2. Analysis area: valley floor (DEM <= 1,800 m) connected to Srinagar only.
3. QA_RADSAT is downloaded for all four scenes with the usual checks, and
   pixels saturated in green or SWIR1 are masked.
4. Pixels with SWIR1 reflectance <= 0 and green > 0 count as water. Pixels
   with both <= 0 stay masked.
5. Cirrus stays unmasked.
6. Flooded = water on 2014-09-10 and not water on 2014-08-25, among pixels
   valid on both dates.
7. The result is not to be adjusted toward any published figure.

**How the rules are implemented (fixed at the same time)**

- R1. The threshold is the value recorded in `data/water_report.json` for
  2014-09-10 (bin edge 0.127). Water = MNDWI > threshold (strictly greater),
  as in step 2.
- R2. Analysis area = DEM <= 1,800 m AND inside the 8-connected patch of such
  pixels that contains 34.08 N, 74.80 E (the patch used in step 2).
- R3. A pixel is saturated if QA_RADSAT bit 2 (band 3, green) or bit 5
  (band 6, SWIR1) is set, per the STAC `classification:bitfields`. QA_RADSAT
  is taken from the same row (036 or 037) that the mosaic used for that pixel
  (`*_source.tif`). Terrain occlusion (bit 11) is reported, not masked.
- R4. Pixels with green <= 0 and SWIR1 > 0 are not covered by rule 4; they
  stay masked, as in step 2 (MNDWI is not meaningful there).
- R5. Valid on a date = in the analysis area AND QA_PIXEL bits 0 (fill),
  1 (dilated cloud), 3 (cloud), 4 (cloud shadow), 5 (snow) all 0 AND not
  saturated (R3) AND not (green <= 0 and SWIR1 <= 0) AND not (green <= 0 and
  SWIR1 > 0). Water on a date = valid AND (MNDWI > threshold OR (SWIR1 <= 0
  and green > 0)).
- R6. Among pixels valid on both dates: flooded = water Sep and not water
  Aug; persistent water = water on both; water lost = water Aug and not
  water Sep. Not assessed = analysis area pixels not valid on at least one
  date.
- R7. Sensitivity variants change only the thresholds; masks are identical.
  (a) each date's own step 2 Otsu value as recorded (2014-08-25: -0.058;
  2014-09-10: 0.127), not recomputed after the new saturation mask;
  (b) 0 on both dates; (c) 0.077 and 0.177 on both dates.
- R8. Area = pixel count x 0.0009 km2 (30 m x 30 m UTM map units).
- R9. Srinagar box: lon 74.65 to 74.95, lat 33.95 to 34.20, the same pixel
  window as the step 1 previews.

**Unflagged clouds on 2014-09-10 (estimate only, masks not changed)**

Computed over analysis-area pixels valid on 2014-09-10. Thresholds are fixed
here, before looking at the result:

- U1. QA_PIXEL's own lower-confidence flag: cloud confidence (bits 8-9) =
  medium (2) while the cloud bit is 0.
- U2. Bright-and-white spectral test with the bands already downloaded:
  green, red, NIR and SWIR1 reflectance all > 0.20, AND
  (max - min) / mean of green, red, NIR < 0.7.
- U3. Shadow zone: pixels reached by moving each U2 pixel away from the sun
  (sun azimuth and elevation from the 2014-09-10 row 036 MTL) for cloud
  heights 500 m to 5,000 m in 250 m steps. Reported as the whole zone (upper
  bound) and the dark part of it (NIR < 0.15 and SWIR1 < 0.10).
- U4. Control: U2 is also applied to 2014-08-25 analysis-area pixels that
  QA_PIXEL flags clear (bit 6); the rate found there indicates how often the
  test fires on clear ground.
- U5. Reported: area of U1, U2, U3 (zone and dark part), and the flooded area
  falling inside each.

**OpenStreetMap water bodies and wetlands**

- O1. Overpass API query for the clip box: ways and relations tagged
  natural=water, natural=wetland, landuse=reservoir or waterway=riverbank.
  The raw response and the download time (UTC) are saved.
- O2. Closed ways and multipolygon relations are turned into polygons; a
  pixel is inside if its centre is inside any polygon. The OSM data are
  current, not from 2014.
- O3. Flooded area is split into inside and outside these polygons. If the
  download fails, the work stops.

**Independent check**

- I1. For each date, the water mask is cross-tabulated against QA_PIXEL bit 7
  (water) on valid analysis-area pixels: agreement, both disagreement
  directions, and Cohen's kappa. Bit 7 is produced by USGS from the same
  scene, so it is a check against a different algorithm, not against
  different data.
