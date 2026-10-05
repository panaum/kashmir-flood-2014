## Step 7 rules: blinded reference-label validation (fixed before any sample was drawn)

These rules were written and committed before the validation sample was drawn
(see the git history of `scripts/step7_rules.md`). The step 3 result
(376.6077 km2 at MNDWI 0.127) remains the primary result. Nothing here
replaces it.

**Nature of the reference.** The reference labels are image-based reference
interpretation of the same Landsat 8 scenes used for the map, made by a human
interpreter who knows the project. They are not field truth and not
independent ground truth, and they cannot resolve features smaller than the
30 m pixels. Every result of this step is reported with that limitation.

### S. Strata (mutually exclusive; a pixel takes the first stratum it matches)

Population: assessed analysis-area pixels (step 3 flood classes 0 to 3 inside
`data/flood/analysis_area.tif`). Not-assessed pixels are never sampled.

- S1 Flooded, built-up: class 1 and built-up proxy.
- S2 Flooded, other: class 1.
- S3 Water in August: classes 2 and 3.
- S4 Dry, inside mapped lakes or wetlands: class 0 and inside the step 3 OSM
  mask `data/osm/water_wetland_mask.tif`.
- S5 Dry, built-up: class 0 and built-up proxy.
- S6 Dry, near flooded ground: class 0 and within 500 m of a flooded (class 1)
  pixel, by exact Euclidean distance between pixel centres.
- S7 Dry, other: class 0.

**Built-up proxy.** A pixel is built-up if the 5 x 5 pixel window (150 m x
150 m) centred on it contains at least 10 Microsoft building-footprint
centroids (`data/msbuildings/footprints.csv`, step 4b; version 2023-04-25,
imagery 2017-2021). This is a present-day layer. **It is used only to decide
where samples are drawn. It is never used as a reference label and never
enters the reference data.**

### P. Sampling

- P1. 56 samples per stratum, 392 in total.
- P2. Random generator: `numpy.random.default_rng(20140910)`. Strata are
  processed in the order S1 to S7. Within each stratum the candidate pixels
  are put in random order with that generator and accepted one by one if the
  pixel centre is at least 300 m from every sample already accepted in any
  stratum. If a stratum cannot reach 56 samples, the work stops.
- P3. After all strata are filled, the 392 samples are put in random order
  with the same generator and numbered V001 to V392 in that order.
- P4. The answer key (sample id, stratum, map class, row, column) is written
  to `data/validation/key.csv`, which is not committed and not shown to the
  interpreter. Its sha256 is committed with the labelling material so the key
  cannot change after labelling.

### B. Blinding of the labelling material

- B1. For each sample the interpreter sees only: the sample id, its
  coordinates, and image chips made from the reflectance bands of the two
  Landsat dates.
- B2. Chips: a detail chip of 51 x 51 pixels (1.53 km) with the centre pixel
  outlined, in two band combinations (SWIR1/NIR/Red and NIR/Red/Green) for
  each date; and a context chip of 201 x 201 pixels (6.03 km) in
  SWIR1/NIR/Red for each date. The scenes have no blue band, so there is no
  true-colour chip. Colour stretch: the step 1 stretch for SWIR1, NIR and Red
  (DATA.md); for Green, the same rule (2nd to 98th percentile of valid valley
  pixels, both dates pooled). The same stretch is used on both dates.
- B3. Never shown, in any file the interpreter opens: map class, stratum,
  MNDWI, water masks, flood classes, OSM layers, building footprints or any
  other classification output. A script check scans the labelling files
  (data and visible text) for such terms and stops if any is found.
- B4. The labelling page runs locally from files; it loads nothing from the
  internet.

### L. Labels the interpreter records

For each sample and each date: water at the centre pixel (`yes`, `no`,
`cannot tell`) and the surface seen (open water; water with vegetation;
built-up; crops or grass; trees; bare soil; wetland or marsh; cloud or
shadow); one overall confidence (high, medium, low); `other_source` (free
text, if other imagery was consulted, with its date); notes.

### C. Cannot-tell responses

- C1. Every `cannot tell` is recorded as such. It is never converted to water
  or not water.
- C2. Counts of `cannot tell` are reported overall, by date and by stratum.
- C3. Primary confusion matrices and accuracy metrics exclude `cannot tell`
  observations for the date concerned; the number excluded is reported next
  to every metric. A flood observation is `cannot tell` if either date is.
  Exclusion assumes that `cannot tell` samples are missing at random within
  their stratum; this assumption is stated.
- C4. If `cannot tell` exceeds 5 % of all samples for a date, or 10 % of any
  stratum, this is flagged, and a sensitivity analysis is reported: all
  metrics and the reference-based area recomputed with every `cannot tell`
  counted as water, and again with every one counted as not water.

### M. Metrics (computed only after the labels are supplied)

Map labels per sample, from the hidden key: water on 25 August = classes 2
or 3; water on 10 September = classes 1 or 2; flooded = class 1. Reference
labels: the interpreter's answers; reference flooded = water on 10 September
and not water on 25 August.

- M1. Because the strata are not the map classes, all estimates use the
  stratified estimator for strata that differ from map classes (Stehman,
  2014). With stratum weights W_h = (stratum area) / (total assessed area),
  n_h usable samples in stratum h, and an indicator y for each sample:
  estimated proportion Y = sum_h W_h * mean_h(y), with variance
  V(Y) = sum_h W_h^2 * s_h^2(y) / n_h.
- M2. Ratios (user's accuracy = precision, producer's accuracy = recall) use
  R = Y / X with indicators y and x, and the linearised variance
  V(R) = (1 / X^2) * sum_h W_h^2 * (s_h^2(y) + R^2 s_h^2(x) - 2 R s_h(x, y)) / n_h.
- M3. Reported per date (water) and for the flood: the estimated error
  matrix (proportions of area), overall accuracy, precision, recall and F1
  (2 P R / (P + R)) for the water / flooded class, each with a 95 %
  confidence interval (plus or minus 1.96 standard errors); also the raw,
  unweighted sample counts by stratum.
- M4. Reference-based flooded area (secondary): A = (total assessed area) x
  Y(reference flooded), with a 95 % confidence interval. It is reported next
  to the primary 376.6077 km2 result and never replaces it.
- M5. Per-stratum results (raw agreement counts) are reported, in particular
  for S1, S4 and S5.
- M6. Low-confidence labels are kept in the primary metrics; the metrics are
  also reported without them as a sensitivity.

**Amendment 1: sample design v2 (committed before the v2 sample was drawn)**

- Change: the validation sample is reduced from 56 samples per stratum (392
  in total, design v1) to **20 samples per stratum (140 in total, design
  v2)**.
- Reason: the v1 sample (drawn in commit `0dc1aea` under the rules of commit
  `cd73ad0`) has not been labelled and has not been used for any metric. A
  140-point sample still represents all seven strata. It is treated as a
  smaller validation and diagnostic sample. **It does not give the same
  statistical precision as 392 samples**: with 20 samples per stratum,
  standard errors and confidence intervals will be wider, and per-stratum
  results will rest on 20 observations each.
- The v1 sample and its material stay in the git history and on disk
  unchanged. They will not be labelled or used for any validation result, so
  that results cannot be chosen between the two designs.
- Unchanged: the seven strata S1 to S7 and their order and definitions,
  including the present-day built-up proxy (used only for stratification,
  never as a reference label); random order within each stratum; the 300 m
  minimum spacing; never sampling not-assessed pixels; the hidden answer key
  with a committed checksum; the blinding rules B1 to B4; the label fields
  (L); cannot-tell handling (C); and the metrics and estimators (M), applied
  with n_h = 20.
- v2 random generator: `numpy.random.default_rng(20140825)`, a new seed; the
  v1 sample and seed (20140910) are not reused. Samples are numbered W001 to
  W140, so they cannot be confused with v1 ids (V001 to V392).
- v2 paths: hidden key `data/validation/v2/key.csv` (not committed);
  labelling material `data/validation/labelling_v2/` (not committed); blind
  sample list, labels template and key checksum in `outputs/validation/v2/`.

**Amendment 2: reference labels collected on the v1 sample (committed before any validation result)**

- What happened: the human labelling was carried out on the original v1
  sample (392 points, ids V001 to V392, rules of commit `cd73ad0`, sample of
  commit `0dc1aea`), not on the v2 sample of amendment 1. This was a mistake:
  the two labelling pages look alike.
- It was discovered before any validation metric, agreement statistic,
  confusion matrix, accuracy estimate or other validation result was
  calculated. Only integrity checks had been run (ids, duplicates, allowed
  values, completeness, counts by stratum).
- The labelling stayed blind: the v1 page and chips show no map class, MNDWI,
  stratum, OSM layer, building outline or validation outcome (checked by the
  B3 blinding check in commit `0dc1aea`), and none was revealed during
  labelling.
- The labels are therefore kept as an independently collected human
  reference set, and amendment 1 is superseded for the choice of sample: the
  validation uses the v1 sample. The v2 sample (commit `fb17be8`) is not
  labelled or used for any result.
- Nothing is rewritten or deleted: the original rules, the v1 sample and its
  key (`data/validation/key.csv`, sha256 1be0d65a...a624 as committed), and
  amendment 1 and its v2 sample stay as they are. The 376.6077 km2 result and
  all earlier results are unchanged.
- Reference set: every v1 row in which August water, September water and
  confidence are all present. Inclusion depends only on whether a row was
  completed, never on its answers or on the map class. No sample is
  rebalanced, replaced or removed on the basis of its labels or its map
  class. The v1 ids were assigned in random order (rule P3) and were
  labelled largely in id order, so the completed rows are a subset of the
  stratified sample with unequal numbers per stratum; the stratified
  estimators of rule M are used with the actual n_h of each stratum.
- The 11 rows that have both water answers but no confidence are completed
  by the interpreter (confidence only) before the reference set is frozen.
  The interim export received on 2026-10-05 is kept, not committed, at
  `data/validation/received/labels_v1_interim_2026-10-05.csv` (sha256
  1408e31efd0e8c1ac0b98d7a3294fa3d252087f7d6fffe2772dd5a600799fbba).
  The final labels file is committed, with its sha256, before any metric is
  computed.
- Surface fields: the surface (cover) fields of rule L were not meaningfully
  collected in this session (blank in almost every completed row). This is
  recorded as a gap. They are not filled in or inferred, and no result that
  depends on them is reported. The water validation uses the August water,
  September water and confidence fields only.

**Amendment 3: clarifications of rule M (committed before any metric was computed)**

Rule M left three points open. They were settled by the project owner before
any validation metric was calculated:

- F1 interval: F1 is computed as the ratio of stratified means F1 = Y / X
  with y = 2 x [map water and reference water] and x = [map water] +
  [reference water] (for the flood: flooded instead of water). This equals
  2 P R / (P + R), and its variance is the M2 linearised ratio variance,
  applied unchanged (the formula holds for any y and x, not only 0/1
  indicators).
- Variance denominator: s_h^2 and s_h(x, y) are the sample variance and
  covariance within stratum h with denominator n_h - 1.
- Interval bounds: 95 % intervals are reported exactly as estimate plus or
  minus 1.96 standard errors, not clipped; every bound outside [0, 1] (or
  below 0 km2 for the area) is flagged.

Implementation details that follow directly from rules S and M: stratum
weights W_h use the stratum areas from rule S on the step 3 rasters (the
same computation that drew the sample), and every sample's stratum in the
hidden key is checked against them; n_h is the number of reference rows in
stratum h usable for the quantity concerned (after cannot-tell exclusion,
rule C); no finite-population correction is applied, as M1 states none.
