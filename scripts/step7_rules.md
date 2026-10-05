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
