# Interpretation note: validation metrics of commit `4e68f07` (diagnostic run)

**These metrics must not be presented as the accuracy of the flood map.**

`outputs/validation/metrics_v1.json` and `outputs/validation/metrics_v1_summary.csv`
were computed under rule M (with amendments 2 and 3 of `scripts/step7_rules.md`)
from the 116 frozen reference rows in `outputs/validation/labels_v1_final.csv`.
They are kept unchanged as a recorded diagnostic validation run.

Why they are not the map's accuracy: the completed reference labels show a
systematic disagreement with the intended operational meaning of "water",
particularly for the 25 August 2014 image. Taken at face value they imply
that about 43 % of the assessed analysis area (about 1,600 km2) was water on
25 August, against 37 km2 in the map and 26 km2 in the USGS water flag, in
an image that is mostly green fields and slopes; in strata the map shows as
dry, roughly 4 in 10 samples were labelled water on that date, and in the
"dry, other" stratum the points labelled water were labelled water on both
dates. This points to a different reading of "water at the centre pixel"
(for example dark or waterlogged vegetated ground) rather than to water
missed by the map. The cause was not investigated point by point; no
individual disagreement point was inspected.

What follows from this: the first validation attempt produced anomalous
human labels. It motivated a prospective, fresh blind validation on the
untouched 140-point W sample, with a clarified definition of "water at the
centre pixel" recorded in amendment 4 before any W label was collected.
The 116 labels were used only to identify this general ambiguity; they were
not used to tune the new protocol in any other way, and they are not
relabelled, modified or reinterpreted.

Unchanged: the frozen 116 labels (sha256 8cf11c56...dc2c), the metrics files
of `4e68f07`, amendment 3 (`e512769`), and the primary flood result of
376.6077 km2.
