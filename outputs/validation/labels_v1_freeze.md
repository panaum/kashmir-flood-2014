# Frozen reference labels, validation design v1 (step 7, amendment 2)

Frozen on 2026-10-05, before any validation metric was calculated.

| | |
|---|---|
| Frozen file | `outputs/validation/labels_v1_final.csv` (exact byte copy of the interpreter's export `labels (2).csv`, not edited) |
| sha256 | `8cf11c56cb3af6333bace9e4748b8c3499c438196e64577c05ad1752d305dc2c` |
| Interim file compared against | `data/validation/received/labels_v1_interim_2026-10-05.csv` (not committed), sha256 `1408e31efd0e8c1ac0b98d7a3294fa3d252087f7d6fffe2772dd5a600799fbba` |
| Sample | v1, ids V001 to V392 (rules `cd73ad0`, sample `0dc1aea`); hidden key sha256 `1be0d65a9446857f9bc374a9a62b0788eac0bae48ae53e249253dd12c657a624` |
| Rules | `scripts/step7_rules.md`, amendment 2 (`b9787f5`) |

## Reference set

A row is completed when August water, September water and confidence are all
present. **The reference set is exactly 116 completed rows** (V001 to V118
range), with no selection by answer or map class.

| S1 | S2 | S3 | S4 | S5 | S6 | S7 | Total |
|---|---|---|---|---|---|---|---|
| 18 | 9 | 21 | 16 | 19 | 17 | 16 | 116 |

Confidence among completed rows: 111 high, 5 medium, 0 low.
"Cannot tell": 0 for 25 August, 0 for 10 September.

## Changes since the interim file

- The 11 intended rows (V002, V004, V006, V007, V008, V009, V013, V014, V015,
  V024, V068) changed only in their confidence field, each from empty to a
  rating.
- No August-water or September-water answer changed in any row. Ids and
  coordinates are unchanged and match the committed blind list.

## V005 anomaly

- V005 contains an unexpected confidence value ("medium").
- V005 was completely empty in the preserved interim file.
- V005 has no August-water and no September-water answer.
- Therefore V005 is not part of the completed reference set, which remains
  exactly 116 rows.
- The V005 confidence is treated as a stray annotation and is excluded from
  all validation analysis. The file was frozen as received; V005 was not
  cleared, completed or otherwise altered.
- The anomaly was discovered before any validation metric was calculated.

## Gaps

- Surface (cover) fields were not meaningfully collected: given for 0 of the
  completed rows on 25 August and 1 on 10 September. They are not filled or
  inferred, and no result depends on them (amendment 2).
- `other_source` and `notes` are empty in all completed rows.

## Statement

No automated map class or prediction was inspected for this freeze; only
the stratum column of the hidden key was read, to count completed rows by
stratum. No validation metric, agreement statistic, confusion matrix,
accuracy estimate, confidence interval or reference-based area has been
calculated.
