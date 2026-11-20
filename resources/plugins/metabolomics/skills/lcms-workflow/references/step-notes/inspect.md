# mb_inspect — input recognition

## Principle

Deterministic header/content sniffing, no guessing. The tool classifies each
file and reports what follow-up steps it supports; the agent decides.

## What it recognizes

- **Feature tables** by vendor headers: MS-DIAL (`Alignment ID`,
  `Average Mz`), MZmine (`row ID`, `m/z`, `retention time`), xcms
  (`mz`/`rt` or `mzmed`/`rtmed`), MetaboAnalyst (second column
  `label`/`class`/`group` with numeric feature columns), else generic.
- **mz/rt aliases** are normalized to `mz`/`rt` metadata columns; they are
  never treated as samples. Non-numeric annotation columns stay metadata.
- **Metadata**: small tables with non-numeric columns and no intensity
  matrix. **MSP/MGF** libraries by content (`BEGIN IONS`, `Num Peaks`).
  **Mappings** (two columns, repeating first column), **compound DBs**
  (id + mass), single-column **m/z / ID lists**.

## Interpretation

- `rt_unit`: median RT > 60 → seconds (xcms convention), else minutes.
  This is an inference — confirm against the vendor export.
- `missing_rate` counts NaN; `zero_rate` counts zeros separately. Both are
  treated as missing by `mb_preprocess`; a high zero rate with low NaN rate
  usually means the exporter already imputed zeros.
- `qc_samples` matches names containing `qc`/`pool` (case-insensitive),
  excluding `blank`. Rename or map samples before QC if the vendor used
  other labels.
- MetaboAnalyst orientation is samples-as-rows: transpose first. The tool
  reports, it does not transform.

## Pitfalls

- A directory scan includes every non-hidden file; previous JSON outputs
  show up as `unknown`. Pass clean input directories.
- Feature IDs are taken from the first column as strings; duplicate IDs
  will collide downstream — check the export.