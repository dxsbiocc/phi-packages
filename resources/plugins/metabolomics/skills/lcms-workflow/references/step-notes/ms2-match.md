# mb_ms2_match — spectral library matching

## Principle

MS2 fragmentation patterns are the strongest routinely available evidence
for small-molecule identity. Matching query spectra against a reference
library by cosine similarity is the standard putative-identification step
(MSnbase/OrgMassSpecR/msPurity practice; MassBank, MoNA, GNPS ecosystems).

## Parsing

- **MSP**: key/value blocks separated by blank lines. Tolerated dialects:
  NIST (`Name:`, `PrecursorMZ:`, `Num Peaks:`), MoNA/RIKEN (`NAME:`,
  `PRECURSORMZ:`), `=` separators, precursor buried in `Comments:`.
  Peak lines: `mz intensity` (whitespace/colon/semicolon separated).
- **MGF**: `BEGIN IONS` … `TITLE=`, `PEPMASS=` (first token), peak lines,
  `END IONS`. Query must be MGF; library may be MSP or MGF (content-sniffed).
- Entries without a precursor m/z cannot pass the precursor filter and are
  counted in `n_library_without_precursor`.

## Scoring

- Precursor gate: |Δ| ≤ ppm × query precursor / 1e6 (default 10 ppm).
- Fragments: sqrt-intensity transform, greedy nearest-pair matching within
  `fragment_tol` Da (default 0.02), each peak used once; cosine = matched
  dot product / (full-vector norms). Unmatched peaks lower the score —
  chimeric spectra are penalized, not cherry-picked.
- Report floor: `min_score` 0.7, `top` 5 hits per query. Common practice
  also asks for ≥3–5 matched fragments; check `matched_fragments`.

## Interpretation

- A hit is **MSI level 2** (putative library match). Level 1 needs an
  authentic standard on the same method (RT + MS2). Keep the level in prose.
- Library provenance is part of the result: name, version, license
  (MassBank CC-BY, MoNA CC-BY-SA, GNPS community libraries, commercial NIST
  needs a license). Confirm with the user before use; cite in the report.
- Isomers share fragments; a high cosine does not distinguish them. RT
  evidence (not implemented here) or standards resolve that.

## Pitfalls

- No MS2 spectra → no matching. Do not infer identity from precursor m/z.
- Different collision energies distort cosine; note the library's energy
  when known. GNPS networking and in-silico tools (SIRIUS/MetFrag) are gaps.