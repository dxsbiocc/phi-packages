# Section patterns

Use the pattern that matches the request. Omit a slot when the source files do not contain it, and mark the omission.

## Methods

One paragraph per step, in run order.

- Material: sample unit, groups, and the inclusion rule, with the file that defines them.
- Processing: tool name and the version recorded by the skill or the environment, plus the parameters that change the result (filters, normalization, model, seed, reference build).
- Statistics: the comparison, the test, and how multiple comparisons were handled, when that was recorded.
- Output: the table or figure path a reader would open.

Do not describe packages that were not used. Do not copy a tutorial from `scanpy` or `scikit-learn` into the methods; describe the run that happened.

## Results

- Start with what was measured and the sample size.
- Then the contrast: groups or conditions, direction, and the estimate as written in the output.
- Then the display: "Figure 1 shows …" only when that file exists. Say what the axes and groups are.
- End with the limit that is visible in the output (missing samples, failed runs, filters that dropped records).

Do not add a mechanism that the table does not measure.

## Figure or table caption

- First sentence: what is plotted or tabulated, and the sample unit.
- Second sentence: encodings (axes, color, panels) and the statistical marks actually drawn.
- Last sentence: the file name, or the table's source path, so the caption can be regenerated.

A caption does not interpret beyond the encodings. Interpretation belongs in results.

## Plain-language summary

Five to eight sentences. No new numbers. Define the sample, what was compared, what was found in the outputs, and what was not measured. Keep identifiers (gene symbols, NCT numbers, file names) as they appear in the source.
