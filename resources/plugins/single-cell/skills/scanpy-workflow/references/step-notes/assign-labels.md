Purpose: write the judged labels back to the object. Principle: the automatic
call is a proposal; the dot plot and the marker tables are the evidence, and
overriding them by hand (or via a mapping JSON) is the analysis, not a failure
of it. Interpretation: report how many clusters got each label and how many
stayed unassigned — unassigned is information about missing evidence. Pitfalls:
do not propagate a coarse label into subtype claims, and do not rename clusters
to match a reference the tissue does not contain.
