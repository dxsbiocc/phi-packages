# MAG function

**Purpose.** Turn MAG gene prediction and functional annotation outputs into
per-bin gene counts, MIMAG structural evidence, and bin × KO / bin × COG
matrices — the "what can these genomes do" layer of the assembly route.

**Principle.** prokka `.txt` summaries give contig/base/CDS/rRNA/tRNA counts;
its `.tsv` feature table names rRNA products, so 5S/16S/23S presence plus
tRNA ≥ 18 yields the MIMAG structural criterion that CheckM tables cannot
show. Prodigal `.faa`/`.gff` (gzipped in nf-core/mag) give CDS counts only.
eggNOG-mapper `.emapper.annotations` rows carry KEGG_ko and COG_category per
query protein; the tool counts hits per bin. Bin attribution is
deterministic: prokka subdirectory name (nf-core layout), prodigal/eggnog
file stem cleaned of `_prodigal` and assembler prefixes, `--bin_map`
(file,bin) overrides everything.

**Interpretation.** `mimag_structural` pass/fail/unknown: unknown means no
prokka feature table for that bin — not a failure. Combine with mg_bin_qc
tiers: completeness/contamination HQ + structural pass is the full MIMAG HQ
claim (at prediction level). KO/COG matrices are functional potential:
hit counts scale with genome size and annotation depth, so compare profiles
(row shares) rather than raw counts across bins of different size.

**Pitfalls.** rRNA/tRNA evidence is prokka's prediction, not barrnap/
tRNAscan-SE confirmation — state the level. When both prokka and prodigal
exist for a bin, prokka's CDS count wins (prokka runs prodigal internally;
summing would double-count). One prodigal file per bin — passing both .faa
and .gff of the same bin takes the max, not the sum. COG categories can be
multi-letter (e.g. "KT"); each letter is counted separately.