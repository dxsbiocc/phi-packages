"""Per-entity normalization shared by DepMap, PRISM, GDSC, and Cell Model Passports.

Gene, cell-line, and compound modules run before identifier mapping. They do not
guess when more than one reference record matches.
"""

from pgx.normalize.cell_lines import CellosaurusIndex, canonical_name
from pgx.normalize.compounds import name_variants, resolve_compounds
from pgx.normalize.genes import HgncIndex, parse_gene_header, resolve_gene

__all__ = [
    "CellosaurusIndex",
    "HgncIndex",
    "canonical_name",
    "name_variants",
    "parse_gene_header",
    "resolve_compounds",
    "resolve_gene",
]
