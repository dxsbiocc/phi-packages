PubMed connects to NCBI E-utilities, PubMed Central (PMC) / Europe PMC open-access JATS XML, and NCBI PubTator 3.0 to search MEDLINE and life-science journal articles, retrieve structured abstracts and open full-text sections, traverse related or citing articles, and extract biomedical entity annotations.

### Capabilities

- Search PubMed/MEDLINE using native MeSH headings and field tags (`[ti]`, `[au]`, `[ta]`, `[pdat]`) as well as PubTator 3.0 entity queries (`@GENE_TP53`, `@DISEASE_`, `@CHEMICAL_`, `@VARIANT_`, `relations:`).
- Retrieve complete bibliographic records, ORCID-linked authors, structured abstracts, and PMID/PMCID/DOI crosswalks.
- Extract Creative Commons or public-domain full-text sections (`abstract`, `methods`, `results`, `discussion`) from PMC / Europe PMC JATS XML.
- Expand from a key paper to related PubMed citations, citing works, or cited references, and batch-resolve raw reference lists.
- Retrieve PubTator 3.0 BioCJSON entity mentions (genes, diseases, chemicals, variants, species, and cell lines) per article.

### Research & Workflow Use Cases

- Build systematic literature reviews with reproducible MeSH and publication-date filters.
- Quote experimental parameters directly from open-access `Methods` sections without downloading PDFs.
- Cross-check variant, drug, or gene relationships annotated across biomedical abstracts via PubTator 3.0.
- Convert mixed lists of PMIDs, PMCIDs, DOIs, and formatted citations into normalized bibliography entries.

### Tools

| Tool | Description |
| --- | --- |
| `usage_terms` | Return usage terms, NCBI E-utilities / PMC / PubTator 3.0 coverage, and rate-limit guidance. |
| `pubmed_search_articles` | Search PubMed biomedical literature with MeSH/field tags or PubTator 3.0 entity queries, date bounds, and article-type filters. |
| `pubmed_get_article` | Retrieve full article metadata, structured abstract, PMID/PMCID/DOI crosswalks, and open-access links. |
| `pubmed_get_fulltext` | Fetch open-access PMC / Europe PMC JATS full text or a specific section when redistribution is permitted. |
| `pubmed_related_articles` | Find similar PubMed citations, citing papers, or referenced works from a seed PMID or DOI. |
| `pubmed_lookup_citations` | Batch-resolve up to 25 PMIDs, PMCIDs, DOIs, or citation strings with identifier conversion. |
| `pubmed_entity_annotations` | Retrieve PubTator 3.0 BioCJSON gene, disease, chemical, variant, and species annotations for an article. |

### Usage Notes

- No sign-in or API key is required. Hosted read-only by [Phi Packages](https://github.com/dxsbiocc/phi-packages) over official public NCBI and Europe PMC APIs.
- Full-text bodies are inlined only when a Creative Commons or public-domain license is confirmed; paywalled or restricted articles return links and metadata only.
- Follow [NCBI Policies and Disclaimers](https://www.ncbi.nlm.nih.gov/home/about/policies/) and cite original authors and NLM resources.
