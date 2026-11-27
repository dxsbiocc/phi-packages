Clinical Trials queries the US National Library of Medicine's official ClinicalTrials.gov REST API v2, providing structured access to global interventional and observational clinical study registrations, eligibility criteria, study arms, endpoints, and posted summary results.

### Capabilities

- Search trials by disease/condition, drug or intervention, target/biomarker keyword, lead sponsor, facility location, recruitment status, phase (`PHASE1`–`PHASE4`), study type, and posted-results availability.
- Retrieve full protocol records by NCT identifier (`NCT...`), including study design, masking/allocation, arms, primary and secondary outcome measures, enrollment counts, and linked PubMed citations.
- Inspect complete inclusion and exclusion eligibility criteria, eligible age bounds, sex, and healthy-volunteer rules.
- Read posted outcome measures, statistical analyses, and serious/other adverse-event summaries from `resultsSection` for completed trials.

### Research & Workflow Use Cases

- Landscape clinical development pipelines for a therapeutic target, drug modality, or oncology indication.
- Compare inclusion/exclusion criteria, endpoints, and sample sizes when designing a new clinical protocol or statistical analysis plan.
- Extract posted primary-endpoint and safety summaries across completed Phase 2/3 trials for meta-analysis and evidence synthesis.

### Tools

| Tool | Description |
| --- | --- |
| `usage_terms` | Return usage terms and NLM ClinicalTrials.gov REST API v2 notes. |
| `ctgov_search_studies` | Search ClinicalTrials.gov studies by condition, intervention, keyword, sponsor, location, status, phase, and posted results. |
| `ctgov_get_study` | Retrieve the full study protocol, design, arms, endpoints, enrollment, locations, and linked PMIDs for one NCT ID. |
| `ctgov_get_eligibility` | Retrieve full inclusion/exclusion criteria, age range, sex, and healthy-volunteer flags for one NCT ID. |
| `ctgov_get_results` | Retrieve posted outcome measures, statistical analyses, and adverse-event summaries for completed trials with results. |

### Usage Notes

- No sign-in or API key is required. Hosted read-only by [Phi Packages](https://github.com/dxsbiocc/phi-packages) over the official [ClinicalTrials.gov API v2](https://clinicaltrials.gov/data-api/about-api).
- Registry records are provided for scientific research and are not medical advice. Always cite `ClinicalTrials.gov` and the corresponding `NCT` identifier.
