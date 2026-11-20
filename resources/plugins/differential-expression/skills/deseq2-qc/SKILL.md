---
name: deseq2-qc
description: >-
  Design-aware preflight of sample metadata, read-only verification of the
  wrapper-produced QC certificate for real DESeq2 results TSV, and
  approval-gated one-sided hypergeometric ORA with BH correction over an
  explicit GMT.
phi:
  environment: phi:python@1
  attachTo: [DifferentialExpression]
  scripts:
    - name: preflight
      description: >-
        Validate user-described study design against sample metadata before
        any wrapper run. Returns needs_input with a structured question
        (header, question, options with label/description), or wrapper_params
        containing only study_design, sample_id_col, contrast_variable,
        reference_level, target_level, and donor_col when paired.
      run: [python, ./scripts/preflight.py]
      approval: read
      timeoutSeconds: 60
      args:
        type: object
        additionalProperties: false
        required: [metadata, design, sample_id_col, condition_col, reference_level, target_level]
        properties:
          metadata:
            type: string
            format: input-path
            description: Sample metadata TSV or CSV.
          design:
            type: string
            enum: [independent, paired, unknown]
            description: >-
              Study design as described by the user. Never inferred. Use
              unknown when the user has not described it.
          sample_id_col:
            type: string
            description: Column holding unique sample identifiers.
          condition_col:
            type: string
            description: Column holding group/condition labels.
          reference_level:
            type: string
            description: Reference (control) level. Must exist in condition_col.
          target_level:
            type: string
            description: Target (treatment) level. Must exist in condition_col.
          donor_col:
            type: string
            description: >-
              Donor/pair identifier column. Required for paired designs;
              donors are never dropped.
    - name: assess
      description: >-
        Read-only verification of the wrapper-produced research_qc sidecar
        against the exact results bytes. Recomputes evaluated/hit counts
        using the thresholds recorded in the certificate (hit when
        padj < padjThreshold and |lfc| >= log2fcThreshold) and checks
        resultMd5 and disposition. Never creates, overwrites, or
        recertifies the sidecar and accepts no new thresholds.
      run: [python, ./scripts/assess.py]
      approval: read
      timeoutSeconds: 120
      args:
        type: object
        additionalProperties: false
        required: [results, qc]
        properties:
          results:
            type: string
            format: input-path
            description: DESeq2 results TSV exactly as produced by the wrapper.
          qc:
            type: string
            format: input-path
            description: The research_qc sidecar produced by the wrapper.
    - name: ora
      description: >-
        Over-representation analysis over the certified results and an
        explicit GMT. Revalidates the certificate and result bytes on every
        call; eligible only. One-sided hypergeometric test, BH across all
        tested terms. Never computes or creates output for zero/invalid/
        missing input. Never overwrites existing unrelated files. Invoke
        only after explicit user write approval.
      run: [python, ./scripts/ora.py]
      approval: write
      timeoutSeconds: 300
      args:
        type: object
        additionalProperties: false
        required: [results, qc, gmt, out]
        properties:
          results:
            type: string
            format: input-path
            description: The same results TSV that was certified.
          qc:
            type: string
            format: input-path
            description: The research_qc sidecar produced by the wrapper.
          gmt:
            type: string
            format: input-path
            description: User-provided GMT gene-set file. No external catalogs.
          out:
            type: string
            format: project-path
            description: >-
              Output TSV path. Must not collide with existing unrelated
              files; choose a fresh path instead of overwriting.
---

# deseq2-qc

Three typed tools: `de_preflight`, `de_assess`, `de_ora`. Arguments are JSON
objects matching each `args` schema. Every scientific branch (ok,
needs_input, skip, blocked, invalid) emits one compact machine-readable JSON
object to stdout with exit 0 — never rows, tables, top terms, or logs.
Exit 1 is reserved for unexpected errors.

- `de_preflight` (read): gate before any wrapper run. `needs_input` carries
  `question = {header, question, options: [{label, description}...]}`; relay
  it verbatim and do not start the wrapper. `ok` carries `wrapper_params`
  with only `study_design`, `sample_id_col`, `contrast_variable`,
  `reference_level`, `target_level`, and `donor_col` when paired. Donors are
  never dropped; there is no unpaired fallback.
- `de_assess` (read): verifies the wrapper-produced `research_qc` sidecar —
  a UTF-8 TSV with exactly one header and one record, columns in order:
  `schemaVersion disposition contrast geneIdColumn padjThreshold
  log2fcThreshold evaluated hits resultMd5`. `schemaVersion` is literally
  `1`; `resultMd5` is the MD5 of the exact results bytes (integrity, not
  security). Hit rule is strict: padj < padjThreshold and
  |lfc| >= log2fcThreshold, exactly as recorded. Disposition: `eligible`
  iff evaluated > 0 and hits > 0; `skip` iff evaluated > 0 and hits = 0;
  `invalid` when evaluated = hits = 0 or the certificate does not match the
  bytes. NA is valid missing; all-NA is invalid. Assess never writes.
- `de_ora` (write, explicit approval): eligible certificates only, explicit
  GMT and required `out` project path, revalidation on every call, atomic
  output, no overwrite of existing unrelated files.