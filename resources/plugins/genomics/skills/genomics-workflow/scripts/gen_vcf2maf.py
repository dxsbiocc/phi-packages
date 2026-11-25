#!/usr/bin/env python3
"""gen_vcf2maf: deterministic converter from VEP / SnpEff / ANNOVAR annotated VCF
(.vcf or .vcf.gz, e.g. from nf-core/sarek) to a GDC-compliant somatic MAF table.

Implemented in pure Python 3 stdlib so it runs identically in plugin:genomics and
lightweight test environments with zero external VCF library dependencies.
"""

from __future__ import annotations

import argparse
import datetime
import gzip
import hashlib
import json
import os
import re
import sys
from pathlib import Path
from typing import Any


# Sequence Ontology (VEP / SnpEff) consequence priority and GDC Variant_Classification mapping.
# Ordered from highest biological impact to lowest so transcript selection picks the most severe coding effect.
SO_PRIORITY: list[tuple[str, str]] = [
    ("transcript_ablation", "Splice_Site"),
    ("splice_acceptor_variant", "Splice_Site"),
    ("splice_donor_variant", "Splice_Site"),
    ("stop_gained", "Nonsense_Mutation"),
    ("frameshift_variant", "FRAME_SHIFT"),
    ("stop_lost", "Nonstop_Mutation"),
    ("start_lost", "Translation_Start_Site"),
    ("initiator_codon_variant", "Translation_Start_Site"),
    ("transcript_amplification", "Intron"),
    ("inframe_insertion", "In_Frame_Ins"),
    ("disruptive_inframe_insertion", "In_Frame_Ins"),
    ("conservative_inframe_insertion", "In_Frame_Ins"),
    ("inframe_deletion", "In_Frame_Del"),
    ("disruptive_inframe_deletion", "In_Frame_Del"),
    ("conservative_inframe_deletion", "In_Frame_Del"),
    ("missense_variant", "Missense_Mutation"),
    ("protein_altering_variant", "Missense_Mutation"),
    ("rare_amino_acid_variant", "Missense_Mutation"),
    ("splice_region_variant", "Splice_Region"),
    ("splice_donor_5th_base_variant", "Splice_Region"),
    ("splice_donor_region_variant", "Splice_Region"),
    ("splice_polypyrimidine_tract_variant", "Splice_Region"),
    ("incomplete_terminal_codon_variant", "Silent"),
    ("start_retained_variant", "Silent"),
    ("stop_retained_variant", "Silent"),
    ("synonymous_variant", "Silent"),
    ("coding_sequence_variant", "Missense_Mutation"),
    ("mature_mirna_variant", "RNA"),
    ("5_prime_utr_variant", "5'UTR"),
    ("5_prime_utr_premature_start_codon_gain_variant", "5'UTR"),
    ("3_prime_utr_variant", "3'UTR"),
    ("non_coding_transcript_exon_variant", "RNA"),
    ("non_coding_transcript_variant", "RNA"),
    ("nc_transcript_variant", "RNA"),
    ("intron_variant", "Intron"),
    ("upstream_gene_variant", "5'Flank"),
    ("downstream_gene_variant", "3'Flank"),
    ("tfbs_ablation", "Targeted_Region"),
    ("tfbs_amplification", "Targeted_Region"),
    ("tf_binding_site_variant", "Targeted_Region"),
    ("regulatory_region_ablation", "Targeted_Region"),
    ("regulatory_region_amplification", "Targeted_Region"),
    ("regulatory_region_variant", "Targeted_Region"),
    ("intergenic_variant", "IGR"),
    ("intergenic_region", "IGR"),
]

SO_RANK = {term: idx for idx, (term, _) in enumerate(SO_PRIORITY)}
SO_TO_MAF = dict(SO_PRIORITY)

ANNOVAR_EXONIC_MAP = {
    "nonsynonymous_snv": "Missense_Mutation",
    "stopgain": "Nonsense_Mutation",
    "stoploss": "Nonstop_Mutation",
    "startloss": "Translation_Start_Site",
    "frameshift_insertion": "Frame_Shift_Ins",
    "frameshift_deletion": "Frame_Shift_Del",
    "frameshift_substitution": "Frame_Shift_Indel",
    "nonframeshift_insertion": "In_Frame_Ins",
    "nonframeshift_deletion": "In_Frame_Del",
    "nonframeshift_substitution": "In_Frame_Ins",
    "synonymous_snv": "Silent",
}

ANNOVAR_FUNC_MAP = {
    "splicing": "Splice_Site",
    "utr5": "5'UTR",
    "utr3": "3'UTR",
    "ncrna_exonic": "RNA",
    "ncrna_intronic": "RNA",
    "ncrna_splicing": "Splice_Site",
    "intronic": "Intron",
    "upstream": "5'Flank",
    "downstream": "3'Flank",
    "intergenic": "IGR",
}

MAF_COLUMNS = [
    "Hugo_Symbol",
    "NCBI_Build",
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Strand",
    "Variant_Classification",
    "Variant_Type",
    "Reference_Allele",
    "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2",
    "Tumor_Sample_Barcode",
    "Matched_Norm_Sample_Barcode",
    "HGVSc",
    "HGVSp_Short",
    "Transcript_ID",
    "Consequence",
    "FILTER",
    "t_depth",
    "t_ref_count",
    "t_alt_count",
    "t_vaf",
    "Context",
]

NORMAL_HINTS = ("normal", "ctrl", "control", "blood", "germline", "pbmc", "buffy")


def emit(payload: dict[str, Any]) -> None:
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    sys.exit(0)


def fail(message: str, code: int = 1) -> None:
    print(json.dumps({"error": message}, indent=2, ensure_ascii=False))
    sys.exit(code)


def open_text(path: str):
    if path.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def write_table_artifact(
    path: str,
    title: str,
    *,
    inputs: tuple[str | None, ...],
    parameters: dict[str, Any],
    rows: int,
    columns: int,
) -> str:
    p = Path(path)
    desc: dict[str, Any] = {
        "contractVersion": "1.0.0",
        "kind": "table",
        "file": p.name,
        "mediaType": "text/tab-separated-values",
        "title": title[:200],
        "table": {"rows": int(rows), "columns": int(columns)},
        "provenance": {
            "createdAt": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "tool": "gen_vcf2maf",
            "skill": "genomics-workflow",
            "parameters": parameters,
        },
    }
    cwd = Path.cwd().resolve()
    prov_in = []
    for raw in inputs:
        if not raw:
            continue
        ip = Path(raw)
        try:
            rel = ip.resolve().relative_to(cwd)
        except ValueError:
            continue
        entry: dict[str, Any] = {"path": str(rel)}
        if ip.is_file():
            entry["sha256"] = hashlib.sha256(ip.read_bytes()).hexdigest()
        prov_in.append(entry)
    if prov_in:
        desc["provenance"]["inputs"] = prov_in
    sidecar = p.with_name(p.name + ".phi-artifact.json")
    tmp = sidecar.with_name(sidecar.name + ".tmp")
    tmp.write_text(json.dumps(desc, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    os.replace(tmp, sidecar)
    return str(p)


class FastaContextReader:
    """Fetch 3-mer genomic context using .fai index if available, else in-memory contig map."""

    def __init__(self, fasta_path: str) -> None:
        self.fasta_path = fasta_path
        self.fai: dict[str, tuple[int, int, int, int]] = {}
        self.seqs: dict[str, str] = {}
        self._fh = None
        fai_path = fasta_path + ".fai"
        if not fasta_path.lower().endswith(".gz") and os.path.isfile(fai_path):
            with open(fai_path, "r", encoding="utf-8") as fh:
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 5:
                        name = parts[0].strip()
                        self.fai[name] = (int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4]))
            self._fh = open(fasta_path, "rb")
        else:
            cur_name = None
            chunks: list[str] = []
            with open_text(fasta_path) as fh:
                for line in fh:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        if cur_name is not None:
                            self.seqs[cur_name] = "".join(chunks).upper()
                        cur_name = line[1:].split()[0]
                        chunks = []
                    else:
                        chunks.append(line)
                if cur_name is not None:
                    self.seqs[cur_name] = "".join(chunks).upper()

    def close(self) -> None:
        if self._fh is not None:
            self._fh.close()

    def _resolve_chrom(self, chrom: str, mapping: dict[str, Any]) -> str | None:
        if chrom in mapping:
            return chrom
        alt = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
        return alt if alt in mapping else None

    def get_triplet(self, chrom: str, pos_1based: int) -> str:
        if pos_1based < 2:
            return ""
        if self._fh is not None and self.fai:
            key = self._resolve_chrom(chrom, self.fai)
            if key is None:
                return ""
            length, offset, line_bases, line_width = self.fai[key]
            if pos_1based + 1 > length:
                return ""
            bases = []
            for p0 in (pos_1based - 2, pos_1based - 1, pos_1based):
                byte_off = offset + (p0 // line_bases) * line_width + (p0 % line_bases)
                self._fh.seek(byte_off)
                bases.append(self._fh.read(1).decode("ascii", errors="ignore").upper())
            tri = "".join(bases)
            return tri if len(tri) == 3 and all(b in "ACGT" for b in tri) else ""
        key = self._resolve_chrom(chrom, self.seqs)
        if key is None:
            return ""
        seq = self.seqs[key]
        if pos_1based + 1 > len(seq):
            return ""
        tri = seq[pos_1based - 2 : pos_1based + 1]
        return tri if len(tri) == 3 and all(b in "ACGT" for b in tri) else ""


def normalize_maf_alleles(pos: int, ref: str, alt: str) -> tuple[int, int, str, str, str]:
    """Convert VCF 1-based POS/REF/ALT into GDC MAF coordinates, alleles, and Variant_Type."""
    r = ref.upper()
    a = alt.upper()
    p = int(pos)

    # Trim identical trailing bases first, then identical leading bases (VCF anchor padding).
    while len(r) > 1 and len(a) > 1 and r[-1] == a[-1]:
        r = r[:-1]
        a = a[:-1]
    while len(r) > 0 and len(a) > 0 and r[0] == a[0]:
        r = r[1:]
        a = a[1:]
        p += 1

    if len(r) == 0 and len(a) > 0:
        # Pure insertion between p-1 and p in GDC MAF convention
        return p - 1, p, "-", a, "INS"
    if len(a) == 0 and len(r) > 0:
        # Pure deletion spanning p .. p + len(r) - 1
        return p, p + len(r) - 1, r, "-", "DEL"
    if len(r) == len(a):
        vtype = {1: "SNP", 2: "DNP", 3: "TNP"}.get(len(r), "ONP")
        return p, p + len(r) - 1, r, a, vtype
    if len(a) > len(r):
        return p, p + len(r) - 1, r, a, "INS"
    return p, p + len(r) - 1, r, a, "DEL"


def resolve_frameshift(vc: str, vtype: str) -> str:
    if vc != "FRAME_SHIFT":
        return vc
    if vtype == "INS":
        return "Frame_Shift_Ins"
    if vtype == "DEL":
        return "Frame_Shift_Del"
    return "Frame_Shift_Del"


def parse_info(info_str: str) -> dict[str, str]:
    if not info_str or info_str == ".":
        return {}
    out: dict[str, str] = {}
    for item in info_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
        elif item:
            out[item] = "true"
    return out


def pick_vep_or_snpeff(
    entries: list[str],
    field_idx: dict[str, int],
    mode: str,
    alt_allele: str,
    vtype: str,
) -> dict[str, str]:
    best_score = (999, 1, 1)
    best: dict[str, str] | None = None

    for raw in entries:
        parts = raw.split("|")
        def get_f(name: str) -> str:
            idx = field_idx.get(name)
            return parts[idx].strip() if idx is not None and idx < len(parts) else ""

        allele_field = get_f("allele")
        # When multi-allelic VCF has multiple ALT entries, match the allele if possible
        if allele_field and allele_field not in (alt_allele, alt_allele[1:] if len(alt_allele) > 1 else alt_allele, "-"):
            if len(entries) > 1 and not alt_allele.endswith(allele_field):
                pass

        cons_raw = get_f("consequence") or get_f("annotation")
        terms = [t.strip().lower() for t in re.split(r"[&,+]", cons_raw) if t.strip()]
        rank = min((SO_RANK.get(t, 900) for t in terms), default=900)
        best_term = min(terms, key=lambda t: SO_RANK.get(t, 900)) if terms else ""
        vc = resolve_frameshift(SO_TO_MAF.get(best_term, "Targeted_Region"), vtype)

        canonical = 0 if get_f("canonical").upper() in ("YES", "1", "TRUE") else 1
        biotype = get_f("biotype") or get_f("transcript_biotype")
        is_coding = 0 if "protein_coding" in biotype.lower() else 1

        symbol = get_f("symbol") or get_f("gene_name") or get_f("gene") or "Unknown"
        tx_id = get_f("feature") or get_f("feature_id")
        hgvsc = get_f("hgvsc") or get_f("hgvs.c")
        hgvsp = get_f("hgvsp") or get_f("hgvs.p")
        if ":" in hgvsc:
            hgvsc = hgvsc.split(":", 1)[1]
        if ":" in hgvsp:
            hgvsp = hgvsp.split(":", 1)[1]
        hgvsp = hgvsp.replace("%3D", "=")

        score = (rank, canonical, is_coding)
        if score < best_score or best is None:
            best_score = score
            best = {
                "Hugo_Symbol": symbol if symbol else "Unknown",
                "Variant_Classification": vc,
                "Consequence": cons_raw,
                "Transcript_ID": tx_id,
                "HGVSc": hgvsc,
                "HGVSp_Short": hgvsp,
                "annotator": mode,
            }

    return best or {
        "Hugo_Symbol": "Unknown",
        "Variant_Classification": "Targeted_Region",
        "Consequence": "",
        "Transcript_ID": "",
        "HGVSc": "",
        "HGVSp_Short": "",
        "annotator": "unannotated",
    }


def parse_annovar(info: dict[str, str], vtype: str) -> dict[str, str]:
    gene = (info.get("Gene.refGene") or info.get("Gene.ensGene") or "Unknown").split(",")[0].split("\\x3b")[0].strip()
    func = (info.get("Func.refGene") or info.get("Func.ensGene") or "").split("\\x3b")[0].strip().lower()
    exonic = (info.get("ExonicFunc.refGene") or info.get("ExonicFunc.ensGene") or "").split("\\x3b")[0].strip().lower()
    aachange = (info.get("AAChange.refGene") or info.get("AAChange.ensGene") or "").split(",")[0]
    tx_id, hgvsc, hgvsp = "", "", ""
    if aachange and ":" in aachange:
        aparts = aachange.split(":")
        if len(aparts) >= 2:
            tx_id = aparts[1]
        for token in aparts:
            if token.startswith("c."):
                hgvsc = token
            elif token.startswith("p."):
                hgvsp = token

    if exonic and exonic in ANNOVAR_EXONIC_MAP:
        vc = ANNOVAR_EXONIC_MAP[exonic]
        if vc == "Frame_Shift_Indel":
            vc = "Frame_Shift_Ins" if vtype == "INS" else "Frame_Shift_Del"
    else:
        vc = ANNOVAR_FUNC_MAP.get(func, "Targeted_Region")

    return {
        "Hugo_Symbol": gene or "Unknown",
        "Variant_Classification": vc,
        "Consequence": exonic or func,
        "Transcript_ID": tx_id,
        "HGVSc": hgvsc,
        "HGVSp_Short": hgvsp,
        "annotator": "annovar",
    }


def extract_sample_depths(
    fmt_keys: list[str], sample_str: str, alt_idx: int, ref_base: str, alt_base: str
) -> dict[str, Any]:
    vals = sample_str.split(":")
    fmap = {k: vals[i] if i < len(vals) else "." for i, k in enumerate(fmt_keys)}
    gt = fmap.get("GT", "")

    # Skip strict homozygous reference calls in multi-sample VCFs
    if gt in ("0/0", "0|0", "./.", ".|.", "."):
        return {"present": False}
    if gt:
        alleles = re.split(r"[/|]", gt)
        if str(alt_idx) not in alleles and len(alleles) >= 1 and all(a.isdigit() or a == "." for a in alleles):
            return {"present": False}

    t_ref, t_alt, t_dp, t_vaf = None, None, None, None

    # Standard GATK/Mutect2/FreeBayes AD + DP + AF
    if "AD" in fmap and fmap["AD"] not in (".", ""):
        ad_parts = fmap["AD"].split(",")
        if len(ad_parts) > alt_idx:
            try:
                t_ref = int(float(ad_parts[0]))
                t_alt = int(float(ad_parts[alt_idx]))
            except ValueError:
                pass
    # Strelka2 SNV tier1 base counts (AU, CU, GU, TU) or indel TAR/TIR
    elif len(ref_base) == 1 and len(alt_base) == 1 and f"{alt_base.upper()}U" in fmap and f"{ref_base.upper()}U" in fmap:
        try:
            t_ref = int(fmap[f"{ref_base.upper()}U"].split(",")[0])
            t_alt = int(fmap[f"{alt_base.upper()}U"].split(",")[0])
        except ValueError:
            pass
    elif "TAR" in fmap and "TIR" in fmap:
        try:
            t_ref = int(fmap["TAR"].split(",")[0])
            t_alt = int(fmap["TIR"].split(",")[0])
        except ValueError:
            pass

    if "DP" in fmap and fmap["DP"] not in (".", ""):
        try:
            t_dp = int(float(fmap["DP"]))
        except ValueError:
            pass
    if t_dp is None and t_ref is not None and t_alt is not None:
        t_dp = t_ref + t_alt

    for af_key in ("AF", "VAF", "FREQ"):
        if af_key in fmap and fmap[af_key] not in (".", ""):
            raw_af = fmap[af_key].split(",")[min(alt_idx - 1, fmap[af_key].count(","))].rstrip("%")
            try:
                val = float(raw_af)
                t_vaf = val / 100.0 if fmap[af_key].endswith("%") or val > 1.0 else val
                break
            except ValueError:
                pass
    if t_vaf is None and t_alt is not None and t_dp and t_dp > 0:
        t_vaf = t_alt / float(t_dp)

    return {
        "present": True,
        "t_depth": str(t_dp) if t_dp is not None else "",
        "t_ref_count": str(t_ref) if t_ref is not None else "",
        "t_alt_count": str(t_alt) if t_alt is not None else "",
        "t_vaf": f"{t_vaf:.4f}" if t_vaf is not None else "",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, help="Annotated somatic VCF (.vcf or .vcf.gz).")
    parser.add_argument("--output_maf", required=True, help="Output GDC MAF TSV path.")
    parser.add_argument("--tumor_sample", default=None, help="Explicit tumor sample column or barcode.")
    parser.add_argument("--normal_sample", default=None, help="Optional matched normal sample column in VCF.")
    parser.add_argument(
        "--include_filtered",
        nargs="?",
        const=True,
        type=lambda x: str(x).lower() not in ("false", "0", "no"),
        default=False,
        help="Include non-PASS VCF records (default false).",
    )
    parser.add_argument(
        "--pass_only",
        nargs="?",
        const=True,
        type=lambda x: str(x).lower() not in ("false", "0", "no"),
        default=None,
        help="CLI alias: retain only PASS / '.' records (default true unless --include_filtered).",
    )
    parser.add_argument("--ref_fasta", default=None, help="Optional reference FASTA to extract 3-mer Context.")
    parser.add_argument("--ncbi_build", default="GRCh38", help="Genome build string (default GRCh38).")
    args = parser.parse_args()
    if args.pass_only is None:
        args.pass_only = not args.include_filtered

    if not os.path.isfile(args.vcf):
        fail(f"VCF file not found: {args.vcf}")
    if args.ref_fasta and not os.path.isfile(args.ref_fasta):
        fail(f"ref_fasta not found: {args.ref_fasta}")

    csq_fields: dict[str, int] = {}
    ann_fields: dict[str, int] = {}
    vcf_samples: list[str] = []
    fasta_reader = FastaContextReader(args.ref_fasta) if args.ref_fasta else None

    out_path = Path(args.output_maf)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    n_vcf_records = 0
    n_filtered_out = 0
    n_maf_rows = 0
    vc_counts: dict[str, int] = {}
    annotators_seen: set[str] = set()
    tumor_samples_seen: set[str] = set()
    warnings: list[str] = []
    target_tumor_cols: list[tuple[int, str]] = []
    normal_col_idx: int | None = None
    matched_norm_barcode = args.normal_sample or "NORMAL"

    try:
        with open_text(args.vcf) as vcf_fh, open(out_path, "w", encoding="utf-8", newline="") as maf_fh:
            maf_fh.write("\t".join(MAF_COLUMNS) + "\n")

            for raw_line in vcf_fh:
                if raw_line.startswith("##"):
                    if raw_line.startswith("##INFO=<ID=CSQ,"):
                        m = re.search(r"Format:\s*([^\"'>]+)", raw_line)
                        if m:
                            csq_fields = {f.strip().lower(): i for i, f in enumerate(m.group(1).split("|"))}
                    elif raw_line.startswith("##INFO=<ID=ANN,"):
                        m = re.search(r"Functional annotations:\s*'([^']+)'", raw_line)
                        if m:
                            ann_fields = {
                                f.strip().lower(): i for i, f in enumerate(m.group(1).split("|"))
                            }
                        else:
                            # Standard SnpEff ANN 16-column spec fallback
                            std_ann = [
                                "allele", "annotation", "annotation_impact", "gene_name", "gene_id",
                                "feature_type", "feature_id", "transcript_biotype", "rank",
                                "hgvs.c", "hgvs.p",
                            ]
                            ann_fields = {k: idx for idx, k in enumerate(std_ann)}
                    continue

                if raw_line.startswith("#CHROM"):
                    header_cols = raw_line.rstrip("\r\n").split("\t")
                    vcf_samples = header_cols[9:] if len(header_cols) > 9 else []
                    if args.normal_sample:
                        if args.normal_sample not in vcf_samples:
                            fail(f"normal_sample '{args.normal_sample}' not in VCF samples: {vcf_samples}")
                        normal_col_idx = vcf_samples.index(args.normal_sample)
                        matched_norm_barcode = args.normal_sample
                    elif len(vcf_samples) == 2 and not args.tumor_sample:
                        # Auto-identify normal column if one sample matches normal naming conventions
                        norm_cands = [
                            i for i, s in enumerate(vcf_samples)
                            if any(h in s.lower() for h in NORMAL_HINTS)
                        ]
                        if len(norm_cands) == 1:
                            normal_col_idx = norm_cands[0]
                            matched_norm_barcode = vcf_samples[normal_col_idx]

                    if args.tumor_sample:
                        if vcf_samples and args.tumor_sample in vcf_samples:
                            target_tumor_cols = [(vcf_samples.index(args.tumor_sample), args.tumor_sample)]
                        elif not vcf_samples:
                            target_tumor_cols = [(-1, args.tumor_sample)]
                        elif len(vcf_samples) == 1:
                            # Allow renaming single VCF sample column to user-supplied tumor barcode
                            target_tumor_cols = [(0, args.tumor_sample)]
                        else:
                            fail(f"tumor_sample '{args.tumor_sample}' not in VCF samples: {vcf_samples}")
                    else:
                        if not vcf_samples:
                            target_tumor_cols = [(-1, Path(args.vcf).stem.replace(".vcf", ""))]
                        else:
                            target_tumor_cols = [
                                (idx, s) for idx, s in enumerate(vcf_samples) if idx != normal_col_idx
                            ]
                    continue

                if not raw_line.strip() or raw_line.startswith("#"):
                    continue

                parts = raw_line.rstrip("\r\n").split("\t")
                if len(parts) < 8:
                    continue
                n_vcf_records += 1
                chrom, pos_s, _vid, ref, alt_str, _qual, flt, info_str = parts[:8]
                if args.pass_only and flt not in ("PASS", "."):
                    n_filtered_out += 1
                    continue

                pos = int(pos_s)
                info = parse_info(info_str)
                fmt_keys = parts[8].split(":") if len(parts) > 8 else []
                sample_fields = parts[9:] if len(parts) > 9 else []

                alts = [a for a in alt_str.split(",") if a and a != "*" and not a.startswith("<")]
                for alt_idx, alt in enumerate(alts, start=1):
                    start_pos, end_pos, maf_ref, maf_alt, vtype = normalize_maf_alleles(pos, ref, alt)

                    if "CSQ" in info and csq_fields:
                        ann = pick_vep_or_snpeff(info["CSQ"].split(","), csq_fields, "vep", alt, vtype)
                    elif "ANN" in info:
                        if not ann_fields:
                            std_ann = [
                                "allele", "annotation", "annotation_impact", "gene_name", "gene_id",
                                "feature_type", "feature_id", "transcript_biotype", "rank",
                                "hgvs.c", "hgvs.p",
                            ]
                            ann_fields = {k: idx for idx, k in enumerate(std_ann)}
                        ann = pick_vep_or_snpeff(info["ANN"].split(","), ann_fields, "snpeff", alt, vtype)
                    elif "Func.refGene" in info or "Gene.refGene" in info or "Func.ensGene" in info:
                        ann = parse_annovar(info, vtype)
                    else:
                        gene_fb = info.get("GENE") or info.get("SYMBOL") or "Unknown"
                        ann = {
                            "Hugo_Symbol": gene_fb,
                            "Variant_Classification": "Targeted_Region",
                            "Consequence": "unannotated",
                            "Transcript_ID": "",
                            "HGVSc": "",
                            "HGVSp_Short": "",
                            "annotator": "unannotated",
                        }
                    annotators_seen.add(ann["annotator"])

                    ctx_3mer = ""
                    if fasta_reader is not None and vtype == "SNP":
                        ctx_3mer = fasta_reader.get_triplet(chrom, start_pos)

                    for s_col_idx, tumor_barcode in target_tumor_cols:
                        if s_col_idx >= 0 and s_col_idx < len(sample_fields):
                            s_info = extract_sample_depths(fmt_keys, sample_fields[s_col_idx], alt_idx, maf_ref, maf_alt)
                            if not s_info["present"]:
                                continue
                        else:
                            # Sites-only VCF fallback: check INFO DP/AF if present
                            info_dp = info.get("DP", "")
                            info_af = info.get("AF", "").split(",")[0] if "AF" in info else ""
                            s_info = {
                                "present": True,
                                "t_depth": info_dp,
                                "t_ref_count": "",
                                "t_alt_count": "",
                                "t_vaf": info_af,
                            }

                        row_vals = [
                            ann["Hugo_Symbol"],
                            args.ncbi_build,
                            chrom,
                            str(start_pos),
                            str(end_pos),
                            "+",
                            ann["Variant_Classification"],
                            vtype,
                            maf_ref,
                            maf_ref,
                            maf_alt,
                            tumor_barcode,
                            matched_norm_barcode,
                            ann["HGVSc"],
                            ann["HGVSp_Short"],
                            ann["Transcript_ID"],
                            ann["Consequence"],
                            flt,
                            s_info["t_depth"],
                            s_info["t_ref_count"],
                            s_info["t_alt_count"],
                            s_info["t_vaf"],
                            ctx_3mer,
                        ]
                        maf_fh.write("\t".join(row_vals) + "\n")
                        n_maf_rows += 1
                        tumor_samples_seen.add(tumor_barcode)
                        vc = ann["Variant_Classification"]
                        vc_counts[vc] = vc_counts.get(vc, 0) + 1
    finally:
        if fasta_reader is not None:
            fasta_reader.close()

    if n_maf_rows == 0:
        warnings.append("0 MAF rows written; check whether all VCF records were filtered out (--pass_only false to include non-PASS)")
    if "unannotated" in annotators_seen:
        warnings.append(
            "VCF records lacked VEP (CSQ), SnpEff (ANN), or ANNOVAR (Func.refGene) annotations; "
            "Variant_Classification defaulted to Targeted_Region — annotate VCF with VEP/SnpEff upstream for coding TMB/oncoplots"
        )
    if args.ref_fasta is None:
        warnings.append("no --ref_fasta supplied: Context column is empty (pass --ref_fasta if running gen_signatures SBS96 downstream)")

    params_meta = {
        "pass_only": args.pass_only,
        "ncbi_build": args.ncbi_build,
        "context_extracted": bool(args.ref_fasta),
    }
    artifact = write_table_artifact(
        str(out_path),
        f"Somatic MAF converted from VCF ({n_maf_rows} variants across {len(tumor_samples_seen)} sample(s))",
        inputs=(args.vcf, args.ref_fasta),
        parameters=params_meta,
        rows=n_maf_rows,
        columns=len(MAF_COLUMNS),
    )

    emit(
        {
            "vcf": os.path.abspath(args.vcf),
            "output_maf": os.path.abspath(out_path),
            "ncbi_build": args.ncbi_build,
            "annotators_detected": sorted(annotators_seen),
            "n_vcf_records": n_vcf_records,
            "n_filtered_out": n_filtered_out,
            "n_maf_variants": n_maf_rows,
            "n_samples": len(tumor_samples_seen),
            "tumor_samples": sorted(tumor_samples_seen),
            "matched_normal": matched_norm_barcode if normal_col_idx is not None or args.normal_sample else None,
            "context_extracted": bool(args.ref_fasta),
            "variant_classifications": vc_counts,
            "artifacts": [artifact],
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
