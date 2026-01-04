#!/usr/bin/env python3
"""
Run GCIC pipeline for a single gene extracted from a large multi-FASTA.

Outputs:
- results/per_gene/<GENE_ID>/* (per-gene outputs)
- results/all_genes_gcic_summary.tsv (single-row summary; OVERWRITTEN each run)

Key change (align with rice definition/output):
- Do NOT compute summary here.
- Read results/per_gene/<GENE_ID>/gene_gcic_summary.tsv produced by GCIC_pipeline.py
  and write results/all_genes_gcic_summary.tsv with rice-compatible columns.
"""

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Dict, Tuple
import pandas as pd
import csv


DEFAULT_MOTIFS = "resources/motifs.txt"
DEFAULT_FAMILY_MAP = "resources/family_map.tsv"

SUMMARY_COLS = [
    "gene_id",
    "chrom",
    "region_start",
    "region_end",
    "has_gcic",
    "gcic_motif_family_count",
    "gcic_total_bp",
    "family_count",
    "family_set",
]


def run(cmd):
    subprocess.run(cmd, check=True)


def load_gene_coordinates(bed_file: str) -> Dict[str, Tuple[str, int, int]]:
    coords = {}
    with open(bed_file) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            gid = parts[3]
            coords[gid] = (chrom, start, end)
    return coords


def extract_gene_fasta(master_fasta: str, gene_id: str) -> str:
    """Return sequence (string) for the gene_id from master_fasta."""
    seq_lines = []
    found = False
    with open(master_fasta) as f:
        for line in f:
            if line.startswith(">"):
                header_gid = line[1:].strip().split()[0]
                if header_gid == gene_id:
                    found = True
                    seq_lines = []
                else:
                    if found:
                        break
                    found = False
            else:
                if found:
                    seq_lines.append(line.strip())
    if not found or not seq_lines:
        raise RuntimeError(f"Gene {gene_id} not found in FASTA: {master_fasta}")
    return "".join(seq_lines)


def write_single_fasta(out_path: Path, gene_id: str, seq: str):
    out_path.write_text(f">{gene_id}\n{seq}\n")


def normalize_has_gcic(val) -> str:
    """Return 'TRUE'/'FALSE' for rice-compatible output."""
    if isinstance(val, bool):
        return "TRUE" if val else "FALSE"
    s = str(val).strip()
    if s.lower() in {"true", "1", "t", "yes"}:
        return "TRUE"
    return "FALSE"


def read_gene_summary(summary_path: Path) -> dict:
    """Read gene_gcic_summary.tsv and normalize to SUMMARY_COLS."""
    if not summary_path.exists() or summary_path.stat().st_size == 0:
        raise RuntimeError(f"Missing per-gene summary: {summary_path}")

    df = pd.read_csv(summary_path, sep="\t")
    if df.empty:
        raise RuntimeError(f"Empty per-gene summary: {summary_path}")

    row = df.iloc[0].to_dict()

    # Validate required columns
    missing = [c for c in SUMMARY_COLS if c not in row]
    if missing:
        raise RuntimeError(
            f"Per-gene summary missing columns {missing}. "
            f"Found columns: {list(df.columns)}"
        )

    out = {k: row.get(k, "") for k in SUMMARY_COLS}
    out["has_gcic"] = normalize_has_gcic(out["has_gcic"])

    # Make ints where appropriate (safe)
    for k in ("region_start", "region_end", "gcic_motif_family_count", "gcic_total_bp", "family_count"):
        try:
            out[k] = int(out[k])
        except Exception:
            pass

    if pd.isna(out.get("family_set")):
        out["family_set"] = ""

    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", choices=["rice", "arabidopsis"], required=True)
    ap.add_argument("--gene-id", required=True)
    ap.add_argument("--master-fasta", required=True)
    ap.add_argument("--bed", required=True)
    ap.add_argument("--gcic-script", default="GCIC_pipeline.py")
    ap.add_argument("--motifs", default=DEFAULT_MOTIFS)
    ap.add_argument("--family-map", dest="family_map", default=DEFAULT_FAMILY_MAP)
    ap.add_argument("--strand-mode", choices=["fwd", "rev", "total"], default="fwd")
    ap.add_argument("--window", type=int, default=None)
    ap.add_argument("--step", type=int, default=None)
    ap.add_argument("--out-root", default="results")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    gene_id = args.gene_id
    species = args.species

    app_root = Path(".").resolve()
    out_root = (app_root / args.out_root).resolve()
    per_gene = out_root / "per_gene" / gene_id
    per_gene.mkdir(parents=True, exist_ok=True)

    coords = load_gene_coordinates(args.bed)
    if gene_id not in coords:
        raise RuntimeError(f"Gene {gene_id} not found in BED: {args.bed}")
    chrom, region_start, region_end = coords[gene_id]

    # Extract and write per-gene FASTA
    seq = extract_gene_fasta(args.master_fasta, gene_id)
    gene_fa = per_gene / f"{gene_id}.fa"
    write_single_fasta(gene_fa, gene_id, seq)

    # Run pipeline
    cmd = [
        sys.executable, args.gcic_script,
        "--species", species,
        "--fasta", str(gene_fa),
        "--gene-id", gene_id,
        "--chrom", str(chrom),
        "--start", str(region_start),
        "--end", str(region_end),
        "--motifs", args.motifs,
        "--family-map", args.family_map,
        "--outdir", str(per_gene),
        "--strand-mode", args.strand_mode
    ]
    if args.window is not None:
        cmd += ["--window", str(args.window)]
    if args.step is not None:
        cmd += ["--step", str(args.step)]

    run(cmd)

    # Read per-gene summary produced by GCIC_pipeline.py
    per_gene_summary = per_gene / "gene_gcic_summary.tsv"
    row = read_gene_summary(per_gene_summary)

    # Write rice-compatible all_genes_gcic_summary.tsv (OVERWRITE each run)
    summary_path = out_root / "all_genes_gcic_summary.tsv"
    with open(summary_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=SUMMARY_COLS, delimiter="\t")
        w.writeheader()
        w.writerow(row)

    print("[INFO] Done.")
    print(f"[INFO] per-gene: {per_gene}")
    print(f"[INFO] summary:  {summary_path}")


if __name__ == "__main__":
    main()
