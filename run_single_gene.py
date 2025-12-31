#!/usr/bin/env python3
"""Run GCIC pipeline for a single gene extracted from a large multi-FASTA.

Outputs:
- results/per_gene/<GENE_ID>/* (per-gene outputs)
- results/all_genes_gcic_summary.tsv (single-row summary; OVERWRITTEN each run)
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

    # Extract and write per-gene FASTA (input for pipeline)
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

    # Build summary from GCIC.islands.multi.tsv
    islands_multi = per_gene / "GCIC.islands.multi.tsv"
    families = set()
    total_bp = 0
    if islands_multi.exists() and islands_multi.stat().st_size > 0:
        df = pd.read_csv(islands_multi, sep="\t")
        if "family" in df.columns:
            families = set(df["family"].dropna().astype(str))
        if "start" in df.columns and "end" in df.columns:
            total_bp = int((df["end"] - df["start"]).clip(lower=0).sum())

    summary_path = out_root / "all_genes_gcic_summary.tsv"
    row = {
        "gene_id": gene_id,
        "species": species,
        "chrom": chrom,
        "region_start": region_start,
        "region_end": region_end,
        "n_families": len(families),
        "total_island_bp": total_bp
    }

    # 重要：毎回上書き（溜めない）
    with open(summary_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(row.keys()), delimiter="\t")
        w.writeheader()
        w.writerow(row)


if __name__ == "__main__":
    main()
