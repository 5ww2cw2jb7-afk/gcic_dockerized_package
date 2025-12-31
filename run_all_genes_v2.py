#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path
import pandas as pd
import csv


DEFAULT_MOTIFS = "resources/motifs.txt"
DEFAULT_FAMILY_MAP = "resources/family_map.tsv"

RICE_FASTA = "gene_regions.5prime.fa"
RICE_BED   = "gene_regions.bed"
AT_FASTA   = "ATgene_regions.5prime.fa"
AT_BED     = "ATgene_regions.bed"


def run(cmd):
    subprocess.run(cmd, check=True)


def load_gene_coordinates(bed_file: str):
    coords = {}
    with open(bed_file) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, gene_id, *_ = line.strip().split("\t")
            coords[gene_id] = (chrom, int(start), int(end))
    return coords


def iter_fasta_records(fasta_path: Path):
    gene_id = None
    seq_lines = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                if gene_id is not None:
                    yield gene_id, "".join(seq_lines)
                gene_id = line[1:].strip().split()[0]
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if gene_id is not None:
            yield gene_id, "".join(seq_lines)


def write_single_fasta(out_fa: Path, gene_id: str, seq: str) -> None:
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    with open(out_fa, "w") as out:
        out.write(f">{gene_id}\n")
        for i in range(0, len(seq), 60):
            out.write(seq[i:i+60] + "\n")


def make_multi_island_fasta(islands_multi_tsv: Path, gene_id: str, seq: str, out_handle, mode: str) -> None:
    df = pd.read_csv(islands_multi_tsv, sep="\t")
    if df.empty:
        return

    for c in ("start", "end", "family"):
        if c not in df.columns:
            return

    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["family"] = df["family"].astype(str)

    grouped = []
    for (s, e), fams in df.groupby(["start", "end"])["family"]:
        fam_set = sorted(set(",".join(fams.astype(str)).split(",")))
        grouped.append((int(s), int(e), fam_set))
    grouped.sort(key=lambda x: (x[0], x[1]))

    for s, e, fam_set in grouped:
        s0 = max(0, s)
        e0 = min(len(seq), e)
        island_seq = seq[s0:e0]
        fam_str = ",".join(fam_set)
        header = f">{gene_id}_{s}-{e} |families={fam_str} |strand_support=. |mode={mode} |src=+ |scan={s}-{e}"
        out_handle.write(header + "\n")
        for i in range(0, len(island_seq), 60):
            out_handle.write(island_seq[i:i+60] + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", choices=["rice", "arabidopsis"], default="rice")
    ap.add_argument("--strand-mode", choices=["fwd", "rev", "total"], default="fwd")
    ap.add_argument("--gcic-script", default="GCIC_pipeline.py")
    ap.add_argument("--motifs", default=DEFAULT_MOTIFS)
    ap.add_argument("--family-map", dest="family_map", default=DEFAULT_FAMILY_MAP)
    ap.add_argument("--gene-id", action="append", help="Run only these gene IDs (repeatable). If omitted, run all.")
    ap.add_argument("--out-root", default="results")
    args = ap.parse_args()

    if args.species == "arabidopsis":
        fasta_path = Path(AT_FASTA)
        bed_path = Path(AT_BED)
    else:
        fasta_path = Path(RICE_FASTA)
        bed_path = Path(RICE_BED)

    if not fasta_path.exists():
        raise SystemExit(f"[ERROR] FASTA not found: {fasta_path}")
    if not bed_path.exists():
        raise SystemExit(f"[ERROR] BED not found: {bed_path}")

    coords = load_gene_coordinates(str(bed_path))
    targets = set(args.gene_id) if args.gene_id else None

    out_root = Path(args.out_root)
    per_gene_dir = out_root / "per_gene"
    per_gene_dir.mkdir(parents=True, exist_ok=True)

    summary_path = out_root / "all_genes_gcic_summary.tsv"
    islands_fa_path = out_root / "all_genes.GCIC.multi.island.fa"

    cols = [
        "gene_id",
        "chrom",
        "region_start",
        "region_end",
        "has_gcic",
        "gcic_motif_family_count",
        "gcic_total_bp",
        "family_set",
    ]

    with open(summary_path, "w", newline="") as fsum, open(islands_fa_path, "w") as ffa:
        w = csv.DictWriter(fsum, fieldnames=cols, delimiter="\t")
        w.writeheader()

        for gid, seq in iter_fasta_records(fasta_path):
            if gid not in coords:
                continue
            if targets is not None and gid not in targets:
                continue

            chrom, region_start, region_end = coords[gid]
            gene_dir = per_gene_dir / gid
            gene_dir.mkdir(parents=True, exist_ok=True)

            gene_fa = gene_dir / f"{gid}.fa"
            write_single_fasta(gene_fa, gid, seq)

            run([
                "python", args.gcic_script,
                "--species", args.species,
                "--fasta", str(gene_fa),
                "--gene-id", gid,
                "--chrom", str(chrom),
                "--start", str(region_start),
                "--end", str(region_end),
                "--motifs", args.motifs,
                "--family-map", args.family_map,
                "--outdir", str(gene_dir),
                "--strand-mode", args.strand_mode
            ])

            islands_multi = gene_dir / "GCIC.islands.multi.tsv"
            families = set()
            total_bp = 0
            if islands_multi.exists() and islands_multi.stat().st_size > 0:
                df = pd.read_csv(islands_multi, sep="\t")
                if not df.empty and all(c in df.columns for c in ["start", "end", "family"]):
                    for (s, e), sub in df.groupby(["start", "end"]):
                        s, e = int(s), int(e)
                        total_bp += (e - s)
                        fams = ",".join(sub["family"].astype(str)).split(",")
                        families.update([x for x in fams if x])

            has_gcic = len(families) > 0
            w.writerow({
                "gene_id": gid,
                "chrom": chrom,
                "region_start": region_start,
                "region_end": region_end,
                "has_gcic": has_gcic,
                "gcic_motif_family_count": len(families),
                "gcic_total_bp": total_bp,
                "family_set": ",".join(sorted(families)),
            })

            if has_gcic and islands_multi.exists() and islands_multi.stat().st_size > 0:
                make_multi_island_fasta(islands_multi, gid, seq, ffa, args.strand_mode)

    print("[INFO] Done.")
    print(f"[INFO] summary: {summary_path}")
    print(f"[INFO] islands:  {islands_fa_path}")


if __name__ == "__main__":
    main()
