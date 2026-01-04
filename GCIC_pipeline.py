#!/usr/bin/env python3
import argparse
import subprocess
import pandas as pd
import os
import sys
from pathlib import Path


def run(cmd):
    subprocess.run(cmd, check=True)


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def read_single_fasta(path: str):
    header = None
    seq_lines = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:].strip().split()[0]
            else:
                seq_lines.append(line)
    if header is None:
        raise ValueError(f"Invalid FASTA (no header): {path}")
    return header, "".join(seq_lines)


def write_fasta(path: str, header: str, seq: str):
    with open(path, "w") as out:
        out.write(f">{header}\n")
        for i in range(0, len(seq), 60):
            out.write(seq[i:i + 60] + "\n")


def infer_species_from_gene_id(gene_id: str) -> str:
    if gene_id and gene_id.upper().startswith("AT"):
        return "arabidopsis"
    return "rice"


def window_step_for_species(species: str):
    # 要件：現状の初期値を維持
    if species == "arabidopsis":
        return 100, 25
    return 200, 50


def resolve_if_relative(p: str, base: Path) -> str:
    pp = Path(p)
    if pp.is_absolute():
        return str(pp)
    cand1 = (base / pp).resolve()
    if cand1.exists():
        return str(cand1)
    return str((Path.cwd() / pp).resolve())


def _parse_family_tokens(fam_value) -> list[str]:
    """Comma-separated family string -> cleaned list of tokens."""
    if pd.isna(fam_value):
        return []
    toks = [x.strip() for x in str(fam_value).split(",")]
    return [t for t in toks if t != ""]


def main(args):
    pkg_root = Path(__file__).resolve().parent
    analysis_dir = pkg_root / "analysis"

    zscore_scan = analysis_dir / "family_zscore_scan.py"
    filter_multi = analysis_dir / "filter_multifamily_islands.py"

    if not zscore_scan.exists():
        raise FileNotFoundError(f"Missing script: {zscore_scan}")
    if not filter_multi.exists():
        raise FileNotFoundError(f"Missing script: {filter_multi}")

    os.makedirs(args.outdir, exist_ok=True)
    prefix = os.path.join(args.outdir, "GCIC")

    species = args.species if args.species else infer_species_from_gene_id(args.gene_id)

    default_window, default_step = window_step_for_species(species)
    window = args.window if args.window is not None else default_window
    step = args.step if args.step is not None else default_step

    motifs_path = resolve_if_relative(args.motifs, pkg_root)
    family_map_path = resolve_if_relative(args.family_map, pkg_root)

    # strand-mode：fwdは無加工、rev/totalは一時FASTAを作る
    fasta_for_scan = args.fasta
    tmp_fa = None

    gene_id_in_fa, seq_fwd = read_single_fasta(args.fasta)
    seq_for_island = seq_fwd  # island FASTA 用（fwd/rev/total と整合させる）

    if args.strand_mode != "fwd":
        tmp_fa = prefix + ".tmp.strand_mode.fa"
        if args.strand_mode == "rev":
            seq_for_island = revcomp(seq_fwd)
        elif args.strand_mode == "total":
            seq_for_island = seq_fwd + revcomp(seq_fwd)
        else:
            raise ValueError(f"Unknown strand-mode: {args.strand_mode}")

        write_fasta(tmp_fa, args.gene_id, seq_for_island)
        fasta_for_scan = tmp_fa

    # ① family z-score scan
    run([
        sys.executable, str(zscore_scan),
        "--fasta", fasta_for_scan,
        "--motifs", motifs_path,
        "--family-map", family_map_path,
        "--window", str(window),
        "--step", str(step),
        "--z", "2.5",
        "--out-prefix", prefix
    ])

    islands_tsv = prefix + ".islands.tsv"
    multi_tsv = prefix + ".islands.multi.tsv"

    # ② multi-family islands 抽出
    if os.path.exists(islands_tsv) and os.path.getsize(islands_tsv) > 0:
        run([sys.executable, str(filter_multi), islands_tsv, multi_tsv])
    else:
        # 空ファイルでも作っておく（下流互換）
        open(multi_tsv, "w").close()

    # ③ GCIC island FASTA を生成（UIで表示するファイル）
    island_fa_path = os.path.join(args.outdir, "GCIC.multi.island.fa")

    if (not os.path.exists(multi_tsv)) or os.path.getsize(multi_tsv) == 0:
        # island が無い場合は空にする
        open(island_fa_path, "w").close()
    else:
        df = pd.read_csv(multi_tsv, sep="\t")
        # 期待カラム：chrom,start,end,family
        required = {"start", "end", "family"}
        if not required.issubset(set(df.columns)):
            open(island_fa_path, "w").close()
        else:
            # start/end は 0-based, endはexclusive（family_zscore_scanの定義）
            df["start"] = df["start"].astype(int)
            df["end"] = df["end"].astype(int)
            df["family"] = df["family"].astype(str)

            # 同じ区間をまとめる（FASTAのヘッダ用にユニーク集合を表示）
            groups = (
                df.groupby(["start", "end"])["family"]
                  .apply(lambda s: ",".join(sorted(set(",".join(s).split(",")))))
                  .reset_index()
            )

            with open(island_fa_path, "w") as out:
                for i, r in groups.iterrows():
                    s = int(r["start"])
                    e = int(r["end"])
                    fams = r["family"]
                    # 安全チェック
                    s2 = max(0, min(s, len(seq_for_island)))
                    e2 = max(0, min(e, len(seq_for_island)))
                    if e2 <= s2:
                        continue
                    island_seq = seq_for_island[s2:e2]
                    header = f"{args.gene_id}|island{i+1}|{s2}-{e2}|families={fams}"
                    out.write(f">{header}\n")
                    for j in range(0, len(island_seq), 60):
                        out.write(island_seq[j:j+60] + "\n")

    # ④ per-gene summary（互換：gene_gcic_summary.tsv）
    # ---- イネと同じ定義に戻す ----
    # gcic_motif_family_count = family の総数（重複込み：usage）
    # family_count            = family の種類数（ユニーク：diversity）
    families_set = set()
    total_occ = 0
    total_bp = 0

    if os.path.exists(multi_tsv) and os.path.getsize(multi_tsv) > 0:
        df = pd.read_csv(multi_tsv, sep="\t")
        if all(c in df.columns for c in ["start", "end", "family"]):
            # island 単位（start,end）でまとめて bp と family を集計
            for (s, e), sub in df.groupby(["start", "end"]):
                s, e = int(s), int(e)
                total_bp += max(0, e - s)

                # sub["family"] には複数行があり得るので、行ごとにトークン化して合算
                for fam in sub["family"]:
                    toks = _parse_family_tokens(fam)
                    total_occ += len(toks)        # ★重複込み（usage）
                    families_set.update(toks)     # ★ユニーク集合（diversity）

    has_gcic = (total_occ > 0) or (len(families_set) > 0)

    out_summary = os.path.join(args.outdir, "gene_gcic_summary.tsv")
    summary = {
        "gene_id": args.gene_id,
        "chrom": args.chrom,
        "region_start": args.start,
        "region_end": args.end,
        "has_gcic": has_gcic,

        # ★イネ定義：重複込み総数（usage）
        "gcic_motif_family_count": int(total_occ),

        # 合計長（現行踏襲）
        "gcic_total_bp": int(total_bp),

        # ★イネ定義：種類数（diversity）
        "family_count": int(len(families_set)),

        # ユニーク集合
        "family_set": ",".join(sorted(families_set))
    }
    pd.DataFrame([summary]).to_csv(out_summary, sep="\t", index=False)

    # 後始末
    if tmp_fa and os.path.exists(tmp_fa):
        try:
            os.remove(tmp_fa)
        except OSError:
            pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="GCIC: Gene-Centered Identification of Cis-regulatory islands")
    parser.add_argument("--species", choices=["rice", "arabidopsis"], default=None)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--gene-id", required=True)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--end", type=int, required=True)
    parser.add_argument("--motifs", required=True)
    parser.add_argument("--family-map", dest="family_map", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--strand-mode", choices=["fwd", "rev", "total"], default="fwd")
    parser.add_argument("--window", type=int, default=None)
    parser.add_argument("--step", type=int, default=None)
    args = parser.parse_args()
    main(args)
