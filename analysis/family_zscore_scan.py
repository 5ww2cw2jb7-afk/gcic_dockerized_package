#!/usr/bin/env python3
import argparse
from pyfaidx import Fasta
import regex as re
import pandas as pd
import numpy as np
from collections import defaultdict
from scipy.stats import zscore

IUPAC = {
    "A":"A","C":"C","G":"G","T":"T",
    "R":"[AG]","Y":"[CT]","S":"[GC]","W":"[AT]",
    "K":"[GT]","M":"[AC]","B":"[CGT]",
    "D":"[AGT]","H":"[ACT]","V":"[ACG]","N":"[ACGT]"
}

def iupac_to_regex(seq):
    return "".join(IUPAC.get(c, c) for c in seq.upper())

def load_motifs(path):
    motifs = []
    with open(path) as f:
        for ln in f:
            if ln.startswith("#") or not ln.strip():
                continue
            name, seq = ln.split()[:2]
            pat = re.compile(f"(?={iupac_to_regex(seq)})", re.I)
            motifs.append((name, pat))
    return motifs

def load_family_map(path):
    d = {}
    df = pd.read_csv(path, sep="\t")
    for _, r in df.iterrows():
        d[r['motif']] = r['family']
    return d

def scan(seq, motifs):
    hits = []
    for name, pat in motifs:
        for m in pat.finditer(seq, overlapped=True):
            hits.append((m.start(), name))
    return hits

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--motifs", required=True)
    ap.add_argument("--family-map", required=True)
    ap.add_argument("--window", type=int, default=200)
    ap.add_argument("--step", type=int, default=50)
    ap.add_argument("--z", type=float, default=2.5)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    motifs = load_motifs(args.motifs)
    fmap = load_family_map(args.family_map)
    fasta = Fasta(args.fasta)

    records = []

    for chrom in fasta.keys():
        seq = str(fasta[chrom]).upper()
        hits = scan(seq, motifs)

        fam_hits = defaultdict(list)
        for pos, m in hits:
            if m in fmap:
                fam_hits[fmap[m]].append(pos)

        for start in range(0, len(seq)-args.window+1, args.step):
            end = start + args.window
            for fam, poss in fam_hits.items():
                c = sum(start <= p < end for p in poss)
                records.append((chrom, start, end, fam, c))

    df = pd.DataFrame(records, columns=["chrom","start","end","family","count"])

    df['z'] = df.groupby('family')['count'].transform(zscore)
    sig = df[df['z'] >= args.z]

    sig.to_csv(f"{args.out_prefix}.family_zscore.tsv", sep="\t", index=False)

    # island merge
    islands = []
    for (chrom, fam), g in sig.groupby(['chrom','family']):
        g = g.sort_values('start')
        cs, ce = None, None
        for _, r in g.iterrows():
            if cs is None:
                cs, ce = r.start, r.end
            elif r.start <= ce:
                ce = r.end
            else:
                islands.append((chrom, cs, ce, fam))
                cs, ce = r.start, r.end
        if cs:
            islands.append((chrom, cs, ce, fam))

    pd.DataFrame(islands, columns=["chrom","start","end","family"])\
      .to_csv(f"{args.out_prefix}.islands.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()
