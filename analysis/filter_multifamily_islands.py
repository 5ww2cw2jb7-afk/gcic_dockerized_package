#!/usr/bin/env python3
import pandas as pd
import sys

if len(sys.argv) != 3:
    print("Usage: python filter_multifamily_islands.py <in.islands.tsv> <out.multi.tsv>")
    sys.exit(1)

in_tsv = sys.argv[1]
out_tsv = sys.argv[2]

# load
df = pd.read_csv(in_tsv, sep="\t")

required_cols = {"chrom", "start", "end", "family"}
if not required_cols.issubset(df.columns):
    print("[FATAL] islands.tsv must have columns:", required_cols)
    sys.exit(1)

# keep islands with >=2 families
multi = (
    df.groupby(["chrom", "start", "end"])
      .filter(lambda x: x["family"].nunique() >= 2)
)

# write
multi.to_csv(out_tsv, sep="\t", index=False)

print(f"[OK] multi-family islands written to {out_tsv}")
print(f"     total islands (before): {df.groupby(['chrom','start','end']).ngroups}")
print(f"     multi-family islands   : {multi.groupby(['chrom','start','end']).ngroups}")
