# Copyright 2025 Xin Huang
#
# GNU General Public License v3.0
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, please see
#
#    https://www.gnu.org/licenses/gpl-3.0.en.html


import os
import pandas as pd

outlier_file = snakemake.input.outliers          
annotation_files = snakemake.input.annotation    
output_file = snakemake.output.annotated_candidates

def chrom_key_from_filename(path: str) -> str:
    base = os.path.basename(path)
    return base.split(".chr", 1)[1].split(".", 1)[0]

# build {chrom: filepath}
anno_map = {chrom_key_from_filename(f): f for f in annotation_files}

# load windows
win = pd.read_csv(outlier_file, sep="\t")
win = win[["Chrom", "Start", "End"]].copy()

# normalize chromosome like 'chr1' -> '1'
win["Chrom"] = win["Chrom"].astype(str).str.replace(r"^chr", "", regex=True)
win["Start"] = pd.to_numeric(win["Start"], errors="coerce").astype("Int64")
win["End"]   = pd.to_numeric(win["End"],   errors="coerce").astype("Int64")
win = win.dropna(subset=["Chrom", "Start", "End"])

filtered = []

for chrom, subw in win.groupby("Chrom"):
    if chrom not in anno_map:
        continue

    ann = pd.read_csv(anno_map[chrom], sep="\t")
    # ensure required columns
    ann["Chr"] = ann["Chr"].astype(str).str.replace(r"^chr", "", regex=True)
    ann["Start"] = pd.to_numeric(ann["Start"], errors="coerce").astype("Int64")
    if "End" not in ann.columns:
        ann["End"] = ann["Start"]
    else:
        ann["End"] = pd.to_numeric(ann["End"], errors="coerce").astype("Int64")
    ann = ann[(ann["Chr"] == chrom) & ann["Start"].notna()]

    if ann.empty:
        continue

    # pre-slice by global bounds
    gmin = int(subw["Start"].min())
    gmax = int(subw["End"].max())
    ann_slice = ann[(ann["Start"] >= gmin) & (ann["Start"] <= gmax)].copy()
    if ann_slice.empty:
        continue

    # window filter: Start in any [Start, End]
    hits = []
    for s0, e0 in zip(subw["Start"].astype(int), subw["End"].astype(int)):
        hits.append(ann_slice[(ann_slice["Start"] >= s0) & (ann_slice["Start"] <= e0)])
    if hits:
        sub = pd.concat(hits, ignore_index=True).drop_duplicates()
        if not sub.empty:
            filtered.append(sub)

# write output
if filtered:
    res = pd.concat(filtered, ignore_index=True)
    res = res.sort_values(by=["Chr", "Start", "End"], ascending=[True, True, True])
    res.to_csv(output_file, sep="\t", index=False)
else:
    # empty file with header only
    pd.DataFrame().to_csv(output_file, sep="\t", index=False)
