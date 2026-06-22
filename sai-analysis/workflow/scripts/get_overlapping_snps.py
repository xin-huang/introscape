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

windows_file = snakemake.input.windows
u_files = snakemake.input.U_snps
q_files = snakemake.input.Q_snps
out_file = snakemake.output.overlapping

def chrom_from_path(path):
    base = os.path.basename(path)
    return base.split(".chr", 1)[1].split(".", 1)[0]

def to_set(s):
    return {x.strip() for x in s.split(",")} if isinstance(s, str) and s else set()

# map chrom -> file
u_map = {chrom_from_path(f): f for f in u_files}
q_map = {chrom_from_path(f): f for f in q_files}

# load windows (assume header with Chrom/Start/End)
win = pd.read_csv(windows_file, sep="\t")
win = win[["Chrom","Start","End"]].copy()
win["Chrom"] = win["Chrom"].astype(str)

rows = []

for chrom, subw in win.groupby("Chrom"):
    if chrom not in u_map or chrom not in q_map:
        continue

    # read only the needed chromosome's U and Q files
    u_df = pd.read_csv(u_map[chrom], sep="\t")
    q_df = pd.read_csv(q_map[chrom], sep="\t")

    # build window -> SNPs index for this chromosome
    u_idx = {(str(r.Chrom), int(r.Start), int(r.End)): str(r.U_SNP) for _, r in u_df.iterrows()}
    q_idx = {(str(r.Chrom), int(r.Start), int(r.End)): str(r.Q_SNP) for _, r in q_df.iterrows()}

    # for each window on this chromosome, take intersection
    for _, r in subw.iterrows():
        key = (str(r["Chrom"]), int(r["Start"]), int(r["End"]))
        if key not in u_idx or key not in q_idx:
            continue
        u_set = to_set(u_idx[key])
        q_set = to_set(q_idx[key])
        inter = u_set & q_set
        if inter:
            rows.append([key[0], key[1], key[2], ",".join(sorted(inter))])

out_df = pd.DataFrame(rows, columns=["Chrom","Start","End","Overlapping_SNPs"])
out_df["Chrom"] = pd.to_numeric(out_df["Chrom"], errors="coerce")
out_df["Start"] = pd.to_numeric(out_df["Start"], errors="coerce")
out_df["End"]   = pd.to_numeric(out_df["End"], errors="coerce")
out_df = out_df.sort_values(by=["Chrom","Start","End"])
out_df.to_csv(out_file, sep="\t", index=False)
