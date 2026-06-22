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


import pandas as pd

f1 = snakemake.input.u_outliers
f2 = snakemake.input.q_outliers
out = snakemake.output.overlapping

df1 = pd.read_csv(f1, sep="\t")
df2 = pd.read_csv(f2, sep="\t")

# take intersection by Chrom, Start, End
common = pd.merge(df1, df2, on=["Chrom","Start","End"])

common.to_csv(out, sep="\t", index=False)
