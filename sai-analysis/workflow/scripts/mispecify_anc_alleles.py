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
import numpy as np

infile = snakemake.input.anc_alleles
outfile = snakemake.output.anc_alleles
prop = float(snakemake.wildcards.prop)
seed = int(snakemake.params.seed)

df = pd.read_csv(infile, sep="\t", header=None,
                 names=["chrom","start","end","allele"])

n = len(df)
k = int(n * prop)

rng = np.random.default_rng(seed)
pick = rng.choice(df.index, size=k, replace=False)
df.loc[pick, "allele"] = "T"

df.to_csv(outfile, sep="\t", header=False, index=False)
