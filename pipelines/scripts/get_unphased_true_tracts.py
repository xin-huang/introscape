import sys
sys.stderr = open(snakemake.log[0], "w")


import pandas as pd
import pyranges as pr

try:
    true_tracts = pd.read_csv(
        snakemake.input.bed, sep="\t", header=None, names=['Chromosome', 'Start', 'End', 'Sample']
    )
except pd.errors.EmptyDataError:
    open(snakemake.output.bed, 'w').close()
else:
    true_tracts['Sample'] = true_tracts.apply(lambda row: "_".join(row['Sample'].split("_")[0:2]), axis=1)
    true_tracts = pr.PyRanges(true_tracts).merge(by='Sample')
    true_tracts.to_csv(snakemake.output.bed, sep="\t", header=False)
