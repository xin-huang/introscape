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


import pysam
import re


def fa2bed(fasta_path: str, output_bed: str) -> None:
    """
    Converts an ancestral allele FASTA file to BED format.

    Parameters
    ----------
    fasta_path : str
        Path to the ancestral allele FASTA file (supports .gz).
    output_bed : str
        Path to the output BED file (gzip compressed).
    """
    fasta = pysam.FastaFile(fasta_path)

    with open(output_bed, "wt") as out:
        for raw_chrom in fasta.references:
            match = re.search(
                r"GRCh\d+:(\d+|X|Y)", raw_chrom
            )  # Extract chromosome number
            if not match:
                print(f"Skipping unrecognized chromosome name: {raw_chrom}")
                continue
            chrom = match.group(1)  # Keep only the number (no "chr" prefix)

            seq = fasta.fetch(raw_chrom).upper()
            for pos, base in enumerate(seq):
                if base in "ACGT":  # Filter out 'N'
                    out.write(f"{chrom}\t{pos}\t{pos+1}\t{base}\n")

    fasta.close()


fa2bed(snakemake.input.fa, snakemake.output.bed)
