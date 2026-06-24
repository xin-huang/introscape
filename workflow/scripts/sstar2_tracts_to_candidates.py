# Copyright 2026 Xin Huang
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

import csv
from pathlib import Path


FIELDNAMES = [
    "analysis",
    "method",
    "dataset",
    "chrom",
    "start",
    "end",
    "score",
    "sample",
    "attributes",
]


def parse_tract_line(line: str, analysis: str, dataset: str) -> dict[str, str]:
    """
    Parse one sstar2 inferred-tract BED line.

    Parameters
    ----------
    line : str
        BED line.
    analysis : str
        Analysis identifier.
    dataset : str
        Dataset identifier.

    Returns
    -------
    dict[str, str]
        Candidate-region row.
    """
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 3:
        raise ValueError(f"Expected at least three BED fields, got: {line!r}")
    sample = fields[3] if len(fields) > 3 else ""
    score = fields[4] if len(fields) > 4 else ""
    attributes = "\t".join(fields[5:]) if len(fields) > 5 else ""
    return {
        "analysis": analysis,
        "method": "sstar2",
        "dataset": dataset,
        "chrom": fields[0],
        "start": fields[1],
        "end": fields[2],
        "score": score,
        "sample": sample,
        "attributes": attributes,
    }


def convert_tracts(
    input_path: str,
    output_path: str,
    analysis: str,
    dataset: str,
) -> None:
    """
    Convert inferred tracts to the common candidate-region table.

    Parameters
    ----------
    input_path : str
        sstar2 BED output path.
    output_path : str
        Candidate-region TSV path.
    analysis : str
        Analysis identifier.
    dataset : str
        Dataset identifier.
    """
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(input_path, "r") as source, output.open("w", newline="") as target:
        writer = csv.DictWriter(target, fieldnames=FIELDNAMES, delimiter="\t")
        writer.writeheader()
        for line in source:
            if not line.strip() or line.startswith("#"):
                continue
            writer.writerow(parse_tract_line(line, analysis, dataset))


convert_tracts(
    str(snakemake.input.tracts),
    str(snakemake.output.candidates),
    str(snakemake.wildcards.analysis),
    str(snakemake.params.dataset),
)
