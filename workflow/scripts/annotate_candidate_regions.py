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
import gzip
import re
from pathlib import Path
from typing import TextIO


ANNOTATED_FIELDNAMES = [
    "analysis",
    "method",
    "dataset",
    "chrom",
    "start",
    "end",
    "score",
    "sample",
    "attributes",
    "gene_id",
    "gene_name",
    "gene_chrom",
    "gene_start",
    "gene_end",
    "gene_biotype",
]
GENE_FIELDNAMES = [
    "analysis",
    "method",
    "dataset",
    "gene_id",
    "gene_name",
    "gene_chrom",
    "gene_start",
    "gene_end",
    "gene_biotype",
]


def open_text(path: str) -> TextIO:
    """
    Open plain text or gzip-compressed text.

    Parameters
    ----------
    path : str
        Input file path.

    Returns
    -------
    Iterable[str]
        Text lines.
    """
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_gtf_attributes(attributes: str) -> dict[str, str]:
    """
    Parse a GTF attribute string.

    Parameters
    ----------
    attributes : str
        Raw GTF attributes.

    Returns
    -------
    dict[str, str]
        Parsed attributes.
    """
    parsed = {}
    for item in attributes.strip().split(";"):
        item = item.strip()
        if not item or " " not in item:
            continue
        key, value = item.split(" ", 1)
        parsed[key] = value.strip().strip('"')
    return parsed


def normalize_chrom(chrom: str) -> str:
    """
    Normalize common chromosome naming styles.

    Parameters
    ----------
    chrom : str
        Chromosome name.

    Returns
    -------
    str
        Normalized chromosome token.
    """
    token = chrom.strip()
    if token.lower().startswith("chr"):
        token = token[3:]
    match = re.match(r"NC_0*(\d+)\.", token)
    if match:
        number = int(match.group(1))
        if number == 23:
            return "X"
        if number == 24:
            return "Y"
        return str(number)
    return token


def read_genes(path: str) -> dict[str, list[dict[str, str]]]:
    """
    Read gene features from a GTF file.

    Parameters
    ----------
    path : str
        GTF file path.

    Returns
    -------
    dict[str, list[dict[str, str]]]
        Gene features keyed by normalized chromosome.
    """
    genes: dict[str, list[dict[str, str]]] = {}
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = parse_gtf_attributes(fields[8])
            gene_id = attrs.get("gene_id", "")
            gene_name = attrs.get("gene", attrs.get("gene_name", gene_id))
            gene = {
                "gene_id": gene_id,
                "gene_name": gene_name,
                "gene_chrom": fields[0],
                "gene_start": str(int(fields[3]) - 1),
                "gene_end": fields[4],
                "gene_biotype": attrs.get("gene_biotype", ""),
            }
            genes.setdefault(normalize_chrom(fields[0]), []).append(gene)
    return genes


def overlaps(candidate: dict[str, str], gene: dict[str, str]) -> bool:
    """
    Return whether one candidate region overlaps one gene.

    Parameters
    ----------
    candidate : dict[str, str]
        Candidate-region row.
    gene : dict[str, str]
        Gene row.

    Returns
    -------
    bool
        True when intervals overlap.
    """
    candidate_start = int(candidate["start"])
    candidate_end = int(candidate["end"])
    gene_start = int(gene["gene_start"])
    gene_end = int(gene["gene_end"])
    return candidate_start < gene_end and candidate_end > gene_start


def annotate_candidates(
    candidates_path: str,
    annotation_path: str,
    annotated_path: str,
    genes_path: str,
) -> None:
    """
    Annotate candidate regions with overlapping genes.

    Parameters
    ----------
    candidates_path : str
        Candidate-region TSV path.
    annotation_path : str
        GTF annotation path.
    annotated_path : str
        Annotated candidate-region TSV path.
    genes_path : str
        Candidate gene TSV path.
    """
    genes_by_chrom = read_genes(annotation_path)
    Path(annotated_path).parent.mkdir(parents=True, exist_ok=True)
    Path(genes_path).parent.mkdir(parents=True, exist_ok=True)
    seen_genes = set()

    with open(candidates_path, "r", newline="") as candidates_file:
        reader = csv.DictReader(candidates_file, delimiter="\t")
        candidates = list(reader)

    with open(annotated_path, "w", newline="") as annotated_file, open(
        genes_path, "w", newline=""
    ) as genes_file:
        annotated_writer = csv.DictWriter(
            annotated_file, fieldnames=ANNOTATED_FIELDNAMES, delimiter="\t"
        )
        gene_writer = csv.DictWriter(
            genes_file, fieldnames=GENE_FIELDNAMES, delimiter="\t"
        )
        annotated_writer.writeheader()
        gene_writer.writeheader()

        for candidate in candidates:
            candidate_genes = genes_by_chrom.get(normalize_chrom(candidate["chrom"]), [])
            matched = [gene for gene in candidate_genes if overlaps(candidate, gene)]
            if not matched:
                annotated_writer.writerow(
                    {
                        **candidate,
                        "gene_id": "",
                        "gene_name": "",
                        "gene_chrom": "",
                        "gene_start": "",
                        "gene_end": "",
                        "gene_biotype": "",
                    }
                )
                continue

            for gene in matched:
                annotated_writer.writerow({**candidate, **gene})
                gene_key = (gene["gene_id"], gene["gene_chrom"])
                if gene_key in seen_genes:
                    continue
                seen_genes.add(gene_key)
                gene_writer.writerow(
                    {
                        "analysis": candidate["analysis"],
                        "method": candidate["method"],
                        "dataset": candidate["dataset"],
                        **gene,
                    }
                )


annotate_candidates(
    str(snakemake.input.candidates),
    str(snakemake.input.annotation),
    str(snakemake.output.annotated),
    str(snakemake.output.genes),
)
