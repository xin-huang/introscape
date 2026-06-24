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

from pathlib import Path
from typing import Iterable

import pandas as pd


def read_metadata(path: str) -> pd.DataFrame:
    """
    Read sample metadata.

    Parameters
    ----------
    path : str
        Metadata file path.

    Returns
    -------
    pandas.DataFrame
        Metadata with canonical ``sample`` and ``population`` columns.
    """
    metadata = pd.read_csv(path, sep=r"\s+", dtype=str)
    column_map = {column.lower(): column for column in metadata.columns}
    sample_column = column_map.get("sample")
    population_column = column_map.get("population")
    if sample_column is None or population_column is None:
        raise ValueError("Metadata must contain Sample and Population columns")
    return metadata.rename(
        columns={sample_column: "sample", population_column: "population"}
    )


def collect_samples(metadata: pd.DataFrame, population: str) -> list[str]:
    """
    Collect samples from one population.

    Parameters
    ----------
    metadata : pandas.DataFrame
        Metadata with ``sample`` and ``population`` columns.
    population : str
        Population identifier.

    Returns
    -------
    list[str]
        Sample identifiers.
    """
    if not population:
        return []
    selected = metadata.loc[metadata["population"] == population, "sample"]
    return selected.dropna().astype(str).tolist()


def write_samples(samples: Iterable[str], path: str) -> None:
    """
    Write one sample identifier per line.

    Parameters
    ----------
    samples : Iterable[str]
        Sample identifiers.
    path : str
        Output sample-list path.
    """
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as handle:
        for sample in samples:
            handle.write(f"{sample}\n")


def write_log(path: str, ref_count: int, tgt_count: int, src_count: int) -> None:
    """
    Write a small sample-list log.

    Parameters
    ----------
    path : str
        Log file path.
    ref_count : int
        Number of reference samples.
    tgt_count : int
        Number of target samples.
    src_count : int
        Number of source samples.
    """
    log_path = Path(path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as handle:
        handle.write(f"reference_samples\t{ref_count}\n")
        handle.write(f"target_samples\t{tgt_count}\n")
        handle.write(f"source_samples\t{src_count}\n")


metadata = read_metadata(str(snakemake.input.metadata))
ref_samples = collect_samples(metadata, str(snakemake.params.ref_population))
tgt_samples = collect_samples(metadata, str(snakemake.params.tgt_population))
src_samples = collect_samples(metadata, str(snakemake.params.src_population))

if not ref_samples:
    raise ValueError(f"No samples found for {snakemake.params.ref_population}")
if not tgt_samples:
    raise ValueError(f"No samples found for {snakemake.params.tgt_population}")

write_samples(ref_samples, str(snakemake.output.ref))
write_samples(tgt_samples, str(snakemake.output.tgt))
write_samples(src_samples, str(snakemake.output.src))
write_log(str(snakemake.log[0]), len(ref_samples), len(tgt_samples), len(src_samples))
